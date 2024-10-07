// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file SolutionMapper/main.C
  \brief Solution field mapping utility for sub-model analysis.
*/

#include "FFrLib/FFrExtractor.H"
#include "FFlLib/FFlInit.H"
#include "FFlLib/FFlLinkHandler.H"
#include "FFlLib/FFlIOAdaptors/FFlReaders.H"
#include "FFlLib/FFlIOAdaptors/FFlFedemWriter.H"
#include "FFlLib/FFlFEParts/FFlShellElementBase.H"
#include "FFlLib/FFlFEParts/FFlNode.H"
#include "FFlLib/FFlElementBase.H"
#include "FFaLib/FFaAlgebra/FFaMath.H"
#include "FFaLib/FFaAlgebra/FFaMat34.H"
#include "FFaLib/FFaCmdLineArg/FFaCmdLineArg.H"
#include "FFaLib/FFaDefinitions/FFaAppInfo.H"
#include "FFaLib/FFaDefinitions/FFaResultDescription.H"
#include "FFaLib/FFaOS/FFaFilePath.H"
#include "FFaLib/FFaOS/FFaTag.H"
#include "Admin/FedemAdmin.H"
#include <fstream>
#include <iomanip>
#include <cstring>

extern void cmdLineArgInitStd(int argc, char** argv);
extern void readOptionFilesStd(const char* program);

//! \cond DO_NOT_DOCUMENT
#define ADDOPTION FFaCmdLineArg::instance()->addOption
#define GETOPTION FFaCmdLineArg::instance()->getValue
//! \endcond


//! \brief Struct with sub-model nodal data.
struct Node
{
  FaVec3          Xn;    //!< Global Cartesian coordinates
  int             elmId; //!< ID of the matching element
  FFlElementBase* elm;   //!< Pointer to the matching element
  double          xi[3]; //!< Parametric coordinates of the matching element
  IntVec          nodes; //!< Global node numbers of the matching element
  DoubleVec       shape; //!< Interpolation function values

  //! \brief The constructor initializes the coordinates.
  Node(const FaVec3& x = FaVec3()) : Xn(x)
  {
    elmId = 0;
    elm = NULL;
    xi[0] = xi[1] = xi[2] = 0.0;
  }

  //! \brief Assignment operator.
  Node& operator=(const FaVec3& x)
  {
    Xn = x;
    elmId = 0;
    elm = NULL;
    xi[0] = xi[1] = xi[2] = 0.0;
    return *this;
  }
};

typedef std::map<int,Node> NodeMap; //!< Convenience type definition.



/*!
  \brief Reads an FE model from file into a FFlLinkHandler object.
  \param[in] fileName File containing the FE model
  \return Pointer to object containing the FE model, or NULL on failure.
*/

FFlLinkHandler* readFEModel (const std::string& fileName)
{
  FFlLinkHandler* feModel = new FFlLinkHandler();
  std::cout <<"\nReading "<< fileName << std::endl;
  if (FFlReaders::instance()->read(fileName,feModel) > 0)
  {
    std::cout <<"Read done."<< std::endl;
    return feModel;
  }

  std::cerr <<" *** Read failed."<< std::endl;
  delete feModel;
  return NULL;
}


/*!
  \brief Searches for nodal points in a global FE model.
  \param ftlFile File containing the global FE model
  \param subFile File containing the FE sub-model
  \param nodFile File containing point coordinates to search for
  \param[out] feModel The global FE model
  \param[out] subModel The FE sub-model
  \param[out] nodes Resulting node-to-element mapping
  \param[out] T_lg Local to global transformation matrix for sub-model
  \return Total number of nodal points to search for, zero if read failure
*/

size_t searchPoints (std::string& ftlFile, FFlLinkHandler*& feModel,
                     std::string& subFile, FFlLinkHandler*& subModel,
                     std::string& nodFile, std::map<int,NodeMap>& nodes,
                     FaMat33& T_lg)
{
  // Evaluate some command-line options
  int iprint = 0, groupIdx = -1;
  DoubleVec trans, rot;
  double nodeTol = 0.0;
  GETOPTION ("debug",iprint);
  GETOPTION ("group",groupIdx);
  GETOPTION ("translate",trans);
  GETOPTION ("rotate",rot);
  GETOPTION ("nodeTol",nodeTol);
  GETOPTION ("offPlaneTol",FFlShellElementBase::offPlaneTol);

  // Set up the local-to-global transformation matrix for the sub-model
  if (!rot.empty())
  {
    for (double& a : rot) a *= M_PI/180.0;
    if (rot.size() < 3) rot.resize(3,0.0);
    T_lg.eulerRotateZYX(rot.data());
  }
  FaMat34 Tlg(T_lg,trans.data());
  std::cout <<"\nLocal-to-global transformation:"<< Tlg << std::endl;

  FFl::initAllReaders();
  FFl::initAllElements();

  if (!subFile.empty()) // Load the FE sub-model
    if (!(subModel = readFEModel(FFaFilePath::checkName(subFile))))
      return 0;

  // Load the sub-model node file
  char cline[128];
  std::cout <<"\nReading "<< FFaFilePath::checkName(nodFile) << std::endl;
  std::ifstream is(nodFile,std::ios::in);
  while (is.getline(cline,128))
    if (cline[0] != '#') // Ignore comment lines
    {
      FaVec3 X;
      FFlNode* node = NULL;
      bool haveCoor = false;
      int c, intPnt = 0, nodeId = 0;
      char* p = strtok(cline," ,");
      char* q = NULL;
      for (c = 0; c < 5 && p; c++, p = strtok(NULL," ,"))
        switch (c) {
        case 0:
          intPnt = strtol(p,&q,10);
          break;
        case 1:
          nodeId = strtol(p,&q,10);
          break;
        default:
          X[c-2] = strtod(p,&q);
          if (q && q[0] == '-') // The value is on the format 1.0-3 (not 1.0e-3)
            X[c-2] *= pow(10.0,atof(q)); // account for the exponent explicitly
          haveCoor = true;
          break;
        }

      // If only one value per line, consider it to be the node Id
      if (c == 1) std::swap(intPnt,nodeId);

      // Store the global nodal coordinates.
      // Use coordinates from the FE sub-model, if provided.
      if (subModel && !(node = subModel->getNode(nodeId)))
        std::cout <<"  ** Node "<< nodeId <<" is not present in the sub-model"
                  <<" (ignored)."<< std::endl;
      else if (node)
        nodes[intPnt][nodeId] = Tlg*node->getPos();
      else if (haveCoor)
        nodes[intPnt][nodeId] = Tlg*X;
    }

  // Print node summary
  size_t nNodes = 0;
  for (const std::pair<const int,NodeMap>& nodeList : nodes)
  {
    nNodes += nodeList.second.size();
    std::cout <<"Interface"<< std::setw(2) << nodeList.first <<":"
              << std::setw(5) << nodeList.second.size() <<" nodes.\n";
  }
  if (nNodes > 0)
    std::cout <<"Total :"<< std::setw(10) << nNodes <<" nodes."<< std::endl;
  else
  {
    std::cerr <<" *** No sub-model nodes specified."<< std::endl;
    return 0;
  }

  // Load the global FE model
  if (!(feModel = readFEModel(FFaFilePath::checkName(ftlFile))))
    return 0;

  // Lambda function searching for a matching node.
  auto&& findMatchingNode = [feModel,nodeTol](const FaVec3& X) -> int
  {
    FFlNode* node = feModel->findClosestNode(X);
    if (node && (X - node->getPos()).length() <= nodeTol)
      return node->getID();

    return 0;
  };

  // Search for matching elements in the global FE model
  // for each nodal point of the sub-model
  std::cout <<"Searching for matching elements ..."
            << (iprint ? "\n" : "  0%") << std::flush;
  int elmGroupID = groupIdx, nodeId;
  size_t inod = 0, ipro = 0, lpro = 0;
  for (std::pair<const int,NodeMap>& nodeList : nodes)
  {
    if (groupIdx == 0) elmGroupID = nodeList.first;
    for (std::pair<const int,Node>& node : nodeList.second)
      if ((nodeId = findMatchingNode(node.second.Xn)))
      {
        // Store ID of the found node
        node.second.elmId = -nodeId;
        node.second.nodes = { nodeId };

        if (iprint > 0)
        {
          // Print search result
          FaVec3 X = feModel->getNode(nodeId)->getPos();
          std::cout <<"   * Node "<< node.first <<": "<< node.second.Xn
                    <<" --> Node "<< nodeId
                    <<": X = "<< X <<"  distance = "
                    << (node.second.Xn-X).length() << std::endl;
        }
        else if ((ipro = 100*(++inod)/nNodes) > lpro) // Print progress update
          std::cout <<"\b\b\b\b" << std::setw(3) << (lpro = ipro)
                    <<"%" << std::flush;
      }
      else if ((node.second.elm = feModel->findPoint(node.second.Xn,
                                                     node.second.xi,
                                                     elmGroupID)))
      {
        // Store element and node IDs of the found element
        node.second.elmId = node.second.elm->getID();
        node.second.nodes.resize(node.second.elm->getNodeCount());
        for (int lnod = 1; lnod <= node.second.elm->getNodeCount(); lnod++)
          node.second.nodes[lnod-1] = node.second.elm->getNodeID(lnod);

        if (iprint > 0)
        {
          // Print search result
          FaVec3 X = node.second.elm->mapping(node.second.xi[0],
                                              node.second.xi[1]);
          std::cout <<"   * Node "<< node.first <<": "<< node.second.Xn
                    <<" --> Element "<< node.second.elmId
                    <<" xi="<< node.second.xi[0] <<" eta="<< node.second.xi[1]
                    <<": X = "<< X <<"  distance = "
                    << (node.second.Xn-X).length() << std::endl;
        }
        else if ((ipro = 100*(++inod)/nNodes) > lpro) // Print progress update
          std::cout <<"\b\b\b\b" << std::setw(3) << (lpro = ipro)
                    <<"%" << std::flush;
      }
      else if (iprint > 0)
        std::cout <<" *** Failed to find matching element for node "
                  << node.first <<": "<< node.second.Xn << std::endl;
      else
        std::cout <<"\n *** Failed to find matching element for node "
                  << node.first <<": "<< node.second.Xn
                  <<"\nSearching continues ..."
                  << std::setw(3) << ipro <<"%" << std::flush;
  }

  if (!iprint)
    std::cout << std::endl;

  return nNodes;
}


/*!
  \brief Writes the nodal mapping container to a binary file.
  \param mapFile Name of file for storage of nodal mapping results
  \param[in] cs Checksum of the FE part the nodal mapping belongs to
  \param[in] nodes Node-to-element mapping container
  \param[in] T_lg Local to global transformation matrix for sub-model
*/

bool writeMappingFile (std::string& mapFile, unsigned int cs,
                       const std::map<int,NodeMap>& nodes,
                       const FaMat33& T_lg)
{
  // Open binary file for output
  FT_FILE fDes = FT_open(FFaFilePath::checkName(mapFile).c_str(),FT_wb);
  if (fDes <= (FT_FILE)0)
  {
    perror(mapFile.c_str());
    return false;
  }

  int i, j, k, nnod = 0;
  bool rotate = !T_lg.isCoincident(FaMat33());

  // Count the number of nodes with matching elements
  for (const std::pair<const int,NodeMap>& nodeList : nodes)
    for (const std::pair<const int,Node>& node : nodeList.second)
      if (node.second.elm || node.second.elmId < 0) nnod++;

  std::cout <<"\nWriting nodal mapping file "<< mapFile << std::endl;

  // Write file header identifying this file type
  const char* fTag = "#FEDEM nodal mapping";
  std::string cNod = ";7.3;" + std::to_string(nnod) + ";\n";
  if (rotate) cNod[3] = '9';
  if (FFaTag::write(fDes,fTag,strlen(fTag),cs,LEN_TAG) < 0 ||
      FT_write(cNod.c_str(),1,cNod.size(),fDes) < 8)
  {
    perror(mapFile.c_str());
    FT_close(fDes);
    return false;
  }

  // Store transformation matrix in case the sub-model is rotated
  size_t nBytes = FT_tell(fDes);
  for (i = 0; i < 3 && rotate; i++)
    nBytes += FT_write(T_lg[i].getPt(),sizeof(double),3,fDes);

  // Loop over all nodes with a matching element
  for (const std::pair<const int,NodeMap>& nodeList : nodes)
    for (const std::pair<const int,Node>& node : nodeList.second)
      if (node.second.elm)
      {
        // Gather integer data for this nodal point
        int nen = node.second.elm->getNodeCount();
        IntVec elmId(4+nen,0);
        elmId[0] = nodeList.first;
        elmId[1] = node.first;
        elmId[2] = node.second.elm->getID();
        elmId[3] = nen;

        // Extract the interpolation coefficients N(xi,eta) from the
        // element by inserting unit displacements at the element nodes.
        // This way we don't have to perform the element-type dependent
        // basis function evaluation when using the written file.
        DoubleVec weights(nen,0.0);
        std::vector<FaVec3> elmW(nen);
        for (i = j = k = 0; i < nen; i++)
        {
          elmId[4+i] = node.second.elm->getNodeID(1+i);
          elmW[i][k] = 1.0;
          if (++k == 3 || i+1 == nen)
          {
            FaVec3 w = node.second.elm->interpolate(node.second.xi,elmW);
            for (k = 0; k < 3 && j < nen; j++, k++)
            {
              weights[j] = w[k];
              elmW[j][k] = 0.0;
            }
            k = 0;
          }
        }
#ifdef FT_DEBUG
        for (int id : elmId) std::cout <<" "<< id;
        std::cout <<": "<< node.second.xi[0] <<" "<< node.second.xi[1] <<",";
        for (double w : weights) std::cout <<" "<< w;
        std::cout << std::endl;
#endif

        // Write data to the binary file
        nBytes += FT_write(elmId.data(),sizeof(int),4+nen,fDes);
        nBytes += FT_write(node.second.xi,sizeof(double),3,fDes);
        nBytes += FT_write(weights.data(),sizeof(double),nen,fDes);
      }
      else if (node.second.elmId < 0)
      {
        // Direct nodal match
        int elmId[5];
        elmId[0] = nodeList.first;
        elmId[1] = node.first;
        elmId[2] = node.second.elmId;
        elmId[3] = 1;
        elmId[4] = node.second.nodes.front();
#ifdef FT_DEBUG
        for (i = 0; i < 5; i++) std::cout <<" "<< elmId[i];
        std::cout << std::endl;
#endif

        // Write data to the binary file
        nBytes += FT_write(elmId,sizeof(int),5,fDes);
      }

  std::cout << std::setw(16) << nBytes <<" bytes written."<< std::endl;
  FT_close(fDes);
  return true;
}


/*!
  \brief Reads the nodal mapping results from binary file.
  \param mapFile Name of file for storage of nodal mapping results
  \param[out] cs Checksum of the FE part the nodal mapping belongs to
  \param[out] nodes Resulting node-to-element mapping
  \param[out] T_lg Local to global transformation matrix for sub-model
  \return Total number of nodes in the \a nodes container
*/

size_t readMappingFile (std::string& mapFile, unsigned int& cs,
                        std::map<int,NodeMap>& nodes, FaMat33& T_lg)
{
  // Open binary file for input
  FT_FILE fDes = FT_open(FFaFilePath::checkName(mapFile).c_str(),FT_rb);
  if (fDes <= (FT_FILE)0)
  {
    perror(mapFile.c_str());
    return 0;
  }

  std::cout <<"\nReading nodal mapping file "<< mapFile << std::endl;

  // Lambda function for writing error message and closing the file.
  auto&& failure = [fDes](const std::string& msg) -> size_t
  {
    std::cerr <<" *** "<< msg << std::endl;
    FT_close(fDes);
    return 0;
  };

  // Read file header and check that it is the right file
  char buf[16];
  float currVer;
  std::string myTag;
  int nnod = FFaTag::read(fDes,myTag,cs,LEN_TAG);
  if (nnod != FFaTag::endian() || myTag.find("#FEDEM nodal mapping") > 0)
    return failure("Invalid mapping file, tag = " + myTag);
  else if (!FT_gets(buf,16,fDes))
    return failure("Error reading file version number");
  else if (sscanf(buf,";%f;%d;",&currVer,&nnod) < 2)
    return failure("Erroneous file version field: " + std::string(buf));
  else if (currVer < 7.3f)
    return failure("Wrong file version: " + std::to_string(currVer));
  else if (fabsf(currVer-7.9f) < 0.1f) // Rotation matrix is stored
  {
    for (int i = 0; i < 3; i++)
      if (FT_read(T_lg[i].getPt(),sizeof(double),3,fDes) < 3)
        return failure("Failed to read rotation matrix");
  }

  int iprint = 0;
  GETOPTION ("debug",iprint);

  int inod, id[4] = {0,0,0,0}; // Read mapping data for each node
  for (inod = 0; inod < nnod; inod++)
    if (FT_read(id,sizeof(int),4,fDes) < 4 || id[3] < 1 || id[3] > 20)
      nnod = 0;
    else
    {
      Node& node = nodes[id[0]][id[1]];
      node.elmId = id[2];
      size_t nen = id[3];
      node.nodes.resize(nen);
      node.shape.resize(nen);
      if (FT_read(node.nodes.data(),sizeof(int),nen,fDes) < nen)
        nnod = 0;
      else if (nen == 1)
        node.shape.front() = 1.0;
      else if (FT_read(node.xi,sizeof(double),3,fDes) < 3 ||
               FT_read(node.shape.data(),sizeof(double),nen,fDes) < nen)
        nnod = 0;

      if (iprint < 1) continue;
      for (int i = 0; i < 4; i++) std::cout <<" "<< id[i];
      for (int jnod : node.nodes) std::cout <<" "<< jnod;
      std::cout <<": xi = "<< node.xi[0] <<" "<< node.xi[1] <<", N =";
      for (double w : node.shape) std::cout <<" "<< w;
      std::cout << std::endl;
    }

  if (nnod > 0)
    std::cout <<"Read done."<< std::endl;
  else
  {
    std::cerr <<" *** Failed to read from mapping file"
              <<"\n     Last node "<< inod <<":";
    for (int i = 0; i < 4; i++) std::cerr <<" "<< id[i];
    std::cerr << std::endl;
  }

  FT_close(fDes);
  return nnod;
}


/*!
  \brief Reads an FE model from file and writes out with external node status.
  \param subFile Name of file with the FE sub-model
  \param model The FE sub-model
  \param nodes List of external nodes
*/

bool writeSubModel (std::string& subFile, FFlLinkHandler*& model, IntVec& nodes)
{
  if (!model) // Load the FE sub-model
    if (!(model = readFEModel(FFaFilePath::checkName(subFile))))
      return false;

  // Flag all interface nodes as external
  int unknown = 0;
  for (int& nodeID : nodes)
  {
    FFlNode* node = model->getNode(nodeID);
    if (node)
    {
      if (node->isRefNode())
      {
        node = model->createAttachableNode(node,node->getPos());
        nodeID = node->getID();
      }
      node->setExternal();
    }
    else if (++unknown <= 10)
      std::cerr <<" *** Node "<< nodeID <<" is not present in the sub-model."
                << std::endl;
  }
  if (unknown > 0)
  {
    std::cerr <<" *** A total of "<< unknown <<" sub-model nodes do not exist"
              << std::endl;
    return false;
  }

  // Write new ftl-file with external nodes
  std::string name = FFaFilePath::getBaseName(subFile,true);
  if (!FFaFilePath::hasPath(subFile) && FFaFilePath::isExtension(subFile,"ftl"))
    name.append("_ext");
  subFile = name + ".ftl";

  FFaAppInfo current;
  std::vector<std::string> meta(4);
  meta[0] = "Fedem version: " + current.version;
  meta[1] = "Created by solution mapping utility, with external node status";
  meta[2] = "This file: " + subFile;
  meta[3] = "Written by: " + current.user + ", " + current.date;

  std::cout <<"\nWriting "<< subFile << std::endl;
  FFlFedemWriter writer(model);
  return writer.write(subFile,true,true,meta);
}


/*!
  \brief Main program for the solution field mapping utility.

  \callgraph
*/

int main (int argc, char** argv)
{
  // Initialize the command-line parser, and
  // define all command-line options for this program
  cmdLineArgInitStd (argc,argv);
  ADDOPTION ("ftlFile","","Name of FE data file for global model");
  ADDOPTION ("frsFile","","Name of results database file for global model");
  ADDOPTION ("subFile","","Name of FE data file for sub-model");
  ADDOPTION ("nodeFile","","Name of nodal coordinate file for sub-model");
  ADDOPTION ("translate",DoubleVec(3,0.0),"Global translation of sub-model");
  ADDOPTION ("rotate",DoubleVec(3,0.0),"Global rotation of sub-model");
  ADDOPTION ("nodeTol",1.0e-6,"Nodal match tolerance");
  ADDOPTION ("offPlaneTol",0.1,"Out-of-plane point match tolerance");
  ADDOPTION ("rotations",false,"Map nodal rotations also");
  ADDOPTION ("mapFile","","Name of file with nodal point mapping results");
  ADDOPTION ("outFile","","Nodal displacement file");
  ADDOPTION ("useDeformation",false,"Read deformation instead of"
             " Total translation and rotation");
  ADDOPTION ("group",-1,"Element group index to search in");
  ADDOPTION ("partId",0,"FE part baseID in results database");
  ADDOPTION ("debug",0,"Debug print switch");

  // Read the command-line option files, if any
  readOptionFilesStd ("fedem_solmap");

  // Print out program name, version and build date
  std::cout <<"\nFedem Solution Mapping utility "
            << FedemAdmin::getVersion() <<" "<< FedemAdmin::getBuildDate();
  if (!FedemAdmin::is64bit()) std::cout <<"   (32bit)";
  std::cout << std::endl;

  FFrExtractor* res = NULL;
  FFlLinkHandler* globalModel = NULL;
  FFlLinkHandler* subModel = NULL;

  // Lambda function that prints an error message and closes down.
  auto&& finished = [&res,&globalModel,&subModel](const char* msg = NULL)
  {
    if (msg) std::cerr <<"\n *** "<< msg <<".\n";
    delete res;
    delete globalModel;
    delete subModel;
    FFrExtractor::releaseMemoryBlocks();
    FFaCmdLineArg::removeInstance();
    FFl::releaseAllReaders();
    FFl::releaseAllElements();
    return msg ? 1 : 0;
  };

  // Evaluate the command-line options
  int baseID = 0;
  bool rotations = false;
  std::string ftlFile, frsFile, subFile, nodFile, mapFile, outFile;
  GETOPTION ("partId",baseID);
  GETOPTION ("rotations",rotations);
  GETOPTION ("ftlFile",ftlFile);
  GETOPTION ("frsFile",frsFile);
  GETOPTION ("subFile",subFile);
  GETOPTION ("nodeFile",nodFile);
  GETOPTION ("mapFile",mapFile);
  GETOPTION ("outFile",outFile);

  FaMat33 Tlg;
  std::map<int,NodeMap> nodes;
  const char* wmsg = NULL;
  unsigned int globCS = 0;
  size_t nNodes = 0;
  if (!nodFile.empty())
  {
    // Establish mapping from the sub-model nodes to global model elements
    if (!(nNodes = searchPoints(ftlFile,globalModel,
                                subFile,subModel,nodFile,nodes,Tlg)))
      return finished("Failure");

    // Write the mapping results to file
    globCS = globalModel->calculateChecksum();
    if (!mapFile.empty() && !writeMappingFile(mapFile,globCS,nodes,Tlg))
      return finished("Failure");
  }
  else if (mapFile.empty())
    return finished("No nodal mapping file specified");
  else if (!(nNodes = readMappingFile(mapFile,globCS,nodes,Tlg)))
    return finished("Failure");

  // Extract the IDs of the external sub-model nodes,
  // i.e., those nodes that found a matching element in the global model
  IntVec disNodes;
  disNodes.reserve(nNodes);
  for (std::pair<const int,NodeMap>& nodeList : nodes)
    for (std::pair<const int,Node>& node : nodeList.second)
      if (!node.second.nodes.empty())
        disNodes.push_back(node.first);

  // Write out an FTL-file with external node status for the sub-model
  if (!subFile.empty() && !writeSubModel(subFile,subModel,disNodes))
    return finished("Failure");

  if (frsFile.empty())
    return finished(nodFile.empty() ? "No frs-file specified" : NULL);

  // Open global model results database file
  res = new FFrExtractor("RDB reader");
  if (res->addFile(FFaFilePath::checkName(frsFile),true))
    std::cout <<"\nOpened result database file "<< frsFile << std::endl;
  else
    return finished("Failed to load specified results database file");

  // Initialize the nodal results reader
  typedef std::pair<FFrEntryBase*,FFrEntryBase*> NodeRes;
  std::map<int,NodeRes> nodeRes;
  for (std::pair<const int,NodeMap>& nodeList : nodes)
    for (std::pair<const int,Node>& node : nodeList.second)
      for (int inod : node.second.nodes)
        nodeRes[inod] = NodeRes(NULL,NULL);

  // Search for the time step variables
  FFrEntryBase* stepVar = res->getTopLevelVar("Time step number");
  FFrEntryBase* timeVar = res->getTopLevelVar("Physical time");
  if (!timeVar || !stepVar)
    return finished("No time step variables in file");

  // Search for the nodal result variables
  size_t nNodeVar = 0;
  FFaResultDescription dDesc("Part",baseID);
  dDesc.varDescrPath = { "Nodes", "", "Dynamic response", "Total translation" };
  dDesc.varRefType   =   "VEC3";
  std::string totRot("Total rotation");

  bool useDeformation = false;
  GETOPTION ("useDeformation",useDeformation);
  if (useDeformation)
  {
    dDesc.varDescrPath.back() = "Translational deformation";
    totRot = "Angular deformation";
  }

  for (std::pair<const int,NodeRes>& node : nodeRes)
  {
    dDesc.varDescrPath[1] = std::to_string(node.first);
    if ((node.second.first = res->search(dDesc)))
      ++nNodeVar;
    else
      std::cout <<"  ** Failed to locate "<< dDesc.getText() << std::endl;
    if (rotations)
    {
      std::swap(dDesc.varDescrPath.back(),totRot);
      if (!(node.second.second = res->search(dDesc)))
        std::cout <<"  ** Failed to locate "<< dDesc.getText() << std::endl;
      std::swap(dDesc.varDescrPath.back(),totRot);
    }
  }

  std::cout <<"   * Successfully located "<< nNodeVar
            <<" nodes in the results database file";
  if (nNodeVar < nodeRes.size())
  {
    std::cout <<" (expected "<< nodeRes.size() <<").";
    wmsg = "Not all nodes were assigned valid displacements.\n     "
           "Review the messages and consider fixing before continuing";
  }
  std::cout << std::endl;

  if (!res->resetRDBPositioning())
    return finished("Invalid results data base file, no time steps stored");

  // Open binary file for output
  FT_FILE fileDes = FT_open(FFaFilePath::checkName(outFile).c_str(),FT_wb);
  if (fileDes <= (FT_FILE)0)
  {
    perror(outFile.c_str());
    return finished("Unable to open output file");
  }

  // Write file header identifying this file type
  const char* fileTag = "#FEDEM nodal displacements";
  const char* fileVer = rotations ? ";7.6;" : ";7.3;";
  if (FFaTag::write(fileDes,fileTag,strlen(fileTag),globCS,LEN_TAG) < 0 ||
      FT_write(fileVer,1,5,fileDes) < 5)
  {
    perror(outFile.c_str());
    FT_close(fileDes);
    return finished("Unable to write file header");
  }

  std::cout <<"\nWriting nodal displacement file "<< outFile << std::endl;

  int iprint = 0;
  GETOPTION ("debug",iprint);

  size_t nBytes = FT_tell(fileDes);

  // Lambda function that writes data to the binary file and checks for error.
  auto&& writeData = [fileDes,&nBytes,&wmsg](const void* data,
                                             size_t n, size_t m = 1)
  {
    int nWrote = FT_write(data,n,m,fileDes);
    if (nWrote > 0)
    {
      nBytes += n*nWrote;
      return true;
    }
    wmsg = "Write failure";
    return false;
  };

  FaMat33 Tgl(Tlg.transpose()); // Global-to-local transformation

  // Loop over all time steps in the RDB-file
  bool firstStep = true;
  do {
    int step;
    double time, data[3];
    std::map<int,FaVec3> nodeDis, nodeRot;
    DoubleVec displ;
    displ.reserve((rotations ? 3 : 6)*nNodes);

    // Read time step data
    if (res->getSingleTimeStepData(stepVar,&step,1) != 1) break;
    if (res->getSingleTimeStepData(timeVar,&time,1) != 1) break;
    std::cout <<"   * Time step "<< step <<" t="<< time << std::endl;

    // Read nodal displacements
    for (const std::pair<const int,NodeRes>& node : nodeRes)
      if (node.second.first)
        if (res->getSingleTimeStepData(node.second.first,data,3) == 3)
        {
          nodeDis[node.first] = Tgl*FaVec3(data);
          if (node.second.second)
            if (res->getSingleTimeStepData(node.second.second,data,3) == 3)
              nodeRot[node.first] = Tgl*FaVec3(data);
        }

    std::cout <<"     Read displacements for "<< nodeDis.size()
              <<" nodes."<< std::endl;
#ifdef FT_DEBUG
    for (const std::pair<const int,FaVec3>& node : nodeDis)
    {
      std::cout << node.first <<": "<< node.second;
      std::map<int,FaVec3>::const_iterator rit = nodeRot.find(node.first);
      if (rit != nodeRot.end()) std::cout <<" "<< rit->second;
      std::cout << std::endl;
    }
#endif

    // Interpolate the nodal displacements onto the sub-model nodes
    for (std::pair<const int,NodeMap>& nodeList : nodes)
      for (std::pair<const int,Node>& node : nodeList.second)
        if (!node.second.nodes.empty())
        {
          FaVec3 dis, rot;
          std::vector<FaVec3> elmDis, elmRot;
          elmDis.reserve(node.second.nodes.size());
          elmRot.reserve(node.second.nodes.size());
          for (int inod : node.second.nodes)
          {
            std::map<int,FaVec3>::const_iterator nit = nodeDis.find(inod);
            if (nit == nodeDis.end()) break;
            elmDis.push_back(nit->second);
            if (rotations)
            {
              nit = nodeRot.find(inod);
              elmRot.push_back(nit == nodeRot.end() ? FaVec3() : nit->second);
            }
          }
          if (elmDis.size() < node.second.nodes.size())
          {
            std::cout <<"  ** Missing nodal results for ";
            if (node.second.elmId > 0)
              std::cout <<"element "<< node.second.elmId
                        <<", node "<< node.first;
            else
              std::cout <<"node "<< node.first <<", it";
            std::cout <<" will be prescribed to zero."<< std::endl;
            if (iprint > 0)
            {
              std::cout <<"     The missing nodes are:";
              for (int inod : node.second.nodes)
                if (nodeDis.find(inod) == nodeDis.end())
                  std::cout <<" "<< inod;
              std::cout << std::endl;
            }
          }
          else if (node.second.elm)
          {
            dis = node.second.elm->interpolate(node.second.xi,elmDis);
            if (rotations)
              rot = node.second.elm->interpolate(node.second.xi,elmRot);
          }
          else if (elmDis.size() == 1)
          {
            dis = elmDis.front();
            if (rotations)
              rot = elmRot.front();
          }
          else
          {
            for (size_t i = 0; i < elmDis.size(); i++)
              dis += node.second.shape[i]*elmDis[i];
            for (size_t j = 0; j < elmRot.size(); j++)
              rot += node.second.shape[j]*elmRot[j];
            if (iprint > 1)
            {
              std::cout <<"\nNode"<< std::setw(7) << node.first <<":";
              for (double w : node.second.shape) std::cout <<" "<< w;
              for (int i = 0; i < 3; i++)
              {
                std::cout <<"\n"<< std::setw(11) << dis[i] <<":";
                for (const FaVec3& v : elmDis) std::cout <<" "<< v[i];
              }
              std::cout << std::endl;
            }
          }
          displ.insert(displ.end(),dis.getPt(),dis.getPt()+3);
          if (rotations)
            displ.insert(displ.end(),rot.getPt(),rot.getPt()+3);
        }

    if (firstStep)
    {
      // Write node numbers to the binary file
      std::string cNodes = std::to_string(disNodes.size()) + ";\n";
      nBytes += FT_write(cNodes.c_str(),1,cNodes.size(),fileDes);
      if (!writeData(disNodes.data(),sizeof(int),disNodes.size()))
        break;
      else
        firstStep = false;
    }
    else if (displ.size() != (rotations ? 6 : 3)*disNodes.size())
    {
      // This should normally not happen
      wmsg = "Inconsistent time history, logic error";
      break;
    }

    // Write nodal results to binary file
    if (!writeData(&step,sizeof(int)) || !writeData(&time,sizeof(double)))
      break;
    else if (!writeData(displ.data(),sizeof(double),displ.size()))
      break;
    else
      std::cout << std::setw(16) << nBytes <<" bytes written."<< std::endl;
  }
  while (res->incrementRDB());
  FT_close(fileDes);

  // Close up and quit
  if (!wmsg) std::cout <<"\nDone."<< std::endl;
  return finished(wmsg);
}
