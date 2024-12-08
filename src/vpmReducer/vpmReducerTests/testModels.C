// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file testModels.C

  \brief Hard-coded basic FE models for testing of FE part reduction.

  \details This file contains a set of classes, each one representing a simple
  FE model which is suitable for testing of the FE part reduction procedures
  without the need of reading an input file. The FE models are created by the
  respective constructors of the FFlLinkHandler sub-classes.

  \author Knut Morten Okstad

  \date 21 Apr 2020
*/

#include <array>

#include "FFaLib/FFaAlgebra/FFaMath.H"
#include "FFaLib/FFaAlgebra/FFaMat33.H"
#include "FFaLib/FFaCmdLineArg/FFaCmdLineArg.H"

#include "FFlLib/FFlLinkHandler.H"
#include "FFlLib/FFlFEParts/FFlNode.H"
#include "FFlLib/FFlFEParts/FFlBeams.H"
#include "FFlLib/FFlFEParts/FFlShells.H"
#include "FFlLib/FFlFEParts/FFlLoads.H"
#include "FFlLib/FFlFEParts/FFlRGD.H"
#include "FFlLib/FFlFEParts/FFlWAVGM.H"
#include "FFlLib/FFlFEParts/FFlBUSH.H"
#include "FFlLib/FFlFEParts/FFlCMASS.H"
#include "FFlLib/FFlFEParts/FFlPBEAMSECTION.H"
#include "FFlLib/FFlFEParts/FFlPBEAMPIN.H"
#include "FFlLib/FFlFEParts/FFlPTHICK.H"
#include "FFlLib/FFlFEParts/FFlPMAT.H"

#ifdef FF_NAMESPACE
using namespace FF_NAMESPACE;
#endif


/*!
  \brief Creates an FE model of a cantilever beam with default properties.
*/

class Cantilever : public FFlLinkHandler
{
public:
  //! \brief The constructor creates the FE model.
  //! \param[in] L Total length of the beam
  //! \param[in] nel Number of beam elements
  //! \param[in] reduce If &gt; 0, the tip node is also marked external
  //! \param[in] twoD If \e true, all internal nodes are considered 3-DOF nodes
  Cantilever(double L, int nel, int reduce, bool twoD)
  {
    FFlAttributeBase* psec = new FFlPBEAMSECTION(1);
    FFlAttributeBase* pmat = new FFlPMAT(1);
    this->addAttribute(psec);
    this->addAttribute(pmat);
    double x = 0.0;
    FFlNode* n1 = new FFlNode(1,x,0.0,0.0,1);
    this->addNode(n1);
    for (int i = 0; i < nel; i++)
    {
      x += L/nel;
      FFlNode* n2 = new FFlNode(2+i,x,0.0,0.0);
      if (twoD)
        n2->setStatus(-28); // fix the local dofs 3,4,5 in all nodes
      this->addNode(n2);
      FFlElementBase* elm = new FFlBEAM2(1+i);
      elm->setNode(1,n1);
      elm->setNode(2,n2);
      elm->setAttribute(psec);
      elm->setAttribute(pmat);
      this->addElement(elm);
      n1 = n2;
    }
    if (reduce)
      n1->setStatus(reduce);
  }
};


/*!
  \brief Creates an FE model of a pinned beam with default properties.
*/

class PinnedBeam : public FFlLinkHandler
{
public:
  //! \brief The constructor creates the FE model.
  //! \param[in] L Total length of the beam
  //! \param[in] nel Number of beam elements
  //! \param[in] ielPin Index for beam element with end release
  PinnedBeam(double L, int nel, int ielPin)
  {
    FFlAttributeBase* psec = new FFlPBEAMSECTION(1);
    FFlAttributeBase* pmat = new FFlPMAT(1);
    FFlPBEAMPIN*      ppin = new FFlPBEAMPIN(1);
    if (ielPin < 0)
      ppin->PA = 456; // Release the rotational DOFs at end 1 of the element
    else
      ppin->PB = 456; // Release the rotational DOFs at end 2 of the element
    this->addAttribute(psec);
    this->addAttribute(pmat);
    this->addAttribute(ppin);
    double x = 0.0;
    FFlElementBase* elm = NULL;
    FFlNode* n1 = new FFlNode(1,x,0.0,0.0,1);
    this->addNode(n1);
    if (ielPin < 0)
      ielPin = -ielPin;
    else if (ielPin == 0)
      ielPin = nel;
    for (int i = 0; i < nel; i++)
    {
      x += L/nel;
      FFlNode* n2 = new FFlNode(2+i,x,0.0,0.0);
      this->addNode(n2);
      elm = new FFlBEAM2(1+i);
      elm->setNode(1,n1);
      elm->setNode(2,n2);
      elm->setAttribute(psec);
      elm->setAttribute(pmat);
      if (i+1 == ielPin)
        elm->setAttribute(ppin);
      this->addElement(elm);
      n1 = n2;
    }
  }
};


/*!
  \brief Creates an FE model of a cantilever shell cylinder.
  \details Either a rigid (RGD) or flexible (WAVGM) spider is also
  created at each end with a supernode at its reference node.
*/

class CylinderShell : public FFlLinkHandler
{
public:
  //! \brief The constructor creates the FE model.
  //! \param[in] L Total length of the cylinder
  //! \param[in] R Radius of cylinder
  //! \param[in] t Shell thickness
  //! \param[in] nL Number of elements in length direction
  //! \param[in] nC Number of elements in circular direction
  //! \param[in] reduce If \e true, the tip node is also marked external
  //! \param[in] flexible If \e true, use a flexible spider at the ends
  CylinderShell(double L, double R, double t, int nL, int nC,
                bool reduce = true, bool flexible = false)
  {
    std::cout <<"   * Generating a cylinder: L="<< L <<" R="<< R <<" t="<< t
              <<" nL="<< nL <<" nC="<< nC << std::endl;

    FFlPTHICK* pthk = new FFlPTHICK(1);
    pthk->thickness = t;
    this->addAttribute(pthk);

    FFlPMAT* pmat = new FFlPMAT(1);
    pmat->youngsModule = 2.1e11;
    pmat->poissonsRatio = 0.3;
    pmat->materialDensity = 7850.0;
    this->addAttribute(pmat);

    int i, j, n = 0;
    FFlNode* n1 = new FFlNode(++n,0.0,0.0,0.0,1);
    this->addNode(n1);
    FFlNode* n2 = new FFlNode(++n,0.0,0.0,L);
    if (reduce)
      n2->setStatus(1);
    this->addNode(n2);
    for (i = 0; i <= nL; i++)
    {
      double z = (double)i*L/(double)nL;
      for (j = 0; j < nC; j++)
      {
        double theta = (double)(2*j)*acos(-1.0)/(double)nC;
        double x = R*cos(theta);
        double y = R*sin(theta);
        this->addNode(new FFlNode(++n,x,y,z));
      }
    }

    int nnod = n;
    FFlElementBase* elm = NULL;
    for (i = n = 0; i < nL; i++)
      for (j = 0; j < nC; j++)
      {
        elm = new FFlQUAD4(++n);
        elm->setAttribute(pthk);
        elm->setAttribute(pmat);
        if (j == nC-1)
        {
          elm->setNode(2,nC*i+3);
          elm->setNode(3,nC*i+3+nC);
        }
        else
        {
          elm->setNode(2,nC*i+j+4);
          elm->setNode(3,nC*i+j+4+nC);
        }
        elm->setNode(1,nC*i+j+3);
        elm->setNode(4,nC*i+j+3+nC);
        this->addElement(elm);
      }

    if (flexible)
    {
      this->addNode(new FFlNode(++nnod,0.0,0.0,0.0));
      elm = new FFlWAVGM(++n);
      elm->setNode(1,nnod);
    }
    else
    {
      elm = new FFlRGD(++n);
      elm->setNode(1,n1);
    }
    for (j = 0; j < nC; j++)
      elm->setNode(2+j,3+j);
    this->addElement(elm);

    if (flexible)
    {
      this->addNode(new FFlNode(++nnod,0.0,0.0,L));
      elm = new FFlWAVGM(++n);
      elm->setNode(1,nnod);
    }
    else
    {
      elm = new FFlRGD(++n);
      elm->setNode(1,n2);
    }
    for (j = 0; j < nC; j++)
      elm->setNode(2+j,nC*nL+3+j);
    this->addElement(elm);

    if (flexible)
    {
      elm = new FFlBUSH(++n);
      elm->setNode(1,n1);
      elm->setNode(2,nnod-1);
      this->addElement(elm);
      elm = new FFlBUSH(++n);
      elm->setNode(1,n2);
      elm->setNode(2,nnod);
      this->addElement(elm);
      elm = new FFlCMASS(++n);
      elm->setNode(1,nnod-1);
      this->addElement(elm);
      elm = new FFlCMASS(++n);
      elm->setNode(1,nnod);
      this->addElement(elm);
    }

    this->resolve();
  }
};


/*!
  \brief Static helper generating a quadrilateral element in a regular mesh.
*/

static FFlElementBase* newQuadElement (int e, int i, int j, int n1, bool parab)
{
  FFlElementBase* elm = NULL;
  if (parab)
  {
    elm = new FFlQUAD8(e);
    elm->setNode(1,n1*(i+0)+j+1);
    elm->setNode(2,n1*(i+0)+j+2);
    elm->setNode(3,n1*(i+0)+j+3);
    elm->setNode(4,n1*(i+1)+j+3);
    elm->setNode(5,n1*(i+2)+j+3);
    elm->setNode(6,n1*(i+2)+j+2);
    elm->setNode(7,n1*(i+2)+j+1);
    elm->setNode(8,n1*(i+1)+j+1);
  }
  else
  {
    elm = new FFlQUAD4(e);
    elm->setNode(1,n1*(i+0)+j+1);
    elm->setNode(2,n1*(i+0)+j+2);
    elm->setNode(3,n1*(i+1)+j+2);
    elm->setNode(4,n1*(i+1)+j+1);
  }

  return elm;
}


/*!
  \brief Static helper generating a triangular element in a regular mesh.
*/

static FFlElementBase* newTriaElement (int e, int i, int j, int n1, bool parab)
{
  FFlElementBase* elm = NULL;
  if (parab)
  {
    elm = new FFlTRI6(e);
    if (e%2 == 0)
    {
      elm->setNode(1,n1*(i+2)+j+3);
      elm->setNode(2,n1*(i+2)+j+2);
      elm->setNode(3,n1*(i+2)+j+1);
      elm->setNode(4,n1*(i+1)+j+2);
      elm->setNode(5,n1*(i+0)+j+3);
      elm->setNode(6,n1*(i+1)+j+3);
    }
    else
    {
      elm->setNode(1,n1*(i+0)+j+1);
      elm->setNode(2,n1*(i+0)+j+2);
      elm->setNode(3,n1*(i+0)+j+3);
      elm->setNode(4,n1*(i+1)+j+2);
      elm->setNode(5,n1*(i+2)+j+1);
      elm->setNode(6,n1*(i+1)+j+1);
    }
  }
  else
  {
    elm = new FFlTRI3(e);
    if (e%2 == 0)
    {
      elm->setNode(1,n1*(i+1)+j+2);
      elm->setNode(2,n1*(i+1)+j+1);
      elm->setNode(3,n1*(i+0)+j+2);
    }
    else
    {
      elm->setNode(1,n1*(i+0)+j+1);
      elm->setNode(2,n1*(i+0)+j+2);
      elm->setNode(3,n1*(i+1)+j+1);
    }
  }

  return elm;
}


/*!
  \brief Creates an FE model of a cantilever shell.
*/

class CantileverShell : public FFlLinkHandler
{
public:
  //! \brief The constructor creates the FE model.
  //! \param[in] L Total length of the cantilever
  //! \param[in] b Width of the cantilever
  //! \param[in] t Shell thickness
  //! \param[in] nL Number of elements in length direction
  //! \param[in] nb Number of elements in width direction
  //! \param[in] reduce If &gt; 0, the tip nodes are also marked external
  CantileverShell(double L, double b, double t, int nL, int nb, int reduce = 1)
  {
    std::cout <<"   * Generating a plate: L="<< L <<" b="<< b <<" t="<< t
              <<" nL="<< nL <<" nb="<< nb << std::endl;

    bool useTriangles = false, useParabolic = false;
    FFaCmdLineArg::instance()->getValue("useTriangles",useTriangles);
    FFaCmdLineArg::instance()->getValue("useParabolic",useParabolic);
    int order = useParabolic ? 2 : 1;

    FaMat33 Tlg;
    DoubleVec alpha(3,0.0);
    FFaCmdLineArg::instance()->getValue("eulerRot",alpha);
    std::cout <<"     alpha =";
    for (double a : alpha) std::cout <<" "<< a;
    std::cout <<" ==> Tlg:"<< Tlg.eulerRotateZYX(FaVec3(alpha.data())*M_PI/180.0);

    FFlPTHICK* pthk = new FFlPTHICK(1);
    pthk->thickness = t;
    this->addAttribute(pthk);

    FFlPMAT* pmat = new FFlPMAT(1);
    pmat->youngsModule = 2.1e11;
    pmat->poissonsRatio = 0.3;
    pmat->materialDensity = 7850.0;
    this->addAttribute(pmat);

    nL *= order;
    nb *= order;
    int nnb = nb+1;

    int i, j, k, e;
    FFlNode* nn = NULL;
    for (i = 0; i <= nL; i++)
      for (j = 0; j <= nb; j++)
      {
        FaVec3 X((double)i*L/(double)nL,(double)j*b/(double)nb,0.0);
        this->addNode(nn = new FFlNode(nnb*i+j+1,Tlg*X));
        if (i == 0 && j%nb == 0)
          nn->setStatus(1);
        else if (i == nL && j%nb == 0 && reduce)
          nn->setStatus(reduce);
      }

    FFlElementBase* elm = NULL;
    for (i = e = 0; i < nL; i += order)
      for (j = 0; j < nb; j += order)
        if (!useTriangles)
        {
          elm = newQuadElement(++e,i,j,nnb,useParabolic);
          elm->setAttribute(pthk);
          elm->setAttribute(pmat);
          this->addElement(elm);
        }
        else for (k = 0; k < 2; k++)
        {
          elm = newTriaElement(++e,i,j,nnb,useParabolic);
          elm->setAttribute(pthk);
          elm->setAttribute(pmat);
          this->addElement(elm);
        }

    this->resolve();
  }
};


/*!
  \brief Creates an FE model of one quarter of a shell cylinder.
*/

class QuartCylShell : public FFlLinkHandler
{
public:
  //! \brief The constructor creates the FE model.
  //! \param[in] L Total length of the cylinder
  //! \param[in] R Radius of cylinder
  //! \param[in] t Shell thickness
  //! \param[in] nL Number of elements in length direction
  //! \param[in] nC Number of elements in circular direction
  //! \param[in] reduce If \e true, the tip nodes are also marked external
  QuartCylShell(double L, double R, double t, int nL, int nC, bool reduce)
  {
    std::cout <<"   * Generating shell cylinder: L="<< L <<" R="<< R <<" t="<< t
              <<" nL="<< nL <<" nC="<< nC << std::endl;

    bool useTriangles = false, useParabolic = false;
    FFaCmdLineArg::instance()->getValue("useTriangles",useTriangles);
    FFaCmdLineArg::instance()->getValue("useParabolic",useParabolic);
    int order = useParabolic ? 2 : 1;

    FFlPTHICK* pthk = new FFlPTHICK(1);
    pthk->thickness = t;
    this->addAttribute(pthk);

    FFlPMAT* pmat = new FFlPMAT(1);
    pmat->youngsModule = 2.1e11;
    pmat->poissonsRatio = 0.3;
    pmat->materialDensity = 7850.0;
    this->addAttribute(pmat);

    nL *= order;
    nC *= order;
    int nnC = nC+1;

    int i, j, k, e;
    FFlNode* nn = NULL;
    for (i = 0; i <= nL; i++)
    {
      double z = (double)i*L/(double)nL;
      for (j = 0; j <= nC; j++)
      {
        double theta = (double)j*acos(0.0)/(double)nC;
        double x = R*cos(theta);
        double y = R*sin(theta);
        this->addNode(nn = new FFlNode(nnC*i+j+1,x,y,z));
        if ((i == 0 || i == nL) && (j == 0 || j == nC))
          if (reduce || j == 0)
            nn->setStatus(1);
      }
    }

    FFlElementBase* elm = NULL;
    for (i = e = 0; i < nL; i += order)
      for (j = 0; j < nC; j += order)
        if (!useTriangles)
        {
          elm = newQuadElement(++e,i,j,nnC,useParabolic);
          elm->setAttribute(pthk);
          elm->setAttribute(pmat);
          this->addElement(elm);
        }
        else for (k = 0; k < 2; k++)
        {
          elm = newTriaElement(++e,i,j,nnC,useParabolic);
          elm->setAttribute(pthk);
          elm->setAttribute(pmat);
          this->addElement(elm);
        }

    this->resolve();
  }
};


/*!
  \brief Static helper calculating a status code from boundary condition flags.
*/

static int BC (const std::array<int,6>& flags)
{
  int c = 0, d = 1;
  for (int f : flags)
  {
    c += f*d;
    d *= 2;
  }

  return -c;
}


/*!
  \brief Creates an FE model of the Pinched diaphragmed cylinder test case.

  \details This class defines the FE model of the Pinched Diaphragmed Cylinder
  benchmark, parameterized only by the number of elements in length- and
  circular direction (all other properties are hard-coded).

  Only one eighth of the cylinder is represented, using symmetry boundary
  conditions. Therefore, the applied load is also only one quarter of the
  equivalent load on the full cylinder.
*/

class PDICylShell : public FFlLinkHandler
{
public:
  //! \brief The constructor creates the FE model.
  //! \param[in] nL Number of elements in length direction
  //! \param[in] nC Number of elements in circular direction
  PDICylShell(int nL, int nC)
  {
    const double L = 600.0; // Total length of the cylinder
    const double R = 300.0; // Radius of cylinder

    bool useTriangles = false, useParabolic = false;
    FFaCmdLineArg::instance()->getValue("useTriangles",useTriangles);
    FFaCmdLineArg::instance()->getValue("useParabolic",useParabolic);
    int order = useParabolic ? 2 : 1;
    nL *= order;
    nC *= order;
    int nnL = nL+1;
    int nnC = nC+1;

    // Define thickness and material properties
    FFlPTHICK* pthk = new FFlPTHICK(1);
    FFlPMAT*   pmat = new FFlPMAT(1);
    pthk->thickness = 3.0;
    pmat->youngsModule = 3.0e6;
    pmat->poissonsRatio = 0.3;
    pmat->materialDensity = 0.0;
    this->addAttribute(pthk);
    this->addAttribute(pmat);

    // Define the external load
    FFlLoadBase* load = new FFlCFORCE(1);
    load->setValue(FaVec3(0.0,-0.25,0.0));
    load->setTarget(nnL*nnC);
    this->addLoad(load);

    int i, j, k, e;
    double dZ = L/(double)(2*nL);
    double dTheta = acos(0.0)/(double)nC;
    std::cout <<"\n   * Generating shell cylinder: "
              << nL <<"x"<< nC << std::endl;

    // Generate the nodal points
    FFlNode* nn = NULL;
    for (i = 0; i <= nL; i++)
      for (j = 0; j <= nC; j++)
      {
        double theta = j*dTheta;
        double x = R*cos(theta);
        double y = R*sin(theta);
        double z = i*dZ;
        this->addNode(nn = new FFlNode(nnC*i+j+1,x,y,z));

        // Define boundary conditions by setting negative nodal status codes
        // i ==  0: ux = uy = rz = 0 (the diaphragm)
        // i == nL: uz = rx = ry = 0 (XY as symmetry plane)
        // j ==  0: uy = rx = rz = 0 (XZ as symmetry plane)
        // j == nC: ux = ry = rz = 0 (YZ as symmetry plane)
        if (i == 0 && j == 0)
          nn->setStatus(BC({1,1,0,1,0,1}));
        else if (i == 0 && j == nC)
          nn->setStatus(BC({1,1,0,0,1,1}));
        else if (i == nL && j == 0)
          nn->setStatus(BC({0,1,1,1,1,1}));
        else if (i == nL && j == nC)
          nn->setStatus(BC({1,0,1,1,1,1}));
        else if (i == 0)
          nn->setStatus(BC({1,1,0,0,0,1}));
        else if (i == nL)
          nn->setStatus(BC({0,0,1,1,1,0}));
        else if (j == 0)
          nn->setStatus(BC({0,1,0,1,0,1}));
        else if (j == nC)
          nn->setStatus(BC({1,0,0,0,1,1}));
      }

    // Generate the element topology
    FFlElementBase* elm = NULL;
    for (i = e = 0; i < nL; i += order)
      for (j = 0; j < nC; j += order)
        if (!useTriangles)
        {
          elm = newQuadElement(++e,i,j,nnC,useParabolic);
          elm->setAttribute(pthk);
          elm->setAttribute(pmat);
          this->addElement(elm);
        }
        else for (k = 0; k < 2; k++)
        {
          elm = newTriaElement(++e,i,j,nnC,useParabolic);
          elm->setAttribute(pthk);
          elm->setAttribute(pmat);
          this->addElement(elm);
        }

    this->resolve();
  }
};


/*!
  \brief Creates an FE model of the Open Pinched Hemisphere test case.

  \details This class defines the FE model of the Open Pinched Hemisphere
  benchmark, parameterized only by the number of elements in azimuth- and
  circular direction (all other properties are hard-coded).

  Only one quarter of the hemisphere is represented, using symmetry boundary
  conditions. Therefore, the applied loads are also only one half of the
  equivalent load on the full hemisphere.
*/

class OPHSphShell : public FFlLinkHandler
{
public:
  //! \brief The constructor creates the FE model.
  //! \param[in] nA Number of elements in azimuth direction
  //! \param[in] nC Number of elements in circular direction
  OPHSphShell(int nA, int nC)
  {
    const double R = 10.0; // Radius of hemisphere

    bool useTriangles = false, useParabolic = false;
    FFaCmdLineArg::instance()->getValue("useTriangles",useTriangles);
    FFaCmdLineArg::instance()->getValue("useParabolic",useParabolic);
    int order = useParabolic ? 2 : 1;
    nA *= order;
    nC *= order;
    int nnC = nC+1;

    // Define thickness and material properties
    FFlPTHICK* pthk = new FFlPTHICK(1);
    FFlPMAT*   pmat = new FFlPMAT(1);
    pthk->thickness = 0.04;
    pmat->youngsModule = 6.825e7;
    pmat->poissonsRatio = 0.3;
    pmat->materialDensity = 0.0;
    this->addAttribute(pthk);
    this->addAttribute(pmat);

    // Define the external loads
    FFlLoadBase* load = new FFlCFORCE(1);
    load->setValue(FaVec3(1.0,0.0,0.0));
    load->setTarget(1);
    this->addLoad(load);
    load = new FFlCFORCE(1);
    load->setValue(FaVec3(0.0,-1.0,0.0));
    load->setTarget(nnC);
    this->addLoad(load);

    int i, j, k, e;
    double dAlpha = acos(0.0)*72.0/(double)(90*nA);
    double dTheta = acos(0.0)/(double)nC;
    std::cout <<"\n   * Generating shell hemisphere: "
              << nA <<"x"<< nC << std::endl;

    // Generate the nodal points
    FFlNode* nn = NULL;
    for (j = 0; j <= nA; j++)
      for (i = 0; i <= nC; i++)
      {
        double theta = i*dTheta;
        double alpha = j*dAlpha;
        double x = R*cos(alpha)*cos(theta);
        double y = R*cos(alpha)*sin(theta);
        double z = R*sin(alpha);
        this->addNode(nn = new FFlNode(nnC*j+i+1,x,y,z));

        // Define boundary conditions by setting negative nodal status codes
        // i ==  0: uy = rx = rz = 0 (XZ as symmetry plane)
        // j == nC: ux = ry = rz = 0 (YZ as symmatry plane)
        if (i == nC/2 && j == 0)
          nn->setStatus(BC({0,0,1,0,0,0}));
        else if (i == 0)
          nn->setStatus(BC({0,1,0,1,0,1}));
        else if (i == nC)
          nn->setStatus(BC({1,0,0,0,1,1}));
      }

    // Generate the element topology
    FFlElementBase* elm = NULL;
    for (i = e = 0; i < nA; i += order)
      for (j = 0; j < nC; j += order)
        if (!useTriangles)
        {
          elm = newQuadElement(++e,i,j,nnC,useParabolic);
          elm->setAttribute(pthk);
          elm->setAttribute(pmat);
          this->addElement(elm);
        }
        else for (k = 0; k < 2; k++)
        {
          elm = newTriaElement(++e,i,j,nnC,useParabolic);
          elm->setAttribute(pthk);
          elm->setAttribute(pmat);
          this->addElement(elm);
        }

    this->resolve();
  }
};


/*!
  \brief Creates an FE model of the Scordelis-Lo Roof test case.

  \details This class defines the FE model of the Scordelis-Lo Roof benchmark,
  parameterized only by the number of elements in length- and circular direction
  (all other properties are hard-coded). One quarter of the roof is represented,
  using symmetry boundary conditions. Apply gravity load i negative y-direction
  (using command-line option -gvec 0.0 -1.0) to run this test model.
*/

class ScoLoRoof : public FFlLinkHandler
{
public:
  //! \brief The constructor creates the FE model.
  //! \param[in] nL Number of elements in length direction
  //! \param[in] nC Number of elements in circular direction
  ScoLoRoof (int nL, int nC)
  {
    const double L = 50.0; // Total length of the roof
    const double R = 25.0; // Radius of cylindric roof
    const double A = 40.0; // Roof angle at the edge

    bool useTriangles = false, useParabolic = false;
    FFaCmdLineArg::instance()->getValue("useTriangles",useTriangles);
    FFaCmdLineArg::instance()->getValue("useParabolic",useParabolic);
    int order = useParabolic ? 2 : 1;
    nL *= order;
    nC *= order;
    int nnC = nC+1;

    // Define thickness and material properties
    FFlPTHICK* pthk = new FFlPTHICK(1);
    FFlPMAT*   pmat = new FFlPMAT(1);
    pthk->thickness = 0.25;
    pmat->youngsModule = 4.32e8;
    pmat->poissonsRatio = 0.0;
    pmat->materialDensity = 360.0;
    this->addAttribute(pthk);
    this->addAttribute(pmat);

    int i, j, k, e;
    double dZ = L/(double)(2*nL);
    double dTheta = acos(0.0)*A/(double)(90*nC);
    std::cout <<"\n   * Generating shell cylinder: "
              << nL <<"x"<< nC << std::endl;

    // Generate the nodal points
    FFlNode* nn = NULL;
    for (i = 0; i <= nL; i++)
      for (j = 0; j <= nC; j++)
      {
        double theta = (1.0-A/90.0)*acos(0.0)+j*dTheta;
        double x = R*cos(theta);
        double y = R*sin(theta);
        double z = i*dZ;
        this->addNode(nn = new FFlNode(nnC*i+j+1,x,y,z));

        // Define boundary conditions by setting negative nodal status codes
        // i ==  0: ux = uy = rz = 0 (the diaphragm)
        // i == nL: uz = rx = ry = 0 (XY as symmetry plane)
        // j == nC: ux = ry = rz = 0 (YZ as symmetry plane)
        if (i == 0 && j == nC)
          nn->setStatus(BC({1,1,0,0,1,1}));
        else if (i == nL && j == nC)
          nn->setStatus(BC({1,0,1,1,1,1}));
        else if (i == 0)
          nn->setStatus(BC({1,1,0,0,0,1}));
        else if (i == nL)
          nn->setStatus(BC({0,0,1,1,1,0}));
        else if (j == nC)
          nn->setStatus(BC({1,0,0,0,1,1}));
      }

    // Generate the element topology
    FFlElementBase* elm = NULL;
    for (i = e = 0; i < nL; i += order)
      for (j = 0; j < nC; j += order)
        if (!useTriangles)
        {
          elm = newQuadElement(++e,i,j,nnC,useParabolic);
          elm->setAttribute(pthk);
          elm->setAttribute(pmat);
          this->addElement(elm);
        }
        else for (k = 0; k < 2; k++)
        {
          elm = newTriaElement(++e,i,j,nnC,useParabolic);
          elm->setAttribute(pthk);
          elm->setAttribute(pmat);
          this->addElement(elm);
        }

    this->resolve();
  }
};


void ffl_setLink (FFlLinkHandler* link);

//! \brief Creates an FE model as a FFlLinkHandler object.
int createFEModel (int iPart, int nel, int nel2,
                   double L, double b, double t, bool twoD, bool solve)
{
  FFlNode::init();
  FFlBEAM2::init();
  FFlTRI3::init();
  FFlTRI6::init();
  FFlQUAD4::init();
  FFlQUAD8::init();
  FFlRGD::init();
  FFlWAVGM::init();
  FFlBUSH::init();
  FFlCMASS::init();
  FFlPBEAMSECTION::init();
  FFlPBEAMPIN::init();
  FFlPTHICK::init();
  FFlPMAT::init();
  FFlCFORCE::init();

  switch (iPart) {
  case 0: ffl_setLink(new Cantilever(L,nel,13,twoD)); break;
  case 1: ffl_setLink(new Cantilever(L,nel,!solve,twoD)); break;
  case 2: ffl_setLink(new PinnedBeam(L,nel,nel2)); break;
  case 3: ffl_setLink(new CylinderShell(L,b,t,nel,nel2,!solve)); break;
  case 4: ffl_setLink(new CylinderShell(L,b,t,nel,nel2,!solve,true)); break;
  case 5: ffl_setLink(new CantileverShell(L,b,t,nel,nel2,!solve)); break;
  case 6: ffl_setLink(new CantileverShell(L,b,t,nel,nel2,13)); break;
  case 7: ffl_setLink(new QuartCylShell(L,b,t,nel,nel2,!solve)); break;
  case 8: ffl_setLink(new PDICylShell(nel,nel2)); break;
  case 9: ffl_setLink(new OPHSphShell(nel,nel2)); break;
  case 10: ffl_setLink(new ScoLoRoof(nel,nel2)); break;
  default:
    std::cerr <<" *** Non-implemented FE part "<< iPart << std::endl;
    return iPart;
  }

  return 0;
}


//! \brief Releases the current FE model from memory.
void clearFEModel ()
{
  ffl_setLink(NULL);
}
