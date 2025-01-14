// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file test_FEsolve.C
  \brief Unit testing for the linear FE solver.
  \details The purpose of this file is to provide some unit/regression tests
  for various FE models, on a low level. It also covers parsing of Nastran bulk
  data files containing the FE model definition. For each model, the linear-
  elastic stiffness matrix is assembled, and solved for some external load
  (defined in the FE data file). Then the response at one node is compared
  against the reference values (u_x, u_y, u_z).
*/

#include "gtest.h"
#include "vpmCommon/FFCmdLineArgInterface.H"
#include <cstdlib>
#include <cmath>


extern "C" {
  //! From FFlLinkHandler_F.C
  int ffl_loadPart(const std::string&);
  //! From reducerInterface.C
  int solvePartDis(int,const int*,double*);
  //! From reducerInterface.C
  void initSolverArgs(int,char**,bool,bool);
  //! From reducerInterface.C
  void freeSingeltons();
}

static std::string srcdir; //!< Full path of the source directory of this test


/*!
  \brief Main program for the FE solver unit test executable.
*/

int main (int argc, char** argv)
{
  // Initialize the google test module.
  // This will remove the gtest-specific values in argv.
  ::testing::InitGoogleTest(&argc,argv);

  // Extract the source directory of the tests
  // to use as prefix for loading the FE data files
  for (int i = 1; i < argc; i++)
    if (!strncmp(argv[i],"--srcdir=",9))
    {
      srcdir = argv[i]+9;
      std::cout <<"Note: Source directory = "<< srcdir << std::endl;
      if (srcdir.back() != '/') srcdir += '/';
      srcdir += "FEparts/"; // sub-folder with FE data files
      while (++i < argc) argv[i-1] = argv[i];
      --argc;
    }

  // Initialize the command-line arguments
  initSolverArgs(argc,argv,true,true);
  ADDOPTION ("gvec",DoubleVec(),"Gravity vector in global coordinates");

  int status = 0;
  if (argc < 2) // Invoke the google test driver
    status = RUN_ALL_TESTS();
  else if (ffl_loadPart(argv[1]) == 1)
  {
    double displ[3];
    int nodeNo = argc > 2 ? atoi(argv[2]) : 1;
    status = solvePartDis(1,&nodeNo,displ);
    std::cout << argv[1] <<": dis("<< nodeNo <<") = "
              << displ[0] <<" "<< displ[1] <<" "<< displ[2] << std::endl;
  }

  // Release heap-allocated singleton objects
  freeSingeltons();
  return status;
}

//! \brief Struct with parameters to instantiate particular units tests over.
struct Case { const char* fileName; int nodeNo; double nodeDis[3]; };

//! \brief Global stream operator to print out a Case instance.
std::ostream& operator<< (std::ostream& os, const Case& c)
{
  return os << c.fileName <<":"
            << std::string(32-strlen(c.fileName),' ')
            <<"U("<< c.nodeNo <<") = "
            << c.nodeDis[0] <<" "<< c.nodeDis[1] <<" "<< c.nodeDis[2];
}

//! \brief Class describing a unit test instance.
class Solve : public testing::Test, public testing::WithParamInterface<Case> {};


/*!
  Creates a parameterized test reading a FE data file, solving it,
  and comparing the nodal response with given reference values.
  GetParam() will return a Case instance, with the actual file name, node number
  and displacement values at that node to compare the response against.
*/

TEST_P(Solve, FEmodel)
{
  double displ[3];
  ASSERT_FALSE(srcdir.empty());
  ASSERT_EQ(ffl_loadPart(srcdir+GetParam().fileName),1);
  ASSERT_EQ(solvePartDis(1,&GetParam().nodeNo,displ),0);
  std::cout << GetParam() <<": "
            << displ[0] <<" "<< displ[1] <<" "<< displ[2] << std::endl;
  double dmx = std::max(fabs(displ[0]),std::max(fabs(displ[1]),fabs(displ[2])));
  double eps = dmx > 10.0 ? dmx*1.0e-5 : 1.0e-4;
  for (int i = 0; i < 3; i++)
    EXPECT_NEAR(displ[i],GetParam().nodeDis[i],eps);
}


/*!
  Instantiate the test over a list of file names,
  node for which to extract displacements,
  and associated reference values to compare with.
*/

INSTANTIATE_TEST_CASE_P(TestBeams, Solve,  //  FE-part file               node  u_x    u_y   u_z
                        testing::Values(Case{ "01_1D_BeamElements.nas",   1251, 0.0, 0.0128, 0.0 },
                                        Case{ "02_1D2D_Hybrid.nas",       1259, 0.0, 0.0129, 0.0 },
                                        Case{ "03_2D_ShellElements.nas", 29358, 0.0, 0.0123, 0.0 },
                                        Case{ "04_3D_VolumeElements.nas", 1257, 0.0, 0.0128, 0.0 },
                                        Case{ "Tripod.nas",               5, 0.0267234, 0.0, 0.0 },
                                        Case{ "L-beam.nas",  3, 0.00806403, 0.0287016, 0.0168831 },
                                        Case{ "PBEAML.nas", 16355,-0.420529, 0.202566, 0.0087845 },
                                        Case{ "CylinderBeam.nas",            6, 0.0, 0.3450, 0.0 },
                                        Case{ "CylinderBeam-generic.nas",    6, 0.0, 0.3450, 0.0 },
                                        Case{ "CylinderBeam-shearY.nas",  6, 0.3434, 0.3450, 0.0 },
                                        Case{ "CylinderBeam-noshear.nas",    6, 0.0, 0.3434, 0.0 }));
