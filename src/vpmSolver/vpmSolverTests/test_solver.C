// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file test_solver.C

  \brief Unit tests for the dynamics solver.

  \details The purpose of this file is to provide basic unit-tests for some
  simple models, as an alternative to the regression tests in solverTests.
  We here focus more on basic functionality of the vpmSolver module, which is
  accessed via the solverInterface API, and uses google test for checking the
  behaviour and also for comparing response variables against reference values.

  \author Knut Morten Okstad, SAP SE

  \date 20 Mar 2020
*/

#include <string>
#include <vector>

#include "../solverInterface.h"
#include "gtest.h"

static std::string srcdir; //!< Full path of the source directory of this test


/*!
  \brief Main program for the dynamics solver unit test executable.
*/

int main (int argc, char** argv)
{
  // Initialize the google test module.
  // This will remove the gtest-specific values in argv.
  ::testing::InitGoogleTest(&argc,argv);

  // Extract the source directory of the tests
  // to use as prefix for loading the input files
  for (int i = 1; i < argc; i++)
    if (!strncmp(argv[i],"--srcdir=",9))
    {
      srcdir = argv[i]+9;
      std::cout <<"Note: Source directory = "<< srcdir << std::endl;
      if (srcdir.back() != '/') srcdir += '/';
      srcdir += "models/"; // sub-folder with input files
      break;
    }

  // Invoke the google test driver
  return RUN_ALL_TESTS();
}


//! \brief Struct with parameters to instantiate particular units tests over.
struct Case
{
  std::string         name;   //!< Name of fsi-file to test (without extension)
  int                 baseId; //!< BaseID of the triad to check the response at
  std::vector<double> refVar; //!< The reference response values to compare with
  bool                linear; //!< If \e true, do a linear static analysis only
};

//! \brief Global stream operator to print out a Case instance.
std::ostream& operator<< (std::ostream& os, const Case& c)
{
  os << c.name <<", "<< c.baseId <<",";
  for (double v : c.refVar) os <<" "<< v;
  return os <<", "<< std::boolalpha << c.linear;
}

//! \brief Class describing a unit test instance.
class Solve : public testing::Test, public testing::WithParamInterface<Case> {};


/*!
  \brief Creates a parameterized solver unit test.
  \details This reads the solver input file, solves the initial equilibrium
  configuration and compares the control response variables against the
  reference values, then (unless it is a linear static analysis only)
  it invokes the time step loop comparing the response after each time step.
  Finally it closes down the solver to clean up memory.

  GetParam() will return a Case instance, with the actual file name, baseID
  of the triad to compare the response at, and its reference values.
*/

TEST_P(Solve, SystemModel)
{
  std::string fsifile = "-fsifile=\"" + srcdir + GetParam().name + ".fsi\"";
  std::string resfile = "-resfile=\"" + GetParam().name;
  if (GetParam().linear) resfile += "-lin";
  resfile += ".res\"";

  // Set up all command-line options for this solver unit test.
  // Rely on the default values where sufficient (see solverInterface.C).
  std::vector<const char*> dynamicOpt = {
    "-initEquilibrium",
    "-tolInitEquil=1.0e-8",
    "-tolDispNorm=1.0e-16",
    "-tolVelNorm=1.0e-15",
    "-tolEnerSum=1.0e-16",
    "-nupdat=5",
    "-timeInc=0.1",
    "-timeEnd=0.3",
    "-allPrimaryVars-",
    "-allSecondaryVars-"
  };
  std::vector<const char*> linelOpt = {
    "-initEquilibrium",
    "-tolInitEquil=0.0",
    "-allPrimaryVars-",
    "-allSecondaryVars-"
  };
  std::vector<const char*>* fco = &dynamicOpt;
  if (GetParam().linear) fco = &linelOpt;
  std::vector<char*> args;
  args.reserve(3+fco->size());
  args.push_back(const_cast<char*>("test_solver"));
  args.push_back(const_cast<char*>(fsifile.c_str()));
  for (const char* opt : *fco)
    args.push_back(const_cast<char*>(opt));
  args.push_back(const_cast<char*>(resfile.c_str()));

  const int                  baseId = GetParam().baseId;
  const std::vector<double>& refVar = GetParam().refVar;

  // Lambda function for checking the response for current time step.
  auto&& checkResponse = [baseId,refVar](size_t ip)
  {
    // Extract response at the specified object
    double var[3];
    int nvar = getStateVar(baseId,var);
    EXPECT_GT(nvar,0);

    // Print response
    std::cout <<"Time="<< getCurrentTime() <<":";
    for (int j = 0; j < nvar; j++)
      std::cout <<" V"<< 1+j <<"="<< var[j];
    std::cout << std::endl;

    // Check it against the reference values
    for (int k = 0; k < nvar && k < 3 && ip+k < refVar.size(); k++)
      EXPECT_NEAR(var[k],refVar[ip+k],1.0e-6);
  };

  // Read input, preprocess the model and setup/solve the initial configuration
  int ierr = solverInit(args.size(),args.data());
  ASSERT_GE(ierr,0);
  checkResponse(0);

  if (!GetParam().linear)
  {
    // Time step loop
    size_t ipr = 3;
    while (solveNext(&ierr))
    {
      ASSERT_EQ(ierr,0);
      checkResponse(ipr);
      ipr += 3;
    }
    ASSERT_EQ(ierr,0);
    checkResponse(ipr);
  }

  // Simulation finished, terminate by closing down external files, etc.
  EXPECT_EQ(solverDone(),0);
}


/*!
  \brief Instantiates the unit tests Solve.SystemModel.
*/

INSTANTIATE_TEST_CASE_P(TestSystemBeam, Solve, // fsifile      baseID    u_x            u_z       linear?
                        testing::Values(Case{ "cantilever-pipe", 17, { 0.344609, 0.0, 9.992951,
                                                                       0.344609, 0.0, 9.992951,
                                                                       0.344609, 0.0, 9.992951 }, false },
                                        Case{ "cantilever-pipe", 17, { 0.345020, 0.0, 10.00000 }, true },
                                     Case{ "cantilever-noshear", 17, { 0.343424, 0.0, 10.00000 }, true }));


/*!
  \brief Creates a solver unit test for linear static analysis.
  \details This reads the solver input file, solves each time step
  and compares the control response variables against the reference values.
  Finally it closes down the solver to clean up memory.
*/

TEST(Solve, Prescribed)
{
  std::string fsifile = "-fsifile=\"" + srcdir + "cantilever-prescribed.fsi\"";

  // Set up all command-line options for this solver unit test.
  // Rely on the default values where sufficient (see solverInterface.C).
  std::vector<const char*> linelOpt = {
    "-numit=-1", // No iterations, pure linear statics
    "-timeInc=0.2",
    "-timeEnd=1.0",
    "-quasiStatic=1.0",
    "-tolDispNorm=1.0e-16",
    "-tolEnerSum=1.0e-16",
    "-allPrimaryVars-",
    "-allSecondaryVars-"
  };
  std::vector<const char*>* fco = &linelOpt;
  std::vector<char*> args;
  args.reserve(2+fco->size());
  args.push_back(const_cast<char*>("test_solver"));
  args.push_back(const_cast<char*>(fsifile.c_str()));
  for (const char* opt : *fco)
    args.push_back(const_cast<char*>(opt));

  //                             --- Triad 25 -----  -- Triad 17 --
  std::vector<double> refVar = { 0.086555, 0.0, 6.0, 0.2, 0.0, 10.0,
                                 0.173111, 0.0, 6.0, 0.4, 0.0, 10.0,
                                 0.259666, 0.0, 6.0, 0.6, 0.0, 10.0,
                                 0.346222, 0.0, 6.0, 0.8, 0.0, 10.0 };

  // Read input, preprocess the model and setup initial configuration
  int ierr = solverInit(args.size(),args.data());
  ASSERT_GE(ierr,0);

  // Time step loop
  bool more = true;
  size_t ip = 0;
  do
  {
    // Solve for next time step
    more = solveNext(&ierr);
    ASSERT_EQ(ierr,0);

    // Extract response at the specified triads, 25 and 17
    double var[6];
    int nvar = getStateVar(25,var);
    int nva2 = getStateVar(17,var+nvar);
    ASSERT_EQ(nvar,3);
    ASSERT_EQ(nva2,3);

    // Print response
    std::cout <<"Time="<< getCurrentTime() <<":";
    for (int j = 0; j < 6; j++)
      std::cout <<" V"<< 1+j <<"="<< var[j];
    std::cout << std::endl;

    // Check it against the reference values
    for (int k = 0; k < 6 && ip < refVar.size(); k++)
      EXPECT_NEAR(var[k],refVar[ip++],1.0e-6);
  }
  while (more);

  // Simulation finished, terminate by closing down external files, etc.
  ASSERT_EQ(solverDone(),0);
}
