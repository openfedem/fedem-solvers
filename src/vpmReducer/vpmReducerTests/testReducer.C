// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file testReducer.C

  \brief Test program for the FE part reducer on some simple models.

  \details This file contains an alternative main program for testing the
  FE reducer on some very basic hard-coded FE parts. The purpose is to provide
  an as simple as possible means to investigate FE model manipulations, and to
  compare different ways of doing things.

  \author Knut Morten Okstad

  \date 5 Mar 2018
*/

#include <cstdlib>
#include <cstring>

#include "FFaLib/FFaCmdLineArg/FFaCmdLineArg.H"
#include "FFaLib/FFaOS/FFaFortran.H"

//! \cond DO_NOT_DOCUMENT
void cmdLineArgInit (int argc, char** argv);
void redirectOutput2Ftn (int whichFile);
int createFEModel (int iPart, int nels, int nel2,
                   double L, double b, double t,
                   bool twoD, bool solve = true);
extern "C" int ffl_loadPart (const std::string& linkFile);
void releaseSingeltons (bool lastPart = true);
void clearFEModel ();

SUBROUTINE (solver,SOLVER) (const int&, const int*, double*, int& iret);
SUBROUTINE (reducer,REDUCER) (int& iret);
SUBROUTINE (smallreducer,SMALLREDUCER) (int& iret);
SUBROUTINE (setdate,SETDATE) (const char*, const int nchar = 0);
SUBROUTINE (setversion,SETVERSION) (const char*, const int nchar = 0);
SUBROUTINE (solvesupel,SOLVESUPEL) (const int& ndim, const int& nfix,
                                    const int& neqs, double* g, int& iret);

#define ADDOPTION FFaCmdLineArg::instance()->addOption
//! \endcond


/*!
  \brief Main program for system-testing of FE model reduction.
*/

int main (int argc, char** argv)
{
  int offs = 0;
  int nels = 4, ncyl = 4;
  double L = 10.0, R = 1.0, t = 0.1;

  if (argc > 1 && argv[1][0] != '-')
  {
    // Extract number of elements and the beam length from
    // the first two command-line options
    offs = 1;
    nels = atoi(argv[1]);
    if (argc > 2 && argv[2][0] != '-')
      L = atof(argv[++offs]);
    if (nels < 1 || L <= 0.0)
      return 1;
  }

  bool slv = false;
  if (argc > 1+offs && !strcmp(argv[1+offs],"-solve"))
  {
    ++offs;
    slv = true; // We are doing a direct solve (no reduction)
  }

  bool small = false;
  if (argc > 1+offs && !strcmp(argv[1+offs],"-small"))
  {
    ++offs;
    small = true; // Use the small test driver
  }

  bool twoD = false;
  if (argc > 1+offs && !strcmp(argv[1+offs],"-2D"))
  {
    ++offs;
    twoD = true; // Consider as 2D problem
  }

  // Check which FE model constructor to use (default is the Cantilever beam)
  int ipart = 1;
  if (argc > 1+offs && !strncmp(argv[1+offs],"-P",2))
    ipart = atoi(argv[++offs]+2);

  if (ipart == 2)
  {
    // Check for element with end release
    ncyl = 0;
    if (argc > 1+offs)
      if (argv[1+offs][0] != '-' || isdigit(argv[1+offs][1]))
        ncyl = atoi(argv[++offs]); // ielpin
  }
  else if (ipart >= 3)
    // Extract cylinder properties from command-line options
    for (int i = 1; 1+offs < argc && argv[1+offs][0] != '-'; i++, offs++)
      switch (i) {
      case 1: L = atof(argv[1+offs]); break;
      case 2: R = atof(argv[1+offs]); break;
      case 3: t = atof(argv[1+offs]); break;
      case 4: nels = atoi(argv[1+offs]); break;
      case 5: ncyl = atoi(argv[1+offs]); break;
      }

  // Initialize command-line parser with the remaining options
  cmdLineArgInit (argc-offs,argv+offs);
  redirectOutput2Ftn (1);

  // Only define command-line options relevant for this test program
  ADDOPTION ("printArray",0,"Debug print switch for enumerated arrays");
  ADDOPTION ("Bmatprecision",2,"B-matrix precision");
  ADDOPTION ("linkId",1,"Link base-ID number");
  ADDOPTION ("linkfile","","Name of link input file");
  ADDOPTION ("ftlout","","Name of link output file in FTL-format");
  ADDOPTION ("samfile","","Name of SAM data file");
  ADDOPTION ("extNodes","","List of external nodes");
  ADDOPTION ("ngen",0,"Number of generalized modes");
  ADDOPTION ("nevred",12,"Number of eigenvalues to compute for reduced system");
  ADDOPTION ("dispfile","","Name of displacement file");
  ADDOPTION ("resfile","","Name of result output file");
  ADDOPTION ("frsfile","","Name of results database file");
  ADDOPTION ("rdbinc",1,"Increment number for the results database files");
  ADDOPTION ("lumpedmass",false,"Use lumped element mass matrices");
  ADDOPTION ("denseSolver",false,"Use the LAPACK equation solver");
  ADDOPTION ("skylineSolver",false,"Use the skyline equation solver");
  ADDOPTION ("tolFactorize",1.0e-12,"Equation solver singularity criterion");
  ADDOPTION ("tolEigval",1.0e-8,"Max acceptable relative error in eigenvalues");
  ADDOPTION ("factorMass",false,"Factorize mass matrix in eigenvalue solver");
  ADDOPTION ("useANDESformulation",true,"Shell formulation option");
  ADDOPTION ("useRotCpl",true,"Use rotation coupling in WAVGM's");
  ADDOPTION ("useEccCpl",true,"Use eccentricity coupling in WAVGM's");
  ADDOPTION ("autoStiffMethod",3,"Method for automatic stiffness computations");
  ADDOPTION ("autoStiffScale",1.0e2,"Scale factor for auto-added springs");
  ADDOPTION ("autoMassScale",1.0e-9,"Scale factor for auto-added masses");
  ADDOPTION ("drillingStiff",1.0e-6,"Fictitious drilling DOF stiffness");
  ADDOPTION ("useParabolic",false,"Use parabolic shell elements");
  ADDOPTION ("useTriangles",false,"Use triangular shell elements");
  ADDOPTION ("gvec",DoubleVec(),"Gravity vector in global coordinates");
  ADDOPTION ("solution",DoubleVec(),"Solution vector to compare with");
  ADDOPTION ("eulerRot",DoubleVec(),"Euler angles for orientation");
  ADDOPTION ("mute",true,"Suppress warnings on undefined options");
  FFaCmdLineArg::instance()->getValue("mute",FFaCmdLineArg::mute);

  F90_NAME(setdate,SETDATE) ("");
  F90_NAME(setversion,SETVERSION) ("");

  // Create the FE model
  int ierr = 0;
  std::string linkFile;
  FFaCmdLineArg::instance()->getValue("linkfile",linkFile);
  if (linkFile.empty())
    ierr = createFEModel(ipart,nels,ncyl,L,R,t,twoD,slv);
  else if ((ierr = ffl_loadPart(linkFile)) < 0)
    std::cerr <<" *** Failed to load FE data file "<< linkFile
              <<" ("<< ierr <<")"<< std::endl;
  else
  {
    ipart = nels = 1;
    ierr = 0;
  }
  if (ierr) return ierr;

  // Invoke the FE reducer or direct solver
  if (slv)
    F90_NAME(solver,SOLVER) (0,NULL,NULL,ierr);
  else if (small)
    F90_NAME(smallreducer,SMALLREDUCER) (ierr);
  else
    F90_NAME(reducer,REDUCER) (ierr);

  // Erase the FE model
  clearFEModel();

  if (ierr || slv || ipart == 2)
  {
    releaseSingeltons();
    return ierr;
  }

  // After successful reduction, test it by solving the superelement system
  int ndim = ipart >= 5 ? 18 + 6*(ipart%2) : (ipart == 0 ? 9 : 12);
  int nfix = ipart >= 5 ? 12 : 6;
  int neqs = ndim + (twoD ? 3 : 6)*(ipart >= 5 ? (nels+1)*(ncyl+1)-4 : nels-1);
  if (ipart == 0)
    neqs += 3;
  else if (ipart == 6)
    neqs += 6;
  DoubleVec g(3,0.0);
  FFaCmdLineArg::instance()->getValue("gvec",g);
  F90_NAME(solvesupel,SOLVESUPEL) (ndim,nfix,neqs,g.data(),ierr);

  releaseSingeltons();
  return ierr;
}
