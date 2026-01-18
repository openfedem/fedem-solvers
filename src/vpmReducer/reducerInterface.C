// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file reducerInterface.C

  \brief This file contains the C++ wrappers of the FEDEM part reducer API.

  \details Each Fortran subroutine of the FE part reducer shared library that is
  supposed to be accessible for outside applications need to have their C++
  wrapper implemented in this file. The subroutines to be wrapped need to be
  declared using the \b SUBROUTINE macro, defined in the FFaFortran.H file.
  They can then be invoked using the \b F90_NAME macro defined in the same file.

  \author Knut Morten Okstad, SAP SE

  \date 4 Apr 2018
*/

#include "FFaLib/FFaCmdLineArg/FFaCmdLineArg.H"
#include "FFaLib/FFaOS/FFaFortran.H"

#if defined(win32) || defined(win64)
#include <windows.h>
//! \brief Macro for C++ binding to function in shared library on Windows
#define DLLexport(ret) extern "C" __declspec(dllexport) ret
#else
//! \brief Macro for C++ binding to function in shared library on Linux
#define DLLexport(ret) extern "C" ret
#endif

//! \cond DO_NOT_DOCUMENT
#define ADDOPTION FFaCmdLineArg::instance()->addOption
#define ADD_PRIVATE_OPTION(a,b,c) ADDOPTION(a,b,c,false)

// Declaration of Fortran subroutines in the reducer library
SUBROUTINE (solver,SOLVER) (const int& n, const int* nodes,
                            double* displ, int& iret);
SUBROUTINE (reducer,REDUCER) (int& iret);
SUBROUTINE (cfemreducer,CFEMREDUCER) (const int& numStiff, int& iret);

extern void cmdLineArgInit (int argc, char** argv);
extern void readOptionFiles (const char* program, int whichFile);
extern void releaseSingeltons (bool lastPart = true);
extern void resetCWD ();


DLLexport(int) initSolverArgs (int argc, char** argv, bool first, bool muted)
{
  if (first)
    cmdLineArgInit(argc,argv);
  else
  {
    // Reuse the existing command-line option definitions,
    // only insert new actual values from the command-line here
    FFaCmdLineArg::init(argc,argv);
    return 0;
  }

  // Define the command-line options here
  ADDOPTION ("printArray",0,"Additional debug print switch for some arrays");
  ADDOPTION ("cachesize",0,"Cache size (KB) to be used by the SPR solver, or"
             "\ncore memory (MB) reserved for numerical data for the GSF solver"
             "\nApplies to the stiffness matrix only when lumped mass is used"
             "\n= 0: Let all numerical data be in core for the GSF solver");
  ADDOPTION ("bufsize_rigid",0,"Buffer size (in DP-words) per rigid element"
             "\n<= 0: Use conservative estimate computed internally");
  ADDOPTION ("linkId",1,"Link base-ID number");
  ADDOPTION ("linkfile","","Name of link input file (must be specified)");
  ADDOPTION ("ftlout","","Name of link output file in FTL-format");
  ADDOPTION ("samfile","","Name of SAM data file");
  ADDOPTION ("extNodes","","List of external nodes"
             "\n(in addition to those specified in the link file)");
  ADDOPTION ("resfile","","Name of result output file");
  ADDOPTION ("frsfile","","Name of results database file");
  ADDOPTION ("rdbinc",1,"Increment number for the results database files");
  ADDOPTION ("neval",0,"Number of eigenvalues/eigenvectors to compute");
  ADDOPTION ("ngen",0,"Number of generalized modes");
  ADDOPTION ("lumpedmass",false,"Use lumped element mass matrices"
             "\nDefault: Use consistent element mass matrices");
  ADDOPTION ("skylineSolver",false,"Use the skyline equation solver"
             "\nDefault: Use the SPR equation solver");
  ADDOPTION ("denseSolver",false,"Use LAPACK dense matrix equation solver"
             "\nDefault: Use the SPR equation solver");
  ADDOPTION ("gsfSolver",0,"Use the GSF equation solver"
             "\n= 0: Use SPR solver for both stiffness and mass matrix"
             "\n= 1: Use GSF solver for stiffness matrix, and SPR for mass"
             "\n= 2: Use GSF solver for both stiffness and mass matrix");
  ADDOPTION ("tolWAVGM",1.0e-4,"Geometric tolerance for WAVGM elements");
  ADDOPTION ("tolFactorize",1.0e-12,"Equation solver singularity criterion"
             "\n(smaller values are less restrictive)"
             "\nThe lowest value allowed is 1.0e-20");
  ADDOPTION ("tolEigval",1.0e-8,"Max acceptable relative error in eigenvalues");
  ADDOPTION ("eigenshift",0.0,"Shift value for eigenvalue analysis"
             "\n(target frequency for generalized DOFs)");
  ADDOPTION ("factorMass",false,"Factorize mass matrix in the eigenvalue solver"
             "\nDefault: Factorize stiffness matrix");
  ADDOPTION ("singularityHandler",1,"Option on how to treat singular matrices"
             "\n= 0: Abort on all occurring singularities"
             "\n= 1: Suppress true zero pivots, abort on reduced-to-zero pivots"
             "\n> 1: Suppress all occurring singularities of any kind");
  ADDOPTION ("autoStiffMethod",3,"Method for automatic stiffness computations"
             "\n= 1: Use k = Min(diag(K)) * 0.1/<tolFactorize>"
             "\n= 2: Use k = Mean(diag(K)) * <autoStiffScale>"
             "\n= 3: Use k = Max(diag(K)) * <autoStiffScale>");
  ADDOPTION ("autoStiffScale",1.0e2,"Scale factor for auto-added springs");
  ADDOPTION ("autoMassScale",1.0e-9,"Scale factor for auto-added masses");
  ADDOPTION ("drillingStiff",1.0e-6,"Fictitious drilling DOF stiffness");

  // The remaining options are ment for development use only
  ADDOPTION ("useRotCpl",true,"Use rotation coupling in WAVGM's",false);
  ADDOPTION ("useEccCpl",true,"Use eccentricity coupling in WAVGM's",false);
  ADDOPTION ("useUnWeightedCG",false,"Do not use weighted WAVGM CG's",false);
  ADDOPTION ("useOldWAVGM",false,"Use the old wavgmConstrEqn routine",false);
  ADDOPTION ("useOldRGD3D",false,"Use the old RGD3D routine (SAM)",false);
  ADDOPTION ("useANDESformulation",true,"Shell formulation option",false);
  ADDOPTION ("useIncompatibleModes",false,"Linear hexahedron option",false);
  ADDOPTION ("oldStyleMass",true,"Always use lumped mass for 3- and "
             "4-noded shells\nas well as for 10-noded tetrahedrons",false);

  FFaCmdLineArg::mute = muted;
  return 0; // No error conditions yet
}


DLLexport(int) solvePartDis (int n, const int* nodes, double* displ)
{
  int ierr = 0;
  readOptionFiles ("fedem_partsol",2);
  F90_NAME(solver,SOLVER) (n,nodes,displ,ierr);

  resetCWD();
  return ierr;
}


DLLexport(int) solvePart ()
{
  // Additional command-line options for the FE part solver
  ADDOPTION ("gvec",DoubleVec(),"Gravity vector in global coordinates");

  int ierr = solvePartDis(0,NULL,NULL);

  releaseSingeltons();
  return ierr;
}


DLLexport(void) freeSingeltons ()
{
  releaseSingeltons();
}


static void initReducerArgs ()
{
  // Additional command-line options for the FE part reducer
  ADDOPTION ("datacheck",false,"Do data check only (exit after data input)");
  ADDOPTION ("Bmatprecision",2,"Storage precision of the B-matrix on disk"
             "\n= 1: Single precision\n= 2: Double precision");
  ADDOPTION ("Bramsize",-1,"In-core size (MB) of displacement recovery matrix"
             "\n<= 0: Store full matrix in core");
  ADD_PRIVATE_OPTION ("dmramsize",-1,"Same as -Bramsize" // Backward compability
                      " but in terms of double words");
  ADDOPTION ("stiffile","","Name of stiffness matrix file");
  ADDOPTION ("massfile","","Name of mass matrix file");
  ADDOPTION ("gravfile","","Name of gravity force vector file");
  ADDOPTION ("dispfile","","Name of displacement vector file");
  ADDOPTION ("loadfile","","Name of load vector file");
  ADDOPTION ("Bmatfile","","Name of B-matrix file");
  ADDOPTION ("eigfile","","Name of eigenvector file");
  ADDOPTION ("nevred",12,"Number of eigenvalues to compute for reduced system");
  ADDOPTION ("nomass",false,"Skip mass matrix reduction");
  ADDOPTION ("nograv",false,"Skip gravity force and mass matrix calculation");

  // Options for non-linear link reduction
  ADD_PRIVATE_OPTION ("forcefile","","Name of force matrix file");
  ADD_PRIVATE_OPTION ("numStatesFile","","Name of file storing the number of "
                      "solved states (for nonlinear links)");
  ADD_PRIVATE_OPTION ("numCfemSolutions",0,"Number of solutions with the "
                      "nonlinear CFEM reducer");
  ADD_PRIVATE_OPTION ("CfemFile","","Name of file with additional CFEM data"
                      "\n(loads and boundary conditions)");

  // The remaining options are ment for development use only
  ADD_PRIVATE_OPTION ("diagmass",false,"Use diagonalized mass matrix"
                      "\n(this will override the -lumpedmass option)");
  ADD_PRIVATE_OPTION ("twoAssemblyLoops",false,"Assemble stiffness and mass in "
                      "separate loops\n(may save memory usage in big models)");
  ADD_PRIVATE_OPTION ("printLowBmatConnections",0,
                      "Print the given number of lowest connections between"
                      "\nexternal and internal degrees of freedom.\nUseful "
                      "for finding disconnected nodes after Guyan reduction.");
}


DLLexport(int) reducePart (bool first, bool last)
{
  if (first) initReducerArgs();

  int ierr = 0, numCfemSolutions = 0;
  readOptionFiles ("fedem_reducer",2);
  FFaCmdLineArg::instance()->getValue("numCfemSolutions",numCfemSolutions);
  if (numCfemSolutions > 0)
    F90_NAME(cfemreducer,CFEMREDUCER) (numCfemSolutions,ierr);
  else
    F90_NAME(reducer,REDUCER) (ierr);

  resetCWD();
  releaseSingeltons(last);
  return ierr;
}

//! \endcond
