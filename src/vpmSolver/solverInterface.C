// SPDX-FileCopyrightText: 2023 SAP SE
//
// SPDX-License-Identifier: Apache-2.0
//
// This file is part of FEDEM - https://openfedem.org
////////////////////////////////////////////////////////////////////////////////

/*!
  \file solverInterface.C

  \brief This file contains the C++ wrappers of the FEDEM dynamics solver API.

  \details Each Fortran subroutine or function of the dynamics solver shared
  library that is supposed to be accessible for outside applications need to
  have their C++ wrapper implemented in this file. The subroutines/functions
  to be wrapped need to be declared using the \b SUBROUTINE, \b INTEGER_FUNCTION
  and \b DOUBLE_FUNCTION macros, which all are defined in the FFaFortran.H file.
  They can then be invoked using the \b F90_NAME macro defined in the same file.

  \author Knut Morten Okstad, Fedem Technology AS

  \date 2 Dec 2016
*/

#include "vpmCommon/FFCmdLineArgInterface.H"
#include "vpmCommon/SignalHandler.H"
#include "FiDeviceFunctions/FiDeviceFunctionFactory.H"
#include "FFrLib/FFrReadOpInit.H"
#ifdef FT_HAS_RECOVERY
#include "FFlLib/FFlInit.H"
#endif
#include "FFaMathExpr/FFaMathExprFactory.H"
#include "FFaLib/FFaDefinitions/FFaMsg.H"
#include "FFaLib/FFaOS/FFaFortran.H"
#include <cstring>

#if defined(win32) || defined(win64)
#include <windows.h>
//! \brief Macro for C++ binding to function in shared library on Windows
#define DLLexport(ret) extern "C" __declspec(dllexport) ret
#else
//! \brief Macro for C++ binding to function in shared library on Linux
#define DLLexport(ret) extern "C" ret
#endif

//! \brief Internal operation status variable.
//! \details This variable is used to record the current operation state by
//! the various functions of the solver interface. It is updated inside those
//! functions only and is used to ensure that they are invoked in a valid order.
static int iop = -100;


/*!
  \brief Helper function to check correct execution order.
  \param[in] prg Name of the calling function
  \param[in] startStep If \e true, we are starting on a new time/load increment
*/
static int checkState (const char* prg, bool startStep = false)
{
  if (iop > 0 && startStep)
  {
    std::cerr <<" *** "<< prg <<": Previous step not finished yet, iop = "<< iop
              << std::endl;
    return -999;
  }
  else if (iop < -9)
  {
    std::cerr <<" *** "<< prg <<": Solver is not initialized."<< std::endl;
    return iop;
  }

  return 0;
}


/*!
  \brief Helper function to release all global heap-allocated objects on exit.
  \param[in] removeFFlSingeltons If \e true, also release FE part objects if any

  \details This function is invoked as the final task of the simulation, and
  in case of error exit, in order to clean up heap memory in case a new solver
  run is invoked later by the same process.
*/
static void releaseGlobalHeapObjects (bool removeFFlSingeltons = true)
{
  if (removeFFlSingeltons)
  {
#ifdef FT_HAS_RECOVERY
    FFl::releaseAllReaders ();
    FFl::releaseAllElements ();
#endif
    FFr::clearReadOps ();
  }
  FiDeviceFunctionFactory::removeInstance ();
  FFaMathExprFactory::removeInstance ();
  FFaCmdLineArg::removeInstance ();
  SignalHandler::reset ();
  FFaMsg::setMessager ();
  iop = -99; // Reset the operation status variable when everything is released
}


//! \cond DO_NOT_DOCUMENT
// Declaration of Fortran subroutines and functions in the solver library
SUBROUTINE (slv_init0,SLV_INIT0) (int& ierr);
SUBROUTINE (slv_inita,SLV_INITA) (const double* xinp, int& nxinp, int& ierr);
#ifdef _NCHAR_AFTER_CHARARG
SUBROUTINE (slv_initb,SLV_INITB) (const char* cfsi, const int nchar, int& ierr);
#else
SUBROUTINE (slv_initb,SLV_INITB) (const char* cfsi, int& ierr, const int nchar);
#endif
#ifdef _NCHAR_AFTER_CHARARG
SUBROUTINE (slv_initc,SLV_INITC) (const char* cfsi, const int nchar,
                                  const double* data, const int& ndat,
                                  int& ierr);
#else
SUBROUTINE (slv_initc,SLV_INITC) (const char* cfsi,
                                  const double* data, const int& ndat,
                                  int& ierr, const int nchar);
#endif
SUBROUTINE (slv_initd,SLV_INITD) (const double* data, const int& ndat,
                                  int& ierr);
SUBROUTINE (slv_restart,SLV_RESTART) (const double* data, const int& ndat,
                                      const int& writeToRDB, int& ierr);
SUBROUTINE (slv_next,SLV_NEXT) (int& iop, const int& finalStep,
                                int& done, int& ierr);
SUBROUTINE (slv_modes,SLV_MODES) (const int& nmod,
                                  double* eval, double* evec,
                                  const int& dofOrder,
                                  const int& useLaPack, int& ierr);
SUBROUTINE (slv_inverse,SLV_INVERSE) (const double* x,
                                      const int* xeqs, const int* feqs,
                                      const int& ndis, const int& nfrs,
                                      int& done, int& ierr);
SUBROUTINE (slv_done,SLV_DONE) (int& ierr);
SUBROUTINE (slv_getsmat,SLV_GETSMAT) (double* Nmat, const int& iop, int& ierr);
SUBROUTINE (slv_getemat,SLV_GETEMAT) (double* Emat, const int& bid,
                                      const int& iop, int& ierr);
SUBROUTINE (slv_rhsvec,SLV_RHSVEC) (double* Rvec, const int& iop, int& ierr);
INTEGER_FUNCTION (slv_geteqn,SLV_GETEQN) (const int& bid, int* meqn);
INTEGER_FUNCTION (slv_getvar,SLV_GETVAR) (const int& bid, double* var);
INTEGER_FUNCTION (slv_syssize,SLV_SYSSIZE) (const int& dofs);
DOUBLE_FUNCTION (slv_gettime,SLV_GETTIME) (const int& tFlag, int& ierr);
SUBROUTINE(slv_settime,SLV_SETTIME) (const double& nextTime, int& ierr);
DOUBLE_FUNCTION (slv_getfunc,SLV_GETFUNC) (const int& uid, const double& x,
                                           int& ierr);
DOUBLE_FUNCTION (slv_getfunction,SLV_GETFUNCTION) (int& ierr,
                                                   const int& uid,
                                                   const char* tag,
                                                   const int nchar);
INTEGER_FUNCTION (slv_getfuncid,SLV_GETFUNCID) (const char* tag,
                                                const int nchar);
INTEGER_FUNCTION (slv_haveresults,SLV_HAVERESULTS) ();
INTEGER_FUNCTION (slv_statesize,SLV_STATESIZE) (const int& posOnly);
INTEGER_FUNCTION (slv_gagessize,SLV_GAGESSIZE) ();
INTEGER_FUNCTION (slv_partsize,SLV_PARTSIZE) (const int& iop, const int& bid);
SUBROUTINE (slv_savestate,SLV_SAVESTATE) (const int& posOnly,
                                          double* data, const int& ndat,
                                          int& ierr);
SUBROUTINE (slv_initgages,SLV_INITGAGES) (const int& iop, double* data,
                                          const int& ndat, int& ierr);
SUBROUTINE (slv_savepart,SLV_SAVEPART) (const int& iop, const int& bid,
                                        double* data, const int& ndat,
                                        int& ierr);
SUBROUTINE (slv_straindisp,SLV_STRAINDISP) (const double* disp,
                                            const int* gageIDs, double* eps,
                                            const int& ndof, const int& ng,
                                            int& ierr);
SUBROUTINE (slv_beamforces,SLV_BEAMFORCES) (const double* disp,
                                            const int* beamIDs, double* forces,
                                            const int& ndof, const int& nBeam,
                                            int& ierr);
SUBROUTINE (slv_reldistance,SLV_RELDISTANCE) (const double* disp,
                                              const int* Ids, double* relDis,
                                              const int& ndof, const int& nId,
                                              int& ierr);
SUBROUTINE (slv_responsevars,SLV_RESPONSEVARS) (const double* disp,
                                                const int* Ids, double* Var,
                                                const int& ndof, const int& nId,
                                                int& ierr);
SUBROUTINE (slv_jointspring,SLV_JOINTSPRING) (double* sprCoeff,
                                              const int& ids, int &ierr);

extern void cmdLineArgInit (int argc, char** argv);
extern void readOptionFiles (const char* program, int whichFile = 1);
extern void resetCWD ();


DLLexport(int) solverInit (int argc, char** argv, const char* cfsi,
                           const double* stateData, const int ndat,
                           const double* gagesData, const int ngda,
                           const double* extInput, int* nxinp)
{
  // Initialize the command-line parser
  cmdLineArgInit (argc,argv);

  // Define all command-line options here
  ADDOPTION ("fsifile","fedem_solver.fsi","Name of solver input file");
  ADDOPTION ("fsi2file","","Name of additional solver input file");
  ADDOPTION ("plugin","","Plugin(s) for user-defined elements and functions");
  ADDOPTION ("skylinesolver",false,"Use skyline solver"
             "\nDefault: Use sparse solver");
  ADDOPTION ("densesolver",false,"Use LAPACK dense matrix equation solver"
             "\nDefault: Use sparse solver");
  ADDOPTION ("pardiso",false,"Use the Pardiso sparse equation solver"
             "\nDefault: Use the SPR sparse equation solver");
  ADDOPTION ("pardisoIndef",false,"Use Pardiso, indefinite system matrices"
             "\nDefault: Use the SPR sparse equation solver");
  ADDOPTION ("nosolveropt",false,"Switch off equation system reordering");
  ADDOPTION ("tolFactorize",DoubleVec({1.0e-12,1.0e-9,1.0e-12}),
             "Linear equation solver singularity critera"
             "\n(1)=time domain simulation"
             "\n(2)=initial static equilibrium and eigenvalue analysis"
             "\n(3)=control system integration"
             "\n(smaller values less restrictive)");
  ADDOPTION ("scaleToKG",1.0,"Scaling factor to SI mass unit [kg]");
  ADDOPTION ("scaleToM",1.0,"Scaling factor to SI length unit [m]");
  ADDOPTION ("scaleToS",1.0,"Scaling factor to SI time unit [s]");
  ADDOPTION ("targetFrequencyRigid",1.0e4,"Target frequency for auto-"
             "stiffness calculation");
  ADDOPTION ("restarttime",-1.0,"Physical time for restart"
             "\n< 0.0: No restart, but regular simulation");
  ADDOPTION ("restartfile","","Response database file(s) to restart from");
  ADDOPTION ("diffraction",0,"Target number of diffraction panels");
  ADDOPTION ("displacementfile","","Name of file with boundary displacements");
  ADDOPTION ("externalfuncfile","","Name of external function value file");

  // Initial static equilibrium parameters
  ADDOPTION ("initEquilibrium",false,"Initial static equilibrium iterations");
  ADDOPTION ("initEqAD",false,"Initial equilibrium with aerodynamic loads");
  ADDOPTION ("stressStiffEqu",false,"Use geometric stiffness for statics");
  ADDOPTION ("tolInitEquil",1.0e-3,"Convergence tolerance for initial"
             " equilibrium iterations");
  ADDOPTION ("limInitEquilStep",1.0,"Initial equilibrium step size limit");
  ADD_PRIVATE_OPTION ("fixSingEqInit",false,"Fix singular equations (if any)"
             "\nduring initial equilibrium iterations");
  ADD_PRIVATE_OPTION ("releaseSingEqInit",false,"Release singular equations"
             " (if any)\nduring initial equilibrium iterations"
             "\n(by adding small ficticious stiffnesses)");

  // Time integration parameters
  ADDOPTION ("timeStart",0.0,"Start time");
  ADDOPTION ("timeEnd",0.0,"Stop time");
  ADDOPTION ("timeInc",0.0,"Initial time increment");
  ADDOPTION ("alphaNewmark",DoubleVec({0.1}),"Time integration parameters");
  ADDOPTION ("alpha1",0.0,"Global mass-proportional damping factor");
  ADDOPTION ("alpha2",0.0,"Global stiffness-proportional damping factor");
  ADDOPTION ("stopGlbDmp",-1.0,"Stop time for structural damping factors");
  ADDOPTION ("quasiStatic",0.0,"Do a quasi-static simulation to this time");
  ADDOPTION ("NewmarkFlag",0,"Newmark time integration option (= cba)"
             "\na > 0: Compute inertia force from residual of previous"
             "\n       increment in the right-hand-side calculation"
             "\nb > 0: Use total solution increment in configuration update"
             "\nc = 1: Use HHT-alpha algorithm equivalent to FENRIS"
             "\nc = 2: Use generalized-alpha algorithm"
             "\n       with interpolated internal forces");
  ADDOPTION ("ignoreIC",false,"Ignore all specified initial conditions");
  ADDOPTION ("stiffDampFiltering",1,
             "Rigid body filtering of stiffness-proportional damping"
             "\n= 0: Deactivated"
             "\n= 1: Use first-order approxiation of deformational velocity"
             "\n= 2: Use second-order approxiation of deformational velocity"
             "\n     (based on Newmark update formula)");
  ADDOPTION ("num_damp_energy_skip",1,"Number of steps without calculation of"
             "\nenergy from stiffness proportional damping");
  ADDOPTION ("centripForceCorr",false,"Use centripital force correction");
  ADDOPTION ("stressStiffDyn",false,"Use geometric stiffness for dynamics");
  ADDOPTION ("stressStiffUpdateSkip",0,"Number of iterations without updating"
             " stress\nstiffening (always updated in predictor step)");
  ADDOPTION ("stressStiffDivergSkip",0,"Number of iterations without stress"
             "\nstiffening on cut-back with same step size");

  // Automatic time-stepping parameters
  ADDOPTION ("autoTimeStep",0,"Time stepping procedure"
             "\n= 0: Fixed time step size"
             "\n= 1: Automatically computed time step size");
  ADDOPTION ("maxInc",0.1,"Maximum time increment");
  ADDOPTION ("minInc",0.001,"Minimum time increment");
  ADDOPTION ("cutbackFactor",1.0,"Time step reduction factor in cut-back");
  ADDOPTION ("cutbackSteps",0,"Number of cut-back steps");
  ADDOPTION ("cutbackSing",false,"Try cut-back when detecting singularities");
  ADDOPTION ("cutbackNegPiv",100,"Try cut-back when detecting negative pivots");

  // Parameters for smooth up-ramping of first load step
  ADDOPTION ("rampSteps",0,"Number of increments in ramp-up stage");
  ADDOPTION ("rampGravity",false,"Ramp up gravity forces also");
  ADDOPTION ("rampData",DoubleVec({1.0,2.0}),"Ramp-up function parameters"
             "\n(1)=max speed during ramp-up stage"
             "\n(2)=total length (in time) of ramp-up stage"
             "\n(3)=rest time after ramp-up before new load (default=0.0)");

  // Newton-Raphson iteration parameters
  ADDOPTION ("maxit",15,"Maximum number of iterations");
  ADDOPTION ("minit",1,"Minimum number of iterations");
  ADDOPTION ("numit",0,"Fixed number of iterations");
  ADDOPTION ("nupdat",0,"Number of iterations with system matrix update");
  ADDOPTION ("maxSeqNoUpdate",100,"Max number of sequential iterations "
             "without system matrix update");
  ADDOPTION ("tolUpdateFactor",0.0,"Convergence criterion scaling factor"
             "\nfor continuing matrix updates");
  ADDOPTION ("lineSearch",false,"Use line search in the nonlinear iterations");
  ADDOPTION ("tolDispNorm",0.0,"Displacement vector convergence tolerance");
  ADDOPTION ("tolDispTra",0.0,"Max displacement tolerance");
  ADDOPTION ("tolDispRot",0.0,"Max rotation tolerance");
  ADDOPTION ("tolDispGen",0.0,"Max generalized DOF tolerance");
  ADDOPTION ("tolVelNorm",0.0,"Velocity vector convergence tolerance");
  ADDOPTION ("tolVelTra", 0.0,"Max velocity tolerance");
  ADDOPTION ("tolVelRot", 0.0,"Max angular velocity tolerance");
  ADDOPTION ("tolVelGen", 0.0,"Max generalized velocity tolerance");
  ADDOPTION ("tolAccNorm",0.0,"Acceleration vector convergence tolerance");
  ADDOPTION ("tolAccTra", 0.0,"Max acceleration tolerance");
  ADDOPTION ("tolAccRot", 0.0,"Max angular acceleration tolerance");
  ADDOPTION ("tolAccGen", 0.0,"Max generalized acceleration tolerance");
  ADDOPTION ("tolResNorm",0.0,"Residual force vector convergence tolerance");
  ADDOPTION ("tolResTra", 0.0,"Max residual force tolerance");
  ADDOPTION ("tolResRot", 0.0,"Max residual torque tolerance");
  ADDOPTION ("tolResGen", 0.0,"Max residual generalized DOF force tolerance");
  ADDOPTION ("tolEnerMax",0.0,"Max energy in a single DOF tolerance");
  ADDOPTION ("tolEnerSum",0.0,"Energy norm convergence tolerance");
  ADD_PRIVATE_OPTION ("tolVelProp",0.0,"Velocity-proportional term to the"
             "\nvelocity vector convergence tolerance");
  ADDOPTION ("monitorWorst",6,"Number of DOFs to monitor on poor convergence");
  ADDOPTION ("monitorIter",2,"Number of iterations to monitor before maxit");
  ADDOPTION ("stopOnDivergence",0,"Number of warnings on possible divergence"
             "\nbefore the dynamics simulation is aborted"
             "\n(0 = no limit)");

  FFaCmdLineArg::additionalHelpText = "\nFor all the convergence tolerance"
    " options, its value is interpreted as follows:"
    "\n = 0.0: This tolerance is ignored"
    "\n > 0.0: This tolerance is in a set of tests"
    " where all must be satisfied"
    "\n < 0.0: This tolerance is in a set of tests"
    " where only one must be satisfied"
    "\n\t(using the absolute value as the actual tolerance value)\n";

  // Eigenvalue calculation parameters
  ADDOPTION ("eiginc",0.0,"Time between each eigenvalue analysis");
  ADDOPTION ("numEigModes",0,"Number of eigenmodes to calculate");
  ADDOPTION ("lancz1",false,"Use the LANCZ1 eigensolver\nDefault: Use LANCZ2");
  ADDOPTION ("damped",false,"Solve the damped eigenproblem using LAPACK");
  ADD_PRIVATE_OPTION ("undamped",false,"Use the LAPACK eigensolver");
  ADD_PRIVATE_OPTION ("symmetric",false,"Use the symmetric LAPACK eigensolver");
  ADDOPTION ("eigenshift",0.0,"Shift value for vibration eigenvalue analysis"
             "\n(negative value captures zero frequency modes)");
  ADDOPTION ("addBC_eigensolver",false,"Use additional BCs on eigensolver");
  ADDOPTION ("factorMass_eigensolver",false,"Factor mass matrix in eigensolver"
             "\nDefault: Factor stiffness matrix");
  ADDOPTION ("stressStiffEig",false,"Use geometric stiffness for eigenvalue"
             " analysis");
  ADDOPTION ("tolEigval",1.0e-8,"Max acceptable relative error in eigenvalues");
  ADDOPTION ("tolEigvector",1.0e-8,"Orthogonality limit for the eigenvectors");
  ADDOPTION ("effModalMass",false,"Compute the effective mass for each mode");
  ADDOPTION ("yamlFile","","YAML file prefix for system mode shape export");

  // Frequency domain analysis parameters
  ADDOPTION ("frequency_domain",false,"Switch for frequency domain solution");
  ADDOPTION ("nrModes",0,"Switch between modal or direct solution, if 0 direct solution");
  ADDOPTION ("sample_freq",100.0,"Defines sampling frequency");
  ADDOPTION ("windowSize",0,"Defines the window size in samples");
  ADDOPTION ("freqfile","","Name of frequency response database file");
  ADD_PRIVATE_OPTION ("sweep_range",DoubleVec({0.0,50.0}),"Sweep frequency range");
  ADD_PRIVATE_OPTION ("sweep_input",IntVec(),"Sweep input DOFs");
  ADD_PRIVATE_OPTION ("freq_output",IntVec(),"Frequency domain analysis output DOFs");
  ADD_PRIVATE_OPTION ("pyplot",false,"Plot graphs directly with python");

  // Control system parameters
  ADDOPTION ("ctrlTol",DoubleVec({0.002,0.002,0.5,1.0e-5}),
             "Control system tolerance parameters"
             "\n(1)=Absolute iteration tolerance"
             "\n(2)=Relative iteration tolerance"
             "\n(3)=Accuracy parameter"
             "\n(4)=Relative perturbation for numerical Jacobian computation");
  ADDOPTION ("delayBuffer",1000,"Initial buffer size for delay elements");

  // General output options
  ADDOPTION ("resfile","fedem_solver.res","Name of result output file");
  ADDOPTION ("frs1file","th_p.frs","Name of primary response database file");
  ADDOPTION ("frs2file","th_s.frs","Name of secondary response database file");
  ADDOPTION ("modesfile","ev_p.frs","Name of primary modes database file");
  ADDOPTION ("ctrlfile","ctrl.frs","Name of control system database file");
  ADDOPTION ("rdbinc",1,"Increment number for the results database files");
  ADDOPTION ("rdblength",0.0,"Maximum time length of results database files");
  ADDOPTION ("VTFfile","","Name of VTF output file");
  ADDOPTION ("printFunc",0,"Option for function output."
             "\nIf 1: Print wave spectrums to result output file"
             "\nIf 2: Print all function evaluations to separate file");
  ADDOPTION ("printTriad",IntVec(),"Print some triads to result output file");
  ADDOPTION ("printinc",0.0,"Time between each print to result output file");
  ADDOPTION ("saveinc2",0.0,"Time between each save of secondary variables");
  ADDOPTION ("saveinc3",0.0,"Time between each save for external recovery");
  ADDOPTION ("saveinc4",0.0,"Time between each save of control system data");
  ADDOPTION ("savestart",0.0,"Time for first save to response database");
  ADDOPTION ("double1",true,"Save primary variables in double precision");
  ADDOPTION ("double2",false,"Save secondary variables in double precision");
  ADDOPTION ("flushinc",0.0,"Time between each database file flush"
             "\n< 0.0: Do not flush results database (let the OS decide)"
             "\n= 0.0: Flush at each time step, no external buffers"
             "\n> 0.0: Flush at specified time interval, use external buffers");

  // Options for controlling output of solution variables
  ADDOPTION ("allPrimaryVars",true,"Output all primary variables");
  ADDOPTION ("allSecondaryVars",false,"Output all secondary variables");
  ADDOPTION ("allRestartVars",false,"Output all variables needed for restart");
  ADDOPTION ("allTriadVars",false,"Output all triad variables");
  ADDOPTION ("allSupelVars",false,"Output all superelement variables");
  ADDOPTION ("allLoadVars",false,"Output all external load variables");
  ADDOPTION ("allSpringVars",false,"Output all spring variables");
  ADDOPTION ("allDamperVars",false,"Output all damper variables");
  ADDOPTION ("allFrictionVars",false,"Output all friction variables");
  ADDOPTION ("allJointVars",false,"Output all joint variables");
  ADDOPTION ("allSystemVars",false,"Output all system variables");
  ADDOPTION ("allControlVars",false,"Output all control variables");
  ADDOPTION ("allLengthVars",false,"Output all length variables");
  ADDOPTION ("allDefVars",false,"Output all deflection variables");
  ADDOPTION ("allVelVars",false,"Output all velocity variables");
  ADDOPTION ("allAccVars",false,"Output all acceleration variables");
  ADDOPTION ("allForceVars",false,"Output all force variables");
  ADDOPTION ("allStiffVars",false,"Output all spring stiffnesses");
  ADDOPTION ("allDampCoeff",false,"Output all damper coefficients");
  ADDOPTION ("allCGVars",false,"Output all centre of gravity variables");
  ADDOPTION ("allHDVars",false,"Output all hydrodynamics variables");
  ADDOPTION ("allGenDOFVars",false,"Output all generalized DOF variables");
  ADDOPTION ("allEnergyVars",false,"Output all energy quantities");
  ADDOPTION ("allAlgorVars",false,"Output all algorithm variables");
  ADDOPTION ("allTireVars",false,"Output all tire variables");
  ADDOPTION ("allEngineVars",false,"Output all engine values");
  ADDOPTION ("allBeamForces",false,"Output all beam sectional forces");
  ADDOPTION ("noBeamForces",false,"Suppress all beam sectional force output");

  // Options for automatic curve export
  ADDOPTION ("curveFile","response.bak.fmm","Name of curve definition file");
  ADDOPTION ("rpcFile","","Get number of repeats, averages, and"
             "\npoints per frame and group, from this RPC-file");
  ADDOPTION ("curvePlotFile","","Name of curve export output file");
  ADDOPTION ("curvePlotType",0,"Format of curve export output file"
             "\n= 0 : ASCII (separate file for each curve)"
             "\n= 1 : DAC, Windows (separate file for each curve)"
             "\n= 2 : DAC, UNIX (separate file for each curve)"
             "\n= 3 : RPC, Windows (all curves in one file)"
             "\n= 4 : RPC, UNIX (all curves in one file)"
             "\n= 5 : ASCII (all curves in one file)");
  ADDOPTION ("curvePlotPrec",0,"Output precision for exported curve data files"
             "\n= 0 : half precision (int*2)"
             "\n= 1 : single precision (real*4)"
             "\n= 2 : double precision (real*8)");

#ifdef FT_HAS_RECOVERY
  // Stress recovery options
  ADDOPTION ("recovery",0,"Recovery option (1=stress, 2=gage, 3=both)");
  ADDOPTION ("frs3file","","Name of stress recovery database file");
  ADDOPTION ("allGages",false,"Output all strain gages");
  ADDOPTION ("partDeformation",1,"Output recovered part deformations to frs"
             "\n(0=off, 1=local, 2=total, 3=local and total");
  ADDOPTION ("partVMStress",1,"Output recovered von Mises stresses on parts"
             "\n(0=off, 1=output to frs, 2=output through state array, 3=both)");
  ADDOPTION ("Bramsize",-1,"In-core size (MB) of displacement recovery matrix"
             "\n< 0: Use the same as in the reducer (default)"
             "\n= 0: Store full matrix in core");
  ADD_PRIVATE_OPTION ("dmramsize",-1,"Same as -Bramsize"
                      " but in terms of double words");
  ADD_PRIVATE_OPTION ("stressForm",0,"General stress formulation option");
  ADD_PRIVATE_OPTION ("ffqStressForm",2,"Stress formulation for the FFQ shell");
  ADD_PRIVATE_OPTION ("fftStressForm",1,"Stress formulation for the FFT shell");
  ADD_PRIVATE_OPTION ("useIncompatibleModes",false,"Linear hexahedron option");
#endif

  // The remaining options are for development use only
  ADD_PRIVATE_OPTION ("constantArcHeight",false,"Use constant arch height"
                      " for the contact elements");
  ADD_PRIVATE_OPTION ("nullifyEvec",false,"Nullify the eccentricity vector"
                      "\nfrom joint position to dependent triad position");
  ADD_PRIVATE_OPTION ("depDirInJointDir",false,"Make dependent joint triads"
                      "\nhave their system directions in joint directions");
  ADD_PRIVATE_OPTION ("ignoreCamEvecF",false,"Ignore moment from eccentricity for cam, force terms");
  ADD_PRIVATE_OPTION ("ignoreCamEvecFgrad",false,"Ignore moment from eccentricity for cam, Newton terms");
  ADD_PRIVATE_OPTION ("saveIter",DoubleVec({0.0,-1.0}),
                      "Save results at every iteration when "
                      "saveIter(1) < time < saveIter(2)");
  ADD_PRIVATE_OPTION ("saveForceContributions",false,"Save stiffness-, damping-,"
                      " inertia-\nand external force contributions on Triads");
  ADD_PRIVATE_OPTION ("centripCorrToFD",false,"Add centripetal moment correction to damping forces");
  ADD_PRIVATE_OPTION ("centripCorrToFI",false,"Add centripetal moment correction to inertia forces");
  ADD_PRIVATE_OPTION ("continueAfterMaxIter",false,"Continue time integration after reaching maximum iterations");
  ADD_PRIVATE_OPTION ("predControl",1,"Control system calculation point (1, 2, or 3)");
  ADD_PRIVATE_OPTION ("corrControl",1,"Control system calculation point (1, 2, or 3)");
  ADD_PRIVATE_OPTION ("nFricSkip",0,"Number of friction calculation skips");
  ADD_PRIVATE_OPTION ("nFricSkipStatic",0,"Number of friction calculation skips in quasi-static stage");
  ADD_PRIVATE_OPTION ("nFricVelUpdateSkip",0,"Number of friction velocity update skips");
  ADD_PRIVATE_OPTION ("tolFricForce",0.0,"Relative friction force tolerance");
  ADD_PRIVATE_OPTION ("fricForm",1,"Friction formulation"
                      "\n= 0 : Original"
                      "\n= 1 : Use time-step dependent predictor"
                      "\n= 2 : Use force-dependent stick friction direction"
                      "\n= 3 : Combination of 1 and 2");
  ADD_PRIVATE_OPTION ("tireIntegrator",3,"Integration scheme for tire variables"
                      "\n= 0 : Forward Euler without iterative correction"
                      "\n= 1 : Backward Euler with iterative correction"
                      "\n= 2 : Improved Euler with iterative correction"
                      "\n= 3 : Basic Theta Method, Theta = 0.505");
  ADD_PRIVATE_OPTION ("tireStaticEquil",1,
                      "Method for incorporating tires in static equilibrium"
                      "\n= 0 : Ignore tires"
                      "\n= 1 : Include tires, but ignore initial velocities"
                      "\n= 2 : Include tires, and use the initial velocities");
  ADD_PRIVATE_OPTION ("printSys",IntVec(),"Equation range for system print");
  ADD_PRIVATE_OPTION ("printSysDofs",false,"Do the system print in DOF-order");
  ADD_PRIVATE_OPTION ("noDiscontinuitySplit",false,"Suppress time step splitting on"
                      "\ncontrol system discontinuities");
  ADD_PRIVATE_OPTION ("zeroControlInit",false,
                      "Set control system initial states to zero,"
                      "\nexcept integrators that are specified by the user");
  ADD_PRIVATE_OPTION ("findAllSingularities",true,"Find all singular equations"
                      " before exiting");
  ADD_PRIVATE_OPTION ("shadowPosAlg",1,
                      "Method for positioning the shadow element C0n"
                      "\n= 1 : Triangle fit based on selected triads"
                      "\n= 2 : Mass-based weighted average of nodal positions"
                      "\n(negative values gives debug print to 'shadowPos.dbg'");
  ADD_PRIVATE_OPTION ("overrideMassScale",-1.0,"Overriding mass scale for all superelements");
  ADD_PRIVATE_OPTION ("overrideMassPropDamp",-1.0,"Overriding mass-proportional damping for all links");
  ADD_PRIVATE_OPTION ("overrideStiffPropDamp",-1.0,"Overriding stiffness-proportional damping for all links");
  ADD_PRIVATE_OPTION ("overrideTireType",-1,"Override tire type for all tires");
  ADD_PRIVATE_OPTION ("JDTalpha",0.5,"lateral, integration damping for John Deere tire");
  ADD_PRIVATE_OPTION ("JDTdampFactor",1.0,"for ILAT /=1; damping for John Deere tire");
  ADD_PRIVATE_OPTION ("JDTstiffFactor",0.0,"for ILAT /=1; stiffness for John Deere tire");
  ADD_PRIVATE_OPTION ("JDTdamping",3,"for ILAT /=1; use finite diff. damping John Deere tire");
  ADD_PRIVATE_OPTION ("JDTstiffness",1,"Stiffness matrix option for John Deere tire dynamics"
                      "\n= 0 : Ignore stiffness"
                      "\n= 1 : Include stiffness"
                      "\n= 2 : Use new stiffness");
  ADD_PRIVATE_OPTION ("resFileFormat",2,"Format on convergence output"
                      "\n= 0 : Suppress all convergence output"
                      "\n= 1 : Use old format (compatible with R4.0 and older)"
                      "\n= 2 : Use new format");
  ADD_PRIVATE_OPTION ("scaledStructDamp",false,"Structural damping matrix based on"
                      "\nscaled mass- and stiffness matrix");
  ADD_PRIVATE_OPTION ("dragForm",2,"Drag formulation");
  ADD_PRIVATE_OPTION ("slamForm",0,"Slam formulation");
  ADD_PRIVATE_OPTION ("noAddedMass",false,"Ignore all added masses");
  ADD_PRIVATE_OPTION ("dumpBeamForces",-1.0,"Print a summary of beam forces");
  ADD_PRIVATE_OPTION ("noFEnodeAtJoints",false,"Disable sectional forces in Joint triads");
  ADD_PRIVATE_OPTION ("noWheelerStretching",false,"Disable depth stretching relative to free wave surface");
  ADD_PRIVATE_OPTION ("saveBeamJoints",false,"Dump sectional forces at beam joints");
  ADD_PRIVATE_OPTION ("dumpWaveNode",0,"Dump wave kinematics at given Triad");
  ADD_PRIVATE_OPTION ("nHDupdat",1,"Number of iterations with hydrodynamics update");
  ADD_PRIVATE_OPTION ("curvature5p",false,"Use 5 point stencil in beam curvature calculations");
  ADD_PRIVATE_OPTION ("HWAFLS",false,"Use hardware module for wave kinematics");
  ADD_PRIVATE_OPTION ("FNV",0,"FNV wave force formulation option");
  ADD_PRIVATE_OPTION ("FNVcutoff",-1.0,"FNV LP cut-off frequency for the kinematics");
  ADD_PRIVATE_OPTION ("FNVlength",0.0,"FNV simulation length");
  ADD_PRIVATE_OPTION ("ignoreAD",false,"Ignore all aerodynamic loads");
  ADD_PRIVATE_OPTION ("ignoreHD",false,"Ignore all hydrodynamic loads");
  ADD_PRIVATE_OPTION ("dumpSupEls",false,"Dump superelements to Fortran file and exit");
  ADD_PRIVATE_OPTION ("consoleDebug",false,"Print solver debug output to console instead of separate file");

  ADD_PRIVATE_OPTION ("subtractExtFromIntForces",false,"subtracting external forces from internal forces");

  // Read the option files if any
  readOptionFiles ("fedem_solver");

  // Initialize the fedem_solver
  int ierr = 0;
  if (cfsi && stateData && ndat > 0)
#ifdef _NCHAR_AFTER_CHARARG
    F90_NAME(slv_initc,SLV_INITC) (cfsi,strlen(cfsi),stateData,ndat,ierr);
#else
    F90_NAME(slv_initc,SLV_INITC) (cfsi,stateData,ndat,ierr,strlen(cfsi));
#endif
  else if (cfsi)
#ifdef _NCHAR_AFTER_CHARARG
    F90_NAME(slv_initb,SLV_INITB) (cfsi,strlen(cfsi),ierr);
#else
    F90_NAME(slv_initb,SLV_INITB) (cfsi,ierr,strlen(cfsi));
#endif
  else if (stateData && ndat > 0)
    F90_NAME(slv_initd,SLV_INITD) (stateData,ndat,ierr);
  else if (extInput && nxinp)
    F90_NAME(slv_inita,SLV_INITA) (extInput,*nxinp,ierr);
  else
    F90_NAME(slv_init0,SLV_INIT0) (ierr);

  if (gagesData && ngda > 0 && ierr == 0)
    F90_NAME(slv_initgages,SLV_INITGAGES) (2,const_cast<double*>(gagesData),
                                           ngda,ierr);

  if (ierr)
  {
    releaseGlobalHeapObjects();
    resetCWD();
    std::cerr <<" *** Solver initialization failure ("<< ierr <<")"<< std::endl;
    return ierr > 0 ? -100-ierr : ierr;
  }

  // Initialization finished, check if we are doing a restart (iop = -3)
  double currentTime = F90_NAME(slv_gettime,SLV_GETTIME) (0,ierr);
  double startTime = F90_NAME(slv_gettime,SLV_GETTIME) (2,ierr);
  iop = currentTime > startTime ? -3 : 0;
  if (!cfsi || (stateData && ndat > 0))
    return 0;

  // In the first run (stateData is NULL) and with model input in the string,
  // return the length of the internal state vector to restart from
  return F90_NAME(slv_statesize,SLV_STATESIZE) (0);
}


DLLexport(int) restartFromState (const double* stateData, const int ndat,
                                 const int writeToRDB)
{
  int ierr = checkState("restartFromState",true);
  if (ierr < 0) return ierr;

  F90_NAME(slv_restart,SLV_RESTART) (stateData,ndat,writeToRDB,ierr);
  if (ierr < 0)
    releaseGlobalHeapObjects();
  else
    iop = -3; // Flag restart

  return ierr;
}


DLLexport(bool) solveNext (int* ierr)
{
  if ((*ierr = checkState("solveNext",true)) < 0)
    return false;

  FiDeviceFunctionFactory::instance()->updateExtFuncFromFile();

  // Solve for next time increment
  int done = 0;
  F90_NAME(slv_next,SLV_NEXT) (iop,0,done,*ierr);
  if (*ierr < 0) releaseGlobalHeapObjects();

  return done == 0;
}


DLLexport(bool) startStep (int* ierr)
{
  if (iop == 0)
    iop = -1;
  else if ((*ierr = checkState("startStep",true)) < 0)
    return false;

  FiDeviceFunctionFactory::instance()->updateExtFuncFromFile();

  // Build the first coefficient matrix of current step and return
  int done = 0;
  F90_NAME(slv_next,SLV_NEXT) (iop,0,done,*ierr);
  if (*ierr < 0) releaseGlobalHeapObjects();

  return done == 0;
}


DLLexport(bool) solveIteration (int* ierr, bool all)
{
  if ((*ierr = checkState("solveIteration")) < 0)
    return false;
  else if (iop < 1)
  {
    std::cerr <<" *** solveIteration: Step not started, iop = "<< iop
              << std::endl;
    *ierr = -999;
    return false;
  }
  else if (all && iop > 2 && iop < 5)
    iop += 12; // Solve all iterations, until convergence

  // Solve current iteration, and assemble the next one
  int done = 0;
  F90_NAME(slv_next,SLV_NEXT) (iop,0,done,*ierr);
  if (*ierr < 0) releaseGlobalHeapObjects();

  return all ? done == 0 : (done == 0 && iop > 0);
}


DLLexport(bool) solveEigenModes (const int nmod,
                                 double* eval, double* evec,
                                 bool dofOrder, int useLaPack, int *ierr)
{
  if ((*ierr = checkState("solveEigenModes") < 0))
    return false;

  F90_NAME(slv_modes,SLV_MODES) (nmod,eval,evec,dofOrder,useLaPack,*ierr);
  if (*ierr < 0) releaseGlobalHeapObjects();

  return *ierr >= 0;
}


DLLexport(bool) solveInverse (const double* x,
                              const int* xeqs, const int* feqs,
                              const int ndis, const int nfrs, int* ierr)
{
  if ((*ierr = checkState("solveInverse",true)) < 0)
    return false;

  int done = 0;
  F90_NAME(slv_inverse,SLV_INVERSE) (x,xeqs,feqs,ndis,nfrs,done,*ierr);
  if (*ierr < 0) releaseGlobalHeapObjects();

  return done == 0 && *ierr == 0;
}


DLLexport(int) solverDone (bool removeSingletons)
{
  int ierr = iop == -99 ? 0 : checkState("solverDone");
  if (ierr < 0) return ierr;

  // Clean up when finished
  F90_NAME(slv_done,SLV_DONE) (ierr);
  releaseGlobalHeapObjects(removeSingletons || ierr);
  resetCWD();

  return ierr;
}


DLLexport(void) solverClose ()
{
  if (iop <= -99)
    releaseGlobalHeapObjects(true);
}


DLLexport(double) getTime (int tFlag, int* ierr)
{
  return F90_NAME(slv_gettime,SLV_GETTIME) (tFlag,*ierr);
}


DLLexport(bool) setTime (double nextTime)
{
  int ierr = 0;
  F90_NAME(slv_settime,SLV_SETTIME) (nextTime,ierr);
  return ierr >= 0;
}


DLLexport(double) evalFunc (int uid, const char* tag, double x, int* ierr)
{
  int lerr = 0;
  int& err = ierr ? *ierr : lerr;

  if (tag)
    return F90_NAME(slv_getfunction,SLV_GETFUNCTION) (err,uid,tag,strlen(tag));

  return F90_NAME(slv_getfunc,SLV_GETFUNC) (uid,x,err);
}


DLLexport(int) getFuncId (const char* tag)
{
  return F90_NAME(slv_getfuncid,SLV_GETFUNCID) (tag,strlen(tag));
}


DLLexport(int) getEquations (int bid, int* meqn)
{
  if (checkState("getEquations") < 0)
    return 0;

  return F90_NAME(slv_geteqn,SLV_GETEQN) (bid,meqn);
}


DLLexport(int) getStateVar (int bid, double* var)
{
  if (checkState("getStateVar") < 0)
    return 0;

  return F90_NAME(slv_getvar,SLV_GETVAR) (bid,var);
}


static bool systemMatrix (int iMat, double* Smat)
{
  int ierr = checkState("systemMatrix");
  if (ierr < 0) return false;

  F90_NAME(slv_getsmat,SLV_GETSMAT) (Smat,iMat,ierr);
  return ierr == 0;
}

DLLexport(bool) getSystemMatrix (double* Smat, int iMat)
{
  return systemMatrix (iMat,Smat);
}

DLLexport(bool) getNewtonMatrix (double* Nmat)
{
  return systemMatrix (0,Nmat);
}

DLLexport(bool) getStiffnessMatrix (double* Kmat)
{
  return systemMatrix (1,Kmat);
}

DLLexport(bool) getMassMatrix (double* Mmat)
{
  return systemMatrix (2,Mmat);
}

DLLexport(bool) getDampingMatrix (double* Cmat)
{
  return systemMatrix (3,Cmat);
}


DLLexport(bool) getElementStiffnessMatrix (double* Kmat, int bid)
{
  int ierr = checkState("getElementStiffnessMatrix");
  if (ierr < 0) return false;

  F90_NAME(slv_getemat,SLV_GETEMAT) (Kmat,bid,1,ierr);
  return ierr == 0;
}


static bool rhsVector (int iVec, double* Rvec)
{
  int ierr = checkState("rhsVector");
  if (ierr < 0) return false;

  F90_NAME(slv_rhsvec,SLV_RHSVEC) (Rvec,iVec,ierr);
  return ierr == 0;
}


DLLexport(bool) getRhsVector (double* Rvec, int iVec)
{
  return rhsVector (iVec < 5 ? 0 : iVec, Rvec);
}

DLLexport(bool) setRhsVector (const double* Rvec)
{
  return rhsVector (3,const_cast<double*>(Rvec));
}

DLLexport(bool) addRhsVector (const double* Rvec)
{
  return rhsVector (4,const_cast<double*>(Rvec));
}


DLLexport(bool) setExtFunc (int funcId, double value)
{
  return FiDeviceFunctionFactory::instance()->setValue (-funcId,0.0,value);
}


DLLexport(bool) solveWindow (const int nStep, int nInc, int nIn, int nOut,
                             const int* fId, const double* times,
                             const double* inputs, double* outputs,
                             int nDat, double* stateData, int* ierr)
{
  if ((*ierr = checkState("solveWindow",true)) < 0)
    return false;

  if (!times) nInc = 0;
  if (!inputs) nIn = 0;
  if (!outputs) nOut = 0;
  if (!stateData) nDat = 0;

  if (nInc > 0 && nInc < nStep)
  {
    *ierr = nInc - nStep;
    std::cerr <<" *** Too short times array specified, "<< nInc <<" < "<< nStep
              << std::endl;
    return false;
  }

  // Loop over the time steps of this time window
  int i, j, done = 0;
  for (i = 0; i < nStep && done == 0 && *ierr == 0; i++)
  {
    if (nInc > 0)
    {
      // Explicit time steps are specified
      F90_NAME(slv_settime,SLV_SETTIME) (times[i],*ierr);
      if (*ierr > 0) *ierr = 0; // Don't abort on warnings
    }

    if (nIn > 0)
    {
      // Update external function values (typically physical sensor readings)
      for (j = 0; j < nIn; j++)
        if (!setExtFunc (1+j,*(inputs+j)))
          --(*ierr);
      inputs += nIn;
    }

    // Invoke the solver to advance the time one step forward
    int finalStep = i+1 < nStep ? 0 : 1;
    if (*ierr == 0)
      F90_NAME(slv_next,SLV_NEXT) (iop,finalStep,done,*ierr);

    if (nOut > 0)
    {
      // Extract the output values
      for (j = 0; j < nOut && *ierr == 0; j++)
        outputs[j] = F90_NAME(slv_getfunc,SLV_GETFUNC) (fId[j],-1.0,*ierr);
      outputs += nOut;
    }
  }

  if (*ierr == 0 && nDat > 0)
  {
    // Save current state
    F90_NAME(slv_savestate,SLV_SAVESTATE) (0,stateData,nDat,*ierr);
    if (*ierr > 0) *ierr = 0; // State array too long - ignore
  }

  if (*ierr < 0) releaseGlobalHeapObjects();

  return done == 0;
}


DLLexport(int) haveResults ()
{
  return F90_NAME(slv_haveresults,SLV_HAVERESULTS) ();
}


DLLexport(int) getSystemSize (bool dofs)
{
  return F90_NAME(slv_syssize,SLV_SYSSIZE) (dofs ? 1 : 0);
}


DLLexport(int) getStateSize ()
{
  return F90_NAME(slv_statesize,SLV_STATESIZE) (0);
}


DLLexport(int) getTransformationStateSize ()
{
  return F90_NAME(slv_statesize,SLV_STATESIZE) (1);
}


DLLexport(int) getGagesSize ()
{
  return F90_NAME(slv_gagessize,SLV_GAGESSIZE) ();
}


DLLexport(int) getPartDeformationStateSize (int bid)
{
  return F90_NAME(slv_partsize,SLV_PARTSIZE) (1,bid);
}


DLLexport(int) getPartStressStateSize (int bid)
{
  return F90_NAME(slv_partsize,SLV_PARTSIZE) (2,bid);
}


DLLexport(bool) saveState (double* stateData, const int ndat)
{
  int ierr = checkState("saveState");
  if (ierr < 0) return false;

  F90_NAME(slv_savestate,SLV_SAVESTATE) (0,stateData,ndat,ierr);
  return ierr >= 0;
}


DLLexport(bool) saveTransformationState (double* stateData, const int ndat)
{
  int ierr = checkState("saveTransformationState");
  if (ierr < 0) return false;

  F90_NAME(slv_savestate,SLV_SAVESTATE) (1,stateData,ndat,ierr);
  return ierr >= 0;
}


DLLexport(bool) saveGages (double* gagesData, const int ngda)
{
  int ierr = checkState("saveGages");
  if (ierr < 0) return false;

  F90_NAME(slv_initgages,SLV_INITGAGES) (1,gagesData,ngda,ierr);
  return ierr >= 0;
}


DLLexport(bool) savePartDeformationState (int bid,
                                          double* stateData, const int ndat)
{
  int ierr = checkState("savePartDeformationState");
  if (ierr < 0) return false;

  F90_NAME(slv_savepart,SLV_SAVEPART) (1,bid,stateData,ndat,ierr);
  return ierr >= 0;
}


DLLexport(bool) savePartStressState (int bid,
                                     double* stateData, const int ndat)
{
  int ierr = checkState("savePartStressState");
  if (ierr < 0) return false;

  F90_NAME(slv_savepart,SLV_SAVEPART) (2,bid,stateData,ndat,ierr);
  return ierr >= 0;
}


DLLexport(bool) getStrainsFromDisp (const double* disp, const int* gageIDs,
                                    double* eps, const int ndof, const int ng)
{
  int ierr = 0;

  F90_NAME(slv_straindisp,SLV_STRAINDISP) (disp,gageIDs,eps,ndof,ng,ierr);
  return ierr >= 0;
}


DLLexport(bool) getBeamForcesFromDisp (const double* disp, const int* beamIDs,
                                       double* forces, const int ndof,
                                       const int nBeam)
{
  int ierr = 0;

  F90_NAME(slv_beamforces,SLV_BEAMFORCES) (disp,beamIDs,forces,ndof,nBeam,ierr);
  return ierr >= 0;
}


DLLexport(bool) getRelDisp (const double* disp, const int* Ids, double* relDis,
                            const int ndof, const int nId)
{
  int ierr = 0;

  F90_NAME(slv_reldistance,SLV_RELDISTANCE) (disp,Ids,relDis,ndof,nId,ierr);
  return ierr >= 0;
}


DLLexport(bool) getRespVars (const double* disp, const int* Ids, double* Var,
                             const int ndof, const int nId)
{
  int ierr = 0;

  F90_NAME(slv_responsevars,SLV_RESPONSEVARS) (disp,Ids,Var,ndof,nId,ierr);
  return ierr >= 0;
}


DLLexport(bool) getJointSprCoeff (double* sprCoeff, int bid)
{
  int ierr = checkState("getJointSprCoeff");
  if (ierr < 0) return false;

  F90_NAME(slv_jointspring,SLV_JOINTSPRING) (sprCoeff,bid,ierr);
  return ierr >= 0;
}


DLLexport(const char*) getFileName (const char* fileOpt, char* fileName, int nc)
{
  fileName[0] = '\0';

  std::string myName;
  FFaCmdLineArg::instance()->getValue (fileOpt,myName);
  if (myName.empty()) return NULL;

  strncat(fileName,myName.c_str(),nc-1);
  return static_cast<int>(myName.size()) < nc ? NULL : fileName;
}

//! \endcond
