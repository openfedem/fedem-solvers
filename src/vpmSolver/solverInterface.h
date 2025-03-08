/* SPDX-FileCopyrightText: 2023 SAP SE
 *
 * SPDX-License-Identifier: Apache-2.0
 *
 * This file is part of FEDEM - https://openfedem.org
 */
/*!
  \file solverInterface.h

  \brief This file defines the API for the FEDEM dynamics solver core.

  \details Each function of the dynamics solver shared library that is
  supposed to be accessible for outside applications need to have their
  declaration in this file.

  \author Knut Morten Okstad, Fedem Technology AS

  \date 20 Dec 2016
*/

#ifndef _SOLVER_INTERFACE_H
#define _SOLVER_INTERFACE_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#else
#include <stdbool.h>
#endif

  /*!
    \brief Initializes a new solver process.
    \param[in] argc Number of command-line arguments passed to the solver
    \param[in] argv List of command-line arguments passed to the solver
    \param[in] cfsi Content of the solver input file describing the model
    \param[in] stateData Complete state vector to restart simulation from
    \param[in] ndat Length of the restart state vector
    \param[in] gagesData Initial strain gage values for restart
    \param[in] ngda Length of the gages array
    \param[in] extInput Initial external function values
    \param[inout] nxinp Length of the extInput array
    \return &ge; 0 : Success
    \return &lt; 0 : An error has occured

    \details
    This function processes the input and sets up necessary data structures
    prior to the time integration loop. It also performs the license checks.
    See the Fedem R7.4 Users Guide, Appendix C.3 for a complete list of all
    command-line arguments that may be specified and their default values.

    If the argument \a cfsi is NULL, the model is assumed defined in the file
    specified via the command-line argument -fsifile instead.

    The arguments \a extInput and \a nxinp are used to transfer initial values
    of the external functions to the solver before the initial equilibrirum
    iterations is performed. The argument \a nxinp is reset to zero on exit,
    if initial equilibrium iterations is not performed, and to a negative valiue
    if the array \a extInput is too short.

    The return value is zero on success, a negative value indicates some error.
    If \a cfsi is NULL and \a stateData is NULL, the return value is the
    required minimum length of the array to be used as state vector for restart.
  */
  int solverInit(int argc, char** argv, const char* cfsi = NULL,
                 const double* stateData = NULL, const int ndat = 0,
                 const double* gagesData = NULL, const int ngda = 0,
                 const double* extInput  = NULL, int* nxinp = NULL);

  /*!
    \brief Restarts a simulation process from the provided state.
    \param[in] stateData Complete state vector to restart simulation from
    \param[in] ndat Length of the restart state vector
    \param[in] writeToRDB Flag for saving response variables (see below)
    \return &nbsp;&nbsp; 0 : Success
    \return &lt; 0 : An error occurred

    \details
    This function re-initializes the mechanism elements with data from the
    provided state array, such that the simulation can continue from that state.

    The value of the input flag \a writeToRDB is interpreted as follows:
    - 0: No results saving
    - 1: Append results to the already opened results database
    - 2: Increment the results database and write to new files
  */
  int restartFromState(const double* stateData, const int ndat,
                       const int writeToRDB = 2);

  /*!
    \brief Executes the simulation for a specified time window.
    \param[in] nStep Number of time/load steps to solve for from current state
    \param[in] nInc Number of explicit time increments
    \param[in] nIn Number of input sensor values
    \param[in] nOut Number of output sensor values
    \param[in] fOut List of user IDs identifying the output sensors in the model
    \param[in] times Explicit times to solver for (dimension nInc)
    \param[in] inputs List of input sensor values (dimension nIn*nStep)
    \param[out] outputs List of output sensor values (dimension nOut*nStep)
    \param[in] nDat Length of the restart state vector
    \param[out] stateData Complete state vector at the end of the time window
    \param[out] ierr A non-zero value indicates some error condition that
                     requires the simulation to be terminated
    \return \e true, unless the end of the simulation has been reached

    \details This function solves the problem for a time/load step window, with
    given values for the external functions (physical sensor readings), and
    extraction of results from another set of general functions in the model.
  */
  bool solveWindow(int nStep, int nInc, int nIn, int nOut, const int* fOut,
                   const double* times, const double* inputs, double* outputs,
                   int nDat, double* stateData, int* ierr);

  /*!
    \brief Returns whether current time step have results to be saved.
  */
  int haveResults();

  /*!
    \brief Returns the required length of the state vector to be is used when
    restarting a simulation from an in-core array.
  */
  int getStateSize();

  /*!
    \brief Returns the required length of the state vector which contains
    translation and rotation information for all triads and superelements.
  */
  int getTransformationStateSize();

  /*!
    \brief Returns the required length of the state vector which contains
    deformation information for the specified FE part.
    \param[in] bid Base ID of the FE part to consider
  */
  int getPartDeformationStateSize(int bid);

  /*!
    \brief Returns the required length of the state vector which contains
    von Mises stress information for the specified FE part.
    \param[in] bid Base ID of the FE part to consider
  */
  int getPartStressStateSize(int bid);

  /*!
    \brief Returns the required length of the initial strain gage vector
    which is used when restarting a simulation from an in-core array.
  */
  int getGagesSize();

  /*!
    \brief Saves current state to the provided array.
    \param[out] stateData Complete state vector
    \param[in] ndat Length of the state vector

    \details This function can be used after a succesful ::solveWindow call,
    or after a ::solveNext call returning \e false indicating the end of the
    simulation, to save current state to an in-core array that can be
    used to restart the simulation later on.
  */
  bool saveState(double* stateData, const int ndat);

  /*!
    \brief Saves current transformation state to the provided core array.
    \param[out] stateData Transformation state vector
    \param[in] ndat Length of the state vector
  */
  bool saveTransformationState(double* stateData, const int ndat);

  /*!
    \brief Saves the initial gage strains to the provided core array.
    \param[out] stateData Gage strain state vector
    \param[in] ndat Length of the state vector
  */
  bool saveGages(double* stateData, const int ndat);

  /*!
    \brief Gets gage strain values from given displacement field.
    \param[in] disp    Displacement vector
    \param[in] gageIDs Array of strain gage ID numbers
    \param[out] eps    Strain tensor at the strain gages
    \param[in] ndof    Length of disp (number of degrees of freedom)
    \param[in] ng      Length of gageIDs (number of gages)
  */
  bool getStrainsFromDisp(const double* disp, const int* gageIDs, double* eps,
                          const int ndof, const int ng);

  /*!
    \brief Gets internal sectional forces from given displacements.
    \param[in] disp    Displacement vector
    \param[in] beamIDs Array of beam IDs
    \param[out] forces Force vector at the sections
    \param[in] ndof    Number of degrees of freedom
    \param[in] nBeam   Number of beams
  */
  bool getBeamForcesFromDisp(const double* disp, const int* beamIDs,
                             double* forces, const int ndof, const int nBeam);

  /*!
    \brief Gets relative distances between triads from given displacements.
    \param[in]  disp   Displacement vector
    \param[in]  Ids    Array of engine IDs
    \param[out] relDis Relative distance change
    \param[in]  ndof   Number of degrees of freedom
    \param[in]  nId    Number of engine IDs
  */
  bool getRelDisp(const double* disp, const int* Ids, double* relDis,
                  const int ndof, const int nId);

  /*!
    \brief Gets response variable values from given displacements.
    \param[in]  disp   Displacement vector
    \param[in]  Ids    Array of engine IDs
    \param[out] resVar Response variable values
    \param[in]  ndof   Number of degrees of freedom
    \param[in]  nId    Number of engine IDs
  */
  bool getRespVars(const double* disp, const int* Ids, double* resVar,
                   const int ndof, const int nId);

  /*!
    \brief Saves the deformation of an FE part to the provided core array.
    \param[in] bid Base ID of the FE part to consider
    \param[out] stateData Deformation state vector
    \param[in] ndat Length of the state vector

    \details This function is used to save current deformation state
    for the specified FE part to the provided core array after a successful
    ::solveNext call, for access by external visualisation modules.
  */
  bool savePartDeformationState(int bid, double* stateData, const int ndat);

  /*!
    \brief Saves the stress state of an FE part to the provided core array.
    \param[in] bid Base ID of the FE part to consider
    \param[out] stateData Stress state vector
    \param[in] ndat Length of the state vector

    \details This function is used to save current von Mises stress state
    for the specified FE part to the provided core array after a successful
    ::solveNext call, for access by external visualisation modules.
  */
  bool savePartStressState(int bid, double* stateData, const int ndat);

  /*!
    \brief Advances the solution one time/load step forward.
    \param[out] ierr Equal to zero on a successful computation, a non-zero
    value indicates some error that requires to simulation to terminate
    \return \e true, unless current time/load step failed to converge
    or the end time of the simulation has been reached
  */
  bool solveNext(int* ierr);

  /*!
    \brief Starts a new time/load step.
    \param[out] ierr Equal to zero on a successful computation, a non-zero
    value indicates some error that requires to simulation to terminate
    \return \e true, unless the simulation has to stop, either because of an
    error condition, or because the end time of the simulation has been reached

    \details This function starts a new time (or load) step, by calculating the
    predicted response and the coefficient matrix and right-hand-side vector of
    the first nonlinear iteration. It has to be followed up by a series of
    ::solveIteration calls in order to continue the simulation,
    but the linear equation system can be manipulated in between.
  */
  bool startStep(int* ierr);

  /*!
    \brief Solves current linearized equation system and updates the state.
    \param[out] ierr Equal to zero on a successful computation, a non-zero
    value indicates some error that requires to simulation to terminate
    \param[in] all If \e true, repeat until convergence
    \return \e true, unless the simulation has to stop, either because of an
    error condition, or (if \a all is \e false) the current time/load step has
    reached convergence

    \details This function solves the current linearized equation system and
    updates all state variables in the model. Then it assembles the linearized
    system of equations for next iteration, unless convergence has been reached.
    if \a all is \e true, this is repeated until convergence is achieved.
  */
  bool solveIteration(int* ierr, bool all = false);

  /*!
    \brief Solves the eigenvalue problem at current state.
    \param[in] nmod Number of eigenmodes to calculate
    \param[out] eval Computed eigenvalues
    \param[out] evec Computed eigenvectors
    \param[in] dofOrder If \e true, return eigenvectors in DOF-order.
    Otherwise, they are returned in equation order
    \param[in] useLaPack Flag usage of LAPACK eigensolvers (1=DSYGVX, 2=DGGEVX)
    \param[out] ierr Equal to zero on a successful computation, a negative value
    indicates memory allocation error or similar issues, positive value indicate
    that the eigenvalue solver failed, but the execution may continue
    \return \e true, unless the simulation has to stop due to an error condition
  */
  bool solveEigenModes(const int nmod, double* eval, double* evec,
                       bool dofOrder, int useLaPack, int* ierr);

  /*!
    \brief Solves the inverse problem at current time/load step.
    \param[in] x Specified displacement values at a set of degrees of freedom
    \param[in] xeqs Equation numbers (1-based) for the specified displacements
    \param[in] feqs Equation numbers (1-based) for the DOFs with unknown
                    external forces
    \param[in] ndis Number of DOFs with specified displacements
    \param[in] nfrs Number of DOFs with unknown external forces
    \param[out] ierr Equal to zero on a successful computation, a non-zero
    value indicates some error that requires the simulation to terminate
    \return \e true, unless the simulation has to stop due to some
    error condition, or the end of the simulation has been reached
  */
  bool solveInverse(const double* x,
                    const int* xeqs, const int* feqs,
                    int ndis, int nfrs, int* ierr);

  /*!
    \brief Closes down current model and cleans up core memory and on disk.
    \param[in] removeSingletons If \e false, singleton objects are not deleted
    \return Zero on success, a non-zero value indicates some error

    \details This function should be invoked only when the time (or load) step
    loop is finished. If the process is going to be restarted later on the same
    model, the \a removeSingletons argument should be \e true, to ensure that
    objects allocated on the heap only on program startup are not deleted.
  */
  int solverDone(bool removeSingletons = true);

  /*!
    \brief Cleans up singleton objects in core memory.
    \details This function needs to be invoked only if the solveDone() call
    was performed with \a removeSingletons = \e false, before finishing.
  */
  void solverClose();

  /*!
    \brief Assigns a sensor value from a physical twin to the simulation model.
    \param[in] funcId ID of the external function to receive the sensor value
    \param[in] value The actual sensor value to be assigned
    \return \e false if \a funcId is out of range, otherwise \e true

    \details This function may be used prior to the ::solveNext call,
    to assign a sensor value from a physical twin to the specified actuator
    or load in the digital twin model.
  */
  bool setExtFunc(int funcId, double value);

  /*!
    \brief Returns the physical time of a time step.
    \param[in] tFlag Flag indicating which time step to return the time for
    - 0 : Return for current time step
    - 1 : Return for next time step
    - 2 : Return the start time of the simulation
    - 3 : Return the stop time of the simulation
    \param ierr Error flag, untouched unless \a tFlag = 1
  */
  double getTime(int tFlag, int* ierr);

  /*!
    \brief Sets the time increment size of the next time step.
    \param[in] nextTime Physical time to calculate next time increment size from
  */
  bool setTime(double nextTime);

  /*!
    \brief Returns the current physical time of the simulation.
  */
  inline double getCurrentTime() { int ierr; return getTime(0,&ierr); }

  /*!
    \brief Returns the physical time of the next step of the simulation.
    \param ierr Untouched on a successful call,
    decremented if the time step size in the model is defined via
    a general function that could not be evaluated.
  */
  inline double getNextTime(int* ierr) { return getTime(1,ierr); }

  /*!
    \brief Evaluates a specified general function in the model.
    \param[in] uid User ID of the general function to be evaluated
    \param[in] tag Alternative function identification
    \param[in] x If zero or positive, this is the function argument value.
    If negative, the value is ignored and the function is evaluated for the
    current state of the sensor object that is defined as the function argument.
    If \a tag is specified, the value is ignored (same behavior as if negative).
    \param[out] ierr Equal to zero on a successful call, a non-zero value
    indicates that the function could not be evaluated
    \return The evaluated function value
  */
  double evalFunc(int uid, const char* tag = NULL,
                  double x = -1.0, int* ierr = NULL);

  /*!
    \brief Returns the current value of the specified general function.
    \param[in] uid User ID of the general function to be evaluated
    \param[out] ie Equal to zero on a successful call, a non-zero value
    indicates that the function could not be evaluated
    \return The evaluated function value
  */
  double getFunc(int uid, int* ie = NULL) { return evalFunc(uid,NULL,-1.0,ie); }

  /*!
    \brief Returns the user ID of a tagged function.
  */
  int getFuncId(const char* tag);

  /*!
    \brief Returns the equation numbers for the DOFs of a specified object.
    \param[in] bid Base ID of the mechanism object to get equations for
    \param[out] meqn List of equation numbers (1-based)
    \return &nbsp;&nbsp; 0 : Non-existing object or object with no DOFs
    \return &gt; 0 : Number of DOFs for the specified object (length of meqn)
  */
  int getEquations(int bid, int* meqn);

  /*!
    \brief Returns the current state variables for a specified object.
    \param[in] bid Base ID of the mechanism object to get statevariables for
    \param[out] var List of variable values
    \return &nbsp;&nbsp; 0 : Non-existing object or object with no DOFs
    \return &gt; 0 : Number of state variables for the specified object

    \details This function is provided to facilitate unit testing where we
    want to assert on certain state variables. For triads, only the updated
    positions are returned (no rotation quantities).
  */
  int getStateVar(int bid, double* var);

  /*!
    \brief Returns the dimension (number of equations) of the system.
    \param[in] dofs If \e true, return the total number of DOFs in the system
  */
  int getSystemSize(bool dofs = false);

  /*!
    \brief Returns current content of the system Newton matrix.
    \param[out] Smat The system matrix in full format
    \param[in] iMat System matrix type flag (0, 1, 2 or 3)
  */
  bool getSystemMatrix(double* Smat, int iMat = 0);

  /*!
    \brief Returns current content of the system Newton matrix.
    \param[out] Nmat The system Newton matrix in full format
  */
  bool getNewtonMatrix(double* Nmat);

  /*!
    \brief Returns current content of the system stiffness matrix.
    \param[out] Kmat The system stiffness matrix in full format
  */
  bool getStiffnessMatrix(double* Kmat);

  /*!
    \brief Returns current content of the system mass matrix.
    \param[out] Mmat The system mass matrix in full format
  */
  bool getMassMatrix(double* Mmat);

  /*!
    \brief Returns current content of the system damping matrix.
    \param[out] Cmat The system damping matrix in full format
  */
  bool getDampingMatrix(double* Cmat);

  /*!
    \brief Returns the content of an element stiffness matrix.
    \param[out] Kmat The element stiffness matrix
    \param[in] bid Base ID of the superelement (beam, reduced FE part or
    generic part) to return the stiffness matrix for
  */
  bool getElementStiffnessMatrix(double* Kmat, int bid);

  /*!
    \brief Returns current content of the system right-hand-side vector.
    \param[out] Rvec The system right-hand-side vector
    \param[in] iVec System vector type flag (0=residual or 5=external)
  */
  bool getRhsVector(double* Rvec, int iVec = 0);

  /*!
    \brief Replaces current content of the system right-hand-side vector.
    \param[in] Rvec The new system right-hand-side vector
  */
  bool setRhsVector(const double* Rvec);

  /*!
    \brief Updates the content of the system right-hand-side vector.
    \param[in] Rvec Vector to be added to the current right-hand-side vector
  */
  bool addRhsVector(const double* Rvec);

  /*!
    \brief Gets spring stiffness coefficients for a joint.
    \param[out] sprCoeff Spring stiffness coefficients
    \param[in] bid Base ID of joint
  */
  bool getJointSprCoeff(double* sprCoeff, int bid);

  /*!
    \brief Gets a file name option from the command-line arguments.
    \param[in] fileOpt Command-line option specifying a file name
    \param[out] fileName Value of the specified command-line option
    \param[in] nc Max length of fileName
  */
  const char* getFileName(const char* fileOpt, char* fileName, int nc);

#ifdef __cplusplus
}
#endif

#endif
