!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!===============================================================================
!> @file solverDriver.f90
!>
!> @brief Global driver routines for the FEDEM Dynamics Solver.
!>
!> @details All subroutines and functions of this file are global wrappers of
!> the subroutines located in the #solvermodule. They are defined in a way to
!> facilitate invocation from C or C++, only with arguments of primitive data
!> types that can be casted to equivalent data types in C/C++.
!>
!> All the subroutines are used by the C++ API defined in solverInterface.C,
!> however, no caller graphs are provided for them since Doxygen seems not to
!> be able to generate cross-language call graphs.
!>
!> @author Knut Morten Okstad, Fedem Technology AS
!>
!> @date 30 Nov 2016
!===============================================================================

!===============================================================================
!> @brief Initializes a new simulation process.
!> @callgraph
subroutine slv_init0 (ierr)
  use solverModule, only : initialize
  implicit none
  integer, intent(out) :: ierr !< Error flag
  call initialize (ierr=ierr)
end subroutine slv_init0

!===============================================================================
!> @brief Initializes a new simulation process.
!> @callgraph
subroutine slv_initA (xinp,nxinp,ierr)
  use kindModule  , only : dp
  use solverModule, only : initialize
  implicit none
  real(dp), intent(in)    :: xinp(*) !< Initial external function values
  integer , intent(inout) :: nxinp   !< Length of the xinp array
  integer , intent(out)   :: ierr    !< Error flag
  call initialize (xinp=xinp,nxinp=nxinp,ierr=ierr)
end subroutine slv_initA

!===============================================================================
!> @brief Initializes a new simulation process.
!> @callgraph
subroutine slv_initB (chfsi,ierr)
  use solverModule, only : initialize
  implicit none
  character(len=*), intent(in)  :: chfsi !< Content of the solver input file
  integer         , intent(out) :: ierr  !< Error flag
  call initialize (chfsi,ierr=ierr)
end subroutine slv_initB

!===============================================================================
!> @brief Initializes a new simulation process and restarts from provided state.
!> @callgraph
subroutine slv_initC (chfsi,data,ndat,ierr)
  use kindModule  , only : dp
  use solverModule, only : initialize
  implicit none
  character(len=*), intent(in)  :: chfsi   !< Content of the solver input file
  real(dp)        , intent(in)  :: data(*) !< Array of state variables
  integer         , intent(in)  :: ndat    !< Length of the data array
  integer         , intent(out) :: ierr    !< Error flag
  call initialize (chfsi,data,ndat,ierr=ierr)
end subroutine slv_initC

!===============================================================================
!> @brief Initializes a new simulation process and restarts from provided state.
!> @callgraph
subroutine slv_initD (data,ndat,ierr)
  use kindModule  , only : dp
  use solverModule, only : initialize
  implicit none
  real(dp), intent(in)  :: data(*) !< Array of state variables
  integer , intent(in)  :: ndat    !< Length of the data array
  integer , intent(out) :: ierr    !< Error flag
  call initialize (stateData=data,ndat=ndat,ierr=ierr)
end subroutine slv_initD

!===============================================================================
!> @brief Restarts a running solver process from the provided state.
!> @callgraph
subroutine slv_restart (data,ndat,writeToRDB,ierr)
  use kindModule  , only : dp
  use solverModule, only : softRestart
  implicit none
  real(dp), intent(in)  :: data(*)    !< Array of state variables
  integer , intent(in)  :: ndat       !< Length of the data array
  integer , intent(in)  :: writeToRDB !< Flag for controlling frs-file output
  integer , intent(out) :: ierr       !< Error flag
  call softRestart (data,ndat,writeToRDB,ierr)
end subroutine slv_restart

!===============================================================================
!> @brief Solves the next time increment by Newton-Raphson iterations.
!> @param iop Operation flag telling what to do, see solvermodule::solvestep
!> @param[in] final If .true., this is assumed to be the final time step
!> @param[out] done If .true., the simulation is finished or failed
!> @param[out] ierr Error flag
!> @callgraph
subroutine slv_next (iop,final,done,ierr)
  use freqResponseModule, only : solveFreqDomain
  use solverModule      , only : solveRampUp, solveStep
  implicit none
  integer, intent(inout) :: iop
  logical, intent(in)    :: final
  logical, intent(out)   :: done
  integer, intent(out)   :: ierr
  if (iop == 0) then
     call solveRampUp (iop,done,ierr)
     if (ierr < 0) return
  end if
  if (iop < 1) then
     call solveFreqDomain (ierr)
     if (ierr < 0) return
  end if
  call solveStep (iop,final,done,ierr)
end subroutine slv_next

!===============================================================================
!> @brief Solves the eigenvalue system at current configuration.
!> @param[in] nMode Number of eigenmodes to compute
!> @param[out] eval Computed eigenvalues
!> @param[out] evec Computed eigenvectors
!> @param[in] dofOrder If .true., return eigenvectors in DOF-order
!> @param[in] useLaPack Flag usage of LAPACK dense matrix eigensolver
!> @param[out] ierr Error flag
!> @callgraph
subroutine slv_modes (nMode,eval,evec,dofOrder,useLaPack,ierr)
  use kindModule  , only : dp
  use solverModule, only : solveModes
  implicit none
  integer , intent(in)  :: nMode, useLaPack
  real(dp), intent(out) :: eval(*), evec(*)
  logical , intent(in)  :: dofOrder
  integer , intent(out) :: ierr
  call solveModes (nMode,eval,evec,dofOrder,useLaPack,ierr)
end subroutine slv_modes

!===============================================================================
!> @brief Solves the inverse problem at the current state.
!> @param[in] x Array with response data at certain DOFs
!> @param[in] xeqs Array of equation numbers of the DOFs with @b x values
!> @param[in] feqs Array of equation numbers of the DOFs for force calculation
!> @param[in] ndis Length of array x and xeqs
!> @param[in] nfrc Length of array feqs
!> @param[out] done If .true., the simulation is finished or failed
!> @param[out] ierr Error flag
!> @callgraph
subroutine slv_inverse (x,xeqs,feqs,ndis,nfrc,done,ierr)
  use kindModule   , only : dp
  use inverseModule, only : solveInverse
  implicit none
  real(dp), intent(in)  :: x(*)
  integer , intent(in)  :: xeqs(*), feqs(*)
  integer , intent(in)  :: ndis, nfrc
  logical , intent(out) :: done
  integer , intent(out) :: ierr
  call solveInverse (x(1:ndis),xeqs(1:ndis),feqs(1:nfrc),done,ierr)
end subroutine slv_inverse

!===============================================================================
!> @brief Terminates the running solver process by deallocating all memory, etc.
!> @callgraph
subroutine slv_done (ierr)
  use solverModule      , only : finalize
  use inverseModule     , only : doneInverse
  use freqResponseModule, only : deallocateFreq
  implicit none
  integer, intent(out) :: ierr !< Error flag
  call doneInverse ()
  call deallocateFreq ()
  call finalize (ierr)
end subroutine slv_done

!===============================================================================
!> @brief Returns the physical time of a specified simulation state.
!> @callgraph
function slv_gettime (tflag,ierr) result(time)
  use kindModule  , only : dp
  use solverModule, only : getTime
  implicit none
  integer, intent(in)    :: tflag !< Indicates which step to return the time for
  integer, intent(inout) :: ierr  !< Error flag
  real(dp) :: time
  call getTime (time,tflag,ierr)
end function slv_gettime

!===============================================================================
!> @brief Sets the time increment size of the next time step.
!> @callgraph
subroutine slv_settime (nextTime,ierr)
  use kindModule  , only : dp
  use solverModule, only : setNewTime
  implicit none
  real(dp), intent(in) :: nextTime !< Time to calculate next time increment from
  integer, intent(out) :: ierr     !< Error flag
  call setNewTime (nextTime,ierr)
end subroutine slv_settime

!===============================================================================
!> @brief Evaluates the specified function.
!> @callgraph
function slv_getfunc (uid,x,ierr) result(value)
  use kindModule  , only : dp
  use solverModule, only : getEngine
  implicit none
  integer , intent(in)    :: uid  !< User ID of the function to evaluate
  real(dp), intent(in)    :: x    !< The function argument value
  integer , intent(inout) :: ierr !< Error flag
  real(dp) :: value
  if (x < 0.0_dp) then
     call getEngine (value,ierr,uid)
  else
     call getEngine (value,ierr,uid,X=x)
  end if
end function slv_getfunc

!===============================================================================
!> @brief Evaluates the (first) function identified by the given tag.
!> @callgraph
function slv_getfunction (ierr,uid,tag) result(value)
  use kindModule  , only : dp
  use solverModule, only : getEngine
  implicit none
  integer         , intent(in)    :: uid  !< User ID of the function to evaluate
  character(len=*), intent(in)    :: tag  !< Tag identifying the function
  integer         , intent(inout) :: ierr !< Error flag
  real(dp) :: value
  call getEngine (value,ierr,uid,tag)
end function slv_getfunction

!===============================================================================
!> @brief Returns the user ID of (first) function identified by the given tag.
!> @callgraph
function slv_getfuncid (tag) result(uid)
  use solverModule, only : getEngineId
  implicit none
  character(len=*), intent(in) :: tag !< Tag identifying the function
  integer :: uid
  uid = getEngineId(tag)
end function slv_getfuncid

!===============================================================================
!> @brief Returns a system matrix in full format.
!> @callgraph
subroutine slv_getSmat (Nmat,iopM,ierr)
  use kindModule  , only : dp
  use solverModule, only : getSystemMatrix
  implicit none
  real(dp), intent(out) :: Nmat(*) !< The system matrix
  integer , intent(in)  :: iopM    !< Flag indicating which matrix to return
  integer , intent(out) :: ierr    !< Error flag
  call getSystemMatrix (Nmat,iopM,ierr)
end subroutine slv_getSmat

!===============================================================================
!> @brief Returns an element matrix.
!> @callgraph
subroutine slv_getEmat (Emat,bid,iopM,ierr)
  use kindModule  , only : dp
  use solverModule, only : getElementMatrix
  implicit none
  real(dp), intent(out) :: Emat(*) !< The element matrix
  integer , intent(in)  :: bid     !< Base ID of the element to consider
  integer , intent(in)  :: iopM    !< Flag indicating which matrix to return
  integer , intent(out) :: ierr    !< Error Flag
  call getElementMatrix (Emat,bid,iopM,ierr)
end subroutine slv_getEmat

!===============================================================================
!> @brief Sets, updates or retrieves a system right-hand-side vector.
!> @callgraph
subroutine slv_RHSvec (Rvec,iopV,ierr)
  use kindModule  , only : dp
  use solverModule, only : getRhsVector, setRhsVector
  implicit none
  real(dp), intent(inout) :: Rvec(*) !< The right-hand-side vector
  integer , intent(in)    :: iopV    !< Operation flag telling what to do
  integer , intent(out)   :: ierr    !< Error Flag
  select case (iopV)
  case (1); call setRhsVector (Rvec,.false.,.false.,ierr)
  case (2); call setRhsVector (Rvec,.true.,.false.,ierr)
  case (3); call setRhsVector (Rvec,.false.,.true.,ierr)
  case (4); call setRhsVector (Rvec,.true.,.true.,ierr)
  case default; call getRhsVector (Rvec,iopV-4,ierr)
  end select
end subroutine slv_RHSvec

!===============================================================================
!> @brief Calculates the strains at strain gages for a given displacement field.
!> @callgraph
subroutine slv_StrainDisp (disp,gageIds,eps,ndof,ng,ierr)
  use kindModule  , only : dp
  use solverModule, only : computeGageStrains
  implicit none
  real(dp), intent(in)  :: disp(*)    !< Displacement vector
  integer , intent(in)  :: gageIds(*) !< Base IDs of strain gages
  real(dp), intent(out) :: eps(*)     !< Strain tensor values at gage positions
  integer , intent(in)  :: ndof       !< Length of displacement vector
  integer , intent(in)  :: ng         !< Length of gageIds array
  integer , intent(out) :: ierr       !< Error flag
  call computeGageStrains (disp,gageIds,ndof,ng,eps,ierr)
end subroutine slv_StrainDisp

!===============================================================================
!> @brief Calculates beam section forces for given displacement vector.
!> @callgraph
subroutine slv_BeamForces (disp,beamIds,forces,ndof,nBeams,ierr)
  use kindModule  , only : dp
  use solverModule, only : computeBeamForces
  implicit none
  real(dp), intent(in)  :: disp(*)    !< Displacement vector
  integer , intent(in)  :: beamIds(*) !< Base IDs of beam elements
  real(dp), intent(out) :: forces(*)  !< Sectional forces at beam ends
  integer , intent(in)  :: ndof       !< Length of displacement vector
  integer , intent(in)  :: nBeams     !< Number of beam elements
  integer , intent(out) :: ierr       !< Error flag
  call computeBeamForces (disp,beamIds,ndof,nBeams,forces,ierr)
end subroutine slv_BeamForces

!===============================================================================
!> @brief Calculates the rel. dist between triads for given displacement vector.
!> @callgraph
subroutine slv_RelDistance (disp,Ids,relDis,ndof,nIds,ierr)
  use kindModule  , only : dp
  use solverModule, only : computeRelativeDistance
  implicit none
  real(dp), intent(in)  :: disp(*)   !< Displacement vector
  integer,  intent(in)  :: Ids(*)    !< User IDs of functions to use
  real(dp), intent(out) :: relDis(*) !< Relative displacement change
  integer,  intent(in)  :: ndof      !< Length of displacement vector
  integer,  intent(in)  :: nIds      !< Length of Ids and relDis arrays
  integer,  intent(out) :: ierr      !< Error flag
  call computeRelativeDistance (disp,Ids,relDis,ndof,nIds,ierr)
end subroutine slv_RelDistance

!===============================================================================
!> @brief Calculates response variable values for given displacement vector.
!> @callgraph
subroutine slv_ResponseVars (disp,Ids,Var,ndof,nIds,ierr)
  use kindModule  , only : dp
  use solverModule, only : computeResponseVars
  implicit none
  real(dp), intent(in)  :: disp(*) !< Displacement vector
  integer,  intent(in)  :: Ids(*)  !< User IDs of functions to use
  real(dp), intent(out) :: Var(*)  !< Response variable values
  integer,  intent(in)  :: ndof    !< Length of displacement vector
  integer,  intent(in)  :: nIds    !< Length of Ids and Var arrays
  integer,  intent(out) :: ierr    !< Error flag
  call computeResponseVars (disp,Ids,Var,ndof,nIds,ierr)
end subroutine slv_ResponseVars

!===============================================================================
!> @brief Gets joint spring stiffness coefficients.
!> @callgraph
subroutine slv_JointSpring (sprCoeff,bid,ierr)
  use kindModule  , only : dp
  use solverModule, only : getJointSpringStiffness
  implicit none
  real(dp), intent(out) :: sprCoeff(*) !< Joint spring stiffness coefficients
  integer,  intent(in)  :: bid         !< Base ID of the spring
  integer,  intent(out) :: ierr        !< Error flag
  call getJointSpringStiffness (sprCoeff,bid,ierr)
end subroutine slv_JointSpring

!===============================================================================
!> @brief Returns the equation numbers affected by a given mechanism object.
!> @callgraph
function slv_geteqn (bid,meqn) result(ndof)
  use solverModule, only : objectEquations
  implicit none
  integer, intent(in)  :: bid     !< Base ID of the mechanism object to consider
  integer, intent(out) :: meqn(*) !< List of equation numbers
  integer :: ndof
  call objectEquations (bid,meqn,ndof)
end function slv_geteqn

!===============================================================================
!> @brief Returns the state variables for a given mechanism object.
!> @callgraph
function slv_getvar (bid,var) result(nvar)
  use kindModule  , only : dp
  use solverModule, only : objectStateVar
  implicit none
  integer , intent(in)  :: bid    !< Base ID of the mechanism object to consider
  real(dp), intent(out) :: var(*) !< List of state variables
  integer :: nvar
  call objectStateVar (bid,var,nvar)
end function slv_getvar

!===============================================================================
!> @brief Returns the dimension of the system matrices and vectors.
!> @callgraph
function slv_syssize (expanded) result(ndim)
  use solverModule, only : systemSize
  implicit none
  logical, intent(in) :: expanded !< If .true, return total number of DOFs
  integer :: ndim
  call systemSize (ndim,expanded)
end function slv_syssize

!===============================================================================
!> @brief Returns the size of the total state array.
!> @callgraph
function slv_statesize (posOnly) result(ndat)
  use solverModule, only : stateVectorSize
  implicit none
  logical, intent(in) :: posOnly !< If .true., consider position matrices only
  integer :: ndat
  call stateVectorSize (posOnly,ndat)
end function slv_statesize

!===============================================================================
!> @brief Returns the size of the strain gage state array.
!> @callgraph
function slv_gagessize () result(ndat)
  use solverModule, only : strainGagesSize
  implicit none
  integer :: ndat
  call strainGagesSize (ndat)
end function slv_gagessize

!===============================================================================
!> @brief Returns the size of the deformation/stress state array for a FE part.
!> @callgraph
function slv_partsize (iopS,bid) result(ndat)
  use solverModule, only : partStateVectorSize
  implicit none
  integer, intent(in) :: iopS !< Operation code telling which state to consider
  integer, intent(in) :: bid  !< Base ID of the FE part to consider
  integer :: ndat
  call partStateVectorSize (iopS,bid,ndat)
end function slv_partsize

!===============================================================================
!> @brief Saves the current state to an array in core.
!> @callgraph
subroutine slv_savestate (posOnly,data,ndat,ierr)
  use kindModule  , only : dp
  use solverModule, only : saveState
  implicit none
  logical , intent(in)  :: posOnly !< If .true., consider position matrices only
  real(dp), intent(out) :: data(*) !< Array receiving the state variables
  integer , intent(in)  :: ndat    !< Length of the state array
  integer , intent(out) :: ierr    !< Error flag
  call saveState (posOnly,data,ndat,ierr)
end subroutine slv_savestate

!===============================================================================
!> @brief Save/restores initial gage strain state to/from an array in core.
!> @callgraph
subroutine slv_initgages (iopS,data,ndat,ierr)
  use kindModule  , only : dp
  use solverModule, only : saveInitGageStrains
  implicit none
  integer , intent(in)  :: iopS    !< Operation code telling what to do
  real(dp), intent(out) :: data(*) !< Array receiving the state variables
  integer , intent(in)  :: ndat    !< Length of the state array
  integer , intent(out) :: ierr    !< Error flag
  call saveInitGageStrains (iopS,data,ndat,ierr)
end subroutine slv_initgages

!===============================================================================
!> @brief Saves current deformation state for an FE part to an array in core.
!> @callgraph
subroutine slv_savepart (iopS,bid,data,ndat,ierr)
  use kindModule  , only : dp
  use solverModule, only : savePartState
  implicit none
  integer , intent(in)  :: iopS    !< Operation code telling what to save
  integer , intent(in)  :: bid     !< Base ID of the FE part to save for
  real(dp), intent(out) :: data(*) !< Array receiving the state variables
  integer , intent(in)  :: ndat    !< Length of the state array
  integer , intent(out) :: ierr    !< Error flag
  call savePartState (iopS,bid,data,ndat,ierr)
end subroutine slv_savepart

!===============================================================================
!> @brief Checks whether current time step have results to be saved.
!> @callgraph
function slv_haveresults ()
  use solverModule, only : haveResults
  implicit none
  integer :: slv_haveresults
  slv_haveresults = haveResults()
end function slv_haveresults
