!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file systemTypeModule.f90
!> @brief System data container.

!!==============================================================================
!> @brief Module with system data containers.

module SystemTypeModule

  use KindModule         , only : dp, i8
  use FunctionTypeModule , only : EngineType, dp
  use NormTypeModule     , only : TestSetType
  use SysMatrixTypeModule, only : SysMatrixType

  implicit none

  !> @brief Data type containing parameters and system matrices of the model.
  type SystemType

     real(dp) :: time     !< Current physical time
     real(dp) :: timeStep !< Current time step size
     real(dp) :: timePrev !< Previous time step size
     real(dp) :: tStart   !< Start time of simulation
     real(dp) :: tEnd     !< Stop time of simulation
     real(dp) :: tQStatic !< Stop time of quasi-static simulation
     real(dp) :: tInc     !< Initial time increment size
     real(dp) :: maxInc   !< Maximum time increment size
     real(dp) :: minInc   !< Minimum time increment size

     !> Indicator for automatic time increments.
     !> - = 0 No, use constant or prescribed time increment sizes
     !> - = 1 Yes, based on local truncation error
     integer  :: varInc
     !> Number of steps to perform after each cut-back
     !> before increasing the time increment size again
     integer  :: nCutStp
     !> Iteration cut-back ratios
     real(dp) :: cutbck(2)

     type(EngineType), pointer :: tIncEngine !< Time step size function
     type(EngineType), pointer :: fxItEngine !< Number of iterations function

     real(dp) :: alpha_m !< Generalized alpha parameter for inertia terms
     real(dp) :: alpha_f !< Generalized alpha parameters for stiffness terms
     real(dp) :: alpha   !< The HHT alpha-parameter for numeric damping
     real(dp) :: beta    !< Newmark time integration parameter
     real(dp) :: gamma   !< Newmark time integration parameter
     real(dp) :: svel    !< Scaling factor for velocity update
     real(dp) :: sacc    !< Scaling factor for acceleration update

     integer  :: maxIt   !< Maximum number of iterations per time step
     integer  :: minIt   !< Minimum number of iterations per time step
     integer  :: fixedIt !< Prescribed number of iterations per time step
     integer  :: nUpdat  !< Number of iterations with tangent matrix update
     !> Maximum number of sequental iterations without tangent matrix updates
     integer  :: maxSequentialNoUpdate
     !> Relative tolerance to determine if tangent matrix update is necessary
     real(dp) :: tolUpdateFactor

     integer(i8) :: nStep    !< Total number of time steps
     integer(i8) :: nIter    !< Total number of iterations updates
     integer(i8) :: nUpdates !< Total number of matrix updates

     integer  :: nUpdaThisStep !< Number of matrix updates for this step
     integer  :: nIterThisStep !< Number of iterations for this step
     real(dp) :: nIterPrevStep !< Number of iterations for previous step

     real(dp) :: equTol !< Initial static equilibrium convergence tolerance
     real(dp) :: equLim !< Initial static equilibrium step size limit

     logical  :: stressStiffIsOn(2) !< Geometric stiffness toggles

     real(dp) :: wDisp(3)  !< Norm scaling factors for Displacements
     real(dp) :: wForce(3) !< Norm scaling factors for Residual forces

     type(TestSetType) :: convergenceSet !< Iteration convergence data

     real(dp), pointer :: urd(:)    !< Velocity in global direction
     real(dp), pointer :: urdd(:)   !< Acceleration in global direction
     real(dp), pointer :: urdp(:)   !< Global velocity in previous step
     real(dp), pointer :: urddp(:)  !< Global acceleration in previous step

     real(dp), pointer :: dk(:)     !< Estimated velocity for predictor step
     real(dp), pointer :: ak(:)     !< Estimated acceleration for predictor step

     real(dp), pointer :: del(:)    !< RHS/Solution vector of the linear system
     real(dp), pointer :: sinc(:,:) !< Iterative solution increment
     real(dp), pointer :: rinc(:)   !< Total increment for current time step

     real(dp), pointer :: RFk(:)    !< Reaction force vector
     real(dp), pointer :: Qk(:)     !< External force vector in system direction
     real(dp), pointer :: Qprev(:)  !< External force vector of previous step
     real(dp), pointer :: FSk(:)    !< Internal stiffness force vector
     real(dp), pointer :: FDk(:)    !< Internal damping force vector
     real(dp), pointer :: FIk(:)    !< Internal inertia force vector
     real(dp), pointer :: FIprev(:)   !< Inertia forces of previous step
     real(dp), pointer :: FIactual(:) !< Inertia forces including the residual
     real(dp), pointer :: residual(:) !< Current residual force vector

     type(SysMatrixType) :: Nmat !< The system Newton matrix

  end type SystemType


  !> @brief Standard routine for writing an object to file.
  interface WriteObject
     module procedure WriteSystemType
  end interface


contains

  !!============================================================================
  !> @brief Standard routine for writing an object to io.
  !>
  !> @param[in] sys System level model data
  !> @param[in] io File unit number to write to
  !> @param[in] complexity Option controlling the amount of output
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 22 Jun 2000

  subroutine WriteSystemType (sys,io,complexity)

    use IdTypeModule       , only : getId
    use SysMatrixTypeModule, only : WriteObject

    type(SystemType), intent(in) :: sys
    integer         , intent(in) :: io
    integer,optional, intent(in) :: complexity

    !! --- Logic section ---

    write(io,'(A)') 'System','{'
    write(io,*) 'tStart   =', sys%tStart
    write(io,*) 'tEnd     =', sys%tEnd
    write(io,*) 'tQStatic =', sys%tQStatic
    write(io,*) 'tInc     =', sys%tInc
    write(io,*) 'maxInc   =', sys%maxInc
    write(io,*) 'minInc   =', sys%minInc
    write(io,*) 'cutback  =', sys%cutbck
    write(io,*) 'nCutStp  =', sys%nCutStp
    write(io,*) 'varInc   =', sys%varInc
    if (associated(sys%tIncEngine)) then
       write(io,*) 'tIncEngine =', trim(getId(sys%tIncEngine%id))
    end if
    if (associated(sys%fxItEngine)) then
       write(io,*) 'fxItEngine =', trim(getId(sys%fxItEngine%id))
    end if
    write(io,*) 'alpha_m  =', sys%alpha_m
    write(io,*) 'alpha_f  =', sys%alpha_f
    write(io,*) 'alpha    =', sys%alpha
    write(io,*) 'beta     =', sys%beta
    write(io,*) 'gamma    =', sys%gamma
    write(io,*) 'dynStressStiff?   =', sys%stressStiffIsOn(1)
    write(io,*) 'equStressStiff?   =', sys%stressStiffIsOn(2)
    write(io,*) 'maxIt    =', sys%maxIt
    write(io,*) 'minIt    =', sys%minIt
    write(io,*) 'fixedIt  =', sys%fixedIt
    write(io,*) 'nUpdat   =', sys%nUpdat
    write(io,*) 'maxSeqNU =', sys%maxSequentialNoUpdate
    write(io,*) 'tolUpdat =', sys%tolUpdateFactor
    write(io,*) 'equTol   =', sys%equTol
    write(io,*) 'equLim   =', sys%equLim
    write(io,'(A,1P3E13.5)') ' wDisp   =', sys%wDisp
    write(io,'(A,1P3E13.5)') ' wForce  =', sys%wForce
    if (present(complexity)) then
       if (complexity > 1) then
          !! Write sparse matrix data structure
          call writeObject (sys%Nmat,io,'Newton matrix',10+complexity)
       end if
    end if
    write(io,'(A)') '}'

  end subroutine WriteSystemType


  !!============================================================================
  !> @brief Initializes the SystemType object.
  !>
  !> @param[out] sys System level model data
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 25 Apr 2002

  subroutine NullifySys (sys)

    use NormTypeModule     , only : nullifyConvSet
    use SysMatrixTypeModule, only : nullifySysMatrix

    type(SystemType), intent(out) :: sys

    !! --- Logic section ---

    sys%time     = 0.0_dp
    sys%timeStep = 0.0_dp
    sys%timePrev = 0.0_dp

    sys%tStart = 0.0_dp
    sys%tEnd   = 0.0_dp
    sys%tQStatic = 0.0_dp
    sys%tInc   = 0.0_dp
    sys%maxInc = 0.0_dp
    sys%minInc = 0.0_dp
    sys%cutbck = 1.0_dp
    sys%nCutStp = 1
    sys%varInc = 0

    nullify(sys%tIncEngine)
    nullify(sys%fxItEngine)

    sys%alpha_m = 0.0_dp
    sys%alpha_f = 0.0_dp
    sys%alpha   = 0.0_dp
    sys%beta    = 0.0_dp
    sys%gamma   = 0.0_dp
    sys%svel    = 0.0_dp
    sys%sacc    = 0.0_dp

    sys%maxIt   = 0
    sys%minIt   = 0
    sys%fixedIt = 0
    sys%nUpdat  = 0

    sys%nStep    = 0_i8
    sys%nIter    = 0_i8
    sys%nUpdates = 0_i8
    sys%nUpdaThisStep = 0
    sys%nIterThisStep = 0
    sys%nIterPrevStep = 0.0_dp

    sys%equTol  = 0.0_dp
    sys%equLim  = 0.0_dp

    sys%stressStiffIsOn = .false.

    sys%wDisp  = 0.0_dp
    sys%wForce = 0.0_dp

    call nullifyConvSet (sys%convergenceSet)

    nullify(sys%urd)
    nullify(sys%urdd)
    nullify(sys%urdp)
    nullify(sys%urddp)
    nullify(sys%dk)
    nullify(sys%ak)
    nullify(sys%del)
    nullify(sys%sinc)
    nullify(sys%rinc)
    nullify(sys%RFk)
    nullify(sys%Qk)
    nullify(sys%Qprev)
    nullify(sys%FSk)
    nullify(sys%FDk)
    nullify(sys%FIk)
    nullify(sys%FIprev)
    nullify(sys%FIactual)
    nullify(sys%residual)

    call nullifySysMatrix (sys%Nmat)

  end subroutine NullifySys


  !!============================================================================
  !> @brief Deallocates the SystemType object.
  !>
  !> @param sys System level model data
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Jan 2017

  subroutine DeallocateSys (sys)

    use NormTypeModule     , only : deallocateConvSet
    use SysMatrixTypeModule, only : deallocateSysMatrix

    type(SystemType), intent(inout) :: sys

    !! --- Logic section ---

    call deallocateConvSet (sys%convergenceSet)
    call deallocateSysMatrix (sys%Nmat)

    if (associated(sys%urd))      deallocate(sys%urd)
    if (associated(sys%urdd))     deallocate(sys%urdd)
    if (associated(sys%urdp))     deallocate(sys%urdp)
    if (associated(sys%urddp))    deallocate(sys%urddp)
    if (associated(sys%dk))       deallocate(sys%dk)
    if (associated(sys%ak))       deallocate(sys%ak)
    if (associated(sys%del))      deallocate(sys%del)
    if (associated(sys%sinc))     deallocate(sys%sinc)
    if (associated(sys%rinc))     deallocate(sys%rinc)
    if (associated(sys%RFk))      deallocate(sys%RFk)
    if (associated(sys%Qk))       deallocate(sys%Qk)
    if (associated(sys%Qprev))    deallocate(sys%Qprev)
    if (associated(sys%FSk))      deallocate(sys%FSk)
    if (associated(sys%FDk))      deallocate(sys%FDk)
    if (associated(sys%FIk))      deallocate(sys%FIk)
    if (associated(sys%FIprev))   deallocate(sys%FIprev)
    if (associated(sys%FIactual)) deallocate(sys%FIactual)
    if (associated(sys%residual)) deallocate(sys%residual)

    call nullifySys (sys)

  end subroutine DeallocateSys


  !!============================================================================
  !> @brief Checks whether current step is quasi-static or dynamic.
  !>
  !> @param[in] sys System level model data
  !> @param[in] time Current simulation time
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 11 Oct 2022

  logical function isQuasiStatic (sys,time)

    type(SystemType)  , intent(in) :: sys
    real(dp), optional, intent(in) :: time

    !! --- Logic section ---

    if (present(time)) then
       isQuasiStatic = time < sys%tQStatic + 0.1_dp*sys%minInc
    else
       isQuasiStatic = sys%time < sys%tQStatic + 0.1_dp*sys%minInc
    end if

  end function isQuasiStatic

end module SystemTypeModule
