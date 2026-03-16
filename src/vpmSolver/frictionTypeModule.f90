!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file frictionTypeModule.f90
!>
!> @brief Friction object data containers.

!!==============================================================================
!> @brief Module with data types joint DOF friction objects.
!>
!> @details The module also contains subroutines for accessing the frictions.

module FrictionTypeModule

  use KindModule        , only : dp
  use IdTypeModule      , only : IdType
  use FunctionTypeModule, only : EngineType
  use SpringTypeModule  , only : SpringBaseType

  implicit none

  !> @brief Friction type names
  character(len=14),parameter :: fricType_p(6) = (/ 'ROT_FRICTION  ', &
       &                                            'TRANS_FRICTION', &
       &                                            'CAM_FRICTION  ', &
       &                                            'BALL_FRICTION ', &
       &                                            'BALL_FRICTION2', &
       &                                            'GENERIC_ENGINE' /)
  !> @cond NO_DOCUMENTATION
  integer, parameter ::                              ROT_FRICTION_p   = 1, &
       &                                             TRANS_FRICTION_p = 2, &
       &                                             CAM_FRICTION_p   = 3, &
       &                                             BALL_JNT_FRICTION_p = 4, &
       &                                             BALL_JNT_FRICTION2_p = 5, &
       &                                             GENERIC_ENGINE_p = 6
  !> @endcond

  real(dp), save :: fricForceTol = 0.0_dp !< Zero tolerance for friction forces


  !> @brief Data type representing a set of friction parameters.
  type FrictionParameterType

     type(IdType) :: id !< General identification data

     integer  :: type             !< Type of friction parameter set
     real(dp) :: typeDepParams(3) !< Joint type dependent parameters
     real(dp) :: CoulombCoeff     !< Coulomb friction coefficient
     real(dp) :: StribeckMagn     !< Magnitude of the Stribeck effect
     real(dp) :: StribeckSpeed    !< Critical speed of the Stribeck effect
     real(dp) :: ViscCoeff        !< Viscous coefficient
     real(dp) :: PrestressLoad    !< Prestress load
     real(dp) :: AsymMagn         !< Magnitude of asymmetries
     real(dp) :: FricAmpl         !< Friction amplitude
     real(dp) :: FricFreq         !< Frequency
     real(dp) :: FricPhase        !< Phase
     real(dp) :: stickStiffness   !< if &gt; 0.0, use spring element

  end type FrictionParameterType


  !> @brief Data type representing a joint DOF friction.
  type FrictionType

     type(FrictionParameterType), pointer :: param !< Friction properties
     type(EngineType)           , pointer :: eng   !< User-defined normal load
     !> Spring with variable yield to handle stick stiffness
     type(SpringBaseType)       , pointer :: spr

     !> Turn friction on/off, used by cam joints and contact elements
     logical :: isActive

     !> @brief Friction initialization flag.
     !> @details The value has the following interpretation:
     !>   - = 0: No initialization (within iteration loop)
     !>   - = 1: Initialize for a new time step
     !>   - = 2: Initialize for start of simulation
     integer :: lInit
     integer :: lSlip(2) !< Stick/Slip indicators

     real(dp) :: Fext  !< Force that should match friction force on stick
     real(dp) :: Fequ  !< Equivalent normal force
     real(dp) :: Fmax  !< Max friction force

     real(dp) :: force       !< Current force
     real(dp) :: forcePrevIt !< Force at previous iteration
     real(dp) :: dF          !< Change in force between last two iterations
     real(dp) :: stiff       !< Relative change in friction force wrt. velocity

     real(dp) :: pos       !< Current position
     real(dp) :: pos0      !< Position at start of time step
     real(dp) :: posPrevIt !< Position at previous iteration

     real(dp) :: vel       !< Current velocity
     real(dp) :: vel0      !< Velocity at start of timestep
     real(dp) :: velPrevIt !< Velocity at previous iteration

     real(dp) :: forcePrev !< Force at the end of last time step
     real(dp) :: posPrev   !< Position at the end of last time step
     real(dp) :: eDmp      !< Energy loss (always positive)

     logical  :: saveVar(3) !< Flags indicating which variables should be saved

  end type FrictionType


  !> @brief Data type representing a friction object pointer.
  !> @details This data type is used to construct arrays of frictions where each
  !> element is a pointer to a friction object, and not the objects themselves.
  type FrictionPtrType
     type(FrictionType), pointer :: p
  end type FrictionPtrType


  !> @brief Returns pointer to object with specified ID.
  interface GetPtrToId
     module procedure GetPtrToIdFrictionParameterType
  end interface

  !> @brief Standard routine for writing an object to file.
  interface WriteObject
     module procedure WriteFrictionType
     module procedure WriteFrictionParameterType
  end interface

  !> @brief Updates the state variables pertaining to previous time step.
  interface UpdateAtConvergence
     module procedure UpdateFrictionAtConvergence
  end interface

  !> @brief Restores the state variables from the last converged time step.
  interface RestoreFromLastStep
     module procedure RestoreFrictionFromLastStep
  end interface


contains

  !!============================================================================
  !> @brief Returns pointer to friction parameters object with specified id.
  !>
  !> @param[in] array Array of frictiontypemodule::frictionparametertype objects
  !>                  to search within
  !> @param[in] id Base ID of the object to search for
  !>
  !> @details If the friction parameters object is not found, NULL is returned.
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 1 Nov 1999

  function GetPtrToIdFrictionParameterType (array,id) result(ptr)

    integer                    , intent(in)         :: id
    type(FrictionParameterType), intent(in), target :: array(:)
    type(FrictionParameterType), pointer            :: ptr

    !! Local variables
    integer :: i

    !! --- Logic section ---

    do i = 1, size(array)
       if (array(i)%id%baseId == id) then
          ptr => array(i)
          return
       end if
    end do

    nullify(ptr)
    write(*,*) '*** GetPtrToIdFrictionParameterType returned nullified, id =',id

  end function GetPtrToIdFrictionParameterType


  !!============================================================================
  !> @brief Standard routine for writing an object to io.
  !>
  !> @param[in] friction The frictiontypemodule::frictiontype object to write
  !> @param[in] io File unit number to write to
  !> @param[in] complexity If present, the value indicates the amount of print
  !>
  !> @author Bjorn Haugen
  !>
  !> @date Sep 1999

  subroutine WriteFrictionType (friction,io,complexity)

    type(FrictionType), intent(in) :: friction
    integer           , intent(in) :: io
    integer, optional , intent(in) :: complexity

    !! --- Logic section ---

    write(io,'(A)') 'Friction','{'
    call writeObject (friction%param,io)

    write(io,*) 'saveVar       =', friction%saveVar
    if (present(complexity)) then
       if (complexity >= 2) then
          write(io,*) 'lSlip         =', friction%lSlip
          write(io,*) 'force         =', friction%force
          write(io,*) 'pos           =', friction%pos
          write(io,*) 'vel           =', friction%vel
          write(io,*) 'stiff         =', friction%stiff
          write(io,*) 'eDmp          =', friction%eDmp
       end if
    end if
    write(io,'(A)') '}'

  end subroutine WriteFrictionType


  !!============================================================================
  !> @brief Standard routine for writing an object to io.
  !>
  !> @param[in] fricData The frictiontypemodule::frictionparametertype
  !>                     object to write
  !> @param[in] io File unit number to write to
  !>
  !> @author Bjorn Haugen
  !>
  !> @date Sep 1999

  subroutine WriteFrictionParameterType (fricData,io)

    use IdTypeModule, only : writeId

    type(FrictionParameterType), intent(in) :: fricData
    integer                    , intent(in) :: io

    !! --- Logic section ---

    write(io,'(A)') 'FrictionData','{'
    call writeId (fricData%id,io)
    write(io,*) 'type          = ', fricType_p(fricData%type)
    write(io,*) 'typeDepParams = ', fricData%typeDepParams
    write(io,*) 'CoulombCoeff  = ', fricData%CoulombCoeff
    write(io,*) 'StribeckMagn  = ', fricData%StribeckMagn
    write(io,*) 'StribeckSpeed = ', fricData%StribeckSpeed
    write(io,*) 'ViscCoeff     = ', fricData%ViscCoeff
    write(io,*) 'PrestressLoad = ', fricData%PrestressLoad
    write(io,*) 'AsymMagn      = ', fricData%AsymMagn
    write(io,*) 'FricAmpl      = ', fricData%FricAmpl
    write(io,*) 'FricFreq      = ', fricData%FricFreq
    write(io,*) 'FricPhase     = ', fricData%FricPhase
    write(io,*) 'StickStifness = ', fricData%stickStiffness
    write(io,'(A)') '}'

  end subroutine WriteFrictionParameterType


  !!============================================================================
  !> @brief Initializes a friction object.
  !>
  !> @param[out] friction The frictiontypemodule::frictionttype
  !>                      object to initialize
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 16 Jun 2006

  subroutine NullifyFriction (friction)

    type(FrictionType), intent(out) :: friction

    !! --- Logic section ---

    nullify(friction%param)
    nullify(friction%eng)
    nullify(friction%spr)

    friction%isActive    = .true.
    friction%lInit       = 2
    friction%lSlip       = 0
    friction%Fext        = 0.0_dp
    friction%Fequ        = 0.0_dp
    friction%Fmax        = 0.0_dp
    friction%force       = 0.0_dp
    friction%forcePrevIt = 0.0_dp
    friction%dF          = 0.0_dp
    friction%pos         = 0.0_dp
    friction%pos0        = 0.0_dp
    friction%posPrevIt   = 0.0_dp
    friction%vel         = 0.0_dp
    friction%vel0        = 0.0_dp
    friction%velPrevIt   = 0.0_dp
    friction%stiff       = 0.0_dp
    friction%forcePrev   = 0.0_dp
    friction%posPrev     = 0.0_dp
    friction%eDmp        = 0.0_dp
    friction%saveVar     = .false.

  end subroutine NullifyFriction


  !!============================================================================
  !> @brief Deallocates an array of friction parameters objects.
  !>
  !> @param frictions The frictiontypemodule::frictionparametertype
  !>                  objects to deallocate
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Jan 2017

  subroutine DeallocateFrictionPrms (frictions)

    use IdTypeModule, only : deallocateId

    type(FrictionParameterType), pointer :: frictions(:)

    !! Local variables
    integer :: i

    !! --- Logic section ---

    do i = 1, size(frictions)
       call deallocateId (frictions(i)%id)
    end do
    deallocate(frictions)
    nullify(frictions)

  end subroutine DeallocateFrictionPrms


  !!============================================================================
  !> @brief Allocates and initializes a friction object.
  !>
  !> @param[out] friction The frictiontypemodule::frictionparameterttype
  !>                      object to allocate/initialize
  !> @param[in] frictionSets All friction parameters objects in the model
  !> @param[in] frictionId Base id of the friction parameters object to use
  !> @param[in] saveVar Flags indicating which variables should be saved
  !> @param[out] err Error flag
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 28 Aug 2002

  subroutine InitializeFriction (friction,frictionSets,frictionId,saveVar,err)

    use reportErrorModule     , only : AllocationError
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getbool

    type(FrictionType)         , pointer     :: friction
    type(FrictionParameterType), intent(in)  :: frictionSets(:)
    integer                    , intent(in)  :: frictionId, saveVar(:)
    integer                    , intent(out) :: err

    !! Local variables
    logical :: allSec, allRest, allFric, allForce, allEnerg

    !! --- Logic section ---

    allocate(friction,STAT=err)
    if (err /= 0) then
       err = AllocationError('InitializeFriction')
       return
    end if

    call NullifyFriction (friction)

    call ffa_cmdlinearg_getbool ('allSecondaryVars',allSec)
    call ffa_cmdlinearg_getbool ('allRestartVars',allRest)
    call ffa_cmdlinearg_getbool ('allFrictionVars',allFric)
    call ffa_cmdlinearg_getbool ('allForceVars',allForce)
    call ffa_cmdlinearg_getbool ('allEnergyVars',allEnerg)

    friction%saveVar(1) = saveVar(1) > 0 .or. allSec .or. allFric .or. allForce
    friction%saveVar(2) = saveVar(2) > 0 .or. allSec .or. allFric .or. allEnerg
    friction%saveVar(3) = allFric .or. allRest ! Additional debug/restart output

    friction%param => GetPtrToId(frictionSets,frictionId)
    if (.not. associated(friction%param)) err = 5

  end subroutine InitializeFriction


  !!============================================================================
  !> @brief Initializes internal friction variables at the start of simulation.
  !>
  !> @param friction The friction object to initialize internal variables for
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 24 Oct 2005

  subroutine UpdateFrictionAtStart (friction)

    type(FrictionType), intent(inout) :: friction

    !! --- Logic section ---

    friction%lInit = 0
    friction%lSlip = 0
    friction%pos0  = friction%pos
    friction%posPrevIt = friction%pos
    friction%velPrevIt = friction%vel

  end subroutine UpdateFrictionAtStart


  !!============================================================================
  !> @brief Updates a friction object after convergence has been achieved.
  !>
  !> @param friction The friction object to update
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 3 Oct 2005

  subroutine UpdateFrictionAtConvergence (friction)

    use SpringTypeModule, only : setYieldLimit

    type(FrictionType), intent(inout) :: friction

    !! --- Logic section ---

    friction%lInit = 1
    friction%pos0  = friction%pos

    if (associated(friction%spr)) then

       !! There is no need to recalculate friction%Fmax here because the
       !! velocities have not been updated since the last updateFrictions call
       !!TODO,kmo: Isn't it better to update this in every iteration also?
       call setYieldLimit (friction%spr,abs(friction%Fmax))

    end if

  end subroutine UpdateFrictionAtConvergence


  !!============================================================================
  !> @brief Restores a friction object from the last converged step.
  !>
  !> @param friction The friction object to restore
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 2 Nov 2008

  subroutine RestoreFrictionFromLastStep (friction)

    type(FrictionType), intent(inout) :: friction

    !! --- Logic section ---

    friction%lInit = 1
    friction%pos   = friction%pos0

  end subroutine RestoreFrictionFromLastStep

end module FrictionTypeModule
