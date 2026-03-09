!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file engineRoutinesModule.f90
!>
!> @brief Subroutines for general function evaluation.

!!==============================================================================
!> @brief Module with subroutines/functions for evaluation of general functions.
!>
!> @details This module contains functions and subroutines for evaluating
!> general functions, which are represented by functiontypemodule::enginetype
!> objects in the model.
!>
!> @note The enginevalue() function may be implicitly recursive
!> if a sensortypemodule::sensortype object defining an argument of
!> the general function is defined on a quantity which again may be controlled
!> by (another) general function. Currently, this applies to sensors on
!> spring- and damper quantities. Due to this, the update subroutines for
!> springs and dampers as well as sensors have all to be within this module,
!> in order to avoid circular use inclusions. No other (non-recursive)
!> subroutines that work on functiontypemodule::enginetype objects should be
!> placed in this module.

module EngineRoutinesModule

  use FunctionTypeModule   , only : EngineType
  use TriadTypeModule      , only : TriadType
  use SpringTypeModule     , only : SpringBaseType
  use DamperTypeModule     , only : DamperBaseType
  use EnvironmentTypeModule, only : EnvironmentType

  implicit none

  !> Public pointer to the environmental data member
  !> in the mechanismtypemodule::mechanismtype object
  type(EnvironmentType), save, public, pointer :: ourEnvir => null()

  !> Private pointer to the corresponding member array
  !> in the mechanismtypemodule::mechanismtype object
  type(EngineType)    , save, private, pointer :: engines(:) => null()
  !> @copydoc engineroutinesmodule::engines
  type(TriadType)     , save, private, pointer :: triads(:)  => null()
  !> @copydoc engineroutinesmodule::engines
  type(SpringBaseType), save, private, pointer :: springs(:) => null()
  !> @copydoc engineroutinesmodule::engines
  type(DamperBaseType), save, private, pointer :: dampers(:) => null()

  !> Flag used for consistent right-hand-side calculation in the predictor step.
  logical, save :: isPredictorStep = .false. !< Equals .true. in first iteration

  private :: EvalArgs, SensorRate


contains

  !!============================================================================
  !> @brief Initialization of private pointers.
  !>
  !> @param[in] eArr All functiontypemodule::enginetype objects in the model
  !> @param[in] tArr All triadtypemodule::triadtype objects in the model
  !> @param[in] sArr All springtypemodule::springbasetype objects in the model
  !> @param[in] dArr All dampertypemodule::damperbasetype objects in the model
  !>
  !> @details This subroutine sets the private pointers of this module to point
  !> to the arrays provided as arguments (which should be those within the
  !> mechanismtypemodule::mechanismtype object). This is to avoid the necessity
  !> to transport these arrays through the enginevalue() and evalargs() calls,
  !> since they are needed by the subroutine updatesensor() which is invoked
  !> from evalargs().
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 6 Jun 2002

  subroutine SetPointersForSensors (eArr,tArr,sArr,dArr)

    type(EngineType)    , intent(in), target :: eArr(:)
    type(TriadType)     , intent(in), target :: tArr(:)
    type(SpringBaseType), intent(in), target :: sArr(:)
    type(DamperBaseType), intent(in), target :: dArr(:)

    !! --- Logic section ---

    engines => eArr
    triads  => tArr
    springs => sArr
    dampers => dArr

  end subroutine SetPointersForSensors


  !!============================================================================
  !> @brief Pre-evaluation of general functions.
  !>
  !> @param engines All general functions in the model
  !> @param ierr Error flag
  !>
  !> @details This subroutine pre-evaluates all functiontypemodule::enginetype
  !> objects that have been flagged for it. It is used in the case that the
  !> order of evaluation matters, in particular for user-defined functions.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 6 Feb 2015

  subroutine preEvaluate (engines,ierr)

    use kindModule             , only : dp
    use explicitFunctionsModule, only : USER_DEFINED_p

    type(EngineType), intent(inout) :: engines(:)
    integer         , intent(inout) :: ierr

    !! Local variables
    integer  :: i
    real(dp) :: f

    !! --- Logic section ---

    do i = 1, size(engines)
       if (associated(engines(i)%func)) then
          if (engines(i)%func%type == USER_DEFINED_p) then
             if (mod(engines(i)%func%intParameters(4),2) > 0) then
                f = EngineValue(engines(i),ierr)
             end if
          end if
       end if
    end do

  end subroutine preEvaluate


  !!============================================================================
  !> @brief Evaluates a general function with optional scaling and offset.
  !>
  !> @param engine Points to the general function to evaluate
  !> @param[in] E1 Scaling factor
  !> @param[in] E0 Offset value
  !> @param ierr Error flag
  !> @param[in] xArg Optional function argument value
  !>
  !> @details
  !> - eVal = E0 + E1*func(x),  if @a func is defined for @a engine
  !> - eVal = E0 + E1*x,        if @a func is not defined for @a engine
  !>
  !> @a x may either be the input @a xArg or sensor-defined.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Karl Erik Thoresen / Bjorn Haugen / Knut Morten Okstad
  !>
  !> @date 6 Jun 2002

  function Evaluate (engine,E1,E0,ierr,xArg) result(eVal)

    use kindModule, only : dp

    type(EngineType) , pointer       :: engine
    real(dp)         , intent(in)    :: E1, E0
    integer          , intent(inout) :: ierr
    real(dp),optional, intent(in)    :: xArg
    real(dp)                         :: eVal

    !! --- Logic section ---

    if (associated(engine)) then
       eVal = E0 + E1 * EngineValue(engine,ierr,xArg)
    else if (present(xArg)) then
       eVal = E0 + E1 * xArg
    else
       eVal = E0
    end if

  end function Evaluate


  !!============================================================================
  !> @brief Evaluates a general function.
  !>
  !> @param engine The general function to evaluate
  !> @param ierr Error flag
  !> @param[in] xArg Optional function argument value
  !>
  !> @details
  !> - EngineValue = func(x),  if @a func is defined for @a engine
  !> - EngineValue = x,        if @a func is not defined for @a engine
  !>
  !> @a x may either be the input @a xArg or sensor-defined.
  !>
  !> @note
  !> This may be an implicit recursive function through the call sequence:
  !> @code
  !>    -> EvalArgs -> updateSensor -> updateSpringBase -> EngineValue
  !>                                   updateDamperBase
  !> @endcode
  !>
  !> @callgraph @callergraph
  !>
  !> @author Karl Erik Thoresen / Bjorn Haugen / Knut Morten Okstad
  !>
  !> @date 6 Jun 2002

  recursive function EngineValue (engine,ierr,xArg) result(eVal)

    use kindModule             , only : dp
    use FunctionTypeModule     , only : FunctionValue
    use IdTypeModule           , only : getId
    use explicitFunctionsModule, only : dbgFunc
    use reportErrorModule      , only : reportError, error_p

    type(EngineType) , intent(inout) :: engine
    integer          , intent(inout) :: ierr
    real(dp),optional, intent(in)    :: xArg

    !! Local variables
    integer  :: i, iDer, lerr, nArg
    real(dp) :: x(10), eVal, ramp

    !! Local stack to keep track of recursive calls
    integer, parameter :: maxStack_p = 10
    integer, save      :: topStack = 0, callStack(maxStack_p)

    !! --- Logic section ---

    eVal = 0.0_dp
    lerr = ierr

    !! Check for infinite recursive engine loop
    do i = 1, topStack
       if (engine%id%baseId == callStack(i)) then
          ierr = ierr - 1
          call reportError (error_p,'Infinite recursive loop for Engine'// &
               &            getId(engine%id),addString='EngineValue')
          return
       end if
    end do

    if (topStack >= maxStack_p) then
       ierr = ierr - 1
       call reportError (error_p,'Too many recursion levels detected for '// &
            &            'Engine'//getId(engine%id),addString='EngineValue')
       return
    end if

    !! Push current engine Id on the call stack
    topStack = topStack + 1
    callStack(topStack) = engine%id%baseId

    !! Define the argument value(s)
    iDer = 0
    nArg = size(engine%args)
    if (associated(engine%func)) then
       if (nArg == 0) then
          !! Check if the argument-less function also should be ramped
          !! assuming it is a time-history function
          if (engine%func%intParameters(2) == 3) nArg = -1
       else if (engine%func%intParameters(2) < 0) then
          iDer = -1 ! Deactivate ramping for this function
       end if
    end if
    call EvalArgs (engine%args,nArg,iDer,xArg,x,ramp,ierr)
    if (ierr < lerr) then
       topStack = topStack - 1 ! Pop the call stack
       call reportError (error_p,'Failed to evaluate arguments for engine'// &
            &                    getId(engine%id),addString='EngineValue')
       return
    end if

    !! Define the engine value
    if (nArg > 1) then
       eVal = ramp * FunctionValue(engine%func,x(1:nArg),ierr)
    else
       eVal = ramp * FunctionValue(engine%func,x(1),ierr)
    end if
    topStack = topStack - 1 ! Pop the call stack

    if (ierr < lerr) then
       call reportError (error_p,'Failed to evaluate Engine'//getId(engine%id),&
            &            addString='EngineValue')
    else if (dbgFunc > 0) then
       write(dbgFunc,"('Engine',A39,':',1PE12.5)") getId(engine%id),eVal
    end if

  end function EngineValue


  !!============================================================================
  !> @brief Evaluates the time-derivative of a general function.
  !>
  !> @param engine The general function to evaluate
  !> @param ierr Error flag
  !> @param[in] tArg Optional function argument value
  !>
  !> @callgraph @callergraph
  !>
  !> @author Trond Arne Svidal                                 @date 13 May 2002
  !> @author Knut Morten Okstad                                @date 17 Nov 2014

  recursive function EngineRate (engine,ierr,tArg) result(rVal)

    use kindModule        , only : dp
    use FunctionTypeModule, only : FunctionDerivative
    use IdTypeModule      , only : getId
    use reportErrorModule , only : reportError, error_p

    type(EngineType) , intent(inout) :: engine
    integer          , intent(inout) :: ierr
    real(dp),optional, intent(in)    :: tArg

    !! Local variables
    integer  :: i, iDer, lerr, nArg
    real(dp) :: x(10), rVal, ramp

    !! --- Logic section ---

    rVal = 0.0_dp
    lerr = ierr

    !! Define the argument value(s)
    iDer = 1
    nArg = size(engine%args)
    if (associated(engine%func)) then
       if (nArg == 0) then
          !! Check if the argument-less function also should be ramped
          !! assuming it is a time-history function
          if (engine%func%intParameters(2) == 3) nArg = -1
       else if (engine%func%intParameters(2) < 0) then
          iDer = -1 ! Deactivate ramping for this function
       end if
    end if
    call EvalArgs (engine%args,nArg,iDer,tArg,x,ramp,ierr)
    if (ierr < lerr) then
       call reportError (error_p,'Failed to evaluate arguments for engine'// &
            &                    getId(engine%id),addString='EngineRate')
       return
    else if (abs(ramp-1.0_dp) < 1.0e-15_dp) then
       ierr = ierr - 1
       call reportError (error_p,'Time-derivative for ramped engines '// &
            &                    'not yet implemented',addString='EngineRate')
       return
    end if

    !! Define the engine rate value
    if (present(tArg)) then
       rVal = FunctionDerivative(engine%func,x(1),1,ierr)
    else if (nArg == 1 .and. associated(engine%args(1)%p)) then
       rVal = FunctionDerivative(engine%func,x(1),1,ierr) &
            * SensorRate(engine%args(1)%p,ierr)
    else
       do i = 1, nArg
          if (associated(engine%args(i)%p)) then
             rVal = rVal + FunctionDerivative(engine%func,x(1:nArg),i,1,ierr) &
                  &      * SensorRate(engine%args(i)%p,ierr)
          end if
       end do
    end if

    if (ierr < lerr) then
       call reportError (error_p,'Failed to evaluate derivative for Engine'// &
            &                    getId(engine%id),addString='EngineRate')
    end if

  end function EngineRate


  !!============================================================================
  !> @brief Evaluates the time-derivative of a sensor value.
  !>
  !> @param sensor The sensor object to evaluate the derivative for
  !> @param ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 17 Nov 2014

  recursive function SensorRate (sensor,ierr) result (rVal)

    use kindModule       , only : dp
    use SensorTypeModule , only : SensorType
    use SensorTypeModule , only : TIME_p, NUM_ITERATIONS_p, ENGINE_p
    use IdTypeModule     , only : getId
    use reportErrorModule, only : internalError, reportError, error_p

    type(SensorType), intent(inout) :: sensor
    integer         , intent(inout) :: ierr

    !! Local variables
    real(dp) :: rVal

    !! --- Logic section ---

    select case (sensor%type)

    case (TIME_p)
       rVal = 1.0_dp

    case (NUM_ITERATIONS_p)
       rVal = 0.0_dp

    case (ENGINE_p)
       rVal = EngineRate(engines(sensor%index),ierr)

    case default
       rVal = 0.0_dp
       ierr = ierr - 1
       call reportError (error_p,'Can not evaluate time-derivative of Sensor'//&
            &            getId(sensor%id),addString='SensorRate')
    end select

  end function SensorRate


  !!============================================================================
  !> @brief Evaluates the argument(s) of a general function.
  !>
  !> @param args The function arguments to evaluate
  !> @param[in] nArg Length of array @a args
  !> @param[in] iDer Derivative order of ramp function (0, 1, or -1 for no ramp)
  !> @param[in] xArg Optional value of (first) function argument
  !> @param[out] x The evaluated argument value(s)
  !> @param[out] rVal Ramp function value or derivative (equals 1.0 if no ramp)
  !> @param ierr Error flag
  !>
  !> @details The output value @a x may either be the input argument @a xArg
  !> or a sensor-defined value. If a ramp function is defined and the argument
  !> type is @a TIME, the argument value @a x is instead set to the end time of
  !> the ramp function, and the ramp function value is returned through @a rVal.
  !>
  !> @note
  !> This may be an implicit recursive subroutine through the call sequence:
  !> @code
  !>    -> updateSensor -> updateSpringBase -> EngineValue -> EvalArgs
  !>                    -> updateDamperBase    EngineRate  /
  !> @endcode
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 22 Aug 2022

  recursive subroutine EvalArgs (args,nArg,iDer,xArg,x,rVal,ierr)

    use kindModule        , only : dp, hugeVal_p
    use FunctionTypeModule, only : FunctionValue, getRampValue
    use SensorTypeModule  , only : SensorPtrType, ourTime, TIME_p

    type(SensorPtrType), intent(inout) :: args(:)
    integer            , intent(in)    :: nArg, iDer
    real(dp), optional , intent(in)    :: xArg
    real(dp)           , intent(out)   :: x(:), rVal
    integer            , intent(inout) :: ierr

    !! Local variables
    integer :: i, iArg

    !! --- Logic section ---

    iArg = 1
    rVal = 1.0_dp
    if (present(xArg)) then
       iArg = 2
       x(1) = xArg
       if (iDer >= 0) then
          call getRampValue (x(1),iDer,rVal,ierr)
       end if
    else if (nArg < 1) then
       if (associated(ourTime) .and. nArg == -1 .and. iDer >= 0) then
          x(1) = ourTime%value
          call getRampValue (x(1),iDer,rVal,ierr)
       else
          x(1) = -hugeVal_p ! Should not be accessed but go bananas if attempted
       end if
    end if

    do i = iArg, nArg
       if (associated(args(i)%p)) then
          call UpdateSensor (args(i)%p,ierr)
          x(i) = args(i)%p%value
          if (args(i)%p%type == TIME_p .and. iDer >= 0) then
             call getRampValue (x(i),iDer,rVal,ierr)
          end if
       else
          x(i) = 0.0_dp
       end if
    end do

  end subroutine EvalArgs


  !!============================================================================
  !> @brief Updates a sensor value depending on the sensor type.
  !>
  !> @param sensor The sensor object to update the value for
  !> @param ierr Error flag
  !>
  !> @details For many sensor types the actual sensor value is a pointer
  !> directly into the object it is evaluating, and this subroutine does nothing
  !> for such sensors.
  !>
  !> @note For some sensors of non-trivial type,
  !> this may be an implicit recursive subroutine through the call sequences:
  !> @code
  !>    -> EngineValue --------------------\
  !>    -> updateSpringBase -> EngineValue -> EvalArgs -> updateSensor
  !>    -> updateDamperBase /
  !> @endcode
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 6 Jun 2002

  recursive subroutine UpdateSensor (sensor,ierr)

    use kindModule       , only : dp
    use SensorTypeModule , only : SensorType, ourTime
    use SensorTypeModule , only : TIME_p, ENGINE_p, CONTROL_p, MATLAB_WS_p
    use SensorTypeModule , only : JOINT_VARIABLE_p, RELATIVE_TRIAD_p, TRIAD_p
    use SensorTypeModule , only : SPRING_AXIAL_p, SPRING_JOINT_p
    use SensorTypeModule , only : DAMPER_AXIAL_p, DAMPER_JOINT_p
    use SensorTypeModule , only : STRAIN_GAGE_p, NUM_ITERATIONS_p
    use SensorTypeModule , only : LENGTH_p, VEL_p, ACC_p, LOCAL_p
    use SensorTypeModule , only : POS_p, REL_POS_p, FORCE_p
    use SensorTypeModule , only : W_SPEED_p, F_VEL_p, F_ACC_p, DYN_P_p
    use WindTurbineRoutinesModule, only : getWindSpeed
    use HydrodynamicsModule, only : getSeaState
    use TriadTypeModule  , only : transVSysToGlob
    use IdTypeModule     , only : getId
    use reportErrorModule, only : internalError, reportError, error_p

    type(SensorType), intent(inout) :: sensor
    integer         , intent(inout) :: ierr

    !! Local variables
    integer  :: lerr
    real(dp) :: work(7)

    !! --- Logic section ---

    lerr = ierr
    select case (sensor%type)

    case ( TIME_p, CONTROL_p, MATLAB_WS_p, JOINT_VARIABLE_p, &
         & STRAIN_GAGE_p, NUM_ITERATIONS_p )
       !! Do nothing, all handled through pointers to the appropriate values

    case (ENGINE_p)

       sensor%value = EngineValue(engines(sensor%index),ierr)

    case (SPRING_AXIAL_p, SPRING_JOINT_p)

       select case (sensor%entity)
       case (LENGTH_p)
          !! Do nothing, the spring length should already be up-to-date
       case default
          !! Update the variables of the associated spring
          call updateSpringBase (springs(sensor%index),ierr)
       end select

    case (DAMPER_AXIAL_p, DAMPER_JOINT_p)

       select case (sensor%entity)
       case (LENGTH_p, VEL_p)
          !! Do nothing, the damper length/velocity should already be up-to-date
       case default
          !! Update the variables of the associated damper
          call updateDamperBase (dampers(sensor%index),ierr)
       end select

    case (TRIAD_p)

       select case (sensor%system)
       case (LOCAL_p)
          !! Evaluate local velocity/acceleration/force component of the triad
          sensor%value = GetLocal(triads(sensor%index),sensor%dof,sensor%entity)
       case default
          !! Do nothing, global position, velocity and acceleration are all
          !! handled through pointers into the associated triad object.
          !! Except for the following quantities:
          select case (sensor%entity)
          case (POS_p)
             if (sensor%dof >= 4 .and. sensor%dof <= 9) then
                !! Evaluate global angular orientation variables
                sensor%value = GetAngle(triads(sensor%index),sensor%dof-3)
             end if
          case (FORCE_p)
             !! Evaluate global force component of the triad
             sensor%value = GetGlobalForce(triads(sensor%index),sensor%dof)
          case (W_SPEED_p)
             !! Evaluate the wind speed at the triad location
             call getWindSpeed (triads(sensor%index)%ur(:,4),ourTime%value, &
                  &             work(1:3),work(4:6),ierr)
             if (ierr < 0) ierr = lerr + ierr
             sensor%value = work(sensor%dof)
          case (F_VEL_p:DYN_P_p)
             !! Evaluate the fluid particle velocity/acceleration/pressure
             call getSeaState (ourEnvir,ourTime%value, &
                  &            triads(sensor%index)%ur(:,4),work,ierr)
             if (ierr < 0) ierr = lerr + ierr
             if (sensor%entity == F_VEL_p) then
                sensor%value = work(sensor%dof)
             else if (sensor%entity == F_ACC_p) then
                sensor%value = work(3+sensor%dof)
             else
                sensor%value = work(7)
             end if
          end select
       end select

    case (RELATIVE_TRIAD_p)

       !! Evaluate relative position/velocity/acceleration between two triads
       sensor%value = GetRelative(triads(sensor%index),triads(sensor%index2), &
            &                     sensor%dof,sensor%entity,sensor%value)

    case default
       ierr = ierr + internalError('UpdateSensor: Invalid sensor type')
    end select

    if (ierr < lerr) then
       call reportError (error_p,'Failed to update Sensor'//getId(sensor%id), &
            &            addString='UpdateSensor')
    end if

  contains

    !!==========================================================================
    !> @brief Evaluates a local velocity, acceleration or force for given triad.
    function GetLocal (triad,dof,entity)

      type(TriadType), intent(in) :: triad
      integer        , intent(in) :: dof, entity

      integer  :: iStart, iDir
      real(dp) :: GetLocal, vec(3)

      select case (dof)
      case (1:3)
         iStart = 1
         iDir   = dof
      case (4:6)
         iStart = 4
         iDir   = dof - 3
      case default
         ierr = ierr + internalError('UpdateSensor: Invalid sensor DOF')
         GetLocal = 0.0_dp
         return
      end select

      select case (entity)
      case (VEL_p)
         vec = triad%urd(iStart:iStart+2)
      case (ACC_p)
         vec = triad%urdd(iStart:iStart+2)
      case (FORCE_p)
         if (triad%isFEnode) then
            GetLocal = triad%supForce(iStart+iDir-1)
            return
         else if (triad%isSEnode) then
            vec = transVSysToGlob(triad,triad%supForce(iStart:iStart+2))
         else
            vec = transVSysToGlob(triad,triad%nodeForce(iStart:iStart+2))
         end if
      case default
         ierr = ierr + internalError('UpdateSensor: Invalid sensor entity')
         GetLocal = 0.0_dp
         return
      end select

      GetLocal = dot_product(triad%ur(:,iDir),vec)

    end function GetLocal

    !!==========================================================================
    !> @brief Evaluates a global force for given triad.
    function GetGlobalForce (triad,dof)

      type(TriadType), intent(in) :: triad
      integer        , intent(in) :: dof

      integer  :: iStart, iDir
      real(dp) :: GetGlobalForce, vec(3)

      select case (dof)
      case (1:3)
         iStart = 1
         iDir   = dof
      case (4:6)
         iStart = 4
         iDir   = dof - 3
      case default
         ierr = ierr + internalError('UpdateSensor: Invalid sensor DOF')
         GetGlobalForce = 0.0_dp
         return
      end select

      if (triad%isFEnode) then
         GetGlobalForce = triad%supForce(iStart+iDir-1)
      else if (triad%isSEnode) then
         vec = transVSysToGlob(triad,triad%supForce(iStart:iStart+2))
         GetGlobalForce = vec(iDir)
      else
         vec = transVSysToGlob(triad,triad%nodeForce(iStart:iStart+2))
         GetGlobalForce = vec(iDir)
      end if

    end function GetGlobalForce

    !!==========================================================================
    !> @brief Evaluates a global rotation angle for given triad.
    function GetAngle (triad,dof)

      use rotationModule, only : FFa_glbEulerZYX, mat_to_vec

      type(TriadType), intent(in) :: triad
      integer        , intent(in) :: dof

      real(dp) :: GetAngle, vec(3)

      if (dof > 0 .and. dof <= 3) then
         !! X, Y or Z-rotation (Euler ZYX)
         call FFa_glbEulerZYX (triad%ur(:,1:3),vec)
         GetAngle = vec(dof)
      else if (dof > 3 .and. dof <= 6) then
         !! X, Y or Z-rotation (Rodrigues parametrization)
         call mat_to_vec (triad%ur(:,1:3),vec)
         GetAngle = vec(dof-3)
      end if

    end function GetAngle

    !!==========================================================================
    !> @brief Evaluates a relative quantity between given two triads.
    function GetRelative (triad1,triad2,dof,entity,oldValue)

      use rotationModule, only : FFa_eulerZYX, deltaRot

      type(TriadType), intent(in) :: triad1, triad2
      integer        , intent(in) :: dof, entity
      real(dp)       , intent(in) :: oldValue

      real(dp)            :: GetRelative, relDis, relPos(3), vec(3)
      real(dp), parameter :: eps_p = 1.0e-12_dp

      select case (entity)

      case (REL_POS_p) ! Relative positions

         select case (dof)

         case (1:3) ! X, Y or Z-distance
            GetRelative = triad2%ur(dof,4) - triad1%ur(dof,4)

         case (4:6) ! X, Y or Z-rotation (Euler ZYX)
            call FFa_eulerZYX (triad1%ur(:,1:3),triad2%ur(:,1:3),relpos)
            GetRelative = relPos(dof-3)

         case (7:9) ! X, Y or Z-rotation (Rodrigues parametrization)
            relPos = deltaRot(triad1%ur(:,1:3),triad2%ur(:,1:3))
            GetRelative = relPos(dof-6)

         case (0)   ! relative distance
            relPos = triad2%ur(:,4) - triad1%ur(:,4)
            GetRelative = sqrt(dot_product(relPos,relPos))

         case default
            ierr = ierr + internalError('UpdateSensor: Invalid sensor DOF')
            GetRelative = 0.0_dp
         end select

      case (VEL_p) ! Relative velocities

         select case (dof)

         case (1:6) ! X, Y or Z-velocity
            GetRelative = triad2%urd(dof) - triad1%urd(dof)

         case (0)   ! relative velocity
            relPos = triad2%ur(:,4) - triad1%ur(:,4)
            relDis = dot_product(relPos,relPos)
            if (relDis < eps_p) then
               GetRelative = oldValue ! Keep the old value
            else
               vec = triad2%urd(1:3) - triad1%urd(1:3)
               GetRelative = dot_product(vec,relPos)/sqrt(relDis)
            end if

         case default
            ierr = ierr + internalError('UpdateSensor: Invalid sensor DOF')
            GetRelative = 0.0_dp
         end select

      case (ACC_p) ! Relative accelerations

         select case (dof)

         case (1:6) ! X, Y or Z-acceleration
            GetRelative = triad2%urdd(dof) - triad1%urdd(dof)

         case (0)   ! relative acceleration
            relPos = triad2%ur(:,4) - triad1%ur(:,4)
            relDis = dot_product(relPos,relPos)
            if (relDis < eps_p) then
               GetRelative = oldValue ! Keep the old value
            else
               vec = triad2%urdd(1:3) - triad1%urdd(1:3)
               GetRelative = dot_product(vec,relPos)/sqrt(relDis)
            end if

         case default
            ierr = ierr + internalError('UpdateSensor: Invalid sensor DOF')
            GetRelative = 0.0_dp
         end select

      case default
         ierr = ierr + internalError('UpdateSensor: Invalid sensor entity')
         GetRelative = 0.0_dp
      end select

    end function GetRelative

  end subroutine UpdateSensor


  !!============================================================================
  !> @brief Updates the spring variables.
  !>
  !> @param spr The spring object to update the variables for
  !> @param ierr Error flag
  !>
  !> @note
  !> This may be an implicit recursive subroutine through the call sequence:
  !> @code
  !>    -> EngineValue -> EvalArgs -> updateSensor -> updateSpringBase
  !> @endcode
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 6 Jun 2002

  recursive subroutine UpdateSpringBase (spr,ierr)

    use kindModule       , only : dp
    use SpringTypeModule , only : checkYieldLimits, springForce, springStiff
    use IdTypeModule     , only : getId
    use reportErrorModule, only : reportError, error_p

    type(SpringBaseType), intent(inout) :: spr
    integer             , intent(inout) :: ierr

    !! Local variables
    integer  :: lerr
    real(dp) :: scale

    !! --- Logic section ---

    if (.not. spr%isActive) then
       if (associated(spr%length)) spr%length = 0.0_dp
       spr%deflection = 0.0_dp
       spr%force      = 0.0_dp
       spr%stiffness  = 0.0_dp
       spr%Estr       = 0.0_dp
       return
    end if

    lerr = ierr

    !! Get the engine value for the stress free length
    if (associated(spr%length0Engine)) then
       spr%length0 = spr%l0 + spr%l1 * EngineValue(spr%length0Engine,ierr)
    else
       spr%length0 = spr%l0
    end if
    if (ierr < lerr) then
       call reportError (error_p,'Failed to evaluate stress-free length '// &
            &                    'engine for Spring'//getId(spr%id), &
            &            addString='UpdateSpringBase')
       return
    end if

    !! Evaluate the spring deflection, force and stiffness
    spr%deflection = spr%length - spr%length0
    spr%force      = springForce(spr,spr%deflection-spr%yieldDeflection,ierr)
    spr%stiffness  = springStiff(spr,spr%deflection-spr%yieldDeflection,ierr)
    if (ierr < lerr) return

    !! Possible scaling of the stiffness and force
    if (spr%deflection >= 0.0_dp) then
       if (associated(spr%stiffScaleEnginePos)) then
          scale = spr%scale1 * EngineValue(spr%stiffScaleEnginePos,ierr)
       else
          scale = spr%scale1
       end if
    else
       if (associated(spr%stiffScaleEngineNeg)) then
          scale = spr%scale1 * EngineValue(spr%stiffScaleEngineNeg,ierr)
       else
          scale = spr%scale1
       end if
    end if

    if (ierr < lerr) then
       call reportError (error_p,'Failed to evaluate stiffness scaling '// &
            &                    'engine for Spring'//getId(spr%id), &
            &            addString='UpdateSpringBase')
       return
    else
       spr%stiffness = spr%stiffness*scale
       spr%force     = spr%force*scale
    end if

    if (associated(spr%yield)) then
       !! Adjust force and stiffness if the spring force is above yield limit
       call checkYieldLimits (spr,spr%yield%pos,spr%yield%neg)
    end if

  end subroutine UpdateSpringBase


  !!============================================================================
  !> @brief Updates the damper variables.
  !>
  !> @param dmp The damper object to update the variables for
  !> @param ierr Error flag
  !>
  !> @note
  !> This may be an implicit recursive subroutine through the call sequence:
  !> @code
  !>    -> EngineValue -> EvalArgs -> updateSensor -> updateDamperBase
  !> @endcode
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 6 Jun 2002

  recursive subroutine UpdateDamperBase (dmp,ierr)

    use kindModule       , only : dp
    use DamperTypeModule , only : DamperBaseType, damperForce, damperCoeff
    use IdTypeModule     , only : getId
    use reportErrorModule, only : reportError, error_p

    type(DamperBaseType), intent(inout) :: dmp
    integer             , intent(inout) :: ierr

    !! Local variables
    integer  :: lerr
    real(dp) :: scale

    !! --- Logic section ---

    if (.not. dmp%isActive) then
       if (associated(dmp%velocity)) dmp%velocity = 0.0_dp
       if (associated(dmp%length))   dmp%length   = 0.0_dp
       dmp%coeff = 0.0_dp
       dmp%force = 0.0_dp
       return
    end if

    lerr = ierr

    !! Evaluate the damping coefficient and force
    dmp%coeff = damperCoeff(dmp,ierr)
    if (.not. isPredictorStep) then
       dmp%force = damperForce(dmp,ierr)
    else if (associated(dmp%velocity)) then
       dmp%force = dmp%coeff*dmp%velocity
    end if
    if (ierr < lerr) return

    !! Possible scaling of the damping coefficient and force
    if (associated(dmp%coeffScaleEngine)) then
       scale = dmp%scale1 * EngineValue(dmp%coeffScaleEngine,ierr)
    else
       scale = dmp%scale1
    end if

    if (ierr < lerr) then
       call reportError (error_p,'Failed to evaluate damper-coefficient '// &
            &                    'scaling engine for Damper'//getId(dmp%id), &
            &            addString='UpdateDamperBase')
    else
       dmp%coeff = dmp%coeff*scale
       dmp%force = dmp%force*scale
    end if

  end subroutine UpdateDamperBase


  !!============================================================================
  !> @brief Update all engines and store the values for saving.
  !>
  !> @param engines All general functions in the model
  !> @param[in] eFlag Update the functions whose @a saveVar equals this value
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 20 Jul 2004

  subroutine UpdateEnginesForSave (engines,eFlag,ierr)

    use FunctionTypeModule, only : EngineType
    use reportErrorModule , only : reportError, debugFileOnly_p

    type(EngineType), intent(inout) :: engines(:)
    integer         , intent(in)    :: eFlag
    integer         , intent(out)   :: ierr

    !! Local variables
    integer :: i

    !! --- Logic section ---

    ierr = 0
    do i = 1, size(engines)
       if (engines(i)%saveVar == eFlag) then
          engines(i)%value = EngineValue(engines(i),ierr)
       end if
    end do
    if (ierr < 0) call reportError (debugFileOnly_p,'UpdateEnginesForSave')

  end subroutine UpdateEnginesForSave


  !!============================================================================
  !> @brief Updates and prints out the current value of all general functions.
  !>
  !> @param[in] istep Time increment counter
  !> @param[in] time Current time
  !> @param engines All general functions in the model
  !> @param[in] lpu File unit number for res-file output
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 12 Jun 2023

  subroutine PrintFuncValues (istep,time,engines,lpu)

    use FunctionTypeModule, only : EngineType, dp, getEngineId
    use reportErrorModule , only : reportError, warning_p

    integer         , intent(in)    :: istep
    real(dp)        , intent(in)    :: time
    type(EngineType), intent(inout) :: engines(:)
    integer         , intent(in)    :: lpu

    !! Local variables
    integer :: i, ierr

    !! --- Logic section ---

    ierr = 0
    do i = 1, size(engines)
       engines(i)%value = EngineValue(engines(i),ierr)
    end do
    if (ierr < 0) then
       call reportError (warning_p,'Function evaluation failed',ierr=ierr)
       if (-ierr <= size(engines)) return
    end if

    write(lpu,600) istep+1, time
    do i = 1, size(engines)
       write(lpu,610) getEngineId(engines(i)), engines(i)%value
    end do

600 format(/'  >>> Function values before step',I6,'  time =',1PE13.5,' <<<')
610 format(5X,A32,' :',1PE13.5)

  end subroutine PrintFuncValues

end module EngineRoutinesModule
