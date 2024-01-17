!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file functionTypeModule.f90
!> @brief Function data container.

!!==============================================================================
!> @brief Module with data types representing function objects.
!>
!> @details The module also contains subroutines for accessing or manipulating
!> the function data, and for evaluating the functions themselves.

module FunctionTypeModule

  use IdTypeModule    , only : IdType, ldesc_p
  use SensorTypeModule, only : SensorPtrType, dp

  implicit none

  !> @brief Data type representing a function shape object.
  type FunctionType
     type(IdType)      :: id                 !< General identification data
     integer           :: type               !< Function type identifier
     integer           :: nArg               !< Number of function arguments
     integer , pointer :: intParameters(:)   !< Type-dependent integer params.
     real(dp), pointer :: realParameters(:)  !< Type-dependent real parameters
  end type FunctionType

  !> @brief Data type representing a general function object.
  type EngineType
     type(IdType)                 :: id      !< General identification data
     type(FunctionType) , pointer :: func    !< Pointer to the function shape
     type(SensorPtrType), pointer :: args(:) !< Function argument sensors
     character(len=ldesc_p)       :: tag     !< Object identification tag
     real(dp)                     :: value   !< Current function value for save
     integer                      :: saveVar !< Flag indicating save to frs
     logical                      :: isUsed  !< Indicates that function is used
  end type EngineType

  !> @brief Data type representing a general function pointer.
  !> @details This data type is used to construct arrays of function objects
  !> where each array element is a pointer to a function object,
  !> and not the functions themselves.
  type EnginePtrType
     type(EngineType), pointer :: p !< Pointer to a general function
  end type EnginePtrType

  !> Ramp function used for scaling of all general functions of time
  type(FunctionType), pointer, save :: ourRamp => null()


  !> @brief Standard routine for writing an object to file.
  interface WriteObject
     module procedure WriteFunctionType
     module procedure WriteEngineType
  end interface

  !> @brief Returns pointer to object with specified ID.
  interface GetPtrToId
     module procedure GetPtrToIdFunction
     module procedure GetPtrToIdEngine
  end interface

  !> @brief Returns the function value.
  interface FunctionValue
     module procedure FunctionValue1
     module procedure FunctionValue2
  end interface

  !> @brief Returns a function derivative.
  interface FunctionDerivative
     module procedure FunctionDerivative1
     module procedure FunctionDerivative2
  end interface

  !> @brief Deallocates a function object.
  interface DeallocateFunctions
     module procedure DeallocateFuncs
     module procedure DeallocateEngines
  end interface


contains

  !!============================================================================
  !> @brief Returns pointer to (first) function shape object with specified ID.
  !>
  !> @param[in] array Array of functiontypemodule::functiontype objects
  !>                  to search within
  !> @param[in] id Base ID of the object to search for
  !>
  !> @details If the function shape is not found, NULL is returned.
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 1 Nov 1999

  function GetPtrToIdFunction (array,id) result(ptr)
    type(FunctionType), pointer :: ptr

    integer           , intent(in)        :: id
    type(FunctionType), intent(in),target :: array(:)

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
    write(*,*) '*** GetPtrToIdFunction returned nullified, baseId =',id

  end function GetPtrToIdFunction


  !!============================================================================
  !> @brief Returns pointer to (first) general function with specified ID.
  !>
  !> @param[in] array Array of functiontypemodule::enginetype objects
  !>                  to search within
  !> @param[in] id Base ID of the object to search for
  !> @param[in] usedByMech If .true., flag the found object as used
  !> @param[out] index The array index of the found object
  !>
  !> @details If the general function is not found, NULL is returned.
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 1 Nov 1999

  function GetPtrToIdEngine (array,id,usedByMech,index) result(ptr)
    type(EngineType), pointer :: ptr

    integer         , intent(in)            :: id
    type(EngineType), intent(in) , target   :: array(:)
    logical         , intent(in) , optional :: usedByMech
    integer         , intent(out), optional :: index

    !! Local variables
    integer :: i

    !! --- Logic section ---

    do i = 1, size(array)
       if (array(i)%id%baseId == id) then
          ptr => array(i)
          if (present(usedByMech)) then
             if (usedByMech) ptr%isUsed = .true.
          else
             ptr%isUsed = .true.
          end if
          if (present(index)) index = i
          return
       end if
    end do

    nullify(ptr)
    write(*,*) '*** GetPtrToIdEngine returned nullified, baseId =',id

  end function GetPtrToIdEngine


  !!============================================================================
  !> @brief Initializes a function shape object.
  !>
  !> @param func The functiontypemodule::functiontype object to initialize
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 10 Sep 1998

  subroutine NullifyFunction (func)

    use IdTypeModule, only : nullifyId

    type(FunctionType), intent(out) :: func

    !! --- Logic section ---

    call nullifyId (func%id)

    func%type = 0
    func%nArg = 1
    nullify(func%intParameters)
    nullify(func%realParameters)

  end subroutine NullifyFunction


  !!============================================================================
  !> @brief Initializes a general function object.
  !>
  !> @param engine The functiontypemodule::enginetype object to initialize
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 10 Sep 1998

  subroutine NullifyEngine (engine)

    use IdTypeModule, only : nullifyId

    type(EngineType), intent(out) :: engine

    !! --- Logic section ---

    call nullifyId (engine%id)

    nullify(engine%func)
    nullify(engine%args)
    engine%tag = ''
    engine%value = 0.0_dp
    engine%saveVar = 0
    engine%isUsed = .false.

  end subroutine NullifyEngine


  !!============================================================================
  !> @brief Deallocates a function shape object.
  !>
  !> @param func The functiontypemodule::functiontype object to deallocate
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Jan 2017

  subroutine DeallocateFunction (func)

    use IdTypeModule, only : deallocateId

    type(FunctionType), intent(inout) :: func

    !! --- Logic section ---

    call deallocateId (func%id)

    if (associated(func%intParameters))  deallocate(func%intParameters)
    if (associated(func%realParameters)) deallocate(func%realParameters)

    call nullifyFunction (func)

  end subroutine DeallocateFunction


  !!============================================================================
  !> @brief Deallocates all function shape objects.
  !>
  !> @param funcs Array of all functiontypemodule::functiontype objects
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Jan 2017

  subroutine deallocateFuncs (funcs)

    type(FunctionType), pointer :: funcs(:)

    !! Local variables
    integer :: i

    !! --- Logic section ---

    if (associated(ourRamp)) then
       call DeallocateFunction (ourRamp)
       deallocate(ourRamp)
       nullify(ourRamp)
    end if

    do i = 1, size(funcs)
       call DeallocateFunction (funcs(i))
    end do

    deallocate(funcs)
    nullify(funcs)

  end subroutine deallocateFuncs


  !!============================================================================
  !> @brief Deallocates a general function object.
  !>
  !> @param engine The functiontypemodule::enginetype object to deallocate
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Jan 2017

  subroutine DeallocateEngine (engine)

    use IdTypeModule, only : deallocateId

    type(EngineType), intent(inout) :: engine

    !! --- Logic section ---

    call deallocateId (engine%id)

    if (associated(engine%args)) deallocate(engine%args)

    call nullifyEngine (engine)

  end subroutine DeallocateEngine


  !!============================================================================
  !> @brief Deallocates all general function objects.
  !>
  !> @param engines Array of all functiontypemodule::enginetype objects
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Jan 2017

  subroutine deallocateEngines (engines)

    type(EngineType), pointer :: engines(:)

    !! Local variables
    integer :: i

    !! --- Logic section ---

    do i = 1, size(engines)
       call DeallocateEngine (engines(i))
    end do

    deallocate(engines)
    nullify(engines)

  end subroutine deallocateEngines


  !!============================================================================
  !> @brief Standard routine for writing an object to file.
  !>
  !> @param[in] func The functiontypemodule::functiontype object to write
  !> @param[in] lpu File unit number to write to
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 8 Sep 1999

  subroutine WriteFunctionType (func,lpu)

    use IdTypeModule           , only : writeId
    use explicitFunctionsModule, only : funcType_p

    type(FunctionType), intent(in) :: func
    integer           , intent(in) :: lpu

    !! --- Logic section ---

    write(lpu,'(A)') 'Function','{'
    call writeId (func%id,lpu)

    write(lpu,*) 'type        = ', trim(funcType_p(func%type))
    write(lpu,*) 'nArg        = ', func%nArg
    if (associated(func%intParameters)) then
       write(lpu,*) 'intParameters  =', func%intParameters
    end if
    if (associated(func%realParameters)) then
       write(lpu,*) 'realParameters =', func%realParameters
    end if

    write(lpu,'(A)') '}'

  end subroutine WriteFunctionType


  !!============================================================================
  !> @brief Standard routine for writing an object to file.
  !>
  !> @param[in] engine The functiontypemodule::enginetype object to write
  !> @param[in] lpu File unit number to write to
  !>
  !> @author Karl Erik Thoresen
  !>
  !> @date 27 Sep 1998

  subroutine WriteEngineType (engine,lpu)

    use IdTypeModule    , only : writeId
    use SensorTypeModule, only : getSensorId

    type(EngineType), intent(in) :: engine
    integer         , intent(in) :: lpu

    !! Local variables
    integer :: i

    !! --- Logic section ---

    write(lpu,'(A)') 'Engine','{'
    call writeId (engine%id,lpu)
    if (engine%tag /= '') write(lpu,*) 'tag         =', engine%tag

    if (associated(engine%args)) then
       do i = 1, size(engine%args)
          if (associated(engine%args(i)%p)) then
             write(lpu,*) 'argSensor   = ', &
                  &       trim(getSensorId(engine%args(i)%p,.true.))
          else
             write(lpu,*) 'argSensor   = NULL'
          end if
       end do
    end if

    if (associated(engine%func)) then
       write(lpu,*) 'func        = ', trim(getFunctionId(engine%func))
    else
       write(lpu,*) 'func        = NULL'
    end if

    write(lpu,'(A)') '}'

  end subroutine WriteEngineType


  !!============================================================================
  !> @brief Returns the id of a function object, including its type.
  !>
  !> @param[in] func The functiontypemodule::functiontype object to get id for
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 8 Feb 2024

  function getFunctionId (func) result(text)

    use IdTypeModule           , only : lId_p, getId
    use explicitFunctionsModule, only : funcType_p, maxFunc_p

    type(FunctionType), intent(in) :: func

    !! Local variables
    character(len=lId_p) :: text

    !! --- Logic section ---

    text = getId(func%id)
    if (func%type > 0 .and. func%type <= maxFunc_p) then
       text = trim(funcType_p(func%type))//text
    else
       text = 'Function'//text(1:lId_p-8)
    end if

  end function getFunctionId


  !!============================================================================
  !> @brief Returns the id of an engine object, including function type.
  !>
  !> @param[in] engine The functiontypemodule::enginetype object to get id for
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 12 Jun 2023

  function getEngineId (engine) result(text)

    use IdTypeModule           , only : lId_p, getId
    use explicitFunctionsModule, only : funcType_p, maxFunc_p

    type(EngineType), intent(in) :: engine

    !! Local variables
    integer              :: itype
    character(len=lId_p) :: text

    !! --- Logic section ---

    text = getId(engine%id)
    if (associated(engine%func)) then
       itype = engine%func%type
       if (itype > 0 .and. itype <= maxFunc_p) then
          text = trim(text)//' '//trim(funcType_p(itype))
       end if
    end if

  end function getEngineId


  !!============================================================================
  !> @brief Evaluates the single-argument function f = f(x).
  !>
  !> @param func Pointer to the function shape to evaluate
  !> @param[in] xArg Function argument value to evaluate the function at
  !> @param ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 6 Nov 1998

  function FunctionValue1 (func,xArg,ierr) result(fVal)

    use explicitFunctionsModule, only : funcValue
    use reportErrorModule      , only : getErrorFile, reportError, error_p
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_getint

    type(FunctionType), pointer       :: func
    real(dp)          , intent(in)    :: xArg
    integer           , intent(inout) :: ierr

    !! Local variables
    real(dp)          :: fVal
    integer           :: lerr, iprint, lpu
    character(len=64) :: errMsg

    !! --- Logic section ---

    fVal = xArg
    if (.not.associated(func)) return
    if (func%id%baseId <= 0) return

    lerr = ierr
    call funcValue (func%id%baseId, func%intParameters, &
         &          func%realParameters, xArg, fVal, ierr)
    call ffa_cmdlinearg_getint ('debug',iprint)
    if (iprint == -3) iprint = 3
    if (iprint > 2 .or. ierr < lerr .or. isNaN(fVal)) then
       if (func%narg > 0) then
          write(errMsg,"(' at x =',1PE12.5)") xArg
       else
          errMsg = ''
       end if
    end if
    if (ierr < lerr) then
       call reportError (error_p,'Unable to evaluate '// &
            &            trim(getFunctionId(func))//errMsg)
       return
    end if
    if (iprint > 2) then
       lpu = getErrorFile()
       write(lpu,600) trim(getFunctionId(func)),trim(errMsg),fVal
600    format(5X,2A,': Value =',1PE12.5)
    end if
    if (isNaN(fVal)) then
       ierr = ierr - 1
       call reportError (error_p,trim(getFunctionId(func))//errMsg// &
            &            ': Returned NaN','Check your input')
    end if

  end function FunctionValue1


  !!============================================================================
  !> @brief Evaluates the multi-argument function f = f(x,y,z,...).
  !>
  !> @param func Pointer to the function shape to evaluate
  !> @param[in] xArg Function argument values to evaluate the function at
  !> @param ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 6 Nov 1998

  function FunctionValue2 (func,xArg,ierr) result(fVal)

    use explicitFunctionsModule, only : funcValue
    use reportErrorModule      , only : getErrorFile, reportError, error_p
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_getint

    type(FunctionType), pointer       :: func
    real(dp)          , intent(in)    :: xArg(:)
    integer           , intent(inout) :: ierr

    !! Local variables
    real(dp)          :: fVal, x
    integer           :: lerr, iprint, lpu
    character(len=64) :: errMsg

    !! --- Logic section ---

    if (size(xArg) > 0) then
       fVal = xArg(1)
    else
       fVal = 0.0_dp
    end if
    if (.not.associated(func)) return
    if (func%id%baseId <= 0) return

    lerr = ierr
    if (size(xArg) > 1) then
       call funcValue (func%id%baseId, func%intParameters, &
            &          func%realParameters, xArg, fVal, ierr)
    else
       x = fVal
       call funcValue (func%id%baseId, func%intParameters, &
            &          func%realParameters, x, fVal, ierr)
    end if
    call ffa_cmdlinearg_getint ('debug',iprint)
    if (iprint == -3) iprint = 3
    if (iprint > 2 .or. ierr < lerr .or. isNaN(fVal)) then
       if (func%narg > 0) then
          write(errMsg,"(' at x =',1P,9E12.5)") xArg
       else
          errMsg = ''
       end if
    end if
    if (ierr < lerr) then
       call reportError (error_p,'Unable to evaluate '// &
            &            trim(getFunctionId(func))//errMsg)
       return
    end if
    if (iprint > 2) then
       lpu = getErrorFile()
       write(lpu,600) trim(getFunctionId(func)),trim(errMsg),fVal
600    format(5X,2A,': Value =',1PE12.5)
    end if
    if (isNaN(fVal)) then
       ierr = ierr - 1
       call reportError (error_p,trim(getFunctionId(func))//errMsg// &
            &            ': Returned NaN','Check your input')
    end if

  end function FunctionValue2


  !!============================================================================
  !> @brief Evaluates the derivative of a single-argument function f(x).
  !>
  !> @param func Pointer to the function shape to evaluate the derivative of
  !> @param[in] xArg Function argument value to evaluate the derivative at
  !> @param[in] derOrder Derivative order (0, 1 or 2)
  !> @param ierr Error flag
  !> @param[in] eps Size of argument interval for numerical derivatives
  !>
  !> @callgraph @callergraph
  !>
  !> @author Trond Arne Svidal
  !>
  !> @date Oct 2002

  function FunctionDerivative1 (func,xArg,derOrder,ierr,eps) result(fDer)

    use explicitFunctionsModule, only : funcDerivative

    type(FunctionType), pointer       :: func
    real(dp)          , intent(in)    :: xArg
    integer           , intent(in)    :: derOrder
    integer           , intent(inout) :: ierr
    real(dp), optional, intent(in)    :: eps

    !! Local variables
    integer             :: lerr
    real(dp)            :: delta, fVal1, fVal2, fVal, fDer
    real(dp), parameter :: eps_p = 1.0e-10_dp

    !! --- Logic section ---

    if (derOrder == 1) then
       fDer = 1.0_dp
    else if (derOrder == 2) then
       fDer = 0.0_dp
    else
       fDer = FunctionValue(func,xArg,ierr)
       return
    end if

    if (.not.associated(func)) return
    if (func%id%baseId <= 0) return

    !! Check if explicit derivative is available, and compute that one if so
    call funcDerivative (func%id%baseId, derOrder, func%intParameters, &
         &               func%realParameters, xArg, fDer, lerr)
    if (lerr == 0) return

    if (present(eps)) then
       delta = max(eps,eps_p) ! Too small can cause numerical problems
    else
       delta = eps_p
    end if

    !! Numerical differentiation based on central difference
    fVal1 = FunctionValue(func,xArg-delta,ierr)
    fVal2 = FunctionValue(func,xArg+delta,ierr)
    if (derOrder == 2) then
       !! Central 2nd derivative based on algebraic Taylor series manipulation
       fVal = FunctionValue(func,xArg,ierr)
       fDer = (fVal1 - 2.0_dp*fVal2 + fVal)/(delta*delta)
    else
       fDer = 0.5_dp*(fVal2 - fVal1)/delta
    end if

  end function FunctionDerivative1


  !!============================================================================
  !> @brief Evaluates the derivative of a multi-argument function f(x,y,z,...).
  !>
  !> @param func Pointer to the function shape to evaluate the derivative of
  !> @param[in] xArg Function argument values to evaluate the derivative at
  !> @param[in] dVar Which argument to calculate derivative with respect to
  !> @param[in] derOrder Derivative order (0, 1 or 2)
  !> @param ierr Error flag
  !> @param[in] eps Size of argument interval for numerical derivatives
  !>
  !> @callgraph @callergraph
  !>
  !> @author Trond Arne Svidal
  !>
  !> @date Oct 2002

  function FunctionDerivative2 (func,xArg,dVar,derOrder,ierr,eps) result(fDer)

    use explicitFunctionsModule, only : funcDerivative
    use reportErrorModule      , only : internalError

    type(FunctionType), pointer       :: func
    real(dp)          , intent(in)    :: xArg(:)
    integer           , intent(in)    :: dVar, derOrder
    integer           , intent(inout) :: ierr
    real(dp), optional, intent(in)    :: eps

    !! Local variables
    integer  :: lerr, intPar(4)
    real(dp) :: fDer

    !! --- Logic section ---

    if (size(xArg) == 1) then
       fDer = FunctionDerivative1(func,xArg(1),derOrder,ierr,eps)
       return
    else if (derOrder == 1 .and. dVar == 1) then
       fDer = 1.0_dp
    else if (derOrder == 1 .or. derOrder == 2) then
       fDer = 0.0_dp
    else
       fDer = FunctionValue(func,xArg,ierr)
       return
    end if

    if (.not.associated(func)) return
    if (func%id%baseId <= 0) return

    !! Compute explicit derivative of the multi-argument function
    intPar(1:3) = func%intParameters(1:3)
    intPar(4) = dVar
    call funcDerivative (func%id%baseId, derOrder, intPar, &
         &               func%realParameters, xArg, fDer, lerr)
    if (lerr == 0) return

    ierr = internalError('FunctionDerivative2: Invalid multi-variable function')

  end function FunctionDerivative2


  !!============================================================================
  !> @brief Integrates a function shape over specified domain.
  !>
  !> @param func Pointer to the function shape to integrate
  !> @param[in] x0 Lower bound of integration domain
  !> @param[in] x1 Upper bound of integration domain
  !> @param[in] intOrder Integration order (0, 1 or 2)
  !> @param ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 6 Nov 1998

  function FunctionIntegral (func,x0,x1,intOrder,ierr) result(fInt)

    use explicitFunctionsModule, only : funcIntegral
    use reportErrorModule      , only : reportError, error_p

    type(FunctionType), pointer       :: func
    real(dp)          , intent(in)    :: x0, x1
    integer           , intent(in)    :: intOrder
    integer           , intent(inout) :: ierr

    !! Local variables
    real(dp), parameter :: epsZero_p = 1.0e-15_dp
    real(dp)            :: fInt, fInt0
    logical             :: identityFunc
    integer             :: lerr
    character(len=64)   :: errMsg

    !! --- Logic section ---

    identityFunc = .true.
    if (associated(func)) then
       if (func%id%baseId > 0) identityFunc = .false.
    end if

    if (identityFunc) then
       if (intOrder == 2) then
          fInt = (x1*x1*x1 - x0*x0*x0) / 6.0_dp
       else if (intOrder == 1) then
          fInt = (x1*x1 - x0*x0) / 2.0_dp
       else
          fInt = x1 - x0
       end if
       return
    end if

    lerr = ierr
    call funcIntegral (func%id%baseId, intOrder, func%intParameters, &
         &             func%realParameters, x1, fInt, ierr)
    if (ierr == lerr .and. abs(x0) > epsZero_p) then
       call funcIntegral (func%id%baseId, intOrder, func%intParameters, &
            &             func%realParameters, x0, fInt0, ierr)
       fInt = fInt - fInt0
    end if
    if (ierr == lerr) return

    if (abs(x0) > epsZero_p) then
       write(errMsg,"(' from x =',1PE12.5,' to x =',E12.5)") x0,x1
    else
       write(errMsg,"(' from x = 0.0 to x =',1PE12.5)") x1
    end if
    call reportError (error_p,'Unable to integrate '// &
         &            trim(getFunctionId(func))//errMsg)

  end function FunctionIntegral


  !!============================================================================
  !> @brief Checks whether a function can be integrated explicitly.
  !>
  !> @param func Pointer to a function shape object
  !> @param[in] order Integration order
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 6 Nov 1998

  function CanIntegrate (func,order)

    use explicitFunctionsModule, only : CanIntegrateExplicitly

    type(FunctionType), pointer    :: func
    integer           , intent(in) :: order
    logical                        :: CanIntegrate

    !! --- Logic section ---

    if (associated(func)) then
       CanIntegrate = CanIntegrateExplicitly(func%type,order)
    else
       CanIntegrate = .true. ! Assumed identity function
    end if

  end function CanIntegrate


  !!============================================================================
  !> @brief Checks if any functions have a control system variable as argument.
  !>
  !> @param engines All general functions in the model
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 27 Oct 2009

  function isCtrlSysUsed (engines)

    use SensorTypeModule, only : ENGINE_p, CONTROL_p

    type(EngineType), intent(in) :: engines(:)
    logical                      :: isCtrlSysUsed

    !! Local variables
    integer :: i

    !! --- Logic section ---

    isCtrlSysUsed = .true.
    do i = 1, size(engines)
       if (engines(i)%isUsed) then
          if (haveCtrlSysArg(engines(i))) return
       end if
    end do
    isCtrlSysUsed = .false.

  contains

    !> @brief Checks if a function have a control system variable as argument.
    recursive function haveCtrlSysArg (engine) result(have)
      type(EngineType), intent(in) :: engine
      logical :: have
      integer :: j
      have = .true.
      do j = 1, size(engine%args)
         if (associated(engine%args(j)%p)) then
            select case (engine%args(j)%p%type)
            case (CONTROL_p)
               return
            case (ENGINE_p)
               if (engine%args(j)%p%index(1) > 0) then
                  if (haveCtrlSysArg(engines(engine%args(j)%p%index(1)))) return
               end if
            end select
         end if
      end do
      have = .false.
    end function haveCtrlSysArg

  end function isCtrlSysUsed


  !!============================================================================
  !> @brief Initializes the ramp-up scaling function (smooth trajectory).
  !>
  !> @param[in] Tstart Start time of trajectory
  !> @param[in] Tend End time of function
  !> @param[in] Tlength Total trajectory time
  !> @param[in] Vmax Maximum speed during trajectory
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 19 Aug 2022

  subroutine initRampFunction (Tstart,Tend,Tlength,Vmax,ierr)

    use explicitFunctionsModule, only : SMOOTH_TRAJ_P
    use kindModule             , only : epsDiv0_p
    use reportErrorModule      , only : noteFileOnly_p, error_p
    use reportErrorModule      , only : reportError, allocationError

    real(dp), intent(in)  :: Tstart, Tend, Tlength, Vmax
    integer , intent(out) :: ierr

    !! Local variables
    character(len=96) :: msg

    !! --- Logic section ---

    if (Tlength < epsDiv0_p) then
       ierr = -1
    else
       ierr = 0
    end if
    if (Tend-Tstart - Tlength < -epsDiv0_p) ierr = ierr - 2
    if (Vmax < epsDiv0_p) ierr = ierr - 4
    if (Vmax*Tlength - 1.0_dp < epsDiv0_p) ierr = ierr - 8
    if (ierr < 0) then
       call reportError (error_p,'Invalid smooth trajectory parameters', &
            &            ierr=ierr, addString='initRampFunction')
       return
    end if

    if (.not. associated(ourRamp)) then
       allocate(ourRamp,stat=ierr)
    end if
    if (ierr == 0) then
       call nullifyFunction (ourRamp)
       allocate(ourRamp%intParameters(1),ourRamp%realParameters(5),stat=ierr)
    end if
    if (ierr /= 0) then
       ierr = allocationError('initRampFunction')
       return
    end if

    !! Define the ramp function (smooth trajectory with amplitude 1.0)
    ourRamp%id%baseId = 1
    ourRamp%id%userId = 1
    ourRamp%id%descr  = 'Smooth ramp function'
    ourRamp%intParameters(1)  = SMOOTH_TRAJ_p
    ourRamp%realParameters(1) = Tstart
    ourRamp%realParameters(2) = Tlength
    ourRamp%realParameters(3) = 2.0_dp*Vmax*Vmax / (Vmax*Tlength - 1.0_dp)
    ourRamp%realParameters(4) = Vmax
    ourRamp%realParameters(5) = Tend
    write(msg,600) ourRamp%realParameters
    call reportError (noteFileOnly_p,'Creating smooth trajectory ramp '// &
         &            'function with parameters:',msg)
600 format('Tstart =',F8.3,' Tlength =',F8.3,' Amax =',F8.3,' Vmax =',F8.3, &
         & ' Tend =',F8.3)

  end subroutine initRampFunction


  !!============================================================================
  !> @brief Evaluates the smooth trajectory ramp function.
  !>
  !> @param t Time to evaluate ramp function at, set to end of ramp on exit
  !> @param[in] derOr Derivative order (0 or 1)
  !> @param[out] rVal The ramp function value
  !> @param ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 19 Aug 2022

  subroutine getRampValue (t,derOr,rVal,ierr)

    real(dp), intent(inout) :: t
    integer , intent(in)    :: derOr
    real(dp), intent(out)   :: rVal
    integer , intent(inout) :: ierr

    !! Local variables
    integer :: lerr

    !! --- Logic section ---

    rVal = 1.0_dp
    lerr = ierr
    if (associated(ourRamp)) then
       if (t < ourRamp%realParameters(5)) then
          if (derOr > 0) then
             rval = FunctionDerivative(ourRamp,t,derOr,ierr)
          else
             rval = FunctionValue(ourRamp,t,ierr)
          end if
          if (ierr == lerr) t = ourRamp%realParameters(5)
       end if
    end if

  end subroutine getRampValue

end module FunctionTypeModule
