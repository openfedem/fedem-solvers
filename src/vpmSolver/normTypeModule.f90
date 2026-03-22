!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file normTypeModule.f90
!>
!> @brief Convergence norm data containers.

!!==============================================================================
!> @brief Module with data types representing convergence norm objects.
!>
!> @details The module also contains subroutines for accessing the norm data.

module NormTypeModule

  use kindModule, only : dp

  implicit none

  integer, parameter :: nNormTypes_p = 4  !< Total number of norm types
  integer, parameter :: iVecNorm_p = 1, & !< Index for L_2 displacement norm
       &                iInfTra_p  = 2, & !< Index for L_inf translation norm
       &                iInfRot_p  = 3, & !< Index for L_inf rotation norm
       &                iInfGen_p  = 4   !< Index for L_inf generalized DOF norm


  !> @brief Data type representing a convergence check definition.
  type TestItemType
     character(len=12) :: title !< Name to be used in res-file headings
     integer  :: code      !< ON/OFF switch for this norm
     integer  :: nIterIncr !< Number of iterations with increasing norm value
     real(dp) :: value     !< The current value of this convergence norm
     real(dp) :: tolerance !< The actual convergence tolerance for this norm
     real(dp) :: absTol    !< Absolute convergence tolerance for this norm
     real(dp) :: relTol    !< Velocity proportional tolerance for this norm
  end type TestItemType

  !> @brief Data type representing a set of convergence checks.
  type TestSetType
     logical :: doingWell    !< If .true., the convergence is good
     logical :: isActive(3)  !< Activation flags (dis/vel/acc, force, energy)
     integer :: startMonitor !< Search for worst DOFs id iter &gt; startMonitor
     integer , pointer  :: worstDOFs(:,:) !< Worst DOFs when convergence is slow
     real(dp), pointer  :: worstEnerg(:)  !< Worst energy when slow convergence
     type(TestItemType) :: disNorms(nNormTypes_p) !< Displacement norms
     type(TestItemType) :: velNorms(nNormTypes_p) !< Velocity norms
     type(TestItemType) :: accNorms(nNormTypes_p) !< Acceleration norms
     type(TestItemType) :: resNorms(nNormTypes_p) !< Force residual norms
     type(TestItemType) :: energyNorms(2)         !< Energy norms
  end type TestSetType

  private :: InitTolerance


contains

  !!============================================================================
  !> @brief Initializes a set of convergence checks object.
  !>
  !> @param[out] convSet The normtypemodule::testsettype object to initialize
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 30 Jun 2004

  subroutine NullifyConvSet (convSet)

    type(TestSetType), intent(out) :: convSet

    !! --- Logic section ---

    convSet%doingWell = .true.
    convSet%isActive = .false.
    convSet%startMonitor = 0
    nullify(convSet%worstDOFs)
    nullify(convSet%worstEnerg)

  end subroutine NullifyConvSet


  !!============================================================================
  !> @brief Deallocates a set of convergence checks object.
  !>
  !> @param convSet The normtypemodule::testsettype object to deallocate
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Jan 2017

  subroutine DeallocateConvSet (convSet)

    type(TestSetType), intent(inout) :: convSet

    !! --- Logic section ---

    if (associated(convSet%worstDOFs))  deallocate(convSet%worstDOFs)
    if (associated(convSet%worstEnerg)) deallocate(convSet%worstEnerg)

    call nullifyConvSet (convSet)

  end subroutine DeallocateConvSet


  !!============================================================================
  !> @brief Initialize the convergence checks from command-line options.
  !>
  !> @param convSet Set of convergence checks
  !> @param[in] maxIt Maximum number of iterations per time step
  !> @param[in] monWorst Number of DOFs to monitor on poor convergence
  !> @param[in] monIter Number of iterations to monitor before @a maxit
  !> @param[in] relTol Velocity proportional tolerance on velocity corrections
  !> @param[out] err Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 30 Jun 2004

  subroutine InitConvChecks (convSet,maxIt,monWorst,monIter,relTol,err)

    use reportErrorModule, only : allocationError, reportError, note_p,warning_p

    type(TestSetType), intent(inout) :: convSet
    integer          , intent(in)    :: maxIt, monWorst, monIter
    real(dp)         , intent(in)    :: relTol
    integer          , intent(out)   :: err

    !! Local variables
    integer :: i

    !! --- Logic section ---

    call initTolerance ('tolDispNorm',' L2(displ)',convSet%disNorms(iVecNorm_p))
    call initTolerance ('tolDispTra','Inf(traDis)',convSet%disNorms(iInfTra_p))
    call initTolerance ('tolDispRot','Inf(angDis)',convSet%disNorms(iInfRot_p))
    call initTolerance ('tolDispGen','Inf(genDis)',convSet%disNorms(iInfGen_p))
    call initTolerance ('tolVelNorm','  L2(veloc)',convSet%velNorms(iVecNorm_p))
    call initTolerance ('tolVelTra' ,'Inf(traVel)',convSet%velNorms(iInfTra_p))
    call initTolerance ('tolVelRot' ,'Inf(angVel)',convSet%velNorms(iInfRot_p))
    call initTolerance ('tolVelGen' ,'Inf(genVel)',convSet%velNorms(iInfGen_p))
    call initTolerance ('tolAccNorm','  L2(accel)',convSet%accNorms(iVecNorm_p))
    call initTolerance ('tolAccTra' ,'Inf(traAcc)',convSet%accNorms(iInfTra_p))
    call initTolerance ('tolAccRot' ,'Inf(angAcc)',convSet%accNorms(iInfRot_p))
    call initTolerance ('tolAccGen' ,'Inf(genAcc)',convSet%accNorms(iInfGen_p))
    call initTolerance ('tolResNorm','    L2(res)',convSet%resNorms(iVecNorm_p))
    call initTolerance ('tolResTra' ,'Inf(traRes)',convSet%resNorms(iInfTra_p))
    call initTolerance ('tolResRot' ,'Inf(angRes)',convSet%resNorms(iInfRot_p))
    call initTolerance ('tolResGen' ,'Inf(genRes)',convSet%resNorms(iInfGen_p))
    call initTolerance ('tolEnerMax','E-norm(max)',convSet%energyNorms(1))
    call initTolerance ('tolEnerSum','E-norm(sum)',convSet%energyNorms(2))

    do i = 1, nNormTypes_p
       if (convSet%disNorms(i)%code > 0) convSet%isActive(1) = .true.
       if (convSet%velNorms(i)%code > 0) convSet%isActive(1) = .true.
       if (convSet%accNorms(i)%code > 0) convSet%isActive(1) = .true.
       if (convSet%resNorms(i)%code > 0) convSet%isActive(2) = .true.
    end do
    if (convSet%energyNorms(1)%code > 0) convSet%isActive(3) = .true.
    if (convSet%energyNorms(2)%code > 0) convSet%isActive(3) = .true.

    if (relTol > 0.0_dp .and. convSet%velNorms(iVecNorm_p)%code > 0) then
       convSet%velNorms(iVecNorm_p)%relTol = relTol
       call reportError (note_p,'Using velocity proportional tolerance '// &
            &            'for the scaled vector norm of velocity corrections')
    else if (.not. any(convSet%isActive)) then
       call reportError (warning_p,'No convergence tests defined')
    end if

    convSet%startMonitor = maxIt - monIter
    allocate(convSet%worstDOFs(monWorst,3), &
         &   convSet%worstEnerg(monWorst), stat=err)
    if (err /= 0) err = allocationError('InitConvChecks')

  end subroutine InitConvChecks


  !!============================================================================
  !> @brief Initializes a convergence check with data from command-line.
  !>
  !> @param[in] tolType Command-line option to get convergence tolerance from
  !> @param[in] normName Name of solution norm defining the convergence check
  !> @param testItem The convergence check object to initialize
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 29 Apr 2003

  subroutine InitTolerance (tolType,normName,testItem)

    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getdouble

    character(len=*)  , intent(in)    :: tolType, normName
    type(TestItemType), intent(inout) :: testItem

    !! Local variables
    real(dp) :: value

    !! --- Logic section ---

    testItem%title     = normName
    testItem%nIterIncr = 0
    testItem%code      = -1
    testItem%value     = 0.0_dp
    testItem%tolerance = 0.0_dp
    testItem%absTol    = 0.0_dp
    testItem%relTol    = 0.0_dp

    call ffa_cmdlinearg_getdouble (tolType,value)
    if (value > 0.0_dp) then
       testItem%absTol = value
       testItem%code   = 2 ! allOf
    else if (value < 0.0_dp) then
       testItem%absTol = -value
       testItem%code   = 1 ! oneOf
    end if

  end subroutine InitTolerance


  !!============================================================================
  !> @brief Checks if the state of the convergence checks set is converged.
  !>
  !> @param convSet Set of convergence checks
  !> @param[in] factor_opt Optional tolerance scaling factor
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 28 Mar 2003

  function HasConverged (convSet,factor_opt)

    type(TestSetType) , intent(in) :: convSet
    real(dp), optional, intent(in) :: factor_opt
    logical                        :: hasConverged

    !! Local variables
    logical :: testsExist(2), oneOfTest, allOfTest
    integer :: iNorm
    integer, parameter :: oneOf_p = 1, allOf_p = 2


    !! --- Logic section ---

    testsExist = .false.
    oneOfTest  = .false.
    allOfTest  = .true.

    do iNorm = 1, nNormTypes_p
       call TestItem (convSet%disNorms(iNorm))
       call TestItem (convSet%velNorms(iNorm))
       call TestItem (convSet%accNorms(iNorm))
       call TestItem (convSet%resNorms(iNorm))
    end do
    call TestItem (convSet%energyNorms(1))
    call TestItem (convSet%energyNorms(2))

    hasConverged = .true.
    if (testsExist(oneOf_p)) hasConverged = hasConverged .and. oneOfTest
    if (testsExist(allOf_p)) hasConverged = hasConverged .and. allOfTest

  contains

    !> @brief Evaluates the state of a convergence check.
    subroutine TestItem (item)
      type(TestItemType), intent(in) :: item
      real(dp) :: tol
      if (present(factor_opt)) then
         tol = item%tolerance*factor_opt
      else
         tol = item%tolerance
      end if
      select case (item%code)
      case(oneOf_p); if (item%value <  tol) oneOfTest = .true.
      case(allOf_p); if (item%value >= tol) allOfTest = .false.
      case default ; return
      end select
      testsExist(item%code) = .true.
    end subroutine TestItem

  end function HasConverged


  !!============================================================================
  !> @brief Check if a given norm is showing sign of possible divergence.
  !>
  !> @param testItem The convergence check object to check for dovergence
  !> @param[in] value Current norm value
  !> @param mayDiverge If .true., a possible divergence is detected
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 19 Mar 2004

  subroutine checkDivergence (testItem,value,mayDiverge)

    type(TestItemType), intent(inout) :: testItem
    real(dp)          , intent(in)    :: value
    logical           , intent(inout) :: mayDiverge

    !! Local variables
    integer, parameter :: maxIterIncrBeforeWarning_p = 2

    !! --- Logic section ---

    if (testItem%code <= 0) return ! Only check norms used in convergence checks

    if (value >= testItem%value .and. testItem%value >= testItem%tolerance) then
       testItem%nIterIncr = testItem%nIterIncr + 1
    else
       testItem%nIterIncr = 0
    end if

    if (testItem%nIterIncr >= maxIterIncrBeforeWarning_p) mayDiverge = .true.

  end subroutine checkDivergence

end module NormTypeModule
