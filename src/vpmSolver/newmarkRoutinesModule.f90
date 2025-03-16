!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file newmarkRoutinesModule.f90
!>
!> @brief Subroutines for dynamic simulation based on Newmark time integration.

!!==============================================================================
!> @brief Module with subroutines for dynamic simulation in the time domain.
!>
!> @details This module contains subroutines implementing the Newmark time
!> integration driver in the FEDEM Dynamics Solver.
!> Both the Hilber-Hughes-Taylor (HHT-&alpha;) and the generalized-&alpha;
!> methods are supported, and which of these algorithms to use is controlled
!> via the third digit of the @a NewmarkFlag option (iopAlg):
!>  - iopAlg = 0: Use HHT-&alpha;
!>  - iopAlg = 1: Use HHT-&alpha; in FENRIS mode
!>  - iopAlg = 2: Use generalized-&alpha;
!>
!> The subroutine newmarkint() is the main driver for the time integration
!> process. It calculates the dynamic equilibrium state of the mechanism
!> through a Newton-Raphson iteration procedure. The time step loop itself
!> is in the caller of this subroutine. All other subroutines in this module
!> are utilities that are invoked by newmarkint() to carry out sub-tasks.
!>
!> All equation numbers below and elsewhere in the comments in the source code
!> are referring to the the Fedem R7.3 Theory Guide.
!> See the Sections 7.1 to 7.4 there for the derivation.

module NewmarkRoutinesModule

  implicit none

  private :: PredictAndUpdate


contains

#ifdef FT_DEBUG
  !> @cond NO_DOCUMENTATION
  subroutine writeArrays (dbg,arr1,name1,arr2,name2,arr3,name3,arr4,name4)
    use kindModule            , only : dp
    use manipMatrixModule     , only : writeObject
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_intValue
    integer                   , intent(in) :: dbg
    real(dp)        , optional, intent(in) :: arr1(:), arr2(:), arr3(:), arr4(:)
    character(len=*), optional, intent(in) :: name1, name2, name3, name4
    if (ffa_cmdlinearg_intValue('debug') == -99) then
       if (present(arr1) .and. present(name1)) then
          call writeObject (arr1,dbg,name1)
       end if
       if (present(arr2) .and. present(name2)) then
          call writeObject (arr2,dbg,name2)
       end if
       if (present(arr3) .and. present(name3)) then
          call writeObject (arr3,dbg,name3)
       end if
       if (present(arr4) .and. present(name4)) then
          call writeObject (arr4,dbg,name4)
       end if
    end if
  end subroutine writeArrays
  !> @endcond
#endif


  !!============================================================================
  !> @brief Starts a time step by computing predicted velocity and acceleration.
  !>
  !! @param[inout] sys System level model data
  !> @param[inout] motions All prescribed motions in the model
  !> @param[in]    iopAlg Algorithm option
  !> @param[out]   ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 3 Jun 2002

  subroutine PredictVelAcc (sys,motions,iopAlg,ierr)

    use SystemTypeModule    , only : SystemType
    use MotionTypeModule    , only : MotionType, dp
    use MotionRoutinesModule, only : updatePrescribedMotions
    use reportErrorModule   , only : reportError, debugFileOnly_p
#ifdef FT_DEBUG
    use dbgUnitsModule      , only : dbgSolve
#endif

    type(SystemType), intent(inout) :: sys
    type(MotionType), intent(inout) :: motions(:)
    integer         , intent(in)    :: iopAlg
    integer         , intent(out)   :: ierr

    !! Local variables
    real(dp) :: AA, AD, DA, DD

    !! --- Logic section ---

#ifdef FT_DEBUG
    write(dbgSolve,"(/'   --- in PredictVelAcc')")
#endif

    call updatePrescribedMotions (sys%time,sys%timeStep,sys%beta,sys%gamma, &
         &                        sys%urd,sys%urdd,motions,.false.,ierr)
    if (ierr < 0) then
       call reportError (debugFileOnly_p,'PredictVelAcc')
       return
    end if

    if (iopAlg > 0) then ! Generalized-alpha method
       !! Calculates  dk + urd  and  ak + urdd  for use in the force predictor
       DD = sys%gamma / sys%beta
       AA = 0.5_dp / sys%beta
    else ! HHT-alpha method
       !! Calculates dk and ak according to eqs. (7.19) and (7.21)
       DD = sys%gamma / sys%beta - 1.0_dp
       AA = 0.5_dp / sys%beta - 1.0_dp
    end if
    DA = (0.5_dp*sys%gamma/sys%beta - 1.0_dp)*sys%timeStep
    AD =  1.0_dp / (sys%beta*sys%timeStep)

    !! Predicted velocity and acceleration
    sys%dk = DD*sys%urd + DA*sys%urdd
    sys%ak = AD*sys%urd + AA*sys%urdd

    !! Initialize the total displacement increment for this time step
    sys%rinc = 0.0_dp

    !! Newmark incrementation scaling factors for velocity and acceleration
    !! see eq. (7.32)
    sys%sacc = 1.0_dp    / (sys%beta*sys%timeStep*sys%timeStep)
    sys%svel = sys%gamma / (sys%beta*sys%timeStep)

#ifdef FT_DEBUG
    call writeArrays (dbgSolve,sys%dk,'sys%dk',sys%ak,'sys%ak')
#endif

  end subroutine PredictVelAcc


  !!============================================================================
  !> @brief Increment the system velocity- and acceleration vectors.
  !>
  !> @param[inout] sys System level model data
  !>
  !> @details Also calculates tolerance variables for use in convergence check.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date Jun 2002

  subroutine IncVelAcc (sys)

    use SystemTypeModule, only : SystemType
#ifdef FT_DEBUG
    use dbgUnitsModule  , only : dbgSolve
#endif

    type(SystemType), intent(inout) :: sys

    !! --- Logic section ---

#ifdef FT_DEBUG
    write(dbgSolve,"(/'   --- in IncVelAcc')")
#endif

    sys%sinc(:,2) = sys%sinc(:,1)*sys%svel ! Iterative velocity correction
    sys%sinc(:,3) = sys%sinc(:,1)*sys%sacc ! Iterative acceleration correction

    sys%rinc = sys%rinc + sys%sinc(:,1)    ! Total displacement increment
    sys%urd  = sys%svel*sys%rinc - sys%dk  ! Updated velocity. eq. (7.20)
    sys%urdd = sys%sacc*sys%rinc - sys%ak  ! Updated acceleration, eq. (7.18)

#ifdef FT_DEBUG
    call writeArrays (dbgSolve,sys%rinc,'sys%rinc', &
         &            sys%urd,'sys%urd',sys%urdd,'sys%urdd')
#endif

  end subroutine IncVelAcc


  !!============================================================================
  !> @brief Updates all mechanism objects.
  !>
  !> @param[in]    sys System level model data
  !> @param[inout] mech Mechanism components of the model
  !> @param[in]    newVelocities Indicates whether new velocities are available.
  !> If not present, Tire and friction objects are not updated.
  !> @param[in]    updateCTI Tire update flag
  !> @param[out]   ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 7 Jun 2002

  subroutine UpdateMechanism (sys,mech,ierr,newVelocities,updateCTI)

    use KindModule                    , only : i8, maxInt_p
    use SystemTypeModule              , only : SystemType
    use MechanismTypeModule           , only : MechanismType
    use TriadTypeModule               , only : clearTriadForces
    use SupElRoutinesModule           , only : UpdateSupEls
    use SupElRoutinesModule           , only : UpdateSeaEnvironment
    use UserdefElRoutinesModule       , only : UpdateUDEs
    use BushingElementRoutinesModule  , only : UpdateBushingElements
    use ContactElementRoutinesModule  , only : UpdateContactElements
    use MasterSlaveJointRoutinesModule, only : UpdateJoints
    use SpringRoutinesModule          , only : UpdateSprings, UpdateSpringYields
    use DamperRoutinesModule          , only : UpdateDampers
    use ForceRoutinesModule           , only : UpdateExternalForces
    use MassRoutinesModule            , only : UpdateMasses
    use TireRoutinesModule            , only : UpdateTires, UpdateTireForces
    use FrictionRoutinesModule        , only : UpdateFrictions
    use profilerModule                , only : startTimer, stopTimer, upd_p
#ifdef FT_DEBUG
    use dbgUnitsModule                , only : dbgSolve, dbgTime, dbgIter
#endif
    use reportErrorModule             , only : reportError, debugFileOnly_p

    type(SystemType)   , intent(in)    :: sys
    type(MechanismType), intent(inout) :: mech
    integer            , intent(out)   :: ierr
    logical, optional  , intent(in)    :: newVelocities
    integer, optional  , intent(in)    :: updateCTI

    !! Local variables
    integer :: iStep
    logical :: updateTyres, newPositions

    !! --- Logic section ---

#ifdef FT_DEBUG
    write(dbgSolve,"(/'   --- in UpdateMechanism')")
    dbgTime = sys%time
    dbgIter = sys%nIterThisStep
#endif

    !! Note: The error flag is intent(inout) in all of the subroutines below
    !! (except for updateSupEls). It is therefore not necessary with a goto 915
    !! statement after each call. Instead, we contine to process all objects
    !! to find all that are in error.

    call startTimer (upd_p)

    if (sys%nStep >= int(maxInt_p,i8)) then
       iStep = 0 ! Avoid integer overflow, only used by FNV
       !! and possibly user-defined elements, so 99.9% harmless here
    else
       iStep = int(sys%nStep)+1
    end if

    if (present(newVelocities)) then
       updateTyres = newVelocities
    else
       updateTyres = .false.
    end if
    newPositions = present(updateCTI)

    if (newPositions) then
       call updateSeaEnvironment (mech%env,mech%triads,mech%sups, &
            &                     sys%time,sys%nIterThisStep,ierr)
       if (ierr < 0) goto 900
    else
       ierr = 0
    end if

    call clearTriadForces (mech%triads)
    call updateJoints (mech%joints,mech%motions,ierr)
    call updateSupEls (mech%sups,mech%supLoads,mech%env, &
         &             sys%beta,sys%gamma,sys%time,sys%timeStep, &
         &             iStep,sys%nIterThisStep,newPositions,ierr)
    call updateUDEs (mech%elms,mech%env,sys%time,sys%timeStep, &
         &           iStep,sys%nIterThisStep,ierr)

    !! Update all spring lengths and damper lengths/velocities first,
    !! since they may be measured by sensors used by other springs/dampers
    call updateSprings (mech%axialSprings,mech%joints,.false.,ierr, &
         &              updateLength=.true.,updateVar=.false.)
    call updateDampers (mech%dampers,sys%timeStep,.false.,ierr, &
         &              updateLV=.true.,updateVar=.false.)

    !! Now all quantities that may be measured by sensors should be up to date,
    !! so go on and update the remaining mechanism variables

    call updateSpringYields (mech%springYields,ierr)
    call updateSprings (mech%axialSprings,mech%joints,.false.,ierr, &
         &              updateLength=.false.,updateVar=.true.)
    call updateDampers (mech%dampers,sys%timeStep,.false.,ierr, &
         &              updateLV=.false.,updateVar=.true.)

    call updateBushingElements (mech%bElems,ierr)
    call updateContactElements (mech%cElems,sys%timeStep,ierr)

    !! Tires must be updated only after mechanism velocities have been updated
    if (updateTyres) call updateTires (sys,mech%tires,ierr,updateCTI)

    call updateExternalForces (mech%forces,2,ierr)
    call updateTireForces     (mech%tires,mech%gravity)
    call updateMasses         (mech%masses,mech%gravity,.false.,ierr)

    if (present(newVelocities)) then
       !! Frictions are updated after all triad force contributions have been
       !! computed. Superelements, axial spring and -dampers, external forces,
       !! additional masses and tires might have such contributions.
       call updateFrictions (mech%joints,mech%cElems,sys%timeStep, &
            &                sys%nIterThisStep,newVelocities,ierr)
    end if

    call stopTimer (upd_p)

900 if (ierr /= 0) call reportError (debugFileOnly_p,'updateMechanism')

  end subroutine UpdateMechanism


  !!============================================================================
  !> @brief Updates mechanism based on predicted velocity and acceleration.
  !>
  !> @param[in]    sam Data for managing system matrix assembly
  !> @param[in]    sys System level model data
  !> @param[inout] mech Mechanism components of the model
  !> @param[in]    iopAlg Algorithm option
  !> @param[out]   ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 7 Jun 2002

  subroutine PredictAndUpdate (sam,sys,mech,iopAlg,ierr)

    use SamModule                 , only : SamType
    use SystemTypeModule          , only : SystemType
    use MechanismTypeModule       , only : MechanismType
    use TriadTypeModule           , only : SetTriadsVelAcc
    use SupElRoutinesModule       , only : SetSupElsVelAcc
    use MasterSlaveJointTypeModule, only : SetJointsVelAcc
    use TireRoutinesModule        , only : SetTireStateVar
    use profilerModule            , only : startTimer, stopTimer, upd_p
#ifdef FT_DEBUG
    use dbgUnitsModule            , only : dbgSolve
#endif
    use reportErrorModule         , only : reportError, debugFileOnly_p

    type(SamType)      , intent(in)    :: sam
    type(SystemType)   , intent(in)    :: sys
    type(MechanismType), intent(inout) :: mech
    integer            , intent(in)    :: iopAlg
    integer            , intent(out)   :: ierr

    !! --- Logic section ---

#ifdef FT_DEBUG
    write(dbgSolve,"(/'   --- in PredictAndUpdate')")
#endif

    call startTimer (upd_p)

    !! Initiate tire state variables for next time step
    call SetTireStateVar (sys,mech%tires)

    !! Transfer the predicted velocities and accelerations down to object level
    call SetTriadsVelAcc (mech%triads,sys%dk,sys%ak)
    call SetSupElsVelAcc (mech%sups,sam,sys%dk,sys%ak)
    call SetJointsVelAcc (mech%joints,sys%dk,sys%ak)

    if (iopAlg > 0) then ! Generalized-alpha method
       !! Adjust dk and ak for use in the force corrector.
       !! In the predictor, urd and urdd were added so we subtract them again.
       sys%dk = sys%dk - sys%urd
       sys%ak = sys%ak - sys%urdd
    end if
#ifdef FT_DEBUG
    call writeArrays (dbgSolve,sys%dk,'sys%dk',sys%ak,'sys%ak')
#endif

    call stopTimer (upd_p)

    !! Now update all mechanism objects (including frictions and tires)
    call updateMechanism (sys,mech,ierr,.true.)

    if (ierr /= 0) call reportError (debugFileOnly_p,'PredictAndUpdate')

  end subroutine PredictAndUpdate


  !!============================================================================
  !> @brief Updates mechanism based on corrected velocity and acceleration.
  !>
  !> @param[in]    sam Data for managing system matrix assembly
  !> @param[inout] sys System level model data
  !> @param[inout] mech Mechanism components of the model
  !> @param[in]    delta Iterative solution vector correction
  !> @param[in]    useTotalInc If .true., use the total solution increment
  !> to update position matrices from the previous configuration
  !> @param[out]   ierr Error flag
  !>
  !> @details Increment all position, velocity and acceleration variables and
  !> bring all mechanism objects up-to-date with the updated configuration.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Karl Erik Thoresen               @date   Sep 1999
  !> @author Bjorn Haugen                     @date   Nov 2001
  !> @author Knut Morten Okstad               @date 4 Jun 2002

  subroutine IncAndUpdate (sam,sys,mech,delta,useTotalInc,ierr)

    use SamModule                     , only : SamType, dp
    use SystemTypeModule              , only : SystemType
    use MechanismTypeModule           , only : MechanismType
    use TriadTypeModule               , only : SetTriadsVelAcc, IncTriadsPos
    use SupElRoutinesModule           , only : SetSupElsVelAcc, IncSupElsGenDofs
    use MasterSlaveJointTypeModule    , only : SetJointsVelAcc
    use MasterSlaveJointRoutinesModule, only : IncJointsVar
    use TriadTypeModule               , only : TransSysVecToGlobal
    use solExtensionModule            , only : csExpand
    use profilerModule                , only : startTimer, stopTimer, upd_p
#ifdef FT_DEBUG
    use dbgUnitsModule                , only : dbgSolve
#endif
    use reportErrorModule             , only : reportError, debugFileOnly_p

    type(SamType)      , intent(in)    :: sam
    type(SystemType)   , intent(inout) :: sys
    type(MechanismType), intent(inout) :: mech
    real(dp)           , intent(in)    :: delta(:)
    logical            , intent(in)    :: useTotalInc
    integer            , intent(out)   :: ierr

    !! Local variables
    integer, parameter :: okCTI_p = 0

    !! --- Logic section ---

#ifdef FT_DEBUG
    write(dbgSolve,"(/'   --- in IncAndUpdate')")
    call writeArrays (dbgSolve,delta,'Solution increment')
#endif

    call startTimer (upd_p)

    !! Expand the equation-ordered delta to DOF-order
    if (sys%nIterThisStep == 0) then ! Predictor step
       call csExpand (sam,delta,sys%sinc(:,1))
    else ! Corrector step, replace the prescribed values, if any, by zero
       call csExpand (sam,delta,sys%sinc(:,1),S2=0.0_dp)
    end if
    call TransSysVecToGlobal (mech%triads,sys%sinc(:,1))

    !! Increment all position, velocity and acceleration variables
    call IncVelAcc (sys)
    if (useTotalInc) then
       call IncTriadsPos (mech%triads,sys%rinc,useTotalInc)
       call IncSupElsGenDofs (mech%sups,sys%rinc,useTotalInc)
       call IncJointsVar (mech%joints,sys%rinc,mech%motions,useTotalInc)
    else
       call IncTriadsPos (mech%triads,sys%sinc(:,1))
       call IncSupElsGenDofs (mech%sups,sys%sinc(:,1))
       call IncJointsVar (mech%joints,sys%sinc(:,1),mech%motions)
    end if

    !! Transfer the updated velocities and accelerations down to object level
    call SetTriadsVelAcc (mech%triads,sys%urd,sys%urdd)
    call SetSupElsVelAcc (mech%sups,sam,sys%urd,sys%urdd)
    call SetJointsVelAcc (mech%joints,sys%urd,sys%urdd)

    call stopTimer (upd_p)

    !! Now update all mechanism objects (including frictions and tires)
    call updateMechanism (sys,mech,ierr,.true.,okCTI_p)

    if (ierr /= 0) call reportError (debugFileOnly_p,'IncAndUpdate')

  end subroutine IncAndUpdate


  !!============================================================================
  !> @brief Advances the dynamics solution one step forward.
  !>
  !> @param[in]    sam Data for managing system matrix assembly
  !> @param[inout] sys System level model data
  !> @param[inout] mech Mechanism components of the model
  !> @param[inout] ctrl Control system data
  !> @param[in]    alpha Global structural damping parameters
  !> @param[in]    extRhs Additional external forces set by external process
  !> @param[inout] iop Control variable defining what to do (see below)
  !> @param[in]    NewmarkFlag Time integration control variable (see below)
  !> @param[inout] resFileFormat Flag for res-file output of convergence history
  !> @param[in]    lpu File unit number for res-file output
  !> @param[out]   ierr Error flag
  !>
  !> @details This subroutine performs different sub-tasks of the solution
  !> process, depending on the value of the control variable @a iop, as follows:
  !> - =  0 : Normal operation, solve the entire step with Newton-iterations
  !> - = -1 : Establish the linearized system of the first iteration
  !>         (the prediction step), and then return to the caller
  !> - = -2 : Same as -1, but also factorize the Newton matrix
  !> - = -3 : Same as 0, but this is the first time step after restart
  !> - =  2 : We are doing Newmark iterations (internal mode)
  !> - =  4 : Do one Newmark iteration on the existing linear system
  !> - =  6 : As 4, and then continue until equilibrium
  !> - = 16 : As 6, but include triangularization of existing equation system
  !>
  !> The @a NewmarkFlag variable is used to select among slightly different
  !> versions of the Newmark time integrator. It is a three-digit value (@e cba)
  !> where each of the digits has the following interpretation:
  !> - a = 1 : Use actual inertia force computed from the residual of the
  !>           previous increment in the right-hand-side calculation,
  !>           see eqs. (7.46), (7.53) and (7.55).
  !> - a = 2 : As for a = 1, but in addition the damping contribution to the
  !>           force predictor, eq. (7.53), is calculated consistently
  !>           for nonlinear dampers, see note ii) below.
  !> - b &gt; 0 : Use total solution increment in the configuration updates
  !> - c = 1 : Use an HHT-&alpha; algorithm equivalent to FENRIS, that is,
  !>           the force predictor is computed from eq. (7.52) and not (7.53).
  !> - c = 2 : Use generalized-&alpha; algorithm, see Section 7.4.3.
  !>
  !> @note i) If a = 0 (default setting), then the last computed inertia force
  !> (FIk) is used instead in the right-hand-side calculation in place of the
  !> eq. (7.46) quantity.
  !> @note ii) a = 2 is equivalent to a = 1 unless the model contains nonlinear
  !> dampers and/or hydrodynamic drag. If a < 2 and nonlinear dampers exist,
  !> the damping force contribution to the predictor is calculated by evaluating
  !> the nonlinear force-velocity curve for the predicted velocity (dk), instead
  !> of multiplying the tangent damping matrix by the velocity according to
  !> eq. (7.53). This is a simplification and not strictly consistent,
  !> but the option a < 2 is retained as default for backward compatibility.
  !> @note iii) The third digit c is equivalent to the variable @e iopAlg
  !> elsewhere in this documentation. c &gt; 0 also implies a = 2.
  !>
  !> @callgraph @callergraph
  !>
  !> @author MTHJ                              @date 1 Jun 1983
  !> @author Ole Ivar Sivertsen                @date 4 Feb 1989
  !> @author Karl Erik Thoresen                @date   Des 1998
  !> @author Bjorn Haugen                      @date   Mar 2002
  !> @author Knut Morten Okstad                @date 9 Jun 2002

  subroutine NewmarkInt (sam,sys,mech,ctrl,alpha,extRhs, &
       &                 iop,NewmarkFlag,resFileFormat,lpu,ierr)

    use KindModule               , only : dp, i8
    use SamModule                , only : SamType
    use SystemTypeModule         , only : SystemType
    use MechanismTypeModule      , only : MechanismType
    use ControlTypeModule        , only : ControlType
    use NormTypeModule           , only : iVecNorm_p, HasConverged
    use IdTypeModule             , only : StrId
    use ControlRoutinesModule    , only : IterateControlSystem
#ifdef FT_HAS_EXTCTRL
    use ExtCtrlSysRoutinesModule , only : IterateExtCtrlSys
#endif
    use SolverRoutinesModule     , only : reportEquationError, meqErr
    use SolverRoutinesModule     , only : reportNegativePivots
    use SolverRoutinesModule     , only : CalculateIterationNorms
    use SolverRoutinesModule     , only : iterationAccelerator
    use SolverRoutinesModule     , only : printConvergence, printSysSolveProcess
    use SolverRoutinesModule     , only : saveStep, restoreLastStep
    use TireRoutinesModule       , only : updateCTIresAtConvergence
    use ForceRoutinesModule      , only : updateExternalForces
    use EnvironmentTypeModule    , only : updateGravity
    use WindTurbineRoutinesModule, only : updateAeroForces
    use SupElRoutinesModule      , only : updateSupElDamping
    use EngineRoutinesModule     , only : preEvaluate, EngineValue
    use EngineRoutinesModule     , only : isPredictorStep
    use solExtensionModule       , only : csSolve
    use addInSysModule           , only : BuildNewtonMat, GetForceVectors
    use profilerModule           , only : startTimer, stopTimer, sol_p, upd_p
    use progressModule           , only : lterm
    use dbgUnitsModule           , only : dbgSolve, dbgShadowPos
    use explicitFunctionsModule  , only : dbgFunc
    use reportErrorModule        , only : reportError, debugFileOnly_p, error_p
    use reportErrorModule        , only : note_p, warning_p, warningFileOnly_p
    use FFaCmdLineArgInterface   , only : ffa_cmdlinearg_getdoubles
    use FFaCmdLineArgInterface   , only : ffa_cmdlinearg_getint
    use FFaCmdLineArgInterface   , only : ffa_cmdlinearg_isTrue

    type(SamType)      , intent(in)    :: sam
    type(SystemType)   , intent(inout) :: sys
    type(MechanismType), intent(inout) :: mech
    type(ControlType)  , intent(inout) :: ctrl
    real(dp)           , intent(in)    :: alpha(2)
    real(dp), optional , intent(in)    :: extRhs(:)
    integer            , intent(inout) :: iop, resFileFormat
    integer            , intent(in)    :: NewmarkFlag, lpu
    integer            , intent(out)   :: ierr

    !! Local variables
    integer, parameter :: nrhs_p = 1
    logical            :: useTotalInc, doTangentUpdate
    integer            :: stressStiffSkip, stressStiffSkipOnDivr
    integer            :: stopOnDivergence, cutbackNegPiv, iopAlg
    integer            :: predControl, corrControl, solveMode
    integer, save      :: nSeqNoUpdates = 0, nNegPivot = 0, ctrlSysMode = 0
    integer, save      :: nWarnPivot = 0, nWarnDiverg = 0, iCutStp = 0
    real(dp)           :: H, SSTIF, SMASS, SDAMP, EPS(1), saveIter(2)
    real(dp), save     :: tolUpdateFactor = 0.0_dp

    !! --- Logic section ---

    call ffa_cmdlinearg_getdoubles ('tolFactorize',EPS,1)
    call ffa_cmdlinearg_getdoubles ('saveIter',saveIter,2)
    call ffa_cmdlinearg_getint('stressStiffUpdateSkip',stressStiffSkip)
    call ffa_cmdlinearg_getint('stressStiffDivergSkip',stressStiffSkipOnDivr)
    call ffa_cmdlinearg_getint('stopOnDivergence',stopOnDivergence)
    call ffa_cmdlinearg_getint('cutbackNegPiv',cutbackNegPiv)
    call ffa_cmdlinearg_getint('predControl',predControl)
    call ffa_cmdlinearg_getint('corrControl',corrControl)
    if (predControl > 1 .or. corrControl > 1 .and. sys%nStep == 0_i8) then
       call reportError (note_p,'Predictor/Corrector control system:'// &
            &            StrId(10*predControl+corrControl))
    end if

    iopAlg      = mod(NewmarkFlag/100,10)
    useTotalInc = mod(NewmarkFlag/10,10) > 0

    if (iop < 1) then
       ctrlSysMode = 3

       !! Calculate initial mechanism forces (for numerical alpha damping)
       if (iop == -3) then
          iop = 0 ! Assuming FIactual is read from restart file
       else if (mod(NewmarkFlag,10) > 0 .or. iopAlg > 0) then
          sys%FIactual = sys%Qk - sys%FDk - sys%FSk
       else ! Ignore remaining residual from previous step
          sys%FIactual = sys%FIk
       end if
       if (iopAlg > 0) then ! Generalized-alpha method
          sys%Qprev  = sys%Qk
          sys%FIprev = sys%FIk
       end if

    end if

100 continue ! re-entry point after iteration cut-back
    if (iop < 1) then

       sys%time = sys%time + sys%timeStep
#ifdef FT_DEBUG
       write(dbgSolve,666) sys%nStep,sys%time,sys%timeStep
#endif
       if (dbgFunc > 0) write(dbgFunc,666) sys%nStep,sys%time,sys%timeStep
666    format(//3X,77('=') /'   --- in NewmarkInt, step, t, dt =',I4,1P2E12.5)

       if (ffa_cmdlinearg_isTrue('rampGravity')) then
          !! Check for ramping of gravity forces
          call updateGravity (mech%env,sys%time,ierr)
          if (ierr < 0) goto 915
       end if

       if (associated(mech%turbine)) then
          !! Calculate loads on wind turbine blades for the new time step
          call updateAeroForces (sys%time,mech%turbine,ierr)
          if (ierr /= 0) goto 915
       end if

       !! Initialize for new time step
       call preEvaluate (mech%engines,ierr)
       call updateExternalForces (mech%forces,1,ierr)
       call PredictVelAcc (sys,mech%motions,iopAlg,ierr)
       if (ierr < 0) goto 915

       call updateSupElDamping (mech%sups,mech%engines,alpha,ierr)
       if (ierr < 0) goto 915

       sys%nIterThisStep = 0
       sys%nUpdaThisStep = 0
       tolUpdateFactor = sys%tolUpdateFactor
       nSeqNoUpdates = 0
       nNegPivot = 0

       if (associated(sys%fxItEngine)) then
          !! A variable fixed number of iterations is used
          solveMode = sys%fixedIt
          sys%fixedIt = nint(EngineValue(sys%fxItEngine,ierr))
          if (ierr < 0) then
             goto 915
          else if (sys%fixedIt < 1 .and. solveMode > 0) then
             call reportError (note_p,'Using iteration tolerances')
          else if (sys%fixedIt /= solveMode) then
             call reportError (note_p,'Using fixed number of iterations'// &
                  &            StrId(sys%fixedIt))
          end if
       end if

    end if

    !! Newton matrix scaling factors for Newmark time integration
    if (iopAlg > 1) then ! Generalized-alpha method, eq. (7.64)
       SSTIF = sys%alpha_f
       SDAMP = sys%alpha_f*sys%svel
       SMASS = sys%alpha_m*sys%sacc
    else ! HHT-alpha method, eq. (7.51)
       SSTIF = sys%alpha + 1.0_dp
       SDAMP = sys%svel*SSTIF
       SMASS = sys%sacc
    end if

    Iterations: do ! --- Iteration loop
       if (iop > 2) goto 200 ! Iteration re-entry after external manipulation

       if (dbgFunc > 0) write(dbgFunc,"('   --- Iter = ',I3)") sys%nIterThisStep

       !!---------------------------------------------------------------
       !!        Decoupled integration, control module first
       !!---------------------------------------------------------------

       call IterateControlSystem (sys,ctrl,ctrlSysMode,sam%mpar,ierr)
       if (ierr < 0) goto 915

       if (sys%nIterThisStep == 0) then ! Predictor step

          isPredictorStep = iopAlg > 0 .or. mod(NewmarkFlag,10) == 2

          if (dbgShadowPos > 0) then
             write(dbgShadowPos,"('====== Time = ',1p,E12.3)") sys%time
          end if

          !! Update mechanism (including frictions) based on predicted
          !! velocity and acceleration, as well as the control output
          call PredictAndUpdate (sam,sys,mech,iopAlg,ierr)
          if (ierr < 0) goto 915

          if (predControl == 2) then
             call IterateControlSystem (sys,ctrl,ctrlSysMode,sam%mpar,ierr)
             if (ierr < 0) goto 915
             call UpdateMechanism (sys,mech,ierr)
             if (ierr < 0) goto 915
          end if

          call GetForceVectors (sys%FSk,sys%FDk,sys%FIk,sys%Qk,sys%RFk, &
               &                sam,mech,sys%nStep,alpha(2),ierr)
          if (ierr < 0) goto 915

          if (present(extRhs)) then
             sys%Qk = sys%Qk + extRhs ! Add forces set by external process
          end if

          !! Force predictor
          if (iopAlg == 2) then ! Generalized-alpha method, eq. (7.65)
             sys%residual = sys%alpha_f*(sys%Qk - sys%Qprev) &
                  &       + sys%alpha_f*sys%FDk + sys%FIactual &
                  &       + sys%alpha_m*sys%FIk - sys%FIprev
          else if (iopAlg == 1) then ! HHT-alpha method, eq. (7.52), FENRIS mode
             sys%residual = SSTIF*(sys%Qk - sys%Qprev + sys%FDk) &
                  &       + sys%FIk - sys%FIprev + sys%FIactual
          else ! HHT-alpha method, eq. (7.53)
             sys%residual = SSTIF*(sys%Qk - sys%FSk + sys%FDk) &
                  &       + sys%FIk - sys%alpha*sys%FIactual
          end if
#ifdef FT_DEBUG
          call writeArrays (dbgSolve,sys%residual,'Force predictor:')
#endif

          ctrlSysMode = 5

       else ! Corrector step

          isPredictorStep = .false.

          if (dbgShadowPos > 0) then
             write(dbgShadowPos,"('------ Iter = ',I3)") sys%nIterThisStep
          end if

          !! Update mechanism to transfer the control output onto the mechanism
          call UpdateMechanism (sys,mech,ierr)
          if (ierr < 0) goto 915

          if (corrControl == 2) then
             call IterateControlSystem (sys,ctrl,ctrlSysMode,sam%mpar,ierr)
             if (ierr < 0) goto 915
             call UpdateMechanism (sys,mech,ierr)
             if (ierr < 0) goto 915
          end if

          call GetForceVectors (sys%FSk,sys%FDk,sys%FIk,sys%Qk,sys%RFk, &
               &                sam,mech,sys%nStep,alpha(2),ierr)
          if (ierr < 0) goto 915

          if (present(extRhs)) then
             sys%Qk = sys%Qk + extRhs ! Add forces set by external process
          end if

          !! Force corrector
          if (iopAlg == 2) then ! Generalized-alpha method, eq. (7.67)
             sys%residual = sys%alpha_f*(sys%Qk - sys%FSk - sys%FDk) &
                  &       + (1.0_dp - sys%alpha_f)*sys%FIactual &
                  &       -           sys%alpha_m *sys%FIk &
                  &       - (1.0_dp - sys%alpha_m)*sys%FIprev
          else ! HHT-alpha method, eq. (7.55)
             sys%residual = SSTIF*(sys%Qk - sys%FSk - sys%FDk) &
                  &       - sys%FIk - sys%alpha*sys%FIactual
          end if
#ifdef FT_DEBUG
          call writeArrays (dbgSolve,sys%residual,'Force corrector:')
#endif

          if (sys%time > saveIter(1) .and. sys%time < saveIter(2)) then
             !! Save intermediate results for current iteration to assess
             !! the convergence, using a pseudo time t + dt*(IT-MAXIT)/MAXIT
             H = sys%time
             sys%time = H + sys%timeStep*(sys%nIterThisStep-sys%maxIt)/sys%maxIt
             call saveStep (sys,mech,ctrl,ierr)
             if (ierr < 0) goto 915
             sys%time = H
          end if

       end if

       sys%del = sys%residual ! Right-hand-side of the equation system
200    continue

       !! Determine if we need to update the tangent matrix in this iteration
       if (sys%nIterThisStep < max(1,sys%nUpdat)) then
          doTangentUpdate = .true. ! Always update in the first nUpdat iters.
       else if (nSeqNoUpdates >= sys%maxSequentialNoUpdate) then
          doTangentUpdate = .true. ! Too many sequentually iters. without update
       else if (sys%tolUpdateFactor <= 1.0_dp) then
          doTangentUpdate = .false. ! No update tolerance is given, tangent OK
       else if (HasConverged(sys%convergenceSet,tolUpdateFactor)) then
          doTangentUpdate = .false. ! Update tolerances were met, tangent OK
       else
          !! The tangent update tolerances were not met, so we need to update.
          !! If this is the first iteration with tangent update after a sequence
          !! of one or more iterations without updates, the tolerances are
          !! sharpened by a factor of 0.1 to ensure we now get closer to
          !! convergence before matrix updates is ommitted again.
          if (.not. doTangentUpdate) tolUpdateFactor = 0.1_dp*tolUpdateFactor
          doTangentUpdate = .true.
       end if

       if (doTangentUpdate .and. iop < 3) then

          !! Assemble the system Newton matrix

          if (sys%nIterThisStep == 0 .and. size(mech%motions) > 0) then
             !! In the predictor step, add contributions to the right-hand-side
             !! of the equation system from the prescribed motions, if any
             call BuildNewtonMat (sys%Nmat, smass, sdamp, sstif, mech, sam, &
                  &               sys%nIterThisStep, stressStiffSkip, &
                  &               sys%stressStiffIsOn(1), alpha(2), ierr, &
                  &               sys%del)
          else
             call BuildNewtonMat (sys%Nmat, smass, sdamp, sstif, mech, sam, &
                  &               sys%nIterThisStep, stressStiffSkip, &
                  &               sys%stressStiffIsOn(1), alpha(2), ierr)
          end if
          if (ierr < 0) goto 915

          if (iop == -1 .or. iop == 2) then
             iop = 4 ! Exit iteration loop for external manipulations
             return
          end if
       end if
       if (doTangentUpdate .and. (iop < 6 .or. iop == 16)) then

          solveMode = 1 ! Triangulate the updated Newton matrix
          call startTimer (sol_p)
          if (size(meqErr) > 1) then
             call csSolve (solveMode, 0, sys%Nmat, sys%del, nrhs_p, lpu, ierr, &
                  &        meqnInErr=meqErr, tolFactorize=EPS(1))
          else
             call csSolve (solveMode, 0, sys%Nmat, sys%del, nrhs_p, lpu, ierr, &
                  &        eqnInErr=meqErr(1), tolFactorize=EPS(1))
          end if
          call stopTimer (sol_p)
          if (ierr < 0) then
             goto 914
          else if (ierr > 0 .and. meqErr(1) > 0) then
             nNegPivot = nNegPivot + ierr
             if (nWarnPivot < 100) then
                nWarnPivot = nWarnPivot + 1
                if (resFileFormat == 1) resFileFormat = 2
                call reportNegativePivots (sam,mech%sups,mech%triads, &
                     &                     mech%joints,ierr,sys%time)
             end if
             if (ierr >= cutbackNegPiv .and. iop < 1) then
                if (sys%timeStep > sys%minInc .and. sys%cutbck(1) < 1.0_dp) then
                   call reportError (note_p,'Trying iteration cut-back')
                   call cutBack (.true.)
                   goto 100
                end if
             end if
          end if

          sys%nUpdates = sys%nUpdates + 1_i8
          sys%nUpdaThisStep = sys%nUpdaThisStep + 1
          nSeqNoUpdates = 0
          if (iop < 0) then
             iop = 6
             return
          end if
       else
          nSeqNoUpdates = nSeqNoUpdates + 1
       end if

       solveMode = 4 ! Calculate the displacement correction
       call startTimer (sol_p)
       call csSolve (solveMode, 0, sys%Nmat, sys%del, nrhs_p, lpu, ierr)
       call stopTimer (sol_p)
       if (ierr < 0) goto 915

       if (ffa_cmdlinearg_isTrue('lineSearch')) then
          call IterationAccelerator(sys%nIterThisStep,sys%residual,sys%del,ierr)
          if (ierr < 0) goto 915
       end if

       !! This also update frictions (the other updateMechanism calls do not)
       call IncAndUpdate (sam,sys,mech,sys%del,useTotalInc,ierr)
       if (ierr < 0) goto 915

#ifdef FT_HAS_EXTCTRL
       if (associated(ctrl%extCtrlSys)) then
          call IterateExtCtrlSys (ctrlSysMode,sys%time-sys%timeStep,sys%time, &
               &                  ctrl%extCtrlSys,ierr)
          if (ierr < 0) goto 915
       end if
#endif

       !! Calculate tolerances and present iteration norms
       call CalculateIterationNorms (sam,sys,sys%nIterThisStep,ierr)
       if (ierr < 0) goto 915

       !! Print iteration information
       if (resFileFormat /= 0) then
          call printConvergence (lpu,resFileFormat-1,sys%nStep+1_i8, &
               &                 sys%nIterThisStep,doTangentUpdate, &
               &                 sys%time,sys%timeStep,sys%del,sys%residual, &
               &                 sys%convergenceSet,sam,mech)
       end if
       if (sys%nStep < 9999999999999_i8) then
          write(LTERM,"(I14)",advance='NO') sys%nStep+1_i8
       else
          write(LTERM,"(14X)",advance='NO')
       end if
       write(LTERM,6001) sys%nIterThisStep,sys%time,sys%timeStep, &
            &            sys%convergenceSet%velNorms(iVecNorm_p)%value, &
            &            sys%convergenceSet%resNorms(iVecNorm_p)%value
6001   format(' |',I6,1X,1P,4(' |',D10.3),' |')

       if (resFileFormat > 1 .and. sys%nIterThisStep == 0) resFileFormat = 1

       if (.not. sys%convergenceSet%doingWell) then
          nWarnDiverg = nWarnDiverg + 1
          if (resFileFormat == 1) resFileFormat = 2
       end if

       if (dbgSolve > 0) then

          !! Print system vector components for current iteration to debug file
          call printSysSolveProcess (dbgSolve,sys%nIterThisStep,sys%time, &
               &                     sys%residual,sys%del, &
               &                     sys%FIk,sys%FDk,sys%FSk,sys%Qk, &
               &                     sys%FIactual,sam,mech)

       end if

       sys%nIterThisStep = sys%nIterThisStep + 1
       sys%nIter         = sys%nIter + 1_i8

       !! Check if fixed number of iterations
       if (sys%fixedIt > 0) then

          !! A fixed number of iterations has been requested
          if (sys%nIterThisStep == sys%fixedIt) exit

       else if (sys%nIterThisStep > sys%maxIt) then

          !! The maximum number of iterations has been reached
          if (iop > 0) then
             ierr = 3
             call MaxIterMessage (error_p,'Aborting')
             goto 915
          else if (nNegPivot > 10 .and. &
               &   stressStiffSkip < stressStiffSkipOnDivr) then
             !! We got many negative pivots during the iterations. This can be
             !! due to too much geometric stiffnesses. Try iteration cut-back
             !! without stress stiffening in the first few iterations.
             call MaxIterMessage (warning_p, &
                  'Trying cut-back without stress stiffening')
             stressStiffSkip = stressStiffSkipOnDivr
             call cutBack (.false.)
             goto 100
          else if (sys%timeStep > sys%minInc .and. sys%cutbck(1) < 1.0_dp) then
             !! Try iteration cut-back with reduced time step size
             call MaxIterMessage (warning_p, &
                  'Trying cut-back with reduced step size')
             call cutBack (.true.)
             goto 100
          else if (ffa_cmdlinearg_isTrue('continueAfterMaxIter')) then
             call MaxIterMessage (warning_p,'Continuing')
             exit
          else
             ierr = 1
             call MaxIterMessage (error_p,'Aborting')
             goto 915
          end if

       else if (sys%nIterThisStep >= max(2,sys%minIt)) then

          !! The minimum number of iterations has now been performed,
          !! now check for convergence in the solution norms
          if (HasConverged(sys%convergenceSet)) exit

       end if

       if (iop == 6 .or. iop == 16) then
          iop = 0
       else if (iop == 4) then
          iop = 2
       end if
    end do Iterations
    if (iop > 0) iop = 0

    !! Time step converged, save CTI state variables
    call startTimer (upd_p)
    call updateCTIresAtConvergence (sys,mech%tires,ierr)
    call stopTimer (upd_p)
    if (ierr < 0) goto 915

    !! Integration step accepted, update mechanism variables
    call GetForceVectors (sys%FSk,sys%FDk,sys%FIk,sys%Qk,sys%RFk, &
         &                sam,mech,sys%nStep,alpha(2),ierr)
    if (ierr < 0) goto 915

    if (present(extRhs)) then
       sys%Qk = sys%Qk + extRhs ! Add forces set by external process
    end if

    !! Try increase time step size again after successful cut-back
    if (sys%cutbck(2) < 1.0_dp) then
       if (iCutStp < sys%nCutStp) then
          iCutStp = iCutStp + 1
       else
          iCutStp = 0
          sys%cutbck(2) = min(1.0_dp,sys%cutbck(2)/sys%cutbck(1))
       end if
    end if

    !! Check if we have had too many warnings on possibly divergence
    if (stopOnDivergence < 1 .or. nWarnDiverg < stopOnDivergence) return

    ierr = 2
    call reportError (error_p,'The maximum number of iterations allowed '// &
         &            'with poor convergence has been reached.')
    goto 915


    !!  E R R O R   R E T U R N

914 continue
    call reportEquationError (sam,mech%sups,mech%triads,mech%joints,ierr)
    if (ierr == -3 .and. iop < 1 .and.ffa_cmdlinearg_isTrue('cutbackSing')) then
       if (nNegPivot > 10 .and. stressStiffSkip < stressStiffSkipOnDivr) then
          !! We got many negative pivots during the iterations. This can be
          !! due to too much geometric stiffnesses. Try iteration cut-back
          !! without stress stiffening in the first few iterations.
          call reportError (note_p, &
               'Trying iteration cut-back without stress stiffening')
          stressStiffSkip = stressStiffSkipOnDivr
          call cutBack (.false.)
          goto 100
       else if (sys%timeStep > sys%minInc .and. sys%cutbck(1) < 1.0_dp) then
          call reportError (note_p,'Trying iteration cut-back')
          call cutBack (.true.)
          goto 100
       end if
    end if

915 continue
    call reportError (debugFileOnly_p,'NewmarkInt')

  contains

    !> @brief Prints a message when the max number of iterations is reached.
    subroutine MaxIterMessage (msgCode,action)
      integer         , intent(in) :: msgCode
      character(len=*), intent(in) :: action
      integer         , save       :: nWarnMaxIter = 0
      integer :: msgUnit
      msgUnit = msgCode
      if (msgCode == warning_p) then
         nWarnMaxIter = nWarnMaxIter + 1
         if (nWarnMaxIter > 10) msgUnit = warningFileOnly_p
      end if
      call reportError (msgUnit, &
           'Max. number of iterations reached without convergence - '//action)
      if (nWarnMaxIter == 10) call reportError (note_p, &
           'No further "Max. number of iterations" messages will be given.', &
           'See fedem_solver.res file for the complete history.')
    end subroutine MaxIterMessage

    !> @brief Resets the current time and time step in case of cut-back.
    subroutine cutBack (reduceStepSize)
      logical, intent(in) :: reduceStepSize
      iCutStp = 0
      ctrlSysMode = 4
      sys%time = sys%time - sys%timeStep
      if (reduceStepSize) then
         sys%cutbck(2) = sys%cutbck(2)*sys%cutbck(1)
         sys%timeStep = max(sys%minInc,sys%cutbck(1)*sys%timeStep)
      end if
      call restoreLastStep (sam,sys,mech)
    end subroutine cutBack

  end subroutine NewmarkInt

end module NewmarkRoutinesModule
