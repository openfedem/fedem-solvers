!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file staticEquilibriumModule.f90
!>
!> @brief Subroutines for quasi-static time-domain simulation.

!!==============================================================================
!> @brief Module with subroutines for quasi-static time-domain simulation.

module staticEquilibriumModule

  implicit none

  private

  public :: staticEquilibrium, staticInt


contains

  !!============================================================================
  !> @brief Performs static equilibrium iterations to find start configuration.
  !>
  !> @param sam Data for managing system matrix assembly
  !> @param sys System level model data
  !> @param mech Mechanism components of the model
  !> @param ctrl Control system data
  !> @param[in] xinp External function values for initial state
  !> @param num_xinp Number of external function values to extract/extracted
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @details This is a rewritten version of the old QUAEQL subroutine by
  !> OIS and others from 1989. It solves for the initial configuration where
  !> all external loads (including gravity) is evaluated at the starting time,
  !> through Newton-Raphson iterations. The purpose is to avoid initial noise
  !> in the subsequent dynamic simulation due to large initial transients.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 30 May 2002

  subroutine staticEquilibrium (sam,sys,mech,ctrl,xinp,num_xinp,lpu,ierr)

    use KindModule              , only : dp, i8
    use SamModule               , only : SamType
    use SystemTypeModule        , only : SystemType
    use MechanismTypeModule     , only : MechanismType
    use ControlTypeModule       , only : ControlType
    use SysMatrixTypeModule     , only : writeObject
    use MotionRoutinesModule    , only : updatePrescribedMotions
    use SolverRoutinesModule    , only : reportEquationError, meqErr
    use SolverRoutinesModule    , only : reportNegativePivots
    use SolverRoutinesModule    , only : printSysSolveProcess, saveStep
    use ControlRoutinesModule   , only : IterateControlSystem
#ifdef FT_HAS_EXTCTRL
    use ExtCtrlSysRoutinesModule, only : InitExtCtrlSys
#endif
    use WindTurbineRoutinesModule,only : updateAeroForces, addInAeroForces
    use NormRoutinesModule      , only : ScaledNorm
    use solExtensionModule      , only : csSolve
    use addInSysModule          , only : BuildStiffMat, GetStaticForceVectors
    use profilerModule          , only : startTimer, stopTimer, sol_p
    use progressModule          , only : lterm, writeProgress
    use manipMatrixModule       , only : writeObject
    use dbgUnitsModule          , only : dbgSolve
    use reportErrorModule       , only : reportError, debugFileOnly_p, error_p
    use reportErrorModule       , only : warning_p, empty_p
    use FFaCmdLineArgInterface  , only : ffa_cmdlinearg_getdoubles
    use FFaCmdLineArgInterface  , only : ffa_cmdlinearg_intValue
    use FFaCmdLineArgInterface  , only : ffa_cmdlinearg_isTrue
    use FiDeviceFunctionInterface,only : fidf_initExtFunc

    type(SamType)      , intent(inout) :: sam
    type(SystemType)   , intent(inout) :: sys
    type(MechanismType), intent(inout) :: mech
    type(ControlType)  , intent(inout) :: ctrl
    real(dp), optional , intent(in)    :: xinp(*)
    integer , optional , intent(inout) :: num_xinp
    integer            , intent(in)    :: lpu
    integer            , intent(out)   :: ierr

    !! Local variables
    logical  :: useAdditionalBCs, notWarned
    integer  :: doIterations, ctrlSysMode, iFactorize, iSolve, iSing, initAero
    integer  :: IT, NINC, dLPU
    real(dp) :: DELN, EPS(2), factor, saveIter(2)

    character(len=128) :: errMsg
    integer, parameter :: stressStiffUpdateSkip = 0

    !! --- Logic section ---

#ifdef FT_DEBUG
    write(dbgSolve,"(//3X,77('=')/'   --- in staticEquilibrium')")
#endif
    if (ffa_cmdlinearg_intValue('debug') /= -99) then
       dLPU = 0
    else if (dbgSolve > 0) then
       dLPU = dbgSolve
    else
       dLPU = lpu
    end if

    ctrlSysMode = 1 ! Initiate the control system
    call IterateControlSystem (sys,ctrl,ctrlSysMode,sam%mpar,ierr)
    if (ierr < 0) goto 915

#ifdef FT_HAS_EXTCTRL
    if (associated(ctrl%extCtrlSys)) then
       !! Initiate the external control systems
       call InitExtCtrlSys (sys%tStart,sys%timeStep,ctrl%extCtrlSys,ierr)
       if (ierr < 0) goto 915
    end if
#endif

    !! Check whether we actually want to do the equilibrium iterations
    if (sys%equLim > 0.0_dp .and. sam%neq > 0) then
       if (sys%equTol > 0.0_dp) then
          doIterations = 2 ! Perform non-linear iterations
          write(lpu,"(/'  Starting initial static equilibrium analysis')")
       else
          doIterations = 1 ! Perform a linear static analysis (single iteration)
          write(lpu,"(/'  Doing linear static analysis')")
       end if
    else
       doIterations = 0
    end if

    if (doIterations > 0) then
       !! Set initial values for the prescribed motions, but only the prescribed
       !! displacements will be accounted for during the initial iterations
       call updatePrescribedMotions (sys%time,0.0_dp,sys%beta,sys%gamma, &
            &                        sys%urd, sys%urdd, mech%motions, &
            &                        doIterations == 1, ierr)
       if (ierr < 0) goto 915
    end if

    !! Update all mechanism objects (initial configuration).
    !! Note(1): This call will also evaluate state-dependent general functions
    !! if they are used to define external loads, prescribed motions, etc.
    !! External functions are however not yet initialized at this point if they
    !! take their values from a file, so they will all evaluate to zero here.
    call staticIncAndUpdate (sam,sys,mech,ierr,doIterations=doIterations>0)
    if (ierr < 0) goto 915

    if (doIterations == 0) then

       if (associated(mech%turbine)) then
          !! Calculate the initial aerodynamic loads on wind turbine blades
          call updateAeroForces (sys%time,mech%turbine,ierr)
          if (ierr < 0) goto 915
       end if

       !! Quasi-static equilibrium iterations are not wanted.
       !! Same comment as the Note(1) above applies here as well.
       call finalStaticUpdate (sam,sys,mech,ierr)
       if (ierr < 0) goto 915

       if (present(num_xinp)) num_xinp = 0
       return

    else if (doIterations == 2 .and. sys%equLim < 1.01_dp*sys%equTol) then

       ierr = -1
       write(errMsg,1011) sys%equLim,sys%equTol
1011   format('The iteration step size limit',1PE12.5,' is too small ', &
            & 'compared to the equilibrium convergence tolerance',E12.5)
       call reportError (error_p,errMsg)
       goto 915

    else if (sys%maxIt < min(sys%minIt,1)) then

       ierr = -1
       write(errMsg,1012) sys%maxIt,sys%minIt
1012   format('maxIt =',i6,'  minIt =',i6)
       call reportError (error_p, &
            'The Max and Min number of iterations are incompatible',errMsg)
       goto 915

    end if

    !! We have go for static equilibrium iterations ...
    if (doIterations == 2) then
       call writeProgress (' --> STATIC EQUILIBRIUM ITERATIONS')
       write(lterm,100)
    else
       call writeProgress (' --> LINEAR STATIC ANALYSIS')
    end if
100 format(/15X,'| It.no. |Max stepsize| Displ. norm| Tol. Limit'/10X,62('-') )

    call ffa_cmdlinearg_getdoubles ('tolFactorize',EPS,2)
    call ffa_cmdlinearg_getdoubles ('saveIter',saveIter,2)

    !! Check singularity handling options
    if (ffa_cmdlinearg_isTrue('fixSingEqInit')) then
       iSing = -1 ! Fix all singular (or nearly singular) DOFs to zero
    else if (ffa_cmdlinearg_isTrue('releaseSingEqInit')) then
       iSing =  1 ! Release all singular DOFs (by adding a small stiffness)
    else
       iSing =  0 ! Abort on singularities
    end if
    notWarned = .true.

    !! Check whether the additional BC's should be enforced
    useAdditionalBCs = sam%ndof1 < sam%neq
    if (useAdditionalBCs) then
       iFactorize = 2
       iSolve = 6
    else
       iFactorize = 1
       iSolve = 4
    end if

    !! Check for inclusion of aerodynamic loads, if any
    if (.not. associated(mech%turbine) .or. doIterations == 1) then
       initAero = -1
    else if (ffa_cmdlinearg_isTrue('initEqAD')) then
       initAero =  1
    else
       initAero =  0
    end if

    !! Initialize the external function values.
    !! TODO: Maybe this should be placed before the first mechanism update call?
    !!       See the Note(1) above.
    if (present(xinp) .and. present(num_xinp)) then
       call fidf_initExtFunc (xinp,num_xinp,ierr)
       if (ierr < 0) then
          call reportError (error_p,'Failed to initialize external functions', &
               &            'The provided input array is too short.',ierr=ierr)
       end if
       num_xinp = ierr
    else
       call fidf_initExtFunc (EPS,0,ierr)
    end if

200 continue
    LoadSteps: do NINC = 0, sys%maxIt-1 ! --- Load increment loop

       Iterations: do IT = 1, sys%maxIt ! --- Equilibrium iteration loop

          call GetStaticForceVectors (sys%FSk,sys%Qk,sys%RFk,sam,mech,ierr)
          if (initAero == 2) then
             call addInAeroForces (sys%Qk,sys%RFk,mech%turbine,sam,ierr)
          end if
          if (ierr < 0) goto 915

          !! Calculate force residual
          call DCOPY (sam%neq,sys%Qk(1),1,sys%residual(1),1)
          call DAXPY (sam%neq,-1.0_dp,sys%FSk(1),1,sys%residual(1),1)
          !! Right-hand-side of equation system
          call DCOPY (sam%neq,sys%residual(1),1,sys%del(1),1)

          if (sys%nIterThisStep == 0 .and. size(mech%motions) > 0) then
             !! In the first iteration, add contributions to the right-hand-side
             !! of the equation system from the prescribed displacements, if any
             call BuildStiffMat (sys%Nmat,mech,sam,it,sys%stressStiffIsOn(2), &
                  &              stressStiffUpdateSkip,ierr,sys%del)
             !! Store the prescribed motion forces in FIactual for printing
             if (dbgSolve > 0) sys%FIactual = sys%del - sys%residual
          else
             call BuildStiffMat (sys%Nmat,mech,sam,it,sys%stressStiffIsOn(2), &
                  &              stressStiffUpdateSkip,ierr)
             if (dbgSolve > 0) sys%FIactual = 0.0_dp
          end if
          if (ierr < 0) goto 915

          if (sys%nUpdates < 5_i8 .and. dLPU > 0) then
             write(dLPU,"(//' >>> DUMP OF SYSTEM STIFFNESS MATRIX <<<<<')")
             write(dLPU,*) '    ndof1       =', sam%ndof1
             write(dLPU,*) '    neq         =', sam%neq
             call writeObject (sys%Nmat,sam%mpar,dLPU)
             write(dLPU,"('')")
             call writeObject (sys%del,dLPU,'     Right-hand-side vector')
          end if

          !! Factor the tangent stiffness matrix
          call startTimer (sol_p)
          if (size(meqErr) > 1) then
             call csSolve (iFactorize, iSing, sys%Nmat, sys%del, 1, lpu, ierr, &
                  &        meqnInErr=meqErr, tolFactorize=EPS(2), &
                  &        neq1=sam%ndof1)
          else
             call csSolve (iFactorize, iSing, sys%Nmat, sys%del, 1, lpu, ierr, &
                  &        eqnInErr=meqErr(1), tolFactorize=EPS(2), &
                  &        neq1=sam%ndof1)
          end if
          call stopTimer (sol_p)
          if (ierr < 0) then
             goto 914
          else if (ierr > 0 .and. iSing == 0) then
             call reportNegativePivots (sam,mech%sups,mech%triads,mech%joints, &
                  &                     ierr,sys%time)
          else if (ierr > 0 .and. notWarned) then
             notWarned = .false.
             call reportEquationError (sam,mech%sups,mech%triads,mech%joints, &
                  -3, 'Analysis continues with modified stiffness matrix')
          end if

          !! Back and forward solve with respect to Right-Hand-Side
          if (useAdditionalBCs) sys%del(sam%ndof1+1:) = 0.0_dp
          call startTimer (sol_p)
          call csSolve (iSolve, iSing, sys%Nmat, sys%del, 1, lpu, ierr, &
               &        neq1=sam%ndof1)
          call stopTimer (sol_p)
          if (ierr < 0) goto 915

          if (dLPU > 0) then
             write(dLPU,"('')")
             call writeObject (sys%del,dLPU,'     Solution vector')
          end if

          DELN = ScaledNorm(sam,sys%del,sys%wDisp,2,.true.,ierr)
          if (ierr < 0) goto 915

          if (doIterations == 1) then
             write(LPU  ,1003) DELN
             write(LTERM,6003) DELN
          else if (DELN > sys%equLim) then
             factor = sys%equLim/DELN ! Scale down the solution increment
             sys%del(1:sam%ndof1) = sys%del(1:sam%ndof1) * factor
             write(LPU  ,1001) NINC, IT, sys%equLim, DELN, sys%equTol, factor
             write(LTERM,6001) NINC, IT, sys%equLim, DELN, sys%equTol, factor
          else
             write(LPU  ,1002) NINC, IT, sys%equLim, DELN, sys%equTol
             write(LTERM,6002) NINC, IT, sys%equLim, DELN, sys%equTol
          end if

1001      format('     NINC =',I4,'  IT =',I3, &
               & ' MAXINC =',1P,D12.4,'  DELN =',D12.4,'  TOL =',D12.4, &
               & '  (inc. scaled down with',D10.2,')')
1002      format('     NINC =',I4,'  IT =',I3, &
               & ' MAXINC =',1P,D12.4,'  DELN =',D12.4,'  TOL =',D12.4)
1003      format('     DELN =',1P,D12.5/)
6001      format(I14,' |',I6,1X,1P,3(' | ',D10.3), &
               & '  (inc. scaled down with',D10.2,')')
6002      format(I14,' |',I6,1X,1P,3(' | ',D10.3))
6003      format(/15X,'Displacement norm: ',1PD12.5)

          !! Increment the static solution and update all mechanism variables
          call staticIncAndUpdate (sam,sys,mech,ierr,sys%del)
          if (ierr < 0) goto 915

          if (dbgSolve > 0) then

             !! Print out system vector components for current iteration
             call printSysSolveProcess (dbgSolve,sys%nIterThisStep,sys%time, &
                  &                     sys%residual,sys%del, &
                  &                     sys%FIk,sys%FDk,sys%FSk,sys%Qk, &
                  &                     sys%FIactual,sam,mech)

          end if

          sys%nUpdates      = sys%nUpdates + 1_i8
          sys%nIter         = sys%nIter + 1_i8
          sys%nIterThisStep = sys%nIterThisStep + 1

          if (doIterations == 1) goto 500 ! Linear analysis, no checking
          if (DELN < sys%equTol .and. IT >= sys%minIt) goto 500 ! Convergence

          if (sys%nIterThisStep < sys%maxIt) then
             if (sys%time > saveIter(1) .and. sys%time < saveIter(2)) then
                !! Save intermediate results for current iteration to assess
                !! the convergence, using a pseudo time t + dt*(IT-MAXIT)/MAXIT
                factor   = sys%time
                sys%time = sys%time &
                     &   + sys%timeStep*(sys%nIterThisStep-sys%maxIt)/sys%maxIt
                call saveStep (sys,mech,ctrl,ierr)
                if (ierr < 0) goto 915
                sys%time = factor
             end if
          end if

          if (DELN > sys%equLim) goto 400 ! Restart the iteration loop

       end do Iterations

       !! The equilibrium iterations did not fully converge, but accept
       write(errMsg,1013) DELN/sys%equTol
1013   format('Current value of DELN/TOL = ',1PE12.5, ' (should be < 1)')
       call reportError (warning_p,'The initial static equilibrium itera'// &
            'tions did not fully converge.',errMsg,'The mechanism may be '// &
            'too far from equilibrium for a dynamics simulation to be started.')
       goto 500

400    continue ! Try restart the iteration loop with a down-scaled load

    end do LoadSteps

    ierr = sys%maxIt
    write(errMsg,1014) sys%maxIt
1014 format('The initial static equilibrium iterations did not converge', &
         ' after',i4,' load iterations')
    call reportError (error_p,errMsg, &
         'Check that the model has sufficient boundary conditions '// &
         'for a quasi-static analysis')
    goto 915

500 continue
    initAero = initAero + 1
    if (initAero == 1 .or. initAero == 2) then
       !! Calculate the initial aerodynamic loads on wind turbine blades
       write(lpu,"(/'  Including aerodynamic forces')")
       call updateAeroForces (sys%time,mech%turbine,ierr)
       if (ierr < 0) goto 915
    end if
    if (initAero == 2) goto 200 ! Reiterate with aerodynamic loads included

    ctrlSysMode = 2 ! Find the initial steady-state of the control system
    call IterateControlSystem (sys,ctrl,ctrlSysMode,sam%mpar,ierr)
    if (ierr < 0) goto 915

    call updatePrescribedMotions (mech%motions,ierr)
    if (ierr < 0) goto 915

    !! Do final configuration update (initial velocity dependent stuff, etc.)
    call finalStaticUpdate (sam,sys,mech,ierr)
    if (ierr < 0) goto 915

    return

    !!  E R R O R   R E T U R N

914 continue
    call reportEquationError (sam,mech%sups,mech%triads,mech%joints,ierr)
    call reportError (error_p,'Initial static equilibrium analysis failed.', &
         'Check that the model has sufficient boundary conditions '// &
         'for a quasi-static analysis.')
    if (ierr == -3) then
       call reportError (empty_p,'Consider Additional Boundary Conditions.')
    end if

915 continue
    call reportError (debugFileOnly_p,'staticEquilibrium')

  end subroutine staticEquilibrium


  !!============================================================================
  !> @brief Advances the quasi-static solution one time step.
  !>
  !> @param sam Data for managing system matrix assembly
  !> @param sys System level model data
  !> @param mech Mechanism components of the model
  !> @param ctrl Control system data
  !> @param[in] extRhs Additional external forces set by some external process
  !> @param iop Control variable defining what to do (see below)
  !> @param resfileFormat Flag for res-file output of convergence history
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @details This subroutine performs different sub-tasks of the solution
  !> process, depending on the value of the control variable iop, as follows:
  !> - =  0 : Normal operation, solve the entire step with Newton-iterations
  !> - = -1 : Establish the linearized system of the first iteration,
  !>          and then return to the calling module
  !> - = -2 : Same as -1, but also factorize the stiffness matrix
  !> - =  1 : We are doing quasi-static iterations (internal mode)
  !> - =  3 : Do one quasi-static iteration on the existing linear system
  !> - =  5 : As 3, and then continue until equilibrium
  !> - = 15 : As 5, but include triangularization of existing equation system
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 30 May 2002

  subroutine staticInt (sam,sys,mech,ctrl,extRhs,iop,resFileFormat,lpu,ierr)

    use KindModule               , only : dp, i8
    use SamModule                , only : SamType
    use SystemTypeModule         , only : SystemType
    use MechanismTypeModule      , only : MechanismType
    use ControlTypeModule        , only : ControlType
    use NormTypeModule           , only : iVecNorm_p, HasConverged
    use IdTypeModule             , only : StrId
    use SysMatrixTypeModule      , only : writeObject
    use ManipMatrixModule        , only : writeObject
    use EnvironmentTypeModule    , only : updateGravity
    use MotionRoutinesModule     , only : updatePrescribedMotions
    use SolverRoutinesModule     , only : reportEquationError, meqErr
    use SolverRoutinesModule     , only : reportNegativePivots
    use SolverRoutinesModule     , only : CalculateIterationNorms
    use SolverRoutinesModule     , only : printConvergence, printSysSolveProcess
    use SolverRoutinesModule     , only : saveStep, restoreLastStep
    use EngineRoutinesModule     , only : EngineValue
    use solExtensionModule       , only : csSolve
    use addInSysModule           , only : BuildStiffMat, GetStaticForceVectors
    use profilerModule           , only : startTimer, stopTimer, sol_p
    use progressModule           , only : lterm
    use fileUtilitiesModule      , only : getDBGfile
    use dbgUnitsModule           , only : dbgSolve
    use reportErrorModule        , only : reportError, debugFileOnly_p, error_p
    use reportErrorModule        , only : note_p, empty_p
    use FFaCmdLineArgInterface   , only : ffa_cmdlinearg_getdoubles
    use FFaCmdLineArgInterface   , only : ffa_cmdlinearg_getint
    use FFaCmdLineArgInterface   , only : ffa_cmdlinearg_isTrue

    type(SamType)      , intent(in)    :: sam
    type(SystemType)   , intent(inout) :: sys
    type(MechanismType), intent(inout) :: mech
    type(ControlType)  , intent(inout) :: ctrl
    real(dp), optional , intent(in)    :: extRhs(:)
    integer            , intent(inout) :: iop, resFileFormat
    integer            , intent(in)    :: lpu
    integer            , intent(out)   :: ierr

    !! Local variables
    integer, parameter :: nrhs_p = 1
    logical            :: doTangentUpdate, pureLinear
    integer, save      :: nSeqNoUpdates = 0, nWarnDiverg = 0, iCutStp = 0
    integer            :: solveMode, stopOnDivergence, cutbackNegPiv
    integer            :: stressStiffUpdateSkip, nFricSkip, iprint, dLPU, lpr
    real(dp)           :: H, EPS(1), saveIter(2)
    real(dp), save     :: tolUpdateFactor = 0.0_dp

    !! --- Logic section ---

    dLPU = lpu
#ifdef FT_DEBUG
    if (iop > 1) dLPU = getDBGfile(35,'inverse.dbg')
    write(dbgSolve,10) sys%time,iop
10  format(//3X,77('=') /'   --- in StaticInt t, iop =',1PE12.5,I3)
#endif

    pureLinear = sys%fixedIt == -1
    call ffa_cmdlinearg_getdoubles ('tolFactorize',EPS,1)
    call ffa_cmdlinearg_getdoubles ('saveIter',saveIter,2)
    call ffa_cmdlinearg_getint('stressStiffUpdateSkip',stressStiffUpdateSkip)
    call ffa_cmdlinearg_getint('nFricSkipStatic',nFricSkip)
    call ffa_cmdlinearg_getint('cutbackNegPiv',cutbackNegPiv)
    call ffa_cmdlinearg_getint('debug',iprint)
    if (iprint > 1) then
       lpr = lpu
    else
       lpr = 0
    end if

    if (iop < 1) then
       !! Reduce the load step size if we are doing iteration cut-back
       h = max(sys%minInc,sys%timeStep*sys%cutbck(2))
       if (h < sys%timeStep) sys%timeStep = h
    end if
100 continue ! re-entry point after iteration cut-back
    if (iop < 1) then
       sys%time = sys%time + sys%timeStep

       if (sys%fixedIt < 0) then
          !! Restore initial configuration in case of no incrementation
          call restoreLastStep (sam,sys,mech)
          sys%nIterThisStep = 0
       end if

       if (ffa_cmdlinearg_isTrue('rampGravity')) then
          !! Check for ramping of gravity forces
          call updateGravity (mech%env,sys%time,ierr)
          if (ierr < 0) goto 915
       end if

       !! Initialize the prescribed motions for new load step
       call updatePrescribedMotions (sys%time,sys%timeStep,sys%beta,sys%gamma, &
            &                        sys%urd, sys%urdd, mech%motions, &
            &                        sys%fixedIt < 0, lpu=lpr, ierr=ierr)
       if (ierr < 0) goto 915

       !! Update all mechanism variables for the new load step
       call staticIncAndUpdate (sam,sys,mech,ierr,noCorotUpd=pureLinear)
       if (ierr < 0) goto 915

       sys%nIterThisStep = 0
       sys%nUpdaThisStep = 0
       tolUpdateFactor = sys%tolUpdateFactor
       nSeqNoUpdates = 0

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

    doTangentUpdate = .true.
    Iterations: do ! --- Iteration loop

       if (iop < 3 .and. sam%neq > 0) then

          !! Compute system force vectors
          call GetStaticForceVectors (sys%FSk,sys%Qk,sys%RFk,sam,mech,ierr)
          if (ierr < 0) goto 915

          if (present(extRhs)) then
             !! Add forces set by external process
             call DAXPY (sam%neq,1.0_dp,extRhs(1),1,sys%Qk(1),1)
          end if
          !! Get the external load vector
          !! (pure linear or multi-loadcase nonlinear analysis)
          call DCOPY (sam%neq,sys%Qk(1),1,sys%residual(1),1)
          if (sys%fixedIt >= 0 .or. sys%nIterThisStep > 0) then
             !! Calculate the force residual (nonlinear analysis)
             call DAXPY (sam%neq,-1.0_dp,sys%FSk(1),1,sys%residual(1),1)
          end if
          !! Right-hand-side of equation system
          call DCOPY (sam%neq,sys%residual(1),1,sys%del(1),1)

       end if
       if (sys%nIterThisStep > 0 .and. iop < 3) then
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

       !! Determine if we need to update the tangent matrix in this iteration
       if (sam%neq < 1) then
          doTangentUpdate = .false. ! No unknowns, all DOFs are prescribed
       else if (sys%nIterThisStep < max(1,sys%nUpdat)) then
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

          !! Assemble the system stiffness matrix

          if (sys%nIterThisStep == 0 .and. size(mech%motions) > 0) then
             !! In the first iteration, add contributions to the right-hand-side
             !! of the equation system from the prescribed displacements, if any
             call BuildStiffMat (sys%Nmat, mech, sam, &
                  &              sys%nIterThisStep, sys%stressStiffIsOn(1), &
                  &              stressStiffUpdateSkip, ierr, sys%del)
             !! Store the prescribed motion forces in FIactual for printing
             if (dbgSolve > 0) sys%FIactual = sys%del - sys%residual
          else
             call BuildStiffMat (sys%Nmat, mech, sam, &
                  &              sys%nIterThisStep, sys%stressStiffIsOn(1), &
                  &              stressStiffUpdateSkip, ierr)
             if (dbgSolve > 0) sys%FIactual = 0.0_dp
          end if
          if (ierr < 0) goto 915

          if (sys%nUpdates < 5_i8 .and. iprint == -99) then
             write(lpu,"(//' >>> DUMP OF SYSTEM STIFFNESS MATRIX <<<<<')")
             write(lpu,*) '    ndof1       =', sam%ndof1
             write(lpu,*) '    neq         =', sam%neq
             call writeObject (sys%Nmat,sam%mpar,lpu)
          end if

          if (abs(iop) == 1) then
             iop = 3
             return
          end if
       end if
       if (doTangentUpdate .and. (iop < 5 .or. iop == 15)) then

          solveMode = 1 ! Triangulate the updated stiffness matrix
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
          else if (ierr > 0) then
             if (resFileFormat == 1) resFileFormat = 2
             call reportNegativePivots (sam,mech%sups,mech%triads,mech%joints, &
                  &                     ierr,sys%time)
             if (ierr >= cutbackNegPiv .and. iop < 1 .and. sys%fixedIt >= 0)then
                if (sys%timeStep > sys%minInc .and. sys%cutbck(1) < 1.0_dp) then
                   call reportError (note_p,'Trying iteration cut-back')
                   call cutBack ()
                   goto 100
                end if
             end if
          end if

          sys%nUpdates = sys%nUpdates + 1_i8
          sys%nUpdaThisStep = sys%nUpdaThisStep + 1
          nSeqNoUpdates = 0
          if (iop < 0) then
             iop = 5
             return
          end if
       else
          nSeqNoUpdates = nSeqNoUpdates + 1
       end if
       if (sam%neq > 0) then

          if (dLPU /= lpu) then
             call writeObject (reshape(sys%del,(/sam%neq,1/)),dLPU, &
                  &            'System right-hand-side vector')
          else if (iprint == -99) then
             call writeObject (sys%del,lpu,'System right-hand-side vector')
          end if

          solveMode = 4 ! Calculate the displacement correction
          call startTimer (sol_p)
          call csSolve (solveMode, 0, sys%Nmat, sys%del, nrhs_p, lpu, ierr)
          call stopTimer (sol_p)
          if (ierr < 0) goto 915

          if (dLPU /= lpu) then
             call writeObject (reshape(sys%del,(/sam%neq,1/)),dLPU, &
                  &            'System solution vector')
          else if (iprint == -99) then
             call writeObject (sys%del,lpu,'System solution vector')
          end if

       end if

       !! Increment the static solution and update all mechanism variables
       if (sys%nIterThisStep < nFricSkip) then
          call staticIncAndUpdate (sam,sys,mech,ierr,sys%del)
       else
          !! Also update frictions (if any), using pseudo velocities
          call staticIncAndUpdate (sam,sys,mech,ierr,sys%del, &
               &                   noCorotUpd=pureLinear,doFrictions=.false.)
       end if
       if (ierr < 0) goto 915

       !! Calculate tolerances and present iteration norms
       call CalculateIterationNorms (sam,sys,sys%nIterThisStep,ierr)
       if (ierr < 0) goto 915

       if (resFileFormat /= 0) then
          !! Print iteration information
          call printConvergence (lpu,resFileFormat-1, &
               &                 sys%nStep+1_i8,sys%nIterThisStep, &
               &                 doTangentUpdate .or. sam%neq < 1, &
               &                 sys%time,sys%timeStep,sys%del,sys%residual, &
               &                 sys%convergenceSet,sam,mech)
       end if
       if (sys%nStep < 9999999999999_i8) then
          write(LTERM,"(I14)",advance='NO') sys%nStep+1_i8
       else
          write(LTERM,"(14X)",advance='NO')
       end if
       write(LTERM,6001) sys%nIterThisStep,sys%time,sys%timeStep, &
            &            sys%convergenceSet%disNorms(iVecNorm_p)%value, &
            &            sys%convergenceSet%resNorms(iVecNorm_p)%value
6001   format(' |',I6,1X,1P,4(' |',D10.3),' |')

       if (resFileFormat > 1 .and. sys%nIterThisStep == 0) resFileFormat = 1

       if (.not. sys%convergenceSet%doingWell) then
          nWarnDiverg = nWarnDiverg + 1
          if (resFileFormat == 1) resFileFormat = 2
       end if

       if (dbgSolve > 0) then

          !! Print out system vector components for current iteration
          call printSysSolveProcess (dbgSolve,sys%nIterThisStep,sys%time, &
               &                     sys%residual,sys%del, &
               &                     sys%FIk,sys%FDk,sys%FSk,sys%Qk, &
               &                     sys%FIactual,sam,mech)

       end if

       sys%nIterThisStep = sys%nIterThisStep + 1
       sys%nIter         = sys%nIter + 1_i8

       !! Check if fixed number of iterations
       if (pureLinear .or. sys%fixedIt > 0) then

          !! A fixed number of iterations has been requested
          if (sys%nIterThisStep == abs(sys%fixedIt)) exit

       else if (sys%nIterThisStep > sys%maxIt) then

          call reportError (error_p, &
               'The maximum number of iterations has been reached')
          if (iop > 0) then
             ierr = 3
             goto 915
          else if (sys%timeStep > sys%minInc .and. sys%cutbck(1) < 1.0_dp) then
             !! Try iteration cut-back with reduced load step size
             call reportError (empty_p,'Trying cut-back with reduced step size')
             call cutBack ()
             goto 100
          else if (ffa_cmdlinearg_isTrue('continueAfterMaxIter')) then
             call reportError (empty_p,'Trying to continue with next load step')
             exit
          else
             ierr = 1
             goto 915
          end if

       else if (sys%nIterThisStep >= max(2,sys%minIt)) then

          !! The minimum number of iterations has now been performed,
          !! now check for convergence in the solution norms
          if (HasConverged(sys%convergenceSet)) exit

       end if

       if (iop == 5 .or. iop == 15) then
          iop = 0
       else if (iop == 3) then
          iop = 1
       end if
    end do Iterations
    if (iop > 0) iop = 0

    !! Integration step accepted, update mechanism variables
    call GetStaticForceVectors (sys%FSk,sys%Qk,sys%RFk,sam,mech,ierr)
    if (ierr < 0) goto 915

    if (present(extRhs)) then
       sys%Qk = sys%Qk + extRhs ! Add forces set by external process
    end if

    if (sys%cutbck(2) < 1.0_dp) then
       !! Try increase load step size again after successful cut-back
       if (iCutStp < sys%nCutStp) then
          iCutStp = iCutStp + 1
       else
          iCutStp = 0
          sys%cutbck(2) = min(1.0_dp,sys%cutbck(2)/sys%cutbck(1))
       end if
    end if

    !! Check if we have had too many warnings on possibly divergence
    call ffa_cmdlinearg_getint('stopOnDivergence',stopOnDivergence)
    if (stopOnDivergence < 1 .or. nWarnDiverg < stopOnDivergence) return

    ierr = 2
    call reportError (error_p,'The maximum number of iterations allowed '// &
         &            'with poor convergence has been reached.')
    goto 915


    !!  E R R O R   R E T U R N

914 continue
    call reportEquationError (sam,mech%sups,mech%triads,mech%joints,ierr)
    if (ierr == -3 .and. iop < 1 .and.ffa_cmdlinearg_isTrue('cutbackSing')) then
       call reportError (note_p,'Trying iteration cut-back')
       call cutBack ()
       goto 100
    else
       call reportError (empty_p,'Check that the model has sufficient '// &
            &            'boundary conditions for a quasi-static analysis.')
    end if

915 continue
    call reportError (debugFileOnly_p,'StaticInt')

  contains

    !> @brief Resets the current time and time step in case of cut-back.
    subroutine cutBack ()
      iCutStp = 0
      sys%time = sys%time - sys%timeStep
      sys%cutbck(2) = sys%cutbck(2)*sys%cutbck(1)
      sys%timeStep = max(sys%minInc,sys%cutbck(1)*sys%timeStep)
      call restoreLastStep (sam,sys,mech)
    end subroutine cutBack

  end subroutine staticInt


  !!============================================================================
  !> @brief Increments all position variables and updates the configuration.
  !>
  !> @param[in] sam Data for managing system matrix assembly
  !> @param sys System level model data
  !> @param mech Mechanism components of the model
  !> @param[out] ierr Error flag
  !> @param[in] delta Iterative solution increment
  !> @param[in] noCorotUpd If .true., the co-rotated frames are not updated
  !> @param[in] doIterations If .false., no initial equilibrium iterations
  !> (used by the tire module only)
  !> @param[in] doFrictions If not present, skip friction update.
  !> If .false., use pseudo-velocities during friction update.
  !>
  !> @details All mechanism objects are brought up-to-date with the updated
  !> configuration while the velocity and acceleration variables are assumed
  !> identically zero.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 30 May 2002

  subroutine staticIncAndUpdate (sam,sys,mech,ierr,delta, &
       &                         noCorotUpd,doIterations,doFrictions)

    use KindModule                    , only : dp, i8, maxInt_p
    use SamModule                     , only : SamType
    use SystemTypeModule              , only : SystemType
    use MechanismTypeModule           , only : MechanismType
    use TriadTypeModule               , only : clearTriadForces
    use TriadTypeModule               , only : TransSysVecToGlobal, IncTriadsPos
    use SupElRoutinesModule           , only : SetSupElsVelAcc, IncSupElsGenDofs
    use SupElRoutinesModule           , only : updateSeaEnvironment
    use SupElRoutinesModule           , only : updateSupElsStatic
    use UserdefElRoutinesModule       , only : updateUDEs
    use BushingElementRoutinesModule  , only : updateBushingElements
    use ContactElementRoutinesModule  , only : updateContactElements
    use MasterSlaveJointRoutinesModule, only : updateJoints, IncJointsVar
    use SpringRoutinesModule          , only : updateSprings, updateSpringYields
    use DamperRoutinesModule          , only : updateDampers
    use ForceRoutinesModule           , only : updateExternalForces
    use MassRoutinesModule            , only : updateMasses
    use TireRoutinesModule            , only : updateTires, updateTireForces
    use FrictionRoutinesModule        , only : updateFrictions
    use solExtensionModule            , only : csExpand
#ifdef FT_DEBUG
    use dbgUnitsModule                , only : dbgSolve, dbgTime, dbgIter
#endif
    use profilerModule                , only : startTimer, stopTimer, upd_p
    use reportErrorModule             , only : reportError, debugFileOnly_p

    type(SamType)      , intent(in)    :: sam
    type(SystemType)   , intent(inout) :: sys
    type(MechanismType), intent(inout) :: mech
    integer            , intent(out)   :: ierr
    real(dp), optional , intent(in)    :: delta(:)
    logical , optional , intent(in)    :: noCorotUpd, doIterations, doFrictions

    !! Local variables
    integer :: iStep, linInc

    !! --- Logic section ---

#ifdef FT_DEBUG
    if (present(delta)) then
       write(dbgSolve,"(/'   --- in staticIncAndUpdate')",advance='NO')
    else
       write(dbgSolve,"(/'   --- in staticUpdate')",advance='NO')
    end if
    if (present(noCorotUpd)) then
       write(dbgSolve,"(', noCorotUpd =',L2)",advance='NO') noCorotUpd
    end if
    if (present(doIterations)) then
       write(dbgSolve,"(', doIterations =',L2)",advance='NO') doIterations
    end if
    if (present(doFrictions)) then
       write(dbgSolve,"(', doFrictions =',L2)",advance='NO') doFrictions
    end if
    write(dbgSolve,"()")
    dbgTime = sys%time
    dbgIter = sys%nIterThisStep
#endif

    call startTimer (upd_p)

    if (sys%nStep > int(maxInt_p,i8)) then
       iStep = 0 ! Avoid overflow, used by user-defined elements only
    else
       iStep = int(sys%nStep)
    end if

    if (present(delta)) then
       if (.not. present(noCorotUpd)) then
          linInc = 0
       else if (noCorotUpd) then
          linInc = 1
       else
          linInc = 0
       end if

       !! Expand the equation-ordered delta to DOF-order
       if (sys%nIterThisStep == 0) then ! Predictor step
          call csExpand (sam,delta,sys%sinc(:,1))
       else ! Corrector step, replace the prescribed values, if any, by zero
          call csExpand (sam,delta,sys%sinc(:,1),S2=0.0_dp)
       end if
       call TransSysVecToGlobal (mech%triads,sys%sinc(:,1))

       sys%rinc = sys%rinc + sys%sinc(:,1) ! Total displacement increment

       !! Increment all position variables
       call IncTriadsPos     (mech%triads,sys%sinc(:,1))
       call IncSupElsGenDofs (mech%sups,sys%sinc(:,1))
       call IncJointsVar     (mech%joints,sys%sinc(:,1),mech%motions)

    else
       linInc = 1
       sys%rinc = 0.0_dp
    end if

    call updateSeaEnvironment (mech%env,mech%triads,mech%sups,sys%time,999,ierr)
    if (ierr < 0) goto 900

    !! Update all mechanism objects (for dampers: update length only)
    call clearTriadForces      (mech%triads)
    call updateJoints          (mech%joints,mech%motions,ierr)
    call updateSupElsStatic    (mech%sups,mech%supLoads,mech%env, &
         &                      sys%time,sys%nIterThisStep,linInc,ierr)
    call updateUDEs            (mech%elms,mech%env,sys%time,-sys%timeStep, &
         &                      iStep,sys%nIterThisStep,ierr)
    call updateDampers         (mech%dampers,sys%timeStep,.false.,ierr, &
         &                      updateVar=.false.)
    call updateSpringYields    (mech%springYields,ierr)
    call updateSprings         (mech%axialSprings,mech%joints,.false.,ierr)
    call updateBushingElements (mech%bElems,ierr,springsOnly=.true.)
    call updateContactElements (mech%cElems,sys%timeStep,ierr)
    if (present(doIterations)) then
       call updateTires (sys,mech%tires,ierr,staticIterations=doIterations)
    else
       call updateTires (sys,mech%tires,ierr,staticIterations=.true.)
    end if
    call updateExternalForces  (mech%forces,0,ierr)
    call updateTireForces      (mech%tires,mech%gravity)
    call updateMasses          (mech%masses,mech%gravity,.false.,ierr)
    if (present(doFrictions)) then
       !! Frictions are updated after all triad force contributions have been
       !! computed. Superelements, axial springs, external forces,
       !! additional masses and tires might have such contributions.
       call updateFrictions (mech%joints,mech%cElems,sys%timeStep, &
            &                sys%nIterThisStep,doFrictions,ierr)
    end if

900 if (ierr < 0) call reportError (debugFileOnly_p,'staticIncAndUpdate')

    call stopTimer (upd_p)

  end subroutine staticIncAndUpdate


  !!============================================================================
  !> @brief Updates the configuration after initial static equilibrium.
  !>
  !> @param[in] sam Data for managing system matrix assembly
  !> @param sys System level model data
  !> @param mech Mechanism components of the model
  !> @param[out] ierr Error flag
  !>
  !> @details This subroutine performs some final updates after the initial
  !> static equilibrium iterations have converged. Typically, this is (initial)
  !> velocity-dependent stuff, like dampers, etc. There is no need to update
  !> those during the iterations themselves which are purely static,
  !> only on the converged configuration.
  !> @note If the subsequent time history simulation is quasi-static, any
  !> non-zero initial conditions are nullified such that all velocity-dependent
  !> quantities will be initialized to zero instead.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 2 Jun 2004

  subroutine finalStaticUpdate (sam,sys,mech,ierr)

    use KindModule                  , only : dp, i8, maxInt_p
    use SamModule                   , only : SamType
    use SystemTypeModule            , only : SystemType, isQuasiStatic
    use MechanismTypeModule         , only : MechanismType
    use TriadTypeModule             , only : clearTriadForces, getTriadsVelAcc
    use MasterSlaveJointTypeModule  , only : updateAtConvergence,getJointsVelAcc
    use SupElTypeModule             , only : updateAtConvergence
    use SpringTypeModule            , only : checkFailure
    use FrictionTypeModule          , only : updateAtConvergence
    use MotionRoutinesModule        , only : InitPrescribedVelAcc
    use HydroDynamicsModule         , only : InitiateHydroDynBodies
    use SupElRoutinesModule         , only : updateSeaEnvironment
    use SupElRoutinesModule         , only : updateSupEls, SetSupElsVelAcc
    use UserdefElRoutinesModule     , only : updateUDEs
    use BushingElementRoutinesModule, only : updateBushingElements
    use SpringRoutinesModule        , only : updateSprings
    use DamperRoutinesModule        , only : updateDampers
    use ForceRoutinesModule         , only : updateExternalForces
    use MassRoutinesModule          , only : updateMasses
    use TireRoutinesModule          , only : updateTireForces
    use SolverRoutinesModule        , only : clearVelAcc
    use profilerModule              , only : startTimer, stopTimer, upd_p
    use reportErrorModule           , only : reportError, debugFileOnly_p
#ifdef FT_DEBUG
    use manipMatrixModule           , only : writeObject
    use dbgUnitsModule              , only : dbgSolve
    use FFaCmdLineArgInterface      , only : ffa_cmdlinearg_intValue
#endif

    type(SamType)      , intent(in)    :: sam
    type(SystemType)   , intent(inout) :: sys
    type(MechanismType), intent(inout) :: mech
    integer            , intent(out)   :: ierr

    !! Local variables
    integer :: i, iStep

    !! --- Logic section ---

#ifdef FT_DEBUG
    write(dbgSolve,"(/'   --- in finalStaticUpdate')")
#endif

    call startTimer (upd_p)

    if (sys%nStep > int(maxInt_p,i8)) then
       iStep = 0 ! To avoid overflow, only used by FNV so harmless here
    else
       iStep = int(sys%nStep)
    end if

    call updateSeaEnvironment (mech%env,mech%triads,mech%sups,sys%time,-1,ierr)
    if (ierr < 0) goto 900

    call clearTriadForces (mech%triads)

    if (isQuasiStatic(sys,sys%time+sys%timeStep)) then
       !! Quasi-static simulation, clear possible initial conditions
       call clearVelAcc (mech%triads,mech%joints,.true.)
    else ! Dynamics simulation
       !! Load initial velocities and/or accelerations into system vectors
       call getTriadsVelAcc (mech%triads,sys%urd,sys%urdd)
       call getJointsVelAcc (mech%joints,sys%urd,sys%urdd)
       !! Extract initial conditions from the prescribed motions
       call InitPrescribedVelAcc (mech%motions,sys%urd,sys%urdd)
       !! Pass the (initial) velocities and accelerations to the superelements
       call SetSupElsVelAcc (mech%sups,sam,sys%urd,sys%urdd)
    end if

    !! Update some velocity-dependent mechanism objects
    call updateSupEls          (mech%sups,mech%supLoads,mech%env, &
         &                      0.0_dp,0.0_dp,sys%time,0.0_dp, &
         &                      iStep,-sys%nIterThisStep,.true.,ierr)
    call updateUDEs            (mech%elms,mech%env,sys%time,sys%timeStep, &
         &                      iStep,-sys%nIterThisStep,ierr)
    call updateSprings         (mech%axialSprings,mech%joints,.false.,ierr, &
         &                      updateLength=.false.,updateVar=.true.)
    call updateDampers         (mech%dampers,0.0_dp,.false.,ierr, &
         &                      updateLV=.false.,updateVar=.true.)
    call updateBushingElements (mech%bElems,ierr,dampersOnly=.true.)
    call updateExternalForces  (mech%forces,0,ierr)
    call updateTireForces      (mech%tires,mech%gravity)
    call updateMasses          (mech%masses,mech%gravity,.false.,ierr)
    if (ierr < 0) goto 900

    !! If failure criterion has been exceeded, set hasFailed for next time step
    do i = 1, size(mech%baseSprings)
       call checkFailure (mech%baseSprings(i))
    end do

    !! Update all previous state variables defining the initial configuration

    sys%urdp  = sys%urd
    sys%urddp = sys%urdd

    do i = 1, size(mech%cElems)
       mech%cElems(i)%cVarPrev = mech%cElems(i)%cVar(1:3,1)
       mech%cElems(i)%wasActive = mech%cElems(i)%isActive
    end do

    do i = 1, size(mech%triads)
       if (mech%triads(i)%nDOFs > 0) then
          mech%triads(i)%urPrev = mech%triads(i)%ur
       end if
    end do

    do i = 1, size(mech%sups)
       call updateAtConvergence (mech%sups(i))
    end do

    call InitiateHydroDynBodies (mech%sups,mech%elms,mech%env,.false.,ierr)
    if (ierr < 0) goto 900

    do i = 1, size(mech%joints)
       call updateAtConvergence (mech%joints(i))
    end do

    do i = 1, size(mech%frictions)
       if (associated(mech%frictions(i)%p)) then
          call updateAtConvergence (mech%frictions(i)%p)
       end if
    end do

900 if (ierr < 0) call reportError (debugFileOnly_p,'finalStaticUpdate')

    call stopTimer (upd_p)

#ifdef FT_DEBUG
    if (ffa_cmdlinearg_intValue('debug') == -99) then
       call writeObject(sys%rinc,dbgSolve,'sys%rinc')
       call writeObject(sys%urd ,dbgSolve,'sys%urd')
       call writeObject(sys%urdd,dbgSolve,'sys%urdd')
    end if
#endif

  end subroutine finalStaticUpdate

end module staticEquilibriumModule
