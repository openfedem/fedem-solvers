!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file solverModule.f90
!>
!> @brief Model data and driver routines for the FEDEM Dynamics Solver.

!!==============================================================================
!> @brief Module with model containers and top lever driver subroutines.
!>
!> @details This module contains all model data and the top-level subroutines of
!> the FEDEM Dynamics solver. It replaces the old top-level driver subroutine of
!> version R7.2 and earlier. The new driver is divided into three main parts:
!>
!> - initialize() - Reads the input files and sets up everything.
!> - solvestep() - Advances the dynamic solution one step forward.
!> - finalize() - Closes down the result database after successful termination.
!>
!> In addition, this module contains several utility subroutines for external
!> manipulation of the solution process through the solver API.

module solverModule

  use kindModule         , only : dp, lfnam_p
  use SamModule          , only : SamType
  use ModesTypeModule    , only : ModesType
  use SystemTypeModule   , only : SystemType, SysMatrixType
  use ControlTypeModule  , only : ControlType
  use MechanismTypeModule, only : MechanismType

  implicit none

  private

  !! Model data containers
  type(SamType)      , save :: sam   !< Data for managing system matrix assembly
  type(ModesType)    , save :: modes !< Data for eigenmodes
  type(SystemType)   , save :: sys   !< System level model data
  type(ControlType)  , save :: ctrl  !< Control system data
  type(MechanismType), save :: mech  !< Mechanism objects of the model
  type(SysMatrixType), save :: Amat  !< System matrix for external communication

  !> @brief Externally set load vector.
  !> @details Assumed constant during Newton iterations.
  real(dp), allocatable, save :: extRhs(:)

  !! Internal global control parameters (undocumented)
  !> @cond NO_DOCUMENTATION
  logical , save :: lDouble(2), lEnergyInt, lSave, anyRes(5)
  integer , save :: IO, iPrint, iSimError, resFilFM, firstStep = 1
  integer , save :: numEigSol, numStiffDampSkip, NewmarkFlag, printFunc
  real(dp), save :: flushInc, lastFlush, bufRat(5)
  real(dp), save :: alpha(2), stopGlbDmp
  real(dp), save :: printInc, lastPrint, eigInc, lastEig
  real(dp), save :: rdbInc, lastRDB, dumpBeams
#ifdef FT_HAS_RECOVERY
  integer , save :: lRec
  real(dp), save :: recInc, lastRec
#endif
  !> @endcond

  character(len=lfnam_p), save :: chmodel    !< Model file name
  character(len=lfnam_p), save :: chnames(5) !< Result database file names
  character(len=lfnam_p), save :: yamlFile   !< Mode shape file name

  character(len=:), allocatable, save :: frsNames !< Recovery result files

  public :: sam, sys, mech, resFilFM ! For frequency response analysis module
  public :: initialize, softRestart, solveStep, solveRampUp, solveModes
  public :: finalize, closeAll
  public :: solveDynamic, getEngine, getEngineId, getTime, setNewTime
  public :: stateVectorSize, saveState
  public :: partStateVectorSize, savePartState
  public :: strainGagesSize, saveInitGageStrains
  public :: getSystemMatrix, getElementMatrix, getRhsVector, setRhsVector
  public :: systemSize, objectEquations, objectStateVar, haveResults
  public :: solverParameters, solveLinEqSystem
  public :: computeGageStrains, computeBeamForces, computeRelativeDistance
  public :: computeResponseVars, getJointSpringStiffness


contains

  !!============================================================================
  !> @brief Initializes the simulation model.
  !>
  !> @param[in] chfsi Text string containing the model definition
  !> @param[in] stateData Current state to restart from
  !> @param[in] ndat Length of the stateData array
  !> @param[in] xinp External function values for initial state
  !> @param nxinp Number of external function values to extract/extracted
  !> @param[out] ierr Error flag
  !>
  !> @details This subroutine reads and sets up a new simulation model,
  !> defines the initial configuration and initializes the result database.
  !> It also handles restart, both from file and from a solution state
  !> provided as argument.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 30 Nov 2016

  subroutine initialize (chfsi,stateData,ndat,xinp,nxinp,ierr)

    use SysMatrixTypeModule       , only : nullifySysMatrix
    use SystemTypeModule          , only : isQuasiStatic
    use SensorTypeModule          , only : numCtrlOut
    use FunctionTypeModule        , only : isCtrlSysUsed
    use FrictionTypeModule        , only : fricForceTol
    use InitiateModule            , only : readSolverData, preprocessSolverData
    use InitiateSystemTypeModule  , only : initTimeStepping
    use InitiateSupElTypeModule   , only : writeSupEls2Ftn
    use RestartModule             , only : restartInit
    use AddInSysModule            , only : getForceVectors
    use TimeStepModule            , only : getInitialTimeStepSize
    use SolverRoutinesModule      , only : openVTF, printMassDistribution
    use SolverRoutinesModule      , only : calculateIterationTolerances, meqErr
    use SolverRoutinesModule      , only : saveStep, updatePreviousState
    use EnergyIntegrationModule   , only : energyInitialization
    use StaticEquilibriumModule   , only : StaticEquilibrium
    use HydroDynamicsModule       , only : diffractionCalc
#ifdef FT_HAS_RECOVERY
    use StressRecoveryModule      , only : initRecovery, writeRecoveryHeaders
    use StressRecoveryModule      , only : writeRosettes2Ftn, gageRecovery
    use StressRecoveryModule      , only : initGageStrains, saveGageResults
    use ProfilerModule            , only : rec_p, rec1_p, sav_p
#endif
    use SaveModule                , only : initiateSaveModule
    use SaveModule                , only : writeSolverHeaders
    use FileUtilitiesModule       , only : getFileName
    use VersionModule             , only : openResFile
    use TimerModule               , only : initTime, startTimer, stopTimer
    use ProfilerModule            , only : ini_p, nProfMod_p
    use ProgressModule            , only : writeProgress, lterm
    use ReportErrorModule         , only : reportError, error_p, note_p
    use ReportErrorModule         , only : allocationError
    use ReportErrorModule         , only : openTerminalOutputFile
    use FFaCmdLineArgInterface    , only : ffa_cmdlinearg_getint
    use FFaCmdLineArgInterface    , only : ffa_cmdlinearg_getbool
    use FFaCmdLineArgInterface    , only : ffa_cmdlinearg_getdouble
    use FFaCmdLineArgInterface    , only : ffa_cmdlinearg_doubleValue
    use FFaCmdLineArgInterface    , only : ffa_cmdlinearg_isTrue
    use FFaCmdLineArgInterface    , only : ffa_cmdlinearg_isSet

    character(len=*), optional, intent(in)    :: chfsi
    real(dp)        , optional, intent(in)    :: stateData(*), xinp(*)
    integer         , optional, intent(in)    :: ndat
    integer         , optional, intent(inout) :: nxinp
    integer                   , intent(out)   :: ierr

    !! Local variables
    character(len=lfnam_p) :: chname
    logical                :: didInitEquil, doTimeStepping
    integer                :: nchar
    real(dp)               :: dt, mass

    !! --- Logic section ---

    call initTime (nProfMod_p)
    call startTimer (ini_p)

    iSimError = 0 ! Error flag for simulation error, but keep results!
    numEigSol = 0 ! Total number of eigenvalue solutions

    !! Open a file for terminal output, if requested.
    !! Default is to use standard out (Fortran unit 6).
    call openTerminalOutputFile (lterm)
    write(lterm,"(/11X,16('='),'> START OF PROGRAM SOLVER <',16('='))")

    !! Open the solver result file and write heading
    call getFileName ('resfile',chname)
    call openResFile (chname,'Dynamics Solver',IO,ierr)
    if (ierr /= 0) return


    !! --- Data input ----------------------------------------------------------

    call initiateSaveModule ()
    call nullifySysMatrix (Amat)

    !! Read model data from the solver input file
    call writeProgress (' --> READ SOLVER INPUT')
    call readSolverData (chfsi,sys,ctrl,mech,chmodel,lEnergyInt,ierr)
    if (ierr /= 0) goto 915

    if (ffa_cmdlinearg_isTrue('dumpSupEls')) then

       !! Write superelement matrices to Fortran source code files
       call writeSupEls2Ftn (mech%sups,ierr)
       if (ierr /= 0) goto 915

    end if


    !! --- Preprocessing of mechanism model and simulation parameters ----------

    !! Get output options from the command-line
    call ffa_cmdlinearg_getint ('debug',iprint)
    call ffa_cmdlinearg_getint ('printFunc',printFunc)
#ifdef FT_HAS_RECOVERY
    call ffa_cmdlinearg_getint ('recovery',lRec)
    if (lRec > 10) lRec = mod(lRec,10)
#endif
    call ffa_cmdlinearg_getint ('resFileFormat',resFilFM)
    if (resFilFM == 1) resFilFM = -1 ! Use the old res-file format

    call ffa_cmdlinearg_getbool ('double1',lDouble(1))
    call ffa_cmdlinearg_getbool ('double2',lDouble(2))
    call ffa_cmdlinearg_getdouble ('dumpBeamForces',dumpBeams)
    call ffa_cmdlinearg_getdouble ('eiginc',eigInc)
    if ( dumpBeams >= 0.0_dp .or. iprint > 0 .or. printFunc == 3 .or. &
         ffa_cmdlinearg_isSet('printTriad') ) then
       call ffa_cmdlinearg_getdouble ('printinc',printInc)
    else
       printInc = -1.0_dp ! Suppress additional res-file print
    end if

    call preprocessSolverData (sam,sys,ctrl,modes,mech,iprint,ierr)
    if (ierr /= 0) goto 915

    if (isQuasiStatic(sys,sys%tStart)) then
       !! The simulation will start with quasi-static load increments
       if (size(mech%tires) > 0) then
          ierr = -1
          call reportError (error_p,'Quasi-static analysis is not '// &
               &            'available for models containing tires')
          goto 915
       else if (isCtrlSysUsed(mech%engines)) then
          ierr = -1
          call reportError (error_p,'Quasi-static analysis is not '// &
               &            'available for models with control systems')
          goto 915
       else if (ctrl%saveVar) then
          ctrl%saveVar = .false. ! Switch off control output when quasi-static
          call reportError (note_p,'Control lines will not have any value '// &
               &            'in Quasi-static simulations')
       end if
    end if

    !! Initialize the convergence tolerances
    call calculateIterationTolerances (sam,sys)

    lastPrint = sys%tStart - printInc
    lastFlush = sys%tStart
    lastEig   = sys%tStart - eigInc ! Do an eigenvalue analysis also at t=tStart
    if (modes%nModes < 1) eigInc = -1.0_dp ! No eigenvalue analysis at all

    !! Compute result files buffer size ratios
    bufRat = 0.0_dp
    call ffa_cmdlinearg_getdouble ('flushinc',flushInc)
    if (flushInc > 0.0_dp) then
       if (associated(sys%tIncEngine) .or. sys%varInc > 0) then
          bufRat(1) = flushInc/sys%minInc
       else
          bufRat(1) = flushInc/sys%tInc
       end if
       call ffa_cmdlinearg_getdouble ('saveinc2',dt)
       if (dt >= sys%minInc) then
          bufRat(2) = flushInc/dt
       else if (dt >= 0.0_dp) then
          bufRat(2) = bufRat(1)
       end if
       if (eigInc >= sys%minInc) then
          bufRat(3) = max(1.0_dp,flushInc/eigInc)
       else if (eigInc >= 0.0_dp) then
          bufRat(3) = max(1.0_dp,bufRat(1))
       end if
       call ffa_cmdlinearg_getdouble ('saveinc4',dt)
       if (dt >= sys%minInc) then
          bufRat(4) = max(1.0_dp,flushInc/dt)
       else if (dt >= 0.0_dp) then
          bufRat(4) = max(1.0_dp,bufRat(2))
       end if
       if (ffa_cmdlinearg_isTrue('frequency_domain')) then
          call ffa_cmdlinearg_getdouble ('sample_freq',dt)
          bufRat(5) = max(1.0_dp,flushInc*dt)
       end if
    end if

    !! Print out some key model parameters
    write(io,100) sam%nel,sam%nnod,sam%mpar(27),sam%neq,sam%mpar(22)
100 format(//5X,'FE MODEL SUMMARY' / 5X,16('-') &
         & //5X,'Number of elements          =',I8 &
         &  /5X,'Number of nodes (total)     =',I8 &
         &  /5X,'Number of external nodes    =',I8 &
         &  /5X,'Number of equations         =',I8 &
         &  /5X,'Number of generalized modes =',I8 / )
    call printMassDistribution (mech%sups,mech%masses,mech%tires, &
         &                      mech%elms,mech%gravity,mass,IO)

    if (size(mech%frictions) > 0) then
       !! Friction force tolerance: <tolFricForce>*mass*g
       fricForceTol = ffa_cmdlinearg_doubleValue('tolFricForce')*mass*9.8_dp * &
            &         ffa_cmdlinearg_doubleValue('scaleToS') * &
            &         ffa_cmdlinearg_doubleValue('scaleToS') / &
            &         ffa_cmdlinearg_doubleValue('scaleToM')
       if (fricForceTol > 0.0_dp) write(io,1000) fricForceTol
1000   format(5X,'Zero force tolerance in friction calculation:',1PE12.5/)
       if (sys%minIt < 3) then
          sys%minIt = 3
          call reportError (note_p,'Setting minimum number of iterations'// &
               &            ' to 3 because this a model with friction')
       end if
    end if

    call stopTimer (ini_p)

    if (ffa_cmdlinearg_isTrue('dumpSupEls')) then

#ifdef FT_HAS_RECOVERY
       if (lRec == 2 .or. lRec == 3) then
          call initRecovery (mech%sups,lRec,rec1_p,iprint,IO,ierr)
          if (ierr < 0) goto 915

          !! Write strain recovery matrices to Fortran source code files
          call writeRosettes2Ftn (mech%sups,ierr)
          if (ierr /= 0) goto 915

       end if
#endif

       call closeAll (ierr)
       if (ierr /= 0) goto 915

       ierr = 1
       write(lterm,"(/11X,16('='),'>  END OF PROGRAM SOLVER  <',16('='))")
       close(IO)
       return

    end if

    if (ffa_cmdlinearg_isTrue('findAllSingularities')) then
       allocate(meqErr(sam%neq),stat=ierr)
    else
       allocate(meqErr(1),stat=ierr)
    end if
    if (ierr /= 0) then
       ierr = allocationError('solver::meqErr')
       goto 915
    end if


    !! --- Initialization of the result database -------------------------------

    call writeProgress (' --> WRITING RESULT DATABASE FILE HEADERS')
    call getFileName ('frs1file',chnames(1))
    call getFileName ('frs2file',chnames(2))
    call getFileName ('modesfile',chnames(3))
    call getFileName ('ctrlfile',chnames(4))
    if (ffa_cmdlinearg_isTrue('frequency_domain')) then
       call getFileName ('freqfile',chnames(5))
    else
       chnames(5) = ''
    end if
    call writeSolverHeaders (chnames,mech,ctrl,modes,bufRat,sam%neq, &
         &                   lDouble,ierr,modelFile=chmodel)
    if (ierr < 0) goto 915

    lSave = .true.
    lastRDB = sys%tStart
    call ffa_cmdlinearg_getdouble ('rdblength',rdbInc)

    call getFileName ('VTFfile',chname)
    if (chname /= '') then
       call openVTF (chname,ierr)
    end if

    !! Get file name prefix for export of system mode shapes to YAML files
    call getFileName ('yamlFile',yamlFile)


    !! --- Initialize for concurrent stress recovery, if requested

#ifdef FT_HAS_RECOVERY
    if (lRec /= 0) then
       call initRecovery (mech%sups,lRec,rec1_p,iprint,IO,ierr)
       if (ierr < 0) goto 915

       if (ffa_cmdlinearg_isSet('frs3file')) then
          !! Allocate a character string long enough to hold the frs-file names
          !! for all superelements that are targeted for stress recovery
          nchar = lfnam_p * numRecoveryParts(mech%sups)
          allocate(character(len=nchar) :: frsNames, stat=ierr)
          if (ierr == 0) then
             call getFileName ('frs3file',frsNames)
          else
             ierr = allocationError('solver::frsNames')
             goto 915
          end if

          call startTimer (sav_p)
          call writeRecoveryHeaders (frsNames,chmodel,mech%sups, &
               &                     bufRat(2),abs(lRec),ierr)
          call stopTimer (sav_p)
          if (ierr < 0) goto 915
       end if

       lastRec = sys%tStart
       call ffa_cmdlinearg_getdouble ('saveinc2',recInc)
    end if
#endif


    !! --- Establish the initial configuration ---------------------------------

    didInitEquil = .false.
    doTimeStepping = .true.
    if (present(stateData) .and. present(ndat)) then
       if (present(nxinp)) nxinp = 0

       !! Restart from the state provided
       call startTimer (ini_p)
       call restartInit (sam,sys,mech,ctrl,stateData,ndat,IO,ierr)
       call stopTimer (ini_p)
       if (ierr /= 0) goto 915

    else if (ffa_cmdlinearg_doubleValue('restarttime') >= 0.0) then
       if (present(nxinp)) nxinp = 0

       !! Restart from a given time step of a previous run
       call startTimer (ini_p)
       call restartInit (sam,sys,mech,ctrl,IO,ierr)
       call stopTimer (ini_p)
       if (ierr /= 0) goto 915

    else

       !! Set initial time and time increment size
       call initTimeStepping (sys,ierr)
       if (ierr < 0) goto 915

       !! Do quasi-static equilibrium iterations to find the start configuration
       call StaticEquilibrium (sam,sys,mech,ctrl,xinp,nxinp,IO,IERR)
       if (ierr /= 0) goto 915

       dt = 0.1_dp*sys%minInc
       doTimeStepping = sys%tEnd < sys%time-dt .or. sys%tEnd > sys%time+dt
       didInitEquil   = sys%nIterThisStep > 0

    end if

    call ffa_cmdlinearg_getint ('diffraction',nchar)
    if (nchar > 0 .and. doTimeStepping) then

       !! Compute hydrodynamic loads from diffraction
       call diffractionCalc (mech%env,mech%env%waveFunc,ierr)
       if (ierr /= 0) goto 915

    end if

    call ffa_cmdlinearg_getdouble ('stopGlbDmp',stopGlbDmp)
    if (stopGlbDmp < 0.0_dp .or. sys%time <= stopGlbDmp) then
       call ffa_cmdlinearg_getdouble ('alpha1',alpha(1))
       call ffa_cmdlinearg_getdouble ('alpha2',alpha(2))
    else
       alpha = 0.0_dp ! No global structural damping
    end if

    !! Compute the initial system force vectors
    call getForceVectors (sys%FSk,sys%FDk,sys%FIk,sys%Qk,sys%RFk, &
         &                sam,mech,sys%nStep,alpha(2),ierr)
    if (ierr /= 0) goto 915

    if (lEnergyInt .and. doTimeStepping) then

       call startTimer (ini_p)
       call writeProgress (' --> ENERGY INITIALIZATION')
       call ffa_cmdlinearg_getint ('num_damp_energy_skip',numStiffDampSkip)
       call energyInitialization (mech,sys%time,ierr)
       call stopTimer (ini_p)
       if (ierr /= 0) goto 915

    end if

    !! Initialize previous response variables based on the initial configuration
    !! in case they are used in inverse methods, etc.
    call updatePreviousState (sam,sys,mech,ierr)
    if (ierr < 0) goto 915

#ifdef FT_HAS_RECOVERY
    if (abs(lRec) == 2 .or. lRec == 3) then

       if (didInitEquil) then
          !! Do gage recovery on the initial equilibrium configuration
          call startTimer (rec_p)
          call gageRecovery (mech%sups,sys%nStep,sys%time,0,IO,ierr)
          call stopTimer (rec_p)
          if (ierr /= 0) goto 915
       end if

       call startTimer (sav_p)
       call saveGageResults (mech%sups,sys%nStep,sys%time,ierr)
       call stopTimer (sav_p)
       if (ierr < 0) goto 915

       if (.not. didInitEquil) then
          !! Reset the zeroInit flag when not doing initial equilibrium,
          !! such that the current (zero) strains are considered being
          !! the initial strains, without further calculations.
          call initGageStrains (0,bufRat,ierr)
       end if

    end if
#endif

    !! Save initial configuration to response database file
    call saveStep (sys,mech,ctrl,ierr)
    if (ierr < 0) goto 915


    !! --- Ready to start the time integration loop ----------------------------

    if (doTimeStepping) then
       call writeProgress (' --> TIME INTEGRATION')
       call ffa_cmdlinearg_getint ('NewmarkFlag',NewmarkFlag)
       write(lterm,200)
    end if
200 format(/10X,'Step | It.no. | Sim. Time | Time Inc. | ', &
         &      'Vel. norm | Res. norm |' /10X,63('-') )

    ierr = 0
    return

915 call closeAll (ierr,'initialize')

  contains

    !> @brief Counts the number of superelements with recovery enabled.
    function numRecoveryParts (sups) result(nrec)
      use SupElTypeModule, only : SupElType
      type(SupElType), intent(in) :: sups(:)
      integer :: i, nrec
      nrec = 0
      do i = 1, size(sups)
         if (associated(sups(i)%rcy)) nrec = nrec + 1
      end do
    end function numRecoveryParts

  end subroutine initialize


  !!============================================================================
  !> @brief Restarts a running simulation model from the provided state.
  !>
  !> @param[in] stateData Current state to restart from
  !> @param[in] ndat Length of the stateData array
  !> @param[in] writeToRDB Control variable for frs-file output
  !> @param[out] ierr Error flag
  !>
  !> @details The value of the argument @a writeToRDB is interpreted as follows:
  !> - &lt; 1 : Suppress all frs-output
  !> -   =  1 : Continue with output to same files as the previous run
  !> -   =  2 : Create a new set of files for this restart run
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 6 Feb 2017

  subroutine softRestart (stateData,ndat,writeToRDB,ierr)

    use RestartModule       , only : restartInit
#ifdef FT_HAS_RECOVERY
    use StressRecoveryModule, only : writeRecoveryHeaders
    use ProfilerModule      , only : sav_p
#endif
    use SaveModule          , only : closeResultFiles, writeSolverHeaders
    use ProfilerModule      , only : startTimer, stopTimer, ini_p
    use ProgressModule      , only : writeProgress, lterm
    use ReportErrorModule   , only : reportError, note_p

    real(dp), intent(in)  :: stateData(*)
    integer , intent(in)  :: ndat, writeToRDB
    integer , intent(out) :: ierr

    !! --- Logic section ---

    write(IO,"('')")
    lSave = writeToRDB > 0
    if (writeToRDB == 2) then
       call writeProgress (' --> INCREMENTING RESULT DATABASE FILES')
       call reportError (note_p,'Incrementing result database files')

       call closeResultFiles (anyRes,ierr)
       if (ierr < 0) goto 915

       lastRDB = sys%time
       lastFlush = sys%time
       call writeSolverHeaders (chnames,mech,ctrl,modes,bufRat,sam%neq, &
            &                   lDouble,ierr,modelFile=chmodel)
       if (ierr < 0) goto 915

#ifdef FT_HAS_RECOVERY
       if (lRec /= 0 .and. allocated(frsNames)) then
          call startTimer (sav_p)
          call writeRecoveryHeaders (frsNames,chmodel,mech%sups, &
               &                     bufRat(2),abs(lRec),ierr)
          call stopTimer (sav_p)
          if (ierr < 0) goto 915
       end if
#endif
    end if

    call startTimer (ini_p)
    call restartInit (sam,sys,mech,ctrl,stateData,ndat,IO,ierr)
    call stopTimer (ini_p)
    if (ierr /= 0) goto 915

    call writeProgress (' --> TIME INTEGRATION')
    write(lterm,200)
200 format(/10X,'Step | It.no. | Sim. Time | Time Inc. | ', &
         &      'Vel. norm | Res. norm |' /10X,63('-') )

    return

915 call closeAll (ierr,'softRestart')

  end subroutine softRestart


  !!============================================================================
  !> @brief Solves for the next time step or iteration.
  !>
  !> @param iop Control variable defining what to do (see below)
  !> @param[in] finalStep If .true., this is the final step of current
  !>                      time window
  !> @param[out] finished If .true., the end of the simulation has been reached,
  !>                      or the time integration failed
  !> @param[out] ierr Error flag
  !>
  !> @details This subroutine performs different sub-tasks of the solution
  !> process, depending on the value of the control variable iop, as follows:
  !> - =  0 : Normal operation, solve the entire step with iterations
  !> - = -1 : Establish the linear system of the first (prediction) iteration,
  !>          and then return to the calling module
  !> - = -2 : Same as -1, but also factorize the Newton matrix
  !> - = -3 : Same as 0, but this is the first time step after restart
  !> - =  1 : We are doing quasi-static iterations (internal mode)
  !> - =  2 : We are doing Newmark iterations (internal mode)
  !> - =  3 : Do one quasi-static iteration on the existing linear system
  !> - =  4 : Do one Newmark iteration on the existing linear system
  !> - =  5 : As 3, and then continue until equilibrium
  !> - =  6 : As 4, and then continue until equilibrium
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 30 Nov 2016

  subroutine solveStep (iop,finalStep,finished,ierr)

    use KindModule             , only : i8
    use TimeStepModule         , only : getTimeStepSize
    use SystemTypeModule       , only : isQuasiStatic
    use SolverRoutinesModule   , only : printResults
    use SolverRoutinesModule   , only : calculateIterationTolerances
    use SolverRoutinesModule   , only : terminateStep, saveStep, dumpStep
    use SolverRoutinesModule   , only : reportEquationError, meqErr
    use FiniteElementModule    , only : reportBeamEndForces, envBeamEndForces
    use EnergyIntegrationModule, only : energyIntegration
    use StaticEquilibriumModule, only : StaticInt
    use NewmarkRoutinesModule  , only : NewmarkInt
    use ModesRoutinesModule    , only : eigenModes, printModes, exportModes
    use EngineRoutinesModule   , only : updateEnginesForSave, printFuncValues
#ifdef FT_HAS_RECOVERY
    use StressRecoveryModule   , only : writeRecoveryHeaders, flushRecoveryFiles
    use StressRecoveryModule   , only : stressRecovery, gageRecovery
    use StressRecoveryModule   , only : saveGageResults
    use ProfilerModule         , only : rec_p, rec2_p
#endif
    use SaveModule             , only : writeSolverHeaders, modeDB
    use SaveModule             , only : flushResultFiles, closeResultFiles
    use ModesTypeModule        , only : writeModesDB
    use TimerModule            , only : startTimer, stopTimer
    use ProfilerModule         , only : sav_p
    use ProgressModule         , only : writeProgress
    use ReportErrorModule      , only : internalError
    use ReportErrorModule      , only : reportError, debugFileOnly_p
    use ReportErrorModule      , only : error_p, warning_p, note_p

    integer, intent(inout) :: iop
    logical, intent(in)    :: finalStep
    logical, intent(out)   :: finished
    integer, intent(out)   :: ierr

    !! Local variables
    real(dp)          :: t
    logical           :: hasReached, linearStatic
    character(len=64) :: errMsg

    !! Statement function
    hasReached(t) = sys%time + 0.1_dp*sys%minInc > t

    !! --- Logic section ---

    if (iSimError /= 0) then
       ierr = internalError('solveStep: Invoked after simulation failure')
       return
    else if (mod(iop,2) == 1) then
       goto 200
    else if (iop > 0) then
       goto 400
    end if

100 continue
    if (printInc >= 0.0_dp .and. hasReached(lastPrint+printInc)) then

       if (printInc > 0.0_dp) then
          lastPrint = sys%time ! Print at fixed time intervals
       else
          printInc = -1.0_dp ! Print at t=tStart only
       end if
       if (printFunc == 3 .and. sys%nStep < 1000000_i8) then
          call printFuncValues (int(sys%nStep),sys%time,mech%engines,IO)
       end if
       if (dumpBeams >= 0.0_dp) then
          call reportBeamEndForces (sys%time,mech%sups,mech%gravity,IO)
       end if
       call printResults (sys%time,mech%triads,mech%sups,iprint,IO)

       if (resFilFM == 1) resFilFM = 2 ! Must write new convergence heading
    end if
    if (dumpBeams >= 0.0_dp .and. sys%time >= dumpBeams) then
       call envBeamEndForces (mech%sups,ierr)
       if (ierr < 0) goto 915
    end if
    if (eigInc >= 0.0_dp .and. hasReached(lastEig+eigInc)) then

       lastEig = sys%time
       numEigSol = numEigSol + 1
       call writeProgress (' --> EIGENVALUE CALCULATIONS')

       call eigenModes (sam,modes,sys,mech,iprint,ierr)
       if (ierr < -1) then
          if (ierr < -10) then
             meqErr(1) = -ierr/10
             call reportEquationError(sam,mech%sups,mech%triads,mech%joints,-5)
          end if
          write(errMsg,"('Eigenvalue analysis failed at T=',1PE12.5)") sys%time
          call reportError (warning_p,errMsg, &
               'No modal results were stored for this time step.', &
               'However, the time integration continues.')
       else if (ierr /= 0) then
          call reportError (error_p,'Eigenvalue analysis failed.', &
               'You can try to rerun the dynamic simulation', &
               'with the eigenvalue analysis switched off.')
          goto 915
       else
          call printModes (sam,modes,sys%time,5*iprint,IO)
          call exportModes (yamlFile,chmodel,modes,mech%triads,mech%sups, &
               &            sys%nStep,sys%time,ierr)
          if (ierr < 0) goto 915
          call startTimer (sav_p)
          call writeModesDB (modes,sys,mech%triads,mech%sups, &
               &             modeDB,lDouble(1),ierr)
          call stopTimer (sav_p)
          if (ierr < 0) goto 915
       end if

       if (resFilFM == 1) resFilFM = 2 ! Must write new convergence heading
    end if

    !! Exit the time step loop ?
    ierr = 0
    finished = sys%tStart <= sys%tEnd .and. hasReached(sys%tEnd)
    if (finished) return

    if (rdbInc > 0.0_dp .and. hasReached(lastRDB+rdbInc)) then

       write(errMsg,"(', at T =',1PE12.5)") sys%time
       call writeProgress (' --> INCREMENTING RESULT DATABASE FILES')
       write(IO,"('')")
       call reportError (note_p,'Incrementing result database files'//errMsg)
       call closeResultFiles (anyRes,ierr)
       if (ierr < 0) goto 915

       lastRDB = sys%time
       lastFlush = sys%time
       call writeSolverHeaders (chnames,mech,ctrl,modes,bufRat,sam%neq, &
            &                   lDouble,ierr,modelFile=chmodel)
       if (ierr < 0) goto 915

#ifdef FT_HAS_RECOVERY
       if (lRec /= 0 .and. allocated(frsNames)) then
          call startTimer (sav_p)
          call writeRecoveryHeaders (frsNames,chmodel,mech%sups, &
               &                     bufRat(2),abs(lRec),ierr)
          call stopTimer (sav_p)
          if (ierr < 0) goto 915
       end if
#endif

    else if (flushInc >= 0.0_dp .and. hasReached(lastFlush+flushInc)) then

       lastFlush = sys%time
       call flushResultFiles (ierr)
       if (ierr < 0) goto 915

#ifdef FT_HAS_RECOVERY
       if (lRec /= 0) then
          call flushRecoveryFiles (sav_p,ierr)
          if (ierr < 0) goto 915
       end if
#endif

    end if

    !! Set the (possibly) velocity proportional tolerances
    call calculateIterationTolerances (sam,sys)

    !! Get new time step size
    sys%timePrev = sys%timeStep
    if (sys%time >= sys%tStart) then
       sys%timeStep = getTimeStepSize(sys,sam,ctrl,1.0_dp,ierr)
       if (ierr < 0) goto 915
    end if

    !! Check if next time step should be solved dynamically
    if (.not. isQuasiStatic(sys,sys%time+sys%timeStep)) goto 300

    !! Advance the quasi-static solution one step
    if (firstStep == 1) then
       firstStep = 2
       write(IO,"(/'  Starting quasi-static analysis')")
    end if
    if (iop == -3) iop = 0
200 if (allocated(extRhs)) then
       call StaticInt (sam,sys,mech,ctrl,extRhs,iop,resFilFM,IO,iSimError)
    else
       call StaticInt (sam,sys,mech,ctrl,iop=iop,resFileFormat=resFilFM, &
            &          lpu=IO,ierr=iSimError)
    end if
    linearStatic = sys%fixedIt < 0
    goto 500

    !! Advance dynamics solution one step through Newmark time integration
300 if (firstStep == 1) then
       write(IO,"(/'  Starting Newmark time integration')")
    else if (firstStep == 2) then
       write(IO,"(/'  Continuing with Newmark time integration')")
       sys%urd = sys%rinc / sys%timePrev ! Estimate starting velocity
    end if
    if (firstStep > 0) then
       firstStep = 0
       write(IO,301) sys%alpha,sys%beta,sys%gamma
       if (NewmarkFlag > 199) then
          write(IO,302) 1.0_dp-sys%alpha_m,1.0_dp-sys%alpha_f
       end if
301    format(5X,'alpha =',F8.4,'  beta =',F8.5,'  gamma =',F8.5)
302    format(5X,'alpha_m =',F9.5,'  alpha_f =',F8.5)
    end if
    if (stopGlbDmp >= 0.0_dp .and. sys%time+sys%timeStep > stopGlbDmp) then
       alpha = 0.0_dp ! Switch off the global structural damping
    end if

400 if (allocated(extRhs)) then
       call NewmarkInt (sam,sys,mech,ctrl,alpha,extRhs, &
            &           iop,NewmarkFlag,resFilFM,IO,iSimError)
    else
       call NewmarkInt (sam,sys,mech,ctrl,alpha, &
            &           iop=iop,NewmarkFlag=NewmarkFlag, &
            &           resFileFormat=resFilFM,lpu=IO,ierr=iSimError)
    end if
    linearStatic = .false.

500 ierr = abs(iSimError)
    finished = iSimError /= 0
    if (finished) goto 916 ! Time integration failure, terminate
    if (iop > 0) return ! Exit for external iteration processing

    sys%nStep = sys%nStep + 1_i8
    if (allocated(extRhs)) then
       extRhs = 0.0_dp
    end if

    if (lEnergyInt) then

       call energyIntegration (mech,sys%time,sys%timeStep, linearStatic, &
            &                  int(numStiffDampSkip,i8) < sys%nStep, ierr)
       if (ierr < 0) goto 915

       mech%Etot = mech%Epot + mech%Ekin + mech%Estr - mech%Einp &
            &    + mech%Edmp - mech%Eext

    end if

    finished = sys%tStart <= sys%tEnd .and. hasReached(sys%tEnd)
    call terminateStep (sam,sys,mech,ctrl, &
         &              linearStatic, finalStep.or.finished, ierr)
    if (ierr < 0) goto 915

#ifdef FT_HAS_RECOVERY
    if (lRec == 1 .or. lRec == 3) then

       if (recInc < 0.0_dp .or. hasReached(lastRec+recInc)) then
          lastRec = sys%time
          call stressRecovery (mech%sups,sys%nStep,lastRec,lSave, &
               &               rec_p,rec2_p,sav_p,resFilFM,iprint,IO,ierr)
          if (ierr /= 0) goto 915
       end if

    end if
    if (abs(lRec) == 2 .or. lRec == 3) then

       call startTimer (rec_p)
       call gageRecovery (mech%sups,sys%nStep,sys%time, &
            &             max(0,iprint)*resFilFM,IO,ierr)
       call stopTimer (rec_p)
       if (ierr /= 0) goto 915

       if (lSave) then
          call startTimer (sav_p)
          call saveGageResults (mech%sups,sys%nStep,sys%time,ierr)
          call stopTimer (sav_p)
          if (ierr < 0) goto 915
       end if

    end if
#endif

    if (lSave) then
       !! Save response variables to result database
       call saveStep (sys,mech,ctrl,ierr,.true.)
       if (ierr < 0) goto 915
    end if

    !! Update engines that are used for inter-process communication
    call updateEnginesForSave (mech%engines,1,ierr)
    if (ierr < 0) goto 915
    if (finished) goto 100 ! Check for eigenvalue calculation at the final time

    return

915 finished = .false.
    iSimError = ierr
    call dumpStep (sys,mech,ctrl,sam,fName='datastructure.inf')
    call closeAll (ierr,'solveStep')
    return

916 call reportError (debugFileOnly_p,'solveStep')

  end subroutine solveStep


  !!============================================================================
  !> @brief Solves the ramp-up stage, if any.
  !>
  !> @param iop Operation flag telling what to do, should be zero
  !> @param[out] finished If .true., the time integration failed
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 19 Aug 2022

  subroutine solveRampUp (iop,finished,ierr)

    use ReportErrorModule, only : reportError, debugFileOnly_p

    integer, intent(inout) :: iop
    logical, intent(out)   :: finished
    integer, intent(out)   :: ierr

    !! --- Logic section ---

    finished = .false.
    do while (.not.finished .and. sys%time < sys%tStart)
       call solveStep (iop,.false.,finished,ierr)
       if (ierr < 0) then
          call reportError (debugFileOnly_p,'solveRampUp')
          return
       end if
    end do

  end subroutine solveRampUp


  !!============================================================================
  !> @brief Solves the eigenvalue system at current configuration.
  !>
  !> @param[in] nModes Number of eigenmodes to solve for
  !> @param[out] eVal Computed eigenvalues (frequencies in [Hz])
  !> @param[out] eVec Mass-normalized eigenvectors
  !> @param[in] dofOrder If .true., expand eigenvectors to DOF-order
  !> @param[in] useLaPack Flag for using LAPack eigenvalue solvers
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 26 Oct 2020

  subroutine solveModes (nModes,eVal,eVec,dofOrder,useLaPack,ierr)

    use kindModule             , only : pi_p
    use initiateModesTypeModule, only : initiateModes
    use ModesRoutinesModule    , only : eigenModes, printModes, exportModes
    use SolverRoutinesModule   , only : reportEquationError, meqErr
    use ProgressModule         , only : writeProgress
    use ReportErrorModule      , only : internalError
    use ReportErrorModule      , only : reportError, debugFileOnly_p, warning_p
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_getint

    integer , intent(in)  :: nModes, useLaPack
    real(dp), intent(out) :: eVal(*), eVec(*)
    logical , intent(in)  :: dofOrder
    integer , intent(out) :: ierr

    !! Local variables
    integer           :: i, j, n
    character(len=64) :: errMsg
    real(dp), pointer :: eigVec(:,:)

    !! --- Logic section ---

    if (iSimError /= 0) then
       ierr = internalError('solveModes: Invoked after simulation failure')
       return
    end if

    if (useLaPack == 1) then
       !! Use the symmetric DSYGVX eigensolver
       call initiateModes (sam,sys,modes,ierr,nModes,5)
    else if (useLaPack == 2) then
       !! Use the non-symmetric DGGEVX eigensolver
       call initiateModes (sam,sys,modes,ierr,nModes,3)
    else
       !! Use the default eigensolver
       call initiateModes (sam,sys,modes,ierr,nModes)
    end if
    if (ierr < 0) goto 915

    numEigSol = numEigSol + 1
    call writeProgress (' --> EIGENVALUE CALCULATIONS')
    call eigenModes (sam,modes,sys,mech,iprint,ierr)
    if (ierr < -1) then
       if (ierr < -10) then
          meqErr(1) = -ierr/10
          call reportEquationError (sam,mech%sups,mech%triads,mech%joints,-5)
       end if
       write(errMsg,"('Eigenvalue analysis failed at T=',1PE12.5)") sys%time
       call reportError (warning_p,errMsg, &
            &            'No modal results were stored for this time step.')
       ierr = -ierr
       return
    else if (ierr < 0) then
       goto 915
    end if

    call printModes (sam,modes,sys%time,5*iprint,IO)
    call exportModes (yamlFile,chmodel,modes,mech%triads,mech%sups, &
         &            sys%nStep,sys%time,ierr)
    if (ierr < 0) goto 915

    if (dofOrder) then
       eigVec => modes%ReVec ! Eigenvectors in DOF-order
    else if (associated(modes%eqVec)) then
       eigVec => modes%eqVec ! Eigenvectors in equation order
    else
       ierr = internalError('solveModes: No equation-ordered eigenvectors')
       return
    end if

    j = 1
    n = size(eigVec,1)
    do i = 1, nModes
       eVal(i) = modes%ReVal(i)*0.5_dp/pi_p
       call DCOPY (n,eigVec(1,i),1,eVec(j),1)
       j = j + n
    end do

    if (resFilFM == 1) resFilFM = 2 ! Must write new convergence heading

    !! Reset from command-line argument
    call ffa_cmdlinearg_getint ('numEigModes',modes%nModes)

    return

915 call reportError (debugFileOnly_p,'solveModes')

  end subroutine solveModes


  !!============================================================================
  !> @brief Terminates the simulation and close the result database.
  !>
  !> @param[out] ierr Error flag
  !>
  !> @details This subroutine is invoked after the time integration has
  !> completed successfully, and there are no more calculations to be performed.
  !> Its main task is to close down the result database, but it also performs
  !> the automatic curve export when that has been requested.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 30 Nov 2016

  subroutine finalize (ierr)

    use KindModule                , only : i8
    use SystemTypeModule          , only : isQuasiStatic
    use SolverRoutinesModule      , only : dumpStep
    use SaveModule                , only : res1DB, res2DB, modeDB, freqDB
    use FiniteElementModule       , only : reportBeamEndForces
#ifdef FT_HAS_RECOVERY
    use StressRecoveryModule      , only : getGageRecoveryFiles
#endif
    use FileUtilitiesModule       , only : getFileName
    use ProfilerModule            , only : startTimer, stopTimer, exp_p
    use ProgressModule            , only : writeProgress, lterm
    use ReportErrorModule         , only : reportError, debugFileOnly_p, error_p
    use FFpBatchExportInterface   , only : ffp_crvexp

    integer, intent(out) :: ierr

    !! Local variables
    integer                :: i
    character(len=64)      :: errMsg
    character(len=lfnam_p) :: chname, crvFile

    !! --- Logic section ---

    if (sys%nIter > 1_i8) then
       write(lterm,"(10X,73('-'))")
       call writeProgress(' --> TIME INTEGRATION DONE. CLOSING RESULT DATABASE')
    else
       call writeProgress(' --> CLOSING RESULT DATABASE')
    end if

#ifdef FT_HAS_RECOVERY
    if (abs(lRec) == 2 .or. lRec == 3) then
       call getGageRecoveryFiles (chnames)
    else
       chnames = ''
    end if
#endif

    call closeAll (ierr)
    if (ierr /= 0) goto 915

    if (iSimError /= 0) then
       !! Simulation failure (divergent model, etc.)
       if (isQuasiStatic(sys)) then
          write(errMsg,1) sys%time
       else
          write(errMsg,2) sys%time
       end if
1      format('Quasi-static Simulation aborted at T =',1PE12.5)
2      format('Dynamic Simulation aborted at T =',1PE12.5)
       call reportError (error_p,errMsg, &
            &            'Results are saved until the last accepted time step.')
    else if (dumpBeams >= 0.0_dp .and. sys%time > dumpBeams) then
       !! Print out beam end force envelopes
       call reportBeamEndForces (-sys%time,mech%sups,mech%gravity,IO)
    end if

    if (iprint > 0) then
       call dumpStep (sys,mech,ctrl,sam,fName='datastructure.inf')
    end if

    if (sys%nUpdates < 100000000_i8 .and. sys%nIter < 1000000000000_i8) then
       write(IO   ,3) sys%nUpdates, sys%nIter
       write(lterm,3) sys%nUpdates, sys%nIter
3      format(//5X,'NUMBER OF MATRIX UPDATES:',I8, &
            &  /5X,'NUMBER OF ITERATIONS:',I12)
    end if


    !! --- Batch curve export --------------------------------------------------

    call getFileName ('curvePlotFile',crvFile)
    if (crvFile == '') goto 999

    if (freqDB%fileName /= '') then
       chname = freqDB%fileName
    else if (res1DB%fileName /= '') then
       chname = res1DB%fileName
       if (res2DB%fileName /= '') then
          chname = trim(chname) //','// res2DB%fileName
       end if
       if (modeDB%fileName /= '') then
          chname = trim(chname) //','// modeDB%fileName
       end if
    else if (res2DB%fileName /= '') then
       chname = res2DB%fileName
       if (modeDB%fileName /= '') then
          chname = trim(chname) //','// modeDB%fileName
       end if
    else if (modeDB%fileName /= '') then
       chname = modeDB%fileName
    else
       chname = ''
    end if
#ifdef FT_HAS_RECOVERY
    if (abs(lRec) == 2 .or. lRec == 3) then
       do i = 1, size(chnames)
          if (chnames(i) /= '') then
             if (chname == '') then
                chname = chnames(i)
             else if (len_trim(chname)+len_trim(chnames(i))+3 <= lfnam_p) then
                chname = trim(chname) //','// chnames(i)
             end if
          end if
       end do
    end if
#endif

    chname = '<'// trim(chname) //'>'
    call writeProgress (' --> EXPORTING RESULTS TO CURVE FILE(S)')
    call startTimer (exp_p)
    call ffp_crvexp (trim(chname),trim(crvFile),trim(chmodel),ierr)
    call stopTimer (exp_p)
    if (ierr == 0) goto 999

    call reportError (error_p,'Auto-export of curve data failed', &
         &            addString='ffp_crvexp')

915 call reportError (debugFileOnly_p,'finalize')
    call dumpStep (sys,mech,ctrl,sam,fName='datastructure.inf')
999 call writeClosure (lterm,ierr)
    call deallocateAll ()

  end subroutine finalize


  !!============================================================================
  !> @brief Closes the result database files and any external modules used.
  !>
  !> @param ierr Error flag
  !> @param[in] abortsub If present, gives the name of the aborting subroutine
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 30 Nov 2016

  subroutine closeAll (ierr,abortsub)

    use SaveModule               , only : closeResultFiles
    use SolverRoutinesModule     , only : closeSolverDB
    use HydroDynamicsModule      , only : closeHydroDyn
    use WindTurbineRoutinesModule, only : closeAeroDyn
    use dbgUnitsModule           , only : closeDbgUnits
#ifdef FT_HAS_RECOVERY
    use StressRecoveryModule     , only : closeRecovery
#endif
    use ProgressModule           , only : lterm
    use reportErrorModule        , only : reportError, debugFileOnly_p

    integer                   , intent(inout) :: ierr
    character(len=*), optional, intent(in)    :: abortsub

    !! Local variables
    integer :: lerr

    !! --- Logic section ---

    if (present(abortsub)) then
       call closeResultFiles (anyRes,lerr)
    else
       call closeSolverDB (mech,anyRes,ierr)
    end if

    call closeHydroDyn (IO)
    if (associated(mech%turbine)) then
       call closeAeroDyn (lerr)
       if (.not.present(abortsub) .and. ierr == 0) ierr = lerr
    end if

#ifdef FT_HAS_RECOVERY
    call closeRecovery (.true.,lerr)
    if (.not.present(abortsub) .and. ierr == 0) ierr = lerr
#endif

    call closeDbgUnits ()

    if (present(abortsub)) then
       call reportError (debugFileOnly_p,abortsub)
       call writeClosure (lterm,ierr)
       call deallocateAll ()
    else if (ierr < 0) then
       call reportError (debugFileOnly_p,'closeAll')
    end if

  end subroutine closeAll


  !!============================================================================
  !> @brief Outputs file and profiling information on program termination.
  !>
  !> @param[in] lcon File unit number for console messages
  !> @param ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 30 Nov 2016

  subroutine writeClosure (lcon,ierr)

    use FileUtilitiesModule, only : getFileName
    use ProfilerModule     , only : reportTiming
    use TimerModule        , only : showTime

    integer, intent(in)    :: lcon
    integer, intent(inout) :: ierr

    !! Local variables
    character(len=lfnam_p) :: chname

    !! --- Logic section ---

    write(lcon,600) ' END OF PROGRAM SOLVER '

    if (ierr == 0) then
       write(lcon,610)
       call writeFileName ('fsifile','SOLVER INPUT FILE')
       call writeFileName ('resfile','RESULT LOG FROM SOLVER')
       if (anyRes(1)) then
          call writeFileName ('frs1file','PRIMARY RESPONSE DATA FROM SOLVER')
       end if
       if (anyRes(2)) then
          call writeFileName ('frs2file','SECONDARY RESPONSE DATA FROM SOLVER')
       end if
       if (anyRes(3)) then
          call writeFileName ('modesfile','EIGENVALUE RESULTS FROM SOLVER')
       end if
       if (anyRes(4)) then
          call writeFileName ('ctrlfile','CONTROL SYSTEM DATA FROM SOLVER')
       end if
       if (anyRes(5)) then
          call writeFileName ('freqfile','FREQUENCY DOMAIN RESULTS FROM SOLVER')
       end if
       call writeFileName ('curvePlotFile','EXPORTED CURVE DATA')
       call writeFileName ('VTFfile','EXPORTED RIGID BODY ANIMATION')
       write(lcon,611)
       if (iSimError /= 0) then
          call getFileName ('resfile',chname)
          if (len_trim(chname) < 43) then
             chname = '             <'//trim(chname)//'>'
          else if (len_trim(chname) < 56) then
             chname = '<'//trim(chname)//'>'
          else
             chname = '<...'//trim(chname(len_trim(chname)-51:))//'>'
          end if
          chname(58:58) = '!'
          write(lcon,650) chname
       end if
    else
       call getFileName ('resfile',chname)
       if (len_trim(chname) < 43) then
          chname = '             <'//trim(chname)//'>'
       else if (len_trim(chname) < 56) then
          chname = '<'//trim(chname)//'>'
       else
          chname = '<...'//trim(chname(len_trim(chname)-51:))//'>'
       end if
       chname(58:58) = '!'
       write(lcon,690) chname
    end if

    call reportTiming (IO,sys%nIter,numEigSol,sys%nStep)
    call showTime (IO)

    if (ierr /= 0 .or. iSimError /= 0) then
       if (ierr == 0) ierr = iSimError
       write(IO,620) 'failed :-('
    else
       write(IO,620) 'successfully completed :-)'
    end if

    close(IO)

600 format(/11X,16('='),'> ',A,' <',16('='))
610 format(/11X,'-----------------------------------------------------------' &
         & /11X,'!                    fedem_solver                         !' &
         & /11X,'!---------------------------------------------------------!' &
         & /11X,'!             F I L E   I N F O R M A T I O N             !' &
         & /11X,'!---------------------------------------------------------!' )
611 format( 11X,'-----------------------------------------------------------'/)
620 format(/ 4X,'Simulation ',A)
650 format(/11X,'-----------------------------------------------------------' &
         & /11x,'!                                                         !' &
         & /11x,'!                  AN ERROR IS DETECTED                   !' &
         & /11x,'!                                                         !' &
         & /11x,'!             ERROR MESSAGES ARE WRITTEN TO               !' &
         & /11X,'!                                                         !' &
         & /11X,'!',A58, &
         & /11X,'!                                                         !' &
         & /11x,'!     RESULTS ARE KEPT UNTIL THE LAST ACCEPTED TIME STEP  !' &
         & /11X,'!                                                         !' &
         & /11X,'-----------------------------------------------------------'/)
690 format(/11X,'-----------------------------------------------------------' &
         & /11X,'!                    fedem_solver                         !' &
         & /11X,'!---------------------------------------------------------!' &
         & /11X,'!            E R R O R   I N F O R M A T I O N            !' &
         & /11X,'!---------------------------------------------------------!' &
         & /11X,'!                                                         !' &
         & /11X,'!                 AN ERROR IS DETECTED                    !' &
         & /11X,'!                                                         !' &
         & /11X,'!             ERROR MESSAGES ARE WRITTEN TO               !' &
         & /11X,'!                                                         !' &
         & /11X,'!',A58, &
         & /11X,'!                                                         !' &
         & /11X,'-----------------------------------------------------------'/)

  contains

    !> @brief Writes a file name with description to the console file unit
    subroutine writeFileName (chFile,chInfo)
      character(len=*), intent(in) :: chFile !< File name
      character(len=*), intent(in) :: chInfo !< File description
      integer :: i, idata
      call getFileName (chFile,chName)
      idata = len_trim(chName)
      if (idata > 22) then
         write(lcon,601) chName(idata-18:),chInfo,(' ',i=len(chInfo),30),'!'
      else if (idata > 0) then
         write(lcon,602) chName,chInfo,(' ',i=len(chInfo),30),'!'
      end if
601   format(11X,'! ...',A19,' : ',A,30A1)
602   format(11X,'! ',A22,' : ',A,30A1)
    end subroutine writeFileName

  end subroutine writeClosure


  !!============================================================================
  !> @brief Deallocates all dynamically allocated data.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Jan 2017

  subroutine deallocateAll ()

    use SamModule              , only : deallocateSAM
    use SystemTypeModule       , only : deallocateSys
    use SysMatrixTypeModule    , only : deallocateSysMatrix
    use ModesTypeModule        , only : deallocateModes
    use ControlTypeModule      , only : deallocateCtrl
    use MechanismTypeModule    , only : deallocateMechanism
    use TimeStepModule         , only : deallocateTimeStep
    use ExplicitFunctionsModule, only : deallocateSplines
    use FiniteElementModule    , only : deallocateBeams
    use AsmExtensionModule     , only : castToInt8
    use ModesRoutinesModule    , only : allocEigenMatrices
    use SolverRoutinesModule   , only : meqErr

    !! --- Logic section ---

    if (allocated(extRhs)) deallocate(extRhs)
    if (allocated(meqErr)) deallocate(meqErr)
    if (allocated(frsNames)) deallocate(frsNames)

    call deallocateTimeStep ()
    call deallocateSplines ()
    call deallocateBeams ()
    nullify(sam%meqn) ! avoid deallocation twice
    call deallocateSys (sys)
    call deallocateSAM (sam)
    call castToInt8 (sam)
    call allocEigenMatrices (modes)
    call deallocateModes (modes)
    call deallocateCtrl (ctrl)
    call deallocateMechanism (mech)
    call deallocateSysMatrix (Amat)

    firstStep = 1 ! In case solver is started again by the same process

  end subroutine deallocateAll


  !!==========================================================================
  !> @brief Returns whether the next step should be solved dynamically or not.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Guenter Glanzer
  !>
  !> @date 10 Jun 2018

  logical function solveDynamic ()

    use SystemTypeModule, only : isQuasiStatic

    !! --- Logic section ---

    solveDynamic = .not. isQuasiStatic(sys)

  end function solveDynamic


  !!==========================================================================
  !> @brief Returns whether current step have results to be saved.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 6 Dec 2024

  integer function haveResults ()

    !! --- Logic section ---

    if (lRec /= 1 .and. lRec /= 3) then
       haveResults = 1 ! Only transformation results (every time step)
    else if (sys%time - lastRec < 1.0e-12_dp) then
       haveResults = 2 ! This step have stress recovery results
    else
       haveResults = 0 ! No results for this time step
    end if

  end function haveResults


  !!============================================================================
  !> @brief Returns the simulation time of the current or next time step.
  !>
  !> @param[out] time The simulation time returned
  !> @param[in] nextStep If .true., return for next step, otherwise current step
  !> @param ierr Error flag
  !>
  !> @details This subroutine does not always work with @a nextStep = .true.
  !> if cut-back is performed, as the time step size then might be decreased.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 5 Dec 2016

  subroutine getTime (time,nextStep,ierr)

    use TimeStepModule, only : getTimeStepSize

    real(dp), intent(out)   :: time
    logical , intent(in)    :: nextStep
    integer , intent(inout) :: ierr

    !! --- Logic section ---

    if (nextStep) then
       time = sys%time + getTimeStepSize(sys,sam,ctrl,1.0_dp,ierr)
    else
       time = sys%time
    end if

  end subroutine getTime


  !!============================================================================
  !> @brief Sets the time increment size to be used for next time step.
  !>
  !> @param[in] nextTime Physical time to calculate time increment size from
  !> @param[out] ierr Error flag
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 22 Jun 2022

  subroutine setNewTime (nextTime,ierr)

    use reportErrorModule, only : reportError, error_p, warning_p, note_p

    real(dp), intent(in)  :: nextTime
    integer , intent(out) :: ierr

    !! Local variables
    real(dp)           :: deltaT
    character(len=128) :: msg

    !! --- Logic section ---

    deltaT = nextTime - sys%time
    if (deltaT < 0.0_dp) then
       ierr = -1
       write(msg,600) sys%time, nextTime
       call reportError (error_p,'Can not go back in time!',msg, &
            &            addString='setNewTime')
    else if (deltaT < sys%minInc) then
       ierr = 1
       write(msg,600) sys%time, nextTime
       call reportError (warning_p,'Time increment too small!',msg, &
            &            'Check your input, using minimum step size instead.')
       sys%tInc = sys%minInc
    else if (deltaT > sys%maxInc) then
       ierr = -2
       write(msg,600) sys%time, nextTime
       call reportError (error_p,'Time increment too large!',msg, &
            &            addString='setNewTime')
    else if (abs(deltaT-sys%tInc) > 0.1_dp*sys%tInc) then
       ierr = 2
       write(msg,600) sys%time, nextTime
       call reportError (note_p,'Time step size changed by more than 10%',msg)
       sys%tInc = deltaT
    else ! All good
       ierr = 0
       sys%tInc = deltaT
    end if

600 format('Current time =',1PE13.5,'  New time =',E13.5)

  end subroutine setNewTime


  !!============================================================================
  !> @brief Returns the current value of the specified engine.
  !>
  !> @param[out] value The function value
  !> @param ierr Error flag
  !> @param[in] userId User ID of the engine to evaluate
  !> @param[in] x Optional function argument value
  !> @param[in] tag Optional tag to identify the engine with
  !>
  !> @details If an engine with @a userId is found, the value of that engine
  !> is returned. Otherwise, the (first) engine with a matching @a tag is used.
  !> If no engine with the given identification is found, zero is returned.
  !> The error flag is decrementend on failure, othwerwise it is not touched.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 5 Dec 2016

  subroutine getEngine (value,ierr,userId,tag,x)

    use EngineRoutinesModule, only : EngineValue
    use reportErrorModule   , only : reportError
    use reportErrorModule   , only : debugFileOnly_p, warningFileOnly_p

    real(dp)        , intent(out)          :: value
    integer         , intent(inout)        :: ierr
    integer         , intent(in)           :: userId
    character(len=*), intent(in), optional :: tag
    real(dp)        , intent(in), optional :: x

    !! Local variables
    integer :: i, jerr
    logical :: found

    !! --- Logic section ---

    if (.not. associated(mech%engines)) then
       !! Invoked before model initialization, or after deallocation
       write(*,*) '*** getEngine invoked outside model scope'
       ierr = -999
    end if

    found = .false.

    jerr = ierr
    do i = 1, size(mech%engines)
       if (mech%engines(i)%id%userId == userId) then
          found = .true.
          value = EngineValue(mech%engines(i),ierr,x)
          exit
       else if (present(tag) .and. .not.found) then
          if (mech%engines(i)%tag == tag) then
             found = .true.
             value = EngineValue(mech%engines(i),ierr,x)
             if (userId <= 0) exit
          end if
       end if
    end do

    if (ierr < jerr) then
       call reportError (debugFileOnly_p,'getEngine',ierr=ierr)
    else if (found) then
       return
    else if (present(tag)) then
       call reportError (warningFileOnly_p,'No Function with tag = '//tag)
    end if

    value = 0.0_dp ! The specified engine does not exist, silently return zero

  end subroutine getEngine


  !!============================================================================
  !> @brief Returns the user Id of the (first) engine with the specified tag.
  !>
  !> @param[in] tag Tag to identify the engine with
  !> @return User ID of tagged function, or -1 if not found
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 25 Jun 2022

  function getEngineId (tag) result(userId)

    use reportErrorModule, only : reportError, error_p

    character(len=*), intent(in) :: tag

    !! Local variables
    integer :: i, userId

    !! --- Logic section ---

    if (.not. associated(mech%engines)) then
       !! Invoked before model initialization, or after deallocation
       write(*,*) '*** getEngineId invoked outside model scope'
       userId = -999
    end if

    do i = 1, size(mech%engines)
       if (mech%engines(i)%tag == tag) then
          userId = mech%engines(i)%id%userId
          return
       end if
    end do

    userId = -1 ! Tag not found
    call reportError (error_p,'No Function with tag = '//tag)

  end function getEngineId


  !!==========================================================================
  !> @brief Returns a system matrix as a full matrix.
  !>
  !> @param[out] Nmat The system matrix
  !> @param[in] iopM Option telling which matrix to return (see below)
  !> @param[out] ierr Error flag
  !>
  !> @details The value of iopM is iterpreted as follows:
  !> - = 1 : Return the system stiffness matrix
  !> - = 2 : Return the system mass matrix
  !> - = 3 : Return the system damping matrix
  !> - &gt; 3 : Return the current system Newton matrix
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 30 Mar 2017

  subroutine getSystemMatrix (Nmat,iopM,ierr)

    use AddInSysModule     , only : BuildStiffMat, BuildDamperMat, BuildMassMat
    use SysMatrixTypeModule, only : shareMatrixStructure, allocateSysMatrix
    use SysMatrixTypeModule, only : convertSysMat
    use reportErrorModule  , only : reportError, debugFileOnly_p

    real(dp), intent(out) :: Nmat(sam%neq,sam%neq)
    integer , intent(in)  :: iopM
    integer , intent(out) :: ierr

    !! --- Logic section ---

    ierr = 0
    if (iopM > 0 .and. iopM <= 3) then
       if (Amat%dim == 0) then
          call shareMatrixStructure (Amat,sys%Nmat)
          call allocateSysMatrix (Amat,ierr)
       end if
       select case (iopM+10*ierr)
       case (1)
          call BuildStiffMat (Amat,mech,sam,0,sys%stressStiffIsOn(2),0,ierr)
       case (2)
          call BuildMassMat (Amat,mech,sam,ierr)
       case (3)
          call BuildDamperMat (Amat,mech,sam,ierr)
       end select
       if (ierr == 0) then
          call convertSysMat (Amat,Nmat,sam%neq,1,1,ierr)
       end if
    else
       !! Return the current Newton matrix
       call convertSysMat (sys%Nmat,Nmat,sam%neq,1,1,ierr)
    end if
    if (ierr < 0) call reportError (debugFileOnly_p,'getNewtonMatrix')

  end subroutine getSystemMatrix


  !!============================================================================
  !> @brief Returns an element matrix for the specified superelement.
  !>
  !> @param[out] Emat The element matrix
  !> @param[in] baseId Base ID of the superelement to return the matrix for
  !> @param[in] iopM Option telling which matrix to return (see below)
  !> @param[out] ierr Error flag
  !>
  !> @details The value of iopM is iterpreted as follows:
  !> - = 1 : Return the element stiffness matrix
  !> - = 2 : Return the element mass matrix
  !> - = 3 : Return the element damping matrix
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 16 Feb 2018

  subroutine getElementMatrix (Emat,baseId,iopM,ierr)

    use SupElTypeModule  , only : SupElType, getPtrToId
    use reportErrorModule, only : reportError, error_p

    real(dp), intent(out) :: Emat(*)
    integer , intent(in)  :: baseId, iopM
    integer , intent(out) :: ierr

    !! Local variables
    type(SupElType), pointer :: supEl

    !! --- Logic section ---

    supEl => getPtrToId(mech%sups,baseId)
    if (.not. associated(supEl)) then
       ierr = baseId
       call reportError (error_p,'No superelement with this baseId',ierr=ierr, &
            &            addString='getElementMatrix')
       return
    end if

    ierr = 0
    select case (iopM)
    case (1); call DCOPY (size(supEl%Kmmat),supEl%Kmmat(1,1),1,Emat(1),1)
    case (2); call DCOPY (size(supEl%Mmat),supEl%Mmat(1,1),1,Emat(1),1)
    case (3); call DCOPY (size(supEl%Cmat),supEl%Cmat(1,1),1,Emat(1),1)
    end select

  end subroutine getElementMatrix


  !!============================================================================
  !> @brief Returns the current right-hand-side vector of the linearized system.
  !>
  !> @param[out] Rvec The right-hand-side vector
  !> @param[in] iopV If equal to 1, return current external load vector instead
  !> @param[out] ierr Error flag
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 30 Mar 2017

  subroutine getRhsVector (Rvec,iopV,ierr)

    real(dp), intent(out) :: Rvec(sam%neq)
    integer , intent(in)  :: iopV
    integer , intent(out) :: ierr

    !! --- Logic section ---

    if (iopV == 1) then
       call DCOPY (sam%neq,sys%Qk(1),1,Rvec(1),1)
    else
       call DCOPY (sam%neq,sys%del(1),1,Rvec(1),1)
    end if
    ierr = 0

  end subroutine getRhsVector


  !!============================================================================
  !> @brief Sets/Updates current right-hand-side vector of linearized system.
  !>
  !> @param Rvec The right-hand-side vector
  !> @param[in] addTo If .true., the provided vector is added to the current
  !> right-hand-side vector in the model, otherwise it replaces the current one
  !> @param[in] keepDuringIterations If .true., the provided vector should
  !> be applied in all iterations of current time step
  !> @param[out] ierr Error flag
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 30 Mar 2017

  subroutine setRhsVector (Rvec,addTo,keepDuringIterations,ierr)

    use reportErrorModule, only : allocationError

    real(dp), intent(in)  :: Rvec(sam%neq)
    logical , intent(in)  :: addTo, keepDuringIterations
    integer , intent(out) :: ierr

    !! --- Logic section ---

    if (keepDuringIterations) then
       if (.not. allocated(extRhs)) then
          allocate(extRhs(size(sys%del)),stat=ierr)
          if (ierr /= 0) then
             ierr = allocationError('setRhsVector')
             return
          end if
       end if
       if (sys%nIterThisStep == 0) then
          call DCOPY (sam%neq,Rvec(1),1,extRhs(1),1)
       end if
    else if (allocated(extRhs)) then
       deallocate(extRhs)
    end if

    if (addTo) then
       call DAXPY (sam%neq,1.0_dp,Rvec(1),1,sys%del(1),1)
    else
       call DCOPY (sam%neq,Rvec(1),1,sys%del(1),1)
    end if

    ierr = 0

  end subroutine setRhsVector


  !!============================================================================
  !> @brief Returns the dimension of the linearized system.
  !>
  !> @param[out] ndim Number of unknowns in the linear equation system
  !> @param[in] expanded If .true., return total number of DOFs in the system
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 30 Mar 2017

  subroutine systemSize (ndim,expanded)

    integer, intent(out)          :: ndim
    logical, intent(in), optional :: expanded

    !! --- Logic section ---

    if (.not. associated(sam%mpar)) then
       !! Invoked before model initialization, or after deallocation
       write(*,*) ' *** systemSize invoked outside model scope'
       ndim = -999
    else if (present(expanded)) then
       if (expanded) then
          ndim = sam%ndof
       else
          ndim = sam%neq
       end if
    else
       ndim = sam%neq
    end if

  end subroutine systemSize


  !!==========================================================================
  !> @brief Returns parameters for the Newmark/HHT-algorithm.
  !>
  !> @param[out] dt Time step size
  !> @param[out] alpha Numerical damping parameter of the HHT-scheme
  !> @param[out] beta Newmark time integration parameter
  !> @param[out] gamma Newmark time integration parameter
  !>
  !> @callergraph
  !>
  !> @author Guenter Glanzer
  !>
  !> @date 10 Jun 2018

  subroutine solverParameters (dt,alpha,beta,gamma)

    real(dp), intent(out) :: dt, alpha, beta, gamma

    !! --- Logic section ---

    dt    = sys%timeStep
    alpha = sys%alpha
    beta  = sys%beta
    gamma = sys%gamma

  end subroutine solverParameters


  !!==========================================================================
  !> @brief Returns the equation numbers for the specified object.
  !>
  !> @param[in] baseId Base ID of the object to get equation numbers for
  !> @param[out] meqn Array of equation numbers
  !> @param[out] ndof Number of degrees of freedom for the object in question
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 30 Mar 2017

  subroutine objectEquations (baseId,meqn,ndof)

    integer, intent(in)  :: baseId
    integer, intent(out) :: meqn(*), ndof

    !! Local variables
    integer :: i, j

    !! --- Logic section ---

    ndof = 0
    do i = 1, size(mech%triads)
       if (mech%triads(i)%id%baseId == baseId) then
          if (associated(mech%triads(i)%eqNumber)) then
             ndof = size(mech%triads(i)%eqNumber)
             meqn(1:ndof) = mech%triads(i)%eqNumber
          end if
          return
       end if
    end do
    do i = 1, size(mech%joints)
       if (mech%joints(i)%id%baseId == baseId) then
          ndof = mech%joints(i)%nJointDOFs
          do j = 1, ndof
             meqn(j) = sam%meqn(mech%joints(i)%jointDofs(j)%sysDOF)
          end do
          return
       end if
    end do

  end subroutine objectEquations


  !!==========================================================================
  !> @brief Returns the current state variables for the specified object.
  !>
  !> @param[in] baseId Base ID of the object to get equation numbers for
  !> @param[out] vars Array of state variables (positions only for Triads)
  !> @param[out] nVar Number of state variables
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 20 Mar 2020

  subroutine objectStateVar (baseId,vars,nVar)

    integer , intent(in)  :: baseId
    real(dp), intent(out) :: vars(*)
    integer , intent(out) :: nVar

    !! Local variables
    integer :: i, j

    !! --- Logic section ---

    nVar = 0
    do i = 1, size(mech%triads)
       if (mech%triads(i)%id%baseId == baseId) then
          if (mech%triads(i)%nDOFs >= 3) then
             nVar = 3
             vars(1:3) = mech%triads(i)%ur(:,4)
          end if
          return
       end if
    end do
    do i = 1, size(mech%joints)
       if (mech%joints(i)%id%baseId == baseId) then
          nVar = mech%joints(i)%nJointDOFs
          do j = 1, nVar
             vars(j) = mech%joints(i)%jointDofs(j)%jVar(1)
          end do
          return
       end if
    end do

  end subroutine objectStateVar


  !!============================================================================
  !> @brief Solves current linear equation system for a set of right-hand-sides.
  !>
  !> @param rhs The right-hand-side vectors of the linear equation system
  !> @param[in] nrhs Number of right-hand-side vectors
  !> @param[out] ierr Error flag
  !>
  !> @details Assuming the Newton matrix already have been factorized.
  !> To be used by the inverse problem algorithm.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 16 Jan 2017

  subroutine solveLinEqSystem (rhs,nrhs,ierr)

    use profilerModule    , only : startTimer, stopTimer, sol_p
    use solExtensionModule, only : csSolve

    real(dp), intent(inout) :: rhs(:)
    integer , intent(in)    :: nrhs
    integer , intent(out)   :: ierr

    !! Local variables
    integer, parameter :: solveMode_p = 4

    !! --- Logic section ---

    call startTimer (sol_p)
    call csSolve (solveMode_p,0,sys%Nmat,rhs,nrhs,IO,ierr)
    call stopTimer (sol_p)

  end subroutine solveLinEqSystem


  !!==========================================================================
  !> @brief Returns the dimension of the state vector.
  !>
  !> @param[in] transOnly If .true., consider position matrices only
  !> @param[out] ndat Required length of the state array
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 16 Jan 2017

  subroutine stateVectorSize (transOnly,ndat)

    use RestartModule, only : getStateSize, getTransformationStateSize

    logical, intent(in)  :: transOnly
    integer, intent(out) :: ndat

    !! --- Logic section ---

    if (.not. associated(mech%triads)) then
       !! Invoked before model initialization, or after deallocation
       write(*,*) ' *** stateVectorSize invoked outside model scope'
       ndat = -999
    else if (transOnly) then
       ndat = getTransformationStateSize(mech)
    else
       ndat = getStateSize(sam,mech,ctrl)
    end if

  end subroutine stateVectorSize


  !!============================================================================
  !> @brief Saves current state to an in-core array.
  !>
  !> @param[in] transOnly If .true., consider position matrices only
  !> @param[out] stateData Array containing all solution dependent variables
  !> @param[in] ndat Length of the stateData array
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 16 Jan 2017

  subroutine saveState (transOnly,stateData,ndat,ierr)

    use RestartModule, only : saveToCore, saveTransformationDataToCore

    logical , intent(in)  :: transOnly
    real(dp), intent(out) :: stateData(*)
    integer , intent(in)  :: ndat
    integer , intent(out) :: ierr

    !! --- Logic section ---

    if (transOnly) then
       call saveTransformationDataToCore (sys,mech,stateData,ndat,ierr)
    else
       call saveToCore (sam,sys,mech,ctrl,stateData,ndat,ierr)
    end if

  end subroutine saveState


  !!============================================================================
  !> @brief Returns the state vector dimension for a FE part.
  !>
  !> @param[in] iopS Option telling for which state, 1: deformation, 2: stresses
  !> @param[in] bid Base ID if the FE part to consider
  !> @param[in] ndat Required length of the state array
  !>
  !> @details The deformation state size equals 4 + 3*(number of nodes) whereas
  !> the stress state size euals 4 + number of nodes. In both cases the first 4
  !> entries are the step number, time, step size and baseID, respectiviely.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Runar Heggelien Refsnaes
  !>
  !> @date 29 Oct 2017

  subroutine partStateVectorSize (iopS,bid,ndat)

#ifdef FT_HAS_RECOVERY
    use StressRecoveryModule, only : getDeformationSize, getStressSize
#endif

    integer, intent(in)  :: iopS, bid
    integer, intent(out) :: ndat

    !! --- Logic section ---

#ifdef FT_HAS_RECOVERY
    if (iopS == 1) then
       ndat = getDeformationSize(bid)
    else
       ndat = getStressSize(bid)
    end if
    if (ndat > 0) ndat = ndat + 4
#else
    ndat = 0
    print *,' ** partStateVectorSize dummy:',bid,iopS
#endif

  end subroutine partStateVectorSize


  !!============================================================================
  !> @brief Saves the deformation/stress state of a FE part to an in-core array.
  !>
  !> @param[in] iopS Option telling what to save, 1: deformation, 2: stress
  !> @param[in] bid Base ID if the FE part to consider
  !> @param[out] data Deformation/stress state array
  !> @param[in] ndat Length of the data array
  !> @param[out] ierr Error flag
  !>
  !> @details If @a iopS equals 1, the von Mises stress value at each node of
  !> the part is saved, otherwise the 3 deformational displacement components
  !> at each node are saved. As the first four values, the current step number,
  !> the current time, the time step size, and the base ID of the FE part,
  !> respectively, are saved.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Runar Heggelien Refsnaes
  !>
  !> @date 31 Oct 2017

  subroutine savePartState (iopS,bid,data,ndat,ierr)

#ifdef FT_HAS_RECOVERY
    use StressRecoveryModule, only : getStress, getDeformation
    use reportErrorModule   , only : reportError, debugFileOnly_p
#endif

    integer , intent(in)  :: iopS, bid
    real(dp), intent(out) :: data(*)
    integer , intent(in)  :: ndat
    integer , intent(out) :: ierr

    !! --- Logic section ---

    data(1) = sys%nStep
    data(2) = sys%time
    data(3) = sys%timeStep
    data(4) = real(bid,dp)

    ierr = 5
#ifdef FT_HAS_RECOVERY
    if (iopS == 1) then
       call getDeformation (bid,data,ndat,ierr)
    else
       call getStress (bid,data,ndat,ierr)
    end if
    if (ierr < 0) then
       call reportError (debugFileOnly_p,'savePartState')
    end if
#else
    print *,' ** savePartState dummy:',bid,iopS,ndat
#endif

  end subroutine savePartState


  !!============================================================================
  !> @brief Returns the dimension of the strain gages vector.
  !>
  !> @param[out] ndat Required length of the strain gage state array
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 18 Dec 2017

  subroutine strainGagesSize (ndat)

#ifdef FT_HAS_RECOVERY
    use StressRecoveryModule, only : getStrainGagesSize
#endif

    integer, intent(out) :: ndat

    !! --- Logic section ---

#ifdef FT_HAS_RECOVERY
    ndat = getStrainGagesSize()
#else
    ndat = 0
#endif

  end subroutine strainGagesSize


  !!============================================================================
  !> @brief Save/restores the initial gage strain to/from an in-core array.
  !>
  !> @param[in] iopS If equal to 1 save to data array, otherwise restore from it
  !> @param data Array with strain gage values
  !> @param[in] ndat Length of the data array
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 18 Dec 2017

  subroutine saveInitGageStrains (iopS,data,ndat,ierr)

#ifdef FT_HAS_RECOVERY
    use StressRecoveryModule, only : initGageStrains
#endif

    integer , intent(in)    :: iopS
    real(dp), intent(inout) :: data(*)
    integer , intent(in)    :: ndat
    integer , intent(out)   :: ierr

    !! --- Logic section ---

#ifdef FT_HAS_RECOVERY
    call initGageStrains (iopS,data(1:ndat),ierr)
#else
    ierr = 0
    print *,' ** saveInitGageStrains dummy:',iopS,ndat,data(1)
#endif

  end subroutine saveInitGageStrains


  !!============================================================================
  !> @brief Computes strain tensor at gage positions for input displacements.
  !>
  !> @param[in] disp    Displacement vector
  !> @param[in] gageIds Base IDs of strain gages to calculate for
  !> @param[in] nDisp   Length of displacement vector (if 0, use internal state)
  !> @param[in] nGage   Length of gageIds array (number of strain gages)
  !> @param[out] eps    Strain tensors at gage positions
  !> @param[out] ierr   Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Guenter Glanzer
  !>
  !> @date 9 Sep 2020

  subroutine computeGageStrains (disp,gageIds,nDisp,nGage,eps,ierr)

    use SolverRoutinesModule, only : restoreLastStep
#ifdef FT_HAS_RECOVERY
    use StressRecoveryModule, only : gageRecovery
    use ReportErrorModule   , only : reportError, warning_p, debugFileOnly_p
#endif

    real(dp), intent(in)  :: disp(*)
    integer , intent(in)  :: gageIds(*), nDisp, nGage
    real(dp), intent(out) :: eps(*)
    integer , intent(out) :: ierr

    !! --- Logic section ---

    if (nDisp > 0) then
       call applySolutionVec (disp)
    end if

#ifdef FT_HAS_RECOVERY
    !! The sups(i)%finit array has to be correctly initialised
    if (nDisp > 0) then
       call gageRecovery (mech%sups,sys%nStep,sys%time,0,IO,ierr)
    else
       !! Use current state
       ierr = 0
    end if
    if (ierr /= 0) then
       call reportError (debugFileOnly_p,'computeGageStrains')
    else
       !! Extract the gage strains
       call getGageStrains (gageIds,eps)
    end if
#else
    ierr = 0
    eps(1:3*nGage) = 0.0_dp
    print *,' ** computeGageStrains dummy:',gageIds(1:nGage)
#endif

    !! Reset all position variables from the stored previous values
    if (nDisp > 0) then
       call restoreLastStep (sam,sys,mech)
    end if

  contains

    !> @brief Updates the configuration by applying given displacement vector.
    subroutine applySolutionVec (disp)

      use TriadTypeModule    , only : TransSysVecToGlobal, IncTriadsPos
      use SupElRoutinesModule, only : IncSupElsGenDofs
      use SupElTypeModule    , only : BuildFinit
      use SolExtensionModule , only : csExpand

      real(dp), intent(in) :: disp(nDisp)

      integer :: i

      !! Expand the equation-ordered solution vector to DOF-order
      call csExpand (sam,disp,sys%sinc(:,1))

      !! Update the sys%sinc(:,1) vector
      call TransSysVecToGlobal (mech%triads,sys%sinc(:,1))

      !! Update the triad(i)%ur(1:3,1:4) matrices
      call IncTriadsPos (mech%triads,sys%sinc(:,1))

      !! Update the sups(i)%genDOFs%ur(:) vectors
      call IncSupElsGenDofs (mech%sups,sys%sinc(:,1))

      !! Calculate the increment in sups(i)%finit for the adapted configuration
      do i = 1, size(mech%sups)
         call BuildFinit (mech%sups(i))
         mech%sups(i)%finit = mech%sups(i)%finit - mech%sups(i)%finitPrev
      end do

    end subroutine applySolutionVec

#ifdef FT_HAS_RECOVERY
    !> @brief Extracts gage strains from the strain rosette objects.
    subroutine getGageStrains (gageIds,epsil)

      use StrainRosetteModule , only : StrainRosetteType
      use StressRecoveryModule, only : getStrainRosette

      integer , intent(in)  :: gageIds(nGage)
      real(dp), intent(out) :: epsil(3,nGage)

      integer                          :: i
      type(StrainRosetteType), pointer :: pGage

      do i = 1, nGage
         pGage => getStrainRosette(gageIds(i))
         if (associated(pGage)) then
            epsil(:,i) = pGage%epsC(1:3)
         else
            epsil(:,i) = 0.0_dp
            call reportError (warning_p,'No strain rosette with this base ID', &
                 &            ierr=gageIds(i))
         end if
      end do

    end subroutine getGageStrains
#endif

  end subroutine computeGageStrains


  !!============================================================================
  !> @brief Computes beam sectional forces for input displacements.
  !>
  !> @param[in] disp    Displacement vector
  !> @param[in] beamIds Base IDs of the beams/triads to calculate for
  !> @param[in] nDofs   Length of displacement vector
  !> @param[in] nBeams  Length of beam identifier array
  !> @param[out] forces Beam sectional forces
  !> @param[out] ierr   Error flag
  !>
  !> @details The array @a beamIds contains triplets of (beamId, triadId, dof),
  !> where triadId is the base ID of one of the two end triads of the beam
  !> element, and dof is a zero-based force component index (-1 implies all).
  !> This subroutine always returns all force components, i.e.,
  !> the dof index is not used here.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Guenter Glanzer
  !>
  !> @date 20 Sep 2020

  subroutine computeBeamForces (disp,beamIds,nDofs,nBeams,forces,ierr)

    use SolverRoutinesModule, only : restoreLastStep
    use reportErrorModule   , only : reportError, warning_p, debugFileOnly_p

    real(dp), intent(in)  :: disp(*)
    integer , intent(in)  :: beamIds(*), nDofs, nBeams
    real(dp), intent(out) :: forces(*)
    integer , intent(out) :: ierr

    !! --- Logic section ---

    !! Calculate internal forces from displacement increment
    call forcesFromDisp (disp,ierr)
    if (ierr < 0) return

    !! Extract the beam section forces
    call getBeamSectionForces (beamIds,forces)

    !! Reset all position variables from the stored previous values
    call restoreLastStep (sam,sys,mech)

  contains

    !> @brief Updates the configuration by applying given displacement vector.
    subroutine forcesFromDisp (disp,ierr)

      use TriadTypeModule    , only : TransSysVecToGlobal, IncTriadsPos
      use TriadTypeModule    , only : clearTriadForces
      use SupElRoutinesModule, only : IncSupElsGenDofs, updateSupElsStatic
      use FiniteElementModule, only : updateBeamSectionForces
      use SolExtensionModule , only : csExpand

      real(dp), intent(in)    :: disp(nDofs)
      integer , intent(inout) :: ierr

      !! Expand the equation-ordered solution vector to DOF-order
      call csExpand (sam,disp,sys%sinc(:,1))

      !! Update the sys%sinc(:,1) vector
      call TransSysVecToGlobal (mech%triads,sys%sinc(:,1))

      !! Update the triad(i)%ur(1:3,1:4) matrices
      call IncTriadsPos (mech%triads,sys%sinc(:,1))

      !! Update the sups(i)%genDOFs%ur(:) vectors
      call IncSupElsGenDofs (mech%sups,sys%sinc(:,1))

      !! All triad forces are set to zero
      call clearTriadForces (mech%triads)

      !! Superelement update (static case).
      !! linInc=2 indicates no corotational frame update,
      !! and the increment finit-finitPrev will be calculated.
      call updateSupElsStatic (mech%sups,mech%supLoads,mech%env, &
           &                   sys%time,sys%nIterThisStep, &
           &                   linInc=2,ierr=ierr)

      !! Beam sectional force update (based on superelement-level quantities)
      call updateBeamSectionForces (mech%sups,ierr)

      if (ierr < 0) call reportError (debugFileOnly_p,'forcesFromDisp')

    end subroutine forcesFromDisp

    !> @brief Extracts sectional forces at triads from the beam elements.
    subroutine getBeamSectionForces (beamIds,forces)

      use SupElTypeModule, only : SupElType, getPtrToId, IsBeam
      use TriadTypeModule, only : TriadType

      integer , intent(in)  :: beamIds(3,nBeams)
      real(dp), intent(out) :: forces(6,nBeams)

      integer                  :: i, j
      type(SupElType), pointer :: pBeam
      type(TriadType), pointer :: pTriad

      forces = 0.0_dp
      do i = 1, nBeams
         pBeam => getPtrToId(mech%sups,beamIds(1,i))
         if (.not. associated(pBeam)) then
            call reportError (warning_p,'No superelement with this base ID', &
                 &            ierr=beamIds(1,i))
         else if (.not. IsBeam(pBeam)) then
            call reportError (warning_p,'This superelement is not a beam', &
                 &            ierr=beamIds(1,i))
         else
            do j = 1, 2
               pTriad => pBeam%triads(j)%p
               if (associated(pTriad%supForce) .and. pTriad%isFEnode) then
                  if (pTriad%id%baseId == beamIds(2,i)) then
                     forces(:,i) = pTriad%supForce
                  end if
               end if
            end do
         end if
      end do

    end subroutine getBeamSectionForces

  end subroutine computeBeamForces


  !!============================================================================
  !> @brief Computes relative distance between 2 triads for input displacements.
  !>
  !> @param[in]  disp   Displacement vector
  !> @param[in]  Ids    User IDs of the functions with relative sensors
  !> @param[out] relDis Relative distances
  !> @param[in]  nDofs  Length of displacement vector
  !> @param[in]  nIds   Length of function identifier array
  !> @param[out] ierr   Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Guenter Glanzer
  !>
  !> @date 22 Oct 2020

  subroutine computeRelativeDistance (disp,Ids,relDis,nDofs,nIds,ierr)

    use TriadTypeModule   , only : TransSysVecToGlobal, IncTriadsPos
    use SolExtensionModule, only : csExpand
    use reportErrorModule , only : reportError, debugFileOnly_p

    real(dp), intent(in)  :: disp(*)
    integer , intent(in)  :: Ids(*), nDofs, nIds
    real(dp), intent(out) :: relDis(*)
    integer , intent(out) :: ierr

    !! Local variables
    integer  :: i
    real(dp) :: relDist

    !! --- Logic section ---

    !! Extract current relative distance between 2 triads
    ierr = 0
    do i = 1, nIds
       call getEngine (relDis(i),ierr,Ids(i))
    end do
    if (ierr < 0) goto 915

    !! Expand the equation-ordered solution vector to DOF-order
    call csExpand (sam,disp(1:nDofs),sys%sinc(:,1))

    !! Update the sys%sinc(:,1) vector
    call TransSysVecToGlobal (mech%triads,sys%sinc(:,1))

    !! Update the triad(i)%ur(1:3,1:4) matrices
    call IncTriadsPos (mech%triads,sys%sinc(:,1))

    !! Extract relative distance change due to the given displacement state
    do i = 1, nIds
       call getEngine (relDist,ierr,Ids(i))
       relDis(i) = relDist - relDis(i)
    end do
    if (ierr < 0) goto 915

    !! Reset triad positions from the stored previous values
    do i = 1, size(mech%triads)
       if (mech%triads(i)%nDOFs > 0) then
          mech%triads(i)%ur = mech%triads(i)%urPrev
       end if
    end do

    return

915 continue
    call reportError (debugFileOnly_p,'computeRelativeDistance',ierr=ierr)

  end subroutine computeRelativeDistance


  !!============================================================================
  !> @brief Computes response variables for input displacements.
  !>
  !> @param[in]  disp  Displacement vector
  !> @param[in]  Ids   User IDs of the functions with response quantities
  !> @param[out] resp  Response variable values
  !> @param[in]  nDofs Length of displacement vector
  !> @param[in]  nIds  Length of function identifier array
  !> @param[out] ierr  Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Guenter Glanzer
  !>
  !> @date 6 June 2023

  subroutine computeResponseVars (disp,Ids,resp,nDofs,nIds,ierr)

    use SolverRoutinesModule, only : restoreLastStep
    use reportErrorModule   , only : reportError, debugFileOnly_p

    real(dp), intent(in)  :: disp(*)
    integer , intent(in)  :: Ids(*), nDofs, nIds
    real(dp), intent(out) :: resp(*)
    integer , intent(out) :: ierr

    !! Local variables
    integer  :: i
    real(dp) :: newValue

    !! --- Logic section ---

    !! Extract current response values
    ierr = 0
    do i = 1, nIds
       call getEngine (resp(i),ierr,Ids(i))
    end do
    if (ierr < 0) goto 915

    !! Calculate variables from displacement vector
    call variablesFromDisp (disp,ierr)
    if (ierr < 0) goto 915

    !! Extract response value changes due to the given displacement state
    do i = 1, nIds
       call getEngine (newValue,ierr,Ids(i))
       resp(i) = newValue - resp(i)
    end do
    if (ierr < 0) goto 915

    !! Reset all position variables from the stored previous values
    call restoreLastStep (sam,sys,mech)

915 continue
    call reportError (debugFileOnly_p,'computeResponseVars',ierr=ierr)

  contains

    !> @brief Updates the configuration by applying given displacement vector.
    subroutine variablesFromDisp (disp,ierr)

      use TriadTypeModule     , only : TransSysVecToGlobal, IncTriadsPos
      use TriadTypeModule     , only : clearTriadForces
      use MasterSlaveJointRoutinesModule, only : updateJoints, IncJointsVar
      use SupElRoutinesModule , only : IncSupElsGenDofs, updateSupElsStatic
      use SpringRoutinesModule, only : updateSprings
      use SolExtensionModule  , only : csExpand

      real(dp), intent(in)    :: disp(nDofs)
      integer , intent(inout) :: ierr

      !! Expand the equation-ordered solution vector to DOF-order
      call csExpand (sam,disp,sys%sinc(:,1))

      !! Update the sys%sinc(:,1) vector
      call TransSysVecToGlobal (mech%triads,sys%sinc(:,1))

      !! Update the triad(i)%ur(1:3,1:4) matrices
      call IncTriadsPos (mech%triads,sys%sinc(:,1))

      !! Update the sups(i)%genDOFs%ur(:) vectors
      call IncSupElsGenDofs (mech%sups,sys%sinc(:,1))

      !! Update all joint DOFs in the model
      call IncJointsVar (mech%joints,sys%sinc(:,1),mech%motions)

      !! All triad forces are set to zero
      call clearTriadForces (mech%triads)

      !! Update joint triad positions and associated constraint variables
      call updateJoints (mech%joints,mech%motions,ierr)

      !! Use Finit as update
      call updateSupElsStatic (mech%sups,mech%supLoads,mech%env, &
           &                   sys%time,sys%nIterThisStep,linInc=2,ierr=ierr)

      !! Update spring variables
      call updateSprings (mech%axialSprings,mech%joints,.false.,ierr)

      if (ierr < 0) call reportError (debugFileOnly_p,'variablesFromDisp')

    end subroutine variablesFromDisp

  end subroutine computeResponseVars


  !!============================================================================
  !> @brief Gets joint spring stiffness coefficients.
  !>
  !> @param[out] sprCoeff Spring stiffness coefficient
  !> @param[in]  bid      Base ID for spring
  !> @param[out] ierr     Error flag
  !>
  !> @callergraph
  !>
  !> @author Guenter Glanzer
  !>
  !> @date 20 June 2023

  subroutine getJointSpringStiffness (sprCoeff,bid,ierr)

    use MasterSlaveJointTypeModule, only : MasterSlaveJointType, GetPtrToId

    real(dp), intent(out)  :: sprCoeff(*)
    integer , intent(in)   :: bid
    integer , intent(out)  :: ierr

    !! Local variables
    type(MasterSlaveJointType), pointer :: joint
    integer                             :: j

    !! --- Logic section ---

    joint => getPtrToId(mech%joints,bid)
    if (.not. associated(joint)) then
       ierr = 1
       return
    else if (.not. associated(joint%springEl)) then
       ierr = 2
       return
    else
       ierr = 0
    end if

    sprCoeff(:size(joint%springEl%spr)) = 1.0_dp ! initialize spring coefficent

    !! Get joint spring stiffness
    do j = 1, size(joint%springEl%spr)
       if (associated(joint%springEl%spr(j)%p)) then
          sprCoeff(j) = joint%springEl%spr(j)%p%stiffness
       end if
    end do

  end subroutine getJointSpringStiffness

end module solverModule
