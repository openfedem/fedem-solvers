!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file initiateModule.f90
!> @brief Initialization of model data from the FEDEM solver input file.

!!==============================================================================
!> @brief Initialization of model data from the FEDEM solver input file.

module initiateModule

  implicit none

contains

  !!============================================================================
  !> @brief Administers reading of model data from the FEDEM solver input file.
  !>
  !> @param[in]  chfsi Text string containing the entire solver model
  !> @param[out] sys System level model data
  !> @param[out] ctrl Control system data
  !> @param[out] mech Mechanism components of the model
  !> @param[out] modelFileName Name of the FEDEM model file
  !> @param[out] doEnergyInt If .true., energy integration is requested
  !> @param[out] err Error flag
  !>
  !> @details This subroutine parses the FEDEM solver input file and initializes
  !> all solver objects based on it. The name of the solver input file is
  !> obtained from the command-line argument `-fsifile` (and `-fsi2file`).
  !> Alternatively, the model can be initialized by passing the entire input
  !> file content through the input variable @a cfsi.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 28 May 2002

  subroutine readSolverData (chfsi,sys,ctrl,mech,modelFileName,doEnergyInt,err)

    use kindModule           , only : dp, lfnam_p
    use IdTypeModule         , only : ReportInputError, getId
    use SystemTypeModule     , only : SystemType
    use ControlTypeModule    , only : ControlType, NullifyCtrl
    use ControlTypeModule    , only : ReadControlSystem
    use MechanismTypeModule  , only : MechanismType, NullifyMechanism
    use MechanismTypeModule  , only : ReadMechanism
    use EnvironmentTypeModule, only : ReadEnvironment, InitiateEnvironment
    use MassTypeModule       , only : InitiateMasses
    use MotionTypeModule     , only : InitiateMotions
    use ForceTypeModule      , only : InitiateForces
    use SupElLoadTypeModule  , only : InitiateLoads

    use headingNameListModule       , only : read_HEADING, modelFile, version
    use inputUtilities              , only : iuCopyToScratch, iuWriteToScratch
    use initiateSystemTypeModule    , only : InitiateSystem
    use initiateTriadTypeModule     , only : ReadTriads
    use initiateSupElTypeModule     , only : ReadSupEls
    use initiateSpringTypeModule    , only : InitiateBaseSprings
    use initiateSpringTypeModule    , only : InitiateSprings
    use initiateSpringTypeModule    , only : InitiateSpringFailures
    use initiateSpringTypeModule    , only : InitiateSpringYields
    use initiateSpringTypeModule    , only : CheckSpringConnections
    use initiateDamperTypeModule    , only : InitiateDampers
    use ContactSurfaceModule        , only : ReadGliderCurves
    use initiateJointTypeModule     , only : ReadJoints, ReadHigherPairs
    use initiateBushingElmTypeModule, only : ReadBushingElements
    use initiateContactElmTypeModule, only : ReadContactElements
    use initiateRoadTypeModule      , only : ReadRoads
    use initiateTireTypeModule      , only : ReadTires
    use initiateUserdefElTypeModule , only : ReadUserdefEls
    use initiateWindTurbineModule   , only : ReadTurbineConfig
    use initiateFunctionTypeModule  , only : InitiateFunctions, InitiateEngines1
    use initiateFrictionTypeModule  , only : InitiateFrictionData
    use initiateFrictionTypeModule  , only : InitFrictionPtrArray
    use initiateSensorTypeModule    , only : InitiateSensors
#ifdef FT_HAS_RECOVERY
    use StressRecoveryModule        , only : ReadStrainRosettes
#endif
    use FileUtilitiesModule         , only : findUnitNumber, getFileName
    use ReportErrorModule           , only : reportError
    use ReportErrorModule           , only : error_p, debugFileOnly_p

    character(len=*), optional, intent(in)  :: chfsi
    type(SystemType)          , intent(out) :: sys
    type(ControlType)         , intent(out) :: ctrl
    type(MechanismType)       , intent(out) :: mech
    character(len=*)          , intent(out) :: modelFileName
    logical                   , intent(out) :: doEnergyInt
    integer                   , intent(out) :: err

    !! Local variables
    integer            :: i, infp, tIncEngine
    real(dp)           :: weight(3)
    character(lfnam_p) :: chname


    !! --- Logic section ---

    call NullifyMechanism (mech) ! Avoid crash on file opening failure, etc.
    call NullifyCtrl (ctrl)

    !! Open a temporary solver input file

    infp = findUnitNumber(50)
    open(infp,IOSTAT=err,STATUS='SCRATCH')
    if (err /= 0) then
       call reportError (error_p,'Unable to open temporary solver input file', &
            &            addString='readSolverData')
       return
    end if

    if (present(chfsi)) then

       !! Write model definition to the temporary solver input file

       call iuWriteToScratch (infp,chfsi,err)
       if (err /= 0) goto 990

    else

       !! Copy the solver input file(s) into one temporary file

       call getFileName ('fsifile',chname)
       call iuCopyToScratch (infp,chname,err)
       if (err /= 0) goto 990

       call getFileName ('fsi2file',chname)
       if (chname /= '') then
          call iuCopyToScratch (infp,chname,err)
          if (err /= 0) goto 990
       end if

    end if
    rewind(infp)


    !! --- Start parsing the solver input file

    call Read_HEADING (infp,err)
    if (err /= 0) then
       call ReportInputError('HEADING')
       goto 990
    end if

    modelFileName = modelFile

    call ReadEnvironment (infp,mech%env,err)
    if (err /= 0) goto 990

    call ReadMechanism (infp,mech,tIncEngine,weight,err)
    if (err /= 0) goto 990

    !! Check if energy calculation is requested
    doEnergyInt = mech%saveVar(2)

    call ReadTriads (infp,mech%triads,err)
    if (err /= 0) goto 990

    call ReadSupEls (infp,mech%env,mech%triads,mech%sups,err)
    if (err /= 0) goto 990

    if (.not. doEnergyInt) then
       !! Check if energy calculation is requested
       do i = 1, size(mech%sups)
          if (mech%sups(i)%saveVar(3)) doEnergyInt = .true.
       end do
    end if

    call InitiateFunctions (infp,mech%gravity,mech%env%seaDepth, &
         &                  mech%env%waveFunc,mech%functions,err)
    if (err /= 0) goto 990

    call InitiateEngines1 (infp,mech%engines,mech%functions,err)
    if (err /= 0) goto 990

    call ReadUserdefEls (infp,mech%triads,mech%engines,mech%elms,err)
    if (err /= 0) goto 990

    call InitiateMasses (infp,mech%masses,mech%triads,mech%engines,err)
    if (err /= 0) goto 990

    call InitiateFrictionData (infp,mech%frictionSets,err)
    if (err /= 0) goto 990

    call InitiateSpringFailures (infp,mech%springFailures,err)
    if (err /= 0) goto 990

    call InitiateSpringYields (infp,mech%engines,mech%springYields,err)
    if (err /= 0) goto 990

    call InitiateBaseSprings (infp,mech%engines,mech%functions, &
         &                    mech%springFailures,mech%springYields, &
         &                    mech%baseSprings,err)
    if (err /= 0) goto 990

    call InitiateSprings (infp,mech%triads, &
         &                mech%baseSprings,mech%axialSprings,err)
    if (err /= 0) goto 990

    call InitiateDampers (infp,mech%triads,mech%engines,mech%functions, &
         &                mech%axialSprings,mech%baseDampers,mech%dampers,err)
    if (err /= 0) goto 990

    call ReadGliderCurves (infp,mech%triads,mech%cSurfs,err)
    if (err /= 0) goto 990

    call ReadJoints (infp,mech%triads,mech%cSurfs, &
         &           mech%baseSprings,mech%baseDampers,mech%dampers, &
         &           mech%engines,mech%joints,mech%frictionSets,err)
    if (err /= 0) goto 990

    call ReadHigherPairs (infp,mech%joints,mech%higherPairs,err)
    if (err /= 0) goto 990

    call InitiateLoads (infp,mech%sups,mech%engines,mech%supLoads,err)
    if (err /= 0) goto 990

    call InitiateForces (infp,mech%triads,mech%joints,mech%sups,mech%engines, &
         &               mech%forces,err)
    if (err /= 0) goto 990

    call getFileName ('displacementfile',chname)
    call InitiateMotions (infp,chname,mech%triads,mech%joints,mech%engines, &
         &                mech%motions,err)
    if (err /= 0) goto 990

    if (.not. doEnergyInt) then
       !! Check if energy calculation is requested
       do i = 1, size(mech%supLoads)
          if (mech%supLoads(i)%saveVar(3)) doEnergyInt = .true.
       end do
       do i = 1, size(mech%forces)
          if (mech%forces(i)%saveVar(3)) doEnergyInt = .true.
       end do
       do i = 1, size(mech%motions)
          if (mech%motions(i)%saveVar(3)) doEnergyInt = .true.
       end do
    end if

    call ReadBushingElements (infp,mech%triads, &
         &                    mech%baseSprings,mech%baseDampers,mech%bElems,err)
    if (err /= 0) goto 990

    call ReadContactElements (infp,mech%triads,mech%cSurfs,mech%frictionSets, &
         &                    mech%baseSprings,mech%baseDampers,mech%cElems,err)
    if (err /= 0) goto 990

    call ReadRoads (infp,mech%functions,mech%engines,mech%roads,err)
    if (err /= 0) goto 990

    call ReadTires (infp,mech%joints,mech%roads,mech%tires,err)
    if (err /= 0) goto 990

    if (.not. doEnergyInt) then
       !! Check if energy calculation is requested
       do i = 1, size(mech%tires)
          if (mech%tires(i)%saveVar(9)) doEnergyInt = .true.
       end do
    end if

    !! Check that all base springs are connected to something
    call CheckSpringConnections (mech%baseSprings,err)
    if (err < 0) goto 990

    !! Establish array of pointers to all frictions
    call InitFrictionPtrArray (mech%joints,mech%cElems,mech%frictions,err)
    if (err /= 0) goto 990

    if (.not. doEnergyInt) then
       !! Check if energy calculation is requested
       do i = 1, size(mech%baseSprings)
          if (mech%baseSprings(i)%saveVar(5)) doEnergyInt = .true.
       end do
       do i = 1, size(mech%baseDampers)
          if (mech%baseDampers(i)%saveVar(5)) doEnergyInt = .true.
       end do
       do i = 1, size(mech%frictions)
          if (mech%frictions(i)%p%saveVar(2)) doEnergyInt = .true.
       end do
    end if

    !! Save input and external energy only if the mechanism has that
    mech%hasEext = size(mech%tires) > 0
    mech%hasEinp = size(mech%supLoads) > 0 .or. &
         &         size(mech%forces) > 0 .or. size(mech%motions) > 0
    if (.not. mech%hasEinp) then
       do i = 1, size(mech%baseSprings)
          if (associated(mech%baseSprings(i)%length0Engine)) then
             mech%hasEinp = .true.
          end if
       end do
    end if

    call InitiateSensors (infp,sys,mech,mech%sensors,err)
    if (err /= 0) goto 990

    call ReadControlSystem (infp,mech%engines,mech%sensors,ctrl,err)
    if (err /= 0) goto 990

    call ReadTurbineConfig (infp,mech%env,mech%turbine, &
         &                  mech%triads,mech%sups,mech%joints,err)
    if (err /= 0) goto 990

    if (associated(mech%turbine)) mech%hasEext = .true.

#ifdef FT_HAS_RECOVERY
    call ReadStrainRosettes (infp,err)
    if (err /= 0) goto 990
#endif

    close(infp)

    !! Get solution algorithm parameters from command-line arguments
    call InitiateSystem (sys,mech%engines,tIncEngine,weight,err)
    if (err /= 0) goto 999

    !! Connect the sea environment functions
    call InitiateEnvironment (mech%env,mech%functions,mech%engines,err)
    if (err /= 0) goto 999

    return

990 close(infp)
999 call reportError (debugFileOnly_p,'readSolverData')

  end subroutine readSolverData


  !!============================================================================
  !> @brief Administers preprocessing of model data for the dynamics solver.
  !>
  !> @param[out] sam Data for managing system matrix assembly
  !> @param      sys System level model data
  !> @param[in]  ctrl Control system data
  !> @param[out] modes Data for eigenmodes
  !> @param      mech Mechanism components of the model
  !> @param[in]  ipsw Print switch; the higher value the more print is produced
  !> @param[out] err Error flag
  !>
  !> @details This subroutine builds up the necessary data structures and
  !> resolves cross references, etc., after the solver input file has been read.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 28 May 2002

  subroutine preprocessSolverData (sam,sys,ctrl,modes,mech,ipsw,err)

    use SamModule          , only : SamType, WriteObject
    use ModesTypeModule    , only : ModesType
    use SystemTypeModule   , only : SystemType, WriteObject
    use ControlTypeModule  , only : ControlType, InitiateControl, WriteObject
    use MechanismTypeModule, only : MechanismType, WriteObject
    use MotionTypeModule   , only : InitiateMotions2

    use SamSolverModule           , only : initSAM, initSAM_dofStatus
    use SamSolverModule           , only : initSAM_topology, initSAM_constraints
    use SamSolverModule           , only : initSAM_preassembly
    use initiateSystemTypeModule  , only : InitiateSystem2
    use initiateModesTypeModule   , only : InitiateModes
    use initiateTriadTypeModule   , only : InitiateTriads2
    use initiateSupElTypeModule   , only : InitiateSupEls2
    use initiateJointTypeModule   , only : InitiateJoints, InitiateHigherPairs
    use initiateTireTypeModule    , only : InitiateTires
    use initiateFunctionTypeModule, only : InitiateEngines2
    use initiateSensorTypeModule  , only : InitiateSensors2
    use initiateUserdefElTypeModule,only : InitiateUserdefEls
    use genericPartModule         , only : InitiateGenericParts
    use EngineRoutinesModule      , only : SetPointersForSensors
    use FNVwaveForceModule        , only : InitiateWaveSpectrumFNV
    use FNVwaveForceModule        , only : InitiateFNVkinematics, calcFNVforces
    use FileUtilitiesModule       , only : findUnitNumber
    use ReportErrorModule         , only : reportError, error_p, debugFileOnly_p

    type(SamType)      , intent(out)   :: sam
    type(ModesType)    , intent(out)   :: modes
    type(SystemType)   , intent(inout) :: sys
    type(ControlType)  , intent(in)    :: ctrl
    type(MechanismType), intent(inout) :: mech
    integer            , intent(in)    :: ipsw
    integer            , intent(out)   :: err

    !! Local variables
    integer :: idata


    !! --- Logic section ---

    idata = 0
    if (ipsw > 0) then
       idata = findUnitNumber(10)
       open(idata,FILE='datastructure.lst',IOSTAT=err)
       if (err /= 0) then
          idata = 0 ! Open failure, silently ignore
       else if (ipsw > 2) then
          write(idata,"('* Dump of solver data structure after data input *'/)")

          !! Dump all input data to datastructure.lst
          call WriteObject (sys,idata,ipsw)
          write(idata,'(/)')

          call WriteObject (mech,idata,-1)
          write(idata,'(/)')

          call WriteObject (ctrl,idata,1)
          write(idata,'(//)')
          call flush(idata)

       end if
    end if

    !! Generate wave spectrum data for FNV force calculation
    call InitiateWaveSpectrumFNV (sys,mech%env,err)
    if (err /= 0) goto 900

    call initSAM (sam,err)
    if (err /= 0) goto 990

    !! Establish the SAM arrays MADOF and MSC (nodal DOF data)
    call initSAM_dofStatus (mech%triads, mech%sups, mech%joints, &
         &                  sam%mpar, sam%madof, sam%msc, sam%dofType, err)
    if (err /= 0) goto 980

    if (sam%ndof < 1) then
       err = 1
       call reportError (error_p,'No DOFs in this model, nothing to solve')
       goto 980
    end if

    !! Establish the SAM arrays MPMNPC and MMNPC (element topology)
    call initSAM_topology (mech%sups, mech%axialSprings, mech%dampers, &
         &                 mech%joints, mech%bElems, mech%cElems, mech%elms, &
         &                 mech%masses, mech%tires, sam%mpar, &
         &                 sam%mmnpc, sam%mpmnpc, err)
    if (err /= 0) goto 980

    if (sam%nel < 1) then
       err = 2
       call reportError (error_p,'No elements in this model, nothing to solve')
       goto 980
    end if

    !! Establish the SAM arrays MPMCEQ, MMCEQ and TTCC (linear dependencies)
    call initSAM_constraints (mech%motions, mech%joints, mech%higherPairs, &
         &                    sam%mpar, sam%madof, sam%msc, &
         &                    sam%mpmceq, sam%mmceq, sam%ttcc, sam%mpreac, err)
    if (err /= 0) goto 980

    !! Initialize the linear equation solver (nodal reordering, preassembly, ..)
    call initSAM_preassembly (sam,sys%Nmat,err)
    if (err /= 0) goto 900

    !! Allocate system matrices and vectors
    call InitiateSystem2 (sam,sys,err)
    if (err /= 0) goto 900

    !! Initialize the eigenvalue solver
    call InitiateModes (sam,sys,modes,err)
    if (err /= 0) goto 900

    !! Initialize control parameters in mpar
    call InitiateControl (sam%mpar,ctrl)

    !! Initialize pointers from triads into system vectors
    call InitiateTriads2 (sam,mech%triads,err)
    if (err /= 0) goto 900

    !! Initialize pointers from superelements into system vectors
    call InitiateSupEls2 (sam,sys,mech%sups,mech%masses,mech%engines,err)
    if (err /= 0) goto 900

    !! Initialize linear dependencies in joints and associated ICs
    call InitiateJoints (mech%joints,mech%forces,mech%motions, &
         &               sam%mpreac,sys%RFk,err)
    if (err /= 0) goto 900

    !! Resolve initial conditions in joints due to higher pairs
    call InitiateHigherPairs (mech%higherPairs)

    !! Initialize generic part velocities (at CoG triads)
    call InitiateGenericParts(mech%sups)

    !! Initialize user-defined elements
    call InitiateUserdefEls (mech%gravity,mech%elms,err)
    if (err /= 0) goto 990

    !! Initialize pointers from prescribed motions to reaction forces
    call InitiateMotions2 (mech%motions,sam%mpreac,sys%RFk,err)
    if (err /= 0) goto 990

    !! Initialize pointers from sensors to measured global variables
    call InitiateSensors2 (mech%sensors,mech%triads,mech%joints,ctrl,err)
    if (err /= 0) goto 900

    !! Initialize pointers from engines to argument sensors/engines
    call InitiateEngines2 (mech%engines,mech%sensors,err)
    if (err /= 0) goto 900

    !! Initialize global pointers to mechanism objects for updateSensor
    call SetPointersForSensors (mech%engines,mech%triads, &
         &                      mech%baseSprings,mech%baseDampers)

    !! Initialize the tire models
    call InitiateTires (sys,mech%tires,mech%roads,err)
    if (err /= 0 .or. sys%tEnd <= sys%tStart) goto 900

    !! Initialize the FNV force histories
    call InitiateFNVkinematics (sys%tEnd-sys%tStart,sys%tInc, &
         &                      mech%env%seaDepth,err)
    if (err /= 0) goto 900

    call calcFNVforces (sys%tEnd-sys%tStart,mech%env%rhow,err)

900 continue
    if (idata < 1) goto 990

    !! Dump all input data to datastructure.lst
    write(idata,"('* Dump of solver data structure after preprocessing *'/)")

    call WriteObject (sys,idata,ipsw)
    write (idata,'(/)')

    call WriteObject (mech,idata,-1)
    write (idata,'(/)')

    call WriteObject (ctrl,idata,1)
    write (idata,'(/)')

980 continue
    if (idata < 1) goto 990

    !! Dump SAM data to datastructure.lst
    call WriteObject (sam,idata,3)

    close(idata)

990 continue
    if (err /= 0) call reportError (debugFileOnly_p,'preprocessSolverData')

  end subroutine preprocessSolverData

end module initiateModule
