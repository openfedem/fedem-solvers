!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file initiateSensorTypeModule.f90
!> @brief Initialization of sensor objects from the solver input file.

!!==============================================================================
!> @brief Initialization of sensor objects from the solver input file.

module InitiateSensorTypeModule

  implicit none

  private :: ReportSensorError


contains

  !> @cond NO_DOCUMENTATION
  subroutine ReportSensorError (sensor,errMsg,errId)
    use SensorTypeModule, only : SensorType
    use IdTypeModule    , only : StrId, ReportInputError
    type(SensorType), intent(in) :: sensor
    character(len=*), intent(in) :: errMsg
    integer         , intent(in) :: errId
    call ReportInputError ('SENSOR',id=sensor%id,msg=errMsg//StrId(errId))
  end subroutine ReportSensorError
  !> @endcond


  !!============================================================================
  !> @brief Allocates the time sensor object.
  !>
  !> @param[in] sys System level model data
  !> @param[out] err Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 22 Aug 2022

  subroutine AllocateTimeSensor (sys,err)

    use SystemTypeModule , only : SystemType
    use SensorTypeModule , only : nullifySensor, ourTime, TIME_p
    use reportErrorModule, only : AllocationError

    type(SystemType), target, intent(in)  :: sys
    integer                 , intent(out) :: err

    !! --- Logic section ---

    err = 0
    if (associated(ourTime)) return

    allocate(ourTime,STAT=err)
    if (err /= 0) then
       err = AllocationError('AllocateTimeSensor')
       return
    end if

    call nullifySensor (ourTime)
    ourTime%type = TIME_p
    ourTime%value => sys%time

  end subroutine AllocateTimeSensor


  !!============================================================================
  !> @brief Initializes sensor objects with data from the solver input file.
  !>
  !> @param[in] infp File unit number for the solver input file
  !> @param[in] sys System level model data
  !> @param[in] mech Mechanism components of the model
  !> @param[out] sensors Array of all sensor objects in the model
  !> @param[out] err Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen                                       @date 7 Sep 1999
  !> @author Knut Morten Okstad                                 @date 6 Jun 2002

  subroutine InitiateSensors (infp,sys,mech,sensors,err)

    use SystemTypeModule          , only : SystemType
    use MechanismTypeModule       , only : MechanismType, EnvironmentType
    use TriadTypeModule           , only : TriadType, GetPtrToId
    use MasterSlaveJointTypeModule, only : MasterSlaveJointType, GetPtrToId
    use ContactElementTypeModule  , only : ContactElementType, GetPtrToId
    use SpringTypeModule          , only : SpringBaseType, GetPtrToId
    use DamperTypeModule          , only : DamperBaseType, GetPtrToId
    use FunctionTypeModule        , only : EngineType, GetPtrToId
    use EngineRoutinesModule      , only : ourEnvir
    use SensorTypeModule          ! using everything from this module
    use IdTypeModule              , only : initId, ReportInputError, ldesc_p
    use inputUtilities            , only : iuGetNumberOfEntries
    use inputUtilities            , only : iuSetPosAtNextEntry, iuCharToInt
    use progressModule            , only : lterm
    use reportErrorModule         , only : AllocationError
    use reportErrorModule         , only : reportError, debugFileOnly_p

    integer                    , intent(in)  :: infp
    type(SystemType)   , target, intent(in)  :: sys
    type(MechanismType), target, intent(in)  :: mech
    type(SensorType)           , pointer     :: sensors(:), tmp(:)
    integer                    , intent(out) :: err

    !! Local variables
    integer :: i, ii, ldof, nSensors, stat
    logical :: needT, nullifyOnError

    type(TriadType)           , pointer :: triad
    type(EngineType)          , pointer :: engine
    type(SpringBaseType)      , pointer, save :: spring
    type(DamperBaseType)      , pointer, save :: damper
    type(MasterSlaveJointType), pointer, save :: joint
    type(ContactElementType)  , pointer, save :: cElem

    character(len=30)  :: type, dofEntity, dofSystem
    character(ldesc_p) :: extDescr
    integer            :: id, extId(10), triadId(4), dof, engineId, &
         &                ctrlVarId, damperId, springId, jointId

    namelist /SENSOR/ id, extId, extDescr, triadId, dof, &
         &            engineId, ctrlVarId, damperId, springId, &
         &            type, dofEntity, dofSystem, jointId

    !! --- Logic section ---

    nSensors = iuGetNumberOfEntries(infp,'&SENSOR',err)
    if (err /= 0 .or. nSensors < 1) goto 900
    write(lterm,*) 'Number of &SENSOR =',nSensors

    allocate(tmp(nSensors),STAT=stat)
    if (stat == 0) then
       call DeallocateSensors (sensors)
       sensors => tmp
    else
       err = AllocationError('InitiateSensors')
       return
    end if

    needT = .false.
    do ii = 1, nSensors

       nullifyOnError = .false.
       call NullifySensor (sensors(ii))
       if (.not. iuSetPosAtNextEntry(infp,'&SENSOR')) then
          err = err - 1
          call ReportInputError ('SENSOR',ii)
          cycle
       end if

       !! Initialize to default
       id = 0; extId = 0; extDescr = ''; triadId = 0; dof = 0
       engineId = 0; ctrlVarId = 0; damperId = 0; springId = 0
       type = ''; dofEntity = ''; dofSystem = ''

       !! Read values
       read(infp,nml=SENSOR,iostat=stat)
       if (stat /= 0) then
          err = err - 1
          call ReportInputError ('SENSOR',ii)
          cycle
       end if

       call initId (sensors(ii)%id,id,extId,extDescr,stat)
       sensors(ii)%type   = iuCharToInt(type,      sensorType_p)
       sensors(ii)%entity = iuCharToInt(dofEntity, sensorEntity_p)
       sensors(ii)%system = iuCharToInt(dofSystem, sensorSystem_p)
       sensors(ii)%dof    = dof

       select case (sensors(ii)%type)

       case (TRIAD_p)

          call allocIndex (sensors(ii)%index,1)
          triad => GetPtrToId(mech%triads,triadId(1),sensors(ii)%index(1))
          if (.not. associated(triad)) then
             err = err - 1
             call ReportSensorError (sensors(ii), &
                  'Non-existing triad, baseId =',triadId(1))
          end if

          if (sensors(ii)%entity <= DYN_P_p) then
             if (sensors(ii)%entity >= W_SPEED_p) needT = .true.
             if (sensors(ii)%entity >= F_VEL_p) ourEnvir => mech%env
          end if

          !! The rest to be initiated in InitiateSensors2


       case (RELATIVE_TRIAD_p)

          if (triadId(4) > 0) then
             call allocIndex (sensors(ii)%index,4)
          else if (triadId(2) > 0) then
             call allocIndex (sensors(ii)%index,2)
          else
             call allocIndex (sensors(ii)%index,0)
          end if
          do i = 1, size(sensors(ii)%index)
             triad => GetPtrToId(mech%triads,triadId(i),sensors(ii)%index(i))
             if (.not. associated(triad)) then
                err = err - 1
                call ReportSensorError (sensors(ii), &
                     'Non-existing triad, baseId =',triadId(i))
             end if
          end do

          select case (sensors(ii)%dof)

          case (0:13) ! 0 is relative, 1-6 are the usual DOFs

          case default
             err = err - 1
             call ReportSensorError (sensors(ii), &
                  'Invalid relative sensor DOF:',sensors(ii)%dof)
          end select


       case (JOINT_VARIABLE_p)

          call allocIndex (sensors(ii)%index,1)
          joint => GetPtrToId(mech%joints,jointId,sensors(ii)%index(1))
          if (associated(joint)) then
             if (associated(joint%slider) .and. sensors(ii)%dof == 3) then
                ldof = 7 ! Slider dof
             else
                ldof = sensors(ii)%dof
             end if

             do dof = 1, size(joint%jointDofs)
                if (joint%jointDofs(dof)%ldof == ldof) then
                   select case (sensors(ii)%entity)

                   case (REL_POS_p)
                      sensors(ii)%value => joint%jointDofs(dof)%jvar(1)

                   case (VEL_p)
                      sensors(ii)%value => joint%jointDofs(dof)%jvar(2)

                   case (ACC_p)
                      sensors(ii)%value => joint%jointDofs(dof)%jvar(3)

                   case (FORCE_p)
                      sensors(ii)%value => joint%jointDofs(dof)%jvar(1) ! tmp

                   case default
                      ldof = 0
                      err = err - 1
                      call ReportSensorError (sensors(ii), &
                           'Invalid sensor entity:',sensors(ii)%entity)
                   end select

                   exit
                end if
             end do

          else

             call allocIndex (sensors(ii)%index,1)
             cElem => GetPtrToId(mech%cElems,jointId,sensors(ii)%index(1))
             if (associated(celem)) then
                ldof = sensors(ii)%dof

                if (ldof >= 1 .and. ldof <= 3) then
                   select case (sensors(ii)%entity)

                   case (REL_POS_p)
                      sensors(ii)%value => cElem%cVar(ldof,1)

                   case (VEL_p)
                      sensors(ii)%value => cElem%cVar(ldof,2)

                   case default
                      ldof = 0
                      err = err - 1
                      call ReportSensorError (sensors(ii), &
                           'Invalid sensor entity:',sensors(ii)%entity)
                   end select
                end if

             else
                err = err - 1
                call ReportSensorError (sensors(ii), &
                     'Non-existing joint, baseId =',jointId)
                if (nullifyOnError) nullify(sensors(ii)%index)
                cycle
             end if

          end if

          if (.not. associated(sensors(ii)%value) .and. ldof > 0) then
             err = err - 1
             call ReportSensorError (sensors(ii), &
                  'Invalid joint sensor DOF:',sensors(ii)%dof)
          end if


       case (TIME_p)

          sensors(ii)%value => sys%time
          ourTime => sensors(ii) ! Used by wind speed and sea state sensors


       case (NUM_ITERATIONS_p)

          sensors(ii)%value => sys%nIterPrevStep


       case (ENGINE_p)

          call allocIndex (sensors(ii)%index,1)
          engine => GetPtrToId(mech%engines,engineId,.false., &
               &               sensors(ii)%index(1))
          if (.not. associated(engine)) then
             err = err - 1
             call ReportSensorError (sensors(ii), &
                  'Non-existing engine, baseId =',engineId)
             cycle
          end if


       case (CONTROL_p)

          sensors(ii)%dof = ctrlVarId

          !! The rest to be initialized in InitiateSensors2


#ifdef FT_HAS_RECOVERY
       case (STRAIN_GAGE_p)

          call allocIndex (sensors(ii)%index,1)
          sensors(ii)%index(1) = engineId

          !! The rest to be initialized in InitiateSensors2
#endif


       case (SPRING_AXIAL_p, SPRING_JOINT_p)

          call allocIndex (sensors(ii)%index,1)
          spring => GetPtrToId(mech%baseSprings,springId,sensors(ii)%index(1))
          if (.not. associated(spring)) then
             err = err - 1
             call ReportSensorError (sensors(ii), &
                  'Non-existing spring, baseId =',springId)
             cycle
          end if

          select case (sensors(ii)%entity)

          case (LENGTH_p)
             sensors(ii)%value => spring%length
          case (DEFL_p)
             sensors(ii)%value => spring%deflection
          case (FORCE_p)
             sensors(ii)%value => spring%force
          case default
             err = err - 1
             call ReportSensorError (sensors(ii), &
                  'Invalid spring sensor entity:',sensors(ii)%entity)
          end select


       case (DAMPER_AXIAL_p, DAMPER_JOINT_p)

          call allocIndex (sensors(ii)%index,1)
          damper => GetPtrToId(mech%baseDampers,damperId,sensors(ii)%index(1))
          if (.not. associated(damper)) then
             err = err - 1
             call ReportSensorError (sensors(ii), &
                  'Non-existing damper, baseId =',damperId)
             cycle
          end if

          select case (sensors(ii)%entity)

          case (LENGTH_p)
             sensors(ii)%value => damper%length
          case (VEL_p)
             sensors(ii)%value => damper%velocity
          case (FORCE_p)
             sensors(ii)%value => damper%force
          case default
             err = err - 1
             call ReportSensorError (sensors(ii), &
                  'Invalid damper sensor entity:',sensors(ii)%entity)
          end select


       case default
          err = err - 1
          call ReportSensorError (sensors(ii), &
               'Invalid sensor type:',sensors(ii)%type)
       end select

       if (nullifyOnError) nullify(sensors(ii)%index)

    end do

    if (err == 0 .and. needT) then
       call allocateTimeSensor (sys,err)
    end if

900 if (err /= 0) call reportError (debugFileOnly_p,'InitiateSensors')

  contains

    !> @cond NO_DOCUMENTAION
    subroutine allocIndex (indices,nInd)
      integer, pointer    :: indices(:)
      integer, intent(in) :: nInd
      integer, save, target :: dummy(1)
      allocate(indices(nInd),STAT=stat)
      if (stat == 0) then
         if (nInd > 0) indices = 0
      else
         err = err + AllocationError('InitiateSensors')
         if (nInd == 1) then
            indices => dummy
            nullifyOnError = .true.
         else
            nullify(indices)
         end if
      end if
    end subroutine allocIndex
    !> @endcond

  end subroutine InitiateSensors


  !!============================================================================
  !> @brief Initializes pointers to the measured quantity for sensor objects.
  !>
  !> @param sensors Array of all sensor objects in the model
  !> @param[in] triads All triads in the model
  !> @param[in] joints Array of all joints in the model
  !> @param[in] ctrl Control system data
  !> @param[out] err Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 7 Mar 2001

  subroutine InitiateSensors2 (sensors,triads,joints,ctrl,err)

    use SensorTypeModule          ! using everything from this module
    use TriadTypeModule           , only : TriadType
    use MasterSlaveJointTypeModule, only : MasterSlaveJointType
    use ControlTypeModule         , only : ControlType
#ifdef FT_HAS_RECOVERY
    use StrainRosetteModule       , only : StrainRosetteType
    use StressRecoveryModule      , only : getStrainRosette
#endif
    use reportErrorModule         , only : AllocationError
    use reportErrorModule         , only : reportError, debugFileOnly_p

    type(SensorType)          , intent(inout)      :: sensors(:)
    type(TriadType)           , intent(in), target :: triads(:)
    type(MasterSlaveJointType), intent(in), target :: joints(:)
    type(ControlType)         , intent(in)         :: ctrl
    integer                   , intent(out)        :: err

    !! Local variables
    integer :: ii, dof, ldof

    type(TriadType)           , pointer :: triad
    type(MasterSlaveJointType), pointer :: joint
#ifdef FT_HAS_RECOVERY
    type(StrainRosetteType), pointer, save :: strGage
#endif

    !! -- Logic section ---

    err = 0
    do ii = 1, size(sensors)

       dof = sensors(ii)%dof
       select case (sensors(ii)%type)

       case (CONTROL_p)

          if (dof < 1 .or. dof > size(ctrl%vreg)) then
             err = err - 1
             call ReportSensorError (sensors(ii), &
                  'Invalid control variable, id =',dof)
          else
             sensors(ii)%value => ctrl%vreg(dof)
          end if

       case (TRIAD_p)

          select case (sensors(ii)%entity)

          case (POS_p)

             select case (sensors(ii)%system)

             case (GLOBAL_p)

                select case (dof)

                case (1:3)
                   sensors(ii)%value => triads(sensors(ii)%index(1))%ur(dof,4)

                case (4:9)
                   continue

                case default
                   err = err - 1
                   call ReportSensorError (sensors(ii), &
                        'Invalid position DOF:',dof)
                end select

             case default
                err = err - 1
                call ReportSensorError (sensors(ii), &
                     'Invalid position coordinate system:',sensors(ii)%system)
             end select

          case (VEL_p)

             select case (dof)

             case (1:6)

                select case (sensors(ii)%system)

                case (GLOBAL_p)
                   sensors(ii)%value => triads(sensors(ii)%index(1))%urd(dof)

                case (LOCAL_p)
                   continue

                case default
                   err = err - 1
                   call ReportSensorError (sensors(ii), &
                        'Invalid velocity coordinate system:', &
                        sensors(ii)%system)
                 end select

             case default
                err = err - 1
                call ReportSensorError (sensors(ii), &
                     'Invalid velocity DOF:',dof)
             end select

          case (ACC_p)

             select case (dof)

             case (1:6)

                select case (sensors(ii)%system)

                case (GLOBAL_p)
                   sensors(ii)%value => triads(sensors(ii)%index(1))%urdd(dof)

                case (LOCAL_p)
                   continue

                case default
                   err = err - 1
                   call ReportSensorError (sensors(ii), &
                        'Invalid acceleration coordinate system:', &
                        sensors(ii)%system)
                end select

             case default
                err = err - 1
                call ReportSensorError (sensors(ii), &
                     'Invalid acceleration DOF:',dof)
             end select

          case (W_SPEED_p:DYN_P_p)
             continue

          case (FORCE_p)

             !! Ensure the force array is allocated (in case it is not plotted)
             triad => triads(sensors(ii)%index(1))
             if (triad%isFEnode .or. triad%isSEnode) then
                if (.not. associated(triad%supForce)) then
                   allocate(triad%supForce(triad%nDOFs),STAT=err)
                   if (err /= 0) then
                      err = AllocationError('InitiateSensors1')
                      return
                   end if
                   triad%supForce = 0.0_dp
                end if
             else if (.not. associated(triad%nodeForce)) then
                err = err - 1
                call ReportSensorError (sensors(ii), &
                     'No nodal force computed for Triad',triad%id%userId)
             end if

          case default
             err = err - 1
             call ReportSensorError (sensors(ii), &
                  'Invalid triad sensor entity:',sensors(ii)%entity)
          end select

       case (JOINT_VARIABLE_p)

          if (sensors(ii)%entity == FORCE_p) then

             joint => joints(sensors(ii)%index(1))
             if (associated(joint%slider) .and. dof == 3) dof = 7 ! Slider dof

             do ldof = 1, size(joint%jointDofs)
                if (joint%jointDofs(ldof)%ldof == dof) then
                   if (associated(joint%jointDofs(ldof)%F)) then
                      sensors(ii)%value => joint%jointDofs(ldof)%F
                   end if
                   exit
                end if
             end do

          end if

#ifdef FT_HAS_RECOVERY
       case (STRAIN_GAGE_p)

          strGage => getStrainRosette(sensors(ii)%index(1))
          if (.not. associated(strGage)) then
             err = err - 1
             call ReportSensorError (sensors(ii), &
                  'Invalid strain rosette ID:',sensors(ii)%index(1))
             cycle
          else if (dof < 1 .or. dof > 4+size(strGage%gages)) then
             err = err - 1
             call ReportSensorError (sensors(ii), &
                  'Invalid strain gage DOF index:',dof)
             cycle
          end if

          select case (sensors(ii)%entity)

          case (STRAIN_p)

             select case (dof)
             case (1:3)
                sensors(ii)%value => strGage%epsP(dof)
             case (4)
                sensors(ii)%value => strGage%epsVM
             case default
                sensors(ii)%value => strGage%gages(dof-4)%epsGage
             end select

          case (STRESS_p)

             select case (dof)
             case (1:3)
                sensors(ii)%value => strGage%sigmaP(dof)
             case (4)
                sensors(ii)%value => strGage%sigmaVM
             case default
                sensors(ii)%value => strGage%gages(dof-4)%sigGage
             end select

          case default
             err = err - 1
             call ReportSensorError (sensors(ii), &
                  'Invalid strain gage sensor entity:',sensors(ii)%entity)

          end select
#endif

       end select

       if (err == 0 .and. .not. associated(sensors(ii)%value)) then
          allocate(sensors(ii)%value,STAT=err)
          if (err /= 0) then
             err = AllocationError('InitiateSensors2')
             return
          end if
          sensors(ii)%allocValue = .true.
          sensors(ii)%value = 0.0_dp
       end if

    end do

    if (err < 0) call reportError (debugFileOnly_p,'InitiateSensors2')

  end subroutine InitiateSensors2

end module InitiateSensorTypeModule
