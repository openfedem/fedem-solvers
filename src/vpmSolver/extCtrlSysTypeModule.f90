!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module ExtCtrlSysTypeModule

  use KindModule        , only : dp, lfnam_p
  use IdTypeModule      , only : IdType, ldesc_p
  use FunctionTypeModule, only : EnginePtrType
  use SensorTypeModule  , only : SensorPtrType

  implicit none

  ! Define element types and enums
  ! ------------------------------
  character(len=6), parameter :: ExtCtrlSysType_p(2) = (/ &
       &                         'MATLAB', &
       &                         'EASY5 '  /)

  integer         , parameter :: MATLAB_p = 1, &
       &                         EASY5_p  = 2


  ! Define the external control system type
  ! ---------------------------------------
  type ExtCtrlSysType
     type(IdType)                    :: id
     integer                         :: sysType
     integer                         :: timeDim
     character(len=lfnam_p)          :: path, simCmd
     type(EnginePtrType)   , pointer :: engineIn(:)
     type(SensorPtrType)   , pointer :: sensorOut(:)
     real(dp)              , pointer :: state(:)
     character(len=ldesc_p), pointer :: match(:)
  end type ExtCtrlSysType


  interface GetPtrToId
     module procedure GetPtrToIdExtCtrlSys
  end interface


contains

  function GetPtrToIdExtCtrlSys (array,id) result(ptr)
    type(ExtCtrlSysType), pointer :: ptr

    !!==========================================================================
    !! Return pointer to (first) element with specified id or NULL if not found.
    !!
    !! Programmer : Torleif Iversen
    !! date/rev   : 16 Aug 2002
    !!==========================================================================

    integer             , intent(in)        :: id
    type(ExtCtrlSysType), intent(in),target :: array(:)

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
    write(*,*) '*** GetPtrToIdExtCtrlSys returned nullified, baseId =',id

  end function GetPtrToIdExtCtrlSys


  subroutine NullifyExtCtrlSys (extCtrlSys)

    !!==========================================================================
    !! Initialize the extCtrlSysType object.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 6 Jan 2003/3.0
    !!==========================================================================

    use IdTypeModule, only : nullifyId

    type(ExtCtrlSysType), intent(out) :: extCtrlSys

    !! --- Logic section ---

    call nullifyId (extCtrlSys%id)
    extCtrlSys%sysType = 0
    extCtrlSys%timeDim = -1
    extCtrlSys%path    = ''
    extCtrlSys%simCmd  = ''

    nullify(extCtrlSys%engineIn)
    nullify(extCtrlSys%sensorOut)
    nullify(extCtrlSys%state)
    nullify(extCtrlSys%match)

  end subroutine NullifyExtCtrlSys


  subroutine DeallocateExtCtrlSys (extCtrlSys)

    !!==========================================================================
    !! Deallocate the extCtrlSysType object.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 23 Jan 2017/1.0
    !!==========================================================================

    use IdTypeModule, only : deallocateId

    type(ExtCtrlSysType), intent(out) :: extCtrlSys

    !! --- Logic section ---

    call deallocateId (extCtrlSys%id)

    if (associated(extCtrlSys%engineIn))  deallocate(extCtrlSys%engineIn)
    if (associated(extCtrlSys%sensorOut)) deallocate(extCtrlSys%sensorOut)
    if (associated(extCtrlSys%state))     deallocate(extCtrlSys%state)
    if (associated(extCtrlSys%match))     deallocate(extCtrlSys%match)

    call nullifyExtCtrlSys (extCtrlSys)

  end subroutine DeallocateExtCtrlSys


  subroutine ReadExtCtrlSysInfo (infp,engines,sensors,extCtrlSys,err)

    !!==========================================================================
    !! Read external control system information from the fsi file,
    !! create external control system types and assign values
    !!
    !! Programmer : Haavar Johan Johnsen
    !! Revised    : Torleif Iversen
    !! date/rev   : June 2002/2.0
    !!==========================================================================

    use FunctionTypeModule, only : EngineType, GetPtrToId
    use SensorTypeModule  , only : SensorType, GetPtrToId, MATLAB_WS_p
    use IdTypeModule      , only : initId, ReportInputError
    use inputUtilities    , only : iuGetNumberOfEntries
    use inputUtilities    , only : iuSetPosAtNextEntry, iuCharToInt
    use reportErrorModule , only : AllocationError
    use reportErrorModule , only : reportError, debugFileOnly_p

    integer             , intent(in)  :: infp
    type(EngineType)    , intent(in)  :: engines(:)
    type(SensorType)    , intent(in)  :: sensors(:)
    type(ExtCtrlSysType), pointer     :: extCtrlSys(:)
    integer             , intent(out) :: err

    !! Local variables
    character(len=ldesc_p) :: errMsg
    integer                :: i, elm, numExtCtrlSys, numInputs, numOutputs, stat
    integer                :: sensorOutID(100)

    ! Set up the namelist for the external control systems
    ! ----------------------------------------------------
    integer                :: id, extId(10), engineInID(100)
    character(len=ldesc_p) :: extDescr, match(100)
    character(len=8)       :: sysType
    character(len=lfnam_p) :: fileName, path
    namelist /EXTERNAL_CONTROL_SYSTEM/ id, extId, extDescr, &
         &   engineInID, sysType, match, fileName, path


    !! --- Logic section ---

    numExtCtrlSys = iuGetNumberOfEntries(infp,'&EXTERNAL_CONTROL_SYSTEM',err)
    if (err /= 0 .or. numExtCtrlSys == 0) return

    ! Allocate space for all external control systems found in the fsi file
    ! ---------------------------------------------------------------------
    allocate(extCtrlSys(numExtCtrlSys),STAT=stat)
    if (stat /= 0) then
       err = AllocationError('ReadExtCtrlSysInfo')
       return
    end if

    do elm = 1, numExtCtrlSys

       ! Step to next external control system
       ! ------------------------------------
       call NullifyExtCtrlSys (extCtrlSys(elm))
       if (.not. iuSetPosAtNextEntry(infp,'&EXTERNAL_CONTROL_SYSTEM')) then
          err = err - 1
          call ReportInputError ('EXTERNAL_CONTROL_SYSTEM',elm)
          cycle
       end if

       ! Initialize variables in the namelist
       ! ------------------------------------
       id = 0; extId = 0; extDescr = ''
       sysType = ''; fileName = ''; path = ''; match = ''; engineInID = 0

       ! Read namelist values from the fsi file
       ! --------------------------------------
       read(infp,nml=EXTERNAL_CONTROL_SYSTEM,iostat=stat)
       if (stat /= 0) then
          err = err - 1
          call ReportInputError ('EXTERNAL_CONTROL_SYSTEM',elm)
          cycle
       end if

       ! Count the number of inputs to the external control system
       ! ---------------------------------------------------------
       numInputs = count(engineInID > 0)

       ! Get sensors of type Matlab
       ! --------------------------
       sensorOutID = 0
       numOutputs = 0
       do i = 1, size(sensors)
          if (sensors(i)%type == MATLAB_WS_p .and. sensors(i)%index(1) == id) then
             numOutputs = numOutputs + 1
             sensorOutID(numOutputs) = sensors(i)%id%baseId
          end if
       end do

       ! Count the number of outputs from the external control system
       ! ------------------------------------------------------------
       numOutputs = count(sensorOutID > 0)

       ! Allocate space for the external control system data
       ! ---------------------------------------------------
       allocate(extCtrlSys(elm)%engineIn(numInputs), &
            &   extCtrlSys(elm)%match(numInputs), &
            &   extCtrlSys(elm)%sensorOut(numOutputs), STAT=stat)
       if (stat /= 0) then
          err = AllocationError('ReadExtCtrlSysInfo')
          return
       end if

       ! Transfer the fsi file values to the created external control system
       ! -------------------------------------------------------------------
       call initId (extCtrlSys(elm)%id,id,extId,extDescr,stat)
       extCtrlSys(elm)%simCmd  = fileName
       extCtrlSys(elm)%sysType = iuCharToInt(sysType,ExtCtrlSysType_p)
       extCtrlSys(elm)%match   = match

       ! Set up a link between the external control system and the engine
       ! ----------------------------------------------------------------
       do i = 1, size(extCtrlSys(elm)%engineIn)
          extCtrlSys(elm)%engineIn(i)%p => GetPtrToId(engines,engineInID(i))
          if (.not. associated(extCtrlSys(elm)%engineIn(i)%p)) then
             err = err - 1
             write(errMsg,*) 'Unconnected input ID ',engineInID(i)
             call ReportInputError ('EXTERNAL_CONTROL_SYSTEM', &
                  &                 elm,extCtrlSys(elm)%id,errMsg)
          end if
       end do

       ! Set up a link between the external control system and the sensor
       ! ----------------------------------------------------------------
       do i = 1, size(extCtrlSys(elm)%sensorOut)
          extCtrlSys(elm)%sensorOut(i)%p => GetPtrToId(sensors,sensorOutID(i))
          if (.not.associated(extCtrlSys(elm)%sensorOut(i)%p)) then
             err = err - 1
             write(errMsg,*) 'Unconnected output ID ',sensorOutID(i)
             call ReportInputError ('EXTERNAL_CONTROL_SYSTEM', &
                  &                 elm,extCtrlSys(elm)%id,errMsg)
          end if
       end do

    end do

    if (err < 0) call reportError (debugFileOnly_p,'ReadExtCtrlSysInfo')

  end subroutine ReadExtCtrlSysInfo

end module ExtCtrlSysTypeModule
