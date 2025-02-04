!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file controlTypeModule.f90
!> @brief Control system data container.

!!==============================================================================
!> @brief Module with data types representing control system object of a model.
!>
!> @details The module also contains subroutines for initializing the control
!> system from the input file, as well as general subroutines and functions for
!> for accessing the control system data.

module ControlTypeModule

  use KindModule          , only : dp
  use IdTypeModule        , only : IdType
  use FunctionTypeModule  , only : EngineType
  use SensorTypeModule    , only : SensorType
#ifdef FT_HAS_EXTCTRL
  use ExtCtrlSysTypeModule, only : ExtCtrlSysType
#endif

  implicit none

  private

  !> @brief Data type representing a control input parameter.
  type CtrlPrm
     integer                   :: var    !< Variable index of the input
     type(EngineType), pointer :: engine !< Function processing the input
     type(SensorType), pointer :: sensor !< Sensor meansuring the input
  end type CtrlPrm

  !> @brief Data type representing the control system of a model.
  type ControlType

     type(CtrlPrm), pointer :: input(:)  !< Control input parameters
     real(dp)     , pointer :: vreg(:)   !< The control state variables
     type(IdType) , pointer :: vregId(:) !< Id for the control line variables
     integer      , pointer :: ivar(:)   !< Control variables for extra lines
     logical                :: saveVar   !< .true. if variables should be saved

     integer           :: mpireg(9)  !< Matrix of pointers into IREG
     integer           :: mprreg(13) !< Matrix of pointers into RREG
     integer , pointer :: ireg(:)    !< Integer buffer for control system
     real(dp), pointer :: rreg(:)    !< Real buffer for control system
     real(dp), pointer :: delay(:)   !< Buffer for delay elements

#ifdef FT_HAS_EXTCTRL
     type(ExtCtrlSysType), pointer :: extCtrlSys(:) !< External control systems
#endif

  end type ControlType


  !> @brief Standard routine for writing an object to file.
  interface WriteObject
     module procedure WriteControlType
  end interface

  public :: ControlType, nullifyCtrl, deallocateCtrl
  public :: ReadControlSystem, InitiateControl, hasControlElements
  public :: WriteObject, writeCtrlSysHeader, writeCtrlSysDB


contains

  !!============================================================================
  !> @brief Initializes the ControlType object.
  !>
  !> @param[out] ctrl Control system data
  !> @param[in] deallocating If .true., the pointers are nullified
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 4 Jul 2004

  subroutine nullifyCtrl (ctrl,deallocating)

    type(ControlType), intent(out) :: ctrl
    logical, optional, intent(in)  :: deallocating

    !! --- Logic section ---

    ctrl%mpireg = 0
    ctrl%mprreg = 0
    ctrl%saveVar = .false.

    if (present(deallocating)) then
       nullify(ctrl%input)
    else
       !! Allocate initially to zero length, such that the size() operator works
       allocate(ctrl%input(0))
    end if

    nullify(ctrl%vreg)
    nullify(ctrl%vregId)
    nullify(ctrl%ivar)
    nullify(ctrl%ireg)
    nullify(ctrl%rreg)
    nullify(ctrl%delay)
#ifdef FT_HAS_EXTCTRL
    nullify(ctrl%extCtrlSys)
#endif

  end subroutine nullifyCtrl


  !!============================================================================
  !> @brief Deallocates the ControlType object.
  !>
  !> @param ctrl Control system data
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Jan 2017

  subroutine deallocateCtrl (ctrl)

    use IdTypeModule        , only : deallocateId
    use AllocationModule    , only : reAllocate
#ifdef FT_HAS_EXTCTRL
    use ExtCtrlSysTypeModule, only : deallocateExtCtrlSys
#endif

    type(ControlType), intent(inout) :: ctrl

    !! Local variables
    integer :: i

    !! --- Logic section ---

    if (associated(ctrl%input)) deallocate(ctrl%input)
    if (associated(ctrl%vreg))  deallocate(ctrl%vreg)
    if (associated(ctrl%vregId)) then
       do i = 1, size(ctrl%vregId)
          call deallocateId (ctrl%vregId(i))
       end do
       deallocate(ctrl%vregId)
    end if
    if (associated(ctrl%ivar))  deallocate(ctrl%ivar)
    if (associated(ctrl%ireg))  deallocate(ctrl%ireg)
    if (associated(ctrl%rreg))  deallocate(ctrl%rreg)
    call reAllocate ('deallocateCtrl',ctrl%delay)

#ifdef FT_HAS_EXTCTRL
    if (associated(ctrl%extCtrlSys)) then
       do i = 1, size(ctrl%extCtrlSys)
          call deallocateExtCtrlSys (ctrl%extCtrlSys(i))
       end do
       deallocate(ctrl%extCtrlSys)
    end if
#endif

    call nullifyCtrl (ctrl,.true.)

  end subroutine deallocateCtrl


  !!============================================================================
  !> @brief Initializes the ControlType object with data from the input file.
  !>
  !> @param[in] infp File unit number for the solver input file
  !> @param[in] engines All general functions in the model
  !> @param[in] sensors All sensors (function argument objects) in the model
  !> @param[out] ctrl Control system data of the model
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date Jul 2000

  subroutine ReadControlSystem (infp,engines,sensors,ctrl,ierr)

    use inputUtilities        , only : iuGetNumberOfEntries
    use allocationModule      , only : reAllocate
    use IdTypeModule          , only : nullifyId
#ifdef FT_HAS_EXTCTRL
    use ExtCtrlSysTypeModule  , only : ReadExtCtrlSysInfo
#endif
    use progressModule        , only : lterm
    use reportErrorModule     , only : reportError, debugFileOnly_p
    use reportErrorModule     , only : AllocationError, GetErrorFile
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getdoubles
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getbool
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getint

    integer          , intent(in)  :: infp
    type(EngineType) , intent(in)  :: engines(:)
    type(SensorType) , intent(in)  :: sensors(:)
    type(ControlType), intent(out) :: ctrl
    integer          , intent(out) :: ierr

    !! Local variables
    logical :: allSec, allRest, allCtrl
    integer :: i, iprint, nCEl, nCIn, nCLin, nCVar, nCV2, nRreg, nIreg, nDelay


    !! --- Logic section ---

    nCVar = iuGetNumberOfEntries(infp,'&CONTROL_VARIABLE',ierr)
    if (ierr /= 0) goto 999
    write(lterm,*) 'Number of &CONTROL_VARIABLE =',nCVar
    nCLin = iuGetNumberOfEntries(infp,'&CONTROL_LINE',ierr)
    if (ierr /= 0) goto 999
    write(lterm,*) 'Number of &CONTROL_LINE =',nCLin
    nCIn  = iuGetNumberOfEntries(infp,'&CONTROL_INPUT',ierr)
    if (ierr /= 0) goto 999
    write(lterm,*) 'Number of &CONTROL_INPUT =',nCIn
    nCEl  = iuGetNumberOfEntries(infp,'&CONTROL_ELEMENT',ierr)
    if (ierr /= 0) goto 999
    write(lterm,*) 'Number of &CONTROL_ELEMENT =',nCEl

    if (nCVar+nCIn+nCEl < 1) goto 900 ! No control system at all

    !! Initialize pointer arrays
    ctrl%mpireg(1) = 1                      ! 1: Control integer parameters
    ctrl%mpireg(2) = ctrl%mpireg(1)+10      ! 2: Variable status codes
    ctrl%mpireg(3) = ctrl%mpireg(2)+nCVar   ! 3: Module where each var. is comp.
    ctrl%mpireg(4) = ctrl%mpireg(3)+nCVar   ! 4: Pivot indicies for BNEW
    ctrl%mpireg(5) = ctrl%mpireg(4)+nCVar   ! 5: Pivot indicies for LNEW
    ctrl%mpireg(6) = ctrl%mpireg(5)+2*nCVar ! 6: Pointers to MTOP mat's in MMTOP
    ctrl%mpireg(7) = ctrl%mpireg(6)+nCEl+1  ! 7: Pointers to RPAR mat's in RRPAR
    ctrl%mpireg(8) = ctrl%mpireg(7)+nCEl+1  ! 8: MMTOP: element topology
    ctrl%mpireg(9) = ctrl%mpireg(8)+25*nCEl ! Assume average # elm-nodes <=25

    nCV2 = nCVar*nCVar
    ctrl%mprreg(1) = 1                       !  1: Control real parameters
    ctrl%mprreg(2) = ctrl%mprreg(1)+10       !  2: Previous control variables
    ctrl%mprreg(3) = ctrl%mprreg(2)+nCVar    !  3: Num. solution method 1
    ctrl%mprreg(4) = ctrl%mprreg(3)+nCVar    !  4: Num. solution method 2
    ctrl%mprreg(5) = ctrl%mprreg(4)+2*nCVar  !  5: Function values
    ctrl%mprreg(6) = ctrl%mprreg(5)+nCVar    !  6: Help vector
    ctrl%mprreg(7) = ctrl%mprreg(6)+4*nCVar  !  7: Local error tolerance
    ctrl%mprreg(8) = ctrl%mprreg(7)+2*nCVar  !  8: Estimated truncation error
    ctrl%mprreg(9) = ctrl%mprreg(8)+nCVar    !  9: Numerical Jacobian
    ctrl%mprreg(10)= ctrl%mprreg(9)+nCV2     ! 10: Iteration matrix for meth. 1
    ctrl%mprreg(11)= ctrl%mprreg(10)+nCV2    ! 11: Iteration matrix for meth. 2
    ctrl%mprreg(12)= ctrl%mprreg(11)+4*nCV2  ! 12: Element parameters RPAR
    ctrl%mprreg(13)= ctrl%mprreg(12)+25*nCEl ! Assume average # parameters <=25

    !! Allocate the control system arrays
    nDelay = 0
    nIreg  = ctrl%mpireg(9)  - 1
    nRreg  = ctrl%mprreg(13) - 1
    deallocate(ctrl%input)
    allocate(ctrl%input(nCIn), ctrl%ireg(nIreg), ctrl%rreg(nRreg), &
         &   ctrl%vreg(nCVar), ctrl%ivar(nCLin), ctrl%vregId(nCVar+nCLin), &
         &   STAT=ierr)
    if (ierr /= 0) then
       ierr = AllocationError('ReadControlSystem')
       return
    end if
    if (nIreg > 0) ctrl%ireg = 0
    if (nCLin > 0) ctrl%ivar = 0
    if (nRreg > 0) ctrl%rreg = 0.0_dp
    if (nCVar > 0) ctrl%vreg = 0.0_dp
    do i = 1, size(ctrl%vregId)
       call NullifyId (ctrl%vregId(i))
    end do

    !! Initialize control system parameters from command-line
    call ffa_cmdlinearg_getdoubles ('ctrlTol',ctrl%rreg(4:7),4)
    call ffa_cmdlinearg_getint ('debug',iprint)

    !! Determine if control variables will be saved (either all or nothing)
    call ffa_cmdlinearg_getbool ('allSecondaryVars',allSec)
    call ffa_cmdlinearg_getbool ('allRestartVars',allRest)
    call ffa_cmdlinearg_getbool ('allControlVars',allCtrl)
    ctrl%saveVar = allSec .or. allRest .or. allCtrl

    !! Read and initialize CONTROL_VARIABLEs
    if (nCVar > 0) then
       call ReadControlVariables (infp, &
            ctrl%ireg(ctrl%mpireg(2):ctrl%mpireg(3)-1), &
            ctrl%vreg, ctrl%vregId, ctrl%ivar, iprint, GetErrorFile(), ierr)
       if (ierr < 0) goto 999
    end if

    !! Read and initialize CONTROL_INPUTs
    if (nCIn > 0) then
       call ReadControlInput (infp, engines, sensors, ctrl%input, &
            ctrl%ireg(ctrl%mpireg(2):ctrl%mpireg(3)-1), ierr)
       if (ierr < 0) goto 999
    end if

    !! Read and initialize the CONTROL_ELEMENTs
    if (nCEl > 0) then
       call ReadControlElements (infp, nCEl, nDelay, &
            ctrl%ireg(ctrl%mpireg(6):ctrl%mpireg(7)-1),   & ! mpmtop
            ctrl%ireg(ctrl%mpireg(7):ctrl%mpireg(8)-1),   & ! mprpar
            ctrl%ireg(ctrl%mpireg(8):ctrl%mpireg(9)-1),   & ! mmtop
            ctrl%rreg(ctrl%mprreg(12):ctrl%mprreg(13)-1), & ! rpar
            iprint, GetErrorFile(), ierr)
       if (ierr < 0) goto 999
    end if

    !! Allocate initial delay buffer (may be reallocated later)
    if (nDelay > 0) then
       call ffa_cmdlinearg_getint ('delayBuffer',ctrl%ireg(10))
       nDelay = ctrl%ireg(10)
    else
       nDelay = 1 ! Still allocate one element, such that ctrl%delay(1) works
    end if
    call reAllocate ('ReadControlSystem',ctrl%delay,nDelay,ierr)
    if (ierr < 0) return
    ctrl%delay = 0.0_dp
900 continue
#ifdef FT_HAS_EXTCTRL
    call ReadExtCtrlSysInfo (infp, engines, sensors, ctrl%extCtrlSys, ierr)
    if (ierr == 0) return
#else
    return
#endif

999 call reportError (debugFileOnly_p,'ReadControlSystem')

  end subroutine ReadControlSystem


  !!============================================================================
  !> @brief Initializes the SAM control array MPAR with control system data.
  !>
  !> @param mpar Matrix of parameters
  !> @param[in] ctrl Control system data
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date Apr 2003

  subroutine InitiateControl (mpar,ctrl)

    integer          , intent(inout) :: mpar(:)
    type(ControlType), intent(in)    :: ctrl

    !! --- Logic section ---

    mpar(45) = size(ctrl%input)                    ! nCIn
    mpar(46) = 0                                   ! nCOut
    if (associated(ctrl%vreg)) then
       mpar(47) = size(ctrl%vreg)                  ! nCVar
    else
       mpar(47) = 0
    end if
    mpar(48) = ctrl%mpireg(7) - ctrl%mpireg(6) - 1 ! nCEl

    !! Set array lengths
    if (associated(ctrl%ireg)) then
       mpar(60) = ctrl%ireg(ctrl%mpireg(6)+mpar(48)) - 1
       mpar(61) = ctrl%ireg(ctrl%mpireg(7)+mpar(48)) - 1
    end if

  end subroutine InitiateControl


  !!============================================================================
  !> @brief Reads data for the control inputs from the solver input file.
  !>
  !> @param[in] infp File unit number for the solver input file
  !> @param[in] engines All general functions in the model
  !> @param[in] sensors All sensors (function argument objects) in the model
  !> @param[out] input Control input parameters
  !> @param[out] mstat Status flags for the control variables
  !> @param[out] ierr Error flag
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date Aug 2000

  subroutine ReadControlInput (infp,engines,sensors,input,mstat,ierr)

    use IdTypeModule      , only : ldesc_p, initId, ReportInputError
    use FunctionTypeModule, only : GetPtrToId
    use SensorTypeModule  , only : GetPtrToId
    use inputUtilities    , only : iuSetPosAtNextEntry
    use reportErrorModule , only : reportError, debugFileOnly_p

    integer         , intent(in)  :: infp
    type(EngineType), intent(in)  :: engines(:)
    type(SensorType), intent(in)  :: sensors(:)
    type(CtrlPrm)   , intent(out) :: input(:)
    integer         , intent(out) :: mstat(:)
    integer         , intent(out) :: ierr

    !! Local variables
    integer      :: idIn, stat
    type(IdType) :: cid

    character(len=ldesc_p) :: extDescr
    integer                :: id, extId(10), iVar, inEngineID, inSensorId
    namelist /CONTROL_INPUT/ id, extId, extDescr, iVar, inEngineID, inSensorId

    !! --- Logic section ---

    ierr = 0
    rewind(infp)

    do idIn = 1, size(input)
       input(idIn)%var = 0
       nullify(input(idIn)%engine)
       nullify(input(idIn)%sensor)
       if (.not. iuSetPosAtNextEntry(infp,'&CONTROL_INPUT')) then
          ierr = ierr - 1
          call ReportInputError ('CONTROL_INPUT',idIn)
          cycle
       end if

       !! Default values
       id=0; extId=0; extDescr=''
       iVar=0; inEngineID=0; inSensorId=0

       read(infp,nml=CONTROL_INPUT,iostat=stat)
       if (stat /= 0) then
          ierr = ierr - 1
          call ReportInputError ('CONTROL_INPUT',idIn)
          cycle
       end if

       call initId (cid,id,extId,extDescr,stat)

       !! Set the input variable number
       if (iVar > 0 .and. iVar <= size(mstat)) then
          input(idIn)%var = iVar
          mstat(iVar) = 1 ! Signals input variable
       else
          ierr = ierr - 1
          call ReportInputError ('CONTROL_INPUT',idIn,cid)
          cycle
       end if

       !! Set the type of input: engine or sensor
       if (inEngineId > 0) then
          input(idIn)%engine => GetPtrToId(engines,inEngineId,.false.)
          if (associated(input(idIn)%engine)) cycle
       else if (inSensorId > 0) then
          input(idIn)%sensor => GetPtrToId(sensors,inSensorId)
          if (associated(input(idIn)%sensor)) cycle
       end if

       ierr = ierr - 1
       call ReportInputError ('CONTROL_INPUT',idIn,cid)

    end do

    if (ierr < 0) call reportError (debugFileOnly_p,'ReadControlInput')

  end subroutine ReadControlInput


  !!============================================================================
  !> @brief Reads data for the control elements from the solver input file.
  !>
  !> @param[in] infp File unit number for the solver input file
  !> @param[in] nCEl Total number of control elements in the model
  !> @param[out] nDelay Number of delay elements in the model
  !> @param[out] mpmtop Matrix of pointers to topology vectors
  !> @param[out] mprpar Matrix of pointers to real data for the control elements
  !> @param[out] mmtop Topology vectors for the control elements
  !> @param[out] rpar Real data for the control elements
  !> @param[in] iprint Print switch; the higher value the more print is produced
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date Aug 2000

  subroutine ReadControlElements (infp,   nCEl,   nDelay, &
       &                          mpmtop, mprpar, mmtop, rpar, &
       &                          iprint, lpu,    ierr )

    use IdTypeModule     , only : ldesc_p, initId, ReportInputError
    use inputUtilities   , only : iuSetPosAtNextEntry
    use reportErrorModule, only : reportError, debugFileOnly_p

    integer , intent(in)  :: infp, nCEl, iprint, lpu
    integer , intent(out) :: nDelay, MPMTOP(:), MPRPAR(:), MMTOP(:)
    real(dp), intent(out) :: RPAR(:)
    integer , intent(out) :: ierr

    !! Local variables
    integer      :: i, iVar, nVar, ix, stat
    type(IdType) :: cid

    logical, external :: checkCtrlParams

    character(len=ldesc_p) :: extDescr
    integer                :: id, extId(10), type, nRealData, variables(100)
    real(dp)               :: realData(100)
    namelist /CONTROL_ELEMENT/ id, extId, extDescr, type, nRealData, &
         &                     realData, variables

    !! --- Logic section ---

    ierr = 0
    nDelay = 0

    if (iprint > 0) write(LPU,"(/1X,'CONTROL ELEMENT PARAMETERS')")

    MPRPAR(1) = 1
    MPMTOP(1) = 1

    do i = 1, nCEl

       if (.not. iuSetPosAtNextEntry(infp,'&CONTROL_ELEMENT')) then
          ierr = ierr - 1
          call ReportInputError ('CONTROL_ELEMENT',i)
          cycle
       end if

       !! Default values
       id=0; extId=0; extDescr=''
       type=0; nRealData=0; realData=0.0_dp; variables=0

       read(infp,nml=CONTROL_ELEMENT,iostat=stat)
       if (stat /= 0) then
          ierr = ierr - 1
          call ReportInputError ('CONTROL_ELEMENT',i)
          cycle
       end if

       call initId (cid,id,extId,extDescr,stat)

       nVar = 0 ! Number of variables (DOFS)
       do iVar = 1, size(variables)
          if (variables(iVar) == 0) then
             nVar = iVar - 1
             exit
          end if
       end do

       !! Check validity of element type number and the number of
       !! parameters and variables for the specified element type

       if (.not. checkCtrlParams(type,nRealData,nVar)) then
          call ReportInputError ('CONTROL_ELEMENT',i,cid)
          if (type > 0) then
             write(LPU,600) nRealData,nVar,type
600          format(10X,'Illegal number of parameters (NO)',I3,', or' &
                  &/10X,'illegal number of connection nodes (IC)',I3, &
                  &/10X,'for control element of type (ETYP)',I3)
          else
             write(LPU,"( 10X,'Undefined type (ETYP)',I3)") -type
          end if
          ierr = ierr - 1
          cycle
       end if

       !! Store type and realData in RPAR array

       ix          = mprpar(i)
       mprpar(i+1) = mprpar(i)+nRealData+1 ! +1 since type is also stored here
       rpar(ix)    = type
       rpar(ix+1:ix+nRealData) = realData(1:nRealData)
       if (iprint > 0) then
          write(lpu,"(3X,'elno.',I3,' type',I3,' realdata',100(1X,E15.5,1X))") &
               &     i,type,rpar(ix+1:ix+nRealData)
       end if

       !! Read module topology data

       ix = mpmtop(i)-1
       mpmtop(i+1) = mpmtop(I)+nVar

       do iVar = 1, nVar
          mmtop(ix+iVar) = abs(variables(iVar))
       end do

       if (iprint > 0) then
          write(LPU,"(12X,'topology ',100(1X,I3,1X))") mmtop(ix+1:ix+nVar)
       end if

       if (type == 11) nDelay = nDelay + 1

    end do

    if (ierr < 0) call reportError (debugFileOnly_p,'ReadControlElements')

  end subroutine ReadControlElements


  !!============================================================================
  !> @brief Reads status codes and initial conditions for the control variables.
  !>
  !> @param[in] infp File unit number for the solver input file
  !> @param[out] mstat Status flags for the control variables
  !> @param[out] vreg The control state variables
  !> @param[out] vregId Id for the control line variables
  !> @param[out] mvar Control variables for extra control lines
  !> @param[in] iprint Print switch; the higher value the more print is produced
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @details This subroutine reads status codes and initial conditions
  !> for all control variables in the model from the solver input file.
  !> It also reads the Id of other control lines referring to each variable.
  !> This is needed such that the control variable value can be written to the
  !> results database for all control lines that refer to it (Bugfix #323).
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !> @date Aug 2000
  !>
  !> @author Knut Morten Okstad
  !> @date 29 Apr 2017

  subroutine ReadControlVariables (infp,mstat,vreg,vregId,mvar,iprint,lpu,ierr)

    use IdTypeModule     , only : ldesc_p,  initId, ReportInputError
    use inputUtilities   , only : iuSetPosAtNextEntry
    use reportErrorModule, only : reportError, debugFileOnly_p

    integer     , intent(in)  :: infp, iprint, lpu
    integer     , intent(out) :: mstat(:), mvar(:)
    real(dp)    , intent(out) :: vreg(:)
    type(IdType), intent(out) :: vregId(:)
    integer     , intent(out) :: ierr

    !! Local variables
    integer :: i, stat

    character(len=ldesc_p) :: extDescr
    integer                :: id, extId(10), iVar, status
    real(dp)               :: initValue
    namelist /CONTROL_VARIABLE/ id, extId, extDescr, iVar, status, initValue
    namelist /CONTROL_LINE/     id, extId, extDescr, iVar

    !! --- Logic section ---

    ierr = 0
    rewind(infp)

    if (iprint > 0) then
       write(lpu,"(/1X,'CONTROL STATE VARIABLE INITIAL VALUE AND STATUS CODE')")
    end if

    do i = 1, size(vreg)

       if (.not. iuSetPosAtNextEntry(infp,'&CONTROL_VARIABLE')) then
          ierr = ierr - 1
          call ReportInputError ('CONTROL_VARIABLE',i)
          cycle
       end if

       !! Default values
       id=0; extId=0; extDescr=''; iVar=0; status=0; initValue=0.0_dp

       read(infp,nml=CONTROL_VARIABLE,iostat=stat)
       if (stat /= 0 .or. iVar < 1 .or. iVar > size(vreg)) then
          ierr = ierr - 1
          call ReportInputError ('CONTROL_VARIABLE',i)
          cycle
       end if

       call initId (vregId(iVar),id,extId,extDescr,stat)
       vreg(iVar)  = initValue
       mstat(iVar) = status
       if (iprint > 0 .and. extId(1) > 0) then
          write(lpu,"(1X,I3,1X,I5,2X,E15.5,5X,I3)") iVar,id,initValue,status
       end if

    end do

    rewind(infp)
    status = size(vreg)

    do i = 1, size(mvar)

       if (.not. iuSetPosAtNextEntry(infp,'&CONTROL_LINE')) then
          ierr = ierr - 1
          call ReportInputError('CONTROL_LINE',i)
          cycle
       end if

       !! Default values
       id=0; extId=0; extDescr=''; iVar=0

       read(infp,nml=CONTROL_LINE,iostat=stat)
       if (stat /= 0 .or. iVar < 1 .or. iVar > size(vreg)) then
          ierr = ierr - 1
          call ReportInputError('CONTROL_LINE',i)
       else
          status = status + 1
          call initId (vregId(status),id,extId,extDescr,stat)
          mvar(i) = iVar
       end if

    end do

    if (ierr < 0) call reportError (debugFileOnly_p,'ReadControlVariables')

  end subroutine ReadControlVariables


  !!============================================================================
  !> @brief Standard routine for writing an object to io.
  !>
  !> @param[in] ctrl Control system data
  !> @param[in] io File unit number to write to
  !> @param[in] complexity If present, the value indicates the amount of print
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 27 Jul 2000

  subroutine WriteControlType (ctrl,io,complexity)

    use IdTypeModule     , only : getId
    use manipMatrixModule, only : writeObject

    type(ControlType), intent(in) :: ctrl
    integer          , intent(in) :: io
    integer, optional, intent(in) :: complexity

    !! Local variables
    integer :: i

    !! --- Logic section ---

    write(io,'(A)') 'Control','{'

    do i = 1, size(ctrl%input)
       if (associated(ctrl%input(i)%engine)) then
          write(io,*) 'input var', ctrl%input(i)%var, &
               ': engineId =', trim(getId(ctrl%input(i)%engine%id))
       else if (associated(ctrl%input(i)%sensor)) then
          write(io,*) 'input var', ctrl%input(i)%var, &
               ': sensorId =', trim(getId(ctrl%input(i)%sensor%id))
       end if
    end do

    if (present(complexity)) then
       if (complexity >= 1) then
          write(io,*) 'size(input)  =', size(ctrl%input)
          write(io,*) 'size(mpireg) =', size(ctrl%mpireg)
          write(io,*) 'size(mprreg) =', size(ctrl%mprreg)
          if (associated(ctrl%ireg)) then
             write(io,*) 'size(ireg)   =', size(ctrl%ireg)
          end if
          if (associated(ctrl%rreg)) then
             write(io,*) 'size(rreg)   =', size(ctrl%rreg)
          end if
          if (associated(ctrl%delay)) then
             write(io,*) 'size(delay)  =', size(ctrl%delay)
          end if
#ifdef FT_HAS_EXTCTRL
          if (associated(ctrl%extCtrlSys)) then
             write(io,*) 'size(extCtrlSys) =', size(ctrl%extCtrlSys)
          end if
#endif
       end if
       if (complexity >= 2) then
          write(io,'(A,15I6)') ' mpireg =', ctrl%mpireg
          write(io,'(A,15I6)') ' mprreg =', ctrl%mprreg
          if (associated(ctrl%ireg)) then
             write(io,*)
             call writeIREG (1,ctrl%mpireg,'ICNT')
             call writeIREG (2,ctrl%mpireg,'MSTAT')
             call writeIREG (3,ctrl%mpireg,'MMOD')
             call writeIREG (4,ctrl%mpireg,'IPVT1')
             call writeIREG (5,ctrl%mpireg,'IPVT2')
             call writeIREG (6,ctrl%mpireg,'MPTOP')
             call writeIREG (7,ctrl%mpireg,'MPRPAR')
             call writeIREG (8,ctrl%mpireg,'MMTOP')
          end if
          if (associated(ctrl%rreg)) then
             write(io,*)
             call writeRREG (1,ctrl%mprreg,'RCNT')
             call writeRREG (2,ctrl%mprreg,'VPREG')
             call writeRREG (3,ctrl%mprreg,'VBREG')
             call writeRREG (4,ctrl%mprreg,'VLREG')
             call writeRREG (5,ctrl%mprreg,'VDREG')
             call writeRREG (6,ctrl%mprreg,'VHREG')
             call writeRREG (7,ctrl%mprreg,'TOL')
             call writeRREG (8,ctrl%mprreg,'ERR')
             call writeRREG (9,ctrl%mprreg,'AJAC')
             call writeRREG(10,ctrl%mprreg,'ABE')
             call writeRREG(11,ctrl%mprreg,'ALO')
             call writeRREG(12,ctrl%mprreg,'RPAR')
          end if
          if (associated(ctrl%delay)) then
             write(io,*)
             call writeObject (ctrl%delay,io,'DELAY')
          end if
       end if
    end if

    write(io,'(A)') '}'

  contains

    !> @brief Writes an indexed integer array to io.
    subroutine writeIREG (i,mp,name)
      integer         , intent(in) :: i, mp(:)
      character(len=*), intent(in) :: name
      call writeObject (ctrl%ireg(mp(i):mp(i+1)-1),io,name)
    end subroutine writeIREG

    !> @brief Writes an indexed real array to io.
    subroutine writeRREG (i,mp,name)
      integer         , intent(in) :: i, mp(:)
      character(len=*), intent(in) :: name
      call writeObject (ctrl%rreg(mp(i):mp(i+1)-1),io,name)
    end subroutine writeRREG

  end subroutine WriteControlType


  !!============================================================================
  !> @brief Checks if control elements are present in the model.
  function hasControlElements (ctrl)
    type(ControlType), intent(in) :: ctrl
    logical :: hasControlElements
    hasControlElements = size(ctrl%input) > 0
  end function hasControlElements


  !!============================================================================
  !> @brief Writes results database headers for the control system.
  !>
  !> @param[in] ctrl Control system data
  !> @param[in] mechId Id of the mechanism object
  !> @param rdb Results database file for control system data
  !>
  !> @details This subroutine saves only those arrays that need to be
  !> restored in a restart. The list of arrays stored *must* match that in
  !> controltypemodule::writectrlsysdb(), and should at least contain
  !> the arrays referred in restartmodule::readcontroldata().
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 25 Oct 2008

  subroutine writeCtrlSysHeader (ctrl,mechId,rdb)

    use IdTypeModule, only : IdType, writeIdHeader, StrId
    use RDBModule   , only : RDBType, idatd, nbi_p, nbd_p, writeVarDef

    type(ControlType), intent(in)    :: ctrl
    type(IdType)     , intent(in)    :: mechId
    type(RDBType)    , intent(inout) :: rdb

    !! Local variables
    integer :: nV
    integer :: iIVar(10) = 0
    integer :: iRVar(15) = 0

    !! --- Logic section ---

    call writeIdHeader ('Mechanism',mechId,idatd)

    call writeIRDef (1,ctrl%mpireg,'ICNT')
    !call writeIRDef (2,ctrl%mpireg,'MSTAT')
    !call writeIRDef (3,ctrl%mpireg,'MMOD')
    !call writeIRDef (4,ctrl%mpireg,'IPVT1')
    !call writeIRDef (5,ctrl%mpireg,'IPVT2')
    !call writeIRDef (6,ctrl%mpireg,'MPTOP')
    !call writeIRDef (7,ctrl%mpireg,'MPRPAR')
    !call writeIRDef (8,ctrl%mpireg,'MMTOP')
    call writeRRDef (1,ctrl%mprreg,'RCNT')
    call writeRRDef (2,ctrl%mprreg,'VPREG')
    !call writeRRDef (3,ctrl%mprreg,'VBREG')
    !call writeRRDef (4,ctrl%mprreg,'VLREG')
    !call writeRRDef (5,ctrl%mprreg,'VDREG')
    call writeRRDef (6,ctrl%mprreg,'VHREG')
    !call writeRRDef (7,ctrl%mprreg,'TOL')
    !call writeRRDef (8,ctrl%mprreg,'ERR')
    !call writeRRDef (9,ctrl%mprreg,'AJAC')
    !call writeRRDef(10,ctrl%mprreg,'ABE')
    !call writeRRDef(11,ctrl%mprreg,'ALO')
    !call writeRRDef(12,ctrl%mprreg,'RPAR')

    if (associated(ctrl%delay) .and. ctrl%ireg(10) > 0) then
       nV = ctrl%ireg(10) ! Initial delay buffer size
       call writeRRDef (13,name='DELAY')
    end if

    write(idatd,"('}')")

  contains

    !> @brief Writes results database header for an indexed integer array.
    subroutine writeIRDef (i,mp,name)
      integer         , intent(in) :: i, mp(:)
      character(len=*), intent(in) :: name
      nV = mp(i+1) - mp(i)
      rdb%nBytes = rdb%nBytes + nV*nbi_p
      call writeVarDef (rdb,iIVar(i),idatd,name,'NONE', &
           &            'VECTOR;('//trim(adjustl(StrId(nV)))//')',nBits=1)
    end subroutine writeIRDef

    !> @brief Writes results database header for an indexed real array.
    subroutine writeRRDef (i,mp,name)
      integer         , intent(in) :: i
      integer,optional, intent(in) :: mp(:)
      character(len=*), intent(in) :: name
      if (present(mp)) nV = mp(i+1) - mp(i)
      rdb%nBytes = rdb%nBytes + nV*nbd_p
      call writeVarDef (rdb,iRVar(i),idatd,name,'NONE', &
           &            'VECTOR;('//trim(adjustl(StrId(nV)))//')',nBits=8*nbd_p)
    end subroutine writeRRDef

  end subroutine writeCtrlSysHeader


  !!============================================================================
  !> @brief Writes control system data to the results database.
  !>
  !> @param rdb Results database file for control system data
  !> @param[in] ctrl Control system data
  !> @param[in] nStep Time step number
  !> @param[in] time Current time
  !> @param[out] ierr Error flag
  !>
  !> @details Only the arrays that need to be restored in a restart are saved.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 25 Oct 2008

  subroutine writeCtrlSysDB (rdb,ctrl,nStep,time,ierr)

    use KindModule       , only : i8
    use RDBModule        , only : RDBType, writeTimeStepDB, writeRDB
    use RDBModule        , only : checkStepSize
    use reportErrorModule, only : reportError, debugFileOnly_p, warning_p

    type(RDBType)    , intent(inout) :: rdb
    type(ControlType), intent(in)    :: ctrl
    integer(i8)      , intent(in)    :: nStep
    real(dp)         , intent(in)    :: time
    integer          , intent(out)   :: ierr

    !! Local variables
    integer :: nBytesPerStep = 0
    integer :: nWarning = 0

    !! --- Logic section ---

    if (nBytesPerStep == 0) nBytesPerStep = -int(rdb%nBytes)

    call writeTimeStepDB (rdb,nStep,time,ierr)

    call writeIREG (1,ctrl%mpireg)
    !call writeIREG (2,ctrl%mpireg)
    !call writeIREG (3,ctrl%mpireg)
    !call writeIREG (4,ctrl%mpireg)
    !call writeIREG (5,ctrl%mpireg)
    !call writeIREG (6,ctrl%mpireg)
    !call writeIREG (7,ctrl%mpireg)
    !call writeIREG (8,ctrl%mpireg)
    call writeRREG (1,ctrl%mprreg)
    call writeRREG (2,ctrl%mprreg)
    !call writeRREG (3,ctrl%mprreg)
    !call writeRREG (4,ctrl%mprreg)
    !call writeRREG (5,ctrl%mprreg)
    call writeRREG (6,ctrl%mprreg)
    !call writeRREG (7,ctrl%mprreg)
    !call writeRREG (8,ctrl%mprreg)
    !call writeRREG (9,ctrl%mprreg)
    !call writeRREG(10,ctrl%mprreg)
    !call writeRREG(11,ctrl%mprreg)
    !call writeRREG(12,ctrl%mprreg)

    if (associated(ctrl%delay) .and. ctrl%ireg(10) > 0) then
       !! We can only save initial portion of the delay buffer
       call writeRDB (rdb,ctrl%delay(1:ctrl%ireg(10)),ierr,.true.)
       if (size(ctrl%delay) > ctrl%ireg(10) .and. nWarning < 10) then
          call reportError (warning_p,'DELAY buffer is saved only with its '// &
               &            'initial size.','A restart of this simulation '// &
               &            'beyond this time might therefore not work.')
          nWarning = nWarning + 1
       end if
    end if

    if (nBytesPerStep < 1 .and. ierr >= 0) then
       !! Check consistency between header and the actual data
       nBytesPerStep = nBytesPerStep + int(rdb%nBytes)
       call checkStepSize (rdb,nBytesPerStep,ierr)
    end if

    if (ierr < 0) call reportError (debugFileOnly_p,'writeCtrlSysDB')

  contains

    !> @brief Writes an indexed integer array to the results database.
    subroutine writeIREG (i,mp)
      integer, intent(in) :: i, mp(:)
      call writeRDB (rdb,ctrl%ireg(mp(i):mp(i+1)-1),ierr)
    end subroutine writeIREG

    !> @brief Writes an indexed real array to the results database.
    subroutine writeRREG (i,mp)
      integer, intent(in) :: i, mp(:)
      call writeRDB (rdb,ctrl%rreg(mp(i):mp(i+1)-1),ierr,.true.)
    end subroutine writeRREG

  end subroutine writeCtrlSysDB

end module ControlTypeModule
