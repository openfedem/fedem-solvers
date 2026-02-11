!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file forceTypeModule.f90
!> @brief External point load object data container.

!!==============================================================================
!> @brief Module with data types representing external point load objects.
!>
!> @details The module also contains subroutines for accessing or manipulating
!> the load data, as well as reading loads from the solver input file.

module ForceTypeModule

  use KindModule                , only : dp, lfnam_p
  use TriadTypeModule           , only : TriadType, IdType
  use MasterSlaveJointTypeModule, only : MasterSlaveJointType
  use SupElTypeModule           , only : SupElType
  use FunctionTypeModule        , only : EngineType

  implicit none

  integer, parameter :: forceVec_p  = -1 !< Force vector type id
  integer, parameter :: momentVec_p = -2 !< Moment vector type id

  !> Force DOF type names
  character(len=9), parameter :: forceType_p(-2:6) = (/ 'moment   ', &
       &                                                'force    ', &
       &                                                'UNDEFINED', &
       &                                                'X-force  ', &
       &                                                'Y-force  ', &
       &                                                'Z-force  ', &
       &                                                'X-moment ', &
       &                                                'Y-moment ', &
       &                                                'Z-moment ' /)

  !> @brief Data type representing an external point load object.
  type ForceType

     type(IdType) :: id !< General identification data

     !> @brief Local DOF the force is acting on [1,6].
     !> @details
     !> - = -1 : Force vector in the direction of @a (v1-v2)
     !> - = -2 : Moment vector in the direction of @a (v1-v2)
     integer :: dof

     !> @brief Flag telling when the force value is to be updated.
     !> @details
     !> - = 1 : Update only in the beginning of each step
     !> - = 2 : Update during Newton iterations
     !> - &gt; 9 : As @a updateFlag-10, but applied in local axes
     integer :: updateFlag

     !> @brief Load type flag.
     !> @details
     !> - = 0 : Normal, time-domain load (or constant)
     !> - = 1 : Frequency domain load function
     !> - = 2 : Frequency domain displacement function
     !> - = 3 : Frequency domain velocity function
     !> - = 4 : Frequency domain acceleration function
     integer :: loadType
     character(len=lfnam_p) :: fftFile !< Name of file with FFT input

     type(TriadType)           , pointer :: triad !< Triad the force acts on
     type(MasterSlaveJointType), pointer :: joint !< Joint the force acts in

     type(EngineType), pointer :: engine !< Function giving current force value
     real(dp)                  :: f0     !< Constant force value
     real(dp)                  :: f1     !< Scaling factor for function value

     !! v1 and v2 are points in the local coordinate system of sup1 and sup2,
     !! defining the start and end point of the force direction vector.
     !! Fixed directions are handled by sup1 and sup2 being NULL.
     real(dp) :: v1(3) !< Start point of the force direction vector
     real(dp) :: v2(3) !< End point of the force direction vector
     type(SupElType), pointer :: sup1 !< Local coordinate system for v1
     type(SupElType), pointer :: sup2 !< Local coordinate system for v2
     type(TriadType), pointer :: trp1 !< Start point of force direction vector
     type(TriadType), pointer :: trp2 !< End point of force direction vector

     real(dp) :: F               !< Current force value
     real(dp) :: Fprev           !< Previous force value
     real(dp) :: forceDir(3)     !< Force direction vector
     real(dp) :: forceDirPrev(3) !< Previous force direction
     real(dp) :: eInp            !< Accumulated input energy from this force

     logical  :: saveVar(3) !< Flags indicating which variables should be saved

  end type ForceType


  !> @brief Standard routine for writing an object to file.
  interface WriteObject
     module procedure WriteForceType
  end interface

  private :: WriteForceType


contains

  !!============================================================================
  !> @brief Initializes the force type with data from the solver input file.
  !>
  !> @param[in] infp File unit number for the solver input file
  !> @param[in] triads Array of all triads in the model
  !> @param[in] joints Array of all joints in the model
  !> @param[in] sups Array of all superelements in the model
  !> @param[in] engines Array of all general functions in the model
  !> @param[out] forces Array of all external point loads in the model
  !> @param[out] err Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Karl Erik Thoresen
  !>
  !> @date 25 Sep 1998

  subroutine InitiateForces (infp,triads,joints,sups,engines,forces,err)

    use IdTypeModule              , only : ldesc_p, initId, ReportInputError
    use TriadTypeModule           , only : GetPtrToId, allocateNodeForce
    use SupElTypeModule           , only : GetPtrToId
    use MasterSlaveJointTypeModule, only : GetPtrToId
    use FunctionTypeModule        , only : GetPtrToId
    use inputUtilities            , only : iuGetNumberOfEntries
    use inputUtilities            , only : iuSetPosAtNextEntry, iuCharToInt
    use progressModule            , only : lterm
    use reportErrorModule         , only : allocationError
    use reportErrorModule         , only : reportError, debugFileOnly_p
    use FFaFilePathInterface      , only : ffa_checkPath
    use FFaCmdLineArgInterface    , only : ffa_cmdlinearg_getbool

    integer                   , intent(in)  :: infp
    type(TriadType)           , intent(in)  :: triads(:)
    type(MasterSlaveJointType), intent(in)  :: joints(:)
    type(SupElType)           , intent(in)  :: sups(:)
    type(EngineType)          , intent(in)  :: engines(:)
    type(ForceType)           , pointer     :: forces(:)
    integer                   , intent(out) :: err

    !! Local variables
    integer :: i, idIn, nLoads, stat
    logical :: doFRA, allSec, allForce, allLoad, allEnergy

    !! Define the LOAD namelist
    integer            :: id, extId(10), lDof, updateFlag, loadType
    integer            :: triadId, jointId, loadEngineId
    integer            :: supEl1Id, supEl2Id, triad1Id, triad2Id
    character(ldesc_p) :: extDescr
    character(lfnam_p) :: fileName
    character(len=20)  :: type
    real(dp)           :: vec1(3), vec2(3), f0, f1
    integer            :: saveVar(3)
    namelist /LOAD/ id, extId, extDescr, type, lDof, updateFlag, loadType, &
         &          triadId, jointId, loadEngineId, supEl1Id, supEl2Id, &
         &          triad1Id, triad2Id, vec1, vec2, f0, f1, saveVar, fileName

    !! --- Logic section ---

    nLoads = iuGetNumberOfEntries(infp,'&LOAD',err)
    if (err /= 0) return
    write(lterm,*) 'Number of &LOAD =',nLoads

    call DeallocateForces (forces)
    allocate(forces(nLoads),STAT=stat)
    if (stat /= 0)then
       stat = AllocationError('InitiateForces 10')
       return
    end if

    call ffa_cmdlinearg_getbool ('frequency_domain',doFRA)
    call ffa_cmdlinearg_getbool ('allSecondaryVars',allSec)
    call ffa_cmdlinearg_getbool ('allLoadVars',allLoad)
    call ffa_cmdlinearg_getbool ('allForceVars',allForce)
    call ffa_cmdlinearg_getbool ('allEnergyVars',allEnergy)

    do idIn = 1, nLoads

       call NullifyForce (forces(idIn))
       if (.not. iuSetPosAtNextEntry(infp,'&LOAD')) then
          err = err - 1
          call ReportInputError ('LOAD',idIn)
          cycle
       end if

       !! Default values
       id=0; extId=0; extDescr=''
       type=''; lDof=0; updateFlag=2; loadType=0; fileName=''
       triadId=0; jointId=0; loadEngineId=0
       supEl1Id=0; supEl2Id=0; triad1Id=0; triad2Id=0
       vec1=0.0_dp; vec2=0.0_dp; f0=0.0_dp; f1=0.0_dp; saveVar=0

       read(infp,nml=LOAD,iostat=stat)
       if (stat /= 0) then
          err = err - 1
          call ReportInputError ('LOAD',idIn)
          cycle
       end if

       call initId (forces(idIn)%id,id,extId,extDescr,stat)

       if (triadId > 0) then ! Force is acting on a triad

          !! Set the pointer to the triad the load is acting on
          forces(idIn)%triad => GetPtrToId(triads,triadId)
          if (.not. associated(forces(idIn)%triad)) then
             err = err - 1
             call ReportInputError ('LOAD',idIn,forces(idIn)%id, &
                  &                 'Invalid triad ID')
             cycle
          end if

          if (lDof > 0 .and. lDof < 7) then
             forces(idIn)%dof = lDof
          else
             forces(idIn)%dof = iuCharToInt(type,forceType_p,startIndex=-2)
          end if
          if (forces(idIn)%dof < -2 .or. forces(idIn)%dof > 6) then
             err = err - 1
             call ReportInputError ('LOAD',idIn,forces(idIn)%id, &
                  &                 'Invalid local DOF in referred triad')
             cycle
          end if

          !! If this triad also is connected to a superelement,
          !! allocate a separate array for storage of the total nodal force
          stat = err
          call allocateNodeForce (forces(idIn)%triad,err)
          if (err < stat) exit

          !! Set the pointer to triad or superelement for the FROM point
          if (supEl1Id > 0) then
             forces(idIn)%sup1 => GetPtrToId(sups,supEl1Id)
             if (.not. associated(forces(idIn)%sup1)) then
                err = err - 1
                call ReportInputError ('LOAD',idIn,forces(idIn)%id, &
                     &                 'Invalid superelement ID for FROM-point')
                cycle
             end if
          else if (triad1Id > 0) then
             forces(idIn)%trp1 => GetPtrToId(triads,triad1Id)
             if (.not. associated(forces(idIn)%trp1)) then
                err = err - 1
                call ReportInputError ('LOAD',idIn,forces(idIn)%id, &
                     &                 'Invalid triad ID for FROM-point')
                cycle
             end if
          end if

          !! Set the pointer to triad or superelement for the TO point
          if (supEl2Id > 0) then
             forces(idIn)%sup2 => GetPtrToId(sups,supEl2Id)
             if (.not. associated(forces(idIn)%sup2)) then
                err = err - 1
                call ReportInputError ('LOAD',idIn,forces(idIn)%id, &
                     &                 'Invalid superelement ID for TO-point')
                cycle
             end if
          else if (triad2Id > 0) then
             forces(idIn)%trp2 => GetPtrToId(triads,triad2Id)
             if (.not. associated(forces(idIn)%trp2)) then
                err = err - 1
                call ReportInputError ('LOAD',idIn,forces(idIn)%id, &
                     &                 'Invalid triad ID for TO-point')
                cycle
             end if
          end if

       else if (jointId > 0) then ! Force is acting on a joint dof

          !! Set the pointer to the joint the load is acting on
          forces(idIn)%joint => GetPtrToId(joints,jointId)
          if (.not. associated(forces(idIn)%joint)) then
             err = err - 1
             call ReportInputError ('LOAD',idIn,forces(idIn)%id, &
                  &                 'Invalid joint ID')
             cycle
          end if

          do i = 1, forces(idIn)%joint%nJointDofs
             if (forces(idIn)%joint%jointDofs(i)%lDof == lDof) then
                forces(idIn)%dof = i
                forces(idIn)%joint%jointDofs(i)%loadIdx = idIn
                exit
             end if
          end do
          if (forces(idIn)%dof == 0) then
             err = err - 1
             call ReportInputError ('LOAD',idIn,forces(idIn)%id, &
                  &                 'Invalid local DOF in referred joint')
             cycle
          end if

       else
          err = err - 1
          call ReportInputError ('LOAD',idIn,forces(idIn)%id)
          cycle
       end if

       !! Set the load engine pointer
       if (loadEngineId > 0) then
          forces(idIn)%engine => GetPtrToId(engines,loadEngineId)
          if (.not. associated(forces(idIn)%engine)) then
             err = err - 1
             call ReportInputError ('LOAD',idIn,forces(idIn)%id, &
                  &                 'Invalid engine ID')
             cycle
          end if
       end if

       forces(idIn)%updateFlag = updateFlag
       if (doFRA .or. loadType > 1) then
          forces(idIn)%loadType = loadType
          if (fileName /= '') then
             call ffa_checkPath (fileName)
             forces(idIn)%fftFile = fileName
          end if
       end if

       forces(idIn)%f0 = f0
       forces(idIn)%f1 = f1

       forces(idIn)%v1 = vec1
       forces(idIn)%v2 = vec2

       !! Determine which secondary triad variables will be saved. Settings in
       !! the solver input file are overruled by command-line arguments, if any.
       do i = 1, size(saveVar)
          forces(idIn)%saveVar(i) = saveVar(i) > 0 .or. allSec .or. allLoad
       end do
       if (allForce)  forces(idIn)%saveVar(1) = .true. ! Global force (vector)
       if (allForce)  forces(idIn)%saveVar(2) = .true. ! Local force (value)
       if (allEnergy) forces(idIn)%saveVar(3) = .true. ! Energy quantities

    end do

    if (err < 0) call reportError (debugFileOnly_p,'InitiateForces')

  end subroutine InitiateForces


  !!============================================================================
  !> @brief Standard routine for writing an object to file.
  !>
  !> @param[in] force The forcetypemodule::forcetype object to write
  !> @param[in] io File unit number to write to
  !> @param[in] complexity If present, the value indicates the amount of print
  !>
  !> @author Karl Erik Thoresen
  !>
  !> @date 27 Sep 1998

  subroutine WriteForceType (force,io,complexity)

    use IdTypeModule, only : writeId

    type(ForceType) , intent(in) :: force
    integer         , intent(in) :: io
    integer,optional, intent(in) :: complexity

    !! --- Logic section ---

    write(io,'(A)') 'Force','{'
    call writeId (force%id,io)

    if (associated(force%triad)) then
       write(io,*) 'dof/type    =', forceType_p(force%dof)
       write(io,*) 'triad%id    =', force%triad%id%baseId,force%triad%id%userId
    else if (associated(force%joint)) then
       write(io,*) 'dof         =', force%dof
       write(io,*) 'joint%id    =', force%joint%id%baseId,force%joint%id%userId
    end if

    if (associated(force%engine)) then
       write(io,*) 'engine%id   =',force%engine%id%baseId,force%engine%id%userId
    end if
    write(io,*) 'f0, f1      =', force%f0, force%f1

    if (associated(force%sup1)) then
       write(io,*) 'sup1(id)    =', force%sup1%id%baseId,force%sup1%id%userId
    else
       write(io,*) 'sup1        = NULL'
    end if
    if (associated(force%sup2)) then
       write(io,*) 'sup2(id)    =', force%sup2%id%baseId,force%sup2%id%userId
    else
       write(io,*) 'sup2        = NULL'
    end if
    if (associated(force%trp1)) then
       write(io,*) 'trp1(id)    =', force%trp1%id%baseId,force%trp1%id%userId
    else
       write(io,*) 'trp1        = NULL'
    end if
    if (associated(force%trp2)) then
       write(io,*) 'trp2(id)    =', force%trp2%id%baseId,force%trp2%id%userId
    else
       write(io,*) 'trp2        = NULL'
    end if

    write(io,*) 'updateFlag  =', force%updateFlag
    write(io,*) 'loadType    =', force%loadType
    if (force%fftFile /= '') then
       write(io,*) 'fftFile     =', trim(force%fftFile)
    end if
    write(io,*) 'saveVar     =', force%saveVar

    if (present(complexity)) then
       if (complexity >= 2) then
          write(io,*) 'v1          =', force%v1
          write(io,*) 'v2          =', force%v2
          write(io,*) 'F           =', force%F
          write(io,*) 'eInp        =', force%eInp
       end if
    end if

    write(io,'(A)') '}'

  end subroutine WriteForceType


  !!============================================================================
  !> @brief Initializes an external point load object.
  !>
  !> @param force The forcetypemodule::forcetype object to initialize
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 14 Mar 2001

  subroutine NullifyForce (force)

    use IdTypeModule, only : nullifyId

    type(ForceType), intent(out) :: force

    !! --- Logic section ---

    call nullifyId (force%id)
    force%dof = 0
    force%updateFlag = 0
    force%loadType = 0
    force%fftFile = ''

    nullify(force%triad)
    nullify(force%joint)
    nullify(force%engine)

    force%f0 = 0.0_dp
    force%f1 = 0.0_dp
    force%v1 = 0.0_dp
    force%v2 = 0.0_dp

    nullify(force%sup1)
    nullify(force%sup2)
    nullify(force%trp1)
    nullify(force%trp2)

    force%F = 0.0_dp
    force%Fprev = 0.0_dp
    force%forceDir = 0.0_dp
    force%forceDirPrev = 0.0_dp
    force%eInp = 0.0_dp
    force%saveVar = .false.

  end subroutine NullifyForce


  !!============================================================================
  !> @brief Deallocates the external point load objects.
  !>
  !> @param forces The forcetypemodule::forcetype objects to deallocate
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Jan 2017

  subroutine DeallocateForces (forces)

    use IdTypeModule, only : deallocateId

    type(ForceType), pointer :: forces(:)

    !! Local variables
    integer :: i

    !! --- Logic section ---

    do i = 1, size(forces)
       call deallocateId (forces(i)%id)
    end do

    deallocate(forces)
    nullify(forces)

  end subroutine DeallocateForces

end module ForceTypeModule
