!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file massTypeModule.f90
!> @brief Point mass object data containers.

!!==============================================================================
!> @brief Module with data types representing point mass objects.
!>
!> @details The module also contains subroutines for accessing the mass data.

module MassTypeModule

  use TriadTypeModule   , only : TriadType, IdType, dp
  use FunctionTypeModule, only : EngineType

  implicit none

  !> @brief Data type representing a point mass object.
  type MassType

     type(IdType) :: id !< General identification data

     integer :: samElNum  !< Element number for sam reference (mpmnpc)
     logical :: addedMass !< if .true., this is an added mass (i.e., no gravity)

     type(TriadType) , pointer :: triad      !< Triad associated with this mass
     type(EngineType), pointer :: massEngine !< Scaling function for mass
     type(EngineType), pointer :: IIEngine   !< Scaling function for intertia

     real(dp) :: m0(3)    !< Constant mass values
     real(dp) :: m1(3)    !< Scalable mass values
     real(dp) :: II0(3,3) !< Constant inertia tensor
     real(dp) :: II1(3,3) !< Scalable inertia tensor

     real(dp) :: mass(3) !< Current mass    = m0  + m1*massEngine
     real(dp) :: II(3,3) !< Current Inertia = II0 + II1*IIEngine

     real(dp) :: Epot0 !< Initial potential energy
     real(dp) :: Epot  !< Potential energy relative to Epot0
     real(dp) :: Ekin  !< Kinetic energy

  end type MassType


  !> @brief Standard routine for writing an object to file.
  interface WriteObject
     module procedure WriteMassType
  end interface

  private :: WriteMassType


contains

  !!============================================================================
  !> @brief Initializes the mass objects with data from the solver input file.
  !>
  !> @param[in] infp File unit number for the solver input file
  !> @param[out] masses Array of all point mass objects in the model
  !> @param[in] triads Array of all triads in the model
  !> @param[in] engines Array of all engines in the model
  !> @param[out] err Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 25 Oct 1999

  subroutine InitiateMasses (infp,masses,triads,engines,err)

    use IdTypeModule          , only : ldesc_p, initId, getId, ReportInputError
    use TriadTypeModule       , only : allocateNodeForce, GetPtrToId
    use FunctionTypeModule    , only : GetPtrToId
    use inputUtilities        , only : iuGetNumberOfEntries, iuSetPosAtNextEntry
    use progressModule        , only : lterm
    use reportErrorModule     , only : allocationError
    use reportErrorModule     , only : reportError, warning_p, debugFileOnly_p
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_isTrue

    integer         , intent(in)  :: infp
    type(MassType)  , pointer     :: masses(:)
    type(TriadType) , intent(in)  :: triads(:)
    type(EngineType), intent(in)  :: engines(:)
    integer         , intent(out) :: err

    !! Local variables
    integer :: idIn, nMass, nSkip, stat

    !! Define the MASS namelist
    integer            :: id, extId(10), triadId, massEngineId, IIEngineId, dof
    logical            :: addedMass
    character(ldesc_p) :: extDescr
    real(dp)           :: mass0(3), mass1(3), II0(3,3), II1(3,3)
    namelist /MASS/ id, extId, extDescr, triadId, massEngineId, IIEngineId, &
         &          dof, addedMass, mass0, mass1, II0, II1

    !! --- Logic section ---

    nSkip = 0
    nMass = iuGetNumberOfEntries(infp,'&MASS',err)
    if (err /= 0) return
    write(lterm,*) 'Number of &MASS =',nMass

    if (ffa_cmdlinearg_isTrue('noAddedMass')) then ! Ignore all added masses
       do idIn = 1, nMass
          if (iuSetPosAtNextEntry(infp,'&MASS')) then
             addedMass=.false.
             read(infp,nml=MASS,iostat=stat)
             if (stat == 0 .and. addedMass) nSkip = nSkip + 1
          end if
       end do
       if (nSkip == nMass) return
       nMass = nMass - nSkip
       rewind(infp)
    end if

    call DeallocateMasses (masses)
    allocate(masses(nMass),STAT=stat)
    if (stat /= 0) then
       err = AllocationError('InitiateMasses')
       return
    end if

    do idIn = 1, nMass

       call NullifyMass (masses(idIn))
100    continue
       if (.not. iuSetPosAtNextEntry(infp,'&MASS')) then
          err = err - 1
          call ReportInputError ('MASS',idIn)
          cycle
       end if

       !! Default values
       id=0; extId=0; extDescr=''; triadId=0; massEngineId=0; IIEngineId=0
       dof=0; mass0=0.0_dp; mass1=0.0_dp; II0=0.0_dp; II1=0.0_dp
       addedMass=.false.

       read(infp,nml=MASS,iostat=stat)
       if (stat /= 0) then
          err = err - 1
          call ReportInputError ('MASS',idIn)
          cycle
       end if
       if (addedMass .and. nSkip > 0) goto 100 ! Skip all added masses

       call initId (masses(idIn)%id,id,extId,extDescr,stat)
       masses(idIn)%addedMass = addedMass

       if (dof < 0 .or. dof > 3) then
          !! Masses are given in separate directions
          masses(idIn)%m0 = mass0
          masses(idIn)%m1 = mass1
       else if (dof == 0) then
          !! The mass is equal in all directions (this is the default)
          masses(idIn)%m0 = mass0(1)
          masses(idIn)%m1 = mass1(1)
       else
          !! The mass is effective only in one direction
          masses(idIn)%m0(dof) = mass0(1)
          masses(idIn)%m1(dof) = mass1(1)
       end if

       masses(idIn)%II0 = II0
       masses(idIn)%II1 = II1

       !! Connect the mass engine
       if (massEngineId > 0) then
          masses(idIn)%massEngine => GetPtrToId(engines,massEngineId)
          if (.not. associated(masses(idIn)%massEngine)) then
             err = err - 1
             call ReportInputError ('MASS',idIn,masses(idIn)%id)
          end if
       end if

       !! Connect the inertia engine
       if (IIEngineId > 0) then
          masses(idIn)%IIEngine => GetPtrToId(engines,IIEngineId)
          if (.not. associated(masses(idIn)%IIEngine)) then
             err = err - 1
             call ReportInputError ('MASS',idIn,masses(idIn)%id)
          end if
       end if

       !! Connect to the triad
       masses(idIn)%triad => GetPtrToId(triads,triadId)
       if (.not. associated(masses(idIn)%triad)) then
          err = err - 1
          call ReportInputError ('MASS',idIn,masses(idIn)%id)
       else if (masses(idIn)%triad%nDOFs < 1) then
          call reportError (warning_p,'Mass on earth-linked triads ignored', &
               &            'Mass Id :'//getId(masses(idIn)%id), &
               &            'Triad Id :'//getId(masses(idIn)%triad%id))
       else
          !! If this triad is also connected to a superelement,
          !! allocate a separate array for storage of the total nodal force
          stat = err
          call allocateNodeForce (masses(idIn)%triad,err)
          if (err < stat) exit
       end if
    end do

    if (err < 0) call reportError (debugFileOnly_p,'InitiateMasses')

  end subroutine InitiateMasses


  !!============================================================================
  !> @brief Standard routine for writing an object to file.
  !>
  !> @param[in] mass The masstypemodule::masstype object to write
  !> @param[in] io File unit number to write to
  !> @param[in] complexity If present, the value indicates the amount of print
  !>
  !> @author Karl Erik Thoresen
  !>
  !> @date 27 Sep 1998

  subroutine WriteMassType (mass,io,complexity)

    use IdTypeModule, only : writeId

    type(MassType)  , intent(in) :: mass
    integer         , intent(in) :: io
    integer,optional, intent(in) :: complexity

    !! Local variables
    integer :: i

    !! --- Logic section ---

    write(io,'(A)') 'Mass','{'
    call writeId (mass%id,io)

    write(io,*) 'samElNum    =', mass%samElNum
    write(io,*) 'addedMass?  =', mass%addedMass

    if (associated(mass%triad)) then
       write(io,*) 'triadId     =', mass%triad%id%baseId
    else
       write(io,*) 'triad       = NULL'
    end if

    if (associated(mass%massEngine)) then
       write(io,*) 'massEngine  =', mass%massEngine%id%baseId
    end if

    if (associated(mass%IIEngine)) then
       write(io,*) 'IIEngine    =', mass%IIEngine%id%baseId
    end if

    write(io,*) 'm0          =', mass%m0
    write(io,*) 'm1          =', mass%m1
    write(io,6) 'II0         =',(mass%II0(i,:),i=1,3)
    write(io,6) 'II1         =',(mass%II1(i,:),i=1,3)

    if (present(complexity)) then
       if (complexity >= 2) then
          write(io,*) 'mass        =', mass%mass
          write(io,6) 'II          =',(mass%II(i,:),i=1,3)
          write(io,*) 'Epot0       =', mass%Epot0
          write(io,*) 'Epot        =', mass%Epot
          write(io,*) 'Ekin        =', mass%Ekin
       end if
    end if

    write(io,'(A)') '}'

6   format(1X,A13,3ES15.6E3/(14X,3ES15.6E3))

  end subroutine WriteMassType


  !!============================================================================
  !> @brief Initializes a point mass object.
  !>
  !> @param mass The masstypemodule::masstype object to initialize
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 30 May 2002

  subroutine NullifyMass (mass)

    use IdTypeModule, only : nullifyId

    type(MassType), intent(out) :: mass

    !! --- Logic section ---

    call nullifyId (mass%id)

    mass%samElNum = 0
    mass%addedMass = .false.

    nullify(mass%triad)
    nullify(mass%massEngine)
    nullify(mass%IIEngine)

    mass%m0    = 0.0_dp
    mass%m1    = 0.0_dp
    mass%mass  = 0.0_dp
    mass%II0   = 0.0_dp
    mass%II1   = 0.0_dp
    mass%II    = 0.0_dp

    mass%ePot0 = 0.0_dp
    mass%ePot  = 0.0_dp
    mass%eKin  = 0.0_dp

  end subroutine NullifyMass


  !!============================================================================
  !> @brief Deallocates the point mass objects.
  !>
  !> @param masses The masstypemodule::masstype objects to deallocate
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Jan 2017

  subroutine DeallocateMasses (masses)

    use IdTypeModule, only : deallocateId

    type(MassType), pointer :: masses(:)

    !! Local variables
    integer :: i

    !! --- Logic section ---

    do i = 1, size(masses)
       call deallocateId (masses(i)%id)
    end do
    deallocate(masses)
    nullify(masses)

  end subroutine DeallocateMasses

end module MassTypeModule
