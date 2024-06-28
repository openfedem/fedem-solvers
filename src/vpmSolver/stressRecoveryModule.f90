!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file stressRecoveryModule.f90
!>
!> @brief Subroutines for stress and/or strain gage recovery during solving.

!!==============================================================================
!> @brief Module with subroutines for recovery within the time integration loop.
!>
!> @details This module contains data and subroutines for performing
!> stress- and/or strain gage recovery during the time integration loop,
!> by calling subroutines from the stress recovery solver.

module StressRecoveryModule

  use SamModule          , only : SamType, dp
  use RDBModule          , only : RDBType
  use StrainRosetteModule, only : StrainElementType
  use DisplacementModule , only : rk

  implicit none

  private

  !> @brief Data type holding recovery data for superelements.
  type RecPartType
     integer           :: baseId !< Base ID of superelement to recover
     integer           :: nRos   !< Number of strain rosettes in superelement
     integer           :: ndim   !< Number of external DOFs in superelement
     type(SamType)     :: sam    !< Assembly management data for superelement
     type(RDBType)     :: rdb    !< Stress results database file
     type(RDBType)     :: rdb2   !< Strain gage results database file
     integer , pointer :: gno(:) !< Group node numbers
     integer , pointer :: map(:) !< Mapping from group DOF to global DOF
     real(rk), pointer :: B(:,:) !< Displacement recovery matrix for node group
     real(rk), pointer :: sv(:)  !< Recovered internal nodal displacements
     real(dp), pointer :: vms(:) !< Recovered von Mises stresses
  end type RecPartType

  !> @brief Superelement recovery data container
  type(RecPartType),   allocatable, save :: part(:)
  !> @brief Strain rosette container
  type(StrainElementType), pointer, save :: strainRosettes(:) => null()

  public :: getStrainRosette, readStrainRosettes, initRecovery, closeRecovery
  public :: writeRosettes2Ftn, writeRecoveryHeaders, flushRecoveryFiles
  public :: getStrainGagesSize, initGageStrains
  public :: stressRecovery, gageRecovery, saveGageResults
  public :: getDeformation, getDeformationSize
  public :: getStress, getStressSize
  public :: getGageRecoveryFiles


contains

  !!============================================================================
  !> @brief Returns pointer to (first) strain rosette object with specified ID.
  !>
  !> @param[in] id Base ID of the object to search for
  !>
  !> @details If the strain rosette object is not found, NULL is returned.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 7 Sep 2016

  function getStrainRosette (id) result(ptr)

    use StrainRosetteModule, only : StrainRosetteType

    type(StrainRosetteType), pointer    :: ptr
    integer                , intent(in) :: id

    !! Local variables
    integer :: i

    !! --- Logic section ---

    if (associated(strainRosettes)) then
       do i = 1, size(strainRosettes)
          if (strainRosettes(i)%id%baseId == id) then
             ptr => strainRosettes(i)%data(1)
             return
          end if
       end do
    end if

    nullify(ptr)
    write(*,*) '*** getStrainRosette returned nullified, BaseId =',id

  end function getStrainRosette


  !!============================================================================
  !> @brief Returns deformation vector of the specified part.
  !>
  !> @param[in] id Base ID of the FE part to return deformations for
  !> @param[in] data Array with deformation vector
  !> @param[in] ndat Size of @a data array
  !> @param idat Running array index, negative on return if error occurred
  !>
  !> @callergraph
  !>
  !> @author Runar Heggelien Refsnaes
  !>
  !> @date 30 Oct 2017

  subroutine getDeformation (id,data,ndat,idat)

    integer , intent(in)    :: id, ndat
    integer , intent(inout) :: idat
    real(dp), intent(out)   :: data(ndat)

    !! Local variables
    integer :: i, j, k

    !! --- Logic section ---

    if (.not.allocated(part)) return ! Do nothing if invoked after deallocation

    do i = 1, size(part)
       if (part(i)%baseId == id) then
          if (idat + 3*part(i)%sam%nnod - 1 > ndat) then
             idat = -i
          else
             do j = 1, part(i)%sam%nnod
                if (part(i)%sam%minex(j) > 0) then
                   do k = 0, 2
                      data(idat+k) = real(part(i)%sv(part(i)%sam%madof(j)+k),dp)
                   end do
                else
                   data(idat:idat+2) = 0.0_dp
                end if
                idat = idat + 3
             end do
          end if
          return
       end if
    end do

  end subroutine getDeformation


  !!============================================================================
  !> @brief Returns the size of deformation vector for the specified part.
  !>
  !> @param[in] id Base ID of the FE part to return deformations for
  !>
  !> @callergraph
  !>
  !> @author Runar Heggelien Refsnaes
  !>
  !> @date 30 Oct 2017

  function getDeformationSize (id) result(ndat)

    integer, intent(in) :: id

    !! Local variables
    integer :: i, ndat

    !! --- Logic section ---

    if (.not.allocated(part)) then
       !! Invoked before allocation or after deallocation
       ndat = -999
       return
    end if

    do i = 1, size(part)
       if (part(i)%baseId == id) then
          ndat = 3*part(i)%sam%nnod
          return
       end if
    end do

    ndat = -1

  end function getDeformationSize


  !!============================================================================
  !> @brief Returns the array of von Mises stresses for the specified part.
  !>
  !> @param[in] id Base ID of the FE part to return stresses for
  !> @param[in] data Array with von Mises stresses
  !> @param[in] ndat Size of @a data array
  !> @param idat Running array index, negative on return if error occurred
  !>
  !> @callergraph
  !>
  !> @author Runar Heggelien Refsnaes
  !>
  !> @date 31 Oct 2017

  subroutine getStress (id,data,ndat,idat)

    integer , intent(in)    :: id, ndat
    integer , intent(inout) :: idat
    real(dp), intent(out)   :: data(ndat)

    !! Local variables
    integer :: i, nvms

    !! --- Logic section ---

    if (.not.allocated(part)) return ! Do nothing if invoked after deallocation

    do i = 1, size(part)
       if (part(i)%baseId == id .and. associated(part(i)%vms)) then
          nvms = size(part(i)%vms)
          if (idat + nvms - 1 > ndat) then
             idat = -i
          else if (nvms > 0) then
             call DCOPY (nvms,part(i)%vms(1),1,data(idat),1)
             idat = idat + nvms
          end if
          return
       end if
    end do

  end subroutine getStress


  !!============================================================================
  !> @brief Returns the size of von Mises stress array for the specified part.
  !>
  !> @param[in] id Base ID of the FE part to return stresses for
  !>
  !> @callergraph
  !>
  !> @author Runar Heggelien Refsnaes
  !>
  !> @date 30 Oct 2017

  function getStressSize (id) result(ndat)

    integer, intent(in) :: id

    !! Local variables
    integer :: i, ndat

    !! --- Logic section ---

    if (.not.allocated(part)) then
       !! Invoked before allocation or after deallocation
       ndat = -999
       return
    end if

    do i = 1, size(part)
       if (part(i)%baseId == id .and. associated(part(i)%vms)) then
          ndat = size(part(i)%vms)
          return
       end if
    end do

    ndat = -1

  end function getStressSize


  !!============================================================================
  !> @brief Initializes all strain rosette elements from the solver input file.
  !>
  !> @param[in] infp File unit number for the solver input file
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 15 Jan 2016

  subroutine readStrainRosettes (infp,ierr)

    use StrainRosetteModule, only : readStrainElements

    integer, intent(in)  :: infp
    integer, intent(out) :: ierr

    !! --- Logic section ---

    call readStrainElements (infp,strainRosettes,err=ierr)

  end subroutine readStrainRosettes


  !!============================================================================
  !> @brief Initiates all strain rosettes with recovery data from superelements.
  !>
  !> @param[in] ipsw Print switch for debug output
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 15 Jan 2016

  subroutine initStrainRosettes (ipsw,lpu,ierr)

    use DisplacementModule     , only : RMatrixPtr, IVectorPtr
    use DisplacementModule     , only : elDispFromSupElDisp
    use StrainRosetteModule    , only : initStrainRosette, checkRosette
    use StrainRosetteModule    , only : initStrainGauges, printRosetteHeading
    use ReportErrorModule      , only : allocationError
    use ReportErrorModule      , only : reportError, debugFileOnly_p
    use FFlLinkHandlerInterface, only : ffl_set

    integer, intent(in)  :: ipsw, lpu
    integer, intent(out) :: ierr

    !! Local variables
    integer                       :: i, j, nRosettes
    type(RMatrixPtr), allocatable :: Hmat(:)
    type(IVectorPtr), allocatable :: globalNodes(:)
    integer         , allocatable :: irec(:)

    !! --- Logic section ---

    nRosettes = size(strainRosettes)
    allocate(Hmat(nRosettes),globalNodes(nRosettes),irec(nRosettes+1),stat=ierr)
    if (ierr /= 0) then
       ierr = allocationError('initStrainRosettes')
       return
    end if

    !! Identify all strain rosettes that are on parts marked for recovery

    irec = 0
    do i = 1, nRosettes
       nullify(Hmat(i)%p)
       nullify(globalNodes(i)%p)
       do j = 1, size(part)
          if (strainRosettes(i)%linkNumber == part(j)%baseId) then
             irec(i) = j
             part(j)%nRos = part(j)%nRos + 1
             exit
          end if
       end do
       if (irec(i) > 0) then
          call ffl_set (irec(i))
          !! Check that the nodal ordering of the strain rosette yields
          !! a proper normal vector and swap the element nodes if not.
          globalNodes(i)%p => strainRosettes(i)%globalNodes
          call checkRosette (strainRosettes(i),globalNodes(i)%p,ierr)
          if (ierr < 0) goto 900
       end if
    end do

    !! Compute the H_el matrices for all strain rosettes

    i = 1
    do while (i <= nRosettes)
       j = i
       if (irec(i) > 0) then
          do j = i, nRosettes
             if (irec(j+1) /= irec(i)) then
                call elDispFromSupElDisp (part(irec(i))%sam,globalNodes(i:j), &
                     &                    Hmat(i:j),ierr,irec(i))
                if (ierr < 0) goto 900
                exit
             end if
          end do
       end if
       i = j + 1
    end do
    deallocate(globalNodes)

    !! Initialize the strain rosettes

    do i = 1, nRosettes
       if (irec(i) > 0) then

          call ffl_set (irec(i))
          call InitStrainRosette (strainRosettes(i), Hmat(i)%p, &
               &                  part(irec(i))%sam%madof, &
               &                  openFiles=ipsw>0, addSigma0=.false., &
               &                  iprint=ipsw, lpu=lpu, ierr=ierr)
          if (ierr /= 0) goto 900

          call InitStrainGauges (strainRosettes(i)%data(1)%gages, &
               &                 strainRosettes(i)%data(1)%alphaGages, &
               &                 strainRosettes(i)%posInGl)

          call PrintRosetteHeading (strainRosettes(i)%data(1), &
               &                    strainRosettes(i)%globalNodes, &
               &                    part(irec(i))%sam%minex, &
               &                    strainRosettes(i)%posInGl, &
               &                    strainRosettes(i)%id%userId, &
               &                    strainRosettes(i)%linkNumber, &
               &                    strainRosettes(i)%lpuAscii)

       end if
    end do
    deallocate(irec)

    do i = 1, nRosettes
       if (associated(Hmat(i)%p)) deallocate(Hmat(i)%p)
    end do
    deallocate(Hmat)

    return

900 call reportError (debugFileOnly_p,'initStrainRosettes')

  end subroutine initStrainRosettes


  !!============================================================================
  !> @brief Initiates all strain rosettes that are on parts marked for recovery.
  !>
  !> @param[out] ierr Error flag
  !>
  !> @details This subroutine is obsolete/not used.
  !> Retained only for the historical reason.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 17 Feb 2017

  subroutine initStrainRosettesFromCore (ierr)

    use IdTypeModule       , only : getId
    use ProgressModule     , only : lterm
    use StrainRosetteModule, only : initStrainGauges
    use ReportErrorModule  , only : allocationError
    use ReportErrorModule  , only : reportError, debugFileOnly_p

    integer, intent(out) :: ierr

    !! Local variables
    integer :: i, j

    !! --- Logic section ---

    do i = 1, size(strainRosettes)
       do j = 1, size(part)
          if (strainRosettes(i)%linkNumber == part(j)%baseId) then

             part(j)%nRos = part(j)%nRos + 1
             allocate(strainRosettes(i)%data(1)%Bcart(3,part(j)%ndim),stat=ierr)
             if (ierr /= 0) then
                ierr = allocationError('initStrainRosettesFromCore')
                return
             end if

             call loadStrainRosetteFromCore (part(j)%baseId,part(j)%nRos, &
                  &                          strainRosettes(i)%data(1)%Bcart, &
                  &                          3*part(j)%ndim,ierr)
             if (ierr == 3*part(j)%ndim) then
                call initStrainGauges (strainRosettes(i)%data(1)%gages, &
                     &                 strainRosettes(i)%data(1)%alphaGages, &
                     &                 strainRosettes(i)%posInGl)
             else
                deallocate(strainRosettes(i)%data(1)%Bcart)
                nullify(strainRosettes(i)%data(1)%Bcart)
             end if
             if (ierr < 0) then
                call reportError (debugFileOnly_p,'initStrainRosettesFromCore')
                return
             end if

             write(lterm,100) trim(getId(strainRosettes(i)%id))
100          format(7X,'Loaded recovery data for Strain Rosette',A, &
                  &    ' from in-core storage')
             exit
          end if
       end do
    end do

  end subroutine initStrainRosettesFromCore


  !!============================================================================
  !> @brief Initializes for superelement stress recovery.
  !>
  !> @param[in] sups All superelements in the model
  !> @param[in] iop Flag telling what type of results to recover, see below
  !> @param[in] trec1 Timer ID for displacement recovery setup
  !> @param[in] ipsw Print switch for debug output
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @details This subroutine performs all necessary preprocessing in order to
  !> recover internal nodal displacements and stresses during time integration.
  !> This involves allocation of additional data structures and loading the
  !> required substructure data from files associated with each superelement.
  !>
  !> What type of results to recover is controlled by the input argument @a iop,
  !> which can have the following values:
  !> - = 1: Deformations and von Mises stresses
  !> - = 2: Strain gages
  !> - = 3: As both 1 and 2
  !> - A negative value implies loading substructure data from in-core arrays
  !>   residing in a shared library, instead of reading from files.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 27 Oct 2015

  subroutine initRecovery (sups,iop,trec1,ipsw,lpu,ierr)

    use SupElTypeModule        , only : SupElType
    use IdTypeModule           , only : getId, strId
    use RDBModule              , only : nullifyRDB
    use SamModule              , only : nullifySAM, writeObject
    use SamStressModule        , only : initiateSAM
    use DisplacementModule     , only : allocateBandEmatrices, openBandEmatrices
    use DisplacementModule     , only : closeBandEmatrix, getGroupBmat
    use FileUtilitiesModule    , only : getDBGfile
    use ManipMatrixModule      , only : writeObject
    use ProgressModule         , only : writeProgress
    use TimerModule            , only : startTimer, stopTimer
    use ReportErrorModule      , only : allocationError, reportError
    use ReportErrorModule      , only : debugFileOnly_p, error_p, note_p
    use FFlLinkHandlerInterface, only : ffl_init, ffl_getElmId
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_getint

    type(SupElType), intent(in)  :: sups(:)
    integer        , intent(in)  :: iop, trec1, ipsw, lpu
    integer        , intent(out) :: ierr

    !! Local variables
    integer           :: i, irec, nnod, lpuDbg, writeVMS
    character(len=16) :: rprec

    !! --- Logic section ---

    if (rk < dp) then
       rprec = 'Single Precision'
    else
       rprec = 'Double Precision'
    end if

    write(lpu,"()")
    call reportError (note_p,'Setting up for stress recovery in '//rprec)

    ierr = 0
    irec = 0
    do i = 1, size(sups)
       if (associated(sups(i)%rcy)) irec = irec + 1
    end do

    allocate(part(irec),STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('initRecovery')
       return
    end if
    if (irec == 0) return ! No stress recovery in this simulation

    do i = 1, irec
       call nullifySAM (part(i)%sam)
       call nullifyRDB (part(i)%rdb)
       call nullifyRDB (part(i)%rdb2)
       nullify(part(i)%gno)
       nullify(part(i)%map)
       nullify(part(i)%B)
       nullify(part(i)%sv)
       nullify(part(i)%vms)
    end do

    if (iop > 0) then
       call allocateBandEmatrices (irec,ierr)
       if (ierr < 0) goto 900
    end if

    call ffa_cmdlinearg_getint ('partVMStress',writeVMS)

    irec = 0
    do i = 1, size(sups)
       if (associated(sups(i)%rcy)) then
          irec = irec + 1
          part(irec)%baseId = sups(i)%id%baseId
          part(irec)%nRos = 0
          part(irec)%ndim = sups(i)%nTotDOFs ! Superelement dimension
          if (iop == -2) cycle ! Strain rosettes have been preprocessed already

          !! Read the FE data file for this superelement
          call writeProgress (' --> Reading FE data files for Part'// &
               &              strId(sups(i)%id%userId))
          call ffl_init (trim(sups(i)%rcy%fileName(4)), &
               &         trim(sups(i)%rcy%elmGroup),ierr)
          if (ierr < 0) then
             call reportError (error_p,'Error loading FE Part data')
             goto 900
          end if

          !! Establish the SAM datastructure
          call writeProgress (' --> Loading substructure recovery files')
          call initiateSAM (sups(i)%rcy%fileName(3),part(irec)%sam,lpu,ierr)
          if (ierr < 0) then
             goto 900
          else if (ipsw > 2) then
             if (sups(i)%nExtNods >= 10) then
                lpuDbg = getDBGfile(20+i)
             else if (part(irec)%sam%nnod < 1000) then
                lpuDbg = lpu
             else
                lpuDbg = 0
             end if
             if (lpuDbg > 0) then
                write(lpuDbg,"(/)")
                call writeObject (part(irec)%sam,lpuDbg,ipsw, &
                     &            'SAM data for Part'//getId(sups(i)%id))
             end if
          end if

          !! Open the B-matrix and the generalized modes files
          call startTimer (trec1)
          call openBandEmatrices (part(irec)%sam%ndof1, &
               &                  part(irec)%sam%ndof2, &
               &                  part(irec)%sam%mpar(22), &
               &                  part(irec)%sam%mpar(31), ierr, &
               &                  sups(i)%rcy%fileName(1), &
               &                  sups(i)%rcy%fileName(2), irec)
          if (ierr < 0) goto 890

          if (mod(sups(i)%rcy%recovery,2) == 1) then
             write(lpu,"()")
             call reportError (note_p,'Stress recovery is performed'// &
                  &            ' for Part'//getId(sups(i)%id),'')

             nnod = part(irec)%sam%mpar(33)
             if (nnod > 0 .and. nnod < part(irec)%sam%nnod) then
                !! Extract the relevant part of the B-matrix
                !! for element group recovery
                call writeProgress(' --> Extracting node group recovery matrix')
                call getGroupBmat (part(irec)%sam,part(irec)%gno, &
                     &             part(irec)%map,part(irec)%B, &
                     &             irec,ipsw,lpu,ierr)
                if (iop < 3 .and. ierr == 0) then
                   !! Release the superelement disk-matrix
                   call closeBandEmatrix (irec,ierr)
                end if
                if (ierr < 0) then
                   goto 890
                else if (ipsw > 3) then
                   call writeObject (part(irec)%gno,lpu,'Recovery node group')
                   call writeObject (part(irec)%map,lpu,'Group DOF mapping')
                   call writeObject (part(irec)%B,lpu,'B-matrix for node group')
                end if
             end if

             !! Allocate arrays for internal displacement and von Mises stress
             allocate(part(irec)%sv(part(irec)%sam%ndof), STAT=ierr)
             if (ierr == 0 .and. writeVMS > 1) then
                nnod = getGroupVMSsize(part(irec)%sam)
                allocate(part(irec)%vms(nnod), STAT=ierr)
             end if
             if (ierr /= 0) then
                call stopTimer (trec1)
                ierr = allocationError('initRecovery: result vectors')
                return
             end if

          end if
          call stopTimer (trec1)
       end if
    end do

    if (abs(iop) == 2 .or. iop == 3) then

       if (associated(strainRosettes)) then
          !! Initialize the strain gages with recovery data
          call writeProgress(' --> Initializing for strain gage recovery')
          if (iop < 0) then
             call initStrainRosettesFromCore (ierr)
          else
             call initStrainRosettes (ipsw,lpu,ierr)
          end if
          if (ierr < 0) goto 900
       end if

       if (iop == 2) then
          !! Release all superelement data if doing strain gages only
          call closeRecovery (.false.,ierr)
          if (ierr < 0) goto 900
       end if

       write(lpu,"()")
       irec = 0
       do i = 1, size(sups)
          if (associated(sups(i)%rcy)) then
             irec = irec + 1
             if (part(irec)%nRos > 0 .and. sups(i)%rcy%recovery/2 == 1) then
                call reportError (note_p,'Gage recovery is performed'// &
                     &            ' for Part'//getId(sups(i)%id))
             end if
          end if
       end do

    end if

    return

890 call stopTimer (trec1)
900 call reportError (debugFileOnly_p,'initRecovery')

  contains

    !> @brief Returns the size of the in-core array with von Mises stresses.
    function getGroupVMSsize (sam) result(vmSize)
      type(SamType), intent(in) :: sam
      integer :: iel, nstrp, vmSize
      vmSize = 0
      do iel = 1, sam%nel
         if (ffl_getElmId(iel) > 0) then
            select case (sam%melcon(iel))
            case (21,23)
               nstrp =  6 ! Three-noded thin shell
            case (22,24)
               nstrp =  8 ! Four-noded thin shell
            case (41)
               nstrp = 10 ! Ten-noded tetrahedron
            case (42)
               nstrp = 15 ! Fifteen-noded wedge
            case (43)
               nstrp = 20 ! Twenty-noded hexahedron
            case (44)
               nstrp =  8 ! Eight-noded hexahedron
            case (45)
               nstrp =  4 ! Four-noded tetrahedron
            case (46)
               nstrp =  6 ! Six-noded wedge
            case default
               nstrp = -3 ! All other element types are ignored
            end select
            vmSize = vmSize + 3+nstrp
         end if
      end do
    end function getGroupVMSsize

  end subroutine initRecovery


  !!============================================================================
  !> @brief Initializes the stress- and strain gage recovery frs-files.
  !>
  !> @param[in] frsFileName List of results database file names
  !> @param[in] modelFileName Name of the FEDEM model file
  !> @param[in] sups All superelements in the model
  !> @param[in] bRat Relative buffer size for the results database files
  !> @param[in] iop Flag telling what type of results to recover,
  !> see stressrecoverymodule::initrecovery
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 3 Feb 2016

  subroutine writeRecoveryHeaders (frsFileName,modelFileName,sups,bRat,iop,ierr)

    use SupElTypeModule        , only : SupElType
    use SaveStressModule       , only : writeStressHeader
    use SaveStrainGageModule   , only : writeStrainGageHeader
    use RDBModule              , only : closeRDBfile
    use ReportErrorModule      , only : reportError, debugFileOnly_p, error_p
    use FFlLinkHandlerInterface, only : ffl_set
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_isTrue
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_getint

    character(len=*), intent(in)  :: frsFileName, modelFileName
    type(SupElType) , intent(in)  :: sups(:)
    real(dp)        , intent(in)  :: bRat
    integer         , intent(in)  :: iop
    integer         , intent(out) :: ierr

    !! Local variables
    integer :: i, irec, ifrs, writeDef, writeVMS
    logical :: lVMS

    !! --- Logic section ---

    ierr = 0
    ifrs = 0
    if (iop == 1 .or. iop == 3) then
       call ffa_cmdlinearg_getint ('partDeformation',writeDef)
       call ffa_cmdlinearg_getint ('partVMStress',writeVMS)
       lVMS = mod(writeVMS,2) > 0
    else
       writeDef = 0
       lVMS = .false.
    end if
    if (writeDef > 0 .or. lVMS) then

       !! Open stress recovery frs-files and write file header
       irec = 0
       do i = 1, size(sups)
          if (associated(sups(i)%rcy)) then
             irec = irec + 1
             if (mod(sups(i)%rcy%recovery,2) == 1) then
                ifrs = ifrs + 1
                call ffl_set (irec)
                call closeRDBfile (part(irec)%rdb,ierr)
                if (ierr /= 0) goto 900
                call writeStressHeader (getFileName(ifrs,frsFileName), &
                     &                  modelFileName,part(irec)%rdb, &
                     &                  part(irec)%sam,sups(i), &
                     &                  iDef=writeDef,lStr=(/lVMS/), &
                     &                  bufRat=min(1.0_dp,bRat),ierr=ierr)
                if (ierr /= 0) goto 900
             end if
          end if
       end do

    end if
    if (iop == 2 .or. iop == 3) then
       if (.not. ffa_cmdlinearg_isTrue('allGages')) return

       !! Open gage recovery frs-files and write file header
       irec = 0
       do i = 1, size(sups)
          if (associated(sups(i)%rcy)) then
             irec = irec + 1
             if (part(irec)%nRos > 0 .and. sups(i)%rcy%recovery/2 == 1) then
                ifrs = ifrs + 1
                call closeRDBfile (part(irec)%rdb2,ierr)
                if (ierr /= 0) goto 900
                call writeStrainGageHeader (getFileName(ifrs,frsFileName), &
                     &                      modelFileName,part(irec)%rdb2, &
                     &                      strainRosettes,bRat,ierr, &
                     &                      partId=sups(i)%id%baseId)
                if (ierr /= 0) goto 900
             end if
          end if
       end do

    end if

    return

900 call reportError (debugFileOnly_p,'writeRecoveryHeaders')

  contains

    !> @brief Extracts one file name from a comma-separated list.
    function getFileName (iNum,nameList) result(fileName)
      use kindModule  , only       :  lfnam_p
      integer         , intent(in) :: iNum
      character(len=*), intent(in) :: nameList
      character(len=lfnam_p)       :: fileName
      integer                      :: i, j, k
      if (nameList(1:1) /= '<') then
         fileName = nameList
      else
         i = 0
         j = 1
         do k = 1, iNum
            i = j + 1
            j = i + index(nameList(i:),',') - 1
            if (j < i) then
               if (k == iNum) then
                  j = len_trim(nameList)
               else
                  fileName = 'noname'
                  call reportError (error_p,'Too few frs-file names specified')
                  return
               end if
            end if
         end do
         fileName = nameList(i+1:j-2)
      end if
    end function getFileName

  end subroutine writeRecoveryHeaders


  !!============================================================================
  !> @brief Closes files and deallocates memory associated with stress recovery.
  !>
  !> @param[in] finished If .true., also deallocate private recovery containers
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 27 Oct 2015

  subroutine closeRecovery (finished,ierr)

    use SamModule              , only : deAllocateSAM
    use StrainRosetteModule    , only : deallocateRosette
    use DisplacementModule     , only : closeBandEmatrices
    use RDBModule              , only : closeRDBfile
    use ReportErrorModule      , only : reportError, debugFileOnly_p
    use FFlLinkHandlerInterface, only : ffl_done

    logical, intent(in)  :: finished
    integer, intent(out) :: ierr

    !! Local variables
    integer :: i

    !! --- Logic section ---

    call closeBandEmatrices (ierr)

    if (finished .and. associated(strainRosettes)) then
       do i = 1, size(strainRosettes)
          call deallocateRosette (strainRosettes(i))
       end do
       deallocate(strainRosettes)
       nullify(strainRosettes)
    end if

    if (allocated(part)) then
       do i = 1, size(part)
          if (associated(part(i)%gno)) then
             deallocate(part(i)%gno)
             nullify(part(i)%gno)
          end if
          if (associated(part(i)%map)) then
             deallocate(part(i)%map)
             nullify(part(i)%map)
          end if
          if (associated(part(i)%B)) then
             deallocate(part(i)%B)
             nullify(part(i)%B)
          end if
          if (associated(part(i)%sv)) then
             deallocate(part(i)%sv)
             nullify(part(i)%sv)
          end if
          if (associated(part(i)%vms)) then
             deallocate(part(i)%vms)
             nullify(part(i)%vms)
          end if
          call deAllocateSAM (part(i)%sam)
          call closeRDBfile (part(i)%rdb,ierr)
          call closeRDBfile (part(i)%rdb2,ierr)
       end do
       if (finished) deallocate(part)
    end if

    call ffl_done

    if (ierr < 0) call reportError (debugFileOnly_p,'closeRecovery')

  end subroutine closeRecovery


  !!============================================================================
  !> @brief Calculates internal displacements and stresses in superelements.
  !>
  !> @param[in] sups All superelements in the model
  !> @param[in] iStep Time increment counter
  !> @param[in] time Current simulation time
  !> @param[in] doSave If .true., save the results to results database file
  !> @param[in] trec Timer ID for recovery operations
  !> @param[in] trec2 Timer ID for displacement recovery
  !> @param[in] tsav Timer ID for file saving operations
  !> @param[in] lPrint If zero, suppress all res-file output
  !> @param[in] jpsw Print switch for debug output
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @details This subroutine is invoked once after each time/load increment
  !> in order to recover internal nodal displacements and von Mises stresses
  !> for the finite elements of the superelements during the time integration.
  !> If recovery for one (more more) element group has been requested,
  !> the internal displacements are calculated only for the nodes connected to
  !> the elements in the specified group(s), unless @a doSave is .false.
  !>
  !> @details for the finite elements
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 27 Oct 2015

  subroutine stressRecovery (sups,iStep,time,doSave, &
       &                     trec,trec2,tsav,lPrint,jpsw,lpu,ierr)

    use KindModule             , only : i8
    use SupElTypeModule        , only : SupElType
    use IdTypeModule           , only : getId
    use ScratchArrayModule     , only : realScratchArray
    use TimerModule            , only : startTimer, stopTimer
    use ReportErrorModule      , only : reportError, debugFileOnly_p
    use FFlLinkHandlerInterface, only : ffl_set
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_getint

    type(SupElType), intent(in)  :: sups(:)
    integer(i8)    , intent(in)  :: iStep
    integer        , intent(in)  :: trec, trec2, tsav, lPrint, jpsw, lpu
    real(dp)       , intent(in)  :: time
    logical        , intent(in)  :: doSave
    integer        , intent(out) :: ierr

    !! Local variables
    integer :: i, irec, ipsw, writeDef, writeVMS
    logical :: writeToFrs

    !! --- Logic section ---

    ipsw = jpsw*min(abs(lPrint),1)
    call ffa_cmdlinearg_getint ('partDeformation',writeDef)
    call ffa_cmdlinearg_getint ('partVMStress',writeVMS)
    writeToFrs = writeDef > 0 .or. mod(writeVMS,2) > 0

    irec = 0
    do i = 1, size(sups)
       if (associated(sups(i)%rcy)) then
          irec = irec + 1
          if (mod(sups(i)%rcy%recovery,2) == 1) then
             call ffl_set (irec)
             if (ipsw > 0) then
                if (iStep < 1000000_i8) then
                   write(lpu,600) 'FE Part'//trim(getId(sups(i)%id)),iStep,time
                else
                   write(lpu,601) 'FE Part'//trim(getId(sups(i)%id)),time
                end if
             end if
600          format('  Doing stress recovery on ',A, ' at time step',I6,1PE12.5)
601          format('  Doing stress recovery on ',A, ' at time',1PE12.5)

             if (doSave .and. writeToFrs) then
                call recoverAndSave (sups(i),part(irec),ierr)
             else
                call recoverNotSave (sups(i),part(irec),ierr)
             end if
             if (ierr < 0) then
                call reportError (debugFileOnly_p,'stressRecovery')
                return
             end if

             !! KMO: Something is very weird. Without this write statement,
             !! the resulting frs-file sometimes becomes corrupted if doing only
             !! stress recovery on one part while doing gage recovery on others
             if (ipsw > 0 .and. part(irec)%rdb%nBytes > 0) then
                write(lpu,"(5X,'Written',I10,' bytes')") part(irec)%rdb%nBytes
             end if
          end if
       end if
    end do
    ierr = 0

  contains

    !> @brief Recovers internal displacements and stresses without saving.
    !> @callgraph @callergraph
    subroutine recoverNotSave (sup,recp,ierr)
      use DisplacementModule  , only : calcIntDisplacements
      use StressRoutinesModule, only : calcStresses

      type(SupElType)  , intent(in)    :: sup  !< Superelement to recover for
      type(RecPartType), intent(inout) :: recp !< Recovery data for superelement
      integer          , intent(out)   :: ierr !< Error flag

      !! Calculate internal displacements for all nodes in the superelement
      call startTimer (trec)
      call startTimer (trec2)
      if (associated(sup%genDOFs)) then
         call calcIntDisplacements (recp%sam,recp%sv,sup%finit, &
              &                     sup%genDOFs%ur,ipsw,lpu,ierr,irec)
      else
         call calcIntDisplacements (recp%sam,recp%sv,sup%finit, &
              &                     (/time/),ipsw,lpu,ierr,irec)
      end if
      call stopTimer (trec2)
      call stopTimer (trec)

      if (ierr >= 0 .and. associated(recp%vms)) then
         !! Calculate von Mises stresses for the specified elements
         call calcStresses (.false.,recp%rdb,recp%sam,sv=recp%sv,vms=recp%vms, &
              &             irec=trec,iprint=ipsw,lpu=lpu,ierr=ierr)
      end if

    end subroutine recoverNotSave

    !> @brief Recovers internal displacements and stresses, and saves to file.
    !> @details If recovery on element group(s) is requested, this subroutine
    !> will calculate internal displacements only for those nodes connected to
    !> the elements in the specified element groups(s).
    !> @callgraph @callergraph
    subroutine recoverAndSave (sup,recp,ierr)
      use DisplacementModule  , only : calcIntDisplacements
      use DisplacementModule  , only : calcTotalDisplacements
      use DisplacementModule  , only : calcGroupDisplacements
      use DisplacementModule  , only : extractExternalDisp
      use StressRoutinesModule, only : calcStresses
      use RDBModule           , only : writeTimeStepDB
      use SaveStressModule    , only : writeDisplacementDB
      use ManipMatrixModule   , only : writeObject

#if FT_HAS_RECOVERY == 1
#define RSCATR SSCATR
#else
#define RSCATR DSCATR
#endif

      type(SupElType)  , intent(in)    :: sup  !< Superelement to recover for
      type(RecPartType), intent(inout) :: recp !< Recovery data for superelement
      integer          , intent(out)   :: ierr !< Error flag

      integer :: ipD, ipT, nDis
      real(rk), pointer :: work(:)

      call startTimer (trec)
      call startTimer (trec2)
      if (associated(recp%gno) .and. associated(recp%B)) then

         !! Recovery on element group(s) has been requested.
         !! Calculate local and total displacements in nodes connected to group.
         !! TODO: Component modes are not accounted for here yet.
         nDis = size(recp%B,1)
         if (writeDef > 1) then
            ipD = 1 + recp%sam%ndof
            call realScratchArray (work,recp%sam%ndof+2*nDis,ierr)
         else
            ipD = 1 + recp%sam%ndof2 ! calcGroupDisplacements uses work(ndof2)
            call realScratchArray (work,recp%sam%ndof2+nDis+1,ierr)
         end if
         if (ierr == 0) then
            ipT = ipD + nDis
            call calcGroupDisplacements (work(ipD:ipT-1),work(ipT:), &
                 &                       recp%sam,sup,recp%B,recp%gno, &
                 &                       ipsw,lpu,ierr)
            call RSCATR (nDIS,recp%map,work(ipD),recp%sv(1),1)
            if (writeDef > 1) then
               call RSCATR (nDIS,recp%map,work(ipT),work(1),1)
            end if
            if (size(recp%gno) > recp%sam%mpar(33) .and. ierr == 0) then
               !! If any of the group nodes are external,
               !! we need to extract displacements for those nodes as well
               call extractExternalDisp (recp%sv,work,recp%sam,sup,recp%gno, &
                    &                    ipsw,lpu,ierr)
            end if
            if (ipsw > 1) then
               call writeObject (recp%sv,lpu,'Expanded displacement vector')
            end if
         end if

      else

         !! Calculate internal displacements for all nodes in the superelement
         if (associated(sup%genDOFs)) then
            call calcIntDisplacements (recp%sam,recp%sv,sup%finit, &
                 &                     sup%genDOFs%ur,ipsw,lpu,ierr,irec)
         else
            call calcIntDisplacements (recp%sam,recp%sv,sup%finit, &
                 &                     (/time/),ipsw,lpu,ierr,irec)
         end if

         if (ierr == 0 .and. writeDef > 1) then
            !! Calculate the total displacement state
            call realScratchArray (work,recp%sam%ndof,ierr)
            if (ierr == 0) then
               call calcTotalDisplacements (sup,recp%sam%madof,recp%sam%minex, &
                    &                       recp%sv,work,ierr)
            end if
         end if

      end if
      call stopTimer (trec2)
      call stopTimer (trec)
      if (ierr < 0) return

      if (recp%rdb%fileNumb >= 0) then
         call startTimer (tsav)

         !! Save time step identification to results database file
         call writeTimeStepDB (recp%rdb,iStep,time,ierr)
         if (ierr >= 0) then

            !! Save expanded displacement state to results database file
            if (writeDef > 1) then
               call writeDisplacementDB (recp%rdb,recp%sam,recp%sv,svtot=work, &
                    &                    ierr=ierr,lDouble2=.false.)
            else if (writeDef > 0) then
               call writeDisplacementDB (recp%rdb,recp%sam,recp%sv, &
                    &                    ierr=ierr,lDouble2=.false.)
            end if

         end if
         call stopTimer (tsav)
         if (ierr < 0) return
      end if

      !! Calculate von Mises stresses and save to results database file
      if (mod(writeVMS,2) > 0 .and. associated(recp%vms)) then
         call calcStresses (.false.,recp%rdb,recp%sam,sv=recp%sv,vms=recp%vms, &
              &             irec=trec,isav=tsav,iprint=ipsw,lpu=lpu,ierr=ierr)
      else if (mod(writeVMS,2) > 0) then
         call calcStresses (.false.,recp%rdb,recp%sam,sv=recp%sv, &
              &             irec=trec,isav=tsav,iprint=ipsw,lpu=lpu,ierr=ierr)
      else if (associated(recp%vms)) then ! calculate without saving
         call calcStresses (.false.,recp%rdb,recp%sam,sv=recp%sv,vms=recp%vms, &
              &             irec=trec,iprint=ipsw,lpu=lpu,ierr=ierr)
      end if

    end subroutine recoverAndSave

  end subroutine stressRecovery


  !!============================================================================
  !> @brief Calculates strain gage results for the superelements.
  !>
  !> @param[in] sups All superelements in the model
  !> @param[in] iStep Time increment counter
  !> @param[in] time Current simulation time
  !> @param[in] lPrint If zero, suppress all res-file output
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 27 Oct 2015

  subroutine gageRecovery (sups,iStep,time,lPrint,lpu,ierr)

    use KindModule        , only : i8
    use SupElTypeModule   , only : SupElType
    use IdTypeModule      , only : getId
    use ManipMatrixModule , only : writeObject
#if FT_HAS_RECOVERY == 1
    use ScratchArrayModule, only : getSingleScratchArray
#endif
    use ReportErrorModule , only : reportError, debugFileOnly_p

    type(SupElType), intent(in)  :: sups(:)
    integer(i8)    , intent(in)  :: iStep
    integer        , intent(in)  :: lPrint, lpu
    real(dp)       , intent(in)  :: time
    integer        , intent(out) :: ierr

    !! Local variables
    integer           :: i, j, irec
    real(rk), pointer :: finit(:)

    !! --- Logic section ---

    ierr = 0
    if (.not. associated(strainRosettes)) return

    irec = 0
    do i = 1, size(sups)
       if (associated(sups(i)%rcy)) then
          irec = irec + 1
          if (part(irec)%nRos > 0 .and. sups(i)%rcy%recovery/2 == 1) then

             if (lPrint /= 0) then
                if (iStep < 1000000_i8) then
                   write(lpu,600) 'FE Part'//trim(getId(sups(i)%id)),iStep,time
                else
                   write(lpu,601) 'FE Part'//trim(getId(sups(i)%id)),time
                end if
             end if
600          format('  Doing gage recovery on ',A, ' at time step',I6,1PE12.5)
601          format('  Doing gage recovery on ',A, ' at time',1PE12.5)

#if FT_HAS_RECOVERY == 1
             finit => getSingleScratchArray(size(sups(i)%finit),ierr)
             if (ierr /= 0) goto 900
             do j = 1, size(finit)
                finit(j) = real(sups(i)%finit(j),rk) ! Cast to single precision
             end do
#else
             finit => sups(i)%finit
#endif
             if (lPrint > 3) then
                call writeObject (finit,lpu,'  ** Nodal deformations (Finit)')
             end if
             do j = 1, size(strainRosettes)
                if (strainRosettes(j)%linkNumber == sups(i)%id%baseId) then
                   call calcGageStrains (strainRosettes(j)%data(1),finit, &
                        &                strainRosettes(j)%lpuAscii,ierr)
                   if (ierr /= 0) goto 900
                   if (lPrint > 3) then
                      write(lpu,"('  ** Strain rosette',A/6X,1P6E13.5)") &
                           &    trim(getId(strainRosettes(j)%id)), &
                           &    strainRosettes(j)%data(1)%epsC, &
                           &    strainRosettes(j)%data(1)%sigmaC
                   end if
                end if
             end do

          end if
       end if
    end do

    return

900 call reportError (debugFileOnly_p,'gageRecovery')

  contains

    !> @brief Convenience wrapper for some #strainrosettemodule subroutines.
    !> @callgraph @callergraph
    subroutine calcGageStrains (rosette,finit,lpu,ierr)
      use StrainRosetteModule, only : StrainRosetteType
      use StrainRosetteModule, only : calcZeroStartRosetteStrains
      use StrainRosetteModule, only : calcRosetteStrains
      use StrainRosetteModule, only : printRosetteStrains
      type(StrainRosetteType), intent(inout) :: rosette
      real(rk)               , intent(in)    :: finit(:)
      integer                , intent(in)    :: lpu
      integer                , intent(out)   :: ierr
      if (rosette%zeroInit) then
         call calcZeroStartRosetteStrains (rosette,finit,ierr)
         if (ierr /= 0) return
      end if
      call calcRosetteStrains (rosette,finit,ierr)
      call printRosetteStrains (rosette,time,lpu)
    end subroutine calcGageStrains

  end subroutine gageRecovery


  !!============================================================================
  !> @brief Save strain gage results to the results database file.
  !>
  !> @param[in] sups All superelements in the model
  !> @param[in] iStep Time increment counter
  !> @param[in] time Current simulation time
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 14 May 2020

  subroutine saveGageResults (sups,iStep,time,ierr)

    use KindModule          , only : i8
    use SupElTypeModule     , only : SupElType
    use SaveStrainGageModule, only : writeStrainGageDB
    use ReportErrorModule   , only : reportError, debugFileOnly_p

    type(SupElType), intent(in)  :: sups(:)
    integer(i8)    , intent(in)  :: iStep
    real(dp)       , intent(in)  :: time
    integer        , intent(out) :: ierr

    !! Local variables
    integer :: i, irec

    !! --- Logic section ---

    ierr = 0
    if (.not. associated(strainRosettes)) return

    irec = 0
    do i = 1, size(sups)
       if (associated(sups(i)%rcy)) then
          irec = irec + 1
          if (part(irec)%nRos > 0 .and. part(irec)%rdb2%fileNumb >= 0) then
             call writeStrainGageDB (part(irec)%rdb2,strainRosettes, &
                  &                  iStep,time,.false.,ierr, &
                  &                  partId=sups(i)%id%baseId)
             if (ierr /= 0) then
                call reportError (debugFileOnly_p,'saveGageResults')
                return
             end if
          end if
       end if
    end do

  end subroutine saveGageResults


  !!============================================================================
  !> @brief Flushes the strain gage results database file to disk.
  !>
  !> @param[in] tsav Timer ID for file saving operations
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 3 Feb 2016

  subroutine flushRecoveryFiles (tsav,ierr)

    use RDBModule        , only : flushRDBfile
    use TimerModule      , only : startTimer, stopTimer
    use ReportErrorModule, only : reportError, debugFileOnly_p

    integer, intent(in)  :: tsav
    integer, intent(out) :: ierr

    !! Local variables
    integer :: i

    !! --- Logic section ---

    ierr = 0

    if (allocated(part)) then
       call startTimer (tsav)
       do i = 1, size(part)
          call flushRDBfile (part(i)%rdb2,ierr)
       end do
       call stopTimer (tsav)
    end if

    if (ierr < 0) call reportError (debugFileOnly_p,'flushRecoveryFiles')

  end subroutine flushRecoveryFiles


  !!============================================================================
  !> @brief Returns the file names of the gage recovery results database files.
  !>
  !> @param[out] chnames List of results database file names
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 13 May 2020

  subroutine getGageRecoveryFiles (chnames)

    character(len=*), intent(out) :: chnames(:)

    !! Local variables
    integer :: i, ifil

    !! --- Logic section ---

    chnames = ''
    if (allocated(part)) then
       ifil = 1
       do i = 1, size(part)
          if (ifil > size(chnames)) return
          chnames(ifil) = part(i)%rdb2%fileName
          if (chnames(ifil) /= '') ifil = ifil + 1
       end do
    end if

  end subroutine getGageRecoveryFiles


  !!============================================================================
  !> @brief Writes a Fortran90 subroutine with the strain recovery matrices.
  !>
  !> @param[in] sups All superelements in the model
  !> @param[out] ierr Error flag
  !>
  !> @details This subroutine is obsolete/not used.
  !> Retained only for the historical reason.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 17 Feb 2017

  subroutine writeRosettes2Ftn (sups,ierr)

    use SupElTypeModule    , only : SupElType
    use IdTypeModule       , only : getId
    use fileUtilitiesModule, only : findUnitNumber
    use progressModule     , only : lterm
    use reportErrorModule  , only : reportError, error_p

    type(SupElType), intent(in)  :: sups(:)
    integer        , intent(out) :: ierr

    !! Local variables
    integer          :: iftn, iinc, i, j, n, irec, iros
    character(len=1) :: c

    !! --- Logic section ---

    if (.not. associated(strainRosettes)) return

    iftn = findUnitNumber(20)
    open(iftn,FILE='loadRecoveryMatrix.f90',STATUS='UNKNOWN',IOSTAT=ierr)
    if (ierr == 0) then
       iinc = findUnitNumber(21)
       open(iinc,FILE='recovery.inc',STATUS='UNKNOWN',IOSTAT=ierr)
    end if
    if (ierr /= 0) then
       call reportError (error_p,'Could not open Fortran output files', &
            &            addString='writeRosettes2Ftn')
       return
    end if

    n = 0
    irec = 0
    do j = 1, size(sups)
       if (.not. associated(sups(j)%rcy)) cycle

       irec = irec + 1
       if (part(irec)%nRos < 1) cycle

       write(lterm,*) 'Writing Fortran code for Recovery on Part'// &
            &         trim(getId(sups(j)%id))
       n = n + 1
       c = i2c(n)
       if (n == 1) then
          write(iftn,600)
       else
          write(iinc,"()")
       end if
       write(iftn,610) c,c,c,c
       write(iftn,"(' else')",advance='NO')
       write(iinc,601) c,sups(j)%id%baseId,c,part(irec)%ndim,c,part(irec)%nRos
       write(iinc,602) c,c,c

       iros = 0
       do i = 1, size(strainRosettes)
          if (strainRosettes(i)%linkNumber == sups(j)%id%baseId) then
             iros = iros + 1
             call writeMatrix (size(strainRosettes(i)%data(1)%Bcart), &
                  &            strainRosettes(i)%data(1)%Bcart, &
                  &            iros == part(irec)%nRos)
          end if
       end do

    end do
    if (n > 0) then
       write(iftn,690)
       write(lterm,*) n,'superelement(s) written to loadRecoveryMatrix.f90'// &
            &         ' and recovery.inc'
    end if

    close(iftn)
    close(iinc)

600 format('subroutine LoadStrainRosetteFromCore (supId,iros,data,ndata,stat)' &
         / ' !DEC$ ATTRIBUTES DLLEXPORT :: LoadStrainRosetteFromCore', &
         / ' integer, intent(in)  :: supId, iros, ndata', &
         / ' real   , intent(out) :: data(ndata)', &
         / ' integer, intent(out) :: stat', &
         / ' include ''recovery.inc''')
601 format(' integer, parameter :: seId',A1,' = ',I8, &
         / ' integer, parameter :: ndim',A1,' = ',I8, &
         / ' integer, parameter :: nros',A1,' = ',I8 )
602 format(' real, parameter :: bmat',A1,'(3*ndim',A1,'*nros',A1,') = (/ &')
610 format(' if (supId == seId',A1,') then', &
         / '    if (iros < 1 .or. iros > nros',A1,') then', &
         / '       stat = 0', &
         / '    else', &
         / '       stat = 3*ndim',A1, &
         / '       data = bmat',A1,'(stat*(iros-1)+1:stat*iros)', &
         / '    end if')
690 format( &
         / '    stat = 0', &
         / ' end if', &
         / 'end subroutine LoadStrainRosetteFromCore')

  contains

    !> @brief Simple conversion of an integer (in range [0,36]) to a character.
    function i2c(i) result(c)
      integer, intent(in) :: i
      character(len=1) :: c
      if (i < 1) then
         c = '0'
      else if (i < 10) then
         c = char(ichar('0')+i)
      else if (i < 37) then
         c = char(ichar('A')+i-10)
      else
         c = '*'
      end if
    end function i2c

    !> @brief Writes matrix data to a formatted file.
    subroutine writeMatrix (ndat,data,lastRos)
      integer , intent(in) :: ndat
      real(rk), intent(in) :: data(ndat)
      logical , intent(in) :: lastRos
      if (ndat > 3) then
         write(iinc,601) data(1:ndat-3)
      end if
      if (lastRos) then
         write(iinc,602) data(ndat-2:)
      else
         write(iinc,601) data(ndat-2:)
      end if
601   format(1P,E16.8,',',E15.8,',',E15.8,', &')
602   format(1P,E16.8,',',E15.8,',',E15.8,' /)')
    end subroutine writeMatrix

  end subroutine writeRosettes2Ftn


  !!============================================================================
  !> @brief Returns the number of strain gage values.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 18 Dec 2017

  function getStrainGagesSize ()

    integer :: getStrainGagesSize

    !! --- Logic section ---

    if (associated(strainRosettes)) then
       getStrainGagesSize = size(strainRosettes)*3
    else
       getStrainGagesSize = 0
    end if

  end function getStrainGagesSize


  !!============================================================================
  !> @brief Saves/restores the initial gage strains to/from an in-core array.
  !>
  !> @param[in] lSave Save to array if equal to 1, restore from it if equal to 2
  !> @param data Array with strain gage values
  !> @param[out] ierr Error flag
  !>
  !> @details If @a lSave is 0, then this subroutine will do nothing,
  !> except for restoring the @a zeroInit flag to .false. for each strain gage.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 18 Dec 2017

  subroutine initGageStrains (lSave,data,ierr)

    use reportErrorModule, only : reportError, note_p, error_p

    integer , intent(in)    :: lSave
    real(dp), intent(inout) :: data(:)
    integer , intent(out)   :: ierr

    !! Local variables
    integer :: i, idat

    !! --- Logic section ---

    ierr = 0
    if (.not. associated(strainRosettes)) return

    if (lSave == 1 .and. size(data) < 3) then
       call reportError (note_p,'No strain gage array provided, saving skipped')
       return
    else if (lSave > 1) then
       ierr = size(data) - size(strainRosettes)*3
    end if

    if (ierr < 0) then
       call reportError (error_p,'Strain gage array is too short', &
            &            addString='initGageStrains')
       return
    end if

    idat = 0
    do i = 1, size(strainRosettes)
       if (lSave == 1) then
          data(idat+1:idat+3) = strainRosettes(i)%data(1)%epsCInit
       else if (lSave > 1) then
          strainRosettes(i)%data(1)%epsCInit = data(idat+1:idat+3)
       end if
       if (lSave == 0 .or. lSave > 1) then
          strainRosettes(i)%data(1)%zeroInit = .false.
       end if
       idat = idat + 3
    end do

  end subroutine initGageStrains

end module StressRecoveryModule
