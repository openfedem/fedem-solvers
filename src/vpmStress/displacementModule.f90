!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file displacementModule.f90
!> @brief Recovery of displacements at internal nodes.

!!==============================================================================
!> @brief Module with real kind parameter to use in the recovery calculations.

module RecKindModule

#if FT_HAS_RECOVERY == 1
  use KindModule, only : sp
#else
  use KindModule, only : dp
#endif

  implicit none

  private

  !> @brief Real kind (precision) to use for all recovery calculations
#if FT_HAS_RECOVERY == 1
  integer, parameter, public :: rk = sp
#else
  integer, parameter, public :: rk = dp
#endif

end module RecKindModule


!!==============================================================================
!> @brief Module with subroutines for recovery of internal displacements.
!>
!> @details This module contains a set of subroutines for calculation of
!> displacements at the internal nodes of a superelement, based on the updated
!> positions of the triads connected to the superelement, as well as the
!> generalized modes (or static gravitation modes). The module is used by all
!> four recovery modules. Some subroutines are also used by the dynamics solver,
!> when (stress and/or gage) recovery during time integration is requested.

module DisplacementModule

  use RecKindModule    , only : rk
#if FT_HAS_RECOVERY == 1
  use sDiskMatrixModule, only : DiskMatrixType
#else
  use DiskMatrixModule , only : DiskMatrixType
#endif

  implicit none

#if FT_HAS_RECOVERY == 1
#define RGEMV  SGEMV
#define RCOPY  SCOPY
#define RAXPY  SAXPY
#define RGEMV  SGEMV
#define RSCATR SSCATR
#else
#define RGEMV  DGEMV
#define RCOPY  DCOPY
#define RAXPY  DAXPY
#define RGEMV  DGEMV
#define RSCATR DSCATR
#endif

  !> @cond NO_DOCUMENTATION
  type RMatrixPtr
     real(rk), pointer :: p(:,:)
  end type RMatrixPtr

  type IVectorPtr
     integer , pointer :: p(:)
  end type IVectorPtr
  !> @endcond

  !> @brief Displacement recovery matrices associated with external nodal DOFs
  type(DiskMatrixType), save, allocatable, private :: BmatDisk(:)
  !> @brief Displacement recovery matrices associated with generalized modes
  type(DiskMatrixType), save, allocatable, private :: EmatDisk(:)

  private :: disExpand, calcTotalNodalDisplacement


contains

  !!============================================================================
  !> @brief Retrieves result database file pointers for system response data.
  !>
  !> @param[in] triads Array of triads to read response variables for
  !> @param[in] sup The superelement to read response variables for
  !> @param[out] pointers Array of result database variable pointers
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 9 Oct 2000

  subroutine readResponsePointers (triads,sup,pointers,ierr)

    use TriadTypeModule      , only : TriadType
    use SupElTypeModule      , only : SupElType
    use PointerKindModule    , only : ptr
    use IdTypeModule         , only : getId
    use reportErrorModule    , only : reportError
    use reportErrorModule    , only : error_p, warning_p, debugFileOnly_p
    use FFrExtractorInterface, only : ffr_findPtr

    type(TriadType), intent(in)  :: triads(:)
    type(SupElType), intent(in)  :: sup
    integer(ptr)   , intent(out) :: pointers(:)
    integer        , intent(out) :: ierr

    !! Local variables
    integer :: i, numGenDof, numFixed, numFound

    !! --- Logic section ---

    numFixed = 0
    pointers = 0_ptr
    do i = 1, size(triads)
       if (triads(i)%nDOFs == 6) then
          call ffr_findPtr ('Position matrix','Triad', &
               &            triads(i)%id%baseId,pointers(i))
       else if (triads(i)%nDOFs == 3) then
          call ffr_findPtr ('Position','Triad', &
               &            triads(i)%id%baseId,pointers(i))
       else
          numFixed = numFixed + 1
       end if
    end do

    if (associated(sup%genDOFs)) then
       numGenDof = sup%genDOFs%nDOFs
    else
       numGenDof = 0
    end if
    i = size(triads) + 1
    call ffr_findPtr ('Position matrix','Part',sup%id%baseId,pointers(i))
    if (numGenDof > 0) then
       call ffr_findPtr ('Generalized displacement','Part', &
            &            sup%id%baseId,pointers(i+1))
    end if

    ierr = 0
    numFound = count(pointers > 0_ptr)
    if (numFound+numFixed == size(triads)+1+min(numGenDof,1)) then
       return ! Found file pointers for all quantities
    else if (numFound == 0) then
       !! No response variables found, most likely OK but give warning
       call reportError (warning_p,'No system-level response variables found.',&
            &            'Stress recovery will be based on local deformations',&
            &            'relative to the modelling configuration of the Part.')
       ierr = 1
       return
    end if

    !! Some (but not all) variables are missing.
    !! This is unlikely and is considered an error.
    do i = 1, size(pointers)
       if (pointers(i) > 0_ptr) cycle

       if (i <= size(triads)) then
          call reportError (error_p,'Cannot find position for Triad'// &
               &            getId(triads(i)%id))
       else if (i-1 == size(triads)) then
          call reportError (error_p,'Cannot find position matrix for Part'// &
               &            getId(sup%id))
       else if (numGenDof > 0) then
          call reportError (error_p, &
               &            'Cannot find generalized displacements for Part'// &
               &            getId(sup%id))
       else
          cycle
       end if
       ierr = ierr - 1
    end do

    call reportError (debugFileOnly_p,'readResponsePointers')

  end subroutine readResponsePointers


  !!============================================================================
  !> @brief Retrieves result database pointer for superelement displacements.
  !>
  !> @param[in] sup The superelement to read displacement response variable for
  !> @param[out] disPtr Result database variable pointer for displacement vector
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 14 Jul 2021

  subroutine readDisplPointer (sup,disPtr,ierr)

    use SupElTypeModule      , only : SupElType
    use PointerKindModule    , only : ptr
    use FFrExtractorInterface, only : ffr_findPtr

    type(SupElType), intent(in)  :: sup
    integer(ptr)   , intent(out) :: disPtr
    integer        , intent(out) :: ierr

    !! --- Logic section ---

    ierr = 0
    call ffr_findPtr ('Vectors|Dynamic response|Displacement', &
         &            'Part',sup%id%baseId,disPtr)

  end subroutine readDisplPointer


  !!============================================================================
  !> @brief Reads position matrix for a superelement from the results database.
  !>
  !> @param[in] sup The superelement to read position matrix for
  !> @param[in] supPtr Result database variable pointer for the superelement
  !> @param[in] iprint Print switch
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 7 Mar 2006

  subroutine readSupElPosition (sup,supPtr,iprint,lpu,ierr)

    use SupElTypeModule      , only : SupElType
    use PointerKindModule    , only : ptr
    use IdTypeModule         , only : getId
    use reportErrorModule    , only : reportError, error_p
    use FFrExtractorInterface, only : ffr_getData

    type(SupElType), intent(inout) :: sup
    integer(ptr)   , intent(in)    :: supPtr
    integer        , intent(in)    :: iprint, lpu
    integer        , intent(out)   :: ierr

    !! Local variables
    integer :: j, k

    !! --- Logic section ---

    !! Read current superelement position
    call ffr_getData (sup%supTr(1,1),12,supPtr,ierr)
    if (ierr /= 0) then
       write(lpu,600) 12,ierr+12
       call reportError (error_p,'Error reading position matrix for Part'// &
            &            getId(sup%id),addString='readSupElPosition')
    else if (iprint > 5) then
       write(lpu,691) 'supTr',((sup%supTr(k,j),j=1,4),k=1,3)
       call flush (lpu)
    end if

600 format(' *** Mismatch between length of wanted array',i3, &
         & ' and actual variable size',i3,' on the results file')
691 format(1X,A5,' =',1P4E13.5/(8X,4E13.5))

  end subroutine readSupElPosition


  !!============================================================================
  !> @brief Reads external nodal displacements for the specified superelement.
  !>
  !> @param[in] triads Array of triads connected to the superelement
  !> @param[in] sup The superelement to read displacements for
  !> @param[in] pointers Array of result database variable pointers
  !> @param[in] iprint Print switch
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @details This subroutine reads the updated positions for the triads
  !> connected to the superelement from the results database, and calculates
  !> the associated local deformational displacements of the superelement.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 9 Oct 2000

  subroutine readSupElDisplacements (triads,sup,pointers,iprint,lpu,ierr)

    use TriadTypeModule      , only : TriadType, writeObject
    use SupElTypeModule      , only : SupElType, writeObject, buildFinit
    use PointerKindModule    , only : ptr
    use IdTypeModule         , only : getId
    use reportErrorModule    , only : reportError, error_p, debugFileOnly_p
    use FFrExtractorInterface, only : ffr_getData

    type(TriadType), intent(inout) :: triads(:)
    type(SupElType), intent(inout) :: sup
    integer(ptr)   , intent(in)    :: pointers(:)
    integer        , intent(in)    :: iprint, lpu
    integer        , intent(out)   :: ierr

    !! Local variables
    integer :: i, k, err

    !! --- Logic section ---

    !! Read current triad positions
    ierr = 0
    do i = 1, size(triads)
       if (triads(i)%nDOFs == 6) then
          call ffr_getData (triads(i)%ur(1,1),12,pointers(i),err)
       else if (triads(i)%nDOFs == 3) then
          call ffr_getData (triads(i)%ur(1,4),3,pointers(i),err)
       else
          cycle
       end if
       if (err /= 0) then
          ierr = ierr + 1
          write(lpu,600) 3*(triads(i)%nDOFs-2),err+3*(triads(i)%nDOFs-2)
          call reportError (error_p,'Error reading position for Triad'// &
               &            getId(triads(i)%id))
       else if (iprint > 3) then
          write(lpu,691) 'Triad'//getId(triads(i)%id),(triads(i)%ur(k,:),k=1,3)
          call flush (lpu)
       end if
    end do

    !! Read current superelement position
    i = size(triads) + 1
    call ffr_getData (sup%supTr(1,1),12,pointers(i),err)
    if (err /= 0) then
       ierr = ierr + 1
       write(lpu,600) 12,err+12
       call reportError (error_p,'Error reading position matrix for Part'// &
            &            getId(sup%id))
    else if (iprint > 3) then
       write(lpu,691) 'Part'//getId(sup%id),(sup%supTr(k,:),k=1,3)
       call flush (lpu)
    end if

    if (sup%genDOFs%nDOFs > 0) then

       !! Read current generalized superelement displacements
       i = size(triads) + 2
       call ffr_getData (sup%genDOFs%ur(1),size(sup%genDOFs%ur),pointers(i),err)
       if (err /= 0) then
          ierr = ierr + 1
          write(lpu,600) size(sup%genDOFs%ur),err+size(sup%genDOFs%ur)
          call reportError (error_p, &
               'Error reading generalized displacements for Part'// &
               getId(sup%id))
       end if

    end if

    if (ierr > 0) then
       call reportError (debugFileOnly_p,'readSupElDisplacements')
    else
       !! Calculate the local deformational displacements for this superelement
       call buildFinit (sup)
    end if

    if (iprint < 5) return

    call WriteObject (sup,lpu,2)
    write(lpu,'(/)')

    do i = 1, size(triads)
       call WriteObject (triads(i),lpu,2)
       write(lpu,'(/)')
    end do

600 format(' *** Mismatch between length of wanted array',i3, &
         & ' and actual variable size',i3,' on the results file')
691 format(1X,A10,' =',1P4E13.5/(13X,4E13.5))

  end subroutine readSupElDisplacements


  !!============================================================================
  !> @brief Calculates the internal displacements due to gravitation forces.
  !>
  !> @param[in] fileName Name of file with static gravitation mode shapes
  !> @param[in] gvec The gravitation vector
  !> @param[out] vii Amplitudes of static gravitation modes
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 19 Feb 2018

  subroutine readStaticDisplacements (fileName,gvec,vii,lpu,ierr)

    use kindModule        , only : dp
    use scratchArrayModule, only : getRealScratchArray
    use binaryDBInterface , only : readDoubleDB
    use reportErrorModule , only : reportError, debugFileOnly_p, error_p

    character(len=*), intent(in)  :: fileName
    real(dp)        , intent(in)  :: gvec(3)
    real(dp)        , intent(out) :: vii(:)
    integer         , intent(in)  :: lpu
    integer         , intent(out) :: ierr

    !! Local variables
    integer           :: ndof1
    real(dp), pointer :: vgi(:)

    !! --- Logic section ---

    ndof1 = size(vii)
    vgi => getRealScratchArray(ndof1*3,ierr)
    if (ierr < 0) then
       call reportError (debugFileOnly_p,'readStaticDisplacements ')
       return
    end if

    call readDoubleDB (fileName,'displacement matrix',ndof1*3,vgi(1),ierr)
    if (ierr < 0) then
       call reportError (error_p,'Failed to read matrix file '//fileName, &
            &            addString='readStaticDisplacements')
       return
    end if

    !! Multiply with the gravitation vector to find the actual displacement
    call DGEMV ('N',ndof1,3,1.0_dp,vgi(1),ndof1,gvec(1),1,1.0_dp,vii(1),1)

    write(lpu,600) trim(fileName),ndof1,minval(vii),maxval(vii)
600 format(/4X,'Reading static gravitation modes from file: ',A &
         & /7X,'ndof1 =',I10,'  range = {',1PE13.5,',',E12.5,' }')

  end subroutine ReadStaticDisplacements


  !!============================================================================
  !> @brief Allocates the B- and E-matrix file objects.
  !>
  !> @param[in] nPart Number of superelements to allocate matrices for
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 13 Jan 2016

  subroutine allocateBandEmatrices (nPart,ierr)

#if FT_HAS_RECOVERY == 1
    use sDiskMatrixModule, only : dmNullify
#else
    use DiskMatrixModule , only : dmNullify
#endif
    use reportErrorModule, only : allocationError

    integer, intent(in)  :: nPart
    integer, intent(out) :: ierr

    !! Local variables
    integer :: i

    !! --- Logic section ---

    allocate(BmatDisk(nPart),EmatDisk(nPart),STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('allocateBandEmatrices')
    else
       do i = 1, nPart
          call dmNullify (BmatDisk(i))
          call dmNullify (EmatDisk(i))
       end do
    end if

  end subroutine allocateBandEmatrices


  !!============================================================================
  !> @brief Opens the B- and E-matrix files associated with a superelement.
  !>
  !> @param[in] ndof1 Number of 1-status (internal) DOFs in superelement
  !> @param[in] ndof2 Number of 2-status (external) DOFs in superelement
  !> @param[in] ngen Number of generalized DOFs in superelement
  !> @param[in] nColSwap Number of columns in one swap section of the B-matrix
  !> @param[out] ierr Error flag
  !> @param[in] Bmatfile Name of the B-matrix file
  !> @param[in] Ematfile Name of the E-matrix file
  !> @param[in] indx Disk-matrix index for the superelement
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 11 Oct 2000

  subroutine openBandEmatrices (ndof1,ndof2,ngen,nColSwap,ierr, &
       &                        Bmatfile,Ematfile,indx)

    use kindModule         , only : lfnam_p
    use FileUtilitiesModule, only : getFileName
#if FT_HAS_RECOVERY == 1
    use sDiskMatrixModule  , only : dmNullify, dmOpen, dmGetSwapSize
#else
    use DiskMatrixModule   , only : dmNullify, dmOpen, dmGetSwapSize
#endif
    use DiskMatrixModule   , only : dmOld_p
    use reportErrorModule  , only : reportError, debugFileOnly_p
    use reportErrorModule  , only : allocationError, internalError

    integer                   , intent(in)  :: ndof1, ndof2, ngen, nColSwap
    integer                   , intent(out) :: ierr
    character(len=*), optional, intent(in)  :: Bmatfile, Ematfile
    integer         , optional, intent(in)  :: indx

    !! Local variables
    character(lfnam_p) :: chname
    integer            :: chksum, nSwap, i

    character(len=30), parameter :: dmTag_p = '#FEDEM generalized modes'

    !! --- Logic section ---

    if (present(indx)) then
       if (.not. allocated(BmatDisk)) then
          ierr = internalError('openBandEmatrices: BmatDisk not allocated')
          return
       else if (ngen > 0 .and. .not.allocated(EmatDisk)) then
          ierr = internalError('openBandEmatrices: EmatDisk not allocated')
          return
       else if (indx < 1 .or. indx > size(BmatDisk)) then
          ierr = internalError('openBandEmatrices: Part index out of range')
          return
       end if
       i = indx
    else
       allocate(BmatDisk(1),STAT=ierr)
       if (ierr == 0 .and. ngen > 0) then
          allocate(EmatDisk(1),STAT=ierr)
       end if
       if (ierr /= 0) then
          ierr = allocationError('openBandEmatrices')
          return
       end if
       call dmNullify (BmatDisk(1))
       if (ngen > 0) then
          call dmNullify (EmatDisk(1))
       end if
       i = 1
    end if

    nSwap = dmGetSwapSize()
    if (present(Bmatfile)) then
       chname = Bmatfile
    else
       call getFileName ('Bmatfile',chname,'_B.fmx')
    end if
    if (nSwap > 0) then
       !! Use the user-provided swap size
       call dmOpen (BmatDisk(i),chname,chksum,ndof1,ndof2,dmOld_p,ierr, &
            &       swapSize=nSwap)
    else if (nSwap < 0) then
       !! Use the same swap size as in the reducer
       call dmOpen (BmatDisk(i),chname,chksum,ndof1,ndof2,dmOld_p,ierr, &
            &       nColSwap=nColSwap)
    else
       !! No swapping. Store the full matrix in core
       call dmOpen (BmatDisk(i),chname,chksum,ndof1,ndof2,dmOld_p,ierr)
    end if
    if (ierr < 0) call reportError (debugFileOnly_p,'openBandEmatrices')
    if (ierr < 0 .or. ngen < 1) return

    if (present(Ematfile)) then
       chname = Ematfile
    else
       call getFileName ('eigfile',chname,'_E.fmx')
    end if
    if (nSwap > 0) then
       !! Use the user-provided swap size
       call dmOpen (EmatDisk(i),chname,chksum,ndof1,ngen,dmOld_p,ierr, &
            &       swapSize=nSwap,wantTag=dmTag_p)
    else
       !! No swapping. Store the full matrix in core
       call dmOpen (EmatDisk(i),chname,chksum,ndof1,ngen,dmOld_p,ierr, &
            &       wantTag=dmTag_p)
    end if
    if (ierr < 0) call reportError (debugFileOnly_p,'openBandEmatrices')

  end subroutine openBandEmatrices


  !!============================================================================
  !> @brief Closes the B- and E-matrix files and deallocates related arrays.
  !>
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 11 Oct 2000

  subroutine closeBandEmatrices (ierr)

#if FT_HAS_RECOVERY == 1
    use sDiskMatrixModule , only : dmClose
#else
    use DiskMatrixModule  , only : dmClose
#endif
    use ScratchArrayModule, only : releaseScratchArrays
    use reportErrorModule , only : reportError, debugFileOnly_p

    integer, intent(out) :: ierr

    !! Local variables
    integer :: i, lerr

    !! --- Logic section ---

    ierr = 0
    call releaseScratchArrays

    !! Close the B-matrix files
    if (allocated(BmatDisk)) then
       do i = 1, size(BmatDisk)
          call dmClose (BmatDisk(i),lerr)
          ierr = ierr + lerr
       end do
       deallocate(BmatDisk)
    end if

    !! Close the E-matrix files
    if (allocated(EmatDisk)) then
       do i = 1, size(EmatDisk)
          call dmClose (EmatDisk(i),ierr)
          ierr = ierr + lerr
       end do
       deallocate(EmatDisk)
    end if

    if (ierr < 0) call reportError (debugFileOnly_p,'closeBandEmatrices')

  end subroutine closeBandEmatrices


  !!============================================================================
  !> @brief Closes the B- and E-matrix files for one particular superelement.
  !>
  !> @param[in] indx Index of the matrix to close
  !> @param[out] ierr Error flag
  !>
  !> @details This subroutine is used by the dynamics solver only, when
  !> stress/displacement recovery during time integration has been requested.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 20 May 2020

  subroutine closeBandEmatrix (indx,ierr)

#if FT_HAS_RECOVERY == 1
    use sDiskMatrixModule, only : dmClose, dmNullify
#else
    use DiskMatrixModule , only : dmClose, dmNullify
#endif
    use reportErrorModule, only : reportError, debugFileOnly_p

    integer, intent(in)  :: indx
    integer, intent(out) :: ierr

    !! Local variables
    integer :: lerr

    !! --- Logic section ---

    ierr = 0

    if (allocated(BmatDisk)) then
       if (indx > 0 .and. indx <= size(BmatDisk)) then
          call dmClose (BmatDisk(indx),lerr)
          call dmNullify (BmatDisk(indx))
          ierr = ierr + lerr
       end if
    end if

    if (allocated(EmatDisk)) then
       if (indx > 0 .and. indx <= size(EmatDisk)) then
          call dmClose (EmatDisk(indx),lerr)
          call dmNullify (EmatDisk(indx))
          ierr = ierr + lerr
       end if
    end if

    if (ierr < 0) call reportError (debugFileOnly_p,'closeBandEmatrix')

  end subroutine closeBandEmatrix


  !!============================================================================
  !> @brief Reads the internal displacements for a superelement from file.
  !>
  !> @param[in] disPtr Result database variable pointer for displacement
  !> @param[in] supId Id of the superelement to read displacements for
  !> @param[in] sam Assembly management data for the superelement
  !> @param[out] sv Local nodal displacements in DOF-order
  !> @param[in] iprint Print switch. If positive, print the local displacements.
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 14 Jul 2021

  subroutine readIntDisplacements (disPtr,supId,sam,sv,iprint,lpu,ierr)

    use IdTypeModule     , only : IdType, getId
    use SamModule        , only : SamType
    use PointerKindModule, only : ptr
    use reportErrorModule, only : reportError, error_p

    integer(ptr) , intent(in)  :: disPtr
    type(IdType) , intent(in)  :: supId
    type(SamType), intent(in)  :: sam
    real(rk)     , intent(out) :: sv(sam%ndof)
    integer      , intent(in)  :: iprint, lpu
    integer      , intent(out) :: ierr

    !! Local variables
    integer :: i

    !! --- Logic section ---

    call ffr_getData (sv(1),sam%ndof,disPtr,ierr)
    if (ierr /= 0) then
       write(lpu,690) sam%ndof,ierr+sam%ndof
       call reportError (error_p,'Error reading displacements for Part'// &
            &            getId(supId),addString='readIntDisplacements')
       if (ierr < 0) ierr = -ierr
    else if (iprint > 1) then
       write(lpu,600)
       do i = 1, sam%nnod
          write(lpu,610) sam%minex(i),sv(sam%madof(i):sam%madof(i+1)-1)
       end do
    end if

600 format( / 5X,'LOCAL DISPLACEMENT IN SUPERELEMENT :' &
         & // 5X,'Node #',4X,'X-disp',7X,'Y-disp',7X,'Z-disp', &
         &                7X,'X-rot',8X,'Y-rot',8X,'Z-rot')
610 format(I9,2X,1P,6E13.5)
690 format(' *** Mismatch between length of wanted array:',i8 &
         & / 5X,'and actual variable size:',i8)

  end subroutine readIntDisplacements


  !!============================================================================
  !> @brief Extracts the internal displacements for a superelement.
  !>
  !> @param[in] sam Assembly management data for the superelement
  !> @param[out] sv Local nodal displacements in DOF-order
  !> @param[in] finit Deformational displacements at the supernodes
  !> @param[in] vg Amplitudes of the generalized modes
  !> @param[in] iprint Print switch. If positive, print the local displacements.
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !> @param[in] indx Disk-matrix index for the superelement
  !> @param[in] vii Amplitudes of static gravitation modes
  !> @param[in] g Gravitation vector
  !>
  !> @details The internal (VI) and external (VE) displacement vectors are
  !> established for the superelement. Then they are rearranged from equation-
  !> order into nodal point DOF-order (SV).
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 11 Oct 2000

  subroutine calcIntDisplacements (sam,sv,finit,vg,iprint,lpu,ierr,indx,vii,g)

    use SamModule         , only : SamType, dp
#if FT_HAS_RECOVERY == 1
    use sDiskMatrixModule , only : dmMatTimesVec
#else
    use DiskMatrixModule  , only : dmMatTimesVec
#endif
    use ScratchArrayModule, only : realScratchArray
    use reportErrorModule , only : reportError, debugFileOnly_p

    type(SamType)    , intent(in)  :: sam
    real(dp)         , intent(in)  :: finit(:), vg(:)
    real(rk)         , intent(out) :: sv(sam%ndof)
    integer          , intent(in)  :: iprint, lpu
    integer          , intent(out) :: ierr
    integer, optional, intent(in)  :: indx
    real(rk),optional, intent(in)  :: vii(:), g(3)

    !! Local variables
    integer           :: i, j, ipVe, ipVi, ndof1, ndof2, neq, ngen
    real(rk), pointer :: work(:)

    !! --- Logic section ---

    ierr  = 0
    neq   = sam%neq
    ndof1 = sam%ndof1
    ndof2 = sam%ndof2
    ngen  = sam%mpar(22)
    call realScratchArray (work,neq+ndof1+ndof2,ierr)
    if (ierr < 0) goto 900

    if (present(indx)) then
       j = indx
    else
       j = 1
    end if

    ipVi = 1 + neq
    ipVe = ipVi + ndof1

    !! Extract the external degrees of freedom
    call RCOPY (ndof2,0.0_rk,0,work(ipVe),1)
#if FT_HAS_RECOVERY == 1
    do i = 1, ndof2
       work(ipVe+sam%dofPosIn2(i)-1) = real(finit(i),rk)
    end do
#else
    call DSCATR (ndof2,sam%dofPosIn2(1),finit(1),work(ipVe),1)
#endif

    !! Calculate internal degrees of freedom
    call dmMatTimesVec (BmatDisk(j),work(ipVe:),work(ipVi:ipVe-1),ierr)

    if (ngen > 0 .and. ierr >= 0) then
       !! Add contributions from the generalized modes to vi
       call dmMatTimesVec (EmatDisk(j),vg,work(ipVi:ipVe-1),ierr,.false.)
    end if
    if (ierr < 0) goto 900

    if (present(g) .and. present(vii)) then
       !! Add contributions from static gravitation forces
       call RGEMV ('N',ndof1,3,1.0_rk,vii(1),ndof1,g(1),1,1.0_rk,work(ipVi),1)
    end if

    !! Add vi(ndof1) and ve(ndof2) into sveq(neq) "equation order"
    call RCOPY (neq,0.0_rk,0,work(1),1)
    call RSCATR (ndof1,sam%meqn1(1),work(ipVi),work(1),1)
    call RSCATR (ndof2,sam%meqn2(1),work(ipVe),work(1),1)

    !! Expand sveq(neq) to sv(ndof) "nodal point order"
    call disExpand (sam,work,sv)

    if (iprint > 1) then

       !! Output to log-file
       write(lpu,600)
       do i = 1, sam%nnod
          write(lpu,610) sam%minex(i),sv(sam%madof(i):sam%madof(i+1)-1)
       end do

600    format( / 5X,'LOCAL DISPLACEMENT IN SUPERELEMENT :' &
            & // 5X,'Node #',4X,'X-disp',7X,'Y-disp',7X,'Z-disp', &
            &                7X,'X-rot',8X,'Y-rot',8X,'Z-rot')
610    format(I9,2X,1P,6E13.5)

    end if

    return

900 call reportError (debugFileOnly_p,'calcIntDisplacements')

  end subroutine calcIntDisplacements


  !!============================================================================
  !> @brief Calculates the total displacements at the nodal points.
  !>
  !> @param[in] sup The superelement to calculate internal displacements for
  !> @param[in] madof Matrix of accumulated DOFs for the superelement
  !> @param[in] minex Matrix of internal to external node number mapping
  !> @param[in] sv Local nodal displacements in DOF-order
  !> @param[out] svTotal Total nodal displacements in DOF-order
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Mar 2018

  subroutine calcTotalDisplacements (sup,madof,minex,sv,svTotal,ierr)

    use SupElTypeModule  , only : SupElType
    use ReportErrorModule, only : reportError, debugFileOnly_p

    type(SupElType), intent(in)  :: sup
    integer        , intent(in)  :: madof(:), minex(:)
    real(rk)       , intent(in)  :: sv(:)
    real(rk)       , intent(out) :: svTotal(:)
    integer        , intent(out) :: ierr

    !! Local variables
    integer :: inod, idof, iend

    !! --- Logic section ---

    call RCOPY (size(svTotal),0.0_rk,0,svTotal(1),1)

    do inod = 1, size(madof)-1
       idof = madof(inod)
       iend = madof(inod+1) - 1
       if (iend > idof .and. iend <= size(svTotal)) then
          call calcTotalNodalDisplacement (sup,inod,minex,sv(idof:iend), &
               &                           svTotal(idof:iend),ierr)
          if (ierr /= 0) then
             call reportError (debugFileOnly_p,'calcTotalDisplacements')
             return
          end if
       end if
    end do

  end subroutine calcTotalDisplacements


  !!============================================================================
  !> @brief Calculates superelement-to-element displacement extraction matrices.
  !>
  !> @param[in] sam Assembly management data for the superelement
  !> @param[in] mnpcPtr Array of pointers to MNPC-arrays for the elements
  !> @param[out] Hel_Ptr Array of element displacement extraction matrices
  !> @param[out] ierr Error flag
  !> @param[in] indx Disk-matrix index for the superelement
  !>
  !> @details This subroutine establishes the transformation matrices that give
  !> the displacement vectors for strain rosette elements in global coordinates
  !> from the superelement displacement vector.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 11 Dec 1999

  subroutine ElDispFromSupElDisp (sam,mnpcPtr,Hel_Ptr,ierr,indx)

    use SamModule         , only : SamType
#if FT_HAS_RECOVERY == 1
    use sDiskMatrixModule , only : dmGetColPtr
#else
    use DiskMatrixModule  , only : dmGetColPtr
#endif
    use ScratchArrayModule, only : realScratchArray
    use reportErrorModule , only : allocationError, reportError, debugFileOnly_p

    type(SamType)   , intent(in)  :: sam
    type(IVectorPtr), intent(in)  :: mnpcPtr(:)
    type(RMatrixPtr), intent(out) :: Hel_Ptr(:)
    integer         , intent(out) :: ierr
    integer,optional, intent(in)  :: indx

    !! Local variables
    integer , pointer :: mnpc(:), meqn1(:), meqn2(:), dofPosIn2(:)
    real(rk), pointer :: Helement(:,:), v1(:), work(:)
    integer           :: ndim, ndof, ndof1, ndof2, neq, nelnod, nMat
    integer           :: i, j, supEldof, inod, gnod, idof, gdof

    !! --- Logic section ---

    ierr  = 0
    ndof1 = sam%ndof1
    ndof2 = sam%ndof2
    ndof  = sam%ndof
    neq   = sam%neq
    ndim  = sam%mpar(24)
    nMat  = size(mnpcPtr)
    call realScratchArray (work,neq+ndof,ierr)
    if (ierr < 0) goto 900

    if (present(indx)) then
       j = indx
    else
       j = 1
    end if

    !! Form the full displacement vector corresponding to each visible
    !! superelement dof. First establish in equation number order,
    !! then expand to all dofs.

    meqn1     => sam%meqn1
    meqn2     => sam%meqn2
    dofPosIn2 => sam%dofPosIn2

    do i = 1, nMat
       idof = 0
       mnpc => mnpcPtr(i)%p
       nelnod = size(mnpc)
       do inod = 1, nelnod
          gnod = mnpc(inod)
          idof = idof + sam%madof(gnod+1)-sam%madof(gnod)
       end do
       allocate(Hel_Ptr(i)%p(idof,ndim),stat=ierr)
       if (ierr /= 0) then
          ierr = allocationError('ElDispFromSupEldisp: H_el')
          return
       end if
       Hel_Ptr(i)%p = 0.0_rk
    end do

    call RCOPY (neq,0.0_rk,0,work(1),1)

    do supEldof = 1, ndim

       if (supEldof <= ndof2) then
          if (ndof1 > 0) then
             v1 => dmGetColPtr(BmatDisk(j),dofPosIn2(supElDof),ierr)
             if (ierr < 0) goto 900
             call RSCATR (ndof1,meqn1(1),v1(1),work(1),1)
          end if
          work(meqn2(supEldof)) = 1.0_rk
       else
          v1 => dmGetColPtr(EmatDisk(j),supElDof-ndof2,ierr)
          if (ierr < 0) goto 900
          call RSCATR (ndof1,meqn1(1),v1(1),work(1),1)
       end if

       call disExpand (sam,work,work(neq+1:))
       if (supEldof <= ndof2) work(meqn2(supEldof)) = 0.0_rk

       do i = 1, nMat
          mnpc => mnpcPtr(i)%p
          Helement => Hel_Ptr(i)%p
          nelnod = size(mnpc)

          idof = 1
          do inod = 1, nelnod
             gnod = mnpc(inod)
             gdof = sam%madof(gnod)
             ndof = sam%madof(gnod+1) - gdof
             call RCOPY (ndof,work(neq+gdof),1,Helement(idof,supEldof),1)
             idof = idof + ndof
          end do
       end do

    end do

    return

900 call reportError (debugFileOnly_p,'ElDispFromSupElDisp')

  end subroutine ElDispFromSupElDisp


  !!============================================================================
  !> @brief Expands and rearranges a displacement vector.
  !>
  !> @param[in] sam Assembly management data for the superelement
  !> @param[in] sveq Displacement vector expressed in equation order
  !> @param[out] svdof Displacement vector expressed in nodal point DOF order
  !>
  !> @details This is the same as the subroutine solextensionmodule::csexpand,
  !> but we don't want the recovery modules to depend on that big module.
  !> The implementation is essentially the same as the subroutine EXPAND
  !> from the SAM library. Notice that the first coefficient of each constraint
  !> equation, i.e., `c0=ttcc(mpmceq(iceq))` is ignored here, under the
  !> assumption that all constraint equations represents coupling elements
  !> (RBE2 and RBE3), and no prescribed motions.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date Feb 2017

  subroutine disExpand (sam,sveq,svdof)

    use SamModule, only : SamType

    type(SamType), intent(in)  :: sam
    real(rk)     , intent(in)  :: sveq(:)
    real(rk)     , intent(out) :: svdof(:)

    !! Local variables
    integer :: ICEQ, IDOF, IEQ, IP

    !! --- Logic section ---

    do IDOF = 1, sam%NDOF
       IEQ  = sam%MEQN(IDOF)
       ICEQ = -IEQ
       if (IEQ .gt. 0 .and. IEQ .le. sam%NEQ) then
          svDOF(IDOF) = svEQ(IEQ)
       else if (ICEQ .gt. 0 .and. ICEQ .le. sam%NCEQ) then
          svDOF(IDOF) = 0.0_rk
          do IP = sam%MPMCEQ(ICEQ)+1, sam%MPMCEQ(ICEQ+1)-1
             if (sam%MMCEQ(IP) .gt. 0 .and. sam%MMCEQ(IP) .le. sam%NDOF) then
                IEQ = sam%MEQN(sam%MMCEQ(IP))
                if (IEQ .gt. 0 .and. IEQ .le. sam%NEQ) then
                   svDOF(IDOF) = svDOF(IDOF) + real(sam%TTCC(IP),rk)*svEQ(IEQ)
                end if
             end if
          end do
       else
          svDOF(IDOF) = 0.0_rk
       end if
    end do

  end subroutine disExpand


  !!============================================================================
  !> @brief Expands and rearranges a displacement vector.
  !>
  !> @param[in] sam Assembly management data for the superelement
  !> @param[in] sveq Displacement vector expressed in status-1 equation order
  !> @param[in] svdof Displacement vector expressed in nodal point DOF order
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date Mar 2018

  subroutine dis1Expand (sam,sveq,svdof,ierr)

    use SamModule         , only : SamType
    use ScratchArrayModule, only : realScratchArray
    use reportErrorModule , only : reportError, debugFileOnly_p

    type(SamType), intent(in)  :: sam
    real(rk)     , intent(in)  :: sveq(:)
    real(rk)     , intent(out) :: svdof(:)
    integer      , intent(out) :: ierr

    !! Local variables
    real(rk), pointer :: work(:)

    !! --- Logic section ---

    call realScratchArray (work,sam%neq,ierr)
    if (ierr < 0) then
       call reportError (debugFileOnly_p,'dis1Expand')
    else
       call RCOPY (sam%neq,0.0_rk,0,work(1),1)
       call RSCATR (sam%ndof1,sam%meqn1(1),sveq(1),work(1),1)
       call disExpand (sam,work,svdof)
    end if

  end subroutine dis1Expand


  !!============================================================================
  !> @brief Extracts a sub-matrix from the B-matrix associated with node groups.
  !>
  !> @param[in] sam Assembly management data for the superelement
  !> @param[out] MNPC List of node numbers for node group to calculate for
  !> @param[out] MGMAP DOF-mapping from node group to global
  !> @param[out] Bsub Displacement recovery matrix for the node group
  !> @param[in] bind Index of the superelement B-matrix to extract from
  !> @param[in] iprint Print switch
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @details This subroutine is used by the dynamics solver only, when
  !> stress/displacement recovery on element groups has been requested.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 18 May 2020

  subroutine getGroupBmat (sam,MNPC,MGMAP,Bsub,bInd,iprint,lpu,ierr)

    use SamModule        , only : SamType, NodeNumFromDofNum
    use IdTypeModule     , only : StrId
#if FT_HAS_RECOVERY == 1
    use sDiskMatrixModule, only : dmSize, dmGetSwapSecPtr
#else
    use DiskMatrixModule , only : dmSize, dmGetSwapSecPtr
#endif
    use reportErrorModule, only : debugFileOnly_p, reportError
    use reportErrorModule, only : allocationError, internalError

    type(SamType)    , intent(in)  :: sam
    integer , pointer, intent(out) :: MNPC(:), MGMAP(:)
    real(rk), pointer, intent(out) :: Bsub(:,:)
    integer          , intent(in)  :: bInd, iprint, lpu
    integer          , intent(out) :: ierr

    !! Local variables
    integer           :: i, iB, ie, ieq, idof, inod, irow, iceq
    integer           :: j, jdof, ldof, lnod, mdof
    integer           :: ndof1, ndof2, ndofG, nSwap, nS
    real(rk), pointer :: Bsec(:,:)

    !! --- Logic section ---

    if (.not. allocated(BmatDisk)) then
       ierr = internalError('getGroupBmat: BmatDisk not allocated')
       return
    else if (bInd < 1 .or. bInd > size(BmatDisk)) then
       ierr = internalError('getGroupBmat: B-matrix index out of range')
       return
    end if

    allocate(MNPC(count(sam%mnnn < 1)),STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('getGroupBmat: MNPC')
       nullify(MNPC)
       nullify(Bsub)
       return
    end if

    !! Find the group nodes
    irow = 0
    ndofG = 0
    do inod = 1, sam%nnod
       if (sam%mnnn(inod) < 0) then ! Internal group node
          irow = irow + 1
          MNPC(irow) = inod
          ndofG = ndofG + sam%madof(inod+1) - sam%madof(inod)
       else if (sam%mnnn(inod) == 0) then ! External group node
          irow = irow + 1
          MNPC(irow) = inod
       end if
    end do

    !! Initialize the group B-matrix
    ndof1 = size(sam%meqn1)
    ndof2 = size(sam%meqn2)
    allocate(MGMAP(ndofG),Bsub(ndofG,ndof2),STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('getGroupBmat: Bsub')
       nullify(Bsub)
       return
    end if
    call RCOPY (ndofG*ndof2,0.0_rk,0,Bsub(1,1),1)

    write(lpu,600) trim(StrId(ndofG)),trim(StrId(ndof2)), &
         &         trim(StrId(ndof1)),trim(StrId(ndof2))
600 format(10X,'Extracting a',A,' x',A,' sub-matrix for displacement recovery',&
         / 10X,'from the',A,' x',A,' global B-matrix of the superelement.')

    !! Loop over swap sections of the B-matrix
    nSwap = int(dmSize(BmatDisk(bInd),3)) ! Number of columns in each section
    do iB = 1, ndof2, nSwap
       nS = min(nSwap,ndof2-iB+1) ! Number of columns in this section

       !! Get pointer to the current B-matrix swap section in core
       Bsec => dmGetSwapSecPtr(BmatDisk(bInd),iB,ierr)
       if (ierr < 0) then
          call reportError (debugFileOnly_p,'getGroupBmat')
          return
       else if (size(Bsec,1) /= ndof1) then
          ierr = internalError('getGroupBmat: Inconsistent SAM data')
          write(lpu,690) size(Bsec,1),ndof1
690       format(9X,'size(Bsec,1) =',I8,' does not equal NDOF1 =',I8)
          return
       end if

       !! Loop over all nodes in current group
       ierr = 0
       jdof = 0 ! Running internal DOF-counter within the node group
       do i = 1, size(MNPC)
          inod = MNPC(i)
          if (sam%mnnn(inod) == 0) cycle ! Ignore the external nodes
          !! Loop over all DOFs in current node
          do idof = sam%MADOF(inod), sam%MADOF(inod+1)-1
             jdof = jdof + 1
             ieq  = sam%MEQN(idof)
             iceq = -ieq
             if (ieq > 0 .and. ieq <= sam%NEQ) then
                !! This is a free internal DOF.
                !! Copy corresponding row of the global B-matrix.
                irow = getInternalDofPos(idof)
                if (irow > 0) then
                   call RCOPY (nS,Bsec(irow,1),ndof1,Bsub(jdof,iB),ndofG)
                else
                   ierr = ierr - 1
                end if
             else if (iceq > 0 .and. iceq <= sam%NCEQ) then
                !! This is a dependent DOF.
                !! Get contributions from the associated independent DOFs.
                do j = sam%MPMCEQ(iceq)+1, sam%MPMCEQ(iceq+1)-1
                   mdof = sam%MMCEQ(j)
                   if (mdof > 0 .and. mdof <= sam%NDOF) then
                      ieq = sam%MEQN(mdof)
                      if (ieq > 0 .and. ieq <= sam%NEQ) then
                         ie = sam%MSC(mdof)
                         if (ie > 0 .and. ie <= ndof2) then
                            !! This is an external independent DOF.
                            !! Add the associated constraint equation
                            !! coefficient directly into the Bsub matrix.
                            Bsub(jdof,ie) = Bsub(jdof,ie) + real(sam%TTCC(j),rk)
                            if (iprint > 0) then
                               ldof = mdof
                               lnod = NodeNumFromDofNum(sam,ldof,.true.)
                               write(lpu,610) mdof,lnod,ldof,idof,inod
610                            format('   * External independent DOF',I8, &
                                    & ' (Node,lDof =',I8,I2, &
                                    & ') for dependent DOF',I8,' (Node =',I8,')')
                            end if
                         else
                            !! This is an internal independent DOF.
                            !! Scale corresponding row of the global B-matrix
                            !! with the constraint equation coefficient and
                            !! add to the jdof'th row of the Bsub matrix.
                            irow = getInternalDofPos(mdof)
                            if (irow > 0) then
                               call RAXPY (nS,real(sam%TTCC(j),rk), &
                                    &      Bsec(irow,1),ndof1, &
                                    &      Bsub(jdof,iB),ndofG)
                            else
                               ierr = ierr - 1
                            end if
                         end if
                      end if
                   end if
                end do
             end if
             MGMAP(jdof) = idof ! Map group DOF to global DOF
          end do
       end do
       if (ierr < 0) then
          ierr = internalError('getGroupBmat: Inconsistent SAM data')
          return
       end if

    end do

  contains

    !> @brief Convenience function with error handling for row-index extraction.
    function getInternalDofPos (idof) result(ipos)
      integer, intent(in) :: idof
      integer :: ipos
      if (idof > 0 .and. idof <= size(sam%dofPosIn1)) then
         ipos = sam%dofPosIn1(idof)
         if (ipos > 0 .and. ipos <= ndof1) return
         write(lpu,691) idof,ipos,ndof1
      else
         write(lpu,692) idof
      end if
      ipos = -idof
691   format(' *** Error: DOF',I8,' was expected to be internal,', &
           & ' but its B-matrix row-index is',I8,' (NDOF1 =',I8,')')
692   format(' *** Error: DOF',I8,' was expected to be internal (out of range)')
    end function getInternalDofPos

  end subroutine getGroupBmat


  !!============================================================================
  !> @brief Extracts the internal displacements for a group of nodes.
  !>
  !> @param[out] vd Internal local nodal displacements
  !> @param[out] vt Internal total nodal displacements
  !> @param[in] sam Assembly management data for the superelement
  !> @param[in] sup The superelement to calculate internal displacements for
  !> @param[in] Bmat Displacement recovery matrix of node group to calculate for
  !> @param[in] nodes List of node numbers for node group to calculate for
  !> @param[in] iprint Print switch. If positive, print the local displacements.
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @details This subroutine is used by the dynamics solver only, when
  !> stress/displacement recovery on element groups has been requested.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 18 May 2020

  subroutine calcGroupDisplacements (vd,vt,sam,sup,Bmat,nodes,iprint,lpu,ierr)

    use SamModule         , only : SamType
    use SupElTypeModule   , only : SupElType
    use ScratchArrayModule, only : realScratchArray
    use ReportErrorModule , only : internalError, reportError, debugFileOnly_p

    real(rk)         , intent(out) :: vd(:), vt(:)
    type(SamType)    , intent(in)  :: sam
    type(SupElType)  , intent(in)  :: sup
    real(rk)         , intent(in)  :: Bmat(:,:)
    integer          , intent(in)  :: nodes(:), iprint, lpu
    integer          , intent(out) :: ierr

    !! Local variables
    integer           :: i, j, jdof, jend
    real(rk), pointer :: work(:)

    !! --- Logic section ---

    ierr = 0
    call realScratchArray (work,sam%ndof2,ierr)
    if (ierr < 0) goto 900

    !! Extract the external degrees of freedom
    call RCOPY (sam%ndof2,0.0_rk,0,work(1),1)
#if FT_HAS_RECOVERY == 1
    do i = 1, sam%ndof2
       work(sam%dofPosIn2(i)) = real(sup%finit(i),rk)
    end do
#else
    call DSCATR (sam%ndof2,sam%dofPosIn2(1),sup%finit(1),work(1),1)
#endif

    !! Calculate internal degrees of freedom
    call RGEMV ('N',size(vd),sam%ndof2,1.0_rk,Bmat(1,1),size(Bmat,1), &
         &      work(1),1,0.0_rk,vd(1),1)

    if (sam%mpar(22) > 0) then
       !! Add contributions from the generalized modes to vi
       !! TODO: Consider this for later...
       ierr = internalError('calcGroupDisplacements: ngen > 0 not implemented')
       return
    end if

    call RCOPY (size(vt),0.0_rk,0,vt(1),1)

    if (iprint > 0) write(lpu,600)

    !! Calculate the total displacements accounting for superelement rotation

    jdof = 1
    do i = 1, size(nodes)
       j = nodes(i)
       if (sam%mnnn(j) < 0) then ! Ignore the external nodes
          jend = jdof-1 + sam%madof(j+1) - sam%madof(j)
          if (iprint > 0) then
             write(lpu,610) sam%minex(j),vd(jdof:jend)
          end if
          if (jend > jdof .and. jend <= size(vt)) then
             call calcTotalNodalDisplacement (sup,j,sam%minex, &
                  &                           vd(jdof:jend),vt(jdof:jend),ierr)
             if (ierr /= 0) goto 900
          end if
          jdof = jend + 1
       end if
    end do

    return

600 format( / 5X,'LOCAL DISPLACEMENT IN SUPERELEMENT :' &
         & // 5X,'Node #',4X,'X-disp',7X,'Y-disp',7X,'Z-disp', &
         &                7X,'X-rot',8X,'Y-rot',8X,'Z-rot')
610 format(I9,2X,1P,6E13.5)

900 call reportError (debugFileOnly_p,'calcGroupDisplacements')

  end subroutine calcGroupDisplacements


  !!============================================================================
  !> @brief Extracts local and total displacements for external group nodes.
  !>
  !> @param sv Local nodal displacements in DOF-order
  !> @param svTotal Total nodal displacements in DOF-order
  !> @param[in] sam Assembly management data for the superelement
  !> @param[in] sup The superelement to extract displacements for
  !> @param[in] nodes List of node numbers for node group to calculate for
  !> @param[in] iprint Print switch. If positive, print the local displacements.
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @details This subroutine is used by the dynamics solver only, when
  !> stress/displacement recovery on element groups has been requested.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 19 Jun 2020

  subroutine extractExternalDisp (sv,svTotal,sam,sup,nodes,iprint,lpu,ierr)

    use SamModule        , only : SamType
    use SupElTypeModule  , only : SupElType
    use ReportErrorModule, only : reportError, debugFileOnly_p

    real(rk)       , intent(inout) :: sv(:)
    real(rk)       , intent(inout) :: svTotal(:)
    type(SamType)  , intent(in)    :: sam
    type(SupElType), intent(in)    :: sup
    integer        , intent(in)    :: nodes(:), iprint, lpu
    integer        , intent(out)   :: ierr

    !! Local variables
    integer :: i, idof, ind2, j, jdof, jend

    !! --- Logic section ---

    ierr = 0
    do i = 1, size(nodes)
       j = nodes(i)
       if (sam%mnnn(j) == 0) then ! Ignore internal nodes
          jdof = sam%madof(j)
          jend = sam%madof(j+1)-1
          do idof = jdof, jend
             ind2 = sam%msc(idof)
             if (ind2 > 0 .and. ind2 <= sam%ndof2) then
                sv(idof) = real(sup%finit(ind2),rk)
             end if
          end do
          if (iprint > 0) then
             write(lpu,610) sam%minex(j),sv(jdof:jend)
          end if
          if (jend > jdof .and. jend <= size(svTotal)) then
             call calcTotalNodalDisplacement (sup,j,sam%minex,sv(jdof:jend), &
                  &                           svTotal(jdof:jend),ierr)
             if (ierr /= 0) then
                call reportError (debugFileOnly_p,'extractExternalDisp')
                return
             end if
          end if
       end if
    end do

610 format(I9,2X,1P,6E13.5)

  end subroutine extractExternalDisp


  !!============================================================================
  !> @brief Calculates the total displacement at a nodal point.
  !>
  !> @param[in] sup The superelement to calculate displacement for
  !> @param[in] inod Internal node number to calculate displacement at
  !> @param[in] minex Matrix of internal to external node number mapping
  !> @param[in] uLoc Local nodal displacement w.r.t. rotated superelement system
  !> @param[out] uTot Total nodal displacement w.r.t. global system
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 May 2020

  subroutine calcTotalNodalDisplacement (sup,inod,minex,uLoc,uTot,ierr)

    use SupElTypeModule        , only : SupElType, dp
    use ManipMatrixModule      , only : matmul34
    use RotationModule         , only : vec_to_mat, deltaRot
    use ReportErrorModule      , only : reportError, error_p
    use FFlLinkHandlerInterface, only : ffl_getNodalCoor

    type(SupElType), intent(in)  :: sup
    integer        , intent(in)  :: inod, minex(:)
    real(rk)       , intent(in)  :: uLoc(:)
    real(rk)       , intent(out) :: uTot(:)
    integer        , intent(out) :: ierr

    !! Local variables
    integer  :: i, nodeNum
    real(dp) :: X0(3), Xn(3), dX(3), dR(3,3), dRot(3)

    !! --- Logic section ---

    nodeNum = -minex(inod)
    if (nodeNum < 0) nodeNum = inod
    call ffl_getNodalCoor (X0(1),X0(2),X0(3),nodeNum,ierr)
    if (ierr /= 0) then
       call reportError (error_p,'Failed to retrieve nodal coordinates', &
            &            addString='calcTotalNodalDisplacement')
       return
    end if

    do i = 1, 3
       dX(i) = real(uLoc(i),dp)
    end do
    Xn = matmul34(sup%supTr,X0+dX) - matmul34(sup%supTrInit,X0)
    do i = 1, 3
       uTot(i) = real(Xn(i),rk)
    end do
    if (size(uLoc) < 6) return

#if FT_HAS_RECOVERY == 1
    do i = 1, 3
       dRot(i) = real(uLoc(3+i),dp)
    end do
    call vec_to_mat (dRot,dR)
#else
    call vec_to_mat (uLoc(4:6),dR)
#endif
    dRot = deltaRot(sup%supTrInit,matmul(dR,sup%supTr))
    do i = 1, 3
       uTot(3+i) = real(dRot(i),rk)
    end do

  end subroutine calcTotalNodalDisplacement

end module DisplacementModule
