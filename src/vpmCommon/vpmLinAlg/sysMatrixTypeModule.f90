!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file sysMatrixTypeModule.f90
!> @brief System matrix data container.

#ifdef FT_HAS_MKL
include 'mkl_pardiso.f90'
#endif

!!==============================================================================
!> @brief Module with kind-parameters for the SPR equation solver.
!> @details This module contains some integer kind parameters that are used by
!> the SPR equation solver to facilitate easy creation of a 64-bit integer
!> version, when needed.

module sprKindModule

  use KindModule, only : dp, i4, i8, nbi_p

  implicit none

#ifdef FT_HAS_SPR_INT8
  integer, parameter :: ik     = i8      !< 64-bit integer kind
  integer, parameter :: nbik_p = 2*nbi_p !< Number of bytes per integer(ik)
#else
  integer, parameter :: ik     = i4      !< 32-bit integer kind
  integer, parameter :: nbik_p = nbi_p   !< Number of bytes per integer(ik)
#endif

  !> @brief Casting an array of integer to/from 64-bit.
  interface castInt
     module procedure castToIK
     module procedure castFromIK
  end interface castInt


contains

  !> @brief Casts an array of integer from default to 64-bit kind.
  !> @param[in] N Array size
  !> @param[in] IIN Input array
  !> @param[out] IOUT Output array
  subroutine castToIK (N,IIN,IOUT)
    integer    , intent(in)  :: N, IIN(N)
    integer(ik), intent(out) :: IOUT(N)
    integer :: i, M
    M = mod(N,7)
    do i = 1, M
       IOUT(i) = int(IIN(i),ik)
    end do
    do i = M+1, N, 7
       IOUT(i  ) = int(IIN(i  ),ik)
       IOUT(i+1) = int(IIN(i+1),ik)
       IOUT(i+2) = int(IIN(i+2),ik)
       IOUT(i+3) = int(IIN(i+3),ik)
       IOUT(i+4) = int(IIN(i+4),ik)
       IOUT(i+5) = int(IIN(i+5),ik)
       IOUT(i+6) = int(IIN(i+6),ik)
    end do
  end subroutine castToIK

  !> @brief Casts an array of integer from 64-bit to default kind.
  !> @param[in] N Array size
  !> @param[in] IIN Input array
  !> @param[out] IOUT Output array
  subroutine castFromIK (N,IIN,IOUT)
    integer    , intent(in)  :: N
    integer(ik), intent(in)  :: IIN(N)
    integer    , intent(out) :: IOUT(N)
    integer :: i, M
    M = mod(N,7)
    do i = 1, M
       IOUT(i) = int(IIN(i))
    end do
    do i = M+1, N, 7
       IOUT(i  ) = int(IIN(i  ))
       IOUT(i+1) = int(IIN(i+1))
       IOUT(i+2) = int(IIN(i+2))
       IOUT(i+3) = int(IIN(i+3))
       IOUT(i+4) = int(IIN(i+4))
       IOUT(i+5) = int(IIN(i+5))
       IOUT(i+6) = int(IIN(i+6))
    end do
  end subroutine castFromIK

end module sprKindModule


!!==============================================================================
!> @brief Module with data types and utility subroutines for system matrices.
!> @details This module contains data types representing system matrices with
!> various storage schemes, and associated utility subroutines for data access.
!> In order of increased memory usage, the following formats are covered:
!> - Diagonal
!> - Sparse
!> - Skyline
!> - Dense

module SysMatrixTypeModule

  use sprKindModule , only : dp, ik
  use FEData        , only : FEDataInput
  use FELinearSolver, only : SAM, CAM, GSFColSup, MessageLevel, Error_Flag
#ifdef FT_HAS_MKL
  use MKL_pardiso   , only : MKL_PARDISO_HANDLE
#endif

  implicit none

#if defined(irix64) || defined(hpux64) || defined(aix64) || defined(win64) || defined(linux64)
  integer, parameter :: nbp_p = 8 !< Number of bytes in a pointer variable
#else
  integer, parameter :: nbp_p = 4 !< Number of bytes in a pointer variable
#endif

  integer, parameter :: diagonalMatrix_p = 0 !< Diagonal format
  integer, parameter :: skylineMatrix_p  = 1 !< Skyline format
  integer, parameter :: sparseMatrix_p   = 2 !< Sparse format using SPR solver
  integer, parameter :: denseMatrix_p    = 3 !< Dense format using LAPACK solver
  integer, parameter :: outOfCore_p      = 4 !< Sparse format using GSF solver
  integer, parameter :: pardiso_p        = 5 !< Sparse format using PARDISO

  !> System matrix type springs.
  character(len=18), parameter :: matrixType_p(0:5) = (/&
       &                         'diagonal matrix   ', &
       &                         'skyline matrix    ', &
       &                         'sparse matrix     ', &
       &                         'dense matrix      ', &
       &                         'out of core matrix', &
       &                         'sparse matrix     '  /)

  integer(ik), parameter :: nspar_p = 60_ik !< Size of control array @a MSPAR

  type(Error_Flag), pointer, save :: GSFinfo !< Error flag for the GSF solver


  !> @brief Data type for skyline matrix storage.
  type SkylineStorageType
     integer          :: eigensolver !< Eigenvalue solver, 1=LANCZ1, 2=LANCZ2
     integer          :: neq         !< Number of EQuations
     integer, pointer :: msky(:)     !< Matrix of SKYline definitions
  end type SkylineStorageType

  !> @brief Data type for sparse matrix storage (SPR solver).
  type SparseStorageType
     integer(ik), pointer :: mspar(:)  !< Matrix of Sparse PARameters
     integer(ik), pointer :: msica(:)  !< Matrix of Storage Information for CA
     integer(ik), pointer :: msifa(:)  !< Matrix of Storage Information for FA
     integer(ik), pointer :: mtrees(:) !< Matrix of elimination assembly TREES
     integer(ik), pointer :: mvarnc(:) !< Matrix of VARiable to Node Correspond.
     real(dp)             :: rinfo(4)  !< FLOPS count for numerical operations
  end type SparseStorageType

  !> @brief Data type for sparse matrix storage (GSF solver).
  type GSFStorageType
     logical :: isFactorized !< Set to .true. when matrix is factorized
     Type(MessageLevel), pointer :: Msg !< Object for error message generation
     Type(FEDataInput) , pointer :: P   !< Object with FE input data
     Type(SAM)         , pointer :: Q   !< Object with symbolic assembly data
     Type(CAM)         , pointer :: S   !< Object with coefficient assembly data
     Type(GSFColSup)   , pointer :: T   !< Object with the factorized matrix
  end type GSFStorageType

  !> @brief Data type for sparse matrix storage (PARDISO solver).
  type PardisoStorageType
     integer :: mtype     !< The matrix type
     integer :: iparm(64) !< Various solution parameters of PARDISO
#ifdef FT_HAS_MKL
     type(MKL_PARDISO_HANDLE) :: pt(64) !< Internal data structure pointers
#endif
     integer, pointer :: ia(:) !< Indices to first column index of each row
     integer, pointer :: ja(:) !< Column indices for the matrix elements
  end type PardisoStorageType

  !> @brief Data type for a system coefficient matrix.
  type SysMatrixType
     logical :: lFirst      !< If .true. this is the 1st time this is assembled
     logical :: isShared    !< Is the data structure shared by another matrix?
     integer :: storageType !< (0=diagonal, 1=skyline, 2=sparse, 3=dense, 4=gsf)
     integer :: dim         !< matrix DIMension, i.e., dim*dim is full storage
     integer                 , pointer :: meqn(:)  !< EQ. Numbers for all dofs
     integer(ik)             , pointer :: meqn8(:) !< 64-bit integer version
     real(dp)                , pointer :: value(:) !< The matrix element VALUES
     integer                 , pointer :: ipiv(:)  !< Dense matrix pivot indices
     type(SkylineStorageType), pointer :: skyline  !< SKYLINE data structure
     type(SparseStorageType) , pointer :: sparse   !< SPARSE data structure
     type(GSFStorageType)    , pointer :: gsf      !< GSF(DNVS) data structure
     type(PardisoStorageType), pointer :: pardiso  !< Pardiso data structure
  end type SysMatrixType


  !> @brief Reallocates a data object.
  interface reAllocate
     module procedure reAllocateSysMat
     module procedure reAllocateSkyline
     module procedure reAllocateSparse
     module procedure reAllocateGSF
     module procedure reAllocatePardiso
  end interface

  !> @brief Performs consistency checking of a data object.
  interface check
     module procedure checkSysMatStorage
  end interface

  !> @brief Standard routine for writing an object to file.
  interface writeObject
     module procedure writeSysMat
     module procedure writeSysMat2
     module procedure writeSkyMat
     module procedure writeSparseMat
     module procedure writeSparseStructure
     module procedure writePardisoMat
  end interface


  private :: reAllocateSysMat, reAllocateSkyline
  private :: reAllocateSparse, reAllocateGSF, reAllocatePardiso
  private :: checkSysMatStorage, checkSkylineStorage
  private :: checkSparseStorage, checkGSFStorage, checkPardisoStorage
  private :: writeSysMat, writeSkyMat, writePardisoMat
  private :: writeSparseMat, writeSparseStructure
  private :: FEDataInput, SAM, CAM, GSFColSup, MessageLevel, Error_Flag


contains

  !!============================================================================
  !> @brief Checks if a system matrix object is allocated or not.
  !>
  !> @param[in] this The sysmatrixtypemodule::sysmatrixtype object to check
  !>
  !> @callergraph

  function isAllocated (this)

    type(SysMatrixType), intent(in) :: this
    logical :: isAllocated

    !! --- Logic section ---

    if (this%storageType == outOfCore_p) then
       isAllocated = associated(this%gsf)
    else
       isAllocated = associated(this%value)
    end if

  end function isAllocated


  !!============================================================================
  !> @brief Initializes a system matrix object.
  !>
  !> @param[out] this The sysmatrixtypemodule::sysmatrixtype object to nullify
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 4 Mar 2003

  subroutine nullifySysMatrix (this)

    type(SysMatrixType), intent(out) :: this

    !! --- Logic section ---

    this%lFirst = .true.
    this%isShared = .false.
    this%storageType = -1
    this%dim = 0
    nullify(this%meqn)
    nullify(this%meqn8)
    nullify(this%value)
    nullify(this%ipiv)
    nullify(this%skyline)
    nullify(this%sparse)
    nullify(this%gsf)
    nullify(this%pardiso)
    nullify(GSFinfo)

  end subroutine nullifySysMatrix


  !!============================================================================
  !> @brief Lets one system matrix object share the data structure of another.
  !>
  !> @param[out] this The sysmatrixtypemodule::sysmatrixtype object
  !> to share data structures from
  !> @param[in] that The sysmatrixtypemodule::sysmatrixtype object
  !> to initialize data structures from
  !> @param[in] needsFactorization If .true., @a this is going to be factorized
  !> @param[in] newType Matrix type of @a this
  !> @param[out] ierr Error flag
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 20 Mar 2003

  subroutine shareMatrixStructure (this,that,needsFactorization,newType,ierr)

    use FELinearSolver   , only : FEAnalyze
    use ReportErrorModule, only : internalError

    type(SysMatrixType), intent(out) :: this
    type(SysMatrixType), intent(in)  :: that
    logical, optional  , intent(in)  :: needsFactorization
    integer, optional  , intent(in)  :: newType
    integer, optional  , intent(out) :: ierr

    !! Local variables
    integer     :: lerr
    integer(ik) :: nspar

    !! --- Logic section ---

    if (present(ierr)) ierr = 0

    this%isShared = .true.
    if (present(newType)) then
       this%storageType = newType
    else
       this%storageType = that%storageType
    end if
    this%dim = that%dim
    this%meqn => that%meqn
    this%meqn8 => that%meqn8
    nullify(this%value)
    nullify(this%ipiv)

    nullify(this%skyline)
    nullify(this%sparse)
    nullify(this%gsf)
    nullify(this%pardiso)

    select case (this%storageType)

    case (skylineMatrix_p)

       this%skyline => that%skyline

    case (sparseMatrix_p)

       if (associated(that%sparse)) then
          nspar = size(that%sparse%mspar,kind=ik)
       else
          nspar = 0_ik
          lerr  = internalError('shareMatrixStructure: Incompatible matrices')
          if (present(ierr)) ierr = lerr
       end if
       call reAllocateSparse ('shareMatrixStructure',this%sparse,nspar,ierr)
       if (nspar > 0_ik .and. associated(this%sparse)) then
          this%sparse%rinfo  =  that%sparse%rinfo
          this%sparse%mspar  =  that%sparse%mspar
          this%sparse%msica  => that%sparse%msica
          this%sparse%msifa  => that%sparse%msifa
          this%sparse%mtrees => that%sparse%mtrees
          this%sparse%mvarnc => that%sparse%mvarnc
       end if

    case (outOfCore_p)

       call reAllocateGSF ('shareMatrixStructure',this%gsf,.false.,ierr)
       if (associated(this%gsf) .and. associated(that%gsf)) then
          this%gsf%Msg => that%gsf%Msg
          this%gsf%P   => that%gsf%P
          this%gsf%Q   => that%gsf%Q
          if (present(needsFactorization)) then
             if (needsFactorization) then
                call FEAnalyze (this%gsf%Msg,this%gsf%Q,this%gsf%T,GSFinfo)
                lerr = checkGSFinfo(GSFinfo)
                if (present(ierr)) ierr = lerr
             end if
          end if
       else if (associated(this%gsf)) then
          lerr = internalError('shareMatrixStructure: Incompatible matrices')
          if (present(ierr)) ierr = lerr
       end if

    case (pardiso_p)

       this%pardiso => that%pardiso

    end select

  end subroutine shareMatrixStructure


  !!============================================================================
  !> @brief Allocates storage for a system matrix.
  !>
  !> @param sysMat The sysmatrixtypemodule::sysmatrixtype object to allocate for
  !> @param ierr Error flag
  !> @param[in] lpu File unit number for res-file output
  !> @param[in] chName Heading for storage requirement print
  !>
  !> @details The storage requirement is optionally written out to @a lpu
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 4 Mar 2003

  subroutine allocateSysMatrix (sysMat,ierr,lpu,chName)

    use KindModule       , only : i4, i8
    use AllocationModule , only : reAllocate, writeStorage
    use ReportErrorModule, only : internalError

    type(SysMatrixType)       , intent(inout) :: sysMat
    integer                   , intent(inout) :: ierr
    integer         , optional, intent(in)    :: lpu
    character(len=*), optional, intent(in)    :: chName

    !! Local variables
    integer(i8) :: nwR, nwI, nwIS
#ifdef FT_DEBUG
    integer(i8) :: info(4)
    integer     :: i
#endif

    type(SparseStorageType), pointer :: spr

    !! --- Logic section ---

    nwI = 3_i8 + size(sysMat%meqn,kind=i8)

    select case (sysMat%storageType)

    case (skylineMatrix_p)

       nwR = int(sysMat%skyline%msky(sysMat%dim),i8)
       nwI = nwI + 2_i8 + size(sysMat%skyline%msky,kind=i8)

    case (sparseMatrix_p)

       spr => sysMat%sparse
       nwR = 4_i8 + int(spr%mspar(8) + spr%mspar(16),i8)
       nwIS = size(spr%mspar,kind=i8) + &
            & size(spr%msica,kind=i8) + size(spr%msifa,kind=i8) + &
            & size(spr%mtrees,kind=i8) + size(spr%mvarnc,kind=i8)
       if (ik > i4) then
          nwI = nwI + 2_i8*nwIS + size(sysMat%meqn8,kind=i8)
       else
          nwI = nwI + nwIS
       end if

    case (diagonalMatrix_p)

       nwR = int(sysMat%dim,i8)

    case (denseMatrix_p)

       nwR = int(sysMat%dim*sysMat%dim,i8)

    case (outOfCore_p)

       nwR = 0_i8 ! In-core array is not used

    case (pardiso_p)

       nwR = size(sysMat%pardiso%ja,kind=i8)
       nwI = nwI + nwR + size(sysMat%pardiso%ia,kind=i8)

    case default

       ierr = internalError('allocateSysMatrix: Invalid matrix type')
       return

    end select

    if (present(chName) .and. present(lpu)) then

       call writeStorage (chName//' ('// &
            &             trim(matrixType_p(sysMat%storageType))//')', &
            &             nwI,nwR,lpu)
#ifdef FT_DEBUG
       if (sysMat%storageType == sparseMatrix_p) then
          do i = 1, 4
             info(i) = int(spr%rinfo(i),i8)
          end do
          write(lpu,600) chName,info
       end if
#endif
    end if

    call reAllocate ('allocateSysMatrix',sysMat%value,nwR,ierr)

#ifdef FT_DEBUG
600 format(/4X,'Number of floating point operations needed for ', A &
         & /7X,'Assembly of FE matrices :',2I11 &
         & /7X,'Factorization           :', I11 &
         & /7X,'Forward and back solve  :', I11 )
#endif

  end subroutine allocateSysMatrix


  !!============================================================================
  !> @brief Deallocates all dynamic arrays for a system matrix object.
  !>
  !> @param this The sysmatrixtypemodule::sysmatrixtype object to deallocate
  !> @param ierr Error flag
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 4 Mar 2003

  subroutine deAllocateSysMatrix (this,ierr)

    use kindModule      , only : i4
    use AllocationModule, only : reAllocate

    type(SysMatrixType), intent(inout) :: this
    integer, optional  , intent(inout) :: ierr

    !! --- Logic section ---

    if (this%isShared) then
       if (associated(this%sparse)) then
          nullify(this%sparse%msica)
          nullify(this%sparse%msifa)
          nullify(this%sparse%mtrees)
          nullify(this%sparse%mvarnc)
       end if
       if (associated(this%gsf)) then
          nullify(this%gsf%Msg)
          nullify(this%gsf%P)
          nullify(this%gsf%Q)
       end if
    end if

    call reAllocate        ('deAllocateSysMatrix',this%value)
    call reAllocate        ('deAllocateSysMatrix',this%ipiv)
    call reAllocateSparse  ('deAllocateSysMatrix',this%sparse)
    call reAllocateGSF     ('deAllocateSysMatrix',this%gsf,ierr=ierr)
    call reAllocatePardiso ('deAllocateSysMatrix',this%pardiso)
    if (this%isShared) return

    call reAllocate        ('deAllocateSysMatrix',this%meqn)
    if (this%storageType == sparseMatrix_p .and. ik > i4) then
       call reAllocate     ('deAllocateSysMatrix',this%meqn8)
    end if
    call reAllocateSkyline ('deAllocateSysMatrix',this%skyline)

  end subroutine deAllocateSysMatrix


  !!============================================================================
  !> @brief Allocates, reallocates or deallocates a system matrix object.
  !>
  !> @param[in] label Text label used for logging of dynamic memory use
  !> @param this The sysmatrixtypemodule::sysmatrixtype object to allocate
  !> @param newAlloc If .true., a new system matrix is allocated
  !> @param ierr Error flag
  !>
  !> @details If reallocated, the existing contents is lost.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 27 Feb 2003

  subroutine reAllocateSysMat (label,this,newAlloc,ierr)

    use KindModule       , only : nbi_p
    use AllocationModule , only : doLogMem, logAllocMem
    use ReportErrorModule, only : allocationError

    character(len=*)   , intent(in)    :: label
    type(SysMatrixType), pointer       :: this
    logical, optional  , intent(in)    :: newAlloc
    integer, optional  , intent(inout) :: ierr

    !! Local variables
    integer :: oldSize, newSize, err

    !! --- Logic section ---

    if (associated(this)) then
       oldSize = 1
       call deAllocateSysMatrix (this)
       deallocate(this)
       nullify(this)
    else
       oldSize = 0
    end if

    if (.not. present(newAlloc)) then
       newSize = 0
    else if (newAlloc) then
       allocate(this,STAT=err)
       if (err == 0) then
          newSize = 1
          call nullifySysMatrix (this)
       else
          newSize = -1
          if (present(ierr)) ierr = ierr + allocationError(label)
       end if
    else
       newSize = 0
    end if

    if (doLogMem) call logAllocMem (label,oldSize,newSize,3*nbi_p+8*nbp_p)

  end subroutine reAllocateSysMat


  !!============================================================================
  !> @brief Allocates, reallocates or deallocates a skyline storage object.
  !>
  !> @param[in] label Text label used for logging of dynamic memory use
  !> @param this The sysmatrixtypemodule::skylinestoragetype object to allocate
  !> @param[in] newAlloc If .true., a new system matrix is allocated
  !> @param ierr Error flag
  !>
  !> @details If reallocated, the existing contents is lost.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 20 Jan 2003

  subroutine reAllocateSkyline (label,this,newAlloc,ierr)

    use KindModule       , only : nbi_p
    use AllocationModule , only : doLogMem, logAllocMem, reAllocate
    use ReportErrorModule, only : allocationError

    character(len=*)        , intent(in)    :: label
    type(SkylineStorageType), pointer       :: this
    logical, optional       , intent(in)    :: newAlloc
    integer, optional       , intent(inout) :: ierr

    !! Local variables
    integer :: oldSize, newSize, err

    !! --- Logic section ---

    if (associated(this)) then
       oldSize = 1
       call reAllocate (label,this%msky)
       deallocate(this)
       nullify(this)
    else
       oldSize = 0
    end if

    if (.not. present(newAlloc)) then
       newSize = 0
    else if (newAlloc) then
       allocate(this,STAT=err)
       if (err == 0) then
          newSize = 1
          this%eigenSolver = 2 ! Using LANCZ2 is default
          this%neq = 0
          nullify(this%msky)
       else
          newSize = -1
          if (present(ierr)) ierr = ierr + allocationError(label)
       end if
    else
       newSize = 0
    end if

    if (doLogMem) call logAllocMem (label,oldSize,newSize,2*nbi_p+nbp_p)

  end subroutine reAllocateSkyline


  !!============================================================================
  !> @brief Allocates, reallocates or deallocates a sparse storage object.
  !>
  !> @param[in] label Text label used for logging of dynamic memory use
  !> @param this The sysmatrixtypemodule::sparsestoragetype object to allocate
  !> @param[in] nspar Size of the control array @a mspar
  !> @param ierr Error flag
  !>
  !> @details If reallocated, the existing contents is lost.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 20 Jan 2003

  subroutine reAllocateSparse (label,this,nspar,ierr)

    use KindModule       , only : nbs_p
    use AllocationModule , only : doLogMem, logAllocMem, reAllocate
    use ReportErrorModule, only : allocationError

    character(len=*)       , intent(in)    :: label
    type(SparseStorageType), pointer       :: this
    integer(ik), optional  , intent(in)    :: nspar
    integer    , optional  , intent(inout) :: ierr

    !! Local variables
    integer :: oldSize, newSize, err

    !! --- Logic section ---

    if (associated(this)) then
       oldSize = 1
       call reAllocate (label,this%mspar)
       call reAllocate (label,this%msica)
       call reAllocate (label,this%msifa)
       call reAllocate (label,this%mtrees)
       call reAllocate (label,this%mvarnc)
       deallocate(this)
       nullify(this)
    else
       oldSize = 0
    end if

    if (.not. present(nspar)) then
       newSize = 0
    else if (nspar > 0_ik) then
       allocate(this,STAT=err)
       if (err == 0) then
          newSize = 1
          this%rinfo = 0.0_dp
          nullify(this%mspar)
          nullify(this%msica)
          nullify(this%msifa)
          nullify(this%mtrees)
          nullify(this%mvarnc)
       else
          newSize = -1
          if (present(ierr)) ierr = ierr + allocationError(label)
       end if
    else
       newSize = 0
    end if

    if (doLogMem) call logAllocMem (label,oldSize,newSize,4*nbs_p+5*nbp_p)
    if (newSize > 0) call reAllocate (label,this%mspar,int(nspar),ierr)

  end subroutine reAllocateSparse


  !!============================================================================
  !> @brief Allocates, reallocates or deallocates a GSF storage object.
  !>
  !> @param[in] label Text label used for logging of dynamic memory use
  !> @param this The sysmatrixtypemodule::gsfstoragetype object to allocate
  !> @param[in] newAlloc If .true., a new system matrix is allocated
  !> @param ierr Error flag
  !>
  !> @details If reallocated, the existing contents is lost.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 8 Sep 2004

  subroutine reAllocateGSF (label,this,newAlloc,ierr)

    use KindModule       , only : nbi_p
    use AllocationModule , only : doLogMem, logAllocMem
    use ReportErrorModule, only : allocationError
    use FELinearSolver   , only : FEDestroy

    character(len=*)    , intent(in)    :: label
    type(GSFStorageType), pointer       :: this
    logical, optional   , intent(in)    :: newAlloc
    integer, optional   , intent(inout) :: ierr

    !! Local variables
    integer :: oldSize, newSize, err

    !! --- Logic section ---

    if (associated(this)) then
       oldSize = 1
       if (present(newAlloc)) then
          if (present(ierr)) then
             call FEDestroy (this%Msg, this%Q, this%S, this%T, GSFinfo)
             ierr = ierr + checkGSFinfo(GSFinfo)
          else
             call FEDestroy (this%Q, this%S, this%T, GSFinfo)
             err = checkGSFinfo(GSFinfo)
          end if
       else
          if (present(ierr)) then
             call FEDestroy (this%Msg, INFO=GSFinfo)
             ierr = ierr + checkGSFinfo(GSFinfo)
          end if
       end if
       if (associated(this%Msg)) deallocate(this%Msg)
       if (associated(this%P))   deallocate(this%P)
       deallocate(this)
       nullify(this)
    else
       oldSize = 0
    end if

    if (.not. present(newAlloc)) then
       newSize = 0
    else
       allocate(this,STAT=err)
       if (err == 0 .and. newAlloc) allocate(this%Msg,STAT=err)
       if (err == 0 .and. newAlloc) allocate(this%P,STAT=err)
       if (err == 0) then
          newSize = 1
          this%isFactorized = .false.
          nullify(this%Q)
          nullify(this%S)
          nullify(this%T)
       else
          newSize = -1
          if (present(ierr)) ierr = ierr + allocationError(label)
       end if
    end if

    if (doLogMem) call logAllocMem (label,oldSize,newSize,nbi_p+5*nbp_p)

  end subroutine reAllocateGSF


  !!============================================================================
  !> @brief Allocates, reallocates or deallocates a Pardiso storage object.
  !>
  !> @param[in] label Text label used for logging of dynamic memory use
  !> @param this The sysmatrixtypemodule::pardisostoragetype object to allocate
  !> @param[in] newAlloc If .true., a new system matrix is allocated
  !> @param ierr Error flag
  !>
  !> @details If reallocated, the existing contents is lost.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 10 Feb 2016

  subroutine reAllocatePardiso (label,this,newAlloc,ierr)

    use KindModule       , only : nbi_p
    use AllocationModule , only : doLogMem, logAllocMem, reAllocate
    use ReportErrorModule, only : allocationError

    character(len=*)        , intent(in)    :: label
    type(PardisoStorageType), pointer       :: this
    logical, optional       , intent(in)    :: newAlloc
    integer, optional       , intent(inout) :: ierr

    !! Local variables
    integer :: oldSize, newSize, i, err

    !! --- Logic section ---

    if (associated(this)) then
       oldSize = 1
       call reAllocate (label,this%ia)
       call reAllocate (label,this%ja)
       deallocate(this)
       nullify(this)
    else
       oldSize = 0
    end if

    if (.not. present(newAlloc)) then
       newSize = 0
    else if (newAlloc) then
       allocate(this,STAT=err)
       if (err == 0) then
          newSize = 1
          this%iparm = 0
#ifdef FT_HAS_MKL
          do i = 1, size(this%pt)
             this%pt(i)%dummy = 0
          end do
#endif
          nullify(this%ia)
          nullify(this%ja)
       else
          newSize = -1
          if (present(ierr)) ierr = ierr + allocationError(label)
       end if
    else
       newSize = 0
    end if

    if (doLogMem) then
       i = 1 + size(this%iparm)
#ifdef FT_HAS_MKL
       i = i + size(this%pt)
#endif
       call logAllocMem (label,oldSize,newSize,i*nbi_p+2*nbp_p)
    end if

  end subroutine reAllocatePardiso


  !!============================================================================
  !> @brief Performs consistency checking of a system matrix storage structure.
  !>
  !> @param[in] sysMat The sysmatrixtypemodule::sysmatrixtype object to check
  !> @param[in] mpar Matrix of parameters
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 26 Feb 2003

  subroutine checkSysMatStorage (sysMat,mpar,lpu,ierr)

    type(SysMatrixType), intent(in)  :: sysMat
    integer            , intent(in)  :: mpar(:), lpu
    integer            , intent(out) :: ierr

    !! --- Logic section ---

    select case (sysMat%storageType)

    case (skylineMatrix_p)

       if (associated(sysMat%skyline)) then
          call checkSkylineStorage (sysMat%skyline,mpar,lpu,ierr)
       else
          ierr = -1
          write(lpu,*) '*** Error: sysMat%skyline is not allocated'
       end if

    case (sparseMatrix_p)

       if (associated(sysMat%sparse)) then
          call checkSparseStorage (sysMat%sparse,mpar,lpu,ierr)
       else
          ierr = -1
          write(lpu,*) '*** Error: sysMat%sparse is not allocated'
       end if

    case (outOfCore_p)

       if (associated(sysMat%gsf)) then
          call checkGSFStorage (sysMat%gsf,lpu,ierr)
       else
          ierr = -1
          write(lpu,*) '*** Error: sysMat%gsf is not allocated'
       end if

    case (pardiso_p)

       if (associated(sysMat%pardiso)) then
          call checkPardisoStorage (sysMat%pardiso,lpu,ierr)
       else
          ierr = -1
          write(lpu,*) '*** Error: sysMat%pardiso is not allocated'
       end if

    case (diagonalMatrix_p, denseMatrix_p)

       ierr = 0 ! Nothing to check here, always OK

    case default

       ierr = -1
       write(lpu,*) '*** Error: Invalid matrix storage type',sysMat%storageType

    end select

    !! Equation numbers
    if (.not. associated(sysMat%meqn)) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: sysMat%meqn is not allocated'
    else if (size(sysMat%meqn) < mpar(7)) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: size(sysMat%meqn) =',size(sysMat%meqn)
       write(lpu,*) '                        nceq =',mpar(7)
    end if

    if (sysMat%storageType == sparseMatrix_p) then
       if (.not. associated(sysMat%meqn8)) then
          ierr = ierr - 1
          write(lpu,*) '*** Error: sysMat%meqn8 is not allocated'
       else if (size(sysMat%meqn8) < mpar(7)) then
          ierr = ierr - 1
          write(lpu,*) '*** Error: size(sysMat%meqn8) =',size(sysMat%meqn)
          write(lpu,*) '                         nceq =',mpar(7)
       end if
    end if

  end subroutine checkSysMatStorage


  !!============================================================================
  !> @brief Performs consistency checking of the skyline data structure.
  !>
  !> @param[in] skyline The sysmatrixtypemodule::skylinestoragetype to check
  !> @param[in] mpar Matrix of parameters
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 24 Sep 1998

  subroutine checkSkylineStorage (skyline,mpar,lpu,ierr)

    type(SkylineStorageType), intent(in)  :: skyline
    integer                 , intent(in)  :: mpar(:), lpu
    integer                 , intent(out) :: ierr

    !! --- Logic section ---

    if (.not. associated(skyline%msky)) then
       ierr = -1
       write(lpu,*) '*** Error: skyline%msky is not allocated'
    else if (size(skyline%msky) < mpar(11)) then
       ierr = -1
       write(lpu,*) '*** Error: size(skyline%msky) =',size(skyline%msky)
       write(lpu,*) '                          neq =',mpar(11)
    else
       ierr = 0
       return
    end if

    write(lpu,*) '*** Error return from checkSkyline, ierr =',ierr

  end subroutine checkSkylineStorage


  !!============================================================================
  !> @brief Performs consistency checking of the sparse data structure.
  !>
  !> @param[in] sparse The sysmatrixtypemodule::sparsestoragetype to check
  !> @param[in] mpar Matrix of parameters
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 24 Sep 1998

  subroutine checkSparseStorage (sparse,mpar,lpu,ierr)

    type(SparseStorageType), intent(in)  :: sparse
    integer                , intent(in)  :: mpar(:), lpu
    integer                , intent(out) :: ierr

    !! --- Logic section ---

    ierr = 0

    if (.not. associated(sparse%mspar)) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: sparse%mspar is not allocated'
    else if (size(sparse%mspar,kind=ik) < nspar_p) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: size(sparse%mspar) =',size(sparse%mspar)
       if (size(sparse%mspar) < 4) goto 999
    end if

    if (.not. associated(sparse%msica)) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: sparse%msica is not allocated'
    else if (size(sparse%msica,kind=ik) < sparse%mspar(2)) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: size(sparse%msica) =',size(sparse%msica)
       write(lpu,*) '              sparse%mspar(2) =',sparse%mspar(2)
    end if

    if (.not. associated(sparse%msifa)) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: sparse%msifa is not allocated'
    else if (size(sparse%msifa,kind=ik) < sparse%mspar(3)) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: size(sparse%msifa) =',size(sparse%msifa)
       write(lpu,*) '              sparse%mspar(3) =',sparse%mspar(3)
    end if

    if (.not. associated(sparse%mtrees)) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: sparse%mtrees is not allocated'
    else if (size(sparse%mtrees,kind=ik) < sparse%mspar(4)) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: size(sparse%mtrees) =',size(sparse%mtrees)
       write(lpu,*) '               sparse%mspar(4) =',sparse%mspar(4)
    end if

    if (.not. associated(sparse%mvarnc)) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: sparse%mvarnc is not allocated'
    else if (size(sparse%mvarnc) < 2*mpar(11)) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: size(sparse%mtrees) =',size(sparse%mtrees)
       write(lpu,*) '                         2*neq =',2*mpar(11)
    end if

999 continue
    if (ierr < 0) write(lpu,*) 'Error return from checkSparse, ierr =',ierr

  end subroutine checkSparseStorage


  !!============================================================================
  !> @brief Performs consistency checking of the sparse data structure.
  !>
  !> @param[in] sparse The sysmatrixtypemodule::pardisostoragetype to check
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 10 Feb 2016

  subroutine checkPardisoStorage (sparse,lpu,ierr)

    type(PardisoStorageType), intent(in)  :: sparse
    integer                 , intent(in)  :: lpu
    integer                 , intent(out) :: ierr

    !! Local variables
    integer :: neq

    !! --- Logic section ---

    ierr = 0

    if (.not. associated(sparse%ia)) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: sparse%ia is not allocated'
    end if
    if (.not. associated(sparse%ja)) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: sparse%ja is not allocated'
    end if

    neq = size(sparse%ia)
    if (ierr == 0 .and. size(sparse%ja) < sparse%ia(neq+1)) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: size(sparse%ja) =',size(sparse%ja)
       write(lpu,*) '          sparse%ia(neq+1) =',sparse%ia(neq+1)
    end if

    if (ierr < 0) write(lpu,*) '*** Error return from checkSparse, ierr =',ierr

  end subroutine checkPardisoStorage


  !!============================================================================
  !> @brief Performs consistency checking of the GSF data structure.
  !>
  !> @param[in] gsf The sysmatrixtypemodule::gsfstoragetype to check
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 13 Sep 2004

  subroutine checkGSFStorage (gsf,lpu,ierr)

    type(GSFStorageType), intent(in)  :: gsf
    integer             , intent(in)  :: lpu
    integer             , intent(out) :: ierr

    !! --- Logic section ---

    ierr = 0

    if (.not. associated(gsf%Msg)) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: gsf%Msg is not allocated'
    end if

    if (.not. associated(gsf%P)) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: gsf%P is not allocated'
    end if

    if (.not. associated(gsf%Q)) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: gsf%Q is not allocated'
    end if

    if (ierr < 0) write(lpu,*) 'Error return from checkGSF, ierr =',ierr

  end subroutine checkGSFStorage


  !!============================================================================
  !> @brief Inspects the GSF error flag and outputs the error messages, if any.
  !>
  !> @param[in] gsf The sysmatrixtypemodule::gsfstoragetype to check
  !> @param INFO Error handling object from the GSF package
  !> @return Negative value if an error was detected, otherwise zero
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 17 Jul 2005

  function checkGSFinfo (INFO,gsf) result(ierr)

    use FELinearSolver   , only : FEErrorHandler
    use ReportErrorModule, only : getErrorFile

    type(Error_Flag)              , pointer    :: INFO
    type(GSFStorageType), optional, intent(in) :: gsf
    integer                                    :: ierr

    ierr = 0
    if (.not. associated(INFO)) return

    !! --- Logic section ---

    if (present(gsf)) then
       call FEErrorHandler (INFO, gsf%P, gsf%Q, getErrorFile())
    else
       call FEErrorHandler (INFO, output=getErrorFile())
    end if

    nullify(INFO)
    ierr = -1

  end function checkGSFinfo


  !!============================================================================
  !> @brief Converts the given system matrix into a full rectangular matrix.
  !>
  !> @param[in] this The sysmatrixtypemodule::sysmatrixtype object to convert
  !> @param[out] fullMat Rectangular matrix representation of @a this
  !> @param[in] ndim Max dimension on @a fullMat
  !> @param[in] startRow First row in @a fullMat to insert reactangular matrix
  !> @param[in] ksa Sign of converted matrix
  !> @param[out] ierr Error flag
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 30 Jul 2003

  subroutine convertSysMat (this,fullMat,ndim,startRow,ksa,ierr)

    use ReportErrorModule, only : getErrorFile, reportError, debugFileOnly_p

    type(SysMatrixType), intent(in)  :: this
    real(dp)           , intent(out) :: fullMat(:,:)
    integer            , intent(in)  :: ndim, startRow, ksa
    integer            , intent(out) :: ierr

    !! Local variables
    integer     :: i, i0, j, mdim, nrow
    integer(ik) :: lerr

    !! --- Logic section ---

    ierr = 0
    i0   = startRow-1
    nrow = size(fullMat,1)
    mdim = min(ndim,this%dim,nrow-i0)
    if (abs(ksa) < 10) fullMat = 0.0_dp

    if (.not. associated(this%value)) return

    select case (this%storageType)

    case (diagonalMatrix_p)

       do i = 1, mdim
          fullMat(i0+i,i) = this%value(i)*sign(1,ksa)
       end do

    case (skylineMatrix_p)

       if (.not. associated(this%skyline)) return
       call SKYCONV (this%value(1), this%skyline%msky(1), &
            &        fullMat(startRow,1), nrow, mdim, sign(11,ksa))

    case (sparseMatrix_p)

       if (.not. associated(this%sparse)) return
       call SPRCONV (this%value(1+this%dim), fullMat(startRow,1), &
            &        this%sparse%mspar(1), this%sparse%mtrees(1), &
            &        this%sparse%msifa(1), int(nrow,ik), int(mdim,ik), &
            &        int(sign(11,ksa),ik), int(getErrorFile(),ik), lerr)
       if (lerr < 0) then
          ierr = int(lerr)
          call reportError (debugFileOnly_p,'convertSysMat')
          return
       end if

    case (denseMatrix_p)

       do j = 1, mdim
          do i = 1, mdim
             fullMat(i0+i,j) = this%value(i+this%dim*(j-1))*sign(1,ksa)
          end do
       end do

    case (pardiso_p)

       do i = 1, mdim
          do j = 1, this%pardiso%ia(i), this%pardiso%ia(i+1)-1
             fullMat(i0+i,j) = this%value(this%pardiso%ja(j))*sign(1,ksa)
          end do
       end do

    end select

  end subroutine convertSysMat


  !!============================================================================
  !> @brief Standard routine for writing an object to file.
  !>
  !> @param this The sysmatrixtypemodule::sysmatrixtype object to write
  !> @param[in] mpar Matrix of parameters
  !> @param[in] io File unit number to write to
  !> @param[in] text If present, write as heading
  !> @param[in] nelL Number of matrix elements to write per line
  !> @param[in] complexity If present, the value indicates the amount of print
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 26 Feb 2003

  subroutine WriteSysMat (this,mpar,io,text,nelL,complexity)

    use ManipMatrixModule, only : writeObject
    use FELinearSolver   , only : GSFDump

    type(SysMatrixType),     intent(inout) :: this
    integer                   , intent(in) :: mpar(:), io
    character(len=*), optional, intent(in) :: text
    integer         , optional, intent(in) :: nelL, complexity

    !! Local variables
    integer :: i, n

    !! --- Logic section ---

    if (present(text)) write(io,'(/A)') text

    if (this%storageType < 0) then
       write(io,*) '    type        = (undefined)'
    else
       write(io,*) '    type        = ',matrixType_p(this%storageType)
    end if
    write(io,*) '    dim         =', this%dim
    if (associated(this%value)) then
       write(io,*) '    size        =', size(this%value)
    end if
    if (associated(this%skyline)) then
       write(io,*) '    eigenSolver =', this%skyline%eigenSolver
    end if

    if (present(complexity)) then
       if (complexity > 2 .and. associated(this%meqn)) then
          n = size(this%meqn)
          write(io,*)
          write(io,*) 'MEQN', mpar(3), n
          write(io,'(2I8)') (i,this%meqn(i),i=1,n)
          write(io,*)
       end if
    end if

    select case (this%storageType)

    case (diagonalMatrix_p)

       if (.not. associated(this%value)) return
       call writeObject (this%value,io,nelLin=nelL)

    case (skylineMatrix_p)

       if (.not. associated(this%skyline)) return
       if (.not. associated(this%value))   return
       call writeObject (this%value,this%skyline,io,nelL,complexity)

    case (sparseMatrix_p)

       if (.not. associated(this%sparse)) return
       if (.not. associated(this%value))  return
       call writeObject (this%value,this%sparse,io,nelL,complexity)

    case (denseMatrix_p)

       n = this%dim
       if (.not. associated(this%value)) return
       call writeObject (reshape(this%value,(/n,n/)),io,nelLin=nelL)

    case (outOfCore_p)

       if (associated(this%gsf)) then
          call GSFDump (this%gsf%Q,this%gsf%S,this%gsf%T,io,GSFInfo)
          i = checkGSFinfo(GSFinfo)
       end if

    case (pardiso_p)

       if (.not. associated(this%pardiso)) return
       if (.not. associated(this%value))   return
       call writeObject (this%value,this%pardiso,io,nelL,complexity)

    end select

  end subroutine WriteSysMat


  !!============================================================================
  !> @brief Standard routine for writing an object to file.
  !>
  !> @param this The sysmatrixtypemodule::sysmatrixtype object to write
  !> @param[in] io File unit number to write to
  !> @param[in] text If present, write as heading
  !> @param[in] complexity If present, the value indicates the amount of print
  !>
  !> @details This subroutine only writes the data structure for sparse matrix.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 26 Feb 2003

  subroutine WriteSysMat2 (this,io,text,complexity)

    type(SysMatrixType)       , intent(in) :: this
    integer                   , intent(in) :: io
    character(len=*), optional, intent(in) :: text
    integer         , optional, intent(in) :: complexity

    !! --- Logic section ---

    if (this%storageType == sparseMatrix_p .and. associated(this%sparse)) then
       call writeSparseStructure (this%sparse,io,complexity,text)
    end if

  end subroutine WriteSysMat2


  !!============================================================================
  !> @brief Writes out a skyline matrix, column by column.
  !>
  !> @param[in] matrix The values of the skyline matrix
  !> @param[in] skyline sysmatrixtypemodule::skylinestoragetype object
  !> associated with the matrix
  !> @param[in] io File unit number to write to
  !> @param[in] nelLin Number of lines per matrix to write
  !> @param[in] complexity If present, the value indicates the amount of print
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 27 Aug 2001

  subroutine WriteSkyMat (matrix,skyline,io,nelLin,complexity)

    use ManipMatrixModule , only : writeObject
    use ScratchArrayModule, only : getRealScratchMatrix

    real(dp)                , intent(in) :: matrix(:)
    type(SkylineStorageType), intent(in) :: skyline
    integer                 , intent(in) :: io
    integer, optional       , intent(in) :: nelLin, complexity

    !! Local variables
    integer           :: i, ieq, ip1, ip2, l, nli, nelL
    real(dp), pointer :: fullMat(:,:)

    !! --- Logic section ---

    if (present(nelLin)) then
       nelL = min(20,nelLin) ! Number of matrix elements to write per line
    else
       nelL = 10
    end if

    if (present(complexity)) then
       if (complexity > 3) then
          nelL = size(skyline%msky)
          write(io,*)
          write(io,*) 'MSKY', skyline%neq, nelL
          write(io,'(2I8)') (i,skyline%msky(i),i=1,nelL)
       end if
    end if

    if (skyline%neq <= nelL) then
       fullMat => getRealScratchMatrix(skyline%neq,skyline%neq,i)
       if (i == 0) then
          !! Write out as full matrix
          call SKYCNV (matrix(1),skyline%msky(1),fullMat(1,1),skyline%neq,1)
          call writeObject (fullMat,io,nelLin=nelL)
          return
       end if
    end if

    ip2 = 0
    do ieq = 1, skyline%neq
       ip1 = ip2 + 1
       ip2 = skyline%msky(ieq)
       nli = 1 + (ip2-ip1)/nelL
       write(io,"(/' +++++ Equation',I8,' +++++')") ieq
       do l = 1, nli
          write(io,'(I9,19I13) ') (ieq+i-ip2,i=ip1,min(ip2,ip1+nelL-1))
          write(io,'(1P20E13.5)') (matrix(i),i=ip1,min(ip2,ip1+nelL-1))
          ip1 = ip1 + nelL
       end do
    end do

  end subroutine WriteSkyMat


  !!============================================================================
  !> @brief Writes out a sparse matrix, row by row.
  !>
  !> @param[in] matrix The values of the sparse matrix
  !> @param[in] sparse sysmatrixtypemodule::sparsestoragetype object
  !> associated with the matrix
  !> @param[in] io File unit number to write to
  !> @param[in] nelLin Number of lines per matrix to write
  !> @param[in] complexity If present, the value indicates the amount of print
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 11 Mar 2003

  subroutine WriteSparseMat (matrix,sparse,io,nelLin,complexity)

    use ManipMatrixModule , only : writeObject
    use ScratchArrayModule, only : getRealScratchMatrix

    real(dp)               , intent(in) :: matrix(:)
    type(SparseStorageType), intent(in) :: sparse
    integer                , intent(in) :: io
    integer, optional      , intent(in) :: nelLin, complexity

    !! Local variables
    integer           :: lerr, nelL
    integer(ik)       :: i, n
    real(dp), pointer :: fullMat(:,:)

    !! --- Logic section ---

    if (present(nelLin)) then
       nelL = min(20,nelLin) ! Number of matrix elements to write per line
    else
       nelL = 10
    end if

    if (present(complexity)) then
       if (complexity > 3) then
          n = size(sparse%mspar)
          write(io,*)
          write(io,*) 'MSPAR', n
          write(io,'(2I8)') (i,sparse%mspar(i),i=1,n)
          if (sparse%mspar(58) > 0_ik) then
             n = sparse%mspar(8)
             write(io,*)
             write(io,*) 'PERI2E, PERE2I', n
             write(io,'(3I8)') (i,sparse%msifa(i),sparse%msifa(n+i),i=1_ik,n)
          end if
       end if
    end if

    n = sparse%mspar(8)
    if (n <= int(nelL,ik)) then
       fullMat => getRealScratchMatrix(int(n),int(n),lerr)
       if (lerr == 0) then
          !! Write out as full matrix
          call SPRCNV (matrix(n+1), fullMat(1,1), sparse%mspar(1), &
               &       sparse%mtrees(1), sparse%msifa(1), n, 1_ik, int(io,ik),i)
          if (i == 0_ik) call writeObject (fullMat,io,nelLin=nelL)
          return
       end if
    end if

    call SPRPRN (matrix(n+1), sparse%mspar(1), sparse%mtrees(1), &
         &       sparse%msifa(1), int(io,ik), i)

  end subroutine WriteSparseMat


  !!==========================================================================
  !> @brief Standard routine for writing an object to file.
  !>
  !> @param[in] this The sysmatrixtypemodule::sparsestoragetype object to write
  !> @param[in] io File unit number to write to
  !> @param[in] complexity If present, the value indicates the amount of print
  !> @param[in] label If present, print out as matrix label
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 25 Apr 2003

  subroutine WriteSparseStructure (this,io,complexity,label)

    use ManipMatrixModule, only : writeObject

    type(SparseStorageType), intent(in) :: this
    integer                , intent(in) :: io
    integer     , optional , intent(in) :: complexity
    character(*), optional , intent(in) :: label

    !! Local variables
    integer(ik) :: XELNOD, XNODEL, XNODES, NODES , PERM  , ELNOD, NODEL , SEPARR
    integer(ik) :: PERI2E, PERE2I, LINDX , XLNZ  , XSPLIT, SPLIT
    integer(ik) :: XSUPER, XSIBL , SIBL  , XLINDX, SUPSUP, ELSEQ, XELSEQ, XBLOCK
    integer(ik) :: NODMAP, SUPMAP, IWUSED

    !! --- Logic section ---

    if (present(label)) then
       write(io,'(A)') 'Sparse '//label,'{'
    else
       write(io,'(A)') 'Sparse','{'
    end if
    if (associated(this%mspar))  write(io,*) 'size(mspar)  =', size(this%mspar)
    if (associated(this%msica))  write(io,*) 'size(msica)  =', size(this%msica)
    if (associated(this%msifa))  write(io,*) 'size(msifa)  =', size(this%msifa)
    if (associated(this%mtrees)) write(io,*) 'size(mtrees) =', size(this%mtrees)
    if (associated(this%mvarnc)) write(io,*) 'size(mvarnc) =', size(this%mvarnc)
    if (present(complexity)) then
       if (associated(this%mspar))    call writeObject (this%mspar ,io,'MSPAR')
       if (complexity == 1) then
          if (associated(this%msica)) call writeObject (this%msica ,io,'MSICA')
          if (associated(this%msifa)) call writeObject (this%msifa ,io,'MSIFA')
          if (associated(this%mtrees))call writeObject (this%mtrees,io,'MTREES')
          if (associated(this%mvarnc))call writeObject (this%mvarnc,io,'MVARNC')
       else if (complexity >= 2 .and. associated(this%mspar)) then
          if (associated(this%msica) .and. this%mspar(1) > 0) then
             XELNOD = 1_ik
             XNODEL = XELNOD + this%mspar(5) + 1_ik
             XNODES = XNODEL + this%mspar(6) + 1_ik
             NODES  = XNODES + this%mspar(6) + 1_ik
             PERM   = NODES  + this%mspar(8)
             ELNOD  = PERM   + this%mspar(6)
             NODEL  = ELNOD  + this%mspar(7)
             SEPARR = NODEL  + this%mspar(7)
             IWUSED = SEPARR + this%mspar(6)
             if (complexity == 3) then
                call writeCon (this%msica,XNODES,NODES,this%mspar(6),'NODES')
                call writeCon (this%msica,XELNOD,ELNOD,this%mspar(5),'ELNOD')
             else
                call writeSubArray (this%msica,XNODES,NODES ,'XNODES')
                call writeSubArray (this%msica,NODES ,PERM  ,'NODES')
                call writeSubArray (this%msica,XELNOD,XNODEL,'XELNOD')
                call writeSubArray (this%msica,ELNOD ,NODEL ,'ELNOD')
             end if
             if (this%mspar(1) > 1) then
                if (complexity == 3) then
                   call writeCon (this%msica,XNODEL,NODEL,this%mspar(6),'NODEL')
                else
                   call writeSubArray (this%msica,XNODEL,XNODES,'XNODEL')
                   call writeSubArray (this%msica,NODEL ,SEPARR,'NODEL')
                end if
             end if
             call writeSubArray (this%msica,PERM  ,ELNOD ,'PERM')
             if (complexity == 3) then
                call writeSet (this%msica,SEPARR,IWUSED,'SEPARR')
             else
                call writeSubArray (this%msica,SEPARR,IWUSED,'SEPARR')
             end if
          end if
          if (associated(this%msifa) .and. this%mspar(1) > 3) then
             PERI2E = 1_ik
             PERE2I = PERI2E + this%mspar(8)
             LINDX  = this%mspar(58)*this%mspar(8) + 1_ik
             XLNZ   = LINDX  + this%mspar(15)
             XSPLIT = XLNZ   + this%mspar(11) + 1_ik
             SPLIT  = XSPLIT + this%mspar(11)
             IWUSED = SPLIT  + this%mspar(22) + 1_ik
             if (this%mspar(58) > 0) then
                call writeSubArray (this%msifa,PERI2E,PERE2I,'PERI2E')
                call writeSubArray (this%msifa,PERE2I,LINDX ,'PERE2I')
             end if
             call writeSubArray (this%msifa,LINDX,XLNZ  ,'LINDX')
             call writeSubArray (this%msifa,XLNZ ,XSPLIT,'XLNZ')
             if (this%mspar(22) > 0) then
                call writeSubArray (this%msifa,XSPLIT,SPLIT ,'XSPLIT')
                call writeSubArray (this%msifa,SPLIT ,IWUSED,'SPLIT')
             end if
          end if
          if (associated(this%mtrees) .and. this%mspar(1) > 2_ik) then
             XSUPER = 1
             XSIBL  = XSUPER + this%mspar(11) + 1_ik
             SIBL   = XSIBL  + this%mspar(11) + 1_ik
             XLINDX = SIBL   + this%mspar(11) + 1_ik
             SUPSUP = XLINDX + this%mspar(11) + 1_ik
             ELSEQ  = SUPSUP + this%mspar(11)
             XELSEQ = ELSEQ  + this%mspar(19)
             XBLOCK = XELSEQ + this%mspar(12) + 1_ik
             IWUSED = XBLOCK + this%mspar(12) + 1_ik
             call writeSubArray (this%mtrees,XSUPER,XSIBL ,'XSUPER')
             call writeSubArray (this%mtrees,XLINDX,SUPSUP,'XLINDX')
             if (complexity == 3) then
                call writeCon (this%mtrees,XSIBL ,SIBL ,this%mspar(11),'SIBL')
                call writeSet (this%mtrees,SUPSUP,ELSEQ,'SUPSUP')
                call writeCon (this%mtrees,XELSEQ,ELSEQ,this%mspar(12),'ELSEQ')
             else
                call writeSubArray (this%mtrees,XSIBL ,SIBL  ,'XSIBL')
                call writeSubArray (this%mtrees,SIBL  ,XLINDX,'SIBL')
                call writeSubArray (this%mtrees,SUPSUP,ELSEQ ,'SUPSUP')
                call writeSubArray (this%mtrees,ELSEQ ,XELSEQ,'ELSEQ')
                call writeSubArray (this%mtrees,XELSEQ,XBLOCK,'XELSEQ')
             end if
             call writeSubArray (this%mtrees,XBLOCK,IWUSED,'XBLOCK')
          end if
          if (associated(this%mvarnc) .and. this%mspar(1) > 3_ik) then
             NODMAP = 1_ik
             SUPMAP = NODMAP + this%mspar(8)
             IWUSED = SUPMAP + this%mspar(8)
             call writeSubArray (this%mvarnc,NODMAP,SUPMAP,'NODMAP')
             call writeSubArray (this%mvarnc,SUPMAP,IWUSED,'SUPMAP')
          end if
       end if
    end if
    write(io,'(A)') '}'

  contains

    !> @brief Convenience wrapper writing a sub-matrix of given matrix.
    subroutine writeSubArray (array,from,to,name)
      integer(ik)     , intent(in) :: array(:), from, to
      character(len=*), intent(in) :: name
      if (to > from) call writeObject (array(from:to-1_ik),io,name)
    end subroutine writeSubArray

    !> @brief Convenience wrapper writing out connectivity arrays.
    subroutine writeCon (array,x,y,n,name)
      integer(ik)     , intent(in) :: array(:), x, y, n
      character(len=*), intent(in) :: name
      integer(ik)                  :: i, j, k
      write(io,'(a)') name
      do i = 1, n
         j = y + array(x+i-1) - 1_ik
         k = y + array(x+i)   - 2_ik
         write(io,"(I8,' :',10I8)") i,array(j:min(k,j+9_ik))
         if (j+9_ik < k) write(io,"(10X,10I8)") array(j+10_ik:k)
      end do
    end subroutine writeCon

    !> @brief Convenience wrapper writing out the index sets of non-zero values.
    subroutine writeSet (array,from,to,name)
      integer(ik)     , intent(in) :: array(:), from, to
      character(len=*), intent(in) :: name
      integer(ik)                  :: i, nSet
      write(io,100) name,name
100   format(a,': the values below are indices for which ',a,'(i) > 0')
      nSet = 0_ik
      do i = from, to-1
         if (array(i) > 0_ik) then
            nSet = nSet + 1
            if (mod(nSet,10) > 0_ik) then
               write(io,'(I8)',advance='no') i-from+1_ik
            else
               write(io,'(I8)') i-from+1_ik
            end if
         end if
      end do
      write(io,"('  size =',I6)") nSet
    end subroutine writeSet

  end subroutine WriteSparseStructure


  !!============================================================================
  !> @brief Writes out a sparse matrix, row by row.
  !>
  !> @param[in] matrix The values of the sparse matrix
  !> @param[in] sparse sysmatrixtypemodule::pardisostoragetype object
  !> associated with the matrix
  !> @param[in] io File unit number to write to
  !> @param[in] nelLin Number of lines per matrix to write
  !> @param[in] complexity If present, the value indicates the amount of print
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 10 Feb 2016

  subroutine WritePardisoMat (matrix,sparse,io,nelLin,complexity)

    use ManipMatrixModule , only : writeObject
    use ScratchArrayModule, only : getRealScratchMatrix

    real(dp)                , intent(in) :: matrix(:)
    type(PardisoStorageType), intent(in) :: sparse
    integer                 , intent(in) :: io
    integer, optional       , intent(in) :: nelLin, complexity

    !! Local variables
    integer           :: i, j, n, ieq, ip1, ip2, nli, nelL
    real(dp), pointer :: fullMat(:,:)

    !! --- Logic section ---

    if (present(nelLin)) then
       nelL = min(20,nelLin) ! Number of matrix elements to write per line
    else
       nelL = 10
    end if

    write(io,*) '    mtype       =', sparse%mtype
    write(io,"('     iparm       =',10I6)") sparse%iparm(1:10)
    write(io,"(18X,10I6)") sparse%iparm(11:)
    do i = 1, size(sparse%iparm)
       if (sparse%iparm(i) > 99999 .or. sparse%iparm(i) < -9999) then
          write(io,"('     iparm(',I2,')   =',I12)") i,sparse%iparm(i)
       end if
    end do

    n = size(sparse%ia) - 1
    if (present(complexity)) then
       if (complexity > 3) then
          write(io,*) 'IA', n+1
          write(io,'(2I8)') (i,sparse%ia(i),i=1,n+1)
       end if
    end if

    if (n <= nelL) then
       fullMat => getRealScratchMatrix(n,n,i)
       if (i == 0) then
          !! Write out as full matrix
          fullMat = 0.0_dp
          do i = 1, n
             do j = sparse%ia(i), sparse%ia(i+1)-1
                fullMat(i,sparse%ja(j)) = matrix(j)
                fullMat(sparse%ja(j),i) = matrix(j)
             end do
          end do
          call writeObject (fullMat,io,nelLin=nelL)
          return
       end if
    end if

    do ieq = 1, n
       ip1 = sparse%ia(ieq)
       ip2 = sparse%ia(ieq+1)-1
       nli = 1 + (ip2-ip1)/nelL
       write(io,"(/' +++++ Equation',I8,' +++++')") ieq
       do i = 1, nli
          write(io,'(I9,19I13) ') (sparse%ja(j),j=ip1,min(ip2,ip1+nelL-1))
          write(io,'(1P20E13.5)') (   matrix(j),j=ip1,min(ip2,ip1+nelL-1))
          ip1 = ip1 + nelL
       end do
    end do

  end subroutine WritePardisoMat


  !!============================================================================
  !> @brief Extracts the diagonal elements from the given system matrix.
  !>
  !> @param this The sysmatrixtypemodule::sysmatrixtype object to extract from
  !> @param[out] diag Diagonal elements of the matrix
  !> @param[out] ierr Error flag
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Mar 2004

  subroutine extractDiagonal (this,diag,ierr)

    use FELinearSolver   , only : FEExtractDiagA
    use ReportErrorModule, only : internalError, getErrorFile
    use ReportErrorModule, only : reportError, debugFileOnly_p

    type(SysMatrixType), intent(inout) :: this
    real(dp)           , intent(out)   :: diag(:)
    integer            , intent(out)   :: ierr

    !! Local variables
    integer     :: i
    integer(ik) :: lerr

    !! --- Logic section ---

    if (.not. associated(this%value) .and. this%storageType /= outOfCore_p) then
       ierr = internalError('extractDiagonal: System matrix not allocated')
       return
    else if (this%dim > size(diag)) then
       ierr = internalError('extractDiagonal: Output array is too small')
       return
    end if

    ierr = 0
    select case (this%storageType)

    case (diagonalMatrix_p)

       diag(1:this%dim) = this%value

    case (skylineMatrix_p)

       if (.not. associated(this%skyline)) then
          ierr = internalError('extractDiagonal: skyline')
          return
       end if

       do i = 1, this%skyline%neq
          diag(i) = this%value(this%skyline%msky(i))
       end do

    case (sparseMatrix_p)

       if (.not. associated(this%sparse)) then
          ierr = internalError('extractDiagonal: sparse')
          return
       end if

       call SPREXD (this%sparse%mspar(1), this%sparse%mtrees(1), &
            &       this%sparse%msifa(1), this%value(1), diag(1), &
            &       int(getErrorFile(),ik), lerr)
       ierr = int(lerr)

    case (outOfCore_p)

       if (.not. associated(this%gsf)) then
          ierr = internalError('extractDiagonal: gsf not allocated')
          return
       end if

       call FEExtractDiagA (this%dim, this%gsf%Q, this%gsf%S, diag, GSFinfo)
       ierr = checkGSFinfo (GSFinfo)

    case (denseMatrix_p)

       do i = 1, this%dim
          diag(i) = this%value(i+this%dim*(i-1))
       end do

    case (pardiso_p)

       do i = 1, this%dim
          diag(i) = this%value(this%pardiso%ja(this%pardiso%ia(i)))
       end do

    case default

       ierr = internalError('extractDiagonal: Invalid matrix type')
       return

    end select

    if (ierr < 0) call reportError (debugFileOnly_p,'extractDiagonal')

  end subroutine extractDiagonal


  !!============================================================================
  !> @brief Saves the given system matrix to temporary file.
  !>
  !> @param this The sysmatrixtypemodule::sysmatrixtype object to save
  !> @param[in] name Name tag of matrix
  !> @param[in] nWord Number of real words in the matrix to save
  !> @param[out] iFile Handle of the file the matrix is written to
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 21 Sep 2004

  subroutine saveSysMat (this,name,nWord,iFile,ierr)

    use KindModule       , only : i8
    use ReportErrorModule, only : reportError, error_p
    use BinaryDBInterface, only : openBinaryDB, writeDoubleDB, tmp_p

    type(SysMatrixType), intent(in)  :: this
    character(len=*)   , intent(in)  :: name
    integer(i8)        , intent(in)  :: nWord
    integer            , intent(out) :: iFile, ierr

    !! --- Logic section ---

    call openBinaryDB ('SysMatrix',tmp_p,iFile,ierr)
    if (ierr /= 0) goto 900

    call writeDoubleDB (iFile,this%value(1),nWord,ierr)
    if (ierr >= 0) return

900 call reportError (error_p,'Can not save system '//name//' on file.', &
         'Check available space on the disk drive for temporary files.')

  end subroutine saveSysMat


  !!============================================================================
  !> @brief Restores the given system matrix from temporary file.
  !>
  !> @param this The sysmatrixtypemodule::sysmatrixtype object to restore
  !> @param[in] name Name tag of matrix
  !> @param[in] nWord Number of real words in the matrix to restore
  !> @param[in] iFile File handle to restore from
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 21 Sep 2004

  subroutine restoreSysMat (this,name,nWord,iFile,ierr)

    use KindModule       , only : i8
    use ReportErrorModule, only : reportError, error_p
    use BinaryDBInterface, only : closeBinaryDB, readDoubleDB

    type(SysMatrixType), intent(inout) :: this
    character(len=*)   , intent(in)    :: name
    integer(i8)        , intent(in)    :: nWord
    integer            , intent(inout) :: iFile
    integer            , intent(out)   :: ierr

    !! --- Logic section ---

    call readDoubleDB (-iFile,this%value(1),abs(nWord),ierr)
    if (ierr < 0) then
       call reportError (error_p,'Can not restore system '//name)
    else if (this%storageType == sparseMatrix_p) then
       this%sparse%mspar(1) = 4_ik
    end if

    if (nWord < 0_i8 .and. ierr >= 0) then
       !! Close (and delete) the temporary file
       call closeBinaryDB (iFile,ierr)
       iFile = -1
    end if

  end subroutine restoreSysMat


  !!============================================================================
  !> @brief Adds matrix elements corresponding to free and constrained DOFs.
  !>
  !> @param[in] EM The element matrix to add
  !> @param[in] MEEN Matrix of element equation numbers
  !> @param[in] MEQN Matrix of equation numbers
  !> @param[in] MPMCEQ Matrix of pointers to constraint equations
  !> @param[in] MMCEQ Matrix of constraint equations
  !> @param[in] TTCC Table of constraint equation coefficients
  !> @param A System matrix coefficients
  !> @param[in] IA Index to to the start of each column in @a A
  !> @param[in] JA Row indices of @a A
  !> @param B Right-hand-side vector
  !> @param[out] IERR Error flag
  !>
  !> @details The elements corresponding to free and constrained DOFs of the
  !> element matrix @a EM is added to the sparse system matrix @a {A,IA,JA}
  !> and the associated right-hand-side vector @a B.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 12 Feb 2016

  subroutine asmSparse (EM,MEEN,MEQN,MPMCEQ,MMCEQ,TTCC,A,IA,JA,B,ierr)

    use reportErrorModule, only : internalError

    real(dp)          , intent(in)    :: EM(:,:), TTCC(:)
    integer           , intent(in)    :: MEEN(:), MEQN(:), MPMCEQ(:), MMCEQ(:)
    integer           , intent(in)    :: IA(:), JA(:)
    real(dp)          , intent(inout) :: A(:)
    real(dp), optional, intent(inout) :: B(:)
    integer           , intent(out)   :: IERR

    !! Local variables
    real(dp), parameter :: eps_p = 1.0e-15_dp

    integer  :: i, j, j1, j2, ip, jp, ieq, jeq, iceq, jceq, ndof, ncol
    real(dp) :: c

    !! --- Logic section ---

    ndof = size(meen)
    ncol = size(EM,2)
    ierr = size(EM,1) - ndof
    if (ierr < 0 .or. (ncol /= 1 .and. ncol < ndof)) then
       ierr = internalError('asmSparse: Too small element matrix')
       return
    end if

    !! Add elements corresponding to free dofs in EM into A

    do i = 1, ndof
       ieq = meen(i)
       if (ieq < 1) cycle

       !! Add contributions to A
       j1 = IA(ieq)
       j2 = IA(ieq+1)-1
       if (ncol > 1) then
          do j = 1, ndof
             jeq = meen(j)
             if (jeq >= ieq) then ! Upper triangle only (assuming symmetry)
                call addToRow (jeq,j1,j2,EM(i,j),ierr)
                if (ierr < 0) return
             end if
          end do
       else ! Diagonal element matrix
          call addToRow (ieq,j1,j2,EM(i,1),ierr)
          if (ierr < 0) return
       end if

    end do

    if (all(meen >= 0)) return ! No prescribed DOFs in this element

    !! Add (appropriately weighted) elements corresponding to constrained
    !! (dependent and prescribed) dofs in EM into A and (optionally) B

    do i = 1, ndof
       iceq = -meen(i)
       if (iceq < 1) cycle

       c = ttcc(mpmceq(iceq))
       if (present(B) .and. abs(c) > eps_p) then

          !! Add contributions to B (right-hand-side)
          if (ncol > 1) then
             do j = 1, ndof
                jeq = meen(j)
                if (jeq > 0) then
                   B(jeq) = B(jeq) - c*EM(i,j)
                else if (jeq < 0) then
                   jceq = -jeq
                   do jp = mpmceq(jceq)+1, mpmceq(jceq+1)-1
                      if (mmceq(jp) > 0) then
                         jeq = meqn(mmceq(jp))
                         B(jeq) = B(jeq) - c*ttcc(jp)*EM(i,j)
                      end if
                   end do
                end if
             end do
          else ! Diagonal element matrix
             do jp = mpmceq(iceq)+1, mpmceq(iceq+1)-1
                if (mmceq(jp) > 0) then
                   jeq = meqn(mmceq(jp))
                   B(jeq) = B(jeq) - c*ttcc(jp)*EM(i,1)
                end if
             end do
          end if

       end if

       !! Add contributions to A
       do ip = mpmceq(iceq)+1, mpmceq(iceq+1)-1
          c  = ttcc(ip)
          if (mmceq(ip) > 0 .and. abs(c) > eps_p) then
             ieq = meqn(mmceq(ip))
             j1 = IA(ieq)
             j2 = IA(ieq+1)-1
             if (ncol > 1) then
                do j = 1, ndof
                   jeq = meen(j)
                   if (jeq >= ieq) then ! Upper triangle only
                      call addToRow (jeq,j1,j2,c*EM(i,j),ierr)
                      if (ierr < 0) return
                   else if (jeq > 0) then
                      call addToRow (ieq,IA(jeq),IA(jeq+1)-1,c*EM(i,j),ierr)
                      if (ierr < 0) return
                   else if (jeq < 0) then
                      jceq = -jeq
                      do jp = mpmceq(jceq)+1, mpmceq(jceq+1)-1
                         if (mmceq(jp) > 0) then
                            jeq = meqn(mmceq(jp))
                            if (jeq >= ieq) then ! Upper triangle only
                               call addToRow (jeq,j1,j2,c*ttcc(jp)*EM(i,j),ierr)
                               if (ierr < 0) return
                            end if
                         end if
                      end do
                   end if
                end do
             else ! Diagonal element matrix
                do jp = mpmceq(iceq)+1, mpmceq(iceq+1)-1
                   if (mmceq(jp) > 0) then
                      jeq = meqn(mmceq(jp))
                      if (jeq >= ieq) then ! Upper triangle only
                         call addToRow (jeq,j1,j2,c*ttcc(jp)*EM(i,1),ierr)
                         if (ierr < 0) return
                      end if
                   end if
                end do
             end if
          end if
       end do

    end do

  contains

    !> @brief Adds one row into the system matrix.
    subroutine addToRow (jeq,j1,j2,value,ierr)
      integer , intent(in)  :: jeq, j1, j2
      real(dp), intent(in)  :: value
      integer , intent(out) :: ierr
      integer :: jmid, jlow, jhigh
      jlow  = j1
      jhigh = j2
      do while (jlow <= jhigh)
         jmid = (jhigh+jlow)/2
         if (JA(jmid) < jeq) then
            jlow = jmid+1
         else if (JA(jmid) > jeq) then
            jhigh = jmid-1
         else
            A(jmid) = A(jmid) + value
            ierr = 0
            return
         end if
      end do
      ierr = internalError('asmSparse: Invalid sparse matrix structure')
    end subroutine addToRow

  end subroutine asmSparse

end module SysMatrixTypeModule
