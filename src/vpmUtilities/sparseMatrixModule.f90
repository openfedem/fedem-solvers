!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file sparseMatrixModule.f90
!> @brief General sparse matrix representation.

!!==============================================================================
!> @brief Module for management of general sparse matrices.
!> @details This module contains a data type and associated utility routines
!> for representing a general sparse matrix stored on classic coordinate form.
!> @note All data members of sparsematrixmodule::sparsematrixtype are private
!> and can therefore be accessed only through the subroutines/functions
!> contained in this module.

module SparseMatrixModule

  use KindModule, only : dp

  implicit none

  !> Initial allocated size for a new sparse matrix.
  integer , parameter, private :: initSize_p   = 10000
  !> Relative size increment when the matrix is too small
  real(dp), parameter, private :: growFactor_p = 1.5_dp

  !> @brief Data type representing a rectangular sparse matrix.
  type SparseMatrixType
     private
     integer           :: nNonZero !< Number of non-zero elements
     integer           :: nRows    !< Number of matrix rows
     integer           :: nCols    !< Number of matrix columns
     real(dp), pointer :: vala(:)  !< Values of non-zero elements
     integer , pointer :: iRow(:)  !< Row indices of the non-zero entries
     integer , pointer :: jCol(:)  !< Column indices of the non-zero entries
  end type SparseMatrixType

  private :: smGrow


contains

  !!============================================================================
  !> @brief Initializes a sparse matrix object.
  !>
  !> @param[out] this The sparsematrixmodule::sparsematrixtype object
  !>                  to be initialized
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date Mar 2003

  subroutine smNullify (this)

    type(SparseMatrixType), intent(out) :: this

    !! --- Logic section ---

    this%nRows    = 0
    this%nCols    = 0
    this%nNonZero = 0

    nullify(this%vala)
    nullify(this%iRow)
    nullify(this%jCol)

  end subroutine smNullify


  !!============================================================================
  !> @brief Allocates a sparse matrix object.
  !>
  !> @param[out] this The sparsematrixmodule::sparsematrixtype object
  !>                  to be allocated
  !> @param[in] nRows Number of matrix rows
  !> @param[in] nCols Number of matrix columns
  !> @param[in] useInitSize If .true., perform initial allocation using the
  !>            sparsematrixmodule::initsize_p value as number of non-zeroes
  !> @param[out] err Error indicator
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date Jan 2003

  subroutine smAllocate (this,nRows,nCols,useInitSize,err)

    use allocationModule, only : reAllocate

    type(SparseMatrixType), intent(out) :: this
    integer               , intent(in)  :: nRows, nCols
    logical               , intent(in)  :: useInitSize
    integer               , intent(out) :: err

    !! --- Logic section ---

    err = 0
    call smNullify (this)

    this%nRows = nRows
    this%nCols = nCols

    if (useInitSize) then
       call reAllocate ('smAllocate',this%vala,initSize_p,err)
       call reAllocate ('smAllocate',this%iRow,initSize_p,err)
       call reAllocate ('smAllocate',this%jCol,initSize_p,err)
    end if

  end subroutine smAllocate


  !!============================================================================
  !> @brief Initializes a sparse matrix with the given data.
  !>
  !> @param this The sparsematrixmodule::sparsematrixtype object to initialize
  !> @param[in] mRow Row indices of the non-zero entries
  !> @param[in] mCol Column indices of the non-zero entries
  !> @param[in] values All non-zero entries of the matrix
  !> @param[out] err Error indicator
  !>
  !> @note The input arrays are referenced directly by the
  !> sparsematrixmodule::sparsematrixtype object (pointer-associated)
  !> and not copied. They must therefore not be accessed from outside the object
  !> after its destruction (by calling sparsematrixmodule::smdeallocate),
  !> as these arrays then would be deallocated as well.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 16 Sep 2004

  subroutine smSetMatrix (this,mRow,mCol,values,err)

    use KindModule      , only : nbi_p, nbd_p
    use allocationModule, only : doLogMem, logAllocMem, reAllocate

    type(SparseMatrixType), intent(inout) :: this
    integer               , pointer       :: mRow(:), mCol(:)
    real(dp)              , pointer       :: values(:)
    integer               , intent(out)   :: err

    !! --- Logic section ---

    err = 0
    if (.not.associated(values)) err = err - 1
    if (.not.associated(mRow))   err = err - 1
    if (.not.associated(mCol))   err = err - 1
    if (err < 0) return

    this%nNonZero = size(values)
    if (this%nNonZero /= size(mRow)) err = err - 1
    if (this%nNonZero /= size(mCol)) err = err - 1
    if (err < 0) return

    call reAllocate ('smSetMatrix',this%vala)
    call reAllocate ('smSetMatrix',this%iRow)
    call reAllocate ('smSetMatrix',this%jCol)
    this%vala => values
    this%iRow => mRow
    this%jCol => mCol

    if (doLogMem) then
       call logAllocMem ('smSetMatrix',0,size(this%vala),nbd_p)
       call logAllocMem ('smSetMatrix',0,size(this%iRow)+size(this%jCol),nbi_p)
    end if

  end subroutine smSetMatrix


  !!============================================================================
  !> @brief Returns the dimension of a sparse matrix.
  !>
  !> @param[in] this The sparsematrixmodule::sparsematrixtype object to return
  !>                 the dimension of
  !> @param[in] iDim Which dimension to return the size for:
  !>            1. Number of rows
  !>            2. Number of columns
  !>            3. Number of non-zeroes
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date Jan 1999

  function smSize (this,iDim)

    type(SparseMatrixType), intent(in) :: this
    integer, optional     , intent(in) :: iDim
    integer                            :: smSize

    !! --- Logic section ---

    if (.not. present(idim)) then
       smSize = this%nRows*this%nCols
    else if (iDim == 1) then
       smSize = this%nRows
    else if (iDim == 2) then
       smSize = this%nCols
    else
       smSize = this%nNonZero
    end if

  end function smSize


  !!============================================================================
  !> @brief Deallocates a sparse matrix.
  !>
  !> @param this The sparsematrixmodule::sparsematrixtype object to deallocate
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date Jan 2003

  subroutine smDeallocate (this)

    use allocationModule, only : reAllocate

    type(SparseMatrixType), intent(inout) :: this

    !! --- Logic section ---

    call reAllocate ('smDeallocate',this%vala)
    call reAllocate ('smDeallocate',this%iRow)
    call reAllocate ('smDeallocate',this%jCol)
    call smNullify (this)

  end subroutine smDeallocate


  !!============================================================================
  !> @brief Transposes a sparse matrix.
  !>
  !> @param this The sparsematrixmodule::sparsematrixtype object to transpose
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date Jan 1999

  subroutine smTranspose (this)

    type(SparseMatrixType), intent(inout) :: this

    !! Local variables
    integer :: tmp, i

    !! --- Logic section ---

    tmp        = this%nRows
    this%nRows = this%nCols
    this%nCols = tmp

    do i = 1, this%nNonZero
       tmp          = this%iRow(i)
       this%iRow(i) = this%jCol(i)
       this%jCol(i) = tmp
    end do

  end subroutine smTranspose


  !!============================================================================
  !> @brief Assigns a matrix element.
  !>
  !> @param this The sparsematrixmodule::sparsematrixtype object
  !>             to receive the new matrix element
  !> @param[in] value The new matrix element value
  !> @param[in] rInd Row index of the new matrix element
  !> @param[in] cInd Column index of the new matrix element
  !> @param[out] err Error indicator
  !> @param[in] lpu Optional file unit number for error- and other messages
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date Jan 1999

  subroutine smSetValue (this,value,rInd,cInd,err,lpu)

    type(SparseMatrixType), intent(inout) :: this
    real(dp)              , intent(in)    :: value
    integer               , intent(in)    :: rInd, cInd
    integer               , intent(out)   :: err
    integer, optional     , intent(in)    :: lpu

    !! Local variables
    integer :: i

    real(dp), parameter :: epsValue_p = 1.0e-30_dp

    !! --- Logic section ---

    err = 0
    if (abs(value) < epsValue_p) return

    !! Check if the indices already exist, and set the associated value if so
    do i = 1, this%nNonZero
       if (this%iRow(i) == rInd .and. this%jCol(i) == cInd) then
          this%vala(i) = value
          return
       end if
    end do

    !! New indices, increment number of nonZeroes and insert new value
    this%nNonZero = this%nNonZero + 1
    if (this%nNonZero > size(this%vala)) then
       call smGrow (this,err,lpu)
       if (err < 0) return
    end if

    i = this%nNonZero
    this%vala(i) = value
    this%iRow(i) = rInd
    this%jCol(i) = cInd

  end subroutine smSetValue


  !!============================================================================
  !> @brief Returns a matrix element.
  !>
  !> @param[in] this The sparsematrixmodule::sparsematrixtype object
  !>                 to return the element for
  !> @param[in] rInd Row index of the matrix element to return
  !> @param[in] cInd Column index of the matrix element to return
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date Jan 1999

  function smGetValue (this,rInd,cInd)

    type(SparseMatrixType), intent(in) :: this
    integer               , intent(in) :: rInd, cInd
    real(dp)                           :: smGetValue

    !! Local variables
    integer :: i

    !! --- Logic section ---

    smGetValue = 0.0_dp

    !! Check if indices represent a nonZero value
    do i = 1, this%nNonZero
       if (this%iRow(i) == rInd .and. this%jCol(i) == cInd) then
          smGetValue = this%vala(i)
          return
       end if
    end do

  end function smGetValue


  !!============================================================================
  !> @brief Returns a matrix column.
  !>
  !> @param[in] this The sparsematrixmodule::sparsematrixtype object
  !>                 to return the column for
  !> @param[in] cInd Index of the matrix column to return
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date Jan 1999

  function smGetCol (this,cInd)

    type(SparseMatrixType), intent(in) :: this
    integer               , intent(in) :: cInd
    real(dp)                           :: smGetCol(this%nRows)

    !! Local variables
    integer :: i

    !! --- Logic section ---

    smGetCol = 0.0_dp

    !! Insert the nonZero elements in the column vector
    do i = 1, this%nNonZero
       if (this%jCol(i) == cInd) smGetCol(this%iRow(i)) = this%vala(i)
    end do

  end function smGetCol


  !!============================================================================
  !> @brief Returns a matrix row.
  !>
  !> @param[in] this The sparsematrixmodule::sparsematrixtype object
  !>                 to return the row for
  !> @param[in] rInd Index of the matrix row to return
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date Jan 1999

  function smGetRow (this,rInd)

    type(SparseMatrixType), intent(in) :: this
    integer               , intent(in) :: rInd
    real(dp)                           :: smGetRow(this%nCols)

    !! Local variables
    integer :: i

    !! --- Logic section ---

    smGetRow = 0.0_dp

    !! Insert the nonZero elements in the row vector
    do i = 1, this%nNonZero
       if (this%iRow(i) == rInd) smGetRow(this%jCol(i)) = this%vala(i)
    end do

  end function smGetRow


  !!============================================================================
  !> @brief Multiplies a sparse matrix with a vector.
  !>
  !> @param[in] this The sparsematrixmodule::sparsematrixtype object (@b A)
  !>                 being the left-hand-size of the multiplication operation
  !> @param[in] xVec The right-hand-side vector (@b x) of the multiplication
  !> @param[out] yVec The resulting vector (@b y)
  !> @param[out] err Error indicator
  !>
  !> @details Performs the operation @b y = @b A * @b x .
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date Jan 1999

  subroutine smMatTimesVec (this,xVec,yVec,err)

    use reportErrorModule, only : internalError

    type(SparseMatrixType), intent(in)  :: this
    real(dp)              , intent(in)  :: xVec(:)
    real(dp)              , intent(out) :: yVec(:)
    integer               , intent(out) :: err

    !! Local variables
    integer :: i, j, k

    !! --- Logic section ---

    if (size(xVec) == this%nCols .and. size(yVec) == this%nRows) then
       err = 0
    else
       err = internalError('smMatTimesVec: Dimension mis-match')
       return
    end if

    yVec = 0.0_dp

    !! Perform the multiplication
    do k = 1, this%nNonZero
       i = this%iRow(k)
       j = this%jCol(k)
       yVec(i) = yVec(i) + this%vala(k)*xVec(j)
    end do

  end subroutine smMatTimesVec


  !!============================================================================
  !> @brief Multiplies the transpose of a sparse matrix with a vector.
  !>
  !> @param[in] this The sparsematrixmodule::sparsematrixtype object (@b A)
  !>                 being the left-hand-size of the multiplication operation
  !> @param[in] xVec The right-hand-side vector (@b x) of the multiplication
  !> @param[out] yVec The resulting vector (@b y)
  !> @param[out] err Error indicator
  !>
  !> @details Performs the operation @b y = <b>A</b><sup>t</sup> * @b x .
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date Jan 1999

  subroutine smMatTransTimesVec (this,xVec,yVec,err)

    use reportErrorModule, only : internalError

    type(SparseMatrixType), intent(in)  :: this
    real(dp)              , intent(in)  :: xVec(:)
    real(dp)              , intent(out) :: yVec(:)
    integer               , intent(out) :: err

    !! Local variables
    integer :: i, j, k

    !! --- Logic section ---

    if (size(xVec) == this%nRows .and. size(yVec) == this%nCols) then
       err = 0
    else
       err = internalError('smMatTransTimesVec: Dimension mis-match')
       return
    end if

    yVec = 0.0_dp

    !! Perform the multiplication
    do k = 1, this%nNonZero
       i = this%jCol(k)
       j = this%iRow(k)
       yVec(i) = yVec(i) + this%vala(k)*xVec(j)
    end do

  end subroutine smMatTransTimesVec


  !!============================================================================
  !> @brief Reallocates a sparse matrix to exactly fit its dimensions.
  !>
  !> @param this The sparsematrixmodule::sparsematrixtype object to reallocate
  !> @param[out] err Error indicator
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date Jan 2003

  subroutine smShrinkToFit (this,err)

    use allocationModule, only : reAllocate

    type(SparseMatrixType), intent(inout) :: this
    integer               , intent(out)   :: err

    !! Local variables
    integer            :: newSize, oldSize
    logical, parameter :: preserveContent_p = .true.

    !! --- Logic section ---

    err     = 0
    oldSize = size(this%vala)
    newSize = this%nNonZero
    call reAllocate ('smShrinkToFit',this%vala,newSize,err,preserveContent_p)
    call reAllocate ('smShrinkToFit',this%iRow,newSize,err,preserveContent_p)
    call reAllocate ('smShrinkToFit',this%jCol,newSize,err,preserveContent_p)

  end subroutine smShrinkToFit


  !!============================================================================
  !> @brief Writes out the size parameters for a sparse matrix.
  !>
  !> @param[in] this The sparsematrixmodule::sparsematrixtype object
  !>                  to write size parameters for
  !> @param[in] name Matrix identifier
  !> @param[in] lpu File unit number to write to
  !> @param[in] writeMatrixElms If .true., write out the matrix elements as well
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date Mar 2003

  subroutine smWrite (this,name,lpu,writeMatrixElms)

    use KindModule         , only : nbd_p, nbi_p, i8
    use AllocationModule   , only : StrBytes
    use ScratchArrayModule , only : getIntegerScratchArray
    use SearchAndSortModule, only : quickSortInt

    type(SparseMatrixType), intent(in) :: this
    character(len=*)      , intent(in) :: name
    integer               , intent(in) :: lpu
    logical, optional     , intent(in) :: writeMatrixElms

    !! Local variables
    integer, parameter :: nelLin = 10
    integer            :: i, j, lastRow, lenLine, nnz
    integer, pointer   :: indices(:)
    character(len=16)  :: cBytes

    !! --- Logic section ---

    write(lpu,600) name, this%nRows, this%nCols
    if (this%nNonZero > 0) then
       cBytes = StrBytes(int(this%nNonZero,i8)*int(nbd_p+2*nbi_p,i8))
       write(lpu,610) this%nNonZero,trim(adjustl(cBytes))
    else
       write(lpu,611)
       return
    end if
    if (this%nRows*this%nCols > 0) then
       write(lpu,620) real(this%nNonZero)/real(this%nRows*this%nCols)
    end if

    if (.not. present(writeMatrixElms)) return
    if (.not. writeMatrixElms) return

    !! Sort the non-zeroes in the printing order (row by row)
    nnz = this%nNonZero
    indices => getIntegerScratchArray(2*nnz,i)
    if (i < 0) return
    do i = 1, nnz
       indices(nnz+i) = this%iRow(i)*this%nCols + this%jCol(i)
    end do
    call quickSortInt (indices(nnz+1:2*nnz),indices(1:nnz))

    !! Write out the matrix elements, up to nelLin elements per row
    lastRow = 0
    lenLine = 0
    do i = 1, nnz
       if (this%iRow(indices(i)) /= lastRow) then
          lenLine = 0
          lastRow = this%iRow(indices(i))
          write(lpu,650) lastRow
       else if (lenLine == nelLin) then
          lenLine = 0
          write(lpu,*)
       end if
       if (lenLine == 0) then
          j = i ! Write the column numbers for current row
          do while (this%iRow(indices(j)) == lastRow .and. lenLine < nelLin)
             lenLine = lenLine + 1
             write(lpu,660,ADVANCE='no') this%jCol(indices(j))
             j = j + 1
             if (j > nnz) exit
          end do
          lenLine = 0
          write(lpu,*)
       end if
       lenLine = lenLine + 1
       write(lpu,670,ADVANCE='no') this%vala(indices(i))
    end do
    write(lpu,*)

600 format(/4X,'Storage requirement for sparse matrix ', A, &
         & /7X,'nRows    = ',I10,'  nCols   = ',I10 )
610 format( 7X,'nNonZero = ',I10,'  (',A,')')
611 format( 7X,'There are no non-zero elements in this matrix')
620 format( 7X,'Relative density = ',1PE12.5 )
650 format(/' +++++ Row',I8,' +++++')
660 format(I9,4X)
670 format(1PE13.5)

  end subroutine smWrite


  !!============================================================================
  !> @brief Expands the allocated space for a sparse matrix by a fixed amount.
  !>
  !> @param this The sparsematrixmodule::sparsematrixtype object to expand
  !> @param[out] err Error indicator
  !> @param[in] lpu Optional file unit number logging the rellocation
  !>
  !> @details The number of times @a n to call this subroutine to grow a matrix
  !> to size @a x with a given @a growfactor is:
  !>
  !>     n = log(x) / log(growfactor)
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date Jan 2003

  subroutine smGrow (this,err,lpu)

    use allocationModule, only : reAllocate

    type(SparseMatrixType), intent(inout) :: this
    integer               , intent(out)   :: err
    integer, optional     , intent(in)    :: lpu

    !! Local variables
    integer            :: newSize, oldSize
    logical, parameter :: preserveContent_p = .true.

    !! --- Logic section ---

    err     = 0
    oldSize = size(this%vala)
    newSize = max(int(oldSize*growFactor_p),initSize_p)
    if (present(lpu)) write(lpu,*) 'SparseMatrixType::smGrow new size :',newSize
    call reAllocate ('smGrow',this%vala,newSize,err,preserveContent_p)
    call reAllocate ('smGrow',this%iRow,newSize,err,preserveContent_p)
    call reAllocate ('smGrow',this%jCol,newSize,err,preserveContent_p)

  end subroutine smGrow

end module SparseMatrixModule
