!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file sDiskMatrixModule.f90
!> @brief General rectangular matrix with disk storage.

!!==============================================================================
!> @brief Module for general rectangular matrices with disk storage,
!> @details This module contains a data type and associated utility routines
!> for representing a general rectangular matrix stored on a disk file.
!> This module is single-precision limited version of @ref diskmatrixmodule.
!> It supports ReadOnly matrices only.

module sDiskMatrixModule

  use KindModule, only : sp, dp, i8

  implicit none

  !> @brief Data type representing a rectangular matrix with disk storage.
  type DiskMatrixType
     integer             :: nRows      !< Number of rows in the entire matrix
     integer             :: nCols      !< Number of columns in the entire matrix
     integer             :: nColSwap   !< Number of columns in one swap section
     integer, private    :: nElSwap    !< Number of elements in one swap section
     integer, private    :: nElastSwap !< Number of elements in the last section
     integer, private    :: nSwap      !< Number of swap sections
     integer, private    :: curSwap    !< Current swap section in core
     integer, private    :: swapFile   !< Swap file number
     integer, private    :: swapCnt(2) !< Number of swap operations performed
     integer(i8),pointer :: fileKey(:) !< File address for all swap sections
     real(sp)   ,pointer :: vala(:,:)  !< Matrix elements for one swap section
     real(dp)   ,pointer :: vald(:)    !< Double precision buffer
  end type DiskMatrixType

  !> @brief Multiplies a disk matrix with a given vector.
  interface dmMatTimesVec
     module procedure dmMatTimesVec_SP
     module procedure dmMatTimesVec_DP
  end interface

  private :: dmSwap, dmGetAddress, dmSingleCast
  private :: dmMatTimesVec_SP, dmMatTimesVec_DP


contains

  !!============================================================================
  !> @brief Returns requested size of one swap section from command-line option.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 31 Jan 2002

  function dmGetSwapSize () result(nw)

    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getint

    !! Local variables
    integer            :: nw, swapMB, swapW
    integer, parameter :: MByte = 262144 ! Number of sp words in one megabyte

    !! --- Logic section ---

    call ffa_cmdlinearg_getint ('Bramsize',swapMB) ! Swap size in MegaBytes
    call ffa_cmdlinearg_getint ('dmramsize',swapW) ! Swap size in sp Words
    if (int(swapMB,i8)*int(MByte,i8) > int(huge(1),i8)) then
       nw = huge(1)        ! = 2 GigaWords
    else if (swapMB >= 0) then
       nw = swapMB * MByte ! = swapMB * 1024^2 / 4
    else
       nw = swapW
    end if

  end function dmGetSwapSize


  !!============================================================================
  !> @brief Initializes a disk matrix object.
  !>
  !> @param[out] this The sdiskmatrixmodule::diskmatrixtype object
  !>                  to be initialized
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 14 Feb 2002

  subroutine dmNullify (this)

    type(DiskMatrixType), intent(out) :: this

    !! --- Logic section ---

    this%nRows      =  0
    this%nCols      =  0
    this%nColSwap   =  0
    this%nElSwap    =  0
    this%nElastSwap =  0
    this%nSwap      =  0
    this%curSwap    =  0
    this%swapFile   = -1
    this%swapCnt    =  0

    nullify(this%fileKey)
    nullify(this%vala)
    nullify(this%vald)

  end subroutine dmNullify


  !!============================================================================
  !> @brief Writes out the content of a disk matrix.
  !>
  !> @param[in] this The sdiskmatrixmodule::diskmatrixtype object to write out
  !> @param[in] lpu File unit number to write to
  !> @param[in] name Matrix identifier
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 18 Mar 2003

  subroutine dmWrite (this,lpu,name)

    use manipMatrixModule, only : writeObject

    type(DiskMatrixType), intent(inout) :: this
    integer             , intent(in)    :: lpu
    character(len=*)    , intent(in)    :: name

    !! Local variables
    integer , parameter :: nelLin = 10
    integer             :: i, ierr
    real(sp), pointer   :: column(:)
    character(len=16)   :: chCol

    !! --- Logic section ---

    if (.not.associated(this%vala) .or. this%nSwap < 1) return

#if defined(win32) || defined(win64)
    if (this%nSwap == 1 .and. this%nelSwap < 64000) then
#else
    if (this%nSwap == 1) then
#endif

       !! Write the whole matrix in one go
       call writeObject (this%vala,lpu,name,nelLin)

    else

       !! Write out the matrix column by column
       do i = 1, this%nCols
          column => dmGetColPtr(this,i,ierr)
          if (ierr < 0) exit
          write(chCol,"('  column',I8)") i
          call writeObject (column,lpu,name//chCol,nelLin)
       end do

    end if

  end subroutine dmWrite


  !!============================================================================
  !> @brief Opens a disk matrix file for ReadOnly operations.
  !>
  !> @param[out] this The sdiskmatrixmodule::diskmatrixtype object to open
  !> @param[in] fileName Name of file to read/write matrix content from/to
  !> @param fileChkSum Check-sum value of the FE part associated with @a this
  !> @param[in] nRows Number of matrix rows
  !> @param[in] nCols Number of matrix columns
  !> @param[in] status File opening mode, must be diskmatrixmodule::dmold_p
  !> @param[out] err Error indicator
  !> @param[in] nColSwap Optional number of matrix columns per swap section
  !> @param[in] swapSize Optional size of swap section in double precision words
  !> @param[in] wantTag Expected file tag for ReadOnly matrices
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 1 Jun 2001

  subroutine dmOpen (this,fileName,fileChkSum,nRows,nCols,status,err, &
       &             nColSwap,swapSize,wantTag)

    use kindModule       , only : nbd_p, nbs_p, nbi_p, maxInt_p
    use allocationModule , only : reAllocate, StrBytes
    use reportErrorModule, only : getErrorFile, reportError, error_p
    use binaryDBInterface, only : openBinaryDB, closeBinaryDB, read_p
    use binaryDBInterface, only : readTagDB, getPositionDB
    use binaryDBInterface, only : readFloatDB, readDoubleDB

    type(DiskMatrixType), intent(out)   :: this
    character(len=*)    , intent(in)    :: fileName
    integer             , intent(inout) :: fileChkSum
    integer             , intent(in)    :: nRows, nCols, status
    integer             , intent(out)   :: err
    integer, optional   , intent(in)    :: nColSwap
    integer, optional   , intent(in)    :: swapSize
    character(len=*), optional, intent(in) :: wantTag

    !! Local variables
    logical           :: storeAsDouble
    integer,parameter :: maxDoubleBuf_p = 2048
    integer           :: i, j, lpu, nElSwap, nDouble, nRead
    integer(i8)       :: nBytes
    character(len=32) :: fileTag
    character(len=64) :: errMsg

    character(len=18), parameter :: dmTag_p = '#FEDEM disk matrix'

    !! --- Logic section ---

    err = 0
    if (status > 0) then
       err = -99
       call reportError (error_p,'Invalid use, for ReadOnly operations only', &
            &            addString='dmOpen')
       return
    else if (associated(this%vala)) then
       return
    end if

    this%nRows    = nRows
    this%nCols    = nCols
    this%nColSwap = nCols ! Default: store full matrix in memory
    storeAsDouble = .false.

    !! Set number of columns in swap section
    if (present(nColSwap)) then
       if (nColSwap > 0) then
          this%nColSwap = adjustSwapSize(nCols,nColSwap)
       end if
    else if (present(swapSize) .and. nRows > 0 .and. nCols > 0) then
       if (swapSize > 0) then
          this%nColSwap = adjustSwapSize(nCols,swapSize/nRows)
       end if
    end if

    !! Check that the swap section size is within 2 Giga words
    nBytes = int(nRows,i8)*int(this%nColSwap,i8)
    if (nBytes > maxInt_p) then
       this%nColSwap = adjustSwapSize(nCols,maxInt_p/nRows)
    end if

    !! Set number of swap sections and number of elements in swap section
    this%nSwap      = 1 + (nCols-1)/max(1,this%nColSwap)
    this%nElSwap    = nRows * this%nColSwap
    this%nElastSwap = int(int(nRows,i8)*int(nCols,i8) - &
         &                int(this%nElSwap,i8)*int(this%nSwap-1,i8))

    !! Report matrix size
    lpu = getErrorFile()
    write(lpu,600) trim(fileName), nRows, nCols, this%nSwap, this%nElSwap
    if (this%nElastSwap < this%nElSwap) write(lpu,601) this%nElastSwap

    !! Open swap file and check the file tag
    call openBinaryDB (fileName,read_p,this%swapFile,err)
    if (err == 0) call readTagDB (this%swapFile,fileTag,fileChkSum,err)

    if (err < 0) then
       call reportError (error_p,'Cannot open disk matrix file '//fileName, &
            &            addString='dmOpen')
       goto 110
    end if

    !! Check whether the matrix was saved in single precision
    i = len_trim(fileTag)
    if (fileTag(i-2:i) == ' SP') then
       i = i - 3
    else
       storeAsDouble = .true.
    end if
    if (present(wantTag)) then
       if (fileTag(1:i) /= wantTag) err = -1
    else
       if (fileTag(1:i) /= dmTag_p) err = -1
    end if
    if (err < 0) then
       call reportError (error_p,'Invalid disk matrix file '//fileName, &
            &            'File tag: '//fileTag,addString='dmOpen')
       goto 100
    end if

    !! Report allocated storage in core
#ifdef FT_HAS_INT4ONLY
    nBytes = 4_i8*int(this%nSwap,i8)
#else
    nBytes = 8_i8*int(this%nSwap,i8)
#endif
    nBytes = nBytes + int(11*nbi_p,i8) + int(nbs_p,i8)*int(this%nElSwap,i8)
    if (storeAsDouble .and. this%nSwap > 1) then
       nBytes = nBytes + int(nbd_p,i8)*int(min(maxDoubleBuf_p,nRows),i8)
    end if
    errMsg = StrBytes(nBytes)
    write(lpu,602) trim(errMsg)

    !! Report allocated storage on file
    if (storeAsDouble) then
       nBytes = 46_i8 + int(nbd_p,i8)*dmSize(this)
    else
       nBytes = 46_i8 + int(nbs_p,i8)*dmSize(this)
    end if
    errMsg = StrBytes(nBytes)
    write(lpu,603) trim(errMsg)
    if (nRows < 1) return

    !! Allocate the in-core swap array
    call reAllocate ('dmOpen',this%fileKey,this%nSwap,err)
    call reAllocate ('dmOpen',this%vala,this%nRows,this%nColSwap,err)
    if (err < 0) goto 110

    if (storeAsDouble) then
       !! Allocate single precision buffer (max size: one matrix column)
       call reAllocate ('dmOpen',this%vald,min(maxDoubleBuf_p,this%nRows),err)
       if (err < 0) goto 100
    end if

    !! Read or write through the file to initialize the file keys
    this%curSwap = 1
    nElSwap = this%nElSwap
    do i = 1, this%nSwap
       if (i == this%nSwap) nElSwap = this%nElastSwap

       call getPositionDB (this%swapFile,this%fileKey(i))
       if (this%fileKey(i) == 0_i8) then
          err = -i
          exit
       end if

       if (storeAsDouble) then
          nDouble = min(size(this%vald),nElSwap)
          do j = 1, nElSwap, max(1,nDouble)
             nRead = min(nDouble,nElSwap+1-j)
             call readDoubleDB (this%swapFile,this%vald(1),nRead,err)
             if (err < 0) exit
             if (i == this%nSwap) then ! Last swap section, cast to single
                call dmSingleCast (this,j,nRead)
             end if
          end do
       else
          call readFloatDB (this%swapFile,this%vala(1,1),nElSwap,err)
       end if
       this%curSwap = i
       if (err < 0) exit

    end do
    if (storeAsDouble .and. this%nSwap < 2) then
       !! No need for double precision buffer when entire matrix is in core
       call reAllocate ('dmOpen',this%vald)
    end if
    if (err >= 0) return

    !! Failure to read/write through the disk matrix file
    if (storeAsDouble) then
       nBytes = 46_i8 + dmSize(this)*int(nbd_p,i8)
    else
       nBytes = 46_i8 + dmSize(this)*int(nbs_p,i8)
    end if
    errMsg = 'The file is smaller than the anticipated '//StrBytes(nBytes)
    call reportError (error_p,'Failed to initialize disk matrix file.', &
         &            errMsg,'File: '//fileName,addString='dmOpen')

100 call closeBinaryDB (this%swapFile,lpu)
110 call reAllocate ('dmOpen',this%fileKey)
    call reAllocate ('dmOpen',this%vala)
    call reAllocate ('dmOpen',this%vald)
    call dmNullify (this)

600 format(/4X,'Opening disk matrix file: ',A, &
         & /4X,'   Number of matrix rows     =',I10, &
         & /4X,'   Number of matrix columns  =',I10, &
         & /4X,'   Number of swap sections   =',I10, &
         & /4X,'   Size of swap sections     =',I10,' words')
601 format( 4X,'   Size of last swap section =',I10,' words')
602 format( 4X,'   Allocated storage         =      ',A)
603 format( 4X,'   Total file size           =      ',A)

  contains

    !!==========================================================================
    !> @brief Adjusts the number of matrix columns in a swap section.
    !> @param[in] nTot Total number of matrix columns
    !> @param[in] nSwap Number of swap sections
    !> @return Number of columns in all swap sections, except for the last one
    !> @details The number of columns is set such that the total number of
    !> columns in the matrix is nearly divisable by the number of swap columns,
    !> and such that the last swap section is as close in size as possisble
    !> to the other sections.
    function adjustSwapSize (nTot,nSwap) result(n)
      integer, intent(in) :: nTot,nSwap
      integer :: n
      if (nSwap >= nTot) then
         n = nTot
      else if (nSwap <= 1) then
         n = 1
      else
         n = nSwap
         do while (n > 1 .and. mod(nTot,n) > 0)
            if (nTot/(n-1) > nTot/n .and. mod(nTot,n-1) > 0) exit
            n = n - 1
         end do
      end if
    end function adjustSwapSize

  end subroutine dmOpen


  !!============================================================================
  !> @brief Returns the size of the given disk matrix.
  !>
  !> @param[in] this The sdiskmatrixmodule::diskmatrixtype object to address
  !> @param[in] iDim Which size parameter to return:
  !> - if not present, return the total size (@a nrows&times;ncols)
  !> - 1: number of matrix rows (@a nrows)
  !> - 2: number of matrix columns (@a ncols)
  !> - for any other value: return the number of columns in swap array
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 29 Jan 1999

  function dmSize (this,iDim)

    type(DiskMatrixType), intent(in) :: this
    integer, optional   , intent(in) :: iDim
    integer(i8)                      :: dmSize

    !! --- Logic section ---

    if (.not. present(idim)) then
       dmSize = int(this%nRows,i8)*int(this%nCols,i8)
    else if (iDim == 1) then
       dmSize = int(this%nRows,i8)
    else if (iDim == 2) then
       dmSize = int(this%nCols,i8)
    else
       dmSize = int(this%nColSwap,i8) ! Number of columns in swap array
    end if

  end function dmSize


  !!============================================================================
  !> @brief Closes the disk matrix file and releases the associated swap array.
  !>
  !> @param this The sdiskmatrixmodule::diskmatrixtype object to close
  !> @param[out] err Error indicator
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 1 Dec 2000

  subroutine dmClose (this,err)

    use allocationModule , only : reAllocate
    use reportErrorModule, only : getErrorFile, reportError, debugFileOnly_p
    use binaryDBInterface, only : closeBinaryDB

    type(DiskMatrixType), intent(inout) :: this
    integer             , intent(out)   :: err

    !! Local variables
    integer :: lpu

    !! --- Logic section ---

    err = 0

    if (this%swapCnt(1) > 0 .or. this%swapCnt(2) > 0) then
       lpu = getErrorFile()
       write(lpu,600) this%swapCnt
    end if

    if (this%nElSwap > 0) then
       !! Close the swap file
       call closeBinaryDB (this%swapFile,err)
       if (err < 0) goto 110
    end if

100 call reAllocate ('dmClose',this%fileKey)
    call reAllocate ('dmClose',this%vala)
    call reAllocate ('dmClose',this%vald)
    call dmNullify (this)
    return

110 call reportError (debugFileOnly_p,'dmClose')
    goto 100

600 format(/4X,'Closing disk matrix file:', &
         & /4X,'   Number write operations performed =',I6, &
         & /4X,'   Number read operations performed  =',I6)

  end subroutine dmClose


  !!============================================================================
  !> @brief Gets swap section number and local array index for matrix element.
  !>
  !> @param[in] this The sdiskmatrixmodule::diskmatrixtype object to address
  !> @param[in] rInd Row index of matrix element to get address of
  !> @param[in] cInd Column index of matrix element to get address of
  !> @param[out] section Index of swap section containing the matrix element
  !> @param[out] index Array index within that swap section for matrix element
  !> @param[out] err Error indicator
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 10 Feb 1999

  subroutine dmGetAddress (this,rInd,cInd,section,index,err)

    use reportErrorModule, only : internalError

    type(DiskMatrixType), intent(inout) :: this
    integer             , intent(in)    :: rInd, cInd
    integer             , intent(out)   :: section, index, err

    !! --- Logic section ---

    if ( rInd <= 0 .or. rInd > this%nRows .or. &
         cInd <= 0 .or. cInd > this%nCols ) then
       section = 0
       index = 0
       err = internalError('dmGetAddress: Invalid matrix indices')
    else
       section = (cInd-1)/max(1,this%nColSwap) + 1
       index = cInd - this%nColSwap*(section-1)
       err = 0
    end if

  end subroutine dmGetAddress


  !!============================================================================
  !> @brief Writes current section to disk and reads a new one into core.
  !>
  !> @param[in] this The sdiskmatrixmodule::diskmatrixtype object to swap for
  !> @param[in] newSection Index of the swap section to read into core
  !> @param[out] err Error indicator
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 1 Dec 2000

  subroutine dmSwap (this,newSection,err)

    use reportErrorModule   , only : internalError, reportError, error_p
    use binaryDBInterface   , only : setPositionDB, readFloatDB, readDoubleDB
    use FFaProfilerInterface, only : ffa_starttimer, ffa_stoptimer

    type(DiskMatrixType), intent(inout) :: this
    integer             , intent(in)    :: newSection
    integer             , intent(out)   :: err

    !! Local variables
    integer :: j, nElSwap, nDouble, nRead

    !! --- Logic section ---

    if (newSection == this%curSwap) then
       err = 0
    else if (newSection > this%nSwap) then
       err = internalError('dmSwap: Invalid swap section number')
    else if (newSection > 0) then
       if (newSection == this%nSwap) then
          nElSwap = this%nElastSwap
       else
          nElSwap = this%nElSwap
       end if

       !! Read specified section into memory
       call ffa_starttimer ('dmSwap:read')
       call setPositionDB (this%swapFile,this%fileKey(newSection),err)
       if (err < 0) goto 100
       if (associated(this%vald)) then
          nDouble = min(size(this%vald),nElSwap)
          do j = 1, nElSwap, max(1,nDouble)
             nRead = min(nDouble,nElSwap+1-j)
             call readDoubleDB (this%swapFile,this%vald(1),nRead,err)
             if (err < 0) goto 120
             call dmSingleCast (this,j,nRead)
          end do
       else
          call readFloatDB (this%swapFile,this%vala(1,1),nElSwap,err)
          if (err < 0) goto 120
       end if
       call ffa_stoptimer ('dmSwap:read')
       this%curSwap = newSection
       this%swapCnt(2) = this%swapCnt(2) + 1
       err = 0
    else
       err = 0
    end if

    return

100 call ffa_stoptimer ('dmSwap:read')
    call reportError (error_p,'Could not set swap file pointer', &
         &            addString='dmSwap')
    return

120 call ffa_stoptimer ('dmSwap:read')
    call reportError (error_p,'Could not read from swap file', &
         &            addString='dmSwap')
    return

  end subroutine dmSwap


  !!============================================================================
  !> @brief Casts the in-core matrix elements from double to single precision.
  !>
  !> @param this The sdiskmatrixmodule::diskmatrixtype object to cast for
  !> @param[in] ielSwap Index to first element in swap section
  !> @param[in] nelSwap Number of elements in swap section to cast for
  !>
  !> @callergraph
  !>
  !> @details It is assumed that the double-precision buffer array
  !> sdiskmatrixmodule::diskmatrixtype::vald is not larger than one column.
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 19 May 2005

  subroutine dmSingleCast (this,ielSwap,nelSwap)

    type(DiskMatrixType), intent(inout) :: this
    integer             , intent(in)    :: ielSwap, nelSwap

    !! Local variables
    integer :: i, nr, ic = 1, lastRow = 0

    !! --- Logic section ---

    if (ielSwap == 1) then
       lastRow = 0
       ic = 1
    end if

    nr = min(nelSwap,this%nRows-lastRow)
    do i = 1, nr
       this%vala(lastRow+i,ic) = real(this%vald(i),sp)
    end do

    lastRow = mod(lastRow+nr,this%nRows)
    if (lastRow == 0) ic = ic + 1
    if (nr == nelSwap) return

    nr = nelSwap - nr
    do i = 1, nr
       this%vala(i,ic) = real(this%vald(nelSwap-nr+i),sp)
    end do

    lastRow = nr

  end subroutine dmSingleCast


  !!============================================================================
  !> @brief Returns the value of a matrix element.
  !>
  !> @param this The sdiskmatrixmodule::diskmatrixtype object to return for
  !> @param[in] rInd Row index of matrix element to return
  !> @param[in] cInd Column index of matrix element to return
  !> @param[out] err Error indicator
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 10 Feb 1999

  function dmGetValue (this,rInd,cInd,err)

    use reportErrorModule, only : reportError, debugFileOnly_p

    type(DiskMatrixType), intent(inout) :: this
    integer             , intent(in)    :: rInd, cInd
    integer             , intent(out)   :: err
    real(sp)                            :: dmGetValue

    !! Local variables
    integer :: index, section

    !! --- Logic section ---

    call dmGetAddress (this,rInd,cInd,section,index,err)
    if (err < 0) goto 100

    call dmSwap (this,section,err)
    if (err < 0) goto 100

    dmGetValue = this%vala(rInd,index)
    return

100 call reportError (debugFileOnly_p,'dmGetValue')
    dmGetValue = 0.0_sp

  end function dmGetValue


  !!============================================================================
  !> @brief Returns a pointer to a matrix containing a swap section.
  !>
  !> @param this The sdiskmatrixmodule::diskmatrixtype object to return for
  !> @param[in] cInd Column index for start of swap section to return for
  !> @param[out] err Error indicator
  !>
  !> @details The pointer to the entire swap section starting at column @a cInd
  !> of the given disk matrix. The function can be used for speedy random access
  !> of the matrix elements.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 13 Sep 2004

  function dmGetSwapSecPtr (this,cInd,err)

    use reportErrorModule, only : reportError, debugFileOnly_p

    type(DiskMatrixType), intent(inout) :: this
    integer             , intent(in)    :: cInd
    integer             , intent(out)   :: err
    real(sp)            , pointer       :: dmGetSwapSecPtr(:,:)

    !! Local variables
    integer :: index, section

    !! --- Logic section ---

    call dmGetAddress (this,1,cInd,section,index,err)
    if (err < 0) goto 100

    call dmSwap (this,section,err)
    if (err < 0) goto 100

    dmGetSwapSecPtr => this%vala
    return

100 call reportError (debugFileOnly_p,'dmGetSwapSecPtr')
    nullify(dmGetSwapSecPtr)

  end function dmGetSwapSecPtr


  !!============================================================================
  !> @brief Returns a pointer to column @a cInd of the given disk matrix.
  !>
  !> @param this The sdiskmatrixmodule::diskmatrixtype object to return for
  !> @param[in] cInd Column index to return pointer for
  !> @param[out] err Error indicator
  !>
  !> @note This function circumvents all protection and should be used with
  !> extreme caution and only for ReadOnly purposes when speed is essential.
  !> Calls to other dm... routines can also change the contents of the memory
  !> being pointed to.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen and Karl Erik Thoresen
  !>
  !> @date 5 May 1999

  function dmGetColPtr (this,cInd,err)

    use reportErrorModule, only : reportError, debugFileOnly_p

    type(DiskMatrixType), intent(inout) :: this
    integer             , intent(in)    :: cInd
    integer             , intent(out)   :: err
    real(sp)            , pointer       :: dmGetColPtr(:)

    !! Local variables
    integer :: index, section

    !! --- Logic section ---

    call dmGetAddress (this,1,cInd,section,index,err)
    if (err < 0) goto 100

    call dmSwap (this,section,err)
    if (err < 0) goto 100

    dmGetColPtr => this%vala(:,index)
    return

100 call reportError (debugFileOnly_p,'dmGetColPtr')
    nullify(dmGetColPtr)

  end function dmGetColPtr


  !!============================================================================
  !> @brief Gets a column of the given disk matrix.
  !>
  !> @param this The sdiskmatrixmodule::diskmatrixtype object to return for
  !> @param[in] cInd Column index to return for
  !> @param[out] col Content of the matrix column
  !> @param[out] err Error indicator
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 10 Feb 1999

  subroutine dmGetCol (this,cInd,col,err)

    use reportErrorModule, only : internalError, reportError, debugFileOnly_p

    type(DiskMatrixType), intent(inout) :: this
    integer             , intent(in)    :: cInd
    real(sp)            , intent(out)   :: col(:)
    integer             , intent(out)   :: err

    !! Local variables
    integer :: index, section

    !! --- Logic section ---

    if (size(col) /= this%nRows) then
       err = internalError('dmGetCol: Matrix dimension mismatch')
       return
    end if

    call dmGetAddress (this,1,cInd,section,index,err)
    if (err < 0) goto 100

    call dmSwap (this,section,err)
    if (err < 0) goto 100

    call SCOPY (this%nRows,this%vala(1,index),1,col(1),1)
    return

100 call reportError (debugFileOnly_p,'dmGetCol')

  end subroutine dmGetCol


  !!============================================================================
  !> @brief Gets a row of the given disk matrix.
  !>
  !> @param this The sdiskmatrixmodule::diskmatrixtype object to return for
  !> @param[in] rInd Row index to return for
  !> @param[out] row Content of the matrix row
  !> @param[out] err Error indicator
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 10 Feb 1999

  subroutine dmGetRow (this,rInd,row,err)

    use reportErrorModule, only : internalError, reportError, debugFileOnly_p

    type(DiskMatrixType), intent(inout) :: this
    integer             , intent(in)    :: rInd
    real(sp)            , intent(out)   :: row(:)
    integer             , intent(out)   :: err

    !! Local variables
    integer :: i, iSect, cInd, index, section

    !! --- Logic section ---

    if (size(row) < this%nCols) then
       err = internalError('dmGetRow: Matrix dimension mismatch')
       return
    end if

    cInd = 1
    do iSect = 1, this%nSwap
       call dmGetAddress (this,rInd,cInd,section,index,err)
       if (err < 0) exit

       call dmSwap (this,section,err)
       if (err < 0) exit

       do i = 1, this%nColSwap
          if (cInd > this%nCols) exit ! Last section may be less than nColSwap
          row(cInd) = this%vala(rInd,index)
          cInd = cInd + 1
          index = index + 1
       end do

    end do

    if (err < 0) call reportError (debugFileOnly_p,'dmGetRow')

  end subroutine dmGetRow


  !!============================================================================
  !> @brief Multiplies a disk matrix with a given single precision vector.
  !>
  !> @param this The sdiskmatrixmodule::diskmatrixtype object (@b A) to multiply
  !> @param[in] xVec The vector (@b x) to multiply matrix @b A with
  !> @param yVec The resulting vector (@b y)
  !> @param[out] err Error indicator
  !> @param[in] doInitialize If not present or .true., the existing content
  !>            of @a yVec is discarded
  !>
  !> @details One of of the following operations are performed
  !>     y = A*x
  !>     y = y + A*x
  !> depending on the value of @a doInitialize.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 10 Feb 1999

  subroutine dmMatTimesVec_SP (this,xVec,yVec,err,doInitialize)

    use reportErrorModule, only : internalError, reportError, debugFileOnly_p

    type(DiskMatrixType), intent(inout) :: this
    real(sp)            , intent(in)    :: xVec(:)
    real(sp)            , intent(inout) :: yVec(:)
    integer             , intent(out)   :: err
    logical, optional   , intent(in)    :: doInitialize

    !! Local variables
    integer :: i, nRows, section, index

    !! --- Logic section ---

    err = 0
    if (this%nRows < 1) return
    if (this%nCols < 1) return

    nRows = this%nRows
    if (size(xVec) < this%nCols .or. size(yVec) /= nRows) then
       err = internalError('dmMatTimesVec: Matrix dimension mismatch')
       return
    end if

    if (.not. present(doInitialize)) then
       call SCOPY (nRows,0.0_sp,0,yVec(1),1)
    else if (doInitialize) then
       call SCOPY (nRows,0.0_sp,0,yVec(1),1)
    end if

    !! Perform the multiplication using SAXPY in order to finish the columns

    do i = 1, this%nCols
       call dmGetAddress (this,1,i,section,index,err)
       if (err < 0) exit

       call dmSwap (this,section,err)
       if (err < 0) exit

       call SAXPY (nRows,xVec(i),this%vala(1,index),1,yVec(1),1)
    end do

    if (err < 0) call reportError (debugFileOnly_p,'dmMatTimesVec')

  end subroutine dmMatTimesVec_SP


  !!============================================================================
  !> @brief Multiplies a disk matrix with a given double precision vector.
  !>
  !> @param this The sdiskmatrixmodule::diskmatrixtype object (@b A) to multiply
  !> @param[in] xVec The vector (@b x) to multiply matrix @b A with
  !> @param yVec The resulting vector (@b y)
  !> @param[out] err Error indicator
  !> @param[in] doInitialize If not present or .true., the existing content
  !>            of @a yVec is discarded
  !>
  !> @details One of of the following operations are performed
  !>     y = A*x
  !>     y = y + A*x
  !> depending on the value of @a doInitialize.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 10 Feb 1999

  subroutine dmMatTimesVec_DP (this,xVec,yVec,err,doInitialize)

    use reportErrorModule, only : internalError, reportError, debugFileOnly_p

    type(DiskMatrixType), intent(inout) :: this
    real(dp)            , intent(in)    :: xVec(:)
    real(sp)            , intent(inout) :: yVec(:)
    integer             , intent(out)   :: err
    logical, optional   , intent(in)    :: doInitialize

    !! Local variables
    integer :: i, nRows, section, index

    !! --- Logic section ---

    err = 0
    if (this%nRows < 1) return
    if (this%nCols < 1) return

    nRows = this%nRows
    if (size(xVec) < this%nCols .or. size(yVec) /= nRows) then
       err = internalError('dmMatTimesVec: Matrix dimension mismatch')
       return
    end if

    if (.not. present(doInitialize)) then
       call SCOPY (nRows,0.0_sp,0,yVec(1),1)
    else if (doInitialize) then
       call SCOPY (nRows,0.0_sp,0,yVec(1),1)
    end if

    !! Perform the multiplication using SAXPY in order to finish the columns

    do i = 1, this%nCols
       call dmGetAddress (this,1,i,section,index,err)
       if (err < 0) exit

       call dmSwap (this,section,err)
       if (err < 0) exit

       call SAXPY (nRows,real(xVec(i),sp),this%vala(1,index),1,yVec(1),1)
    end do

    if (err < 0) call reportError (debugFileOnly_p,'dmMatTimesVec')

  end subroutine dmMatTimesVec_DP


  !!============================================================================
  !> @brief Multiplies the transpose of a disk matrix with a given vector.
  !>
  !> @param this The sdiskmatrixmodule::diskmatrixtype object (@b A) to multiply
  !> @param[in] xVec The vector (@b x) to multiply matrix @b A with
  !> @param yVec The resulting vector (@b y)
  !> @param[out] err Error indicator
  !>
  !> @details The following operation is performed
  !>     y = A^t*x
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 10 Feb 1999

  subroutine dmMatTransTimesVec (this,xVec,yVec,err)

    use reportErrorModule, only : internalError, reportError, debugFileOnly_p

    type(DiskMatrixType), intent(inout) :: this
    real(sp)            , intent(in)    :: xVec(:)
    real(sp)            , intent(out)   :: yVec(:)
    integer             , intent(out)   :: err

    !! Local variables
    integer :: i, nRows, section, index

    real(sp), external :: SDOT

    !! --- Logic section ---

    err = 0
    if (this%nRows < 1) return
    if (this%nCols < 1) return

    nRows = this%nRows
    if (size(xVec) /= nRows .or. size(yVec) < this%nCols) then
       err = internalError('dmMatTransTimesVec: Matrix dimension mismatch')
       return
    end if

    !! Perform the multiplication using dot-products
    !! in order to work in a column-like fashion

    do i = 1, this%nCols
       call dmGetAddress (this,1,i,section,index,err)
       if (err < 0) exit

       call dmSwap (this,section,err)
       if (err < 0) exit

       yVec(i) = SDOT(nRows,this%vala(1,index),1,xVec(1),1)
    end do

    if (err < 0) call reportError (debugFileOnly_p,'dmMatTransTimesVec')

  end subroutine dmMatTransTimesVec


  !!============================================================================
  !> @brief Finds the absolute maximum value for each row in a disk matrix.
  !>
  !> @param this The sdiskmatrixmodule::diskmatrixtype object to address
  !> @param[out] maxRowVal Absolute maximum row values
  !> @param[out] err Error indicator
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 25 Mar 2003

  subroutine dmFindAbsMaxRowValues (this,maxRowVal,err)

    use reportErrorModule, only : internalError, reportError, debugFileOnly_p

    type(DiskMatrixType), intent(inout) :: this
    real(sp)            , intent(out)   :: maxRowVal(:)
    integer             , intent(out)   :: err

    !! Local variables
    integer  :: i, j, nRows, section, index
    real(sp) :: colVal

    !! --- Logic section ---

    err = 0
    if (this%nRows < 1) return
    if (this%nCols < 1) return

    nRows = this%nRows
    if (size(maxRowVal) < nRows) then
       err = internalError('dmFindAbsMaxRowValues: Matrix dimension mismatch')
       return
    end if

    !! Search for largest row value column by column

    do i = 1, this%nCols
       call dmGetAddress (this,1,i,section,index,err)
       if (err < 0) exit

       call dmSwap (this,section,err)
       if (err < 0) exit

       if (i == 1) then
          maxRowVal(1:nRows) = this%vala(:,index)
       else
          do j = 1, nRows
             colVal = this%vala(j,index)
             if (abs(colVal) > abs(maxRowVal(j))) maxRowVal(j) = colVal
          end do
       end if
    end do

    if (err < 0) call reportError (debugFileOnly_p,'dmFindAbsMaxRowValues')

  end subroutine dmFindAbsMaxRowValues

end module sDiskMatrixModule
