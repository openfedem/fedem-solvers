!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file binaryDBInterface.f90
!> @brief Fortran interface for the C functions for doing binary IO.

!!==============================================================================
!> @brief Fortran interface for the C functions for doing binary IO.

module binaryDBInterface

  implicit none

  integer, parameter :: tmp_p = 0    !< Opening mode for temporary files
  integer, parameter :: read_p = 1   !< Read file opening mode
  integer, parameter :: write_p = 2  !< Write file opening mode
  integer, parameter :: append_p = 3 !< Append file opening mode
  integer, parameter :: update_p = 4 !< Update file opening mode
  integer, parameter :: reawri_p = 5 !< Read/write file opening mode

  private :: closeBinary_DB, closeBinaryDB1, closeBinaryDB2
  private :: readDoubleD4, readDoubleD8, writeDoubleD4, writeDoubleD8
  private :: readDoubleDB4, readDoubleDB8, readDoubleDBfile
  private :: writeDoubleDB4, writeDoubleDB8, writeDoubleDBfile

  !> @brief Reads a double precision array from the specified binary file.
  interface readDoubleDB
     module procedure readDoubleDB4, readDoubleDB8
     module procedure readDoubleDBfile
  end interface

  !> @brief Writes a double precision array to the specified binary file.
  interface writeDoubleDB
     module procedure writeDoubleDB4, writeDoubleDB8
     module procedure writeDoubleDBfile
  end interface

  !> @brief Closes the specified binary file.
  interface closeBinaryDB
     module procedure closeBinaryDB1, closeBinaryDB2
  end interface


  interface

     !> @brief Opens a binary direct access file for read or write.
     !> @param[in] fileName Name of file to open
     !> @param[in] fileType Type of file to open
     !> -  0 : Temporary file with no name, automatically deleted when closed
     !> -  1 : Existing file, opened for read-only
     !> -  2 : New file, opened for write
     !> -  3 : Existing file, opened for append
     !> -  4 : New file, opened for read and write
     !> -  5 : Existing file, opened for read and write
     !> @param[out] fileNum Assigned file handle (non-negative number)
     !> @param[out] status Exit status (negative on error)
     !> -  0 : Everything is OK
     !> - -1 : Could not open the specified file
     !> - -2 : Too many simultaneously opened files
     subroutine openBinaryDB (fileName, fileType, fileNum, status)
       character, intent(in)  :: filename*(*)
       integer  , intent(in)  :: fileType
       integer  , intent(out) :: fileNum, status
     end subroutine openBinaryDB

     !> @brief Closes the specified binary file.
     !> @param[in] fileNum File handle for the file to close
     !> @param[in] forceDelete If @e .true., force deleting the file after close
     !> @param[out] status Exit status (negative on error)
     subroutine closeBinary_DB (fileNum, forceDelete, status)
       integer, intent(in)  :: fileNum
       logical, intent(in)  :: forceDelete
       integer, intent(out) :: status
     end subroutine closeBinary_DB

     !> @brief Deletes the named file.
     !> @param[in] fileName Name of file to delete
     subroutine deleteDB (fileName)
       character, intent(in) :: filename*(*)
     end subroutine deleteDB

     !> @brief Allocates an in-core buffer for the specified binary file.
     !> @param[in] fileNum File handle for the file to allocate for
     !> @param[in] nBytes Number of bytes in the buffer
     !> @param[out] status Exit status (negative on error)
     subroutine setBufDB (fileNum, nBytes, status)
       integer, intent(in)  :: fileNum, nBytes
       integer, intent(out) :: status
     end subroutine setBufDB

     !> @brief Flushes the in-core buffer of the specified binary file to disk.
     !> @param[in] fileNum File handle for the file to flush
     !> @param[out] status Exit status (negative on error)
     subroutine flushBinaryDB (fileNum, status)
       integer, intent(in)  :: fileNum
       integer, intent(out) :: status
     end subroutine flushBinaryDB

     !> @brief Writes an integer array to the specified binary file.
     !> @param[in] fileNum File handle for the file to write to
     !> @param[in] data The integer array to write
     !> @param[in] nData Number of words to write
     !> @param[out] status Exit status
     !> (negative on error, otherwise number of bytes written)
     subroutine writeIntDB (fileNum, data, nData, status)
       integer, intent(in)  :: fileNum, nData
       integer, intent(in)  :: data
       integer, intent(out) :: status
     end subroutine writeIntDB

     !> @brief Writes a single precision array to the specified binary file.
     !> @param[in] fileNum File handle for the file to write to
     !> @param[in] data The real array to write
     !> @param[in] nData Number of words to write
     !> @param[out] status Exit status
     !> (negative on error, otherwise number of bytes written)
     subroutine writeFloatDB (fileNum, data, nData, status)
       use KindModule, only  :  sp
       integer , intent(in)  :: fileNum, nData
       real(sp), intent(in)  :: data
       integer , intent(out) :: status
     end subroutine writeFloatDB

     !> @brief Writes a double precision array to the specified binary file.
     !> @param[in] fileNum File handle for the file to write to
     !> @param[in] data The real array to write
     !> @param[in] nData Number of words to write
     !> @param[out] status Exit status
     !> (negative on error, otherwise number of bytes written)
     subroutine writeDoubleD4 (fileNum, data, nData, status)
       use KindModule, only  :  dp
       integer , intent(in)  :: fileNum, nData
       real(dp), intent(in)  :: data
       integer , intent(out) :: status
     end subroutine writeDoubleD4

     !> @brief Writes a big double precision array to the specified binary file.
     !> @param[in] fileNum File handle for the file to write to
     !> @param[in] data The real array to write
     !> @param[in] nData Number of words to write
     !> @param[out] status Exit status
     !> (negative on error, otherwise number of bytes written)
     subroutine writeDoubleD8 (fileNum, data, nData, status)
       use KindModule, only     :  i8, dp
       integer    , intent(in)  :: fileNum
       integer(i8), intent(in)  :: nData
       real(dp)   , intent(in)  :: data
       integer    , intent(out) :: status
     end subroutine writeDoubleD8

     !> @brief Writes a character string to the specified binary file.
     !> @param[in] fileNum File handle for the file to write to
     !> @param[in] data The character string to write
     !> @param[out] status Exit status
     !> (negative on error, otherwise number of bytes written)
     subroutine writeCharDB (fileNum, data, status)
       integer  , intent(in)  :: fileNum
       character, intent(in)  :: data*(*)
       integer  , intent(out) :: status
     end subroutine writeCharDB

     !> @brief Writes a file tag and checksum to the specified binary file.
     !> @param[in] fileNum File handle for the file to write to
     !> @param[in] tag The file tag to write
     !> @param[in] chksum The file checksum to write
     !> @param[out] status Exit status
     !> (negative on error, otherwise number of bytes written)
     subroutine writeTagDB (fileNum, tag, chksum, status)
       integer  , intent(in)  :: fileNum, chksum
       character, intent(in)  :: tag*(*)
       integer  , intent(out) :: status
     end subroutine writeTagDB

     !> @brief Reads an integer array from a specified binary file.
     !> @param[in] fileNum File handle for the file to read from
     !> @param[in] data The integer array to read into
     !> @param[in] nData Number of words to read
     !> @param[out] status Exit status
     !> (negative on error, otherwise number of bytes read)
     subroutine readIntDB (fileNum, data, nData, status)
       integer, intent(in)  :: fileNum, nData
       integer, intent(out) :: data
       integer, intent(out) :: status
     end subroutine readIntDB

     !> @brief Reads a single precision array from a specified binary file.
     !> @param[in] fileNum File handle for the file to read from
     !> @param[in] data The real array to read into
     !> @param[in] nData Number of words to read
     !> @param[out] status Exit status
     !> (negative on error, otherwise number of bytes read)
     subroutine readFloatDB (fileNum, data, nData, status)
       use KindModule, only  :  sp
       integer , intent(in)  :: fileNum, nData
       real(sp), intent(out) :: data
       integer , intent(out) :: status
     end subroutine readFloatDB

     !> @brief Reads a double precision array from a specified binary file.
     !> @param[in] fileNum File handle for the file to read from
     !> @param[in] data The real array to read into
     !> @param[in] nData Number of words to read
     !> @param[out] status Exit status
     !> (negative on error, otherwise number of bytes read)
     subroutine readDoubleD4 (fileNum, data, nData, status)
       use KindModule, only  :  dp
       integer , intent(in)  :: fileNum, nData
       real(dp), intent(out) :: data
       integer , intent(out) :: status
     end subroutine readDoubleD4

     !> @brief Reads a big double precision array from a specified binary file.
     !> @param[in] fileNum File handle for the file to read from
     !> @param[in] data The real array to read into
     !> @param[in] nData Number of words to read
     !> @param[out] status Exit status
     !> (negative on error, otherwise number of bytes read)
     subroutine readDoubleD8 (fileNum, data, nData, status)
       use KindModule, only     :  i8, dp
       integer    , intent(in)  :: fileNum
       integer(i8), intent(in)  :: nData
       real(dp)   , intent(out) :: data
       integer    , intent(out) :: status
     end subroutine readDoubleD8

     !> @brief Reads a character string from a specified binary file.
     !> @param[in] fileNum File handle for the file to read from
     !> @param[in] data The character string to read into
     !> @param[in] nData Number of characters to read
     !> @param[out] status Exit status
     !> (negative on error, otherwise number of bytes read)
     subroutine readCharDB (fileNum, data, nData, status)
       integer  , intent(in)  :: fileNum, nData
       character, intent(out) :: data*(*)
       integer  , intent(out) :: status
     end subroutine readCharDB

     !> @brief Reads the file tag and checksum from a specified binary file.
     !> @param[in] fileNum File handle for the file to read from
     !> @param[out] tag The character string to read the file tag into
     !> @param[out] chksum The file checksum
     !> @param[out] status Exit status
     !> (negative on error, otherwise number of bytes read)
     subroutine readTagDB (fileNum, tag, chksum, status)
       integer  , intent(in)  :: fileNum
       character, intent(out) :: tag*(*)
       integer  , intent(out) :: chksum, status
     end subroutine readTagDB

     !> @brief Sets the read position for a specified binary file.
     !> @param[in] fileNum File handle for the file to set position for
     !> @param[in] pos The file position to set
     !> @param[out] status Exit status (negative on error, otherwise zero)
     subroutine setPositionDB (fileNum, pos, status)
       use KindModule, only     :  i8
       integer    , intent(in)  :: fileNum
       integer(i8), intent(in)  :: pos
       integer    , intent(out) :: status
     end subroutine setPositionDB

     !> @brief Gets the current position for a specified binary file.
     !> @param[in] fileNum File handle for the file to get position for
     !> @param[out] pos The current file position (zero on error)
     subroutine getPositionDB (fileNum, pos)
       use KindModule, only     :  i8
       integer    , intent(in)  :: fileNum
       integer(i8), intent(out) :: pos
     end subroutine getPositionDB

     !> @brief Writes a character string at a specified file location.
     !> @param[in] fileNum File handle for the file to write to
     !> @param[in] afterTag Tag identifying the file position to write to
     !> @param[in] data The character string to write
     !> @param[out] status Exit status
     !> (negative on error, otherwise number of bytes written)
     !>
     !> @details The actual location is identified by the string @a afterTag.
     !> Any existing data at the specified location will be overwritten such
     !> that the total file size will not change, unless at the end of the file.
     !> If the @a afterTag is not found, the file is not touched.
     subroutine putCharDB (fileNum, afterTag, data, status)
       integer  , intent(in)  :: fileNum
       character, intent(in)  :: afterTag*(*), data*(*)
       integer  , intent(out) :: status
     end subroutine putCharDB

     !> @brief Copies data from one binary file to another.
     !> @param[in] tofile File handle for the file to copy to
     !> @param[in] fromfile File handle for the file to copy from
     !> @param[in] nbytes Number of bytes to copy
     !> @param[out] status Exit status
     !> (negative on error, otherwise number of bytes copied)
     subroutine copyBinaryDB (tofile, fromfile, nbytes, status)
       integer  , intent(in)  :: tofile, fromfile, nbytes
       integer  , intent(out) :: status
     end subroutine copyBinaryDB

  end interface


contains

  !!============================================================================
  !> @brief Reads a double precision array from a specified binary file.
  !> @param[in] fileNum File handle for the file to read from
  !> @param[in] data The real array to read into
  !> @param[in] nData Number of words to read
  !> @param[out] status Exit status
  subroutine readDoubleDB4 (fileNum, data, nData, status)
    use KindModule, only  :  dp
    integer , intent(in)  :: fileNum, nData
    real(dp), intent(out) :: data
    integer , intent(out) :: status
    call readDoubleD4 (fileNum,data,nData,status)
  end subroutine readDoubleDB4

  !> @brief Reads a big double precision array from a specified binary file.
  !> @param[in] fileNum File handle for the file to read from
  !> @param[in] data The real array to read into
  !> @param[in] nData Number of words to read
  !> @param[out] status Exit status
  subroutine readDoubleDB8 (fileNum, data, nData, status)
    use KindModule, only     :  i8, dp
    integer    , intent(in)  :: fileNum
    integer(i8), intent(in)  :: nData
    real(dp)   , intent(out) :: data
    integer    , intent(out) :: status
    call readDoubleD8 (fileNum,data,nData,status)
  end subroutine readDoubleDB8


  !!============================================================================
  !> @brief Writes a double precision array to the specified binary file.
  !> @param[in] fileNum File handle for the file to write to
  !> @param[in] data The real array to write
  !> @param[in] nData Number of words to write
  !> @param[out] status Exit status
  subroutine writeDoubleDB4 (fileNum, data, nData, status)
    use KindModule, only  :  dp
    integer , intent(in)  :: fileNum, nData
    real(dp), intent(in)  :: data
    integer , intent(out) :: status
    call writeDoubleD4 (fileNum,data,nData,status)
  end subroutine writeDoubleDB4

  !> @brief Writes a big double precision array to the specified binary file.
  !> @param[in] fileNum File handle for the file to write to
  !> @param[in] data The real array to write
  !> @param[in] nData Number of words to write
  !> @param[out] status Exit status
  subroutine writeDoubleDB8 (fileNum, data, nData, status)
    use KindModule, only     :  i8, dp
    integer    , intent(in)  :: fileNum
    integer(i8), intent(in)  :: nData
    real(dp)   , intent(in)  :: data
    integer    , intent(out) :: status
    call writeDoubleD8 (fileNum,data,nData,status)
  end subroutine writeDoubleDB8


  !!============================================================================
  !> @brief Closes the specified binary file.
  !> @param[in] fileNum File handle for the file to close
  !> @param[out] status Exit status (negative on error)
  subroutine closeBinaryDB1 (fileNum, status)
    integer, intent(in)  :: fileNum
    integer, intent(out) :: status
    call closeBinary_DB (fileNum,.false.,status)
  end subroutine closeBinaryDB1

  !> @brief Closes the specified binary file.
  !> @param[in] fileNum File handle for the file to close
  !> @param[in] forceDelete If @e .true., force deleting the file after close
  !> @param[out] status Exit status (negative on error)
  subroutine closeBinaryDB2 (fileNum, forceDelete, status)
    integer, intent(in)  :: fileNum
    logical, intent(in)  :: forceDelete
    integer, intent(out) :: status
    call closeBinary_DB (fileNum,forceDelete,status)
  end subroutine closeBinaryDB2


  !!============================================================================
  !> @brief Reads a double precision matrix from a named binary file into core.
  !>
  !> @param[in] fileName Name of the file to read
  !> @param[in] dtype File tag identifying the correct data type
  !> @param[in] ndata Number of words to read
  !> @param[in] data The real array to read into
  !> @param ierr Error flag
  !> @param jerr Optional error flag for the read operation
  !>
  !> @author Knut Morten Okstad
  !>
  !> date 21 Dec 2016

  subroutine readDoubleDBfile (fileName,dtype,ndata,data,ierr,jerr)

    use kindModule       , only : dp
    use reportErrorModule, only : reportError, error_p

    character(len=*), intent(in)    :: fileName, dtype
    integer         , intent(in)    :: ndata
    real(dp)        , intent(out)   :: data
    integer         , intent(inout) :: ierr
    integer,optional, intent(inout) :: jerr

    !! Local variables
    integer           :: ifile, chksum, status
    character(len=32) :: ctag

    !! --- Logic section ---

    call openBinaryDB (fileName,read_p,ifile,status)
    if (status < 0) then
       call reportError (error_p,'Failed to open binary file '//fileName, &
            &            addString='readDoubleDBfile')
       ierr = ierr + status
       return
    end if

    call readTagDB (ifile,ctag,chksum,status)
    if (status < 0) then
       call reportError (error_p,'Failed to extract file tag from '//fileName, &
            &            addString='readDoubleDBfile')
       ierr = ierr + status
    else if (ctag /= '#FEDEM '//dtype) then
       call reportError (error_p,trim(fileName)// &
            &            ' is not a '//dtype //' file','tag='//trim(ctag), &
            &            addString='readDoubleDBfile')
       ierr = ierr - 1
    else
       call readDoubleD4 (ifile,data,ndata,status)
       if (status < 0) then
          if (present(jerr)) then
             jerr = jerr + status
          else
             ierr = ierr + status
          end if
       end if
    end if

    call closeBinaryDB (ifile,status)

  end subroutine readDoubleDBfile


  !!============================================================================
  !> @brief Writes a double precision matrix to a named binary file.
  !>
  !> @param[in] fileName Name of the file to write
  !> @param[in] dtype File tag identifying the data type to write
  !> @param[in] chksum File checksum
  !> @param[in] data The real array to write
  !> @param ierr Error flag
  !> @param sbuf Single precision buffer used for type casting
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 3 Jan 2005

  subroutine writeDoubleDBfile (fileName,dtype,chksum,data,ierr,sbuf)

    use kindModule, only : sp, dp

    character(len=*) , intent(in)  :: fileName, dtype
    integer          , intent(in)  :: chksum
    real(dp)         , intent(in)  :: data(:,:)
    integer          , intent(out) :: ierr
    real(sp),optional, intent(out) :: sbuf(:)

    !! Local variables
    integer :: i, j, ifile, nvec

    !! --- Logic section ---

    call openBinaryDB (fileName,write_p,ifile,ierr)
    if (ierr < 0) return

    if (present(sbuf)) then

       !! Store as single-precision (half file size)
       nvec = size(data,1)
       ierr = size(sbuf) - nvec
       if (ierr < 0) return

       call writeTagDB (ifile,'#FEDEM '//dtype//' SP',chksum,ierr)
       if (ierr < 0) return

       do i = 1, size(data,2)
          do j = 1, nvec
             sbuf(j) = real(data(j,i),sp)
          end do
          call writeFloatDB (ifile,sbuf(1),nvec,ierr)
          if (ierr < 0) return
       end do

    else

       call writeTagDB (ifile,'#FEDEM '//dtype,chksum,ierr)
       if (ierr < 0) return

       call writeDoubleDB (ifile,data(1,1),size(data),ierr)
       if (ierr < 0) return

    end if

    call closeBinaryDB (ifile,ierr)

  end subroutine writeDoubleDBfile

end module binaryDBInterface
