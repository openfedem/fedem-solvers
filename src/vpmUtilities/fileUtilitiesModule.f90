!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file fileUtilitiesModule.f90
!> @brief Subroutines for file handling.

!!==============================================================================
!> @brief Module with subroutines for file handling.
!> @details This module contains various utility subroutines and functions
!> for management of file unit numbers that mostly are used for debug output.
!> You have to configure doxygen with the option ENABLED_SECTIONS = FULL_DOC
!> to extract the detailed documentation of the symbols in this module.
!>
!> @author Bjorn Haugen
!> @date Sep 2002

module FileUtilitiesModule
  !> @cond FULL_DOC

  implicit none

  !> Array of unit numbers for the currently open debug files.
  integer, pointer, save, private :: dbg_files(:) => null()


contains

  !!============================================================================
  !> @brief Returns a file unit number to a debug file.
  !> @details If @a fileName is not provided, the file is instead named
  !> "debug_<idx>.dbg". The unit number is also stored in @a dbg_files(idx).
  !> The array #dbg_files is automatically reallocated if the @a idx value
  !> provided is larger than the current size. If the named file already exists,
  !> it is overwritten unless @a newfile is specified and equals .true.
  !> In the latter case a new file name is constructed.
  function GetDBGfile (idx,fileName,newfile) result(file)

    use progressModule, only : lterm

    integer         , intent(in)           :: idx
    character(len=*), intent(in), optional :: fileName
    logical         , intent(in), optional :: newfile

    !! Local variables
    character(len=64) :: dbgName
    character(len=8)  :: ext
    integer, pointer  :: old_dbg_files(:)
    integer           :: file, nFiles, i, j

    !! --- Logic section ---

    if (associated(dbg_files)) then
       nFiles = size(dbg_files)
    else
       nFiles = 0
    end if

    !! Reallocate if dbg_files is not large enough
    if (idx > nFiles) then
       if (nFiles > 0) write(lterm,600) nFiles,idx
       old_dbg_files => dbg_files
       allocate(dbg_files(max(10,idx)),stat=file)
       if (file == 0) then
          if (nFiles > 0) dbg_files(1:nFiles) = old_dbg_files(1:nFiles)
          dbg_files(nFiles+1:) = 0
       else
          nullify(dbg_files)
       end if
       if (associated(old_dbg_files)) deallocate(old_dbg_files)
    else
       file = 0
    end if

    !! Open a new file, if necessary
    if (idx < 1 .or. file /= 0) then
       file = 0 ! Use stderr if idx is invalid
       return
    else if (dbg_files(idx) == 0) then
       dbg_files(idx) = findUnitNumber(80)
       if (present(fileName)) then
          dbgName = fileName
       else
          write(dbgName,"(i4)") idx
          dbgName = 'debug_'//trim(adjustl(dbgName))//'.dbg'
       end if
       if (present(newfile)) then
          i = index(dbgName,'.')
          if (newfile .and. i > 0) then
             j = 1 ! Create new filename, to avoid overwriting existing one
             ext = dbgName(i:)
             dbgName(i:) = '001'//ext
             do while (isFile(dbgName))
                j = j + 1
                write(dbgName(i:i+2),"(i3.3)") j
             end do
          end if
       end if
       write(lterm,610) trim(dbgName),idx,dbg_files(idx)
       open(dbg_files(idx),FILE=trim(dbgName),STATUS='UNKNOWN')
    end if

    !! Return the file pointer
    file = dbg_files(idx)

600 format('DBG file: Reallocating dgb_files array:',i4,'  -->',i4)
610 format('DBG file: Opening ',a,', index =',i4,', Fortran unit =',i4)

  end function GetDBGfile


  !!============================================================================
  !> @brief Closes and deallocates all debug files.
  !> @author Knut Morten Okstad
  !> @date Jan 2017
  subroutine CloseDBGfiles ()

    integer :: i

    if (associated(dbg_files)) then
       do i = 1, size(dbg_files)
          if (dbg_files(i) > 6) close(dbg_files(i))
       end do
       deallocate(dbg_files)
       nullify(dbg_files)
    end if

  end subroutine CloseDBGfiles


  !!============================================================================
  !> @brief Checks if a named file exists.
  function isFile (fileName) result(exists)

    character(len=*), intent(in) :: fileName
    logical                      :: exists

    if (fileName /= '') then
       inquire(FILE=filename,EXIST=exists)
    else
       exists = .false.
    end if

  end function isFile


  !!============================================================================
  !> @brief Returns a FORTRAN file unit number currently not in use.
  function findUnitNumber (iStart) result(iUnit)

    integer, intent(in) :: istart
    integer             :: iUnit
    logical             :: exists, named

    iUnit = iStart-1
    named = .true.
    do while (named)
       iUnit = iUnit + 1
       inquire(iunit,EXIST=exists,NAMED=named)
    end do

  end function findUnitNumber


  !!============================================================================
  !> @brief Prints all open files with unit number in the range [iStart,iEnd].
  subroutine printOpenFiles (iStart,iEnd,iOut)

    integer, intent(in) :: istart, iEnd, iOut
    integer             :: iUnit
    logical             :: exists, named
    character(len=128)  :: fileName

    iUnit = iStart-1
    named = .true.
    do iUnit = iStart, iEnd
       inquire(iunit,EXIST=exists,NAMED=named,NAME=fileName)
       if (named) then
          write(iOut,"(' File unit',I6,' exists, name: ',A)") iUnit, fileName
       end if
    end do

  end subroutine printOpenFiles

  !> @endcond
  !!============================================================================
  !> @brief Extracts a file name from a command-line option.
  !>
  !> @param[in]  file The command-line option to extract from
  !> @param[out] name The value of the command-line option
  !> @param[in]  suffix File extension to use when using default file
  !>
  !> @details If needed the string value is converted to UNIX or Windows format.
  !> If the specified command-line option is not used, the default file name is
  !> `<partname><suffix>` where `<partname>` is the name of the FE data file
  !> without the file extension, as provided with the `-linkfile` option.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 20 Sep 2000

  subroutine getFileName (file,name,suffix)

    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getstring
    use FFaFilePathInterface  , only : ffa_checkpath

    character(len=*),           intent(in)  :: file
    character(len=*),           intent(out) :: name
    character(len=*), optional, intent(in)  :: suffix

    !! Local variables
    integer :: idot

    !! --- Logic section ---

    call ffa_cmdlinearg_getstring (file,name)
    if (name == '' .and. present(suffix)) then
       !! File name not given, use default
       call ffa_cmdlinearg_getstring ('linkfile',name)
       if (name == '') then
          name = 'fedem'//suffix
          return
       end if
       idot = index(name,'.',.TRUE.)
       if (idot < 1) idot = len_trim(name)+1
       name(idot:) = suffix
    end if

    !! Unix/Windows pathname conversion
    call ffa_checkpath (name)

  end subroutine getFileName

end module FileUtilitiesModule
