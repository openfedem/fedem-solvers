!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file reportErrorModule.f90
!> @brief Subroutines for error message handling.

!!==============================================================================
!> @brief Module with subroutines for printing error messages, etc.
!> @details This module contains various utility subroutines and functions
!> for printing error messages to the console unit or a designated log-file.
!> You have to configure doxygen with the option ENABLED_SECTIONS = FULL_DOC
!> to extract the detailed documentation of the symbols in this module.

module ReportErrorModule
  !> @cond FULL_DOC

  implicit none

  integer, save, private :: errorFile = 915 !< Default error file unit number

  !! Legal errCode values in ReportError          Formatting of message:
  integer, parameter :: empty_p   = 0          !< "          <message>"
  integer, parameter :: note_p    = 1          !< "Note    : <message>"
  integer, parameter :: warning_p = 2          !< "Warning : <message>"
  integer, parameter :: error_p   = 3          !< "Error   : <message>"
  integer, parameter :: debug_p   = 4          !< "Debug   : <message>"
  integer, parameter :: noteFileOnly_p    = -1 !< "Note    : <message>"
  integer, parameter :: warningFileOnly_p = -2 !< "Warning : <message>"
  integer, parameter :: errorFileOnly_p   = -3 !< "Error   : <message>"
  integer, parameter :: debugFileOnly_p   = -4 !< "Debug   : <message>"

  !> Error message prefices
  character(len=10), parameter :: preString_p(0:4) = (/ '          ', &
       &                                                'Note :    ', &
       &                                                'Warning : ', &
       &                                                'Error :   ', &
       &                                                'Debug :   ' /)


contains

  !!============================================================================
  !> @brief Sets the file unit number to be used for error messages.

  subroutine SetErrorFile (lpuErr)

    integer, intent(in) :: lpuErr !< File unit number to use for error messages

    errorFile = lpuErr

  end subroutine SetErrorFile


  !!============================================================================
  !> @brief Returns the file unit number to be used for messages.
  !> @param[in] whichOne Which type of message
  !> - = 0 : Console (stderr)
  !> - > 0 : Error file (default)

  integer function GetErrorFile (whichOne)

    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_isTrue

    integer, intent(in), optional :: whichOne

    if (.not. present(whichOne)) then
       GetErrorFile = errorFile ! Return the errorFile unit number
    else if (whichOne > 0) then
       GetErrorFile = errorFile ! Return the errorFile unit number
    else if (ffa_cmdlinearg_isTrue('consolemsg')) then
       GetErrorFile =  6 ! Return the unit number for console output
    else
       GetErrorFile = -1 ! Suppress all error messages to console
    end if

  end function GetErrorFile


  !!============================================================================
  !> @brief Writes an error message to the console, and to the error file.
  !> @details The error file is either \a lpu (if present) or \a errorFile.
  !> The error message may be divided into (up to) four lines.

  subroutine ReportError (errCode,string1,string2,string3,string4, &
       &                  lpu,ierr,addString)

    integer     , intent(in)           :: errCode
    character(*), intent(in)           :: string1
    character(*), intent(in), optional :: string2, string3, string4, addString
    integer     , intent(in), optional :: lpu, ierr

    integer :: lpuErr
#if defined(irix) || defined(irix64)
    integer :: istat
#endif

    if (errCode >= 0) then
       !! Write error message to console
       lpuErr = GetErrorFile(0)
       if (lpuErr >= 0) then
          write(lpuErr,100) preString_p(errCode), trim(string1)
          if (present(string2)) write(lpuErr,101) trim(string2)
          if (present(string3)) write(lpuErr,101) trim(string3)
          if (present(string4)) write(lpuErr,101) trim(string4)
       end if
    end if

    !! Write error message to unit lpu if present, otherwise write to errorFile
    if (present(lpu)) then
       lpuErr = lpu
    else
       lpuErr = errorFile
    end if

    if (present(ierr)) then
       write(lpuErr,110) preString_p(abs(errCode)), trim(string1), ierr
    else
       write(lpuErr,100) preString_p(abs(errCode)), trim(string1)
    end if
    if (present(string2)) write(lpuErr,101) trim(string2)
    if (present(string3)) write(lpuErr,101) trim(string3)
    if (present(string4)) write(lpuErr,101) trim(string4)

    !! If present, write the additional string to the error file
    if (present(addString)) write(lpuErr,100) preString_p(0),trim(addString)
#if defined(irix) || defined(irix64)
    call flush(lpuErr,istat)
#else
    call flush(lpuErr)
#endif

100 format(A10,A)
101 format(10X,A)
110 format(A10,A,', IERR =',I6)

  end subroutine ReportError


  !!============================================================================
  !> @brief Opens a file for terminal output.

  subroutine OpenTerminalOutputFile (iunit)

    use FileUtilitiesModule   , only : findUnitNumber
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getint

    integer, intent(out) :: iunit !< File unit number for console output

    integer :: ierr

    call ffa_cmdlinearg_getint ('terminal',iunit)
    if (iunit < 0) then
       !! Terminal output is not wanted, so write to a scratch file
       !! which automatically will be deleted upon program completion
       iunit = findUnitNumber(8)
       open(iunit,STATUS='SCRATCH',IOSTAT=ierr)
       if (ierr /= 0) iunit = 0
    else if (iunit > 6) then
       !! Write terminal output to 'progress_info.res' in current directory
       iunit = findUnitNumber(iunit)
       open(iunit,FILE='progress_info.res',STATUS='UNKNOWN',IOSTAT=ierr)
       if (ierr /= 0) iunit = 0
    end if

  end subroutine OpenTerminalOutputFile


  !!============================================================================
  !> @brief Generates a memory allocation error message.
  !> @param[in] errorString Internal identification, e.g., a subroutine name

  integer function AllocationError (errorString)

    character(len=*), intent(in) :: errorString

    call ReportError (error_p,'Unable to allocate memory.', &
         'The problem may be too big to be solved on your computer. Free up'// &
         ' some memory by','shutting down applications and then try again.' // &
         ' If not already in use, you may also try','the -Bramsize option'  // &
         ' to reduce the memory need (reducer and recovery processes only).',  &
         addString=errorString)

    AllocationError = -1

  end function AllocationError


  !!============================================================================
  !> @brief Generates an internal error message.
  !> @param[in] errorString Internal identification, e.g., a subroutine name
  !>
  !> @details This function shall be used only to mark logic program errors
  !> that normally should not happen (bugs), and when no other sensible
  !> feedback message can be given to the user. When compiled in debug mode,
  !> the function triggers a segmentation fault to facilitate debugging.

  integer function InternalError (errorString)

    character(len=*), intent(in) :: errorString

    call ReportError (error_p,'Internal error', &
         'You may have encountered a bug in the program', &
         'Please give details to developers@openfedem.org', &
         addString=errorString)

    InternalError = -1

#ifdef FT_DEBUG
#if defined(win32) || defined(win64)
    call abort('due to internal error')
#else
    call flush() ! Flush all open files before aborting
    call abort()
#endif
#endif

  end function InternalError

  !> @endcond
end module ReportErrorModule
