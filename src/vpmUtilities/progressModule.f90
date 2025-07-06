!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file progressModule.f90
!> @brief Subroutines for reporting progress info during computations.

!!==============================================================================
!> @brief Module with subroutines for reporting progress during computations.
!> @details This module contains various utility subroutines and functions
!> for printing progress messages to the console during lengthy computations.
!> You have to configure doxygen with the option ENABLED_SECTIONS = FULL_DOC
!> to extract the detailed documentation of the symbols in this module.
!>
!> @author Knut Morten Okstad
!> @date 9 Feb 2001

module ProgressModule
  !> @cond FULL_DOC

  use KindModule, only : sp

  implicit none

  integer , save          :: lterm   !< File unit number for terminal output
  real(sp), save, private :: last(2) !< Used to log incremental time consumption

  interface writeProgress
     module procedure writeMessage
     module procedure writeProsent
     module procedure writeProsentDouble
  end interface

  private :: getTimeUsage, writeMessage, writeProsent, writeProsentDouble


contains

  !!============================================================================
  !> @brief Returns incremental time consumption since the previous invokation
  !> @param[out] time Time measurements
  !> - time(1) = wall time since previous invokation of this function
  !> - time(2) = percentage CPU-time usage in time(1)
  !> - time(3) = wall time since the start of the program
  !> @param[in] hundreds If present, also return times smaller that 0.01s
  !>
  !> @author Knut Morten Okstad
  !> @date 11 Sep 2006

  function getTimeUsage (time,hundreds) result(changed)

    logical               :: changed
    real(sp), intent(out) :: time(3)
    logical , intent(in)  :: hundreds
    real(sp), external    :: CPUSEC, CLKSEC

    !! Local variables
    real(sp) :: currentCPU, currentCLK

    !! --- Logic section ---

    currentCPU = CPUSEC(0.0_sp)
    currentCLK = CLKSEC(0.0_sp)
    time(1) = currentCLK - last(2)
    time(3) = currentCLK
    if (time(1) > 1.0_sp .or. (time(1) >= 0.005_sp .and. hundreds)) then
       time(2) = 100.0_sp * (currentCPU-last(1)) / time(1)
       last(1) = currentCPU
       last(2) = currentCLK
       changed = .true.
    else
       time(2) = 0.0_sp
       changed = .false.
    end if

  end function getTimeUsage


  !!============================================================================
  !> @brief Writes progress information for a do-loop to the terminal.
  !> @details A message is written for each 1% completion of the loop,
  !> as long as the wall time usage since the last message is 1 sec or more.
  !> A message is written for each 10% completion of the loop, as long as the
  !> wall time usage since the last message is 0.01 sec or more.

  subroutine writeProsent (loopVar,loopEnd)

    integer, intent(in) :: loopVar, loopEnd

    !! Local variables
    real(sp) :: time(3), prosent, lastProsent = 1.0_sp, lastWrite = 0.0_sp
    logical  :: endLoop

    !! --- Logic section ---

    if (loopEnd <= 0) return

    prosent = 100.0_sp * (real(loopVar,sp))/(real(loopEnd,sp))
    endLoop = loopVar >= loopEnd

    if (prosent >= lastProsent .or. endLoop) then
       if (getTimeUsage(time,prosent >= lastWrite+10.0_sp .or. endLoop)) then
          if (time(2) < 99.5_sp) then
             write(lterm,600) nint(prosent), time(1), nint(time(2))
          else
             write(lterm,601) nint(prosent), time(1)
          end if
          call flush (lterm)
          lastWrite   = real(nint(prosent),sp)
          lastProsent = lastWrite + 1.0_sp
       else
          lastProsent = lastProsent + 1.0_sp
       end if
       if (endLoop) then
          lastWrite   = 0.0_sp
          lastProsent = 1.0_sp
       end if
    end if

600 format(I20,'% completed  =>',F8.2,I4,'% CPU time')
601 format(I20,'% completed  =>',F8.2)

  end subroutine writeProsent


  !!============================================================================
  !> @brief Writes progress information for a do-loop to the terminal.
  !> @details Same as writeProsent except that the loop variables now are
  !> double precision reals.

  subroutine writeProsentDouble (current,total)

    use kindModule, only : dp, epsDiv0_p

    real(dp), intent(in) :: current, total

    !! Local variables
    real(sp) :: time(3), prosent, lastProsent = 1.0_sp, lastWrite = 0.0_sp
    logical  :: endLoop

    !! --- Logic section ---

    if (total <= epsDiv0_p) return

    prosent = 100.0_sp * real(current/total,sp)
    endLoop = current >= total

    if (prosent >= lastProsent .or. endLoop) then
       if (getTimeUsage(time,prosent >= lastWrite+10.0_sp .or. endLoop)) then
          if (time(2) < 99.5_sp) then
             write(lterm,600) nint(prosent), time(1), nint(time(2))
          else
             write(lterm,601) nint(prosent), time(1)
          end if
          call flush (lterm)
          lastWrite   = real(nint(prosent),sp)
          lastProsent = lastWrite + 1.0_sp
       else
          lastProsent = lastProsent + 1.0_sp
       end if
       if (endLoop) then
          lastWrite   = 0.0_sp
          lastProsent = 1.0_sp
       end if
    end if

600 format(I20,'% completed  =>',F8.2,I4,'% CPU time')
601 format(I20,'% completed  =>',F8.2)

  end subroutine writeProsentDouble


  !!============================================================================
  !> @brief Write a progress message with optional memory usage to the terminal.
  !>

  subroutine writeMessage (msg,reportMemoryUsage,newLine)

    use FFaProfilerInterface, only : ffa_getMemUsage

    character(len=*), intent(in) :: msg
    logical,optional, intent(in) :: reportMemoryUsage, newLine

    !! Local variables
    real(sp) :: time(3), usage(4)
    logical  :: firstCall = .true.

    !! --- Logic section ---

    if (firstCall) then
       last = 0.0_sp
       firstCall = .false.
       write(lterm,600)
    end if

    if (.not. present(newLine)) then
       write(lterm,"()")
    else if (newLine) then
       write(lterm,"()")
    end if

    if (getTimeUsage(time,.true.)) then
       if (time(2) < 99.5_sp) then
          write(lterm,601) msg, time(1), nint(time(2)), time(3)
       else
          write(lterm,602) msg, time(1), time(3)
       end if
    else
       write(lterm,603) msg, time(3)
    end if

    if (present(reportMemoryUsage)) then
       if (reportMemoryUsage) then
          call ffa_getMemUsage (usage)
          if (usage(3) > 0.0) write(lterm,610) 'Work',usage(1),usage(1)/usage(3)
          if (usage(4) > 0.0) write(lterm,610) 'Page',usage(2),usage(2)/usage(4)
       end if
    end if
    call flush (lterm)

600 format(/T57,'Incremental Rel. Total'/T57,' Wall time  CPU  Wall t')
601 format(2X,A,T57,'=>',F8.2,I4,'%',F8.1)
602 format(2X,A,T57,'=>',F8.2,5X,F8.1)
603 format(2X,A,T57,'=>',13X,F8.1)
610 format(7X,A,' memory:',F9.2,' MB',2PF7.2,'% of peak')

  end subroutine writeMessage

  !> @endcond
end module ProgressModule
