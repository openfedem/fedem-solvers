!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file dbgUnitsModule.f90
!> @brief File unit numbers for debug print.

!!==============================================================================
!> @brief Module with file unit numbers for debug print.
!> @details This module contains all file unit numbers used for debug print of
!> various tasks. They are gathered in a separate module to better control their
!> deallocation on program termination, without actual termination of the
!> running process itself.

module dbgUnitsModule

  use KindModule, only : dp

  implicit none

  integer, save :: dbgCorot = 0 !< For corotational update calculations
  integer, save :: dbgJoint = 0 !< For joint update calculations
  integer, save :: dbgCurve = 0 !< For contact curve calculations
  integer, save :: dbgFric  = 0 !< For friction calculations
  integer, save :: dbgSolve = 0 !< For general system level calculations
  integer, save :: dbgTire = 0  !< For tire update calculations
  integer, save :: dbgRoad = 0  !< For road evaluations
  integer, save :: dbgRot = 0   !< For mass matrix correction calculations
  integer, save :: dbgUDE = 0   !< For user-defined element calculations
  integer, save :: dbgAD = 0    !< For aerodynamc load calculations
  integer, save :: dbgHD = 0    !< For hydrodynamc load calculations

  integer , save :: dbgIter = 0      !< Current iteration counter
  real(dp), save :: dbgTime = 0.0_dp !< Current physical time


contains

  !!============================================================================
  !> @brief Closes all debug print files.
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 1 Feb 2017

  subroutine closeDbgUnits

    use FileUtilitiesModule, only : closeDBGfiles

    call closeDBGfiles ()
    call closeUnit (dbgCorot)
    call closeUnit (dbgJoint)
    call closeUnit (dbgCurve)
    call closeUnit (dbgFric)
    call closeUnit (dbgSolve)
    call closeUnit (dbgTire)
    call closeUnit (dbgRoad)
    call closeUnit (dbgRot)
    call closeUnit (dbgUDE)
    call closeUnit (dbgAD)
    call closeUnit (dbgHD)

  contains
    !> @cond NO_DOCUMENTAION
    subroutine closeUnit (iunit)
      integer, intent(inout) :: iunit
      if (iunit > 0) then
         close(iunit)
         iunit = 0
      end if
    end subroutine closeUnit
    !> @endcond

  end subroutine closeDBGUnits

end module dbgUnitsModule
