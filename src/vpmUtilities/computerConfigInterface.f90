!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file computerConfigInterface.f90
!> @brief Fortran interface for functions for extracting computer configuration.
!> @details This file contains a module with Fortran interface definitions for
!> global functions for extraction of the computer configuration and user name.
!> See the file computerConfig.C for the actual implementation.

!!==============================================================================
!> @brief Fortran interface for functions for extracting computer configuration.

module computerConfigInterface

  implicit none

  interface

     !> @brief Returns a string identifying the computer the program is run on.
     subroutine getComputerConfig (cid)
       character, intent(out) :: cid*(*)
     end subroutine getComputerConfig

     !> @brief Returns the user name of the running process owner.
     subroutine getUserName (cuser)
       character, intent(out) :: cuser*(*)
     end subroutine getUserName

  end interface

end module computerConfigInterface
