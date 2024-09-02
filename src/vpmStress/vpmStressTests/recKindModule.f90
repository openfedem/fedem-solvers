!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file recKindModule.f90
!> @brief Real kind parameter to use in the recovery calculations.

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
