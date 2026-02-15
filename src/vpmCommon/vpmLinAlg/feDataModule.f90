!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file feDataModule.f90
!> @brief FE data interface for the GSF equation solver.

!!==============================================================================
!> @brief FE data interface for the GSF equation solver.
!> @details This module defines a data type that is used to transport the
!> the FE data from FEDEM into the DNVS _Linear Algebra Modules_ subroutines.
!> @note The module *must* be named FEData and the data type FEDataInput,
!> as these names are also used within the DNVS library. However, the contents
!> of the data type is arbitrary and can be application specific.

module FEData

  use SamModule, only : SamType

  implicit none

  type FEDataInput
#ifdef __GNUC__
     integer :: idx = 0 !< Index to the SAM data instance
#else
     type(SamType), pointer :: sam => null() !< Pointer to the SAM data instance
#endif
  end type FEDataInput

end module FEData
