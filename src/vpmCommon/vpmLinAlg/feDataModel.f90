!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file feDataModel.f90
!> @brief SAM data container for the system matrices.

!!==============================================================================
!> @brief SAM data container for the system matrices.
!> @details This module is used to handle a global container of SAM data for
!> the (up to) two system matrices of the FEDEM model. This container is then
!> referred by the call-back subroutines of the GSF equation solver,
!> to extract the FE data required to build its internal data structure.

module FEModel

  use SamModule, only : SamType

  !> @brief Data type representing a pointer to a sammodule::samtype instance.
  type SamPtrType
    type(SamType), pointer :: p => null() !< Pointer to the SAM instance
  end type SamPtrType

  type(SamPtrType), save :: ourSams(2) !< Global SAM data container

end module FEModel
