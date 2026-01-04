!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file feDataRoutines.f90
!> @brief Call-backs for the GSF equation solver.

!> @cond FULL_DOC
!!==============================================================================
!> @brief Module with kind-parameters for the GSF equation solver call-backs.
!> @details The two kind parameters defined in this module have to match the
!> corresponding parameters of the kind_values module of the GSF library.

module FEKind
  integer, parameter :: is = SELECTED_INT_KIND(9)
  integer, parameter :: wp = SELECTED_REAL_KIND(15,307)
end module FEKind
!> @endcond


!!=============================================================================
!> @brief Extracts the element connectivity table for the FE model.

subroutine GetConnectivityFEData (T, nels, nnpc, xnpc, npc)
  !DEC$ ATTRIBUTES DLLEXPORT :: GetConnectivityFEData

  use FEKind, only : is
  use FEData, only : FEDataInput

  implicit none

  type(FEDataInput), intent(in)  :: T
  integer(is)      , intent(out) :: nels, nnpc
  integer(is)      , pointer     :: xnpc(:), npc(:)

  nels =  T%sam%nel
  nnpc =  T%sam%nmmnpc
  xnpc => T%sam%mpmnpc
  npc  => T%sam%mmnpc

end subroutine GetConnectivityFEData


!!=============================================================================
!> @brief Extracts the constraint equations of the linear couplings.

subroutine GetConstraintsFEData (T, nceq, nnceq, xceq, ceq, tcc)
  !DEC$ ATTRIBUTES DLLEXPORT :: GetConstraintsFEData

  use FEKind, only : is, wp
  use FEData, only : FEDataInput

  implicit none

  type(FEDataInput), intent(in)        :: T
  integer(is)      , intent(out)       :: nceq, nnceq
  integer(is)      , pointer           :: xceq(:), ceq(:)
  real(wp)         , pointer, optional :: tcc(:)

  nceq  =  T%sam%nceq
  nnceq =  T%sam%nmmceq
  xceq  => T%sam%mpmceq
  ceq   => T%sam%mmceq
  if (present(tcc)) then
     tcc => T%sam%ttcc
  end if

end subroutine GetConstraintsFEData


!!=============================================================================
!> @brief Extracts the nodal DOF information.

subroutine GetPartitionFEData (T, nnod, ndof, madof, msc, minex)
  !DEC$ ATTRIBUTES DLLEXPORT :: GetPartitionFEData

  use FEKind, only : is
  use FEData, only : FEDataInput

  implicit none

  type(FEDataInput), intent(in)  :: T
  integer(is)      , intent(out) :: nnod, ndof
  integer(is)      , pointer     :: madof(:), msc(:), minex(:)

  nnod  =  T%sam%nnod
  ndof  =  T%sam%ndof
  madof => T%sam%madof
  msc   => T%sam%msc
  minex => T%sam%minex

end subroutine GetPartitionFEData
