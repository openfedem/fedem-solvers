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
#ifdef __GNUC__
  use FEModel, only : ourSams
#endif

  implicit none

  type(FEDataInput), intent(in)  :: T
  integer(is)      , intent(out) :: nels, nnpc
  integer(is)      , pointer     :: xnpc(:), npc(:)

#ifdef __GNUC__
  nels =  ourSams(T%idx)%p%nel
  nnpc =  ourSams(T%idx)%p%nmmnpc
  xnpc => ourSams(T%idx)%p%mpmnpc
  npc  => ourSams(T%idx)%p%mmnpc
#else
  nels =  T%sam%nel
  nnpc =  T%sam%nmmnpc
  xnpc => T%sam%mpmnpc
  npc  => T%sam%mmnpc
#endif

end subroutine GetConnectivityFEData


!!=============================================================================
!> @brief Extracts the constraint equations of the linear couplings.

subroutine GetConstraintsFEData (T, nceq, nnceq, xceq, ceq, tcc)
  !DEC$ ATTRIBUTES DLLEXPORT :: GetConstraintsFEData

  use FEKind, only : is, wp
  use FEData, only : FEDataInput
#ifdef __GNUC__
  use FEModel, only : ourSams
#endif

  implicit none

  type(FEDataInput), intent(in)        :: T
  integer(is)      , intent(out)       :: nceq, nnceq
  integer(is)      , pointer           :: xceq(:), ceq(:)
  real(wp)         , pointer, optional :: tcc(:)

#ifdef __GNUC__
  nceq  =  ourSams(T%idx)%p%nceq
  nnceq =  ourSams(T%idx)%p%nmmceq
  xceq  => ourSams(T%idx)%p%mpmceq
  ceq   => ourSams(T%idx)%p%mmceq
  if (present(tcc)) then
     tcc => ourSams(T%idx)%p%ttcc
  endif
#else
  nceq  =  T%sam%nceq
  nnceq =  T%sam%nmmceq
  xceq  => T%sam%mpmceq
  ceq   => T%sam%mmceq
  if (present(tcc)) then
     tcc => T%sam%ttcc
  end if
#endif

end subroutine GetConstraintsFEData


!!=============================================================================
!> @brief Extracts the nodal DOF information.

subroutine GetPartitionFEData (T, nnod, ndof, madof, msc, minex)
  !DEC$ ATTRIBUTES DLLEXPORT :: GetPartitionFEData

  use FEKind, only : is
  use FEData, only : FEDataInput
#ifdef __GNUC__
  use FEModel, only : ourSams
#endif

  implicit none

  type(FEDataInput), intent(in)  :: T
  integer(is)      , intent(out) :: nnod, ndof
  integer(is)      , pointer     :: madof(:), msc(:), minex(:)

#ifdef __GNUC__
  nnod  =  ourSams(T%idx)%p%nnod
  ndof  =  ourSams(T%idx)%p%ndof
  madof => ourSams(T%idx)%p%madof
  msc   => ourSams(T%idx)%p%msc
  minex => ourSams(T%idx)%p%minex
#else
  nnod  =  T%sam%nnod
  ndof  =  T%sam%ndof
  madof => T%sam%madof
  msc   => T%sam%msc
  minex => T%sam%minex
#endif

end subroutine GetPartitionFEData
