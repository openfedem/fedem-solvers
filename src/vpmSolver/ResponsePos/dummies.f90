!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!! Dummy subroutines

subroutine CTIO (ID, NVAR, VAR)
  INTEGER ID, NVAR
  DOUBLE PRECISION VAR(NVAR)
  print*,' ** CTIO dummy: ',ID,NVAR,VAR(1)
end subroutine CTIO

function checkCtrlParams (elmType,nRealData,nVar)
  integer, intent(inout) :: elmType
  integer, intent(in)    :: nRealData, nVar
  logical :: checkCtrlParams
  print*,' ** checkCtrlParams dummy: ',elmType,nRealData,nVar
  checkCtrlParams = .false.
end function checkCtrlParams

module FFlLinkHandlerInterface; contains
  subroutine ffl_getnodalcoor (x,y,z,inod,ierr)
    use kindModule, only : dp
    real(dp) :: x, y, z
    integer  :: inod, ierr
    ierr = -1
    print*,' ** ffl_getnodalcoor dummy: ',x,y,z,inod
  end subroutine ffl_getnodalcoor
end module FFlLinkHandlerInterface

module FiUserElmInterface; contains
  subroutine fi_ude4 (eId,eType,idx,iwork,rwork,name,nvar)
    use kindModule, only : dp
    integer         , intent(in)  :: eId, eType, idx, iwork
    real(dp)        , intent(in)  :: rwork
    character(len=*), intent(out) :: name
    integer         , intent(out) :: nvar
    name = ''
    nvar = 0
    print*,' ** fi_ude4 dummy: ',eId,eType,idx,iwork,rwork
  end subroutine fi_ude4
  subroutine fi_ude5 (eId,eType,idx,iwork,rwork,value,nvar)
    use kindModule, only : dp
    integer , intent(in)  :: eId, eType, idx, iwork
    real(dp), intent(in)  :: rwork
    real(dp), intent(out) :: value
    integer , intent(out) :: nvar
    value = 0.0_dp
    nvar = 0
    print*,' ** fi_ude5 dummy: ',eId,eType,idx,iwork,rwork
  end subroutine fi_ude5
end module FiUserElmInterface

module hydroDynamicsModule; contains
  function getCalculatedFluidMotion (triad,ele,vel,acc)
    use TriadTypeModule, only : TriadType, dp
    logical :: getCalculatedFluidMotion
    type(TriadType)   , intent(in) :: triad
    real(dp),optional , intent(out) :: ele, vel(3), acc(3)
    getCalculatedFluidMotion = .false.
    if (present(ele)) ele = 0.0_dp
    if (present(vel)) vel = 0.0_dp
    if (present(acc)) acc = 0.0_dp
    print*,' ** getCalculatedFluidMotion dummy: ',triad%id%baseId
  end function getCalculatedFluidMotion
end module hydroDynamicsModule

module finiteElementModule; contains
  function GetSectionalForces (sup,iend)
    use SupElTypeModule, only : SupElType, dp
    type(SupElType), intent(in) :: sup
    integer, intent(in) :: iend
    real(dp) :: GetSectionalForces(6)
    GetSectionalForces = 0.0_dp
    print*,' ** getSectionalForces dummy: ',sup%id%baseId,iend
  end function GetSectionalForces
end module finiteElementModule

module SysMatrixTypeModule
  implicit none
  type SysMatrixType
     integer :: dim
  end type SysMatrixType
  interface writeObject
    module procedure writeSysMat
  end interface
  private :: writeSysMat
contains
  subroutine nullifySysMatrix (this)
    type(SysMatrixType), intent(out) :: this
    this%dim = 0
  end subroutine nullifySysMatrix
  subroutine deallocateSysMatrix (this)
    type(SysMatrixType), intent(inout) :: this
    print*,' ** deallocateSysMat dummy: ',this%dim
  end subroutine deallocateSysMatrix
  subroutine writeSysMat (this,io,text,complexity)
    type(SysMatrixType), intent(in)           :: this
    integer            , intent(in)           :: io
    character(len=*)   , intent(in), optional :: text
    integer            , intent(in), optional :: complexity
    print*,' ** writeSysMat dummy: ',this%dim,io
    if (present(complexity) .and. present(text)) print*,text
  end subroutine writeSysMat
end module SysMatrixTypeModule

module RecKindModule
  use kindModule, only : sp
  integer, parameter, public :: rk = sp
end module RecKindModule
