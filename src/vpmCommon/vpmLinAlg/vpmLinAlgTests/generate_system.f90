!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!!==============================================================================
!> @brief Utilities for generation of simple FE mesh data structures (SAM).
!>
!> @author Knut Morten Okstad, SAP SE
!>
!> @date 16 Aug 2021

!!==============================================================================
!> @brief Generates a simple structured FE model for unit testing.
!>
!> @param[in] ipart Model identifier (see testModels.C)
!> @param[in] nel Number of elements in each direction
!> @param[in] geom Geometry parameters (length, width, etc.)
!> @param[out] sam Data for managing system matrix assembly
!> @param[in] ipsw Print switch for debug output
!> @param[in] lpu File unit number for debug output
!> @param[out] ierr Error flag

subroutine generate_model (ipart,nel,geom,sam,ipsw,lpu,ierr)

  use SamModule       , only : SamType, dp
  use TimerModule     , only : initTime
  use SamReducerModule, only : initiateSAM
  use AllocationModule, only : reAllocate

  implicit none

  integer      , intent(in)  :: ipart, nel(2)
  real(dp)     , intent(in)  :: geom(3)
  type(SamType), intent(out) :: sam
  integer      , intent(in)  :: ipsw, lpu
  integer      , intent(out) :: ierr

  !! Local variables
  logical , parameter :: linearSolve = .true.
  integer , parameter :: partId = 1, ngen = 10
  real(dp), pointer   :: tncoor(:,:) => null()
  integer             :: i, ndim, nenod

  !! --- Logic section ---

  call make_femodel (ipart,nel(1),nel(2),geom(1),geom(2),geom(3),ierr)
  if (ierr < 0) then
     write(lpu,"(' *** Failure generating FE model, ierr =',I6)") ierr
     return
  end if

  call initTime ()
  call initiateSAM (sam,tncoor,linearSolve,ipsw,lpu,ierr)
  call reAllocate ('generate_model',tncoor)
  if (ierr /= 0) return

  nenod = 0
  ndim = ngen
  do i = 1, sam%nnod
     if (sam%mnnn(i) > 1) then
        nenod = nenod + 1
        ndim = ndim + sam%madof(i+1)-sam%madof(i)
     end if
  end do

  sam%mpar(18) = partId
  sam%mpar(19) = ngen
  sam%mpar(22) = ngen
  sam%mpar(24) = ndim
  sam%mpar(29) = nenod

  write(lpu,1) sam%nel,sam%nnod,sam%mpar(29),sam%mpar(22), &
       &       sam%mpar(17),sam%mpar(27),sam%mpar(21),sam%mpar(24)

1 format(//4X,'FE MODEL SUMMARY' / 4X,16('-') &
       & //4X,'Number of elements             =',I8 &
       &  /4X,'Number of nodes (total)        =',I8 &
       &  /4X,'Number of external nodes       =',I8 &
       &  /4X,'Number of generalized modes    =',I8 &
       &  /4X,'Number of materials            =',I8 &
       &  /4X,'Number of geometric properties =',I8 &
       &  /4X,'Number of load cases           =',I8 &
       &  /4X,'Dimension of reduced matrices  =',I8 / )

end subroutine generate_model


!!==============================================================================
!> @brief Preprocesses the FE data structures.
!>
!> @param sam Data for managing system matrix assembly
!> @param[out] sysK System stiffness matrix
!> @param[out] sysM System mass matrix
!> @param[in] lumpedMass If .true., use a lumped mass matrix formulation
!> @param[in] ipsw Print switch for debug output
!> @param[in] lpu File unit number for debug output
!> @param[out] ierr Error flag

subroutine preproc_system (sam,sysK,sysM,lumpedMass,ipsw,lpu,ierr)

  use SamModule          , only : SamType, writeObject, writeSamSize
  use SysMatrixTypeModule, only : SysMatrixType
  use SysMatrixTypeModule, only : sparseMatrix_p, diagonalMatrix_p
  use SysMatrixTypeModule, only : nullifySysMatrix, shareMatrixStructure
  use SysMatrixTypeModule, only : check, writeObject
  use AsmExtensionModule , only : castToInt8, csAllocPointerArrays
  use AsmExtensionModule , only : csRenumber, csPreassemble
  use AsmExtensionModule , only : csBeginAssembly, csEndAssembly
  use ProgressModule     , only : writeProgress
  use ReportErrorModule  , only : reportError, error_p, debugFileOnly_p

  implicit none

  type(SamType)      , intent(inout) :: sam
  type(SysMatrixType), intent(out)   :: sysK, sysM
  logical            , intent(in)    :: lumpedMass
  integer            , intent(in)    :: ipsw, lpu
  integer            , intent(out)   :: ierr

  !! Local variables
  logical, parameter :: linearSolve = .true., factorMass = .false.
  character(len=128) :: errMsg

  !! --- Logic section ---

  call writeProgress (' --> Creating data structures for equation solver')
  call castToInt8 (sam,ierr)
  if (ierr < 0) goto 990

  call csAllocPointerArrays (sam,sysK,sparseMatrix_p,lpu,ierr,ipsw)
  if (ierr < 0) goto 991

  call csRenumber (sam,sysK,lpu,ierr)
  if (ierr < 0) goto 991

  call writeProgress (' --> Pre-assembling the stiffness matrix')
  call csPreassemble (sam,sysK,.not.linearSolve,lpu,ierr,ipsw)
  if (ierr < 0) goto 991

  call check (sysK,sam%mpar,lpu,ierr)
  if (ierr < 0) goto 991

  if (sam%mpar(19) < 1 .and. sam%mpar(30) < 1) then

     ! No system mass matrix needed
     call nullifySysMatrix (sysM)

  else if (.not. lumpedMass) then

     ! Let the mass- and stiffness matrices share datastructures
     call shareMatrixStructure (sysM,sysK,factorMass,ierr=ierr)

  else if (sam%nmmceq > 0 .and. lumpedMass) then

     ! Lumped mass is requested for a FE model that has constraint equations.
     ! The system mass matrix will then have a different sparsity pattern
     ! than the system stiffness matrix.
     call csAllocPointerArrays (sam,sysM,-sysK%storageType,lpu,ierr,ipsw)
     if (ierr < 0) goto 900

     call csRenumber (sam,sysM,lpu,ierr)
     if (ierr < 0) goto 900

     call writeProgress (' --> Pre-assembling the lumped mass matrix')
     call csPreassemble (sam,sysM,.false.,lpu,ierr,ipsw,sam%meqn)
     if (ierr < 0) goto 900

     call check (sysM,sam%mpar,lpu,ierr)

  else

     ! Lumped mass is requested for a FE model with no constraint equations.
     ! The system mass matrix will then be a pure diagonal matrix.
     call shareMatrixStructure (sysM,sysK,factorMass,diagonalMatrix_p)

  end if

  call writeSAMsize (sam,lpu)

  if (sam%mpar(22) > sam%ndof1) then
     ierr = sam%ndof1 - sam%mpar(22)
     write(errMsg,690) sam%mpar(22), sam%ndof1
690  format('NGEN =',I8,',  NDOF1 =',I8)
     call reportError (error_p,'The number of generalized DOFs (NGEN)', &
          'is larger than the number of internal nodal DOFs (NDOF1).', &
          errMsg,'This is not allowed, reduce NGEN to less than NDOF1.')
  end if

  if (ierr >= 0) call csBeginAssembly (sysK,'Stiffness matrix',lpu,ierr)
  if (ierr >= 0) call csBeginAssembly (sysM,'Mass matrix',lpu,ierr)

900 continue
  if (ipsw > 1 .and. ierr <= 0) then
     call writeObject (sam,lpu,4)
     call writeObject (sysK,lpu,' === System stiffness matrix',10,3)
     call writeObject (sysM,lpu,' === System mass matrix',10,3)
  else if (ipsw >= -7 .and. ipsw <= -5 .and. ierr == 0) then
     call writeObject (sam,lpu,5)
  end if
  if (ierr == 0) then
     call writeProgress (' --> Done preprocessing.')
  else
     call reportError (debugFileOnly_p,'preproc_system')
  end if
  return

990 call nullifySysMatrix (sysK)
991 call nullifySysMatrix (sysM)
  goto 900

end subroutine preproc_system
