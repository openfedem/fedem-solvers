!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file asmExtensionModule.f90
!> @brief System matrix assembly extension module.

#if defined(FT_HAS_INT4_ONLY) || defined(FT_USE_INT8)
#define FT_SPR_INT_DEFAULT 1
#elif !defined(FT_HAS_SPR_INT8)
#define FT_SPR_INT_DEFAULT 1
#endif

!!==============================================================================
!> @brief System matrix assembly extension module.
!>
!> @details This module contains subroutines that administer the preprocessing
!> steps of a linear system of equations, and the assembly of element matrices
!> and associated right-hand-side vectors into system matrices and vectors.
!> The module defines an interface between the FEDEM solver applications
!> and the various linear algebra packages that are being used.

module AsmExtensionModule

  use sprKindModule, only : ik

  implicit none

  !> Data type for 64-bit version of some integer arrays in sammodule::samtype.
  type Sam8Type
     !> Copy of sammodule::samtype::mpar
     integer(ik), pointer :: mpar(:)   => null()
     !> Copy of sammodule::samtype::mpmnpc
     integer(ik), pointer :: mpmnpc(:) => null()
     !> Copy of sammodule::samtype::mmnpc
     integer(ik), pointer :: mmnpc(:)  => null()
     !> Copy of sammodule::samtype::madof
     integer(ik), pointer :: madof(:)  => null()
     !> Copy of sammodule::samtype::msc
     integer(ik), pointer :: msc(:)    => null()
     !> Copy of sammodule::samtype::mpmceq
     integer(ik), pointer :: mpmceq(:) => null()
     !> Copy of sammodule::samtype::mmceq
     integer(ik), pointer :: mmceq(:)  => null()
     !> Copy of sammodule::samtype::minex
     integer(ik), pointer :: minex(:)  => null()
  end type Sam8Type

  type(Sam8Type), private, save :: mySam !< 64-bit integer version of SAM arrays

  logical, private, save :: bigModel !< Flag for extra progress logging

  private :: csExpandElMat, csGetElmEq, writeMEQN
#ifdef FT_HAS_MKL
  private :: csPreAssSparse
#endif


contains

  !!============================================================================
  !> @brief Writes out the assigned equation numbers.
  !>
  !> @param[in] prefix Message prefix
  !> @param[in] madof Matrix of accumulated DOFs
  !> @param[in] meqn Matrix of equation numbers
  !> @param[in] lpu File unit number for res-file output
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 28 May 2018

  subroutine writeMEQN (prefix,madof,meqn,lpu)

    character(len=*), intent(in) :: prefix
    integer         , intent(in) :: madof(:), meqn(:), lpu

    integer :: i, ndof

    !! --- Logic section ---

    ndof = madof(size(madof))-1
    write(lpu,"(A)") prefix//':MEQN'
    do i = 1, size(madof)-1
       if (madof(i+1) > madof(i)) then
          if (ndof > 999999) then
             write(lpu,8) i, madof(i), madof(i+1)-1, meqn(madof(i):madof(i+1)-1)
          else if (ndof > 99999) then
             write(lpu,7) i, madof(i), madof(i+1)-1, meqn(madof(i):madof(i+1)-1)
          else if (ndof > 9999) then
             write(lpu,6) i, madof(i), madof(i+1)-1, meqn(madof(i):madof(i+1)-1)
          else if (ndof > 999) then
             write(lpu,5) i, madof(i), madof(i+1)-1, meqn(madof(i):madof(i+1)-1)
          else if (ndof > 99) then
             write(lpu,4) i, madof(i), madof(i+1)-1, meqn(madof(i):madof(i+1)-1)
          else
             write(lpu,3) i, madof(i), madof(i+1)-1, meqn(madof(i):madof(i+1)-1)
          end if
       end if
    end do

3   format(I3,',',I3,' -',I3,' :',6I3)
4   format(I4,',',I4,' -',I4,' :',6I4)
5   format(I5,',',I5,' -',I5,' :',6I5)
6   format(I6,',',I6,' -',I6,' :',6I6)
7   format(I7,',',I7,' -',I7,' :',6I7)
8   format(I8,',',I8,' -',I8,' :',6I8)

  end subroutine writeMEQN


  !!============================================================================
  !> @brief Creates a 64-bit integer version of the SAM-arrays needed by SPR.
  !>
  !> @param[in] samData Data for managing system matrix assembly
  !> @param[out] ierr Error flag
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 24 May 2017

  subroutine castToInt8 (samData,ierr)

    use sprKindModule   , only : castInt
    use SamModule       , only : SamType
    use allocationModule, only : reAllocate

    type(SamType), target, intent(in)  :: samData
    integer,     optional, intent(out) :: ierr

    character(len=30), parameter :: prgnam_p = 'AsmExtensionModule::castToInt8'

    !! --- Logic section ---

    if (present(ierr)) ierr = 0
#ifdef FT_SPR_INT_DEFAULT
    !! The integer kinds i4 and i8 are the same, and ik can be either of them,
    !! or the SPR equation solver uses the default (32bit) integer kind
    if (present(ierr)) then
       mySam%mpar   => samData%mpar
       mySam%mpmnpc => samData%mpmnpc
       mySam%mmnpc  => samData%mmnpc
       mySam%madof  => samData%madof
       mySam%msc    => samData%msc
       mySam%mpmceq => samData%mpmceq
       mySam%mmceq  => samData%mmceq
       mySam%minex  => samData%minex
    else
       nullify(mySam%mpar)
       nullify(mySam%mpmnpc)
       nullify(mySam%mmnpc)
       nullify(mySam%madof)
       nullify(mySam%msc)
       nullify(mySam%mpmceq)
       nullify(mySam%mmceq)
       nullify(mySam%minex)
    end if
#else
    !! The SPR equation solver uses 64-bit integer,
    !! so create duplicates of the SAM array which it needs.
    if (present(ierr)) then
       call reAllocate (prgnam_p,mySam%mpar  ,size(samData%mpar)  ,ierr)
       call reAllocate (prgnam_p,mySam%mpmnpc,size(samData%mpmnpc),ierr)
       call reAllocate (prgnam_p,mySam%mmnpc ,size(samData%mmnpc) ,ierr)
       call reAllocate (prgnam_p,mySam%madof ,size(samData%madof) ,ierr)
       call reAllocate (prgnam_p,mySam%msc   ,size(samData%msc)   ,ierr)
       call reAllocate (prgnam_p,mySam%mpmceq,size(samData%mpmceq),ierr)
       if (associated(samData%mmceq)) then
          call reAllocate (prgnam_p,mySam%mmceq,size(samData%mmceq),ierr)
       end if
       call reAllocate (prgnam_p,mySam%minex ,size(samData%minex) ,ierr)
       if (ierr < 0) return
       call castInt (size(samData%mpar)  ,samData%mpar  ,mySam%mpar)
       call castInt (size(samData%mpmnpc),samData%mpmnpc,mySam%mpmnpc)
       call castInt (size(samData%mmnpc) ,samData%mmnpc ,mySam%mmnpc)
       call castInt (size(samData%madof) ,samData%madof ,mySam%madof)
       call castInt (size(samData%msc)   ,samData%msc   ,mySam%msc)
       call castInt (size(samData%mpmceq),samData%mpmceq,mySam%mpmceq)
       if (associated(samData%mmceq)) then
          call castInt (size(samData%mmceq),samData%mmceq,mySam%mmceq)
       end if
       call castInt (size(samData%minex) ,samData%minex ,mySam%minex)
    else
       call reAllocate (prgnam_p,mySam%mpar  )
       call reAllocate (prgnam_p,mySam%mpmnpc)
       call reAllocate (prgnam_p,mySam%mmnpc )
       call reAllocate (prgnam_p,mySam%madof )
       call reAllocate (prgnam_p,mySam%msc   )
       call reAllocate (prgnam_p,mySam%mpmceq)
       call reAllocate (prgnam_p,mySam%mmceq )
       call reAllocate (prgnam_p,mySam%minex )
    end if
#endif

  end subroutine castToInt8


  !!============================================================================
  !> @brief Initializes pointer arrays for the equation solvers.
  !>
  !> @param samData Data for managing system matrix assembly
  !> @param[out] sysMat The system matrix to initialize data structure for
  !> @param[in] matrixType Type of system matrix data structure
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !> @param[in] ipsw Print switch
  !> @param[in] cacheSize Sparse matricx cache size
  !> @param[in] needsFactorization If .true., this matrix will be factorized
  !> @param[in] masterGSF GSF data structure of sibling system matrix
  !>
  !> @details The symbolic assembly stage is also performed
  !> for the sparse matrix equation solvers.
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !> @date Sep 1998
  !>
  !> @author Knut Morten Okstad
  !> @date Feb 2003

  subroutine csAllocPointerArrays (samData,sysMat,matrixType,lpu,ierr, &
       &                           ipsw,cacheSize,needsFactorization,masterGSF)

    use sprKindModule      , only : i4, castInt
    use SamModule          , only : SamType, LumpedSAM
    use SysMatrixTypeModule, only : SysMatrixType, GSFInfo, checkGSFInfo
    use SysMatrixTypeModule, only : GSFStorageType, SparseStorageType
    use SysMatrixTypeModule, only : matrixType_p, nspar_p, diagonalMatrix_p
    use SysMatrixTypeModule, only : skylineMatrix_p, sparseMatrix_p
    use SysMatrixTypeModule, only : denseMatrix_p, outOfCore_p, pardiso_p
    use SysMatrixTypeModule, only : reAllocate, nullifySysMatrix
    use allocationModule   , only : reAllocate
    use progressModule     , only : writeProgress
    use FELinearSolver     , only : FEInitialize, FEAnalyze, FEMatrixOrder
    use FELinearSolver     , only : FEGetActiveOrder, FEGetSuperelementOrder
    use FELinearSolver     , only : FEPermuteIndices, SetLinearSolverType
    use FELinearSolver     , only : SetMemorySize, SetOutputUnit
    use FELinearSolver     , only : SetScratchFileLocation, WriteMessageData
    use scratchArrayModule , only : getIntegerScratchArray
    use reportErrorModule  , only : note_p, error_p, reportError, internalError

    type(SamType)       , target  , intent(inout) :: samData
    type(SysMatrixType) ,           intent(out)   :: sysMat
    integer             ,           intent(in)    :: matrixType, lpu
    integer             ,           intent(out)   :: ierr
    integer             , optional, intent(in)    :: ipsw, cacheSize
    logical             , optional, intent(in)    :: needsFactorization
    type(GSFStorageType), optional, pointer       :: masterGSF

    !! Local variables
    logical            :: lumpedElMats, firstMatrix
    integer(ik)        :: lerr
    integer            :: nWork, n0
    integer, pointer   :: iWork(:)
    integer, parameter :: GSFFactOpt = 2 ! Standard block LDLt factorization
    character(len=256) :: tmpDir

    type(SparseStorageType), pointer :: sparse

    character(len=40), parameter :: prgnam_p = &
         'AsmExtensionModule::csAllocPointerArrays'

    !! --- Logic section ---

    call nullifySysMatrix (sysMat)

    if (samData%ndof <= 0) then
       ierr = internalError(prgnam_p//': NDOF is not set yet')
       return
    else
       ierr = 0
    end if

    if (matrixType == outOfCore_p + sparseMatrix_p) then
       sysMat%storageType = outOfCore_p
    else
       sysMat%storageType = abs(matrixType)
    end if
    lumpedElMats = matrixType < 0 .and. matrixType /= -pardiso_p
    firstMatrix = samData%neq == 0
    bigModel = samData%nnod > 100000

    if (lumpedElMats) then
       !! All element matrix contributions are lumped (diagonalized), and
       !! the only off-diagonal terms may result from constraint equations.
       !! This is allowed only for the second system matrix, the first system
       !! matrix defined must always be consistent.
       if (firstMatrix .or. .not.associated(samData%meqn)) then
          ierr = internalError(prgnam_p//': NEQ and/or MEQN are not set yet')
          return
       else if (size(samData%meqn) /= samData%ndof) then
          ierr = internalError(prgnam_p//': Invalid dimension on MEQN')
          return
       else if (sysMat%storageType <= skylineMatrix_p) then
          ierr = internalError(prgnam_p//': Lumping is not implemented for '// &
               &               matrixType_p(sysMat%storageType))
          return
       end if
    end if

    select case (sysMat%storageType)

    case (diagonalMatrix_p, denseMatrix_p)

       call reAllocate (prgnam_p,sysMat%meqn,samData%ndof,ierr)
       if (ierr < 0) return

    case (skylineMatrix_p)

       if (samData%neq > 0) then
          sysMat%dim = samData%neq  ! Use NEQ to allocate pointer arrays
       else
          sysMat%dim = samData%ndof ! Use NDOF to allocate pointer arrays
       end if

       call reAllocate (prgnam_p,sysMat%meqn,samData%ndof,ierr)
       call reAllocate (prgnam_p,sysMat%skyline,.true.,ierr)
       if (ierr < 0) return

       call reAllocate (prgnam_p,sysMat%skyline%msky,sysMat%dim,ierr)
       if (ierr < 0) return

    case (sparseMatrix_p)

       call reAllocate (prgnam_p,sysMat%meqn,samData%ndof,ierr)
#ifdef FT_SPR_INT_DEFAULT
       sysMat%meqn8 => sysMat%meqn
#else
       call reAllocate (prgnam_p,sysMat%meqn8,samData%ndof,ierr)
#endif
       call reAllocate (prgnam_p,sysMat%sparse,nspar_p,ierr)
       if (ierr < 0) return

       sparse => sysMat%sparse
       nWork = samData%nnod + 2*samData%ndof
       if (ik > i4) nWork = nWork*2 ! We are using 64-bit integers
       iWork => getIntegerScratchArray(nWork,ierr)
       if (ierr < 0) goto 999

       if (lumpedElMats) then

          !! All element matrix contributions are lumped (diagonalized), and
          !! the only off-diagonal terms may result from constraint equations.
          !! Create a separate sparse structure accounting for these only,
          !! using MEQN for the previously defined matrix.
          call ICOPY (samData%ndof,samData%meqn(1),1,sysMat%meqn(1),1)
          if (ik > i4) call castInt (samData%ndof,sysMat%meqn,sysMat%meqn8)
          !TODO,kmo: This MEQN should be before the MMD reordering, yes?
          call SPRSEL (nspar_p,mySam%mpar(1),sparse%mspar(1))

          if (mySam%mpar(20) == 0_ik) then
             !! Compute size of the largest element
             call SPRMXD (mySam%madof(1), mySam%mpmnpc(1), &
                  &       mySam%mmnpc(1), mySam%mpar(1))
          end if

       else

          !! Preprocessing for consistent element matrices
          n0 = count(mySam%msc == 0)
          if (bigModel) call writeProgress ('     SPRPRP ...',newline=.false.)
          call SPRPRP (mySam%madof(1), mySam%minex(1),  mySam%mpmnpc(1), &
               &       mySam%mmnpc(1), mySam%mpmceq(1), mySam%mmceq(1), &
               &       mySam%msc(1), nspar_p, int(lpu,ik), mySam%mpar(1), &
               &       sparse%mspar(1), sysMat%meqn8(1), iWork(1), lerr)
          if (bigModel) call writeProgress ('     SPRPRP done',newline=.false.)
          if (lerr < 0_ik) goto 998

          n0 = count(mySam%msc == 0) - n0
          if (n0 > 0) then
             !! The MSC array was modified, need to save the updated version
             if (ik > i4) call castInt (samData%ndof,mySam%msc,samData%msc)
             write(lpu,"('Note :',I6,' unconnected internal DOFs detected')") n0
          end if

       end if

       sysMat%dim = int(sparse%mspar(8))
       if (sysMat%dim < 7) then
          !! Very few unknowns in this system ==> most DOFs are prescribed.
          !! Switch to dense matrix format and release the unneccesary arrays.
          sysMat%storageType = denseMatrix_p
          call reAllocate (prgnam_p,sysMat%sparse)
          if (ik > i4) then
             call reAllocate (prgnam_p,sysMat%meqn8)
          end if
          if (.not.associated(samData%meqn)) then
             samData%meqn => sysMat%meqn
          end if
          if (sysMat%dim > 0) then
             write(tmpDir,"('Only',I2,' free DOFs in this system')") sysMat%dim
          else
             tmpDir = 'No free DOFs in this system'
          end if
          call reportError (note_p,trim(tmpDir)// &
               &            ', switching to Dense matrix format')
          return
       end if

       if (.not. present(cacheSize)) then
          sparse%mspar(9) = 0_ik
       else if (cacheSize > 0) then
          write(lpu,"('Note : Setting sparse-solver cachesize =',I8)") cacheSize
          sparse%mspar(9) = int(cacheSize,ik)
       else
          sparse%mspar(9) = 0_ik
       end if

       call reAllocate (prgnam_p,sparse%msica,int(sparse%mspar(35)),ierr)
       if (ierr < 0) return

       !! Perform symbolic assembly
       if (bigModel) call writeProgress ('     SPRSAS ...',newline=.false.)
       call SPRSAS (mySam%mpar(1),  mySam%mpmnpc(1), mySam%mmnpc(1), &
            &       mySam%madof(1), mySam%msc(1),    mySam%mpmceq(1), &
            &       mySam%mmceq(1), sysMat%meqn8(1), sparse%mspar(1), &
            &       sparse%msica(1), iWork(1), nspar_p, int(lpu,ik), lerr)
       if (bigModel) call writeProgress ('     SPRSAS done',newline=.false.)
       if (lerr < 0) goto 998

       if (ik > i4) call castInt (samData%ndof,sysMat%meqn8,sysMat%meqn)

       !! Reallocation of msica with the correct size
       call reAllocate (prgnam_p,sparse%msica,int(sparse%mspar(2)),ierr,.true.)
       if (ierr < 0) return

    case (outOfCore_p)

       call reAllocate (prgnam_p,sysMat%gsf,.true.,ierr)
       if (ierr < 0) return

       call SetOutputUnit (sysMat%gsf%Msg, lpu)
       call SetLinearSolverType (sysMat%gsf%Msg, GSFFactOpt)
       if (present(cacheSize)) then
          call SetMemorySize (sysMat%gsf%Msg, cacheSize*131072)
       end if
#if defined(win32) || defined(win64)
       call getenv ('TMP',tmpDir)
       if (tmpDir == '') tmpDir = 'C:'
#else
       call getenv ('TMPDIR',tmpDir)
       if (tmpDir == '') tmpDir = '/var/tmp'
#endif
       call SetScratchFileLocation (sysMat%gsf%Msg, tmpDir)
       if (firstMatrix) call WriteMessageData (sysMat%gsf%Msg)

       if (lumpedElMats) then

          !! All element matrix contributions are lumped (diagonalized), and
          !! the only off-diagonal terms may result from constraint equations.
          !! Create a separate SAM structure accounting for these only.

          call LumpedSAM (samData,sysMat%gsf%P%sam,ierr)
          if (ierr /= 0) goto 999

       else
          sysMat%gsf%P%sam => samData
       end if

       !! Initialize the out-of-core matrix management system (this needs to be
       !! done only once, regardless of whether we have one or two matrices)
       if (firstMatrix) then
          call FEInitialize (sysMat%gsf%Msg, GSFinfo)
          ierr = checkGSFinfo(GSFinfo)
          if (ierr /= 0) goto 999
       end if

       !! Perform symbolic assembly, matrix analysis, etc.
       if (.not. present(masterGSF)) then
          call FEAnalyze (sysMat%gsf%Msg, sysMat%gsf%P, &
               &          sysMat%gsf%Q,   sysMat%gsf%T, GSFinfo)
       else if (needsFactorization) then
          !! Let current matrix have the same ordering as given by masterGSF
          call FEAnalyze (sysMat%gsf%Msg, sysMat%gsf%P, masterGSF%Q, &
               &          sysMat%gsf%Q,   sysMat%gsf%T, GSFinfo)
       else
          !! Let current matrix have the same ordering as given by masterGSF and
          !! ommit the matrix analysis since this matrix will not be factorized
          call FEAnalyze (sysMat%gsf%Msg, sysMat%gsf%P, masterGSF%Q, &
               &          sysMat%gsf%Q,   GSFinfo)
       end if
       ierr = checkGSFinfo(GSFinfo)
       if (ierr /= 0) goto 999

       if (firstMatrix) then

          !! Allocate the matrix of equation numbers
          call reAllocate (prgnam_p,samData%meqn,samData%ndof,ierr)
          if (ierr < 0) return

          !! Get the total number of equations and those of status 1 and 2
          samData%neq = FEMatrixOrder(sysMat%gsf%Q)
          call FEGetSuperelementOrder (sysMat%gsf%Q,samData%ndof1,samData%ndof2)

          !! Get the matrix of equation numbers for the active DOFs
          !! The numbers will here be in the external, application-defined,
          !! order and thus the same for both matrices
          call FEGetActiveOrder (sysMat%gsf%Q,samData%meqn,GSFinfo)
          ierr = checkGSFinfo(GSFinfo)
          if (ierr /= 0) goto 999

          if (matrixType == outOfCore_p + sparseMatrix_p) then
             call reAllocate (prgnam_p,sysMat%meqn,samData%ndof,ierr)
             if (ierr < 0) return

             !! Make a copy of the matrix of equation numbers for this matrix
             !! and permute them to the internal, GSF-defined, order. This will
             !! only be used to establish the internal permutation from the
             !! SPR-order to GSF-order that is needed during eigenvalue analysis

             call ICOPY (samData%ndof,samData%meqn(1),1,sysMat%meqn(1),1)
             call FEPermuteIndices ('T',sysMat%gsf%Q,sysMat%meqn,GSFinfo)
             ierr = checkGSFinfo(GSFinfo)
             if (ierr /= 0) goto 999

          else
             !! Matrix of equation numbers for the first matrix
             sysMat%meqn => samData%meqn
          end if

       else
          !! Matrix of equation numbers for the second matrix
          sysMat%meqn => samData%meqn
       end if

       !! Store the dimension of this matrix
       sysMat%dim = samData%neq

    case (pardiso_p)

       call reAllocate (prgnam_p,sysMat%meqn,samData%ndof,ierr)
       call reAllocate (prgnam_p,sysMat%pardiso,.true.,ierr)
       if (ierr < 0) return

       if (matrixType > 0) then
          sysMat%pardiso%mtype =  2 ! Real, symmetric, positive definite matrix
       else
          sysMat%pardiso%mtype = -2 ! Real, symmetric, indefinite matrix
       end if

    end select

    if (.not.associated(samData%meqn) .and. associated(sysMat%meqn)) then
       samData%meqn => sysMat%meqn
    end if
    if (associated(samData%meqn) .and. present(ipsw)) then
       if (ipsw > 1 .and. sysMat%storageType == sparseMatrix_p) then
          call writeMEQN ('csAllocPointerArrays',samData%madof,samData%meqn,lpu)
       end if
    end if

    return

998 ierr = int(lerr)
999 call reportError (error_p, &
         'Failed to initialize system matrix pointer arrays',addString=prgnam_p)

  end subroutine csAllocPointerArrays


  !!============================================================================
  !> @brief Performs nodal renumbering for the equation solver.
  !>
  !> @param samData Data for managing system matrix assembly
  !> @param sysMat The system matrix to perform renumbering for
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !> @date Sep 1998
  !>
  !> @author Knut Morten Okstad
  !> @date Jan 2003

  subroutine csRenumber (samData,sysMat,lpu,ierr)

    use kindModule         , only : i4
    use SamModule          , only : SamType
    use SysMatrixTypeModule, only : SysMatrixType
    use SysMatrixTypeModule, only : skylineMatrix_p, sparseMatrix_p
    use allocationModule   , only : reAllocate
    use progressModule     , only : writeProgress
    use scratchArrayModule , only : getIntegerScratchArray
    use reportErrorModule  , only : error_p, reportError

    type(SamType)      , intent(inout) :: samData
    type(SysMatrixType), intent(inout) :: sysMat
    integer            , intent(in)    :: lpu
    integer            , intent(out)   :: ierr

    !! Local variables
    integer, parameter :: ipsw = 1, iflag = 1
    integer(ik)        :: lerr
    integer            :: nWork, nstart, neff
    integer, pointer   :: iWork(:)

    character(len=30), parameter :: prgnam_p = 'AsmExtensionModule::csRenumber'

    !! --- Logic section ---

    ierr = 0

    select case (sysMat%storageType)

    case (skylineMatrix_p)

       nWork =  swsWorkSize(samdata%mpmnpc)
       iWork => getIntegerScratchArray(nWork,ierr)
       if (ierr < 0) goto 999

       call reAllocate (prgnam_p,samData%mnnn,samData%nnod,ierr)
       if (ierr < 0) return

       nstart = 0 ! Automatically determine the starting node
       call SWSNUM (samData%mpmnpc(1), samData%mmnpc(1), samData%minex(1), &
            &       samData%msc(1),    samData%madof(1), iWork(1), &
            &       samData%mnnn(1),   samData%nnod,     samData%nel, &
            &       nWork, ipsw, lpu, iflag, nstart, neff, ierr)
       if (ierr < 0) goto 999

    case (sparseMatrix_p)

       nWork = int(sysMat%sparse%mspar(37))
       if (ik > i4) nWork = nWork*2 ! We are using 64-bit integers
       iWork => getIntegerScratchArray(nWork,ierr)
       if (ierr < 0) goto 999

       if (bigModel) call writeProgress ('     SPRRNM ...',newline=.false.)
       call SPRRNM (sysMat%sparse%mspar(1), sysMat%sparse%msica(1), &
            &       iWork(1), int(lpu,ik), lerr)
       if (bigModel) call writeProgress ('     SPRRNM done',newline=.false.)
       if (lerr < 0) goto 998

    end select

    return

998 ierr = int(lerr)
999 call reportError (error_p,'Failed to optimize equation system ordering', &
         &            addString=prgnam_p)

  contains

    !> @brief Returns the work space size required by SAM subroutine SWSNUM.
    function swsWorkSize (mpmnpc)
      integer, intent(in) :: mpmnpc(:)
      integer             :: i, nn, nw1, nw2, swsWorkSize
      nw1 = samData%nnod
      do i = 1, size(mpmnpc)-1
         nn = mpmnpc(i+1) - mpmnpc(i)
         nw1 = nw1 + nn*(nn-1)
      end do
      nw2 = 4*samData%nnod + 2 + 2*samData%nel
      swsWorkSize = max(nw1,nw2)*2
    end function swsWorkSize

  end subroutine csRenumber


  !!============================================================================
  !> @brief Performs pre-assembly for the given system matrix.
  !>
  !> @param samData Data for managing system matrix assembly
  !> @param sysMat The system matrix to perform per-assembly for
  !> @param[in] doSupelReduction If .true., we are doing superelement reduction
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !> @param[in] ipsw Print switch
  !> @param[in] meqnK Matrix of equation numbers for system stiffness matrix
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !> @date Sep 1998
  !>
  !> @author Knut Morten Okstad
  !> @date Jan 2003

  subroutine csPreassemble (samData,sysMat,doSupelReduction,lpu,ierr,ipsw,meqnK)

    use sprKindModule      , only : i4, castInt
    use SamModule          , only : SamType
    use SysMatrixTypeModule, only : SysMatrixType, nspar_p
    use SysMatrixTypeModule, only : SkylineStorageType, SparseStorageType
    use SysMatrixTypeModule, only : matrixType_p, diagonalMatrix_p
    use SysMatrixTypeModule, only : skylineMatrix_p, sparseMatrix_p
    use SysMatrixTypeModule, only : denseMatrix_p, outOfCore_p, pardiso_p
#ifdef FT_HAS_MKL
    use MKL_pardiso        , only : pardiso
#ifndef win32
    use MKL_pardiso        , only : pardisoinit
#endif
#endif
    use allocationModule   , only : reAllocate
    use scratchArrayModule , only : getIntegerScratchArray, releaseScratchArrays
    use reportErrorModule  , only : error_p, reportError

    type(SamType)      , intent(inout) :: samData
    type(SysMatrixType), intent(inout) :: sysMat
    logical            , intent(in)    :: doSupelReduction
    integer            , intent(in)    :: lpu
    integer            , intent(out)   :: ierr
    integer, optional  , intent(in)    :: ipsw, meqnK(:)

    !! Local variables
    integer(ik)          :: lerr
    integer              :: nWork
    integer    , pointer :: iWork(:)
    integer(ik), pointer :: meqnK8(:)

    type(SkylineStorageType), pointer :: skyline
    type(SparseStorageType) , pointer :: sparse

    character(len=33), parameter :: prgnam_p = &
         'AsmExtensionModule::csPreAssemble'

    !! --- Logic section ---

    select case (sysMat%storageType)

    case (skylineMatrix_p)

       skyline => sysMat%skyline
       nWork   =  prassWorkSize(samData%mpmnpc,samData%mmnpc,samData%madof)
       iWork   => getIntegerScratchArray(nWork,ierr)
       if (ierr < 0) goto 999

       call PRASSX (samData%madof(1), samData%minex(1), samData%mpmnpc(1), &
            &       samData%mmnpc(1), samData%msc(1),   samData%mpmceq(1), &
            &       samData%mmceq(1), samData%mnnn(1),  lpu, &
            &       samData%mpar(1),  iWork(1),         sysMat%meqn(1), &
            &       skyline%msky(1),  ierr)
       if (ierr < 0) go to 999

       skyline%eigenSolver = 2
       skyline%neq = samData%neq
       sysMat%dim  = samData%neq

       !! Reallocation of msky with the correct size
       call reAllocate (prgnam_p,skyline%msky,skyline%neq,ierr,.true.)
       if (ierr < 0) return

    case (sparseMatrix_p)

       sparse => sysMat%sparse
       nWork = int(sparse%mspar(38))
       if (ik > i4) nWork = nWork*2 ! We are using 64-bit integers
       iWork => getIntegerScratchArray(nWork,ierr)
       if (ierr < 0) goto 999

       call reAllocate (prgnam_p,sparse%mtrees,int(sparse%mspar(36)),ierr)
       if (ierr < 0) return

       call SPRTRS (sparse%mspar(1), sparse%msica(1), sparse%mtrees(1), &
            &       sysMat%meqn8(1), nspar_p,         sparse%mspar(2), &
            &       sparse%mspar(36), mySam%mpar(3),  sparse%mspar(38), &
            &       iWork(1),        sparse%rinfo(1), int(lpu,ik), lerr)
       if (lerr < 0_ik) goto 998

       call reAllocate (prgnam_p,sparse%msifa,int(sparse%mspar(3)),ierr)
       if (ierr < 0) return

       if (sparse%mspar(51) > 0_ik .and. present(meqnK)) then
#ifdef FT_SPR_INT_DEFAULT
          call SPRABC (sparse%mspar(1), sparse%msifa(1), sysMat%meqn(1), &
               &       meqnK(1),        mySam%mpar(3),   int(lpu,ik), lerr)
#else
          nullify(meqnK8)
          call reAllocate (prgnam_p,meqnK8,size(meqnK),ierr)
          if (ierr < 0) return
          call castInt (size(meqnK),meqnK,meqnK8)
          call SPRABC (sparse%mspar(1), sparse%msifa(1), sysMat%meqn8(1), &
               &       meqnK8(1),       mySam%mpar(3),   int(lpu,ik), lerr)
          call reAllocate (prgnam_p,meqnK8)
#endif
          if (lerr < 0_ik) goto 998
       end if

       nWork = 3*int(sparse%mspar(6)) + 2*int(sparse%mspar(8)) + 1
       if (ik > i4) nWork = nWork*2 ! We are using 64-bit integers
       iWork => getIntegerScratchArray(nWork,ierr)
       if (ierr < 0) goto 999

       call SPRSMB (sparse%mspar(1), sparse%msica(1), sparse%mtrees(1), &
            &       sparse%msifa(1), sysMat%meqn8(1), iWork(1), &
            &       int(lpu,ik), lerr)
       if (lerr < 0) goto 998

       !! Reallocation of msifa with the correct size
       call reAllocate (prgnam_p,sparse%msifa,int(sparse%mspar(3)),ierr,.true.)
       call reAllocate (prgnam_p,sparse%mvarnc,2*int(sparse%mspar(8)),ierr)
       if (ierr < 0) return

       call SPRPMP (sparse%msica(1), sparse%mtrees(1), sysMat%meqn8(1), &
            &       mySam%mpar(1), sparse%mspar(1),  sparse%mvarnc(1))

       if (ik > i4) then
          call castInt (size(samData%mpar),mySam%mpar,samData%mpar)
          call castInt (samData%ndof,sysMat%meqn8,sysMat%meqn)
       end if

    case (diagonalMatrix_p, denseMatrix_p)

       call SYSEQ (samData%msc(1),  samData%mpmceq(1), samData%mmceq(1), lpu, &
            &      samData%mpar(1), sysMat%meqn(1),    ierr)
       if (ierr < 0) goto 999

       sysMat%dim = samData%neq

    case (pardiso_p)

#ifdef FT_HAS_MKL
       call SYSEQ (samData%msc(1),  samData%mpmceq(1), samData%mmceq(1), lpu, &
            &      samData%mpar(1), sysMat%meqn(1),    ierr)
       if (ierr < 0) goto 999

       sysMat%dim = samData%neq

       !! Compute the sparsity pattern of the system matrix
       call csPreAssSparse (samData,sysMat%pardiso,ierr)
       if (ierr < 0) goto 999

#ifndef win32
       !! Initialize default control parameters for MKL PARDISO
       call pardisoinit (sysMat%pardiso%pt, sysMat%pardiso%mtype, &
            &            sysMat%pardiso%iparm)
#endif

       sysMat%pardiso%iparm(6)  = 1 ! Store solution in the RHS vector
       sysMat%pardiso%iparm(27) = 1 ! Invoke the matrix checker

       !! Perform symbolic factorization (phase=11)
       call pardiso (sysMat%pardiso%pt, 1, 1, &
            &        sysMat%pardiso%mtype, 11, samData%neq, &
            &        sysMat%value, sysMat%pardiso%ia, sysMat%pardiso%ja, &
            &        sysMat%ipiv, 1, sysMat%pardiso%iparm, 0, &
            &        sysMat%value, sysMat%value, ierr)
       if (ierr < 0) then
          write(lpu,"('Error :   PARDISO returned IERR = ',I6)") ierr
          goto 999
       end if
#else
       ierr = -999
       write(lpu,"('Error :   MKL PARDISO is not available')")
       goto 999
#endif

    end select

    if (present(ipsw)) then
       if (ipsw > 2) then
          call writeMEQN ('csPreassemble',samData%madof,samData%meqn,lpu)
       end if
    end if

    if (doSupelReduction) then

       !! Establish equation system partitioning into external and internal DOFs

       iWork => getIntegerScratchArray(samData%neq,ierr)
       if (ierr < 0) goto 999

       call csCreateEqPartition (sysMat%storageType /= outOfCore_p, &
            &                    samData,iWork,ierr)
       if (ierr < 0) return

    end if

    call releaseScratchArrays
    return

998 ierr = int(lerr)
999 call reportError (error_p, &
         'Failed to pre-assemble system '//matrixType_p(sysMat%storageType), &
         addString=prgnam_p)

  contains

    !> @brief Returns the work space size required by SAM subroutine PRASSX.
    function prassWorkSize (mpmnpc,mmnpc,madof) result(maxElDofs)
      integer, intent(in) :: mpmnpc(:), mmnpc(:), madof(:)
      integer             :: iel, elDofs, ip, gNod, maxElDofs
      maxElDofs = samData%nnod
      do iEl = 1, size(mpmnpc)-1
         elDofs = 0
         do ip = mpmnpc(iel), mpmnpc(iel+1)-1
            gNod = mmnpc(ip)
            elDofs = elDofs + madof(gnod+1) - madof(gnod)
         end do
         if (elDofs > maxElDofs) maxElDofs = elDofs
      end do
    end function prassWorkSize

  end subroutine csPreassemble


  !!============================================================================
  !> @brief Partitions the equation system into external and internal DOFs.
  !>
  !> @param[in] is2orderedLast If .true., the status-2 DOFs are orered last
  !> @param samData Data for managing system matrix assembly
  !> @param[out] eqnStatus Status flags for the equations
  !> @param[out] ierr Error flag
  !>
  !> @callergraph
  !>
  !> @brief Knut Morten Okstad
  !>
  !> @date 16 Jan 2003

  subroutine csCreateEqPartition (is2orderedLast,samData,eqnStatus,ierr)

    use SamModule        , only : SamType, invertEqPartition
    use allocationModule , only : reAllocate
    use reportErrorModule, only : note_p, debugFileOnly_p
    use reportErrorModule, only : reportError, internalError

    logical      , intent(in)    :: is2OrderedLast
    type(SamType), intent(inout) :: samData
    integer      , intent(out)   :: eqnStatus(:), ierr

    !! Local variables
    integer :: ii, ie, ieq, idof, status

    character(len=39), parameter :: prgnam_p = &
         'AsmExtensionModule::csCreateEqPartition'

    !! --- Logic section ---

    ierr = 0
    call reAllocate (prgnam_p,samData%meqn1,samData%ndof1,ierr)
    call reAllocate (prgnam_p,samData%meqn2,samData%ndof2,ierr)
    if (ierr < 0) return

    !! Find the status code for all equation numbers

    ie = 0
    eqnStatus = 0
    do idof = 1, samData%ndof
       status = samData%msc(idof)
       if (status == 1 .or. status == 2) then
          ieq = samData%meqn(idof)
          eqnStatus(ieq) = status
          if (status == 2 .and. ieq < samData%ndof1 .and. is2OrderedLast) then
             ie = ie + 1
          end if
       end if
    end do
    if (ie > 0) then
       call reportError (note_p,'This part has more than one internal/'// &
            'external DOF block','in the system stiffness- and mass matrices.')
    end if

    !! Establish list of equation numbers with status 1 and 2
    ii = 0
    ie = 0
    do ieq = 1, samData%neq
       if (eqnStatus(ieq) == 1) then
          ii = ii + 1
          if (ii <= samData%ndof1) samData%meqn1(ii) = ieq
       else if (eqnStatus(ieq) == 2) then
          ie = ie + 1
          if (ie <= samData%ndof2) samData%meqn2(ie) = ieq
       end if
    end do

    if (ii /= samData%ndof1) then
       ierr = internalError(prgnam_p//': Sum ii /= ndof1')
    else if (ie /= samData%ndof2) then
       ierr = internalError(prgnam_p//': Sum ie /= ndof2')
    else if (samData%ndof2 > 0) then
       !! Form the control array dofPosIn2 (status 2 eq. numbers in DOF-order)
       call invertEqPartition (.false.,.true.,samData,ierr)
       if (ierr < 0) call reportError (debugFileOnly_p,prgnam_p)
    end if

  end subroutine csCreateEqPartition


  !!============================================================================
  !> @brief Initializes the element assembly.
  !>
  !> @param sysMat The system matrix to initialize element assembly for
  !> @param[in] matName Matrix name, used for storage requirement print
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @callergraph
  !>
  !> @brief Knut Morten Okstad
  !>
  !> @date 8 Sep 2004

  subroutine csBeginAssembly (sysMat,matName,lpu,ierr)

    use SysMatrixTypeModule, only : SysMatrixType, outOfCore_p, dp
    use SysMatrixTypeModule, only : allocateSysMatrix, checkGSFinfo, GSFinfo
    use FELinearSolver     , only : FEBeginAssembly
    use reportErrorModule  , only : debugFileOnly_p, reportError

    type(SysMatrixType), intent(inout) :: sysMat
    character(len=*)   , intent(in)    :: matName
    integer            , intent(in)    :: lpu
    integer            , intent(out)   :: ierr

    !! --- Logic section ---

    ierr = 0
    call allocateSysMatrix (sysMat,ierr,lpu,matName)
    if (ierr /= 0) goto 999

    if (sysMat%storageType == outOfCore_p) then
       call FEBeginAssembly (sysMat%gsf%Q, sysMat%gsf%S, GSFinfo)
       ierr = checkGSFinfo(GSFinfo)
       if (ierr /= 0) goto 999
    else
       sysMat%value = 0.0_dp
    end if

    return

999 continue
    call reportError (debugFileOnly_p,'AsmExtensionModule::csBeginAssembly')

  end subroutine csBeginAssembly


  !!============================================================================
  !> @brief Finalizes the element assembly.
  !>
  !> @param sysMat The system matrix to finalize element assembly for
  !> @param[out] ierr Error flag
  !>
  !> @callergraph
  !>
  !> @brief Knut Morten Okstad
  !>
  !> @date 8 Sep 2004

  subroutine csEndAssembly (sysMat,ierr)

    use SysMatrixTypeModule, only : SysMatrixType, outOfCore_p
    use SysMatrixTypeModule, only : allocateSysMatrix, checkGSFinfo, GSFinfo
    use FELinearSolver     , only : FEEndAssembly
    use reportErrorModule  , only : debugFileOnly_p, reportError

    type(SysMatrixType), intent(inout) :: sysMat
    integer            , intent(out)   :: ierr

    !! --- Logic section ---

    ierr = 0
    if (sysMat%storageType == outOfCore_p) then
       sysMat%gsf%isFactorized = .false.
       if (associated(sysMat%gsf%T)) then
          call FEEndAssembly (sysMat%gsf%Msg, sysMat%gsf%Q, &
               &              sysMat%gsf%S, sysMat%gsf%T, GSFinfo)
       else
          call FEEndAssembly (sysMat%gsf%Q, sysMat%gsf%S)
       end if
       ierr = checkGSFinfo(GSFinfo)
       if (ierr /= 0) goto 999
    end if

    return

999 continue
    call reportError (debugFileOnly_p,'AsmExtensionModule::csEndAssembly')

  end subroutine csEndAssembly


  !!============================================================================
  !> @brief Expands an element matrix to comply with sammodule::samtype::madof.
  !>
  !> @param[in] elMat Element stiffness matrix as returned from element routine
  !> @param[out] expMat Expanded elemetn matrix
  !> @param[in] mpmnpc Matrix of pointers to MNPCs
  !> @param[in] mmnpc Matrix of nodal point correspondances
  !> @param[in] madof Matrix of accumulated DOFs
  !> @param[in] iel Element index
  !> @param[in] nedof Number of rows and columns in @a elMat
  !> @param[out] ierr Error flag
  !>
  !> @details The element matrix @a elMat is expanded such that it matches
  !> the number of DOFs per node as defined by the sam array MADOF.
  !>
  !> @callergraph
  !>
  !> @brief Knut Morten Okstad
  !>
  !> @date 2 Jun 2005

  subroutine csExpandElMat (elMat, expMat, mpmnpc, mmnpc, madof, &
       &                    iel, nedof, ierr)

    use scratchArrayModule, only : getRealScratchMatrix, dp
    use reportErrorModule , only : debugFileOnly_p, reportError, internalError

    real(dp), intent(in), target :: elMat(:,:)
    real(dp), pointer            :: expMat(:,:)
    integer , intent(in)         :: mpmnpc(:), mmnpc(:), madof(:), iel, nedof
    integer , intent(out)        :: ierr

    !! Local variables
    integer :: ien, jen, ip1, ip2, jp1, jp2, inod, jnod, ipmnpc
    integer :: nelnod, neldof, nndof, ndof, nd1

    !! --- Logic section ---

    ipmnpc = mpmnpc(iel)
    nelnod = mpmnpc(iel+1) - ipmnpc
    neldof = 0
    do ien = 1, nelnod
       inod = mmnpc(ipmnpc+ien-1)
       neldof = neldof + madof(inod+1) - madof(inod)
    end do

    if (nedof == neldof) then
       !! No expansion needed, the element dimension matches the size in madof
       expMat => elMat
       ierr = 0
       return
    else if (nedof > neldof) then
       ierr = internalError('csExpandElMat: Element matrix is too big')
       return
    end if

    expMat => getRealScratchMatrix(neldof,neldof,ierr)
    if (ierr < 0) then
       call reportError (debugFileOnly_p,'csExpandElMat')
       return
    end if
    expMat = 0.0_dp

    nndof = nedof/nelnod
    ip1 = 1
    ip2 = 1
    nd1 = nndof-1
    do ien = 1, nelnod
       jp1 = 1
       jp2 = 1
       inod = mmnpc(ipmnpc+ien-1)
       ndof = madof(inod+1) - madof(inod)
       if (ndof < nndof) then
          ierr = internalError('csExpandElMat: Inconsistency in MADOF')
          return
       end if
       do jen = 1, nelnod
          expMat(ip1:ip1+nd1,jp1:jp1+nd1) = elMat(ip2:ip2+nd1,jp2:jp2+nd1)
          jnod = mmnpc(ipmnpc+jen-1)
          jp1 = jp1 + madof(jnod+1) - madof(jnod)
          jp2 = jp2 + nndof
       end do
       ip1 = ip1 + ndof
       ip2 = ip2 + nndof
    end do

  end subroutine csExpandElMat


  !!============================================================================
  !> @brief Adds an element matrix into the system matrix.
  !>
  !> @param[in] samData Data for managing system matrix assembly
  !> @param[in] iel Element index
  !> @param[in] nedof Number of DOFs in element
  !> @param[in] elMat The element matrix to add
  !> @param sysMat The system matrix to add into
  !> @param[out] ierr Error flag
  !> @param sysRhs System right-hand-side vector
  !>
  !> @details This subroutine does the same as asmextensionmodule::csaddem,
  !> but uses assumed-size declaration of @a elMat instead.
  !> It should therefore be used if the calling subroutine declares @a elMat as
  !> a 2D array @a elMat(N,N), where @a N &gt; @a nedof,
  !> and @a elMat(1:nedof,1:nedof) is passed as argument.
  !> Calling asmextensionmodule::csaddem directly will not work in this case.
  !>
  !> @callergraph
  !>
  !> @brief Knut Morten Okstad
  !>
  !> @date 28 Oct 2002

  subroutine csAddElM (samData,iel,nedof,elMat,sysMat,ierr,sysRhs)

    use SamModule          , only : SamType, dp
    use SysMatrixTypeModule, only : SysMatrixType

    type(SamType)      , intent(in)    :: samData
    integer            , intent(in)    :: iel, nedof
    real(dp)           , intent(in)    :: elMat(nedof,nedof)
    type(SysMatrixType), intent(inout) :: sysMat
    integer            , intent(out)   :: ierr
    real(dp), optional , intent(inout) :: sysRhs(:)

    !! --- Logic section ---

    call csAddEM (samData,iel,nedof,elMat,sysMat,ierr,sysRhs)

  end subroutine csAddElM


  !!============================================================================
  !> @brief Adds an element matrix into the system matrix.
  !>
  !> @param[in] samData Data for managing system matrix assembly
  !> @param[in] iel Element index
  !> @param[in] nedof Number of DOFs in element
  !> @param[in] elMat The element matrix to add
  !> @param sysMat The system matrix to add into
  !> @param[out] ierr Error flag
  !> @param sysRhs System right-hand-side vector
  !> @param[in] lpu File unit number for debug print
  !> @param[in] ipsw Print switch
  !>
  !> @callergraph
  !>
  !> @note The argument @a elMat is here declared assumed-shape (KMO 21.10.02).
  !>
  !> @author Bjorn Haugen
  !> @date Sep 1998
  !>
  !> @author Knut Morten Okstad
  !> @date Mar 2003

  subroutine csAddEM (samData,iel,nedof,elMat,sysMat,ierr,sysRhs,lpu,ipsw)

    use kindModule         , only : dp, i4
    use SamModule          , only : SamType
    use SysMatrixTypeModule, only : SysMatrixType, GSFinfo, checkGSFinfo
    use SysMatrixTypeModule, only : SkylineStorageType, SparseStorageType
    use SysMatrixTypeModule, only : matrixType_p
    use SysMatrixTypeModule, only : skylineMatrix_p, sparseMatrix_p
    use SysMatrixTypeModule, only : denseMatrix_p, outOfCore_p, pardiso_p
    use SysMatrixTypeModule, only : asmSparse
    use FELinearSolver     , only : FEAssembleElement
    use scratchArrayModule , only : getIntegerScratchArray
    use TimerModule        , only : startTimer, stopTimer
    use reportErrorModule  , only : internalError, getErrorFile
    use reportErrorModule  , only : reportError, error_p

    type(SamType)      , intent(in)            :: samData
    integer            , intent(in)            :: iel, nedof
    real(dp)           , intent(in)            :: elMat(:,:)
    type(SysMatrixType), intent(inout)         :: sysMat
    integer            , intent(out)           :: ierr
    real(dp), optional , intent(inout), target :: sysRhs(:)
    integer , optional , intent(in)            :: lpu, ipsw

    !! Local variables
    integer(ik)       :: lerr
    integer           :: nrhs
    integer , pointer :: iWork(:)
    real(dp), pointer :: rhs(:), tmpMat(:,:)

    type(SkylineStorageType), pointer :: skyline
    type(SparseStorageType) , pointer :: sparse

    character(len=27), parameter :: prgnam_p = 'AsmExtensionModule::csAddEM'

    !! --- Logic section ---

    if (size(elMat,1) /= nedof .or. size(elMat,2) < nedof) then
       ierr = internalError(prgnam_p//': Invalid element matrix dimension')
       write(getErrorFile(),"(39X,3I6)") size(elMat,1),size(elMat,2),nedof
       return
    else if (present(sysRhs)) then
       nrhs = 1
       rhs => sysRhs ! contributions from prescribed DOFs will be added
    else
       nrhs = 0
       rhs => sysMat%value ! dummy, but should point to something
    end if

    call startTimer (16)

    !! Note: The assembly routines below (except for FEAssembleElement) uses
    !! samData%mpar(18) to determine how the element matrices are organized.
    !! = 0 : Elements with a varying number of nodal DOFs (NNDOF) are allowed
    !!       (e.g., superelements with component modes in the Dynamics Solver).
    !!       For each node I it is assumed that NNDOF = MADOF(I+I) - MADOF(I).
    !! > 0 : All element nodes are assumed to have NNDOF = NEDOF/NENOD dofs.
    !!       However, this may be less than MADOF(I+1) - MADOF(I) for that node.
    !!       This is the case for, e.g., rigid elements in the FE Part Reducer.

    select case (sysMat%storageType)

    case (skylineMatrix_p)

       skyline => sysMat%skyline
       iWork   => getIntegerScratchArray(nedof,ierr)
       if (ierr < 0) goto 999

       call ADDEM1 (elMat(1,1),       samData%ttcc(1),   samData%mpar(1),   &
            &       samData%madof(1), sysMat%meqn(1),    samData%mpmnpc(1), &
            &       samData%mmnpc(1), samData%mpmceq(1), samData%mmceq(1),  &
            &       skyline%msky(1),  iel,               nedof,             &
            &       getErrorFile(),   nrhs,              sysMat%value(1),   &
            &       rhs(1),           iWork(1),          ierr)
       if (ierr < 0) goto 999

    case (denseMatrix_p)

       iWork => getIntegerScratchArray(nedof,ierr)
       if (ierr < 0) goto 999

       call ADDEM2 (elMat(1,1),       samData%ttcc(1),   samData%mpar(1),   &
            &       samData%madof(1), sysMat%meqn(1),    samData%mpmnpc(1), &
            &       samData%mmnpc(1), samData%mpmceq(1), samData%mmceq(1),  &
            &       iel,              nedof,             samData%neq,       &
            &       getErrorFile(),   nrhs,              sysMat%value(1),   &
            &       rhs(1),           iWork(1),          ierr)
       if (ierr < 0) goto 999

    case (sparseMatrix_p)

       sparse => sysmat%sparse
       if (ik > i4) then ! We are using 64-bit integer
          iWork => getIntegerScratchArray(int(sparse%mspar(39))*2,ierr)
       else
          iWork => getIntegerScratchArray(int(sparse%mspar(39)),ierr)
       end if
       if (ierr < 0) goto 999

       call SPRADM (elMat(1,1),      samData%ttcc(1),  mySam%mpar(1),    &
            &       sparse%mspar(1), mySam%madof(1),   sysMat%meqn8(1),  &
            &       mySam%mpmnpc(1), mySam%mmnpc(1),   mySam%mpmceq(1),  &
            &       mySam%mmceq(1),  sparse%msica(1),  sparse%mtrees(1), &
            &       sparse%msifa(1), sparse%mvarnc(1), sysMat%value(1),  &
            &       rhs(1),          iWork(1),         int(iel,ik),      &
            &       int(nedof,ik), int(getErrorFile(),ik), int(nrhs,ik), lerr)
       if (lerr < 0_ik) goto 998

    case (outOfCore_p)

       !! Expand the element matrix, if necessary, to comply with madof
       call csExpandElMat (elMat, tmpMat, &
            &              samData%mpmnpc, samData%mmnpc, samData%madof, &
            &              iel, nedof, ierr)
       if (ierr < 0) goto 999

       if (present(IPSW) .and. present(LPU)) then
          if (IPSW == -7) then ! Export to ACD test program
             write(LPU,"(/'EK',I8,':',I4)") iel, size(tmpMat,1)
             write(LPU,"(1P6E22.15)") tmpMat
          end if
       end if

       call FEAssembleElement (1, iel, tmpMat, sysMat%gsf%P, sysMat%gsf%Q, &
            &                  sysMat%gsf%S, nrhs, rhs, GSFinfo)
       ierr = checkGSFinfo(GSFinfo)
       if (ierr < 0) goto 999

    case (pardiso_p)

       iWork => getIntegerScratchArray(nedof,ierr)
       if (ierr < 0) goto 999

       call csGetElmEq (samData%madof, samData%meqn, samData%mpmnpc, &
            &           samData%mmnpc, samData%mpar(18), iel, ierr, iWork)
       ierr = -abs(ierr-nedof)
       if (ierr < 0) goto 999

       call asmSparse (elMat, iWork, samData%meqn, &
            &          samData%mpmceq, samData%mmceq, samData%ttcc, &
            &          sysMat%value, sysMat%pardiso%IA, sysMat%pardiso%JA, &
            &          sysRhs, ierr)
       if (ierr < 0) goto 999

    case default

       ierr = internalError(prgnam_p//': Invalid matrix type')

    end select

    call stopTimer (16)
    return

998 ierr = int(lerr)
999 call stopTimer (16)
    call reportError (error_p, &
         &            'Failed to add element matrix into system '// &
         &            matrixType_p(sysMat%storageType),addString=prgnam_p)

  end subroutine csAddEM


  !!============================================================================
  !> @brief Adds a diagonal element matrix into the system matrix.
  !>
  !> @param[in] samData Data for managing system matrix assembly
  !> @param[in] iel Element index
  !> @param[in] nedof Number of DOFs in element
  !> @param[in] elMat The diagonal element matrix to add
  !> @param sysMat The system matrix to add into
  !> @param[out] ierr Error flag
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date Mar 2003

  subroutine csAddED (samData,iel,nedof,elMat,sysMat,ierr)

    use kindModule         , only : dp, i4
    use SamModule          , only : SamType
    use SysMatrixTypeModule, only : SysMatrixType, GSFinfo, checkGSFinfo
    use SysMatrixTypeModule, only : SparseStorageType
    use SysMatrixTypeModule, only : matrixType_p, diagonalMatrix_p
    use SysMatrixTypeModule, only : sparseMatrix_p, outOfCore_p, pardiso_p
    use SysMatrixTypeModule, only : asmSparse
    use FELinearSolver     , only : FEAssembleElement
    use scratchArrayModule , only : getIntegerScratchArray
    use TimerModule        , only : startTimer, stopTimer
    use reportErrorModule  , only : internalError, getErrorFile
    use reportErrorModule  , only : reportError, error_p

    type(SamType)      , intent(in)    :: samData
    integer            , intent(in)    :: iel, nedof
    real(dp)           , intent(in)    :: elMat(:)
    type(SysMatrixType), intent(inout) :: sysMat
    integer            , intent(out)   :: ierr

    !! Local variables
    integer(ik)        :: lerr
    integer            :: i, j, ipnod, nndof
    integer, parameter :: nrhs_p = 0
    integer, pointer   :: iWork(:)

    type(SparseStorageType), pointer :: sparse

    character(len=27), parameter :: prgnam_p = 'AsmExtensionModule::csAddED'

    !! --- Logic section ---

    if (size(elMat) /= nedof) then
       ierr = internalError(prgnam_p//': Invalid element matrix dimension')
       write(getErrorFile(),"(39X,2I6)") size(elMat),nedof
       return
    end if

    call startTimer (16)

    select case (sysMat%storageType)

    case (diagonalMatrix_p)

       iWork => getIntegerScratchArray(nedof,ierr)
       if (ierr < 0) goto 999

       call ADDED (elMat(1),         samData%ttcc(1),   samData%mpar(1),   &
            &      samData%madof(1), sysMat%meqn(1),    samData%mpmnpc(1), &
            &      samData%mmnpc(1), samData%mpmceq(1), samData%mmceq(1),  &
            &      iel,              nedof,             getErrorFile(),    &
            &      nrhs_p,           sysMat%value(1),   sysMat%value(1),   &
            &      iWork(1),         ierr)
       if (ierr < 0) goto 999

    case (sparseMatrix_p)

       sparse => sysmat%sparse
       if (ik > i4) then ! We are using 64-bit integer
          iWork => getIntegerScratchArray(int(sparse%mspar(39))*2,ierr)
       else
          iWork => getIntegerScratchArray(int(sparse%mspar(39)),ierr)
       end if
       if (ierr < 0) goto 999

       if (sparse%mspar(51) /= 1_ik) then
          ierr = internalError(prgnam_p//': Invalid matrix type')
          return
       end if

       call SPRADM (elMat(1),        samData%ttcc(1),  mySam%mpar(1),  &
            &       sparse%mspar(1), mySam%madof(1),   sysMat%meqn8(1),  &
            &       mySam%mpmnpc(1), mySam%mmnpc(1),   mySam%mpmceq(1),  &
            &       mySam%mmceq(1),  sparse%msica(1),  sparse%mtrees(1), &
            &       sparse%msifa(1), sparse%mvarnc(1), sysMat%value(1),  &
            &       sysMat%value(1), iWork(1),         int(iel,ik),      &
            &       int(nedof,ik),   int(getErrorFile(),ik), 0_ik, lerr)
       if (lerr < 0_ik) goto 998

    case (outOfCore_p)

       !TODO,kmo: This works only when NEDOF = NNDOF*NENOD and NNDOF = constant
       ipnod = samData%mpmnpc(iel)
       nndof = nedof/(samData%mpmnpc(iel+1)-ipnod)
       do i = 1, nedof
          j = samData%madof(samData%mmnpc(ipnod+(i-1)/NNDOF)) + mod(i-1,NNDOF)
          call FEAssembleElement (1, j, reshape(elMat(i:i),(/1,1/)), &
               &                  sysMat%gsf%P, sysMat%gsf%Q, &
               &                  sysMat%gsf%S, nrhs_p, sysMat%value, GSFinfo)
          ierr = checkGSFinfo(GSFinfo)
          if (ierr < 0) goto 999
       end do

    case (pardiso_p)

       iWork => getIntegerScratchArray(nedof,ierr)
       if (ierr < 0) goto 999

       call csGetElmEq (samData%madof, samData%meqn, samData%mpmnpc, &
            &           samData%mmnpc, samData%mpar(18), iel, nndof, iWork)
       ierr = -abs(nndof-nedof)
       if (ierr < 0) goto 999

       call asmSparse (reshape(elMat,(/nedof,1/)), iWork, samData%meqn, &
            &          samData%mpmceq, samData%mmceq, samData%ttcc, &
            &          sysMat%value, sysMat%pardiso%IA, sysMat%pardiso%JA, &
            &          ierr=ierr)
       if (ierr < 0) goto 999

    case default

       ierr = internalError(prgnam_p//': Invalid matrix type')

    end select

    call stopTimer (16)
    return

998 ierr = int(lerr)
999 call stopTimer (16)
    call reportError (error_p, &
         &            'Failed to add element matrix into system '// &
         &            matrixType_p(sysMat%storageType),addString=prgnam_p)

  end subroutine csAddED


  !!============================================================================
  !> @brief Adds an element vector into a system vector.
  !>
  !> @param[in] samData Data for managing system matrix assembly
  !> @param[in] iel Element index
  !> @param[in] nedof Number of DOFs in element
  !> @param[in] eV The element vector to add
  !> @param sysV The system vector to add into
  !> @param sysReac Vector of reaction forces
  !> @param[out] ierr Error flag
  !>
  !> @details This subroutine takes into account free and dependent DOFs.
  !> Contributions from fixed or prescribed DOFs are added into @a sysReac.
  !>
  !> @callergraph
  !>
  !> @author Karl Erik Thoresen
  !> @date Dec 1998
  !>
  !> @author Knut Morten Okstad
  !> @date Oct 2005

  subroutine csAddEV (samData,iel,nedof,eV,sysV,sysReac,ierr)

    use SamModule         , only : SamType, dp
    use scratchArrayModule, only : getIntegerScratchArray
    use reportErrorModule , only : getErrorFile, reportError, error_p

    type(SamType), intent(in)    :: samData
    integer      , intent(in)    :: iel, nedof
    real(dp)     , intent(in)    :: eV(:)
    real(dp)     , intent(inout) :: sysV(:), sysReac(:)
    integer      , intent(out)   :: ierr

    !! Local variables
    integer, pointer :: meen(:)

    !! --- Logic section ---

    if (associated(samData%mpreac)) then

       call ADDEV1 (eV(1)            , samData%ttcc(1)  , samData%mpar(1)  , &
            &       samData%madof(1) , samData%meqn(1)  , samData%mpmnpc(1), &
            &       samData%mmnpc(1) , samData%mpmceq(1), samData%mmceq(1) , &
            &       samData%mpreac(1), abs(iel), nedof  , sign(1,iel)      , &
            &       getErrorFile()   , sysV(1)          , sysReac(1), ierr)

    else

       meen => getIntegerScratchArray(nedof,ierr)
       if (ierr < 0) goto 999

       call ADDEV2 (eV(1)            , samData%ttcc(1)  , samData%mpar(1)  , &
            &       samData%madof(1) , samData%meqn(1)  , samData%mpmnpc(1), &
            &       samData%mmnpc(1) , samData%mpmceq(1), samData%mmceq(1) , &
            &       abs(iel)         , nedof            , getErrorFile()   , &
            &       sysV(1)          , meen(1)          , ierr)

    end if
    if (ierr < 0) goto 999

    return

999 continue
    call reportError (error_p, &
         &            'Failed to add element vector into system vector', &
         &            addString='AsmExtensionModule::csAddEV')

  end subroutine csAddEV


  !!============================================================================
  !> @brief Adds a nodal vector into a system vector.
  !>
  !> @param[in] samData Data for managing system matrix assembly
  !> @param[in] inod Node index
  !> @param[in] nodV The nodal vector to add
  !> @param sysV The system vector to add into
  !> @param sysReac Vector of reaction forces
  !> @param[out] ierr Error flag
  !>
  !> @details This subroutine takes into account free and dependent DOFs.
  !> Contributions from fixed or prescribed DOFs are added into @a sysReac.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date Oct 2005

  subroutine csAddNV (samData,inod,nodV,sysV,sysReac,ierr)

    use SamModule        , only : SamType, dp
    use reportErrorModule, only : getErrorFile, reportError, error_p

    type(SamType), intent(in)    :: samData
    integer      , intent(in)    :: inod
    real(dp)     , intent(in)    :: nodV(:)
    real(dp)     , intent(inout) :: sysV(:), sysReac(:)
    integer      , intent(out)   :: ierr

    !! --- Logic section ---

    if (associated(samData%mpreac)) then

       call ADDNV1(nodV(1)         , samData%ttcc(1)  , samData%mpar(1)  , &
            &      samData%madof(1), samData%meqn(1)  , samData%mpmceq(1), &
            &      samData%mmceq(1), samData%mpreac(1), &
            &      abs(inod)       , sign(1,inod)     , getErrorFile()   , &
            &      sysV(1)         , sysReac(1)       , ierr)

    else

       call ADDNV (nodV(1)         , samData%ttcc(1)  , samData%mpar(1)  , &
            &      samData%madof(1), samData%meqn(1)  , samData%mpmceq(1), &
            &      samData%mmceq(1), abs(inod)        , getErrorFile()   , &
            &      sysV(1)         , ierr)

    end if

    if (ierr < 0) then
       call reportError (error_p, &
            &            'Failed to add node vector into system vector', &
            &            addString='AsmExtensionModule::csAddNV')
    else
       ierr = 0
    end if

  end subroutine csAddNV


  !!============================================================================
  !> @brief Extracts an element vector from a system vector.
  !>
  !> @param[in] samData Data for managing system matrix assembly
  !> @param[in] iel Element index
  !> @param[out] elV Extracted element vector.
  !> @param[in] sysV The system vector to extract from
  !> @param[out] nedof Total number of DOFs in element
  !>
  !> @callergraph
  !>
  !> @author Karl Erik Thoresen
  !>
  !> @date Dec 1998

  subroutine csExtractElV (samData,iel,elV,sysV,nedof)

    use SamModule        , only : SamType, dp
    use reportErrorModule, only : internalError

    type(SamType)   , intent(in)  :: samData
    integer         , intent(in)  :: iel
    real(dp)        , intent(out) :: elV(:)
    real(dp)        , intent(in)  :: sysV(:)
    integer,optional, intent(out) :: nedof

    !! Local variables
    integer :: NELDOF

    !! --- Logic section ---

    call EXTEV (sysV(1),samData%madof(1),samData%mpmnpc(1),samData%mmnpc(1), &
         &      iel,elV(1),NELDOF)

    if (NELDOF > size(elV)) then
       NELDOF = internalError('AsmExtensionModule::csExtractElV: '// &
            &                 'Overflow of local array elV')
    end if

    if (present(nedof)) nedof = NELDOF

  end subroutine csExtractElV


  !!============================================================================
  !> @brief Finds the matrix of element equation numbers for an element.
  !>
  !> @param[in] MADOF Matrix of accumulated DOFs
  !> @param[in] MEQN Matrix of equation numbers
  !> @param[in] MPMNPC Matrix of pointers to MNPCs
  !> @param[in] MMNPC Matrix of nodal point correspondances
  !> @param[in] IPAR18 Value of @a mpar(18)
  !> @param[in] IEL Element index
  !> @param[out] NEDOF Total number of DOFs in element
  !> @param[out] MEEN Matrix of equation numbers associated with the element
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 12 Feb 2016

  subroutine csGetElmEq (MADOF,MEQN,MPMNPC,MMNPC,IPAR18,IEL,NEDOF,MEEN)

    integer          , intent(in)  :: MADOF(:), MEQN(:), MPMNPC(:), MMNPC(:)
    integer          , intent(in)  :: IPAR18, IEL
    integer          , intent(out) :: NEDOF
    integer, optional, intent(out) :: MEEN(:)

    !! Local variables
    integer :: i, ii, in, is, ie, IP, NENOD, NNDOF

    !! --- Logic section ---

    IP    = MPMNPC(IEL)
    NENOD = MPMNPC(IEL+1) - IP
    if (IPAR18 > 0) then
       !! Assume that NNDOF = MADOF(I+I) - MADOF(I) for node I
       NNDOF = madof(size(madof))
    else if (present(MEEN)) then
       !! Assume that NEDOF = NNDOF*NENOD, and NNDOF is constant
       NNDOF = size(MEEN)/NENOD
    else
       NEDOF = -1
       return
    end if

    NEDOF = 0
    do ii = 1, NENOD
       in = MMNPC(ip+ii-1)
       is = MADOF(in)
       ie = MIN(is+NNDOF,MADOF(in+1)) - 1
       do i = is, ie
          NEDOF = NEDOF + 1
          if (present(MEEN)) then
             MEEN(NEDOF) = MEQN(i)
          end if
       end do
    end do

  end subroutine csGetElmEq


#ifdef FT_HAS_MKL
  !!============================================================================
  !> @brief Calculates the system matrix sparsity pattern for an FE model.
  !>
  !> @param[in] samData Data for managing system matrix assembly
  !> @param sparse sparse Sparse data structure for the Pardiso equation solver
  !> @param[out] ierr Error flag
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 12 Feb 2016

  subroutine csPreAssSparse (samData,sparse,IERR)

    use SamModule          , only : SamType
    use SysMatrixTypeModule, only : PardisoStorageType
    use SparseMatrixUtils  , only : spr_init, spr_preass
    use SparseMatrixUtils  , only : spr_getnnz, spr_getpattern
    use allocationModule   , only : reAllocate
    use scratchArrayModule , only : getIntegerScratchArray
    use reportErrorModule  , only : reportError, debugFileOnly_p

    type(SamType)           , intent(in)    :: samData
    type(PardisoStorageType), intent(inout) :: sparse
    integer                 , intent(out)   :: ierr

    !! Local variables
    integer          :: iel, nedof
    integer, pointer :: iWork(:)

    character(len=39), parameter :: prgnam_p = &
         'AsmExtensionModule::csPreAssSparse'

    !! --- Logic section ---

    ierr = 0
    call spr_init (samData%neq)

    do iel = 1, samData%nel

       call csGetElmEq (samData%madof, samData%meqn, samData%mpmnpc, &
            &           samData%mmnpc, 1, iel, nedof)
       if (nedof < 1) cycle

       iWork => getIntegerScratchArray(nedof,ierr)
       if (ierr < 0) then
          call reportError (debugFileOnly_p,prgnam_p)
          return
       end if

       call csGetElmEq (samData%madof, samData%meqn, samData%mpmnpc, &
            &           samData%mmnpc, samData%mpar(18), iel, nedof, iWork)
       call spr_preass (nedof, iWork(1), samData%meqn(1), &
            &           samData%mpmceq(1), samData%mmceq(1))

    end do

    call reAllocate (prgnam_p,sparse%IA,samData%neq+1,ierr)
    call reAllocate (prgnam_p,sparse%JA,spr_getnnz(),ierr)
    if (ierr < 0) return

    call spr_getpattern (sparse%IA(1),sparse%JA(1))

  end subroutine csPreAssSparse
#endif

end module AsmExtensionModule
