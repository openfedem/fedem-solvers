!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file solExtensionModule.f90
!> @brief Subroutines for solution of linear equation systems.

!!==============================================================================
!> @brief Module with subroutines for solution of linear equation systems.
!>
!> @details This module contains subroutines for solving the linear system of
!> equations that arise in Fedem solver applications. This includes both direct
!> solution for a set of right-hand-side vectors, the reduction of a system of
!> equations such that only external DOFs are retained, as well as eigenvalue
!> solutions. The module defines an interface between the Fedem solver
!> applications and the various 3rd-party linear algebra packages that
!> are being used.

module SolExtensionModule

  use kindModule, only : dp

  implicit none

  real(dp), parameter, private :: epsilon_p = 1.0e-12_dp !< Default tolerance

  !> @brief Solves the linear system of equations Ax = B.
  !> @sa solextensionmodule::cssolveb.
  interface csSolve
     module procedure csSolveA
     module procedure csSolveB
  end interface

  private :: reportSingularities, csSolveA, csSolveB


contains

  !!============================================================================
  !> @cond FULL_DOC
  !> @brief Prints an error message on detection of singular equation systems.

  subroutine reportSingularities (iopSing,diag,neq,nsing,meqnInErr,eqnInErr,err)
    use sprKindModule, only : ik
    integer    , intent(in)    :: iopSing
    real(dp)   , intent(in)    :: diag(:)
    integer(ik), intent(in)    :: neq
    integer(ik), intent(out)   :: nsing
    integer(ik), intent(inout) :: err
    integer(ik), intent(out), optional :: meqnInErr(:), eqnInErr
    integer(ik) :: nrsing
    nsing = 0_ik
    nrsing = 0_ik
    if (present(meqnInErr)) then
       !! Report all singular equations encountered
       call SPRF3S (neq,nsing,nrsing,meqnInErr(1),diag(1))
    else if (present(eqnInErr)) then
       nsing = -1_ik ! Report only the first singular equation
       call SPRF3S (neq,nsing,nrsing,eqnInErr,diag(1))
    end if
    if (abs(iopSing) > 1) then
       err = nsing ! Ignore all singularities
    else if (iopSing /= 0 .and. nsing == nrsing) then
       err = nsing ! All singularities are true ones, ignore them all
    end if
  end subroutine reportSingularities
  !> @endcond


  !!============================================================================
  !> @brief Solves the linear system of equations Ax = B.
  !> @details This version accepts the right-hand-side vectors @a rhs to be
  !> provided as a 2D neq&times;nrhs matrix instead of a 1D array.
  !> @sa solextensionmodule::cssolveb for detailed parameter description.

  subroutine csSolveA (iop,iopSing,sysMat,rhs,lpu,ierr, &
       &               eqnInErr,meqnInErr,tolFactorize,scaleOnSing,neq1)
    use sprKindModule      , only : ik
    use SysMatrixTypeModule, only : SysMatrixType
    integer              , intent(in)    :: iop, iopSing
    type(SysMatrixType)  , intent(inout) :: sysMat
    real(dp)             , intent(inout) :: rhs(:,:)
    integer              , intent(in)    :: lpu
    integer              , intent(out)   :: ierr
    integer(ik), optional, intent(out)   :: eqnInErr, meqnInErr(:)
    real(dp)   , optional, intent(in)    :: tolFactorize, scaleOnSing
    integer    , optional, intent(in)    :: neq1
    integer           :: nrhs
    real(dp), pointer :: rVec(:)
    nrhs = size(rhs,2)
    rVec => castToVec(rhs,size(rhs,1)*nrhs)
    call csSolveB (iop,iopSing,sysMat,rVec,nrhs,lpu,ierr, &
         &         eqnInErr,meqnInErr,tolFactorize,scaleOnSing,neq1)
  contains
    !> @cond FULL_DOC
    !> @brief Matrix-to-vector type casting.
    function castToVec (array,ndim)
      integer , intent(in)         :: ndim
      real(dp), intent(in), target :: array(ndim)
      real(dp), pointer            :: castToVec(:)
      castToVec => array
    end function castToVec
    !> @endcond
  end subroutine csSolveA


  !!============================================================================
  !> @brief Solves the linear system of equations Ax = B.
  !>
  !> @param[in] iop Calculation option (see below)
  !> @param[in] iopSing Singularity handling options (see below)
  !> @param sysMat The coefficient matrix @b A
  !> @param rhs The right-hand-side vector(s) @b B, solution on output
  !> @param[in] nrhs Number of right-hand-side vectors
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !> @param[out] eqnInErr Equation number of the first singularity detected
  !> @param[out] meqnInErr Equation numbers of all singularities detected
  !> @param[in] tolFactorize Singularity tolerance
  !> @param[in] scaleOnSing Scaling factor for singular pivot elements
  !> @param[in] neq1 Number of equations to solve for
  !>
  !> @details Basic equations: Ax = B, LDL'x = B, LDy = B, L'x = y.
  !> The value on the input argument @a iop defines which solution step to do:
  !> - = 1 : Factorization
  !>         (compute L and D in A = LDL')
  !> - = 2 : Factorization and forward reduction
  !>         (compute L and D in A = LDL' and y in LDy = B)
  !> - = 3 : Factorization, forward reduction and back substitution
  !>         (compute L and D in A = LDL', y in LDy = B and x in L'x = y)
  !> - = 4 : Forward reduction and back substitution
  !>         (compute y in LDy = B and x in L'x = y)
  !> - = 5 : Forward reduction only
  !>         (compute y in LDy = B)
  !> - = 6 : Back substitution only
  !>         (compute x in L'x = y)
  !>
  !> The singularity handling flag @a iopSing is interpreted as follows:
  !> - &gt; 1 : Do a solve with modified matrix to override singularities
  !> -   =  1 : Same as for &gt; 1, but override only the true singularities
  !> -   =  0 : Abort if the matrix is singular
  !> - &lt; 0 : Same as for &gt; 0, but suppress all the singular equations
  !>
  !> @author Bjorn Haugen
  !> @date Sep 1998
  !>
  !> @author Knut Morten Okstad
  !> @date Jan 2003

  subroutine csSolveB (iop,iopSing,sysMat,rhs,nrhs,lpu,ierr, &
       &               eqnInErr,meqnInErr,tolFactorize,scaleOnSing,neq1)

    use sprKindModule      , only : i4, ik
    use SysMatrixTypeModule, only : SysMatrixType, checkGSFinfo, GSFinfo
    use SysMatrixTypeModule, only : SparseStorageType, denseMatrix_p
    use SysMatrixTypeModule, only : skylineMatrix_p, sparseMatrix_p, outOfCore_p
#ifdef FT_HAS_MKL
    use SysMatrixTypeModule, only : pardiso_p
    use MKL_pardiso        , only : pardiso
#endif
    use FELinearSolver     , only : FEFactorize, FESolve
    use manipMatrixModule  , only : MatrixPointToArray_real
    use scratchArrayModule , only : getIntegerScratchArray, getRealScratchArray
    use scratchArrayModule , only : getRealScratchMatrix, releaseScratchArrays
    use allocationModule   , only : reAllocate
    use reportErrorModule  , only : reportError, internalError
    use reportErrorModule  , only : debugFileOnly_p, error_p, errorFileOnly_p
    use FFaProfilerInterface,only : ffa_startTimer, ffa_stopTimer

    integer            , intent(in)    :: iop, iopSing
    type(SysMatrixType), intent(inout) :: sysMat
    real(dp)           , intent(inout) :: rhs(:)
    integer            , intent(in)    :: nrhs, lpu
    integer            , intent(out)   :: ierr
    integer(ik), optional, intent(out) :: eqnInErr, meqnInErr(:)
    real(dp), optional , intent(in)    :: tolFactorize, scaleOnSing
    integer , optional , intent(in)    :: neq1

    !! Local variables
    integer           :: nWork, nnd = 0, ndof1, neq
#ifdef FT_HAS_MKL
    integer           :: phase
#endif
    integer(ik)       :: lerr, neqs, nneg
    integer , pointer :: iWork(:)
    real(dp), pointer :: rWork(:), xWork(:,:), rhMat(:,:)
    real(dp)          :: tolerance(4)

    type(SparseStorageType), pointer :: sparse

    character(len=27), parameter :: prgnam_p = 'SolExtensionModule::csSolve'

    external LMVY18, LMMY18

    !! --- Logic section ---

    call ffa_startTimer ('csSolve')

    lerr = 0_ik
    if (present(meqnInErr) .and. iop <= 3) meqnInErr(1) = 0_ik

    if (present(tolFactorize)) then
       tolerance(1) = tolFactorize
    else
       tolerance(1) = epsilon_p
    end if
    if (present(scaleOnSing)) then
       tolerance(4) = scaleOnSing
    else
       tolerance(4) = 10.0_dp
    end if

    if (present(neq1)) then
       ndof1 = neq1
    else
       ndof1 = sysMat%dim
    end if

    select case (sysMat%storageType)

    case (skylineMatrix_p) ! Use the SAM skyline solver

       neq = sysMat%skyline%neq
       if (present(meqnInErr) .or. iopSing /= 0) then

          !! Allocate real work array
          rWork => getRealScratchArray(2*neq,ierr)
          if (ierr < 0) goto 999

          !! Modify the matrix during the factorization,
          !! if needed, to avoid singularities
          call skySolve (sysMat%value(1), rhs(1), sysMat%skyline%msky(1), &
               &         tolerance(1), tolerance(4), neq, ndof1, nrhs, lpu, &
               &         sign(iop,iopSing), rWork(1), nnd, ierr)
          if (ierr == -3) then
             call reportSingularities (abs(iopSing),rWork,int(neq,ik),nneg, &
                  &                    meqnInErr,eqnInErr,lerr)
             nnd = int(nneg)
          end if

       else

          !! Abort on the first singularity encountered
          call SKYSOL (sysMat%value(1), rhs(1), sysMat%skyline%msky(1), &
               &       tolerance(1), neq, ndof1, nrhs, lpu, iop, nnd, ierr)
          if (ierr <= -2 .and. present(eqnInErr)) eqnInErr = int(nnd,ik)

       end if
       if (ierr < 0 .and. lerr == 0_ik) goto 998

    case (denseMatrix_p) ! Use the LAPack dense matrix solver

       neq = sysMat%dim

       if (iop <= 3) then ! Factorization

          call reAllocate (prgnam_p,sysMat%ipiv,neq,ierr)
          if (ierr < 0) goto 999

          call DGETRF (ndof1,ndof1,sysMat%value(1),neq,sysMat%ipiv(1),ierr)
          if (ierr /= 0) then
             if (ierr > 0) then
                if (present(meqnInErr)) then
                   meqnInErr(1) = int(ierr,ik)
                else if (present(eqnInErr)) then
                   eqnInErr = int(ierr,ik)
                end if
                ierr = -4
             else
                ierr = ierr - 10
             end if
             call reportError (errorFileOnly_p,'LAPACK::DGETRF')
             goto 999
          end if

       end if
       if (iop >= 2 .and. iop < 6) then ! Forward reduction

          call DLASWP (nrhs,rhs(1),neq,1,ndof1,sysMat%ipiv(1),1)
          call DTRSM ('L','L','N','U',ndof1,nrhs, &
               &      1.0_dp,sysMat%value(1),neq,rhs(1),neq)

       end if
       if (iop >= 3 .and. iop /= 5) then ! Back substitution

          call DTRSM ('L','U','N','N',ndof1,nrhs, &
               &      1.0_dp,sysMat%value(1),neq,rhs(1),neq)

       end if

    case (sparseMatrix_p) ! Use the SAM sparse matrix solver (in-core only)

       sparse => sysMat%sparse
       if (ndof1 == sysMat%dim-int(sparse%mspar(54))) then
          sparse%mspar(57) = 1_ik ! We're going to solve for the status-1 DOFs
       else if (ndof1 == sysMat%dim) then
          sparse%mspar(57) = 0_ik ! We're going to solve for all DOFs
       else
          ierr = internalError(prgnam_p//': Invalid input arguments')
          goto 990
       end if
       neqs = sparse%mspar(8)

       if (iop <= 3) then ! Factorization

          !! Allocate integer work array
          nWork = int(max(sparse%mspar(13),sparse%mspar(14)+1_ik))
          if (ik > i4) nWork = nWork*2 ! We are using 64-bit integers
          iWork => getIntegerScratchArray(nWork,ierr)
          if (ierr < 0) goto 999

          !! Allocate real work array
          neq   = int(neqs)
          nWork = int(sparse%mspar(17))
          if (present(meqnInErr) .or. iopSing /= 0) then
             nWork = nWork + 2*neq
             !! Ensure that the work array is large enough also for the SPRFS1
             !! call, such that it will not be reallocated before that call.
             !! Because we may need the 2*neq first values for the SPRF3R call.
             if (iopSing < 0) nWork = max(nWork,int(sparse%mspar(14)))
          end if
          rWork => getRealScratchArray(nWork,ierr)
          if (ierr < 0) goto 999

          !! Perform the factorization
          if (present(meqnInErr) .or. iopSing /= 0) then
             !! Modify the matrix during the factorization,
             !! if needed, to avoid singularities
             call SPRFC3 (sparse%mspar(1), sparse%mtrees(1), sparse%msifa(1), &
               &          sysMat%value(1), tolerance(1), iWork(1), &
               &          rWork(1+2*neq), rWork(1), int(lpu,ik), lerr, &
               &          LMVY18, LMMY18)
             if (lerr == -3_ik .and. sparse%mspar(27) < 0_ik) then
                call reportSingularities (abs(iopSing),rWork,neqs,nneg, &
                     &                    meqnInErr,eqnInErr,lerr)
                nnd = int(nneg)
             else if (lerr == 3_ik .and. present(meqnInErr)) then
                !! Report all negative pivots encountered
                lerr = 0_ik
                call SPRF3T (neqs,lerr,meqnInErr(1),rWork(1))
             end if
          else
             !! Abort on the first singularity encountered
             call SPRFC1 (sparse%mspar(1), sparse%mtrees(1), sparse%msifa(1), &
                  &       sysMat%value(1), tolerance(1), iWork(1), rWork(1), &
                  &       int(lpu,ik), lerr, LMVY18, LMMY18)
             if (lerr == -3_ik .and. present(eqnInErr)) then
                eqnInErr = sparse%mspar(27)
             end if
          end if
          if (lerr < 0_ik) goto 997

       end if
       sparse%mspar(57) = 2_ik*sparse%mspar(57)
       if (iop >= 2) then ! Forward reduction and backward substitution

          !! Allocate integer work array
          nWork = int(sparse%mspar(13))
          if (ik > i4) nWork = nWork*2 ! We are using 64-bit integers
          iWork => getIntegerScratchArray(nWork,ierr)
          if (ierr < 0) goto 999

          !! Allocate real work array
          nWork = int(sparse%mspar(14))
          rWork => getRealScratchArray(nWork,ierr)
          if (ierr < 0) goto 999

          if (iop <= 5) then

             if (iopSing < 0 .and. nnd > 0) then
                !! Suppress the singular (or near singular) equations.
                !! CAUTION: This call will only work when rWork has not been
                !! altered after the SPRFC3 call. Thus, this can not be used in
                !! e.g. a modified Newton procedure, where the matrix is not
                !! factorized before each solve.
                call SPRF3R (neqs,int(nrhs,ik),rWork(1),rhs(1),int(lpu,ik))
             end if

             !! Perform the forward reduction
             call SPRFS1 (sparse%mspar(1), sparse%mtrees(1), sparse%msifa(1), &
                  &       sysMat%value(1), rhs(1), int(sysMat%dim,ik), &
                  &       int(nrhs,ik), iWork(1), rWork(1), int(lpu,ik), lerr)
             if (lerr < 0_ik) goto 997

          end if
          if (iop >= 3 .and. iop /= 5) then

             !! Perform the back substitution
             call SPRBS1 (sparse%mspar(1), sparse%mtrees(1), sparse%msifa(1), &
                  &       sysMat%value(1), rhs(1), int(sysMat%dim,ik), &
                  &       int(nrhs,ik), rWork(1), int(lpu,ik), lerr)
             if (lerr < 0_ik) goto 997

             nnd = 0 ! Reset nnd to prevent using rWork above when not allocated
             call releaseScratchArrays

          end if
       end if
       sparse%mspar(57) = 0_ik

    case (outOfCore_p) ! Use the DnV GSF solver with out-of-core possibility

       if (ndof1 /= sysMat%dim) then
          ierr = internalError(prgnam_p//': Status-1 DOF solve not implemented')
          goto 990
       end if

       if (iop <= 3 .and. .not. sysMat%gsf%isFactorized) then ! Factorization

          call ffa_startTimer ('FEFactorize')
          call FEFactorize (sysMat%gsf%Msg, sysMat%gsf%S, sysMat%gsf%T, GSFinfo)
          call ffa_stopTimer ('FEFactorize')
          ierr = checkGSFinfo(GSFinfo,sysMat%gsf)
          if (ierr < 0) goto 998
          if (iop <= 1) goto 990

          sysMat%gsf%isFactorized = .true.
       end if

       !! Allocate real work array
       neq   =  sysMat%dim
       rhMat => MatrixPointToArray_real(rhs,neq,nrhs)
       xWork => getRealScratchMatrix(neq,nrhs,ierr)
       if (ierr < 0) goto 999

       call ffa_startTimer ('FESolve')
       if (iop == 2 .or. iop == 5) then ! Forward reduction

          call FESolve (sysMat%gsf%Msg, 'FPEA', neq, nrhs, &
               &        sysMat%gsf%Q, sysMat%gsf%T, rhMat, xWork, GSFinfo)

       else if (iop == 6) then ! Back substitution

          call FESolve (sysMat%gsf%Msg, 'BPXA', neq, nrhs, &
               &        sysMat%gsf%Q, sysMat%gsf%T, rhMat, xWork, GSFinfo)

       else ! Forward reduction and back substitution

          call FESolve (sysMat%gsf%Msg, 'SPBA', neq, nrhs, &
               &        sysMat%gsf%Q, sysMat%gsf%T, rhMat, xWork, GSFinfo)

       end if
       call ffa_stopTimer ('FESolve')
       ierr = checkGSFinfo(GSFinfo,sysMat%gsf)
       if (ierr < 0) goto 998

       call DCOPY (neq*nrhs,xWork(1,1),1,rhs(1),1)

#ifdef FT_HAS_MKL
    case (pardiso_p)

       if (iop <= 1) then
          phase = 22
       else if (iop == 2 .or. iop == 3) then
          phase = 23
       else if (iop == 5) then
          phase = 331
       else if (iop == 6) then
          phase = 333
       else
          phase = 33
       end if

       !! Allocate real work array
       neq = sysMat%dim
       rWork => getRealScratchArray(neq*nrhs,ierr)
       if (ierr < 0) goto 999

       call pardiso (sysMat%pardiso%pt, 1, 1, &
            &        sysMat%pardiso%mtype, phase, neq, &
            &        sysMat%value, sysMat%pardiso%ia, sysMat%pardiso%ja, &
            &        sysMat%ipiv, nrhs, sysMat%pardiso%iparm, 0, &
            &        rhs, rWork, ierr)
       if (ierr < 0) goto 998

       if (sysMat%pardiso%iparm(6) == 0 .and. iop > 1) then
          call DCOPY (neq*nrhs,rWork(1),1,rhs(1),1)
       end if
#endif

    case default

       ierr = internalError(prgnam_p//': Illegal matrix type')
       goto 990

    end select

    if (lerr > 0_ik) then
       ierr = int(lerr)
    else
       ierr = 0
    end if

990 call ffa_stopTimer ('csSolve')
    return

997 ierr = int(lerr)
998 call reportError (error_p,'Failed to solve equation system',IERR=ierr)
999 call reportError (debugFileOnly_p,prgnam_p)
    goto 990

  end subroutine csSolveB


  !!============================================================================
  !> @brief Performs the superelement reduction for the equation system.
  !>
  !> @param[in] iopSing Singularity handling options (see below)
  !> @param sysMat The system stiffness matrix
  !> @param[in] sysRhs The system right-hand-side vector(s)
  !> @param[out] supElMat The reduced stiffness matrix, @b K
  !> @param[out] supElRhs The reduced right-hand-side vector(s), @b f
  !> @param[in] meqn1 Matrix of status 1 (internal) equation numbers
  !> @param[in] meqn2 Matrix of status 2 (external) equation numbers
  !> @param[in] nrhs Number of right-hand-side vectors
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !> @param[out] Bmatrix Displacement recovery matrix
  !> @param[out] eqnInErr Equation number of the first singularity detected
  !> @param[out] meqnInErr Equation numbers of all singularities detected
  !> @param[in] tolFactorize Singularity tolerance
  !> @param[in] scaleOnSing Scaling factor for singular pivot elements
  !>
  !> @details A static condensation of the internal (i) degrees of freedom
  !> is performed, retaining only the external (e) degrees of freedom, i.e.,
  !> @verbatim
  !>    |Kii Kie|*|vi| = |fi|    ====>   K*ve = f
  !>    |Kei Kee| |ve|   |fe|
  !>     (sysMat)      (sysRhs)
  !> @endverbatim
  !> where
  !> @verbatim
  !>    K = Kee - Kei*(Kii^-1)*Kie    (returned in supElMat)
  !>    f = fe  - Kei*(Kii^-1)*fi     (returned in supElRhs)
  !> @endverbatim
  !>
  !> If the @a Bmatrix is present, it is also returned
  !> @verbatim
  !>    B = -(Kii^-1)*Kie
  !> @endverbatim
  !>
  !> The singularity handling flag @a iopSing is interpreted as follows:
  !> - &gt; 1 : Do a solve with modified matrix to override singularities
  !> -   =  1 : Same as for &gt; 1, but override only the true singularities
  !> -   =  0 : Abort if the matrix is singular
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !> @date Nov 1998
  !>
  !> @author Knut Morten Okstad
  !> @date Jan 2003

  subroutine csSolveSupEl (iopSing,sysMat,sysRhs,supElMat,supElRhs, &
       &                   meqn1,meqn2,nrhs,lpu,ierr, &
       &                   Bmatrix,eqnInErr,meqnInErr,tolFactorize,scaleOnSing)

    use sprKindModule      , only : i4, ik, castInt
    use SysMatrixTypeModule, only : SysMatrixType, checkGSFinfo, GSFinfo
    use SysMatrixTypeModule, only : SparseStorageType, GSFStorageType
    use SysMatrixTypeModule, only : extractDiagonal
    use SysMatrixTypeModule, only : skylineMatrix_p, sparseMatrix_p, outOfCore_p
    use DiskMatrixModule   , only : DiskMatrixType
    use DiskMatrixModule   , only : dmSetValue, dmGetSwapSecPtr, dmSize
    use FELinearSolver     , only : FEFindIsolatedZeroPivots
    use FELinearSolver     , only : FEInsertPivots, FEPermuteIndices
    use FELinearSolver     , only : FEFactorize, FESolve, FEGetSuperelementOrder
    use FELinearSolver     , only : FEExtractL22, FEExtractB12, FEExtractS
    use scratchArrayModule , only : getIntegerScratchArray, getRealScratchArray
    use scratchArrayModule , only : getRealScratchMatrix, releaseScratchArrays
    use progressModule     , only : writeProgress
    use reportErrorModule  , only : reportError, internalError
    use reportErrorModule  , only : debugFileOnly_p, error_p, warning_p
    use FFaProfilerInterface,only : ffa_startTimer, ffa_stopTimer

    integer             , intent(in)              :: iopSing
    type(SysMatrixType) , intent(inout)           :: sysMat
    real(dp)            , intent(in)              :: sysRhs(:,:)
    real(dp)            , intent(out)             :: supElMat(:,:)
    real(dp)            , intent(out)  , target   :: supElRhs(:,:)
    integer             , intent(in)              :: meqn1(:), meqn2(:)
    integer             , intent(in)              :: nrhs, lpu
    integer             , intent(out)             :: ierr
    type(DiskMatrixType), intent(inout), optional :: Bmatrix
    integer(ik)         , intent(out)  , optional :: eqnInErr, meqnInErr(:)
    real(dp)            , intent(in)   , optional :: tolFactorize, scaleOnSing

    !! Local variables
    integer             :: i, iop, idof, jdof, nWork, neq1, neq2, neq, nnd, nnn
    integer(ik)         :: lerr, neqs, nneg
    integer , pointer   :: iWork(:)
    real(dp), pointer   :: rWork(:), xWork(:,:), vecNeq(:)
    real(dp)            :: tolerance(4)
    real(dp), parameter :: zeroTol_p = 1.0e-127_dp

    type(SparseStorageType), pointer :: sparse

    character(len=32), parameter :: prgnam_p = &
         'SolExtensionModule::csSolveSupEl'

    external LMVY18, LMMY18

    !! --- Logic section ---

    call ffa_startTimer ('csSolveSupEl')

    lerr = 0_ik
    if (present(meqnInErr)) meqnInErr(1) = 0_ik

    neq1 = size(meqn1)
    neq2 = size(meqn2)
    neq  = neq1 + neq2
    neqs = neq
    if (neq /= sysMat%dim) then
       ierr = internalError(prgnam_p//': Invalid input arguments NEQ1, NEQ2')
       goto 990
    end if

    if (present(tolFactorize)) then
       tolerance(1) = tolFactorize
    else
       tolerance(1) = epsilon_p
    end if
    if (present(scaleOnSing)) then
       tolerance(4) = scaleOnSing
    else
       tolerance(4) = 10.0_dp
    end if

    select case (sysMat%storageType)

    case (skylineMatrix_p) ! Use the SAM skyline solver

       nWork = nrhs
       if (present(meqnInErr) .or. iopSing > 0) nWork = nWork + 2
       if (nWork > 0) then
          xWork => getRealScratchMatrix(neq,nWork,ierr)
          if (ierr < 0) goto 999
       else
          xWork => supElRhs ! Unused, but should point to something
       end if

       if (nrhs > 0) call DCOPY (neq*nrhs,sysRhs(1,1),1,xWork(1,1),1)

       iop = 1 + min(1,nrhs)
       if (present(meqnInErr) .or. iopSing > 0) then

          !! Modify the matrix during the factorization,
          !! if needed, to avoid singularities
          call skySolve (sysMat%value(1), xWork(1,1), sysMat%skyline%msky(1),&
               &         tolerance(1), tolerance(4), neq, neq1, nrhs, &
               &         lpu, iop, xWork(1,nrhs+1), nnd, ierr)
          if (ierr == -3) then
             call reportSingularities (iopSing,xWork(:,nrhs+1),neqs,nneg, &
                  &                    meqnInErr,eqnInErr,lerr)
             nnd = int(nneg)
          end if
       else
          !! Abort on the first singularity encountered
          call SKYSOL (sysMat%value(1), xWork(1,1), sysMat%skyline%msky(1), &
               &       tolerance(1), neq, neq1, nrhs, &
               &       lpu, iop, nnd, ierr)
          if (ierr <= -2 .and. present(eqnInErr)) eqnInErr = int(nnd,ik)
       end if
       if (ierr < 0 .and. lerr == 0_ik) then
          goto 998
       else if (ierr == 3) then
          call warningOnNegativePivots (nnd)
       end if

       call EXSEM (sysMat%value(1), sysMat%skyline%msky(1), &
            &      neq, neq2, lpu, supElMat(1,1), ierr)
       if (ierr < 0) goto 999

       if (nrhs > 0) supElRhs = xWork(neq1+1:neq,1:nrhs)
       call writeProgress(0.1_dp*neq2+nrhs,real(nrhs+neq2,dp))

       if (present(Bmatrix)) then

          !! Compute the B-matrix for retracking of internal dofs

          vecNeq => getRealScratchArray(neq,ierr)
          if (ierr < 0) goto 999

          iop = 6
          do idof = 1, neq2

             vecNeq = 0.0_dp
             vecNeq(meqn2(idof)) = 1.0_dp
             call SKYSOL (sysMat%value(1), vecNeq(1), sysMat%skyline%msky(1), &
                  &       tolerance(1), neq, neq1, 1, lpu, iop, nnd, ierr)
             if (ierr < 0) goto 998

             do i = 1, neq1
                call dmSetValue (Bmatrix, vecNeq(meqn1(i)), i, idof, ierr)
                if (ierr < 0) goto 999
             end do

             call writeProgress(0.1_dp*neq2+nrhs+0.9_dp*idof,real(nrhs+neq2,dp))
          end do

       end if

    case (sparseMatrix_p) ! Use the SAM sparse matrix solver

       sparse => sysMat%sparse
       if (neqs /= sparse%mspar(8) .or. int(neq2,ik) /= sparse%mspar(54)) then
          ierr = internalError(prgnam_p//': Invalid input arguments NEQ, NEQ2')
          goto 990
       end if

       if (sparse%mspar(1) < 5_ik) then ! Check if matrix already is factorized

          nWork = int(max(sparse%mspar(13),sparse%mspar(14)+1_ik))
          if (ik > i4) nWork = nWork*2 ! We are using 64-bit integers
          iWork => getIntegerScratchArray(nWork,ierr)
          if (ierr < 0) goto 999

          nWork = int(sparse%mspar(17))
          if (present(meqnInErr) .or. iopSing > 0) nWork = nWork + 2*neq
          rWork => getRealScratchArray(nWork,ierr)
          if (ierr < 0) goto 999

          !! Perform the factorization
          sparse%mspar(57) = 1_ik
          if (present(meqnInErr) .or. iopSing > 0) then
             !! Modify the matrix during the factorization,
             !! if needed, to avoid singularities
             call SPRFC3 (sparse%mspar(1), sparse%mtrees(1), sparse%msifa(1), &
                  &       sysMat%value(1), tolerance(1), iWork(1), &
                  &       rWork(1+2*neq), rWork(1), int(lpu,ik), lerr, &
                  &       LMVY18, LMMY18)
             if (lerr == -3_ik .and. sparse%mspar(27) < 0_ik) then
                call reportSingularities (iopSing,rWork,neqs,nneg, &
                     &                    meqnInErr,eqnInErr,lerr)
                nnd = int(nneg)
             end if
          else
             !! Abort on the first singularity encountered
             call SPRFC1 (sparse%mspar(1), sparse%mtrees(1), sparse%msifa(1), &
                  &       sysMat%value(1), tolerance(1), iWork(1), rWork(1), &
                  &       int(lpu,ik), lerr, LMVY18, LMMY18)
             if (lerr == -3_ik .and. present(eqnInErr)) then
                eqnInErr = sparse%mspar(27)
             end if
          end if
#ifdef FT_DEBUG
          write(lpu,600) tolerance(1:3)
600       format(/5X,'Factorization threshold :',1PE12.5, &
               & /5X,'Max diagonal decay      :',1PE12.5, &
               & /5X,'TRACE(A)/MIN(D(I))      :',1PE12.5 )
#endif
          if (lerr < 0_ik) then
             goto 997
          else if (lerr == 3_ik) then
             call warningOnNegativePivots (int(sparse%mspar(26)))
          end if

       end if

       if (ik > i4) then ! We are using 64-bit integers
          iWork => getIntegerScratchArray(neq*2,ierr)
       else
          iWork => getIntegerScratchArray(neq,ierr)
       end if
       if (ierr < 0) goto 999

       !! Extract the superelement matrix
       call SPREXK (sparse%mspar(1), sparse%mtrees(1), sparse%msifa(1), &
            &       sysMat%value(neq+1), supElMat(1,1), iWork(1))
       call writeProgress(1.0_dp,10.0_dp)

       sparse%mspar(57) = 2_ik
       if (nrhs > 0) then

          nWork = int(max(neq,sparse%mspar(13)))
          if (ik > i4) nWork = nWork*2 ! We are using 64-bit integers
          iWork => getIntegerScratchArray(nWork,ierr)
          if (ierr < 0) goto 999

          nWork = neq*nrhs + int(sparse%mspar(14))
          rWork => getRealScratchArray(nWork,ierr)
          if (ierr < 0) goto 999

          call DCOPY (neq*nrhs,sysRhs(1,1),1,rWork(1),1)

          !! Perform the forward reduction
          if (sparse%mspar(1) == 5_ik .or. sparse%mspar(1) == 7_ik) then
             call SPRFS1 (sparse%mspar(1), sparse%mtrees(1), sparse%msifa(1), &
                  &       sysMat%value(1), rWork(1), neqs, int(nrhs,ik), &
                  &       iWork(1), rWork(neq*nrhs+1), int(lpu,ik), lerr)
          else
             call SPRFS2 (sparse%mspar(1), sparse%mtrees(1), sparse%msifa(1), &
                  &       sysMat%value(1), rWork(1), neqs, int(nrhs,ik), &
                  &       iWork(1), rWork(neq+nrhs+1), int(lpu,ik), lerr)
          end if
          if (lerr < 0_ik) goto 997

          !! Extract the superelement right-hand-side vectors
          call SPREXS (sparse%mspar(1), sparse%mtrees(1), sparse%msifa(1), &
               &       int(nrhs,ik), rWork(1), supElRhs(1,1), iWork(1))
          call writeProgress(0.1_dp*neq2+nrhs,real(nrhs+neq2,dp))

       end if
       if (present(Bmatrix)) then

          !! Compute the B-matrix for retracking of internal dofs

          nWork = int(sparse%mspar(14)) + neq
          rWork => getRealScratchArray(nWork,ierr)
          if (ierr < 0) goto 999

          vecNeq => rWork(sparse%mspar(14)+1_ik:sparse%mspar(14)+neqs)

          !! Perform the back substitution for each external dof
          do idof = 1, neq2

             vecNeq = 0.0_dp
             vecNeq(meqn2(idof)) = 1.0_dp
             if (sparse%mspar(1) == 5_ik .or. sparse%mspar(1) == 7_ik) then
                call SPRBS1 (sparse%mspar(1),sparse%mtrees(1),sparse%msifa(1), &
                     &       sysMat%value(1), vecNeq(1), int(sysMat%dim,ik), &
                     &       1_ik, rWork(1), int(lpu,ik), lerr)
             else
                call SPRBS2 (sparse%mspar(1),sparse%mtrees(1),sparse%msifa(1), &
                     &       sysMat%value(1), vecNeq(1), int(sysMat%dim,ik), &
                     &       1_ik, rWork(1), int(lpu,ik), lerr)
             end if
             if (lerr < 0_ik) goto 997

             do i = 1, neq1
                call dmSetValue (Bmatrix, vecNeq(meqn1(i)), i, idof, ierr)
                if (ierr < 0) goto 999
             end do

             call writeProgress(0.1_dp*neq2+nrhs+0.9_dp*idof,real(nrhs+neq2,dp))
          end do

       end if

    case (outOfCore_p) ! Use the DnV GSF solver

       if (.not. sysMat%gsf%isFactorized) then

          if (present(meqnInErr) .or. iopSing > 0) then

             !! Check for zero and negative pivots before factorization
             rWork => getRealScratchArray(sysMat%dim,ierr)
             if (ierr < 0) goto 999

             nnn = 0
             nullify(iWork)
             call FEFindIsolatedZeroPivots (sysMat%gsf%S,iWork,GSFinfo)
             ierr = checkGSFinfo(GSFinfo,sysMat%gsf)
             if (ierr < 0) goto 998

             if (associated(iWork)) then
                if (iopSing > 0) then
                   !! Isolated zero pivots were found,
                   !! fix them by inserting 1.0 on the diagonal
                   rWork(1:size(iWork)) = 1.0_dp
                   call FEInsertPivots (sysMat%gsf%S,iWork,rWork,GSFinfo)
                   ierr = checkGSFinfo(GSFinfo,sysMat%gsf)
                   if (ierr < 0) goto 998
                end if
                if (present(meqnInErr)) then
                   !! Find external equation numbers of the zero pivots found
                   call FEPermuteIndices ('F',sysMat%gsf%Q,iWork,GSFinfo)
                   ierr = checkGSFinfo(GSFinfo,sysMat%gsf)
                   if (ierr < 0) goto 998
                   nnn = size(iWork)
                   nnd = min(nnn,size(meqnInErr))
                   call castInt (nnd,iWork,meqnInErr)
                end if
                deallocate(iWork)
             end if

             call extractDiagonal (sysMat,rWork,ierr)
             if (ierr < 0) goto 999

             nnd = 0
             do i = 1, sysMat%dim
                if (rWork(i) < -zeroTol_p) then
                   nnd = nnd + 1
#ifdef FT_DEBUG
                   write(lpu,"(5X,'Negative pivot:',I8,1PE12.5)") i,rWork(i)
#endif
                else if (present(meqnInErr) .and. rWork(i) <= zeroTol_p) then
                   if (count(meqnInErr == int(i,ik)) == 0) then
                      nnn = nnn + 1
                      if (nnn <= size(meqnInErr)) meqnInErr(nnn) = int(i,ik)
                   end if
                end if
             end do
             lerr = int(nnn,ik)
             if (present(meqnInErr)) then
                if (nnn+1 <= size(meqnInErr)) meqnInErr(nnn+1) = -lerr
             end if
             if (nnd > 0) call warningOnNegativePivots (nnd)

#ifdef FT_DEBUG
             iop = 1
             idof = 1
             jdof = 1
             tolerance = rWork(1)
             do i = 2, sysMat%dim
                if (rWork(i) < tolerance(1)) then
                   idof = i
                   tolerance(1) = rWork(i)
                else if (rWork(i) > tolerance(3)) then
                   jdof = i
                   tolerance(3) = rWork(i)
                else if (abs(rWork(i)) < abs(tolerance(2))) then
                   iop = i
                   tolerance(2) = rWork(i)
                end if
             end do
             write(lpu,"(/5X,'Smallest pivot:',I8,1PE12.5)") idof,tolerance(1)
             write(lpu,"( 5X,'Zeroest  pivot:',I8,1PE12.5)") iop ,tolerance(2)
             write(lpu,"( 5X,'Largest  pivot:',I8,1PE12.5)") jdof,tolerance(3)
             do i = 1, sysMat%dim
                if (abs(rWork(i)) < epsilon_p) then
                   write(lpu,"(5X,'Small    pivot:',I8,1PE12.5)") i,rWork(i)
                end if
             end do
#endif
          end if

          !! Perform the factorization
          call ffa_startTimer ('FEFactorize')
          call FEFactorize (sysMat%gsf%Msg, sysMat%gsf%S, sysMat%gsf%T, GSFinfo)
          call ffa_stopTimer ('FEFactorize')
          ierr = checkGSFinfo(GSFinfo,sysMat%gsf)
          if (ierr < 0 .and. lerr > 0_ik) ierr = -3
          if (ierr < 0) goto 998

          sysMat%gsf%isFactorized = .true.
       end if

       !! Extract the superelement matrix
       call ffa_startTimer ('FEExtractL22')
       call FEExtractL22 ('A', sysMat%gsf%Q, sysMat%gsf%T, supElMat, GSFinfo)
       call ffa_stopTimer ('FEExtractL22')
       ierr = checkGSFinfo(GSFinfo,sysMat%gsf)
       if (ierr < 0) goto 998

       call writeProgress(1.0_dp,10.0_dp)

       if (nrhs > 0) then

          xWork => getRealScratchMatrix(neq,nrhs,ierr)
          if (ierr < 0) goto 999

          !! Perform the forward reduction
          call ffa_startTimer ('FESolve')
          call FESolve (sysMat%gsf%Msg, 'FPEA', neq, nrhs, &
               &        sysMat%gsf%Q, sysMat%gsf%T, sysRhs, xWork, &
               &        GSFinfo, .true.)
          call ffa_stopTimer ('FESolve')
          ierr = checkGSFinfo(GSFinfo,sysMat%gsf)
          if (ierr < 0) goto 998

          !! Extract the superelement right-hand-side vectors
          call ffa_startTimer ('FEExtractS')
          call FEExtractS (sysMat%gsf%Q, sysMat%gsf%T, nrhs, &
               &           xWork, supElRhs, GSFinfo)
          call ffa_stopTimer ('FEExtractS')
          ierr = checkGSFinfo(GSFinfo,sysMat%gsf)
          if (ierr < 0) goto 998

          call writeProgress(0.1_dp*neq2+nrhs,real(nrhs+neq2,dp))

       end if
       if (present(Bmatrix)) then

          iWork => getIntegerScratchArray(neq,ierr)
          if (ierr < 0) goto 999

          !! Compute the B-matrix for retracking of internal dofs
          nWork = int(dmSize(Bmatrix,3)) ! Number of columns in swap section
          do idof = 1, neq2, nWork

             xWork => dmGetSwapSecPtr(Bmatrix,idof,ierr)
             if (ierr < 0) goto 999

             jdof = min(idof+nWork-1,neq2)
             call ffa_startTimer ('FEExtractB12')
             call FEExtractB12 (sysMat%gsf%Msg, sysMat%gsf%Q, sysMat%gsf%T, &
                  &             xWork, idof, jdof, iWork, GSFinfo)
             call ffa_stopTimer ('FEExtractB12')
             ierr = checkGSFinfo(GSFinfo,sysMat%gsf)
             if (ierr < 0) goto 998

             call writeProgress(0.1_dp*neq2+nrhs+0.9_dp*jdof,real(nrhs+neq2,dp))
          end do

       end if

    case default

       ierr = internalError(prgnam_p//': Illegal matrix type')
       goto 990

    end select

    if (lerr > 0_ik) ierr = int(lerr)

990 call releaseScratchArrays
    call ffa_stopTimer ('csSolveSupEl')
    return

997 ierr = int(lerr)
998 call reportError (error_p,'Failed to reduce the equation system',IERR=ierr)
999 call reportError (debugFileOnly_p,prgnam_p)
    goto 990

  contains

    !> @cond FULL_DOC
    !> @brief Prints a warning when negative pivot elements are detected.
    subroutine warningOnNegativePivots (numneg)
      integer, intent(in) :: numneg
      character(len=32)   :: errMsg
      write(errMsg,"('There are',i6,' negative pivots.')") numneg
      call reportError (warning_p, &
           'The stiffness matrix is not positive definite.',errMsg, &
           'This can be caused by too many ill-shaped finite elements.', &
           'Using this FE part may lead to an instable dynamics simulation.')
    end subroutine warningOnNegativePivots
    !> @endcond

  end subroutine csSolveSupEl


  !!============================================================================
  !> @brief Solves the eigenvalue problem (@b A + &lambda;@b B) * @b u = @b 0.
  !>
  !> @param aMat System matrix @b A
  !> @param bMat System matrix @b B
  !> @param[in] mip Matrix of input parameters, see LANCZ2 or SPRLAN in SAM
  !> @param[out] mop Matrix of output parameters, see LANCZ2 or SPRLAN
  !> @param[in] meqn1 Matrix of status 1 (internal) equation numbers
  !> @param[in] neq Size of the equation system to solve
  !> @param[in] nEigVal Number of eigenvalues to return
  !> @param[in] nEigVec Number of eigenvectors to return
  !> @param[in] tolerance Array of tolerances, see LANCZ2 or SPRLAN
  !> @param[in] shift Shift value
  !> @param[out] eigValues Computed eigenvalues
  !> @param[out] eigVectors Computed eigenvectors
  !> @param[in] ipsw Print switch for debug output
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] meqnInErr Equation numbers of all singularities detected
  !> @param[out] ierr Error flag
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !> @date Oct 1998
  !>
  !> @author Knut Morten Okstad
  !> @date Jan 2003

  subroutine csLanczosEigenSolve (aMat,bMat,mip,mop,meqn1,neq,nEigVal,nEigVec, &
       &                          tolerance,shift,eigValues,eigVectors, &
       &                          ipsw,lpu,meqnInErr,ierr)

    use sprKindModule      , only : i4, ik, castInt
    use SysMatrixTypeModule, only : SysMatrixType, checkGSFinfo, GSFinfo
    use SysMatrixTypeModule, only : SparseStorageType, GSFStorageType
    use SysMatrixTypeModule, only : diagonalMatrix_p, skylineMatrix_p
    use SysMatrixTypeModule, only : sparseMatrix_p, outOfCore_p
    use FELinearSolver     , only : FEGetSuperelementOrder
    use FELanczosModule    , only : FELanczosSolve
    use scratchArrayModule , only : getIntegerScratchArray, getRealScratchArray
    use scratchArrayModule , only : releaseScratchArrays
    use reportErrorModule  , only : reportError, internalError
    use reportErrorModule  , only : debugFileOnly_p, error_p, warning_p
    use FFaProfilerInterface,only : ffa_startTimer, ffa_stopTimer

    type(SysMatrixType), intent(inout) :: aMat, bMat
    integer            , intent(in)    :: mip(:), meqn1(:)
    integer            , intent(in)    :: neq, nEigVal, nEigVec
    integer            , intent(out)   :: mop(:)
    real(dp)           , intent(inout) :: tolerance(3)
    real(dp)           , intent(out)   :: eigValues(:), eigVectors(:,:)
    real(dp)           , intent(in)    :: shift
    integer            , intent(in)    :: ipsw, lpu
    integer(ik)        , intent(out)   :: meqnInErr(:)
    integer            , intent(out)   :: ierr

    !! Local variables
    logical           :: illegalComb
    integer           :: i, aType, bType, nWork, maxlan, neqAll, neq1, neq2
    integer(ik)       :: lerr, mmip(10), mmop(10)
    integer , pointer :: iWork(:)
    real(dp), pointer :: rWork(:)

    type(SparseStorageType), pointer :: sprA, sprB
    type(GSFStorageType)   , pointer :: gsfA, gsfB

    character(len=39), parameter :: prgnam_p = &
         'SolExtensionModule::csLanczosEigenSolve'

    external LMVY18, LMMY18, CMVY18, CMMY18

    !! --- Logic section ---

    call ffa_startTimer ('csLanczosEigenSolve')

    ierr   = 0
    aType  = aMat%storageType
    bType  = bMat%storageType
    neqAll = aMat%dim
    maxlan = min(size(eigValues),size(eigVectors,2),neq)
    if (bType == aType .or. bType == diagonalMatrix_p) then
       illegalComb = .false.
    else if (aType == outOfCore_p .and. bType == sparseMatrix_p) then
       illegalComb = .false.
    else
       illegalComb = .true.
    end if
    if (illegalComb) then
       ierr = internalError(prgnam_p//': Incompatible system matrices')
    else if (maxlan < nEigVal .or. neqAll < neq .or. bMat%dim /= neqAll) then
       ierr = internalError(prgnam_p//': Invalid input arguments')
    else if (size(eigVectors,1) < neq) then
       ierr = internalError(prgnam_p//': Invalid dimensions on eigVectors')
    else if (size(mip) < 7 .or. size(mop) < 10) then
       ierr = internalError(prgnam_p//': Invalid dimensions on MIP and/or MOP')
    else if (aType == outOfCore_p .and. bType == diagonalMatrix_p) then
       if (neq < neqAll .and. size(meqn1) /= neq) then
          ierr = internalError(prgnam_p//': Invalid dimensions on MEQN1')
       end if
    end if

    meqnInErr = 0_ik
    if (ierr < 0) return

    select case (aType)

    case (skylineMatrix_p) ! Use the SAM skyline solver

       nWork =  2*maxlan
       iWork => getIntegerScratchArray(nWork,ierr)
       if (ierr < 0) goto 999

       nWork =  2*neq + maxlan*(maxlan+8)
       rWork => getRealScratchArray(nWork,ierr)
       if (ierr < 0) goto 999

       if (aMat%skyline%eigensolver == 1) then
          call LANCZ1 (aMat%value(1), bMat%value(1), tolerance(1),       &
               &       aMat%skyline%msky(1), mip(1), mop(1),             &
               &       eigValues(1), rWork(1), eigVectors(1,1),          &
               &       rWork(1+maxlan), rWork(1+maxlan+2*neq), iWork(1), &
               &       shift, neq, nEigVal, nEigVec, maxlan, lpu, ipsw, ierr)
       else
          call LANCZ2 (aMat%value(1), bMat%value(1), tolerance(1),       &
               &       aMat%skyline%msky(1), mip(1), mop(1),             &
               &       eigValues(1), rWork(1), eigVectors(1,1),          &
               &       rWork(1+maxlan), rWork(1+maxlan+2*neq), iWork(1), &
               &       shift, neq, nEigVal, nEigVec, maxlan, lpu, ipsw, ierr)
       end if
       if (ierr < 0) then
          if (mop(1) < 0) meqnInErr(1) = int(-mop(1),ik)
          goto 998
       end if

       !! Expand the eigVectors array from (neq,nEigVec) to (neqAll,nEigVec)
       neqAll = size(eigVectors,1)
       if (neq < neqAll) call DMEXP (neq,nEigVec,neqAll,nEigVec,eigVectors(1,1))

    case (sparseMatrix_p) ! Use the SAM sparse matrix solver

       sprA => aMat%sparse
       if (associated(bMat%sparse)) then
          sprB => bMat%sparse
       else
          sprB => aMat%sparse
       end if
       neq2 = int(sprA%mspar(54))
       neq1 = neqAll - neq2

       if (neq == neqAll) then
          sprA%mspar(57) = 0_ik
          sprB%mspar(57) = 0_ik
       else if (neq == neq1) then
          sprA%mspar(57) = 1_ik
          sprB%mspar(57) = 1_ik
       else
          ierr = internalError(prgnam_p//': Invalid input arguments')
          return
       end if

       nWork =  max(2*maxlan+int(sprA%mspar(13)),neqAll)
       if (ik > i4) nWork = 2*nWork
       iWork => getIntegerScratchArray(nWork,ierr)
       if (ierr < 0) goto 999

       nWork =  max(maxlan*(maxlan+7), int(sprA%mspar(14)+sprA%mspar(17)))
       nWork =  nWork + maxlan + 2*neq
       rWork => getRealScratchArray(nWork,ierr)
       if (ierr < 0) goto 999

       call castInt (size(mip),mip,mmip)
       call SPRLAN (aMat%value(1), bMat%value(1),  tolerance(1), &
            &       sprA%mspar(1), sprA%mtrees(1), sprA%msifa(1), &
            &       sprB%mspar(1), sprB%mtrees(1), sprB%msifa(1), &
            &       mmip(1), mmop(1), meqnInErr(1), &
            &       eigValues(1), rWork(1), eigVectors(1,1), &
            &       rWork(1+maxlan), rWork(1+maxlan+2*neq), iWork(1), shift, &
            &       int(neq,ik), int(nEigVal,ik), int(nEigVec,ik), &
            &       int(maxlan,ik), int(lpu,ik), int(ipsw,ik), lerr, &
            &       LMVY18, LMMY18, CMVY18, CMMY18)
       call castInt (size(mop),mmop,mop)
       ierr = int(lerr)
       if (ierr < 0) goto 998

       sprA%mspar(57) = 0_ik
       sprB%mspar(57) = 0_ik

       !! Expand the eigVectors array from (neq,nEigVec) to (neqAll,nEigVec)
       neqAll = size(eigVectors,1)
       if (neq < neqAll) call DMEXP (neq,nEigVec,neqAll,nEigVec,eigVectors(1,1))

    case (outOfCore_p) ! Use the DnV GSF solver

       gsfA => aMat%gsf
       if (associated(bMat%gsf)) then
          gsfB => bMat%gsf
       else
          gsfB => aMat%gsf
       end if

       nullify(sprB)
       nullify(rWork)
       if (bType == diagonalMatrix_p) then
          call FEGetSuperelementOrder (gsfA%Q, neq1, neq2)
          if (neq == neqAll) then
             rWork => bMat%value
          else if (neq == neq1 .and. neq1+neq2 == neqAll) then
             !! Squeeze out the status 2 equations from the diagonal matrix
             rWork => getRealScratchArray(neq1,ierr)
             if (ierr < 0) goto 999
             do i = 1, neq1
                rWork(i) = bMat%value(meqn1(i))
             end do
          else
             ierr = internalError(prgnam_p//': Invalid input arguments')
             return
          end if
       else if (bType == sparseMatrix_p) then
          call FEGetSuperelementOrder (gsfA%Q, neq1, neq2)
          sprB => bMat%sparse
          if (neq == neqAll) then
             sprB%mspar(57) = 0_ik
          else if (neq == neq1 .and. neq1+neq2 == neqAll) then
             sprB%mspar(57) = 1_ik
          else
             ierr = internalError(prgnam_p//': Invalid input arguments')
             return
          end if
       end if

       if (abs(shift) > epsilon(1.0_dp)) then
          call reportError (warning_p,'Non-zero shift is not supported for '// &
               &            'the GSF equation solver (ignored)')
       end if

       call castInt (size(mip),mip,mmip)
       call ffa_startTimer ('FELanczosSolve')
       if (bType == sparseMatrix_p) then
          !! Matrix A is GSF and Matrix B is SPR
          call FELanczosSolve (gsfA%Msg, gsfA%Q, gsfA%S, gsfA%T, &
               &               gsfB%Msg, gsfB%Q, gsfB%S, gsfB%T, &
               &               sprB%mspar, sprB%mtrees, sprB%msifa, &
               &               valueB=bMat%value, TOL=tolerance, &
               &               EVL=eigValues, EVC=eigVectors, &
               &               MIP=mmip, MOP=mmop, N=neq, &
               &               NEVAL=nEigVal, NEVEC=nEigVec, MAXLAN=maxlan, &
               &               MSING=meqnInErr, LPU=lpu, IPSW=ipsw, &
               &               INFO=GSFinfo)
       else
          !! Matrix A is GSF and Matrix B is either GSF or diagonal
          call FELanczosSolve (gsfA%Msg, gsfA%Q, gsfA%S, gsfA%T, &
               &               gsfB%Msg, gsfB%Q, gsfB%S, gsfB%T, &
               &               diagB=rWork, TOL=tolerance, &
               &               EVL=eigValues, EVC=eigVectors, &
               &               MIP=mmip, MOP=mmop, N=neq, &
               &               NEVAL=nEigVal, NEVEC=nEigVec, MAXLAN=maxlan, &
               &               MSING=meqnInErr, LPU=lpu, IPSW=ipsw, &
               &               INFO=GSFinfo)
       end if
       call ffa_stopTimer ('FELanczosSolve')
       call castInt (size(mop),mmop,mop)
       ierr = checkGSFinfo(GSFinfo,gsfA)
       if (ierr < 0) goto 998

       ierr = count(meqnInErr > 0_ik)
       if (mip(1) == 1)                      gsfA%isFactorized = .true.
       if (bType == outOfCore_p) then
          if (mip(2) == 1 .and. mip(6) == 2) gsfB%isFactorized = .true.
       else if (bType == sparseMatrix_p) then
          sprB%mspar(57) = 0_ik
       end if

    case default

       ierr = internalError(prgnam_p//': Invalid matrix type')

    end select

990 call releaseScratchArrays
    call ffa_stopTimer ('csLanczosEigenSolve')
    return

998 if (ipsw < 0) then
       call reportError (error_p,'Failed to solve eigenvalue problem', &
            &            'Rerun with command-line option -debug=1 to get', &
            &            'lower-level messages from the eigenvalue solver.', &
            &            IERR=ierr)
    else
       call reportError (error_p,'Failed to solve eigenvalue problem',IERR=ierr)
    end if
999 call reportError (debugFileOnly_p,prgnam_p)
    goto 990

  end subroutine csLanczosEigenSolve


  !!============================================================================
  !> @brief Expands and rearranges a solution vector.
  !>
  !> @param[in] samData Assembly management data
  !> @param[in] sveq Solution vector expressed in equation order
  !> @param[out] svdof Solution vector expressed in nodal point DOF order
  !> @param[in] S1 Scaling factor for the free DOFs (default 1.0)
  !> @param[in] S2 Scaling factor for the specified DOFs (default equal to S1)
  !>
  !> @details This subroutine expands and rearranges a system vector expressed
  !> in equation order, into a vector corresponding to nodal point order.
  !>
  !> @callergraph
  !>
  !> @author Karl Erik Thoresen
  !> @date Sep 1999
  !>
  !> @author Knut Morten Okstad
  !> @date Jul 2001

  subroutine csExpand (samData,sveq,svdof,S1,S2)

    use kindModule, only : dp
    use SamModule , only : SamType

    type(SamType)    , intent(in)  :: samData
    real(dp)         , intent(in)  :: sveq(:)
    real(dp)         , intent(out) :: svdof(:)
    real(dp),optional, intent(in)  :: S1, S2

    !! Local variables
    integer  :: neq, ndof
    real(dp) :: F1, F2

    !! --- Logic section ---

    if (present(S1)) then
       F1 = S1
    else
       F1 = 1.0_dp
    end if
    if (present(S2)) then
       F2 = S2
    else
       F2 = F1
    end if

    ndof  = min(size(svdof),samData%ndof)
    neq   = min(size(sveq),samData%neq)
    svdof = 0.0_dp

    call EXPAND (sveq(1)         , samData%ttcc(1), samData%mpmceq(1), &
         &       samData%mmceq(1), samData%meqn(1), F1, F2, ndof, neq, svdof(1))

  end subroutine csExpand

end module SolExtensionModule
