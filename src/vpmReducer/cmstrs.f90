!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file cmstrs.f90
!> @brief Subroutines for Component Mode Synthesis (CMS) transformation.

!!==============================================================================
!> @brief Module with subroutines for doing Component Mode Synthesis.

module CmstrsModule

  implicit none

  private :: reportSingEqns


contains

  !!============================================================================
  !> @brief Solves the eigenvalue problem of the FE model.
  !>
  !> @param[in] sam Data for managing system matrix assembly
  !> @param SSM System stiffness matrix
  !> @param SMM System mass matrix
  !> @param[in] nEigVal Number of eigenvalues to compute
  !> @param[in] factorMass If .true., factorize the mass matrix
  !> @param[in] iopSing Option for handling matrix singularities
  !> @param[in] tolFactorize Zero pivot tolerance
  !> @param[in] tolEigval Zero tolerance on eigenvalues
  !> @param[in] eigShift Eigenvalue shift
  !> @param[out] EVL Computed eigenvalues
  !> @param[out] PHI Associated eigenvectors
  !> @param[in] IPSW Print switch
  !> @param[in] LPU File unit number for res-file output
  !> @param[out] IERR Error flag
  !>
  !> @details The @a nEigVal lowest eigenvalues and associated eigenvectors are
  !> computed using the Lanczos eigenvalue solver. The LAPack solver is used if
  !> the system matrices are dense/full, or if @a tolEigval is non-positive.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Ketil Aamnes
  !> @date Oct 1988
  !>
  !> @author Knut Morten Okstad
  !> @date Mar 2003

  subroutine EIGVAL (sam, SSM, SMM, nEigVal, factorMass, iopSing, &
       &             tolFactorize, tolEigval, eigShift, &
       &             EVL, PHI, IPSW, LPU, IERR)

    use sprKindModule      , only : ik, nbik_p
    use KindModule         , only : dp, nbd_p
    use SamModule          , only : SamType
    use SysMatrixTypeModule, only : SysMatrixType, diagonalMatrix_p
    use SysMatrixTypeModule, only : outOfCore_p, sparseMatrix_p, denseMatrix_p
    use MatExtensionModule , only : csGetSub11
    use SolExtensionModule , only : csLanczosEigenSolve
    use DenseMatrixModule  , only : solveEigenvalues
    use ManipMatrixModule  , only : writeObject
    use AllocationModule   , only : reAllocate, logAllocMem, doLogMem
    use ProgressModule     , only : writeProgress
    use ReportErrorModule  , only : AllocationError, reportError
    use ReportErrorModule  , only : note_p, warning_p, debugFileOnly_p

    type(SamType)      , intent(in)    :: sam
    type(SysMatrixType), intent(inout) :: ssm, smm
    integer            , intent(in)    :: nEigVal, iopSing, IPSW, LPU
    logical            , intent(in)    :: factorMass
    real(dp)           , intent(in)    :: tolFactorize, tolEigval, eigShift
    real(dp)           , pointer       :: EVL(:), PHI(:,:)
    integer            , intent(out)   :: IERR

    !! Local variables
    integer                  :: i, ndof1, maxlan, mip(7), mop(10)
    integer(ik), allocatable :: meqErr(:)
    real(dp)                 :: tolerances(3)
    real(dp), allocatable    :: lambda(:), phiFull(:,:), K11(:,:), M11(:,:)
    character(len=128)       :: errMsg

    !! --- Logic section ---

    ndof1 = sam%ndof1
    if (ndof1 < 1 .or. nEigVal < 1) return

    call writeProgress(' --> EIGVAL   ; Execute eigenvalue solution')

    if (tolEigval <= 0.0_dp .or. ssm%storageType == denseMatrix_p) then

       !! Use the LAPack eigenvalue solver (full matrices) for small systems
       write(lpu,'(/)')
       call reportError (note_p,'Using the LAPack eigenvalue solver')

       if (smm%storageType == diagonalMatrix_p) then
          i = 1
       else
          i = sam%ndof1
       end if
       allocate(K11(ndof1,ndof1),M11(ndof1,i),meqErr(0), &
            &   lambda(nEigVal),phiFull(ndof1,nEigVal),STAT=ierr)
       if (ierr /= 0) then
          ierr = AllocationError('EIGVAL')
          return
       else if (doLogMem) then
          call logAllocMem ('EIGVAL',0,(ndof1+i+nEigVal)*ndof1+nEigVal,nbd_p)
       end if

       call csGetSub11 (ssm,sam%meqn1,K11,ierr)
       if (ierr == 0) then
          call csGetSub11 (smm,sam%meqn1,M11,ierr)
          if (ierr == 0) then
             call solveEigenvalues (K11,M11,lambda,phiFull, &
                       &            ndof1,nEigVal,nEigVal,ierr)
          end if
       end if

       if (doLogMem) then
          call logAllocMem ('EIGVAL',size(K11)+size(M11),0,nbd_p)
       end if
       deallocate(K11,M11)
       if (ierr < 0) goto 900

       goto 100

    end if

    mip(1) =  1   ! Stores aMat on entry
    if (ssm%storageType==outOfCore_p .and. smm%storageType==sparseMatrix_p) then
       mip(2) = 3 ! Stores bMat on entry (in-core sparse)
    else if (smm%storageType == diagonalMatrix_p) then
       mip(2) = 2 ! Stores diagonal of bMat on entry
    else
       mip(2) = 1 ! Stores bMat on entry
    end if
    mip(3) = -1   ! Pseudo random start vector generated by sprlan
    mip(4) =  1   ! Full reorthogonalization (single pass)
    mip(5) =  2   ! Frequency of solving small triangular eigenvalue problem
    if (factorMass) then
       mip(6) = 2 ! Transform problem through factorization of bMat
    else
       mip(6) = 1 ! Transform problem through factorization of aMat
    end if
    if (iopSing > 0) then
       mip(7) = 1+iopSing ! Modify matrix diagonal on singularities and continue
    else
       mip(7) = 1 ! Modify matrix diagonal on singularities, but abort
    end if

    tolerances(1) = tolEigval
    tolerances(2) = tolFactorize
    tolerances(3) = 1.0e-8_dp

    maxlan = min(ndof1,3*nEigVal+12)
    allocate(lambda(maxlan),phiFull(ndof1,maxlan),STAT=ierr)
    if (ierr /= 0) then
       ierr = AllocationError('EIGVAL')
       return
    else if (doLogMem) then
       call logAllocMem ('EIGVAL',0,size(lambda)+size(phiFull),nbd_p)
    end if
    allocate(meqErr(ndof1),STAT=ierr)
    if (ierr /= 0) then
       ierr = AllocationError('EIGVAL')
       return
    else if (doLogMem) then
       call logAllocMem ('EIGVAL',0,size(meqErr),nbik_p)
    end if

    call csLanczosEigenSolve (ssm, smm, mip, mop, sam%meqn1, ndof1, &
         &                    nEigVal, nEigVal, tolerances, eigShift, &
         &                    lambda, phiFull, ipsw-1, lpu, meqErr, ierr)

    if (ierr /= 0) then
       if (ssm%storageType == outOfCore_p) then
          !! The GSF solver is being used
          if (meqErr(1) > 0) then
             call reportSingEqns ('stiffness- and/or mass matrix',sam,meqErr)
          end if
       else if (smm%storageType == sparseMatrix_p) then
          !! Both ssm and smm are sparse
          if (ssm%sparse%mspar(27) /= 0 .and. smm%sparse%mspar(27) /= 0) then
             call reportSingEqns ('stiffness- and mass matrix',sam,meqErr)
          else if (ssm%sparse%mspar(27) /= 0) then
             call reportSingEqns ('stiffness matrix',sam,meqErr)
          else if (smm%sparse%mspar(27) /= 0) then
             call reportSingEqns ('mass matrix',sam,meqErr)
          end if
       else if (ssm%storageType == sparseMatrix_p) then
          !! Only ssm is sparse (smm is diagonal)
          if (ssm%sparse%mspar(27) /= 0) then
             call reportSingEqns ('stiffness- and/or mass matrix',sam,meqErr)
          end if
       else if (ierr == -2) then
          call reportSingEqns ('stiffness matrix',sam,meqErr)
       else if (ierr == -3) then
          call reportSingEqns ('mass matrix',sam,meqErr)
       end if
       if (ierr < 0) then
          if (ndof1 < 1000) then
             call reportError (note_p,'The LANCZOS eigenvalue solver failed'// &
                  &            ' for a system with less than 1000 unknowns.', &
                  &            'Try switching to the LAPack solver instead'// &
                  &            ' (specify negative Eigenvalue tolerance).')
          end if
          goto 900
       else if (meqErr(1) > 0) then
          call reportError (warning_p,'The above singularities have been '// &
               'circumvented by replacing the zero (or nearly zero)', &
               'pivots by some non-zero values (see res-file for details).')
       end if
       ierr = 0
    end if

100 continue
    call reAllocate ('EIGVAL',evl,nEigVal,ierr)
    call reAllocate ('EIGVAL',phi,ndof1,nEigVal,ierr)
    if (ierr < 0) goto 900

    do i = 1, nEigVal
       if (lambda(i) < 0.0_dp) then
          ierr = ierr + 1
          if (ierr == 1) write(lpu,*)
          write(lpu,610) i, lambda(i)
610       format('  ** Warning: Eigenvalue no.',I4,' is negative:',1PE13.5)
          evl(i) = -sqrt(-lambda(i))
       else
          evl(i) = sqrt(lambda(i))
       end if
    end do

    if (ierr > 0) then
       write(errMsg,620) ierr
620    format(I4,' negative eigenvalue(s) associated with', &
            &    ' component mode shapes were detected.')
       call reportError (warning_p,errMsg, &
            &            'Check that the FE-model is consistent.')
    end if

    call DCOPY (ndof1*nEigVal,phiFull(1,1),1,phi(1,1),1)

    if (54.ge.sam%mpar(25).and.54.le.sam%mpar(26)) then
       call writeObject(evl,lpu,'Angular Eigenfrequencies')
    end if
    if (56.ge.sam%mpar(25).and.56.le.sam%mpar(26)) then
       call writeObject(phi,lpu,'Eigenvectors')
    end if

900 continue
    if (ierr < 0) call reportError (debugFileOnly_p,'EIGVAL')
    if (doLogMem) then
       call logAllocMem ('EIGVAL',size(lambda)+size(phiFull),0,nbd_p)
       call logAllocMem ('EIGVAL',size(meqErr),0,nbik_p)
    end if
    deallocate(lambda,phiFull,meqErr)

  end subroutine EIGVAL


  !!============================================================================
  !> @brief Establishes the full gravitational force vectors.
  !>
  !> @param[in] sam Data for managing system matrix assembly
  !> @param smm System mass matrix
  !> @param[out] grav Gravitational force vectors
  !> @param[out] ierr Error flag
  !>
  !> @details The gravitatinal force vectors are calculated from the unit
  !> system mass matrix and unit accelerations:
  !> @verbatim
  !>    |g| = [M] |a|
  !> @endverbatim
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date Feb 2018

  subroutine GRAVfull (sam, smm, grav, ierr)

    use SamModule          , only : SamType, dp
    use SysMatrixTypeModule, only : SysMatrixType
    use MatExtensionModule , only : csPremult
    use ReportErrorModule  , only : reportError, debugFileOnly_p

    type(SamType)      , intent(in)    :: sam
    type(SysMatrixType), intent(inout) :: smm
    real(dp)           , intent(out)   :: grav(:,:)
    integer            , intent(out)   :: ierr

    !! Local variables
    integer :: idir, idof, ieq, inod

    !! --- Logic section ---

    call DCOPY (size(grav),0.0_dp,0,grav(1,1),1)
    do inod = 1, sam%nnod
       idof = sam%madof(inod) - 1
       do idir = 1, 3
          ieq  = sam%meqn(idof+idir)
          if (ieq > 0) grav(ieq,idir) = 1.0_dp
       end do
    end do

    call csPremult (smm,grav,ierr)
    if (ierr < 0) call reportError (debugFileOnly_p,'GRAVfull')

  end subroutine GRAVfull


  !!============================================================================
  !> @brief Computes the gravitation force vector for the superelement DOFs.
  !>
  !> @param[in] sam Data for managing system matrix assembly
  !> @param smm System mass matrix
  !> @param BmatDisk Displacement recovery matrix
  !> @param[in] phi Computed eigenvectors (component mode shapes)
  !> @param[out] gravg Gravitational force vectors
  !> @param[in] ndim Superelement dimension
  !> @param[in] ngen Number of component modes
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @details
  !> First, compute the full gravitation force vector (g) based on unit
  !> acceleration vectors (a) and the mass matrix (M):
  !> @verbatim
  !>       g = M*a   or    |g_i| = | M_ii M_ie ||a_i|
  !>                       |g_e|   | M_ei M_ee ||a_e|
  !> @endverbatim
  !>
  !> Note that the internal and external DOFs can be mixed, i.e.,
  !> not necessarily having the nice 2*2 block structure indicated above.
  !>
  !> Then transform this to the reduced DOF set for the superelement through
  !> the virtual work equation:
  !> @verbatim
  !>       dvs'*gs = dv'*g  where v = H*vs  or |v_i| = | B  phi ||vs_e|
  !>                                           |v_e|   | I   0  ||vs_g|
  !> @endverbatim
  !>
  !> This gives the transformation of the force vector as
  !> @verbatim
  !>       gs = H'*g   or   |gs_e| = |  B'  I ||g_i|
  !>                        |gs_g| = | phi' 0 ||g_e|
  !> @endverbatim
  !>
  !> where (gs) is the superelement force vector, _e indicates external DOFs,
  !> _i indicates internal DOFs, and _g indicates generalized DOFs.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Ketil Aamnes
  !> @date Oct 1988
  !>
  !> @author Bjorn Haugen
  !> @date Nov 1998
  !>
  !> @author Knut Morten Okstad
  !> @date Mar 2003

  subroutine GRAV (sam, smm, BmatDisk, phi, gravg, ndim, ngen, lpu, ierr)

    use KindModule         , only : dp, nbd_p
    use SamModule          , only : SamType
    use SysMatrixTypeModule, only : SysMatrixType
    use DiskMatrixModule   , only : DiskMatrixType, dmGetColPtr
    use ManipMatrixModule  , only : writeObject
    use AllocationModule   , only : logAllocMem, doLogMem
    use ProgressModule     , only : writeProgress
    use ReportErrorModule  , only : AllocationError
    use ReportErrorModule  , only : reportError, debugFileOnly_p

    type(SamType)       , intent(in)    :: sam
    type(SysMatrixType) , intent(inout) :: smm
    type(DiskMatrixType), intent(inout) :: BmatDisk
    integer             , intent(in)    :: ndim, ngen, lpu
    real(dp)            , intent(in)    :: phi(:,:)
    real(dp)            , intent(out)   :: gravg(:,:)
    integer             , intent(out)   :: ierr

    !! Local variables
    integer               :: idir, idof, nnod, ndof1, ndof2, neq
    real(dp), allocatable :: gFull(:,:), gSupEl(:,:), vec_neq(:)
    real(dp), pointer     :: vec_ndof1(:)
    real(dp), external    :: DDOT

    !! --- Logic section ---

    call writeProgress(' --> GRAV     ; Compute gravitation forces')

    call DCOPY (size(gravg),0.0_dp,0,gravg(1,1),1)
    nnod  = sam%nnod
    ndof1 = sam%ndof1
    ndof2 = sam%ndof2
    neq   = sam%neq

    allocate(gFull(neq,3),gSupEl(ndim,3),vec_neq(neq),STAT=ierr)
    if (ierr /= 0) then
       ierr = AllocationError('GRAV')
       return
    else if (doLogMem) then
       call logAllocMem ('GRAV',0,4*neq+3*ndim,nbd_p)
    end if

    !! Establish the full gravitational force vectors
    call GRAVfull (sam,smm,gFull,ierr)
    if (ierr < 0) goto 900

    if (ndof1 > 0) then

       !! Establish the supEl gravitational force vector for the external dofs

       call DCOPY (neq,0.0_dp,0,vec_neq(1),1)
       do idof = 1, ndof2
          vec_ndof1 => dmGetColPtr(BmatDisk,idof,ierr)
          if (ierr < 0) goto 900

          call DSCATR (ndof1,sam%meqn1(1),vec_ndof1(1),vec_neq(1),1)
          vec_neq(sam%meqn2(idof)) = 1.0_dp
          do idir = 1, 3
             gSupEl(idof,idir) = DDOT(neq,vec_neq(1),1,gFull(1,idir),1)
          end do
          vec_neq(sam%meqn2(idof)) = 0.0_dp
       end do

       !! Establish the contribution to the generalized dofs, if any

       do idof = 1, ngen
          call DSCATR (ndof1,sam%meqn1(1),phi(1,idof),vec_neq(1),1)
          do idir = 1, 3
             gSupEl(ndof2+idof,idir) = DDOT(neq,vec_neq(1),1,gFull(1,idir),1)
          end do
       end do

    else

       !! All nodal DOFs are external

       do idir = 1, 3
          call DGATHR (ndof2,sam%meqn2(1),gFull(1,idir),gSupEl(1,idir),1)
       end do

    end if

    !! Write the gravitation vectors in dof ordering

    do idir = 1, 3
       call DGATHR (ndof2,sam%dofPosIn2(1),gSupEl(1,idir),gravg(1,idir),1)
       if (ngen < 1) cycle
       call DCOPY (ngen,gSupEl(ndof2+1,idir),1,gravg(ndof2+1,idir),1)
    end do

    if (88 >= sam%mpar(25) .and. 88 <= sam%mpar(26)) then
       call writeObject(gravg,lpu,'GRAV')
    end if

    call writeProgress('     Done Computing gravitation forces')

900 continue
    if (doLogMem) then
       call logAllocMem ('GRAV',4*neq+3*ndim,0,nbd_p)
    end if
    deallocate(gFull,gSupEl,vec_neq)

    if (ierr < 0) call reportError (debugFileOnly_p,'GRAV')

  end subroutine GRAV


  !!============================================================================
  !> @brief Performs CMS-transformation on the system matrices.
  !>
  !> @param[in] iopSing Option for handling matrix singularities
  !> @param sam Data for managing system matrix assembly
  !> @param bmatDisk Displacement recovery matrix
  !> @param SMII System mass matrix associated with internal DOFs
  !> @param[in] SMIE System mass matrix coupling the internal and external DOFs
  !> @param[out] SMMAT Reduced mass matrix
  !> @param[in] SMEE System mass matrix associated with external DOFs
  !> @param SKII System stiffness matrix associated with internal DOFs
  !> @param[in] SKEE System stiffness matrix associated with external DOFs
  !> @param[out] SKMAT Reduced stiffness matrix
  !> @param[in] PHI Computed eigenvectors (component mode shapes)
  !> @param[in] fullRHS System right-hand-side vectors
  !> @param[out] RHS Reduced right-hand-side vectors
  !> @param[out] VGI Internal displacements due to gravitational forces
  !> @param[in] tolFactorize Zero pivot tolerance (singular system matrix)
  !> @param iStiff File handle for temporary storage of system stiffness matrix
  !> @param[in] NGEN Number of component modes
  !> @param[in] NRHS Number of right-hand-side vectors
  !> @param[in] LPU File unit number for res-file output
  !> @param[out] IERR Error flag
  !>
  !> @details CMS : "Component Mode Synthesis".
  !>
  !> @callgraph @callergraph
  !>
  !> @author Ketil Aamnes
  !> @date Oct 1988
  !>
  !> @author Knut Morten Okstad
  !> @date Mar 2003

  subroutine CMSTRS (iopSing, sam, bmatDisk, &
       &             SMII, SMIE, SMEE, SMMAT, SKII, SKEE, SKMAT, &
       &             PHI, fullRHS, RHS, VGI, tolFactorize, &
       &             iStiff, NGEN, NRHS, LPU, IERR)

    use sprKindModule      , only : ik, nbik_p
    use KindModule         , only : dp, i8, nbd_p
    use SamModule          , only : SamType
    use SysMatrixTypeModule, only : SysMatrixType, diagonalMatrix_p, outOfCore_p
    use DiskMatrixModule   , only : DiskMatrixType
    use SparseMatrixModule , only : SparseMatrixType
    use SysMatrixTypeModule, only : restoreSysMat
    use DiskMatrixModule   , only : dmWrite, dmSetReadOnly, dmGetColPtr
    use SparseMatrixModule , only : smWrite, smSize, smMatTransTimesVec
    use SolExtensionModule , only : csSolveSupEl, csSolve
    use MatExtensionModule , only : csPremult, csTransform
    use ManipMatrixModule  , only : writeObject
    use AllocationModule   , only : logAllocMem, doLogMem
    use ProgressModule     , only : writeProgress
    use ReportErrorModule  , only : AllocationError
    use ReportErrorModule  , only : reportError, warning_p, debugFileOnly_p

    integer               , intent(in)    :: iopSing
    type(SamType)         , intent(inout) :: sam
    type(DiskMatrixType)  , intent(inout) :: bmatDisk
    type(SysMatrixType)   , intent(inout) :: smii, skii
    type(SparseMatrixType), intent(in)    :: smie
    integer               , intent(inout) :: iStiff
    integer               , intent(in)    :: NGEN, NRHS, LPU
    real(dp)              , intent(in)    :: tolFactorize
    real(dp)              , intent(in)    :: SMEE(:,:), SKEE(:,:)
    real(dp), target      , intent(in)    :: PHI(:,:)
    real(dp)              , intent(in)    :: fullRHS(:,:)
    real(dp)              , intent(out)   :: SMMAT(:,:), SKMAT(:,:)
    real(dp)              , intent(out)   :: RHS(:,:), VGI(:,:)
    integer               , intent(out)   :: IERR

    !! Local variables
    integer                  :: i, j, NDOF1, NDOF2, NEQ1
    integer(i8)              :: nStiff
    integer(ik), allocatable :: meqErr(:)
    real(dp)                 :: dummy(1,1)
    real(dp),    allocatable :: SM11(:,:), SM12(:,:), SM22(:,:)
    real(dp),    allocatable :: SK11(:,:), SK22(:,:)
    real(dp),    allocatable :: RHS1(:,:), RHS2(:,:), gFull(:,:)
    real(dp),    allocatable :: vec_ndof2(:), vec2_neq(:)
    real(dp),    pointer     :: vec_neq(:), vec_ndof1(:), phiFull(:,:)
    real(dp),    external    :: DDOT

    !! --- Logic section ---

    ierr  = 0
    NDOF1 = sam%ndof1
    NDOF2 = sam%ndof2
    call writeProgress(' --> CMSTRS   ; Execute CMS-transformation')

    if (NDOF1 < 1) then

       !! Special build of SMMAT and SKMAT if NDOF1=0 (i.e., no internal DOFs)

       do I = 1, NDOF2
          do J = 1, NDOF2
             SKMAT(J,I) = SKEE(sam%dofPosIn2(J),sam%dofPosIn2(I))
          end do
       end do
       call padRedmatrix (SKMAT,NDOF2)

       if (size(smee) > 0) then
          if (smii%storageType == diagonalMatrix_p) then
             call DCOPY (size(smee),0.0_dp,0,SMMAT(1,1),1)
             do I = 1, NDOF2
                SMMAT(I,I) = SMEE(sam%dofPosIn2(I),1)
             end do
          else
             do I = 1, NDOF2
                do J = 1, NDOF2
                   SMMAT(J,I) = SMEE(sam%dofPosIn2(J),sam%dofPosIn2(I))
                end do
             end do
             call padRedmatrix (SMMAT,NDOF2)
          end if
       end if

       !! Special build of RHS if NDOF1=0 (i.e., no internal DOFs)

       do I = 1, NRHS
          call DGATHR (NDOF2,sam%dofPosIn2(1),fullRHS(1,I),RHS(1,I),1)
       end do

       if (sam%mpar(25) < 83 .and. sam%mpar(26) > 80) then
          call writeMatrices (sam%mpar,lpu)
       end if

       return

    end if

    !! === CMS TRANSFORMATION OF THE STIFFNESS MATRIX ===

    !! Compute  K22 = Kee + B'*Kie
    !! and compute the external-to-internal influence matrix BMAT
    !! ie B = -(Kii)^-1*Kie  or  Solve Kii*B = -Kie  with respect to B

    call writeProgress('     Compute K22 = Kee + B`*Kie ,  B = -(Kii)^-1*Kie')

    allocate(meqErr(ndof1),sk22(ndof2,ndof2),STAT=ierr)
    if (ierr /= 0) then
       ierr = AllocationError('CMSTRS: sk22')
       return
    else if (doLogMem) then
       call logAllocMem ('CMSTRS',0,size(meqErr),nbik_p)
       call logAllocMem ('CMSTRS',0,size(sk22),nbd_p)
    end if
    if (nrhs > 0) then
       allocate(rhs2(ndof2,nrhs),STAT=ierr)
       if (ierr /= 0) then
          ierr = AllocationError('CMSTRS: rhs2')
          return
       else if (doLogMem) then
          call logAllocMem ('CMSTRS',0,size(rhs2),nbd_p)
       end if

       call csSolveSupEl (iopSing, skii, fullRHS, sk22, rhs2, &
            &             sam%meqn1, sam%meqn2, nrhs, lpu, ierr, &
            &             bmatDisk, meqnInErr=meqErr, tolFactorize=tolFactorize)
    else
       call csSolveSupEl (iopSing, skii, fullRHS, sk22, dummy, &
            &             sam%meqn1, sam%meqn2, nrhs, lpu, ierr, &
            &             bmatDisk, meqnInErr=meqErr, tolFactorize=tolFactorize)
    endif

    if (ierr /= 0) then
       if (ierr > 0 .or. ierr == -3) then
          call reportSingEqns ('stiffness matrix',sam,meqErr)
       end if
       if (ierr < 0) goto 900
       call reportError (warning_p,'The above singularities have been '// &
            'circumvented by replacing the zero (or nearly zero)', &
            'pivots by some small non-zero values (see res-file for details).')
    end if

    if (doLogMem) then
       call logAllocMem ('CMSTRS',size(meqErr),0,nbik_p)
    end if
    deallocate(meqErr)

    call writeProgress('     Done computing K22 and Bmat')

    !! Add K22 into SKMAT
    do j = 1, ndof2
       do i = 1, ndof2
          skmat(i,j) = sk22(sam%dofPosIn2(i),sam%dofPosIn2(j))
       end do
    end do
    call padRedMatrix (skmat,ndof2)

    !! Set read only protection on the bmatDisk matrix
    !! (this will save time since swap only reads from file)
    call dmSetReadOnly (bmatDisk,ierr)
    if (ierr < 0) goto 900

    if (size(VGI) > 0) then

       if (iStiff >= 0) then
          !! Restore the un-factorized stiffness matrix from temporary file
          nStiff = size(skii%value,kind=i8)
          call restoreSysMat (skii,'stiffness matrix',nStiff,iStiff,ierr)
          if (ierr < 0) goto 900
       end if

       allocate(gFull(sam%neq,3),STAT=ierr)
       if (ierr /= 0) then
          ierr = AllocationError('CMSTRS')
          return
       else if (doLogMem) then
          call logAllocMem ('CMSTRS',0,sam%neq*3,nbd_p)
       end if

       !! Establish the full gravitational force vectors
       call GRAVfull (sam,smii,gFull,ierr)
       if (ierr < 0) goto 900

       !! Solve for the internal DOFs due to unit gravitational forces
       call csSolve (2,iopSing,skii,gFull,lpu,ierr,neq1=ndof1)
       if (ierr < 0) goto 900

       do i = 1, 3
          call DSCATR (ndof2,sam%meqn2(1),0.0_dp,gFull(1,i),0)
       end do

       call csSolve (6,iopSing,skii,gFull,lpu,ierr,neq1=ndof1)
       if (ierr < 0) goto 900

       do i = 1, 3
          call DGATHR (ndof1,sam%meqn1(1),gFull(1,i),VGi(1,i),1)
       end do

       call writeProgress('     Done computing [vgi] = (Kii)^-1 * [gx,gy,gz]')

    end if

    if (smii%storageType == outOfCore_p) then
       neq1 = NDOF1   ! Use only the internal DOFs
    else
       neq1 = sam%neq ! Use full system
    end if

    if (neq1 > NDOF1 .or. (ngen > 0 .and. nrhs > 0)) then
       allocate(vec_neq(neq1),STAT=ierr)
       if (ierr /= 0) then
          ierr = AllocationError('CMSTRS: vec_neq')
          return
       else if (doLogMem) then
          call logAllocMem ('CMSTRS',0,size(vec_neq),nbd_p)
       end if
    else
       nullify(vec_neq)
    end if

    if (ngen > 0) then

       !! Compute K11 = PHI'*Kii*PHI

       if (iStiff >= 0) then
          !! Restore the un-factorized stiffness matrix from temporary file
          nStiff = size(skii%value,kind=i8)
          call restoreSysMat (skii,'stiffness matrix',nStiff,iStiff,ierr)
          if (ierr < 0) goto 900
       end if

       if (neq1 == NDOF1) then
          phiFull => phi ! Using internal DOFs only
       else
          allocate(phiFull(neq1,ngen),STAT=ierr)
          if (ierr /= 0) then
             ierr = AllocationError('CMSTRS: phiFull')
             return
          else if (doLogMem) then
             call logAllocMem ('CMSTRS',0,size(phiFull),nbd_p)
          end if
          call DCOPY (size(phiFull),0.0_dp,0,phiFull(1,1),1)
          do i = 1, ngen
             call DSCATR (NDOF1,sam%meqn1(1),phi(1,i),phiFull(1,i),1)
          end do
       end if

       allocate(sk11(ngen,ngen),STAT=ierr)
       if (ierr /= 0) then
          ierr = AllocationError('CMSTRS: sk11')
          return
       else if (doLogMem) then
          call logAllocMem ('CMSTRS',0,size(sk11),nbd_p)
       end if
       call csTransform (skii,phiFull,sk11,ierr)
       if (ierr < 0) goto 900

       !! Add K11 into SKMAT
       do i = 1, ngen
          call DCOPY (ngen,sk11(1,i),1,skmat(ndof2+1,ndof2+i),1)
       end do

       if (nrhs > 0) then

          !! Compute RHS1 = PHI'*RHSi

          allocate(rhs1(ngen,nrhs),STAT=ierr)
          if (ierr /= 0) then
             ierr = AllocationError('CMSTRS: rhs1')
             return
          else if (doLogMem) then
             call logAllocMem ('CMSTRS',0,size(rhs1),nbd_p)
          end if

          if (neq1 == NDOF1) then
             !! phiFull is of dimension NDOF1*NGEN, whereas fullRHS is NEQ*NRHS
             !! Must extract the NDOF1 internal DOFs from fullRHS
             do i = 1, nrhs
                call DGATHR (ndof1,sam%meqn1(1),fullRHS(1,i),vec_neq(1),1)
                call DGEMV ('T',ndof1,ngen,1.0_dp,phiFull(1,1),ndof1, &
                     &      vec_neq(1),1,0.0_dp,rhs1(1,i),1)
             end do
          else
             !! phiFull is of dimension NEQ*NGEN, but with external DOFs zero
             call DGEMM ('T','N',ngen,nrhs,neq1,1.0_dp,phiFull(1,1),neq1, &
                  &      fullRHS(1,1),neq1,0.0_dp,rhs1(1,1),ngen)
          end if

       end if

    end if

    if (size(smee) < 1) goto 800 ! Mass matrix reduction not needed


    !! === CMS-TRANSFORMATION OF THE MASS MATRIX ===

    if (smSize(smie,3) > 0) then
       allocate(vec_ndof2(ndof2),STAT=ierr)
       if (ierr /= 0) then
          ierr = AllocationError('CMSTRS: vec_ndof2')
          return
       else if (doLogMem) then
          call logAllocMem ('CMSTRS',0,size(vec_ndof2),nbd_p)
       end if
    end if

    if (neq1 > NDOF1) then
       call DCOPY (neq1,0.0_dp,0,vec_neq(1),1)
    end if

    if (ngen > 0) then

       allocate(sm11(ngen,ngen),sm12(ngen,ndof2),STAT=ierr)
       if (ierr /= 0) then
          ierr = AllocationError('CMSTRS: sm11 and sm12')
          return
       else if (doLogMem) then
          call logAllocMem ('CMSTRS',0,size(sm11)+size(sm12),nbd_p)
       end if

       !! Compute M11 = Phi'*Mii*Phi
       call csTransform (smii,phiFull,sm11,ierr)
       if (ierr < 0) goto 900

       !! Compute Mii*Phi (stored in phiFull)
       call csPremult (smii,phiFull,ierr)
       if (ierr < 0) goto 900

       !! Compute M21 = M12' = B'*Mii*Phi
       do i = ndof2, 1, -1 ! Loop backwards to minimize read operations
          vec_ndof1 => getBcolumn(i,ierr)
          if (ierr < 0) goto 900
          do j = 1, ngen
             sm12(j,i) = DDOT(neq1,vec_ndof1(1),1,phiFull(1,j),1)
          end do
       end do

       if (neq1 > NDOF1) then
          if (doLogMem) then
             call logAllocMem ('CMSTRS',size(phiFull),0,nbd_p)
          end if
          deallocate(phiFull)
       end if

       if (smSize(smie,3) > 0) then
          !! Compute M12 = M12 + Phi'*Mie
          do i = 1, ngen
             call smMatTransTimesVec (smie,phi(:,i),vec_ndof2,ierr)
             if (ierr < 0) goto 900
             call DAXPY (ndof2,1.0_dp,vec_ndof2(1),1,sm12(i,1),ngen)
          end do
       end if

    end if

    call writeProgress('     Compute M22 = B`*Mii*B + Mei*B + B`*Mie + Mee')

    allocate(sm22(ndof2,ndof2),vec2_neq(neq1),STAT=ierr)
    if (ierr /= 0) then
       ierr = AllocationError('CMSTRS: sm22 and vec2_neq')
       return
    else if (doLogMem) then
       call logAllocMem ('CMSTRS',0,size(sm22)+size(vec2_neq),nbd_p)
    end if

    !! Compute M22 = B'*Mii*B
    do j = 1, ndof2
       vec_ndof1 => getBcolumn(j,ierr)
       if (ierr < 0) goto 900

       call csPremult (smii,vec_ndof1,vec2_neq,1,ierr)
       if (ierr < 0) goto 900

       sm22(j,j) = DDOT(neq1,vec_ndof1(1),1,vec2_neq(1),1)

       do i = j-1, 1, -1 ! Loop backwards to minimize read operations

          vec_ndof1 => getBcolumn(i,ierr)
          if (ierr < 0) goto 900

          sm22(i,j) = DDOT(neq1,vec_ndof1(1),1,vec2_neq(1),1)
          sm22(j,i) = sm22(i,j)

       end do

       call writeProgress(j,ndof2)
    end do

    call writeProgress('     Done computing  B`*Mii*B')

    if (smii%storageType == diagonalMatrix_p) then

       !! Add M22 = M22 + Mee
       do i = 1, ndof2
          sm22(i,i) = sm22(i,i) + smee(i,1)
       end do

    else if (smSize(smie,3) > 0) then

       !! Add M22 = M22 + Mei*B + B'*Mie
       do i = 1, ndof2
          vec_ndof1 => dmGetColPtr(bmatDisk,i,ierr)
          if (ierr < 0) goto 900

          call smMatTransTimesVec (smie,vec_ndof1,vec_ndof2,ierr)
          if (ierr < 0) goto 900

          !! M22 += Mei*B (= Mie'*B )
          call DAXPY (ndof2,1.0_dp,vec_ndof2(1),1,sm22(1,i),1)
          !! M22 += B'*Mie (= (Mie'*B)')
          call DAXPY (ndof2,1.0_dp,vec_ndof2(1),1,sm22(i,1),ndof2)
       end do

       !! Add M22 = M22 + Mee
       call DAXPY (ndof2*ndof2,1.0_dp,smee(1,1),1,sm22(1,1),1)

    else

       !! Add M22 = M22 + Mee
       call DAXPY (ndof2*ndof2,1.0_dp,smee(1,1),1,sm22(1,1),1)

    end if

    call writeProgress('     Done computing  M22')

    !! Add M22 into SMMAT
    do j = 1, ndof2
       do i = 1, ndof2
          smmat(i,j) = sm22(sam%dofPosIn2(i),sam%dofPosIn2(j))
       end do
    end do

    if (ngen > 0) then

       !! Add M12 into SMMAT
       j = ndof2+1
       do i = 1, ndof2
          call DCOPY (ngen,sm12(1,sam%dofPosIn2(i)),1,smmat(j,i),1)
          call DCOPY (ngen,sm12(1,sam%dofPosIn2(i)),1,smmat(i,j),size(smmat,1))
       end do

       !! Add M11 into SMMAT
       do i = 1, ngen
          call DCOPY (ngen,sm11(1,i),1,smmat(j,ndof2+i),1)
       end do

    end if
    call padRedMatrix (smmat,ndof2+ngen)

800 continue
    if (allocated(rhs2)) then

       !! Add RHS2 into RHS
       do i = 1, nrhs
          call DGATHR (ndof2,sam%dofPosIn2(1),rhs2(1,i),rhs(1,i),1)
       end do

    end if
    if (allocated(rhs1)) then

       !! Add RHS1 into RHS
       do i = 1, nrhs
          call DCOPY (ngen,rhs1(1,i),1,rhs(ndof2+1,1),1)
       end do

    end if

900 continue
    call writeMatrices (sam%mpar,lpu)
    if (allocated(vec_ndof2)) then
       if (doLogMem) call logAllocMem ('CMSTRS',size(vec_ndof2),0,nbd_p)
       deallocate(vec_ndof2)
    end if
    if (associated(vec_neq)) then
       if (doLogMem) call logAllocMem ('CMSTRS',size(vec_neq),0,nbd_p)
       deallocate(vec_neq)
    end if
    if (allocated(vec2_neq)) then
       if (doLogMem) call logAllocMem ('CMSTRS',size(vec2_neq),0,nbd_p)
       deallocate(vec2_neq)
    end if
    if (allocated(sk11)) then
       if (doLogMem) call logAllocMem ('CMSTRS',size(sk11),0,nbd_p)
       deallocate(sk11)
    end if
    if (allocated(sm11)) then
       if (doLogMem) call logAllocMem ('CMSTRS',size(sm11),0,nbd_p)
       deallocate(sm11)
    end if
    if (allocated(sm12)) then
       if (doLogMem) call logAllocMem ('CMSTRS',size(sm12),0,nbd_p)
       deallocate(sm12)
    end if
    if (allocated(sm22)) then
       if (doLogMem) call logAllocMem ('CMSTRS',size(sm22),0,nbd_p)
       deallocate(sm22)
    end if
    if (allocated(sk22)) then
       if (doLogMem) call logAllocMem ('CMSTRS',size(sk22),0,nbd_p)
       deallocate(sk22)
    end if
    if (allocated(rhs1)) then
       if (doLogMem) call logAllocMem ('CMSTRS',size(rhs1),0,nbd_p)
       deallocate(rhs1)
    end if
    if (allocated(rhs2)) then
       if (doLogMem) call logAllocMem ('CMSTRS',size(rhs2),0,nbd_p)
       deallocate(rhs2)
    end if

    if (ierr < 0) call reportError (debugFileOnly_p,'CMSTRS')

  contains

    !> @brief Pads a reduced matrix with zeroes if larger than @a ndim.
    subroutine padRedMatrix (A,ndim)
      real(dp), intent(inout) :: A(:,:)
      integer , intent(in)    :: ndim
      integer :: j, n
      n = size(A,1) - NDIM
      if (n > 0) then
         do j = 1, NDIM
            call DCOPY (n,0.0_dp,0,A(NDIM+1,j),1)
         end do
      end if
      n = size(A,2) - NDIM
      if (n > 0) then
         call DCOPY (size(A,1)*n,0.0_dp,0,A(1,NDIM+1),1)
      end if
    end subroutine padRedMatrix

    !> @brief Extracts a column from the B-matrix expanded to length @a neq1.
    function getBcolumn (icol,ierr) result(vec_ndof1)
      integer , intent(in)  :: icol
      integer , intent(out) :: ierr
      real(dp), pointer     :: vec_ndof1(:)
      vec_ndof1 => dmGetColPtr(bmatDisk,icol,ierr)
      if (ierr == 0 .and. neq1 > NDOF1) then
         call DSCATR (NDOF1,sam%meqn1(1),vec_ndof1(1),vec_neq(1),1)
         vec_ndof1 => vec_neq
      end if
    end function getBcolumn

    !> @brief Debug print of some matrices used in the CMS-reduction.
    subroutine writeMatrices (mpar,lpu)
      integer, intent(in) :: mpar(:), lpu
      if (MPAR(25) > 99 .or. MPAR(26) <= 0) return
      if (size(smee) > 0) then
         if (51>=MPAR(25).and.51<=MPAR(26)) call writeObject(SMEE,LPU,'SMee')
      end if
      if (52>=MPAR(25).and.52<=MPAR(26)) call writeObject(SKEE,LPU,'SKee')
      if (size(smee) > 0) then
         if (53>=MPAR(25).and.53<=MPAR(26)) call smWrite(SMIE,'SMie',LPU,.true.)
      end if
      if (allocated(sk11)) then
         if (79>=mpar(25).and.79<=mpar(26)) call writeObject(sk11,lpu,'SK11')
      end if
      if (allocated(sm11)) then
         if (76>=mpar(25).and.76<=mpar(26)) call writeObject(sm11,lpu,'SM11')
      end if
      if (allocated(sm12)) then
         if (77>=mpar(25).and.77<=mpar(26)) call writeObject(sm12,lpu,'SM12')
      end if
      if (allocated(sm22)) then
         if (78>=mpar(25).and.78<=mpar(26)) call writeObject(sm22,lpu,'SM22')
      end if
      if (allocated(sk22)) then
         if (80>=mpar(25).and.80<=mpar(26)) call writeObject(sk22,lpu,'SK22')
      end if
      if (size(smee) > 0) then
         if (81>=MPAR(25).and.81<=MPAR(26)) call writeObject(SMMAT,LPU,'SMMAT')
      end if
      if (82>=MPAR(25).and.82<=MPAR(26)) call writeObject(SKMAT,LPU,'SKMAT')
      if (NRHS > 0) then
         if (83>=MPAR(25).and.83<=MPAR(26)) call writeObject(RHS,LPU,'RHS')
      end if
      if (size(VGI) > 0) then
         if (84>=MPAR(25).and.84<=MPAR(26)) call writeObject(VGI,LPU,'VGI')
      end if
      if (85>=MPAR(25).and.85<=MPAR(26)) call dmWrite(bmatDisk,LPU,'BMAT')
    end subroutine writeMatrices

  end subroutine CMSTRS


  !!============================================================================
  !> @brief Computes and prints the inertia properties of the superelement.
  !>
  !> @param[in] M Reduced dense superelement mass matrix
  !> @param[in] UNDPOS Position matrices for the supernodes
  !> @param[in] MEDOF Matrix of External DOF indices for the supernodes
  !> @param[in] NDIM  Dimension on the superelement matrices
  !> @param[in] NENOD Number of supernodes
  !> @param[out] MASS Total mass of the superelement
  !> @param[in]  LPU File unit number for res-file output
  !> @param[out] IERR Error flag
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date July 2002

  subroutine JCMS (M,UNDPOS,MEDOF,NDIM,NENOD,MASS,LPU,IERR)

    use KindModule       , only : dp, epsDiv0_p, nbd_p
    use AllocationModule , only : logAllocMem, doLogMem
    use ReportErrorModule, only : AllocationError

    integer , intent(in)  :: NDIM, NENOD, MEDOF(:), LPU
    real(dp), intent(in)  :: M(:,:), UNDPOS(:,:)
    real(dp), intent(out) :: MASS
    integer , intent(out) :: IERR

    !! Local variables
    integer               :: i, j, n
    real(dp)              :: CR(3), Io(6,6), Ic(6,6)
    real(dp), allocatable :: R(:,:)

    !! --- Logic section ---

    allocate(R(ndim,6),STAT=ierr)
    if (ierr /= 0) then
       ierr = AllocationError('JCMS')
       return
    else if (doLogMem) then
       call logAllocMem ('JCMS',0,size(R),nbd_p)
    end if
    call DCOPY (ndim*6,0.0_dp,0,R(1,1),1)

    !! Calculate Inertia matrix about origo
    do i = 1, nenod
       j = medof(i) - 1
       n = medof(i+1) - medof(i)

       R(j+1,1) = 1.0_dp
       R(j+2,2) = 1.0_dp
       R(j+3,3) = 1.0_dp

       R(j+2,4) = -undpos(3,i)
       R(j+3,4) =  undpos(2,i)
       if (n >= 4) R(j+4,4) = 1.0_dp

       R(j+1,5) =  undpos(3,i)
       R(j+3,5) = -undpos(1,i)
       if (n >= 5) R(j+5,5) = 1.0_dp

       R(j+1,6) = -undpos(2,i)
       R(j+2,6) =  undpos(1,i)
       if (n >= 6) R(j+6,6) = 1.0_dp
    end do
    Io = matmul(transpose(R),matmul(M,R))

    !! Calculate the Mass
    Mass = (Io(1,1)+Io(2,2)+Io(3,3)) / 3.0_dp

    !! Calculate the mass centre
    if (mass > epsDiv0_p) then
       CR(1) = Io(2,6)/mass
       CR(2) = Io(3,4)/mass
       CR(3) = Io(1,5)/mass
    else
       CR = 0.0_dp
    end if

    !! Calculate Inertia matrix about the mass centre
    do i = 1, nenod
       j = medof(i) - 1
       R(j+2,4) = -(undpos(3,i) - CR(3))
       R(j+3,4) =  (undpos(2,i) - CR(2))
       R(j+1,5) =  (undpos(3,i) - CR(3))
       R(j+3,5) = -(undpos(1,i) - CR(1))
       R(j+1,6) = -(undpos(2,i) - CR(2))
       R(j+2,6) =  (undpos(1,i) - CR(1))
    end do
    Ic = matmul(transpose(R),matmul(M,R))

    if (doLogMem) then
       call logAllocMem ('JCMS',size(R),0,nbd_p)
    end if
    deallocate(R)

    write(LPU,600) Mass, CR, &
         &         Io(4,4), Ic(4,4), &
         &         Io(5,5), Ic(5,5), &
         &         Io(6,6), Ic(6,6), &
         &         Io(4,5), Ic(4,5), &
         &         Io(5,6), Ic(5,6), &
         &         Io(6,4), Ic(6,4)

600 format(///4X,'MASS AND INERTIA PROPERTIES OF THE REDUCED SYSTEM' &
         &   /4X,'-------------------------------------------------' &
         &  //4X,'MASS      = ',1PE14.5, &
         &   /4X,'CG X_COOR = ',1PE14.5, &
         &   /4X,'CG Y_COOR = ',1PE14.5, &
         &   /4X,'CG Z_COOR = ',1PE14.5, &
         &  //4X,'               About origin      About CG' &
         &   /4X,'IXX       = ',1PE14.5, E18.5, &
         &   /4X,'IYY       = ',1PE14.5, E18.5, &
         &   /4X,'IZZ       = ',1PE14.5, E18.5, &
         &   /4X,'IXY       = ',1PE14.5, E18.5, &
         &   /4X,'IYZ       = ',1PE14.5, E18.5, &
         &   /4X,'IZX       = ',1PE14.5, E18.5 )

  end subroutine JCMS


  !!============================================================================
  !> @brief Computes and prints the eigenvalues of the reduced system.
  !>
  !> @param sk Dense superelement stiffness matrix
  !> @param sm Dense superelement mass matrix
  !> @param Bmat Matrix relating the internal and external nodal DOFs
  !> @param phi Component mode shapes
  !> @param[out] lambda Computed eigenvalues
  !> @param[out] Zequ Computed eigenvectors
  !> @param[in]  meqn1 Matrix of status 1 (internal) equation numbers
  !> @param[in]  meqn2 Matrix of status 2 (external) equation numbers
  !> @param[in]  n Dimension of the eigenvalue problem
  !> @param[in]  ndof1 Number of internal nodal DOFs (status 1 equations)
  !> @param[in]  ndof2 Number of external nodal DOFs (status 2 equations)
  !> @param[in]  nVal Number of eigenvalues to solve for
  !> @param[in]  ipsw Print switch for debug output
  !> @param[in]  lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @details The associated eigenvectors are optionally computed
  !> and expanded to all internal nodes of the superelement.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date Oct 2002

  subroutine EIGCMS (sk,sm,Bmat,phi,lambda,Zequ,meqn1,meqn2, &
       &             n,ndof1,ndof2,nVal,ipsw,lpu,ierr)

    use KindModule       , only : dp, pi_p
    use DiskMatrixModule , only : DiskMatrixType, dmMatTimesVec
    use DenseMatrixModule, only : solveEigenvalues
    use ProgressModule   , only : writeProgress
    use ManipMatrixModule, only : writeObject
    use ReportErrorModule, only : AllocationError, reportError
    use ReportErrorModule, only : warning_p, debugFileOnly_p

    integer             , intent(in)    :: n, ndof1, ndof2, nVal, ipsw, lpu
    integer             , intent(in)    :: meqn1(:), meqn2(:)
    real(dp)            , intent(inout) :: sk(:,:), sm(:,:)
    type(DiskMatrixType), intent(inout) :: Bmat
    real(dp)            , intent(in)    :: phi(:,:)
    real(dp)            , intent(out)   :: lambda(:), Zequ(:,:)
    integer             , intent(out)   :: ierr

    !! Local variables
    integer               :: i, nVec
    real(dp), allocatable :: Z(:,:), Work(:)

    !! --- Logic section ---

    call writeProgress(' --> EIGCMS   ; Compute eigenvalues of reduced system')

    if (size(Zequ) > 0) then
       nVec = nVal
       allocate(Z(N,nVec),Work(ndof1),STAT=ierr)
    else
       nVec = 0
       allocate(Z(1,1),Work(0),STAT=ierr)
    end if
    if (ierr /= 0) then
       ierr = AllocationError('EIGCMS: Eigenvectors')
       return
    end if

    call solveEigenvalues (SK,SM,lambda,Z,n,nVal,nVec,ierr)
    if (ierr > 0) then
       call reportError (warning_p, &
            'A dynamics simulation using this superelement may fail.')
    end if
    if (ierr /= 0) goto 100

    write(lpu,600)
    do i = 1, nVal
       if (lambda(i) < 0.0_dp) then
          lambda(i) = sqrt(-lambda(i))
          write(lpu,610) i,-lambda(i)*0.5_dp/pi_p
       else
          lambda(i) = sqrt(lambda(i))
          write(lpu,611) i,lambda(i)*0.5_dp/pi_p
       end if
    end do

    if (nVec < 1) goto 100

    if (ipsw > 0) then
       write(lpu,*)
       call writeObject(Z(:,1:nVec),lpu,'    CORRESPONDING EIGENVECTORS')
    end if

    call DCOPY (size(Zequ),0.0_dp,0,Zequ(1,1),1)
    do i = 1, min(nVal,size(Zequ,2))

       !! Calculate eigenvector components for the internal degrees of freedom
       call dmMatTimesVec (Bmat,Z(1:ndof2,i),Work,ierr)
       if (ierr < 0) goto 100

       !! Add contributions from the generalized modes to vi (stored in Work)
       if (N > ndof2) then
          !! Work = Work + matmul(phi,Z(ndof2+1:,i))
          !! Using BLAS routine instead of matmul to avoid stack overflow
          call DGEMV ('N',ndof1,N-ndof2,1.0_dp,phi(1,1),ndof1, &
               &      Z(ndof2+1,i),1,1.0_dp,Work(1),1)
       end if

       !! Add vi(ndof1) and Z(ndof2) into Zequ(neq) "equation order"
       if (ndof1 > 0) call DSCATR (ndof1,meqn1(1),Work(1),Zequ(1,i),1)
       if (ndof2 > 0) call DSCATR (ndof2,meqn2(1), Z(1,i),Zequ(1,i),1)

    end do

100 continue
    deallocate(Z,work)
    if (ierr /= 0) call reportError (debugFileOnly_p,'EIGCMS')

600 format(//4X,'EIGEN FREQUENCIES OF THE REDUCED SYSTEM', &
         &  /4X,'---------------------------------------', &
         & //4X,'MODE  FREQUENCY [Hz]')
610 format(I7,1PE17.8,'*i')
611 format(I7,1PE17.8)

  end subroutine EIGCMS


  !!============================================================================
  !> @brief Prints an error message on singular system system.
  !>
  !> @param[in] mattype Type of system matrix (stiffness or mass)
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[in] meqErr List of singular equation numbers
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad

  subroutine reportSingEqns (mattype,sam,meqErr)

    use sprKindModule    , only : ik
    use SamModule        , only : SamType, NodeNumFromEqNum
    use ReportErrorModule, only : reportError, error_p, empty_p

    character(len=*), intent(in) :: mattype
    type(SamType)   , intent(in) :: sam
    integer(ik)     , intent(in) :: meqErr(:)

    !! Local variables
    integer           :: i, j, k
    character(len=64) :: nodMsg

    !! --- Logic section ---

    call reportError (error_p,'Singular or ill-conditioned '//mattype// &
         &            ' detected at the following node(s):')
    do i = 1, size(meqErr)
       j = int(meqErr(i))
       if (j <= 0) return

       k = NodeNumFromEqNum(sam,j,.true.)
       write(nodMsg,"('Node',i9,' local DOF',i3,' ( eqn. # ',i9,' )')") &
            &       k, j, meqErr(i)
       call reportError (empty_p,nodMsg)
    end do

  end subroutine reportSingEqns

end module CmstrsModule
