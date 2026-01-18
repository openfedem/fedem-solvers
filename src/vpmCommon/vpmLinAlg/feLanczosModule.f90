!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module FELanczosModule

  use kindModule, only : wp => dp

  implicit none

  private

  real(wp), allocatable, save :: work1(:,:), work2(:,:)

  public :: FELanczosSolve


contains

  subroutine FELanczosSolve (MsgA,structA,coeffA,factorA, &
       &                     MsgB,structB,coeffB,factorB, &
       &                     MPARB,MTREEB,MSIFB,valueB,diagB, &
       &                     TOL,EVL,EVC,MIP,MOP, &
       &                     N,NEVAL,NEVEC,MAXLAN,MSING,LPU,IPSW,INFO)

    !!==========================================================================
    !!    TASK :  To determine the NEVAL smallest eigenvalues Lambda_i
    !!            of the generalized, symmetric eigenproblem
    !!                 (A - Lambda_i*B)*q_i = 0
    !!            The routine also returns NLV (=MOP(3)) ortogonal Lanczos
    !!            vectors (in V) or, optionally, approximations to NEVEC
    !!            eigenvectors (corresponding to the NEVEC first eigenvalues).
    !!            values).  The eigenvectors, which are B-orthonormal,
    !!            are returned as the first NEVEC vectors in V.
    !!    A is a symmetric matrix stored in sparse (compressed) format.
    !!    B is either a symmetric sparse (KSB=1) or a diagonal (KSB=2) matrix.
    !!    For both A and B the original matrices or their appropriate factors
    !!    are accepted as input (in coeffA/factorA and coeffB/factorB).
    !!    A truncated Lanczos method is used to transform the problem to
    !!    tridiagonal, special form. The reduced, tridiagonal eigenproblem
    !!    is solved by implicit QL transformation.
    !!    Various reorthogonalization schemes may be specified.
    !!
    !!    This is the GSF version of SPRLAN (out-of-core matrices).
    !!    The matrix B may optionally be represented by an SPR data structure,
    !!    when KSB=MIP(2)=3.
    !!
    !! Programmer : Knut Morten Okstad                   date/rev : Jan 2005/1.0
    !!==========================================================================

    use SprKindModule , only : ik
    use FELinearSolver, only : MessageLevel, SAM, CAM, GSFColSup, Error_Flag
    use FELinearSolver, only : FEMatrixOrder, FEExtractDiagA, FESolve
    use FELinearSolver, only : FEPermuteVectors, FEPermuteIndices
    use FELinearSolver, only : FESetError, SetLinearSolverType
    use ProgressModule, only : writeProgress

    use FFaProfilerInterface, only : ffa_startTimer, ffa_stopTImer

    type(MessageLevel), pointer       :: MsgA, MsgB
    type(SAM)         , pointer       :: structA, structB
    type(CAM)         , pointer       :: coeffA, coeffB
    type(GSFColSup)   , pointer       :: factorA, factorB
    integer(ik), optional, intent(in) :: MPARB(:), MTREEB(:), MSIFB(:)
    real(wp), optional, intent(inout) :: valueB(:), diagB(:)
    real(wp)          , intent(inout) :: TOL(3)
    real(wp)          , intent(out)   :: EVL(:), EVC(:,:)
    integer(ik)       , intent(in)    :: MIP(:)
    integer           , intent(in)    :: N,NEVAL,NEVEC,MAXLAN,LPU,IPSW
    integer(ik)       , intent(out)   :: MOP(:), MSING(:)
    type(Error_Flag)  , pointer       :: INFO

    !! Local variables
    character(len=1)         :: AllInt
    logical                  :: solveAllInOneGo, reStart
    integer                  :: I,IERR,IFLAG,IPA,IPB,IPD,IPK,IPKM1,IPKP1,IPS, &
         &                      JORT,K,KEVEX,KM1,LLP,NEVEXT,NLV,NSDIG,NVC
    integer(ik)              :: J,KORT,KEX,KPD,KSA,KSB,KSR,LERR,LLLP, &
         &                      NACEV,NEQA,NEQB,NGDEV
    real(wp)                 :: ALFA,AUX,BETA,RAN,RLPR,SIGMA,SMALL,TRSH,VPMAX
    integer    , external    :: IMP
    real(wp)   , external    :: DDOT
    integer(ik), allocatable :: LWA(:,:), IPDIAG(:)
    real(wp)   , allocatable :: RWA(:,:), RWE(:,:), EVERR(:), R(:)

    !! --- Logic section ---

    IF (IPSW >= 0) THEN
       LLP = LPU
    ELSE
       LLP = -1
    END IF
    LLLP = int(LLP,ik)

#if FT_DEBUG > 1
    IF (LLP > 0) WRITE(LLP,*) 'ENTERING FELanczosSolve'
#endif

    IF (IPSW > 0 .AND. NEVAL > 0) THEN
       CALL SPRLP1 (0_ik,MIP(1),0.0_wp,int(N,ik), &
            &       int(NEVAL,ik),int(NEVEC,ik),int(MAXLAN,ik),int(LPU,ik))
    END IF

    !! --- Check input parameters ----------------------------------------------

    IERR = 0
    KSA  = MIP(1)
    KSB  = MIP(2)
    KSR  = ABS(MIP(3))
    KORT = MIP(4)
    KEX  = MIP(5)
    KPD  = MIP(6)
    NEQA = int(FEMatrixOrder(structA),ik)
    NEQB = int(FEMatrixOrder(structB),ik)
    CALL SPRLCK (KSA,KSB,MIP(3),KORT,KEX,KPD,int(N,ik),NEQA,NEQB, &
         &       int(NEVAL,ik),int(NEVEC,ik),int(MAXLAN,ik),LLLP,LERR)
    IF (LERR < 0_ik) GOTO 989

    !! --- Set/initiate parameters, variables and pointers ---------------------

    NEVEXT = 0
    NACEV  = 0_ik
    NGDEV  = 0_ik
    NLV    = 0
    NSDIG  = IMP(3)
    RLPR   = 10.0_wp**REAL(-NSDIG,wp)
    TRSH   = SQRT(RLPR)
    SMALL  = 10.0_wp**REAL(-NSDIG+NSDIG/5,wp)
    RAN    = 0.1234567891234567_wp
    MOP    = 0_ik
    J      = 0_ik

    if (N == FEMatrixOrder(structA)) then
       AllInt = 'A'
    else
       AllInt = 'I'
    end if

    !! Default tolerance values
    IF (TOL(1) <= 0.0_wp) TOL(1) = 10.0_wp**REAL(-NSDIG*2/3)
    IF (TOL(2) <= 0.0_wp) TOL(2) = TOL(1)
    IF (TOL(3) <  1.0_wp) TOL(3) = 10.0_wp

    allocate(R(N),work1(N,NEVEC),work2(N,1),stat=IERR)
    IF (IERR /= 0) GOTO 980

    IF (KSB == 2_ik) THEN
       !! Reorder the diagonal mass matrix to internal order
       call FEPermuteVectors (AllInt,'T',N,structA,diagB,INFO)
       IF (associated(INFO)) GOTO 990
    END IF

    !! Start vector
    IF (KSR == 1_ik) THEN
       call FERandomVec (R,RAN)
    ELSE IF (KSR == 3_ik) THEN
       call DCOPY (N,1.0_wp,0,R(1),1)
    ELSE IF (KSR == 2_ik) THEN
       IF (KSB == 1_ik) THEN
          !! Extract the diagonal from the GSF-matrix B
          call FEExtractDiagA (coeffB,R,INFO)
          IF (associated(INFO)) THEN
             CALL SPRLER (20_ik,J,J,J,LLLP,LERR)
             GOTO 989
          END IF
       ELSE IF (KSB == 3_ik) THEN
          !! Extract the diagonal from the SPR-matrix B
          allocate(IPDIAG(N),stat=IERR)
          IF (IERR /= 0) GOTO 980
          CALL SPRDG1 (MPARB(1),MTREEB(1),MSIFB(1), &
               &       IPDIAG(1),LLLP,LERR)
          IF (LERR .LT. 0_ik) THEN
             CALL SPRLER (20_ik,J,J,J,LLLP,LERR)
             GOTO 989
          END IF
          DO I = 1, N
             R(I) = valueB(IPDIAG(I))
          END DO
          deallocate(IPDIAG)
       ELSE
          !! Matrix B is diagonal
          call DCOPY (N,diagB(1),1,R(1),1)
       END IF
    END IF

    !! Factorize the coefficient matrices
    call writeProgress ('     Factoring coefficient matrices')
    call FELFac (INFO,LPU,KSA,KSB,int(KPD), &
         &       N,NEVAL,TOL(2),MOP(6),MSING, &
         &       MsgA,structA,coeffA,factorA, &
         &       MsgB,structB,coeffB,factorB,diagB, &
         &       MPARB,MTREEB,MSIFB,valueB)
    IF (associated(INFO)) GOTO 990

    IF (KSA == 1_ik .AND. NEVAL == 0) GOTO 999 ! exit if NEVAL=0

    !! Pointers in RWA
    IPA   = 1
    IPB   = IPA + 1
    IPD   = IPB + 1
    IPS   = IPD + 1
    IPK   = IPS + 1
    IPKM1 = IPK + 1
    IPKP1 = IPKM1 + 1

    allocate(LWA(MAXLAN,2),RWA(MAXLAN,7),EVERR(MAXLAN),stat=IERR)
    IF (IERR /= 0) GOTO 980

    LWA = 0_ik
    call DCOPY (MAXLAN*7,0.0_wp,0,RWA(1,1),1)
    call DCOPY (MAXLAN,0.0_wp,0,EVL(1),1)

    !!--------------------------------------------------------------------------
    !! Generate Lanczos vectors and reduce the inverse eigenproblem to
    !! tridiagonal form, of dimension K by K, and if  solveAllInOneGo = .FALSE.
    !! solve this when K = NEVAL, NEVAL+INT(3*NEVAL/2) and then for each KEX'th
    !! step until NEVAL eigenvalues are accepted (twice!)
    !! If  solveAllInOneGo = .TRUE.  tridiagonalization is completed and the
    !! tridiagonal problem is solved only once (i.e. complete solution)
    !!
    !! K = 1, 2, ......, M    where M is less than or equal to MAXLAN
    !!--------------------------------------------------------------------------

    KEVEX = NEVAL
    solveAllInOneGo = MAXLAN == N .AND. N < 25

    IF (MIP(3) > 0) THEN ! premultiply start vector by matrix H

       IF (KPD == 1_ik) THEN
          IFLAG = 5
       ELSE
          IFLAG = 3
       END IF
       if (abs(KSB) == 2_ik) then
          !! Matrix A is GSF and matrix B is diagonal
          call FESMM (INFO,LLP,IFLAG,N,AllInt,R,R,MsgA,structA,factorA, &
               &      diagB=diagB)
       else if (abs(KSB) == 3_ik) then
          !! Matrix A is GSF and matrix B is SPR
          call FESMM (INFO,LLP,IFLAG,N,AllInt,R,R,MsgA,structA,factorA, &
               &      MPARB=MPARB,MTREEB=MTREEB,MSIFB=MSIFB,valueB=valueB)
       else
          !! Both matrices are GSF
          call FESMM (INFO,LLP,IFLAG,N,AllInt,R,R,MsgA,structA,factorA, &
               &      structB,coeffB,factorB)
       end if
       IF (associated(INFO)) THEN
          CALL SPRLER (24_ik,J,J,J,LLLP,LERR)
          GOTO 989
       END IF

    END IF

    !! Length of R
    AUX = DDOT(N,R(1),1,R(1),1)
    IF (AUX < SMALL) THEN
       call FERandomVec (R,RAN)
       AUX = DDOT(N,R(1),1,R(1),1)
    END IF
    BETA = SQRT(AUX)
    IF (KPD == 1_ik) THEN
       IFLAG = 6
    ELSE
       IFLAG = 4
    END IF
    RWA(1,IPB) = 1.0_wp

    reStart = .false.
    IF (KORT == -2_ik) THEN
       JORT = -1
    ELSE
       JORT = 0
    END IF

    call writeProgress ('     Starting Lanczos loop')
    do K = 1, MAXLAN ! ====================================== Start Lanczos loop
       J = int(K,ik)
       KM1 = K - 1
       call ffa_startTimer ('LanczosLoop')
       IF (.not. solveAllInOneGo) KEVEX = KEVEX-1

       !! Form Lanczos vector
       AUX = 1.0_wp / BETA
       call DMVV (N,R(1),1,AUX,0,EVC(1,K))

       IF (KM1 > 0 .AND. .not.reStart) THEN

          !! --- Reorthogonalization - if necessary and/or if specified --------

          IF (KORT > 0_ik) THEN
             MOP(8) = MOP(8) + 1_ik
             DO I = 1, int(KORT)
                call FEMGSOrth (EVC(1:N,:),EVC(1:N,K),KM1)
                IF (IPSW > 1) WRITE (LPU,6010) K
                MOP(9) = MOP(9) + int(KM1,ik)
             END DO
          ELSE IF (KORT == -1_ik .AND. K > 2) THEN
             IF (JORT == 1) THEN
                call FEMGSOrth (EVC(1:N,:),EVC(1:N,K),KM1)
                IF (IPSW > 1) WRITE (LPU,6010) K
                MOP(8) = MOP(8) + 1_ik
                MOP(9) = MOP(9) + int(KM1,ik)
                call DCOPY (K,RLPR,0,RWA(1,IPK),1)
                call DCOPY (K,RLPR,0,RWA(1,IPKM1),1)
                JORT = 0
             ELSE
                call FEReOrl (RWA(1:KM1,IPA),RWA(1:KM1,IPB),RWA(1:KM1,IPK), &
                     &        RWA(1:KM1,IPKM1),RWA(1:KM1,IPKP1), &
                     &        BETA,RLPR,TRSH,KM1,VPMAX)
                IF (IPSW > 1) WRITE (LPU,6020) K,VPMAX
                IF (VPMAX > TRSH) THEN
                   call FEMGSOrth (EVC(1:N,:),EVC(1:N,K),KM1)
                   IF (IPSW > 1) WRITE (LPU,6010) K
                   MOP(8) = MOP(8) + 1_ik
                   MOP(9) = MOP(9) + int(KM1,ik)
                   JORT   = 1
                END IF
             END IF
          ELSE IF (KORT == -2_ik) THEN
             JORT = JORT+1
             IF (JORT > 0) THEN
                call FEMGSOrth (EVC(1:N,:),EVC(1:N,K),KM1)
                IF (IPSW > 1) WRITE (LPU,6010) K
                MOP(8) = MOP(8) + 1_ik
                MOP(9) = MOP(9) + int(KM1,ik)
             END IF
             IF (JORT == 2) JORT = -2
          END IF

       END IF

       !! Prepare new R
       if (abs(KSB) == 2_ik) then
          !! Matrix A is GSF and matrix B is diagonal
          call FESMM (INFO,LLP,IFLAG,N,AllInt,EVC(1:N,K),R, &
               &      MsgA,structA,factorA,diagB=diagB)
       else if (abs(KSB) == 3_ik) then
          !! Matrix A is GSF and matrix B is SPR
          call FESMM (INFO,LLP,IFLAG,N,AllInt,EVC(1:N,K),R, &
               &      MsgA,structA,factorA, &
               &      MPARB=MPARB,MTREEB=MTREEB,MSIFB=MSIFB,valueB=valueB)
       else
          !! Both matrices are GSF
          call FESMM (INFO,LLP,IFLAG,N,AllInt,EVC(1:N,K),R, &
               &      MsgA,structA,factorA,structB,coeffB,factorB)
       end if
       IF (associated(INFO)) THEN
          CALL SPRLER (24_ik,J,J,J,LLLP,LERR)
          GOTO 989
       END IF

       IF (KM1 > 0 .AND. .not.reStart) THEN
          call DAXPY (N,-BETA,EVC(1,KM1),1,R(1),1)
       END IF
       reStart = .false.

       !! Element alpha_k of tridiagonal form
       ALFA = DDOT(N,EVC(1,K),1,R(1),1)
       RWA(K,IPA) = ALFA

       NLV = K

       IF (KEVEX == 0 .OR. K == MAXLAN) THEN

          !! --- Solve reduced, tridiagonal eigenproblem -----------------------

          IF (NEVEXT > 0) LWA(:,1) = LWA(:,2)

          CALL SPRLTD (NLV,LWA(1,2),LWA(1,1),EVERR(1),EVL(1), &
               &       RWA(1,IPA),RWA(1,IPB),RWA(1,IPD),RWA(1,IPS), &
               &       TOL(1),NGDEV,NACEV,int(IPSW,ik),int(LPU,ik),LERR)
          IF (LERR < 0_ik) GOTO 989

          !! Next eigenvalue extraction
          NEVEXT = NEVEXT+1
          IF (NEVEXT == 1 .AND. NEVAL > 1) THEN
             KEVEX = NEVAL/2
          ELSE
             KEVEX = int(KEX)
          END IF

          call writeProgress (int(NGDEV),NEVAL)
          IF (int(NGDEV) >= NEVAL .OR. K == MAXLAN) exit ! loop exit
       END IF

       !! --- Next (unscaled) Lanczos vector -----------------------------------
       call DAXPY (N,-ALFA,EVC(1,K),1,R(1),1)

       !! Length of R
       BETA = sqrt(DDOT(N,R(1),1,R(1),1))
       AUX  = abs(RWA(1,IPA))

       IF (BETA < RLPR*AUX) THEN

          !! ------ Breakdown; generate new (random) start vector and make it
          !!        ortonormal to previous Lanczos vectors, and start again
          call FERandomVec (R,RAN)
          AUX = 1.0_wp / sqrt(DDOT(N,R(1),1,R(1),1))
          call DMVV (N,R(1),1,AUX,0,EVC(1,K))
          call FEMGSOrth (EVC(1:N,:),R,K)
          BETA = 1.0_wp
          RWA(K,IPB) = 0.0_wp

          MOP(10) = MOP(10) + 1_ik
          IF (MOP(10) > 10_ik) THEN
             CALL SPRLER (31_ik,J,J,J,LLLP,LERR)
             cycle
          END IF
          reStart = .true.
          cycle
       END IF

       RWA(K+1,IPB) = BETA
       call writeProgress (K,MAXLAN)

       call ffa_stopTimer ('LanczosLoop')
    end do ! =============================================== End of Lanczos loop

    call ffa_stopTimer ('LanczosLoop')
    if (int(NGDEV) < NEVAL) call writeProgress (NEVAL,NEVAL)
    call writeProgress ('     Done')

    SMALL = RLPR*ABS(EVL(1))

    IF (NGDEV == 0_ik) THEN
       CALL SPRLER (32_ik,int(MAXLAN,ik),J,J,LLLP,LERR)
       GOTO 989
    ELSE IF (int(NGDEV) < NEVAL) THEN
       CALL SPRLER (1_ik,NGDEV,int(MAXLAN,ik),int(NEVAL,ik),LLLP,LERR)
    END IF

    !! Restore eigenvalues Lambda_i
    DO K = 1, NLV
       IF (ABS(EVL(K)) > SMALL) THEN
          EVL(K) = 1.0_wp/EVL(K)
       ELSE
          IF (K-1 < int(NGDEV)) NGDEV = int(K-1,ik)
          IF (K-1 < NEVAL) THEN
             CALL SPRLER (3_ik,int(K-1,ik),J,J,LLLP,LERR)
             EXIT
          END IF
       END IF
    END DO

    MOP(1) = NGDEV
    MOP(2) = NACEV
    MOP(3) = int(NLV,ik)
    MOP(5) = MIN(int(NEVEC,ik),NGDEV)
    MOP(7) = NEVEXT

    !! Print out eigenvalues
    IF (IPSW > 0) THEN
       CALL SPRLP2 (EVL(1),EVERR(1),TOL(1),LWA(1,2), &
            &       MOP(1),int(NLV,ik),int(LPU,ik))
    END IF
    IF (NEVEC <= 0) GOTO 999

    !! --- Determine (B-orthonormal) approx. to the NVC first eigenvectors
    NVC = int(MOP(5))
    allocate(RWE(NLV,NLV),stat=IERR)
    IF (IERR /= 0) GOTO 980

    !! Re-solve small eigenproblem
    CALL QLS3DV (RWA(1,IPA),RWA(1,IPB),RWE(1,1),NLV,LLP,IERR)
    IF (IERR < 0) THEN
       CALL SPRLER (25_ik,int(NLV,ik),J,J,LLLP,LERR)
       GOTO 989
    END IF

    !! Rearrange ordering
    CALL EVARR (RWA(1,IPA),RWE(1,1),NLV,NLV,NVC,4)

    !! --- Restore eigenvectors of special problem -----------------------------

    call writeProgress ('     Computing eigenvectors')

    IF (KPD == 1_ik) THEN

       !! Scale eigenvectors
       SIGMA = minval(RWA(1:NVC,IPA))
       IF (SIGMA < TRSH) THEN
          SIGMA = ABS(SIGMA) + TRSH
          call DAXPY (NVC,1.0_wp,SIGMA,0,RWA(1,IPA),1)
       ELSE
          SIGMA = 0.0_wp
       END IF
       DO K = 1, NVC
          call DSCAL (NLV,1.0_wp/SQRT(RWA(K,IPA)),RWE(1,K),1)
          call DAXPY (NVC,-1.0_wp,SIGMA,0,RWA(1,IPA),1)
       END DO

    END IF

    call DGEMM ('N','N',N,NVC,NLV,1.0_wp,EVC(1,1),size(EVC,1), &
         &      RWE(1,1),NLV,0.0_wp,work1(1,1),N)

    !! --- Restore desired eigenvectors ----------------------------------------

    IF (KPD == 1_ik) THEN

#if FT_DEBUG > 1
       IF (LLP > 0) WRITE(LLP,*) 'CALLING  FESolve(A,Action=BPX'//AllInt//')'
#endif
       call ffa_startTimer ('FESolve')
       call FESolve (MsgA,'BPX'//AllInt,N,NVC,structA,factorA,work1,EVC,INFO)
       call ffa_stopTimer ('FESolve')
       IF (associated(INFO)) THEN
          CALL SPRLER (27_ik,J,J,J,LLLP,LERR)
          GOTO 989
       END IF

    ELSE IF (KSB == -1_ik) THEN

       !! Matrix B is GSF
#if FT_DEBUG > 1
       IF (LLP > 0) WRITE(LLP,*) 'CALLING  FESolve(B,Action=BPX'//AllInt//')'
#endif
       call SetLinearSolverType (MsgB,1)
       call ffa_startTimer ('FESolve')
       call FESolve (MsgB,'BPX'//AllInt,N,NVC,structB,factorB,work1,EVC,INFO)
       call ffa_stopTimer ('FESolve')
       IF (associated(INFO)) THEN
          CALL SPRLER (27_ik,J,J,J,LLLP,LERR)
          GOTO 989
       END IF
       call SetLinearSolverType (MsgA,4-int(KSA))

    ELSE IF (KSB == -3_ik) THEN

       !! Matrix B is SPR
       if (size(RWA,kind=ik) < MPARB(14)) then
          deallocate(RWA)
          allocate(RWA(MPARB(14),1_ik),stat=IERR)
          if (IERR /= 0) goto 980
       end if
       CALL SPRBS2 (MPARB(1),MTREEB(1),MSIFB(1),valueB(1),work1(1,1), &
            &       int(N,ik),int(NVC,ik),RWA(1,1),LLLP,LERR)
       IF (LERR < 0_ik) THEN
          CALL SPRLER (27_ik,J,J,J,LLLP,LERR)
          GOTO 989
       END IF

       call FEPermuteVectors (AllInt,'F',N,NVC,structA,work1,EVC,INFO)

    ELSE IF (MOP(6) == 0_ik) THEN

       !! Matrix B is diagonal
       DO K = 1, N
          R(K) = 1.0_wp / diagB(K)
       END DO
       DO K = 1, NVC
          call DMVV (N,work1(1,K),1,R(1),1,EVC(1,K))
       END DO

       call FEPermuteVectors (AllInt,'F',N,NVC,structA,EVC,INFO)

    ELSE

       DO K = 1, NVC
          call DCOPY (N,work1(1,K),1,R(1),1)
          call FESMM (INFO,LLP,2,N,AllInt,R,EVC(1:N,K), &
               &      MsgA,structA,factorA,diagB=diagB)
          IF (associated(INFO)) THEN
             CALL SPRLER (26_ik,int(K,ik),J,J,LLLP,LERR)
             GOTO 989
          END IF
          call DSCAL (N,1.0_wp/RWA(K,IPA),EVC(1,K),1)
       END DO

       call FEPermuteVectors (AllInt,'F',N,NVC,structA,EVC(1:N,1:NVC),INFO)

    END IF

    goto 999

    !! ------------------------------------------------ Exit -------------------

980 CONTINUE
    call FESetError (INFO,'FELanczosSolve',-1)
    goto 999
989 CONTINUE
    IERR = int(LERR)
990 CONTINUE
    call FESetError (INFO,'FELanczosSolve')
999 CONTINUE
    if (allocated(R))     deallocate(R)
    if (allocated(LWA))   deallocate(LWA)
    if (allocated(RWA))   deallocate(RWA)
    if (allocated(RWE))   deallocate(RWE)
    if (allocated(EVERR)) deallocate(EVERR)
    if (allocated(work1)) deallocate(work1)
    if (allocated(work2)) deallocate(work2)
#if FT_DEBUG > 1
    IF (LLP > 0) WRITE(LLP,*) 'LEAVING  FELanczosSolve'
#endif

6010 FORMAT (5X,'Lanczos step',I5,' :  reorthogonalize')
6020 FORMAT (5X,'Lanczos step',I5,' :  max. inner product =',1PE12.4)

  contains

    subroutine FEMGSOrth (G,V,M)
      real(wp), intent(in)    :: G(:,:)
      real(wp), intent(inout) :: V(:)
      integer , intent(in)    :: M
      integer                 :: j
      real(wp)                :: x
#if FT_DEBUG > 1
      IF (LLP > 0) WRITE(LLP,*) 'ENTERING FEMGSOrth',M
#endif
      do j = 1, M
         x = DDOT(N,V(1),1,G(1,j),1)
         call DAXPY (N,-x,G(1,j),1,V(1),1)
      end do
      x = 1.0_wp / sqrt(DDOT(N,V(1),1,V(1),1))
      call DSCAL (N,x,V(1),1)
    end subroutine FEMGSOrth

    subroutine FERandomVec (V,x)
      real(wp), intent(out)   :: V(:)
      real(wp), intent(inout) :: x
      real(wp), parameter     :: pi = 3.141592653589793_wp
      integer                 :: i, n
      n = size(V)
      do i = 1, n
         if (abs(x) < 1.0e-15_wp) x = 0.1_wp * pi
         x = x*x / real(n,wp)
         do while (x < 1.0_wp)
            x = 10.0_wp * x
         end do
         x = x - real(int(x),wp)
         if (x > 0.5_wp) x = x - 1.0_wp
         v(i) = x
      end do
    end subroutine FERandomVec

    subroutine FEReOrl (ALFA,BETA,VVK,VVKM1,VVKP1,BTAKP1,RLPR,TRSH,K,VPMAX)
      real(wp), intent(in)    :: ALFA(:), BETA(:)
      real(wp), intent(inout) :: VVK(:), VVKM1(:), VVKP1(:)
      real(wp), intent(in)    :: BTAKP1, RLPR, TRSH
      integer , intent(in)    :: K
      real(wp), intent(out)   :: VPMAX
      integer                 :: i
      real(wp)                :: x

#if FT_DEBUG > 1
      IF (LLP > 0) WRITE(LLP,*) 'ENTERING FEReOrl',K
#endif

      VVK(K)     = 1.0_wp
      VVKM1(K-1) = 1.0_wp
      IF (K == 2) VVK(1) = RLPR
      VPMAX = RLPR

      X = ((ALFA(1)-ALFA(K))*VVK(1) + BETA(2)*VVK(2) - BETA(K)*VVKM1(1))/BTAKP1
      IF (abs(X) > VPMAX) VPMAX = abs(X)
      VVKP1(1) = X

      IF (K > 2) THEN
         do I = 2, K-1
            X = (  BETA(I)*VVK(I-1) + (ALFA(I)-ALFA(K))*VVK(I) &
                 + BETA(I+1)*VVK(I+1) - BETA(K)*VVKM1(I) )/ BTAKP1
            IF (abs(X) > VPMAX) VPMAX = abs(X)
            VVKP1(I) = X
         end do
      END IF

      VVKP1(K) = RLPR
      IF (VPMAX <= TRSH) THEN ! prepare for next step
         VVKM1 = VVK
         VVK   = VVKP1
      END IF

    end subroutine FEReOrl

  end subroutine FELanczosSolve


  subroutine FELFac (INFO,LPU,KSA,KSB,KPD,N,NEVAL,EPS,NN,MSING, &
       &             MsgA,structA,coeffA,factorA, &
       &             MsgB,structB,coeffB,factorB,diagB, &
       &             MPARB,MTREEB,MSIFB,valueB)

    !!==========================================================================
    !! Administer the matrix factorizations needed for the Lanczos eigensolver.
    !!
    !! Programmer : Knut Morten Okstad                   date/rev : Sep 2004/1.0
    !!==========================================================================

    use SprKindModule , only : ik
    use FELinearSolver, only : MessageLevel, SAM, CAM, GSFColSup, Error_Flag
    use FELinearSolver, only : SetLinearSolverType, FEFactorize
    use FELinearSolver, only : FEErrorHandler, FESetError, FEPivotError

    use FFaProfilerInterface, only : ffa_startTimer, ffa_stopTImer

    type(Error_Flag)  , pointer       :: INFO
    integer           , intent(in)    :: LPU, KPD, N, NEVAL
    integer(ik)       , intent(inout) :: KSA, KSB
    real(wp)          , intent(in)    :: EPS
    integer(ik)       , intent(out)   :: NN, MSING(:)
    type(MessageLevel), pointer       :: MsgA, MsgB
    type(SAM)         , pointer       :: structA, structB
    type(CAM)         , pointer       :: coeffA, coeffB
    type(GSFColSup)   , pointer       :: factorA, factorB
    integer(ik), optional, intent(in) :: MPARB(:), MTREEB(:), MSIFB(:)
    real(wp), optional, intent(inout) :: diagB(:), valueB(:)

    !! Local variables
    integer                  :: NZP, KTYPE
    integer(ik)              :: I, LERR
    real(wp)                 :: TPR(3)
    real(wp)   , parameter   :: epsZero_p = 1.0e-16_wp
    integer(ik), allocatable :: LWA(:)
    real(wp)   , allocatable :: RWA(:)

    external CMVY18, CMMY18

    !! --- Logic section ---

    call ffa_startTimer ('FELFac')
#if FT_DEBUG > 1
    IF (LPU > 0) WRITE(LPU,*) 'ENTERING FELFac',KSA,KSB,KPD
#endif

    NN   = 0_ik
    NZP  = 0
    LERR = 0_ik

    IF (KSA == 1_ik) THEN ! Factorize matrix A

       call FEFindAndFixZeroPivots (structA,coeffA,NZP,MSING,INFO)
       if (associated(INFO)) then
          LERR = -1_ik
          goto 900
       end if

       IF (KPD == 1 .AND. NEVAL > 0) THEN
          KSA = 3_ik
       ELSE
          KSA = 2_ik
       END IF
       KTYPE = 4 - int(KSA)
#if FT_DEBUG > 1
       IF (LPU > 0) WRITE(LPU,*) 'CALLING  FEFactorize(A)',KTYPE
#endif
       call SetLinearSolverType (MsgA,KTYPE)
       call ffa_startTimer ('FEFactorize')
       call FEFactorize (MsgA,coeffA,factorA,INFO)
       call ffa_stopTimer ('FEFactorize')
       if (associated(INFO)) then
          if (FEPivotError(INFO, KSA == 3_ik)) LERR = -3_ik
          CALL SPRLER (21_ik,I,I,I,int(LPU,ik),LERR)
          if (LERR /= -2_ik) goto 900
       end if
    ELSE
       KTYPE = 4 - int(KSA)
    END IF

    IF (KSB > 0_ik .AND. KPD == 2) THEN ! Factorize matrix B

       IF (KSB == 1_ik) THEN ! Matrix B is GSF

          call FEFindAndFixZeroPivots (structB,coeffB,NZP,MSING,INFO)
          if (associated(INFO)) then
             LERR = -1_ik
             goto 900
          end if

          KSB = -1_ik
#if FT_DEBUG > 1
          IF (LPU > 0) WRITE(LPU,*) 'CALLING  FEFactorize(B)',1
#endif
          call SetLinearSolverType (MsgB,1)
          call ffa_startTimer ('FEFactorize')
          call FEFactorize (MsgB,coeffB,factorB,INFO)
          call ffa_stopTimer ('FEFactorize')
          if (associated(INFO)) then
             if (FEPivotError(INFO)) LERR = -3_ik
             CALL SPRLER (22_ik,KSB,I,I,int(LPU,ik),LERR)
          end if

          !! Reset solver type for matrix A in case MsgA and MsgB are shared
          call SetLinearSolverType (MsgA,KTYPE)

       ELSE IF (KSB == 3_ik) THEN ! Matrix B is SPR

          allocate(LWA(MPARB(14)+1_ik),RWA(MPARB(17)),stat=LERR)
          if (LERR /= 0_ik) then
             LERR = -1_ik
             call FESetError (INFO,'FELFac',-1)
             goto 900
          end if

          KSB = -3
          TPR(1) = EPS
          CALL SPRFC2 (MPARB(1),MTREEB(1),MSIFB(1),valueB(1),TPR(1), &
               &       LWA(1),RWA(1),int(LPU,ik),LERR,CMVY18,CMMY18)
          IF (LERR < 0_ik) CALL SPRLER (22_ik,int(KSB,ik),I,I,int(LPU,ik),LERR)

          IF (LERR < 0_ik .AND. MPARB(27) > 0_ik .AND. NZP < size(MSING)) THEN
             MSING(NZP+1) = MPARB(27)
          END IF

          deallocate(LWA,RWA)

       ELSE ! Matrix B is diagonal

          KSB = -2
          DO I = 1, int(N,ik)
             IF (diagB(I) > epsZero_p) THEN
                diagB(I) = SQRT(diagB(I))
             ELSE IF (diagB(I) >= -epsZero_p) THEN
                NN = NN + 1_ik
             ELSE
                LERR = -3_ik
                CALL SPRLER (22_ik,int(KSB,ik),I,I,int(LPU,ik),LERR)
                EXIT
             END IF
          END DO

       END IF

    ELSE IF (ABS(KSB) == 2) THEN ! B is diagonal and needs not to be factorized

       DO I = 1, int(N,ik)
          IF (abs(diagB(I)) <= epsZero_p) NN = NN + 1_ik
       END DO

    END IF

900 if (LERR < 0_ik) then
       call FESetError (INFO,'FELFac')
    else if (associated(INFO)) then
       call FEErrorHandler (INFO, output=LPU)
       nullify(INFO)
    end if

#if FT_DEBUG > 1
    IF (LPU > 0) WRITE(LPU,*) 'LEAVING  FELFac',KSA,KSB
#endif
    call ffa_stopTimer ('FELFac')

  end subroutine FELFac


  subroutine FEFindAndFixZeroPivots (struct,coeff,NZERO,MSING,INFO)

    !!==========================================================================
    !! Repair isolated zero pivot elements in the matrix prior to factorization.
    !!
    !! Programmer : Knut Morten Okstad                   date/rev : Sep 2005/1.0
    !!==========================================================================

    use SprKindModule , only : ik
    use FELinearSolver, only : SAM, CAM, Error_Flag
    use FELinearSolver, only : FEPermuteIndices, FESetError
    use FELinearSolver, only : FEFindIsolatedZeroPivots, FEInsertPivots

    type(SAM)       , pointer       :: struct
    type(CAM)       , pointer       :: coeff
    integer         , intent(inout) :: NZERO
    integer(ik)     , intent(inout) :: MSING(:)
    type(Error_Flag), pointer       :: INFO

    !! Local variables
    integer          :: i, NZST
    integer(ik)      :: iEq
    integer, pointer :: iWork(:)

    !! --- Logic section ---

    nullify(iWork)
    call FEFindIsolatedZeroPivots (coeff,iWork,INFO)
    if (Get_ErrStat()) return
    if (.not. associated(iWork)) return

    !! Isolated zero pivots were found,
    !! fix them by inserting 1.0 on the diagonal
    call DCOPY (size(iWork),1.0_wp,0,work1(1,1),1)
    call FEInsertPivots (coeff,iWork,work1(:,1),INFO)
    if (Get_ErrStat()) return

    !! Convert to external equation numbers
    call FEPermuteIndices ('F',struct,iWork,INFO)
    if (Get_ErrStat()) return

    !! Store the singular equation numbers in MSING, making sure the same
    !! equation is not counted twice (in case both matrices are singular)
    NZST = min(NZERO,size(MSING))
    do i = 1, size(iWork)
       iEq = int(iWork(i),ik)
       if (NZST > 0) then
          if (count(MSING(1:NZST) == iEq) > 0) cycle
       end if
       NZERO = NZERO + 1
       if (NZERO > size(MSING)) cycle
       NZST = NZERO
       MSING(NZST) = iEq
    end do

    deallocate(iWork)

  contains

    function Get_ErrStat ()
      logical :: Get_ErrStat
      Get_ErrStat = associated(INFO)
      if (Get_ErrStat) then
         call FESetError (INFO,'FEFindAndFixZeroPivots')
      end if
    end function Get_ErrStat

  end subroutine FEFindAndFixZeroPivots


  subroutine FESMM (INFO,LPU,IFLAG,Ndim,AllInt,X,Y, &
       &            MsgA,structA,factorA, &
       &            structB,coeffB,factorB,diagB, &
       &            MPARB,MTREEB,MSIFB,valueB)

    !!==========================================================================
    !! To determine:
    !!
    !!    X := (A-1)*UbT*X         IFLAG = 1
    !!    Y  = (A-1)*UbT*X         IFLAG = 2
    !!    X := Ub*(A-1)*UbT*X      IFLAG = 3
    !!    Y  = Ub*(A-1)*UbT*X      IFLAG = 4
    !!    X := (Ua-T)*B*(Ua-1)*X   IFLAG = 5
    !!    Y  = (Ua-T)*B*(Ua-1)*X   IFLAG = 6
    !!
    !!  where X and Y are vectors,
    !!        A is a symmetric sparse matrix, represented by
    !!        its factors La*Da*LaT, stored in (structA,factorA),
    !!        Ua is the Cholesky factor of A, stored in (structA,factorA),
    !!        B is either a symmetric sparse matrix, stored in (structB,coeffB),
    !!        a symmetric sparse in-core matrix, stored in (MPARB,...,valueB),
    !!        or a diagonal matrix, stored in diagB.
    !!        Ub is the Cholesky factor of B, stored in (structB,factorB), or in
    !!        (MPARB,...,valueB), or in diagB if B is a diagonal matrix.
    !!
    !! Programmer : Knut Morten Okstad                   date/rev : Jan 2005/1.0
    !!==========================================================================

    use SprKindModule , only : ik
    use FELinearSolver, only : MessageLevel, SAM, CAM, GSFColSup, Error_Flag
    use FELinearSolver, only : FESolve, FETRMM, FEPRM, FESetError

    use FFaProfilerInterface, only : ffa_startTimer, ffa_stopTImer

    type(Error_Flag)         , pointer       :: INFO
    integer                  , intent(in)    :: LPU, IFLAG, Ndim
    character(len=1)         , intent(in)    :: AllInt
    real(wp)                 , intent(inout) :: X(:)
    real(wp)                 , intent(out)   :: Y(:)
    type(MessageLevel)       , pointer       :: MsgA
    type(SAM)                , pointer       :: structA
    type(GSFColSup)          , pointer       :: factorA
    type(SAM)      , optional, pointer       :: structB
    type(CAM)      , optional, pointer       :: coeffB
    type(GSFColSup), optional, pointer       :: factorB
    real(wp)       , optional, intent(in)    :: diagB(:), valueB(:)
    integer(ik)    , optional, intent(in)    :: MPARB(:), MTREEB(:), MSIFB(:)

    !! Local variables
    integer(ik)       :: IERR
    real(wp), pointer :: xMat(:,:)

    !! --- Logic section ---

    call ffa_startTimer ('FESMM')
#if FT_DEBUG > 1
    IF (LPU > 0) WRITE(LPU,*) 'ENTERING FESMM, IFLAG =',IFLAG
    IF (LPU > 0) WRITE(LPU,"(5X,'X =',1P6E12.4/(8X,6E12.4))") X
#endif

    IERR = 0_ik
    select case (IFLAG)

    case(1:4) ! Compute (A-1)*UbT*X

       if (present(diagB)) then
          call DMVV (Ndim,diagB(1),1,X(1),1,work1(1,1))
       else if (present(structB) .and. present(factorB)) then
#if FT_DEBUG > 1
          IF (LPU > 0) WRITE(LPU,*) 'CALLING  FETRMM(B,'//AllInt//',N,N)'
#endif
          xMat => MatrixPointToArray(X,Ndim,1)
          call ffa_startTimer ('FETRMM')
          call FETRMM (AllInt,'N','N','N',Ndim,1,structB,1.0_wp,factorB, &
               &       xMat,work1,INFO)
          call ffa_stopTimer ('FETRMM')
          if (associated(INFO)) goto 90
       else if (present(MPARB) .and. present(valueB)) then
          if (.not.present(MTREEB) .or. .not.present(MTREEB)) then
             call FESetError (INFO,'FESMM',-99)
             goto 100
          end if
          call SPRTMM (valueB(MPARB(8)+1),X(1),X(1),work1(1,1), &
               &       MPARB(1),MTREEB(1),MSIFB(1),int(Ndim,ik), &
               &       1_ik,1_ik,22_ik,int(LPU,ik),IERR)
          if (IERR < 0_ik) goto 90
       else
          call FESetError (INFO,'FESMM',-99)
          goto 100
       end if

#if FT_DEBUG > 1
       IF (LPU > 0) WRITE(LPU,*) 'CALLING  FESolve(A,Action=SNB'//AllInt//')'
#endif
       call ffa_startTimer ('FESolve')
       call FESolve (MsgA,'SNB'//AllInt,Ndim,1,structA,factorA,work1,work2,INFO)
       call ffa_stopTimer ('FESolve')
       if (associated(INFO)) goto 90

       select case (IFLAG)

       case(1) ! Store result in X (overwriting the input)

          call DCOPY (Ndim,work2(1,1),1,X(1),1)

       case(2) ! Store result in Y

          call DCOPY (Ndim,work2(1,1),1,Y(1),1)

       case(3:4) ! Compute Ub*(A-1)*UbT*X

          if (present(diagB)) then
             call DMVV (Ndim,diagB(1),1,work2(1,1),1,work1(1,1))
          else if (present(structB) .and. present(factorB)) then
#if FT_DEBUG > 1
             IF (LPU > 0) WRITE(LPU,*) 'CALLING  FETRMM(B,'//AllInt//',T,N)'
#endif
             call ffa_startTimer ('FETRMM')
             call FETRMM (AllInt,'T','N','N',Ndim,1,structB,1.0_wp,factorB, &
                  &       work2,work1,INFO)
             call ffa_stopTimer ('FETRMM')
             if (associated(INFO)) goto 90
          else
             call SPRTMM (valueB(MPARB(8)+1),work2(1,1),X(1),work1(1,1), &
                  &       MPARB(1),MTREEB(1),MSIFB(1),int(Ndim,ik), &
                  &       1_ik,1_ik,11_ik,int(LPU,ik),IERR)
             if (IERR < 0_ik) goto 90
          end if

       end select

    case(5:6) ! Compute (Ua-T)*B*(Ua-1)*X

#if FT_DEBUG > 1
       IF (LPU > 0) WRITE(LPU,*) 'CALLING  FESolve(A,Action=BNB'//AllInt//')'
#endif
       xMat => MatrixPointToArray(X,Ndim,1)
       call ffa_startTimer ('FESolve')
       call FESolve (MsgA,'BNB'//AllInt,Ndim,1,structA,factorA,xMat,work1,INFO)
       call ffa_stopTimer ('FESolve')
       if (associated(INFO)) goto 90

       if (present(diagB)) then
          call DMVV (Ndim,diagB(1),1,work1(1,1),1,work2(1,1))
       else if (present(structB) .and. present(coeffB)) then
#if FT_DEBUG > 1
          IF (LPU > 0) WRITE(LPU,*) 'CALLING  FEPRM(B,'//AllInt//')'
#endif
          call ffa_startTimer ('FEPRM')
          call FEPRM (AllInt,'N',Ndim,1,structB,1.0_wp,coeffB,work1, &
               &      0.0_wp,work2,INFO)
          call ffa_stopTimer ('FEPRM')
          if (associated(INFO)) goto 90
       else if (present(MPARB) .and. present(valueB)) then
          if (.not.present(MTREEB) .or. .not.present(MTREEB)) then
             call FESetError (INFO,'FESMM',-98)
             goto 100
          end if
          call SPRPRM (valueB(MPARB(8)+1),work1(1,1),X(1),work2(1,1), &
               &       MPARB(1),MTREEB(1),MSIFB(1),int(Ndim,ik), &
               &       1_ik,1_ik,88_ik,int(LPU,ik),IERR)
          if (IERR < 0_ik) goto 90
       else
          call FESetError (INFO,'FESMM',-98)
          goto 100
       end if

#if FT_DEBUG > 1
       IF (LPU > 0) WRITE(LPU,*) 'CALLING  FESolve(A,Action=FNB'//AllInt//')'
#endif
       call ffa_startTimer ('FESolve')
       call FESolve (MsgA,'FNB'//AllInt,Ndim,1,structA,factorA,work2,work1,INFO)
       call ffa_stopTimer ('FESolve')
       if (associated(INFO)) goto 90

    case default
       call FESetError (INFO,'FESMM',-97)
       goto 100
    end select

    select case (IFLAG)
    case(3,5) ! Store result in X (overwriting the input)
       call DCOPY (Ndim,work1(1,1),1,X(1),1)
       X = work1(:,1)
    case(4,6) ! Store result in Y
       call DCOPY (Ndim,work1(1,1),1,Y(1),1)
    end select

 90 if (associated(INFO) .or. IERR < 0_ik) call FESetError (INFO,'FESMM')
#if FT_DEBUG > 1
    IF (LPU > 0) WRITE(LPU,"(5X,'Y =',1P6E12.4/(8X,6E12.4))") work1(:,1)
    IF (LPU > 0) WRITE(LPU,*) 'LEAVING  FESMM'
#endif
100 call ffa_stopTimer ('FESMM')

  contains

    function MatrixPointToArray (array,rows,columns)
      integer , intent(in)         :: rows, columns
      real(wp), intent(in), target :: array(rows,columns)
      real(wp), pointer            :: MatrixPointToArray(:,:)
      MatrixPointToArray => array
    end function MatrixPointToArray

  end subroutine FESMM

end module FELanczosModule
