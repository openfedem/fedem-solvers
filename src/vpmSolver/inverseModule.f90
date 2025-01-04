!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module inverseModule

  !!============================================================================
  !! This module contains driver subroutines for solving inverse problems.
  !!============================================================================

  use KindModule, only : dp

  implicit none

  private

  real(dp), allocatable :: U(:,:)     ! displacement (global direction)
  real(dp), allocatable :: Ud(:,:)    ! velocity     (global direction)
  real(dp), allocatable :: Udd(:,:)   ! acceleration (global direction)
  real(dp), allocatable :: cfa(:,:)   ! vector for the HHT-alpha method
  real(dp), allocatable :: F_inv(:,:) ! inverse force vector (needed for output)
  real(dp), allocatable :: Qvec(:)    ! constant force vector (e.g. gravity)
  real(dp), allocatable :: V0(:)      ! temporary array
  real(dp), allocatable :: V1(:)      ! temporary array
  real(dp), allocatable :: V2(:)      ! temporary array
  real(dp), allocatable :: cFact(:)   ! vector for the HHT - alpha method
  real(dp), allocatable :: CV(:)      ! damping part (holds C*V1)
  real(dp), allocatable :: MA(:)      ! acceleration part (holds M*V2)
  real(dp), allocatable :: KU(:)      ! stiffness part (holds K*V0)
  real(dp), allocatable :: Un(:)      ! displacement increment vector
  real(dp), allocatable :: Udn(:)     ! velocity increment vector
  real(dp), allocatable :: Uddn(:)    ! acceleration increment vector
  real(dp), allocatable :: Csc(:)     ! dynamic part of the rhs
  real(dp), allocatable :: Uc(:)      ! displacements caused by dynamic forces
  real(dp), allocatable :: M(:,:)     ! mass matrix
  real(dp), allocatable :: C(:,:)     ! damping matrix
  real(dp), allocatable :: K(:,:)     ! stiffness matrix
#ifdef FT_DEBUG
  real(dp), allocatable :: KH(:,:)    ! newton matrix
  integer , allocatable :: ipiv(:)    ! Pivot indices of permuation matrix P
#endif
  real(dp), allocatable :: FH(:)      ! rhs vector for equation (KH * Un = FH)

  public :: solveInverse, doneInverse


contains

  subroutine prepare_inv_mem (ndim, incs, err)

    !!==========================================================================
    !! memory management - memory allocation
    !!
    !! Programmer : Guenter Glanzer                   date/rev : 10 Jul 2018/1.0
    !!==========================================================================

    use ReportErrorModule, only : allocationError

    integer, intent(in)  :: ndim, incs
    integer, intent(out) :: err

    if (allocated(U)) then
       err = 0
       return
    end if

    allocate(U    (ndim,incs), &
         &   Ud   (ndim,incs), &
         &   Udd  (ndim,incs), &
         &   cfa  (ndim,incs), &
         &   F_inv(ndim,incs), &
         &   Qvec (ndim),      &
         &   V0   (ndim),      &
         &   V1   (ndim),      &
         &   V2   (ndim),      &
         &   cFact(ndim),      &
         &   CV   (ndim),      &
         &   MA   (ndim),      &
         &   KU   (ndim),      &
         &   Un   (ndim),      &
         &   Udn  (ndim),      &
         &   Uddn (ndim),      &
         &   Csc  (ndim),      &
         &   Uc   (ndim),      &
         &   M    (ndim,ndim), &
         &   C    (ndim,ndim), &
         &   K    (ndim,ndim), &
         &   FH   (ndim), STAT=err)

#ifdef FT_DEBUG
    if (err == 0) then
       allocate(KH(ndim,ndim),ipiv(ndim), STAT=err)
    end if
#endif

    if (err /= 0) then
       err = AllocationError('prepare mem in inverseModule')
       return
    end if

    !! initialisation
    U     = 0.0_dp;  Ud   = 0.0_dp;  Udd = 0.0_dp;  cfa   = 0.0_dp
    F_inv = 0.0_dp;  Qvec = 0.0_dp;  FH  = 0.0_dp
    M     = 0.0_dp;  C    = 0.0_dp;  K   = 0.0_dp
    V0    = 0.0_dp;  V1   = 0.0_dp;  V2  = 0.0_dp;  cFact = 0.0_dp;  CV = 0.0_dp
    MA    = 0.0_dp;  KU   = 0.0_dp;  Un  = 0.0_dp;  Udn = 0.0_dp
    Uddn  = 0.0_dp;  Csc  = 0.0_dp;  Uc  = 0.0_dp

  end subroutine prepare_inv_mem


  subroutine doneInverse ()

    !!==========================================================================
    !! memory management - deallocation
    !!
    !! Programmer : Guenter Glanzer                   date/rev : 10 Jul 2018/1.0
    !!==========================================================================

    if (allocated(U))     deallocate(U)
    if (allocated(Ud))    deallocate(Ud)
    if (allocated(Udd))   deallocate(Udd)
    if (allocated(cfa))   deallocate(cfa)

    if (allocated(F_inv)) deallocate(F_inv)
    if (allocated(Qvec))  deallocate(Qvec)

    if (allocated(K))     deallocate(K)
    if (allocated(C))     deallocate(C)
    if (allocated(M))     deallocate(M)

    if (allocated(V0))    deallocate(V0)
    if (allocated(V1))    deallocate(V1)
    if (allocated(V2))    deallocate(V2)

    if (allocated(cFact)) deallocate(cFact)
    if (allocated(CV))    deallocate(CV)
    if (allocated(MA))    deallocate(MA)
    if (allocated(KU))    deallocate(KU)

    if (allocated(Un))    deallocate(Un)
    if (allocated(Udn))   deallocate(Udn)
    if (allocated(Uddn))  deallocate(Uddn)
    if (allocated(Csc))   deallocate(Csc)
    if (allocated(Uc))    deallocate(Uc)
    if (allocated(FH))    deallocate(FH)

#ifdef FT_DEBUG
    if (allocated(KH))    deallocate(KH)
    if (allocated(ipiv))  deallocate(ipiv)
#endif

  end subroutine doneInverse


  subroutine solveInverse (x,x_eqs,f_eqs,finished,ierr)

    !!==========================================================================
    !! Solves the inverse problem for next time step.
    !!
    !! Programmer : Knut Morten Okstad                date/rev : 23 Jan 2018/1.0
    !!==========================================================================

    use SolverModule     , only : systemSize, solveStep, solveDynamic
    use SolverModule     , only : getRhsVector, setRhsVector
    use ReportErrorModule, only : reportError, debugFileOnly_p
    use ReportErrorModule, only : allocationError

    real(dp), intent(in)  :: x(:)
    integer , intent(in)  :: x_eqs(:), f_eqs(:)
    logical , intent(out) :: finished
    integer , intent(out) :: ierr

    !! Local variables
    integer               :: iop, ndim
    integer , save        :: iDynStep = 0
    real(dp), allocatable :: Rvec(:)

    !! --- Logic section ---

    !! Allocate right-hand-side vector
    call systemSize (ndim)
    allocate(Rvec(ndim),stat=ierr)
    if (ierr /= 0) then
       ierr = allocationError('solveInverse')
       return
    end if

    iop = -2 ! Build and factor the first coefficient matrix of current step
    call solveStep (iop,.false.,finished,ierr)
    if (ierr < 0) goto 900

    if (finished) then
       if (ierr > 0) then
          goto 910 ! Time integration failure
       else
          goto 890 ! End of time range
       end if
    end if

    !! Extract the current right-hand-side vector (external forces)
    call getRhsVector (Rvec,1,ierr)
    if (ierr < 0) goto 900

    !! Solve the inverse problem
    if (solveDynamic()) then
       iDynStep = iDynStep + 1
       if (iDynStep == 1) then
          !! Solve the first dynamic step quasi-statically
          call coreInverse (x,x_eqs,f_eqs,Rvec,ierr)
       else
          !! Solve the dynamic inverse problem
          call coreInverseDyn (x,x_eqs,f_eqs,Rvec,iDynStep-1,ierr)
       end if
    else
       !! Still in the quasi-static time domain
       call coreInverse (x,x_eqs,f_eqs,Rvec,ierr)
    end if
    if (ierr < 0) goto 900

    !! Update the right-hand-side vector
    call setRhsVector (Rvec,.true.,.true.,ierr)
    if (ierr < 0) goto 900

    !! Do the equilibrium iterations
    call solveStep (iop,.false.,finished,ierr)
    if (ierr < 0) goto 900

890 deallocate(Rvec)
    return

900 call doneInverse ()
910 call reportError (debugFileOnly_p,'solveInverse')
    goto 890

  end subroutine solveInverse


  subroutine coreInverseDyn (x,x_eqs,f_eqs,F,iStep,ierr)

    !!==========================================================================
    !! Solves the dynamic inverse problem for next time step.
    !!
    !! Programmer : Guenter Glanzer                  date/rev : 10 June 2018/1.0
    !!==========================================================================

    use SolverModule     , only : systemSize, solveLinEqSystem
    use SolverModule     , only : setRhsVector, getSystemMatrix
    use ReportErrorModule, only : reportError, debugFileOnly_p

    real(dp), intent(in)    :: x(:)
    real(dp), intent(inout) :: F(:)
    integer , intent(in)    :: x_eqs(:), f_eqs(:), iStep
    integer , intent(out)   :: ierr

    !! Local variables
    integer  :: ndim, n, ns, i0, i1
    real(dp) :: alphaH, a0, a1, a2, a3, a4, a5
    real(dp) :: uc0(size(x_eqs))

    !! --- Logic section ---

    call inv_init (alphaH,a0,a1,a2,a3,a4,a5)

    !! allocate right-hand-side vector
    call systemSize (ndim)
    call prepare_inv_mem (ndim,2,ierr)
    if (ierr < 0) goto 900

    !! system matrices and newton matrix
    call getSystemMatrix (M,2,ierr) ! 2 returns the mass matrix
    if (ierr < 0) goto 900
    call getSystemMatrix (K,1,ierr) ! 1 returns the stiffness matrix
    if (ierr < 0) goto 900
    call getSystemMatrix (C,3,ierr) ! 3 returns damping matrix
    if (ierr < 0) goto 900

#ifdef FT_DEBUG
    !! building newton matrix
    KH = a0*M + (1.0_dp+alphaH)*(a1*C + K)
#endif

    !! save constant force vector
    Qvec(1:ndim) = F(1:ndim)

    !! HHT-alpha/newmark algorithm
    !! displacements/velocities and accelerations are stored in a ndim*2 matrix
    i0 = mod(iStep,2) + 1
    i1 = mod(iStep+1,2) + 1

    V0 = U(:,i0)
    V1 = a1 * U(:,i0) + a4 * Ud(:,i0) + a5 * Udd(:,i0)
    V2 = a0 * U(:,i0) + a2 * Ud(:,i0) + a3 * Udd(:,i0)

    CV = matmul(C,V1)
    MA = matmul(M,V2)
    KU = matmul(K,V0)

    cFact = cfa(:,i0)

    !! inverse part

    Csc = (1.0_dp+alphaH)*CV + MA - alphaH*cFact

#ifdef FT_DEBUG
    call solveEquation (.true.,KH,Csc,Uc,ndim,ierr)
#else
    call solveLinEqSystem (Csc,1,ierr)
    Uc = Csc
#endif
    if (ierr < 0) goto 900

    !! reduce measurements by dynamic proportion
    ns = size(x_eqs)
    do n = 1, ns
       uc0(n) = x(n) - Uc(x_eqs(n))
    end do

    call coreInverse (uc0,x_eqs,f_eqs,F,ierr)
    if (ierr < 0) goto 900

    !! correction is necessary,
    !! because HHT algorithm needs the total force on the structure
    F = F + Qvec

    !! scaling for fedem input
    F = (1.0_dp/(1.0_dp+alphaH))*F

    FH = (1.0_dp+alphaH)*(CV+F) + MA - alphaH*cFact

#ifdef FT_DEBUG
    call solveEquation (.false.,KH,FH,Un,ndim,ierr)
#else
    call solveLinEqSystem (FH,1,ierr)
    Un = FH
#endif
    if (ierr < 0) goto 900

    !! update accelerations and velocities
    Uddn = a0*Un - V2
    Udn  = a1*Un - V1

    cFact = F - matmul(C,Udn) - matmul(K,Un)

    !! store displacements, velocities and accelerations
    U(:,i1)   = Un
    Ud(:,i1)  = Udn
    Udd(:,i1) = Uddn
    cfa(:,i1) = cFact

    !! store F for output
    F_inv(:,i1) = F

#ifdef FT_DEBUG
    !! write solution to file
    call write_results_to_file (iStep,i1,ndim)
#endif

    !! shift generalized force vector by the constant force vector Qvec
    !! (e.g. fedem specification)
    F = F - Qvec

    return

900 call reportError (debugFileOnly_p,'coreInverseDyn')

#ifdef FT_DEBUG
  contains

    subroutine solveEquation (doFactor,A,b,x,n,ierr)
      use ReportErrorModule, only : error_p
      !! Solves a linear system of equations using LAPACK subroutines
      logical,  intent(in)    :: doFactor
      real(dp), intent(inout) :: A(:,:)
      real(dp), intent(in)    :: b(:)
      real(dp), intent(out)   :: x(:)
      integer,  intent(in)    :: n
      integer,  intent(out)   :: ierr
      call DCOPY (n,b(1),1,x(1),1)
      if (doFactor) then
         call DGESV (n,1,A(1,1),size(A,1),ipiv(1),x(1),size(x),ierr)
      else
         call DGETRS ('N',n,1,A(1,1),size(A,1),ipiv(1),x(1),size(x),ierr)
      end if
      if (ierr /= 0) then
         call reportError (error_p,'LAPACK::DGESV returned IERR = ',ierr=ierr)
      end if
    end subroutine solveEquation
#endif

  end subroutine coreInverseDyn


  subroutine inv_init (alphaH, a0, a1, a2, a3, a4, a5)

    !!==========================================================================
    !! Parameters for the HHT-alpha/Newmark scheme
    !!
    !! Programmer : Guenter Glanzer                   date/rev : 10 Jul 2018/1.0
    !!==========================================================================

    use SolverModule, only : solverParameters

    real(dp), intent(out) :: alphaH, a0, a1, a2, a3, a4, a5

    real(dp) :: dt, beta, gamma

    call solverParameters (dt, alphaH, beta, gamma)

    a0 = 1.0_dp/(beta*dt*dt)
    a1 =  gamma/(beta*dt)
    a2 = 1.0_dp/(beta*dt)
    a3 = 0.5_dp/(beta) - 1.0_dp
    a4 =  gamma/(beta) - 1.0_dp
    a5 = 0.5_dp*dt*(gamma/beta - 2.0_dp)

  end subroutine inv_init


#ifdef FT_DEBUG
  subroutine write_results_to_file (seqNr, row_nr, numcols)

    !!==========================================================================
    !! Displacements and inverse forces are written to file (visualisation)
    !! file names: disp.txt
    !!              inv.txt
    !!
    !! Programmer : Guenter Glanzer                date/rev : 10 Jul 2018/1.0
    !!==========================================================================

    use FileUtilitiesModule, only : findUnitNumber

    integer, intent(in) :: seqNr, row_nr, numcols

    character(len=*), parameter :: disp_file = "disp.txt"
    character(len=*), parameter :: inv_file  = "inv.txt"

    integer           :: j, stat, lpu
    character(len=16) :: access

    if (seqNr == 1) then
       access = 'sequential'
    else
       access = 'append'
    end if

    !! displacement part
    lpu = findUnitNumber(35)
    open(lpu,file=disp_file,status='unknown',access=access,iostat=stat)
    if (stat /= 0) then
       write(*,*) 'Open ', disp_file, ' failed with iostat = ', stat
       return
    end if

    write(lpu,*) (U(j,row_nr),",",j=1,numcols)
    close(lpu)

    !! inverse force part
    open(lpu,file=inv_file,status='unknown',access=access,iostat=stat)
    if (stat /= 0) then
       write(*,*) 'Open ', inv_file, ' failed with iostat = ', stat
       return
    end if

    write(lpu,*) (F_inv(j,row_nr),",",j=1,numcols)
    close(lpu)

  end subroutine write_results_to_file
#endif


  subroutine coreInverse (x,x_eqs,f_eqs,F,ierr)

    !!==========================================================================
    !! Solves the linear equation system of the inverse problem.
    !!
    !! Programmer : Knut Morten Okstad                date/rev : 23 Jan 2018/1.0
    !!==========================================================================

    use SolverModule     , only : solveLinEqSystem
    use ReportErrorModule, only : reportError, error_p, debugFileOnly_p
    use ReportErrorModule, only : allocationError, internalError
#ifdef FT_DEBUG
    use manipMatrixModule  , only : writeObject
    use fileUtilitiesModule, only : getDBGfile
#endif

    real(dp), intent(in)    :: x(:)
    integer , intent(in)    :: x_eqs(:), f_eqs(:)
    real(dp), intent(inout) :: F(:)
    integer , intent(out)   :: ierr

    !! Local variables
    integer               :: i, j, ldR, ndis, ndim, nfrc
    real(dp)              :: Lwork
    real(dp), allocatable :: Rhs(:), c(:,:), R(:), Work(:)
    character(len=64)     :: errMsg

#ifdef FT_DEBUG
    integer :: lpu
    lpu = getDBGfile(35,'inverse.dbg')
#endif

    !! --- Logic section ---

    ierr = 0
    ndis = size(x_eqs)
    nfrc = size(f_eqs)
    if (ndis < 1) then
       ierr = internalError('coreInverse: No displacement variables')
       return
    else if (nfrc < 1) then
       ierr = internalError('coreInverse: No force variables')
       return
    else if (size(x) /= ndis) then
       ierr = internalError('coreInverse: Inconsistent input')
       return
    end if

    !! Check that the given equation numbers are in valid range
    ndim = size(F)
    call checkEqns (x_eqs,ierr)
    call checkEqns (f_eqs,ierr)
    if (ierr < 0) return

    i = -1 ! Calculate required size of the Work array
    ldR = max(ndis,nfrc)
    call DGELS ('N',ndis,nfrc,1,Lwork,ndis,Lwork,ldR,Lwork,i,ierr)
    if (ierr /= 0) goto 900

    !! Allocate work arrays
    allocate(Rhs(ndim+ndim*nfrc),c(ndis,nfrc),R(ldR),Work(int(Lwork)),stat=ierr)
    if (ierr /= 0) then
       ierr = allocationError('coreInverse')
       return
    end if

    !! Solve the linear system of equations [K]*{x} = {F0,F1,F2,...,Fn}
    !! where F0 is the current external load (gravity, etc.), and Fi,i=1,n
    !! are unit load vectors corresponding to the unknown applied loads. The
    !! stiffness matrix [K]=sys%Nmat is assumed to be factorized at this point.
    call DCOPY (ndim,F(1),1,Rhs(1),1)
    call DCOPY (ndim*nfrc,0.0_dp,0,Rhs(ndim+1),1)
    do i = 1, nfrc
       Rhs(ndim*i+f_eqs(i)) = 1.0_dp
    end do
#ifdef FT_DEBUG
    call writeObject (reshape(Rhs,(/ndim,1+nfrc/)),lpu,'Inverse::RHS')
#endif
    call solveLinEqSystem (Rhs,1+nfrc,ierr)
    if (ierr < 0) goto 910
#ifdef FT_DEBUG
    call writeObject (reshape(Rhs,(/ndim,1+nfrc/)),lpu,'Inverse::Sol')
#endif

    !! Establish and solve the inverse problem using least squares
    do i = 1, ndis
       do j = 1, nfrc
          c(i,j) = Rhs(ndim*j+x_eqs(i))
       end do
       R(i) = x(i) - Rhs(x_eqs(i))
    end do
#ifdef FT_DEBUG
    call writeObject (c,lpu,'Inverse::c')
    call writeObject (reshape(R(1:ndis),(/ndis,1/)),lpu,'Inverse::R')
#endif
    call DGELS ('N',ndis,nfrc,1,c(1,1),ndis,R(1),ldR,Work(1),size(Work),ierr)
    if (ierr /= 0) goto 900

    !! Calculate the resulting loads
    call DCOPY (ndim,0.0_dp,0,F(1),1)
    do i = 1, nfrc
       F(f_eqs(i)) = R(i)
    end do
#ifdef FT_DEBUG
    call writeObject (reshape(R(1:nfrc),(/nfrc,1/)),lpu,'Inverse::S')
    call writeObject (reshape(F,(/ndim,1/)),lpu,'Inverse::F')
#endif

890 deallocate(Rhs,c,R,Work)
    return

900 write(errMsg,"('Error return from LAPACK::DGELS, INFO =',I6)") ierr
    call reportError (error_p,'Failed to solve the inverse problem', &
         &            addString=errMsg)

910 call reportError (debugFileOnly_p,'coreInverse')
    if (i > 0) goto 890

  contains

    subroutine checkEqns (meqn,ierr)
      integer, intent(in)    :: meqn(:)
      integer, intent(inout) :: ierr
      do i = 1, size(meqn)
         if (meqn(i) < 1 .or. meqn(i) > ndim) then
            write(errMsg,"('Equation no.',I6,' is out of range')") meqn(i)
            ierr = ierr + internalError('coreInverse: '//errMsg)
         end if
      end do
    end subroutine checkEqns

  end subroutine coreInverse

end module inverseModule
