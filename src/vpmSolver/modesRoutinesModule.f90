!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file modesRoutinesModule.f90
!> @brief Subroutines for eigenvalue analysis.

!!==============================================================================
!> @brief Module with subroutines for the system-level eigenvalue analysis.
!>
!> @author Knut Morten Okstad
!> @date 23 Mar 2000

module ModesRoutinesModule

  use KindModule, only : dp

  implicit none

  private

  !> @cond NO_DOCUMENTATION
  !! Internal work arrays for LAPACK eigensolver (undocumented)
  integer , save              :: Lwork = 0, N = 0
  integer , save, allocatable :: indx(:), Iwork(:)
  logical , save, allocatable :: Bwork(:)
  real(dp), save, allocatable :: alphaR(:), alphaI(:), beta(:), Work(:)
  real(dp), save, allocatable :: A(:,:), B(:,:)
  real(dp), save, pointer     :: V(:,:)
  real(dp), save, allocatable :: Lscale(:), Rscale(:), rcondE(:), rcondV(:)
  !> @endcond

  public :: allocEigenMatrices, eigenModes, printModes, exportModes


contains

  !!============================================================================
  !> @brief Allocates internal work arrays used by the eigenvalue solver.
  !>
  !> @param modes Eigenmode data container
  !> @param[in] neq Number of equations, i.e., dimension of system matrices
  !> @param[in] neqEig Number of active equation during the eigenvalue analysis
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph

  subroutine allocEigenMatrices (modes,neq,neqEig,ierr)

    use ModesTypeModule  , only : ModesType
    use reportErrorModule, only : AllocationError, InternalError

    type(ModesType)  , intent(inout) :: modes
    integer, optional, intent(in)    :: neq, neqEig
    integer, optional, intent(out)   :: ierr

    !! Local variables
    integer  :: idum, jerr, maxLan
    real(dp) :: rdum, Nwork

    !! --- Logic section ---

    if (present(ierr)) then
       ierr = 0
    end if

    if (present(neq)) then
       if (modes%solver == 4 .or. modes%solver == 6) then
          idum = modes%solver/2
          if (N == idum*neq) return
          N = idum*neq
       else
          if (N == neq) return
          N = neq
       end if
    else
       N = 0
       Lwork = 0
    end if

    if (allocated(A))      deallocate(A)
    if (allocated(B))      deallocate(B)
    if (allocated(alphaR)) deallocate(alphaR)
    if (allocated(alphaI)) deallocate(alphaI)
    if (allocated(beta))   deallocate(beta)
    if (allocated(indx))   deallocate(indx)
    if (allocated(Lscale)) deallocate(Lscale)
    if (allocated(Rscale)) deallocate(Rscale)
    if (allocated(rcondE)) deallocate(rcondE)
    if (allocated(rcondV)) deallocate(rcondV)
    if (allocated(Iwork))  deallocate(Iwork)
    if (allocated(Bwork))  deallocate(Bwork)
    if (allocated(Work))   deallocate(Work)

    if (associated(V)) then
       if (.not.associated(V,modes%eqVec)) deallocate(V)
       nullify(V)
    end if
    if (associated(modes%eqVec)) then
       deallocate(modes%eqVec)
       nullify(modes%eqVec)
    end if

    if (N < 1) return ! We are deallocating only

    if (modes%solver <= 2) then

       if (present(neqEig)) then
          maxLan = neqEig
       else
          maxLan = neq
       end if
       maxLan = min(maxLan,12+3*modes%nModes) ! Safe upper bound

       allocate(alphaR(N),modes%eqVec(N,maxLan),stat=jerr)
       if (jerr /= 0) goto 915
       return

    else if (modes%solver == 4 .or. modes%solver == 6) then

       allocate(V(N,N),stat=jerr)
       if (jerr /= 0) goto 915

    else

       allocate(modes%eqVec(N,N),stat=jerr)
       if (jerr /= 0) goto 915
       V => modes%eqVec

    end if

    if (modes%solver <= 4 .or. modes%solver == 6) then
       allocate(alphaR(N),alphaI(N),A(N,N),B(N,N),beta(N),indx(N), &
            &   Lscale(N),Rscale(N),rcondE(N),rcondV(N), &
            &   Iwork(N+6),Bwork(N),stat=jerr)
    else
       allocate(alphaR(N),A(N,N),B(N,N),indx(N),Iwork(5*N),stat=jerr)
    end if
    if (jerr /= 0) goto 915

    Lwork = -1 ! Find the optimal size of the scratch array Lwork
    if (modes%solver <= 4 .or. modes%solver == 6) then
       call DGGEVX ('B','N','V','B', &
            &       N,A(1,1),N,B(1,1),N,alphaR(1),alphaI(1),beta(1), &
            &       V(1,1),N,V(1,1),N,idum,idum,Lscale(1),Rscale(1),rdum,rdum, &
            &       rcondE(1),rcondV(1),Nwork,Lwork,Iwork(1),Bwork(1),jerr)
    else
       call DSYGVX (1,'V','I','U', &
            &       N,A(1,1),N,B(1,1),N,rdum,rdum, &
            &       1,modes%nModes,rdum,idum,alphaR(1),V(1,1),N, &
            &       Nwork,Lwork,Iwork(1),indx(1),jerr)
    end if
    if (jerr /= 0) then
       if (present(ierr)) ierr = InternalError('allocEigenMatrices')
       return
    end if

    Lwork = INT(Nwork)
    allocate(Work(Lwork),stat=jerr)
    if (jerr == 0) return

915 continue
    if (present(ierr)) ierr = AllocationError('allocEigenMatrices')

  end subroutine allocEigenMatrices


  !!============================================================================
  !> @brief Builds the full matrices of the generalized eigenvalue problem.
  !>
  !> @param modes Eigenmode data container
  !> @param[in] neq Number of equations, i.e., dimension of system matrices
  !> @param[in] M System mass matrix
  !> @param[in] C System damping matrix
  !> @param[in] K System stiffness matrix
  !> @param[in] Q System steady-state error elimination matrix
  !> @param[out] ierr Error flag
  !>
  !> @details The generalized eigenvalue problem is on the form
  !>          @b A @b z = &lambda; @b B @b z
  !>
  !> Without damping:
  !>
  !>          [K]*{z} = lambda*[M]*{z}
  !>
  !> With damping (2n state-space form):
  !>
  !>          | K  0 |*{z} = lambda * | -C -M |*{z}
  !>          | 0 -M |                | -M  0 |
  !>
  !> With damping and/or steady-state error elimination (3n state-space form):
  !>
  !>          | Q  0  0 |*{z} = lambda * | -K -C -M |*{z}
  !>          | 0  K  0 |                |  K  0  0 |
  !>          | 0  0  M |                |  0  M  0 |
  !>
  !> @note The 3n state-space method does not work quite yet.
  !> No eigenvalues are obtained for unknown reasons.
  !> However, the @b A and @b B matrices seems to be correct.
  !> The method has been proven valid in MATLAB.
  !> But the @b A and @b B matrices obtained through this code yield
  !> the same eigenvalue results in MATLAB.
  !> Can it be due to singular or very sparse matrices?
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !> @date 22 Mar 2000
  !>
  !> @author Magne Bratland
  !> @date 8 Oct 2010

  subroutine fillEigenSystem (modes,neq,M,C,K,Q,ierr)

    use ModesTypeModule    , only : ModesType
    use SysMatrixTypeModule, only : SysMatrixType, convertSysMat
    use reportErrorModule  , only : reportError, debugFileOnly_p

    type(ModesType)    , intent(inout) :: modes
    integer            , intent(in)    :: neq
    type(SysMatrixType), intent(in)    :: M, C, K, Q
    integer            , intent(out)   :: ierr

    !! Local variables
    integer :: i, j, n2, n3

    !! --- Logic section ---

    call allocEigenMatrices (modes,neq,ierr)
    if (ierr < 0) goto 915

    A = 0.0_dp
    B = 0.0_dp
    n2 = 2*neq
    n3 = 3*neq

    !! The A matrix

    if (modes%solver == 6) then
       !! A(1:n,1:n) = Q
       call convertSysMat (Q,A,neq,1,11,ierr)
       if (ierr < 0) goto 915
       !! A(n+1:2n,n+1:2n) = K
       call convertSysMat (K,A(:,neq+1:),neq,neq+1,11,ierr)
       if (ierr < 0) goto 915
       !! A(2n+1:3n,2n+1:3n) = M
       call convertSysMat (M,A(:,n2+1:),neq,n2+1,11,ierr)
       if (ierr < 0) goto 915
    else
       !! A(1:n,1:n) = K
       call convertSysMat (K,A,neq,1,11,ierr)
       if (ierr < 0) goto 915
    end if
    if (modes%solver == 4) then
       !! A(n+1:2n,n+1:2n) = -M
       call convertSysMat (M,A(:,neq+1:),neq,neq+1,-11,ierr)
       if (ierr < 0) goto 915
    end if

    !! The B matrix

    if (modes%solver == 6) then
       !! B(1:n,1:n) = -K
       B(1:neq,1:neq) = -A(neq+1:n2,neq+1:n2)
       !! B(n+1:2n,1:n) = K
       B(neq+1:n2,1:neq) = A(neq+1:n2,neq+1:n2)
       !! B(1:n,n+1:2n) = -C
       call convertSysMat (C,B(:,neq+1:),neq,1,-11,ierr)
       if (ierr < 0) goto 915
       !! B(2n+1:3n,n+1:2n) = M
       B(n2+1:n3,neq+1:n2) = A(n2+1:n3,n2+1:n3)
       !! B(1:n,2n+1:3n) = -M
       B(1:neq,n2+1:n3) = -A(n3+1:n3,n2+1:n3)
    else if (modes%solver == 4) then
       !! B(n+1:2n,1:n) = M
       B(neq+1:n2,1:neq) = A(neq+1:n2,neq+1:n2)
       !! B(1:n,n+1:2n) = M
       B(1:neq,neq+1:n2) = A(neq+1:n2,neq+1:n2)
       !! B(1:n,1:n) = -C
       call convertSysMat (C,B,neq,1,-11,ierr)
       if (ierr < 0) goto 915
    else
       !! B(1:n,1:n) = M
       call convertSysMat (M,B,neq,1,11,ierr)
       if (ierr < 0) goto 915
    end if

    if (modes%solver == 5) then
       !! Symmetrice the dense stiffness matrix, A
       do j = 2, neq
          do i = 1, j-1
             A(i,j) = 0.5_dp*(A(i,j)+A(j,i))
             A(j,i) = A(i,j)
          end do
       end do
    end if

    return

915 call reportError (debugFileOnly_p,'fillEigenSystem')

  end subroutine fillEigenSystem


  !!============================================================================
  !> @brief Sorts a real array in ascending order.
  !>
  !> @param[in] V Array of real values to sort
  !> @param ind Array of indices defining the sorted order
  !>
  !> @callergraph

  subroutine bubbelSort (V,ind)

    real(dp), intent(in)    :: V(:)
    integer , intent(inout) :: ind(:)

    !! Local variables
    integer             :: i, j, k
    real(dp), parameter :: eps = 1.0e-8_dp

    !! --- Logic section ---

    do i = 1, size(ind)
       k = 0
       do j = size(ind), i+1, -1
          if (V(ind(j-1)) > V(ind(j))+V(ind(j))*eps) then
             k = ind(j)
             ind(j) = ind(j-1)
             ind(j-1) = k
          end if
       end do
       if (k == 0) return
    end do

  end subroutine bubbelSort


  !!============================================================================
  !> @brief Extracts the real eigenvalues and associated eigenvectors.
  !>
  !> @param[in] sam Data for managing system matrix assembly
  !> @param modes Eigenmode data
  !> @param[out] normFactor Scaling factors for normalizing the eigenmode shapes
  !> @param[out] modeOrder Order of the modes with respect to eigenfrequency
  !> @param[in] iprint Print switch for additional output
  !> @param[in] io File unit number to write to
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph

  subroutine findRealEigenValues (sam,modes,normFactor,modeOrder,iprint,io,ierr)

    use SamModule         , only : SamType
    use ModesTypeModule   , only : ModesType
    use KindModule        , only : epsDiv0_p
    use SolExtensionModule, only : csExpand
    use MatExtensionModule, only : csTransform
    use reportErrorModule , only : reportError, debugFileOnly_p

    type(SamType)  , intent(in)    :: sam
    type(ModesType), intent(inout) :: modes
    real(dp)       , intent(out)   :: normFactor(:)
    integer        , intent(out)   :: modeOrder(:)
    integer        , intent(in)    :: iprint, io
    integer        , intent(out)   :: ierr

    !! Local variables
    integer  :: i, j
    real(dp) :: wmax

    !! --- Logic section ---

    ierr = 0
    wmax = 1.0_dp
    do i = 1, N
       indx(i) = i
       if (.not. allocated(beta)) then
          if (i > modes%nModes .and. modes%solver > 4) then
             work(i:) = 0.0_dp
             indx(i:) = N+1
             exit
          end if
          work(i) = alphaR(i)
       else if (beta(i) <= -epsDiv0_p) then
          write(io,"(' *** Eigenvalue #',I6,' is invalid, beta = ',1PE13.5)") &
               i, beta(i)
          work(i) = 0.0_dp
          indx(i) = N+1
          ierr = ierr - 1
       else if (beta(i) < epsDiv0_p) then
          write(io,"(' *** Eigenvalue #',I6,' is infinite (ignored)')") i
          work(i) = 0.0_dp
          indx(i) = N+1
       else if (abs(alphaI(i)) > epsDiv0_p) then
          write(io,"(' *** Eigenvalue #',I6,' is complex, alphaI =',1PE13.5)") &
               i, alphaI(i)
          work(i) = 0.0_dp
          indx(i) = N+1
          ierr = ierr - 1
       else
          work(i) = alphaR(i)/beta(i)
       end if
       wmax = max(wmax,work(i))
    end do
    work(N+1) = wmax+wmax
    if (ierr < 0) goto 915

    !! Sort the modes in ascending order of eigenfrequency
    call bubbelSort (work,indx)

    do i = 1, modes%nModes
       j = indx(i)
       if (j > N) return

       !! Compute the real eigenvalue
       if (allocated(beta)) then
          modes%ReVal(i) = sqrt(alphaR(j)/beta(j))
       else
          modes%ReVal(i) = sqrt(alphaR(j))
       end if

       !! Orthogonalize the eigenvector with respect to the mass matrix
       call csTransform (modes%Mmat,modes%eqVec(:,j),wmax,ierr)
       if (ierr /= 0) goto 915

       !! Expand from equation-based to nodal-based ordering
       wmax = 1.0_dp/sqrt(wmax)
       normFactor(i) = wmax ! By LIM. Used in ModalMasses
       modeOrder(i)  = j    ! To get the ordering and scaling of the modes right
       call csExpand (sam,modes%eqVec(:,j),modes%ReVec(:,i),wmax,0.0_dp)
       if (iprint > 0 .and. allocated(rcondE)) then
          write(io,6000) i, 1.0_dp, modes%ReVal(i)**2.0_dp, rcondE(j), rcondV(j)
       end if

    end do
    return

915 continue
    call reportError (debugFileOnly_p,'findRealEigenValues')

6000 format(' Mode',I4,' :  M =',1P,E12.5,'  K =',E12.5, &
          & '  cond(E) =',E12.5,'  cond(V) =',E12.5)

  end subroutine findRealEigenValues


  !!============================================================================
  !> @brief Extracts the complex eigenvalues and associated eigenvectors.
  !>
  !> @param[in] sam Data for managing system matrix assembly
  !> @param modes Eigenmode data
  !> @param[in] iprint Print switch for additional output
  !> @param[in] io File unit number to write to
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph

  subroutine findEigenValues (sam,modes,iprint,io,ierr)

    use SamModule         , only : SamType
    use ModesTypeModule   , only : ModesType
    use KindModule        , only : epsDiv0_p
    use SolExtensionModule, only : csExpand
    use MatExtensionModule, only : csTransform
    use reportErrorModule , only : reportError, debugFileOnly_p

    type(SamType)  , intent(in)    :: sam
    type(ModesType), intent(inout) :: modes
    integer        , intent(in)    :: iprint, io
    integer        , intent(out)   :: ierr

    !! Local variables
    integer  :: i, j, k, l, nmod, ndof, dof1, dof2
    real(dp) :: cr, ci, reVal, imVal, wmax

    real(dp), parameter :: overdamped = 10.0_dp

    !! --- Logic section ---

    !! To whom it may concern:
    !! From explanation of DGGEVX routine
    !! Note:   the   quotients   ALPHAR(j)/BETA(j)    and
    !! ALPHAI(j)/BETA(j)  may  easily over- or underflow,
    !! and BETA(j) may even  be  zero.   Thus,  the  user
    !! should   avoid   naively   computing   the   ratio
    !! ALPHA/BETA. However, ALPHAR  and  ALPHAI  will  be
    !! always  less  than  and  usually  comparable  with
    !! norm(A) in magnitude, and BETA  always  less  than
    !! and usually comparable with norm(B).
    !! Magne Bratland, 15 Nov 2010.

    ierr = 0
    wmax = 1.0_dp
    do i = 1, N
       indx(i) = i
       if (beta(i) <= -epsDiv0_p) then
          write(io,"(' *** Eigenvalue #',I6,' is invalid, beta=',1PE13.5)") &
               i, beta(i)
          work(i) = 0.0_dp
          indx(i) = N+1
          ierr = ierr - 1
       else if (beta(i) < epsDiv0_p) then
          write(io,"(' *** Eigenvalue #',I6,' is infinite (ignored)')") i
          work(i) = 0.0_dp
          indx(i) = N+1
       else if (abs(alphaI(i)) < epsDiv0_p) then
          work(i) = alphaR(i)/beta(i)
       else
          work(i) = abs(alphaI(i)/beta(i))
       end if
       wmax = max(wmax,work(i))
    end do
    work(N+1) = wmax+wmax

    !! Sort the modes in ascending order of eigenfrequency
    call bubbelSort (work,indx)

    if (modes%solver == 6) then
       !! Using the 3n state-space solution.
       !! This method will yield at least one purely real eigenvalue
       !! in addition to the solution of the 2n state-space method.
       !! The purely real eigenvalue can be discarded.
       !! Should there be a note about when the real value are positive,
       !! that the system is unstable, perhaps?
       ndof = N/3
       dof1 = ndof+1
       dof2 = ndof*2
    else
       ndof = N/2
       dof1 = 1
       dof2 = ndof
    end if

    k = 0
    nmod = 0
    do i = 1, N
       j = indx(i)
       k = k + 1
       modes%nModes = nmod
       if (j > N) return

       !! Compute the complex eigenvalue
       reVal = alphaR(j)/beta(j)
       imVal = alphaI(j)/beta(j)

       !! Orthogonalize the eigenvectors with respect to the mass matrix
       !! and expand from equation-based to nodal-based ordering

       if (abs(alphaI(j)) < epsDiv0_p) then

          !! Real eigenvalue
          if (nmod == size(modes%ReVal)) return

          if (modes%solver == 6) then
             write(io,"('  ** Ignoring purely real mode',2I6,1P3E13.5)") &
                  &      i, j, reVal, rcondE(j), rcondV(j)
          else if (reVal < -overdamped) then
             write(io,"('  ** Ignoring overdamped real mode',2I6,1P3E13.5)") &
                  &      i, j, reVal, rcondE(j), rcondV(j)
          else
             nmod = nmod + 1
             modes%ReVal(nmod) = reVal
             modes%ImVal(nmod) = imVal
             call csTransform (modes%Mmat,V(1:ndof,j),cr,ierr)
             if (ierr /= 0) goto 915
             wmax = 1.0_dp/sqrt(cr)
             call csExpand (sam,V(1:ndof,j),modes%ReVec(:,nmod),wmax,0.0_dp)
             modes%ImVec(:,nmod) = 0.0_dp
             if (iprint > 0) then
                write(io,6000) nmod, 1.0_dp, &
                     &        -reVal*2.0_dp, reVal*reVal, &
                     &         rcondE(j), rcondV(j)
#ifdef FT_DEBUG
             else
                write(io,"('   * Real mode',20X,2I6,1P3E13.5)") &
                     &      i, j, reVal, rcondE(j), rcondV(j)
#endif
             end if
          end if
          k = 0

       else if (mod(k,2) == 1) then

          !! Complex-conjugate eigenvalues
          if (nmod == size(modes%ReVal)) return

          l = indx(i+1)
          nmod = nmod + 1
          modes%ReVal(nmod) = reVal
          modes%ImVal(nmod) = imVal
          call csTransform (modes%Mmat,V(dof1:dof2,j),cr,ierr)
          if (ierr /= 0) goto 915
          call csTransform (modes%Mmat,V(dof1:dof2,l),ci,ierr)
          if (ierr /= 0) goto 915
          wmax = 1.0_dp/sqrt(cr+ci)
          call csExpand (sam,V(dof1:dof2,j),modes%ReVec(:,nmod),wmax,0.0_dp)
          call csExpand (sam,V(dof1:dof2,l),modes%ImVec(:,nmod),wmax,0.0_dp)
          if (iprint > 0) then
             write(io,6000) nmod, 1.0_dp, &
                  &        -reVal*2.0_dp, reVal*reVal + imVal*imVal, &
                  &         rcondE(j), rcondV(j)
#ifdef FT_DEBUG
          else
             write(io,"('   * Complex mode',17X,2I6,1P2E13.5,'*i ',2E13.5)") &
                  &      i, j, reVal, imVal, rcondE(j), rcondV(j)
#endif
          end if

       else if (modes%ImVal(nmod) < 0.0_dp) then

          !! Store only the complex mode with positive imaginary eigenvalue
          modes%ReVal(nmod) = reVal
          modes%ImVal(nmod) = imVal
          modes%ImVec(:,nmod) = -modes%ImVec(:,nmod)
       end if

    end do
    return

915 continue
    call reportError (debugFileOnly_p,'findEigenValues')

6000 format(' Mode',I4,' :  M =',1P,E12.5,'  C =',E12.5,'  K =',E12.5, &
          & '  cond(E) =',E12.5,'  cond(V) =',E12.5)

  end subroutine findEigenValues


  !!============================================================================
  !> @brief Prints the eigenvalues and eigenvectors to unit IO.
  !>
  !> @param[in] sam Data for managing system matrix assembly
  !> @param modes Eigenmode data
  !> @param[in] time Current simulation time
  !> @param[in] iprint Print switch for additional output
  !> @param[in] io File unit number to write to
  !>
  !> @callergraph

  subroutine printModes (sam,modes,time,iprint,io)

    use KindModule     , only : pi_p
    use SamModule      , only : SamType
    use ModesTypeModule, only : ModesType

    type(SamType)  , intent(in) :: sam
    type(ModesType), intent(in) :: modes
    real(dp)       , intent(in) :: time
    integer        , intent(in) :: iprint, IO

    !! Local variables
    integer  :: i, j, k, Nprint
    real(dp) :: scale, effMass(6)
    real(dp) :: hypot, phase, x, y

    !! Statement functions
    hypot(x,y) = SQRT(x*x+y*y)
    phase(x,y) = MOD(360.0_dp+ATAN2(y,x)*180.0_dp/pi_p,360.0_dp)

    !! --- Logic section ---

    !! Write heading
    write(IO,1000) time
    select case (modes%solver)
    case (0) ; write(IO,1001) 'SAM::SPRLAN'
    case (1) ; write(IO,1001) 'SAM::LANCZ1'
    case (2) ; write(IO,1001) 'SAM::LANCZ2'
    case(3:4); write(IO,1001) 'LAPACK::DGGEVX'
    case (5) ; write(IO,1001) 'LAPACK::DSYGVX'
    case (6) ; write(IO,1001) 'LAPACK::DGGEVX'
    end select

    !! Write eigenvalues
    scale = 0.5_dp / pi_p
    if (modes%solver == 4 .or. modes%solver == 6) then
       write(IO,1002) (I,I=1,5)
       write(IO,1003) (modes%ReVal(I),modes%ImVal(I)*scale,I=1,modes%nModes)
    else
       write(IO,1004) (I,I=1,10)
       write(IO,1005) (modes%ReVal(I)*scale,I=1,modes%nModes)
    end if

    if (associated(modes%effMass)) then

       !! Write effective modal masses
       effMass = 0.0_dp
       write(IO,1100)
       write(IO,1101)
       do K = 1, modes%nModes
          write(IO,1102) K, modes%ReVal(K)*scale, modes%effMass(K,:)
          effMass = effMass + modes%effMass(K,:)
       end do
       write(IO,1103) effMass
    end if

    if (iprint > 9) then

       !! Write eigenvectors
       if (iprint < 15) then
          Nprint = 1
       else
          Nprint = modes%nModes
       end if

       do K = 1, Nprint
          if (modes%solver == 4 .or. modes%solver == 6) then
             write(IO,2001) K
          else
             write(IO,2002) K
          end if
          do I = 1, sam%nnod
             if (modes%solver == 4 .or. modes%solver == 6) then
                !! Complex eigenvectors are written as magnitude and phase
                write(IO,2003) I,(hypot(modes%ReVec(J,K),modes%ImVec(J,K)), &
                     &            phase(modes%ReVec(J,K),modes%ImVec(J,K)), &
                     &            J = sam%madof(I),sam%madof(I+1)-1)
             else
                write(IO,2004) I,(modes%ReVec(J,K), &
                     &            J = sam%madof(I),sam%madof(I+1)-1)
             end if
          end do
       end do

    end if

    write(IO,2005)

1000 format(/130('*')/3X,'EIGENVALUE ANALYSIS AT TIME (S):',1PE12.5)
1001 format(130('*')///3X,'EIGENFREQUENCIES [Hz] (using ',A,') :')
1002 format(/5(3X,'LAMBDA & EIGFREQ',I2,4X),' ....'/)
1003 format(1P,5(1X,2E11.4,'  '))
1004 format(/10(2X,'EIGFREQ',I2,2X)/)
1005 format(1P,10E13.5)

1100 format(//3X,'EFFECTIVE MASSES/INERTIAS FOR THE CALCULATED EIGENVALUES:'/)
1101 format(3X,'Mode #',5X, 'EigFreq [Hz]', &
          & 9X,'X',12X,'Y',12X,'Z',12X,'Rx',11X,'Ry',11X,'Rz'/)
1102 format(1P,4X,I3,5X,E13.5,2X,'|',6E13.5)
1103 format(/'TOTAL EFFECTIVE MASSES:',5X,6E13.5/)

2001 format(//3X,'EIGENVECTOR :',I2 // ' Node No', &
          & 6X,'F1',19X,'F2',19X,'F3',19X,'F4',19X,'F5',19X,'F6'/)
2002 format(//3X,'EIGENVECTOR :',I2 // '    Node No',&
          & 9X,'F1',13X,'F2',13X,'F3',13X,'F4',13X,'F5',13X,'F6'/)
2003 format(I5,6(1PE15.7,0PF6.1))
2004 format(I8,3X,1P,8E15.7)
2005 format(130('*'))

  end subroutine printModes


  !!============================================================================
  !> @brief Calculates damped or undamped eigenvalues and eigenvectors.
  !>
  !> @param[in] sam Data for managing system matrix assembly
  !> @param modes Eigenmode data
  !> @param sys System level model data
  !> @param mech Mechanism components of the model
  !> @param[in] ctrl Control system data
  !> @param pCS Data for coupled control system and structure modal analysis
  !> @param[in] iprint Print switch for additional output
  !> @param[out] ierr Error flag
  !>
  !> @details This is the main driver of the system-level eigenvalue analysis.
  !> The parameter modestypemodule::modestype::solver determines which actual
  !> subroutine to be used, as follows:
  !> - 0 : Use SPRLAN (undamped systems only)
  !> - 1 : Use LANCZ1 (undamped systems only)
  !> - 2 : Use LANCZ2 (undamped systems only)
  !> - 3 : Use LAPACK::DGGEVX (undamped system)
  !> - 4 : Use LAPACK::DGGEVX (damped system, 2n state-space form)
  !> - 5 : Use LAPACK::DSYGVX (undamped, symmetric system)
  !> - 6 : Use LAPACK::DGGEVX (3n state-space form)
  !>
  !> @callgraph @callergraph

  subroutine eigenModes (sam,modes,sys,mech,ctrl,pCS,iprint,ierr)

    use sprKindModule      , only : ik
    use SamModule          , only : SamType
    use ModesTypeModule    , only : ModesType
    use SystemTypeModule   , only : SystemType
    use MechanismTypeModule, only : MechanismType
    use ControlTypeModule  , only : ControlType
    use ControlStructModule, only : ControlStructType, writeObject
    use ControlStructModule, only : BuildStructControlJacobi
    use TriadTypeModule    , only : hasLocalDirections, transSysToGlob
    use AsmExtensionModule , only : csAddEM
    use SolExtensionModule , only : csLanczosEigenSolve, csExpand
    use MatExtensionModule , only : csTransform, csCopyMat
    use AddInSysModule     , only : BuildStiffMat, BuildDamperMat, BuildMassMat
    use profilerModule     , only : startTimer, stopTimer, eig_p
    use reportErrorModule  , only : allocationError, getErrorFile
    use reportErrorModule  , only : reportError, debugFileOnly_p, warning_p
    use reportErrorModule  , only : errorFileOnly_p, warning_p, note_p, empty_p
#ifdef FT_DEBUG
    use SysMatrixTypeModule, only : writeObject
    use manipMatrixModule  , only : writeObject
    use fileUtilitiesModule, only : getDBGfile
#endif

    type(SamType)       , intent(in)    :: sam
    type(ModesType)     , intent(inout) :: modes
    type(SystemType)    , intent(inout) :: sys
    type(MechanismType) , intent(inout) :: mech
    type(ControlType)   , intent(inout) :: ctrl
    type(ControlStructType), pointer    :: pCS
    integer             , intent(in)    :: iprint
    integer             , intent(out)   :: ierr

    !! Local variables
    logical               :: complexModes, mIsDestroyed
    integer               :: i, idof, ldof, io, j, mop(10), neqEig
    integer(ik)           :: meqErr(2)
    real(dp)              :: Abnorm, Bbnorm
    integer , allocatable :: modeOrder(:)
    real(dp), allocatable :: normFactors(:)
    character(len=64)     :: emsg
    integer, parameter    :: iter = -1, zeroStressStiffUpdateSkip = 0
    integer, save         :: mip(7) = (/ 1, 1, -1, 2, 2, 2, 0 /)
    logical, save         :: isFirstNote = .true.
    real(dp), external    :: DLAMCH

    !! --- Logic section ---

    complexModes = modes%solver == 4 .or. modes%solver == 6

    call BuildStiffMat (modes%Kmat,mech,sam,iter,modes%stressStiffIsOn, &
         &              zeroStressStiffUpdateSkip,ierr)
    if (ierr /= 0) goto 915

    call BuildDamperMat (modes%Cmat,mech,sam,ierr)
    if (ierr /= 0) goto 915

    call BuildMassMat (modes%Mmat,mech,sam,ierr)
    if (ierr /= 0) goto 915

    if (associated(pCS)) then
       !! Call routine for df/dx for controller
       call BuildStructControlJacobi (pCS,ctrl,sys,sam%mpar,ierr)
       if (ierr < 0) goto 915
       !! Note: If the controller contains non-collocated sensors and actuators,
       !! the matrices will be unsymmetric. However, these will by the existing
       !! code be faulty symmetrizied
       if (modes%solver == 6) then
          if (IAND(ierr,1) == 0) then
             call reportError (warning_p,'The Q-matrix is identically zero '// &
                  &            'since there are no position controllers.', &
                  &            'Maybe use the 2n state-space solver instead?')
          end if
          call csAddEM (sam,pCS%samElNum,pCS%nDofs,pCS%SSEEMat,modes%Qmat,ierr)
          if (ierr /= 0) goto 915
       end if
       call csAddEM (sam,pCS%samElNum,pCS%nDofs,pCS%stiffMat,modes%Kmat,ierr)
       if (ierr /= 0) goto 915
       call csAddEM (sam,pCS%samElNum,pCS%nDofs,pCS%dampMat,modes%Cmat,ierr)
       if (ierr /= 0) goto 915
       call csAddEM (sam,pCS%samElNum,pCS%nDofs,pCS%massMat,modes%Mmat,ierr)
       if (ierr /= 0) goto 915
    end if

    !! Make a copy of the mass matrix since csLanczosEigenSolve may destroy it
    mIsDestroyed = modes%solver==1 .or. (modes%solver<3 .and. modes%factorMass)
    if (mIsDestroyed) then
       call csCopyMat (sys%Nmat,modes%Mmat,ierr)
       if (ierr /= 0) goto 915
    end if

#ifdef FT_DEBUG
    if (mod(iprint,100) == -99) then
       io = getDBGfile(3,'modes.dbg')
       if (iprint == -199) then
          j = 3
       else
          j = 0
       end if
       write(io,"(/'=== Eigenvalue analysis at time =',1pe12.5)") sys%time
       write(io,"('TolEigval, TolFactor, TolEigvec =',1p3e12.5)") modes%tol
       if (associated(pCS)) call writeObject (pCS,io,5)
       call writeObject (modes%Qmat,io,'Qmat in Eigenvalue analysis',12)
       call writeObject (modes%Kmat,io,'Kmat in Eigenvalue analysis',12,j)
       call writeObject (modes%Cmat,io,'Cmat in Eigenvalue analysis',12,j)
       call writeObject (modes%Mmat,io,'Mmat in Eigenvalue analysis',12,j)
    end if
#endif

    !! neqEig = sam%ndof1 = sam%neq - sam%ndof2  if "-addBC_eigenSolver" is set.
    !! sam%ndof2 is the number of degrees of freedom that should be fixed during
    !! eigenvalue calculations only. These equations are always ordered last in
    !! the system matrices such that we may solve for the first ndof1 equations
    !! only, simply by using ndof1 in place of neq in the equation solver calls.

    if (modes%addBC) then
       neqEig = sam%ndof1
    else
       neqEig = sam%neq
    end if

    call startTimer (eig_p)

    if (modes%solver <= 2) then

       call allocEigenMatrices (modes,sam%neq,neqEig,ierr)
       if (ierr /= 0) then
          call stopTimer (eig_p)
          goto 915
       end if

       !! Solve the real symmetric eigenvalue problem using the LANCZOS solver

       if (modes%factorMass) then
          mip(6) = 2
       else
          mip(6) = 1
       end if
       io = getErrorFile()
       call csLanczosEigenSolve (modes%Kmat, modes%Mmat, MIP, MOP, sam%meqn, &
            &                    neqEig, modes%nModes, modes%nModes, &
            &                    modes%tol, modes%shift, alphaR, modes%eqVec, &
            &                    min(1,iprint/2), IO, meqErr, IERR)
       if (ierr < 0) then
          if (ierr >= -5) ierr = ierr - 10*int(meqErr(1))
          if (modes%solver /= 1 .and. neqEig < 50 .and. isFirstNote) then
             isFirstNote = .false.
             if (modes%solver == 2) then
                emsg = 'LANCZ2 failed'
             else
                emsg = 'SPRLAN failed'
             end if
             call reportError (note_p,'The eigenvalue solver '// trim(emsg) // &
                  &            ' for a system with less than 50 unknowns.')
             if (modes%solver == 2) then
                emsg = ' (-lancz1)'
             else
                emsg = ' (-skylinesolver -nosolveropt -lancz1)'
             end if
             call reportError (empty_p,'Try switching to the LANCZ1 solver '// &
                  'instead'//emsg,'Or use the LAPack eigensolver (-undamped)')
          end if
          call stopTimer (eig_p)
          goto 915
       else if (ierr > 0) then
          write(emsg,"(' at time =',1pe12.5)") sys%time
          if (modes%solver > 0) then
             call reportError (warning_p,'Not all eigenvectors were accepted'//&
                  &            emsg)
          else
             call reportError (warning_p,'Negative definite stiffness matrix.',&
                  &            'Detected in eigenvalue analysis '//emsg)
          end if
          ierr = 0
       end if

       do j = 1, modes%nModes
          modes%ReVal(j) = dsqrt(abs(alphaR(j)))
          if (alphaR(j) < 0.0_dp) ierr = ierr + 1
          call csExpand (sam,modes%eqVec(:,j),modes%ReVec(:,j),1.0_dp,0.0_dp)
          if (iprint > 0) then
             write(io,6000) j, 1.0_dp, modes%ReVal(j)**2.0_dp
          end if
       end do
       if (ierr > 0) then
          write(emsg,"(I3,' eigenvalues are negative!')") ierr
          call reportError (warning_p,emsg)
          write(IO,"(2X,10(2X,'OMEGA^2',I2,2X))") (i,i=1,10)
          write(IO,"(2X,10E13.5)") alphaR(1:modes%nModes)
          ierr = 0
       end if

    else

       !! Establish the first-order generalized eigenvalue problem
       call fillEigenSystem (modes,neqEig,modes%Mmat, &
            &                modes%Cmat,modes%Kmat,modes%Qmat,ierr)
       if (ierr /= 0) then
          call stopTimer (eig_p)
          goto 915
       end if
#ifdef FT_DEBUG
       if (mod(iprint,100) == -99 .and. size(A,1) <= 10) then
          call writeObject (A,io,'Resulting A-matrix',10)
          call writeObject (B,io,'Resulting B-matrix',10)
       end if
#endif

    end if
    if (modes%solver == 5) then

       !! Solve the generalized symmetric eigenvalue problem using LAPACK
       Abnorm = 2.0_dp*DLAMCH('S')
       call DSYGVX (1, 'V', 'I', 'U', &
            &       N, A(1,1), N, B(1,1), N, Work(1), Work(1), &
            &       1, modes%nModes, Abnorm, j, alphaR(1), V(1,1), N, &
            &       Work(1), Lwork, Iwork(1), indx(1), ierr)
       if (ierr /= 0) then
          write(emsg,"('Error from LAPACK::DSYGVX, INFO =',I6)") ierr
          call reportError (errorFileOnly_p,emsg)
          ierr = -2
          call stopTimer (eig_p)
          goto 915
       end if

    else if (modes%solver >= 3) then

       !! Solve the generalized non-symmetric eigenvalue problem using LAPACK
       call DGGEVX ('B', 'N', 'V', 'B', &
            &       N, A(1,1), N, B(1,1), N, alphaR(1), alphaI(1), beta(1), &
            &       V(1,1), N, V(1,1), N, i, j, Lscale(1), Rscale(1), &
            &       Abnorm, Bbnorm, rcondE(1), rcondV(1), &
            &       Work(1), Lwork, Iwork(1), Bwork(1), ierr)
       if (ierr /= 0) then
          write(emsg,"('Error from LAPACK::DGGEVX, INFO =',I6)") ierr
          call reportError (errorFileOnly_p,emsg)
          ierr = -2
          call stopTimer (eig_p)
          goto 915
       end if
#ifdef FT_DEBUG
       if (mod(iprint,100) == -99) then
          write(io,600)
600       format(/'Mode   alphaR       alphaI         beta', &
               &  '       Re(eigval)   Im(eigval)')
          do j = 1, N
             if (beta(j) > 1.0e-16_dp) then
                write(io,"(I4,1P5E13.5)") j, alphaR(j), alphaI(j), beta(j), &
                     &                    alphaR(j)/beta(j), alphaI(j)/beta(j)
             else
                write(io,"(I4,1P3E13.5)") j, alphaR(j), alphaI(j), beta(j)
             end if
          end do
       end if
#endif

    end if
    if (modes%solver >= 3) then
       io = getErrorFile()
       if (complexModes) then

          !! Extract complex eigenvalues and expand the associated eigenvectors
          call findEigenValues (sam,modes,iprint,io,ierr)
          if (modes%nModes < size(modes%ReVal)) then
             modes%ReVal(modes%nModes+1:)   = 0.0_dp
             modes%ImVal(modes%nModes+1:)   = 0.0_dp
             modes%ReVec(:,modes%nModes+1:) = 0.0_dp
             modes%ImVec(:,modes%nModes+1:) = 0.0_dp
          end if

       else

          allocate(modeOrder(modes%nModes),normFactors(modes%nModes),stat=ierr)
          if (ierr /= 0) then
             ierr = allocationError('eigenModes')
             call stopTimer (eig_p)
             return
          end if

          !! Extract real eigenvalues and expand the associated eigenvectors
          call findRealEigenValues (sam,modes,normFactors,modeOrder, &
               &                    iprint,io,ierr)

       end if
    end if

    call stopTimer (eig_p)
    if (ierr /= 0) goto 915

    !! Compute the damping ratios for each eigenmode
    do j = 1, modes%nModes
       if (complexModes) then
          modes%dampRat(j) = -modes%ReVal(j)*2.0_dp
       else if (modes%solver >= 3) then
          call csTransform (modes%Cmat,V(:,modeOrder(j)),modes%dampRat(j),ierr)
          if (ierr < 0) goto 915
          modes%dampRat(j) = 0.5_dp*modes%dampRat(j)*normFactors(j) &
               &           / modes%ReVal(j)
       else
          call csTransform (modes%Cmat,modes%eqVec(:,j),modes%dampRat(j),ierr)
          if (ierr < 0) goto 915
          modes%dampRat(j) = 0.5_dp*modes%dampRat(j) / modes%ReVal(j)
       end if
    end do

    !! Transform all eigenvectors to global directions
    do i = 1, size(mech%triads)
       if (mech%triads(i)%nDOFs == 0) cycle
       if (.not. hasLocalDirections(mech%triads(i))) cycle

       idof = mech%triads(i)%sysDOF
       ldof = idof + mech%triads(i)%nDOFs-1
       do j = 1, modes%nModes
          call transSysToGlob (mech%triads(i),modes%ReVec(idof:ldof,j))
          if (complexModes) then
             call transSysToGlob (mech%triads(i),modes%ImVec(idof:ldof,j))
          end if
       end do

    end do

    if (associated(modes%effMass)) then
       !! Restore the unfactored mass matrix, if needed
       if (mIsDestroyed) then
          call csCopyMat (modes%Mmat,sys%Nmat,ierr)
          if (ierr /= 0) goto 915
       end if
       !! Calculate effective masses for each mode
       call modalMasses (mech%sups,modes,sam,V,normFactors,modeOrder,ierr)
       if (ierr /= 0) goto 915
    end if
    if (allocated(normFactors)) deallocate(normFactors)
    if (allocated(modeOrder))   deallocate(modeOrder)

    return

915 continue
    call reportError (debugFileOnly_p,'eigenModes')

6000 format(' Mode',I4,' :  M =',1P,E12.5,'  K =',E12.5)

  end subroutine eigenModes


  !!============================================================================
  !> @brief Calculates effective modal masses for the eigenmodes.
  !>
  !> @param[in] sups All superelements in the model
  !> @param modes Eigenmode data
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[in] eigVectors The eigenvectors
  !> @param[in] normFactors Scaling factors for normalizing the eigenmode shapes
  !> @param[in] modeOrder Order of the modes with respect to eigenfrequency
  !> @param[out] ierr Error flag
  !>
  !> @details This subroutine is used to judge the significance of each
  !> vibration mode for each translation- and rotation component.
  !> Modes with relatively high effective masses can be readily excited by base
  !> exitation in the actual direction, while modes with lower effective masses
  !> may not be that critical. This is a commonly used feature, e.g. in NASTRAN.
  !>
  !> The calculation is as follows:
  !>
  !>   [temp] = [&Phi;]^t[M]*[R]
  !>
  !>   M_eff(i,j) = temp(i,j)^2
  !>
  !> where the eigenvector matrix, [&Phi;], is transposed and assumed to be
  !> mass-normalized, [M] is the compact mass matrix, and [R] is the influence
  !> matrix representing a rigid body unit displacement of the mechanism.
  !> Finally, [M_eff] is a nModes&times;6 matrix of the effective modal masses.
  !> The results are stored in the array modestypemodule::modestype::effmass.
  !>
  !> @todo It works now for LANCZ2 and for the -undamped solvers (although I
  !> encountered something that looks like a bug in the -undamped calculation of
  !> the eigenvectors). Furthermore, local orientation of triads are not (yet)
  !> taken into account, but should normally not be critical as long as the
  !> local triad masses are the same in all coordinate directions (for inertias
  !> it will induce errors, however). Component modes seems to be treated
  !> correctly, but has not been tested too much. Some more elegance could also
  !> be put into the creation of unit displacement vectors, just to make sure
  !> that all kinds of joints and boundary conditions are treated correctly.
  !> In fact, there might be a better and easier solution to this than
  !> the approach I have used here.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Leif Ivar Myklebust
  !>
  !> @date 2 Jul 2003

  subroutine modalMasses (sups,modes,sam,eigVectors,normFactors,modeOrder,ierr)

    use SamModule         , only : SamType
    use ModesTypeModule   , only : ModesType
    use SupElTypeModule   , only : SupElType
    use MatExtensionModule, only : csPremult
    use reportErrorModule , only : allocationError, reportError, debugFileOnly_p

    type(SupElType), intent(in)    :: sups(:)
    type(ModesType), intent(inout) :: modes
    type(SamType)  , intent(in)    :: sam
    real(dp)       , intent(in)    :: eigVectors(:,:), normFactors(:)
    integer        , intent(in)    :: modeOrder(:)
    integer        , intent(out)   :: ierr

    !! Local variables
    integer               :: i, imod
    real(dp), allocatable :: ttcc_0(:), R_mat(:,:), eig(:,:), modal_part(:,:)

    !! --- Logic Section ---

    allocate(ttcc_0(sam%nmmceq), R_mat(sam%neq,6), &
         &   eig(sam%neq,modes%nModes), modal_part(modes%nModes,6), stat=ierr)
    if (ierr /= 0) then
       ierr = allocationError('modalMasses')
       return
    end if

    R_mat = 0.0_dp

    !! Check if there are any joints that are not attached to ground.
    !! If so, the factors in TTCC should be set to zero in order to keep
    !! the joint rigid while establishing the unit displacement vectors.
    ttcc_0 = sam%ttcc
    do i = 1, sam%nceq
       if (sam%mpmceq(i+1)-sam%mpmceq(i) > 7) then
          !! The joint is not connected to ground
          ttcc_0(sam%mpmceq(i)+1:sam%mpmceq(i+1)-1) = 0.0_dp
       end if
    end do

    !! (1) Keep the eigenvectors on equation form, but scale them
    !! (for modes%solver = 3) to keep things mass-normalized
    if (modes%solver >= 3) then
       do imod = 1, modes%nModes
          eig(:,imod) = eigVectors(:,modeOrder(imod)) * normFactors(imod)
       end do
    else
       eig = eigVectors(:,1:modes%nModes)
    end if

    !! (2) Establish the influence vectors in all 6 directions
    !! These are the vectors that represent the displacement of the masses
    !! resulting from static application of a unit ground displacement
    do i = 1, size(sups)
       call addUnitDisplVector (sam, sups(i), ttcc_0, R_mat, ierr)
       if (ierr /= 0) goto 915
    end do

    !! (3) Multiply the mass matrix with the influence matrix
    call csPremult (modes%Mmat, R_mat, ierr)
    if (ierr /= 0) goto 915

    !! (4) Establish the participation matrix by premultiplying
    !! with the transposed eigenvectors
    modal_part = matmul(transpose(eig),R_mat)

    !! (5) Square the participation matrix, and voila!
    !! If the eigenvectors were not mass normalized, we would here divide each
    !! term with the corresponding generalized mass.
    do imod = 1, modes%nModes
       do i = 1, 6
          modes%EffMass(imod,i) = modal_part(imod,i)*modal_part(imod,i)
       end do
    end do

900 deallocate(R_mat, ttcc_0, eig, modal_part)
    return

915 call reportError (debugFileOnly_p,'modalMasses')
    goto 900

  end subroutine modalMasses


  !!============================================================================
  !> @brief Adds a superelement unit displacement vector into the system vector.
  !>
  !> @param[in] samData Data for managing system matrix assembly
  !> @param[in] sup The superelement to add unit displacement vector for
  !> @param[in] ttcc_0 Modified table of constraint equation coefficients
  !> @param sysV System displacement vector
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Leif Ivar Myklebust
  !>
  !> @date 2 Jul 2003

  subroutine addUnitDisplVector (samData,sup,ttcc_0,sysV,ierr)

    use SamModule                 , only : SamType
    use SupElTypeModule           , only : SupElType
    use MassMatrixCorrectionModule, only : mmcRigAccelVectors
    use reportErrorModule         , only : allocationError, getErrorFile
    use reportErrorModule         , only : reportError, debugFileOnly_p

    type(SamType)  , intent(in)    :: samData
    type(SupElType), intent(in)    :: sup
    real(dp)       , intent(in)    :: ttcc_0(:)
    real(dp)       , intent(inout) :: sysV(:,:)
    integer        , intent(out)   :: ierr

    !! Local variables
    real(dp), parameter   :: eps_p = 1.0e-16_dp
    real(dp), parameter   :: Center(3) = (/ 0.0_dp, 0.0_dp, 0.0_dp /)
    integer               :: i, k
    integer , allocatable :: meen(:)
    real(dp), allocatable :: R(:,:), sysV_temp(:)

    !! --- Logic section ---

    allocate(meen(sup%nTotDofs), R(sup%nTotDofs,6), &
         &   sysV_temp(samdata%neq), stat=ierr)
    if (ierr /= 0) then
       ierr = allocationError('addUnitDisplVector')
       return
    end if

    !! Establish rigid body acceleration vectors about global origin
    call mmcRigAccelVectors (sup, Center, R)

    do k = 1, 6 ! Loop through all directions TX TY TZ RX RY RZ
       sysV_temp = sysV(:,k)
       call ADDEV (R(1,k)          , ttcc_0(1)        , samData%mpar(1)  , &
            &      samData%madof(1), samData%meqn(1)  , samData%mpmnpc(1), &
            &      samData%mmnpc(1), samData%mpmceq(1), samData%mmceq(1) , &
            &      sup%samElNum    , sup%nTotDofs     , getErrorFile()   , &
            &      sysV_temp(1)    , meen(1)          , ierr)
       if (ierr /= 0) exit

       !! To avoid double counting of dofs
       do i = 1, samdata%neq
          if (abs(sysV_temp(i)) > eps_p .and. abs(sysV(i,k)) <= eps_p) then
             sysV(i,k) = sysV_temp(i)
          end if
       end do
    end do

    deallocate(meen, R, sysV_temp)

    if (ierr /= 0) call reportError (debugFileOnly_p,'addUnitDisplVector')

  end subroutine addUnitDisplVector


  !!============================================================================
  !> @brief Exports mode shapes with associated beam section forces to YAML.
  !>
  !> @param[in] yamlFile Name of the YAML-file top write
  !> @param[in] modelFile Name of the model file
  !> @param[in] modes Eigenmode data
  !> @param triads All triads in the model
  !> @param sups All superelements in the model
  !> @param[in] iStep Time increment counter
  !> @param[in] time Current simulation time
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 21 Jun 2017

  subroutine exportModes (yamlFile,modelFile,modes,triads,sups,iStep,time,ierr)

    use KindModule         , only : i8
    use ModesTypeModule    , only : ModesType
    use TriadTypeModule    , only : TriadType, IncTriadsPos
    use SupElTypeModule    , only : SupElType, IsBeam
    use fileUtilitiesModule, only : findUnitNumber
#ifdef FT_DEBUG
    use dbgUnitsModule     , only : dbgSolve
#endif
    use kindModule         , only : lfnam_p
    use reportErrorModule  , only : reportError
    use reportErrorModule  , only : note_p, error_p, debugFileOnly_p

    character(len=*), intent(in)    :: yamlFile, modelFile
    type(ModesType) , intent(in)    :: modes
    type(TriadType) , intent(inout) :: triads(:)
    type(SupElType) , intent(inout) :: sups(:)
    integer(i8)     , intent(in)    :: iStep
    real(dp)        , intent(in)    :: time
    integer         , intent(out)   :: ierr

    !! Local variables
    integer  :: i, j, iDbg, iYaml, iYlen
    real(dp) :: scale

    character(len=lfnam_p) :: fileName

    !! --- Logic section ---

    ierr = 0
    if (yamlFile == '') return

    iDbg = 0
#ifdef FT_DEBUG
    if (iStep < 1000_i8) then
       iDbg = dbgSolve
       write(iDbg,"(/'   --- in exportModes, step,t =',I3,1PE12.5)") iStep,time
    end if
#endif

    if (iStep > 0_i8) then
       write(fileName,"(I20)") iStep
       fileName = trim(yamlFile)//'_'//trim(adjustl(fileName))//'.yaml'
    else
       fileName = trim(yamlFile)//'.yaml'
    end if

    iYaml = findUnitNumber(10)
    iYlen = len_trim(fileName)
    open(iYaml,FILE=fileName(1:iYlen),STATUS='UNKNOWN',IOSTAT=ierr)
    if (ierr /= 0) then
       call reportError (error_p,'Unable to open output file '//fileName)
       return
    else
       call reportError (note_p,'Exporting mode shapes to file '//fileName)
    end if

    call writeYAMLheader (iYaml,'fedem_solver',modelFile,iStep,time)
    do i = 1, modes%nModes
       if (iDbg > 0) write(iDbg,"(/5X,'- Mode',I3,':')") i
       scale = getScale(modes%ReVec(:,i))
       call writeYAMLdispl (iYaml,triads,i,modes%ReVal(i),modes%ReVec(:,i))
       call incTriadsPos (triads,modes%ReVec(:,i)*scale)
       call writeYAMLforces (iYaml,sups,scale,iDbg,ierr)
       do j = 1, size(triads)
          if (triads(j)%nDOFs > 0) triads(j)%ur = triads(j)%urPrev
       end do
       do j = 1, size(sups)
          if (IsBeam(sups(j))) sups(j)%supTr = sups(j)%supTrPrev
       end do
       if (ierr < 0) exit
    end do

    close(iYaml)
    if (ierr < 0) call reportError (debugFileOnly_p,'exportModes')

  contains

    !> @brief Calculates scaling factor for an eigenvector.
    function getScale (eigVec) result(scale)
      real(dp), intent(in) :: eigVec(:)
      real(dp), parameter  :: maxRot_p = 0.1_dp

      real(dp) :: maxRot, scale
      integer  :: i, j, k, idof

      !! Find the largest angular component in the eigenVector
      maxRot = 0.0_dp
      do i = 1, size(sups)
         if (IsBeam(sups(i))) then
            do j = 1, sups(i)%nExtNods
               idof = sups(i)%triads(j)%p%sysDOF + 3
               do k = 4, sups(i)%triads(j)%p%nDOFs
                  if (abs(eigVec(idof)) > maxRot) maxRot = abs(eigVec(idof))
                  idof = idof + 1
               end do
            end do
         end if
      end do

      !! Scale this eigenvector such that no rotational components are larger
      !! than maxRot_p, such that we are safely within the linear regime when
      !! calculating the internal forces
      if (maxRot > maxRot_p) then
         scale = maxRot_p/maxRot
      else
         scale = 1.0_dp
      end if

    end function getScale

  end subroutine exportModes


  !!============================================================================
  !> @brief Writes the header of a YAML frequency response file.
  !>
  !> @param[in] iYaml File unit number of the YAML-file to write
  !> @param[in] prog Program name
  !> @param[in] modelFile Name of the model file
  !> @param[in] iStep Time increment counter
  !> @param[in] time Current simulation time
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 21 Jun 2017

  subroutine writeYAMLheader (iYaml,prog,modelFile,istep,time)

    use kindModule             , only : i8
    use versionModule          , only : progVer, buildDate, getCurrentDate
    use computerConfigInterface, only : getUserName, getComputerConfig

    integer         , intent(in) :: iYaml
    character(len=*), intent(in) :: prog, modelFile
    integer(i8)     , intent(in) :: iStep
    real(dp)        , intent(in) :: time

    !! Local variables
    character(len=16) :: user
    character(len=80) :: chid
    character(len=32) :: date

    !! --- Logic section ---

    call getUserName (user)
    call getComputerConfig (chid)
    call getCurrentDate (date)

    write(iYaml,600) 1.0
    write(iYaml,601) trim(prog), trim(progVer), trim(buildDate)
    write(iYaml,602) trim(modelFile)
    write(iYaml,603) trim(user)
    write(iYaml,604) trim(chid)
    write(iYaml,605) trim(date)
    if (iStep > 99999_i8) then
       write(iYaml,608) time
    else if (iStep > 0_i8) then
       write(iYaml,606) iStep, time
    end if
    write(iYaml,607)

600 format('# YAML formatted input file for frequency response.' &
         //'Format version:',F4.1)
601 format(/'Solver version: ',A,1X,A,' (Build date ',A,')')
602 format(/'Model input file: ',A)
603 format(/'User: ',A)
604 format(/'Computer: ',A)
605 format(/'Date and time: ',A)
606 format(/'Simulation time step:',I6,1PE12.5)
608 format(/'Simulation time:',1PE12.5)
607 format(/'Modes:')

  end subroutine writeYAMLheader


  !!============================================================================
  !> @brief Writes the mode shape eigenvector to a YAML frequency response file.
  !>
  !> @param[in] iYaml File unit number of the YAML-file to write
  !> @param[in] triads All triads in the model
  !> @param[in] iMod Mode shape index
  !> @param[in] eigVal Associated eigenvalue
  !> @param[in] eigVec Associated eigenvector
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 21 Jun 2017

  subroutine writeYAMLdispl (iYaml,triads,iMod,eigVal,eigVec)

    use kindModule     , only : pi_p
    use TriadTypeModule, only : TriadType

    integer        , intent(in) :: iYaml, iMod
    type(TriadType), intent(in) :: triads(:)
    real(dp)       , intent(in) :: eigVal, eigVec(:)

    !! Local variables
    integer :: i, idof

    !! --- Logic section ---

    if (iMod > 1) write(iYaml,"('')")
    write(iYaml,600) iMod, eigVal*0.5_dp/pi_p, &
         &           max(maxval(eigVec),-minval(eigVec))
    do i = 1, size(triads)
       if (triads(i)%nDOFs >= 6) then
          idof = triads(i)%sysDOF
          write(iYaml,601) triads(i)%id%baseId, eigVec(idof:idof+5)
       end if
    end do

600 format('- Name: Mode',I3 &
         /'  Frequency:',1PE12.5,' # [Hz]' &
         /'  Amplitude:',1PE12.5 &
         /'  Shape: # TriadID, dx, dy, dz, rx, ry, rz [m, m, m, rad, rad, rad]')
601 format('  - [',I6,1P,6(',',E12.5),']')

  end subroutine writeYAMLdispl


  !!============================================================================
  !> @brief Writes beam sectional forces to a YAML frequency response file.
  !>
  !> @param[in] iYaml File unit number of the YAML-file to write
  !> @param sups All superelements in the model
  !> @param[in] eigScale Scaling factor for the eigenvector
  !> @param[in] iDbg File unit number for debug output
  !> @param[in] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 21 Jun 2017

  subroutine writeYAMLforces (iYaml,sups,eigScale,iDbg,ierr)

    use SupElTypeModule  , only : SupElType, IsBeam, getSupElId
    use SupElTypeModule  , only : UpdateSupElCorot, BuildFinit
    use manipMatrixModule, only : writeObject
    use reportErrorModule, only : reportError, debugFileOnly_p

    integer        , intent(in)    :: iYaml, iDbg
    type(SupElType), intent(inout) :: sups(:)
    real(dp)       , intent(in)    :: eigScale
    integer        , intent(out)   :: ierr

    !! Local variables
    integer :: i

    !! --- Logic section ---

    ierr = 0
    write(iYaml,600) eigScale
    do i = 1, size(sups)
       if (IsBeam(sups(i))) then

          !! Update co-rotated superelement coordinate system
          call UpdateSupElCorot (sups(i),0,ierr)
          if (ierr < 0) cycle

          !! Calculate current deformation vector in co-rotated system
          call BuildFinit (sups(i))

          !! Calculate the stiffness forces (stress resultants) in the beam
          sups(i)%FS = matmul(sups(i)%KmMat,sups(i)%finit)
          write(iYaml,601) sups(i)%id%baseId,sups(i)%FS(1:12)

          if (iDbg > 0) then
             write(iDbg,"('')")
             call writeObject (sups(i)%supTr,iDbg, &
                  'Updated coordinate system for '//getSupElId(sups(i)))
             call writeObject (sups(i)%finit,iDbg,'Nodal deformations')
             call writeObject (sups(i)%FS,iDbg,'Resulting stiffness forces')
          end if

       end if
    end do
    if (ierr < 0) call reportError (debugFileOnly_p,'writeYAMLforces')

600 format(/'  Eigenvector scale:',1PE12.5, &
         & /'  Stress resultants in system beams: # BeamID,', &
         & ' N1, Vy1, Vz1, Mx1, Myy1, Mzz1, N2, Vy2, Vz2, Mx2, Myy2, Mzz2', &
         & ' [N N N Nm Nm Nm N N N Nm Nm Nm]')
601 format('  - [',I6,1P,12(',',E12.5),']')

  end subroutine writeYAMLforces

end module ModesRoutinesModule
