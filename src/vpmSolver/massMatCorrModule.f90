!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file massMatCorrModule.f90
!>
!> @brief Subroutines for mass matrix correction due to high-speed rotations.

!!==============================================================================
!> @brief Module with subroutines for superelement mass matrix correction.
!>
!> @details This module contains some subroutines for calculating mass matrix
!> corrections for superelements what undergo high-speed rotations.

module MassMatrixCorrectionModule

  implicit none

  private :: mmcMassMatrixCorrection


contains

  !!============================================================================
  !> @brief Generates the 6 rigid body acceleration vectors for a superelement.
  !>
  !> @param[in] supel The superelement to generate the acceleration vectors for
  !> @param[in] rotCenter The rotation center of the superelement
  !> @param[out] rigVdd Rigid body acceleration vectors of the superelement
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen                                         @date Oct 1999

  subroutine mmcRigAccelVectors (supel,rotCenter,rigVdd)

    use SupElTypeModule, only : SupElType, dp

    type(SupElType), intent(in)  :: supel
    real(dp)       , intent(in)  :: rotCenter(3)
    real(dp)       , intent(out) :: rigVdd(:,:)

    !! Local variables
    integer  :: i, j, ip
    real(dp) :: xc, yc, zc

    !! --- Logic section ---

    rigVdd = 0.0_dp

    xc = rotCenter(1)
    yc = rotCenter(2)
    zc = rotCenter(3)

    do i = 1, supel%nExtNods
       ip = supel%triads(i)%firstDOF - 1

       !! Rigid translation in X-direction
       rigVdd(ip+1,1) = 1.0_dp

       !! Rigid translation in Y-direction
       rigVdd(ip+2,2) = 1.0_dp

       !! Rigid translation in Z-direction
       rigVdd(ip+3,3) = 1.0_dp

       !! Rigid rotation about an X-axis through the mass center
       rigVdd(ip+2,4) = -(supel%TrUndeformed(3,4,i) - zc)
       rigVdd(ip+3,4) =  (supel%TrUndeformed(2,4,i) - yc)

       !! Rigid rotation about an Y-axis through the mass center
       rigVdd(ip+1,5) =  (supel%TrUndeformed(3,4,i) - zc)
       rigVdd(ip+3,5) = -(supel%TrUndeformed(1,4,i) - xc)

       !! Rigid rotation about an Z-axis through the mass center
       rigVdd(ip+1,6) = -(supel%TrUndeformed(2,4,i) - yc)
       rigVdd(ip+2,6) =  (supel%TrUndeformed(1,4,i) - xc)

       if (supel%triads(i)%p%nDOFs >= 6) then
          do j = 4, 6
             rigVdd(ip+j,j) = 1.0_dp
          end do
       end if

    end do

  end subroutine mmcRigAccelVectors


  !!============================================================================
  !> @brief Calculates the rigid body mass properties for a superelement.
  !>
  !> @param supel The superelement to calculate mass properties for
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen                                      @date  1 Oct 1997
  !> @author Karl Erik Thoresen                                @date  4 Jun 1998
  !> @author Bjorn Haugen                                      @date 14 Oct 1999
  !> @author Knut Morten Okstad                                @date 21 Jan 2004

  subroutine mmcRigidMassProperties (supel,ierr)

    use KindModule        , only : dp, epsDiv0_p
    use SupElTypeModule   , only : SupElType
    use scratchArrayModule, only : getRealScratchMatrix
    use reportErrorModule , only : reportError, debugFileOnly_p

    type(SupElType), intent(inout) :: supel
    integer        , intent(out)   :: ierr

    !! Local variables
    real(dp)          :: RTMR(6,6)
    real(dp), pointer :: R(:,:)

    !! --- Logic section ---

    R => getRealScratchMatrix(supel%nTotDofs,6,ierr)
    if (ierr < 0) then
       call reportError (debugFileOnly_p,'mmcRigidMassProperties')
       return
    end if

    !! Establish rigid body acceleration vectors about origin
    supel%posMassCenter = 0.0_dp
    call mmcRigAccelVectors (supel,supel%posMassCenter,R)

    !! Compute RTMR - characteristic inertia matrix about origin
    RTMR = matmul(transpose(R),matmul(supel%Mmat,R))

    supel%mass = (RTMR(1,1) + RTMR(2,2) + RTMR(3,3))/3.0_dp ! Total mass

    !! Compute mass center from the inertia-matrix coupling terms
    !! between translation and rotation
    if (supel%mass > epsDiv0_p) then
       supel%posMassCenter(1) = RTMR(2,6)/supel%mass
       supel%posMassCenter(2) = RTMR(3,4)/supel%mass
       supel%posMassCenter(3) = RTMR(1,5)/supel%mass
    end if

    !! Establish rigid body acceleration vectors about mass center
    call mmcRigAccelVectors (supel,supel%posMassCenter,R)

    !! Compute characteristic inertia matrix about mass center
    supel%inertia = matmul(transpose(R),matmul(supel%Mmat,R))

  end subroutine mmcRigidMassProperties


  !!============================================================================
  !> @brief Initializes the matrix  (R^t*R)^-1*R^t  for mass matrix correction.
  !>
  !> @param supel The superelement to calculate the rtr_r matrix for
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Karl Erik Thoresen / Bjorn Haugen                 @date 1990s
  !> @author Knut Morten Okstad                                @date 21 Jan 2004

  subroutine mmcInitMassMatrixCorrection (supel,ierr)

    use SupElTypeModule   , only : SupElType, dp
    use IdTypeModule      , only : getId
    use DenseMatrixModule , only : solveAxB
    use scratchArrayModule, only : getRealScratchMatrix
    use reportErrorModule , only : AllocationError, reportError, debugFileOnly_p

    type(SupElType), intent(inout) :: supel
    integer        , intent(out)   :: ierr

    !! Local variables
    integer           :: ndofs
    real(dp)          :: RTR(6,6)
    real(dp), pointer :: R(:,:)

    !! --- Logic section ---

    ndofs = supel%nTotDOFs
    if (associated(supel%genDOFs)) then
       ndofs = ndofs - supel%genDOFs%nDOFs
    end if
    allocate(supel%rtr_rt(6,ndofs),STAT=ierr)
    if (ierr /= 0) then
       ierr = AllocationError('mmcInitMassMatrixCorrection')
       return
    end if

    R => getRealScratchMatrix(ndofs,6,ierr)
    if (ierr < 0) goto 900

    !! Establish rigid body acceleration vectors about mass center
    call mmcRigAccelVectors (supel,supel%posMassCenter,R)

    !! Multiply R^t*R
    supel%rtr_rt = transpose(R)
    RTR = matmul(supel%rtr_rt,R)

    !! Solve for the correction matrix = (R^t*R)^-1 * R^t
    call solveAxB (RTR,supel%rtr_rt,ierr)
    if (ierr == 0) return

900 call reportError (debugFileOnly_p,'mmcInitMassMatrixCorrection'// &
         &            ' for FE Part'//getId(supel%id))

  end subroutine mmcInitMassMatrixCorrection


  !!============================================================================
  !> @brief Calculates the error in the mass matrix.
  !>
  !> @param[in] supel The superelement to calculate the mass matrix error for
  !> @param[in] elVd Superelement velocity vector in local axis directions
  !> @param[out] momentError Torque error estimate for the superelement
  !> @param[out] forceError Force error estimate for the superelement
  !> @param[out] radiusVec Vector from axis of rotation to center of gravity
  !> @param[out] ErrNorm Error norms
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Karl Erik Thoresen / Bjorn Haugen                 @date 1990s
  !> @author Knut Morten Okstad                                @date 21 Jan 2004

  subroutine mmcMassMatrixCorrection (supel,elVd,momentError,forceError, &
       &                              radiusVec,ErrNorm,ierr)

    use SupElTypeModule   , only : SupElType, dp
    use manipMatrixModule , only : cross_product
    use scratchArrayModule, only : getRealScratchMatrix
    use rotationModule    , only : spin
    use reportErrorModule , only : reportError, debugFileOnly_p

    type(SupElType)  , intent(in)  :: supel
    real(dp)         , intent(in)  :: elVd(:)
    real(dp)         , intent(out) :: momentError(3), forceError(3)
    real(dp),optional, intent(out) :: radiusVec(3), ErrNorm(2)
    integer ,optional, intent(out) :: ierr

    !! Local variables
    integer             :: i, j, ndofs, nnods
    real(dp), parameter :: eps_p = tiny(1.0_dp)*1000.0_dp
    real(dp), pointer   :: rigVec(:,:), rotAcc(:)
    real(dp)            :: CR(3),spinO(3,3),spinO2(3,3),a(6),RTMRot(6), &
         &                 Energy,dE,n1(3),n2(3),n3(3),Uct,omega,omegaVec(3), &
         &                 Fcorrect(3),Tcorrect(3),Lc(3),len,radius

    !! --- Logic section ---

    momentError = 0.0_dp
    forceError  = 0.0_dp
    if (present(radiusVec)) radiusVec = 0.0_dp
    if (present(ErrNorm))   ErrNorm   = 0.0_dp
    if (present(ierr))      ierr      = 0

    !! Get local scratch array
    nnods = supel%nExtNods
    ndofs = supel%nTotDOFs
    if (associated(supel%genDOFs)) then
       ndofs = ndofs - supel%genDOFs%nDOFs
    end if
    rigVec => getRealScratchMatrix(ndofs,7,i)
    if (i < 0) then
       if (present(ierr)) ierr = i
       call reportError (debugFileOnly_p,'mmcMassMatrixCorrection')
       return
    end if
    rotAcc => rigVec(:,7)

    !! Find the rigid body velocity vector for the element (store in a)
    a = matmul(supel%rtr_rt,elVd(1:ndofs))

    omegaVec = a(4:6)                            ! Vector along rotation axis
    omega = sqrt(dot_product(omegaVec,omegaVec)) ! Angular velocity
    if (omega < eps_p) return                    ! No corrections necessary

    n1 = omegaVec/omega                          ! Normalize the vector
    n3 = cross_product(n1,a)
    len = sqrt(dot_product(n3,n3))
    if (len > eps_p) n3 = n3/len     ! Scale to unity

    n2     = cross_product(n3,n1)    ! Already unity
    Uct    = dot_product(a(1:3),n2)  ! Velocity normal to rotation axis
    radius = Uct/omega               ! Distance from rotation axis to CG
    if (present(radiusVec)) radiusVec = -radius*n3 ! Vector from rot. axis to CG

    CR = supel%posMassCenter + radius*n3 ! Point on the axis of rotation

    call spin (omegaVec,spinO)     ! Compute spin of the angular velocity vector
    spinO2 = matmul(spinO,spinO)

    !! Rigid body acceleration vectors about mass-center
    call mmcRigAccelVectors (supel,CR,rigVec(:,1:6))

    !! Establish element acceleration vector for angular velocity omega about CR
    rotAcc = 0.0_dp
    do i = 1, nnods
       j = supel%triads(i)%firstDOF
       !! Nodal acceleration vector for angular velocity omega about CR
       rotAcc(j:j+2) = matmul(spinO2,supel%TrUndeformed(:,4,i)-CR)
    end do

    !! Compute RTMRot = R'*M*a
    RTMRot = matmul(matmul(supel%Mmat(1:ndofs,1:ndofs),rotAcc),rigVec(:,1:6))

    !! RTMRot now contains the force vector at CR given by the rotation a
    !! Calculate the centrifugal moments around the center of rotation

    !! Compute angular momentum about the mass center
    Lc = matmul(supel%inertia(4:6,4:6),omegaVec)

    !! Compute correct torque as rate of change of angular momentum
    Tcorrect = matmul(spinO,Lc)

    !! Compute correct reaction force as : f = m*a = m*( w x ( w x r ) )
    Fcorrect = matmul(spinO2,-n3)*radius*supel%mass

    !! Torqe error estimate as finite element torqe minus rigid body torque
    momentError = RTMRot(4:6) - Tcorrect

    !! Correction in resulting force as finite element resultant force
    !! minus correct force from rigid body dynamics
    forceError = RTMRot(1:3) - Fcorrect

    if (.not.present(ErrNorm)) return

    !! Calculate an error norm
    Energy = (dot_product(Lc,omegaVec) + &
         &    dot_product(a(1:3),a(1:3))*supel%mass)*0.5_dp
    dE     =  dot_product(MomentError,omegaVec)

    if (abs(Energy) > eps_p) ErrNorm(1) = dE/Energy
    ErrNorm(2) = sqrt(dot_product(MomentError,MomentError))

  end subroutine mmcMassMatrixCorrection


  !!============================================================================
  !> @brief Checks the error in the mass matrix and print warning if increasing.
  !>
  !> @param supel The superelement to check mass matrix error for
  !> @param[in] uld Superelement velocity vector in local axis directions
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Karl Erik Thoresen / Bjorn Haugen                 @date 1990s
  !> @author Knut Morten Okstad                                @date 21 Jan 2004

  subroutine mmcMassMatrixWarning (supel,uld,ierr)

    use SupElTypeModule       , only : SupElType, dp
    use IdTypeModule          , only : getId
    use fileUtilitiesModule   , only : getDBGfile
    use dbgUnitsModule        , only : dbgRot
    use reportErrorModule     , only : reportError, debugFileOnly_p, warning_p
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_intValue

    type(SupElType), intent(inout) :: supel
    real(dp)       , intent(in)    :: uld(:)
    integer        , intent(out)   :: ierr

    !! Local variables
    real(dp), parameter :: minRotError_p = 0.1_dp
    real(dp)            :: MomentError(3), forceError(3), errNorm(2)

    !! --- Logic section ---

    if (dbgRot == 0 .and. ffa_cmdlinearg_intValue('debug') > 0) then
       !! Open the centripetal error file
       dbgRot = getDBGfile(4,'centr_err.lst')
       write(dbgRot,100)
100    format(/'  Before/without centripetal moment correction '&
            & /'  -------------------------------------------- '&
            &//'  Part   dE/Ekin        Tx error       Ty error       Tz error')
    end if

    call mmcMassMatrixCorrection (supel, uld, momentError, forceError, &
         &                        ErrNorm=errNorm, ierr=ierr)
    if (ierr < 0) then
       call reportError (debugFileOnly_p,'mmcMassMatrixWarning')
       return
    end if

    if (dbgRot > 0) then
       write(dbgRot,'(i3,4e15.3)') supel%id%userId,errNorm(1),MomentError
    end if

    if (supel%massCorrFlag > 0 .or. supel%massCorrFlag <= -5) return

    !! Print warning if moment error is increasing in 5 subsequent time steps

    if (abs(errNorm(1)) > max(supel%RotError,minRotError_p)) then

       !! This error is increasing
       supel%RotError = abs(errNorm(1))
       supel%massCorrFlag = supel%massCorrFlag - 1

       if (supel%massCorrFlag == -5) then
          !! The error has been increasing during the 5 last time steps
          call reportError (warning_p, &
               'Indication of a poor centripetal force on the rotating Part'// &
               getId(supel%id), 'Try distributing the triads better '// &
               'or use centripetal moment correction.')
       end if

    else
       !! There was not detected any error increase in this time step
       supel%RotError = 0.0_dp
       supel%massCorrFlag = 0
    end if

  end subroutine mmcMassMatrixWarning


  !!============================================================================
  !> @brief
  !> Returns a correction force vector that accounts for the mass matrix error.
  !>
  !> @param[in] sup The superelement to calculate for correction vector for
  !>
  !> @callgraph @callergraph
  !>
  !> @author Karl Erik Thoresen / Bjorn Haugen                 @date 1990s
  !> @author Knut Morten Okstad                                @date 21 Jan 2004

  function mmcGetMassTorqueCorrection (sup) result(Qk)

    use SupElTypeModule  , only : SupElType, dp
    use manipMatrixModule, only : cross_product

    type(SupElType), intent(in) :: sup
    real(dp)                    :: Qk(sup%nTotDofs)

    !! Local variabels
    integer  :: i, j, n, nnods, nrnods
    real(dp) :: momentError(3), forceError(3), radiusVec(3)
    real(dp) :: fVec(6), fMom(3), rVec(3)

    !! --- Logic section ---

    Qk = 0.0_dp
    call mmcMassMatrixCorrection (sup,sup%uld,momentError,forceError,radiusVec)

    !! Count number of nodes with rotational dofs ( 6 dofs )
    !! Find the average nodal position
    nnods  = sup%nExtNods
    nrnods = 0
    rVec = 0.0_dp
    do i = 1, nnods
       if (sup%triads(i)%p%nDOFs >= 6) nrnods = nrnods + 1
       rVec = rVec + sup%TrUndeformed(:,4,i) - sup%posMassCenter
    end do
    rVec = radiusVec + rVec/dble(nnods)
    fMom = cross_product(rVec,forceError)

    !! Correct for moment error by subtracting torque from all nodes with
    !! rotational dofs
    fVec(1:3) =  forceError/dble(nnods)
    fVec(4:6) = (momentError-fMom)/dble(nrnods)
    do i = 1, nnods
       j = sup%triads(i)%firstDOF
       n = sup%triads(i)%p%nDOFs
       Qk(j:j+n-1) = fVec(1:n)
    end do

  end function mmcGetMassTorqueCorrection

end module MassMatrixCorrectionModule
