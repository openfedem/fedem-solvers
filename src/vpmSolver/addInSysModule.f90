!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file addInSysModule.f90
!>
!> @brief Subroutines for system-level assembly.

!!==============================================================================
!> @brief Module with subroutines for system-level assembly.
!>
!> @details This module contains subroutines for assembling the system
!> left-hand-side matrices and associated right-hand-side force vectors
!> of the linearized equation system.

module AddInSysModule

  implicit none

contains

  !!============================================================================
  !> @brief Builds the system Newton matrix.
  !>
  !> @param Nmat The system Newton matrix
  !> @param[in] scaleM Scaling factor for the system mass matrix
  !> @param[in] scaleC Scaling factor for the system damping matrix
  !> @param[in] scaleK Scaling factor for the system stiffness matrix
  !> @param mech Mechanism components of the model
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[in] iter Iteration counter
  !> @param[in] stressStiffUpdateSkip Number of initial iterations without
  !> stress stiffening in each time step
  !> @param[in] globalStressStiffIsOn Global stress stiffening flag
  !> @param[in] alpha2 Stiffness-proportional damping factor for axial springs
  !> @param[out] ierr Error flag
  !> @param Rhs Right-hand-side force vector
  !>
  !> @details The Newton matrix is a linear combination of the mass-, damping-
  !> and stiffness matrices, Nmat = scaleM*M + scaleC*C + scaleK*K.
  !> If @a Rhs is provided, contributions due to prescribed motions,
  !> if any, will be added.
  !>
  !> @author Karl Erik Thoresen                                   @date Feb 1999
  !> @author Bjorn Haugen                                         @date Des 2000
  !> @author Knut Morten Okstad                                   @date Jul 2002

  subroutine BuildNewtonMat (Nmat, scaleM, scaleC, scaleK, mech, sam, iter, &
       &                     stressStiffUpdateSkip, globalStressStiffIsOn,  &
       &                     alpha2, ierr, Rhs)

    use SamModule                   , only : SamType, dp
    use SysMatrixTypeModule         , only : SysMatrixType
    use MechanismTypeModule         , only : MechanismType
    use SupElRoutinesModule         , only : BuildSupNewtonMat, addInSupMat
    use UserdefElRoutinesModule     , only : BuildUDENewtonMat, addInUDEMat
    use BushingElementRoutinesModule, only : addInBushingElementMat
    use ContactElementRoutinesModule, only : addInContactElementStiffMat
    use ContactElementRoutinesModule, only : addInContactElementDamperMat
    use SpringRoutinesModule        , only : addInSpringStiffMat
    use DamperRoutinesModule        , only : addInDamperMat
    use MassRoutinesModule          , only : addInMassMat
    use TireRoutinesModule          , only : addInTireMassMat
    use TireRoutinesModule          , only : addInTireStiffMat
    use TireRoutinesModule          , only : addInTireDamperMat
    use IdTypeModule                , only : getId
    use profilerModule              , only : startTimer, stopTimer, asm_p
    use reportErrorModule           , only : warning_p, debugFileOnly_p
    use reportErrorModule           , only : reportError
#ifdef FT_DEBUG
    use dbgUnitsModule        , only : dbgSolve
    use SysMatrixTypeModule   , only : writeObject
    use manipMatrixModule     , only : writeObject
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getint
#endif

    type(SysMatrixType), intent(inout) :: Nmat
    real(dp)           , intent(in)    :: scaleM, scaleC, scaleK
    type(MechanismType), intent(inout) :: mech
    type(SamType)      , intent(in)    :: sam
    integer            , intent(in)    :: iter, stressStiffUpdateSkip
    logical            , intent(in)    :: globalStressStiffIsOn
    real(dp)           , intent(in)    :: alpha2
    integer            , intent(out)   :: ierr
    real(dp), optional , intent(inout) :: Rhs(:)

    !! Local variables
    integer             :: i
    real(dp)            :: scaleD
    logical             :: reComputeSupMat, newTanStiff
    logical , parameter :: lDynamics = .true.
    real(dp), parameter :: eps_p = 1.0e-16_dp
    real(dp), save      :: oldScaleM = 0.0_dp
    real(dp), save      :: oldScaleK = 0.0_dp
    real(dp), save      :: oldScaleC = 0.0_dp
    logical , save      :: firstMsg = .true.

#ifdef FT_DEBUG
    integer :: iprint
    if (dbgSolve > 0) then
       call ffa_cmdlinearg_getint ('debug',iprint)
       write(dbgSolve,"(/'   --- in BuildNewtonMat, iter =',I3)") iter
       write(dbgSolve,"('       scaleM,scaleC,scaleK :',1p3e12.4)") &
            &                   scaleM,scaleC,scaleK
    else
       iprint = 0
    end if
    if (iprint == -99 .and. present(Rhs)) then
       write(dbgSolve,*)
       call writeObject (Rhs,dbgSolve,'Incoming system RHS vector')
       write(dbgSolve,*)
    end if
#endif

    !! --- Logic section ---

    call startTimer (asm_p)

    ierr = 0
    Nmat%value = 0.0_dp

    if (Nmat%lFirst) then
       reComputeSupMat = .true.
       newTanStiff = .true.
    else
       !! Check if factors have changed, recompute Newton matrices if so
       reComputeSupMat = abs(scaleM-oldScaleM) + &
            &            abs(scaleK-oldScaleK) + &
            &            abs(scaleC-oldScaleC) > 1.0e-16_dp
       newTanStiff = iter == 0 .or. iter > stressStiffUpdateSkip
    end if

    if (reComputeSupMat) then
       oldScaleM = scaleM
       oldScaleK = scaleK
       oldScaleC = scaleC
    end if

    do i = 1, size(mech%sups)

       if ( reComputeSupMat .or. mech%sups(i)%stressStiffFlag(1) == 1 .or. &
            abs(mech%sups(i)%massScl(1)-mech%sups(i)%massScl(2)) > eps_p .or. &
            abs(mech%sups(i)%stifScl(1)-mech%sups(i)%stifScl(2)) > eps_p ) then

          call BuildSupNewtonMat (newTanStiff, scaleM, scaleC, scaleK, &
               &                  mech%sups(i), ierr)
          if (ierr < 0) goto 900

       end if

       call addInSupMat (mech%sups(i)%Nmat, Nmat, mech%sups(i), sam, ierr, Rhs)
       if (ierr < 0) goto 900

    end do

#ifdef FT_DEBUG
    if (iprint == -99 .and. size(mech%sups) > 0) then
       write(dbgSolve,100) 'superelement'
       call writeObject (Nmat,sam%mpar,dbgSolve)
       if (present(Rhs)) call writeObject (Rhs,dbgSolve,'System RHS vector')
    end if
100 format(/'System Newton matrix after ',A,' assembly')
#endif

    newTanStiff = globalStressStiffIsOn .and. newTanStiff

    do i = 1, size(mech%elms)

       call BuildUDENewtonMat (scaleM, scaleC, scaleK, mech%elms(i), ierr)
       if (ierr < 0) goto 900

       call addInUDEMat (mech%elms(i)%Nmat, Nmat, mech%elms(i), sam, ierr, Rhs)
       if (ierr < 0) goto 900

    end do

#ifdef FT_DEBUG
    if (iprint == -99 .and. size(mech%elms) > 0) then
       write(dbgSolve,100) 'userdefined element'
       call writeObject (Nmat,sam%mpar,dbgSolve)
       if (present(Rhs)) call writeObject (Rhs,dbgSolve,'System RHS vector')
    end if
#endif

    do i = 1, size(mech%bElems)
       call addInBushingElementMat (scaleK, scaleC, Nmat, &
            &                       mech%bElems(i), sam, ierr, Rhs)
       if (ierr < 0) goto 900
    end do

#ifdef FT_DEBUG
    if (iprint == -99 .and. size(mech%bElems) > 0) then
       write(dbgSolve,100) 'bushing element'
       call writeObject (Nmat,sam%mpar,dbgSolve)
       if (present(Rhs)) call writeObject (Rhs,dbgSolve,'System RHS vector')
    end if
#endif

    do i = 1, size(mech%cElems)
       call addInContactElementStiffMat (newTanStiff, scaleK, Nmat, &
            &                            mech%cElems(i), sam, ierr, Rhs)
       if (ierr < 0) goto 900
       call addInContactElementDamperMat (newTanStiff, scaleC, scaleK, Nmat, &
            &                             mech%cElems(i), sam, ierr, Rhs)
       if (ierr < 0) goto 900
    end do

#ifdef FT_DEBUG
    if (iprint == -99 .and. size(mech%cElems) > 0) then
       write(dbgSolve,100) 'contact element'
       call writeObject (Nmat,sam%mpar,dbgSolve)
       if (present(Rhs)) call writeObject (Rhs,dbgSolve,'System RHS vector')
    end if
#endif

    do i = 1, size(mech%axialSprings)
       scaleD = scaleC*(mech%axialSprings(i)%alpha2 + alpha2)
       if (newTanStiff .and. abs(scaleD) > 1.0e-16_dp) then
          call addInSpringStiffMat (.true., scaleK, Nmat, &
               &                    mech%axialSprings(i), sam, ierr, Rhs)
          if (ierr < 0) goto 900
          if (.not.associated(mech%axialSprings(i)%spr(1)%p%length0Engine)) then
             call addInSpringStiffMat (.false., scaleD, Nmat, &
                  &                    mech%axialSprings(i), sam, ierr, Rhs)
             if (ierr < 0) goto 900
          else if (firstMsg) then
             call reportError (warning_p,'Stiffness-proportional damping '//&
                  'for Axial spring'//trim(getId(mech%axialSprings(i)%id)), &
                  'is deactivated due to stress-free length change.', &
                  'Use Axial an damper instead.')
             firstMsg = .false.
          end if
       else
          call addInSpringStiffMat (newTanStiff, scaleK+scaleD, Nmat, &
               &                    mech%axialSprings(i), sam, ierr, Rhs)
          if (ierr < 0) goto 900
       end if
    end do
    do i = 1, size(mech%joints)
       if (associated(mech%joints(i)%springEl)) then
          call addInSpringStiffMat (newTanStiff, scaleK, Nmat, &
               &                    mech%joints(i)%springEl, sam, ierr, Rhs)
          if (ierr < 0) goto 900
       end if
       if (associated(mech%joints(i)%sprFric)) then
          call addInSpringStiffMat (newTanStiff, scaleK, Nmat, &
               &                    mech%joints(i)%sprFric, sam, ierr, Rhs)
          if (ierr < 0) goto 900
       end if
    end do

#ifdef FT_DEBUG
    if (iprint == -99 .and. size(mech%axialSprings) > 0) then
       write(dbgSolve,100) 'spring'
       call writeObject (Nmat,sam%mpar,dbgSolve)
       if (present(Rhs)) call writeObject (Rhs,dbgSolve,'System RHS vector')
    end if
#endif

    do i = 1, size(mech%dampers)
       call addInDamperMat (newTanStiff, scaleC, scaleK, Nmat, &
            &               mech%dampers(i), sam, ierr, Rhs)
       if (ierr < 0) goto 900
    end do

#ifdef FT_DEBUG
    if (iprint == -99 .and. size(mech%dampers) > 0) then
       write(dbgSolve,100) 'damper'
       call writeObject (Nmat,sam%mpar,dbgSolve)
       if (present(Rhs)) call writeObject (Rhs,dbgSolve,'System RHS vector')
    end if
#endif

    do i = 1, size(mech%masses)
       call addInMassMat (scaleM, Nmat, mech%masses(i), sam, ierr, Rhs)
       if (ierr < 0) goto 900
    end do

#ifdef FT_DEBUG
    if (iprint == -99 .and. size(mech%masses) > 0) then
       write(dbgSolve,100) 'additional mass'
       call writeObject (Nmat,sam%mpar,dbgSolve)
       if (present(Rhs)) call writeObject (Rhs,dbgSolve,'System RHS vector')
    end if
#endif

    do i = 1, size(mech%tires)
       call addInTireStiffMat (lDynamics, scaleK, Nmat, &
            &                  mech%tires(i), sam, ierr, Rhs)
       if (ierr < 0) goto 900
       call addInTireDamperMat (scaleC, Nmat, mech%tires(i), sam, ierr, Rhs)
       if (ierr < 0) goto 900
       call addInTireMassMat (scaleM, Nmat, mech%tires(i), sam, ierr, Rhs)
       if (ierr < 0) goto 900
    end do

#ifdef FT_DEBUG
    if (iprint == -99 .and. size(mech%tires) > 0) then
       write(dbgSolve,100) 'tire'
       call writeObject (Nmat,sam%mpar,dbgSolve)
       if (present(Rhs)) call writeObject (Rhs,dbgSolve,'System RHS vector')
    end if
#endif

    Nmat%lFirst = .false.

900 call stopTimer (asm_p)
    if (ierr < 0) call reportError (debugFileOnly_p,'BuildNewtonMat')

  end subroutine BuildNewtonMat


  !!============================================================================
  !> @brief Builds the system mass matrix.
  !>
  !> @param Mmat The system mass matrix
  !> @param[in] mech Mechanism components of the model
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[out] ierr Error flag
  !>
  !> @author Karl Erik Thoresen                                   @date Aug 1999
  !> @author Knut Morten Okstad                                   @date Jul 2002

  subroutine BuildMassMat (Mmat, mech, sam, ierr)

    use SamModule              , only : SamType, dp
    use SysMatrixTypeModule    , only : SysMatrixType
    use MechanismTypeModule    , only : MechanismType
    use SupElRoutinesModule    , only : addInSupMat
    use UserdefElRoutinesModule, only : addInUDEMat
    use MassRoutinesModule     , only : addInMassMat
    use TireRoutinesModule     , only : addInTireMassMat
    use profilerModule         , only : startTimer, stopTimer, asm_p
    use reportErrorModule      , only : reportError, debugFileOnly_p

    type(SysMatrixType), intent(inout) :: Mmat
    type(MechanismType), intent(in)    :: mech
    type(SamType)      , intent(in)    :: sam
    integer            , intent(out)   :: ierr

    !! Local variables
    integer             :: i
    real(dp), parameter :: scaleI = 1.0_dp

    !! --- Logic section ---

    call startTimer (asm_p)

    ierr = 0
    Mmat%value = 0.0_dp

    do i = 1, size(mech%sups)
       if (associated(mech%sups(i)%Mmat)) then
          call addInSupMat (mech%sups(i)%Mmat, Mmat, mech%sups(i), sam, ierr, &
               &            scale=mech%sups(i)%massScl(1))
          if (ierr < 0) goto 900
       end if
       if (associated(mech%sups(i)%Mamat)) then
          call addInSupMat (mech%sups(i)%Mamat, Mmat, mech%sups(i), sam, ierr)
          if (ierr < 0) goto 900
       end if
    end do

    do i = 1, size(mech%elms)
       call addInUDEMat (mech%elms(i)%Mmat, Mmat, mech%elms(i), sam, ierr)
       if (ierr < 0) goto 900
    end do

    do i = 1, size(mech%masses)
       call addInMassMat (scaleI, Mmat, mech%masses(i), sam, ierr)
       if (ierr < 0) goto 900
    end do

    do i = 1, size(mech%tires)
       call addInTireMassMat (scaleI, Mmat, mech%tires(i), sam, ierr)
       if (ierr < 0) goto 900
    end do

900 call stopTimer (asm_p)
    if (ierr < 0) call reportError (debugFileOnly_p,'BuildMassMat')

  end subroutine BuildMassMat


  !!============================================================================
  !> @brief Builds the system damping matrix.
  !>
  !> @param Cmat The system damping matrix
  !> @param[in] mech Mechanism components of the model
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[out] ierr Error flag
  !>
  !> @author Knut Morten Okstad                                   @date Jul 2002

  subroutine BuildDamperMat (Cmat, mech, sam, ierr)

    use SamModule                   , only : SamType, dp
    use SysMatrixTypeModule         , only : SysMatrixType
    use MechanismTypeModule         , only : MechanismType
    use SupElRoutinesModule         , only : addInSupMat
    use UserdefElRoutinesModule     , only : addInUDEMat
    use BushingElementRoutinesModule, only : addInBushingElementMat
    use ContactElementRoutinesModule, only : addInContactElementDamperMat
    use DamperRoutinesModule        , only : addInDamperMat
    use TireRoutinesModule          , only : addInTireDamperMat
    use profilerModule              , only : startTimer, stopTimer, asm_p
    use reportErrorModule           , only : reportError, debugFileOnly_p

    type(SysMatrixType), intent(inout) :: Cmat
    type(MechanismType), intent(in)    :: mech
    type(SamType)      , intent(in)    :: sam
    integer            , intent(out)   :: ierr

    !! Local variables
    integer             :: i
    real(dp), parameter :: scaleK = 0.0_dp, scaleC = 1.0_dp

    !! --- Logic section ---

    call startTimer (asm_p)

    ierr = 0
    Cmat%value = 0.0_dp

    do i = 1, size(mech%sups)
       if (associated(mech%sups(i)%Cmat)) then
          !! We here assume that only stiffness-proportional damping is used,
          !! unless the stiffness-scaling factor is equal to one
          call addInSupMat (mech%sups(i)%Cmat, Cmat, mech%sups(i), sam, ierr, &
               &            scale=mech%sups(i)%stifScl(1))
          if (ierr < 0) goto 900
       end if
       if (associated(mech%sups(i)%Cdmat)) then
          call addInSupMat (mech%sups(i)%Cdmat, Cmat, mech%sups(i), sam, ierr)
          if (ierr < 0) goto 900
       end if
    end do

    do i = 1, size(mech%elms)
       if (associated(mech%elms(i)%Cmat)) then
          call addInUDEMat (mech%elms(i)%Cmat, Cmat, mech%elms(i), sam, ierr)
          if (ierr < 0) goto 900
       end if
    end do

    do i = 1, size(mech%bElems)
       call addInBushingElementMat (scaleK, scaleC, Cmat, &
            &                       mech%bElems(i), sam, ierr)
       if (ierr < 0) goto 900
    end do

    do i = 1, size(mech%cElems)
       call addInContactElementDamperMat (.false., scaleC, scaleK, Cmat, &
            &                             mech%cElems(i),sam,ierr)
       if (ierr < 0) goto 900
    end do

    do i = 1, size(mech%dampers)
       call addInDamperMat (.false., scaleC, scaleK, Cmat, &
            &               mech%dampers(i), sam, ierr)
       if (ierr < 0) goto 900
    end do

    do i = 1, size(mech%tires)
       call addInTireDamperMat (scaleC, Cmat, mech%tires(i), sam, ierr)
       if (ierr < 0) goto 900
    end do

900 call stopTimer (asm_p)
    if (ierr < 0) call reportError (debugFileOnly_p,'BuildDamperMat')

  end subroutine BuildDamperMat


  !!============================================================================
  !> @brief Builds the system stiffness matrix.
  !>
  !> @param Kmat The system stiffness matrix
  !> @param mech Mechanism components of the model
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[in] iter Iteration counter (&lt; 0 for eigenvalue analysis)
  !> @param[in] globalStressStiffIsOn Global stress stiffening flag
  !> @param[in] stressStiffUpdateSkip Number of initial iterations without
  !> stress stiffening in each time step
  !> @param[out] ierr Error flag
  !> @param Rhs Right-hand-side force vector
  !>
  !> @details If @a Rhs is provided, contributions due to prescribed motions,
  !> if any, will be added.
  !>
  !> @author Karl Erik Thoresen                                   @date Aug 1999
  !> @author Knut Morten Okstad                                   @date Jul 2002

  subroutine BuildStiffMat (Kmat, mech, sam, iter, globalStressStiffIsOn, &
       &                    stressStiffUpdateSkip, ierr, Rhs)

    use SamModule                   , only : SamType, dp
    use SysMatrixTypeModule         , only : SysMatrixType
    use MechanismTypeModule         , only : MechanismType
    use SupElRoutinesModule         , only : CompTanStiff, addInSupMat
    use UserdefElRoutinesModule     , only : addInUDEMat
    use BushingElementRoutinesModule, only : addInBushingElementMat
    use ContactElementRoutinesModule, only : addInContactElementStiffMat
    use SpringRoutinesModule        , only : addInSpringStiffMat
    use TireRoutinesModule          , only : addInTireStiffMat
    use profilerModule              , only : startTimer, stopTimer, asm_p
    use reportErrorModule           , only : reportError, debugFileOnly_p

    type(SysMatrixType), intent(inout) :: Kmat
    type(MechanismType), intent(inout) :: mech
    type(SamType)      , intent(in)    :: sam
    integer            , intent(in)    :: iter
    logical            , intent(in)    :: globalStressStiffIsOn
    integer            , intent(in)    :: stressStiffUpdateSkip
    integer            , intent(out)   :: ierr
    real(dp), optional , intent(inout) :: Rhs(:)

    !! Local variables
    integer             :: i
    logical             :: newTanStiff, lDynamics
    real(dp), parameter :: scaleK = 1.0_dp, scaleC = 0.0_dp

    !! --- Logic section ---

    call startTimer (asm_p)

    ierr = 0
    Kmat%value = 0.0_dp

    newTanStiff = iter <= 0 .or. iter > stressStiffUpdateSkip

    do i = 1, size(mech%sups)
       if ( (iter >= 0 .and. mech%sups(i)%stressStiffFlag(2) == 1) .or. &
            (iter <  0 .and. mech%sups(i)%stressStiffFlag(3) == 1) ) then

          if (newTanStiff) then
             call CompTanStiff (mech%sups(i), ierr)
             if (ierr < 0) goto 900
          end if

          call addInSupMat (mech%sups(i)%KtMat, Kmat, &
               &            mech%sups(i), sam, ierr, Rhs, &
               &            mech%sups(i)%stifScl(1))
       else
          call addInSupMat (mech%sups(i)%KmMat, Kmat, &
               &            mech%sups(i), sam, ierr, Rhs, &
               &            mech%sups(i)%stifScl(1))
       end if
       if (ierr < 0) goto 900
    end do

    newTanStiff = newTanStiff .and. globalStressStiffIsOn

    do i = 1, size(mech%elms)
       call addInUDEMat (mech%elms(i)%Kmat, Kmat, &
            &            mech%elms(i), sam, ierr, Rhs)
       if (ierr < 0) goto 900
    end do

    do i = 1, size(mech%bElems)
       call addInBushingElementMat (scaleK, scaleC, Kmat, &
            &                       mech%bElems(i), sam, ierr, Rhs)
       if (ierr < 0) goto 900
    end do

    do i = 1, size(mech%cElems)
       call addInContactElementStiffMat (newTanStiff, scaleK, Kmat, &
            &                            mech%cElems(i), sam, ierr, Rhs)
       if (ierr < 0) goto 900
    end do

    do i = 1, size(mech%axialSprings)
       call addInSpringStiffMat (newTanStiff, scaleK, Kmat, &
            &                    mech%axialSprings(i), sam, ierr, Rhs)
       if (ierr < 0) goto 900
    end do
    do i = 1, size(mech%joints)
       if (associated(mech%joints(i)%springEl)) then
          call addInSpringStiffMat (newTanStiff, scaleK, Kmat, &
               &                    mech%joints(i)%springEl, sam, ierr, Rhs)
          if (ierr < 0) goto 900
       end if
       if (associated(mech%joints(i)%sprFric)) then
          call addInSpringStiffMat (newTanStiff, scaleK, Kmat, &
               &                    mech%joints(i)%sprFric, sam, ierr, Rhs)
          if (ierr < 0) goto 900
       end if
    end do

    lDynamics = iter < 0 ! Use same stiffness in eigen-analysis as in dynamics
    do i = 1, size(mech%tires)
       call addInTireStiffMat (lDynamics, scaleK, Kmat, &
            &                  mech%tires(i), sam, ierr, Rhs)
       if (ierr < 0) goto 900
    end do

900 call stopTimer (asm_p)
    if (ierr < 0) call reportError (debugFileOnly_p,'BuildStiffMat')

  end subroutine BuildStiffMat


  !!============================================================================
  !> @brief Builds the system force vectors.
  !>
  !> @param[out] FSk Internal stiffness forces
  !> @param[out] FDk Internal damping forces
  !> @param[out] FIk Internal inertia forces
  !> @param[out] Qk External forces
  !> @param[out] RFk Reaction forces
  !> @param[in] sam Data for managing system matrix assembly
  !> @param mech Mechanism components of the model
  !> @param[in] istep Time increment counter
  !> @param[in] alpha2 Stiffness-proportional damping factor for axial springs
  !> @param[out] ierr Error flag
  !>
  !> @details this subroutine puts forces from all mechanism components
  !> into the associated system force vectors.
  !>
  !> @author Ole Ivar Sivertsen                                @date 19 Jan 1989
  !> @author Karl Erik Thoresen                                @date 1999
  !> @author Bjorn Haugen                                      @date 2001
  !> @author Knut Morten Okstad                                @date Jun 2002

  subroutine GetForceVectors (FSk, FDk, FIk, Qk, RFk, sam, mech, &
       &                      istep, alpha2, ierr)

    use KindModule                  , only : dp, i8, maxInt_p
    use SamModule                   , only : SamType
    use MechanismTypeModule         , only : MechanismType
    use SupElRoutinesModule         , only : addInSupForces
    use UserdefElRoutinesModule     , only : addInUDEForces
    use ForceRoutinesModule         , only : addInExternalForces
    use SpringRoutinesModule        , only : addInSpringForces
    use DamperRoutinesModule        , only : addInDamperForces
    use BushingElementRoutinesModule, only : addInBushingElementForces
    use ContactElementRoutinesModule, only : addInContactElementForces
    use FrictionRoutinesModule      , only : addInFrictionForces
    use FrictionRoutinesModule      , only : findForceInFriction
    use MassRoutinesModule          , only : addInInertiaForces
    use MassRoutinesModule          , only : addInGravitationForces
    use TireRoutinesModule          , only : addInTireForces
    use TireRoutinesModule          , only : addInTireInertiaForces
    use WindTurbineRoutinesModule   , only : addInAeroForces
    use FNVwaveForceModule          , only : addInFNVwaveForces
    use profilerModule              , only : startTimer, stopTimer, asm_p
    use reportErrorModule           , only : reportError, debugFileOnly_p
#ifdef FT_DEBUG
    use dbgUnitsModule        , only : dbgSolve
    use manipMatrixModule     , only : WriteObject
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getint
#endif

    real(dp)           , intent(out)   :: FSk(:), FDk(:), FIk(:), Qk(:), RFk(:)
    type(SamType)      , intent(in)    :: sam
    type(MechanismType), intent(inout) :: mech
    integer(i8)        , intent(in)    :: istep
    real(dp)           , intent(in)    :: alpha2
    integer            , intent(out)   :: ierr

    !! Local variables
    integer :: i
    logical :: haveSprings

#ifdef FT_DEBUG
    integer :: iprint
    if (dbgSolve > 0) then
       call ffa_cmdlinearg_getint ('debug',iprint)
       write(dbgSolve,"(/'   --- in GetForceVectors')")
    else
       iprint = 0
    end if
#endif

    !! --- Logic section ---

    call startTimer (asm_p)

    FSk = 0.0_dp
    FDk = 0.0_dp
    FIk = 0.0_dp
    Qk  = 0.0_dp
    RFk = 0.0_dp

    call addInSupForces (sam, mech%sups, FSk, FDk, FIk, Qk, RFk, ierr)
#ifdef FT_DEBUG
    if (iprint == -99 .and. size(mech%sups) > 0) then
       call writeObject (FSk,dbgSolve,'FS after superelement assembly')
       call writeObject (FDk,dbgSolve,'FD after superelement assembly')
       call writeObject (FIk,dbgSolve,'FI after superelement assembly')
       call writeObject (Qk ,dbgSolve,'Q after superelement assembly')
       call writeObject (RFk,dbgSolve,'RF after superelement assembly')
    end if
#endif

    call addInUDEForces (sam, mech%elms, FSk, FDk, FIk, Qk, RFk, ierr)
#ifdef FT_DEBUG
    if (iprint == -99 .and. size(mech%elms) > 0) then
       call writeObject (FSk,dbgSolve,'FS after userdefined element assembly')
       call writeObject (FDk,dbgSolve,'FD after userdefined element assembly')
       call writeObject (FIk,dbgSolve,'FI after userdefined element assembly')
       call writeObject (Qk ,dbgSolve,'Q after userdefined element assembly')
       call writeObject (RFk,dbgSolve,'RF after userdefined element assembly')
    end if
#endif

    call addInExternalForces (Qk, RFk, mech%forces, sam, ierr)
    if (istep < int(maxInt_p,i8)) then
       call addInFNVwaveForces (Qk, RFk, mech%env%Tsea, sam, int(istep)+1, ierr)
    end if
#ifdef FT_DEBUG
    if (iprint == -99 .and. size(mech%forces) > 0) then
       call writeObject (Qk ,dbgSolve,'Q after external force assembly')
       call writeObject (RFk,dbgSolve,'RF after external force assembly')
    end if
#endif

    haveSprings = size(mech%axialSprings) > 0
    do i = 1, size(mech%axialSprings)
       call addInSpringForces (FSk, RFk, mech%axialSprings(i), sam, ierr)
       if (.not. associated(mech%axialSprings(i)%spr(1)%p%length0Engine)) then
          call addInSpringForces (FDk, RFk, mech%axialSprings(i), sam, &
               &                  mech%axialSprings(i)%alpha2 + alpha2, ierr)
       end if
    end do
    do i = 1, size(mech%joints)
       if (associated(mech%joints(i)%springEl)) then
          haveSprings = .true.
          call addInSpringForces (FSk, RFk, mech%joints(i)%springEl, sam, ierr)
       end if
       if (associated(mech%joints(i)%sprFric)) then
          haveSprings = .true.
          call addInSpringForces (FSk, RFk, mech%joints(i)%sprFric, sam, ierr)
       end if
    end do

#ifdef FT_DEBUG
    if (iprint == -99 .and. haveSprings) then
       call writeObject (FSk,dbgSolve,'FS after spring assembly')
       call writeObject (RFk,dbgSolve,'RF after spring assembly')
    end if
#endif

    do i = 1, size(mech%dampers)
       call addInDamperForces (FDk, RFk, mech%dampers(i), sam, ierr)
    end do
#ifdef FT_DEBUG
    if (iprint == -99 .and. size(mech%dampers) > 0) then
       call writeObject (FDk,dbgSolve,'FD after damper assembly')
       call writeObject (RFk,dbgSolve,'RF after damper assembly')
    end if
#endif

    do i = 1, size(mech%bElems)
       call addInBushingElementForces (FSk, RFk, mech%bElems(i), sam, ierr, &
            &                          addSprings=.true.)
       call addInBushingElementForces (FDk, RFk, mech%bElems(i), sam, ierr, &
            &                          addDampers=.true.)
    end do
#ifdef FT_DEBUG
    if (iprint == -99 .and. size(mech%bElems) > 0) then
       call writeObject (FSk,dbgSolve,'FS after bushing element assembly')
       call writeObject (FDk,dbgSolve,'FD after bushing element assembly')
       call writeObject (RFk,dbgSolve,'RF after bushing element assembly')
    end if
#endif

    do i = 1, size(mech%cElems)
       call addInContactElementForces (FSk, RFk, mech%cElems(i), sam, ierr, &
            &                          addSprings=.true.)
       call addInContactElementForces (FDk, RFk, mech%cElems(i), sam, ierr, &
            &                          addDampers=.true.)
    end do
#ifdef FT_DEBUG
    if (iprint == -99 .and. size(mech%cElems) > 0) then
       call writeObject (FSk,dbgSolve,'FS after contact element assembly')
       call writeObject (FDk,dbgSolve,'FD after contact element assembly')
       call writeObject (RFk,dbgSolve,'RF after contact element assembly')
    end if
#endif

    do i = 1, size(mech%masses)
       call addInInertiaForces (FIk, RFk, mech%masses(i), sam, ierr)
       call addInGravitationForces (Qk, RFk, mech%masses(i), mech%gravity, &
            &                       sam, ierr)
    end do
#ifdef FT_DEBUG
    if (iprint == -99 .and. size(mech%masses) > 0) then
       call writeObject (FIk,dbgSolve,'FI after additional mass assembly')
       call writeObject (Qk ,dbgSolve,'Q after additional mass assembly')
       call writeObject (RFk,dbgSolve,'RF after additional mass assembly')
    end if
#endif

    do i = 1, size(mech%tires)
       call addInTireInertiaForces (FIk, RFk, mech%tires(i), sam, ierr)
       call addInTireForces (Qk, RFk, mech%tires(i), mech%gravity, sam, ierr)
    end do
#ifdef FT_DEBUG
    if (iprint == -99 .and. size(mech%tires) > 0) then
       call writeObject (FIk,dbgSolve,'FI after tire assembly')
       call writeObject (Qk ,dbgSolve,'Q after tire assembly')
       call writeObject (RFk,dbgSolve,'RF after tire assembly')
    end if
#endif

    if (associated(mech%turbine)) then
       call addInAeroForces (Qk, RFk, mech%turbine, sam, ierr)
#ifdef FT_DEBUG
       if (iprint == -99) then
          call writeObject (Qk ,dbgSolve,'Q after aero force assembly')
          call writeObject (RFk,dbgSolve,'RF after aero force assembly')
       end if
#endif
    end if

    if (size(mech%frictions) > 0) then
       do i = 1, size(mech%joints)
          call findForceInFriction (mech%joints(i), sam, RFk, Qk, Fsk, FDk, FIk)
          call addInFrictionForces (mech%joints(i), sam, Qk, RFk, ierr)
       end do
       do i = 1, size(mech%cElems)
          !! Temporarily removed until findForceInFriction is implemented
          !call findForceInFriction (mech%cElems(i), sam, RFk, Qk, Fsk, FDk,FIk)
          call addInContactElementForces (Qk, RFk, mech%cElems(i), sam, ierr, &
               &                          addFriction=.true.)
       end do
#ifdef FT_DEBUG
       if (iprint == -99) then
          call writeObject (Qk ,dbgSolve,'Q after friction assembly')
          call writeObject (RFk,dbgSolve,'RF after friction assembly')
       end if
#endif
    end if

    call stopTimer (asm_p)
    if (ierr < 0) call reportError (debugFileOnly_p,'GetForceVectors')

  end subroutine GetForceVectors


  !!============================================================================
  !> @brief Builds the quasti-static system force vectors.
  !>
  !> @param[out] FSk Internal stiffness forces
  !> @param[out] Qk External forces
  !> @param[out] RFk Reaction forces
  !> @param[in] sam Data for managing system matrix assembly
  !> @param mech Mechanism components of the model
  !> @param[out] ierr Error flag
  !>
  !> @details this subroutine puts static forces from all mechanism components
  !> into the associated system force vectors.
  !>
  !> @author Knut Morten Okstad                                   @date Jun 2002

  subroutine GetStaticForceVectors (FSk, Qk, RFk, sam, mech, ierr)

    use SamModule                   , only : SamType, dp
    use MechanismTypeModule         , only : MechanismType
    use TireTypeModule              , only : STI_p
    use SupElRoutinesModule         , only : addInStaticSupForces
    use UserdefElRoutinesModule     , only : addInStaticUDEForces
    use ForceRoutinesModule         , only : addInExternalForces
    use SpringRoutinesModule        , only : addInSpringForces
    use BushingElementRoutinesModule, only : addInBushingElementForces
    use ContactElementRoutinesModule, only : addInContactElementForces
    use FrictionRoutinesModule      , only : addInFrictionForces
    use FrictionRoutinesModule      , only : findForceInFriction
    use MassRoutinesModule          , only : addInGravitationForces
    use TireRoutinesModule          , only : addInTireForces
    use profilerModule              , only : startTimer, stopTimer, asm_p
    use reportErrorModule           , only : reportError, debugFileOnly_p
#ifdef FT_DEBUG
    use dbgUnitsModule        , only : dbgSolve
    use manipMatrixModule     , only : WriteObject
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getint
#endif
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_intValue

    real(dp)           , intent(out)   :: FSk(:), Qk(:), RFk(:)
    type(SamType)      , intent(in)    :: sam
    type(MechanismType), intent(inout) :: mech
    integer            , intent(out)   :: ierr

    !! Local variables
    integer :: i
    logical :: haveSprings

#ifdef FT_DEBUG
    integer :: iprint
    if (dbgSolve > 0) then
       call ffa_cmdlinearg_getint ('debug',iprint)
       write(dbgSolve,"(/'   --- in GetStaticForceVectors')")
    else
       iprint = 0
    end if
#endif

    !! --- Logic section ---

    call startTimer (asm_p)

    FSk = 0.0_dp
    Qk  = 0.0_dp
    RFk = 0.0_dp

    call addInStaticSupForces (sam, mech%sups, FSk, Qk, RFk, ierr)
#ifdef FT_DEBUG
    if (iprint == -99 .and. size(mech%sups) > 0) then
       call writeObject (FSk,dbgSolve,'FS after superelement assembly')
       call writeObject (Qk ,dbgSolve,'Q after superelement assembly')
       call writeObject (RFk,dbgSolve,'RF after superelement assembly')
    end if
#endif

    call addInStaticUDEForces (sam, mech%elms, FSk, Qk, RFk, ierr)
#ifdef FT_DEBUG
    if (iprint == -99 .and. size(mech%elms) > 0) then
       call writeObject (FSk,dbgSolve,'FS after userdefined element assembly')
       call writeObject (Qk ,dbgSolve,'Q after userdefined element assembly')
       call writeObject (RFk,dbgSolve,'RF after userdefined element assembly')
    end if
#endif

    call addInExternalForces (Qk, RFk, mech%forces, sam, ierr)
#ifdef FT_DEBUG
    if (iprint == -99 .and. size(mech%forces) > 0) then
       call writeObject (Qk ,dbgSolve,'Q after external force assembly')
       call writeObject (RFk,dbgSolve,'RF after external force assembly')
    end if
#endif

    haveSprings = size(mech%axialSprings) > 0
    do i = 1, size(mech%axialSprings)
       call addInSpringForces (FSk, RFk, mech%axialSprings(i), sam, ierr)
    end do
    do i = 1, size(mech%joints)
       if (associated(mech%joints(i)%springEl)) then
          haveSprings = .true.
          call addInSpringForces (FSk, RFk, mech%joints(i)%springEl, sam, ierr)
       end if
       if (associated(mech%joints(i)%sprFric)) then
          haveSprings = .true.
          call addInSpringForces (FSk, RFk, mech%joints(i)%sprFric, sam, ierr)
       end if
    end do
#ifdef FT_DEBUG
    if (iprint == -99 .and. haveSprings) then
       call writeObject (FSk,dbgSolve,'FS after spring assembly')
       call writeObject (RFk,dbgSolve,'RF after spring assembly')
    end if
#endif

    do i = 1, size(mech%bElems)
       call addInBushingElementForces (FSk, RFk, mech%bElems(i), sam, ierr, &
            &                          addSprings=.true.)
    end do
#ifdef FT_DEBUG
    if (iprint == -99 .and. size(mech%bElems) > 0) then
       call writeObject (FSk,dbgSolve,'FS after bushing element assembly')
       call writeObject (RFk,dbgSolve,'RF after bushing element assembly')
    end if
#endif

    do i = 1, size(mech%cElems)
       call addInContactElementForces (FSk, RFk, mech%cElems(i), sam, ierr, &
            &                          addSprings=.true.)
       call addInContactElementForces (Qk, RFk, mech%cElems(i), sam, ierr, &
            &                          addFriction=.true.)
    end do
#ifdef FT_DEBUG
    if (iprint == -99 .and. size(mech%cElems) > 0) then
       call writeObject (FSk,dbgSolve,'FS after contact element assembly')
       call writeObject (RFk,dbgSolve,'RF after contact element assembly')
    end if
#endif

    do i = 1, size(mech%masses)
       call addInGravitationForces (Qk, RFk, mech%masses(i), mech%gravity, &
            &                       sam, ierr)
    end do
#ifdef FT_DEBUG
    if (iprint == -99 .and. size(mech%masses) > 0) then
       call writeObject (Qk ,dbgSolve,'Q after additional mass assembly')
       call writeObject (RFk,dbgSolve,'RF after additional mass assembly')
    end if
#endif

    if (ffa_cmdlinearg_intValue('tireStaticEquil') > 0) then
       do i = 1, size(mech%tires)
          if (mech%tires(i)%api == STI_p) then
             call addInTireForces (Qk, RFk, mech%tires(i),mech%gravity,sam,ierr)
          end if
       end do
#ifdef FT_DEBUG
       if (iprint == -99 .and. size(mech%tires) > 0) then
          call writeObject (Qk,dbgSolve,'Q after tire assembly')
          call writeObject (RFk,dbgSolve,'RF after tire assembly')
       end if
#endif
    end if

    if (size(mech%frictions) > 0) then
       do i = 1, size(mech%joints)
          call findForceInFriction (mech%joints(i), sam, RFk, Qk, Fsk)
          call addInFrictionForces (mech%joints(i), sam, Qk, RFk, ierr)
       end do
       do i = 1, size(mech%cElems)
          !! Temporarily removed until findForceInFriction is implemented
          !call findForceInFriction (mech%cElems(i), sam, RFk, Qk, Fsk)
          call addInContactElementForces (Qk, RFk, mech%cElems(i), sam, ierr, &
               &                          addFriction=.true.)
       end do
#ifdef FT_DEBUG
       if (iprint == -99) then
          call writeObject (Qk ,dbgSolve,'Q after friction assembly')
          call writeObject (RFk,dbgSolve,'RF after friction assembly')
       end if
#endif
    end if

    call stopTimer (asm_p)
    if (ierr < 0) call reportError (debugFileOnly_p,'GetStaticForceVectors')

  end subroutine GetStaticForceVectors

end module AddInSysModule
