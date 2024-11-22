!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file supElRoutinesModule.f90
!>
!> @brief Subroutines for superelement calculations.

!!==============================================================================
!> @brief Module with subroutines for superelement calculations.
!>
!> @details This module contains a set of subroutines for performing various
!> computation tasks on the supeltypemodule::supeltype objects in the model
!> during the dynamic or quasi-static simulation.

module SupElRoutinesModule

  implicit none

  private :: calcBuoyancyForces, calcMorisonForces
  private :: updateSupElLoad, scaledMatmul


contains

  !!============================================================================
  !> @brief Extracts local velocities and accelerations for the superelements.
  !>
  !> @param sups All superelements in the model
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[in] velGlobal Global velocity vector
  !> @param[in] accGlobal Global acceleration vector
  !>
  !> @details This subroutine extracts the local velocities and accelerations
  !> (@a uld and @a uldd, respectively) from the given global velocity and
  !> acceleration vectors for all superelements.
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 16 Nov 2001

  subroutine SetSupElsVelAcc (sups,sam,velGlobal,accGlobal)

    use SupElTypeModule   , only : SupElType, dp
    use SamModule         , only : SamType
    use AsmExtensionModule, only : csExtractElv

    type(SupElType), intent(inout) :: sups(:)
    type(SamType)  , intent(in)    :: sam
    real(dp)       , intent(in)    :: velGlobal(:), accGlobal(:)

    !! Local variables
    integer  :: i, j, n
    real(dp) :: trTr(3,3)

    !! --- Logic section ---

    do i = 1, size(sups)

       call csExtractElv (sam,sups(i)%samElNum,sups(i)%uld ,velGlobal)
       call csExtractElv (sam,sups(i)%samElNum,sups(i)%uldd,accGlobal)

       !! Transform to the corotated superelement coordinate system
       trTr = transpose(sups(i)%supTr(:,1:3))

       do j = 1, sups(i)%nExtNods
          if (sups(i)%triads(j)%p%nDOFs < 3) cycle

          n = sups(i)%triads(j)%firstDOF
          sups(i)%uld (n:n+2) = matmul(trTr,sups(i)%uld (n:n+2))
          sups(i)%uldd(n:n+2) = matmul(trTr,sups(i)%uldd(n:n+2))

          if (sups(i)%triads(j)%p%nDOFs < 6) cycle

          sups(i)%uld (n+3:n+5) = matmul(trTr,sups(i)%uld (n+3:n+5))
          sups(i)%uldd(n+3:n+5) = matmul(trTr,sups(i)%uldd(n+3:n+5))

       end do
    end do

  end subroutine SetSupElsVelAcc


  !!============================================================================
  !> @brief Calculates deformational velocities and accelerations.
  !>
  !> @param sup The superelement to calculate velocity/acceleration for
  !> @param[in] beta Newmark time integration parameter
  !> @param[in] gamma Newmark time integration parameter
  !> @param[in] h Current time increment size
  !>
  !> @details This subroutine integrates the deformational velocities and
  !> accelerations based on the Newmark integration parameters and the state
  !> at the previous time step.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 6 Jul 2006

  subroutine BuildFinitVelAcc (sup,beta,gamma,h)

    use SupElTypeModule, only : SupElType, dp

    type(SupElType), intent(inout) :: sup
    real(dp)       , intent(in)    :: beta, gamma, h

    !! Local variables
    integer  :: nDOF
    real(dp) :: c1, c2, c3

    !! --- Logic section ---

    !! Find the deformational displacement increment
    nDOF = size(sup%vld)
    call DCOPY (nDOF,sup%finit(1),1,sup%vld(1),1)
    call DAXPY (nDOF,-1.0_dp,sup%finitPrev(1),1,sup%vld(1),1)

    if (associated(sup%vldd)) then ! Newmark integration

       !! Integrate the deformational acceleration
       c1 = 1.0_dp/(beta*h*h)
       c2 = 1.0_dp/(beta*h)
       c3 = 1.0_dp - 0.5_dp/beta
       !! vldd = c1*vld - c2*vldPrev + c3*vlddPrev
       call DAXPY (nDOF, c1,sup%vld(1)     ,1,sup%vldd(1),1)
       call DAXPY (nDOF,-c2,sup%vldPrev(1) ,1,sup%vldd(1),1)
       call DAXPY (nDOF, c3,sup%vlddPrev(1),1,sup%vldd(1),1)

       !! Integrate the deformational velocity
       c1 = gamma/(beta*h)
       c2 = 1.0_dp - gamma/beta
       c3 = 0.5_dp*h*gamma/beta - h
       !! vld = c1*vld + c2*vldPrev - c3*vlddPrev
       call DSCAL (nDOF, c1,sup%vld(1)     ,1)
       call DAXPY (nDOF, c2,sup%vldPrev(1) ,1,sup%vld(1),1)
       call DAXPY (nDOF,-c3,sup%vlddPrev(1),1,sup%vld(1),1)

    else ! Use first order approximation (this should be R4.1 compatible)

       !! vld = vld / h
       call DSCAL (nDOF,1.0_dp/h,sup%vld(1),1)

    end if

  end subroutine BuildFinitVelAcc


  !!============================================================================
  !> @brief Increments the generalized DOFs for all superelements.
  !>
  !> @param sups All superelements in the model
  !> @param[in] solinc Current global solution increment
  !> @param[in] useTotalInc If .true., update from previous solution state
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 16 Nov 2001

  subroutine IncSupElsGenDofs (sups,solinc,useTotalInc)

    use SupElTypeModule, only : SupElType, dp

    type(SupElType)  , intent(inout) :: sups(:)
    real(dp)         , intent(in)    :: solinc(:)
    logical, optional, intent(in)    :: useTotalInc

    !! Local variables
    logical :: fromPrevious
    integer :: i, j, n

    !! --- Logic section ---

    if (present(useTotalInc)) then
       fromPrevious = useTotalInc
    else
       fromPrevious = .false.
    end if

    do i = 1, size(sups)
       if (associated(sups(i)%genDOFs)) then
          j = sups(i)%genDOFs%sysDOF
          n = j+sups(i)%genDOFs%nDOFs-1
          if (fromPrevious) then
             sups(i)%genDOFs%ur = sups(i)%genDOFs%urPrev + solinc(j:n)
          else
             sups(i)%genDOFs%ur = sups(i)%genDOFs%ur + solinc(j:n)
          end if
       end if
    end do

  end subroutine IncSupElsGenDofs


  !!============================================================================
  !> @brief Updates all superelements in the model based on the computed state.
  !>
  !> @param sups All superelements in the model
  !> @param supLoads All superelement loads in the model
  !> @param env Environmental data
  !> @param[in] beta Newmark time integration parameter
  !> @param[in] gamma Newmark time integration parameter
  !> @param[in] time Current physical time
  !> @param[in] timeStep Current time increment size
  !> @param[in] istep Time increment counter
  !> @param[in] iter Iteration counter
  !> @param[in] newPositions If .true., the position variables have been updated
  !> @param ierr Error flag
  !>
  !> @details This subroutine updates the superelement objects:
  !> - the deflection vector, @a Finit, is established
  !> - superelement forces are calculated and transformed to system directions
  !> - the total nodal forces are calculated in system directions
  !>
  !> @callgraph @callergraph
  !>
  !> @author Karl Erik Thoresen                                   @date Sep 1999
  !> @author Knut Morten Okstad                                   @date Nov 2001
  !> @author Knut Morten Okstad                                   @date Feb 2002
  !> @author Knut Morten Okstad                                   @date Apr 2008
  !> @author Knut Morten Okstad                                   @date Aug 2011

  subroutine UpdateSupEls (sups,supLoads,env,beta,gamma,time,timeStep, &
       &                   istep,iter,newPositions,ierr)

    use SupElTypeModule           , only : SupElType, dp
    use SupElLoadTypeModule       , only : SupElLoadType
    use EnvironmentTypeModule     , only : EnvironmentType
    use SupElTypeModule           , only : UpdateSupElCorot, BuildFinit, IsBeam
    use SupElTypeModule           , only : GetRsupToSys
    use TriadTypeModule           , only : UpdateNodeForce
    use massMatrixCorrectionModule, only : mmcGetMassTorqueCorrection
    use HydrodynamicsModule       , only : getDragForces, getDiffractionForces
    use EngineRoutinesModule      , only : isPredictorStep
#ifdef FT_DEBUG
    use IdTypeModule              , only : getId
    use manipMatrixModule         , only : writeObject
    use dbgUnitsModule            , only : dbgSolve
#endif
    use dbgUnitsModule            , only : dbgShadowPos
    use reportErrorModule         , only : reportError, debugFileOnly_p
    use FFaCmdLineArgInterface    , only : ffa_cmdlinearg_getint

    type(SupElType)      , intent(inout) :: sups(:)
    type(SupElLoadType)  , intent(inout) :: supLoads(:)
    type(EnvironmentType), intent(inout) :: env
    real(dp)             , intent(in)    :: beta, gamma, time, timeStep
    integer              , intent(in)    :: istep, iter
    logical              , intent(in)    :: newPositions
    integer              , intent(inout) :: ierr

    !! Local variables
    integer  :: i, j, n, m, lerr, stat, dragForm
    real(dp) :: RsupToSys(3,3), sclD
    real(dp), target  :: dummy(6)
    real(dp), pointer :: f(:,:), fn(:)

    !! --- Logic section ---

#ifdef FT_DEBUG
    write(dbgSolve,600) time,istep,iter
600 format(/'   --- in UpdateSupEls, time =',1PE12.5,' istep =',I6,' iter =',I3)
#endif

    call ffa_cmdlinearg_getint ('dragForm',dragForm)

    !! Scaling factor for Morison drag forces
    if (isPredictorStep) then
       sclD = 1.0_dp
    else
       sclD = 0.5_dp ! because the damping matrix Cdmat is linear wrt. velocity
    end if

    lerr = ierr
    do i = 1, size(sups)
       m = size(sups(i)%KmMat,1)
       n = size(sups(i)%KmMat,2)

       if (newPositions) then
          stat = 0

          !! Update the co-rotated superelement coordinate system
          call UpdateSupElCorot (sups(i),dbgShadowPos,stat)
          ierr = ierr + stat
          if (stat < 0) cycle

          !! Calculate superelement deformation vector in co-rotated system
          call BuildFinit (sups(i))
       end if

       if (associated(sups(i)%Mmat)) then
          !! Calculate the inertia forces in the superelement, Fi = M*uldd
          call scaledMatmul (m,n,sups(i)%massScl(1), &
               &             sups(i)%Mmat,sups(i)%uldd,sups(i)%Fi)
       else
          call DCOPY (m,0.0_dp,0,sups(i)%Fi(1),1)
       end if

       if (associated(sups(i)%Cmat)) then
          !! Calculate the damping forces in the superelement, Fd = C*uld.
          !! We here assume that only stiffness-proportional damping is used,
          !! unless the stiffness-scaling factor is equal to one.
          call scaledMatmul (m,n,sups(i)%stifScl(1), &
               &             sups(i)%Cmat,sups(i)%uld,sups(i)%Fd)
       else
          call DCOPY (m,0.0_dp,0,sups(i)%Fd(1),1)
       end if

       if (timeStep > 0.0_dp .and. associated(sups(i)%vld)) then
          j = size(sups(i)%vld)

          !! Calculate current deformational velocity for this superelement
          call BuildFinitVelAcc (sups(i),beta,gamma,timeStep)

          !! Filter out rigid body modes from the stiffness-proportional part
          !! Fd = Fd - alpha2*K*(uld-vld)
          !! using sups(i)%Q and sups(i)%Fs as temporaries below
          call DCOPY (j,sups(i)%uld(1),1,sups(i)%Q(1),1)
          call DAXPY (j,-1.0_dp,sups(i)%vld(1),1,sups(i)%Q(1),1)
          call scaledMatmul (j,j,sups(i)%stifScl(1), &
               &             sups(i)%KmMat,sups(i)%Q,sups(i)%Fs,m)
          call DAXPY (j,-sups(i)%kDmpFactor,sups(i)%Fs(1),1,sups(i)%Fd(1),1)
       end if

       !! Calculate the stiffness forces in the superelement, Fs = K*Finit
       call scaledMatmul (m,n,sups(i)%stifScl(1), &
            &             sups(i)%KmMat,sups(i)%Finit,sups(i)%Fs)

       if (associated(sups(i)%fg)) then
          !! Calculate the gravitational forces which is an external force
          call scaledMatmul (m,3,sups(i)%massScl(1),sups(i)%fg, &
               &             matmul(env%gravity,sups(i)%supTr(:,1:3)),sups(i)%Q)
       else
          call DCOPY (m,0.0_dp,0,sups(i)%Q(1),1)
       end if

       if (env%rhow > 0.0_dp .and. associated(sups(i)%hydyn)) then
          if (IsBeam(sups(i))) then
             !! Calculate Morison forces and add to the external forces
             call calcMorisonForces (sups(i),env,time,istep,abs(iter),ierr)
             if (dragForm < 1) then
                !! Add damping and inertia forces due to the structural motion
                !! Not using BLAS here, since these are 12x12 matrices only
                sups(i)%FD = sups(i)%FD + matmul(sups(i)%Cdmat,sups(i)%uld)*sclD
                sups(i)%FI = sups(i)%FI + matmul(sups(i)%Mamat,sups(i)%uldd)
             end if
          else
             !! Calculate buoyancy forces and add to the external forces
             call calcBuoyancyForces (sups(i),env,time,abs(iter),ierr)
             if (.not. isPredictorStep) then
                !! Calculate drag and slam forces and add to the damping forces
                call getDragForces (sups(i)%FD,sups(i),sups(i)%hydyn,env, &
                     &              time,timeStep,ierr)
             end if
          end if
          !! Calculate diffraction forces and add to the external forces
          call getDiffractionForces (sups(i)%Q,sups(i),sups(i)%hydyn, &
               &                     env,time,ierr)
       end if

       !! Calculate the mass matrix correction. The correction has no terms
       !! in the mass matrix and is therefore added as an external force only.
       select case (sups(i)%massCorrFlag)
       case (1)
          sups(i)%Q  = sups(i)%Q  + mmcGetMassTorqueCorrection(sups(i))
       case (2)
          sups(i)%FD = sups(i)%FD - mmcGetMassTorqueCorrection(sups(i))
       case (3)
          sups(i)%FI = sups(i)%FI - mmcGetMassTorqueCorrection(sups(i))
       end select

    end do
    if (ierr < lerr) goto 900

    !! Calculate superelement loads and add to the external force vector
    do i = 1, size(supLoads)
       call updateSupElLoad (supLoads(i)%sup%Q,supLoads(i)%sup%S, &
            &                supLoads(i),ierr)
       if (ierr < lerr) goto 900
    end do

#ifdef FT_DEBUG
    do i = 1, size(sups)
       write(dbgSolve,"()")
       call writeObject (sups(i)%supTr,dbgSolve,'Updated coordinate system '// &
            &            'for Superelement'//trim(getId(sups(i)%id)))
       call writeObject (sups(i)%finit,dbgSolve,'Nodal deformations')
       if (timeStep > 0.0_dp .and. associated(sups(i)%vld)) then
          call writeObject (sups(i)%vld,dbgSolve,'Nodal deformational velocity')
       end if
       call writeObject (sups(i)%uld,dbgSolve,'Nodal velocity')
       call writeObject (sups(i)%uldd,dbgSolve,'Nodal acceleration')
       call writeObject (sups(i)%FS,dbgSolve,'Stiffness forces')
       call writeObject (sups(i)%FD,dbgSolve,'Damping forces')
       call writeObject (sups(i)%FI,dbgSolve,'Inertia forces')
       call writeObject (sups(i)%Q,dbgSolve,'External forces')
    end do
#endif

    !! Transform the superelement forces to the system coordinates
    !! and update the sectional forces as well as the total nodal forces

    do i = 1, size(sups)

       do j = 1, sups(i)%nExtNods
          if (sups(i)%triads(j)%p%nDOFs < 3) cycle

          n = sups(i)%triads(j)%firstDOF
          m = n-1 + sups(i)%triads(j)%p%nDOFs
          RsupToSys = GetRsupToSys(sups(i),j)
          sups(i)%FS(n:n+2) = matmul(RsupToSys,sups(i)%FS(n:n+2))
          sups(i)%FD(n:n+2) = matmul(RsupToSys,sups(i)%FD(n:n+2))
          sups(i)%FI(n:n+2) = matmul(RsupToSys,sups(i)%FI(n:n+2))
          sups(i)%Q (n:n+2) = matmul(RsupToSys,sups(i)%Q (n:n+2))
          if (sups(i)%triads(j)%p%nDOFs >= 6) then
             sups(i)%FS(n+3:n+5) = matmul(RsupToSys,sups(i)%FS(n+3:n+5))
             sups(i)%FD(n+3:n+5) = matmul(RsupToSys,sups(i)%FD(n+3:n+5))
             sups(i)%FI(n+3:n+5) = matmul(RsupToSys,sups(i)%FI(n+3:n+5))
             sups(i)%Q (n+3:n+5) = matmul(RsupToSys,sups(i)%Q (n+3:n+5))
          end if

          if (associated(sups(i)%triads(j)%p%supForce)) then
             fn => sups(i)%triads(j)%p%supForce
          else
             fn => dummy(1:sups(i)%triads(j)%p%nDOFs)
          end if
          fn = sups(i)%FS(n:m) + sups(i)%FD(n:m) + sups(i)%FI(n:m) &
               -sups(i)%Q(n:m)

          call updateNodeForce (sups(i)%triads(j)%p,fn)

          if (associated(sups(i)%triads(j)%p%force)) then
             m = n + sups(i)%triads(j)%p%nDOFs - 1
             f => sups(i)%triads(j)%p%force
             f(:,1) = f(:,1) + sups(i)%FS(n:m)
             f(:,2) = f(:,2) + sups(i)%FD(n:m)
             f(:,3) = f(:,3) + sups(i)%FI(n:m)
             f(:,4) = f(:,4) + sups(i)%Q (n:m)
          end if

       end do

    end do

    return

900 call reportError (debugFileOnly_p,'UpdateSupEls')

  end subroutine UpdateSupEls


  !!============================================================================
  !> @brief Updates all superelements in the model based on the computed state.
  !>
  !> @param sups All superelements in the model
  !> @param supLoads All superelement loads in the model
  !> @param env Environmental data
  !> @param[in] time Current physical time
  !> @param[in] iter Iteration counter
  !> @param[in] linInc Option for linear (1) or incremental (2) deformations
  !> @param ierr Error flag
  !>
  !> @details This subroutine updates the superelement objects statically.
  !> It is used in place of supelroutinesmodule::updatesupels
  !> in a quasi-static simulation.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 5 June 2002

  subroutine UpdateSupElsStatic (sups,supLoads,env,time,iter,linInc,ierr)

    use SupElTypeModule      , only : SupElType, dp
    use SupElLoadTypeModule  , only : SupElLoadType
    use EnvironmentTypeModule, only : EnvironmentType
    use SupElTypeModule      , only : UpdateSupElCorot, BuildFinit, GetRsupToSys
    use TriadTypeModule      , only : UpdateNodeForce
#ifdef FT_DEBUG
    use IdTypeModule         , only : getId
    use manipMatrixModule    , only : writeObject
    use dbgUnitsModule       , only : dbgSolve
#endif
    use dbgUnitsModule       , only : dbgShadowPos
    use reportErrorModule    , only : reportError, debugFileOnly_p

    type(SupElType)      , intent(inout) :: sups(:)
    type(SupElLoadType)  , intent(inout) :: supLoads(:)
    type(EnvironmentType), intent(inout) :: env
    real(dp)             , intent(in)    :: time
    integer              , intent(in)    :: iter
    integer, optional    , intent(in)    :: linInc
    integer              , intent(inout) :: ierr

    !! Local variables
    logical  :: useCorotUpd, useFinitInc
    integer  :: i, j, n, m, lerr, stat
    real(dp) :: RsupToSys(3,3)
    real(dp), target  :: dummy(6)
    real(dp), pointer :: f(:,:), fn(:)

    !! --- Logic section ---

#ifdef FT_DEBUG
    write(dbgSolve,600) time,iter
600 format(/'   --- in UpdateSupElsStatic, time =',1PE12.5,' iter =',I3)
#endif

    if (present(linInc)) then
       useCorotUpd = linInc == 0
       useFinitInc = linInc > 1
    else
       useCorotUpd = .true.
       useFinitInc = .false.
    end if

    lerr = ierr
    do i = 1, size(sups)
       m = size(sups(i)%KmMat,1)
       n = size(sups(i)%KmMat,2)

       if (useCorotUpd) then
          stat = 0 ! Update the co-rotated superelement coordinate system
          call UpdateSupElCorot (sups(i),dbgShadowPos,stat)
          ierr = ierr + stat
          if (stat < 0) cycle
       end if

       !! Calculate current superelement deformation vector in co-rotated system
       call BuildFinit (sups(i))
       if (useFinitInc) then ! Use incremental deformation
          sups(i)%finit = sups(i)%finit - sups(i)%finitPrev
       end if

       !! The inertia and damping forces are zero in a statics solution
       call DCOPY (m,0.0_dp,0,sups(i)%Fi(1),1)
       call DCOPY (m,0.0_dp,0,sups(i)%Fd(1),1)

       !! Calculate the stiffness forces in the superelement, Fs = K*Finit
       call scaledMatmul (m,n,sups(i)%stifScl(1), &
            &             sups(i)%KmMat,sups(i)%Finit,sups(i)%Fs)

       if (associated(sups(i)%fg) .and. .not. useFinitInc) then
          !! Calculate the gravitational forces which is an external force
          call scaledMatmul (m,3,sups(i)%stifScl(1),sups(i)%fg, &
               &             matmul(env%gravity,sups(i)%supTr(:,1:3)),sups(i)%Q)
       else
          call DCOPY (m,0.0_dp,0,sups(i)%Q(1),1)
       end if

       if (env%rhow > 0.0_dp .and. associated(sups(i)%hydyn)) then
          !! Calculate buoyancy forces and add to the external forces
          call calcBuoyancyForces (sups(i),env,time,iter,ierr)
       end if

    end do
    if (ierr < lerr) goto 900

    !! Calculate superelement loads and add to the external force vector
    do i = 1, size(supLoads)
       call updateSupElLoad (supLoads(i)%sup%Q,supLoads(i)%sup%S, &
            &                supLoads(i),ierr)
       if (ierr < lerr) goto 900
    end do

#ifdef FT_DEBUG
    do i = 1, size(sups)
       if (useCorotUpd) then
          write(dbgSolve,"()")
          call writeObject (sups(i)%supTr,dbgSolve,'Updated co-rotated '// &
              &            'system for Superelement'//trim(getId(sups(i)%id)))
       else
          write(dbgSolve,"(/7X,'for Superelement',A/)") trim(getId(sups(i)%id))
       end if
       call writeObject (sups(i)%finit,dbgSolve,'Nodal deformations')
       call writeObject (sups(i)%FS,dbgSolve,'Resulting stiffness forces')
       call writeObject (sups(i)%Q,dbgSolve,'External forces')
    end do
#endif

    !! Transform the superelement forces to the system coordinates
    !! and update the sectional forces as well as the total nodal forces

    do i = 1, size(sups)

       do j = 1, sups(i)%nExtNods
          if (sups(i)%triads(j)%p%nDOFs < 3) cycle

          n = sups(i)%triads(j)%firstDOF
          m = n-1 + sups(i)%triads(j)%p%nDOFs
          RsupToSys = GetRsupToSys(sups(i),j)
          sups(i)%FS(n:n+2) = matmul(RsupToSys,sups(i)%FS(n:n+2))
          sups(i)%Q (n:n+2) = matmul(RsupToSys,sups(i)%Q (n:n+2))
          if (sups(i)%triads(j)%p%nDOFs >= 6) then
             sups(i)%FS(n+3:n+5) = matmul(RsupToSys,sups(i)%FS(n+3:n+5))
             sups(i)%Q (n+3:n+5) = matmul(RsupToSys,sups(i)%Q (n+3:n+5))
          end if

          if (associated(sups(i)%triads(j)%p%supForce)) then
             fn => sups(i)%triads(j)%p%supForce
          else
             fn => dummy(1:sups(i)%triads(j)%p%nDOFs)
          end if
          fn = sups(i)%FS(n:m) - sups(i)%Q(n:m)
          call updateNodeForce (sups(i)%triads(j)%p,fn)

          if (associated(sups(i)%triads(j)%p%force)) then
             f => sups(i)%triads(j)%p%force
             f(:,1) = f(:,1) + sups(i)%FS(n:m)
             f(:,4) = f(:,4) + sups(i)%Q (n:m)
          end if

       end do

    end do

    return

900 call reportError (debugFileOnly_p,'UpdateSupElsStatic')

  end subroutine UpdateSupElsStatic


  !!============================================================================
  !> @brief Adds the superelement loads to the system external load vector.
  !>
  !> @param Q System external load vector
  !> @param S Normalized superelement load vectors
  !> @param supLoad The superelement load object to assemble
  !> @param ierr Error flag
  !>
  !> @details This subroutine updates the amplitude of the superelement load and
  !> adds its contribution to the external force vector of the global system.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten okstad
  !>
  !> @date 8 Apr 2008

  subroutine UpdateSupElLoad (Q,S,supLoad,ierr)

    use kindModule          , only : dp
    use SupElLoadTypeModule , only : SupElLoadType
    use IdTypeModule        , only : getId
    use EngineRoutinesModule, only : Evaluate, UpdateSensor
    use reportErrorModule   , only : reportError, error_p

    real(dp)           , intent(inout) :: Q(:)
    real(dp)           , intent(in)    :: S(:,:)
    type(SupElLoadType), intent(inout) :: supLoad
    integer            , intent(inout) :: ierr

    !! Local variables
    integer  :: lerr
    real(dp) :: xArg

    !! --- Logic section ---

    lerr = ierr
    if (.not. associated(supLoad%Engine)) then
       xArg = 0.0_dp
    else if (.not. associated(supLoad%Engine%args(1)%p)) then
       xArg = -supLoad%delay
    else
       call UpdateSensor (supLoad%Engine%args(1)%p,ierr)
       xArg = supLoad%Engine%args(1)%p%value - supLoad%delay
    end if

    !! Get the engine value for this load in current state
    supLoad%F = Evaluate(supLoad%Engine,supLoad%f1,supLoad%f0,ierr,xArg)
    if (ierr < lerr) then
       call reportError (error_p,'Failed to evaluate engine '// &
            &                    'for superelement load'//getId(supLoad%id), &
            &            addString='updateSupElLoad')
    else
       !! Update the external force vector, Q
       call DAXPY (size(Q),supLoad%F,S(1,supLoad%loadCase),1,Q(1),1)
    end if

  end subroutine UpdateSupElLoad


  !!============================================================================
  !> @brief Updates the superelement damping matrices.
  !>
  !> @param sups All superelements in the model
  !> @param engs All general functions in the model
  !> @param[in] alpha Global structural damping parameters
  !> @param[out] ierr Error flag
  !>
  !> @details Effective with time-dependent structural damping coefficients.
  !> This subroutine also updates the time-dependent stiffness- and mass
  !> scaling factors.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 1 Dec 2010

  subroutine UpdateSupElDamping (sups,engs,alpha,ierr)

    use SupElTypeModule     , only : SupElType, dp
    use FunctionTypeModule  , only : EngineType
    use EngineRoutinesModule, only : EngineValue
    use reportErrorModule   , only : reportError, debugFileOnly_p

    type(SupElType) , intent(inout) :: sups(:)
    type(EngineType), intent(inout) :: engs(:)
    real(dp)        , intent(in)    :: alpha(2)
    integer         , intent(out)   :: ierr

    !! Local variables
    integer             :: i, n
    logical             :: newAlpha, newDmp
    real(dp)            :: dmpFac
    real(dp), save      :: prevAlpha(2) = 0.0_dp
    real(dp), parameter :: eps_p = 1.0e-16_dp

    !! --- Logic section ---

    newAlpha  = any(abs(alpha-prevAlpha) > eps_p)
    prevAlpha = alpha

    ierr = 0
    do i = 1, size(sups)
       if (sups(i)%dmpSclIdx > 0) then
          sups(i)%dmpScl(2) = sups(i)%dmpScl(1)
          sups(i)%dmpScl(1) = EngineValue(engs(sups(i)%dmpSclIdx),ierr)
          newDmp = abs(sups(i)%dmpScl(1) - sups(i)%dmpScl(2)) > eps_p
       else
          newDmp = .false.
       end if
       if (sups(i)%massSclIdx > 0) then
          sups(i)%massScl(2) = sups(i)%massScl(1)
          sups(i)%massScl(1) = EngineValue(engs(sups(i)%massSclIdx),ierr)
       end if
       if (sups(i)%stifSclIdx > 0) then
          sups(i)%stifScl(2) = sups(i)%stifScl(1)
          sups(i)%stifScl(1) = EngineValue(engs(sups(i)%stifSclIdx),ierr)
       end if
       if (ierr == 0 .and. (newAlpha .or. newDmp)) then

          !! Update the damping coefficients and compute new damping matrix
          n = size(sups(i)%Cmat)
          call DCOPY (n,0.0_dp,0,sups(i)%Cmat(1,1),1)

          dmpFac = sups(i)%mDmp0*sups(i)%dmpScl(1) + alpha(1)
          if (abs(dmpFac) > eps_p) then
             sups(i)%mDmpFactor = dmpFac
             call DAXPY (n,dmpFac,sups(i)%Mmat(1,1),1,sups(i)%Cmat(1,1),1)
          else
             sups(i)%mDmpFactor = 0.0_dp
          end if

          dmpFac = sups(i)%kDmp0*sups(i)%dmpScl(1) + alpha(2)
          if (abs(dmpFac) > eps_p) then
             sups(i)%kDmpFactor = dmpFac
             call DAXPY (n,dmpFac,sups(i)%KmMat(1,1),1,sups(i)%Cmat(1,1),1)
          else
             sups(i)%kDmpFactor = 0.0_dp
          end if

       end if
    end do

    if (ierr < 0) call reportError (debugFileOnly_p,'UpdateSupElDamping')

  end subroutine UpdateSupElDamping


  !!============================================================================
  !> @brief Updates the current sea state.
  !>
  !> @param env Environmental data
  !> @param[in] triads All triads in the model
  !> @param[in] sups All superelements in the model
  !> @param[in] time Current physical time
  !> @param[in] iter Iteration counter
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Oct 2008

  subroutine UpdateSeaEnvironment (env,triads,sups,time,iter,ierr)

    use EnvironmentTypeModule , only : EnvironmentType, dp
    use TriadTypeModule       , only : TriadType
    use SupElTypeModule       , only : SupElType
    use HydrodynamicsModule   , only : initFluidMotions
    use FunctionTypeModule    , only : FunctionValue
    use EngineRoutinesModule  , only : EngineValue
    use reportErrorModule     , only : reportError, debugFileOnly_p
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_intValue

    type(EnvironmentType), intent(inout) :: env
    type(TriadType)      , intent(in)    :: triads(:)
    type(SupElType)      , intent(in)    :: sups(:)
    real(dp)             , intent(in)    :: time
    integer              , intent(in)    :: iter
    integer              , intent(out)   :: ierr

    !! --- Logic section ---

    ierr = 0
    if (env%rhow <= 0.0_dp) return

    if (iter < ffa_cmdlinearg_intValue('nHDupdat')) then
       call initFluidMotions (time,env,triads,sups,ierr)
    end if

    if (associated(env%waveFunc)) then
       env%seaLevel = env%sea0 + FunctionValue(env%waveFunc,time,ierr)
    end if

    if (associated(env%cScaleEng)) then
       env%currScale = EngineValue(env%cScaleEng,ierr)
    end if

    if (associated(env%wScaleEng)) then
       if (associated(env%wScaleEng,env%cScaleEng)) then
          env%seaScale = env%currScale
       else
          env%seaScale = EngineValue(env%wScaleEng,ierr)
       end if
    end if

    if (associated(env%hScaleEng)) then
       env%hdfScale = EngineValue(env%hScaleEng,ierr)
    end if

    if (ierr < 0) call reportError (debugFileOnly_p,'UpdateSeaEnvironment')

  end subroutine UpdateSeaEnvironment


  !!============================================================================
  !> @brief Calculates Morison force contributions for a two-noded beam element.
  !>
  !> @param sup The beam element to calculate Morison forces for
  !> @param env Environmental data
  !> @param[in] time Current simulation time
  !> @param[in] istep Time increment counter
  !> @param[in] iter Iteration counter
  !> @param ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 24 Sep 2019

  subroutine calcMorisonForces (sup,env,time,istep,iter,ierr)

    use SupElTypeModule       , only : SupElType, dp
    use EnvironmentTypeModule , only : EnvironmentType
    use HydrodynamicsModule   , only : getMorisonForces
    use IdTypeModule          , only : getId

    type(SupElType)      , intent(inout) :: sup
    type(EnvironmentType), intent(inout) :: env
    real(dp)             , intent(in)    :: time
    integer              , intent(in)    :: istep, iter
    integer              , intent(inout) :: ierr

    !! Local variables
    integer :: err

    !! --- Logic section ---

    call getMorisonForces ('Beam element'//getId(sup%id),sup%triads,sup%supTr, &
                           sup%uld,sup%uldd,sup%Q,sup%MaMat,sup%CdMat, &
                           sup%hydyn,env,time,istep,iter,err)
    if (err < 0) ierr = ierr + err

  end subroutine calcMorisonForces


  !!============================================================================
  !> @brief Calculates buoyancy forces and associated load correction stiffness.
  !>
  !> @param sup Superelement to calculate buoyancy forces for
  !> @param env Environmental data
  !> @param[in] time Current simulation time
  !> @param[in] iter Iteration counter
  !> @param ierr Error flag
  !>
  !> @details For a (partly) submerged body.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 2 Jul 2008

  subroutine calcBuoyancyForces (sup,env,time,iter,ierr)

    use SupElTypeModule      , only : SupElType, HydroDynType, dp
    use EnvironmentTypeModule, only : EnvironmentType
    use HydroDynamicsModule  , only : getBuoyancyForces
    use IdTypeModule         , only : getId
    use rotationModule       , only : EccExpand
    use profilerModule       , only : startTimer, stopTimer, hyd_p
#ifdef FT_DEBUG
    use dbgUnitsModule       , only : dbgHD
#endif
    use reportErrorModule    , only : reportError, debugFileOnly_p

    type(SupElType)      , intent(inout) :: sup
    type(EnvironmentType), intent(inout) :: env
    real(dp)             , intent(in)    :: time
    integer              , intent(in)    :: iter
    integer              , intent(inout) :: ierr

    !! Local variables
    integer  :: err
    real(dp) :: g, eVec(3), kL(6,6)

    !! --- Logic section ---

    if (sup%rigidFlag > 0) then
       call getBuoyancyForces ('Superelement'//getId(sup%id),sup%triads, &
            &                  sup%supTr,sup%hydyn,env,g,time,iter,ierr=err)
    else
       call getBuoyancyForces ('Beam element'//getId(sup%id),sup%triads, &
            &                  sup%supTr,sup%hydyn,env,g,time,iter,sup%Q,err)
    end if
    if (err < 0) then
       ierr = ierr + err
       call reportError (debugFileOnly_p,'calcBuoyancyForces')
       return
    else if (sup%rigidFlag <= 0) then
       return
    end if

    call startTimer (hyd_p)

    !! Buoyancy force vector in the center of gravity triad
    eVec = sup%hydyn%C0b(:,1) - sup%Trundeformed(:,4,1)
    call EccExpand (eVec,env%rhow*g*sup%hydyn%Vb(1)*sup%hydyn%wn,sup%hydyn%B)
#ifdef FT_DEBUG
    write(dbgHD,"(7X,'B  =',1P6E13.5)") sup%hydyn%B
    write(dbgHD,"(11X,1P6E13.5)") matmul(sup%supTr(:,1:3),sup%hydyn%B(1:3))
#endif

    !! Add to the superelement external load vector
    if (associated(sup%Bgp)) then
       sup%Q = sup%Q + matmul(sup%hydyn%B,sup%Bgp)
    else
       sup%Q(1:6) = sup%Q(1:6) + sup%hydyn%B
    end if

    if (sup%hydyn%As(1) > 0.0_dp .and. associated(sup%KlMat)) then
       !! Load-correction stiffness
       call EccExpand (eVec,-env%rhow*g*sup%hydyn%As(1)*sup%hydyn%wn,kL)
#ifdef FT_DEBUG
       write(dbgHD,"(7X,'kL =',1P6E13.5/(11X,1P6E13.5))") (kL(err,:),err=1,6)
#endif
       if (associated(sup%Bgp)) then
          sup%KlMat = matmul(transpose(sup%Bgp),matmul(kL,sup%Bgp))
       else
          sup%KlMat = 0.0_dp
          sup%KlMat(1:6,1:6) = kL
       end if
    else if (associated(sup%KlMat)) then
       sup%KlMat = 0.0_dp
    end if

    call stopTimer (hyd_p)

  end subroutine calcBuoyancyForces


  !!============================================================================
  !> @brief Adds superelement forces into corresponding system force vectors.
  !>
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[in] sups All superelements in the model
  !> @param FSk Internal stiffness force vector
  !> @param FDk Internal damping force vector
  !> @param FIk Internal internal force vector
  !> @param Qk  External force vector (gravitation loads)
  !> @param RFk Reaction forces
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date December 2001

  subroutine AddInSupForces (sam,sups,FSk,FDk,FIk,Qk,RFk,ierr)

    use SamModule         , only : SamType
    use SupElTypeModule   , only : SupElType, dp
    use IdTypeModule      , only : getId
    use AsmExtensionModule, only : csAddEV
    use profilerModule    , only : startTimer, stopTimer, sup2_p
    use reportErrorModule , only : reportError, error_p

    type(SamType)  , intent(in)    :: sam
    type(SupElType), intent(in)    :: sups(:)
    real(dp)       , intent(inout) :: FSk(:), FDk(:), FIk(:), Qk(:), RFk(:)
    integer        , intent(out)   :: ierr

    !! Local variables
    integer :: i, iel, err, nedof

    !! --- Logic section ---

    call startTimer (sup2_p)

    ierr = 0
    do i = 1, size(sups)
       iel   = sups(i)%samElNum
       nedof = sups(i)%nTotDofs
       call csAddEV (sam, iel, nedof, sups(i)%FS, FSk, RFk, err); ierr=ierr+err
       call csAddEV (sam, iel, nedof, sups(i)%FD, FDk, RFk, err); ierr=ierr+err
       call csAddEV (sam, iel, nedof, sups(i)%FI, FIk, RFk, err); ierr=ierr+err
       call csAddEV (sam,-iel, nedof, sups(i)%Q , Qk,  RFk, err); ierr=ierr+err
       if (ierr /= 0) then
          call stopTimer (sup2_p)
          call reportError (error_p,'Inconsistent control information.', &
               &            'Encountered during assembly of superelement'// &
               &            getId(sups(i)%id),addString='AddInSupForces')
          return
       end if
    end do

    call stopTimer (sup2_p)

  end subroutine AddInSupForces


  !!============================================================================
  !> @brief Adds superelement forces into corresponding system force vectors.
  !>
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[in] sups All superelements in the model
  !> @param FSk Internal stiffness force vector
  !> @param Qk External force vector (gravitation loads)
  !> @param RFk Reaction forces
  !> @param[out] ierr Error flag
  !>
  !> @details This subroutine only considers the static forces, the stiffness
  !> and external forces.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 5 June 2002

  subroutine AddInStaticSupForces (sam,sups,FSk,Qk,RFk,ierr)

    use SamModule         , only : SamType
    use SupElTypeModule   , only : SupElType, dp
    use IdTypeModule      , only : getId
    use AsmExtensionModule, only : csAddEV
    use profilerModule    , only : startTimer, stopTimer, sup2_p
    use reportErrorModule , only : reportError, error_p

    type(SamType)  , intent(in)    :: sam
    type(SupElType), intent(in)    :: sups(:)
    real(dp)       , intent(inout) :: FSk(:)
    real(dp)       , intent(inout) :: Qk(:)
    real(dp)       , intent(inout) :: RFk(:)
    integer        , intent(out)   :: ierr

    !! Local variables
    integer :: i, iel, err, nedof

    !! --- Logic section ---

    call startTimer (sup2_p)

    ierr = 0
    do i = 1, size(sups)
       iel   = sups(i)%samElNum
       nedof = sups(i)%nTotDofs
       call csAddEV (sam, iel, nedof, sups(i)%FS, FSk, RFk, err); ierr=ierr+err
       call csAddEV (sam,-iel, nedof, sups(i)%Q , Qk,  RFk, err); ierr=ierr+err
       if (ierr /= 0) then
          call stopTimer (sup2_p)
          call reportError (error_p,'Inconsistent control information.', &
               &            'Encountered during assembly of superelement'// &
               &            getId(sups(i)%id),addString='AddInStaticSupForces')
          return
       end if
    end do

    call stopTimer (sup2_p)

  end subroutine AddInStaticSupForces


  !!============================================================================
  !> @brief Adds a superelement matrix into the equivalent system matrix.
  !>
  !> @param[in] supMat The superelement matrix to add into the system
  !> @param sysMat System coefficient matrix (stiffness, mass, ...)
  !> @param[in] sup The superelement the provided matrix is associated with
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[out] err Error flag
  !> @param sysRhs System right-hand-side vector
  !> @param[in] scale Optional scaling factor for the matrix to add
  !>
  !> @details This subroutine transforms the superelement matrix @a supMat
  !> to the system directions for each DOF associated with it, and then adds it
  !> into the system matrix @a sysMat. If any of the superelement DOFs are
  !> prescribed, the associated force contributions are added into the system
  !> force vector @a sysRhs, if present.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Karl Erik Thoresen
  !>
  !> @date Januar 1999

  subroutine AddInSupMat (supMat,sysMat,sup,sam,err,sysRhs,scale)

    use SamModule          , only : SamType
    use SysMatrixTypeModule, only : SysMatrixType
    use SupElTypeModule    , only : SupElType, GetRsupToSys, dp
    use AsmExtensionModule , only : csAddEM
    use scratchArrayModule , only : getRealScratchMatrix
    use profilerModule     , only : startTimer, stopTimer, sup1_p
    use reportErrorModule  , only : getErrorFile, reportError, debugFileOnly_p

    real(dp)           , intent(in)    :: supMat(:,:)
    type(SysMatrixType), intent(inout) :: sysMat
    type(SupElType)    , intent(in)    :: sup
    type(SamType)      , intent(in)    :: sam
    integer            , intent(out)   :: err
    real(dp), optional , intent(inout) :: sysRhs(:)
    real(dp), optional , intent(in)    :: scale

    !! Local variables
    integer           :: i, n, nDim, lpu
    real(dp)          :: RsupToSys(3,3), tmp(3)
    real(dp), pointer :: sMat(:,:)

    !! --- Logic section ---

    call startTimer (sup1_p)

    lpu = getErrorFile()
    nDim = size(supMat,1)
    sMat => getRealScratchMatrix(nDim,nDim,err)
    if (err < 0) goto 900

    call DCOPY (nDim*nDim,supMat(1,1),1,sMat(1,1),1)
    if (present(scale)) then
       call DSCAL (nDim*nDim,scale,sMat(1,1),1)
    end if

    do i = 1, sup%nExtNods
       n = sup%triads(i)%firstDOF

       if (sup%triads(i)%p%nDOFs >= 3) then
          RsupToSys = transpose(GetRsupToSys(sup,i))
          call MATTRA (RsupToSys(1,1),sMat(1,1),tmp(1),nDim,3,n,lpu,err)
          if (err < 0) goto 900
       end if

       if (sup%triads(i)%p%nDOFs >= 6) then
          call MATTRA (RsupToSys(1,1),sMat(1,1),tmp(1),nDim,3,n+3,lpu,err)
          if (err < 0) goto 900
       end if

    end do

    !! Add this matrix to the system matrix
    call csAddEM (sam,sup%samElNum,sup%nTotDofs,sMat,sysMat,err,sysRhs)
    if (err == 0) goto 999

900 call reportError (debugFileOnly_p,'AddInSupMat')
999 call stopTimer (sup1_p)

  end subroutine AddInSupMat


  !!============================================================================
  !> @brief Computes the tangential superelement stiffness matrix.
  !>
  !> @param sup The superelement to calculate tangent stiffness for
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Dag Rune Christensen                              @date 23 Oct 1997
  !> @author Bjorn Haugen                                      @date  3 Oct 2003
  !> @author Knut Morten Okstad                                @date 10 Apr 2019

  subroutine CompTanStiff (sup,ierr)

    use SupElTypeModule  , only : SupElType, dp
    use corotUtilModule  , only : addKgrToKm
    use reportErrorModule, only : reportError, debugFileOnly_p

    type(SupElType), intent(inout) :: sup
    integer        , intent(out)   :: ierr

    !! Local variables
    integer :: m, n

    !! --- Logic section ---

    m = size(sup%KmMat,1)
    n = size(sup%KmMat,2)

    !! Calculate internal forces, Fs = Km*Finit
    call scaledMatmul (m,n,1.0_dp,sup%KmMat,sup%Finit,sup%FS)

    !! Add geometric stiffness to material stiffness
    call addKgrToKm (sup,sup%KtMat,sup%KmMat,sup%FS,ierr)
    if (ierr < 0) then
       call reportError (debugFileOnly_p,'CompTanStiff')

    else if (associated(sup%KlMat)) then
       !! Add load correction stiffness
       call DAXPY (m*n,1.0_dp,sup%KlMat(1,1),1,sup%KtMat(1,1),1)
    end if

  end subroutine CompTanStiff


  !!============================================================================
  !> @brief Computes the superelement Newton matrix.
  !>
  !> @param[in] newTangent If .true., the tangential stiffness is updated
  !> @param[in] scaleM Mass matrix scaling factor
  !> @param[in] scaleC Damping matrix scaling factor
  !> @param[in] scaleK Stiffness matrix scaling factor
  !> @param sup The superelement to calculate tangent stiffness for
  !> @param[out] ierr Error flag
  !>
  !> @details The Newton matrix is a linear combination of the stiffness-,
  !> damping- (if any) and stiffness matrices of the superelement.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 10 Apr 2019

  subroutine BuildSupNewtonMat (newTangent,scaleM,scaleC,scaleK,sup,ierr)

    use SupElTypeModule  , only : SupElType, dp
    use reportErrorModule, only : reportError, debugFileOnly_p

    logical        , intent(in)    :: newTangent
    real(dp)       , intent(in)    :: scaleM, scaleC, scaleK
    type(SupElType), intent(inout) :: sup
    integer        , intent(out)   :: ierr

    !! Local variables
    integer :: n

    !! --- Logic section ---

    n = size(sup%Nmat)
    if (associated(sup%Mmat)) then

       !! Structural mass matrix contribution
       call DCOPY (n,sup%Mmat(1,1),1,sup%Nmat(1,1),1)
       call DSCAL (n,scaleM*sup%massScl(1),sup%Nmat(1,1),1)

    else
       call DCOPY (n,0.0_dp,0,sup%Nmat(1,1),1)
    end if

    if (associated(sup%Cmat)) then

       !! Structural damping matrix contribution.
       !! We here assume that only stiffness-proportional damping is used,
       !! unless the stiffness-scaling factor is equal to one.
       call DAXPY (n,scaleC*sup%stifScl(1),sup%Cmat(1,1),1,sup%Nmat(1,1),1)

    end if
    if (sup%stressStiffFlag(1) == 1) then

       if (newTangent) then
          call CompTanStiff (sup,ierr)
          if (ierr < 0) then
             call reportError (debugFileOnly_p,'BuildSupNewtonMat')
             return
          end if
       end if

       !! Tangential stiffness matrix contribution
       call DAXPY (n,scaleK*sup%stifScl(1),sup%Ktmat(1,1),1,sup%Nmat(1,1),1)

    else

       !! Material stiffness matrix contribution
       call DAXPY (n,scaleK*sup%stifScl(1),sup%Kmmat(1,1),1,sup%Nmat(1,1),1)

    end if
    if (associated(sup%Mamat)) then

       !! Mass matrix contributions due to hydrodynamic added mass
       call DAXPY (n,scaleM,sup%Mamat(1,1),1,sup%Nmat(1,1),1)

    end if
    if (associated(sup%Cdmat)) then

       !! Damping matrix contributions due to hydrodynamic drag
       call DAXPY (n,scaleC,sup%Cdmat(1,1),1,sup%Nmat(1,1),1)

    end if

  end subroutine BuildSupNewtonMat


  !!============================================================================
  !> @brief Calculates the scaled matrix-vector product Y = &alpha;*A*X.
  !>
  !> @param[in] m Number of rows in matrix @b A
  !> @param[in] n Number of columns in matrix @b A and length of vector @b X.
  !> @param[in] alpha Scaling factor, &alpha;
  !> @param[in] A The matrix to multiply with
  !> @param[in] X The vector to multiply with
  !> @param[out] Y The output vector
  !> @param[in] ldA Leading dimension of matrix @b A
  !>
  !> @details This is just a convenience wrapper over the BLAS subroutine DGEMV.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 21 Oct 2019

  subroutine scaledMatmul (m,n,alpha,A,X,Y,ldA)

    use kindModule, only : dp

    integer , intent(in)  :: m, n
    real(dp), intent(in)  :: alpha, A(:,:), X(:)
    real(dp), intent(out) :: Y(:)
    integer , intent(in), optional :: ldA

    !! --- Logic section ---

    if (present(ldA)) then
       call DGEMV ('N',m,n,alpha,A(1,1),ldA,X(1),1,0.0_dp,Y(1),1)
    else
       call DGEMV ('N',m,n,alpha,A(1,1),m,X(1),1,0.0_dp,Y(1),1)
    end if

  end subroutine scaledMatmul

end module SupElRoutinesModule
