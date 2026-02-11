!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file forceRoutinesModule.f90
!>
!> @brief Subroutines for external point load calculations.

!!==============================================================================
!> @brief Module with subroutines for external point load calculations.
!>
!> @details This module contains a set of subroutines for performing various
!> computation tasks on the forcetypemodule::forcetype objects in the model
!> (external point loads) during the dynamic or quasi-static simulation.

module ForceRoutinesModule

  implicit none

  private

  public :: updateExternalForces, addInExternalForces


contains

  !!============================================================================
  !> @brief Updates the external forces (point loads).
  !>
  !> @param forces All external point loads in the model
  !> @param[in] iFlag Option telling how to update the forces (see below)
  !> @param ierr Error flag
  !>
  !> @details The way the forces are updated are governed by the value of the
  !> argument @a iFlag as follows:
  !> - iFlag < 0: Update all forces, but not the nodal force vectors (restart)
  !> - iFlag = 0: Update all forces
  !> - iFlag = 1: Update forces that depends on the previous configuration only
  !> - iFlag = 2: Update forces that depends on current (unknown) configuration
  !>
  !> @callgraph @callergraph
  !>
  !> @author Karl Erik Thoresen
  !> @date Feb 1999
  !>
  !> @author Knut Morten Okstad
  !> @date Feb 2018

  subroutine updateExternalForces (forces,iFlag,ierr)

    use ForceTypeModule     , only : ForceType, forceVec_p, momentVec_p, dp
    use IdTypeModule        , only : getId
    use TriadTypeModule     , only : updateNodeForce
    use EngineRoutinesModule, only : Evaluate
    use manipMatrixModule   , only : matmul34
    use reportErrorModule   , only : reportError
    use reportErrorModule   , only : error_p, warning_p, debugFileOnly_p

    type(ForceType), intent(inout) :: forces(:)
    integer        , intent(in)    :: iFlag
    integer        , intent(inout) :: ierr

    !! Local variables
    integer  :: i, err, lerr
    real(dp) :: fDir(3), eV(6), dLen

    !! --- Logic section ---

    lerr = ierr
    do i = 1, size(forces)
       if (forces(i)%loadType > 0) cycle ! Skip frequency-domain loads

       !! Get the engine value for this force in current state
       if (iFlag <= 0 .or. iFlag == mod(forces(i)%updateFlag,10)) then
          forces(i)%F = Evaluate(forces(i)%engine, &
               &                 forces(i)%f1,forces(i)%f0,ierr)
       end if
       if (ierr < lerr) then
          lerr = ierr
          call reportError (error_p,'Failed to evaluate force '// &
               &                    'engine for Force'//getId(forces(i)%id))
          cycle
       end if

       if (forces(i)%dof == forceVec_p .or. forces(i)%dof == momentVec_p) then

          !! Set up a vector describing the force direction
          if (associated(forces(i)%sup1)) then
             fDir = matmul34(forces(i)%sup1%supTr,forces(i)%v1)
          else if (associated(forces(i)%trp1)) then
             fDir = matmul34(forces(i)%trp1%ur,forces(i)%v1)
          else
             fDir = forces(i)%v1
          end if
          if (associated(forces(i)%sup2)) then
             fDir = matmul34(forces(i)%sup2%supTr,forces(i)%v2) - fDir
          else if (associated(forces(i)%trp2)) then
             fDir = matmul34(forces(i)%trp2%ur,forces(i)%v2)
          else
             fDir = forces(i)%v2 - fDir
          end if

          !! Normalize the direction vector
          dLen = sqrt(dot_product(fDir,fDir))
          if (dLen <= tiny(fDir(1))*10.0_dp) then
             call reportError (warning_p,'Force'//getId(forces(i)%id), &
                  'Could not find direction for this force because '// &
                  'the points defining the direction coincide.', &
                  'Applying force = 0.0 and continuing.')
             forces(i)%F        = 0.0_dp
             forces(i)%forceDir = 0.0_dp
          else
             forces(i)%forceDir = fDir / dLen
          end if

       end if
       if ((iFlag == 0 .or. iFlag == 2) .and. associated(forces(i)%triad)) then

          !! Establish a force vector for this triad
          call calcExtTriadForce (forces(i),eV,err)
          if (err < 0) ierr = ierr - 1

          !! Update the nodal forces of this triad
          call updateNodeForce (forces(i)%triad,-eV)

       end if

    end do

    if (ierr < lerr) call reportError (debugFileOnly_p,'updateExternalForces')

  end subroutine updateExternalForces


  !!============================================================================
  !> @brief Assembles system force vector contributions from the external loads.
  !>
  !> @param Qk System external force vector
  !> @param RFk System reaction forces associated with constrained DOFs
  !> @param[in] forces All external point loads in the model
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Karl Erik Thoresen
  !>
  !> @date Feb 1999

  subroutine addInExternalForces (Qk,RFk,forces,sam,ierr)

    use ForceTypeModule   , only : ForceType, dp
    use SamModule         , only : SamType
    use AsmExtensionModule, only : csAddNV
    use reportErrorModule , only : internalError, reportError, debugFileOnly_p

    type(ForceType), intent(in)    :: forces(:)
    type(SamType)  , intent(in)    :: sam
    real(dp)       , intent(inout) :: Qk(:), RFk(:)
    integer        , intent(inout) :: ierr

    !! Local variables
    integer  :: i, err, lerr, nDOFs, samNodNum
    real(dp) :: eV(10)

    !! --- Logic section ---

    lerr = ierr
    do i = 1, size(forces)
       if (forces(i)%loadType > 0) cycle ! Skip frequency-domain loads

       if (associated(forces(i)%triad)) then

          nDOFs = forces(i)%triad%nDOFs
          if (nDOFs < 1) then
             cycle
          else if (nDOFs > size(eV)) then
             ierr = ierr + internalError('addInExternalForces 10')
             cycle
          end if

          !! Establish a force vector for this triad
          call calcExtTriadForce (forces(i),eV,err)
          if (err < 0) ierr = ierr - 1

          samNodNum = forces(i)%triad%samNodNum

       else if (associated(forces(i)%joint)) then

          nDOFs = forces(i)%joint%nJointDOFs
          if (nDOFs < 1) then
             cycle
          else if (nDOFs > size(eV)) then
             ierr = ierr + internalError('addInExternalForces 20')
             cycle
          end if

          !! Establish a force vector for this joint
          eV = 0.0_dp
          eV(forces(i)%dof) = forces(i)%F

          samNodNum = forces(i)%joint%samNodNum

       else
          cycle
       end if

       call csAddNV (sam, -samNodNum, eV(1:nDOFs), Qk, RFk, err)
       if (err /= 0) ierr = ierr - 1

    end do

    if (ierr < lerr) call reportError (debugFileOnly_p,'addInExternalForces')

  end subroutine addInExternalForces


  !!============================================================================
  !> @brief Calculates system force vector contributions from an external load.
  !>
  !> @param[in] force The external point load object to calculate forces for
  !> @param[out] eV Nodal force vector
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 30 Jan 2008

  subroutine calcExtTriadForce (force,eV,ierr)

    use ForceTypeModule  , only : ForceType, forceVec_p, momentVec_p, dp
    use TriadTypeModule  , only : transGlobToSys, transTriadToSys
    use reportErrorModule, only : internalError

    type(ForceType), intent(in)  :: force
    real(dp)       , intent(out) :: eV(6)
    integer        , intent(out) :: ierr

    !! --- Logic section ---

    ierr = 0
    eV = 0.0_dp

    select case (force%dof)

    case (forceVec_p)
       eV(1:3) = force%F * force%forceDir
       call transGlobToSys (force%triad,eV(1:3))

    case (momentVec_p)
       eV(4:6) = force%F * force%forceDir
       call transGlobToSys (force%triad,eV(4:6))

    case (1:6) ! The force is acting directly on a nodal dof
       eV(force%dof) = force%F
       if (force%updateFlag > 9) then
          !! Apply the load in local triad axes.
          !! For backward compatibility with earlier versions (pre R7.3).
          call transTriadToSys (force%triad,eV(1:force%triad%nDOFs))
       end if

    case default
       ierr = internalError('calcExtTriadForce')

    end select

  end subroutine calcExtTriadForce

end module ForceRoutinesModule
