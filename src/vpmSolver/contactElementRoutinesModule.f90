!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file contactElementRoutinesModule.f90
!> @brief Subroutines for contact element calculations.

!!==============================================================================
!> @brief Module with subroutines for contact element calculations.

module ContactElementRoutinesModule

  use KindModule              , only : dp, epsDiv0_p
  use ContactElementTypeModule, only : ContactElementType

  implicit none

  private

  real(dp), parameter :: epsZero_p = 1.0e-15_dp !< Zero tolerance

  public :: addInContactElementStiffMat, addInContactElementDamperMat
  public :: addInContactElementForces, updateContactElements


contains

#ifdef FT_DEBUG
  !> @cond NO_DOCUMENTATION
  subroutine writeMat (mType,eM,elmId)
    use IdTypeModule     , only : IdType, getId
    use manipMatrixModule, only : writeObject
    use dbgUnitsModule   , only : dbgSolve
    character(len=*), intent(in) :: mType
    real(dp)        , intent(in) :: eM(:,:)
    type(IdType)    , intent(in) :: elmId
    if (dbgSolve > 0) then
       if (size(eM,2) > 1) then
          call writeObject (eM,dbgSolve,trim(mType)// &
               &            ' matrix for Contact element'//trim(getId(elmId)))
       else
          call writeObject (reshape(eM,(/size(eM,1)/)),dbgSolve,trim(mType)// &
               &            ' vector for Contact element'//trim(getId(elmId)))
       end if
    end if
  end subroutine writeMat
  !> @endcond
#endif


  !!============================================================================
  !> @cond FULL_DOC
  !> @brief Utility for looping over the current active curve point indices.
  !> @callergraph
  !> @author Knut Morten Okstad
  !> @date 16 Jun 2006
  !> @endcond

  function getActivePointIndex (cElem,cpIndex)

    type(ContactElementType), intent(in) :: cElem
    integer                 , intent(in) :: cpIndex
    integer                              :: getActivePointIndex

    !! --- Logic section ---

    getActivePointIndex = 0
    if (cpIndex == cElem%lastDomain) then
       getActivePointIndex = 1
    else if (.not. cElem%isAtCorner) then
       if (cpIndex == 1 + mod(cElem%lastDomain,size(cElem%cSurf%CPoints))) then
          getActivePointIndex = 2
       end if
    end if

  end function getActivePointIndex


  !!============================================================================
  !> @brief Calculates the stiffness matrix for each contact spring.
  !>
  !> @param[in] includeStressStiff If .true., geometric stiffness is included
  !> @param[in] scaleK Stiffness scaling factor
  !> @param Nmat System matrix to add stiffness contributions into
  !> @param[in] cElem The contact element to calculate stiffness matrix for
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[out] ierr Error flag
  !> @param Rhs System right-hand-side vector associated with @a Nmat
  !>
  !> @details The element stiffness matrix is added into
  !> the system Newton matrix, multiplied by the factor @a scaleK.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !> @date Oct 2001
  !>
  !> @author Knut Morten Okstad
  !> @date Sep 2005

  subroutine addInContactElementStiffMat (includeStressStiff, scaleK, &
       &                                  Nmat, cElem, sam, ierr, Rhs)

    use SamModule               , only : SamType
    use SysMatrixTypeModule     , only : SysMatrixType
    use ContactElementTypeModule, only : ignoreEccN
    use TriadTypeModule         , only : transGlobToSys
    use AsmExtensionModule      , only : csAddEM
    use scratchArrayModule      , only : getRealScratchMatrix
    use manipMatrixModule       , only : dyadic_product
    use rotationModule          , only : EccExpand
    use reportErrorModule       , only : reportError, debugFileOnly_p

    logical                 , intent(in)    :: includeStressStiff
    real(dp)                , intent(in)    :: scaleK
    type(SysMatrixType)     , intent(inout) :: Nmat
    type(ContactElementType), intent(in)    :: cElem
    type(SamType)           , intent(in)    :: sam
    integer                 , intent(out)   :: ierr
    real(dp), optional      , intent(inout) :: Rhs(:)

    !! Local variables
    integer           :: i, j, ip, nDOFs, nDOFm
    logical           :: hasStiff(3)
    real(dp)          :: stf, stiffInG(3,3), rotStiff(3,3), sFull(6,6)
    real(dp)          :: KgNom, KgMat(3,3)
    real(dp), pointer :: eM(:,:)

    !! --- Logic section ---

    ierr = 0
    if (.not. cElem%isActive) return

    stiffInG = 0.0_dp
    rotStiff = 0.0_dp
    hasStiff = .false.

    !! --- Contact point stiffness from contact springs

    do i = 1, 3
       if (associated(cElem%springs(i)%p)) then

          !! Material stiffness contribution from spring
          stf = cElem%springs(i)%p%stiffness * scaleK
          if (abs(stf) > epsZero_p) hasStiff(1) = .true.
          stiffInG = stiffInG + dyadic_product(cElem%CPosInG(:,i)) * stf

       end if
       if (associated(cElem%springs(3+i)%p)) then

          !! Material stiffness contribution from rotational spring
          stf = cElem%springs(3+i)%p%stiffness * scaleK
          if (abs(stf) > epsZero_p) hasStiff(2) = .true.
          rotStiff = rotStiff + dyadic_product(cElem%CPosInG(:,i)) * stf

       end if
    end do

    if (includeStressStiff .and. associated(cElem%springs(1)%p)) then

       !! Geometric stiffness contribution from normal spring
       KgNom = 0.0_dp
       if (abs(cElem%springs(1)%p%length) > epsDiv0_p) then
          KgNom = cElem%springs(1)%p%force / cElem%springs(1)%p%length
          if (abs(KgNom) > cElem%springs(1)%p%stiffness) then
             !! Use the material stiffness as an upper limit to avoid
             !! numerical stability issues when the length is close to zero
             KgNom = sign(cElem%springs(1)%p%stiffness,KgNom) * scaleK
          else
             KgNom = KgNom * scaleK
          end if
       else if (abs(cElem%springs(1)%p%force) <= epsZero_p) then
          !! Both the length and force is zero, use l'Hopital's rule assuming
          !! that force = stiff*length (i.e. length0 must be zero too)
          if (abs(cElem%springs(1)%p%length0) <= epsZero_p) then
             KgNom = cElem%springs(1)%p%stiffness * scaleK
          end if
       end if

       if (abs(KgNom) > epsZero_p) then
          if (cElem%radius > 0.0_dp) then
             !! Radial contact, add geometric stiffness in the local Y-direction
             hasStiff(3) = .true.
             KgMat = KgNom * dyadic_product(cElem%CPosInG(:,2))
          else if (cElem%isAtCorner) then
             KgMat = 0.0_dp
          end if
          if (cElem%isAtCorner) then
             !! Corner point, add geometric stiffness in the local Z-direction
             hasStiff(3) = .true.
             KgMat = KgMat + KgNom * dyadic_product(cElem%CPosInG(:,3))
          end if
       end if

    end if

    if (associated(cElem%friction)) then
       if (associated(cElem%friction%spr)) then

          !! Material stiffness contribution from friction spring
          stf = cElem%friction%spr%stiffness * scaleK
          if (abs(stf) > epsZero_p) hasStiff(1) = .true.
          stiffInG = stiffInG + dyadic_product(cElem%CPosInG(:,3)) * stf

       end if
    end if

    if (.not. any(hasStiff)) return

    !! Allocate scratch array for the element matrix, eM
    eM => getRealScratchMatrix(cElem%nDOFs,cElem%nDOFs,ierr)
    if (ierr /= 0) goto 900

    !! --- Build element stiffness matrix in global system

    eM    = 0.0_dp
    sFull = 0.0_dp

    !! Take into account possible eccentricity
    !!    | Fm | = |    I    :  0 || Kt : 0  || I : Spin(e)' |
    !!    | Mm |   | Spin(e) :  I ||  0 : Kr || 0 :   I      |

    !! Triad1 (follower) stiffness
    nDOFs = cElem%triad1%nDOFs
    if (ignoreEccN) then
       sFull(1:3,1:3) = stiffInG
    else if (nDOFs > 0) then
       call EccExpand (cElem%eVec(:,2),stiffInG,sFull)
    end if
    if (hasStiff(3)) sFull(1:3,1:3) = sFull(1:3,1:3) + KgMat
    if (hasStiff(2)) sFull(4:6,4:6) = sFull(4:6,4:6) + rotStiff

    if (nDOFs > 0) eM(1:nDOFs,1:nDOFs) = sFull(1:nDOFs,1:nDOFs)
    ip = 1 + nDOFs

    !! Contributions to curve point DOFs
    if (.not. ignoreEccN) then
       call EccExpand (cElem%eVec(:,1),stiffInG,sFull)
       if (hasStiff(3)) sFull(1:3,1:3) = sFull(1:3,1:3) + KgMat
       if (hasStiff(2)) sFull(4:6,4:6) = sFull(4:6,4:6) + rotStiff
    end if

    do j = 1, size(cElem%CSurf%CPoints)
       i = getActivePointIndex(cElem,j)
       nDOFm = cElem%CSurf%CPoints(j)%triad%nDOFs
       if (i > 0 .and. nDOFm > 0) then
          stf = cElem%coeff(i)
          eM(ip:ip+nDOFm-1, ip:ip+nDOFm-1) = sFull(1:nDOFm,1:nDOFm) * stf*stf
          if (nDOFs > 0) then
             eM( 1:   nDOFs  , ip:ip+nDOFm-1) = -sFull(1:nDOFs,1:nDOFm) * stf
             eM(ip:ip+nDOFm-1,  1:   nDOFs  ) = -sFull(1:nDOFm,1:nDOFs) * stf
          end if
       end if
       ip = ip + nDOFm
    end do

    !! --- Transform to system directions

    ip = 1
    if (nDOFs > 0) then
       call transGlobToSys (cElem%triad1, eM, ip, nDOFs)
       ip = ip + nDOFs
    end if

    do j = 1, size(cElem%CSurf%CPoints)
       i = getActivePointIndex(cElem,j)
       nDOFm = cElem%CSurf%CPoints(j)%triad%nDOFs
       if (i > 0 .and. nDOFm > 0) then
          call transGlobToSys (cElem%CSurf%CPoints(j)%triad, eM, ip, nDOFm)
       end if
       ip = ip + nDOFm
    end do

    !! --- Add to system matrix

#ifdef FT_DEBUG
    call writeMat ('Stiffness',eM,cElem%id)
#endif
    call csAddEM (sam,cElem%samElnum,cElem%nDOFs,eM,Nmat,ierr,Rhs)
    if (ierr == 0) return

900 call reportError (debugFileOnly_p,'addInContactElementStiffMat')

  end subroutine addInContactElementStiffMat


  !!============================================================================
  !> @brief Calculates the damping matrix for each contact damper.
  !>
  !> @param[in] includeStressStiff If .true., geometric stiffness is included
  !> @param[in] scaleC Damping scaling factor
  !> @param[in] scaleK Stiffness scaling factor
  !> @param Nmat System matrix to add damping contributions into
  !> @param[in] cElem The contact element to calculate damping matrix for
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[out] ierr Error flag
  !> @param Rhs System right-hand-side vector associated with @a Nmat
  !>
  !> @details The element damping matrix are added into
  !> the system Newton matrix, multiplied by the factor @a scaleC.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !> @date Oct 2001
  !>
  !> @author Knut Morten Okstad
  !> @date Sep 2005

  subroutine addInContactElementDamperMat (includeStressStiff, scaleC, scaleK, &
       &                                   Nmat, cElem, sam, ierr, Rhs)

    use SamModule               , only : SamType
    use SysMatrixTypeModule     , only : SysMatrixType
    use ContactElementTypeModule, only : ignoreEccN
    use TriadTypeModule         , only : transGlobToSys
    use AsmExtensionModule      , only : csAddEM
    use scratchArrayModule      , only : getRealScratchMatrix
    use manipMatrixModule       , only : dyadic_product
    use rotationModule          , only : EccExpand
    use reportErrorModule       , only : reportError, debugFileOnly_p

    logical                 , intent(in)    :: includeStressStiff
    real(dp)                , intent(in)    :: scaleC, scaleK
    type(SysMatrixType)     , intent(inout) :: Nmat
    type(ContactElementType), intent(in)    :: cElem
    type(SamType)           , intent(in)    :: sam
    integer                 , intent(out)   :: ierr
    real(dp), optional      , intent(inout) :: Rhs(:)

    !! Local variables
    integer           :: i, j, ip, nDOFs, nDOFm
    logical           :: hasDamp(3)
    real(dp)          :: dmp, dampInG(3,3), rotDamp(3,3), dFull(6,6)
    real(dp)          :: KgNom, KgMat(3,3)
    real(dp), pointer :: eM(:,:)

    !! --- Logic section ---

    ierr = 0
    if (.not. cElem%isActive) return

    dampInG = 0.0_dp
    rotDamp = 0.0_dp
    hasDamp = .false.

    !! --- Contact point damping from contact dampers

    do i = 1, 3
       if (associated(cElem%dampers(i)%p)) then

          !! Damping contribution from damper
          dmp = cElem%dampers(i)%p%coeff * scaleC
          if (abs(dmp) > epsZero_p) hasDamp(1) = .true.
          dampInG = dampInG + dyadic_product(cElem%DPosInG(:,i)) * dmp

       end if
       if (associated(cElem%dampers(3+i)%p)) then

          !! Damping contribution from rotational damper
          dmp = cElem%dampers(3+i)%p%coeff * scaleC
          if (abs(dmp) > epsZero_p) hasDamp(2) = .true.
          rotDamp = rotDamp + dyadic_product(cElem%DPosInG(:,i)) * dmp

       end if
    end do

    if (includeStressStiff .and. associated(cElem%dampers(1)%p)) then

       !! Geometric stiffness contribution from normal damper
       KgNom = 0.0_dp
       if (abs(cElem%dampers(1)%p%length) > epsDiv0_p) then
          KgNom = cElem%dampers(1)%p%force / cElem%dampers(1)%p%length
          if (abs(KgNom) > cElem%springs(1)%p%stiffness) then
             !! When the follower passes through the contact surface with a
             !! significant velocity the geometric stiffness from the damper
             !! may become non-physically high, therefore we use the actual
             !! material stiffness of the associated contact spring as a limit
             KgNom = sign(cElem%springs(1)%p%stiffness,KgNom) * scaleK
          else
             KgNom = KgNom * scaleK
          end if
       else if (abs(cElem%dampers(1)%p%force) > epsZero_p) then
          !! Zero damper length, use the material spring stiffness instead
          KgNom = cElem%springs(1)%p%stiffness * scaleK
          KgNom = sign(KgNom,cElem%dampers(1)%p%force)
       end if

       if (abs(KgNom) > epsZero_p) then
          if (cElem%radius > 0.0_dp .and. .not.cElem%nonRadialDamping) then
             !! Radial contact, add geometric stiffness in the local Y-direction
             hasDamp(3) = .true.
             KgMat = KgNom * dyadic_product(cElem%DPosInG(:,2))
          else if (cElem%isAtCorner) then
             KgMat = 0.0_dp
          end if
          if (cElem%isAtCorner) then
             !! Corner point, add geometric stiffness in the local Z-direction
             hasDamp(3) = .true.
             KgMat = KgMat + KgNom * dyadic_product(cElem%DPosInG(:,3))
          end if
       end if

    end if

    if (.not. any(hasDamp)) return

    !! Allocate scratch array for the element matrix, eM
    eM => getRealScratchMatrix(cElem%nDOFs,cElem%nDOFs,ierr)
    if (ierr /= 0) goto 900

    !! --- Build element damping matrix in global system

    eM    = 0.0_dp
    dFull = 0.0_dp

    !! Take into account possible eccentricity
    !!    | Fm | = |    I    :  0 || Ct : 0  || I : Spin(e)' |
    !!    | Mm |   | Spin(e) :  I ||  0 : Cr || 0 :   I      |

    !! Triad1 (follower) damping
    nDOFs = cElem%triad1%nDOFs
    if (ignoreEccN) then
       dFull(1:3,1:3) = dampInG
    else if (nDOFs > 0) then
       call EccExpand (cElem%eVec(:,2),dampInG,dFull)
    end if
    if (hasDamp(3)) dFull(1:3,1:3) = dFull(1:3,1:3) + KgMat
    if (hasDamp(2)) dFull(4:6,4:6) = dFull(4:6,4:6) + rotDamp

    if (nDOFs > 0) eM(1:nDOFs,1:nDOFs) = dFull
    ip = 1 + nDOFs

    !! Contributions to curve point DOFs
    if (.not. ignoreEccN) then
       call EccExpand (cElem%eVec(:,1),dampInG,dFull)
       if (hasDamp(3)) dFull(1:3,1:3) = dFull(1:3,1:3) + KgMat
       if (hasDamp(2)) dFull(4:6,4:6) = dFull(4:6,4:6) + rotDamp
    end if

    do j = 1, size(cElem%CSurf%CPoints)
       i = getActivePointIndex(cElem,j)
       nDOFm = cElem%CSurf%CPoints(j)%triad%nDOFs
       if (i > 0 .and. nDOFm > 0) then
          dmp = cElem%coeff(i)
          eM(ip:ip+nDOFm-1, ip:ip+nDOFm-1) = dFull(1:nDOFm,1:nDOFm) * dmp*dmp
          if (nDOFs > 0) then
             eM( 1:   nDOFs  , ip:ip+nDOFm-1) = -dFull(1:nDOFs,1:nDOFm) * dmp
             eM(ip:ip+nDOFm-1,  1:   nDOFs  ) = -dFull(1:nDOFm,1:nDOFs) * dmp
          end if
       end if
       ip = ip + nDOFm
    end do

    !! --- Transform to system directions

    ip = 1
    if (nDOFs > 0) then
       call transGlobToSys (cElem%triad1, eM, ip, nDOFs)
       ip = ip + nDOFs
    end if

    do j = 1, size(cElem%CSurf%CPoints)
       i = getActivePointIndex(cElem,j)
       nDOFm = cElem%CSurf%CPoints(j)%triad%nDOFs
       if (i > 0 .and. nDOFm > 0) then
          call transGlobToSys (cElem%CSurf%CPoints(j)%triad, eM, ip, nDOFm)
       end if
       ip = ip + nDOFm
    end do

    !! --- Add to system matrix

#ifdef FT_DEBUG
    call writeMat ('Damping',eM,cElem%id)
#endif
    call csAddEM (sam,cElem%samElnum,cElem%nDOFs,eM,Nmat,ierr,Rhs)
    if (ierr == 0) return

900 call reportError (debugFileOnly_p,'addInContactElementDamperMat')

  end subroutine addInContactElementDamperMat


  !!============================================================================
  !> @brief Calculates system force vector terms from a contact element.
  !>
  !> @param F System right-hand-side vector
  !> @param RF Reaction for vector
  !> @param[in] cElem The contact element to calculate forces for
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[out] ierr Error flag
  !> @param[in] addSprings If .true., add forces from contact springs
  !> @param[in] addDampers If .true., add forces from contact dampers
  !> @param[in] addFriction If .true., add forces from contact friction
  !>
  !> @details The force terms are added into the system vector @a F
  !> and the reaction force vectos @a RF.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date Oct 2001

  subroutine addInContactElementForces (F, RF, cElem, sam, ierr, &
       &                                addSprings, addDampers, addFriction)

    use SamModule               , only : SamType
    use ContactElementTypeModule, only : ignoreEccF
    use TriadTypeModule         , only : transGlobToSys
    use AsmExtensionModule      , only : csAddEV
    use scratchArrayModule      , only : getRealScratchArray
    use manipMatrixModule       , only : cross_product
    use reportErrorModule       , only : reportError, debugFileOnly_p

    real(dp)                , intent(inout) :: F(:), RF(:)
    type(ContactElementType), intent(in)    :: cElem
    type(SamType)           , intent(in)    :: sam
    integer                 , intent(inout) :: ierr
    logical, optional       , intent(in)    :: addSprings,addDampers,addFriction

    !! Local variables
    integer           :: i, j, ip, nDOFs, nDOFm, err
    logical           :: hasForce(3), fricOnly, fricFromSpring
    real(dp)          :: forceInG(3), momInG(3), momM(3), momS(3), fFull(6)
    real(dp), pointer :: eV(:)
#ifdef FT_DEBUG
    character(len=16) :: forceType
#endif

    !! --- Logic section ---

    if (.not. cElem%isActive) return

    forceInG = 0.0_dp
    momInG   = 0.0_dp
    hasForce = .false.
    fricOnly = .true.

    fricFromSpring = .false.
    if (associated(cElem%friction)) then
       if (associated(cElem%friction%spr)) fricFromSpring = .true.
    end if

    do i = 1, 3
       j = 3 + i

       if (present(addSprings)) then
          if (addSprings .and. associated(cElem%springs(i)%p)) then

             !! Force contribution from spring
             forceInG = forceInG + cElem%springs(i)%p%force * cElem%CPosInG(:,i)
             if (abs(cElem%springs(i)%p%force) > epsZero_p) hasForce(1) = .true.
             fricOnly = .false.

          end if
          if (addSprings .and. associated(cElem%springs(j)%p)) then

             !! Force contribution from rotational spring
             momInG = momInG + cElem%springs(j)%p%force * cElem%CPosInG(:,i)
             if (abs(cElem%springs(j)%p%force) > epsZero_p) hasForce(1) = .true.
             fricOnly = .false.

          end if
       end if

       if (present(addDampers)) then
          if (addDampers .and. associated(cElem%dampers(i)%p)) then

             !! Force contribution from damper
             forceInG = forceInG + cElem%dampers(i)%p%force * cElem%DPosInG(:,i)
             if (abs(cElem%dampers(i)%p%force) > epsZero_p) hasForce(2) = .true.
             fricOnly = .false.

          end if
          if (addDampers .and. associated(cElem%dampers(j)%p)) then

             !! Force contribution from rotational damper
             momInG = momInG + cElem%dampers(j)%p%force * cElem%DPosInG(:,i)
             if (abs(cElem%dampers(j)%p%force) > epsZero_p) hasForce(2) = .true.
             fricOnly = .false.

          end if
       end if

    end do

    if (present(addSprings)) then
       if (addSprings .and. fricFromSpring) then

          !! Force from friction spring
          forceInG = forceInG + cElem%friction%force * cElem%CPosInG(:,3)
          if (abs(cElem%friction%force) > epsZero_p) hasForce(1) = .true.
          fricOnly = .false.

       end if
    end if

    if (present(addFriction) .and. .not.fricFromSpring) then
       if (addFriction .and. associated(cElem%friction)) then

          !! Force contribution from friction
          !! The friction routines calculate the friction force to be positive
          !! for positive velocity, but since it is resisting motion it should
          !! have a negative value here, when added to the external force vector
          if (fricOnly) then
             forceInG = -cElem%friction%force * cElem%CPosInG(:,3)
          else
             forceInG = forceInG + cElem%friction%force * cElem%CPosInG(:,3)
          end if
          if (abs(cElem%friction%force) > epsZero_p) hasForce(3) = .true.

       else
          fricOnly = .false.
       end if
    end if

    if (.not. any(hasForce)) return

    !! Allocate scratch array for the element vector, eV
    eV => getRealScratchArray(cElem%nDOFs,err)
    if (err /= 0) goto 900

    !! --- Build element force vector in global system

    eV = 0.0_dp

    !! Take into account possible eccentricity
    !!    | Fm | = |    I    :  0 || F |
    !!    | Mm |   | Spin(e) :  I || M |

    if (ignoreEccF) then
       momM = momInG
       momS = momInG
    else
       momM = momInG + cross_product(cElem%eVec(:,1),forceInG)
       momS = momInG + cross_product(cElem%eVec(:,2),forceInG)
    end if

    !! Add contribution to triad1 (follower)
    fFull(1:3) = forceInG
    fFull(4:6) = momS
    nDOFs = cElem%triad1%nDOFs
    if (nDOFs > 0) then
       eV(1:nDOFs) = fFull(1:nDOFs)
       call transGlobToSys (cElem%triad1,eV(1:nDOFs))
    end if
    ip = 1 + nDOFs

    !! Add contribution to curve point DOFs
    fFull(4:6) = momM
    do j = 1, size(cElem%CSurf%CPoints)
       i = getActivePointIndex(cElem,j)
       nDOFm = cElem%CSurf%CPoints(j)%triad%nDOFs
       if (i > 0 .and. nDOFm > 0) then
          eV(ip:ip+nDOFm-1) = -fFull(1:nDOFm) * cElem%coeff(i)
          call transGlobToSys (cElem%CSurf%CPoints(j)%triad,eV(ip:ip+nDOFm-1))
       end if
       ip = ip + nDOFm
    end do

    !! --- Add the contact element forces to the system internal force vector

#ifdef FT_DEBUG
    if (count(hasForce) > 1) then
       forceType = 'Force'
    else if (hasForce(1)) then
       forceType = 'Stiffness force'
    else if (hasForce(2)) then
       forceType = 'Damping force'
    else
       forceType = 'Friction force'
    end if
    call writeMat (forceType,reshape(eV,(/cElem%nDOFs,1/)),cElem%id)
#endif
    if (fricOnly) then
       call csAddEV (sam,-cElem%samElnum, cElem%nDOFs, eV, F, RF, err)
    else
       call csAddEV (sam, cElem%samElnum, cElem%nDOFs, eV, F, RF, err)
    end if
    if (err == 0) return

900 ierr = ierr - 1
    call reportError (debugFileOnly_p,'addInContactElementForces')

  end subroutine addInContactElementForces


  !!============================================================================
  !> @brief Updates all contact elements.
  !>
  !> @param[in] cElems All contact elements of the model
  !> @param[in] timeStep Time increment size
  !> @param[out] ierr Error flag
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date May 2002

  subroutine updateContactElements (cElems,timeStep,ierr)

#ifdef FT_DEBUG
    use fileUtilitiesModule   , only : getDBGfile
    use dbgUnitsModule        , only : dbgCurve
#endif
    use IdTypeModule          , only : getId
    use reportErrorModule     , only : reportError, warning_p, debugFileOnly_p
    use reportErrorModule     , only : getErrorFile
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getint

    type(ContactElementType), intent(inout) :: cElems(:)
    real(dp)                , intent(in)    :: timeStep
    integer                 , intent(inout) :: ierr

    !! Local variables
    integer           :: i, lerr, iprint, lpu
    real(dp)          :: domainThickness
    character(len=64) :: msg1, msg2

    !! --- Logic section ---

    lpu = getErrorFile()
    call ffa_cmdlinearg_getint ('debug',iprint)

    lerr = ierr
    do i = 1, size(cElems)
       if (cElems(i)%lastDomain < 1) cycle ! Ignore this contact element

#ifdef FT_DEBUG
       if (size(cElems) < 4 .or. cElems(i)%id%userId == -iprint) then
          dbgCurve = getDBGfile(7,'contactCurve.dbg')
       else
          dbgCurve = 0
       end if
#endif
       call updateContactPosition (cElems(i),cElems(i)%cSurf,lpu)
       call updateContactSprDmp (cElems(i),timeStep,ierr)
       if (.not. cElems(i)%isActive .or. ierr < lerr) cycle

       if (cElems(i)%radius > 0.0_dp) then
          domainThickness = cElems(i)%radius
       else
          domainThickness = cElems(i)%thickness
       end if

       if (cElems(i)%springs(1)%p%length0 >= domainThickness) then
          cElems(i)%nWarning = cElems(i)%nWarning + 1
          if (cElems(i)%nWarning > 5) cycle
          write(msg1,600) cElems(i)%springs(1)%p%length0
          write(msg2,601) domainThickness
600       format('The stress-free length of the Tx spring,',1PE12.5)
601       format('is larger than the contact domain Thickness,',1PE12.5)
          call reportError (warning_p, &
               'Contact element'//trim(getId(cElems(i)%id)),msg1,msg2, &
               'This might cause numerical instabilities. '// &
               'Consider to increase the Thickness.')
       end if

    end do

    if (ierr < lerr) call reportError (debugFileOnly_p,'updateContactElements')

  end subroutine updateContactElements


  !!============================================================================
  !> @brief Updates the local coordinates of the contact element follower.
  !>
  !> @param cElem The contact element to update contact position for
  !> @param cSurf Contact surface associated with @a cElem
  !> @param[in] lpu File unit number for res-file output
  !>
  !> @details The local coordinates are updates with respect to the
  !> contact surface. Then the associated velocities are computed.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !> @date Oct 2001
  !>
  !> @author Knut Morten Okstad
  !> @date Sep 2005

  subroutine updateContactPosition (cElem,cSurf,lpu)

    use ContactSurfaceModule, only : GliderCurveType
    use CurveTypeModule     , only : GetPosCornerSection, GetPosOnCircSection
    use IdTypeModule        , only : getId
    use manipMatrixModule   , only : matmul34
#ifdef FT_DEBUG
    use dbgUnitsModule      , only : dbgCurve, dbgIter, dbgTime
#endif

    type(ContactElementType), intent(inout) :: cElem
    type(GliderCurveType)   , intent(inout) :: cSurf
    integer                 , intent(in)    :: lpu

    !! Local variables
    logical  :: wasActive, wasAtCorner
    integer  :: i, j, iInc, iPrevChk, iPrevDom, lastDomain, nCurvePoints
    real(dp) :: Xsi, XsiPrev, dXPrev(6), CPosPrev(3,4), FPosInG(3,4)
#ifdef FT_DEBUG
    real(dp) :: h, hMax, teta, lenV, secV(3)
#endif

    !! --- Logic section ---

    FPosInG = matmul34(cElem%triad1%ur,cElem%FPosInT) ! global follower position

#ifdef FT_DEBUG
    if (dbgCurve > 0) then
       write(dbgCurve,"(/'===== time =',1pe12.5,' iter =',i3)") dbgTime,dbgIter
       write(dbgCurve,"( '      Contact Element',A)") trim(getId(cElem%id))
       write(dbgCurve,"( '      Follower',1p3e13.5)") FPosInG(:,4)
    end if
#endif

    !! Initialization
    wasActive        = cElem%isActive
    wasAtCorner      = cElem%isAtCorner
    lastDomain       = cElem%lastDomain
    nCurvePoints     = size(cSurf%CPoints)
    cElem%isActive   = .false.
    cElem%isAtCorner = .false.
    cElem%coeff      = 0.0_dp

    !! Check curve section domains, starting with the last one in contact
    XsiPrev  = 0.0_dp
    dXPrev   = 0.0_dp
    CPosPrev = 0.0_dp
    iPrevDom = 0
    iPrevChk = 0
    iInc     = 0
    i        = lastDomain
    do while (.not. cElem%isActive)
       j = 1 + mod(i,nCurvePoints)
#ifdef FT_DEBUG
       if (dbgCurve > 0) write(dbgCurve,"(' ---- Checking domain:',2I4)") i, j
#endif

       call GetPosOnCircSection (FPosInG, cElem%cVarPrev, &
            &                    cSurf%CPoints(i), cSurf%CPoints(j), &
            &                    cElem%thickness, cElem%width, cElem%radius, &
            &                    cElem%isActive, cElem%CPosInG, &
            &                    Xsi, cElem%cVar(:,1))

       if (cElem%isActive) then

          !! We are in contact with current domain
          if (iPrevDom > 0) then

             !! We are in contact with the previously checked domain too,
             !! now choose the one having the smallest dX value
             if (abs(cElem%cVar(1,1)) < abs(dXprev(1))) then
                cElem%lastDomain = i ! current domain is closest
             else
                cElem%lastDomain = iPrevDom ! previous domain is closest
                cElem%CPosInG    = CPosPrev
                cElem%cVar(:,1)  = dXPrev
                Xsi              = XsiPrev
             end if

          else if (cSurf%isLooping .or. nCurvePoints > 2) then

             !! This is the first domain that we detect contact for, we need to
             !! check either the preceeding or the succeeding domain too,
             !! to find out which one actually is the closest
             iPrevDom = i
             i = iPrevChk
             if (Xsi >= 0.5_dp) then
                if (cSurf%isLooping .or. iPrevDom < nCurvePoints-1) then
                   !! Check the succeeding domain next
                   i = 1 + mod(iPrevDom,nCurvePoints)
                end if
             else
                if (cSurf%isLooping .or. iPrevDom > 1) then
                   !! Check the preceeding domain next
                   i = 1 + mod(nCurvePoints+iPrevDom-2,nCurvePoints)
                end if
             end if

             if (i == iPrevChk) then
                !! No more domains to check, so stick with the one we just found
                cElem%lastDomain = iPrevDom
             else
                !! Continue and check domain "i"
                cElem%isActive = .false.
                CPosPrev       = cElem%CPosInG
                dXPrev         = cElem%cVar(:,1)
                XsiPrev        = Xsi
             end if

          end if

       else if (iPrevDom > 0) then

          !! We are in contact with the previously checked domain,
          !! but not with the current one, exit search
          cElem%isActive   = .true.
          cElem%lastDomain = iPrevDom
          cElem%CPosInG    = CPosPrev
          cElem%cVar(:,1)  = dXPrev
          Xsi              = XsiPrev

       else

          !! Not yet in contact, continue to check neighboring domains
          !! Alternate between succeeding and preceeding domains,
          !! moving gradually farther from the starting domain
          iPrevChk = i
          if (iInc <= 0) then
             iInc = 1-iInc ! Try a succeeding domain next
          else
             iInc = -iInc  ! Try a preceeding domain next
          end if
          if (cSurf%isLooping) then
             i = 1 + mod(nCurvePoints+cElem%lastDomain+iInc-1,nCurvePoints)
             if (i == cElem%lastDomain) exit ! Exit search, no contact
          else
             i = cElem%lastDomain + iInc
             if (i >= nCurvePoints) then     ! End of non-looping curve
                iInc = -iInc                 ! Try a preceeding domain next
                i = cElem%lastDomain + iInc
                if (i < 1) exit              ! Exit search, no contact
             else if (i < 1) then            ! Start of non-looping curve
                iInc = 1-iInc                ! Try a succeeding domain next
                i = cElem%lastDomain + iInc
                if (i > nCurvePoints-1) exit ! Exit search, no contact
             end if
          end if

       end if

    end do

    !! Check corner section domains, starting with the one previously in contact
    iInc = 0
    if (cSurf%isLooping) then
       i = cElem%lastDomain
    else
       i = min(max(2,cElem%lastDomain),nCurvePoints-1)
       j = i
    end if
    do while (.not. cElem%isActive)
#ifdef FT_DEBUG
       if (dbgCurve > 0) write(dbgCurve,"(' ---- Checking Corner point:',i4)") i
#endif

       call GetPosCornerSection (FPosInG, cElem%cVarPrev, &
            &                    cSurf%CPoints(i), &
            &                    cElem%thickness, cElem%width, cElem%radius, &
            &                    cElem%isActive, cElem%CPosInG, &
            &                    cElem%cVar(:,1))

       if (cElem%isActive) then

          cElem%lastDomain = i
          cElem%isAtCorner = .true.
          Xsi = 0.0_dp

       else

          !! Not yet in contact, continue to check neighboring domains
          !! Alternate between succeeding and preceeding domains,
          !! moving gradually farther from the starting domain
          if (iInc <= 0) then
             iInc = 1-iInc ! Try a succeeding domain next
          else
             iInc = -iInc  ! Try a preceeding domain next
          end if
          if (cSurf%isLooping) then
             i = 1 + mod(nCurvePoints+cElem%lastDomain+iInc-1,nCurvePoints)
             if (i == cElem%lastDomain) exit ! Exit search, no contact
          else
             i = j + iInc
             if (i >= nCurvePoints) then     ! End of non-looping curve
                iInc = -iInc                 ! Try a preceeding domain next
                i = j + iInc
                if (i < 2) exit              ! Exit search, no contact
             else if (i <= 1) then           ! Start of non-looping curve
                iInc = 1-iInc                ! Try a succeeding domain next
                i = j + iInc
                if (i > nCurvePoints-1) exit ! Exit search, no contact
             end if
          end if

       end if

    end do

    if (cElem%isActive) then

       i = cElem%lastDomain
       if (cElem%isAtCorner) then
          j = i
       else
          j = 1 + mod(i,nCurvePoints)
       end if

       !! Set the nodal coefficients
       cElem%coeff(1) = 1.0_dp - Xsi
       cElem%coeff(2) = Xsi

       !! Update contact point velocities and eccentricity vectors
       call updateContactVelEcc (cElem,cSurf%CPoints(i),cSurf%CPoints(j))

       if (cSurf%isLooping .and. wasActive) then
          !! Adjust the previous slide variable when looping contact curve
          call checkLooping (cSurf%loopLength,cElem%cVar(3,1),cElem%cVarPrev(3))
       end if

       if (.not. wasActive) then
          !! Initiate the stress free length of the friction spring (if any),
          !! such that it is stress free in the first iteration after activation
          if (associated(cElem%friction)) then
             if (associated(cElem%friction%spr)) then
                cElem%friction%spr%l0 = cElem%cVar(3,1)
             end if
          end if
       end if

#ifdef FT_DEBUG
       if (cSurf%CPoints(i)%curvature > epsDiv0_p) then
          secV = cSurf%CPoints(j)%triad%ur(:,4) - cSurf%CPoints(i)%triad%ur(:,4)
          lenV = sqrt(dot_product(secV,secV))
          teta = asin(0.5_dp*lenV*cSurf%CPoints(i)%curvature)
          hMax = (1.0_dp-cos(teta)) / cSurf%CPoints(i)%curvature
       else
          hMax = 0.0_dp
       end if
       h = sqrt(dot_product(cElem%eVec(:,1),cElem%eVec(:,1)))
       if (dbgCurve > 0) then
          write(dbgCurve,"('Found domain i,j, Xsi, dX,dY,dZ :',2I3,1P7E12.3)") &
               i, j, xsi, cElem%cVar(:,1)
          write(dbgCurve,"(41X,'velocity :',1P6E12.3)") cElem%cVar(:,2)
          write(dbgCurve,"('eVecGlSecant :',1P3E11.3,' h,hMax :',2E10.3)") &
               cElem%eVec(:,1),h,hMax
          write(dbgCurve,"('eVecFollower :',1P3E11.3)") &
               cElem%eVec(:,2)
       end if
#endif
    else
       cElem%cVar = 0.0_dp
    end if

    if (cElem%isActive .eqv. wasActive) then
#ifdef FT_DEBUG
       !! In debug mode, report also when follower moves to a new segment
       if (wasAtCorner .eqv. cElem%isAtCorner) then
          if (lastDomain == cElem%lastDomain) return
       end if
#else
       return
#endif
    end if

    !! Output the contact status change to res-file
    if (cElem%isActive) then
       if (wasActive) then
          write(lpu,600) trim(getId(cElem%id)),'changed segment',FPosInG(:,4)
       else
          write(lpu,600) trim(getId(cElem%id)),'activated',FPosInG(:,4)
       end if
#ifdef FT_DEBUG
       write(lpu,610) 1, trim(getId(cSurf%CPoints(i)%triad%id)), &
            &         cSurf%CPoints(i)%triad%ur(:,4), cElem%coeff(1)
       if (cElem%isAtCorner) then
          write(lpu,"(14X,A)") 'Corner point'
       else
          write(lpu,610) 2, trim(getId(cSurf%CPoints(j)%triad%id)), &
               &         cSurf%CPoints(j)%triad%ur(:,4), cElem%coeff(2)
       end if
#endif
    else
       write(lpu,600) trim(getId(cElem%id)),'deactivated',FPosInG(:,4)
    end if

600 format(5X,'Contact Element',A,1X,A,' at follower location',1P3E11.3)
#ifdef FT_DEBUG
610 format(5X,'Point',I2,': Triad',A,1X,1P3E11.3,' coeff =',0PF5.2)
#endif

  contains

    !> @brief Adjusts the previous slide variable when looping contact curve.
    subroutine checkLooping (loopLength,Z,Zprev)
      real(dp), intent(in)    :: loopLength, Z
      real(dp), intent(inout) :: Zprev
      if (Zprev < Z - 0.5_dp*loopLength) then
         Zprev = Zprev + loopLength
      else if (Zprev > Z + 0.5_dp*loopLength) then
         Zprev = Zprev - loopLength
      end if
    end subroutine checkLooping

  end subroutine updateContactPosition


  !!============================================================================
  !> @brief Updates the contact point velocities and eccentricty vectors.
  !>
  !> @param cElem The contact element to update velocity and eccentricity for
  !> @param[in] CP1 First point of the glider curve segment
  !> @param[in] CP2 Second point of the glider curve segment
  !>
  !> @details The contact point velocity is computed from the follower triad.
  !> The eccentricity vectors are computed from the glider triad positions.
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !> @date Oct 2001
  !>
  !> @author Knut Morten Okstad
  !> @date Sep 2005 / 2.0

  subroutine updateContactVelEcc (cElem, CP1, CP2)

    use CurvePointTypeModule, only : CurvePointType

    type(ContactElementType), intent(inout) :: cElem
    type(CurvePointType)    , intent(in)    :: CP1, CP2

    !! Local variables
    integer  :: nDOFs
    real(dp) :: v1(3), secVec(3), velVec(6), lenVec

    !! --- Logic section ---

    !! Update slider variable (curvilinear position along the glider curve)
    if (cElem%isAtCorner) then
       cElem%cVar(3,1) = CP1%slideVar
    else
       if (cElem%cSurf%isLooping .and. CP1%slideVar > CP2%slideVar) then
          lenVec = cElem%cSurf%loopLength
       else
          lenVec = CP2%slideVar
       end if
       cElem%cVar(3,1) = cElem%coeff(1)*CP1%slideVar + cElem%coeff(2)*lenVec
    end if

    !! Define relative velocity of the follower
    !TODO,bh: Base this on velocity of circle section origin?
    velVec = 0.0_dp
    nDOFs  = cElem%triad1%nDOFs
    if (nDOFs > 0) then
       velVec(1:nDOFs) = cElem%triad1%urd
    end if
    nDOFs = CP1%triad%nDOFs
    if (nDOFs > 0) then
       velVec(1:nDOFs) = velVec(1:nDOFs) - cElem%coeff(1)*CP1%triad%urd
    end if
    nDOFs = CP2%triad%nDOFs
    if (nDOFs > 0) then
       velVec(1:nDOFs) = velVec(1:nDOFs) - cElem%coeff(2)*CP2%triad%urd
    end if

    !! Transform to local coordinates
    cElem%cVar(1:3,2) = matmul(velVec(1:3),cElem%CPosInG(:,1:3))
    cElem%cVar(4:6,2) = matmul(velVec(4:6),cElem%CPosInG(:,1:3))

    !! Eccentricity from glider-secant to contact point
    !! NB: normal to secant to minimize moments
    !TODO,bh: Handle position relative to triad
    v1     = cElem%CPosInG(:,4) - CP1%triad%ur(:,4)
    secVec = CP2%triad%ur(:,4)  - CP1%triad%ur(:,4)
    lenVec = sqrt(dot_product(secVec,secVec))
    if (lenVec > epsDiv0_p) then
       secVec = secVec/lenVec
       cElem%eVec(:,1) = v1 - secVec*dot_product(v1,secVec)
    else
       cElem%eVec(:,1) = v1
    end if

    !! Eccentricity from follower triad to contact point
    !TODO,kmo: Handle position relative to triad(?)
    cElem%eVec(:,2) = cElem%CPosInG(:,4) - cElem%triad1%ur(:,4)

!!$!bh, begin test new
!!$    velVec = cElem%triad1%urd(1:3) &
!!$         & - cElem%coeff(1) * CP1%triad%urd(1:3) &
!!$         & - cElem%coeff(2) * CP2%triad%urd(1:3) &
!!$         & - cElem%coeff(1) * cross_product(CP1%triad%urd(4:6), &
!!$         &                                  cElem%eVec(:,1)) &
!!$         & - cElem%coeff(2) * cross_product(CP2%triad%urd(4:6), &
!!$         &                                  cElem%eVec(:,1))
!!$!bh, end test new

  end subroutine updateContactVelEcc


  !!============================================================================
  !> @brief Updates the contact springs and dampers.
  !>
  !> @param cElem The contact element to update
  !> @param[in] timeStep Time increment size
  !> @param ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !> @date Oct 2001
  !>
  !> @author Knut Morten Okstad
  !> @date Sep 2005 / 2.0

  subroutine updateContactSprDmp (cElem,timeStep,ierr)

    use CurveTypeModule     , only : rotateAboutZ
    use DamperTypeModule    , only : updateDamperVelocity
    use EngineRoutinesModule, only : UpdateSpringBase, UpdateDamperBase
    use IdTypeModule        , only : getId
#ifdef FT_DEBUG
    use dbgUnitsModule      , only : dbgCurve
#endif
    use reportErrorModule   , only : reportError, error_p

    type(ContactElementType), intent(inout) :: cElem
    real(dp)                , intent(in)    :: timeStep
    integer                 , intent(inout) :: ierr

    !! Local variables
    integer  :: i, lerr
    real(dp) :: def(2), velDmp(2), velInG(3)

    !! --- Logic section ---

    if (cElem%isActive .and. cElem%nonRadialDamping) then
       !! Contact variables are in the radial system, but the dampers should be
       !! in the local rectangular system, must transform the orientation
       !! matrix and the deflection and velocity variables accordingly
       cElem%DPosInG = cElem%CPosInG(:,1:3)
       call rotateAboutZ (cElem%DPosInG,-cElem%cVar(2,1))
       def(1) = cElem%cVar(1,1) * cos(cElem%cVar(2,1))
       def(2) = cElem%cVar(1,1) * sin(cElem%cVar(2,1))
       velInG = matmul(cElem%CPosInG(:,1:3),cElem%cVar(1:3,2))
       velDmp = matmul(velInG,cElem%DPosInG(:,1:2))
       do i = 1, 2
          if (associated(cElem%dampers(i)%p)) then
             cElem%dampers(i)%p%length   = def(i)
             cElem%dampers(i)%p%velocity = velDmp(i)
          end if
       end do
    end if

    lerr = ierr
    do i = 1, 6
       if (associated(cElem%springs(i)%p)) then

          !! Update spring states
          cElem%springs(i)%p%isActive = cElem%isActive
          call UpdateSpringBase (cElem%springs(i)%p,ierr)
#ifdef FT_DEBUG
          if (dbgCurve > 0) then
             write(dbgCurve,"('    Spring',I2,' force,deflection:',1P2E12.3)") &
                  i, cElem%springs(i)%p%force, cElem%springs(i)%p%deflection
          end if
#endif

       end if
       if (associated(cElem%dampers(i)%p)) then

          !! Update damper states
          cElem%dampers(i)%p%isActive = cElem%isActive
          call updateDamperVelocity (cElem%dampers(i)%p,timeStep)
          call UpdateDamperBase (cElem%dampers(i)%p,ierr)
#ifdef FT_DEBUG
          if (dbgCurve > 0) then
             write(dbgCurve,"('    Damper',I2,' force,velocity:',1P2E12.3)") &
                  i, cElem%dampers(i)%p%force, cElem%dampers(i)%p%velocity
          end if
#endif

       end if
    end do

    if (ierr < lerr) then
       call reportError (error_p,'Failed to update contact element'// &
            &            getId(cElem%id),addString='updateContactSprDmp')
    end if

  end subroutine updateContactSprDmp

end module ContactElementRoutinesModule
