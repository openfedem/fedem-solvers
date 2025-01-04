!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file elStressModule.f90
!> @brief Calculation of element stresses and strains.

!!==============================================================================
!> @brief Module with subroutines for element stress and strain calculation.

module ElstressModule

  use kindModule, only : dp

  implicit none

  logical , parameter :: lStiffProj = .true.    !< Activate stiffness projection
  real(dp), parameter :: sqrt3_p = sqrt(3.0_dp) !< For Gauss-point extrapolation

  private

  public :: extractEV, ElStress


contains

  !!============================================================================
  !> @brief Extracts element displacements from superelement displacements.
  !>
  !> @param[in] IEL Local element index
  !> @param[in] sam Assembly management data for the superelement
  !> @param[in] SV Nodal displacements of the superelement
  !> @param[out] EV Nodal displacements of the element
  !> @param[out] NEDOF Number of DOFs in element
  !>
  !> @details This subroutine does the same as the SAM subroutine EXTEV/EXTEV2,
  !> but needs to be reimplemented to account for different real kind of the
  !> input argument @a SV.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 11 Oct 2000

  subroutine extractEV (IEL,sam,SV,EV,NEDOF)

    use SamModule    , only : SamType
    use recKindModule, only : rk

    integer      , intent(in)  :: IEL
    type(SamType), intent(in)  :: sam
    real(rk)     , intent(in)  :: SV(:)
    real(dp)     , intent(out) :: EV(:)
    integer      , intent(out) :: NEDOF

    !! Local variables
    integer :: iedof, in, ip, js, je, nd, NNDOF

    !! --- Logic section ---

    NEDOF = 0

    select case (sam%melcon(iel))
    case (11,21:24,31:32) ! Beam and shell elements with rotational DOFs
       NNDOF = 6
    case (41:46) ! Continuum elements with translational DOFs only
       NNDOF = 3
    case default
       return ! All other elements have no DOFs
    end select

    do ip = sam%mpmnpc(IEL), sam%mpmnpc(IEL+1)-1
       in = sam%mmnpc(ip)
       js = sam%madof(in)
       nd = sam%madof(in+1) - js
       if (nd > NNDOF) nd = NNDOF
       je = js + nd-1
       iedof = NEDOF + 1
       NEDOF = NEDOF + nd
       if (NEDOF <= size(EV)) then
          EV(iedof:NEDOF) = SV(js:je)
       end if
    end do

    if (NEDOF > size(EV)) NEDOF = size(EV) - NEDOF

  end subroutine extractEV


  !!============================================================================
  !> @brief Computes FE stresses and strains for the specified element.
  !>
  !> @param[in] IEL Local element index
  !> @param[in] IELTYP Element type flag
  !> @param[in] V Nodal displacements of the element
  !> @param[out] S Stress resultants
  !> @param[out] E Stress resultant-conjugate strains
  !> @param[out] Sigma stresses
  !> @param[out] Epsil strains
  !> @param[out] F Element nodal forces
  !> @param[out] NENOD Number of element nodes
  !> @param[out] NSTRP Number of stress calculation points
  !> @param[in] IPSW Print switch
  !> @param[in] LPU File unit number for res-file output
  !> @param[out] IERR Error flag
  !>
  !> @details This subroutine computes beam section forces and shell stress
  !> resultants at element nodes, and stresses and strains at predefined results
  !> points (for beams and shells) or at the element nodes (for solids).
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 11 Oct 2000

  subroutine ElStress (IEL, IELTYP, V, S, E, Sigma, Epsil, F, &
       &               NENOD, NSTRP, IPSW, LPU, IERR)

    use ReportErrorModule      , only : error_p, warning_p, reportError
    use FFlLinkHandlerInterface, only : ffl_getElmId

    real(dp)          , intent(inout) :: V(:)
    real(dp)          , intent(out)   :: S(:), E(:), Sigma(:), Epsil(:)
    real(dp), optional, intent(out)   :: F(:)
    integer           , intent(in)    :: IEL, IELTYP, IPSW, LPU
    integer           , intent(out)   :: NENOD, NSTRP, IERR

    !! Local variables
    integer           :: n, elmNo, iprint
    character(len=64) :: errmsg

    !! --- Logic section ---

    ierr  = 0
    nenod = 0
    nstrp = 0
    elmNo = ffl_getElmId(iel)
    if (elmno < 1) return

    if (IPSW < 0 .and. -IPSW == elmNo) then
       IPRINT = 99 ! Debugging of individual elements
    else
       IPRINT = IPSW
    end if

    select case (IELTYP)

    case (11) ! Two-noded beam
       nenod = 2
       nstrp = 0
       call STR11 (iel,v,F,S,E,Sigma,Epsil,IPRINT-5,LPU,IERR)

    case (21) ! Three-noded thin shell
       nenod = 3
       nstrp = 6
       call STR21 (iel,v,F,S,E,Sigma,Epsil,IPRINT-5,LPU,IERR)

    case (22) ! Four-noded thin shell
       nenod = 4
       nstrp = 8
       call STR22 (iel,v,F,S,E,Sigma,Epsil,IPRINT-5,LPU,IERR)

    case (23) ! Three-noded thin ANDES shell
       nenod = 3
       nstrp = 6
       call STR23 (iel,v,F,S,E,Sigma,Epsil,LPU,IERR)

    case (24) ! Four-noded thin ANDES shell
       nenod = 4
       nstrp = 8
       call STR24 (iel,v,F,S,E,Sigma,Epsil,LPU,IERR)

    case (31) ! Six-noded triangular thick shell
       nenod = 6
       nstrp = 12
       call STR31 (iel,v,S,E,Sigma,Epsil,IPRINT-5,LPU,IERR)

    case (32) ! Eight-noded quadrilateral thick shell
       nenod = 8
       nstrp = 16
       call STR32 (iel,v,S,E,Sigma,Epsil,IPRINT-5,LPU,IERR)

    case (41) ! Ten-noded tetrahedron
       nenod = 10
       nstrp = 10
       call STR41 (iel,v,Sigma,Epsil,IPRINT-5,LPU,IERR)

    case (42) ! Fifteen-noded wedge
       nenod = 15
       nstrp = 15
       call STR42 (iel,v,Sigma,Epsil,IPRINT-5,LPU,IERR)

    case (43) ! Twenty-noded hexahedron
       nenod = 20
       nstrp = 20
       call STR43 (iel,v,Sigma,Epsil,IPRINT-5,LPU,IERR)

    case (44) ! Eight-noded hexahedron
       nenod = 8
       nstrp = 8
       call STR44 (iel,v,Sigma,Epsil,IPRINT-5,LPU,IERR)

    case (45) ! Four-noded tetrahedron
       nenod = 4
       nstrp = 4
       call STR45 (iel,v,Sigma,Epsil,IPRINT-5,LPU,IERR)

    case (46) ! Six-noded wedge
       nenod = 6
       nstrp = 6
       call STR46 (iel,v,Sigma,Epsil,IPRINT-5,LPU,IERR)

    case default ! Silently ignore all the other element types
       return

    end select

    if (ierr /= 0) then
       write(errmsg,"('Stress calculation failed for element',2i8)") iel,elmNo
       if (ierr > 0) then
          call reportError (warning_p,errmsg,'Continuing with the next element')
       else
          call reportError (error_p,errmsg,addString='ElStress')
       end if
       return
    end if

    !! The element routines return the 'engineering' shear strain (gamma_xy)
    !! but on the stress results file we need the tensorial quantity

    do n = 1, nstrp
       select case (IELTYP)
       case (21:24) ! Thin shell elements
          Epsil(3*n) = Epsil(3*n)*0.5_dp
       case (31:32,41:46) ! Thick shell and solid elements
          Epsil(6*n-2) = Epsil(6*n-2)*0.5_dp
          Epsil(6*n-1) = Epsil(6*n-1)*0.5_dp
          Epsil(6*n)   = Epsil(6*n)  *0.5_dp
       end select
    end do

    if (iprint == 0) return

    !! Print out stress results to the log-file

    select case (IELTYP)

    case (11) ! Beam elements
       write(LPU,610) elmNo,(n,S(6*n-5:6*n),n=1,nenod)
       write(LPU,611)       (n,E(6*n-5:6*n),n=1,nenod)
       write(LPU,612)  (n,Sigma(n),Epsil(n),n=1,nstrp)

    case (21:24) ! Thin shell elements
       write(LPU,620) elmNo,(n,S(6*n-5:6*n),n=1,nenod)
       write(LPU,621)       (n,E(6*n-5:6*n),n=1,nenod)
       write(LPU,622)   (n,Sigma(3*n-2:3*n),n=1,nstrp)
       write(LPU,623)   (n,Epsil(3*n-2:3*n),n=1,nstrp)

    case (31:32) ! Thick shell elements
       write(LPU,630) elmNo,(n,Sigma(6*n-5:6*n),n=1,nstrp)
       write(LPU,631)       (n,Epsil(6*n-5:6*n),n=1,nstrp)

    case (41:46) ! Solid elements
       write(LPU,640) elmNo,(n,Sigma(6*n-5:6*n),n=1,nstrp)
       write(LPU,641)       (n,Epsil(6*n-5:6*n),n=1,nstrp)

    case default
       return

    end select

    if (iprint < 2 .or. .not.present(F)) return

    select case (IELTYP)
    case (41:46)
       write(LPU,642) elmNo,(n,F(3*n-2:3*n),n=1,nenod)
    case default
       write(LPU,643) elmNo,(n,F(6*n-5:6*n),n=1,nenod)
    end select

610 format(//5X,'SECTIONAL FORCES AND STRESSES FOR BEAM ELEMENT',I8,' :' &
         & //5X,'Node #',4X,'Nx',11X,'Vy',11X,'Vz',11X,'Mx',11X,'My',11X,'Mz' &
         &  /1P,(I9,2X,6E13.5))
611 format( /5X,'Node #',4X,'epsilon_x',5X,'gamma_y',6X,'gamma_z', &
         &               6X,'kappa_x',6X,'kappa_y',6X,'kappa_z' &
         &  /1P,(I9,2X,6E13.5))
612 format( /5X,'Point#',4X,'sigma_xx',2X,'epsilon_xx' &
         &  /1P,(I9,2X,2E13.5))

620 format(//5X,'STRESS RESULTANTS AND STRESSES FOR SHELL ELEMENT',I8,' :' &
         & //5X,'Node #',4X,'Nxx',10X,'Nyy',10X,'Nxy', &
         &              10X,'Mxx',10X,'Myy',10X,'Mxy' &
         &  /1P,(I9,2X,6E13.5))
621 format( /5X,'Node #',3X,'epsilon_xx',3X,'epsilon_yy',3X,'epsilon_xy', &
         &               4X,'kappa_xx',5X,'kappa_yy',5X,'kappa_xy' &
         &  /1P,(I9,2X,6E13.5))
622 format( /5X,'Point#',4X,'sigma_xx',5X,'sigma_yy',5X,'sigma_xy' &
         &  /1P,(I9,2X,3E13.5))
623 format( /5X,'Point#',3X,'epsilon_xx',3X,'epsilon_yy',3X,'epsilon_xy' &
         &  /1P,(I9,2X,3E13.5))

630 format(//5X,'STRESSES FOR SHELL ELEMENT',I8,' :' &
         & //5X,'Point#',4X,'sigma_xx',5X,'sigma_yy',5X,'sigma_zz', &
         &               5X,'sigma_xy',5X,'sigma_xz',5X,'sigma_yz' &
         &  /1P,(I9,2X,6E13.5))
631 format( /5X,'Point#',3X,'epsilon_xx',3X,'epsilon_yy',3X,'epsilon_zz', &
         &               3X,'epsilon_xy',3X,'epsilon_xz',3X,'epsilon_yz' &
         &  /1P,(I9,2X,6E13.5))

640 format(//5X,'STRESSES FOR SOLID ELEMENT',I8,' :' &
         & //5X,'Node #',4X,'sigma_xx',5X,'sigma_yy',5X,'sigma_zz', &
         &               5X,'sigma_xy',5X,'sigma_xz',5X,'sigma_yz' &
         &  /1P,(I9,2X,6E13.5))
641 format( /5X,'Node #',3X,'epsilon_xx',3X,'epsilon_yy',3X,'epsilon_zz', &
         &               3X,'epsilon_xy',3X,'epsilon_xz',3X,'epsilon_yz' &
         &  /1P,(I9,2X,6E13.5))

642 format(//5X,'NODAL FORCES FOR SOLID ELEMENT',I8,' :' &
         & //5X,'Node #',6X,'F_x',10X,'F_y',10X,'F_z' &
         &  /1P,(I9,2X,3E13.5))
643 format(//5X,'NODAL FORCES FOR ELEMENT',I8,' :' &
         & //5X,'Node #',6X,'F_x',10X,'F_y',10X,'F_z', &
         &              10X,'M_x',10X,'M_y',10X,'M_z' &
         &  /1P,(I9,2X,6E13.5))

  end subroutine ElStress


  !!============================================================================
  !> @brief Utility to retrieve element data needed for stress recovery.
  !>
  !> @param[in] PRG Name of calling subroutine
  !> @param[in] IEL Local element index
  !> @param[out] XG Global X-coordinates of the element nodes
  !> @param[out] YG Global Y-coordinates of the element nodes
  !> @param[out] ZG Global Z-coordinates of the element nodes
  !> @param[out] EMOD Youngs modulus
  !> @param[out] RNY Poisson's ratio
  !> @param[out] THK Nodal shell thicknesses
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 29 Jan 2021

  subroutine getElementData (PRG,IEL,XG,YG,ZG,EMOD,RNY,THK,ierr)

    use FFlLinkHandlerInterface, only : ffl_getcoor, ffl_getmat, ffl_getthick
    use ReportErrorModule      , only : internalError

    character(len=*)  , intent(in)  :: PRG
    integer           , intent(in)  :: IEL
    real(dp)          , intent(out) :: XG(:), YG(:), ZG(:), EMOD, RNY
    real(dp), optional, intent(out) :: THK(:)
    integer           , intent(out) :: ierr

    !! Local variables
    real(dp) :: dummyRho

    !! --- Logic section ---

    call ffl_getcoor (XG,YG,ZG,IEL,ierr)
    if (ierr < 0) then
       ierr = internalError(PRG//': Retrieving element nodal coordinates')
       return
    end if

    call ffl_getmat (EMOD,RNY,dummyRho,IEL,ierr)
    if (ierr < 0) then
       ierr = internalError(PRG//': Retrieving element material data')
       return
    end if

    if (present(THK)) then
       call ffl_getthick (THK,IEL,ierr)
       if (ierr < 0) then
          ierr = internalError(PRG//': Retrieving shell element thickness')
       end if
    end if

  end subroutine getElementData


  !> @cond FULL_DOC
  !!****************************************************************************
  !> @brief Computes FE stresses and strains for a 2-noded beam element.

  subroutine STR11 (iel,EV,EF,SF,SS,sigma,epsil,IPSW,LPU,IERR)

    use ReportErrorModule      , only : warning_p, reportError, internalError
    use FFlLinkHandlerInterface, only : ffl_getCoor
    use FFlLinkHandlerInterface, only : ffl_getBeamSection, ffl_getPinFlags

    integer, parameter :: nedof = 12

    integer , intent(in)            :: iel, IPSW, LPU
    real(dp), intent(in)            :: EV(nedof) !< Element nodal displacements
    real(dp), intent(out), optional :: EF(nedof) !< Element nodal forces
    real(dp), intent(out)           :: SF(6,2)   !< Element nodal section forces
    real(dp), intent(out)           :: SS(6,2)   !< Force-conjugate strains
    real(dp), intent(out)           :: sigma(1)  !< Stresses at result points
    real(dp), intent(out)           :: epsil(1)  !< Strains at result points
    integer , intent(out)           :: IERR

    !! Local variables
    integer  :: IPA, IPB
    real(dp) :: XG(5), YG(5), ZG(5), BSEC(14)
    real(dp) :: Bl, EK(nedof,nedof), Ex, Ey, Ez, Sg(6), SM(3), SN(3), T(3,3)

    !! --- Logic section ---

    call ffl_getCoor (XG,YG,ZG,IEL,ierr)
    if (ierr < 0) then
       ierr = internalError('STR11: Retrieving element nodal coordinates')
       return
    end if

    call ffl_getBeamSection (BSEC,IEL,ierr)
    if (ierr < 0) then
       ierr = internalError('STR11: Retrieving beam element data')
       return
    end if

    call ffl_getPinFlags (IPA,IPB,IEL,ierr)
    if (ierr < 0) then
       ierr = internalError('STR11: Retrieving beam pin flags')
       return
    end if

    !! --- Compute the element stiffness matrix

    call BEAM31 (EK(1,1),XG(1),YG(1),ZG(1),BSEC(2),BSEC(9),BSEC(11),BSEC(13), &
         &       BSEC(14),IPA,IPB,IEL,LPU,IPSW,IERR)
    if (ierr /= 0) then
       ierr = 1
       call reportError (warning_p,'Failure calculating beam stiffness matrix')
       return
    end if

    !! --- Calculate the element nodal forces

    if (present(EF)) then
       EF = matmul(EK,EV)
    end if

    !! --- Calculate the local sectional forces

    ! Global forces at local node 1
    if (present(EF)) then
       Sg = EF(1:6)
    else
       Sg = matmul(EK(1:6,:),EV)
    end if

    ! Eccentricity at local node 1
    Ex = Xg(4) - Xg(1)
    Ey = Yg(4) - Yg(1)
    Ez = Zg(4) - Zg(1)

    ! Transform eccentric moments to element end 1
    Sg(4) = Sg(4) - Ez*Sg(2) + Ey*Sg(3)
    Sg(5) = Sg(5) + Ez*Sg(1) - Ex*Sg(3)
    Sg(6) = Sg(6) - Ey*Sg(1) + Ex*Sg(2)

    ! Transform to local coordinates
    call DCOS30(T,Xg,Yg,Zg)
    SN = matmul(T,Sg(1:3))
    SM = matmul(T,Sg(4:6))

    ! Sectional forces at element end 1
    SF(1,1) = -SN(1)
    SF(2,1) =  SN(2)
    SF(3,1) =  SN(3)
    SF(4,1) = -SM(1) + BSEC(10)*SN(3) - BSEC(11)*SN(2)
    SF(5,1) =  SM(2)
    SF(6,1) =  SM(3)

    ! Element length
    Ex = XG(2) - XG(1)
    Ey = YG(2) - YG(1)
    Ez = ZG(2) - ZG(1)
    Bl = sqrt(Ex*Ex + Ey*Ey + Ez*Ez)

    ! Sectional forces at element end 2
    SF(:,2) = SF(:,1)
    SF(5,2) = SF(5,2) + SN(3)*Bl
    SF(6,2) = SF(6,2) + SN(2)*Bl

    ! Sectional force conjugate strains ...  later !

    SS = 0.0_dp

    !! --- Calculate element stresses and conjugate strains

    ! Compute axial stress and strain at result points ...  later !

    IERR = 0
    Sigma(1) = 0.0_dp
    Epsil(1) = 0.0_dp

  end subroutine STR11


  !!****************************************************************************
  !> @brief Computes FE stresses and strains for a 3-noded shell element.

  subroutine STR21 (iel,EV,EF,SR,SS,sigma,epsil,IPSW,LPU,IERR)

    use IsoMatModule                  , only : isoMat2D, isoMat2Dinv
    use StrainAndStressUtilitiesModule, only : getShellElementAxes
    use StrainAndStressUtilitiesModule, only : getShellStressTrans
    use ReportErrorModule             , only : warning_p, reportError
    use FFaTensorTransformsInterface  , only : tratensor
    use FFaCmdLineArgInterface        , only : ffa_cmdlinearg_getint

    integer, parameter :: nenod = 3, nedof = 6*nenod

    integer , intent(in)            :: iel, IPSW, LPU
    real(dp), intent(in)            :: EV(nedof) !< Element nodal displacements
    real(dp), intent(out), optional :: EF(nedof) !< Element nodal forces
    real(dp), intent(out) :: SR(6,nenod)         !< Stress resultants at nodes
    real(dp), intent(out) :: SS(6,nenod)         !< SR-conjugate strains
    real(dp), intent(out) :: sigma(3,2*nenod)    !< Stresses at result points
    real(dp), intent(out) :: epsil(3,2*nenod)    !< Strains at result points
    integer , intent(out) :: IERR

    !! Local variables
    real(dp), parameter :: alpha = 1.5_dp ! Hardcoded from the old EMAT
    real(dp), parameter :: beta  = 0.5_dp ! Hardcoded from the old EMAT

    integer  :: i, stressForm
    real(dp) :: THK(nenod), XG(nenod), YG(nenod), ZG(nenod), EK(nedof,nedof)
    real(dp) :: VBL(9), VML(9), VX(3), VZ(3), T_str(2,2), EMOD, RNY, E(3,3)
    real(dp) :: AKB(9,9), AKM(7,9), ESMB(3,9), ESMM(3,9), RBF(3), RMF(3)

    !! --- Logic section ---

    call getElementData ('STR21',IEL,XG,YG,ZG,EMOD,RNY,THK,IERR)
    if (ierr < 0) return

    call isoMat2D (EMOD,RNY,E)

    !! --- Compute the element stiffness matrix

    call ffa_cmdlinearg_getint ('fftStressForm',stressForm)
    if (stressForm == 1) then
       CALL FTSA31 (EK(1,1),AKM(1,1),AKB(1,1),E(1,1), &
            &       XG(1),YG(1),ZG(1),THK(1),1,IEL,LPU,IPSW,IERR)
    else
       CALL FTS31 (EK(1,1),AKM(1,1),AKB(1,1),E(1,1), &
            &      XG(1),YG(1),ZG(1),THK(1),ALPHA,BETA,IEL,LPU,IPSW,IERR)
    end if
    if (ierr /= 0) then
       ierr = 1
       call reportError (warning_p,'Failure calculating shell stiffness matrix')
       return
    end if

    !! --- Calculate the element nodal forces

    if (present(EF)) then
       EF = matmul(EK,EV)
    end if

    !! --- Compute the 2D stress transformation matrix
    call getShellElementAxes (nenod,XG,YG,ZG,VX,VML,VZ,LPU,ierr)
    if (ierr /= 0) return

    call getShellStressTrans (VX,VZ,T_str,LPU,ierr)
    if (ierr /= 0) return

    !! --- Calculate stress resultants at the element center

    if (stressForm == 1) then
       CALL FTSA32 (ESMM(1,1),ESMB(1,1),AKM(1,1),AKB(1,1),XG(1),YG(1),ZG(1))
    else
       CALL FTS32 (ESMM(1,1),ESMB(1,1),AKM(1,1),AKB(1,1),XG(1),YG(1),ZG(1), &
            &      E(1,1),ALPHA,IEL,LPU,IERR)
    end if
    if (ierr /= 0) then
       ierr = 1
       call reportError (warning_p,'Failure calculating shell stress matrices')
       return
    end if

    CALL FTS38 (VML(1),VBL(1),EV(1),XG(1),YG(1),ZG(1))
    RMF = matmul(ESMM,VML)
    RBF = matmul(ESMB,VBL)

    !! --- Transform to stress output coordinate system
    call tratensor (2,RMF,T_str)
    call tratensor (2,RBF,T_str)

    !! --- Nodal stress resultants

    SR(1:3,1) = RMF
    SR(4:6,1) = RBF
    SR(:,2) = SR(:,1)
    SR(:,3) = SR(:,1)

    !! --- Nodal stress resultant conjugate strains ...  later !

    SS = 0.0_dp

    !! --- Nodal stresses

    do i = 1, nenod
       !! Top surface
       sigma(:,i)       = (RMF + RBF*6.0_dp/THK(i))/THK(i)
       !! Bottom surface
       sigma(:,nenod+i) = (RMF - RBF*6.0_dp/THK(i))/THK(i)
    end do

    !! --- Nodal strains

    call isoMat2Dinv (EMOD,RNY,E)
    epsil = matmul(E,sigma)

  end subroutine STR21


  !!****************************************************************************
  !> @brief Computes FE stresses and strains for a 4-noded shell element.

  subroutine STR22 (iel,EV,EF,SR,SS,sigma,epsil,IPSW,LPU,IERR)

    use IsoMatModule                  , only : isoMat2D
    use PmatModule                    , only : pMatStiff
    use StrainAndStressUtilitiesModule, only : getShellElementAxes
    use StrainAndStressUtilitiesModule, only : getShellStressTrans
    use ReportErrorModule             , only : warning_p, reportError
    use FFaCmdLineArgInterface        , only : ffa_cmdlinearg_getint

    integer, parameter :: nenod = 4, nedof = 6*nenod

    integer , intent(in)            :: iel, IPSW, LPU
    real(dp), intent(inout)         :: EV(nedof) !< Element nodal displacements
    real(dp), intent(out), optional :: EF(nedof) !< Element nodal forces
    real(dp), intent(out) :: SR(6,nenod)         !< Stress resultants at nodes
    real(dp), intent(out) :: SS(6,nenod)         !< SR-conjugate strains
    real(dp), intent(out) :: sigma(3,2*nenod)    !< Stresses at result points
    real(dp), intent(out) :: epsil(3,2*nenod)    !< Strains at result points
    integer , intent(out) :: IERR

    !! Local variables
    real(dp), parameter :: alpha = 1.5_dp ! Hardcoded from the old EMAT
    real(dp), parameter :: beta  = 0.9_dp ! Hardcoded from the old EMAT

    integer  :: stressForm
    real(dp) :: THK(nenod), XG(nenod), YG(nenod), ZG(nenod)
    real(dp) :: EMOD, RNY, E(3,3), T_el(3,3), T_str(2,2)
    real(dp) :: PMAT(nedof,nedof), EK(nedof,nedof)

    !! --- Logic section ---

    call getElementData ('STR22',IEL,XG,YG,ZG,EMOD,RNY,THK,IERR)
    if (ierr < 0) return

    call isoMat2D (EMOD,RNY,E)

    if (lStiffProj) then

       !! Compute projection matrix
       call pMatStiff (xg,yg,zg,PMAT,LPU,IERR)
       if (IERR < 0) then
          IERR = 1
          call reportError (warning_p,'Failure computing projection matrix')
          return
       end if

       !! Compute deformational displacement vector through projection
       ev = matmul(pmat,ev)

    end if
    if (present(EF)) then

       !! Compute the element stiffness matrix
       CALL FQS31 (EK(1,1),E(1,1),XG(1),YG(1),ZG(1),THK(1),ALPHA,BETA, &
            &      IEL,LPU,IPSW,IERR,EMOD,RNY)
       IF (IERR < 0) THEN
          IERR = 1
          call reportError (warning_p,'Failure computing stiffness matrix')
          return
       END IF

       !! Compute the element nodal forces
       if (lStiffProj) then
          EF = matmul(transpose(PMAT),matmul(EK,EV))
       else
          EF = matmul(EK,EV)
       end if

    end if

    !! Compute the global-to-local transformation matrix
    call getShellElementAxes (nenod,XG,YG,ZG,T_el(1,:),T_el(2,:),T_el(3,:), &
         &                    LPU,ierr)
    if (ierr /= 0) return

    !! Compute the 2D stress transformation matrix
    call getShellStressTrans (T_el(1,:),T_el(3,:),T_str,LPU,ierr)
    if (ierr /= 0) return

    call ffa_cmdlinearg_getint ('ffqStressForm',stressForm)
    if (stressForm > 0) then

       !! Extrapolation from Gauss points (either 1x1 or 2x2)
       call STR22a (stressForm,XG,YG,ZG,THK,E,T_el,T_str,EV, &
            &       SR,SS,sigma,epsil)

    else

       !! Direct evaluation in nodes (using Femlib routines)
       call STR22b (EMOD,RNY,THK,XG,YG,ZG,E,T_str,EV, &
            &       SR,SS,sigma,epsil,IEL,IPSW,LPU,IERR)
       if (IERR < 0) IERR = 1

    end if

  end subroutine STR22

  !!****************************************************************************
  !> @brief Computes FE stresses and strains for a 4-noded shell element.

  subroutine STR22a (nGauss,XG,YG,ZG,THK,Cmat,T_el,T_str,EV,SR,SS,sigma,epsil)

    use StrainAndStressUtilitiesModule, only : StrainDispQuad4
    use FFaTensorTransformsInterface  , only : tratensor

    integer , parameter   :: nenod = 4, nedof = 6*nenod
    integer , intent(in)  :: nGauss
    real(dp), intent(in)  :: XG(:), YG(:), ZG(:), THK(:)
    real(dp), intent(in)  :: Cmat(:,:), T_el(:,:), T_str(:,:), EV(:)
    real(dp), intent(out) :: SR(:,:), SS(:,:), sigma(:,:), epsil(:,:)

    !! Local variables
    integer  :: i, j, i1, i2, j1, j2, n
    real(dp) :: hHalf, vld(nedof), B_L(3,nedof), B_U(3,nedof)
    real(dp) :: epsGU(3,2,2), epsGL(3,2,2), epsU(3), epsL(3)
    real(dp) :: sigU(3), sigL(3), sigBen(3), sigMem(3), gauss(2)

    integer, parameter :: iClose(4) = (/1,2,2,1/)
    integer, parameter :: iFar(4)   = (/2,1,1,2/)
    integer, parameter :: jClose(4) = (/1,1,2,2/)
    integer, parameter :: jFar(4)   = (/2,2,1,1/)

    !! Factors to extrapolate from Gauss points to nodal values
    real(dp), parameter :: f1 = 0.5_dp + 0.5_dp*sqrt3_p
    real(dp), parameter :: f2 = 0.5_dp - 0.5_dp*sqrt3_p

    !! --- Logic section ---

    if (nGauss == 1) then
       gauss(1) = 0.0_dp
       gauss(2) = 1.0_dp
    else
       gauss(1) = -1.0_dp / sqrt3_p
       gauss(2) =  1.0_dp / sqrt3_p
    end if

    hHalf = sum(THK) / real(2*nenod)
    do i = 1, nGauss
       do j = 1, nGauss

          call StrainDispQuad4 (6,XG,YG,ZG, T_el, gauss(i),gauss(j), hHalf, B_U)
          call StrainDispQuad4 (6,XG,YG,ZG, T_el, gauss(i),gauss(j),-hHalf, B_L)

          ! Form local coordinate displacement vector
          if (i == 1 .and. j == 1) then
             do i1 = 1, 8
                j1 = 3*i1-2
                j2 = 3*i1
                vld(j1:j2) = matmul(T_el,ev(j1:j2))
             end do
          end if

          ! Strain at gauss point, upper and lower surface
          epsGU(:,i,j) = matmul(B_U,vld)
          epsGL(:,i,j) = matmul(B_L,vld)

          ! Transform to stress output coordinate system
          call traStrain (epsGU(:,i,j),T_str)
          call traStrain (epsGL(:,i,j),T_str)

       end do
    end do

    do n = 1, nenod

       ! Extrapolate gauss-point strains to element node
       if (nGauss == 1) then
          epsU = epsGU(:,1,1)
          epsL = epsGL(:,1,1)
       else
          i1 = iClose(n)
          j1 = jClose(n)
          i2 = iFar(n)
          j2 = jFar(n)
          epsU = f1*epsGU(:,i1,j1) + f2*epsGU(:,i2,j2)
          epsL = f1*epsGL(:,i1,j1) + f2*epsGL(:,i2,j2)
       end if

       epsil(:,n)       = epsU
       epsil(:,nenod+n) = epsL

       ! Extrapolated nodal stresses
       sigU = matmul(Cmat,epsU)
       sigL = matmul(Cmat,epsL)

       sigma(:,n)       = sigU
       sigma(:,nenod+n) = sigL

       sigMem = (sigU + sigL)*0.5_dp ! Membrane part of stress
       sigBen = (sigU - sigL)*0.5_dp ! Bending part of stress

       ! Nodal stress resultants
       SR(1:3,n) = sigMem*thk(i)
       SR(4:6,n) = sigBen*thk(i)*thk(i)/6.0_dp

       ! Nodal stress resultant conjugate strains ...  later !
       SS(:,n) = 0.0_dp

    end do

  contains

    !> @brief Transforms a 2D strain tensor.
    subroutine traStrain (eps,T_str)
      real(dp), intent(inout) :: eps(3)
      real(dp), intent(in)    :: T_str(2,2)
      eps(3) = 0.5_dp * eps(3) ! to tensorial shear strain, epsilon_xy
      call tratensor (2,eps,T_str)
      eps(3) = 2.0_dp * eps(3) ! back to enginering shear strain, gamma_xy
    end subroutine traStrain

  end subroutine STR22a

  !!****************************************************************************
  !> @brief Computes FE stresses and strains for a 4-noded shell element.

  subroutine STR22b (EMOD,RNY,TH,XG,YG,ZG,Cmat,T_str,EV, &
       &             SR,SS,sigma,epsil,IEL,IPSW,IW,IERR)

    use manipMatrixModule           , only : invert33
    use FFaTensorTransformsInterface, only : tratensor

    integer , parameter   :: nenod = 4
    real(dp), intent(in)  :: EMOD, RNY, TH(:), XG(:), YG(:), ZG(:)
    real(dp), intent(in)  :: Cmat(:,:), T_str(:,:), EV(:)
    real(dp), intent(out) :: SR(:,:), SS(:,:), sigma(:,:), epsil(:,:)
    integer , intent(in)  :: IEL, IPSW, IW
    integer , intent(out) :: IERR

    !! Local variables
    integer  :: n, m
    real(dp) :: Cinv(3,3)

    !! --- Logic section ---

    CALL FQS32 (EMOD,RNY,Cmat(1,1),TH(1),XG(1),YG(1),ZG(1),EV(1), &
         &      SIGMA(1,1),SR(1,1),IEL,IW,IPSW,IERR)
    if (IERR < 0) return

    do n = 1, nenod
       m = nenod + n

       ! Transform to stress output coordinate system
       call tratensor (2,sigma(:,n),T_str)
       call tratensor (2,sigma(:,m),T_str)
       call tratensor (2,SR(1:3,n),T_str)
       call tratensor (2,SR(4:6,n),T_str)

    end do

    ! Stress-conjugate strains
    Cinv = invert33(Cmat,IW,IERR)
    if (IERR == 0) epsil = matmul(Cinv,sigma)

    ! Nodal stress resultant conjugate strains ...  later !
    SS = 0.0_dp

  end subroutine STR22b


  !!****************************************************************************
  !> @brief Computes FE stresses and strains for a 3-noded ANDES shell element.

  subroutine STR23 (iel,EV,EF,SR,SS,sigma,epsil,LPU,IERR)

    use IsoMatModule                  , only : isoMat2D, isoMat2Dinv
    use StrainAndStressUtilitiesModule, only : getShellElementAxes
    use StrainAndStressUtilitiesModule, only : getShellStressTrans
    use ReportErrorModule             , only : warning_p, reportError
    use FFaTensorTransformsInterface  , only : tratensor
    use FFaCmdLineArgInterface        , only : ffa_cmdlinearg_getint

    integer, parameter :: nenod = 3, nedof = 6*nenod

    integer , intent(in)            :: iel, lpu
    real(dp), intent(in)            :: EV(nedof) !< Element nodal displacements
    real(dp), intent(out), optional :: EF(nedof) !< Element nodal forces
    real(dp), intent(out) :: SR(6,nenod)         !< Stress resultants at nodes
    real(dp), intent(out) :: SS(6,nenod)         !< SR-conjugate strains
    real(dp), intent(out) :: sigma(3,2*nenod)    !< Stresses at result points
    real(dp), intent(out) :: epsil(3,2*nenod)    !< Strains at result points
    integer , intent(out) :: IERR

    !! Local variables
    integer  :: i
    real(dp) :: THK(nenod), XG(nenod), YG(nenod), ZG(nenod), EK(nedof,nedof)
    real(dp) :: VBL(9), VML(9), VX(3), VZ(3), T_str(2,2), EMOD, RNY, E(3,3)
    real(dp) :: AKB(9,9), AKM(7,9), ESMB(3,9), ESMM(3,9), RBF(3), RMF(3)

    !! --- Logic section ---

    call getElementData ('STR23',IEL,XG,YG,ZG,EMOD,RNY,THK,IERR)
    if (ierr < 0) return

    call isoMat2D (EMOD,RNY,E)

    !! --- Compute the element stiffness matrix

    CALL FTSA31 (EK(1,1),AKM(1,1),AKB(1,1),E(1,1), &
         &       XG(1),YG(1),ZG(1),THK(1),1,IEL,LPU,0,IERR)
    if (ierr /= 0) then
       ierr = 1
       call reportError (warning_p,'Failure calculating shell stiffness matrix')
       return
    end if

    !! --- Calculate the element nodal forces

    if (present(EF)) then
       EF = matmul(EK,EV)
    end if

    !! --- Compute the 2D stress transformation matrix
    call getShellElementAxes (nenod,XG,YG,ZG,VX,VML,VZ,LPU,ierr)
    if (ierr /= 0) return

    call getShellStressTrans (VX,VZ,T_str,LPU,ierr)
    if (ierr /= 0) return

    !! --- Calculate stress resultants at the element center

    CALL FTSA32 (ESMM(1,1),ESMB(1,1),AKM(1,1),AKB(1,1),XG(1),YG(1),ZG(1))
    if (ierr /= 0) then
       ierr = 1
       call reportError (warning_p,'Failure calculating shell stress matrices')
       return
    end if

    CALL FTS38 (VML(1),VBL(1),EV(1),XG(1),YG(1),ZG(1))
    RMF = matmul(ESMM,VML)
    RBF = matmul(ESMB,VBL)

    !! --- Transform to stress output coordinate system
    call tratensor (2,RMF,T_str)
    call tratensor (2,RBF,T_str)

    !! --- Nodal stress resultants

    SR(1:3,1) = RMF
    SR(4:6,1) = RBF
    SR(:,2) = SR(:,1)
    SR(:,3) = SR(:,1)

    !! --- Nodal stress resultant conjugate strains ...  later !

    SS = 0.0_dp

    !! --- Nodal stresses

    do i = 1, nenod
       !! Top surface
       sigma(:,i)       = (RMF + RBF*6.0_dp/THK(i))/THK(i)
       !! Bottom surface
       sigma(:,nenod+i) = (RMF - RBF*6.0_dp/THK(i))/THK(i)
    end do

    !! --- Nodal strains

    call isoMat2Dinv (EMOD,RNY,E)
    epsil = matmul(E,sigma)

  end subroutine STR23


  !!****************************************************************************
  !> @brief Computes FE stresses and strains for a 4-noded ANDES shell element.

  subroutine STR24 (iel,EV,EF,SR,SS,sigma,epsil,LPU,IERR)

    use IsoMatModule                  , only : isoMat2D
    use PmatModule                    , only : pMatStiff
    use StrainAndStressUtilitiesModule, only : getShellElementAxes
    use StrainAndStressUtilitiesModule, only : getShellStressTrans
    use ReportErrorModule             , only : warning_p, reportError

    integer, parameter :: nenod = 4, nedof = 6*nenod

    integer , intent(in)            :: iel, lpu
    real(dp), intent(inout)         :: EV(nedof) !< Element nodal displacements
    real(dp), intent(out), optional :: EF(nedof) !< Element nodal forces
    real(dp), intent(out) :: SR(6,nenod)         !< Stress resultants at nodes
    real(dp), intent(out) :: SS(6,nenod)         !< SR-conjugate strains
    real(dp), intent(out) :: sigma(3,2*nenod)    !< Stresses at result points
    real(dp), intent(out) :: epsil(3,2*nenod)    !< Strains at result points
    integer , intent(out) :: IERR

    !! Local variables
    real(dp), parameter :: alpha = 1.5_dp ! Hardcoded from the old EMAT
    real(dp), parameter :: beta  = 0.9_dp ! Hardcoded from the old EMAT

    real(dp) :: THK(nenod), XG(nenod), YG(nenod), ZG(nenod)
    real(dp) :: EMOD, RNY, Cmat(3,3), T_el(3,3), T_str(2,2)
    real(dp) :: PMAT(nedof,nedof), EK(nedof,nedof)

    !! --- Logic section ---

    call getElementData ('STR24',IEL,XG,YG,ZG,EMOD,RNY,THK,IERR)
    if (ierr < 0) return

    call isoMat2D (EMOD,RNY,Cmat)

    !! Compute projection matrix
    call pMatStiff (xg,yg,zg,PMAT,LPU,IERR)
    if (IERR < 0) then
       IERR = 1
       call reportError (warning_p,'Failure computing projection matrix')
       return
    end if

    !! Compute deformational displacement vector through projection
    ev = matmul(pmat,ev)

    if (present(EF)) then

       !! Compute the element stiffness matrix
       !! TODO: Replace by real ANDES stiffness
       CALL FQS31 (EK(1,1),Cmat(1,1),XG(1),YG(1),ZG(1),THK(1),ALPHA,BETA, &
            &      IEL,LPU,0,IERR,EMOD,RNY)
       IF (IERR < 0) THEN
          IERR = 1
          call reportError (warning_p,'Failure computing stiffness matrix')
          return
       END IF

       !! Compute the element nodal forces
       EF = matmul(transpose(PMAT),matmul(EK,EV))

    end if

    !! Compute the global-to-local transformation matrix
    call getShellElementAxes (nenod,XG,YG,ZG,T_el(1,:),T_el(2,:),T_el(3,:), &
         &                    LPU,ierr)
    if (ierr /= 0) return

    !! Compute the 2D stress transformation matrix
    call getShellStressTrans (T_el(1,:),T_el(3,:),T_str,LPU,ierr)
    if (ierr /= 0) return

    !! Extrapolation from 2x2 Gauss points
    call STR22a (2,XG,YG,ZG,THK,Cmat,T_el,T_str,EV,SR,SS,sigma,epsil)

  end subroutine STR24


  !!****************************************************************************
  !> @brief Computes FE stresses and strains for a 6-noded shell element.

  subroutine STR31 (iel,EV,SR,SS,sigma,epsil,IPSW,LPU,IERR)

    use IsoMatModule                , only : isoMat2Dinv
    use FFaTensorTransformsInterface, only : tratensor
    use ReportErrorModule           , only : warning_p, reportError

    integer, parameter :: nenod = 6, nedof = 5*nenod

    integer , intent(in)    :: iel, IPSW, LPU
    real(dp), intent(inout) :: EV(6*nenod)      !< Element nodal displacements
    real(dp), intent(out)   :: SR(6,nenod)      !< Stress resultants at nodes
    real(dp), intent(out)   :: SS(6,nenod)      !< SR-conjugate strains
    real(dp), intent(out)   :: sigma(6,2*nenod) !< Stresses at result points
    real(dp), intent(out)   :: epsil(6,2*nenod) !< Strains at result points
    integer , intent(out)   :: IERR

    !! Local variables
    integer  :: i, k, ip, ICODIR, LIN
    real(dp) :: THK(nenod), XG(nenod), YG(nenod), ZG(nenod)
    real(dp) :: EMOD, RNY, L1(3), L2(3), ZETA, RL1(nenod), RL2(nenod)
    real(dp) :: LAMBI(nenod,3,3), LAMP(3,3), SIG(5,5*nenod), Einv(5,5)

    !! --- Logic section ---

    call getElementData ('STR31',IEL,XG,YG,ZG,EMOD,RNY,THK,ierr)
    if (ierr < 0) return

    !! Set up the inverse of the constitutive matrix defined in SCTS32
    call isoMat2Dinv (EMOD,RNY,Einv)
    Einv(4,4) = Einv(3,3)*1.2_dp
    Einv(5,5) = Einv(4,4)

    ICODIR = 1
    LIN    = 1
    call SCTS30 (XG,YG,ZG,IEL,ICODIR,RL1,RL2,LAMBI,IPSW,LPU,IERR)
    if (IERR < 0) then
       IERR = 1
       call reportError (warning_p,'Failure computing nodal direction cosines')
       return
    end if

    !! Transform nodal rotations to local directions
    do i = 1, nenod
       EV(5*i-4:5*i-2) = EV(6*i-5:6*i-3)
       EV(5*i-1:5*i) = matmul(LAMBI(i,1:2,:),EV(6*i-2:6*i))
    end do

    !! Area coordinates of the stress samping points (edge mid-points)
    L1 = 0.5_dp; L1(2) = 0.0_dp
    L2 = 0.5_dp; L2(3) = 0.0_dp

    !! Loop over top and bottom surface
    ip = 0
    do k = 1, 2
       ZETA = real(3-2*k,dp)

       !! Loop over the stress sampling points
       do i = 1, 3
          call SCTS32 (SIG(1,1),LAMP(1,1),XG(1),YG(1),ZG(1),THK(1), &
               &       IEL,EMOD,RNY,L1(i),L2(i),ZETA, &
               &       LAMBI(1,1,1),ICODIR,LIN,IPSW,LPU,IERR)
          if (IERR < 0) then
             IERR = 1
             call reportError (warning_p,'Failure computing stress matrix')
             return
          end if

          !! Calculate local stress- and strain vectors at this point
          sigma(1:5,ip+3+i) = matmul(SIG,EV(1:nedof))
          epsil(1:5,ip+3+i) = matmul(Einv,sigma(1:5,ip+3+i))

          sigma(4:6,ip+3+i) = sigma(3:5,ip+3+i)
          sigma(  3,ip+3+i) = 0.0_dp ! Assume zero through-the-thickness stress
          epsil(4:6,ip+3+i) = epsil(3:5,ip+3+i)
          epsil(  3,ip+3+i) = 0.0_dp ! Assume zero through-the-thickness strain

          !! Transform to global coordinate system
          call tratensor (3,sigma(:,ip+3+i),LAMP)
          call tratensor (3,epsil(:,ip+3+i),LAMP)

       end do

       !! Extrapolate the corner nodes
       sigma(:,ip+1) = sigma(:,ip+4) + sigma(:,ip+6) - sigma(:,ip+5)
       sigma(:,ip+2) = sigma(:,ip+5) + sigma(:,ip+4) - sigma(:,ip+6)
       sigma(:,ip+3) = sigma(:,ip+6) + sigma(:,ip+5) - sigma(:,ip+4)
       epsil(:,ip+1) = epsil(:,ip+4) + epsil(:,ip+6) - epsil(:,ip+5)
       epsil(:,ip+2) = epsil(:,ip+5) + epsil(:,ip+4) - epsil(:,ip+6)
       epsil(:,ip+3) = epsil(:,ip+6) + epsil(:,ip+5) - epsil(:,ip+4)

       ip = ip + nenod
    end do

    !! Nodal stress resultants and conjugate strains ...  maybe later !
    SR = 0.0_dp
    SS = 0.0_dp

    ierr = 0

  end subroutine STR31


  !!****************************************************************************
  !> @brief Computes FE stresses and strains for a 8-noded shell element.

  subroutine STR32 (iel,EV,SR,SS,sigma,epsil,IPSW,LPU,IERR)

    use IsoMatModule                , only : isoMat2Dinv
    use FFaTensorTransformsInterface, only : tratensor
    use ReportErrorModule           , only : warning_p, reportError

    integer, parameter :: nenod = 8, nedof = 5*nenod

    integer , intent(in)    :: iel, IPSW, LPU
    real(dp), intent(inout) :: EV(6*nenod)      !< Element nodal displacements
    real(dp), intent(out)   :: SR(6,nenod)      !< Stress resultants at nodes
    real(dp), intent(out)   :: SS(6,nenod)      !< SR-conjugate strains
    real(dp), intent(out)   :: sigma(6,2*nenod) !< Stresses at result points
    real(dp), intent(out)   :: epsil(6,2*nenod) !< Strains at result points
    integer , intent(out)   :: IERR

    !! Local variables
    integer  :: i, j, k, ip, ICODIR
    real(dp) :: THK(nenod), XG(nenod), YG(nenod), ZG(nenod)
    real(dp) :: EMOD, RNY, XI, ETA, ZETA, XII(nenod), ETI(nenod)
    real(dp) :: LAMBI(nenod,3,3), LAMP(3,3), SIG(5,nedof)
    real(dp) :: SigPt(6,2,2), EpsPt(6,2,2), Einv(5,5)

    !! Factors to extrapolate from Gauss points to nodal values
    real(dp), parameter :: f1 = 0.5_dp + 0.5_dp*sqrt3_p
    real(dp), parameter :: f2 = 0.5_dp - 0.5_dp*sqrt3_p

    !! --- Logic section ---

    call getElementData ('STR32',IEL,XG,YG,ZG,EMOD,RNY,THK,ierr)
    if (ierr < 0) return

    !! Set up the inverse of the constitutive matrix defined in SCQS32
    call isoMat2Dinv (EMOD,RNY,Einv)
    Einv(4,4) = Einv(3,3)*1.2_dp
    Einv(5,5) = Einv(4,4)

    ICODIR = 1
    call SCQS30 (XG(1),YG(1),ZG(1),IEL,ICODIR, &
         &       XII(1),ETI(1),LAMBI(1,1,1),IPSW,LPU,IERR)
    if (IERR < 0) then
       IERR = 1
       call reportError (warning_p,'Failure computing nodal direction cosines')
       return
    end if

    !! Transform nodal rotations to local directions
    do i = 1, nenod
       EV(5*i-4:5*i-2) = EV(6*i-5:6*i-3)
       EV(5*i-1:5*i) = matmul(LAMBI(i,1:2,:),EV(6*i-2:6*i))
    end do

    !! Loop over top and bottom surface
    ip = 0
    do k = 1, 2
       ZETA = real(3-2*k,dp)

       !! Loop over the stress sampling points (Barlow points)
       do j = 1, 2
          ETA = real(2*j-3,dp)/sqrt3_p
          do i = 1, 2
             XI = real(2*i-3,dp)/sqrt3_p

             !! Evaluate the stress matrix
             call SCQS32 (SIG(1,1),LAMP(1,1),XG(1),YG(1),ZG(1),THK(1), &
                  &       IEL,EMOD,RNY,XI,ETA,ZETA, &
                  &       XII(1),ETI(1),LAMBI(1,1,1),ICODIR,IPSW,LPU,IERR)
             if (IERR < 0) then
                IERR = 1
                call reportError (warning_p,'Failure computing stress matrix')
                return
             end if

             !! Calculate local stress- and strain vectors at this point
             SigPt(1:5,i,j) = matmul(SIG,EV(1:nedof))
             EpsPt(1:5,i,j) = matmul(Einv,SigPt(1:5,i,j))

             SigPt(4:6,i,j) = SigPt(3:5,i,j)
             SigPt(  3,i,j) = 0.0_dp ! Assume zero through-the-thickness stress
             EpsPt(4:6,i,j) = EpsPt(3:5,i,j)
             EpsPt(  3,i,j) = 0.0_dp ! Assume zero through-the-thickness strain

             !! Transform to global coordinate system
             call tratensor (3,SigPt(:,i,j),LAMP)
             call tratensor (3,EpsPt(:,i,j),LAMP)

           end do
       end do

       !! Extrapolate the corner nodes
       sigma(:,ip+1) = f1*SigPt(:,1,1) + f2*SigPt(:,2,2)
       sigma(:,ip+3) = f1*SigPt(:,2,1) + f2*SigPt(:,1,2)
       sigma(:,ip+5) = f1*SigPt(:,2,2) + f2*SigPt(:,1,1)
       sigma(:,ip+7) = f1*SigPt(:,1,2) + f2*SigPt(:,2,1)
       epsil(:,ip+1) = f1*EpsPt(:,1,1) + f2*EpsPt(:,2,2)
       epsil(:,ip+3) = f1*EpsPt(:,2,1) + f2*EpsPt(:,1,2)
       epsil(:,ip+5) = f1*EpsPt(:,2,2) + f2*EpsPt(:,1,1)
       epsil(:,ip+7) = f1*EpsPt(:,1,2) + f2*EpsPt(:,2,1)
       !! Interpolate the mid-side nodes
       do i = 2, 8, 2
          sigma(:,ip+i) = 0.5_dp*(sigma(:,ip+i-1) + sigma(:,ip+mod(i,8)+1))
          epsil(:,ip+i) = 0.5_dp*(epsil(:,ip+i-1) + epsil(:,ip+mod(i,8)+1))
       end do

       ip = ip + nenod
    end do

    !! Nodal stress resultants and conjugate strains ...  maybe later !
    SR = 0.0_dp
    SS = 0.0_dp

    ierr = 0

  end subroutine STR32


  !!****************************************************************************
  !> @brief Computes FE stresses and strains for a 10-noded tetrahedron element.

  subroutine STR41 (iel, v, sigma, epsil, IPSW, LPU, IERR)

    use IsoMatModule          , only : isoMat3Dinv
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getint

    integer, parameter :: nenod = 10, nedof = 3*nenod, nstrp = 10

    integer , intent(in)  :: iel, IPSW, LPU
    real(dp), intent(in)  :: v(nedof)       !< Displacements at nodal points
    real(dp), intent(out) :: sigma(6,nenod) !< Stresses at nodal points
    real(dp), intent(out) :: epsil(6,nenod) !< Strains at nodal points
    integer , intent(out) :: IERR

    !! Local variables
    integer  :: n, stressForm
    real(dp) :: EMOD, RNY, Einv(6,6)
    real(dp) :: XG(nenod), YG(nenod), ZG(nenod), SIGG(6,nstrp)

    real(dp), parameter :: alpha_p = 1.927051062810166_dp
    real(dp), parameter :: beta_p = -0.309017015969668_dp

    !! --- Logic section ---

    call getElementData ('STR41',IEL,XG,YG,ZG,EMOD,RNY,IERR=ierr)
    if (ierr < 0) return
    !
    call isoMat3Dinv (EMOD,RNY,Einv)
    !
    call ffa_cmdlinearg_getint ('stressForm',stressForm)
    if (stressForm == 0) then
       n = 10 ! Direct nodal evaluation
    else
       n = 4 ! Extrapolate from 4 Gauss points
    end if
    !
    ! ------ ITET32 CALCULATES STRESSES IN THE SELECTED GAUSS-POINTS
    !        (NSTRP) OR DIRECTLY IN THE NODAL POINTS
    !        FOR A 10-NODED VOLUME ELEMENT (ISOPARAMETRIC TETRAHEDRON)
    !
    CALL ITET32 (SIGG(1,1),V(1),V(1),V(1),V(1),V(1), &
         &       XG(1),YG(1),ZG(1),EMOD,RNY,N,0,0,IEL,LPU,IPSW,IERR)
    if (ierr /= 0) then
       ierr = 1
       return
    end if
    !
    if (stressForm == 0) then
       sigma = SIGG
    else
       !! Extrapolate to the corner nodes (1,3,5,10)
       sigma(:,1)  = alpha_p*SIGG(:,1) + beta_p*(SIGG(:,2)+SIGG(:,3)+SIGG(:,4))
       sigma(:,3)  = alpha_p*SIGG(:,2) + beta_p*(SIGG(:,1)+SIGG(:,3)+SIGG(:,4))
       sigma(:,5)  = alpha_p*SIGG(:,3) + beta_p*(SIGG(:,1)+SIGG(:,2)+SIGG(:,4))
       sigma(:,10) = alpha_p*SIGG(:,4) + beta_p*(SIGG(:,1)+SIGG(:,2)+SIGG(:,3))
       !! Interpolate the mid-edge nodes (2,4,6,7,8,9)
       sigma(:,2)  = 0.5_dp*(sigma(:,1) + sigma(:,3))
       sigma(:,4)  = 0.5_dp*(sigma(:,3) + sigma(:,5))
       sigma(:,6)  = 0.5_dp*(sigma(:,5) + sigma(:,1))
       sigma(:,7)  = 0.5_dp*(sigma(:,1) + sigma(:,10))
       sigma(:,8)  = 0.5_dp*(sigma(:,3) + sigma(:,10))
       sigma(:,9)  = 0.5_dp*(sigma(:,5) + sigma(:,10))
    end if
    !
    ! --- COMPUTE NODAL STRAINS
    !
    epsil = matmul(Einv,sigma)

  end subroutine STR41


  !!****************************************************************************
  !> @brief Computes FE stresses and strains for a 15-noded wedge element.

  subroutine STR42 (iel, v, sigma, epsil, IPSW, LPU, IERR)

    use IsoMatModule          , only : isoMat3Dinv
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getint

    integer, parameter :: nenod = 15, nedof = 3*nenod, nstrp = 7, nstrpZ = 3

    integer , intent(in)  :: iel, IPSW, LPU
    real(dp), intent(in)  :: v(nedof)       !< Displacements at nodal points
    real(dp), intent(out) :: sigma(6,nenod) !< Stresses at nodal points
    real(dp), intent(out) :: epsil(6,nenod) !< Strains at nodal points
    integer , intent(out) :: IERR

    !! Local variables
    integer  :: n, nZ, stressForm
    real(dp) :: EMOD, RNY, Einv(6,6), XG(nenod), YG(nenod), ZG(nenod)
    real(dp) :: SIGG(6,nstrp*nstrpZ)

    real(dp), parameter :: zm1_p = 0.5_dp*sqrt3_p - 0.5_dp
    real(dp), parameter :: zp1_p = zm1_p + 1.0_dp

    !! --- Logic section ---

    call getElementData ('STR42',IEL,XG,YG,ZG,EMOD,RNY,IERR=ierr)
    if (ierr < 0) return
    !
    call isoMat3Dinv (EMOD,RNY,Einv)
    !
    call ffa_cmdlinearg_getint ('stressForm',stressForm)
    if (stressForm == 0) then
       n  = 6 ! Direct nodal evaluation
       nZ = 3
    else
       n  = 3 ! Extrapolate from 3x2 Gauss points giving a linear variation
       nZ = 2
    end if
    !
    ! ------ IPRI32 CALCULATES STRESSES IN THE SELECTED GAUSS-POINTS
    !        (NSTRP*NSTRPZ) OR DIRECTLY IN THE NODAL POINTS
    !        FOR A 15-NODED VOLUME ELEMENT (ISOPARAMETRIC PRISM).
    !
    CALL IPRI32 (SIGG(1,1),V(1),V(1),V(1),V(1),V(1), &
         &       XG(1),YG(1),ZG(1),EMOD,RNY,n,nZ,0,0,IEL,LPU,IPSW,IERR)
    if (ierr /= 0) then
       ierr = 1
       return
    end if
    !
    if (stressForm == 0) then
       !! Extract the nodal stresses
       sigma(:,1:7)   = SIGG(:,1:7)
       sigma(:,8)     = SIGG(:,9)
       sigma(:,9)     = SIGG(:,11)
       sigma(:,10:15) = SIGG(:,13:18)
    else
       !! Extrapolate to the corner edges of the prism
       epsil(:,1)  = SIGG(:,3) + SIGG(:,1) - SIGG(:,2)
       epsil(:,2)  = SIGG(:,2) + SIGG(:,1) - SIGG(:,3)
       epsil(:,3)  = SIGG(:,2) + SIGG(:,3) - SIGG(:,1)
       epsil(:,4)  = SIGG(:,6) + SIGG(:,4) - SIGG(:,5)
       epsil(:,5)  = SIGG(:,5) + SIGG(:,4) - SIGG(:,6)
       epsil(:,6)  = SIGG(:,5) + SIGG(:,6) - SIGG(:,4)
       !! Extrapolate to the top and bottom triangular faces
       sigma(:,1)  = zp1_p*epsil(:,1) - zm1_p*epsil(:,4)
       sigma(:,3)  = zp1_p*epsil(:,2) - zm1_p*epsil(:,5)
       sigma(:,5)  = zp1_p*epsil(:,3) - zm1_p*epsil(:,6)
       sigma(:,10) = zp1_p*epsil(:,4) - zm1_p*epsil(:,1)
       sigma(:,12) = zp1_p*epsil(:,5) - zm1_p*epsil(:,2)
       sigma(:,14) = zp1_p*epsil(:,6) - zm1_p*epsil(:,3)
       !! Interpolate the mid-edge nodes
       sigma(:,2)  = 0.5_dp*(sigma(:,1)  + sigma(:,3))
       sigma(:,4)  = 0.5_dp*(sigma(:,3)  + sigma(:,5))
       sigma(:,6)  = 0.5_dp*(sigma(:,5)  + sigma(:,1))
       sigma(:,7)  = 0.5_dp*(sigma(:,1)  + sigma(:,10))
       sigma(:,8)  = 0.5_dp*(sigma(:,3)  + sigma(:,12))
       sigma(:,9)  = 0.5_dp*(sigma(:,5)  + sigma(:,14))
       sigma(:,11) = 0.5_dp*(sigma(:,10) + sigma(:,12))
       sigma(:,13) = 0.5_dp*(sigma(:,12) + sigma(:,14))
       sigma(:,15) = 0.5_dp*(sigma(:,14) + sigma(:,10))
    end if
    !
    ! --- COMPUTE NODAL STRAINS
    !
    epsil = matmul(Einv,sigma)

  end subroutine STR42


  !!****************************************************************************
  !> @brief Computes FE stresses and strains for a 20-noded hexahedron element.

  subroutine STR43 (iel, v, sigma, epsil, IPSW, LPU, IERR)

    use IsoMatModule          , only : isoMat3Dinv
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getint

    integer, parameter :: nenod = 20, nedof = 3*nenod, nstrp = 27

    integer , intent(in)  :: iel, IPSW, LPU
    real(dp), intent(in)  :: v(nedof)       !< Displacements at nodal points
    real(dp), intent(out) :: sigma(6,nenod) !< Stresses at nodal points
    real(dp), intent(out) :: epsil(6,nenod) !< Strains at nodal points
    integer , intent(out) :: IERR

    !! Local variables
    integer  :: i, j, k, l, n, stressForm
    real(dp) :: EMOD, RNY, Einv(6,6), x, y, z
    real(dp) :: XG(nenod), YG(nenod), ZG(nenod), SIGG(6,nstrp)

    real(dp), parameter :: zero_p = 0.0_dp
    real(dp), parameter :: one_p = sqrt3_p
    real(dp), parameter :: Xinod_p(60) = (/-one_p,-one_p,-one_p, &
         &                                 zero_p,-one_p,-one_p, &
         &                                  one_p,-one_p,-one_p, &
         &                                  one_p,zero_p,-one_p, &
         &                                  one_p, one_p,-one_p, &
         &                                 zero_p, one_p,-one_p, &
         &                                 -one_p, one_p,-one_p, &
         &                                 -one_p,zero_p,-one_p, &
         &                                 -one_p,-one_p,zero_p, &
         &                                  one_p,-one_p,zero_p, &
         &                                  one_p, one_p,zero_p, &
         &                                 -one_p, one_p,zero_p, &
         &                                 -one_p,-one_p, one_p, &
         &                                 zero_p,-one_p, one_p, &
         &                                  one_p,-one_p, one_p, &
         &                                  one_p,zero_p, one_p, &
         &                                  one_p, one_p, one_p, &
         &                                 zero_p, one_p, one_p, &
         &                                 -one_p, one_p, one_p, &
         &                                 -one_p,zero_p, one_p /)
    real(dp), parameter :: Xin_p(3,20) = reshape(Xinod_p,(/3,20/))

    !! --- Logic section ---

    call getElementData ('STR43',IEL,XG,YG,ZG,EMOD,RNY,IERR=ierr)
    if (ierr < 0) return
    !
    call isoMat3Dinv (EMOD,RNY,Einv)
    !
    call ffa_cmdlinearg_getint ('stressForm',stressForm)
    if (stressForm == 0) then
       n = 0 ! Direct nodal evaluation
    else
       n = 2 ! Extrapolate from 2x2x2 Gauss points
    end if
    !
    ! ------ IHEX32 CALCULATES STRESSES IN THE SELECTED GAUSS-POINTS
    !        (NSTRXI*NSTRET*NSTRZE) OR DIRECTLY IN THE NODAL POINTS
    !        FOR A 20-NODED VOLUME ELEMENT (ISOPARAMETRIC HEXAHEDRON).
    !
    CALL IHEX32 (SIGG(1,1),V(1),V(1),V(1),V(1),V(1), &
         &       XG(1),YG(1),ZG(1),EMOD,RNY,n,n,n,0,0,IEL,LPU,IPSW,IERR)
    if (ierr /= 0) then
       ierr = 1
       return
    end if
    !
    if (stressForm == 0) then
       !! Extract the nodal stresses
       sigma(:,1)  = SIGG(:,1)
       sigma(:,2)  = SIGG(:,2)
       sigma(:,3)  = SIGG(:,3)
       sigma(:,4)  = SIGG(:,6)
       sigma(:,5)  = SIGG(:,9)
       sigma(:,6)  = SIGG(:,8)
       sigma(:,7)  = SIGG(:,7)
       sigma(:,8)  = SIGG(:,4)
       sigma(:,9)  = SIGG(:,10)
       sigma(:,10) = SIGG(:,12)
       sigma(:,11) = SIGG(:,18)
       sigma(:,12) = SIGG(:,16)
       sigma(:,13) = SIGG(:,19)
       sigma(:,14) = SIGG(:,20)
       sigma(:,15) = SIGG(:,21)
       sigma(:,16) = SIGG(:,24)
       sigma(:,17) = SIGG(:,27)
       sigma(:,18) = SIGG(:,26)
       sigma(:,19) = SIGG(:,25)
       sigma(:,20) = SIGG(:,22)
    else
       !! Extrapolate nodal stresses from the 2x2x2 Gauss points
       sigma = 0.0_dp
       do n = 1, nenod
          l = 0
          do k = -1, 1, 2
             z = 1.0_dp + k*Xin_p(3,n)
             do j = -1, 1, 2
                y = 1.0_dp + j*Xin_p(2,n)
                do i = -1, 1, 2
                   x = 1.0_dp + i*Xin_p(1,n)
                   l = l + 1
                   sigma(:,n) = sigma(:,n) + SIGG(:,l)*x*y*z*0.125_dp
                end do
             end do
          end do
       end do
    end if
    !
    ! --- COMPUTE NODAL STRAINS
    !
    epsil = matmul(Einv,sigma)

  end subroutine STR43


  !!****************************************************************************
  !> @brief Computes FE stresses and strains for a 8-noded hexahedron element.

  subroutine STR44 (iel, v, sigma, epsil, IPSW, LPU, IERR)

    use IsoMatModule          , only : isoMat3D, isoMat3Dinv
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getint
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_isTrue

    integer, parameter :: NENOD = 8, NEDOF = 3*NENOD, NIDOF = 9, NIP = 3

    integer , intent(in)  :: iel, IPSW, LPU
    real(dp), intent(in)  :: v(nedof)       !< Displacements at nodal points
    real(dp), intent(out) :: sigma(6,nenod) !< Stresses at nodal points
    real(dp), intent(out) :: epsil(6,nenod) !< Strains at nodal points
    integer , intent(out) :: IERR

    !! Local variables
    integer  :: i, n, iop, lop, stressForm
    real(dp) :: EMOD, RNY, E(36), Einv(6,6), detJ, vol, x, y, z, sig(6)
    real(dp) :: XG(nenod), YG(nenod), ZG(nenod), Vi(nidof)
    real(dp) :: EK(nedof,nedof), A1(nedof,nidof), A2(nidof,nidof)
    real(dp) :: B(6,nidof), SI(6,nedof)

    real(dp), parameter :: XI(8)   = (/-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0,-1.0/)
    real(dp), parameter :: ETA(8)  = (/-1.0,-1.0, 1.0, 1.0,-1.0,-1.0, 1.0, 1.0/)
    real(dp), parameter :: ZETA(8) = (/-1.0,-1.0,-1.0,-1.0, 1.0, 1.0, 1.0, 1.0/)

    !! --- Logic section ---

    call getElementData ('STR44',IEL,XG,YG,ZG,EMOD,RNY,IERR=ierr)
    if (ierr < 0) return
    !
    call isoMat3D (EMOD,RNY,E)
    call isoMat3Dinv (EMOD,RNY,Einv)
    !
    if (ffa_cmdlinearg_isTrue('useIncompatibleModes')) then
       IOP = 1 ! Internal DOFs are included ==> HEX11 element
       CALL HEXA31 (EK(1,1),A1(1,1),A2(1,1),E(1),XG(1),YG(1),ZG(1), &
            &       NIP,IOP,IEL,LPU,IPSW,IERR)
       IF (IERR /= 0) THEN
          ierr = 1
          return
       END IF
       vi = matmul(v,A1) ! Extract the internal DOF values
    else
       IOP = 0 ! No internal DOFs are included ==> HEX8 element
    end if

    !!============================================
    !! Stress formulation option for HEX8/HEX11:
    !! = 0 : Direct nodal evaluation
    !! = 1 : Volume-averaged stress (constant)
    !! > 1 : Extrapolation from 2x2x2 Gauss points
    !!============================================
    call ffa_cmdlinearg_getint ('stressForm',stressForm)
    if (stressForm == 0) then
       LOP = 0
    else
       LOP = NIP-1
       sig = 0.0_dp
       vol = 0.0_dp
    end if

    do n = 1, nenod

       !! Evaluate the stress-displacement matrices, SI and B
       CALL HEXA32 (SI(1,1),B(1,1),A1(1,1),A2(1,1),detJ, &
            &       E(1),XG(1),YG(1),ZG(1),XI(n),ETA(n),ZETA(n),IOP,LOP,IERR)
       if (ierr /= 0) then
          ierr = 1
          return
       end if

       sigma(:,n) = matmul(SI,v)
       if (iop == 1) then
          !! Add contributions from the internal DOFs
          sigma(:,n) = sigma(:,n) + matmul(B,vi)
       end if

       if (stressForm == 0) then
          !! Evaluate nodal strains
          epsil(:,n) = matmul(Einv,sigma(:,n))
       else if (stressForm == 1) then
          !! Accumulate the volume-averaged stress
          sig = sig + sigma(:,n)*detJ
          vol = vol + detJ
       end if

    end do

    if (stressForm == 1) then
       !! Compute volume-averaged stresses and assign to each element node
       sigma(:,1) = sig/vol
       epsil(:,1) = matmul(Einv,sigma(:,1))
       do n = 2, nenod
          sigma(:,n) = sigma(:,1)
          epsil(:,n) = epsil(:,1)
       end do
    else if (stressForm >= 2) then
       !! Extrapolate nodal stresses from the 2x2x2 Gauss points
       si = 0.0_dp
       do n = 1, nenod
          x = -XI(n)*sqrt3_p
          y = -ETA(n)*sqrt3_p
          z = -ZETA(n)*sqrt3_p
          do i = 1, 8
             call LINHEX (i,x,y,z,rny,0)
             si(:,n) = si(:,n) + sigma(:,i)*rny
          end do
       end do
       sigma = si(:,1:nenod)
       epsil = matmul(Einv,sigma)
    end if

  end subroutine STR44


  !!****************************************************************************
  !> @brief Computes FE stresses and strains for a 4-noded tetrahedron element.

  subroutine STR45 (iel, v, sigma, epsil, IPSW, LPU, IERR)

    use IsoMatModule , only : isoMat3D
    use CSTetraModule, only : CSTetStrain

    integer, parameter :: nenod = 4, nedof = 3*nenod

    integer , intent(in)  :: iel, IPSW, LPU
    real(dp), intent(in)  :: v(nedof)       !< Displacements at nodal points
    real(dp), intent(out) :: sigma(6,nenod) !< Stresses at nodal points
    real(dp), intent(out) :: epsil(6,nenod) !< Strains at nodal points
    integer , intent(out) :: IERR

    !! Local variables
    integer  :: n
    real(dp) :: EMOD, RNY, E(6,6), XG(nenod), YG(nenod), ZG(nenod)

    !! --- Logic section ---

    call getElementData ('STR45',IEL,XG,YG,ZG,EMOD,RNY,IERR=ierr)
    if (ierr < 0) return

    call isoMat3D (EMOD,RNY,E)

    call CSTetStrain (XG,YG,ZG,V,epsil(:,1),RNY,LPU,IPSW,IERR)
    if (ierr /= 0) then
       ierr = 1
       return
    end if

    sigma(:,1) = matmul(E,epsil(:,1))

    do n = 2, nenod
       sigma(:,n) = sigma(:,1)
       epsil(:,n) = epsil(:,1)
    end do

  end subroutine STR45


  !!****************************************************************************
  !> @brief Computes FE stresses and strains for a 6-noded wedge element.

  subroutine STR46 (iel, v, sigma, epsil, IPSW, LPU, IERR)

    use IsoMatModule          , only : isoMat3D
    use Ipri6Module           , only : Ipri6Strain
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getint

    integer, parameter :: nenod = 6, nedof = 3*nenod

    integer , intent(in)  :: iel, IPSW, LPU
    real(dp), intent(in)  :: v(nedof)       !< Displacements at nodal points
    real(dp), intent(out) :: sigma(6,nenod) !< Stresses at nodal points
    real(dp), intent(out) :: epsil(6,nenod) !< Strains at nodal points
    integer , intent(out) :: IERR

    !! Local variables
    integer  :: stressForm
    real(dp) :: EMOD, RNY, E(6,6), XG(nenod), YG(nenod), ZG(nenod)

    !! --- Logic section ---

    call getElementData ('STR46',IEL,XG,YG,ZG,EMOD,RNY,IERR=ierr)
    if (ierr < 0) return

    call isoMat3D (EMOD,RNY,E)

    !! Legal choices: 0,1,2,3 (see subroutine Ipri6Strain for details)
    call ffa_cmdlinearg_getint ('stressForm',stressForm)
    if (stressForm == 1) stressForm = 3

    call Ipri6Strain (XG,YG,ZG,V,epsil,IEL,stressForm,LPU,IPSW,IERR)
    if (ierr < 0) then
       ierr = 1
    else if (ierr == 0) then
       sigma = matmul(E,epsil)
    end if

  end subroutine STR46

  !> @endcond
end module ElstressModule
