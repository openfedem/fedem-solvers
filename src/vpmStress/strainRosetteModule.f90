!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module StrainRosetteModule

  use KindModule  , only : dp, sp
  use IdTypeModule, only : IdType

  implicit none

#if FT_HAS_RECOVERY == 1
  integer, parameter :: rk = sp !< Single precision real kind
#else
  integer, parameter :: rk = dp !< Double precision real kind
#endif

  type StrainGageType
     integer  :: lpu_FiDF       !< Device function file for writing
     integer  :: fatHandle      !< Handle to fatigue calculation data
     real(dp) :: Teps_NfromC(3) !< Strain gage transformation
     real(dp) :: epsGage        !< Current gage strain (time-dependent)
     real(dp) :: sigGage        !< Current gage stress (time-dependent)
  end type StrainGageType


  type StrainRosetteType

     !! Constant (time-independent) quantities
     real(rk)            , pointer :: Bcart(:,:) !< Strain-displacement matrix
     type(StrainGageType), pointer :: gages(:)   !< Data for all gages

     integer  :: fatHandle   !< Handle to fatigue calculation data
     real(dp) :: gateValue   !< Gate value [MPa] for PVX and rainflow analysis
     real(dp) :: snCurve(4)  !< Parameters defining the S-N curve to use

     real(dp) :: alphaGages  !< Angle between the gages
     real(dp) :: zPos        !< Position along the local Z-axis
     real(dp) :: Cmat(3,3)   !< Constitutive matrix
     real(dp) :: epsCInit(3) !< Initial Cartesian strains
     real(dp) :: sigmaC0(3)  !< Initial residual stresses
     logical  :: zeroInit    !< Make gage strains of the first time step zero?

     !! Time dependent quantities
     real(dp) :: epsC(3)     !< Cartesian strains
     real(dp) :: epsP(3)     !< Principal strain (max, min, signed abs max)
     real(dp) :: gammaMax    !< Maximum shear strain
     real(dp) :: epsVM       !< von Mises strain
     real(dp) :: sigmaC(3)   !< Cartesian stresses
     real(dp) :: sigmaP(3)   !< Principal stresses (max, min, signed abs max)
     real(dp) :: tauMax      !< Maximum shear stress
     real(dp) :: sigmaVM     !< von Mises stress
     real(dp) :: alpha1      !< Angle of max principle strain/stress
     real(dp) :: alphaGamma  !< Angle of max shear

  end type StrainRosetteType


  type StrainElementType

     type(IdType) :: id    !< General identification data

     integer :: elmNumber  !< Which element this strain rosette is defined on
     integer :: linkNumber !< Which FE part this strain rosette is attached to
     integer :: lpuASCII   !< File for ASCII print

     integer , pointer :: globalNodes(:) !< Nodes of strain-measuring element
     real(dp)          :: posInGl(3,4)   !< Global rosette position matrix
     real(dp), pointer :: X0(:,:)        !< Nodal coordinates
     real(dp), pointer :: T0(:,:)        !< Initial element coordinate system
     real(rk), pointer :: disp(:,:)      !< Nodal deformational displacements
     real(rk), pointer :: vgii(:,:)      !< Nodal gravitational displacements
     real(dp), pointer :: ur(:)          !< Global displacements and rotations

     type(StrainRosetteType), pointer :: data(:) !< Surface strain rosettes

  end type StrainElementType


  interface nullifyRosette
     module procedure nullifyStrainRosette
     module procedure nullifyStrainElement
  end interface

  private :: nullifyStrainRosette, nullifyStrainElement


contains

  subroutine nullifyStrainRosette (rosette)

    !!==========================================================================
    !! Initialize the StrainRosetteType object.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 30 May 2001/1.0
    !!==========================================================================

    type(StrainRosetteType), intent(out) :: rosette

    !! --- Logic section ---

    nullify(rosette%Bcart)
    nullify(rosette%gages)

    rosette%fatHandle  = 0
    rosette%gateValue  = 0.0_dp
    rosette%snCurve    = 0.0_dp
    rosette%alphaGages = 0.0_dp
    rosette%zPos       = 0.0_dp
    rosette%Cmat       = 0.0_dp
    rosette%epsCInit   = 0.0_dp
    rosette%sigmaC0    = 0.0_dp
    rosette%zeroInit   = .false.
    rosette%epsC       = 0.0_dp
    rosette%epsP       = 0.0_dp
    rosette%gammaMax   = 0.0_dp
    rosette%epsVM      = 0.0_dp
    rosette%sigmaC     = 0.0_dp
    rosette%sigmaP     = 0.0_dp
    rosette%tauMax     = 0.0_dp
    rosette%sigmaVM    = 0.0_dp
    rosette%alpha1     = 0.0_dp
    rosette%alphaGamma = 0.0_dp

  end subroutine nullifyStrainRosette


  subroutine nullifyStrainElement (rosette)

    !!==========================================================================
    !! Initialize the StrainElementType object.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 6 Feb 2004/1.0
    !!==========================================================================

    use IdTypeModule, only : nullifyId

    type(StrainElementType), intent(out) :: rosette

    !! --- Logic section ---

    call nullifyId (rosette%id)

    rosette%elmNumber  = 0
    rosette%linkNumber = 0
    rosette%lpuASCII   = 0
    rosette%posInGl    = 0.0_dp
    nullify(rosette%globalNodes)
    nullify(rosette%X0)
    nullify(rosette%T0)
    nullify(rosette%disp)
    nullify(rosette%vgii)
    nullify(rosette%ur)
    nullify(rosette%data)

  end subroutine nullifyStrainElement


  subroutine deallocateRosette (rosette)

    !!==========================================================================
    !! Deallocate the StrainElementType object.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 23 Jan 2017/1.0
    !!==========================================================================

    type(StrainElementType), intent(inout) :: rosette

    !! Local variables
    integer :: i

    !! --- Logic section ---

    if (associated(rosette%id%assId)) then
       deallocate(rosette%id%assId)
    end if

    if (associated(rosette%globalNodes)) then
       deallocate(rosette%globalNodes)
    end if

    if (associated(rosette%X0)) then
       deallocate(rosette%X0)
    end if

    if (associated(rosette%T0)) then
       deallocate(rosette%T0)
    end if

    if (associated(rosette%disp)) then
       deallocate(rosette%disp)
    end if

    if (associated(rosette%vgii)) then
       deallocate(rosette%vgii)
    end if

    if (associated(rosette%ur)) then
       deallocate(rosette%ur)
    end if

    if (associated(rosette%data)) then
       do i = 1, size(rosette%data)
          if (associated(rosette%data(i)%Bcart)) then
             deallocate(rosette%data(i)%Bcart)
          end if
          if (associated(rosette%data(i)%gages)) then
             deallocate(rosette%data(i)%gages)
          end if
       end do
       deallocate(rosette%data)
    end if

    call nullifyStrainElement (rosette)

  end subroutine deallocateRosette


  subroutine evaluateStrainGages (gages,epsC,sigC)

    !!==========================================================================
    !! Evaluate the given strain gages based on the cartesian strain state epsC
    !! and the cartesian stress state sigC.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 30 May 2001/1.0
    !!==========================================================================

    type(StrainGageType), intent(inout) :: gages(:)
    real(dp)            , intent(in)    :: epsC(:), sigC(:)

    !! Local variables
    integer :: i

    !! --- Logic section ---

    do i = 1, size(gages)
       gages(i)%epsGage = dot_product(gages(i)%Teps_NfromC,epsC)
       gages(i)%sigGage = dot_product(gages(i)%Teps_NfromC,sigC)
    end do

  end subroutine evaluateStrainGages


  subroutine calcRosetteStrains (rosette,displ,ierr)

    !!==========================================================================
    !! Evaluate the given strain rosette based on the given displacement vector.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 30 May 2001/2.0
    !!==========================================================================

    use StrainAndStressUtilitiesModule, only : PrincipleStrains2D
    use StrainAndStressUtilitiesModule, only : PrincipleStresses2D
    use reportErrorModule             , only : internalError

    type(StrainRosetteType), intent(inout) :: rosette
    real(rk)               , intent(in)    :: displ(:)
    integer                , intent(out)   :: ierr

    !! --- Logic section ---

    if (.not. associated(rosette%Bcart)) then
       ierr = internalError('calcRosetteStrains: No B-matrix for rosette')
       return
    else if (size(displ) /= size(rosette%Bcart,2)) then
       ierr = internalError('calcRosetteStrains: Invalid displacement vector')
       return
    else
       ierr = 0
    end if

    !! Compute the Cartesian strain and stress
    rosette%epsC   = matmul(rosette%Bcart,displ) + rosette%epsCInit
    rosette%sigmaC = matmul(rosette%Cmat,rosette%epsC) + rosette%sigmaC0

    if (associated(rosette%gages)) then
       !! Compute the strain-gage strains
       call evaluateStrainGages (rosette%gages,rosette%epsC,rosette%sigmaC)
    end if

    !! Compute principle strain, max shear strain, angle of principle strain
    call PrincipleStrains2D (rosette%epsC, &
         &                   rosette%epsP(1), rosette%epsP(2), &
         &                   rosette%gammaMax, &
         &                   rosette%alpha1, rosette%alphaGamma)

    !! Signed abs max strain
    if (abs(rosette%epsP(1)) > abs(rosette%epsP(2))) then
       rosette%epsP(3) = rosette%epsP(1)
    else
       rosette%epsP(3) = rosette%epsP(2)
    end if

    !! von Mises strain
    rosette%epsVM = sqrt(rosette%epsP(1)*rosette%epsP(1) &
         &             + rosette%epsP(2)*rosette%epsP(2) &
         &             - rosette%epsP(1)*rosette%epsP(2))

    !! Compute principle stress and max shear stress
    call PrincipleStresses2D (rosette%sigmaC, &
         &                    rosette%sigmaP(1), rosette%sigmaP(2), &
         &                    rosette%tauMax)

    !! Signed abs max stress
    if (abs(rosette%sigmaP(1)) > abs(rosette%sigmaP(2))) then
       rosette%sigmaP(3) = rosette%sigmaP(1)
    else
       rosette%sigmaP(3) = rosette%sigmaP(2)
    end if

    !! von Mises stress
    rosette%sigmaVM = sqrt(rosette%sigmaP(1)*rosette%sigmaP(1) &
         &               + rosette%sigmaP(2)*rosette%sigmaP(2) &
         &               - rosette%sigmaP(1)*rosette%sigmaP(2))

  end subroutine calcRosetteStrains


  subroutine calcZeroStartRosetteStrains (rosette,displ,ierr)

    !!==========================================================================
    !! Evaluate the initial strain of a strain rosette.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 30 May 2001/2.0
    !!==========================================================================

    use reportErrorModule, only : internalError

    type(StrainRosetteType), intent(inout) :: rosette
    real(rk)               , intent(in)    :: displ(:)
    integer                , intent(out)   :: ierr

    !! --- Logic section ---

    if (.not. associated(rosette%Bcart)) then
       ierr = internalError('calcZeroStartRosetteStrains: No B-matrix')
    else if (size(displ) /= size(rosette%Bcart,2)) then
       ierr = internalError('calcZeroStartRosetteStrains: Invalid displacement')
    else
       ierr = 0
       rosette%epsCInit = -matmul(rosette%Bcart,displ)
       rosette%zeroInit = .false. ! Compute initial strains only once
    end if

  end subroutine CalcZeroStartRosetteStrains


  subroutine PrintRosetteStrains (rosette,time,lpu)

    !!==========================================================================
    !! Print the current state of a strain rosette to ASCII file.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 25 May 2001 / 2.0
    !!==========================================================================

    use kindModule, only : pi_p

    type(StrainRosetteType), intent(in) :: rosette
    real(dp)               , intent(in) :: time
    integer                , intent(in) :: lpu

    !! Local variables
    integer  :: i, nGages
    real(dp) :: alpha1, alphaGamma

    if (lpu <= 0) return

    !! --- Logic section ---

    if (associated(rosette%gages)) then
       nGages = size(rosette%gages)
    else
       nGages = 0
    end if

    !! Note: only multiply with 18 (and not 180) due to the 1P format
    alpha1     = rosette%alpha1     * 18.0_dp/pi_p
    alphaGamma = rosette%alphaGamma * 18.0_dp/pi_p

    write(lpu,"(1P,E14.6E3,6E11.3E2,2F11.5,12E11.3E2)") time, rosette%epsC, &
         rosette%epsP(1:2), rosette%gammaMax, alpha1, alphaGamma, &
         rosette%sigmaC, rosette%sigmaVM, rosette%sigmaP(1:2), rosette%tauMax, &
         (rosette%gages(i)%epsGage, i=1,nGages)

  end subroutine PrintRosetteStrains


  subroutine PrintInitialStrain (rosette,time,lpu)

    !!==========================================================================
    !! Print the initial state of a strain rosette to ASCII file.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 25 May 2001 / 2.0
    !!==========================================================================

    use StrainAndStressUtilitiesModule, only : PrincipleStrains2D
    use StrainAndStressUtilitiesModule, only : PrincipleStresses2D

    type(StrainRosetteType), intent(inout) :: rosette
    real(dp)               , intent(in)    :: time
    integer                , intent(in)    :: lpu

    !! --- Logic section ---

    if (lpu <= 0) return

    !! Compute the Cartesian strain and stress
    rosette%epsC   = -rosette%epsCInit
    rosette%sigmaC = matmul(rosette%Cmat,rosette%epsC) + rosette%sigmaC0

    if (associated(rosette%gages)) then
       !! Compute the strain-gage strains
       call evaluateStrainGages (rosette%gages,rosette%epsC,rosette%sigmaC)
    end if

    !! Compute principle strain, max shear strain, angle of principle strain
    call PrincipleStrains2D (rosette%epsC, &
         &                   rosette%epsP(1), rosette%epsP(2), &
         &                   rosette%gammaMax, &
         &                   rosette%alpha1, rosette%alphaGamma)

    !! Compute principle stress and max shear stress
    call PrincipleStresses2D (rosette%sigmaC, &
         &                    rosette%sigmaP(1), rosette%sigmaP(2), &
         &                    rosette%tauMax)

    write(lpu,"('# Initial strain:')")
    call PrintRosetteStrains (rosette,time,lpu)
    write(lpu,"('# Strain readings (= True strain - Initial strain):')")

  end subroutine PrintInitialStrain


  subroutine PrintRosetteHeading (rosette, globalNodes, minex, posInGl, &
       &                          rosId, feId, lpu)

    !!==========================================================================
    !! Print ASCII file heading for a given strain rosette.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 25 May 2001 / 2.0
    !!==========================================================================

    use kindModule, only : pi_p

    type(StrainRosetteType), intent(in) :: rosette
    integer                , intent(in) :: globalNodes(:), minex(:)
    real(dp)               , intent(in) :: posInGl(:,:)
    integer                , intent(in) :: rosId, feId, lpu

    !! Local variables
    integer  :: i, nGages
    real(dp) :: alpha

    !! --- Logic section ---

    if (lpu <= 0) return

    if (associated(rosette%gages)) then
       nGages = size(rosette%gages)
    else
       nGages = 0
    end if

    !! Write information about the strain gage
    write(lpu,"('#')")
    write(lpu,"('#  Strain rosette identifier:',I8,'   Part :',I8)") rosId, feId
    write(lpu,"('#  Global nodes             :',9I8)") minex(globalNodes)
    write(lpu,"('#  Number of gages          :',I8)")  nGages
    if (nGages > 1) then
       alpha = rosette%alphaGages * 180.0_dp/pi_p
       write(lpu,"('#  Angle between gages      :',F12.3)") alpha
    end if

    write(lpu,"('#  Position                 :',1P,3E12.3E3)") posInGl(:,4)
    write(lpu,"('#  Position along Z-axis    :',1P, E12.3E3)") rosette%zPos
    write(lpu,"('#  X-direction unit vector  :',1P,3E12.3E3)") posInGl(:,1)
    write(lpu,"('#  Y-direction unit vector  :',1P,3E12.3E3)") posInGl(:,2)
    write(lpu,"('#  Z-direction unit vector  :',1P,3E12.3E3)") posInGl(:,3)

    !! Write the heading for the results
    write(lpu,"('#')")
    write(lpu,ADVANCE='YES',FMT="('#',I9,20I11)") (i,i=1,16+nGages)
    write(lpu,ADVANCE='NO', FMT="('#',A13,20A11)") 'time','eps_x','eps_y', &
         'gamma_xy','eps1','eps2','gammaMax','alpha1','alphaGamm', &
         'sigma_x','sigma_y','tau_xy','von_Mises','sig1','sig2','tau_max'
    do i = 1, nGages
       write(lpu,ADVANCE='NO',FMT="(A10,I1)") 'eps_gage',i
    end do
    write(lpu,*)

  end subroutine PrintRosetteHeading


  subroutine calcElmCoordSystem (id, globalNodes, posInGl, X_el, T_el, &
       &                         lpu, ierr, useElCoordSys)

    !!==========================================================================
    !! Calulate nodal coordinates and transformation matrix for the
    !! strain-gage element. Reorder the nodes and repeat calculation if the
    !! initial ordering gives the wrong Z-axis direction.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 25 May 2001 / 1.0
    !!==========================================================================

    use kindModule                    , only : epsDiv0_p
    use IdTypeModule                  , only : getId
    use StrainAndStressUtilitiesModule, only : getShellElementAxes
    use rotationModule                , only : orthoNorm3
    use manipMatrixModule             , only : cross_product
    use reportErrorModule             , only : reportError, error_p
    use FFlLinkHandlerInterface       , only : ffl_getNodalCoor

    type(IdType)     , intent(in)    :: id
    integer          , intent(in)    :: globalNodes(:), lpu
    real(dp)         , intent(inout) :: posInGl(:,:)
    real(dp)         , intent(out)   :: X_el(:,:), T_el(:,:)
    integer          , intent(out)   :: ierr
    logical, optional, intent(in)    :: useElCoordSys

    !! Local variables
    integer :: i, nElNodes

    !! --- Logic section ---

    ierr = 0
    nElNodes = size(globalNodes)

    !! Get global coordinates of the element nodes
    do i = 1, nElNodes
       call ffl_getNodalCoor (X_el(1,i),X_el(2,i),X_el(3,i),globalNodes(i),ierr)
       if (ierr /= 0) then
          call reportError (error_p,'Unable to get element coordinates for '// &
               &            'Rosette'//getId(id),addString='calcElmCoordSystem')
          return
       end if
    end do

    !! Get the transformation matrix from global to element coordinate system
    call getShellElementAxes (nElNodes, X_el(1,:), X_el(2,:), X_el(3,:), &
         &                    T_el(1,:), T_el(2,:), T_el(3,:), &
         &                    lpu, ierr, useElCoordSys)
    if (ierr /= 0) goto 910
    if (.not.present(useElCoordSys)) return ! Do not recompute matrix posInGl

    !! Calculate the strain gage position
    do i = 1, 3
       posInGl(i,4) = sum(X_el(i,:)) / nElNodes
    end do

    if (useElCoordSys) then
       !! Use element coordinate system the for X- and Z-direction definitions
       posInGl(:,1:3) = transpose(T_el)
       return
    end if

    !! Check that the element Z-axis points in the correct direction
    if (dot_product(T_el(3,:),posInGl(:,3)) < epsDiv0_p) goto 900

    !! Finish initializing the orientation part of the rosettes location
    posInGl(:,3) = T_el(3,:) ! Same Z-axis as the element
    posInGl(:,2) = cross_product(posInGl(:,3),posInGl(:,1))
    posInGl(:,1) = cross_product(posInGl(:,2),posInGl(:,3))
    call orthoNorm3 (posInGl(:,1:3))
    return

900 ierr = -1
910 call reportError (error_p,'Could not calculate coordinate system for '// &
         &            'Rosette'//getId(id),'Check the rosette definition.', &
         &            addString='calcElmCoordSystem')

  end subroutine calcElmCoordSystem


  subroutine InitStrainRosette (rosette, H_el, madof, addSigma0, sigma0, &
       &                        iprint, lpu, ierr, &
       &                        useElCoordSys, openFiles, calcDisp, vgii)

    !!==========================================================================
    !! Compute strain-displacement matrices for all strain rosettes.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 25 May 2001 / 2.0
    !!==========================================================================

    use StrainAndStressUtilitiesModule, only : StrainDispCST, StrainDispQuad4
    use IdTypeModule                  , only : getId
    use fileUtilitiesModule           , only : findUnitNumber
    use manipMatrixModule             , only : writeObject
    use reportErrorModule             , only : allocationError, internalError
    use reportErrorModule             , only : reportError
    use reportErrorModule             , only : error_p, debugFileOnly_p
    use FFaTensorTransformsInterface  , only : tratensor

    type(StrainElementType), intent(inout) :: rosette
    real(rk)               , intent(in)    :: H_el(:,:)
    integer                , intent(in)    :: madof(:)
    logical                , intent(in)    :: addSigma0
    real(dp), optional     , intent(inout) :: sigma0(:)
    integer                , intent(in)    :: iprint, lpu
    integer                , intent(out)   :: ierr
    logical , optional     , intent(in)    :: useElCoordSys, openFiles, calcDisp
    real(rk), optional     , intent(in)    :: vgii(:)

    !! Local variables
    integer               :: i, j, iros, idof, jdof
    integer               :: n, ndim, nElDof, nElNodes, nNDof
    real(dp)              :: b(3), c(2,2), Teps(3,3), T_el(3,3)
    character(len=16)     :: chId
    character(len=64)     :: chName
    real(dp), allocatable :: X_el(:,:), B_el(:,:)
    real(rk), allocatable :: bscr(:,:)

    !! --- Logic section ---

    ndim     = size(H_el,2)
    nElDof   = size(H_el,1)
    nElNodes = size(rosette%globalNodes)

    !! Find max number of nodal DOFs in the element

    nNDof = 0
    do i = 1, nElNodes
       n = madof(rosette%globalNodes(i)+1) - madof(rosette%globalNodes(i))
       nNDof = max(nNDof,n)
    end do

    allocate(X_el(3,nElNodes),B_el(3,nElNodes*nNDof),bscr(3,nElDof),STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('InitStrainRosette 1')
       return
    end if

    !! Calculate position and transformation matrices

    call calcElmCoordSystem (rosette%id, rosette%globalNodes, rosette%posInGl, &
         &                   X_el, T_el, lpu, ierr, useElCoordSys)
    if (ierr /= 0) then
       call reportError (debugFileOnly_p,'InitStrainRosette')
       return
    end if

    c = matmul(T_el(1:2,:),rosette%posInGl(:,1:2))
    Teps(1,1) = c(1,1)*c(1,1)
    Teps(1,2) = c(2,1)*c(2,1)
    Teps(1,3) = c(1,1)*c(2,1)
    Teps(2,1) = c(1,2)*c(1,2)
    Teps(2,2) = c(2,2)*c(2,2)
    Teps(2,3) = c(1,2)*c(2,2)
    Teps(3,1) = 2.0_dp*c(1,1)*c(1,2)
    Teps(3,2) = 2.0_dp*c(2,2)*c(2,1)
    Teps(3,3) = c(1,1)*c(2,2) + c(1,2)*c(2,1)

    !! Transform the given 3D residual stress tensor into the element plane

    if (addSigma0 .and. present(sigma0)) then
       if (size(sigma0) >= 6) then
          call tratensor (3,sigma0,T_el)
          sigma0(3) = sigma0(4)
       else
          ierr = internalError('InitStrainRosette: sigma0 vector too small')
          return
       end if
    end if

    do iros = 1, size(rosette%data)

       !! Form strain-displacement matrix at the strain-gage position

       select case (nElNodes)

       case (3) !! ====== Triangular element ======

          call StrainDispCST (nNDof, X_el(1,:), X_el(2,:), X_el(3,:), T_el, &
               &              rosette%data(iros)%zPos, B_el)

       case (4) !! ====== Quadrilateral element ======

          call StrainDispQuad4 (nNDof, X_el(1,:), X_el(2,:), X_el(3,:), T_el, &
               &                0.0_dp,0.0_dp, rosette%data(iros)%zPos, B_el)

       case default

          ierr = internalError('InitStrainRosette')
          return

       end select

       allocate(rosette%data(iros)%Bcart(3,ndim),STAT=ierr)
       if (ierr /= 0) then
          ierr = allocationError('InitStrainRosette 2')
          return
       end if

       !! Transform strain-displacement matrix to global dof directions.
       !! B_el is also compressed to the correct size if the element is
       !! connected to both 3-DOF and 6-DOF nodes.
       idof = 1
       jdof = 1
       do i = 1, nElNodes
          n = madof(rosette%globalNodes(i)+1) - madof(rosette%globalNodes(i))
          do j = 0, n-1, 3
             B_el(:,idof+j:idof+j+2) = matmul(B_el(:,jdof+j:jdof+j+2),T_el)
          end do
          idof = idof + n
          jdof = jdof + nNDof
       end do

       !! Transform the strain-displacement matrix to rosette directions
       do i = 1, nElDof
          b = matmul(Teps,B_el(:,i))
          do j = 1, 3
             bscr(j,i) = real(b(j),rk)
          end do
       end do

       !! Establish the strain displacement matrix for Cartesian strain in
       !! rosette coordinate system, with respect to the superelement DOFs
       rosette%data(iros)%Bcart = matmul(bscr,H_el)
       if (iprint > 3) then
          call writeObject (rosette%data(iros)%Bcart,lpu, &
               &            'B-matrix for Strain Rosette'//getId(rosette%id))
       end if

       !! Transform the 2D residual stress tensor to rosette directions
       if (addSigma0 .and. present(sigma0)) then
          rosette%data(iros)%sigmaC0 = matmul(Teps,sigma0(1:3))
       end if

       !! Calculate initial strains due to gravitational displacements
       if (present(vgii)) then
          do i = 1, nElNodes
             call calcStrain (rosette%data(iros)%epsCInit, &
                  &           madof(rosette%globalNodes(i)), &
                  &           madof(rosette%globalNodes(i)+1))
          end do
       end if

    end do

    if (size(rosette%data) == 1 .and. present(calcDisp)) then
       if (calcDisp) then
          allocate(rosette%X0(3,nElNodes),rosette%disp(3,nElNodes), &
               &   rosette%T0(3,3),rosette%ur(6),STAT=ierr)
          if (ierr /= 0) then
             ierr = allocationError('InitStrainRosette 3')
             return
          end if
          rosette%X0 = X_el
          rosette%T0 = transpose(T_el)
          rosette%disp = 0.0_rk
          rosette%ur = 0.0_dp
          if (present(vgii)) then
             allocate(rosette%vgii(3,nElNodes),STAT=ierr)
             if (ierr /= 0) then
                ierr = allocationError('InitStrainRosette 3a')
                return
             end if
             do i = 1, nElNodes
                idof = madof(rosette%globalNodes(i))
                rosette%vgii(:,i) = vgii(idof:idof+2)
             end do
          end if

       end if
    end if

    deallocate(X_el,B_el,bscr)

    if (present(openFiles)) then
       if (.not.openFiles) then
          rosette%lpuASCII = 0
          return
       end if
    end if

    !! Open the strain gage file for ASCII output
    rosette%lpuASCII = findUnitNumber(20)
    write(chId,'(I12)') rosette%id%userId
    chName = 'rosette'//trim(adjustl(chId))//'.asc'
    open(rosette%lpuASCII,FILE=trim(chName),STATUS='REPLACE',IOSTAT=ierr)
    if (ierr /= 0) then
       call reportError (error_p,'Could not open rosette output file '//chName,&
            &            addString='InitStrainRosette')
    end if

  contains

    subroutine calcStrain (eps,idof1,idof2)
      real(dp), intent(inout) :: eps(3)
      integer , intent(in)    :: idof1, idof2
      do idof = idof1, idof2-1
         do j = 1, 3
            eps(j) = eps(j) + real(bscr(j,1+idof-idof1)*vgii(idof),dp)
         end do
      end do
    end subroutine calcStrain

  end subroutine InitStrainRosette


  subroutine CalcRosetteDisplacements (rosette, madof, H_el, displ, supTr, &
       &                               lpu, ierr)

    !!==========================================================================
    !! Compute nodal displacements in a strain rosette.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 31 May 2007 / 1.0
    !!==========================================================================

    use StrainAndStressUtilitiesModule, only : getShellElementAxes
    use rotationModule                , only : FFa_glbEulerZYX
    use manipMatrixModule             , only : matmul34
    use reportErrorModule             , only : reportError, debugFileOnly_p

    type(StrainElementType), intent(inout) :: rosette
    integer                , intent(in)    :: madof(:), lpu
    real(rk)               , intent(in)    :: H_el(:,:), displ(:)
    real(dp)               , intent(in)    :: supTr(:,:)
    integer                , intent(out)   :: ierr

    !! Local variables
    integer  :: i, j, n, ldof, nElNodes
    real(dp) :: Xn(3,4), posR(3), Tn(3,3), Tinc(3,3)

    !! --- Logic section ---

    nElNodes = size(rosette%globalNodes)

    posR = 0.0_dp
    ldof = 1
    do i = 1, nElNodes
       n = rosette%globalNodes(i)
       rosette%disp(:,i) = matmul(H_el(ldof:ldof+2,:),displ)
       if (associated(rosette%vgii)) then
          rosette%disp(:,i) = rosette%disp(:,i) + rosette%vgii(:,i)
       end if
       do j = 1, 3
          posR(j) = posR(j) + real(rosette%disp(j,i),dp)
       end do
       ldof = ldof + madof(n+1)-madof(n)
    end do

    !! Updated global strain gage position
    posR = rosette%posInGl(:,4) + posR/real(nElNodes,dp)
    rosette%ur(1:3) = matmul34(supTr,posR)

    !! Updated strain gage orientation
    Xn(:,1:nElNodes) = rosette%X0 + rosette%disp
    call getShellElementAxes (nElNodes, Xn(1,:), Xn(2,:), Xn(3,:), &
         &                    Tn(:,1), Tn(:,2), Tn(:,3), lpu, ierr)
    if (ierr /= 0) then
       call reportError (debugFileOnly_p,'CalcRosetteDisplacements')
       return
    end if

    !! Calculate Euler angles representing the updated global orientation
    Tinc = matmul(Tn,transpose(rosette%T0))
    call FFa_glbEulerZYX (matmul(Tinc,supTr(:,1:3)),rosette%ur(4:6))

  end subroutine CalcRosetteDisplacements

end module StrainRosetteModule
