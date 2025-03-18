!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file finiteElementModule.f90
!> @brief Finite element object data containers and subroutines.

!!==============================================================================
!> @brief Module with data types representing finite element objects.
!>
!> @details This module contains data and subroutines for handling basic
!> finite elements on the system level in a mechanism model.
!> The finite elements are represented by the same supeltypemodule::supeltype
!> data type as the superelements, but with some additional functionalities
!> for treatment of two-noded beams.

module FiniteElementModule

  use TriadTypeModule, only : TriadType, IdType, dp
  use SupElTypeModule, only : SupElType, SupElPtrType

  implicit none

  private

  !> @brief Data type for beam cross section properties.
  type BeamPropertyType
     type(IdType) :: id      !< General identification data
     real(dp) :: Rho         !< Mass density
     real(dp) :: Ix          !< Polar area moment
     real(dp) :: EP(6)       !< Material and cross section parameters
     real(dp) :: CA(2)       !< Shear reduction factors, As = CA*A
     real(dp) :: XS(2)       !< Shear center offset w.r.t. neutral axis
     real(dp) :: Zvec(3)     !< Local Z-definition vector
     real(dp) :: Morison(11) !< Morison parameters (Ca, Cm, Cd, Dd, Db)
     real(dp) :: Dint        !< Diameter of internal fluid body
     real(dp) :: rhoInt      !< Mass density of internal fluid
  end type BeamPropertyType

  !> @brief Data type representing a node in a linked list of beam strings.
  type BeamType
     type(SupElPtrType), pointer :: elm(:) !< Ordered list of beam elements
     type(BeamType)    , pointer :: next   !< Pointer to next beam string node
  end type BeamType

  !> @brief Data type representing a beam joint.
  type BeamJointType
     type(TriadType)   , pointer :: triad   !< The triad the connects the beams
     type(SupElPtrType), pointer :: legs(:) !< List of connected beam elements
  end type BeamJointType

  !> Total number of beam elements in the model
  integer            , public     , save :: nBeamEl
  !> Mass density of the marine growth layer
  real(dp)           , public     , save :: massGrowth
  !> Mass density of the internal fluid
  real(dp)           , public     , save :: massIntFluid
  !> Maximum and minimum beam sectional forces over the time history
  real(dp)           , allocatable, save :: beamEnvelope(:,:,:)
  !> Head of the linked list of beam strings
  type(BeamType)     , pointer    , save :: beams
  !> All beam joints in the model
  type(BeamJointType), allocatable, save :: joints(:)

  !> @brief Returns pointer to object with specified ID.
  interface GetPtrToId
     module procedure GetPtrToIdBeamProp
  end interface

  public :: BeamPropertyType, BeamJointType, BeamType, beams
  public :: ReadElementProps, CalcElementMatrices, GetPtrToId
  public :: reportBeamEndForces, envBeamEndForces, GetSectionalForces
  public :: updateBeamSectionForces, writeBeamJointForces
  public :: findBeams, findBeamJoints, deallocateBeams


contains

  !!============================================================================
  !> @brief Returns pointer to (first) beam property with specified ID.
  !>
  !> @param[in] array Array of finiteelementmodule::beampropertytype objects
  !>                  to search within
  !> @param[in] id Base ID of the object to search for
  !>
  !> @details If the beam property is not found, NULL is returned.
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 26 Sep 2008

  function GetPtrToIdBeamProp (array,id) result(ptr)
    type(BeamPropertyType), pointer :: ptr

    integer               , intent(in)         :: id
    type(BeamPropertyType), intent(in), target :: array(:)

    !! Local variables
    integer :: i

    !! --- Logic section ---

    do i = 1, size(array)
       if (array(i)%id%baseId == id) then
          ptr => array(i)
          return
       end if
    end do

    nullify(ptr)
    write(*,*) '*** GetPtrToIdBeamProp returned nullified, baseId =',id

  end function GetPtrToIdBeamProp


  !!============================================================================
  !> @brief Read finite element properties from the solver input file.
  !>
  !> @param[in] infp File unit number for the solver input file
  !> @param bProp Array of all beam properties in the model
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 26 Sep 2008

  subroutine ReadElementProps (infp,bProp,ierr)

    use idTypeModule     , only : ldesc_p, initId, nullifyId, ReportInputError
    use inputUtilities   , only : iuGetNumberOfEntries, iuSetPosAtNextEntry
    use progressModule   , only : lterm
    use reportErrorModule, only : AllocationError, reportError, debugFileOnly_p

    integer               , intent(in)  :: infp
    type(BeamPropertyType), pointer     :: bProp(:)
    integer               , intent(out) :: ierr

    !! Local variables
    integer :: i, nElmP, stat

    !! Define the ELEMENT_PROPERTY namelist
    integer                :: id, extId(10)
    character(len=ldesc_p) :: extDescr
    real(dp)               :: geometry(8), orientation(3), &
         &                    material(4), hydyn(11), rho_int, D_int
    namelist /ELEMENT_PROPERTY/ id, extId, extDescr, geometry, orientation, &
         &                      material, hydyn, rho_int, D_int

    !! --- Logic section ---

    nElmP = iuGetNumberOfEntries(infp,'&ELEMENT_PROPERTY',ierr)
    if (ierr /= 0) return
    write(lterm,*) 'Number of &ELEMENT_PROPERTY =',nElmP

    allocate(bProp(nElmP),STAT=stat)
    if (stat /= 0) then
       ierr = AllocationError('ReadElementProps')
       return
    end if

    do i = 1, nElmP
       call nullifyId (bProp(i)%id)

       if (.not. iuSetPosAtNextEntry(infp,'&ELEMENT_PROPERTY')) then
          ierr = ierr - 1
          call ReportInputError ('ELEMENT_PROPERTY',i)
          cycle
       end if

       id=0; extId=0; extDescr=''
       geometry=0.0_dp; orientation=0.0_dp; material=0.0_dp; hydyn=0.0_dp
       rho_int=0.0_dp; D_int=0.0
       read(infp,nml=ELEMENT_PROPERTY,iostat=stat)
       if (stat /= 0) then
          ierr = ierr - 1
          call ReportInputError ('ELEMENT_PROPERTY',i)
          cycle
       end if

       call initId (bProp(i)%id,id,extId,extDescr,stat)
       bProp(i)%rho     = material(1)   ! Mass density
       bProp(i)%EP(1:2) = material(2:3) ! Stiffness moduli (E,G)
       bProp(i)%EP(3:6) = geometry(1:4) ! Cross section parameters (A,Iy,Iz,It)
       if (material(2) < 0.0_dp) then   ! Generic cross section data is given
          bProp(i)%Ix   = material(4)   ! Torsional mass inertia
          bProp(i)%CA   = geometry(5:6) ! Shear stiffness (G*Asy,G*Asz)
       else                             ! A, Iy, Iz, etc. are given
          bProp(i)%Ix   = geometry(2) + geometry(3) ! Polar area moment (Ip)
          !! Note that if Iy=geometry(2) and Iz=geometry(3) both are zero, Ix is
          !! assumed equal to the torsional constant It=geometry(4) instead.
          !! This is correct for circular cross sections only.
          if (bProp(i)%Ix <= 0.0_dp) then
             bProp(i)%Ix = geometry(4)
          end if
          bProp(i)%CA   = geometry(5:6) ! Shear reduction factors (kappa_x[yz])
       end if
       bProp(i)%XS      = geometry(7:8) ! Shear centre
       bProp(i)%Zvec    = orientation   ! Local Z-axis
       bProp(i)%Morison = hydyn         ! Hydrodynamic quantities
       bProp(i)%rhoInt  = rho_int       ! Density of internal fluid
       if (rho_int > 0.0_dp) then
          bProp(i)%Dint = D_int         ! Diameter of internal fluid body
       else
          bProp(i)%Dint = 0.0_dp
       end if
    end do

    if (ierr < 0) call reportError (debugFileOnly_p,'ReadElementProps')

  end subroutine ReadElementProps


  !!============================================================================
  !> @brief Calculates the element matrices for a basic finite element.
  !>
  !> @param elem The finite element to calculate element matrices for
  !> @param[in] bProp Beam properties associated with the element
  !> @param[in] env Environmental data
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 26 Sep 2008

  subroutine CalcElementMatrices (elem,bProp,env,ierr)

    use EnvironmentTypeModule , only : EnvironmentType
    use IdTypeModule          , only : getId
    use SupElTypeModule       , only : IsBeam
    use reportErrorModule     , only : getErrorFile, reportError, error_p
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getint

    type(SupElType)       , intent(inout) :: elem
    type(BeamPropertyType), intent(in)    :: bProp
    type(EnvironmentType) , intent(in)    :: env
    integer               , intent(out)   :: ierr

    !! Local variables
    integer :: lpu, iprint

    !! --- Logic section ---

    lpu = getErrorFile()
    call ffa_cmdlinearg_getint ('debug',iprint)

    if (IsBeam(elem)) then
       call CalcBeamElementMatrices (elem,bProp,env,lpu,iprint,ierr)
    else
       ierr = -1
       call reportError (error_p,'Invalid finite element')
    end if
    if (ierr /= 0) then
       call reportError (error_p,'Detected for Element'//getId(elem%id), &
            &            addString='CalcElementMatrices')
    end if

  end subroutine CalcElementMatrices


  !!============================================================================
  !> @brief Calculates the element matrices for a beam finite element.
  !>
  !> @param beam The beam element to calculate element matrices for
  !> @param[in] bProp Beam properties associated with the element
  !> @param[in] env Environmental data
  !> @param[in] lpu File unit number for res-file output
  !> @param[in] ipsw Print switch for debug output
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 26 Sep 2008

  subroutine CalcBeamElementMatrices (beam,bProp,env,lpu,ipsw,ierr)

    use KindModule           , only : pi_p, epsDiv0_p
    use EnvironmentTypeModule, only : EnvironmentType, glob2sea
    use SupElTypeModule      , only : InitiateHyDyn
    use reportErrorModule    , only : AllocationError, reportError, error_p
    use manipMatrixModule    , only : trans1V

    type(SupElType)       , intent(inout) :: beam
    type(BeamPropertyType), intent(in)    :: bProp
    type(EnvironmentType) , intent(in)    :: env
    integer               , intent(in)    :: lpu, ipsw
    integer               , intent(out)   :: ierr

    !! Local variables
    real(dp) :: Xg(5,3), Xl(3), Tlg(3,3), bL, Area, rhoM, rhoG, M
    real(dp) :: Z1, Z2, MorPar(11)
    integer  :: necc

    !! --- Logic section ---

    if (associated(beam%eccVec)) then
       necc = size(beam%eccVec,2)
    else
       necc = 0
    end if

    Xg(4,:) = beam%TrUndeformed(:,4,1)
    Xg(5,:) = beam%TrUndeformed(:,4,2)
    if (necc >= 2) then
       !! Eccentric beam ends w.r.t. the triads
       Xg(1,:) = Xg(4,:) + beam%eccVec(:,1)
       Xg(2,:) = Xg(5,:) + beam%eccVec(:,2)
    else
       Xg(1:2,:) = Xg(4:5,:)
    end if
    if (.not. any(abs(bProp%Zvec) > epsDiv0_p)) then
       Tlg = trans1V(Xg(2,:)-Xg(1,:),lpu,ierr)
       if (ierr < 0) then
          call reportError (error_p,'Degenerated beam element (zero length).', &
               &            addString='CalcBeamElementMatrices')
          return
       end if
       Xg(3,:) = Xg(1,:) + Tlg(:,3)
    else
       Xg(3,:) = Xg(1,:) + bProp%Zvec
    end if

    if (size(beam%KmMat,1) /= 12 .or. size(beam%KmMat,2) /= 12) then
       ierr = -1
       call reportError (error_p,'Invalid dimension on beam element matrix.', &
            'Please note that beam triads can not be attached to ground.', &
            'Fix all DOFs in the triad instead and try again.', &
            addString='CalcBeamElementMatrices')
       return
    end if

    !!==========================================================================
    !! Invoke Femlib routine for calculation of material stiffness matrix.
    !!
    !! SUBROUTINE BEAM31 (EK,X,Y,Z,EP,CA,XS,EFFLEN,PHI,
    !!                    IPINA,IPINB,ID,IW,IPSW,IERR)
    !! EK(1:12,1:12) - Stiffness matrix in global coordinate axes
    !! X(1:5), Y(1:5), Z(1:5) - Coordinates of beam element
    !!       X(1),Y(1),Z(1) - Global coordinates of neutral axis in end 1
    !!       X(2),Y(2),Z(2) - Global coordinates of neutral axis in end 2
    !!       X(3),Y(3),Z(3) - Global coordinates of point in local XZ-plane
    !!       X(4),Y(4),Z(4) - Global coordinates of node 1
    !!       X(5),Y(5),Z(5) - Global coordinates of node 2
    !! EP(1:6) - Material and cross section parameters
    !!       EP(1) = E - Young's modulus
    !!       EP(2) = G - Shear modulus
    !!       EP(3) = A - Cross section area
    !!       EP(4) = Iy - Area moment about local Y-axis
    !!       EP(5) = Iz - Area moment about local Z-axis
    !!       EP(6) = Ix - Torsional constant (polar area moment)
    !! CA(1:2) - Shear reduction factors in local Y- and Z-directions
    !! XS(1:2) - Shear centre Y,Z-coordinates (relative to neutral axis)
    !! EFFLEN - Effective beam length (0.0: computed from coordinates)
    !! PHI - Angle (in degrees) between local Y-axis and principal Y-Axis
    !! IPINA, IPINB - End release flags (0: No end release here)
    !! ID - Element ID used for error messages, etc.
    !! IW - Output unit for debug print
    !! IPSW - Print switch
    !! IERR - Error flag
    !!==========================================================================

    call BEAM31 (beam%KmMat(1,1),Xg(1,1),Xg(1,2),Xg(1,3), &
         &       bProp%EP(1),bProp%CA(1),bProp%XS(1),0.0_dp,0.0_dp,0,0, &
         &       beam%id%userId,lpu,ipsw,ierr)
    if (ierr < 0) then
       call reportError (error_p,'Failed to compute beam stiffness', &
            &            addString='CalcBeamElementMatrices')
       return
    else if (any(isNan(beam%KmMat))) then
       call reportError (error_p,'Beam element stiffness matrix contains NaN', &
            &            'Check the material/geometry properties', &
            &            addString='CalcBeamElementMatrices')
       write(lpu,100) 'EP',bProp%EP
       write(lpu,100) 'CA',bProp%CA
       write(lpu,100) 'XS',bProp%XS
100    format(10X,A2,' =',1P9E13.5)
       ierr = -9
       return
    end if

    !! Find the effective mass densities for this beam element
    rhoM = bProp%rho ! For mass matrix
    rhoG = bProp%rho ! For gravity forces

    if (associated(beam%hydyn)) then
       Xl = Xg(2,:) - Xg(1,:)
       bL = sqrt(Xl(1)*Xl(1) + Xl(2)*Xl(2) + Xl(3)*Xl(3))
       if (beam%hydyn%bodyIndex == -2 .and. abs(bProp%Dint) > epsDiv0_p) then

          !! Cross section area of the internal fluid body
          Area = 0.25_dp*pi_p*bProp%Dint*bProp%Dint
          !! Add total mass from internal fluid in this element
          M = bProp%rhoInt*Area
          massIntFluid = massIntFluid + M*bL

          !! Modify the mass density to account for internal fluid
          if (bProp%EP(1) < 0.0_dp .or. bProp%EP(3) <= epsDiv0_p) then
             rhoM = rhoM + M             ! Mass per unit length
          else
             rhoM = rhoM + M/bProp%EP(3) ! Mass per unit volume
          end if
          if (bProp%Dint > 0.0_dp) then
             rhoG = rhoM ! Include the internal fluid mass as gravity load too
          end if

       end if
       if (env%tG > 0.0_dp .and. env%zg(1) < env%zg(2)) then

          !! Check for marine growth on this element
          Xl = glob2sea(env,beam%triads(1)%p%ur(:,4))
          Z1 = Xl(3) ! Z-coordinate of end 1 with respect to MSL
          Xl = glob2sea(env,beam%triads(2)%p%ur(:,4))
          Z2 = Xl(3) ! Z-coordinate of end 2 with respect to MSL
          if (Z1 > Z2) then
             Z2 = Z1
             Z1 = Xl(3)
          end if
          Area = hasGrowth(Z1,Z2,env%zg(1),env%zg(2))
          if (Area > 1.0e-4_dp) then
             !! This element has marine growth,
             !! modify the drag and buoyancy diameters accordingly
             MorPar = bprop%Morison
             MorPar(4) = MorPar(4) + env%tG + env%tG
             MorPar(5) = MorPar(5) + env%tG + env%tG
             call initiateHyDyn (beam%hydyn,MorPar)
             !! Cross section area of the marine growth
             Area = pi_p*(bProp%Morison(5)+env%tG)*env%tG * Area
             !! Add total mass from marine growth on this element
             M = env%rhoG*Area
             massGrowth = massGrowth + M*bL

             !! Modify mass density to account for marine growth
             if (bProp%EP(1) < 0.0_dp .or. bProp%EP(3) <= epsDiv0_p) then
                rhoM = rhoM + M
                rhoG = rhoG + M
             else
                rhoM = rhoM + M/bProp%EP(3)
                rhoG = rhoG + M/bProp%EP(3)
             end if
          end if

       end if
    end if

    if (bProp%EP(1) <= epsDiv0_p .and. rhoG > epsDiv0_p) then
       Area = 1.0_dp      ! Element mass is given as mass per unit length
       z1   = bProp%Ix/rhoG
       z2   = bProp%Ix/rhoM
    else
       Area = bProp%EP(3) ! Element mass is given as mass per unit volume
       z1   = bProp%Ix
       z2   = bProp%Ix
    end if
    if (necc >= 4) then
       !! Eccentric mass w.r.t. the triads
       Xg(1,:) = Xg(1,:) + beam%eccVec(:,3)
       Xg(2,:) = Xg(2,:) + beam%eccVec(:,4)
    end if

    !!==========================================================================
    !! Invoke Femlib routine for calculation of consistent mass matrix.
    !!
    !! SUBROUTINE BEAM35 (EM,X,Y,Z,A,RHO,RIX,IOP,IPINA,IPINB,IW,IPSW,IERR)
    !! EM(1:12,1:12) - Mass matrix in global coordinate axes
    !! X(1:5), Y(1:5), Z(1:5) - Coordinates of beam element
    !!       X(1),Y(1),Z(1) - Global coordinates of CoG in end 1
    !!       X(2),Y(2),Z(2) - Global coordinates of CoG in end 2
    !!       X(3),Y(3),Z(3) - Direction vector for local Z-axis
    !!       X(4),Y(4),Z(4) - Global coordinates of node 1
    !!       X(5),Y(5),Z(5) - Global coordinates of node 2
    !! A - Cross section area
    !! RHO - Mass density
    !! RIX - Polar area moment
    !! IOP - Mass matrix option (3: Consistent mass matrix)
    !! IPINA, IPINB - End release flags (0: No end release here)
    !! IW - Output unit for debug print
    !! IPSW - Print switch
    !! IERR - Error flag
    !!==========================================================================

    call BEAM35 (beam%Mmat(1,1),Xg(1,1),Xg(1,2),Xg(1,3),Area,rhoG,z1,3,0,0, &
         &       lpu,ipsw,ierr)

    beam%fg = beam%Mmat(:,1:3) + beam%Mmat(:,7:9)

    if (abs(rhoM-rhoG) > epsDiv0_p) then
       !! Recalculate the mass matrix including internal fluid mass also
       call BEAM35 (beam%Mmat(1,1),Xg(1,1),Xg(1,2),Xg(1,3),Area,rhoM,z2,3,0,0, &
            &       lpu,ipsw,ierr)
    end if

    !! Decouple the nodeForce from the supForce array again, since each triad
    !! now may be connected to several elements. The nodeForce (if needed)
    !! will then be the sum of contributions from all the connected elements,
    !! and any other connected objects (springs/dampers/masses/tires/forces).
    !! The supForce will now contain the section force instead, calculated
    !! from the element referring to that triad as its first local node.
    if (associated(beam%triads(1)%p%nodeForce,beam%triads(1)%p%supForce)) then
       nullify(beam%triads(1)%p%nodeForce)
    end if
    if (associated(beam%triads(2)%p%nodeForce,beam%triads(2)%p%supForce)) then
       nullify(beam%triads(2)%p%nodeForce)
    end if

    !! Signal that these triads now represent FE nodes with sectional forces
    !! in local (element) directions stored in the supForce array
    beam%triads(1)%p%isFEnode = .true.
    beam%triads(2)%p%isFEnode = .true.

    !! Initialize the beam bending stiffness that is used in the alternative
    !! bending moment calculation based on curvatures derived from the
    !! deformed beam string configuration (see updateBeamChainQuantities).
    allocate(beam%EI(2),STAT=ierr)
    if (ierr /= 0) then
       ierr = AllocationError('CalcBeamElementMatrices')
    else if (bProp%EP(1) < 0.0_dp) then
       beam%EI = bProp%EP(4:5)
    else
       beam%EI = bProp%EP(1)*bProp%EP(4:5)
    end if

  contains

    function hasGrowth (z1,z2,zg1,zg2)
      real(dp), intent(in) :: z1, z2, zg1, zg2
      real(dp)             :: hasGrowth
      if (z1 >= zg2 .or. z2 <= zg1) then
         hasGrowth = 0.0_dp
      else if (z1 < zg1 .and. z2 <= zg2) then
         hasGrowth = (z2-zg1)/(z2-z1)
      else if (z1 >= zg1 .and. z2 > zg2) then
         hasGrowth = (zg2-z1)/(z2-z1)
      else if (z1 < zg1 .and. z2 > zg2) then
         hasGrowth = (zg2-zg1)/(z2-z1)
      else
         hasGrowth = 1.0_dp
      end if
    end function hasGrowth

  end subroutine CalcBeamElementMatrices


  !!============================================================================
  !> @brief Detects all beams in the model, that is, chains of beam elements.
  !>
  !> @param[in] sups All superelements in the model
  !> @param[out] ierr Error flag
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 12 Oct 2011

  subroutine findBeams (sups,ierr)

    use IdTypeModule     , only : getId, StrId
    use SupElTypeModule  , only : IsBeam
    use progressModule   , only : lterm
    use reportErrorModule, only : allocationError, reportError, note_p

    type(SupElType), intent(in), target :: sups(:)
    integer        , intent(out)        :: ierr

    !! Local variables
    integer                 :: i, n2, ilast, istart, nel, nBeams
    type(BeamType), pointer :: current

    !! --- Logic section ---

    ierr = 0
    nBeamEl = 0
    nullify(beams)
    if (size(sups) < 1) then
       return
    else if (IsBeam(sups(1))) then
       nBeamEl = 1
       istart = 1
       ilast = 1
       nel = 1
    else
       istart = 0
       ilast = 0
       nel = 0
    end if

    n2 = 2
    nBeams = 0
    do i = 2, size(sups)
       if (IsBeam(sups(i))) then
          nBeamEl = nBeamEl + 1
          if (istart < 1) then ! First element in a new beam
             istart = i
             ilast = i
             nel = 1
             cycle
          else if (i > ilast+1) then ! We have reached the end of a beam
             if (nel > 1) then ! Only consider beams of two or more elements
                call newBeam
                if (ierr /= 0) then
                   ierr = allocationError('findBeams')
                   return
                end if
             end if
             n2 = 2
             !! Must check if the current and/or previous element
             !! are parts of a new beam before continuing the loop
             if (i == ilast+2 .and. IsBeam(sups(ilast+1))) then
                !! Starting a new beam from the previous element
                istart = ilast+1
                ilast = istart
                nel = 1
             else if (IsBeam(sups(i))) then
                !! Starting a new beam from the current element
                istart = i
                ilast = i
                nel = 1
                cycle
             else
                istart = 0
                cycle
             end if
          end if
          if (associated(sups(i)%triads(1)%p,sups(i-1)%triads(n2)%p)) then
             !! Add a new element to current beam
             nel = nel + 1
             ilast = i
             n2 = 2
          else if (associated(sups(i)%triads(2)%p,sups(i-1)%triads(n2)%p)) then
             !! Add a new element to current beam (swapped nodes)
             nel = nel + 1
             ilast = i
             n2 = 1
          end if
       end if
    end do

    if (istart > 0 .and. nel > 1) then
       call newBeam
       if (ierr /= 0) then
          ierr = allocationError('findBeams')
       end if
    end if

    if (nBeams > 0) then
       write(lterm,*) nBeams,' beams detected (',nBeamEl,' elements)'
    end if

    current => beams
    do while (associated(current))
       call reportError (note_p,'Structural beam: Beam'// &
            &            trim(getId(current%elm(1)%p%id))//' -- Beam'// &
            &            trim(getId(current%elm(size(current%elm))%p%id))// &
            &            ','//trim(StrId(size(current%elm)))//' elements')
       current => current%next
    end do

  contains

    !> @brief Create a new node in the linked list of beams.
    subroutine newBeam
      integer :: j, k
      if (associated(beams)) then
         !! Allocate a new entry in the linked list of beams
         allocate(current%next,stat=ierr)
         if (ierr /= 0) return
         current => current%next
      else
         !! This is the first beam encountered
         allocate(beams,stat=ierr)
         if (ierr /= 0) return
         current => beams
      end if
      allocate(current%elm(nel),stat=ierr)
      if (ierr /= 0) return
      nullify(current%next)
      do j = 1, nel
         current%elm(j)%p => sups(istart+j-1)
         allocate(current%elm(j)%p%scoord,stat=ierr)
         if (ierr /= 0) return
         current%elm(j)%p%scoord = 0.0_dp
         do k = 1, 2
            if (.not. associated(current%elm(j)%p%triads(k)%p%scoord)) then
               allocate(current%elm(j)%p%triads(k)%p%scoord, &
                    &   current%elm(j)%p%triads(k)%p%kappa(3), &
                    &   current%elm(j)%p%triads(k)%p%Moment(2),stat=ierr)
               if (ierr /= 0) return
               current%elm(j)%p%triads(k)%p%scoord = 0.0_dp
               current%elm(j)%p%triads(k)%p%kappa  = 0.0_dp
               current%elm(j)%p%triads(k)%p%Moment = 0.0_dp
            end if
         end do
      end do
      nbeams = nbeams + 1
    end subroutine newBeam

  end subroutine findBeams


  !!============================================================================
  !> @brief Detects all beam joints in the model.
  !>
  !> @param[in] triads All triads in the model
  !> @param sups All superelements in the model
  !> @param[out] ierr Error flag
  !>
  !> @details Beam joints are triads connected to three (or more) beam elements.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 3 Nov 2010

  subroutine findBeamJoints (triads,sups,ierr)

    use TriadTypeModule       , only : GetPtrToId
    use SupElTypeModule       , only : IsBeam
    use progressModule        , only : lterm
    use reportErrorModule     , only : allocationError
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getbool

    type(TriadType), intent(in)    :: triads(:)
    type(SupElType), intent(inout) :: sups(:)
    integer        , intent(out)   :: ierr

    !! Local variables
    integer                      :: i, j, k, nTriad, nSup, nJoint
    logical                      :: allBeamForces, noBeamForces
    integer        , allocatable :: nLeg(:)
    type(TriadType), pointer     :: triad

    !! --- Logic section ---

    nTriad = size(triads)
    allocate(nLeg(nTriad),STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('findBeamJoints 1')
       return
    end if

    !! Find the number beam elements (legs) connected to each triad
    nLeg = 0
    nSup = size(sups)
    do i = 1, nSup
       if (IsBeam(sups(i))) then
          triad => GetPtrToId(triads,sups(i)%triads(1)%p%id%baseId,j)
          if (associated(triad)) nLeg(j) = nLeg(j) + 1
          triad => GetPtrToId(triads,sups(i)%triads(2)%p%id%baseId,j)
          if (associated(triad)) nLeg(j) = nLeg(j) + 1
       end if
    end do

    !! Triads with 3 or more connections defines a joint
    nJoint = count(nLeg > 2)
    if (nJoint < 1) goto 100

    write(lterm,*) nJoint,' beam joint(s) detected'

    allocate(joints(nJoint),STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('findBeamJoints 2')
       return
    end if

    !! Establish the joint topology objects and transform the nLeg array
    !! into a triad-to-joint index map
    j = 1
    do i = 1, nJoint
       do while (nLeg(j) < 3)
          nLeg(j) = 0
          j = j + 1
       end do
       allocate(joints(i)%legs(nLeg(j)),STAT=ierr)
       if (ierr /= 0) then
          ierr = allocationError('findBeamJoints 3')
          return
       end if
       do k = 1, nLeg(j)
          nullify(joints(i)%legs(k)%p)
       end do
       nullify(joints(i)%triad)
       nLeg(j) = i
       j = j + 1
    end do
    if (j <= nTriad) nLeg(j:nTriad) = 0

    !! Set up the beam joint topology and modify the rigidFlag:
    !! rigidFlag = -1: Zero or one element connection at each end
    !! rigidFlag = -2: Two or more connections and end 1, one or less at end 2
    !! rigidFlag = -3: Two or more connections and end 2, one or less at end 1
    !! rigidFlag = -4: Two or more connections and both ends
    !! rigidFlag = -5: We don't care, flag all beams for sectional force output
100 call ffa_cmdlinearg_getbool ('allBeamForces',allBeamForces)
    call ffa_cmdlinearg_getbool ('noBeamForces',noBeamForces)
    do i = 1, nSup
       if (IsBeam(sups(i))) then
          if (nJoint > 0) then
             triad => GetPtrToId(triads,sups(i)%triads(1)%p%id%baseId,j)
             if (associated(triad) .and. nLeg(j) > 0) then
                call initJointLeg (joints(nLeg(j)),sups(i))
                sups(i)%rigidFlag = -2
             else
                sups(i)%rigidFlag = -1
             end if
             triad => GetPtrToId(triads,sups(i)%triads(2)%p%id%baseId,j)
             if (associated(triad) .and. nLeg(j) > 0) then
                call initJointLeg (joints(nLeg(j)),sups(i))
                sups(i)%rigidFlag = sups(i)%rigidFlag - 2
             end if
          end if
          if (allBeamForces) then
             sups(i)%rigidFlag = -5
          else if (noBeamForces .or. nJoint < 1) then
             sups(i)%rigidFlag = -1
          end if
          if (sups(i)%rigidFlag < -1) then
             sups(i)%saveVar(2) = .true.
          end if
       end if
    end do

    deallocate(nLeg)

  contains

    !> @brief Adds a beam element to a beam joint,
    subroutine initJointLeg (joint,beam)
      type(BeamJointType)    , intent(inout) :: joint
      type(SupElType), target, intent(in)    :: beam
      if (associated(joint%triad)) then
         do k = 2, size(joint%legs)
            if (.not. associated(joint%legs(k)%p)) then
               joint%legs(k)%p => beam
               return
            end if
         end do
      else
         joint%triad => triad
         joint%legs(1)%p => beam
      end if
    end subroutine initJointLeg

  end subroutine findBeamJoints


  !!============================================================================
  !> @brief Deallocates the linked lists of beam elements.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Jan 2017

  subroutine deallocateBeams ()

    !! Local variables
    integer                 :: i
    type(BeamType), pointer :: current, p

    !! --- Logic section ---

    if (allocated(beamEnvelope)) then
       deallocate(beamEnvelope)
    end if

    if (allocated(joints)) then
       do i = 1, size(joints)
          deallocate(joints(i)%legs)
       end do
       deallocate(joints)
    end if

    current => beams
    do while (associated(current))
       p => current%next
       deallocate(current%elm)
       deallocate(current)
       current => p
    end do

    nullify(beams)
    nBeamEl = 0

  end subroutine deallocateBeams


  !!============================================================================
  !> @brief Returns the sectional forces at one end of a beam element.
  !>
  !> @param[in] sup The beam element to calculated sectional forces for
  !> @param[in] iend Which end to calculate sectional forces at (1 or 2)
  !>
  !> @details The forces are calculated in the local coordinate system of the
  !> two-noded beam element.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 16 Sep 2009

  function GetSectionalForces (sup,iend) result(secf)

    use SupElTypeModule,        only : GetRsupToSys
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_isTrue

    type(SupElType), intent(in) :: sup
    integer        , intent(in) :: iend

    !! Local variables
    integer  :: idof, jdof
    real(dp) :: secf(6), T_lg(3,3)

    !! --- Logic section ---

    secf = 0.0_dp
    if (iend < 1 .or. iend > size(sup%triads)) return
    idof = sup%triads(iend)%firstDOF
    jdof = idof + sup%triads(iend)%p%nDOFs-1
    if (ffa_cmdlinearg_isTrue('subtractExtFromIntForces')) then
       secf = sup%Fs(idof:jdof) - sup%Q(idof:jdof)
    else
       secf = sup%Fs(idof:jdof)
    end if
    T_lg = GetRsupToSys(sup,iend)
    secf(1:3) = matmul(secf(1:3),T_lg)
    secf(4:6) = matmul(secf(4:6),T_lg)

  end function GetSectionalForces


  !!============================================================================
  !> @brief Reports all beam sectional forces at current configuration.
  !>
  !> @param[in] t Current simulation time,
  !> &ge;0: Report forces at current time step,
  !> &lt;0: Report current minimum and maximum forces (envelope).
  !> @param[in] sups All superelements in the model
  !> @param[in] g Gravitation vector
  !> @param[in] lpu File unit number for res-file output
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 13 Nov 2008

  subroutine reportBeamEndForces (t,sups,g,lpu)

    use kindModule         , only : epsDiv0_p
    use SupElTypeModule    , only : IsBeam
    use IdTypeModule       , only : getId
    use FileUtilitiesModule, only : findUnitNumber

    real(dp)          , intent(in) :: t
    type(SupElType)   , intent(in) :: sups(:)
    real(dp), optional, intent(in) :: g(3)
    integer           , intent(in) :: lpu

    !! Local variables
    integer  :: i, j, iop, iun
    real(dp) :: n(3), z, sf(2,6)

    !! --- Logic section ---

    write(lpu,600) abs(t)
    if (t >= 0.0_dp) then
       iop = 1
    else if (allocated(beamEnvelope)) then
       iop = 2
    else
       return
    end if

    n = 0.0_dp
    if (present(g)) then
       !! Print out height in the (opposite) direction of the gravity vector
       z = sqrt(g(1)*g(1) + g(2)*g(2) + g(3)*g(3))
       if (z > epsDiv0_p) then
          n = -g/z
       else
          n(3) = 1.0_dp
       end if
       if (iop == 1) then
          write(lpu,601) 'Height'
       else
          write(lpu,602) 'Height'
       end if
    else
       !! Print out accumulated position coordinate along the beam
       if (iop == 1) then
          write(lpu,601) 'Position'
       else
          write(lpu,602) 'Position'
       end if
    end if

    j = 0
    z = 0.0_dp
    do i = 1, size(sups)
       if (IsBeam(sups(i))) then
          !! Print forces at the first triad of each beam element
          if (iop == 1) then
             sf(1,:) = -GetSectionalForces(sups(i),1)
          else
             sf = transpose(beamEnvelope(1:6,:,i))
          end if
          if (present(g)) z = dot_product(n,sups(i)%triads(1)%p%ur(:,4))
          write(lpu,612) getId(sups(i)%id), z, sf(1:iop,:)
          if (.not. present(g)) then
             n = sups(i)%triads(2)%p%ur(:,4) - sups(i)%triads(1)%p%ur(:,4)
             z = z + sqrt(n(1)*n(1) + n(2)*n(2) + n(3)*n(3))
          end if
          j = i
       end if
    end do
    if (j > 0) then
       !! Print forces at the second triad of the last beam element
       if (iop == 1) then
          sf(1,:) = GetSectionalForces(sups(j),2)
       else
          sf = transpose(beamEnvelope(1:6,:,j+1))
       end if
       if (present(g)) z = dot_product(n,sups(j)%triads(2)%p%ur(:,4))
       write(lpu,612) getId(sups(j)%id), z, sf(1:iop,:)
    end if

    if (iop == 1) return

    iun = findUnitNumber(10)
    open(iun,FILE='riser_envelope.asc')

    write(lpu,620)
    if (present(g)) then
       write(lpu,622) 'Height'
       write(iun,623) 'Height'
    else
       write(lpu,622) 'Position'
       write(iun,623) 'Position'
    end if

    j = 0
    z = 0.0_dp
    do i = 1, size(sups)
       if (IsBeam(sups(i))) then
          !! Print position at the first triad of each beam element
          if (present(g)) z = dot_product(n,sups(i)%triads(1)%p%ur(:,4))
          write(lpu,611) getId(sups(i)%id), z, (beamEnvelope(iop,:,i),iop=7,9)
          write(iun,613) z, (beamEnvelope(iop,:,i),iop=1,9)
          if (.not. present(g)) then
             n = sups(i)%triads(2)%p%ur(:,4) - sups(i)%triads(1)%p%ur(:,4)
             z = z + sqrt(n(1)*n(1) + n(2)*n(2) + n(3)*n(3))
          end if
          j = i
       end if
    end do

    !! Print position at the second triad of the last beam element
    if (present(g)) z = dot_product(n,sups(j)%triads(2)%p%ur(:,4))
    write(lpu,611) getId(sups(j)%id), z, (beamEnvelope(iop,:,j+1),iop=7,9)
    write(iun,613) z, (beamEnvelope(iop,:,j+1),iop=1,9)

    close(iun)
    deallocate(beamEnvelope)

600 format(//5X,'>>> Summary of beam sectional forces, at t =',1PE12.5,' <<<' &
         &  /5X,60('='))
601 format( /5X,'Beam identification',11X,A8, &
         &   8X,'Nx',11X,'Ny',11X,'Nz',11X,'Mx',11X,'My',11X,'Mz')
602 format( /5X,'Beam identification',11X,A8, &
         &   6X,'min Nx',7X,'max Nx',7X,'min Ny',7X,'max Ny', &
         &   7X,'min Nz',7X,'max Nz',7X,'min Mx',7X,'max Mx', &
         &   7X,'min My',7X,'max My',7X,'min Mz',7X,'max Mz')
611 format(  5X,A29,F9.3,2X,1P,6E13.5)
612 format(  5X,A29,F9.3,2X,1P,12E13.5)
613 format(  F10.3,1P,18E13.5)

620 format( /5X,'>>> Envelopes for nodal positions along the beam <<<' &
         &  /5X,52('='))
622 format( /5X,'Beam identification',11X,A8, &
         &   6X,'min X',8X,'max X',8X,'min Y',8X,'max Y',8X,'min Z',8X,'max Z')
623 format('# Beam envelopes'/'# ',A8, &
         &   4X,'min_Nx',7X,'max_Nx',7X,'min_Ny',7X,'max_Ny', &
         &   7X,'min_Nz',7X,'max_Nz',7X,'min_Mx',7X,'max_Mx', &
         &   7X,'min_My',7X,'max_My',7X,'min_Mz',7X,'max_Mz', &
         &   7X,'min_X',8X,'max_X',8X,'min_Y',8X,'max_Y',8X,'min_Z',8X,'max_Z')

  end subroutine reportBeamEndForces


  !!============================================================================
  !> @brief Calculates the max/min beam sectional forces and positions.
  !>
  !> @param[in] sups All superelements in the model
  !> @param[out] ierr Error flag
  !>
  !> @details The calculation is performed over the entire time history
  !> of the dynamics simulation and is stored in the @a beamEnvelope array.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 27 Apr 2009

  subroutine envBeamEndForces (sups,ierr)

    use kindModule       , only : hugeVal_p
    use SupElTypeModule  , only : IsBeam
    use reportErrorModule, only : allocationError

    type(SupElType), intent(in)  :: sups(:)
    integer        , intent(out) :: ierr

    !! Local variables
    integer  :: i, j
    real(dp) :: sf(9)

    !! --- Logic section ---

    ierr = 0
    if (nBeamEl < 1) return ! No beam elements in this model

    if (.not. allocated(beamEnvelope)) then
       !! Allocate array for accumulation of max/min forces and positions
       allocate(beamEnvelope(9,2,size(sups)+1), STAT=ierr)
       if (ierr /= 0) then
          ierr = allocationError('envBeamEndForces')
          return
       end if
       beamEnvelope(:,1,:) =  hugeVal_p
       beamEnvelope(:,2,:) = -hugeVal_p
    end if

    j = 0
    do i = 1, size(sups)
       if (IsBeam(sups(i))) then
          !! Calculate sectional forces and nodal position
          !! at the first triad of each beam element
          sf(1:6) = -GetSectionalForces(sups(i),1)
          sf(7:9) = sups(i)%triads(1)%p%ur(:,4)
          beamEnvelope(:,1,i) = min(beamEnvelope(:,1,i),sf)
          beamEnvelope(:,2,i) = max(beamEnvelope(:,2,i),sf)
          j = i
       end if
    end do

    if (j > 0) then
       !! Calculate sectional forces and nodal position
       !! at the second triad of the last beam element
       sf(1:6) = GetSectionalForces(sups(j),2)
       sf(7:9) = sups(j)%triads(2)%p%ur(:,4)
       beamEnvelope(:,1,j+1) = min(beamEnvelope(:,1,j+1),sf)
       beamEnvelope(:,2,j+1) = max(beamEnvelope(:,2,j+1),sf)
    else
       !! This model did not have any beam elements.
       !! Release the allocated array to avoid unneccesary loops later.
       deallocate(beamEnvelope)
       nBeamEl = 0
    end if

  end subroutine envBeamEndForces


  !!============================================================================
  !> @brief Calculates the current beam sectional forces.
  !>
  !> @param sups All superelements in the model
  !> @param[out] ierr Error flag
  !>
  !> @details The running curve coordinate along the beams is alos updated.
  !> This coordinate runs as long as the beam elements are in consequtive order
  !> with no other superelements in between.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 14 Dec 2012

  subroutine updateBeamSectionForces (sups,ierr)

    use kindModule            , only : epsDiv0_p
    use SupElTypeModule       , only : IsBeam
    use IdTypeModule          , only : getId
    use manipMatrixModule     , only : cross_product
    use reportErrorModule     , only : reportError, error_p, debugFileOnly_p
    use reportErrorModule     , only : internalError
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getbool

    type(SupElType), intent(inout) :: sups(:)
    integer        , intent(out)   :: ierr

    !! Local variables
    integer  :: i, nBeam, n1
    logical  :: use5points
    real(dp) :: dl, s0, s1, s2, s3, T(3), N(3), S(3), Trot(3,3), Tcrv(3,3)
    real(dp) :: dX(3), x0(3), x1(3), x2(3), x3(3), x4(3), tolTan

    type(TriadType), pointer :: p0, p1, p2, p3, p4
    type(BeamType) , pointer :: current

    !! --- Logic section ---

    call ffa_cmdlinearg_getbool ('curvature5p',use5points)

    !! Calculation of curve lengths and sectional forces along the beam

    ierr = 0
    nBeam = 0
    tolTan = 1.0_dp
    current => beams
    do while (associated(current))
       p2 => current%elm(1)%p%triads(2)%p
       p3 => current%elm(2)%p%triads(1)%p
       p4 => current%elm(2)%p%triads(2)%p
       do i = 1, size(current%elm)
          nBeam = nBeam + 1
          current%elm(i)%p%rigidFlag = -10 ! Flag this as a beam chain element
          p1 => current%elm(i)%p%triads(1)%p
          p2 => current%elm(i)%p%triads(2)%p
          dX = p2%ur(:,4) - p1%ur(:,4)
          dl = sqrt(dX(1)*dX(1) + dX(2)*dX(2) + dX(3)*dX(3))
          if (dl <= epsDiv0_p) then
             ierr = ierr - 1
             call reportError (error_p,'Degenerated Beam'// &
                  &            getId(current%elm(i)%p%id))
          else if (dl < tolTan) then
             tolTan = dl
          end if
          !! Check for swapped nodes
          n1 = 1
          if (i == 1) then
             if (associated(p1,p3)) then
                p1 => p2
                p2 => p3
                n1 = 2
             else if (associated(p1,p4)) then
                p1 => p2
                p2 => p4
                n1 = 2
             end if
          else if (associated(p2,p4)) then
             p2 => p1
             p1 => p4
             n1 = 2
          end if
          p2%scoord = p1%scoord + dl
          current%elm(i)%p%scoord = p1%scoord + 0.5_dp*dl
          if (associated(p1%supForce) .and. p1%isFEnode) then
             p1%supForce = -GetSectionalForces(current%elm(i)%p,n1)
          end if
          p3 => p1
          p4 => p2
       end do
       if (associated(p2%supForce) .and. p2%isFEnode) then
          n1 = 3 - n1
          p2%supForce = GetSectionalForces(current%elm(size(current%elm))%p,n1)
       end if
       current => current%next
    end do

    !! Calculation of curvature directly from the beam shape

    if (ierr == 0) then
       current => beams
       tolTan = tolTan*epsDiv0_p*epsDiv0_p
    end if
    do while (associated(current))
       p1 => current%elm(1)%p%triads(1)%p
       p2 => current%elm(1)%p%triads(2)%p
       p3 => current%elm(2)%p%triads(1)%p
       p4 => current%elm(2)%p%triads(2)%p
       !! Check if the first element has swapped nodes
       if (associated(p1,p3)) then
          p1 => p2
          p2 => p3
       else if (associated(p1,p4)) then
          p1 => p2
          p2 => p4
       end if
       p0 => p1
       do i = 2, size(current%elm)
          if (associated(current%elm(i)%p%triads(1)%p,p2)) then
             p3 => current%elm(i)%p%triads(2)%p
          else if (associated(current%elm(i)%p%triads(2)%p,p2)) then
             !! Swapped nodes for this element
             p3 => current%elm(i)%p%triads(1)%p
          else
             !! Logic error, should not happen
             ierr = internalError('updateBeamSectionForces 1')
             return
          end if
          x0 = p0%ur(:,4)
          x1 = p1%ur(:,4)
          x2 = p2%ur(:,4)
          x3 = p3%ur(:,4)
          if (use5points .and. i < size(current%elm)) then
             if (associated(current%elm(i+1)%p%triads(1)%p,p3)) then
                p4 => current%elm(i+1)%p%triads(2)%p
             else if (associated(current%elm(i+1)%p%triads(2)%p,p3)) then
                !! Swapped nodes for next element
                p4 => current%elm(i+1)%p%triads(1)%p
             else
                !! Logic error, should not happen
                ierr = internalError('updateBeamSectionForces 2')
                return
             end if
          end if
          x4 = p4%ur(:,4)
          s0 = p1%scoord - p0%scoord
          s1 = p2%scoord - p1%scoord
          s2 = p3%scoord - p2%scoord
          s3 = p4%scoord - p3%scoord

          !! Curve tangent in global coordinates
          dX = s1*s1*x3 + (s2*s2-s1*s1)*x2 - s2*s2*x1
          dl = sqrt(dX(1)*dX(1) + dX(2)*dX(2) + dX(3)*dX(3))
          if (dl <= tolTan) then
             ierr = ierr - 1
             call reportError (error_p,'Unable to calculate curvature at '// &
                  &            'Triad'//getId(p2%id),'Degenerated beam')
          else
             T = dX / dl
          end if
          !! Curvature in global coordinates
          if (use5points .and. i > 2 .and. i < size(current%elm)) then
             !! We have four equal(?) segments, use 5-point stencil
             !! See https://en.wikipedia.org/wiki/Five-point_stencil
             dl = s0+s1+s2+s3
             dX = (16.0_dp*(x3+x1) - 30.0_dp*x2 - x4 - x0) / (0.75_dp*dl*dl)
          else
             !! Use generalized 3-point stencil
             dX = 2.0_dp * (s2*x1 - (s1+s2)*x2 + s1*x3) / ((s1*s2)*(s1+s2))
          end if
          dl = sqrt(dX(1)*dX(1) + dX(2)*dX(2) + dX(3)*dX(3))
          if (dl <= epsDiv0_p) then
             !! The beam is perfectly straight at this point ==> zero curvature
             p2%kappa = 0.0_dp
             p2%Moment = 0.0_dp
          else
             N = dX / dl
             S = cross_product(N,T)
             !! Rotate the element axes about the S-axis such that the
             !! local X-axis coincides with the tangent direction, T
             Tcrv(:,1) = T
             Tcrv(:,2) = S
             Tcrv(:,3) = N
             Trot(1,:) = current%elm(i)%p%supTr(:,1)
             Trot(2,:) = S
             Trot(3,:) = cross_product(Trot(1,:),S)
             Trot = matmul(matmul(Tcrv,Trot),current%elm(i)%p%supTr(:,1:3))
             !! Find curvature in local element directions
             p2%kappa = matmul(dX,Trot)
             p2%kappa(1) = 0.0_dp
             p2%Moment = current%elm(i)%p%EI*p2%kappa(2:3)
          end if
          p0 => p1
          p1 => p2
          p2 => p3
       end do
       current => current%next
    end do

    if (ierr < 0) call reportError (debugFileOnly_p,'updateBeamSectionForces')

    if (nBeam >= size(sups)) return ! Pure beamstring model

    !! Calculation for the remaining beams that are not in a chain.

    do i = 1, size(sups)
       if (IsBeam(sups(i)) .and. sups(i)%rigidFlag > -10) then
          p1 => sups(i)%triads(1)%p ! First end of this element
          p2 => sups(i)%triads(2)%p ! Second end of this element
          if (associated(p1%supForce) .and. p1%isFEnode) then
             p1%supForce = -GetSectionalForces(sups(i),1)
          end if
          if (associated(p2%supForce) .and. p2%isFEnode) then
             p2%supForce = GetSectionalForces(sups(i),2)
          end if
       end if
    end do

  end subroutine updateBeamSectionForces


  !!============================================================================
  !> @brief Writes beam joint forces at current configuration to ASCII files.
  !>
  !> @param[in] t Current simulation time
  !> @param[out] ierr Error flag
  !>
  !> @details One file is written for each Triad with name `Triad_<userId>.asc`.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 3 Nov 2010

  subroutine writeBeamJointForces (t,ierr)

    use IdTypeModule       , only : getId, StrId
    use FileUtilitiesModule, only : findUnitNumber

    real(dp), intent(in)  :: t
    integer , intent(out) :: ierr

    !! Local variables
    logical           :: append
    integer           :: i, j, iunit, nLegs
    real(dp)          :: sf(6)
    character(len=32) :: fname

    !! --- Logic section ---

    ierr = 0
    if (.not. allocated(joints)) return

    iunit = findUnitNumber(10)

    do i = 1, size(joints)

       nLegs = size(joints(i)%legs)
       fname = 'Triad_'//trim(adjustl(StrId(joints(i)%triad%id%userId)))//'.asc'
       inquire(FILE=fname,EXIST=append)

       if (append) then
          !! Append to existing file
          open(iunit,FILE=fname,POSITION='APPEND')
       else
          !! Open a new file and write header
          open(iunit,FILE=fname)
          write(iunit,600,ADVANCE='NO') trim(getId(joints(i)%triad%id))
          do j = 1, nLegs
             write(iunit,601,ADVANCE='NO') trim(getId(joints(i)%legs(j)%p%id))
          end do
          write(iunit,610,ADVANCE='NO')
          do j = 1, nLegs
             write(iunit,611,ADVANCE='NO')
          end do
          write(iunit,*)
       end if

       !! Write beam end force for each leg
       write(iunit,620,ADVANCE='NO') t
       do j = 1, nLegs
          call getLegForce (joints(i)%legs(j)%p,joints(i)%triad,sf,ierr)
          write(iunit,621,ADVANCE='NO') sf(1),sf(5),sf(6)
       end do
       write(iunit,*)
       close(iunit)

    end do

600 format('# Sectional forces in Beam Joint',A / '#',T16)
601 format('Leg',A,T43)
610 format(/'# Time ')
611 format(T13,'Nx',T27,'My',T41,'Mz')
620 format(1P,E12.5)
621 format(1P,3E14.5)

  contains

    !> @brief Extracts the sectional force at a beam end.
    subroutine getLegForce (beam,triad,sf,ierr)
      type(SupElType), intent(in)  :: beam
      type(TriadType), pointer     :: triad
      real(dp)       , intent(out) :: sf(6)
      integer        , intent(out) :: ierr
      if (associated(beam%triads(1)%p,triad)) then
         sf = -GetSectionalForces(beam,1)
      else if (associated(beam%triads(2)%p,triad)) then
         sf = GetSectionalForces(beam,2)
      else
         ierr = ierr - 1
      end if
    end subroutine getLegForce

  end subroutine writeBeamJointForces

end module FiniteElementModule
