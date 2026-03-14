!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file windTurbineRoutinesModule.f90
!>
!> @brief Subroutines for wind turbine calculations.

!!==============================================================================
!> @brief Module with subroutines for wind turbine calculations.

module windTurbineRoutinesModule

  implicit none

  private

  public :: AeroInput, updateAeroForces, addInAeroForces, closeAeroDyn
  public :: getWindSpeed, getHubWindSpeed, getTipWindSpeed, getBladeDeflections


contains

  !!============================================================================
  !> @brief Sets up the information needed for, and initializes AeroDyn.
  !>
  !> @param env Environmental data
  !> @param turb Wind turbine configuration object
  !> @param[in] ADfile Name of the AeroDyn input file
  !> @param[in] NumBl Number of turbine blades
  !> @param[in] CompAero If .true., Aerodynamic forces are computed
  !> @param[in] CompNoise If .true., turbulent inflow noise is accounted for
  !> @param[in] SumPrint If .true., the resultant aerodynamic forces are printed
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 11 Oct 2010 / 2.0

  subroutine AeroInput (env,turb,ADFile,NumBl, &
       &                CompAero,CompNoise,SumPrint,ierr)

    use EnvironmentTypeModule, only : EnvironmentType
    use WindTurbineTypeModule, only : TurbineConfig
#if FT_HAS_AERODYN > 12
    use AeroDyn              , only : AD_InitOptions, AD_Init
    use versionModule        , only : progVer
#ifdef FT_DEBUG
    use fileUtilitiesModule  , only : getDBGfile
    use dbgUnitsModule       , only : dbgAD
#endif
    use reportErrorModule    , only : debugFileOnly_p, error_p
    use reportErrorModule    , only : getErrorFile, allocationError, reportError
    use NWTC_Library         , only : NWTC_Init
#else
    use reportErrorModule    , only : internalError
#endif

    !! Arguments
    type(EnvironmentType)      , intent(inout) :: env
    type(TurbineConfig), target, intent(inout) :: turb
    character(len=*)           , intent(in)    :: ADFile
    integer                    , intent(in)    :: NumBl
    logical                    , intent(in)    :: CompAero, CompNoise, SumPrint
    integer                    , intent(out)   :: ierr

    !! Local variables
#if FT_HAS_AERODYN > 12
    integer              :: nADbnods
    type(AD_InitOptions) :: ADoptions
#endif

    !! --- Logic section ---

    ierr = 0
#if FT_HAS_AERODYN > 12
#ifdef FT_DEBUG
    dbgAD = getDBGfile(11,'aerodynamics.dbg')
#endif

    !! Set up the AeroDyn parameters
    ADoptions%ADInputFile = ADFile
    ADoptions%OutRootName = 'fedem_aerodyn'
    ADoptions%WrSumFile   = SumPrint

    !! Define the AeroDyn turbine configuration
    call convertTurbineConfig (turb,turb%ADinterface,ierr)
    if (ierr < 0) then
       call reportError (debugFileOnly_p,'AeroInput')
       return
    else if (NumBl /= size(turb%ADinterface%Blade)) then
       ierr = -1
       call reportError (error_p,'Inconsistent number of blades', &
            &            addString='AeroInput')
       return
    end if

    !! Initialize AeroDyn, returning relative markers
    call NWTC_init ('fedem_solver',progVer,.false.,getErrorFile())
    turb%ADmarkers = AD_Init(ADoptions, turb%ADinterface, ierr)
    if (ierr /= 0) then
       ierr = -1
       call reportError (debugFileOnly_p,'AD_Init')
       call reportError (debugFileOnly_p,'AeroInput')
       return
    else if (.not. allocated(turb%ADmarkers%Blade)) then
       ierr = -1
       call reportError (error_p,'No AeroDyn blade nodes detected', &
            &            addString='AeroInput')
       return
    end if

    !! Get the number of blade nodes from the returned data structure
    nADbnods = size(turb%ADmarkers%Blade,1)
    if (nADbnods /= size(turb%node,2)) then
       ierr = -2
       call reportError (error_p,'AeroDyn blade node mismatch', &
            &            addString='AeroInput')
       return
    end if

    !! Allocate variables for AeroDyn forces
    allocate(turb%ADoptions%SetMulTabLoc(nADbnods,NumBl),STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('AeroInput: ADoptions%SetMulTabLoc')
       return
    end if

    turb%ADoptions%SetMulTabLoc  = .false.
    turb%ADoptions%LinearizeFlag = .false.

    !! Set Fedem variables based on AeroDyn inputs
    call Set_Fedem_Params (env,turb,CompAero,CompNoise,ierr)
    if (ierr < 0) call reportError (debugFileOnly_p,'AeroInput')
#else
    print *,' ** AeroInput dummy: ',env%rhoA,turb%id%userId,ADFile,NumBl, &
         &                          CompAero,CompNoise,SumPrint
    if (CompAero) ierr = internalError('AeroInput: AeroDyn is not available')
#endif

  end subroutine AeroInput


#if FT_HAS_AERODYN > 12
  !!============================================================================
  !> @brief Sets Fedem variables based on AeroDyn inputs.
  !>
  !> @param env Environmental data
  !> @param turbine Wind turbine configuration object
  !> @param[in] CompAero If .true., Aerodynamic forces are computed
  !> @param[in] CompNoise If .true., turbulent inflow noise is accounted for
  !> @param[out] ierr Error flag
  !>
  !> @details This subroutine is called at the start of the simulation to set up
  !> simulation variables for Fedem based on the AeroDyn parameters.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 11 Oct 2010

  subroutine Set_Fedem_Params (env,turbine,CompAero,CompNoise,ierr)

    use kindModule           , only : dp, pi_p
    use EnvironmentTypeModule, only : EnvironmentType
    use WindTurbineTypeModule, only : TurbineConfig
    use IdTypeModule         , only : getId
    use reportErrorModule    , only : allocationError
    use reportErrorModule    , only : getErrorFile, reportError
    use reportErrorModule    , only : error_p, warning_p, debugFileOnly_p
    use AeroDyn              , only : AD_GetConstant, ReKi
    use Element              , only : RELM
    use Blade                , only : DR

    !! Arguments
    type(EnvironmentType), intent(inout) :: env
    type(TurbineConfig)  , intent(inout) :: turbine
    logical              , intent(in)    :: CompAero, CompNoise
    integer              , intent(out)   :: ierr

    !! Local variables
    integer            :: i, lpu
    real(dp)           :: HHub
    character(len=128) :: MSG

    !! --- Logic section ---

    if (size(turbine%bladeElm) /= size(RELM)) then
       ierr = -1
       write(MSG,600) size(turbine%bladeElm), size(RELM)
600    format('Number of blades in Fedem:',I4,'  in AeroDyn:',I4)
       call reportError (error_p,'The number of blade elements in the '// &
            &            'AeroDyn input file has to match the Fedem model.', &
            &            MSG,'Update the AeroDyn input file and try again.', &
            &            addString='Set_Fedem_Params')
       return
    end if

    !! Override some blade geometry parameters read from the AeroDyn input file
    !! with what we derived from the Fedem turbine model. Only the blade Chord
    !! and Twist angle, as well as the airfoil file index parameters of the
    !! element data read will be used.
    do i = 1, size(RELM)
       RELM(i) = real(turbine%bladeElm(i)%R,ReKi)
       DR(i)   = real(turbine%bladeElm(i)%DR,ReKi)
    end do

    env%rhoa = AD_GetConstant('AirDensity',ierr)
    turbine%RefHH = AD_GetConstant('RefHt',i)
    if (ierr /= 0 .or. i /= 0) then
       ierr = -1
       call reportError (debugFileOnly_p,'Set_Fedem_Params')
       return
    end if

    HHub = turbine%PtfmRef + turbine%hub%ur(3,4) - turbine%tower%ur(3,4)

    lpu = getErrorFile()
    write(lpu,610) trim(getId(turbine%id)), &
         &         turbine%blade(1)%precone*180.0_dp/pi_p, HHub, turbine%TipRad

610 format(/5X,'>>> Computed parameters for Turbine',A &
         & /9X,'Precone angle =',F6.2,' [deg]' &
         & /9X,'Reference hub height =',1PE12.5 &
         & /9X,'Tip radius =',E12.5)

    if (CompNoise) then

       !! Let's compute the turbulence intensity and average wind speed
       !! for the turbulent inflow noise calculation:

       !InpPosition = (/ 0.0, 0.0, HHub /)
       !call Noise_CalcTI(0.0_ReKi, TMax, DT, InpPosition)
       !KinViscosity = AD_GetConstant('KinVisc',ierr)
       ierr = -1
       call reportError (error_p,'Noise_CalcTI is not implemented.', &
            &            'Set CompNoise to .false. and try again.', &
            &            addString='Set_Fedem_Params')

    end if
    if (CompAero) then

       !! Let's see if the hub-height in AeroDyn and Fedem are within 10%
       if (abs(HHub-turbine%RefHH) > 0.1*HHub) then
          write(MSG,620) turbine%RefHH, HHub
          call reportError (warning_p,MSG)
620       format('The AeroDyn input reference hub height (',1PE12.5,' ) and ' &
                 'the Fedem hub height (',E12.5,' ) differ by more than 10%')
       end if

    end if

  end subroutine Set_Fedem_Params


  !!============================================================================
  !> @brief Converts a Fedem wind turbine configuration to AeroDyn's format.
  !>
  !> @param[in] Fedem Wind turbine configuration object, Fedem format
  !> @param[out] AD Wind turbine configuration object, AeroDyn format
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 15 Oct 2010

  subroutine convertTurbineConfig (Fedem,AD,ierr)

    use WindTurbineTypeModule, only : TurbineConfig
    use AeroDyn              , only : AeroConfig, ReKi
    use reportErrorModule    , only : allocationError, reportError, error_p

    !! Arguments
    type(TurbineConfig), intent(in)  :: Fedem
    type(AeroConfig)   , intent(out) :: AD
    integer            , intent(out) :: ierr

    !! --- Logic section ---

    ierr = 0

    if (.not. associated(Fedem%tower)) then
       ierr = ierr - 1
       call reportError(error_p,'No tower triad specified')
    end if
    if (.not. associated(Fedem%nacelle)) then
       ierr = ierr - 1
       call reportError(error_p,'No nacelle triad specified')
    end if
    if (.not. associated(Fedem%hub)) then
       ierr = ierr - 1
       call reportError(error_p,'No hub triad specified')
    end if
    if (ierr < 0) return

    AD%BladeLength = real(Fedem%TipRad - Fedem%HubRad, ReKi)

    allocate(AD%Blade(size(Fedem%blade)),STAT=ierr)
    if (ierr == 0) then
       call updateTurbineConfig (Fedem,AD)
    else
       ierr = allocationError('convertTurbineConfig')
    end if

    !! Initialize the AeroDyn markers not yet in use
    call updateADMarker (AD%TailFin)
    call updateADMarker (AD%Substructure)
    call updateADMarker (AD%Foundation)

  end subroutine convertTurbineConfig


  !!============================================================================
  !> @brief Updates the AeroDyn wind turbine configuration.
  !>
  !> @param[in] Fedem Wind turbine configuration object, Fedem format
  !> @param[out] AD Wind turbine configuration object, AeroDyn format
  !>
  !> @note AeroDyn requires the inverse orientation matrix compared to Fedem
  !> (global-to-local vs. local-to-global). Therefore the transpose().
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 15 Oct 2010

  subroutine updateTurbineConfig (Fedem,AD)

    use WindTurbineTypeModule, only : TurbineConfig
    use AeroDyn              , only : AeroConfig, ReKi

    !! Arguments
    type(TurbineConfig), intent(in)    :: Fedem
    type(AeroConfig)   , intent(inout) :: AD

    !! Local variables
    integer :: i

    !! --- Logic section ---

    call updateADMarker (AD%Tower,Fedem%tower)
    call updateADMarker (AD%Nacelle,Fedem%nacelle)
    if (associated(Fedem%hubEl)) then
       call updateADMarker (AD%Hub,Fedem%hub,Fedem%hubEl)
    else
       call updateADMarker (AD%Hub,Fedem%hub)
       !! The Hub triad has a different coordinate system; x,y,z --> z,x,y
       AD%Hub%Orientation(1,:) = real(Fedem%hub%ur(:,3),ReKi)
       AD%Hub%Orientation(2,:) = real(Fedem%hub%ur(:,1),ReKi)
       AD%Hub%Orientation(3,:) = real(Fedem%hub%ur(:,2),ReKi)
    end if

    !! The rotor furl marker is new, and not yet present in the Fedem model.
    !! Use instead the Hub triad for positioning and the Shaft for orientation.
    call updateADMarker (AD%RotorFurl,Fedem%hub)
    AD%RotorFurl%Orientation(1,:) = real(Fedem%shaft%ur(:,3),ReKi)
    AD%RotorFurl%Orientation(2,:) = real(Fedem%shaft%ur(:,1),ReKi)
    AD%RotorFurl%Orientation(3,:) = real(Fedem%shaft%ur(:,2),ReKi)
    if (Fedem%shaft%nDOFs >= 6) then
       AD%RotorFurl%RotationVel   = real(Fedem%shaft%urd(4:6),ReKi)
    else
       AD%RotorFurl%RotationVel   = 0.0_ReKi
    end if

    do i = 1, size(Fedem%blade)
       call updateADMarker (AD%Blade(i),Fedem%blade(i)%triad)
    end do

  end subroutine updateTurbineConfig


  !!============================================================================
  !> @brief Updates an AeroDyn marker based on the corresponding Fedem triad.
  !>
  !> @param ADmarker A marker associated with the AeroDyn turbine configuration
  !> @param[in] triad Triad corresponding to the specified marker
  !> @param[in] sup Beam element defining the position/orientation of the marker
  !> @param[in] offset Eccentricity of the AeroDyn marker w.r.t the pitch axis
  !>
  !> @note AeroDyn requires the inverse orientation matrix compared to Fedem
  !> (global-to-local vs. local-to-global). Therefore the transpose().
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 15 Oct 2010

  subroutine updateADMarker (ADmarker,triad,sup,offset)

    use AeroDyn        , only : Marker, ReKi
    use TriadTypeModule, only : TriadType, dp
    use SupElTypeModule, only : SupElType

    !! Arguments
    type(Marker)   , intent(inout)        :: ADmarker
    type(TriadType), intent(in), optional :: triad
    type(SupElType), intent(in), optional :: sup
    real(dp)       , intent(in), optional :: offset(2)

    !! Local variables
    real(dp) :: omega(3), vl(3)

    !! --- Logic section ---

    if (present(sup)) then
       ADmarker%Position    = real(sup%supTr(:,4),ReKi)
       ADmarker%Orientation = real(transpose(sup%supTr(:,1:3)),ReKi)
       if (present(triad) .and. triad%nDOFs >= 6) then
          ADmarker%TranslationVel = real(triad%urd(1:3),ReKi)
          ADmarker%RotationVel    = real(triad%urd(4:6),ReKi)
       else
          ADmarker%TranslationVel = 0.0_ReKi
          ADmarker%RotationVel    = 0.0_ReKi
       end if
    else if (present(triad)) then
       if (present(offset)) then
          !! The Aerodyn marker have an eccentricity from the pitch axis
          vl = triad%ur(:,4) + matmul(triad%ur(:,2:3),offset)
       else
          vl = triad%ur(:,4)
       end if
       ADmarker%Position = real(vl,ReKi)
       ADmarker%Orientation = real(transpose(triad%ur(:,1:3)),ReKi)
       if (triad%nDOFs >= 6) then
          ADmarker%TranslationVel = real(triad%urd(1:3),ReKi)
          ADmarker%RotationVel    = real(triad%urd(4:6),ReKi)
          if (present(offset)) then
             if (any(abs(offset) > 1.0e-8_dp)) then
                !! Add velocity contribution due to eccentricity
                omega = matmul(triad%urd(4:6),triad%ur(:,1:3))
                vl(1) =  omega(2)*offset(2) - omega(3)*offset(1)
                vl(2) = -omega(1)*offset(2)
                vl(3) =  omega(1)*offset(1)
                ADmarker%TranslationVel = ADmarker%TranslationVel &
                     +                    real(matmul(triad%ur(:,1:3),vl),ReKi)
             end if
          end if
       end if
    else
       ADmarker%Position       = 0.0_ReKi
       ADmarker%Orientation    = 0.0_ReKi
       ADmarker%TranslationVel = 0.0_ReKi
       ADmarker%RotationVel    = 0.0_ReKi
    end if

  end subroutine updateADMarker

#endif

  !!============================================================================
  !> @brief Gets the aerodynamic forces at given time on the given wind turbine.
  !>
  !> @param[in] time Current physical time
  !> @param turb Wind turbine configuration object
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 1 Dec 2009

  subroutine updateAeroForces (time,turb,ierr)

    use WindTurbineTypeModule, only : TurbineConfig, dp
#if FT_HAS_AERODYN > 12
    use AeroDyn              , only : AD_CalculateLoads, AD_Output, ReKi
#ifdef FT_DEBUG
    use AeroDyn              , only : writeConfig, writeMarkers, writeLoads
    use dbgUnitsModule       , only : dbgAD
#endif
#endif
    use profilerModule       , only : startTimer, stopTimer, aed_p
    use reportErrorModule    , only : reportError, debugFileOnly_p

    !! Arguments
    real(dp)           , intent(in)    :: time
    type(TurbineConfig), intent(inout) :: turb
    integer            , intent(out)   :: ierr

    !! Local variables
    integer       :: i, j
    logical, save :: lFirst = .true.

    !! --- Logic section ---

    if (size(turb%node,1) < 1) return ! Aerodynamics has been turned off

    call startTimer (aed_p)

#if FT_HAS_AERODYN > 12
    !! Update the AeroDyn wind turbine markers
    call updateTurbineConfig (turb,turb%ADInterface)
    do i = 1, size(turb%node,1)
       do j = 1, size(turb%node,2)
          call updateADMarker (turb%ADmarkers%Blade(j,i),turb%node(i,j)%triad, &
               &               offset=turb%bladeElm(j)%ACoffset)
       end do
    end do
#ifdef FT_DEBUG
    write(dbgAD,"(//3X,77('=')/'   --- in updateAeroForces t =',1PE12.5/)") time
    call writeConfig (turb%ADInterface,dbgAD)
    call writeMarkers (turb%ADmarkers,dbgAD)
#endif

    !! Invoke the AeroDyn main driver
    turb%ADloads = AD_CalculateLoads(real(time,ReKi), &
         &                           turb%ADmarkers, &
         &                           turb%ADInterface, &
         &                           turb%ADoptions,ierr)
    if (ierr /= 0) goto 900

#ifdef FT_DEBUG
    call writeLoads (turb%ADloads,dbgAD)
#endif
    do i = 1, size(turb%node,1)
       do j = 1, size(turb%node,2)
          if (lFirst) turb%node(i,j)%aeroForcePrev = 0.0_dp
          !! AeroDyn delivers force/moment per unit length.
          !! In Fedem, we need the force resultants only.
          call scale (turb%ADloads%Blade(j,i)%Force,turb%bladeElm(j)%DR)
          call scale (turb%ADloads%Blade(j,i)%Moment,turb%bladeElm(j)%DR)
          !! Extract Normal/Tangent loads and Pitch moment (for plotting only)
          call getADforces (turb%ADInterface%Blade(i), &
               &            turb%ADmarkers%Blade(j,i), &
               &            turb%ADloads%Blade(j,i), turb%node(i,j)%aeroForce)
       end do
    end do
    call AD_Output ! Write out AeroDyn element data, if requested

#else
    ierr = -999
    do i = 1, size(turb%node,1)
       do j = 1, size(turb%node,2)
          turb%node(i,j)%aeroForce = time ! Dummy assignment
       end do
    end do
    goto 900
#endif

    !! Turn off the lFirst flag after we've looped through all the blades once
100 lFirst = .false.
    call stopTimer (aed_p)
    return

900 call reportError (debugFileOnly_p,'updateAeroForces')
    goto 100

#if FT_HAS_AERODYN > 12
  contains

    !> @brief Scales a real vector by given scalar.
    subroutine scale (vec,a)
      real(ReKi), intent(inout) :: vec(:)
      real(dp)  , intent(in)    :: a
      vec = vec*real(a,ReKi)
    end subroutine scale

    !> @brief Extracts the normal- and tangent loads and the pitch moment.
    subroutine getADforces (bl,nd,ADload,force)
      use AeroDyn , only : Marker, Load
      type(Marker), intent(in)  :: bl, nd
      type(Load)  , intent(in)  :: ADload
      real(dp)    , intent(out) :: force(3)
      real(dp) :: pitch, cp, sp
      !! This calculation of current pitch angle of the element
      !! is taken from the function AD_CalculateLoads in AeroDyn.f90
      pitch = -atan2(-dot_product(bl%Orientation(1,:),nd%Orientation(2,:)), &
           &          dot_product(bl%Orientation(1,:),nd%Orientation(1,:)))
      sp = sin(pitch)
      cp = cos(pitch)
      !! Undo the transformation done within AeroDyn such that we can output
      !! the same force components as FAST (and Fedem built with AeroDyn 12.58)
      force(1) = cp*ADload%Force(1) + sp*ADload%Force(2)
      force(2) = sp*ADload%Force(1) - cp*ADload%Force(2)
      force(3) = ADload%Moment(3)
    end subroutine getADforces
#endif

  end subroutine updateAeroForces


  !!============================================================================
  !> @brief Calculates system force vector contributions from a wind turbine.
  !>
  !> @param Q External force vector
  !> @param RF Reaction forces
  !> @param[in] turbine Wind turbine configuration object
  !> @param[in] sam Data for managing system matrix assembly
  !> @param ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 1 Dec 2009

  subroutine addInAeroForces (Q,RF,turbine,sam,ierr)

    use SamModule            , only : SamType, dp
    use WindTurbineTypeModule, only : TurbineConfig
#if FT_HAS_AERODYN > 12
    use TriadTypeModule      , only : transTriadToSys
#endif
    use AsmExtensionModule   , only : csAddNV
    use reportErrorModule    , only : reportError, debugFileOnly_p

    !! Arguments
    real(dp)           , intent(inout) :: Q(:), RF(:)
    type(TurbineConfig), intent(in)    :: turbine
    type(SamType)      , intent(in)    :: sam
    integer            , intent(inout) :: ierr

    !! Local variables
    integer  :: i, j, err, lerr
    real(dp) :: eV(6)

    !! --- Logic section ---

    lerr = ierr
    do i = 1, size(turbine%node,1)
       do j = 1, size(turbine%node,2)
#if FT_HAS_AERODYN > 12
          !! The AeroDyn forces now refer to the local triad coordinate system
          eV(1:3) = turbine%ADloads%Blade(j,i)%Force
          eV(4:6) = turbine%ADloads%Blade(j,i)%Moment
          eV(4) = eV(4) + turbine%bladeElm(j)%ACoffset(1)*eV(3) &
               &        - turbine%bladeElm(j)%ACoffset(2)*eV(2)
          eV(5) = eV(5) + turbine%bladeElm(j)%ACoffset(2)*eV(1)
          eV(6) = eV(6) - turbine%bladeElm(j)%ACoffset(1)*eV(1)
          call transTriadToSys (turbine%node(i,j)%triad,eV)
#endif
          call csAddNV (sam,-turbine%node(i,j)%triad%samNodNum,eV,Q,RF,err)
          if (err /= 0) ierr = ierr - 1
       end do
    end do

    if (ierr < lerr) call reportError (debugFileOnly_p,'addInAeroForces')

  end subroutine addInAeroForces


  !!============================================================================
  !> @brief Evaluates the wind speed at the hub center.
  !>
  !> @param[in] time Current physical time
  !> @param[in] turbine Wind turbine configuration object
  !> @param[out] dws Hub wind speed including disturbance
  !> @param[out] uws Undisturbed hub wind speed
  !> @param[out] ierr Error flag
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 9 Feb 2011

  subroutine getHubWindSpeed (time,turbine,dws,uws,ierr)

    use WindTurbineTypeModule, only : TurbineConfig, dp
#if FT_HAS_AERODYN > 12
    use AeroDyn , only : AD_GetUndisturbedWind, ReKi
    use AeroSubs, only : AD_WindVelocityWithDisturbance
#endif

    !! Arguments
    real(dp)           , intent(in)  :: time
    type(TurbineConfig), intent(in)  :: turbine
    real(dp)           , intent(out) :: dws(3), uws(3)
    integer            , intent(out) :: ierr

    !! --- Logic section ---

    ierr = 0
    dws = 0.0_dp
    uws = 0.0_dp
    if (size(turbine%node,1) < 1) return ! Aerodynamics has been turned off

#if FT_HAS_AERODYN > 12
    dws = AD_WindVelocityWithDisturbance(turbine%ADinterface%Hub%Position,ierr)
    uws = AD_GetUndisturbedWind(real(time,ReKi), &
         &                      turbine%ADinterface%Hub%Position,ierr)
#else
    print *,' ** getHubWindSpeed dummy: ',time
#endif

  end subroutine getHubWindSpeed


  !!============================================================================
  !> @brief Evaluates the wind speed at the tip of the given blade.
  !>
  !> @param[in] iBlade Index of the blade to evaluate the wind for
  !> @param[in] time Current physical time
  !> @param[in] turbine Wind turbine configuration object
  !> @param[out] dws Blade tip wind speed including disturbance
  !> @param[out] uws Undisturbed blade tip wind speed
  !> @param[out] ierr Error flag
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 9 Feb 2011

  subroutine getTipWindSpeed (iBlade,time,turbine,dws,uws,ierr)

    use WindTurbineTypeModule, only : TurbineConfig, dp

    !! Arguments
    integer            , intent(in)  :: iBlade
    real(dp)           , intent(in)  :: time
    type(TurbineConfig), intent(in)  :: turbine
    real(dp)           , intent(out) :: dws(3), uws(3)
    integer            , intent(out) :: ierr

    !! --- Logic section ---

    dws = 0.0_dp
    uws = 0.0_dp
    if (size(turbine%node,1) < 1) then
       ierr = 0 ! Aerodynamics has been turned off
    else if (iBlade > 0 .and. iBlade <= size(turbine%blade)) then
       call getWindSpeed (turbine%blade(iBlade)%tip%ur(:,4),time,dws,uws,ierr)
    else
       ierr = -1
    end if

  end subroutine getTipWindSpeed


  !!============================================================================
  !> @brief Evaluates the wind speed at a given location center.
  !>
  !> @param[in] pos Spatial point to evaluate the wind at
  !> @param[in] time Current physical time
  !> @param[out] dws Wind speed including disturbance
  !> @param[out] uws Undisturbed wind speed
  !> @param[out] ierr Error flag
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 5 Nov 2013

  subroutine getWindSpeed (pos,time,dws,uws,ierr)

    use kindModule, only : dp
#if FT_HAS_AERODYN > 12
    use AeroDyn , only : AD_GetUndisturbedWind, ReKi
    use AeroSubs, only : AD_WindVelocityWithDisturbance
#endif

    !! Arguments
    real(dp), intent(in)  :: pos(3), time
    real(dp), intent(out) :: dws(3), uws(3)
    integer , intent(out) :: ierr

    !! Local variables
#if FT_HAS_AERODYN > 12
    real(ReKi) :: XYZ(3)
#endif

    !! --- Logic section ---

    ierr = 0

#if FT_HAS_AERODYN > 12
    XYZ = real(pos,ReKi) ! cast to Real
    dws = AD_WindVelocityWithDisturbance(XYZ,ierr)
    uws = AD_GetUndisturbedWind(real(time,ReKi),XYZ,ierr)
#else
    dws = 0.0_dp
    uws = 0.0_dp
    print *,' ** getWindSpeed dummy: ',pos,time
#endif

  end subroutine getWindSpeed


  !!============================================================================
  !> @brief Evaluates the deflection along the blades.
  !>
  !> @param turbine Wind turbine configuration object
  !>
  !> @details The deflections are calculated relative to the blade root.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 9 Feb 2011

  subroutine getBladeDeflections (turbine)

    use WindTurbineTypeModule, only : TurbineConfig

    !! Arguments
    type(TurbineConfig), intent(inout) :: turbine

    !! Local variables
    integer :: i, j

    !! --- Logic section ---

    do i = 1, size(turbine%node,1)
       do j = 1, size(turbine%node,2)
          call calcBladeDeflection (turbine%node(i,j)%triad, &
               &                    turbine%blade(i)%cone, turbine%hub, &
               &                    turbine%bladeElm(j)%R)
       end do
       call calcBladeDeflection (turbine%blade(i)%tip, &
            &                    turbine%blade(i)%cone, turbine%hub, &
            &                    turbine%TipRad)
    end do

  contains

    !> @brief Evaluates the blade deflection at a given triad.
    subroutine calcBladeDeflection (triad,blade,hub,R)
      use TriadTypeModule, only : TriadType, transGlobToTriad, dp
      use rotationModule , only : deltaRot
      type(TriadType), intent(inout) :: triad
      type(TriadType), intent(in)    :: blade, hub
      real(dp)       , intent(in)    :: R
      triad%ur_def(1:3) = triad%ur(:,4) - hub%ur(:,4) - blade%ur(:,3)*R
      triad%ur_def(4:6) = deltaRot(blade%ur,triad%ur)
      call transGlobToTriad (blade,triad%ur_def(1:6))
      if (size(triad%ur_def) >= 12) then
         triad%ur_def(7:12) = triad%urd
         call transGlobToTriad (blade,triad%ur_def(7:12))
      end if
    end subroutine calcBladeDeflection

  end subroutine getBladeDeflections


  !!============================================================================
  !> @brief Terminates the aerodynamics module.
  !>
  !> @param[out] ierr Error flag
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 9 Feb 2011

  subroutine closeAeroDyn (ierr)

#if FT_HAS_AERODYN > 12
    use AeroDyn, only : AD_Terminate
#endif

    !! Arguments
    integer, intent(out) :: ierr

    !! --- Logic section ---

#if FT_HAS_AERODYN > 12
    call AD_Terminate (ierr)
#else
    ierr = 0
#endif

  end subroutine closeAeroDyn

end module windTurbineRoutinesModule
