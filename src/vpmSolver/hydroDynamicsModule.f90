!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file hydroDynamicsModule.f90
!> @brief Calculation of hydrodynamic element loads.

!!==============================================================================
!> @brief Module with subroutines for hydrodynamic load calculations.
!>
!> @details This module deals with the calculation of hydrodynamic loads
!> (buoyancy, inertia loads due to added mass, damping forces due to drag, etc.)
!> for beam- and superelements. It also manages the evaluation of the fluid
!> particle motions, which are needed in the hydrodynamic load calculations.
!>
!> Some of the subroutines and functions of this module are not documented.
!> You have to configure doxygen with the option ENABLED_SECTIONS = FULL_DOC
!> to extract detailed documentation of those subroutines and functions.

module HydroDynamicsModule

  use kindModule, only : sp, dp

  implicit none

  private

  real(dp), allocatable, save :: waterMotion(:,:) !< Fluid particle motions
  integer , allocatable, save :: calcWMotion(:)   !< Have calculated motions?

  real(dp), parameter :: eps_p = 1.0e-16_dp !< Zero tolerance

  !> @cond FULL_DOC
  real(sp), save :: tp(2)        !< Internal timing variables
  real(sp), save :: wavTime(2,2) !< Timing of wave calculation module (HW & SW)
  integer , save :: wavCall(2)   !< Number of wave calculation invokations

  logical, save :: useHWAFLS !< Flag for hardware wave evaluation
  logical, save :: useFNV    !< Flag for nonlinear wave force calculations
  !> @endcond

  !> @brief Updates hydrodynamics quantities after a time step is converged.
  interface UpdateAtConvergence
     module procedure UpdateHydroDynamicsAtConvergence
  end interface

  public :: initFluidMotions, getCalculatedFluidMotion, closeHydroDyn
  public :: InitiateHydroDynBodies, UpdateAtConvergence
  public :: getWaveElevation, getMorisonForces, getBuoyancyForces, getDragForces
  public :: getSeaState, diffractionCalc, getDiffractionForces


contains

  !> @cond FULL_DOC
  !> @brief Global-to-local transformation of a point vector.
  function glob2loc (Tlg,xg) result(xl)
    use manipMatrixModule, only : matmul34, invert34
    real(dp), intent(in) :: Tlg(3,4), xg(3)
    real(dp)             :: xl(3)
    xl = matmul34(invert34(Tlg),xg)
  end function glob2loc

  !> @brief Local-to-global transformation of a point vector.
  function loc2glob (Tlg,xl) result(xg)
    use manipMatrixModule, only : matmul34
    real(dp), intent(in) :: Tlg(3,4), xl(3)
    real(dp)             :: xg(3)
    xg = matmul34(Tlg,xl)
  end function loc2glob

  !> @brief Starts a timer.
  subroutine startTW ()
    real(sp), external   :: CLKSEC, CPUSEC
    tp(1) = CLKSEC(0.0_sp)
    tp(2) = CPUSEC(0.0_sp)
  end subroutine startTW

  !> @brief Stops a timer.
  subroutine stopTW (i)
    integer , intent(in) :: i
    real(sp), external   :: CLKSEC, CPUSEC
    wavTime(i,1) = wavTime(i,1) + CLKSEC(tp(1))
    wavTime(i,2) = wavTime(i,2) + CPUSEC(tp(2))
  end subroutine stopTW
  !> @endcond


  !!============================================================================
  !> @brief Initializes the fluid particle motion cache for all triads in water.
  !>
  !> @param[in] time Current simulation time
  !> @param[in] env Environmental data
  !> @param[in] triads All triads in the model
  !> @param[in] sups All superelements in the model
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 30 Jun 2010

  subroutine initFluidMotions (time,env,triads,sups,ierr)

    use EnvironmentTypeModule  , only : EnvironmentType
    use TriadTypeModule        , only : TriadType
    use SupElTypeModule        , only : SupElType
    use DiffractionModule      , only : NemohInit
    use explicitFunctionsModule, only : WAVE_SINUS_p
    use profilerModule         , only : startTimer, stopTimer, wav_p, hyd_p
    use reportErrorModule      , only : allocationError
    use reportErrorModule      , only : debugFileOnly_p, reportError
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_isTrue
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_intValue

    real(dp)             , intent(in)  :: time
    type(EnvironmentType), intent(in)  :: env
    type(TriadType)      , intent(in)  :: triads(:)
    type(SupElType)      , intent(in)  :: sups(:)
    integer              , intent(out) :: ierr

    !! Local variables
    integer  :: i, nnod
    real(dp) :: X(3)

    !! --- Logic section ---

    call startTimer (hyd_p)

    useHWAFLS = .false.
    useFNV    = .false.
    if (associated(env%waveFunc)) then
       if (env%waveFunc%type == WAVE_SINUS_p) then
          useHWAFLS = ffa_cmdlinearg_isTrue('HWAFLS')
          useFNV    = ffa_cmdlinearg_intValue('FNV') > 0
       end if
    end if

    if (.not. allocated(calcWMotion)) then
       nnod = 0
       do i = 1, size(triads)
          if (triads(i)%inFluid > nnod) nnod = triads(i)%inFluid
       end do
       allocate(calcWMotion(nnod),waterMotion(11,nnod),STAT=ierr)
       if (ierr /= 0) then
          ierr = allocationError('initFluidMotion')
          return
       end if
       calcWMotion = 0
       waterMotion = 0.0_dp
       wavTime = 0.0_sp
       wavCall = 0
       if (ffa_cmdlinearg_intValue('diffraction') > 0) then
          call NemohInit (sups,ierr)
          if (ierr /= 0) call reportError (debugFileOnly_p,'initFluidMotions')
       end if
       if (useHWAFLS) then
          i = env%waveFunc%intParameters(3)
          call HWAFLS_init (nnod,i,env%seaDepth,env%waveFunc%realParameters)
       end if
    else
       calcWMotion = 0
       waterMotion(1:10,:) = 0.0_dp ! retain the 11th column
       ierr = 0
    end if

    if (useHWAFLS .and. time >= 0.0_dp) then
       call startTimer (wav_p)
       do i = 1, size(triads)
          nnod = triads(i)%inFluid
          if (nnod > 0) then
             X = glob2loc(env%Tsea,triads(i)%ur(:,4))
             call HWAFLS_setNode (nnod,X(1),X(3))
          end if
       end do
       call startTW ()
       call HWAFLS_evaluate (time)
       call stopTW (1)
       call stopTimer (wav_p)
    end if

    call stopTimer (hyd_p)

  end subroutine initFluidMotions


  !!============================================================================
  !> @brief Returns the latest calculated fluid particle motion at a triad.
  !>
  !> @param[in] triad The triad object to get particle motion data for
  !> @param[out] elev Sea surface elevation at the triad location
  !> @param[out] fvel Water particle velocity at the triad location
  !> @param[out] facc Water particle acceleration at the triad location
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 29 Nov 2010

  function getCalculatedFluidMotion (triad,elev,fvel,facc)

    use TriadTypeModule, only : TriadType

    type(TriadType)  , intent(in)  :: triad
    real(dp),optional, intent(out) :: elev, fvel(3), facc(3)
    logical                        :: getCalculatedFluidMotion

    !! --- Logic section ---

    getCalculatedFluidMotion = triad%inFluid > 0

    if (allocated(waterMotion) .and. getCalculatedFluidMotion) then
       if (present(fvel)) fvel = waterMotion(4:6,triad%inFluid)
       if (present(facc)) facc = waterMotion(7:9,triad%inFluid)
       if (present(elev)) elev = waterMotion(10,triad%inFluid)
    end if

  end function getCalculatedFluidMotion


  !!============================================================================
  !> @brief Returns the water surface normal vector in local coordinate system.
  !>
  !> @param[out] normal Unit normal vector for the water surface
  !> @param[in] gravity Global gravitation vector
  !> @param[in] Tlg Local to global transformation matrix (optional)
  !> @return The gravitation constant
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 27 Mar 2009

  function waterSurfaceNormal (normal,gravity,Tlg) result(g)

    use kindModule       , only : epsDiv0_p
    use reportErrorModule, only : reportError, error_p

    real(dp), intent(out)          :: normal(3)
    real(dp), intent(in)           :: gravity(3)
    real(dp), intent(in), optional :: Tlg(:,:)

    !! Local variables
    real(dp) :: g

    !! --- Logic section ---

    g = gravity(1)*gravity(1) + gravity(2)*gravity(2) + gravity(3)*gravity(3)
    if (g < epsDiv0_p) then
       g = 0.0_dp
       call reportError (error_p,'Can not do hydrodynamics in zero gravity')
       return
    end if

    g = sqrt(g)
    if (present(Tlg)) then
       normal = matmul(-gravity/g,Tlg(1:3,1:3))
    else
       normal = -gravity/g
    end if

  end function waterSurfaceNormal


  !!============================================================================
  !> @brief Evaluates the wave profile, velocity and acceleration at a point.
  !>
  !> @param[in] waveFunc The wave function to evaluate for
  !> @param[in] waveTheory Flag indicating which wave theory to use
  !> @param[in] Twave Wave coordinate system (axis directions and origin)
  !> @param[in] g Gravitation constant
  !> @param[in] depth Water depth (zero means infinite)
  !> @param[in] Xg Global coordinates of the evaluation point
  !> @param[in] time Current simulation time
  !> @param[in] scale Scaling factor
  !> @param[out] wave Wave kinematics quantities at the evaluation point.
  !> wave(:,1) = Projection of point coordinates onto the current sea surface.
  !> wave(:,2) = Water particle velocity.
  !> wave(:,3) = Water particle acceleration.
  !> @param[out] dynp Dynamic pressure (optional)
  !> @param stat Status flag (negative on error exit).
  !> 1 = The point is below or on the water surface.
  !> 2 = The point is above the water surface (no kinematics).
  !!
  !> @details The calculated wave kinematics quantities are referring to
  !> the global coordinate system, unless when @a stat is negative on entry.
  !> In the latter case the wave coordinate system is used.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 25 Mar 2009

  subroutine evaluateWave (waveFunc,waveTheory,Twave,g,depth,Xg,time, &
       &                   scale,wave,dynp,stat)

    use FunctionTypeModule     , only : FunctionType
    use explicitFunctionsModule, only : WAVE_SINUS_p, WAVE_STREAMLINE_p
    use explicitFunctionsModule, only : WAVE_STOKES5_p, WAVE_EMBEDDED_p
    use explicitFunctionsModule, only : SINUSOIDAL_p, COMPL_SINUS_p
    use explicitFunctionsModule, only : DELAYED_COMPL_SINUS_p, USER_DEFINED_p
    use waveFunctionsModule    , only : stokes2Wave, stokes5Wave, streamWave
    use waveFunctionsModule    , only : deepWave, deepWaves, embeddedWave
    use waveFunctionsModule    , only : finiteDepthWave, finiteDepthWaves, twoPi
    use waveFunctionsModule    , only : userDefinedWave
    use dbgUnitsModule         , only : dbgHD ! Note that this use statement
    !! has always to be there (not only when FT_DEBUG is defined), or else the
    !! dependendency generator in cmake fails on Linux when building debug mode
    use reportErrorModule      , only : reportError, error_p
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_getbool

    type(FunctionType), intent(in)    :: waveFunc
    integer           , intent(in)    :: waveTheory
    real(dp)          , intent(in)    :: Twave(3,4), g,depth, Xg(3), time, scale
    real(dp)          , intent(out)   :: wave(3,3)
    real(dp), optional, intent(out)   :: dynp
    integer           , intent(inout) :: stat

    !! Local variables
    integer           :: n, m, nDir
    logical           :: noWStretch, atSurface
    real(dp)          :: X(4)
    real(dp), target  :: sFunc(6)
    real(dp), pointer :: rFunc(:)

    !! --- Logic section ---

    n = 0
    m = 3
    nDir = 1
    wave = 0.0_dp
    select case (waveFunc%type)
    case (WAVE_SINUS_p) ! Irregular waves
       n    = waveFunc%intParameters(3)
       nDir = waveFunc%intParameters(4)
       m    = size(waveFunc%realParameters)/(n*nDir)
       rFunc => waveFunc%realParameters
    case (SINUSOIDAL_p) ! Regular waves
       n = 1
       sFunc(1) =  waveFunc%realParameters(3)
       sFunc(2) =  waveFunc%realParameters(1) * twoPi
       sFunc(3) = -waveFunc%realParameters(2) * twoPi
       rFunc => sFunc(1:3)
    case (COMPL_SINUS_p, DELAYED_COMPL_SINUS_p) ! Two-component sine wave
       n = 2
       sFunc(1) =  waveFunc%realParameters(5)
       sFunc(2) =  waveFunc%realParameters(1) * twoPi
       sFunc(3) = -waveFunc%realParameters(3) * twoPi
       sFunc(4) =  waveFunc%realParameters(6)
       sFunc(5) =  waveFunc%realParameters(2) * twoPi
       sFunc(6) = -waveFunc%realParameters(4) * twoPi
       rFunc => sFunc
    case (WAVE_STOKES5_p, WAVE_STREAMLINE_p) ! Nonlinear wave functions
       n = (size(waveFunc%realParameters)-3)/2
    case (WAVE_EMBEDDED_p) ! Irregular waves with embedded streamline wave
       n = waveFunc%intParameters(3)
       m = size(waveFunc%realParameters)/n
       rFunc => waveFunc%realParameters
    case (USER_DEFINED_p) ! User-defined wave function (plug-in)
       n = max(1,size(waveFunc%realParameters))
       nDir = waveFunc%nArg
    end select
    if (n < 1 .or. nDir < 1) then
       stat = -1
       call reportError (error_p,'Invalid wave profile function', &
            &            addString='evaluateWave')
       return
    else if (waveFunc%type == USER_DEFINED_p .and. nDir /= 4) then
       stat = -nDir
       call reportError (error_p,'Invalid user-defined wave profile function', &
            &            'The number of function arguments should be 4', &
            &            addString='evaluateWave')
       return
    end if

    X(1:3) = glob2loc(Twave,Xg) ! Evaluation point in the wave coordinate system
    X(4) = time

    atSurface = abs(stat) == 2
    call ffa_cmdlinearg_getbool ('noWheelerStretching',noWStretch)

    if (waveFunc%type == WAVE_EMBEDDED_p) then

       !! Nonlinear streamline wave theory embedded in an irreagular sea state
       call embeddedWave (waveFunc%intParameters, waveFunc%realParameters, &
            &             g, depth, X(1), X(3), time, atSurface, noWStretch, &
            &             wave(3,1), wave(1,2), wave(3,2), wave(1,3), wave(3,3))

    else if (waveFunc%type == WAVE_STREAMLINE_p) then

       !! Nonlinear streamline wave theory
       call streamWave (waveFunc%realParameters, depth, &
            &           X(1), X(3), time, atSurface, wave(3,1), &
            &           wave(1,2), wave(3,2), wave(1,3), wave(3,3))

    else if (waveFunc%type == WAVE_STOKES5_p) then

       !! 5th order Stokes wave theory
       call stokes5Wave (waveFunc%realParameters, depth, &
            &            X(1), X(3), time, atSurface, wave(3,1), &
            &            wave(1,2), wave(3,2), wave(1,3), wave(3,3))

    else if (waveFunc%type == USER_DEFINED_p) then

       !! User-defined wave function (plug-in)
       call userDefinedWave (waveFunc%id%baseId, &
            &                waveFunc%intParameters, &
            &                waveFunc%realParameters, g, depth, &
            &                X, atSurface, wave(3,1), wave(:,2), wave(:,3), m)
       if (m < 0) then
          stat = m
          return
       end if

    else if (waveTheory == 2) then

       !! 2nd order Stokes wave theory
       call stokes2Wave (reshape(rFunc,(/m,n/)), g, depth, &
            &            X(1), X(3), time, atSurface, noWStretch, wave(3,1), &
            &            wave(1,2), wave(3,2), wave(1,3), wave(3,3))

    else if (depth > 0.0_dp) then

       !! Airy (linear) wave theory, finite depth
       if (nDir > 1) then
          call finiteDepthWaves (reshape(rFunc,(/m,n,nDir/)), g, depth, &
               &                 X(1),X(2),X(3), time, atSurface, noWStretch, &
               &                 wave(3,1), wave(1,2), wave(2,2), wave(3,2), &
               &                 wave(1,3), wave(2,3), wave(3,3), dynp)
       else !! No wave spreading
          call finiteDepthWave (reshape(rFunc,(/m,n/)), g, depth, &
               &                X(1), X(3), time, atSurface, noWStretch, &
               &                wave(3,1), wave(1,2), wave(3,2), &
               &                wave(1,3), wave(3,3), dynp)
       end if

    else

       !! Airy (linear) wave theory, deep water
       if (nDir > 1) then
          call deepWaves (reshape(rFunc,(/m,n,nDir/)), g, &
               &          X(1), X(2), X(3), time, atSurface, noWStretch, &
               &          wave(3,1), wave(1,2), wave(2,2), wave(3,2), &
               &          wave(1,3), wave(2,3), wave(3,3), dynp)
       else
          call deepWave (reshape(rFunc,(/m,n/)), g, X(1), X(3), time, &
               &         atSurface, noWStretch, wave(3,1), &
               &         wave(1,2), wave(3,2), wave(1,3), wave(3,3), dynp)
       end if

    end if
#ifdef FT_DEBUG
    if (nDir > 1) then
       write(dbgHD,"(7X,'x, y, z, eta, v, a =',1P10E13.5)",advance='NO') &
            X(1:3), wave(3,1), wave(:,2), wave(:,3)
    else
       write(dbgHD,"(7X,'x, z, eta, v, a =',1P7E13.5)",advance='NO') &
            X(1), X(3), wave(3,1), wave(1,2), wave(3,2), wave(1,3), wave(3,3)
    end if
    if (present(dynp)) then
       write(dbgHD,"(', p =',1PE13.5)") dynp
    else
       write(dbgHD,*)
    end if
#else
    m = dbgHD ! Dummy statement
#endif

    m = stat
    if (atSurface .or. X(3) <= wave(3,1)) then
       stat = 1 ! This point is on or below the water surface
    else
       stat = 2 ! This point is above the water surface
    end if
    if (m < 0) return

    !! Scale and transform to global coordinate system
    X(3) = wave(3,1)*scale
    wave(:,1) = loc2glob(Twave,X(1:3))
    if (stat == 2) then
       wave(:,2:3) = wave(:,2:3)*scale
    else
       wave(:,2) = matmul(Twave(:,1:3),wave(:,2)*scale)
       wave(:,3) = matmul(Twave(:,1:3),wave(:,3)*scale)
    end if

  end subroutine evaluateWave


  !!============================================================================
  !> @brief Evaluates the sea current velocity at the given point and time.
  !>
  !> @param[in] currFunc The sea current velocity function to evaluate for
  !> @param[in] dirFunc The sea current direction function
  !> @param[in] Tsea Coordinate system for the sea current
  !> @param[in] Xg Global coordinates of the evaluation point
  !> @param[in] time Current simulation time
  !> @param[in] scale Scaling factor
  !> @param[out] cvel Sea current velocity in global axis directions
  !> @param stat Status flag (negative on error exit).
  !> 1 = The point is below or on the water surface.
  !> 2 = The point is above the water surface (no kinematics).
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 27 Oct 2009

  subroutine evaluateCurrent (currFunc,dirFunc,Tsea,Xg,time,scale,cvel,stat)

    use FunctionTypeModule, only : FunctionType, FunctionValue
#ifdef FT_DEBUG
    use dbgUnitsModule    , only : dbgHD
#endif
    use reportErrorModule , only : reportError, debugFileOnly_p

    type(FunctionType), pointer       :: currFunc, dirFunc
    real(dp)          , intent(in)    :: Tsea(3,4), Xg(3), time, scale
    real(dp)          , intent(out)   :: cvel(3)
    integer           , intent(inout) :: stat

    !! Local variables
    real(dp) :: Xloc(3), curr, cdir

    !! --- Logic section ---

    Xloc = glob2loc(Tsea,Xg)
    if (stat == 0) then ! Check if this point is in the water
       cvel = 0.0_dp
       if (Xloc(3) > 0.0_dp) then
          stat = 2 ! This point is above the (calm) sea surface
       else if (.not. associated(currFunc)) then
          stat = 1 ! This point is in still water (no waves, no current)
       end if
       if (stat > 0) return
    else if (.not. associated(currFunc)) then
       cvel = 0.0_dp
       return
    end if

    stat = 0
    curr = scale*FunctionValue(currFunc,getArg(currFunc%nArg),stat)
    if (associated(dirFunc)) then
       cdir = FunctionValue(dirFunc,getArg(dirFunc%nArg),stat)
    else
       cdir = 0.0_dp
    end if
    if (stat < 0) then
       call reportError (debugFileOnly_p,'evaluateCurrent')
       return
    end if

    cvel(1) = curr*cos(cdir)
    cvel(2) = curr*sin(cdir)
    cvel(3) = 0.0_dp
#ifdef FT_DEBUG
    write(dbgHD,"(7X,'z, current(dir,x,y) =',1P4E13.5)") Xloc(3),cdir,cvel(1:2)
#endif
    cvel = matmul(Tsea(:,1:3),cvel)
    stat = 1

  contains

    !> @brief Returns the function arguments depending on its dimensionality.
    function getArg (numArg) result(x)
      integer, intent(in) :: numArg
      real(dp)            :: x(numArg)
      if (numArg < 1) return
      if (numArg < 3) then
         x(1) = Xloc(3)
      else
         x(1:3) = Xloc
      end if
      if (numArg == 2 .or. numArg >= 4) then
         x(numArg) = time
      else if (numArg > 4) then
         x(5:numArg) = 0.0_dp
      end if
    end function getArg

  end subroutine evaluateCurrent


  !!============================================================================
  !> @brief Evaluates the sea state at the given point and time.
  !>
  !> @param[in] env Environmental data
  !> @param[in] g Gravitation constant
  !> @param[in] time Current simulation time
  !> @param[in] istep Time increment counter
  !> @param[in] inod Nodal number of sea kinematics evaluation point
  !> @param[in] x Global coordinates of sea kinematics evaluation point
  !> @param[out] waterMotion Sea kinematics state at current point.
  !> waterMotion(1:3) = Projection of point x onto the current sea surface.
  !> waterMotion(4:6) = Particle velocity at point x in global coordinates.
  !> waterMotion(7:9) = Particle acceleration at point x in global coordinates.
  !> @param stat Status flag (negative on error exit)
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 9 Sep 2011

  subroutine evaluateSea (env,g,time,istep,inod,x,waterMotion,stat)

    use EnvironmentTypeModule, only : EnvironmentType
    use FNVwaveForceModule   , only : getWaveKinematics
    use profilerModule       , only : startTimer, stopTimer, wav_p
#ifdef FT_DEBUG
    use dbgUnitsModule       , only : dbgHD
#endif
    use reportErrorModule    , only : debugFileOnly_p, reportError

    type(EnvironmentType), intent(in)    :: env
    real(dp)             , intent(in)    :: g, time, x(3)
    integer              , intent(in)    :: istep, inod
    real(dp)             , intent(out)   :: waterMotion(:)
    integer              , intent(inout) :: stat

    !! Local variables
    integer  :: istat
    real(dp) :: cvel(3)

    !! --- Logic section ---

    istat = stat
    if (associated(env%waveFunc)) then
       call startTimer (wav_p)
       if (useFNV .and. istat == 0) then
          !! Use pre-evaluated wave kinematics
          call getWaveKinematics (istep,inod,waterMotion,stat)
       else if (useHWAFLS .and. istat == 0) then
          !! Evaluate wave kinematics in hardware
          call startTW ()
          call HWAFLS_getNode (inod,waterMotion(1))
          call stopTW (1)
          wavCall(1) = wavCall(1) + 1
          stat = 3
#ifdef FT_DEBUG
          write(dbgHD,"(7X,'x, z, eta, v, a =',1P7E13.5)") X(1), X(3), &
               &      waterMotion(3:4), waterMotion(6:7), waterMotion(9)
#endif
       else
          !! Evaluate wave kinematics in software (possibly at updated location)
          call startTW ()
          call evaluateWave (env%waveFunc,env%waveTheory,env%Tsea,g, &
               &             env%seaDepth,x,time,env%seaScale,waterMotion, &
               &             stat=stat)
          call stopTW (2)
          wavCall(2) = wavCall(2) + 1
       end if
       call stopTimer (wav_p)
       if (stat == 3) then
          cvel = glob2loc(env%Tsea,x)
          if (cvel(3) <= waterMotion(3)) then
             stat = 1 ! Point is in the sea
          else
             stat = 2 ! Point is above the sea
          end if
          cvel(3) = waterMotion(3)
          waterMotion(1:3) = loc2glob(env%Tsea,cvel)*env%seaScale
          waterMotion(4:9) = waterMotion(4:9)*env%seaScale
          if (stat == 1) then
             waterMotion(4:6) = matmul(env%Tsea(:,1:3),waterMotion(4:6))
             waterMotion(7:9) = matmul(env%Tsea(:,1:3),waterMotion(7:9))
          end if
       end if
    else
       waterMotion = 0.0_dp ! calm sea
    end if

    if (stat == 0 .or. stat == 1 .or. (istat == 2 .and. stat > 0)) then
       if (istat == 2) stat = 2
       call evaluateCurrent (env%currFunc,env%cDirFunc,env%Tsea, &
            &                x,time,env%currScale,cvel,stat)
       waterMotion(4:6) = waterMotion(4:6) + cVel
    end if

    if (stat < 0) call reportError (debugFileOnly_p,'evaluateSea')

  end subroutine evaluateSea


  !!============================================================================
  !> @brief Returns the wave height at the given point and time.
  !>
  !> @param[in] env Environmental data
  !> @param[in] Xg Global coordinates of sea kinematics evaluation point
  !> @param[in] time Current simulation time
  !> @param[out] stat Status flag (negative on error exit)
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 8 Feb 2013

  function getWaveElevation (env,Xg,time,stat)

    use EnvironmentTypeModule, only : EnvironmentType
#ifdef FT_DEBUG
    use dbgUnitsModule       , only : dbgHD
#endif

    type(EnvironmentType), intent(in)  :: env
    real(dp)             , intent(in)  :: time, Xg(3)
    integer              , intent(out) :: stat

    !! Local variables
    real(dp) :: Xw(3), wave(9), getWaveElevation

    !! --- Logic section ---

#ifdef FT_DEBUG
    write(dbgHD,"(/'   --- in getWaveElevation, X, time :',1P4E13.5)") Xg,time
#endif

    if (associated(env%waveFunc)) then
       !! Modify the point coordinate to ensure we are above the current
       !! sea leavel (to avoid uneccesary kinematics evaluations)
       wave(1:3) = glob2loc(env%Tsea,Xg);
       wave(3) = 100.0_dp
       Xw = loc2glob(env%Tsea,wave(1:3))
       stat = -1 ! Get the wave height, in sea coordinate system
       call evaluateWave (env%waveFunc,env%waveTheory,env%Tsea, &
            &             sqrt(sum(env%gravity*env%gravity)), &
            &             env%seaDepth,Xw,time,env%seaScale,wave,stat=stat)
       getWaveElevation = wave(3)
    else
       stat = 0
       getWaveElevation = 0.0_dp
    end if

  end function getWaveElevation


  !!============================================================================
  !> @brief Evaluates the sea state at the given point and time.
  !>
  !> @param[in] env Environmental data
  !> @param[in] time Current simulation time
  !> @param[in] Xg Global coordinates of sea kinematics evaluation point
  !> @param[out] sea Sea kinematics state at current point.
  !> sea(1:3) = Particle velocity at point @a Xg in global coordinates.
  !> sea(4:6) = Particle acceleration at point @a Xg in global coordinates.
  !> sea(7) = Dynamic pressure at point @a Xg in global coordinates.
  !> @param[out] stat Status flag (negative on error exit)
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 29 Jun 2015

  subroutine getSeaState (env,time,Xg,sea,stat)

    use EnvironmentTypeModule, only : EnvironmentType
    use profilerModule       , only : startTimer, stopTimer, wav_p
#ifdef FT_DEBUG
    use dbgUnitsModule       , only : dbgHD
#endif
    use reportErrorModule    , only : debugFileOnly_p, reportError

    type(EnvironmentType), intent(in)  :: env
    real(dp)             , intent(in)  :: time, Xg(3)
    real(dp)             , intent(out) :: sea(:)
    integer              , intent(out) :: stat

    !! Local variables
    real(dp) :: wave(9), dynp, cvel(3)

    !! --- Logic section ---

#ifdef FT_DEBUG
    write(dbgHD,"(/'   --- in getSeaState, X, time :',1P4E13.5)") Xg,time
#endif

    stat = 0
    if (associated(env%waveFunc)) then
       call startTimer (wav_p)
       call startTW ()
       call evaluateWave (env%waveFunc,env%waveTheory,env%Tsea, &
            &             sqrt(sum(env%gravity*env%gravity)), &
            &             env%seaDepth,Xg,time,env%seaScale,wave,dynp,stat)
       call stopTW (2)
       wavCall(2) = wavCall(2) + 1
       call stopTimer (wav_p)
       sea(1:6) = wave(4:9)
       sea(7) = env%rhow*dynp
    else
       sea = 0.0_dp ! calm sea
    end if

    if (stat == 0 .or. stat == 1) then
       call evaluateCurrent (env%currFunc,env%cDirFunc,env%Tsea, &
            &                Xg,time,env%currScale,cvel,stat)
       sea(1:3) = sea(1:3) + cVel
    end if

    if (stat < 0) call reportError (debugFileOnly_p,'getSeaState')

  end subroutine getSeaState


  !!============================================================================
  !> @brief Initiates all hydrodynamic bodies in the model.
  !>
  !> @param sups All superelements in the model
  !> @param elms All user-defined elements in the model
  !> @param[in] env Environmental data
  !> @param[in] restart If .true., this is a restart simulation
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 12 Aug 2008

  subroutine InitiateHydroDynBodies (sups,elms,env,restart,ierr)

    use SupElTypeModule       , only : SupElType
    use UserdefElTypeModule   , only : UserdefElType
    use EnvironmentTypeModule , only : EnvironmentType
    use DiffractionModule     , only : NemohMesh
    use manipMatrixModule     , only : matmul34, invert34
    use reportErrorModule     , only : reportError, debugFileOnly_p
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getint

    type(SupElType)      , intent(inout) :: sups(:)
    type(UserdefElType)  , intent(inout) :: elms(:)
    type(EnvironmentType), intent(in)    :: env
    logical              , intent(in)    :: restart
    integer              , intent(out)   :: ierr

    !! Local variables
    integer :: i, lerr, ndpanel

    !! --- Logic section ---

    call ffa_cmdlinearg_getint ('diffraction',ndpanel)

    ierr = 0
    do i = 1, size(sups)
       if (associated(sups(i)%hydyn)) then

          if (restart) then
             call updateHydroDynBody (sups(i)%hydyn,sups(i)%triads, &
                  &                   sups(i)%supTr,env%seaLevel, &
                  &                   env%gravity,ierr=lerr)
             if (lerr < 0) ierr = ierr - 1
          end if

          call updateAtConvergence (sups(i)%hydyn,sups(i)%supTr,lerr)
          if (lerr < 0) ierr = ierr - 1

          if (sups(i)%hydyn%bodyIndex >= 0) then
             !! Prepare for diffraction analysis
             call NemohMesh ('NemohDB',sups(i)%id,sups(i)%hydyn%bodyIndex, &
                  &          matmul34(invert34(env%Tsea),sups(i)%supTr), &
                  &          sups(i)%posMassCenter,nfobj=ndpanel,ierr=lerr)
             if (lerr < 0) ierr = ierr - 1
          end if

       end if
    end do
    do i = 1, size(elms)
       if (associated(elms(i)%hydyn)) then

          if (restart) then
             call updateHydroDynBody (elms(i)%hydyn,elms(i)%triads, &
                  &                   elms(i)%Tlg,env%seaLevel, &
                  &                   env%gravity,ierr=lerr)
             if (lerr < 0) ierr = ierr - 1
          end if

          call updateAtConvergence (elms(i)%hydyn,elms(i)%Tlg,lerr)
          if (lerr < 0) ierr = ierr - 1

       end if
    end do

    if (ierr < 0) call reportError (debugFileOnly_p,'InitiateHydroDynBodies')

  end subroutine InitiateHydroDynBodies


  !!============================================================================
  !> @brief Updates the hydrodynamic body quantities.
  !>
  !> @param hydyn Data for hydrodynamic force calculation
  !> @param[in] triads Triads on the superelement to update hydrodynamics for
  !> @param[in] supTr Superelement position matrix
  !> @param[in] sLev Current sea level
  !> @param[in] gravity Global gravitation vector
  !> @param[out] g Gravitation constant
  !> @param[out] wb Buoyancy weight factors for two-noded beams
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 12 Aug 2008

  subroutine updateHydroDynBody (hydyn,triads,supTr,sLev,gravity,g,wb,ierr)

    use SupElTypeModule        , only : HydroDynType
    use TriadTypeModule        , only : TriadPtrType
    use reportErrorModule      , only : reportError, debugFileOnly_p
    use FFaBodyHandlerInterface, only : ffa_partial_volume

    type(HydroDynType), intent(inout) :: hydyn
    type(TriadPtrType), intent(in)    :: triads(:)
    real(dp)          , intent(in)    :: supTr(:,:), sLev, gravity(3)
    real(dp),optional , intent(out)   :: g, wb(2)
    integer           , intent(out)   :: ierr

    !! Local variables
    real(dp) :: z0, Lb

    !! --- Logic section ---

    ierr = -1
    Lb = waterSurfaceNormal(hydyn%wn,gravity,supTr)
    if (Lb <= 0.0_dp) then
       goto 100
    else if (present(g)) then
       g = Lb
    end if

    !! Find the local height of the water surface w.r.t. link origin
    z0 = sLev + dot_product(gravity/Lb,supTr(:,4))

    if (hydyn%bodyIndex >= 0) then
       !! Compute buoyancy volume and center, and area of the waterline plane,
       !! for a (possibly) partly submerged body
       call ffa_partial_volume (hydyn%bodyIndex, hydyn%wn, z0, &
            &                   hydyn%Vb(1), hydyn%As(1), &
            &                   hydyn%C0b(:,1), hydyn%C0As(:,1), ierr)
       if (present(wb) .and. size(triads) == 2) then
          wb = 0.5_dp !TODO: Compute wb from C0b and the nodal coordinates
       end if
    else if (size(triads) == 2) then
       !! Compute buoyancy volume and center for a two-noded beam
       call getBeamLength (triads(1)%p, triads(2)%p, hydyn%wn, z0, z0, &
            &              supTr, wb, Lb, hydyn%C0b(:,1))
       hydyn%Vb(1) = Lb * hydyn%As(1)
       ierr = 0
    end if

100 if (ierr < 0) call reportError (debugFileOnly_p,'updateHydroDynBody')

  end subroutine updateHydroDynBody


  !!============================================================================
  !> @brief Computes the buoyancy length and center for a two-noded beam.
  !>
  !> @param[in] triad1 Triad art first end of the beam
  !> @param[in] triad2 Triad art second end of the beam
  !> @param[in] normal Normal vector of the water surface in global coordinates
  !> @param[in] h1 Height of the water surface along the normal vector for end 1
  !> @param[in] h2 Height of the water surface along the normal vector for end 2
  !> @param[in] Tlg Local to global transformation matrix for the beam element
  !> @param[out] weight Buoyancy weight factors
  !> @param[out] Lb Length of the submerged part of the beam (buoyancy length)
  !> @param[out] C0b Centre of buoyancy
  !> @param[out] iEnd Which end is below the water surface, if partly submerged.
  !> iEnd < 0 : This element is completely submerged.
  !> iEnd = 0 : This element is completely in the air.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 29 Sep 2008

  subroutine getBeamLength (triad1,triad2,normal,h1,h2,Tlg,weight,Lb,C0b,iEnd)

    use TriadTypeModule, only : TriadType

    type(TriadType)   , intent(in)  :: triad1, triad2
    real(dp)          , intent(in)  :: normal(3), h1, h2
    real(dp), optional, intent(in)  :: Tlg(3,4)
    real(dp), optional, intent(out) :: weight(2)
    real(dp),           intent(out) :: Lb
    real(dp), optional, intent(out) :: C0b(3)
    integer , optional, intent(out) :: iEnd

    !! Local variables
    real(dp) :: z1, z2, xi, x1(3), x2(3)

    !! --- Logic section ---

    if (present(Tlg)) then
       x1 = glob2loc(Tlg,triad1%ur(:,4)) ! Local coordinates of node 1
       x2 = glob2loc(Tlg,triad2%ur(:,4)) ! Local coordinates of node 2
    else
       x1 = triad1%ur(:,4) ! Global coordinates of node 1
       x2 = triad2%ur(:,4) ! Global coordinates of node 2
    end if
    Lb = sqrt(sum((x2-x1)*(x2-x1)))  ! Current length of the beam
    z1 = dot_product(x1,normal) - h1 ! Height above surface for node 1
    z2 = dot_product(x2,normal) - h2 ! Height above surface for node 2
    xi = z1 / (z1-z2)                ! Surface intersection, relative to node 1
    if (z1 <= 0.0_dp .and. z2 <= 0.0_dp) then ! beam is completely submerged
       if (present(weight)) weight = 0.5_dp
       if (present(iEnd)) iEnd = -2
    else if (z1 >= 0.0_dp .and. z2 >= 0.0_dp) then ! beam is above water
       Lb = 0.0_dp
       if (present(weight)) weight = 0.0_dp
       if (present(iEnd)) iEnd = 0
    else if (z1 < 0.0_dp) then ! only first end is submerged
       Lb = Lb * xi
       if (present(weight)) then
          weight(2) = 0.5_dp*xi
          weight(1) = 1.0_dp - weight(2)
       end if
       if (present(iEnd)) iEnd = 1
    else !   z2 < 0.0_dp       , only second end is submerged
       Lb = Lb * (1.0_dp-xi)
       if (present(weight)) then
          weight(1) = 0.5_dp - 0.5_dp*xi
          weight(2) = 1.0_dp - weight(1)
       end if
       if (present(iEnd)) iEnd = 2
    end if

    if (present(C0b) .and. present(weight)) then
       C0b = weight(1)*x1 + weight(2)*x2 ! Buoyancy center
    end if

  end subroutine getBeamLength


  !!============================================================================
  !> @brief Calculates Morison force contributions for a two-noded beam element.
  !>
  !> @param[in] beamId ID string of the beam element used in feedback messages
  !> @param[in] triads Triads connected to the beam elements
  !> @param[in] supTr Position matrix for the beam element
  !> @param[in] urd Nodal velocities in local coordinates
  !> @param[in] urdd Nodal accelerations in local coordinates
  !> @param Q External nodal forces with added mass, drag and buoyancy terms
  !> @param[out] eMa Element added mass matrix
  !> @param[out] eCd Element damping matrix due to drag
  !> @param hydyn Data for hydrodynamic force calculation
  !> @param env Environmental data
  !> @param[in] time Current simulation time
  !> @param[in] istep Time increment counter
  !> @param[in] iter Iteration counter
  !> @param[out] ierr Error flag
  !>
  !> @brief The inertia- and damping forces due to added mass and drag
  !> are calculated for a two-noded element assuming a circular cross section.
  !> The corresponding left-hand side matrix contributions are also calculated.
  !> Buoyancy forces are also calculated, but no load-correction stiffness.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 27 Mar 2009

  subroutine getMorisonForces (beamId,triads,supTr,urd,urdd,Q,eMa,eCd, &
       &                       hydyn,env,time,istep,iter,ierr)

    use kindModule            , only : epsDiv0_p
    use TriadTypeModule       , only : TriadPtrType
    use SupElTypeModule       , only : HydroDynType
    use EnvironmentTypeModule , only : EnvironmentType
    use FunctionTypeModule    , only : FunctionValue
    use IdTypeModule          , only : getId, StrId
    use rotationModule        , only : EccExpand
    use profilerModule        , only : startTimer, stopTimer, hyd_p
    use fileUtilitiesModule   , only : getDBGfile
#ifdef FT_DEBUG
    use dbgUnitsModule        , only : dbgHD
#if FT_DEBUG > 2
    use manipMatrixModule     , only : writeObject
#endif
#endif
    use reportErrorModule     , only : getErrorFile, internalError, reportError
    use reportErrorModule     , only : warning_p, debugFileOnly_p
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getint
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_isTrue

    character(len=*)     , intent(in)    :: beamId
    type(TriadPtrType)   , intent(in)    :: triads(:)
    real(dp)             , intent(in)    :: supTr(:,:), urd(:), urdd(:)
    real(dp)             , intent(inout) :: Q(:)
    real(dp)             , intent(out)   :: eMa(:,:), eCd(:,:)
    type(HydroDynType)   , intent(inout) :: hydyn
    type(EnvironmentType), intent(inout) :: env
    real(dp)             , intent(in)    :: time
    integer              , intent(in)    :: istep, iter
    integer              , intent(out)   :: ierr

    !! Local variables
    integer  :: i, j, k, n1, n2, iEnd, dragForm, dumpWave, dbgW
    real(dp) :: g, xi, Ca56, Ca(3), Cm(3), Cd(4), x1(3), x2(3), e(3,2), wb(2)
    real(dp) :: Lb, Vb, uld(3), uld4, uldd(6,2), vr(4,2), vel(3,2), acc(3,2)
    real(dp) :: Fd(3,2), Fm(3,2), Fs(3,2), CdMat(3,3,2), Cs
    real(dp) :: eta(2), surfaceMotion(9), entryTime, decayTime, rhowEff, v_rel

    !! --- Logic section ---

    eMa = 0.0_dp
    eCd = 0.0_dp
    if (.not. allocated(calcWMotion)) then
       ierr = 1 ! Fluid particle motion cache not initialized yet, ignore
       return
    end if

    call startTimer (hyd_p)

    call ffa_cmdlinearg_getint ('dumpWaveNode',dumpWave)
    call ffa_cmdlinearg_getint ('dragForm',dragForm)

    !! dragForm=0: Old drag formulation, equivalent to R7.0.3  and earlier.
    !! dragForm=1: New formulation based on force resultant in the normal plane
    !!             of the element, and linear variation of the force intensity
    !!             along the beam element (see Fedem theory guide, Appendix B).
    !! dragForm=2: As dragForm=1, but instead of linear interpolation along the
    !!             element, a piece-wise constant force intensity is assumed.
    !! dragForm=3: Drag and added mass forces are as for dragForm=2, whereas the
    !!             the damping matrix contributions are as for dragForm=0.
    !! dragForm=4: As dragForm=1, but with zero damping matrix contributions.
    !! dragForm=5: As dragForm=2, but with zero damping matrix contributions.

    hydyn%B = 0.0_dp
    hydyn%M = 0.0_dp
    hydyn%D = 0.0_dp
    hydyn%S = 0.0_dp

    !! Find gravity constant (g) and water surface normal (wn)
    ierr = -1
    g = waterSurfaceNormal(hydyn%wn,env%gravity)
    if (g <= 0.0_dp) goto 900

#ifdef FT_DEBUG
    if (dbgHD == 0) dbgHD = getDBGfile(5,'hydrodynamics.dbg')
    write(dbgHD,"(/'   --- in getMorisonForces, ',A)",advance='NO') trim(beamId)
    do i = 1, 2
       j = triads(i)%p%inFluid
       if (j > 0) write(dbgHD,"(I3)",advance='NO') calcWMotion(j)
    end do
    write(dbgHD,"(/7X,'time, iter      =',1PE13.5,I4)") time, iter
    write(dbgHD,"(7X,'rhow, g, normal =',1P5E13.5)") env%rhow, g, hydyn%wn
#endif


    !!===========================!!
    !! Sea kinematics evaluation !!
    !!===========================!!

    !! Evaluate the sea state (wave and current) at both nodes
    eta  = 0.0_dp
    vel  = 0.0_dp
    acc  = 0.0_dp
    vr   = 0.0_dp
    uldd = 0.0_dp
    iEnd = 0
    do i = 1, 2
       j = triads(i)%p%inFluid
       if (j < 1) cycle ! dry triad

       if (calcWMotion(j) == 0) then ! Evaluate only once per node per step
          call evaluateSea (env,g,time,istep,j,triads(i)%p%ur(:,4), &
               &            waterMotion(1:9,j),calcWMotion(j))
          ierr = calcWMotion(j)
          if (ierr < 0) goto 900

          if (dumpWave == triads(i)%p%id%userid) then
             dbgW = getDBGfile(8,'Triad'//trim(adjustl(StrId(dumpWave)))// &
                  &            '_wave.asc')
             write(dbgW,600) time,iter,waterMotion(1:9,j)
600          format(1PE12.5,I3,3(1X,3E13.5))
          end if
       end if

       !! Current sea surface height
       if (associated(env%waveFunc)) then
          eta(i) = dot_product(waterMotion(1:3,j),hydyn%wn)
       else
          eta(i) = env%seaLevel
       end if
       waterMotion(10,j) = eta(i)

       if (abs(calcWMotion(j)) == 1) then
          !! Local fluid particle velocity and acceleration
          vel(:,i) = matmul(waterMotion(4:6,j),supTr(:,1:3))
          acc(:,i) = matmul(waterMotion(7:9,j),supTr(:,1:3))
          n1 = triads(i)%p%nDOFs
          if (n1 > 0) then
             !! Relative fluid particle velocity w.r.t. structural motion
             k = triads(i)%firstDOF
             vr(1:3,i) = vel(:,i) - urd(k:k+2)
             if (n1 > 3) vr(4,i) = -urd(k+3)
             !! Structural acceleration
             uldd(1:n1,i) = urdd(k:k+n1-1)
          end if
          iEnd = iEnd + 10
       end if

       if (calcWMotion(j) == 1 .and. waterMotion(11,j) <= 0.0_dp) then
          waterMotion(11,j) = time ! Log the water entry time of this node
#ifdef FT_DEBUG
          write(dbgHD,"(7X,'Triad',A,' is entering the water')") &
               &      trim(getId(triads(i)%p%id))
#endif
       else if (calcWMotion(j) > 1 .and. abs(waterMotion(11,j)) > eps_p) then
          waterMotion(11,j) = 0.0_dp ! This node is no longer in the water
#ifdef FT_DEBUG
          write(dbgHD,"(7X,'Triad',A,' is leaving the water')") &
               &      trim(getId(triads(i)%p%id))
#endif
       end if
       iEnd = iEnd + 1

    end do

#ifdef FT_DEBUG
    write(dbgHD,"(7X,'wave profile    =',1P2E13.5)") eta
    if (iEnd > 10 .or. any(abs(vel) > eps_p) .or. any(abs(acc) > eps_p)) then
       if (any(abs(vel) > eps_p)) then
          write(dbgHD,"(7X,'fluid velocity  =',1P6E13.5)") vel
       end if
       if (any(abs(acc) > eps_p)) then
          write(dbgHD,"(4X,'fluid acceleration =',1P6E13.5)") acc
       end if
       write(dbgHD,"(5X,'relative velocity =',1P10E13.5)") vr
    end if
    write(dbgHD,"(18X,'iEnd =',I3)") iEnd
#endif

    Vb = hydyn%Vb(1)
    if (iEnd > 10) then ! at least one node is in the water
       !! Calculate the intersection with the current water surface, if any,
       !! and the length of the submerged part of the beam element
       call getBeamLength (triads(1)%p, triads(2)%p, hydyn%wn, &
            &              eta(1), eta(2), weight=wb, Lb=Lb, iEnd=iEnd)
       if (iEnd > 0) then
          env%seaLevel = eta(1)*wb(1) + eta(2)*wb(2) ! update current sea level
       end if
       hydyn%Vb(1) = Lb * hydyn%As(1)
    else if (iEnd == 1) then
       iEnd = 0
       hydyn%Vb(1) = 0.0_dp
       call reportError (warning_p,trim(beamId)// &
            ' has one end grounded and the other end above water.', &
            'No hydrodynamic forces will be calculated for this element.', &
            'Consider refining the model if this is not correct.')
    else
       iEnd = 0 ! both nodes are either grounded or above the water, OK
       hydyn%Vb(1) = 0.0_dp
    end if

    ierr = 0
    if (iEnd == 0) then
       if (Vb > 0.0_dp) then
          write(getErrorFile(),610) trim(beamId),'leaving'
       end if
       goto 910
    else if (abs(Vb) <= eps_p .and. hydyn%As(1) > 0.0_dp) then
       write(getErrorFile(),610) trim(beamId),'entering'
    end if
610 format(5X,A,' is ',A,' the water')

    Lb = 0.5_dp * Lb ! Half length of the submerged part
#ifdef FT_DEBUG
    write(dbgHD,"(7X,'Vb, Lb, wb iEnd =',1P4E13.5,I4)") hydyn%Vb(1),Lb,wb,iEnd
#endif
    if (ffa_cmdlinearg_isTrue('ignoreHD')) goto 910

    e = 0.0_dp
    if (iEnd > 0) then
       iEnd = 3 - iEnd ! the local node that is above water
       xi = 2.0_dp*wb(2) - 2.0_dp + real(iEnd,dp)
       x1 = triads(1)%p%ur(:,4)
       x2 = triads(2)%p%ur(:,4)

       !! Eccentricity transformation to the beam ends
       if (iEnd == 2) then
          e(:,2) = (1.0_dp-xi)*glob2loc(supTr,x2-x1)
       else if (iEnd == 1) then
          e(:,1) = -xi*glob2loc(supTr,x2-x1)
       end if

       !! Since one beam end is above the water, we need to re-evaluate the
       !! water kinematics at the intersection point to get correct forces
       ierr = 2
       call evaluateSea (env,g,time,istep,0,x1+xi*(x2-x1),surfaceMotion,ierr)
       if (ierr < 0) goto 900

       if (dumpWave == triads(iEnd)%p%id%userid) then
          dbgW = getDBGfile(8,'Triad'//trim(adjustl(StrId(dumpWave)))// &
               &            '_wave.asc')
          write(dbgW,600) time,iter,surfaceMotion
       end if

       !! Local fluid particle velocity and acceleration
       vel(:,iEnd) = matmul(surfaceMotion(4:6),supTr(:,1:3))
       acc(:,iEnd) = matmul(surfaceMotion(7:9),supTr(:,1:3))

       !! Relative fluid particle velocity w.r.t. structural motion
       i = triads(1)%firstDOF
       j = triads(2)%firstDOF
       n1 = triads(1)%p%nDOFs
       n2 = triads(2)%p%nDOFs
       if (n1 >= 3 .and. n2 >= 3) then
          uld = urd(i:i+2)*(1.0_dp-xi) + urd(j:j+2)*xi
       else if (n1 >= 3) then
          uld = urd(i:i+2)*(1.0_dp-xi)
       else if (n2 >= 3) then
          uld = urd(j:j+2)*xi
       else
          uld = 0.0_dp
       end if
       if (n1 > 3 .and. n2 > 3) then
          uld4 = urd(i+3)*(1.0_dp-xi) + urd(j+3)*xi
       else if (n1 > 3) then
          uld4 = urd(i+3)*(1.0_dp-xi)
       else if (n2 > 3) then
          uld4 = urd(j+3)*xi
       else
          uld4 = 0.0_dp
       end if
       vr(1:3,iEnd) = vel(:,iEnd) - uld
       vr(4,iEnd) = -uld4

       !! Structural acceleration
       uldd(:,iEnd) = 0.0_dp
       if (n1 > 0) then
          uldd(1:n1,iEnd) = urdd(i:i+n1-1)*(1.0_dp-xi)
       end if
       if (n2 > 0) then
          uldd(1:n2,iEnd) = uldd(1:n2,iEnd) + urdd(j:j+n2-1)*xi
       end if

       j = triads(iEnd)%p%inFluid
       if (calcWMotion(j) > 0 .and. abs(waterMotion(11,j)) <= eps_p) then
          waterMotion(11,j) = -time
#ifdef FT_DEBUG
          write(dbgHD,"(7X,'Triad',A,' is not yet in the water')") &
               &      trim(getId(triads(iEnd)%p%id))
#endif
       end if

#ifdef FT_DEBUG
       write(dbgHD,"(7X,'reevaluated surface kinematics for end',I3)") iEnd
       write(dbgHD,"(7X,'fluid velocity  =',1P3E13.5)") vel(:,iEnd)
       write(dbgHD,"(4X,'fluid acceleration =',1P3E13.5)") acc(:,iEnd)
       write(dbgHD,"(5X,'relative velocity =',1P5E13.5)") vr(:,iEnd)
#endif
    end if

    !!==================================!!
    !! End of sea kinematics evaluation !!
    !!==================================!!


    !! Effective fluid density (possibly scaled to ramp up hydrodynamic forces)
    rhowEff = env%rhow*env%hdfScale

    !! Set up diagonal mass and damping matrices
    Ca(1)   = rhowEff * hydyn%Morison(6) * Lb
    Ca(2:3) = rhowEff * hydyn%Morison(1) * Lb
    Ca56    = Ca(2) * Lb*Lb/52.5_dp ! rotational added mass
    Cm(1)   = rhowEff * hydyn%Morison(7) * Lb
    Cm(2:3) = rhowEff * hydyn%Morison(2) * Lb
    Cd(1)   = rhowEff * hydyn%Morison(4) * Lb
    Cd(2:3) = rhowEff * hydyn%Morison(3) * Lb
    Cd(4)   = rhowEff * hydyn%Morison(5) * Lb
    Cs      = rhowEff * hydyn%Morison(8) * Lb
#if FT_DEBUG > 1
    write(dbgHD,"(7X,'Ca =',1P3E13.5)") Ca(1:2),Ca56
    write(dbgHD,"(7X,'Cm =',1P2E13.5)") Cm(1:2)
    write(dbgHD,"(7X,'Cd =',1P3E13.5)") Cd(1:2),Cd(4)
    write(dbgHD,"(7X,'Cs =',1P2E13.5)") Cs
#endif

    !! Buoyancy force
    hydyn%wn = matmul(hydyn%wn,supTr(:,1:3))
    hydyn%B(1:3) = env%rhow*g * hydyn%Vb(1)*hydyn%wn


    !!==================================================!!
    !! Evaluate nodal hydrodynamics force contributions !!
    !! and the associated tangent matrix contributions  !!
    !!==================================================!!

    Fs = 0.0_dp
    CdMat = 0.0_dp
    do i = 1, 2
       j = triads(i)%firstDOF

       !! Find the relative cross-flow velocity magnitude
       v_rel = sqrt(vr(2,i)*vr(2,i)+vr(3,i)*vr(3,i))

       if (dragForm > 0) then
          !! Added mass forces
          Fm(:,i) = Cm*acc(:,i) - Ca*uldd(1:3,i)
          !! Drag force in axial direction
          Fd(1,i) = Cd(1)*abs(vr(1,i))*vr(1,i)
          !! Normal drag forces in element coordinate system
          Fd(2:3,i) = Cd(2)*v_rel*vr(2:3,i)
       else
          !! Added mass forces
          Fm(:,i) = Cm*acc(:,i)
          !! Old drag formulation (R7.0 compatible, incorrect for skew flows)
          Fd(:,i) = Cd(1:3)*abs(vr(1:3,i))*vel(:,i)
       end if

       if (Cs > 0.0_dp .and. triads(i)%p%inFluid > 0) then
          !! Check if a slamming contribution should be added
          entryTime = abs(waterMotion(11,triads(i)%p%inFluid))
          decayTime = entryTime + hydyn%Morison(9)
          if (time >= max(entryTime,eps_p) .and. time <= decayTime) then
             if (dragForm > 0) then
                Fs(2:3,i) = Cs*v_rel*vr(2:3,i)
             else
                Fs(2:3,i) = Cs*abs(vr(2:3,i))*vr(2:3,i)
             end if
          end if
       end if

       !! Store the drag and added mass forces (for plotting only)
       hydyn%D(1:3) = hydyn%D(1:3) + Fd(:,i)
       hydyn%M(1:3) = hydyn%M(1:3) + Fm(:,i)
       hydyn%S(1:3) = hydyn%S(1:3) + Fs(:,i)

       if (dragForm == 1 .or. dragForm == 2) then
          !! Drag damping matrix. It is defined via the derivative of the
          !! drag force (Fd) with respect to the structure velocity (uld).
          !! See eqs. (B.20)-(B.23) in the Fedem R7.1 theory guide.
          CdMat(1,1,i) = Cd(1)*2.0_dp*abs(vr(1,i))
          if (v_rel > epsDiv0_p) then
             CdMat(2,2,i) = Cd(2)*(vr(2,i)*vr(2,i)/v_rel + v_rel)
             CdMat(3,3,i) = Cd(2)*(vr(3,i)*vr(3,i)/v_rel + v_rel)
             CdMat(2,3,i) = Cd(2)*(vr(2,i)*vr(3,i)/v_rel)
             CdMat(3,2,i) = CdMat(2,3,i)
          end if
       end if

       n1 = triads(i)%p%nDOFs
       if (n1 < 3) cycle ! Grounded triad, no contributions

#ifdef FT_DEBUG
       write(dbgHD,"(7X,'G',I1,' =',1P3E13.5)") i,Q(j:j+2) ! gravity force
       write(dbgHD,"(7X,'B',I1,' =',1P3E13.5)") i,hydyn%B(1:3)*wb(i)
#endif

       !! Add external force vector due to buoyancy
       Q(j:j+2) = Q(j:j+2) + hydyn%B(1:3)*wb(i)

       if (dragForm /= 1 .and. dragForm /= 4) then
          !! Add external force vector due to fluid motion,
          !! assuming piece-wise constant variation along the beam
          Q(j:j+2) = Q(j:j+2) + Fm(:,i) + Fd(:,i) + Fs(:,i)
          if (n1 >= 6) then
             Q(j+4:j+5) = Q(j+4:j+5) - Ca56*uldd(5:6,i) ! Added mass moments
          end if
       end if

       if (dragForm == 2) then
          if (n1 >= 6) then
             call EccExpand (e(:,i),CdMat(:,:,i),eCd(j:j+5,j:j+5))
          else
             eCd(j:j+2,j:j+2) = CdMat(:,:,i)
          end if
       else if (dragForm < 1 .or. dragForm == 3) then
          !! Drag damping matrix. It is defined via the derivative of the
          !! drag force (Fd) with respect to the structure velocity (uld).
          !! But since Fd is a quadratic function of the velocity (vr)
          !! we get the factor 2.0 in the damping matrix expression here.
          x1 = Cd(1:3)*2.0_dp*abs(vr(1:3,i)) ! diagonal damping matrix
          if (n1 >= 6) then
             call EccExpand (e(:,i),x1,eCd(j:j+5,j:j+5))
          else
             do k = 0, 2
                eCd(j+k,j+1) = x1(1+k)
             end do
          end if
       end if
       if (dragForm < 1 .or. dragForm == 2 .or. dragForm == 3) then
          !! Drag damping matrix term due to axial spin velocity
          if (n1 > 3) then
             eCd(j+3,j+3) = eCd(j+3,j+3) + Cd(4)*2.0_dp*abs(vr(4,i))
          end if
       end if

       !! Added mass matrix
       if (n1 >= 6) then
          call EccExpand (e(:,i),Ca,eMa(j:j+5,j:j+5))
          do k = 4, 5
             eMa(j+k,j+k) = eMa(j+k,j+k) + Ca56
          end do
       else
          do k = 0,2
             eMa(j+k,j+k) = Ca(1+k)
          end do
       end if

    end do

    if (dragForm == 1 .or. dragForm == 4) then
       !! Assuming linear variation along the beam
       do i = 1, 2
          n1 = triads(i)%p%nDOFs
          if (n1 < 3) cycle ! Grounded triad, no contributions

          j = triads(i)%firstDOF
          k = 3 - i

          !! Add external force vector due to fluid motion
          Q(j:j+2) = Q(j:j+2) + &
               &    ((Fm(:,i) + Fd(:,i) + Fs(:,i))*2.0_dp + &
               &      Fm(:,k) + Fd(:,k) + Fs(:,k))/3.0_dp
          if (n1 >= 6) then
             Q(j+4:j+5) = Q(j+4:j+5) & ! Added mass moments
                  &     - Ca56*(uldd(5:6,i)*2.0_dp + uldd(5:6,k))/3.0_dp
          end if
          if (dragForm == 1) then
             !! Drag damping matrix
             eCd(j:j+2,j:j+2) = (2.0_dp*CdMat(:,:,i) + CdMat(:,:,k))/3.0_dp
             if (n1 >= 6) then
                call EccExpand (e(:,i),eCd(j:j+2,j:j+2),eCd(j:j+5,j:j+5))
                eCd(j+3,j+3) = eCd(j+3,j+3) &
                     + Cd(4)*2.0_dp*(2.0_dp*abs(vr(4,i))+abs(vr(4,k)))/3.0_dp
             end if
          end if
       end do
    end if

#ifdef FT_DEBUG
    write(dbgHD,"(7X,'Fm =',1P6E13.5)") Fm
    write(dbgHD,"(7X,'Fd =',1P6E13.5)") Fd
    if (any(abs(Fs) > eps_p)) then
       write(dbgHD,"(7X,'Fs =',1P6E13.5)") Fs
    end if
    do i = 1, 2
       n1 = triads(i)%p%nDOFs
       if (n1 > 0) then
          j = triads(i)%firstDOF
          write(dbgHD,"(7X,'Q',I1,' =',1P3E13.5)") i,Q(j:j+2)
#if FT_DEBUG > 2
          call writeObject (eMa(j:j+n1-1,j:j+n1-1),dbgHD,'MaMat')
          call writeObject (eCd(j:j+n1-1,j:j+n1-1),dbgHD,'CdMat')
#endif
       end if
    end do
#endif

    ierr = 0
    goto 910

900 call reportError (debugFileOnly_p,'getMorisonForces')
910 call stopTimer (hyd_p)
    do i = 1, 2
       j = triads(i)%p%inFluid
       if (j > 0) calcWMotion(j) = -abs(calcWMotion(j))
    end do

  end subroutine getMorisonForces


  !!============================================================================
  !> @brief Calculates buoyancy force resultant for a superelement or beam.
  !>
  !> @param[in] supId ID string of the superelement used in feedback messages
  !> @param[in] triads Triads connected to the superelement
  !> @param[in] supTr Position matrix for the superelement
  !> @param hydyn Data for hydrodynamic force calculation
  !> @param env Environmental data
  !> @param[out] g Gravitation constant
  !> @param[in] time Current simulation time
  !> @param[in] iter Iteration counter
  !> @param Q External nodal forces for beams, including buoyancy on output
  !> @param[out] ierr Error flag
  !>
  !> @details The contributions to the external load vector form the hydrostatic
  !> buoyancy may optionally be calculated, but only for two-noded elements.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 2 Jul 2008

  subroutine getBuoyancyForces (supId,triads,supTr,hydyn,env,g,time,iter,Q,ierr)

    use TriadTypeModule      , only : TriadPtrType
    use SupElTypeModule      , only : HydroDynType
    use EnvironmentTypeModule, only : EnvironmentType
    use profilerModule       , only : startTimer, stopTimer, hyd_p
#ifdef FT_DEBUG
    use fileUtilitiesModule  , only : getDBGfile
    use dbgUnitsModule       , only : dbgHD
#endif
    use reportErrorModule    , only : reportError, debugFileOnly_p
    use reportErrorModule    , only : internalError, getErrorFile

    character(len=*)     , intent(in)    :: supId
    type(TriadPtrType)   , intent(in)    :: triads(:)
    real(dp)             , intent(in)    :: supTr(:,:)
    type(HydroDynType)   , intent(inout) :: hydyn
    type(EnvironmentType), intent(inout) :: env
    real(dp)             , intent(out)   :: g
    real(dp)             , intent(in)    :: time
    integer              , intent(in)    :: iter
    real(dp), optional   , intent(inout) :: Q(:)
    integer              , intent(out)   :: ierr

    !! Local variables
    integer  :: i, j
    real(dp) :: As, Vb, wb(2)

    !! --- Logic section ---

    call startTimer (hyd_p)

    hydyn%B = 0.0_dp
    hydyn%M = 0.0_dp
    hydyn%D = 0.0_dp
    hydyn%S = 0.0_dp

    As = hydyn%As(1)
    Vb = hydyn%Vb(1)

#ifdef FT_DEBUG
    if (dbgHD == 0) dbgHD = getDBGfile(5,'hydrodynamics.dbg')
#endif
    if (associated(env%waveFunc)) then
       !! Get sea surface elevation at the origin of this part
       env%seaLevel = env%sea0 + getWaveElevation(env,supTr(:,4),time,ierr)
       if (ierr < 0) goto 900
    end if

    !! Update the hydrodynamic body quantities for current configuration
    call updateHydroDynBody (hydyn,triads,supTr, &
         &                   env%seaLevel,env%gravity,g,wb,ierr)
    if (ierr < 0) goto 900

    if (abs(hydyn%Vb(1)) <= eps_p) then
       if (Vb <= 0.0_dp) goto 800 ! The entire body is above the water surface
       write(getErrorFile(),610) trim(supId),'leaving the water'
    else if (hydyn%Vb(1) > 0.0_dp .and. abs(Vb) <= eps_p) then
       write(getErrorFile(),610) trim(supId),'entering the water'
    else if (abs(hydyn%As(1)) <= eps_p .and. As > 0.0_dp) then
       write(getErrorFile(),610) trim(supId),'totally submerged'
    else if (hydyn%As(1) > 0.0_dp .and. abs(As) <= eps_p) then
       write(getErrorFile(),610) trim(supId),'partially submerged'
    end if
610 format(5X,A,' is ',A)
#ifdef FT_DEBUG
    write(dbgHD,"(/'   --- in getBuoyancyForces, body',I4)") hydyn%bodyIndex
    write(dbgHD,"(7X,'time, iter      =',1PE13.5,I4)") time, iter
    write(dbgHD,"(7X,'rhow, g, normal =',1P5E13.5)") env%rhow, g, hydyn%wn
    write(dbgHD,"(7X,'Vb, C0b         =',1P8E13.5)") hydyn%Vb,hydyn%C0b
    write(dbgHD,"(37X,1P3E13.5)") loc2glob(supTr,hydyn%C0b(:,1))
    if (hydyn%As(1) > 0.0_dp) then
       write(dbgHD,"(7X,'As, C0As        =',1P8E13.5)") hydyn%As, hydyn%C0As
       write(dbgHD,"(50X,1P3E13.5)") loc2glob(supTr,hydyn%C0As(:,1))
    end if
#else
    i = iter ! Dummy statement
#endif

    if (present(Q) .and. size(triads) == 2) then

       !! Two-noded beam, calculate the buoyancy force at each node
       hydyn%B(1:3) = env%rhow*g*hydyn%Vb(1)*hydyn%wn
       do i = 1, 2
          if (triads(i)%p%nDOFs >= 3) then
             j = triads(i)%firstDOF
#ifdef FT_DEBUG
             write(dbgHD,"(7X,'G',I1,' =',1P3E13.5)") i,Q(j:j+2)
             write(dbgHD,"(7X,'B',I1,' =',1P3E13.5)") i,hydyn%B(1:3)*wb(i)
             if (j == 1) then
                write(dbgHD,"(11X,1P3E13.5)") matmul(supTr(:,1:3),hydyn%B(1:3))
             end if

#endif
             Q(j:j+2) = Q(j:j+2) + hydyn%B(1:3)*wb(i)
#ifdef FT_DEBUG
             write(dbgHD,"(7X,'Q',I1,' =',1P3E13.5)") i,Q(j:j+2)
#endif
          end if
       end do

    else if (present(Q)) then

       !! For FE parts we would need to distribute the buoyancy force over the
       !! connected triads, depending on how much volume to associate with each
       ierr = internalError('getBuoyancyForces: Not implemented for FE parts')

    end if

800 call stopTimer (hyd_p)
    return

900 call reportError (debugFileOnly_p,'getBuoyancyForces')
    goto 800

  end subroutine getBuoyancyForces


  !!============================================================================
  !> @brief Calculates drag (and slam) forces for a (partly) submerged body.
  !>
  !> @param Fd Damping forces at the centre of gravity due to drag and slam
  !> @param[in] sup Superelement to calculate damping forces for
  !> @param hydyn Data for hydrodynamic force calculation
  !> @param[in] env Environmental data
  !> @param[in] time Current simulation time
  !> @param[in] dt Time increment size
  !> @param ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 8 Jul 2008

  subroutine getDragForces (Fd,sup,hydyn,env,time,dt,ierr)

    use SupElTypeModule        , only : SupElType, HydroDynType
    use EnvironmentTypeModule  , only : EnvironmentType
    use rotationModule         , only : EccExpand
    use profilerModule         , only : startTimer, stopTimer, hyd_p
#ifdef FT_DEBUG
    use dbgUnitsModule         , only : dbgHD
#endif
    use reportErrorModule      , only : reportError, debugFileOnly_p
    use reportErrorModule      , only : internalError
    use FFaBodyHandlerInterface, only : ffa_total_volume
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_getint

    real(dp)             , intent(inout) :: Fd(:)
    type(SupElType)      , intent(in)    :: sup
    type(HydroDynType)   , intent(inout) :: hydyn
    type(EnvironmentType), intent(in)    :: env
    real(dp)             , intent(in)    :: time, dt
    integer              , intent(inout) :: ierr

    !! Local variables
    logical  :: haveDrag, haveSlam
    integer  :: i, j, err, slamForm
    real(dp) :: Vol, Vb, dAdt, u_n, uld(6), F(6), eVec(3)

    !! --- Logic section ---

    call startTimer (hyd_p)

    !! Compute total volume of the body
    call ffa_total_volume (hydyn%bodyIndex, Vol, uld, err)
    if (err < 0) goto 900

    !! Normalized buoyancy volume (0: completely above, 1: fully submerged)
    Vb = hydyn%Vb(1)/Vol
    if (Vb <= 0.0_dp) goto 910

#ifdef FT_DEBUG
    write(dbgHD,"(/'   --- in getDragForces, body',I4)") hydyn%bodyIndex
    write(dbgHD,"(7X,'Vol, Vb(ratio)  =',1P5E13.5)") Vol, Vb
#endif

    if (sup%rigidFlag > 0) then

       haveDrag = .false.
       haveSlam = .false.

       !! Find the overall velocity in the local link coordinate system
       uld = sup%triads(1)%p%urd
       do i = 2, size(sup%triads)
          uld = uld + sup%triads(i)%p%urd
       end do
       uld(1:3) = matmul(uld(1:3)/size(sup%triads),sup%supTr(:,1:3))
       uld(4:6) = matmul(uld(4:6)/size(sup%triads),sup%supTr(:,1:3))

       if (  abs(hydyn%slamPar(1)) > eps_p .and. &
            (hydyn%slamPar(3) <= eps_p .or. time < hydyn%slamPar(3)) ) then
          u_n = dot_product(hydyn%wn,uld(1:3)) ! velocity in normal direction
          if (u_n < 0.0_dp) then
             call ffa_cmdlinearg_getint ('slamForm',slamForm)
             if (slamForm < 2) then
                dAdt = hydyn%As(1) - hydyn%As(2) ! waterline area increment
                haveSlam = dAdt > hydyn%As(1)*hydyn%slamPar(2)
             else
                dAdt = hydyn%Vb(1) - hydyn%Vb(2) ! buoyancy volume increment
                haveSlam = dAdt > hydyn%Vb(1)*hydyn%slamPar(2)
             end if
          end if
       end if

       if (associated(hydyn%dragPar)) then

          !! Local drag force vector at the buoyancy center
          call computeLocalDragForce (F,hydyn%dragPar,uld)
#ifdef FT_DEBUG
          write(dbgHD,"(7X,'uld             =',1P7E13.5)") uld
          write(dbgHD,"(7X,'Fdrag           =',1P6E13.5)") F
#endif

          !! Drag force vector in the center of gravity triad
          eVec = hydyn%C0b(:,1) - sup%Trundeformed(:,4,1)
          call EccExpand (eVec,F,hydyn%D)
          hyDyn%D(4:6) = hyDyn%D(4:6) + F(4:6)
#ifdef FT_DEBUG
          write(dbgHD,"(7X,'D  =',1P6E13.5)") hydyn%D
          write(dbgHD,"(11X,1P3E13.5)") matmul(sup%supTr(:,1:3),hydyn%D(1:3))
#endif
       end if

       !! Check for local triad contributions too (unevenly distributed drag)
       do i = 2, size(sup%triads)
          if ( associated(sup%triads(i)%p%dragPar) .and. &
               sup%triads(i)%p%ur(3,4) < 0.0_dp ) then

             !! This triad is below the water surface, find the
             !! local drag force attacking at this triad position
             uld(1:3) = matmul(sup%triads(i)%p%urd(1:3),sup%supTr(:,1:3))
             call computeCoupledDragForce (F,sup%triads(i)%p%dragPar,uld)
#ifdef FT_DEBUG
             write(dbgHD,"(7X,'uld_',I1,'           =',1P3E13.5)") i,uld(1:3)
             write(dbgHD,"(7X,'Fdrag_',I1,'         =',1P3E13.5)") i,F(1:3)
#endif

             !! Transform to the center of gravity triad
             eVec = sup%Trundeformed(:,4,i) - sup%Trundeformed(:,4,1)
             call EccExpand (eVec,F(1:3),F)
             hyDyn%D = hyDyn%D + F
#ifdef FT_DEBUG
             write(dbgHD,"(7X,'D_',I1,'             =',1P6E13.5)") i,F
#endif
          end if
       end do

       if (haveDrag) then
          !! Add to the superelement damping force vector
          if (associated(sup%Bgp)) then
             Fd = Fd - matmul(hydyn%D,sup%Bgp)
          else
             Fd(1:6) = Fd(1:6) - hydyn%D
          end if
       end if

       if (haveSlam) then

          !! Find the slam contribution, it is a force attacking in the center
          !! of the added waterline area in the direction of the surface normal
          if (slamForm == 2) then
             eVec = (hydyn%Vb(1)*loc2glob(sup%supTr,hydyn%C0b(:,1)) - &
                  &  hydyn%Vb(2)*loc2glob(sup%supTrPrev,hydyn%C0b(:,2))) / dAdt
             hydyn%C0s = glob2loc(sup%supTr,eVec)
          else if (slamForm == 1) then
             eVec = (hydyn%As(1)*loc2glob(sup%supTr,hydyn%C0As(:,1)) - &
                  &  hydyn%As(2)*loc2glob(sup%supTrPrev,hydyn%C0As(:,2))) / dAdt
             hydyn%C0s = glob2loc(sup%supTr,eVec)
          else if (slamForm == -1 .and. hydyn%bodyIndex >= 0) then
             call ffa_inc_area (hydyn%bodyIndex, hydyn%wn, sup%supTr, &
                  &             dAdt, hydyn%C0s, err)
             if (err < 0) goto 900
          else
             hydyn%C0s = hydyn%C0As(:,1)
          end if

          dAdt = dAdt / dt ! waterline area rate
          F(1:3) = env%rhow*abs(hydyn%slamPar(1))*dAdt*dAdt * hydyn%wn
#ifdef FT_DEBUG
          write(dbgHD,"(7X,'dAdt            =',1P,E13.5)") dAdt
          write(dbgHD,"(7X,'C0(dA)          =',1P3E13.5)") hydyn%C0s
          write(dbgHD,"(24X,1P3E13.5)") loc2glob(sup%supTr,hydyn%C0s)
          write(dbgHD,"(7X,'Fslam           =',1P3E13.5)") F(1:3)
#endif

          !! Transform to the center of gravity triad
          eVec = hydyn%C0s - sup%Trundeformed(:,4,1)
          call EccExpand (eVec,F(1:3),hydyn%S)
#ifdef FT_DEBUG
          write(dbgHD,"(7X,'S  =',1P6E13.5)") hydyn%S
          write(dbgHD,"(11X,1P3E13.5)") matmul(sup%supTr(:,1:3),hydyn%S(1:3))
#endif

          !! Add to the superelement damping force vector
          if (associated(sup%Bgp)) then
             Fd = Fd - matmul(hydyn%S,sup%Bgp)
          else
             Fd(1:6) = Fd(1:6) - hydyn%S
          end if

       end if

    else

       !! For FE parts we would need to distribute the drag force over the
       !! connected triads, depending on how much volume to associate with each
       ierr = internalError('getDragForces: Not implemented for FE parts')

    end if
    goto 910

900 ierr = ierr + err
    call reportError (debugFileOnly_p,'getDragForces')
910 call stopTimer (hyd_p)

  contains

    subroutine computeLocalDragForce (Fdrag,dragParam,uld)
      real(dp), intent(out)   :: Fdrag(:)
      real(dp), intent(in)    :: dragParam(:,:)
      real(dp), intent(inout) :: uld(:)
      Fdrag = 0.0_dp
      do j = 1, size(dragParam,2)
         if (dragParam(1,j) > 0.0_dp) then
            if (dragParam(3,j) > 0.0_dp .and. abs(uld(j)) > dragParam(3,j)) then
               !! Cut-off for very high velocities
               uld(j) = sign(dragParam(3,j),uld(j))
            end if
            Fdrag(j) = env%rhow*sign(dragParam(1,j)*uld(j)*uld(j),-uld(j))
            if (dragParam(2,j) > 0.0_dp) then
               !! Account for reduced drag for only a partly submerged body
               Fdrag(j) = Fdrag(j)*max(Vb,1.0_dp-dragParam(2,j))
            end if
            haveDrag = .true.
         end if
      end do
    end subroutine computeLocalDragForce

    subroutine computeCoupledDragForce (Fdrag,dragMat,uld)
      real(dp), intent(out)   :: Fdrag(3)
      real(dp), intent(in)    :: dragMat(3,3)
      real(dp), intent(inout) :: uld(3)
      Fdrag = 0.0_dp
      do j = 1, 3
         if (dragMat(j,j) > 0.0_dp) then
            Fdrag = Fdrag + Vb*env%rhow*dragMat(:,j)*sign(uld(j)*uld(j),-uld(j))
            haveDrag = .true.
         end if
      end do
    end subroutine computeCoupledDragForce

  end subroutine getDragForces


  !!============================================================================
  !> @brief Updates the HydroDynType object after convergence has been achieved.
  !>
  !> @param hydyn Data for hydrodynamic force calculation
  !> @param[in] supTr Position matrix for the superelement to update for
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 11 Aug 2008

  subroutine UpdateHydroDynamicsAtConvergence (hydyn,supTr,ierr)

    use SupElTypeModule        , only : HydroDynType
    use FFaBodyHandlerInterface, only : ffa_save_intersection
    use reportErrorModule      , only : reportError, debugFileOnly_p

    type(HydroDynType), intent(inout) :: hydyn
    real(dp)          , intent(in)    :: supTr(3,4)
    integer           , intent(out)   :: ierr

    !! --- Logic section ---

    ierr = 0
    hydyn%Vb(2) = hydyn%Vb(1)
    hydyn%C0b(:,2) = hydyn%C0b(:,1)
    if (hyDyn%bodyIndex < 0) return

    !! Save data for the current intersection plane for calculation of
    !! centroid of the area increment (used for slam force calculations)
    hydyn%As(2) = hydyn%As(1)
    hydyn%C0As(:,2) = hydyn%C0As(:,1)
    call ffa_save_intersection (hyDyn%bodyIndex,supTr,ierr)
    if (ierr < 0) then
       call reportError (debugFileOnly_p,'UpdateHydroDynamicsAtConvergence')
    end if

  end subroutine UpdateHydroDynamicsAtConvergence


  !!============================================================================
  !> @brief Closes the hydrodynamics module and report some timings.
  !>
  !> @param[in] lpu File unit number for res-file output
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 9 Jan 2012

  subroutine closeHydroDyn (lpu)

    use DiffractionModule, only : NemohClose

    integer, intent(in) :: lpu

    !! --- Logic section ---

    call HWAFLS_close ()
    call NemohClose ()

    if (any(wavCall > 0)) write(lpu,600)
    if (wavCall(1) > 0) write(lpu,610) 'hardware',wavCall(1),wavTime(1,:)
    if (wavCall(2) > 0) write(lpu,610) 'software',wavCall(2),wavTime(2,:)

    if (allocated(waterMotion)) deallocate(waterMotion)
    if (allocated(calcWMotion)) deallocate(calcWMotion)

600 format(//47X,'# calls Wall time  CPU time' &
         &  / 4X,'+',39('-'),'+',3(9('-'),'+'))
610 format(5X,'Wave kinematics evaluations in ',A,I10,2F10.2)

  end subroutine closeHydroDyn


  !!============================================================================
  !> @brief Performs diffraction analysis using Nemoh.
  !>
  !> @param[in] env Environmental data
  !> @param[in] waveFunc The wave function to use in the diffraction analysis
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 06 Mar 2015

  subroutine diffractionCalc (env,waveFunc,ierr)

    use EnvironmentTypeModule  , only : EnvironmentType
    use FunctionTypeModule     , only : FunctionType
    use explicitFunctionsModule, only : WAVE_SINUS_p, SINUSOIDAL_p
    use explicitFunctionsModule, only : COMPL_SINUS_p, DELAYED_COMPL_SINUS_p
    use waveFunctionsModule    , only : twoPi
    use DiffractionModule      , only : NemohAnalysis
    use reportErrorModule      , only : debugFileOnly_p, error_p, reportError

    type(EnvironmentType), intent(in)  :: env
    type(FunctionType)   , intent(in)  :: waveFunc
    integer              , intent(out) :: ierr

    !! Local variables
    integer           :: n, m
    real(dp), target  :: sFunc(6)
    real(dp), pointer :: rFunc(:)

    !! --- Logic section ---

    n = 0
    m = 3
    select case (waveFunc%type)
    case (WAVE_SINUS_p) ! Irregular waves
       n = waveFunc%intParameters(3)
       m = size(waveFunc%realParameters)/n
       rFunc => waveFunc%realParameters
    case (SINUSOIDAL_p) ! Regular waves
       n = 1
       sFunc(1) =  waveFunc%realParameters(3)
       sFunc(2) =  waveFunc%realParameters(1) * twoPi
       sFunc(3) = -waveFunc%realParameters(2) * twoPi
       rFunc => sFunc(1:3)
    case (COMPL_SINUS_p, DELAYED_COMPL_SINUS_p) ! Two-component sine wave
       n = 2
       sFunc(1) =  waveFunc%realParameters(5)
       sFunc(2) =  waveFunc%realParameters(1) * twoPi
       sFunc(3) = -waveFunc%realParameters(3) * twoPi
       sFunc(4) =  waveFunc%realParameters(6)
       sFunc(5) =  waveFunc%realParameters(2) * twoPi
       sFunc(6) = -waveFunc%realParameters(4) * twoPi
       rFunc => sFunc
    end select
    if (n < 1) then
       ierr = -1
       call reportError (error_p,'Invalid wave function for diffraction', &
            &            addString='diffractionCalc')
       return
    end if

    !! Perform the diffraction analysis
    call NemohAnalysis ('NemohDB',env,reshape(rFunc(1:m*n),(/m,n/)),ierr=ierr)
    if (ierr < 0) call reportError (debugFileOnly_p,'diffractionCalc')

  end subroutine diffractionCalc


  !!============================================================================
  !> @brief Extracts the diffraction force at a given time for a superelement.
  !>
  !> @param Q External forces at centre of gravity including diffraction effects
  !> @param[in] sup Superelement to extract diffraction forces
  !> @param hydyn Data for hydrodynamic force calculation
  !> @param[in] env Environmental data
  !> @param[in] time Current simulation time
  !> @param ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 16 Mar 2015

  subroutine getDiffractionForces (Q,sup,hydyn,env,time,ierr)

    use SupElTypeModule      , only : SupElType, HydroDynType
    use EnvironmentTypeModule, only : EnvironmentType
    use DiffractionModule    , only : getExcitationForce

    real(dp)             , intent(inout) :: Q(:)
    type(SupElType)      , intent(in)    :: sup
    type(HydroDynType)   , intent(inout) :: hydyn
    type(EnvironmentType), intent(in)    :: env
    real(dp)             , intent(in)    :: time
    integer              , intent(inout) :: ierr

    !! Local variables
    integer  :: lerr
    real(dp) :: Fex(6), Tsl(3,3)

    !! --- Logic section ---

    call getExcitationForce (hydyn%bodyIndex,time,Fex,lerr)
    if (lerr > 0) return ! No diffraction force for this superelement

    !! Transform from sea coordinate system to local axes of the superelement
    Tsl = matmul(transpose(sup%supTr(:,1:3)),env%Tsea(:,1:3))
    hydyn%D(1:3) = matmul(Tsl,Fex(1:3))
    hydyn%D(4:6) = matmul(Tsl,Fex(4:6))

    Q(1:6) = Q(1:6) + hydyn%D ! Add to external forces at the centre of gravity

    ierr = ierr + lerr

  end subroutine getDiffractionForces

end module HydroDynamicsModule
