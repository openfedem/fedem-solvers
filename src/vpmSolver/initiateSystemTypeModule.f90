!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file initiateSystemTypeModule.f90
!> @brief Initialization of system data from the solver input file.

!!==============================================================================
!> @brief Initialization of system data from the solver input file.

module initiateSystemTypeModule

  use kindModule, only : dp

  implicit none

  real(dp), private, parameter :: epsT_p = 1.0e-12_dp !< Zero tolerance for time


contains

  !!============================================================================
  !> @brief Initializes system level model data from the command-line arguments.
  !>
  !> @param[out] sys System level model data
  !> @param[in] engines All general functions in the model
  !> @param[in] tIncEngineID Base ID of the function returning time step size
  !> @param[in] wFac Weighting factors for the DOFs of different kind,
  !> wFac(1)=translations, wFac(2)=rotation, wFac(3)=generalized DOFs
  !> @param[out] err Error flag
  !>
  !> @callergraph
  !>
  !> @brief Knut Morten Okstad
  !>
  !> @date 28 May 2002

  subroutine InitiateSystem (sys,engines,tIncEngineId,wFac,err)

    use SystemTypeModule      , only : SystemType, nullifySys, dp
    use FunctionTypeModule    , only : EngineType
    use IdTypeModule          , only : getId
    use FunctionTypeModule    , only : GetPtrToId
    use fileUtilitiesModule   , only : getDBGfile
    use dbgUnitsModule        , only : dbgSolve
#ifdef FT_DEBUG
    use progressModule        , only : lterm
#endif
    use reportErrorModule     , only : reportError, error_p, warning_p, note_p
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_isTrue
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getdouble
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getdoubles
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getint
#ifndef FT_DEBUG
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getints
#endif

    type(SystemType), intent(out) :: sys
    type(EngineType), intent(in)  :: engines(:)
    integer         , intent(in)  :: tIncEngineId
    real(dp)        , intent(in)  :: wFac(3)
    integer         , intent(out) :: err

    !! Local variables
    integer  :: NewmarkFlag, printSys(2)
    real(dp) :: alpha(2)

    !! --- Logic section ---

    err = 0
    call nullifySys (sys)

    !! Newmark time integration parameters
    call ffa_cmdlinearg_getint ('NewmarkFlag',NewmarkFlag)
    call ffa_cmdlinearg_getdoubles ('alphaNewmark',alpha,2)
    if (NewmarkFlag < 200) then ! Using the HHT-alpha method
       call checkRange ('alpha',alpha(1),0.0_dp,1.0_dp/3.0_dp)
       sys%alpha = -alpha(1) ! NB: Negative value
    else ! Using the generalized-alpha method
       call checkRange ('alpha_f',alpha(1),0.0_dp,0.5_dp)
       call checkRange ('alpha_m',alpha(2),3.0_dp*alpha(1)-1.0_dp,alpha(1))
       !! Notice that the variable sys%alpha_f and sys%alpha_m
       !! here are defined as 1.0 - (values from command-line option)
       sys%alpha_f = 1.0_dp - alpha(1)
       sys%alpha_m = 1.0_dp - alpha(2)
       sys%alpha = alpha(2) - alpha(1)
    end if
    sys%beta  = (1.0_dp-sys%alpha)*(1.0_dp-sys%alpha)/4.0_dp
    sys%gamma =  0.5_dp-sys%alpha

    !! General time stepping parameters
    call ffa_cmdlinearg_getdouble ('timeStart',sys%tStart)
    call ffa_cmdlinearg_getdouble ('timeEnd',sys%tEnd)
    call ffa_cmdlinearg_getdouble ('timeInc',sys%tInc)
    call ffa_cmdlinearg_getdouble ('maxInc',sys%maxInc)
    call ffa_cmdlinearg_getdouble ('minInc',sys%minInc)
    call ffa_cmdlinearg_getdouble ('quasiStatic',sys%tQStatic)
    call ffa_cmdlinearg_getdouble ('cutbackFactor',sys%cutbck(1))
    call ffa_cmdlinearg_getint ('cutbackSteps',sys%nCutStp)
    call ffa_cmdlinearg_getint ('autoTimeStep',sys%varInc)
    if (tIncEngineId > 0) then
       sys%tIncEngine => GetPtrToId(engines,tIncEngineId)
       if (associated(sys%tIncEngine)) then
          call reportError (note_p,'Time step size is governed by Engine'// &
               &            getId(sys%tIncEngine%id))
          if (sys%minInc < epsT_p) then
             sys%minInc = 1.0e-4_dp
             if (abs(sys%tEnd-sys%tStart) > epsT_p) then
                call reportError (note_p,'Minimum step size is set to 0.0001')
             end if
          end if
       else
          err = -1
          call reportError (error_p,'Invalid time step engine specification', &
               &            addString='InitiateSystem')
       end if
    else if (sys%tInc < 1.0e-6_dp) then
       sys%tInc = 1.0e-6_dp
       if (abs(sys%tEnd-sys%tStart) > epsT_p) then
          call reportError (warning_p,'The initial time step is too small,'// &
               &                      ' it is therefore reset to 1.0e-6')
       end if
    else if (sys%tInc > 0.8_dp*sys%maxInc .and. sys%varInc < 1) then
       sys%maxInc = 1.25_dp*sys%tInc
    end if

    !! Newton-Raphson iteration parameters
    call ffa_cmdlinearg_getbool ('stressStiffDyn',sys%stressStiffIsOn(1))
    call ffa_cmdlinearg_getint ('maxit',sys%maxIt)
    call ffa_cmdlinearg_getint ('minit',sys%minIt)
    call ffa_cmdlinearg_getint ('numit',sys%fixedIt)
    if (sys%fixedIt < -3) then
       sys%fxItEngine => GetPtrToId(engines,-sys%fixedIt)
       if (associated(sys%fxItEngine)) then
          call reportError (note_p,'Number of iterations is governed by '// &
               &            'Engine'//getId(sys%fxItEngine%id))
       else
          call reportError (warning_p,'Invalid iteration engine specification',&
               &            'Using convergence tolerances instead.')
       end if
       sys%fixedIt = 0
    end if

    !! Criteria for tangent matrix update during iterations
    call ffa_cmdlinearg_getint ('nupdat',sys%nUpdat)
    call ffa_cmdlinearg_getint ('maxSeqNoUpdate',sys%maxSequentialNoUpdate)
    call ffa_cmdlinearg_getdouble ('tolUpdateFactor',sys%tolUpdateFactor)

    !! Quasi-static equilibrium parameters
    if (ffa_cmdlinearg_isTrue('initEquilibrium')) then
       call ffa_cmdlinearg_getbool ('stressStiffEqu',sys%stressStiffIsOn(2))
       call ffa_cmdlinearg_getdouble ('tolInitEquil',sys%equTol)
       call ffa_cmdlinearg_getdouble ('limInitEquilStep',sys%equLim)
    end if

    !! Weight factors for norm calculation (see function ScaledNorm)
    sys%wDisp  = wFac ! The weights were read from the MECHANISM namelist
    sys%wForce = 0.0_dp
    if (wFac(1) > 0.0_dp) sys%wForce(1) = 1.0_dp/wFac(1)
    if (wFac(2) > 0.0_dp) sys%wForce(2) = 1.0_dp/wFac(2)
    if (wFac(3) > 0.0_dp) sys%wForce(3) = 1.0_dp/wFac(3)

#ifdef FT_DEBUG
    if (ffa_cmdlinearg_isTrue('consoleDebug')) then
       !! System debug print to console unit
       dbgSolve = lterm
    else
       printSys = 1
    end if
#else
    call ffa_cmdlinearg_getints ('printSys',printSys,2)
#endif
    if (printSys(2) > 0) then
       !! System debug print to separate output file
       dbgSolve = getDBGfile(13,'solveProcess.dbg')
    end if

  contains

    !> @brief Checks whether an alpha-parameter is in the valid range.
    subroutine checkRange (name,val,vmin,vmax)
      character(len=*), intent(in)    :: name
      real(dp)        , intent(inout) :: val
      real(dp)        , intent(in)    :: vmin,vmax
      character(len=128)  :: msg1, msg2
      real(dp)            :: oldVal
      real(dp), parameter :: eps_p = 1.0e-8_dp
      oldVal = val
      if (val < vmin-eps_p) then
         val = vmin
      else if (val > vmax+eps_p) then
         val = vmax
      else
         return
      end if
      write(msg1,601) name,oldVal,vmin,vmax
      write(msg2,602) name,val
      call reportError (warning_p,msg1,msg2)
601   format('The specified parameter ',A,' =',F8.4, &
           & ' is outside the stable range [',F8.4,',',F8.4,' ].')
602   format('Resetting ',A,' to ',F8.4)
    end subroutine checkRange

  end subroutine InitiateSystem


  !!============================================================================
  !> @brief Allocates the system vectors and Newton matrix.
  !>
  !> @param[in] sam Data for managing system matrix assembly
  !> @param sys System level model data
  !> @param[out] err Error flag
  !>
  !> @details This subroutine allocates all the system vectors (displacement,
  !> velocity, acceleration and force vectors) and the system Newton matrix.
  !> The dynamic convergence tolerances are also initialized.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 28 May 2003

  subroutine InitiateSystem2 (sam,sys,err)

    use SamModule             , only : SamType
    use SystemTypeModule      , only : SystemType, dp
    use SysMatrixTypeModule   , only : allocateSysMatrix
    use NormTypeModule        , only : initConvChecks
    use reportErrorModule     , only : reportError, debugFileOnly_p
    use reportErrorModule     , only : getErrorFile, allocationError
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_intValue
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_doubleValue

    type(SamType)   , intent(in)    :: sam
    type(SystemType), intent(inout) :: sys
    integer         , intent(out)   :: err

    !! --- Logic section ---

    allocate(sys%urd     (sam%ndof), &
         &   sys%urdd    (sam%ndof), &
         &   sys%urdp    (sam%ndof), &
         &   sys%urddp   (sam%ndof), &
         &   sys%dk      (sam%ndof), &
         &   sys%ak      (sam%ndof), &
         &   sys%sinc    (sam%ndof,3), &
         &   sys%rinc    (sam%ndof), &
         &   sys%RFk     (max(1,sam%mpar(17))), &
         &   sys%Qk      (max(1,sam%neq)), &
         &   sys%FSk     (max(1,sam%neq)), &
         &   sys%FDk     (max(1,sam%neq)), &
         &   sys%FIk     (max(1,sam%neq)), &
         &   sys%FIactual(max(1,sam%neq)), &
         &   sys%residual(max(1,sam%neq)), &
         &   sys%del     (max(1,sam%neq)), STAT=err)
    if (err /= 0) then
       err = AllocationError('InitiateSystem2')
       return
    end if

    sys%urd=0.0_dp; sys%urdd=0.0_dp
    sys%urdp=0.0_dp; sys%urddp=0.0_dp; sys%dk=0.0_dp; sys%ak=0.0_dp
    sys%sinc=0.0_dp; sys%rinc=0.0_dp; sys%RFk=0.0_dp
    sys%Qk=0.0_dp; sys%FSk=0.0_dp; sys%FDk=0.0_dp; sys%FIk=0.0_dp

    select case (mod(ffa_cmdlinearg_intValue('NewmarkFlag')/100,10)) ! iopAlg
    case (1:2)
       allocate(sys%Qprev(sam%neq),sys%FIprev(sam%neq),STAT=err)
       if (err /= 0) then
          err = AllocationError('InitiateSystem2')
          return
       end if
       sys%Qprev=0.0_dp; sys%FIprev=0.0_dp
    end select

    if (sam%neq > 0) then
       call allocateSysMatrix (sys%Nmat,err,getErrorFile(), &
            &                  'system Newton matrix')
       if (err /= 0) goto 915
    end if
    if (abs(sys%tEnd-sys%tStart) <= epsT_p) return ! No time incrementation

    !! Initialize the dynamic convergence tolerance parameters
    call initConvChecks (sys%convergenceSet, sys%maxit, &
         &               min(ffa_cmdlinearg_intValue('monitorWorst'),sam%neq), &
         &               ffa_cmdlinearg_intValue('monitorIter'), &
         &               ffa_cmdlinearg_doubleValue('tolVelProp'), err)
    if (err == 0) return

915 call reportError (debugFileOnly_p,'InitiateSystem2')

  end subroutine InitiateSystem2


  !!============================================================================
  !> @brief Initializes the time stepping before starting integration loop.
  !>
  !> @param sys System level model data
  !> @param[out] err Error flag
  !>
  !> @details In case a dynamic ramp-up is specified, the actual start time is
  !> shifted back the duration of the ramp-up, such that the start of the actual
  !> simulation will happen at the same time as without ramp-up. The number of
  !> ramp-up increments is used to define the time increment size during the
  !> ramp-up stage. The ramp function is a smooth trajectory function.
  !> @sa functiontypemodule::initrampfunction
  !>
  !> @callgraph @callergraph
  !>
  !> @brief Knut Morten Okstad
  !>
  !> @date 19 Aug 2022

  subroutine InitTimeStepping (sys,err)

    use SystemTypeModule        , only : SystemType, dp
    use FunctionTypeModule      , only : initRampFunction
    use InitiateSensorTypeModule, only : allocateTimeSensor
    use TimeStepModule          , only : getInitialTimeStepSize
    use reportErrorModule       , only : reportError, debugFileOnly_p
    use FFaCmdLineArgInterface  , only : ffa_cmdlinearg_getint
    use FFaCmdLineArgInterface  , only : ffa_cmdlinearg_getdoubles

    type(SystemType), intent(inout) :: sys
    integer         , intent(out)   :: err

    !! Local variables
    integer  :: nRampSteps
    real(dp) :: rampData(3), rampLength, rampStart, rampStop, velMaxRamp

    !! --- Logic section ---

    !! Set initial time and time increment size
    sys%time     = sys%tStart
    sys%timeStep = getInitialTimeStepSize(sys,err)
    if (err < 0) goto 915

    call ffa_cmdlinearg_getint ('rampSteps',nRampSteps)
    if (nRampSteps < 1) return ! No up-ramping

    !! We are using dynamic ramp-up of the time-dependent input functions
    call ffa_cmdlinearg_getdoubles ('rampData',rampData,3)
    velMaxRamp = rampData(1) ! Max normalized speed during ramp-up
    rampLength = rampData(2) ! Total length (in time) of ramp-up stage
    !! Define the start and stop time of the ramp-up stage
    rampStart = sys%tStart - (rampLength + max(rampData(3),0.0_dp))
    rampStop  = sys%tStart + sys%timeStep
    !! Create the ramp function (smooth trajectory)
    call initRampFunction (rampStart,rampStop,rampLength,velMaxRamp,err)
    if (err < 0) goto 915

    !! Set initial time and time increment size for ramp-up stage
    sys%time     = rampStart
    sys%timeStep = (sys%tStart - rampStart) / real(nRampSteps,dp)

    !! Ensure a time sensor exists for the ramp-up function
    call allocateTimeSensor (sys,err)
    if (err == 0) return

915 call reportError (debugFileOnly_p,'InitTimeStepping')

  end subroutine InitTimeStepping

end module initiateSystemTypeModule
