!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file controlRoutinesModule.f90
!> @brief Subroutine for invocation of the regulation system solver.

!!==============================================================================
!> @brief Module with a subroutine that wraps the regulation system solver.

module ControlRoutinesModule

  implicit none

  private :: SetControlInput

contains

  !!============================================================================
  !> @brief Sets input values for control variables from engines and sensors.

  subroutine SetControlInput (ctrl,ierr)

    use ControlTypeModule   , only : ControlType
    use EngineRoutinesModule, only : EngineValue
    use IdTypeModule        , only : getId
    use reportErrorModule   , only : reportError, error_p, debugFileOnly_p

    type(ControlType), intent(inout) :: ctrl
    integer          , intent(out)   :: ierr

    !! Local variables
    integer :: i, j, lerr

    !! --- Logic section ---

    ierr = 0
    do i = 1, size(ctrl%input)
       j = ctrl%input(i)%var

       lerr = ierr
       if (associated(ctrl%input(i)%engine)) then
          ctrl%vreg(j) = EngineValue(ctrl%input(i)%engine,ierr)
       end if

       if (ierr < lerr) then
          call reportError (error_p,'Unable to set input value '// &
               &            'for control variable'//getId(ctrl%vregId(j)))
       end if

    end do

    if (ierr < 0) call reportError (debugFileOnly_p,'SetControlInput')

  end subroutine SetControlInput


  !!============================================================================
  !> @brief Iterates the control module.

  subroutine IterateControlSystem (sys,ctrl,MODE,MSIM,IERR)

    use ControlTypeModule, only : ControlType
    use SystemTypeModule , only : SystemType
    use IdTypeModule     , only : StrId
    use kindModule       , only : dp
    use profilerModule   , only : startTimer, stopTimer, ctrl_p
    use allocationModule , only : reAllocate
    use reportErrorModule, only : reportError, error_p, warning_p
    use reportErrorModule, only : noteFileOnly_p, debugFileOnly_p

    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getint
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getbool
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getdoubles

    type(SystemType) , intent(in)    :: sys
    type(ControlType), intent(inout) :: ctrl
    integer          , intent(in)    :: MODE, MSIM(:)
    integer          , intent(out)   :: IERR

    !! Local variables
    integer  :: IPRINT, IOP, NDDIM, ISTAT(10), IAERR(2)
    logical  :: zCtrl0, noDcSplit
    real(dp) :: tolFact(3)

    !! --- Logic section ---

    IERR = 0
    if (maxval(MSIM(45:48)) < 1) return ! No control system at all

    call startTimer (ctrl_p)

    !! Establish the VREG vector with sensor and engine input
    call SetControlInput (ctrl,ierr)
    if (ierr < 0) goto 900

    if (MSIM(48) < 1) goto 900 ! No control elements

    call ffa_cmdlinearg_getint ('debug',iprint)
    call ffa_cmdlinearg_getbool ('zeroControlInit',zCtrl0)
    call ffa_cmdlinearg_getbool ('noDiscontinuitySplit',noDcSplit)
    call ffa_cmdlinearg_getdoubles ('tolFactorize',tolFact,3)
    IOP   = MODE
    IAERR = 0
    ISTAT = 0
100 continue
    call RSFED (IOP, sys%time+sys%timeStep, sys%timeStep, MSIM(47), MSIM(48), &
         &      ctrl%vreg(1), ctrl%mpireg(1), ctrl%ireg(1), &
         &      ctrl%mprreg(1), ctrl%rreg(1), ctrl%delay(1), size(ctrl%delay), &
         &      IAERR,ISTAT,zCtrl0,.not.noDcSplit,tolFact(3),iprint)

    if (IAERR(1) == 999 .and. IOP > 0) then
       IOP = -MODE
       NDDIM = int(1.2*IAERR(2))
       call reportError (noteFileOnly_p,'DELAY buffer size increased from '// &
            &            StrId(size(ctrl%delay)) //' to '// StrId(NDDIM))
       call reAllocate ('IterateControlSystem',ctrl%delay,NDDIM,ierr,.true.)
       if (ierr == 0) goto 100
    else if (IAERR(1) == -7) then
       call OUTMES (sys%time+sys%timeStep,IAERR,ISTAT)
       call reportError (warning_p, &
            &            'Difficulty with discontinuity in control system')
    else if (IAERR(1) < 0) then
       ierr = IAERR(1)
       call OUTMES (sys%time+sys%timeStep,IAERR,ISTAT)
       call reportError (error_p,'Failure iterating control system')
    end if

900 if (ierr < 0) call reportError (debugFileOnly_p,'IterateControlSystem')
    call stopTimer (ctrl_p)

  end subroutine IterateControlSystem

end module ControlRoutinesModule
