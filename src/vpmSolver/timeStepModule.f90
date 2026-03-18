!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file timeStepModule.f90
!>
!> @brief Subroutines/functions for automatic time stepping.

!!==============================================================================
!> @brief Module with subroutines/functions for automatic time stepping.
!>
!> @details The auto-time stepping is based in local trunction errors
!> calculated from the last three acceleration states. This is not much in use
!> any longer and is only kept for historical reasons. The module also handles
!> prescribed time step size defined by a general function (time step engine).

module TimeStepModule

  use KindModule, only : dp

  implicit none

  private

  real(dp), pointer    , save :: a1(:)  !< Bottom of acceleration stack
  real(dp), pointer    , save :: a2(:)  !< Middle of acceleration stack
  real(dp), pointer    , save :: a3(:)  !< Top of acceleration stack
  real(dp), allocatable, save :: err(:) !< Internal work array

  public :: pushAccelStack, getInitialTimeStepSize, getTimeStepSize
  public :: deallocateTimeStep


contains

  !!============================================================================
  !> @brief Pushes the given vector onto the acceleration stack.
  !>
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[in] acc Acceleration vector to be pushed
  !> @param[out] ierr Error flag
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 17 Jul 2000

  subroutine pushAccelStack (sam,acc,ierr)

    use SamModule        , only : SamType
    use reportErrorModule, only : allocationError

    type(SamType), intent(in)  :: sam
    real(dp)     , intent(in)  :: acc(:)
    integer      , intent(out) :: ierr

    !! Local variables
    integer           :: i, j
    real(dp), pointer :: a0(:)

    !! --- Logic section ---

    ierr = 0
    if (.not. allocated(err)) then
       allocate(a1(sam%neq),a2(sam%neq),a3(sam%neq),err(sam%neq),STAT=ierr)
       if (ierr /= 0) then
          ierr = allocationError('pushAccelStack')
          return
       end if
    end if

    a0 => a1
    do i = 1,sam%ndof
       j = sam%meqn(i)
       if (j > 0) a0(j) = acc(i)
    end do

    a1 => a2
    a2 => a3
    a3 => a0

  end subroutine pushAccelStack


  !!============================================================================
  !> @brief Calculates the next time step based on local truncation errors.
  !>
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[in] sys System level model data
  !> @param[in] errlim Truncation error limit
  !> @param[in] HP Previous time step size
  !> @param[in] H Current time step size
  !> @param[out] HX Next time step size
  !> @param[out] istat Status flag:
  !>   - = 0 : The time step in unchanged
  !>   - > 0 : The time step is increased or decreased
  !>   - < 0 : Error
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 17 Jul 2000

 subroutine getAutoStep (sam,sys,errlim,HP,H,HX,istat)

    use SamModule         , only : SamType
    use SystemTypeModule  , only : SystemType
    use NormRoutinesModule, only : ScaledNorm
    use reportErrorModule , only : reportError, debugFileOnly_p

    type(SamType)   , intent(in)  :: sam
    type(SystemType), intent(in)  :: sys
    real(dp)        , intent(in)  :: errlim, HP, H
    real(dp)        , intent(out) :: HX
    integer         , intent(out) :: istat

    !! Local variables
    real(dp), parameter :: EPS    = 1.0e-12_dp
    real(dp), parameter :: PFAC   = 0.85_dp
    real(dp), parameter :: RO     = 1.2_dp
    real(dp), parameter :: MAXINC = 2.0_dp

    integer  :: i
    real(dp) :: c1, c2, c3, error1, error2, errDis, errVel, dk, vk

    !! --- Logic Section ---

    HX = H
    istat = 0
    if (.not. allocated(err)) return

    !! Evaluate the two first terms of the displacement truncation error

    c1 = 0.5_dp * (sys%beta-1.0_dp/ 6.0_dp) * H*H
    c2 = 0.5_dp * (sys%beta-1.0_dp/12.0_dp) * H*H*H/(H+HP)
    c3 = c2 * H/HP

    do i = 1,sam%neq
       error1 = c1 * (a3(i)-a1(i))
       error2 = c2 * (a3(i)-a2(i)) - c3*(a2(i)-a1(i))
       err(i) = (abs(error1) + abs(error2))/(errlim*H)
    end do

    errDis = ScaledNorm(sam,err,sys%wDisp,doExpand=.true.,ierr=istat)
    if (istat < 0) goto 915

    !! Evaluate the two first terms of the velocity truncation error

    c1 = 0.5_dp * (sys%gamma-1.0_dp/2.0_dp) * H
    c2 = 0.5_dp * (sys%gamma-1.0_dp/3.0_dp) * H*H /(H+HP)
    c3 = c2 * h/hp

    do i = 1,sam%neq
       error1 = c1 * (a3(i)-a1(i))
       error2 = c2 * (a3(i)-a2(i)) - c3*(a2(i)-a1(i))
       err(i) = (abs(error1) + abs(error2))/(errlim*H)
    end do

    errVel = ScaledNorm(sam,err,sys%wDisp,doExpand=.true.,ierr=istat)
    if (istat < 0) goto 915

    !! Determine the time step correction, reduction or increase

    DK = (1.0_dp/errDis)**(1.0_dp/3.0_dp)
    if (abs(sys%gamma-0.5_dp) < eps) then
       VK = (1.0_dp/errVel)**(1.0_dp/3.0_dp)
    else
       VK = (1.0_dp/errVel)**(1.0_dp/2.0_dp)
    end if
    DK = MAX(DK,VK)

    if (errDis > RO .or. DK > 1.0_dp/PFAC) then

       istat = 1
       HX = MIN(PFAC*H*DK,MAXINC*H)

       if (HX < sys%minInc) HX = sys%minInc
       if (HX > sys%maxInc) HX = sys%maxInc

    end if

    return

915 continue
    call reportError (debugFileOnly_p,'getAutoStep')

  end subroutine getAutoStep


  !!============================================================================
  !> @brief Returns the length of the initial time step.
  !>
  !> @param sys System level model data
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 21 Aug 2007

  function getInitialTimeStepSize (sys,ierr) result(H0)

    use SystemTypeModule    , only : SystemType
    use EngineRoutinesModule, only : EngineValue
    use reportErrorModule   , only : reportError, debugFileOnly_p

    real(dp)                        :: H0
    type(SystemType), intent(inout) :: sys
    integer         , intent(out)   :: ierr

    !! --- Logic section ---

    ierr = 0
    if (associated(sys%tIncEngine)) then

       !! Prescribed time step length variation
       H0 = EngineValue(sys%tIncEngine,ierr)
       if (ierr < 0) then
          H0 = sys%minInc
          call reportError (debugFileOnly_p,'getInitialTimeStepSize')
       else if (H0 < sys%minInc) then
          !! Adjust the time step length if it is less than the minimum length
          H0 = sys%minInc
       end if

    else

       !! Fixed step length
       H0 = sys%tInc

    end if

  end function getInitialTimeStepSize


  !!============================================================================
  !> @brief Returns the length of the next time step.
  !>
  !> @param sys System level model data
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[in] ctrl Control system data
  !> @param[in] errlim Truncation error limit
  !> @param ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 26 Feb 2002

  function getTimeStepSize (sys,sam,ctrl,errlim,ierr) result(HX)

    use KindModule          , only : i8
    use SystemTypeModule    , only : SystemType
    use SamModule           , only : SamType
    use ControlTypeModule   , only : ControlType
    use EngineRoutinesModule, only : EngineValue
    use reportErrorModule   , only : reportError, debugFileOnly_p

    real(dp)                         :: HX
    type(SystemType) , intent(inout) :: sys
    type(SamType)    , intent(in)    :: sam
    type(ControlType), intent(in)    :: ctrl
    real(dp)         , intent(in)    :: errlim
    integer          , intent(inout) :: ierr

    !! Local variables
    integer :: lerr

    !! --- Logic section ---

    if (associated(sys%tIncEngine)) then

       !! Prescribed time step length variation
       lerr = 0
       HX = EngineValue(sys%tIncEngine,lerr)
       if (lerr < 0) goto 915

       !! Adjust the time step length if it is less than the minimum step length
       if (HX < sys%minInc) HX = sys%minInc

    else if (sys%varInc > 0 .and. sys%nStep > 2_i8) then

       !! Automatically adjusted time step based on local truncation error
       call getAutoStep (sam,sys,errlim,sys%timePrev,sys%timeStep,HX,lerr)
       if (lerr < 0) goto 915

       if (associated(ctrl%rreg)) then
          !! Time step governed by some control system limit?
          if (ctrl%rreg(2) > sys%time .and. ctrl%rreg(2) < sys%time+HX) then
             HX = ctrl%rreg(2) - sys%time
          end if
       end if

       !! Adjust the time step according to the max and min step length
       if (HX < sys%minInc) HX = sys%minInc
       if (HX > sys%maxInc) HX = sys%maxInc

    else

       !! Fixed step length
       HX = sys%tInc

    end if

    return

915 continue
    ierr = ierr - lerr
    call reportError (debugFileOnly_p,'getTimeStepSize')

  end function getTimeStepSize


  !!============================================================================
  !> @brief Deallocates the internal buffers for automatic time stepping.
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Jan 2017

  subroutine deallocateTimeStep ()

    !! --- Logic section ---

    if (associated(a1)) deallocate(a1)
    if (associated(a2)) deallocate(a2)
    if (associated(a3)) deallocate(a3)
    if (allocated(err)) deallocate(err)

    nullify(a1)
    nullify(a2)
    nullify(a3)

  end subroutine deallocateTimeStep

end module TimeStepModule
