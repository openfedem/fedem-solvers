!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file initiateFrictionTypeModule.f90
!> @brief Initialization of friction objects from the solver input file.

!!==============================================================================
!> @brief Initialization of friction objects from the solver input file.

module InitiateFrictionTypeModule

  implicit none

contains

  !!============================================================================
  !> @brief Initializes friction parameter sets with data from the input file.
  !>
  !> @param[in] infp File unit number for the solver input file
  !> @param[out] frictionData Array of all friction parameter objects
  !> @param[out] err Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date Sep 1998

  subroutine InitiateFrictionData (infp,frictionData,err)

    use IdTypeModule      , only : ldesc_p, initId, nullifyId, ReportInputError
    use FrictionTypeModule, only : FrictionParameterType, fricType_p, dp
    use FrictionTypeModule, only : DeallocateFrictionPrms
    use inputUtilities    , only : iuGetNumberOfEntries, iuSetPosAtNextEntry
    use inputUtilities    , only : iuCharToInt
    use progressModule    , only : lterm
    use reportErrorModule , only : AllocationError

    integer                    , intent(in)  :: infp
    type(FrictionParameterType), pointer     :: frictionData(:)
    integer                    , intent(out) :: err

    !! Local variables
    integer :: idIn, nSets, stat

    !! Define the FRICTION_SET namelist
    character(len=30)  :: type
    character(ldesc_p) :: extDescr
    integer  :: id, extId(10)
    real(dp) :: typeDepParams(10), CoulombCoeff, StribeckMagn, StribeckSpeed, &
         &      ViscCoeff, PrestressLoad, AsymMagn, &
         &      FricAmpl, FricFreq, FricPhase, StickStiffness
    namelist /FRICTION_SET/ id, extId, extDescr, type, typeDepParams, &
         &                  CoulombCoeff, StribeckMagn, StribeckSpeed, &
         &                  ViscCoeff, PrestressLoad, AsymMagn, &
         &                  FricAmpl, FricFreq, FricPhase, StickStiffness

    !! --- Logic section ---

    nSets = iuGetNumberOfEntries(infp,'&FRICTION_SET',err)
    if (err /= 0) return
    write(lterm,*) 'Number of &FRICTION_SET =', nSets

    call DeallocateFrictionPrms (frictionData)
    allocate(frictionData(nSets),STAT=stat)
    if (stat /= 0) then
       err = AllocationError('InitiateFrictionData 10')
       return
    end if

    do idIn = 1, nSets
       call nullifyId (frictionData(idIn)%id)
       if (.not. iuSetPosAtNextEntry(infp,'&FRICTION_SET')) then
          err = err - 1
          call ReportInputError ('FRICTION_SET',idIn)
          cycle
       end if

       !! Default values
       id = 0; extId = 0; extDescr = ''; type = ''; typeDepParams = 0.0_dp
       CoulombCoeff = 0.0_dp; StribeckMagn = 0.0_dp; StribeckSpeed = 0.0_dp
       ViscCoeff = 0.0_dp; PrestressLoad = 0.0_dp; AsymMagn = 0.0_dp
       FricAmpl = 0.0_dp; FricFreq = 0.0_dp; FricPhase = 0.0_dp
       StickStiffness = 0.0_dp

       read(infp,nml=FRICTION_SET,iostat=stat)
       if (stat /= 0) then
          err = err - 1
          call ReportInputError ('FRICTION_SET',idIn)
          cycle
       end if

       call initId (frictionData(idIn)%id,id,extId,extDescr,stat)
       frictionData(idIn)%type = iuCharToInt(type,fricType_p,stat)
       if (stat < 0) then
          err = err - 1
          call ReportInputError ('FRICTION_SET',idIn,frictionData(idIn)%id)
          cycle
       end if

       frictionData(idIn)%typeDepParams = typeDepParams(1:3)
       frictionData(idIn)%CoulombCoeff  = CoulombCoeff
       frictionData(idIn)%StribeckMagn  = StribeckMagn
       frictionData(idIn)%StribeckSpeed = StribeckSpeed
       frictionData(idIn)%ViscCoeff     = ViscCoeff
       frictionData(idIn)%PrestressLoad = PrestressLoad
       frictionData(idIn)%AsymMagn      = AsymMagn
       frictionData(idIn)%FricAmpl      = FricAmpl
       frictionData(idIn)%FricFreq      = FricFreq
       frictionData(idIn)%FricPhase     = FricPhase

    end do

  end subroutine InitiateFrictionData


  !!============================================================================
  !> @brief Creates an array of pointers to all frictions.
  !>
  !> @param[in] joints Array of all joints in the model
  !> @param[in] cElems Array of all contact elements in the model
  !> @param[out] frictions Array of pointers to all friction objects
  !> @param[out] err Error flag
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date Sep 2000

  subroutine InitFrictionPtrArray (joints,cElems,frictions,err)

    use MasterSlaveJointTypeModule, only : MasterSlaveJointType
    use ContactElementTypeModule  , only : ContactElementType
    use FrictionTypeModule        , only : FrictionPtrType
    use reportErrorModule         , only : AllocationError

    type(MasterSlaveJointType), intent(in), target :: joints(:)
    type(ContactElementType)  , intent(in), target :: cElems(:)
    type(FrictionPtrType)     , pointer            :: frictions(:)
    integer                   , intent(out)        :: err

    !! Local variables
    integer :: i, j, nFric

    !! --- Logic section ---

    nFric = 0
    do i = 1, size(joints)
       do j = 1, joints(i)%nJointDOFs
          if (associated(joints(i)%jointDofs(j)%friction)) nFric = nFric + 1
       end do
    end do

    do i = 1, size(cElems)
       if (associated(cElems(i)%friction)) nFric = nFric + 1
    end do

    deallocate(frictions)
    allocate(frictions(nFric),STAT=err)
    if (err /= 0) then
       err = AllocationError('InitFrictionPtrArray')
       return
    else if (nFric == 0) then
       return
    end if

    nFric = 0
    do i = 1, size(joints)
       do j = 1, joints(i)%nJointDOFs
          if (associated(joints(i)%jointDofs(j)%friction)) then
             nFric = nFric + 1
             frictions(nFric)%p => joints(i)%jointDofs(j)%friction
          end if
       end do
    end do

    do i = 1, size(cElems)
       if (associated(cElems(i)%friction)) then
          nFric = nFric + 1
          frictions(nFric)%p => cElems(i)%friction
       end if
    end do

  end subroutine InitFrictionPtrArray

end module InitiateFrictionTypeModule
