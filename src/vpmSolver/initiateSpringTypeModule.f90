!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file initiateSpringTypeModule.f90
!> @brief Initialization of spring objects from the solver input file.

!!==============================================================================
!> @brief Initialization of spring objects from the solver input file.

module InitiateSpringTypeModule

  implicit none

contains

  !!============================================================================
  !> @brief Initializes spring failure objects with data from the input file.
  !>
  !> @param[in] infp File unit number for the solver input file
  !> @param[out] springFailures Array of all spring failure objects in the model
  !> @param[out] err Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date Apr 2005

  subroutine InitiateSpringFailures (infp,springFailures,err)

    use kindModule       , only : hugeVal_p, dp
    use IdTypeModule     , only : ldesc_p, initId, ReportInputError
    use SpringTypeModule , only : SpringFailureType
    use SpringTypeModule , only : DeallocateFailures, NullifyFailure
    use inputUtilities   , only : iuGetNumberOfEntries, iuSetPosAtNextEntry
    use progressModule   , only : lterm
    use reportErrorModule, only : AllocationError, reportError, debugFileOnly_p

    integer                , intent(in)  :: infp
    type(SpringFailureType), pointer     :: springFailures(:)
    integer                , intent(out) :: err

    !! Local Variables
    integer :: idIn, nSpringFailures, stat

    !! Define the SPRING_FAILURE namelist
    character(ldesc_p) :: extDescr
    integer  :: id, extId(10)
    logical  :: compFailure
    real(dp) :: deflectionMax, deflectionMin, forceMax, forceMin
    namelist /SPRING_FAILURE/ id, extId, extDescr, compFailure, &
         &                    deflectionMax, deflectionMin, forceMax, forceMin

    !! --- Logic section ---

    nSpringFailures = iuGetNumberOfEntries(infp,'&SPRING_FAILURE',err)
    if (err /= 0) return
    write(lterm,*) 'Number of &SPRING_FAILURE =',nSpringFailures

    call DeallocateFailures (springFailures)
    allocate(springFailures(nSpringFailures),STAT=stat)
    if (stat /= 0) then
       err = AllocationError('InitiateSpringFailures')
       return
    end if

    do idIn = 1, nSpringFailures
       call NullifyFailure (springFailures(idIn))

       if (.not. iuSetPosAtNextEntry(infp,'&SPRING_FAILURE')) then
          err = err - 1
          call ReportInputError ('SPRING_FAILURE',idIn)
          cycle
       end if

       !! Default values
       id=0; extId=0; extDescr=''
       compFailure=.false.
       deflectionMax=-hugeVal_p; deflectionMin=hugeVal_p
       forceMax=-hugeVal_p; forceMin=hugeVal_p

       read(infp,nml=SPRING_FAILURE,iostat=stat)
       if (stat /= 0) then
          err = err - 1
          call ReportInputError ('SPRING_FAILURE',idIn)
          cycle
       end if

       call initId (springFailures(idIn)%id,id,extId,extDescr,stat)

       springFailures(idIn)%compFailure = compFailure

       springFailures(idIn)%deflRange(1)  = deflectionMin
       springFailures(idIn)%deflRange(2)  = deflectionMax

       springFailures(idIn)%deflActive(1) = deflectionMin <  hugeVal_p
       springFailures(idIn)%deflActive(2) = deflectionMax > -hugeVal_p

       springFailures(idIn)%forceRange(1)  = forceMin
       springFailures(idIn)%forceRange(2)  = forceMax

       springFailures(idIn)%forceActive(1) = forceMin <  hugeVal_p
       springFailures(idIn)%forceActive(2) = forceMax > -hugeVal_p

       if ( all(springFailures(idIn)%deflActive) .and. &
            deflectionMax <= deflectionMin ) then
          err = err - 1
          call ReportInputError ('SPRING_FAILURE',idIn,springFailures(idIn)%id,&
               &                 msg='deflectionMax <= deflectionMin')
       end if
       if ( all(springFailures(idIn)%forceActive) .and. &
            forceMax <= forceMin ) then
          err = err - 1
          call ReportInputError ('SPRING_FAILURE',idIn,springFailures(idIn)%id,&
               &                 msg='forceMax <= forceMin')
       end if

    end do

    if (err < 0) call reportError (debugFileOnly_p,'InitiateSpringFailures')

  end subroutine InitiateSpringFailures


  !!============================================================================
  !> @brief Initializes spring yield objects with data from the input file.
  !>
  !> @param[in] infp File unit number for the solver input file
  !> @param[in] engines Array of all engines in the model
  !> @param[out] springYields Array of all spring yield objects in the model
  !> @param[out] err Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date Apr 2005

  subroutine InitiateSpringYields (infp,engines,springYields,err)

    use kindModule        , only : hugeVal_p, dp
    use IdTypeModule      , only : ldesc_p, initId, StrId, ReportInputError
    use FunctionTypeModule, only : EngineType, GetPtrToId
    use SpringTypeModule  , only : SpringYieldType
    use SpringTypeModule  , only : DeallocateYields, NullifyYield
    use inputUtilities    , only : iuGetNumberOfEntries, iuSetPosAtNextEntry
    use progressModule    , only : lterm
    use reportErrorModule , only : AllocationError, reportError, debugFileOnly_p

    integer              , intent(in)  :: infp
    type(EngineType)     , target      :: engines(:)
    type(SpringYieldType), pointer     :: springYields(:)
    integer              , intent(out) :: err

    !! Local Variables
    integer :: idIn, nSpringYields, stat

    !! Define the SPRING_YIELD namelist
    character(ldesc_p) :: extDescr
    integer  :: id, extId(10)
    integer  :: yieldForceMaxEngine, yieldForceMinEngine
    real(dp) :: yieldForceMax, yieldForceMin, yieldDeflectionAbsMax
    namelist /SPRING_YIELD/ id, extId, extDescr, &
         &                  yieldForceMaxEngine, yieldForceMinEngine, &
         &                  yieldForceMax, yieldForceMin, yieldDeflectionAbsMax

    !! --- Logic section ---

    nSpringYields = iuGetNumberOfEntries(infp,'&SPRING_YIELD',err)
    if (err /= 0) return
    write(lterm,*) 'Number of &SPRING_YIELD =',nSpringYields

    call DeallocateYields (springYields)
    allocate(springYields(nSpringYields),STAT=stat)
    if (stat /= 0) then
       err = AllocationError('InitiateSpringYields')
       return
    end if

    do idIn = 1, nSpringYields
       call NullifyYield (springYields(idIn))

       if (.not. iuSetPosAtNextEntry(infp,'&SPRING_YIELD')) then
          err = err - 1
          call ReportInputError ('SPRING_YIELD',idIn)
          cycle
       end if

       !! Default values
       id=0; extId=0; extDescr=''
       yieldForceMaxEngine=0; yieldForceMinEngine=0
       yieldForceMax=-hugeVal_p; yieldForceMin=hugeVal_p
       yieldDeflectionAbsMax=0.0_dp

       read(infp,nml=SPRING_YIELD,iostat=stat)
       if (stat /= 0) then
          err = err - 1
          call ReportInputError ('SPRING_YIELD',idIn)
          cycle
       end if

       call initId (springYields(idIn)%id,id,extId,extDescr,stat)

       if (yieldForceMaxEngine > 0) then
          springYields(idIn)%pos%active = .true.
          springYields(idIn)%pos%ysf(2) = 1.0_dp
          springYields(idIn)%pos%engine => GetPtrToId(engines, &
               &                                      yieldForceMaxEngine)
          if (.not. associated(springYields(idIn)%pos%engine)) then
             err = err - 1
             call ReportInputError ('SPRING_YIELD',idIn,springYields(idIn)%id, &
                  msg='Non-existing engine, id ='//StrId(yieldForceMaxEngine))
          end if
       else if (yieldForceMax > -hugeVal_p) then
          springYields(idIn)%pos%active = .true.
          springYields(idIn)%pos%ysf(1) = yieldForceMax
       end if

       if (yieldForceMinEngine > 0) then
          springYields(idIn)%neg%active = .true.
          springYields(idIn)%neg%ysf(2) = 1.0_dp
          springYields(idIn)%neg%engine => GetPtrToId(engines, &
               &                                      yieldForceMinEngine)
          if (.not. associated(springYields(idIn)%neg%engine)) then
             err = err - 1
             call ReportInputError ('SPRING_YIELD',idIn,springYields(idIn)%id, &
                  msg='Non-existing engine, id ='//StrId(yieldForceMinEngine))
          end if
       else if (yieldForceMin < hugeVal_p) then
          springYields(idIn)%neg%active = .true.
          springYields(idIn)%neg%ysf(1) = yieldForceMin
       end if

       if ( springYields(idIn)%pos%active .and. &
            springYields(idIn)%neg%active .and. &
            yieldForceMaxEngine <= 0 .and. yieldForceMinEngine <= 0 .and. &
            yieldForceMax <= yieldForceMin ) then
          err = err - 1
          call ReportInputError ('SPRING_YIELD',idIn,springYields(idIn)%id, &
               &                 msg='yieldForceMax <= yieldForceMin')
       end if

       springYields(idIn)%deflectionMax    = yieldDeflectionAbsMax
       springYields(idIn)%hasDeflectionMax = yieldDeflectionAbsMax > 0.0_dp

    end do

    if (err < 0) call reportError (debugFileOnly_p,'InitiateSpringYields')

  end subroutine InitiateSpringYields


  !!============================================================================
  !> @brief Initializes spring elements with data from the solver input file.
  !>
  !> @param[in] infp File unit number for the solver input file
  !> @param[in] triads Array of all triads in the model
  !> @param[in] baseSprings Array of all spring base objects in the model
  !> @param[out] springs Array of all spring elements in the model
  !> @param[out] err Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Karl Erik Thoresen and Bjorn Haugen
  !> @date Nov 1998
  !>
  !> @author Knut Morten Okstad
  !> @date Jul 2002

  subroutine InitiateSprings (infp,triads,baseSprings,springs,err)

    use IdTypeModule      , only : ldesc_p
    use IdTypeModule      , only : initId, getId, StrId, ReportInputError
    use TriadTypeModule   , only : dp, TriadType, allocateNodeForce, GetPtrToId
    use SpringTypeModule  , only : SpringBaseType, SpringType
    use SpringTypeModule  , only : GetPtrToId, DeallocateSprings, NullifySpring
    use inputUtilities    , only : iuGetNumberOfEntries, iuSetPosAtNextEntry
    use progressModule    , only : lterm
    use reportErrorModule , only : AllocationError, reportError
    use reportErrorModule , only : warning_p, debugFileOnly_p

    integer             , intent(in)  :: infp
    type(TriadType)     , intent(in)  :: triads(:)
    type(SpringBaseType), intent(in)  :: baseSprings(:)
    type(SpringType)    , pointer     :: springs(:)
    integer             , intent(out) :: err

    !! Local Variables
    integer, parameter :: maxTriads_p = 10
    integer  :: i, idIn, nSprings, nSpr, nTriad, stat
    real(dp) :: fDir(3), sLen(6)

    !! Define the SPRING_ELEMENT namelist
    character(ldesc_p) :: extDescr
    integer  :: id, extId(10), triadIDs(maxTriads_p), springBaseId(7)
    real(dp) :: alpha2, couplStiff(15)
    namelist /SPRING_ELEMENT/ id, extId, extDescr, triadIDs, springBaseId, &
         &                    alpha2, couplStiff

    !! --- Logic section ---

    nSprings = iuGetNumberOfEntries(infp,'&SPRING_ELEMENT',err)
    if (err /= 0) return
    write(lterm,*) 'Number of &SPRING_ELEMENT =',nSprings

    call DeallocateSprings (springs)
    allocate(springs(nSprings),STAT=stat)
    if (stat /= 0) then
       err = AllocationError('InitiateSprings 10')
       return
    end if

    do idIn = 1, nSprings

       call NullifySpring (springs(idIn))
       if (.not. iuSetPosAtNextEntry(infp,'&SPRING_ELEMENT')) then
          err = err - 1
          call ReportInputError ('SPRING_ELEMENT',idIn)
          cycle
       end if

       !! Default values
       id=0; extId=0; extDescr=''; triadIDs=0; springBaseId=0
       alpha2=0.0_dp; couplStiff=0.0_dp
       read(infp,nml=SPRING_ELEMENT,iostat=stat)
       if (stat /= 0) then
          err = err - 1
          call ReportInputError ('SPRING_ELEMENT',idIn)
          cycle
       end if

       allocate(springs(idIn)%id,STAT=stat)
       if (stat /= 0) then
          err = AllocationError('InitiateSprings 11')
          return
       end if

       call initId (springs(idIn)%id,id,extId,extDescr,stat)

       nSpr = 6
       do while (nSpr > 0)
          if (springBaseID(nSpr) > 0) exit
          nSpr = nSpr - 1
       end do
       if (nSpr == 0) then
          err = err - 1
          call ReportInputError ('SPRING_ELEMENT',idIn,springs(idIn)%id, &
               msg='No base springs specified')
          cycle
       end if

       nTriad = count(triadIDs > 0)
       if (nTriad == 0) then
          err = err - 1
          call ReportInputError ('SPRING_ELEMENT',idIn,springs(idIn)%id, &
               msg='No triads specified')
          cycle
       end if

       allocate(springs(idIn)%spr(nSpr),springs(idIn)%triads(nTriad),STAT=stat)
       if (stat /= 0) then
          err = AllocationError('InitiateSprings 15')
          return
       end if

       stat = err
       do i = 1, nSpr
          if (springBaseId(i) > 0) then
             springs(idIn)%spr(i)%p => GetPtrToId(baseSprings,springBaseId(i))
          else
             nullify(springs(idIn)%spr(i)%p)
             cycle
          end if
          if (.not. associated(springs(idIn)%spr(i)%p)) then
             err = err - 1
             call ReportInputError ('SPRING_ELEMENT',idIn,springs(idIn)%id, &
                  msg='Non-existing base spring, id ='//StrId(springBaseId(i)))
          else if (springs(idIn)%spr(i)%p%dof >= 0) then
             err = err - 1
             call ReportInputError ('SPRING_ELEMENT',idIn,springs(idIn)%id, &
                  msg='Reference to a base spring that already is in use')
          end if
       end do

       do i = 1, nTriad
          springs(idIn)%triads(i)%p => GetPtrToId(triads,triadIDs(i))
          if (associated(springs(idIn)%triads(i)%p)) then
             springs(idIn)%nDOFs = springs(idIn)%nDOFs + &
                  &                springs(idIn)%triads(i)%p%nDOFs
             !! Set samNodNum temporarily to flag that the triads should receive
             !! proper node numbers later, also when the triad is grounded
             springs(idIn)%triads(i)%p%samNodNum = 1
          else
             err = err - 1
             call ReportInputError ('SPRING_ELEMENT',idIn,springs(idIn)%id, &
                  msg='Non-existing triad, id ='//StrId(triadIDs(i)))
          end if
       end do
       if (err < stat) cycle

       if (nTriad == 1) then

          !! Spring to ground
          sLen(1:3) = springs(idIn)%triads(1)%p%ur(:,4)
          call FFa_glbEulerZYX (springs(idIn)%triads(1)%p%ur(:,1:3),sLen(4:6))

       else if (nTriad == 2 .and. springBaseId(7) > 0) then

          !! Global spring
          sLen(1:3) = springs(idIn)%triads(2)%p%ur(:,4) &
               &    - springs(idIn)%triads(1)%p%ur(:,4)
          call FFa_eulerZYX (springs(idIn)%triads(1)%p%ur(:,1:3), &
               &             springs(idIn)%triads(2)%p%ur(:,1:3),sLen(4:6))

       else if (nSpr == 1 .and. associated(springs(idIn)%spr(1)%p)) then

          !! Axial spring
          springs(idIn)%alpha2 = alpha2

          allocate(springs(idIn)%spr(1)%p%length, &
               &   springs(idIn)%forceDir(4,nTriad-1),STAT=stat)
          if (stat /= 0) then
             err = AllocationError('InitiateSprings 20')
             return
          end if

          slen = 0.0_dp
          do i = 2, nTriad
             fDir = springs(idIn)%triads(i)%p%ur(:,4) &
                  - springs(idIn)%triads(i-1)%p%ur(:,4)
             slen(1) = slen(1) + sqrt(dot_product(fDir,fDir))
          end do
          springs(idIn)%spr(1)%p%allocatedLength = .true.
          springs(idIn)%spr(1)%p%length = slen(1)
          springs(idIn)%spr(1)%p%dof = 0
          springs(idIn)%forceDir = 0.0_dp
          nSpr = -1

       else

          err = err - 1
          call ReportInputError ('SPRING_ELEMENT',idIn,springs(idIn)%id)
          cycle

       end if
       if (nSpr > 0) then

          !! Initializations for global springs and springs to ground
          do i = 1, nSpr
             if (.not. associated(springs(idIn)%spr(i)%p)) cycle

             allocate(springs(idIn)%spr(i)%p%length,STAT=stat)
             if (stat /= 0) then
                err = AllocationError('InitiateSprings 30')
                return
             end if

             springs(idIn)%spr(i)%p%allocatedLength = .true.
             springs(idIn)%spr(i)%p%dof = i
             springs(idIn)%spr(i)%p%length = sLen(i)
             if (triadIDs(2) == 0) then
                springs(idIn)%spr(i)%p%l0 = sLen(i) !hack
             end if

          end do
          do i = 15, 1, -1
             if (abs(couplStiff(i)) <= 1.0e-16_dp) cycle

             !! We have some explicit coupling terms for this spring
             allocate(springs(idIn)%couplStiff(i),STAT=stat)
             if (stat /= 0) then
                err = AllocationError('InitiateSprings 35')
                return
             end if
             springs(idIn)%couplStiff = couplStiff(1:i)
             exit

          end do

       end if

       !! If the spring triads also are connected to a superelement,
       !! allocate separate arrays for storage of the total nodal force
       stat = err
       do i = 1, nTriad
          call allocateNodeForce (springs(idIn)%triads(i)%p,err)
          if (err < stat) exit
       end do

       !! Consistency check: Deactivate all completely grounded springs
       if (springs(idIn)%nDOFs < 1) then
          do i = 1, abs(nSpr)
             if (associated(springs(idIn)%spr(i)%p)) then
                springs(idIn)%spr(i)%p%isActive = .false.
             end if
          end do
          call reportError (warning_p, &
               'Spring element'//trim(getId(springs(idIn)%id))// &
               ' is only connected to grounded Triads (ignored).')
       end if

    end do

    if (err < 0) call reportError (debugFileOnly_p,'InitiateSprings')

  end subroutine InitiateSprings


  !!============================================================================
  !> @brief Initializes spring base objects with data from the input file.
  !>
  !> @param[in] infp File unit number for the solver input file
  !> @param[in] engines Array of all engines in the model
  !> @param[in] functions Array of all functions in the model
  !> @param[in] springFailures Array of all spring failure objects in the model
  !> @param[in] springYields Array of all spring yield objects in the model
  !> @param[out] springs Array of all spring base objects in the model
  !> @param[out] err Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !> @date Apr 2002
  !>
  !> @author Knut Morten Okstad
  !> @date 1 Jul 2002

  subroutine InitiateBaseSprings (infp,engines,functions,springFailures, &
       &                          springYields,springs,err)

    use IdTypeModule          , only : ldesc_p
    use IdTypeModule          , only : initId, getId, StrId, ReportInputError
    use FunctionTypeModule    , only : EngineType, FunctionType, dp
    use FunctionTypeModule    , only : FunctionDerivative
    use SpringTypeModule      , only : SpringFailureType, SpringYieldType
    use SpringTypeModule      , only : SpringBaseType, GetPtrToId
    use SpringTypeModule      , only : DeallocateSprings, NullifySpring
    use inputUtilities        , only : iuGetNumberOfEntries, iuSetPosAtNextEntry
    use progressModule        , only : lterm
    use reportErrorModule     , only : allocationError, getErrorFile
    use reportErrorModule     , only : reportError, debugFileOnly_p
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getbool

    integer                , intent(in)  :: infp
    type(EngineType)       , intent(in)  :: engines(:)
    type(FunctionType)     , intent(in)  :: functions(:)
    type(SpringFailureType), intent(in)  :: springFailures(:)
    type(SpringYieldType)  , intent(in)  :: springYields(:)
    type(SpringBaseType)   , pointer     :: springs(:)
    integer                , intent(out) :: err

    !! Local Variables
    integer :: i, idIn, nSprings, stat
    logical :: allSec, allSpring, allStiff
    logical :: allDef, allLength, allForce, allEnergy

    !! Define the SPRING_BASE namelist
    character(ldesc_p) :: extDescr
    integer  :: id, extId(10), lengthEngineId, forceFuncId
    integer  :: stiffFuncId, stiffScaleEngineId
    integer  :: stiffScaleEnginePosId, stiffScaleEngineNegId, saveVar(5)
    integer  :: springFailureId, springYieldId, unLoadType
    real(dp) :: l0, l1, s0, s1
    namelist /SPRING_BASE/ id, extId, extDescr, &
         &                 lengthEngineId, stiffFuncId, forceFuncId, &
         &                 l0, l1, s0, s1, stiffScaleEngineId, &
         &                 stiffScaleEnginePosId, stiffScaleEngineNegId, &
         &                 springFailureId, springYieldId, unLoadType, saveVar

    !! --- Logic section ---

    nSprings = iuGetNumberOfEntries(infp,'&SPRING_BASE',err)
    if (err /= 0) return
    write(lterm,*) 'Number of &SPRING_BASE =',nSprings

    call DeallocateSprings (springs)
    allocate(springs(nSprings),STAT=stat)
    if (stat /= 0) then
       err = AllocationError('initiateSpringBaseType 10')
       return
    end if

    call ffa_cmdlinearg_getbool ('allSecondaryVars',allSec)
    call ffa_cmdlinearg_getbool ('allSpringVars',allSpring)
    call ffa_cmdlinearg_getbool ('allStiffVars',allStiff)
    call ffa_cmdlinearg_getbool ('allLengthVars',allLength)
    call ffa_cmdlinearg_getbool ('allDefVars',allDef)
    call ffa_cmdlinearg_getbool ('allForceVars',allForce)
    call ffa_cmdlinearg_getbool ('allEnergyVars',allEnergy)

    do idIn = 1, nSprings

       call NullifySpring (springs(idIn))
       if (.not. iuSetPosAtNextEntry(infp,'&SPRING_BASE')) then
          err = err - 1
          call ReportInputError ('SPRING_BASE',idIn)
          cycle
       end if

       !! Default values
       id=0; extId=0; extDescr=''
       lengthEngineId=0; stiffFuncId=0; forceFuncId=0
       stiffScaleEngineId=0; stiffScaleEnginePosId=0; stiffScaleEngineNegId=0
       springFailureId=0; springYieldId=0; unLoadType=0
       l0=0.0_dp; l1=0.0_dp; s0=0.0_dp; s1=0.0_dp; saveVar=0

       read(infp,nml=SPRING_BASE,iostat=stat)
       if (stat /= 0) then
          err = err - 1
          call ReportInputError ('SPRING_BASE',idIn)
          cycle
       end if

       call initId (springs(idIn)%id,id,extId,extDescr,stat)
       springs(idIn)%l0 = l0
       springs(idIn)%l1 = l1
       springs(idIn)%s0 = s0
       springs(idIn)%s1 = s1

       if (stiffScaleEnginePosId < 1) stiffScaleEnginePosId = stiffScaleEngineId
       if (stiffScaleEngineNegId < 1) stiffScaleEngineNegId = stiffScaleEngineId

       !! Connect the specified engines and functions to this spring
       call ConnectFunctions (lengthEngineId,stiffFuncId,forceFuncId, &
            &                 stiffScaleEnginePosId,stiffScaleEngineNegId, &
            &                 springs(idIn),stat)
       if (stat < 0) then
          err = err - 1
          call ReportInputError ('SPRING_BASE',idIn,springs(idIn)%id)
          cycle
       end if

       !! Connect to failure criterion, if any
       if (springFailureId > 0) then
          springs(idIn)%failure => GetPtrToId(springFailures,springFailureId)
          if (.not. associated(springs(idIn)%failure)) then
             err = err - 1
             call ReportInputError ('SPRING_BASE',idIn,springs(idIn)%id, &
                  msg='Non-existing spring failure id='//StrId(springFailureId))
             cycle
          end if
       end if

       !! Connect to yield criterion, if any
       if (springYieldId > 0) then
          springs(idIn)%yield => GetPtrToId(springYields,springYieldId)
          if (.not. associated(springs(idIn)%yield)) then
             err = err - 1
             call ReportInputError ('SPRING_BASE',idIn,springs(idIn)%id, &
                  msg='Non-existing spring yield id='//StrId(springYieldId))
             cycle
          end if
       end if

       if (unLoadType > 0 .and. associated(springs(idIn)%forceFunction)) then
          springs(idIn)%unLoadType = unLoadType
          if (abs(springs(idIn)%s0) <= 1.0e-16_dp) then
             !! Find the (initial) tangent/secant stiffness by differentiating
             !! the force-deflection function at deflection=0.0
             springs(idIn)%s0 = FunctionDerivative(springs(idIn)%forceFunction,&
                  &                                0.0_dp,1,err)
          end if
          stat = getErrorFile()
          write(stat,600) trim(getId(springs(idIn)%id)), springs(idIn)%s0
600       format(5X,'Initial tangent stiffness for Spring',A,' :',1PE12.5)
       end if

       !! Determine which secondary spring variables will be saved. Settings in
       !! the solver input file are overruled by command-line arguments, if any.
       do i = 1, size(saveVar)
          springs(idIn)%saveVar(i) = saveVar(i) > 0 .or. allSec .or. allSpring
       end do
       if (allStiff)  springs(idIn)%saveVar(1) = .true. ! Spring stiffness
       if (allLength) springs(idIn)%saveVar(2) = .true. ! Spring length
       if (allDef)    springs(idIn)%saveVar(3) = .true. ! Spring deflection
       if (allForce)  springs(idIn)%saveVar(4) = .true. ! Spring force
       if (allEnergy) springs(idIn)%saveVar(5) = .true. ! Spring energy

    end do

    if (err < 0) call reportError (debugFileOnly_p,'InitiateBaseSprings')

  contains

    !> @brief Connects functions and engines to the given spring base object.
    subroutine ConnectFunctions (lengthId,stiffId,forceId, &
         &                       posScaleId,negScaleId,spr,err)

      use FunctionTypeModule, only : GetPtrToId

      integer             , intent(in)  :: lengthId, stiffId, forceId
      integer             , intent(in)  :: posScaleId, negScaleId
      type(SpringBaseType), intent(out) :: spr
      integer             , intent(out) :: err

      err = 0

      if (lengthId > 0) then
         !! Connect the length engine
         spr%length0Engine => GetPtrToId(engines,lengthId)
         if (.not. associated(spr%length0Engine)) err = err - 1
      end if

      if (stiffId > 0) then
         !! Connect the stiffness-deflection function
         spr%stiffnessFunction => GetPtrToId(functions,stiffId)
         if (.not. associated(spr%stiffnessFunction)) err = err - 1
      else if (forceId > 0) then
         !! Connect the force-deflection function
         spr%forceFunction => GetPtrToId(functions,forceId)
         if (.not. associated(spr%forceFunction)) err = err - 1
      end if

      if (posScaleId > 0) then
         !! Connect the stiffness scale engine for tension
         spr%stiffScaleEnginePos => GetPtrToId(engines,posScaleId)
         if (.not. associated(spr%stiffScaleEnginePos)) err = err - 1
      end if

      if (negScaleId > 0) then
         !! Connect the stiffness scale engine for compression
         spr%stiffScaleEngineNeg => GetPtrToId(engines,negScaleId)
         if (.not. associated(spr%stiffScaleEngineNeg)) err = err - 1
      end if

    end subroutine ConnectFunctions

  end subroutine InitiateBaseSprings


  !!============================================================================
  !> @brief Deactivates all base springs that are not referred, for consistency.
  !>
  !> @param springs Array of all spring base objects in the model
  !> @param[out] stat Number of spring base objects that are deactivated
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 5 Jul 2004

  subroutine CheckSpringConnections (springs,stat)

    use IdTypeModule     , only : getId
    use SpringTypeModule , only : SpringBaseType
    use reportErrorModule, only : reportError, warning_p, warningFileOnly_p

    type(SpringBaseType), intent(inout) :: springs(:)
    integer             , intent(out)   :: stat

    !! Local variables
    integer :: i

    !! --- Logic section ---

    stat = 0
    do i = 1, size(springs)
       if (springs(i)%dof < 0) then
          stat = stat + 1
          springs(i)%isActive = .false.
          call reportError (warningFileOnly_p, &
               'Base Spring'//trim(getId(springs(i)%id))// &
               ' is not referenced.')
       end if
    end do

    if (stat > 0) then
       call reportError (warning_p,'The mechanism may be inconsistent.', &
            &            'Some spring objects are not referenced', &
            &            'and are ignored by the Dynamics Solver', &
            &            addString='CheckSpringConnections')
    end if

  end subroutine CheckSpringConnections

end module InitiateSpringTypeModule
