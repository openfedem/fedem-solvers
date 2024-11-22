!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file initiateSupElTypeModule.f90
!> @brief Initialization of superelements from the solver input file.

!!==============================================================================
!> @brief Initialization of superelements from the solver input file.
module initiateSupElTypeModule

  implicit none

  private

  public :: ReadSupEls, InitiateSupEls2, WriteSupEls2Ftn


contains

  !!============================================================================
  !> @brief Initializes superelements with data from the input file.
  !>
  !> @param[in] infp File unit number for the solver input file
  !> @param[in] env Environmental data
  !> @param[in] triads Array of all triads in the model
  !> @param[out] sups Array of all superelements in the model
  !> @param[out] err Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 1 Jun 2000

  subroutine ReadSupEls (infp,env,triads,sups,err)

    use EnvironmentTypeModule  , only : EnvironmentType
    use TriadTypeModule        , only : TriadType
    use TriadTypeModule        , only : GetPtrToId, AllocateTriadForces
    use SupElTypeModule        , only : SupElType, getSupElId
    use SupElTypeModule        , only : nullifySupEl, deallocateSupEls
    use SupElTypeModule        , only : InitiateHyDyn, IsBeam
    use SupElTypeModule        , only : TransformOffset, updateSupElCorot
    use IdTypeModule           , only : initId, deallocateId, getId, StrId
    use IdTypeModule           , only : ReportInputError
    use SupElNamelistModule ! Defines all variables in the SUP_EL namelist
    use GenericPartModule      , only : ReadGenericPart
    use FiniteElementModule    , only : massGrowth, massIntFluid
    use FiniteElementModule    , only : BeamPropertyType, GetPtrToId
    use FiniteElementModule    , only : ReadElementProps, CalcElementMatrices
    use FiniteElementModule    , only : FindBeams, FindBeamJoints
    use FNVwaveForceModule     , only : findColumns
    use inputUtilities         , only : iuGetNumberOfEntries
    use inputUtilities         , only : iuSetPosAtNextEntry
    use progressModule         , only : lterm
    use fileUtilitiesModule    , only : getDBGfile
    use dbgUnitsModule         , only : dbgShadowPos
    use reportErrorModule      , only : debugFileOnly_p, note_p, warning_p
    use reportErrorModule      , only : reportError, getErrorFile
    use reportErrorModule      , only : allocationError, internalError
    use FFaBodyHandlerInterface, only : FFa_initbody
    use FFaFilePathInterface   , only : FFa_checkPath
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_getbool
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_getint
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_getdouble
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_intValue

    integer              , intent(in)  :: infp
    type(EnvironmentType), intent(in)  :: env
    type(TriadType)      , intent(in)  :: triads(:)
    type(SupElType)      , pointer     :: sups(:), tmp(:)
    integer              , intent(out) :: err

    !! Local variables
    integer :: i, idIn, j, nDOF, nDOFs, nSupEl, nFnod, stat, lerr
    integer :: forceFNV, stiffDampFiltering, defaultShadowPosAlg, lpuDbg
    logical :: allPri, allSec, allRest, allSupel
    logical :: allCG, allHD, allGenDOF, allEnergy
    logical :: centripForceCorrIsOn(3), stressStiffIsOn(3)
    logical :: scaledStructDamp, scaleD, haveDamping(0:2), hasScaling
    real(dp) :: oMassScale, oMassDmp, oStiffDmp, glbAlpha1, glbAlpha2, stopGD
    real(dp), parameter :: eps_p = 1.0e-16_dp
    type(BeamPropertyType), pointer :: beamProp(:), bProp

    !! Statement function
    hasScaling(oMassScale) = abs(oMassScale-1.0_dp)*abs(oMassScale) > eps_p

    !! --- Logic section ---

    nullify(beamProp)
    call ReadElementProps (infp,beamProp,err)
    if (err < 0) goto 900

    nSupEl = iuGetNumberOfEntries(infp,'&SUP_EL',err)
    if (err /= 0 .or. nSupEl < 1) goto 900
    write(lterm,*) 'Number of &SUP_EL =',nSupEl

    allocate(tmp(nSupEl),STAT=stat)
    if (stat == 0) then
       call deallocateSupEls (sups)
       sups => tmp
    else
       err = AllocationError('ReadSupEls 1')
       goto 910
    end if

    call ffa_cmdlinearg_getbool ('allPrimaryVars',allPri)
    call ffa_cmdlinearg_getbool ('allSecondaryVars',allSec)
    call ffa_cmdlinearg_getbool ('allRestartVars',allRest)
    call ffa_cmdlinearg_getbool ('allSupelVars',allSupel)
    call ffa_cmdlinearg_getbool ('allCGVars',allCG)
    call ffa_cmdlinearg_getbool ('allHDVars',allHD)
    call ffa_cmdlinearg_getbool ('allGenDOFVars',allGenDOF)
    call ffa_cmdlinearg_getbool ('allEnergyVars',allEnergy)
    call ffa_cmdlinearg_getbool ('centripForceCorr',centripForceCorrIsOn(1))
    call ffa_cmdlinearg_getbool ('centripCorrToFD',centripForceCorrIsOn(2))
    call ffa_cmdlinearg_getbool ('centripCorrToFI',centripForceCorrIsOn(3))
    call ffa_cmdlinearg_getbool ('stressStiffDyn',stressStiffIsOn(1))
    call ffa_cmdlinearg_getbool ('stressStiffEqu',stressStiffIsOn(2))
    call ffa_cmdlinearg_getbool ('stressStiffEig',stressStiffIsOn(3))

    call ffa_cmdlinearg_getdouble ('alpha1',glbAlpha1)
    call ffa_cmdlinearg_getdouble ('alpha2',glbAlpha2)
    call ffa_cmdlinearg_getdouble ('stopGlbDmp',stopGD)
    if ((glbAlpha1 > 0.0_dp .or. glbAlpha2 > 0.0_dp) .and. stopGD > 0.0_dp) then
       scaledStructDamp = .true. ! Rayleigh damping based on scaled matrices
    else
       call ffa_cmdlinearg_getbool ('scaledStructDamp',scaledStructDamp)
    end if

    call ffa_cmdlinearg_getdouble ('overrideMassScale',oMassScale)
    call ffa_cmdlinearg_getdouble ('overrideMassPropDamp',oMassDmp)
    call ffa_cmdlinearg_getdouble ('overrideStiffPropDamp',oStiffDmp)

    call ffa_cmdlinearg_getint ('FNV',forceFNV)
    call ffa_cmdlinearg_getint ('stiffDampFiltering',stiffDampFiltering)
    call ffa_cmdlinearg_getint ('shadowPosAlg',defaultShadowPosAlg)
    if (defaultShadowPosAlg < 0) dbgShadowPos = getDBGfile(14,'shadowPos.dbg')

    massGrowth = 0.0_dp
    massIntFluid = 0.0_dp

    nFnod = 0
    do idIn = 1, nSupEl
       lerr = err

       call nullifySupEl (sups(idIn))
       if (.not. iuSetPosAtNextEntry(infp,'&SUP_EL',stat)) then
          err = err - 1
          call ReportInputError('SUP_EL',idIn)
          if (stat < 0) then
             exit
          else
             cycle
          end if
       end if

       call read_SUP_EL (infp,stat)
       if (stat /= 0) then
          err = err - 1
          call ReportInputError('SUP_EL',idIn)
          cycle
       end if

       call initId (sups(idIn)%id,id,extId,extDescr,stat)
       sups(idIn)%nExtNods = numTriads
       sups(idIn)%nLoadCase = numLoadCase
       if (numTriads == 2 .and. elPropId > 0) then
          write(lterm,*) '   Beam', trim(getId(sups(idIn)%id)), &
               &         '   elPropId =', elPropId
       else
          write(lterm,*) '   Part', trim(getId(sups(idIn)%id)), &
               &         '   nTriads =', numTriads, &
               &         '   nGenDOFs =', numGenDOFs
       end if
       if (shadowPosAlg /= 0) then
          sups(idIn)%shadowPosAlg = shadowPosAlg
       else
          sups(idIn)%shadowPosAlg = abs(defaultShadowPosAlg)
       end if

       if (sups(idIn)%shadowPosAlg > 0) then
          !! Initiate the superelement-wise flags with the global flag values
          !! Override global values if corresponding superelement flags are set
          do i = 1, 3
             if (stressStiffFlag(i) >= 0) then
                sups(idIn)%stressStiffFlag(i) = stressStiffFlag(i)
             else if (stressStiffIsOn(i)) then
                sups(idIn)%stressStiffFlag(i) = 1
             end if
          end do
          if (massCorrFlag >= 0) then
             sups(idIn)%massCorrFlag = massCorrFlag
          else if (centripForceCorrIsOn(1)) then
             sups(idIn)%massCorrFlag = 1
          end if
          if (sups(idIn)%massCorrFlag == 1) then
             if (centripForceCorrIsOn(2)) then
                sups(idIn)%massCorrFlag = 2
             else if (centripForceCorrIsOn(3)) then
                sups(idIn)%massCorrFlag = 3
             end if
          end if
       end if

       !! Override superelement scaling factors by global values, if set
       if (oMassScale >= 0.0_dp) massScale = oMassScale
       if (oMassDmp   >= 0.0_dp) alpha1    = oMassDmp
       if (oStiffDmp  >= 0.0_dp) alpha2    = oStiffDmp

       !! Rayleigh structural damping coefficients
       haveDamping(1) = alpha1 > 0.0_dp
       haveDamping(2) = alpha2 > 0.0_dp
       if (numGenDOFs > 0) then
          if (any(alpha3(1:numGenDOFs) > 0.0_dp)) then
             haveDamping(1) = .true.
          else
             alpha3(1:numGenDOFs) = alpha1
          end if
          if (any(alpha4(1:numGenDOFs) > 0.0_dp)) then
             haveDamping(2) = .true.
          else
             alpha4(1:numGenDOFs) = alpha2
          end if
       end if
       if (stopGD < 0.0_dp) then
          if (glbAlpha1 > 0.0_dp) haveDamping(1) = .true.
          if (glbAlpha2 > 0.0_dp) haveDamping(2) = .true.
       end if
       haveDamping(0) = haveDamping(1) .or. haveDamping(2)

       if (bodyFile /= '' .and. env%rhow > 0.0_dp) then
          !! Hydrodynamics loads are requested
          allocate(sups(idIn)%hydyn, STAT=stat)
          if (stat /= 0) then
             err = AllocationError('ReadSupEls 1a')
             goto 910
          end if
          call initiateHyDyn (sups(idIn)%hydyn)
       end if

       !! --- Initialization of triad connectivities and positions -------------

       if (numTriads > maxTriads_p) then
          err = internalError('ReadSupEls: Too many superelement triads'// &
               &              StrId(numTriads))
          goto 910
       end if

       allocate(sups(idIn)%triads(numTriads), &
            &   sups(idIn)%TrUndeformed(3,4,numTriads), STAT=stat)
       if (stat /= 0) then
          err = AllocationError('ReadSupEls 2')
          goto 910
       end if

       do i = 1, numTriads
          nullify(sups(idIn)%triads(i)%p)
          sups(idIn)%triads(i)%firstDOF = 0
       end do
       sups(idIn)%TrUndeformed = 0.0_dp

       if (maxval(nodeIds(1:numTriads)) > 0) then
          !! External FE node numbers are provided
          allocate(sups(idIn)%nodeId(numTriads), STAT=stat)
          if (stat /= 0) then
             err = AllocationError('ReadSupEls 2a')
             goto 910
          end if
          sups(idIn)%nodeId = nodeIds(1:numTriads)
       end if

       do i = 1, numTriads

          if (.not. iuSetPosAtNextEntry(infp,'&TRIAD_UNDPOS')) then
             err = err - 1
             lerr = lerr - 1
             call ReportTriadError(sups(idIn)%id,i)
             cycle
          end if

          call read_TRIAD_UNDPOS (infp,stat)
          if (stat /= 0) then
             err = err - 1
             lerr = lerr - 1
             call ReportTriadError(sups(idIn)%id,i)
             cycle
          else if (supElId /= id) then
             err = err - 1
             lerr = lerr - 1
             call ReportTriadError(sups(idIn)%id,triadId,supElId)
             cycle
          end if

          !! The UNDPOS records are not necessarily in the same order
          !! as the superelement triads. Therefore we must do a search.
          do j = 1, numTriads
             if (triadId == triadIds(j)) goto 100
          end do
          err = err - 1
          lerr = lerr - 1
          call ReportTriadError(sups(idIn)%id,triadId,supElId, &
               'This triad is not connected to this superelement')
          cycle

100       continue
          sups(idIn)%triads(j)%p => GetPtrToId(triads,triadId)
          sups(idIn)%TrUndeformed(:,:,j) = transpose(undPosInSupElSystem)

          !! Check for eccentricity vectors
          if (any(abs(eccVec) > eps_p) .or. any(abs(eccMass) > eps_p)) then
             if (.not. associated(sups(idIn)%eccVec)) then
                allocate(sups(idIn)%eccVec(3,2*numTriads), STAT=stat)
                if (stat /= 0) then
                   err = AllocationError('ReadSupEls 2b')
                   goto 910
                end if
                sups(idIn)%eccVec = 0.0_dp
             end if
             sups(idIn)%eccVec(:,j) = eccVec
             sups(idIn)%eccVec(:,numTriads+j) = eccMass
          end if

       end do

       !! Count the number of DOFs in the superelement
       nDof = 0
       do i = 1, numTriads
          if (associated(sups(idIn)%triads(i)%p)) then
             sups(idIn)%triads(i)%firstDOF = nDof + 1
             nDof = nDof + sups(idIn)%triads(i)%p%nDOFs
          else
             err = err - 1
          end if
       end do
       sups(idIn)%nTotDOFs = nDof + numGenDOFs
       if (err < lerr) then
          call ReportInputError('SUP_EL',idIn,sups(idIn)%id, &
               'Non-existing triad(s) referred')
          cycle
       end if

       if (ffa_cmdlinearg_intValue('debug') < 2) then
          lpuDbg = 0
       else if (numTriads < 10) then
          lpuDbg = getErrorFile()
          write(lpuDbg,"(/)")
       else
          lpuDbg = getDBGfile(20+idIn, &
               &              'supEl'//trim(adjustl(StrId(idIn)))//'.dbg')
       end if

       !! Set up the local superelement coordinate system
       sups(idIn)%supTr = transpose(supPos)
       if (sups(idIn)%shadowPosAlg == 1) then
          do i = 1, numTriads
             if (triadIds(i) == refTriad1Id) sups(idIn)%refTriad(1) = i
             if (triadIds(i) == refTriad2Id) sups(idIn)%refTriad(2) = i
             if (triadIds(i) == refTriad3Id) sups(idIn)%refTriad(3) = i
          end do
          err = err - count(sups(idIn)%refTriad == 0)

          call TransformOffset (sups(idIn),offset1,1)
          call TransformOffset (sups(idIn),offset2,2)
          call TransformOffset (sups(idIn),offset3,3)
          stat = 1 ! Initiate the co-rotated superelement coordinate system
          call updateSupElCorot (sups(idIn),dbgShadowPos,stat)
          if (stat < 0) err = err - 1
       end if

       !! --- Allocation of superelement matrices and vectors ------------------

       nDof = sups(idIn)%nTotDOFs
       allocate(sups(idIn)%Q        (nDof), &
            &   sups(idIn)%FS       (nDof), &
            &   sups(idIn)%FD       (nDof), &
            &   sups(idIn)%FI       (nDof), &
            &   sups(idIn)%finit    (nDof), &
            &   sups(idIn)%finitPrev(nDof), &
            &   sups(idIn)%uld      (nDof), &
            &   sups(idIn)%uldd     (nDof), &
            &   STAT=stat)
       if (stat /= 0) then
          err = AllocationError('ReadSupEls 3')
          goto 910
       end if

       sups(idIn)%Q         = 0.0_dp
       sups(idIn)%FS        = 0.0_dp
       sups(idIn)%FD        = 0.0_dp
       sups(idIn)%FI        = 0.0_dp
       sups(idIn)%finit     = 0.0_dp
       sups(idIn)%finitPrev = 0.0_dp
       sups(idIn)%uld       = 0.0_dp
       sups(idIn)%uldd      = 0.0_dp

       do i = 1, numTriads
          nDOFs = sups(idIn)%triads(i)%p%nDOFs
          if (nDOFs > 0) then
             lerr = err
             j = sups(idIn)%triads(i)%firstDOF
             sups(idIn)%triads(i)%p%ur_def => sups(idIn)%finit(j:j+nDOFs-1)
             call allocateTriadForces (sups(idIn)%triads(i)%p,err)
             if (err < lerr) goto 900
          end if
       end do

       nDOFs = nDof - numGenDOFs
       if (alpha2 > 0.0_dp .and. stiffDampFiltering > 0) then
          !! Additional array(s) for stiffness-proportional damping calculation
          allocate(sups(idIn)%vld(nDOFs), STAT=stat)
          if (stat /= 0) then
             err = AllocationError('ReadSupEls 3a')
             goto 910
          end if
          if (stiffDampFiltering > 1) then
             allocate(sups(idIn)%vldPrev (nDOFs), &
                  &   sups(idIn)%vldd    (nDOFs), &
                  &   sups(idIn)%vlddPrev(nDOFs), &
                  &   STAT=stat)
             if (stat /= 0) then
                err = AllocationError('ReadSupEls 3b')
                goto 910
             end if
             sups(idIn)%vldPrev  = 0.0_dp
             sups(idIn)%vldd     = 0.0_dp
             sups(idIn)%vlddPrev = 0.0_dp
          end if
          sups(idIn)%vld = 0.0_dp
       end if

       if (numGenDOFs > 0) then
          !! Additional arrays for component modes
          allocate(sups(idIn)%genDOFs, STAT=stat)
          if (stat /= 0) then
             err = AllocationError('ReadSupEls 4')
             goto 910
          end if
          sups(idIn)%genDOFs%nDOFs = numGenDOFs
          sups(idIn)%genDOFs%firstDOF = nDOFs + 1
          allocate(sups(idIn)%genDOFs%BC(numGenDOFs), &
               &   sups(idIn)%genDOFs%alpha1(numGenDOFs), &
               &   sups(idIn)%genDOFs%alpha2(numGenDOFs), &
               &   sups(idIn)%genDOFs%ur(numGenDOFs), &
               &   sups(idIn)%genDOFs%urPrev(numGenDOFs), &
               &   sups(idIn)%genDOFs%energy(3,numGenDOFs), STAT=stat)
          nullify(sups(idIn)%genDOFs%urd)
          nullify(sups(idIn)%genDOFs%urdd)
          nullify(sups(idIn)%genDOFs%sysDOF)
          if (stat /= 0) then
             err = AllocationError('ReadSupEls 5')
             goto 910
          end if
          sups(idIn)%genDOFs%BC     = BC(1:numGenDOFs)
          sups(idIn)%genDOFs%ur     = 0.0_dp
          sups(idIn)%genDOFs%urPrev = 0.0_dp
          sups(idIn)%genDOFs%energy = 0.0_dp
       end if

       allocate(sups(idIn)%Nmat (nDof,nDof), &
            &   sups(idIn)%Mmat (nDof,nDof), &
            &   sups(idIn)%KmMat(nDof,nDof), &
            &   sups(idIn)%fg   (nDof,3), &
            &   STAT=stat)
       if (stat /= 0) then
          err = AllocationError('ReadSupEls 6')
          goto 910
       end if

       if (any(sups(idIn)%stressStiffFlag > 0)) then
          !! Tangent stiffness is allocated only if stress stiffening is active
          allocate(sups(idIn)%KtMat(nDof,nDof), STAT=stat)
          if (stat /= 0) then
             err = AllocationError('ReadSupEls 6a')
             goto 910
          end if
       end if

       if (numLoadCase > 0) then
          !! External load vectors
          allocate(sups(idIn)%S(nDof,numLoadCase), STAT=stat)
          if (stat /= 0) then
             err = AllocationError('ReadSupEls 6c')
             goto 910
          end if
       end if

       !! --- Input or generation of the constant superelement matrices --------

       nullify(bProp)
       do i = 1, size(inputFiles)
          if (inputFiles(i) /= '') call FFa_checkPath (inputFiles(i))
       end do
       if (numStates > 0) then
          call ReadNonlinSupElMatrices (inputFiles,sups(idIn),numStates,stat)
       else if (inputFiles(1) /= '') then
          if (inputFiles(7) /= '' .and. recoveryFlag > 0) then
             allocate(sups(idIn)%rcy, STAT=stat)
             if (stat /= 0) then
                err = AllocationError('ReadSupEls 7')
                goto 910
             end if
             sups(idIn)%rcy%fileName = inputFiles(4:7)
             sups(idIn)%rcy%recovery = recoveryFlag
             if (mod(recoveryFlag,2) == 0) then
                sups(idIn)%rcy%elmGroup = '' ! Only strain gage recovery
             else if (ffa_cmdlinearg_intValue('recovery') > 10) then
                sups(idIn)%rcy%elmGroup = '' ! Ignore the element groups
             else
                sups(idIn)%rcy%elmGroup = elmGroups ! Groups for stress recovery
             end if
          end if
          call ReadSupElMatrices (inputFiles,sups(idIn),stat)
       else if (elPropId > 0) then
          bProp => GetPtrToId(beamProp,elPropId)
          if (associated(bProp)) then
             sups(idIn)%rigidFlag = -elPropId ! Beam element
             if (associated(sups(idIn)%hydyn)) then
                call initiateHyDyn (sups(idIn)%hydyn,bProp%Morison)
                if (bodyFile == 'FILLED') then
                   sups(idIn)%hydyn%bodyIndex = -2 ! Indicate internal fluid
                end if
             end if
             call CalcElementMatrices (sups(idIn),bProp,env,stat)
          else
             call ReportInputError('SUP_EL',idIn,sups(idIn)%id, &
                  'Non-existing beam element property referred')
             err = err - 1
          end if
       else
          call ReadGenericPart (infp,sups(idIn),stat)
       end if
       if (stat /= 0) then
          err = err - 1
          cycle
       end if

       if (.not. associated(sups(idIn)%Mmat)) then
          haveDamping = .false. ! No damping in pure (quasi-)static analysis
          sups(idIn)%massCorrFlag = 0 ! Deactivate mass-matrix correction
       else if (haveDamping(0) .or. glbAlpha1+glbAlpha2 > 0.0_dp) then
          !! Damping matrix is allocated only if structural damping is active
          allocate(sups(idIn)%Cmat(nDof,nDof), STAT=stat)
          if (stat /= 0) then
             err = AllocationError('ReadSupEls 6b')
             goto 910
          end if
       end if

       scaleD = scaledStructDamp .or. stiffEngineId > 0 .or. .not.haveDamping(0)
       if (scaleD) then

          !! Apply stiffness and/or mass scaling factors
          !! before computing the structural damping matrix, if any
          call scaleMatrices (sups(idIn)%KmMat,sups(idIn)%Mmat,sups(idIn)%fg)
          if (haveDamping(1) .and. hasScaling(massScale)) then
             call reportError (note_p,'Mass-proportional damping for '// &
                  &            trim(getSupElId(sups(idIn))) // &
                  &            ' is based on the scaled mass matrix')
          end if
          if (haveDamping(2) .and. hasScaling(stiffScale)) then
             call reportError (note_p,'Stiffness-proportional damping for '// &
                  &            trim(getSupElId(sups(idIn))) // &
                  &            ' is based on the scaled stiffness matrix')
          end if

       end if
       if (haveDamping(0)) then

          !! Compute the structural damping matrix
          call computeStructDamping (sups(idIn)%Cmat, &
               &                     sups(idIn)%Mmat,sups(idIn)%KmMat, &
               &                     glbAlpha1+alpha1,glbAlpha2+alpha2)

       end if
       if (.not. scaleD) then

          !! Apply stiffness and/or mass scaling factors. The structural damping
          !! above is now based on the unscaled mass- and stiffness matrices.
          call scaleMatrices (sups(idIn)%KmMat,sups(idIn)%Mmat,sups(idIn)%fg)

          !! Modify the structural damping parameters such that they yield
          !! consistent damping forces when applied to the scaled mass- and
          !! stiffness matrices later (see UpdateSupEls and SupElEnergies)
          if (hasScaling(massScale)) then
             alpha1 = alpha1 / massScale
             call DSCAL (numGenDOFs,1.0_dp/massScale,alpha3(1),1)
          end if
          if (hasScaling(stiffScale)) then
             alpha2 = alpha2 / stiffScale
             call DSCAL (numGenDOFs,1.0_dp/stiffScale,alpha4(1),1)
          end if

       end if
       if (strDmpEngineId > 0 .and. (alpha1 > 0.0_dp .or. alpha2 > 0.0_dp)) then
          !! Time-dependent structural damping
          sups(idIn)%dmpSclIdx = strDmpEngineId
          sups(idIn)%mDmp0 = alpha1
          sups(idIn)%kDmp0 = alpha2
       else
          sups(idIn)%mDmpFactor = alpha1
          sups(idIn)%kDmpFactor = alpha2
       end if
       if (numGenDOFs > 0) then
          sups(idIn)%genDOFs%alpha1 = alpha3(1:numGenDOFs)
          sups(idIn)%genDOFs%alpha2 = alpha4(1:numGenDOFs)
       end if

       if (massEngineId > 0) then

          !! Time-dependent mass scaling
          if (haveDamping(1)) then
             call reportError (note_p,'The time-dependent mass scaling '// &
                  'assigned to '//getSupElId(sups(idIn)), &
                  'is not accounted for when calculating the '// &
                  'mass-proportional damping contributions.')
          end if
          sups(idIn)%massSclIdx = massEngineId

       end if
       if (stiffEngineId > 0) then

          !! Time-dependent stiffness scaling
          if (haveDamping(1)) then
             call reportError (warning_p,'Time-dependent stiffness scaling '// &
                  'can not be combined with mass-proportional damping.', &
                  'The time-dependent scaling is therefore switched off for '//&
                  getSupElId(sups(idIn)))
          else
             sups(idIn)%stifSclIdx = stiffEngineId
          end if

       end if

       if (associated(sups(idIn)%hydyn)) then

          !! Currently no geometry-based buoyancy for beams
          if (associated(bProp)) bodyFile = 'NONE'

          if (bodyFile /= 'NONE') then
             !! Read body data for buoyancy and drag/slam calculation
             call FFa_checkPath (bodyFile)
             call FFa_initbody (trim(bodyFile),sups(idIn)%hydyn%bodyIndex)
             if (sups(idIn)%hydyn%bodyIndex < 0) err = err - 1
          else if (sups(idIn)%rigidFlag >= 0) then
             err = err - 1
             call ReportInputError('SUP_EL',idIn,sups(idIn)%id, &
                  'Hydrodynamics for superelements without '// &
                  'a geometry description is not supported.')
          else if (numTriads /= 2) then
             err = err - 1
             call ReportInputError('SUP_EL',idIn,sups(idIn)%id, &
                  'Hydrodynamics not supported for'//trim(StrId(numTriads))// &
                  '-noded flexible elements')
          end if

          !! Drag and slam parameters
          nDOFs = 0
          do i = 1, 6
             if (dragParams(1,i) > 0.0_dp) nDOFs = i
          end do
          if (nDOFs > 0 .and. sups(idIn)%rigidFlag >= 0) then
             allocate(sups(idIn)%hydyn%dragPar(3,nDOFs), STAT=stat)
             if (stat /= 0) then
                err = AllocationError('ReadSupEls 8')
                goto 910
             end if
             sups(idIn)%hydyn%dragPar = dragParams(:,1:nDOFs)
          end if
          sups(idIn)%hydyn%slamPar = slamParams

          if (any(sups(idIn)%stressStiffFlag > 0) .and. bodyFile /= 'NONE') then
             !! Load correction stiffness for buoyancy forces
             allocate(sups(idIn)%KlMat(nDof,nDof), STAT=stat)
             if (stat /= 0) then
                err = AllocationError('ReadSupEls 8a')
                goto 910
             end if
             sups(idIn)%KlMat = 0.0_dp
          end if

          if (associated(bProp)) then
             !! Added mass and damping from Morison equation
             allocate(sups(idIn)%MaMat(nDof,nDof), &
                  &   sups(idIn)%CdMat(nDof,nDof), STAT=stat)
             if (stat /= 0) then
                err = AllocationError('ReadSupEls 8b')
                goto 910
             end if
             sups(idIn)%MaMat = 0.0_dp
             sups(idIn)%CdMat = 0.0_dp
             !! Mark all triads that need hydrodynamics calculation
             do j = 1, sups(idIn)%nExtNods
                if (sups(idIn)%triads(j)%p%inFluid == 0) then
                   nFnod = nFnod + 1
                   sups(idIn)%triads(j)%p%inFluid = nFnod
                end if
             end do
          end if

       end if

       !! Determine if primary variables (position matrix) will be saved
       sups(idIn)%savePos = allPri .or. allRest .or. allSupel .or. savePos > 0

       !! Determine which secondary supel variables will be saved. Settings in
       !! the solver input file are overruled by command-line arguments, if any.
       do i = 1, size(saveVar)
          sups(idIn)%saveVar(i) = saveVar(i) > 0 .or. allSec .or. allSupel
       end do
       if (allRest)   sups(idIn)%saveVar(1:2) = .true. ! Restart variables
       if (allCG)     sups(idIn)%saveVar(1) = .true. ! Centre of gravity
       if (IsBeam(sups(idIn))) then
          sups(idIn)%saveVar(2) = saveVar(2) > 0     ! End sectional forces
       else if (allGenDOF) then
          sups(idIn)%saveVar(2) = .true.             ! Generalized DOFs (d,v,a)
       end if
       if (allEnergy) sups(idIn)%saveVar(3) = .true. ! Energy quantities
       if (allHD)     sups(idIn)%saveVar(4) = .true. ! Hydrodynamics quantities

    end do

    !! Check if we have any beams and beam joints
    if (err == 0) call findBeams (sups,err)
    if (err == 0 .and. forceFNV > 0) call findColumns(env%Tsea,beamProp,err)
    if (err == 0) call findBeamJoints (triads,sups,err)

900 if (err /= 0) call reportError (debugFileOnly_p,'ReadSupEls')
910 if (associated(beamProp)) then
       do i = 1, size(beamProp)
          call deallocateId (beamProp(i)%id)
       end do
       deallocate(beamProp)
    end if

  contains

    !> @brief Subroutine for debug print of a superelement matrix.
    subroutine writeSupMatrix (label,sup,A)
      use manipMatrixModule, only : writeObject

      character(len=*), intent(in) :: label
      type(SupElType) , intent(in) :: sup
      real(dp)        , intent(in) :: A(:,:)

      call writeObject (A,lpuDbg,label//' for '//getSupElId(sup),eps=1.0e-4_dp)
      write(lpuDbg,"()")

    end subroutine writeSupMatrix

    !> @brief Computes structural damping matrix assuming Rayleigh damping.
    subroutine computeStructDamping (C,M,K,alpha1,alpha2)

      real(dp), intent(out) :: C(:,:)
      real(dp), intent(in)  :: M(:,:), K(:,:), alpha1, alpha2

      call DCOPY (size(C),0.0_dp,0,C(1,1),1)

      if (abs(alpha1) > eps_p) then
         call DAXPY (size(C),alpha1,M(1,1),1,C(1,1),1)
      end if

      if (abs(alpha2) > eps_p) then
         call DAXPY (size(C),alpha2,K(1,1),1,C(1,1),1)
      end if

      !! Account for possibly individual damping factors on the generalized DOFs
      !! assuming that the sub-matrices associated with these DOFs are diagonal

      nDof = size(C,1) - numGenDOFs
      do i = 1, numGenDOFs
         j = nDof + i
         if (abs(alpha3(i)-alpha1) > eps_p) then
            C(j,j)      = C(j,j)      + (alpha3(i)-alpha1)*M(j,j)
            C(1:nDof,j) = C(1:nDof,j) + (alpha3(i)-alpha1)*M(1:nDof,j)
            C(j,1:nDof) = C(j,1:nDof) + (alpha3(i)-alpha1)*M(j,1:nDof)
         end if
         if (abs(alpha4(i)-alpha2) > eps_p) then
            C(j,j) = C(j,j) + (alpha4(i)-alpha2)*K(j,j)
            !! The coupling matrix between the triad DOFs and generalized DOFs
            !! is assumed always zero (results from the CMS reduction)
         end if
      end do

      if (lpuDbg < 1) return

      !! Debug output of structural damping matrix
      call writeSupMatrix ('Damping matrix',sups(idIn),C)

    end subroutine computeStructDamping

    !> @brief Applies scaling factors on the superelement matrices.
    subroutine scaleMatrices (Km,M,fg)

      real(dp), pointer :: Km(:,:), M(:,:), fg(:,:)

      if (abs(stiffScale-1.0_dp) > eps_p) then
         if (associated(Km)) then
            call DSCAL (size(Km),stiffScale,Km(1,1),1)
         end if
      end if
      if (abs(massScale-1.0_dp) > eps_p) then
         if (associated(M)) then
            call DSCAL (size(M),massScale,M(1,1),1)
         end if
         if (associated(fg)) then
            call DSCAL (size(fg),massScale,fg(1,1),1)
         end if
      end if

      if (lpuDbg < 1) return

      !! Debug output of reduced FE matrices

      if (associated(Km)) call writeSupMatrix ('Stiffness matrix',sups(idIn),Km)
      if (associated(M))  call writeSupMatrix ('Mass matrix',sups(idIn),M)
      if (associated(fg)) call writeSupMatrix ('Gravity forces',sups(idIn),fg)

    end subroutine scaleMatrices

  end subroutine ReadSupEls


  !!============================================================================
  !> @brief Initializes more superelement arrays after system initialization.
  !>
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[in] sys System level model data
  !> @param sups All superelements in the model
  !> @param[in] masses All point masses in the model
  !> @param[in] engines All general functions in the model
  !> @param[out] err Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 3 Aug 2000

  subroutine InitiateSupEls2 (sam,sys,sups,masses,engines,err)

    use SamModule                 , only : SamType, dp
    use SystemTypeModule          , only : SystemType
    use SupElTypeModule           , only : SupElType
    use MassTypeModule            , only : MassType
    use FunctionTypeModule        , only : EngineType, GetPtrToId
    use IdTypeModule              , only : ReportInputError
    use MassMatrixCorrectionModule, only : mmcRigidMassProperties
    use MassMatrixCorrectionModule, only : mmcInitMassMatrixCorrection
    use CorotUtilModule           , only : formShadowPosGrad
    use dbgUnitsModule            , only : dbgShadowPos
    use reportErrorModule         , only : reportError, debugFileOnly_p

    type(SamType)   , intent(in)    :: sam
    type(SystemType), intent(in)    :: sys
    type(SupElType) , intent(inout) :: sups(:)
    type(MassType)  , intent(in)    :: masses(:)
    type(EngineType), intent(in)    :: engines(:)
    integer         , intent(out)   :: err

    !! Local variables
    integer :: idIn, inod, nSupEls, firstDOF, lastDOF, lerr

    !! --- Logic section ---

    err = 0
    lerr = 0

    nSupEls = size(sups)
    do idIn = 1, nSupEls

       if (associated(sups(idIn)%genDOFs)) then
          !! Initialize pointers to the generalized DOFs of this superelement
          inod     = sups(idIn)%genDOFs%samNodNum
          firstDOF = sam%madof(inod)
          lastDOF  = sam%madof(inod+1)-1
          sups(idIn)%genDOFs%sysDOF => sam%madof(inod)
          sups(idIn)%genDOFs%urd    => sys%urd(firstDOF:lastDOF)
          sups(idIn)%genDOFs%urdd   => sys%urdd(firstDOF:lastDOF)
       end if

       if (associated(sups(idIn)%Mmat)) then
          !! Compute rigid body mass properties
          call mmcRigidMassProperties (sups(idIn),err)
          if (err /= 0) exit
       end if

       if (sups(idIn)%massCorrFlag > 0 .or. sups(idIn)%saveVar(3)) then
          !! Initiate matrices needed for mass matrix correction
          call mmcInitMassMatrixCorrection (sups(idIn),err)
          if (err /= 0) exit
       end if

       if ( any(sups(idIn)%stressStiffFlag > 0) .or. &
            sups(idIn)%shadowPosAlg > 1 ) then
          !! Establish gradients of shadow element position
          call formShadowPosGrad (sups(idIn),dbgShadowPos,err)
          if (err < 0) exit
          if (err > 0) lerr = lerr + 1
       end if

       sups(idIn)%supTrPrev = sups(idIn)%supTr
       sups(idIn)%supTrInit = sups(idIn)%supTr
       sups(idIn)%addedMass = hasAddedMass(sups(idIn))

       if (sups(idIn)%dmpSclIdx > 0) then
          if (.not. associated(GetPtrToId(engines,sups(idIn)%dmpSclIdx, &
               &                          index=sups(idIn)%dmpSclIdx))) then
             lerr = lerr - 1
             call ReportInputError('SUP_EL',idIn,sups(idIn)%id, &
                  &                'Non-existing engine referred')
          end if
       end if
       if (sups(idIn)%massSclIdx > 0) then
          if (.not. associated(GetPtrToId(engines,sups(idIn)%massSclIdx, &
               &                          index=sups(idIn)%massSclIdx))) then
             lerr = lerr - 1
             call ReportInputError('SUP_EL',idIn,sups(idIn)%id, &
                  &                'Non-existing engine referred')
          end if
       end if
       if (sups(idIn)%stifSclIdx > 0) then
          if (.not. associated(GetPtrToId(engines,sups(idIn)%stifSclIdx, &
               &                          index=sups(idIn)%stifSclIdx))) then
             lerr = lerr - 1
             call ReportInputError('SUP_EL',idIn,sups(idIn)%id, &
                  &                'Non-existing engine referred')
          end if
       end if

    end do
    if (lerr > 0 .and. err >= 0) err = -lerr

    if (err /= 0) call reportError (debugFileOnly_p,'InitiateSupEls2')

  contains

    !> @brief Checks if a superelement has additional masses on triads.
    function hasAddedMass (supel)

      type(SupElType), intent(in) :: supel
      logical :: hasAddedMass
      integer :: i, id, j, nMasses, nTriads

      hasAddedMass = .true.

      nMasses = size(masses)
      nTriads = supel%nExtNods
      do i = 1, nMasses
         if (abs(sum(masses(i)%m0))+abs(sum(masses(i)%m1)) > 1.0e-16_dp) then
            id = masses(i)%triad%id%baseId
            do j = 1, nTriads
               if (supel%triads(j)%p%id%baseId == id) return
            end do
         end if
      end do

      hasAddedMass = .false.

    end function hasAddedMass

  end subroutine InitiateSupEls2


  !!============================================================================
  !> @brief Reads into core the reduced superelement matrices from file.
  !>
  !> @param[in] fileNames List of files containing the superelement matrices
  !> @param sup The superelement to read matrices for
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 26 Sep 2000

  subroutine ReadSupElMatrices (fileNames,sup,ierr)

    use SupElTypeModule  , only : SupElType
    use reportErrorModule, only : reportError, error_p, note_p, debugFileOnly_p

    type(SupElType) , intent(inout) :: sup
    character(len=*), intent(in)    :: fileNames(4)
    integer         , intent(out)   :: ierr

    !! Local variables
    integer :: lerr(3), ndim

    !! --- Logic section ---

    lerr = 0
    ndim = sup%nTotDOFs

    !! Read reduced material stiffness matrix from file
    call readMatrix (fileNames(1),sup%id%baseId,'stiffness matrix', &
         &           ndim*ndim,sup%KmMat(1,1),lerr)
    if (lerr(3) > 0) then
       deallocate(sup%KmMat)
       nullify(sup%KmMat)
       call reportError (error_p,'No superelement stiffness provided.')
       lerr(1) = -100
    end if

    !! Read reduced mass matrix from file
    call readMatrix (fileNames(2),sup%id%baseId,'mass matrix', &
         &           ndim*ndim,sup%Mmat(1,1),lerr)
    if (lerr(3) > 0) then
       deallocate(sup%Mmat)
       nullify(sup%Mmat)
    end if

    !! Read gravity force vectors from file
    call readMatrix (fileNames(3),sup%id%baseId,'gravity force vectors', &
         &           ndim*3,sup%fg(1,1),lerr)
    if (lerr(3) > 0) then
       deallocate(sup%fg)
       nullify(sup%fg)
    end if

    if (associated(sup%S)) then
       !! Read load vectors from file
       call readMatrix (fileNames(4),sup%id%baseId,'load vectors', &
            &           ndim*sup%nLoadCase,sup%S(1,1),lerr)
       if (lerr(3) > 0) then
          deallocate(sup%S)
          nullify(sup%S)
       end if
    end if

    if (lerr(2) < 0) then
       call reportError (note_p,'The above errors can be due to '// &
            'an inconsistency in the number of external DOFs', &
            'in the reduced FE matrices compared to the number of triads', &
            '(and component modes) specified in the mechanism for that part.', &
            'Please verify that the reduced FE matrices are up to date '// &
            'with respect to the mechanism model.')
    else if (mod(lerr(1),100) < 0) then
       call reportError (note_p,'The above errors can be due to very long '// &
            'pathnames for the matrix files.','Try to simplify the '// &
            'pathnames if they are long or have many levels.')
    end if

    if (lerr(1) < 0) then
       call reportError (debugFileOnly_p,'ReadSupElMatrices')
    end if
    ierr = lerr(1)

  end subroutine ReadSupElMatrices


  !!============================================================================
  !> @brief Reads into core the reduced nonlinear superelement matrices.
  !>
  !> @param[in] fileNames List of files containing the superelement matrices
  !> @param sup The superelement to read matrices for
  !> @param[in] numData Number of states in the nonlinear reduction
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Leif Ivar Myklebust and Bjorn Haugen
  !>
  !> @date 11 Sep 2006

  subroutine ReadNonlinSupElMatrices (fileNames,sup,numData,ierr)

    use SupElTypeModule  , only : SupElType
    use reportErrorModule, only : allocationError, reportError
    use reportErrorModule, only : error_p, note_p, debugFileOnly_p

    type(SupElType) , intent(inout) :: sup
    character(len=*), intent(in)    :: fileNames(5)
    integer         , intent(in)    :: numData
    integer         , intent(out)   :: ierr

    !! Local variables
    integer :: lerr(3), ndim

    !! --- Logic section ---

    lerr = 0
    ndim = sup%nTotDOFs

    if (numData < 2) then
       ierr = -1
       call reportError (error_p,'Too few nonlinear states', &
            &            addString='ReadNonlinSupElMatrices')
       return
    end if

    allocate(sup%nonlin,STAT=ierr)
    if (ierr == 0) then
       allocate(sup%nonlin%force(ndim,numData), &
            &   sup%nonlin%disp (ndim,numData), &
            &   sup%nonlin%stiff(ndim,ndim,numData), &
            &   STAT=ierr)
    end if
    if (ierr /= 0) then
       ierr = AllocationError('ReadNonlinSupElMatrices')
       return
    end if

    !! Read reduced nonlinear material stiffness matrix from file
    call readMatrix (fileNames(1),sup%id%baseId,'stiffness matrix', &
         &           ndim*ndim*numData,sup%nonlin%stiff(1,1,1),lerr)
    !! Read reduced nonlinear force vectors from file
    call readMatrix (fileNames(5),sup%id%baseId,'force vectors', &
         &           ndim*numData,sup%nonlin%force(1,1),lerr)
    !! Read reduced nonlinear displacement vectors from file
    call readMatrix (fileNames(4),sup%id%baseId,'displacement vectors', &
         &           ndim*numData,sup%nonlin%disp(1,1),lerr)
    !! Read reduced mass matrix from file
    call readMatrix (fileNames(2),sup%id%baseId,'mass matrix', &
         &           ndim*ndim,sup%Mmat(1,1),lerr)
    !! Read gravity force vectors from file
    call readMatrix (fileNames(3),sup%id%baseId,'gravity force vectors', &
         &           ndim*3,sup%fg(1,1),lerr)

    if (lerr(2) < 0) then
       call reportError (note_p,'The above errors can be due to '// &
            'an inconsistency in the number of external DOFs', &
            'in the reduced FE matrices compared to the number of triads', &
            '(and component modes) specified in the mechanism for that part.', &
            'Please verify that the reduced FE matrices are up to date '// &
            'with respect to the mechanism model.')
    end if
    if (lerr(1) < 0) then
       call reportError (debugFileOnly_p,'ReadNonlinSupElMatrices')
    end if
    ierr = lerr(1)

  end subroutine ReadNonlinSupElMatrices


  !!============================================================================
  !> @brief Reads a matrix from a binary file into core.
  !>
  !> @param[in] fileName The file containing the matrix to read
  !> @param[in] supId Base ID of the superelement the matrix is associated with
  !> @param[in] dtype Type of matrix to read (stiffness, mass, etc.)
  !> @param[in] ndata Size of the matrix to read
  !> @param[out] data The matrix data
  !> @param ierr Error flags
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 26 Sep 2000

  subroutine ReadMatrix (fileName,supId,dtype,ndata,data,ierr)

    use kindModule       , only : dp
    use progressModule   , only : lterm
    use binaryDBInterface, only : readDoubleDB
    use reportErrorModule, only : reportError, error_p, debugFileOnly_p

    character(len=*), intent(in)    :: fileName, dtype
    integer         , intent(in)    :: supId, ndata
    real(dp)        , intent(out)   :: data
    integer         , intent(inout) :: ierr(:)

    !! Local variables
    integer :: lerr(2), stat2
    logical :: haveFile

    !! --- Logic section ---

    if (fileName == '') then
       haveFile = .false.
    else
       inquire(FILE=fileName,EXIST=haveFile)
    end if

    if (haveFile) then
       stat2 = 0
    else
       !! Check if the superelement matrix is hard-coded in a shared library
       call loadSupelMatrixFromCore (supId,dtype,data,ndata,stat2)
    end if

    lerr = ierr(1:2)
    if (stat2 == ndata) then
       write(lterm,100) dtype
100    format(7X,'Loading superelement ',A,' from in-core storage')
    else if (fileName == '') then
       ierr(3) = 1 ! No matrix provided, not necessarily an error
       return
    else if (stat2 >= 0) then
       !! Read the superelement matrix from file
       call readDoubleDB (fileName,dtype,ndata,data,ierr(1),ierr(2))
    end if

    ierr(3) = 0
    if (ierr(2) < lerr(2) .or. stat2 < 0) then
       call reportError (error_p,'Failed to read matrix file '//fileName, &
            &            addString='ReadMatrix')
       ierr(1) = ierr(1) - 1
    else if (ierr(1) < lerr(1)) then
       call reportError (debugFileOnly_p,'ReadMatrix')
    end if

  end subroutine ReadMatrix


  !!============================================================================
  !> @brief Writes a Fortran90 subroutine with hard-coded superelement matrices.
  !>
  !> @param sups All superelements in the model
  !> @param[out] ierr Error flag
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 2 Jan 2017

  subroutine WriteSupEls2Ftn (sups,ierr)

    use SupElTypeModule    , only : SupElType, getSupElId
    use fileUtilitiesModule, only : findUnitNumber
    use progressModule     , only : lterm
    use reportErrorModule  , only : reportError, error_p

    type(SupElType), intent(in)  :: sups(:)
    integer        , intent(out) :: ierr

    !! Local variables
    integer          :: iftn, iinc, idIn, n, nDof
    character(len=1) :: c

    !! --- Logic section ---

    iftn = findUnitNumber(20)
    open(iftn,FILE='loadSupelMatrix.f90',STATUS='UNKNOWN',IOSTAT=ierr)
    if (ierr == 0) then
       iinc = findUnitNumber(21)
       open(iinc,FILE='supels.inc',STATUS='UNKNOWN',IOSTAT=ierr)
    end if
    if (ierr /= 0) then
       call reportError (error_p,'Could not open Fortran output files', &
            &            addString='writeSupEls2Ftn')
       return
    end if

    n = 0
    do idIn = 1, size(sups)
       if (sups(idIn)%rigidFlag == 0) then
          write(lterm,*) 'Writing Fortran code for '// &
               &         trim(getSupElId(sups(idIn)))
          n = n + 1
          c = i2c(n)
          if (n == 1) then
             write(iftn,600)
          else
             write(iinc,"()")
          end if
          nDof = sups(idIn)%nTotDOFs
          write(iinc,601) c, sups(idIn)%id%baseId, c, nDof
          if (associated(sups(idIn)%KmMat)) then
             write(iinc,602) 's',c,c,'ndim'//c
             call writeMatrix (size(sups(idIn)%KmMat),sups(idIn)%KmMat)
          end if
          if (associated(sups(idIn)%Mmat)) then
             write(iinc,602) 'm',c,c,'ndim'//c
             call writeMatrix (size(sups(idIn)%Mmat),sups(idIn)%Mmat)
          end if
          if (associated(sups(idIn)%fg)) then
             write(iinc,602) 'g',c,c,'3'
             call writeMatrix (size(sups(idIn)%fg),sups(idIn)%fg)
          end if
          write(iftn,610) c
          if (associated(sups(idIn)%KmMat)) then
             write(iftn,611) 's', 'smat'//c, 'smat'//c
          end if
          if (associated(sups(idIn)%Mmat)) then
             write(iftn,611) 'm', 'mmat'//c, 'mmat'//c
          end if
          if (associated(sups(idIn)%fg)) then
             write(iftn,611) 'g', 'gmat'//c, 'gmat'//c
          end if
          write(iftn,612,advance='NO')
       end if
    end do
    if (n > 0) then
       write(iftn,690)
       write(lterm,*) n,'superelement(s) written to loadSupelMatrix.f90 and'// &
            &         ' supels.inc'
       write(lterm,*) 'Compile this into libvpmSupEls.so on Linux using:'
       write(lterm,*) '$ gfortran -fPIC -shared loadSupelMatrix.f90 '// &
            &         '-o libvpmSupEls.so'
    end if

    close(iftn)
    close(iinc)

600 format('subroutine LoadSupelMatrixFromCore (supId,dtype,data,ndata,stat)' &
         / ' !DEC$ ATTRIBUTES DLLEXPORT    :: LoadSupelMatrixFromCore', &
         / ' integer         , intent(in)  :: supId, ndata', &
         / ' character(len=1), intent(in)  :: dtype', &
         / ' double precision, intent(out) :: data(ndata)', &
         / ' integer         , intent(out) :: stat', &
         / ' include ''supels.inc''')
601 format(' integer, parameter :: seId',A1,' = ',I8, &
         / ' integer, parameter :: ndim',A1,' = ',I8 )
602 format(' double precision, parameter :: ',A1,'mat',A1, &
         & '(ndim',A1,'*',A,') = (/ &')
610 format(' if (supId == seId',A1,') then', &
         / '    select case (dtype)')
611 format('       case (''',A1,'''); data = ',A5,'; stat = size(',A5,')')
612 format('    end select', &
         / ' else')
690 format( &
         / '    stat = 0', &
         / ' end if', &
         / 'end subroutine LoadSupelMatrixFromCore')

  contains

    !> @brief Converts an integer (less that 37) to a character.
    function i2c(i) result(c)
      integer, intent(in) :: i
      character(len=1) :: c
      if (i < 1) then
         c = '0'
      else if (i < 10) then
         c = char(ichar('0')+i)
      else if (i < 37) then
         c = char(ichar('A')+i-10)
      else
         c = '*'
      end if
    end function i2c

    !> @brief Writes a matrix as a Fortran data statement.
    subroutine writeMatrix (ndat,data)
      integer         , intent(in) :: ndat
      double precision, intent(in) :: data(ndat)
      integer :: i, idat, ncont
      integer, parameter :: nValPrLine = 89 ! Set this based on compiler limit
      ncont = (ndat-1) / nValPrLine ! Number of continuation lines
      idat = 1
      do i = 1, ncont
         write(iinc,603,advance='NO') data(idat:idat+nValPrLine-1)
         write(iinc,"('&')")
         idat = idat + nValPrLine
      end do
      if (ndat > 1) then
         write(iinc,603,advance='NO') data(idat:ndat-1)
      end if
      write(iinc,604) data(ndat)
603   format(1P,99(E22.15,','))
604   format(1PE22.15,' /)')
    end subroutine writeMatrix

  end subroutine WriteSupEls2Ftn

end module initiateSupElTypeModule
