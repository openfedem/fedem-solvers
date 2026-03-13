!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file initiateTriadTypeModule.f90
!> @brief Initialization of triads from the solver input file.

!!==============================================================================
!> @brief Initialization of triads from the solver input file.

module initiateTriadTypeModule

  implicit none

contains

  !!============================================================================
  !> @brief Initializes triads with data from the input file.
  !>
  !> @param[in] infp File unit number for the solver input file
  !> @param[out] triads Array of all triads in the model
  !> @param[out] err Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 21 Jul 2000

  subroutine ReadTriads (infp,triads,err)

    use TriadTypeModule       , only : TriadType, NullifyTriad, DeallocateTriads
    use IdTypeModule          , only : initId, ReportInputError
    use inputUtilities        , only : iuGetNumberOfEntries, iuSetPosAtNextEntry
    use rotationModule        , only : orthonorm3
    use progressModule        , only : lterm
    use reportErrorModule     , only : AllocationError
    use reportErrorModule     , only : reportError, debugFileOnly_p
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getbool
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getint
    use TriadNamelistModule ! Defines all variables in the TRIAD namelist

    integer        , intent(in)  :: infp
    type(TriadType), pointer     :: triads(:)
    integer        , intent(out) :: err

    !! Local variables
    integer :: i, idIn, nTriads, stat
    logical :: allPri, allSec, allRest, allTriad
    logical :: allDef, allVel, allAcc, allForce, allHD
    logical :: needAddBC, ignoreIC, saveForceContribs

    !! --- Logic section ---

    nTriads = iuGetNumberOfEntries(infp,'&TRIAD',err)
    if (err /= 0) goto 900
    write(lterm,*) 'Number of &TRIAD =',nTriads

    call DeallocateTriads (triads)
    allocate(triads(nTriads),STAT=stat)
    if (stat /= 0) then
       err = AllocationError('ReadTriads 1')
       return
    end if

    !! Check if we need to account for the additional boundary conditions (BC=2)
    call ffa_cmdlinearg_getbool ('initEquilibrium',needAddBC)
    if (.not. needAddBC) then
       call ffa_cmdlinearg_getint ('numEigModes',stat)
       if (stat > 0) call ffa_cmdlinearg_getbool ('addBC_eigensolver',needAddBC)
    end if

    call ffa_cmdlinearg_getbool ('ignoreIC',ignoreIC)
    call ffa_cmdlinearg_getbool ('allPrimaryVars',allPri)
    call ffa_cmdlinearg_getbool ('allSecondaryVars',allSec)
    call ffa_cmdlinearg_getbool ('allRestartVars',allRest)
    call ffa_cmdlinearg_getbool ('allTriadVars',allTriad)
    call ffa_cmdlinearg_getbool ('allDefVars',allDef)
    call ffa_cmdlinearg_getbool ('allVelVars',allVel)
    call ffa_cmdlinearg_getbool ('allAccVars',allAcc)
    call ffa_cmdlinearg_getbool ('allHDVars',allHD)
    call ffa_cmdlinearg_getbool ('allForceVars',allForce)
    call ffa_cmdlinearg_getbool ('saveForceContributions',saveForceContribs)

    do idIn = 1, nTriads

       call NullifyTriad (triads(idIn))
       if (.not. iuSetPosAtNextEntry(infp,'&TRIAD')) then
          err = err - 1
          call ReportInputError('TRIAD',idIn)
          cycle
       end if

       call read_TRIAD (infp,stat)
       if (stat /= 0) then
          err = err - 1
          call ReportInputError('TRIAD',idIn)
          cycle
       end if

       call initId (triads(idIn)%id,id,extId,extDescr,stat)
       triads(idIn)%nDOFs = nDOFs
       triads(idIn)%ur = transpose(ur)
       call orthonorm3 (triads(idIn)%ur)

       if (nDOFS <= 0) cycle ! Ground linked triad

       allocate(triads(idIn)%urd(nDOFs), &
            &   triads(idIn)%urdd(nDOFs), &
            &   triads(idIn)%BC(nDOFs), &
            &   STAT=stat)
       if (stat /= 0) then
          err = AllocationError('ReadTriads 2')
          return
       end if

       if (ignoreIC) then
          triads(idIn)%urd  = 0.0_dp
          triads(idIn)%urdd = 0.0_dp
       else
          triads(idIn)%urd  = urd(1:nDOFs)
          triads(idIn)%urdd = urdd(1:nDOFs)
       end if

       if (needAddBC) then
          triads(idIn)%BC = BC(1:nDOFs)
       else
          triads(idIn)%BC = 1
          !! We still want the always fixed BCs (if any)
          do i = 1, nDOFs
             if (BC(i) <= 0) triads(idIn)%BC(i) = 0
          end do
       end if

       if (sysDir == 1) then
          !! Boundary conditions apply in initial local directions
          allocate(triads(idIn)%sysDirInG(3,3),STAT=stat)
          if (stat /= 0) then
             err = AllocationError('ReadTriads 2a')
             return
          end if
          triads(idIn)%sysDirInG = triads(idIn)%ur(:,1:3)
       else if (sysDir == 2) then
          !! Boundary conditions apply in updated local directions
          triads(idIn)%sysDirInG => triads(idIn)%urPrev(:,1:3)
       else if (sysDir == 3) then
          triads(idIn)%sysDirInG => triads(idIn)%ur(:,1:3)
       end if

       if (any(abs(dragParams) > 1.0e-15_dp)) then
          allocate(triads(idIn)%dragPar(3,3),STAT=stat)
          if (stat /= 0) then
             err = AllocationError('ReadTriads 4')
             return
          end if
          triads(idIn)%dragPar = dragParams
       end if

       !! Determine if primary variables (position matrix) will be saved
       triads(idIn)%savePos = allPri .or. allRest .or. allTriad .or. savePos > 0

       !! Determine which secondary triad variables will be saved. Settings in
       !! the solver input file are overruled by command-line arguments, if any.
       do i = 1, size(saveVar)
          triads(idIn)%saveVar(i) = saveVar(i) > 0 .or. allSec .or. allTriad
       end do
       if (allRest)  triads(idIn)%saveVar(1:3) = .true. ! Restart variables
       if (allVel)   triads(idIn)%saveVar(1) = .true. ! Global velocities
       if (allAcc)   triads(idIn)%saveVar(2) = .true. ! Global accelerations
       if (allForce) triads(idIn)%saveVar(3) = .true. ! Global forces
       if (allVel)   triads(idIn)%saveVar(4) = .true. ! Local velocities
       if (allAcc)   triads(idIn)%saveVar(5) = .true. ! Local acceleration
       if (allForce) triads(idIn)%saveVar(6) = .true. ! Local forces
       if (allDef)   triads(idIn)%saveVar(7) = .true. ! Global deformations
       triads(idIn)%saveVar(8) = saveForceContribs    ! For debugging only
       if (allHD)    triads(idIn)%saveVar(9) = .true. ! Hydrodynamics quantities

    end do

900 if (err < 0) call reportError (debugFileOnly_p,'ReadTriads')

  end subroutine ReadTriads


  !!============================================================================
  !> @brief Initializes more triad arrays after system initialization.
  !>
  !> @param[in] sam Data for managing system matrix assembly
  !> @param triads All triads in the model
  !> @param[out] err Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 21 Jul 2000

  subroutine initiateTriads2 (sam,triads,err)

    use SamModule      , only : SamType
    use TriadTypeModule, only : TriadType, transSysToGlob

    type(SamType)  , intent(in)    :: sam
    type(TriadType), intent(inout) :: triads(:)
    integer        , intent(out)   :: err

    !! Local variables
    integer :: idIn, iNod, nTriads, nDOF, firstDOF, lastDOF

    !! --- Logic section ---

    err = 0

    nTriads = size(triads)
    do idIn = 1, nTriads
       nDOF = triads(idIn)%nDOFs
       iNod = triads(idIn)%samNodNum
       if (nDOF < 1 .or. iNod < 1) cycle

       firstDOF = sam%madof(iNod)
       lastDOF  = sam%madof(iNod+1)-1

       triads(idIn)%sysDOF   => sam%madof(iNod)
       triads(idIn)%eqNumber => sam%meqn(firstDOF : lastDOF)
       triads(idIn)%urPrev   =  triads(idIn)%ur

       !! Initial conditions are assumed specified in the system directions
       call transSysToGlob (triads(idIn),triads(idIn)%urd)
       call transSysToGlob (triads(idIn),triads(idIn)%urdd)
    end do

  end subroutine InitiateTriads2

end module initiateTriadTypeModule
