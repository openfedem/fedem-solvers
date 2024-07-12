!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file initiateStressModule.f90
!> @brief Initialization of triads and superelements from the solver input file.

!!==============================================================================
!> @brief Initialization of triads and superelements from the solver input file.

module initiateStressModule

  implicit none

  private :: readSupEls, readTriads, initSupEls


contains

  !!============================================================================
  !> @brief Initializes a superelement and connected triads from the input file.
  !>
  !> @param[in] iSup Base ID of the superelement to be initialized
  !> @param[out] gvec Gravity vector
  !> @param[out] sup The superelement that that was initialized
  !> @param[out] triads Array of triads connected to the superelement
  !> @param[out] modelFileName Name of model file
  !> @param[in] iprint Print switch
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 2 Oct 2000

  subroutine readSolverData (iSup,gvec,sup,triads,modelFileName,iprint,lpu,ierr)

    use kindModule               , only : dp, lfnam_p
    use TriadTypeModule          , only : TriadType, writeObject
    use SupElTypeModule          , only : SupElType, writeObject
    use IdTypeModule             , only : ReportInputError
    use HeadingNameListModule    , only : Read_HEADING, modelFile
    use EnvironmentNameListModule, only : Read_ENVIRONMENT, gravity
    use inputUtilities           , only : iuCopyToScratch, iuSetPosAtNextEntry
    use fileUtilitiesModule      , only : findUnitNumber, getFileName
    use reportErrorModule        , only : reportError, error_p, debugFIleOnly_p

    integer         , intent(in)  :: iSup, iprint, lpu
    real(dp)        , intent(out) :: gvec(3)
    type(SupElType) , intent(out) :: sup
    type(TriadType) , pointer     :: triads(:)
    character(len=*), intent(out) :: modelFileName
    integer         , intent(out) :: ierr

    !! Local variables
    integer, pointer   :: triadIds(:)
    integer            :: i, infp
    character(lfnam_p) :: chname

    !! --- Logic section ---

    infp = findUnitNumber(50)
    open(infp,IOSTAT=ierr,STATUS='SCRATCH')
    if (ierr /= 0) then
       call reportError (error_p,'Unable to open temporary solver input file')
       goto 990
    end if

    call getFileName ('fsifile',chname)
    call iuCopyToScratch (infp,chname,ierr)
    if (ierr /= 0) goto 990

    rewind(infp)

    call Read_HEADING (infp,ierr)
    if (ierr /= 0) then
       call ReportInputError ('HEADING')
       goto 990
    end if

    modelFileName = modelFile

    if (iuSetPosAtNextEntry(infp,'&ENVIRONMENT')) then

       call Read_ENVIRONMENT (infp,ierr)
       if (ierr /= 0) then
          call ReportInputError ('ENVIRONMENT')
          goto 990
       end if

       gvec = gravity

    else ! No environment record, assume zero gravity
       rewind(infp)
       gvec = 0.0_dp
    end if

    allocate(triadIds(0),STAT=ierr)
    if (ierr == 0) call readSupEls (infp,iSup,triadIds,sup,ierr)
    if (ierr == 0) call readTriads (infp,triadIds,triads,ierr)
    if (ierr == 0) call initSupEls (triadIds,triads,sup,ierr)
    if (associated(triadIds)) deallocate(triadIds)

    close(infp)

    if (iprint < 3) goto 990

    call WriteObject (sup,lpu,1)
    write(lpu,'(/)')

    do i = 1, size(triads)
       call WriteObject (triads(i),lpu,1)
       write(lpu,'(/)')
    end do

990 continue
    if (ierr /= 0) call reportError (debugFileOnly_p,'readSolverData')

  end subroutine readSolverData


  !!============================================================================
  !> @brief Initializes triads with data from the input file.
  !>
  !> @param[in] infp File unit number for the solver input file
  !> @param[in] triadIds Array of user Ids of the triads to consider
  !> @param[out] triads Array of initialized triad objects
  !> @param[out] err Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 29 Sep 2000

  subroutine readTriads (infp,triadIds,triads,err)

    use TriadTypeModule    , only : TriadType, nullifyTriad
    use TriadNamelistModule, only : read_TRIAD, id, extId, extDescr, nDOFs, ur
    use IdTypeModule       , only : initId, ReportInputError
    use inputUtilities     , only : iuGetNumberOfEntries, iuSetPosAtNextEntry
    use progressModule     , only : lterm
    use reportErrorModule  , only : allocationError
    use reportErrorModule  , only : reportError, debugFileOnly_p

    integer         , intent(in)  :: infp, triadIds(:)
    type(TriadType) , pointer     :: triads(:)
    integer         , intent(out) :: err

    !! Local variables
    integer :: i, idIn, indexTr, nTriads, nTriadIds, stat

    !! --- Logic section ---

    nTriads = iuGetNumberOfEntries(infp,'&TRIAD',err)
    if (err /= 0) goto 900

    write(lterm,"(15X,'Number of &TRIAD  =',I6)") nTriads

    nTriadIds = min(nTriads,size(triadIds))
    allocate(triads(nTriadIds),STAT=stat)
    if (stat /= 0) then
       err = AllocationError('readTriads')
       return
    end if

    indexTr = 0
    do idIn = 1, nTriads

       if (.not. iuSetPosAtNextEntry(infp,'&TRIAD')) then
          err = err - 1
          call ReportInputError ('TRIAD',idIn)
          cycle
       end if

       call read_TRIAD (infp,stat)
       if (stat /= 0) then
          err = err - 1
          call ReportInputError ('TRIAD',idIn)
          cycle
       end if

       !! Check if this triad is within the list of triadIds
       do i = 1, nTriadIds
          if (triadIds(i) == id) goto 100
       end do
       cycle

100    indexTr = indexTr + 1
       call nullifyTriad (triads(indexTr))
       call initId (triads(indexTr)%id,id,extId,extDescr,stat)
       triads(indexTr)%nDOFs = nDOFs
       triads(indexTr)%ur = transpose(ur)

    end do
    do while (indexTr < nTriadIds)
       indexTr = indexTr + 1
       call nullifyTriad (triads(indexTr))
    end do

900 if (err < 0) call reportError (debugFileOnly_p,'readTriads')

  end subroutine readTriads


  !!============================================================================
  !> @brief Initializes the superelement with data from the input file.
  !>
  !> @param[in] infp File unit number for the solver input file
  !> @param[in] baseId Base Id of the superelement to consider
  !> @param[out] triads Array of user Ids of the triads to consider
  !> @param[out] sup The superelement object to consider
  !> @param[out] err Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 29 Sep 2000

  subroutine readSupEls (infp,baseId,triads,sup,err)

    use SupElTypeModule    , only : SupElType, nullifySupEl
    use SupElNamelistModule, only : read_SUP_EL, read_TRIAD_UNDPOS
    use SupElNamelistModule, only : reportTriadError
    use SupElNamelistModule, only : id, extId, extDescr, supPos
    use SupElNamelistModule, only : maxTriads_p, numTriads, triadIds, numGenDOFs
    use SupElNamelistModule, only : supElId, triadId, undPosInSupElSystem
    use IdTypeModule       , only : initId, ReportInputError
    use inputUtilities     , only : iuGetNumberOfEntries, iuSetPosAtNextEntry
    use progressModule     , only : lterm
    use reportErrorModule  , only : allocationError
    use reportErrorModule  , only : reportError, error_p, debugFileOnly_p

    integer        , intent(in)  :: infp, baseId
    integer        , pointer     :: triads(:)
    type(SupElType), intent(out) :: sup
    integer        , intent(out) :: err

    !! Local variables
    integer           :: i, idIn, j, nSupEl, stat
    character(len=64) :: errmsg

    !! --- Logic section ---

    call nullifySupEl (sup)
    nSupEl = iuGetNumberOfEntries(infp,'&SUP_EL',err)
    if (err /= 0) goto 900

    write(lterm,"(15X,'Number of &SUP_EL =',I6)") nSupEl

    do idIn = 1, nSupEl

       if (.not. iuSetPosAtNextEntry(infp,'&SUP_EL')) then
          err = err - 1
          call ReportInputError ('SUP_EL',idIn)
          cycle
       end if

       call read_SUP_EL (infp,stat)
       if (stat /= 0) then
          err = err - 1
          call ReportInputError ('SUP_EL',idIn)
       else if (id == baseId) then
          goto 200
       end if

    end do

    err = -999
    write(errMsg,*) 'baseID',baseId,' was not found on the solver input file'
    call reportError (error_p,errMsg,addString='readSupEls')
    return

200 call initId (sup%id,id,extId,extDescr,stat)
    sup%nExtNods = numTriads
    sup%supTr = transpose(supPos)
    sup%supTrInit = sup%supTr

    if (numTriads > maxTriads_p) then
       err = err - 1
       write(*,*) '*** To many triads',numTriads,' for superelement',id
       write(*,*) '    Array "triadIds" is overwritten! Watch out...!'
    end if
    if (err < 0) goto 900

    if (associated(triads)) deallocate(triads)
    allocate(triads(numTriads), &
         &   sup%genDOFs, &
         &   sup%triads(numTriads), &
         &   sup%TrUndeformed(3,4,numTriads), STAT=err)
    if (err /= 0) then
       err = AllocationError('readSupEls')
       return
    end if

    sup%genDOFs%nDOFs = numGenDOFs
    sup%genDOFs%firstDOF = 0
    sup%genDOFs%samNodNum = 0
    nullify(sup%genDOFs%sysDOF)
    nullify(sup%genDOFs%BC)
    nullify(sup%genDOFs%alpha1)
    nullify(sup%genDOFs%alpha2)
    nullify(sup%genDOFs%ur)
    nullify(sup%genDOFs%urPrev)
    nullify(sup%genDOFs%urd)
    nullify(sup%genDOFs%urdd)

    do i = 1, numTriads

       if (.not. iuSetPosAtNextEntry(infp,'&TRIAD_UNDPOS')) then
          err = err - 1
          call ReportTriadError (sup%id,i)
          cycle
       end if

       call read_TRIAD_UNDPOS (infp,stat)
       if (stat /= 0) then
          err = err - 1
          call ReportTriadError (sup%id,i)
          cycle
       else if (supElId /= id) then
          err = err - 1
          call ReportTriadError (sup%id,triadId,supElId)
          cycle
       end if

       !! The UNDPOS records are not necessarily in the same order
       !! as the superelement triads. Therefore we must do a search.
       do j = 1, numTriads
          if (triadId == triadIds(j)) goto 100
       end do
       err = err - 1
       call ReportTriadError (sup%id,triadId,supElId, &
            'This Triad is not connected to this Part')
       cycle

100    continue
       triads(j) = triadIds(j)
       sup%TrUndeformed(:,:,j) = transpose(undPosInSupElSystem)

    end do

900 if (err < 0) call reportError (debugFileOnly_p,'readSupEls')

  end subroutine readSupEls


  !!============================================================================
  !> @brief Initializes the superelement with data from triads.
  !>
  !> @param[in] infp File unit number for the solver input file
  !> @param[in] triadIds Array of user Ids of the triads to consider
  !> @param[in] triads Array of triad objects to consider
  !> @param sup The superelement object to consider
  !> @param[out] err Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 29 Sep 2000

  subroutine initSupEls (triadIds,triads,sup,err)

    use TriadTypeModule  , only : TriadType, GetPtrToId, dp
    use SupElTypeModule  , only : SupElType
    use IdTypeModule     , only : getId
    use progressModule   , only : lterm
    use reportErrorModule, only : allocationError, reportError, debugFileOnly_p

    integer        , intent(in)    :: triadIds(:)
    type(TriadType), intent(in)    :: triads(:)
    type(SupElType), intent(inout) :: sup
    integer        , intent(out)   :: err

    !! Local variables
    integer :: i, nDof, stat

    !! --- Logic section ---

    err = 0
    write(lterm,"(15X,'Part',A,'; nTriads =',I4,'  nGenDofs =',I4)") &
         &      trim(getId(sup%id)), sup%nExtNods, sup%genDOFs%nDOFs

    nDof = 0
    do i = 1, sup%nExtNods
       sup%triads(i)%p => GetPtrToId(triads,triadIds(i))
       if (associated(sup%triads(i)%p)) then
          sup%triads(i)%firstDOF = nDof + 1
          nDof = nDof + sup%triads(i)%p%nDOFs
       else
          err = err - 1
       end if
    end do
    sup%nTotDofs = nDof + sup%genDOFs%nDOFs
    sup%genDOFs%firstDOF = nDof + 1

    allocate(sup%finit(sup%nTotDofs), STAT=stat)
    if (stat /= 0) then
       err = AllocationError('initSupEls: sup%finit')
       return
    end if

    if (sup%genDOFs%nDOFs > 0) then
       allocate(sup%genDOFs%ur(sup%genDOFs%nDOFs), STAT=stat)
       if (stat /= 0) then
          err = AllocationError('initSupEls: sup%genDOFs%ur')
          return
       end if
       sup%genDOFs%ur = 0.0_dp
    end if

    if (err < 0) call reportError (debugFileOnly_p,'initSupEls')

  end subroutine initSupEls

end module initiateStressModule
