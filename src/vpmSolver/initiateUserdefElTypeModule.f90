!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file initiateUserdefElTypeModule.f90
!> @brief Initialization of user-defined elements from the solver input file.

!!==============================================================================
!> @brief Initialization of user-defined elements from the solver input file.

module initiateUserdefElTypeModule

  implicit none

contains

  !!============================================================================
  !> @brief Initializes user-defined elements with data from the input file.
  !>
  !> @param[in] infp File unit number for the solver input file
  !> @param[in] triads Array of all triads in the model
  !> @param[in] engines All general functions in the model
  !> @param[out] elms Array of all user-defined elements in the model
  !> @param[out] err Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 31 Mar 2014

  subroutine ReadUserdefEls (infp,triads,engines,elms,err)

    use IdTypeModule       , only : ldesc_p, initId, StrId, ReportInputError
    use TriadTypeModule    , only : TriadType, GetPtrToId, dp
    use FunctionTypeModule , only : EngineType, GetPtrToId
    use UserdefElTypeModule, only : UserdefElType, NullifyUserdefEl
    use UserdefElTypeModule, only : DeallocateUdElms
    use SupElTypeModule    , only : InitiateHyDyn
    use inputUtilities     , only : iuGetNumberOfEntries, iuSetPosAtNextEntry
    use reportErrorModule  , only : allocationError, internalError
    use reportErrorModule  , only : debugFileOnly_p, reportError
    use IdTypeModule       , only : nullifyId
    use progressModule     , only : lterm

    integer            , intent(in)  :: infp
    type(TriadType)    , intent(in)  :: triads(:)
    type(EngineType)   , intent(in)  :: engines(:)
    type(UserdefElType), pointer     :: elms(:)
    integer            , intent(out) :: err

    !! Define the USER_EL namelist
    integer, parameter :: maxTriads_p = 1000
    integer, parameter :: maxIpar_p = 100, maxRpar_p = 100, maxVar_p = 10
    integer            :: id, extId(10), eType, numTriads, triadIds(maxTriads_p)
    integer            :: nipar, nrpar, nvar, ipar(maxIpar_p), idVar(maxVar_p)
    real(dp)           :: stiffScale, massScale, alpha1, alpha2
    real(dp)           :: rpar(4+maxRpar_p), morison(7)
    character(ldesc_p) :: extDescr

    namelist /USER_EL/ id, extId, extDescr, eType, numTriads, triadIds, &
         &             stiffScale, massScale, alpha1, alpha2, &
         &             nipar, ipar, nrpar, rpar, nvar, idVar, morison

    !! Local variables
    logical  :: doHydro
    integer  :: i, idIn, nDOF, nUserEl, nFnod, stat
    real(dp) :: MorPar(11)

    !! --- Logic section ---

    nUserEl = iuGetNumberOfEntries(infp,'&USER_EL',err)
    if (err /= 0) goto 900
    write(lterm,*) 'Number of &USER_EL =',nUserEl

    call DeallocateUdElms (elms)
    allocate(elms(nUserEl),STAT=stat)
    if (stat /= 0) then
       err = AllocationError('ReadUserdefEls 1')
       return
    end if

    nFnod = 0
    do idIn = 1, size(triads)
       if (triads(idIn)%inFluid > nFnod) nFnod = triads(idIn)%inFluid
    end do

    do idIn = 1, nUserEl

       call nullifyUserdefEl (elms(idIn))
       if (.not. iuSetPosAtNextEntry(infp,'&USER_EL',stat)) then
          err = err - 1
          call ReportInputError ('USER_EL',idIn)
          if (stat < 0) then
             exit
          else
             cycle
          end if
       end if

       !! Default values
       id=0; extId=0; extDescr=''
       eType=0; numTriads=0; triadIds=0
       stiffScale=1.0_dp; massScale=1.0_dp; alpha1=0.0_dp; alpha2=0.0_dp
       nipar=0; ipar=0; nrpar=0; rpar=0.0_dp; nvar=0; idVar=0; morison=0.0_dp

       read(infp,nml=USER_EL,iostat=stat)
       if (stat /= 0) then
          err = err - 1
          call ReportInputError ('USER_EL',idIn)
          cycle
       end if

       call initId (elms(idIn)%id,id,extId,extDescr,stat)
       elms(idIn)%type = eType

       if (numTriads > maxTriads_p) then
          err = internalError('ReadUserdefEls: Too many element triads'// &
               &              StrId(numTriads))
          return
       else if (nipar > maxIpar_p) then
          err = internalError('ReadUserdefEls: Too many integer parameters'// &
               &              StrId(nipar))
          return
       else if (nrpar > maxRpar_p) then
          err = internalError('ReadUserdefEls: Too many Real parameters'// &
               &              StrId(nrpar))
          return
       end if

       deallocate(elms(idIn)%triads,elms(idIn)%TrUndeformed,elms(idIn)%engines)
       allocate(elms(idIn)%triads(numTriads), &
            &   elms(idIn)%TrUndeformed(3,4,numTriads), &
            &   elms(idIn)%iwork(3+nipar+numTriads), &
            &   elms(idIn)%rwork(4+nrpar), &
            &   elms(idIn)%engines(nvar),STAT=stat)
       if (stat /= 0) then
          err = AllocationError('ReadUserdefEls 2')
          return
       end if

       doHydro = morison(4) > 0.0_dp ! Perform hydrodynamic force calculation

       !! Initialize the triad connections and
       !! determine dimension of the element matrices
       nDOF = 0
       stat = err
       do i = 1, numTriads
          elms(idIn)%triads(i)%p => GetPtrToId(triads,triadIds(i))
          if (associated(elms(idIn)%triads(i)%p)) then
             elms(idIn)%TrUndeformed(:,:,i) = elms(idIn)%triads(i)%p%ur
             elms(idIn)%triads(i)%firstDOF = nDOF + 1
             nDOF = nDOF + elms(idIn)%triads(i)%p%nDOFs
             if (doHydro .and. elms(idIn)%triads(i)%p%inFluid == 0) then
                nFnod = nFnod + 1
                elms(idIn)%triads(i)%p%inFluid = nFnod
             end if
          else
             err = err - 1
             elms(idIn)%triads(i)%firstDOF = 0
          end if
       end do
       elms(idIn)%nTotDOFs = nDOF
       if (err < stat) then
          call ReportInputError ('USER_EL',idIn,elms(idIn)%id, &
               &                 'Invalid triad connection')
       end if

       !! Initialize engine connections for any state-dependent variables
       stat = err
       do i = 1, nvar
          if (idVar(i) > 0) then
             elms(idIn)%engines(i)%p => GetPtrToId(engines,idVar(i))
             if (.not. associated(elms(idIn)%engines(i)%p)) err = err - 1
          else
             nullify(elms(idIn)%engines(i)%p)
          end if
       end do
       if (err < stat) then
          call ReportInputError ('USER_EL',idIn,elms(idIn)%id, &
               &                 'Invalid engine connection')
       end if

       !! Initialize the work arrays. Notice that the first three locations of
       !! iwork and the first four locations of rwork have fixed interpretation.
       !! The rest will depend on the element plugin.
       elms(idIn)%iwork(1) = 3+nipar+numTriads ! Total length of iwork
       elms(idIn)%iwork(2) = 4+nrpar           ! Total length of rwork
       elms(idIn)%iwork(3) = nvar   ! Number of state-dependent variables
       elms(idIn)%rwork(1) = alpha1 ! Mass-proportional damping factor
       elms(idIn)%rwork(2) = alpha2 ! Stiffness-proportional damping factor
       elms(idIn)%rwork(3) = stiffScale ! Constant stiffness scaling factor
       elms(idIn)%rwork(4) = massScale  ! Constant mass scaling factor
       if (nipar > 0) elms(idIn)%iwork(4:3+nipar) = ipar(1:nipar)
       if (nrpar > 0) elms(idIn)%rwork(5:4+nrpar) = rpar(1:nrpar)
       !! Store the firstDOF indices at the end of the integer work array, such
       !! that user-defined elements can handle a varying number triad DOFs
       do i = 1, numTriads
          elms(idIn)%iwork(3+nipar+i) = elms(idIn)%triads(i)%firstDOF
       end do

       !! Allocate and initialize the element matrices
       allocate(elms(idIn)%Q   (nDof), &
            &   elms(idIn)%Fs  (nDof), &
            &   elms(idIn)%Fd  (nDof), &
            &   elms(idIn)%Fi  (nDof), &
            &   elms(idIn)%Kmat(nDof,nDof), &
            &   elms(idIn)%Cmat(nDof,nDof), &
            &   elms(idIn)%Mmat(nDof,nDof), &
            &   elms(idIn)%Nmat(nDof,nDof), &
            &   STAT=stat)
       if (stat /= 0) then
          err = AllocationError('ReadUserdefEls 3')
          return
       end if

       elms(idIn)%Q    = 0.0_dp
       elms(idIn)%Fs   = 0.0_dp
       elms(idIn)%Fd   = 0.0_dp
       elms(idIn)%Fi   = 0.0_dp
       elms(idIn)%Kmat = 0.0_dp
       elms(idIn)%Cmat = 0.0_dp
       elms(idIn)%Mmat = 0.0_dp
       elms(idIn)%Nmat = 0.0_dp

       if (doHydro) then

          !! Hydrodynamic parameters are specified
          allocate(elms(idIn)%hydyn, STAT=stat)
          if (stat /= 0) then
             err = AllocationError('ReadUserdefEls 4')
             return
          end if

          MorPar(1:4) = morison(1:4) ! Ca, Cm, Cd and drag-diameter (D)
          MorPar(5) = morison(4) ! Assume equal buoyancy and drag diameters
          MorPar(6) = morison(5) ! Axial drag coefficient
          MorPar(7) = 0.0_dp
          MorPar(8) = morison(6) ! Axial added mass coefficient (Ca)
          MorPar(9) = morison(7) ! Axial added mass coefficient (Cm)
          MorPar(10:11) = 0.0_dp
          call initiateHyDyn (elms(idIn)%hydyn,MorPar)

       end if

    end do

900 if (err < 0) call reportError (debugFileOnly_p,'ReadUserdefEls')

  end subroutine ReadUserdefEls


  !!============================================================================
  !> @brief Initializes the user-defined elements after reading the input file.
  !>
  !> @param[in] gravity Gravitation vector
  !> @param elms Array of all user-defined elements in the model
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 31 Mar 2014

  subroutine InitiateUserdefEls (gravity,elms,ierr)

    use UserDefElTypeModule    , only : UserdefElType, dp
    use UserDefElRoutinesModule, only : InitializeUDE
    use FiUserElmInterface     , only : Fi_UDE_Init
    use IdTypeModule           , only : ldesc_p
#ifdef FT_DEBUG
    use fileUtilitiesModule    , only : getDBGfile
    use dbgUnitsModule         , only : dbgUDE
#endif
    use reportErrorModule      , only : getErrorFile
    use reportErrorModule      , only : reportError, debugFileOnly_p
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_getstring

    real(dp)           , intent(in)    :: gravity(3)
    type(UserdefElType), intent(inout) :: elms(:)
    integer            , intent(out)   :: ierr

    !! Local variables
    integer            :: i, lpu
    character(ldesc_p) :: signature
    character(len=520) :: plugin

    !! --- Logic section ---

    ierr = 0
    if (size(elms) < 1) return

    call ffa_cmdlinearg_getstring ('plugin',plugin)
    call Fi_UDE_Init (plugin,gravity(1),signature,ierr)
    if (ierr < 0) goto 900

    lpu = getErrorFile()
    if (size(elms) == 1) then
       write(lpu,601) trim(signature)
    else
       write(lpu,600) size(elms),trim(signature)
    end if
600 format(/5X,'This model contains',I6,' user-defined elements.'/5X,A)
601 format(/5X,'This model contains one user-defined element.'/5X,A)

#ifdef FT_DEBUG
    dbgUDE = getDBGfile(66,'userelm.dbg')
#endif

    do i = 1, size(elms)
       call InitializeUDE (elms(i)%id,elms(i)%type,elms(i)%nTotDofs, &
            &              elms(i)%triads,elms(i)%iwork,elms(i)%rwork,ierr)
       if (ierr < 0) goto 900
    end do

    return

900 call reportError (debugFileOnly_p,'InitiateUserdefEls')

  end subroutine InitiateUserdefEls

end module initiateUserdefElTypeModule
