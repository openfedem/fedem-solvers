!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file samReducerModule.f90
!> @brief Initialisation of the SAM data structure for the FE part reducer.

!!==============================================================================
!> @brief Initialisation of the SAM data structure for the FE part reducer.

module SamReducerModule

  implicit none

  private :: BeamPinConstraints


contains

  !!============================================================================
  !> @brief Establishes constraint equations associated with beam pin flags.
  !>
  !> @param[in] madof Matrix of accumulated DOFs
  !> @param[in] mpbeam Matrix of IDs for beam elements with pin flags
  !> @param[in] mpmnpc Matrix of pointers to MNPCs
  !> @param     mmnpc Matrix of nodal point correspondances
  !> @param[in] msc Matrix of status codes
  !> @param     minex Matrix of internal to external node number mapping
  !> @param     nceq Number of constraint equations
  !> @param     nmmceq Number of elements in MMCEQ
  !> @param     mpmceq Matrix of pointers to MCEQs
  !> @param     mmceq Matrix of constraint equation definitions
  !> @param     ttcc Table of constraint equation coefficients
  !> @param[in] nxnod Number of extra nodes to account for beam pin flags
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @details This subroutine appends the constraint equations accounting for
  !> beam pin flags to the linear coupling information arrays, MPMCEQ, MMCEQ
  !> and TTCC of sammodule::samtype. The element topology array (MMNPC)
  !> is also modified for the affected beam elements.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 16 Sep 2002

  subroutine BeamPinConstraints (madof, mpbeam, mpmnpc, mmnpc, msc, minex, &
       &                         nceq, nmmceq, mpmceq, mmceq, ttcc, &
       &                         nxnod, lpu, ierr)

    use kindModule             , only : dp, epsDiv0_p
    use manipMatrixModule      , only : cross_product
    use reportErrorModule      , only : internalError, reportError
    use reportErrorModule      , only : warning_p, debugFileOnly_p
    use FFlLinkHandlerInterface, only : ffl_getCoor, ffl_getPinFlags
    use FFlLinkHandlerInterface, only : ffl_getElmId

    integer , intent(in)    :: madof(:), mpbeam(:), mpmnpc(:), msc(:)
    integer , intent(inout) :: mmnpc(:), minex(:)
    integer , intent(inout) :: nceq, nmmceq, mpmceq(:), mmceq(:)
    real(dp), intent(inout) :: ttcc(:)
    integer , intent(in)    :: nxnod, lpu
    integer , intent(out)   :: ierr

    !! Local variables
    integer  :: i, j, k, iel, iceq, ipc, ipMN, neccpb, nmceq
    integer  :: maxSiz, iDof, iNod, depDof, depNod, depLoc, PF(2)
    real(dp) :: C, e(3), xg(5), yg(5), zg(5), Tgl(3,3)
    real(dp), parameter :: zero_p = 1.0e-7_dp

    !! --- Logic section ---

    maxSiz = min(size(mmceq),size(ttcc))
    depNod = size(minex) - nxnod
    iceq   = nceq+1
    neccpb = 0

    do iel = 1, size(mpbeam)

       !! Get the beam pin flags
       call ffl_getPinFlags (PF(1),PF(2),mpbeam(iel),ierr)
       if (ierr < 0) goto 990

       !! Get element nodal coordinates
       call ffl_getCoor (xg,yg,zg,mpbeam(iel),ierr)
       if (ierr < 0) goto 990

       !! Compute global-to-local transformation matrix for this element
       !! The below is similar as in the Femlib subroutine DCOS30, but is
       !! reproduced here to avoid dependency to Femlib from this module.
       e(1) = Xg(2) - Xg(1)
       e(2) = Yg(2) - Yg(1)
       e(3) = Zg(2) - Zg(1)
       Tgl(1,:) = e / sqrt(e(1)*e(1) + e(2)*e(2) + e(3)*e(3))
       e(1) = Xg(3) - Xg(1)
       e(2) = Yg(3) - Yg(1)
       e(3) = Zg(3) - Zg(1)
       e    = cross_product(e,Tgl(1,:))
       Tgl(2,:) = e / sqrt(e(1)*e(1) + e(2)*e(2) + e(3)*e(3))
       Tgl(3,:) = cross_product(Tgl(1,:),Tgl(2,:))

       ipMN = mpmnpc(mpbeam(iel)) - 1
       do i = 1, 2
          if (PF(i) < 1) cycle

          !!TODO: Remove this and replace by correct eccentricity handling
          e(1) = xg(3+i) - xg(i)
          e(2) = yg(3+i) - yg(i)
          e(3) = zg(3+i) - zg(i)
          if (e(1)*e(1) + e(2)*e(2) + e(3)*e(3) > epsDiv0_p) then
             neccpb = neccpb + 1
             write(lpu,600) ffl_getElmId(mpbeam(iel)), i, PF(i), e
600          format('  ** Warning: Beam element', I10, &
                  & ' has end release at the eccentric end', I2 &
                  / '              pin flag =', I6, &
                  & ',  eccentricity vector =', 1P3E13.5 )
          end if

          !! Substitute the independent node with the new dependent node
          !! in the element connectivity for this beam element
          depNod = depNod + 1
          iNod = mmnpc(ipMN+i)
          mmnpc(ipMN+i) = depNod
          minex(depNod) = -iNod ! Identify the new dependent node

          do depDof = madof(depNod), madof(depNod+1)-1
             if (msc(depDof) > 0) cycle ! This dependent DOF is released (free)

             !! Add a constraint equation in TTCC for current dependent DOF
             depLoc = depDof - madof(depNod)
             iDof   = madof(iNod) + 3*(depLoc/3) - 1
             depLoc = 1 + mod(depLoc,3)
             nmmceq = nmmceq + 1
             if (nmmceq > maxSiz) goto 910
             mmceq(nmmceq) = depDof
             ttcc(nmmceq)  = 0.0_dp

             !! The dependent DOF is defined in the element coordinate system.
             !! Thus, it might be coupled to (up to) three independent DOFs.
             do j = 1, 3
                C = Tgl(depLoc,j)
                iDof = iDof + 1
                if (abs(C) <= zero_p) cycle ! Ignore terms much less than one

                if (msc(iDof) > 0) then

                   !! Add a constraint equation forcing current dependent DOF
                   !! to be equal the current independent DOF
                   nmmceq = nmmceq + 1
                   if (nmmceq > maxSiz) goto 910
                   mmceq(nmmceq) = iDof
                   ttcc(nmmceq)  = C

                else

                   !! Current DOF is dependent in another constraint equation
                   do k = 1, nceq
                      if (mmceq(mpmceq(k)) == iDof) goto 100
                   end do
                   goto 900

                   !! Make a copy the RHS of the other contraint equation
                   !! scaled by the coefficient C
100                continue
                   ipc   = mpmceq(k)
                   nmceq = mpmceq(k+1) - ipc - 1
                   if (nmmceq+nmceq > maxSiz) goto 910
                   mmceq(nmmceq+1:nmmceq+nmceq) = mmceq(ipc+1:ipc+nmceq)
                   ttcc(nmmceq:nmmceq+nmceq)    = C*ttcc(ipc:ipc+nmceq)
                   nmmceq = nmmceq + nmceq

                end if

             end do

             iceq = iceq + 1
             if (iceq > size(mpmceq)) goto 920
             mpmceq(iceq) = nmmceq + 1

          end do
       end do
    end do

    nceq = iceq-1
    if (neccpb > 0) then
       call reportError (warning_p, &
            'This part contains eccentric beam elements with end release(s).', &
            'This may be incorrectly accounted for in the current version.', &
            'Please contact Fedem support.', &
            addString='BeamPinConstraints')
    end if

    return

900 ierr = internalError('BeamPinConstraints: Invalid SAM datastructure')
    return
910 ierr = internalError('BeamPinConstraints: mmceq/ttcc too small')
    return
920 ierr = internalError('BeamPinConstraints: mpmceq too small')
    return
990 call reportError (debugFileOnly_p,'BeamPinConstraints')

  end subroutine BeamPinConstraints


  !!============================================================================
  !> @brief Initializes the SAM data structure for fedem_reducer.
  !>
  !> @param[out] sam Data for managing system matrix assembly
  !> @param[out] tncoor Table of nodal coordinates of the FE model
  !> @param[in] linearSolve If .true., we are doing a direct linear solve
  !> @param[in] ipsw Print switch for debug output
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 20 Sep 2000

  subroutine initiateSAM (sam,tncoor,linearSolve,ipsw,lpu,ierr)

    use SamModule              , only : SamType, dp, version_p
    use SamModule              , only : nullifySam, check, writeObject
    use RigidModule            , only : MultiPointConstraints, Rigid3D
    use RigidModule            , only : checkRgdElements, checkRbarElements
    use RigidModule            , only : krigid, krbar, kwavgm
    use ProgressModule         , only : writeProgress
    use allocationModule       , only : reAllocate, logAllocMem, doLogMem
    use reportErrorModule      , only : allocationError, reportError
    use reportErrorModule      , only : error_p, warning_p, debugFileOnly_p
    use FFlLinkHandlerInterface, only : ffl_getSize, ffl_getNodes, ffl_getTopol
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_isTrue

    logical      , intent(in)  :: linearSolve
    integer      , intent(in)  :: ipsw, lpu
    type(SamType), intent(out) :: sam
    real(dp)     , pointer     :: tncoor(:,:)
    integer      , intent(out) :: ierr

    !! Local variables
    character(len=128)    :: errMsg
    integer               :: imcex, imceq, ipceq, itceq
    integer               :: ncex, nmcex, nceq, nmceq
    integer               :: npbeam, nrgd, nrbar, nwavgm
    integer , allocatable :: mpbeam(:), mprgd(:), mprbar(:), mpwavgm(:)
    integer , allocatable :: iscr(:), meldep(:)
    real(dp), allocatable :: rscr(:)

    !! --- Logic section ---

    call nullifySam (sam)

    !! Initialize the Matrix of Parameters
    call reAllocate ('initiateSAM',sam%mpar,50,ierr)
    if (ierr < 0) return

    sam%mpar = 0
    sam%mpar(50) = version_p

    sam%nnod   => sam%mpar(1)  !! N. of active nodes (nanod)
    sam%nel    => sam%mpar(2)  !! N. of elements
    sam%ndof   => sam%mpar(3)  !! N. of degrees of freedom
    sam%ndof1  => sam%mpar(4)  !! N. of 1-status dofs
    sam%ndof2  => sam%mpar(5)  !! N. of 2-status dofs
    sam%nspdof => sam%mpar(6)  !! N. of specified Dofs
    sam%nceq   => sam%mpar(7)  !! N. of constraint equations
    sam%nsdof  => sam%mpar(8)  !! N. of supressed dofs
    sam%npdof  => sam%mpar(9)  !! N. of prescribed dofs
    sam%nddof  => sam%mpar(10) !! N. of dependent dofs
    sam%neq    => sam%mpar(11) !! N. of equations
    sam%nsky   => sam%mpar(12) !! N. of elements in SKYline array
    sam%nmmnpc => sam%mpar(15) !! Number of elements in mmnpc
    sam%nmmceq => sam%mpar(16) !! Number of elements in mmceq

    !! Get some model size parameters
    call ffl_getSize (sam%nnod, sam%nel, sam%ndof, sam%nmmnpc, sam%mpar(17), &
         &            sam%mpar(23), npbeam, nrgd, nrbar, nwavgm, &
         &            sam%mpar(27), sam%mpar(28), ierr)
    if (ierr < 0) goto 900

    !! Allocate some global nodal arrays
    call writeProgress (' --> Creating nodal arrays')
    call reAllocate ('initiateSAM',tncoor   ,sam%nnod,3,ierr)
    call reAllocate ('initiateSAM',sam%madof,sam%nnod+1,ierr)
    call reAllocate ('initiateSAM',sam%minex,sam%nnod  ,ierr)
    call reAllocate ('initiateSAM',sam%mnnn ,sam%nnod  ,ierr)
    call reAllocate ('initiateSAM',sam%msc  ,sam%ndof  ,ierr)
    if (ierr < 0) return

    !! Get the global nodal data
    call ffl_getNodes (sam%nnod, sam%ndof, &
         &             sam%madof(1), sam%minex(1), sam%mnnn(1), sam%msc(1), &
         &             tncoor(1,1), tncoor(1,2), tncoor(1,3), ierr)
    if (ierr /= 0) then
       call reportError (error_p,'Failure creating global nodal arrays', &
            &            addString='ffl_getNodes')
       goto 900
    end if

    !! Allocate some global topology arrays
    call writeProgress (' --> Creating element topology arrays')
    call reAllocate ('initiateSAM',sam%mpmnpc,sam%nel+1 ,ierr)
    call reAllocate ('initiateSAM',sam%mmnpc ,sam%nmmnpc,ierr)
    call reAllocate ('initiateSAM',sam%melcon,sam%nel   ,ierr)
    if (ierr < 0) return

    allocate(mpbeam(npbeam) , &
         &   mpwavgm(nwavgm), &
         &   mprbar(nrbar)  , &
         &   mprgd(nrgd)    , &
         &   meldep(nrgd)   , STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('initiateSAM 2')
       return
    else if (doLogMem) then
       call logAllocMem ('initiateSAM',0,npbeam+nwavgm+nrbar+2*nrgd,4)
    end if

    !! Get the global topology data.
    !! Note that ffl_getTopol may return a smaller sam%nmmnpc than ffl_getSize.
    !! Thus, we cannot use size(sam%mmnpc) when writing sam%mmnpc in saveSAM.
    call ffl_getTopol (sam%nel, sam%nmmnpc, &
         &             sam%melcon(1), sam%mmnpc(1), sam%mpmnpc(1), &
         &             mpbeam, mprgd, mprbar, mpwavgm, ierr)
    if (ierr < 0) then
       call reportError (error_p,'Failure creating global topology arrays', &
            &            addString='ffl_getTopol')
       goto 900
    end if

    !! Fix all external DOFs to directly solve for some external load
    if (linearSolve) call fixExternalDOFs (sam%msc)

    if (nrgd > 0) then

       !! Check RGD elements for common nodes and establish chain connectivity
       call checkRgdElements (mprgd, mprbar, mpwavgm, meldep, &
            &                 sam%minex, sam%mpmnpc, sam%mmnpc, sam%mnnn, ierr)
       if (ierr /= 0) goto 900

    end if
    if (nrbar > 0) then

       !! Check RBAR elements for common nodes
       call checkRbarElements (mprbar, mprgd, &
            &                  sam%minex, sam%mpmnpc, sam%mmnpc, &
            &                  sam%madof, sam%msc, sam%mnnn, lpu, ierr)
       if (ierr /= 0) goto 900

    end if

    ipceq = 1
    imceq = 1
    itceq = 1
    if (sam%mpar(28) > 0 .or. sam%mpar(23) > 0) then

       !! Allocate scratch arrays since we don't know the exact size yet
       call writeProgress (' --> Creating constraint equations')
       call estimateConstrSize (sam%mpar,sam%melcon,sam%mpmnpc,nceq,nmceq)
       imceq = ipceq + 2 + nceq
       allocate(iscr(imceq+nmceq-1),rscr(nmceq),STAT=ierr)
       if (ierr /= 0) then
          ierr = allocationError('initiateSAM 3')
          return
       else if (doLogMem) then
          call logAllocMem ('initiateSAM',0,size(iscr),4)
          call logAllocMem ('initiateSAM',0,size(rscr),8)
       end if

       iscr(ipceq) = 1

    else
       allocate(iscr(0),rscr(0)) ! Only to suppress compiler warning
       nceq  = 0
       nmceq = 0
    end if

    if (nwavgm > 0) then

       !! Define the linear couplings due to interpolation constraint elements
       !! and any other explicitly defined multi-point constraints
       call MultiPointConstraints (sam%minex, sam%madof, mpwavgm, &
            &                      sam%mpmnpc, sam%mmnpc, &
            &                      tncoor(:,1), tncoor(:,2), tncoor(:,3), &
            &                      sam%mpar, sam%msc, &
            &                      iscr(ipceq:imceq-1), iscr(imceq:), &
            &                      rscr(itceq:), ipsw, lpu, ierr)
       if (ierr == -4) then
          call reportError (error_p, &
               'Insufficient scratch storage for constraint elements.', &
               'Use command-line option -bufsize_rigid to increase the size.')
          goto 900
       else if (ierr < 0) then
          call reportError (error_p,'Error associated with constraint elements')
          goto 900
       else if (ierr > 0) then
          write(errMsg,600) ierr
600       format(I8,' interpolation constraint elements (RBE3) are not fully', &
               &    ' accounted for.')
          call reportError (warning_p,adjustl(errMsg), &
               'Execution continues, but the FE model might be inconsistent.', &
               'Please inspect the details in the fedem_reducer.res file and', &
               'correct the model before using it in a dynamics simulation.')
       end if

    end if

    if (doLogMem) call logAllocMem ('initiateSAM',nwavgm,0,4)
    deallocate(mpwavgm)

    if (nrgd > 0 .or. nrbar > 0) then

       !! Define the linear couplings due to rigid elements
       ncex  = sam%nceq
       nmcex = sam%nmmceq
       ipceq = 2 + ncex
       itceq = 1 + nmcex
       imcex = imceq
       imceq = imcex + nmcex
       if (ffa_cmdlinearg_isTrue('useOldRGD3D')) then

          !! Original SAM routine
          call RGD3D (sam%minex(1), sam%madof(1), sam%melcon(1),         &
               &      sam%mpmnpc(1), sam%mmnpc(1), iscr(1), iscr(imcex), &
               &      tncoor(1,1), tncoor(1,2), tncoor(1,3), rscr(1),    &
               &      sam%mpar(1), sam%msc(1), iscr(ipceq), iscr(imceq), &
               &      rscr(itceq), ncex, krigid, 0, 0, 0, nmceq-nmcex,   &
               &      lpu, ierr)

       else

          !! This new version accounts for RGD (or RBE2) component definitions,
          !! chained RGDs, and RBARs with independent DOFs on both nodes
          call Rigid3D (mprgd, mprbar, meldep, &
               &        sam%minex, sam%madof, sam%mpmnpc, sam%mmnpc, &
               &        iscr(1:ncex+1), iscr(imcex:imcex+nmcex-1), &
               &        tncoor(:,1), tncoor(:,2), tncoor(:,3), rscr(1:nmcex), &
               &        sam%mpar, sam%msc, iscr(ipceq:imceq-1), iscr(imceq:), &
               &        rscr(itceq:), ipsw, lpu, ierr)

       end if
       if (ierr == -4) then
          call reportError (error_p, &
               'Insufficient scratch storage allocated for rigid elements.', &
               'Use command-line option -bufsize_rigid to increase the size.')
          goto 900
       else if (ierr < 0) then
          call reportError (error_p,'Error associated with rigid elements')
          goto 900
       end if

    end if

    if (doLogMem) call logAllocMem ('initiateSAM',nrbar+2*nrgd,0,4)
    deallocate(mprgd,mprbar,meldep)

    if (npbeam > 0) then

       !! Beam with pin flags exists. Modify the element topology for the
       !! affected elements and establish the associated constraint equations.
       call BeamPinConstraints (sam%madof, mpbeam, &
            &                   sam%mpmnpc, sam%mmnpc, sam%msc, sam%minex, &
            &                   sam%nceq, sam%nmmceq, iscr(ipceq:imceq-1), &
            &                   iscr(imceq:), rscr(itceq:), sam%mpar(23), &
            &                   lpu, ierr)
       if (ierr < 0) goto 900

    end if

    if (doLogMem) call logAllocMem ('initiateSAM',npbeam,0,4)
    deallocate(mpbeam)

    if (sam%nmmceq > 0) then

       !! Now allocate and store the linear coupling arrays permanently
       call reAllocate ('initiateSAM',sam%mpmceq,sam%nceq+1,ierr)
       call reAllocate ('initiateSAM',sam%mmceq ,sam%nmmceq,ierr)
       call reAllocate ('initiateSAM',sam%ttcc  ,sam%nmmceq,ierr)
       if (ierr < 0) return

       sam%mpmceq = iscr(ipceq:ipceq+sam%nceq)
       sam%mmceq  = iscr(imceq:imceq+sam%nmmceq-1)
       sam%ttcc   = rscr(itceq:itceq+sam%nmmceq-1)
       sam%nspdof = sam%nspdof + sam%nceq
       sam%nddof  = sam%nceq
       if (doLogMem) then
          call logAllocMem ('initiateSAM',size(iscr),0,4)
          call logAllocMem ('initiateSAM',size(rscr),0,8)
       end if
       deallocate(iscr,rscr)

    else

       !! Allocate length-1 arrays in the event of no linear couplings.
       !! This is needed to avoid run-time error when compiled with array bound
       !! checking, since the address of the first array element is passed to
       !! Fortran-77 routines in order to prevent stack overflow issues.
       call reAllocate ('initiateSAM',sam%mpmceq,1,ierr)
       call reAllocate ('initiateSAM',sam%mmceq ,1,ierr)
       call reAllocate ('initiateSAM',sam%ttcc  ,1,ierr)
       if (ierr < 0) return

       sam%mpmceq(1) = 1
       sam%mmceq(1)  = 0
       sam%ttcc(1)   = 0.0_dp

    end if

    sam%mpar(18) = 1 ! Flag that we now are in the reducer
    call check (sam,lpu,ierr)

900 continue
    if (ierr /= 0 .or. ipsw > 2) call writeObject (sam,lpu,2)
    if (ierr /= 0) call reportError (debugFileOnly_p,'initiateSAM')

  contains

    !> @brief Translates the external status code (2) to fixed (0).
    elemental subroutine fixExternalDOFs (isc)
      integer, intent(inout) :: isc
      if (isc == 2) isc = 0
    end subroutine fixExternalDOFs

    !> @brief Finds a conservative estimate on the constraint array sizes.
    subroutine estimateConstrSize (mpar,melcon,mpmnpc,nceq,nmceq)
      use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getint
      integer, intent(in)  :: mpar(:), melcon(:), mpmnpc(:)
      integer, intent(out) :: nceq, nmceq
      integer :: iel, nrbuf, ncex, nmcex

      nceq  = 5*mpar(23)
      nmceq = 4*nceq
      ncex  = 0
      nmcex = 0
      nrbuf = 1
      do iel = 1, mpar(2)
         if (melcon(iel) == krbar) then
            nrbuf = 2
            nceq  = nceq  +  6
            nmceq = nmceq + 21
         else if (melcon(iel) == krigid) then
            nrbuf = 2
            nceq  = nceq  +  6*(mpmnpc(iel+1)-mpmnpc(iel)-1)
            nmceq = nmceq + 21*(mpmnpc(iel+1)-mpmnpc(iel)-1)
         else if (melcon(iel) == kwavgm) then
            ncex  = ncex  + 6
            nmcex = nmcex + 6 + 24*(mpmnpc(iel+1)-mpmnpc(iel)-1)
         end if
      end do

      !! If we have both rigid/rbar and wavgm elements, the required space for
      !! the wavgm elements is doubled, because these arrays are copied once
      nceq  = nceq  + nrbuf*ncex
      nmceq = nmceq + nrbuf*nmcex

      !! For backward compability
      call ffa_cmdlinearg_getint ('bufsize_rigid',nrbuf)
      if (nrbuf > 0) nmceq = 20*mpar(23) + nrbuf*mpar(28)

    end subroutine estimateConstrSize

  end subroutine initiateSAM


  !!============================================================================
  !> @brief Writes arrays of the SAM data structure to binary file.
  !>
  !> @param[in] filename Name of binary file to store the SAM data in
  !> @param[in] chksum Checksum value for the FE part owing the @a sam object
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 29 Sep 2000

  subroutine saveSAM (filename,chksum,sam,ierr)

    use SamModule        , only : SamType
    use reportErrorModule, only : reportError, debugFileOnly_p
    use BinaryDBInterface, only : openBinaryDB, closeBinaryDB, writeTagDB
    use BinaryDBInterface, only : writeIntDB, writeDoubleDB, write_p

    character*(*), intent(in)  :: filename
    integer      , intent(in)  :: chksum
    type(SamType), intent(in)  :: sam
    integer      , intent(out) :: ierr

    !! Local variables
    integer :: ifile

    !! --- Logic section ---

    call openBinaryDB (filename,write_p,ifile,ierr)
    if (ierr < 0) goto 990

    call writeTagDB (ifile,'#SAM data',chksum,ierr)
    if (ierr < 0) goto 990

    call writeIntDB (ifile,size(sam%mpar),1,ierr)
    if (ierr < 0) goto 990

    call writeIntDB (ifile,sam%mpar(1),size(sam%mpar),ierr)
    if (ierr < 0) goto 990

    call writeIntDB (ifile,sam%madof(1),size(sam%madof),ierr)
    if (ierr < 0) goto 990

    call writeIntDB (ifile,sam%minex(1),size(sam%minex),ierr)
    if (ierr < 0) goto 990

    call writeIntDB (ifile,sam%mnnn(1),size(sam%mnnn),ierr)
    if (ierr < 0) goto 990

    call writeIntDB (ifile,sam%msc(1),size(sam%msc),ierr)
    if (ierr < 0) goto 990

    call writeIntDB (ifile,sam%mpmnpc(1),size(sam%mpmnpc),ierr)
    if (ierr < 0) goto 990

    call writeIntDB (ifile,sam%mmnpc(1),sam%nmmnpc,ierr)
    if (ierr < 0) goto 990

    call writeIntDB (ifile,sam%melcon(1),size(sam%melcon),ierr)
    if (ierr < 0) goto 990

    if (sam%nceq > 0) then
       call writeIntDB (ifile,sam%mpmceq(1),size(sam%mpmceq),ierr)
       if (ierr < 0) goto 990

       call writeIntDB (ifile,sam%mmceq(1),size(sam%mmceq),ierr)
       if (ierr < 0) goto 990

       call writeDoubleDB (ifile,sam%ttcc(1),size(sam%ttcc),ierr)
       if (ierr < 0) goto 990
    end if

    if (associated(sam%meqn)) then
       call writeIntDB (ifile,sam%meqn(1),size(sam%meqn),ierr)
       if (ierr < 0) goto 990
    end if

    if (associated(sam%meqn1)) then
       call writeIntDB (ifile,sam%meqn1(1),size(sam%meqn1),ierr)
       if (ierr < 0) goto 990
    end if

    if (associated(sam%meqn2)) then
       call writeIntDB (ifile,sam%meqn2(1),size(sam%meqn2),ierr)
       if (ierr < 0) goto 990
    end if

    if (associated(sam%dofPosIn2)) then
       call writeIntDB (ifile,sam%dofPosIn2(1),size(sam%dofPosIn2),ierr)
       if (ierr < 0) goto 990
    end if

    call closeBinaryDB (ifile,ierr)
    if (ierr == 0) return

990 call reportError (debugFileOnly_p,'saveSAM')

  end subroutine saveSAM

end module samReducerModule
