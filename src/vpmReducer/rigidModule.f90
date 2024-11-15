!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file rigidModule.f90
!>
!> @brief Linear multi-point constraint handling.

!!==============================================================================
!> @brief Module with subroutines for setting up linear multi-point constraints.
!>
!> @details This module contains a collection of subroutines for processing the
!> rigid coupling elements (RGD, RBAR, etc.) for generation of the equivalent
!> constraint equations (CE) of the algebraic system.

module RigidModule

  implicit none

  private

  integer, parameter, public :: krigid = 61 !< Type id for rigid elements
  integer, parameter, public :: krbar  = 62 !< Type id for rigid bars
  integer, parameter, public :: kwavgm = 63 !< Type id for constraint elements

  public :: MultiPointConstraints, checkRgdElements, checkRbarElements, Rigid3D


contains

  !> @cond NO_DOCUMENTATION
  subroutine reportRigidError (numErr, msg1, msg2, msg3)

    use reportErrorModule, only : reportError, error_p, errorFileOnly_p

    character(len=*), intent(in)           :: msg1
    character(len=*), intent(in), optional :: msg2, msg3
    integer         , intent(inout)        :: numErr

    if (numErr < 10) then
       call reportError (error_p,msg1,msg2,msg3)
    else
       call reportError (errorFileOnly_p,msg1,msg2,msg3)
    end if
    numErr = numErr + 1

  end subroutine reportRigidError
  !> @endcond

  !!============================================================================
  !> @brief Generates constraint equation arrays in the SAM datastructure.
  !>
  !> @param[in] minex Matrix of internal to external node number
  !> @param[in] madof Matrix of accumulated DOFs
  !> @param[in] mpam Matrix of WAVGM element indices
  !> @param[in] mpmnpc Matrix of pointers to MNPCs
  !> @param[in] mmnpc Matrix of nodal point correspondances
  !> @param[in] txc Table of x-coordinates for the nodal points
  !> @param[in] tyc Table of Y-coordinates for the nodal points
  !> @param[in] tzc Table of Z-coordinates for the nodal points
  !> @param mpar Matrix of parameters
  !> @param msc Matrix of status codes
  !> @param[out] mpmceq Matrix of pointers to MCEQs
  !> @param[out] mmceq Matrix of constraint equation definitions
  !> @param[out] ttcc Table of constraint equation coefficients
  !> @param[in] ipsw Print switch
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @details This subroutine generates the constraint information arrays,
  !> MPMCEQ, MMCEQ and TTCC, that account for implicitly defined constraints
  !> via the interpolation constraint elements (WAVGM or RBE3) and any other
  !> explicitly defined linear multi-point constraints.
  !>
  !> Status codes in MSC are updated and the final number of constraint
  !> equations (NCEQ) is stored in MPAR(7). The actual number of elements in
  !> MMCEQ and TTCC containing relevant information (NMMCEQ) is also
  !> determined and returned in MPAR(16).
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 15 Aug 2002

  subroutine MultiPointConstraints (MINEX,MADOF,MPAM,MPMNPC,MMNPC,TXC,TYC,TZC, &
       &                            MPAR,MSC,MPMCEQ,MMCEQ,TTCC,IPSW,LPU,IERR)

    use kindModule             , only : dp, epsDiv0_p
    use WAvgMotionModule       , only : wavgmConstrEqn
    use allocationModule       , only : doLogMem, logAllocMem
    use reportErrorModule      , only : allocationError, reportError
    use reportErrorModule      , only : warningFileOnly_p, debugFileOnly_p
    use FFlLinkHandlerInterface, only : ffl_getElmId, ffl_getwavgm
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_isTrue
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_getdouble

    integer , intent(in)    :: minex(:), madof(:), mpam(:), mpmnpc(:), mmnpc(:)
    integer , intent(in)    :: ipsw, lpu
    real(dp), intent(in)    :: txc(:), tyc(:), tzc(:)
    integer , intent(inout) :: mpar(:), msc(:)
    integer , intent(out)   :: mpmceq(:), mmceq(:)
    real(dp), intent(out)   :: ttcc(:)
    integer , intent(out)   :: ierr

    !! Local variables
    character(len=128)    :: errMsg
    integer               :: i, j, k, l, n, iel, xel, iam, iceq, ipcur, iwarn
    integer               :: depDof, refC, ignoredC(3), indC(6), nenod, nW
    logical               :: refDofs(6)
    real(dp)              :: sumW, epsX, X0(3), tolX(3)
    integer , allocatable :: mDofs(:)
    real(dp), allocatable :: deltaX(:,:), work(:), weight(:), omega(:)

    !! --- Logic section ---

    call ffa_cmdlinearg_getdouble ('tolWAVGM',epsX)

    !! Find the largest WAVGM element in terms of number of independent nodes
    nenod = 0
    do iam = 1, size(mpam)
       iel = mpam(iam)
       nenod = max(nenod,mpmnpc(iel+1)-mpmnpc(iel)-1)
    end do
    call ffl_getwavgm (nW,indC,X0,-1,ierr)
    if (ierr < 0) goto 900

    !! Allocate local work arrays
    allocate(mDofs(6*nenod), STAT=ierr)
    if (ierr /= 0) then
       ierr = AllocationError('MultiPointConstraints')
       return
    else if (doLogMem) then
       call logAllocMem ('MultiPointConstraints',0,6*nenod,4)
    end if
    allocate(deltaX(3,1+nenod),work(nenod), &
         &   weight(nW),omega(6*nenod), STAT=ierr)
    if (ierr /= 0) then
       ierr = AllocationError('MultiPointConstraints')
       return
    else if (doLogMem) then
       call logAllocMem ('MultiPointConstraints',0,3+10*nenod+nW,8)
    end if

    !! Process all weighted averaged motion elements (RBE3)
    iceq  = 1
    ipcur = 1
    iwarn = 0
    do iam = 1, size(mpam)

       iel = mpam(iam)
       xel = ffl_getElmId(iel)
       nenod = mpmnpc(iel+1) - mpmnpc(iel) - 1

       if ( ffa_cmdlinearg_isTrue('useOldWAVGM') .or. &     !TODO: Retained for
            ffa_cmdlinearg_isTrue('useUnWeightedCG') ) then !TODO: testing only

          !! Compute element centroid and offset from each node to the centroid
          do i = 1, nenod
             n = mmnpc(mpmnpc(iel)+i)
             X0(1) = X0(1) + txc(n)
             X0(2) = X0(2) + tyc(n)
             X0(3) = X0(3) + tzc(n)
          end do
          X0 = X0 / nenod
          do i = 1, 1+nenod
             n = mmnpc(mpmnpc(iel)+i-1)
             deltaX(1,i) = txc(n) - X0(1)
             deltaX(2,i) = tyc(n) - X0(2)
             deltaX(3,i) = tzc(n) - X0(3)
          end do

          !! Compute geometric tolerances
          do j = 1, 3
             tolX(j) = max(maxval(-deltaX(j,1:1+nenod))*epsX, &
                  &        maxval( deltaX(j,1:1+nenod))*epsX, epsDiv0_p)
          end do
          if (ipsw > 3) then
             write(lpu,600) xel
             write(lpu,610) tolX, X0, deltaX(:,1), (i,deltaX(:,1+i),i=1,nenod)
          end if

       else if (ipsw > 3) then !TODO: This will later be replaced by an "if"

          write(lpu,600) xel
          write(lpu,620) (minex(mmnpc(j)),j=mpmnpc(iel),mpmnpc(iel)+nenod)

       end if

       !! Get element properties
       call ffl_getwavgm (refC,indC,weight,iel,ierr)
       if (ierr == 0) then
          !! No properties defined - use default settings
          refC = 123456
          indC(1:3) = -1
          indC(4:6) = 0
       else if (ierr < 0) then
          goto 900
       end if

       if (ffa_cmdlinearg_isTrue('useOldWAVGM')) then !TODO: Remove this later
          !! Normalize the weighting factors
          j = indC(1)
          if (j > 0) then
             k = j + nenod - 1
             sumW = sum(weight(j:k))
             if (abs(sumW) > epsDiv0_p) weight(j:k) = weight(j:k) / sumW
          end if
          do i = 2, 6
             j = indC(i)
             if (j < 1 .or. j == indC(i-1)) cycle
             k = j + nenod - 1
             sumW = sum(weight(j:k))
             if (abs(sumW) > epsDiv0_p) weight(j:k) = weight(j:k) / sumW
          end do
          call procDofComp (refC,refDofs,iel,lpu,ierr)
       else
          call procDofComp (abs(refC),refDofs,iel,lpu,ierr)
       end if
       if (ierr < 0) goto 900

       ignoredC = 0
       do i = 1, 6
          if (.not. refDofs(i)) cycle

          !! Find global DOF number for current dependent DOF
          depDof = madof(mmnpc(mpmnpc(iel))) + i-1
          if (depDof >= madof(mmnpc(mpmnpc(iel))+1)) then
             ignoredC(1) = 10*ignoredC(1) + i ! Ignore invalid dependent DOF
             cycle
          else if (msc(depDof) /= 1) then
             ignoredC(2) = 10*ignoredC(2) + i ! This DOF is aready constrained
             cycle
          end if

          if (refC >= 0) then
             !! Compute constraint coefficients from the WAVGM definitions
             call wavgmConstrEqn (xel,i,nenod,indC,txc,tyc,tzc, &
                  &               mmnpc(mpmnpc(iel):mpmnpc(iel+1)-1), &
                  &               epsX,tolX,deltaX,work,weight,omega,ipsw,lpu)
          else ! Explicit constraints (assuming 6 DOFs per independent node)
             call explConstrEqn (6*nenod,indC(i),weight,omega)
          end if

          !! Find global DOF numbers for the independent DOFs
          mDofs = 0
          do j = 1, nenod
             n = mmnpc(mpmnpc(iel)+j)
             do k = 1, madof(n+1) - madof(n)
                l = 6*j-6 + k
                if (abs(omega(l)) > epsDiv0_p) mDofs(l) = madof(n) + k-1
             end do
          end do
          if (ipsw > 3) then
             write(lpu,630) depDof, mDofs(1:6)
             write(lpu,631) omega(1:6)
             do l = 7, 6*nenod, 6
                write(lpu,632) mDofs(l:l+5)
                write(lpu,633) omega(l:l+5)
             end do
          end if

          if (maxval(mDofs) > 0) then

             !! Insert constraint equation into the SAM data structure
             call defineConstraintEqn (depDof,mDofs(1:6*nenod),omega, &
                  &                    mpmceq(1:1),mmceq,ttcc, & ! dummy args
                  &                    msc,mpmceq,mmceq,ttcc, &
                  &                    iceq,ipcur,lpu,ierr)
             if (ierr < 0) goto 900

          else
             ignoredC(3) = 10*ignoredC(3) + i ! Dependent DOF has no coupling
          end if

       end do

       if (ignoredC(1) > 0) then
          write(errMsg,690) ignoredC(1),xel
          call reportError (warningFileOnly_p,errMsg, &
               'No other structural elements are connected to these DOFs.')
       end if
       if (ignoredC(2) > 0) then
          write(errMsg,690) ignoredC(2),xel
          call reportError (warningFileOnly_p,errMsg, &
               'These DOFs are either external (or fixed), '// &
               'or dependent in previously defined constraint elements.')
       end if
       if (ignoredC(3) > 0) then
          write(errMsg,690) ignoredC(3),xel
          call reportError (warningFileOnly_p,errMsg, &
               'No associated independent DOFs.')
       end if
       if (maxval(ignoredC) > 0) iwarn = iwarn + 1
    end do
    ierr = iwarn

    if (doLogMem) then
       call logAllocMem ('MultiPointConstraints',size(mDofs),0,4)
       call logAllocMem ('MultiPointConstraints', &
            &            size(deltaX)+size(work)+size(weight)+size(omega),0,8)
    end if
    deallocate(mDofs,deltaX,work,weight,omega)

    !! Restore MSC, set NCEQ and NMMCEQ and complete MPMCEQ
    do i = 1, mpar(3)
       if (msc(i) < 0) msc(i) = 0
    end do
    mpar(7)  = iceq-1
    mpar(16) = ipcur-1
    mpmceq(iceq) = ipcur

    return

600 format(/5X,'>>> Processing WAVGM element',I8,' <<<')
610 format(/19X,'X-coor',7X,'Y-coor',7X,'Z-coor', &
         & /5X,'Tolerance     ',1P3E13.5, &
         & /5X,'Centroid      ',1P3E13.5, &
         & /5X,'Dependent node',1P3E13.5, &
         &/(5X,'Node      ',I4, 1P3E13.5))
620 format(/9X,'Ref node',I8 / 9X,'Nodes   ',10I8 /(17X,10I8))
630 format(/5X,'Dependent DOF   :',I8 / 5X,'Independent DOFs:',6I8)
631 format( 5X,'Weights         :',6F8.3)
632 format(22X,6I8)
633 format(22X,6F8.3)
690 format('Dependent DOF(s)',I7,' of constraint element',I8,' are ignored.')

900 continue
    call reportError (debugFileOnly_p,'MultiPointConstraints')

contains

  !> @brief Copies coefficients of an explicitly defined constraint equation.
  !> @details Assuming here that all independent nodes have six DOFs.
  !> See also FFlSesamReader::readLinearDependencies(), where the
  !> weight matrix is populated from the SESAM BLDEP records.
  subroutine explConstrEqn (nedof,indC,weight,omega)
    integer , intent(in)  :: nedof, indC
    real(dp), intent(in)  :: weight(:)
    real(dp), intent(out) :: omega(:)
    if (indC > 0 .and. indC+nedof-1 <= size(weight)) then
       call DCOPY (nedof,weight(indC),1,omega(1),1)
    else ! Assume equal weight on all DOFs
       call DCOPY (nedof,1.0_dp,0,omega(1),1)
    end if
  end subroutine explConstrEqn

  end subroutine MultiPointConstraints


  !!============================================================================
  !> @brief Checks the topology of the RGD elements.
  !>
  !> @param[in] mprgd Matrix of RGD element indices
  !> @param[in] mprbar Matrix of RBAR element indices
  !> @param[in] mpam Matrix of WAVGM element indices
  !> @param[out] meldep RGD element chain definition (see below)
  !> @param[in] minex Matrix of internal to external node number
  !> @param[in] mpmnpc Matrix of pointers to MNPCs
  !> @param mmnpc Matrix of nodal point correspondances
  !> @param nodeStatus Nodal status flags
  !> @param[out] numErr Error flag
  !>
  !> @details The RGD elements are checked against the following rules:
  !>  - Only one node per RGD element may be on a triad (external DOFs).
  !>  - If a dependent node (and only one) is on a triad, the independent node
  !>    of that element cannot be a dependent in another constraint element.
  !>  - Common dependent nodes is allowed in two RGD elements only if they also
  !>    share a common independent node.
  !>  - An independent node of a RGD-element may be the dependent node in only
  !>    one other element.
  !>
  !> If a dependent node (and only one) is on a triad, and the second condition
  !> above is satisfied, that dependent node and the independent node of the
  !> element are swapped.
  !>
  !> This subroutine also establishes the topology of chained RGD elements.
  !> The output array @a meldep is a mapping @a ielA-->ielB, where RGD element
  !> @a ielB is the element in which the reference node of RGD element @a ielA
  !> is a dependent node. If the reference node of element @a ielA is
  !> not dependent, then @a ielB=0.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 21 Aug 2003

  subroutine checkRgdElements (mprgd, mprbar, mpam, meldep, minex, &
       &                       mpmnpc, mmnpc, nodeStatus, numErr)

    use reportErrorModule      , only : reportError, note_p, error_p
    use FFlLinkHandlerInterface, only : ffl_getElmId

    integer, intent(in)    :: mprgd(:), mprbar(:), mpam(:)
    integer, intent(in)    :: minex(:), mpmnpc(:)
    integer, intent(inout) :: mmnpc(:), nodeStatus(:)
    integer, intent(out)   :: meldep(:), numErr

    !! Local variables
    integer :: i1, i2, iel1, iel2, inod1, inod2, jnod1, jnod2, nrgd
    character(len=128) :: errMsg1, errMsg2

    !! --- Logic section ---

    nrgd   = size(mprgd)
    numErr = 0
    meldep = 0

    !! Increment the nodeStatus by 10 for all dependent nodes
    do iel1 = 1, nrgd
       do i1 = mpmnpc(mprgd(iel1))+1, mpmnpc(mprgd(iel1)+1)-1
          inod1 = mmnpc(i1)
          nodeStatus(inod1) = nodeStatus(inod1) + 10
       end do
    end do
    do iel1 = 1, size(mprbar)
       do i1 = mpmnpc(mprbar(iel1)), mpmnpc(mprbar(iel1)+1)-1
          inod1 = mmnpc(i1)
          nodeStatus(inod1) = nodeStatus(inod1) + 10
       end do
    end do
    do iel1 = 1, size(mpam)
       inod1 = mmnpc(mpmnpc(mpam(iel1)))
       nodeStatus(inod1) = nodeStatus(inod1) + 10
    end do

    !! Check that not more than one RGD-node has triad status (external node).
    !! If one dependent node has triad status, try to swap it with the
    !! independent node of that element.
    do iel1 = 1, nrgd
       i2 = 0
       inod1 = 0
       do i1 = mpmnpc(mprgd(iel1)), mpmnpc(mprgd(iel1)+1)-1
          if (mod(nodeStatus(mmnpc(i1)),10) == 2) then
             if (i2 == 0) then
                i2 = i1-mpmnpc(mprgd(iel1))+1 ! The first external node
                inod1 = mmnpc(i1)
             else if (i2 > 0) then
                i2 = -2 ! Second external node
             else
                i2 = i2-1 ! More than two external nodes
             end if
          end if
       end do
       if (i2 < 0) then
          !! More than one external node in the element
          write(errMsg1,600) ffl_getElmId(mprgd(iel1)), -i2
          call reportRigidError(numErr,errMsg1)
       else if (i2 > 1) then
          !! A dependent node has triad status,
          !! swap it with the independent node, if possible
          jnod1 = mmnpc(mpmnpc(mprgd(iel1)))
          if (nodeStatus(jnod1) < 10) then
             !! Swap inod1 with jnod1 in this elements connectivity table
             mmnpc(mpmnpc(mprgd(iel1)))      = inod1
             mmnpc(mpmnpc(mprgd(iel1))+i2-1) = jnod1
             nodeStatus(inod1) = nodeStatus(inod1) - 10
             nodeStatus(jnod1) = nodeStatus(jnod1) + 10
             write(errMsg1,610) ffl_getElmId(mprgd(iel1)), minex(inod1)
             write(errMsg2,611) minex(jnod1)
             call reportError (note_p,errMsg1,errMsg2)
          else
             !! Unable to swap, because the independent node is a dependent
             !! in another constraint element
             write(errMsg1,610) ffl_getElmId(mprgd(iel1)), minex(inod1)
             write(errMsg2,612) minex(jnod1)
             call reportRigidError(numErr,errMsg1,errMsg2, &
                  'because the latter node is a dependent in another element')
          end if
       end if
    end do

    !! Reset the nodeStatus array
    do inod1 = 1, size(nodeStatus)
       if (nodeStatus(inod1) > 9) nodeStatus(inod1) = mod(nodeStatus(inod1),10)
    end do

    !! Check for common dependent nodes
    do iel1 = 1, nrgd
       jnod1 = mmnpc(mpmnpc(mprgd(iel1)))
       do i1 = mpmnpc(mprgd(iel1))+1, mpmnpc(mprgd(iel1)+1)-1
          inod1 = mmnpc(i1)

          do iel2 = iel1+1, nrgd
             jnod2 = mmnpc(mpmnpc(mprgd(iel2)))
             do i2 = mpmnpc(mprgd(iel2))+1, mpmnpc(mprgd(iel2)+1)-1
                inod2 = mmnpc(i2)
                !! Error if common dependent node and not same independent node
                if (inod2 == inod1 .and. jnod1 /= jnod2) then
                   write(errMsg1,620) ffl_getElmId(mprgd(iel1)), &
                        &             ffl_getElmId(mprgd(iel2)),minex(inod1)
                   call reportRigidError(numErr,errMsg1)
                end if
             end do
          end do

          !! Check connectivity
          do iel2 = 1, nrgd
             jnod2 = mmnpc(mpmnpc(mprgd(iel2)))
             if (jnod2 == inod1) meldep(iel2) = iel1
          end do

       end do
    end do

    if (numErr > 0) then
       write(errMsg1,690) numErr
       call reportError (error_p,errMsg1,addString='checkRgdElements')
    end if

600 format('RGD element',i8,' has a triad (external DOFs) in',i3,' nodes')
610 format('RGD element',i8,' has a triad (external DOFs) on dependent node',i8)
611 format('This dependent node is swapped with the independent node',i8)
612 format('Unable to swap this dependent node with independent node',i8)
620 format('RGD elements',i8,' and',i8, ' have common dependent node',i8)
690 format('Found',i6,' errors associated with RGD elements')

  end subroutine checkRgdElements


  !!============================================================================
  !> @brief Checks consistency of RBAR elements.
  !>
  !> @param[in] mprbar Matrix of RBAR element indices
  !> @param[in] mprgd Matrix of RGD element indices
  !> @param[in] minex Matrix of internal to external node number
  !> @param[in] mpmnpc Matrix of pointers to MNPCs
  !> @param[in] mmnpc Matrix of nodal point correspondances
  !> @param[in] madof Matrix of accumulated DOFs
  !> @param[in] msc Matrix of status codes
  !> @param[in] nodeStatus Nodal status flags
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] numErr Error flag
  !>
  !> @details This subroutine checks that RBAR and RGD elements have
  !> no common nodes, and that no RBAR nodes are external.
  !>
  !> The subroutine also checks DOF components of RBAR elements to ensure that:
  !> - No DOF is defined to be both dependent and independent.
  !> - Independent DOFs jointly define any rigid body motion.
  !> - No rotational independent DOFs exist on 3-DOF nodes.
  !>
  !> If rotational dependent DOFs exist on 3-DOF nodes, a warning is issued.
  !>
  !> @callergraph
  !>
  !> @author Tommy Jorstad
  !>
  !> @date 29 Sep 2002

  subroutine checkRbarElements (mprbar, mprgd, minex, mpmnpc, mmnpc, &
       &                        madof, msc, nodeStatus, lpu, numErr)

    use reportErrorModule      , only : reportError
    use reportErrorModule      , only : error_p, warning_p, warningFileOnly_p
    use FFlLinkHandlerInterface, only : ffl_getElmId

    integer, intent(in)  :: mprgd(:), mprbar(:), minex(:), mpmnpc(:), mmnpc(:)
    integer, intent(in)  :: madof(:), msc(:), nodeStatus(:), lpu
    integer, intent(out) :: numErr

    !! Local variables
    integer            :: i, j, k, j1, j2, e1, e2
    integer            :: rgdComp(4), nrgd, nrbar, numWarns, ierr
    integer            :: nodeA, nodeB, nodeC, nodeD, nndofA, nndofB, fdA, fdB
    logical            :: md1(6), md2(6), sd1(6), sd2(6), dofsOk
    character(len=128) :: errMsg

    !! --- Logic section ---

    nrgd  = size(mprgd)
    nrbar = size(mprbar)
    numErr = 0
    numWarns = 0

    do j1 = 1, nrbar
       e1 = mprbar(j1)
       nodeA = mmnpc(mpmnpc(e1))
       nodeB = mmnpc(mpmnpc(e1)+1)

       !! Check for common RGD/RBAR nodes
       do j2 = 1, nrgd
          e2 = mprgd(j2)
          do i = mpmnpc(e2), mpmnpc(e2+1)-1
             nodeC = mmnpc(i)
             if (nodeA == nodeC .or. nodeB == nodeC) then
                write(errMsg,600) ffl_getElmId(e2),minex(nodeC),ffl_getElmId(e1)
                call reportRigidError(numErr,errMsg)
             end if
          end do
       end do

       !! Check for common RBAR/RBAR nodes
       do j2 = j1+1, nrbar
          e2 = mprbar(j2)
          nodeC = mmnpc(mpmnpc(e2))
          nodeD = mmnpc(mpmnpc(e2)+1)
          if (nodeA == nodeC .or. nodeA == nodeD .or. &
              nodeB == nodeC .or. nodeB == nodeD) then
             write(errMsg,640) ffl_getElmId(e2),ffl_getElmId(e1)
             call reportRigidError(numErr,errMsg)
          end if
       end do

       !! Check that the RBAR nodes do not have triad status
       if (nodeStatus(nodeA) == 2 .or. nodeStatus(nodeB) == 2) then
          if (nodeStatus(nodeA) == 2) then
             write(errMsg,610) ffl_getElmId(e1),minex(nodeA)
          else if (nodeStatus(nodeB) == 2) then
             write(errMsg,610) ffl_getElmId(e1),minex(nodeB)
          end if
          call reportRigidError(numErr,errMsg)
       end if

       !! Set up RBAR component definition
       call ffl_getRgdDofComp (rgdComp,e1,ierr)
       dofsOk = (ierr == 0)

       call procRBarDef (rgdComp,md1,md2,sd1,sd2,e1,lpu,ierr)
       if (ierr /= 0) dofsOk = .false.

       !! There must be 6 independent DOFs present
       if (count(md1) + count(md2) /= 6) dofsOk = .false.

       !! No DOF may be both independent and dependent on the same node
       do i = 1, 6
          if (md1(i) .and. sd1(i)) dofsOk = .false.
          if (md2(i) .and. sd2(i)) dofsOk = .false.
       end do

       !! The independent DOFs must represent any motion
       do i = 1, 3
          if (.not. (md1(i) .or. md2(i))) dofsOk = .false.
       end do

       !! If no independent x-rotation  ->  y- or z-translation must be fixed
       !! If no independent y-rotation  ->  z- or x-translation must be fixed
       !! If no independent z-rotation  ->  x- or y-translation must be fixed
       do i = 4, 6
          if (.not. (md1(i) .or. md2(i))) then
             j = 1 + mod(i,3)
             if (.not. (md1(j) .and. md2(j))) then
                k = 1 + mod(j,3)
                if (.not. (md1(k) .and. md2(k))) dofsOk = .false.
             end if
          end if
       end do

       if (.not. dofsOk) then
          write(errMsg,620) ffl_getElmId(e1), rgdComp
          call reportRigidError(numErr,errMsg)
       end if

       !! Check for rotational DOFs on 3-DOF nodes
       fdA = madof(nodeA)
       nndofA = madof(nodeA+1) - fdA
       if (nndofA == 6) then
          if (msc(fdA+3)+msc(fdA+4)+msc(fdA+5) == 0) nndofA = 3
       else
          nndofA = 3
       end if
       fdB = madof(nodeB)
       nndofB = madof(nodeB+1) - fdB
       if (nndofB == 6) then
          if (msc(fdB+3)+msc(fdB+4)+msc(fdB+5) == 0) nndofB = 3
       else
          nndofB = 3
       end if

       do i = 4, 6
          if (nndofA == 3 .and. (md1(i) .or. md2(i))) then
             ierr = -1
             write(errMsg,630) 'independent', i, ffl_getElmId(e1)
             call reportRigidError(numErr,errMsg)
          end if
          if (nndofB == 3 .and. (sd2(i) .or. sd1(i))) then
             write(errMsg,630) 'dependent', i, ffl_getElmId(e1)
             if (numWarns < 10) then
                call reportError (warning_p,errMsg)
             else
                call reportError (warningFileOnly_p,errMsg)
             end if
             numWarns = numWarns + 1
          end if
       end do
    end do

    if (numErr > 0) then
       write(errMsg,690) numErr, numWarns
       call reportError (error_p,errMsg,addString='checkRbarElements')
    end if

600 format('RGD element',i8,' has common node',i8,' with RBAR element',i8)
610 format('RBAR element',i8, &
         & ' has a triad (external DOFs) on dependent node',i8)
620 format('Invalid DOF components in RBAR element',i8,' :',4i8)
630 format('Invalid ',a,' DOF',i2,' on 3-DOF node in RBAR element',i8)
640 format('RBAR elements',i8,' and',i8,' have common nodes')
690 format('Found',i6,' errors and',i6, &
         & ' warnings associated with RBAR elements')

  end subroutine checkRbarElements


  !!============================================================================
  !> @brief Generates constraint equation arrays in the SAM datastructure.
  !>
  !> @param[in] mprgd Matrix of RGD element indices
  !> @param[in] mprbar Matrix of RBAR element indices
  !> @param meldep RGD element couplings
  !> @param[in] minex Matrix of internal to external node number
  !> @param[in] madof Matrix of accumulated DOFs
  !> @param[in] mpmnpc Matrix of pointers to MNPCs
  !> @param[in] mmnpc Matrix of nodal point correspondances
  !> @param[in] mpmcex Matrix of pointers to external MCEQs
  !> @param[in] mmcex Matrix of external constraint equation definitions
  !> @param[in] txc Table of x-coordinates for the nodal points
  !> @param[in] tyc Table of Y-coordinates for the nodal points
  !> @param[in] tzc Table of Z-coordinates for the nodal points
  !> @param[in] ttccx Table of external constraint equation coefficients
  !> @param mpar Matrix of parameters
  !> @param msc Matrix of status codes
  !> @param[out] mpmceq Matrix of pointers to MCEQs
  !> @param[out] mmceq Matrix of constraint equation definitions
  !> @param[out] ttcc Table of constraint equation coefficients
  !> @param[in] ipsw Print switch
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @details This subroutine generates new versions of the constraint
  !> information arrays MPMCEQ, MMCEQ and TTCC, that account for both implicitly
  !> defined constraints via rigid and semi-rigid elements (RGD and RBAR) and
  !> explicitly defined constraints recorded in arrays MPMCEX, MMCEX and TTCCX.
  !>
  !> The RGD element is a completely rigid or semi-rigid 3D element with 6 DOFs,
  !> an arbitrary number of nodes and all independent DOFs associated with the
  !> first (reference) node. The DOFs of the remaining element nodes may then
  !> be coupled to the independent node, making them dependent DOFs.
  !> However, not all DOFs in these nodes need not be dependent.
  !>
  !> The RBAR element is a completely rigid or semi-rigid two-noded element with
  !> independent DOFs in both nodes, and an arbitrary number of dependent DOFs.
  !>
  !> Status codes in MSC are updated and the final number of constraint
  !> equations (NCEQ) is stored in MPAR(7). The actual number of elements in
  !> MMCEQ and TTCC containing relevant information (NMMCEQ) is also
  !> determined and returned in MPAR(16).
  !>
  !> This subroutine is an extended version of the SAM routine RGD3D.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 18 Jul 2002

  subroutine Rigid3D (mprgd,mprbar,meldep,minex,madof, &
       &              mpmnpc,mmnpc,mpmcex,mmcex,txc,tyc,tzc,ttccx, &
       &              mpar,msc,mpmceq,mmceq,ttcc,ipsw,lpu,ierr)

    use kindModule             , only : dp
    use scratchArrayModule     , only : getIntegerScratchArray
    use reportErrorModule      , only : reportError, debugFileOnly_p
    use FFlLinkHandlerInterface, only : ffl_getElmId

    integer , intent(in)    :: mprgd(:), mprbar(:), minex(:), mpmnpc(:),mmnpc(:)
    integer , intent(in)    :: madof(:), mpmcex(:), mmcex(:)
    real(dp), intent(in)    :: txc(:), tyc(:), tzc(:), ttccx(:)
    integer , intent(inout) :: meldep(:), mpar(:), msc(:)
    integer , intent(out)   :: mpmceq(:), mmceq(:)
    real(dp), intent(out)   :: ttcc(:)
    integer , intent(in)    :: ipsw, lpu
    integer , intent(out)   :: ierr

    !! Local variables
    logical          :: slvDofs1(6), mstDofs1(6), slvDofs2(6), mstDofs2(6)
    integer          :: i, iceq, ipcur, idof, iel, jp, mstNod
    integer          :: nrgd, nrbar, nm, ncex, maxsiz, rgdComp(4)
    integer, pointer :: elDepth(:)

    !! --- Logic section ---

    ierr  = 0
    nrgd  = size(mprgd)
    nrbar = size(mprbar)
    ncex  = size(mpmcex)-1

    !! Set status codes (in MSC) for explicitly constrained DOFs to -1
    !! for prescribed DOFs, and to -2 for dependent DOFs
    do i = 1, ncex
       jp = mpmcex(i)
       idof = mmcex(jp)
       if (idof < 1 .or. idof > mpar(3)) then
          call RGDERR (1,i,idof,idof,lpu,ierr)
          goto 900
       end if
       nm = mpmcex(i+1) - jp - 1
       if (nm == 0) then
          msc(idof) = -1
       else if (nm > 0) then
          msc(idof) = -2
       else
          call RGDERR (7,i,nm,nm,lpu,ierr)
          goto 900
       end if
    end do

    elDepth => getIntegerScratchArray(nrgd,ierr)
    if (ierr < 0) goto 900

    !! For each RGD element check status of each independent and dependent DOF,
    !! and increase legal positive status code by 10
    nm = 0
    do i = 1, nrgd
       iel = mprgd(i)
       mstDofs1 = .true. ! All DOFs are free for an RBE2 reference node
       slvDofs1 = .false.

       !! Check independent DOFs
       mstNod = mmnpc(mpmnpc(iel))
       !! Note that krigid is used to specify the RGD kind
       call chkAndIncIndepDofs (mstNod,krigid,mstDofs1,iel,ierr)
       if (ierr < 0) goto 900

       !! Get dependent DOF components
       call ffl_getRgdDofComp (rgdComp,iel,ierr)
       if (ierr < 0) goto 900

       call procDofComp (rgdComp(1),slvDofs1,iel,lpu,ierr)
       if (ierr < 0) goto 900

       !! Check the dependent DOFs
       do jp = mpmnpc(iel)+1, mpmnpc(iel+1)-1
          call chkAndIncDepDofs (mmnpc(jp),slvDofs1,iel,ierr)
          if (ierr < 0) goto 900
       end do

       !! Calculate the chain depth for this RGD element
       elDepth(i) = calcDepth(i)
       if (elDepth(i) < 0) then
          write(lpu,*)
          nm = nm + 1
       else if (elDepth(i) > 1 .and. ipsw > 0) then
          write(lpu,600,ADVANCE='no') ffl_getElmId(iel)
          jp = meldep(i)
          do while (jp > 0)
             write(lpu,"(I8)",ADVANCE='no') ffl_getElmId(mprgd(jp))
             jp = meldep(jp)
          end do
          write(lpu,*)
       end if

    end do
    ierr = -nm
    if (ierr < 0) goto 900

    !! For each RBAR element check status of each independent and dependent DOF,
    !! and increase legal positive status code by 10
    do i = 1, nrbar
       iel = mprbar(i)

       !! Get component definition
       call ffl_getRgdDofComp (rgdComp,iel,ierr)
       if (ierr < 0) goto 900

       call procRBarDef (rgdComp,mstDofs1,mstDofs2,slvDofs1,slvDofs2, &
            &            iel,lpu,ierr)
       if (ierr < 0) goto 900

       !! Check independent DOFs (RBAR may have independent DOFs on both nodes)
       if (any(mstDofs1)) then
          jp   = mpmnpc(iel)
          mstNod = mmnpc(jp)
          call chkAndIncIndepDofs (mstNod,krbar,mstDofs1,iel,ierr)
          if (ierr < 0) goto 900
       end if
       if (any(mstDofs2)) then
          jp   = mpmnpc(iel) + 1
          mstNod = mmnpc(jp)
          call chkAndIncIndepDofs (mstNod,krbar,mstDofs2,iel,ierr)
          if (ierr < 0) goto 900
       end if

       !! Check dependent DOFs
       jp = mpmnpc(iel)
       call chkAndIncDepDofs (mmnpc(jp),slvDofs1,iel,ierr)
       if (ierr < 0) goto 900
       call chkAndIncDepDofs (mmnpc(jp+1),slvDofs2,iel,ierr)
       if (ierr < 0) goto 900
    end do

    !! Initialize current dependent (ICEQ) and pointer in MMCEQ/TTCC (IPCUR)
    iceq  = 1
    ipcur = 1

    !! =========================================================================
    !! rgdCE establishes the constraint equation for each dependent DOF of each
    !! RGD element in MPMCEQ/MMCEQ/TTCC.
    !! =========================================================================

    if (nrgd > 0) then
       call rgdCE (mprgd,elDepth,meldep,msc,mpmcex,mmcex,ttccx, &
            &      mpmceq,mmceq,mpmnpc,mmnpc,ttcc,madof,txc,tyc,tzc,minex, &
            &      iceq,ipcur,lpu,ierr)
       if (ierr < 0) goto 900
    end if

    !! =========================================================================
    !! rbarCE establishes the constraint equation for each dependent DOF of each
    !! RBAR element in MPMCEQ/MMCEQ/TTCC.
    !! =========================================================================

    if (nrbar > 0) then
       call rbarCE (mprbar,msc,mpmcex,mmcex,ttccx,mpmceq,mmceq,mpmnpc,mmnpc, &
            &       ttcc,madof,txc,tyc,tzc,minex,iceq,ipcur,lpu,ierr)
       if (ierr < 0) goto 900
    end if

    !! =========================================================================
    !! Append explicitly defined constraint equations to MPMCEQ/MMCEQ/TTCC, and
    !! expand/contract on any specified independent DOFs such that all
    !! independent DOFs included in the final version of the constraint arrays
    !! are free DOFs.
    !! =========================================================================

    if (ncex > 0) then
       maxsiz = min(size(ttcc),size(mmceq))
       call RGDXPL (mpar(1),minex(1),madof(1), &
            &       mpmcex(1),mmcex(1),ttccx(1),msc(1), &
            &       mpmceq(1),mmceq(1),ttcc(1), &
            &       ncex,maxsiz,lpu,iceq,ipcur,ierr)
       if (ierr < 0) goto 900
    end if

    !! Restore MSC, set NCEQ and NMMCEQ and complete MPMCEQ
    do i = 1, mpar(3)
       if (msc(i) < 0) then
          msc(i) = 0
       else if (msc(i) > 10) then
          msc(i) = msc(i) - 10
       end if
    end do
    mpar(7)  = iceq-1
    mpar(16) = ipcur-1
    mpmceq(iceq) = ipcur

    return

600 format('  ** Warning: Chained RGD elements:',I8)

900 continue
    call reportError (debugFileOnly_p,'Rigid3D')

  contains

    !> @brief Calculates chain depth of an RGD element based on array MELDEP.
    !> @details The function is invoked recursively until the end of a chain
    !> is detected. If RGD elements are connected in a loop, a negative value
    !> will be returned.
    recursive function calcDepth(ep) result (depth)
      integer, intent(in) :: ep
      integer             :: depth
      if (ep < 1 .or. ep > size(meldep)) then
         depth = 0 ! Invalid element number, ep
      else if (meldep(ep) == 0) then
         depth = 1 ! No other RGD element connected to this one
      else if (meldep(ep) > 0) then
         meldep(ep) = -meldep(ep) ! Negate current element entry to detect loops
         depth = 1 + calcDepth(-meldep(ep)) ! Check next element in the chain
         meldep(ep) = -meldep(ep) ! Restore positive values in meldep
         if (depth > 0) return    ! Return chain depth from current element ep
         !! A loop was detected, print all RGD elements in the loop
         write(lpu,"(I8)",ADVANCE='no') ffl_getElmId(mprgd(ep))
         depth = depth - 2
      else
         !! We have detected a loop, abort to avoid infinite recursion
         depth = -1
         write(lpu,610,ADVANCE='no') ffl_getElmId(mprgd(ep))
      end if
610   format(' *** Error: Looped RGD elements:',I8)
    end function calcDepth

    !> @brief Checks for illegal status codes on rigid independent DOFs.
    subroutine chkAndIncIndepDofs (inod,kind,dofs,iel,ierr)
      integer, intent(in)  :: inod, kind, iel
      logical, intent(in)  :: dofs(6)
      integer, intent(out) :: ierr

      integer :: i, isc, idof, nndof

      ierr  = 0
      idof  = madof(inod)
      nndof = madof(inod+1) - idof

      if (kind == krbar) then
         !! RBAR independent nodes may be 3-DOF nodes
         !! if none of the rotational DOFs are free
         if (nndof == 3 .and. any(dofs(4:6))) then
            write(lpu,600) iel,minex(inod)
         end if
      else if (kind == krigid) then
         !! RBE2 independent nodes must have all 6 DOFs
         if (nndof /= 6) then
            call RGDERR (13,nndof,minex(inod),iel,lpu,ierr)
            return
         end if
      end if

      do i = 1, nndof
         if (dofs(i)) then
            isc = msc(idof)
            if (isc > 0 .and. isc < 11) then
               msc(idof) = isc + 10
            else if (isc == -2) then
               call RGDERR (5,i,minex(inod),i,lpu,ierr)
               return
            end if
         end if
         idof = idof + 1
      end do

600   format('  ** Warning: Rotational independent DOFs for element',I8, &
           & ' used on 3-DOF node',I8)

    end subroutine chkAndIncIndepDofs

    !> @brief Checks for illegal status codes on rigid dependent DOFs.
    subroutine chkAndIncDepDofs (inod,dofs,iel,ierr)
      integer, intent(in)  :: inod, iel
      logical, intent(in)  :: dofs(6)
      integer, intent(out) :: ierr

      integer :: i, isc, idof, nndof

      idof  = madof(inod )
      nndof = madof(inod+1) - idof
      if (nndof /= 6 .and. nndof /= 3) then
         call RGDERR (14,nndof,minex(inod),iel,lpu,ierr)
         return
      end if

      do i = 1, nndof
         if (dofs(i)) then
            isc = msc(idof)
            if (isc > 0) then
               msc(idof) = isc + 10
            else if (isc == -1) then
               call RGDERR (10,1,minex(inod),i,lpu,ierr)
               return
            end if
         end if
         idof = idof + 1
      end do

    end subroutine chkAndIncDepDofs

  end subroutine Rigid3D


  !!============================================================================
  !> @brief Construct the constraint equations for RGD (or RBE2) elements.
  !>
  !> @param[in] mprgd Matrix of RGD element indices
  !> @param[in] elDepth Depth of each RGD element in chains
  !> @param[in] meldep RGD element couplings
  !> @param msc Matrix of status codes
  !> @param[in] mpmcex Matrix of pointers to external MCEQs
  !> @param[in] mmcex Matrix of external constraint equation definitions
  !> @param[in] ttccx Table of external constraint equation coefficients
  !> @param mpmceq Matrix of pointers to MCEQs
  !> @param mmceq Matrix of constraint equation definitions
  !> @param[in] mpmnpc Matrix of pointers to MNPCs
  !> @param[in] mmnpc Matrix of nodal point correspondances
  !> @param ttcc Table of constraint equation coefficients
  !> @param[in] madof Matrix of accumulated DOFs
  !> @param[in] txc Table of x-coordinates for the nodal points
  !> @param[in] tyc Table of Y-coordinates for the nodal points
  !> @param[in] tzc Table of Z-coordinates for the nodal points
  !> @param[in] minex Matrix of internal to external node number
  !> @param iceq Running constraint equation counter
  !> @param ipcur Running index into MMCEQ
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @details Chained RGD elements are detected and handled correctly
  !> (chains are traversed by depth such that all elements of depth i
  !> are processed prior to any element of depth i+1).
  !>
  !> @callgraph @callergraph
  !>
  !> @author Tommy Jorstad
  !>
  !> @date 18 Jul 2002

  subroutine rgdCE (mprgd,elDepth,meldep,msc,mpmcex,mmcex, &
       &            ttccx,mpmceq,mmceq,mpmnpc,mmnpc,ttcc,madof, &
       &            txc,tyc,tzc,minex,iceq,ipcur,lpu,ierr)

    use kindModule             , only : dp, epsDiv0_p
    use allocationModule       , only : reAllocate
    use reportErrorModule      , only : internalError
    use reportErrorModule      , only : reportError, debugFileOnly_p
    use FFlLinkHandlerInterface, only : ffl_getRgdDofComp

    integer , intent(in)    :: mprgd(:), elDepth(:), meldep(:)
    integer , intent(in)    :: mpmcex(:), mmcex(:), mpmnpc(:), mmnpc(:)
    integer , intent(in)    :: madof(:), minex(:), lpu
    real(dp), intent(in)    :: ttccx(:), txc(:), tyc(:), tzc(:)
    integer , intent(inout) :: msc(:), mpmceq(:), mmceq(:), ipcur, iceq
    real(dp), intent(inout) :: ttcc(:)
    integer , intent(out)   :: ierr

    !! Local variables
    integer             :: i, j, k, ip, jp, iel, idof, iDepth, iIns
    integer             :: maxDepth, mnumb, cesiz
    integer             :: nrgd, nndof, nceq, mstNod, slvNod, fMstDof, fSlvDof
    integer             :: depDof, status, rgdComp(4), ceInd(6)
    real(dp), pointer   :: c(:)
    integer , pointer   :: m(:)
    logical             :: foundCE, mstDofs(6), slvDofs(6), usedMst(6)

    integer , parameter :: dofMst(3,6) = reshape((/ 1, 5, 6, &
         &                                          2, 4, 6, &
         &                                          3, 4, 5, &
         &                                          4,-1,-1, &
         &                                          5,-1,-1, &
         &                                          6,-1,-1 /), (/ 3,6 /) )

    !! --- Logic section ---

    ierr = 0
    nrgd = size(mprgd)
    maxDepth = maxval(elDepth)
    nullify(c)
    nullify(m)

    !! Initial allocation of local dynamics arrays (may be increased below)
    call reAllocate ('rgdCE',m,3,ierr)
    call reAllocate ('rgdCE',c,3,ierr)
    if (ierr < 0) goto 900

    do iDepth = 1, maxDepth ! Loop over all element depths
       do i = 1, nrgd       ! Loop over all elements with this depth
          if (elDepth(i) /= iDepth) cycle

          iel     = mprgd(i)
          mstNod  = mmnpc(mpmnpc(iel))
          fMstDof = madof(mstNod)
          nndof   = madof(mstNod+1) - fMstDof

          !! Get DOF components for element's dependent DOFs
          call ffl_getRgdDofComp (rgdComp,iel,ierr)
          if (ierr < 0) goto 900

          call procDofComp (rgdComp(1),slvDofs,iel,lpu,ierr)
          if (ierr < 0) goto 900

          !! Check if the reference node of the element is a dependent node
          if (meldep(i) == 0) then
             !! Not a dependent node, so all DOFs are indepenedent
             mstDofs = .true.
          else
             !! Dependent node, so we check which DOFs are dependent DOFs
             mstDofs(1:nndof) = msc(fMstDof:fMstDof+nndof-1) > 0

             !! Then find which DOFs they depend on
             usedMst = .false.
             do j = 1, 6
                if (slvDofs(j)) then
                   do k = 1, 3
                      if (dofMst(k,j) > 0) usedMst(dofMst(k,j)) = .true.
                   end do
                end if
             end do

             !! If some DOFs depend on other DOFs that are dependent themselves,
             !! we need to find the constraint equation indices of those DOFs
             ceInd = 0
             do j = 1, 6
                if (usedMst(j) .and. .not. mstDofs(j)) then
                   idof = fMstDof + j - 1
                   foundCE = .false.
                   do k = 1, iceq
                      if (mmceq(mpmceq(k)) == idof) then
                         foundCE = .true.
                         ceInd(j) = k
                         exit
                      end if
                   end do
                   if (.not. foundCE) then
                      ierr = internalError('rgdCE: invalid SAM data structure')
                      return
                   end if
                end if
             end do
          end if

          !! At this stage,
          !! if the reference node is not dependent, then mstDofs is all true,
          !! if it is dependent, then mstDofs indicates which DOFs are free
          !! and ceInd contains the equation indices to find the CEs of the
          !! DOFs that are dependent, but still needed (there is no point in
          !! finding indices for DOFs that are not used by the dependent nodes).

          !! Now we loop over this element's dependent DOFs to set CEs
          do jp = mpmnpc(iel)+1, mpmnpc(iel+1)-1
             slvNod  = mmnpc(jp)
             fSlvDof = madof(slvNod)
             nndof   = madof(slvNod+1) - fSlvDof

             !! Check number of dofs on this node (both 6 and 3 are allowed)
             if (nndof == 6) then
                if (msc(fSlvDof+3)+msc(fSlvDof+4)+msc(fSlvDof+5) == 0) nndof = 3
             else
                nndof = 3
             end if

             !! Loop over the dependent DOFs on node slvNod
             do j = 1, nndof
                if (.not. slvDofs(j)) cycle

                depDof = fSlvDof + j - 1
                status = msc(depDof)

                !! Set up "local" constraining independent DOFs
                nceq = 3
                do k = 1, 3
                   if (dofMst(k,j) > 0) then
                      m(k) = fMstDof + dofMst(k,j) - 1
                   else
                      m(k) = 0
                   end if
                end do

                !! Set up "local" constraint coeffisients
                c(1) = 1.0_dp
                if (j == 1) then
                   c(2) = tzc(slvNod) - tzc(mstNod)
                   c(3) = tyc(mstNod) - tyc(slvNod)
                else if (j == 2) then
                   c(2) = tzc(mstNod) - tzc(slvNod)
                   c(3) = txc(slvNod) - txc(mstNod)
                else if (j == 3) then
                   c(2) = tyc(slvNod) - tyc(mstNod)
                   c(3) = txc(mstNod) - txc(slvNod)
                end if

                !! Loop over the reference node DOFs, if they are dependent,
                !! substitute their current entry by its governing CE
                do k = 1, 3
                   mnumb = dofMst(k,j)
                   !! Skip if DOF is not used in a coupling, or if it is free
                   if (mnumb <= 0)     cycle
                   if (mstDofs(mnumb)) cycle

                   !! Find size of CE to insert for this constrained DOF
                   if (ceInd(mnumb)+1 == iceq) then
                      cesiz = ipcur - mpmceq(ceInd(mnumb)) - 1
                   else
                      cesiz = mpmceq(ceInd(mnumb)+1) - mpmceq(ceInd(mnumb)) - 1
                   end if

                   if (nceq+cesiz-1 > size(c)) then
                      !! Reallocate the DOF and coefficients arrays
                      call reAllocate ('rgdCE',m,nceq+cesiz-1,ierr,.true.)
                      call reAllocate ('rgdCE',c,nceq+cesiz-1,ierr,.true.)
                      if (ierr < 0) goto 900
                   end if

                   !! iIns is starting index for insertion of the replacing CE
                   iIns = nceq+k-3

                   !! If not last mstdof, then move last part of m/c backwards
                   if (iIns < nceq) then
                      m(iIns+cesiz:nceq+cesiz-1) = m(iIns+1:nceq)
                      c(iIns+cesiz:nceq+cesiz-1) = c(iIns+1:nceq)
                   end if

                   !! Then m/c gets the CE of the DOF instead of the DOF itself
                   ip = mpmceq(ceInd(mnumb))
                   m(iIns:iIns+cesiz-1) =         mmceq(ip+1:ip+cesiz)
                   c(iIns:iIns+cesiz-1) = c(iIns)*(ttcc(ip+1:ip+cesiz))

                   nceq = nceq + cesiz-1
                end do

                !! Ignore coupling to DOFs with negligible coefficients
                !! (may yield potential ill-conditioning)
                do k = 1, nceq
                   if (abs(c(k)) < epsDiv0_p) m(k) = 0
                end do

                !! Lists m and c are now ok for this DOF, so store this CE
                call defineConstraintEqn (depDof,m(1:nceq),c, &
                     &                    mpmcex,mmcex,ttccx,msc, &
                     &                    mpmceq,mmceq,ttcc,iceq,ipcur,lpu,ierr)
                if (ierr < 0) then
                   call RGDERR (-2,minex(mstNod),minex(slvNod),i,lpu,i)
                   goto 900
                end if

                if (status == -2) then
                   nceq = size(mpmcex)-1
                   call RGDSWP (minex(1),msc(1),mpmcex(1),mmcex(1),ttccx(1), &
                        &       nceq,slvNod,depDof,lpu,ierr)
                   if (ierr < 0) goto 900
                end if

             end do !! dependent node DOF loop
          end do !! element dependent node loop

       end do !! element loop
    end do !! depth loop

    !! Release the local dynamic arrays
    call reAllocate ('rgdCE',m)
    call reAllocate ('rgdCE',c)

    return

900 continue
    call reportError (debugFileOnly_p,'rgdCE')

  end subroutine rgdCE


  !!============================================================================
  !> @brief Constructs the constraint equations for RBAR elements.
  !>
  !> @param[in] mprbar Matrix of RBAR element indices
  !> @param msc Matrix of status codes
  !> @param[in] mpmcex Matrix of pointers to external MCEQs
  !> @param[in] mmcex Matrix of external constraint equation definitions
  !> @param[in] ttccx Table of external constraint equation coefficients
  !> @param mpmceq Matrix of pointers to MCEQs
  !> @param mmceq Matrix of constraint equation definitions
  !> @param[in] mpmnpc Matrix of pointers to MNPCs
  !> @param[in] mmnpc Matrix of nodal point correspondances
  !> @param ttcc Table of constraint equation coefficients
  !> @param[in] madof Matrix of accumulated DOFs
  !> @param[in] txc Table of x-coordinates for the nodal points
  !> @param[in] tyc Table of Y-coordinates for the nodal points
  !> @param[in] tzc Table of Z-coordinates for the nodal points
  !> @param[in] minex Matrix of internal to external node number
  !> @param iceq Running constraint equation counter
  !> @param ipcur Running index into MMCEQ
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Tommy Jorstad
  !>
  !> @date 18 Jul 2002

  subroutine rbarCE (mprbar,msc,mpmcex,mmcex,ttccx,mpmceq,mmceq,mpmnpc,mmnpc, &
       &             ttcc,madof,txc,tyc,tzc,minex,iceq,ipcur,lpu,ierr)

    use kindModule             , only : dp
    use reportErrorModule      , only : reportError, error_p, debugFileOnly_p
    use FFlLinkHandlerInterface, only : ffl_getRgdDofComp

    integer , intent(in)    :: mprbar(:), mpmcex(:), mmcex(:)
    integer , intent(in)    :: mpmnpc(:), mmnpc(:), madof(:), minex(:), lpu
    real(dp), intent(in)    :: ttccx(:), txc(:), tyc(:), tzc(:)
    integer , intent(inout) :: msc(:), mpmceq(:), mmceq(:), ipcur, iceq
    real(dp), intent(inout) :: ttcc(:)
    integer , intent(out)   :: ierr

    !! Local variables
    real(dp), parameter :: eps_p = 1.0e-8_dp
    integer             :: i, j, k, iel, nndof, status, depDof, ncex, nrbar
    integer             :: slvNod, mstNod, fSlvDof, fMstDof, rgdComp(4), m(3)
    logical             :: mstDofs1(6), slvDofs1(6), mstDofs2(6), slvDofs2(6)
    logical             :: sdA(6), mdA(6), mdB(6)
    real(dp)            :: c(3), dx, dy, dz, tol
    character(len=128)  :: errMsg, errMsg2

    !! --- Logic section ---

    ierr  = 0
    ncex  = size(mpmcex)-1
    nrbar = size(mprbar)

    do i = 1, nrbar
       iel = mprbar(i)

       !! Get component definition
       call ffl_getRgdDofComp (rgdComp,iel,ierr)
       if (ierr < 0) goto 900

       call procRBarDef (rgdComp,mstDofs1,mstDofs2,slvDofs1,slvDofs2, &
            &            iel,lpu,ierr)
       if (ierr < 0) goto 900

       !! Set up CEs (dependent node corresponds to node A)
       mstNod = mmnpc(mpmnpc(iel))
       slvNod = mmnpc(mpmnpc(iel)+1)
       mdB    = mstDofs1
       mdA    = mstDofs2
       sdA    = slvDofs2

       do j = 1, 2

          fSlvDof = madof(slvNod)
          fMstDof = madof(mstNod)

          dx  = txc(mstNod) - txc(slvNod)
          dy  = tyc(mstNod) - tyc(slvNod)
          dz  = tzc(mstNod) - tzc(slvNod)
          tol = sqrt(dx*dx + dy*dy + dz*dz)*eps_p

          nndof = madof(slvNod+1) - fSlvDof
          if (nndof == 6) then
             if (msc(fSlvDof+3)+msc(fSlvDof+4)+msc(fSlvDof+5) == 0) nndof = 3
          else
             nndof = 3
          end if

          !! Loop over dependent DOFs, processing only those needed
          do k = 1, nndof
             if (.not. sdA(k)) cycle

             m      = 0
             c(1)   = 1.0_dp
             c(2:3) = 0.0_dp

             select case (k)
             case (1)
                if (mdB(1)) m(1) = fMstDof
                if (mdB(5)) m(2) = fMstDof + 4
                if (mdB(6)) m(3) = fMstDof + 5
                c(2) = -dz
                c(3) =  dy
             case (2)
                if (mdB(2)) m(1) = fMstDof + 1
                if (mdB(4)) m(2) = fMstDof + 3
                if (mdB(6)) m(3) = fMstDof + 5
                c(2) =  dz
                c(3) = -dx
             case (3)
                if (mdB(3)) m(1) = fMstDof + 2
                if (mdB(4)) m(2) = fMstDof + 3
                if (mdB(5)) m(3) = fMstDof + 4
                c(2) = -dy
                c(3) =  dx
             case (4)
                if (mdB(4)) then
                   m(1) = fMstDof + 3
                else if (mdA(2) .and. mdB(2) .and. abs(dz) > tol) then
                   m(1) = fSlvDof + 1
                   m(2) = fMstDof + 1
                   c(1) = -1.0_dp / abs(dz)
                   c(2) = -c(1)
                else if (mdA(3) .and. mdB(3) .and. abs(dy) > tol) then
                   m(1) = fSlvDof + 2
                   m(2) = fMstDof + 2
                   c(1) = 1.0_dp / abs(dy)
                   c(2) = -c(1)
                end if
             case (5)
                if (mdB(5)) then
                   m(1) = fMstDof + 4
                else if (mdA(1) .and. mdB(1) .and. abs(dz) > tol) then
                   m(1) = fSlvDof
                   m(2) = fMstDof
                   c(1) = 1.0_dp / abs(dz)
                   c(2) = -c(1)
                else if (mdA(3) .and. mdB(3) .and. abs(dx) > tol) then
                   m(1) = fSlvDof + 2
                   m(2) = fMstDof + 2
                   c(1) = -1.0_dp / abs(dx)
                   c(2) = -c(1)
                end if
             case (6)
                if (mdB(6)) then
                   m(1) = fMstDof + 5
                else if (mdA(1) .and. mdB(1) .and. abs(dy) > tol) then
                   m(1) = fSlvDof
                   m(2) = fMstDof
                   c(1) = -1.0_dp / abs(dy)
                   c(2) = -c(1)
                else if (mdA(2) .and. mdB(2) .and. abs(dx) > tol) then
                   m(1) = fSlvDof + 1
                   m(2) = fMstDof + 1
                   c(1) = 1.0_dp / abs(dx)
                   c(2) = -c(1)
                end if
             end select

             if (m(1) == 0) then
                write(errMsg,600) k, slvNod
                write(errMsg2,601) iel
                call reportError (error_p,errMsg,errMsg2)
                ierr = -1
                goto 900
             end if

             depDof = fSlvDof+k-1
             status = msc(depDof)
             call defineConstraintEqn (depDof,m,c,mpmcex,mmcex,ttccx, &
                  &                    msc,mpmceq,mmceq,ttcc,iceq,ipcur, &
                  &                    lpu,ierr)
             if (ierr < 0) then
                call RGDERR (-2,minex(mstNod),minex(slvNod),k,lpu,k)
                goto 900
             end if

             if (status == -2) then
                call RGDSWP (minex(1),msc(1),mpmcex(1),mmcex(1),ttccx(1),ncex, &
                     &       slvNod,depDof,lpu,ierr)
                if (ierr < 0) goto 900
             end if
          end do

          mstNod = mmnpc(mpmnpc(iel)+1)
          slvNod = mmnpc(mpmnpc(iel))
          mdB    = mstDofs2
          mdA    = mstDofs1
          sdA    = slvDofs1
       end do
    end do

    return

600 format('Could not establish constraint equation for DOF ',I1,', node',I8)
601 format('on RBAR element',I8,' due to the initial element geometry')

900 continue
    call reportError (debugFileOnly_p,'rbarCE')

  end subroutine rbarCE


  !!============================================================================
  !> @brief Establishes the constraint equation for a given dependent DOF.
  !>
  !> @param[in] iDofS The dependent DOF to generate constriant equation for
  !> @param[in] iDofM Array of independent DOFs to be coupled to DOF @a iDofS
  !> @param[in] C Array of constraint equation coefficients
  !> @param[in] mpmcex Matrix of pointers to external MCEQs
  !> @param[in] mmcex Matrix of external constraint equation definitions
  !> @param[in] ttccx Table of external constraint equation coefficients
  !> @param msc Matrix of status codes
  !> @param mpmceq Matrix of pointers to MCEQs
  !> @param mmceq Matrix of constraint equation definitions
  !> @param ttcc Table of constraint equation coefficients
  !> @param iceq Running constraint equation counter
  !> @param ipcur Running index into MMCEQ
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @details The dependent DOF (IDOFS) is coupled to the independent DOFs
  !> (DOF numbers in array IDOFM) through the following equation:
  !>
  !> @code
  !>    D-IDOFS = C1*(D-IDOFM1) + C2*(D-IDOFM2) + ... + C<n>*(D-IDOFM<n>)
  !> @endcode
  !>
  !> The equation is stored in the MMCEQ and TTCC arrays, from index IPCUR
  !> onwards. The IPCUR and ICEQ indices are then updated.
  !>
  !> This implementation is based on the SAM routine RGDCE3.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 22 Jul 2002

  subroutine defineConstraintEqn (iDofS,iDofM,C,mpmcex,mmcex,ttccx, &
       &                          msc,mpmceq,mmceq,ttcc,iceq,ipcur,lpu,ierr)

    use kindModule       , only : dp
    use reportErrorModule, only : reportError, debugFileOnly_p

    integer , intent(in)    :: iDofS, iDofM(:), mpmcex(:), mmcex(:), lpu
    real(dp), intent(in)    :: C(:), ttccx(:)
    integer , intent(inout) :: msc(:), mpmceq(:), mmceq(:), iceq, ipcur
    real(dp), intent(inout) :: ttcc(:)
    integer , intent(out)   :: ierr

    !! Local variables
    integer  :: i, j, ipceq, idof, ndofs, maxSiz
    logical  :: firstDof
    real(dp) :: C0

    !! --- Logic section ---

    ndofs  = size(iDofM)
    maxSiz = min(size(mmceq),size(ttcc))
    if (ipcur + ndofs > maxSiz) goto 900
    maxSiz = size(mpmceq)
    if (iceq > maxSiz) goto 900

    ierr         = 0
    firstDof     = .true.
    mpmceq(iceq) = 0
    mmceq(ipcur) = 0

    do i = 1, ndofs
       idof = iDofM(i)
       if (idof <= 0) cycle

       if (msc(idof) > 0) then

          if (firstDof .or. msc(iDofS) == 0) then
             !! Store dependent DOF and the zero constraint coefficient
             mpmceq(iceq) = ipcur
             mmceq(ipcur) = iDofS
             ttcc(ipcur)  = 0.0_dp
             ipcur        = ipcur + 1
          end if

          !! Store the i'th independent DOF and its associated coefficient
          mmceq(ipcur) = idof
          ttcc(ipcur)  = C(i)
          ipcur        = ipcur + 1
          firstDof     = .false.

          !! Mark this DOF as dependent
          msc(iDofS) = -3

       else if (msc(idof) == 0) then

          !! If the first independent DOF is fixed,
          !! the dependent DOF is also marked as fixed
          if (firstDof) msc(iDofS) = 0

       else if (msc(idof) == -1) then

          !! This independent DOF has an explicit prescribed value, get it
          C0 = 0.0_dp
          do j = 1, size(mpmcex)-1
             if (mmcex(mpmcex(j)) == idof) then
                C0 = ttccx(mpmcex(j))
                exit
             end if
          end do

          if (mpmceq(iceq) > 0) then
             ipceq        = mpmceq(iceq)
          else
             ipceq        = ipcur
             mpmceq(iceq) = ipcur
             ttcc(iceq)   = 0.0_dp
             ipcur        = ipcur + 1
          end if

          !! Add prescribed value of the independent DOF to this dependent DOF
          mmceq(ipceq) = iDofS
          ttcc(ipceq)  = ttcc(ipceq) + C(i)*C0
          firstDof     = .false.

          !! Mark the dependent DOF as having an (explicit) prescribed value
          msc(iDofS) = -1

       end if

    end do

    if (mpmceq(iceq) > 0) iceq = iceq + 1

    return

900 continue
    call RGDERR (4,maxSiz,0,0,lpu,ierr)
    call reportError (debugFileOnly_p,'defineConstraintEqn')

  end subroutine defineConstraintEqn


  !!============================================================================
  !> @brief Processes an m-set number of RGD and RBAR elements.
  !>
  !> @param[in] dofNumbs The m-set number
  !> @param[in] iel Element index
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] dofs Array of boolean DOF flags (on/off)
  !> @param[out] ierr Error flag
  !>
  !> @details If @a dofNumbs is zero, all DOFs are flagged as .true.
  !>
  !> @callergraph
  !>
  !> @author Tommy Jorstad
  !>
  !> @date 18 Jul 2002

  subroutine procDofComp (dofNumbs,dofs,iel,lpu,ierr)

    integer, intent(in)  :: dofNumbs, iel, lpu
    logical, intent(out) :: dofs(6)
    integer, intent(out) :: ierr

    !! Local variables
    integer :: i, iDof, div, dofIds

    !! --- Logic section ---

    div    = 100000
    dofs   = .false.
    dofIds = dofNumbs
    if (dofNumbs > 0) then
       do i = 1, 6
          iDof = dofIds/div
          if (iDof > 0 .and. iDof < 7) then
             dofs(iDof) = .true.
          else if (iDof /= 0) then
             ierr = -1
             write(lpu,600) iel,dofNumbs
             return
          end if
          dofIds = dofIds - iDof*div
          div = div/10
       end do
    else if (dofNumbs == 0) then
       dofs = .true.
    else if (dofNumbs < 0) then
       ierr = -2
       write(lpu,600) dofNumbs,iel
    end if

600 format(5X,'Invalid DOF components',I8,'  for rigid element',I8)

  end subroutine procDofComp


  !!============================================================================
  !> @brief Processes the component definition of an RBAR element.
  !>
  !> @param[in] rgdComp Dependent/independent DOFs at each node
  !> @param[in] iel Element index
  !> @param[out] d1 Array of independent DOFs at node 1
  !> @param[out] d2 Array of independent DOFs at node 2
  !> @param[out] d3 Array of dependent DOFs at node 1
  !> @param[out] d4 Array of dependent DOFs at node 2
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Tommy Jorstad
  !>
  !> @date 18 Jul 2002

  subroutine procRBarDef (rgdComp,d1,d2,d3,d4,iel,lpu,ierr)

    integer, intent(in)  :: rgdComp(4), iel, lpu
    logical, intent(out) :: d1(6), d2(6), d3(6), d4(6)
    integer, intent(out) :: ierr

    !! --- Logic section ---

    ierr = 0

    if (rgdComp(1) > 0) then
       call procDofComp (rgdComp(1),d1,iel,lpu,ierr)
    else
       d1 = .false.
    end if

    if (rgdComp(2) > 0) then
       call procDofComp (rgdComp(2),d2,iel,lpu,ierr)
    else
       d2 = .false.
    end if

    if (rgdComp(3) > 0) then
       call procDofComp (rgdComp(3),d3,iel,lpu,ierr)
    else if (rgdComp(4) > 0) then
       d3 = .false.
    end if

    if (rgdComp(4) > 0) then
       call procDofComp (rgdComp(4),d4,iel,lpu,ierr)
    else if (rgdComp(3) > 0) then
       d4 = .false.
    else
       d3 = .not. d1
       d4 = .not. d2
    end if

  end subroutine procRBarDef

end module RigidModule
