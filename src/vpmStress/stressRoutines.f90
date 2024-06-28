!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file stressRoutines.f90
!> @brief Calculation and storage of superelement stresses.

!!==============================================================================
!> @brief Module with subroutine for superelement stress calculation/storage.

module StressRoutinesModule

  implicit none

contains

  !!============================================================================
  !> @brief Calculates element stresses and stores on results database file.
  !>
  !> @param[in] addStress0 Include residual stresses from external source
  !> @param rdb The results database file to write stress results to
  !> @param[in] sam Assembly management data for the superelement
  !> @param[in] vtf VTF file handle
  !> @param[in] vtfId VTF result block ID
  !> @param[in] isup Based ID of the superelement
  !> @param[in] sv Superelement displacement vector
  !> @param[out] sf Superelement nodal forces
  !> @param[out] vms von Mises stress state array
  !> @param[in] irec Recovery timer handle
  !> @param[in] isav Saving timer handle
  !> @param[in] iprint Print switch
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @details The von Mises stress is also written to the VTF-file,
  !> if such output is specified. The von Mises stress can also be returned
  !> in the state array @a vms, for use by external visualization tools.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 11 Oct 2000

  subroutine calcStresses (addStress0,rdb,sam,vtf,vtfId,isup,sv,sf,vms, &
       &                   irec,isav,iprint,lpu,ierr)

    use KindModule                    , only : dp, i8, hugeVal_p
    use RecKindModule                 , only : rk
    use RDBModule                     , only : RDBType, flushRDBfile
    use SamModule                     , only : SamType
    use ElStressModule                , only : extractEV, elStress
    use ResStressModule               , only : addResStress
    use SaveStressModule              , only : saveElmOrder, writeStressHeader
    use SaveStressModule              , only : writeStressDB, writeStrMeasureDB
    use SaveVTFModule2                , only : VTFA_RESMAP_ELEMENT
    use SaveVTFModule2                , only : VTFA_RESMAP_ELEMENT_NODE
    use SaveVTFModule2                , only : VTF_RESSET_ABSMAX
    use SaveVTFModule2                , only : VTFFResBlockCreate
    use SaveVTFModule2                , only : VTFFResBlockDelete
    use SaveVTFModule2                , only : VTFFResBlockSetMapToBlockID
    use SaveVTFModule2                , only : VTFFResBlockSetNumRes
    use SaveVTFModule2                , only : VTFFFileWriteBlock
    use SaveVTFModule2                , only : writeStressVTF
    use TimerModule                   , only : startTimer, stopTimer
    use ReportErrorModule             , only : internalError, reportError
    use ReportErrorModule             , only : error_p, debugFileOnly_p
    use FFaTensorTransformsInterface  , only : vonMises, princval, maxshearvalue
    use FFaCmdLineArgInterface        , only : ffa_cmdlinearg_getbool
    use FFaCmdLineArgInterface        , only : ffa_cmdlinearg_isTrue

    logical          , intent(in)    :: addStress0
    type(RDBType)    , intent(inout) :: rdb
    type(SamType)    , intent(in)    :: sam
    integer ,optional, intent(in)    :: vtf, vtfId, isup, irec, isav
    real(rk),optional, intent(in)    :: sv(:)
    integer          , intent(in)    :: iprint, lpu
    real(dp),optional, intent(out)   :: sf(:)
    real(dp),optional, intent(out)   :: vms(:)
    integer          , intent(out)   :: ierr

    !! Local variables
    integer, parameter :: iReSet = VTF_RESSET_ABSMAX ! hard-coded for now
    integer, parameter :: maxstp = 50, maxnod = 10, nStrRes = 8
    logical, parameter :: lSS  = .false. ! No stress resultant conjugate strains
    logical     :: lDouble, lForce, lSave, lStoreVMS
    logical     :: lSR, lStress, lStrain, lStrRes(nStrRes)
    integer     :: i, j, iel, jel, idat, iBlock, lerr, iReMap, vtfFile
    integer     :: n, ncmp0, ncmp1, ncmp2, nenod, nstrp
    integer(i8) :: nBytes, nWrote
    real(dp)    :: EV(6*maxnod), EF(6*maxnod), SR(6*maxnod), SS(6*maxnod)
    real(dp)    :: Stress(6*maxstp), Strain(6*maxstp), resMat(nStrRes,maxstp)

    !! --- Logic section ---

    if (present(isup)) then
       call ffa_cmdlinearg_getbool ('double',lDouble)
       call ffa_cmdlinearg_getbool ('SR',lSR)
       call ffa_cmdlinearg_getbool ('stress',lStress)
       call ffa_cmdlinearg_getbool ('strain',lStrain)
       call ffa_cmdlinearg_getbool ('vmStress',lStrRes(1))
       call ffa_cmdlinearg_getbool ('vmStrain',lStrRes(5))
       call ffa_cmdlinearg_getbool ('maxPStress',lStrRes(2))
       call ffa_cmdlinearg_getbool ('maxPStrain',lStrRes(6))
       call ffa_cmdlinearg_getbool ('minPStress',lStrRes(3))
       call ffa_cmdlinearg_getbool ('minPStrain',lStrRes(7))
       call ffa_cmdlinearg_getbool ('maxSStress',lStrRes(4))
       call ffa_cmdlinearg_getbool ('maxSStrain',lStrRes(8))
       lSave = lSR .or. lStress .or. lStrain .or. any(lStrRes)
       lStoreVMS = .false.
    else ! Solver mode, recover von Mises stress only
       lDouble = .false.
       lSR = .false.
       lStress = .false.
       lStrain = .false.
       lStrRes(1) = .true.
       lStrRes(2:) = .false.
       lSave = present(isav)
       if (present(vms)) then
          lStoreVMS = size(vms) > 0
       else
          lStoreVMS = .false.
       end if
    end if
    if (present(SF)) then
       call ffa_cmdlinearg_getbool ('nodalForces',lForce)
       if (lForce) SF = 0.0_dp
    else
       lForce = .false.
    end if
    if (present(vtf) .and. present(vtfId) .and. present(isup) .and. lSave) then
       vtfFile = vtf
    else
       vtfFile = -1
    end if

    if (vtfFile >= 0) then
       if (present(isav)) call startTimer (isav)

       if (ffa_cmdlinearg_isTrue('VTFavgelm')) then
          iReMap = VTFA_RESMAP_ELEMENT ! Export averaged element results
       else
          iReMap = VTFA_RESMAP_ELEMENT_NODE ! Export element node results
       end if

       !! Create a scalar VTF block for the von Mises stress
       iBlock = VTFFResBlockCreate(vtfId,1,iReMap,0)
       if (iBlock < 0) then
          ierr = iBlock
          call reportError (error_p,'Error creating scalar result block', &
               &            ierr=ierr)
          goto 500
       end if
       call VTFFResBlockSetMapToBlockID (iBlock,isup)
       ierr = VTFFResBlockSetNumRes(iBlock,sam%nmmnpc) - 1
       if (ierr < 0) then
          call reportError (error_p,'Error allocating scalar result block', &
               &            ierr=ierr)
          goto 500
       end if

       if (present(isav)) call stopTimer (isav)
    end if

    idat = 0 ! Index counter for the von Mises Stress vector
    nBytes = rdb%nBytes
    do jel = 1, sam%nel
       if (allocated(saveElmOrder)) then
          iel = saveElmOrder(jel)
       else
          iel = jel
       end if
       if (present(irec)) call startTimer (irec)

       SR = 0.0_dp
       SS = 0.0_dp
       Stress = 0.0_dp
       Strain = 0.0_dp
       nenod = 0
       nstrp = 0
       ncmp0 = numberOfComponents(sam%melcon(iel),0)
       ncmp1 = numberOfComponents(sam%melcon(iel),1)
       ncmp2 = numberOfComponents(sam%melcon(iel),2)

       if (present(sv)) then

          !! Extract element nodal displacements
          call extractEV (iel,sam,sv,ev,lerr)
          if (lerr < 0) then
             call reportError (error_p,'Too many element nodes', &
                  &            addString='extractEV')
          end if

          !! Compute element results
          if (lForce .and. lerr > 0) then
             call elStress (iel,sam%melcon(iel),EV,SR,SS,Stress,Strain,EF, &
                  &         nenod,nstrp,iprint,lpu,lerr)
          else if (lerr > 0) then
             call elStress (iel,sam%melcon(iel),EV,SR,SS,Stress,Strain, &
                  &         nenod=nenod,nstrp=nstrp, &
                  &         ipsw=iprint,lpu=lpu,ierr=lerr)
          end if

       else
          lerr = 0
       end if

       if (addStress0 .and. lerr >= 0) then
          !! Add residual stresses to the stress tensor
          call addResStress (sam,Stress,iel,nenod,nstrp,iprint,lpu,lerr)
       end if

       if (present(irec)) call stopTimer (irec)
       if (lerr < 0) then
          ierr = lerr
          goto 510
       else if (ncmp1*nstrp > 6*maxstp) then
          ierr = internalError('calcStresses: Overflow of arrays Stress,Strain')
          return
       else if (ncmp0*nenod > 6*maxnod) then
          ierr = internalError('calcStresses: Overflow of arrays EF,SR,SS')
          return
       end if

       if (lForce .and. nenod > 0) then
          !! Add element nodal forces into the global vector
          n = sam%mpmnpc(iel)
          call addEF (SF,EF,sam%madof,sam%mmnpc(n:n+nenod-1),nenod,ncmp0)
       end if

       if (ncmp2 > 0 .and. nenod > 0 .and. lSave) then

          if (present(isav)) call startTimer (isav)
          if (lerr > 0) then
             !! Stress calculation failed for this element, write hugeVal
             if (lSR) SR = hugeVal_p
             if (lSS) SS = hugeVal_p
          end if

          !! Save stress resultants and conjugate strains on database file
          call writeStressDB (rdb,lSR,lSS,lDouble,ncmp2,nenod,SR,SS,ierr)
          if (ierr < 0) then
             call reportError (error_p, &
                  &            'Error saving element stress resultants')
             goto 500
          else if (iprint > 6) then
             nWrote = rdb%nBytes - nBytes
             nBytes = rdb%nBytes
             write(lpu,600) iel,nWrote,nBytes
          end if
          if (present(isav)) call stopTimer (isav)

       end if
       if (ncmp1 > 0 .and. nstrp > 0) then

          if (any(lStrRes) .and. nstrp > maxstp) then
             ierr = internalError('calcStresses: Overflow of array resMat')
             return
          end if

          if (lerr > 0) then
             !! Stress calculation failed, write hugeVal
             if (lStress)      Stress = hugeVal_p
             if (lStrain)      Strain = hugeVal_p
             if (any(lStrRes)) resMat = hugeVal_p
          else
             if (present(irec)) call startTimer (irec)

             !! Calculate von Mises and/or principal values
             n = numberOfComponents(sam%melcon(iel),3)
             j = 1
             resMat = 0.0_dp
             do i = 1, nstrp
                if (lStrRes(1) .or. vtfFile >= 0) then
                   resMat(1,i) = vonMises(n,Stress(j))
                end if
                if (any(lStrRes(2:4))) then
                   call princval (n,Stress(j),resMat(2:4,i))
                   call maxshearvalue (n,resMat(2:4,i),SS(1))
                   resMat(3,i) = resMat(1+n,i)
                   resMat(4,i) = SS(1)
                end if
                if (lStrRes(5)) then
                   resMat(5,i) = vonMises(n,Strain(j))
                end if
                if (any(lStrRes(6:8))) then
                   call princval (n,Strain(j),resMat(6:8,i))
                   call maxshearvalue (n,resMat(6:8,i),SS(1))
                   resMat(7,i) = resMat(5+n,i)
                   resMat(8,i) = SS(1)
                end if
                j = j + ncmp1
             end do

             if (present(irec)) call stopTimer (irec)
          end if

          if (lSave) then
             if (present(isav)) call startTimer (isav)

             !! Save stress/strain tensors and derived values on database file
             call writeStrMeasureDB (rdb,lStress,lStrain,lStrRes,lDouble, &
                  &                  ncmp1,nStrRes,nstrp,Stress,Strain, &
                  &                  resMat,ierr)
             if (ierr < 0) then
                call reportError (error_p,'Error saving element stresses')
                goto 500
             else if (iprint > 6) then
                nWrote = rdb%nBytes - nBytes
                nBytes = rdb%nBytes
                write(lpu,600) iel,nWrote,nBytes
             end if

             if (vtfFile >= 0) then
                if (lerr == 0) then
                   !! Write von Mises stresses to the VTF-file
                   nstrp = nstrp / nenod
                   call writeStressVTF (iBlock,nenod,nstrp,iReSet,iReMap, &
                        &               resMat(1,:),ierr)
                else
                   !! Write undefined value to the VTF-file for failed elements
                   call writeStressVTF (iBlock,nenod,2,0,iReMap,resMat,ierr)
                end if
                if (ierr < 0) goto 500
             end if

             if (present(isav)) call stopTimer (isav)
          end if

          if (lStoreVMS) then
             !! Store von Mises stresses in core array for all elements
             vms(idat+1) = iel  !     ffl_getElmId(iel) ?
             vms(idat+2) = nenod
             vms(idat+3) = nstrp
             vms(idat+4:idat+3+nstrp) = resMat(1,1:nstrp)
             idat = idat + 3 + nstrp
          end if

       else if (ncmp1 > 0 .and. vtfFile >= 0) then

          !! Write undefined value to the VTF-file for the stress-less elements
          nenod = sam%mpmnpc(iel+1) - sam%mpmnpc(iel)
          call writeStressVTF (iBlock,nenod,2,0,iReMap,resMat,ierr)
          if (ierr < 0) goto 510

       end if

    end do
    if (lSave) then
       if (present(isav)) call startTimer (isav)

       !! Time step finished, flush the results database file to disk
       call flushRDBfile (rdb,ierr)
       if (ierr < 0) goto 500

       if (vtfFile >= 0) then

          !! Write the scalar VTF block to file
          ierr = VTFFFileWriteBlock(vtfFile,iBlock) - 1
          if (ierr < 0) then
             call reportError (error_p,'Error writing scalar result block', &
                  &            ierr=ierr)
             goto 500
          end if
          call VTFFResBlockDelete (iBlock)

       end if

       if (present(isav)) call stopTimer (isav)
    end if
    if (lForce .and. iprint > 0) then

       !! Output the nodal forces to log-file
       call printNodalForces (sam%nnod,sam%madof,sam%minex,SF,lpu)

    end if

    return

500 if (present(isav)) call stopTimer (isav)
510 call reportError (debugFileOnly_p,'calcStresses')

600 format('Element',I8,' : wrote',I6,' bytes. Totally written',I10,' bytes.')

  contains

    !> @brief Returns the number of result components for an element type.
    function numberOfComponents (ieltyp,restyp) result(nComp)
      integer, intent(in) :: ieltyp, restyp
      integer             :: nComp
      select case (restyp)
      case (0) ! Nodal forces
         select case (ieltyp)
         case (11);    nComp = 6 ! Beam element
         case (21:24); nComp = 6 ! Thin shell element
         case (41:46); nComp = 3 ! Solid element
         case default; nComp = 0
         end select
      case (1) ! Stresses and strains
         select case (ieltyp)
         case (11);    nComp = 1 ! Beam element
         case (21:24); nComp = 3 ! Thin shell element
         case (31:32); nComp = 6 ! Thick shell element
         case (41:46); nComp = 6 ! Solid element
         case default; nComp = 0
         end select
      case (2) ! Stress resultants and conjugate strains
         select case (ieltyp)
         case (11);    nComp = 6 ! Beam element
         case (21:24); nComp = 6 ! Thin shell element
         case (31:32); nComp = 6 ! Thick shell element
         case default; nComp = 0
         end select
      case (3) ! Stress/strain tensor dimension
         select case (ieltyp)
         case (11);    nComp = 1 ! Beam element
         case (21:24); nComp = 2 ! Thin shell element
         case (31:32); nComp = 3 ! Thick shell element
         case (41:46); nComp = 3 ! Solid element
         case default; nComp = 0
         end select
      case default
         nComp = 0
      end select
    end function numberOfComponents

    !> @brief Adds element nodal forces into the global vector.
    subroutine addEF (SF,EF,madof,mnpc,nenod,nndof)
      integer , intent(in)    :: nenod,nndof,madof(:),mnpc(nenod)
      real(dp), intent(in)    :: EF(nndof,nenod)
      real(dp), intent(inout) :: SF(:)
      integer                 :: i, j, k
      do i = 1, nenod
         j = madof(mnpc(i))
         k = j + nndof - 1
         SF(j:k) = SF(j:k) + EF(:,i)
      end do
    end subroutine addEF

    !> @brief Prints nodal forces to unit LPU, omitting nodes with zero forces.
    subroutine printNodalForces (nnod,madof,minex,sf,lpu)
      integer , intent(in) :: nnod, madof(:), minex(:), lpu
      real(dp), intent(in) :: sf(:)
      integer              :: i, j, k
      real(dp)             :: tolf(6)
      real(dp), parameter  :: epsTol_p = 1.0e-8_dp

      tolf = 0.0_dp
      do i = 1, size(madof)-1
         do j = madof(i), madof(i+1)-1
            k = j+1 - madof(i)
            tolf(k) = max(tolf(k),abs(sf(j)))
         end do
      end do
      tolf = tolf * epsTol_p

      write(lpu,600)
      do i = 1, nnod
         do j = madof(i), madof(i+1)-1
            if (abs(sf(j)) > tolf(j+1-madof(i))) goto 100
         end do
         cycle
100      write(lpu,610,advance='NO') minex(i)
         do j = madof(i), madof(i+1)-1
            if (abs(sf(j)) > tolf(j+1-madof(i))) then
               write(lpu,620,advance='NO') sf(j)
            else
               write(lpu,630,advance='NO')
            end if
         end do
         write(lpu,*)
      end do
      write(lpu,650) tolf

600   format(//5X,'NODAL FORCES IN SUPERELEMENT :', &
           &      ' (values less than Tolerance are not printed)' &
           & //5X,'Node #    Fx',11X,'Fy',11X,'Fz',11X,'Mx',11X,'My',11X,'Mz' )
610   format(I9,2X)
620   format(1P,E13.5)
630   format(13(' '))
650   format('  Tolerance',1P6E13.5/)
    end subroutine printNodalForces

  end subroutine calcStresses

end module StressRoutinesModule

