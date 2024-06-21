!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

subroutine gage (ierr)

  use KindModule                , only : dp, i8, epsDiv0_p, lfNam_p
  use SamModule                 , only : SamType, deAllocateSAM, writeObject
  use SamStressModule           , only : initiateSAM
  use DisplacementModule        , only : openBandEmatrices, closeBandEmatrices
  use DisplacementModule        , only : readSupElDisplacements, readSolverData
  use DisplacementModule        , only : readResponsePointers, getFileName
  use DisplacementModule        , only : readStaticDisplacements, dis1Expand
  use DisplacementModule        , only : ElDispFromSupElDisp
  use DisplacementModule        , only : RMatrixPtr, IVectorPtr
  use IdTypeModule              , only : getId
  use TriadTypeModule           , only : TriadType
  use SupElTypeModule           , only : SupElType
  use StrainGageModule          , only : InitStrainGages, ReadStrainGageData
  use StrainGageModule          , only : PrintStrainGages_FiDF
  use StrainGageModule          , only : addFatiguePoints, reportDamage
  use StrainGageModule          , only : deallocateFatigue
  use StrainRosetteModule       , only : StrainElementType
  use StrainRosetteModule       , only : InitStrainRosette
  use StrainRosetteModule       , only : PrintRosetteHeading
  use StrainRosetteModule       , only : PrintInitialStrain, PrintRosetteStrains
  use StrainRosetteModule       , only : calcRosetteDisplacements
  use StrainRosetteModule       , only : calcZeroStartRosetteStrains
  use StrainRosetteModule       , only : calcRosetteStrains, deallocateRosette
  use RDBModule                 , only : RDBType, flushRDBfile, closeRDBfile
  use SaveStrainGageModule      , only : writeStrainGageHeader,writeStrainGageDB
  use ManipMatrixModule         , only : writeObject
  use TimerModule               , only : initTime, showTime
  use VersionModule             , only : openResFile
  use ProgressModule            , only : lterm, writeProgress
  use PointerKindModule         , only : ptr
  use ReportErrorModule         , only : note_p, error_p, debugFileOnly_p
  use ReportErrorModule         , only : openTerminalOutputFile
  use ReportErrorModule         , only : allocationError, reportError
  use FiDeviceFunctionInterface , only : FiDF_CloseAll
  use FFaCmdLineArgInterface    , only : ffa_cmdlinearg_getint
  use FFaCmdLineArgInterface    , only : ffa_cmdlinearg_getbool
  use FFaCmdLineArgInterface    , only : ffa_cmdlinearg_getdouble
  use FFaCmdLineArgInterface    , only : ffa_cmdlinearg_isSet
  use FFlLinkHandlerInterface   , only : ffl_init, ffl_done
  use FFrExtractorInterface     , only : ffr_init, ffr_done, ffr_getnextstep

  implicit none

  integer, intent(out)      :: ierr
  logical                   :: nullStartStrains, writeAsciiFiles, useFSIformat
  logical                   :: lDeformation, lGrav
  integer                   :: iFatigue, iprint, isup, i, lpu, nRosettes
  integer(i8)               :: iStep
  integer(ptr), allocatable :: rPointers(:)
  real(dp)                  :: currTime, lastTime, startTime, stopTime, tInc
  real(dp)                  :: flushInc, lastFlush, bufRatio, tIncDac
  real(dp)                  :: toMPaScale, g(3)
  real(dp), allocatable     :: vgii(:), vii(:)
  type(RDBType)             :: rdb
  type(SamType)             :: samSup
  type(SupElType)           :: sup
  type(TriadType), pointer  :: triads(:)
  character(len=lfnam_p)    :: chname, modelFileName

  type(StrainElementType), pointer :: strainRosettes(:)
  type(RMatrixPtr),    allocatable :: H_all(:)
  type(IVectorPtr),    allocatable :: globalNodes_all(:)

  !! --- Logic section ---

  call initTime ()
  call openTerminalOutputFile (lterm)
  write(lterm,6000) 'START OF PROGRAM GAGE'

  ! Open the gage result file and write heading
  call getFileName ('resfile',chname,'_gage.res')
  call openResFile (chname,'Strain Gage Recovery',lpu,ierr)
  if (ierr /= 0) return

  call ffa_cmdlinearg_getint ('debug',iprint)
  call ffa_cmdlinearg_getint ('fatigue',iFatigue)
  call ffa_cmdlinearg_getbool ('nullify_start_rosettestrains',nullStartStrains)
  call ffa_cmdlinearg_getbool ('writeAsciiFiles',writeAsciiFiles)
  call ffa_cmdlinearg_getbool ('deformation',lDeformation)
  call ffa_cmdlinearg_getdouble ('stressToMPaScale',toMPaScale)
  call ffa_cmdlinearg_getdouble ('statm',startTime)
  call ffa_cmdlinearg_getdouble ('stotm',stopTime)
  call ffa_cmdlinearg_getdouble ('tinc',tInc)
  call ffa_cmdlinearg_getdouble ('dac_sampleinc',tIncDac)
  if (tIncDac >= 0.0_dp .and. tIncDac < tInc) tIncDac = tInc

  ! Compute result file buffer size ratio
  call ffa_cmdlinearg_getdouble ('flushinc',flushInc)
  if (flushInc <= 0.0_dp) then
     bufRatio = 0.0_dp
  else if (tInc > epsDiv0_p) then
     bufRatio = flushInc/tInc
  else
     bufRatio = min(1000.0_dp,flushInc/0.0001_dp) ! We don't know the time step
  end if

  ! --- Read the link file

  call writeProgress (' --> Reading link files')
  call ffl_init (0,0,.false.,ierr)
  if (ierr /= 0) then
     call reportError (error_p,'Error initializing the link object')
     goto 900
  end if

  ! --- Establish the SAM datastructure

  call getFileName ('samfile',chname,'.fsm')
  call initiateSAM (chname,samSup,lpu,ierr)
  if (ierr /= 0) goto 900

  ! --- Read superelement data from the solver input file

  call ffa_cmdlinearg_getint ('linkId',isup)
  if (isup < 1) isup = samSup%mpar(18) ! Fetch baseID from reducer if not given
  write(lterm,"(/11X,'--> Process Link; baseID (isup) =',I6)") isup
  call readSolverData (isup,g,sup,triads,modelFileName,iprint,lpu,ierr)
  if (ierr /= 0) goto 900

  if (iprint > 4) then
     call WriteObject (samSup,lpu,3)
     write(lpu,"(/)")
  end if

  ! --- Read definition for possible strain gages on this link

  call writeProgress (' --> Initializing strain rosettes')
  call getFileName ('rosfile',chname)
  useFSIformat = chname(len_trim(chname)-3:) == '.fsi'
  call readStrainGageData (chname,useFSIformat,sup%id,strainRosettes,ierr)
  if (ierr /= 0) goto 900

  if (.not. associated(strainRosettes)) then
     call reportError (note_p,'No strain rosettes on this link')
     deallocate(triads)
     call deAllocateSAM (samSup)
     call ffl_done
     goto 900
  end if

  ! --- Open the B-matrix and the generalized modes files

  call openBandEmatrices (samSup%ndof1,samSup%ndof2, &
       &                  samSup%mpar(22),samSup%mpar(31),ierr)
  if (ierr < 0) goto 900

  nRosettes = size(strainRosettes)
  allocate(globalNodes_all(nRosettes),H_all(nRosettes),stat=ierr)
  if (ierr /= 0) then
     ierr = allocationError('globalNodes_all and H_all')
     goto 900
  end if

  do i = 1, nRosettes
     globalNodes_all(i)%p => strainRosettes(i)%globalNodes
  end do

  ! --- Compute the H_el matrices for all strain rosettes

  call ElDispFromSupElDisp (samSup,globalNodes_all,H_all,ierr)
  if (ierr < 0) goto 900

  deallocate(globalNodes_all)

  ! --- Read in static displacements due to gravitation forces, if any

  lGrav = sqrt(g(1)*g(1)+g(2)*g(2)+g(3)*g(3)) > 1.0e-8_dp
  if (lGrav .and. ffa_cmdlinearg_isSet('dispfile')) then
     call getFileName ('dispfile',chname,'_V.fmx')
     inquire(FILE=chname,EXIST=lGrav)
  else
     lGrav = .false.
  end if
  if (lGrav) then
     allocate(vgii(samSup%ndof),vii(samSup%ndof1),STAT=ierr)
     if (ierr /= 0) then
        ierr = allocationError('gravitation deflection vector')
        goto 900
     end if

     call readStaticDisplacements (chname,g,vii,lpu,ierr)
     if (ierr < 0) goto 900

     call dis1Expand (samSup,vii,vgii,ierr)
     if (ierr < 0) goto 900

     deallocate(vii)
  end if

  ! --- Initialize the strain rosettes

  do i = 1, nRosettes

     if (useFSIformat) then
        ! New input format, assume the position matrix is given explicitly
        if (lGrav) then
           call InitStrainRosette (strainRosettes(i),H_all(i)%p,samSup%madof, &
                &                  .false.,calcDisp=.true.,vgii=vgii, &
                &                  openFiles=writeAsciiFiles, &
                &                  lpu=lpu,iprint=iprint,ierr=ierr)
        else
           call InitStrainRosette (strainRosettes(i),H_all(i)%p,samSup%madof, &
                &                  .false.,calcDisp=.true., &
                &                  openFiles=writeAsciiFiles, &
                &                  lpu=lpu,iprint=iprint,ierr=ierr)
        end if
     else
        ! Old input format, must recompute the rosette position matrix
        if (lGrav) then
           call InitStrainRosette (strainRosettes(i),H_all(i)%p,samSup%madof, &
                &                  .false.,useElCoordSys=.false.,vgii=vgii, &
                &                  openFiles=writeAsciiFiles, &
                &                  lpu=lpu,iprint=iprint,ierr=ierr)
        else
           call InitStrainRosette (strainRosettes(i),H_all(i)%p,samSup%madof, &
                &                  .false.,useElCoordSys=.false., &
                &                  openFiles=writeAsciiFiles, &
                &                  lpu=lpu,iprint=iprint,ierr=ierr)
        end if
     end if
     if (ierr /= 0) goto 900

     call InitStrainGages (strainRosettes(i)%id%userId, &
          &                strainRosettes(i)%data(1)%gages, &
          &                strainRosettes(i)%data(1)%alphaGages, &
          &                strainRosettes(i)%posInGl, tIncDac, ierr)
     if (ierr < 0) goto 900

     call PrintRosetteHeading (strainRosettes(i)%data(1), &
          &                    strainRosettes(i)%globalNodes, samSup%minex, &
          &                    strainRosettes(i)%posInGl, &
          &                    strainRosettes(i)%id%userId, &
          &                    strainRosettes(i)%linkNumber, &
          &                    strainRosettes(i)%lpuAscii)

  end do

  if (.not.useFSIformat) then
     do i = 1, nRosettes
        deallocate(H_all(i)%p)
     end do
     deallocate(H_all)
  end if

  if (lGrav) deallocate(vgii)

  call closeBandEmatrices (ierr)
  if (ierr < 0) goto 900

  call ffl_done

  ! --- Open the solver results database

  call writeProgress (' --> Reading solver result files')
  call ffr_init ('frsfile',ierr)
  if (ierr /= 0) then
     call reportError (error_p,'Error opening results database')
     goto 900
  end if

  allocate(rPointers(size(triads)+2),STAT=ierr)
  if (ierr /= 0) then
     ierr = allocationError('item group pointers')
     goto 900
  end if

  call readResponsePointers (triads,sup,rPointers,ierr)
  if (ierr /= 0) goto 900

  ! --- Initialize the strain gage results database

  call writeProgress (' --> Writing result database headers')
  call getFileName ('rdbfile',chname,'.frs')
  call writeStrainGageHeader (chname,modelFileName,rdb,strainRosettes, &
       &                      bufRatio,ierr)
  if (ierr /= 0) goto 900

  ! --- Time loop

  lastTime = stopTime
  currTime = startTime - 1.0_dp
  lastFlush = currTime
  call writeProgress (' --> Starting time loop')
  do while (ffr_getnextstep(startTime,stopTime,tInc,currTime,iStep))

     if (iprint > 0 .and. iStep < 100000000_i8) then
        write(lterm,6010) currTime
        write(lpu,6011) iStep,currTime
     end if

     ! Read system response data for this time step
     call readSupElDisplacements (triads,sup,rPointers,iprint,lpu,ierr)
     if (ierr /= 0) goto 900

     if (iprint > 3) then
        call writeObject (sup%finit,lpu,'  ** Nodal deformations (Finit)')
     end if

     ! Calculate strain gage results
     if (iprint > 0) write(lterm,6030)
     do i = 1, nRosettes

        if (strainRosettes(i)%data(1)%zeroInit .or. nullStartStrains) then

           call CalcZeroStartRosetteStrains (strainRosettes(i)%data(1), &
                &                            sup%finit,ierr)
           if (ierr /= 0) goto 900

           call PrintInitialStrain (strainRosettes(i)%data(1),currTime, &
                &                   strainRosettes(i)%lpuAscii)

        end if

        if (associated(strainRosettes(i)%disp)) then
           call CalcRosetteDisplacements (strainRosettes(i),samSup%madof, &
                &                         H_all(i)%p,sup%finit,sup%supTr, &
                &                         lpu,ierr)
           if (ierr /= 0) goto 900
           if (iprint > 1) then
              write(lpu,6031) trim(getId(strainRosettes(i)%id)), &
                   &          strainRosettes(i)%disp
           end if
        end if

        call CalcRosetteStrains (strainRosettes(i)%data(1),sup%finit,ierr)
        if (ierr /= 0) goto 900

        call PrintRosetteStrains (strainRosettes(i)%data(1),currTime, &
             &                    strainRosettes(i)%lpuAscii)

        call PrintStrainGages_FiDF (strainRosettes(i)%data(1),currTime)

        if (iprint > 1) then
           write(lpu,"('  ** Strain rosette',A/6X,1P6E13.5)") &
                &    trim(getId(strainRosettes(i)%id)), &
                &    strainRosettes(i)%data(1)%epsC, &
                &    strainRosettes(i)%data(1)%sigmaC
        end if

        if (iFatigue > 0) then
           call AddFatiguePoints (strainRosettes(i)%data(1),currTime,toMPaScale)
        end if

     end do
     nullStartStrains = .false.

     call writeStrainGageDB (rdb,strainRosettes,iStep,currTime, &
          &                  lDeformation,ierr)
     if (ierr < 0) goto 900

     if (flushInc >= 0.0_dp .and. currTime >= lastFlush+flushInc) then
        lastFlush = currTime
        call flushRDBfile (rdb,ierr)
        if (ierr < 0) goto 900
     end if

     if (iprint < 1) then
        lastTime = currTime
        call writeProgress (currTime-startTime,stopTime-startTime)
     end if

  end do
  if (iStep < -1_i8) then
     ierr = -1
     write(chname,"('time =',1PE12.5)") currTime
     call reportError (error_p,'Error searching for results after '//chname)
     goto 900
  else if (iprint < 1 .and. lastTime < stopTime) then
     call writeProgress (stopTime-startTime,stopTime-startTime)
  end if

  ! --- Finito, deallocate and close down

  if (allocated(H_all)) then
     do i = 1, size(H_all)
        deallocate(H_all(i)%p)
     end do
     deallocate(H_all)
  end if

  if (iFatigue > 0) then
     call writeProgress (' --> Time loop done. Performing fatigue calculation')
     call reportDamage (strainRosettes,iFatigue,lpu)
     call writeProgress (' --> Closing database files')
  else
     call writeProgress (' --> Time loop done. Closing database files')
  end if

  do i = 1, size(strainRosettes)
     if (strainRosettes(i)%lpuAscii > 0) close(strainRosettes(i)%lpuAscii)
     call deallocateFatigue (strainRosettes(i)%data(1))
     call deallocateRosette (strainRosettes(i))
  end do
  deallocate(strainRosettes)
  deallocate(triads,rPointers)
  call deAllocateSAM (samSup)
  call closeRDBfile (rdb,ierr)
  call FiDF_CloseAll
  call ffr_done

900 continue
  write(lterm,6000) ' END OF PROGRAM GAGE '
  if (ierr /= 0) then
     call reportError (debugFileOnly_p,'gage')
     call getFileName ('resfile',chname,'_gage.res')
     write(lterm,6100) chname
     call showTime (lpu)
     write(lpu,6200) 'failed :-('
  else
     call showTime (lpu)
     write(lpu,6200) 'successfully completed :-)'
  end if
  call closeLogFiles ()

6000 format(/11X,16('='),'> ',A,' <',16('='))
6010 format(/11X,'--> ......Simulation time :',1PD13.5)
6011 format(/120('*')//5X,'TIME STEP',I8,' : INTEGRATION TIME T =',1PD13.5)
6030 format( 11X,'--> .........................Strain gage results')
6031 format(/5X,'Nodal deformations for Strain Rosette', A / (5X,1P3E13.5))
6100 format(/11X,'An error is detected. Error messages are written to ',A)
6200 format(/ 4X,'Strain gage calculation ',A)

end subroutine gage
