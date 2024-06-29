!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

subroutine fpp (ierr)

  use KindModule             , only : sp, dp, i8, lfnam_p
  use SamModule              , only : SamType, writeObject, deAllocateSAM
  use SamStressModule        , only : initiateSAM
  use InitiateStressModule   , only : readSolverData
  use DisplacementModule     , only : IVectorPtr, RMatrixPtr
  use DisplacementModule     , only : getFileName, ElDispFromSupElDisp
  use DisplacementModule     , only : openBandEmatrices, closeBAndEMatrices
  use DisplacementModule     , only : readResponsePointers
  use DisplacementModule     , only : readSupElDisplacements
  use DisplacementModule     , only : readStaticDisplacements, dis1Expand
  use TriadTypeModule        , only : TriadType
  use SupElTypeModule        , only : SupElType
  use RDBModule              , only : RDBType, writeTimeStepDB, closeRDBFile
#ifdef FT_HAS_FPPINTERFACE
  use FPPModule              , only : fppInit, fppOpen, fppClose
  use FPPModule              , only : fppProcess, fppProcessLast
#endif
  use FatigueModule          , only : fatigueInit, fatigueAddPoint
  use FatigueModule          , only : fatigueDamage
  use StrainCoatModule       , only : StrainCoatType, StrainElementType, iBuf
  use StrainCoatModule       , only : initiateStrainCoats, deallocateCoat
  use StrainCoatModule       , only : printStrainCoatInput, calcAngleData
  use StrainCoatModule       , only : printStrainCoatData, calcStrainCoatData
  use StrainRosetteModule    , only : initStrainRosette
  use StrainRosetteModule    , only : printRosetteHeading
  use StrainRosetteModule    , only : printRosetteStrains
  use SaveStrainCoatModule   , only : initiateStrainCoatHeader
  use SaveStrainCoatModule   , only : finalizeStrainCoatHeader
  use SaveStrainCoatModule   , only : writeElementsHeader, writeElementsDB
  use SaveStrainCoatModule   , only : writeHistoryHeader, writeHistoryDB
  use ResStressModule        , only : readResStress, getSolidSurfaceStress
  use TimerModule            , only : initTime, showTime
  use VersionModule          , only : openResFile
  use ProgressModule         , only : writeProgress, lterm
  use AllocationModule       , only : profile => doLogMem
  use PointerKindModule      , only : ptr
  use ReportErrorModule      , only : openTerminalOutputFile
  use ReportErrorModule      , only : allocationError, reportError
  use ReportErrorModule      , only : error_p, warning_p, debugFileOnly_p
  use FFaProfilerInterface   , only : ffa_initProfiler, ffa_reportTimer
  use FFaProfilerInterface   , only : ffa_startTimer, ffa_stopTimer
  use FFaCmdLineArgInterface , only : ffa_cmdlinearg_isSet
  use FFaCmdLineArgInterface , only : ffa_cmdlinearg_getint
  use FFaCmdLineArgInterface , only : ffa_cmdlinearg_getbool
  use FFaCmdLineArgInterface , only : ffa_cmdlinearg_getfloat
  use FFaCmdLineArgInterface , only : ffa_cmdlinearg_getdouble
  use FFlLinkHandlerInterface, only : ffl_init, ffl_done
  use FFlLinkHandlerInterface, only : ffl_getnostrc, ffl_getnostrcmat
  use FFrExtractorInterface  , only : ffr_init, ffr_done, ffr_getnextstep

  implicit none

  integer, intent(out)      :: ierr
  logical, parameter        :: useElementCoordSys_p = .true.
  logical                   :: lGrav, writeRosetteFiles, moreElements
  logical                   :: writeHistory, oldRange, addStress0, haveSNdata
  integer                   :: lpu, iprint, i, j, isup, iSurface, fppType
  integer                   :: nStrainCoat, nStrainCoatTotal, bufSiz, angBinSize
  integer(i8)               :: iStp
  integer(ptr), allocatable :: rPointers(:)
  real(sp)                  :: pvxGate
  real(dp)                  :: biAxialGate, toMPaScale, sigma3D(6), g(3)
  real(dp)                  :: currTime, lastTime, startTime, stopTime, tInc
  real(dp), allocatable     :: vgii(:), vii(:)
  type(RDBType)             :: rdb
  type(SamType)             :: samSup
  type(SupElType)           :: sup
  type(TriadType), pointer  :: triads(:)
  character(len=lfnam_p)    :: chname, modelFileName
  character(len=32)         :: errMsg

  type(StrainElementType), pointer :: rosette
  type(StrainCoatType),    pointer :: strainCoats(:)
  type(RMatrixPtr),    allocatable :: H_all(:)
  type(IVectorPtr),    allocatable :: globalNodes_all(:)

  !! --- Logic section ---

  call initTime ()
  call ffa_initprofiler ('fedem_fpp profiler')
  call ffa_starttimer ('fedem_fpp')
  call openTerminalOutputFile (lterm)
  write(lterm,6000) 'START OF PROGRAM FPP'

  ! Open the fpp result file and write heading
  call getFileName ('resfile',chname,'_fpp.res')
  call openResFile (chname,'Damage Recovery',lpu,ierr)
  if (ierr /= 0) return

  ! Get output options from the command-line
  call ffa_cmdlinearg_getint ('debug',iprint)
  call ffa_cmdlinearg_getint ('surface',iSurface)
  call ffa_cmdlinearg_getint ('BufSizeInc',bufSiz)
  call ffa_cmdlinearg_getint ('angleBins',angBinSize)
  call ffa_cmdlinearg_getint ('HistDataType',fppType)
  call ffa_cmdlinearg_getbool ('oldRange',oldRange)
  call ffa_cmdlinearg_getbool ('writeHistory',writeHistory)
  call ffa_cmdlinearg_getfloat ('PVXGate',pvxGate)
  call ffa_cmdlinearg_getdouble ('biAxialGate',biAxialGate)
  call ffa_cmdlinearg_getdouble ('stressToMPaScale',toMPaScale)
  call ffa_cmdlinearg_getdouble ('statm',startTime)
  call ffa_cmdlinearg_getdouble ('stotm',stopTime)
  call ffa_cmdlinearg_getdouble ('tinc',tInc)
  writeRosetteFiles = iprint > 1

  ! Import residual stresses from external file, if specified
  addStress0 = .false.
  if (ffa_cmdlinearg_isSet('resStressFile')) then
     call writeProgress (' --> Importing residual stresses',profile)
     call readResStress (ierr)
     if (ierr < 0) then
        call reportError (error_p,'Error reading external stress file')
        goto 900
     else
        addStress0 = ierr == 0
     end if
  end if

  ! --- Read the link file

  call writeProgress (' --> Reading link files',profile)
  call ffa_starttimer ('ffl_init')
  call ffl_init (.true.,ierr)
  call ffa_stoptimer ('ffl_init')
  if (ierr /= 0) then
     call reportError (error_p,'Error initializing the link object')
     goto 900
  end if

  nStrainCoatTotal = ffl_getnostrc()
  if (nStrainCoatTotal < 1) then
     call reportError (error_p,'Link does not contain any strain coat elements')
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
  write(lterm,"(15X,'Number of strain coats =',I6)") nStrainCoatTotal
  call readSolverData (isup,g,sup,triads,modelFileName,iprint,lpu,ierr)
  if (ierr /= 0) goto 900

  if (iprint > 4) then
     call WriteObject (samSup,lpu,3)
     write(lpu,"(/)")
  end if

  ! --- Open the B-matrix and the generalized modes files

  call openBandEmatrices (samSup%ndof1,samSup%ndof2, &
       &                  samSup%mpar(22),samSup%mpar(31),ierr)
  if (ierr < 0) goto 900

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

  ! --- Open the solver results database

  call writeProgress (' --> Reading solver result files',profile)
  call ffa_starttimer ('ffr_init')
  call ffr_init ('frsfile',ierr)
  call ffa_stoptimer ('ffr_init')
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

  ! --- Initialize the strain coat results database

  if (.not. writeHistory) then ! Write summary frs-file

     call writeProgress (' --> Initializing result database headers')
     call initiateStrainCoatHeader (modelFileName,rdb,sup,ierr)
     if (ierr /= 0) goto 900

     ! The binary data is written to a temporary file
     call writeTimeStepDB (rdb,1_i8,stopTime,ierr,.true.)
     if (ierr < 0) goto 900

  end if

  ! --- Allocate strain coat elements

  call ffa_cmdlinearg_getint ('blockSize',nStrainCoat)
  if (nStrainCoat > nStrainCoatTotal) then
     nStrainCoat = nStrainCoatTotal
  else if (nStrainCoat < 1) then
     nStrainCoat = min(2000,nStrainCoatTotal)
     write(errMsg,"(' Reset to',I5)") nStrainCoat
     call reportError (warning_p,'Invalid value on option -blockSize.'//errMsg)
  end if

  allocate(strainCoats(nStrainCoat),STAT=ierr)
  if (ierr /= 0) then
     ierr = allocationError('strain coat elements')
     goto 900
  end if

  ! --- Open and initialize the fpp-file

  if (fppType > 0) then

     call fatigueInit (ierr)
     if (ierr < 0) goto 900

#ifdef FT_HAS_FPPINTERFACE
  else if (fppType < 0) then

     call writeProgress (' --> Initializing fpp-file')
     call getFileName ('fppfile',chname,'.fpp')
     call fppOpen (chname,sup%id,stopTime,ffl_getnostrcmat(),ierr)
     if (ierr < 0) goto 900
#endif

  end if

  ! --- Element block loop

  haveSNdata = .false.
  moreElements = .true.
  do while (moreElements)

     call writeProgress (' --> Initializing strain coats',profile)
     call ffa_starttimer ('initStrainCoats')
     call InitiateStrainCoats (isup,iSurface,fppType,bufSiz,angBinSize, &
          &                    strainCoats,nStrainCoat,ierr)
     call ffa_stoptimer ('initStrainCoats')
     if (ierr < 0) goto 900
     if (ierr > 0) moreElements = .false.
     if (nStrainCoat < 1) exit

     write(lterm,6010) strainCoats(1)%tensorRosette%id%baseId, &
          &            strainCoats(nStrainCoat)%tensorRosette%id%baseId, &
          &            nStrainCoatTotal

     allocate(globalNodes_all(nStrainCoat),H_all(nStrainCoat),stat=ierr)
     if (ierr /= 0) then
        ierr = allocationError('globalNodes_all and H_all')
        goto 900
     end if

     do i = 1, nStrainCoat
        globalNodes_all(i)%p => strainCoats(i)%tensorRosette%globalNodes
     end do

     ! --- Compute the H_el matrices for all strain rosettes

     call writeProgress (' --> Computing transformation matrices',profile)
     call ffa_starttimer ('elDispFromSupElDisp')
     call ElDispFromSupElDisp (samSup,globalNodes_all,H_all,ierr)
     call ffa_stoptimer ('elDispFromSupElDisp')
     if (ierr < 0) goto 900

     deallocate(globalNodes_all)

     ! --- Initialize the strain rosettes

     call writeProgress (' --> Computing strain-displ matrices',profile)
     do i = 1, nStrainCoat

        if (strainCoats(i)%fppType > 0) haveSNdata = .true.

        if (addStress0) then
           ! Get surface residual stress tensor of the underlying solid element
           rosette => strainCoats(i)%tensorRosette
           call getSolidSurfaceStress (samSup, &
                &                      rosette%elmNumber, &
                &                      rosette%globalNodes, &
                &                      sigma3D, ierr)
           if (ierr < 0) goto 900
        else
           ierr = 0
        end if

        call ffa_starttimer ('initStrainRosette')
        if (lGrav) then
           call InitStrainRosette (strainCoats(i)%tensorRosette, H_all(i)%p, &
                &                  samSup%madof, ierr>0, sigma3D, 0,lpu, ierr, &
                &                  useElementCoordSys_p, writeRosetteFiles, &
                &                  vgii=vgii)
        else
           call InitStrainRosette (strainCoats(i)%tensorRosette, H_all(i)%p, &
                &                  samSup%madof, ierr>0, sigma3D, 0,lpu, ierr, &
                &                  useElementCoordSys_p, writeRosetteFiles)
        end if
        call ffa_stoptimer ('initStrainRosette')
        if (ierr < 0) goto 900

        if (ierr > 0) then
           ! Degenerated strain coat element, just ignore it
           call deAllocateCoat(strainCoats(i))
           allocate(strainCoats(i)%tensorRosette%data(0))
           allocate(strainCoats(i)%results(0))
           allocate(strainCoats(i)%fppProcessorHandle(0))
        else if (writeRosetteFiles) then
           rosette => strainCoats(i)%tensorRosette
           call PrintRosetteHeading (rosette%data(1), &
                &                    rosette%globalNodes, samSup%minex, &
                &                    rosette%posInGl, rosette%id%userId, &
                &                    rosette%linkNumber, rosette%lpuAscii)
        end if

        deallocate(H_all(i)%p)
     end do
     deallocate(H_all)

     if (iprint > 0) then

        ! Write strain coat data nicely to ASCII results file
        call printStrainCoatInput (strainCoats(1:nStrainCoat),lpu)

     end if

     if (writeHistory) then ! Write history frs-files

        call writeProgress (' --> Writing result database headers')
        call getFileName ('rdbfile',chname,'.frs')
        call writeHistoryHeader (chname,modelFileName,rdb,sup, &
             &                   strainCoats(1:nStrainCoat),ierr)
        if (ierr /= 0) goto 900

     end if
#ifdef FT_HAS_FPPINTERFACE
     if (fppType < 0) then

        call writeProgress (' --> Initializing fpp processors',profile)
        call fppInit (strainCoats(1:nStrainCoat),ierr)
        if (ierr < 0) goto 900

     end if
#endif

     ! --- Time loop

     iBuf = 0
     lastTime = stopTime
     currTime = startTime - 1.0_dp
     call writeProgress (' --> Starting time loop',profile)
     do while (ffr_getnextstep(startTime,stopTime,tInc,currTime,iStp))

        iBuf = iBuf + 1
        if (iprint > 2 .and. iStp < 10000000_i8) write(lpu,6011) iStp,currTime

        ! Read system response data for this time step
        call ffa_starttimer ('readSupElDisp')
        call readSupElDisplacements (triads,sup,rPointers,iprint,lpu,ierr)
        call ffa_stoptimer ('readSupElDisp')
        if (ierr /= 0) goto 900

        ! Calculate strain coat results
        do i = 1, nStrainCoat
           if (.not. associated(strainCoats(i)%tensorRosette%data)) cycle

           call ffa_starttimer ('calcStrainCoatData')
           call CalcStrainCoatData (strainCoats(i), sup%finit, &
                &                   biAxialGate, toMPaScale, ierr)
           call ffa_stoptimer ('calcStrainCoatData')
           if (ierr /= 0) goto 900

           if (writeRosetteFiles) then
              rosette => strainCoats(i)%tensorRosette
              call PrintRosetteStrains (rosette%data(1), currTime, &
                   &                    rosette%lpuAscii)
           end if
        end do

        if (writeHistory) then

           call writeTimeStepDB (rdb,iStp,currTime,ierr)
           if (ierr < 0) goto 900

           call writeHistoryDB (rdb,strainCoats(1:nStrainCoat),ierr)
           if (ierr < 0) goto 900

        end if
        if (fppType > 0) then

           ! Accumulate time histories for fatigue calculation
           call fatigueAddPoint (strainCoats(1:nStrainCoat),currTime)

#ifdef FT_HAS_FPPINTERFACE
        else if (fppType < 0 .and. iBuf == bufSiz) then

           ! Calculate fpp results
           call fppProcess (strainCoats(1:nStrainCoat),iBuf,ierr)
           if (ierr < 0) goto 900

           iBuf = 0
#endif
        end if

        lastTime = currTime
        call writeProgress (currTime-startTime,stopTime-startTime)

     end do
     if (iStp < -1_i8) then
        ierr = -1
        write(chname,"('time =',1PE12.5)") currTime
        call reportError (error_p,'Error searching for results after '//chname)
        goto 900
     else if (lastTime < stopTime) then
        call writeProgress (stopTime-startTime,stopTime-startTime)
     end if

     ! --- Calculate most popular angle and angle spread

     do i = 1, nStrainCoat
        do j = 1, size(strainCoats(i)%results)
           call calcAngleData (strainCoats(i)%results(j)%angBin, &
                &              strainCoats(i)%results(j)%popAngle, &
                &              strainCoats(i)%results(j)%angSpread, &
                &              strainCoats(i)%results(j)%strRange, oldRange)
        end do
     end do

     ! --- Close all FPP processors

     call writeProgress (' --> Time loop done. Closing fpp processors',profile)
     if (writeHistory) then

        call closeRDBfile (rdb,ierr)
        if (ierr < 0) goto 900

     end if
     if (fppType > 0) then

        ! Calculate damage
        call fatigueDamage (strainCoats(1:nStrainCoat),real(pvxGate,dp),ierr)
        if (ierr < 0) goto 900

#ifdef FT_HAS_FPPINTERFACE
     else if (fppType < 0) then

        ! Calculate fpp results
        call fppProcessLast (strainCoats(1:nStrainCoat),iBuf,biAxialGate,ierr)
        if (ierr < 0) goto 900
#endif

     end if
     if (iprint > 0) then

        ! Write strain coat results nicely to ASCII results file
        call printStrainCoatData (strainCoats(1:nStrainCoat),lpu)

     end if
     if (.not. writeHistory) then

        ! --- Write strain coat summary to database file

        call writeProgress (' --> Writing to result database',profile)
        call writeElementsHeader (rdb,strainCoats(1:nStrainCoat),.false.,ierr)
        if (ierr /= 0) goto 900

        call writeElementsDB (rdb,strainCoats(1:nStrainCoat), &
             &                stopTime-startTime,ierr)
        if (ierr < 0) goto 900

     end if

     do i = 1, nStrainCoat
        j = strainCoats(i)%tensorRosette%lpuAscii
        if (j > 0) close(j)
        call deallocateCoat (strainCoats(i))
     end do

  end do

  ! --- Finito, deallocate and close down

  call writeProgress (' --> Block loop done. Closing database files',profile)

  if (.not. haveSNdata .and. fppType > 0) then
     call reportError (warning_p,'None of the element groups processed '// &
          &            'were assigned an S-N curve.', &
          &            'Damage and Life contour plots are therefore not '// &
          &            'created for this part.')
  end if

  if (.not. writeHistory) then
     call getFileName ('rdbfile',chname,'.frs')
     call finalizeStrainCoatHeader (chname,rdb,ierr)
     call closeRDBfile (rdb,ierr)
     if (ierr < 0) goto 900
  end if

#ifdef FT_HAS_FPPINTERFACE
  if (fppType < 0) then
     call fppClose (ierr)
  end if
#endif

  deallocate(strainCoats)
  deallocate(triads,rPointers)
  if (lGrav) deallocate(vgii)

  call deAllocateSAM (samSup)
  call closeBandEmatrices (ierr)
  call ffl_done
  call ffr_done

900 continue
  if (ierr /= 0) then
     call writeProgress (' --> Done, with errors',profile)
  else
     call writeProgress (' --> Done',profile)
  end if
  write(lterm,6000) ' END OF PROGRAM FPP '
  if (ierr /= 0) then
     call reportError (debugFileOnly_p,'fpp')
     call getFileName ('resfile',chname,'_fpp.res')
     write(lterm,6100) chname
     call showTime (lpu)
     write(lpu,6200) 'failed :-('
  else
     call showTime (lpu)
     write(lpu,6200) 'successfully completed :-)'
  end if
  call closeLogFiles ()
  call ffa_stoptimer ('fedem_fpp')
  call ffa_reporttimer

6000 format(/11X,16('='),'> ',A,' <',16('='))
6010 format( 15X,'Processing elements',I8,'  to',I8,'  of total',I8,' elements')
6011 format(/120('*')//5X,'TIME STEP',I8,' : INTEGRATION TIME T =',1PD13.5)
6100 format(/11X,'An error is detected. Error messages are written to ',A)
6200 format(/ 4X,'Strain coat calculation ',A)

end subroutine fpp
