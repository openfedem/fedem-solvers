!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

subroutine stress (ierr)

  use KindModule             , only : dp, i8, maxInt_p, lfNam_p
  use SamModule              , only : SamType, deAllocateSAM, writeObject
  use SamStressModule        , only : initiateSAM
  use InitiateStressModule   , only : readSolverData
  use DisplacementModule     , only : openBandEmatrices, closeBandEmatrices
  use DisplacementModule     , only : readSupElPosition, readSupElDisplacements
  use DisplacementModule     , only : readResponsePointers, readDisplPointer
  use DisplacementModule     , only : readIntDisplacements, calcIntDisplacements
  use DisplacementModule     , only : calcTotalDisplacements, getFileName
  use TriadTypeModule        , only : TriadType
  use SupElTypeModule        , only : SupElType, nullifySupel
  use RDBModule              , only : RDBType, writeTimeStepDB, closeRDBfile
  use SaveStressModule       , only : initWriteDisp, writeStressHeader
  use SaveStressModule       , only : writeDisplacementDB, saveElmOrder
  use SaveNastranModule2     , only : writeSupElDeformation
  use SaveNastranModule2     , only : findFENodeNumsForTriads
  use SaveVTFModule2         , only : VTFFFileSetOutputDebugError
  use SaveVTFModule2         , only : VTFFStateInfoBlockCreate
  use SaveVTFModule2         , only : VTFFFileAppendFile
  use SaveVTFModule2         , only : VTFFFileCloseFile, VTFFFileDelete
  use SaveVTFModule2         , only : writeTimeStepVTF, writeDummyStressVTF
  use SaveVTFModule2         , only : writeSupElTransformVTF, writeFooterVTF
  use SaveVTFModule2         , only : writeDisplacementVTF
  use ResStressModule        , only : readResStress
  use StressRoutinesModule   , only : calcStresses
  use TimerModule            , only : initTime, showTime
  use VersionModule          , only : openResFile
  use ProgressModule         , only : lterm, writeProgress
  use PointerKindModule      , only : ptr
  use ReportErrorModule      , only : note_p, error_p, debugFileOnly_p
  use ReportErrorModule      , only : openTerminalOutputFile
  use ReportErrorModule      , only : allocationError, reportError
  use BinaryDBInterface      , only : readDoubleDB
  use FFaCmdLineArgInterface , only : ffa_cmdlinearg_getint
  use FFaCmdLineArgInterface , only : ffa_cmdlinearg_getbool
  use FFaCmdLineArgInterface , only : ffa_cmdlinearg_getdouble
  use FFaCmdLineArgInterface , only : ffa_cmdlinearg_getstring
  use FFaCmdLineArgInterface , only : ffa_cmdlinearg_isSet
  use FFaCmdLineArgInterface , only : ffa_cmdlinearg_isTrue
  use FFaCmdLineArgInterface , only : ffa_cmdlinearg_intValue
  use FFlLinkHandlerInterface, only : ffl_init, ffl_done
  use FFlLinkHandlerInterface, only : ffl_export_vtf, ffl_elmorder_vtf
  use FFrExtractorInterface  , only : ffr_init, ffr_done
  use FFrExtractorInterface  , only : ffr_setposition, ffr_getnextstep

  implicit none

  integer, intent(out)     :: ierr
  integer(i8)              :: iStep
  integer                  :: lpu, iprint, isup
  integer                  :: vtfFile, vtfInfo, idVTF(2)
  integer(ptr),allocatable :: rPointers(:)
  integer(ptr)             :: disPtr
  logical                  :: lDeformation, lDumpDef, lGrav, addStress0, lVTFAt0
  real(dp)                 :: currTime, startTime, stopTime, tInc, dsVTF, grv(3)
  real(dp), allocatable    :: sv(:), sf(:), vii(:)
  type(RDBType)            :: rdb
  type(SamType)            :: samSup
  type(SupElType)          :: sup
  type(TriadType), pointer :: triads(:)
  character(len=lfnam_p)   :: chname, modelFileName

  !! --- Logic section ---

  call initTime ()
  call openTerminalOutputFile (lterm)
  write(lterm,6000) 'START OF PROGRAM STRESS'

  ! Open the stress results file and write heading
  call getFileName ('resfile',chname,'_stress.res')
  call openResFile (chname,'Stress Recovery',lpu,ierr)
  if (ierr /= 0) return

  call ffa_cmdlinearg_getint ('debug',iprint)
  call ffa_cmdlinearg_getbool ('dumpDefNas',lDumpDef)
  call ffa_cmdlinearg_getbool ('deformation',lDeformation)
  call ffa_cmdlinearg_getdouble ('statm',startTime)
  call ffa_cmdlinearg_getdouble ('stotm',stopTime)
  call ffa_cmdlinearg_getdouble ('tinc',tInc)

  ! Import residual stresses from external file, if specified
  addStress0 = .false.
  if (ffa_cmdlinearg_isSet('resStressFile')) then
     if (startTime > 0.0_dp) then
        call writeProgress (' --> Importing residual stresses')
        call readResStress (ierr)
        if (ierr < 0) then
           call reportError (error_p,'Error reading external stress file')
           goto 900
        else
           addStress0 = ierr == 0
        end if
     else
        call reportError (note_p,'Residual stresses are ignored when '// &
             &            'start time <= zero')
     end if
  end if

  ! --- Read the link file

  call writeProgress (' --> Reading link files')
  call ffl_init (0,0,.true.,ierr)
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
  write(lterm,"(/11X,'--> Process Part; baseID (isup) =',I6)") isup
  if (ffa_cmdlinearg_isSet('fsifile')) then
     call readSolverData (isup,grv,sup,triads,modelFileName,iprint,lpu,ierr)
     if (ierr /= 0) goto 900
  else
     allocate(triads(0)) ! No supernodes, assuming direct solution is available
     call nullifySupEl (sup)
     call ffa_cmdlinearg_getstring ('linkfile',sup%id%descr)
     sup%id%baseId = isup
  end if

  if (iprint > 4) then
     call WriteObject (samSup,lpu,3)
     write(lpu,"(/)")
  end if

  if (lDumpDef) then
     ! Find FE node numbers for the triads connected to this superelement
     call findFENodeNumsForTriads (samSup,sup,ierr)
     if (ierr < 0) goto 900
  end if

  ! --- Open the B-matrix and the generalized modes files

  if (ffa_cmdlinearg_isSet('Bmatfile')) then
     call openBandEmatrices (samSup%ndof1,samSup%ndof2, &
          &                  samSup%mpar(22),samSup%mpar(31),ierr)
     if (ierr < 0) goto 900
  end if

  allocate(sv(samSup%ndof), STAT=ierr)
  if (ierr /= 0) then
     ierr = allocationError('displacement vector')
     goto 900
  end if
  if (lDeformation .or. ffa_cmdlinearg_isTrue('nodalForces')) then
     allocate(sf(samSup%ndof), STAT=ierr)
     if (ierr /= 0) then
        ierr = allocationError('force vector')
        goto 900
     end if
  end if

  ! --- Read in static displacements due to gravitation forces, if any

  lGrav = sqrt(grv(1)*grv(1)+grv(2)*grv(2)+grv(3)*grv(3)) > 1.0e-8_dp
  if (lGrav .and. ffa_cmdlinearg_isSet('dispfile')) then
     call getFileName ('dispfile',chname,'_V.fmx')
     inquire(FILE=chname,EXIST=lGrav)
  else
     lGrav = .false.
  end if
  if (lGrav) then
     allocate(vii(samSup%ndof1*3), STAT=ierr)
     if (ierr /= 0) then
        ierr = allocationError('gravitation deflection vector')
        goto 900
     end if
     call readDoubleDB (chname,'displacement matrix',samSup%ndof1*3,vii(1),ierr)
     if (ierr < 0) then
        call reportError (error_p,'Failed to read matrix file '//chname, &
             &            addString='readDoubleDB')
        goto 900
     end if
  end if

  ! --- Output part geometry to the VTF-file, if required

  chname = ''
  if (ffa_cmdlinearg_isSet('VTFfile')) then
     call ffa_cmdlinearg_getstring ('VTFfile',chname)
  end if

  if (chname /= '') then
     call writeProgress (' --> Writing part geometry to VTF-file')
     call ffl_export_vtf (trim(chname),trim(sup%id%descr),sup%id%baseId,ierr)
     if (ierr /= 0) then
        call reportError (error_p,'Error writing to VTF-file '//chname)
        goto 900
     end if

     ! Reopen the VTF-file for fortran output
     call VTFFFileSetOutputDebugError (1)
     vtfFile = VTFFFileAppendFile(trim(chname))
     if (vtfFile < 0) then
        ierr = vtfFile
        call reportError (error_p,'Error re-opening VTF-file '//chname)
        goto 900
     end if

     ! Initialize VTF result block id numbers
     idVTF = ffa_cmdlinearg_intValue('VTFoffset')
     if (ffa_cmdlinearg_intValue('VTFparts') > 0) then
        vtfInfo = VTFFStateInfoBlockCreate()
     else
        vtfInfo = -1
     end if

     allocate(saveElmOrder(samSup%nel), STAT=ierr)
     if (ierr /= 0) then
        ierr = allocationError('element permutation vector')
        goto 900
     end if

     call ffl_elmorder_vtf (saveElmOrder(1),ierr)
     if (ierr > 0) then
        call reportError (note_p,'The element processing order is changed '// &
             &            'to optimize the VTF export')
     else
        deallocate(saveElmOrder)
     end if

     call ffa_cmdlinearg_getdouble ('VTFdscale',dsVTF)
     lVTFAt0 = startTime > 0.0_dp .and. ffa_cmdlinearg_isTrue('VTFinit')

  else
     vtfFile = -1
     lVTFAt0 = .false.
  end if

  ! --- Open the solver results database

  call writeProgress (' --> Reading solver result files')
  call ffr_init ('frsfile',ierr)
  if (ierr /= 0) then
     call reportError (error_p,'Error opening results database')
     goto 900
  end if

  allocate(rPointers(size(triads)+2), STAT=ierr)
  if (ierr /= 0) then
     ierr = allocationError('item group pointers')
     goto 900
  end if

  call readResponsePointers (triads,sup,rPointers,ierr)
  if (ierr < 0) then
     goto 900
  else if (ierr > 0) then
     deallocate(rPointers)
  end if

  call readDisplPointer (sup,disPtr,ierr)
  if (ierr < 0) then
     goto 900
  else if (disPtr > 0_ptr) then
     lDeformation = .false.
     lDumpDef = .false.
  end if

  ! --- Initialize the stress results database

  call initWriteDisp (samSup)

  call writeProgress (' --> Writing result database headers')
  call getFileName ('rdbfile',chname,'.frs')
  if (lDeformation) then
     call writeStressHeader (chname,modelFileName,rdb,samSup,sup,3,ierr)
  else
     call writeStressHeader (chname,modelFileName,rdb,samSup,sup,ierr=ierr)
  end if
  if (ierr /= 0) goto 900

  ! --- Write residual stress state read from external file, if any

  isup = 0
  iStep = 0_i8
  currTime = 0.0_dp
  if (lVTFAt0) then
     isup = vtfFile
     idVTF(2) = idVTF(2) + 2
  end if
  if (addStress0) then

     write(lterm,6010) currTime
     if (iprint > 0 .and. iStep < 100000000_i8) write(lpu,6011) iStep,currTime

     if (rdb%fileNumb >= 0) then
        ! Save time step identification to results database file
        call writeTimeStepDB (rdb,iStep,currTime,ierr)
        if (ierr < 0) goto 900
     end if

     if (lDeformation) then
        sv = 0.0_dp ! Save a zero displacement state to results database file
        call writeDisplacementDB (rdb,samSup,sv,svtot=sv,ierr=ierr)
        if (ierr < 0) goto 900
     end if

     ! Calculate and save residual stresses to results database file
     write(lterm,6030)
     call calcStresses (addStress0,rdb,samSup,isup,idVTF(2),sup%id%baseId, &
          &             iprint=iprint,lpu=lpu,ierr=ierr)
     if (ierr < 0) goto 900

  else if (lVTFAt0) then

     ! Write a dummy scalar block for superelements without residual stresses
     call writeDummyStressVTF (vtfFile,idVTF(2),sup%id%baseId,ierr)
     if (ierr < 0) goto 900

  end if
  if (lVTFAt0) then

     ! Position the result files at t=0.0
     call ffr_setposition (0.0_dp,currTime,iStep)
     if (iStep < 0_i8) then
        ierr = -1
        call reportError (error_p,'Error searching for results at time = 0.0')
        goto 900
     end if

     if (allocated(rPointers)) then
        ! Read superelement position at the initial state
        call readSupElPosition (sup,rPointers(size(triads)+1),iprint,lpu,ierr)
        if (ierr /= 0) goto 900
     end if

     ! Write time step info to VTF-file
     call writeTimeStepVTF (vtfInfo,int(iStep),currTime,ierr)
     if (ierr < 0) goto 900

     ! Write superelement transformations to VTF-file
     idVTF(1) = idVTF(1) + 1
     call writeSupElTransformVTF (vtfFile,idVTF(1),sup,ierr)
     if (ierr < 0) goto 900

  end if

  ! --- Time loop

  currTime = startTime - 1.0_dp
  call writeProgress (' --> Starting time loop')
  do while (ffr_getnextstep(startTime,stopTime,tInc,currTime,iStep))

     write(lterm,6010) currTime
     if (iprint > 0 .and. iStep < 100000000_i8) write(lpu,6011) iStep,currTime

     if (rdb%fileNumb >= 0) then
        ! Save time step identification to results database file
        call writeTimeStepDB (rdb,iStep,currTime,ierr)
        if (ierr < 0) goto 900
     end if

     if (allocated(rPointers)) then
        ! Read system response data for this time step
        call readSupElDisplacements (triads,sup,rPointers,iprint,lpu,ierr)
        if (ierr /= 0) goto 900
     end if

     if (lDumpDef .and. iStep < int(maxInt_p,i8)) then
        ! Write superelement deformations to Nastran bulk-data file
        call writeSupElDeformation (int(iStep),sup,ierr)
        if (ierr /= 0) goto 900
     end if

     if (vtfFile >= 0 .and. iStep < int(maxInt_p,i8)) then
        ! Write time step info to VTF-file
        call writeTimeStepVTF (vtfInfo,int(iStep),currTime,ierr)
        if (ierr < 0) goto 900
        ! Write superelement transformations to VTF-file
        idVTF(1) = idVTF(1) + 1
        call writeSupElTransformVTF (vtfFile,idVTF(1),sup,ierr)
        if (ierr < 0) goto 900
     end if

     ! Read or calculate the internal displacements
     write(lterm,6020)
     if (disPtr > 0_ptr) then
        call readIntDisplacements (disPtr,sup%id,samSup,sv,iprint,lpu,ierr)
     else if (lGrav) then
        call calcIntDisplacements (samSup,sv,sup%finit,sup%genDOFs%ur, &
             &                     iprint,lpu,ierr,vii=vii, &
             &                     g=matmul(grv,sup%supTr(:,1:3)))
     else
        call calcIntDisplacements (samSup,sv,sup%finit,sup%genDOFs%ur, &
             &                     iprint,lpu,ierr)
     end if
     if (ierr < 0) goto 900

     if (lDeformation) then
        ! Calculate total displacement state
        call calcTotalDisplacements (sup,samSup%madof,samSup%minex,sv,sf,ierr)
        if (ierr < 0) goto 900

        ! Save expanded displacement state to results database file
        call writeDisplacementDB (rdb,samSup,sv,svtot=sf,ierr=ierr)
        if (ierr < 0) goto 900

        if (vtfFile >= 0) then
           ! Write expanded displacement state to VTF-file
           idVTF(2) = idVTF(2) + 1
           call writeDisplacementVTF (vtfFile,idVTF(2),sup,samSup,sv,dsVTF,ierr)
           if (ierr < 0) goto 900
           idVTF(2) = idVTF(2) + 1
        end if

     else if (vtfFile >= 0) then
        idVTF(2) = idVTF(2) + 2
     end if

     ! Calculate and save stresses to results database file and VTF
     write(lterm,6030)
     call calcStresses (addStress0,rdb,samSup,vtfFile,idVTF(2),sup%id%baseId, &
          &             sv,sf,iprint=iprint,lpu=lpu,ierr=ierr)
     if (ierr < 0) goto 900

  end do
  if (iStep < -1_i8) then
     ierr = -1
     write(chname,"('time =',1PE12.5)") currTime
     call reportError (error_p,'Error searching for results after '//chname)
     goto 900
  end if

  ! --- Finito, deallocate and close down

  call writeProgress (' --> Time loop done. Closing database files')
  deallocate(triads)
  if (allocated(rPointers)) deallocate(rPointers)
  if (allocated(sv))        deallocate(sv)
  if (allocated(sf))        deallocate(sf)
  if (allocated(vii))       deallocate(vii)
  call deAllocateSAM (samSup)
  if (vtfFile >= 0) then
     if (allocated(saveElmOrder)) deallocate(saveElmOrder)
     call ffa_cmdlinearg_getint ('VTFparts',isup)
     if (isup > 0) then
        ! Write result definition blocks for all parts to the VTF-file
        call ffa_cmdlinearg_getint ('VTFoffset',iprint)
        call writeFooterVTF (vtfFile,vtfInfo,isup,iprint,idVTF(1)-iprint, &
             &               lVTFAt0,lDeformation,'von Mises stress',ierr)
        if (ierr < 0) goto 900
     end if
     ierr = VTFFFileCloseFile(vtfFile) - 1
     call VTFFFileDelete (vtfFile)
     if (ierr < 0) goto 900
  end if
  call closeBandEmatrices (ierr)
  call closeRDBfile (rdb,ierr)
  call ffl_done
  call ffr_done

900 continue
  write(lterm,6000) ' END OF PROGRAM STRESS '
  if (ierr /= 0) then
     call reportError (debugFileOnly_p,'stress')
     call getFileName ('resfile',chname,'_stress.res')
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
6020 format( 11X,'--> .........................Displacements')
6030 format( 11X,'--> .........................Stress results')
6100 format(/11X,'An error is detected. Error messages are written to ',A)
6200 format(/ 4X,'Stress calculation ',A)

end subroutine stress
