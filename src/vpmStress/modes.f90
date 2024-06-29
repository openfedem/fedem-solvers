!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

subroutine modes (nStep,nMode,timeStep,modeNum,processThisMode,ierr)

  use KindModule             , only : dp, i8, maxInt_p, lfnam_p
  use TriadTypeModule        , only : TriadType
  use SupElTypeModule        , only : SupElType
  use SamModule              , only : SamType, writeObject, deallocateSAM
  use SamStressModule        , only : initiateSAM
  use InitiateStressModule   , only : readSolverData
  use DisplacementModule     , only : openBandEmatrices, closeBandEmatrices
  use DisplacementModule     , only : readResponsePointers, calcIntDisplacements
  use DisplacementModule     , only : readSupElDisplacements
  use IdTypeModule           , only : StrId
  use RDBModule              , only : RDBType, openRDBfile, closeRDBfile
  use RDBModule              , only : writeTimeStepDB
  use SaveStressModule       , only : initWriteDisp, writeDisplacementDB
  use SaveStressModule       , only : writeModeHeader, writeModesHeader
  use SaveVTFModule2         , only : writeSupElTransformVTF, writeFooterVTF
  use SaveVTFModule2         , only : writeDisplacementVTF, writetimestepVTF
  use SaveVTFModule2         , only : VTFFFileAppendFile, VTFFFileCloseFile
  use SaveVTFModule2         , only : VTFFFileDelete
  use SaveVTFModule2         , only : VTFFFileSetOutputDebugError
  use SaveVTFModule2         , only : VTFFStateInfoBlockCreate
  use ModesRoutinesModule2   , only : readModesPointers
  use ModesRoutinesModule2   , only : readFrequencies, readSupelModes
  use ModesRoutinesModule2   , only : calcStrainEnergyDensity
  use FileUtilitiesModule    , only : getFileName
  use TimerModule            , only : initTime, showTime
  use VersionModule          , only : openResFile
  use ProgressModule         , only : lterm, writeProgress
  use PointerKindModule      , only : ptr
  use ReportErrorModule      , only : allocationError, reportError
  use ReportErrorModule      , only : openTerminalOutputFile
  use ReportErrorModule      , only : note_p, error_p, debugFileOnly_p
  use BinaryDBInterface      , only : readDoubleDB
  use FFaCmdLineArgInterface , only : ffa_cmdlinearg_getint
  use FFaCmdLineArgInterface , only : ffa_cmdlinearg_getbool
  use FFaCmdLineArgInterface , only : ffa_cmdlinearg_isTrue
  use FFaCmdLineArgInterface , only : ffa_cmdlinearg_isSet
  use FFaCmdLineArgInterface , only : ffa_cmdlinearg_intValue
  use FFlLinkHandlerInterface, only : ffl_init, ffl_done, ffl_export_vtf
  use FFrExtractorInterface  , only : ffr_init, ffr_done, ffr_setposition

  implicit none

  integer , intent(in)       :: nStep, nMode
  real(dp), intent(in)       :: timeStep(nStep)
  integer , intent(in)       :: modeNum(nMode)
  logical , intent(inout)    :: processThisMode(nStep,nMode)
  integer , intent(out)      :: ierr
  integer                    :: lpu, iprint, i, j, isup, istep, lerr
  integer                    :: iComp, nedf, iExp, vtfInfo, idVTF(2)
  integer(i8)                :: gotStep
  integer      , allocatable :: vtfIds(:,:), vtfFile(:)
  integer(ptr) , allocatable :: rPointers(:), mPointers(:,:)
  logical                    :: lComplex, lEnergyDensity, lGrav, multiFiles
  real(dp)                   :: gotTime, dsVTF, grv(3)
  real(dp)     , allocatable :: sv(:), vii(:), freq(:), eigFinit(:,:,:)
  type(RDBType), allocatable :: rdb(:)
  type(SamType)              :: samSup
  type(SupElType)            :: sup
  type(TriadType), pointer   :: triads(:)
  character(len=lfnam_p)     :: chname, modelFileName
  character(len=64)          :: errMsg

  !! --- Logic section ---

  call initTime ()
  call openTerminalOutputFile (lterm)
  write(lterm,6000) 'START OF PROGRAM MODES'

  ! Open the modes results file and write heading
  call getFileName ('resfile',chname,'_modes.res')
  call openResFile (chname,'Modal Recovery',lpu,ierr)
  if (ierr /= 0) return

  call ffa_cmdlinearg_getint ('debug',iprint)
  call ffa_cmdlinearg_getint ('linkId',isup)
  call ffa_cmdlinearg_getbool ('damped',lComplex)
  call ffa_cmdlinearg_getbool ('energy_density',lEnergyDensity)
  if (lComplex) then
     iComp = 2
  else
     iComp = 1
  end if

  ! --- Read the link file

  call writeProgress (' --> Reading link files')
  call ffl_init (.false.,ierr)
  if (ierr /= 0) then
     call reportError (error_p,'Error initializing the link object')
     goto 900
  end if

  ! --- Establish the SAM datastructure

  call getFileName ('samfile',chname,'.fsm')
  call initiateSAM (chname,samSup,lpu,ierr)
  if (ierr /= 0) goto 900

  ! --- Read superelement data from the solver input file

  if (isup < 1) isup = samSup%mpar(18) ! Fetch baseID from reducer if not given
  write(lterm,"(/11X,'--> Process Part; baseID (isup) =',I6)") isup
  call readSolverData (isup,grv,sup,triads,modelFileName,iprint,lpu,ierr)
  if (ierr /= 0) goto 900

  if (iprint > 4) then
     call WriteObject (samSup,lpu,3)
     write(lpu,"(/)")
  end if

  ! --- Open the B-matrix and the generalized modes files

  call openBandEmatrices (samSup%ndof1,samSup%ndof2, &
       &                  samSup%mpar(22),samSup%mpar(31),ierr)
  if (ierr < 0) goto 900

  allocate(sv(samSup%ndof*iComp),eigFinit(sup%nTotDofs,iComp,nMode), STAT=ierr)
  if (ierr /= 0) then
     ierr = allocationError('displacement- and eigenvectors')
     goto 900
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

  allocate(vtfFile(nMode), STAT=ierr)
  if (ierr /= 0) then
     ierr = allocationError('VTF-file handles')
     goto 900
  end if

  vtfFile = -1
  vtfInfo = -1
  iExp    = 0
  chname  = ''

  if (ffa_cmdlinearg_isSet('VTFfile')) then
     call ffa_cmdlinearg_getstring ('VTFfile',chname)
  end if

  if (chname /= '') then
     if (ffa_cmdlinearg_isTrue('VTFexpress')) then
        iExp = index(chname,'.',.TRUE.)
        if (iExp < 1) iExp = len_trim(chname)
     end if
     if (iExp > 0 .and. nStep > 1) then
        call reportError (note_p, &
             'When writing express VTF-files, only the mode shapes', &
             'specified for the first time step are exported')
     end if
     if (lComplex) then
        call reportError (note_p, &
             'VTF output of complex mode shapes is not supported.', &
             'Only the imaginary part of the mode shape vector is exported.')
     end if

     do i = 1, nMode
        if (iExp > 0 .and. nMode > 1) then
           ! Export each mode shape to separate VTF-files
           chname(iExp:) = '_'//trim(adjustl(StrId(modeNum(i))))//'.vtf'
        else if (i > 1) then
           ! Use the same file for all mode shapes
           vtfFile(i) = vtfFile(i-1)
           cycle
        end if

        call writeProgress (' --> Writing part geometry to VTF-file')
        call ffl_export_vtf (trim(chname),trim(sup%id%descr),sup%id%baseId,ierr)
        if (ierr /= 0) then
           call reportError (error_p,'Error writing to VTF-file '//chname)
           goto 900
        end if

        ! Reopen the VTF-file for fortran output
        call VTFFFileSetOutputDebugError (1)
        vtfFile(i) = VTFFFileAppendFile(trim(chname))
        if (vtfFile(i) < 0) then
           ierr = vtfFile(i)
           call reportError (error_p,'Error re-opening VTF-file '//chname)
           goto 900
        end if

     end do

     ! Initialize VTF result block id numbers
     idVTF = ffa_cmdlinearg_intValue('VTFoffset')
     if (ffa_cmdlinearg_intValue('VTFparts') > 0) then
        allocate(vtfIds(nStep,nMode),freq(nMode), STAT=ierr)
        if (ierr /= 0) then
           ierr = allocationError('VTF block id-numbers')
           goto 900
        end if
        vtfIds = 0
        freq = 0.0_dp
        if (iExp == 0) then
           vtfInfo = VTFFStateInfoBlockCreate()
        end if
     end if

     call ffa_cmdlinearg_getdouble ('VTFdscale',dsVTF)

  end if

  ! --- Open the solver results database

  call writeProgress (' --> Reading solver result files')
  call ffr_init ('frsfile',ierr)
  if (ierr /= 0) then
     call reportError (error_p,'Error opening results database')
     goto 900
  end if

  allocate(rPointers(size(triads)+2), &
       &   mPointers(size(triads)+1,nMode), STAT=ierr)
  if (ierr /= 0) then
     ierr = allocationError('item group pointers')
     goto 900
  end if

  call readResponsePointers (triads,sup,rPointers,ierr)
  if (ierr /= 0) goto 900

  call readModesPointers (sup,nMode,modeNum,processThisMode,mPointers,ierr)
  if (ierr >= nMode) then
     call reportError (error_p,'Could not find any eigenmodes to expand')
     goto 900
  else if (ierr < 0) then
     goto 900
  end if

  ! --- Initialize the modes results database

  call initWriteDisp (samSup,multiFiles)
  if (.not. multiFiles) then
     multiFiles = lEnergyDensity .or. .not.all(processThisMode)
  end if

  if (multiFiles) then
     allocate(rdb(nMode+1), STAT=ierr)
  else
     allocate(rdb(1), STAT=ierr)
  end if
  if (ierr /= 0) then
     ierr = allocationError('result database files')
     goto 900
  end if

  call writeProgress (' --> Writing result database headers')
  call getFileName ('rdbfile',chname,'.frs')

  if (multiFiles) then

     ! File for the dynamic response
     call writeModeHeader (chname,modelFileName,0,.false., &
          &                rdb(1),samSup,sup,ierr)
     if (ierr /= 0) goto 900
     call closeRDBfile (rdb(1),ierr,.true.)
     if (ierr < 0) goto 900

     ! One additional file for each mode
     do i = 1, nMode
        if (any(processThisMode(:,i))) then
           call writeModeHeader (chname,modelFileName,modeNum(i),lComplex, &
                &                rdb(i+1),samSup,sup,ierr)
           if (ierr /= 0) goto 900
           call closeRDBfile (rdb(i+1),ierr,.true.)
           if (ierr < 0) goto 900
        else
           write(errMsg,"('Eigenmode',i3,' will not be expanded')") modeNum(i)
           call reportError (note_p,errMsg)
        end if
     end do

  else

     ! Write all mode shapes to the same file
     call writeModesHeader (chname,modelFileName,modeNum,lComplex, &
          &                 rdb(1),samSup,sup,ierr)
     if (ierr /= 0) goto 900

  end if

  ! --- Time step loop

  lerr = 0
  call writeProgress (' --> Starting time loop')
  do iStep = 1, nStep

     ! Position the result files at current time step
     call ffr_setposition (timeStep(iStep),gotTime,gotStep)
     if (gotStep < 0_i8) then
        lerr = lerr + 1
        write(chname,"('time =',1PE12.5)") timeStep(iStep)
        call reportError (error_p,'Error searching for results at '//chname)
        cycle
     end if
     write(lterm,6010) gotTime
     write(lpu,6011) gotStep,gotTime

     ! Read system response data for this time step
     call readSupElDisplacements (triads,sup,rPointers,iprint,lpu,ierr)
     if (ierr /= 0) goto 900

     if (iStep == 1 .and. allocated(freq)) then
        ! Read system eigenfrequencies for this time step
        call readFrequencies (nMode,modeNum,processThisMode(iStep,:),freq,ierr)
        if (ierr /= 0) goto 900
     end if

     ! Read system modes data for this time step
     call readSupElModes (iComp,nMode,processThisMode(iStep,:),mPointers,sup, &
          &               sv,eigFinit,ierr)
     if (ierr /= 0) goto 900

     idVTF(1) = idVTF(1) + 1
     if (vtfFile(1) >= 0 .and. gotStep <= int(maxInt_p,i8)) then
        ! Write time step info to VTF-file
        call writeTimeStepVTF (vtfInfo,int(gotStep),gotTime,ierr)
        if (ierr < 0) goto 900
        ! Write superelement transformations to VTF-file(s)
        call writeSupElTransformVTF (vtfFile(1),idVTF(1),sup,ierr)
        if (ierr < 0) goto 900
        if (iExp > 0 .and. nMode > 1) then
           do i = 2, nMode
              call writeSupElTransformVTF (vtfFile(i),idVTF(1),sup,ierr)
              if (ierr < 0) goto 900
           end do
        end if
     end if

     ! Calculate internal displacements
     write(lterm,6020)
     if (lGrav) then
        call calcIntDisplacements (samSup,sv,sup%finit,sup%genDOFs%ur, &
             &                     iprint,lpu,ierr,vii=vii, &
             &                     g=matmul(grv,sup%supTr(:,1:3)))
     else
        call calcIntDisplacements (samSup,sv,sup%finit,sup%genDOFs%ur, &
             &                     iprint,lpu,ierr)
     end if
     if (ierr < 0) goto 900

     ! Save time step identification and the expanded displacement state
     ! to the deformation results database file
     if (multiFiles) then
        call openRDBfile (rdb(1),ierr=ierr)
        if (ierr < 0) goto 900
     end if
     call writeTimeStepDB (rdb(1),gotStep,gotTime,ierr)
     if (ierr < 0) goto 900
     call writeDisplacementDB (rdb(1),samSup,sv,ierr)
     if (ierr < 0) goto 900
     if (multiFiles) then
        call closeRDBfile (rdb(1),ierr)
        if (ierr < 0) goto 900
     end if

     idVTF(2) = idVTF(2) + 1
     if (vtfFile(1) >= 0 .and. iExp == 0) then
        ! Write expanded displacement state to VTF-file
        call writeDisplacementVTF (vtfFile(1),idVTF(2),sup,samSup,sv,dsVTF,ierr)
        if (ierr < 0) goto 900
     end if

     ! --- Eigenmode loop

     nedf = sup%nTotDofs - sup%genDOFs%nDOFs
     write(lterm,6030,advance='NO')
     do i = 1, nMode
        if (.not. processThisMode(iStep,i)) cycle

        write(lterm,"(I4)",advance='NO') modeNum(i)
        if (iprint > 0) write(lpu,6031) i
        do j = 1, iComp

           ! Calculate internal displacements corresponding to this eigenmode
           call calcIntDisplacements (samSup,sv(1+samSup%ndof*(j-1)), &
                &                     eigFinit(:,j,i), &
                &                     eigFinit(nedf+1:sup%nTotDofs,j,i), &
                &                     iprint,lpu,ierr)
           if (ierr < 0) goto 900

        end do

        if (multiFiles) then

           ! Save time step identification and the expanded mode shape
           ! to the results database file associated with current mode
           call openRDBfile (rdb(i+1),ierr=ierr)
           if (ierr < 0) goto 900
           call writeTimeStepDB (rdb(i+1),gotStep,gotTime,ierr)
           if (ierr < 0) goto 900
           call writeDisplacementDB (rdb(i+1),samSup,sv,ierr,iComp)
           if (ierr < 0) goto 900

           if (lEnergyDensity) then
              ! Calculate strain energy density for this mode and save to rdb
              call calcStrainEnergyDensity (rdb(i+1),samSup,sv,iprint,lpu,ierr)
              if (ierr < 0) goto 900
           end if

           call closeRDBfile (rdb(i+1),ierr)
           if (ierr < 0) goto 900

        else

           ! Save the expanded mode shape to the results database file
           call writeDisplacementDB (rdb(1),samSup,sv,ierr,iComp)
           if (ierr < 0) goto 900

        end if

        if (vtfFile(i) >= 0) then
           ! Write the expanded mode shape to VTF-file
           if (iExp == 0) idVTF(2) = idVTF(2) + 1
           if (allocated(vtfIds)) vtfIds(iStep,i) = idVTF(2)
           call writeDisplacementVTF (vtfFile(i),idVTF(2),sup,samSup, &
                &                     sv(1+samSup%ndof*(iComp-1):),dsVTF,ierr)
           if (ierr < 0) goto 900
        end if

     end do

     if (iExp > 0) exit ! Only write one time step when express files
  end do

  ! --- Finito, deallocate and close down

  call closeRDBfile (rdb(1),ierr)
  if (ierr < 0) goto 900

  call writeProgress (' --> Done time loop. Closing database files')
  deallocate(triads,sv,eigFinit,rPointers,mPointers,rdb)
  if (lGrav) deallocate(vii)
  call deAllocateSAM (samSup)
  if (vtfFile(1) >= 0) then
     if (allocated(vtfIds) .and. lerr == 0) then
        ! Write result definition blocks for all parts to the VTF-file
        call ffa_cmdlinearg_getint ('VTFoffset',istep)
        call ffa_cmdlinearg_getint ('VTFparts',isup)
        if (iExp > 0) then
           do i = 1, nMode
              call writeFooterVTF (vtfFile(i),istep,isup, &
                   &               modeNum(i),vtfIds(1,i),freq(i),ierr)
              if (ierr < 0) goto 900
           end do
        else
           call writeFooterVTF (vtfFile(1),vtfInfo,istep,isup,nStep,nMode, &
                &               modeNum,vtfIds,freq,ierr)
           if (ierr < 0) goto 900
        end if
     end if
     if (allocated(vtfIds)) deallocate(vtfIds)
     if (allocated(freq)) deallocate(freq)
     if (iExp > 0) then
        do i = 1, nMode
           ierr = VTFFFileCloseFile(vtfFile(i)) - 1
           call VTFFFileDelete (vtfFile(i))
           if (ierr < 0) goto 900
        end do
     else
        ierr = VTFFFileCloseFile(vtfFile(1)) - 1
        call VTFFFileDelete (vtfFile(1))
        if (ierr < 0) goto 900
     end if
  end if
  call closeBandEmatrices (ierr)
  call ffl_done
  call ffr_done
  if (ierr >= 0) ierr = lerr

900 continue
  write(lterm,6000) ' END OF PROGRAM MODES '
  if (ierr /= 0) then
     call reportError (debugFileOnly_p,'modes')
     call getFileName ('resfile',chname,'_modes.res')
     write(lterm,6100) chname
     call showTime (lpu)
     write(lpu,6200) 'failed :-('
  else
     call showTime (lpu)
     write(lpu,6200) 'successfully completed :-)'
  end if
  call closeLogFiles ()

6000 format(//11X,16('='),'> ',A,' <',16('='))
6010 format(//11X,'--> ......Simulation time :',1PD13.5)
6011 format(/120('*')//5X,'TIME STEP',I8,' : INTEGRATION TIME T =',1PD13.5)
6020 format( 11X,'--> .........................Displacements')
6030 format( 11X,'--> .........................Modes results')
6031 format(//5X,9('='),'> MODE',I3,' <',75('=') )
6100 format(/11X,'An error is detected. Error messages are written to ',A)
6200 format(/ 4X,'Modal expansion ',A)

end subroutine modes
