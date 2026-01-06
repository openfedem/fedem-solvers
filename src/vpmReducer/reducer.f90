!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file reducer.f90
!> @brief Main driver for the FEDEM FE part reducer.

!!==============================================================================
!> @brief Main driver for the FEDEM FE part reducer.
!>
!> @param[out] ierr Error flag
!>
!> @callgraph
!>
!> @author Knut Morten Okstad
!>
!> @date 15 Oct 2000

subroutine reducer (ierr)

  use KindModule                , only : dp, sp, i8, nbs_p, nbd_p, lfnam_p
  use KindModule                , only : epsDiv0_p, pi_p
  use SamModule                 , only : SamType, deallocateSAM
  use SamModule                 , only : NodeNumFromEqNum
  use SamReducerModule          , only : saveSAM
  use SaveReducerModule         , only : writeModes
  use SysMatrixTypeModule       , only : SysMatrixType, outOfCore_p
  use SysMatrixTypeModule       , only : sparseMatrix_p, diagonalMatrix_p
  use SysMatrixTypeModule       , only : saveSysMat, restoreSysMat
  use SysMatrixTypeModule       , only : deallocateSysMatrix
  use SparseMatrixModule        , only : SparseMatrixType, smDeallocate
  use DiskMatrixModule          , only : DiskMatrixType, dmSize, dmGetSwapSize
  use DiskMatrixModule          , only : dmNullify, dmOpen, dmClose
  use DiskMatrixModule          , only : dmFindMinConnections
  use AsmExtensionModule        , only : csBeginAssembly, csEndAssembly
  use AsmExtensionModule        , only : csRepermute, castToInt8
  use InputReducerModule        , only : readReducerData
  use InaddModule               , only : INADD, extractSubMat
  use CmstrsModule              , only : CMSTRS, EIGVAL, EIGCMS, JCMS, GRAV
  use FileUtilitiesModule       , only : getFileName
  use TimerModule               , only : initTime, showTime
  use VersionModule             , only : openResFile
  use ProgressModule            , only : lterm, writeProgress
  use AllocationModule          , only : doLogMem, logAllocMem, reAllocate
  use ScratchArrayModule        , only : releaseScratchArrays
  use ReportErrorModule         , only : empty_p, note_p, warning_p, error_p
  use ReportErrorModule         , only : debugFileOnly_p, openTerminalOutputFile
  use ReportErrorModule         , only : allocationError, reportError
  use BinaryDBInterface         , only : writeDoubleDB, closeBinaryDB
  use FFaProfilerInterface      , only : ffa_initProfiler, ffa_reportTimer
  use FFaProfilerInterface      , only : ffa_startTimer, ffa_stopTimer
  use FFaProfilerInterface      , only : ffa_printMemStatus
  use FFaCmdLineArgInterface    , only : ffa_cmdlinearg_getint
  use FFaCmdLineArgInterface    , only : ffa_cmdlinearg_getbool
  use FFaCmdLineArgInterface    , only : ffa_cmdlinearg_getdouble
  use FFaCmdLineArgInterface    , only : ffa_cmdlinearg_isSet
  use FFlLinkHandlerInterface   , only : ffl_done, ffl_calcs, ffl_getcs
  use FFlLinkHandlerInterface   , only : ffl_addcs_int, ffl_addcs_double

  implicit none

  integer , intent(out)  :: ierr
  logical                :: dataChk, noGrav, noRedMass, diagMass, lumpedMass
  logical                :: factorMass, calcGdisp, twoLoops
  integer                :: i, imass, istiff, iopSing, lpu, iprint
  integer                :: cs, ndim, nenod, neval, ngen, nevred, nrhs, mlc(100)
  integer                :: iNod, lDof, lowBconn, Bprec
  integer(i8)            :: nMass, nStiff
  real(dp), parameter    :: tolMass_p = 1.0e-6_dp
  real(sp)               :: autoBramRatio
  real(dp)               :: sMass, rMass, tolEigval, tolFactorize, eigenShift
  real(dp)               :: minAval(100)
  real(dp), allocatable  :: smee(:,:), skee(:,:), rhs(:,:), vgi(:,:)
  real(dp), allocatable  :: rval(:), rvec(:,:), gravec(:,:), sm(:,:), sk(:,:)
  real(sp), allocatable  :: work(:)
  real(dp), pointer      :: tenc(:,:), eval(:), evec(:,:), csRhs(:,:)
  real(dp), pointer      :: tmpVec(:,:), tmpRhs(:,:)
  real(dp), target       :: dummy(1,1)
  integer , pointer      :: medof(:)
  type(SamType)          :: sam
  type(SysMatrixType)    :: csStiffMat, csMassMat
  type(SparseMatrixType) :: smieSparse, skieSparse
  type(DiskMatrixType)   :: BmatDisk
  character(len=lfnam_p) :: chName
  character(len=64)      :: errMsg


  !! --- Logic section ---

  call initTime ()
  call ffa_initprofiler ('fedem_reducer profiler')
  call ffa_starttimer ('fedem_reducer')
  call openTerminalOutputFile (lterm)
  write(lterm,6000) 'START OF PROGRAM REDUCER'
  call ffa_printMemStatus (lterm)

  ierr = 0
  iMass = -1
  iStiff = -1
  nMass = 0_i8
  nStiff = 0_i8
  nullify(tenc)
  nullify(eval)
  nullify(evec)
  nullify(csRhs)
  nullify(medof)
  call dmNullify (BmatDisk)

  ! Open the reducer result file and write heading
  call getFileName ('resfile',chname,'.res')
  call openResFile (chname,'FE model Reducer',lpu,ierr)
  if (ierr /= 0) return

  call ffa_cmdlinearg_getint ('debug',iprint)
  call ffa_cmdlinearg_getint ('singularityHandler',iopSing)
  call ffa_cmdlinearg_getint ('Bmatprecision',Bprec)
  call ffa_cmdlinearg_getbool ('datacheck',dataChk)
  call ffa_cmdlinearg_getbool ('twoAssemblyLoops',twoLoops)
  call ffa_cmdlinearg_getbool ('nograv',noGrav)
  call ffa_cmdlinearg_getbool ('nomass',noRedMass)
  call ffa_cmdlinearg_getbool ('diagmass',diagMass)
  call ffa_cmdlinearg_getbool ('lumpedmass',lumpedMass)
  call ffa_cmdlinearg_getbool ('factorMass',factorMass)
  call ffa_cmdlinearg_getdouble ('tolEigval',tolEigval)
  call ffa_cmdlinearg_getdouble ('tolFactorize',tolFactorize)
  call ffa_cmdlinearg_getdouble ('eigenshift',eigenShift)


  ! --- Read the link file and establish the SAM datastructure

  call writeProgress (' --> Reading FE data file')
  call readReducerData (sam, csStiffMat, csMassMat, sMass, rMass, &
       &                tenc, medof, mlc, noGrav, diagMass, lumpedMass, &
       &                .false., iprint, lpu, ierr)
  if (dataChk .or. ierr /= 0) then
     call reAllocate ('reducer',tenc)
     call reAllocate ('reducer',medof)
     if (ierr /= 0 .or. (iprint /= -7 .and. iprint /= -6)) goto 900
  else if (diagMass .and. sam%nmmceq > 0) then
     call reportError (warning_p,'This FE part has constraint equations.', &
          'Forcing a diagonal mass matrix may yield poor convergence', &
          'and/or inaccurate results in the dynamics simulation.', &
          'Switching to lumped- or consistent element mass is recommended.')
  end if
  call ffa_printMemStatus (lterm)

  call ffa_cmdlinearg_getint ('printArray',sam%mpar(25))
  if (sam%mpar(25) < 100) then
     sam%mpar(26) = sam%mpar(25)
  else
     sam%mpar(26) = mod(sam%mpar(25),100)
     sam%mpar(25) =     sam%mpar(25)/100
  end if

  ! --- Calculate the link checksum

  call ffl_calcs (cs)
  write(lterm,1) cs
1 format(7X,'FE part checksum :',I12)
  call ffl_addcs_int (sam%mpar(22))           ! ngen
  call ffl_addcs_int (sam%mpar(19))           ! neval
  call ffl_addcs_int (2)                      ! 1=skyline, 2=sparse
  call ffl_addcs_double (tolEigval)
  call ffl_addcs_double (tolFactorize)
  if (diagMass) then
     call ffl_addcs_int (1)
  else if (lumpedMass) then
     call ffl_addcs_int (2)
  end if
  if (.not.factorMass .and. sam%mpar(22) > 0) call ffl_addcs_int (1)
  call ffl_getcs (cs,ierr)
  if (ierr /= 0) goto 890
  write(lterm,1) cs

  if (tolFactorize < 1.0e-20_dp) then
     tolFactorize = 1.0e-20_dp
     call reportError (note_p,'A factorization tolerance lower than 1e-20'// &
          ' is not allowed.','The factorization tolerance is reset to 1e-20')
  end if


  ! --- Assemble the system load vectors, if any

  if (ffa_cmdlinearg_isSet('loadfile')) then
     nrhs = sam%mpar(21)
  else
     nrhs = 0 ! Don't compute reduced load vectors
  end if
  if (nrhs > 0) then

     allocate(csRhs(sam%neq,nrhs),STAT=ierr)
     if (ierr /= 0) then
        ierr = allocationError('reducer: System right-hand-side vectors')
        goto 890
     else if (doLogMem) then
        call logAllocMem ('reducer',0,size(csRhs),nbd_p)
     end if

     do i = 1, nrhs
        call INADD (sam, csRhs(:,i), mlc(i), iprint, lpu, ierr)
        if (ierr < 0) goto 890
     end do

     tmpRhs => csRhs
  else
     tmpRhs => dummy
     dummy = 0.0_dp
  end if

  if (twoLoops .or. noGrav) then

     ! --- Assemble stiffness and mass matrices separately, stiffness first

     if (.not.dataChk) then
        call writeProgress (' --> Allocating substructure stiffness matrix')
        call csBeginAssembly (csStiffMat,'system stiffness matrix',lpu,ierr)
        if (ierr /= 0) goto 890
        call ffa_printMemStatus (lterm)
     end if

     ! --- Now assemble the substructure stiffness matrix

     call INADD (sam, csStiffMat, csMassMat, sMass, rMass, &
          &      diagMass.or.lumpedMass, dataChk, 1, iprint, lpu, ierr)
     if (ierr < 0) then
        call reportError (error_p,'Can not build substructure stiffness matrix')
        goto 890
     end if

     if (.not.dataChk) then
        call csEndAssembly (csStiffMat,ierr)
        if (ierr /= 0) goto 890
        call ffa_printMemStatus (lterm)
     end if

     if (noGrav) then
        call ffl_done ! Assembly finished, release the link object
        goto 2 ! Skip gravity force calculation and mass matrix reduction
     end if

  end if
  if (twoLoops) then

     ! --- Then do the mass matrix

     if (.not.dataChk) then
        call writeProgress (' --> Allocating substructure mass matrix')
        call csBeginAssembly (csMassMat,'system mass matrix',lpu,ierr)
        if (ierr /= 0) goto 890
        call ffa_printMemStatus (lterm)
     end if

     ! --- Now assemble the substructure mass matrix

     call INADD (sam, csStiffMat, csMassMat, sMass, rMass, &
          &      diagMass.or.lumpedMass, dataChk, 2, iprint, lpu, ierr)
     if (ierr < 0) then
        call reportError (error_p,'Can not build substructure mass matrix')
        goto 890
     end if

     call ffl_done ! Assembly finished, release the link object
     if (dataChk) goto 900

     call csEndAssembly (csMassMat,ierr)
     if (ierr /= 0) goto 900

  else ! Assemble the stiffness and mass in the same loop

     ! --- Allocate substructure stiffness- and mass matrices

     if (.not.dataChk) then
        call writeProgress (' --> Allocating substructure matrices')
        call csBeginAssembly (csStiffMat,'system stiffness matrix',lpu,ierr)
        if (ierr /= 0) goto 890
        call csBeginAssembly (csMassMat,'system mass matrix',lpu,ierr)
        if (ierr /= 0) goto 890
        call ffa_printMemStatus (lterm)
     end if

     ! --- Now assemble the substructure mass- and stiffness matrices

     call INADD (sam, csStiffMat, csMassMat, sMass, rMass, &
          &      diagMass.or.lumpedMass, dataChk, 3, iprint, lpu, ierr)
     if (ierr < 0) then
        call reportError (error_p,'Can not build substructure matrices')
        goto 890
     end if

     call ffl_done ! Assembly finished, release the link object
     if (dataChk) goto 900

     call csEndAssembly (csStiffMat,ierr)
     if (ierr /= 0) goto 900
     call csEndAssembly (csMassMat,ierr)
     if (ierr /= 0) goto 900

  end if

2 continue
  call ffa_printMemStatus (lterm)

  ! --- Store the substructure mass- and stiffness matrices on temporary files

  calcGdisp = ffa_cmdlinearg_isSet('dispfile')
  if (calcGdisp .or. sam%mpar(19) > 0) then
     if (csStiffMat%storageType /= outOfCore_p) then
        nStiff = size(csStiffMat%value,kind=i8)
     end if
  end if
  if (sam%mpar(19) > 0) then
     if (csMassMat%storageType /= outOfCore_p .and. factorMass) then
        nMass = size(csMassMat%value,kind=i8)
     end if
     autoBramRatio = 0.5_sp
  else
     autoBramRatio = 0.9_sp
  end if

  if (nStiff > 0_i8) then
     call writeProgress ('     Saving substructure matrices on temporary files')
     call saveSysMat (csStiffMat,'stiffness matrix',nStiff,iStiff,ierr)
     if (ierr < 0) goto 900
     if (nMass > 0_i8) call writeProgress (1,2)
  end if

  if (nMass > 0_i8) then
     call saveSysMat (csMassMat,'mass matrix',nMass,iMass,ierr)
     if (ierr < 0) goto 900
     if (nStiff > 0_i8) call writeProgress (2,2)
  end if


  ! --- Allocate submatrices ie and ee

  if (noGrav) noRedMass = .true. ! No gravity also implies no mass matrix
  if (noRedMass) then
     ndim = 0
     i = 0
  else if (csMassMat%storageType == diagonalMatrix_p) then
     ndim = sam%ndof2
     i = 1
  else
     ndim = sam%ndof2
     i = ndim
  end if
  allocate(smee(ndim,i),skee(sam%ndof2,sam%ndof2),STAT=ierr)
  if (ierr /= 0) then
     ierr = allocationError('reducer: Submatrices smee and skee')
     goto 900
  else if (doLogMem) then
     call logAllocMem ('reducer',0,size(smee)+size(skee),nbd_p)
  end if

  ! --- Extract submatrices ie and ee from the substructure matrices

  call extractSubMat (sam, csStiffMat, csMassMat, &
       &              smieSparse, smee, skieSparse, skee, iprint, lpu, ierr)
  if (ierr /= 0) goto 900
  call ffa_printMemStatus (lterm)

  ! --- Allocate disk swap matrix

  write(lpu,*)
  call getFileName ('Bmatfile',chname,'_B.fmx')
  call dmOpen (BmatDisk,chname,cs,sam%ndof1,sam%ndof2,Bprec,ierr, &
       &       swapSize=dmGetSwapSize(sam%ndof1,sam%ndof2,autoBramRatio))
  if (ierr < 0) goto 900
  call ffa_printMemStatus (lterm)

  sam%mpar(31) = int(dmSize(BmatDisk,3)) ! Number of columns pr swap section

  ! --- Store all SAM-arrays on file for fedem_stress

  call getFileName ('samfile',chname,'_SAM.fsm')
  call saveSAM (chname,cs,sam,ierr)
  if (ierr < 0) then
     call reportError (error_p,'Can not save SAM-data to file '//chname)
     goto 900
  end if


  neval = sam%mpar(19)
  ngen  = sam%mpar(22)
  ndim  = sam%mpar(24)
  nenod = sam%mpar(29)

  if (neval > 0) then

     ! --- Compute eigenvalues and eigenvectors of the substructure

     call EIGVAL (sam, csStiffMat, csMassMat, neval, factorMass, &
          &       iopSing, tolFactorize, tolEigval, eigenShift, &
          &       eval, evec, iprint, lpu, ierr)
     if (ierr < 0) then
        call reportError (error_p,'Can not perform eigenvalue calculation')
        goto 900
     end if
     call ffa_printMemStatus (lterm)

     write(lpu,6010) (i,eval(i)*0.5_dp/pi_p,i=1,neval)

     call getFileName ('eigfile',chname,'_E.fmx')
     if (Bprec == 1) then
        ! Allocate work array for single precision casting
        allocate(work(size(evec,1)),STAT=ierr)
        if (ierr /= 0) then
           ierr = allocationError('reducer: SP work array')
           goto 900
        else if (doLogMem) then
           call logAllocMem ('reducer',0,size(work),nbs_p)
        end if
        call writeDoubleDB (chname,'generalized modes',cs,evec,ierr,work)
        if (doLogMem) then
           call logAllocMem ('reducer',size(work),0,nbs_p)
        end if
        deallocate(work)
     else
        call writeDoubleDB (chname,'generalized modes',cs,evec,ierr)
     end if
     if (ierr /= 0) then
        call reportError (error_p,'Can not write eigenvectors to file')
        goto 900
     end if

     call getFileName ('frsfile',chname)
     if (chname /= '') then
        ! Save the component mode shapes to frs-file for plotting
        call writeModes (chname,sam,sam%mpar(18),'Component modes', &
             &           neval,eval,evec,ierr)
        if (ierr /= 0) then
           call reportError (warning_p,'Failed to write mode shapes file', &
                &                      'Execution continues...')
        end if
     end if

     if (csStiffMat%storageType == sparseMatrix_p) then
        ! Restore in-core sparse matrix only if shift has been applied
        ! Otherwise, re-use the factorization already computed
        if (csStiffMat%sparse%mspar(1) < 9) nStiff = -nStiff
     end if

     if (nStiff > 0_i8) then
        ! Restore the un-factorized stiffness matrix from temporary file
        call restoreSysMat (csStiffMat,'stiffness matrix',nStiff,iStiff,ierr)
        if (ierr < 0) goto 900
     end if

     if (nMass > 0_i8) then
        ! Restore the un-factorized mass matrix from temporary file
        call restoreSysMat (csMassMat,'mass matrix',-nMass,iMass,ierr)
        if (ierr < 0) goto 900
     end if

     if ( csStiffMat%storageType == outOfCore_p .and. &
          csMassMat%storageType == sparseMatrix_p ) then
        ! Recompute the internal permutation arrays for the SPR mass matrix
        ! now based on the external equation ordering
        call csRepermute (csMassMat, sam%meqn, lpu, ierr)
        if (ierr < 0) goto 900
     end if

     tmpVec => evec
  else
     tmpVec => dummy
     dummy = 0.0_dp
  end if


  ! --- Allocate reduced mass/stiffness matrices, and load vectors

  if (noRedMass) then
     i = 0
  else
     i = ndim
  end if
  allocate(sm(i,i),sk(ndim,ndim),rhs(ndim,nrhs),STAT=ierr)
  if (ierr == 0 .and. calcGdisp) then
     allocate(vgi(sam%ndof1,3),STAT=ierr)
  else
     allocate(vgi(0,0))
  end if
  if (ierr /= 0) then
     ierr = allocationError('reducer: Reduced substructure matrices')
     goto 900
  else if (doLogMem) then
     call logAllocMem ('reducer',0,size(sm)+size(sk)+size(rhs)+size(vgi),nbd_p)
  end if

  ! --- Perform the CMS-transformation

  call CMSTRS (iopSing, sam, BmatDisk, &
       &       csMassMat, smieSparse, smee, sm, csStiffMat, skee, sk, &
       &       tmpVec, tmpRhs, rhs, vgi, tolFactorize, &
       &       iStiff, ngen, nrhs, lpu, ierr)
  if (ierr /= 0) then
     call reportError (error_p,'Can not perform CMS-transformation')
     goto 900
  end if
  call ffa_printMemStatus (lterm)

  call ffa_cmdlinearg_getint ('printLowBmatConnections',lowBconn)
  if (lowBconn > 0) then
     call dmFindMinConnections (BmatDisk,mlc(1:lowBconn),minAval,ierr)
     if (ierr < 0) goto 900
     write(lpu,6300)
     do i = 1, lowBconn
        lDof = mlc(i)
        iNod = NodeNumFromEqNum(sam,lDof,.true.)
        write(lpu,6310) iNod,lDof,minAval(mlc(i))
     end do
  end if

  if (nStiff /= 0_i8) then
     ! Close (and delete) the temporary stiffness matrix file
     call closeBinaryDB (iStiff,ierr)
     iStiff = -1
  end if

  if (nrhs > 0) then
     if (doLogMem) then
        call logAllocMem ('reducer',size(csRhs),0,nbd_p)
     end if
     deallocate(csRhs)
  end if
  if (doLogMem) then
     call logAllocMem ('reducer',size(skee)+size(smee),0,nbd_p)
  end if
  deallocate(skee,smee)
  if (sam%ndof1 > 0) then
     call smDeallocate (skieSparse)
     call smDeallocate (smieSparse)
  end if

  if (noRedMass) then
     sMass = 0.0_dp
     rMass = 0.0_dp
  else

     ! --- Calculate mass properties about link origin and mass center

     call JCMS (sm,tenc,medof,ndim,nenod,rMass,lpu,ierr)
     if (ierr /= 0) then
        call reportError (error_p,'Can not compute inertia properties')
        goto 900
     end if

     call reAllocate ('reducer',tenc)
     call reAllocate ('reducer',medof)

  end if

  if (abs(sMass) <= epsDiv0_p) then
     call reportError (note_p,'This FE part has no mass')
  else if (abs(rMass-sMass) > tolMass_p*sMass) then
     write(errMsg,"('Mass deviation = ',1PE12.5,'  (',2PF6.2,'%)')") &
          &        rMass-sMass, abs(rMass-sMass)/sMass
     call reportError (warning_p,'The total mass computed from the '// &
          'reduced substructure matrices','deviates from the mass '// &
          'computed during finite element assembly.',errMsg)
     call reportError (empty_p,'This could be caused by poor constraint '// &
          'equations from rigid and/or interpolation constraint elements.', &
          'Please check the FE model or contact Fedem support with details.')
  end if

  ! --- Store the reduced substructure matrices for fedem_solver

  call getFileName ('stiffile',chname,'_S.fmx')
  call writeDoubleDB (chname,'stiffness matrix',cs,sk,ierr)
  if (ierr /= 0) then
     call reportError (error_p,'Can not save reduced stiffness matrix to file')
     goto 900
  end if

  if (.not. noRedMass) then

     call getFileName ('massfile',chname,'_M.fmx')
     call writeDoubleDB (chname,'mass matrix',cs,sm,ierr)
     if (ierr /= 0) then
        call reportError (error_p,'Can not save reduced mass matrix to file')
        goto 900
     end if

  end if
  if (calcGdisp) then

     call getFileName ('dispfile',chname,'_V.fmx')
     call writeDoubleDB (chname,'displacement matrix',cs,vgi,ierr)
     if (ierr /= 0) then
        call reportError (error_p,'Can not save displacement matrix to file')
        goto 900
     end if

  end if
  if (nrhs > 0) then

     call getFileName ('loadfile',chname,'_L.fmx')
     call writeDoubleDB (chname,'load vectors',cs,rhs,ierr)
     if (ierr /= 0) then
        call reportError (error_p,'Can not save reduced load vectors to file')
        goto 900
     end if

  end if

  call ffa_cmdlinearg_getint ('nevred',nevred)
  if (nevred > 0 .and. sMass > epsDiv0_p .and. rMass > epsDiv0_p) then

     ! --- Compute eigenvalues of the reduced substructure

     call getFileName ('frsfile',chname)
     if (nevred > ndim) nevred = ndim
     if (chname /= '') then
        allocate(rVal(ndim),rVec(sam%neq,nevred),STAT=ierr)
     else
        allocate(rVal(ndim),rVec(0,0),STAT=ierr)
     end if
     if (ierr /= 0) then
        ierr = allocationError('reducer: Eigenvalues (and expanded '// &
             &                 'eigenvectors) of reduced substructure')
        goto 900
     else if (doLogMem) then
        call logAllocMem ('reducer',0,size(rVal)+size(rVec),nbd_p)
     end if

     call EIGCMS (sk,sm,BmatDisk,tmpVec,rVal,rVec,sam%meqn1,sam%meqn2, &
          &       ndim,sam%ndof1,sam%ndof2,nevred,iprint,lpu,ierr)
     if (ierr < 0) goto 900

     if (chname /= '' .and. ierr == 0) then
        ! Save the expanded mode shapes to frs-file for plotting
        call writeModes (chname,sam,sam%mpar(18),'Free-free reduced modes', &
             &           nevred,rVal,rVec,ierr)
        if (ierr /= 0) then
           call reportError (warning_p,'Failed to write mode shapes file', &
                &                      'Execution continues...')
        end if
     end if

     if (doLogMem) then
        call logAllocMem ('reducer',size(rVal)+size(rVec),0,nbd_p)
     end if
     deallocate(rVal,rVec)

  end if

  if (doLogMem) then
     call logAllocMem ('reducer',size(sk)+size(sm)+size(rhs)+size(vgi),0,nbd_p)
  end if
  deallocate(sk,sm,rhs,vgi)

  if (noGrav) goto 900


  ! --- Calculate and store gravitational forces

  allocate(gravec(ndim,3),STAT=ierr)
  if (ierr /= 0) then
     ierr = allocationError('reducer: Gravity force vectors')
     goto 900
  else if (doLogMem) then
     call logAllocMem ('reducer',0,size(gravec),nbd_p)
  end if

  call GRAV (sam,csMassMat,BmatDisk,tmpVec,gravec,ndim,ngen,lpu,ierr)
  if (ierr /= 0) then
     call reportError (error_p,'Can not compute gravity forces')
     goto 900
  end if

  call getFileName ('gravfile',chname,'_G.fmx')
  call writeDoubleDB (chname,'gravity force vectors',cs,gravec,ierr)
  if (ierr /= 0) then
     call reportError (error_p,'Can not save gravity force vectors to file')
     goto 900
  end if

  if (doLogMem) then
     call logAllocMem ('reducer',size(gravec),0,nbd_p)
  end if
  deallocate(gravec)

  if (neval > 0) then
     call reAllocate ('reducer',eval)
     call reAllocate ('reducer',evec)
  end if
  goto 900


  ! --- Finito

890 continue
  call ffl_done (.true.)
900 continue
  call writeProgress (' --> Closing files and cleaning up')
  if (csStiffMat%storageType == outOfCore_p) then
     nullify(csStiffMat%meqn) ! Points to sam%meqn, avoid deallocation twice
     if (csMassMat%storageType == outOfCore_p) then
        nullify(csMassMat%meqn) ! Points to sam%meqn, avoid deallocation twice
     end if
  else
     nullify(sam%meqn) ! Points to csStiffMat%meqn, avoid deallocation twice
  end if
  i = 0
  if (iStiff >= 0) call closeBinaryDB (iStiff,i)
  if (iMass  >= 0) call closeBinaryDB (iMass,i)
  call dmClose (BmatDisk,i)
  call deAllocateSysMatrix (csStiffMat,i)
  call deAllocateSysMatrix (csMassMat)
  call deAllocateSAM (sam)
  call castToInt8 (sam)
  call releaseScratchArrays
  call writeProgress ('     Finished')

  write(lterm,6000) ' END OF PROGRAM REDUCER '
  call ffa_printMemStatus (lterm)

  if (ierr == 0 .and. i /= 0) ierr = i
  if (ierr /= 0) then
     call reportError (debugFileOnly_p,'reducer')
     call getFileName ('resfile',chname,'.res')
     write(lterm,6100) trim(chname)
     call showTime (lpu)
     write(lpu,6200) 'failed :-('
  else
     call showTime (lpu)
     write(lpu,6200) 'successfully completed :-)'
  end if

  call closeLogFiles ()
  call ffa_stoptimer ('fedem_reducer')
  call ffa_reporttimer

6000 format(/11X,16('='),'> ',A,' <',16('='))
6010 format(//4X,'COMPONENT MODE FREQUENCIES' &
          &  /4X,'--------------------------' &
          & //4X,'MODE  FREQUENCY [Hz]'/ 1P,(I7,E17.8))
6100 format(/11X,'An error is detected. Error messages are written to ',A)
6200 format(/ 4X,'Model reduction ',A)
6300 format(//1X,62('-') &
          & /' --- Internal node DOFs with lowest connection to external DOFs' &
          & /1X,62('-'))
6310 format(' Node',I7,' DOF',I3,' Max connection value ',1PE12.5)

end subroutine reducer
