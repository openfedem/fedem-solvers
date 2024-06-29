!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

subroutine cfemReducer (numStiff,ierr)
#ifdef FT_HAS_CFEM

  use KindModule                , only : dp, i8, nbd_p, lfnam_p
  use KindModule                , only : epsDiv0_p, pi_p
  use SamModule                 , only : SamType, deallocateSAM
  use SamReducerModule          , only : saveSAM
  use SysMatrixTypeModule       , only : SysMatrixType, diagonalMatrix_p
  use SysMatrixTypeModule       , only : saveSysMat, deallocateSysMatrix
  use SparseMatrixModule        , only : SparseMatrixType, smDeallocate
  use DiskMatrixModule          , only : DiskMatrixType, dmSize, dmGetSwapSize
  use DiskMatrixModule          , only : dmNullify, dmOpen, dmClose, dmNdp_p
  use AsmExtensionModule        , only : csBeginAssembly, csEndAssembly
  use InaddModule               , only : INADD, extractSubMat
  use CfemFedemModule           , only : readReducerDataAndSolveCfem
  use CmstrsModule              , only : CMSTRS, EIGVAL, EIGCMS, JCMS, GRAV
  use FileUtilitiesModule       , only : getFileName
  use TimerModule               , only : initTime, showTime
  use VersionModule             , only : openResFile
  use ProgressModule            , only : lterm, writeProgress
  use FileUtilitiesModule       , only : findUnitNumber
  use AllocationModule          , only : doLogMem, logAllocMem, reAllocate
  use ScratchArrayModule        , only : releaseScratchArrays
  use ReportErrorModule         , only : empty_p, note_p, warning_p, error_p
  use ReportErrorModule         , only : debugFileOnly_p, openTerminalOutputFile
  use ReportErrorModule         , only : allocationError, reportError
  use BinaryDBInterface         , only : getPositionDB, setPositionDB
  use BinaryDBInterface         , only : readDoubleDB, writeDoubleDB
  use BinaryDBInterface         , only : closeBinaryDB
  use FFaCmdLineArgInterface    , only : ffa_cmdlinearg_getint
  use FFaCmdLineArgInterface    , only : ffa_cmdlinearg_getbool
  use FFaCmdLineArgInterface    , only : ffa_cmdlinearg_getdouble
  use FFlLinkHandlerInterface   , only : ffl_done, ffl_calcs, ffl_getcs
  use FFlLinkHandlerInterface   , only : ffl_addcs_int, ffl_addcs_double

  implicit none

  integer , intent(in)   :: numStiff
  integer , intent(out)  :: ierr
  logical                :: diagMass, lumpedMass, factorMass
  integer                :: i, imass, istiff, iopSing, lpu, iprint
  integer                :: cs, ndim, nenod, neval, ngen, nevred
  real(dp), parameter    :: tolMass_p = 1.0e-6_dp
  real(dp)               :: sMass, rMass, tolEigval, tolFactorize, eigenShift
  real(dp), allocatable  :: smee(:,:), skee(:,:), vgi(:,:)
  real(dp), allocatable  :: rval(:), rvec(:,:), gravec(:,:), sm(:,:), sk(:,:)
  real(dp), pointer      :: tenc(:,:), eval(:), evec(:,:)
  integer , pointer      :: medof(:)
  type(SamType)          :: sam
  type(SysMatrixType)    :: csStiffMat, csMassMat
  type(SparseMatrixType) :: smieSparse, skieSparse
  type(DiskMatrixType)   :: BmatDisk
  character(len=lfnam_p) :: chName
  character(len=64)      :: errMsg

  ! for nonlinear reduction, by LIM
  integer(i8)            :: currPos
  integer                :: nstif, active_stiff, iUnit
  real(dp), allocatable  :: force_red(:,:), displ_red(:,:), K_red(:,:), M_red(:,:)
  real(dp), pointer      :: stiff_f2c(:), force_disp_f2c(:,:)


  !! --- Logic section ---

  call initTime ()
  call openTerminalOutputFile (lterm)
  write(lterm,6000) 'START OF PROGRAM REDUCER'

  ierr = 0
  iMass = -1
  iStiff = -1
  nullify(tenc)
  nullify(medof)
  call dmNullify (BmatDisk)

  ! Open the reducer result file and write heading
  call getFileName ('resfile',chname,'.res')
  call openResFile (chname,'Non-linear FE model Reducer',lpu,ierr)
  if (ierr /= 0) return

  call ffa_cmdlinearg_getint ('debug',iprint)
  call ffa_cmdlinearg_getint ('singularityHandler',iopSing)
  call ffa_cmdlinearg_getbool ('diagmass',diagMass)
  call ffa_cmdlinearg_getbool ('lumpedmass',lumpedMass)
  call ffa_cmdlinearg_getbool ('factorMass',factorMass)
  call ffa_cmdlinearg_getdouble ('tolEigval',tolEigval)
  call ffa_cmdlinearg_getdouble ('tolFactorize',tolFactorize)
  call ffa_cmdlinearg_getdouble ('eigenshift',eigenShift)


  ! --- Read the link file and establish the SAM datastructure

  call writeProgress (' --> Reading FE data file,'// &
       &              ' solving for nonlinear solutions sets')
  call getFileName ('CfemFile',chName,'.dat')
  call readReducerDataAndSolveCfem (sam, csStiffMat, csMassMat, tenc, medof, &
       &                            diagMass, lumpedMass, iprint, lpu, ierr, &
       &                            stiff_f2c, force_disp_f2c, &
       &                            numStiff, sam%mpar(40), chName)
  if (ierr /= 0) then
     call reAllocate ('reducer',tenc)
     call reAllocate ('reducer',medof)
     goto 900
  else if (diagMass .and. sam%nmmceq > 0) then
     call reportError (warning_p,'This FE part has constraint equations.', &
          'Forcing a diagonal mass matrix may yield poor convergence', &
          'and/or inaccurate results in the dynamics simulation.', &
          'Switching to lumped- or consistent element mass is recommended.')
  end if

  call ffa_cmdlinearg_getint ('printArray',sam%mpar(25))
  if (sam%mpar(25) < 100) then
     sam%mpar(26) = sam%mpar(25)
  else
     sam%mpar(26) = mod(sam%mpar(25),100)
     sam%mpar(25) =     sam%mpar(25)/100
  end if

  ! --- Calculate the link checksum

  call ffl_calcs (ierr)
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
  if (ierr /= 0) goto 900

  if (tolFactorize < 1.0e-20_dp) then
     tolFactorize = 1.0e-20_dp
     call reportError (note_p,'A factorization tolerance lower than 1e-20'// &
          ' is not allowed.','The factorization tolerance is reset to 1e-20')
  end if


  ! --- Allocate substructure stiffness- and mass matrices

  call writeProgress (' --> Allocating substructure matrices')
  call csBeginAssembly (csStiffMat,'system stiffness matrix',lpu,ierr)
  if (ierr /= 0) goto 900
  call csBeginAssembly (csMassMat,'system mass matrix',lpu,ierr)
  if (ierr /= 0) goto 900

  ! --- Now assemble the substructure mass matrix

  call INADD (sam, csStiffMat, csMassMat, sMass, rMass, &
       &      diagMass.or.lumpedMass, .false., 2, iprint, lpu, ierr)
  if (ierr /= 0) then
     call reportError (error_p,'Can not build substructure matrices')
     goto 900
  end if

  call ffl_done ! Assembly finished, release the link object

  call csEndAssembly (csStiffMat,ierr)
  if (ierr /= 0) goto 900
  call csEndAssembly (csMassMat,ierr)
  if (ierr /= 0) goto 900

  ! --- Store the substructure mass- and stiffness matrices on temporary files

  call writeProgress ('     Saving substructure matrices to temporary files')

  currPos = size(csStiffMat%value,kind=i8)*int(sam%mpar(40),i8)
  eval => csStiffMat%value
  csStiffMat%value => stiff_f2c
  call saveSysMat (csStiffMat,'stiffness matrices',currPos,iStiff,ierr)
  if (ierr < 0) goto 900
  csStiffMat%value => eval
  deallocate(stiff_f2c)
  call writeProgress (1,2)

  currPos = size(csMassMat%value,kind=i8)
  call saveSysMat (csMassMat,'mass matrix',currPos,iMass,ierr)
  if (ierr < 0) goto 900
  call writeProgress (2,2)


  neval = sam%mpar(19)
  ngen  = sam%mpar(22)
  ndim  = sam%mpar(24)
  nenod = sam%mpar(29)
  nstif = sam%mpar(40)

  nullify(eval)
  nullify(evec)

  ! --- Allocate reduced mass/stiffness matrices

  allocate(sm(ndim,ndim),sk(ndim,ndim),vgi(0,0),STAT=ierr)
  if (ierr /= 0) then
     ierr = allocationError('reducer: Reduced substructure matrices')
     goto 900
  else if (doLogMem) then
     call logAllocMem ('reducer',0,size(sm)+size(sk)+size(vgi),nbd_p)
  end if

  ! --- Allocate submatrices ie and ee

  if (csMassMat%storageType == diagonalMatrix_p) then
     i = 1
  else
     i = sam%ndof2
  end if
  allocate(smee(sam%ndof2,i),skee(sam%ndof2,sam%ndof2),STAT=ierr)
  if (ierr /= 0) then
     ierr = allocationError('reducer: Submatrices smee and skee')
     goto 900
  else if (doLogMem) then
     call logAllocMem ('reducer',0,size(smee)+size(skee),nbd_p)
  end if

  ! Only keep force-displacement vectors for external DOFs for nonlinear results
  allocate(force_red(ndim,nstif), &
       &   displ_red(ndim,nstif), &
       &   K_red(ndim*ndim,nstif), &
       &   M_red(ndim*ndim,nstif), STAT=ierr)
  if (ierr /= 0) then
     ierr = allocationError('reducer: force_red, displ_red, K_red and M_red')
     goto 900
  else if (doLogMem) then
     call logAllocMem ('reducer',0,2*size(force_red)+2*size(K_red),nbd_p)
  end if
  K_red = 0.0_dp
  M_red = 0.0_dp


  ! --- Loop through nonlinear stiffness system

  do active_stiff = 1, nstif

     write(errMsg,"(' ==> CMSTRS for stiffness number',i3,' ===')") active_stiff
     call writeProgress (errMsg)

     if (active_stiff == 1) then
        ! Rewind, and read the first stiffness matrix from file
        call readDoubleDB (-iStiff,csStiffMat%value(1),size(csStiffMat%value),ierr)
        if (ierr < 0) then
           call reportError (error_p,'Can not restore system stiffness matrix')
           goto 900
        end if
     else
        ! Read stiffness matrix number "active_stiff" from file
        call getPositionDB (iStiff,currPos)
        call readDoubleDB (iStiff,csStiffMat%value(1),size(csStiffMat%value),ierr)
        if (currPos <= 0_i8 .or. ierr < 0) then
           call reportError (error_p,'Can not restore system stiffness matrix')
           goto 900
        end if
     end if

     ! --- Extract submatrices ie and ee from the substructure matrices

     call extractSubMat (sam, csStiffMat, csMassMat, &
          &              smieSparse, smee, skieSparse, skee, iprint, lpu, ierr)
     if (ierr /= 0) goto 900

     ! --- Allocate disk swap matrix

     call getFileName ('Bmatfile',chname,'_B.fmx')
     call dmOpen (BmatDisk,chname,cs,sam%ndof1,sam%ndof2,dmNdp_p,ierr, &
          &       swapSize=dmGetSwapSize())
     if (ierr < 0) goto 900

     sam%mpar(31) = int(dmSize(BmatDisk,3))

     if (active_stiff == 1) then

        ! --- Store all SAM-arrays on file for fedem_stress

        call getFileName ('samfile',chname,'_SAM.fsm')
        call saveSAM (chname,cs,sam,ierr)
        if (ierr < 0) then
           call reportError (error_p,'Can not save SAM-data to file '//chname)
           goto 900
        end if

     end if

     if (neval > 0) then

        ! --- Compute eigenvalues and eigenvectors of the substructure

        call EIGVAL (sam, csStiffMat, csMassMat, neval, factorMass, &
             &       iopSing, tolFactorize, tolEigval, eigenShift, &
             &       eval, evec, iprint, lpu, ierr)
        if (ierr < 0) then
           call reportError (error_p,'Can not perform eigenvalue calculation')
           goto 900
        end if

        write(lpu,6010) (i,eval(i)*0.5_dp/pi_p,i=1,neval)

        call getFileName ('eigfile',chname,'_E.fmx')
        call writeDoubleDB (chname,'generalized modes',cs,evec,ierr)
        if (ierr /= 0) then
           call reportError (error_p,'Can not write eigenvectors to file')
           goto 900
        end if

        if (active_stiff > 1) then
           ! Read matrix number "active_stiff" from file
           call setPositionDB (iStiff,currPos,ierr)
           call readDoubleDB (iStiff,csStiffMat%value(1),size(csStiffMat%value),ierr)
        else
           ! Rewind, and read the first matrix from file
           call readDoubleDB (-iStiff,csStiffMat%value(1),size(csStiffMat%value),ierr)
        end if
        if (ierr < 0) then
           call reportError (error_p,'Can not restore system stiffness matrix')
           goto 900
        end if

        call readDoubleDB (-iMass,csMassMat%value(1),size(csMassMat%value),ierr)
        if (ierr < 0) then
           call reportError (error_p,'Can not restore system mass matrix')
           goto 900
        end if

     end if

     ! --- Perform the CMS-transformation

     call CMSTRS (iopSing, sam, BmatDisk, &
          &       csMassMat, smieSparse, smee, sm, csStiffMat, skee, sk, &
          &       evec, evec, evec, vgi, tolFactorize, &
          &       iStiff, ngen, 0, lpu, ierr)
     if (ierr /= 0) then
        call reportError (error_p,'Can not perform CMS-transformation')
        goto 900
     end if

     ! --- Calculate mass properties about link origin and mass center

     call JCMS (sm,tenc,medof,ndim,nenod,rMass,lpu,ierr)
     if (ierr /= 0) then
        call reportError (error_p,'Can not compute inertia properties')
        goto 900
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

     K_red(:,active_stiff) = reshape(sk,(/ndim*ndim/))
     M_red(:,active_stiff) = reshape(sm,(/ndim*ndim/))
     ! Reduced force- and displacement vectors
     i = (active_stiff-1)*sam%ndof + sam%ndof1
     force_red(1:sam%ndof2,active_stiff) = force_disp_f2c(i+1:i+sam%ndof2,1)
     displ_red(1:sam%ndof2,active_stiff) = force_disp_f2c(i+1:i+sam%ndof2,2)

     call ffa_cmdlinearg_getint ('nevred',nevred)
     if (nevred > 0 .and. smass > epsDiv0_p .and. rmass > epsDiv0_p) then

        ! --- Compute eigenvalues of the reduced substructure

        if (nevred > ndim) nevred = ndim
        allocate(rVal(ndim),rVec(0,0),STAT=ierr)
        if (ierr /= 0) then
           ierr = allocationError('reducer: Eigenvalues of reduced system')
           goto 900
        else if (doLogMem) then
           call logAllocMem ('reducer',0,size(rVal)+size(rVec),nbd_p)
        end if

        call EIGCMS (sk,sm,BmatDisk,evec,rVal,rVec,sam%meqn1,sam%meqn2, &
             &       ndim,sam%ndof1,sam%ndof2,nevred,iprint,lpu,ierr)
        if (ierr < 0) goto 900

        if (doLogMem) then
           call logAllocMem ('reducer',size(rVal)+size(rVec),0,nbd_p)
        end if
        deallocate(rVal,rVec)

     end if


     ! --- Calculate and store gravitational forces

     allocate(gravec(ndim,3),STAT=ierr)
     if (ierr /= 0) then
        ierr = allocationError('reducer: Gravity force vectors')
        goto 900
     else if (doLogMem) then
        call logAllocMem ('reducer',0,size(gravec),nbd_p)
     end if

     call GRAV (sam,csMassMat,BmatDisk,evec,gravec,ndim,ngen,lpu,ierr)
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

     call dmClose (BmatDisk,ierr)
     if (ierr < 0) goto 900

  end do

  call reAllocate ('reducer',tenc)
  call reAllocate ('reducer',medof)

  deallocate(force_disp_f2c)

  if (doLogMem) then
     call logAllocMem ('reducer',size(skee)+size(smee),0,nbd_p)
     call logAllocMem ('reducer',size(sk)+size(sm)+size(vgi),0,nbd_p)
  end if
  deallocate(skee,smee)
  deallocate(sk,sm,vgi)

  if (sam%ndof1 > 0) then
     call smDeallocate (skieSparse)
     call smDeallocate (smieSparse)
  end if


  ! --- Store the reduced substructure matrices for fedem_solver

  call getFileName ('stiffile',chname,'_S.fmx')
  call writeDoubleDB (chname,'stiffness matrix',cs,K_red,ierr)
  if (ierr /= 0) then
     call reportError (error_p,'Can not save reduced stiffness matrix to file')
     goto 900
  end if

  call getFileName ('massfile',chname,'_M.fmx')
  call writeDoubleDB (chname,'mass matrix',cs,M_red,ierr)
  if (ierr /= 0) then
     call reportError (error_p,'Can not save reduced mass matrix to file')
     goto 900
  end if

  ! --- Store nonlinear force and displacement vectors

  call getFileName ('forcefile',chname,'_F.fmx')
  call writeDoubleDB (chname,'force vectors',cs,force_red,ierr)
  if (ierr /= 0) then
     call reportError (error_p,'Can not save NL force vectors to file')
     goto 900
  end if

  call getFileName ('dispfile',chname,'_D.fmx')
  call writeDoubleDB (chname,'displacement vectors',cs,displ_red,ierr)
  if (ierr /= 0) then
     call reportError (error_p,'Can not save NL disp vectors to file')
     goto 900
  end if

  if (doLogMem) then
     call logAllocMem ('reducer',size(K_red)+size(M_red)+ &
          &                      size(force_red)+size(displ_red),0,nbd_p)
  end if
  deallocate(K_red,M_red,force_red,displ_red)

  call getFileName ('numStatesFile',chname,'_numStates.txt')
  iUnit = findUnitNumber(30)
  open(iUnit,FILE=chname,STATUS='UNKNOWN')
  write(iUnit,*) nstif
  close(iUnit)


  ! --- Finito

900 continue
  call writeProgress (' --> Closing files and cleaning up')
  call ffl_done
  nullify(sam%meqn) ! Points to csStiffMat%meqn, avoid deallocation twice
  i = 0
  if (iStiff >= 0) call closeBinaryDB (iStiff,i)
  if (iMass  >= 0) call closeBinaryDB (iMass,i)
  call deAllocateSysMatrix (csStiffMat,i)
  call deAllocateSysMatrix (csMassMat)
  call deAllocateSAM (sam)
  call releaseScratchArrays
  call writeProgress ('     Finished')

  write(lterm,6000) ' END OF PROGRAM REDUCER '

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

6000 format(/11X,16('='),'> ',A,' <',16('='))
6010 format(//4X,'COMPONENT MODE FREQUENCIES' &
          &  /4X,'--------------------------' &
          & //4X,'MODE  FREQUENCY [Hz]'/ 1P,(I7,E17.8))
6100 format(/11X,'An error is detected. Error messages are written to ',A)
6200 format(/ 4X,'Model reduction ',A)

#else
  !! Dummy version
  integer, intent(in)  :: numStiff
  integer, intent(out) :: ierr
  print *,'*** cfemReducer dummy:',numStiff
  ierr = -1
#endif

end subroutine cfemReducer
