!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file smallReducer.f90
!> @brief Simplified driver for the FE part reducer.
!>
!> @details This driver does essentially the same as the main driver reducer()
!> except for no component modes option. It operates on full matrices only.
!> This is used to verify the FE part reduction procedure on smaller models.
!>
!> @author Knut Morten Okstad
!>
!> @date 26 Feb 2018

!!==============================================================================
!> @brief Test driver for the FE part reducer operating on full matrices only.

subroutine smallReducer (ierr)

  use KindModule             , only : dp, lfnam_p
  use SamModule              , only : SamType, deallocateSAM
  use SysMatrixTypeModule    , only : SysMatrixType, deallocateSysMatrix
  use AsmExtensionModule     , only : csBeginAssembly, csEndAssembly, castToInt8
  use MatExtensionModule     , only : csGetSub11, csGetSub12, csGetSub22
  use InputReducerModule     , only : readReducerData
  use SamReducerModule       , only : saveSAM
  use InaddModule            , only : INADD
  use CmstrsModule           , only : GRAVfull
  use FileUtilitiesModule    , only : getFileName
  use TimerModule            , only : initTime, showTime
  use VersionModule          , only : openResFile
  use ProgressModule         , only : lterm, writeProgress
  use ScratchArrayModule     , only : releaseScratchArrays
  use ManipMatrixModule      , only : writeObject
  use ReportErrorModule      , only : error_p, debugFileOnly_p, reportError
  use ReportErrorModule      , only : allocationError, openTerminalOutputFile
  use BinaryDBInterface      , only : writeDoubleDB
  use FFaCmdLineArgInterface , only : ffa_cmdlinearg_getint
  use FFlLinkHandlerInterface, only : ffl_done

  implicit none

  integer , intent(out)  :: ierr
  integer                :: i, lpu, iprint, ndim, neqi, mlc(100)
  real(dp)               :: sMass, rMass
  integer , allocatable  :: ipiv(:)
  real(dp), allocatable  :: sm(:,:), smii(:,:), smee(:,:), smie(:,:), BMtmp(:,:)
  real(dp), allocatable  :: sk(:,:), skii(:,:), skee(:,:), skie(:,:), Bmat(:,:)
  real(dp), allocatable  :: gFull(:,:), gravec(:,:), vgi(:,:)
  real(dp), pointer      :: tenc(:,:)
  integer , pointer      :: medof(:)
  type(SamType)          :: sam
  type(SysMatrixType)    :: csStiffMat, csMassMat
  character(len=lfnam_p) :: chName

  !! --- Logic section ---

  call initTime ()
  call openTerminalOutputFile (lterm)
  write(lterm,6000) 'START OF PROGRAM small REDUCER'

  nullify(tenc)
  nullify(medof)

  ! Open the reducer result file and write heading
  call getFileName ('resfile',chname,'.res')
  call openResFile (chname,'FE model Reducer',lpu,ierr)
  if (ierr /= 0) return

  call ffa_cmdlinearg_getint ('debug',iprint)

  ! --- Read the link file and establish the SAM datastructure

  call writeProgress (' --> Reading FE data file')
  call readReducerData (sam, csStiffMat, csMassMat, sMass, rMass, &
       &                tenc, medof, mlc, .false., .false., .false., .false., &
       &                iprint, lpu, ierr)
  if (ierr /= 0) goto 900

  ! --- Allocate substructure stiffness- and mass matrices

  deallocate(tenc,medof)
  call writeProgress (' --> Allocating substructure matrices')
  call csBeginAssembly (csStiffMat,'system stiffness matrix',lpu,ierr)
  if (ierr /= 0) goto 890
  call csBeginAssembly (csMassMat,'system mass matrix',lpu,ierr)
  if (ierr /= 0) goto 890

  ! --- Now assemble the substructure mass- and stiffness matrices

  call INADD (sam, csStiffMat, csMassMat, sMass, rMass, &
       &      .false., .false., 3, iprint, lpu, ierr)
  if (ierr /= 0) then
     call reportError (error_p,'Can not build substructure matrices')
     goto 890
  end if

  call ffl_done (.true.) ! Assembly finished, release the link object

  call csEndAssembly (csStiffMat,ierr)
  if (ierr /= 0) goto 900
  call csEndAssembly (csMassMat,ierr)
  if (ierr /= 0) goto 900

  ! --- Store all SAM-arrays on file

  call getFileName ('samfile',chname,'_SAM.fsm')
  call saveSAM (chname,1,sam,ierr)
  if (ierr < 0) then
     call reportError (error_p,'Can not save SAM-data to file '//chname)
     goto 900
  end if

  ! --- Allocate submatrices

  allocate(smii(sam%ndof1,sam%ndof1),skii(sam%ndof1,sam%ndof1), &
       &   smee(sam%ndof2,sam%ndof2),skee(sam%ndof2,sam%ndof2), &
       &   smie(sam%ndof1,sam%ndof2),skie(sam%ndof1,sam%ndof2), &
       &   Bmat(sam%ndof1,sam%ndof2),ipiv(sam%ndof1), STAT=ierr)
  if (ierr /= 0) then
     ierr = allocationError('reducer: Submatrices ii, ee and ie')
     goto 900
  end if

  ! --- Extract submatrices ii, ie and ee from the substructure matrices

  call csGetSub11 (csStiffMat,sam%meqn1,skii,ierr)
  if (ierr /= 0) goto 900
  call csGetSub11 (csMassMat,sam%meqn1,smii,ierr)
  if (ierr /= 0) goto 900

  call csGetSub12 (csStiffMat,sam%meqn1,sam%meqn2,skie,ierr)
  if (ierr /= 0) goto 900
  call csGetSub12 (csMassMat,sam%meqn1,sam%meqn2,smie,ierr)
  if (ierr /= 0) goto 900

  call csGetSub22 (csStiffMat,sam%meqn2,skee,ierr)
  if (ierr /= 0) goto 900
  call csGetSub22 (csMassMat,sam%meqn2,smee,ierr)
  if (ierr /= 0) goto 900

  if (iprint > 2) then
     call writeObject (SKii,lpu,'SKii')
     call writeObject (SMii,lpu,'SMii')
     call writeObject (SMie,lpu,'SMie')
     call writeObject (SKie,lpu,'SKie')
     call writeObject (SMee,lpu,'SMee')
     call writeObject (SKee,lpu,'SKee')
  end if

  ! --- Calculate the B-matrix, B = Kii^-1 * Kie

  neqi = sam%ndof1
  ndim = sam%ndof2
  call DCOPY (neqi*ndim,skie(1,1),1,Bmat(1,1),1)
  call DGESV (neqi,ndim,skii(1,1),neqi,ipiv(1),Bmat(1,1),neqi,ierr)
  if (ierr /= 0) then
     call reportError (error_p,'DGESV returned IERR = ',ierr=ierr)
     goto 900
  end if

  ! --- Allocate reduced mass/stiffness matrices

  allocate(sk(ndim,ndim),sm(ndim,ndim),BMtmp(neqi,ndim),STAT=ierr)
  if (ierr /= 0) then
     ierr = allocationError('reducer: Reduced substructure matrices')
     goto 900
  end if

  ! --- Perform the CMS-transformation, sK = Kee + B^t*Kie

  call DSCAL (ndim*neqi,-1.0_dp,Bmat(1,1),1)
  call DCOPY (ndim*ndim,skee(1,1),1,sk(1,1),1)       ! SK = SKee
  call DGEMM ('T','N',ndim,ndim,neqi,1.0_dp,Bmat(1,1),neqi, &
       &      skie(1,1),neqi,1.0_dp,sk(1,1),ndim)    ! SK += B^t*SKie

  ! --- Perform the CMS-transformation, sM = Mee + Mie^t*B + B^t*Mie + B^t*Mii*B

  call DCOPY (ndim*ndim,smee(1,1),1,sm(1,1),1)       ! SM = SMee
  call DGEMM ('T','N',ndim,ndim,neqi,1.0_dp,smie(1,1),neqi, &
       &      Bmat(1,1),neqi,1.0_dp,sm(1,1),ndim)    ! SM += SMie^t*B
  call DCOPY (neqi*ndim,SMie(1,1),1,BMtmp(1,1),1)    ! BMtmp = SMie
  call DGEMM ('N','N',neqi,ndim,neqi,1.0_dp,smii(1,1),neqi, &
       &      Bmat(1,1),neqi,1.0_dp,BMtmp(1,1),neqi) ! BMtmp += SMii*B
  call DGEMM ('T','N',ndim,ndim,neqi,1.0_dp,Bmat(1,1),neqi, &
       &      BMtmp(1,1),neqi,1.0_dp,sm(1,1),ndim)   ! SM += B^t*BMtmp

  ! --- Calculate gravitational forces

  allocate(gFull(sam%neq,3),gravec(ndim,3),vgi(neqi,3),STAT=ierr)
  if (ierr /= 0) then
     ierr = allocationError('reducer: Gravity force vectors')
     goto 900
  end if

  call GRAVfull (sam,csMassMat,gFull,ierr)
  if (ierr /= 0) then
     call reportError (error_p,'Can not compute gravity forces')
     goto 900
  end if

  if (iprint > 2) then
     call writeObject (gFull,lpu,'Qg')
  end if

  do i = 1, neqi
     vgi(i,:)    = gFull(sam%meqn1(i),:)
  end do
  do i = 1, ndim
     gravec(i,:) = gFull(sam%meqn2(i),:)
  end do
  call DGEMM ('T','N',ndim,3,neqi,1.0_dp,Bmat(1,1),neqi, &
       &      vgi(1,1),neqi,1.0_dp,gravec(1,1),ndim) ! g += B^t*vgi

  ! --- Solve for the internal displacement due to gravity, vgi = Kii^-1 * Qgi

  call DGETRS ('N',neqi,3,skii(1,1),neqi,ipiv(1),vgi(1,1),neqi,ierr)
  if (ierr /= 0) then
     call reportError (error_p,'DGETRS returned IERR = ',ierr=ierr)
     goto 900
  end if

  if (iprint > 1) then
     call writeObject (SM,lpu,'SMMAT')
     call writeObject (SK,lpu,'SKMAT')
     call writeObject (VGI,lpu,'VGI')
     call writeObject (Bmat,lpu,'BMAT',10)
     call writeObject (gravec,lpu,'GRAV')
  end if

  ! --- Store the reduced substructure matrices for fedem_solver

  call writeDoubleDB ('fedem_S.fmx','stiffness matrix',i,sk,ierr)
  if (ierr /= 0) then
     call reportError (error_p,'Can not save reduced stiffness matrix to file')
     goto 900
  end if

  call writeDoubleDB ('fedem_M.fmx','mass matrix',i,sm,ierr)
  if (ierr /= 0) then
     call reportError (error_p,'Can not save reduced mass matrix to file')
     goto 900
  end if

  call writeDoubleDB ('fedem_V.fmx','displacement matrix',i,vgi,ierr)
  if (ierr /= 0) then
     call reportError (error_p,'Can not save displacement matrix to file')
     goto 900
  end if

  call writeDoubleDB ('fedem_B.fmx','disk matrix',i,Bmat,ierr)
  if (ierr /= 0) then
     call reportError (error_p,'Can not save disk matrix to file')
     goto 900
  end if

  call writeDoubleDB ('fedem_G.fmx','gravity force vectors',i,gravec,ierr)
  if (ierr /= 0) then
     call reportError (error_p,'Can not save gravity force vectors to file')
     goto 900
  end if

  deallocate(skii,skee,skie,smii,smee,smie)
  deallocate(Bmat,BMtmp,ipiv,sk,sm,gFull,gravec,vgi)
  goto 900

  ! --- Finito

890 continue
  call ffl_done (.true.)
900 continue
  call writeProgress (' --> Closing files and cleaning up')
  nullify(sam%meqn) ! Points to csStiffMat%meqn, avoid deallocation twice
  call deAllocateSysMatrix (csStiffMat)
  call deAllocateSysMatrix (csMassMat)
  call deAllocateSAM (sam)
  call castToInt8 (sam)
  call releaseScratchArrays
  call writeProgress ('     Finished')

  write(lterm,6000) ' END OF PROGRAM small REDUCER '

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
6100 format(/11X,'An error is detected. Error messages are written to ',A)
6200 format(/ 4X,'Model reduction ',A)

end subroutine smallReducer
