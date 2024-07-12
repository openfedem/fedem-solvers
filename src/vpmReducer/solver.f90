!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file solver.f90
!> @brief Driver for the FEDEM linear FE part solver.

!!==============================================================================
!> @brief Driver for the FEDEM linear FE part solver.
!>
!> @param[in] n Number of nodes to return displacements for
!> @param[in] nodes List of nodes to return displacements for
!> @param[out] displ Array of nodal displacements
!> @param[out] ierr Error flag
!>
!> @details This program may be used to solve the linear FE equations on a part
!> subjected to some external load. Its main purpose is for development use,
!> to simplify the verification of FE formulations, constraint handling, etc.
!> The program can also solve for the natural frequencies (eigenvalues).
!>
!> @callgraph
!>
!> @author Knut Morten Okstad
!>
!> @date 11 Oct 2002

subroutine solver (n,nodes,displ,ierr)

  use sprKindModule             , only : dp, ik
  use KindModule                , only : lfnam_p, pi_p
  use SamModule                 , only : SamType
  use SamModule                 , only : NodeNumFromEqNum, deAllocateSAM
  use SamReducerModule          , only : saveSAM
  use SaveReducerModule         , only : writeModes, writeDeformation
  use SysMatrixTypeModule       , only : SysMatrixType, outOfCore_p
  use SysMatrixTypeModule       , only : deAllocateSysMatrix
  use AsmExtensionModule        , only : csBeginAssembly, csEndAssembly
  use AsmExtensionModule        , only : castToInt8
  use SolExtensionModule        , only : csSolve
  use InputReducerModule        , only : readReducerData
  use InaddModule               , only : INADD
  use CmstrsModule              , only : EIGVAL
  use FileUtilitiesModule       , only : getFileName
  use TimerModule               , only : initTime, showTime
  use VersionModule             , only : openResFile
  use ProgressModule            , only : lterm, writeProgress
  use ManipMatrixModule         , only : writeObject
  use ScratchArrayModule        , only : releaseScratchArrays
  use IdTypeModule              , only : StrId
  use ReportErrorModule         , only : openTerminalOutputFile
  use ReportErrorModule         , only : allocationError, reportError, note_p
  use ReportErrorModule         , only : error_p, warning_p, debugFileOnly_p
  use BinaryDBInterface         , only : closeBinaryDB
  use FFaCmdLineArgInterface    , only : ffa_cmdlinearg_getint
  use FFaCmdLineArgInterface    , only : ffa_cmdlinearg_getbool
  use FFaCmdLineArgInterface    , only : ffa_cmdlinearg_getdouble
  use FFaCmdLineArgInterface    , only : ffa_cmdlinearg_getdoubles
  use FFaCmdLineArgInterface    , only : ffa_cmdlinearg_isSet
  use FFlLinkHandlerInterface   , only : ffl_done

  implicit none

  integer, intent(in)    :: n
  integer, intent(inout) :: nodes(n)
  real(dp), intent(out)  :: displ(3,n)
  integer , intent(out)  :: ierr
  logical                :: lumpedMass, factorMass
  integer                :: i, iopAdd, iopSing, lpu, iprint
  integer                :: j, k, ieq, neval, nMax, nrhs, mlc(100)
  integer(ik)            :: eqnInErr
  real(dp)               :: dMax, sMass, rMass, grav(3)
  real(dp)               :: tolEigval, tolFactorize, eigenShift
  real(dp), allocatable  :: rhs(:,:), Qg(:,:)
  real(dp), pointer      :: tenc(:,:), eval(:), evec(:,:)
  integer , pointer      :: medof(:)
  type(SamType)          :: sam
  type(SysMatrixType)    :: csStiffMat, csMassMat
  character(len=lfnam_p) :: chName

  !! --- Logic section ---

  call initTime ()
  call openTerminalOutputFile (lterm)
  write(lterm,6000) 'START OF PROGRAM SOLVER'

  grav = 0.0_dp
  nullify(tenc)
  nullify(eval)
  nullify(evec)
  nullify(medof)

  ! Open the solver result file and write heading
  call getFileName ('resfile',chname,'.res')
  call openResFile (chname,'FE model Solver',lpu,ierr)
  if (ierr /= 0) return

  call ffa_cmdlinearg_getint ('debug',iprint)
  call ffa_cmdlinearg_getint ('singularityHandler',iopSing)
  call ffa_cmdlinearg_getbool ('lumpedmass',lumpedMass)
  call ffa_cmdlinearg_getdouble ('tolFactorize',tolFactorize)
  call ffa_cmdlinearg_getdoubles ('gvec',grav,3)
  if (grav(1)*grav(1) + grav(2)*grav(2) + grav(3)*grav(3) > 1.0e-8_dp) then
     nrhs = 1
  else
     nrhs = 0
  end if

  ! --- Read the link file and establish the SAM datastructure

  call writeProgress (' --> Reading FE data file')
  call readReducerData (sam, csStiffMat, csMassMat, sMass, rMass, &
       &                tenc, medof, mlc, .false., .false., lumpedMass, &
       &                .false., iprint, lpu, ierr, nrhs)
  if (ierr /= 0) goto 900

  call ffa_cmdlinearg_getint ('printArray',sam%mpar(25))
  if (sam%mpar(25) < 100) then
     sam%mpar(26) = sam%mpar(25)
  else
     sam%mpar(26) = mod(sam%mpar(25),100)
     sam%mpar(25) =     sam%mpar(25)/100
  end if

  if (tolFactorize < 1.0e-20_dp) then
     tolFactorize = 1.0e-20_dp
     call reportError (note_p,'A factorization tolerance lower than 1e-20'// &
          ' is not allowed.','The factorization tolerance is reset to 1e-20')
  end if

  ! --- Assemble the system load vectors, if any

  neval = sam%mpar(19)
  if (neval > 0) then
     nrhs = 0
  else
     nrhs = max(sam%mpar(21),sam%mpar(30))
  end if
  if (nrhs > 0) then

     allocate(rhs(sam%neq,nrhs),STAT=ierr)
     if (ierr /= 0) then
        ierr = allocationError('Right-hand-side vectors')
        goto 890
     end if

     rhs = 0.0_dp
     do i = 1, sam%mpar(21)
        call writeProgress (' --> INADD    ; Load case '//StrId(mlc(i)))
        call INADD (sam, rhs(:,i), mlc(i), iprint, lpu, ierr)
        if (ierr < 0) goto 890
     end do

  end if

  if (neval > 0) then
     iopAdd = 3 ! Assemble both stiffness and mass
  else if (sam%mpar(30) > 0) then
     iopAdd = -3 ! Assemble stiffness and gravitation forces
  else
     iopAdd = 1 ! Assemble stiffness only
  end if

  ! --- Allocate substructure stiffness- and mass matrices

  call writeProgress (' --> Allocating substructure matrices')
  call csBeginAssembly (csStiffMat,'system stiffness matrix',lpu,ierr)
  if (ierr /= 0) goto 890
  if (iopAdd == 3) then
     call csBeginAssembly (csMassMat,'system mass matrix',lpu,ierr)
     if (ierr /= 0) goto 890
  end if

  ! --- Now assemble the substructure mass- and stiffness matrices

  if (sam%mpar(30) > 0) then

     ! Also assemble gravitation force vectors, {gx,gy,gz} = [M]*{ax,ay,az},
     ! where ax,ay,az are unit accelerations vectors in each global direction.
     ! Notice that we have to do this multiplication on the element level and
     ! assemble the resulting vectors, and not perform the matrix-vector product
     ! on the system level as in the reducer, due to the boundary conditions.
     allocate(Qg(sam%neq,3),STAT=ierr)
     if (ierr /= 0) then
        ierr = allocationError('Gravitation force vectors')
        return
     end if

     call INADD (sam, csStiffMat, csMassMat, sMass, rMass, &
          &      lumpedMass, .false., iopAdd, iprint, lpu, ierr, Qg)
  else
     call INADD (sam, csStiffMat, csMassMat, sMass, rMass, &
          &      lumpedMass, .false., iopAdd, iprint, lpu, ierr)
  end if
  if (ierr < 0) then
     call reportError (error_p,'Can not build substructure matrices')
     goto 890
  end if

  call ffl_done (.true.) ! Assembly finished, release the link object

  call csEndAssembly (csStiffMat,ierr)
  if (ierr /= 0) goto 900
  if (iopAdd == 3) then
     call csEndAssembly (csMassMat,ierr)
     if (ierr /= 0) goto 900
  end if

  ! --- Calculate gravity load

  if (sam%mpar(30) > 0) then

     ! Multiply with the specified gravitation constant
     ! and add to each right-hand-side vector
     do i = 1, nrhs
        call DGEMV ('N',sam%neq,3,1.0_dp,Qg(1,1),sam%neq, &
             &      grav(1),1,1.0_dp,rhs(1,i),1)
     end do
     deallocate(Qg)

  end if

  ! --- Store all SAM-arrays on file for fedem_stress

  if (ffa_cmdlinearg_isSet('samfile')) then
     ieq = sam%ndof1; sam%ndof1 = 0 ! Since meqn1 is not stored, ndof1 must be 0
     call getFileName ('samfile',chname,'_SAM.fsm')
     call saveSAM (chname,0,sam,ierr)
     if (ierr < 0) then
        call reportError (error_p,'Can not save SAM-data to file '//chname)
        goto 900
     end if
     sam%ndof1 = ieq
  end if

  call getFileName ('frsfile',chname)

  if (neval > 0) then

     ! --- Compute eigenvalues and eigenvectors of the substructure

     call ffa_cmdlinearg_getbool ('factorMass',factorMass)
     call ffa_cmdlinearg_getdouble ('tolEigval',tolEigval)
     call ffa_cmdlinearg_getdouble ('eigenshift',eigenShift)
     call EIGVAL (sam, csStiffMat, csMassMat, neval, factorMass, &
          &       iopSing, tolFactorize, tolEigval, eigenShift, &
          &       eval, evec, iprint, lpu, ierr)
     if (ierr < 0) then
        call reportError (error_p,'Can not perform eigenvalue calculation')
        goto 900
     end if

     write(lpu,6010) (i,eval(i)*0.5_dp/pi_p,i=1,neval)

     if (chname /= '') then
        ! Save the mode shapes to frs-file for plotting
        call writeProgress (' --> Saving results')
        call writeModes (chname,sam,sam%mpar(18),'Eigenmodes', &
             &           neval,eval,evec,ierr)
        if (ierr /= 0) then
           call reportError (warning_p,'Failed to write mode shapes file', &
                &                      'Execution continues...')
        end if
     end if

     deallocate(eval,evec)

  else if (nrhs > 0) then

     if (iprint > 2) then
        call writeObject (rhs,lpu,'Right-hand-side vector(s)')
     end if

     ! --- Direct solve for FE nodal displacements

     call writeProgress(' --> Solving FE equations')

     call csSolve (3,iopSing,csStiffMat,rhs,lpu,ierr, &
          &        eqnInErr,tolFactorize=tolFactorize)
     if (ierr == -3) then
        j = int(eqnInErr)
        i = NodeNumFromEqNum(sam,j)
        if (i > 0) then
           call reportError (error_p, &
                'Singularity encountered in node '//StrId(sam%minex(i)))
        else
           write(lpu,"(10X,'eqnInErr =',I8)") eqnInErr
        end if
     else if (ierr < 0) then
        call reportError (error_p,'Linear solve failed, ierr =',ierr=ierr)
     else

        if (iprint > 2) then
           call writeObject (rhs,lpu,'Solution vector(s)')
        end if

        ! Extract displacement values for control point output
        ! and print out the max value in each direction
        write(lpu,"()")
        displ = 0.0_dp
        do k = 1, 3
           nMax = 0
           dMax = 0.0_dp
           do i = 1, sam%nnod
              if (sam%madof(i)+k <= sam%madof(i+1)) then
                 ieq = sam%meqn(sam%madof(i)+k-1)
                 if (ieq > 0 .and. ieq <= sam%neq) then
                    if (nMax == 0 .or. abs(rhs(ieq,1)) > dMax) then
                       nMax = sam%minex(i)
                       dMax = rhs(ieq,1)
                    end if
                    do j = 1, n
                       if (nodes(j) == sam%minex(i)) displ(k,j) = rhs(ieq,1)
                    end do
                 end if
              end if
           end do
           write(lpu,600) char(ichar('W')+k),dMax,nMax
600        format('    Max ',A1,'-displacement:',1PE12.5,'  in node',I8)
        end do
        write(lpu,"()")

        if (chname /= '') then
           ! Save the displacement vector(s) to frs-file for plotting
           call writeProgress (' --> Saving results')
           if (iprint > 0) iprint = lpu
           call writeDeformation (chname,sam,sam%mpar(18),nrhs,mlc,rhs, &
                &                 iprint,ierr)
        end if

     end if

     deallocate(rhs)

  end if
  ierr = 0
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
  call deAllocateSysMatrix (csStiffMat,i)
  call deAllocateSysMatrix (csMassMat)
  call deAllocateSAM (sam)
  call castToInt8 (sam)
  call releaseScratchArrays
  call writeProgress ('     Finished')

  write(lterm,6000) ' END OF PROGRAM SOLVER '

  if (ierr == 0 .and. i /= 0) ierr = i
  if (ierr /= 0) then
     call reportError (debugFileOnly_p,'solver')
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
6010 format(/5X,'MODE  FREQUENCY [Hz]'/ 1P,(I8,E17.8))
6100 format(/11X,'An error is detected. Error messages are written to ',A)
6200 format(/ 4X,'Linear solver ',A)

end subroutine solver
