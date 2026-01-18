!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file inputReducerModule.f90
!> @brief Data input and preprocessing for fedem_reducer.

!!==============================================================================
!> @brief Module with subroutines for input and preprocessing of FE model data.

module InputReducerModule

  implicit none

contains

  !!============================================================================
  !> @brief Administers data input and preprocessing for fedem_reducer.
  !>
  !> @param[out] sam Data for managing system matrix assembly
  !> @param[out] sysK System stiffness matrix
  !> @param[out] sysM System mass matrix
  !> @param[out] mass Total mass of the FE part
  !> @param[out] rotMass Total rotational mass of the FE part
  !> @param[out] tenc Table of External Node Coordinates
  !> @param[out] medof Matrix of External DOFs
  !> @param[out] mlc Matrix of Load Case identifiers
  !> @param[in] noMass If .true., the system mass matrix is not needed
  !> @param[in] diagMass If .true., use a diagonalized system mass matrix
  !> @param[in] lumpedMass If .true., use a lumped mass matrix formulation
  !> @param[in] demoEdition If .true., a demo license is used
  !> @param[in] ipsw Print switch for debug output
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !> @param[in] directSolve If present, we are doing a direct linear solve
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 20 Sep 2000

  subroutine readReducerData (sam,sysK,sysM,mass,rotMass,tenc,medof,mlc, &
       &                      noMass,diagMass,lumpedMass,demoEdition, &
       &                      ipsw,lpu,ierr,directSolve)

    use SamModule              , only : SamType, dp
    use SamModule              , only : nullifySAM, writeSAMSize, writeObject
    use SamReducerModule       , only : initiateSAM
    use SysMatrixTypeModule    , only : SysMatrixType, diagonalMatrix_p
    use SysMatrixTypeModule    , only : skylineMatrix_p, sparseMatrix_p
    use SysMatrixTypeModule    , only : denseMatrix_p, outofCore_p
    use SysMatrixTypeModule    , only : nullifySysMatrix, shareMatrixStructure
    use SysMatrixTypeModule    , only : check, writeObject
    use AsmExtensionModule     , only : castToInt8, csAllocPointerArrays
    use AsmExtensionModule     , only : csRenumber, csPreassemble
    use ManipMatrixModule      , only : writeObject
    use AllocationModule       , only : reAllocate
    use ProgressModule         , only : writeProgress
    use ReportErrorModule      , only : reportError
    use ReportErrorModule      , only : note_p, error_p, debugFileOnly_p
    use FFlLinkHandlerInterface, only : ffl_init, ffl_done
    use FFlLinkHandlerInterface, only : ffl_getLoadCases, ffl_massProp
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_getint
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_getbool
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_isSet
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_isTrue
    use FFaProfilerInterface   , only : ffa_getPhysMem

    type(SamType)      , intent(out) :: sam
    type(SysMatrixType), intent(out) :: sysK, sysM
    real(dp)           , intent(out) :: mass, rotMass
    real(dp)           , pointer     :: tenc(:,:)
    integer            , pointer     :: medof(:)
    logical            , intent(in)  :: noMass, diagMass,lumpedMass, demoEdition
    integer            , intent(in)  :: ipsw, lpu
    integer            , intent(out) :: mlc(:), ierr
    integer,  optional , intent(in)  :: directSolve

    !! Local variables
    logical            :: linearSolve, factorMass
    integer            :: i, nenod, neval, nlc, ngen, ndim, matrixType, nFreeMem
    real(dp)           :: cg(3), inertia(6)
    real(dp), pointer  :: tncoor(:,:)
    character(len=128) :: errMsg

    !! --- Logic section ---

    nullify(tncoor)
    call ffl_init (0,0,ierr)
    if (ierr /= 0) then
       call nullifySAM (sam)
       call reportError (error_p,'Failure parsing FE data file', &
            &            addString='ffl_reducer_init')
       goto 990
    end if

    linearSolve = present(directSolve)

    ! --- Initialize the SAM datastructure
    call initiateSAM (sam,tncoor,linearSolve,ipsw,lpu,ierr)
    if (ierr /= 0) goto 990

    call ffa_cmdlinearg_getint ('linkId',sam%mpar(18))
    call ffa_cmdlinearg_getint ('ngen',ngen)
    if (ngen > 0) then
       call ffa_cmdlinearg_getint ('neval',neval)
       call ffa_cmdlinearg_getbool ('factorMass',factorMass)
       sam%mpar(19) = max(neval,ngen)
    else
       factorMass = .false.
    end if

    ! --- Compute the number of external nodes and the reduced matrix dimension
    nenod = 0
    ndim = ngen
    do i = 1, sam%nnod
       if (sam%mnnn(i) > 1) then
          nenod = nenod + 1
          ndim = ndim + sam%madof(i+1)-sam%madof(i)
       end if
    end do

    ! --- Get the number of load cases
    nlc = size(mlc)
    call ffl_getloadcases (mlc(1),nlc)

    sam%mpar(21) = nlc
    sam%mpar(22) = ngen
    sam%mpar(24) = ndim
    sam%mpar(29) = nenod
    if (present(directSolve)) then
       sam%mpar(30) = directSolve
       if (nlc < 1) mlc(1) = 1
    end if

    write(lpu,600) sam%nel,sam%nnod,sam%mpar(29),sam%mpar(22), &
         &         sam%mpar(17),sam%mpar(27),sam%mpar(21),sam%mpar(24)

    ! --- Compute and report mass properties of the unreduced FE model
    call writeProgress (' --> Computing mass properties')
    call ffl_massprop (mass,cg,inertia)
    write(lpu,610) mass,cg,inertia(1),inertia(4:5), &
         &         inertia(2),inertia(6),inertia(3)
    call flush(lpu)

    rotMass = sum(inertia(1:3)) / 3.0_dp

600 format(//4X,'FE MODEL SUMMARY' / 4X,16('-') &
         & //4X,'Number of elements             =',I8 &
         &  /4X,'Number of nodes (total)        =',I8 &
         &  /4X,'Number of external nodes       =',I8 &
         &  /4X,'Number of generalized modes    =',I8 &
         &  /4X,'Number of materials            =',I8 &
         &  /4X,'Number of geometric properties =',I8 &
         &  /4X,'Number of load cases           =',I8 &
         &  /4X,'Dimension of reduced matrices  =',I8 / )
610 format(  4X,'Total mass                     =',1PE13.5 &
         &  /4X,'Center of gravity              =',3E13.5 &
         &  /4X,'Moments of inertia (tensor)    =',3E13.5 &
         &  /49X,2E13.5 /38X,'symm.',19X,E13.5 / )

    if (.not. (linearSolve .or. noMass)) then

       ! --- Store away nodal coordinates of the external nodes (needed in JCMS)
       call reAllocate ('readReducerData',tenc,3,nenod,ierr)
       call reAllocate ('readReducerData',medof,nenod+1,ierr)
       if (ierr < 0) goto 990

       nenod = 0
       medof(1) = 1
       do i = 1, sam%nnod
          if (sam%mnnn(i) > 1) then
             nenod = nenod + 1
             tenc(:,nenod) = tncoor(i,:)
             medof(nenod+1) = medof(nenod) + sam%madof(i+1)-sam%madof(i)
          end if
       end do

    end if

    ! --- Check if any of the nodes lack stiffness or boundary conditions
    call writeProgress (' --> Checking model consistency')
    call CheckSingularNodes (sam%nel,sam%nnod,sam%minex, &
         &                   sam%mpmnpc,sam%mmnpc,sam%madof,sam%msc,ierr)
    if (ierr > 0) then
       if (ipsw > 0) then
          call writeObject(sam%minex,lpu,'MINEX')
          call writeObject(sam%mpmnpc,lpu,'MPMNPC')
          call writeObject(sam%mmnpc,lpu,'MMNPC')
          call writeObject(sam%madof,lpu,'MADOF')
          call writeObject(sam%msc,lpu,'MSC')
       end if
       goto 990
    end if

    if (associated(sam%mmceq)) then

       ! --- Check linear couplings
       call CheckLinearCouplings (sam%mpmceq,sam%mmceq,sam%msc,ierr)
       if (ierr > 0) then
          if (ipsw > 0) then
             call writeObject(sam%mpmceq,lpu,'MPMCEQ')
             call writeObject(sam%mmceq,lpu,'MMCEQ')
             call writeObject(sam%msc,lpu,'MSC')
          end if
          goto 990
       end if
    end if

    ! --- Allocate the pointer structure for the linear equation solver
    call writeProgress (' --> Creating data structures for equation solver')
    if (ffa_cmdlinearg_isTrue('skylineSolver')) then
       matrixType = skylineMatrix_p
    else if (ffa_cmdlinearg_isTrue('denseSolver')) then
       matrixType = denseMatrix_p
    else
       call ffa_cmdlinearg_getint ('gsfSolver',matrixType)
       if (demoEdition .or. nenod == sam%nnod .or. matrixType < 1) then
          matrixType = sparseMatrix_p
       else if (matrixType == 1) then
          matrixType = outOfCore_p + sparseMatrix_p
       else
          matrixType = outOfCore_p
       end if
    end if

    ngen = 0
    if (.not.demoEdition .and. ffa_cmdlinearg_isSet('cachesize')) then
       call ffa_cmdlinearg_getint ('cachesize',ngen)
    else if (matrixType >= outOfCore_p) then
       nFreeMem = ffa_getPhysMem(.false.)
       if (nFreeMem > 0) then
          ngen = min(32,nFreeMem/4)
          write(errmsg,650) ngen, 100*ngen/nFreeMem
650       format('Setting GSF solver RAM size to',I4,'MB (',I3, &
               & '% of available memory).')
          call reportError (note_p,errmsg,'Use option -cachesize to override.')
       end if
    end if

    if (matrixType == sparseMatrix_p) then
       call castToInt8 (sam,ierr)
       if (ierr < 0) goto 990
    else if (matrixType >= outOfCore_p) then
       write(lpu,"(/' BEGIN EQUATION SOLVER INFO')")
    end if

    call csAllocPointerArrays (sam,sysK,matrixType,lpu,ierr,ipsw, &
         &                     cacheSize=ngen)
    if (ierr < 0) goto 991

    ! --- Optimize the equation system
    call csRenumber (sam,sysK,lpu,ierr)
    if (ierr < 0) goto 991

    ! --- Preassembly
    call writeProgress (' --> Pre-assembling the stiffness matrix')
    call csPreassemble (sam,sysK,.not.linearSolve,lpu,ierr,ipsw)
    if (ierr < 0) goto 991

    ! --- Check for consistency
    call check (sysK,sam%mpar,lpu,ierr)
    if (ierr < 0) goto 991

    ! --- Create the mass matrix data structure

    if (noMass .or. (linearSolve .and. max(sam%mpar(19),sam%mpar(30)) < 1)) then

       ! No system mass matrix needed
       call nullifySysMatrix (sysM)

    else if (.not.diagMass .and. .not.lumpedMass) then

       ! Let the mass- and stiffness matrices share datastructures
       call shareMatrixStructure (sysM,sysK,factorMass,ierr=ierr)
       if (ierr < 0) goto 900

    else if (sam%nmmceq > 0 .and. lumpedMass) then

       ! Lumped mass is requested for a FE model that has constraint equations.
       ! The system mass matrix will then have a different sparsity pattern
       ! than the system stiffness matrix.

       if (matrixType == outOfCore_p+sparseMatrix_p) then
          ! Use the SPR solver for the mass matrix
          call castToInt8 (sam,ierr)
          if (ierr < 0) goto 990
          call csAllocPointerArrays (sam,sysM,-sparseMatrix_p,lpu,ierr,ipsw)
       else
          ! Use the same solver for the lumped mass as used for the stiffness
          call csAllocPointerArrays (sam,sysM,-sysK%storageType,lpu,ierr, &
               &                     ipsw,ngen,factorMass,sysK%gsf)
       end if
       if (ierr < 0) goto 900

       ! --- Optimize the equation system
       call csRenumber (sam,sysM,lpu,ierr)
       if (ierr < 0) goto 900

       ! --- Preassembly
       call writeProgress (' --> Pre-assembling the lumped mass matrix')
       if (matrixType == outOfCore_p+sparseMatrix_p .and. sam%mpar(22) > 0) then
          call csPreassemble (sam,sysM,.false.,lpu,ierr,ipsw,sysK%meqn)
       else
          call csPreassemble (sam,sysM,.false.,lpu,ierr,ipsw,sam%meqn)
       end if
       if (ierr < 0) goto 900

       ! --- Check for consistency
       call check (sysM,sam%mpar,lpu,ierr)
       if (ierr < 0) goto 900

    else

       ! Lumped mass is requested for a FE model with no constraint equations.
       ! The system mass matrix will then be a pure diagonal matrix.
       call shareMatrixStructure (sysM,sysK,factorMass,diagonalMatrix_p)

    end if

    if (matrixType == outOfCore_p+sparseMatrix_p) then
       ! Reset the stiffness matrix equation numbers to the external order
       call reAllocate ('inputReducer',sysK%meqn)
       sysK%meqn => sam%meqn
       ! Reset the mass matrix equation numbers also if sharing datastructures
       if (sysM%isShared) sysM%meqn => sam%meqn
    end if

    if (matrixType >= outOfCore_p) write(lpu,"(' END EQUATION SOLVER INFO'/)")

    if (ipsw > 0) then

       ! --- Print out substructure data to res-file
       call SUBPRIN (sam%mpar(18),sam%mpar,sam%mpmnpc,sam%mmnpc,sam%melcon, &
            &        sam%minex,tncoor(:,1),tncoor(:,2),tncoor(:,3),sam%madof, &
            &        sam%msc,sam%mpmceq,sam%mmceq,sam%ttcc,sam%mnnn,LPU)

    end if

    call writeSAMsize (sam,lpu)

    if (sam%ndof2 < 1 .and. .not.linearSolve) then
       ierr = -1
       call reportError (error_p,'This FE part has no external DOFs.', &
            'Check that it is properly attached to other model members.')
    else if (sam%mpar(22) > sam%ndof1) then
       ierr = -1
       write(errMsg,690) sam%mpar(22), sam%ndof1
690    format('NGEN =',I8,',  NDOF1 =',I8)
       call reportError (error_p,'The number of generalized DOFs (NGEN)', &
            'is larger than the number of internal nodal DOFs (NDOF1).', &
            errMsg,'This is not allowed, reduce NGEN to less than NDOF1.')
    end if

900 continue
    if (ipsw > 1 .and. ierr <= 0) then
       call writeObject (sam,lpu,4)
       call writeObject (sysK,sam%mpar,lpu,' === System stiffness matrix',10,3)
       call writeObject (sysM,sam%mpar,lpu,' === System mass matrix',10,3)
    else if (ipsw >= -7 .and. ipsw <= -5 .and. ierr == 0) then
       call writeObject (sam,lpu,5)
    end if
    call reAllocate ('readReducerData',tncoor)
    if (ierr == 0) return
    call ffl_done (.true.)
    call reportError (debugFileOnly_p,'readReducerData')
    return

990 call nullifySysMatrix (sysK)
991 call nullifySysMatrix (sysM)
    goto 900

  end subroutine readReducerData


  !!============================================================================
  !> @brief Checks the DOF status correctness for all linear couplings.
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen

  subroutine checkLinearCouplings (mpmceq,mmceq,msc,ierr)

    use ReportErrorModule, only : reportError, error_p, errorFileOnly_p

    integer, intent(in)  :: mpmceq(:), mmceq(:), msc(:)
    integer, intent(out) :: ierr

    !! Local variables
    integer            :: i, ieq, iStart, iEnd, inDep, iDep, nceq
    character(len=128) :: errMsg

    !! --- Logic section ---

    ierr = 0
    nceq = size(mpmceq)-1

    do ieq = 1, nceq
       iStart = mpmceq(ieq)
       iEnd   = mpmceq(ieq+1)-1
       iDep   = mmceq(iStart)
       if (msc(iDep) == 1 .or. msc(iDep) == 2) then
          ierr = ierr + 1
          write(errMsg,600) msc(iDep),'dependent',iDep,ieq
          call reportError (errorFileOnly_p,errMsg)
       end if
       do i = iStart+1, iEnd
          inDep = mmceq(i)
          if (msc(inDep) /= 1 .and. msc(inDep) /= 2) then
             ierr = ierr + 1
             write(errMsg,600) msc(iDep),'independent',inDep,ieq
             call reportError (errorFileOnly_p,errMsg)
          end if
       end do
    end do

    if (ierr > 0) then
       write(errMsg,690) ierr
       call reportError (error_p,errMsg,addString='checkLinearCouplings')
    end if

600 format('Invalid status code',I2,' for ',A,' DOF',I8, &
         & ' in constraint equation',I8)
690 format('Found',I4,' errors associated with linear couplings')

  end subroutine checkLinearCouplings


  !!============================================================================
  !> @brief Checks that all nodes are connected to an element.
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen

  subroutine checkSingularNodes (nel,nnod,minex,mpmnpc,mmnpc,madof,msc,numSing)

    use ScratchArrayModule, only : getIntegerScratchArray
    use ReportErrorModule , only : reportError, error_p, errorFileOnly_p

    integer, intent(in)  :: nel, nnod
    integer, intent(in)  :: minex(:), mpmnpc(:), mmnpc(:), madof(:), msc(:)
    integer, intent(out) :: numSing

    !! Local variables
    integer            :: iel, inod, ip
    integer, pointer   :: nConnected(:)
    character(len=128) :: errMsg

    !! --- Logic section ---

    nConnected => getIntegerScratchArray(nnod,numSing)
    if (numSing < 0) return

    !! Flag all nodes with specified DOFs as connected,
    !! and nodes with one or more free DOFs as unconnected
    do inod = 1, nnod
       nConnected(inod) = 1 ! assume the node is constrained
       do ip = madof(inod), madof(inod+1)-1
          if (msc(ip) == 1 .or. msc(ip) == 2) nConnected(inod) = 0
       end do
    end do

    !! Flag all element nodes as connected
    do iel = 1, nel
       do ip = mpmnpc(iel), mpmnpc(iel+1)-1
          inod = mmnpc(ip)
          nConnected(inod) = nConnected(inod) + 1
          inod = -minex(inod)
          if (inod > 0 .and. inod <= nnod) then
             !! This element node is coupled to another node via beam pin
             !! constraint equations. Flag that the associated independent node
             !! is connected in case no other elements uses that node.
             nConnected(inod) = nConnected(inod) + 1
          end if
       end do
    end do

    !! Report all nodes that do not have stiffness or boundary conditions
    numSing = 0
    do inod = 1, nnod
       if (nConnected(inod) < 1) then
          write(errMsg,600) minex(inod)
          if (numSing < 10) then
             call reportError (error_p,errMsg)
          else
             call reportError (errorFileOnly_p,errMsg)
          end if
          numSing = numSing + 1
       end if
    end do

    if (numSing > 0) then
       write(errMsg,690) numSing
       call reportError (error_p,errMsg,addString='checkSingularNodes')
    end if

600 format('Node',I8,' has DOFs without stiffness or constraints')
690 format('Found',I4,' nodes without stiffness or boundary conditions')

  end subroutine checkSingularNodes


  !!============================================================================
  !> @brief Prints out data for a substructure.
  !>
  !> @callergraph
  !>
  !> @author Ketil Aamnes
  !> @date 14 Oct 1988
  !>
  !> @author Knut Morten Okstad
  !> @date 10 Sep 2002

  subroutine SUBPRIN (LINK,MPAR,MPMNPC,MNPC,MEKN,MINEX, &
       &              X,Y,Z,MADOF,MSC,MPMCEQ,MCEQ,COEFF,INODE,LPU)

    use KindModule             , only : dp
    use FFlLinkHandlerInterface, only : ffl_getElmId

    integer , intent(in) :: LINK, MPAR(:), MPMNPC(:), MNPC(:), MEKN(:), MINEX(:)
    integer , intent(in) :: MADOF(:), MSC(:), MPMCEQ(:), MCEQ(:), INODE(:), LPU
    real(dp), intent(in) :: COEFF(:), X(:), Y(:), Z(:)

    !! Local variables
    logical :: lFirst
    integer :: I,ICOMPM,ICOMPS,J,INODM,INODS,NMST

    !! --- Logic section ---

    write(LPU,600) LINK,MPAR(1:13)
    write(LPU,601) MPAR(14:28)

    ! --- Print out element data

    write(LPU,610)
    do I = 1, MPAR(2)
       if (MPMNPC(I+1)-MPMNPC(I) <= 10) then
          write(LPU,611) I,ffl_getElmId(I),MEKN(I),MNPC(MPMNPC(I):MPMNPC(I+1)-1)
       else
          write(LPU,612) I,ffl_getElmId(I),MEKN(I),MNPC(MPMNPC(I):MPMNPC(I+1)-1)
       end if
    end do

    ! --- Print out nodal data

    write(LPU,630)
    do I = 1, MPAR(1)
       write(LPU,631) I,MINEX(I),INODE(I),X(I),Y(I),Z(I), &
            &         MSC(MADOF(I):MADOF(I+1)-1)
    end do

    ! --- Print out prescribed displacements

    lFirst = .true.
    do I = 1, MPAR(7)
       if (MPMCEQ(I+1)-MPMCEQ(I) == 1) then
          if (lFirst) write(LPU,640)
          lFirst = .false.
          ICOMPS = MCEQ(MPMCEQ(I))
          INODS  = nodeNumFromEqNum(ICOMPS)
          write(LPU,641) INODS,ICOMPS,COEFF(MPMCEQ(I))
       end if
    end do

    ! --- Print out dependent displacements

    lFirst = .true.
    do I = 1, MPAR(7)
       NMST = MPMCEQ(I+1)-MPMCEQ(I) - 1
       if (NMST > 0) then
          if (lFirst) write(LPU,650)
          lFirst = .false.
          ICOMPS = MCEQ(MPMCEQ(I))
          ICOMPM = MCEQ(MPMCEQ(I)+1)
          INODS  = nodeNumFromEqNum(ICOMPS)
          INODM  = nodeNumFromEqNum(ICOMPM)
          write(LPU,651) INODS,INODM,ICOMPS,ICOMPM,COEFF(MPMCEQ(I):MPMCEQ(I)+1)
          do J = 2, NMST
             ICOMPM = MCEQ(MPMCEQ(I)+J)
             INODM  = nodeNumFromEqNum(ICOMPM)
             write(LPU,652) INODM,ICOMPM,COEFF(MPMCEQ(I)+J)
          end do
       end if
    end do

600 format(/4X,'SUBSTRUCTURE TYPE NUMBER :',I8 /4X,34('=') &
         &//4X,' 1. NNOD   =',I10 &
         & /4X,' 2. NEL    =',I10 &
         & /4X,' 3. NDOF   =',I10 &
         & /4X,' 4. NDOF1  =',I10 &
         & /4X,' 5. NDOF2  =',I10 &
         & /4X,' 6. NSPDOF =',I10 &
         & /4X,' 7. NCEQ   =',I10 &
         & /4X,' 8. NSDOF  =',I10 &
         & /4X,' 9. NPDOF  =',I10 &
         & /4X,'10. NDDOF  =',I10 &
         & /4X,'11. NEQ    =',I10 &
         & /4X,'12. NSKY   =',I10 &
         & /4X,'13. LBW    =',I10 )
601 format( 4X,'14. NODLBW =',I10 &
         & /4X,'15. NEMNPC =',I10 &
         & /4X,'16. NEMCEQ =',I10 &
         & /4X,'17. NEGR   =',I10 &
         & /4X,'18. linkId =',I10 &
         & /4X,'19. NEVAL  =',I10 &
         & /4X,'20. MAXNIV =',I10 &
         & /4X,'21. MAXIT  =',I10 &
         & /4X,'22. NGEN   =',I10 &
         & /4X,'23. NBEL   =',I10 &
         & /4X,'24. NDIM   =',I10 &
         & /4X,'25. IPRIN1 =',I10 &
         & /4X,'26. IPRIN2 =',I10 &
         & /4X,'27. NPROP  =',I10 &
         & /4X,'28. NLMPC  =',I10 )

610 format(/4X,'ELEMENT DATA' &
         & /4X,'------------' &
         //'  ELEMENT NUMBER   ELEMENT' &
         &/' INTERNAL EXTERNAL PROPERTY   INTERNAL NODE NUMBERS' /)
611 format(2I9,I5,2X,10I8 )
612 format(2I9,I5,2X,10I8 / (25X,10I8) )

630 format(/4X,'NODE DATA' &
         & /4X,'---------' &
         &//4X,'NODE NUMBER' &
         &/' INTERNAL EXTERNAL  IST',6X,'XCOOR',9X,'YCOOR',9X,'ZCOOR', &
         & '     IX  IY  IZ  IRX IRY IRZ'/)
631 format(2I9,I4,1X,1P3E14.6,6I4)

640 format(/4X,'PRESCRIBED DISPLACEMENTS' &
         & /4X,'------------------------' &
         &//4X,'NODE    COMPONENT    DISPLACEMENT' &
         & /4X,'NUMBER'/)
641 format(2I9,5X,1PE12.3)

650 format(/4X,'DEPENDENT DISPLACEMENTS' &
         & /4X,'-----------------------')
651 format(//4X,'DEPENDENT NODE :',I6,    11X,'INDEPENDENT NODE :',I6, &
         &  /4X,'COMPONENT      :',I6,    11X,'COMPONENT        :',I6, &
         &  /4X,'COEFF          :',1PE12.3,5X,'COEFFICIENT      :',1PE12.3)
652 format(/37X,'INDEPENDENT NODE :',I6 &
         & /37X,'COMPONENT   :',I6 &
         & /37X,'COEFFICIENT :',1PE12.3)

  contains

    !> @brief Returns the node number and local DOF for a global DOF.
    function nodeNumFromEqNum (IDOF) result(INOD)
      integer, intent(inout) :: IDOF
      integer                :: INOD
      do INOD = 1, size(MADOF)-1
         if (IDOF >= MADOF(INOD) .and. IDOF < MADOF(INOD+1)) then
            IDOF = IDOF+1 - MADOF(INOD)
            return
         end if
      end do
    end function nodeNumFromEqNum

  end subroutine SUBPRIN

end module InputReducerModule
