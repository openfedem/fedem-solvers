!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file saveReducerModule.f90
!>
!> @brief Save of FE model reduction data to result database files.

!!==============================================================================
!> @brief Module with subroutines for saving of FE model reduction results.
!>
!> @details This module contains a collection of subroutines for writing results
!> data from the FE model reduction process or direct linear solve
!> to the results database.

module saveReducerModule

  implicit none

  private

  integer, save :: idVal, idVec(2)
  integer, save :: idDis, id3Dof

  public :: writeModes, writeDeformation


contains

  !!============================================================================
  !> @brief Administers writing of results database for fedem_reducer.
  !>
  !> @param[in] chname Name of frs-file receiving the results data
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[in] partId Base ID of the FE part
  !> @param[in] modeType Name of eigenmode item group
  !> @param[in] nMode Number of eigenmodes
  !> @param[in] ev Eigenvalues
  !> @param[in] sv Eigenvectors
  !> @param[out] ierr Error flag
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 4 Oct 2002

  subroutine writeModes (chname,sam,partId,modeType,nMode,ev,sv,ierr)

    use KindModule            , only : dp, i8, lfnam_p
    use RDBModule             , only : RDBType, idatd, openRDBfile, closeRDBfile
    use RDBModule             , only : openHeaderFiles, nullifyRDB, writeRDB
    use RDBModule             , only : writeTimeStepHeader, writeTimeStepDB
    use SamModule             , only : SamType
    use IdTypeModule          , only : IdType, nullifyId, writeIdHeader
    use ReportErrorModule     , only : reportError, debugFileOnly_p
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getstring

    character*(*), intent(in)  :: chname, modeType
    type(SamType), intent(in)  :: sam
    integer      , intent(in)  :: partId, nMode
    real(dp)     , intent(in)  :: ev(:), sv(:,:)
    integer      , intent(out) :: ierr

    !! Local variables
    integer, parameter :: nBits = 32
    type(RDBType)      :: rdb
    type(IdType)       :: supId
    character(lfnam_p) :: fileName

    !! --- Logic section ---

    idVal = 0
    idVec = 0
    call nullifyRDB (rdb)
    call nullifyId (supId)
    supId%baseId = partId

    !! Open the temporary header files

    call ffa_cmdlinearg_getstring ('linkfile',fileName)
    call openHeaderFiles ('fedem_reducer',ierr, &
         &                'modes data base file',mName=fileName)
    if (ierr /= 0) goto 190

    !! Write modes results definitions

    call writeTimeStepHeader (rdb)
    call writeIdHeader ('Part',supId,idatd,.true.)
    write(idatd,610) modeType
    call writeModesHeader (rdb,nMode,nBits)
    call writeVectorHeader (rdb,sam,nMode,nBits,ierr)
    if (ierr /= 0) goto 190
    write(idatd,"('  ]'/'}')")

    !! Open the results database file and copy the file header to it
    call openRDBfile (rdb,ierr,chname,'#FEDEM modal data')
    if (ierr < 0) goto 190

    !! Write time step identification
    call writeTimeStepDB (rdb,0_i8,0.0_dp,ierr)
    if (ierr < 0) goto 190

    !! Write eigenvalues
    call writeRDB (rdb,ev(1:nMode),ierr)
    if (ierr < 0) goto 190

    !! Write mode shape data
    call writeDisplacementDB (rdb,sam,nMode,sv,ierr=ierr)
    if (ierr < 0) goto 190

    !! Close results data base file
    ierr = 0
    call closeRDBfile (rdb,ierr)
    if (ierr == 0) return

190 call reportError (debugFileOnly_p,'writeModes')

610 format('  [;"',a,'";')

  end subroutine writeModes


  !!============================================================================
  !> @brief Administers writing of results database for fedem_partsol.
  !>
  !> @param[in] chname Name of frs-file receiving the results data
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[in] partId Base ID of the FE part
  !> @param[in] nlc Number of load cases
  !> @param[in] mlc Load case identifiers
  !> @param[in] sv Displacement vector in equation order
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 19 Mar 2018

  subroutine writeDeformation (chname,sam,partId,nlc,mlc,sv,lpu,ierr)

    use KindModule            , only : dp, lfnam_p
    use SamModule             , only : SamType
    use RDBModule             , only : RDBType, idatd
    use RDBModule             , only : openRDBfile, closeRDBfile, nullifyRDB
    use RDBModule             , only : openHeaderFiles, writeTimeStepHeader
    use IdTypeModule          , only : IdType, nullifyId, writeIdHeader
    use ReportErrorModule     , only : reportError, debugFileOnly_p
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getstring

    character*(*), intent(in)  :: chname
    type(SamType), intent(in)  :: sam
    integer      , intent(in)  :: partId, nlc, mlc(:), lpu
    real(dp)     , intent(in)  :: sv(:,:)
    integer      , intent(out) :: ierr

    !! Local variables
    integer, parameter :: nBits = 32
    type(RDBType)      :: rdb
    type(IdType)       :: supId
    character(lfnam_p) :: fileName

    !! --- Logic section ---

    idVec = 0
    idDis = 0
    id3Dof = 0
    call nullifyRDB (rdb)
    call nullifyId (supId)
    supId%baseId = partId

    !! Open the temporary header files

    call ffa_cmdlinearg_getstring ('linkfile',fileName)
    call openHeaderFiles ('fedem_reducer',ierr, &
         &                'response data base file',mName=fileName)
    if (ierr /= 0) goto 190

    !! Write displacement results definitions

    call writeTimeStepHeader (rdb)
    call writeIdHeader ('Part',supId,idatd,.true.)
    call writeVectorHeader (rdb,sam,-1,2*nBits)
    call writeNodesHeader (rdb,sam,nBits)
    write(idatd,"('}')")

    !! Open the results database file and copy the file header to it
    call openRDBfile (rdb,ierr,chname,'#FEDEM response data')
    if (ierr < 0) goto 190

    !! Write displacement vectors
    if (lpu > 0) then
       call writeDisplacementDB (rdb,sam,nlc,sv,mlc,lpu,ierr)
    else
       call writeDisplacementDB (rdb,sam,nlc,sv,mlc,ierr=ierr)
    end if
    if (ierr < 0) goto 190

    !! Close results data base file
    ierr = 0
    call closeRDBfile (rdb,ierr)
    if (ierr == 0) return

190 call reportError (debugFileOnly_p,'writeDeformation')

  end subroutine writeDeformation


  !!============================================================================
  !> @brief Writes modal result definitions to the temporary header files.
  !>
  !> @param rdb Results database file for modes output
  !> @param[in] nMode Number of mode shapes
  !> @param[in] nbit Number of bits per real word
  !>
  !> @author Knut Morten Okstad
  !>
  !> @callergraph
  !>
  !> @date 4 Oct 2002

  subroutine writeModesHeader (rdb,nMode,nbit)

    use RDBModule, only : rDBType, idatd, ivard

    type(RDBType), intent(inout) :: rdb
    integer      , intent(in)    :: nMode, nbit

    !! Local variables
    integer :: iMode

    !! --- Logic section ---

    write(idatd,610)
    if (idVal == 0) then
       rdb%nVar = rdb%nVar + 1
       idVal = rdb%nVar
       write(ivard,611) rdb%nVar,'Eigenvalue','ANGLE/TIME',nbit
    end if
    do iMode = 1, nMode
       write(idatd,612) iMode,idVal
    end do
    write(idatd,"('    ]')")
    rdb%nBytes = rdb%nBytes + 3*nMode*nbit/8

610 format('    [;"Eigenvalues";')
611 format('<',i1,';"',a,'";',a,';FLOAT;',i2,';SCALAR>')
612 format('      [;"Mode',i3,'";<',i1,'>]')

  end subroutine writeModesHeader


  !!============================================================================
  !> @brief Writes displacement vector definitions to temporary header files.
  !>
  !> @param rdb Results database file for displacement output
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[in] nMode Number of mode shapes/flag (see below)
  !> @param[in] nbit Number of bits per real word
  !> @param[out] ierr Error flag
  !>
  !> @details If @a nMode &gt; 0, we are writing eigenmode shape vectors
  !> and its value indicates the number of mode shapes to write.
  !> If @a nMode equals zero, we are writing the displacement response for
  !> visualization use (translational DOF components only). If @a nMode &lt; 0
  !> the entire displacement vector is stored (including rotational DOFs).
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 4 Oct 2002

  subroutine writeVectorHeader (rdb,sam,nMode,nbit,ierr)

    use RDBModule        , only : RDBType, idatd, ivard
    use SamModule        , only : SamType
    use ReportErrorModule, only : reportError, warning_p

    type(RDBType)    , intent(inout) :: rdb
    type(SamType)    , intent(in)    :: sam
    integer          , intent(in)    :: nMode, nbit
    integer, optional, intent(out)   :: ierr

    !! Local variables
    integer :: i, j, ndof

    !! --- Logic section ---

    j = 1
    ndof = 3*sam%nnod
    if (nMode < 0) then
       j = 2 ! Write the whole displacement vector
       ndof = sam%ndof
    else if (present(ierr)) then
       j = 1
       ierr = 0 ! Check that all nodes have (at least) 3 DOFs
       do i = 1, sam%nnod
          if (sam%madof(i+1) < sam%madof(i)+3) ierr = ierr + 1
       end do
       if (ierr > 0) then
          call reportError (warning_p,'Not all nodes have (at least) 3 DOFs')
          return
       end if
    end if

    if (idVec(j) == 0) then
       rdb%nVar = rdb%nVar + 1
       idVec(j) = rdb%nVar
       if (nMode < 0) then
          write(ivard,610) rdb%nVar,'Displacement','LENGTH',nbit,ndof
       else
          write(ivard,610) rdb%nVar,'Translational deformation', &
               &           'LENGTH',nbit,ndof
       end if
    end if

    if (nMode > 0) then
       write(idatd,611)
       do i = 1, nMode
          write(idatd,612) i,idVec(j)
       end do
       write(idatd,"('    ]')")
    else
       write(idatd,613) idVec(j)
    end if

    rdb%nBytes = rdb%nBytes + ndof*nbit/8

610 format('<',i1,';"',a,'";',a,';FLOAT;',i2,';VECTOR;(',i8,')>')
611 format('    [;"Vectors";')
612 format('      [;"Mode',i3,'";<',i1,'>]')
613 format('  [;"Vectors";'/'    [;"Dynamic response";<',i1,'>]'/'  ]')

  end subroutine writeVectorHeader


  !!============================================================================
  !> @brief Writes nodal translation definitions to the temporary header files.
  !>
  !> @param rdb Results database file for displacement output
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[in] nbit Number of bits per real word
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 4 Oct 2002

  subroutine writeNodesHeader (rdb,sam,nbit)

    use RDBModule, only : RDBType, idatd, iitem, ivard
    use SamModule, only : SamType

    type(RDBType), intent(inout) :: rdb
    type(SamType), intent(in)    :: sam
    integer      , intent(in)    :: nbit

    !! Local variables
    integer :: i, ndof

    !! --- Logic section ---

    ndof = 0

    write(idatd,620)
    do i = 1, sam%nnod
       !! Check that the node has (at least) 3 DOFs
       if (sam%madof(i+1) < sam%madof(i)+3) cycle
       if (idDis == 0) then
          rdb%nVar = rdb%nVar + 1
          idDis = rdb%nVar
          write(ivard,621) rdb%nVar,'Translational deformation','LENGTH',nbit
       end if
       if (id3DOF == 0) then
          rdb%nIG = rdb%nIG + 1
          id3DOF  = rdb%nIG
          write(iitem,622) id3DOF,idDis
       end if
       write(idatd,623) sam%minex(i),id3DOF
       ndof = ndof + 3
    end do
    write(idatd,"('  ]')")

    rdb%nBytes = rdb%nBytes + ndof*nbit/8

620 format('  [;"Nodes";')
621 format('<',i1,';"',a,'";',a,';FLOAT;',i2,';VEC3;(3);(("d_x","d_y","d_z"))>')
622 format('[',i1,';"Dynamic response";<',i1,'>]')
623 format('    [;',i8,';[',i1,']]')

  end subroutine writeNodesHeader


  !!============================================================================
  !> @brief Writes nodal displacements to the results database file.
  !>
  !> @param rdb Results database file for displacement output
  !> @param[in] sam Data for managing system matrix assembly
  !> @param[in] nComp Number of load cases or mode shapes
  !> @param[in] sv Displacement vector in equation order
  !> @param[in] mlc Load case identifiers
  !> @param[in] lpu File unit number for res-file output
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 4 Oct 2002

  subroutine writeDisplacementDB (rdb,sam,nComp,sv,mlc,lpu,ierr)

    use KindModule        , only : dp, i8, epsDiv0_p
    use RDBModule         , only : RDBType, writeRDB, writeTimeStepDB
    use SamModule         , only : SamType
    use SolExtensionModule, only : csExpand
    use ScratchArrayModule, only : getRealScratchArray
    use ReportErrorModule , only : internalError, reportError, debugFileOnly_p

    type(RDBType)    , intent(inout) :: rdb
    type(SamType)    , intent(in)    :: sam
    integer          , intent(in)    :: nComp
    real(dp)         , intent(in)    :: sv(:,:)
    integer, optional, intent(in)    :: mlc(:), lpu
    integer          , intent(out)   :: ierr

    !! Local variables
    integer           :: i, j, iCmp, nDim
    real(dp)          :: maxN
    real(dp), pointer :: rscr(:), sveq(:), svdof(:)

    !! --- Logic section ---

    nDim = size(sv,1)
    if (nDim == sam%ndof1) then
       rscr => getRealScratchArray(sam%neq+sam%ndof,ierr)
       if (ierr < 0) return
       sveq  => rscr(1:sam%neq)
       svdof => rscr(1+sam%neq:sam%neq+sam%ndof)
    else if (nDim == sam%neq) then
       svdof => getRealScratchArray(sam%ndof,ierr)
       if (ierr < 0) goto 190
    else
       ierr = internalError('writeDisplacementDB: Invalid vector dimension')
       return
    end if

    !! Write translations only for efficient loading of animation

    do iCmp = 1, nComp

       !! Expand displacement vector from equation order to DOF order
       if (nDim == sam%neq) then
          call csExpand (sam,sv(:,iCmp),svdof,1.0_dp,0.0_dp)
       else
          sveq = 0.0_dp
          do i = 1, sam%ndof1
             sveq(sam%meqn1(i)) = sv(i,iCmp)
          end do
          call csExpand (sam,sveq,svdof,1.0_dp,0.0_dp)
       end if

       if (.not. present(mlc)) then
          !! Normalize such that the largest nodal displacement equals 1.0
          maxN = 0.0_dp
          do i = 1, sam%nnod
             j = sam%madof(i)
             if (sam%madof(i+1) > j+2) then
                maxN = max(maxN,dot_product(svdof(j:j+2),svdof(j:j+2)))
             end if
          end do
          if (maxN > epsDiv0_p) svdof = svdof/sqrt(maxN)
       else if (present(lpu) .and. iCmp <= size(mlc)) then
          !! Debug output of displacement vector
          write(lpu,600) mlc(iCmp)
          do i = 1, sam%nnod
             if (sam%madof(i) < sam%madof(i+1)) then
                write(lpu,601) sam%minex(i),svdof(sam%madof(i):sam%madof(i+1)-1)
             end if
          end do
       end if

       if (present(mlc)) then
          !! Write load case identification
          if (iCmp <= size(mlc)) then
             call writeTimeStepDB (rdb,int(iCmp,i8),real(mlc(iCmp),dp),ierr)
          else
             call writeTimeStepDB (rdb,int(iCmp,i8),real(iCmp,dp),ierr)
          end if
          !! Write the whole displacement vector to database file
          call writeRDB (rdb,svdof,ierr,.true.)
          if (ierr < 0) goto 190
       end if

       !! Write nodal translations to database file
       do i = 1, sam%nnod
          j = sam%madof(i)
          if (sam%madof(i+1) < j+3) cycle ! Skip DOF-less node
          call writeRDB (rdb,svdof(j:j+2),ierr)
          if (ierr < 0) goto 190
       end do

    end do

    return

190 call reportError (debugFileOnly_p,'writeDisplacementDB')

600 format(/5X,'NODAL DISPLACEMENTS FOR LOAD CASE',I6,':' &
         &//5X,'Node #',4X,'X-disp',7X,'Y-disp',7X,'Z-disp', &
         &              7X,'X-rot ',7X,'Y-rot ',7X,'Z-rot')
601 format(I9,2X,1P,6E13.5)

  end subroutine writeDisplacementDB

end module saveReducerModule
