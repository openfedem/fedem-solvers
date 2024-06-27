!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file saveStressModule.f90
!> @brief Save of stress recovery variables to result database files.

!!==============================================================================
!> @brief Module with subroutines for saving of stress recovery results.
!>
!> @details This module contains a collection of subroutines for writing
!> stress recovery results data of an FE part to the results database.

module saveStressModule

  implicit none

  private

  !> @brief Variable- and item group indices for result database files.
  type HeaderId
     integer :: idDis(3),idRot(3)   ! Translation and rotation variables
     integer :: idStress(3)         ! Stress tensor variables
     integer :: idStrain(3)         ! Strain tensor variables
     integer :: idSR1,idSR2,idCurva ! Stress resultant and curvature variables
     integer :: idSF1,idSF2         ! Sectional force variables (beam elements)
     integer :: idSS1,idSS2         ! Sectional strain variables (beam elements)
     integer :: idStressVm1,idStrainVm1 !Von Mises stress and strain
     integer :: idStressMaxP1,idStrainMaxP1 !Maximum principal stess and strain
     integer :: idStressMinP1,idStrainMinP1 !Minimum principal stress and strain
     integer :: idStressMaxS1,idStrainMaxS1 !Max shear stress and strain
     integer :: idEnergyDens        ! Modes scaled strain energy density
     integer :: idBEAM
     integer :: idTRI3, idQUAD4
     integer :: idTRI6, idQUAD8
     integer :: idTET4, idTET10
     integer :: idPEN6, idPEN15
     integer :: idHEX8, idHEX20
  end type HeaderId

  logical, save :: lNodes = .true. !< Flag displacement output as nodal results
  integer, save :: lVector = 0     !< Flag displacement output as vector results

  !> @brief Processing order for the elements.
  !> @details Used only when direct VTF export also is done in the same run.
  integer, allocatable, save, public :: saveElmOrder(:)

  public :: writeStressHeader, writeModesHeader, writeModeHeader
  public :: initWriteDisp, writeDisplacementDB, writeStressDB, writeStrMeasureDB


contains

  !!============================================================================
  !> @brief Determines the output format for nodal deformations.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 22 Jun 2008

  subroutine initWriteDisp (sam,writeNodes)

    use SamModule             , only : SamType
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_isTrue

    type(SamType)    , intent(in)  :: sam
    logical, optional, intent(out) :: writeNodes

    !! Local variables
    integer :: i

    !! --- Logic section ---

    if (ffa_cmdlinearg_isTrue('write_vector')) then
       !! Check that all nodes have (at least 3 DOFs)
       lVector = 1
       do i = 1, sam%nnod
          if (sam%madof(i+1) < sam%madof(i)+3) lVector = 0
       end do

       if (lVector > 0) then
          lNodes = ffa_cmdlinearg_isTrue('write_nodes')
       end if
    end if

    if (present(writeNodes)) then
       writeNodes = lNodes
    end if

  end subroutine initWriteDisp


  !!============================================================================
  !> @brief Administers writing of results database header for fedem_stress.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 20 Nov 2000

  subroutine writeStressHeader (resFileName,modelFileName,rdb,sam,sup, &
       &                        iDef,ierr,lStr,bufRat)

    use kindModule            , only : lfnam_p, dp, i8
    use RDBModule             , only : RDBType, idatd, nullifyRDB, openRDBfile
    use RDBModule             , only : openHeaderFiles, writeTimeStepHeader
    use SamModule             , only : SamType
    use SupElTypeModule       , only : SupElType
    use IdTypeModule          , only : writeIdHeader
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getbool
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getstring

    character(len=*) , intent(in)  :: resFileName, modelFileName
    type(RDBType)    , intent(out) :: rdb
    type(SamType)    , intent(in)  :: sam
    type(SupElType)  , intent(in)  :: sup
    integer          , intent(out) :: ierr
    integer, optional, intent(in)  :: iDef
    logical, optional, intent(in)  :: lStr(:)
    real(dp),optional, intent(in)  :: bufRat

    ! Local variables
    type(HeaderId)     :: header
    integer            :: nBits
    logical            :: lDouble, lStrRes(12)
    logical, parameter :: lDef = .false. ! No element result point deformations
    logical, parameter :: lSS  = .false. ! No stress resultant conjugate strains
    character(lfnam_p) :: linkFileName

    ! --- Logic section ---

    ierr = 0
    call nullifyRDB (rdb)
    header%idDis     = 0
    header%idRot     = 0
    header%idStress  = 0
    header%idStrain  = 0
    header%idSR1     = 0
    header%idSR2     = 0
    header%idCurva   = 0
    header%idSF1     = 0
    header%idSF2     = 0
    header%idSS1     = 0
    header%idSS2     = 0
    header%idStressVm1 = 0
    header%idStrainVm1 = 0
    header%idStressMaxP1 = 0
    header%idStrainMaxP1 = 0
    header%idStressMinP1 = 0
    header%idStrainMinP1 = 0
    header%idStressMaxS1 = 0
    header%idStrainMaxS1 = 0
    header%idBEAM  = 0
    header%idTRI3  = 0
    header%idQUAD4 = 0
    header%idTRI3  = 0
    header%idQUAD4 = 0
    header%idTET4  = 0
    header%idTET10 = 0
    header%idPEN6  = 0
    header%idPEN15 = 0
    header%idHEX8  = 0
    header%idHEX20 = 0

    if (.not. present(lStr)) then
       call ffa_cmdlinearg_getbool ('SR',lStrRes(10))
       call ffa_cmdlinearg_getbool ('stress',lStrRes(11))
       call ffa_cmdlinearg_getbool ('strain',lStrRes(12))
       call ffa_cmdlinearg_getbool ('vmStress',lStrRes(1))
       call ffa_cmdlinearg_getbool ('vmStrain',lStrRes(5))
       call ffa_cmdlinearg_getbool ('maxPStress',lStrRes(2))
       call ffa_cmdlinearg_getbool ('maxPStrain',lStrRes(6))
       call ffa_cmdlinearg_getbool ('minPStress',lStrRes(3))
       call ffa_cmdlinearg_getbool ('minPStrain',lStrRes(7))
       call ffa_cmdlinearg_getbool ('maxSStress',lStrRes(4))
       call ffa_cmdlinearg_getbool ('maxSStrain',lStrRes(8))
       lStrRes(9) = .false. ! No scaled strain energy density yet
       if (present(iDef) .or. any(lStrRes)) then
          call ffa_cmdlinearg_getstring ('linkfile',linkFileName)
          call ffa_cmdlinearg_getbool ('double',lDouble)
       else
          return ! No result output requested
       end if
    else if (associated(sup%rcy)) then
       !! Solver mode, write only deformation and some scalar stress measures
       lStrRes = .false.
       lStrRes(1:size(lStr)) = lStr
       if (iDef > 0 .or. any(lStrRes)) then
          linkFileName = sup%rcy%fileName(4)
          lDouble = .false.
       else
          return ! No result output requested
       end if
    else
       return
    end if

    if (lDouble) then
       nBits = 64 ! Write all results in double precision
    else
       nBits = 32 ! Write all results in single precision
    end if

    ! Open the temporary header files
    call openHeaderFiles ('fedem_stress',ierr, &
         &                'response data base file',modelFileName,linkFileName)
    if (ierr /= 0) return

    ! Write stress results definitions

    call writeTimeStepHeader (rdb)
    call writeIdHeader ('Part',sup%id,idatd,.true.)
    call writeNodesHeader (header,rdb,sam,'Dynamic response',nBits,iDef,.false.)
    call writeElementsHeader (header,rdb,sam,nBits,lStrRes(10),lSS, &
         &                    lStrRes(11),lStrRes(12),lStrRes(1:9),lDef)
    write(idatd,"('}')")

    rdb%nBstep = int(rdb%nBytes)
    if (present(bufRat)) then
       rdb%bufSize = int(rdb%nBytes*bufRat)
       rdb%nBytes  = 0_i8
    end if

    ! Open the results database file and copy the file header to it

    call openRDBfile (rdb,ierr,resFileName,'#FEDEM response data')

  end subroutine writeStressHeader


  !!============================================================================
  !> @brief Administers writing of results database header for fedem_modes.
  !>
  !> @details This version writes one file for each mode shape.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 20 Nov 2000

  subroutine writeModeHeader (resFileName,modelFileName,iMode,lCompl, &
       &                      rdb,sam,sup,ierr)

    use kindModule            , only : lfnam_p
    use RDBModule             , only : RDBType, idatd, nullifyRDB, openRDBfile
    use RDBModule             , only : openHeaderFiles, writeTimeStepHeader
    use SamModule             , only : SamType
    use SupElTypeModule       , only : SupElType
    use IdTypeModule          , only : writeIdHeader
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getbool
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getstring

    character(len=*), intent(in)  :: resFileName, modelFileName
    integer         , intent(in)  :: iMode
    logical         , intent(in)  :: lCompl
    type(RDBType)   , intent(out) :: rdb
    type(SamType)   , intent(in)  :: sam
    type(SupElType) , intent(in)  :: sup
    integer         , intent(out) :: ierr

    ! Local variables
    type(HeaderId)     :: header
    integer            :: nBits
    logical            :: lStrRes(9), lDouble
    character(len=16)  :: igName
    character(lfnam_p) :: linkFileName

    ! --- Logic section ---

    call nullifyRDB (rdb)
    header%idDis = 0
    header%idRot = 0
    header%idEnergyDens = 0
    header%idBEAM  = 0
    header%idTRI3  = 0
    header%idQUAD4 = 0
    header%idTRI6  = 0
    header%idQUAD8 = 0
    header%idTET4  = 0
    header%idTET10 = 0
    header%idPEN6  = 0
    header%idPEN15 = 0
    header%idHEX8  = 0
    header%idHEX20 = 0
    if (iMode > 0) then
       write(igName,"('Mode',i3)") iMode
    else
       igName = 'Dynamic response'
    end if

    lStrRes = .false.
    call ffa_cmdlinearg_getbool ('energy_density',lStrRes(9))
    call ffa_cmdlinearg_getbool ('double',lDouble)
    if (lDouble) then
       nBits = 64 ! Write all results in double precision
    else
       nBits = 32 ! Write all results in single precision
    end if

    ! Open the temporary header files

    call ffa_cmdlinearg_getstring ('linkfile',linkFileName)
    call openHeaderFiles ('fedem_modes',ierr, &
         &                'modes data base file',modelFileName,linkFileName)
    if (ierr /= 0) return

    ! Write modes results definitions

    call writeTimeStepHeader (rdb)
    call writeIdHeader ('Part',sup%id,idatd,.true.)
    call writeNodesHeader (header,rdb,sam,trim(igName),nBits,1,lCompl)
    if (iMode > 0) then
       call writeElementsHeader (header,rdb,sam,nBits,.false.,.false., &
            &                    .false.,.false.,lStrRes,.false.)
    end if
    write(idatd,"('}')")

    ! Open the results database file and copy the file header to it

    call openRDBfile (rdb,ierr,resFileName,'#FEDEM modal data')

  end subroutine writeModeHeader


  !!============================================================================
  !> @brief Administers writing of results database header for fedem_modes.
  !>
  !> @details This version writes all mode shapes as a vector to the same file.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 20 Nov 2000

  subroutine writeModesHeader (resFileName,modelFileName,iModes,lCompl, &
       &                       rdb,sam,sup,ierr)

    use kindModule            , only : lfnam_p, i8
    use RDBModule             , only : RDBType, idatd, nullifyRDB, openRDBfile
    use RDBModule             , only : openHeaderFiles, writeTimeStepHeader
    use SamModule             , only : SamType
    use SupElTypeModule       , only : SupElType
    use IdTypeModule          , only : writeIdHeader
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getbool
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getstring

    character(len=*), intent(in)  :: resFileName, modelFileName
    integer         , intent(in)  :: iModes(:)
    logical         , intent(in)  :: lCompl
    type(RDBType)   , intent(out) :: rdb
    type(SamType)   , intent(in)  :: sam
    type(SupElType) , intent(in)  :: sup
    integer         , intent(out) :: ierr

    ! Local variables
    type(HeaderId)     :: header
    integer            :: i, nBits
    logical            :: lDouble
    character(len=16)  :: igName
    character(lfnam_p) :: linkFileName

    ! --- Logic section ---

    call nullifyRDB (rdb)
    header%idDis(3) = 0
    header%idRot(3) = 0

    lVector = 2 ! Write all mode shape vectors to the same file

    call ffa_cmdlinearg_getbool ('double',lDouble)
    if (lDouble) then
       nBits = 64 ! Write all results in double precision
    else
       nBits = 32 ! Write all results in single precision
    end if

    ! Open the temporary header files

    call ffa_cmdlinearg_getstring ('linkfile',linkFileName)
    call openHeaderFiles ('fedem_modes',ierr, &
         &                'modes data base file',modelFileName,linkFileName)
    if (ierr /= 0) return

    ! Write modes results definitions

    call writeTimeStepHeader (rdb)
    call writeIdHeader ('Part',sup%id,idatd,.true.)
    write(idatd,610)
    call writeNodesHeader (header,rdb,sam,'Dynamic response',nBits,1,.false.)
    do i = 1, size(iModes)
       write(igName,"('Mode',i3)") iModes(i)
       call writeNodesHeader (header,rdb,sam,trim(igName),nBits,1,lCompl)
    end do
    write(idatd,"('  ]'/'}')")

    rdb%nBstep = int(rdb%nBytes)
    rdb%nBytes = 0_i8

    ! Open the results database file and copy the file header to it

    call openRDBfile (rdb,ierr,resFileName,'#FEDEM modal data')

610 format('  [;"Vectors";')

  end subroutine writeModesHeader


  !!==========================================================================
  !> @brief Writes nodal result definitions to the temporary header files.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 20 Nov 2000

  subroutine writeNodesHeader (header,rdb,sam,igName,nbit,iDef,lCompl)

    use RDBModule, only : RDBType, ivard, iitem, idatd
    use SamModule, only : SamType

    type(HeaderId)  , intent(inout) :: header
    type(RDBType)   , intent(inout) :: rdb
    type(SamType)   , intent(in)    :: sam
    character(len=*), intent(in)    :: igName
    integer         , intent(in)    :: nbit
    integer,optional, intent(in)    :: iDef
    logical         , intent(in)    :: lCompl

    !! Local variables
    integer :: i, id3DOF, id6DOF, nRot

    !! --- Logic section ---

    if (.not. present(iDef)) return

    if (iDef > 0 .and. lNodes) then

       !! Write out deformations as nodal data, including rotations

       id3DOF = 0
       id6DOF = 0
       write(idatd,600)
       do i = 1, sam%nnod
          if (sam%mpar(33) > 0 .and. sam%mnnn(i) > 0) cycle ! outside group
          if (sam%madof(i+1) < sam%madof(i)+3) cycle ! Fixed node, no results

          if (header%idDis(1) == 0) then
             rdb%nVar = rdb%nVar + 1
             header%idDis(1) = rdb%nVar
             write(ivard,601) rdb%nVar,'Translational deformation', &
                  &           'LENGTH',nbit,'VEC3','"d_x","d_y","d_z"'
          end if
          if (header%idDis(2) == 0 .and. iDef > 1 .and. .not.lCompl) then
             rdb%nVar = rdb%nVar + 1
             header%idDis(2) = rdb%nVar
             write(ivard,601) rdb%nVar,'Total translation', &
                  &           'LENGTH',nbit,'VEC3','"u_x","u_y","u_z"'
          end if
          if (sam%madof(i+1) > sam%madof(i)+5) then
             if (header%idRot(1) == 0) then
                rdb%nVar = rdb%nVar + 1
                header%idRot(1) = rdb%nVar
                write(ivard,601) rdb%nVar,'Angular deformation','ANGLE', &
                     &           nbit,'ROT3','"theta_x","theta_y","theta_z"'
             end if
             if (header%idRot(2) == 0 .and. iDef > 1 .and. .not.lCompl) then
                rdb%nVar = rdb%nVar + 1
                header%idRot(2) = rdb%nVar
                write(ivard,601) rdb%nVar,'Total rotation','ANGLE', &
                     &           nbit,'ROT3','"theta_x","theta_y","theta_z"'
             end if
             if (id6DOF == 0) then
                rdb%nIG = rdb%nIG + 1
                id6DOF  = rdb%nIG
                if (lCompl) then
                   write(iitem,602) id6DOF,igName, &
                        &           header%idDis(1),header%idRot(1), &
                        &           header%idDis(1),header%idRot(1)
                else if (iDef > 1) then
                   write(iitem,614) id6DOF,igName, &
                        &           header%idDis(1),header%idRot(1), &
                        &           header%idDis(2),header%idRot(2)
                else
                   write(iitem,603) id6DOF,igName, &
                        &           header%idDis(1),header%idRot(1)
                end if
             end if
             write(idatd,606) sam%minex(i),id6DOF
             if (lCompl .or. iDef > 1) then
                rdb%nBytes = rdb%nBytes + 12*nbit/8
             else
                rdb%nBytes = rdb%nBytes + 6*nbit/8
             end if
          else
             if (id3DOF == 0) then
                rdb%nIG = rdb%nIG + 1
                id3DOF  = rdb%nIG
                if (lCompl) then
                   write(iitem,604) id3DOF,igName, &
                        &           header%idDis(1),header%idDis(1)
                else if (iDef > 1) then
                   write(iitem,603) id3DOF,igName,header%idDis(1:2)
                else
                   write(iitem,605) id3DOF,igName,header%idDis(1)
                end if
             end if
             write(idatd,606) sam%minex(i),id3DOF
             if (lCompl .or. iDef > 1) then
                rdb%nBytes = rdb%nBytes + 6*nbit/8
             else
                rdb%nBytes = rdb%nBytes + 3*nbit/8
             end if
          end if

       end do
       write(idatd,"('  ]')")

    end if
    if (iDef == 1 .and. lVector > 0) then

       !! All nodes have (at least) three translational DOFs.
       !! Write out deformations as a big 3*nNode vector.
       !! The rotations (if any) are written as a separate vector.

       nRot = 0
       do i = 1, sam%nnod
          if (sam%minex(i) < 0) cycle ! Skip internal beam nodes
          if (sam%madof(i+1) < sam%madof(i)+6) cycle ! No rotation for this node
          nRot = nRot + 3
          if (lCompl) then
             rdb%nBytes = rdb%nBytes + 3*nbit/4
          else
             rdb%nBytes = rdb%nBytes + 3*nbit/8
          end if
       end do

       if (lCompl) then
          rdb%nBytes = rdb%nBytes + 3*(sam%nnod-sam%mpar(23))*nbit/4
       else
          rdb%nBytes = rdb%nBytes + 3*(sam%nnod-sam%mpar(23))*nbit/8
       end if

       if (header%idDis(3) == 0) then
          rdb%nVar = rdb%nVar + 1
          header%idDis(3) = rdb%nVar
          write(ivard,611) rdb%nVar,'Translational deformation', &
               &                    'LENGTH',nbit,3*(sam%nnod-sam%mpar(23))
       end if
       if (nRot > 0 .and. header%idRot(3) == 0) then
          rdb%nVar = rdb%nVar + 1
          header%idRot(3) = rdb%nVar
          write(ivard,611) rdb%nVar,'Angular deformation','LENGTH',nbit,nRot
       end if

       if (lVector == 1) write(idatd,610)

       if (nRot > 0) then
          if (lCompl) then
             write(idatd,618) igName,header%idDis(3),header%idRot(3), &
                  &                  header%idDis(3),header%idRot(3)
          else
             write(idatd,619) igName,header%idDis(3),header%idRot(3)
          end if
       else
          if (lCompl) then
             write(idatd,616) igName,header%idDis(3),header%idDis(3)
          else
             write(idatd,617) igName,header%idDis(3)
          end if
       end if

       if (lVector == 1) write(idatd,"('  ]')")

    end if

600 format('  [;"Nodes";')
601 format('<',i3,';"',a,'";',a,';FLOAT;',i2,';',a,';(3);((',a,'))>')
602 format('[',i3,';"',a,'";[;"Re";<',i3,'><',i3,'>][;"Im";<',i3,'><',i3,'>]]')
603 format('[',i3,';"',a,'";<',i3,'><',i3,'>]')
604 format('[',i3,';"',a,'";[;"Re";<',i3,'>][;"Im";<',i3,'>]]')
605 format('[',i3,';"',a,'";<',i3,'>]')
606 format('    [;',i8,';[',i3,']]')
610 format('  [;"Vectors";')
611 format('<',i3,';"',a,'";',a,';FLOAT;',i2,';VECTOR;(',i8,')>')
614 format('[',i3,';"',a,'";<',i3,'><',i3,'><',i3,'><',i3,'>]')
616 format('    [;"',a,'";[;"Re";<',i3,'>][;"Im";<',i3,'>]]')
617 format('    [;"',a,'";<',i3,'>]')
618 format('    [;"',a,'";[;"Re";<',i3,'><',i3,'>][;"Im";<',i3,'><',i3,'>]]')
619 format('    [;"',a,'";<',i3,'><',i3,'>]')

  end subroutine writeNodesHeader


  !!============================================================================
  !> @brief Writes element result definitions to the temporary header files.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 20 Nov 2000

  subroutine writeElementsHeader (header,rdb,sam,nbit, &
         &                        lSR,lSS,lStress,lStrain,lStrRes,lDef)

    use RDBModule              , only : RDBType, idatd
    use SamModule              , only : SamType
    use FFlLinkHandlerInterface, only : ffl_getElmId

    type(headerId), intent(inout) :: header
    type(RDBType) , intent(inout) :: rdb
    type(SamType) , intent(in)    :: sam
    integer       , intent(in)    :: nbit
    logical       , intent(in)    :: lSR, lSS, lStress, lStrain
    logical       , intent(in)    :: lStrRes(:), lDef

    ! Local variables
    integer :: elmno, iel, jel, itemGroup, nPoint = 0

    ! --- Logic section ---

    if (.not. (lSR .or. lSS .or. lStress .or. lStrain .or. any(lStrRes))) return

    write(idatd,600)

    do jel = 1, sam%nel
       if (allocated(saveElmOrder)) then
          iel = saveElmOrder(jel)
       else
          iel = jel
       end if

       elmno = ffl_getElmId(iel)
       if (elmno <= 0) cycle

       select case(sam%melcon(iel))

       case(11) !! Two-noded beam
          if (.not. (lSR .or. lSS .or. nPoint > 0)) cycle

          if (header%idBEAM == 0) then
             call writeBeamHeader (header,rdb,header%idBEAM,'BEAM2',2,nPoint, &
                  &                nbit,lSR,lSS,lStress,lStrain,lDef,lStrRes)
          end if
          itemGroup = header%idBEAM

       case(21,23) !! Three-noded thin shell
          if (header%idTRI3 == 0) then
             call writeShellHeader (header,rdb,header%idTRI3,'TRI3',3, &
                  &                 nbit,lSR,lSS,lStress,lStrain,lDef,lStrRes)
          end if
          itemGroup = header%idTRI3

       case(22,24) !! Four-noded thin shell
          if (header%idQUAD4 == 0) then
             call writeShellHeader (header,rdb,header%idQUAD4,'QUAD4',4, &
                  &                 nbit,lSR,lSS,lStress,lStrain,lDef,lStrRes)
          end if
          itemGroup = header%idQUAD4

       case(31) !! Six-noded thick shell
          if (header%idTRI6 == 0) then
             call writeShellHeader (header,rdb,header%idTRI6,'TRI6',6, &
                  &                 nbit,lSR,lSS,lStress,lStrain,lDef,lStrRes)
          end if
          itemGroup = header%idTRI6

       case(32) !! Eight-noded thin shell
          if (header%idQUAD8 == 0) then
             call writeShellHeader (header,rdb,header%idQUAD8,'QUAD8',8, &
                  &                 nbit,lSR,lSS,lStress,lStrain,lDef,lStrRes)
          end if
          itemGroup = header%idQUAD8

       case(41) !! Ten-noded tetrahedron
          if (header%idTET10 == 0) then
             call writeSolidHeader (header,rdb,header%idTET10,'TET10',10, &
                  &                 nbit,lStress,lStrain,lDef,lStrRes)
          end if
          itemGroup = header%idTET10

       case(42) !! Fiften-noded wedge
          if (header%idPEN15 == 0) then
             call writeSolidHeader (header,rdb,header%idPEN15,'WEDG15',15, &
                  &                 nbit,lStress,lStrain,lDef,lStrRes)
          end if
          itemGroup = header%idPEN15

       case(43) !! Twenty-noded hexahedron
          if (header%idHEX20 == 0) then
             call writeSolidHeader (header,rdb,header%idHEX20,'HEX20',20, &
                  &                 nbit,lStress,lStrain,lDef,lStrRes)
          end if
          itemGroup = header%idHEX20

       case(44) !! Eight-noded hexahedron
          if (header%idHEX8 == 0) then
             call writeSolidHeader (header,rdb,header%idHEX8,'HEX8',8, &
                  &                 nbit,lStress,lStrain,lDef,lStrRes)
          end if
          itemGroup = header%idHEX8

       case(45) !! Four-noded tetrahedron
          if (header%idTET4 == 0) then
             call writeSolidHeader (header,rdb,header%idTET4,'TET4',4, &
                  &                 nbit,lStress,lStrain,lDef,lStrRes)
          end if
          itemGroup = header%idTET4

       case(46) !! Six-noded wedge
          if (header%idPEN6 == 0) then
             call writeSolidHeader (header,rdb,header%idPEN6,'WEDG6',6, &
                  &                 nbit,lStress,lStrain,lDef,lStrRes)
          end if
          itemGroup = header%idPEN6

       case default
          cycle

       end select

       write(idatd,601) elmno,itemGroup

    end do
    write(idatd,"('  ]')")

600 format('  [;"Elements";')
601 format('    [;',i8,';[',i3,']]')

  end subroutine writeElementsHeader


  !!============================================================================
  !> @brief Writes beam element item group definition to temporary header files.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 20 Nov 2000

  subroutine writeBeamHeader (header,rdb,itemGroup,eltype,nelnod,npoint,nbit, &
       &                      lSR,lSS,lStress,lStrain,lDef,lStrRes)

    use RDBModule, only : RDBType, ivard, iitem

    type(headerId)  , intent(inout) :: header
    type(RDBType)   , intent(inout) :: rdb
    integer         , intent(inout) :: itemGroup
    character(len=*), intent(in)    :: eltype
    integer         , intent(in)    :: nelnod, npoint, nbit
    logical         , intent(in)    :: lSR, lSS, lStress, lStrain, lDef
    logical         , intent(in)    :: lStrRes(:)

    ! Local variables
    integer           :: i
    character(len=64) :: chVars

    ! --- Logic section ---

    rdb%nIG   = rdb%nIG + 1
    itemGroup = rdb%nIG
    write(iitem,600) itemGroup,trim(eltype)

    if (lSR .or. lSS) then

       i = 1
       if (lSR) then
          if (header%idSF1 == 0) then
             rdb%nVar = rdb%nVar + 1
             header%idSF1 = rdb%nVar
             write(ivard,601) rdb%nVar,'Beam sectional force', &
                  &           'FORCE',nbit,'VEC3',3,'"N","V_y","V_z"'
          end if
          write(chVars(i:),"('<',i3,'>')") header%idSF1
          i = i + 5
          if (header%idSF2 == 0) then
             rdb%nVar = rdb%nVar + 1
             header%idSF2 = rdb%nVar
             write(ivard,601) rdb%nVar,'Beam sectional moment', &
                  &           'FORCE*LENGTH',nbit,'VEC3',3,'"M_x","M_y","M_z"'
          end if
          write(chVars(i:),"('<',i3,'>')") header%idSF2
          i = i + 5
       end if

       if (lSS) then
          if (header%idSS1 == 0) then
             rdb%nVar = rdb%nVar + 1
             header%idSS1 = rdb%nVar
             write(ivard,601) rdb%nVar,'Strain','NONE',nbit,'VEC3',3, &
                  &           '"epsilon_xx","epsilon_xy","epsilon_xz"'
          end if
          write(chVars(i:),"('<',i3,'>')") header%idSS1
          i = i + 5
          if (header%idSS2 == 0) then
             rdb%nVar = rdb%nVar + 1
             header%idSS2 = rdb%nVar
             write(ivard,601) rdb%nVar,'Curvature','NONE/LENGTH', &
                  &           nbit,'VEC3',3,'"kappa_xx","kappa_xy","kappa_xz"'
          end if
          write(chVars(i:),"('<',i3,'>')") header%idSS2
          i = i + 5
       end if

       write(iitem,603) 'Element nodes','Basic'
       write(iitem,602) (i,trim(chVars),i=1,nelnod)
       write(iitem,"('    ]'/'  ]')")

    end if
    if ((lStress .or. lStrain .or. any(lStrRes)) .and. npoint > 0) then

       i = 1
       if (lStress) then
          call writeStressTensorDef (header%idStress(1),rdb%nVar,1,nbit)
          write(chVars(i:),"('<',i3,'>')") header%idStress(1)
          i = i + 5
       end if

       if (lStrain) then
          call writeStrainTensorDef (header%idStrain(1),rdb%nVar,1,nbit)
          write(chVars(i:),"('<',i3,'>')") header%idStrain(1)
          i = i + 5
       end if

       if (lDef) then
          if (header%idDis(1) == 0) then
             rdb%nVar = rdb%nVar + 1
             header%idDis(1) = rdb%nVar
             write(ivard,601) rdb%nVar,'Translational deformation', &
                  &           'LENGTH',nbit,'VEC3',3,'"d_x","d_y","d_z"'
          end if
          write(chVars(i:),"('<',i3,'>')") header%idDis(1)
          i = i + 5
       end if

       ! ---- Stress measure section
       if (lStrRes(1)) then
          if (header%idStressVm1 == 0) then
             rdb%nVar = rdb%nVar + 1
             header%idStressVm1 = rdb%nVar
             write(ivard,605) rdb%nVar,'Von Mises stress','FORCE/AREA',nbit
          end if
          write(chVars(i:),"('<',i3,'>')") header%idStressVm1
          i = i + 5
       end if

       if (lStrRes(2)) then
          if (header%idStressMaxP1 == 0) then
             rdb%nVar = rdb%nVar + 1
             header%idStressMaxP1 = rdb%nVar
             write(ivard,605) rdb%nVar,'Max principal stress','FORCE/AREA',nbit
          end if
          write(chVars(i:),"('<',i3,'>')") header%idStressMaxP1
          i = i + 5
       end if

       if (lStrRes(3)) then
          if (header%idStressMinP1 == 0) then
             rdb%nVar = rdb%nVar + 1
             header%idStressMinP1 = rdb%nVar
             write(ivard,605) rdb%nVar,'Min principal stress','FORCE/AREA',nbit
          end if
          write(chVars(i:),"('<',i3,'>')") header%idStressMinP1
          i = i + 5
       end if

       if (lStrRes(4)) then
          if (header%idStressMaxS1 == 0) then
             rdb%nVar = rdb%nVar + 1
             header%idStressMaxS1 = rdb%nVar
             write(ivard,605) rdb%nVar,'Max shear stress','FORCE/AREA',nbit
          end if
          write(chVars(i:),"('<',i3,'>')") header%idStressMaxS1
          i = i + 5
       end if

       ! ---- Strain measure section
       if (lStrRes(5)) then
          if (header%idStrainVm1 == 0) then
             rdb%nVar = rdb%nVar + 1
             header%idStrainVm1 = rdb%nVar
             write(ivard,605) rdb%nVar,'Von Mises strain','NONE',nbit
          end if
          write(chVars(i:),"('<',i3,'>')") header%idStrainVm1
          i = i + 5
       end if

       if (lStrRes(6)) then
          if (header%idStrainMaxP1 == 0) then
             rdb%nVar = rdb%nVar + 1
             header%idStrainMaxP1 = rdb%nVar
             write(ivard,605) rdb%nVar,'Max principal strain','NONE',nbit
          end if
          write(chVars(i:),"('<',i3,'>')") header%idStrainMaxP1
          i = i + 5
       end if

       if (lStrRes(7)) then
          if (header%idStrainMinP1 == 0) then
             rdb%nVar = rdb%nVar + 1
             header%idStrainMinP1 = rdb%nVar
             write(ivard,605) rdb%nVar,'Min principal strain','NONE',nbit
          end if
          write(chVars(i:),"('<',i3,'>')") header%idStrainMinP1
          i = i + 5
       end if

       if (lStrRes(8)) then
          if (header%idStrainMaxS1 == 0) then
             rdb%nVar = rdb%nVar + 1
             header%idStrainMaxS1 = rdb%nVar
             write(ivard,605) rdb%nVar,'Max shear strain','NONE',nbit
          end if
          write(chVars(i:),"('<',i3,'>')") header%idStrainMaxS1
          i = i + 5
       end if

       if (lStrRes(9)) then
          if (header%idEnergyDens == 0) then
             rdb%nVar = rdb%nVar + 1
             header%idEnergyDens = rdb%nVar
             write(ivard,605) rdb%nVar,'Scaled strain energy density', &
                  &           'FORCE/AREA',nbit
          end if
          write(chVars(i:),"('<',i3,'>')") header%idEnergyDens
          i = i + 5
       end if

       ! ---- Write item group
       write(iitem,603) 'Evaluation points','Basic'
       write(iitem,602) (i,trim(chVars),i=1,npoint)
       write(iitem,"('    ]'/'  ]')")

    end if
    write(iitem,"(']')")

600 format('[',i3,';"',a,'";')
601 format('<',i3,';"',a,'";',a,';FLOAT;',i2,';',a,';(',i1,');((',a,'))>')
602 format('      [;',i2,';',a,']')
603 format('  [;"',a,'";'/'    [;"',a,'";')
605 format('<',i3,';"',a,'";',a,';FLOAT;',i2,';SCALAR>')

  end subroutine writeBeamHeader


  !!============================================================================
  !> @brief Writes shell item group definition to temporary header files.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 20 Nov 2000

  subroutine writeShellHeader (header,rdb,itemGroup,eltype,nelnod,nbit, &
       &                       lSR,lSS,lStress,lStrain,lDef,lStrRes)

    use RDBModule, only : RDBType, ivard, iitem

    type(headerId)  , intent(inout) :: header
    type(RDBType)   , intent(inout) :: rdb
    integer         , intent(inout) :: itemGroup
    character(len=*), intent(in)    :: eltype
    integer         , intent(in)    :: nelnod, nbit
    logical         , intent(in)    :: lSR, lSS, lStress, lStrain, lDef
    logical         , intent(in)    :: lStrRes(:)

    ! Local variables
    integer           :: i, j
    character(len=64) :: chVars

    ! --- Logic section ---

    rdb%nIG   = rdb%nIG + 1
    itemGroup = rdb%nIG
    write(iitem,600) itemGroup,trim(eltype)

    if (lSR .or. lSS) then

       i = 1
       if (lSR) then
          if (header%idSR1 == 0) then
             rdb%nVar = rdb%nVar + 1
             header%idSR1 = rdb%nVar
             write(ivard,601) rdb%nVar,'Shell stress resultant force', &
                  &           'FORCE/LENGTH',nbit,'TENSOR2',3, &
                  &           '"n_xx","n_yy","n_xy"'
          end if
          write(chVars(i:),"('<',i3,'>')") header%idSR1
          i = i + 5
          if (header%idSR2 == 0) then
             rdb%nVar = rdb%nVar + 1
             header%idSR2 = rdb%nVar
             write(ivard,601) rdb%nVar,'Shell stress resultant moment', &
                  &           'FORCE*LENGTH/LENGTH',nbit,'TENSOR2',3, &
                  &           '"m_xx","m_yy","m_xy"'
          end if
          write(chVars(i:),"('<',i3,'>')") header%idSR2
          i = i + 5
       end if

       if (lSS) then
          call writeStrainTensorDef (header%idStrain(2),rdb%nVar,2,nbit)
          write(chVars(i:),"('<',i3,'>')") header%idStrain(2)
          i = i + 5
          if (header%idCurva == 0) then
             rdb%nVar = rdb%nVar + 1
             header%idCurva = rdb%nVar
             write(ivard,601) rdb%nVar,'Curvature', &
                  &           'NONE/LENGTH',nbit,'TENSOR2',3, &
                  &           '"kappa_xx","kappa_yy","kappa_xy"'
          end if
          write(chVars(i:),"('<',i3,'>')") header%idCurva
          i = i + 5
       end if

       write(iitem,603) 'Basic'
       write(iitem,602) (i,trim(chVars),i=1,nelnod)
       write(iitem,"('    ]')")

    end if
    if (lStress .or. lStrain .or. any(lStrRes)) then

       if (nelnod > 4) then
          j = 3
       else
          j = 2
       end if
       i = 1
       if (lStress) then
          call writeStressTensorDef (header%idStress(j),rdb%nVar,j,nbit)
          write(chVars(i:),"('<',i3,'>')") header%idStress(j)
          i = i + 5
       end if

       if (lStrain) then
          call writeStrainTensorDef (header%idStrain(j),rdb%nVar,j,nbit)
          write(chVars(i:),"('<',i3,'>')") header%idStrain(j)
          i = i + 5
       end if

       if (lDef) then
          if (header%idDis(1) == 0) then
             rdb%nVar = rdb%nVar + 1
             header%idDis(1) = rdb%nVar
             write(ivard,601) rdb%nVar,'Translational deformation', &
                  &           'LENGTH',nbit,'VEC3',3,'"d_x","d_y","d_z"'
          end if
          write(chVars(i:),"('<',i3,'>')") header%idDis(1)
          i = i + 5
       end if

       ! ---- Stress measure section
       if (lStrRes(1)) then
          if (header%idStressVm1 == 0) then
             rdb%nVar = rdb%nVar + 1
             header%idStressVm1 = rdb%nVar
             write(ivard,605) rdb%nVar,'Von Mises stress','FORCE/AREA',nbit
          end if
          write(chVars(i:),"('<',i3,'>')") header%idStressVm1
          i = i + 5
       end if

       if (lStrRes(2)) then
          if (header%idStressMaxP1 == 0) then
             rdb%nVar = rdb%nVar + 1
             header%idStressMaxP1 = rdb%nVar
             write(ivard,605) rdb%nVar,'Max principal stress','FORCE/AREA',nbit
          end if
          write(chVars(i:),"('<',i3,'>')") header%idStressMaxP1
          i = i + 5
       end if

       if (lStrRes(3)) then
          if (header%idStressMinP1 == 0) then
             rdb%nVar = rdb%nVar + 1
             header%idStressMinP1 = rdb%nVar
             write(ivard,605) rdb%nVar,'Min principal stress','FORCE/AREA',nbit
          end if
          write(chVars(i:),"('<',i3,'>')") header%idStressMinP1
          i = i + 5
       end if

       if (lStrRes(4)) then
          if (header%idStressMaxS1 == 0) then
             rdb%nVar = rdb%nVar + 1
             header%idStressMaxS1 = rdb%nVar
             write(ivard,605) rdb%nVar,'Max shear stress','FORCE/AREA',nbit
          end if
          write(chVars(i:),"('<',i3,'>')") header%idStressMaxS1
          i = i + 5
       end if

       ! ---- Strain measure section
       if (lStrRes(5)) then
          if (header%idStrainVm1 == 0) then
             rdb%nVar = rdb%nVar + 1
             header%idStrainVm1 = rdb%nVar
             write(ivard,605) rdb%nVar,'Von Mises strain','NONE',nbit
          end if
          write(chVars(i:),"('<',i3,'>')") header%idStrainVm1
          i = i + 5
       end if

       if (lStrRes(6)) then
          if (header%idStrainMaxP1 == 0) then
             rdb%nVar = rdb%nVar + 1
             header%idStrainMaxP1 = rdb%nVar
             write(ivard,605) rdb%nVar,'Max principal strain','NONE',nbit
          end if
          write(chVars(i:),"('<',i3,'>')") header%idStrainMaxP1
          i = i + 5
       end if

       if (lStrRes(7)) then
          if (header%idStrainMinP1 == 0) then
             rdb%nVar = rdb%nVar + 1
             header%idStrainMinP1 = rdb%nVar
             write(ivard,605) rdb%nVar,'Min principal strain','NONE',nbit
          end if
          write(chVars(i:),"('<',i3,'>')") header%idStrainMinP1
          i = i + 5
       end if

       if (lStrRes(8)) then
          if (header%idStrainMaxS1 == 0) then
             rdb%nVar = rdb%nVar + 1
             header%idStrainMaxS1 = rdb%nVar
             write(ivard,605) rdb%nVar,'Max shear strain','NONE',nbit
          end if
          write(chVars(i:),"('<',i3,'>')") header%idStrainMaxS1
          i = i + 5
       end if

       if (lStrRes(9)) then
          if (header%idEnergyDens == 0) then
             rdb%nVar = rdb%nVar + 1
             header%idEnergyDens = rdb%nVar
             write(ivard,605) rdb%nVar,'Scaled strain energy density', &
                  &           'FORCE/AREA',nbit
          end if
          write(chVars(i:),"('<',i3,'>')") header%idEnergyDens
          i = i + 5
       end if

       ! ---- Write item group
       write(iitem,603) 'Top'
       write(iitem,602) (i,trim(chVars),i=1,nelnod)
       write(iitem,"('    ]')")
       write(iitem,603) 'Bottom'
       write(iitem,602) (i,trim(chVars),i=1,nelnod)
       write(iitem,"('    ]')")

    end if
    write(iitem,"('  ]'/']')")

600 format('[',i3,';"',a,'";'/'  [;"Element nodes";')
601 format('<',i3,';"',a,'";',a,';FLOAT;',i2,';',a,';(',i1,');((',a,'))>')
602 format('      [;',i2,';',a,']')
603 format('    [;"',a,'";')
605 format('<',i3,';"',a,'";',a,';FLOAT;',i2,';SCALAR>')

  end subroutine writeShellHeader


  !!============================================================================
  !> @brief Writes solid item group definition to temporary header files.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 20 Nov 2000

  subroutine writeSolidHeader (header,rdb,itemGroup,eltype,nelnod,nbit, &
       &                       lStress,lStrain,lDef,lStrRes)

    use RDBModule, only : RDBType, ivard, iitem

    type(headerId)  , intent(inout) :: header
    type(RDBType)   , intent(inout) :: rdb
    integer         , intent(inout) :: itemGroup
    character(len=*), intent(in)    :: eltype
    integer         , intent(in)    :: nelnod, nbit
    logical         , intent(in)    :: lStress, lStrain, lDef
    logical         , intent(in)    :: lStrRes(:)

    ! Local variables
    integer           :: i
    character(len=64) :: chVars

    ! --- Logic section ---

    rdb%nIG   = rdb%nIG + 1
    itemGroup = rdb%nIG
    write(iitem,600) itemGroup,trim(eltype)

    i = 1
    if (lStress) then
       call writeStressTensorDef (header%idStress(3),rdb%nVar,3,nbit)
       write(chVars(i:),"('<',i3,'>')") header%idStress(3)
       i = i + 5
    end if

    if (lStrain) then
       call writeStrainTensorDef (header%idStrain(3),rdb%nVar,3,nbit)
       write(chVars(i:),"('<',i3,'>')") header%idStrain(3)
       i = i + 5
    end if

    if (lDef) then
       if (header%idDis(1) == 0) then
          rdb%nVar = rdb%nVar + 1
          header%idDis(1) = rdb%nVar
          write(ivard,601) rdb%nVar,'Translational deformation', &
               &           'LENGTH',nbit,'VEC3',3,'"d_x","d_y","d_z"'
       end if
       write(chVars(i:),"('<',i3,'>')") header%idDis(1)
       i = i + 5
    end if

    ! ---- Stress measure section
    if (lStrRes(1)) then
       if (header%idStressVm1 == 0) then
          rdb%nVar = rdb%nVar + 1
          header%idStressVm1 = rdb%nVar
          write(ivard,605) rdb%nVar,'Von Mises stress','FORCE/AREA',nbit
       end if
       write(chVars(i:),"('<',i3,'>')") header%idStressVm1
       i = i + 5
    end if

    if (lStrRes(2)) then
       if (header%idStressMaxP1 == 0) then
          rdb%nVar = rdb%nVar + 1
          header%idStressMaxP1 = rdb%nVar
          write(ivard,605) rdb%nVar,'Max principal stress','FORCE/AREA',nbit
       end if
       write(chVars(i:),"('<',i3,'>')") header%idStressMaxP1
       i = i + 5
    end if

    if (lStrRes(3)) then
       if (header%idStressMinP1 == 0) then
          rdb%nVar = rdb%nVar + 1
          header%idStressMinP1 = rdb%nVar
          write(ivard,605) rdb%nVar,'Min principal stress','FORCE/AREA',nbit
       end if
       write(chVars(i:),"('<',i3,'>')") header%idStressMinP1
       i = i + 5
    end if

    if (lStrRes(4)) then
       if (header%idStressMaxS1 == 0) then
          rdb%nVar = rdb%nVar + 1
          header%idStressMaxS1 = rdb%nVar
          write(ivard,605) rdb%nVar,'Max shear stress','FORCE/AREA',nbit
       end if
       write(chVars(i:),"('<',i3,'>')") header%idStressMaxS1
       i = i + 5
    end if

    ! ---- Strain measure section
    if (lStrRes(5)) then
       if (header%idStrainVm1 == 0) then
          rdb%nVar = rdb%nVar + 1
          header%idStrainVm1 = rdb%nVar
          write(ivard,605) rdb%nVar,'Von Mises strain','NONE',nbit
       end if
       write(chVars(i:),"('<',i3,'>')") header%idStrainVm1
       i = i + 5
    end if

    if (lStrRes(6)) then
       if (header%idStrainMaxP1 == 0) then
          rdb%nVar = rdb%nVar + 1
          header%idStrainMaxP1 = rdb%nVar
          write(ivard,605) rdb%nVar,'Max principal strain','NONE',nbit
       end if
       write(chVars(i:),"('<',i3,'>')") header%idStrainMaxP1
       i = i + 5
    end if

    if (lStrRes(7)) then
       if (header%idStrainMinP1 == 0) then
          rdb%nVar = rdb%nVar + 1
          header%idStrainMinP1 = rdb%nVar
          write(ivard,605) rdb%nVar,'Min principal strain','NONE',nbit
       end if
       write(chVars(i:),"('<',i3,'>')") header%idStrainMinP1
       i = i + 5
    end if

    if (lStrRes(8)) then
       if (header%idStrainMaxS1 == 0) then
          rdb%nVar = rdb%nVar + 1
          header%idStrainMaxS1 = rdb%nVar
          write(ivard,605) rdb%nVar,'Max shear strain','NONE',nbit
       end if
       write(chVars(i:),"('<',i3,'>')") header%idStrainMaxS1
       i = i + 5
    end if

    if (lStrRes(9)) then
       if (header%idEnergyDens == 0) then
          rdb%nVar = rdb%nVar + 1
          header%idEnergyDens = rdb%nVar
          write(ivard,605) rdb%nVar,'Scaled strain energy density', &
               &           'FORCE/AREA',nbit
       end if
       write(chVars(i:),"('<',i3,'>')") header%idEnergyDens
       i = i + 5
    end if

    ! ---- Write item group
    write(iitem,602) (i,trim(chVars),i=1,nelnod)
    write(iitem,"('    ]'/'  ]'/']')")

600 format('[',i3,';"',a,'";'/'  [;"Element nodes";'/'    [;"Basic";')
601 format('<',i3,';"',a,'";',a,';FLOAT;',i2,';',a,';(',i1,');((',a,'))>')
602 format('      [;',i2,';',a,']')
605 format('<',i3,';"',a,'";',a,';FLOAT;',i2,';SCALAR>')

  end subroutine writeSolidHeader


  !!============================================================================
  !> @brief Writes stress tensor definition to the temporary header files.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 1 Feb 2021

  subroutine writeStressTensorDef (idStress,nVar,iDim,nbit)

    use RDBModule, only : ivard

    integer, intent(inout) :: idStress, nVar
    integer, intent(in)    :: iDim, nbit

    !! --- Logic section ---

    if (idStress > 0) return

    nVar     = nVar + 1
    idStress = nVar
    select case (iDim)
    case (1)
       write(ivard,601) nVar,'Stress','FORCE/AREA',nbit,'TENSOR1',1,'"sigma_xx"'
    case (2)
       write(ivard,601) nVar,'Stress','FORCE/AREA',nbit,'TENSOR2',3, &
            &           '"sigma_xx","sigma_yy","sigma_xy"'
    case (3)
       write(ivard,601) nVar,'Stress','FORCE/AREA',nbit,'TENSOR3',6, &
            &           '"sigma_xx","sigma_yy","sigma_zz",'// &
            &           '"sigma_xy","sigma_xz","sigma_yz"'
    end select

601 format('<',i3,';"',a,'";',a,';FLOAT;',i2,';',a,';(',i1,');((',a,'))>')

  end subroutine writeStressTensorDef


  !!============================================================================
  !> @brief Writes strain tensor definition to the temporary header files.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 1 Feb 2021

  subroutine writeStrainTensorDef (idStrain,nVar,iDim,nbit)

    use RDBModule, only : ivard

    integer, intent(inout) :: idStrain, nVar
    integer, intent(in)    :: iDim, nbit

    !! --- Logic section ---

    if (idStrain > 0) return

    nVar     = nVar + 1
    idStrain = nVar
    select case (iDim)
    case (1)
       write(ivard,601) nVar,'Strain','NONE',nbit,'TENSOR1',1,'"epsilon_xx"'
    case (2)
       write(ivard,601) nVar,'Strain','NONE',nbit,'TENSOR2',3, &
            &           '"epsilon_xx","epsilon_yy","epsilon_xy"'
    case (3)
       write(ivard,601) nVar,'Strain','NONE',nbit,'TENSOR3',6, &
            &           '"epsilon_xx","epsilon_yy","epsilon_zz",'// &
            &           '"epsilon_xy","epsilon_xz","epsilon_yz"'
    end select

601 format('<',i3,';"',a,'";',a,';FLOAT;',i2,';',a,';(',i1,');((',a,'))>')

  end subroutine writeStrainTensorDef


  !!============================================================================
  !> @brief Writes nodal displacements to the results database file.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 27 Nov 2000

  subroutine writeDisplacementDB (rdb,sam,sv,ierr,iComp,svtot,lDouble2)

    use RecKindModule         , only : rk
    use RDBModule             , only : RDBType, writeRDB
    use SamModule             , only : SamType
    use reportErrorModule     , only : reportError, debugFileOnly_p
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getbool

    type(RDBType)     , intent(inout) :: rdb
    type(SamType)     , intent(in)    :: sam
    real(rk)          , intent(in)    :: sv(:)
    integer           , intent(out)   :: ierr
    integer,  optional, intent(in)    :: iComp
    real(rk), optional, intent(in)    :: svtot(:)
    logical,  optional, intent(in)    :: lDouble2

    !! Local variables
    integer :: i, j, k, lStart, nComp
    logical :: lDouble

    !! --- Logic section ---

    if (present(lDouble2)) then
       lDouble = lDouble2
    else
       call ffa_cmdlinearg_getbool ('double',lDouble)
    end if
    if (present(iComp)) then
       nComp = iComp
    else
       nComp = 1
    end if

    if (lNodes) then

       do i = 1, sam%nnod
          if (sam%mpar(33) > 0 .and. sam%mnnn(i) > 0) cycle ! outside group
          j = sam%madof(i)
          if (sam%madof(i+1) > j+2) then
             if (sam%madof(i+1) > j+5) then
                k = j+5
             else
                k = j+2
             end if
             do lStart = 0, sam%ndof*(nComp-1), sam%ndof
                call writeRDB (rdb,sv(lStart+j:lStart+k),ierr,lDouble)
             end do
             if (present(svtot)) then
                call writeRDB (rdb,svtot(j:k),ierr,lDouble)
             end if
             if (ierr < 0) exit
          end if
       end do

    end if
    if (lVector > 0) then

       !! Write translations only, for efficient loading of animation
       do lStart = 0, (nComp-1)*sam%ndof, sam%ndof
          do i = 1, sam%nnod
             if (sam%minex(i) < 0) cycle ! Skip internal beam nodes
             j = sam%madof(i)
             call writeRDB (rdb,sv(lStart+j:lStart+j+2),ierr,lDouble)
             if (ierr < 0) exit
          end do
          !! Write rotations as separate vectors
          do i = 1, sam%nnod
             if (sam%minex(i) < 0) cycle ! Skip internal beam nodes
             j = sam%madof(i)+3
             if (j+3 > sam%madof(i+1)) cycle ! No rotation in this node
             call writeRDB (rdb,sv(lStart+j:lStart+j+2),ierr,lDouble)
             if (ierr < 0) exit
          end do
       end do

    end if

    if (ierr < 0) call reportError (debugFileOnly_p,'writeDisplacementDB')

  end subroutine writeDisplacementDB


  !!============================================================================
  !> @brief Writes element stresses and strains to the results database file.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 27 Nov 2000

  subroutine writeStressDB (rdb,lStress,lStrain,lDouble, &
       &                    nComp,nStrP,Stress,Strain,ierr)

    use RDBModule        , only : RDBType, writeRDB
    use KindModule       , only : dp
    use reportErrorModule, only : reportError, debugFileOnly_p

    type(RDBType), intent(inout) :: rdb
    logical      , intent(in)    :: lStress, lStrain, lDouble
    integer      , intent(in)    :: nComp, nStrP
    real(dp)     , intent(in)    :: Stress(nComp,nStrP), Strain(nComp,nStrP)
    integer      , intent(out)   :: ierr

    ! Local variables
    integer :: i

    ! --- Logic section ---

    if (lStress .and. lStrain) then

       ! Save both stresses and strains
       do i = 1, nStrP
          call writeRDB (rdb,Stress(:,i),ierr,lDouble)
          call writeRDB (rdb,Strain(:,i),ierr,lDouble)
       end do

    else if (lStress) then

       ! Save stresses only
       call writeRDB (rdb,Stress,ierr,lDouble)

    else if (lStrain) then

       ! Save strains only
       call writeRDB (rdb,Strain,ierr,lDouble)

    end if

    if (ierr < 0) call reportError (debugFileOnly_p,'writeStressDB')

  end subroutine writeStressDB


  !!============================================================================
  !> @brief Writes derived stress/strain measures to the results database file.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 28 Oct 2002

  subroutine writeStrMeasureDB (rdb,lStress,lStrain,lStrRes,lDouble, &
       &                        nComp,nRes,nStrP,Stress,Strain,resMat,ierr)

    use RDBModule        , only : RDBType, writeRDB
    use KindModule       , only : dp
    use reportErrorModule, only : reportError, debugFileOnly_p

    type(RDBType), intent(inout) :: rdb
    logical      , intent(in)    :: lStress, lStrain, lDouble, lStrRes(:)
    integer      , intent(in)    :: nComp, nRes, nStrp
    real(dp)     , intent(in)    :: resMat(nRes,nStrP)
    real(dp)     , intent(in)    :: Stress(nComp,nStrP), Strain(nComp,nStrP)
    integer      , intent(out)   :: ierr

    ! Local variables
    integer :: i, j

    ! --- Logic section ---

    if (.not. any(lStrRes)) then

       ! No derived measures requested, save tensors only
       call writeStressDB (rdb,lStress,lStrain,lDouble, &
            &              nComp,nStrP,Stress,Strain,ierr)

    else

       do i = 1, nStrP
          if (lStress) then

             ! Save the stress tensor
             call writeRDB (rdb,Stress(:,i),ierr,lDouble)

          end if
          if (lStrain) then

             ! Save the strain tensor
             call writeRDB (rdb,Strain(:,i),ierr,lDouble)

          end if
          do j = 1, nRes
             if (lStrRes(j)) then

                ! Save derived quantity
                call writeRDB (rdb,resMat(j,i),ierr,lDouble)

             end if
          end do
       end do

    end if

    if (ierr < 0) call reportError (debugFileOnly_p,'writeStrMeasureDB')

  end subroutine writeStrMeasureDB

end module saveStressModule
