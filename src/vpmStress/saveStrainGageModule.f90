!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module saveStrainGageModule

  implicit none

  private :: writeRosetteHeader


contains

  subroutine writeStrainGageHeader (resFileName,modelFileName,rdb,rosettes, &
       &                            bufRatio,ierr,initRDB,partId)

    !!==========================================================================
    !! Administer writing of results database header for all strain gages.
    !!
    !! Programmer : Tommy Stokmo Jørstad
    !! date/rev   : 27 Sept 2001/1.0
    !!==========================================================================

    use KindModule            , only : dp, i8, lfnam_p
    use RDBModule             , only : RDBType, nullifyRDB, openRDBfile
    use RDBModule             , only : openHeaderFiles, writeTimeStepHeader
    use StrainRosetteModule   , only : StrainElementType
    use reportErrorModule     , only : reportError, debugFileOnly_p
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getbool
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getstring

    character(len=*)       , intent(in)  :: resFileName, modelFileName
    type(RDBType)          , intent(out) :: rdb
    type(StrainElementType), intent(in)  :: rosettes(:)
    real(dp)               , intent(in)  :: bufRatio
    integer                , intent(out) :: ierr
    logical, optional      , intent(in)  :: initRDB
    integer, optional      , intent(in)  :: partId

    !! Local variables
    logical            :: lDef
    character(lfnam_p) :: linkFileName

    !! --- Logic section ---

    if (.not.present(initRDB)) then
       call nullifyRDB (rdb)
    else if (initRDB) then
       call nullifyRDB (rdb)
    end if

    !! Open the temporary header files
    if (present(partId)) then
       lDef = .false.
       call openHeaderFiles ('fedem_solver',ierr,'response data base file', &
            &                modelFileName)
    else
       call ffa_cmdlinearg_getbool ('deformation',lDef)
       call ffa_cmdlinearg_getstring ('linkfile',linkFileName)
       call openHeaderFiles ('fedem_gage',ierr,'response data base file', &
            &                modelFileName,linkFileName)
    end if
    if (ierr < 0) goto 500

    !! Write results data definitions to header files
    call writeTimeStepHeader (rdb)
    call writeRosetteHeader (rdb,rosettes,lDef,partId)

    !! Set buffer size
    rdb%nBstep  = int(rdb%nBytes)
    rdb%bufSize = int(rdb%nBytes*bufRatio)
    rdb%nBytes  = 0_i8

    !! Open the results database file and copy the file header to it
    call openRDBfile (rdb,ierr,resFileName,'#FEDEM response data')

500 if (ierr < 0) call reportError (debugFileOnly_p,'writeStrainGageHeader')

  end subroutine writeStrainGageHeader


  subroutine writeRosetteHeader (rdb,rosettes,lDef,partId)

    !!==========================================================================
    !! Administer writing of results database header for strain rosettes.
    !!
    !! Programmer : Tommy Stokmo Jørstad
    !! date/rev   : 27 Sept 2001/1.0
    !!==========================================================================

    use RDBModule          , only : RDBType, idatd, nbs_p
    use RDBModule          , only : writeVarDef, writeItgDef, idatd
    use IdTypeModule       , only : writeIdHeader, StrId
    use StrainRosetteModule, only : StrainElementType

    type(RDBType)          , intent(inout) :: rdb
    type(StrainElementType), intent(in)    :: rosettes(:)
    logical                , intent(in)    :: lDef
    integer, optional      , intent(in)    :: partId

    !! Local variables
    integer :: i, j, iFile, nGages, nNodes
    integer :: idPos, idAngle, idDef
    integer :: idEpsC, idSigmaC
    integer :: idAlpha1, idAlphaGamma
    integer :: idEpsGage, idSigGage
    integer :: idNode, idStrainGages(3)

    !! --- Logic section ---

    idPos         = 0
    idAngle       = 0
    idDef         = 0
    idEpsC        = 0
    idSigmaC      = 0
    idAlpha1      = 0
    idAlphaGamma  = 0
    idEpsGage     = 0
    idSigGage     = 0
    idStrainGages = 0

    do i = 1, size(rosettes)
       if (present(partId)) then
          if (rosettes(i)%linkNumber /= partId) cycle
       end if

       nGages = size(rosettes(i)%data(1)%gages)
       rdb%nBytes = rdb%nBytes + (8+2*nGages)*nbs_p
       if ( lDef .and. associated(rosettes(i)%disp) .and. &
            associated(rosettes(i)%globalNodes) ) then
          nNodes = size(rosettes(i)%globalNodes)
          rdb%nBytes = rdb%nBytes + 3*nNodes*nbs_p
       else
          nNodes = 0
       end if
       if (associated(rosettes(i)%ur)) then
          rdb%nBytes = rdb%nBytes + 6*nbs_p
       end if

       call writeIdHeader ('Strain rosette',rosettes(i)%id,idatd)

       call writeVarDef (rdb,idAlpha1,idatd, &
            &            'Angle of maximum principal strain/stress','ANGLE')
       call writeVarDef (rdb,idAlphaGamma,idatd, &
            &            'Angle of maximum shear','ANGLE')
       call writeVarDef (rdb,idEpsC,idatd, &
            &            'Strain tensor','NONE','TENSOR2', &
            &            '"epsilon_x","epsilon_y","epsilon_xy"')
       call writeVarDef (rdb,idSigmaC,idatd, &
            &            'Stress tensor','FORCE/AREA','TENSOR2', &
            &            '"sigma_x","sigma_y","sigma_xy"')

       do j = 1, nGages
          call writeItGDef (rdb,idStrainGages(j),iFile, &
               &            'Gage '//char(ichar('0')+j))
          call writeVarDef (rdb,idEpsGage,iFile,'Gage strain')
          call writeVarDef (rdb,idSigGage,iFile,'Gage stress')
          if (iFile > 0) write(iFile,"(']')")
       end do

       do j = 1, nNodes
          idNode = 0
          call writeItGDef (rdb,idNode,iFile, &
               &            'Node'//trim(StrId(rosettes(i)%globalNodes(j))))
          call writeVarDef (rdb,idDef,iFile,'Deformation', &
               &            'LENGTH','VEC3','"d_x","d_y","d_z"')
          write(iFile,"(']')")
       end do

       if (associated(rosettes(i)%ur)) then
          call writeVarDef (rdb,idPos,idatd,'Position', &
               &            'LENGTH','VEC3','"x","y","z"')
          call writeVarDef (rdb,idAngle,idatd,'Euler angles', &
               &            'ANGLE','ROT3','"eps_x","eps_y","eps_z"')
       end if

       write(idatd,"('}')")

    end do

  end subroutine writeRosetteHeader


  subroutine writeStrainGageDB (respDB,rosettes,iStep,time,lDef,ierr,partId)

    !!==========================================================================
    !! Administer writing of results database for strain gages and rosettes.
    !!
    !! Programmer : Tommy Stokmo Jørstad
    !! date/rev   : 27 Sept 2001/1.0
    !!==========================================================================

    use KindModule         , only : dp, i8
    use RDBModule          , only : RDBType, checkStepSize
    use RDBModule          , only : writeTimeStepDB, writeRDB
    use StrainRosetteModule, only : StrainElementType
    use reportErrorModule  , only : reportError, debugFileOnly_p

    type(RDBType)          , intent(inout) :: respDB
    type(StrainElementType), intent(in)    :: rosettes(:)
    integer(i8)            , intent(in)    :: iStep
    real(dp)               , intent(in)    :: time
    logical                , intent(in)    :: lDef
    integer                , intent(out)   :: ierr
    integer, optional      , intent(in)    :: partId

    !! Local variables
    integer  :: i, j, nBytesPerStep = 0
    real(dp) :: epsTensorial(3)

    !! --- Logic section ---

    if (nBytesPerStep == 0) nBytesPerStep = -int(respDB%nBytes)

    call writeTimeStepDB (respDB,iStep,time,ierr)
    if (ierr < 0) goto 500

    do i = 1, size(rosettes)
       if (present(partId)) then
          if (rosettes(i)%linkNumber /= partId) cycle
       end if

       !! Replace the engineering shear strain by the corresponding tensorial
       !! quantity before writing the results, such that the UI is able to
       !! compute correct von Mises and principal values
       epsTensorial(1:2) = rosettes(i)%data(1)%epsC(1:2)
       epsTensorial(3) = 0.5_dp * rosettes(i)%data(1)%epsC(3)
       call writeRDB (respDB,rosettes(i)%data(1)%alpha1,ierr)
       call writeRDB (respDB,rosettes(i)%data(1)%alphaGamma,ierr)
       call writeRDB (respDB,epsTensorial,ierr)
       call writeRDB (respDB,rosettes(i)%data(1)%sigmaC,ierr)
       if (ierr < 0) goto 500

       do j = 1, size(rosettes(i)%data(1)%gages)
          call writeRDB (respDB,rosettes(i)%data(1)%gages(j)%epsGage,ierr)
          call writeRDB (respDB,rosettes(i)%data(1)%gages(j)%sigGage,ierr)
          if (ierr < 0) goto 500
       end do

       if (lDef .and. associated(rosettes(i)%disp)) then
          call writeRDB (respDB,rosettes(i)%disp,ierr)
          if (ierr < 0) goto 500
       end if

       if (associated(rosettes(i)%ur)) then
          call writeRDB (respDB,rosettes(i)%ur,ierr)
          if (ierr < 0) goto 500
       end if

    end do

    if (nBytesPerStep < 1) then
       !! Check consistency between header and the actual data
       nBytesPerStep = nBytesPerStep + int(respDB%nBytes)
       call checkStepSize (respDB,nBytesPerStep,ierr)
       if (ierr < 0) goto 500
    end if

    ierr = 0
    return

500 call reportError (debugFileOnly_p,'writeStrainGageDB')

  end subroutine writeStrainGageDB

end module saveStrainGageModule
