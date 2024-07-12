!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module StrainGageModule

  implicit none

  private :: ReadStrainGageOldData


contains

  subroutine ReadStrainGageData (chName,useFSIformat,linkId,strainRosettes,ierr)

    !!==========================================================================
    !! Administer the input of strain gage data from file.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 6 Feb 2004 / 2.0
    !!==========================================================================

    use StrainRosetteModule, only : StrainElementType, ReadStrainElements
    use IdTypeModule       , only : IdType
    use inputUtilities     , only : iuCopyToScratch
    use fileUtilitiesModule, only : findUnitNumber
    use reportErrorModule  , only : reportError, error_p, debugFileOnly_p

    character(len=*)       , intent(in)  :: chName
    logical                , intent(in)  :: useFSIformat
    type(IdType)           , intent(in)  :: linkId
    type(StrainElementType), pointer     :: strainRosettes(:)
    integer                , intent(out) :: ierr

    !! Local variables
    integer :: infp

    !! --- Logic section ---

    infp = findUnitNumber(50)
    open(infp,IOSTAT=ierr,STATUS='SCRATCH')
    if (ierr /= 0) then
       call reportError (error_p,'Unable to open temporary solver input file')
       goto 900
    end if

    call iuCopyToScratch (infp,chname,ierr)
    if (ierr /= 0) goto 900

    rewind(infp)
    if (useFSIformat) then
       call ReadStrainElements (infp,strainRosettes,linkId%baseId,ierr)
    else
       call ReadStrainGageOldData (infp,strainRosettes,linkId%userId,ierr)
    end if
    close(infp)

900 continue
    if (ierr /= 0) call reportError (debugFileOnly_p,'ReadStrainGageData')

  end subroutine ReadStrainGageData


  subroutine ReadStrainGageOldData (inFile,strainRosettes,thisLinkNumber,err)

    !!==========================================================================
    !! Read strain gage data from the rosette definition file.
    !!
    !! Programmer : Bjorn Haugen
    !! date/rev   : 10 Dec 1998 / 1.0
    !!==========================================================================

    use kindModule             , only : dp, pi_p
    use StrainRosetteModule    , only : StrainElementType
    use StrainRosetteModule    , only : nullifyRosette, checkRosette
    use StrainRosetteModule    , only : singleGage_p, doubleGage_p
    use StrainRosetteModule    , only : tripleGage1_p, tripleGage2_p
    use isoMatModule           , only : isoMat2D
    use reportErrorModule      , only : allocationError, getErrorFile
    use reportErrorModule      , only : reportError, errorFileOnly_p, error_p
    use FFlLinkHandlerInterface, only : ffl_ext2int
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_getint

    integer                , intent(in)  :: inFile, thisLinkNumber
    type(StrainElementType), pointer     :: strainRosettes(:)
    integer                , intent(out) :: err

    !! Local variables
    integer  :: i, id, idIn, lpu, rosetteType, linkNumber, rdbinc
    integer  :: numberOfRosettes, numberOfNodes, numberOfGages
    real(dp) :: alpha, Emod, Poisson, length, xVector(3), zVector(3)

    integer, parameter :: ncharcommand = 3
    integer, parameter :: linelength   = 132 ! Max. line length FII can handle
    character(len=128) :: errmsg

    !! External functions
    logical , external :: lnxfii, linfii, lflfii
    integer , external :: innfii, ierfii
    real(dp), external :: dflfii

    !! --- Logic section ---

    err = 0
    lpu = getErrorFile()
    errmsg = ''
    call ffa_cmdlinearg_getint ('rdbinc',rdbinc) ! For creation of dummy baseID

    !! Open FII (Fortrain Input Interpretor)
    rewind(inFile)
    call fiiopn(lineLength,ncharCommand,inFile,lpu,-1,err)
    call fiisst('#')

    !! Count number of straingages and strainRosettes for current link
    numberOfRosettes = 0
    do while (readNextRosette(id,rosetteType,linkNumber))
       if (linkNumber == thisLinkNumber) numberOfRosettes = numberOfRosettes + 1
    end do
    if (err < 0 .or. numberOfRosettes < 1) goto 900

    !! Allocate new array
    allocate(strainRosettes(numberOfRosettes), STAT=err)
    if (err /= 0) then
       err = allocationError('ReadStrainGageOldData 1')
       return
    end if

    !! Now actually read the straingage definitions
    rewind(inFile)
    call fiiopn(lineLength,ncharCommand,inFile,lpu,-1,err)
    call fiisst('#')

    idIn = 0
    do while (readNextRosette(id,rosetteType,linkNumber))
       if (linkNumber /= thisLinkNumber) cycle ! Skip for all other links

       idIn = idIn + 1
       call nullifyRosette (strainRosettes(idIn))

       strainRosettes(idIn)%id%userId = id
       !! Create a unique dummy base id to be used on the frs-file
       strainRosettes(idIn)%id%baseId = idIn + linkNumber*1000 + rdbinc*10000000
       strainRosettes(idIn)%linkNumber = linkNumber

       numberOfNodes = readInteger()
       if (err < 0) goto 900

       allocate(strainRosettes(idIn)%globalNodes(numberOfNodes), &
            &   strainRosettes(idIn)%data(1), STAT=err)
       if (err /= 0) then
          err = allocationError('ReadStrainGageOldData 2')
          return
       end if

       call nullifyRosette (strainRosettes(idIn)%data(1))

       do i = 1, numberOfNodes
          strainRosettes(idIn)%globalNodes(i) = readInteger()
       end do
       strainRosettes(idIn)%data(1)%zPos = readDP()
       do i = 1,3
          xVector(i) = readDP()
       end do
       do i = 1,3
          zVector(i) = readDP()
       end do
       if (err < 0) goto 900

       !! Normalize rosette direction definition
       length = sqrt(sum(xVector**2))
       if (length < 1000.0_dp*tiny(xVector(1))) then
          write(errmsg,"('Undefined x-direction for rosette number:',I8)") id
          goto 900
       else
          xVector = xVector/length
       end if

       length = sqrt(sum(zVector**2))
       if (length < 1000.0_dp*tiny(zVector(1))) then
          write(errmsg,"('Undefined z-direction for rosette number:',I8)") id
          goto 900
       else
          zVector = zVector/length
       end if

       strainRosettes(idIn)%posInGl(:,1) = xVector
       strainRosettes(idIn)%posInGl(:,2) = 0.0_dp
       strainRosettes(idIn)%posInGl(:,3) = zVector
       strainRosettes(idIn)%posInGl(:,4) = 0.0_dp

       Emod    = readDP()
       Poisson = readDP()
       if (err < 0) goto 900

       call isoMat2D (Emod,Poisson,strainRosettes(idIn)%data(1)%Cmat)

       !! Check that the nodal ordering of the strain rosette yields
       !! a proper normal vector and swap the element nodes if not.
       call checkRosette (strainRosettes(idIn), &
            &             strainRosettes(idIn)%globalNodes,err)
       if (err < 0) goto 900

       select case (rosetteType)

       case (singleGage_p)
          numberOfGages = 1
          alpha = 0.0_dp

       case (doubleGage_p)
          numberOfGages = 2
          alpha = pi_p/2.0_dp ! 90 degrees between the gages

       case (tripleGage1_p)
          numberOfGages = 3
          alpha = pi_p/3.0_dp ! 60 degrees between the gages

       case (tripleGage2_p)
          numberOfGages = 3
          alpha = pi_p/4.0_dp ! 45 degrees between the gages

       case default
          write(errmsg,"('Undefined rosette-type:',I8)") rosetteType
          goto 900

       end select

       strainRosettes(idIn)%data(1)%alphaGages = alpha
       allocate(strainRosettes(idIn)%data(1)%gages(numberOfGages),STAT=err)
       if (err /= 0) then
          err = allocationError('ReadStrainGageOldData 3')
          return
       end if

    end do

900 if (errmsg /= '') then
       err = err-1
       call reportError (error_p,errmsg,addString='ReadStrainGageOldData')
    else if (err < 0) then
       call reportError (error_p,'Failed to read rosette definition file', &
            &            addString='ReadStrainGageOldData')
    end if

  contains

    function readInteger()
      integer :: readInteger
      if (linfii()) then
         readInteger = innfii()
      else
         readInteger = 0
         err = err - 1
         call reportError (errorFileOnly_p,'Failure reading an integer value')
      end if
    end function readInteger

    function readDP()
      real(dp) :: readDP
      if (lflfii()) then
         readDP = dflfii()
      else
         readDP = 0.0_dp
         err = err - 1
         call reportError (errorFileOnly_p,'Failure reading a real value')
      end if
    end function readDP

    function readNextRosette (rId,rType,linkId)
      integer, intent(out) :: rId,rType,linkId
      logical :: readNextRosette

      readNextRosette = .false.
      if (err < 0) return

      call fiirln()
      err = ierfii()
      if (err < 0) return

      if (lnxfii('END')) return
      if (lnxfii('EOF')) return

      rId    = readInteger()
      rType  = readInteger()
      linkId = readInteger()
      if (err == 0) readNextRosette = .true.

    end function readNextRosette

  end subroutine ReadStrainGageOldData


  subroutine PrintStrainGages_FiDF (rosette,time)

    !!==========================================================================
    !! Write gage strains to DAC files.
    !!
    !! Programmer : Bjorn Haugen
    !! date/rev   : 10 Dec 1998 / 1.0
    !!==========================================================================

    use StrainRosetteModule      , only : StrainRosetteType, dp
    use FiDeviceFunctionInterface, only : FiDF_SetValue

    type(StrainRosetteType), intent(in) :: rosette
    real(dp)               , intent(in) :: time

    !! Local variables
    integer :: i, ifile

    !! --- Logic section ---

    do i = 1, size(rosette%gages)
       iFile = rosette%gages(i)%lpu_FiDF
       if (iFile < 1) cycle

       !! Write in microstrains, i.e., multiply with 1.0e6
       call FiDF_SetValue (iFile,time,rosette%gages(i)%epsGage*1.0e6_dp)
    end do

  end subroutine PrintStrainGages_FiDF


  subroutine InitStrainGages (rosId,gages,alphaGages,posInGl,dT_FiDF,ierr)

    !!==========================================================================
    !! Initiate the strain-gage data and open the device function file
    !! for each strain gage in the strain rosette.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 25 May 2001 / 2.0
    !!==========================================================================

    use kindModule               , only : dp, epsDiv0_p
    use StrainRosetteModule      , only : StrainGageType, InitStrainGauges
    use reportErrorModule        , only : reportError, error_p, warning_p
    use FiDeviceFunctionInterface, only : FiDF_OpenWrite, FiDF_SetStep
    use FiDeviceFunctionInterface, only : FiDF_SetXAxis, FiDF_SetYAxis
    use FiDeviceFunctionInterface, only : DACFile_p

    integer             , intent(in)    :: rosId
    type(StrainGageType), intent(inout) :: gages(:)
    real(dp)            , intent(in)    :: alphaGages, posInGl(:,:), dT_FiDF
    integer             , intent(out)   :: ierr

    !! Local variables
    integer :: i

    character(len=1)  :: chGage
    character(len=16) :: chId
    character(len=64) :: chName

    character(len=7), parameter :: chrUnit1  = 'Seconds'
    character(len=2), parameter :: chrUnit2  = 'ue'
    character(len=4), parameter :: chrTitle1 = 'Time'
    character(len=6), parameter :: chrTitle2 = 'Strain'

    !! --- Logic section ---

    ierr = 0
    call initStrainGauges (gages,alphaGages,posInGl)
    if (dT_FiDF < 0.0_dp) return ! Suppress DAC-file output

    !! Initialize a DAC-file for each strain gage
    do i = 1, size(gages)

       write(chId,'(I12)') rosId
       write(chGage,'(I1)') i
       chName = 'rosette'//trim(adjustl(chId))//'_gage'//chGage//'.dac'
       call FiDF_OpenWrite (trim(chName),DACFile_p,gages(i)%lpu_FiDF,ierr)
       if (ierr /= 0) then
          call reportError (error_p,'Unable to open DAC-file '//chName, &
               &            addString='InitStrainGages')
          return
       end if

       !! Initialize axis information
       if (dT_FiDF < epsDiv0_p) then
          call reportError (warning_p,'Zero step-size for DAC-file '//chName)
       else
          call FiDF_SetStep (gages(i)%lpu_FiDF,dT_FiDF)
       end if
       call FiDF_SetXAxis (gages(i)%lpu_FiDF,chrTitle1,ChrUnit1)
       call FiDF_SetYAxis (gages(i)%lpu_FiDF,chrTitle2,ChrUnit2)

    end do

  end subroutine InitStrainGages


  subroutine addFatiguePoints (rosette,time,toMPaScale)

    !!==========================================================================
    !! Add current stress state of the strain gages to the fatigue calculator.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 30 May 2007/1.0
    !!==========================================================================

    use StrainRosetteModule, only : StrainRosetteType, dp
    use FFpFatigueInterface, only : ffp_addPoint

    type(StrainRosetteType), intent(inout) :: rosette
    real(dp)               , intent(in)    :: time, toMPaScale

    !! Local variables
    integer :: i

    !! --- Logic section ---

    call ffp_addPoint (rosette%fatHandle, time, rosette%sigmaP(1)*toMPaScale)

    do i = 1, size(rosette%gages)
       call ffp_addPoint (rosette%gages(i)%fatHandle, time, &
            &             rosette%gages(i)%sigGage*toMPaScale)
    end do

  end subroutine addFatiguePoints


  subroutine calcGageDamage (gages,gateValue,snData,damage,printCyc)

    !!==========================================================================
    !! Compute damage in the strain gages based on the given S-N curve data.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 30 May 2007/1.0
    !!==========================================================================

    use StrainRosetteModule, only : StrainGageType, dp
    use FFpFatigueInterface, only : ffp_getDamage

    type(StrainGageType), intent(in)  :: gages(:)
    real(dp)            , intent(in)  :: gateValue, snData(:)
    real(dp)            , intent(out) :: damage(:)
    integer, optional   , intent(in)  :: printCyc

    !! Local variables
    integer :: i

    !! --- Logic section ---

    do i = 1, size(gages)
       damage(1+i) = ffp_getDamage(gages(i)%fatHandle,gateValue,snData,printCyc)
    end do

  end subroutine calcGageDamage


  subroutine calcGageNCycles (gages,low,high,nCycles)

    !!==========================================================================
    !! Compute the number of stress cycles in the range [low,high].
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 4 Jul 2007/1.0
    !!==========================================================================

    use StrainRosetteModule, only : StrainGageType, dp
    use FFpFatigueInterface, only : ffp_getNumCycles

    type(StrainGageType), intent(in)  :: gages(:)
    real(dp)            , intent(in)  :: low, high
    integer             , intent(out) :: nCycles(:)

    !! Local variables
    integer :: i

    !! --- Logic section ---

    do i = 1, size(gages)
       nCycles(1+i) = ffp_getNumCycles(gages(i)%fatHandle,low,high)
    end do

  end subroutine calcGageNCycles


  subroutine reportDamage (rosettes,iFatigue,lpu)

    !!==========================================================================
    !! Calculate and report the damage in all strain rosettes.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 30 May 2007/1.0
    !!==========================================================================

    use StrainRosetteModule   , only : StrainElementType, dp
    use IdTypeModule          , only : getID
    use FFpFatigueInterface   , only : ffp_getDamage, ffp_getNumCycles
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getdouble

    type(StrainElementType), intent(in) :: rosettes(:)
    integer                , intent(in) :: iFatigue, lpu

    !! Local variables
    integer           :: i, j, n
    real(dp)          :: gate, s0, s1, binSize, sndata(4), damage(4)
    integer           :: nCycles(4)
    character(len=48) :: cCycles

    !! --- Logic section ---

    nCycles = -1
    call ffa_cmdlinearg_getdouble('binSize',binSize)

    do i = 1, size(rosettes)
       n = size(rosettes(i)%data(1)%gages)

       !! Use default values from command-line arguments if not given
       gate = rosettes(i)%data(1)%gateValue
       sndata = rosettes(i)%data(1)%snCurve
       if (gate <= 0.0_dp)      call ffa_cmdlinearg_getdouble('gate',gate)
       if (sndata(1) <= 0.0_dp) call ffa_cmdlinearg_getdouble('loga1',sndata(1))
       if (sndata(2) <= 0.0_dp) call ffa_cmdlinearg_getdouble('loga2',sndata(2))
       if (sndata(3) <= 0.0_dp) call ffa_cmdlinearg_getdouble('m1',sndata(3))
       write(lpu,600) gate,sndata

       !! Calculate damage
       if (iFatigue == 2) then ! Print debug pvx and cycles files
          damage(1) = ffp_getDamage(rosettes(i)%data(1)%fatHandle,gate,sndata,1)
          call calcGageDamage (rosettes(i)%data(1)%gages,gate,sndata,damage,1)
       else
          damage(1) = ffp_getDamage(rosettes(i)%data(1)%fatHandle,gate,sndata)
          call calcGageDamage (rosettes(i)%data(1)%gages,gate,sndata,damage)
       end if
       write(lpu,610) getId(rosettes(i)%id),damage(1:1+n)

       !! Count the stress cycles
       s0 = 0.0_dp
       s1 = binSize
       nCycles(1) = ffp_getNumCycles(rosettes(i)%data(1)%fatHandle,s0,s1)
       call calcGageNCycles (rosettes(i)%data(1)%gages,s0,s1,nCycles)
       do while (any(nCycles >= 0))
          if (any(nCycles > 0)) then
             do j = 1, 1+n
                if (nCycles(j) >= 0) then
                   write(cCycles(12*j-11:12*j),"(I12)") nCycles(j)
                else
                   cCycles(12*j-11:12*j) = ' '
                end if
             end do
             write(lpu,620) s0,s1,trim(cCycles)
          end if
          s0 = s1
          s1 = s0 + binSize
          nCycles(1) = ffp_getNumCycles(rosettes(i)%data(1)%fatHandle,s0,s1)
          call calcGageNCycles (rosettes(i)%data(1)%gages,s0,s1,nCycles)
       end do

       write(lpu,690)
    end do

600 format(/5X,'===== Computed damage in strain rosette =====' &
         &/11X,'gate value :',1PE12.5, &
         &/11X,'S-N curve  : log(a1) =',0PF7.3,' log(a2) =',F7.3, &
         &     ' m1 =',F7.3,' m2 =',F7.3, &
         &//5X,'Strain Rosette',22X,'Max princ.  Gage 1      Gage 2  ...')
610 format( 4X,A36,1P4E12.5 /)
620 format( 5X,'Stress cycles ',F7.2,' -',F7.2,A)
690 format( 5X,'============================================='/ )

  end subroutine reportDamage


  subroutine deallocateFatigue (rosette)

    !!==========================================================================
    !! Release all fatigue data associated with the given strain rosette.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 30 May 2007/1.0
    !!==========================================================================

    use StrainRosetteModule, only : StrainRosetteType
    use FFpFatigueInterface, only : ffp_releaseData

    type(StrainRosetteType), intent(inout) :: rosette

    !! Local variables
    integer :: i

    !! --- Logic section ---

    call ffp_releaseData (rosette%fatHandle)
    do i = 1, size(rosette%gages)
       call ffp_releaseData (rosette%gages(i)%fatHandle)
    end do

  end subroutine deallocateFatigue

end module StrainGageModule
