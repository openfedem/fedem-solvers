!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module StrainCoatModule

  use KindModule         , only : sp, dp
  use StrainRosetteModule, only : StrainElementType

  implicit none

  type AngleBinType
     integer  :: nVal
     real(dp) :: sigMax, sigMin
     real(dp) :: epsMax, epsMin
  end type AngleBinType

  type AngleBinPtr
     type(AngleBinType), pointer :: p
  end type AngleBinPtr


  type CoatResultType
     real(sp)         , pointer :: fppValue(:)
     type(AngleBinPtr), pointer :: angBin(:)
     integer  :: SNcurve(2)
     real(dp) :: sCF, fatValue, damage
     real(dp) :: sigMax, sigMin, tauMax, vmsMax
     real(dp) :: epsMax, epsMin, gammaMax, vmeMax
     real(dp) :: strRange(2)
     real(dp) :: biAxialSum, biAxialSqr
     real(dp) :: popAngle, angSpread
     integer  :: nBiAxial
     integer  :: matGroup
     integer  :: resultSet
  end type CoatResultType


  type StrainCoatType
     type(StrainElementType)       :: tensorRosette
     type(CoatResultType), pointer :: results(:)
     integer                       :: fppType
     integer             , pointer :: fppProcessorHandle(:)
  end type StrainCoatType


  integer, save :: iBuf = 0 !! Time step buffer counter

  integer, parameter :: maxCoats_p = 3
  integer, parameter :: maxNodes_p = 9

  !! Result set names for frs-file output
  character(len=6), parameter :: resultSet_p(0:maxCoats_p) = (/ &
       &                         'Basic ', &
       &                         'Bottom', &
       &                         'Mid   ', &
       &                         'Top   ' /)

  !! Element type names for frs-file output
  character(len=6), parameter :: elmType_p(3:maxNodes_p) = (/ &
       &                         'STRCT3', &
       &                         'STRCQ4', &
       &                         '      ', &
       &                         'STRCT6', &
       &                         '      ', &
       &                         'STRCQ8', &
       &                         'STRCQ9' /)


contains

  subroutine deallocateCoat (coat)

    !!==========================================================================
    !! Release the allocated memory of a StrainCoatType object.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 30 May 2001/1.0
    !!==========================================================================

    use StrainRosetteModule, only : deallocateRosette

    type(StrainCoatType), intent(inout) :: coat

    !! Local variables
    integer :: i, j

    !! --- Logic section ---

    if (associated(coat%results)) then
       do i = 1, size(coat%results)
          if (associated(coat%results(i)%fppValue)) then
             deallocate(coat%results(i)%fppValue)
          end if
          if (associated(coat%results(i)%angBin)) then
             do j = 1, size(coat%results(i)%angBin)
                if (associated(coat%results(i)%angBin(j)%p)) then
                   deallocate(coat%results(i)%angBin(j)%p)
                end if
             end do
             deallocate(coat%results(i)%angBin)
          end if
       end do
       deallocate(coat%results)
    end if

    if (associated(coat%fppProcessorHandle)) then
       deallocate(coat%fppProcessorHandle)
    end if

    call deallocateRosette (coat%tensorRosette)

  end subroutine deallocateCoat


  subroutine nullifyCoat (coat)

    !!==========================================================================
    !! Initialize the StrainCoatType object.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 30 May 2001/1.0
    !!==========================================================================

    use StrainRosetteModule, only : nullifyRosette

    type(StrainCoatType), intent(out) :: coat

    !! --- Logic section ---

    coat%fppType = 0
    nullify(coat%results)
    nullify(coat%fppProcessorHandle)
    call nullifyRosette(coat%tensorRosette)

  end subroutine nullifyCoat


  subroutine nullifyResults (coat)

    !!==========================================================================
    !! Initialize the CoatResultType object.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 30 May 2001/1.0
    !!==========================================================================

    use KindModule, only : hugeVal_p

    type(CoatResultType), intent(out) :: coat

    !! --- Logic section ---

    nullify(coat%fppValue)
    nullify(coat%angBin)

    coat%sigMax     = 0.0_dp
    coat%tauMax     = 0.0_dp
    coat%vmsMax     = 0.0_dp
    coat%epsMax     = 0.0_dp
    coat%gammaMax   = 0.0_dp
    coat%vmeMax     = 0.0_dp
    coat%sigMin     = hugeVal_p
    coat%epsMin     = hugeVal_p
    coat%SNcurve    = 0
    coat%sCF        = 1.0_dp
    coat%fatValue   = 0.0_dp
    coat%damage     = 0.0_dp
    coat%strRange   = 0.0_dp
    coat%biAxialSum = 0.0_dp
    coat%biAxialSqr = 0.0_dp
    coat%popAngle   = 0.0_dp
    coat%angSpread  = 0.0_dp
    coat%nBiAxial   = 0
    coat%matGroup   = 0
    coat%resultSet  = 0

  end subroutine nullifyResults


  subroutine initiateStrainCoats (linkId,iSurface,fppType,bufSize,angBinSize, &
       &                          strainCoats,nStrainCoats,ierr)

    !!==========================================================================
    !! Initiate all strain coat elements with data from the link object.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 30 May 2001/1.0
    !!==========================================================================

    use isoMatModule           , only : isoMat2D
    use StrainRosetteModule    , only : nullifyRosette
    use reportErrorModule      , only : internalError, allocationError
    use reportErrorModule      , only : reportError, error_p
    use FFlLinkHandlerInterface, only : ffl_getstraincoat, ffl_ext2int

    integer             , intent(in)  :: linkId, iSurface, fppType
    integer             , intent(in)  :: bufSize, angBinSize
    type(StrainCoatType), intent(out) :: strainCoats(:)
    integer             , intent(out) :: nStrainCoats
    integer             , intent(out) :: ierr

    !! Local variables
    logical, parameter :: elmID_p    = .false.
    integer, parameter :: maxnod_p   = 8
    integer, parameter :: maxpoint_p = 3
    integer            :: i, j, k, id, ip, elmId
    integer            :: nNodes, nPoints, nSurf, iCoat = 0
    integer            :: nodes(maxnod_p), mNum(maxpoint_p), rset(maxpoint_p)
    integer            :: SNcurve(2,maxpoint_p)
    real(dp)           :: E(maxpoint_p), nu(maxpoint_p), Zpos(maxpoint_p)
    real(dp)           :: SCF(maxpoint_p)
    character(len=64)  :: errMsg

    !! --- Logic section ---

    nStrainCoats = size(strainCoats)

    do i = 1, nStrainCoats

       call ffl_getstraincoat (id,nNodes,nPoints,nodes, &
            &                  mNum,E,nu,Zpos,rSet,SCF,SNcurve,elmId,ierr)
       if (ierr < 0) then
          ierr = internalError('ffl_getStrainCoat')
          return
       else if (ierr > 0) then
          nStrainCoats = i - 1 ! No more strain coat elements in this link
          return
       else if (nNodes > maxnod_p .or. nPoints > maxpoint_p) then
          write(errMsg,"('Too many nodes or points:',2i6)") nNodes, nPoints
          call reportError (error_p,'Invalid strain coat definition.',errMsg, &
               &            addString='initiateStrainCoats')
          ierr = -1
          return
       end if

       if (iSurface > 0) then
          nSurf = 0
          do j = 1, nPoints
             if (iSurface == rSet(j)) nSurf = nSurf + 1
          end do
       else
          nSurf = nPoints
       end if

       iCoat = iCoat + 1
       call nullifyCoat (strainCoats(i))
       strainCoats(i)%tensorRosette%id%userId = id
       strainCoats(i)%tensorRosette%id%baseId = iCoat
       strainCoats(i)%tensorRosette%elmNumber = ffl_ext2int(elmID_p,elmId)
       strainCoats(i)%tensorRosette%linkNumber = linkId
       if (fppType < 0 .or. minval(SNcurve(:,1:nPoints)) >= 0) then
          strainCoats(i)%fppType = fppType
       else
          strainCoats(i)%fppType = 0
       end if

       allocate(strainCoats(i)%tensorRosette%globalNodes(nNodes), &
            &   strainCoats(i)%tensorRosette%data(nSurf), &
            &   strainCoats(i)%results(nSurf), &
            &   STAT=ierr)
       if (ierr /= 0) then
          ierr = allocationError('initiateStrainCoats 1')
          return
       end if

       do j = 1, nNodes
          strainCoats(i)%tensorRosette%globalNodes(j) = nodes(j)
       end do

       if (strainCoats(i)%fppType /= 0) then
          allocate(strainCoats(i)%fppProcessorHandle(nSurf), STAT=ierr)
          if (ierr /= 0) then
             ierr = allocationError('initiateStrainCoats 1a')
             return
          end if
          if (nSurf > 0) strainCoats(i)%fppProcessorHandle = 0
       end if

       ip = 0
       do j = 1, nPoints
          if (iSurface > 0 .and. iSurface /= rSet(j)) cycle
          ip = ip + 1
          call nullifyRosette (strainCoats(i)%tensorRosette%data(ip))
          call nullifyResults (strainCoats(i)%results(ip))
          call isoMat2D (E(j),nu(j),strainCoats(i)%tensorRosette%data(ip)%Cmat)
          strainCoats(i)%tensorRosette%data(ip)%zPos = Zpos(j)
          strainCoats(i)%results(ip)%matGroup = mNum(j)
          strainCoats(i)%results(ip)%resultSet = rSet(j)
          strainCoats(i)%results(ip)%SNcurve = SNcurve(:,j)
          strainCoats(i)%results(ip)%sCF = SCF(j)
          if (strainCoats(i)%fppType < 0) then
             allocate(strainCoats(i)%results(ip)%fppValue(bufSize), &
                  &   strainCoats(i)%results(ip)%angBin(angBinSize-1),STAT=ierr)
          else
             allocate(strainCoats(i)%results(ip)%angBin(angBinSize-1),STAT=ierr)
          end if
          if (ierr /= 0) then
             ierr = allocationError('initiateStrainCoats 2')
             return
          end if
          do k = 1, angBinSize-1
             nullify(strainCoats(i)%results(ip)%angBin(k)%p)
          end do
       end do

    end do

  end subroutine initiateStrainCoats


  subroutine calcStrainCoatData (strainCoat,displ,biAxialGate,toMPaScale,ierr)

    !!==========================================================================
    !! Update the strain coat element based on the nodal displacement vector.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 29 May 2001/2.0
    !!==========================================================================

    use KindModule         , only : pi_p
    use StrainRosetteModule, only : calcRosetteStrains
    use reportErrorModule  , only : allocationError
    use reportErrorModule  , only : reportError, debugFileOnly_p

    type(StrainCoatType), intent(inout) :: strainCoat
    real(dp)            , intent(in)    :: displ(:), biAxialGate, toMPaScale
    integer             , intent(out)   :: ierr

    !! Local variables
    integer  :: i
    real(dp) :: fppValue

    !! --- Logic section ---

    ierr = 0
    do i = 1, size(strainCoat%tensorRosette%data)

       call calcRosetteStrains (strainCoat%tensorRosette%data(i),displ,ierr)
       if (ierr /= 0) then
          call reportError (debugFileOnly_p,'calcStrainCoatData')
          return
       end if

       !! Keep track of the largest values
       call updateMax (strainCoat%results(i)%epsMax, &
            &          strainCoat%tensorRosette%data(i)%epsP(1))
       call updateMin (strainCoat%results(i)%epsMin, &
            &          strainCoat%tensorRosette%data(i)%epsP(2))
       call updateMax (strainCoat%results(i)%sigMax, &
            &          strainCoat%tensorRosette%data(i)%sigmaP(1))
       call updateMin (strainCoat%results(i)%sigMin, &
            &          strainCoat%tensorRosette%data(i)%sigmaP(2))
       call updateMax (strainCoat%results(i)%gammaMax, &
            &          strainCoat%tensorRosette%data(i)%gammaMax)
       call updateMax (strainCoat%results(i)%tauMax, &
            &          strainCoat%tensorRosette%data(i)%tauMax)
       call updateMax (strainCoat%results(i)%vmeMax, &
            &          strainCoat%tensorRosette%data(i)%epsVM)
       call updateMax (strainCoat%results(i)%vmsMax, &
            &          strainCoat%tensorRosette%data(i)%sigmaVM)
       call updateAngBin (strainCoat%tensorRosette%data(i)%alpha1, &
            &             strainCoat%tensorRosette%data(i)%sigmaP, &
            &             strainCoat%tensorRosette%data(i)%epsP, &
            &             strainCoat%results(i)%angBin)
       if (ierr /= 0) then
          ierr = allocationError('calcStrainCoatData')
          return
       end if

       !! Update the max stress/strain ranges
       strainCoat%results(i)%strRange(1) = strainCoat%results(i)%sigMax &
            &                            - strainCoat%results(i)%sigMin
       strainCoat%results(i)%strRange(2) = strainCoat%results(i)%epsMax &
            &                            - strainCoat%results(i)%epsMin

       !! Accumulate bi-axiality ratio
       if (strainCoat%tensorRosette%data(i)%sigmaP(3) > biAxialGate) then
          call updateBiAxial (strainCoat%results(i)%biAxialSum, &
               &              strainCoat%results(i)%biAxialSqr, &
               &              strainCoat%tensorRosette%data(i)%sigmaP)
          strainCoat%results(i)%nBiAxial = strainCoat%results(i)%nBiAxial + 1
       end if

       select case (strainCoat%fppType)

       case (-1,1) ! Signed max principle stress (MPa)
          fppValue = strainCoat%tensorRosette%data(i)%sigmaP(3)*toMPaScale

       case (-2) ! Signed max principle strain (Micro strain)
          fppValue = strainCoat%tensorRosette%data(i)%epsP(3)*1.0e6_dp

       case default
          cycle

       end select
       if (strainCoat%fppType > 0) then
          strainCoat%results(i)%fatValue = fppValue * strainCoat%results(i)%sCF
       else
          strainCoat%results(i)%fppValue(iBuf) = real(fppValue,sp)
       end if

    end do

  contains

    subroutine updateMax (current,newVal)
      real(dp), intent(inout) :: current
      real(dp), intent(in)    :: newVal
      if (newVal > current) current = newVal
    end subroutine updateMax

    subroutine updateMin (current,newVal)
      real(dp), intent(inout) :: current
      real(dp), intent(in)    :: newVal
      if (newVal < current) current = newVal
    end subroutine updateMin

    subroutine updateBiAxial (current,current2,Pval)
      real(dp), intent(inout) :: current,current2
      real(dp), intent(in)    :: Pval(:)
      real(dp) :: biaxial
      if (abs(Pval(1)) > abs(Pval(2))) then
         biaxial = Pval(2)/Pval(1)
      else
         biaxial = Pval(1)/Pval(2)
      end if
      current  = current  + biaxial
      current2 = current2 + biaxial*biaxial
    end subroutine updateBiAxial

    subroutine updateAngBin (angle,sig,eps,angBin)
      real(dp)         , intent(in)    :: angle ! in the interval [-pi/2,pi/2]
      real(dp)         , intent(in)    :: sig(2), eps(2)
      type(AngleBinPtr), intent(inout) :: angBin(:)
      integer :: iAng, jAng, nBin
      nBin = size(angBin)
      iAng = nint((angle/pi_p+0.5_dp)*nBin) ! bin of largest principal value
      jAng = nint((angle/pi_p+1.0_dp)*nBin) ! bin of smallest principal value
      if (iAng < 1)    iAng = nBin          ! first and last bin are equivalent
      if (jAng > nBin) jAng = jAng - nBin   ! first and last bin are equivalent

      if (.not. associated(angBin(iAng)%p)) then
         allocate(angBin(iAng)%p,STAT=ierr)
         if (ierr /= 0) return
         angBin(iAng)%p%nVal   = 1
         angBin(iAng)%p%sigMax = sig(1)
         angBin(iAng)%p%sigMin = sig(1)
         angBin(iAng)%p%epsMax = eps(1)
         angBin(iAng)%p%epsMin = eps(1)
      else
         angBin(iAng)%p%nVal   = angBin(iAng)%p%nVal + 1
         angBin(iAng)%p%sigMax = max(angBin(iAng)%p%sigMax,sig(1))
         angBin(iAng)%p%sigMin = min(angBin(iAng)%p%sigMin,sig(1))
         angBin(iAng)%p%epsMax = max(angBin(iAng)%p%epsMax,eps(1))
         angBin(iAng)%p%epsMin = min(angBin(iAng)%p%epsMin,eps(1))
      end if

      if (.not. associated(angBin(jAng)%p)) then
         allocate(angBin(jAng)%p,STAT=ierr)
         if (ierr /= 0) return
         angBin(jAng)%p%nVal   = 0
         angBin(jAng)%p%sigMax = sig(2)
         angBin(jAng)%p%sigMin = sig(2)
         angBin(jAng)%p%epsMax = eps(2)
         angBin(jAng)%p%epsMin = eps(2)
      else
         angBin(jAng)%p%sigMax = max(angBin(jAng)%p%sigMax,sig(2))
         angBin(jAng)%p%sigMin = min(angBin(jAng)%p%sigMin,sig(2))
         angBin(jAng)%p%epsMax = max(angBin(jAng)%p%epsMax,eps(2))
         angBin(jAng)%p%epsMin = min(angBin(jAng)%p%epsMin,eps(2))
      end if
    end subroutine updateAngBin

  end subroutine calcStrainCoatData


  subroutine calcAngleData (angBin,popAng,angSpd,sRange,useOldRange)

    !!==========================================================================
    !! Calculate the biggest gap in angle distribution, set angSpd to 180-gap.
    !! Also determine the most popular angle, popAng.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 25 February 2002/2.0
    !!==========================================================================

    type(AngleBinPtr), intent(inout) :: angBin(:)
    real(dp)         , intent(out)   :: popAng, angSpd, sRange(2)
    logical          , intent(in)    :: useOldRange

    !! Local variables
    real(dp) :: binSize
    integer  :: i, iGap, maxGap, firstGap, nBins, mVal

    !! --- Logic section ---

    nBins = size(angBin)
    binSize = 180.0_dp/nBins
    if (.not. useOldRange) sRange = 0.0_dp

    !! Find the stress and strain ranges, and the most popular angle
    iGap = 0
    mVal = 0
    do i = 1, nBins
       if (associated(angBin(i)%p)) then
          if (.not. useOldRange) then
             sRange(1) = max(sRange(1),angBin(i)%p%sigMax-angBin(i)%p%sigMin)
             sRange(2) = max(sRange(2),angBin(i)%p%epsMax-angBin(i)%p%epsMin)
          end if
          if (angBin(i)%p%nVal > mVal) then
             iGap = i
             mVal = angBin(i)%p%nVal
          else if (angBin(i)%p%nVal == 0) then
             deallocate(angBin(i)%p)
             nullify(angBin(i)%p)
          end if
       end if
    end do
    popAng = iGap*binSize - 90.0_dp

    !! Find the angle spread
    firstGap = 0
    maxGap = 0
    iGap = 0
    do i = 1, nBins
       if (associated(angBin(i)%p)) then
          if (firstGap == 0) then
             firstGap = i
          else if (iGap > 0) then
             maxGap = max(maxGap,i-iGap+1)
             iGap = 0
          end if
       else if (iGap == 0 .and. firstGap > 0) then
          iGap = i
       end if
    end do

    if (iGap > 0) firstGap = firstGap + nBins-iGap+1
    if (firstGap > maxGap) maxGap = firstGap

    angSpd = 180.0_dp - maxGap*binSize

  end subroutine calcAngleData


  subroutine printStrainCoatInput (strainCoats,lpu)

    !!==========================================================================
    !! Print strain coat input data to ASCII results file.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 19 June 2008/1.0
    !!==========================================================================

    type(StrainCoatType), intent(in) :: strainCoats(:)
    integer             , intent(in) :: lpu

    !! Local variables
    integer           :: i, j
    character(len=80) :: cfmt

    !! --- Logic section ---

    write(lpu,600)
    do i = 1, size(strainCoats)
       cfmt = "(2I8,A10,1PE14.5,I10,I11,I3,0PF12.3)"
       if (size(strainCoats(i)%results) > 1) then
          cfmt(36:) = "/(14X,A12,1PE14.5,I10,I11,I3,0PF12.3))"
       end if
       write(lpu,cfmt)   strainCoats(i)%tensorRosette%id%baseId, &
            &            strainCoats(i)%tensorRosette%id%userId, &
            (resultSet_p(strainCoats(i)%results(j)%resultSet), &
            &            strainCoats(i)%tensorRosette%data(j)%zPos, &
            &            strainCoats(i)%results(j)%matGroup, &
            &            strainCoats(i)%results(j)%SNcurve, &
            &            strainCoats(i)%results(j)%sCF, &
            &   j=1,size(strainCoats(i)%results))
    end do
    write(lpu,620)

600 format(//5X,'STRAIN COAT PROPERTY SUMMARY' /5X,74('-') &
          & /5X,'Strain Coat  Result set  Z-position  Material Group  ', &
          &     'S-N curve  SCF factor')
620 format(4X,75('-')/)

  end subroutine printStrainCoatInput


  subroutine printStrainCoatData (strainCoats,lpu)

    !!==========================================================================
    !! Print strain coat results data to ASCII results file.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 5 June 2001/1.0
    !!==========================================================================

    type(StrainCoatType), intent(in) :: strainCoats(:)
    integer             , intent(in) :: lpu

    !! Local variables
    integer            :: i, j
    character(len=112) :: cfmt

    !! --- Logic section ---

    write(lpu,600)
    do i = 1, size(strainCoats)
       cfmt = "(2I8,A10,3X,1P,2(1X,3E13.5),1X,2E13.5,1X,0P,2F13.5)"
       if (size(strainCoats(i)%results) > 1) then
          cfmt(51:) = "/(14X,A12,3X,1P,2(1X,3E13.5),1X,2E13.5,1X,0P,2F13.5))"
       end if
       write(lpu,cfmt)   strainCoats(i)%tensorRosette%id%baseId, &
            &            strainCoats(i)%tensorRosette%id%userId, &
            (resultSet_p(strainCoats(i)%results(j)%resultSet), &
            &            strainCoats(i)%results(j)%sigMax, &
            &            strainCoats(i)%results(j)%tauMax, &
            &            strainCoats(i)%results(j)%vmsMax, &
            &            strainCoats(i)%results(j)%epsMax, &
            &            strainCoats(i)%results(j)%gammaMax, &
            &            strainCoats(i)%results(j)%vmeMax, &
            &   BiAxMean(strainCoats(i)%results(j)), &
            & BiAxStdDev(strainCoats(i)%results(j)), &
            &            strainCoats(i)%results(j)%popAngle, &
            &            strainCoats(i)%results(j)%angSpread, &
            &   j=1,size(strainCoats(i)%results))
       if (strainCoats(i)%fppType > 0) then
          write(lpu,612) (strainCoats(i)%results(j)%damage, &
               & j=1,size(strainCoats(i)%results))
       end if
    end do
    write(lpu,620)

600 format(//5X,'STRAIN COAT RECOVERY SUMMARY' /5X,24('-'), &
         &      '   --------------- Stress --------------', &
         &      '   --------------- Strain --------------', &
         &      '   ------ Biaxiality ------   - Principal dir. angle -' &
         &  /5X,'Strain Coat  Result set  ', &
         &      '    Max P1      Max shear  Max von Mises', &
         &      '    Max P1      Max shear  Max von Mises ', &
         &      '    Mean      Std. Dev.      Most popular    spread' )
612 format(20X,'Damage =  ',1P3E13.5)
620 format(4X,159('-')/)

  end subroutine printStrainCoatData


  function BiAxMean (coat)

    !!==========================================================================
    !! Returns biaxiality mean for the given strain coat element.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 1 March 2002/1.0
    !!==========================================================================

    type(CoatResultType), intent(in) :: coat
    real(dp)                         :: BiAxMean

    !! --- Logic section ---

    BiAxMean = coat%biAxialSum / max(1,coat%nBiAxial)

  end function BiAxMean


  function BiAxStdDev (coat)

    !!==========================================================================
    !! Returns biaxiality standard deviation for the given strain coat element.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 1 March 2002/1.0
    !!==========================================================================

    type(CoatResultType), intent(in) :: coat
    real(dp)                         :: BiAxStdDev, BiAxMean, dvar, dnum

    !! --- Logic section ---

    BiAxStdDev = 0.0_dp
    dnum = coat%nBiAxial
    if (dnum > 1.0_dp) then
       BiAxMean = coat%biAxialSum/dnum
       dvar     = coat%biAxialSqr/dnum - BiAxMean*BiAxMean
       if (dvar > 0.0_dp) then
          BiAxStdDev = sqrt(dvar*dnum/(dnum-1.0_dp))
       end if
    end if

  end function BiAxStdDev

end module StrainCoatModule
