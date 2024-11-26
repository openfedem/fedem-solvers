!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file initiateFunctionTypeModule.f90
!> @brief Initialization of functions from the solver input file.

!!==============================================================================
!> @brief Module with a namelist for reading function data.
module FunctionNamelistModule

  use KindModule  , only : dp, lfnam_p
  use IdTypeModule, only : ldesc_p

  implicit none

  !! Define the FUNCTION namelist
  integer, parameter :: maxRealData_p = 10000   !< Max number of data values
  integer            :: id                      !< Base ID of the function
  integer            :: extId(10)               !< User ID path of the function
  integer            :: realDataSize            !< Size of real data array
  integer            :: extrapolationType       !< Extrapolation type flag
  integer            :: channel                 !< Column index for data file
  integer            :: nArg                    !< Number of function arguments
  integer            :: waveId                  !< Base Id for wave function
  integer            :: waveDir                 !< Wave direction angle [deg]
  integer            :: seed                    !< Random seed for wave function
  real(dp)           :: realData(maxRealData_p) !< Real data array
  character(ldesc_p) :: extDescr                !< User description
  character(len=40)  :: type                    !< Function type flag
  character(lfnam_p) :: fileName                !< Name of function data file
  character(len=500) :: expression              !< Explicit function expression

  namelist /FUNCTION/ id, extId, extDescr, type, realDataSize, realData, &
       &              extrapolationType, fileName, channel, &
       &              expression, nArg, waveId, waveDir, seed

contains

  !> @cond FULL_DOC
  !> @brief Reads the FUNCTION namelist from the given file unit.
  subroutine read_FUNCTION (infp,stat)
    integer, intent(in)  :: infp
    integer, intent(out) :: stat

    !! Default values
    id=0; extId=0; extDescr=''
    type=''; realDataSize=0; extrapolationType=0
    fileName=''; channel=0; expression=''; nArg=1
    waveId=0; waveDir=0; seed=0

    read(infp,nml=FUNCTION,iostat=stat)

  end subroutine read_FUNCTION
  !> @endcond

end module FunctionNamelistModule


!!==============================================================================
!> @brief Initialization of function shapes from the solver input file.
module initiateFunctionTypeModule

  use kindModule, only : lfnam_p

  implicit none

  private

  !> @cond NO_DOCUMENTATION
  integer, allocatable, save :: argSensor(:)
  integer,              save :: maxChannel = 0

  type IntVec
     integer, pointer :: values(:)
  end type IntVec

  type(IntVec), allocatable, save :: argSensors(:)
  !> @endcond

  !> @brief Data type for temporary linked list of RAO data.
  type RAOType
     character(len=lfnam_p) :: RAOfile     !< Name of RAO data file
     integer                :: waveDir     !< Wave direction angle [deg]
     integer                :: waveId      !< Base Id for wave function
     integer                :: motionId(6) !< Base Ids for resulting motions
     type(RAOType), pointer :: next        !< Pointer to next RAO
  end type RAOType

  public :: InitiateFunctions, InitiateEngines1, InitiateEngines2


contains

  !!============================================================================
  !> @brief Initializes the function type with data from the solver input file.
  !>
  !> @param[in] infp File unit number for the solver input file
  !> @param[in] gravity Gravitation vector
  !> @param[in] seaDepth Water depth for wave function
  !> @param[in] waveFunc Pointer to the wave function to be used
  !> @param[out] functions Array of all function shapes in the model
  !> @param[out] err Error flag
  !>
  !> @details This subroutine also processes the irregular wave function,
  !> generating time history parameters from the wave spectrum, etc.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 7 Sep 1999

  subroutine InitiateFunctions (infp,gravity,seaDepth,waveFunc,functions,err)

    use FunctionTypeModule       , only : FunctionType, NullifyFunction
    use FunctionTypeModule       , only : DeallocateFunctions
    use explicitFunctionsModule  , only : spline_p, nSpln_p, funcType_p
    use explicitFunctionsModule  , only : WAVE_SINUS_p, WAVE_EMBEDDED_p
    use explicitFunctionsModule  , only : WAVE_STOKES5_p, WAVE_STREAMLINE_p
    use explicitFunctionsModule  , only : SINUSOIDAL_p, CONSTANT_p
    use explicitFunctionsModule  , only : DEVICE_FUNCTION_p, MATH_EXPRESSION_p
    use explicitFunctionsModule  , only : USER_DEFINED_p, dbgFunc
    use waveFunctionsModule      , only : initFUNC4, initFUNC7, initFUNC8
    use waveFunctionsModule      , only : initFUNC9, waveNumber, checkDepth
    use kindModule               , only : dp, pi_p, epsDiv0_p
    use IdTypeModule             , only : initId, getId, ReportInputError
    use manipMatrixModule        , only : MatrixPointToArray_real
    use inputUtilities           , only : iuGetNumberOfEntries
    use inputUtilities           , only : iuSetPosAtNextEntry, iuCharToInt
    use progressModule           , only : lterm
    use reportErrorModule        , only : error_p, warning_p, note_p
    use reportErrorModule        , only : debugFileOnly_p, reportError
    use reportErrorModule        , only : AllocationError, getErrorFile
    use fileUtilitiesModule      , only : getDBGfile
    use FiDeviceFunctionInterface, only : FiDF_Open
    use FFaFilePathInterface     , only : ffa_checkPath
    use FFaCmdLineArgInterface   , only : ffa_cmdlinearg_getbool
    use FFaCmdLineArgInterface   , only : ffa_cmdlinearg_getstring
    use FFaCmdLineArgInterface   , only : ffa_cmdlinearg_intValue
    use FFaMathExprInterface     , only : ffame_create
    use FFaUserFuncInterface     , only : ffauf_init, ffauf_getnopar
    use FFaUserFuncInterface     , only : ffauf_getflag
    use FunctionNameListModule ! Defines all variables in the FUNCTION namelist

    integer           , intent(in)    :: infp
    real(dp)          , intent(in)    :: gravity(3)
    real(dp)          , intent(inout) :: seaDepth
    type(FunctionType), pointer       :: waveFunc, functions(:)
    integer           , intent(out)   :: err

    !! Local variables
    integer, parameter :: max_embsw_p = 20
    integer            :: i, j, k, nFunc, nWave, noSpectrum, fileId, stat, lpu
    integer            :: printFunc, intDataSize, intData(4+max_embsw_p)
    real(dp)           :: H3, T1, Ga, W0, dW, g, sfunc(3,max_embsw_p), ramp(2)
    real(dp), pointer  :: waveData(:,:), tmpW(:)
    character(len=128) :: err_msg
    character(len=520) :: plugin

    type(RAOType), pointer :: RAOlist, tmp

    !! --- Logic section ---

    nFunc = iuGetNumberOfEntries(infp,'&FUNCTION',err)
    if (err /= 0) return
    write(lterm,*) 'Number of &FUNCTION =',nFunc

    call DeallocateFunctions (functions)
    allocate(functions(nFunc),STAT=stat)
    if (stat /= 0) then
       err = AllocationError('InitiateFunctions 1')
       return
    end if

    printFunc = ffa_cmdlinearg_intValue('printFunc')
    if (nFunc > 0) then
       if (printFunc == 2) then
          dbgFunc = getDBGfile(12,'functions.txt',.true.)
#ifdef FT_DEBUG
       else
          dbgFunc = getDBGfile(12,'functions.dbg')
#endif
       end if
    end if

    call ffa_cmdlinearg_getstring ('plugin',plugin)

    !! Gravitation constant
    g = sqrt(sum(gravity*gravity))

    nullify(RAOlist)

    if (associated(waveFunc) .and. ffa_cmdlinearg_intValue('FNV') > 0) then
       !! Do not generate JONSWAP spectrum for this wave function (FNV is used)
       noSpectrum = waveFunc%id%baseId
       printFunc = 0
    else
       noSpectrum = 0
    end if

    do i = 1, nFunc
       nullify(waveData)
       nullify(tmpW)

       call NullifyFunction (functions(i))
       if (.not. iuSetPosAtNextEntry(infp,'&FUNCTION')) then
          err = err - 1
          call ReportInputError ('FUNCTION',i)
          cycle
       end if

       call read_FUNCTION (infp,stat)
       if (stat /= 0) then
          err = err - 1
          call ReportInputError ('FUNCTION',i)
          cycle
       end if

       call initId (functions(i)%id,id,extId,extDescr,stat)

       stat = err
       functions(i)%type = iuCharToInt(type,funcType_p,err)
       if (err < stat) call ReportInputError ('FUNCTION',i,functions(i)%id)

       stat = 0
       nWave = 0
       intDataSize = 3
       if (functions(i)%type == DEVICE_FUNCTION_p) then
          intDataSize = 5
          if (realDataSize < 3) then
             realData(1+realDataSize:3) = 0.0_dp
             if (realDataSize < 2) realData(2) = 1.0_dp
             realDataSize = 3
          end if
       else if (functions(i)%type == WAVE_SINUS_p) then
          if (id == noSpectrum) then
             !! FNV method is used, do not generate JONSWAP spectrum here
             nWave = 0
             realDataSize = 5
          else
             !! Wave function composed of nWave sinusoidal components
             nWave = realDataSize/3
          end if
          if (waveId == 0 .and. nWave > 0) then
             if (seaDepth > 0.0_dp .and. g > epsDiv0_p) then
                !! Precompute wave numbers for each component from depth and g
                waveData => MatrixPointToArray_real(realData,4,nWave)
                realDataSize = 4*nWave
             else
                !! The wave number is taken as omega^2/g for all components
                waveData => MatrixPointToArray_real(realData,3,nWave)
             end if
          end if
          intDataSize = 4
       else if (functions(i)%type == WAVE_EMBEDDED_p) then
          !! Wave function composed of nWave sinusoidal components
          !! + embedded streamline waves
          nWave = realDataSize/3
          waveData => MatrixPointToArray_real(realData,4,nWave)
          realDataSize = 4*nWave
          intDataSize = 4+channel
       else if (functions(i)%type == USER_DEFINED_p) then
          intDataSize = 4 ! Add space for evaluation flag
       end if
       if (nWave > 0 .and. realDataSize > maxRealData_p) then
          allocate(tmpW(realDataSize),STAT=stat)
          if (stat /= 0) then
             err = AllocationError('InitiateFunctions 2')
             return
          else
             waveData => MatrixPointToArray_real(tmpW,realDataSize/nWave,nWave)
          end if
       else if (realDataSize > maxRealData_p .or. realDataSize < 0) then
          err = err - 1
          write(err_msg,600) realDataSize, maxRealData_p
600       format('Number of parameters =',I6,'  Max allowed =',I6)
          call reportError (error_p,'Too many real parameters in Function'// &
               &            getId(functions(i)%id),err_msg)
          cycle
       else if (intDataSize > size(intData)) then
          err = err - 1
          write(err_msg,600) intDataSize, size(intData)
          call reportError (error_p,'Too many integer parameters in Function'//&
               &            getId(functions(i)%id),err_msg)
          cycle
       end if

       !! ===================================
       !! Preprocessing of wave function data
       !! ===================================

       if (functions(i)%type > 0 .and. functions(i)%type < WAVE_SINUS_p) then

          if (nArg > 1 .and. g > epsDiv0_p) then
             !! Sinusoidal functions with wave number need the gravity constant
             realDataSize = realDataSize + 1
             realData(realDataSize) = g

             if (seaDepth > 0.0_dp) then
                !! Check that the finite water depth is not too large
                if (functions(i)%type == SINUSOIDAL_p) then
                   W0 = realData(1)*2.0_dp*pi_p
                else
                   W0 = min(realData(1),realData(2))*2.0_dp*pi_p
                end if
                if (checkDepth(seaDepth,waveNumber(W0,g,seaDepth))) then
                   write(err_msg,610) seaDepth,functions(i)%type
610                format('The specified water depth',1PE12.5, &
                        ' is too high for finite depth wave calculation.',I6)
                   call reportError (note_p,err_msg, &
                        'Changing to infinite water depth formulation.')
                   seaDepth = 0.0_dp
                end if
             end if
          end if

       else if (functions(i)%type == WAVE_SINUS_p) then

          if (fileName /= '') then
             !! Wave function from file
             call ffa_checkPath (fileName)
             if (waveId > 0) then
                !! RAO function, store data temporarily in a linked list
                if (.not.insertRAO(fileName,channel,waveDir,waveId,id)) return
             else if (associated(waveData)) then
                !! Read wave component data from ASCII file
                call initFUNC4 (fileName,g,seaDepth,waveData,stat,seed)
             end if
             waveDir = 1 ! Number of wave directions (always 1)
          else if (associated(waveData)) then
             !! Generate wave spectrum from statistical parameters
             H3 = realData(1) ! Significant wave height
             T1 = realData(2) ! Mean/Peak wave period
             W0 = realData(3) ! Lowest circular frequency in spectrum
             dW = realData(4) ! Frequency increment betweeen wave components
             Ga = realData(5) ! Spectral peakedness
             if (extrapolationType > 4) then
                !! Use new implementation, based on the DNV RP 205
                j = mod(extrapolationType,10) - 4
                if (j == 1 .or. j == 2) Ga = 1.0_dp ! Pierson-Moskowitz
                if (j == 1 .or. j == 3) W0 = -W0 ! Constant frequency intervals
                call initFUNC4 (H3,T1,Ga,W0,dW,g,seaDepth,waveData,stat, &
                     &          seed,waveDir,extrapolationType/10)
                if (stat > 0 .and. stat < 4) then
                   call reportError (warning_p,'Check definition of Wave '// &
                        &            'Function'//getId(functions(i)%id))
                end if
             else if (extrapolationType > 0) then
                !! Old implementation, based on O. Faltinsen's book
                call initFUNC4 (extrapolationType,H3,T1,Ga,W0,dW,g,seaDepth, &
                     &          waveData,stat,seed)
             else
                stat = -1
             end if
          end if

          if (seaDepth > 0.0_dp .and. stat > 3) then
             write(err_msg,610) seaDepth,stat-3
             call reportError (note_p,err_msg, &
                  'Changing to infinite water depth formulation.')
             seaDepth = 0.0_dp
          end if

       else if (functions(i)%type == WAVE_STOKES5_p .or. &
            &   functions(i)%type == WAVE_STREAMLINE_p) then

          !! Solve the nonlinear streamline wave equation
          T1 = realData(1) ! Wave period
          H3 = realData(2) ! Wave height
          W0 = realData(3) ! Current velocity
          if (functions(i)%type == WAVE_STOKES5_p) then
             call initFUNC7 (H3,T1,W0,g,seaDepth,realData,stat)
          else
             call initFUNC8 (H3,T1,W0,g,seaDepth,realData,stat)
          end if
          realDataSize = 2*stat+1
          if (stat > 0 .and. printFunc == 1) then
             lpu = getErrorFile()
             write(lpu,620) trim(getId(functions(i)%id))
             write(lpu,"(I8,1P2E13.5)") (j,realData(3+j),realData(2+stat+j), &
                  &                      j=1,stat-1)
             write(lpu,"('')")
620          format(/5X,'>>> Wave parameters for Function',A,' <<<', &
                  & /5X,'  i     amp_i        c_(i+1)')
          end if

       else if (functions(i)%type == WAVE_EMBEDDED_p) then

          !! Embedded streamline waves
          H3 = realData(1) ! Significant wave height
          T1 = realData(2) ! Mean wave period
          W0 = realData(3) ! Lowest circular frequency in spectrum
          dW = realData(4) ! Frequency increment betweeen wave components
          Ga = realData(5) ! Spectral peakedness
          ramp = realData(6:7) ! Transition parameters
          sfunc(:,1:channel) = reshape(realData(8:7+3*channel),(/3,channel/))
          call initFUNC9 (extrapolationType,4,nWave,H3,T1,Ga,W0,dw,g,seaDepth, &
               &          ramp,sfunc(:,1:channel),intData,realData,stat,seed)
          if (stat > maxRealData_p) then
             stat = -stat
             write(err_msg,600) -stat, maxRealData_p
             call reportError (error_p,'Too many wave components. '//err_msg)
          else if (stat > 0) then
             realDataSize = stat
          end if

       end if
       if (stat < 0) then
          err = err - 1
          call reportError (error_p,'Failed to initialize Wave Function'// &
               &            getId(functions(i)%id),ierr=stat)
          cycle
       else if (associated(waveData) .and. printFunc == 1) then
          lpu = getErrorFile()
          nWave = nWave/waveDir
          do k = 1, waveDir
             write(lpu,640,advance='NO') trim(getId(functions(i)%id))
             if (waveDir > 1) then
                write(lpu,641) 90.0_dp*real(k+k-waveDir-1,dp)/real(waveDir+1,dp)
             else
                write(lpu,"(' <<<')")
             end if
             write(lpu,643)
             do j = 1, nWave
                write(lpu,"(I8,1P4E13.5)") j,waveData(:,j+nWave*(k-1))
             end do
          end do
          write(lpu,"('')")
640       format(/5X,'>>> Wave spectrum for Function',A)
641       format(', theta =',F8.2,' <<<')
643       format( 5X,'  i      A_i       omega_i      epsilon_i      k_i')
       end if

       allocate(functions(i)%intParameters(intDataSize),STAT=stat)
       if (stat /= 0) then
          err = AllocationError('InitiateFunctions 3')
          return
       end if
       if (associated(tmpW)) then
          functions(i)%realParameters => tmpW
       else
          allocate(functions(i)%realParameters(realDataSize),STAT=stat)
          if (stat /= 0) then
             err = AllocationError('InitiateFunctions 4')
             return
          else if (realDataSize > 0) then
             functions(i)%realParameters = realData(1:realDataSize)
          end if
       end if

       functions(i)%intParameters(1) = functions(i)%type
       functions(i)%intParameters(2) = extrapolationType
       functions(i)%intParameters(3) = 0

       select case (functions(i)%type)
       case (SPLINE_p)
          !! Store index to static spline data of this function
          functions(i)%intParameters(3) = splineIndex()

       case (SINUSOIDAL_p)
          !! Flag that this is a general sine function and not a wave function.
          !! See explicitFunctionsModule::FUNC1 for the interpretation.
          functions(i)%intParameters(3) = 1

       case (WAVE_SINUS_p)
          if (nWave > 0) then
            !! Store the number of wave components of this function
            functions(i)%intParameters(3) = nWave
         else
            !! Store the random seed (for use in FNV calculations)
            functions(i)%intParameters(3) = seed
         end if
         !! Store the number of wave directions of this function
         functions(i)%intParameters(4) = waveDir

       case (WAVE_EMBEDDED_p)
          !! Store data for embedded streamline waves for this function
          functions(i)%intParameters(3:) = intData(3:intDataSize)

       case (CONSTANT_p)
          !! No arguments needed for this function
          functions(i)%nArg = 0

       case (DEVICE_FUNCTION_p)

          !! Device function (DAC/ASCII-file or RPC-channel)
          fileId = channel
          call ffa_checkPath (fileName)
          call FiDF_Open (trim(fileName),fileId,stat)
          if (stat /= 0) then
             err = err - 1
             call reportError (error_p, &
                  'Could not open external device file: '//fileName,ierr=stat)
          end if
          functions(i)%intParameters(3) = fileId
          functions(i)%intParameters(4) = channel
          functions(i)%intParameters(5) = int(realData(3))
          if (fileId == 0) functions(i)%nArg = 0

       case (MATH_EXPRESSION_p)

          !! User-defined math expression
          call ffame_create (nArg,expression,id,stat)
          if (stat < 0) then
             err = err - 1
             call reportError (error_p, &
                  'Malformed function expression: '//expression,ierr=stat)
          end if
          functions(i)%intParameters(3) = id
          functions(i)%nArg = nArg

       case (USER_DEFINED_p)

          !! User-defined function (plug-in)
          stat = ffauf_init(plugin,channel,expression)
          if (stat < 0) then
             err = err - 1
             call reportError (error_p, &
                  'Failed to initialize user-defined function',ierr=stat)
          else if (expression /= '') then
             call reportError (note_p,expression)
          end if
          functions(i)%intParameters(3) = channel
          functions(i)%intParameters(4) = ffauf_getflag(channel)
          functions(i)%nArg = max(0,stat)

          stat = ffauf_getnopar(channel)
          if (realDataSize < stat) then
             err = err - 1
             write(err_msg,660) realDataSize, stat
660          format('Too few function parameters',I4,'. Expected',I4,'.')
             call reportError (error_p, &
                  'Failed to initialize user-defined function.',err_msg)
          end if

       end select

    end do

    !! Now create motions functions from the specified RAO tables
    do while (associated(RAOlist) .and. err >= 0)
       call createRAOMotions (RAOlist,functions,err)
       tmp => RAOlist
       RAOlist => tmp%next
       deallocate(tmp)
    end do

    if (err < 0) call reportError (debugFileOnly_p,'InitiateFunctions')

  contains

    !> @brief Returns the next spline index in the cyclic range [1,nSpln_p].
    function SplineIndex ()
      integer, save :: index = 1
      integer       :: SplineIndex
      SplineIndex = index
      if (index < nSpln_p) then
         index = index + 1
      else
         index = 1
      end if
    end function SplineIndex

    !> @brief Inserts data for an RAO function into a linked list.
    function insertRAO (RAOfile,RAOdof,waveDir,waveId,motionId)
      character(len=*), intent(in) :: RAOfile
      integer         , intent(in) :: RAOdof, waveDir, waveId, motionId
      type(RAOType)   , pointer    :: curr, prev
      logical                      :: insertRAO
      insertRAO = .true.
      if (RAOdof < 1 .or. RAOdof > 6) then
         call ReportInputError ('FUNCTION',i,functions(i)%id)
         err = err - 1
         return
      end if
      nullify(prev)
      curr => RAOlist
      do while (associated(curr))
         if ( curr%RAOfile == RAOfile .and. &
              curr%waveDir == waveDir .and. &
              curr%waveId  == waveId ) goto 100
         prev => curr
         curr => prev%next
      end do
      allocate(curr,STAT=stat)
      if (stat /= 0) then
         err = AllocationError('InitiateFunctions 3')
         insertRAO = .false.
         return
      else if (associated(prev)) then
         prev%next => curr
      else
         RAOlist => curr
      end if
      nullify(curr%next)
      curr%RAOfile = fileName
      curr%waveDir = waveDir
      curr%waveId = waveId
      curr%motionId = 0
100   curr%motionId(RAOdof) = motionId
    end function insertRAO

  end subroutine InitiateFunctions


  !!============================================================================
  !> @brief Creates the real data for the RAO motion functions.
  !>
  !> @param[in] RAO Data for an RAO motion function
  !> @param functions Array of all function shapes in the model
  !> @param err Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 4 Mar 2008

  subroutine createRAOMotions (RAO,functions,err)

    use kindModule             , only : dp, pi_p
    use FunctionTypeModule     , only : FunctionType, GetPtrToId
    use explicitFunctionsModule, only : WAVE_SINUS_p, SINUSOIDAL_p
    use explicitFunctionsModule, only : COMPL_SINUS_p, DELAYED_COMPL_SINUS_p
    use reportErrorModule      , only : reportError, error_p, debugFileOnly_p
    use reportErrorModule      , only : AllocationError
    use FiRAOTableInterface    , only : FiConvertWaveData
    use FiRAOTableInterface    , only : FiExtractMotion, FiReleaseMotion

    type(RAOType)     , intent(in)    :: RAO
    type(FunctionType), intent(inout) :: functions(:)
    integer           , intent(inout) :: err

    !! Local variables
    integer                     :: idof, ierr, nRows, nComp
    real(dp)                    :: sFunc(6)
    character(len=64)           :: msg
    type(FunctionType), pointer :: waveFunc, motionFunc

    !! --- Logic section ---

    ierr = err
    waveFunc => GetPtrToId(functions,RAO%waveId)
    if (.not. associated(waveFunc)) then
       err = err - 1
       call reportError (error_p,'Non-existing wave function referred')
    else if (waveFunc%type == WAVE_SINUS_p) then
       nComp = waveFunc%intParameters(3)
       nRows = size(waveFunc%realParameters)/nComp
       call FiConvertWaveData (trim(RAO%RAOfile),RAO%waveDir, &
            &                  nRows,nComp,waveFunc%realParameters(1),err)
    else if (waveFunc%type == SINUSOIDAL_p) then
       nComp = 1
       nRows = 3
       sFunc(1) =  waveFunc%realParameters(3)
       sFunc(2) =  waveFunc%realParameters(1) * 2.0_dp*pi_p
       sFunc(3) = -waveFunc%realParameters(2) * 2.0_dp*pi_p
       call FiConvertWaveData (trim(RAO%RAOfile),RAO%waveDir, &
            &                  nRows,nComp,sFunc(1),err)
    else if (waveFunc%type == COMPL_SINUS_p .or. &
         &   waveFunc%type == DELAYED_COMPL_SINUS_p) then
       nComp = 2
       nRows = 3
       sFunc(1) =  waveFunc%realParameters(5)
       sFunc(2) =  waveFunc%realParameters(1) * 2.0_dp*pi_p
       sFunc(3) = -waveFunc%realParameters(3) * 2.0_dp*pi_p
       sFunc(4) =  waveFunc%realParameters(6)
       sFunc(5) =  waveFunc%realParameters(2) * 2.0_dp*pi_p
       sFunc(6) = -waveFunc%realParameters(4) * 2.0_dp*pi_p
       call FiConvertWaveData (trim(RAO%RAOfile),RAO%waveDir, &
            &                  nRows,nComp,sFunc(1),err)
    else
       err = err - 1
       call reportError (error_p,'Invalid Wave function type')
    end if
    if (err == ierr-2) then
       write(msg,"('Failed to read data for wave angle',I4)") RAO%waveDir
       call reportError (error_p,trim(msg)//' from RAO-file '//RAO%RAOfile)
    end if

    do idof = 1, 6
       if (RAO%motionId(idof) > 0 .and. err == ierr) then
          motionFunc => GetPtrToId(functions,RAO%motionId(idof))
          if (.not. associated(motionFunc)) then
             err = err - 1
             call reportError (error_p,'Non-existing motion function referred')
          else if (motionFunc%type /= WAVE_SINUS_p) then
             err = err - 1
             call reportError (error_p,'Invalid RAO motion function type')
          else
             motionFunc%intParameters(3) = nComp
             if (size(motionFunc%realParameters) /= 3*nComp) then
                deallocate(motionFunc%realParameters)
                allocate(motionFunc%realParameters(3*nComp),STAT=ierr)
                if (ierr /= 0) then
                   err = err + AllocationError('createRAOMotions')
                   return
                end if
             end if
             call FiExtractMotion (idof-1,motionFunc%realParameters(1),err)
          end if
       end if
    end do

    call FiReleaseMotion ()

    if (err < ierr) call reportError (debugFileOnly_p,'createRAOMotions')

  end subroutine createRAOMotions


  !!============================================================================
  !> @brief Initializes the engine type with data from the solver input file.
  !>
  !> @param[in] infp File unit number for the solver input file
  !> @param[out] engines Array of all general functions in the model
  !> @param[in] functions Array of all function shapes in the model
  !> @param[out] err Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 7 Sep 1999

  subroutine InitiateEngines1 (infp,engines,functions,err)

    use FunctionTypeModule    , only : FunctionType, EngineType, GetPtrToId
    use FunctionTypeModule    , only : DeallocateEngines, NullifyEngine
    use IdTypeModule          , only : ldesc_p, initId, ReportInputError
    use inputUtilities        , only : iuGetNumberOfEntries, iuSetPosAtNextEntry
    use progressModule        , only : lterm
    use explicitFunctionsModule,only : USER_DEFINED_p, DEVICE_FUNCTION_p
    use reportErrorModule     , only : reportError, debugFileOnly_p
    use reportErrorModule     , only : AllocationError
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_isTrue

    integer           , intent(in)  :: infp
    type(EngineType)  , pointer     :: engines(:)
    type(FunctionType), intent(in)  :: functions(:)
    integer           , intent(out) :: err

    !! Local variables
    integer :: i, j, nEng, stat

    !! Define the ENGINE namelist
    character(ldesc_p) :: extDescr, tag
    integer            :: id, extId(10), functionId, nArg, argSensorId(10)
    namelist /ENGINE/ id, extId, extDescr, tag, functionId, nArg, argSensorId

    !! --- Logic section ---

    nEng = iuGetNumberOfEntries(infp,'&ENGINE',err)
    if (err /= 0) return
    write(lterm,*) 'Number of &ENGINE =', nEng

    call DeallocateEngines (engines)
    allocate(engines(nEng),argSensor(nEng),argSensors(nEng),STAT=stat)
    if (stat /= 0) then
       err = AllocationError('InitateEngines')
       return
    end if

    argSensor = 0

    do i = 1, nEng

       nullify(argSensors(i)%values)
       call NullifyEngine (engines(i))
       if (.not. iuSetPosAtNextEntry(infp,'&ENGINE')) then
          err = err - 1
          call ReportInputError ('ENGINE',i)
          cycle
       end if

       !! Default values
       id=0; extId=0; extDescr=''; tag=''
       functionId=0; nArg=0; argSensorId=0

       read(infp,nml=ENGINE,iostat=stat)
       if (stat /= 0) then
          err = err - 1
          call ReportInputError ('ENGINE',i)
          cycle
       end if

       call initId (engines(i)%id,id,extId,extDescr,stat)
       engines(i)%tag = tag

       if (functionId > 0) then
          engines(i)%func => GetPtrToId(functions,functionId)
          if (.not. associated(engines(i)%func)) then
             err = err - 1
             call ReportInputError ('ENGINE',i,engines(i)%id)
             cycle
          else if (nArg < engines(i)%func%nArg) then
             err = err - 1
             write(extDescr,"('nArg =',I2,', expected at least',I2)") &
                  &             nArg, engines(i)%func%nArg
             call ReportInputError ('ENGINE',i,engines(i)%id, &
                  &                 'Too few function arguments '//extDescr)
          end if
       end if

       !! Temporarily handle the sensor arguments.
       !! Actually connecting the pointers is done in InitiateEngines2
       !! after all engines have been read and initialized.
       if (nArg > 1) then
          allocate(argSensors(i)%values(nArg),STAT=stat)
          if (stat /= 0) then
             err = AllocationError('InitateEngines1')
             return
          end if
          do j = 1, nArg
             argSensors(i)%values(j) = argSensorId(j)
          end do
       else if (nArg == 1) then
          argSensor(i) = argSensorId(1)
       end if

       if (ffa_cmdlinearg_isTrue('allEngineVars')) then
          engines(i)%saveVar = 2 ! Evaluate this engine at each saved time step
       end if
       if (associated(engines(i)%func)) then
          if (engines(i)%func%type == USER_DEFINED_p) then
             if (mod(engines(i)%func%intParameters(4),4) > 1) then
                engines(i)%saveVar = 1 ! Evaluate engine after each time step
             end if
          else if (engines(i)%func%type == DEVICE_FUNCTION_p) then
             if (engines(i)%func%nArg == 0) then
                maxChannel = max(maxChannel,engines(i)%func%intParameters(4))
             end if
          end if
       end if

    end do

    if (err < 0) call reportError (debugFileOnly_p,'InitiateEngines1')

  end subroutine InitiateEngines1


  !!============================================================================
  !> @brief Connects the engines to their argument sensors.
  !>
  !> @param engines Array of all general functions in the model
  !> @param[in] sensors All sensors (function argument objects) in the model
  !> @param[out] err Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 20 Jun 2002

  subroutine InitiateEngines2 (engines,sensors,err)

    use FunctionTypeModule       , only : EngineType, GetPtrToId
    use SensorTypeModule         , only : SensorType, GetPtrToId
    use IdTypeModule             , only : ReportInputError
    use explicitFunctionsModule  , only : DEVICE_FUNCTION_p
    use progressModule           , only : lterm
    use reportErrorModule        , only : allocationError
    use reportErrorModule        , only : reportError, debugFileOnly_p, error_p
    use FiDeviceFunctionInterface, only : FiDF_extfunc
    use FFaFilePathInterface     , only : ffa_checkPath
    use FFaCmdLineArgInterface   , only : ffa_cmdlinearg_getstring

    type(EngineType), intent(inout) :: engines(:)
    type(SensorType), intent(in)    :: sensors(:)
    integer         , intent(out)   :: err

    !! Local variables
    integer :: i, j, stat, sensorId, extFunc(maxChannel)
    logical :: multiArg
    character(len=lfnam_p) :: fileName
    character(len=1024)    :: labels

    !! --- Logic section ---

    err = 0

    do i = 1, size(engines)

       multiArg = .false.
       sensorId = argSensor(i)
       if (associated(argSensors(i)%values)) then
          multiArg = .true.
          allocate(engines(i)%args(size(argSensors(i)%values)),STAT=stat)
       else if (sensorId /= 0) then
          allocate(engines(i)%args(1),STAT=stat)
       else
          allocate(engines(i)%args(0),STAT=stat)
       end if
       if (stat /= 0) then
          err = allocationError('InitiateEngines2')
          return
       end if

       stat = err
       do j = 1, size(engines(i)%args)
          if (multiArg) sensorId = argSensors(i)%values(j)
          if (sensorId == 0) then
             nullify(engines(i)%args(j)%p) ! Assume zero-valued argument
          else
             engines(i)%args(j)%p => GetPtrToId(sensors,sensorId)
             if (.not. associated(engines(i)%args(j)%p)) then
                err = err - 1
                call ReportInputError ('ENGINE',i,engines(i)%id)
             end if
          end if
       end do
       if (err < stat) then
          deallocate(engines(i)%args)
          nullify(engines(i)%args)
       end if

       if (multiArg) deallocate(argSensors(i)%values)

       if (sensorId == 0 .and. associated(engines(i)%func)) then
          if (engines(i)%func%type == DEVICE_FUNCTION_p) then
             extFunc(engines(i)%func%intParameters(4)) = i
          end if
       end if

    end do

    deallocate(argSensor,argSensors)

    call ffa_cmdlinearg_getstring ('externalfuncfile',fileName)
    if (fileName /= '' .and. maxChannel > 0 .and. err == 0) then
       !! Read the external function values from file.
       !! Find the labels which are used to identify the columns.
       labels = '<"'//engines(extFunc(1))%tag
       do i = 2, maxChannel
          labels = trim(labels)//'","'//engines(extFunc(i))%tag
       end do
       write(lterm,*) 'Read external function values from file ',trim(fileName)
       if (len_trim(labels) > maxChannel) then
          labels = trim(labels)//'">'
          write(lterm,*) 'Using channel labels ',trim(labels)
       else
          labels = ''
       end if
       call ffa_checkPath (fileName)
       call FiDF_extfunc (err,trim(fileName),trim(labels))
       if (err < 0) then
          call reportError (error_p,'Failed to initialize external functions')
       end if
    end if

    if (err < 0) call reportError (debugFileOnly_p,'InitiateEngines2')

  end subroutine InitiateEngines2

end module initiateFunctionTypeModule
