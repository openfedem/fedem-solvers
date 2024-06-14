!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file rdbModule.f90
!> @brief Data type and subroutines for writing FEDEM result files.

!!==============================================================================
!> @brief Module with data type and subroutines for writing FEDEM result files.

module rdbModule

  use KindModule, only : lfnam_p, nbi_p, nbs_p, nbd_p, sp, i8

  implicit none

  private

  !> @brief Data type for a binary result database file.
  type RDBType
     character(lfnam_p) :: fileName !< Name of the results database file
     integer            :: fileNumb !< File number used by the binaryDB routines
     integer            :: bufSize  !< Size of file buffer (in bytes)
     integer(i8)        :: nBytes   !< Byte counter for debugging
     integer            :: nBstep   !< Number of bytes written per time step
     integer            :: nVar,nIG !< Number of variable- and item group defs
     logical            :: tmpBin   !< Write binary data to a temporary file
     logical            :: empty    !< .true. until some binary data is written
  end type RDBType

  integer, save, private :: incr  = -1 !< File increment number
  integer, save, public  :: ivard = 65 !< Temporary variable description file
  integer, save, public  :: iitem = 66 !< Temporary item group description file
  integer, save, public  :: idatd = 67 !< Temporary data description file

  !!TODO,kmo: Check if the below still is an issue or if this can be refactored.
  !! Note that since ivard, iitem and idatd are global static variables
  !! (and not members of RDBType) we can only write one file header at once.
  !! Trying to put these parameters as variables in RDBType only results in a
  !! compiler error on SGI on write statements like "write(rdb%ivard,600) ..."

  real(sp), allocatable, save :: work(:) !< For double-to-single casting

  public :: nbi_p, nbs_p, nbd_p
  public :: RDBType, nullifyRDB
  public :: openHeaderFiles, writeTimeStepHeader
  public :: openRDBfile, flushRDBfile, closeRDBfile
  public :: writeRDB, writeTimeStepDB
  public :: writeVarDef, writeItGDef
  public :: checkStepSize

  !> @brief Writes a result quantity to the specified results database file.
  interface writeRDB
     module procedure writeBoolVarDB
     module procedure writeIntVarDB
     module procedure writeIntArrayDB
     module procedure writeFloatArray1DB
     module procedure writeFloatArray2DB
     module procedure writeDoubleVarDB
     module procedure writeDoubleArray1DB
     module procedure writeDoubleArray2DB
  end interface


contains

  !!============================================================================
  !> @cond FULL_DOC
  !> @brief Prints an error message on write failure.
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 22 Mar 2001

  subroutine reportWriteError (prg,rdb)

    use reportErrorModule, only : reportError, error_p

    character(len=*), intent(in) :: prg
    type(RDBType)   , intent(in) :: rdb

    !! Local variables
    character(len=64) :: addMsg

    !! --- Logic section ---

    write(addMsg,"('Number of bytes successfully written:',I10)") rdb%nBytes
    call reportError (error_p,'Error writing to binary file '// &
         &            trim(rdb%fileName),addMsg,addString=prg)

  end subroutine reportWriteError
  !> @endcond


  !!============================================================================
  !> @brief Checks that the estimated number of bytes per time is correct.
  !>
  !> @param[in] rdb The result database file to consider
  !> @param[in] nByteStep Estimated number of bytes per time step for this file
  !> @param[out] ierr Error flag
  !>
  !> @details An error message is issued if the actual amount written per step
  !> does not match the estimated time step size.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 22 Mar 2001

  subroutine checkStepSize (rdb,nByteStep,ierr)

    use reportErrorModule, only : reportError, error_p

    type(RDBType), intent(in)  :: rdb
    integer      , intent(in)  :: nByteStep
    integer      , intent(out) :: ierr

    !! Local variables
    character(len=64) :: addMsg

    !! --- Logic section ---

    if (rdb%bufSize > 0) then
       ierr = -abs(nByteStep - rdb%nBstep)
    else
       ierr = 0
    end if
    if (ierr == 0) return

    write(addMsg,600) rdb%nBstep, nByteStep
    call reportError (error_p, &
         'You have detected an inconsistency in the results management.', &
         'Please give details to developers@openfedem.org','This error '// &
         'can be bypassed by switching to automatic flushing of result files', &
         '(specify -flushinc -1 as an additional solver option).', &
         addString=trim(addMsg)//',   File: '//rdb%fileName)

600 format('Estimated step size:',I8,',   Actual size:',I8)

  end subroutine checkStepSize


  !!============================================================================
  !> @brief Initializes a rdbmodule::rdbtype object.
  !>
  !> @param[out] rdb The result database file to consider
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 22 Mar 2001

  subroutine nullifyRDB (rdb)

    type(RDBType), intent(out) :: rdb

    !! --- Logic section ---

    rdb%fileName = ''
    rdb%fileNumb = -1
    rdb%bufSize  = 0
    rdb%nBytes   = 0_i8
    rdb%nBstep   = 0
    rdb%nVar     = 0
    rdb%nIG      = 0
    rdb%tmpBin   = .false.
    rdb%empty    = .true.

  end subroutine nullifyRDB


  !!============================================================================
  !> @brief Opens the temporary header files and write the meta data section.
  !>
  !> @param[in] prog Name of the program module that writes the file to open
  !> @param[out] ierr Error flag
  !> @param[in] comment Additional description to write in the meta data section
  !> @param[in] mFile Name of the model file to which the written file belong
  !> @param[in] mName Name of the model to which the written file belong
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 20 Nov 2000

  subroutine openHeaderFiles (prog,ierr,comment,mFile,mName)

    use versionModule          , only : progVer, buildDate, getCurrentDate
    use fileUtilitiesModule    , only : findUnitNumber
    use reportErrorModule      , only : reportError, error_p
    use computerConfigInterface, only : getUserName, getComputerConfig

    character(len=*)         , intent(in)  :: prog
    integer                  , intent(out) :: ierr
    character(len=*),optional, intent(in)  :: comment, mFile, mName

    !! Local variables
    character(len=16) :: user
    character(len=96) :: chid
    character(len=32) :: date

    !! --- Logic section ---

    ivard = findUnitNumber(65)
    open(ivard,iostat=ierr,status='SCRATCH')
    if (ierr == 0) then
       iitem = findUnitNumber(66)
       open(iitem,iostat=ierr,status='SCRATCH')
       if (ierr == 0) then
          idatd = findUnitNumber(67)
          open(idatd,iostat=ierr,status='SCRATCH')
          if (ierr == 0) then
             call getUserName (user)
             call getComputerConfig (chid)
             call getCurrentDate (date)
             if (present(mFile)) then
                write(ivard,600) 'AssociatedModelFileName', trim(mFile)
             end if
             if (present(mName)) then
                write(ivard,600) 'ModelName              ', trim(mName)
             end if
             if (present(comment)) then
                write(ivard,600) 'InformationText        ', trim(comment)
             end if
             write(ivard,600) 'User                   ', trim(user)
             write(ivard,600) 'Computer               ', trim(chid)
             write(ivard,600) 'DateTime               ', trim(date)
             write(ivard,600) 'UsedTime               ', '00:00:00.00'
             write(ivard,600) 'Module                 ', trim(prog)
             write(ivard,600) 'ModuleVersion          ', &
                  &           trim(progVer)//' '//trim(buildDate)
             write(ivard,601) 'VARIABLES:'
             write(idatd,601) 'DATABLOCKS:'
             return
          end if
       end if
    end if

    ierr = -1
    call reportError (error_p,'Could not open temporary header files', &
         &            addString='openHeaderFiles')

600 format(1x,a,' = ',a,';')
601 format(a)

  end subroutine openHeaderFiles


  !!============================================================================
  !> @brief Opens a binary results database file.
  !>
  !> @param[in] rdb The result database file to open
  !> @param[out] ierr Error flag
  !> @param[in] fileName Absolute path to the result database file
  !> @param[in] fileTag Tag to write in the first line identifying the file type
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 27 Nov 2000

  subroutine openRDBfile (rdb,ierr,fileName,fileTag)

    use kindModule            , only : maxInt_p
    use versionModule         , only : getCurrentDate
    use reportErrorModule     , only : reportError, error_p, noteFileOnly_p
    use binaryDBInterface     , only : openBinaryDB, closeBinaryDB, setBufDB
    use binaryDBInterface     , only : writeTagDB, copyBinaryDB, flushBinaryDB
    use binaryDBInterface     , only : write_p, append_p
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getint

    type(RDBType), intent(inout)        :: rdb
    integer      , intent(out)          :: ierr
    character*(*), intent(in), optional :: fileName, fileTag

    !! Local variables
    integer           :: i, tmpFile, tmpBytes
    character(len=32) :: ch

    !! --- Logic section ---

    if (rdb%tmpBin) then
       tmpFile = rdb%fileNumb
       if (rdb%nBytes > int(maxInt_p,i8)) then
          call reportError (error_p,'Can not copy temporary RDB-file '// &
               &            'larger than 2GB',addString='openRDBfile')
          ierr = -int(min(int(maxInt_p,i8),rdb%nBytes-int(maxInt_p,i8)))
          return
       else
          tmpBytes = int(rdb%nBytes)
       end if
    end if

    if (present(fileName)) then

       !! Create a new file name using the file increment number
       if (incr < 0) then
          call ffa_cmdlinearg_getint ('rdbinc',incr)
       else if (incr > 0) then
          incr = incr + 1
       end if
       if (incr > 0) then
          write(ch,'(i6)') incr
       else
          call getCurrentDate (ch,.true.)
       end if
       i = index(fileName,'.',.TRUE.)
       if (i > 1) then
          rdb%fileName = fileName(1:i-1)//'_'//trim(adjustl(ch))//fileName(i:)
       else
          rdb%fileName = trim(fileName)//'_'//trim(adjustl(ch))
       end if

       !! Open a new database file for writing
       call openBinaryDB (rdb%fileName,write_p,rdb%fileNumb,ierr)

    else

       !! Open an existing database file for append
       call openBinaryDB (rdb%fileName,append_p,rdb%fileNumb,ierr)

    end if

    if (ierr < 0) then
       call reportError (error_p,'Could not open RDB-file '//rdb%fileName, &
            &            addString='openRDBfile')
       return
    end if

    if (present(fileTag)) then

       !! Write file tag and copy the ASCII file header to the binary file
       call writeTagDB (rdb%fileNumb,fileTag,0,ierr)
       call copyHeaderToBinaryFile (rdb%fileNumb,ierr)
       if (ierr < 0) then
          call reportError (error_p,'Could not write header to RDB-file '// &
               &            rdb%fileName,addString='openRDBfile')
          return
       end if

       rdb%nBytes = ierr
       write(ch,"(3I10)") rdb%nBytes,rdb%nBstep,rdb%bufSize
       call reportError (noteFileOnly_p,'Opening RDB-file: '//rdb%fileName, &
            &            '   Header size   ='//ch(1:10)//' Bytes', &
            &            '   Data per step ='//ch(11:20)//' Bytes', &
            &            '   Buffer size   ='//ch(21:30)//' Bytes')

       !! Flush file header to disk
       call flushBinaryDB (rdb%fileNumb,ierr)
       if (ierr < 0) then
          call reportError (error_p,'Could not flush RDB-file '// &
               &            rdb%fileName,addString='openRDBfile')
          return
       end if

    end if

    if (rdb%tmpBin) then

       !! Binary data have already been written to a temporary file.
       !! Now copy it into the permanent file.
       !! Note that the file is still considered empty if only time step data
       !! has been written. The rdb%empty flag must therefore be set explicitly
       !! to .false. if more binary data is to be written later. If not,
       !! the file will be deleted when closed.

       if (tmpBytes > 0) then

          call copyBinaryDB (rdb%fileNumb,tmpFile,tmpBytes,ierr)
          if (ierr < 0) then
             call reportError (error_p,'Could not copy temporary RDB-file', &
                  &            addString='openRDBfile')
             return
          end if

          rdb%nBytes = rdb%nBytes + ierr
          rdb%empty  = tmpBytes <= nbi_p + nbd_p ! Only time step data?
       end if
       rdb%tmpBin = .false.

       call closeBinaryDB (tmpFile,ierr)
       if (ierr < 0) then
          call reportError (error_p,'Could not close temporary RDB-file', &
               &            addString='openRDBfile')
          return
       end if

    end if

    !! Allocate a buffer for this file, if requested
    call setBufDB (rdb%fileNumb,rdb%bufSize,ierr)
    if (ierr < 0) then
       call reportError (error_p,'Could not allocate RDB-file buffer', &
            &            addString='openRDBfile')
    end if

  end subroutine openRDBfile


  !!============================================================================
  !> @brief Copies contents of the temporary header files to the binary file.
  !>
  !> @param[in] ifile File handle (0-based index) for the binary file
  !> @param nBytes Number of bytes written to the file, negative on error
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 27 Nov 2000

  subroutine copyHeaderToBinaryFile (ifile,nBytes)

    use binaryDBInterface, only : writeCharDB

    integer, intent(in)    :: ifile
    integer, intent(inout) :: nBytes

    !! Local variables
    integer            :: ierr
    character(len=256) :: chline ! NB: Max length of a header line is 256 !

    !! --- Logic section ---

    ierr = nBytes
    rewind(ivard)
    do while (ierr >= 0)
       read(ivard,'(a)',end=100) chline
       call writeCharDB (ifile,trim(chline)//char(13),ierr)
       nBytes = nBytes + ierr
    end do
100 close(ivard)

    rewind(iitem)
    do while (ierr >= 0)
       read(iitem,'(a)',end=200) chline
       call writeCharDB (ifile,trim(chline)//char(13),ierr)
       nBytes = nBytes + ierr
    end do
200 close(iitem)

    rewind(idatd)
    do while (ierr >= 0)
       read(idatd,'(a)',end=300) chline
       call writeCharDB (ifile,trim(chline)//char(13),ierr)
       nBytes = nBytes + ierr
    end do
300 close(idatd)

    if (ierr >= 0) then
       call writeCharDB (ifile,'DATA:',ierr)
       nBytes = nBytes + ierr
    end if

    if (ierr < 0) nBytes = ierr

  end subroutine copyHeaderToBinaryFile


  !!============================================================================
  !> @brief Writes a variable definition.
  !>
  !> @param rdb The result database file to write to
  !> @param iVar Variable identifier (1-based index)
  !> @param[in] iFile File unit number to write a reference to the variable to
  !> @param[in] cId Variable description
  !> @param[in] cDim Dimension of the variable
  !> @param[in] cType Type of the variable
  !> @param[in] cComp1 Description of first dimension components
  !> @param[in] cComp2 Description of second dimension components
  !> @param[in] nBits Number of bits in each variable component (32 or 64)
  !>
  !> @details This subroutine writes a variable definition to the variable
  !> definition section of the results file header, and optionally a reference
  !> to it in the item group- or data definition section.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 30 Apr 2002

  subroutine writeVarDef (rdb,iVar,iFile,cId,cDim,cType,cComp1,cComp2,nBits)

    type(RDBType)   , intent(inout)       :: rdb
    integer         , intent(inout)       :: iVar
    integer         , intent(in)          :: iFile
    character(len=*), intent(in)          :: cId
    character(len=*), intent(in),optional :: cDim, cType, cComp1, cComp2
    integer         , intent(in),optional :: nBits

    !! Local variables
    integer :: iRef, nB

    !! --- Logic section ---

    if (iVar == 0) then ! This variable has not yet been defined
       rdb%nVar = rdb%nVar + 1
       iRef = rdb%nVar
    else
       iRef = iVar
    end if

    if (iVar == 0) then ! Write the variable definition
       iVar = iRef
       if (.not. present(nBits)) then
          nB = 32 ! Default is single precision
       else
          nB = nBits
       end if
       select case (iVar)
       case ( 0: 9); write(ivard,601,advance='NO') iVar
       case (10:99); write(ivard,602,advance='NO') iVar
       case default; write(ivard,603,advance='NO') iVar
       end select
       if (.not. present(cDim)) then
          write(ivard,610) cId,'NONE',nB,'SCALAR'
       else if (.not. present(cType)) then
          write(ivard,610) cId,cDim,nB,'SCALAR'
       else if (cType == 'NUMBER' .or. nB < 8) then
          write(ivard,620) cId,cDim,8*nbi_p,cType
       else if (.not. present(cComp1)) then
          write(ivard,610) cId,cDim,nB,cType
       else if (.not. present(cComp2)) then
          write(ivard,630) cId,cDim,nB,cType,cComp1
       else
          write(ivard,640) cId,cDim,nB,cType,cComp1,cComp2
       end if
    end if

    if (iFile > 0) then ! Write a reference to the variable
       select case (iVar)
       case ( 0: 9); write(iFile,601,advance='NO') iVar
       case (10:99); write(iFile,602,advance='NO') iVar
       case default; write(iFile,603,advance='NO') iVar
       end select
       write(iFile,"('>')",advance='NO')
    end if

601 format('<',i1)
602 format('<',i2)
603 format('<',i3)
610 format(';"',a,'";',a,';FLOAT;',i2,';',a,'>')
620 format(';"',a,'";',a,';INT;',i2,';',a,'>')
630 format(';"',a,'";',a,';FLOAT;',i2,';',a,';(3);((',a,'))>')
640 format(';"',a,'";',a,';FLOAT;',i2,';',a,';(3,4);((',a,'),(',a,'))>')

  end subroutine writeVarDef


  !!============================================================================
  !> @brief Writes the beginning of an item group definition.
  !>
  !> @param rdb The result database file to write to
  !> @param itemGroup Item group identifier (1-based index)
  !> @param[out] iFile File unit number for the item group definition,
  !> if zero, this item group definition has already been written
  !> @param[in] cId Item group description
  !>
  !> @details This subroutine writes the beginning of an item group definition
  !> to the item group definition section of the results file header,
  !> and a reference to it in the data section.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 30 Apr 2002

  subroutine writeItGDef (rdb,itemGroup,iFile,cId)

    type(RDBType)   , intent(inout) :: rdb
    integer         , intent(inout) :: itemGroup
    integer         , intent(out)   :: iFile
    character(len=*), intent(in)    :: cId

    !! --- Logic section ---

    if (itemGroup == 0) then ! This item group has not yet been defined
       rdb%nIG   = rdb%nIG + 1
       itemGroup = rdb%nIG
       iFile = iitem
    else if (itemGroup < 0) then ! Item group defined, but not written yet
       itemGroup = -itemGroup
       iFile = iitem
    else
       iFile = 0
    end if

    if (iFile > 0) then
       !! Write the item group definition
       select case (itemGroup)
       case ( 0: 9); write(iFile,601,advance='NO') itemGroup
       case (10:99); write(iFile,602,advance='NO') itemGroup
       case default; write(iFile,603,advance='NO') itemGroup
       end select
       write(iFile,600,advance='NO') cId
    end if

    !! Write a reference to the item group
    select case (itemGroup)
    case ( 0: 9); write(idatd,601,advance='NO') itemGroup
    case (10:99); write(idatd,602,advance='NO') itemGroup
    case default; write(idatd,603,advance='NO') itemGroup
    end select
    write(idatd,"(']')",advance='NO')

600 format(';"',a,'";')
601 format('[',i1)
602 format('[',i2)
603 format('[',i3)

  end subroutine writeItGDef


  !!============================================================================
  !> @brief Writes time step definition to the temporary header files.
  !>
  !> @param rdb The result database file to write to
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 10 Jan 2003

  subroutine writeTimeStepHeader (rdb)

    type(RDBType), intent(inout) :: rdb

    !! Local variables
    integer :: idTimeStep, idPhysTime

    !! --- Logic section ---

    idTimeStep = 0
    idPhysTime = 0
    call writeVarDef (rdb,idTimeStep,idatd,'Time step number','NONE','NUMBER')
    call writeVarDef (rdb,idPhysTime,idatd,'Physical time','TIME',nBits=64)
    write(idatd,*)

    rdb%nBytes = rdb%nBytes + nbi_p + nbd_p

  end subroutine writeTimeStepHeader


  !!============================================================================
  !> @brief Writes current time step to the specified results database file.
  !>
  !> @param rdb The result database file to write to
  !> @param[in] iStep Time step number
  !> @param[in] time Current time
  !> @param[out] ierr Error flag
  !> @param[in] openTmp If present and .true., a temporary binary file is
  !> opened if the binary results database file is not opened yet
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 20 Nov 2000

  subroutine writeTimeStepDB (rdb,iStep,time,ierr,openTmp)

    use kindModule            , only : dp, i8, maxInt_p
    use reportErrorModule     , only : getErrorFile, internalError, reportError
    use reportErrorModule     , only : error_p, debugFileOnly_p
    use binaryDBInterface     , only : openBinaryDB, tmp_p
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getint

    type(RDBType)   , intent(inout) :: rdb
    integer(i8)     , intent(in)    :: iStep
    real(dp)        , intent(in)    :: time
    integer         , intent(out)   :: ierr
    logical,optional, intent(in)    :: openTmp

    !! Local variables
    integer :: lpu, iprint, jStep

    !! --- Logic section ---

    if (rdb%fileNumb < 0) then
       if (present(openTmp)) rdb%tmpBin = openTmp
       if (.not. rdb%tmpBin) then
          ierr = internalError('writeTimeStepDB: RDB-file has not been opened')
          return
       else
          !! Open a (buffered) temporary database file for writing
          call openBinaryDB ('dummy',tmp_p,rdb%fileNumb,ierr)
          if (ierr == 0) call setBufDB (rdb%fileNumb,rdb%bufSize,ierr)
          if (ierr < 0) then
             call reportError (error_p,'Could not open temporary RDB-file', &
                  &            addString='writeTimeStepDB')
             return
          end if
       end if
    end if

    call ffa_cmdlinearg_getint ('debug',iprint)
    if (iprint > 5 .and. rdb%nBytes > 0_i8 .and. .not.rdb%tmpBin) then
       lpu = getErrorFile()
       write(lpu,600) trim(rdb%fileName),rdb%nBytes
       rdb%nBytes = 0_i8
    end if

    !! TODO: Remove this when the results database (including FFrLib)
    !!       fully supports 64-bit integer time step counters.
    if (iStep > int(maxInt_p,i8)) then
       call reportError (error_p,'The time step counter is exhausted.', &
            '64-bit counter is not yet implemented in the result database.', &
            'Deactivate all result output and restart to continue.', &
            addString='writeTimeStepDB')
       ierr = -99
       return
    else
       jStep = int(iStep) ! Cast to 32-bit integer
    end if

    ierr = 0
    call writeRDB (rdb,jStep,ierr)
    call writeRDB (rdb,time,ierr,.true.)
    if (ierr < 0) then
       call reportError (debugFileOnly_p,'writeTimeStepDB')
    else if (.not. rdb%tmpBin) then
       rdb%empty = .false.
    end if

600 format('  RDB-file ',a,' :',i10,' bytes written')

  end subroutine writeTimeStepDB


  !!============================================================================
  !> @brief Writes a bolean variable to the specified results database file.
  !>
  !> @param rdb The result database file to write to
  !> @param[in] var The variable value to write
  !> @param ierr Error flag
  !>
  !> @callgraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 16 Jun 2003

  subroutine writeBoolVarDB (rdb,var,ierr)

    type(RDBType), intent(inout) :: rdb
    logical      , intent(in)    :: var
    integer      , intent(inout) :: ierr

    !! --- Logic section ---

    if (var) then
       call writeIntVarDB (rdb,1,ierr)
    else
       call writeIntVarDB (rdb,0,ierr)
    end if

  end subroutine writeBoolVarDB


  !!============================================================================
  !> @brief Writes an integer variable to the specified results database file.
  !>
  !> @param rdb The result database file to write to
  !> @param[in] var The variable value to write
  !> @param ierr Error flag
  !>
  !> @callgraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 6 Dec 2000

  subroutine writeIntVarDB (rdb,var,ierr)

    use binaryDBInterface, only : writeIntDB

    type(RDBType), intent(inout) :: rdb
    integer      , intent(in)    :: var
    integer      , intent(inout) :: ierr

    !! --- Logic section ---

    if (ierr < 0) return

    call writeIntDB (rdb%fileNumb,var,1,ierr)
    if (ierr > 0) then
       rdb%nBytes = rdb%nBytes + ierr
    else if (ierr < 0) then
       call reportWriteError ('writeIntVarDB',rdb)
    end if

  end subroutine writeIntVarDB


  !!============================================================================
  !> @brief Writes an integer array to the specified results database file.
  !>
  !> @param rdb The result database file to write to
  !> @param[in] array The array to write
  !> @param ierr Error flag
  !>
  !> @callgraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 25 Oct 2008

  subroutine writeIntArrayDB (rdb,array,ierr)

    use binaryDBInterface, only : writeIntDB

    type(RDBType), intent(inout) :: rdb
    integer      , intent(in)    :: array(:)
    integer      , intent(inout) :: ierr

    !! --- Logic section ---

    if (ierr < 0) return

    call writeIntDB (rdb%fileNumb,array(1),size(array),ierr)
    if (ierr > 0) then
       rdb%nBytes = rdb%nBytes + ierr
    else if (ierr < 0) then
       call reportWriteError ('writeIntArrayDB',rdb)
    end if

  end subroutine writeIntArrayDB


  !!============================================================================
  !> @brief Writes a 1D float array to the specified results database file.
  !>
  !> @param rdb The result database file to write to
  !> @param[in] array The array to write
  !> @param ierr Error flag
  !> @param[in] writeAsDouble Currently unused (gives warning)
  !>
  !> @callgraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 10 Jan 2003

  subroutine writeFloatArray1DB (rdb,array,ierr,writeAsDouble)

    use binaryDBInterface, only : writeFloatDB
    use reportErrorModule, only : reportError, warningFileOnly_p

    type(RDBType)   , intent(inout) :: rdb
    real(sp)        , intent(in)    :: array(:)
    integer         , intent(inout) :: ierr
    logical,optional, intent(in)    :: writeAsDouble

    !! --- Logic section ---

    if (ierr < 0) return

    if (present(writeAsDouble)) then
       if (writeAsDouble) then
          call reportError (warningFileOnly_p,'Ignoring writeAsDouble request',&
               &            addString='writeFloatArray1DB')
       end if
    end if

    call writeFloatDB (rdb%fileNumb,array(1),size(array),ierr)
    if (ierr > 0) then
       rdb%nBytes = rdb%nBytes + ierr
    else if (ierr < 0) then
       call reportWriteError ('writeFloatArray1DB',rdb)
    end if

  end subroutine writeFloatArray1DB


  !!============================================================================
  !> @brief Writes a 2D float array to the specified results database file.
  !>
  !> @param rdb The result database file to write to
  !> @param[in] array The array to write
  !> @param ierr Error flag
  !>
  !> @callgraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 10 Feb 2017

  subroutine writeFloatArray2DB (rdb,array,ierr)

    use binaryDBInterface, only : writeFloatDB

    type(RDBType), intent(inout) :: rdb
    real(sp)     , intent(in)    :: array(:,:)
    integer      , intent(inout) :: ierr

    !! --- Logic section ---

    if (ierr < 0) return

    call writeFloatDB (rdb%fileNumb,array(1,1),size(array),ierr)
    if (ierr > 0) then
       rdb%nBytes = rdb%nBytes + ierr
    else if (ierr < 0) then
       call reportWriteError ('writeFloatArray2DB',rdb)
    end if

  end subroutine writeFloatArray2DB


  !!============================================================================
  !> @brief Writes a double variable to the specified results database file.
  !>
  !> @param rdb The result database file to write to
  !> @param[in] var The variable value to write
  !> @param ierr Error flag
  !> @param[in] writeAsDouble If present and .true., write as 64-bit real,
  !> otherwise cast to 32-bit real
  !>
  !> @callgraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 6 Dec 2000

  subroutine writeDoubleVarDB (rdb,var,ierr,writeAsDouble)

    use kindModule       , only : dp
    use binaryDBInterface, only : writeDoubleDB, writeFloatDB

    type(RDBType)   , intent(inout) :: rdb
    real(dp)        , intent(in)    :: var
    integer         , intent(inout) :: ierr
    logical,optional, intent(in)    :: writeAsDouble

    !! Local variables
    logical :: writeDoubles

    !! --- Logic section ---

    if (ierr < 0) return

    if (present(writeAsDouble)) then
       writeDoubles = writeAsDouble
    else
       writeDoubles = .false. ! Write in single precision is the default
    end if

    if (writeDoubles) then
       call writeDoubleDB (rdb%fileNumb,var,1,ierr)
    else
       call writeFloatDB (rdb%fileNumb,real(var,sp),1,ierr)
    end if

    if (ierr > 0) then
       rdb%nBytes = rdb%nBytes + ierr
    else if (ierr < 0) then
       call reportWriteError ('writeDoubleVarDB',rdb)
    end if

  end subroutine writeDoubleVarDB


  !!============================================================================
  !> @brief Writes a 1D double array to the specified results database file.
  !>
  !> @param rdb The result database file to write to
  !> @param[in] array The array to write
  !> @param ierr Error flag
  !> @param[in] writeAsDouble If present and .true., write as 64-bit reals,
  !> otherwise cast to 32-bit reals
  !>
  !> @callgraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 6 Dec 2000

  subroutine writeDoubleArray1DB (rdb,array,ierr,writeAsDouble)

    use kindModule       , only : dp
    use binaryDBInterface, only : writeDoubleDB, writeFloatDB
    use reportErrorModule, only : allocationError

    type(RDBType)   , intent(inout) :: rdb
    real(dp)        , intent(in)    :: array(:)
    integer         , intent(inout) :: ierr
    logical,optional, intent(in)    :: writeAsDouble

    !! Local variables
    integer :: i, ndat
    logical :: writeDoubles

    !! --- Logic section ---

    if (ierr < 0) return

    if (present(writeAsDouble)) then
       writeDoubles = writeAsDouble
    else
       writeDoubles = .false. ! Write in single precision is the default
    end if

    if (writeDoubles) then
       call writeDoubleDB (rdb%fileNumb,array(1),size(array),ierr)
    else
       !! Type-cast array into single precision using scratch array work
       ndat = size(array)
       if (.not. allocated(work)) then
          allocate(work(ndat),STAT=ierr)
       else if (size(work) < ndat) then
          deallocate(work)
          allocate(work(ndat),STAT=ierr)
       else
          ierr = 0
       end if
       if (ierr /= 0) then
          ierr = allocationError('writeDoubleArray1DB')
          return
       end if
       do i = 1, ndat
          work(i) = real(array(i),sp)
       end do
       call writeFloatDB (rdb%fileNumb,work(1),ndat,ierr)
    end if

    if (ierr > 0) then
       rdb%nBytes = rdb%nBytes + ierr
    else if (ierr < 0) then
       call reportWriteError ('writeDoubleArray1DB',rdb)
    end if

  end subroutine writeDoubleArray1DB


  !!============================================================================
  !> @brief Writes a 2D double array to the specified results database file.
  !>
  !> @param rdb The result database file to write to
  !> @param[in] array The array to write
  !> @param[in] writeAsDouble If present and .true., write as 64-bit reals,
  !> otherwise cast to 32-bit reals
  !> @param ierr Error flag
  !>
  !> @callgraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 6 Dec 2000

  subroutine writeDoubleArray2DB (rdb,array,ierr,writeAsDouble)

    use kindModule       , only : dp
    use binaryDBInterface, only : writeDoubleDB, writeFloatDB
    use reportErrorModule, only : allocationError

    type(RDBType)   , intent(inout) :: rdb
    real(dp)        , intent(in)    :: array(:,:)
    integer         , intent(inout) :: ierr
    logical,optional, intent(in)    :: writeAsDouble

    !! Local variables
    integer :: i, j, ndat, dim1
    logical :: writeDoubles

    !! --- Logic section ---

    if (ierr < 0) return

    if (present(writeAsDouble)) then
       writeDoubles = writeAsDouble
    else
       writeDoubles = .false. ! Write in single precision is the default
    end if

    if (writeDoubles) then
       call writeDoubleDB (rdb%fileNumb,array(1,1),size(array),ierr)
    else
       !! Type-cast array into single precision using scratch array work
       ndat = size(array)
       if (.not. allocated(work)) then
          allocate(work(ndat),STAT=ierr)
       else if (size(work) < ndat) then
          deallocate(work)
          allocate(work(ndat),STAT=ierr)
       else
          ierr = 0
       end if
       if (ierr /= 0) then
          ierr = allocationError('writeDoubleArray2DB')
          return
       end if
       dim1 = size(array,1)
       do i = 1, size(array,2)
          do j = 1, dim1
             work(j+(i-1)*dim1) = real(array(j,i),sp)
          end do
       end do
       call writeFloatDB (rdb%fileNumb,work(1),ndat,ierr)
    end if

    if (ierr > 0) then
       rdb%nBytes = rdb%nBytes + ierr
    else if (ierr < 0) then
       call reportWriteError ('writeDoubleArray2DB',rdb)
    end if

  end subroutine writeDoubleArray2DB


  !!============================================================================
  !> @brief Flushes the specified results database file to disk.
  !>
  !> @param rdb The result database file to flush
  !> @param ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 10 Jan 2003

  subroutine flushRDBfile (rdb,ierr)

    use reportErrorModule, only : reportError, error_p
    use binaryDBInterface, only : flushBinaryDB

    type(RDBType), intent(inout) :: rdb
    integer      , intent(inout) :: ierr

    !! --- Logic section ---

    if (rdb%fileNumb < 0 .or. ierr < 0) return

    call flushBinaryDB (rdb%fileNumb,ierr)
    if (ierr < 0) then
       call reportError (error_p,'Could not flush RDB-file '//rdb%fileName, &
            &            addString='flushRDBfile')
    end if

  end subroutine flushRDBfile


  !!============================================================================
  !> @brief Closes the specified results database file.
  !>
  !> @param rdb The result database file to close
  !> @param ierr Error flag
  !> @param[in] retainIfEmpty If present and .true., don't delete empty file
  !> @param[out] haveData If present, set to .true. unless the file is empty
  !>
  !> @details The file is deleted if it is empty when closed (i.e., no binary
  !> data has been written to it) unless @a retainIfEmpty is .true.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 20 Nov 2000

  subroutine closeRDBfile (rdb,ierr,retainIfEmpty,haveData)

    use timerModule      , only : getCurrentTime
    use reportErrorModule, only : getErrorFile, reportError, error_p, warning_p
    use binaryDBInterface, only : openBinaryDB, closeBinaryDB, putCharDB
    use binaryDBInterface, only : reawri_p

    type(RDBType)   , intent(inout) :: rdb
    integer         , intent(inout) :: ierr
    logical,optional, intent(in)    :: retainIfEmpty
    logical,optional, intent(out)   :: haveData

    !! Local variables
    logical           :: deleteIt
    integer           :: lerr, lpu, hours, minutes, seconds
    real(sp)          :: usedTime(2), hundreds
    character(len=12) :: chUsed
    real(sp),external :: cpusec

    !! --- Logic section ---

    if (present(haveData)) haveData = .false.
    if (rdb%fileNumb < 0) return ! already closed or never opened

    if (allocated(work)) deallocate(work)

    if (.not. present(retainIfEmpty)) then
       deleteIt = rdb%empty
    else if (retainIfEmpty) then
       deleteIt = .false.
    else
       deleteIt = rdb%empty
    end if

    lpu = getErrorFile()
    if (deleteIt) then
       write(lpu,600) trim(rdb%fileName)
    else if (rdb%nBytes > 0_i8) then
       write(lpu,601) trim(rdb%fileName),rdb%nBytes
    end if

    call closeBinaryDB (rdb%fileNumb,deleteIt,lerr)
    if (lerr < 0) then
       ierr = ierr + lerr
       call reportError (error_p,'Could not close RDB-file '//rdb%fileName, &
            &            addString='closeRDBfile')
       return
    else if (deleteIt) then
       rdb%fileNumb = -99
       rdb%nBytes = 0_i8
       return
    end if

    !! Re-open the file in read-write mode
    call openBinaryDB (rdb%fileName,reawri_p,rdb%fileNumb,lerr)
    if (lerr < 0) then
       ierr = ierr + lerr
       call reportError (error_p,'Could not re-open RDB-file '//rdb%fileName, &
            &            addString='closeRDBfile')
       return
    end if

    !! Write the used CPU-time into the file-info section
    usedTime = getCurrentTime(1)
    hours    = int(usedTime(2))/3600
    minutes  = int(usedTime(2))/60 - 60*hours
    seconds  = int(usedTime(2)) - 60*minutes - 3600*hours
    hundreds = usedTime(2) - int(usedTime(2))
    write(chUsed,"(i3,':',i2.2,':',i2.2,f3.2)") hours,minutes,seconds,hundreds
    call putCharDB (rdb%fileNumb,'UsedTime                =',chUsed,lerr)
    if (lerr /= 0) then ! warning if lerr > 0, error if lerr < 0
       lpu = warning_p*max(min(lerr,1),0) + error_p*max(min(-lerr,1),0)
       call reportError (lpu,'Could not insert used CPU time '// &
            &            'into RDB-file '//rdb%fileName, &
            &            addString='closeRDBfile')
       ierr = ierr + min(0,lerr)
    end if

    call closeBinaryDB (rdb%fileNumb,lerr)
    if (lerr < 0) then
       ierr = ierr + lerr
       call reportError (error_p,'Could not close RDB-file '//rdb%fileName, &
            &            addString='closeRDBfile')
    end if

    rdb%fileNumb = -99
    rdb%nBytes = 0_i8
    rdb%empty = .true.

    if (present(haveData)) haveData = .true.

600 format('  RDB-file ',a,' : deleted (empty)')
601 format('  RDB-file ',a,' :',i11,' bytes written')

  end subroutine closeRDBfile

end module rdbModule
