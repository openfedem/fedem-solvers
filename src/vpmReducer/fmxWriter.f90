!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file fmxWriter.f90
!> @brief Shared library wrapper for reading and writing FEDEM fmx files.
!> @author Knut Morten Okstad, SAP SE
!> @date 22 Nov 2019

!> @brief Module encapsulation of the IO_FMX subroutine.

module FMXwriter

  implicit none

contains

  !> @brief Read/writes a double precision array from/to the specified fmx-file.
  !> @param[in] prefix File name prefix
  !> @param[in] itype Matrix type (1=stiffness, 2=mass, 3=gravity forces)
  !> @param data Array with matrix content
  !> @param[in] nval Size of the data array
  !> @param[in] chkSum Write to file if present, otherwise read from file
  !> @param[out] ierr Error flag (negative value indicates error, otherwise OK)

  subroutine io_FMX (prefix,itype,data,nval,chkSum,ierr)

    use kindModule       , only : dp, lfnam_p
    use binaryDBInterface, only : writeDoubleDB, readDoubleDB

    character(len=*) , intent(in)    :: prefix
    integer          , intent(in)    :: itype, nval
    real(dp)         , intent(inout) :: data(*)
    integer, optional, intent(in)    :: chkSum
    integer          , intent(out)   :: ierr

    character(len=lfnam_p) :: fileName
    character(len=32)      :: fileTag

    select case (itype)
    case (1)
       fileName = trim(prefix)//'_S.fmx'
       fileTag  = 'stiffness matrix'
    case (2)
       fileName = trim(prefix)//'_M.fmx'
       fileTag  = 'mass matrix'
    case (3)
       fileName = trim(prefix)//'_G.fmx'
       fileTag  = 'gravity force vectors'
    case (4)
       fileName = trim(prefix)//'_B.fmx'
       fileTag  = 'disk matrix'
    case default
       print *,'*** FMX: Unknown file type',itype
       ierr = -max(1,itype)
       return
    end select

    ierr = 0
    if (present(chkSum)) then
       print *,'  * Writing '//trim(fileTag)//' file: '//trim(fileName),nval
       call writeDoubleDB (fileName,fileTag,chkSum, &
            &              reshape(data(1:nval),(/nval,1/)),ierr)
    else
       print *,'  * Reading '//trim(fileTag)//' file: '//trim(fileName),nval
       call readDoubleDB (fileName,fileTag,nval,data(1),ierr)
    end if
    if (ierr < 0) print *,'*** Failure:',ierr

  end subroutine io_FMX

end module FMXwriter


!> @brief Writes a double precision array to the specified fmx-file.
!> @param[in] prefix File name prefix
!> @param[in] itype Matrix type (1=stiffness, 2=mass, 3=gravity forces)
!> @param[in] data Array with matrix content
!> @param[in] nval Size of the data array
!> @return Error status. Negative value indicates an error, otherwise OK

function writeFMX (prefix,itype,data,nval) result(ierr)
!DEC$ ATTRIBUTES DLLEXPORT :: WRITEFMX
  use FMXwriter, only : io_FMX
  character(len=*), intent(in)    :: prefix
  integer         , intent(in)    :: itype, nval
  double precision, intent(inout) :: data(*)
  integer         , parameter     :: chkSum = 0
  integer                         :: ierr
  call io_FMX (prefix,itype,data,nval,chkSum,ierr)
end function writeFMX


!> @brief Reads a double precision array from the specified fmx-file.
!> @param[in] prefix File name prefix
!> @param[in] itype Matrix type (1=stiffness, 2=mass, 3=gravity forces)
!> @param[out] data Array with matrix content
!> @param[in] nval Size of the data array
!> @return Error status. Negative value indicates an error, otherwise OK

function readFMX (prefix,itype,data,nval) result(ierr)
!DEC$ ATTRIBUTES DLLEXPORT :: READFMX
  use FMXwriter, only : io_FMX
  character(len=*), intent(in)    :: prefix
  integer         , intent(in)    :: itype, nval
  double precision, intent(inout) :: data(*)
  integer                         :: ierr
  call io_FMX (prefix,itype,data,nval,ierr=ierr)
end function readFMX


!> @brief Initializes the Fortran file unit to be used for error messages.

subroutine initFMX ()
!DEC$ ATTRIBUTES DLLEXPORT :: INITFMX
  use reportErrorModule, only : setErrorFile
  call setErrorFile (6)
end subroutine initFMX


!> @brief Reads some SAM data arrays from the specified fsm-file.
!> @param[in] prefix File name prefix
!> @param ndof Size of the @a meqn array
!> @param ndof1 Size of the @a meqn1 array
!> @param ndof2 Size of the @a meqn2 array
!> @param nceq Size of the @ mpmceq array
!> @param nmmceq Size of the @a mmceq and @a ttcc arrays
!> @param[out] meqn Matrix of equation numbers for all DOFs
!> @param[out] meqn1 Matrix of status 1 equation numbers
!> @param[out] meqn2 Matrix of status 2 equation numbers
!> @param[out] mmceq Matrix of matrices of constraint equations
!> @param[out] mpmceq Matrix of pointers to matrices of constraint equations
!> @param[out] ttcc Table of tables of constraint equation coefficients
!> @return Error status. Negative value indicates an error, otherwise OK
!>
!> @details This function should be invoked twice. In the first call,
!> @a prefix should point to the file to read from and the size parameters
!> @a ndof, @a ndof1 ... @a nmmceq are extracted. The file remains open.
!> Before the second call, the arrays @a meqn, @a meqn1 ... @a ttcc need to be
!> allocated to the lenghts returned by the first call. These arrays are then
!> read in the second call of this function (with a blank @a prefix argument).

function readFSM (prefix,ndof,ndof1,ndof2,nceq,nmmceq, &
     &            meqn,meqn1,meqn2,mmceq,mpmceq,ttcc) result(ierr)
  !DEC$ ATTRIBUTES DLLEXPORT :: READFSM
  use kindModule       , only : dp, lfnam_p
  use binaryDBInterface, only : openBinaryDB, closeBinaryDB, read_p
  use binaryDBInterface, only : readTagDB, readIntDB, readDoubleDB
  use reportErrorModule, only : reportError, error_p

  implicit none

  character(len=*), intent(in)    :: prefix
  integer         , intent(inout) :: ndof, ndof1, ndof2, nceq, nmmceq
  integer         , intent(out)   :: meqn(*), meqn1(*), meqn2(*)
  integer         , intent(out)   :: mmceq(*), mpmceq(*)
  real(dp)        , intent(out)   :: ttcc(*)

  character(len=lfnam_p) :: fileName
  character(len=32)      :: fileTag
  integer, save          :: nnod, nel, nmmnpc, ifile = -1
  integer, allocatable   :: itmp(:)
  integer :: ierr, npar

  if (ifile < 0 .and. len(prefix) > 0) then
     fileName = trim(prefix)//'_SAM.fsm'
     call openBinaryDB (fileName,read_p,ifile,ierr)
     if (ierr < 0) then
        call reportError (error_p,'Error opening SAM-file '//fileName)
        return
     end if

     call readTagDB (ifile,fileTag,npar,ierr)
     if (ierr < 0) then
        call reportError (error_p,'Error reading file tag from '//fileName)
        return
     else if (fileTag /= '#SAM data') then
        call reportError (error_p,'File "'//trim(fileName)// &
             &            '" is not a SAM data file, tag='//fileTag)
        ierr = -2
        return
     else
        print *,'  * Reading SAM data file: '//trim(fileName)
     end if

     !! Initialize the Matrix of Parameters
     call readIntDB (ifile,npar,1,ierr)
     if (ierr < 0) then
        call reportError (error_p,'Error reading MPAR size from '//fileName)
        return
     end if

     allocate(itmp(npar),stat=ierr)
     call readIntDB (ifile,itmp(1),npar,ierr)
     if (ierr < 0) then
        call reportError (error_p,'Error reading MPAR from '//fileName)
        return
     end if

     !! Extract the required array size parameters
     nnod   = itmp(1)
     nel    = itmp(2)
     ndof   = itmp(3)
     ndof1  = itmp(4)
     ndof2  = itmp(5)
     nceq   = itmp(7)
     nmmnpc = itmp(15)
     nmmceq = itmp(16)
     deallocate(itmp)
     return
  else if (ifile < 0) then
     ierr = ifile
     return
  end if

  !! Allocate temporary array for reading beyond the arrays not needed
  npar = max(ndof,nnod+1,nel+1,nmmnpc)
  allocate(itmp(npar),stat=ierr)

  !! Read until we get to the arrays that we actually want
  if (ierr == 0) call readIntDB (ifile,itmp(1),nnod+1,ierr)
  if (ierr >= 0) call readIntDB (ifile,itmp(1),nnod,ierr)
  if (ierr >= 0) call readIntDB (ifile,itmp(1),nnod,ierr)
  if (ierr >= 0) call readIntDB (ifile,itmp(1),ndof,ierr)
  if (ierr >= 0) call readIntDB (ifile,itmp(1),nel+1,ierr)
  if (ierr >= 0) call readIntDB (ifile,itmp(1),nmmnpc,ierr)
  if (ierr >= 0) call readIntDB (ifile,itmp(1),nel,ierr)
  deallocate(itmp)

  !! Here we are...
  if (nceq > 0) then
     if (ierr >= 0) call readIntDB (ifile,mpmceq(1),nceq+1,ierr)
     if (ierr >= 0) call readIntDB (ifile,mmceq(1),nmmceq,ierr)
     if (ierr >= 0) call readDoubleDB (ifile,ttcc(1),nmmceq,ierr)
  end if
  if (ierr >= 0) call readIntDB (ifile,meqn(1),ndof,ierr)
  if (ierr >= 0) call readIntDB (ifile,meqn1(1),ndof1,ierr)
  if (ierr >= 0) call readIntDB (ifile,meqn2(1),ndof2,ierr)

  if (ierr < 0) then
     call reportError (error_p,'Error reading from SAM-file '//fileName)
     return
  end if

  call closeBinaryDB (ifile,ierr)
  if (ierr < 0) call reportError (error_p,'Error closing SAM-file '//fileName)

  ifile = -1

end function readFSM
