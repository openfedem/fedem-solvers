!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

subroutine solveSupEl (ndim,nfix,neq,g,ierr)

  !!============================================================================
  !! Test program for verifying the solution using superelements.
  !!
  !! Programmer : Knut Morten Okstad
  !! date/rev   : 23 Feb 2018/1.0
  !!============================================================================

  use kindModule         , only : dp, lfnam_p
  use fileUtilitiesModule, only : getFileName
  use reportErrorModule  , only : setErrorFile

  implicit none

  integer , intent(in)   :: ndim, nfix, neq
  real(dp), intent(in)   :: g(3)
  integer , intent(out)  :: ierr
  integer , parameter    :: lpu = 10
  character(len=lfnam_p) :: chName

  call getFileName ('resfile',chname,'.res')
  open(lpu,FILE=chname,STATUS='UNKNOWN',POSITION='APPEND',IOSTAT=ierr)
  if (ierr /= 0) then
     write(*,*) '*** Unable to open res-file '//trim(chname)
     return
  end if

  call setErrorFile (lpu)
  write(lpu,"(//'System level solve verification')")

  if (neq > ndim) then
     call solveSupEl1 (ndim,nfix,neq,g,lpu,ierr)
  else
     call solveSupEl2 (ndim,nfix,g,lpu,ierr)
  end if

  close(lpu)

end subroutine solveSupEl


subroutine solveSupEl1 (ndim,nfix,neq,g,lpu,ierr)

  !!============================================================================
  !! Test program for verifying the solution using superelements.
  !!
  !! Programmer : Knut Morten Okstad
  !! date/rev   : 23 Feb 2018/1.0
  !!============================================================================

  use kindModule       , only : dp
  use binaryDBInterface, only : readDoubleDB
  use manipMatrixModule, only : writeObject

  implicit none

  integer , intent(in)  :: ndim, nfix, neq, lpu
  real(dp), intent(in)  :: g(3)
  integer , intent(out) :: ierr

  !! Rely on stack allocation for the local arrays
  integer  :: ipiv(ndim-nfix)
  real(dp) :: Kmat(ndim,ndim), Fg(ndim,3)
  real(dp) :: A(ndim-nfix,ndim-nfix), R(ndim-nfix)
  real(dp) :: B(neq-ndim,ndim), vii(neq-ndim), vei(neq-ndim), vgi(neq-ndim,3)

  character(len=11), parameter :: file1 = 'fedem_S.fmx'
  character(len=11), parameter :: file2 = 'fedem_G.fmx'
  character(len=11), parameter :: file3 = 'fedem_V.fmx'
  character(len=11), parameter :: file4 = 'fedem_B.fmx'

  !! Read reduction matrices from files
  ierr = 0
  call readDoubleDB (file1,'stiffness matrix',ndim*ndim,Kmat(1,1),ierr)
  call readDoubleDB (file2,'gravity force vectors',ndim*3,Fg(1,1),ierr)
  call readDoubleDB (file3,'displacement matrix',(neq-ndim)*3,vgi(1,1),ierr)
  call readDoubleDB (file4,'disk matrix',(neq-ndim)*ndim,B(1,1),ierr)
  if (ierr < 0) then
     write(*,*) '*** Failure reading matrix files.'
     return
  end if

  call writeObject (Kmat,lpu,'Kmat')
  call writeObject (Fg,lpu,'Fg')
  call writeObject (vgi,lpu,'v_gi')
  call writeObject (B,lpu,'Bmat')

  !! Establish the top-level equation system [A]*{v} = {R} where {R}
  !! consists of the static gravitation load, {R} = [Fg]*{gx,gy,gz}.
  !! Enforce boundary conditions by eliminating the first <nfix> DOFs.
  A = Kmat(nfix+1:,nfix+1:)
  R = matmul(Fg(nfix+1:,:),g)
  write(lpu,6) 'Right-hand-side vector',R

  !! Solve for displacements, {v} = [A]-1 * {R}
  call DGESV (ndim-nfix,1,A(1,1),ndim-nfix,ipiv(1),R(1),ndim-nfix,ierr)
  if (ierr /= 0) write(lpu,*) '*** DGESV returned INFO =',ierr
  write(lpu,6) 'Solution vector, v_e',R
6 format(/A,':'/(1P6E13.5))

  !! Calculate displacements for the internal DOFs
  vei = matmul(B(:,nfix+1:),R)
  vii = matmul(vgi,g)

  write(lpu,6) 'v^e_i',vei
  write(lpu,6) 'v^i_i',vii
  write(lpu,6) 'v_i = v^e_i + v^i_i',vii + vei

end subroutine solveSupEl1


subroutine solveSupEl2 (ndim,nfix,g,lpu,ierr)

  !!============================================================================
  !! Test program for verifying the solution using superelements.
  !!
  !! Programmer : Knut Morten Okstad
  !! date/rev   : 6 Oct 2022/1.0
  !!============================================================================

  use kindModule            , only : dp, lfnam_p
  use fileUtilitiesModule   , only : getFileName
  use binaryDBInterface     , only : readDoubleDB
  use manipMatrixModule     , only : writeObject
  use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getdoubles

  implicit none

  integer , intent(in)  :: ndim, nfix, lpu
  real(dp), intent(in)  :: g(3)
  integer , intent(out) :: ierr

  !! Rely on stack allocation for the local arrays
  integer  :: ipiv(ndim-nfix), i
  real(dp) :: Kmat(ndim,ndim), Fg(ndim,3), S(6)
  real(dp) :: A(ndim-nfix,ndim-nfix), R(ndim-nfix)

  character(len=lfnam_p) :: file1, file2
  call getFileName ('stiffile',file1,'_S.fmx')
  call getFileName ('gravfile',file2,'_G.fmx')

  !! Read reduction matrices from files
  ierr = 0
  call readDoubleDB (file1,'stiffness matrix',ndim*ndim,Kmat(1,1),ierr)
  call readDoubleDB (file2,'gravity force vectors',ndim*3,Fg(1,1),ierr)
  if (ierr < 0) then
     write(*,*) '*** Failure reading matrix files.'
     return
  end if

  call writeObject (Kmat,lpu,'Kmat')
  call writeObject (Fg,lpu,'Fg')

  !! Establish the top-level equation system [A]*{v} = {R} where {R}
  !! consists of the static gravitation load, {R} = [Fg]*{gx,gy,gz}.
  !! Enforce boundary conditions by eliminating the first <nfix> DOFs.
  A = Kmat(nfix+1:,nfix+1:)
  R = matmul(Fg(nfix+1:,:),g)
  write(lpu,6) 'Right-hand-side vector',R

  !! Solve for displacements, {v} = [A]-1 * {R}
  call DGESV (ndim-nfix,1,A(1,1),ndim-nfix,ipiv(1),R(1),ndim-nfix,ierr)
  if (ierr /= 0) write(lpu,*) '*** DGESV returned INFO =',ierr
  write(lpu,6) 'Solution vector, v_e',R
6 format(/A,':'/(1P6E13.5))

  !! Verify the solution
  call ffa_cmdlinearg_getdoubles ('solution',S,6)
  do i = 1, 6
     if (abs(S(i)) > 0.0_dp .and. abs(R(i)-S(i)) > 1.0e-9_dp) then
        write(*,*) '*** Discrepancy in solution',i,R(i),R(i)-S(i)
        ierr = ierr + 1
     end if
  end do

end subroutine solveSupEl2
