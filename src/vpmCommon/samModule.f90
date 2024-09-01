!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file samModule.f90
!> @brief Assembly of FE matrices into system matrices.

!!==============================================================================
!> @brief This module contains a data structure for assembly of FE matrices.
!>
!> @details The names and meanings of (most of) the data members of the
!> sammodule::samtype structure are adopted from Kolbein Bell's pionering work
!> on the field. See his reports on the SAM library for a thorough elaboration.

module SamModule

  use KindModule, only : dp

  implicit none

  integer, parameter :: version_p = 3 !< Version tag

  !> @brief Data structure for management of FE matrix assembly.
  type SamType

     integer, pointer :: mpar(:)   !< Matrix of parameters
     integer, pointer :: mpmnpc(:) !< Matrix of pointers to MNPCs
     integer, pointer :: mmnpc(:)  !< Matrix of nodal point correspondances
     integer, pointer :: madof(:)  !< Matrix of accumulated DOFs
     integer, pointer :: msc(:)    !< Matrix of status codes
     integer, pointer :: dofType(:)!< Type of DOF (traDof_p, rotDof_p, genDof_p)
     integer, pointer :: mpmceq(:) !< Matrix of pointers to MCEQs
     integer, pointer :: mmceq(:)  !< Matrix of constraint equation definitions
     real(dp),pointer :: ttcc(:)   !< Table of constraint equation coefficients
     integer, pointer :: mpreac(:) !< Matrix of pointers to reaction forces
     integer, pointer :: minex(:)  !< Matrix of internal to external node number
     integer, pointer :: mnnn(:)   !< Matrix of new nodal numbers (optimization)
     integer, pointer :: melcon(:) !< Matrix of element number conversions
     integer, pointer :: meqn(:)   !< Matrix of equation numbers for all DOFs
     integer, pointer :: meqn1(:)  !< Matrix of status 1 equation numbers
     integer, pointer :: meqn2(:)  !< Matrix of status 2 equation numbers
     integer, pointer :: dofPosIn1(:) !< DOF-positions in the 11-matrix
     integer, pointer :: dofPosIn2(:) !< DOF-positions in the 22-matrix

     !! These point to specific elements in mpar
     integer, pointer :: nnod      !< Number of active nodes (NANOD)
     integer, pointer :: nel       !< Number of elements
     integer, pointer :: ndof      !< Number of degrees of freedom
     integer, pointer :: ndof1     !< Number of 1-status dofs
     integer, pointer :: ndof2     !< Number of 2-status dofs
     integer, pointer :: nspdof    !< Number of specified dofs
     integer, pointer :: nceq      !< Number of constraint equations
     integer, pointer :: nsdof     !< Number of supressed dofs
     integer, pointer :: npdof     !< Number of prescribed dofs
     integer, pointer :: nddof     !< Number of dependent dofs
     integer, pointer :: neq       !< Number of equations
     integer, pointer :: nsky      !< Number of elements in the skyline array
     integer, pointer :: nmmnpc    !< Number of elements in mmnpc
     integer, pointer :: nmmceq    !< Number of elements in mmceq

  end type SamType

  !> @brief Standard routine for writing an object to file.
  interface writeObject
     module procedure writeSAM
  end interface

  private :: writeSAM


contains

  !!============================================================================
  !> @brief Nullifies all pointers to enable safe testing after allocation.
  !>
  !> @param[out] samData The sammodule::samtype object to nullify
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 24 Sep 1998

  subroutine nullifySAM (samData)

    type(SamType), intent(out) :: samData

    !! --- Logic section ---

    nullify(samData%mpar)
    nullify(samData%mpmnpc)
    nullify(samData%mmnpc)
    nullify(samData%madof)
    nullify(samData%msc)
    nullify(samData%dofType)
    nullify(samData%mpmceq)
    nullify(samData%mmceq)
    nullify(samData%ttcc)
    nullify(samData%mpreac)
    nullify(samData%minex)
    nullify(samData%mnnn)
    nullify(samData%melcon)
    nullify(samData%meqn)
    nullify(samData%meqn1)
    nullify(samData%meqn2)
    nullify(samData%dofPosIn1)
    nullify(samData%dofPosIn2)

    nullify(samData%nnod)
    nullify(samData%nel)
    nullify(samData%ndof)
    nullify(samData%ndof1)
    nullify(samData%ndof2)
    nullify(samData%nspdof)
    nullify(samData%nceq)
    nullify(samData%nsdof)
    nullify(samData%npdof)
    nullify(samData%nddof)
    nullify(samData%neq)
    nullify(samData%nsky)
    nullify(samData%nmmnpc)
    nullify(samData%nmmceq)

  end subroutine nullifySAM


  !!============================================================================
  !> @brief Deallocates all dynamic arrays within the SAM structure.
  !>
  !> @param samData The sammodule::samtype object to deallocate
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 14 Feb 2003

  subroutine deallocateSAM (samData)

    use AllocationModule, only : reAllocate

    type(SamType), intent(inout) :: samData

    !! --- Logic section ---

    if (associated(samData%minex)) then
       !! Avoid double free if mnnn points to minex
       if (associated(samData%mnnn,samData%minex)) nullify(samData%mnnn)
    end if

    call reAllocate ('deallocateSAM',samData%mpar)
    call reAllocate ('deallocateSAM',samData%mpmnpc)
    call reAllocate ('deallocateSAM',samData%mmnpc)
    call reAllocate ('deallocateSAM',samData%madof)
    call reAllocate ('deallocateSAM',samData%msc)
    call reAllocate ('deallocateSAM',samData%dofType)
    call reAllocate ('deallocateSAM',samData%mpmceq)
    call reAllocate ('deallocateSAM',samData%mmceq)
    call reAllocate ('deallocateSAM',samData%ttcc)
    call reAllocate ('deallocateSAM',samData%mpreac)
    call reAllocate ('deallocateSAM',samData%minex)
    call reAllocate ('deallocateSAM',samData%mnnn)
    call reAllocate ('deallocateSAM',samData%melcon)
    call reAllocate ('deallocateSAM',samData%meqn)
    call reAllocate ('deallocateSAM',samData%meqn1)
    call reAllocate ('deallocateSAM',samData%meqn2)
    call reAllocate ('deallocateSAM',samData%dofPosIn1)
    call reAllocate ('deallocateSAM',samData%dofPosIn2)

  end subroutine deallocateSAM


  !!============================================================================
  !> @brief Performs consistency checking of the SAM data structure.
  !>
  !> @param[in] samData The sammodule::samtype object to check consistency of
  !> @param[in] lpu File unit number for error messages
  !> @param[out] ierr Error flag
  !> @param[in] checkEquations If present and .true., also check the meqn array
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 24 Sep 1998

  subroutine check (samData,lpu,ierr,checkEquations)

    use ReportErrorModule, only : internalError

    type(SamType)   , intent(in)  :: samData
    integer         , intent(in)  :: lpu
    integer         , intent(out) :: ierr
    logical,optional, intent(in)  :: checkEquations

    !! Local variables
    integer :: nel, nnod, ndof, nceq

    !! --- Logic section ---

    ierr = 0
    if (.not. associated(samData%mpar)) then
       ierr = -1
       write(lpu,*) '*** Error: samData%mpar is not allocated'
       goto 999
    else if (size(samData%mpar) < 30) then
       ierr = -1
       write(lpu,*) '*** Error: size(samData%mpar) =',size(samData%mpar)
       goto 999
    end if

    !! Check alias pointers within samData%mpar
    if (.not.associated(samData%nnod,   samData%mpar( 1))) ierr = ierr -1
    if (.not.associated(samData%nel,    samData%mpar( 2))) ierr = ierr -2
    if (.not.associated(samData%ndof,   samData%mpar( 3))) ierr = ierr -4
    if (.not.associated(samData%ndof1,  samData%mpar( 4))) ierr = ierr -8
    if (.not.associated(samData%ndof2,  samData%mpar( 5))) ierr = ierr -16
    if (.not.associated(samData%nspdof, samData%mpar( 6))) ierr = ierr -32
    if (.not.associated(samData%nceq ,  samData%mpar( 7))) ierr = ierr -64
    if (.not.associated(samData%nsdof,  samData%mpar( 8))) ierr = ierr -128
    if (.not.associated(samData%npdof,  samData%mpar( 9))) ierr = ierr -256
    if (.not.associated(samData%nddof,  samData%mpar(10))) ierr = ierr -512
    if (.not.associated(samData%neq,    samData%mpar(11))) ierr = ierr -1024
    if (.not.associated(samData%nsky,   samData%mpar(12))) ierr = ierr -2048
    if (.not.associated(samData%nmmnpc, samData%mpar(15))) ierr = ierr -16384
    if (.not.associated(samData%nmmceq, samData%mpar(16))) ierr = ierr -32768
    if (ierr < 0) then
       write(lpu,*) '*** Warning: Inconsistency in samData, ierr =',ierr
       ierr = 0
    end if

    nnod = samData%mpar(1)
    nel  = samData%mpar(2)
    ndof = samData%mpar(3)
    nceq = samData%mpar(7)

    !! Element connection tables
    if (.not. associated(samData%mpmnpc)) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: samData%mpmnpc is not allocated'
       goto 100
    else if (size(samData%mpmnpc) < nel+1) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: size(samData%mpmnpc) =',size(samData%mpmnpc)
       write(lpu,*) '                          nel+1 =',nel+1
       goto 100
    end if

    if (.not. associated(samData%mmnpc)) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: samData%mmnpc is not allocated'
    else if (size(samData%mmnpc) < samData%mpmnpc(nel+1)-1) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: size(samData%mmnpc) =',size(samData%mmnpc)
       write(lpu,*) '               mpmnpc(nel+1)-1 =',samData%mpmnpc(nel+1)-1
    end if

    !! Dofs pr. node and dof status codes
100 if (.not. associated(samData%madof)) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: samData%madof is not allocated'
    else if (size(samData%madof) < nnod) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: size(samData%madof) =',size(samData%madof)
       write(lpu,*) '                          nnod =',nnod
    end if

    if (.not. associated(samData%msc)) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: samData%msc is not allocated'
    else if (size(samData%msc) < ndof) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: size(samData%msc) =',size(samData%msc)
       write(lpu,*) '                        ndof =',ndof
    end if

    if (.not. associated(samData%dofType)) then
       if (samData%mpar(18) == 0) then ! Only the solver needs this array
          ierr = ierr - 1
          write(lpu,*) '*** Error: samData%dofType is not allocated'
       end if
    else if (size(samData%dofType) < ndof) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: size(samData%dofType) =',size(samData%dofType)
       write(lpu,*) '                            ndof =',ndof
    end if

    !! Constraint equations
    if (.not. associated(samData%mpmceq)) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: samData%mpmceq is not allocated'
       goto 200
    else if (size(samData%mpmceq) < nceq+1) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: size(samData%mpmceq) =',size(samData%mpmceq)
       write(lpu,*) '                         nceq+1 =',nceq+1
       goto 200
    else if (size(samData%mpmceq) <= 1) then
       goto 200
    end if

    if (.not. associated(samData%mmceq)) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: samData%mmceq is not allocated'
    else if (size(samData%mmceq) < samData%mpmceq(nceq+1)-1) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: size(samData%mmceq) =',size(samData%mmceq)
       write(lpu,*) '              mpmceq(nceq+1)-1 =',samData%mpmceq(nceq+1)-1
    end if

    if (.not. associated(samData%ttcc)) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: samData%ttcc is not allocated'
    else if (size(samData%ttcc) < samData%mpmceq(nceq+1)-1) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: size(samData%ttcc) =',size(samData%ttcc)
       write(lpu,*) '             mpmceq(nceq+1)-1 =',samData%mpmceq(nceq+1)-1
    end if

    !! External/internal node numbers + renumbering
200 if (.not. associated(samData%minex)) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: samData%minex is not allocated'
    else if (size(samData%minex) < nnod) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: size(samData%minex) =',size(samData%minex)
       write(lpu,*) '                          nnod =',nnod
    end if

    if (.not. associated(samData%mnnn)) then
       if (samData%mpar(18) /= 0) then ! The solver does not need this array
          ierr = ierr - 1
          write(lpu,*) '*** Error: samData%mnnn is not allocated'
       end if
    else if (size(samData%mnnn) < nnod) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: size(samData%mnnn) =',size(samData%mnnn)
       write(lpu,*) '                         nnod =',nnod
    end if

    !! External/internal element numbers
    if (.not. associated(samData%melcon)) then
       if (samData%mpar(18) /= 0) then ! The solver does not need this array
          ierr = ierr - 1
          write(lpu,*) '*** Error: samData%melcon is not allocated'
       end if
    else if (size(samData%melcon) < nel) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: size(samData%melcon) =',size(samData%melcon)
       write(lpu,*) '                            nel =',nel
    end if

    if (.not. present(checkEquations)) goto 999
    if (.not. checkEquations) goto 999

    !! Equation numbers
    if (.not. associated(samData%meqn)) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: samData%meqn is not allocated'
    else if (size(samData%meqn) < nceq) then
       ierr = ierr - 1
       write(lpu,*) '*** Error: size(samData%meqn) =',size(samData%meqn)
       write(lpu,*) '                         nceq =',nceq
    end if

999 continue
    if (ierr < 0) ierr = internalError('check: Inconsistent SAM data structure')

  end subroutine check


  !!============================================================================
  !> @brief Returns internal node number and local dof-number for a global dof.
  !>
  !> @param[in] sam The sammodule::samtype object to consider
  !> @param idof Global DOF number on input, local DOF number on output
  !> @param extNum If present and .true., return external node number instead
  !> @return Internal node number in range [1,nnod]
  !> @return -1 if @a idof is out of range
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 29 Jun 1999

  function NodeNumFromDofNum (sam,idof,extNum)

    use ReportErrorModule, only : internalError

    type(SamType)   , intent(in)    :: sam
    integer         , intent(inout) :: idof
    logical,optional, intent(in)    :: extNum
    integer                         :: NodeNumFromDofNum

    !! Local variables
    integer           :: inod
    character(len=64) :: errMsg

    !! --- Logic section ---

    do inod = 1, size(sam%madof)-1
       if (idof >= sam%madof(inod) .and. idof < sam%madof(inod+1)) then
          idof = idof+1 - sam%madof(inod) ! Local DOF within this node
          if (.not.associated(sam%minex) .or. .not.present(extNum)) then
             NodeNumFromDofNum = inod
          else if (extNum) then
             NodeNumFromDofNum = sam%minex(inod)
          else
             NodeNumFromDofNum = inod
          end if
          return
       end if
    end do

    write(errMsg,"('NodeNumFromDofNum: Invalid DOF number:',I10)") idof
    NodeNumFromDofNum = internalError(errMsg)

  end function NodeNumFromDofNum


  !!============================================================================
  !> @brief Returns internal node number and local dof-number for an equation.
  !>
  !> @param[in] sam The sammodule::samtype object to consider
  !> @param ieq Equation number on input, local DOF number on output
  !> @param extNum If present and .true., return external node number instead
  !> @return Internal node number in range [1,nnod]
  !> @return -1 if @a ieq is out of range
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 29 Jun 1999

  function NodeNumFromEqNum (sam,ieq,extNum)

    use ReportErrorModule, only : internalError

    type(SamType)   , intent(in)    :: sam
    integer         , intent(inout) :: ieq
    logical,optional, intent(in)    :: extNum
    integer                         :: NodeNumFromEqNum

    !! Local variables
    integer           :: idof
    character(len=64) :: errMsg

    !! --- Logic section ---

    do idof = 1, size(sam%meqn)
       if (ieq == sam%meqn(idof)) then
          ieq = idof
          NodeNumFromEqNum = NodeNumFromDofNum(sam,ieq,extNum)
          return
       end if
    end do

    write(errMsg,"('NodeNumFromEqNum: Invalid equation number:',I10)") ieq
    NodeNumFromEqNum = internalError(errMsg)

  end function NodeNumFromEqNum


  !!============================================================================
  !> @brief Writes out allocated storage for the SAM data structure.
  !>
  !> @param[in] sam The sammodule::samtype object to write for
  !> @param[in] lpu File unit number to print to
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 7 Mar 2003

  subroutine WriteSAMsize (sam,lpu)

    use KindModule      , only : i8
    use AllocationModule, only : writeStorage

    type(SamType), intent(in) :: sam
    integer      , intent(in) :: lpu

    !! Local variables
    integer(i8) :: nwI, nwR

    !! --- Logic section ---

    nwI = 0_i8
    nwR = 0_i8
    if (associated(sam%mpar))    nwI = nwI + size(sam%mpar,kind=i8)
    if (associated(sam%mpmnpc))  nwI = nwI + size(sam%mpmnpc,kind=i8)
    if (associated(sam%mmnpc))   nwI = nwI + size(sam%mmnpc,kind=i8)
    if (associated(sam%madof))   nwI = nwI + size(sam%madof,kind=i8)
    if (associated(sam%msc))     nwI = nwI + size(sam%msc,kind=i8)
    if (associated(sam%dofType)) nwI = nwI + size(sam%dofType,kind=i8)
    if (associated(sam%mpmceq))  nwI = nwI + size(sam%mpmceq,kind=i8)
    if (associated(sam%mmceq))   nwI = nwI + size(sam%mmceq,kind=i8)
    if (associated(sam%ttcc))    nwR = nwR + size(sam%ttcc,kind=i8)
    if (associated(sam%mpreac))  nwI = nwI + size(sam%mpreac,kind=i8)
    if (associated(sam%minex))   nwI = nwI + size(sam%minex,kind=i8)
    if (associated(sam%mnnn) .and. .not.associated(sam%mnnn,sam%minex)) then
       nwI = nwI + size(sam%mnnn,kind=i8)
    end if
    if (associated(sam%melcon))  nwI = nwI + size(sam%melcon,kind=i8)
    if (associated(sam%meqn))    nwI = nwI + size(sam%meqn,kind=i8)
    if (associated(sam%meqn1))   nwI = nwI + size(sam%meqn1,kind=i8)
    if (associated(sam%meqn2))   nwI = nwI + size(sam%meqn2,kind=i8)
    if (associated(sam%dofPosIn1)) nwI = nwI + size(sam%dofPosIn1,kind=i8)
    if (associated(sam%dofPosIn2)) nwI = nwI + size(sam%dofPosIn2,kind=i8)
    call writeStorage ('the SAM data structure',nwI,nwR,lpu)

  end subroutine WriteSAMsize


  !!============================================================================
  !> @brief Standard routine for writing an object to file.
  !>
  !> @param[in] sam The sammodule::samtype object to print
  !> @param[in] io File unit number to write to
  !> @param[in] complexity If present, the value indicates the amount of print
  !> @param[in] title If present, write as heading, otherwise 'SAM' is written
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 17 Oct 2000

  subroutine WriteSAM (sam,io,complexity,title)

    type(SamType)   , intent(in) :: sam
    integer         , intent(in) :: io
    integer,optional, intent(in) :: complexity
    character(len=*), intent(in), optional :: title

    !! Local variables
    integer :: i, j, k, n, start, end

    !! --- Logic section ---

    if (.not. associated(sam%mpar)) return

    if (present(title)) then
       write(io,'(A)') trim(title),'{'
    else
       write(io,'(A)') 'SAM','{'
    end if
    write(io,600) sam%mpar(1:16),size(sam%mpar),sam%mpar(17:)

    if (present(complexity)) then
       if (complexity > 1 .and. complexity /= 5) then

          if (associated(sam%mpmnpc)) then
             n = size(sam%mpmnpc)
             write(io,*)
             write(io,*) 'MPMNPC', sam%mpar(2), n
             write(io,'(2I8)') (i,sam%mpmnpc(i),i=1,n)
             if (associated(sam%mmnpc)) then
                write(io,*)
                write(io,*) 'MMNPC', sam%mpar(15), size(sam%mmnpc)
                do i = 1,sam%mpar(2)
                   j = sam%mpmnpc(i)
                   k = sam%mpmnpc(i+1)-1
                   write(io,'(I8," :",20I8/(10X,20I8))') i, sam%mmnpc(j:k)
                end do
             end if
          end if

          if (associated(sam%madof)) then
             n = size(sam%madof)
             write(io,*)
             write(io,*) 'MADOF', sam%mpar(1), n
             write(io,'(2I8)') (i,sam%madof(i),i=1,n)
             if (associated(sam%msc)) then
                write(io,*)
                write(io,*) 'MSC', sam%mpar(3), size(sam%msc)
                do i = 1,sam%mpar(1)
                   j = sam%madof(i)
                   k = sam%madof(i+1)-1
                   if (any(sam%msc(j:k) < 0)) then
                      write(io,'(I8," :",20I4/(10X,20I4))') i, sam%msc(j:k)
                   else
                      write(io,'(I8," :",20I2/(10X,20I2))') i, sam%msc(j:k)
                   end if
                end do
             end if
             if (associated(sam%dofType)) then
                write(io,*)
                write(io,*) 'dofType', sam%mpar(3), size(sam%dofType)
                do i = 1,sam%mpar(1)
                   j = sam%madof(i)
                   k = sam%madof(i+1)-1
                   write(io,'(I8," :",20I2/(10X,20I2))') i, sam%dofType(j:k)
                end do
             end if
          end if

          if (associated(sam%mpmceq)) then
             n = size(sam%mpmceq)
             write(io,*)
             write(io,*) 'MPMCEQ', sam%mpar(7), n
             write(io,'(2I8)') (i,sam%mpmceq(i),i=1,n)
             if (associated(sam%mmceq) .and. associated(sam%ttcc)) then
                write(io,*)
                write(io,*) 'MMCEQ/TTCC', size(sam%mmceq), size(sam%ttcc)
                do i = 1,sam%mpar(7)
                   j = sam%mpmceq(i)
                   k = sam%mpmceq(i+1)-1
                   end = min(j+11,k)
                   write(io,'(I8," : ",12I8)') i,sam%mmceq(j:end)
                   write(io,'(11X,12F8.3)') sam%ttcc(j:end)
                   do start = j+12,k,12
                      end = min(start+11,k)
                      write(io,'(11X,12I8)') sam%mmceq(start:end)
                      write(io,'(11X,12F8.3)') sam%ttcc(start:end)
                   end do
                end do
             end if
          end if

          if (associated(sam%mpreac)) then
             n = size(sam%mpreac)
             write(io,*)
             write(io,*) 'mpreac', sam%mpar(3), n
             write(io,'(2I8)') (i,sam%mpreac(i),i=1,n)
          end if

          if (associated(sam%minex) .and. associated(sam%mnnn)) then
             n = size(sam%minex)
             write(io,*)
             write(io,*) 'MINEX/MNNN',sam%mpar(1),size(sam%minex),size(sam%mnnn)
             write(io,'(3I8)') (i,sam%minex(i),sam%mnnn(i),i=1,sam%mpar(1))
          else if (associated(sam%minex)) then
             n = size(sam%minex)
             write(io,*)
             write(io,*) 'MINEX',sam%mpar(1),size(sam%minex)
             write(io,'(2I8)') (i,sam%minex(i),i=1,sam%mpar(1))
          end if

          if (associated(sam%melcon)) then
             n = size(sam%melcon)
             write(io,*)
             write(io,*) 'MELCON', sam%mpar(2), n
             write(io,'(2I8)') (i,sam%melcon(i),i=1,n)
          end if

       end if
       if (complexity > 2 .and. complexity /= 5) then

          if (associated(sam%meqn)) then
             n = size(sam%meqn)
             write(io,*)
             write(io,*) 'MEQN', sam%mpar(3), n
             write(io,'(2I8)') (i,sam%meqn(i),i=1,n)
          end if

          if (associated(sam%meqn1)) then
             n = size(sam%meqn1)
             write(io,*)
             write(io,*) 'MEQN1', sam%mpar(4), n
             write(io,'(2I8)') (i,sam%meqn1(i),i=1,n)
          end if

          if (associated(sam%meqn2)) then
             n = size(sam%meqn2)
             write(io,*)
             write(io,*) 'MEQN2', sam%mpar(5), n
             write(io,'(2I8)') (i,sam%meqn2(i),i=1,n)
          end if

          if (associated(sam%dofPosIn1)) then
             n = size(sam%dofPosIn1)
             write(io,*)
             write(io,*) 'dofPosIn1', sam%mpar(3), n
             write(io,'(2I8)') (i,sam%dofPosIn1(i),i=1,n)
          end if

          if (associated(sam%dofPosIn2)) then
             n = size(sam%dofPosIn2)
             write(io,*)
             write(io,*) 'dofPosIn2', sam%mpar(5), n
             write(io,'(2I8)') (i,sam%dofPosIn2(i),i=1,n)
          end if

       end if
       if (complexity == 5) then ! Alternative format for external debugging

          if (associated(sam%mpmnpc)) then
             write(io,*)
             write(io,*) 'MPMNPC:', size(sam%mpmnpc)
             write(io,'(10I8)') sam%mpmnpc
          end if
          if (associated(sam%mmnpc)) then
             write(io,*)
             write(io,*) 'MMNPC:', size(sam%mmnpc)
             write(io,'(10I8)') sam%mmnpc
          end if
          if (associated(sam%madof)) then
             write(io,*)
             write(io,*) 'MADOF:', size(sam%madof)
             write(io,'(10I8)') sam%madof
          end if
          if (associated(sam%msc)) then
             write(io,*)
             write(io,*) 'MSC:', size(sam%msc)
             write(io,'(10I8)') sam%msc
          end if
          if (associated(sam%dofType)) then
             write(io,*)
             write(io,*) 'dofType:', size(sam%dofType)
             write(io,'(10I8)') sam%dofType
          end if
          if (associated(sam%mpmceq)) then
             write(io,*)
             write(io,*) 'MPMCEQ:', size(sam%mpmceq)
             write(io,'(10I8)') sam%mpmceq
          end if
          if (associated(sam%mmceq)) then
             write(io,*)
             write(io,*) 'MMCEQ:', size(sam%mmceq)
             write(io,'(10I8)') sam%mmceq
          end if
          if (associated(sam%ttcc)) then
             write(io,*)
             write(io,*) 'TTCC:', size(sam%ttcc)
             write(io,'(1P6E22.15)') sam%ttcc
          end if
          if (associated(sam%mpreac)) then
             write(io,*)
             write(io,*) 'MPREAC:', size(sam%mpreac)
             write(io,'(10I8)') sam%mpreac
          end if
          if (associated(sam%minex)) then
             write(io,*)
             write(io,*) 'MINEX:',size(sam%minex)
             write(io,'(10I8)') sam%minex
          end if
          if (associated(sam%mnnn)) then
             write(io,*)
             write(io,*) 'MNNN:',size(sam%mnnn)
             write(io,'(10I8)') sam%mnnn
          end if
          if (associated(sam%melcon)) then
             write(io,*)
             write(io,*) 'MELCON:', size(sam%melcon)
             write(io,'(10I8)') sam%melcon
          end if
          if (associated(sam%meqn)) then
             n = size(sam%meqn)
             write(io,*)
             write(io,*) 'MEQN:', size(sam%meqn)
             write(io,'(10I8)') sam%meqn
          end if
          if (associated(sam%meqn1)) then
             write(io,*)
             write(io,*) 'MEQN1:', size(sam%meqn1)
             write(io,'(10I8)') sam%meqn1
          end if
          if (associated(sam%meqn2)) then
             write(io,*)
             write(io,*) 'MEQN2:', size(sam%meqn2)
             write(io,'(10I8)') sam%meqn2
          end if
          if (associated(sam%dofPosIn1)) then
             write(io,*)
             write(io,*) 'dofPosIn1:', size(sam%dofPosIn1)
             write(io,'(10I8)') sam%dofPosIn1
          end if
          if (associated(sam%dofPosIn2)) then
             write(io,*)
             write(io,*) 'dofPosIn2:', size(sam%dofPosIn2)
             write(io,'(10I8)') sam%dofPosIn2
          end if

       end if
    end if
    write(io,'(A/)') '}'

600 format (' mpar(1)  = nnod   =',I8, &
         & /' mpar(2)  = nel    =',I8, &
         & /' mpar(3)  = ndof   =',I8, &
         & /' mpar(4)  = ndof1  =',I8, &
         & /' mpar(5)  = ndof2  =',I8, &
         & /' mpar(6)  = nspdof =',I8, &
         & /' mpar(7)  = nceq   =',I8, &
         & /' mpar(8)  = nsdof  =',I8, &
         & /' mpar(9)  = npdof  =',I8, &
         & /' mpar(10) = nddof  =',I8, &
         & /' mpar(11) = neq    =',I8, &
         & /' mpar(12) = nsky   =',I8, &
         & /' mpar(13) = lbw    =',I8, &
         & /' mpar(14) = nodlbw =',I8, &
         & /' mpar(15) = nmmnpc =',I8, &
         & /' mpar(16) = nmmceq =',I8, &
         & /' mpar(17:',I2,') =',4I8, &
         &/('              ',10I8))

  end subroutine WriteSAM


  !!============================================================================
  !> @brief Duplicates the SAM object to emulate a model with 1-DOF elements.
  !>
  !> @param[in]  orgSam The sammodule::samtype object to be duplicated
  !> @param[out] newSam Pointer to the new sammodule::samtype object
  !> @param[out] ierr   Error flag (negative on allocation error)
  !>
  !> @details All the arrays of the @a orgSam object are shared with the
  !> @a newSam object (that is, the array pointers points to the same physical
  !> address), except for @a madof which is re-initialized with 1 DOF per node.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 10 Sep 2004

  subroutine LumpedSAM (orgSam,newSam,ierr)

    use AllocationModule , only : reAllocate
    use ReportErrorModule, only : allocationError

    type(SamType), intent(in)  :: orgSam
    type(SamType), pointer     :: newSam
    integer      , intent(out) :: ierr

    !! Local variables
    integer :: i

    !! --- Logic section ---

    allocate(newSam,stat=ierr)
    if (ierr /= 0) then
       ierr = allocationError('LumpedSAM')
       return
    end if

    call nullifySam (newSam)

    !! Initialize the Matrix of Parameters
    call reAllocate ('LumpedSAM',newSam%mpar,size(orgSam%mpar),ierr)
    if (ierr < 0) return

    newSam%mpar   =  orgSam%mpar
    newSam%nnod   => newSam%mpar(1)  ! Number of active nodes (nanod)
    newSam%nel    => newSam%mpar(2)  ! Number of elements
    newSam%ndof   => newSam%mpar(3)  ! Number of degrees of freedom
    newSam%ndof1  => newSam%mpar(4)  ! Number of 1-status dofs
    newSam%ndof2  => newSam%mpar(5)  ! Number of 2-status dofs
    newSam%nspdof => newSam%mpar(6)  ! Number of specified dofs
    newSam%nceq   => newSam%mpar(7)  ! Number of constraint equations
    newSam%nsdof  => newSam%mpar(8)  ! Number of supressed dofs
    newSam%npdof  => newSam%mpar(9)  ! Number of prescribed dofs
    newSam%nddof  => newSam%mpar(10) ! Number of dependent dofs
    newSam%neq    => newSam%mpar(11) ! Number of equations
    newSam%nsky   => newSam%mpar(12) ! Number of elements in SKYline array
    newSam%nmmnpc => newSam%mpar(15) ! Number of elements in mmnpc
    newSam%nmmceq => newSam%mpar(16) ! Number of elements in mmceq

    newSam%nnod   = orgSam%ndof
    newSam%nel    = orgSam%ndof
    newSam%nmmnpc = orgSam%ndof

    call reAllocate ('LumpedSAM',newSam%madof,newSam%nnod+1,ierr)
    if (ierr < 0) return

    newSam%madof = (/ (i,i=1,newSam%nnod+1) /)

    newSam%minex  => newSam%madof
    newSam%mpmnpc => newSam%madof
    newSam%mmnpc  => newSam%madof
    newSam%mnnn   => newSam%madof
    newSam%melcon => newSam%madof

    newSam%msc    => orgSam%msc
    newSam%mpmceq => orgSam%mpmceq
    newSam%mmceq  => orgSam%mmceq
    newSam%ttcc   => orgSam%ttcc

  end subroutine LumpedSAM


  !!============================================================================
  !> @brief Forms the control arrays @a dofPosIn1 and/or @a dofPosIn2 in SAM.
  !>
  !> @param lp1 If .true., create the @a dofPosIn1 array
  !> @param lp2 If .true., create the @a dofPosIn2 array
  !> @param sam The sammodule::samtype object to consider
  !> @param[out] ierr Error flag
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 18 May 2020

  subroutine invertEqPartition (lp1,lp2,sam,ierr)

    use ScratchArrayModule, only : getIntegerScratchArray
    use AllocationModule  , only : reAllocate
    use ReportErrorModule , only : internalError

    logical      , intent(in)    :: lp1, lp2
    type(SamType), intent(inout) :: sam
    integer      , intent(out)   :: ierr

    !! Local variables
    integer          :: i1, i2, ieq, idof
    integer, pointer :: ieqdof(:)

    !! --- Logic section ---

    ierr = 0
    ieqdof => sam%meqn ! To avoid silly compiler warning from gfortran

    if (lp1) then
       if (associated(sam%meqn1)) then
          ieqdof => getIntegerScratchArray(sam%neq,ierr)
          call reAllocate ('invertEqPartition',sam%dofPosIn1,sam%ndof,ierr)
       else
          ierr = internalError('invertEqPartition: MEQN1 not defined')
       end if
    end if
    if (lp2) then
       if (associated(sam%meqn2)) then
          call reAllocate ('invertEqPartition',sam%dofPosIn2,sam%ndof2,ierr)
       else
          ierr = ierr + internalError('invertEqPartition: MEQN2 not defined')
       end if
    end if
    if (ierr < 0) return

    if (lp1) then

       !! First, invert the DOF-to-equation mapping, MEQN -> IEQDOF
       call ICOPY (sam%neq,-1,0,ieqdof(1),1)
       do idof = 1, sam%ndof
          ieq  = sam%meqn(idof)
          if (ieq > sam%neq) then
             ierr = internalError('invertEqPartition: Invalid MEQN array')
             return
          else if (ieq > 0) then
             ieqdof(ieq) = idof
          end if
       end do

       !! Then use the temporary IEQDOF array to compute dofPosIn1.
       !! This should be much faster than using the intrinsic function
       !! findloc for large models, as we here avoid the nested do-loop.
       call ICOPY (sam%ndof,-1,0,sam%dofPosIn1(1),1)
       do i1 = 1, sam%ndof1
          ieq = sam%meqn1(i1)
          if (ieq > 0 .and. ieq <= sam%neq) then
             idof = ieqdof(ieq)
             if (idof > 0) then
                sam%dofPosIn1(idof) = i1
                cycle
             end if
          end if
          ierr = internalError('invertEqPartition: Invalid MEQN1 array')
          return
       end do

    end if
    if (lp2) then

       !! Since ndof2 usually is much smaller than ndof and ndof1,
       !! we just use the findloc intrinsic here to invert the MEQN2 array
       i2 = 0
       do idof = 1, sam%ndof
          if (sam%msc(idof) == 2) then
             i2 = i2 + 1
             sam%dofPosIn2(i2) = findloc(sam%meqn2,sam%meqn(idof),1)
          end if
       end do

    end if

#if defined(__GNUC__) && (__GNUC__ < 9)
    !! The findloc intrinsic is a Fortran 2008 feature.
    !! Not available in older gfortran implementations.
    !! Remove this when migrating to gfortran 9.0
  contains
    !> @cond NO_DOCUMENTATION
    function findloc (array,val,dim) result(i)
      integer, intent(in) :: array(:), val, dim
      integer :: i
      do i = dim, size(array)
         if (array(i) == val) return
      end do
      i = 0
    end function findloc
    !> @endcond
#endif

  end subroutine invertEqPartition

end module SamModule
