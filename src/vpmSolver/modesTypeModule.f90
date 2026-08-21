!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file modesTypeModule.f90
!> @brief Eigenmodes data container.

!!==============================================================================
!> @brief Module with a data type representing the eigenmodes of the model.
!>
!> @details The module also contains subroutines for writing modal results
!> to the results database, as well as general subroutines and functions for
!> for accessing the modal data.

module ModesTypeModule

  use KindModule         , only : dp
  use SysMatrixTypeModule, only : SysMatrixType

  implicit none

  !> @brief Data type containing parameters and matrices for modal analysis.
  type ModesType

     integer  :: solver !< Which eigensolver to use
     integer  :: nModes !< Number of eigenmodes to calculate
     real(dp) :: shift  !< Shift value
     real(dp) :: tol(3) !< Tolerances for factorization and eigenvalues/vectors
     logical  :: addBC  !< .true. if the additional BCs should be applied
     logical  :: factorMass !< .true. if factorization of mass matrix in LANCZ2
     logical  :: stressStiffIsOn !< .true. if geometric stiffness should be used

     type(SysMatrixType) :: Kmat !< System stiffness matrix
     type(SysMatrixType) :: Cmat !< System damping matrix
     type(SysMatrixType) :: Mmat !< System mass matrix
     type(SysMatrixType) :: Qmat !< System steady-state error elimination matrix

     real(dp), pointer :: ReVal(:)   !< Real part of computed eigenvalues
     real(dp), pointer :: ImVal(:)   !< Imaginary part of computed eigenvalues
     real(dp), pointer :: ReVec(:,:) !< Real part of computed eigenvectors
     real(dp), pointer :: ImVec(:,:) !< Imaginary part of computed eigenvectors
     real(dp), pointer :: eqVec(:,:) !< Eigenvectors (real) in equation order
     real(dp), pointer :: dampRat(:) !< Computed damping ratios

     real(dp), pointer :: effMass(:,:) !< Effective modal masses

  end type ModesType


contains

  !!============================================================================
  !> @brief Initializes the ModesType object.
  !>
  !> @param[out] modes Eigenmode data
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 7 Jun 2002

  subroutine nullifyModes (modes)

    use SysMatrixTypeModule, only : nullifySysMatrix

    type(ModesType), intent(out) :: modes

    !! --- Logic section ---

    modes%solver = 0
    modes%nModes = 0
    modes%shift  = 0.0_dp
    modes%tol    = 0.0_dp
    modes%addBC  = .false.
    modes%factorMass = .false.
    modes%stressStiffIsOn = .false.

    nullify(modes%ReVal)
    nullify(modes%ImVal)
    nullify(modes%ReVec)
    nullify(modes%ImVec)
    nullify(modes%eqVec)
    nullify(modes%dampRat)
    nullify(modes%effMass)

    call nullifySysMatrix (modes%Kmat)
    call nullifySysMatrix (modes%Cmat)
    call nullifySysMatrix (modes%Mmat)
    call nullifySysMatrix (modes%Qmat)

  end subroutine nullifyModes


  !!============================================================================
  !> @brief Deallocates the ModesType object.
  !>
  !> @param modes Eigenmode data
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Jan 2017

  subroutine deallocateModes (modes)

    use SysMatrixTypeModule, only : deallocateSysMatrix

    type(ModesType), intent(inout) :: modes

    !! --- Logic section ---

    if (associated(modes%ReVal))   deallocate(modes%ReVal)
    if (associated(modes%ImVal))   deallocate(modes%ImVal)
    if (associated(modes%ReVec))   deallocate(modes%ReVec)
    if (associated(modes%ImVec))   deallocate(modes%ImVec)
    if (associated(modes%eqVec))   deallocate(modes%eqVec)
    if (associated(modes%dampRat)) deallocate(modes%dampRat)
    if (associated(modes%effMass)) deallocate(modes%effMass)

    call deallocateSysMatrix (modes%Kmat)
    call deallocateSysMatrix (modes%Cmat)
    call deallocateSysMatrix (modes%Mmat)

    call nullifyModes (modes)

  end subroutine deallocateModes


  !!============================================================================
  !> @brief Writes results database headers for the modal results.
  !>
  !> @param[in] modes Eigenmode data
  !> @param[in] mechId Id of the mechanism object
  !> @param[in] triads All triads in the model
  !> @param[in] sups All superelements in the model
  !> @param rdb Results database file for control system data
  !> @param[in] nBits Size of the real words of the results data (32 or 64)
  !> @param ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date January 2001

  subroutine writeModesHeader (modes,mechId,triads,sups,rdb,nBits,ierr)

    use IdTypeModule     , only : IdType, writeIdHeader
    use TriadTypeModule  , only : TriadType
    use SupElTypeModule  , only : SupElType
    use RDBModule        , only : RDBType, iitem, ivard, idatd, nbs_p
    use reportErrorModule, only : internalError

    type(ModesType) , intent(in)    :: modes
    type(IdType)    , intent(in)    :: mechId
    type(TriadType) , intent(in)    :: triads(:)
    type(SupElType) , intent(in)    :: sups(:)
    type(RDBType)   , intent(inout) :: rdb
    integer         , intent(in)    :: nBits
    integer,optional, intent(inout) :: ierr

    !! Local variables
    integer :: i, j, idReEVal, idImEVal, idEVal, idFreq, idDamp
    integer :: idTra, idRot, id3dof, id6dof
    logical :: complexModes


    !! --- Logic section ---

    complexModes = modes%solver == 4 .or. modes%solver == 6

    if (present(ierr)) then
       ! Check that the variable- and item group reference numbers will not be
       ! more than a one-digit number. If they will, we get a format overflow.
       i = rdb%nVar
       j = rdb%nIG
       if (j > 8 .or. i > 5 .or. (complexModes .and. i > 4)) then
          ierr = ierr + internalError('writeModesHeader: format overflow')
          return
       end if
    end if

    rdb%nIG = rdb%nIG + 1
    idEVal  = rdb%nIG
    write(iitem,602) idEVal,'Eigenvalues'

    if (complexModes) then

       ! Write headers for complex eigenvalues
       rdb%nVar = rdb%nVar + 1
       idReEVal = rdb%nVar
       write(ivard,600) idReEVal,'Re','ANGLE/TIME',nBits
       rdb%nVar = rdb%nVar + 1
       idImEVal = rdb%nVar
       write(ivard,600) idImEVal,'Im','ANGLE/TIME',nBits
       rdb%nVar = rdb%nVar + 1
       idFreq   = rdb%nVar
       write(ivard,600) idFreq,'Eigenfrequency','NONE/TIME',8*nbs_p
       rdb%nVar = rdb%nVar + 1
       idDamp   = rdb%nVar
       write(ivard,600) idDamp,'Damping ratio','NONE',8*nbs_p
       do j = 1, size(modes%ReVal)
          write(iitem,603) j,'<',idReEVal,'><',idImEVal, &
               &           '><',idFreq,'><',idDamp,'>]'
       end do
       rdb%nBytes = rdb%nBytes + size(modes%ReVal)*(2*nBits/8 + 2*nbs_p)

    else

       ! Write headers for real eigenvalues
       rdb%nVar = rdb%nVar + 1
       idReEVal = rdb%nVar
       write(ivard,600) idReEVal,'Eigenvalue','ANGLE/TIME',nBits
       rdb%nVar = rdb%nVar + 1
       idFreq   = rdb%nVar
       write(ivard,600) idFreq,'Eigenfrequency','NONE/TIME',8*nbs_p
       rdb%nVar = rdb%nVar + 1
       idDamp   = rdb%nVar
       write(ivard,600) idDamp,'Damping ratio','NONE',8*nbs_p
       do j = 1, size(modes%ReVal)
          write(iitem,603) j,'<',idReEVal,'><',idFreq,'><',idDamp,'>]'
       end do
       rdb%nBytes = rdb%nBytes + size(modes%ReVal)*(nBits/8 + 2*nbs_p)

    end if

    write(iitem,"(']')")
    call writeIdHeader ('Mechanism',mechId,idatd)
    write(idatd,606) idEVal

    ! Write headers for the eigenvectors
    rdb%nVar = rdb%nVar + 1
    idTra    = rdb%nVar
    write(ivard,601) idTra,'Translational deformation','LENGTH',nBits, &
         &           'VEC3','"d_x","d_y","d_z"'
    rdb%nVar = rdb%nVar + 1
    idRot    = rdb%nVar
    write(ivard,601) idRot,'Angular deformation','ANGLE',nBits, &
         &           'ROT3','"theta_x","theta_y","theta_z"'

    id3dof = 0
    id6dof = 0
    do i = 1, size(triads)
       if (triads(i)%id%baseId < 0) then ! Dummy CG triad for generic parts
          cycle
       else if (triads(i)%nDOFS == 3) then

          ! 3-dof triads, translations only
          if (id3dof == 0) call write3DofGroup (ierr)
          call writeIdHeader ('Triad',triads(i)%id,idatd)
          write(idatd,606) id3dof

          if (complexModes) then
             rdb%nBytes = rdb%nBytes + size(modes%ReVal)*6*nBits/8
          else
             rdb%nBytes = rdb%nBytes + size(modes%ReVal)*3*nBits/8
          end if

       else if (triads(i)%nDOFS == 6) then

          ! 6-dof triads, translations and rotation angles
          if (id6dof == 0) call write6DofGroup (ierr)
          call writeIdHeader ('Triad',triads(i)%id,idatd)
          write(idatd,606) id6dof

          if (complexModes) then
             rdb%nBytes = rdb%nBytes + size(modes%ReVal)*12*nBits/8
          else
             rdb%nBytes = rdb%nBytes + size(modes%ReVal)*6*nBits/8
          end if

       end if
    end do

    do i = 1, size(sups)
       if (associated(sups(i)%genDOFs)) then

          ! Superelement with generalized DOFs, no predefined item groups
          ! because the number of generalized DOFs may vary from supel to supel.
          call writeIdHeader ('Part',sups(i)%id,idatd)
          write(idatd,"(a)") '[;"Eigenvectors";'
          do j = 1, size(modes%ReVal)
             if (complexModes) then
                write(idatd,605) j,nBits,sups(i)%genDOFs%nDOFs, &
                     &             nBits,sups(i)%genDOFs%nDOFs
             else
                write(idatd,604) j,nBits,sups(i)%genDOFs%nDOFs
             end if
          end do
          write(idatd,"(']}')")

          if (complexModes) then
             rdb%nBytes = rdb%nBytes + size(modes%ReVal) &
                  &                  * sups(i)%genDOFs%nDOFs*2*nBits/8
          else
             rdb%nBytes = rdb%nBytes + size(modes%ReVal) &
                  &                  * sups(i)%genDOFs%nDOFs*nBits/8
          end if

       end if
    end do

600 format('<',i1,';"',a,'";',a,';FLOAT;',i2,';SCALAR>')
601 format('<',i1,';"',a,'";',a,';FLOAT;',i2,';',a,';(3);((',a,'))>')
602 format('[',i1,';"',a,'";')
603 format('  [;"Mode',i3,'";',a,4(i1,a))
604 format('  [;"Mode',i3,'";', &
         & '<;"Generalized DOFs";LENGTH;FLOAT;',i2,';VECTOR;(',i3,')>]')
605 format('  [;"Mode',i3,'";[;"Re";', &
         & '<;"Generalized DOFs";LENGTH;FLOAT;',i2,';VECTOR;(',i3,')>]' &
         / '              [;"Im";', &
         & '<;"Generalized DOFs";LENGTH;FLOAT;',i2,';VECTOR;(',i3,')>]]')
606 format('[',i1,']}')

  contains

    !> @brief Writes item group definition for a 3-DOF eigenvector.
    subroutine write3DofGroup(ierr)
      integer,optional,intent(inout) :: ierr

      if (present(ierr) .and. rdb%nIG > 8) then
         ierr = ierr + internalError('writeModesHeader: format overflow')
         return
      end if

      rdb%nIG = rdb%nIG + 1
      id3dof  = rdb%nIG
      write(iitem,602) id3dof,'Eigenvectors'
      do j = 1, size(modes%ReVal)
         if (complexModes) then
            write(iitem,603) j,'[;"Re";<',idTra,'>][;"Im";<',idTra,'>]'
         else
            write(iitem,603) j,'<',idTra,'>]'
         end if
      end do
      write(iitem,"(']')")

602   format('[',i1,';"',a,'";')
603   format('  [;"Mode',i3,'";',a,i1,a,i1,a)

    end subroutine write3DofGroup

    !> @brief Writes item group definition for a 6-DOF eigenvector.
    subroutine write6DofGroup(ierr)
      integer,optional,intent(inout) :: ierr

      if (present(ierr)) then
         if (rdb%nIG > 8 .or. (complexModes .and. rdb%nIG > 6)) then
            ierr = ierr + internalError('writeModesHeader: format overflow')
            return
         end if
      end if

      if (complexModes) then
         rdb%nIG = rdb%nIG + 1
         write(iitem,602) rdb%nIG,'Re','<',idTra,'><',idRot,'>]'
         rdb%nIG = rdb%nIG + 1
         write(iitem,602) rdb%nIG,'Im','<',idTra,'><',idRot,'>]'
      end if
      rdb%nIG = rdb%nIG + 1
      id6dof  = rdb%nIG
      write(iitem,602) id6dof,'Eigenvectors'
      do j = 1, size(modes%ReVal)
         if (complexModes) then
            write(iitem,603) j,'[',id6dof-2,'][',id6dof-1,']]'
         else
            write(iitem,603) j,'<',idTra,'><',idRot,'>]'
         end if
      end do
      write(iitem,"(']')")

602   format('[',i1,';"',a,'";',a,i1,a,i1,a)
603   format('  [;"Mode',i3,'";',a,i1,a,i1,a)

    end subroutine write6DofGroup

  end subroutine writeModesHeader


  !!============================================================================
  !> @brief Writes modal results to database file.
  !>
  !> @param[in] modes Eigenmode data
  !> @param[in] sys System level model data
  !> @param[in] triads All triads in the model
  !> @param[in] sups All superelements in the model
  !> @param rdb Results database file for control system data
  !> @param[in] writeAsDouble If .true., save real words as double precision
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date January 2001

  subroutine writeModesDB (modes,sys,triads,sups,rdb,writeAsDouble,ierr)

    use KindModule            , only : pi_p
    use SystemTypeModule      , only : SystemType
    use TriadTypeModule       , only : TriadType
    use SupElTypeModule       , only : SupElType
    use RDBModule             , only : RDBType, writeRDB, writeTimeStepDB
    use RDBModule             , only : checkStepSize
    use reportErrorModule     , only : reportError
    use reportErrorModule     , only : empty_p, warning_p, debugFileOnly_p
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_getint

    type(ModesType) , intent(in)    :: modes
    type(SystemType), intent(in)    :: sys
    type(TriadType) , intent(in)    :: triads(:)
    type(SupElType) , intent(in)    :: sups(:)
    type(RDBType)   , intent(inout) :: rdb
    logical         , intent(in)    :: writeAsDouble
    integer         , intent(out)   :: ierr

    !! Local variables
    integer           :: i, j, k, m, nModes, nTriad, nSupel, nBytesPerStep = 0
    logical           :: complexModes
    real(dp)          :: freq
    character(len=64) :: errMsg

    !! --- Logic section ---

    complexModes = modes%solver == 4 .or. modes%solver == 6

    if (nBytesPerStep == 0) nBytesPerStep = -int(rdb%nBytes)

    nModes = size(modes%ReVal)
    nTriad = size(triads)
    nSupEl = size(sups)

    !! Check if the requested number of eigenmodes were found
    call ffa_cmdlinearg_getint('numEigModes',m)
    if (modes%nModes < m) then
       write(errMsg,600) m,modes%nModes
600    format('Requested',I4,' eigenmodes. Only',I4,' were found.')
       call reportError (warning_p,errMsg)
       if (modes%nModes < nModes) then
          write(errMsg,601) nModes,nModes-modes%nModes
601       format('Saving',I4,' eigenmodes, the last',I4,' mode(s) being zero.')
          call reportError (empty_p,errMsg)
       end if
    end if

    call writeTimeStepDB (rdb,sys%nStep,sys%time,ierr)
    if (ierr < 0) goto 500

    !! Write the eigenvalues to database file
    do i = 1, nModes
       call writeRDB (rdb,modes%ReVal(i),ierr,writeAsDouble)
       if (complexModes) then
          call writeRDB (rdb,modes%ImVal(i),ierr,writeAsDouble)
          freq = modes%ImVal(i) * 0.5_dp/pi_p
       else
          freq = modes%ReVal(i) * 0.5_dp/pi_p
       end if
       call writeRDB (rdb,freq,ierr)
       call writeRDB (rdb,modes%dampRat(i),ierr)
       if (ierr < 0) goto 500
    end do

    !! Write the triad components of the eigenvectors to database file
    do i = 1, nTriad
       if (triads(i)%id%baseId < 0) then ! Dummy CG triad for generic parts
          cycle
       else if (triads(i)%nDOFs == 3 .or. triads(i)%nDOFs == 6) then
          j = triads(i)%sysDOF
          k = triads(i)%nDOFs + j-1
          do m = 1, nModes
             call writeRDB (rdb,modes%ReVec(j:k,m),ierr,writeAsDouble)
             if (complexModes) then
                call writeRDB (rdb,modes%ImVec(j:k,m),ierr,writeAsDouble)
             end if
             if (ierr < 0) goto 500
          end do
       end if
    end do

    !! Write the generalized DOF-components, if any,
    !! of the eigenvectors to database file
    do i = 1, nSupEl
       if (associated(sups(i)%genDOFs)) then
          j = sups(i)%genDOFs%sysDOF
          k = sups(i)%genDOFs%nDOFs + j-1
          do m = 1, nModes
             call writeRDB (rdb,modes%ReVec(j:k,m),ierr,writeAsDouble)
             if (complexModes) then
                call writeRDB (rdb,modes%ImVec(j:k,m),ierr,writeAsDouble)
             end if
             if (ierr < 0) goto 500
          end do
       end if
    end do

    if (nBytesPerStep < 1) then
       !! Check consistency between header and the actual data
       nBytesPerStep = nBytesPerStep + int(rdb%nBytes)
       call checkStepSize (rdb,nBytesPerStep,ierr)
       if (ierr < 0) goto 500
    end if

    return

500 call reportError (debugFileOnly_p,'writeModesDB')

  end subroutine writeModesDB

end module ModesTypeModule
