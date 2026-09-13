!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file samSolverModule.f90
!> @brief Initialisation of the SAM data structure for the dynamics solver.

!!==============================================================================
!> @brief Initialisation of the SAM data structure for the dynamics solver.

module SamSolverModule

  implicit none

  integer, parameter, private :: traDOF_p = 1 !< Translational DOF type flag
  integer, parameter, private :: rotDOF_p = 2 !< Rotational DOF type flag
  integer, parameter, private :: genDOF_p = 3 !< Generalized DOF type flag


contains

  !!============================================================================
  !> @brief Initializes the SAM data structure for fedem_solver.
  !>
  !> @param[out] sam Data for managing system matrix assembly
  !> @param[out] ierr Error flag
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 21 Jul 2000

  subroutine initSAM (sam,ierr)

    use SamModule       , only : SamType, nullifySAM, version_p
    use allocationModule, only : reAllocate

    type(SamType), intent(out) :: sam
    integer      , intent(out) :: ierr

    !! --- Logic section ---

    call nullifySAM (sam)
    call reAllocate ('initSAM',sam%mpar,80,ierr)
    if (ierr < 0) return

    !! Set pointers to key parameters in MPAR
    sam%nnod   => sam%mpar(1)
    sam%nel    => sam%mpar(2)
    sam%ndof   => sam%mpar(3)
    sam%ndof1  => sam%mpar(4)
    sam%ndof2  => sam%mpar(5)
    sam%nspdof => sam%mpar(6)
    sam%nceq   => sam%mpar(7)
    sam%nsdof  => sam%mpar(8)
    sam%npdof  => sam%mpar(9)
    sam%nddof  => sam%mpar(10)
    sam%neq    => sam%mpar(11)
    sam%nsky   => sam%mpar(12)
    sam%nmmnpc => sam%mpar(15)
    sam%nmmceq => sam%mpar(16)

    sam%mpar = 0
    sam%mpar(80) = version_p ! mark SAM structure with current version

  end subroutine initSAM


  !!============================================================================
  !> @brief Initializes some nodal/dof arrays in SAM.
  !>
  !> @param triads All triads in the model
  !> @param sups All superelements in the model
  !> @param joints All joints in the model
  !> @param mpar Matrix of parameters
  !> @param[out] madof Matrix of accumulated DOFs
  !> @param[out] msc Matrix of status codes
  !> @param[out] dofType Type flag for each DOF
  !> @param[out] ierr Error flag
  !>
  !> @details This subroutine initializes the SAM arrays @a MADOF and @a MSC
  !> based on Triad, SupEl and Joint data. It also sets the type flag
  !> for each DOF (translational, rotational or generalized).
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 21 Jul 2000

  subroutine initSAM_dofStatus (triads,sups,joints,mpar,madof,msc,dofType,ierr)

    use TriadTypeModule           , only : TriadType
    use SupElTypeModule           , only : SupElType
    use MasterSlaveJointTypeModule, only : MasterSlaveJointType
    use allocationModule          , only : reAllocate

    type(TriadType)           , intent(inout) :: triads(:)
    type(SupElType)           , intent(inout) :: sups(:)
    type(MasterSlaveJointType), intent(inout) :: joints(:)
    integer                   , intent(inout) :: mpar(:)
    integer                   , pointer       :: madof(:), msc(:), dofType(:)
    integer                   , intent(out)   :: ierr

    !! Local variables
    integer :: i, idIn, inod, l, ndof, ngen, nnod, nTriads, nSupEls, nJoints

    !! --- Logic section ---

    ierr = 0
    ndof = 0
    ngen = 0
    nnod = 0

    !! Count nodes and DOFs, and assign SAM node numbers to the triads,
    !! superelements (those having generalized DOFs) and joints

    !! SAM nodes on the triads
    nTriads = size(triads)
    do idIn = 1, nTriads
       if (triads(idIn)%nDOFs > 0) then
          nnod = nnod + 1
          ndof = ndof + triads(idIn)%nDOFs
          triads(idIn)%samNodNum = nnod
       else if (triads(idIn)%samNodNum > 0) then
          nnod = nnod + 1
          ngen = ngen + 1
          triads(idIn)%samNodNum = nnod
       end if
    end do

    mpar(27) = nnod - ngen ! number of external nodes
    ngen = 0

    !! SAM nodes on the elements for generalized DOFs
    nSupEls = size(sups)
    do idIn = 1, nSupEls
       if (associated(sups(idIn)%genDOFs)) then
          nnod = nnod + 1
          ngen = ngen + sups(idIn)%genDOFs%nDOFs
          ndof = ndof + sups(idIn)%genDOFs%nDOFs
          sups(idIn)%genDOFs%samNodNum = nnod
       end if
    end do

    mpar(26) = nSupEls
    mpar(22) = ngen

    !! SAM nodes on the joints for handling joint DOFs
    nJoints = size(joints)
    do idIn = 1, nJoints
       if (joints(idIn)%nJointDOFs > 0) then
          nnod = nnod + 1
          ndof = ndof + joints(idIn)%nJointDOFs
          joints(idIn)%samNodNum = nnod
       end if
    end do

    mpar(1) = nnod
    mpar(3) = ndof

    call reAllocate ('initSAM_dofStatus',madof,nnod+1,ierr)
    call reAllocate ('initSAM_dofStatus',msc,ndof,ierr)
    call reAllocate ('initSAM_dofStatus',dofType,ndof,ierr)
    if (ierr < 0) return

    !! Initialize MADOF, MSC and dofType

    madof(1) = 1
    do idIn = 1, nTriads
       inod = triads(idIn)%samNodNum
       if (triads(idIn)%nDOFs > 0) then
          madof(inod+1) = madof(inod) + triads(idIn)%nDOFs
          !! Set status codes for the triad DOFs
          msc(madof(inod) : madof(inod+1)-1) = triads(idIn)%BC
          !! Set the DOF type
          dofType(madof(inod) : madof(inod)+2) = traDof_p
          if (madof(inod+1)-1 > madof(inod)+2) then
             dofType(madof(inod)+3 : madof(inod+1)-1) = rotDof_p
          end if
       else if (inod > 0) then
          madof(inod+1) = madof(inod)
       end if
    end do

    do idIn = 1, nSupEls
       if (associated(sups(idIn)%genDOFs)) then
          inod = sups(idIn)%genDOFs%samNodNum
          madof(inod+1) = madof(inod) + sups(idIn)%genDOFs%nDOFs
          !! Set status codes for the generalized DOFs
          msc(madof(inod) : madof(inod+1)-1) = sups(idIn)%genDOFs%BC
          !! Set the DOF type
          dofType(madof(inod) : madof(inod+1)-1) = genDof_p
       end if
    end do

    do idIn = 1, nJoints
       inod = joints(idIn)%samNodNum
       if (inod > 0) then
          madof(inod+1) = madof(inod) + joints(idIn)%nJointDOFs
          !! Set status codes for the joint DOFs
          msc(madof(inod) : madof(inod+1)-1) = joints(idIn)%BC
          !! Set the DOF type
          do i = 1, joints(idIn)%nJointDOFs
             l = joints(idIn)%jointDofs(i)%lDof
             if (l == 7) l = 3 ! slider DOF
             if (l >= 1 .and. l <= 3) then
                dofType(madof(inod)+i-1) = traDof_p
             else
                dofType(madof(inod)+i-1) = rotDof_p
             end if
          end do
       else ! No joint DOFs (rigid joint), set samNodNum to a valid MADOF index
          joints(idIn)%samNodNum = nnod+1
       end if
    end do

  end subroutine initSAM_dofStatus


  !!============================================================================
  !> @brief Initializes element topology arrays in SAM.
  !>
  !> @param sups All superelements in the model
  !> @param springs All axial springs in the model
  !> @param dampers All axial dampers in the model
  !> @param joints All joints in the model
  !> @param bElems All bushing elements in the model
  !> @param cElems All contact elements in the model
  !> @param uElems All user-defined elements in the model
  !> @param masses All masses in the model
  !> @param tires All tires in the model
  !> @param ctrlStructEl Data for coupled control system modal analysis
  !> @param mpar Matrix of parameters
  !> @param[out] mmnpc Matrix of nodal point correspondances
  !> @param[out] mpmnpc Matrix of pointers to MNPCs
  !> @param[out] ierr Error flag
  !>
  !> @details This subroutine initializes the SAM arrays @a MMNPC and @a MPMNPC
  !> based on SupEl, Spring, Damper, Joint and Tire data, and based on Bushing-,
  !> Contact-, User-defined- and Mass-element data.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 21 Jul 2000

  subroutine initSAM_topology (sups,springs,dampers,joints, &
       &                       bElems,cElems,uElems,masses,tires,ctrlStructEl, &
       &                       mpar,mmnpc,mpmnpc,ierr)

    use SupElTypeModule           , only : SupElType
    use SpringTypeModule          , only : SpringType, isActive
    use DamperTypeModule          , only : DamperType
    use MasterSlaveJointTypeModule, only : MasterSlaveJointType
    use BushingElementTypeModule  , only : BushingElementType
    use ContactElementTypeModule  , only : ContactElementType
    use UserdefElTypeModule       , only : UserdefElType
    use MassTypeModule            , only : MassType
    use TireTypeModule            , only : TireType
    use ControlStructModule       , only : ControlStructType
    use MasterSlaveJointTypeModule, only : getJointId, GetPtrToOwner
    use allocationModule          , only : reAllocate
    use reportErrorModule         , only : AllocationError
    use reportErrorModule         , only : reportError, error_p, empty_p

    type(SupElType)           , intent(inout) :: sups(:)
    type(SpringType)          , intent(inout) :: springs(:)
    type(DamperType)          , intent(inout) :: dampers(:)
    type(MasterSlaveJointType), intent(inout) :: joints(:)
    type(BushingElementType)  , intent(inout) :: bElems(:)
    type(ContactElementType)  , intent(inout) :: cElems(:)
    type(UserdefElType)       , intent(inout) :: uElems(:)
    type(MassType)            , intent(inout) :: masses(:)
    type(TireType)            , intent(inout) :: tires(:)
    type(ControlStructType)   , pointer       :: ctrlStructEl
    integer                   , intent(inout) :: mpar(:)
    integer                   , pointer       :: mmnpc(:), mpmnpc(:)
    integer                   , intent(out)   :: ierr

    !! Local variables
    integer :: i, idIn, nels, nnpc, nSupEls, nBElems, nCElems, nUElems
    integer :: nSpring, nDamper, nMasses, nTires
    logical :: haveSpr(2)
    logical                   , allocatable :: isUsed(:)
    type(MasterSlaveJointType), pointer     :: joint

    !! --- Logic section ---

    ierr = 0
    nnpc = 0

    !! Count number of entries in the element topology array, MMNPC
    nSupEls = size(sups)
    do idIn = 1, nSupEls
       nnpc = nnpc + sups(idIn)%nExtNods
       if (associated(sups(idIn)%genDOFs)) nnpc = nnpc + 1
    end do

    nSpring = 0
    do idIn = 1, size(springs)
       if (.not. isActive(springs(idIn))) cycle ! ignore unused springs
       nSpring = nSpring + 1
       nnpc = nnpc + size(springs(idIn)%triads)
    end do

    nDamper = 0
    do idIn = 1, size(dampers)
       if (.not. dampers(idIn)%dmp%isActive) cycle ! ignore unused dampers
       nDamper = nDamper + 1
       if (associated(dampers(idIn)%triad2)) then
          nnpc = nnpc + 2
       else
          nnpc = nnpc + 1
       end if
    end do

    nBElems = size(bElems)
    do idIn = 1, nBElems
       nnpc = nnpc + 1 + size(bElems(idIn)%spokes)
    end do

    nCElems = size(cElems)
    do idIn = 1, nCElems
       nnpc = nnpc + 1 + size(cElems(idIn)%CSurf%CPoints)
    end do

    nUElems = size(uElems)
    do idIn = 1, nUElems
       nnpc = nnpc + size(uElems(idIn)%triads)
    end do

    nMasses = size(masses)
    nnpc    = nnpc + nMasses

    nTires  = size(tires)
    nnpc    = nnpc + nTires

    nels = nSupEls + nSpring + nDamper + nMasses + nTires &
         + nBElems + nCElems + nUElems

    !! Joints are counted as elements only if they have a multi-DOF spring
    do idIn = 1, size(joints)
       if ( associated(joints(idIn)%springEl) .or. &
            associated(joints(idIn)%sprFric) ) then
          nels = nels + 1
          nnpc = nnpc + 1
       end if
    end do

    !! Add control structure interaction as an element, if it exists
    if (associated(ctrlStructEl)) then
       nels = nels + 1
       nnpc = nnpc + size(ctrlStructEl%samMNPC)
    end if

    mpar(2)  = nels
    mpar(15) = nnpc

    call reAllocate ('initSAM_topology',mmnpc,nnpc,ierr)
    call reAllocate ('initSAM_topology',mpmnpc,nels+1,ierr)
    if (ierr < 0) return

    allocate(isUsed(mpar(1)),STAT=ierr)
    if (ierr /= 0) then
       ierr = AllocationError('initSAM_topology')
       return
    end if

    !! Initialize MMNPC and MPMNPC, and assign SAM element numbers to the
    !! superelements, springs and dampers, bushing-, contact- and mass-elements

    mpmnpc(1) = 1
    nels      = 0
    isUsed    = .false.

    !! Superelements, number of nodes equals the number of triads
    !! + one internal node to handle generalized DOFs, if any
    do idIn = 1, nSupEls
       nels = nels + 1
       nnpc = sups(idIn)%nExtNods
       sups(idIn)%samElNum = nels
       do i = 1, nnpc
          mmnpc(mpmnpc(nels)+i-1) = sups(idIn)%triads(i)%p%samNodNum
          isUsed(sups(idIn)%triads(i)%p%samNodNum) = .true.
       end do
       if (associated(sups(idIn)%genDOFs)) then
          mmnpc(mpmnpc(nels)+nnpc) = sups(idIn)%genDOFs%samNodNum
          isUsed(sups(idIn)%genDOFs%samNodNum) = .true.
          nnpc = nnpc + 1
       end if
       mpmnpc(nels+1) = mpmnpc(nels) + nnpc
    end do

    !! Axial (or global) springs, number of nodes equals the number of triads
    do idIn = 1, size(springs)
       if (.not. isActive(springs(idIn))) cycle ! ignore unused springs

       nels = nels + 1
       nnpc = size(springs(idIn)%triads)
       springs(idIn)%samElNum = nels
       do i = 1, nnpc
          mmnpc(mpmnpc(nels)+i-1) = springs(idIn)%triads(i)%p%samNodNum
          isUsed(springs(idIn)%triads(i)%p%samNodNum) = .true.
       end do
       mpmnpc(nels+1) = mpmnpc(nels) + nnpc

    end do

    !! Joint springs, one element node for each spring
    do idIn = 1, size(joints)
       if (joints(idIn)%samNodNum > mpar(1)) cycle

       isUsed(joints(idIn)%samNodNum) = .true.
       haveSpr(1) = associated(joints(idIn)%springEl)
       haveSpr(2) = associated(joints(idIn)%sprFric)
       if (.not. any(haveSpr)) cycle

       nels = nels + 1
       if (haveSpr(1)) joints(idIn)%springEl%samElNum = nels
       if (haveSpr(2)) joints(idIn)%sprFric%samElNum = nels

       !! Connect to the internal node of the joint
       mmnpc(mpmnpc(nels)) = joints(idIn)%samNodNum
       mpmnpc(nels+1) = mpmnpc(nels) + 1

    end do

    !! Dampers, number of nodes equals the number of triads for axial dampers
    !! and one for joint dampers
    do idIn = 1, size(dampers)
       if (.not. dampers(idIn)%dmp%isActive) cycle ! ignore unused dampers

       nels = nels + 1
       dampers(idIn)%samElNum = nels

       if (associated(dampers(idIn)%triad2)) then

          nnpc = 2 ! Axial damper, connect to the two triad nodes
          mmnpc(mpmnpc(nels))   = dampers(idIn)%triad1%samNodNum
          mmnpc(mpmnpc(nels)+1) = dampers(idIn)%triad2%samNodNum
          isUsed(dampers(idIn)%triad1%samNodNum) = .true.
          isUsed(dampers(idIn)%triad2%samNodNum) = .true.

       else if (associated(dampers(idIn)%triad1)) then

          nnpc = 1 ! Damper to ground, connect to the single triad node
          mmnpc(mpmnpc(nels)) = dampers(idIn)%triad1%samNodNum
          isUsed(dampers(idIn)%triad1%samNodNum) = .true.

       else

          !! Joint damper, connect to the internal node of the joint
          joint => GetPtrToOwner(joints,dampers(idIn)%dmp)
          if (.not. associated(joint)) then
             nnpc = 0
          else if (joint%samNodNum > mpar(1)) then
             nnpc = 0
          else
             nnpc = 1
             mmnpc(mpmnpc(nels)) = joint%samNodNum
             isUsed(joint%samNodNum) = .true.
          end if

       end if

       mpmnpc(nels+1) = mpmnpc(nels) + nnpc

    end do

    !! Bushing elements, number of nodes equals the number of spokes
    !! + centre triad
    do idIn = 1, nBElems
       nels = nels + 1
       nnpc = 1 + size(bElems(idIn)%spokes)
       bElems(idIn)%samElNum = nels
       mmnpc(mpmnpc(nels)) = bElems(idIn)%centre%samNodNum
       isUsed(bElems(idIn)%centre%samNodNum) = .true.
       do i = 1, nnpc-1
          mmnpc(mpmnpc(nels)+i) = bElems(idIn)%spokes(i)%triad%samNodNum
          isUsed(bElems(idIn)%spokes(i)%triad%samNodNum) = .true.
       end do
       mpmnpc(nels+1) = mpmnpc(nels) + nnpc
    end do

    !! Contact elements, number of nodes equals the number of curve points
    !! + contact (follower) triad
    do idIn = 1, nCElems
       nels = nels + 1
       nnpc = 1 + size(cElems(idIn)%CSurf%CPoints)
       cElems(idIn)%samElNum = nels
       mmnpc(mpmnpc(nels)) = cElems(idIn)%triad1%samNodNum
       isUsed(cElems(idIn)%triad1%samNodNum) = .true.
       do i = 1, nnpc-1
          mmnpc(mpmnpc(nels)+i) = cElems(idIn)%CSurf%CPoints(i)%triad%samNodNum
          isUsed(cElems(idIn)%CSurf%CPoints(i)%triad%samNodNum) = .true.
       end do
       mpmnpc(nels+1) = mpmnpc(nels) + nnpc
    end do

    !! User-defined elements, number of nodes equals the number of triads
    do idIn = 1, nUElems
       nels = nels + 1
       nnpc = size(uElems(idIn)%triads)
       uElems(idIn)%samElNum = nels
       do i = 1, nnpc
          mmnpc(mpmnpc(nels)+i-1) = uElems(idIn)%triads(i)%p%samNodNum
          isUsed(uElems(idIn)%triads(i)%p%samNodNum) = .true.
       end do
       mpmnpc(nels+1) = mpmnpc(nels) + nnpc
    end do

    !! Triad masses, number of nodes is 1
    do idIn = 1, nMasses
       nels = nels + 1
       nnpc = 1
       masses(idIn)%samElNum = nels
       mmnpc(mpmnpc(nels)) = masses(idIn)%triad%samNodNum
       isUsed(masses(idIn)%triad%samNodNum) = .true.
       mpmnpc(nels+1) = mpmnpc(nels) + nnpc
    end do

    !! Tires, number of nodes is 1
    do idIn = 1, nTires
       nels = nels + 1
       nnpc = 1
       tires(idIn)%samElNum = nels
       mmnpc(mpmnpc(nels)) = tires(idIn)%RimTriad%samNodNum
       isUsed(tires(idIn)%RimTriad%samNodNum) = .true.
       mpmnpc(nels+1) = mpmnpc(nels) + nnpc
    end do

    !! Control structure interaction (for eigenvalue analysis, primarily)
    if (associated(ctrlStructEl)) then
       nels = nels + 1
       nnpc = size(ctrlStructEl%samMNPC)
       ctrlStructEl%samElNum = nels
       mpmnpc(nels+1) = mpmnpc(nels) + nnpc
       mmnpc(mpmnpc(nels):mpmnpc(nels+1)-1) = ctrlStructEl%samMNPC
    end if

    !! Consistency check: All dependent nodes must be connected to something
    if (count(isUsed) < mpar(1)) then
       do idIn = 1, size(joints)
          if (.not. isUsed(joints(idIn)%STriad%samNodNum)) then
             ierr = ierr + 1
             call reportError (error_p,'The dependent triad in '// &
                  trim(getJointId(joints(idIn)))//' is not connected.')
          end if
       end do
       if (ierr > 0) then
          call reportError (empty_p, &
               'The model is therefore inconsistent and must be corrected.', &
               addString='initSAM_topology')
       end if
    end if

    deallocate(isUsed)

  end subroutine initSAM_topology


  !!============================================================================
  !> @brief Initializes linear coupling arrays in SAM.
  !>
  !> @param motions All prescribed motions in the model
  !> @param joints All joints in the model
  !> @param higherPairs All higher pair objects in the model
  !> @param mpar Matrix of parameters
  !> @param[in] madof Matrix of accumulated DOFs
  !> @param[in] msc Matrix of status codes
  !> @param[out] mpmceq Matrix of pointers to MCEQs
  !> @param[out] mmceq Matrix of constraint equation definitions
  !> @param[out] ttcc Table of constraint equation coefficients
  !> @param[out] mpreac Matrix of pointers to reaction forces
  !> @param[out] ierr Error flag
  !>
  !> @details This subroutine initializes the SAM arrays @a MPMCEQ, @a MMCEQ and
  !> @a TTCC based on Motion and Joint data. The array @a TTCC is here
  !> only initialized to zero (except for the higher pair constraints).
  !> The numerical values of the constraint equations are calculated by
  !> the subroutines updateDependencies and updatePrescribedMotions.
  !> The array @a MPREAC is also initialized by this subroutine
  !> based on Motion data and the fixed DOFs.
  !>
  !> @callergraph
  !>
  !> @author Bjorn Haugen
  !> @date 16 Oct 2000
  !>
  !> @author Knut Morten Okstad
  !> @date 13 Oct 2004

  subroutine initSAM_constraints (motions,joints,higherPairs,mpar,madof,msc, &
       &                          mpmceq,mmceq,ttcc,mpreac,ierr)

    use MotionTypeModule          , only : MotionType, dp
    use MasterSlaveJointTypeModule, only : MasterSlaveJointType, HigherPairType
    use MasterSlaveJointTypeModule, only : getJointId, GetNumberOfMasterDofs
    use IdTypeModule              , only : getId
    use allocationModule          , only : reAllocate
    use reportErrorModule         , only : InternalError
    use reportErrorModule         , only : reportError, debugFileOnly_p, error_p

    type(MotionType)          , intent(inout) :: motions(:)
    type(MasterSlaveJointType), intent(inout), target :: joints(:)
    type(HigherPairType)      , intent(inout) :: higherPairs(:)
    integer                   , intent(inout) :: mpar(:), msc(:)
    integer                   , intent(in)    :: madof(:)
    integer                   , pointer       :: mpmceq(:), mmceq(:), mpreac(:)
    real(dp)                  , pointer       :: ttcc(:)
    integer                   , intent(out)   :: ierr

    !! Local variables
    integer           :: i, j, k, lerr, iceq, sDof, mDof, ip, ipR, ipMceq
    integer           :: nMotion, nJoints, nHPs, nTriads, nDeps, nceq, nmmceq
    character(len=16) :: chDof

    !! --- Logic section ---

    nMotion = size(motions)
    nJoints = size(joints)
    nHPs    = size(higherPairs)

    !! Ensure that the given DOF status codes are valid
    if (any(msc < 0) .or. any(msc > 2)) then
       ierr = InternalError('initSAM_constraints: Invalid MSC array')
       return
    end if

    !! Mark the dependent DOFs of all joints
    ierr = 0
    do i = 1, nJoints
       sDof = madof(joints(i)%STriad%samNodNum)
       do j = 1, size(joints(i)%slaveDofs)
          if (msc(sDof) > 0) then
             msc(sDof) = -i - nMotion
          else
             ierr = ierr + 1
             write(chDof,"('System DOF',I6)") sDof
             call reportError (error_p,'Invalid dependent DOF specified for '//&
                  getJointId(joints(i)),chDof//' is already marked as fixed')
          end if
          sDof = sDof + 1
       end do
    end do

    !! Count the number of constraint equations and constraint coefficients
    nceq   = nMotion + nHPs
    nmmceq = nMotion + 2*nHPs

    do i = 1, nJoints
       nDeps   = size(joints(i)%slaveDofs)
       nTriads = countIndependentDofs(joints(i),ierr)
       nceq    = nceq   + nDeps
       nmmceq  = nmmceq + nDeps*(1+nTriads)
    end do

    !! Initialize some key numbers in MPAR
    mpar(7)  = nceq
    mpar(16) = nmmceq

    !! Allocate the pointer array and coefficient table
    lerr = 0
    iceq = 0
    call reAllocate ('initSAM_constraints',mpmceq,nceq+1,lerr)
    call reAllocate ('initSAM_constraints',mmceq,max(1,nmmceq),lerr)
    call reAllocate ('initSAM_constraints',ttcc ,max(1,nmmceq),lerr)
    if (lerr < 0) goto 900

    mpmceq = 0
    mmceq = 0
    ttcc = 0.0_dp

    if (nMotion > 0 .or. any(msc == 0)) then
       !! Allocate reaction force pointer array for fixed and prescribed DOFs
       call reAllocate ('initSAM_constraints',mpreac,mpar(3),lerr)
       if (lerr < 0) goto 900
    end if

    if (associated(mpreac)) mpreac = 0

    ipR = 1
    if (any(msc == 0)) then
       !! Set pointers into reaction force array for the fixed DOFs
       do i = 1, mpar(3)
          if (msc(i) == 0) then
             mpReac(i) = ipR
             ipR = ipR + 1
          end if
       end do
    end if

    !! Initialize the linear couplings according to prescribed motions
    ip   = 1
    do i = 1, nMotion
       if (associated(motions(i)%triad)) then
          sDof = madof(motions(i)%triad%samNodNum) + motions(i)%dof - 1
       else if (associated(motions(i)%joint)) then
          sDof = madof(motions(i)%joint%samNodNum) + motions(i)%dof - 1
       else
          ierr = InternalError('initSAM_constraints: Incomplete motion(s)')
          goto 900
       end if

       if (msc(sDof) < 1) then
          ierr = ierr + 1
          write(chDof,"('System DOF',I6)") sDof
          call reportError (error_p,'Invalid DOF specification for Motion'// &
               &            getId(motions(i)%id),chDof// &
               &            ' is already marked as suppressed')
          cycle
       end if

       !! Add a constraint equation for the prescribed DOF
       iceq = iceq + 1
       mpmceq(iceq) = ip
       msc(sDof) = -iceq ! update status code to suppressed (prescribed)
       mmceq(ip) = sDof

       if (associated(motions(i)%triad)) then
          motions(i)%triad%BC(motions(i)%dof) = -iceq ! mark DOF as prescribed
       end if

       !! Set short-cut pointers to system DOF number and the prescribed value
       motions(i)%sDof = sDof
       motions(i)%cc => ttcc(ip)

       ip = ip + 1

       !! Set pointer into reaction force array for the prescribed DOF
       mpreac(sDof) = ipR
       ipR = ipR + 1
    end do

    !! Store the number of reaction force components in MPAR
    mpar(17) = ipR-1
    ipR = 1

    !! Initialize the linear couplings according to joint constraints
    do i = 1, nJoints

       lerr = ierr
       sDof = madof(joints(i)%STriad%samNodNum)
       do j = 1, size(joints(i)%slaveDofs)

          !! Add a constraint equation for the dependent DOF
          iceq = iceq + 1
          mpmceq(iceq) = ip
          mmceq(ip) = sDof
          ip = ip + 1

          !! Find the independent DOFs of this constraint equation
          !! and store them in MMCEQ
          call connectIndependentDofs (joints(i))
          if (ierr < 0) goto 900

          !! Set short-cut pointers for constraint equations and coefficients
          ipMceq = mpmceq(iceq)
          joints(i)%slaveDofs(j)%mceq => mmceq(ipMceq:ip-1)
          joints(i)%slaveDofs(j)%tcc  => ttcc (ipMceq:ip-1)

          sDof = sDof + 1
       end do

       !! Set start index to MCEQ of this constraint for each independent triad
       ipMceq = 2
       do k = 1, size(joints(i)%JMTriads)
          joints(i)%JMTriads(k)%ipMceq = ipMceq
          if (joints(i)%JMTriads(k)%nDOFs < 1) cycle
          mDof = madof(joints(i)%JMTriads(k)%triad%samNodNum)
          if (msc(mdof) >= -nMotion) then
             ipMceq = ipMceq + joints(i)%JMTriads(k)%nDOFs
          else if (.not. associated(joints(i)%chain)) then
             !! The independent joint triad is dependent in another joint
             j = -msc(mdof)-nMotion
             if (size(joints(j)%slaveDofs) /= joints(i)%JMTriads(k)%nDOFs) then
                ierr = ierr + 1
                call reportError (error_p,trim(getJointId(joints(i)))// &
                     ' can not be chained with',getJointId(joints(j)), &
                     'because the number of dependent DOFs of the latter ', &
                     'does not match the number of independent triad DOFs '// &
                     'of the former')
                exit
             else if (j >= i) then
                ierr = ierr + 1
                call reportError (error_p,trim(getJointId(joints(i)))// &
                     ' can not be chained with ',getJointId(joints(j)), &
                     'because the latter is ordered after the former')
                exit
             end if
             joints(i)%chain => joints(j)
             ipMceq = ipMceq + GetNumberOfMasterDofs(joints(i)%chain)
          else
             ierr = ierr + 1
             call reportError (error_p,'Invalid '//getJointId(joints(i)), &
                  'Only one joint triad can be dependent in another joint')
             exit
          end if
       end do
       if (ierr > lerr) cycle

       !! Set system DOF and index to MCEQ of this contraint for each joint DOF
       mDof = madof(joints(i)%samNodNum)
       do k = 1, joints(i)%nJointDOFs
          joints(i)%jointDofs(k)%sysDof = mDof
          joints(i)%jointDofs(k)%ipMceq = ipMceq
          ipMceq = ipMceq + 1
          mDof = mDof + 1
       end do

       !! Sanity check
       do j = 1, size(joints(i)%slaveDofs)
          if (ipMceq-1 /= size(joints(i)%slaveDofs(j)%mceq)) then
             ierr = InternalError('initSAM_constraints: Logic error 3')
             goto 900
          end if
       end do

    end do

    !! Add linear couplings according to higher pairs
    do i = 1, nHPs
       if (ip+1 > size(mmceq)) then
          ierr = InternalError('initSAM_constraints: Logic error 4')
          goto 900
       end if

       !! Get the system DOFs for the input and output joints
       call getHPSystemDOFs (higherPairs(i)%slaveJoint, &
            &                higherPairs(i)%masterJoint, &
            &                higherPairs(i)%slaveJointDof, &
            &                higherPairs(i)%masterJointDof, &
            &                sDof,mDof)
       if (sDof < 1 .or. mDof < 1) cycle

       !! Add one constraint equation for each higher pair
       iceq = iceq + 1
       mpmceq(iceq) = ip

       !! Set output (dependent) part of the constraint equation
       msc(sDof) = -iceq ! set the status code to suppressed
       mmceq(ip) = sDof
       ip = ip + 1

       !! Set input (independent) part of the constraint equation
       mmceq(ip) = mDof
       ttcc(ip) = higherPairs(i)%coeff
       ip = ip + 1

    end do

    iceq = iceq + 1
    mpmceq(iceq) = ip ! pointer to last (dummy) equation

    !! Restore zero entries in MSC for all suppressed DOFs
    do i = 1, size(msc)
       if (msc(i) < 0) msc(i) = 0
    end do

    if (ierr > 0) call reportError (debugFileOnly_p,'initSAM_constraints')

900 if (iceq < nceq+1) then
       mpmceq(iceq+1:) = 0 ! To avoid crash later due to uninitialized indices
    end if
    if (lerr < 0 .and. ierr >= 0) ierr = lerr

  contains

    !> @brief Counts the number of free independent DOFs in a joint.
    !> @details This is a recursive function accounting for possible chaining.
    recursive function countIndependentDofs (joint,ierr) result (nDofs)
      type(MasterSlaveJointType), intent(in)    :: joint
      integer                   , intent(inout) :: ierr
      integer :: nDofs            ! return variable
      integer :: jIdx, k, l, mDof ! local stack variables
      integer :: recursionLev = 0 ! static variable

      !! The joint variables are independent DOFs
      nDofs = joint%nJointDOFs

      if (recursionLev < 10) then
         recursionLev = recursionLev + 1
      else
         ierr = ierr + 1 ! avoid infinite recursion in case of modeling error
         call reportError (error_p,'Possibly infinite joint chain '// &
              &            'detected for '//getJointId(joint), &
              &            'The maximum allowed length of a joint chain is 10')
         return
      end if

      !! Check the joint triad DOFs
      do k = 1, size(joint%JMTriads)
         if (joint%JMTriads(k)%nDOFs < 1) cycle

         jIdx = 0
         mDof = madof(joint%JMTriads(k)%triad%samNodNum)
         do l = 1, joint%JMTriads(k)%nDOFs
            if (msc(mDof) < -nMotion) then
               !! This independent triad DOF is dependent in another joint
               !! so count the independent DOFs of that joint instead
               if (jIdx == 0) then
                  if (l == 1) then
                     jIdx = -msc(mDof) - nMotion
                     nDofs = nDofs + countIndependentDofs(joints(jIdx),ierr)
                  else
                     goto 900
                  end if
               end if
            else if (jIdx == 0) then
               nDofs = nDofs + 1
            else
               goto 900
            end if
            mDof = mDof + 1
         end do

      end do

      goto 999

900   ierr = ierr + 1
      call reportError (error_p,'Invalid chaining for '//getJointId(joint), &
           'Some DOFs in Triad'//trim(getId(joint%JMTriads(k)%triad%id))// &
           'are dependent in another joint whereas others are free/suppressed')
999   recursionLev = recursionLev - 1

    end function countIndependentDofs

    !> @brief Builds the linear coupling arrays @a mmceq/mpmceq for a joint.
    !> @details This is a recursive subroutine accounting for possible chaining.
    recursive subroutine connectIndependentDofs (joint)
      type(MasterSlaveJointType), intent(in) :: joint

      integer :: jIdx, k, l, mDof ! local stack variables

      !! DOFs at all independent joint triads
      do k = 1, size(joint%JMTriads)
         if (joint%JMTriads(k)%nDOFs < 1) cycle

         jIdx = 0
         mDof = madof(joint%JMTriads(k)%triad%samNodNum)
         do l = 1, joint%JMTriads(k)%nDOFs
            if (msc(mDof) < -nMotion) then
               !! This independent triad DOF is dependent in another joint
               !! so connect to the independent DOFs of the other joint instead
               if (jIdx == 0 .and. l == 1) then
                  jIdx = -msc(mDof) - nMotion
                  call connectIndependentDofs (joints(jIdx))
                  if (ierr < 0) return
               end if
            else if (ip > size(mmceq)) then
               ierr = InternalError('initSAM_constraints: Logic error 1')
               return
            else if (msc(mDof) > 0) then
               !! Free independent DOF
               mmceq(ip) = mDof
               ip = ip + 1
            else
               !! Fixed or prescribed independent DOF
               mmceq(ip) = -mDof
               ip = ip + 1
            end if
            mDof = mDof + 1
         end do

      end do

      !! Joint DOFs
      mDof = madof(joint%samNodNum)
      do k = 1, joint%nJointDOFs
         if (ip > size(mmceq)) then
            ierr = InternalError('initSAM_constraints: Logic error 2')
            return
         else if (msc(mDof) > 0) then
            !! Free joint DOF
            mmceq(ip) = mDof
         else if (msc(mDof) < -nMotion) then
            !! Should never get here
            ierr = InternalError('initSAM_constraints: Corrupted MSC')
            return
         else
            !! Fixed or prescribed joint DOF
            mmceq(ip) = -mDof
         end if
         ip   = ip + 1
         mDof = mDof + 1
      end do

    end subroutine connectIndependentDofs

    !> @brief Gets the system DOFs for a higher pair object.
    subroutine getHPSystemDOFs (out,inp,locSDOF,locMDOF,sysSDOF,sysMDOF)
      type(MasterSlaveJointType), intent(inout) :: out
      type(MasterSlaveJointType), intent(in)    :: inp
      integer                   , intent(in)    :: locSDOF, locMDOF
      integer                   , intent(out)   :: sysSDOF, sysMDOF

      !! Get output system DOF and pointer to the dependent DOF in joint MCEQ
      if (locSDOF > 0 .and. locSDOF <= out%nJointDOFs) then
         sysSDOF = out%jointDofs(locSDOF)%sysDOF
         ipMceq  = out%jointDofs(locSDOF)%ipMceq
         out%jointDofs(locSDOF)%coeff = higherPairs(i)%coeff
      else
         sysSDOF = 0
         ierr    = ierr + 1
         call reportError (error_p,'Invalid output DOF specification'// &
              &            ' for Higher Pair'//getId(higherPairs(i)%id))
         return
      end if

      if (msc(sysSDOF) < 1) then
         write(chDof,"('System DOF',I6)") sysSDOF
         sysSDOF = 0
         ierr    = ierr + 1
         call reportError (error_p,'Invalid output DOF specification'// &
              &            ' for Higher Pair'//getId(higherPairs(i)%id), &
              &            chDof//' is already marked as suppressed')
         return
      end if

      !! Get the input system DOF
      if (locMDOF > 0 .and. locMDOF <= inp%nJointDOFs) then
         sysMDOF = inp%jointDofs(locMDOF)%sysDOF
      else
         sysMDOF = 0
         ierr    = ierr + 1
         call reportError (error_p,'Invalid input DOF specification'// &
              &            ' for Higher Pair'//getId(higherPairs(i)%id))
         return
      end if

      if (msc(sysMDOF) < 1) then
         write(chDof,"('System DOF',I6)") sysMDOF
         sysMDOF = 0
         ierr    = ierr + 1
         call reportError (error_p,'Invalid input DOF specification'// &
              &            ' for Higher Pair'//getId(higherPairs(i)%id), &
              &            chDof//' is marked as suppressed or prescribed', &
              &            'This is currently not supported')
         return
      end if

      !! Replace the independent DOFs in the joint constraint equations
      do j = 1, size(out%slaveDofs)
         out%slaveDofs(j)%mceq(ipMceq) = sysMDOF
      end do

    end subroutine getHPSystemDOFs

  end subroutine initSAM_constraints


  !!============================================================================
  !> @brief Does final initialization of SAM (nodal reordering and preassembly).
  !>
  !> @param sam Data for managing system matrix assembly
  !> @param[out] sysMat System matrix to perform preassembly for
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 26 Jul 2000

  subroutine initSAM_preassembly (sam,sysMat,ierr)

    use SamModule             , only : SamType, check
    use SysMatrixTypeModule   , only : SysMatrixType, denseMatrix_p, pardiso_p
    use SysMatrixTypeModule   , only : skylineMatrix_p, sparseMatrix_p
    use AsmExtensionModule    , only : castToInt8, csAllocPointerArrays
    use AsmExtensionModule    , only : csRenumber, csPreAssemble
    use allocationModule      , only : reAllocate
    use reportErrorModule     , only : getErrorFile
    use reportErrorModule     , only : reportError, debugFileOnly_p
    use FFaCmdLineArgInterface, only : ffa_cmdlinearg_isTrue

    type(SamType)      , intent(inout) :: sam
    type(SysMatrixType), intent(out)   :: sysMat
    integer            , intent(out)   :: ierr

    !! Local variables
    integer :: i, lpu, matrixType

    !! --- Logic section ---

    call reAllocate ('initSAM_preassembly',sam%minex,sam%nnod,ierr)
    if (ierr < 0) return

    do i = 1, sam%nnod
       sam%minex(i) = i
    end do

    if (ffa_cmdlinearg_isTrue('skylinesolver')) then
       matrixType = skylineMatrix_p
    else if (ffa_cmdlinearg_isTrue('densesolver')) then
       matrixType = denseMatrix_p
    else if (ffa_cmdlinearg_isTrue('pardiso')) then
       matrixType = pardiso_p
    else if (ffa_cmdlinearg_isTrue('pardisoIndef')) then
       matrixType = -pardiso_p
    else
       matrixType = sparseMatrix_p
       call castToInt8 (sam,ierr)
       if (ierr < 0) goto 900
    end if
    lpu = getErrorFile()
    call csAllocPointerArrays (sam,sysMat,matrixType,lpu,ierr)
    if (ierr < 0) goto 900

    if (ffa_cmdlinearg_isTrue('nosolveropt')) then
       if (matrixType == skylineMatrix_p) sam%mnnn => sam%minex
    else
       call csRenumber (sam,sysMat,lpu,ierr)
       if (ierr < 0) goto 900
    end if

    call csPreassemble (sam,sysMat,.false.,lpu,ierr)
    if (ierr < 0) goto 900

    sam%meqn => sysMat%meqn
    call check (sam,lpu,ierr,.true.)
    if (ierr == 0) return

900 continue
    call reportError (debugFileOnly_p,'initSAM_preassembly')

  end subroutine initSAM_preassembly

end module SamSolverModule
