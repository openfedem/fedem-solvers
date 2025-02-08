!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file initiateJointTypeModule.f90
!> @brief Initialization of joint objects from the solver input file.

!!==============================================================================
!> @brief Initialization of joint objects from the solver input file.

module initiateJointTypeModule

  implicit none

contains

  !!============================================================================
  !> @brief Initializes joint objects with data from the solver input file.
  !>
  !> @param[in] infp File unit number for the solver input file
  !> @param[in] triads Array of all triads in the model
  !> @param[in] gCurves Array of all glider curve objects in the model
  !> @param[in] springs Array of all base spring objects in the model
  !> @param[in] dampers Array of all base damper objects in the model
  !> @param[in] dmpElms Array of all damper elements in the model
  !> @param[in] engines Array of all engines in the model
  !> @param joints Array of all joint objects in the model
  !> @param[in] frictionSets Array of all friction parameter sets in the model
  !> @param[out] err Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Karl Erik Thoresen / Bjorn Haugen
  !> @date Nov 1998
  !>
  !> @author Bjorn Haugen
  !> @date Oct 2000
  !>
  !> @author Knut Morten Okstad
  !> @date Jul 2002

  subroutine ReadJoints (infp,triads,gCurves,springs,dampers,dmpElms,engines, &
       &                 joints,frictionSets,err)

    use IdTypeModule              , only : IdType, ldesc_p
    use IdTypeModule              , only : initId, getId, ReportInputError
    use TriadTypeModule           , only : TriadType, GetPtrToId, dp
    use SpringTypeModule          , only : SpringBaseType, GetPtrToId
    use SpringTypeModule          , only : nullifySpring
    use DamperTypeModule          , only : DamperBaseType, DamperType
    use DamperTypeModule          , only : GetPtrToId
    use FunctionTypeModule        , only : EngineType, GetPtrToId
    use FrictionTypeModule        , only : FrictionParameterType
    use ContactSurfaceModule      , only : GliderCurveType, GetPtrToId
    use MasterSlaveJointTypeModule, only : MasterSlaveJointType, getJointId
    use MasterSlaveJointTypeModule, only : DeallocateJoints, NullifyJoint
    use MasterSlaveJointTypeModule, only : rotParamTypes_p
    use MasterSlaveJointTypeModule, only : FREE_p, PRISMATIC_p, CAM_p, AXIAL_p
    use inputUtilities            , only : iuGetNumberOfEntries
    use inputUtilities            , only : iuSetPosAtNextEntry, iuCharToInt
    use rotationModule            , only : vec_to_mat, orthonorm3
    use manipMatrixModule         , only : matmul34, invert34
    use searchAndSortModule       , only : quickSortInt
    use progressModule            , only : lterm
    use reportErrorModule         , only : AllocationError
    use reportErrorModule         , only : reportError, debugFileOnly_p, error_p
    use reportErrorModule         , only : warningFileOnly_p, warning_p, note_p
    use FFaCmdLineArgInterface    , only : ffa_cmdlinearg_getbool

    integer                    , intent(in)    :: infp
    type(TriadType)            , intent(in)    :: triads(:)
    type(GliderCurveType)      , intent(in)    :: gCurves(:)
    type(SpringBaseType)       , intent(in)    :: springs(:)
    type(DamperBaseType)       , intent(in)    :: dampers(:)
    type(DamperType)           , intent(inout) :: dmpElms(:)
    type(EngineType)           , intent(in)    :: engines(:)
    type(MasterSlaveJointType) , pointer       :: joints(:)
    type(FrictionParameterType), intent(in)    :: frictionSets(:)
    integer                    , intent(out)   :: err

    !! Local variables
    integer :: i, j, idIn, iMat, iDof, lDof, nJoints, nSDOFs, nJTriads, maxJTD
    integer :: springCpl(2), stat, lerr
    logical :: allSec, allRest, allJoint, allLen, allVel, allAcc, allForce
    logical :: ignoreIC, nullifyEvec, depDirInJointDir, haveSpringInterConn
    logical :: needAddBC, noFEnodeAtJoints, compFailure
    real(dp):: PJG(3,4), PJGinv(3,4), PTGinv(3,4)
    type(MasterSlaveJointType), pointer :: joint

    !! Type strings defining stiffness couplings
    character(len=4), parameter :: cplTypes_p(5) =(/ 'NONE', &
         &                                           'XY  ', &
         &                                           'YZ  ', &
         &                                           'ZX  ', &
         &                                           'XYZ ' /)

    !! Type enums defining stiffness couplings
    integer, parameter :: NONE_p = 1, &
         &                XY_p   = 2, &
         &                YZ_p   = 3, &
         &                ZX_p   = 4, &
         &                XYZ_p  = 5

    !! Define the MASTERSLAVEJOINT namelist
    integer, parameter :: MAX_JVAR = 10
    integer  :: id, extId(10), type, version, slaveId, masterId, &
         &      nJointVars, JointVarDefs(2,MAX_JVAR), BC(MAX_JVAR), &
         &      springId(MAX_JVAR), damperId(MAX_JVAR), &
         &      frictionSetId(MAX_JVAR), frictionSpringId(MAX_JVAR), &
         &      frictionEngineId(MAX_JVAR), saveVar(5,MAX_JVAR)
    character(ldesc_p) :: extDescr
    character(len=20)  :: rotParam, tranSpringCpl, rotSpringCpl
    real(dp) :: initPosInGlobal(4,3), JVarInitVal(MAX_JVAR), &
         &      JVarInitVel(MAX_JVAR), JVarInitAcc(MAX_JVAR), &
         &      camThickness

    namelist /MASTERSLAVEJOINT/ id, extId, extDescr, type, version, &
         &    camThickness, slaveId, masterId, initPosInGlobal, &
         &    nJointVars, JointVarDefs, JVarInitVal, JVarInitVel, JVarInitAcc, &
         &    BC, springId, damperId, frictionSetId, &
         &    frictionSpringId, frictionEngineId, saveVar, &
         &    rotParam, tranSpringCpl, rotSpringCpl

    !! --- Logic section ---

    nJoints = iuGetNumberOfEntries(infp,'&MASTERSLAVEJOINT',err)
    if (err /= 0) goto 900
    write(lterm,*) 'Number of &MASTERSLAVEJOINT =',nJoints

    call DeallocateJoints (joints)
    allocate(joints(nJoints),STAT=stat)
    if (stat /= 0) then
       err = AllocationError('ReadJoints 10')
       return
    end if

    !! Check if we need to account for the additional boundary conditions (BC=2)
    call ffa_cmdlinearg_getbool ('initEquilibrium',needAddBC)
    if (.not. needAddBC) then
       call ffa_cmdlinearg_getint ('numEigModes',stat)
       if (stat > 0) call ffa_cmdlinearg_getbool ('addBC_eigensolver',needAddBC)
    end if

    call ffa_cmdlinearg_getbool ('ignoreIC',ignoreIC)
    call ffa_cmdlinearg_getbool ('nullifyEvec',nullifyEvec)
    call ffa_cmdlinearg_getbool ('depDirInJointDir',depDirInJointDir)
    call ffa_cmdlinearg_getbool ('noFEnodeAtJoints',noFEnodeAtJoints)
    call ffa_cmdlinearg_getbool ('allSecondaryVars',allSec)
    call ffa_cmdlinearg_getbool ('allRestartVars',allRest)
    call ffa_cmdlinearg_getbool ('allJointVars',allJoint)
    call ffa_cmdlinearg_getbool ('allLengthVars',allLen)
    call ffa_cmdlinearg_getbool ('allVelVars',allVel)
    call ffa_cmdlinearg_getbool ('allAccVars',allAcc)
    call ffa_cmdlinearg_getbool ('allForceVars',allForce)

    do idIn = 1, nJoints

       call NullifyJoint (joints(idIn))
       if (.not. iuSetPosAtNextEntry(infp,'&MASTERSLAVEJOINT')) then
          err = err - 1
          call ReportInputError ('MASTERSLAVEJOINT',idIn)
          cycle
       end if

       !! Default values
       id=0; extId=0; extDescr=''; type=0; version=0
       slaveId=0; masterId=0; nJTriads=1
       initPosInGlobal=0.0_dp; nJointVars=0; JointVarDefs=0; BC=1
       JVarInitVal=0.0_dp; JVarInitVel=0.0_dp; JVarInitAcc=0.0_dp
       springId=0; damperId=0; frictionSetId=0
       frictionSpringId=0; frictionEngineId=0; saveVar=0
       camThickness=-1.0_dp
       rotParam='FOLLOWER_AXIS'; tranSpringCpl='NONE'; rotSpringCpl='NONE'

       read(infp,nml=MASTERSLAVEJOINT,iostat=stat)
       if (stat /= 0) then
          err = err - 1
          call ReportInputError ('MASTERSLAVEJOINT',idIn)
          cycle
       end if

       if (nJointVars > MAX_JVAR) then
          err = -nJointVars
          call ReportInputError ('MASTERSLAVEJOINT',idIn, &
               msg='Too many joint variables')
          return ! Abort parsing due to possible overwrite of local array
       end if

       if (ignoreIC) then
          JVarInitVel = 0.0_dp
          JVarInitAcc = 0.0_dp
       end if

       joint => joints(idIn)
       call initId (joint%id,id,extId,extDescr,stat)
       joint%type = type
       joint%version = version

       if (type == FREE_p .and. version < 0) then
          call reportError (note_p,'Rotation to translation coupling '// &
               &            'set to zero for '//getJointId(joint))
       end if

       if (type >= PRISMATIC_p .and. type <= CAM_p) then
          joint%nJointDOFs = nJointVars + 1 ! Account for additional slider dof
       else
          joint%nJointDOFs = nJointVars
       end if

       !! Check rotation formulation and spring inter-connection
       joint%rotParam = iuCharToInt(rotParam,rotParamTypes_p,stat)
       springCpl(1) = iuCharToInt(tranSpringCpl,cplTypes_p,stat)
       springCpl(2) = iuCharToInt(rotSpringCpl,cplTypes_p,stat)
       if (stat < 0) then
          err = err - 1
          call ReportInputError ('MASTERSLAVEJOINT',idIn,joint%id)
       end if

       !! Check if this joint has any joint springs
       haveSpringInterConn = .false.
       if (any(springId(1:joint%nJointDOFs) > 0)) then
          allocate(joint%springEl,STAT=stat)
          if (stat /= 0) then
             err = AllocationError('ReadJoints 14')
             return
          end if
          call NullifySpring (joint%springEl,joint%id,joint%nJointDOFs)
          call AllocJointSpring (joint%springEl,joint%nJointDOFs,stat)
          if (stat /= 0) then
             err = AllocationError('ReadJoints 15')
             return
          end if
       end if

       !! Check if this joint has any friction springs
       if (any(frictionSpringId(1:joint%nJointDOFs) > 0)) then
          allocate(joint%sprFric,STAT=stat)
          if (stat /= 0) then
             err = AllocationError('ReadJoints 16')
             return
          end if
          call NullifySpring (joint%sprFric,joint%id,joint%nJointDOFs)
          allocate(joint%sprFric%spr(joint%nJointDOFs),STAT=stat)
          if (stat /= 0) then
             err = AllocationError('ReadJoints 17')
             return
          end if
       end if

       !! Dependent triad, assumed to be aligned with joint direction initially
       joint%STriad => GetPtrToId(triads,slaveId)
       if (.not. associated(joint%STriad)) then
          err = err - 1
          call ReportInputError ('MASTERSLAVEJOINT',id=joint%id, &
               msg='Invalid dependent triad Id')
          cycle
       else if (joint%STriad%nDOFs < 1) then
          err = err - 1
          call ReportInputError ('MASTERSLAVEJOINT',id=joint%id, &
               msg='The dependent triad is grounded')
       else if (joint%STriad%dependent) then
          err = err - 1
          call ReportInputError ('MASTERSLAVEJOINT',id=joint%id, &
               msg='The dependent triad is used in another joint')
       else
          joint%STriad%dependent = .true.
          if (noFEnodeAtJoints .and. joint%STriad%isFEnode) then
             call reportError (note_p,'Deactivating sectional force output '// &
                  &            'for dependent Triad'//getId(joint%STriad%id))
             joint%STriad%isFEnode = .false.
          end if
       end if

       if (type == AXIAL_p) then
          nSDOFs = 1 ! Only one dependent DOF in axial joints
          maxJTD = 3 ! Only translational independent DOFs are included
       else
          nSDOFs = joint%STriad%nDOFs
          maxJTD = 6
       end if

       if (type >= PRISMATIC_p .and. type <= CAM_p) then
          !! Handle slider DOF and glider curve for poinit-to-path joints
          allocate(joint%slider,STAT=stat)
          if (stat /= 0) then
             err = AllocationError('ReadJoints 18')
             return
          end if
          joint%slider%master => GetPtrToId(gCurves,masterId)
          if (.not. associated(joint%slider%master)) then
             err = err - 1
             call ReportInputError ('MASTERSLAVEJOINT',id=joint%id, &
                  &             msg='Invalid master curve Id')
             cycle
          end if
          nJTriads = size(joint%slider%master%CPoints)
          if (type < CAM_p .and. version > 2) then
             if (nJTriads == 3) then
                call reportError (note_p,'Using quadratic interpolation for '//&
                     &            getJointId(joint))
             else if (nJTriads > 3) then
                call reportError (note_p,'Using cubic interpolation for '// &
                     &            getJointId(joint))
             end if
          end if
       end if

       !! Allocate arrays of joint triads, joint DOFs and dependent DOFs
       allocate(joint%JMTriads(nJTriads), &
            &   joint%jointDofs(joint%nJointDOFs), &
            &   joint%BC(joint%nJointDOFs), &
            &   joint%slaveDofs(nSDOFs), STAT=stat)
       if (stat /= 0) then
          err = AllocationError('ReadJoints 20')
          return
       end if

       if (needAddBC) then
          joints(idIn)%BC = BC(1:joint%nJointDOFs)
       else
          joints(idIn)%BC = 1
          !! We still want the always fixed BCs (if any)
          do i = 1, joint%nJointDOFs
             if (BC(i) <= 0) joints(idIn)%BC(i) = 0
          end do
       end if

       !! Set pointers to joint triads
       if (associated(joint%slider)) then
          do i = 1, nJTriads
             joint%JMTriads(i)%triad => joint%slider%master%CPoints(i)%triad
          end do
       else
          joint%JMTriads(1)%triad => GetPtrToId(triads,masterId)
          if (.not. associated(joint%JMTriads(1)%triad)) then
             err = err - 1
             call ReportInputError ('MASTERSLAVEJOINT',id=joint%id, &
                  &             msg='Invalid master triad Id')
             !! Must clean up to avoid crash on dangling pointers later
             deallocate(joint%JMTriads)
             nullify(joint%JMTriads)
             cycle
          end if
       end if

       if (nJTriads == 1) then
          !! If this joint connects to beam elements, let the joint triads
          !! share the running coordinate along the beam, such that beam
          !! diagram plots can be generated for the entire beam structure
          if ( associated(joint%STriad%scoord) .and. &
               associated(joint%JMTriads(1)%triad%scoord) ) then
             deallocate(joint%STriad%scoord)
             joint%STriad%scoord => joint%JMTriads(1)%triad%scoord
          end if
       end if

       !! Initialize the joint DOFs
       lerr = err
       do iDof = 1, nJointVars

          call InitJointDof (joint%jointDofs(iDof),JointVarDefs(1,iDof), &
               &             JointVarDefs(2,iDof),iDof,joint%nJointDOFs, &
               &             jVarInitVal(iDof),jVarInitVel(iDof), &
               &             jVarInitAcc(iDof),springId(iDof),damperId(iDof), &
               &             frictionSetId(iDof),frictionSpringId(iDof), &
               &             frictionEngineId(iDof),saveVar(:,iDof),stat)

          if (stat > 0) then
             err = err - 1
             call ReportInputError ('MASTERSLAVEJOINT',id=joint%id, &
                  msg='Invalid joint spring-, damper- or friction Id')
          else if (stat < 0) then
             err = stat
             return
          end if
       end do
       if (err < lerr) cycle

       if (associated(joint%slider)) then

          iDof = joint%nJointDOFs
          joint%slider%slideDof  => joint%jointDofs(iDof)
          joint%slider%tangent   =  0.0_dp
          joint%slider%rotGrad   =  0.0_dp
          joint%slider%thickness =  camThickness

          call InitJointDof (joint%slider%slideDof,7,0,iDof,joint%nJointDOFs, &
               &             jVarInitVal(iDof),jVarInitVel(iDof), &
               &             jVarInitAcc(iDof),springId(iDof),damperId(iDof), &
               &             frictionSetId(iDof),frictionSpringId(iDof), &
               &             frictionEngineId(iDof),saveVar(:,iDof),stat)

          if (stat > 0) then
             err = err - 1
             call ReportInputError ('MASTERSLAVEJOINT',id=joint%id, &
                  msg='Invalid slider spring-, damper- or friction Id')
             cycle
          else if (stat < 0) then
             err = stat
             return
          end if

          !! Allocate array for constraint coefficients for the joint triads
          allocate(joint%slider%coeff(nJTriads), STAT=stat)
          if (stat /= 0) then
             err = AllocationError('ReadJoints 30')
             return
          end if
          joint%slider%coeff = 0.0_dp

       end if

       if (associated(joint%springEl)) then

          !! Initialize the spring element for the joint
          compFailure = .false.
          do i = 1, size(joint%springEl%spr)

             !! The ordinary joint springs, nullify or point to actual spring
             joint%springEl%spr(i)%p => joint%jointDofs(i)%spring

             !! Check if the entire joint should fail simultaneously
             if (associated(joint%springEl%spr(i)%p)) then
                if (associated(joint%springEl%spr(i)%p%failure)) then
                   if (joint%springEl%spr(i)%p%failure%compFailure) then
                      compFailure = .true.
                   end if
                end if
             end if

          end do
          if (compFailure) then
             allocate(joint%springEl%hasFailed, STAT=stat)
             if (stat /= 0) then
                err = AllocationError('ReadJoints 35')
                return
             end if
             joint%springEl%hasFailed = .false.
          end if

       end if
       if (haveSpringInterConn) then

          !! Set up matrix of spring interconnectivity
          !! Note the ldof / idof difference

          select case (springCpl(1))
          case (XY_p)
             call InterConnect (1,2,JointVarDefs(1,:),joint%springEl%lConn)
          case (YZ_p)
             call InterConnect (2,3,JointVarDefs(1,:),joint%springEl%lConn)
          case (ZX_p)
             call InterConnect (3,1,JointVarDefs(1,:),joint%springEl%lConn)
          case (XYZ_p)
             call InterConnect (1,2,JointVarDefs(1,:),joint%springEl%lConn)
             call InterConnect (2,3,JointVarDefs(1,:),joint%springEl%lConn)
             call InterConnect (3,1,JointVarDefs(1,:),joint%springEl%lConn)
          end select

          select case (springCpl(2))
          case (XY_p)
             call InterConnect (4,5,JointVarDefs(1,:),joint%springEl%lConn)
          case (YZ_p)
             call InterConnect (5,6,JointVarDefs(1,:),joint%springEl%lConn)
          case (ZX_p)
             call InterConnect (6,4,JointVarDefs(1,:),joint%springEl%lConn)
          case (XYZ_p)
             call InterConnect (4,5,JointVarDefs(1,:),joint%springEl%lConn)
             call InterConnect (5,6,JointVarDefs(1,:),joint%springEl%lConn)
             call InterConnect (6,4,JointVarDefs(1,:),joint%springEl%lConn)
          end select

          stat = 0
          do i = 1, size(joint%springEl%lConn,1)
             do j = 1, size(joint%springEl%lConn,2)
                if (joint%springEl%lConn(i,j)) then
                   if (.not. associated(joint%springEl%spr(j)%p)) stat = stat+1
                end if
             end do
          end do
          if (stat > 0) then
             call reportError (error_p,'Missing spring definitions for Joint'//&
                  &            trim(getId(joint%id))//' with interconnections')
             err = err - stat
          end if

       end if

       if (associated(joint%sprFric)) then

          !! Initialize the friction spring element for the joint
          do i = 1, size(joint%sprFric%spr)
             !! Point to actual spring of the friction, or nullify
             if (associated(joint%jointDofs(i)%friction)) then
                joint%sprFric%spr(i)%p => joint%jointDofs(i)%friction%spr
             else
                nullify(joint%sprFric%spr(i)%p)
             end if
          end do

       end if

       if (joint%nJointDOFs > 0) then
          allocate(joint%dofOutputOrder(joint%nJointDOFs), STAT=stat)
          if (stat /= 0) then
             err = AllocationError('ReadJoints 36')
             return
          end if
          if (nJTriads > 1) JointVarDefs(1,joint%nJointDOFs) = 7
          call quickSortInt (JointVarDefs(1,1:joint%nJointDOFs), &
               &             joint%dofOutputOrder)
       end if

       do i = 1, nJTriads
          joint%JMTriads(i)%nDOFs = min(maxJTD,joint%JMTriads(i)%triad%nDOFs)
          if (associated(joint%slider)) then
             joint%JMTriads(i)%JPosInT => joint%slider%master%CPoints(i)%CPosInT
          else
             allocate(joint%JMTriads(i)%JPosInT(3,4),STAT=stat)
             if (stat /= 0) then
                err = AllocationError('ReadJoints 40')
                return
             end if
          end if
          if (noFEnodeAtJoints .and. joint%JMTriads(i)%triad%isFEnode) then
             call reportError (note_p,'Deactivating sectional force output '// &
                  &            'for Triad'//getId(joint%JMTriads(i)%triad%id))
             joint%JMTriads(i)%triad%isFEnode = .false.
          end if
       end do

       !! Allocate the intermediate position matrices and associated arrays
       joint%nMats = max(0,maxval(JointVarDefs(2,1:nJointVars)))
       allocate(joint%PosFromJVars(3,4,joint%nMats), &
            &   joint%thetaVec(3,joint%nMats), &
            &   joint%numRot(joint%nMats), &
            &   STAT=stat)
       if (stat /= 0) then
          err = AllocationError('ReadJoints 50')
          return
       end if

       !! Set the initial value of the position matrices
       joint%PosFromJVars = 0.0_dp
       joint%thetaVec     = 0.0_dp
       joint%numRot       = 0
       do iDof = 1, nJointVars
          ldof = joint%jointDofs(iDof)%ldof
          iMat = joint%jointDofs(iDof)%iMat
          if (ldof <= 3) then
             joint%PosFromJVars(ldof,4,iMat) = joint%jointDofs(iDof)%jVar(1)
          else
             joint%thetaVec(ldof-3,iMat) = joint%jointDofs(iDof)%jVar(1)
             joint%numRot(iMat)          = joint%numRot(iMat) + 1
          end if
       end do
       do j = 1, joint%nMats
          call vec_to_mat (joint%thetaVec(:,j),joint%PosFromJVars(:,1:3,j))
       end do

       !! Position of Joint in Global system
       PJG = transpose(initPosInGlobal) ! due to reading row by row
       call orthonorm3 (PJG(:,1:3))     ! due to possible round-off
       joint%JPosInG = PJG              ! Initialize to input joint position

       !! Initial position of the dependent joint triad relative to Joint
       joint%STPos0InJ = matmul34(invert34(PJG),joint%STriad%ur)
       do j = joint%nMats, 1, -1
          PJGinv = invert34(joint%PosFromJVars(:,:,j))
          joint%STPos0InJ = matmul34(PJGinv,joint%STPos0InJ)
       end do

       !! Set the eccentricity vector to the dependent triad to zero
       if (nullifyEvec) joint%STPos0InJ(:,4) = 0.0_dp

       !! Make the dependent triad have system directions along the joint
       !! orientation. Not necessary with new joint formulation,
       !! except for axial joints.
       if (type == AXIAL_p .or. depDirInJointDir) then
          if (.not. associated(joint%STriad%sysDirInG)) then
             joint%STriad%sysDirInG => joint%JPosInG(:,1:3)
          else if (type == AXIAL_p) then
             call reportError (error_p,'Dependent Triad'// &
                  &            trim(getId(joint%STriad%id))//' of '// &
                  &            trim(getJointId(joint))// &
                  &            ' already has local system directions')
             err = err -1
          end if
       end if

       if (nJTriads == 1) then
          !! Position of Joint in joint Triad system
          PTGinv = invert34(joint%JMTriads(1)%triad%ur)
          joint%JMTriads(1)%JPosInT = matmul34(PTGinv,joint%JPosInG)
       end if

    end do

    !! Consistency check: Deactivate joint dampers that are not referred
    stat = 0
    do i = 1, size(dmpElms)
       if (.not. associated(dmpElms(i)%triad1) .and. dmpElms(i)%nDOFs < 1) then
          stat = stat + 1
          dmpElms(i)%dmp%isActive = .false.
          call reportError (warningFileOnly_p, &
               'Joint Damper'//trim(getId(dmpElms(i)%dmp%id))// &
               ' is not referenced by any joint.')
       end if
    end do
    if (stat > 0) then
       call reportError (warning_p,'The mechanism may be inconsistent.', &
            &            'Some joint springs/dampers are not referenced', &
            &            'and are ignored by the Dynamics Solver', &
            &            addString='ReadJoints')
    end if

900 if (err < 0) call reportError (debugFileOnly_p,'ReadJoints')

  contains

    !> @brief Initializes the inter-connected spring data for a joint.
    subroutine AllocJointSpring (joint,nDOFs,jerr)

      use SpringTypeModule, only : SpringType

      type(SpringType), intent(inout)  :: joint
      integer         , intent(in)     :: nDOFs
      integer         , intent(out)    :: jerr

      allocate(joint%spr(nDOFs),STAT=jerr)
      if (jerr /= 0 .or. all(springCpl == NONE_p)) return

      allocate(joint%lConn(2,nDOFs),joint%stiffMat(nDOFs,nDOFs),STAT=jerr)
      if (jerr /= 0) return

      joint%lConn = .false.
      joint%stiffMat = 0.0_dp
      haveSpringInterConn = .true.

    end subroutine AllocJointSpring


    !> @brief Defines a spring coupling between two joint DOFs.
    subroutine InterConnect (iDof,jDof,lDofArr,lConn)

      integer, intent(in)  :: iDof, jDof, lDofArr(:)
      logical, intent(out) :: lConn(:,:)

      do i = 1, size(lConn,2)
         if (lDofArr(i) == iDof) then
            lConn(1+(iDof-1)/3,i) = .true.
         else if (lDofArr(i) == jDof) then
            lConn(1+(jDof-1)/3,i) = .true.
         end if
      end do

    end subroutine InterConnect


    !> @brief Initializes a joint DOF.
    subroutine InitJointDof (jointDof,lDof,iMat,iDof,nJointDOFs, &
         &                   jVarInitVal,jVarInitVel,jVarInitAcc, &
         &                   springId,damperId,frictionId,frictionSpringId, &
         &                   frictionEngineId,saveVar,ierr)

      use MasterSlaveJointTypeModule, only : JointDofType, nullifyJointDof
      use FrictionTypeModule        , only : InitializeFriction
      use FrictionTypeModule        , only : BALL_JNT_FRICTION_p

      type(JointDofType), intent(out), target :: jointDof

      integer , intent(in)  :: lDof, iMat, iDof, nJointDOFs
      real(dp), intent(in)  :: jVarInitVal, jVarInitVel, jVarInitAcc
      integer , intent(in)  :: springId, damperId, frictionId, frictionSpringId
      integer , intent(in)  :: frictionEngineId, saveVar(:)
      integer , intent(out) :: ierr

      integer :: dIdx

      ierr = 0
      call nullifyJointDof (jointDof)
      jointDof%iMat    = iMat
      jointDof%lDof    = lDof
      jointDof%jVar(1) = jVarInitVal
      jointDof%jVar(2) = jVarInitVel
      jointDof%jVar(3) = jVarInitAcc

      !! Check if this joint dof is fixed. In that case, ignore all saveVar
      !! settings since then it will be no non-zero results for this DOF at all.
      !! And it won't make sense with any spring, damper or friction, either.
      if (joint%BC(iDof) == 0) then
         !! This joint dof is fixed. Only the reaction force is non-zero.
         jointDof%fixed = .true.
         jointDof%saveVar(4) = saveVar(4) > 0 .or.allSec.or.allJoint.or.allForce
         return
      end if

      if (springId > 0) then

         !! Set the spring between joint position and the dependent triad
         jointDof%spring => GetPtrToId(springs,springId)
         if (.not. associated(jointDof%spring)) then
            ierr = 1
            return
         else if (jointDof%spring%dof >= 0) then
            ierr = 2
            return
         end if

         jointDof%spring%dof = iDof
         !! Connect the spring length to the joint variable
         jointDof%spring%length => jointDof%jVar(1)

      end if
      if (damperId > 0) then

         !! Set the damper between joint position and the dependent triad
         jointDof%damper => GetPtrToId(dampers,damperId,dIdx)
         if (.not. associated(jointDof%damper)) then
            ierr = 3
            return
         else if (.not. associated(dmpElms(dIdx)%dmp,jointDof%damper)) then
            ierr = 4
            return
         end if

         jointDof%damper%isActive = .true.
         jointDof%damper%dof = iDof
         dmpElms(dIdx)%nDOFs = nJointDOFs

         !! Connect the damper length to joint variable
         jointDof%damper%length => jointDof%jVar(1)

         if (associated(jointDof%damper%spr)) then
            if (associated(jointDof%spring)) then
               !! Connect to associated spring for deformational damping
               jointDof%damper%spr => jointDof%spring
            else
               nullify(jointDof%damper%spr)
               call reportError (note_p, &
                    'Switching off deformational damping for Joint Damper'// &
                    trim(getId(jointDof%damper%id))// &
                    ' since no associated spring')
            end if
         end if

         if (associated(jointDof%damper%spr)) then
            !! Allocate velocity variable
            allocate(jointDof%damper%velocity,STAT=stat)
            if (stat /= 0) then
               ierr = AllocationError('ReadJoints 60')
               return
            end if
         else
            !! Connect the damper velocity to the joint variable
            !! for non-deformational dampers
            jointDof%damper%velocity => jointDof%jVar(2)
         end if

      end if

      if (frictionId > 0) then

         !! Set the friction between joint position and the dependent triad
         call InitializeFriction (jointDof%friction,frictionSets, &
              &                   frictionId,saveVar(4:5),ierr)
         if (ierr /= 0) return

         if (jointDof%friction%param%type == BALL_JNT_FRICTION_p) then
            if (.not. associated(joint%multiDofFriction)) then
               call InitializeFriction (joint%multiDofFriction,frictionSets, &
                    &                   frictionId,saveVar(4:5),ierr)
               if (ierr /= 0) return
            end if
         end if

         if (frictionSpringId > 0) then

            !! Set the spring between joint position and the dependent triad
            jointDof%friction%spr => GetPtrToId(springs,frictionSpringId)
            if (.not. associated(jointDof%friction%spr)) then
               ierr = 1
               return
            else if (jointDof%friction%spr%dof >= 0) then
               ierr = 2
               return
            end if

            jointDof%friction%spr%dof    = iDof
            !! Connect the spring length to the joint variable
            jointDof%friction%spr%length => jointDof%jVar(1)
            !! Ensure zero deflection initially
            jointDof%friction%spr%l0     =  jointDof%jVar(1)

         end if

         if (frictionEngineId > 0) then
            jointDof%friction%eng => GetPtrToId(engines,frictionEngineId)
         end if

      end if

      !! Determine which joint variables will be saved. Settings in the
      !! solver input file are overruled by command-line arguments, if any.
      do j = 1, 3
         jointDof%saveVar(j) = saveVar(j) > 0 .or. allSec.or.allRest.or.allJoint
      end do
      if (allLen) jointDof%saveVar(1) = .true.
      if (allVel) jointDof%saveVar(2) = .true.
      if (allAcc) jointDof%saveVar(3) = .true.

    end subroutine InitJointDof

  end subroutine ReadJoints


  !!============================================================================
  !> @brief Initializes higher pairs with data from the solver input file.
  !>
  !> @param[in] infp File unit number for the solver input file
  !> @param[in] joints Array of all joint objects in the model
  !> @param higherPairs Array of all higher pair objects in the model
  !> @param[out] err Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 7 Dec 2000

  subroutine ReadHigherPairs (infp,joints,higherPairs,err)

    use IdTypeModule              , only : ldesc_p, initId, ReportInputError
    use MasterSlaveJointTypeModule, only : MasterSlaveJointType, GetPtrToId, dp
    use MasterSlaveJointTypeModule, only : HigherPairType, DeallocateHigherPairs
    use inputUtilities            , only : iuGetNumberOfEntries
    use inputUtilities            , only : iuSetPosAtNextEntry
    use progressModule            , only : lterm
    use reportErrorModule         , only : AllocationError
    use reportErrorModule         , only : reportError, debugFileOnly_p

    integer                   , intent(in)  :: infp
    type(MasterSlaveJointType), intent(in)  :: joints(:)
    type(HigherPairType)      , pointer     :: higherPairs(:)
    integer                   , intent(out) :: err

    !! Local variables
    integer :: idIn, nHPs, stat

    !! Define the HIGHER_PAIR namelist
    integer :: id, extId(10)
    integer :: slaveJoint, slaveJointDof, masterJoint, masterJointDof
    real(dp):: coeff
    character(ldesc_p) :: extDescr
    namelist /HIGHER_PAIR/ id, extId, extDescr, &
         &    slaveJoint, slaveJointDof, masterJoint, masterJointDof, coeff

    !! --- Logic section ---

    nHPs = iuGetNumberOfEntries(infp,'&HIGHER_PAIR',err)
    if (err /= 0) return
    write(lterm,*) 'Number of &HIGHER_PAIR =',nHPs

    call DeallocateHigherPairs (higherPairs)
    allocate(higherPairs(nHPs),STAT=stat)
    if (stat /= 0) then
       err = AllocationError('ReadHigherPairs')
       return
    end if

    do idIn = 1, nHPs

       if (.not. iuSetPosAtNextEntry(infp,'&HIGHER_PAIR')) then
          err = err - 1
          call ReportInputError ('HIGHER_PAIR',idIn)
          cycle
       end if

       !! Default values
       id=0; extId=0; extDescr=''; coeff=0.0_dp
       slaveJoint=0; slaveJointDof=-1; masterJoint=0; masterJointDof=-1
       read(infp,nml=HIGHER_PAIR,iostat=stat)
       if (stat /= 0) then
          err = err - 1
          call ReportInputError ('HIGHER_PAIR',idIn)
          cycle
       end if

       call initId (higherPairs(idIn)%id,id,extId,extDescr,stat)

       higherPairs(idIn)%coeff = coeff
       higherPairs(idIn)%slaveJointDof  = slaveJointDof
       higherPairs(idIn)%masterJointDof = masterJointDof
       higherPairs(idIn)%slaveJoint  => GetPtrToId(joints,slaveJoint)
       higherPairs(idIn)%masterJoint => GetPtrToId(joints,masterJoint)
       if (.not. associated(higherPairs(idIn)%slaveJoint)) then
          err = err - 1
          call ReportInputError ('HIGHER_PAIR',id=higherPairs(idIn)%id, &
               &                 msg='Invalid joint Id')
       end if
       if (.not. associated(higherPairs(idIn)%masterJoint)) then
          err = err - 1
          call ReportInputError ('HIGHER_PAIR',id=higherPairs(idIn)%id, &
               &                 msg='Invalid joint Id')
       end if
       if (err == 0) then
          call InitiateHP (higherPairs(idIn)%slaveJoint,slaveJointDof, &
               &           higherPairs(idIn)%masterJoint,masterJointDof)
       end if

    end do

    if (err < 0) call reportError (debugFileOnly_p,'ReadHigherPairs')

  contains

    !> @brief Computes initial values for joint variables via a higher pair.
    subroutine InitiateHP (out,oDof,inp,iDof)

      type(MasterSlaveJointType), intent(inout) :: out
      type(MasterSlaveJointType), intent(in)    :: inp
      integer                   , intent(in)    :: oDof, iDof

      !! Initialize the output joint velocity and acceleration from that of
      !! the input joint. Note that the transmission ratio is not applied here.
      !! That has to be done after the initial conditions on the joint variables
      !! have been resolved, in InitiateHigherPairs.
      if (oDof > 0 .and. oDof <= out%nJointDOFs) then
         if (iDof > 0 .and. iDof <= inp%nJointDOFs) then
            out%jointDofs(oDof)%jVar(2:3) = inp%jointDofs(iDof)%jVar(2:3)
         end if
      end if

    end subroutine InitiateHP

  end subroutine ReadHigherPairs


  !!============================================================================
  !> @brief Computes constraint coefficient and joint variable start values.
  !>
  !> @param joints Array of all joint objects in the model
  !> @param[in] forces Array of all external point load objects in the model
  !> @param[in] motions Array of all prescribed motion objects in the model
  !> @param[out] mpreac Matrix of pointers to reaction forces
  !> @param[in] RF System reaction forces associated with constrained DOFs
  !> @param[out] ierr Error flag
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 24 Oct 2002

  subroutine InitiateJoints (joints,forces,motions,mpreac,RF,ierr)

    use ForceTypeModule               , only : ForceType
    use MotionTypeModule              , only : MotionType
    use MasterSlaveJointTypeModule    , only : MasterSlaveJointType, dp
    use MasterSlaveJointRoutinesModule, only : moveJointInPosition
    use MasterSlaveJointRoutinesModule, only : updateDependencies
    use MasterSlaveJointRoutinesModule, only : initiateJointVars
#ifdef FT_DEBUG
    use dbgUnitsModule                , only : dbgJoint, dbgCurve
    use fileUtilitiesModule           , only : getDBGfile
#endif
    use reportErrorModule             , only : reportError, debugFileOnly_p
    use reportErrorModule             , only : allocationError

    type(MasterSlaveJointType), intent(inout) :: joints(:)
    type(ForceType),   target , intent(in)    :: forces(:)
    type(MotionType)          , intent(in)    :: motions(:)
    integer,                    pointer       :: mpreac(:)
    real(dp),          target , intent(in)    :: RF(:)
    integer                   , intent(out)   :: ierr

    !! Local variables
    logical :: haveFriction
    integer :: i, j, ipl, ipr, err

    !! --- Logic section ---

#ifdef FT_DEBUG
    if (size(joints) > 0) dbgJoint = getDBGfile(6,'joint.dbg')
#endif

    ierr = 0
    do i = 1, size(joints)

#ifdef FT_DEBUG
       if (size(joints(i)%JMTriads) > 2) then
          dbgCurve = getDBGfile(7,'gliderCurve.dbg')
       end if
#endif

       !! Assign pointers to the reaction force or applied load for each dof
       haveFriction = .false.
       do j = 1, joints(i)%nJointDOFs
          ipl = joints(i)%jointDofs(j)%loadIdx
          ipr = joints(i)%jointDofs(j)%sysDOF
          if (.not.associated(mpreac)) then
             ipr = 0
          else if (ipr <= size(mpreac)) then
             ipr = mpreac(ipr)
          else
             ipr = 0
          end if
          if (ipr > 0 .and. ipr <= size(RF)) then
             joints(i)%jointDofs(j)%F => RF(ipr)
          else if (ipl > 0 .and. ipl <= size(forces)) then
             joints(i)%jointDofs(j)%F => forces(ipl)%F
          end if
          if (associated(joints(i)%jointDofs(j)%friction)) haveFriction = .true.
       end do

       call moveJointInPosition (joints(i),err)
       if (err < 0) ierr = ierr + err

       call updateDependencies (joints(i),motions)
       call initiateJointVars (joints(i),err)
       if (err < 0) ierr = ierr + err

       !! Ensure the nodeForce array is allocated for joints with friction
       if (haveFriction .and. .not.associated(joints(i)%STriad%nodeForce)) then
          allocate(joints(i)%STriad%nodeForce(joints(i)%STriad%nDOFs),STAT=err)
          if (err /= 0) then
             ierr = allocationError('InitiateJoints')
             return
          end if
       end if

    end do

    if (ierr < 0) call reportError (debugFileOnly_p,'InitiateJoints')

  end subroutine InitiateJoints


  !!============================================================================
  !> @brief Computes initial values for joint variables via higher pairs.
  !>
  !> @param higherPairs Array of all higher pair objects in the model
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 2 Oct 2018

  subroutine InitiateHigherPairs (higherPairs)

    use MasterSlaveJointTypeModule, only : HigherPairType
    use MasterSlaveJointTypeModule, only : MasterSlaveJointType, dp

    type(HigherPairType), pointer :: higherPairs(:)

    !! Local variables
    integer :: i

    !! --- Logic section ---

    do i = 1, size(higherPairs)
       call InitiateHP (higherPairs(i)%slaveJoint, &
            &           higherPairs(i)%slaveJointDof, &
            &           higherPairs(i)%masterJoint, &
            &           higherPairs(i)%masterJointDof, &
            &           higherPairs(i)%coeff)
    end do

  contains

    !> @brief Computes initial values for joint variables via a higher pair.
    subroutine InitiateHP (out,oDof,inp,iDof,c)

      type(MasterSlaveJointType), intent(inout) :: out
      type(MasterSlaveJointType), intent(in)    :: inp
      integer                   , intent(in)    :: oDof, iDof
      real(dp)                  , intent(in)    :: c

      !! Initialize the output joint velocity and acceleration from the
      !! corresponding values of the input joint, taking into account
      !! the transmission ratio (c).
      if (oDof > 0 .and. oDof <= out%nJointDOFs) then
         if (iDof > 0 .and. iDof <= inp%nJointDOFs) then
            out%jointDofs(oDof)%jVar(2:3) = c*inp%jointDofs(iDof)%jVar(2:3)
         end if
      end if

    end subroutine InitiateHP

  end subroutine InitiateHigherPairs

end module initiateJointTypeModule
