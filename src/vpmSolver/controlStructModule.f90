!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file controlStructModule.f90
!> @brief Coupled control system and structure modal analysis.

!!==============================================================================
!> @brief Module with data types and subroutines for coupled modal analysis.
!>
!> @details This module contains some data types and associated subroutines for
!> conducting coupled modal analysis of mechanim models with control systems.

module ControlStructModule

  use ControlTypeModule, only : CtrlPrm
  use SensorTypeModule , only : SensorType, dp
  use ForceTypeModule  , only : ForceType

  implicit none


  !> @brief Data type representing a control input gradient.
  type StructSensorGradType
     !> Which control input element this sensor is connected to.
     !> - @a j_index i.e. second (column) index in Jacobi_CinToCout(:,:)
     !> - @a i_index i.e. first (row) index in Jacobi_SensorToCIn(:,:)
     integer :: iCin
     !> Which VREG this is connected to (not necessarily equal to j_index)
     integer :: whichVreg

     integer :: lNode(2)    !< Local node number for the element matrices
     integer :: dofStart(2) !< Local dof start for this node
     integer :: nDofs(2)    !< number of dofs for this node

     !> Pointer to the control parameter which has structural input
     type(CtrlPrm)   , pointer :: pCtrlPrm
     type(SensorType), pointer :: pSensor  !< Pointer to the actual sensor
  end type StructSensorGradType


  !> @brief Data type representing a control output gradient.
  type ForceControlGradType
     !> Which control output element this force is connected to
     !> - @a i_index (first)  index in Jacobi_CinToCout(:,:)
     !> - @a j_index (second) index in Jacobi_ForceToCout(:,:)
     integer :: iCout
     !> which VREG this is connected to (not necessarily equal to j_index)
     integer :: whichVreg

     integer :: lNode    !< Local node number for the element matrices
     integer :: dofStart !< Local dof start for this node
     integer :: nDofs    !< number of dofs for this node

     !> Pointer to the force with the control output as force magnitude
     type(ForceType), pointer :: pForce
  end type ForceControlGradType


  !> @brief Data type representing the control-system/structure coupling.
  type ControlStructType

     integer :: ctrlSysEigFlag !< Flag telling which perturbation method to use

     integer :: samElNum !< Element number for SAM reference
     integer :: nDOFs    !< Total number of DOFs in this element

     integer, pointer :: samMNPC(:) !< Matrix of Nodal Point Correspondance

     !! CONTROL PERTURBATION

     integer, pointer :: whichVregIn(:)  !< vreg indices for structural sensors
     integer, pointer :: whichVregOut(:) !< vreg indices for structural loads

     type(StructSensorGradType), pointer :: structSensors(:) !< Sensor side
     type(ForceControlGradType), pointer :: controlForces(:) !< Force side


     !! The tangent matrices.
     !! These must have a lookup table to internal DOFs from both the input side
     !! (sensor) and output side (force) where the case of the input being equal
     !! to the output is properly handled (i.e., they have the same DOF).
     !! The dimension of the matrices are "Number of unique inputs and outputs".

     !> @brief Table of controller properties.
     !> @details Dimension:
     !> (no. of outputs from controller) *
     !> (no. of controller properties) *
     !> (no. of inputs to controller)
     !> - @a ctrlProps(:,1,:) = @b Q
     !> - @a ctrlProps(:,2,:) = @b K
     !> - @a ctrlProps(:,3,:) = @b C
     !> - @a ctrlProps(:,4,:) = @b M
     real(dp), pointer :: ctrlProps(:,:,:)

     !! System equation: M*(d2x/dt2) + C*(dx/dt) + K*X + Q*int(x)dt = F
     real(dp), pointer :: massMat(:,:)  !< mass matrix @b M
     real(dp), pointer :: dampMat(:,:)  !< damping matrix @b C
     real(dp), pointer :: stiffMat(:,:) !< stiffness matrix @b K
     real(dp), pointer :: SSEEMat(:,:)  !< steady-state error matrix @b Q

     !> Rows of the sensor gradient matrix
     real(dp), pointer :: Grad_CinWrtSensor(:,:) !< Dimension: (nCin,nDofs)
     !> Columns of the force gradient matrix
     real(dp), pointer :: Grad_ForceWrtCout(:,:) !< Dimension: (nDofs,nCout)

     integer :: SSEEMatIsNonSymmetric  !< If .true., the @b Q matrix is nonsym.
     integer :: stiffMatIsNonSymmetric !< If .true., the @b K matrix is nonsym.
     integer :: dampMatIsNonSymmetric  !< If .true., the @b C matrix is nonsym.
     integer :: massMatIsNonSymmetric  !< If .true., the @b M matrix is nonsym

  end type ControlStructType


  !> @brief Standard routine for writing an object to file.
  interface WriteObject
     module procedure WriteControlStructType
  end interface


contains

  !!============================================================================
  !> @brief Standard routine for writing an object to io.
  !>
  !> @param[in] pCS Data for coupled control system and structure modal analysis
  !> @param[in] io File unit number to write to
  !> @param[in] complexity Indicates the amount of print
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 19 Feb 2024

  subroutine WriteControlStructType (pCS,io,complexity)

    use IdTypeModule     , only : getId
    use manipMatrixModule, only : writeObject

    type(ControlStructType), intent(in) :: pCS
    integer                , intent(in) :: io
    integer                , intent(in) :: complexity

    !! Local variables
    integer :: i

    !! --- Logic section ---

    write(io,"(/A/A)") 'Control-Structure coupling','{'

    write(io,*) 'ctrlSysEigFlag =', pCS%ctrlSysEigFlag
    write(io,*) 'samElNum       =', pCS%samElNum
    write(io,*) 'nDOFs          =', pCS%nDOFs
    if (associated(pCS%samMNPC)) then
       write(io,*) 'samMNPC        =', pCS%samMNPC
    end if
    if (associated(pCS%whichVRegIn)) then
       write(io,*) 'whichVRegIn    =', pCS%whichVRegIn
    end if
    if (associated(pCS%whichVRegOut)) then
       write(io,*) 'whichVRegOut   =', pCS%whichVRegOut
    end if

    if (complexity > 1 .and. associated(pCS%structSensors)) then
       do i = 1, size(pCS%structSensors)
          call writeStructSensor (i,pCS%structSensors(i))
       end do
    end if
    if (complexity > 1 .and. associated(pCS%controlForces)) then
       do i = 1, size(pCS%controlForces)
          call writeForceControl (i,pCS%controlForces(i))
       end do
    end if

    if (complexity > 2 .and. associated(pCS%massMat)) then
       write(io,"()")
       call writeObject (pCS%massMat,io,'Mass matrix')
    end if
    if (complexity > 2 .and. associated(pCS%dampMat)) then
       write(io,"()")
       call writeObject (pCS%dampMat,io,'Damping matrix')
    end if
    if (complexity > 2 .and. associated(pCS%stiffMat)) then
       write(io,"()")
       call writeObject (pCS%stiffMat,io,'Stiffness matrix')
    end if
    if (complexity > 2 .and. associated(pCS%SSEEMat)) then
       write(io,"()")
       call writeObject (pCS%SSEEMat,io,'Steady-state error elimination matrix')
    end if

    write(io,"(A)") '}'

  contains

    !> @brief Writes out a controlstructmodule::structsensorgradtype object.
    subroutine writeStructSensor (idx,ssg)
      integer                   , intent(in) :: idx
      type(StructSensorGradType), intent(in) :: ssg
      write(io,"(/' Structural sensor',I3,I4)") idx, ssg%iCin
      write(io,*) ' whichVreg  =', ssg%whichVreg
      write(io,*) ' lNode      =', ssg%lNode
      write(io,*) ' dofStart   =', ssg%dofStart
      write(io,*) ' nDOFs      =', ssg%nDofs
      if (associated(ssg%pCtrlPrm)) then
         write(io,*) ' var        =',ssg%pCtrlPrm%var
         if (associated(ssg%pCtrlPrm%engine)) then
            write(io,*) ' engine(id) =', ssg%pCtrlPrm%engine%id%userId
         end if
      end if
      if (associated(ssg%pSensor)) then
         write(io,*) ' sensor(id) =', ssg%pSensor%id%userId
      end if
    end subroutine writeStructSensor

    !> @brief Writes out a controlstructmodule::forcecontrolgradtype object.
    subroutine writeForceControl (idx,fcg)
      integer                   , intent(in) :: idx
      type(ForceControlGradType), intent(in) :: fcg
      write(io,"(/' Force control',I3,I4)") idx, fcg%iCout
      write(io,*) ' whichVreg  =', fcg%whichVreg
      write(io,*) ' lNode      =', fcg%lNode
      write(io,*) ' dofStart   =', fcg%dofStart
      write(io,*) ' nDOFs      =', fcg%nDofs
      if (associated(fcg%pForce)) then
         write(io,*) ' force(id)  =', fcg%pForce%id%userId
      end if
    end subroutine writeForceControl

  end subroutine WriteControlStructType


  !!============================================================================
  !> @brief Initializes the control struct data type.
  !>

  subroutine InitiateControlStruct (pCS,inputs,triads,joints,forces, &
       &                            numVreg,numNod,ierr)

    use TriadTypeModule           , only : TriadType
    use MasterSlaveJointTypeModule, only : MasterSlaveJointType
    use SensorTypeModule          , only : TRIAD_p, RELATIVE_TRIAD_p
    use SensorTypeModule          , only : JOINT_VARIABLE_p
    use ForceTypeModule           , only : ForceType
    use ControlTypeModule         , only : getSensor
    use ReportErrorModule         , only : allocationError
    use ReportErrorModule         , only : reportError, debugFileOnly_p

    type(ControlStructType)   , intent(out)        :: pCS
    type(CtrlPrm)             , intent(in), target :: inputs(:)
    type(TriadType)           , intent(in), target :: triads(:)
    type(MasterSlaveJointType), intent(in), target :: joints(:)
    type(ForceType)           , intent(in), target :: forces(:)
    integer                   , intent(in)         :: numVreg, numNod
    integer                   , intent(out)        :: ierr

    !! Local variables

    type(TriadType)           , pointer :: pTriad1, pTriad2
    type(MasterSlaveJointType), pointer :: pJoint
    type(StructSensorGradType), pointer :: pS
    type(ForceControlGradType), pointer :: pF

    integer, allocatable :: lNode_from_SAM_node(:), local_MADOF(:)
    integer, allocatable :: iCin_from_allVreg(:), iCout_from_allVreg(:)
    integer, pointer     :: tmpIdx(:)

    integer :: i, n, iCin, iCout, iSAM, nElNodes, whichVreg
    integer :: numStructCtrlParams, numControlForces, numVregIn, numVregOut

    !! --- Logic section ---

    ierr = 0
    nullify(pCS%structSensors)
    nullify(pCS%controlForces)

    !! Find the control input parameters that have structural sensors
    !! (triad- and joint dofs only)

    call getCtrlParamsWithStructSensors (inputs,numStructCtrlParams,tmpIdx)
    if (numStructCtrlParams < 0) goto 915
    if (numStructCtrlParams == 0) return

    !! Initialize the input side

    allocate(lNode_from_SAM_node(numNod), iCin_from_AllVreg(numVreg), &
         &   pCS%structSensors(numStructCtrlParams), STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('InitiateControlStruct')
       return
    end if

    lNode_from_SAM_node = 0
    iCin_from_AllVreg  = 0


    !! SENSOR side initialization

    do i = 1, numStructCtrlParams
       pS          => pCS%structSensors(i)
       pS%pCtrlPrm => inputs(tmpIdx(i))
       pS%pSensor  => getSensor(pS%pCtrlPrm)

       pS%iCin = 0  !! i.e. not set yet
       pS%whichVreg = pS%pCtrlPrm%var
       iCin_from_AllVreg(pS%whichVreg) = 1 !! flag that this vreg is part of the cIn side

       select case (pS%pSensor%type)
       case (TRIAD_p)
          pTriad1 => triads(pS%pSensor%index(1))
          lNode_from_SAM_node(pTriad1%samNodNum) = 1

       case (RELATIVE_TRIAD_p)
          pTriad1 => triads(pS%pSensor%index(1))
          pTriad2 => triads(pS%pSensor%index(2))
          lNode_from_SAM_node(pTriad1%samNodNum) = 1
          lNode_from_SAM_node(pTriad2%samNodNum) = 1

       case (JOINT_VARIABLE_p)
          pJoint => joints(pS%pSensor%index(1))
          lNode_from_SAM_node(pJoint%samNodNum) = 1
       end select

    end do

    !! Count number of vreg with structural input
    numVregIn = 0
    do i = 1, numVreg
       if (iCin_from_AllVreg(i) > 0) then
          numVregIn = numVregIn + 1
          iCin_from_AllVreg(i) = numVregIn
       end if
    end do

    !! Set which vreg in (with structural input) are to be perturbed
    allocate(pCS%whichVregIn(numVregIn), STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('InitiateControlStruct 2')
       return
    end if

    !! Set which iCin this sensor is connected to
    do i = 1, numStructCtrlParams
       whichVreg                 = pCS%structSensors(i)%whichVreg
       iCin                      = iCin_from_AllVreg(whichVreg)
       pCS%structSensors(i)%iCin = iCin
       pCS%whichVregIn(iCin)     = whichVreg
    end do

    !! Clean up some scratch space
    deallocate(iCin_from_AllVreg,tmpIdx)

    !! FORCES side initialization
    !! Find the forces whose sensor is a controller variable

    call getCtrlOutForces (forces, numControlForces, tmpIdx)
    if (numControlForces < 0) goto 915

    allocate(pCS%controlForces(numControlForces), &
         &   iCout_from_AllVreg(numVreg), STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('InitiateControlStruct 3')
       return
    end if

    iCout_from_AllVreg = 0

    do i = 1, numControlForces
       pF         => pCS%controlForces(i)
       pF%pForce  => forces(tmpIdx(i))

       pF%iCout = 0 !! i.e. not yet set
       pF%whichVreg = pF%pForce%engine%args(1)%p%dof
       iCout_from_AllVreg(pF%whichVreg) = 1  !! flag this as one to read from

       if      (associated(pF%pForce%triad)) then

          lNode_from_SAM_node(pF%pForce%triad%samNodNum) = 1

       else if (associated(pF%pForce%joint)) then

          lNode_from_SAM_node(pF%pForce%joint%samNodNum) = 1

       else
          !! Error, actually
       end if

    end do

    !! Count number of vreg which are to be read
    numVregOut = 0
    do i = 1, numVreg
       if (iCout_from_AllVreg(i) > 0) then
          numVregOut = numVregOut + 1
          iCout_from_AllVreg(i) = numVregOut
       end if
    end do

    allocate(pCS%whichVregOut(numVregOut), STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('InitiateControlStruct 3a')
       return
    end if

    !! Set which vreg in are to be read

    do i = 1, numControlForces
       whichVreg                  = pCS%controlForces(i)%whichVreg
       iCout                      = iCout_from_AllVreg(whichVreg)
       pCS%controlForces(i)%iCout = iCout
       pCS%whichVregOut(iCout)    = whichVreg
    end do

    !! Clean up some scratch space
    deallocate(iCout_from_AllVreg,tmpIdx)

    !! Count number of element nodes
    nElNodes = count(lNode_from_SAM_node > 0)

    allocate(pCS%samMNPC(nElNodes), local_MADOF(nElNodes+1), STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('InitiateControlStruct 4')
       return
    end if

    !! Set the local element node numbers

    nelNodes = 0
    do iSam = 1, numNod
       if ( lNode_from_SAM_node(iSam) > 0 ) then
          nELNodes = nElNodes + 1
          lNode_from_SAM_node(iSam) = nELNodes
          pCS%samMNPC(nELNodes) = iSam
       end if
    end do


    !! Set the local node association and dofStart for all forces and sensors

    local_MADOF = 0
    do i = 1, numStructCtrlParams
       pS => pCS%structSensors(i)
       pS%lNode    = 0
       pS%nDOFs    = 0
       pS%dofStart = 0

       select case (pS%pSensor%type)
       case (TRIAD_p)
          pTriad1 => triads(pS%pSensor%index(1))
          pS%lNode(1) = lNode_from_SAM_node(pTriad1%samNodNum)
          pS%nDOFs(1) = pTriad1%nDOFs

       case (RELATIVE_TRIAD_p)
          pTriad1 => triads(pS%pSensor%index(1))
          pTriad2 => triads(pS%pSensor%index(2))
          pS%lNode(1) = lNode_from_SAM_node(pTriad1%samNodNum)
          pS%nDOFs(1) = pTriad1%nDOFs
          pS%lNode(2) = lNode_from_SAM_node(pTriad2%samNodNum)
          pS%nDOFs(2) = pTriad2%nDOFs

       case (JOINT_VARIABLE_p)
          pJoint => joints(pS%pSensor%index(1))
          pS%lNode(1) = lNode_from_SAM_node(pJoint%samNodNum)
          pS%nDOFs(1) = pJoint%nJointDOFs
       end select

       do n = 1, 2
          if (pS%lNode(n) > 0) local_MADOF(pS%lNode(n)) = pS%nDOFs(n)
       end do
    end do

    do i = 1, numControlForces
       pF => pCS%controlForces(i)

       if      (associated(pF%pForce%triad)) then
          pF%lNode = lNode_from_SAM_node(pF%pForce%triad%samNodNum)
          pF%nDOFs = pF%pForce%triad%nDOFs
       else if (associated(pF%pForce%joint)) then
          pF%lNode = lNode_from_SAM_node(pF%pForce%joint%samNodNum)
          pF%nDOFs = pF%pForce%joint%nJointDOFs
       else
          pF%lNode = 0
          pF%nDOFs = 0
       end if

       if (pF%lNode > 0) local_MADOF(pF%lNode) = pF%nDOFs
    end do

    !! Now accumulate the dofStart and total number of dofs for this element

    pCS%nDOFs = 0
    do i = 1, nElNodes
       n = local_MADOF(i)
       local_MADOF(i) = pCS%nDOFs + 1
       pCS%nDOFs = pCS%nDOFs + n
    end do
    local_MADOF(nElNodes+1) = pCS%nDOFs + 1

    !! Also store the dofStart in the forces and sensors

    do i = 1, numStructCtrlParams
       pS => pCS%structSensors(i)
       if (pS%lNode(1) > 0) pS%dofStart(1) = local_MADOF(pS%lNode(1))
       if (pS%lNode(2) > 0) pS%dofStart(2) = local_MADOF(pS%lNode(2))
    end do

    do i = 1, numControlForces
       pF => pCS%controlForces(i)
       pCS%controlForces(i)%dofStart = local_MADOF(pF%lNode)
    end do

    !! Clean up some scratch space
    deallocate(lNode_from_SAM_node,local_MADOF)

    allocate(pCS%Grad_CinWrtSensor(numVregIn,pCS%nDofs), &
         &   pCS%Grad_ForceWrtCout(pCS%nDofs,numVregOut), &
         &   pCS%ctrlProps(numVregOut,4,numVregIn), &
         &   pCS%massMat(pCS%nDofs,pCS%nDofs), &
         &   pCS%dampMat(pCS%nDofs,pCS%nDofs), &
         &   pCS%stiffMat(pCS%nDofs,pCS%nDofs), &
         &   pCS%SSEEMat(pCS%nDofs,pCS%nDofs), STAT=ierr)
    if (ierr /= 0) ierr = allocationError('InitiateControlStruct 5')

    return

915 ierr = -1
    call reportError (debugFileOnly_p,'InitiateControlStruct')

  end subroutine InitiateControlStruct


  !!============================================================================
  !> @brief Finds all control input parameters coupled to structural DOFs.

  subroutine getCtrlParamsWithStructSensors (ctrlParams, numStructCtrlParams, &
       &                                     ctrlParamsWithStructSensors)

    use ReportErrorModule, only : allocationError, reportError, warning_p

    type(CtrlPrm), target, intent(in)  :: ctrlParams(:)
    integer,               intent(out) :: numStructCtrlParams
    integer,      pointer, intent(out) :: ctrlParamsWithStructSensors(:)

    !! Local variables
    integer :: i, iPrm

    !! --- Logic section ---

    numStructCtrlParams = 0
    do i = 1, size(ctrlParams)
       if (hasStructSensor(ctrlParams(i),.false.)) then
          numStructCtrlParams = numStructCtrlParams + 1
       end if
    end do
    if (numStructCtrlParams == 0) then
       call reportError (warning_p,'No structural inputs in control system.', &
            'Coupled control system modal analysis is therefore switched off.')
       return
    end if

    allocate(ctrlParamsWithStructSensors(numStructCtrlParams),stat=iPrm)
    if (iPrm /= 0) then
       numStructCtrlParams = allocationError('getCtrlParamsWithStructSensors')
       return
    end if

    iPrm = 0
    do i = 1, size(ctrlParams)
       if (hasStructSensor(ctrlParams(i),.true.)) then
          iPrm = iPrm + 1
          ctrlParamsWithStructSensors(iPrm) = i
       end if
    end do

  contains

    !> @brief Checks if a control parameter is coupled to a structural DOF.
    logical function hasStructSensor (prm,notify)
      use ControlTypeModule, only : getSensor
      use SensorTypeModule , only : SensorType, sensorType_p
      use SensorTypeModule , only : TRIAD_P, RELATIVE_TRIAD_p
      use SensorTypeModule , only : JOINT_VARIABLE_p, TIME_p
      type(CtrlPrm), intent(in) :: prm
      logical      , intent(in) :: notify
      type(SensorType), pointer :: sensor
      sensor => getSensor(prm)
      if (associated(sensor)) then
         hasStructSensor = sensor%type == TRIAD_P .or. &
              &            sensor%type == RELATIVE_TRIAD_p .or. &
              &            sensor%type == JOINT_VARIABLE_p
         if (.not.hasStructSensor .and. sensor%type/=TIME_p .and. notify) then
            call reportError (warning_p,'Unsupported sensor type '// &
                 &            trim(sensorType_p(sensor%type))//' (ignored)', &
                 &            addString='getCtrlParamsWithStructSensors')
         end if
      else ! Logic error, should not happen
         hasStructSensor = .false.
      end if
    end function hasStructSensor

  end subroutine getCtrlParamsWithStructSensors


  !!============================================================================
  !> @brief Finds all control out forces.

  subroutine getCtrlOutForces (forces,numCtrlOutForces,ctrlOutForces)

    use ForceTypeModule  , only : ForceType
    use ReportErrorModule, only : allocationError, reportError, warning_p

    type(ForceType), target, intent(in)  :: forces(:)
    integer,                 intent(out) :: numCtrlOutForces
    integer,        pointer, intent(out) :: ctrlOutForces(:)

    !! Local variables
    integer :: i, iCtrl

    !! --- Logic section ---

    numCtrlOutForces = 0
    do i = 1, size(forces)
       if (isCtrlOutForce(forces(i))) then
          numCtrlOutForces = numCtrlOutForces + 1
       end if
    end do
    if (numCtrlOutForces == 0) then
       call reportError (warning_p,'No control output forces in the model')
       return
    end if

    allocate(ctrlOutForces(numCtrlOutForces),stat=iCtrl)
    if (iCtrl /= 0) then
       numCtrlOutForces = allocationError('getCtrlOutForces')
       return
    end if

    iCtrl = 0
    do i = 1, size(forces)
       if (isCtrlOutForce(forces(i))) then
          iCtrl = iCtrl + 1
          ctrlOutForces(iCtrl) = i
       end if
    end do

  contains

    !> @brief Checks if the source of a force is a control system variable.
    logical function isCtrlOutForce (force)
      use SensorTypeModule, only : CONTROL_p
      type(ForceType), intent(in) :: force
      isCtrlOutForce = .false.
      if (.not. associated(force%engine)) return
      if (.not. associated(force%engine%args)) return
      if (size(force%engine%args) < 1) return
      isCtrlOutForce = force%engine%args(1)%p%type == CONTROL_p
    end function isCtrlOutForce

  end subroutine getCtrlOutForces


  !!============================================================================
  !> @brief Computes the gradient matrices of forces w.r.t. controller inputs.
  !>

  subroutine BuildStructControlJacobi (pCS,ctrl,sys,msim,ierr)

    use ControlTypeModule   , only : ControlType
    use SystemTypeModule    , only : SystemType
    use ForceRoutinesModule , only : calcExtForceGradient
    use EngineRoutinesModule, only : SensorGradient
    use ReportErrorModule   , only : reportError, debugFileOnly_p

    type(ControlStructType), intent(inout) :: pCS
    type(ControlType)      , intent(in)    :: ctrl
    type(SystemType)       , intent(inout) :: sys
    integer                , intent(in)    :: msim(:)
    integer                , intent(out)   :: ierr

    !! Local variables
    integer  :: i, iNode, iCin, iCout, iStart, nStep, numVregIn, numVregOut
    real(dp) :: sGrad(12), fGrad(6)

    !! --- Logic section ---

    ierr = 0
    numVregIn = size(pCS%whichVregIn)
    numVregOut = size(pCS%whichVregOut)

    !! Compute the force gradients w.r.t control outputs

    pCS%Grad_ForceWrtCout = 0.0_dp

    do i = 1, size(pCS%controlForces)
       call calcExtForceGradient (pCS%controlForces(i)%pForce, fGrad, ierr)
       if (ierr < 0) goto 915

       !! Insert into Grad_ForceWrtCout
       iStart = pCS%controlForces(i)%dofStart
       iCout  = pCS%controlForces(i)%iCout
       call DAXPY (pCS%controlForces(i)%nDofs, 1.0_dp, &
            &      fGrad(1), 1, pCS%Grad_ForceWrtCout(iStart,iCout), 1)
    end do

    !! Compute the control input gradients w.r.t displacement dofs

    pCS%Grad_CinWrtSensor = 0.0_dp

    do i = 1, size(pCS%structSensors)
       call SensorGradient (pCS%structSensors(i)%pCtrlPrm%engine, sGrad, ierr)
       if (ierr < 0) goto 915

       do iNode = 1, 2
          if (pCS%structSensors(i)%lNode(iNode) > 0) then
             !! Insert into Grad_CinWrtSensor
             iStart = pCS%structSensors(i)%dofStart(iNode)
             iCin   = pCS%structSensors(i)%iCin
             call DAXPY (pCS%structSensors(i)%nDofs(iNode), 1.0_dp, &
                  &      sGrad(iNode*6-5), 1, &
                  &      pCS%Grad_CinWrtSensor(iCin,iStart), numVregIn)
          end if
       end do
    end do


    !! TODO, Magne: Add controller gradients here and build the full gradients
    pCS%ctrlProps = 0.0_dp

    select case (pCS%ctrlSysEigFlag)
    case (1) ! nPerturb = 3: P, I and D gains
       call EstimateControllerProperties01 (sys, ctrl, msim, &
            &                               pCS%whichVregIn, pCS%whichVregOut, &
            &                               pCS%ctrlProps, ierr)

    case (2) ! nPerturb = 4
       call EstimateControllerProperties02 (sys, ctrl, msim, &
            &                               pCS%whichVregIn, pCS%whichVregOut, &
            &                               pCS%ctrlProps, ierr)

    case (3) ! nPerturb = 6
       call EstimateControllerProperties03 (sys, ctrl, msim, &
            &                               pCS%whichVregIn, pCS%whichVregOut, &
            &                               pCS%ctrlProps, ierr)

    case (4:8) ! nPerturb = 5
       nStep = pCS%ctrlSysEigFlag-3 ! nStep = 1...5
       call EstimateControllerProperties04 (sys, ctrl, msim, &
            &                               pCS%whichVregIn, pCS%whichVregOut, &
            &                               nStep, pCS%ctrlProps, ierr)

    case (9:11) ! nPerturb = 5
       nStep = 10**(pCS%ctrlSysEigFlag-8) ! nStep = 10,100,1000
       call EstimateControllerProperties04 (sys, ctrl, msim, &
            &                               pCS%whichVregIn, pCS%whichVregOut, &
            &                               nStep, pCS%ctrlProps, ierr)

    case (500) ! nPerturb = 1
       call EstimateControllerProperties500 (sys, ctrl, msim, &
            &                                pCS%whichVregIn, pCS%whichVregOut, &
            &                                pCS%ctrlProps, ierr)

    case default ! Error
       if (pCS%ctrlSysEigFlag > 0) then
          ierr = -pCS%ctrlSysEigFlag
       else
          ierr = -999
       end if
    end select
    if (ierr < 0) goto 915

    nStep = ierr

    !! Steady-state error elimination matrix (Q)
    call dMMM (pCS%SSEEMat,pCS%ctrlProps(:,1,:),ierr)
    if (ierr < 0) goto 915

    !! Stiffness matrix (K)
    call dMMM (pCS%stiffMat,pCS%ctrlProps(:,2,:),ierr)
    if (ierr < 0) goto 915

    !! Damping matrix (C)
    call dMMM (pCS%dampMat,pCS%ctrlProps(:,3,:),ierr)
    if (ierr < 0) goto 915

    !! Mass matrix (M)
    call dMMM (pCS%massMat,pCS%ctrlProps(:,4,:),ierr)
    if (ierr < 0) goto 915

    !! Check for symmetry
    pCS%SSEEMatIsNonSymmetric  = isNonSymmetric(pCS%SSEEMat)
    pCS%stiffMatIsNonSymmetric = isNonSymmetric(pCS%stiffMat)
    pCS%dampMatIsNonSymmetric  = isNonSymmetric(pCS%dampMat)
    pCS%massMatIsNonSymmetric  = isNonSymmetric(pCS%massMat)

!    !! Symmetrize
!    !! TODO,Magne: Create input option on this
!    do i = 1,pCS%nDofs
!       do j = i,pCS%nDofs
!            pCS%SSEEMat(i,j)  = (pCS%SSEEMat(i,j) + pCS%SSEEMat(j,i)) * 0.5_dp
!            pCS%SSEEMat(j,i)  =  pCS%SSEEMat(i,j)
!            pCS%stiffMat(i,j) = (pCS%stiffMat(i,j) + pCS%stiffMat(j,i)) * 0.5_dp
!            pCS%stiffMat(j,i) =  pCS%stiffMat(i,j)
!            pCS%dampMat(i,j)  = (pCS%dampMat(i,j) + pCS%dampMat(j,i)) * 0.5_dp
!            pCS%dampMat(j,i)  =  pCS%dampMat(i,j)
!            pCS%massMat(i,j)  = (pCS%massMat(i,j) + pCS%massMat(j,i)) * 0.5_dp
!            pCS%massMat(j,i)  =  pCS%massMat(i,j)
!       end do
!    end do

    ierr = nStep ! Pass the return status from EstimateControllerProperties

    return

915 call reportError (debugFileOnly_p,'BuildStructControlJacobi')

  contains

    !> @brief Performs the matrix-matrix multiplication Cmat = A*ctrlProp*B.
    !> @details A = pCS%Grad_ForceWrtCout and B = pCS%Grad_CinWrtSensor.
    !> Using LAPACK::DGEMM('N','N',m,n,k,alpha,A,ldA,B,ldB,beta,C,ldC)
    subroutine dmmm (Cmat,ctrlProp,ierr)
      use scratchArrayModule, only : getRealScratchArray
      real(dp), intent(out) :: Cmat(:,:)
      real(dp), intent(in)  :: ctrlProp(:,:)
      integer , intent(out) :: ierr
      real(dp), pointer     :: rWork(:)
      rWork => getRealScratchArray(pCS%nDofs*numVregIn,ierr)
      if (ierr == 0) then
         !! rWork(nDofs,nVregIn) = Grad_ForceWrtCout(nDofs,nVregOut)
         !!                      * ctrlProp(nVregOut,nVregIn)
         call DGEMM ('N','N', pCS%nDofs, numVregIn, numVregOut, &
              &       1.0_dp, pCS%Grad_ForceWrtCout(1,1), pCS%nDofs, &
              &               ctrlProp(1,1), numVregOut, &
              &       0.0_dp, rWork(1), pCS%nDofs)
         !! Cmat(nDofs,nDofs) = -rWork(nDofs,nVregIn)
         !!                   * Grad_CinWrtSensor(nVregIn,nDofs)
         call DGEMM ('N','N', pCS%nDofs,  pCS%nDofs, numVregIn,&
              &      -1.0_dp, rWork(1), pCS%nDofs, &
              &               pCS%Grad_CinWrtSensor(1,1), numVregIn, &
              &       0.0_dp, Cmat(1,1), pCS%nDofs)
      end if
    end subroutine dmmm

    !> @brief Checks if the matrix @b A has non-symmetric terms.
    function isNonSymmetric (A)
      real(dp), intent(in) :: A(:,:)
      integer :: i, j, isNonSymmetric
      isNonSymmetric = 0
      do i = 1, size(A,1)-1
         do j = i+1, size(A,2)
            if (abs(A(i,j)-A(j,i)) > 1.0e-15_dp) then
               isNonSymmetric = 1
               return
            end if
         end do
      end do
    end function isNonSymmetric

  end subroutine BuildStructControlJacobi


  !!============================================================================
  !> @brief Estimates the controller gradient properties by perturbation.
  !>
  !> @param sys System level model data
  !> @param ctrl Control system data
  !> @param[in] msim Matrix of simulation parameters
  !> @param[in] whichVregIn Which vreg inputs to perturb
  !> @param[in] whichVregout Which vreg outputs to read variation from
  !> @param[out] ctrlProps 3D table for storing controller properties
  !> (no. of outputs, no. of controller properties, no. of inputs)
  !> - ctrlProps(:,1,:) = Q
  !> - ctrlProps(:,2,:) = K
  !> - ctrlProps(:,3,:) = C
  !> - ctrlProps(:,4,:) = M
  !> @param[out] ierr Error flag
  !>
  !> @details
  !> We use a perturbation method, simular to the Matrix Stiffness Method /
  !> (Virtual) Displacement Method / Unit Load Method to find the equivalent
  !> mechanical properties of the controller. These controller properties will
  !> be added to existing matrices when conducting modal/eigenvalue analysis.
  !>
  !> In this routine, the controller is limited to be of type PID.
  !> The equation for the system is:
  !>
  !>     M*x'' + C*x' + K*x + Q*int(x)dt = F
  !>
  !> or
  !>
  !>     M*(d2x/dt2) + C*(dx/dt) + K*x + Q*int(x)dt = F
  !>
  !> where M is the mass, C is the damping, K is the stiffness,
  !> and Q is the steady state error elimination.
  !>
  !> Controller input = y, controller output = u.
  !>
  !> The values of interest are:
  !> - du/dy: The change du in output from controller with respect to
  !>          the change dy in input to controller.
  !>          du/dy = proportional gain, Kp.
  !> - du/(int(dy)dt): The change du in output from controller with respect to
  !>                    the change int(dy)dt in input to controller.
  !>                    du/(int(dy)dt) = integral gain, Ki.
  !> - du/(dy/dt): The change du in output from controller with respect to
  !>               the change dy/dt in input to controller.
  !>               du/(dy/dt) = derivative gain, Kd.
  !>
  !> Working order:
  !> -# Do an initial perturbation on the controller with dy = 0 and dt &ne; 0.
  !>    This is to insure dy/dt = 0.
  !> -# Get the initial values y0 and u0 for the controller.
  !> -# Establish dy(j) and dt(j). j = number of perturbations.
  !>    For a PID-controller: j = 1...3.
  !> -# Calculate d(int(dy(j))dt) and d(dy/dt).
  !> -# Calculate y(j) and t(j).
  !> -# Iterate the controller with these new values for the input y(j)
  !>    and time t(j) and save the reaction from the controller u(j)
  !>    due to the change in the input.
  !> -# Calculate du(j) based on u0 and u(j).
  !> -# Calculate Kp, Ki and Kd.
  !> -# Based on sensor type (position, velocity or acceleration),
  !>    calculate Q, K, C and M.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Magne Bratland
  !>
  !> @date 25 Nov 2010

  subroutine EstimateControllerProperties01 (sys, ctrl, msim, &
       &                                     whichVregIn, whichVregOut, &
       &                                     ctrlProps, ierr)

    use SystemTypeModule , only : SystemType, dp
    use ControlTypeModule, only : ControlType, copyCtrl
    use SensorTypeModule , only : POS_p, VEL_p, ACC_p
    use DenseMatrixModule, only : solveAxB
    use ReportErrorModule, only : allocationError, reportError
    use ReportErrorModule, only : debugFileOnly_p

    type(SystemType) , intent(inout) :: sys
    type(ControlType), intent(in)    :: ctrl
    integer          , intent(in)    :: msim(:), whichVregIn(:), whichVregOut(:)
    real(dp)         , intent(out)   :: ctrlProps(:,:,:)
    integer          , intent(out)   :: ierr

    !! Local variables

    integer, parameter :: nPerturb = 3 !< Number of controller perturbations

    type(ControlType), pointer :: ctrlCopy => null()
    type(ControlType), pointer :: ctrlCopyCopy => null()

    real(dp) :: orgTime, orgTimeStep !< original time and timestep values
    real(dp) :: dt0, dt              !< quazi-initial and incremental time steps
    real(dp) :: dy                   !< incremental controller input step
    real(dp) :: dintydt              !< integral of y with respect to time
    real(dp) :: y0                   !< initial values of the controller inputs
    real(dp), allocatable :: uy0(:)  !< controller output when input is y0
    real(dp), allocatable :: uy(:)   !< controller output when input is y=y0+dy
    real(dp) :: dyMatrix(nPerturb,nPerturb) !< matrix of perturbation parameters
    real(dp), allocatable :: duTable(:,:)   !< table for storing du = u(y)-u(y0)

    integer :: i, j, iInput, numVregIn, numVregOut, stat

    !! --- Logic section ---

    ctrlProps = 0.0_dp

    numVregIn  = size(whichVregIn)
    numVregOut = size(whichVregOut)

    allocate(uy0(numVregOut),uy(numVregOut), &
         &   duTable(nPerturb,numVregOut), STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('EstimateControllerProperties01')
       return
    end if

    !! Store the initial values of the time
    orgTime = sys%time
    orgTimeStep = sys%timeStep

    dt0 = orgTimeStep*0.1_dp   !! TODO Magne: Change this value?
    sys%time = orgTime + dt0

    stat = 0
    do i = 1, numVregIn
       iInput = whichVregIn(i)

       call copyCtrl (ctrl,ctrlCopy,ierr)
       if (ierr < 0) goto 915

       dy = 0.0_dp ! Initial perturbation with dy = 0
       call PerturbController (sys,ctrlCopy,msim,iInput,dt0,dy, &
            &                  numVregOut,whichVregOut,uy0,ierr)
       if (ierr < 0) goto 915

       !! Save the value of y0
       y0 = abs(ctrlCopy%vreg(iInput))

       !! Perturbation
       do j = 1, nPerturb
          !! Establish dy-matrix
          dt = orgTimeStep*0.1_dp*j
          dy = dt
          dintydt = (y0 + 0.5_dp*dy)*dt

          dyMatrix(j,1) = dintydt
          dyMatrix(j,2) = dy
          dyMatrix(j,3) = dy/dt

          !! Reset current controller (ctrlCopyCopy) to original state (ctrlCopy)
          call copyCtrl (ctrlCopy,ctrlCopyCopy,ierr)
          if (ierr < 0) goto 915

          !! Perturb
          call PerturbController (sys,ctrlCopyCopy,msim,iInput,dt,dy, &
               &                  numVregOut,whichVregOut,uy,ierr)
          if (ierr < 0) goto 915

          !! Calculate du=uy-uy0 and store in a table
          duTable(j,:) = uy - uy0
       end do

       !! Calculate controller properties
       call solveAxB (dyMatrix,duTable,ierr)
       if (ierr < 0) goto 915

       select case (ctrlCopy%input(iInput)%engine%args(1)%p%entity)
       case (POS_p) ! find dintydt(j), dy(j), ddydt(j)
          ctrlProps(:,1,i) = duTable(1,:)
          ctrlProps(:,2,i) = duTable(2,:)
          ctrlProps(:,3,i) = duTable(3,:)
          stat = IOR(stat,1+2+4)
       case (VEL_p) ! find dintydt(j), dy(j), ddydt(j)
          ctrlProps(:,2,i) = duTable(1,:)
          ctrlProps(:,3,i) = duTable(2,:)
          ctrlProps(:,4,i) = duTable(3,:)
          stat = IOR(stat,2+4+8)
       case (ACC_p) ! find dintydt(j), dy(j)
          ctrlProps(:,3,i) = duTable(1,:)
          ctrlProps(:,4,i) = duTable(2,:)
          stat = IOR(stat,4+8)
       case default
          !! Error
       end select
    end do
    ierr = stat

900 continue

    !! Final reset time
    sys%time = orgTime
    sys%timeStep = orgTimeStep

    call deallocateCtrlCopy (ctrlCopy)
    call deallocateCtrlCopy (ctrlCopyCopy)
    deallocate(uy0,uy,duTable)
    return

915 call reportError (debugFileOnly_p,'EstimateControllerProperties01')
    goto 900

  end subroutine EstimateControllerProperties01


  subroutine EstimateControllerProperties02 (sys, ctrl, msim, &
       &                                     whichVregIn, whichVregOut, &
       &                                     ctrlProps, ierr)

    !!==========================================================================
    !! Purpose:
    !! Use a perturbation method, simular to the Matrix Stiffness Method /
    !! (Virtual) Displacement Method / Unit Load Method to find
    !! the equivalent mechanical properties of the controller. These controller properties
    !! will be added to existing matrices when conducting modal analysis / eigenvalue analysis.
    !! There are a total of 4 values of interest:
    !! mass M, damping C, stiffness K and steady state error elimination Q.
    !! M*x'' + C*x' + K*x + Q*int(x)dt = F or M*(d2x/dt2) + C*(dx/dt) + K*X + Q*int(x)dt = F
    !! Controller input = y, controller output = u.
    !! The method is depending on sensor input (position x, velocity x' or acceleration x'')
    !! - If the sensor input is position:
    !!   du = Q*d(int(y)dt) + K*dy + C*d(dy/dt) + M*d(d2y/dt2)
    !! - If the sensor input is velocity:
    !!   du = Q*d(int(int(y)dt)dt) + K*d(int(y)dt)+ C*dy + M*d(dy/dt)
    !! - If the sensor input is acceleration:
    !!   du = Q*d(int(int(int(y)dt)dt)dt) + K*d(int(int(y)dt)dt) + C*d(int(y)dt)+ M*dy
    !!
    !! Working order:
    !! Per controller input, estimate all controller outputs:
    !! 1) Do one initial perturbation on the controller with dy = 0 and dt /= 0.
    !!    This is to insure dy0/dt = 0.
    !! 2) Get the initial values y0 and u0 for the controller.
    !! 3) Establish dy(j) and dt(j). j = number of perturbations. j = 1...4.
    !! 4) Depending on sensor input:
    !!    Position: Calculate d(int(y(j))dt), d(dy(j)/dt(j)) and d(d2y(j)/dt(j)2).
    !!    Velocity: Calculate d(int(int(y(j))dt)dt), d(int(y(j))dt) and d(dy(j)/dt(j))
    !!    Acceleration: Calculate d(int(int(int(y(j))dt)dt)dt), d(int(int(y(j))dt)dt)
    !!    and d(int(y(j))dt)
    !! 5) Calculate y(j) and t(j).
    !! 6) Iterate the controller with these new values for the input y(j) and time t(j)
    !!    and save the reaction from the controller u(j) due to the change in the input.
    !! 7) Calculate du(j) based on u0 and u(j).
    !! 8) Calculate Q, K, C and M.
    !!
    !! Programmer : Magne Bratland
    !! date/rev   : 01 June 2010 / 1.0
    !!==========================================================================

    use SystemTypeModule , only : SystemType, dp
    use ControlTypeModule, only : ControlType, copyCtrl
    use SensorTypeModule , only : POS_p, VEL_p, ACC_p
    use DenseMatrixModule, only : solveAxB
    use ReportErrorModule, only : allocationError, reportError
    use ReportErrorModule, only : debugFileOnly_p

    type(SystemType)   , intent(inout) :: sys
    type(ControlType)  , intent(in)    :: ctrl
    integer            , intent(in)    :: msim(:)
    integer,             intent(in)    :: whichVregIn(:)    !! Which vreg in to perturb
    integer,             intent(in)    :: whichVregOut(:)   !! Which vreg out to read variation from
    real(dp),            intent(out)   :: ctrlProps(:,:,:)  !! table for storing controller properties
    !						  		                        !!(no. of outputs from controller,
    !				 				                        !! no. of controller properties,
    !		   			   			                        !! no. of inputs to controller)
    !                                                       !! ctrlProps(:,1,:) = Q
    !                                                       !! ctrlProps(:,2,:) = K
    !                                                       !! ctrlProps(:,3,:) = C
    !                                                       !! ctrlProps(:,4,:) = M

    integer            , intent(inout) :: ierr

    ! Local variables

    integer, parameter :: nPerturb = 4 !< Number of controller perturbations

    type(ControlType), pointer :: ctrlCopy => null()
    type(ControlType), pointer :: ctrlCopyCopy => null()

    real(dp) :: orgTime                         !! original/initial value for the time
    real(dp) :: orgTimeStep                     !! original/initial value for the time step
    real(dp) :: dt0                             !! quazi-initial time step
    real(dp) :: dy0                             !! quazi-initial input perturbation, dy0 = 0
    real(dp) :: dt              !! incremental time step
    real(dp) :: dy              !! incremental controller input step
    real(dp) :: dy1             !! incremental controller input step no. 1
    !                                           !! to be used when deriving d(d2y/dt2)
    real(dp) :: dy2             !! incremental controller input step no. 2
    !                                           !! to be used when deriving d(d2y/dt2)
    real(dp) :: y1              !! y1 = y0+dy1
    real(dp) :: ddydt           !! d(dy/dt), 1st derivative of y with respect to time
    real(dp) :: dd2ydt2         !! d(d2y/dt2), 2nd derivative of y with respect to time
    real(dp) :: dintydt         !! definite integral of y with respect to time
    real(dp) :: dintintydt      !! definite double integral of y with respect to time
    real(dp) :: dintintintydt   !! definite triple integral of y with respect to time
    real(dp) :: y0              !! initial values of the controller inputs
    real(dp), allocatable :: uy0(:)             !! u(y0), output from controller when input is y0
    real(dp), allocatable :: uy1(:)             !! u(y1)
    real(dp), allocatable :: uy(:)              !! u(y), output from controller when input is y = y0+dy
    real(dp) :: dyMatrix(nPerturb,nPerturb)     !! matrix of perturbation parameters
    real(dp), allocatable :: duTable(:,:)       !! table for storing du = u(y)-u(y0)

    integer :: i, j, numVregIn, numVregOut

    !! --- Logic section ---

    ctrlProps = 0.0_dp

    numVregIn  = size(whichVregIn)
    numVregOut = size(whichVregOut)

    allocate(uy0(numVregOut),uy1(numVregOut),uy(numVregOut), &
         &   duTable(nPerturb,numVregOut), STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('EstimateControllerProperties02')
       return
    end if

    !! Store the initial values of the time
    orgTime = sys%time
    orgTimeStep = sys%timeStep

    !! Make copy of the controller
    call copyCtrl (ctrl,ctrlCopy,ierr)
    if (ierr < 0) goto 915

    !! Do one time perturbation
    dt0 = sys%timeStep*0.1_dp   !! TODO Magne: Change this value?
    dy0 = 0.0_dp
    sys%time = sys%time + dt0

    !! Start perturbation with time step dt0
    do i = 1, numVregIn
       call PerturbController (sys,ctrlCopy,msim,whichVregIn(i),dt0,dy0, &
            &                  numVregOut,whichVregOut,uy,ierr)
       if (ierr < 0) goto 915
    end do

    !! Save the value of u(y0) in an array
    do i = 1, numVregOut
       uy0(i) = ctrlCopy%vreg(whichVregOut(i))
    end do

    !! Depending on sensor input, there are 3 various ways to perturb the system
    do i = 1, numVregIn
       !! Save the value of y0
       y0 = ctrlCopy%vreg(whichVregIn(i))
       select case (ctrlCopy%input(i)%engine%args(1)%p%entity)
       case (POS_p)
          do j = 1, nPerturb
             !! Establish dy-matrix
             dt  = orgTimeStep*0.1_dp*j
             dy1 = dt
             dy2 = -dy1
             y1  = y0 + dy1
             dintydt = (y1  + 0.5*dp*dy2)*dt
             ddydt   = (dy2 -        dy1)/dt
             dd2ydt2 = (dy2 - 2.0_dp*dy1)/(dt*dt)

             dyMatrix(j,1) = dintydt
             dyMatrix(j,2) = dy2
             dyMatrix(j,3) = ddydt
             dyMatrix(j,4) = dd2ydt2

             !! Reset current controller (ctrlCopyCopy) to original state (ctrlCopy)
             call copyCtrl (ctrlCopy,ctrlCopyCopy,ierr)
             if (ierr < 0) goto 915

             !! To derive d(d2y/dt2), the system has to be perturbed three times (two + initial)

             !! First perturbation
             call PerturbController (sys,ctrlCopyCopy,msim,whichVregIn(i),dt,dy1, &
                  &                  numVregOut,whichVregOut,uy,ierr)
             if (ierr < 0) goto 915

             !! Save the value of u(y1) in an array
             uy1 = uy

             !! Second perturbation
             call PerturbController (sys,ctrlCopyCopy,msim,whichVregIn(i),dt,dy2, &
                  &                  numVregOut,whichVregOut,uy,ierr)
             if (ierr < 0) goto 915

             !! Calculate du=uy-uy1 and store in a table
             duTable(j,:) = uy - uy1
          end do

       case (VEL_p)
          do j = 1, nPerturb
             !! Establish dy-matrix
             dt = orgTimeStep*0.1_dp*j
             dy = dt
             dintydt    = (y0        + dy/2.0_dp)*dt
             dintintydt = (y0/2.0_dp + dy/6.0_dp)*dt*dt

             dyMatrix(j,1) = dintintydt
             dyMatrix(j,2) = dintydt
             dyMatrix(j,3) = dy
             dyMatrix(j,4) = dy/dt

             !! Reset current controller (ctrlCopyCopy) to original state (ctrlCopy)
             call copyCtrl (ctrlCopy,ctrlCopyCopy,ierr)
             if (ierr < 0) goto 915

             !! Perturb
             call PerturbController (sys,ctrlCopyCopy,msim,whichVregIn(i),dt,dy, &
                  &                  numVregOut,whichVregOut,uy,ierr)
             if (ierr < 0) goto 915

             !! Calculate du=uy-uy0 and store in a table
             duTable(j,:) = uy - uy0
          end do

       case (ACC_p)
          do j = 1, nPerturb
             !! Establish dy-matrix
             dt = orgTimeStep*0.1_dp*j
             dy = dt
             dintydt       = (y0        + dy/2.0_dp )*dt
             dintintydt    = (y0/2.0_dp + dy/6.0_dp )*dt*dt
             dintintintydt = (y0/6.0_dp + dy/24.0_dp)*dt*dt*dt

             dyMatrix(j,1) = dintintintydt
             dyMatrix(j,2) = dintintydt
             dyMatrix(j,3) = dintydt
             dyMatrix(j,4) = dy

             !! Reset current controller (ctrlCopyCopy) to original state (ctrlCopy)
             call copyCtrl (ctrlCopy,ctrlCopyCopy,ierr)
             if (ierr < 0) goto 915

             !! Perturb
             call PerturbController (sys,ctrlCopyCopy,msim,whichVregIn(i),dt,dy, &
                  &                  numVregOut,whichVregOut,uy,ierr)
             if (ierr < 0) goto 915

             !! Calculate du=uy-uy0 and store in a table
             duTable(j,:) = uy - uy0
          end do

       case default
          !! Error
          goto 915
       end select

       !! Calculate controller properties
       call solveAxB (dyMatrix,duTable,ierr)
       if (ierr < 0) goto 915

       ctrlProps(:,:,i) = transpose(duTable)
    end do

900 continue

    !! Final reset time
    sys%time = orgTime
    sys%timeStep = orgTimeStep

    call deallocateCtrlCopy (ctrlCopy)
    call deallocateCtrlCopy (ctrlCopyCopy)
    deallocate(uy0,uy1,uy,duTable)
    return

915 call reportError (debugFileOnly_p,'EstimateControllerProperties02')
    goto 900

  end subroutine EstimateControllerProperties02


  subroutine EstimateControllerProperties03 (sys, ctrl, msim, &
       &                                     whichVregIn, whichVregOut, &
       &                                     ctrlProps, ierr)

    !!==========================================================================
    !! Purpose:
    !! Use a perturbation method, simular to the Matrix Stiffness Method /
    !! (Virtual) Displacement Method / Unit Load Method to find
    !! the equivalent mechanical properties of the controller. These controller properties
    !! will be added to existing matrices when conducting modal analysis / eigenvalue analysis.
    !! There are a total of 4 values of interest:
    !! mass M, damping C, stiffness K and steady state error elimination Q.
    !! M*x'' + C*x' + K*x + Q*int(x)dt = F or M*(d2x/dt2) + C*(dx/dt) + K*X + Q*int(x)dt = F
    !! Controller input = y, controller output = u.
    !! The method is depending on sensor input (position x, velocity x' or acceleration x'')
    !! - If the sensor input is position:
    !!   du = Q*d(int(y)dt) + K*dy + C*d(dy/dt) + M*d(d2y/dt2)
    !! - If the sensor input is velocity:
    !!   du = Q*d(int(int(y)dt)dt) + K*d(int(y)dt)+ C*dy + M*d(dy/dt)
    !! - If the sensor input is acceleration:
    !!   du = Q*d(int(int(int(y)dt)dt)dt) + K*d(int(int(y)dt)dt) + C*d(int(y)dt)+ M*dy
    !! This is a just-in-case-algorithm; all possible parameters of interest will be derived:
    !! - d(int(int(int(y(j))dt)dt)dt)
    !! - d(int(int(y(j))dt)dt)
    !! - d(int(y(j))dt
    !! - dy(j)
    !! - d(dy(j)/dt(j))
    !! - d(d2y(j)/dt(j)2)
    !!
    !! Working order:
    !! Per controller input, estimate all controller outputs:
    !! 1) Do one initial perturbation on the controller with dy = 0 and dt /= 0.
    !!    This is to insure dy0/dt = 0.
    !! 2) Get the initial values y0 and u0 for the controller.
    !! 3) Establish dy(j) and dt(j). j = number of perturbations. j = 1...6.
    !! 4) Calculate d(int(int(int(y(j))dt)dt)dt), d(int(int(y(j))dt)dt), d(int(y(j))dt)
    !!    d(dy(j)/dt(j)) and d(d2y(j)/dt(j)2).
    !! 5) Calculate y(j) and t(j).
    !! 6) Iterate the controller with these new values for the input y(j) and time t(j)
    !!    and save the reaction from the controller u(j) due to the change in the input.
    !! 7) Calculate du(j) based on u0 and u(j).
    !! 8) Calculate Q, K, C and M.
    !!
    !! Programmer : Magne Bratland
    !! date/rev   : 01 June 2010 / 1.0
    !!==========================================================================

    use SystemTypeModule , only : SystemType, dp
    use ControlTypeModule, only : ControlType, copyCtrl
    use SensorTypeModule , only : POS_p, VEL_p, ACC_p
    use DenseMatrixModule, only : solveAxB
    use ReportErrorModule, only : allocationError, reportError
    use ReportErrorModule, only : debugFileOnly_p

    type(SystemType)   , intent(inout) :: sys
    type(ControlType)  , intent(in)    :: ctrl
    integer            , intent(in)    :: msim(:)
    integer,             intent(in)    :: whichVregIn(:)    !! Which vreg in to perturb
    integer,             intent(in)    :: whichVregOut(:)   !! Which vreg out to read variation from
    real(dp),            intent(out)   :: ctrlProps(:,:,:)  !! table for storing controller properties
    !						  		                        !!(no. of outputs from controller,
    !				 				                        !! no. of controller properties,
    !		   			   			                        !! no. of inputs to controller)
    !                                                       !! ctrlProps(:,1,:) = Q
    !                                                       !! ctrlProps(:,2,:) = K
    !                                                       !! ctrlProps(:,3,:) = C
    !                                                       !! ctrlProps(:,4,:) = M

    integer            , intent(inout) :: ierr

    ! Local variables

    integer, parameter :: nPerturb = 6 !< Number of controller perturbations

    type(ControlType), pointer :: ctrlCopy => null()
    type(ControlType), pointer :: ctrlCopyCopy => null()

    real(dp) :: orgTime                          !! original/initial value for the time
    real(dp) :: orgTimeStep                      !! original/initial value for the time step
    real(dp) :: dt0                              !! quazi-initial time step
    real(dp) :: dy0                              !! quazi-initial input perturbation, dy0 = 0
    real(dp) :: dt               !! incremental time step
    real(dp) :: dy1              !! incremental controller input step no. 1
    !                                            !! to be used when deriving d(d2y/dt2)
    real(dp) :: dy2              !! incremental controller input step no. 2
    !                                            !! to be used when deriving d(d2y/dt2)
    real(dp) :: y1               !! y1 = y0+dy1
    real(dp) :: ddydt            !! d(dy/dt), 1st derivative of y
    !										     !! with respect to time t
    real(dp) :: dd2ydt2          !! d(d2y/dt2), 2nd derivative of y
    !										     !! with respect to time t
    real(dp) :: dintydt          !! definite integral of y with respect to time t
    real(dp) :: dintintydt       !! definite double integral of y with respect to time t
    real(dp) :: dintintintydt    !! definite triple integral of y with respect to time t
    real(dp) :: y0               !! initial values of the controller inputs
    real(dp), allocatable :: uy0(:)              !! u(y0), output from controller when input is y0
    real(dp), allocatable :: uy1(:)              !! u(y1)
    real(dp), allocatable :: uy(:)               !! u(y), output from controller when input is y = y0+dy
    real(dp) :: dyMatrix(nPerturb,nPerturb)      !! matrix of perturbation parameters
    real(dp), allocatable :: duTable(:,:)        !! table for storing du-results

    integer :: i, j, numVregIn, numVregOut

    !! --- Logic section ---

    ctrlProps = 0.0_dp

    numVregIn  = size(whichVregIn)
    numVregOut = size(whichVregOut)

    allocate(uy0(numVregOut),uy1(numVregOut),uy(numVregOut), &
         &   duTable(nPerturb,numVregOut), STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('EstimateControllerProperties03')
       return
    end if

    !! Store the initial values of the time
    orgTime = sys%time
    orgTimeStep = sys%timeStep

    !! Make copy of controller
    call copyCtrl (ctrl,ctrlCopy,ierr)
    if (ierr < 0) goto 915

    !! Do one time perturbation
    dt0 = sys%timeStep*0.1_dp   !! TODO Magne: Change this value?
    dy0 = 0.0_dp
    sys%time = sys%time + dt0

    !! Start perturbation with time step dt0
    do i = 1, numVregIn
       call PerturbController (sys,ctrlCopy,msim,whichVregIn(i),dt0,dy0, &
            &                  numVregOut,whichVregOut,uy0,ierr)
       if (ierr < 0) goto 915
    end do

    !! Save the value of u(y0) in an array
    do i = 1, numVregOut
       uy0(i) = ctrlCopy%vreg(whichVregOut(i))
    end do

    !! New way: do to many perturbations and exlude redundant results. Do then
    !! 6 three-step perturbations

    do i = 1, numVregIn
       !! Save the value of y0
       y0 = ctrlCopy%vreg(whichVregIn(i))
       do j = 1, nPerturb
          !! Establish dy-matrix
          dt  = orgTimeStep*0.1_dp*j
          dy1 = dt
          dy2 = -dy1
          y1  = y0 + dy1
          dintydt       = (y1        + dy2/2.0_dp )*dt
          dintintydt    = (y1/2.0_dp + dy2/6.0_dp )*dt*dt
          dintintintydt = (y1/6.0_dp + dy2/24.0_dp)*dt*dt*dt
          ddydt   = (dy2 -        dy1)/dt
          dd2ydt2 = (dy2 - 2.0_dp*dy1)/(dt*dt)

          dyMatrix(j,1) = dintintintydt
          dyMatrix(j,2) = dintintydt
          dyMatrix(j,3) = dintydt
          dyMatrix(j,4) = dy2
          dyMatrix(j,5) = ddydt
          dyMatrix(j,6) = dd2ydt2

          !! Reset current controller (ctrlCopyCopy) to original state (ctrlCopy)
          call copyCtrl (ctrlCopy,ctrlCopyCopy,ierr)
          if (ierr < 0) goto 915

          !! To derive d(d2y/dt2), the system has to be perturbed three times (two + initial)

          !! First perturbation
          call PerturbController (sys,ctrlCopyCopy,msim,whichVregIn(i),dt,dy1, &
               &                  numVregOut,whichVregOut,uy1,ierr)
          if (ierr < 0) goto 915

          !! Second perturbation
          call PerturbController (sys,ctrlCopyCopy,msim,whichVregIn(i),dt,dy2, &
               &                  numVregOut,whichVregOut,uy,ierr)
          if (ierr < 0) goto 915

          !! Calculate du=uy-uy1 and store in a table
          duTable(j,:) = uy - uy1
       end do

       !! Calculate controller properties
       call solveAxB (dyMatrix,duTable,ierr)
       if (ierr < 0) goto 915

       select case (ctrlCopy%input(i)%engine%args(1)%p%entity)
       case (POS_p) ! dintydt(j), dy(j), ddydt(j) and dd2ydt2(j)
          ctrlProps(:,1,i) = duTable(3,:)
          ctrlProps(:,2,i) = duTable(4,:)
          ctrlProps(:,3,i) = duTable(5,:)
          ctrlProps(:,4,i) = duTable(6,:)
       case (VEL_p) ! dintintydt(j), dintydt(j), dy(j), ddydt(j)
          ctrlProps(:,1,i) = duTable(2,:)
          ctrlProps(:,2,i) = duTable(3,:)
          ctrlProps(:,3,i) = duTable(4,:)
          ctrlProps(:,4,i) = duTable(5,:)
       case (ACC_p) ! dintintintydt(j), dintintydt(j), dintydt(j), dy(j)
          ctrlProps(:,1,i) = duTable(1,:)
          ctrlProps(:,2,i) = duTable(2,:)
          ctrlProps(:,3,i) = duTable(3,:)
          ctrlProps(:,4,i) = duTable(4,:)
       case default
          !! Error
       end select

    end do

900 continue

    !! Final reset time
    sys%time = orgTime
    sys%timeStep = orgTimeStep

    call deallocateCtrlCopy (ctrlCopy)
    call deallocateCtrlCopy (ctrlCopyCopy)
    deallocate(uy0,uy1,uy,duTable)
    return

915 call reportError (debugFileOnly_p,'EstimateControllerProperties03')
    goto 900

  end subroutine EstimateControllerProperties03


  subroutine EstimateControllerProperties04 (sys, ctrl, msim, &
       &                                     whichVregIn, whichVregOut, &
       &                                     nStep, ctrlProps, ierr)

    !!==========================================================================
    !! Purpose:
    !! Use a perturbation method, simular to the Matrix Stiffness Method /
    !! (Virtual) Displacement Method / Unit Load Method to find
    !! the equivalent mechanical properties of the controller. These controller properties
    !! will be added to existing matrices when conducting modal analysis / eigenvalue analysis.
    !! In this routine, the controller is limited to be of type PID. This algorithm has additional
    !! steps in between perturbations to check for accuracy.
    !! The equation for the system is:
    !! M*x'' + C*x' + K*x + Q*int(x)dt = F or M*(d2x/dt2) + C*(dx/dt) + K*X + Q*int(x)dt = F
    !! where M is mass, C is damping, K is stiffness and Q is steady state error elimination
    !! Controller input = y, controller output = u.
    !! The values of interest are:
    !! du/dy: Change du in output from controller with respect to
    !!        change dy in input to controller.
    !!        du/dy = proportional gain, Kp
    !! du/(int(dy)dt): Change du in output from controller with respect to
    !!                 change int(dy)dt in input to controller.
    !!                 du/(int(dy)dt) = integral gain, Ki
    !! du/(dy/dt): Change du in output from controller with respect to
    !!             change dy/dt in input to controller.
    !!             du/(dy/dt) = derivative gain, Kd
    !!
    !! Working order:
    !! 1) Do one initial perturbation on the controller with dy = 0 and dt /= 0.
    !!    This is to insure dy/dt = 0.
    !! 2) Get the initial values y0 and u0 for the controller.
    !! 3) Establish dy(j) and dt(j). j = number of perturbations.
    !!    For a PID-controller: j = 1...3.
    !! 4) Calculate d(int(dy(j)dt) and d(dy/dt).
    !! 5) Calculate y(j) and t(j).
    !! 6) Iterate the controller with these new values for the input y(j) and time t(j)
    !!    and save the reaction from the controller u(j) due to the change in the input.
    !! 7) Calculate du(j) based on u0 and u(j).
    !! 8) Calculate Kp, Ki and Kd.
    !! 9) Based on sensor type (position, velocity or acceleration), calculate Q, K, C and M.
    !!
    !! Programmer : Magne Bratland
    !! date/rev   : 03 June 2010 / 1.0
    !!==========================================================================

    use SystemTypeModule , only : SystemType, dp
    use ControlTypeModule, only : ControlType, copyCtrl
    use SensorTypeModule , only : POS_p, VEL_p, ACC_p
    use DenseMatrixModule, only : solveAxB
    use ReportErrorModule, only : allocationError, reportError
    use ReportErrorModule, only : debugFileOnly_p

    type(SystemType)   , intent(inout) :: sys
    type(ControlType)  , intent(in)    :: ctrl
    integer            , intent(in)    :: msim(:)
    integer,             intent(in)    :: whichVregIn(:)    !! Which vreg in to perturb
    integer,             intent(in)    :: whichVregOut(:)   !! Which vreg out to read variation from
    integer,             intent(in)    :: nStep             !! number of steps in between perturbations
    real(dp),            intent(out)   :: ctrlProps(:,:,:)  !! table for storing controller properties
    !						  		                        !!(no. of outputs from controller,
    !				 				                        !! no. of controller properties,
    !		   			   			                        !! no. of inputs to controller)
    !                                                       !! ctrlProps(:,1,:) = Q
    !                                                       !! ctrlProps(:,2,:) = K
    !                                                       !! ctrlProps(:,3,:) = C
    !                                                       !! ctrlProps(:,4,:) = M

    integer            , intent(inout) :: ierr

    ! Local variables

    integer, parameter :: nPerturb = 5 !< Number of controller perturbations

    type(ControlType), pointer :: ctrlCopy => null()
    type(ControlType), pointer :: ctrlCopyCopy => null()

    real(dp) :: orgTime                          !! original/initial value for the time
    real(dp) :: orgTimeStep                      !! original/initial value for the time step
    real(dp) :: dt0                              !! quazi-initial time step
    real(dp) :: dy0                              !! quazi-initial input perturbation, dy0 = 0
    real(dp) :: ddt                              !! incremental steps of dt
    real(dp) :: ddy                              !! incremental steps of dy
    real(dp) :: dt               !! incremental time step
    real(dp) :: dy               !! incremental controller input step
    real(dp) :: dintydt          !! definite integral of y with respect to time
    real(dp) :: dintintydt       !! definite double integral of y with respect to time
    real(dp) :: dintintintydt    !! definite triple integral of y with respect to time
    real(dp) :: y0               !! initial values of the controller inputs
    real(dp), allocatable :: uy0(:)              !! u(y0), output from controller when input is y0
    real(dp), allocatable :: uy(:)               !! u(y), output from controller when input is y = y0+dy
    real(dp) :: dyMatrix(nPerturb,nPerturb)      !! matrix of perturbation parameters
    real(dp), allocatable :: duTable(:,:)        !! table for storing du = u(y)-u(y0)

    integer :: i, ii, j, numVregIn, numVregOut

    !! --- Logic section ---

    ctrlProps = 0.0_dp

    numVregIn  = size(whichVregIn)
    numVregOut = size(whichVregOut)

    allocate(uy0(numVregOut),uy(numVregOut), &
         &   duTable(nPerturb,numVregOut), STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('EstimateControllerProperties04')
       return
    end if

    !! Store the initial values of the time
    orgTime = sys%time
    orgTimeStep = sys%timeStep

    !! Make copy of controller
    call copyCtrl (ctrl,ctrlCopy,ierr)
    if (ierr < 0) goto 915

    !! Do one initial time perturbation
    dt0 = sys%timeStep*0.1_dp   !! TODO Magne: Change this value?
    dy0 = 0.0_dp
    sys%time = sys%time + dt0

    !! Start perturbation with time step dt0
    do i = 1, numVregIn
       !! Perturb
       ddt = dt0/real(nStep,dp)
       do ii = 1, nStep
          call PerturbController (sys,ctrlCopy,msim,whichVregIn(i),ddt,dy0, &
               &                  numVregOut,whichVregOut,uy,ierr)
          if (ierr < 0) goto 915
       end do
    end do

    !! Save the value of u(y0) in an array
    do i = 1, numVregOut
       uy0(i) = ctrlCopy%vreg(whichVregOut(i))
    end do

    !! New way: do to many perturbations and exlude redundant results. Do then
    !! 6 three-step perturbations

    do i = 1, numVregIn
       !! Save the value of y0
       y0 = ctrlCopy%vreg(whichVregIn(i))
       do j = 1, nPerturb
          !! Establish dy-matrix
          dt            = orgTimeStep*0.1_dp*j
          dy            = dt
          dintydt       = (y0        + dy/2.0_dp )*dt
          dintintydt    = (y0/2.0_dp + dy/6.0_dp )*dt*dt
          dintintintydt = (y0/6.0_dp + dy/24.0_dp)*dt*dt*dt

          dyMatrix(j,1) = dintintintydt
          dyMatrix(j,2) = dintintydt
          dyMatrix(j,3) = dintydt
          dyMatrix(j,4) = dy
          dyMatrix(j,5) = dy/dt

          !! Reset current controller (ctrlCopyCopy) to original state (ctrlCopy)
          call copyCtrl (ctrlCopy,ctrlCopyCopy)
          if (ierr < 0) goto 915

          !! First perturbation
          ddt = dt/real(nStep,dp)
          ddy = dy/real(nStep,dp)
          do ii = 1, nStep
             call PerturbController (sys,ctrlCopyCopy,msim,whichVregIn(i),ddt,ddy, &
                  &                  numVregOut,whichVregOut,uy,ierr)
             if (ierr < 0) goto 915
          end do

          !! Calculate du=uy-uy0 and store in a table
          duTable(j,:) = uy - uy0
       end do

       !! Calculate controller properties
       call solveAxB (dyMatrix,duTable,ierr)
       if (ierr < 0) goto 915

       select case (ctrlCopy%input(i)%engine%args(1)%p%entity)
       case (POS_p) ! find dintydt(j), dy(j), ddydt(j)
          ctrlProps(:,1,i) = duTable(3,:)
          ctrlProps(:,2,i) = duTable(4,:)
          ctrlProps(:,3,i) = duTable(5,:)
       case (VEL_p) ! find dintintydt(j), dintydt(j), dy(j), ddydt(j)
          ctrlProps(:,1,i) = duTable(2,:)
          ctrlProps(:,2,i) = duTable(3,:)
          ctrlProps(:,3,i) = duTable(4,:)
          ctrlProps(:,4,i) = duTable(5,:)
       case (ACC_p) ! find dintintintydt(j),dintintydt(j),dintydt(j), dy(j)
          ctrlProps(:,1,i) = duTable(1,:)
          ctrlProps(:,2,i) = duTable(2,:)
          ctrlProps(:,3,i) = duTable(3,:)
          ctrlProps(:,4,i) = duTable(4,:)
       case default
          !! Error
       end select

    end do

900 continue

    !! Final reset time
    sys%time = orgTime
    sys%timeStep = orgTimeStep

    call deallocateCtrlCopy (ctrlCopy)
    call deallocateCtrlCopy (ctrlCopyCopy)
    deallocate(uy0,uy,duTable)
    return

915 call reportError (debugFileOnly_p,'EstimateControllerProperties04')
    goto 900

  end subroutine EstimateControllerProperties04


  !!============================================================================
  !> @brief Perturbation method without initial perturbation.
  !> @details Perturbation w.r.t. a single integral.
  !> @a nPerturb = 1 in this subroutine, therefore no loop necessary.
  !>
  !> @author Magne Bratland
  !>
  !> @date 7 July 2010

  subroutine EstimateControllerProperties500 (sys, ctrl, msim, &
       &                                      whichVregIn, whichVregOut, &
       &                                      ctrlProps, ierr)

    use SystemTypeModule , only : SystemType, dp
    use ControlTypeModule, only : ControlType, copyCtrl
    use SensorTypeModule , only : POS_p, VEL_p, ACC_p
    use ReportErrorModule, only : allocationError, reportError
    use ReportErrorModule, only : debugFileOnly_p

    type(SystemType) , intent(inout) :: sys
    type(ControlType), intent(in)    :: ctrl
    integer          , intent(in)    :: msim(:)
    integer          , intent(in)    :: whichVregIn(:)   !! Which vreg in to perturb
    integer          , intent(in)    :: whichVregOut(:)  !! Which vreg out to read variation from
    real(dp)         , intent(out)   :: ctrlProps(:,:,:) !! table for storing controller properties
    integer          , intent(inout) :: ierr

    !! Local variables

    type(ControlType), pointer :: ctrlCopy => null()

    real(dp) :: y0          ! initial values of the controller inputs
    real(dp), allocatable :: uy0(:) ! u(y0), output from controller when input to controller is y0
    real(dp), allocatable :: uy(:)  ! u(y), output from controller when input to controller is y = y0+dy
    real(dp) :: orgTime     ! original/initial value for the time
    real(dp) :: orgTimeStep ! original/initial value for the time step
    real(dp) :: dt          ! incremental time step
    real(dp) :: dy          ! incremental controller input step
    real(dp) :: ddydt       ! d(dy/dt), 1st derivative of y with respect to time t
    real(dp) :: dintydt     ! definite integral of y with respect to time t

    integer :: i, iEntity, numVregIn, numVregOut

    !! --- Logic section ---

    ctrlProps = 0.0_dp

    numVregIn  = size(whichVregIn)
    numVregOut = size(whichVregOut)

    allocate(uy0(numVregOut),uy(numVregOut), STAT=ierr)
    if (ierr /= 0) then
       ierr = allocationError('EstimateControllerProperties500')
       return
    end if

    !! Store the initial time values
    orgTime = sys%time
    orgTimeStep = sys%timeStep

    !! Do one initial time perturbation
    sys%time = sys%time + sys%timeStep

    !! Save the value of u(y0) in an array
    do i = 1, numVregOut
       uy0(i) = ctrl%vreg(whichVregOut(i))
    end do

    !! Make copy of controller
    call copyCtrl (ctrl,ctrlCopy,ierr)
    if (ierr < 0) goto 915

    !! Perturbation
    do i = 1, numVregIn
       y0 = ctrlCopy%vreg(whichVregIn(i))

       !! Establish dy-matrix
       dt      = orgTimeStep
       dy      = dt
       ddydt   = dy/dt
       dintydt = (y0 + dy/2.0_dp)*dt

       !! The perturbation sequence
       !! Reset current controller (ctrlCopy) to original state (ctrl)
       call copyCtrl (ctrl,ctrlCopy,ierr)
       if (ierr < 0) goto 915

       !! Perturb the controller
       call PerturbController (sys,ctrlCopy,msim,whichVregIn(i),dt,dy, &
            &                  numVregOut,whichVregOut,uy,ierr)
       if (ierr < 0) goto 915

       iEntity = ctrlCopy%input(i)%engine%args(1)%p%entity
       select case (iEntity)
       case (POS_p,VEL_p,ACC_p)
          !! find dintydt(j), dy(j), ddydt(j)
          ctrlProps(:,iEntity,i) = (1.0_dp/dintydt) * (uy - uy0)
       case default ! Error
       end select

    end do

900 continue

    !! Final reset time
    sys%time = orgTime
    sys%timeStep = orgTimeStep

    call deallocateCtrlCopy (ctrlCopy)
    deallocate(uy0,uy)
    return

915 call reportError (debugFileOnly_p,'EstimateControllerProperties500')
    goto 900

  end subroutine EstimateControllerProperties500


  !!============================================================================
  !> @brief Perturbs one of the inputs of the control system.
  !>
  !> @param sys System level model data
  !> @param ctrl Control system data
  !> @param[in] msim Matrix of simulation parameters
  !> @param[in] iPert Which input to perturb
  !> @param[in] dt Incremental time step
  !> @param[in] dy Incremental step for use in du/dy
  !> @param[in] numVregOut Number of outputs from the controller to read
  !> @param[in] whichVregOut Which outputs from the controller to read
  !> @param[out] uy Perturbed output from controller, u(y) = u(y0+dy)
  !> @param[out] ierr Error flag
  !>
  !> @details This subroutine perturbs one of the control inputs and calculates
  !> the reaction in all of the control outputs.
  !>
  !> Working order:
  !> 1) Change the input y for the controller from y0 to y = y0+dy,
  !>    where dy is a small number. Change the time from t0 to t = t0+dt,
  !>    where dt is a small number.
  !> 2) Run this new input through the controller to get the reaction
  !>    (output) u(y) from the controller due to the change in the input.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Magne Bratland
  !>
  !> @date 18 Sept 2009

  subroutine PerturbController (sys,ctrl,msim,iPert,dt,dy, &
       &                        numVregOut,whichVregOut,uy,ierr)

    use SystemTypeModule     , only : SystemType, dp
    use ControlTypeModule    , only : ControlType
    use ControlRoutinesModule, only : IterateControlSystem

    type(SystemType) , intent(inout) :: sys
    type(ControlType), intent(inout) :: ctrl
    integer          , intent(in)    :: msim(:), iPert
    real(dp)         , intent(in)    :: dt, dy
    integer          , intent(in)    :: numVregOut, whichVregOut(:)
    real(dp)         , intent(out)   :: uy(:)
    integer          , intent(out)   :: ierr

    !! Local variables
    integer, parameter :: ctrlSysMode = 3 !< controller integration mode

    !! --- Logic section ---

    !! Change time step to dt
    sys%timeStep = dt

    !! Change input from y0 to y (y = y0+dy)
    ctrl%vreg(ctrl%input(iPert)%var) = ctrl%vreg(ctrl%input(iPert)%var) + dy

    !! Iterate the controller with the new values for y (y = y0+dy)
    !! and t (t = t0+dt), to derive the value of u(y)
    call IterateControlSystem (sys,ctrl,ctrlSysMode,msim,ierr,setInput=.false.)

    !! Save the value of u(y) in an array
    call DGATHR (numVregOut,whichVregOut(1),ctrl%vreg(1),uy(1),1)

  end subroutine PerturbController


  !!============================================================================
  !> @brief Deallocates a control system copy.
  !>
  !> @param ctrl Pointer to controltypemodule::controltype object to deallocate.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 19 Jan 2024

  subroutine deallocateCtrlCopy (ctrl)

    use ControlTypeModule, only : ControlType, deallocateCtrl

    type(ControlType), pointer :: ctrl

    !! --- Logic section ---

    if (associated(ctrl)) then
       call deallocateCtrl (ctrl,.false.)
       deallocate(ctrl)
    end if

  end subroutine deallocateCtrlCopy

end module ControlStructModule
