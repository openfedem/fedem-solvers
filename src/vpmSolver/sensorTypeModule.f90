!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file sensorTypeModule.f90
!> @brief Sensor data container.

!!==============================================================================
!> @brief Module with data types representing sensor objects
!> (general function arguments).
!>
!> @details The module also contains subroutines for accessing or manipulating
!> the sensor data.

module SensorTypeModule

  use KindModule  , only : dp
  use IdTypeModule, only : IdType

  implicit none

  !> Sensor type strings, i.e., what type of object the sensor is measuring
  character(len=14), parameter:: sensorType_p(12) =(/ 'TIME          ', &
       &                                              'ENGINE        ', &
       &                                              'CONTROL       ', &
       &                                              'TRIAD         ', &
       &                                              'DAMPER_AXIAL  ', &
       &                                              'SPRING_AXIAL  ', &
       &                                              'JOINT_VARIABLE', &
       &                                              'SPRING_JOINT  ', &
       &                                              'DAMPER_JOINT  ', &
       &                                              'RELATIVE_TRIAD', &
       &                                              'STRAIN_GAGE   ', &
       &                                              'NUM_ITERATIONS' /)

  integer, parameter :: TIME_p           =  1 !< Time sensor enum value
  integer, parameter :: ENGINE_p         =  2 !< Engine sensor enum value
  integer, parameter :: CONTROL_p        =  3 !< Control sensor enum value
  integer, parameter :: TRIAD_p          =  4 !< Triad sensor enum value
  integer, parameter :: DAMPER_AXIAL_p   =  5 !< Axial sensor damper enum value
  integer, parameter :: SPRING_AXIAL_p   =  6 !< Axial sensor spring enum value
  integer, parameter :: JOINT_VARIABLE_p =  7 !< Joint sensor enum value
  integer, parameter :: SPRING_JOINT_p   =  8 !< Joint spring sensor enum value
  integer, parameter :: DAMPER_JOINT_p   =  9 !< Joint damper sensor enum value
  integer, parameter :: RELATIVE_TRIAD_p = 10 !< Relative sensor enum value
  integer, parameter :: STRAIN_GAGE_p    = 11 !< Strain gage sensor enum value
  integer, parameter :: NUM_ITERATIONS_p = 12 !< Iterations sensor enum value

  !> Sensor entity strings, i.e., which result quantity the sensor is measuring
  character(len=7), parameter :: sensorEntity_p(14)= (/ 'POS    ', &
       &                                                'VEL    ', &
       &                                                'ACC    ', &
       &                                                'DEFL   ', &
       &                                                'FORCE  ', &
       &                                                'LENGTH ', &
       &                                                'REL_POS', &
       &                                                'ANGLE  ', &
       &                                                'W_SPEED', &
       &                                                'F_VEL  ', &
       &                                                'F_ACC  ', &
       &                                                'DYN_P  ', &
       &                                                'STRAIN ', &
       &                                                'STRESS ' /)

  integer, parameter :: POS_p     = 1 !< Position sensor enum value
  integer, parameter :: VEL_p     = 2 !< Velocity sensor enum value
  integer, parameter :: ACC_p     = 3 !< Acceleration sensor enum value
  integer, parameter :: DEFL_p    = 4 !< Deflection sensor enum value
  integer, parameter :: FORCE_p   = 5 !< Force sensor enum value
  integer, parameter :: LENGTH_p  = 6 !< Length sensor enum value
  integer, parameter :: REL_POS_p = 7 !< Relative position sensor enum value
  integer, parameter :: ANGLE_p   = 8 !< Relative angle sensor enum value
  integer, parameter :: W_SPEED_p = 9 !< Wind speed sensor enum value
  integer, parameter :: F_VEL_p  = 10 !< Fluid velocity sensor enum value
  integer, parameter :: F_ACC_p  = 11 !< Fluid acceleration sensor enum value
  integer, parameter :: DYN_P_p  = 12 !< Dynamic pressure sensor enum value
  integer, parameter :: STRAIN_p = 13 !< Strain sensor enum value
  integer, parameter :: STRESS_p = 14 !< Stress sensor enum value

  !> Coordinate system type strings
  character(len=6), parameter :: sensorSystem_p(2) = (/ 'GLOBAL', &
       &                                                'LOCAL ' /)

  integer, parameter :: GLOBAL_p = 1 !< Global coordinate system enum value
  integer, parameter :: LOCAL_p  = 2 !< Local coordinate system enum value


  !> @brief Data type representing a sensor object.
  type SensorType
     type(IdType) :: id     !< General identification data
     integer      :: type   !< Object type to be measured
     integer      :: entity !< Result quantity type to be measured
     integer      :: system !< Global/local flag for measured quantity
     integer      :: dof    !< Local dof index of measured quantity
     integer , pointer :: index(:) !< Indices to the measured object(s)
     real(dp), pointer :: value !< Pointer to current sensor value
     logical      :: allocValue !< If .true., the @a value is allocated
  end type SensorType

  !> @brief Data type representing a sensor pointer.
  !> @details This data type is used to construct arrays of sensor objects
  !> where each array element is a pointer to a sensor object,
  !> and not the sensors themselves.
  type SensorPtrType
     type(SensorType), pointer :: p !< Pointer to a sensor
  end type SensorPtrType

  !> Pointer to the one and only time sensor of the model
  type(SensorType), pointer, save :: ourTime => null()


  !> @brief Returns pointer to object with specified ID.
  interface GetPtrToId
     module procedure GetPtrToIdSensor
  end interface

  !> @brief Standard routine for writing an object to file.
  interface WriteObject
     module procedure WriteSensorType
  end interface


contains

  !!============================================================================
  !> @brief Returns pointer to (first) sensor object with specified ID.
  !>
  !> @param[in] array Array of sensortypemodule::sensortype objects
  !>                  to search within
  !> @param[in] id Base ID of the object to search for
  !>
  !> @details If the sensor is not found, NULL is returned.
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 1 Nov 1999

  function GetPtrToIdSensor (array,id) result(ptr)
    type(SensorType), pointer :: ptr

    integer         , intent(in)        :: id
    type(SensorType), intent(in),target :: array(:)

    !! Local variables
    integer :: i

    !! --- Logic section ---

    do i = 1, size(array)
       if (array(i)%id%baseId == id) then
          ptr => array(i)
          return
       end if
    end do

    nullify(ptr)
    write(*,*) '*** GetPtrToIdSensor returned nullified, baseId =',id

  end function GetPtrToIdSensor


  !!============================================================================
  !> @brief Standard routine for writing an object to file.
  !>
  !> @param[in] sensor The sensortypemodule::sensortype object to write
  !> @param[in] lpu File unit number to write to
  !> @param[in] complexity If present, the value indicates the amount of print
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 10 Sep 1998

  subroutine WriteSensorType (sensor,lpu,complexity)

    use IdTypeModule, only : writeId

    type(SensorType), intent(in) :: sensor
    integer         , intent(in) :: lpu
    integer,optional, intent(in) :: complexity

    !! Local variables
    integer :: i, idx

    !! --- Logic section ---

    write(lpu,'(A)') 'Sensor','{'
    call writeId (sensor%id,lpu)

    if (sensor%type > 0 .and. sensor%type <= size(sensorType_p)) then
       write(lpu,*) 'type        = ', sensorType_p(sensor%type)
    else
       write(lpu,*) 'type        = UNDEFINED'
    end if
    if (sensor%entity > 0 .and. sensor%entity <= size(sensorEntity_p)) then
       write(lpu,*) 'entity      = ', sensorEntity_p(sensor%entity)
    end if
    if (sensor%system > 0 .and. sensor%system <= size(sensorSystem_p)) then
       write(lpu,*) 'system      = ', sensorSystem_p(sensor%system)
    end if
    if (sensor%dof > 0) then
       write(lpu,*) 'dof         =', sensor%dof
    end if

    if (associated(sensor%index)) then
       if (size(sensor%index) > 0) then
          idx = sensor%index(1)
       else
          idx = 0
       end if
       select case (sensor%type)

       case (TRIAD_p, RELATIVE_TRIAD_p)
          do i = 1, size(sensor%index)
             if (sensor%index(i) > 0) then
                write(lpu,*) ' triad'//char(ichar('0')+i)//' #    =', &
                     &       sensor%index(i)
             end if
          end do

       case (JOINT_VARIABLE_p)
          if (idx  > 0) write(lpu,*) 'joint #     =', idx

       case (SPRING_AXIAL_p, SPRING_JOINT_p)
          if (idx  > 0) write(lpu,*) 'spring #    =', idx

       case (DAMPER_AXIAL_p, DAMPER_JOINT_p)
          if (idx  > 0) write(lpu,*) 'damper #    =', idx

       case (ENGINE_p)
          if (idx  > 0) write(lpu,*) 'engine #    =', idx

       case (STRAIN_GAGE_p)
          if (idx  > 0) write(lpu,*) 'strainGage # =', idx

       end select
    end if

    if (present(complexity)) then
       if (complexity >= 2 .and. associated(sensor%value)) then
          write(lpu,*) 'value       =', sensor%value
       end if
    end if

    write(lpu,'(A)') '}'

  end subroutine WriteSensorType


  !!============================================================================
  !> @brief Returns the id of a sensor object, including its type.
  !>
  !> @param[in] sensor The sensortypemodule::sensortype object to get id for
  !> @param[in] typeAsPrefix If .true., use the sensor type as prefix
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 8 Feb 2024

  function getSensorId (sensor,typeAsPrefix) result(text)

    use IdTypeModule, only : lId_p, getId

    type(SensorType) , intent(in) :: sensor
    logical, optional, intent(in) :: typeAsPrefix

    !! Local variables
    character(len=lId_p) :: text

    !! --- Logic section ---

    text = getId(sensor%id)
    if (sensor%type < 1 .or. sensor%type > size(sensorType_p)) then
       text = 'Sensor'//text(1:lId_p-6)
    else if (.not. present(typeAsPrefix)) then
       text = 'Sensor'//text(1:lId_p-6)
    else if (typeAsPrefix) then
       text = trim(sensorType_p(sensor%type))//text
    else
       text = 'Sensor'//trim(text)//' of type '//sensorType_p(sensor%type)
    end if

  end function getSensorId


  !!============================================================================
  !> @brief Initializes a sensor object.
  !>
  !> @param sensor The sensortypemodule::sensortype object to initialize
  !>
  !> @callgraph @callergraph
  !>
  !> @author Bjorn Haugen
  !>
  !> @date 10 Sep 1998

  subroutine NullifySensor (sensor)

    use IdTypeModule, only : nullifyId

    type(SensorType), intent(out) :: sensor

    !! --- Logic section ---

    call nullifyId (sensor%id)

    sensor%type   = 0
    sensor%entity = 0
    sensor%system = 0
    sensor%dof    = 0
    nullify(sensor%index)
    nullify(sensor%value)
    sensor%allocValue = .false.

  end subroutine NullifySensor


  !!============================================================================
  !> @brief Deallocates a sensor object.
  !>
  !> @param sensor The sensortypemodule::sensortype object to deallocate
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Jan 2017

  subroutine DeallocateSensor (sensor)

    use IdTypeModule, only : deallocateId

    type(SensorType), intent(inout) :: sensor

    !! --- Logic section ---

    call deallocateId (sensor%id)

    if (associated(sensor%index)) deallocate(sensor%index)
    if (sensor%allocValue) deallocate(sensor%value)

    call nullifySensor (sensor)

  end subroutine DeallocateSensor


  !!============================================================================
  !> @brief Deallocates all sensor objects.
  !>
  !> @param sensors Array of all sensortypemodule::sensortype objects
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Jan 2017

  subroutine deallocateSensors (sensors)

    type(SensorType), pointer :: sensors(:)

    !! Local variables
    integer :: i

    !! --- Logic section ---

    do i = 1, size(sensors)
       if (sensors(i)%type == TIME_p .and. associated(ourTime,sensors(i))) then
          nullify(ourTime)
       end if
       call DeallocateSensor (sensors(i))
    end do

    deallocate(sensors)
    nullify(sensors)

    if (associated(ourTime)) then
       deallocate(ourTime)
       nullify(ourTime)
    end if

  end subroutine deallocateSensors


  !!============================================================================
  !> @brief Returns the number of control output sensors in the model.
  !>
  !> @param sensors Array of all sensortypemodule::sensortype objects
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 5 Oct 2006

  function numCtrlOut (sensors) result(nOut)

    type(SensorType), intent(in) :: sensors(:)

    !! Local variables
    integer :: i, nOut

    !! --- Logic section ---

    nOut = 0
    do i = 1, size(sensors)
       if (sensors(i)%type == CONTROL_p) nOut = nOut + 1
    end do

  end function numCtrlOut

end module SensorTypeModule
