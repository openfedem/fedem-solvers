!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file supElTypeModule.f90
!> @brief Superelement object data container.

!!==============================================================================
!> @brief Module with data types representing superelement objects.
!>
!> @details The module also contains subroutines for accessing or manipulating
!> the superelement data.

module SupElTypeModule

  use KindModule     , only : dp, lfnam_p
  use TriadTypeModule, only : TriadPtrType, IdType

  implicit none

  !> @brief Data type for the generalized DOFs associated with component modes.
  type GeneralizedDofs

     integer :: nDOFs     !< Number of generalized DOFs (NGEN)
     integer :: samNodNum !< Node number for SAM reference
     integer :: firstDOF  !< Index to first generalized DOF in superelement

     integer , pointer :: sysDOF !< System DOF number for first generalized DOF

     integer , pointer :: BC(:)     !< Boundary condition codes (0=fixed 1=free)
     real(dp), pointer :: alpha1(:) !< Mass-proportional damping factors
     real(dp), pointer :: alpha2(:) !< Stiffness-proportional damping factors
     real(dp), pointer :: ur(:)     !< Current displacement of generalized DOFs
     real(dp), pointer :: urPrev(:) !< Previous isplacement of generalized DOFs
     real(dp), pointer :: urd(:)    !< Velocity of generalized DOFs
     real(dp), pointer :: urdd(:)   !< Acceleration of generalized DOFs
     real(dp), pointer :: energy(:,:) !< Energy contributions from each DOF

  end type GeneralizedDofs


  !> @brief Data type for the nonlinear force-displacement representation.
  type NonlinForceStiffType

     real(dp), pointer :: force(:,:)   !< Forces for all points
     real(dp), pointer :: disp(:,:)    !< Displacements for all points
     real(dp), pointer :: stiff(:,:,:) !< Stiffness for all points

  end type NonlinForceStiffType


  !> @brief Data type for the hydrodynamic force calculation.
  type HydroDynType

     integer :: bodyIndex !< Handle to volume data for hydrodynamics calculation

     real(dp), pointer :: dragPar(:,:) !< Input parameters for drag calculation

     real(dp) :: slamPar(3) !< Slamming parameters
     real(dp) :: Morison(9) !< Morison equation coefficients for beams
     real(dp) :: wn(3)      !< Current normal vector of water surface
     real(dp) :: Vb(2)      !< Buoyancy volume
     real(dp) :: C0b(3,2)   !< Buoyancy center
     real(dp) :: As(2)      !< Waterline area
     real(dp) :: C0As(3,2)  !< Waterline center
     real(dp) :: C0s(3)     !< Slam attack point
     real(dp) :: B(6)       !< Buoyancy forces in the center of gravity
     real(dp) :: M(3)       !< Added mass forces in the center of gravity
     real(dp) :: D(6)       !< Drag forces in the center of gravity
     real(dp) :: S(6)       !< Slam forces in the center of gravity

  end type HydroDynType


  !> @brief Data type for the integrated stress recovery.
  type RecoveryType
     integer                :: recovery    !< Flag for stress/gage recovery
     character(len=lfnam_p) :: fileName(4) !< File names for stress recovery
     character(len=128)     :: elmGroup    !< Element groups to do recovery for
  end type RecoveryType


  !> @brief Data type representing a superelement object.
  type SupElType

     type(IdType) :: id !< General identification data

     integer :: samElNum  !< Element number for SAM reference
     integer :: nExtNods  !< Number of external nodes (NENOD)
     integer :: nTotDOFs  !< Number of total DOFs (NDIM)
     integer :: nLoadCase !< Number of load cases (NLC)

     !!=== Data for controlling the co-rotated coordinate system's position ===
     integer :: shadowPosAlg !< Flag for co-rotated position update algorithm

     !! shadowPosAlg = 1 : Based on max triangle triads
     integer  :: refTriad(3)    !< Nodal indices of the three reference triads
     real(dp) :: offset(3,3)    !< Reference point offsets from triads
     real(dp) :: triRelSup(3,4) !< Position of triangle relative to superelement

     !! shadowPosalg = 2 : Based on weighted average of displacement being zero
     !> A 6*NEDOF matrix which gives the shadow element displacements as
     !> weighted average of all element DOFs. When multiplied with a "proper"
     !> deformational displacement vector, all averages should be zero.
     real(dp), pointer :: shadowPosGrad(:,:)

     real(dp) :: supTr(3,4)     !< Superelement position, current time step
     real(dp) :: supTrPrev(3,4) !< Superelement position, previous timestep
     real(dp) :: supTrInit(3,4) !< Superelement position, initial configuration

     !! Triad arrays
     type(TriadPtrType), pointer :: triads(:)  !< All triads on the superelement
     integer           , pointer :: nodeId(:)  !< External FE node numbers
     real(dp), pointer :: TrUndeformed(:,:,:)  !< Undeformed position matrices
     real(dp), pointer :: eccVec(:,:)          !< Nodal eccentricity vectors

     type(GeneralizedDofs), pointer :: genDOFs !< Data for the generalized DOFs

     logical  :: addedMass        !< Additional masses on associated triads?
     real(dp) :: mass             !< Superelement mass
     real(dp) :: posMassCenter(3) !< Mass center position in superelement system
     real(dp) :: inertia(6,6)     !< Rigid body inertia properties
     !                            !! about mass center in superelement system

     real(dp) :: ePot0 !< Initial potential energy
     real(dp) :: ePot  !< Current potential energy (relative to initial)
     real(dp) :: eKin  !< Current kinetic energy
     real(dp) :: eStr  !< Current strain energy
     real(dp) :: eDmp  !< Energy loss from damping

     integer  :: rigidFlag          !< 0: Flexible, 1: Generic, 2: Rigid
     integer  :: stressStiffFlag(3) !< Include stress stiffening?

     integer  :: massCorrFlag !< Perform moment correction for spin?
     real(dp) :: RotError     !< Current moment error, if the correction is off

     real(dp) :: mDmpFactor !< Mass proportional damping      (ALPHA1)
     real(dp) :: kDmpFactor !< Stiffness proportional damping (ALPHA2)
     real(dp) :: mDmp0      !< Initial mass proportional damping factor
     real(dp) :: kDmp0      !< Initial stiffness proportional damping factor
     real(dp) :: dmpScl(2)  !< Current and previous scaling of mDmp0 and kDmp0
     integer  :: dmpSclIdx  !< Index to structural damping scaling function

     real(dp) :: massScl(2) !< Current and previous mass scaling factors
     integer  :: massSclIdx !< Index to mass scaling function
     real(dp) :: stifScl(2) !< Current and previous stiffness scaling factors
     integer  :: stifSclIdx !< Index to stiffness scaling function

     type(NonlinForceStiffType), pointer :: nonlin !< Nonlinear link data
     type(HydroDynType)        , pointer :: hydyn  !< Hydrodynamics data
     type(RecoveryType)        , pointer :: rcy    !< Data for stress recovery

     real(dp), pointer :: scoord !< Running coordinate along beams
     real(dp), pointer :: EI(:)  !< Bending stiffness for beams

     real(dp), pointer :: Nmat(:,:)  !< Newton matrix
     real(dp), pointer :: Mmat(:,:)  !< Structural mass matrix
     real(dp), pointer :: Mamat(:,:) !< Virtual added mass matrix
     real(dp), pointer :: Cmat(:,:)  !< Structural damping matrix
     real(dp), pointer :: Cdmat(:,:) !< Drag damping matrix
     real(dp), pointer :: KmMat(:,:) !< Material stiffness matrix
     real(dp), pointer :: KlMat(:,:) !< Load correction stiffness matrix
     real(dp), pointer :: KtMat(:,:) !< Tangent stiffness matrix

     real(dp), pointer :: Bgp(:,:)    !< B-matrix for generic part CoG triad
     real(dp), pointer :: rtr_rt(:,:) !< The matrix ((R^t*R)^-1)*R^t
     !                                !! for mass matrix correction routines

     real(dp), pointer :: Q(:)  !< External Forces
     real(dp), pointer :: FS(:) !< Forces related to the stiffness matrix
     real(dp), pointer :: FD(:) !< Forces related to the damping matrix
     real(dp), pointer :: FI(:) !< Forces related to the inertia matrix

     real(dp), pointer :: uld(:)       !< Velocity in local direction
     real(dp), pointer :: uldd(:)      !< Acceleration in local direction
     real(dp), pointer :: vld(:)       !< Deformational velocity
     real(dp), pointer :: vldPrev(:)   !< Deformational velocity, previous step
     real(dp), pointer :: vldd(:)      !< Deformational acceleration
     real(dp), pointer :: vlddPrev(:)  !< Deformational acceleration, previous
     real(dp), pointer :: finit(:)     !< Deformational displacements
     real(dp), pointer :: finitPrev(:) !< Deformational displacements, previous
     real(dp), pointer :: fg(:,:)      !< Gravitational force vectors
     real(dp), pointer :: S(:,:)       !< External load vectors

     logical :: savePos    !< Flag indicating whether position should be saved
     logical :: saveVar(4) !< Flags indicating which variables should be saved

  end type SupElType


  !> @brief Data type representing a superelement pointer.
  !> @details This data type is used to construct arrays of superelements where
  !> each element is a pointer to an object, and not the objects themselves.
  type SupElPtrType
     type(SupElType), pointer :: p
  end type SupElPtrType


  !> @brief Returns pointer to object with specified ID.
  interface GetPtrToId
     module procedure GetPtrToIdSupEl
  end interface

  !> @brief Standard routine for writing an object to file.
  interface WriteObject
     module procedure WriteSupElType
  end interface

  !> @brief Updates the state variables pertaining to previous time step.
  interface updateAtConvergence
     module procedure updatePreviousValues
  end interface

  !> @brief Restores the state variables from the last converged time step.
  interface restoreFromLastStep
     module procedure restorePreviousValues
  end interface

  private :: GetPtrToIdSupEl, WriteSupElType
  private :: updatePreviousValues, restorePreviousValues


contains

  !!============================================================================
  !> @brief Returns pointer to (first) superelement with specified ID.
  !>
  !> @param[in] array Array of supeltypemodule::supeltype objects
  !>                  to search within
  !> @param[in] id Base ID of the object to search for
  !> @param[in] userId If @e .true., search for a user ID instead
  !>
  !> @details If the superelement is not found, NULL is returned.
  !>
  !> @author Bjorn Haugen / Knut Morten Okstad
  !>
  !> @date 13 Jan 2010

  function GetPtrToIdSupEl (array,id,userId) result(ptr)

    integer        , intent(in)           :: id
    type(SupElType), intent(in), target   :: array(:)
    logical        , intent(in), optional :: userId
    type(SupElType), pointer              :: ptr

    !! Local variables
    integer :: i
    logical :: searchUserId

    !! --- Logic section ---

    if (present(userId)) then
       searchUserId = userId
    else
       searchUserId = .false.
    end if

    do i = 1, size(array)
       if (searchUserId) then
          if (array(i)%id%userId /= id) cycle
       else
          if (array(i)%id%baseId /= id) cycle
       end if
       ptr => array(i)
       return
    end do

    nullify(ptr)
    if (searchUserId) then
       write(*,*) '*** GetPtrToIdSupEl returned nullified, userId =',id
    else
       write(*,*) '*** GetPtrToIdSupEl returned nullified, baseId =',id
    end if

  end function GetPtrToIdSupEl


  !!============================================================================
  !> @brief Standard routine for writing an object to file.
  !>
  !> @param[in] sup The supeltypemodule::supeltype object to write
  !> @param[in] io File unit number to write to
  !> @param[in] complexity If present, the value indicates the amount of print
  !>
  !> @author Karl Erik Thoresen
  !>
  !> @date 27 Sep 1998

  subroutine WriteSupElType (sup,io,complexity)

    use IdTypeModule, only : writeId

    type(SupElType) , intent(in) :: sup
    integer         , intent(in) :: io
    integer,optional, intent(in) :: complexity

    !! Local variables
    integer :: i, j

    !! --- Logic section ---

    write(io,'(A)') 'Superelement','{'
    call writeId (sup%id,io)

    write(io,*) 'samElNum    =', sup%samElNum
    write(io,*) 'nExtNod     =', sup%nExtNods
    write(io,*) 'nTotDofs    =', sup%nTotDofs
    write(io,*) 'nLoadCase   =', sup%nLoadCase

    if (associated(sup%hydyn)) then
       write(io,*) 'bodyIndex   =', sup%hydyn%bodyIndex
       if (associated(sup%hydyn%dragPar)) then
          write(io,1) 'dragPar     =', sup%hydyn%dragPar
       endif
       write(io,3) 'slamPar     =', sup%hydyn%slamPar
       write(io,3) 'Morison     =', sup%hydyn%Morison
    end if

    if (associated(sup%genDOFs)) then
       write(io,*) 'nGenDOFs    =', sup%genDOFs%nDOFs
       write(io,*) 'samNodNum   =', sup%genDOFs%samNodNum
       write(io,*) 'firstDOF    =', sup%genDOFs%firstDOF
       if (associated(sup%genDOFs%sysDOF)) then
          write(io,*) 'sysGenDOF   =', sup%genDOFs%sysDOF
       end if
    end if

    write(io,*) 'shadowPosAlg  =', sup%shadowPosAlg
    if (sup%shadowPosAlg == 1) then
       write(io,*) 'refTriad(idx) =', sup%refTriad
       write(io,1) 'offset        =', (sup%offset(i,:),i=1,3)
       write(io,2) 'triRelSup     =', (sup%triRelSup(i,:),i=1,3)
    end if
    write(io,2) 'supTr         =', (sup%supTr(i,:),i=1,3)
    if (present(complexity)) then
       if (complexity >= 2) then
          write(io,2) 'supTrPrev     =', (sup%supTrPrev(i,:),i=1,3)
          write(io,2) 'supTrInit     =', (sup%supTrInit(i,:),i=1,3)
       end if
    end if
    if (associated(sup%scoord)) then
       write(io,3) 'scoord        =', sup%scoord
    end if
    if (associated(sup%EI)) then
       write(io,3) 'EIy,EIz       =', sup%EI
    end if

    write(io,*) 'addedMass?    =', sup%addedMass
    write(io,3) 'mass          =', sup%mass
    write(io,3) 'posMassCenter =', sup%posMassCenter
    write(io,6) 'inertia       =', (sup%inertia(i,:),i=1,6)
    write(io,*) 'rigid?        =', sup%rigidFlag
    write(io,*) 'stressStiff?  =', sup%stressStiffFlag
    write(io,*) 'massCorr?     =', sup%massCorrFlag
    write(io,3) 'mDmpFactor    =', sup%mDmpFactor
    write(io,3) 'kDmpFactor    =', sup%kDmpFactor
    write(io,*) 'dmpScale(id)  =', sup%dmpSclIdx
    write(io,*) 'massScale(id) =', sup%massSclIdx
    write(io,*) 'stifScale(id) =', sup%stifSclIdx
    write(io,*) 'savePos       =', sup%savePos
    write(io,*) 'saveVar       =', sup%saveVar

    if (present(complexity)) then
       if (complexity >= 1) then

          if (associated(sup%triads)) then
             write(io,9) 'triads(id)', (sup%triads(j)%p%id%baseId, &
                  &                     j=1,size(sup%triads))
             write(io,9) 'dofStart  ', (sup%triads(j)%firstDOF, &
                  &                     j=1,size(sup%triads))
          end if

          if (associated(sup%nodeId)) then
             write(io,9) 'nodes(id) ', sup%nodeId
          end if

          if (associated(sup%TRundeformed)) then
             do i = 1, size(sup%TRundeformed,3)
                write(io,"(' TRundeformed, supNode',i3,' ='/(4ES15.6E3))") i, &
                     &  (sup%TRundeformed(j,:,i),j=1,3)
             end do
          end if

          if (associated(sup%eccVec)) then
             do i = 1, size(sup%eccVec,2)
                write(io,"(' eccVec(',i3,')  =',3ES15.6E3)") i,sup%eccVec(:,i)
             end do
          end if

          if (associated(sup%shadowPosGrad)) then
             write(io,*) 'shadowPosGrad ='
             do i = 1, size(sup%shadowPosGrad,2)
                write(io,4) sup%shadowPosGrad(:,i)
             end do
          end if

          if (associated(sup%fg)) then
             do j = 1, size(sup%fg,2)
                write(io,*) 'fg'//char(ichar('W')+j)//'           ='
                write(io,4) sup%fg(:,j)
             end do
          end if
          if (associated(sup%S)) then
             do j = 1, size(sup%S,2)
                write(io,"(' S',i3,'          =')") j
                write(io,4) sup%S(:,j)
             end do
          end if

          if (associated(sup%genDOFs)) then
             if (associated(sup%genDOFs%BC)) then
                write(io,"(A,30I2/(16X,30I2))")' BCgen        =',sup%genDOFs%BC
             end if
             if (associated(sup%genDOFs%alpha1)) then
                write(io,*) 'alpha1        ='
                write(io,4) sup%genDOFs%alpha1
             end if
             if (associated(sup%genDOFs%alpha2)) then
                write(io,*) 'alpha2        ='
                write(io,4) sup%genDOFs%alpha2
             end if
          end if

       end if
       if (complexity >= 2) then

          if (associated(sup%hydyn)) then
             write(io,3) 'wn            =', sup%hydyn%wn
             write(io,3) 'As            =', sup%hydyn%As
             write(io,3) 'Vb            =', sup%hydyn%Vb
             write(io,3) 'C0s           =', sup%hydyn%C0s
             write(io,3) 'C0b           =', sup%hydyn%C0b
             write(io,3) 'B             =', sup%hydyn%B
             write(io,3) 'M             =', sup%hydyn%M
             write(io,3) 'D             =', sup%hydyn%D
             write(io,3) 'S             =', sup%hydyn%S
          end if
          if (associated(sup%genDOFs)) then
             if (associated(sup%genDOFs%ur)) then
                write(io,*) 'urGen         ='
                write(io,4) sup%genDOFs%ur
             end if
             if (associated(sup%genDOFs%urPrev)) then
                write(io,*) 'urGenPrev     ='
                write(io,4) sup%genDOFs%urPrev
             end if
             if (associated(sup%genDOFs%urd)) then
                write(io,*) 'urdGen        ='
                write(io,4) sup%genDOFs%urd
             end if
             if (associated(sup%genDOFs%urdd)) then
                write(io,*) 'urddGen       ='
                write(io,4) sup%genDOFs%urdd
             end if
          end if

          if (associated(sup%uld)) then
             write(io,*) 'uld           ='
             write(io,4) sup%uld
          end if
          if (associated(sup%uldd)) then
             write(io,*) 'uldd          ='
             write(io,4) sup%uldd
          end if

          if (associated(sup%vld)) then
             write(io,*) 'vld           ='
             write(io,4) sup%vld
          end if
          if (associated(sup%vldd)) then
             write(io,*) 'vldd          ='
             write(io,4) sup%vldd
          end if

          if (associated(sup%finit)) then
             write(io,*) 'finit         ='
             write(io,4) sup%finit
          end if

          write(io,3) 'RotError      =', sup%RotError
          write(io,3) 'ePot0         =', sup%ePot0
          write(io,3) 'ePot          =', sup%ePot
          write(io,3) 'eKin          =', sup%eKin
          write(io,3) 'eStr          =', sup%eStr
          write(io,3) 'eDmp          =', sup%eDmp

       end if
       if (complexity >= 3 .or. complexity < 0) then

          write(io,*)
          if (associated(sup%shadowPosGrad)) then
             write(io,*) 'size(shadowPosGrad) =', size(sup%shadowPosGrad,1) &
                  &                             , size(sup%shadowPosGrad,2)
          end if
          write(io,*) 'size(triads)        =', size(sup%triads)
          write(io,*) 'size(TrUndeformed)  =', size(sup%TrUndeformed,1) &
               &                             , size(sup%TrUndeformed,2) &
               &                             , size(sup%TrUndeformed,3)
          if (associated(sup%eccVec)) then
             write(io,*) 'size(eccVec)  =', size(sup%eccVec,1) &
                  &                       , size(sup%eccVec,2)
          end if
          if (associated(sup%Nmat)) then
             write(io,*) 'size(Nmat)          =', size(sup%Nmat,1) &
                  &                             , size(sup%Nmat,2)
          end if
          if (associated(sup%Mmat)) then
             write(io,*) 'size(Mmat)          =', size(sup%Mmat,1) &
                  &                             , size(sup%Mmat,2)
          end if
          if (associated(sup%Mamat)) then
             write(io,*) 'size(Mamat)         =', size(sup%Mamat,1) &
                  &                             , size(sup%Mamat,2)
          end if
          if (associated(sup%Cmat)) then
             write(io,*) 'size(Cmat)          =', size(sup%Cmat,1) &
                  &                             , size(sup%Cmat,2)
          end if
          if (associated(sup%Cdmat)) then
             write(io,*) 'size(Cdmat)         =', size(sup%Cdmat,1) &
                  &                             , size(sup%Cdmat,2)
          end if
          if (associated(sup%KmMat)) then
             write(io,*) 'size(KmMat)         =', size(sup%KmMat,1) &
                  &                             , size(sup%KmMat,2)
          end if
          if (associated(sup%KlMat)) then
             write(io,*) 'size(KlMat)         =', size(sup%KlMat,1) &
                  &                             , size(sup%KlMat,2)
          end if
          if (associated(sup%KtMat)) then
             write(io,*) 'size(KtMat)         =', size(sup%KtMat,1) &
                  &                             , size(sup%KtMat,2)
          end if
          if (associated(sup%Bgp)) then
             write(io,*) 'size(Bgp)           =', size(sup%Bgp,1) &
                  &                             , size(sup%Bgp,2)
          end if
          if (associated(sup%rtr_rt)) then
             write(io,*) 'size(rtr_rt)        =', size(sup%rtr_rt,1) &
                  &                             , size(sup%rtr_rt,2)
          end if
          if (associated(sup%Q)) then
             write(io,*) 'size(Q)             =', size(sup%Q)
          end if
          if (associated(sup%FS)) then
             write(io,*) 'size(FS)            =', size(sup%FS)
          end if
          if (associated(sup%FD)) then
             write(io,*) 'size(FD)            =', size(sup%FD)
          end if
          if (associated(sup%FI)) then
             write(io,*) 'size(FI)            =', size(sup%FI)
          end if
          if (associated(sup%genDOFs)) then
             if (associated(sup%genDOFs%BC)) then
                write(io,*) 'size(BC)            =', size(sup%genDOFs%BC)
             end if
             if (associated(sup%genDOFs%alpha1)) then
                write(io,*) 'size(alpha1)        =', size(sup%genDOFs%alpha1)
             end if
             if (associated(sup%genDOFs%alpha2)) then
                write(io,*) 'size(alpha2)        =', size(sup%genDOFs%alpha2)
             end if
             if (associated(sup%genDOFs%urd)) then
                write(io,*) 'size(urdGen)        =', size(sup%genDOFs%urd)
             end if
             if (associated(sup%genDOFs%urdd)) then
                write(io,*) 'size(urddGen)       =', size(sup%genDOFs%urdd)
             end if
          end if
          if (associated(sup%uld)) then
             write(io,*) 'size(uld)           =', size(sup%uld)
          end if
          if (associated(sup%uldd)) then
             write(io,*) 'size(uldd)          =', size(sup%uldd)
          end if
          if (associated(sup%vld)) then
             write(io,*) 'size(vld)           =', size(sup%vld)
          end if
          if (associated(sup%vldd)) then
             write(io,*) 'size(vldd)          =', size(sup%vldd)
          end if
          if (associated(sup%finit)) then
             write(io,*) 'size(finit)         =', size(sup%finit)
          end if
          if (associated(sup%fg)) then
             write(io,*) 'size(fg)            =', size(sup%fg,1), size(sup%fg,2)
          end if
          if (associated(sup%S)) then
             write(io,*) 'size(S)             =', size(sup%S,1), size(sup%S,2)
          end if

       end if
    end if

    write(io,'(A)') '}'

1   format(1X,A,3(/3ES15.6E3))
2   format(1X,A,3(/4ES15.6E3))
3   format(1X,A,  10ES15.6E3 )
4   format(        6ES15.6E3 )
6   format(1X,A,6(/6ES15.6E3))
9   format(1X,A,'    =',10I9/(16X,10I9))

  end subroutine WriteSupElType


  !!============================================================================
  !> @brief Initializes a superelement object.
  !>
  !> @param sup The supeltypemodule::supeltype object to initialize
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 13 Oct 2000

  subroutine NullifySupEl (sup)

    use IdTypeModule, only : nullifyId

    type(SupElType), intent(out) :: sup

    !! --- Logic section ---

    call nullifyId (sup%id)

    sup%samElNum = 0
    sup%nExtNods = 0
    sup%nTotDOFs = 0
    sup%nLoadCase = 0

    sup%shadowPosAlg = 0
    sup%refTriad = 0
    nullify(sup%shadowPosGrad)
    sup%offset = 0.0_dp
    sup%triRelSup = 0.0_dp
    sup%supTr = 0.0_dp
    sup%supTrPrev = 0.0_dp
    sup%supTrInit = 0.0_dp

    nullify(sup%triads)
    nullify(sup%nodeId)
    nullify(sup%TrUndeformed)
    nullify(sup%eccVec)
    nullify(sup%genDOFs)

    sup%addedMass = .false.
    sup%mass = 0.0_dp
    sup%posMassCenter = 0.0_dp
    sup%inertia = 0.0_dp
    sup%ePot0 = 0.0_dp
    sup%ePot = 0.0_dp
    sup%eKin = 0.0_dp
    sup%eStr = 0.0_dp
    sup%eStr = 0.0_dp
    sup%eDmp = 0.0_dp
    sup%rigidFlag = 0
    sup%stressStiffFlag = 0
    sup%massCorrFlag = 0
    sup%RotError = 0.0_dp
    sup%mDmpFactor = 0.0_dp
    sup%kDmpFactor = 0.0_dp
    sup%mDmp0 = 0.0_dp
    sup%kDmp0 = 0.0_dp
    sup%dmpScl = 1.0_dp
    sup%massScl = 1.0_dp
    sup%stifScl = 1.0_dp
    sup%dmpSclIdx = 0
    sup%massSclIdx = 0
    sup%stifSclIdx = 0
    sup%savePos = .false.
    sup%saveVar = .false.

    nullify(sup%Nmat)
    nullify(sup%Mmat)
    nullify(sup%Mamat)
    nullify(sup%Cmat)
    nullify(sup%Cdmat)
    nullify(sup%KmMat)
    nullify(sup%KlMat)
    nullify(sup%KtMat)
    nullify(sup%Bgp)
    nullify(sup%rtr_rt)
    nullify(sup%Q)
    nullify(sup%FS)
    nullify(sup%FD)
    nullify(sup%FI)
    nullify(sup%uld)
    nullify(sup%uldd)
    nullify(sup%vld)
    nullify(sup%vldPrev)
    nullify(sup%vldd)
    nullify(sup%vlddPrev)
    nullify(sup%finit)
    nullify(sup%finitPrev)
    nullify(sup%fg)
    nullify(sup%S)
    nullify(sup%nonlin)
    nullify(sup%hydyn)
    nullify(sup%rcy)
    nullify(sup%scoord)
    nullify(sup%EI)

  end subroutine NullifySupEl


  !!============================================================================
  !> @brief Initiates the hydrodynamics quantities for a superelement. 
  !>
  !> @param hydyn Hydrodynamics data container
  !> @param[in] MorPar Morison input parameters
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 25 Mar 2009

  subroutine InitiateHyDyn (hydyn,MorPar)

    use KindModule, only : pi_p

    type(HydroDynType), intent(out) :: hydyn
    real(dp), optional, intent(in)  :: MorPar(:)

    !! --- Logic section ---

    hydyn%bodyIndex = -1
    nullify(hydyn%dragPar)
    hydyn%slamPar = 0.0_dp
    hydyn%Morison = 0.0_dp

    hydyn%wn   = 0.0_dp
    hydyn%Vb   = 0.0_dp
    hydyn%C0b  = 0.0_dp
    hydyn%As   = 0.0_dp
    hydyn%C0As = 0.0_dp
    hydyn%C0s  = 0.0_dp
    hydyn%B    = 0.0_dp
    hydyn%M    = 0.0_dp
    hydyn%D    = 0.0_dp
    hydyn%S    = 0.0_dp

    if (present(MorPar)) then
       !! Initiate coefficients of the Morison equation
       hydyn%As(1) = 0.25_dp*pi_p*MorPar(4)*MorPar(4) ! Ad = pi*Dd^2/4
       hydyn%Morison(1) = hydyn%As(1)*MorPar(1)       ! Ca*Ad
       hydyn%Morison(2) = hydyn%As(1)*MorPar(2)       ! Cm*Ad
       hydyn%Morison(3) = 0.5_dp*MorPar(4)*MorPar(3)  ! Cd*Dd/2
       hydyn%Morison(4) = 0.5_dp*pi_p*MorPar(4)*MorPar(6) ! Cd(axial)*pi*Dd
       hydyn%Morison(5) = 0.5_dp*pi_p*MorPar(4)*MorPar(7) ! Cd(spin)*pi*Dd
       hydyn%Morison(6) = hydyn%As(1)*MorPar(8)       ! Ca(axial)*Ad
       hydyn%Morison(7) = hydyn%As(1)*MorPar(9)       ! Cm(axial)*Ad
       hydyn%Morison(8) = hydyn%As(1)*MorPar(10)      ! Cs*Ad
       hydyn%Morison(9) = MorPar(11)                  ! Slam decay time
       hydyn%As(1) = 0.25_dp*pi_p*MorPar(5)*MorPar(5) ! Ab = pi*Db^2/4
    end if

  end subroutine InitiateHyDyn


  !!============================================================================
  !> @brief Transforms a reference triad offset vector.
  !>
  !> @param sup The supeltypemodule::supeltype object to transform for
  !> @param[in] offset The reference point offset vector to transform
  !> @param[in] index Which reference point (1, 2 or 3) to transform for
  !>
  !> @details The offset vector is transformed from the superelement coordinate
  !> system to the coordinate system of the reference triad itself.
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 13 Nov 2015

  subroutine TransformOffset (sup,offset,index)

    type(SupElType), intent(inout) :: sup
    real(dp)       , intent(in)    :: offset(3)
    integer        , intent(in)    :: index

    !!  offset          = ur^t * supTr * offset = ((supTr*offset)^t * ur)^t
    sup%offset(:,index) = matmul(matmul(sup%supTr(:,1:3),offset), &
         &                       sup%triads(sup%refTriad(index))%p%ur(:,1:3))

  end subroutine TransformOffset


  !!============================================================================
  !> @brief Deallocates a superelement object.
  !>
  !> @param sup The supeltypemodule::supeltype object to deallocate
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 23 Jan 2017

  subroutine DeallocateSupEl (sup)

    use IdTypeModule, only : deallocateId

    type(SupElType), intent(inout) :: sup

    !! --- Logic section ---

    call deallocateId (sup%id)

    if (associated(sup%shadowPosGrad)) deallocate(sup%shadowPosGrad)
    if (associated(sup%triads))        deallocate(sup%triads)
    if (associated(sup%nodeId))        deallocate(sup%nodeId)
    if (associated(sup%TrUndeformed))  deallocate(sup%TrUndeformed)
    if (associated(sup%eccVec))        deallocate(sup%eccVec)
    if (associated(sup%genDOFs)) then
       deallocate(sup%genDOFs%BC,sup%genDOFs%alpha1,sup%genDOFs%alpha2)
       deallocate(sup%genDOFs%ur,sup%genDOFs%urPrev,sup%genDOFs%energy)
       deallocate(sup%genDOFs)
    end if
    if (associated(sup%Nmat))      deallocate(sup%Nmat)
    if (associated(sup%Mmat))      deallocate(sup%Mmat)
    if (associated(sup%Mamat))     deallocate(sup%Mamat)
    if (associated(sup%Cmat))      deallocate(sup%Cmat)
    if (associated(sup%Cdmat))     deallocate(sup%Cdmat)
    if (associated(sup%Kmmat))     deallocate(sup%Kmmat)
    if (associated(sup%Klmat))     deallocate(sup%Klmat)
    if (associated(sup%Ktmat))     deallocate(sup%Ktmat)
    if (associated(sup%Bgp))       deallocate(sup%Bgp)
    if (associated(sup%rtr_rt))    deallocate(sup%rtr_rt)
    if (associated(sup%Q))         deallocate(sup%Q)
    if (associated(sup%FS))        deallocate(sup%FS)
    if (associated(sup%FD))        deallocate(sup%FD)
    if (associated(sup%FI))        deallocate(sup%FI)
    if (associated(sup%uld))       deallocate(sup%uld)
    if (associated(sup%uldd))      deallocate(sup%uldd)
    if (associated(sup%vld))       deallocate(sup%vld)
    if (associated(sup%vldPrev))   deallocate(sup%vldPrev)
    if (associated(sup%vldd))      deallocate(sup%vldd)
    if (associated(sup%vlddPrev))  deallocate(sup%vlddPrev)
    if (associated(sup%finit))     deallocate(sup%finit)
    if (associated(sup%finitPrev)) deallocate(sup%finitPrev)
    if (associated(sup%fg))        deallocate(sup%fg)
    if (associated(sup%S))         deallocate(sup%S)
    if (associated(sup%nonlin)) then
       deallocate(sup%nonlin%force,sup%nonlin%disp,sup%nonlin%stiff)
       deallocate(sup%nonlin)
    end if
    if (associated(sup%hydyn)) then
       if (associated(sup%hydyn%dragPar)) deallocate(sup%hydyn%dragPar)
       deallocate(sup%hydyn)
    end if
    if (associated(sup%rcy))    deallocate(sup%rcy)
    if (associated(sup%scoord)) deallocate(sup%scoord)
    if (associated(sup%EI))     deallocate(sup%EI)

    call nullifySupEl (sup)

  end subroutine DeallocateSupEl


  !!============================================================================
  !> @brief Deallocates an array of superelement objects.
  !>
  !> @param sups The supeltypemodule::supeltype objects to deallocate
  !>
  !> @callgraph @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> date 23 Jan 2017

  subroutine deallocateSupEls (sups)

    type(SupElType), pointer :: sups(:)

    !! Local variables
    integer :: i

    !! --- Logic section ---

    do i = 1, size(sups)
       call DeallocateSupEl (sups(i))
    end do

    deallocate(sups)
    nullify(sups)

  end subroutine deallocateSupEls


  !!============================================================================
  !> @brief Updates the co-rotated system for a superelement.
  !>
  !> @param sup The supeltypemodule::supeltype object to update for
  !> @param dbgUnit File unit number for debug print out
  !> @param stat Status flag, negative value on output indicates an error
  !>
  !> @details If @a stat equals 1 on input, the relative position matrix
  !> for the triangle is initiated (used only for @a shadowPosAlg == 1).
  !>
  !> @callgraph @callergraph
  !>
  !> @author Karl Erik Thoresen
  !> @date Sep 1999
  !>
  !> @author Bjorn Haugen
  !> @date Nov 2003

  subroutine updateSupElCorot (sup,dbgUnit,stat)

    use IdTypeModule     , only : getId
    use manipMatrixModule, only : trans3p, matmul34, invert34
    use rotationModule   , only : vec_to_mat
    use reportErrorModule, only : reportError, error_p, getErrorFile

    type(SupElType), intent(inout) :: sup
    integer        , intent(in)    :: dbgUnit
    integer        , intent(inout) :: stat

    !! Local variables
    integer             :: j, error
    integer , parameter :: maxIter_p = 10
    real(dp)            :: pnt(3,3), triTr(3,4)
    real(dp), parameter :: rotTol_p = 1.0e-10_dp, traTol_p = 1.0e-12_dp
    real(dp)            :: urRigDisp(3,4), rotErr, traErr, rCG(3), dispCG(6)

    !! --- Logic section ---

    if (sup%shadowPosAlg == 1) then

       !! Find the three points defining the triangle
       do j = 1, 3
          pnt(:,j) = matmul34(sup%triads(sup%refTriad(j))%p%ur,sup%offset(:,j))
       end do

       !! Calculate the coordinate system for the triangle
       triTr = trans3P(pnt(:,1),pnt(:,2),pnt(:,3),getErrorFile(),error)
       if (error /= 0) then
          stat = -1
          call reportError (error_p, &
               'Cannot calculate co-rotated reference system '// &
               'for superelement'//getId(sup%id), &
               'Check reference triad positions and offset vectors', &
               addString='updateSupElCorot')
       else if (stat == 1) then
          !! Calculate relative position matrix for the triangle
          sup%triRelSup = matmul34(invert34(triTr),sup%supTr)
       else
          !! Update the co-rotated superelement transformation matrix
          sup%supTr = matmul34(triTr,sup%triRelSup)
       end if

    else if (sup%shadowPosAlg == 2) then

       if (dbgUnit > 0) then
          write(dbgUnit,"('====== Superelement',A)") trim(getId(sup%id))
       end if

       do j = 1, maxIter_p

          call BuildFinit (sup)
          dispCG = matmul(sup%shadowPosGrad,sup%finit)
          traErr = sqrt(dot_product(dispCG(1:3),dispCG(1:3)))
          rotErr = sqrt(dot_product(dispCG(4:6),dispCG(4:6)))

          if (dbgUnit > 0) then
             write(dbgUnit,600) j, dispCG, traErr, rotErr
 600         format('iter=',I2,'  dispCG=',1P,6E11.3,'  traErr,rotErr=',2E11.3)
          end if
          if (rotErr <= rotTol_p .and. traErr <= traTol_p) exit

          !! Establish the rotation matrix
          call vec_to_mat (dispCG(4:6),urRigDisp(:,1:3))

          !! Rotated position of center of gravity
          rCG = matmul(urRigDisp(:,1:3),sup%posMassCenter)

          !! Origin of new shadow element coordinate system measured
          !! in current shadow element coordinate system
          urRigDisp(:,4) = dispCG(1:3) + sup%posMassCenter - rCG

          !! Increment the C0n shadow element position
          sup%supTr = matmul34(sup%supTr,urRigDisp)

       end do

    end if

  end subroutine updateSupElCorot


  !!============================================================================
  !> @brief Updates the state variables pertaining to the previous time step.
  !>
  !> @param sup The supeltypemodule::supeltype object to update state for
  !>
  !> @details This subroutine is invoked once after convergence as been reached.
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 28 Oct 2008

  subroutine updatePreviousValues (sup)

    type(SupElType), intent(inout) :: sup

    !! --- Logic section ---

    sup%supTrPrev = sup%supTr
    sup%finitPrev = sup%finit
    if (associated(sup%genDOFs))  sup%genDOFs%urPrev = sup%genDOFs%ur
    if (associated(sup%vldPrev))  sup%vldPrev  = sup%vld
    if (associated(sup%vlddPrev)) sup%vlddPrev = sup%vldd

  end subroutine updatePreviousValues


  !!============================================================================
  !> @brief Restores the state variables from the last converged time step.
  !>
  !> @param sup The supeltypemodule::supeltype object to restore state for
  !>
  !> @details This subroutine is invoked when doing iteration cut-back.
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 28 Oct 2008

  subroutine restorePreviousValues (sup)

    type(SupElType), intent(inout) :: sup

    !! --- Logic section ---

    sup%supTr = sup%supTrPrev
    sup%finit = sup%finitPrev
    if (associated(sup%genDOFs))  sup%genDOFs%ur = sup%genDOFs%urPrev
    if (associated(sup%vldPrev))  sup%vld  = sup%vldPrev
    if (associated(sup%vlddPrev)) sup%vldd = sup%vlddPrev

    if (associated(sup%hydyn)) then
       sup%hydyn%Vb(1) = sup%hydyn%Vb(2)
       sup%hydyn%C0b(:,1) = sup%hydyn%C0b(:,2)
       if (sup%hyDyn%bodyIndex >= 0) then
          sup%hydyn%As(1) = sup%hydyn%As(2)
          sup%hydyn%C0As(:,1) = sup%hydyn%C0As(:,2)
       end if
    end if

  end subroutine restorePreviousValues


  !!============================================================================
  !> @brief Establishes the deformation vector (finit) for a superelement.
  !>
  !> @param sup The supeltypemodule::supeltype object to obtain deformation for
  !>
  !> @callergraph
  !>
  !> @author Karl Erik Thoresen
  !>
  !> @date 3 Sep 1998

  subroutine BuildFinit (sup)

    use manipMatrixModule, only : matmul34, invert34

    type(SupElType), intent(inout) :: sup

    !! Local variables
    integer  :: i, j, n
    real(dp) :: invSupTr(3,4) !< The inverse of the superelement position matrix
    real(dp) :: urLocal(3,4)  !< Position matrix in superelement coordinates
    real(dp) :: dR(3,3)       !< Skew-symmetric rotation tensor

    !! --- Logic section ---

    invSupTr = invert34(sup%supTr)

    do i = 1, size(sup%triads)
       n = sup%triads(i)%p%nDOFs
       if (n < 3) cycle

       urLocal = matmul34(invSuptr,sup%triads(i)%p%ur)

       !! add in translational dofs
       j = sup%triads(i)%firstDOF
       sup%finit(j:j+2) = urLocal(:,4) - sup%Trundeformed(:,4,i)

       if (n >= 6) then

          !! add in rotational dofs
          dR = matmul(urLocal(:,1:3),sup%Trundeformed(:,1:3,i))
          !bh call mat_to_vec(dR,sup%finit(j+3:j+5))
          sup%finit(j+3) = dR(3,2)
          sup%finit(j+4) = dR(1,3)
          sup%finit(j+5) = dR(2,1)

       end if

    end do

    if (.not. associated(sup%genDOFs)) return
    if (sup%genDOFs%nDOFs < 1) return

    !! add in generalized dofs
    i = sup%genDOFs%firstDOF
    j = i + sup%genDOFs%nDOFs-1
    sup%finit(i:j) = sup%genDOFs%ur

  end subroutine BuildFinit


  !!============================================================================
  !> @brief Returns the transformation matrix to system directions for a triad.
  !>
  !> @param[in] sup The supeltypemodule::supeltype object to transform for
  !> @param[in] n 1-based index of the superelement triad to transform for
  !>
  !> @details The returned 3x3 matrix is to be used to transform vectors from
  !> the superelement directions to system directions for the @a n'th connected
  !> triad of the superelement.
  !>
  !> @callergraph
  !>
  !> @author Karl Erik Thoresen
  !>
  !> @date Dec 1998

  function GetRsupToSys (sup,n) result(RsupToSys)

    type(SupElType), intent(in) :: sup
    integer        , intent(in) :: n
    real(dp)                    :: RsupToSys(3,3)

    !! --- Logic section ---

    if (associated(sup%triads(n)%p%sysDirInG)) then
       RsupToSys = matmul(transpose(sup%triads(n)%p%sysDirInG),sup%supTr(:,1:3))
    else
       RsupToSys = sup%supTr(:,1:3)
    end if

  end function GetRsupToSys


  !!============================================================================
  !> @brief Gets current velocity and acceleration from all superelements.
  !>
  !> @param[in] sups All supeltypemodule::supeltype objects in the model
  !> @param[out] velGlobal System velocity vector
  !> @param[out] accGlobal System acceleration vector
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 1 Nov 2005

  subroutine GetSupElsVelAcc (sups,velGlobal,accGlobal)

    type(SupElType), intent(in)  :: sups(:)
    real(dp)       , intent(out) :: velGlobal(:), accGlobal(:)

    !! Local variables
    integer :: i, j, k

    !! --- Logic section ---

    do i = 1, size(sups)
       if (.not. associated(sups(i)%genDOFs)) cycle

       j = sups(i)%genDOFs%sysDOF
       k = j + sups(i)%genDOFs%nDOFs-1

       velGlobal(j:k) = sups(i)%genDOFs%urd
       accGlobal(j:k) = sups(i)%genDOFs%urdd

    end do

  end subroutine GetSupElsVelAcc


  !!============================================================================
  !> @brief Returns whether a superelement is a beam element or not.
  !>
  !> @param[in] sup The supeltypemodule::supeltype object to check
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 2 Jul 2013

  function IsBeam (sup)

    type(SupElType), intent(in) :: sup
    logical                     :: IsBeam

    !! --- Logic section ---

    IsBeam = sup%nExtNods == 2 .and. sup%rigidFlag < 0

  end function IsBeam


  !!============================================================================
  !> @brief Returns the full ID of a superelement.
  !>
  !> @param[in] sup The supeltypemodule::supeltype object to get the ID for
  !>
  !> @details The full ID of an object consists for of the type name of the
  !> object, the user ID and the user description.
  !>
  !> @callergraph
  !>
  !> @author Knut Morten Okstad
  !>
  !> @date 15 Sep 2017

  function GetSupElId (sup)

    use IdTypeModule, only : lId_p, getId

    type(SupElType), intent(in) :: sup
    character(len=5+lId_p)      :: GetSupElId

    !! --- Logic section ---

    if (IsBeam(sup)) then
       GetSupElId = 'Beam'//getId(sup%id)
    else
       GetSupElId = 'Part'//getId(sup%id)
    end if

  end function GetSupElId

end module SupElTypeModule


!!==============================================================================
!> @brief Module with a namelists for reading superelement data.

module SupElNamelistModule

  use KindModule  , only : dp, lfnam_p
  use IdTypeModule, only : ldesc_p

  implicit none

  !! Define the SUP_EL namelist
  integer, parameter :: maxTriads_p = 10000 !< Max number of triads
  integer, parameter :: maxGenDofs_p = 1000 !< Max number of component modes

  integer :: id                    !< Base ID of the triad
  integer :: extId(10)             !< User ID path of the superelement
  integer :: numLoadCase           !< Number of load cases
  integer :: numStates             !< Number of states in non-linear reduction
  integer :: numGenDOFs            !< Number of component modes
  integer :: numTriads             !< Number of connected triads
  integer :: triadIds(maxTriads_p) !< Base ID of all connected triads
  integer :: nodeIds(maxTriads_p)  !< Node numbers associated with the triads
  integer :: elPropId              !< ID of the (beam) element properties
  integer :: shadowPosAlg          !< Flag for co-rotated positioning algorithm
  integer :: refTriad1Id           !< Base ID of first reference triad
  integer :: refTriad2Id           !< Base ID of second reference triad
  integer :: refTriad3Id           !< Base ID of third reference triad
  integer :: stressStiffFlag(3)    !< Flags for enabling of geometric stiffness
  integer :: massCorrFlag          !< Flag for enabling mass correction
  integer :: recoveryFlag          !< Flag for enabling stress recovery
  integer :: BC(maxGenDofs_p)      !< Boundary conditions for component modes
  integer :: strDmpEngineId        !< Base ID for structural damping function
  integer :: massEngineId          !< Base ID for mass scaling function
  integer :: stiffEngineId         !< Base ID for stiffness scaling function

  integer :: savePos    !< Flag indicating whether position should be saved
  integer :: saveVar(3) !< Flags indicating which variables should be saved

  character(ldesc_p) :: extDescr      !< User description
  character(lfnam_p) :: inputFiles(7) !< Superelement matrix files
  character(lfnam_p) :: bodyFile      !< Geometry file for hydrodynamics
  character(len=128) :: elmGroups     !< Element groups to do recovery for

  real(dp) :: supPos(4,3) !< Initial position matrix (undeformed state)
  real(dp) :: offset1(3)  !< Offset vector for reference point wrt. refTriad1Id
  real(dp) :: offset2(3)  !< Offset vector for reference point wrt. refTriad2Id
  real(dp) :: offset3(3)  !< Offset vector for reference point wrt. refTriad3Id
  real(dp) :: stiffScale  !< Stiffness scaling factor
  real(dp) :: massScale   !< Mass scaling factor

  real(dp) :: dragParams(3,6) !< Drag parameters
  real(dp) :: slamParams(3)   !< Slamming parameters

  real(dp) :: alpha1               !< Mass proportional damping factor
  real(dp) :: alpha2               !< Stiffness proportional damping factor
  real(dp) :: alpha3(maxGenDofs_p) !< Mass-proportional damping for comp. modes
  real(dp) :: alpha4(maxGenDofs_p) !< Stiff-proportional damping for comp. modes

  namelist /SUP_EL/ id, extId, extDescr, &
       &            inputFiles, bodyFile, elmGroups, elPropId, &
       &            numLoadCase, numStates, numGenDOFs, &
       &            numTriads, triadIds, nodeIds, &
       &            shadowPosAlg, supPos, &
       &            refTriad1Id, offset1, &
       &            refTriad2Id, offset2, &
       &            refTriad3Id, offset3, &
       &            stressStiffFlag, massCorrFlag, recoveryFlag, &
       &            stiffEngineId, massEngineId, stiffScale, massScale, &
       &            dragParams, slamParams, &
       &            strDmpEngineId, alpha1, alpha2, alpha3, alpha4, BC, &
       &            savePos, saveVar

  !! Define the TRIAD_UNDPOS namelist
  integer :: supElId !< Superelement base ID
  integer :: triadId !< Triad base ID

  real(dp) :: undPosInSupElSystem(4,3) !< Undeformed relative triad position
  real(dp) :: eccVec(3)                !< Eccentricity vector for stiffness
  real(dp) :: eccMass(3)               !< Eccentricity vector for mass

  namelist /TRIAD_UNDPOS/ supElId, triadId, undPosInSupElSystem, eccVec, eccMass

contains

  !> @cond FULL_DOC
  !> @brief Reads the SUP_EL namelist from the given file unit.
  subroutine read_SUP_EL (infp,stat)

    integer, intent(in)  :: infp
    integer, intent(out) :: stat

    !! Default values
    id=0; extId=0; extDescr=''
    inputFiles=''; bodyFile=''; elmGroups=''
    numLoadCase=0; numStates=0; numGenDOFs=0
    numTriads=0; triadIds=0; nodeIds=0; elPropId=0
    shadowPosAlg=0; supPos=0.0_dp
    refTriad1Id=0; offset1=0.0_dp
    refTriad2Id=0; offset2=0.0_dp
    refTriad3Id=0; offset3=0.0_dp
    stressStiffFlag=-1; massCorrFlag=-1; recoveryFlag=-1
    stiffScale=1.0_dp; massScale=1.0_dp; dragParams=0.0_dp; slamParams=0.0_dp
    alpha1=0.0_dp; alpha2=0.0_dp; alpha3=0.0_dp; alpha4=0.0_dp
    massEngineId=0; stiffEngineId=0; strDmpEngineId=0
    BC=1; savePos=0; saveVar=0

    read(infp,nml=SUP_EL,iostat=stat)

  end subroutine read_SUP_EL

  !> @brief Reads the TRIAD_UNDPOS namelist from the given file unit.
  subroutine read_TRIAD_UNDPOS (infp,stat)

    integer, intent(in)  :: infp
    integer, intent(out) :: stat

    !! Default values
    supElId=0; triadId=0
    undPosInSupElSystem=0.0_dp
    eccVec=0.0_dp; eccMass=0.0_dp

    read(infp,nml=TRIAD_UNDPOS,iostat=stat)

  end subroutine read_TRIAD_UNDPOS

  !> @brief Produces an error message related to reading TRIAD_UNDPOS namelist.
  subroutine reportTriadError (id,tid,sid,addMsg)

    use IdTypeModule     , only : IdType, getId
    use reportErrorModule, only : reportError, error_p

    type(IdType)              , intent(in) :: id
    integer                   , intent(in) :: tid
    integer         , optional, intent(in) :: sid
    character(len=*), optional, intent(in) :: addMsg

    character(len=128) :: errMsg

    if (present(sid)) then
       write(errMsg,*) 'Invalid TRIAD_UNDPOS specification for Part', &
            &          trim(getId(id)),': supelId =',sid,' triadId =',tid
    else
       write(errMsg,*) 'Could not read TRIAD_UNDPOS namelist number',tid, &
            &          'for Part',trim(getId(id))
    end if

    if (present(addMsg)) then
       call reportError (error_p,adjustl(errMsg),addMsg)
    else
       call reportError (error_p,adjustl(errMsg))
    end if

  end subroutine ReportTriadError
  !> @endcond

end module SupElNamelistModule
