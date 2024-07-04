!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module InaddModule

  implicit none

  private

  integer, save :: unitBad(3) = 0 !< File unit for groups of bad elements

  interface INADD
     module procedure lhsINADD
     module procedure rhsINADD
  end interface

  public :: INADD, extractSubMat


contains

  SUBROUTINE lhsINADD (sam, SSM, SMM, SMASS, RMASS, diagMass, dataCheck, &
       &               IOPadd, IPSW, LPU, IERR, gFull)

    !******************************************************************
    !
    !  F E D E M  (FINITE ELEMENT DYNAMICS IN ELASTIC MECHANISMS)
    !
    !  P R E P R O S E S S O R :   INADD                    S I N T E F
    !
    !
    !     DENNE SUBRUTINEN FORETAR BEREGNING OG INNADDERING AV ELEMENT-
    !     MATRISER (STIVHET/MASSE) I SYSTEMMATRISER FOR EN SUBSTRUKTUR.
    !
    !     IOP = 1 : Assemble stiffness matrix only
    !     IOP = 2 : Assemble mass matrix only
    !     IOP = 3 : Assemble both stiffness and mass
    !
    !     PROGRAMMERT AV : KETIL AAMNES
    !     DATO / VERSJON : 88-10-14 / 1.0
    !
    !******************************************************************

    use KindModule             , only : dp, epsDiv0_p
    use SamModule              , only : SamType
    use SysMatrixTypeModule    , only : SysMatrixType, writeObject
    use AsmExtensionModule     , only : csAddEM, csAddED
    use MatExtensionModule     , only : csGetDiag
    use FileUtilitiesModule    , only : findUnitNumber
    use ManipMatrixModule      , only : writeObject
    use ProgressModule         , only : writeProgress
    use ReportErrorModule      , only : internalError
    use ReportErrorModule      , only : reportError, empty_p, note_p
    use ReportErrorModule      , only : warning_p, error_p, debugFileOnly_p
    use FFlLinkHandlerInterface, only : ffl_getElmId, ffl_getSpring
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_getint
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_getdouble

    type(SamType)      , intent(in)    :: sam
    type(SysMatrixType), intent(inout) :: ssm, smm
    real(dp)           , intent(inout) :: SMASS
    real(dp)           , intent(in)    :: RMASS
    logical            , intent(in)    :: diagMass, dataCheck
    integer            , intent(in)    :: IOPadd, IPSW, LPU
    integer            , intent(out)   :: IERR
    real(dp), optional , intent(out)   :: gFull(:,:)

    !! Local variables
    integer , parameter :: MAXDOF = 120 ! Size of the largest element type
    integer             :: iel, jpsw, iop, iopt, nedof, nErr(3), nMass, nBush
    integer             :: elmShape(2)
    logical             :: addStiff, addMass
    real(dp), parameter :: steel_p = 2.0e11_dp ! Default 'rigid part' stiffness
    real(dp)            :: scaleMass, scaleStiff, kt(3), kr(3), mt(3), mr(3)
    real(dp)            :: EK(MAXDOF**2), EM(MAXDOF**2), EMASS
    character(len=64)   :: errMsg

    !! --- Logic section ---

    call writeProgress(' --> INADD    ; Build system matrices')

    !! Initialization

    iop   = abs(IOPadd)
    jpsw  = IPSW-10
    ierr  = 0
    nErr  = 0
    nMass = 0
    nBush = 0
    addStiff = mod(IOP,2) == 1
    addMass  = IOP > 1
    if (addMass) then
       SMASS = 0.0_dp
       if (present(gFull)) gFull = 0.0_dp
    end if

    !! Start assembly loop over all finite elements

    do iel = 1, sam%nel
       if (IPSW < 0) then ! Debugging of individual elements
          if (-IPSW == ffl_getElmId(iel)) jpsw = 99
       end if

       select case(sam%melcon(iel))
       case(11) !! Two-noded beam
          CALL ADD11(IOP,iel,nedof,ek,em,emass,diagMass,jpsw,lpu,ierr)

       case(21) !! Three-noded thin shell FFT3
          CALL ADD21(IOP,iel,nedof,ek,em,emass,diagMass,jpsw,lpu,ierr)

       case(22) !! Four-noded thin shell FFQ4
          CALL ADD22(IOP,iel,nedof,ek,em,emass,diagMass,jpsw,lpu,ierr)

       case(23) !! Three-noded thin shell ANDEST3
          CALL ADD_ANDES3(IOP,iel,nedof,ek,em,emass,diagMass,lpu,ierr)

       case(24) !! Four-noded thin shell ANDESQ4
          CALL ADD_ANDES4(IOP,iel,nedof,ek,em,emass,diagMass,lpu,ierr)

       case(31) !! Six-noded triangular thick shell
          CALL ADD31(IOP,iel,nedof,ek,em,emass,diagMass,jpsw,lpu,ierr)

       case(32) !! Eight-noded quadrilateral thick shell
          CALL ADD32(IOP,iel,nedof,ek,em,emass,diagMass,jpsw,lpu,ierr)

       case(41) !! Ten-noded tetrahedron
          CALL ADD41(IOP,iel,nedof,ek,em,emass,diagMass,jpsw,lpu,ierr)

       case(42) !! Fifteen-noded wedge
          CALL ADD42(IOP,iel,nedof,ek,em,emass,diagMass,jpsw,lpu,ierr)

       case(43) !! Twenty-noded hexahedron
          CALL ADD43(IOP,iel,nedof,ek,em,emass,diagMass,jpsw,lpu,ierr)

       case(44) !! Eight-noded hexahedron
          CALL ADD44(IOP,iel,nedof,ek,em,emass,diagMass,jpsw,lpu,ierr)

       case(45) !! Four-noded tetrahedron
          CALL ADD45(IOP,iel,nedof,ek,em,emass,diagMass,jpsw,lpu,ierr)

       case(46) !! Six-noded wedge
          CALL ADD46(IOP,iel,nedof,ek,em,emass,diagMass,jpsw,lpu,ierr)

       case(50) !! One-noded mass element with no properties
          if (addMass) nMass = nMass + 1
          goto 100

       case(51) !! One-noded mass element (no stiffness)
          if (.not. addMass) goto 100
          CALL ADD51(sam,iel,nedof,em,emass,diagMass,jpsw,lpu,ierr)

       case(61:63) !! Constraint elements (neither stiffness nor mass)
          goto 100

       case(70) !! Two-noded spring with no properties
          if (addStiff) nBush = nBush + 1
          goto 100

       case(71) !! Two-noded spring element (no mass)
          if (.not. addStiff) goto 100
          call ffl_getSpring (ek,nedof,iel,ierr)
          if (jpsw > -5 .and. ierr == 0) then
             write(lpu,"(//5X,'Stiffness matrix for spring element',I10/)") &
                  &    ffl_getElmId(iel)
             call writeObject (reshape(ek,(/nedof,nedof/)),lpu,eps=1.0e-4_dp)
          end if

       case(72) !! Two-noded eccentric spring (no mass)
          if (.not. addStiff) goto 100
          CALL ADD72(iel,nedof,ek,jpsw,lpu,ierr)

       case default
          write(errMsg,*) sam%melcon(iel)
          ierr = internalError('INADD: Unknown element type:'//errMsg)
          return

       end select

       if (ierr < 0) then
          write(errMsg,"('for element',I10,'  of type',I3)") &
               &        ffl_getElmId(iel), sam%melcon(iel)
          call ReportError (empty_p,errMsg)
          nErr(1) = nErr(1) + 1
          goto 100 ! Go on and try all elements before aborting the program
       else if (ierr > 0) then
          if (unitBad(3) == 0) then
             unitBad(3) = findUnitNumber(10)
             open(unitBad(3),FILE='bad_elements.ftl',STATUS='UNKNOWN')
             write(unitBad(3),"('GROUP{3 ')",ADVANCE='no')
          end if
          write(unitBad(3),'(I8)',ADVANCE='no') ffl_getElmId(iel)
       end if

       elmShape = (/ nedof, nedof /)

       if (IPSW /= 0 .and. sam%melcon(iel) /= 51 .and. addStiff) then
          call checkMatrix (ffl_getElmId(iel),sam%melcon(iel), &
               &            reshape(ek,elmShape),abs(IPSW),lpu,nErr(2:3))
       end if

       if (IPSW == -6) then
          !! Export to ACD test program (for the GSF equation solver)
          if (sam%melcon(iel) /= 51 .and. addStiff) then
             write(lpu,"(/'EK',I8,':',I4)") iel, nedof
             write(lpu,"(1P6E22.15)") ek(1:nedof*nedof)
          end if
          if (sam%melcon(iel) < 60 .and. IOP == 3) then
             write(lpu,"(/'EM',I8,':',I4)") iel, nedof
             if (diagMass) then
                write(lpu,"(1P6E22.15)") em(1:nedof)
             else
                write(lpu,"(1P6E22.15)") em(1:nedof*nedof)
             end if
          end if
       end if

       if (jpsw > -8 .and. sam%melcon(iel) < 60 .and. addMass) then
          if (jpsw > 0) write(lpu,*)
          write(lpu,"(5X,'Element',I10,'  Type',I3,' :  Mass =',1PE13.5)") &
               &     ffl_getElmId(iel), sam%melcon(iel), EMASS
          if (jpsw > 0) then
             if (diagMass) then
                call writeObject (em(1:nedof),lpu, &
                     &            '     Lumped element mass matrix (diag):', &
                     &            eps=1.0e-4_dp)
             else
                call writeObject (reshape(em,elmShape),lpu, &
                     &            '     Element mass matrix:', &
                     &            eps=1.0e-4_dp)
             end if
          end if
       end if

       if (sam%melcon(iel) /= 51 .and. .not. dataCheck .and. addStiff) then
#ifdef FT_DEBUG
          !! Check if the element stiffness matrices have not-a-number values.
          !! This can be due to uninitialized variables, division by zero, etc.
          iopt = count(isNan(ek(1:nedof*nedof)))
          if (iopt > 0) then
             write(errMsg,"('Ignoring element',I10,I3,' due to',I4,' NaN')") &
                  &       ffl_getElmId(iel), sam%melcon(iel), iopt
             call ReportError (error_p,errMsg)
             nErr(1) = nErr(1) + 1
             goto 100 ! Go on and try all elements before aborting the program
          end if
#endif

          !! Assemble the substructure stiffness matrix
          call csAddEM (sam,iel,nedof,reshape(ek,elmShape),ssm,ierr, &
               &        lpu=lpu,ipsw=IPSW)
          if (ierr < 0) goto 890

       end if
       if (sam%melcon(iel) < 60 .and. addMass) then

          if (present(gFull) .and. .not.diagMass) then
             call csAddGravity (iel,reshape(em,elmShape),gFull)
          end if

          !! Calculate the total mass of the substructure
          SMASS = SMASS + EMASS
          if (dataCheck) goto 100
          if (IOPadd < 2) goto 100

          !! Assemble the substructure mass matrix
          if (diagMass) then
             call csAddED (sam,iel,nedof,em(1:nedof),smm,ierr)
          else
             call csAddEM (sam,iel,nedof,reshape(em,elmShape),smm,ierr)
          end if
          if (ierr < 0) goto 890

       end if

100    call WriteProgress (iel,sam%nel)
       if (IPSW < 0 .and. jpsw > 0) jpsw = IPSW-10

    end do
    if (nErr(1) > 0) then
       ierr = -nErr(1)
       write(lpu,"()")
       write(errMsg,"(I8,' erroneous elements detected')") nErr(1)
       goto 900
    end if
    if (nErr(2) > 0) then
       ierr = nErr(2)
       write(lpu,"()")
       write(errMsg,"(I8,' elements with')") nErr(2)
       call ReportError (warning_p,trim(adjustl(errMsg))// &
            &            ' negative diagonal stiffness terms detected')
    end if
    if (nErr(3) > 0) then
       ierr = ierr + nErr(3)
       write(lpu,"()")
       write(errMsg,"(I8,' elements with')") nErr(3)
       call ReportError (warning_p,trim(adjustl(errMsg))// &
            &            ' non-symmetric stiffness matrix detected')
    end if

    if (unitBad(1) > 0) then
       write(unitBad(1),600) 'NEGATIVE DEFINITE'
       close(unitBad(1))
    end if
    if (unitBad(2) > 0) then
       write(unitBad(2),600) 'NON-SYMMETRIC'
       close(unitBad(2))
    end if
    if (unitBad(3) > 0) then
       write(unitBad(3),600) 'BAD'
       close(unitBad(3))
    end if
600 format(' {NAME "',A,' ELEMENTS"}}')

    if (addMass) then
       WRITE(LPU,605) SMASS
605    FORMAT(/4X,'TOTAL MASS OF SUBSTRUCTURE :',1PD14.5)
    end if

    if (nBush+nMass > 0) then

       call writeProgress('     Property-less spring/mass elements')

       !! Property-less spring and/or mass elements are present in the model.
       !! These are typically elements that have been added automatically by
       !! the Fedem GUI as a connection between dependent nodes and Triads.
       !! We need now to assign some stiffness and mass values to these
       !! elements that emulate a rigid and mass-less connection.

       write(errMsg,"(I4,A,I4,A)") nBush,' spring elements and', &
            &                      nMass,' mass elements without properties'
       call ReportError (note_p,trim(adjustl(errMsg))//' were found.', &
            'Automatically computed stiffness/mass values will be assigned.')

       if (nBush > 0) then
          !! Find characteristic stiffness coefficients
          call csGetDiag (sam,ssm,kt,kr,ierr)
          if (ierr < 0) goto 910
          call ffa_cmdlinearg_getint('autoStiffMethod',iopt)
          call ffa_cmdlinearg_getdouble('autoStiffScale',scaleStiff)
          if (kt(3) < epsDiv0_p .or. kr(3) < epsDiv0_p) then
             if (iopt == 1) iopt = 2
             write(errMsg,"('The value',1PE12.5,' is used instead.')") steel_p
             call ReportError (warning_p,'This part lacks translational '// &
                  'and/or rotational FE stiffness completely.',errMsg)
             if (kt(3) < epsDiv0_p) kt = steel_p
             if (kr(3) < epsDiv0_p) kr = steel_p !TODO: scale with L^2 ?
          else if (iopt == 1) then
             call ffa_cmdlinearg_getdouble('tolFactorize',scaleStiff)
             scaleStiff = 0.1_dp / scaleStiff
          end if
          kt(1) = kt(iopt)*scaleStiff
          kr(1) = kr(iopt)*scaleStiff
          write(errMsg,"('Added stiffness terms:',1P2E12.5)") kt(1),kr(1)
          call ReportError (empty_p,errMsg)
       end if

       if (IOPadd < 2) then
          nMass = 0 ! No mass matrix assembly
       else if (nMass > 0 .and. SMASS > epsDiv0_p) then
          !! Find characteristic mass coefficients
          call csGetDiag (sam,smm,mt,mr,ierr)
          if (ierr < 0) goto 910
          call ffa_cmdlinearg_getdouble('autoMassScale',scaleMass)
          if (mr(3) < epsDiv0_p) then
             mr(3) = rMass / sam%nnod
             write(errMsg,"('The value',1PE12.5,' is used instead.')") mr(3)
             call ReportError (warning_p,'This part lacks '// &
                  'rotational FE mass completely.',errMsg)
          end if
          mt(1) = mt(3)*scaleMass
          mr(1) = mr(3)*scaleMass
          write(errMsg,"('Added mass terms     :',1P2E12.5)") mt(1),mr(1)
          call ReportError (empty_p,errMsg)
       else if (nMass > 0) then
          nMass = 0
          call ReportError (warning_p,'Mass-less part. No masses are added.')
       end if

       do iel = 1, sam%nel
          select case(sam%melcon(iel))

          case(50) !! One-noded mass element
             if (nMass == 0) cycle
             CALL ADD50 (sam,iel,nedof,em,mt(1),mr(1),diagMass,jpsw,lpu,ierr)
             if (ierr < 0) goto 910
             if (diagMass) then
                call csAddED (sam,iel,nedof,em(1:nedof),smm,ierr)
             else
                elmShape = (/ nedof, nedof /)
                call csAddEM (sam,iel,nedof,reshape(em,elmShape),smm,ierr)
                if (present(gFull)) then
                   call csAddGravity (iel,reshape(em,elmShape),gFull)
                end if
             end if
             if (ierr < 0) goto 890

          case(70) !! Two-noded spring
             if (nBush == 0) cycle
             CALL ADD70 (iel,nedof,ek,kt(1),kr(1),jpsw,lpu,ierr)
             if (ierr < 0) goto 910
             elmShape = (/ nedof, nedof /)
             call csAddEM (sam,iel,nedof,reshape(ek,elmShape),ssm,ierr)
             if (ierr < 0) goto 890

          end select
       end do

    end if

    if (sam%mpar(25) <= 41 .and. sam%mpar(26) >= 41 .and. addStiff) then
       call writeObject (ssm,sam%mpar,lpu, &
            &            ' >>> Assembled stiffness matrix <<<', &
            &            complexity=IPSW)
    end if
    if (sam%mpar(25) <= 42 .and. sam%mpar(26) >= 42 .and. IOPadd > 1) then
       call writeObject (smm,sam%mpar,lpu, &
            &            ' >>> Assembled mass matrix <<<', &
            &            complexity=IPSW)
    end if
    if (sam%mpar(25) <= 44 .and. sam%mpar(26) >= 44 .and. present(gFull)) then
       call writeObject (gFull,lpu,'Gravitation force vectors')
    end if

    return

890 write(errMsg,"('Matrix assembly failed for element',I10)") ffl_getElmId(iel)
900 call ReportError (error_p,adjustl(errMsg))
910 call ReportError (debugFileOnly_p,'INADD')

  contains

    subroutine csAddGravity (iel,EM,g)
      integer, intent(in)     :: iel
      real(dp), intent(in)    :: EM(:,:)
      real(dp), intent(inout) :: g(:,:)

      integer  :: i, j, inpc, idof, jdof, ldof, ieq, nndof

      ldof = 1
      nndof = nedof / (sam%mpmnpc(iel+1)-sam%mpmnpc(iel))
      do inpc = sam%mpmnpc(iel), sam%mpmnpc(iel+1)-1
         idof = sam%madof(sam%mmnpc(inpc))
         do i = 1, nndof
            ieq = sam%meqn(idof)
            if (ieq > 0) then
               do jdof = 1, nedof
                  j = 1 + mod(jdof-1,nndof)
                  if (j <= 3) then
                     g(ieq,j) = g(ieq,j) + EM(ldof,jdof)
                  end if
               end do
            end if
            ldof = ldof + 1
            idof = idof + 1
         end do
      end do

    end subroutine csAddGravity

  END SUBROUTINE lhsINADD


  subroutine rhsINADD (sam, RHS, ILC, IPSW, LPU, IERR)

    !!==========================================================================
    !! Assemble the system load vector for the external load case ILC.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 2 Apr 2008 / 1.0
    !!==========================================================================

    use SamModule              , only : SamType, dp
    use AsmExtensionModule     , only : csAddEV, csAddNV
    use ProgressModule         , only : writeProgress
    use ReportErrorModule      , only : reportError
    use ReportErrorModule      , only : empty_p, error_p, debugFileOnly_p
    use FFlLinkHandlerInterface, only : ffl_getNoLoad, ffl_getLoad, ffl_getElmId

    type(SamType), intent(in)  :: sam
    real(dp)     , intent(out) :: RHS(:)
    integer      , intent(in)  :: ILC, IPSW, LPU
    integer      , intent(out) :: IERR

    !! Local variables
    integer, parameter :: MAXDOF = 120 ! Size of the largest element type
    integer, parameter :: MAXNOD = 20
    integer            :: iel, face, iLoad, nLoad, nedof, nErr
    real(dp)           :: ES(MAXDOF), P(3*MAXNOD)
    character(len=64)  :: errMsg

    !! --- Logic section ---

    ierr  = 0
    nErr  = 0
    rhs   = 0.0_dp
    nLoad = ffl_getNoLoad(ilc)
    if (nLoad < 1) return ! empty load case

    write(errMsg,"(I8)") ilc
    call writeProgress (' --> INADD    ; Build system load vector '// &
         &              trim(adjustl(errMsg)))

    !! Loop over all loads in current load case
    iLoad = 0
    call ffl_getLoad (ilc,iel,face,p(1))
    do while (iel > 0)
       iLoad = iLoad + 1
       if (face > 0) then

          !! Element load
          select case (sam%melcon(iel))

          case(21) !! Three-noded thin shell
             call LOAD21 (iel,p,nedof,es,ipsw,lpu,ierr)

          case(22) !! Four-noded thin shell
             call LOAD22 (iel,p,nedof,es,ipsw,lpu,ierr)

          case(23) !! Three-noded ANDES shell (TODO: separate load routine?)
             call LOAD21 (iel,p,nedof,es,ipsw,lpu,ierr)

          case(24) !! Four-noded ANDES shell (TODO: separate load routine?)
             call LOAD22 (iel,p,nedof,es,ipsw,lpu,ierr)

          case(41) !! Ten-noded tetrahedron
             call LOAD41 (iel,face,p,nedof,es,ipsw,lpu,ierr)

          case(42) !! Fifteen-noded wedge
             call LOAD42 (iel,face,p,nedof,es,ipsw,lpu,ierr)

          case(43) !! Twenty-noded hexahedron
             call LOAD43 (iel,face,p,nedof,es,ipsw,lpu,ierr)

          case(44) !! Eight-noded hexahedron
             call LOAD44 (iel,face,p,nedof,es,ipsw,lpu,ierr)

          case(45) !! Four-noded tetrahedron
             call LOAD45 (iel,face,p,nedof,es,ipsw,lpu,ierr)

          case(46) !! Six-noded wedge
             call LOAD46 (iel,face,p,nedof,es,ipsw,lpu,ierr)

          case default
             ierr = -1
             call ReportError (error_p,'Element load not available')
          end select

          if (ierr < 0) then
             write(errMsg,"('for element',I10,'  of type',I3)") &
                  &        ffl_getElmId(iel), sam%melcon(iel)
             call ReportError (empty_p,errMsg)
             nErr = nErr + 1
          else
             !! Add the element load vector into the right-hand-side vector
             call csAddEV (sam,iel,nedof,es,rhs,rhs,ierr)
             if (ierr < 0) then
                write(errMsg,"('Load vector assembly failed for element',I10)")&
                     &        ffl_getElmId(iel)
                goto 900
             end if
          end if

       else if (face < 0) then

          !! Nodal point load
          if (face == -1) then
             es(1:3) = p(1:3) ! force only
             es(4:6) = 0.0_dp ! zero out the moment part
          else if (face == -2) then
             es(4:6) = p(1:3) ! moment only
             es(1:3) = 0.0_dp ! zero out the force part
          else
             es(1:6) = p(1:6) ! force and moment
          end if

          !! Add the load directly into the right-hand-side vector
          call csAddNV (sam,iel,es,rhs,rhs,ierr)
          if (ierr < 0) then
             write(errMsg,"('Load vector assembly failed for node',I10)") &
                  &        sam%minex(iel)
             goto 900
          end if

       end if

       call WriteProgress (iLoad,nLoad)
       call ffl_getLoad (ilc,iel,face,p(1))
    end do
    if (nErr > 0) then
       ierr = -nErr
       write(errMsg,"(I8,' erroneous element loads detected')") nErr
       goto 900
    end if

    return

900 call ReportError (error_p,adjustl(errMsg))
    call ReportError (debugFileOnly_p,'INADD')

  end subroutine rhsINADD


  subroutine checkMatrix (iel,ieltyp,ek,ipsw,lpu,nErr)

    !!==========================================================================
    !! Check that the matrix is symmetric and have non-negative diagonal terms.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 10 Mar 2004 / 1.0
    !!==========================================================================

    use kindModule         , only : dp
    use fileUtilitiesModule, only : findUnitNumber

    integer , intent(in)    :: iel, ieltyp, ipsw, lpu
    real(dp), intent(in)    :: ek(:,:)
    integer , intent(inout) :: nErr(2)

    !! Local variables
    real(dp), parameter :: eps_p = 1.0e-14_dp
    integer             :: i, j, nTerms

    !! --- Logic section ---

    nTerms = 0
    do i = 1, size(ek,1)
       if (ek(i,i) < 0.0_dp) then
          nTerms = nTerms + 1
          if (ipsw > 2) write(lpu,690) i,ek(i,i)
       end if
    end do
    if (nTerms > 0) then
       nErr(1) = nErr(1) + 1
       write(lpu,691) nTerms,iel,ieltyp
       nTerms = 0
       if (unitBad(1) == 0) then
          unitBad(1) = findUnitNumber(10)
          open(unitBad(1),FILE='negpiv_elements.ftl',STATUS='UNKNOWN')
          write(unitBad(1),"('GROUP{1 ')",ADVANCE='no')
       end if
       write(unitBad(1),"(I8)",ADVANCE='no') iel
    end if

    do i = 1, size(ek,1)
       do j = i+1, size(ek,1)
          if (abs(ek(i,j)-ek(j,i)) > eps_p*max(abs(ek(i,j)),abs(ek(j,i)))) then
             nTerms = nTerms + 1
             if (ipsw > 2) write(lpu,692) i,j,ek(i,j),ek(j,i),ek(i,j)-ek(j,i)
          end if
       end do
    end do
    if (nTerms > 0) then
       nErr(2) = nErr(2) + 1
       write(lpu,693) nTerms,iel,ieltyp
       if (unitBad(2) == 0) then
          unitBad(2) = findUnitNumber(10)
          open(unitBad(2),FILE='nonsym_elements.ftl',STATUS='UNKNOWN')
          write(unitBad(2),"('GROUP{2 ')",ADVANCE='no')
       end if
       write(unitBad(2),'(I8)',ADVANCE='no') iel
    end if

690 format(5X,'Dof',i3,': Diag = ',1PE12.5)
691 format(I6,' negative diagonal term(s) found in element',I8,' of type',I4)
692 format(5X,'Dofs',2i3,': Upper, Lower, diff =',1P2E18.10,E13.5)
693 format(I6,' non-symmetric term(s) found in element',I8,' of type',I4)

  end subroutine checkMatrix


  subroutine extractSubMat (sam, SSM, SMM, SMIE, SMEE, SKIE, SKEE, &
       &                    IPSW, LPU, IERR, extractSKie)

    !!==========================================================================
    !! Extract submatrices 12 and 22 from the system mass- and stiffness matrix.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 3 Oct 2002 / 1.0 (extracted from INADD above)
    !!==========================================================================

    use SamModule          , only : SamType, dp
    use SysMatrixTypeModule, only : SysMatrixType, diagonalMatrix_p, outOfCore_p
    use SparseMatrixModule , only : SparseMatrixType, smNullify, smAllocate
    use SparseMatrixModule , only : smShrinkToFit, smWrite
    use MatExtensionModule , only : csGetSub12, csGetSub22
    use ProgressModule     , only : writeProgress
    use ReportErrorModule  , only : internalError, reportError, debugFileOnly_p

    type(SamType)         , intent(in)    :: sam
    type(SysMatrixType)   , intent(inout) :: SSM, SMM
    type(SparseMatrixType), intent(out)   :: SMIE, SKIE
    real(dp)              , intent(out)   :: SMEE(:,:), SKEE(:,:)
    integer               , intent(in)    :: IPSW, LPU
    integer               , intent(out)   :: IERR
    logical, optional     , intent(in)    :: extractSKie ! Usually not needed

    !! --- Logic section ---

    if (sam%ndof2 <= 0) then

       ierr = internalError('extractSubMat: NDOF2 is zero (no ext. DOFs)')
       return

    else if (sam%ndof1 <= 0) then

       call writeProgress('     Extracting submatrices Kee and Mee')

       call smNullify (skie)
       call smNullify (smie)

    else

       call writeProgress('     Extracting submatrices Kie,Mie and Kee,Mee')

       if (.not. present(extractSKie)) then
          call smNullify (skie)
       else if (.not. extractSKie) then
          call smNullify (skie)
       else

          !! Extract submatrix Kie from system matrix SK

          call smAllocate (skie,sam%ndof1,sam%ndof2, &
               &           ssm%storageType /= outOfCore_p,ierr)
          if (ierr /= 0) goto 900

          call csGetSub12 (ssm,skie,ierr)
          if (ierr /= 0) goto 900

          call smShrinkToFit (skie,ierr)
          call smWrite (skie,'K_ie',lpu, ipsw > 5 .and. ierr == 0)
          if (ierr /= 0) goto 900

       end if
       if (smm%storageType == diagonalMatrix_p .or. size(smee) < 1) then
          call smNullify (smie)
       else

          !! Extract submatrix Mie from system matrix SM

          call smAllocate (smie,sam%ndof1,sam%ndof2, &
               &           smm%storageType /= outOfCore_p,ierr)
          if (ierr /= 0) goto 900

          call csGetSub12 (smm,smie,ierr)
          if (ierr /= 0) goto 900

          call smShrinkToFit (smie,ierr)
          call smWrite (smie,'M_ie',lpu, ipsw > 5 .and. ierr == 0)
          if (ierr /= 0) goto 900

       end if

    end if

    !! Extract submatrix Kee from system matrix SK

    call csGetSub22 (ssm,sam%meqn2,skee,ierr)
    if (ierr /= 0) goto 900

    !! Extract submatrix Mee from system matrix SM

    if (size(smee) > 0) then
       call csGetSub22 (smm,sam%meqn2,smee,ierr)
       if (ierr /= 0) goto 900
    end if

    call writeProgress('     Done extracting submatrices')

    return

900 call ReportError (debugFileOnly_p,'extractSubMat')

  end subroutine extractSubMat


  function dsum (n,A,inc)
    use KindModule, only : dp
    integer , intent(in) :: n, inc
    real(dp), intent(in) :: A(n)
    real(dp) :: dsum
    integer  :: i
    dsum = 0.0_dp
    do i = 1, n, inc
       dsum = dsum + A(i)
    end do
  end function dsum

  function lumpSolidMass (n,MC,ML) result(eMass)
    use KindModule, only : dp
    integer , intent(in) :: n
    real(dp), intent(inout) :: MC(n*n)
    real(dp), intent(inout) :: ML(n)
    real(dp) :: eMass
    integer  :: i, j, k
    j = 1
    k = 0
    do i = 1, n
       k = k + n
       ML(i) = sum(MC(j:k))
       j = j + n
    end do
    eMass = dsum(n,ML,3)
  end function lumpSolidMass

  subroutine expand56 (n,V5,V6)
    use KindModule, only  : dp
    integer , intent(in)  :: n
    real(dp), intent(in)  :: V5(:)
    real(dp), intent(out) :: V6(:)
    integer :: i, j
    j = 1
    do i = 1, n, 6
       V6(i:i+4) = V5(j:j+4)
       V6(i+5) = 0.0_dp
       j = j + 5
    end do
  end subroutine expand56


  ! ######################################################################
  ! ######################################################  ADD11      ###
  ! ######################################################################

  SUBROUTINE ADD11(IOP, IEL, NEDOF, EK, EM, EMASS, diagMass, IPSW, LPU, IERR)

    ! **********************************************************************
    !
    !   F E D E M  (FINITE ELEMENT DYNAMICS IN ELASTIC MECHANISMS)
    !
    !   P R E P R O S E S S O R :   ADD11                        S I N T E F
    !
    !
    !     DENNE SUBRUTINEN BEREGNER ELEMENTMATRISER (STIVHET/MASSE)
    !     FOR BJELKEELEMENTET BEAM (ELTYPE = 11).
    !     DEN BEREGNER OGS] ELEMENTETS MASSE.
    !
    !     PROGRAMMERT AV : KETIL AAMNES
    !     DATO / VERSJON : 88-10-14 / 1.0
    !
    ! **********************************************************************

    use KindModule             , only : dp, epsDiv0_p
    use ReportErrorModule      , only : reportError, error_p
    use FFlLinkHandlerInterface, only : ffl_getCoor, ffl_getElmId, ffl_getNSM
    use FFlLinkHandlerInterface, only : ffl_getPinFlags, ffl_getBeamSection

    integer , parameter   :: neldof = 12
    integer , intent(in)  :: iop, iel, ipsw, lpu
    logical , intent(in)  :: diagMass
    real(dp), intent(out) :: ek(neldof,neldof), em(neldof,neldof), emass
    integer , intent(out) :: nedof, ierr

    !! Local variables
    integer           :: IPA, IPB
    real(dp)          :: XG(5), YG(5), ZG(5), BSEC(14), NSMASS
    character(len=64) :: errMsg

    !! --- Logic section ---

    nedof = neldof
    em    = 0.0_dp
    ek    = 0.0_dp

    ! --- HENTER N\DVENDIGE ELEMENTDATA

    call ffl_getCoor (XG,YG,ZG,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element nodal coordinates')
       return
    end if

    call ffl_getBeamSection (BSEC,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving beam element data')
       return
    end if

    call ffl_getNSM (NSMass,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element non-structural mass')
       return
    end if

    ! Adjust material density for non-structural mass
    if (BSEC(4) > epsDiv0_p) BSEC(1) = BSEC(1) + NSMass/BSEC(4)

    ! Get the beam pin flags, if any
    call ffl_getPinFlags (IPA,IPB,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving beam pin flags')
       return
    end if

    ! --- BEREGNER STIVHETSMATRISE FOR ELEMENTET

    if (ipsw > -7) write(lpu,991) ffl_getElmId(IEL),BSEC

    if (mod(iop,2) == 0) goto 100
    CALL BEAM31(EK(1,1),XG(1),YG(1),ZG(1),BSEC(2),BSEC(9),BSEC(11),BSEC(13), &
         &      BSEC(14),IPA,IPB,IEL,LPU,IPSW,IERR)
    if (ierr /= 0) then
       if (BSEC(3) <= 1.0e-16_dp) then
          write(errMsg,"('A zero shear modulus G =',1PE12.5)") BSEC(3)
          call ReportError (error_p,trim(errMsg)//' was detected for a beam.', &
               'This is not allowed. If you want to model a beam with no '// &
               'shear stiffness,','set the transverse shear-reduction '// &
               'factors for the associated beam property to zero instead.')
       end if
       call ReportError (error_p,'computing beam stiffness matrix')
       return
    end if

    ! --- BEREGNER MASSEMATRISE FOR ELEMENTET

100 if (iop < 2) return
    CALL BEAM35(EM(1,1),XG(1),YG(1),ZG(1),BSEC(4),BSEC(1),BSEC(8),2, &
         &      IPA,IPB,LPU,IPSW,IERR)
    EMASS = EM(1,1)+EM(7,1)

    if (ipsw > -7) write(lpu,992) EM(:,1)

    if (.not. diagMass) then
       CALL BEAM35(EM(1,1),XG(1),YG(1),ZG(1),BSEC(4),BSEC(1),BSEC(8),3, &
            &      IPA,IPB,LPU,IPSW,IERR)
    end if

991 FORMAT(/5X,'Beam element',I10/8X,'BSEC  =',1P,6E13.5/(15X,6E13.5))
992 FORMAT( 8X,'LMASS =',1P,12E13.5)

  END SUBROUTINE ADD11


  ! ######################################################################
  ! ######################################################  ADD21      ###
  ! ######################################################################

  SUBROUTINE ADD21(IOP, IEL, NEDOF, EK, EM, EMASS, diagMass, IPSW, LPU, IERR)

    ! **********************************************************************
    !
    !   F E D E M  (FINITE ELEMENT DYNAMICS IN ELASTIC MECHANISMS)
    !
    !   P R E P R O S E S S O R :   ADD21                        S I N T E F
    !
    !
    !     DENNE SUBRUTINEN BEREGNER ELEMENTMATRISER (STIVHET/MASSE)
    !     FOR TYNNSKALLELEMENTET FTS (ELTYPE = 21).
    !     DEN BEREGNER OGS] ELEMENTETS MASSE.
    !
    !     PROGRAMMERT AV : KETIL AAMNES
    !     DATO / VERSJON : 88-10-14 / 1.0
    !
    ! **********************************************************************

    use KindModule             , only : dp, epsDiv0_p
    use IsoMatModule           , only : isoMat2D
    use ReportErrorModule      , only : reportError, error_p
    use FFlLinkHandlerInterface, only : ffl_getCoor, ffl_getMat, ffl_getThick
    use FFlLinkHandlerInterface, only : ffl_getNSM

    integer , parameter   :: neldof = 18
    integer , intent(in)  :: iop, iel, ipsw, lpu
    logical , intent(in)  :: diagMass
    real(dp), intent(out) :: ek(neldof,neldof), em(neldof,neldof), emass
    integer , intent(out) :: nedof, ierr

    !! Local variables
    real(dp), parameter :: alpha = 1.5_dp ! Hardcoded from the old EMAT
    real(dp), parameter :: beta  = 0.5_dp ! Hardcoded from the old EMAT

    real(dp) :: EMV(neldof), XG(3), YG(3), ZG(3), E(3,3)
    real(dp) :: TH(3), HH(3,9), AKB(9,9), EMOD, RNY, RHO, NSMASS

    !! --- Logic section ---

    nedof = neldof
    em    = 0.0_dp
    ek    = 0.0_dp

    ! --- HENTER N\DVENDIGE ELEMENTDATA

    call ffl_getcoor (XG,YG,ZG,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element nodal coordinates')
       return
    end if

    call ffl_getmat (EMOD,RNY,RHO,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element material data')
       return
    end if

    call ffl_getthick (TH,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element thickness')
       return
    end if

    call ffl_getNSM (NSMass,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element non-structural mass')
       return
    end if

    ! Adjust material density for non-structural mass
    if (TH(1) > epsDiv0_p) RHO = RHO + NSMass/TH(1)

    if (mod(iop,2) == 0) goto 100
    call isoMat2D (EMOD,RNY,E)

    ! --- BEREGNER ELEMENTETS STIVHETSMATRISE

    CALL FTS31(EK(1,1),HH(1,1),AKB(1,1),E(1,1),XG(1),YG(1),ZG(1),TH(1), &
         &     ALPHA,BETA,IEL,LPU,IPSW,IERR)
    IF (IERR /= 0) THEN
       call ReportError (error_p,'computing shell stiffness matrix')
       return
    END IF

    ! --- BEREGNER ELEMENTETS MASSEMATRISE

100 if (iop < 2) return
    CALL FTS35(EMV(1),XG(1),YG(1),ZG(1),TH(1),RHO)

    if (diagMass) then
       call DCOPY(NEDOF,EMV(1),1,EM(1,1),1)
    else
       call DCOPY(NEDOF,EMV(1),1,EM(1,1),NEDOF+1)
    end if

    ! --- BEREGNER ELEMENTETS MASSE

    EMASS = dsum(NEDOF,EMV,6)

  END SUBROUTINE ADD21


  ! ######################################################################
  ! ######################################################  ADD22      ###
  ! ######################################################################

  SUBROUTINE ADD22(IOP, IEL, NEDOF, EK, EM, EMASS, diagMass, IPSW, LPU, IERR)

    ! **********************************************************************
    !
    !   F E D E M  (FINITE ELEMENT DYNAMICS IN ELASTIC MECHANISMS)
    !
    !   P R E P R O S E S S O R :   ADD22                        S I N T E F
    !
    !
    !     DENNE SUBRUTINEN BEREGNER ELEMENTMATRISER (STIVHET/MASSE)
    !     FOR TYNNSKALLELEMENTET FQS (ELTYPE = 22).
    !     DEN BEREGNER OGS] ELEMENTETS MASSE.
    !
    !     PROGRAMMERT AV : TERJE R\LV]G
    !     DATO / VERSJON : 88-10-14 / 1.0
    !
    ! **********************************************************************

    use KindModule             , only : dp, epsDiv0_p
    use IsoMatModule           , only : isoMat2D
    use PmatModule             , only : pMatStiff
    use ReportErrorModule      , only : reportError, error_p
    use FFlLinkHandlerInterface, only : ffl_getCoor, ffl_getMat, ffl_getThick
    use FFlLinkHandlerInterface, only : ffl_getNSM

    integer , parameter   :: neldof = 24
    integer , intent(in)  :: iop, iel, ipsw, lpu
    logical , intent(in)  :: diagMass
    real(dp), intent(out) :: ek(neldof,neldof), em(neldof,neldof), emass
    integer , intent(out) :: nedof, ierr

    !! Local variables
    real(dp), parameter :: alpha = 1.5_dp ! Hardcoded from the old EMAT
    real(dp), parameter :: beta  = 0.9_dp ! Hardcoded from the old EMAT

    real(dp) :: EMV(neldof), XG(4), YG(4), ZG(4), E(3,3), TH(4)
    real(dp) :: EMOD, RNY, RHO, NSMASS, PMAT(neldof,neldof)

    !! --- Logic section ---

    nedof = neldof
    em    = 0.0_dp
    ek    = 0.0_dp

    ! --- HENTER N\DVENDIGE ELEMENTDATA

    call ffl_getcoor (XG,YG,ZG,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element nodal coordinates')
       return
    end if

    call ffl_getmat (EMOD,RNY,RHO,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element material data')
       return
    end if

    call ffl_getthick (TH,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element thickness')
       return
    end if

    call ffl_getNSM (NSMass,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element non-structural mass')
       return
    end if

    ! Adjust material density for non-structural mass
    if (TH(1) > epsDiv0_p) RHO = RHO + NSMass/TH(1)

    if (mod(iop,2) == 0) goto 100
    call isoMat2D (EMOD,RNY,E)

    ! --- BEREGNER ELEMENTETS STIVHETSMATRISE

    CALL FQS31(EM(1,1),E(1,1),XG(1),YG(1),ZG(1),TH(1),ALPHA,BETA,IEL,LPU, &
         &     IPSW,IERR,EMOD,RNY)
    IF (IERR < 0) THEN
       call ReportError (error_p,'computing shell stiffness matrix')
       return
    END IF

    ! Compute projection matrix that restores stress free rigid body motions
    CALL PMATSTIFF(xg,yg,zg,PMAT,LPU,IERR)
    IF (IERR < 0) THEN
       call ReportError (error_p,'computing projection matrix')
       return
    END IF
    ! Project the stiffness matrix: K = P'*K*P
    CALL SYMTRA(EM(1,1),PMAT(1,1),EK(1,1),EMV(1),NEDOF,NEDOF)

    ! --- BEREGNER ELEMENTETS MASSEMATRISE

100 if (iop < 2) return
    CALL FQS35(EMV(1),XG(1),YG(1),ZG(1),TH(1),RHO)
    if (diagMass) then
       call DCOPY(NEDOF,EMV(1),1,EM(1,1),1)
    else
       EM = 0.0_dp
       call DCOPY(NEDOF,EMV(1),1,EM(1,1),NEDOF+1)
    end if

    ! --- BEREGNER ELEMENTETS MASSE (SUMMERER BIDRAG I X-RETNING)

    EMASS = dsum(NEDOF,EMV,6)

  END SUBROUTINE ADD22

  ! ######################################################################
  ! ######################################################  ADD_ANDES3 ###
  ! ######################################################################

  SUBROUTINE ADD_ANDES3(IOP, IEL, NEDOF, &
       &                EK, EM, EMASS, diagMass, LPU, IERR)

    ! **********************************************************************
    !   Stiffness and mass matrix for ANDEST3 element.
    !
    !   Coded by Bjorn Haugen
    ! **********************************************************************

    use KindModule             , only : dp, epsDiv0_p
    use Andes3ShellModule      , only : Andes3shell_stiffmat
    use CompositeTypeModule    , only : CompositeType, Cmatrix
    use ManipMatrixModule      , only : trans3p
    use ReportErrorModule      , only : reportError, error_p
    use FFlLinkHandlerInterface, only : ffl_getCoor, ffl_getElmId, ffl_getNSM
    use FFlLinkHandlerInterface, only : ffl_getPMATSHELL, ffl_getThick
    use FFlLinkHandlerInterface, only : ffl_hasAttribute, ffl_getAttributeID

    integer , parameter   :: neldof = 18
    integer , intent(in)  :: iop, iel, lpu
    logical , intent(in)  :: diagMass
    real(dp), intent(out) :: ek(neldof,neldof), em(neldof,neldof), emass
    integer , intent(out) :: nedof, ierr

    !! Local variables
    integer , parameter :: ltype  = 2
    real(dp), parameter :: alpha  = 1.5_dp
    real(dp), parameter :: alphaH = 0.5_dp

    integer  :: i, j, compID
    real(dp) :: EMV(neldof), XG(3), YG(3), ZG(3), P1(3), P2(3), P3(3)
    real(dp) :: TH(3), EMOD, E2, G12, G1Z, G2Z, RNY, RHO, NSMASS

    type(CompositeType) :: laminate

    real(dp) :: Kmat(18,18), Cmat(6,6), thetaMaterial
    real(dp) :: xl(3), yl(3), Trel(3,4), Tel(3,3), TelT(3,3)
    real(dp) :: X21, Y21, Z21, X31, Y31, Z31, SL21, SL31, COSG, SING

    !! --- Logic section ---

    nedof = neldof
    em    = 0.0_dp
    ek    = 0.0_dp

    ! --- HENTER N\DVENDIGE ELEMENTDATA

    call ffl_getcoor (XG,YG,ZG,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element nodal coordinates')
       return
    end if

    X21 = XG(2)-XG(1)
    Y21 = YG(2)-YG(1)
    Z21 = ZG(2)-ZG(1)
    X31 = XG(3)-XG(1)
    Y31 = YG(3)-YG(1)
    Z31 = ZG(3)-ZG(1)
    SL21 = sqrt(X21*X21 + Y21*Y21 + Z21*Z21)
    SL31 = sqrt(X31*X31 + Y31*Y31 + Z31*Z31)
    COSG = (X31*X21 + Y31*Y21 + Z31*Z21)/(SL31*SL21)
    if (COSG > -1.0_dp .and. COSG < 1.0_dp) then
       SING = sqrt(1.0_dp - COSG*COSG)
    else
       SING = 0.0_dp
    end if

    XL(1) = 0.0_dp
    XL(2) = SL21
    XL(3) = SL31*COSG
    YL(1) = 0.0_dp
    YL(2) = 0.0_dp
    YL(3) = SL31*SING

    !!=== Element transformation matrix

    P1(1)=XG(1)  !defining points instead of x, y and z vectors
    P1(2)=YG(1)
    P1(3)=ZG(1)

    P2(1)=XG(2)
    P2(2)=YG(2)
    P2(3)=ZG(2)

    P3(1)=XG(3)
    P3(2)=YG(3)
    P3(3)=ZG(3)

    trel = trans3P(P1, P2, P3, lpu, ierr)
    if (ierr < 0) then
       write(lpu,600) ffl_getElmId(iel)
       ierr = -ierr
600    format('  ** The 3-noded shell element',I6, &
            & ' has zero area and is ignored.' &
            / '     Check your FE model for consistency.')
       return
    end if

    TelT = trel(1:3,1:3)
    Tel  = transpose(TelT)

    if (ffl_hasAttribute('PMAT',iel,ierr)) then

       !!=== Standard isotropic material and thickness

       if (ierr < 0) then
          call ReportError (error_p,'Error retrieving element material data')
          return
       end if

       call ffl_getmat (EMOD,RNY,RHO,IEL,ierr)
       if (ierr < 0) then
          call ReportError (error_p,'Error retrieving element material data')
          return
       end if

       call ffl_getthick (TH,IEL,ierr)
       if (ierr < 0) then
          call ReportError (error_p,'Error retrieving element thickness')
          return
       end if

       !Form resulting shell resultant form constitutive matrix here

       laminate%nPly = 1
       laminate%PID  = 0
       allocate(laminate%T(1),laminate%theta(1),laminate%pMat(1),STAT=ierr)

       laminate%T(1)         = (TH(1)+TH(2)+TH(3)) / 3.0_dp
       laminate%Z0           = -0.5_dp*laminate%T(1)
       laminate%theta(1)     = 0.0_dp
       laminate%pMat(1)%E1   = EMOD
       laminate%pMat(1)%E2   = EMOD
       laminate%pMat(1)%NU12 = RNY
       laminate%pMat(1)%G12  = 0.5_dp * EMOD / (1.0_dp+RNY)
       laminate%pMat(1)%RHO  = RHO

       thetaMaterial = 0.0_dp

    else if (ffl_hasAttribute('PMATSHELL',iel,ierr)) then

       !!=== Single layer anisotropic material

       if (ierr < 0) then
          call ReportError (error_p,'Error retrieving shell material data')
          return
       end if

       compID = ffl_getAttributeID('PMATSHELL',iel)
       call ffl_getPMATSHELL (compID, EMOD, E2, RNY, G12, G1Z, G2Z, RHO, ierr)
       if (ierr < 0) then
          call ReportError (error_p,'Error retrieving shell material data')
          return
       end if

       call ffl_getthick (TH,IEL,ierr)
       if (ierr < 0) then
          call ReportError (error_p,'Error retrieving element thickness')
          return
       end if

       !Form resulting shell resultant form constitutive matrix here

       laminate%nPly = 1
       laminate%PID  = compID
       allocate(laminate%T(1),laminate%theta(1),laminate%pMat(1),STAT=ierr)

       laminate%T(1)         = (TH(1)+TH(2)+TH(3)) / 3.0_dp
       laminate%Z0           = -0.5_dp*laminate%T(1)

       laminate%theta(1)     = 0.0_dp
       laminate%pMat(1)%E1   = EMOD
       laminate%pMat(1)%E2   = E2
       laminate%pMat(1)%NU12 = RNY
       laminate%pMat(1)%G12  = G12
       laminate%pMat(1)%RHO  = RHO

       if (abs(G12) > epsDiv0_p) then
          laminate%pMat(1)%G12 = G12
       else
          laminate%pMat(1)%G12 = 0.25_dp * (EMOD+E2) / (1.0_dp+RNY)
       end if

       thetaMaterial = 0.0_dp

    else

       !!=== Composite material, TODO later...
       ierr = -999
       call ReportError (error_p,'Composite material not yet implemented')
       return
    end if

    call ffl_getNSM (NSMass,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element non-structural mass')
       return
    end if

    ! Adjust material density for non-structural mass
    if (laminate%T(1) > epsDiv0_p) RHO = RHO + NSMass/laminate%T(1)

    if (mod(iop,2) == 0) goto 100

    ! --- BEREGNER ELEMENTETS STIVHETSMATRISE

    call Cmatrix (laminate,thetaMaterial,Cmat)
    call Andes3shell_stiffmat (xl,yl,Cmat,alpha,alphaH,ltype,Kmat,LPU,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'computing shell stiffness matrix')
       return
    end if

    do i = 1,18,3
       do j = 1,18,3
          EK(i:i+2,j:j+2) = matmul(TelT,matmul(Kmat(i:i+2,j:j+2),Tel))
       end do
    end do

    ! --- BEREGNER ELEMENTETS MASSEMATRISE

100 deallocate(laminate%T,laminate%theta,laminate%pMat)
    if (iop < 2) return

    CALL FTS35(EMV(1),XG(1),YG(1),ZG(1),TH(1),RHO)

    if (diagMass) then
       call DCOPY(NEDOF,EMV(1),1,EM(1,1),1)
    else
       call DCOPY(NEDOF,EMV(1),1,EM(1,1),NEDOF+1)
    end if

    ! --- BEREGNER ELEMENTETS MASSE

    EMASS = dsum(NEDOF,EMV,6)

  END SUBROUTINE ADD_ANDES3


  ! ######################################################################
  ! ######################################################  ADD_ANDES4 ###
  ! ######################################################################

  SUBROUTINE ADD_ANDES4(IOP, IEL, NEDOF, &
       &                EK, EM, EMASS, diagMass, LPU, IERR)

    ! **********************************************************************
    !   Stiffness and mass matrix for ANDESQ4 element
    !
    !   Coded by Bjorn Haugen
    ! **********************************************************************

    use KindModule             , only : dp, epsDiv0_p
    use IsoMatModule           , only : isoMat2D
    use Andes4ShellModule      , only : Andes4shell_stiffmat
    use ManipMatrixModule      , only : trans3p
    use PmatModule             , only : pMatStiff
    use ReportErrorModule      , only : reportError, error_p
    use FFlLinkHandlerInterface, only : ffl_getCoor, ffl_getElmId
    use FFlLinkHandlerInterface, only : ffl_getMat, ffl_getNSM, ffl_getThick

    integer , parameter   :: neldof = 24
    integer , intent(in)  :: iop, iel, lpu
    logical , intent(in)  :: diagMass
    real(dp), intent(out) :: ek(neldof,neldof), em(neldof,neldof), emass
    integer , intent(out) :: nedof, ierr

    !! Local variables
    integer , parameter :: ltype = 2
    real(dp), parameter :: alpha = 1.5_dp
    real(dp), parameter :: beta  = 0.9_dp

    integer  :: i, j
    real(dp) :: EMV(neldof), XG(4), YG(4), ZG(4), E(3,3), TH(4)
    real(dp) :: EMOD, RNY, RHO, NSMASS
    real(dp) :: coorG(3,4), coorL(3,4), PMAT(neldof,neldof)
    real(dp) :: TelT(3,3), Tel(3,3), trel(3,4)
    real(dp) :: Xl(4), Yl(4), Zl(4), tmp(3)
    real(dp) :: Cmat(6,6), thickness, fBen, Kmat(24,24)

    !! --- Logic section ---

    nedof = neldof
    em    = 0.0_dp
    ek    = 0.0_dp

    ! --- HENTER N\DVENDIGE ELEMENTDATA

    call ffl_getcoor (XG,YG,ZG,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element nodal coordinates')
       return
    end if

    coorG(1,:) = xg
    coorG(2,:) = yg
    coorG(3,:) = zg

    !!=== Element transformation matrix

    tmp  = 0.5_dp*(coorG(:,3) + coorG(:,4))
    trel = trans3P(coorG(:,1), coorG(:,2), tmp, lpu, ierr)
    if (ierr < 0) then
       write(lpu,600) ffl_getElmId(iel)
       ierr = -ierr
600    format('  ** The 4-noded shell element',I6, &
            & ' has zero area and is ignored.' &
            / '     Check your FE model for consistency.')
       return
    end if

    TelT = trel(1:3,1:3)
    Tel  = transpose(TelT)

    do i = 1,4
       coorG(1,i) = xg(i) - xg(1)
       coorG(2,i) = yg(i) - yg(1)
       coorG(3,i) = zg(i) - zg(1)
    end do

    do i = 1,4
       coorL(:,i) = matmul(Tel,coorG(:,i))
       Xl(i) = dot_product(Tel(1,:),coorG(:,i))
       Yl(i) = dot_product(Tel(2,:),coorG(:,i))
       Zl(i) = dot_product(Tel(3,:),coorG(:,i))
    end do

    call ffl_getmat (EMOD,RNY,RHO,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element material data')
       return
    end if

    call ffl_getthick (TH,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element thickness')
       return
    end if

    thickness = (TH(1)+TH(2)+TH(3)+TH(4))/4.0_dp

    call ffl_getNSM (NSMass,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element non-structural mass')
       return
    end if

    ! Adjust material density for non-structural mass
    if (thickness > epsDiv0_p) RHO = RHO + NSMass/thickness

    if (mod(iop,2) == 0) goto 100

    call isoMat2D (EMOD,RNY,E)
    fBen = thickness*thickness*thickness / 12.0_dp
    do i = 1,3
       do j = 1,3
          Cmat(i,j)     = E(i,j) * thickness
          Cmat(i+3,j+3) = E(i,j) * fBen
       end do
    end do
    Cmat(1:3,4:6) = 0.0_dp
    Cmat(4:6,1:3) = 0.0_dp

    ! --- BEREGNER ELEMENTETS STIVHETSMATRISE

    ! Zl has to be close to zero, i.e. best fit in X-Y plane
    call Andes4shell_stiffmat (Xl,Yl,Zl,Cmat,alpha,beta,ltype,Kmat,LPU,IERR)
    IF (IERR < 0) THEN
       call ReportError (error_p,'computing shell stiffness matrix')
       return
    END IF

    ! Transform to global
    do i = 1,24,3
       do j = 1,24,3
          Kmat(i:i+2,j:j+2) = matmul(TelT,matmul(Kmat(i:i+2,j:j+2),Tel))
       end do
    end do

    ! Compute projection matrix that restores stress free rigid body motions
    CALL PMATSTIFF(coorG,PMAT,LPU,IERR)
    IF (IERR < 0) THEN
       call ReportError (error_p,'computing projection matrix')
       return
    END IF
    ! Project the stiffness matrix: K = P'*K*P
    CALL SYMTRA(Kmat(1,1),PMAT(1,1),EK(1,1),EMV(1),NEDOF,NEDOF)

    ! --- BEREGNER ELEMENTETS MASSEMATRISE

100 if (iop < 2) return
    CALL FQS35(EMV(1),XG(1),YG(1),ZG(1),TH(1),RHO)
    if (diagMass) then
       call DCOPY(NEDOF,EMV(1),1,EM(1,1),1)
    else
       call DCOPY(NEDOF,EMV(1),1,EM(1,1),NEDOF+1)
    end if

    ! --- BEREGNER ELEMENTETS MASSE (SUMMERER BIDRAG I X-RETNING)

    EMASS = dsum(NEDOF,EMV,6)

  END SUBROUTINE ADD_ANDES4


  ! ######################################################################
  ! ######################################################  ADD31      ###
  ! ######################################################################

  SUBROUTINE ADD31(IOP, IEL, NEDOF, EK, EM, EMASS, diagMass, IPSW, LPU, IERR)

    ! **********************************************************************
    !
    !   F E D E M  (FINITE ELEMENT DYNAMICS IN ELASTIC MECHANISMS)
    !
    !   P R E P R O S E S S O R :   ADD31                        S I N T E F
    !
    !
    !     DENNE SUBRUTINEN BEREGNER ELEMENTMATRISER (STIVHET/MASSE)
    !     FOR SKALLELEMENTET SCTS (ELTYPE = 31).
    !     DEN BEREGNER OGS] ELEMENTETS MASSE.
    !
    !     PROGRAMMERT AV : KETIL AAMNES
    !     DATO / VERSJON : 88-10-14 / 1.0
    !
    ! **********************************************************************

    use KindModule             , only : dp
    use ReportErrorModule      , only : reportError, error_p
    use FFlLinkHandlerInterface, only : ffl_getCoor, ffl_getMat, ffl_getThick
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_getdouble

    integer , parameter   :: neldof = 36
    integer , intent(in)  :: iop, iel, ipsw, lpu
    logical , intent(in)  :: diagMass
    real(dp), intent(out) :: ek(neldof,neldof), em(neldof,neldof), emass
    integer , intent(out) :: nedof, ierr

    !! Local variables
    integer  :: LIN, IINT, ICODIR, IMAS
    real(dp) :: EML(30,30), EMLU(30), EPS6
    real(dp) :: XG(6), YG(6), ZG(6), TH(6), LAMBI(6,3,3), EMOD, RNY, RHO

    !! --- Logic section ---

    nedof = neldof
    em    = 0.0_dp
    ek    = 0.0_dp
    call ffa_cmdlinearg_getdouble ('drillingStiff',EPS6)

    ! --- HENTER N\DVENDIGE ELEMENTDATA

    call ffl_getcoor (XG,YG,ZG,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element nodal coordinates')
       return
    end if

    call ffl_getmat (EMOD,RNY,RHO,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element material data')
       return
    end if

    call ffl_getthick (TH,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element thickness')
       return
    end if

    ! --- BEREGNER ELEMENTETS STIVHETSMATRISE

    ICODIR = 1
    IINT   = 9
    if (mod(iop,2) == 0) goto 100
    LIN    = 1
    CALL SCTS31(EML(1,1),XG(1),YG(1),ZG(1),TH(1), &
         &      IEL,EMOD,RNY,ICODIR,LAMBI(1,1,1),IINT,LIN,IPSW,LPU,IERR)
    IF (IERR < 0) THEN
       call ReportError (error_p,'Error computing thick shell stiffness matrix')
       return
    END IF

    ! Transform the rotational DOFs to global axes
    call congruenceTrans56 (EML,EK,LAMBI,6,EPS6)

    ! --- BEREGNER ELEMENTETS MASSEMATRISE

100 if (iop < 2) return
    IMAS   = 1
    CALL SCTS35(EML(1,1),EMLU(1),XG(1),YG(1),ZG(1),TH(1), &
         &      IEL,RHO,IINT,LAMBI(1,1,1),ICODIR,IMAS,IPSW,LPU,IERR)
    IF (IERR < 0) THEN
       call ReportError (error_p,'Error computing thick shell mass matrix')
       return
    END IF

    EMASS = dsum(size(emlu),EMLU,5)

    if (diagMass) then
       call expand56 (NEDOF,EMLU,EM(:,1))
    else ! Transform the rotational DOFs to global axes
       call congruenceTrans56 (EML,EM,LAMBI,6,0.0_dp)
    end if

  END SUBROUTINE ADD31


  ! ######################################################################
  ! ######################################################  ADD32      ###
  ! ######################################################################

  SUBROUTINE ADD32(IOP, IEL, NEDOF, EK, EM, EMASS, diagMass, IPSW, LPU, IERR)

    ! **********************************************************************
    !
    !   F E D E M  (FINITE ELEMENT DYNAMICS IN ELASTIC MECHANISMS)
    !
    !   P R E P R O S E S S O R :   ADD32                        S I N T E F
    !
    !
    !     DENNE SUBRUTINEN BEREGNER ELEMENTMATRISER (STIVHET/MASSE)
    !     FOR SKALLELEMENTET SCQS (ELTYPE = 32).
    !     DEN BEREGNER OGS] ELEMENTETS MASSE.
    !
    !     PROGRAMMERT AV : KETIL AAMNES
    !     DATO / VERSJON : 88-10-14 / 1.0
    !
    ! **********************************************************************

    use KindModule             , only : dp
    use ReportErrorModule      , only : reportError, error_p
    use FFlLinkHandlerInterface, only : ffl_getCoor, ffl_getMat, ffl_getThick
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_getdouble

    integer , parameter   :: neldof = 48
    integer , intent(in)  :: iop, iel, ipsw, lpu
    logical , intent(in)  :: diagMass
    real(dp), intent(out) :: ek(neldof,neldof), em(neldof,neldof), emass
    integer , intent(out) :: nedof, ierr

    !! Local variables
    integer  :: IINT, ICODIR, IMAS
    real(dp) :: EML(40,40), EMLU(40), EPS6
    real(dp) :: XG(8), YG(8), ZG(8), TH(8), LAMBI(8,3,3), EMOD, RNY, RHO

    !! --- Logic section ---

    nedof = neldof
    em    = 0.0_dp
    ek    = 0.0_dp
    call ffa_cmdlinearg_getdouble ('drillingStiff',EPS6)

    ! --- HENTER N\DVENDIGE ELEMENTDATA

    call ffl_getcoor (XG,YG,ZG,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element nodal coordinates')
       return
    end if

    call ffl_getmat (EMOD,RNY,RHO,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element material data')
       return
    end if

    call ffl_getthick (TH,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element thickness')
       return
    end if

    ! --- BEREGNER ELEMENTETS STIVHETSMATRISE

    ICODIR = 1
    IINT   = 2
    if (mod(iop,2) == 0) goto 100
    CALL SCQS31(EML(1,1),XG(1),YG(1),ZG(1),TH(1), &
         &      IEL,EMOD,RNY,ICODIR,IINT,IINT,LAMBI(1,1,1),IPSW,LPU,IERR)
    IF (IERR < 0) THEN
       call ReportError (error_p,'Error computing thick shell stiffness matrix')
       return
    END IF

    ! Transform the rotational DOFs to global axes
    call congruenceTrans56 (EML,EK,LAMBI,8,EPS6)

    ! --- BEREGNER ELEMENTETS MASSEMATRISE

100 if (iop < 2) return
    IINT   = 4
    IMAS   = 1
    CALL SCQS35(EML(1,1),EMLU(1),XG(1),YG(1),ZG(1),TH(1), &
         &      IEL,RHO,IINT,IINT,LAMBI(1,1,1),ICODIR,IMAS,IPSW,LPU,IERR)
    IF (IERR < 0) THEN
       call ReportError (error_p,'Error computing thick shell mass matrix')
       return
    END IF

    EMASS = dsum(size(EMLU),EMLU,5)

    if (diagMass) then
       call expand56 (NEDOF,EMLU,EM(:,1))
    else ! Transform the rotational DOFs to global axes
       call congruenceTrans56 (EML,EM,LAMBI,8,0.0_dp)
    end if

  END SUBROUTINE ADD32


  ! ######################################################################
  ! ######################################################  ADD41      ###
  ! ######################################################################

  SUBROUTINE ADD41(IOP, IEL, NEDOF, EK, EM, EMASS, diagMass, IPSW, LPU, IERR)

    ! **********************************************************************
    !
    !   F E D E M  (FINITE ELEMENT DYNAMICS IN ELASTIC MECHANISMS)
    !
    !   P R E P R O S E S S O R :   ADD41                        S I N T E F
    !
    !
    !     DENNE SUBRUTINEN BEREGNER ELEMENTMATRISER (STIVHET/MASSE)
    !     FOR VOLUMELEMENTET ITET (ELTYPE = 41).
    !     DEN BEREGNER OGS] ELEMENTETS MASSE.
    !
    !     PROGRAMMERT AV : KETIL AAMNES
    !     DATO / VERSJON : 88-10-14 / 1.0
    !
    ! **********************************************************************

    use KindModule             , only : dp
    use ReportErrorModule      , only : reportError, error_p
    use FFlLinkHandlerInterface, only : ffl_getCoor, ffl_getMat
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_isTrue

    integer , parameter   :: NENOD = 10
    integer , intent(in)  :: iop, iel, ipsw, lpu
    logical , intent(in)  :: diagMass
    real(dp), intent(out) :: ek(*), em(*), emass
    integer , intent(out) :: nedof, ierr

    !! Local variables
    integer  :: INT, LUMP
    real(dp) :: EMOD, RNY, RHO, EMO(495)
    real(dp) :: EMLU(3*NENOD), XG(NENOD), YG(NENOD), ZG(NENOD)

    !! --- Logic section ---

    nedof = 3*NENOD
    eMass = 0.0_dp
    call DCOPY (NEDOF*NEDOF,0.0_dp,0,EK(1),1)
    call DCOPY (NEDOF*NEDOF,0.0_dp,0,EM(1),1)

    ! --- HENTER N\DVENDIGE ELEMENTDATA

    call ffl_getcoor (XG,YG,ZG,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element nodal coordinates')
       return
    end if
    !
    call ffl_getmat (EMOD,RNY,RHO,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element material data')
       return
    end if

    ! --- BEREGNER ELEMENTETS STIVHETSMATRISE

    if (mod(iop,2) == 0) goto 100
    INT  = 5
    CALL ITET31(EMO(1),XG(1),YG(1),ZG(1),EMOD,RNY,IEL,INT,LPU,IPSW,IERR)
    IF (IERR /= 0) THEN
       call ReportError (error_p,'computing solid element stiffness matrix')
       return
    END IF

    CALL EXPSOLID(EK,EMO,NENOD)

    ! --- BEREGNER ELEMENTETS MASSEMATRISE

100 if (iop < 2) return
    INT  = 5
    LUMP = 1 ! Compute lumped mass matrix
    CALL ITET35(EMO(1),EMLU(1),XG(1),YG(1),ZG(1),RHO,INT,LUMP,IEL,LPU,IPSW,IERR)
    IF (IERR /= 0) THEN
       call ReportError (error_p,'computing solid element mass matrix')
       return
    END IF

    eMass = dsum(NEDOF,EMLU,3)

    if (diagMass) then
       call DCOPY(NEDOF,EMLU(1),1,EM(1),1)
    else if (ffa_cmdlinearg_isTrue('oldStyleMass')) then
       call DCOPY(NEDOF,EMLU(1),1,EM(1),NEDOF+1)
    else
       LUMP = 0 ! Compute consistent mass matrix
       CALL ITET35(EMO(1),EMLU(1),XG(1),YG(1),ZG(1),RHO,INT,LUMP, &
            &      IEL,LPU,IPSW,IERR)
       IF (IERR /= 0) THEN
          call ReportError (error_p,'computing solid element mass matrix')
          return
       END IF
       CALL EXPSOLID(EM,EMO,NENOD)
    end if

  END SUBROUTINE ADD41


  ! ######################################################################
  ! ######################################################  ADD42      ###
  ! ######################################################################

  SUBROUTINE ADD42(IOP, IEL, NEDOF, EK, EM, EMASS, diagMass, IPSW, LPU, IERR)

    ! **********************************************************************
    !
    !   F E D E M  (FINITE ELEMENT DYNAMICS IN ELASTIC MECHANISMS)
    !
    !   P R E P R O S E S S O R :   ADD42                        S I N T E F
    !
    !
    !     DENNE SUBRUTINEN BEREGNER ELEMENTMATRISER (STIVHET/MASSE)
    !     FOR VOLUMELEMENTET IPRI (ELTYPE = 42).
    !     DEN BEREGNER OGS] ELEMENTETS MASSE.
    !
    !     PROGRAMMERT AV : KETIL AAMNES
    !     DATO / VERSJON : 88-10-14 / 1.0
    !
    ! **********************************************************************

    use KindModule             , only : dp
    use ReportErrorModule      , only : reportError, error_p
    use FFlLinkHandlerInterface, only : ffl_getCoor, ffl_getMat

    integer , parameter   :: NENOD = 15
    integer , intent(in)  :: iop, iel, ipsw, lpu
    logical , intent(in)  :: diagMass
    real(dp), intent(out) :: ek(*), em(*), emass
    integer , intent(out) :: nedof, ierr

    !! Local variables
    integer  :: INT1, INT2, LUMP
    real(dp) :: EMOD, RNY, RHO, EMO(1080)
    real(dp) :: EMLU(3*NENOD), XG(NENOD), YG(NENOD), ZG(NENOD)

    !! --- Logic section ---

    nedof = 3*NENOD
    eMass = 0.0_dp
    call DCOPY (NEDOF*NEDOF,0.0_dp,0,EK(1),1)
    call DCOPY (NEDOF*NEDOF,0.0_dp,0,EM(1),1)

    ! --- HENTER N\DVENDIGE ELEMENTDATA

    call ffl_getcoor (XG,YG,ZG,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element nodal coordinates')
       return
    end if

    call ffl_getmat (EMOD,RNY,RHO,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element material data')
       return
    end if

    ! --- BEREGNER ELEMENTETS STIVHETSMATRISE

    if (mod(iop,2) == 0) goto 100
    INT1 = 7
    INT2 = 4
    EMO  = 0.0_dp
    CALL IPRI31(EMO(1),XG(1),YG(1),ZG(1),EMOD,RNY,IEL,INT1,INT2,LPU,IPSW,IERR)
    IF (IERR /= 0) THEN
       call ReportError (error_p,'computing solid element stiffness matrix')
       return
    END IF

    CALL EXPSOLID(EK,EMO,NENOD)

    ! --- BEREGNER ELEMENTETS MASSEMATRISE

100 if (iop < 2) return
    INT2 = 4
    LUMP = 1
    EMO  = 0.0_dp
    CALL IPRI35(EMO(1),EMLU(1),XG(1),YG(1),ZG(1),RHO,INT2,LUMP,IEL, &
         &      LPU,IPSW,IERR)
    IF (IERR /= 0) THEN
       call ReportError (error_p,'computing solid element mass matrix')
       return
    END IF

    eMass = dsum(NEDOF,EMLU,3)

    if (diagMass) then
       call DCOPY(NEDOF,EMLU(1),1,EM(1),1)
    else
       LUMP = 0
       EMO  = 0.0_dp
       CALL IPRI35(EMO(1),EMLU(1),XG(1),YG(1),ZG(1),RHO,INT2,LUMP,IEL, &
            &      LPU,IPSW,IERR)
       IF (IERR /= 0) THEN
          call ReportError (error_p,'computing solid element mass matrix')
          return
       END IF
       CALL EXPSOLID(EM,EMO,NENOD)
    end if

  END SUBROUTINE ADD42


  ! ######################################################################
  ! ######################################################  ADD43      ###
  ! ######################################################################

  SUBROUTINE ADD43(IOP, IEL, NEDOF, EK, EM, EMASS, diagMass, IPSW, LPU, IERR)

    ! **********************************************************************
    !
    !   F E D E M  (FINITE ELEMENT DYNAMICS IN ELASTIC MECHANISMS)
    !
    !   P R E P R O S E S S O R :   ADD43                        S I N T E F
    !
    !
    !     DENNE SUBRUTINEN BEREGNER ELEMENTMATRISER (STIVHET/MASSE)
    !     FOR VOLUMELEMENTET IHEX (ELTYPE = 43).
    !     DEN BEREGNER OGS] ELEMENTETS MASSE.
    !
    !     PROGRAMMERT AV : KETIL AAMNES
    !     DATO / VERSJON : 88-10-14 / 1.0
    !
    ! **********************************************************************

    use KindModule             , only : dp
    use ReportErrorModule      , only : reportError, error_p
    use FFlLinkHandlerInterface, only : ffl_getCoor, ffl_getMat

    integer , parameter   :: NENOD = 20
    integer , intent(in)  :: iop, iel, ipsw, lpu
    logical , intent(in)  :: diagMass
    real(dp), intent(out) :: ek(*), em(*), emass
    integer , intent(out) :: nedof, ierr

    !! Local variables
    integer  :: INT, LUMP
    real(dp) :: EMOD, RNY, RHO, EMO(1890)
    real(dp) :: EMLU(3*NENOD), XG(NENOD), YG(NENOD), ZG(NENOD)

    !! --- Logic section ---

    nedof = 3*NENOD
    eMass = 0.0_dp
    call DCOPY (NEDOF*NEDOF,0.0_dp,0,EK(1),1)
    call DCOPY (NEDOF*NEDOF,0.0_dp,0,EM(1),1)

    ! --- HENTER N\DVENDIGE ELEMENTDATA

    call ffl_getcoor (XG,YG,ZG,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element nodal coordinates')
       return
    end if

    call ffl_getmat (EMOD,RNY,RHO,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element material data')
       return
    end if

    ! --- BEREGNER ELEMENTETS STIVHETSMATRISE

    if (mod(iop,2) == 0) goto 100
    INT  = 3
    CALL IHEX31(EMO(1),XG(1),YG(1),ZG(1),EMOD,RNY,IEL,INT,LPU,IPSW,IERR)
    IF (IERR /= 0) THEN
       call ReportError (error_p,'computing solid element stiffness matrix')
       return
    END IF

    CALL EXPSOLID(EK,EMO,NENOD)

    ! --- BEREGNER ELEMENTETS MASSEMATRISE

100 if (iop < 2) return
    INT  = 3
    LUMP = 1
    EMO  = 0.0_dp
    CALL IHEX35(EMO(1),EMLU(1),XG(1),YG(1),ZG(1),RHO,INT,INT,INT,LUMP, &
         &      IEL,LPU,IPSW,IERR)
    IF (IERR /= 0) THEN
       call ReportError (error_p,'computing solid element mass matrix')
       return
    END IF

    eMass = dsum(NEDOF,EMLU,3)

    if (diagMass) then
       call DCOPY(NEDOF,EMLU(1),1,EM(1),1)
    else
       INT  = 3
       LUMP = 0
       EMO  = 0.0_dp
       CALL IHEX35(EMO(1),EMLU(1),XG(1),YG(1),ZG(1),RHO,INT,INT,INT,LUMP, &
            &      IEL,LPU,IPSW,IERR)
       IF (IERR /= 0) THEN
          call ReportError (error_p,'computing solid element mass matrix')
          return
       END IF
       CALL EXPSOLID(EM,EMO,NENOD)
    end if

  END SUBROUTINE ADD43


  ! ######################################################################
  ! ######################################################  ADD44      ###
  ! ######################################################################

  SUBROUTINE ADD44(IOP, IEL, NEDOF, EK, EM, EMASS, diagMass, IPSW, LPU, IERR)

    ! **********************************************************************
    !
    !   F E D E M  (FINITE ELEMENT DYNAMICS IN ELASTIC MECHANISMS)
    !
    !   P R E P R O S E S S O R :   ADD44                        S I N T E F
    !
    !
    !     DENNE SUBRUTINEN BEREGNER ELEMENTMATRISER (STIVHET/MASSE)
    !     FOR VOLUMELEMENTET HEXA (ELTYPE = 44).
    !     DEN BEREGNER OGS] ELEMENTETS MASSE.
    !
    !     PROGRAMMERT AV : KETIL AAMNES
    !     DATO / VERSJON : 88-10-14 / 1.0
    !
    ! **********************************************************************

    use KindModule             , only : dp
    use IsoMatModule           , only : isoMat3D
    use ReportErrorModule      , only : reportError, error_p
    use FFlLinkHandlerInterface, only : ffl_getCoor, ffl_getMat
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_isTrue

    integer , parameter   :: NENOD = 8, NIP = 3
    integer , intent(in)  :: iop, iel, ipsw, lpu
    logical , intent(in)  :: diagMass
    real(dp), intent(out) :: ek(*), em(*), emass
    integer , intent(out) :: nedof, ierr

    !! Local variables
    integer  :: JOP
    real(dp) :: A1(3*NENOD,9), A2(9,9), EMLU(3*NENOD), EMTR(6,6)
    real(dp) :: XG(NENOD), YG(NENOD), ZG(NENOD), EMOD, RNY, RHO

    !! --- Logic section ---

    nedof = 3*NENOD
    eMass = 0.0_dp
    call DCOPY (NEDOF*NEDOF,0.0_dp,0,EK(1),1)
    call DCOPY (NEDOF*NEDOF,0.0_dp,0,EM(1),1)

    ! --- HENTER N\DVENDIGE ELEMENTDATA

    call ffl_getcoor (XG,YG,ZG,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element nodal coordinates')
       return
    end if

    call ffl_getmat (EMOD,RNY,RHO,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element material data')
       return
    end if

    if (mod(iop,2) == 0) goto 100
    call isoMat3D (EMOD,RNY,EMTR)

    ! --- BEREGNER ELEMENTETS STIVHETSMATRISE

    if (ffa_cmdlinearg_isTrue('useIncompatibleModes')) then
       JOP = 1
    else
       JOP = 0
    end if
    CALL HEXA31(EK(1),A1(1,1),A2(1,1),EMTR(1,1),XG(1),YG(1),ZG(1), &
         &      NIP,JOP,IEL,LPU,IPSW,IERR)
    IF (IERR /= 0) THEN
       call ReportError (error_p,'computing solid element stiffness matrix')
       return
    END IF

    ! --- BEREGNER ELEMENTETS MASSEMATRISE

100 if (iop < 2) return
    CALL HEXA35(EM(1),XG(1),YG(1),ZG(1),RHO,LPU,IERR)
    IF (IERR /= 0) THEN
       call ReportError (error_p,'computing solid element mass matrix')
       return
    END IF

    if (diagMass) then
       eMass = lumpSolidMass(NEDOF,EM,EM)
    else
       eMass = lumpSolidMass(NEDOF,EM,EMLU)
    end if

  END SUBROUTINE ADD44


  ! ######################################################################
  ! ######################################################  EXPSOLID   ###
  ! ######################################################################

  SUBROUTINE EXPSOLID (EKF,EKO,NENOD)

    ! **********************************************************************
    !
    !   F E D E M  (FINITE ELEMENT DYNAMICS IN ELASTIC MECHANISMS)
    !
    !   P R E P R O S E S S O R :   EXPSOLID                     S I N T E F
    !
    !
    !     DENNE SUBRUTINEN BEREGNER DEN FULLE ELEMENTMATRISEN (STIVHET/
    !     MASSE) FOR ET VOLUMELEMENT MED NENOD ANTALL KNUTEPUNKTER. DEN
    !     OPPRINNELIGE ELEMENTMATRISEN ER LAGRET SOM ET ENDIMENSJONALT
    !     ARRAY.
    !
    !     PROGRAMMERT AV : KETIL AAMNES
    !     DATO / VERSJON : 88-10-14 / 1.0
    !
    ! **********************************************************************

    use KindModule, only : dp

    integer , intent(in)  :: NENOD
    real(dp), intent(in)  :: EKO(:)
    real(dp), intent(out) :: EKF(NENOD*3,NENOD*3)

    !! Local variables
    integer :: I, J, K

    !! --- Logic section ---

    EKF = 0.0_dp

    ! --- OPPBYGGING AV \VRE TRIANGEL AV ELEMENTSTIVHETSMATRISEN EKF

    K = 1
    DO I = 1, NENOD*3-2, 3
       DO J = I, NENOD*3-2, 3
          EKF(I  ,J  ) = EKO(K  )
          EKF(I+1,J  ) = EKO(K+1)
          EKF(I+2,J  ) = EKO(K+2)
          EKF(I  ,J+1) = EKO(K+3)
          EKF(I+1,J+1) = EKO(K+4)
          EKF(I+2,J+1) = EKO(K+5)
          EKF(I  ,J+2) = EKO(K+6)
          EKF(I+1,J+2) = EKO(K+7)
          EKF(I+2,J+2) = EKO(K+8)
          K = K + 9
       END DO
    END DO

    ! --- OPPBYGGING AV NEDRE TRIANGEL - SYMMETRI

    DO I = 1, NENOD*3
       DO J = 1, NENOD*3
          EKF(J,I) = EKF(I,J)
       END DO
    END DO

  END SUBROUTINE EXPSOLID


  ! ######################################################################
  ! ######################################################  ADD45      ###
  ! ######################################################################

  SUBROUTINE ADD45(IOP, IEL, NEDOF, EK, EM, EMASS, diagMass, IPSW, LPU, IERR)

    ! **********************************************************************
    !
    !   F E D E M  (FINITE ELEMENT DYNAMICS IN ELASTIC MECHANISMS)
    !
    !   P R E P R O S E S S O R :   ADD45
    !
    !     The routine computes element stiffness- and mass matrices for
    !     the volume element CSTET (Constant Strain Tetrahedron,
    !     ELTYPE = 45). The routine also computes the element mass.
    !
    !     PROGRAMMERT AV : Bjorn Haugen
    !     DATO / VERSJON : 97-9-4 / 1.0
    !
    ! **********************************************************************

    use KindModule             , only : dp
    use IsoMatModule           , only : isoMat3D
    use CSTetraModule          , only : CSTETSTIFF, CSTETMASS
    use ReportErrorModule      , only : reportError, error_p
    use FFlLinkHandlerInterface, only : ffl_getCoor, ffl_getMat

    integer , intent(in)  :: iop, iel, ipsw, lpu
    logical , intent(in)  :: diagMass
    real(dp), intent(out) :: ek(*), em(*), emass
    integer , intent(out) :: nedof, ierr

    !! Local variables
    real(dp) :: EMLU(12), EMTR(6,6), XG(4),YG(4),ZG(4), EMOD,RNY,RHO,VOLUME

    !! --- Logic section ---

    nedof = 12
    eMass = 0.0_dp
    call DCOPY (NEDOF*NEDOF,0.0_dp,0,EK(1),1)
    call DCOPY (NEDOF*NEDOF,0.0_dp,0,EM(1),1)

    ! --- Get necessary element data

    call ffl_getcoor (XG,YG,ZG,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element nodal coordinates')
       return
    end if

    call ffl_getmat (EMOD,RNY,RHO,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element material data')
       return
    end if

    if (mod(iop,2) == 0) goto 100
    call isoMat3D (EMOD,RNY,EMTR)

    ! --- Calculate stiffness

    CALL CSTETSTIFF(XG,YG,ZG,EMTR,EK(1:nedof*nedof),VOLUME,LPU,IPSW,IERR)
    IF (IERR /= 0) THEN
       call ReportError (error_p,'computing solid element stiffness matrix')
       return
    END IF

    ! --- Calculate mass

100 if (iop < 2) return
    CALL CSTETMASS(XG,YG,ZG,RHO,EM(1:nedof*nedof),VOLUME,LPU,IPSW,IERR)
    IF (IERR /= 0) THEN
       call ReportError (error_p,'computing solid element mass matrix')
       return
    END IF

    if (diagMass) then
       eMass = lumpSolidMass(NEDOF,EM,EM)
    else
       eMass = lumpSolidMass(NEDOF,EM,EMLU)
    end if

  END SUBROUTINE ADD45


  ! ######################################################################
  ! ######################################################  ADD46      ###
  ! ######################################################################

  SUBROUTINE ADD46(IOP, IEL, NEDOF, EK, EM, EMASS, diagMass, IPSW, LPU, IERR)

    ! **********************************************************************
    !
    !   F E D E M  (FINITE ELEMENT DYNAMICS IN ELASTIC MECHANISMS)
    !
    !   P R E P R O S E S S O R :   ADD46
    !
    !     The routine computes the element stiffness- and mass matrices for
    !     the volume element IPRI6 (Isoparametric PRIsm element with 6
    !     nodes, ELTYPE = 46). The routine also computes the element mass.
    !
    !     PROGRAMMERT AV : Trond Arne Svidal
    !     DATO / VERSJON : 99-11-3 / 1.0
    !
    ! **********************************************************************

    use KindModule             , only : dp
    use IsoMatModule           , only : isoMat3D
    use Ipri6Module            , only : IPRI6STIFF, IPRI6MASS
    use ReportErrorModule      , only : reportError, error_p
    use FFlLinkHandlerInterface, only : ffl_getCoor, ffl_getMat

    integer , intent(in)  :: iop, iel, ipsw, lpu
    logical , intent(in)  :: diagMass
    real(dp), intent(out) :: ek(*), em(*), emass
    integer , intent(out) :: nedof, ierr

    !! Local variables
    real(dp) :: EMLU(18), EMTR(6,6), XG(6), YG(6), ZG(6), EMOD, RNY, RHO

    !! --- Logic section ---

    nedof = 18
    eMass = 0.0_dp
    call DCOPY (NEDOF*NEDOF,0.0_dp,0,ek(1),1)
    call DCOPY (NEDOF*NEDOF,0.0_dp,0,em(1),1)

    ! --- Get necessary element data

    call ffl_getcoor (XG,YG,ZG,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element nodal coordinates')
       return
    end if

    call ffl_getmat (EMOD,RNY,RHO,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element material data')
       return
    end if

    if (mod(iop,2) == 0) goto 100
    call isoMat3D (EMOD,RNY,EMTR)

    ! --- Calculate stiffness

    CALL IPRI6STIFF(XG, YG, ZG, EMTR, EK(1:nedof*nedof), IEL, LPU, IPSW, IERR)
    IF (IERR /= 0) THEN
       call ReportError (error_p,'computing solid element stiffness matrix')
       return
    END IF

    ! --- Calulate mass

100 if (iop < 2) return
    CALL IPRI6MASS(XG, YG, ZG, RHO, EM(1:nedof*nedof), IEL, LPU, IPSW, IERR)
    IF (IERR /= 0) THEN
       call ReportError (error_p,'computing solid element mass matrix')
       return
    END IF

    if (diagMass) then
       eMass = lumpSolidMass(NEDOF,EM,EM)
    else
       eMass = lumpSolidMass(NEDOF,EM,EMLU)
    end if

  END SUBROUTINE ADD46


  ! ######################################################################
  ! ######################################################  ADD50      ###
  ! ######################################################################

  SUBROUTINE ADD50(sam, IEL, NEDOF, EM, Mt, Mr, diagMass, IPSW, LPU, IERR)

    ! **********************************************************************
    !
    !   F E D E M  (FINITE ELEMENT DYNAMICS IN ELASTIC MECHANISMS)
    !
    !   P R E P R O S E S S O R :   ADD50
    !
    !     Concentrated mass at a grid point (ELTYPE = 50).
    !
    !     PROGRAMMERT AV : Knut Morten Okstad
    !     DATO / VERSJON : 2004-03-23 / 1.0
    !
    ! **********************************************************************

    use SamModule              , only : SamType, dp
    use ManipMatrixModule      , only : writeObject
    use FFlLinkHandlerInterface, only : ffl_getElmId

    type(SamType), intent(in)  :: sam
    integer      , intent(in)  :: iel, ipsw, lpu
    logical      , intent(in)  :: diagMass
    real(dp)     , intent(in)  :: Mt, Mr
    real(dp)     , intent(out) :: em(*)
    integer      , intent(out) :: nedof, ierr

    !! Local variables
    integer :: i, mnode

    !! --- Logic section ---

    ierr  = 0
    mnode = sam%mmnpc(sam%mpmnpc(iel))
    nedof = min(6,sam%madof(mnode+1) - sam%madof(mnode))

    if (ipsw > -7) write(lpu,991) ffl_getElmId(iel),Mt,Mr

    if (diagMass) then
       em(1:3) = Mt
       em(4:6) = Mr
    else
       em(1:nedof*nedof) = 0.0_dp
       do i = 1, 3
          em(i+nedof*(i-1)) = Mt
       end do
       do i = 4, nedof
          em(i+nedof*(i-1)) = Mr
       end do
    end if

    if (ipsw > -5) then
       write(lpu,994) ffl_getElmId(iel)
       if (diagMass) then
          call writeObject (em(1:nedof),lpu)
       else
          call writeObject (reshape(em(1:nedof*nedof),(/nedof,nedof/)),lpu)
       end if
    end if

991 format(/5X,'Mass element',I10/8X,'M  =',1P,2E13.5)
994 format(//5X,'Mass matrix for mass element',I10/)

  END SUBROUTINE ADD50


  ! ######################################################################
  ! ######################################################  ADD51      ###
  ! ######################################################################

  SUBROUTINE ADD51(sam, IEL, NEDOF, EM, EMASS, diagMass, IPSW, LPU, IERR)

    ! **********************************************************************
    !
    !   F E D E M  (FINITE ELEMENT DYNAMICS IN ELASTIC MECHANISMS)
    !
    !   P R E P R O S E S S O R :   ADD51
    !
    !     Concentrated mass at a grid point (ELTYPE = 51).
    !
    !     PROGRAMMERT AV : Knut Morten Okstad
    !     DATO / VERSJON : 2000-11-13 / 1.0
    !
    ! **********************************************************************

    use SamModule              , only : SamType, dp
    use ManipMatrixModule      , only : writeObject
    use ReportErrorModule      , only : reportError, error_p, warning_p
    use FFlLinkHandlerInterface, only : ffl_getMass, ffl_getElmId

    type(SamType), intent(in)  :: sam
    integer      , intent(in)  :: iel, ipsw, lpu
    logical      , intent(in)  :: diagMass
    real(dp)     , intent(out) :: em(*), emass
    integer      , intent(out) :: nedof, ierr

    !! Local variables
    integer            :: i, j, k, mnode
    real(dp)           :: M(36), W(6), Work(36)
    character(len=128) :: errMsg

    !! --- Logic section ---

    eMass = 0.0_dp

    call ffl_getMass (em,iel,nedof,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving concentrated node mass')
       return
    end if

    mnode = sam%mmnpc(sam%mpmnpc(iel))
    nedof = min(nedof,sam%madof(mnode+1) - sam%madof(mnode))

    if (nedof < 6) then

       ! Compress the 6x6 mass matrix down to <nedof>x<nedof>
       k = 0
       do j = 1, 6
          do i = 1, 6
             if (i <= nedof .and. j <= nedof) then
                k = k + 1
                em(k) = em(i+6*j-6)
             else
                eMass = eMass + abs(em(i+6*j-6))
             end if
          end do
       end do

       if (eMass > 0.0_dp) then
          write(errMsg,600) sam%minex(mnode), ffl_getElmId(iel), nedof
          call ReportError (warning_p,errMsg,'Mass contributions to the '// &
               &            'non-existing DOFs (inertia terms) are ignored.')
       end if

    end if

    eMass = (em(1)+em(nedof+2)+em(2*nedof+3))/3.0_dp

    ! Check that the matrix is positive definite by computing its eigenvalues
    call DCOPY (nedof*nedof,em(1),1,M(1),1)
    call DSYEV ('N','U',nedof,M(1),nedof,w(1),work(1),size(work),ierr)
    if (ierr /= 0) then
       write(errMsg,610) ffl_getElmId(iel)
       call ReportError (warning_p,'Eigenvalue solver failed to converge', &
            &            errMsg,addString='DSYEV call in ADD51')
       if (ierr > 0) write(lpu,611) ierr
       ierr = 0
    else if (w(1) <= 0.0_dp) then
       write(errMsg,620) ffl_getElmId(iel)
       call ReportError (warning_p,errMsg)
       write(lpu,621) w(1:nedof)
    end if

    if (diagMass) then
       do i = 1, nedof
          em(i) = sum(em(nedof*(i-1)+1:nedof*i))
       end do
    end if

    if (ipsw > -5) then
       write(lpu,994) ffl_getElmId(iel)
       if (diagMass) then
          call writeObject (em(1:nedof),lpu)
       else
          call writeObject (reshape(em(1:nedof),(/nedof,nedof/)),lpu)
       end if
    end if

600 format('Node',I10,'  which is connected to mass element',I10, &
         & '  has only',I2,' DOFs.')
610 format('The element matrix of mass element',I10,' is invalid')
611 format(I11,' off-diagonal term(s) of an intermediate tridiagonal form', &
         &/10X,'did not converge to zero')
620 format('The element matrix of mass element',I10,' is not positive definite')
621 format(10X,'Eigenvalues :',1P6E13.5)
994 format(//5X,'Mass matrix for mass element',I10/)

  END SUBROUTINE ADD51


  ! ######################################################################
  ! ######################################################  ADD70      ###
  ! ######################################################################

  SUBROUTINE ADD70(IEL, NEDOF, EK, Kt, Kr, IPSW, LPU, IERR)

    ! **********************************************************************
    !
    !   F E D E M  (FINITE ELEMENT DYNAMICS IN ELASTIC MECHANISMS)
    !
    !   P R E P R O S E S S O R :   ADD70
    !
    !
    !     Two-noded spring element with eccentricities (ELTYPE = 70).
    !
    !     PROGRAMMERT AV : Knut Morten Okstad
    !     DATO / VERSJON : 2004-03-23 / 1.0
    !
    ! **********************************************************************

    use KindModule             , only : dp
    use ManipMatrixModule      , only : writeObject
    use ReportErrorModule      , only : reportError, error_p
    use FFlLinkHandlerInterface, only : ffl_getCoor, ffl_getElmId

    integer , parameter   :: neldof = 12
    integer , intent(in)  :: iel, ipsw, lpu
    real(dp), intent(in)  :: Kt, Kr
    real(dp), intent(out) :: ek(neldof,neldof)
    integer , intent(out) :: nedof, ierr

    !! Local variables
    integer  :: i
    real(dp) :: XG(3), YG(3), ZG(3), S(6)

    !! --- Logic section ---

    nedof = neldof
    ek    = 0.0_dp

    ! Get necessary element data

    call ffl_getCoor (XG,YG,ZG,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element nodal coordinates')
       return
    end if

    if (ipsw > -7) write(lpu,991) ffl_getElmId(IEL),Kt,Kr

    ! Set up the local stiffness matrix for the spring element

    do i = 1, 3
       EK(  i,  i) =  Kt
       EK(6+i,6+i) =  Kt
       EK(  i,6+i) = -Kt
       EK(6+i,  i) = -Kt
    end do
    do i = 4, 6
       EK(  i,  i) =  Kr
       EK(6+i,6+i) =  Kr
       EK(  i,6+i) = -Kr
       EK(6+i,  i) = -Kr
    end do

    ! Transform from actual spring position to the nodal points

    S(1) = XG(2) - XG(1)
    S(2) = YG(2) - YG(1)
    S(3) = ZG(2) - ZG(1)
    S(4) = XG(3) - XG(1)
    S(5) = YG(3) - YG(1)
    S(6) = ZG(3) - ZG(1)
    if (ipsw > -7) write(lpu,993) S
    call TRIX30 (EK,S(1),S(2),S(3),S(4),S(5),S(6))

    if (ipsw > -5) then
       write(LPU,994) ffl_getElmId(IEL)
       call writeObject (EK,LPU,nellIn=nedof,eps=1.0e-4_dp)
    end if

991 format(/5X,'Spring element',I10/8X,'K  =',1P,2E13.5)
993 format(8X,'E1 =',1P,3E13.5 / 8X,'E2 =',1P,3E13.5)
994 format(//5X,'Stiffness matrix for spring element',I10/)

  END SUBROUTINE ADD70


  ! ######################################################################
  ! ######################################################  ADD72      ###
  ! ######################################################################

  SUBROUTINE ADD72(IEL, NEDOF, EK, IPSW, LPU, IERR)

    ! **********************************************************************
    !
    !   F E D E M  (FINITE ELEMENT DYNAMICS IN ELASTIC MECHANISMS)
    !
    !   P R E P R O S E S S O R :   ADD72
    !
    !
    !     Two-noded spring element with eccentricities (ELTYPE = 72).
    !
    !     PROGRAMMERT AV : Knut Morten Okstad
    !     DATO / VERSJON : 2003-09-01 / 1.0
    !
    ! **********************************************************************

    use KindModule             , only : dp
    use ManipMatrixModule      , only : writeObject
    use ReportErrorModule      , only : reportError, error_p
    use FFlLinkHandlerInterface, only : ffl_getCoor, ffl_getElmId
    use FFlLinkHandlerInterface, only : ffl_getElCoorSys, ffl_getBush

    integer , parameter   :: neldof = 12
    integer , intent(in)  :: iel, ipsw, lpu
    real(dp), intent(out) :: ek(neldof,neldof)
    integer , intent(out) :: nedof, ierr

    !! Local variables
    integer  :: i
    real(dp) :: XG(3), YG(3), ZG(3), S(6), T(3,3)

    !! --- Logic section ---

    nedof = neldof
    ek    = 0.0_dp

    ! Get necessary element data

    call ffl_getCoor (XG,YG,ZG,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element nodal coordinates')
       return
    end if

    call ffl_getElCoorSys (T,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element coordinate system')
       return
    end if

    call ffl_getBush (S,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving bushing element data')
       return
    end if
    if (ipsw > -7) write(lpu,991) ffl_getElmId(IEL),S

    ! Set up the local stiffness matrix for the spring element

    do i = 1, 6
       EK(  i,  i) =  S(i)
       EK(6+i,6+i) =  S(i)
       EK(  i,6+i) = -S(i)
       EK(6+i,  i) = -S(i)
    end do

    ! Transform to global coordinates

    if (ipsw > -7) write(lpu,992) (T(i,:),i=1,3)
    call MPRO30 (EK,T,NEDOF)

    ! Transform from actual spring position to the nodal points

    S(1) = XG(2) - XG(1)
    S(2) = YG(2) - YG(1)
    S(3) = ZG(2) - ZG(1)
    S(4) = XG(3) - XG(1)
    S(5) = YG(3) - YG(1)
    S(6) = ZG(3) - ZG(1)
    if (ipsw > -7) write(lpu,993) S
    call TRIX30 (EK,S(1),S(2),S(3),S(4),S(5),S(6))

    if (ipsw > -5) then
       write(LPU,994) ffl_getElmId(IEL)
       call writeObject (EK,LPU,nellIn=nedof,eps=1.0e-4_dp)
    end if

991 format(/5X,'Bushing element',I10/8X,'K  =',1P,6E13.5)
992 format(8X,'T  =',1P,3E13.5 / (12X,1P,3E13.5) )
993 format(8X,'E1 =',1P,3E13.5 / 8X,'E2 =',1P,3E13.5)
994 format(//5X,'Stiffness matrix for bushing element',I10/)

  END SUBROUTINE ADD72


  subroutine LOAD21 (IEL, P, NEDOF, ES, IPSW, LPU, IERR)

    !!==========================================================================
    !! This routine computes the element load vector for the thin shell element
    !! FTS (3-noded Triangle, ELTYPE = 21), due to a pressure load.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 4 Apr 2008 / 1.0
    !!==========================================================================

    use KindModule             , only : dp
    use ManipMatrixModule      , only : writeObject
    use ReportErrorModule      , only : reportError, error_p, debugFileOnly_p
    use FFlLinkHandlerInterface, only : ffl_getCoor, ffl_getElmId

    integer , parameter   :: NENOD = 3
    integer , intent(in)  :: iel, ipsw, lpu
    real(dp), intent(in)  :: p(3,NENOD)
    real(dp), intent(out) :: es(6*NENOD)
    integer , intent(out) :: nedof, ierr

    !! Local variables
    real(dp) :: XG(NENOD), YG(NENOD), ZG(NENOD)

    !! --- Logic section ---

    nedof = 6*NENOD

    !! Get global nodal coordinates for the element
    call ffl_getcoor (XG,YG,ZG,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element nodal coordinates')
       call ReportError (debugFileOnly_p,'LOAD21')
    end if

    !! Calculate the element load vector
    call FTS33 (ES,XG,YG,ZG,P(1,:),P(2,:),P(3,:),1)

    if (ipsw >= 2) then
       write(LPU,900) ffl_getElmId(IEL)
       call writeObject (ES,LPU,nellIn=6,eps=1.0e-4_dp)
       call writeObject (P,LPU,'Pressure intensities')
    end if

900 format(//5X,'Load vector for shell element',I10/)

  end subroutine LOAD21


  subroutine LOAD22 (IEL, P, NEDOF, ES, IPSW, LPU, IERR)

    !!==========================================================================
    !! This routine computes the element load vector for the thin shell element
    !! FQS (4-noded Quadrilateral, ELTYPE = 22), due to a pressure load.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 4 Apr 2008 / 1.0
    !!==========================================================================

    use KindModule             , only : dp
    use ManipMatrixModule      , only : writeObject
    use ReportErrorModule      , only : reportError, error_p, debugFileOnly_p
    use FFlLinkHandlerInterface, only : ffl_getCoor, ffl_getElmId

    integer , parameter   :: NENOD = 4
    integer , intent(in)  :: iel, ipsw, lpu
    real(dp), intent(in)  :: p(3,NENOD)
    real(dp), intent(out) :: es(6*NENOD)
    integer , intent(out) :: nedof, ierr

    !! Local variables
    real(dp) :: XG(NENOD), YG(NENOD), ZG(NENOD)

    !! --- Logic section ---

    nedof = 6*NENOD

    !! Get global nodal coordinates for the element
    call ffl_getcoor (XG,YG,ZG,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element nodal coordinates')
       call ReportError (debugFileOnly_p,'LOAD22')
    end if

    !! Calculate the element load vector
    call FQS33 (ES,XG,YG,ZG,P(1,:),P(2,:),P(3,:),1)

    if (ipsw >= 2) then
       write(LPU,900) ffl_getElmId(IEL)
       call writeObject (ES,LPU,nellIn=6,eps=1.0e-4_dp)
       call writeObject (P,LPU,'Pressure intensities')
    end if

900 format(//5X,'Load vector for shell element',I10/)

  end subroutine LOAD22


  subroutine LOAD41 (IEL, IFACE, P, NEDOF, ES, IPSW, LPU, IERR)

    !!==========================================================================
    !! This routine computes the element load vector for the volume element
    !! ITET (10-noded Isoparametric Tetrahedron, ELTYPE = 41), due to
    !! a pressure load on one of its faces.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 3 Apr 2008 / 1.0
    !!==========================================================================

    use KindModule             , only : dp
    use ReportErrorModule      , only : reportError, error_p, debugFileOnly_p
    use FFlLinkHandlerInterface, only : ffl_getCoor

    integer , parameter   :: NENOD = 10
    integer , intent(in)  :: iel, iface, ipsw, lpu
    real(dp), intent(in)  :: p(3,6)
    real(dp), intent(out) :: es(3*NENOD)
    integer , intent(out) :: nedof, ierr

    !! Local variables
    integer  :: i, ISURF
    real(dp) :: XG(NENOD), YG(NENOD), ZG(NENOD), PG(3,NENOD)

    !! Face topology, this must be identical to similar array in FFlTET10::init
    integer, parameter :: m_face(24) = (/ 1, 6, 5, 4, 3, 2, &
         &                                3, 4, 5, 9,10, 8, &
         &                                1, 7,10, 9, 5, 6, &
         &                                1, 2, 3, 8,10, 7 /)
    integer, parameter :: mface(6,4) = reshape(m_face,(/6,4/))

    !! --- Logic section ---

    nedof = 3*NENOD

    if (iface < 1 .or. iface > size(mface,2)) then
       ierr = -1
       call ReportError (error_p,'Invalid face index', ierr=iface)
       goto 990
    end if

    !! Get global nodal coordinates for the element
    call ffl_getcoor (XG,YG,ZG,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element nodal coordinates')
       goto 990
    end if

    !! Expand the face pressure intensity array to a 3*NENOD array
    PG = 0.0_dp
    do i = 1, size(mface,1)
       PG(:,mface(i,iface)) = p(:,i)
    end do

    !! Calculate the element load vector
    ISURF = 10**(4-iface)
    call ITET33 (ES,XG,YG,ZG,PG,ISURF,2,0.0_dp,1,0,4,IEL,LPU,IPSW,IERR)
    if (ierr == 0) return

990 call ReportError (debugFileOnly_p,'LOAD41')

  end subroutine LOAD41


  subroutine LOAD42 (IEL, IFACE, P, NEDOF, ES, IPSW, LPU, IERR)

    !!==========================================================================
    !! This routine computes the element load vector for the volume element
    !! IPRI (15-noded Isoparametric Prism, ELTYPE = 42), due to
    !! a pressure load on one of its faces.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 3 Apr 2008 / 1.0
    !!==========================================================================

    use KindModule             , only : dp
    use ReportErrorModule      , only : reportError, error_p, debugFileOnly_p
    use FFlLinkHandlerInterface, only : ffl_getCoor

    integer , parameter   :: NENOD = 15
    integer , intent(in)  :: iel, iface, ipsw, lpu
    real(dp), intent(in)  :: p(3,8)
    real(dp), intent(out) :: es(3*NENOD)
    integer , intent(out) :: nedof, ierr

    !! Local variables
    integer  :: i, ISURF
    real(dp) :: XG(NENOD), YG(NENOD), ZG(NENOD), PG(3,NENOD)

    !! Face topology, this must be identical to similar array in FFlWEDG15::init
    integer, parameter :: m_face(40) = (/ 1, 2, 3, 8,12,11,10, 7, &
         &                                3, 4, 5, 9,14,13,12, 8, &
         &                                1, 7,10,15,14, 9, 5, 6, &
         &                                1, 6, 5, 4, 3, 2, 0, 0, &
         &                               10,11,12,13,14,15, 0, 0 /)
    integer, parameter :: mface(8,5) = reshape(m_face,(/8,5/))

    !! --- Logic section ---

    nedof = 3*NENOD

    if (iface < 1 .or. iface > size(mface,2)) then
       ierr = -1
       call ReportError (error_p,'Invalid face index', ierr=iface)
       goto 990
    end if

    !! Get global nodal coordinates for the element
    call ffl_getcoor (XG,YG,ZG,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element nodal coordinates')
       goto 990
    end if

    !! Expand the face pressure intensity array to a 3*NENOD array
    PG = 0.0_dp
    do i = 1, size(mface,1)
       if (mface(i,iface) > 0) PG(:,mface(i,iface)) = p(:,i)
    end do

    !! Calculate the element load vector
    ISURF = 10**(5-iface)
    call IPRI33 (ES,XG,YG,ZG,PG,ISURF,2,0.0_dp,1,2,2,4,IEL,LPU,IPSW,IERR)
    if (ierr == 0) return

990 call ReportError (debugFileOnly_p,'LOAD42')

  end subroutine LOAD42


  subroutine LOAD43 (IEL, IFACE, P, NEDOF, ES, IPSW, LPU, IERR)

    !!==========================================================================
    !! This routine computes the element load vector for the volume element
    !! IHEX (20-noded Isoparametric Hexahedron, ELTYPE = 43), due to
    !! a pressure load on one of its faces.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 3 Apr 2008 / 1.0
    !!==========================================================================

    use KindModule             , only : dp
    use ReportErrorModule      , only : reportError, error_p, debugFileOnly_p
    use FFlLinkHandlerInterface, only : ffl_getCoor

    integer , parameter   :: NENOD = 20
    integer , intent(in)  :: iel, iface, ipsw, lpu
    real(dp), intent(in)  :: p(3,8)
    real(dp), intent(out) :: es(3*NENOD)
    integer , intent(out) :: nedof, ierr

    !! Local variables
    integer  :: i, ISURF
    real(dp) :: XG(NENOD), YG(NENOD), ZG(NENOD), PG(3,NENOD)

    !! Face topology, this must be identical to similar array in FFlHEX20::init
    integer, parameter :: m_face(48) = (/ 1, 2, 3,10,15,14,13, 9, &
         &                                 3, 4, 5,11,17,16,15,10, &
         &                                 5, 6, 7,12,19,18,17,11, &
         &                                 1, 9,13,20,19,12, 7, 8, &
         &                                13,14,15,16,17,18,19,20, &
         &                                 1, 8, 7, 6, 5, 4, 3, 2 /)
    integer, parameter :: mface(8,6) = reshape(m_face,(/8,6/))

    !! --- Logic section ---

    nedof = 3*NENOD

    if (iface < 1 .or. iface > size(mface,2)) then
       ierr = -1
       call ReportError (error_p,'Invalid face index', ierr=iface)
       goto 990
    end if

    !! Get global nodal coordinates for the element
    call ffl_getcoor (XG,YG,ZG,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element nodal coordinates')
       goto 990
    end if

    !! Expand the face pressure intensity array to a 3*NENOD array
    PG = 0.0_dp
    do i = 1, size(mface,1)
       if (mface(i,iface) > 0) PG(:,mface(i,iface)) = p(:,i)
    end do

    !! Calculate the element load vector
    ISURF = 10**(6-iface)
    call IHEX33 (ES,XG,YG,ZG,PG,ISURF,2,0.0_dp,1,2,2,2,IEL,LPU,IPSW,IERR)
    if (ierr == 0) return

990 call ReportError (debugFileOnly_p,'LOAD43')

  end subroutine LOAD43


  subroutine LOAD44 (IEL, IFACE, P, NEDOF, ES, IPSW, LPU, IERR)

    !!==========================================================================
    !! This routine computes the element load vector for the volume element
    !! HEXA (8-noded trilinear Hexahedron, ELTYPE = 44), due to
    !! a (constant) pressure load on one of its faces.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 3 Apr 2008 / 1.0
    !!==========================================================================

    use KindModule             , only : dp
    use IsoMatModule           , only : isoMat3D
    use ReportErrorModule      , only : reportError, error_p, debugFileOnly_p
    use FFlLinkHandlerInterface, only : ffl_getCoor, ffl_getmat
    use FFaCmdLineArgInterface , only : ffa_cmdlinearg_isTrue

    integer , parameter   :: NENOD = 8, NIP = 3
    integer , intent(in)  :: iel, iface, ipsw, lpu
    real(dp), intent(in)  :: p(3,4)
    real(dp), intent(out) :: es(3*NENOD)
    integer , intent(out) :: nedof, ierr

    !! Local variables
    integer  :: IOP
    real(dp) :: XG(NENOD), YG(NENOD), ZG(NENOD), PX, PY, PZ
    real(dp) :: EMOD, RNY, E(36), EK(24,24), A1(24,9), A2(9,9)

    !! --- Logic section ---

    nedof = 3*NENOD

    !! Get global nodal coordinates for the element
    call ffl_getcoor (XG,YG,ZG,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element nodal coordinates')
       goto 990
    end if

    if (ffa_cmdlinearg_isTrue('useIncompatibleModes')) then
       !! Incompatible modes are used; need to evaluate the stiffness matrix
       !! and the associated static condensation matrices A1 and A2
       !! in order to compute a consistent load vector
       call ffl_getmat (EMOD,RNY,E(1),IEL,ierr)
       if (ierr < 0) then
          call ReportError (error_p,'Error retrieving element material data')
          goto 990
       end if
       call isoMat3D (EMOD,RNY,E)
       IOP = 1 ! Internal DOFs are included
       call HEXA31 (EK(1,1),A1(1,1),A2(1,1),E(1),XG(1),YG(1),ZG(1), &
            &       NIP,IOP,IEL,LPU,IPSW,IERR)
       IF (IERR /= 0) THEN
          call ReportError (error_p,'Error computing element stiffness matrix')
          goto 990
       END IF
    else
       IOP = 0 ! No internal DOFs are included
    end if

    !! Calculate the element load vector, assuming constant pressure
    !! (the HEXA33 routine does not support a linearly varying pressure)
    PX = 0.25_dp*sum(p(1,:))
    PY = 0.25_dp*sum(p(2,:))
    PZ = 0.25_dp*sum(p(3,:))
    call HEXA33 (ES(1),A2(1,1),A1(1,1),XG(1),YG(1),ZG(1),PX,PY,PZ,IFACE,IOP)
    return

990 call ReportError (debugFileOnly_p,'LOAD44')

  end subroutine LOAD44


  subroutine LOAD45 (IEL, IFACE, P, NEDOF, ES, IPSW, LPU, IERR)

    !!==========================================================================
    !! This routine computes the element load vector for the volume element
    !! CSTET (Constant Strain Tetrahedron, ELTYPE = 45), due to
    !! a pressure load on one of its faces.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 3 Apr 2008 / 1.0
    !!==========================================================================

    use KindModule             , only : dp
    use CSTetraModule          , only : CSTetPload
    use ReportErrorModule      , only : reportError, error_p, debugFileOnly_p
    use FFlLinkHandlerInterface, only : ffl_getCoor

    integer , intent(in)  :: iel, iface, ipsw, lpu
    real(dp), intent(in)  :: p(:)
    real(dp), intent(out) :: es(:)
    integer , intent(out) :: nedof, ierr

    !! Local variables
    real(dp) :: XG(4), YG(4), ZG(4)

    !! --- Logic section ---

    nedof = 12

    !! Get global nodal coordinates for the element
    call ffl_getcoor (XG,YG,ZG,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element nodal coordinates')
       goto 990
    end if

    !! Calculate the element load vector
    call CSTetPload (XG,YG,ZG,P,IFACE,ES,LPU,IPSW,ierr)
    if (ierr == 0) return

990 call ReportError (debugFileOnly_p,'LOAD45')

  end subroutine LOAD45


  subroutine LOAD46 (IEL, IFACE, P, NEDOF, ES, IPSW, LPU, IERR)

    !!==========================================================================
    !! This routine computes the element load vector for the volume element
    !! IPRI6 (Isoparametric PRIsm element with 6 nodes, ELTYPE = 46), due to
    !! a pressure load on one of its faces.
    !!
    !! Programmer : Knut Morten Okstad
    !! date/rev   : 3 Apr 2008 / 1.0
    !!==========================================================================

    use KindModule             , only : dp
    use Ipri6Module            , only : Ipri6Pload
    use ReportErrorModule      , only : reportError, error_p, debugFileOnly_p
    use FFlLinkHandlerInterface, only : ffl_getCoor

    integer , intent(in)  :: iel, iface, ipsw, lpu
    real(dp), intent(in)  :: p(:)
    real(dp), intent(out) :: es(:)
    integer , intent(out) :: nedof, ierr

    !! Local variables
    real(dp) :: XG(6), YG(6), ZG(6)

    !! --- Logic section ---

    nedof = 18

    !! Get global nodal coordinates for the element
    call ffl_getcoor (XG,YG,ZG,IEL,ierr)
    if (ierr < 0) then
       call ReportError (error_p,'Error retrieving element nodal coordinates')
       goto 990
    end if

    !! Calculate the element load vector
    call Ipri6Pload (XG,YG,ZG,P,IFACE,ES,LPU,IPSW,ierr)
    if (ierr == 0) return

990 call ReportError (debugFileOnly_p,'LOAD46')

  end subroutine LOAD46

end module InaddModule
