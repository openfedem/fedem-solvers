!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module DiffractionModule

  use KindModule, only : dp, sp

  implicit none

  private

  integer          , save        :: nDOF
  integer          , allocatable :: np(:), nf(:)
  real(sp)         , allocatable :: CoG(:,:), wave(:,:)
  character(len=32), allocatable :: cBody(:)

  real(sp)   , pointer :: A(:,:,:) => null()
  real(sp)   , pointer :: B(:,:,:) => null()
  complex(sp), pointer :: Fe(:,:)  => null()

  public :: NemohInit, NemohClose, NemohMesh, NemohAnalysis, getExcitationForce


contains

  function BodyName (id)

    !!==========================================================================
    !! Create a unique body file name based on a given superelement description.
    !!
    !! Programmer : Knut Morten Okstad              date/rev : 27 Feb 2015 / 1.0
    !!==========================================================================

    use IdTypeModule, only : IdType

    type(IdType), intent(in) :: id
    character(len=18)        :: BodyName

    !! Local variables
    integer          :: i, j
    character(len=8) :: baseId
    character(len=1) :: char

    !! --- Logic section ---

    if (id%descr .eq. '') then
       BodyName = 'Body_'
       j = 5
    else
       !! Create a string from the description that can be used as a filename
       j = 1
       do i = 1, len_trim(id%descr)
          char = id%descr(i:i)
          if ( (char .ge. 'a' .and. char .le. 'z') .or. &
               (char .ge. 'A' .and. char .le. 'Z') ) then
             BodyName(j:j) = char
          else if (char .eq. '_' .or. char .eq. ' ') then
             BodyName(j:j) = '_'
          else
             cycle
          end if
          j = j + 1
          if (j .eq. 10) exit
       end do
    end if

    !! Add the base id to create a unique filename
    write(baseId,'(I8)') id%baseId
    BodyName(j:) = adjustl(baseId)

  end function BodyName


  subroutine NemohInit (sups,ierr)

    !!==========================================================================
    !! Initializes the internal data structures for Nemoh analysis.
    !!
    !! Programmer : Knut Morten Okstad              date/rev : 27 Feb 2015 / 1.0
    !!==========================================================================

    use IdTypeModule     , only : getId
    use SupElTypeModule  , only : SupElType
    use reportErrorModule, only : allocationError, reportError, note_p, error_p

    type(SupElType), intent(in)  :: sups(:)
    integer        , intent(out) :: ierr

    !! Local variables
    integer :: i, nBN

    !! --- Logic section ---

    nBN = 0
    nDOF = 0
    do i = 1, size(sups)
       if (associated(sups(i)%hydyn)) then
          if (sups(i)%hydyn%bodyIndex >= nBN) nBN = sups(i)%hydyn%bodyIndex + 1
          if (sups(i)%hydyn%bodyIndex >= 0) then
             nDOF = nDOF + 6
             call reportError (note_p,'Diffraction analysis is performed '// &
                  &            'for Superelement'//getId(sups(i)%id))
          end if
       end if
    end do
    if (nDOF == 0) then
       ierr = 1 ! No wet bodies in this model
       call reportError (error_p,'No superelements with hydrodynamics are '// &
            &            'present in this model.', &
            &            'Diffraction analysis is NOT performed.', &
            &            addString='NemohInit')
       return
    end if

    allocate(cBody(nBN),CoG(3,nBN),np(nBN),nf(nBN),stat=ierr)
    if (ierr /= 0) then
       ierr = allocationError('NemohInit')
    else
       CoG = 0.0_sp
       np = 0
       nf = 0
    end if

  end subroutine NemohInit


  subroutine NemohClose

    !!==========================================================================
    !! Deallocates the global static arrays used in the Nemoh analysis.
    !!
    !! Programmer : Knut Morten Okstad              date/rev : 22 Nov 2015 / 1.0
    !!==========================================================================

    use FFaBodyHandlerInterface, only : ffa_erase_bodies

    !! --- Logic section ---

    call ffa_erase_bodies()

    nDOF = 0
    if (allocated(cBody)) deAllocate(cBody)
    if (allocated(np))    deAllocate(np)
    if (allocated(nf))    deAllocate(nf)
    if (allocated(CoG))   deAllocate(CoG)
    if (allocated(wave))  deAllocate(wave)
    if (associated(A))    deallocate(A)
    if (associated(B))    deallocate(B)
    if (associated(Fe))   deallocate(Fe)
    nullify(A)
    nullify(B)
    nullify(Fe)

  end subroutine NemohClose


  subroutine NemohMesh (prefix,sId,bodyIdx,supTr,Xc,rho,g,nfobj,ierr)

    !!==========================================================================
    !! Invokes the Nemoh meshing routine for a superelement geometry.
    !!
    !! Programmer : Knut Morten Okstad              date/rev : 27 Feb 2015 / 1.0
    !!==========================================================================

    use IdTypeModule           , only : IdType, getId
    use manipMatrixModule      , only : matmul34
#ifdef FT_HAS_NEMOH
    use Nemoh                  , only : Mesh
    use progressModule         , only : lterm
#endif
    use reportErrorModule      , only : error_p, debugFileOnly_p
    use reportErrorModule      , only : reportError, allocationError
    use FFaBodyHandlerInterface, only : ffa_get_nofaces, ffa_get_face

    character(len=*)  , intent(in)  :: prefix
    type(IdType)      , intent(in)  :: sId
    integer           , intent(in)  :: bodyIdx, nfobj
    real(dp)          , intent(in)  :: supTr(3,4), Xc(3)
    real(sp), optional, intent(in)  :: rho, g
    integer           , intent(out) :: ierr

    !! Local variables
    integer  :: i, j, k, nFace
    real(dp) :: face(3,4), v(3)
    real(sp), allocatable :: Xf(:,:,:)

    !! --- Logic section ---

    if (.not.allocated(cBody) .or. bodyIdx >= size(cBody)) then
       call reportError (error_p,'Invalid body index',addString='NemohMesh')
       ierr = -1
       return
    end if

    !! Get number of faces in this body
    call ffa_get_nofaces (bodyIdx,nFace)
    if (nFace < 4) then
       ierr = -2
       call reportError (error_p,'Invalid geometry for superelement'// &
            &            getId(sId),addString='NemohMesh')
    end if

    allocate(Xf(nFace,4,3),stat=ierr)
    if (ierr /= 0) then
       ierr = allocationError('NemohMesh')
       return
    end if

    !! Get coordinates of the face vertices
    do i = 1, nFace
       call ffa_get_face (bodyIdx,i-1,face,ierr)
       if (ierr < 0) goto 900
       if (ierr == 4) then
          do j = 1, 4
             v = matmul34(supTr,face(:,j))
             do k = 1, 3
                Xf(i,j,k) = real(v(k),sp)
             end do
          end do
       else
          ierr = -3
          call reportError (error_p,'Triangular faces are not supported, '// &
               &            'superelement'//getId(sId),addString='NemohMesh')
          return
       end if
    end do

#ifdef FT_HAS_NEMOH
    write(lterm,*)
    write(lterm,*) '--> Diffraction analysis for Superelement'//trim(getId(sId))

    !! Invoke the Nemoh meshing module for this body
    i = 1 + bodyIdx
    CoG(:,i) = matmul34(supTr,Xc)
    cBody(i) = BodyName(sId)
    call Mesh (prefix,trim(cBody(i)),Xf,nfobj,CoG(:,i),rho,g,nf(i),np(i),ierr)
#else
    print *,' ** NemohMesh dummy: ',trim(prefix),nfobj,Xc,rho,g
#endif
    if (ierr == 0) goto 990

900 call reportError (debugFileOnly_p,'NemohMesh')
990 deallocate(Xf)

  end subroutine NemohMesh


  subroutine NemohAnalysis (prefix,env,waveData,beta,M,omega,Kext,Cext,RAO, &
    &                       ipsw,ierr)

    !!==========================================================================
    !! Invokes the Nemoh preprocessor and solver for the floating bodies.
    !!
    !! Programmer : Knut Morten Okstad              date/rev : 04 Mar 2015 / 1.0
    !!==========================================================================

    use EnvironmentTypeModule, only : EnvironmentType
#ifdef FT_HAS_NEMOH
    use Nemoh, only : PreProcess, NemohSolve, PostProcess
    use Nemoh, only : writeNemohFile, getHydroStatics, getResults
    use reportErrorModule, only : error_p, reportError
#endif
    use reportErrorModule, only : allocationError, internalError, getErrorFile

    character(len=*)     , intent(in)    :: prefix
    type(EnvironmentType), intent(in)    :: env
    real(dp)   , optional, intent(in)    :: waveData(:,:)
    real(sp)   , optional, intent(in)    :: beta(:), M(:,:)
    real(sp)   , optional, intent(inout) :: omega(:)
    real(sp)   , optional, intent(in)    :: Kext(:), Cext(:)
    complex(sp), optional, intent(out)   :: RAO(:,:,:)
    integer    , optional, intent(in)    :: ipsw
    integer              , intent(out)   :: ierr

    !! Local variables
    integer  :: i, j, lpu, nFreq, iFreq
#ifdef FT_HAS_NEMOH
    real(sp) :: depth, rho, g
#endif
    real(sp), allocatable :: C(:,:)

    !! --- Logic section ---

    ierr = 0
    if (present(RAO)) RAO = 0.0_sp
    if (nDOF == 0 .or. .not.allocated(cBody)) return

    if (present(waveData)) then
       iFreq = max(size(waveData,1),3)
       nFreq = size(waveData,2)
    else if (present(omega)) then
       iFreq = 1
       nFreq = size(omega)
    else
       ierr = internalError('NemohAnalysis: No wave frequences provided')
       return
    end if
    if (present(M) .and. present(RAO)) then
       allocate(wave(nFreq,iFreq),C(nDOF,nDOF),stat=ierr)
    else
       allocate(wave(nFreq,iFreq),stat=ierr)
    end if
    if (ierr /= 0) then
       ierr = allocationError('NemohAnalysis')
       return
    end if

    if (present(waveData)) then
       do i = 1, iFreq
          do j = 1, nFreq
             wave(j,i) = real(waveData(i,j),sp)
          end do
       end do
       if (iFreq > 2) iFreq = 2
       if (present(omega) .and. nFreq > 0) then
          do j = 1, nFreq
             omega(j) = real(waveData(j,iFreq),sp)
          end do
          if (nFreq < size(omega)) then
             omega(nFreq+1:) = 0.0_sp
          end if
       end if
    else if (present(omega)) then
       wave(:,1) = omega
    end if

    lpu = 0
    if (present(M) .and. present(RAO) .and. present(ipsw)) then
       if (ipsw > 0) lpu = getErrorFile()
    end if
#ifdef FT_HAS_NEMOH
    if (present(M) .and. present(RAO)) then
       j = 1
       do i = 1, size(cBody)
          if (np(i) > 0) then
             call getHydroStatics (prefix,cBody(i),KH=C(j:j+5,j:j+5),ierr=ierr)
             j = j + 6
          end if
       end do
       if (lpu > 0) then
          call writeMat (C,lpu,'Hydrostatic stiffness matrix (from Nemoh)')
          call writeMat (M,lpu,'Structure mass matrix (from Fedem)')
       end if
    end if

    depth = env%seaDepth
    rho = env%rhoW
    g = sqrt(sum(env%gravity*env%gravity))
    call writeNemohFile (prefix,cBody,np,nf,CoG,depth,0.0_sp,0.0_sp,rho,g)
    call PreProcess (prefix,wave(:,iFreq),beta,ierr)
    if (ierr == 0) then
       call NemohSolve (prefix,ierr)
    end if
    if (ierr == 0) then
       if (present(M) .and. present(RAO)) then
          call PostProcess (prefix,C,M,Kext,Cext,RAO,.true.,lpu,ierr)
       else
          call PostProcess (prefix,ierr=ierr)
       end if
    end if
    if (ierr == 0 .and. .not.present(RAO)) then
       call getResults (prefix,A,B,Fe,ierr)
    end if
    if (ierr /= 0) then
       call reportError (error_p,'Diffraction analysis failed', &
            &            addString='NemohAnalysis')
    end if
#else
    if (present(ipsw)) then
       print *,' ** NemohAnalysis dummy: ',trim(prefix),env%rhoW, &
            &                              present(beta),present(Kext), &
            &                              present(Cext)
    end if
#endif
    deallocate(C)

#ifdef FT_HAS_NEMOH
  contains

    subroutine WriteMat (array,lpu,name)

      real(sp)        , intent(in)           :: array(:,:)
      integer         , intent(in), optional :: lpu
      character(len=*), intent(in), optional :: name

      integer , parameter :: nell = 6 ! Number of values to write per line
      real(sp), parameter :: epsZ = 1.0E-8
      integer  :: i, j, k, m, n

      write(lpu,'(/5X,A)') name

      m = size(array,1)
      n = size(array,2)
      do j = 1, n, nell
         write(lpu,'(1X,20I13)') (i,i=j,min(j+nell-1,n))
         do i = 1, m
            write(lpu,'(I6)',ADVANCE='no') i
            do k = j, min(j+nell-1,n)
               if (abs(array(i,k)) > epsZ) then
                  write(lpu,'(1P,E13.5)',ADVANCE='no') array(i,k)
               else
                  write(lpu,'(F9.1,4X)',ADVANCE='no') 0.0
               end if
            end do
            write(lpu,*)
         end do
      end do

    end subroutine WriteMat
#endif

  end subroutine NemohAnalysis


  subroutine getExcitationForce (ibody,t,Fex,ierr)

    !!==========================================================================
    !! Get excitation force from diffraction for a submerged body.
    !!
    !! Programmer : Knut Morten Okstad              date/rev : 16 Mar 2015 / 1.0
    !!==========================================================================

    integer , intent(in)  :: ibody
    real(dp), intent(in)  :: t
    real(dp), intent(out) :: Fex(6)
    integer , intent(out) :: ierr

    !! Local variables
    integer  :: i, j
    real(dp) :: A, phi

    !! --- Logic section ---

    ierr = 1
    Fex = 0.0_dp
    if (.not.associated(Fe) .or. .not.allocated(wave)) return
    if (ibody < 0 .or. ibody*6+6 > size(Fe,1) .or. size(wave,2) < 3) return

    ierr = 0
    do i = 1, size(wave,1)
       do j = 1, 6
          A = abs(Fe(ibody*6+j,i))
          phi = wave(i,3) + imag(Fe(ibody*6+j,i))
          Fex(j) = Fex(j) + wave(i,1)*A*sin(wave(i,2)*t-phi)
       end do
    end do

  end subroutine getExcitationForce

end module DiffractionModule
