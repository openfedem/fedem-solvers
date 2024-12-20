!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

!> @file manipMatrixModule.f90
!> @brief Utilities for manipulation of integer and real matrices.

!!==============================================================================
!> @brief Module with subroutines for manipulation of integer and real matrices.
!> @details This module contains various utility subroutines (and functions)
!> for manipulation of matrices, such as specialiced multiplication routines,
!> printing to file, etc.
!>
!> @brief Karl Erik Thoresen
!> @date 27 Sep 1998

module ManipMatrixModule

  implicit none

  !> @brief Standard routine for writing an object to file.
  interface WriteObject
     module procedure WriteIntArray
#if !defined(FT_HAS_INT4_ONLY) && !defined(FT_USE_INT8)
     module procedure WriteInt8Array
#endif
     module procedure WriteRealArray
     module procedure WriteRealMatrix
     module procedure WriteDoubleArray
     module procedure WriteDoubleMatrix
  end interface

  !> @brief Matrix multiplication involving diagonal matrices.
  interface diagMatMul
     module procedure diagMatMulMat
     module procedure diagMatMulVec
  end interface

  !> @brief Matrix multiplication involving 3&times;4 position matrices.
  interface matmul34
     module procedure matmul34Vec
     module procedure matmul34Mat
  end interface

  !> @brief Dyadic product between two vectors.
  interface dyadic_product
     module procedure traVmulV
     module procedure traUmulV
  end interface

  private :: WriteIntArray, WriteRealArray, WriteRealMatrix
  private :: WriteDoubleArray, WriteDoubleMatrix
  private :: diagMatMulMat, diagMatMulVec, matmul34Vec, matmul34Mat
  private :: traVmulV, traUmulV


contains

  !!============================================================================
  !> @brief Returns the cross product between the vectors y and z.
  !>
  !> @author Karl Erik Thoresen
  !> @date Jan 2000

  function cross_product (y,z)

    use KindModule, only : dp

    real(dp), intent(in) :: y(3), z(3)
    real(dp)             :: cross_product(3)

    !! --- Logic section ---

    cross_product(1) = y(2)*z(3) - y(3)*z(2)
    cross_product(2) = y(3)*z(1) - y(1)*z(3)
    cross_product(3) = y(1)*z(2) - y(2)*z(1)

  end function cross_product


  !!============================================================================
  !> @brief Returns the dyadic vector product `x`<sup>T</sup>`x`.
  !>
  !> @author Karl Erik Thoresen
  !> @date Aug 1999

  function traVmulV (x)

    use KindModule, only : dp

    real(dp), intent(in) :: x(:)
    real(dp)             :: traVmulV(size(x),size(x))

    !! Local variables
    integer :: i, j, n

    !! --- Logic section ---

    n = size(x)
    do i = 1, n
       do j = 1, n
          traVmulV(i,j) = x(i)*x(j)
       end do
    end do

  end function traVmulV


  !!============================================================================
  !> @brief Returns the dyadic vector product `x`<sup>T</sup>`y`.
  !>
  !> @author Knut Morten Okstad
  !> @date Jun 2004

  function traUmulV (x,y)

    use KindModule, only : dp

    real(dp), intent(in) :: x(:), y(:)
    real(dp)             :: traUmulV(size(x),size(y))

    !! Local variables
    integer :: i, j, m, n

    !! --- Logic section ---

    m = size(x)
    n = size(y)
    do i = 1, m
       do j = 1, n
          traUmulV(i,j) = x(i)*y(j)
       end do
    end do

  end function traUmulV


  !!============================================================================
  !> @brief Performs the transformation
  !> `T = A`<sup>T</sup>&lowast;`D`&lowast;`A`
  !> where `D` is a diagonal matrix.
  !>
  !> @author Knut Morten Okstad
  !> @date Aug 2001

  subroutine diagTransform (D,A,T,ierr)

    use KindModule, only : dp

    real(dp), intent(in)  :: D(:), A(:,:)
    real(dp), intent(out) :: T(:,:)
    integer , intent(out) :: ierr

    !! Local variables
    integer  :: i, j, k, m, n
    real(dp) :: s

    !! --- Logic section ---

    n = size(D)
    m = size(A,2)
    if (size(A,1) /= n .or. size(T,1) /= m .or. size(T,2) /= m) then
       ierr = -1
       return
    end if

    ierr = 0
    do i = 1, m
       do j = i, m
          s = 0.0_dp
          do k = 1, n
             s = s + A(k,i)*D(k)*A(k,j)
          end do
          T(i,j) = s
          if (j > i) T(j,i) = s
       end do
    end do

  end subroutine diagTransform


  !!============================================================================
  !> @brief Performs the matrix multiplication
  !> `B = c1`&lowast;`B + c2`&lowast;`D`&lowast;`A`.
  !> @details `D` is a diagonal matrix, whereas `B` and `A` are full matrices.
  !> The value of the parameters `c1` and `c2` are determined from the input
  !> variable `iflag`, as follows:
  !> - iflag =  1 : c1 = 0, c2 =  1
  !> - iflag = -1 : c1 = 0, c2 = -1
  !> - iflag =  2 : c1 = 1, c2 =  1
  !> - iflag = -2 : c1 = 1, c2 = -1
  !>
  !> @author Knut Morten Okstad
  !> @date Jan 2003

  subroutine diagMatMulMat (D,A,B,iflag,ierr)

    use KindModule, only : dp

    real(dp), intent(in)  :: D(:), A(:,:)
    real(dp), intent(out) :: B(:,:)
    integer , intent(in)  :: iflag
    integer , intent(out) :: ierr

    !! Local variables
    integer :: i, n

    !! --- Logic section ---

    n = size(D)
    if (size(A,1) /= n .or. size(B,1) /= n .or. size(A,2) /= size(B,2)) then
       ierr = -1
       return
    end if

    ierr = 0
    do i = 1, n
       select case (iflag)
       case ( 1); B(i,:) =          D(i)*A(i,:)
       case (-1); B(i,:) =        - D(i)*A(i,:)
       case ( 2); B(i,:) = B(i,:) + D(i)*A(i,:)
       case (-2); B(i,:) = B(i,:) - D(i)*A(i,:)
       end select
    end do

  end subroutine diagMatMulMat


  !!============================================================================
  !> @brief Performs the matrix-vector multiplication
  !> `Y = c1`&lowast;`Y + c2`&lowast;`D`&lowast;`X`.
  !> @details `D` is a diagonal matrix, whereas `Y` and `X` are vectors.
  !> The value of the parameters `c1` and `c2` are determined from the
  !> input variable `iflag`, as follows:
  !> - iflag =  1 : c1 = 0, c2 =  1
  !> - iflag = -1 : c1 = 0, c2 = -1
  !> - iflag =  2 : c1 = 1, c2 =  1
  !> - iflag = -2 : c1 = 1, c2 = -1
  !>
  !> @author Knut Morten Okstad
  !> @date Jan 2003

  subroutine diagMatMulVec (D,X,Y,iflag,ierr)

    use KindModule, only : dp

    real(dp), intent(in)  :: D(:), X(:)
    real(dp), intent(out) :: Y(:)
    integer , intent(in)  :: iflag
    integer , intent(out) :: ierr

    !! Local variables
    integer :: i, n

    !! --- Logic section ---

    n = size(D)
    if (size(X) /= n .or. size(Y) /= n) then
       ierr = -1
       return
    end if

    ierr = 0
    do i = 1, n
       select case (iflag)
       case ( 1); Y(i) =        D(i)*X(i)
       case (-1); Y(i) =      - D(i)*X(i)
       case ( 2); Y(i) = Y(i) + D(i)*X(i)
       case (-2); Y(i) = Y(i) - D(i)*X(i)
       end select
    end do

  end subroutine diagMatMulVec


  !!============================================================================
  !> @brief Returns the matrix product `a`&lowast;`b`.
  !> @details `a` and `b` are both 3&times;4 matrices with the assumed dimension
  !> 4&times;4, and the 4th row being `[0 0 0 1]`.
  !>
  !> @author Karl Erik Thoresen
  !> @date Feb 1999

  function matmul34Mat (a,b)

    use KindModule, only : dp

    real(dp), intent(in) :: a(3,4), b(3,4)
    real(dp)             :: matmul34Mat(3,4)

    !! --- Logic section ---

    matmul34Mat(1:3,1:3) = matmul(a(1:3,1:3),b(1:3,1:3))
    matmul34Mat(1:3,4)   = matmul(a(1:3,1:3),b(1:3,4)) +  a(1:3,4)

  end function matmul34Mat


  !!============================================================================
  !> @brief Returns the matrix-vector product `a`&lowast;`b`.
  !> @details `a` is a 3&times;4 matrices with the assumed dimension
  !> 4&times;4 and the 4th row being `[0 0 0 1]`,
  !> whereas `b` is a vector with size 3 and assumed size 4.
  !>
  !> @author Karl Erik Thoresen
  !> @date Feb 1999

  function matmul34Vec (a,b)

    use KindModule, only : dp

    real(dp), intent(in) :: a(3,4), b(3)
    real(dp)             :: matmul34Vec(3)

    !! --- Logic section ---

    matmul34Vec = matmul(a(:,1:3),b) + a(:,4)

  end function matmul34Vec


  !!============================================================================
  !> @brief Returns the inverse of the 2&times;2 matrix `a`.
  !>
  !> @author Knut Morten Okstad
  !> @date Nov 2000

  function invert22 (a,lpu,ierr) result(b)

    use KindModule, only : dp, epsDiv0_p, hugeVal_p

    real(dp), intent(in)            :: a(2,2)
    integer , intent(in) , optional :: lpu
    integer , intent(out), optional :: ierr
    real(dp)                        :: b(2,2)

    !! Local variables
    real(dp) :: det

    !! --- Logic section ---

    det = a(1,1)*a(2,2) - a(2,1)*a(1,2)

    if (abs(det) < epsDiv0_p) then
       if (present(lpu)) then
          write(lpu,*) '*** Cannot invert singular 2x2 matrix'
          call writeObject(a,lpu)
       else
          write(*,*) '*** Cannot invert singular 2x2 matrix'
       end if
       if (present(ierr)) ierr = -1
       b = hugeVal_p
       return
    end if

    b(1,1) =  a(2,2)/det
    b(2,2) =  a(1,1)/det
    b(1,2) = -a(1,2)/det
    b(2,1) = -a(2,1)/det

    if (present(ierr)) ierr = 0

  end function invert22


  !!============================================================================
  !> @brief Returns the inverse of the 3&times;3 matrix `a`.
  !>
  !> @author Knut Morten Okstad
  !> @date Nov 2000

  function invert33 (a,lpu,ierr) result(b)

    use KindModule, only : dp, epsDiv0_p, hugeVal_p

    real(dp), intent(in)            :: a(3,3)
    integer , intent(in) , optional :: lpu
    integer , intent(out), optional :: ierr
    real(dp)                        :: b(3,3)

    !! Local variables
    real(dp) :: det

    !! --- Logic section ---

    det = a(1,1)*(a(2,2)*a(3,3) - a(3,2)*a(2,3)) &
         -a(1,2)*(a(2,1)*a(3,3) - a(3,1)*a(2,3)) &
         +a(1,3)*(a(2,1)*a(3,2) - a(3,1)*a(2,2))

    if (abs(det) < epsDiv0_p) then
       if (present(lpu)) then
          write(lpu,*) '*** Cannot invert singular 3x3 matrix'
          call writeObject(a,lpu)
       else
          write(*,*) '*** Cannot invert singular 3x3 matrix'
       end if
       if (present(ierr)) ierr = -1
       b = hugeVal_p
       return
    end if

    b(1,1) =  (a(2,2)*a(3,3) - a(3,2)*a(2,3))/det
    b(1,2) = -(a(1,2)*a(3,3) - a(3,2)*a(1,3))/det
    b(1,3) =  (a(1,2)*a(2,3) - a(2,2)*a(1,3))/det
    b(2,1) = -(a(2,1)*a(3,3) - a(3,1)*a(2,3))/det
    b(2,2) =  (a(1,1)*a(3,3) - a(3,1)*a(1,3))/det
    b(2,3) = -(a(1,1)*a(2,3) - a(2,1)*a(1,3))/det
    b(3,1) =  (a(2,1)*a(3,2) - a(3,1)*a(2,2))/det
    b(3,2) = -(a(1,1)*a(3,2) - a(3,1)*a(1,2))/det
    b(3,3) =  (a(1,1)*a(2,2) - a(2,1)*a(1,2))/det

    if (present(ierr)) ierr = 0

  end function invert33


  !!============================================================================
  !> @brief Returns the inverse of the 3&times;4 transformation matrix `a`.
  !> @details The matrix `a` is assumed to have the implicit dimension 4&times;4
  !> with the 4th line being equal to `[0 0 0 1]`.
  !>
  !> @author Knut Morten Okstad
  !> @date Nov 2000

  function invert34 (a) result(b)

    use KindModule, only : dp

    real(dp), intent(in) :: a(3,4)
    real(dp)             :: b(3,4)

    !! --- Logic section ---

    b(:,1:3) = transpose(a(:,1:3))
    b(:,4) = -matmul(b(:,1:3),a(:,4))

  end function invert34


  !!============================================================================
  !> @brief Returns a 3&times;4 transformation matrix calculated from 3 points.
  !>
  !> @author Knut Morten Okstad
  !> @date 12 Jan 2001

  function trans3P (P1,P2,P3,lpu,ierr) result(tr)

    use KindModule, only : dp, epsDiv0_p, hugeVal_p

    real(dp), intent(in)            :: P1(3), P2(3), P3(3)
    integer , intent(in) , optional :: lpu
    integer , intent(out), optional :: ierr
    real(dp)                        :: tr(3,4)

    !! Local variables
    real(dp) :: vlength

    !! --- Logic section ---

    tr(:,1) = P2 - P1
    tr(:,2) = P3 - P1
    tr(:,3) = cross_product(tr(:,1),tr(:,2))
    tr(:,4) = P1

    vlength = len3sqr(tr(:,1))
    if (vlength > epsDiv0_p) then
       tr(:,1) = tr(:,1) / sqrt(vlength)
    else
       if (present(ierr)) ierr = -1
       if (present(lpu)) then
          write(lpu,*) '*** trans3P: Point 1 and 2 are coinciding'
          write(lpu,"(14X,'P1 =',1P3E13.5)") P1
          write(lpu,"(14X,'P2 =',1P3E13.5)") P2
          write(lpu,"(11X,'P2-P1 =',1P3E13.5)") P2-P1
       else
          write(*,*) '*** trans3P: Point 1 and 2 are coinciding'
       end if
       tr = hugeVal_p
       return
    end if

    vlength = len3sqr(tr(:,3))
    if (vlength > epsDiv0_p) then
       tr(:,3) = tr(:,3) / sqrt(vlength)
       tr(:,2) = cross_product(tr(:,3),tr(:,1))
    else
       if (len3sqr(tr(:,2)) <= epsDiv0_p) then
          if (present(ierr)) ierr = -2
          if (present(lpu)) then
             write(lpu,*) '*** trans3P: Point 1 and 3 are coinciding'
             write(lpu,"(14X,'P1 =',1P3E13.5)") P1
             write(lpu,"(14X,'P3 =',1P3E13.5)") P3
          else
             write(*,*) '*** trans3P: Point 1 and 3 are coinciding'
          end if
       else if (len3sqr(P3-P2) <= epsDiv0_p) then
          if (present(ierr)) ierr = -3
          if (present(lpu)) then
             write(lpu,*) '*** trans3P: Point 2 and 3 are coinciding'
             write(lpu,"(14X,'P2 =',1P3E13.5)") P2
             write(lpu,"(14X,'P3 =',1P3E13.5)") P3
          else
             write(*,*) '*** trans3P: Point 2 and 3 are coinciding'
          end if
       else
          if (present(ierr)) ierr = -4
          if (present(lpu)) then
             write(lpu,*) '*** trans3P: The three points are on a straight line'
             write(lpu,"(14X,'P1 =',1P3E13.5)") P1
             write(lpu,"(14X,'P2 =',1P3E13.5)") P2
             write(lpu,"(14X,'P3 =',1P3E13.5)") P3
          else
             write(*,*) '*** trans3P: The three points are on a straight line'
          end if
       end if
       if (present(lpu)) then
          write(lpu,"(11X,'P2-P1 =',1P3E13.5)") P2-P1
          write(lpu,"(11X,'P3-P1 =',1P3E13.5)") P3-P1
          write(lpu,"(11X,'P3-P2 =',1P3E13.5)") P3-P2
          write(lpu,"(' (P2-P1)x(P3-P1) =',1P3E13.5)") tr(:,3)
       end if
       tr = hugeVal_p
       return
    end if

    if (present(ierr)) ierr = 0

  contains

    !> brief Calculates the square of the length of a 3D vector.
    real(dp) function len3sqr (v)
      real(dp), intent(in) :: v(3)
      len3sqr = v(1)*v(1) + v(2)*v(2) + v(3)*v(3)
    end function len3sqr

  end function trans3P


  !!============================================================================
  !> @brief Returns a 3&times;3 transformation matrix based on the given
  !> X-direction vector `V1`.
  !>
  !> @author Knut Morten Okstad
  !> @date 28 Jun 2002

  function trans1V (V1,lpu,ierr) result(T)

    use KindModule, only : dp, epsDiv0_p, hugeVal_p

    real(dp), intent(in)            :: V1(3)
    integer , intent(in) , optional :: lpu
    integer , intent(out), optional :: ierr
    real(dp)                        :: T(3,3)

    !! Local variables
    real(dp) :: sqrLen

    !! --- Logic section ---

    sqrLen = V1(1)*V1(1) + V1(2)*V1(2) + V1(3)*V1(3)
    if (abs(sqrLen) > epsDiv0_p) then
       T(:,1) = V1 / sqrt(sqrLen)
    else
       if (present(lpu)) then
          write(lpu,*) '*** trans1V: The X-direction vector has zero length'
          write(lpu,"(14X,'V1 = ',1P3E13.5)") V1
       else
          write(*,*) '*** trans1V: The X-direction vector has zero length'
       end if
       if (present(ierr)) ierr = -1
       T = hugeVal_p
       return
    end if

    if (abs(T(3,1)) > abs(T(2,1))) then
       ! Define V2 by projecting the global Y-axis onto the plane, then V3=V1xV2
       T(1,2) = -T(2,1)*T(1,1)
       T(2,2) =  T(1,1)*T(1,1) + T(3,1)*T(3,1)
       T(3,2) = -T(2,1)*T(3,1)
       T(:,2) =  T(:,2) / sqrt(T(1,2)*T(1,2) + T(2,2)*T(2,2) + T(3,2)*T(3,2))
       T(:,3) =  cross_product(T(:,1),T(:,2))
    else
       ! Define V3 by projecting the global Z-axis onto the plane, then V2=V3xV1
       T(1,3) = -T(3,1)*T(1,1)
       T(2,3) = -T(3,1)*T(2,1)
       T(3,3) =  T(1,1)*T(1,1) + T(2,1)*T(2,1)
       T(:,3) =  T(:,3) / sqrt(T(1,3)*T(1,3) + T(2,3)*T(2,3) + T(3,3)*T(3,3))
       T(:,2) =  cross_product(T(:,3),T(:,1))
    end if

    if (present(ierr)) ierr = 0

  end function trans1V


  !!============================================================================
  !> @brief Writes out and integer array in a formatted way.
  !> @param[in] array The array to print out
  !> @param[in] lpu File unit number to print to
  !> @param[in] text Heading to be printed before the matrix values

  subroutine WriteIntArray (array,lpu,text)

    integer     , intent(in)           :: array(:)
    integer     , intent(in)           :: lpu
    character(*), intent(in), optional :: text

    !! Local variables
    integer :: i, n

    !! --- Logic section ---

    if (present(text)) write(lpu,'(A)') trim(text)

    n = size(array)
    do i = 1, n, 10
       write(lpu,"(i6,' :',10I6)") i,array(i:min(i+9,n))
    end do

  end subroutine WriteIntArray


#if !defined(FT_HAS_INT4_ONLY) && !defined(FT_USE_INT8)
  !!============================================================================
  !> @brief Writes out an integer array in a formatted way.
  !> @param[in] array The array to print out
  !> @param[in] lpu File unit number to print to
  !> @param[in] text Heading to be printed before the matrix values

  subroutine WriteInt8Array (array,lpu,text)

    use KindModule, only : i8

    integer(i8) , intent(in)           :: array(:)
    integer     , intent(in)           :: lpu
    character(*), intent(in), optional :: text

    !! Local variables
    integer :: i, n

    !! --- Logic section ---

    if (present(text)) write(lpu,'(A)') text

    n = size(array)
    do i = 1, n, 10
       write(lpu,"(i6,' :',10I6)") i,array(i:min(i+9,n))
    end do

  end subroutine WriteInt8Array
#endif


  !!============================================================================
  !> @brief Writes out a single precision array in a formatted way.
  !> @param[in] array The array to print out
  !> @param[in] lpu File unit number to print to
  !> @param[in] text Heading to be printed before the array elements
  !> @param[in] nellIn Max number of values to print per line
  !> @param[in] eps Zero tolerance, truncate terms with smaller absolute value

  subroutine WriteRealArray (array,lpu,text,nellIn,eps)

    use KindModule, only : sp, epsDiv0_p

    real(sp)    , intent(in)           :: array(:)
    integer     , intent(in)           :: lpu
    integer     , intent(in), optional :: nellIn
    character(*), intent(in), optional :: text
    real(sp)    , intent(in), optional :: eps

    !! Local variables
    integer  :: i, j, n, nell
    real(sp) :: epsZero

    !! --- Logic section ---

    if (present(nellIn)) then
       nell = min(20,nellIn) ! Number of array elements to write per line
    else
       nell = 6
    end if

    if (present(eps)) then
       epsZero = eps
    else
       epsZero = epsDiv0_p
    end if

    if (present(text)) write(lpu,'(A)') trim(text)

    n = size(array)
    do i = 1, n, nell
       write(lpu,'(i6)',ADVANCE='no') i
       do j = i, min(i+nell-1,n)
          if (abs(array(j)) > epsZero) then
             write(lpu,'(ES15.6E3)',ADVANCE='no') array(j)
          else
             write(lpu,'(F11.1,4X)',ADVANCE='no') 0.0
          end if
       end do
       write(lpu,*)
    end do

  end subroutine WriteRealArray


  !!============================================================================
  !> @brief Writes out a 2D single precision array in a formatted way.
  !> @param[in] array The array to print out
  !> @param[in] lpu File unit number to print to
  !> @param[in] text Heading to be printed before the matrix elements
  !> @param[in] nellIn Max number of values to print per line
  !> @param[in] eps Zero tolerance, truncate terms with smaller absolute value

  subroutine WriteRealMatrix (array,lpu,text,nellIn,eps)

    use KindModule, only : sp, epsDiv0_p

    real(sp)    , intent(in)           :: array(:,:)
    integer     , intent(in)           :: lpu
    integer     , intent(in), optional :: nellIn
    character(*), intent(in), optional :: text
    real(sp)    , intent(in), optional :: eps

    !! Local variables
    integer  :: i, j, k, m, n, nell
    real(sp) :: epsZero

    !! --- Logic section ---

    if (present(nellIn)) then
       nell = min(20,nellIn) ! Number of matrix elements to write per line
    else
       nell = 6
    end if

    if (present(eps)) then
       epsZero = eps
    else
       epsZero = epsDiv0_p
    end if

    if (present(text)) write(lpu,'(A)') trim(text)

    m = size(array,1)
    n = size(array,2)
    do j = 1, n, nell
       write(lpu,'(1x,20i15)') (i,i=j,min(j+nell-1,n))
       do i = 1, m
          write(lpu,'(i6)',ADVANCE='no') i
          do k = j, min(j+nell-1,n)
             if (abs(array(i,k)) > epsZero) then
                write(lpu,'(ES15.6E3)',ADVANCE='no') array(i,k)
             else
                write(lpu,'(F11.1,4X)',ADVANCE='no') 0.0
             end if
          end do
          write(lpu,*)
       end do
    end do

  end subroutine WriteRealMatrix


  !!============================================================================
  !> @brief Writes out a double precision array in a formatted way.
  !> @param[in] array The array to print out
  !> @param[in] lpu File unit number to print to
  !> @param[in] text Heading to be printed before the array elements
  !> @param[in] nellIn Max number of values to print per line
  !> @param[in] eps Zero tolerance, truncate terms with smaller absolute value

  subroutine WriteDoubleArray (array,lpu,text,nellIn,eps)

    use KindModule, only : dp, epsDiv0_p

    real(dp)    , intent(in)           :: array(:)
    integer     , intent(in)           :: lpu
    integer     , intent(in), optional :: nellIn
    character(*), intent(in), optional :: text
    real(dp)    , intent(in), optional :: eps

    !! Local variables
    integer  :: i, j, n, nell
    real(dp) :: epsZero

    !! --- Logic section ---

    if (present(nellIn)) then
       nell = min(20,nellIn) ! Number of array elements to write per line
    else
       nell = 6
    end if

    if (present(eps)) then
       epsZero = eps
    else
       epsZero = epsDiv0_p
    end if

    if (present(text)) write(lpu,'(A)') trim(text)

    n = size(array)
    do i = 1, n, nell
       write(lpu,'(i6)',ADVANCE='no') i
       do j = i, min(i+nell-1,n)
          if (abs(array(j)) > epsZero) then
             write(lpu,'(ES15.6E3)',ADVANCE='no') array(j)
          else
             write(lpu,'(F11.1,4X)',ADVANCE='no') 0.0
          end if
       end do
       write(lpu,*)
    end do

  end subroutine WriteDoubleArray


  !!============================================================================
  !> @brief Writes out a 2D double precision array in a formatted way.
  !> @param[in] array The array to print out
  !> @param[in] lpu File unit number to print to
  !> @param[in] text Heading to be printed before the matrix elements
  !> @param[in] nellIn Max number of values to print per line
  !> @param[in] eps Zero tolerance, truncate terms with smaller absolute value

  subroutine WriteDoubleMatrix (array,lpu,text,nellIn,eps)

    use KindModule, only : dp, epsDiv0_p

    real(dp)    , intent(in)           :: array(:,:)
    integer     , intent(in)           :: lpu
    integer     , intent(in), optional :: nellIn
    character(*), intent(in), optional :: text
    real(dp)    , intent(in), optional :: eps

    !! Local variables
    integer  :: i, j, k, m, n, nell
    real(dp) :: epsZero

    !! --- Logic section ---

    if (present(nellIn)) then
       nell = min(20,nellIn) ! Number of matrix elements to write per line
    else
       nell = 6
    end if

    if (present(eps)) then
       epsZero = eps
    else
       epsZero = epsDiv0_p
    end if

    if (present(text)) write(lpu,'(A)') trim(text)

    m = size(array,1)
    n = size(array,2)
    do j = 1, n, nell
       write(lpu,'(1x,20i15)') (i,i=j,min(j+nell-1,n))
       do i = 1, m
          write(lpu,'(i6)',ADVANCE='no') i
          do k = j, min(j+nell-1,n)
             if (abs(array(i,k)) > epsZero) then
                write(lpu,'(ES15.6E3)',ADVANCE='no') array(i,k)
             else
                write(lpu,'(F11.1,4X)',ADVANCE='no') 0.0
             end if
          end do
          write(lpu,*)
       end do
    end do

  end subroutine WriteDoubleMatrix


  !!============================================================================
  !> @brief Unifies the given matrix @b A.
  !> @details If @b A is not square, the greatest possible square matrix will be
  !> unified and the rest will be zeroed.

  subroutine unify (A)

    use KindModule, only : dp

    real(dp), intent(inout) :: A(:,:)

    !! Local variables
    integer :: i, n

    !! --- Logic section ---

    A = 0.0_dp
    n = min(size(A,1),size(A,2))
    do i = 1, n
       A(i,i) = 1.0_dp
    end do

  end subroutine unify


  !!============================================================================
  !> @brief Returns a matrix pointer to an integer array.
  !> @details @a array can be a one-dimensional array as actual argument,
  !> but is a matrix as dummy argument.

  function MatrixPointToArray_int (array,rows,columns)

    integer, intent(in)         :: rows, columns
    integer, intent(in), target :: array(rows,columns)
    integer, pointer            :: MatrixPointToArray_int(:,:)

    !! --- Logic section ---

    MatrixPointToArray_int => array

  end function MatrixPointToArray_int


  !!============================================================================
  !> @brief Returns a matrix pointer to a real array.
  !> @details @a array can be a one-dimensional array as actual argument,
  !> but is a matrix as dummy argument.

  function MatrixPointToArray_real (array,rows,columns)

    use KindModule, only : dp

    integer , intent(in)         :: rows, columns
    real(dp), intent(in), target :: array(rows,columns)
    real(dp), pointer            :: MatrixPointToArray_real(:,:)

    !! --- Logic section ---

    MatrixPointToArray_real => array

  end function MatrixPointToArray_real

end module ManipMatrixModule
