module functions
  implicit none
contains

function binsearch_2(n,a)
  implicit none
  integer(2) :: binsearch_2,left,right,mid
  integer(2) :: n,a(:)

  left=1
  right=size(a)
  if ( n >= a(right) ) then
    binsearch_2=right
    return
  end if
  do while ( .true. )
    mid=(left+right)/2
    if ( left == mid .or. n == a(mid) ) exit
    if ( n > a(mid) ) then
      left=mid
    else
      right=mid
    end if
  end do
  binsearch_2=mid
  return
end function binsearch_2

function binsearch_8(n,a)
  implicit none
  integer(2) :: binsearch_8,left,right,mid
  integer(8) :: n,a(:)

  left=1
  right=size(a)
  if ( n >= a(right) ) then
    binsearch_8=right
    return
  end if
  do while ( .true. )
    mid=(left+right)/2
    if ( left == mid .or. n == a(mid) ) exit
    if ( n > a(mid) ) then
      left=mid
    else
      right=mid
    end if
  end do
  binsearch_8=mid
  return
end function binsearch_8

function nchoosek(n,k)
  implicit none
  integer(2) :: k,n,i,ia,ib
  integer(8) :: nchoosek
  real(8) :: A,B

  if ( k < n-k ) then
    ia=n-k+1
    ib=k
  else
    ia=k+1
    ib=n-k
  end if
  A=1;B=1
  do i=ia,n
    A=A*i
  end do
  do i=2,ib
    B=B*i
  end do
  nchoosek=A/B+0.5
  return
end function nchoosek

function complement(a,m)
  implicit none
  integer(2) :: i,j,m
  integer(2) :: a(:),c(m)
  integer(2) :: complement(m)

  c=1
  c(a)=0
  j=0
  complement=0
  do i=1,m
    if ( c(i) /= 0 ) then
      j=j+1
      complement(j)=i
    end if
  end do
  return
end function complement

function width(i_1,i_2,i_4,i_8)
  integer(1) :: width
  integer(1),optional :: i_1
  integer(2),optional :: i_2
  integer(4),optional :: i_4
  integer(8),optional :: i_8
  character(len=20) :: w

  if ( present(i_1) ) write(w,'(I20)') i_1
  if ( present(i_2) ) write(w,'(I20)') i_2
  if ( present(i_4) ) write(w,'(I20)') i_4
  if ( present(i_8) ) write(w,'(I20)') i_8
  width=len_trim(adjustl(w))
end function width

function norm(A)
  implicit none
  real(8) :: norm
  real(8) :: A(3)

  norm=dsqrt(A(1)**2+A(2)**2+A(3)**2)
  return
end function norm

function cross(A,B)
  implicit none
  real(8) :: cross(3)
  real(8) :: A(3),B(3)

  cross(1)=A(2)*B(3)-A(3)*B(2)
  cross(2)=A(3)*B(1)-A(1)*B(3)
  cross(3)=A(1)*B(2)-A(2)*B(1)
  return
end function cross

function trace(A)
  implicit none
  integer(1) :: trace,A(:,:)

  trace=A(1,1)+A(2,2)+A(3,3)
  return
end function trace

function det(A)
  implicit none
  integer(1) :: det
  integer(1) :: A(3,3)

  det=    A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
  det=det-A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
  det=det+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
  return
end function det

function det2(A)
  implicit none
  real(8) :: det2,A(2,2)

  det2=A(1,1)*A(2,2)-A(1,2)*A(2,1)
  return
end function det2

function det3(A)
  implicit none
  real(8) :: det3,A(3,3)

  det3=     A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))
  det3=det3-A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))
  det3=det3+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
  return
end function det3

function inverse(A)
  implicit none
  real(8) :: inverse(3,3)
  real(8) :: A(3,3)
  integer(1) :: ind(2,3),i,j

  ind(:,1)=(/ 2, 3 /)
  ind(:,2)=(/ 1, 3 /)
  ind(:,3)=(/ 1, 2 /)
  do i=1,3
    do j=1,3
      inverse(j,i)=(-1)**(i+j)*det2(A(ind(:,i),ind(:,j)))
    end do
  end do
  inverse=inverse/det3(A)
  return
end function inverse

recursive subroutine sort(A,left,right)
  implicit none
  integer(2),intent(in out) :: A(:)
  integer(2),intent(in) :: left,right
  integer(2) :: temp,t
  integer(2) :: i,j

  if ( left > right ) then
    return
  end if
  temp=A(left)
  i=left
  j=right
  do while ( i /= j )
    do while ( A(j) >= temp .and. i < j )
      j=j-1
    end do
    do while ( A(i) <= temp .and. i < j )
      i=i+1
    end do
    if ( i < j ) then
      t=A(i)
      A(i)=A(j)
      A(j)=t
    end if
  end do
  A(left)=A(i)
  A(i)=temp
  call sort(A,left,i-1_2)
  call sort(A,i+1_2,right)
end subroutine sort

end module functions
