module groups
  implicit none
contains

subroutine grouping(group,eqamat,k,na,a,x,symbol,natom,site)
  use functions
  use structure
  implicit none
  real(8) :: a(:,:),x(:,:)
  integer(1) :: site
  integer(2) :: eqamat(:,:),natom(:),na,k,ib,ie
  integer(2) :: no,ng,i,j,o,ic,order(na),deg(na)
  integer(2),allocatable :: group(:)
  character(len=2) :: symbol(:)
  logical(1) :: occ(na)

  no=size(eqamat,2)
  occ=.true.
  ng=0
  deg=0
  j=0
  do i=1,na
    if ( occ(i) .eqv. .false. ) cycle
    ng=ng+1
    do o=1,no
      ic=eqamat(i,o)
      if ( occ(ic) .eqv. .true. ) then
        j=j+1
        order(j)=ic
        deg(ng)=deg(ng)+1
        occ(ic)=.false.
      end if
    end do
  end do

  ib=1
  ie=0
  do i=1,ng
    ie=ie+deg(i)
    call sort(order(ib:ie),1_2,deg(i))
    ib=ib+deg(i)
  end do
  call system('rm SPOSCAR_NEW 2> /dev/null')
  do i=1,ng-1
    j=sum(deg(1:i))
    if ( order(j) /= j ) then
      call reorder(order,eqamat,natom,x,site)
      write(*,'(A)') '  NOTE: The order of atomic positions has been changed !'
      call writepos(a,x,symbol,natom,'SPOSCAR_NEW')
      exit
    end if
  end do

  o=na-k+1
  do i=2,ng
    if ( sum(deg(1:i-1)) >= o ) exit
  end do
  ng=i-1
  allocate(group(ng+1))
  group(1)=1
  do i=2,ng+1
    group(i)=sum(deg(1:i-1))+1
  end do
  return
end subroutine grouping

subroutine reorder(order,eqamat,natom,x,site)
  implicit none
  real(8) :: x(:,:)
  integer(1) :: site
  integer(2) :: order(:),eqamat(:,:),natom(:)
  integer(2),allocatable :: iorder(:)
  integer(2) :: no,na,i,ib,ie

  no=size(eqamat,2)
  na=size(order)
  allocate(iorder(na))
  forall(i=1:na)
    iorder(order(i))=i
  end forall
  do i=1,no
    eqamat(:,i)=iorder(eqamat(order,i))
  end do
  ib=sum(natom(1:site-1))+1
  ie=ib+natom(site)-1
  x(:,ib:ie)=x(:,order+ib-1)
  return
end subroutine reorder

end module groups
