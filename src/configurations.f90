module configurations
  use functions
  implicit none
contains

subroutine irrconfig(iconf,deg,eqamat,k,lpro,fast,cmax)
  use progress
  use stdout
  implicit none
  type :: datalink
    integer(2),allocatable :: iconf(:)
    integer(2) :: deg
    type(datalink),pointer :: next
  end type datalink
  integer(2),allocatable :: iconf(:,:),deg(:)
  integer(2) :: eqamat(:,:),k(:)
  integer(2),allocatable :: aconf(:),oconf(:)
  integer(8),allocatable :: ncm(:),mc(:),wc(:),E(:,:),ne(:,:)
  integer(2),allocatable :: lconf(:),m(:),neqa(:)
  integer(1) :: nk
  integer(2) :: na,no,ns,j,l,nl
  integer(4) :: n
  integer(8) :: i,nc,ic,se,o,cmax
  logical(1),allocatable :: occ(:)
  logical(1) :: lpro,fast,ldeg
  type(cls_prog) :: prog
  type(datalink),pointer :: p,head,next

  nk=size(k)-1
  na=sum(k)
  ns=na-k(nk+1)
  no=size(eqamat,1)
  allocate(lconf(ns))
  allocate(ncm(nk))
  allocate(wc(nk))
  allocate(m(nk))
  call initial(lconf,nc,ncm,wc,mc,m,k)
  if ( nc < 0 ) call stderr('The number of atomic configurations is out-of-range !')
  call matrixE(E,k,m)
  call neqatom(neqa,ne,eqamat,E,na,ncm(1),wc(1),k(1),fast,ldeg)
  call stdout_3(nc)
  allocate(aconf(ns))
  allocate(oconf(ns))
  allocate(head)
  nullify(head%next)
  p=>head
  n=0
!  fast=.false.
  if ( fast ) then
    call fast_mode
  else
    call std_mode
  end if
  call stdout_4(n)
  allocate(iconf(ns,n))
  allocate(deg(n))
  nullify(p%next)
  i=0
  p=>head
  do while ( associated(p%next) )
    p=>p%next
    i=i+1
    iconf(:,i)=p%iconf
    deg(i)=p%deg
  end do
  p=>head
  do while ( associated(p) )
    next=>p%next
    deallocate(p)
    p=>next
  end do
  if ( fast ) then
    if ( ldeg ) then
      o=0
      do i=1,n
        o=o+deg(i)
      end do
      deg=nint(deg*(dble(nc)/o))
    else
      deg=-1
    end if
  end if
  if ( cmax > 0 .and. cmax < n ) call randconf(iconf,deg,cmax)
  return
contains

subroutine fast_mode
  write(*,'(A)') '( Fast Mode ) ...'
  o=0
  se=sum(ne(2,:))
  nl=size(neqa)
  call prog%set(num=se)
  do l=1,nl
    allocate(occ(ne(2,l)))
    occ=.true.
    do i=1,ne(2,l)
      o=o+1
      if ( occ(i) .eqv. .false. ) cycle
      ic=i+ne(1,l)-1
      call int2conf(ic,aconf,lconf,E,mc,na,ns,ncm,nk,k,m)
      n=n+1
      if ( lpro .and. mod(n,100) == 0 ) call prog%put(o)
      allocate(p%next)
      p=>p%next
      allocate(p%iconf(ns))
      p%iconf=aconf
      p%deg=0
      do j=1,no
        oconf(1:k(1))=eqamat(j,aconf(1:k(1)))
        if ( all(oconf(1:k(1)) /= neqa(l)) ) cycle
        oconf(k(1)+1:)=eqamat(j,aconf(k(1)+1:))
        call sort(oconf(1:k(1)),1_2,k(1))
        call conf2int(ic,oconf,lconf,E,wc,na,ns,nc,nk,k)
        ic=ic-ne(1,l)+1
        if ( occ(ic) .eqv. .true. ) then
          p%deg=p%deg+1
          occ(ic)=.false.
        end if
      end do
    end do
    deallocate(occ)
  end do
  if ( lpro ) call prog%put(se)
  return
end subroutine fast_mode

subroutine std_mode
  write(*,'(A)') '( Standard Mode ) ...'
  allocate(occ(nc))
  occ=.true.
  call prog%set(num=nc)
  do i=1,nc
    if ( occ(i) .eqv. .false. ) cycle
    call int2conf(i,aconf,lconf,E,mc,na,ns,ncm,nk,k,m)
    n=n+1
    if ( lpro .and. mod(n,100) == 0 ) call prog%put(i)
    allocate(p%next)
    p=>p%next
    allocate(p%iconf(ns))
    p%iconf=aconf
    p%deg=0
    do j=1,no
      oconf=eqamat(j,aconf)
      call sort(oconf(1:k(1)),1_2,k(1))
      call conf2int(ic,oconf,lconf,E,wc,na,ns,nc,nk,k)
      if ( occ(ic) .eqv. .true. ) then
        p%deg=p%deg+1
        occ(ic)=.false.
      end if
    end do
  end do
  if ( lpro ) call prog%put(nc)
  return
end subroutine std_mode

end subroutine irrconfig

subroutine initial(lconf,nc,ncm,wc,mc,m,k)
  implicit none
  integer(8),allocatable :: mc(:)
  integer(1) :: i,nk
  integer(2) :: lconf(:),k(:),m(:),j,ibegin,iend
  integer(8) :: nc,ncm(:),wc(:)

  nk=size(k)-1
  m(1)=sum(k)
  do i=1,nk-1
    m(i+1)=m(i)-k(i)
  end do
  do i=1,nk
    ncm(i)=nchoosek(m(i),k(i))
  end do
  ibegin=1
  iend=0
  do i=1,nk
    iend=iend+k(i)
    lconf(ibegin:iend)=getlconf(m(i),k(i))
    ibegin=ibegin+k(i)
  end do
  nc=1
  do i=1,nk
    nc=nc*ncm(i)
  end do
  wc(1)=nc/ncm(1)
  do i=2,nk
    wc(i)=wc(i-1)/ncm(i)
  end do
  allocate(mc(nk-1))
  do i=1,nk-1
    mc(i)=1
    do j=i+1,nk
      mc(i)=mc(i)*ncm(j)
    end do
  end do
  return
end subroutine initial

subroutine matrixE(E,k,m)
  implicit none
  integer(8),allocatable :: E(:,:)
  integer(8),allocatable :: D(:,:)
  integer(2) :: k(:),i,j,m(:),ni,nj
  integer(1) :: nk

  nk=size(k)-1
  ni=maxval(m-k(1:nk))+1
  nj=maxval(k(1:nk))
  allocate(D(ni,nj))
  allocate(E(ni,nj))
  D(1,:)=0
  do i=2,ni
    do j=1,nj
      D(i,j)=nchoosek(nj-j+i-2_2,nj-j)
    end do
  end do
  do i=1,ni
    E(i,1)=sum(D(1:i,1))
  end do
  E(:,2:)=D(:,1:nj-1)
  return
end subroutine matrixE

subroutine neqatom(neqa,ne,eqamat,E,na,nc,wc,k,fast,ldeg)
  implicit none
  integer(2) :: eqamat(:,:)
  integer(8) :: E(:,:),wc,nc,m
  integer(8),allocatable :: ne(:,:)
  integer(2) :: no,na,i,j,ic,neqa_t(na),deg(na),k,n,lc,nj
  integer(2),allocatable :: neqa(:)
  logical(1) :: occ(na),fast,ldeg

  no=size(eqamat,1)
  occ=.true.
  n=0
  neqa_t=0
  deg=0
  do i=1,na
    if ( occ(i) .eqv. .false. ) cycle
    n=n+1
    neqa_t(n)=i
    do j=1,no
      ic=eqamat(j,i)
      if ( occ(ic) .eqv. .true. ) then
        deg(n)=deg(n)+1
        occ(ic)=.false.
      end if
    end do
  end do

  ldeg=.false.
  if ( n == 1 ) then
    fast=.true.
    ldeg=.true.
  elseif ( fast ) then
    do i=1,n-1
      if ( neqa_t(i+1) /= sum(deg(1:i))+1 ) then
        fast=.false.
        exit
      end if
    end do
  end if

  if ( fast ) then
    lc=na-k+1
    do i=1,n
      if ( neqa_t(i) > lc ) exit
    end do
    n=i-1
    allocate(neqa(n))
    neqa=neqa_t(1:n)

    allocate(ne(2,n))
    ne(1,1)=1
    m=0
    nj=size(E,2)-k
    do j=1,k
      m=m+E(lc-1,nj+j)
    end do
    ne(2,1)=nc-m

    do i=2,n
      m=0
      do j=1,k
        m=m+E(lc-neqa(i)+1,nj+j)
      end do
      ne(1,i)=nc-m
      m=0
      do j=1,k
        m=m+E(lc-neqa(i),nj+j)
      end do
      if ( lc == neqa(i) ) m=-1
      ne(2,i)=nc-m
    end do
    ne(2,:)=ne(2,:)-ne(1,:)
    ne(1,:)=(ne(1,:)-1)*wc+1
    ne(2,:)=ne(2,:)*wc
  end if
  return
end subroutine neqatom

subroutine int2conf(ic,aconf,lconf,E,mc,na,ns,ncm,nk,k,m)
  implicit none
  integer(2) :: aconf(:),lconf(:)
  integer(1) :: nk,i
  integer(2) :: m(:),k(:),ns,na,ibegin,iend
  integer(8) :: ic,mc(:),ik(nk),n,ncm(:),E(:,:)
  integer(2) :: label(na),cconf(ns)

  n=ic
  do i=1,nk-1
    ik(i)=ceiling(dble(n)/mc(i))
    n=n-(ik(i)-1)*mc(i)
  end do
  ik(nk)=n
  ibegin=1
  iend=0
  do i=1,nk
    iend=iend+k(i)
    cconf(ibegin:iend)=combconf(ik(i),ncm(i),E,m(i),k(i),lconf(ibegin:iend))
    ibegin=ibegin+k(i)
  end do
  aconf(1:k(1))=cconf(1:k(1))
  ibegin=k(1)+1
  iend=k(1)
  n=k(1)
  do i=2,nk
    iend=iend+k(i)
    label=complement(aconf(1:n),na)
    aconf(ibegin:iend)=label(cconf(ibegin:iend))
    ibegin=ibegin+k(i)
    n=n+k(i)
  end do
  return
end subroutine int2conf

subroutine conf2int(ic,aconf,lconf,E,wc,na,ns,nc,nk,k)
  implicit none
  integer(2) :: aconf(:),lconf(:)
  integer(1) :: nk,i
  integer(2) :: ns,na,k(:)
  integer(8) :: ic,nc,wc(:),m,E(:,:)
  integer(2) :: cconf(ns),conf(ns),j
  integer(2) :: n,nj,ibegin,iend

  cconf(1:k(1))=aconf(1:k(1))
  ibegin=k(1)+1
  iend=k(1)
  do i=2,nk
    iend=iend+k(i)
    cconf(ibegin:iend)=getcconf(aconf(1:iend),na,k(i),ibegin)
    ibegin=ibegin+k(i)
  end do
  conf=lconf-cconf+1
  ic=0
  ibegin=0
  nj=size(E,2)
  do i=1,nk
    m=0
    n=nj-k(i)
    do j=1,k(i)
      m=m+E(conf(ibegin+j),n+j)
    end do
    ic=ic+wc(i)*m
    ibegin=ibegin+k(i)
  end do
  ic=nc-ic
  return
end subroutine conf2int

function getcconf(a,m,n,k)
  implicit none
  integer(2) :: i,j,k,l,m,n
  integer(2) :: getcconf(n)
  integer(2) :: a(:),c(m)

  c=1
  c(a(1:k-1))=0
  c(a(k:))=2
  j=0
  l=0
  do i=1,maxval(a(k:))
    j=j+c(i)
    if ( c(i) == 2 ) then
      l=l+1
      j=j-1
      getcconf(l)=j
    end if
  end do
  return
end function getcconf

function getlconf(n,k)
  implicit none
  integer(2) :: i,n,k,getlconf(k)

  forall(i=1:k) getlconf(i)=n-k+i
  return
end function getlconf

function combconf(ic,nc,E,n,k,lconf)
  implicit none
  integer(8) :: ic,nc
  integer(8),target :: E(:,:)
  integer(2) :: n,k,ni,nj,i
  integer(2) :: combconf(k),lconf(:)

  ic=nc-ic
  ni=n-k+1
  nj=size(E,2)-k
  do i=1,k-1
    combconf(i)=binsearch(ic,E(1:ni,nj+i))
    ni=combconf(i)
    ic=ic-E(combconf(i),nj+i)
  end do
  combconf(i)=ic+1
  combconf=lconf-combconf+1
  return
end function combconf

subroutine randconf(iconf,deg,cmax)
  implicit none
  integer(2),allocatable :: iconf(:,:),deg(:)
  integer(2),allocatable :: iconf_t(:,:),deg_t(:)
  integer(8),allocatable :: a(:)
  integer(8) :: nc,ns,cmax

  nc=size(iconf,2)
  ns=size(iconf,1)
  allocate(iconf_t(ns,nc))
  allocate(deg_t(nc))
  iconf_t=iconf
  deg_t=deg
  deallocate(iconf)
  deallocate(deg)
  allocate(iconf(ns,cmax))
  allocate(deg(cmax))
  allocate(a(cmax))
  call random(a,nc)
  iconf=iconf_t(:,a)
  deg=deg_t(a)
  return
end subroutine randconf

end module configurations
