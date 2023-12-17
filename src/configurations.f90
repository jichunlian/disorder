module configurations
  use functions
  implicit none
contains

subroutine irrconfig(iconf,deg,eqamat,group,k,lpro)
  use progress
  use stdout
  implicit none
  type :: datalink
    integer(2),allocatable :: iconf(:)
    integer(2) :: deg
    type(datalink),pointer :: next
  end type datalink
  integer(2) :: eqamat(:,:),group(:),k(:)
  integer(2),allocatable :: aconf(:),oconf(:),lconf(:),iconf(:,:)
  integer(2),allocatable :: deg(:),degm(:),ib(:),ie(:),m(:),kg(:)
  integer(8),allocatable :: E(:,:),ncm(:),nct(:),ncp(:),il(:)
  integer(1) :: nk
  integer(2) :: na,ns,no,ng,g,o
  integer(4) :: nirr,i
  integer(8) :: i2,i3,i4,i5,nc,np,ic
  logical(1),allocatable :: occ2(:),occ3(:),occ4(:)
  logical(1),allocatable :: occ5(:),valid(:,:)
  logical(1) :: lpro
  type(cls_prog) :: prog
  type(datalink),pointer :: p,head,next

  call preprocessing
  call eliminating
  call postprocessing
  return
contains

subroutine preprocessing
  nk=size(k)-1
  na=sum(k)
  ns=na-k(nk+1)
  no=size(eqamat,2)
  ng=size(group)-1
  allocate(ncm(nk))
  allocate(nct(ng))
  allocate(kg(ng))
  allocate(ncp(nk-1))
  allocate(degm(nk))
  allocate(m(nk))
  allocate(ib(nk))
  allocate(ie(nk))
  allocate(aconf(ns))
  allocate(oconf(ns))
  allocate(il(nk))
  allocate(lconf(ns))
  allocate(valid(no,nk))
  call confnum(nc,ncm,nct,ncp,ib,ie,m,k,group,ng,kg)
  call stdout_3(nc)
  call matrixE(E,m,k)
  call lastconf(lconf,m,k)
  if ( nk > 1 ) allocate(occ3(ncm(2)))
  if ( nk > 2 ) allocate(occ4(ncm(3)))
  if ( nk > 3 ) allocate(occ5(ncm(4)))
  allocate(head)
  nullify(head%next)
  p=>head
  np=0
  do g=1,ng
    np=np+nct(g)
  end do
  do g=2,nk
    np=np*ncm(g)
  end do
  if ( lpro ) call prog%set(num=np)
  nirr=0
  return
end subroutine preprocessing

subroutine eliminating
  do g=1,ng
    allocate(occ2(nct(g)))
    occ2=.true.
    do i2=1,nct(g)
      if ( occ2(i2) .eqv. .false. ) cycle
      call binary
      if ( nk == 1 ) then
        il=i2
        call trialclose(degm(1))
        cycle
      end if
      occ3=.true.
      do i3=1,ncm(2)
        if ( occ3(i3) .eqv. .false. ) cycle
        call multinary(2_1,i3,occ3,degm(2),ib(2),ie(2))
        degm(2)=degm(1)*degm(2)
        if ( nk == 2 ) then
          il=(/ i2, i3 /)
          call trialclose(degm(2))
          cycle
        end if
        occ4=.true.
        do i4=1,ncm(3)
          if ( occ4(i4) .eqv. .false. ) cycle
          call multinary(3_1,i4,occ4,degm(3),ib(3),ie(3))
          degm(3)=degm(2)*degm(3)
          if ( nk == 3 ) then
            il=(/ i2, i3, i4 /)
            call trialclose(degm(3))
            cycle
          end if
          occ5=.true.
          do i5=1,ncm(4)
            if ( occ5(i5) .eqv. .false. ) cycle
            call multinary(4_1,i5,occ5,degm(4),ib(4),ie(4))
            degm(4)=degm(3)*degm(4)
            if ( nk == 4 ) then
              il=(/ i2, i3, i4, i5 /)
              call trialclose(degm(4))
              cycle
            end if
          end do
        end do
      end do
    end do
    deallocate(occ2)
  end do
  return
end subroutine eliminating

subroutine postprocessing
  if ( lpro ) call prog%put(np)
  call stdout_4(nirr)
  allocate(iconf(ns,nirr))
  allocate(deg(nirr))
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
  call degeneracy_correction
  return
end subroutine postprocessing

subroutine binary
  valid(:,1)=.false.
  aconf(1)=0
  call int2cconf(i2,aconf(2:k(1)),nct(g),E,m(1)-group(g),k(1)-1_2,lconf(2:k(1))-group(g))
  aconf(1:k(1))=aconf(1:k(1))+group(g)
  degm(1)=0
  do o=1,no
    oconf(1:k(1))=eqamat(aconf(1:k(1)),o)
    if ( all(oconf(1:k(1)) /= group(g)) ) cycle
    call sort(oconf(1:k(1)),1_2,k(1))
    call cconf2int(ic,oconf(2:k(1))-group(g),lconf(2:k(1))-group(g),nct(g),E,k(1)-1_2)
    if ( ic == i2 ) valid(o,1)=.true.
    if ( occ2(ic) .eqv. .true. ) then
      degm(1)=degm(1)+1
      occ2(ic)=.false.
    end if
  end do
  return
end subroutine binary

subroutine multinary(n,im,occ,deg,ib,ie)
  integer(1) :: n
  integer(8) :: im
  logical(1) :: occ(:)
  integer(2) :: deg,ib,ie

  valid(:,n)=.false.
  call int2cconf(im,aconf(ib:ie),ncm(n),E,m(n),k(n),lconf(ib:ie))
  call cconf2aconf(aconf(ib:ie),aconf(1:ib-1),na)
  deg=0
  do o=1,no
    if ( valid(o,n-1) .eqv. .false. ) cycle
    oconf(ib:ie)=eqamat(aconf(ib:ie),o)
    call sort(oconf(ib:ie),1_2,k(n))
    call aconf2cconf(oconf(ib:ie),aconf(1:ib-1),na)
    call cconf2int(ic,oconf(ib:ie),lconf(ib:ie),ncm(n),E,k(n))
    if ( ic == im ) valid(o,n)=.true.
    if ( occ(ic) .eqv. .true. ) then
      deg=deg+1
      occ(ic)=.false.
    end if
  end do
  return
end subroutine multinary

subroutine trialclose(deg)
  integer(2) :: deg,i

  nirr=nirr+1
  if ( lpro .and. mod(nirr,100) == 0 ) then
    do i=2,g
      il(1)=il(1)+nct(i-1)
    end do
    ic=il(nk)
    do i=1,nk-1
      ic=ic+(il(i)-1)*ncp(i)
    end do
    call prog%put(ic)
  end if
  allocate(p%next)
  p=>p%next
  allocate(p%iconf(ns))
  p%iconf=aconf
  p%deg=deg
  return
end subroutine trialclose

subroutine degeneracy_correction
  integer(2) :: k1

  if ( group(2) == na+1 ) then
    deg=nint(deg*(dble(na)/k(1)))
  else
    do i=1,nirr
      do g=1,ng
        if ( iconf(1,i) == group(g) ) exit
      end do
      k1=binsearch_2(group(g+1)-1_2,iconf(1:k(1),i))
      deg(i)=nint(deg(i)*(dble(kg(g))/k1))
    end do
  end if
  return
end subroutine degeneracy_correction

end subroutine irrconfig

subroutine confnum(nc,ncm,nct,ncp,ib,ie,m,k,group,ng,kg)
  implicit none
  integer(2) :: i,j,nk,ng,kg(:),group(:)
  integer(2) :: k(:),m(:),ib(:),ie(:)
  integer(8) :: nc,ncm(:),nct(:),ncp(:)

  nk=size(k)-1
  m(1)=sum(k)
  do i=1,nk-1
    m(i+1)=m(i)-k(i)
  end do
  do i=1,nk
    ncm(i)=nchoosek(m(i),k(i))
  end do
  do i=1,ng
    kg(i)=group(i+1)-group(i)
  end do
  j=m(1)-1
  do i=1,ng
    nct(i)=nchoosek(j,k(1)-1_2)
    j=j-kg(i)
  end do
  do i=1,nk-1
    ncp(i)=1
    do j=i+1,nk
      ncp(i)=ncp(i)*ncm(j)
    end do
  end do
  nc=1
  do i=1,nk
    nc=nc*ncm(i)
  end do
  ib(1)=1
  ie(1)=k(1)
  do i=1,nk-1
    ib(i+1)=ib(i)+k(i)
    ie(i+1)=ie(i)+k(i+1)
  end do
  return
end subroutine confnum

subroutine matrixE(E,m,k)
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

subroutine lastconf(lconf,m,k)
  implicit none
  integer(2) :: lconf(:),m(:),k(:)
  integer(2) :: ib,j
  integer(1) :: i,nk

  nk=size(k)-1
  ib=0
  do i=1,nk
    forall(j=1:k(i)) lconf(ib+j)=m(i)-k(i)+j
    ib=ib+k(i)
  end do
end subroutine lastconf

subroutine int2cconf(ic,cconf,nc,E,m,k,lconf)
  implicit none
  integer(8) :: ic,nc,n
  integer(8),target :: E(:,:)
  integer(2) :: m,k,ni,nj,i
  integer(2) :: cconf(k),lconf(:)

  n=nc-ic
  ni=m-k+1
  nj=size(E,2)-k
  do i=1,k-1
    cconf(i)=binsearch_8(n,E(1:ni,nj+i))
    ni=cconf(i)
    n=n-E(cconf(i),nj+i)
  end do
  cconf(i)=n+1
  cconf=lconf-cconf+1
  return
end subroutine int2cconf

subroutine cconf2aconf(conf,fconf,na)
  implicit none
  integer(2) :: conf(:),fconf(:)
  integer(2) :: na
  integer(2) :: label(na)

  label=complement(fconf,na)
  conf=label(conf)
  return
end subroutine cconf2aconf

subroutine aconf2cconf(conf,fconf,na)
  implicit none
  integer(2) :: i,j,l,na
  integer(2) :: conf(:),fconf(:),c(na)

  c=1
  c(fconf)=0
  c(conf)=2
  j=0
  l=0
  do i=1,maxval(conf)
    j=j+c(i)
    if ( c(i) == 2 ) then
      l=l+1
      j=j-1
      conf(l)=j
    end if
  end do
  return
end subroutine aconf2cconf

subroutine cconf2int(ic,cconf,lconf,nc,E,k)
  implicit none
  integer(2) :: cconf(:),lconf(:),k
  integer(8) :: ic,nc,E(:,:)
  integer(2) :: n,j,conf(k)

  ic=0
  n=size(E,2)-k
  conf=lconf-cconf+1
  do j=1,k
    ic=ic+E(conf(j),n+j)
  end do
  ic=nc-ic
  return
end subroutine cconf2int

end module configurations
