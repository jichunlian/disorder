module outfiles
  use functions
  implicit none
contains

subroutine output(eqamat,leqa,spgmat,lspg,iconf,deg,k,lcfg,a,x,atom,natom,symbols,site,lpos)
  implicit none
  integer(2) :: eqamat(:,:),iconf(:,:),deg(:)
  integer(2) :: k(:),natom(:)
  integer(1) :: site
  real(8) :: spgmat(:,:,:),a(:,:),x(:,:)
  logical(1) :: leqa,lspg,lcfg,lpos
  character(len=2) :: atom(:),symbols(:)

  write(*,'(A,$)') '  Writing output files :'
  call system('rm EQAMAT 2> /dev/null')
  if ( leqa ) then
    write(*,'(A,$)') '  EQAMAT'
    call outeqamat(eqamat)
  end if

  call system('rm SPGMAT 2> /dev/null')
  if ( lspg ) then
    write(*,'(A,$)') '  SPGMAT'
    call outspgmat(spgmat)
  end if

  call system('rm CFGMAT 2> /dev/null')
  if ( lcfg ) then
    write(*,'(A,$)') '  CFGMAT'
    call outconfig(iconf,deg,k)
  end if

  call system('rm -r poscar/ 2> /dev/null')
  if ( lpos ) then
    write(*,'(A,$)') '  POSCAR'
    call outposcar(a,x,atom,natom,symbols,k,site,iconf)
  end if
  write(*,'(/)')

  return
end subroutine output

subroutine outeqamat(eqamat)
  implicit none
  integer(2) :: eqamat(:,:)
  integer(2) :: no,na,i
  character(len=20) :: fm

  no=size(eqamat,1)
  na=size(eqamat,2)
  write(fm,*) na
  open(11,file='EQAMAT')
  do i=1,no
    write(11,'('//fm//'I4)') eqamat(i,:)
  end do
  close(11)
  return
end subroutine outeqamat

subroutine outspgmat(spgmat)
  implicit none
  real(8) :: spgmat(:,:,:)
  integer(2) :: nr,nt,no,i,j
  character(len=20) :: fm1,fm2

  nr=size(spgmat,3)
  nt=size(spgmat,2)-3
  no=nr*nt

  write(fm1,*) width(i_2=nr)
  write(fm2,*) width(i_2=nt)

  open(12,file='SPGMAT')
  write(12,'(I6,A3,I'//fm1//',A3,I'//fm2//')') no,' = ',nr,' x ',nt
  write(fm1,*) nt
  do i=1,nr
    write(12,*)
    do j=1,3
      write(12,'(3I4,'//fm1//'F12.6)') int(spgmat(j,1:3,i)),spgmat(j,4:,i)
    end do
  end do
  close(12)
  return
end subroutine outspgmat

subroutine outconfig(iconf,deg,k)
  implicit none
  integer(2) :: iconf(:,:),deg(:),k(:)
  integer(1) :: j,nk
  integer(8) :: i,nc
  integer(2) :: ibegin,iend
  character(len=20) :: fm

  nc=size(iconf,2)
  nk=size(k)
  open(13,file='CFGMAT')
  do i=1,nc
    write(13,'(I10,I6\)') i,deg(i)
    ibegin=1
    iend=0
    do j=1,nk-1
      iend=iend+k(j)
      write(fm,*) k(j)
      write(13,'(2X,'//fm//'I4\)') iconf(ibegin:iend,i)
      ibegin=ibegin+k(j)
    end do
    write(13,*)
  end do
  close(13)
  return
end subroutine outconfig

subroutine outposcar(a,x,atom,natom,symbols,k,site,iconf)
  use structure
  implicit none
  real(8) :: a(3,3),x(:,:)
  integer(2) :: natom(:),k(:),iconf(:,:)
  integer(1) :: site
  character(len=2) :: atom(:),symbols(:)
  real(8),allocatable :: x0(:,:),x1(:,:),x2(:,:),x_new(:,:)
  character(len=2),allocatable :: atom_new(:)
  integer(2),allocatable :: natom_new(:),nx(:),nx0(:),nx1(:),nx2(:),aconf(:)
  integer(1) :: ntype,ntype_new,nd,nk,nKw,j
  integer(2) :: na,na0,na1,na2,na_new,n
  integer(8) :: nc,i
  character(len=20) :: fm,tail,poscar

  na0=natom(site)
  allocate(nx0(na0))
  allocate(x0(3,na0))
  forall(i=1:na0) nx0(i)=i
  nx0=nx0+sum(natom(1:site-1))
  x0=x(:,nx0)

  ntype=size(natom)
  na=sum(natom)
  na1=na-na0
  allocate(nx1(na1))
  allocate(nx(na))
  allocate(x1(3,na1))
  nx=complement(nx0,na)
  nx1=nx(1:na1)
  x1=x(:,nx1)

  nk=size(k)
  ntype_new=ntype+nk-1
  if ( any(symbols == 'Kw') ) then
    ntype_new=ntype_new-1
  end if
  allocate(natom_new(ntype_new))
  allocate(atom_new(ntype_new))
  j=0
  do i=1,ntype
    if ( i /= site ) then
      j=j+1
      natom_new(j)=natom(i)
      atom_new(j)=atom(i)
    end if
  end do
  nKw=0
  do i=1,nk
    if ( symbols(i) /= 'Kw' ) then
      j=j+1
      natom_new(j)=k(i)
      atom_new(j)=symbols(i)
    else
      nKw=i
    end if
  end do

  na_new=sum(natom_new)
  allocate(x_new(3,na_new))
  na2=na_new-na1
  allocate(x2(3,na2))
  allocate(nx2(na2))

  nc=size(iconf,2)
  nd=width(i_8=nc)
  write(fm,*) nd
  call system('mkdir poscar 2> /dev/null')
  na=sum(k)
  allocate(aconf(na))
  do i=1,nc
    if ( nKw == nk ) then
      nx2=iconf(i,:)
    else
      aconf=complement(iconf(:,i),na)
      n=na-k(nk)
      aconf(n+1:)=aconf(1:k(nk))
      aconf(1:n)=iconf(:,i)
      if ( nkW /= 0 ) then
        n=sum(k(1:nKw))
        aconf(n-k(nKw)+1:na2)=aconf(n+1:)
      end if
      nx2=aconf(1:na2)
    end if
    x2=x0(:,nx2)
    x_new(:,1:na1)=x1
    x_new(:,na1+1:)=x2
    write(tail,'(I'//fm//'.'//fm//')') i
    tail=trim(tail)
    poscar='poscar/POSCAR-'//tail//''
    call writepos(a,x_new,atom_new,natom_new,poscar)
  end do
  return
end subroutine outposcar

end module outfiles
