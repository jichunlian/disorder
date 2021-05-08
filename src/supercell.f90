program supercell
  use structure
  use functions
  implicit none
  real(8) :: a(3,3),la,lb,lc,sina,cosa,cosb,cosr
  real(8),allocatable :: x(:,:),xs(:,:),T(:,:),Pd(:,:)
  integer(1) :: M(3,3),bd(2,3),vertex(3,8)
  integer(2),allocatable :: natom(:)
  integer(1),allocatable :: Pc(:,:)
  character(len=2),allocatable :: atom(:)
  integer(2) :: i,j,k,n,nt,np,na

  read(*,*) M(1,:)
  read(*,*) M(2,:)
  read(*,*) M(3,:)
  M=transpose(M) 
  nt=abs(det(M))
  allocate(T(3,nt))
  vertex(:,1)=0
  vertex(:,2:4)=M
  vertex(:,5)=M(:,1)+M(:,2)
  vertex(:,6)=M(:,1)+M(:,3)
  vertex(:,7)=M(:,2)+M(:,3)
  vertex(:,8)=M(:,1)+M(:,2)+M(:,3)
  do i=1,3
    bd(1,i)=maxval(vertex(i,:))
    bd(2,i)=minval(vertex(i,:))
  end do
  np=(bd(1,1)-bd(2,1))*(bd(1,2)-bd(2,2))*(bd(1,3)-bd(2,3))
  allocate(Pc(3,np))
  allocate(Pd(3,np))
  n=0
  do i=bd(2,1),bd(1,1)-1
    do j=bd(2,2),bd(1,2)-1
      do k=bd(2,3),bd(1,3)-1
        n=n+1
        Pc(:,n)=(/ i, j, k /)
      end do
    end do
  end do
  Pd=matmul(inverse(dble(M)),dble(Pc))
  n=0
  do i=1,np
    if ( all(Pd(:,i) > -0.01) .and. all(Pd(:,i) < 0.99) ) then
      n=n+1
      T(:,n)=Pd(:,i)
    end if
  end do
  call readpos(a,x,atom,natom,'POSCAR')
  x=matmul(a,x)
  a=matmul(a,dble(M))
  x=matmul(inverse(a),x)
  na=sum(natom)
  natom=natom*nt
  allocate(xs(3,na*nt))
  n=0
  do i=1,na
    do j=1,nt
      n=n+1
      xs(:,n)=x(:,i)+T(:,j)
    end do
  end do
  xs=xs-floor(xs)
  la=norm(a(:,1))
  lb=norm(a(:,2))
  lc=norm(a(:,3))
  sina=norm(cross(a(:,1),a(:,2)))/(la*lb)
  cosa=dot_product(a(:,1),a(:,2))/(la*lb)
  cosb=dot_product(a(:,1),a(:,3))/(la*lc)
  cosr=dot_product(a(:,2),a(:,3))/(lb*lc)
  a(1,1)=la
  a(2,1)=0
  a(3,1)=0
  a(1,2)=lb*cosa
  a(2,2)=lb*sina
  a(3,2)=0
  a(1,3)=lc*cosr
  a(2,3)=(lc*cosb-a(1,3)*cosa)/sina
  a(3,3)=dsqrt(lc**2-a(1,3)**2-a(2,3)**2)
  call writepos(a,xs,atom,natom,'SPOSCAR','SPOSCAR') 
end program supercell
