module symmetry
  use stdout
  use functions
  implicit none
contains

subroutine eqamatrix(eqamat,spgmat,a,x,natom,symprec,site)
  implicit none
  integer(2),allocatable :: eqamat(:,:)
  real(8) :: a(3,3),x(:,:),symprec
  integer(2) :: natom(:),satom
  integer(1) :: site
  integer(2) :: i,j,n,nr,nt,no,na,s
  integer(2),allocatable :: atom(:)
  integer(1),allocatable :: rot(:,:,:)
  integer(2),allocatable :: eqamat1(:,:)
  integer(2),allocatable :: eqamat2(:,:)
  integer(2),allocatable :: eqamat3(:,:)
  real(8),allocatable :: tra1(:,:)
  real(8),allocatable :: tra2(:,:)
  real(8),allocatable :: spgmat(:,:,:)
  character(len=20) :: cls,sys

  call rotation(eqamat1,rot,tra1,a,x,natom,symprec)
  call translation(eqamat2,tra2,a,x,natom,symprec)
  call spgmatrix(spgmat,rot,tra1,tra2)
  call pointgroup(rot,cls,sys)
  nr=size(eqamat1,1)
  nt=size(eqamat2,1)
  satom=sum(natom)
  no=nr*nt
  call stdout_2(no,nr,nt,cls,sys)
  allocate(eqamat3(no,satom))
  n=0
  do i=1,nr
    do j=1,nt
      n=n+1
      eqamat3(n,:)=eqamat2(j,eqamat1(i,:))
    end do
  end do
  na=natom(site)
  allocate(atom(na))
  s=sum(natom(1:site-1))
  forall (i=1:na)
    atom(i)=s+i
  end forall
  allocate(eqamat(no,na))
  eqamat=eqamat3(:,atom)-s
  return
end subroutine eqamatrix

subroutine spgmatrix(spgmat,rot,tra1,tra2)
  implicit none
  real(8),allocatable :: spgmat(:,:,:)
  real(8) :: tra1(:,:),tra2(:,:)
  real(8),allocatable :: tra(:,:)
  integer(2) :: nr,nt,i,j
  integer(1) :: rot(:,:,:)

  nr=size(tra1,2)
  nt=size(tra2,2)
  allocate(tra(3,nt))
  allocate(spgmat(3,nt+3,nr))
  do i=1,nr
    do j=1,nt
      tra(:,j)=tra1(:,i)+tra2(:,j)
    end do
    spgmat(:,1:3,i)=dble(rot(:,:,i))
    spgmat(:,4:,i)=tra
  end do
  return
end subroutine spgmatrix

subroutine rotation(eqamat,rot,tra,a,x,natom,symprec)
  implicit none
  integer(2),allocatable :: eqamat(:,:)
  real(8) :: a(3,3),x(:,:),symprec
  integer(2) :: natom(:),satom
  real(8) :: G(3,3),dG(3,3),x_t(3),T(3),dx(3),c(3)
  real(8),allocatable :: tra(:,:),tra_t(:,:)
  integer(1) :: W11,W12,W13,W21,W22,W23,W31,W32,W33
  integer(1) :: W(3,3),rot_t(3,3,48)
  integer(2) :: i,j,k,m,n,m1,m2,n1,n2,col,num
  integer(1),allocatable :: rot(:,:,:)
  logical(1) :: discard

  G=matmul(transpose(a),a)
  n=0
  do W11=-1,1
    do W12=-1,1
      do W13=-1,1
        do W21=-1,1
          do W22=-1,1
            do W23=-1,1
              do W31=-1,1
                do W32=-1,1
                  do W33=-1,1
                    W(1,:)=(/ W11, W12, W13 /)
                    W(2,:)=(/ W21, W22, W23 /)
                    W(3,:)=(/ W31, W32, W33 /)
                    if ( abs(det(W)) == 1 ) then
                      dG=dabs(G-matmul(matmul(transpose(W),G),W))
                      if ( all(dG < symprec) ) then
                        n=n+1
                        rot_t(:,:,n)=W
                      end if
                    end if
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end do

  allocate(tra_t(3,n))
  m=minloc(natom,1)
  n1=natom(m)
  m2=sum(natom(1:m))
  m1=m2-n1+1
  satom=sum(natom)
  num=0
  do i=1,n
    W=rot_t(:,:,i)
    do j=m1,m2
      T=x(:,j)-matmul(W,x(:,m1))
      k=1
      n1=1
      n2=natom(1)
      do m=1,satom
        x_t=matmul(W,x(:,m))+T
        discard=.true.
        do n=n1,n2
          dx=x_t-x(:,n)
          c=dabs(matmul(a,(dx-nint(dx))))
          if ( all(c < symprec) ) then
            discard=.false.
            exit
          end if
        end do
        if ( discard == .true. ) then
          exit
        end if        
        if ( m == sum(natom(1:k)) .and. m /= satom ) then
          k=k+1
          n1=n2+1
          n2=n2+natom(k)
        end if
      end do
      if ( discard == .false. ) then
        num=num+1
        rot_t(:,:,num)=W
        tra_t(:,num)=T
        exit
      end if
    end do
  end do

  allocate(rot(3,3,num))
  allocate(tra(3,num))
  rot=rot_t(:,:,1:num)
  tra=tra_t(:,1:num)
  allocate(eqamat(num,satom))
  do i=1,num
    W=rot(:,:,i)
    T=tra(:,i)
    col=1
    do j=1,satom
      x_t=matmul(W,x(:,j))+T
      do k=1,satom
        dx=x_t-x(:,k)
        c=dabs(matmul(a,(dx-nint(dx))))
        if ( all(c < symprec) ) then
          eqamat(i,col)=k
          col=col+1
          exit
        end if
      end do
    end do
  end do
  return
end subroutine rotation

subroutine translation(eqamat,tra,a,x,natom,symprec)
  integer(2),allocatable :: eqamat(:,:)
  real(8) :: a(3,3),x(:,:),symprec
  integer(2) :: natom(:),satom
  real(8) :: x_t(3),T(3),dx(3),c(3)
  real(8),allocatable :: tra(:,:),tra_t(:,:)
  integer(2) :: m,m1,m2,n,n1,n2,i,j,k,num,col
  logical(1) :: discard

  m=minloc(natom,1)
  n=natom(m)
  allocate(tra_t(3,n))
  m2=sum(natom(1:m))
  m1=m2-n+1
  num=0
  satom=sum(natom)
  do i=m1,m2
    T=x(:,i)-x(:,m1)
    k=1
    n1=1
    n2=natom(1)
    do m=1,satom
      x_t=x(:,m)+T
      discard=.true.
      do n=n1,n2
        dx=x_t-x(:,n)
        c=dabs(matmul(a,(dx-nint(dx))))
        if ( all(c < symprec) ) then
          discard=.false.
          exit
        end if
      end do
      if ( discard == .true. ) then
        exit
      end if        
      if ( m == sum(natom(1:k)) .and. m /= satom ) then
        k=k+1
        n1=n2+1
        n2=n2+natom(k)
      end if
    end do
    if ( discard == .false. ) then
      num=num+1
      tra_t(:,num)=T
    end if
  end do

  allocate(tra(3,num))
  tra=tra_t(:,1:num)
  allocate(eqamat(num,satom))
  do i=1,num
    T=tra(:,i)
    col=1
    do j=1,satom
      x_t=x(:,j)+T
      do k=1,satom
        dx=x_t-x(:,k)
        c=dabs(matmul(a,(dx-nint(dx))))
        if ( all(c < symprec) ) then
          eqamat(i,col)=k
          col=col+1
          exit
        end if
      end do
    end do
  end do
  return
end subroutine translation

subroutine pointgroup(rot,cls,sys)
  implicit none
  integer(1) :: rot(:,:,:)
  integer(1) :: rottyp(10),nrt(32,10)
  integer(1) :: trW,detW
  integer(2) :: i,n
  character(len=20) :: cls,sys 

  n=size(rot,3)
  rottyp=0
  do i=1,n
    trW=trace(rot(:,:,i))
    detW=det(rot(:,:,i))
    if     ( trW == -2 .and. detW == -1 ) then
        rottyp(1)=rottyp(1)+1
    elseif ( trW == -1 .and. detW == -1 ) then
        rottyp(2)=rottyp(2)+1
    elseif ( trW ==  0 .and. detW == -1 ) then
        rottyp(3)=rottyp(3)+1
    elseif ( trW ==  1 .and. detW == -1 ) then
        rottyp(4)=rottyp(4)+1
    elseif ( trW == -3 .and. detW == -1 ) then
        rottyp(5)=rottyp(5)+1
    elseif ( trW ==  3 .and. detW ==  1 ) then
        rottyp(6)=rottyp(6)+1
    elseif ( trW == -1 .and. detW ==  1 ) then
        rottyp(7)=rottyp(7)+1
    elseif ( trW ==  0 .and. detW ==  1 ) then
        rottyp(8)=rottyp(8)+1
    elseif ( trW ==  1 .and. detW ==  1 ) then
        rottyp(9)=rottyp(9)+1
    elseif ( trW ==  2 .and. detW ==  1 ) then
        rottyp(10)=rottyp(10)+1
    end if
  end do

   nrt(1,:)=(/ 0, 0, 0, 0, 0, 1, 0, 0, 0, 0 /)
   nrt(2,:)=(/ 0, 0, 0, 0, 1, 1, 0, 0, 0, 0 /)
   nrt(3,:)=(/ 0, 0, 0, 0, 0, 1, 1, 0, 0, 0 /)
   nrt(4,:)=(/ 0, 0, 0, 1, 0, 1, 0, 0, 0, 0 /)
   nrt(5,:)=(/ 0, 0, 0, 1, 1, 1, 1, 0, 0, 0 /)
   nrt(6,:)=(/ 0, 0, 0, 0, 0, 1, 3, 0, 0, 0 /)
   nrt(7,:)=(/ 0, 0, 0, 2, 0, 1, 1, 0, 0, 0 /)
   nrt(8,:)=(/ 0, 0, 0, 3, 1, 1, 3, 0, 0, 0 /)
   nrt(9,:)=(/ 0, 0, 0, 0, 0, 1, 1, 0, 2, 0 /)
  nrt(10,:)=(/ 0, 2, 0, 0, 0, 1, 1, 0, 0, 0 /)
  nrt(11,:)=(/ 0, 2, 0, 1, 1, 1, 1, 0, 2, 0 /)
  nrt(12,:)=(/ 0, 0, 0, 0, 0, 1, 5, 0, 2, 0 /)
  nrt(13,:)=(/ 0, 0, 0, 4, 0, 1, 1, 0, 2, 0 /)
  nrt(14,:)=(/ 0, 2, 0, 2, 0, 1, 3, 0, 0, 0 /)
  nrt(15,:)=(/ 0, 2, 0, 5, 1, 1, 5, 0, 2, 0 /)
  nrt(16,:)=(/ 0, 0, 0, 0, 0, 1, 0, 2, 0, 0 /)
  nrt(17,:)=(/ 0, 0, 2, 0, 1, 1, 0, 2, 0, 0 /)
  nrt(18,:)=(/ 0, 0, 0, 0, 0, 1, 3, 2, 0, 0 /)
  nrt(19,:)=(/ 0, 0, 0, 3, 0, 1, 0, 2, 0, 0 /)
  nrt(20,:)=(/ 0, 0, 2, 3, 1, 1, 3, 2, 0, 0 /)
  nrt(21,:)=(/ 0, 0, 0, 0, 0, 1, 1, 2, 0, 2 /)
  nrt(22,:)=(/ 2, 0, 0, 1, 0, 1, 0, 2, 0, 0 /)
  nrt(23,:)=(/ 2, 0, 2, 1, 1, 1, 1, 2, 0, 2 /)
  nrt(24,:)=(/ 0, 0, 0, 0, 0, 1, 7, 2, 0, 2 /)
  nrt(25,:)=(/ 0, 0, 0, 6, 0, 1, 1, 2, 0, 2 /)
  nrt(26,:)=(/ 2, 0, 0, 4, 0, 1, 3, 2, 0, 0 /)
  nrt(27,:)=(/ 2, 0, 2, 7, 1, 1, 7, 2, 0, 2 /)
  nrt(28,:)=(/ 0, 0, 0, 0, 0, 1, 3, 8, 0, 0 /)
  nrt(29,:)=(/ 0, 0, 8, 3, 1, 1, 3, 8, 0, 0 /)
  nrt(30,:)=(/ 0, 0, 0, 0, 0, 1, 9, 8, 6, 0 /)
  nrt(31,:)=(/ 0, 6, 0, 6, 0, 1, 3, 8, 0, 0 /)
  nrt(32,:)=(/ 0, 6, 8, 9, 1, 1, 9, 8, 6, 0 /)

  do i=1,32
    if ( rottyp(1)  == nrt(i,1) .and. &
         rottyp(2)  == nrt(i,2) .and. &
         rottyp(3)  == nrt(i,3) .and. &
         rottyp(4)  == nrt(i,4) .and. &
         rottyp(5)  == nrt(i,5) .and. &
         rottyp(6)  == nrt(i,6) .and. &
         rottyp(7)  == nrt(i,7) .and. &
         rottyp(8)  == nrt(i,8) .and. &
         rottyp(9)  == nrt(i,9) .and. &
         rottyp(10) == nrt(i,10) ) then
      exit
    end if
  end do

  select case (i)
    case (1)
      sys='Triclinic'
      cls='C_1 (1)'
    case (2)
      sys='Triclinic'
      cls='C_i = S_2 (-1)'
    case (3)
      sys='Monoclinic'
      cls='C_2 (2)'
    case (4)
      sys='Monoclinic'
      cls='C_1h = C_s (m)'
    case (5)
      sys='Monoclinic'
      cls='C_2h (2/m)'
    case (6)
      sys='Orthorhombic'
      cls='D_2 = V (222)'
    case (7)
      sys='Orthorhombic'
      cls='C_2v (mm2)'
    case (8)
      sys='Orthorhombic'
      cls='D_2h = V_h (mmm)'
    case (9)
      sys='Tetragonal'
      cls='C_4 (4)'
    case (10)
      sys='Tetragonal'
      cls='S_4 (-4)'
    case (11)
      sys='Tetragonal'
      cls='C_4h (4/m)'
    case (12)
      sys='Tetragonal'
      cls='D_4 (422)'
    case (13)
      sys='Tetragonal'
      cls='C_4v (4mm)'
    case (14)
      sys='Tetragonal'
      cls='D_2d = V_d (-42m)'
    case (15)
      sys='Tetragonal'
      cls='D_4h (4/mmm)'
    case (16)
      sys='Trigonal'
      cls='C_3 (3)'
    case (17)
      sys='Trigonal'
      cls='C_3i = S_6 (-3)'
    case (18)
      sys='Trigonal'
      cls='D_3 (32)'
    case (19)
      sys='Trigonal'
      cls='C_3v (3m)'
    case (20)
      sys='Trigonal'
      cls='D_3d (-3m)'
    case (21)
      sys='Hexagonal'
      cls='C_6 (6)'
    case (22)
      sys='Hexagonal'
      cls='C_3h (-6)'
    case (23)
      sys='Hexagonal'
      cls='C_6h (6/m)'
    case (24)
      sys='Hexagonal'
      cls='D_6 (622)'
    case (25)
      sys='Hexagonal'
      cls='C_6v (6mm)'
    case (26)
      sys='Hexagonal'
      cls='D_3h (-6m2)'
    case (27)
      sys='Hexagonal'
      cls='D_6h (6/mmm)'
    case (28)
      sys='Cubic'
      cls='T (23)'
    case (29)
      sys='Cubic'
      cls='T_h (m-3)'
    case (30)
      sys='Cubic'
      cls='O (432)'
    case (31)
      sys='Cubic'
      cls='T_d (-43m)'
    case (32)
      sys='Cubic'
      cls='O_h (m-3m)'
    case default
      call stderr("ERROR : Can't Find Point Group !")
  end select
  return
end subroutine pointgroup

end module symmetry 
