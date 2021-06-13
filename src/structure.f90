module structure
  implicit none
contains

subroutine readpos(a,x,symbol,natom,poscar)
  use functions
  implicit none
  real(8) :: a(3,3)
  real(8),allocatable :: x(:,:)
  integer(1) :: ntype,io_error
  integer(2) :: satom,i
  integer(2),allocatable :: natom(:)
  real(8) :: scal
  character(len=2),allocatable :: symbol(:)
  character(len=1) :: coordinate
  character(len=*) :: poscar

  open(10,file=poscar)
  read(10,*)
  read(10,*) scal
  read(10,*) a
  a=scal*a
  read(10,*)
  ntype=1
  io_error=0
  do while (io_error == 0)
    allocate(natom(ntype))
    read(10,*,iostat=io_error) natom
    ntype=ntype+1
    deallocate(natom)
    backspace(10)
  end do
  ntype=ntype-2
  allocate(symbol(ntype))
  allocate(natom(ntype))
  backspace(10)
  backspace(10)
  read(10,*) symbol
  read(10,*) natom
  satom=sum(natom)
  allocate(x(3,satom))
  read(10,*) coordinate
  if ( coordinate == 'S' ) then
    read(10,*) coordinate
  end if
  do i=1,satom
    read(10,*) x(:,i)
  end do
  close(10)
  if ( coordinate == 'C' ) then
    x=matmul(inverse(a),x)
  end if
  return
end subroutine readpos

subroutine writepos(a,x,symbol,natom,poscar)
  implicit none
  real(8) :: a(3,3),x(:,:)
  integer(2) :: natom(:),satom,i
  integer(1) :: ntype
  character(len=2) :: symbol(:),fm
  character(len=*) :: poscar

  open(11,file=poscar)
  write(11,'(A)') trim(poscar)
  write(11,'(A)') '   1.00000000000000'
  write(11,'(1X,3F22.16)') a(:,1)
  write(11,'(1X,3F22.16)') a(:,2)
  write(11,'(1X,3F22.16)') a(:,3)
  ntype=size(natom)
  write(fm,'(I2)') ntype
  write(11,'('//fm//'A6)') symbol
  write(11,'('//fm//'I6)') natom
  write(11,'(A6)') 'Direct'
  satom=sum(natom)
  do i=1,satom
    write(11,'(3F20.16)') x(:,i)
  end do
  close(11)
  return
end subroutine writepos

end module structure
