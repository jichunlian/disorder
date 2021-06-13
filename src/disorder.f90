program disorder
  use structure
  use symmetry
  use groups
  use configurations
  use outfiles
  use stdout
  implicit none
  real(8) :: a(3,3),prec
  real(8),allocatable :: x(:,:)
  integer(2),allocatable :: natom(:),k(:),subs(:)
  integer(2),allocatable :: eqamat(:,:),group(:)
  integer(2),allocatable :: iconf(:,:),deg(:)
  real(8),allocatable :: spgmat(:,:,:)
  integer(1) :: site,nsub,io_err,i
  real(4) :: time0,time1
  logical(1) :: alive,leqa,lspg,lcfg,lpos,lpro
  character(len=2),allocatable :: symbols(:),symb(:),atom(:)
  namelist /input/ nsub,subs,symb,prec,site,leqa,lspg,lcfg,lpos,lpro

  call cpu_time(time0)
  call stdout_0

! Setting default value
  nsub=2
  site=1
  prec=1D-5
  lcfg=.true.
  lspg=.false.
  leqa=.false.
  lpos=.false.
  lpro=.false.
! End setting

  allocate(subs(5))
  allocate(symb(5))
  inquire(file='INDSOD',exist=alive)
  if ( alive .eqv. .false. ) call stderr('The INDSOD file does not exist !')
  open(10,file='INDSOD')
  read(10,nml=input,iostat=io_err)
  close(10)
  if ( io_err /= 0 ) call stderr('The error occurred while reading the INDSOD file !')
  if ( nsub > 5 .or. nsub < 2 ) call stderr('NSUB needs to be between 2 and 5 !')
  if ( prec > 1D-2 ) call stderr('PREC is recommended no greater than 1D-2 !')
  allocate(k(nsub))
  allocate(symbols(nsub))
  k=subs(1:nsub)
  if ( any(k < 1) ) call stderr('The elements of SUBS need to be greater than 0 !')
  symbols=symb(1:nsub)
  do i=1,nsub
    if ( ichar(symbols(i)(1:1)) == 0 ) call stderr('The number of SYMB is mismatch !')
  end do
  inquire(file='SPOSCAR',exist=alive)
  if ( alive .eqv. .false. ) call stderr('The SPOSCAR file does not exist !')
  call readpos(a,x,atom,natom,'SPOSCAR')
  if ( site > size(natom) .or. site < 1 ) call stderr('The value of SITE is not allowed !')
  if ( sum(k) /= natom(site) ) call stderr('The sum of SUBS is incorrect !')
  call stdout_1(natom,atom,site,k,symbols)
  call eqamatrix(eqamat,spgmat,a,x,natom,prec,site)
  call grouping(group,eqamat,k(1),natom(site),a,x,atom,natom,site)
  call irrconfig(iconf,deg,eqamat,group,k,lpro)
  call output(eqamat,leqa,spgmat,lspg,iconf,deg,k,lcfg,a,x,atom,natom,symbols,site,lpos) 
  call cpu_time(time1)
  call stdout_5(nint(time1-time0))
end program disorder
