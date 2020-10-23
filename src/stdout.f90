module stdout
  use functions
  implicit none

contains

subroutine stdout_0
  implicit none

  write(*,*) "     _ _                   _           "
  write(*,*) "  __| (_)___  ___  _ __ __| | ___ _ __ "
  write(*,*) " / _` | / __|/ _ \| '__/ _` |/ _ \ '__|"
  write(*,*) "| (_| | \__ \ (_) | | | (_| |  __/ |   "
  write(*,*) " \__,_|_|___/\___/|_|  \__,_|\___|_|  "
  call system('echo -e "      `date +%F\ %H:%M:%S`\c"')
  write(*,*) "      0.5.1"
  write(*,*)
  return
end subroutine stdout_0

subroutine stdout_1(natom,atom,site,k,symbols)
  implicit none
  integer(2) :: natom(:),k(:),na
  integer(1) :: site,ntype,i
  character(len=2) :: atom(:),symbols(:)
  character(len=20) :: fm1,fm2

  write(*,'(A)') '  Reading INDSOD and SPOSCAR ...'
  ntype=size(natom)
  na=sum(natom)
  write(fm1,*) width(i_1=ntype)
  write(fm2,*) width(i_2=na)
  write(*,'(A,I'//fm1//',A,I'//fm2//',A,$)') '  Found ',ntype,' types and ',na,' atoms ('
  do i=1,ntype
    write(fm1,*) width(i_2=natom(i))
    write(fm2,*) len_trim(adjustl(atom(i)))
    write(*,'(1X,I'//fm1//',1X,A'//fm2//',$)') natom(i),atom(i)
  end do
  write(*,'(A)') ' )'
  write(fm1,*) len_trim(adjustl(atom(site)))
  write(*,'(A,A'//fm1//',A,$)') '  The ',atom(site),' site will be substituted by'
  do i=1,size(k)
    write(fm1,*) width(i_2=k(i))
    write(fm2,*) len_trim(adjustl(symbols(i)))
    write(*,'(1X,I'//fm1//',1X,A'//fm2//',$)') k(i),symbols(i)
  end do
  write(*,'(/)')
  return
end subroutine stdout_1

subroutine stdout_2(no,nr,nt,cls,sys)
  implicit none
  integer(2) :: no,nr,nt
  character(len=20) :: cls,sys
  character(len=20) :: fm1,fm2,fm3

  write(fm1,*) width(i_2=no)
  write(fm2,*) width(i_2=nr)
  write(fm3,*) width(i_2=nt)
  write(*,'(A)') '  Searching Space Group Operations ...'
  write(*,'(A,I'//fm1//',A,I'//fm2//',A,I'//fm3//',A)') '  Found ',&
       no,' operations ( ',nr,' rotations and ',nt,' pure translations )'
  write(fm1,*) len_trim(sys)
  write(fm2,*) len_trim(cls)
  if ( sys(1:1) == 'O' ) then
    write(*,'(A,$)') '  The supercell is an '
  else
    write(*,'(A,$)') '  The supercell is a '
  end if
  write(*,'(A'//fm1//',A,A'//fm2//',/)') sys,' lattice with a point group of ',cls
  write(*,'(A)') '  Preprocessing Work Related To Combinatorics ...'
  return
end subroutine stdout_2

subroutine stdout_3(nc)
  implicit none
  integer(8) :: nc
  character(len=20) fm1

  write(fm1,*) width(i_8=nc)
  write(*,'(A,I'//fm1//',A,/)') '  Found ',nc,' atomic configurations'
  write(*,'(A,$)') '  Eliminating Duplicate Configurations '
  return
end subroutine stdout_3

subroutine stdout_4(n)
  implicit none
  integer(4) :: n
  character(len=20) fm1

  write(fm1,*) width(i_4=n)
  write(*,'(A,I'//fm1//',A,/)') '  Found ',n,' irreducible configurations'
  return
end subroutine stdout_4

subroutine stdout_5(time)
  implicit none
  integer(4) :: time
  integer(1) :: h,m,s

  h=time/3600
  m=(time-h*3600)/60
  s=time-h*3600-m*60
  write(*,'(A,I2.2,A,I2.2,A,I2.2,A,/)') '  Program finished !    Elapsed time : ',&
       h,' hour ',m,' min ',s,' sec'

  write(*,'(A)') '+----------------------------------------------------------+'
  write(*,'(A)') '|                     * How to Site ? *                    |'
  write(*,'(A)') '| Please Cite The Following Article When You Use disorder: |'
  write(*,'(A)') '| [1] J.-C. Lian, H.-Y. Wu, W.-Q. Huang, W. Hu, and G.-F.  |'
  write(*,'(A)') '|     Huang, Phys. Rev. B, 102, 134209 (2020).             |'
  write(*,'(A)') '+----------------------------------------------------------+'
  write(*,*)
  return
end subroutine stdout_5

subroutine stderr(massage)
  implicit none
  character(len=*) massage

  write(*,'(/,2X,A)') massage
  write(*,*) '   _____ ____  ____   ___  ____  '
  write(*,*) '  | ____|  _ \|  _ \ / _ \|  _ \ '
  write(*,*) '  |  _| | |_) | |_) | | | | |_) |'
  write(*,*) '  | |___|  _ <|  _ <| |_| |  _ < '
  write(*,*) '  |_____|_| \_\_| \_\\___/|_| \_\'
  write(*,*)
  stop
end subroutine stderr

end module stdout
