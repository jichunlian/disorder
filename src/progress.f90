module progress
  implicit none
  type :: cls_prog
    integer(8)  :: i,num
    integer(1) :: lens=50
    character(len=1) :: head='#',tail='-'
  contains
    procedure :: set
    procedure :: put
  end Type

contains

subroutine set(prog,num,lens,head,tail)
  class( cls_prog ) :: prog
  integer(8) :: num
  integer(1),optional :: lens
  character(len=1),optional :: head,tail

  prog%num=num
  if ( present(lens) ) prog%lens=lens
  if ( present(head) ) prog%head=head
  if ( present(tail) ) prog%tail=tail
end subroutine set
  
subroutine put(prog,i)
  class( cls_prog ) :: prog
  integer(8) :: i
  character(len=1) :: br
  integer(1) :: nh

  prog%i=i
  if ( prog%i > prog%num ) prog%i=prog%num    
  nh=nint(dble(prog%i)*prog%lens/dble(prog%num))
  if ( prog%i < prog%num ) then
    br = char(13)
  else
    br = char(10)
  end if
  write(*,'(4A,F6.2,2A,$)') '  [',repeat(prog%head,nh),&
         repeat(prog%tail,prog%lens-nh),'] ',prog%i*100.D0/prog%num,'%',br
end subroutine put
  
end module progress
