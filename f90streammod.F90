module f90streammod
  !use iso_c_binding
  use, intrinsic :: iso_c_binding
  implicit none

  integer,private,save :: cout = 6

  public :: f90stream_get_cout_unit
  public :: f90stream_set_cout_unit
  public :: f90stream_flush_cout
  private 

  interface
     subroutine f90stream_flush_cout() bind(c,name='f90stream_flush_cout')
       implicit none
     end subroutine f90stream_flush_cout
  end interface

contains  

  subroutine f90stream_cout(n,cstring) &
       & bind(c,name='f90stream_cout')
    !use iso_c_binding
    use, intrinsic :: iso_c_binding
    integer(kind=c_int),intent(in),value :: n
    character(kind=c_char,len=1),dimension(n),intent(in) :: cstring
    character(len=16) fmt
    if ( n > 0 ) then
       write(fmt,'(A,I13,A)')"(",n,"A)"
       write(cout,fmt,ADVANCE="no")cstring
    end if
  end subroutine f90stream_cout

  function f90stream_get_cout_unit() result(i) &
       & bind(c,name='f90stream_get_cout_unit')
    !use iso_c_binding
    use, intrinsic :: iso_c_binding
    integer(kind=c_int) :: i
    i = cout
  end function f90stream_get_cout_unit

  subroutine f90stream_set_cout_unit(i) &
       & bind(c,name='f90stream_set_cout_unit')
    !use iso_c_binding
    use, intrinsic :: iso_c_binding
    integer(kind=c_int),intent(in),value :: i
    call f90stream_flush_cout()
    cout = i
  end subroutine f90stream_set_cout_unit

end module f90streammod

