module JSON_Parser
  !> @file json-parser.f90
  !! Multiline doxygen supported comment
  !! in Fortran  


  implicit none
  contains


  !> @brief Example function
  !! @param[in]  a  Example parameter
  !! @test Returns double input integer for a number of test cases
  function foo(a)
    integer, intent(in) :: a
    integer :: foo
    print *, a
    foo = 2*a
  end function
end module