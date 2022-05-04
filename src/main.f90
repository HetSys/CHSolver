program main
  use globals
  use json_parser
  implicit none
  real(dp) :: CH_params(6)
  character(:), allocatable :: init
  integer :: level
  call initialise()

  call read_json("default", CH_params, init, level)


end program main