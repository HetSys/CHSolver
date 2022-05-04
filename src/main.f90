program main
  use globals
  use json_parser
  use setup
  implicit none
  real(dp) :: CH_params(6)
  character(:), allocatable :: init
  integer :: level
  real(dp), allocatable :: T(:)
  call initialise()

  call read_json("default", CH_params, init, level)
  

  call log_tspace(0.0_DP, 10.0_DP, 10, T)
  print *, T


end program main