program main
  use globals
  use json_parser
  use setup
  use solvers

  implicit none

  real(dp) :: CH_params(6) ! equation parameters
  character(:), allocatable :: init ! initial condition type
  integer :: level ! grid level
  real(dp), allocatable :: Tout(:) ! output times
  real(dp), pointer, contiguous :: c0(:,:) ! initial concentration

  ! set up logging
  call initialise()

  ! input parameters from JSON file
  call read_json("default", CH_params, init, level)

  ! set output times
  call log_tspace(0.0_dp, 5.0_dp, 10, Tout)

  ! initial concentration
  call setup_grid(c0, level, init)

  ! call solver
  call solver_1(Tout, c0, CH_params, SOLVER_FD2)

  ! clean up
  deallocate(c0)
  deallocate(Tout)

end program main
