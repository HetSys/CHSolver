program main
  use globals
  use json_parser
  use setup
  use solvers
  use command_line

  implicit none

  real(dp) :: CH_params(6) ! equation parameters
  character(120) :: fname, run_name, outdir
  character(:), allocatable :: init
  integer :: level ! grid level
  real(dp), allocatable :: Tout(:) ! output times
  real(dp), pointer, contiguous :: c0(:,:) ! initial concentration
  logical :: errors, all_params_fnd

  ! Lin/Logspace params
  integer :: space_selected, num_outputs
  real(dp) :: start_val, stop_val

  ! Default fname, run_name, and outdir
  fname = "input-data.json"
  run_name = "default"
  outdir = "./out"

  ! set up logging
  call initialise()

  ! Get JSON filename, run name, and output directory
  call get_io_commands(fname, run_name, outdir, all_params_fnd)

  ! input parameters from JSON file
  if (.NOT. all_params_fnd) call read_json(trim(fname), trim(run_name), CH_params, init, level, Tout)

  ! Grab any overriding input params from the command line
  call get_input_commands(CH_params, level, init, Tout)

  ! Grab overriding linspace or logspace args
  call get_lin_log_args(space_selected, start_val, stop_val, num_outputs)

  if (space_selected == LINSPACE_SELECTED) call lin_tspace(start_val, stop_val, num_outputs, Tout)
  if (space_selected == LOGSPACE_SELECTED) call log_tspace(start_val, stop_val, num_outputs, Tout)

  ! Validate input params
  call validate_params(CH_params, init, level, Tout, errors)

  if (errors) then
    call logger%fatal("main", "Errors found in input params, aborting.")
    stop
  end if

  ! initial concentration
  call setup_grid(c0, level, init)

  ! call solver
  call solver_1(Tout, c0, CH_params, SOLVER_FD2)

  ! clean up
  deallocate(c0)
  deallocate(Tout)

end program main
