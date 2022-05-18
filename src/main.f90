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
  integer :: grid_res ! grid level
  real(dp), allocatable :: Tout(:) ! output times
  real(dp), pointer, contiguous :: c0(:,:) ! initial concentration
  logical :: errors, all_params_fnd
  integer :: ierr

  ! Lin/Logspace params
  integer :: space_selected, num_outputs
  real(dp) :: start_val, stop_val

  ! Default fname, run_name, and outdir
  fname = "input-data.json"
  run_name = "default"
  outdir = "./out"

  !!! START  JSON, CLI, and LOGGING
  
  ! CLI 
  call initialise()

  ! Get JSON filename, run name, and output directory
  call get_io_commands(fname, run_name, outdir, all_params_fnd)

  ! input parameters from JSON file
  if (.NOT. all_params_fnd) call read_json(trim(fname), trim(run_name), CH_params, init, grid_res, Tout)

  ! Grab any overriding input params from the command line
  call get_input_commands(CH_params, grid_res, init, Tout)

  ! Grab overriding linspace or logspace args
  call get_lin_log_args(space_selected, start_val, stop_val, num_outputs)

  if (space_selected == LINSPACE_SELECTED) call lin_tspace(start_val, stop_val, num_outputs, Tout)
  if (space_selected == LOGSPACE_SELECTED) call log_tspace(start_val, stop_val, num_outputs, Tout)

  ! Validate input params
  call validate_params(CH_params, init, grid_res, Tout, errors)

  if (errors) then
    call logger%fatal("main", "Errors found in input params, aborting.")
    stop
  end if

  !!! END JSON, CLI, and LOGGING


  !!! START GRID AND HDF5 SETUP 
  ! initial concentration
  call setup_grid(c0, grid_res, init)

  call output_init(outdir, [2, grid_res], CH_params, ierr)

  !!! END GRID AND HDF5 SETUP 


 
  ! call solver
  call solver_1(Tout, c0, CH_params, SOLVER_FD2, ierr)


  call output_final(ierr)

  ! clean up
  deallocate(c0)
  deallocate(Tout)

end program main
