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
  real(dp) :: t0, dt
  real(dp), allocatable :: Tout(:), updated_Tout(:) ! output times
  real(dp), dimension(:,:), allocatable :: c0, c1! initial concentration
  logical :: errors, all_params_fnd
  integer :: ierr
  integer :: st_Tout

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

  ! !!! END GRID AND HDF5 SETUP



  ! call solver

   call setup_grid(c0, grid_res, ch_params, init)
   call output_init(outdir, [2, grid_res], CH_params, ierr)

   call solver_1(Tout, c0, CH_params, SOLVER_PS, ierr)
   print*, c0

  !call chkpnt_init(outdir, ch_params, t0, ierr)
  !allocate(c0(c_dims(1), c_dims(2)), c1(c_dims(1), c_dims(2)))
  !call read_hdf5_chkpnt(c0, c1, dt, ierr)

  !do st_Tout = 1, size(Tout)
  !  if (Tout(st_Tout) > t0) then
  !    exit
  !  end if
  !end do

  !allocate(updated_Tout(size(Tout) - st_Tout))

  !updated_Tout = Tout(st_Tout:)

  !call solver_2(t0, updated_Tout, c0, c1, dt, CH_params, SOLVER_FD2, ierr)



  call output_final(ierr)



  ! clean up
  deallocate(c0)
  deallocate(Tout)

end program main
