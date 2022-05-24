program main
  use globals
  use json_parser
  use setup
  use solvers
  use command_line
  use comms

  implicit none

  real(dp) :: CH_params(6) ! equation parameters
  character(120) :: fname, run_name, outdir
  character(:), allocatable :: init
  integer :: grid_res ! grid level
  real(dp), allocatable :: Tout(:) ! output times
  real(dp), dimension(:,:), allocatable :: c0 ! initial concentration
  logical :: errors, all_params_fnd
  integer :: ierr, n

  ! lin/logspace params
  integer :: space_selected, num_outputs
  real(dp) :: start_val, stop_val

  ! start up MPI comms
  call comms_init()

  ! default fname, run_name, and outdir
  fname = "input-data.json"
  run_name = "default"
  outdir = "./out"

  ! Set up on rank 0 only
  if (myrank == 0) then
    ! CLI
    call initialise()

    ! Get JSON filename, run name, and output directory
    call get_io_commands(fname, run_name, outdir, all_params_fnd)

    ! input parameters from JSON file
    if (.NOT. all_params_fnd) then
      call read_json(trim(fname), trim(run_name), CH_params, init, grid_res, Tout)
    endif

    ! Grab any overriding input params from the command line
    call get_input_commands(CH_params, grid_res, init, Tout)

    ! Grab overriding linspace or logspace args
    ! TODO: should this not be in get_input_commands
    call get_lin_log_args(space_selected, start_val, stop_val, num_outputs)

    if (space_selected == LINSPACE_SELECTED) then
      call lin_tspace(start_val, stop_val, num_outputs, Tout)
    endif
    if (space_selected == LOGSPACE_SELECTED) then
      call log_tspace(start_val, stop_val, num_outputs, Tout)
    endif

    ! Validate input params
    call validate_params(CH_params, init, grid_res, Tout, errors)

    ! abort on error
    if (errors) then
      call logger%fatal("main", "Errors found in input params, aborting.")
      call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
    end if
  else
    grid_res = 0
  endif

  ! send setup to other procs
  n = 2**grid_res
  call broadcast_setup(CH_params, Tout, c0, n)

  if (myrank == 0) then
    ! initial concentration
    call setup_grid(c0, grid_res, ch_params, init)

    ! set up HDF5 outputting
    call output_init(outdir, [2, grid_res], CH_params, ierr)
  endif

  ! call solver (on all procs)
  call solver_1(Tout, c0, CH_params, SOLVER_FD2, ierr)

  if (myrank == 0) then
    ! finish HDF5 outputting
    call output_final(ierr)
  endif

  ! clean up (on all procs)
  deallocate(c0)
  deallocate(Tout)

end program main
