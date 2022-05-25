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
  real(dp), dimension(:,:), allocatable :: c0, c1 ! initial concentration
  logical :: errors, all_params_fnd
  integer :: ierr, n

  ! lin/logspace params
  integer :: space_selected, num_outputs
  real(dp) :: start_val, stop_val

  logical :: do_restart
  integer :: checkpoint_number, st_tout
  real(dp) :: restart_time, dt
  real(dp), allocatable :: updated_Tout(:) ! output times

  ! MPI timing
  real(dp) :: t0

  integer :: selected_solver

  ! start up MPI comms
  call comms_init()

  selected_solver = 2 ! Default to ps solver
  do_restart = .false.
  t0 = 0.0_dp
  if (myrank == 0) then
    t0 = MPI_Wtime()
  endif

  ! default fname, run_name, and outdir
  fname = "input-data.json"
  run_name = "default"
  outdir = "./out"

  ! Set up on rank 0 only
  if (myrank == 0) then
    ! CLI
    call initialise()

    ! Get solver selection
    call get_selected_solver(selected_solver)

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

    ! check if restart required
    restart_time = -1.0_dp
    checkpoint_number = -1
    call get_checkpoint_commands(checkpoint_number, restart_time, do_restart)

    if (do_restart) then
      if (checkpoint_number > 0) then
        call chkpnt_init(outdir, ch_params, t0, ierr, n_chkpnt=checkpoint_number)
      else if (restart_time >= 0.0_dp) then
        call chkpnt_init(outdir, ch_params, t0, ierr, start_before_time=restart_time)
      endif
      allocate(c0(c_dims(1), c_dims(2)))
      allocate(c1(c_dims(1), c_dims(2)))
      call read_hdf5_chkpnt(c0, c1, dt, ierr)
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
  call broadcast_setup(CH_params, Tout, c0, n, do_restart, selected_solver)

  ! call solver (on all procs)
  if (.not. do_restart) then
    if (myrank == 0) then
      ! initial concentration
      call setup_grid(c0, grid_res, ch_params, init)

      ! set up HDF5 outputting
      call output_init(outdir, [2, grid_res], CH_params, ierr)
    endif
    
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call solver_1(Tout, c0, CH_params, selected_solver, ierr)
  else
    ! ensure that Tout(1) >= t0
    if (myrank == 0) then
      do st_tout=1,size(tout)
        if (tout(st_tout) > t0) then
          exit
        endif
      enddo
    endif
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    allocate(updated_tout(size(tout)-st_tout))
    updated_tout = Tout(st_tout:)
    call solver_2(t0, updated_tout, c0, c1, dt, CH_params, selected_solver, ierr)
  endif

  if (myrank == 0) then
    ! finish HDF5 outputting
    call output_final(ierr)
  endif

  ! clean up (on all procs)
  deallocate(c0)
  deallocate(Tout)

  if (myrank == 0) then
    call logger%info("main", "total time: "//to_string(MPI_Wtime()-t0))
  endif

  ! shut down comms
  call comms_final()

end program main
