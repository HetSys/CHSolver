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
  real(dp) :: t0, dt, restart_time
  real(dp), allocatable :: Tout(:), updated_Tout(:) ! output times
  real(dp), dimension(:,:), allocatable :: c0, c1! initial concentration
  logical :: errors, all_params_fnd, do_restart
  integer :: ierr
  integer :: st_Tout, checkpoint_number

  ! Lin/Logspace params
  integer :: space_selected, num_outputs
  real(dp) :: start_val, stop_val

  ! Solver Selection
  integer :: selected_solver

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


  ! Grab overriding restarts args

  restart_time = -99.0_DP
  checkpoint_number = -99

  call get_checkpoint_commands(checkpoint_number, restart_time, do_restart)
  
  if(do_restart) then
    if (checkpoint_number > 0) then
      print *, "Number"
      call chkpnt_init(outdir, ch_params, t0, ierr, n_chkpnt=checkpoint_number)
    else if (restart_time >= 0) then
      print *, "Time"
      call chkpnt_init(outdir, ch_params, t0, ierr, start_before_time=restart_time)
    end if
  end if



  ! Validate input params
  call validate_params(CH_params, init, grid_res, Tout, errors)

  if (errors) then
    call logger%fatal("main", "Errors found in input params, aborting.")
    stop 1
  end if

  !!! END JSON, CLI, and LOGGING

  ! SOLVER SELECTION
  call get_selected_solver(selected_solver)

  ! GRID INITIALISATION

  if (.NOT. do_restart) then
   call setup_grid(c0, grid_res, ch_params, init)
   call output_init(outdir, [2, grid_res], CH_params, ierr)

   
   call solver_1(Tout, c0, CH_params, selected_solver, ierr)
  else

  ! call solver_1(Tout, c0, CH_params, SOLVER_FD2, ierr)

    call chkpnt_init(outdir, ch_params, t0, ierr)
    allocate(c0(c_dims(1), c_dims(2)), c1(c_dims(1), c_dims(2)))
    call read_hdf5_chkpnt(c0, c1, dt, ierr)

    do st_Tout = 1, size(Tout)
      if (Tout(st_Tout) > t0) then
        exit
      end if 
    end do

    allocate(updated_Tout(size(Tout) - st_Tout))

    updated_Tout = Tout(st_Tout:)
    call solver_2(t0, updated_Tout, c0, c1, dt, CH_params, selected_solver, ierr)
  end if

  
  
  call output_final(ierr)



  ! clean up
  deallocate(c0)
  deallocate(Tout)

end program main
