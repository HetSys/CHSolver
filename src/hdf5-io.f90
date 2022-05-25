  !> @brief Module containing all the subroutines for HDF5 outputting and checkpointing, used within main and solver code.
  !! @details A new HDF5 file (.chkpnt) is created for each output time requested, with a metadata.dat file storing the parameters used by the solver, as well as the times outputted at. 
module hdf5_io
  use hdf5
  use globals
  use iso_fortran_env

  implicit none

  !Format for Meta Data Records
  character(*), parameter, private :: metaformat = "(I5,4X,F15.5)", &
                                      metadata_fname = "/metadata.dat", &
                                      metadata2_fname = "/metadata2.dat"

  integer :: cur_chkpnt
  character(len=:), allocatable :: folder

  !Concentration Metadata
  integer(hid_t) :: c_dset_id !dataset id for c
  integer(hid_t) :: c_prev_dset_id !dataset id for c_prev
  integer(hid_t) :: c_dspace_id ! dataspace identifier for 2/3D
  integer(hsize_t), allocatable, dimension(:) :: c_dims !dimensions for dataspace for 2/3D
  
  !dt Metadata
  integer(hid_t) :: dt_dset_id !dataset id for c
  integer(hid_t) :: dt_dspace_id ! dataspace identifier for 3/4D
  integer(hsize_t), dimension(1) :: dt_dims

  contains


  !> @brief Creates a folder to store the trajectory and inside the folder a metadata file is created. HDF5 variables are allocated to fit grid lengths/ranks.
  !! @param[in]  foldername Folder name to store the checkpoint files and metadata.
  !! @param[in]  grid_params (1) Number of dimensions, (2) 2 raised to the power of this number will give the grid resolution.
  !! @param[in]  sys_params Parameters used by the solver, to be placed in metadata.
  !! @param[out]  error error code
  subroutine output_init(foldername, grid_params, sys_params, error)
    character(*), intent(in) :: foldername
    integer, intent(in), dimension(2) :: grid_params
    real(dp), intent(in), dimension(6) :: sys_params
    integer, intent(out) :: error
    
    integer :: iu

    allocate(c_dims(grid_params(1)))

    folder = foldername 
    c_dims = int(2**grid_params(2), hsize_t)
    dt_dims = 1

    call execute_command_line("rm "//trim(foldername)//"/*.chkpnt "//trim(foldername)//"/metadata.dat 2> /dev/null", wait=.true.)
    call execute_command_line("mkdir -p "//trim(foldername), wait=.true.)

    open(newunit=iu, file=trim(foldername)//metadata_fname, status="new")

    write(iu, "('grid_params',1X,I5,1X,I5)") grid_params(1), grid_params(2)
    write(iu, "('system_params',1X,F15.5,1X,F15.5,1X,F15.5)") sys_params(1), sys_params(2), sys_params(3)
    write(iu, "('system_params',1X,F15.5,1X,F15.5,1X,F15.5)") sys_params(4), sys_params(5), sys_params(6)
    write(iu, *)

    close(iu)

    cur_chkpnt = 0

    call h5open_f(error)

  end subroutine
  



  !> @brief Called from within solvers to output c, c_prev and dt to a HDF5 checkpoint file.
  !! @param[in]  c Concentration at current timestep.
  !! @param[in]  c_prev Concentration at previous timestep.
  !! @param[in]  dt Difference in time between the current and previous timesteps.
  !! @param[in]  t The time at the current step.
  !! @param[out]  error error code

  subroutine write_to_traj(c, c_prev, time, dt, error)
    real(dp), intent(in), dimension(:, :) :: c, c_prev
    real(dp), intent(in) :: time, dt
    integer, intent(out) :: error

    integer(hid_t) ::  file_id
    integer :: iu, i


    cur_chkpnt = cur_chkpnt + 1

    
    call h5fcreate_f(trim(folder)//"/"//TRIM(adjustl(to_string(cur_chkpnt)))//".chkpnt", h5f_acc_trunc_f, file_id, error)
    
    call h5screate_simple_f(2, c_dims, c_dspace_id, error)
    call h5screate_simple_f(1, dt_dims, dt_dspace_id, error)
    
    call h5dcreate_f(file_id,"c", h5t_native_double, c_dspace_id, c_dset_id, error)
    call h5dcreate_f(file_id,"c_prev", h5t_native_double, c_dspace_id, c_prev_dset_id, error)
    call h5dcreate_f(file_id,"dt", h5t_native_double, dt_dspace_id, dt_dset_id, error)
        
    call h5dwrite_f(c_dset_id, h5t_native_double, c, c_dims, error, c_dspace_id)
    call h5dwrite_f(c_prev_dset_id, h5t_native_double, c_prev, c_dims, error, c_dspace_id)
    call h5dwrite_f(dt_dset_id, h5t_native_double, dt, dt_dims, error, dt_dspace_id)

    call h5dclose_f(c_dset_id,error)
    call h5dclose_f(c_prev_dset_id,error)
    call h5dclose_f(dt_dset_id,error)

    call h5sclose_f(c_dspace_id,error)
    call h5sclose_f(dt_dspace_id,error)

    call h5fclose_f(file_id, error)

    open(newunit=iu, file=trim(folder)//metadata_fname, status="old")
    
    do i = 1, cur_chkpnt + 3
      read(iu, *)
    end do

    write(iu, metaformat) cur_chkpnt, time
    
    close(iu)


  end subroutine

  
  !> @brief Closes interfaces.
  !! @param[out]  error error code
  subroutine output_final(error)
    integer, intent(out) :: error

    integer :: iu, iu2, i, current
    real(dp) :: time
    character(128) :: buffer

    open(newunit=iu, file=trim(folder)//metadata_fname, status="old")
    open(newunit=iu2, file=trim(folder)//metadata2_fname, status="new")

    do i =1, 3
      read(iu, "(A)") buffer
      write(iu2, "(A)") buffer
    end do

    read(iu, *)
    write(iu2, *)

    do i = 1, cur_chkpnt
      read(iu, metaformat) current, time
      write(iu2, metaformat) current, time
    end do

    close(iu, status="delete")
    close(iu2)

    call rename(trim(folder)//metadata2_fname, trim(folder)//metadata_fname)


    call h5close_f(error)

  end subroutine

  !> @brief Reads in a HDF5 checkpoint file, used in restarting
  !! @param[in]  c Concentration at current timestep.
  !! @param[in]  c_prev Concentration at previous timestep.
  !! @param[in]  dt Difference in time between the current and previous timesteps.
  !! @param[out]  error error code
  subroutine read_hdf5_chkpnt(c, c_prev, dt, err)
    real(dp), dimension(:, :), intent(inout), allocatable :: c, c_prev
    real(dp), intent(out) :: dt
    integer, intent(out) :: err
    
    integer(hid_t) :: file_id

    ! allocate(c(c_dims(1), c_dims(2)), c_prev(c_dims(1), c_dims(2)))

    call h5fopen_f(trim(folder)//"/"//TRIM(adjustl(to_string(cur_chkpnt)))//".chkpnt", &
     h5f_acc_rdonly_f, file_id, err)

    call h5dopen_f(file_id, "c", c_dset_id, err)
    call h5dopen_f(file_id, "c_prev", c_prev_dset_id, err)
    call h5dopen_f(file_id, "dt", dt_dset_id, err)

    call h5dread_f(c_dset_id, h5t_native_double, c, c_dims, err)
    call h5dread_f(c_prev_dset_id, h5t_native_double, c_prev, c_dims, err)
    call h5dread_f(dt_dset_id, h5t_native_double, dt, dt_dims, err)

    call h5dclose_f(c_dset_id,err)
    call h5dclose_f(c_prev_dset_id,err)
    call h5dclose_f(dt_dset_id,err)

  end subroutine 

  !> @brief Initialises a restart from a specified checkpoint.
  !! @param[in]  chkpnt_folder Name of folder to store checkpoints and metadatai
  !! @param[out] n grid length in checkpoints
  !! @param[out]  sys_params Parameters used by the solver, to be placed in metadata.
  !! @param[out]  current_time Time of the system at the checkpoitn selected
  !! @param[out]  error error code
  !! @param[in]  n_chkpnt (Optional) Restart from the n-th checkpoint. Default is latest checkpoint.
  !! @param[in]  start_before_time (Optional) Selects latest checkpoint with time before start_before_time. n_chkpnt supercedes this, if both are supplied, in the case of conflict.

  subroutine chkpnt_init(chkpnt_folder, n,  sys_params, current_time, error, n_chkpnt, start_before_time)
    character(*), intent(in) :: chkpnt_folder
    integer, intent(out) :: n
    integer, optional, intent(in) :: n_chkpnt
    integer, intent(out) :: error
    real(dp), optional, intent(in) :: start_before_time
    real(dp), intent(out) :: current_time
    real(dp), intent(out) , dimension(6):: sys_params

    integer :: grid_rank, grid_res

    real(dp) :: time

    integer :: iu, iu2, pos, tot_chkpnt, i, current, res_chkpnt, ios
    character(128) :: buffer

    current = 1


    open(newunit=iu, file=trim(chkpnt_folder)//trim(metadata_fname), status="old")

    ios = 0
    tot_chkpnt = 0

    do while (ios == 0)
      read(iu, *, iostat=ios)
      tot_chkpnt = tot_chkpnt + 1 
    end do

    tot_chkpnt = tot_chkpnt - 5

    close(iu)

    open(newunit=iu, file=trim(chkpnt_folder)//trim(metadata_fname), status="old")

    read(iu, "(A)") buffer
    pos = scan(buffer, " ")
    buffer = buffer(pos+1:)
    read(buffer, "(I5,1X,I5)") grid_rank, grid_res

    read(iu, "(A)") buffer
    pos = scan(buffer, " ")
    buffer = buffer(pos+1:)
    read(buffer, "(F15.5,1X,F15.5,1X,F15.5)") sys_params(1), sys_params(2), sys_params(3)
    
    read(iu, "(A)") buffer
    pos = scan(buffer, " ")
    buffer = buffer(pos+1:)
    read(buffer, "(F15.5,1X,F15.5,1X,F15.5)") sys_params(4), sys_params(5), sys_params(6)
    
    read(iu, *)
    

    ! Start from a numbered checkpoint (res_chkpnt) provided
    if(present(n_chkpnt)) then

      res_chkpnt = n_chkpnt


      if(present(start_before_time)) then
        !call logger%error("start_from_chkpnt", "Only specify which checkpoint or before which time, defaulting to checkpint")
      end if

        if (n_chkpnt > tot_chkpnt) then
          !call logger%error("start_from_chkpnt", "checkpoint specified was never reached, defaulting to final checkpoint")
          res_chkpnt = tot_chkpnt
        end if

        do i =1, res_chkpnt-1
          read(iu, *)
        end do
        read(iu, metaformat) res_chkpnt, current_time


    ! Start from the latest checkpoint before a specified time (start_before_time)
    else if (present(start_before_time)) then



      do while (current .ne. tot_chkpnt)

        read(iu, *) current, current_time
        

        if (current_time .ge. start_before_time) then
        
          close(iu)
          open(newunit=iu, file=trim(chkpnt_folder)//metadata_fname, status="old")

          do i =1,current+3
            read(iu, *)
          end do

          read(iu, metaformat) res_chkpnt, current_time
          exit

        end if

      end do

      if (tot_chkpnt == current) then
        call logger%error("start_from_chkpnt", "Time specified was never reached, choosing the last chkpnt")
        res_chkpnt = tot_chkpnt
      end if


    ! No optional parameter specified, defaulting to last checkpoint
    else 
      res_chkpnt = tot_chkpnt
      do i=1,res_chkpnt-1
        read(iu, *)
      end do
      read(iu, metaformat) res_chkpnt, current_time
    end if


    close(iu)


    open(newunit=iu, file=trim(chkpnt_folder)//metadata_fname, status="old")
    open(newunit=iu2, file=trim(chkpnt_folder)//metadata2_fname, status="new")


    write(iu2, "('grid_params',1X,I5,1X,I5)") grid_rank, grid_res
    write(iu2, "('system_params',1X,F15.5,1X,F15.5,1X,F15.5)") sys_params(1), sys_params(2), sys_params(3)
    write(iu2, "('system_params',1X,F15.5,1X,F15.5,1X,F15.5)") sys_params(4), sys_params(5), sys_params(6)

    do i =1, 4
      read(iu, *)
    end do

    write(iu2, *)

    do i = 1, res_chkpnt
      read(iu, metaformat) current, time
      write(iu2, metaformat) current, time
    end do

    close(iu, status="delete")
    close(iu2)


    call rename(trim(chkpnt_folder)//metadata2_fname, trim(chkpnt_folder)//metadata_fname )

    do i = res_chkpnt+1, tot_chkpnt 
      open(newunit=iu, file=trim(chkpnt_folder)//"/"//trim(adjustl(to_string(i)))//".chkpnt", status="old")
      close(iu, status="delete")
    end do

    cur_chkpnt = res_chkpnt

    if (allocated(c_dims)) then
      deallocate(c_dims)
    end if

    allocate(c_dims(grid_rank))
    n = 2**grid_res
    c_dims = int(n, hsize_t)
    dt_dims = 1

    folder = chkpnt_folder

    call h5open_f(error)

  end subroutine 

end module
