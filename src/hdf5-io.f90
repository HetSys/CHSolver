! Metadata file format (as of now)
!   Start of file:  "nchkpnts {(I5)}"
!   Buffer Line:
!   Records Line: "{(I5)}  {(F15.5)}" representing chkpnt and time when chkpnt is taken
! 

module hdf5_io
  use hdf5
  use globals
  use iso_fortran_env

  implicit none

  !Format for Meta Data Records
  character(*), parameter :: metaformat = "(I5,4X,F15.5)"

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
  integer(hsize_t), allocatable, dimension(:) :: dt_dims

  interface write_to_traj
    module procedure write_to_traj_2D
    module procedure write_to_traj_3D
  end interface write_to_traj

  contains

  !This function should create a folder to store the trajectory 
  !In the folder, a metadata file should be created
  !Alongside this, hdf5 variables should be allocated to fit grid lengths/ranks
  subroutine output_init(foldername, grid_params, sys_params, error)
    character(*), intent(in) :: foldername
    integer, intent(in), dimension(2) :: grid_params
    real(dp), intent(in), dimension(6) :: sys_params
    integer, intent(out) :: error
    
    integer :: iu

    allocate(c_dims(grid_params(1)))
    allocate(dt_dims(1))

    folder = foldername 
    c_dims = int(2**grid_params(2), hsize_t)
    dt_dims = 1

    call execute_command_line("rm -r "//trim(foldername), wait=.true.)
    call execute_command_line("mkdir "//trim(foldername), wait=.true.)

    open(newunit=iu, file=trim(foldername)//"/metadata", status="new")

    write(iu, "('chkpnts',1X,I5)") 0
    write(iu, "('grid_params',1X,I5,1X,I5)") grid_params(1), grid_params(2)
    write(iu, "('system_params',1X,F15.5,1X,F15.5,1X,F15.5)") sys_params(1), sys_params(2), sys_params(3)
    write(iu, "('system_params',1X,F15.5,1X,F15.5,1X,F15.5)") sys_params(4), sys_params(5), sys_params(6)
    write(iu, *)

    close(iu)

    cur_chkpnt = 0

    call h5open_f(error)

  end subroutine

  subroutine write_to_traj_2D(c, c_prev, time, dt, error)
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
    
    ! print *, c
    
    call h5dwrite_f(c_dset_id, h5t_native_double, c, c_dims, error, c_dspace_id)
    call h5dwrite_f(c_prev_dset_id, h5t_native_double, c_prev, c_dims, error, c_dspace_id)
    call h5dwrite_f(dt_dset_id, h5t_native_double, dt, dt_dims, error, dt_dspace_id)

    call h5dclose_f(c_dset_id,error)
    call h5dclose_f(c_prev_dset_id,error)
    call h5dclose_f(dt_dset_id,error)

    call h5sclose_f(c_dspace_id,error)
    call h5sclose_f(dt_dspace_id,error)

    call h5fclose_f(file_id, error)

    open(newunit=iu, file=trim(folder)//"/metadata", status="old")
    
    do i = 1, cur_chkpnt + 4
      read(iu, *)
    end do

    write(iu, metaformat) cur_chkpnt, time
    
    close(iu)


  end subroutine

  subroutine write_to_traj_3D(c, c_prev, time, dt, error)
    real(dp), intent(in), dimension(:, :, :) :: c, c_prev
    real(dp), intent(in) :: time, dt 
    integer, intent(out) :: error

    integer(hid_t) ::  file_id
    integer :: iu, i

    cur_chkpnt = cur_chkpnt + 1

    
    call h5fcreate_f(trim(folder)//"/"//TRIM(adjustl(to_string(cur_chkpnt)))//".chkpnt", h5f_acc_trunc_f, file_id, error)
    
    call h5screate_simple_f(3, c_dims, c_dspace_id, error)
    call h5screate_simple_f(1, dt_dims, dt_dspace_id, error)
    
    call h5dcreate_f(file_id,"c", h5t_native_double, c_dspace_id, c_dset_id, error)
    call h5dcreate_f(file_id,"c_prev", h5t_native_double, c_dspace_id, c_prev_dset_id, error)
    call h5dcreate_f(file_id,"dt", h5t_native_double, dt_dspace_id, dt_dset_id, error)
    
    call h5dwrite_f(c_dset_id, h5t_native_real, c, c_dims, error, c_dspace_id)
    call h5dwrite_f(c_prev_dset_id, h5t_native_real, c_prev, c_dims, error, c_dspace_id)
    call h5dwrite_f(dt_dset_id, h5t_native_real, dt, dt_dims, error, dt_dspace_id)

    call h5dclose_f(c_dset_id,error)
    call h5dclose_f(c_prev_dset_id,error)
    call h5dclose_f(dt_dset_id,error)

    call h5sclose_f(c_dspace_id,error)
    call h5sclose_f(dt_dspace_id,error)

    call h5fclose_f(file_id, error)

    open(newunit=iu, file=trim(folder)//"/metadata", status="old")
    
    do i = 1, cur_chkpnt + 4
      read(iu, *)
    end do

    write(iu, metaformat) cur_chkpnt, time
    
    close(iu)

  end subroutine

  subroutine output_final(error)
    integer, intent(out) :: error

    integer :: iu, iu2, i, current
    real(dp) :: time
    character(128) :: buffer

    open(newunit=iu, file=trim(folder)//"/metadata", status="old")
    open(newunit=iu2, file=trim(folder)//"/metadata2", status="new")


    read(iu, *)
    write(iu2, "('chkpnts',1X,I5)") cur_chkpnt

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

    call rename(trim(folder)//"/metadata2", trim(folder)//"/metadata" )


    call h5close_f(error)

  end subroutine


  ! subroutine read_hdf5_slice(filename, c, c_prev, dt)
  !   character(*), intent(in) :: filename
  !   real(dp), dimension(:, :) :: c, c_prev
  !   real(dp) :: dt
  !   integer 
      



  ! end subroutine 

  subroutine continue_from_chkpnt(chkpnt_folder, sys_params, current_time, error, n_chkpnt, start_before_time)
    character(*), intent(in) :: chkpnt_folder
    integer, optional, intent(in) :: n_chkpnt
    integer, intent(out) :: error
    real(dp), optional, intent(in) :: start_before_time
    real(dp), intent(out) :: current_time
    real(dp), intent(out) , dimension(6):: sys_params

    integer :: grid_rank, grid_len

    real(dp) :: time

    integer :: iu, iu2, pos, tot_chkpnt, i, current, res_chkpnt
    character(128) :: buffer

    current = 1

 
    open(newunit=iu, file=chkpnt_folder//"/metadata", status="old")

    read(iu, "(A)") buffer
    pos = scan(buffer, " ")
    buffer = buffer(pos+1:)
    read(buffer, "(I5)") tot_chkpnt

    read(iu, "(A)") buffer
    pos = scan(buffer, " ")
    buffer = buffer(pos+1:)
    read(buffer, "(I5,1X,I5)") grid_rank, grid_len

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
          open(newunit=iu, file=trim(chkpnt_folder)//"/metadata", status="old")

          do i =1,current+3
            read(iu, *)
          end do

          read(iu, metaformat) res_chkpnt, current_time
          exit

        end if

      end do

      if (tot_chkpnt == current) then
        !call logger%error("start_from_chkpnt", "Time specified was never reached, choosing the last chkpnt")
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


    open(newunit=iu, file=trim(chkpnt_folder)//"/metadata", status="old")
    open(newunit=iu2, file=trim(chkpnt_folder)//"/metadata2", status="new")


    read(iu, *)
    write(iu2, "('chkpnts',1X,I5)") n_chkpnt

    write(iu2, "('grid_params',1X,I5,1X,I5)") grid_rank, grid_len
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


    call rename(trim(chkpnt_folder)//"/metadata2", trim(chkpnt_folder)//"/metadata" )

    do i = res_chkpnt+1, tot_chkpnt 
      open(newunit=iu, file=trim(chkpnt_folder)//"/"//trim(adjustl(to_string(i)))//".chkpnt", status="old")
      close(iu, status="delete")
    end do

    cur_chkpnt = res_chkpnt

    if (allocated(c_dims)) then
      deallocate(c_dims)
    end if

    allocate(c_dims(grid_rank))
    c_dims = int(grid_len, hsize_t)

    folder = chkpnt_folder



    call h5open_f(error)

  end subroutine 

end module