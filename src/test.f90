program main

  use comms
  implicit none

  integer :: err
  real(dp) :: grid(4,4)
  real(dp), allocatable :: mpi_grid(:, :)
  integer :: i, j

  call comms_init()


  grid(1:2, 1:2) = 1
  grid(3:4, 1:2) = 2
  grid(1:2, 3:4) = 3
  grid(3:4, 3:4) = 4


  call grid_scatter(grid, 4, mpi_grid)

  ! do i = 1, 2
  !   do j = 1, 2

  !     mpi_grid(i, j) = mpi_grid(i-1,j) + mpi_grid(i+1,j) + mpi_grid(i,j-1) + mpi_grid(i,j+1)
  !     mpi_grid(i,j) = mpi_grid(i,j)/4

  !   end do
  ! end do

  if (myrank .ne. 0) goto 30
  
  do i = 1,4
    do j = 1,4
      write(*, "(F5.2, 1X)", advance="no") grid(i, j)
    end do
    write(*,*)
  end do
  write(*, *)
  write(*, *) myrank
  do i = 0,3
    do j = 0,3
      write(*, "(F5.2, 1X)", advance="no") mpi_grid(i, j)
    end do
    write(*,*)
  end do
  write(*, *)

  30 call mpi_barrier(mpi_comm_world, mpi_err)
 
  
  if (myrank .ne. 1) goto 31

  write(*, *) myrank
  do i = 0,3
    do j = 0,3
      write(*, "(F5.2, 1X)", advance="no") mpi_grid(i, j)
    end do
    write(*,*)
  end do
  write(*, *)

  31 call mpi_barrier(mpi_comm_world, mpi_err)


  
  if (myrank .ne. 2) goto 32

  write(*, *) myrank
  do i = 0,3
    do j = 0,3
      write(*, "(F5.2, 1X)", advance="no") mpi_grid(i, j)
    end do
    write(*,*)
  end do
  write(*, *)

  32 call mpi_barrier(mpi_comm_world, mpi_err)

    
  if (myrank .ne. 3) goto 33

  write(*, *) myrank
  do i = 0,3
    do j = 0,3
      write(*, "(F5.2, 1X)", advance="no") mpi_grid(i, j)
    end do
    write(*,*)
  end do
  write(*, *)

  33 call mpi_barrier(mpi_comm_world, mpi_err)

  call grid_gather(grid, 4, mpi_grid)


  if (myrank .ne. 0) goto 40
  
  do i = 1,4
    do j = 1,4
      write(*, "(F5.2, 1X)", advance="no") grid(i, j)
    end do
    write(*,*)
  end do
  write(*, *)


  40 call comms_final()

end program
