module comms

  use solver_utils
  use mpi
  implicit none

  integer :: nproc, myrank, mpi_err    
  integer :: nproc_row, mycoords(2)
  integer :: cart_comm, neigh(4) !left, right, down, up

  contains

  subroutine comms_init()

    ! Assert that nrank

    call mpi_init(mpi_err)
    call mpi_comm_size(mpi_comm_world, nproc, mpi_err)


    call ilog2(nproc, nproc_row)
    nproc_row = nproc_row/2
    nproc_row = 2**nproc_row 

    call mpi_cart_create(mpi_comm_world, 2, [nproc_row, nproc_row], &
     [.true., .true.], .true., cart_comm, mpi_err)

    call mpi_comm_rank(mpi_comm_world, myrank, mpi_err)

    call mpi_cart_coords(cart_comm, myrank, 2, mycoords, mpi_err)


    call mpi_cart_shift(cart_comm, 1, 1, neigh(1), neigh(2), mpi_err)
    call mpi_cart_shift(cart_comm, 0, 1, neigh(4), neigh(3), mpi_err)

  end subroutine comms_init


  ! Send Initial Parameters Everywhere
  subroutine broadcast_setup(CH_params, grid_res)
    real(dp), intent(in), dimension(6) :: CH_params
    integer, intent(in) :: grid_res

    call mpi_bcast(CH_params, 6, mpi_float, 0, mpi_comm_world, mpi_err)
    call mpi_bcast(grid_res, 1, mpi_integer, 0, mpi_comm_world, mpi_err)


  end subroutine broadcast_setup

  subroutine grid_scatter(grid, grid_res, mpi_grid)
    integer, intent(in) :: grid_res
    real(dp), intent(in), dimension(grid_res,grid_res) :: grid    
    real(dp), intent(out), allocatable, dimension(:, :) :: mpi_grid

    integer :: rank, mpi_res
    integer :: rankcoords(2)
    integer :: dpsize
    integer :: subgrid_basic, subgrid
    integer(kind=mpi_address_kind) :: extent
    integer, dimension(nproc) :: displs, counts

    integer :: i, j

    mpi_res = int(grid_res/nproc_row)
    allocate(mpi_grid(0:mpi_res+1, 0:mpi_res+1))

    call mpi_type_create_subarray(2, [grid_res, grid_res], [mpi_res, mpi_res], [0, 0],&
     mpi_order_fortran, mpi_double_precision, subgrid_basic, mpi_err)

    call mpi_type_size(mpi_double_precision, dpsize, mpi_err)
    extent = mpi_res*dpsize
    call mpi_type_create_resized(subgrid_basic, int(0, mpi_address_kind), extent, subgrid, mpi_err)
    call mpi_type_commit(subgrid, mpi_err)

    do rank = 0, nproc-1

      call mpi_cart_coords(cart_comm, rank, 2, rankcoords, mpi_err)

      displs(rank+1) = rankcoords(2) + nproc_row*mpi_res*rankcoords(1)
      counts = 1
    end do
    mpi_grid = 0.0_dp

    call MPI_scatterv(grid, counts, displs, subgrid, &
     mpi_grid(1:mpi_res,1:mpi_res), mpi_res*mpi_res, mpi_double_precision, &
     0, mpi_comm_world, mpi_err)

    call send_edge(mpi_res+2, mpi_grid(mpi_res, :), mpi_grid(0, :), "d")
    call send_edge(mpi_res+2, mpi_grid(1, :), mpi_grid(mpi_res+1, :), "u")

    ! call send_edge(mpi_res+2, mpi_grid(:, mpi_res), mpi_grid(:, 0), "r")
    ! call send_edge(mpi_res+2, mpi_grid(:, 1), mpi_grid(:, mpi_res+1), "l")

    call MPI_Type_free(subgrid,mpi_err)

  end subroutine grid_scatter

  subroutine grid_gather(grid, grid_res, mpi_grid)
    integer, intent(in) :: grid_res
    real(dp), intent(in), dimension(grid_res,grid_res) :: grid    
    real(dp), intent(out), dimension(0:int(grid_res/nproc_row)+1, 0:int(grid_res/nproc_row)+1) :: mpi_grid

    integer :: rank, mpi_res
    integer :: rankcoords(2)
    integer :: dpsize
    integer :: subgrid_basic, subgrid
    integer(kind=mpi_address_kind) :: extent
    integer, dimension(nproc) :: displs, counts

    integer :: i, j

    mpi_res = int(grid_res/nproc_row)

    call mpi_type_create_subarray(2, [grid_res, grid_res], [mpi_res, mpi_res], [0, 0],&
     mpi_order_fortran, mpi_double_precision, subgrid_basic, mpi_err)

    call mpi_type_size(mpi_double_precision, dpsize, mpi_err)
    extent = mpi_res*dpsize
    call mpi_type_create_resized(subgrid_basic, int(0, mpi_address_kind), extent, subgrid, mpi_err)
    call mpi_type_commit(subgrid, mpi_err)

    do rank = 0, nproc-1

      call mpi_cart_coords(cart_comm, rank, 2, rankcoords, mpi_err)

      displs(rank+1) = rankcoords(2) + nproc_row*mpi_res*rankcoords(1)
      counts = 1
    end do

    call MPI_gatherv(mpi_grid(1:mpi_res,1:mpi_res), mpi_res*mpi_res, mpi_double_precision, &
    grid, counts, displs, subgrid, &
     0, mpi_comm_world, mpi_err)


  end subroutine grid_gather


  subroutine send_edge(n, sent, recv, direction)
    integer, intent(in) :: n
    real(dp), dimension(n), intent(in) :: sent
    real(dp), dimension(n), intent(out) :: recv
    character, intent(in) :: direction

    integer :: dir_int1, dir_int2

    select case(direction)

    case("l")
      dir_int1 = 4
      dir_int2 = 3

    case("r")
      dir_int1 = 3
      dir_int2 = 4

    case("u")
      dir_int1 = 1
      dir_int2 = 2

    case("d")
      dir_int1 = 2
      dir_int2 = 1

    end select

    call MPI_Sendrecv(sent, n, mpi_double_precision, neigh(dir_int1), 1003, recv, n, &
    mpi_double_precision, neigh(dir_int2), 1003, cart_comm, mpi_status_ignore, mpi_err)


  end subroutine send_edge


  subroutine comms_final()
    call mpi_finalize(mpi_err)
  end subroutine comms_final




end module comms
