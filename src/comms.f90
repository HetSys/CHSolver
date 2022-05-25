!> @brief Module wrapping MPI boilerplate into useful functions.
module comms

  use solver_utils
  use mpi
  implicit none

  integer :: nproc, myrank, mpi_err    
  integer :: nproc_row, mycoords(2)
  integer :: cart_comm
  integer :: neigh(4) !left, right, down, up
  integer :: daig_neigh(4) !ul, ur, dl, dr
  contains

  !> @brief Sets up MPI communications
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

  !> @brief Sends initial parameters to all processes
  !! @param[inout] CH_params        equation and domain parameters
  !! @param[inout] Tout             output times
  !! @param[inout] grid_res         grid resolution
  !! @param[inout] do_restart       whether we are restarting from a checkpoint
  !! @param[inout] to               initial time
  !! @param[inout] selected_solver  solver type
  subroutine broadcast_setup(CH_params, Tout, grid, grid_res, do_restart, t0, selected_solver)
    real(dp), intent(inout), dimension(6) :: CH_params
    real(dp), intent(inout), allocatable, dimension(:) :: Tout
    real(dp), dimension(:,:), allocatable :: grid
    logical, intent(inout) :: do_restart
    integer, intent(inout) :: grid_res, selected_solver
    real(dp), intent(inout) :: t0
    integer :: Tout_size

    if (myrank == 0) then
      Tout_size = size(Tout)
    else
      Tout_size = 0
    endif

    ! broadcast sizes
    call mpi_bcast(grid_res, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
    call mpi_bcast(Tout_size, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
    call mpi_bcast(t0, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_err)

    call mpi_bcast(do_restart, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpi_err)
    call mpi_bcast(selected_solver, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)


    ! allocate on non-master procs
    if (myrank /= 0) then
      allocate(Tout(Tout_size))
      allocate(grid(grid_res,grid_res))
    endif

    ! broadcast data
    call mpi_bcast(CH_params, 6, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_err)
    call mpi_bcast(Tout, size(Tout), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpi_err)
  end subroutine broadcast_setup

  !> @brief Scatters a grid from rank 0 to mpi_grids on all processes
  !! @param[in] grid       global grid
  !! @param[in] grid_res   global grid resolution
  !! @param[out] mpi_grid  local grid
  !! @param[in] mpi_res    local grid resolution
  subroutine grid_scatter(grid, grid_res, mpi_grid, mpi_res)
    integer, intent(in) :: grid_res
    integer, intent(in) :: mpi_res
    real(dp), intent(in), dimension(grid_res,grid_res) :: grid
    real(dp), intent(out), dimension(0:mpi_res+1, 0:mpi_res+1) :: mpi_grid

    integer :: rank
    integer :: rankcoords(2)
    integer :: dpsize
    integer :: subgrid_basic, subgrid
    integer(kind=mpi_address_kind) :: extent
    integer, dimension(nproc) :: displs, counts

    integer :: req1, req2, req3, req4
    integer :: req5, req6, req7, req8


    ! mpi_res = int(grid_res/nproc_row)
    ! allocate(mpi_grid(0:mpi_res+1, 0:mpi_res+1))

    call mpi_type_create_subarray(2, [grid_res, grid_res], [mpi_res, mpi_res], [0, 0],&
     mpi_order_fortran, mpi_double_precision, subgrid_basic, mpi_err)

    call mpi_type_size(mpi_double_precision, dpsize, mpi_err)
    extent = int(mpi_res*dpsize,mpi_address_kind)
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

     call send_edge(mpi_res,  mpi_grid(mpi_res, 1:mpi_res), "d", req1)
     call send_edge(mpi_res,  mpi_grid(1, 1:mpi_res), "u", req2)
     call send_edge(mpi_res,  mpi_grid(1:mpi_res, mpi_res), "r", req3)
     call send_edge(mpi_res,  mpi_grid(1:mpi_res, 1), "l", req4)
     call send_corner(mpi_grid(1, 1), "ul", req5)
     call send_corner(mpi_grid(1, mpi_res), "ur", req6)
     call send_corner(mpi_grid(mpi_res, 1), "dl", req7)
     call send_corner(mpi_grid(mpi_res, mpi_res), "dr", req8)

     call mpi_wait(req1, mpi_status_ignore, mpi_err)
     call mpi_wait(req2, mpi_status_ignore, mpi_err)
     call mpi_wait(req3, mpi_status_ignore, mpi_err)
     call mpi_wait(req4, mpi_status_ignore, mpi_err)
     call mpi_wait(req5, mpi_status_ignore, mpi_err)
     call mpi_wait(req6, mpi_status_ignore, mpi_err)
     call mpi_wait(req7, mpi_status_ignore, mpi_err)
     call mpi_wait(req8, mpi_status_ignore, mpi_err)

     call recv_edge(mpi_res,  mpi_grid(0, 1:mpi_res), "d")
     call recv_edge(mpi_res,  mpi_grid(mpi_res+1, 1:mpi_res), "u")
     call recv_edge(mpi_res,  mpi_grid(1:mpi_res, 0), "r")
     call recv_edge(mpi_res,  mpi_grid(1:mpi_res, mpi_res+1), "l")
     call recv_corner(mpi_grid(0, 0), "dr")
     call recv_corner(mpi_grid(0, mpi_res+1), "dl")
     call recv_corner(mpi_grid(mpi_res+1, 0), "ur")
     call recv_corner(mpi_grid(mpi_res+1, mpi_res+1), "ul")

    call MPI_Type_free(subgrid,mpi_err)
  end subroutine grid_scatter

  !> @brief Gathers a grid to rank 0 from mpi_grids on all processes
  !! @param[in] grid         global grid
  !! @param[inout] grid_res  global grid resolution
  !! @param[in] mpi_grid     local grid
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

    mpi_res = int(grid_res/nproc_row)

    call mpi_type_create_subarray(2, [grid_res, grid_res], [mpi_res, mpi_res], [0, 0],&
     mpi_order_fortran, mpi_double_precision, subgrid_basic, mpi_err)

    call mpi_type_size(mpi_double_precision, dpsize, mpi_err)
    extent = int(mpi_res*dpsize,mpi_address_kind)
    call mpi_type_create_resized(subgrid_basic, int(0, mpi_address_kind), extent, subgrid, mpi_err)
    call mpi_type_commit(subgrid, mpi_err)

    do rank = 0, nproc-1

      call mpi_cart_coords(cart_comm, rank, 2, rankcoords, mpi_err)

      displs(rank+1) = rankcoords(2) + nproc_row*mpi_res*rankcoords(1)
      counts = 1 ! TODO: is this right
    end do

    ! if (myrank == 0) print *, "9", counts, displs
    ! if (myrank == 0) then
    !   print *, "9.5"
    !   do i = 1,4
    !     do j = 1,4
    !       write(*, "(F5.2, 1X)", advance="no") mpi_grid(i, j)
    !     end do
    !     write(*,*)
    !   end do
    ! end if
    ! if (myrank == 0) then
    !   print *, "9.75"
    !   do i = 1,8
    !     do j = 1,8
    !       write(*, "(F5.2, 1X)", advance="no") grid(i, j)
    !     end do
    !     write(*,*)
    !   end do
    !   print *, subgrid
    ! end if
    call MPI_gatherv(mpi_grid(1:mpi_res,1:mpi_res), mpi_res*mpi_res, mpi_double_precision, &
    grid, counts, displs, subgrid, &
     0, mpi_comm_world, mpi_err)

    ! if (myrank == 0) print *, "10"
    call MPI_Type_free(subgrid,mpi_err)
    ! if (myrank == 0) print *, "11"
  end subroutine grid_gather

  !> @brief Sends a corner to the process in the corresponding direction
  !! @details Direction will always correspond to the perspective of sending proc
  !! Please call mpi_wait(req, mpi_ignore_status, mpi_err) to guarantee a sent
  !! request is complete (see scatter for an example).
  !! @param[in] val        value of the corner
  !! @param[in] direction  direction to neighbour
  !! @param[out] req       request ID
  subroutine send_corner(val, direction, req)
    real(dp), intent(in) :: val
    character(2), intent(in) :: direction
    integer, intent(out) :: req

    integer :: dir(2), recv_coords(2), recv_rank
    integer :: err
    dir = [0, 0]

    select case(direction)

    case("ul")
      dir = [-1, -1]

    case("ur")
      dir = [-1, 1]

    case("dl")
      dir = [1, -1]

    case("dr")
      dir = [1, 1]

    case default
      call logger%fatal("recv_edge", "invalid direction: not in ul, ur, dl, dr")
      call MPI_Abort(MPI_COMM_WORLD, 1, err)

    end select

    recv_coords(1) = mod(mycoords(1) + dir(1), nproc_row)
    recv_coords(2) = mod(mycoords(2) + dir(2), nproc_row)

    call mpi_cart_rank(cart_comm, recv_coords, recv_rank, mpi_err)

    call mpi_isend(val, 1, mpi_double_precision, &
    recv_rank, 1003, cart_comm, req, mpi_err)
  end subroutine


  !> @brief Receives a corner from the process in the corresponding direction
  !! @details Direction will always correspond to the perspective of sending proc
  !! Please call mpi_wait(req, mpi_ignore_status, mpi_err) to guarantee a sent
  !! request is complete (see scatter for an example).
  !! @param[out] val       value of the corner
  !! @param[in] direction  direction to neighbour
  subroutine recv_corner(val, direction)
    real(dp), intent(out) :: val
    character(2), intent(in) :: direction

    integer :: dir(2), source_coords(2), source_rank
    integer :: err
    dir = [0, 0]

    select case(direction)

    case("ul")
      dir = [-1, -1]

    case("ur")
      dir = [-1, 1]

    case("dl")
      dir = [1, -1]

    case("dr")
      dir = [1, 1]

    case default
      call logger%fatal("recv_edge", "invalid direction: not in ul, ur, dl, dr")
      call MPI_Abort(MPI_COMM_WORLD, 1, err)

    end select

    source_coords(1) = mod(mycoords(1) - dir(1), nproc_row)
    source_coords(2) = mod(mycoords(2) - dir(2), nproc_row)

    call mpi_cart_rank(cart_comm, source_coords, source_rank, mpi_err)

    call mpi_recv(val, 1, mpi_double_precision, &
    source_rank, 1003, cart_comm, mpi_status_ignore, mpi_err)
  end subroutine

  !> @brief Sends an edge to the process in the corresponding direction
  !! @details Direction will always correspond to the perspective of sending proc
  !! Please call mpi_wait(req, mpi_ignore_status, mpi_err) to guarantee a sent
  !! request is complete (see scatter for an example).
  !! @param[in] n          length of the edge
  !! @param[in] edge       values along the edge
  !! @param[in] direction  direction to neighbour
  !! @param[out] req       request ID
  subroutine send_edge(n, edge, direction, req)
    integer, intent(in) :: n
    real(dp), dimension(n), intent(in) :: edge
    character, intent(in) :: direction
    integer, intent(out) :: req

    integer :: dir_int
    integer :: err

    dir_int = 0

    select case(direction)

    case("l")
      dir_int = 4

    case("r")
      dir_int = 3

    case("u")
      dir_int = 1

    case("d")
      dir_int = 2

    case default
      call logger%fatal("recv_edge", "invalid direction: not in u, d, l, r")
      call MPI_Abort(MPI_COMM_WORLD, 1, err)

    end select

    call mpi_isend(edge, n, mpi_double_precision, &
     neigh(dir_int), 1003, cart_comm, req, mpi_err)
  end subroutine send_edge

  !> @brief Receives an edge from the process in the corresponding direction
  !! @details Direction will always correspond to the perspective of sending proc
  !! Please call mpi_wait(req, mpi_ignore_status, mpi_err) to guarantee a sent
  !! request is complete (see scatter for an example).
  !! @param[in] n          length of the edge
  !! @param[in] edge       values along the edge
  !! @param[in] direction  direction to neighbour
  subroutine recv_edge(n, edge, direction)
    integer, intent(in) :: n
    real(dp), dimension(n), intent(out) :: edge
    character, intent(in) :: direction

    integer :: dir_int
    integer :: err
    dir_int = 0

    select case(direction)

    case("l")
      dir_int = 3

    case("r")
      dir_int = 4

    case("u")
      dir_int = 2

    case("d")
      dir_int = 1

    case default
      call logger%fatal("recv_edge", "invalid direction: not in u, d, l, r")
      call MPI_Abort(MPI_COMM_WORLD, 1, err)

    end select

    call mpi_recv(edge, n, mpi_double_precision, &
     neigh(dir_int), 1003, cart_comm, mpi_status_ignore, mpi_err)
  end subroutine recv_edge

  !> @brief Shuts down MPI communications
  subroutine comms_final()
    call mpi_finalize(mpi_err)
  end subroutine comms_final
end module comms
