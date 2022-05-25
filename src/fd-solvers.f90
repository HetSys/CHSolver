module multigrid
  use globals

  implicit none
  save

  !> @brief multigrid internal grid type
  !!
  !! This is just a wrapper around an array of reals. An array of these forms a
  !! multigrid
  type :: t_grid
    real(dp), dimension(:,:), allocatable :: grid
  end type

  contains


  !> @brief Allocates storage for a multigrid with a given level
  !!
  !! @param mg[inout]  multigrid
  !! @param level[in]  maximum level
  subroutine multigrid_alloc(mg, level)
    implicit none
    type(t_grid), dimension(:), allocatable, intent(inout) :: mg
    integer, intent(in) :: level
    integer :: i

    allocate(mg(1:level))

    do i=1,level
      allocate(mg(i)%grid(0:(1+2**i),0:(1+2**i)))
    enddo
  end subroutine multigrid_alloc


  !> @brief Deallocates storage for a multigrid with a given level
  !!
  !! @param mg[inout]  multigrid
  !! @param level[in]  maximum level
  subroutine multigrid_dealloc(mg, level)
    implicit none
    type(t_grid), dimension(:), allocatable, intent(inout) :: mg
    integer, intent(in) :: level
    integer :: i

    do i=1,level
      deallocate(mg(i)%grid)
    enddo

    deallocate(mg)
  end subroutine multigrid_dealloc
end module multigrid


module fd_solvers
  use globals
  use multigrid
  use solver_utils
  use hdf5_io
  use comms

  implicit none
  save

  ! solver parameters
  real(dp), parameter :: tau = 2.0_dp
  integer, parameter :: niter = 11
  integer, parameter :: local_min = 2

  private tau, niter, invert, laplacian, compute_g, vcycle, smooth, restrict, prolongate
  public solver_ufds2t2

  contains


  !> @brief solves the dimensionless CH equation
  !!
  !! @param[in] Tout  output times
  !! @param[in] c0    initial concentration
  !! @param[in] eps2  dimensionless number
  subroutine solver_ufds2t2(t_0, Tout, CH_params, c0, eps2, errors, c1, dt_in)
    implicit none
    real(dp), intent(in) :: t_0
    real(dp), intent(in) :: Tout(:)
    real(dp), intent(in) :: eps2
    real(dp), intent(in), dimension(6) :: CH_params
    real(dp), dimension(:,:), allocatable, intent(in) :: c0
    real(dp), dimension(:,:), allocatable, intent(in), optional :: c1
    real(dp), intent(in), optional :: dt_in
    integer :: errors

    integer :: N ! grid size
    integer :: level ! grid level
    real(dp) :: dx ! grid spacing
    real(dp) :: dt, dt_min, dt_max ! timesteps
    real(dp) :: dt0, dt1, dt_out ! timesteps
    real(dp) :: a0, a1, a2, b0, b1 ! time constants
    real(dp) :: t, t_out,  tmax
    real(dp) :: t0, t1 ! adaptive timestep
    real(dp), dimension(2,2) :: A ! smoothing matrix
    real(dp) :: eps
    integer :: it, i ! iterators
    logical :: outflag
    character(len=32) :: msg ! logging message
    integer :: req1, req2, req3, req4 ! MPI comms
    logical :: double_start

    ! grid storage (local)
    real(dp), dimension(:,:), allocatable :: phi, psi, g, b, phi_prev, g_prev, work
    type(t_grid), dimension(:), allocatable :: E1, E2, R1, R2
    real(dp), dimension(:,:), allocatable :: c

    ! grid storage (global - rank 0 only)
    real(dp), dimension(:,:), allocatable :: phi_global, psi_global
    type(t_grid), dimension(:), allocatable :: E1_global, E2_global, R1_global, R2_global
    real(dp), dimension(:,:), allocatable :: c_global, c_prev_global
    integer :: N_global ! grid size
    integer :: level_global ! grid level


    ! ======================================================================== !
    !   SETUP                                                                  !
    ! ======================================================================== !
    ! set variables
    N = size(c0,1) / nproc_row
    call ilog2(N,level)
    dx = 1.0_dp/(real(N*nproc_row,dp))
    dt_min = 2.5_dp * eps2
    dt_max = dx
    dt = dt_min
    tmax = maxval(tout)
    t = t_0
    it = 1
    eps = sqrt(eps2)
    t0 = 10.0_dp * eps + t_0
    t1 = 20.0_dp * eps + t_0
    double_start = .false.

    call logger%info("solver_ufds2t2", "eps:"//to_string(eps))

    ! set global variables
    call ilog2(nproc_row,level_global)
    level_global = level_global + level
    N_global = N * nproc_row

    24 format(A, F7.3) ! output message

    ! allocate storage
    allocate(phi(0:N+1,0:N+1))
    allocate(psi(0:N+1,0:N+1))
    allocate(g(0:N+1,0:N+1))
    allocate(b(0:N+1,0:N+1))
    allocate(c(0:N+1,0:N+1))
    allocate(phi_prev(0:N+1,0:N+1))
    allocate(g_prev(0:N+1,0:N+1))
    allocate(work(0:N+1,0:N+1))


    ! allocate multigrid storage
    call multigrid_alloc(E1, level)
    call multigrid_alloc(E2, level)
    call multigrid_alloc(R1, level)
    call multigrid_alloc(R2, level)

    
    ! allocate global storage
    allocate(phi_global(N_global,N_global))
    allocate(psi_global(N_global,N_global))
    allocate(c_global(N_global,N_global))
    allocate(c_prev_global(N_global,N_global))
    call multigrid_alloc(E1_global, level_global-level+local_min)
    call multigrid_alloc(E2_global, level_global-level+local_min)
    call multigrid_alloc(R1_global, level_global-level+local_min)
    call multigrid_alloc(R2_global, level_global-level+local_min)

    ! check if optionals are present
    if (present(c1)) then
      if (.not. present(dt_in)) then
        call logger%warning("solver_ufds2t2", "no timestep provided, defaulting to first order")
      else
        double_start = .true.
      endif
    else
      if (present(dt_in)) then
        call logger%warning("solver_ufds2t2", "no paired concentration provided, defaulting to first order")
      endif
    endif

    ! send c0 (global) to phi (local)
    if (nproc > 1) then
        call grid_scatter(c0, N_global, phi, N)
    else
      phi(1:N,1:N) = c0
    end if

    ! set coupled variable
    call laplacian(phi, psi, dx, N)
    psi = tau*phi - eps2*psi

    ! swap all four edges of psi
    if (nproc > 1) then
      call send_edge(n, psi(1,1:N), "u", req1)
      call send_edge(n, psi(N,1:N), "d", req2)
      call send_edge(n, psi(1:N,1), "l", req3)
      call send_edge(n, psi(1:N,N), "r", req4)
      call mpi_wait(req1, mpi_status_ignore, mpi_err)
      call mpi_wait(req2, mpi_status_ignore, mpi_err)
      call mpi_wait(req3, mpi_status_ignore, mpi_err)
      call mpi_wait(req4, mpi_status_ignore, mpi_err)
      call recv_edge(n, psi(N+1,1:N), "u")
      call recv_edge(n, psi(0,1:N), "d")
      call recv_edge(n, psi(1:N,N+1), "l")
      call recv_edge(n, psi(1:N,0), "r")
    end if

    ! output if required
    if (tout(it)-t_0 < epsilon(tout(it))) then
      ! gather grid to rank 0
      if (nproc > 1) then
        call grid_gather(c_global, N_global, phi)
        call grid_gather(c_prev_global, N_global, phi_prev)
      else
        c_global = phi(1:N,1:N)
        c_prev_global = phi_prev(1:N,1:N)
      end if

      if (myrank == 0) then
        write(msg, 24) "Initial output at t=  0.000"
        call logger%info("solver_ufds2t2", msg)
        dt_out = dt
        t_out = t
        call dimensionalise(CH_params, c_global, t_out)
        call dimensionalise(CH_params, c_prev_global, dt_out)

        call write_to_traj(c_global, c_prev_global, t_out, dt_out, errors)
      endif

      it = it + 1
    endif


    ! ======================================================================== !
    !   FIRST TIMESTEP (first order)                                           !
    ! ======================================================================== !
    if (.not. double_start) then
      ! restrict timestep if we would otherwise exceed an output time
      if (t + dt + epsilon(t) > Tout(it)) then
        dt = tout(it) - t
        outflag = .true.
      else
        outflag = .false.
      endif
      t = t + dt
      dt0 = dt ! store current timestep

      ! compute RHS
      call compute_g(g, phi, dx, N, work)
      b = phi/dt + g

      ! store current variables
      phi_prev = phi
      g_prev = g

      ! set first order smoothing matrix
      A(1,1) = 1.0_dp/dt
      A(1,2) = 4.0_dp/(dx*dx)
      A(2,1) = -(tau+4.0_dp*eps2/(dx*dx))
      A(2,2) = 1.0_dp
      call invert(A)

      ! solve system with repeated v-cycles
      do i=1,niter
        ! compute residuals TODO: maybe use MPI
        call laplacian(psi, R1(level)%grid, dx, n)
        R1(level)%grid = R1(level)%grid - phi/dt + b
        call laplacian(phi, R2(level)%grid, dx, n)
        R2(level)%grid = tau*phi-eps2*R2(level)%grid - psi

        E1(level)%grid = 0.0_dp
        E2(level)%grid = 0.0_dp

        ! TODO: finish iteration conditional on the size of the residual

        ! perform a single v-cycle
        if (nproc > 1) then
          call vcycle(A, E1, E2, R1, R2, E1_global, E2_global, R1_global, R2_global, &
                      eps2, N, dx, level)
        else
          call vcycle0(A, E1, E2, R1, R2, eps2, N, dx, level)
        endif

        ! update with errors
        phi = phi + E1(level)%grid
        psi = psi + E2(level)%grid

        ! send edges between ranks
        if (nproc > 1) then
          call send_edge(n, phi(1,1:N), "u", req1)
          call send_edge(n, phi(N,1:N), "d", req2)
          call send_edge(n, phi(1:N,1), "l", req3)
          call send_edge(n, phi(1:N,N), "r", req4)

          call mpi_wait(req1, mpi_status_ignore, mpi_err)
          call mpi_wait(req2, mpi_status_ignore, mpi_err)
          call mpi_wait(req3, mpi_status_ignore, mpi_err)
          call mpi_wait(req4, mpi_status_ignore, mpi_err)

          call recv_edge(n, phi(N+1,1:N), "u")
          call recv_edge(n, phi(0,1:N), "d")
          call recv_edge(n, phi(1:N,N+1), "l")
          call recv_edge(n, phi(1:N,0), "r")

          call send_edge(n, psi(1,1:N), "u", req1)
          call send_edge(n, psi(N,1:N), "d", req2)
          call send_edge(n, psi(1:N,1), "l", req3)
          call send_edge(n, psi(1:N,N), "r", req4)

          call mpi_wait(req1, mpi_status_ignore, mpi_err)
          call mpi_wait(req2, mpi_status_ignore, mpi_err)
          call mpi_wait(req3, mpi_status_ignore, mpi_err)
          call mpi_wait(req4, mpi_status_ignore, mpi_err)

          call recv_edge(n, psi(N+1,1:N), "u")
          call recv_edge(n, psi(0,1:N), "d")
          call recv_edge(n, psi(1:N,N+1), "l")
          call recv_edge(n, psi(1:N,0), "r")
        endif
      enddo

      ! output if required
      if (outflag) then
        ! gather grid to rank 0
        if (nproc > 1) then
          call grid_gather(c_global, N_global, phi)
          call grid_gather(c_prev_global, N_global, phi_prev)
        else
          c_global = phi(1:N,1:N)
          c_prev_global = phi_prev(1:N,1:N)
        end if

        if (myrank == 0) then
          write(msg, 24) "Output at t=", t
          call logger%info("solver_ufds2t2", msg)
          dt_out = dt
          t_out = t
          ! c = phi(1:N,1:N)
          ! c_prev = phi_prev(1:N,1:N)
          call dimensionalise(CH_params, c_global, t_out)
          call dimensionalise(CH_params, c_prev_global, dt_out)

          call write_to_traj(c_global, c_prev_global, t_out, dt_out, errors)
        endif

        it = it + 1
      endif
    else
      ! compute RHS
      call compute_g(g, phi, dx, N, work)

      ! move c0 to previous timestep
      phi_prev = phi
      g_prev = g

      !! REAPEAT FOR C1
      ! send c0 (global) to phi (local)
      if (nproc > 1) then
        call grid_scatter(c1, N_global, phi, N)
      else
        phi(1:N,1:N) = c1
      end if

      ! set coupled variable
      call laplacian(phi, psi, dx, N)
      psi = tau*phi - eps2*psi

      ! swap all four edges of psi
      if (nproc > 1) then
        call send_edge(n, psi(1,1:N), "u", req1)
        call send_edge(n, psi(N,1:N), "d", req2)
        call send_edge(n, psi(1:N,1), "l", req3)
        call send_edge(n, psi(1:N,N), "r", req4)
        call mpi_wait(req1, mpi_status_ignore, mpi_err)
        call mpi_wait(req2, mpi_status_ignore, mpi_err)
        call mpi_wait(req3, mpi_status_ignore, mpi_err)
        call mpi_wait(req4, mpi_status_ignore, mpi_err)
        call recv_edge(n, psi(N+1,1:N), "u")
        call recv_edge(n, psi(0,1:N), "d")
        call recv_edge(n, psi(1:N,N+1), "l")
        call recv_edge(n, psi(1:N,0), "r")
      end if

      ! set timesteps
      dt0 = dt_in
      dt_min = dt_in
    endif


    ! ======================================================================== !
    !   REMAINING TIMESTEPS (second order)                                     !
    ! ======================================================================== !
    do while (t < tmax)
      ! set current timestep TODO: condition on curvature rather than time
      if (t < t0) then
        dt = dt_min
      else if (t < t1) then
        dt = ((t1-t) * dt_min + (t-t0) * dt_max)/(t1-t0)
      else
        dt = dt_max
      endif

      ! restrict timestep if we would otherwise exceed an output time
      if (t + dt + epsilon(t) >= tout(it)) then
        dt = tout(it) - t
        outflag = .true.
      else
        outflag = .false.
      endif
      t = t + dt

      ! timestep variables
      dt1 = dt+dt0
      a0 = dt*dt/(dt0*dt1)
      a1 = -dt1/dt0
      a2 = (dt0+dt+dt)/dt1
      b0 = -dt/dt0
      b1 = dt1/dt0
      dt0 = dt ! store current timestep

      ! compute RHS
      call compute_g(g, phi, dx, N, work)
      b = -a1/dt*phi - a0/dt*phi_prev + b1*g + b0*g_prev

      ! store current variables
      phi_prev = phi
      g_prev = g

      ! set second order smoothing matrix
      A(1,1) = a2/dt
      A(1,2) = 4.0_dp/(dx*dx)
      A(2,1) = -(tau+4.0_dp*eps2/(dx*dx))
      A(2,2) = 1.0_dp
      call invert(A)

      ! solve system with repeated v-cycles
      do i=1,niter
        ! compute residuals
        call laplacian(psi, R1(level)%grid, dx, n)
        R1(level)%grid = R1(level)%grid - a2/dt*phi + b
        call laplacian(phi, R2(level)%grid, dx, n)
        R2(level)%grid = tau*phi-eps2*R2(level)%grid - psi

        E1(level)%grid = 0.0_dp
        E2(level)%grid = 0.0_dp

        ! TODO: finish iteration conditional on the size of the residual

        ! perform a single v-cycle
        if (nproc > 1) then
          call vcycle(A, E1, E2, R1, R2, E1_global, E2_global, R1_global, R2_global, &
                      eps2, N, dx, level)
        else
          call vcycle0(A, E1, E2, R1, R2, eps2, N, dx, level)
        endif

        ! update with errors
        phi = phi + E1(level)%grid
        psi = psi + E2(level)%grid

        ! send edges between ranks
        if (nproc > 1) then
          call send_edge(n, phi(1,1:N), "u", req1)
          call send_edge(n, phi(N,1:N), "d", req2)
          call send_edge(n, phi(1:N,1), "l", req3)
          call send_edge(n, phi(1:N,N), "r", req4)

          call mpi_wait(req1, mpi_status_ignore, mpi_err)
          call mpi_wait(req2, mpi_status_ignore, mpi_err)
          call mpi_wait(req3, mpi_status_ignore, mpi_err)
          call mpi_wait(req4, mpi_status_ignore, mpi_err)

          call recv_edge(n, phi(N+1,1:N), "u")
          call recv_edge(n, phi(0,1:N), "d")
          call recv_edge(n, phi(1:N,N+1), "l")
          call recv_edge(n, phi(1:N,0), "r")

          call send_edge(n, psi(1,1:N), "u", req1)
          call send_edge(n, psi(N,1:N), "d", req2)
          call send_edge(n, psi(1:N,1), "l", req3)
          call send_edge(n, psi(1:N,N), "r", req4)

          call mpi_wait(req1, mpi_status_ignore, mpi_err)
          call mpi_wait(req2, mpi_status_ignore, mpi_err)
          call mpi_wait(req3, mpi_status_ignore, mpi_err)
          call mpi_wait(req4, mpi_status_ignore, mpi_err)

          call recv_edge(n, psi(N+1,1:N), "u")
          call recv_edge(n, psi(0,1:N), "d")
          call recv_edge(n, psi(1:N,N+1), "l")
          call recv_edge(n, psi(1:N,0), "r")
        endif

        ! print size of residual
        ! print *, i, "r1 max = ", maxval(abs(R1(level)%grid))
      enddo

      ! output if required
      if (outflag) then
        ! gather grid to rank 0
        if (nproc > 1) then
          call grid_gather(c_global, N_global, phi)
          call grid_gather(c_prev_global, N_global, phi_prev)
        else
          c_global = phi(1:N,1:N)
          c_prev_global = phi_prev(1:N,1:N)
        end if

        if (myrank == 0) then
          write(msg, 24) "Output at t=", t
          call logger%info("solver_ufds2t2", msg)
          dt_out = dt
          t_out = t
          call dimensionalise(CH_params, c_global, t_out)
          call dimensionalise(CH_params, c_prev_global, dt_out)

          call write_to_traj(c_global, c_prev_global, t_out, dt_out, errors)
        endif

        it = it + 1
      endif
    enddo


    ! ========================================================================== !
    !   CLEAN UP                                                                 !
    ! ========================================================================== !
    ! deallocate storage
    deallocate(phi)
    deallocate(psi)
    deallocate(g)
    deallocate(b)
    deallocate(c)
    deallocate(phi_prev)
    deallocate(g_prev)
    deallocate(work)

    ! deallocate multigrid storage
    call multigrid_dealloc(E1, level)
    call multigrid_dealloc(E2, level)
    call multigrid_dealloc(R1, level)
    call multigrid_dealloc(R2, level)

    ! deallocate global storage
    deallocate(phi_global)
    deallocate(psi_global)
    deallocate(c_global)
    deallocate(c_prev_global)
    call multigrid_dealloc(E1_global, level_global-level+local_min)
    call multigrid_dealloc(E2_global, level_global-level+local_min)
    call multigrid_dealloc(R1_global, level_global-level+local_min)
    call multigrid_dealloc(R2_global, level_global-level+local_min)
  end subroutine solver_ufds2t2


  !> @brief inverts a non-singular 2x2 matrix A in place
  subroutine invert(A)
    implicit none
    real(dp), dimension(2,2), intent(inout) :: A
    real(dp) :: det, tmp

    det = A(1,1)*A(2,2) - A(1,2)*A(2,1)
    tmp = A(1,1)

    A(1,1) = A(2,2)
    A(1,2) = -A(1,2)
    A(2,1) = -A(2,1)
    A(2,2) = tmp

    A = A / det
  end subroutine invert


  !> @brief computes the periodic Laplacian for a grid in vector form
  !!
  !! @param x   input grid
  !! @param y   output grid, y = Lx
  !! @param dx  grid spacing
  !! @param n   grid size
  subroutine laplacian(x, y, dx, n)
    implicit none
    integer, intent(in) :: n
    real(dp), dimension(0:n+1,0:n+1), intent(in) :: x
    real(dp), dimension(0:n+1,0:n+1), intent(inout) :: y
    real(dp), intent(in) :: dx
    integer :: i, j
    real(dp) :: dx2_ ! interim constants

    dx2_ = 1.0_dp / (dx*dx)

    ! interior
    do j=1,n
      do i=1,n
        y(i,j) = dx2_*(x(i+1,j) + x(i-1,j) + x(i,j+1) + x(i,j-1) - 4*x(i,j))
      enddo
    enddo
  end subroutine laplacian


  !> @brief computes the periodic Laplacian for a grid in vector form
  !!
  !! @param g    output g vector
  !! @param phi  input phi vector
  !! @param dx   grid spacing
  !! @param tau  solver parameter
  !! @param n    grid size
  !! @param work allocated work vector (same size as g)
  subroutine compute_g(g, phi, dx, n, work)
    implicit none
    integer, intent(in) :: n
    real(dp), dimension(0:n+1,0:n+1), intent(in) :: phi
    real(dp), dimension(0:n+1,0:n+1), intent(inout) :: work
    real(dp), dimension(0:n+1,0:n+1), intent(inout) :: g
    real(dp), intent(in) :: dx

    work = phi * (phi*phi - (1+tau))

    call laplacian(work, g, dx, n)
  end subroutine compute_g


  !> @brief performs a single iteration of the multigrid v-cycle
  !!
  !! @param A      2x2 smoothing matrix (solution matrix for a 1x1 system)
  !! @param E1     error multigrid for the first variable
  !! @param E2     error multigrid for the second variable
  !! @param R1     residual multigrid for the first variable
  !! @param R2     residual multigrid for the second variable
  !! @param eps2   PDE parameter
  !! @param N      grid size
  !! @param dx     grid spacing
  !! @param level  grid level
  subroutine vcycle(A, E1, E2, R1, R2, E1g, E2g, R1g, R2g, eps2, N, dx, level)
    implicit none
    integer, intent(in) :: N, level
    real(dp), dimension(2,2), intent(in) :: A
    type(t_grid), dimension(:), intent(inout) :: E1, E2, R1, R2
    type(t_grid), dimension(:), intent(inout) :: E1g, E2g, R1g, R2g
    real(dp), intent(in) :: eps2, dx

    integer :: l, nl, l_global, nl_global
    real(dp) :: dxl

    nl = n
    dxl = dx

    ! go up, smoothing and restricting
    do l=level,(local_min+1),-1
      E1(l)%grid = 0.0_dp
      E2(l)%grid = 0.0_dp

      call smooth(A, E1(l)%grid, E2(l)%grid, R1(l)%grid, R2(l)%grid, eps2, nl, dxl)
      call restrict(R1(l)%grid, R2(l)%grid, R1(l-1)%grid, R2(l-1)%grid, nl)

      nl = nl/2;
      dxl = dxl * 2.0_dp
    enddo

    ! ! smooth at level local_min TODO: remove bad workaround
    ! E1(l)%grid = 0.0_dp
    ! E2(l)%grid = 0.0_dp
    ! do i=1,5
    !   call smooth(A, E1(local_min)%grid, E2(local_min)%grid, R1(local_min)%grid, R2(local_min)%grid, &
    !               eps2, nl, dxl)
    ! enddo

    ! gather onto rank 0
    call ilog2(nproc_row, l_global)
    l_global = l_global + l
    nl_global = nl*nproc_row
    call grid_gather(E1g(l_global)%grid(1:nl_global,1:nl_global), nl_global, E1(l)%grid)
    call grid_gather(E2g(l_global)%grid(1:nl_global,1:nl_global), nl_global, E2(l)%grid)
    call grid_gather(R1g(l_global)%grid(1:nl_global,1:nl_global), nl_global, R1(l)%grid)
    call grid_gather(R2g(l_global)%grid(1:nl_global,1:nl_global), nl_global, R2(l)%grid)

    ! final vcycles on rank 0
    if (myrank == 0) then
      ! print *, nproc, nproc_row, l_global
      call vcycle0(A, E1g, E2g, R1g, R2g, eps2, nl*nproc_row, dxl, l_global)
    endif

    ! scatter back to all ranks
    call grid_scatter(E1g(l_global)%grid(1:nl_global,1:nl_global), nl_global, E1(l)%grid, nl)
    call grid_scatter(E2g(l_global)%grid(1:nl_global,1:nl_global), nl_global, E2(l)%grid, nl)
    call grid_scatter(R1g(l_global)%grid(1:nl_global,1:nl_global), nl_global, R1(l)%grid, nl)
    call grid_scatter(R2g(l_global)%grid(1:nl_global,1:nl_global), nl_global, R2(l)%grid, nl)

    ! go down, smoothing and prolongating
    do l=local_min,level-1
      call prolongate(E1(l+1)%grid, E2(l+1)%grid, E1(l)%grid, E2(l)%grid, nl)

      nl = nl*2;
      dxl = dxl * 0.5_dp

      call smooth(A, E1(l+1)%grid, E2(l+1)%grid, R1(l+1)%grid, R2(l+1)%grid, eps2, nl, dxl)
    enddo
  end subroutine vcycle


  !> @brief performs a single red/black smooth
  !!
  !! @param A      2x2 smoothing matrix (inverse solution matrix for a 1x1 system)
  !! @param E1     error array for the first variable at level
  !! @param E2     error array for the second variable at level
  !! @param R1     residual array for the first variable at level
  !! @param R2     residual array for the second variable at level
  !! @param eps2   PDE parameter
  !! @param N      grid size
  !! @param dx     grid spacing
  subroutine smooth(A, E1, E2, R1, R2, eps2, N, dx)
    implicit none
    integer, intent(in) :: n
    real(dp), dimension(2,2), intent(in) :: A
    real(dp), dimension(0:n+1,0:n+1), intent(inout) :: E1, E2, R1, R2
    real(dp), intent(in) :: eps2, dx

    real(dp), dimension(2) :: rhs
    real(dp) :: dx2_
    integer :: i, j, shift
    integer :: req1, req2, req3, req4

    dx2_ = 1.0_dp / (dx*dx)

    ! send edges
    call send_edge(n, E1(1,1:N), "u", req1)
    call send_edge(n, E1(N,1:N), "d", req2)
    call send_edge(n, E1(1:N,1), "l", req3)
    call send_edge(n, E1(1:N,N), "r", req4)

    call send_edge(n, E2(1,1:N), "u", req1)
    call send_edge(n, E2(N,1:N), "d", req2)
    call send_edge(n, E2(1:N,1), "l", req3)
    call send_edge(n, E2(1:N,N), "r", req4)


    ! ================ !
    ! SMOOTH RED NODES !
    ! ================ !
    ! interior (while the edges are sent)
    do j=2,n-1
      shift = mod(j,2)
      do i=2+shift,n-1+shift,2
        ! compute RHS
        rhs(1) = R1(i,j) + dx2_ * (E2(i+1,j) + E2(i-1,j) + E2(i,j+1) + E2(i,j-1))
        rhs(2) = R2(i,j) - eps2*dx2_ * (E1(i+1,j) + E1(i-1,j) + E1(i,j+1) + E1(i,j-1))

        ! solve for new errors
        E1(i,j) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
        E2(i,j) = A(2,1)*rhs(1) + A(2,2)*rhs(2)
      enddo
    enddo

    ! receive edges
    call mpi_wait(req1, mpi_status_ignore, mpi_err)
    call mpi_wait(req2, mpi_status_ignore, mpi_err)
    call mpi_wait(req3, mpi_status_ignore, mpi_err)
    call mpi_wait(req4, mpi_status_ignore, mpi_err)
    call recv_edge(n, E1(N+1,1:N), "u")
    call recv_edge(n, E1(0,1:N), "d")
    call recv_edge(n, E1(1:N,N+1), "l")
    call recv_edge(n, E1(1:N,0), "r")

    call mpi_wait(req1, mpi_status_ignore, mpi_err)
    call mpi_wait(req2, mpi_status_ignore, mpi_err)
    call mpi_wait(req3, mpi_status_ignore, mpi_err)
    call mpi_wait(req4, mpi_status_ignore, mpi_err)
    call recv_edge(n, E2(N+1,1:N), "u")
    call recv_edge(n, E2(0,1:N), "d")
    call recv_edge(n, E2(1:N,N+1), "l")
    call recv_edge(n, E2(1:N,0), "r")

    ! left edge
    do i=3,n,2
      rhs(1) = R1(i,1) + dx2_ * (E2(i+1,1) + E2(i-1,1) + E2(i,2) + E2(i,0))
      rhs(2) = R2(i,1) - eps2*dx2_ * (E1(i+1,1) + E1(i-1,1) + E1(i,2) + E1(i,0))
      E1(i,1) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
      E2(i,1) = A(2,1)*rhs(1) + A(2,2)*rhs(2)
    enddo

    ! right edge
    do i=2,n,2
      rhs(1) = R1(i,n) + dx2_ * (E2(i+1,n) + E2(i-1,n) + E2(i,n+1) + E2(i,n-1))
      rhs(2) = R2(i,n) - eps2*dx2_ * (E1(i+1,n) + E1(i-1,n) + E1(i,n+1) + E1(i,n-1))
      E1(i,n) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
      E2(i,n) = A(2,1)*rhs(1) + A(2,2)*rhs(2)
    enddo

    ! top edge
    do j=1,n,2
      rhs(1) = R1(1,j) + dx2_ * (E2(2,j) + E2(0,j) + E2(1,j+1) + E2(1,j-1))
      rhs(2) = R2(1,j) - eps2*dx2_ * (E1(2,j) + E1(0,j) + E1(1,j+1) + E1(1,j-1))
      E1(1,j) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
      E2(1,j) = A(2,1)*rhs(1) + A(2,2)*rhs(2)
    enddo

    ! bottom edge
    do j=2,n,2
      rhs(1) = R1(n,j) + dx2_ * (E2(n+1,j) + E2(n-1,j) + E2(n,j+1) + E2(n,j-1))
      rhs(2) = R2(n,j) - eps2*dx2_ * (E1(n+1,j) + E1(n-1,j) + E1(n,j+1) + E1(n,j-1))
      E1(n,j) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
      E2(n,j) = A(2,1)*rhs(1) + A(2,2)*rhs(2)
    enddo


    ! ================== !
    ! SMOOTH BLACK NODES !
    ! ================== !
    ! interior and edges
    do j=1,n
      shift = mod(j,2)
      do i=1+shift,n+shift,2
        ! compute RHS
        rhs(1) = R1(i,j) + dx2_ * (E2(i+1,j) + E2(i-1,j) + E2(i,j+1) + E2(i,j-1))
        rhs(2) = R2(i,j) - eps2*dx2_ * (E1(i+1,j) + E1(i-1,j) + E1(i,j+1) + E1(i,j-1))

        ! solve for new errors
        E1(i,j) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
        E2(i,j) = A(2,1)*rhs(1) + A(2,2)*rhs(2)
      enddo
    enddo
  end subroutine smooth


  !> @brief restricts (ie coarsens) to level-1 from level
  !!
  !! @param R1     residual multigrid for the first variable
  !! @param R2     residual multigrid for the second variable
  !! @param N      grid size
  subroutine restrict(R1f, R2f, R1c, R2c, N)
    implicit none
    integer, intent(in) :: n
    real(dp), dimension(0:n+1,0:n+1), intent(inout) :: R1f, R2f
    real(dp), dimension(0:n/2+1,0:n/2+1), intent(inout) :: R1c, R2c

    integer :: i, j, im, jm, nc
    nc = n/2

    ! use a simple average across the grid
    do j=1,nc
      jm = 2*j-1
      do i=1,nc
        im = 2*i-1

        ! fine -> coarse
        R1c(i,j) = 0.25_dp*(R1f(im,jm) + R1f(im+1,jm) + R1f(im,jm+1) + R1f(im+1,jm+1))
        R2c(i,j) = 0.25_dp*(R2f(im,jm) + R2f(im+1,jm) + R2f(im,jm+1) + R2f(im+1,jm+1))
      enddo
    enddo
  end subroutine restrict


  !> @brief prolongates (ie refines) from level to level+1 (bilinear interpolation)
  !!
  !! @param E1     error multigrid for the first variable
  !! @param E2     error multigrid for the second variable
  !! @param N      grid size
  subroutine prolongate(E1f, E2f, E1c, E2c, N)
    implicit none
    integer, intent(in) :: n
    real(dp), dimension(0:n+1,0:n+1), intent(inout) :: E1c, E2c
    real(dp), dimension(0:2*n+1,0:2*n+1), intent(inout) :: E1f, E2f

    integer :: Nf
    integer :: i, j, if, jf
    real(dp), parameter :: w1 = 0.5625_dp
    real(dp), parameter :: w2 = 0.1875_dp
    real(dp), parameter :: w3 = 0.0625_dp

    integer :: req01, req02, req03, req04, req05, req06, req07, req08
    integer :: req11, req12, req13, req14, req15, req16, req17, req18

    Nf = 2*N

    ! send edges and corners
    call send_edge(Nf, E1f(1,1:Nf), "u", req01)
    call send_edge(Nf, E1f(Nf,1:Nf), "d", req02)
    call send_edge(Nf, E1f(1:Nf,1), "l", req03)
    call send_edge(Nf, E1f(1:Nf,Nf), "r", req04)
    call send_corner(E1f(1, 1), "ul", req05)
    call send_corner(E1f(1, Nf), "ur", req06)
    call send_corner(E1f(Nf, 1), "dl", req07)
    call send_corner(E1f(Nf, Nf), "dr", req08)

    call send_edge(Nf, E2f(1,1:Nf), "u", req11)
    call send_edge(Nf, E2f(Nf,1:Nf), "d", req12)
    call send_edge(Nf, E2f(1:Nf,1), "l", req13)
    call send_edge(Nf, E2f(1:Nf,Nf), "r", req14)
    call send_corner(E2f(1, 1), "ul", req15)
    call send_corner(E2f(1, Nf), "ur", req16)
    call send_corner(E2f(Nf, 1), "dl", req17)
    call send_corner(E2f(Nf, Nf), "dr", req18)

    ! interior
    do j=2,n-1
      jf = 2*j-1
      do i=2,n-1
        if = 2*i-1

        ! largest contribution to nearest
        E1f(if,jf) = E1f(if,jf) + w1*E1c(i,j)
        E1f(if+1,jf) = E1f(if+1,jf) + w1*E1c(i,j)
        E1f(if,jf+1) = E1f(if,jf+1) + w1*E1c(i,j)
        E1f(if+1,jf+1) = E1f(if+1,jf+1) + w1*E1c(i,j)

        ! lesser contribution to intermediate
        E1f(if-1,jf) = E1f(if-1,jf) + w2*E1c(i,j)
        E1f(if-1,jf+1) = E1f(if-1,jf+1) + w2*E1c(i,j)
        E1f(if,jf-1) = E1f(if,jf-1) + w2*E1c(i,j)
        E1f(if,jf+2) = E1f(if,jf+2) + w2*E1c(i,j)
        E1f(if+1,jf-1) = E1f(if+1,jf-1) + w2*E1c(i,j)
        E1f(if+1,jf+2) = E1f(if+1,jf+2) + w2*E1c(i,j)
        E1f(if+2,jf) = E1f(if+2,jf) + w2*E1c(i,j)
        E1f(if+2,jf+1) = E1f(if+2,jf+1) + w2*E1c(i,j)

        ! least contribution to furthest
        E1f(if-1,jf-1) = E1f(if-1,jf-1) + w3*E1c(i,j)
        E1f(if-1,jf+2) = E1f(if-1,jf+2) + w3*E1c(i,j)
        E1f(if+2,jf-1) = E1f(if+2,jf-1) + w3*E1c(i,j)
        E1f(if+2,jf+2) = E1f(if+2,jf+2) + w3*E1c(i,j)

        ! largest contribution to nearest
        E2f(if,jf) = E2f(if,jf) + w1*E2c(i,j)
        E2f(if+1,jf) = E2f(if+1,jf) + w1*E2c(i,j)
        E2f(if,jf+1) = E2f(if,jf+1) + w1*E2c(i,j)
        E2f(if+1,jf+1) = E2f(if+1,jf+1) + w1*E2c(i,j)

        ! lesser contribution to intermediate
        E2f(if-1,jf) = E2f(if-1,jf) + w2*E2c(i,j)
        E2f(if-1,jf+1) = E2f(if-1,jf+1) + w2*E2c(i,j)
        E2f(if,jf-1) = E2f(if,jf-1) + w2*E2c(i,j)
        E2f(if,jf+2) = E2f(if,jf+2) + w2*E2c(i,j)
        E2f(if+1,jf-1) = E2f(if+1,jf-1) + w2*E2c(i,j)
        E2f(if+1,jf+2) = E2f(if+1,jf+2) + w2*E2c(i,j)
        E2f(if+2,jf) = E2f(if+2,jf) + w2*E2c(i,j)
        E2f(if+2,jf+1) = E2f(if+2,jf+1) + w2*E2c(i,j)

        ! least contribution to furthest
        E2f(if-1,jf-1) = E2f(if-1,jf-1) + w3*E2c(i,j)
        E2f(if-1,jf+2) = E2f(if-1,jf+2) + w3*E2c(i,j)
        E2f(if+2,jf-1) = E2f(if+2,jf-1) + w3*E2c(i,j)
        E2f(if+2,jf+2) = E2f(if+2,jf+2) + w3*E2c(i,j)
      enddo
    enddo

    ! receive edges and corners
    call mpi_wait(req01, mpi_status_ignore, mpi_err)
    call mpi_wait(req02, mpi_status_ignore, mpi_err)
    call mpi_wait(req03, mpi_status_ignore, mpi_err)
    call mpi_wait(req04, mpi_status_ignore, mpi_err)
    call mpi_wait(req05, mpi_status_ignore, mpi_err)
    call mpi_wait(req06, mpi_status_ignore, mpi_err)
    call mpi_wait(req07, mpi_status_ignore, mpi_err)
    call mpi_wait(req08, mpi_status_ignore, mpi_err)
    call recv_edge(Nf, E1f(Nf+1,1:Nf), "u")
    call recv_edge(Nf, E1f(0,1:Nf), "d")
    call recv_edge(Nf, E1f(1:Nf,N+1), "l")
    call recv_edge(Nf, E1f(1:Nf,0), "r")
    call recv_corner(E1f(0, 0), "dr")
    call recv_corner(E1f(0, Nf+1), "dl")
    call recv_corner(E1f(Nf+1, 0), "ur")
    call recv_corner(E1f(Nf+1, Nf+1), "ul")

    call mpi_wait(req11, mpi_status_ignore, mpi_err)
    call mpi_wait(req12, mpi_status_ignore, mpi_err)
    call mpi_wait(req13, mpi_status_ignore, mpi_err)
    call mpi_wait(req14, mpi_status_ignore, mpi_err)
    call mpi_wait(req15, mpi_status_ignore, mpi_err)
    call mpi_wait(req16, mpi_status_ignore, mpi_err)
    call mpi_wait(req17, mpi_status_ignore, mpi_err)
    call mpi_wait(req18, mpi_status_ignore, mpi_err)
    call recv_edge(Nf, E2f(Nf+1,1:Nf), "u")
    call recv_edge(Nf, E2f(0,1:Nf), "d")
    call recv_edge(Nf, E2f(1:Nf,Nf+1), "l")
    call recv_edge(Nf, E2f(1:Nf,0), "r")
    call recv_corner(E2f(0, 0), "dr")
    call recv_corner(E2f(0, Nf+1), "dl")
    call recv_corner(E2f(Nf+1, 0), "ur")
    call recv_corner(E2f(Nf+1, Nf+1), "ul")

    ! left
    do j=2,n-1
      jf = 2*j-1
      do i=1,1
        if = 2*i-1

        ! largest contribution to nearest
        E1f(if,jf) = E1f(if,jf) + w1*E1c(i,j)
        E1f(if+1,jf) = E1f(if+1,jf) + w1*E1c(i,j)
        E1f(if,jf+1) = E1f(if,jf+1) + w1*E1c(i,j)
        E1f(if+1,jf+1) = E1f(if+1,jf+1) + w1*E1c(i,j)

        ! lesser contribution to intermediate
        E1f(if-1,jf) = E1f(if-1,jf) + w2*E1c(i,j)
        E1f(if-1,jf+1) = E1f(if-1,jf+1) + w2*E1c(i,j)
        E1f(if,jf-1) = E1f(if,jf-1) + w2*E1c(i,j)
        E1f(if,jf+2) = E1f(if,jf+2) + w2*E1c(i,j)
        E1f(if+1,jf-1) = E1f(if+1,jf-1) + w2*E1c(i,j)
        E1f(if+1,jf+2) = E1f(if+1,jf+2) + w2*E1c(i,j)
        E1f(if+2,jf) = E1f(if+2,jf) + w2*E1c(i,j)
        E1f(if+2,jf+1) = E1f(if+2,jf+1) + w2*E1c(i,j)

        ! least contribution to furthest
        E1f(if-1,jf-1) = E1f(if-1,jf-1) + w3*E1c(i,j)
        E1f(if-1,jf+2) = E1f(if-1,jf+2) + w3*E1c(i,j)
        E1f(if+2,jf-1) = E1f(if+2,jf-1) + w3*E1c(i,j)
        E1f(if+2,jf+2) = E1f(if+2,jf+2) + w3*E1c(i,j)

        ! largest contribution to nearest
        E2f(if,jf) = E2f(if,jf) + w1*E2c(i,j)
        E2f(if+1,jf) = E2f(if+1,jf) + w1*E2c(i,j)
        E2f(if,jf+1) = E2f(if,jf+1) + w1*E2c(i,j)
        E2f(if+1,jf+1) = E2f(if+1,jf+1) + w1*E2c(i,j)

        ! lesser contribution to intermediate
        E2f(if-1,jf) = E2f(if-1,jf) + w2*E2c(i,j)
        E2f(if-1,jf+1) = E2f(if-1,jf+1) + w2*E2c(i,j)
        E2f(if,jf-1) = E2f(if,jf-1) + w2*E2c(i,j)
        E2f(if,jf+2) = E2f(if,jf+2) + w2*E2c(i,j)
        E2f(if+1,jf-1) = E2f(if+1,jf-1) + w2*E2c(i,j)
        E2f(if+1,jf+2) = E2f(if+1,jf+2) + w2*E2c(i,j)
        E2f(if+2,jf) = E2f(if+2,jf) + w2*E2c(i,j)
        E2f(if+2,jf+1) = E2f(if+2,jf+1) + w2*E2c(i,j)

        ! least contribution to furthest
        E2f(if-1,jf-1) = E2f(if-1,jf-1) + w3*E2c(i,j)
        E2f(if-1,jf+2) = E2f(if-1,jf+2) + w3*E2c(i,j)
        E2f(if+2,jf-1) = E2f(if+2,jf-1) + w3*E2c(i,j)
        E2f(if+2,jf+2) = E2f(if+2,jf+2) + w3*E2c(i,j)
      enddo
    enddo

    ! right
    do j=2,n-1
      jf = 2*j-1
      do i=n,n
        if = 2*i-1

        ! largest contribution to nearest
        E1f(if,jf) = E1f(if,jf) + w1*E1c(i,j)
        E1f(if+1,jf) = E1f(if+1,jf) + w1*E1c(i,j)
        E1f(if,jf+1) = E1f(if,jf+1) + w1*E1c(i,j)
        E1f(if+1,jf+1) = E1f(if+1,jf+1) + w1*E1c(i,j)

        ! lesser contribution to intermediate
        E1f(if-1,jf) = E1f(if-1,jf) + w2*E1c(i,j)
        E1f(if-1,jf+1) = E1f(if-1,jf+1) + w2*E1c(i,j)
        E1f(if,jf-1) = E1f(if,jf-1) + w2*E1c(i,j)
        E1f(if,jf+2) = E1f(if,jf+2) + w2*E1c(i,j)
        E1f(if+1,jf-1) = E1f(if+1,jf-1) + w2*E1c(i,j)
        E1f(if+1,jf+2) = E1f(if+1,jf+2) + w2*E1c(i,j)
        E1f(if+2,jf) = E1f(if+2,jf) + w2*E1c(i,j)
        E1f(if+2,jf+1) = E1f(if+2,jf+1) + w2*E1c(i,j)

        ! least contribution to furthest
        E1f(if-1,jf-1) = E1f(if-1,jf-1) + w3*E1c(i,j)
        E1f(if-1,jf+2) = E1f(if-1,jf+2) + w3*E1c(i,j)
        E1f(if+2,jf-1) = E1f(if+2,jf-1) + w3*E1c(i,j)
        E1f(if+2,jf+2) = E1f(if+2,jf+2) + w3*E1c(i,j)

        ! largest contribution to nearest
        E2f(if,jf) = E2f(if,jf) + w1*E2c(i,j)
        E2f(if+1,jf) = E2f(if+1,jf) + w1*E2c(i,j)
        E2f(if,jf+1) = E2f(if,jf+1) + w1*E2c(i,j)
        E2f(if+1,jf+1) = E2f(if+1,jf+1) + w1*E2c(i,j)

        ! lesser contribution to intermediate
        E2f(if-1,jf) = E2f(if-1,jf) + w2*E2c(i,j)
        E2f(if-1,jf+1) = E2f(if-1,jf+1) + w2*E2c(i,j)
        E2f(if,jf-1) = E2f(if,jf-1) + w2*E2c(i,j)
        E2f(if,jf+2) = E2f(if,jf+2) + w2*E2c(i,j)
        E2f(if+1,jf-1) = E2f(if+1,jf-1) + w2*E2c(i,j)
        E2f(if+1,jf+2) = E2f(if+1,jf+2) + w2*E2c(i,j)
        E2f(if+2,jf) = E2f(if+2,jf) + w2*E2c(i,j)
        E2f(if+2,jf+1) = E2f(if+2,jf+1) + w2*E2c(i,j)

        ! least contribution to furthest
        E2f(if-1,jf-1) = E2f(if-1,jf-1) + w3*E2c(i,j)
        E2f(if-1,jf+2) = E2f(if-1,jf+2) + w3*E2c(i,j)
        E2f(if+2,jf-1) = E2f(if+2,jf-1) + w3*E2c(i,j)
        E2f(if+2,jf+2) = E2f(if+2,jf+2) + w3*E2c(i,j)
      enddo
    enddo

    ! top
    do j=1,1
      jf = 2*j-1
      do i=1,n
        if = 2*i-1

        ! largest contribution to nearest
        E1f(if,jf) = E1f(if,jf) + w1*E1c(i,j)
        E1f(if+1,jf) = E1f(if+1,jf) + w1*E1c(i,j)
        E1f(if,jf+1) = E1f(if,jf+1) + w1*E1c(i,j)
        E1f(if+1,jf+1) = E1f(if+1,jf+1) + w1*E1c(i,j)

        ! lesser contribution to intermediate
        E1f(if-1,jf) = E1f(if-1,jf) + w2*E1c(i,j)
        E1f(if-1,jf+1) = E1f(if-1,jf+1) + w2*E1c(i,j)
        E1f(if,jf-1) = E1f(if,jf-1) + w2*E1c(i,j)
        E1f(if,jf+2) = E1f(if,jf+2) + w2*E1c(i,j)
        E1f(if+1,jf-1) = E1f(if+1,jf-1) + w2*E1c(i,j)
        E1f(if+1,jf+2) = E1f(if+1,jf+2) + w2*E1c(i,j)
        E1f(if+2,jf) = E1f(if+2,jf) + w2*E1c(i,j)
        E1f(if+2,jf+1) = E1f(if+2,jf+1) + w2*E1c(i,j)

        ! least contribution to furthest
        E1f(if-1,jf-1) = E1f(if-1,jf-1) + w3*E1c(i,j)
        E1f(if-1,jf+2) = E1f(if-1,jf+2) + w3*E1c(i,j)
        E1f(if+2,jf-1) = E1f(if+2,jf-1) + w3*E1c(i,j)
        E1f(if+2,jf+2) = E1f(if+2,jf+2) + w3*E1c(i,j)

        ! largest contribution to nearest
        E2f(if,jf) = E2f(if,jf) + w1*E2c(i,j)
        E2f(if+1,jf) = E2f(if+1,jf) + w1*E2c(i,j)
        E2f(if,jf+1) = E2f(if,jf+1) + w1*E2c(i,j)
        E2f(if+1,jf+1) = E2f(if+1,jf+1) + w1*E2c(i,j)

        ! lesser contribution to intermediate
        E2f(if-1,jf) = E2f(if-1,jf) + w2*E2c(i,j)
        E2f(if-1,jf+1) = E2f(if-1,jf+1) + w2*E2c(i,j)
        E2f(if,jf-1) = E2f(if,jf-1) + w2*E2c(i,j)
        E2f(if,jf+2) = E2f(if,jf+2) + w2*E2c(i,j)
        E2f(if+1,jf-1) = E2f(if+1,jf-1) + w2*E2c(i,j)
        E2f(if+1,jf+2) = E2f(if+1,jf+2) + w2*E2c(i,j)
        E2f(if+2,jf) = E2f(if+2,jf) + w2*E2c(i,j)
        E2f(if+2,jf+1) = E2f(if+2,jf+1) + w2*E2c(i,j)

        ! least contribution to furthest
        E2f(if-1,jf-1) = E2f(if-1,jf-1) + w3*E2c(i,j)
        E2f(if-1,jf+2) = E2f(if-1,jf+2) + w3*E2c(i,j)
        E2f(if+2,jf-1) = E2f(if+2,jf-1) + w3*E2c(i,j)
        E2f(if+2,jf+2) = E2f(if+2,jf+2) + w3*E2c(i,j)
      enddo
    enddo

    ! bottom
    do j=n,n
      jf = 2*j-1
      do i=1,n
        if = 2*i-1

        ! largest contribution to nearest
        E1f(if,jf) = E1f(if,jf) + w1*E1c(i,j)
        E1f(if+1,jf) = E1f(if+1,jf) + w1*E1c(i,j)
        E1f(if,jf+1) = E1f(if,jf+1) + w1*E1c(i,j)
        E1f(if+1,jf+1) = E1f(if+1,jf+1) + w1*E1c(i,j)

        ! lesser contribution to intermediate
        E1f(if-1,jf) = E1f(if-1,jf) + w2*E1c(i,j)
        E1f(if-1,jf+1) = E1f(if-1,jf+1) + w2*E1c(i,j)
        E1f(if,jf-1) = E1f(if,jf-1) + w2*E1c(i,j)
        E1f(if,jf+2) = E1f(if,jf+2) + w2*E1c(i,j)
        E1f(if+1,jf-1) = E1f(if+1,jf-1) + w2*E1c(i,j)
        E1f(if+1,jf+2) = E1f(if+1,jf+2) + w2*E1c(i,j)
        E1f(if+2,jf) = E1f(if+2,jf) + w2*E1c(i,j)
        E1f(if+2,jf+1) = E1f(if+2,jf+1) + w2*E1c(i,j)

        ! least contribution to furthest
        E1f(if-1,jf-1) = E1f(if-1,jf-1) + w3*E1c(i,j)
        E1f(if-1,jf+2) = E1f(if-1,jf+2) + w3*E1c(i,j)
        E1f(if+2,jf-1) = E1f(if+2,jf-1) + w3*E1c(i,j)
        E1f(if+2,jf+2) = E1f(if+2,jf+2) + w3*E1c(i,j)

        ! largest contribution to nearest
        E2f(if,jf) = E2f(if,jf) + w1*E2c(i,j)
        E2f(if+1,jf) = E2f(if+1,jf) + w1*E2c(i,j)
        E2f(if,jf+1) = E2f(if,jf+1) + w1*E2c(i,j)
        E2f(if+1,jf+1) = E2f(if+1,jf+1) + w1*E2c(i,j)

        ! lesser contribution to intermediate
        E2f(if-1,jf) = E2f(if-1,jf) + w2*E2c(i,j)
        E2f(if-1,jf+1) = E2f(if-1,jf+1) + w2*E2c(i,j)
        E2f(if,jf-1) = E2f(if,jf-1) + w2*E2c(i,j)
        E2f(if,jf+2) = E2f(if,jf+2) + w2*E2c(i,j)
        E2f(if+1,jf-1) = E2f(if+1,jf-1) + w2*E2c(i,j)
        E2f(if+1,jf+2) = E2f(if+1,jf+2) + w2*E2c(i,j)
        E2f(if+2,jf) = E2f(if+2,jf) + w2*E2c(i,j)
        E2f(if+2,jf+1) = E2f(if+2,jf+1) + w2*E2c(i,j)

        ! least contribution to furthest
        E2f(if-1,jf-1) = E2f(if-1,jf-1) + w3*E2c(i,j)
        E2f(if-1,jf+2) = E2f(if-1,jf+2) + w3*E2c(i,j)
        E2f(if+2,jf-1) = E2f(if+2,jf-1) + w3*E2c(i,j)
        E2f(if+2,jf+2) = E2f(if+2,jf+2) + w3*E2c(i,j)
      enddo
    enddo
  end subroutine prolongate


  ! ========================================================================== !
  !   SERIAL SUBROUTINES                                                       !
  ! ========================================================================== !
  subroutine laplacian0(x, y, dx, n)
    implicit none
    integer, intent(in) :: n
    real(dp), dimension(0:n+1,0:n+1), intent(in) :: x
    real(dp), dimension(0:n+1,0:n+1), intent(inout) :: y
    real(dp), intent(in) :: dx
    integer :: i, j
    real(dp) :: dx2_ ! interim constants

    dx2_ = 1.0_dp / (dx*dx)

    ! interior
    do j=2,n-1
      do i=2,n-1
        y(i,j) = dx2_*(x(i+1,j) + x(i-1,j) + x(i,j+1) + x(i,j-1) - 4*x(i,j))
      enddo
    enddo

    ! left/right
    do j=2,n-1
      y(1,j) = dx2_*(x(2,j) + x(n,j) + x(1,j+1) + x(1,j-1) - 4*x(1,j))
      y(n,j) = dx2_*(x(1,j) + x(n-1,j) + x(n,j+1) + x(n,j-1) - 4*x(n,j))
    enddo

    ! top/bottom
    do i=2,n-1
      y(i,1) = dx2_*(x(i+1,1) + x(i-1,1) + x(i,2) + x(i,1) - 4*x(i,1))
      y(i,n) = dx2_*(x(i+1,n) + x(i-1,n) + x(i,1) + x(i,n-1) - 4*x(i,n))
    enddo

    ! corners
    y(1,1) = dx2_*(x(2,1) + x(n,1) + x(1,2) + x(1,n) - 4*x(1,1))
    y(n,1) = dx2_*(x(1,1) + x(n-1,1) + x(n,2) + x(n,n) - 4*x(n,1))
    y(1,n) = dx2_*(x(2,n) + x(n,n) + x(1,1) + x(1,n-1) - 4*x(1,n))
    y(n,n) = dx2_*(x(1,n) + x(n-1,n) + x(n,1) + x(n,n-1) - 4*x(n,n))
  end subroutine laplacian0

  subroutine compute_g0(g, phi, dx, n, work)
    implicit none
    integer, intent(in) :: n
    real(dp), dimension(0:n+1,0:n+1), intent(in) :: phi
    real(dp), dimension(0:n+1,0:n+1), intent(inout) :: work
    real(dp), dimension(0:n+1,0:n+1), intent(inout) :: g
    real(dp), intent(in) :: dx

    work = phi * (phi*phi - (1+tau))

    call laplacian0(work, g, dx, n)
  end subroutine compute_g0

  subroutine vcycle0(A, E1, E2, R1, R2, eps2, N, dx, level)
    implicit none
    real(dp), dimension(2,2), intent(in) :: A
    type(t_grid), dimension(:), intent(inout) :: E1, E2, R1, R2
    real(dp), intent(in) :: eps2, dx
    integer, intent(in) :: N, level

    integer :: l, nl
    real(dp) :: dxl

    nl = n
    dxl = dx

    ! go up, smoothing and restricting
    do l=level,2,-1
      E1(l)%grid = 0.0_dp
      E2(l)%grid = 0.0_dp

      call smooth0(A, E1(l)%grid, E2(l)%grid, R1(l)%grid, R2(l)%grid, eps2, nl, dxl)
      call restrict0(R1(l)%grid, R2(l)%grid, R1(l-1)%grid, R2(l-1)%grid, nl)

      nl = nl/2;
      dxl = dxl * 2.0_dp
    enddo

    ! smooth at level 1
    E1(l)%grid = 0.0_dp
    E2(l)%grid = 0.0_dp
    call smooth0(A, E1(1)%grid, E2(1)%grid, R1(1)%grid, R2(1)%grid, eps2, nl, dxl)
    call smooth0(A, E1(1)%grid, E2(1)%grid, R1(1)%grid, R2(1)%grid, eps2, nl, dxl)
    call smooth0(A, E1(1)%grid, E2(1)%grid, R1(1)%grid, R2(1)%grid, eps2, nl, dxl)
    call smooth0(A, E1(1)%grid, E2(1)%grid, R1(1)%grid, R2(1)%grid, eps2, nl, dxl)

    ! go down, smoothing and prolongating
    do l=1,level-1
      call prolongate0(E1(l+1)%grid, E2(l+1)%grid, E1(l)%grid, E2(l)%grid, nl)

      nl = nl*2;
      dxl = dxl * 0.5_dp

      call smooth0(A, E1(l+1)%grid, E2(l+1)%grid, R1(l+1)%grid, R2(l+1)%grid, eps2, nl, dxl)
    enddo
  end subroutine vcycle0

  subroutine smooth0(A, E1, E2, R1, R2, eps2, N, dx)
    implicit none
    integer, intent(in) :: n
    real(dp), dimension(2,2), intent(in) :: A
    real(dp), dimension(0:n+1,0:n+1), intent(inout) :: E1, E2, R1, R2
    real(dp), intent(in) :: eps2, dx

    real(dp), dimension(2) :: rhs
    real(dp) :: dx2_
    integer :: i, j, shift

    dx2_ = 1.0_dp / (dx*dx)

    ! ================ !
    ! SMOOTH RED NODES !
    ! ================ !
    ! interior
    do j=2,n-1
      shift = mod(j,2)
      do i=2+shift,n-1+shift,2
        ! compute RHS
        rhs(1) = R1(i,j) + dx2_ * (E2(i+1,j) + E2(i-1,j) + E2(i,j+1) + E2(i,j-1))
        rhs(2) = R2(i,j) - eps2*dx2_ * (E1(i+1,j) + E1(i-1,j) + E1(i,j+1) + E1(i,j-1))

        ! solve for new errors
        E1(i,j) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
        E2(i,j) = A(2,1)*rhs(1) + A(2,2)*rhs(2)
      enddo
    enddo

    ! left/right
    do j=2,n-1,2
      rhs(1) = R1(1,j+1) + dx2_ * (E2(2,j+1) + E2(n,j+1) + E2(1,j+2) + E2(1,j))
      rhs(2) = R2(1,j+1) - eps2*dx2_ * (E1(2,j+1) + E1(n,j+1) + E1(1,j+2) + E1(1,j))
      E1(1,j+1) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
      E2(1,j+1) = A(2,1)*rhs(1) + A(2,2)*rhs(2)

      rhs(1) = R1(n,j) + dx2_ * (E2(1,j) + E2(n-1,j) + E2(n,j+1) + E2(n,j-1))
      rhs(2) = R2(n,j) - eps2*dx2_ * (E1(1,j) + E1(n-1,j) + E1(n,j+1) + E1(n,j-1))
      E1(n,j) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
      E2(n,j) = A(2,1)*rhs(1) + A(2,2)*rhs(2)
    enddo

    ! top/bottom
    do i=2,n-1,2
      rhs(1) = R1(i+1,1) + dx2_ * (E2(i+2,1) + E2(i,1) + E2(i+1,2) + E2(i+1,n))
      rhs(2) = R2(i+1,1) - eps2*dx2_ * (E1(i+2,1) + E1(i,1) + E1(i+1,2) + E1(i+1,n))
      E1(i+1,1) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
      E2(i+1,1) = A(2,1)*rhs(1) + A(2,2)*rhs(2)

      rhs(1) = R1(i,n) + dx2_ * (E2(i+1,n) + E2(i-1,n) + E2(i,1) + E2(i,n-1))
      rhs(2) = R2(i,n) - eps2*dx2_ * (E1(i+1,n) + E1(i-1,n) + E1(i,1) + E1(i,n-1))
      E1(i,n) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
      E2(i,n) = A(2,1)*rhs(1) + A(2,2)*rhs(2)
    enddo

    ! upper left corner
    rhs(1) = R1(1,1) + dx2_ * (E2(2,1) + E2(n,1) + E2(1,2) + E2(1,n))
    rhs(2) = R2(1,1) - eps2*dx2_ * (E1(2,1) + E1(n,1) + E1(1,2) + E1(1,n))
    E1(1,1) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
    E2(1,1) = A(2,1)*rhs(1) + A(2,2)*rhs(2)

    ! lower right corner
    rhs(1) = R1(n,n) + dx2_ * (E2(1,n) + E2(n-1,n) + E2(n,1) + E2(n,n-1))
    rhs(2) = R2(n,n) - eps2*dx2_ * (E1(1,n) + E1(n-1,n) + E1(n,1) + E1(n,n-1))
    E1(n,n) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
    E2(n,n) = A(2,1)*rhs(1) + A(2,2)*rhs(2)


    ! ================== !
    ! SMOOTH BLACK NODES !
    ! ================== !
    ! interior
    do j=2,n-1
      shift = mod(j-1,2)
      do i=2+shift,n-1+shift,2
        ! compute RHS
        rhs(1) = R1(i,j) + dx2_ * (E2(i+1,j) + E2(i-1,j) + E2(i,j+1) + E2(i,j-1))
        rhs(2) = R2(i,j) - eps2*dx2_ * (E1(i+1,j) + E1(i-1,j) + E1(i,j+1) + E1(i,j-1))

        ! solve for new errors
        E1(i,j) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
        E2(i,j) = A(2,1)*rhs(1) + A(2,2)*rhs(2)
      enddo
    enddo

    ! left/right
    do j=2,n-1,2
      rhs(1) = R1(1,j) + dx2_ * (E2(2,j) + E2(n,j) + E2(1,j+1) + E2(1,j-1))
      rhs(2) = R2(1,j) - eps2*dx2_ * (E1(2,j) + E1(n,j) + E1(1,j+1) + E1(1,j-1))
      E1(1,j) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
      E2(1,j) = A(2,1)*rhs(1) + A(2,2)*rhs(2)

      rhs(1) = R1(n,j+1) + dx2_ * (E2(1,j+1) + E2(n-1,j+1) + E2(n,j+2) + E2(n,j))
      rhs(2) = R2(n,j+1) - eps2*dx2_ * (E1(1,j+1) + E1(n-1,j+1) + E1(n,j+2) + E1(n,j))
      E1(n,j+1) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
      E2(n,j+1) = A(2,1)*rhs(1) + A(2,2)*rhs(2)
    enddo

    ! top/bottom
    do i=2,n-1,2
      rhs(1) = R1(i,1) + dx2_ * (E2(i+1,1) + E2(i-1,1) + E2(i,2) + E2(i,n))
      rhs(2) = R2(i,1) - eps2*dx2_ * (E1(i+1,1) + E1(i-1,1) + E1(i,2) + E1(i,n))
      E1(i,1) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
      E2(i,1) = A(2,1)*rhs(1) + A(2,2)*rhs(2)

      rhs(1) = R1(i+1,n) + dx2_ * (E2(i+2,n) + E2(i,n) + E2(i+1,1) + E2(i+1,n-1))
      rhs(2) = R2(i+1,n) - eps2*dx2_ * (E1(i+2,n) + E1(i,n) + E1(i+1,1) + E1(i+1,n-1))
      E1(i+1,n) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
      E2(i+1,n) = A(2,1)*rhs(1) + A(2,2)*rhs(2)
    enddo

    ! upper right
    rhs(1) = R1(1,n) + dx2_ * (E2(2,n) + E2(n,n) + E2(1,1) + E2(1,n-1))
    rhs(2) = R2(1,n) - eps2*dx2_ * (E1(2,n) + E1(n,n) + E1(1,1) + E1(1,n-1))
    E1(1,n) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
    E2(1,n) = A(2,1)*rhs(1) + A(2,2)*rhs(2)

    ! lower left
    rhs(1) = R1(n,1) + dx2_ * (E2(1,1) + E2(n-1,1) + E2(n,2) + E2(n,n))
    rhs(2) = R2(n,1) - eps2*dx2_ * (E1(1,1) + E1(n-1,1) + E1(n,2) + E1(n,n))
    E1(n,1) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
    E2(n,1) = A(2,1)*rhs(1) + A(2,2)*rhs(2)
  end subroutine smooth0

  subroutine restrict0(R1f, R2f, R1c, R2c, N)
    implicit none
    integer, intent(in) :: n
    real(dp), dimension(0:n+1,0:n+1), intent(inout) :: R1f, R2f
    real(dp), dimension(0:n/2+1,0:n/2+1), intent(inout) :: R1c, R2c

    integer :: i, j, im, jm, nc
    nc = n/2

    ! use a simple average across the grid
    do j=1,nc
      jm = 2*j-1
      do i=1,nc
        im = 2*i-1

        ! fine -> coarse
        R1c(i,j) = 0.25_dp*(R1f(im,jm) + R1f(im+1,jm) + R1f(im,jm+1) + R1f(im+1,jm+1))
        R2c(i,j) = 0.25_dp*(R2f(im,jm) + R2f(im+1,jm) + R2f(im,jm+1) + R2f(im+1,jm+1))
      enddo
    enddo
  end subroutine restrict0

  subroutine prolongate0(E1f, E2f, E1c, E2c, N)
    implicit none
    integer, intent(in) :: n
    real(dp), dimension(0:n+1,0:n+1), intent(inout) :: E1c, E2c
    real(dp), dimension(0:2*n+1,0:2*n+1), intent(inout) :: E1f, E2f

    integer :: Nf
    integer :: i, j, if, jf
    real(dp), parameter :: w1 = 0.5625_dp
    real(dp), parameter :: w2 = 0.1875_dp
    real(dp), parameter :: w3 = 0.0625_dp

    Nf = 2*N

    ! interior
    do j=2,N-1
      jf = 2*j-1
      do i=2,N-1
        if = 2*i-1

        ! largest contribution to nearest
        E1f(if,jf) = E1f(if,jf) + w1*E1c(i,j)
        E1f(if+1,jf) = E1f(if+1,jf) + w1*E1c(i,j)
        E1f(if,jf+1) = E1f(if,jf+1) + w1*E1c(i,j)
        E1f(if+1,jf+1) = E1f(if+1,jf+1) + w1*E1c(i,j)

        ! lesser contribution to intermediate
        E1f(if-1,jf) = E1f(if-1,jf) + w2*E1c(i,j)
        E1f(if-1,jf+1) = E1f(if-1,jf+1) + w2*E1c(i,j)
        E1f(if,jf-1) = E1f(if,jf-1) + w2*E1c(i,j)
        E1f(if,jf+2) = E1f(if,jf+2) + w2*E1c(i,j)
        E1f(if+1,jf-1) = E1f(if+1,jf-1) + w2*E1c(i,j)
        E1f(if+1,jf+2) = E1f(if+1,jf+2) + w2*E1c(i,j)
        E1f(if+2,jf) = E1f(if+2,jf) + w2*E1c(i,j)
        E1f(if+2,jf+1) = E1f(if+2,jf+1) + w2*E1c(i,j)

        ! least contribution to furthest
        E1f(if-1,jf-1) = E1f(if-1,jf-1) + w3*E1c(i,j)
        E1f(if-1,jf+2) = E1f(if-1,jf+2) + w3*E1c(i,j)
        E1f(if+2,jf-1) = E1f(if+2,jf-1) + w3*E1c(i,j)
        E1f(if+2,jf+2) = E1f(if+2,jf+2) + w3*E1c(i,j)

        ! largest contribution to nearest
        E2f(if,jf) = E2f(if,jf) + w1*E2c(i,j)
        E2f(if+1,jf) = E2f(if+1,jf) + w1*E2c(i,j)
        E2f(if,jf+1) = E2f(if,jf+1) + w1*E2c(i,j)
        E2f(if+1,jf+1) = E2f(if+1,jf+1) + w1*E2c(i,j)

        ! lesser contribution to intermediate
        E2f(if-1,jf) = E2f(if-1,jf) + w2*E2c(i,j)
        E2f(if-1,jf+1) = E2f(if-1,jf+1) + w2*E2c(i,j)
        E2f(if,jf-1) = E2f(if,jf-1) + w2*E2c(i,j)
        E2f(if,jf+2) = E2f(if,jf+2) + w2*E2c(i,j)
        E2f(if+1,jf-1) = E2f(if+1,jf-1) + w2*E2c(i,j)
        E2f(if+1,jf+2) = E2f(if+1,jf+2) + w2*E2c(i,j)
        E2f(if+2,jf) = E2f(if+2,jf) + w2*E2c(i,j)
        E2f(if+2,jf+1) = E2f(if+2,jf+1) + w2*E2c(i,j)

        ! least contribution to furthest
        E2f(if-1,jf-1) = E2f(if-1,jf-1) + w3*E2c(i,j)
        E2f(if-1,jf+2) = E2f(if-1,jf+2) + w3*E2c(i,j)
        E2f(if+2,jf-1) = E2f(if+2,jf-1) + w3*E2c(i,j)
        E2f(if+2,jf+2) = E2f(if+2,jf+2) + w3*E2c(i,j)
      enddo
    enddo

    ! left edge
    do j=2,N-1
      jf = 2*j-1

      ! largest contribution to nearest
      E1f(1,jf) = E1f(1,jf) + w1*E1c(1,j)
      E1f(2,jf) = E1f(2,jf) + w1*E1c(1,j)
      E1f(1,jf+1) = E1f(1,jf+1) + w1*E1c(1,j)
      E1f(2,jf+1) = E1f(2,jf+1) + w1*E1c(1,j)

      ! lesser contribution to intermediate
      E1f(Nf,jf) = E1f(Nf,jf) + w2*E1c(1,j)
      E1f(Nf,jf+1) = E1f(Nf,jf+1) + w2*E1c(1,j)
      E1f(1,jf-1) = E1f(1,jf-1) + w2*E1c(1,j)
      E1f(1,jf+2) = E1f(1,jf+2) + w2*E1c(1,j)
      E1f(2,jf-1) = E1f(2,jf-1) + w2*E1c(1,j)
      E1f(2,jf+2) = E1f(2,jf+2) + w2*E1c(1,j)
      E1f(3,jf) = E1f(3,jf) + w2*E1c(1,j)
      E1f(3,jf+1) = E1f(3,jf+1) + w2*E1c(1,j)

      ! least contribution to furthest
      E1f(Nf,jf-1) = E1f(Nf,jf-1) + w3*E1c(1,j)
      E1f(Nf,jf+2) = E1f(Nf,jf+2) + w3*E1c(1,j)
      E1f(3,jf-1) = E1f(3,jf-1) + w3*E1c(1,j)
      E1f(3,jf+2) = E1f(3,jf+2) + w3*E1c(1,j)

      ! largest contribution to nearest
      E2f(1,jf) = E2f(1,jf) + w1*E2c(1,j)
      E2f(2,jf) = E2f(2,jf) + w1*E2c(1,j)
      E2f(1,jf+1) = E2f(1,jf+1) + w1*E2c(1,j)
      E2f(2,jf+1) = E2f(2,jf+1) + w1*E2c(1,j)

      ! lesser contribution to intermediate
      E2f(Nf,jf) = E2f(Nf,jf) + w2*E2c(1,j)
      E2f(Nf,jf+1) = E2f(Nf,jf+1) + w2*E2c(1,j)
      E2f(1,jf-1) = E2f(1,jf-1) + w2*E2c(1,j)
      E2f(1,jf+2) = E2f(1,jf+2) + w2*E2c(1,j)
      E2f(2,jf-1) = E2f(2,jf-1) + w2*E2c(1,j)
      E2f(2,jf+2) = E2f(2,jf+2) + w2*E2c(1,j)
      E2f(3,jf) = E2f(3,jf) + w2*E2c(1,j)
      E2f(3,jf+1) = E2f(3,jf+1) + w2*E2c(1,j)

      ! least contribution to furthest
      E2f(Nf,jf-1) = E2f(Nf,jf-1) + w3*E2c(1,j)
      E2f(Nf,jf+2) = E2f(Nf,jf+2) + w3*E2c(1,j)
      E2f(3,jf-1) = E2f(3,jf-1) + w3*E2c(1,j)
      E2f(3,jf+2) = E2f(3,jf+2) + w3*E2c(1,j)
    enddo

    ! right edge
    do j=2,N-1
      jf = 2*j-1

      ! largest contribution to nearest
      E1f(Nf-1,jf) = E1f(Nf-1,jf) + w1*E1c(N,j)
      E1f(Nf,jf) = E1f(Nf,jf) + w1*E1c(N,j)
      E1f(Nf-1,jf+1) = E1f(Nf-1,jf+1) + w1*E1c(N,j)
      E1f(Nf,jf+1) = E1f(Nf,jf+1) + w1*E1c(N,j)

      ! lesser contribution to intermediate
      E1f(Nf-2,jf) = E1f(Nf-2,jf) + w2*E1c(N,j)
      E1f(Nf-2,jf+1) = E1f(Nf-2,jf+1) + w2*E1c(N,j)
      E1f(Nf-1,jf-1) = E1f(Nf-1,jf-1) + w2*E1c(N,j)
      E1f(Nf-1,jf+2) = E1f(Nf-1,jf+2) + w2*E1c(N,j)
      E1f(Nf,jf-1) = E1f(Nf,jf-1) + w2*E1c(N,j)
      E1f(Nf,jf+2) = E1f(Nf,jf+2) + w2*E1c(N,j)
      E1f(1,jf) = E1f(1,jf) + w2*E1c(N,j)
      E1f(1,jf+1) = E1f(1,jf+1) + w2*E1c(N,j)

      ! least contribution to furthest
      E1f(Nf-2,jf-1) = E1f(Nf-2,jf-1) + w3*E1c(N,j)
      E1f(Nf-2,jf+2) = E1f(Nf-2,jf+2) + w3*E1c(N,j)
      E1f(1,jf-1) = E1f(1,jf-1) + w3*E1c(N,j)
      E1f(1,jf+2) = E1f(1,jf+2) + w3*E1c(N,j)

      ! largest contribution to nearest
      E2f(Nf-1,jf) = E2f(Nf-1,jf) + w1*E2c(N,j)
      E2f(Nf,jf) = E2f(Nf,jf) + w1*E2c(N,j)
      E2f(Nf-1,jf+1) = E2f(Nf-1,jf+1) + w1*E2c(N,j)
      E2f(Nf,jf+1) = E2f(Nf,jf+1) + w1*E2c(N,j)

      ! lesser contribution to intermediate
      E2f(Nf-2,jf) = E2f(Nf-2,jf) + w2*E2c(N,j)
      E2f(Nf-2,jf+1) = E2f(Nf-2,jf+1) + w2*E2c(N,j)
      E2f(Nf-1,jf-1) = E2f(Nf-1,jf-1) + w2*E2c(N,j)
      E2f(Nf-1,jf+2) = E2f(Nf-1,jf+2) + w2*E2c(N,j)
      E2f(Nf,jf-1) = E2f(Nf,jf-1) + w2*E2c(N,j)
      E2f(Nf,jf+2) = E2f(Nf,jf+2) + w2*E2c(N,j)
      E2f(1,jf) = E2f(1,jf) + w2*E2c(N,j)
      E2f(1,jf+1) = E2f(1,jf+1) + w2*E2c(N,j)

      ! least contribution to furthest
      E2f(Nf-2,jf-1) = E2f(Nf-2,jf-1) + w3*E2c(N,j)
      E2f(Nf-2,jf+2) = E2f(Nf-2,jf+2) + w3*E2c(N,j)
      E2f(1,jf-1) = E2f(1,jf-1) + w3*E2c(N,j)
      E2f(1,jf+2) = E2f(1,jf+2) + w3*E2c(N,j)
    enddo

    ! top edge
    do i=2,N-1
      if = 2*i-1

      ! largest contribution to nearest
      E1f(if,1) = E1f(if,1) + w1*E1c(i,1)
      E1f(if+1,1) = E1f(if+1,1) + w1*E1c(i,1)
      E1f(if,2) = E1f(if,2) + w1*E1c(i,1)
      E1f(if+1,2) = E1f(if+1,2) + w1*E1c(i,1)

      ! lesser contribution to intermediate
      E1f(if-1,1) = E1f(if-1,1) + w2*E1c(i,1)
      E1f(if-1,2) = E1f(if-1,2) + w2*E1c(i,1)
      E1f(if,Nf) = E1f(if,Nf) + w2*E1c(i,1)
      E1f(if,3) = E1f(if,3) + w2*E1c(i,1)
      E1f(if+1,Nf) = E1f(if+1,Nf) + w2*E1c(i,1)
      E1f(if+1,3) = E1f(if+1,3) + w2*E1c(i,1)
      E1f(if+2,1) = E1f(if+2,1) + w2*E1c(i,1)
      E1f(if+2,2) = E1f(if+2,2) + w2*E1c(i,1)

      ! least contribution to furthest
      E1f(if-1,Nf) = E1f(if-1,Nf) + w3*E1c(i,1)
      E1f(if-1,3) = E1f(if-1,3) + w3*E1c(i,1)
      E1f(if+2,Nf) = E1f(if+2,Nf) + w3*E1c(i,1)
      E1f(if+2,3) = E1f(if+2,3) + w3*E1c(i,1)

      ! largest contribution to nearest
      E2f(if,1) = E2f(if,1) + w1*E2c(i,1)
      E2f(if+1,1) = E2f(if+1,1) + w1*E2c(i,1)
      E2f(if,2) = E2f(if,2) + w1*E2c(i,1)
      E2f(if+1,2) = E2f(if+1,2) + w1*E2c(i,1)

      ! lesser contribution to intermediate
      E2f(if-1,1) = E2f(if-1,1) + w2*E2c(i,1)
      E2f(if-1,2) = E2f(if-1,2) + w2*E2c(i,1)
      E2f(if,Nf) = E2f(if,Nf) + w2*E2c(i,1)
      E2f(if,3) = E2f(if,3) + w2*E2c(i,1)
      E2f(if+1,Nf) = E2f(if+1,Nf) + w2*E2c(i,1)
      E2f(if+1,3) = E2f(if+1,3) + w2*E2c(i,1)
      E2f(if+2,1) = E2f(if+2,1) + w2*E2c(i,1)
      E2f(if+2,2) = E2f(if+2,2) + w2*E2c(i,1)

      ! least contribution to furthest
      E2f(if-1,Nf) = E2f(if-1,Nf) + w3*E2c(i,1)
      E2f(if-1,3) = E2f(if-1,3) + w3*E2c(i,1)
      E2f(if+2,Nf) = E2f(if+2,Nf) + w3*E2c(i,1)
      E2f(if+2,3) = E2f(if+2,3) + w3*E2c(i,1)
    enddo

    ! bottom edge
    do i=2,N-1
      if = 2*i-1

      ! largest contribution to nearest
      E1f(if,Nf-1) = E1f(if,Nf-1) + w1*E1c(i,N)
      E1f(if+1,Nf-1) = E1f(if+1,Nf-1) + w1*E1c(i,N)
      E1f(if,Nf) = E1f(if,Nf) + w1*E1c(i,N)
      E1f(if+1,Nf) = E1f(if+1,Nf) + w1*E1c(i,N)

      ! lesser contribution to intermediate
      E1f(if-1,Nf-1) = E1f(if-1,Nf-1) + w2*E1c(i,N)
      E1f(if-1,Nf) = E1f(if-1,Nf) + w2*E1c(i,N)
      E1f(if,Nf-2) = E1f(if,Nf-2) + w2*E1c(i,N)
      E1f(if,1) = E1f(if,1) + w2*E1c(i,N)
      E1f(if+1,Nf-2) = E1f(if+1,Nf-2) + w2*E1c(i,N)
      E1f(if+1,1) = E1f(if+1,1) + w2*E1c(i,N)
      E1f(if+2,Nf-1) = E1f(if+2,Nf-1) + w2*E1c(i,N)
      E1f(if+2,Nf) = E1f(if+2,Nf) + w2*E1c(i,N)

      ! least contribution to furthest
      E1f(if-1,Nf-2) = E1f(if-1,Nf-2) + w3*E1c(i,N)
      E1f(if-1,1) = E1f(if-1,1) + w3*E1c(i,N)
      E1f(if+2,Nf-2) = E1f(if+2,Nf-2) + w3*E1c(i,N)
      E1f(if+2,1) = E1f(if+2,1) + w3*E1c(i,N)

      ! largest contribution to nearest
      E2f(if,Nf-1) = E2f(if,Nf-1) + w1*E2c(i,N)
      E2f(if+1,Nf-1) = E2f(if+1,Nf-1) + w1*E2c(i,N)
      E2f(if,Nf) = E2f(if,Nf) + w1*E2c(i,N)
      E2f(if+1,Nf) = E2f(if+1,Nf) + w1*E2c(i,N)

      ! lesser contribution to intermediate
      E2f(if-1,Nf-1) = E2f(if-1,Nf-1) + w2*E2c(i,N)
      E2f(if-1,Nf) = E2f(if-1,Nf) + w2*E2c(i,N)
      E2f(if,Nf-2) = E2f(if,Nf-2) + w2*E2c(i,N)
      E2f(if,1) = E2f(if,1) + w2*E2c(i,N)
      E2f(if+1,Nf-2) = E2f(if+1,Nf-2) + w2*E2c(i,N)
      E2f(if+1,1) = E2f(if+1,1) + w2*E2c(i,N)
      E2f(if+2,Nf-1) = E2f(if+2,Nf-1) + w2*E2c(i,N)
      E2f(if+2,Nf) = E2f(if+2,Nf) + w2*E2c(i,N)

      ! least contribution to furthest
      E2f(if-1,Nf-2) = E2f(if-1,Nf-2) + w3*E2c(i,N)
      E2f(if-1,1) = E2f(if-1,1) + w3*E2c(i,N)
      E2f(if+2,Nf-2) = E2f(if+2,Nf-2) + w3*E2c(i,N)
      E2f(if+2,1) = E2f(if+2,1) + w3*E2c(i,N)
    enddo

    ! top left corner
    ! largest contribution to nearest
    E1f(1,1) = E1f(1,1) + w1*E1c(1,1)
    E1f(2,1) = E1f(2,1) + w1*E1c(1,1)
    E1f(1,2) = E1f(1,2) + w1*E1c(1,1)
    E1f(2,2) = E1f(2,2) + w1*E1c(1,1)

    ! lesser contribution to intermediate
    E1f(Nf,1) = E1f(Nf,1) + w2*E1c(1,1)
    E1f(Nf,2) = E1f(Nf,2) + w2*E1c(1,1)
    E1f(1,Nf) = E1f(1,Nf) + w2*E1c(1,1)
    E1f(1,3) = E1f(1,3) + w2*E1c(1,1)
    E1f(2,Nf) = E1f(2,Nf) + w2*E1c(1,1)
    E1f(2,3) = E1f(2,3) + w2*E1c(1,1)
    E1f(3,1) = E1f(3,1) + w2*E1c(1,1)
    E1f(3,2) = E1f(3,2) + w2*E1c(1,1)

    ! least contribution to furthest
    E1f(Nf,Nf) = E1f(Nf,Nf) + w3*E1c(1,1)
    E1f(Nf,3) = E1f(Nf,3) + w3*E1c(1,1)
    E1f(3,Nf) = E1f(3,Nf) + w3*E1c(1,1)
    E1f(3,3) = E1f(3,3) + w3*E1c(1,1)

    ! largest contribution to nearest
    E2f(1,1) = E2f(1,1) + w1*E2c(1,1)
    E2f(2,1) = E2f(2,1) + w1*E2c(1,1)
    E2f(1,2) = E2f(1,2) + w1*E2c(1,1)
    E2f(2,2) = E2f(2,2) + w1*E2c(1,1)

    ! lesser contribution to intermediate
    E2f(Nf,1) = E2f(Nf,1) + w2*E2c(1,1)
    E2f(Nf,2) = E2f(Nf,2) + w2*E2c(1,1)
    E2f(1,Nf) = E2f(1,Nf) + w2*E2c(1,1)
    E2f(1,3) = E2f(1,3) + w2*E2c(1,1)
    E2f(2,Nf) = E2f(2,Nf) + w2*E2c(1,1)
    E2f(2,3) = E2f(2,3) + w2*E2c(1,1)
    E2f(3,1) = E2f(3,1) + w2*E2c(1,1)
    E2f(3,2) = E2f(3,2) + w2*E2c(1,1)

    ! least contribution to furthest
    E2f(Nf,Nf) = E2f(Nf,Nf) + w3*E2c(1,1)
    E2f(Nf,3) = E2f(Nf,3) + w3*E2c(1,1)
    E2f(3,Nf) = E2f(3,Nf) + w3*E2c(1,1)
    E2f(3,3) = E2f(3,3) + w3*E2c(1,1)


    ! top right corner
    ! largest contribution to nearest
    E1f(Nf-1,1) = E1f(Nf-1,1) + w1*E1c(N,1)
    E1f(Nf,1) = E1f(Nf,1) + w1*E1c(N,1)
    E1f(Nf-1,2) = E1f(Nf-1,2) + w1*E1c(N,1)
    E1f(Nf,2) = E1f(Nf,2) + w1*E1c(N,1)

    ! lesser contribution to intermediate
    E1f(Nf-2,1) = E1f(Nf-2,1) + w2*E1c(N,1)
    E1f(Nf-2,2) = E1f(Nf-2,2) + w2*E1c(N,1)
    E1f(Nf-1,Nf) = E1f(Nf-1,Nf) + w2*E1c(N,1)
    E1f(Nf-1,3) = E1f(Nf-1,3) + w2*E1c(N,1)
    E1f(Nf,Nf) = E1f(Nf,Nf) + w2*E1c(N,1)
    E1f(Nf,3) = E1f(Nf,3) + w2*E1c(N,1)
    E1f(1,1) = E1f(1,1) + w2*E1c(N,1)
    E1f(1,2) = E1f(1,2) + w2*E1c(N,1)

    ! least contribution to furthest
    E1f(Nf-2,Nf) = E1f(Nf-2,Nf) + w3*E1c(N,1)
    E1f(Nf-2,3) = E1f(Nf-2,3) + w3*E1c(N,1)
    E1f(1,Nf) = E1f(1,Nf) + w3*E1c(N,1)
    E1f(1,3) = E1f(1,3) + w3*E1c(N,1)

    ! largest contribution to nearest
    E2f(Nf-1,1) = E2f(Nf-1,1) + w1*E2c(N,1)
    E2f(Nf,1) = E2f(Nf,1) + w1*E2c(N,1)
    E2f(Nf-1,2) = E2f(Nf-1,2) + w1*E2c(N,1)
    E2f(Nf,2) = E2f(Nf,2) + w1*E2c(N,1)

    ! lesser contribution to intermediate
    E2f(Nf-2,1) = E2f(Nf-2,1) + w2*E2c(N,1)
    E2f(Nf-2,2) = E2f(Nf-2,2) + w2*E2c(N,1)
    E2f(Nf-1,Nf) = E2f(Nf-1,Nf) + w2*E2c(N,1)
    E2f(Nf-1,3) = E2f(Nf-1,3) + w2*E2c(N,1)
    E2f(Nf,Nf) = E2f(Nf,Nf) + w2*E2c(N,1)
    E2f(Nf,3) = E2f(Nf,3) + w2*E2c(N,1)
    E2f(1,1) = E2f(1,1) + w2*E2c(N,1)
    E2f(1,2) = E2f(1,2) + w2*E2c(N,1)

    ! least contribution to furthest
    E2f(Nf-2,Nf) = E2f(Nf-2,Nf) + w3*E2c(N,1)
    E2f(Nf-2,3) = E2f(Nf-2,3) + w3*E2c(N,1)
    E2f(1,Nf) = E2f(1,Nf) + w3*E2c(N,1)
    E2f(1,3) = E2f(1,3) + w3*E2c(N,1)


    ! bottom left corner
    ! largest contribution to nearest
    E1f(1,Nf-1) = E1f(1,Nf-1) + w1*E1c(1,N)
    E1f(2,Nf-1) = E1f(2,Nf-1) + w1*E1c(1,N)
    E1f(1,Nf) = E1f(1,Nf) + w1*E1c(1,N)
    E1f(2,Nf) = E1f(2,Nf) + w1*E1c(1,N)

    ! lesser contribution to intermediate
    E1f(Nf,Nf-1) = E1f(Nf,Nf-1) + w2*E1c(1,N)
    E1f(Nf,Nf) = E1f(Nf,Nf) + w2*E1c(1,N)
    E1f(1,Nf-2) = E1f(1,Nf-2) + w2*E1c(1,N)
    E1f(1,1) = E1f(1,1) + w2*E1c(1,N)
    E1f(2,Nf-2) = E1f(2,Nf-2) + w2*E1c(1,N)
    E1f(2,1) = E1f(2,1) + w2*E1c(1,N)
    E1f(3,Nf-1) = E1f(3,Nf-1) + w2*E1c(1,N)
    E1f(3,Nf) = E1f(3,Nf) + w2*E1c(1,N)

    ! least contribution to furthest
    E1f(Nf,Nf-2) = E1f(Nf,Nf-2) + w3*E1c(1,N)
    E1f(Nf,1) = E1f(Nf,1) + w3*E1c(1,N)
    E1f(3,Nf-2) = E1f(3,Nf-2) + w3*E1c(1,N)
    E1f(3,1) = E1f(3,1) + w3*E1c(1,N)

    ! largest contribution to nearest
    E2f(1,Nf-1) = E2f(1,Nf-1) + w1*E2c(1,N)
    E2f(2,Nf-1) = E2f(2,Nf-1) + w1*E2c(1,N)
    E2f(1,Nf) = E2f(1,Nf) + w1*E2c(1,N)
    E2f(2,Nf) = E2f(2,Nf) + w1*E2c(1,N)

    ! lesser contribution to intermediate
    E2f(Nf,Nf-1) = E2f(Nf,Nf-1) + w2*E2c(1,N)
    E2f(Nf,Nf) = E2f(Nf,Nf) + w2*E2c(1,N)
    E2f(1,Nf-2) = E2f(1,Nf-2) + w2*E2c(1,N)
    E2f(1,1) = E2f(1,1) + w2*E2c(1,N)
    E2f(2,Nf-2) = E2f(2,Nf-2) + w2*E2c(1,N)
    E2f(2,1) = E2f(2,1) + w2*E2c(1,N)
    E2f(3,Nf-1) = E2f(3,Nf-1) + w2*E2c(1,N)
    E2f(3,Nf) = E2f(3,Nf) + w2*E2c(1,N)

    ! least contribution to furthest
    E2f(Nf,Nf-2) = E2f(Nf,Nf-2) + w3*E2c(1,N)
    E2f(Nf,1) = E2f(Nf,1) + w3*E2c(1,N)
    E2f(3,Nf-2) = E2f(3,Nf-2) + w3*E2c(1,N)
    E2f(3,1) = E2f(3,1) + w3*E2c(1,N)


    ! bottom right corner
    ! largest contribution to nearest
    E1f(Nf-1,Nf-1) = E1f(Nf-1,Nf-1) + w1*E1c(N,N)
    E1f(Nf,Nf-1) = E1f(Nf,Nf-1) + w1*E1c(N,N)
    E1f(Nf-1,Nf) = E1f(Nf-1,Nf) + w1*E1c(N,N)
    E1f(Nf,Nf) = E1f(Nf,Nf) + w1*E1c(N,N)

    ! lesser contribution to intermediate
    E1f(Nf-2,Nf-1) = E1f(Nf-2,Nf-1) + w2*E1c(N,N)
    E1f(Nf-2,Nf) = E1f(Nf-2,Nf) + w2*E1c(N,N)
    E1f(Nf-1,Nf-2) = E1f(Nf-1,Nf-2) + w2*E1c(N,N)
    E1f(Nf-1,1) = E1f(Nf-1,1) + w2*E1c(N,N)
    E1f(Nf,Nf-2) = E1f(Nf,Nf-2) + w2*E1c(N,N)
    E1f(Nf,1) = E1f(Nf,1) + w2*E1c(N,N)
    E1f(1,Nf-1) = E1f(1,Nf-1) + w2*E1c(N,N)
    E1f(1,Nf) = E1f(1,Nf) + w2*E1c(N,N)

    ! least contribution to furthest
    E1f(Nf-2,Nf-2) = E1f(Nf-2,Nf-2) + w3*E1c(N,N)
    E1f(Nf-2,1) = E1f(Nf-2,1) + w3*E1c(N,N)
    E1f(1,Nf-2) = E1f(1,Nf-2) + w3*E1c(N,N)
    E1f(1,1) = E1f(1,1) + w3*E1c(N,N)

    ! largest contribution to nearest
    E2f(Nf-1,Nf-1) = E2f(Nf-1,Nf-1) + w1*E2c(N,N)
    E2f(Nf,Nf-1) = E2f(Nf,Nf-1) + w1*E2c(N,N)
    E2f(Nf-1,Nf) = E2f(Nf-1,Nf) + w1*E2c(N,N)
    E2f(Nf,Nf) = E2f(Nf,Nf) + w1*E2c(N,N)

    ! lesser contribution to intermediate
    E2f(Nf-2,Nf-1) = E2f(Nf-2,Nf-1) + w2*E2c(N,N)
    E2f(Nf-2,Nf) = E2f(Nf-2,Nf) + w2*E2c(N,N)
    E2f(Nf-1,Nf-2) = E2f(Nf-1,Nf-2) + w2*E2c(N,N)
    E2f(Nf-1,1) = E2f(Nf-1,1) + w2*E2c(N,N)
    E2f(Nf,Nf-2) = E2f(Nf,Nf-2) + w2*E2c(N,N)
    E2f(Nf,1) = E2f(Nf,1) + w2*E2c(N,N)
    E2f(1,Nf-1) = E2f(1,Nf-1) + w2*E2c(N,N)
    E2f(1,Nf) = E2f(1,Nf) + w2*E2c(N,N)

    ! least contribution to furthest
    E2f(Nf-2,Nf-2) = E2f(Nf-2,Nf-2) + w3*E2c(N,N)
    E2f(Nf-2,1) = E2f(Nf-2,1) + w3*E2c(N,N)
    E2f(1,Nf-2) = E2f(1,Nf-2) + w3*E2c(N,N)
    E2f(1,1) = E2f(1,1) + w3*E2c(N,N)
  end subroutine prolongate0
end module fd_solvers
