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

    allocate(mg(0:level))

    do i=0,level
      allocate(mg(i)%grid((2**i),(2**i)))
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

    do i=0,level
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

  implicit none
  save

  ! solver parameters
  real(dp), parameter :: tau = 2.0_dp
  integer, parameter :: niter = 11

  private tau, niter, invert, laplacian, compute_g, vcycle, smooth, restrict, prolongate
  public solver_ufds2t2

  contains


  !> @brief solves the dimensionless CH equation
  !!
  !! @param[in] Tout  output times
  !! @param[in] c0    initial concentration
  !! @param[in] eps2  dimensionless number
  subroutine solver_ufds2t2(Tout, CH_params, c0, eps2, errors)
    implicit none
    real(dp), intent(in) :: Tout(:)
    real(dp), intent(in) :: eps2
    real(dp), intent(in), dimension(6) :: CH_params
    real(dp), dimension(:,:), allocatable, intent(in) :: c0
    integer :: errors

    integer :: N ! grid size
    integer :: level ! grid level
    real(dp) :: dx ! grid spacing
    real(dp) :: dt, dt0, dt1, dt_out ! timesteps
    real(dp) :: a0, a1, a2, b0, b1 ! time constants
    real(dp) :: t, t_out,  tmax
    real(dp), dimension(2,2) :: A ! smoothing matrix
    real(dp) :: eps
    integer :: it, i, j ! iterators
    logical :: outflag
    character(len=48) :: msg ! logging message

    ! grid storage
    real(dp), dimension(:,:), allocatable :: phi, psi, g, b, phi_prev, g_prev, work
    real(dp), dimension(:,:), allocatable :: c, c_prev
    type(t_grid), dimension(:), allocatable :: E1, E2, R1, R2


    ! ======================================================================== !
    !   SETUP                                                                  !
    ! ======================================================================== !
    ! set variables
    N = size(c0,1)
    call ilog2(N,level)
    dx = 1.0_dp/(real(N,dp))
    dt = 2.5_dp * eps2
    tmax = maxval(tout)
    t = 0.0_dp
    it = 1
    eps = sqrt(eps2)

    24 format(A, F7.3) ! output message

    ! allocate storage
    allocate(phi(N,N))
    allocate(psi(N,N))
    allocate(g(N,N))
    allocate(b(N,N))
    allocate(phi_prev(N,N))
    allocate(g_prev(N,N))
    allocate(work(N,N))
    allocate(c(N, N))
    allocate(c_prev(N, N))

    ! allocate multigrid storage
    call multigrid_alloc(E1, level)
    call multigrid_alloc(E2, level)
    call multigrid_alloc(R1, level)
    call multigrid_alloc(R2, level)

    ! initial condition
    do j=1,N
      do i=1,N
        phi(i,j) = c0(i,j)
      enddo
    enddo

    ! set coupled variable
    call laplacian(phi, psi, dx, N)
    psi = tau*phi - eps2*psi

    !! TO OUT 
    ! output if required
    if (tout(it) < epsilon(tout(it))) then
      write(msg, 24) "Initial condition output at t=  0.000"
      call logger%info("solver_ufds2t2", msg)
      dt_out = dt
      t_out = t
      c = phi
      c_prev = phi_prev
      call dimensionalise(CH_params, c, t_out)
      call dimensionalise(CH_params, c_prev, dt_out)

      call write_to_traj(c, c_prev, t_out, dt_out, errors)
      it = it + 1
    endif


    ! ========================================================================== !
    !   FIRST TIMESTEP (first order)                                             !
    ! ========================================================================== !
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
      ! compute residuals
      call laplacian(psi, R1(level)%grid, dx, n)
      R1(level)%grid = R1(level)%grid - phi/dt + b
      call laplacian(phi, R2(level)%grid, dx, n)
      R2(level)%grid = tau*phi-eps2*R2(level)%grid - psi

      ! TODO: finish iteration conditional on the size of the residual

      ! perform a single v-cycle
      call vcycle(A, E1, E2, R1, R2, eps2, N, dx, level)

      ! update with errors
      phi = phi + E1(level)%grid
      psi = psi + E2(level)%grid
    enddo

    !! TO OUT 
    ! conditionally output
    if (outflag) then
      write(msg, 24) "Output at t=", t
      call logger%info("solver_ufds2t2", msg)

      dt_out = dt
      t_out = t
      c = phi
      c_prev = phi_prev
      call dimensionalise(CH_params, c, t_out)
      call dimensionalise(CH_params, c_prev, dt_out)

      call write_to_traj(c, c_prev, t_out, dt_out, errors)

      it = it + 1
    endif


    ! ========================================================================== !
    !   REMAINING TIMESTEPS (second order)                                       !
    ! ========================================================================== !
    do while (t < tmax)
      ! set current timestep TODO: condition on curvature rather than time
      if (t < 100.0_dp * eps) then
        dt = 2.5_dp * eps2
      else
        dt = 2.5_dp * eps2
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

        ! TODO: finish iteration conditional on the size of the residual

        ! perform a single v-cycle
        call vcycle(A, E1, E2, R1, R2, eps2, N, dx, level)

        ! update with errors
        phi = phi + E1(level)%grid
        psi = psi + E2(level)%grid

        ! print size of residual
        ! print *, i, "r1 max = ", maxval(abs(R1(level)%grid))
      enddo

      !! TO OUT 
      ! conditionally output
      if (outflag) then
        write(msg, 24) "Output at t=", t
        call logger%info("solver_ufds2t2", msg)
        dt_out = dt
        t_out = t
        c = phi
        c_prev = phi_prev
        call dimensionalise(CH_params, c, t_out)
        call dimensionalise(CH_params, c_prev, dt_out)
  
        call write_to_traj(c, c_prev, t_out, dt_out, errors)
        
        it = it + 1
      endif
    enddo

    ! ========================================================================== !
    !   CLEAN UP                                                                 !
    ! ========================================================================== !
    deallocate(phi)
    deallocate(psi)
    deallocate(g)
    deallocate(b)
    deallocate(phi_prev)
    deallocate(g_prev)
    deallocate(work)
    call multigrid_dealloc(E1, level)
    call multigrid_dealloc(E2, level)
    call multigrid_dealloc(R1, level)
    call multigrid_dealloc(R2, level)
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
    real(dp), dimension(:,:), allocatable, intent(in) :: x
    real(dp), dimension(:,:), allocatable, intent(inout) :: y
    real(dp), intent(in) :: dx
    integer, intent(in) :: n
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
    real(dp), dimension(:,:), allocatable, intent(in) :: phi
    real(dp), dimension(:,:), allocatable, intent(inout) :: work
    real(dp), dimension(:,:), allocatable, intent(inout) :: g
    real(dp), intent(in) :: dx
    integer, intent(in) :: n

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
  subroutine vcycle(A, E1, E2, R1, R2, eps2, N, dx, level)
    implicit none
    real(dp), dimension(2,2), intent(in) :: A
    type(t_grid), dimension(:), allocatable, intent(inout) :: E1, E2, R1, R2
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

      call smooth(A, E1(l)%grid, E2(l)%grid, R1(l)%grid, R2(l)%grid, eps2, nl, dxl)
      call restrict(R1(l)%grid, R2(l)%grid, R1(l-1)%grid, R2(l-1)%grid, nl)

      nl = nl/2;
      dxl = dxl * 2.0_dp
    enddo

    ! smooth at level 1
    E1(l)%grid = 0.0_dp
    E2(l)%grid = 0.0_dp
    call smooth(A, E1(1)%grid, E2(1)%grid, R1(1)%grid, R2(1)%grid, eps2, nl, dxl)
    call smooth(A, E1(1)%grid, E2(1)%grid, R1(1)%grid, R2(1)%grid, eps2, nl, dxl)
    call smooth(A, E1(1)%grid, E2(1)%grid, R1(1)%grid, R2(1)%grid, eps2, nl, dxl)
    call smooth(A, E1(1)%grid, E2(1)%grid, R1(1)%grid, R2(1)%grid, eps2, nl, dxl)

    ! go down, smoothing and prolongating
    do l=1,level-1
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
    real(dp), dimension(2,2), intent(in) :: A
    real(dp), dimension(:,:), allocatable, intent(inout) :: E1, E2, R1, R2
    real(dp), intent(in) :: eps2, dx
    integer, intent(in) :: N

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
  end subroutine smooth


  !> @brief restricts (ie coarsens) to level-1 from level
  !!
  !! @param R1     residual multigrid for the first variable
  !! @param R2     residual multigrid for the second variable
  !! @param N      grid size
  subroutine restrict(R1f, R2f, R1c, R2c, N)
    implicit none
    real(dp), dimension(:,:), allocatable, intent(inout) :: R1f, R2f, R1c, R2c
    integer, intent(in) :: N

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
    real(dp), dimension(:,:), allocatable, intent(inout) :: E1f, E2f, E1c, E2c
    integer, intent(in) :: N

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
  end subroutine prolongate
end module fd_solvers
