module multigrid
  use globals

  implicit none
  save

  !> @brief multigrid internal grid type
  !!
  !! This is just a wrapper around an array of reals. An array of these forms a
  !! multigrid
  type :: t_grid
    real(dp), dimension(:), allocatable :: grid
  end type

  contains


  !> @brief Allocates storage for a multigrid with a given level
  !!
  !! @param mg[inout]  multigrid
  !! @param level[in]  maximum level
  subroutine multigrid_alloc(mg, level)
    implicit none
    type(t_grid), dimension(:), allocatable, intent(inout) :: mg(:)
    integer, intent(in) :: level
    integer :: i

    allocate(mg(0:level))

    do i=0,level
      allocate(mg(i)%grid(2**(2*i)))
    enddo
  end subroutine multigrid_alloc


  !> @brief Deallocates storage for a multigrid with a given level
  !!
  !! @param mg[inout]  multigrid
  !! @param level[in]  maximum level
  subroutine multigrid_dealloc(mg, level)
    implicit none
    type(t_grid), dimension(:), allocatable, intent(inout) :: mg(:)
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
    real(dp), dimension(:), allocatable :: phi, psi, g, b, phi_prev, g_prev, work
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
    allocate(phi(N*N))
    allocate(psi(N*N))
    allocate(g(N*N))
    allocate(b(N*N))
    allocate(phi_prev(N*N))
    allocate(g_prev(N*N))
    allocate(work(N*N))
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
        phi(i+n*(j-1)) = c0(i,j)
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
      do j=1,n
        do i=1,n
          c(j,i) = phi(i+n*(j-1))
          c_prev(j,i) = phi_prev(i+n*(j-1))
        enddo
      enddo
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
      print *, ""
      print *, ""
      call vcycle_flat(A, E1, E2, R1, R2, eps2, N, dx, level)
      stop

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
      do j=1,n
        do i=1,n
          c(j,i) = phi(i+n*(j-1))
          c_prev(j,i) = phi_prev(i+n*(j-1))
        enddo
      enddo
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
        do j=1,n
          do i=1,n
            c(j,i) = phi(i+n*(j-1))
            c_prev(j,i) = phi_prev(i+n*(j-1))
          enddo
        enddo
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
    real(dp), dimension(:), allocatable, intent(in) :: x
    real(dp), dimension(:), allocatable, intent(inout) :: y
    real(dp), intent(in) :: dx
    integer, intent(in) :: n
    integer :: i, j, ij
    real(dp) :: dx2_ ! interim constants

    dx2_ = 1.0_dp / (dx*dx)

    ! interior
    do j=2,n-1
      do i=2,n-1
        ij = i+n*(j-1)
        y(ij) = dx2_*(x(ij+1) + x(ij-1) + x(ij+n) + x(ij-n) - 4*x(ij))
      enddo
    enddo

    ! left/right
    do j=2,n-1
      ij = 1+n*(j-1)
      y(ij) = dx2_*(x(ij+1) + x(ij+n-1) + x(ij+n) + x(ij-n) - 4*x(ij))

      ij = n+n*(j-1)
      y(ij) = dx2_*(x(ij-n+1) + x(ij-1) + x(ij+n) + x(ij-n) - 4*x(ij))
    enddo

    ! top/bottom
    do i=2,n-1
      ij = i+n*(1-1)
      y(ij) = dx2_*(x(ij+1) + x(ij-1) + x(ij+n) + x(ij+n*(n-1)) - 4*x(ij))

      ij = i+n*(n-1)
      y(ij) = dx2_*(x(ij+1) + x(ij-1) + x(ij-n*(n-1)) + x(ij-n) - 4*x(ij))
    enddo

    ! corners
    ij = 1+n*(1-1)
    y(ij) = dx2_*(x(ij+1) + x(ij+n-1) + x(ij+n) + x(ij+n*(n-1)) - 4*x(ij))

    ij = 1+n*(n-1)
    y(ij) = dx2_*(x(ij+1) + x(ij+n-1) + x(ij-n*(n-1)) + x(ij-n) - 4*x(ij))

    ij = n+n*(1-1)
    y(ij) = dx2_*(x(ij-n+1) + x(ij-1) + x(ij+n) + x(ij+n*(n-1)) - 4*x(ij))

    ij = n+n*(n-1)
    y(ij) = dx2_*(x(ij-n+1) + x(ij-1) + x(ij-n*(n-1)) + x(ij-n) - 4*x(ij))
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
    real(dp), dimension(:), allocatable, intent(in) :: phi
    real(dp), dimension(:), allocatable, intent(inout) :: work
    real(dp), dimension(:), allocatable, intent(inout) :: g
    real(dp), intent(in) :: dx
    integer, intent(in) :: n
    integer :: i

    do i=1,n*n
      work(i) = phi(i) * (phi(i)*phi(i) - (1+tau))
    enddo

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
  recursive subroutine vcycle(A, E1, E2, R1, R2, eps2, N, dx, level)
    implicit none
    real(dp), dimension(2,2), intent(in) :: A
    type(t_grid), dimension(:), allocatable, intent(inout) :: E1, E2, R1, R2
    real(dp), intent(in) :: eps2, dx
    integer, intent(in) :: N, level

    print "(A, I0, A, I0)", "call vcycle level ", level, ", N = ", N
    E1(level)%grid = 0.0_dp
    E2(level)%grid = 0.0_dp

    ! start smooth
    call smooth(A, E1(level)%grid, E2(level)%grid, R1(level)%grid, R2(level)%grid, eps2, N, dx)

    ! if the level we're at is greater than 0 then a smooth won't solve it
    if (level > 1) then
      ! restrict (level -> level-1)
      call restrict(R1, R2, N, level)

      ! recurse
      call vcycle(A, E1, E2, R1, R2, eps2, N/2, dx*2.0_dp, level-1)

      ! prolongate (level-1 -> level) and correct
      call prolongate(E1, E2, N/2, level-1)
    else
      ! extra smooths on the 2x2 grid
      call smooth(A, E1(level)%grid, E2(level)%grid, R1(level)%grid, R2(level)%grid, eps2, N, dx)
      call smooth(A, E1(level)%grid, E2(level)%grid, R1(level)%grid, R2(level)%grid, eps2, N, dx)
    endif

    ! start smooth
    call smooth(A, E1(level)%grid, E2(level)%grid, R1(level)%grid, R2(level)%grid, eps2, N, dx)
  end subroutine vcycle

  !! flat version of vcycle
  subroutine vcycle_flat(A, E1, E2, R1, R2, eps2, N, dx, level)
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
      print "(A, I0, A, I0)", "call vcycle level ", l, ", N = ", nl
      E1(l)%grid = 0.0_dp
      E2(l)%grid = 0.0_dp

      call smooth(A, E1(l)%grid, E2(l)%grid, R1(l)%grid, R2(l)%grid, eps2, nl, dxl)
      call restrict(R1, R2, nl, l)

      nl = nl/2;
      dxl = dxl * 2.0_dp
    enddo

    ! smooth at level 1
    print "(A, I0, A, I0)", "call vcycle level ", l, ", N = ", nl
    E1(l)%grid = 0.0_dp
    E2(l)%grid = 0.0_dp
    call smooth(A, E1(1)%grid, E2(1)%grid, R1(1)%grid, R2(1)%grid, eps2, nl, dxl)
    call smooth(A, E1(1)%grid, E2(1)%grid, R1(1)%grid, R2(1)%grid, eps2, nl, dxl)
    call smooth(A, E1(1)%grid, E2(1)%grid, R1(1)%grid, R2(1)%grid, eps2, nl, dxl)
    call smooth(A, E1(1)%grid, E2(1)%grid, R1(1)%grid, R2(1)%grid, eps2, nl, dxl)

    ! go down, smoothing and prolongating
    do l=1,level
      call prolongate(E1, E2, nl, l)

      nl = nl*2;
      dxl = dxl * 0.5_dp

      call smooth(A, E1(l)%grid, E2(l)%grid, R1(l)%grid, R2(l)%grid, eps2, nl, dxl)
    enddo
  end subroutine vcycle_flat


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
  !! @param level  grid level
  subroutine smooth(A, E1, E2, R1, R2, eps2, N, dx)
    implicit none
    real(dp), dimension(2,2), intent(in) :: A
    real(dp), dimension(:), allocatable, intent(inout) :: E1, E2, R1, R2
    real(dp), intent(in) :: eps2, dx
    integer, intent(in) :: N

    real(dp), dimension(2) :: rhs
    real(dp) :: dx2_
    integer :: i, j, ij, shift

    call ilog2(N, i)
    print "(A, I0, A, I0)", "call smooth level ", i, ", N = ", N

    dx2_ = 1.0_dp / (dx*dx)

    ! ================ !
    ! SMOOTH RED NODES !
    ! ================ !
    ! interior
    do j=2,n-1
      shift = mod(j,2)
      do i=2+shift,n-1+shift,2
        ij = i+n*(j-1)

        ! compute RHS
        rhs(1) = R1(ij) + dx2_ * (E2(ij+1) + E2(ij-1) + E2(ij+n) + E2(ij-n))
        rhs(2) = R2(ij) - eps2*dx2_ * (E1(ij+1) + E1(ij-1) + E1(ij+n) + E1(ij-n))

        ! solve for new errors
        E1(ij) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
        E2(ij) = A(2,1)*rhs(1) + A(2,2)*rhs(2)
      enddo
    enddo

    ! left/right
    do j=2,n-1,2
      ij = 1+n*j
      rhs(1) = R1(ij) + dx2_ * (E2(ij+1) + E2(ij+n-1) + E2(ij+n) + E2(ij-n))
      rhs(2) = R2(ij) - eps2*dx2_ * (E1(ij+1) + E1(ij+n-1) + E1(ij+n) + E1(ij-n))
      E1(ij) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
      E2(ij) = A(2,1)*rhs(1) + A(2,2)*rhs(2)

      ij = n+n*(j-1)
      rhs(1) = R1(ij) + dx2_ * (E2(ij-n+1) + E2(ij-1) + E2(ij+n) + E2(ij-n))
      rhs(2) = R2(ij) - eps2*dx2_ * (E1(ij-n+1) + E1(ij-1) + E1(ij+n) + E1(ij-n))
      E1(ij) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
      E2(ij) = A(2,1)*rhs(1) + A(2,2)*rhs(2)
    enddo

    ! top/bottom
    do i=2,n-1,2
      ij = i+1+n*(1-1)
      rhs(1) = R1(ij) + dx2_ * (E2(ij+1) + E2(ij-1) + E2(ij+n) + E2(ij+n*(n-1)))
      rhs(2) = R2(ij) - eps2*dx2_ * (E1(ij+1) + E1(ij-1) + E1(ij+n) + E1(ij+n*(n-1)))
      E1(ij) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
      E2(ij) = A(2,1)*rhs(1) + A(2,2)*rhs(2)

      ij = i+n*(n-1)
      rhs(1) = R1(ij) + dx2_ * (E2(ij+1) + E2(ij-1) + E2(ij-n*(n-1)) + E2(ij-n))
      rhs(2) = R2(ij) - eps2*dx2_ * (E1(ij+1) + E1(ij-1) + E1(ij-n*(n-1)) + E1(ij-n))
      E1(ij) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
      E2(ij) = A(2,1)*rhs(1) + A(2,2)*rhs(2)
    enddo

    ! upper left corner
    ij = 1+n*(1-1)
    rhs(1) = R1(ij) + dx2_ * (E2(ij+1) + E2(ij+n-1) + E2(ij+n) + E2(ij+n*(n-1)))
    rhs(2) = R2(ij) - eps2*dx2_ * (E1(ij+1) + E1(ij+n-1) + E1(ij+n) + E1(ij+n*(n-1)))
    E1(ij) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
    E2(ij) = A(2,1)*rhs(1) + A(2,2)*rhs(2)

    ! lower right corner
    ij = n+n*(n-1)
    rhs(1) = R1(ij) + dx2_ * (E2(ij-n+1) + E2(ij-1) + E2(ij-n*(n-1)) + E2(ij-n))
    rhs(2) = R2(ij) - eps2*dx2_ * (E1(ij-n+1) + E1(ij-1) + E1(ij-n*(n-1)) + E1(ij-n))
    E1(ij) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
    E2(ij) = A(2,1)*rhs(1) + A(2,2)*rhs(2)


    ! ================== !
    ! SMOOTH BLACK NODES !
    ! ================== !
    ! interior
    do j=2,n-1
      shift = mod(j-1,2)
      do i=2+shift,n-1+shift,2
        ij = i+n*(j-1)

        ! compute RHS
        rhs(1) = R1(ij) + dx2_ * (E2(ij+1) + E2(ij-1) + E2(ij+n) + E2(ij-n))
        rhs(2) = R2(ij) - eps2*dx2_ * (E1(ij+1) + E1(ij-1) + E1(ij+n) + E1(ij-n))

        ! solve for new errors
        E1(ij) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
        E2(ij) = A(2,1)*rhs(1) + A(2,2)*rhs(2)
      enddo
    enddo

    ! left/right
    do j=2,n-1,2
      ij = 1+n*(j-1)
      rhs(1) = R1(ij) + dx2_ * (E2(ij+1) + E2(ij+n-1) + E2(ij+n) + E2(ij-n))
      rhs(2) = R2(ij) - eps2*dx2_ * (E1(ij+1) + E1(ij+n-1) + E1(ij+n) + E1(ij-n))
      E1(ij) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
      E2(ij) = A(2,1)*rhs(1) + A(2,2)*rhs(2)

      ij = n+n*j
      rhs(1) = R1(ij) + dx2_ * (E2(ij-n+1) + E2(ij-1) + E2(ij+n) + E2(ij-n))
      rhs(2) = R2(ij) - eps2*dx2_ * (E1(ij-n+1) + E1(ij-1) + E1(ij+n) + E1(ij-n))
      E1(ij) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
      E2(ij) = A(2,1)*rhs(1) + A(2,2)*rhs(2)
    enddo

    ! top/bottom
    do i=2,n-1,2
      ij = i+n*(1-1)
      rhs(1) = R1(ij) + dx2_ * (E2(ij+1) + E2(ij-1) + E2(ij+n) + E2(ij+n*(n-1)))
      rhs(2) = R2(ij) - eps2*dx2_ * (E1(ij+1) + E1(ij-1) + E1(ij+n) + E1(ij+n*(n-1)))
      E1(ij) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
      E2(ij) = A(2,1)*rhs(1) + A(2,2)*rhs(2)

      ij = i+1+n*(n-1)
      rhs(1) = R1(ij) + dx2_ * (E2(ij+1) + E2(ij-1) + E2(ij-n*(n-1)) + E2(ij-n))
      rhs(2) = R2(ij) - eps2*dx2_ * (E1(ij+1) + E1(ij-1) + E1(ij-n*(n-1)) + E1(ij-n))
      E1(ij) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
      E2(ij) = A(2,1)*rhs(1) + A(2,2)*rhs(2)
    enddo

    ! upper right
    ij = 1+n*(n-1)
    rhs(1) = R1(ij) + dx2_ * (E2(ij+1) + E2(ij+n-1) + E2(ij-n*(n-1)) + E2(ij-n))
    rhs(2) = R2(ij) - eps2*dx2_ * (E1(ij+1) + E1(ij+n-1) + E1(ij-n*(n-1)) + E1(ij-n))
    E1(ij) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
    E2(ij) = A(2,1)*rhs(1) + A(2,2)*rhs(2)

    ! lower left
    ij = n+n*(1-1)
    rhs(1) = R1(ij) + dx2_ * (E2(ij-n+1) + E2(ij-1) + E2(ij+n) + E2(ij+n*(n-1)))
    rhs(2) = R2(ij) - eps2*dx2_ * (E1(ij-n+1) + E1(ij-1) + E1(ij+n) + E1(ij+n*(n-1)))
    E1(ij) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
    E2(ij) = A(2,1)*rhs(1) + A(2,2)*rhs(2)
  end subroutine smooth


  !> @brief restricts (ie coarsens) to level-1 from level
  !!
  !! @param R1     residual multigrid for the first variable
  !! @param R2     residual multigrid for the second variable
  !! @param N      grid size
  !! @param level  grid level
  subroutine restrict(R1, R2, N, level)
    implicit none
    type(t_grid), dimension(:), allocatable, intent(inout) :: R1, R2
    integer, intent(in) :: N, level

    integer :: i, j, ij, im, jm, ijm

    print "(A, I0, A, I0)", "call restrict level ", level, ", N = ", N

    ! use a simple average across the grid
    do j=1,2**(level-1)
      jm = 2*j-1
      do i=1,2**(level-1)
        im = 2*i-1

        ij = i+(n/2)*(j-1)
        ijm = im+n*(jm-1)

        ! fine -> coarse
        R1(level-1)%grid(ij) = 0.25_dp*(  R1(level)%grid(ijm) &
                                        + R1(level)%grid(ijm+1) &
                                        + R1(level)%grid(ijm+N) &
                                        + R1(level)%grid(ijm+N+1) )


        R2(level-1)%grid(ij) = 0.25_dp*(  R2(level)%grid(ijm) &
                                        + R2(level)%grid(ijm+1) &
                                        + R2(level)%grid(ijm+N) &
                                        + R2(level)%grid(ijm+N+1) )
      enddo
    enddo
  end subroutine restrict


  !> @brief prolongates (ie refines) from level to level+1 (bilinear interpolation)
  !!
  !! @param E1     error multigrid for the first variable
  !! @param E2     error multigrid for the second variable
  !! @param N      grid size
  !! @param level  grid level
  subroutine prolongate(E1, E2, N, level)
    implicit none
    type(t_grid), dimension(:), allocatable, intent(inout) :: E1, E2
    integer, intent(in) :: N, level

    integer :: Nf
    integer :: i, j, ij, if, jf, ijf
    real(dp), parameter :: w1 = 0.5625_dp
    real(dp), parameter :: w2 = 0.1875_dp
    real(dp), parameter :: w3 = 0.0625_dp

    print "(A, I0, A, I0)", "call prolongate level ", level, ", N = ", N

    Nf = 2*N

    ! interior
    do j=2,N-1
      jf = 2*j-1
      do i=2,N-1
        if = 2*i-1

        ij = i+N*(j-1)
        ijf = if+Nf*(jf-1)

        ! largest contribution to nearest
        E1(level+1)%grid(ijf) = E1(level+1)%grid(ijf) + w1*E1(level)%grid(ij)
        E1(level+1)%grid(ijf+1) = E1(level+1)%grid(ijf+1) + w1*E1(level)%grid(ij)
        E1(level+1)%grid(ijf+Nf) = E1(level+1)%grid(ijf+Nf) + w1*E1(level)%grid(ij)
        E1(level+1)%grid(ijf+Nf+1) = E1(level+1)%grid(ijf+Nf+1) + w1*E1(level)%grid(ij)

        ! lesser contribution to intermediate
        E1(level+1)%grid(ijf-1) = E1(level+1)%grid(ijf-1) + w2*E1(level)%grid(ij)
        E1(level+1)%grid(ijf+Nf-1) = E1(level+1)%grid(ijf+Nf-1) + w2*E1(level)%grid(ij)
        E1(level+1)%grid(ijf-Nf) = E1(level+1)%grid(ijf-Nf) + w2*E1(level)%grid(ij)
        E1(level+1)%grid(ijf+2*Nf) = E1(level+1)%grid(ijf+2*Nf) + w2*E1(level)%grid(ij)
        E1(level+1)%grid(ijf+1-Nf) = E1(level+1)%grid(ijf+1-Nf) + w2*E1(level)%grid(ij)
        E1(level+1)%grid(ijf+1+2*Nf) = E1(level+1)%grid(ijf+1+2*Nf) + w2*E1(level)%grid(ij)
        E1(level+1)%grid(ijf+2) = E1(level+1)%grid(ijf+2) + w2*E1(level)%grid(ij)
        E1(level+1)%grid(ijf+2+Nf) = E1(level+1)%grid(ijf+2+Nf) + w2*E1(level)%grid(ij)

        ! least contribution to furthest
        E1(level+1)%grid(ijf-1-Nf) = E1(level+1)%grid(ijf-1-Nf) + w3*E1(level)%grid(ij)
        E1(level+1)%grid(ijf-1+2*Nf) = E1(level+1)%grid(ijf-1+2*Nf) + w3*E1(level)%grid(ij)
        E1(level+1)%grid(ijf+2-Nf) = E1(level+1)%grid(ijf+2-Nf) + w3*E1(level)%grid(ij)
        E1(level+1)%grid(ijf+2+2*Nf) = E1(level+1)%grid(ijf+2+2*Nf) + w3*E1(level)%grid(ij)

        ! largest contribution to nearest
        E2(level+1)%grid(ijf) = E2(level+1)%grid(ijf) + w1*E2(level)%grid(ij)
        E2(level+1)%grid(ijf+1) = E2(level+1)%grid(ijf+1) + w1*E2(level)%grid(ij)
        E2(level+1)%grid(ijf+Nf) = E2(level+1)%grid(ijf+Nf) + w1*E2(level)%grid(ij)
        E2(level+1)%grid(ijf+Nf+1) = E2(level+1)%grid(ijf+Nf+1) + w1*E2(level)%grid(ij)

        ! lesser contribution to intermediate
        E2(level+1)%grid(ijf-1) = E2(level+1)%grid(ijf-1) + w2*E2(level)%grid(ij)
        E2(level+1)%grid(ijf+Nf-1) = E2(level+1)%grid(ijf+Nf-1) + w2*E2(level)%grid(ij)
        E2(level+1)%grid(ijf-Nf) = E2(level+1)%grid(ijf-Nf) + w2*E2(level)%grid(ij)
        E2(level+1)%grid(ijf+2*Nf) = E2(level+1)%grid(ijf+2*Nf) + w2*E2(level)%grid(ij)
        E2(level+1)%grid(ijf+1-Nf) = E2(level+1)%grid(ijf+1-Nf) + w2*E2(level)%grid(ij)
        E2(level+1)%grid(ijf+1+2*Nf) = E2(level+1)%grid(ijf+1+2*Nf) + w2*E2(level)%grid(ij)
        E2(level+1)%grid(ijf+2) = E2(level+1)%grid(ijf+2) + w2*E2(level)%grid(ij)
        E2(level+1)%grid(ijf+2+Nf) = E2(level+1)%grid(ijf+2+Nf) + w2*E2(level)%grid(ij)

        ! least contribution to furthest
        E2(level+1)%grid(ijf-1-Nf) = E2(level+1)%grid(ijf-1-Nf) + w3*E2(level)%grid(ij)
        E2(level+1)%grid(ijf-1+2*Nf) = E2(level+1)%grid(ijf-1+2*Nf) + w3*E2(level)%grid(ij)
        E2(level+1)%grid(ijf+2-Nf) = E2(level+1)%grid(ijf+2-Nf) + w3*E2(level)%grid(ij)
        E2(level+1)%grid(ijf+2+2*Nf) = E2(level+1)%grid(ijf+2+2*Nf) + w3*E2(level)%grid(ij)
      enddo
    enddo

    ! left edge
    i = 1
    if = 2*i-1
    do j=2,N-1
      jf = 2*j-1

      ij = i+N*(j-1)
      ijf = if+Nf*(jf-1)

      ! largest contribution to nearest
      E1(level+1)%grid(ijf) = E1(level+1)%grid(ijf) + w1*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+1) = E1(level+1)%grid(ijf+1) + w1*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+Nf) = E1(level+1)%grid(ijf+Nf) + w1*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+Nf+1) = E1(level+1)%grid(ijf+Nf+1) + w1*E1(level)%grid(ij)

      ! lesser contribution to intermediate
      E1(level+1)%grid(ijf-1+Nf) = E1(level+1)%grid(ijf-1+Nf) + w2*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+2*Nf-1) = E1(level+1)%grid(ijf+2*Nf-1) + w2*E1(level)%grid(ij)
      E1(level+1)%grid(ijf-Nf) = E1(level+1)%grid(ijf-Nf) + w2*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+2*Nf) = E1(level+1)%grid(ijf+2*Nf) + w2*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+1-Nf) = E1(level+1)%grid(ijf+1-Nf) + w2*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+1+2*Nf) = E1(level+1)%grid(ijf+1+2*Nf) + w2*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+2) = E1(level+1)%grid(ijf+2) + w2*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+2+Nf) = E1(level+1)%grid(ijf+2+Nf) + w2*E1(level)%grid(ij)

      ! least contribution to furthest
      E1(level+1)%grid(ijf-1) = E1(level+1)%grid(ijf-1) + w3*E1(level)%grid(ij)
      E1(level+1)%grid(ijf-1+3*Nf) = E1(level+1)%grid(ijf-1+3*Nf) + w3*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+2-Nf) = E1(level+1)%grid(ijf+2-Nf) + w3*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+2+2*Nf) = E1(level+1)%grid(ijf+2+2*Nf) + w3*E1(level)%grid(ij)

      ! largest contribution to nearest
      E2(level+1)%grid(ijf) = E2(level+1)%grid(ijf) + w1*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+1) = E2(level+1)%grid(ijf+1) + w1*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+Nf) = E2(level+1)%grid(ijf+Nf) + w1*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+Nf+1) = E2(level+1)%grid(ijf+Nf+1) + w1*E2(level)%grid(ij)

      ! lesser contribution to intermediate
      E2(level+1)%grid(ijf-1+Nf) = E2(level+1)%grid(ijf-1+Nf) + w2*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+2*Nf-1) = E2(level+1)%grid(ijf+2*Nf-1) + w2*E2(level)%grid(ij)
      E2(level+1)%grid(ijf-Nf) = E2(level+1)%grid(ijf-Nf) + w2*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+2*Nf) = E2(level+1)%grid(ijf+2*Nf) + w2*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+1-Nf) = E2(level+1)%grid(ijf+1-Nf) + w2*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+1+2*Nf) = E2(level+1)%grid(ijf+1+2*Nf) + w2*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+2) = E2(level+1)%grid(ijf+2) + w2*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+2+Nf) = E2(level+1)%grid(ijf+2+Nf) + w2*E2(level)%grid(ij)

      ! least contribution to furthest
      E2(level+1)%grid(ijf-1) = E2(level+1)%grid(ijf-1) + w3*E2(level)%grid(ij)
      E2(level+1)%grid(ijf-1+3*Nf) = E2(level+1)%grid(ijf-1+3*Nf) + w3*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+2-Nf) = E2(level+1)%grid(ijf+2-Nf) + w3*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+2+2*Nf) = E2(level+1)%grid(ijf+2+2*Nf) + w3*E2(level)%grid(ij)
    enddo

    ! right edge
    i = N
    if = 2*i-1
    do j=2,N-1
      jf = 2*j-1

      ij = i+N*(j-1)
      ijf = if+Nf*(jf-1)

      ! largest contribution to nearest
      E1(level+1)%grid(ijf) = E1(level+1)%grid(ijf) + w1*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+1) = E1(level+1)%grid(ijf+1) + w1*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+Nf) = E1(level+1)%grid(ijf+Nf) + w1*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+Nf+1) = E1(level+1)%grid(ijf+Nf+1) + w1*E1(level)%grid(ij)

      ! lesser contribution to intermediate
      E1(level+1)%grid(ijf-1) = E1(level+1)%grid(ijf-1) + w2*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+Nf-1) = E1(level+1)%grid(ijf+Nf-1) + w2*E1(level)%grid(ij)
      E1(level+1)%grid(ijf-Nf) = E1(level+1)%grid(ijf-Nf) + w2*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+2*Nf) = E1(level+1)%grid(ijf+2*Nf) + w2*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+1-Nf) = E1(level+1)%grid(ijf+1-Nf) + w2*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+1+2*Nf) = E1(level+1)%grid(ijf+1+2*Nf) + w2*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+2-Nf) = E1(level+1)%grid(ijf+2-Nf) + w2*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+2) = E1(level+1)%grid(ijf+2) + w2*E1(level)%grid(ij)

      ! least contribution to furthest
      E1(level+1)%grid(ijf-1-Nf) = E1(level+1)%grid(ijf-1-Nf) + w3*E1(level)%grid(ij)
      E1(level+1)%grid(ijf-1+2*Nf) = E1(level+1)%grid(ijf-1+2*Nf) + w3*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+2-2*Nf) = E1(level+1)%grid(ijf+2-2*Nf) + w3*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+2*Nf) = E1(level+1)%grid(ijf+2+Nf) + w3*E1(level)%grid(ij)

      ! largest contribution to nearest
      E2(level+1)%grid(ijf) = E2(level+1)%grid(ijf) + w1*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+1) = E2(level+1)%grid(ijf+1) + w1*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+Nf) = E2(level+1)%grid(ijf+Nf) + w1*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+Nf+1) = E2(level+1)%grid(ijf+Nf+1) + w1*E2(level)%grid(ij)

      ! lesser contribution to intermediate
      E2(level+1)%grid(ijf-1) = E2(level+1)%grid(ijf-1) + w2*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+Nf-1) = E2(level+1)%grid(ijf+Nf-1) + w2*E2(level)%grid(ij)
      E2(level+1)%grid(ijf-Nf) = E2(level+1)%grid(ijf-Nf) + w2*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+2*Nf) = E2(level+1)%grid(ijf+2*Nf) + w2*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+1-Nf) = E2(level+1)%grid(ijf+1-Nf) + w2*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+1+2*Nf) = E2(level+1)%grid(ijf+1+2*Nf) + w2*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+2-Nf) = E2(level+1)%grid(ijf+2-Nf) + w2*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+2) = E2(level+1)%grid(ijf+2) + w2*E2(level)%grid(ij)

      ! least contribution to furthest
      E2(level+1)%grid(ijf-1-Nf) = E2(level+1)%grid(ijf-1-Nf) + w3*E2(level)%grid(ij)
      E2(level+1)%grid(ijf-1+2*Nf) = E2(level+1)%grid(ijf-1+2*Nf) + w3*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+2-2*Nf) = E2(level+1)%grid(ijf+2-2*Nf) + w3*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+2*Nf) = E2(level+1)%grid(ijf+2+Nf) + w3*E2(level)%grid(ij)
    enddo

    ! top edge
    j = 1
    jf = 2*j-1
    do i=2,N-1
      if = 2*i-1

      ij = i+N*(j-1)
      ijf = if+Nf*(jf-1)

      ! largest contribution to nearest
      E1(level+1)%grid(ijf) = E1(level+1)%grid(ijf) + w1*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+1) = E1(level+1)%grid(ijf+1) + w1*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+Nf) = E1(level+1)%grid(ijf+Nf) + w1*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+Nf+1) = E1(level+1)%grid(ijf+Nf+1) + w1*E1(level)%grid(ij)

      ! lesser contribution to intermediate
      E1(level+1)%grid(ijf-1) = E1(level+1)%grid(ijf-1) + w2*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+Nf-1) = E1(level+1)%grid(ijf+Nf-1) + w2*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+Nf*(Nf-1)) = E1(level+1)%grid(ijf+Nf*(Nf-1)) + w2*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+2*Nf) = E1(level+1)%grid(ijf+2*Nf) + w2*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+1+Nf*(Nf-1)) = E1(level+1)%grid(ijf+1+Nf*(Nf-1)) + w2*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+1+2*Nf) = E1(level+1)%grid(ijf+1+2*Nf) + w2*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+2) = E1(level+1)%grid(ijf+2) + w2*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+2+Nf) = E1(level+1)%grid(ijf+2+Nf) + w2*E1(level)%grid(ij)

      ! least contribution to furthest
      E1(level+1)%grid(ijf-1+Nf*(Nf-1)) = E1(level+1)%grid(ijf-1+Nf*(Nf-1)) + w3*E1(level)%grid(ij)
      E1(level+1)%grid(ijf-1+2*Nf) = E1(level+1)%grid(ijf-1+2*Nf) + w3*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+2+Nf*(Nf-1)) = E1(level+1)%grid(ijf+2+Nf*(Nf-1)) + w3*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+2+2*Nf) = E1(level+1)%grid(ijf+2+2*Nf) + w3*E1(level)%grid(ij)

      ! largest contribution to nearest
      E2(level+1)%grid(ijf) = E2(level+1)%grid(ijf) + w1*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+1) = E2(level+1)%grid(ijf+1) + w1*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+Nf) = E2(level+1)%grid(ijf+Nf) + w1*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+Nf+1) = E2(level+1)%grid(ijf+Nf+1) + w1*E2(level)%grid(ij)

      ! lesser contribution to intermediate
      E2(level+1)%grid(ijf-1) = E2(level+1)%grid(ijf-1) + w2*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+Nf-1) = E2(level+1)%grid(ijf+Nf-1) + w2*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+Nf*(Nf-1)) = E2(level+1)%grid(ijf+Nf*(Nf-1)) + w2*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+2*Nf) = E2(level+1)%grid(ijf+2*Nf) + w2*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+1+Nf*(Nf-1)) = E2(level+1)%grid(ijf+1+Nf*(Nf-1)) + w2*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+1+2*Nf) = E2(level+1)%grid(ijf+1+2*Nf) + w2*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+2) = E2(level+1)%grid(ijf+2) + w2*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+2+Nf) = E2(level+1)%grid(ijf+2+Nf) + w2*E2(level)%grid(ij)

      ! least contribution to furthest
      E2(level+1)%grid(ijf-1+Nf*(Nf-1)) = E2(level+1)%grid(ijf-1+Nf*(Nf-1)) + w3*E2(level)%grid(ij)
      E2(level+1)%grid(ijf-1+2*Nf) = E2(level+1)%grid(ijf-1+2*Nf) + w3*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+2+Nf*(Nf-1)) = E2(level+1)%grid(ijf+2+Nf*(Nf-1)) + w3*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+2+2*Nf) = E2(level+1)%grid(ijf+2+2*Nf) + w3*E2(level)%grid(ij)
    enddo

    ! bottom edge
    j = N
    jf = 2*j-1
    do i=2,N-1
      if = 2*i-1

      ij = i+N*(j-1)
      ijf = if+Nf*(jf-1)

      ! largest contribution to nearest
      E1(level+1)%grid(ijf) = E1(level+1)%grid(ijf) + w1*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+1) = E1(level+1)%grid(ijf+1) + w1*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+Nf) = E1(level+1)%grid(ijf+Nf) + w1*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+Nf+1) = E1(level+1)%grid(ijf+Nf+1) + w1*E1(level)%grid(ij)

      ! lesser contribution to intermediate
      E1(level+1)%grid(ijf-1) = E1(level+1)%grid(ijf-1) + w2*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+Nf-1) = E1(level+1)%grid(ijf+Nf-1) + w2*E1(level)%grid(ij)
      E1(level+1)%grid(ijf-Nf) = E1(level+1)%grid(ijf-Nf) + w2*E1(level)%grid(ij)
      E1(level+1)%grid(if) = E1(level+1)%grid(if) + w2*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+1-Nf) = E1(level+1)%grid(ijf+1-Nf) + w2*E1(level)%grid(ij)
      E1(level+1)%grid(if+1) = E1(level+1)%grid(if+1) + w2*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+2) = E1(level+1)%grid(ijf+2) + w2*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+2+Nf) = E1(level+1)%grid(ijf+2+Nf) + w2*E1(level)%grid(ij)

      ! least contribution to furthest
      E1(level+1)%grid(ijf-1-Nf) = E1(level+1)%grid(ijf-1-Nf) + w3*E1(level)%grid(ij)
      E1(level+1)%grid(if-1) = E1(level+1)%grid(if-1) + w3*E1(level)%grid(ij)
      E1(level+1)%grid(ijf+2-Nf) = E1(level+1)%grid(ijf+2-Nf) + w3*E1(level)%grid(ij)
      E1(level+1)%grid(if+2) = E1(level+1)%grid(if+2) + w3*E1(level)%grid(ij)

      ! largest contribution to nearest
      E2(level+1)%grid(ijf) = E2(level+1)%grid(ijf) + w1*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+1) = E2(level+1)%grid(ijf+1) + w1*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+Nf) = E2(level+1)%grid(ijf+Nf) + w1*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+Nf+1) = E2(level+1)%grid(ijf+Nf+1) + w1*E2(level)%grid(ij)

      ! lesser contribution to intermediate
      E2(level+1)%grid(ijf-1) = E2(level+1)%grid(ijf-1) + w2*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+Nf-1) = E2(level+1)%grid(ijf+Nf-1) + w2*E2(level)%grid(ij)
      E2(level+1)%grid(ijf-Nf) = E2(level+1)%grid(ijf-Nf) + w2*E2(level)%grid(ij)
      E2(level+1)%grid(if) = E2(level+1)%grid(if) + w2*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+1-Nf) = E2(level+1)%grid(ijf+1-Nf) + w2*E2(level)%grid(ij)
      E2(level+1)%grid(if+1) = E2(level+1)%grid(if+1) + w2*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+2) = E2(level+1)%grid(ijf+2) + w2*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+2+Nf) = E2(level+1)%grid(ijf+2+Nf) + w2*E2(level)%grid(ij)

      ! least contribution to furthest
      E2(level+1)%grid(ijf-1-Nf) = E2(level+1)%grid(ijf-1-Nf) + w3*E2(level)%grid(ij)
      E2(level+1)%grid(if-1) = E2(level+1)%grid(if-1) + w3*E2(level)%grid(ij)
      E2(level+1)%grid(ijf+2-Nf) = E2(level+1)%grid(ijf+2-Nf) + w3*E2(level)%grid(ij)
      E2(level+1)%grid(if+2) = E2(level+1)%grid(if+2) + w3*E2(level)%grid(ij)
    enddo

    ! top left corner
    j = 1
    i = 1
    jf = 2*j-1
    if = 2*i-1
    ij = i+N*(j-1)
    ijf = if+Nf*(jf-1)

    ! largest contribution to nearest
    E1(level+1)%grid(ijf) = E1(level+1)%grid(ijf) + w1*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+1) = E1(level+1)%grid(ijf+1) + w1*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+Nf) = E1(level+1)%grid(ijf+Nf) + w1*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+Nf+1) = E1(level+1)%grid(ijf+Nf+1) + w1*E1(level)%grid(ij)

    ! lesser contribution to intermediate
    E1(level+1)%grid(ijf-1+Nf) = E1(level+1)%grid(ijf-1+Nf) + w2*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+2*Nf-1) = E1(level+1)%grid(ijf+2*Nf-1) + w2*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+Nf*(Nf-1)) = E1(level+1)%grid(ijf+Nf*(Nf-1)) + w2*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+2*Nf) = E1(level+1)%grid(ijf+2*Nf) + w2*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+1+Nf*(Nf-1)) = E1(level+1)%grid(ijf+1+Nf*(Nf-1)) + w2*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+1+2*Nf) = E1(level+1)%grid(ijf+1+2*Nf) + w2*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+2) = E1(level+1)%grid(ijf+2) + w2*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+2+Nf) = E1(level+1)%grid(ijf+2+Nf) + w2*E1(level)%grid(ij)

    ! least contribution to furthest
    E1(level+1)%grid(ijf-1+Nf*Nf) = E1(level+1)%grid(ijf-1+Nf*Nf) + w3*E1(level)%grid(ij)
    E1(level+1)%grid(ijf-1+3*Nf) = E1(level+1)%grid(ijf-1+3*Nf) + w3*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+2+Nf*(Nf-1)) = E1(level+1)%grid(ijf+2+Nf*(Nf-1)) + w3*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+2+2*Nf) = E1(level+1)%grid(ijf+2+2*Nf) + w3*E1(level)%grid(ij)

    ! largest contribution to nearest
    E2(level+1)%grid(ijf) = E2(level+1)%grid(ijf) + w1*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+1) = E2(level+1)%grid(ijf+1) + w1*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+Nf) = E2(level+1)%grid(ijf+Nf) + w1*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+Nf+1) = E2(level+1)%grid(ijf+Nf+1) + w1*E2(level)%grid(ij)

    ! lesser contribution to intermediate
    E2(level+1)%grid(ijf-1+Nf) = E2(level+1)%grid(ijf-1+Nf) + w2*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+2*Nf-1) = E2(level+1)%grid(ijf+2*Nf-1) + w2*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+Nf*(Nf-1)) = E2(level+1)%grid(ijf+Nf*(Nf-1)) + w2*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+2*Nf) = E2(level+1)%grid(ijf+2*Nf) + w2*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+1+Nf*(Nf-1)) = E2(level+1)%grid(ijf+1+Nf*(Nf-1)) + w2*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+1+2*Nf) = E2(level+1)%grid(ijf+1+2*Nf) + w2*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+2) = E2(level+1)%grid(ijf+2) + w2*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+2+Nf) = E2(level+1)%grid(ijf+2+Nf) + w2*E2(level)%grid(ij)

    ! least contribution to furthest
    E2(level+1)%grid(ijf-1+Nf*Nf) = E2(level+1)%grid(ijf-1+Nf*Nf) + w3*E2(level)%grid(ij)
    E2(level+1)%grid(ijf-1+3*Nf) = E2(level+1)%grid(ijf-1+3*Nf) + w3*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+2+Nf*(Nf-1)) = E2(level+1)%grid(ijf+2+Nf*(Nf-1)) + w3*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+2+2*Nf) = E2(level+1)%grid(ijf+2+2*Nf) + w3*E2(level)%grid(ij)


    ! top right corner
    j = 1
    i = N
    jf = 2*j-1
    if = 2*i-1
    ij = i+N*(j-1)
    ijf = if+Nf*(jf-1)

    ! largest contribution to nearest
    E1(level+1)%grid(ijf) = E1(level+1)%grid(ijf) + w1*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+1) = E1(level+1)%grid(ijf+1) + w1*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+Nf) = E1(level+1)%grid(ijf+Nf) + w1*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+Nf+1) = E1(level+1)%grid(ijf+Nf+1) + w1*E1(level)%grid(ij)

    ! lesser contribution to intermediate
    E1(level+1)%grid(ijf-1) = E1(level+1)%grid(ijf-1) + w2*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+Nf-1) = E1(level+1)%grid(ijf+Nf-1) + w2*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+Nf*(Nf-1)) = E1(level+1)%grid(ijf+Nf*(Nf-1)) + w2*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+2*Nf) = E1(level+1)%grid(ijf+2*Nf) + w2*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+1+Nf*(Nf-1)) = E1(level+1)%grid(ijf+1+Nf*(Nf-1)) + w2*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+1+2*Nf) = E1(level+1)%grid(ijf+1+2*Nf) + w2*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+2-Nf) = E1(level+1)%grid(ijf+2-Nf) + w2*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+2) = E1(level+1)%grid(ijf+2) + w2*E1(level)%grid(ij)

    ! least contribution to furthest
    E1(level+1)%grid(ijf-1+Nf*(Nf-1)) = E1(level+1)%grid(ijf-1+Nf*(Nf-1)) + w3*E1(level)%grid(ij)
    E1(level+1)%grid(ijf-1+2*Nf) = E1(level+1)%grid(ijf-1+2*Nf) + w3*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+2+Nf*(Nf-2)) = E1(level+1)%grid(ijf+2+Nf*(Nf-2)) + w3*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+2+Nf) = E1(level+1)%grid(ijf+2+Nf) + w3*E1(level)%grid(ij)

    ! largest contribution to nearest
    E2(level+1)%grid(ijf) = E2(level+1)%grid(ijf) + w1*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+1) = E2(level+1)%grid(ijf+1) + w1*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+Nf) = E2(level+1)%grid(ijf+Nf) + w1*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+Nf+1) = E2(level+1)%grid(ijf+Nf+1) + w1*E2(level)%grid(ij)

    ! lesser contribution to intermediate
    E2(level+1)%grid(ijf-1) = E2(level+1)%grid(ijf-1) + w2*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+Nf-1) = E2(level+1)%grid(ijf+Nf-1) + w2*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+Nf*(Nf-1)) = E2(level+1)%grid(ijf+Nf*(Nf-1)) + w2*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+2*Nf) = E2(level+1)%grid(ijf+2*Nf) + w2*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+1+Nf*(Nf-1)) = E2(level+1)%grid(ijf+1+Nf*(Nf-1)) + w2*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+1+2*Nf) = E2(level+1)%grid(ijf+1+2*Nf) + w2*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+2-Nf) = E2(level+1)%grid(ijf+2-Nf) + w2*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+2) = E2(level+1)%grid(ijf+2) + w2*E2(level)%grid(ij)

    ! least contribution to furthest
    E2(level+1)%grid(ijf-1+Nf*(Nf-1)) = E2(level+1)%grid(ijf-1+Nf*(Nf-1)) + w3*E2(level)%grid(ij)
    E2(level+1)%grid(ijf-1+2*Nf) = E2(level+1)%grid(ijf-1+2*Nf) + w3*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+2+Nf*(Nf-2)) = E2(level+1)%grid(ijf+2+Nf*(Nf-2)) + w3*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+2+Nf) = E2(level+1)%grid(ijf+2+Nf) + w3*E2(level)%grid(ij)


    ! bottom left corner
    j = N
    i = 1
    jf = 2*j-1
    if = 2*i-1
    ij = i+N*(j-1)
    ijf = if+Nf*(jf-1)

    ! largest contribution to nearest
    E1(level+1)%grid(ijf) = E1(level+1)%grid(ijf) + w1*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+1) = E1(level+1)%grid(ijf+1) + w1*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+Nf) = E1(level+1)%grid(ijf+Nf) + w1*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+Nf+1) = E1(level+1)%grid(ijf+Nf+1) + w1*E1(level)%grid(ij)

    ! lesser contribution to intermediate
    E1(level+1)%grid(ijf-1+Nf) = E1(level+1)%grid(ijf-1+Nf) + w2*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+2*Nf-1) = E1(level+1)%grid(ijf+2*Nf-1) + w2*E1(level)%grid(ij)
    E1(level+1)%grid(ijf-Nf) = E1(level+1)%grid(ijf-Nf) + w2*E1(level)%grid(ij)
    E1(level+1)%grid(if) = E1(level+1)%grid(if) + w2*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+1-Nf) = E1(level+1)%grid(ijf+1-Nf) + w2*E1(level)%grid(ij)
    E1(level+1)%grid(if+1) = E1(level+1)%grid(if+1) + w2*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+2) = E1(level+1)%grid(ijf+2) + w2*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+2+Nf) = E1(level+1)%grid(ijf+2+Nf) + w2*E1(level)%grid(ij)

    ! least contribution to furthest
    E1(level+1)%grid(ijf-1) = E1(level+1)%grid(ijf-1) + w3*E1(level)%grid(ij)
    E1(level+1)%grid(if-1+Nf) = E1(level+1)%grid(if-1+Nf) + w3*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+2-Nf) = E1(level+1)%grid(ijf+2-Nf) + w3*E1(level)%grid(ij)
    E1(level+1)%grid(if+2) = E1(level+1)%grid(if+2) + w3*E1(level)%grid(ij)

    ! largest contribution to nearest
    E2(level+1)%grid(ijf) = E2(level+1)%grid(ijf) + w1*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+1) = E2(level+1)%grid(ijf+1) + w1*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+Nf) = E2(level+1)%grid(ijf+Nf) + w1*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+Nf+1) = E2(level+1)%grid(ijf+Nf+1) + w1*E2(level)%grid(ij)

    ! lesser contribution to intermediate
    E2(level+1)%grid(ijf-1+Nf) = E2(level+1)%grid(ijf-1+Nf) + w2*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+2*Nf-1) = E2(level+1)%grid(ijf+2*Nf-1) + w2*E2(level)%grid(ij)
    E2(level+1)%grid(ijf-Nf) = E2(level+1)%grid(ijf-Nf) + w2*E2(level)%grid(ij)
    E2(level+1)%grid(if) = E2(level+1)%grid(if) + w2*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+1-Nf) = E2(level+1)%grid(ijf+1-Nf) + w2*E2(level)%grid(ij)
    E2(level+1)%grid(if+1) = E2(level+1)%grid(if+1) + w2*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+2) = E2(level+1)%grid(ijf+2) + w2*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+2+Nf) = E2(level+1)%grid(ijf+2+Nf) + w2*E2(level)%grid(ij)

    ! least contribution to furthest
    E2(level+1)%grid(ijf-1) = E2(level+1)%grid(ijf-1) + w3*E2(level)%grid(ij)
    E2(level+1)%grid(if-1+Nf) = E2(level+1)%grid(if-1+Nf) + w3*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+2-Nf) = E2(level+1)%grid(ijf+2-Nf) + w3*E2(level)%grid(ij)
    E2(level+1)%grid(if+2) = E2(level+1)%grid(if+2) + w3*E2(level)%grid(ij)


    ! bottom right corner
    j = N
    i = N
    jf = 2*j-1
    if = 2*i-1
    ij = i+N*(j-1)
    ijf = if+Nf*(jf-1)

    ! largest contribution to nearest
    E1(level+1)%grid(ijf) = E1(level+1)%grid(ijf) + w1*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+1) = E1(level+1)%grid(ijf+1) + w1*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+Nf) = E1(level+1)%grid(ijf+Nf) + w1*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+Nf+1) = E1(level+1)%grid(ijf+Nf+1) + w1*E1(level)%grid(ij)

    ! lesser contribution to intermediate
    E1(level+1)%grid(ijf-1) = E1(level+1)%grid(ijf-1) + w2*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+Nf-1) = E1(level+1)%grid(ijf+Nf-1) + w2*E1(level)%grid(ij)
    E1(level+1)%grid(ijf-Nf) = E1(level+1)%grid(ijf-Nf) + w2*E1(level)%grid(ij)
    E1(level+1)%grid(if) = E1(level+1)%grid(if) + w2*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+1-Nf) = E1(level+1)%grid(ijf+1-Nf) + w2*E1(level)%grid(ij)
    E1(level+1)%grid(if+1) = E1(level+1)%grid(if+1) + w2*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+2-Nf) = E1(level+1)%grid(ijf+2-Nf) + w2*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+2) = E1(level+1)%grid(ijf+2) + w2*E1(level)%grid(ij)

    ! least contribution to furthest
    E1(level+1)%grid(ijf-1-Nf) = E1(level+1)%grid(ijf-1-Nf) + w3*E1(level)%grid(ij)
    E1(level+1)%grid(if-1) = E1(level+1)%grid(if-1) + w3*E1(level)%grid(ij)
    E1(level+1)%grid(ijf+2-2*Nf) = E1(level+1)%grid(ijf+2-2*Nf) + w3*E1(level)%grid(ij)
    E1(level+1)%grid(if+2-Nf) = E1(level+1)%grid(if+2-Nf) + w3*E1(level)%grid(ij)

    ! largest contribution to nearest
    E2(level+1)%grid(ijf) = E2(level+1)%grid(ijf) + w1*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+1) = E2(level+1)%grid(ijf+1) + w1*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+Nf) = E2(level+1)%grid(ijf+Nf) + w1*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+Nf+1) = E2(level+1)%grid(ijf+Nf+1) + w1*E2(level)%grid(ij)

    ! lesser contribution to intermediate
    E2(level+1)%grid(ijf-1) = E2(level+1)%grid(ijf-1) + w2*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+Nf-1) = E2(level+1)%grid(ijf+Nf-1) + w2*E2(level)%grid(ij)
    E2(level+1)%grid(ijf-Nf) = E2(level+1)%grid(ijf-Nf) + w2*E2(level)%grid(ij)
    E2(level+1)%grid(if) = E2(level+1)%grid(if) + w2*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+1-Nf) = E2(level+1)%grid(ijf+1-Nf) + w2*E2(level)%grid(ij)
    E2(level+1)%grid(if+1) = E2(level+1)%grid(if+1) + w2*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+2-Nf) = E2(level+1)%grid(ijf+2-Nf) + w2*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+2) = E2(level+1)%grid(ijf+2) + w2*E2(level)%grid(ij)

    ! least contribution to furthest
    E2(level+1)%grid(ijf-1-Nf) = E2(level+1)%grid(ijf-1-Nf) + w3*E2(level)%grid(ij)
    E2(level+1)%grid(if-1) = E2(level+1)%grid(if-1) + w3*E2(level)%grid(ij)
    E2(level+1)%grid(ijf+2-2*Nf) = E2(level+1)%grid(ijf+2-2*Nf) + w3*E2(level)%grid(ij)
    E2(level+1)%grid(if+2-Nf) = E2(level+1)%grid(if+2-Nf) + w3*E2(level)%grid(ij)
  end subroutine prolongate


  ! TODO: remove, for testing only
  subroutine temp_output_data(b, fname)
    implicit none
    real(dp), dimension(:), allocatable, intent(in) :: b
    character(*), intent(in) :: fname
    integer :: i
    integer :: fid

    open(newunit=fid, file=fname, form='formatted')
    do i=1,ubound(b,1)
      write(fid, *) b(i)
    enddo
    close(fid)
  end subroutine temp_output_data
end module fd_solvers
