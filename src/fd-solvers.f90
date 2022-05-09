module multigrid
  use globals

  implicit none
  save

  !> @brief multigrid internal grid type
  !!
  !! This is just a wrapper around an array of reals. An array of these forms a
  !! multigrid
  type :: t_grid
    real(dp), pointer, contiguous :: grid(:)
  end type

  contains


  !> @brief Allocates storage for a multigrid with a given level
  !!
  !! @param mg[inout]  multigrid
  !! @param level[in]  maximum level
  subroutine multigrid_alloc(mg, level)
    implicit none
    type(t_grid), pointer, contiguous, intent(inout) :: mg(:)
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
    type(t_grid), pointer, contiguous, intent(inout) :: mg(:)
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
  subroutine solver_ufds2t2(Tout, c0, eps2)
    implicit none
    real(dp), intent(in) :: Tout(:)
    real(dp), intent(in) :: eps2
    real(dp), pointer, contiguous, intent(in) :: c0(:,:)

    integer :: N ! grid size
    integer :: level ! grid level
    real(dp) :: dx ! grid spacing
    real(dp) :: dt, dt0, dt1 ! timesteps
    real(dp) :: a0, a1, a2, b0, b1 ! time constants
    real(dp) :: t, tmax
    real(dp), dimension(2,2) :: A ! smoothing matrix
    real(dp) :: eps
    integer :: it, i, j ! iterators
    logical :: outflag
    character(len=48) :: msg ! logging message

    ! grid storage
    real(dp), pointer, contiguous :: phi(:), psi(:), g(:), b(:), &
                                     phi_prev(:), g_prev(:), work(:)
    type(t_grid), pointer, contiguous :: E1(:), E2(:), R1(:), R2(:)

    ! temp variables
    character(len=128) :: fname


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

    23 format(A, I0.5, A) ! output file format
    24 format(A, F7.3) ! output message

    ! allocate storage
    allocate(phi(N*N))
    allocate(psi(N*N))
    allocate(g(N*N))
    allocate(b(N*N))
    allocate(phi_prev(N*N))
    allocate(g_prev(N*N))
    allocate(work(N*N))

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

    ! output if required
    if (tout(it) < epsilon(tout(it))) then
      write(msg, 24) "Initial condition output at t=  0.000"
      call logger%info("solver_ufds2t2", msg)
      write(fname, 23) "out/phi-", it, ".dat"
      call temp_output_data(phi, fname)
      write(fname, 23) "out/psi-", it, ".dat"
      call temp_output_data(psi, fname)

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
    call compute_g(g, phi, dx, tau, N, work)
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

    ! conditionally output
    if (outflag) then
      write(msg, 24) "Output at t=", t
      call logger%info("solver_ufds2t2", msg)
      write(fname, 23) "out/phi-", it, ".dat"
      call temp_output_data(phi, fname)
      write(fname, 23) "out/psi-", it, ".dat"
      call temp_output_data(psi, fname)

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
        dt = dx
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
      call compute_g(g, phi, dx, tau, N, work)
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

      ! conditionally output
      if (outflag) then
        write(msg, 24) "Output at t=", t
        call logger%info("solver_ufds2t2", msg)
        write(fname, 23) "out/phi-", it, ".dat"
        call temp_output_data(phi, fname)
        write(fname, 23) "out/psi-", it, ".dat"
        call temp_output_data(psi, fname)

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
    real(dp), pointer, contiguous, intent(in) :: x(:)
    real(dp), pointer, contiguous, intent(out) :: y(:)
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
  subroutine compute_g(g, phi, dx, tau, n, work)
    implicit none
    real(dp), pointer, contiguous, intent(in) :: phi(:)
    real(dp), pointer, contiguous, intent(out) :: work(:)
    real(dp), pointer, contiguous, intent(out) :: g(:)
    real(dp), intent(in) :: dx, tau
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
    type(t_grid), pointer, contiguous, intent(inout) :: E1(:), E2(:), R1(:), R2(:)
    real(dp), intent(in) :: eps2, dx
    integer, intent(in) :: N, level

    E1(level)%grid = 0.0_dp
    E2(level)%grid = 0.0_dp

    ! start smooth
    call smooth(A, E1, E2, R1, R2, eps2, N, dx, level)

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
      call smooth(A, E1, E2, R1, R2, eps2, N, dx, level)
      call smooth(A, E1, E2, R1, R2, eps2, N, dx, level)
    endif

    ! start smooth
    call smooth(A, E1, E2, R1, R2, eps2, N, dx, level)
  end subroutine vcycle


  !> @brief performs a single red/black smooth
  !!
  !! @param A      2x2 smoothing matrix (inverse solution matrix for a 1x1 system)
  !! @param E1     error multigrid for the first variable
  !! @param E2     error multigrid for the second variable
  !! @param R1     residual multigrid for the first variable
  !! @param R2     residual multigrid for the second variable
  !! @param eps2   PDE parameter
  !! @param N      grid size
  !! @param dx     grid spacing
  !! @param level  grid level
  subroutine smooth(A, E1, E2, R1, R2, eps2, N, dx, level)
    implicit none
    real(dp), dimension(2,2), intent(in) :: A
    type(t_grid), pointer, contiguous, intent(inout) :: E1(:), E2(:), R1(:), R2(:)
    real(dp), intent(in) :: eps2, dx
    integer, intent(in) :: N, level

    real(dp), dimension(2) :: rhs
    real(dp) :: dx2_
    integer :: i, j, ij, shift
    real(dp), pointer, contiguous :: e1g(:), e2g(:), r1g(:), r2g(:)

    dx2_ = 1.0_dp / (dx*dx)

    e1g => E1(level)%grid
    e2g => E2(level)%grid
    r1g => R1(level)%grid
    r2g => R2(level)%grid

    ! ================ !
    ! SMOOTH RED NODES !
    ! ================ !
    ! interior
    do j=2,n-1
      shift = mod(j,2)
      do i=2+shift,n-1+shift,2
        ij = i+n*(j-1)

        ! compute RHS
        rhs(1) = r1g(ij) + dx2_ * (e2g(ij+1)+e2g(ij-1)+e2g(ij+n)+e2g(ij-n))
        rhs(2) = r2g(ij) - eps2*dx2_ * (e1g(ij+1)+e1g(ij-1)+e1g(ij+n)+e1g(ij-n))

        ! solve for new errors
        e1g(ij) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
        e2g(ij) = A(2,1)*rhs(1) + A(2,2)*rhs(2)
      enddo
    enddo

    ! left/right
    do j=2,n-1,2
      ij = 1+n*j
      rhs(1) = r1g(ij) + dx2_ * (e2g(ij+1)+e2g(ij+n-1)+e2g(ij+n)+e2g(ij-n))
      rhs(2) = r2g(ij) - eps2*dx2_ * (e1g(ij+1)+e1g(ij+n-1)+e1g(ij+n)+e1g(ij-n))
      e1g(ij) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
      e2g(ij) = A(2,1)*rhs(1) + A(2,2)*rhs(2)

      ij = n+n*(j-1)
      rhs(1) = r1g(ij) + dx2_ * (e2g(ij-n+1)+e2g(ij-1)+e2g(ij+n)+e2g(ij-n))
      rhs(2) = r2g(ij) - eps2*dx2_ * (e1g(ij-n+1)+e1g(ij-1)+e1g(ij+n)+e1g(ij-n))
      e1g(ij) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
      e2g(ij) = A(2,1)*rhs(1) + A(2,2)*rhs(2)
    enddo

    ! top/bottom
    do i=2,n-1,2
      ij = i+1+n*(1-1)
      rhs(1) = r1g(ij) + dx2_ * (e2g(ij+1)+e2g(ij-1)+e2g(ij+n)+e2g(ij+n*(n-1)))
      rhs(2) = r2g(ij) - eps2*dx2_ * (e1g(ij+1)+e1g(ij-1)+e1g(ij+n)+e1g(ij+n*(n-1)))
      e1g(ij) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
      e2g(ij) = A(2,1)*rhs(1) + A(2,2)*rhs(2)

      ij = i+n*(n-1)
      rhs(1) = r1g(ij) + dx2_ * (e2g(ij+1)+e2g(ij-1)+e2g(ij-n*(n-1))+e2g(ij-n))
      rhs(2) = r2g(ij) - eps2*dx2_ * (e1g(ij+1)+e1g(ij-1)+e1g(ij-n*(n-1))+e1g(ij-n))
      e1g(ij) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
      e2g(ij) = A(2,1)*rhs(1) + A(2,2)*rhs(2)
    enddo

    ! upper left corner
    ij = 1+n*(1-1)
    rhs(1) = r1g(ij) + dx2_ * (e2g(ij+1)+e2g(ij+n-1)+e2g(ij+n)+e2g(ij+n*(n-1)))
    rhs(2) = r2g(ij) - eps2*dx2_ * (e1g(ij+1)+e1g(ij+n-1)+e1g(ij+n)+e1g(ij+n*(n-1)))
    e1g(ij) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
    e2g(ij) = A(2,1)*rhs(1) + A(2,2)*rhs(2)

    ! lower right corner
    ij = n+n*(n-1)
    rhs(1) = r1g(ij) + dx2_ * (e2g(ij-n+1)+e2g(ij-1)+e2g(ij-n*(n-1))+e2g(ij-n))
    rhs(2) = r2g(ij) - eps2*dx2_ * (e1g(ij-n+1)+e1g(ij-1)+e1g(ij-n*(n-1))+e1g(ij-n))
    e1g(ij) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
    e2g(ij) = A(2,1)*rhs(1) + A(2,2)*rhs(2)


    ! ================== !
    ! SMOOTH BLACK NODES !
    ! ================== !
    ! interior
    do j=2,n-1
      shift = mod(j-1,2)
      do i=2+shift,n-1+shift,2
        ij = i+n*(j-1)

        ! compute RHS
        rhs(1) = r1g(ij) + dx2_ * (e2g(ij+1)+e2g(ij-1)+e2g(ij+n)+e2g(ij-n))
        rhs(2) = r2g(ij) - eps2*dx2_ * (e1g(ij+1)+e1g(ij-1)+e1g(ij+n)+e1g(ij-n))

        ! solve for new errors
        e1g(ij) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
        e2g(ij) = A(2,1)*rhs(1) + A(2,2)*rhs(2)
      enddo
    enddo

    ! left/right
    do j=2,n-1,2
      ij = 1+n*(j-1)
      rhs(1) = r1g(ij) + dx2_ * (e2g(ij+1)+e2g(ij+n-1)+e2g(ij+n)+e2g(ij-n))
      rhs(2) = r2g(ij) - eps2*dx2_ * (e1g(ij+1)+e1g(ij+n-1)+e1g(ij+n)+e1g(ij-n))
      e1g(ij) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
      e2g(ij) = A(2,1)*rhs(1) + A(2,2)*rhs(2)

      ij = n+n*j
      rhs(1) = r1g(ij) + dx2_ * (e2g(ij-n+1)+e2g(ij-1)+e2g(ij+n)+e2g(ij-n))
      rhs(2) = r2g(ij) - eps2*dx2_ * (e1g(ij-n+1)+e1g(ij-1)+e1g(ij+n)+e1g(ij-n))
      e1g(ij) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
      e2g(ij) = A(2,1)*rhs(1) + A(2,2)*rhs(2)
    enddo

    ! top/bottom
    do i=2,n-1,2
      ij = i+n*(1-1)
      rhs(1) = r1g(ij) + dx2_ * (e2g(ij+1)+e2g(ij-1)+e2g(ij+n)+e2g(ij+n*(n-1)))
      rhs(2) = r2g(ij) - eps2*dx2_ * (e1g(ij+1)+e1g(ij-1)+e1g(ij+n)+e1g(ij+n*(n-1)))
      e1g(ij) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
      e2g(ij) = A(2,1)*rhs(1) + A(2,2)*rhs(2)

      ij = i+1+n*(n-1)
      rhs(1) = r1g(ij) + dx2_ * (e2g(ij+1)+e2g(ij-1)+e2g(ij-n*(n-1))+e2g(ij-n))
      rhs(2) = r2g(ij) - eps2*dx2_ * (e1g(ij+1)+e1g(ij-1)+e1g(ij-n*(n-1))+e1g(ij-n))
      e1g(ij) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
      e2g(ij) = A(2,1)*rhs(1) + A(2,2)*rhs(2)
    enddo

    ! upper right
    ij = 1+n*(n-1)
    rhs(1) = r1g(ij) + dx2_ * (e2g(ij+1)+e2g(ij+n-1)+e2g(ij-n*(n-1))+e2g(ij-n))
    rhs(2) = r2g(ij) - eps2*dx2_ * (e1g(ij+1)+e1g(ij+n-1)+e1g(ij-n*(n-1))+e1g(ij-n))
    e1g(ij) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
    e2g(ij) = A(2,1)*rhs(1) + A(2,2)*rhs(2)

    ! lower left
    ij = n+n*(1-1)
    rhs(1) = r1g(ij) + dx2_ * (e2g(ij-n+1)+e2g(ij-1)+e2g(ij+n)+e2g(ij+n*(n-1)))
    rhs(2) = r2g(ij) - eps2*dx2_ * (e1g(ij-n+1)+e1g(ij-1)+e1g(ij+n)+e1g(ij+n*(n-1)))
    e1g(ij) = A(1,1)*rhs(1) + A(1,2)*rhs(2)
    e2g(ij) = A(2,1)*rhs(1) + A(2,2)*rhs(2)
  end subroutine smooth


  !> @brief restricts (ie coarsens) to level-1 from level
  !!
  !! @param R1     residual multigrid for the first variable
  !! @param R2     residual multigrid for the second variable
  !! @param N      grid size
  !! @param level  grid level
  subroutine restrict(R1, R2, N, level)
    implicit none
    type(t_grid), pointer, contiguous, intent(inout) :: R1(:), R2(:)
    integer, intent(in) :: N, level

    integer :: i, j, ij, im, jm, ijm

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
    type(t_grid), pointer, contiguous, intent(inout) :: E1(:), E2(:)
    integer, intent(in) :: N, level

    real(dp), pointer, contiguous :: e1c(:), e2c(:), e1f(:), e2f(:)
    integer :: Nf
    integer :: i, j, ij, if, jf, ijf
    real(dp), parameter :: w1 = 0.5625_dp
    real(dp), parameter :: w2 = 0.1875_dp
    real(dp), parameter :: w3 = 0.0625_dp

    Nf = 2*N
    e1c => E1(level)%grid
    e2c => E2(level)%grid
    e1f => E1(level+1)%grid
    e2f => E2(level+1)%grid

    ! interior
    do j=2,N-1
      jf = 2*j-1
      do i=2,N-1
        if = 2*i-1

        ij = i+N*(j-1)
        ijf = if+Nf*(jf-1)

        ! largest contribution to nearest
        e1f(ijf) = e1f(ijf) + w1*e1c(ij)
        e1f(ijf+1) = e1f(ijf+1) + w1*e1c(ij)
        e1f(ijf+Nf) = e1f(ijf+Nf) + w1*e1c(ij)
        e1f(ijf+Nf+1) = e1f(ijf+Nf+1) + w1*e1c(ij)

        ! lesser contribution to intermediate
        e1f(ijf-1) = e1f(ijf-1) + w2*e1c(ij)
        e1f(ijf+Nf-1) = e1f(ijf+Nf-1) + w2*e1c(ij)
        e1f(ijf-Nf) = e1f(ijf-Nf) + w2*e1c(ij)
        e1f(ijf+2*Nf) = e1f(ijf+2*Nf) + w2*e1c(ij)
        e1f(ijf+1-Nf) = e1f(ijf+1-Nf) + w2*e1c(ij)
        e1f(ijf+1+2*Nf) = e1f(ijf+1+2*Nf) + w2*e1c(ij)
        e1f(ijf+2) = e1f(ijf+2) + w2*e1c(ij)
        e1f(ijf+2+Nf) = e1f(ijf+2+Nf) + w2*e1c(ij)

        ! least contribution to furthest
        e1f(ijf-1-Nf) = e1f(ijf-1-Nf) + w3*e1c(ij)
        e1f(ijf-1+2*Nf) = e1f(ijf-1+2*Nf) + w3*e1c(ij)
        e1f(ijf+2-Nf) = e1f(ijf+2-Nf) + w3*e1c(ij)
        e1f(ijf+2+2*Nf) = e1f(ijf+2+2*Nf) + w3*e1c(ij)

        ! largest contribution to nearest
        e2f(ijf) = e2f(ijf) + w1*e2c(ij)
        e2f(ijf+1) = e2f(ijf+1) + w1*e2c(ij)
        e2f(ijf+Nf) = e2f(ijf+Nf) + w1*e2c(ij)
        e2f(ijf+Nf+1) = e2f(ijf+Nf+1) + w1*e2c(ij)

        ! lesser contribution to intermediate
        e2f(ijf-1) = e2f(ijf-1) + w2*e2c(ij)
        e2f(ijf+Nf-1) = e2f(ijf+Nf-1) + w2*e2c(ij)
        e2f(ijf-Nf) = e2f(ijf-Nf) + w2*e2c(ij)
        e2f(ijf+2*Nf) = e2f(ijf+2*Nf) + w2*e2c(ij)
        e2f(ijf+1-Nf) = e2f(ijf+1-Nf) + w2*e2c(ij)
        e2f(ijf+1+2*Nf) = e2f(ijf+1+2*Nf) + w2*e2c(ij)
        e2f(ijf+2) = e2f(ijf+2) + w2*e2c(ij)
        e2f(ijf+2+Nf) = e2f(ijf+2+Nf) + w2*e2c(ij)

        ! least contribution to furthest
        e2f(ijf-1-Nf) = e2f(ijf-1-Nf) + w3*e2c(ij)
        e2f(ijf-1+2*Nf) = e2f(ijf-1+2*Nf) + w3*e2c(ij)
        e2f(ijf+2-Nf) = e2f(ijf+2-Nf) + w3*e2c(ij)
        e2f(ijf+2+2*Nf) = e2f(ijf+2+2*Nf) + w3*e2c(ij)
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
      e1f(ijf) = e1f(ijf) + w1*e1c(ij)
      e1f(ijf+1) = e1f(ijf+1) + w1*e1c(ij)
      e1f(ijf+Nf) = e1f(ijf+Nf) + w1*e1c(ij)
      e1f(ijf+Nf+1) = e1f(ijf+Nf+1) + w1*e1c(ij)

      ! lesser contribution to intermediate
      e1f(ijf-1+Nf) = e1f(ijf-1+Nf) + w2*e1c(ij)
      e1f(ijf+2*Nf-1) = e1f(ijf+2*Nf-1) + w2*e1c(ij)
      e1f(ijf-Nf) = e1f(ijf-Nf) + w2*e1c(ij)
      e1f(ijf+2*Nf) = e1f(ijf+2*Nf) + w2*e1c(ij)
      e1f(ijf+1-Nf) = e1f(ijf+1-Nf) + w2*e1c(ij)
      e1f(ijf+1+2*Nf) = e1f(ijf+1+2*Nf) + w2*e1c(ij)
      e1f(ijf+2) = e1f(ijf+2) + w2*e1c(ij)
      e1f(ijf+2+Nf) = e1f(ijf+2+Nf) + w2*e1c(ij)

      ! least contribution to furthest
      e1f(ijf-1) = e1f(ijf-1) + w3*e1c(ij)
      e1f(ijf-1+3*Nf) = e1f(ijf-1+3*Nf) + w3*e1c(ij)
      e1f(ijf+2-Nf) = e1f(ijf+2-Nf) + w3*e1c(ij)
      e1f(ijf+2+2*Nf) = e1f(ijf+2+2*Nf) + w3*e1c(ij)

      ! largest contribution to nearest
      e2f(ijf) = e2f(ijf) + w1*e2c(ij)
      e2f(ijf+1) = e2f(ijf+1) + w1*e2c(ij)
      e2f(ijf+Nf) = e2f(ijf+Nf) + w1*e2c(ij)
      e2f(ijf+Nf+1) = e2f(ijf+Nf+1) + w1*e2c(ij)

      ! lesser contribution to intermediate
      e2f(ijf-1+Nf) = e2f(ijf-1+Nf) + w2*e2c(ij)
      e2f(ijf+2*Nf-1) = e2f(ijf+2*Nf-1) + w2*e2c(ij)
      e2f(ijf-Nf) = e2f(ijf-Nf) + w2*e2c(ij)
      e2f(ijf+2*Nf) = e2f(ijf+2*Nf) + w2*e2c(ij)
      e2f(ijf+1-Nf) = e2f(ijf+1-Nf) + w2*e2c(ij)
      e2f(ijf+1+2*Nf) = e2f(ijf+1+2*Nf) + w2*e2c(ij)
      e2f(ijf+2) = e2f(ijf+2) + w2*e2c(ij)
      e2f(ijf+2+Nf) = e2f(ijf+2+Nf) + w2*e2c(ij)

      ! least contribution to furthest
      e2f(ijf-1) = e2f(ijf-1) + w3*e2c(ij)
      e2f(ijf-1+3*Nf) = e2f(ijf-1+3*Nf) + w3*e2c(ij)
      e2f(ijf+2-Nf) = e2f(ijf+2-Nf) + w3*e2c(ij)
      e2f(ijf+2+2*Nf) = e2f(ijf+2+2*Nf) + w3*e2c(ij)
    enddo

    ! right edge
    i = N
    if = 2*i-1
    do j=2,N-1
      jf = 2*j-1

      ij = i+N*(j-1)
      ijf = if+Nf*(jf-1)

      ! largest contribution to nearest
      e1f(ijf) = e1f(ijf) + w1*e1c(ij)
      e1f(ijf+1) = e1f(ijf+1) + w1*e1c(ij)
      e1f(ijf+Nf) = e1f(ijf+Nf) + w1*e1c(ij)
      e1f(ijf+Nf+1) = e1f(ijf+Nf+1) + w1*e1c(ij)

      ! lesser contribution to intermediate
      e1f(ijf-1) = e1f(ijf-1) + w2*e1c(ij)
      e1f(ijf+Nf-1) = e1f(ijf+Nf-1) + w2*e1c(ij)
      e1f(ijf-Nf) = e1f(ijf-Nf) + w2*e1c(ij)
      e1f(ijf+2*Nf) = e1f(ijf+2*Nf) + w2*e1c(ij)
      e1f(ijf+1-Nf) = e1f(ijf+1-Nf) + w2*e1c(ij)
      e1f(ijf+1+2*Nf) = e1f(ijf+1+2*Nf) + w2*e1c(ij)
      e1f(ijf+2-Nf) = e1f(ijf+2-Nf) + w2*e1c(ij)
      e1f(ijf+2) = e1f(ijf+2) + w2*e1c(ij)

      ! least contribution to furthest
      e1f(ijf-1-Nf) = e1f(ijf-1-Nf) + w3*e1c(ij)
      e1f(ijf-1+2*Nf) = e1f(ijf-1+2*Nf) + w3*e1c(ij)
      e1f(ijf+2-2*Nf) = e1f(ijf+2-2*Nf) + w3*e1c(ij)
      e1f(ijf+2*Nf) = e1f(ijf+2+Nf) + w3*e1c(ij)

      ! largest contribution to nearest
      e2f(ijf) = e2f(ijf) + w1*e2c(ij)
      e2f(ijf+1) = e2f(ijf+1) + w1*e2c(ij)
      e2f(ijf+Nf) = e2f(ijf+Nf) + w1*e2c(ij)
      e2f(ijf+Nf+1) = e2f(ijf+Nf+1) + w1*e2c(ij)

      ! lesser contribution to intermediate
      e2f(ijf-1) = e2f(ijf-1) + w2*e2c(ij)
      e2f(ijf+Nf-1) = e2f(ijf+Nf-1) + w2*e2c(ij)
      e2f(ijf-Nf) = e2f(ijf-Nf) + w2*e2c(ij)
      e2f(ijf+2*Nf) = e2f(ijf+2*Nf) + w2*e2c(ij)
      e2f(ijf+1-Nf) = e2f(ijf+1-Nf) + w2*e2c(ij)
      e2f(ijf+1+2*Nf) = e2f(ijf+1+2*Nf) + w2*e2c(ij)
      e2f(ijf+2-Nf) = e2f(ijf+2-Nf) + w2*e2c(ij)
      e2f(ijf+2) = e2f(ijf+2) + w2*e2c(ij)

      ! least contribution to furthest
      e2f(ijf-1-Nf) = e2f(ijf-1-Nf) + w3*e2c(ij)
      e2f(ijf-1+2*Nf) = e2f(ijf-1+2*Nf) + w3*e2c(ij)
      e2f(ijf+2-2*Nf) = e2f(ijf+2-2*Nf) + w3*e2c(ij)
      e2f(ijf+2*Nf) = e2f(ijf+2+Nf) + w3*e2c(ij)
    enddo

    ! top edge
    j = 1
    jf = 2*j-1
    do i=2,N-1
      if = 2*i-1

      ij = i+N*(j-1)
      ijf = if+Nf*(jf-1)

      ! largest contribution to nearest
      e1f(ijf) = e1f(ijf) + w1*e1c(ij)
      e1f(ijf+1) = e1f(ijf+1) + w1*e1c(ij)
      e1f(ijf+Nf) = e1f(ijf+Nf) + w1*e1c(ij)
      e1f(ijf+Nf+1) = e1f(ijf+Nf+1) + w1*e1c(ij)

      ! lesser contribution to intermediate
      e1f(ijf-1) = e1f(ijf-1) + w2*e1c(ij)
      e1f(ijf+Nf-1) = e1f(ijf+Nf-1) + w2*e1c(ij)
      e1f(ijf+Nf*(Nf-1)) = e1f(ijf+Nf*(Nf-1)) + w2*e1c(ij)
      e1f(ijf+2*Nf) = e1f(ijf+2*Nf) + w2*e1c(ij)
      e1f(ijf+1+Nf*(Nf-1)) = e1f(ijf+1+Nf*(Nf-1)) + w2*e1c(ij)
      e1f(ijf+1+2*Nf) = e1f(ijf+1+2*Nf) + w2*e1c(ij)
      e1f(ijf+2) = e1f(ijf+2) + w2*e1c(ij)
      e1f(ijf+2+Nf) = e1f(ijf+2+Nf) + w2*e1c(ij)

      ! least contribution to furthest
      e1f(ijf-1+Nf*(Nf-1)) = e1f(ijf-1+Nf*(Nf-1)) + w3*e1c(ij)
      e1f(ijf-1+2*Nf) = e1f(ijf-1+2*Nf) + w3*e1c(ij)
      e1f(ijf+2+Nf*(Nf-1)) = e1f(ijf+2+Nf*(Nf-1)) + w3*e1c(ij)
      e1f(ijf+2+2*Nf) = e1f(ijf+2+2*Nf) + w3*e1c(ij)

      ! largest contribution to nearest
      e2f(ijf) = e2f(ijf) + w1*e2c(ij)
      e2f(ijf+1) = e2f(ijf+1) + w1*e2c(ij)
      e2f(ijf+Nf) = e2f(ijf+Nf) + w1*e2c(ij)
      e2f(ijf+Nf+1) = e2f(ijf+Nf+1) + w1*e2c(ij)

      ! lesser contribution to intermediate
      e2f(ijf-1) = e2f(ijf-1) + w2*e2c(ij)
      e2f(ijf+Nf-1) = e2f(ijf+Nf-1) + w2*e2c(ij)
      e2f(ijf+Nf*(Nf-1)) = e2f(ijf+Nf*(Nf-1)) + w2*e2c(ij)
      e2f(ijf+2*Nf) = e2f(ijf+2*Nf) + w2*e2c(ij)
      e2f(ijf+1+Nf*(Nf-1)) = e2f(ijf+1+Nf*(Nf-1)) + w2*e2c(ij)
      e2f(ijf+1+2*Nf) = e2f(ijf+1+2*Nf) + w2*e2c(ij)
      e2f(ijf+2) = e2f(ijf+2) + w2*e2c(ij)
      e2f(ijf+2+Nf) = e2f(ijf+2+Nf) + w2*e2c(ij)

      ! least contribution to furthest
      e2f(ijf-1+Nf*(Nf-1)) = e2f(ijf-1+Nf*(Nf-1)) + w3*e2c(ij)
      e2f(ijf-1+2*Nf) = e2f(ijf-1+2*Nf) + w3*e2c(ij)
      e2f(ijf+2+Nf*(Nf-1)) = e2f(ijf+2+Nf*(Nf-1)) + w3*e2c(ij)
      e2f(ijf+2+2*Nf) = e2f(ijf+2+2*Nf) + w3*e2c(ij)
    enddo

    ! bottom edge
    j = N
    jf = 2*j-1
    do i=2,N-1
      if = 2*i-1

      ij = i+N*(j-1)
      ijf = if+Nf*(jf-1)

      ! largest contribution to nearest
      e1f(ijf) = e1f(ijf) + w1*e1c(ij)
      e1f(ijf+1) = e1f(ijf+1) + w1*e1c(ij)
      e1f(ijf+Nf) = e1f(ijf+Nf) + w1*e1c(ij)
      e1f(ijf+Nf+1) = e1f(ijf+Nf+1) + w1*e1c(ij)

      ! lesser contribution to intermediate
      e1f(ijf-1) = e1f(ijf-1) + w2*e1c(ij)
      e1f(ijf+Nf-1) = e1f(ijf+Nf-1) + w2*e1c(ij)
      e1f(ijf-Nf) = e1f(ijf-Nf) + w2*e1c(ij)
      e1f(if) = e1f(if) + w2*e1c(ij)
      e1f(ijf+1-Nf) = e1f(ijf+1-Nf) + w2*e1c(ij)
      e1f(if+1) = e1f(if+1) + w2*e1c(ij)
      e1f(ijf+2) = e1f(ijf+2) + w2*e1c(ij)
      e1f(ijf+2+Nf) = e1f(ijf+2+Nf) + w2*e1c(ij)

      ! least contribution to furthest
      e1f(ijf-1-Nf) = e1f(ijf-1-Nf) + w3*e1c(ij)
      e1f(if-1) = e1f(if-1) + w3*e1c(ij)
      e1f(ijf+2-Nf) = e1f(ijf+2-Nf) + w3*e1c(ij)
      e1f(if+2) = e1f(if+2) + w3*e1c(ij)

      ! largest contribution to nearest
      e2f(ijf) = e2f(ijf) + w1*e2c(ij)
      e2f(ijf+1) = e2f(ijf+1) + w1*e2c(ij)
      e2f(ijf+Nf) = e2f(ijf+Nf) + w1*e2c(ij)
      e2f(ijf+Nf+1) = e2f(ijf+Nf+1) + w1*e2c(ij)

      ! lesser contribution to intermediate
      e2f(ijf-1) = e2f(ijf-1) + w2*e2c(ij)
      e2f(ijf+Nf-1) = e2f(ijf+Nf-1) + w2*e2c(ij)
      e2f(ijf-Nf) = e2f(ijf-Nf) + w2*e2c(ij)
      e2f(if) = e2f(if) + w2*e2c(ij)
      e2f(ijf+1-Nf) = e2f(ijf+1-Nf) + w2*e2c(ij)
      e2f(if+1) = e2f(if+1) + w2*e2c(ij)
      e2f(ijf+2) = e2f(ijf+2) + w2*e2c(ij)
      e2f(ijf+2+Nf) = e2f(ijf+2+Nf) + w2*e2c(ij)

      ! least contribution to furthest
      e2f(ijf-1-Nf) = e2f(ijf-1-Nf) + w3*e2c(ij)
      e2f(if-1) = e2f(if-1) + w3*e2c(ij)
      e2f(ijf+2-Nf) = e2f(ijf+2-Nf) + w3*e2c(ij)
      e2f(if+2) = e2f(if+2) + w3*e2c(ij)
    enddo

    ! top left corner
    j = 1
    i = 1
    jf = 2*j-1
    if = 2*i-1
    ij = i+N*(j-1)
    ijf = if+Nf*(jf-1)

    ! largest contribution to nearest
    e1f(ijf) = e1f(ijf) + w1*e1c(ij)
    e1f(ijf+1) = e1f(ijf+1) + w1*e1c(ij)
    e1f(ijf+Nf) = e1f(ijf+Nf) + w1*e1c(ij)
    e1f(ijf+Nf+1) = e1f(ijf+Nf+1) + w1*e1c(ij)

    ! lesser contribution to intermediate
    e1f(ijf-1+Nf) = e1f(ijf-1+Nf) + w2*e1c(ij)
    e1f(ijf+2*Nf-1) = e1f(ijf+2*Nf-1) + w2*e1c(ij)
    e1f(ijf+Nf*(Nf-1)) = e1f(ijf+Nf*(Nf-1)) + w2*e1c(ij)
    e1f(ijf+2*Nf) = e1f(ijf+2*Nf) + w2*e1c(ij)
    e1f(ijf+1+Nf*(Nf-1)) = e1f(ijf+1+Nf*(Nf-1)) + w2*e1c(ij)
    e1f(ijf+1+2*Nf) = e1f(ijf+1+2*Nf) + w2*e1c(ij)
    e1f(ijf+2) = e1f(ijf+2) + w2*e1c(ij)
    e1f(ijf+2+Nf) = e1f(ijf+2+Nf) + w2*e1c(ij)

    ! least contribution to furthest
    e1f(ijf-1+Nf*Nf) = e1f(ijf-1+Nf*Nf) + w3*e1c(ij)
    e1f(ijf-1+3*Nf) = e1f(ijf-1+3*Nf) + w3*e1c(ij)
    e1f(ijf+2+Nf*(Nf-1)) = e1f(ijf+2+Nf*(Nf-1)) + w3*e1c(ij)
    e1f(ijf+2+2*Nf) = e1f(ijf+2+2*Nf) + w3*e1c(ij)

    ! largest contribution to nearest
    e2f(ijf) = e2f(ijf) + w1*e2c(ij)
    e2f(ijf+1) = e2f(ijf+1) + w1*e2c(ij)
    e2f(ijf+Nf) = e2f(ijf+Nf) + w1*e2c(ij)
    e2f(ijf+Nf+1) = e2f(ijf+Nf+1) + w1*e2c(ij)

    ! lesser contribution to intermediate
    e2f(ijf-1+Nf) = e2f(ijf-1+Nf) + w2*e2c(ij)
    e2f(ijf+2*Nf-1) = e2f(ijf+2*Nf-1) + w2*e2c(ij)
    e2f(ijf+Nf*(Nf-1)) = e2f(ijf+Nf*(Nf-1)) + w2*e2c(ij)
    e2f(ijf+2*Nf) = e2f(ijf+2*Nf) + w2*e2c(ij)
    e2f(ijf+1+Nf*(Nf-1)) = e2f(ijf+1+Nf*(Nf-1)) + w2*e2c(ij)
    e2f(ijf+1+2*Nf) = e2f(ijf+1+2*Nf) + w2*e2c(ij)
    e2f(ijf+2) = e2f(ijf+2) + w2*e2c(ij)
    e2f(ijf+2+Nf) = e2f(ijf+2+Nf) + w2*e2c(ij)

    ! least contribution to furthest
    e2f(ijf-1+Nf*Nf) = e2f(ijf-1+Nf*Nf) + w3*e2c(ij)
    e2f(ijf-1+3*Nf) = e2f(ijf-1+3*Nf) + w3*e2c(ij)
    e2f(ijf+2+Nf*(Nf-1)) = e2f(ijf+2+Nf*(Nf-1)) + w3*e2c(ij)
    e2f(ijf+2+2*Nf) = e2f(ijf+2+2*Nf) + w3*e2c(ij)


    ! top right corner
    j = 1
    i = N
    jf = 2*j-1
    if = 2*i-1
    ij = i+N*(j-1)
    ijf = if+Nf*(jf-1)

    ! largest contribution to nearest
    e1f(ijf) = e1f(ijf) + w1*e1c(ij)
    e1f(ijf+1) = e1f(ijf+1) + w1*e1c(ij)
    e1f(ijf+Nf) = e1f(ijf+Nf) + w1*e1c(ij)
    e1f(ijf+Nf+1) = e1f(ijf+Nf+1) + w1*e1c(ij)

    ! lesser contribution to intermediate
    e1f(ijf-1) = e1f(ijf-1) + w2*e1c(ij)
    e1f(ijf+Nf-1) = e1f(ijf+Nf-1) + w2*e1c(ij)
    e1f(ijf+Nf*(Nf-1)) = e1f(ijf+Nf*(Nf-1)) + w2*e1c(ij)
    e1f(ijf+2*Nf) = e1f(ijf+2*Nf) + w2*e1c(ij)
    e1f(ijf+1+Nf*(Nf-1)) = e1f(ijf+1+Nf*(Nf-1)) + w2*e1c(ij)
    e1f(ijf+1+2*Nf) = e1f(ijf+1+2*Nf) + w2*e1c(ij)
    e1f(ijf+2-Nf) = e1f(ijf+2-Nf) + w2*e1c(ij)
    e1f(ijf+2) = e1f(ijf+2) + w2*e1c(ij)

    ! least contribution to furthest
    e1f(ijf-1+Nf*(Nf-1)) = e1f(ijf-1+Nf*(Nf-1)) + w3*e1c(ij)
    e1f(ijf-1+2*Nf) = e1f(ijf-1+2*Nf) + w3*e1c(ij)
    e1f(ijf+2+Nf*(Nf-2)) = e1f(ijf+2+Nf*(Nf-2)) + w3*e1c(ij)
    e1f(ijf+2+Nf) = e1f(ijf+2+Nf) + w3*e1c(ij)

    ! largest contribution to nearest
    e2f(ijf) = e2f(ijf) + w1*e2c(ij)
    e2f(ijf+1) = e2f(ijf+1) + w1*e2c(ij)
    e2f(ijf+Nf) = e2f(ijf+Nf) + w1*e2c(ij)
    e2f(ijf+Nf+1) = e2f(ijf+Nf+1) + w1*e2c(ij)

    ! lesser contribution to intermediate
    e2f(ijf-1) = e2f(ijf-1) + w2*e2c(ij)
    e2f(ijf+Nf-1) = e2f(ijf+Nf-1) + w2*e2c(ij)
    e2f(ijf+Nf*(Nf-1)) = e2f(ijf+Nf*(Nf-1)) + w2*e2c(ij)
    e2f(ijf+2*Nf) = e2f(ijf+2*Nf) + w2*e2c(ij)
    e2f(ijf+1+Nf*(Nf-1)) = e2f(ijf+1+Nf*(Nf-1)) + w2*e2c(ij)
    e2f(ijf+1+2*Nf) = e2f(ijf+1+2*Nf) + w2*e2c(ij)
    e2f(ijf+2-Nf) = e2f(ijf+2-Nf) + w2*e2c(ij)
    e2f(ijf+2) = e2f(ijf+2) + w2*e2c(ij)

    ! least contribution to furthest
    e2f(ijf-1+Nf*(Nf-1)) = e2f(ijf-1+Nf*(Nf-1)) + w3*e2c(ij)
    e2f(ijf-1+2*Nf) = e2f(ijf-1+2*Nf) + w3*e2c(ij)
    e2f(ijf+2+Nf*(Nf-2)) = e2f(ijf+2+Nf*(Nf-2)) + w3*e2c(ij)
    e2f(ijf+2+Nf) = e2f(ijf+2+Nf) + w3*e2c(ij)


    ! bottom left corner
    j = N
    i = 1
    jf = 2*j-1
    if = 2*i-1
    ij = i+N*(j-1)
    ijf = if+Nf*(jf-1)

    ! largest contribution to nearest
    e1f(ijf) = e1f(ijf) + w1*e1c(ij)
    e1f(ijf+1) = e1f(ijf+1) + w1*e1c(ij)
    e1f(ijf+Nf) = e1f(ijf+Nf) + w1*e1c(ij)
    e1f(ijf+Nf+1) = e1f(ijf+Nf+1) + w1*e1c(ij)

    ! lesser contribution to intermediate
    e1f(ijf-1+Nf) = e1f(ijf-1+Nf) + w2*e1c(ij)
    e1f(ijf+2*Nf-1) = e1f(ijf+2*Nf-1) + w2*e1c(ij)
    e1f(ijf-Nf) = e1f(ijf-Nf) + w2*e1c(ij)
    e1f(if) = e1f(if) + w2*e1c(ij)
    e1f(ijf+1-Nf) = e1f(ijf+1-Nf) + w2*e1c(ij)
    e1f(if+1) = e1f(if+1) + w2*e1c(ij)
    e1f(ijf+2) = e1f(ijf+2) + w2*e1c(ij)
    e1f(ijf+2+Nf) = e1f(ijf+2+Nf) + w2*e1c(ij)

    ! least contribution to furthest
    e1f(ijf-1) = e1f(ijf-1) + w3*e1c(ij)
    e1f(if-1+Nf) = e1f(if-1+Nf) + w3*e1c(ij)
    e1f(ijf+2-Nf) = e1f(ijf+2-Nf) + w3*e1c(ij)
    e1f(if+2) = e1f(if+2) + w3*e1c(ij)

    ! largest contribution to nearest
    e2f(ijf) = e2f(ijf) + w1*e2c(ij)
    e2f(ijf+1) = e2f(ijf+1) + w1*e2c(ij)
    e2f(ijf+Nf) = e2f(ijf+Nf) + w1*e2c(ij)
    e2f(ijf+Nf+1) = e2f(ijf+Nf+1) + w1*e2c(ij)

    ! lesser contribution to intermediate
    e2f(ijf-1+Nf) = e2f(ijf-1+Nf) + w2*e2c(ij)
    e2f(ijf+2*Nf-1) = e2f(ijf+2*Nf-1) + w2*e2c(ij)
    e2f(ijf-Nf) = e2f(ijf-Nf) + w2*e2c(ij)
    e2f(if) = e2f(if) + w2*e2c(ij)
    e2f(ijf+1-Nf) = e2f(ijf+1-Nf) + w2*e2c(ij)
    e2f(if+1) = e2f(if+1) + w2*e2c(ij)
    e2f(ijf+2) = e2f(ijf+2) + w2*e2c(ij)
    e2f(ijf+2+Nf) = e2f(ijf+2+Nf) + w2*e2c(ij)

    ! least contribution to furthest
    e2f(ijf-1) = e2f(ijf-1) + w3*e2c(ij)
    e2f(if-1+Nf) = e2f(if-1+Nf) + w3*e2c(ij)
    e2f(ijf+2-Nf) = e2f(ijf+2-Nf) + w3*e2c(ij)
    e2f(if+2) = e2f(if+2) + w3*e2c(ij)


    ! bottom right corner
    j = N
    i = N
    jf = 2*j-1
    if = 2*i-1
    ij = i+N*(j-1)
    ijf = if+Nf*(jf-1)

    ! largest contribution to nearest
    e1f(ijf) = e1f(ijf) + w1*e1c(ij)
    e1f(ijf+1) = e1f(ijf+1) + w1*e1c(ij)
    e1f(ijf+Nf) = e1f(ijf+Nf) + w1*e1c(ij)
    e1f(ijf+Nf+1) = e1f(ijf+Nf+1) + w1*e1c(ij)

    ! lesser contribution to intermediate
    e1f(ijf-1) = e1f(ijf-1) + w2*e1c(ij)
    e1f(ijf+Nf-1) = e1f(ijf+Nf-1) + w2*e1c(ij)
    e1f(ijf-Nf) = e1f(ijf-Nf) + w2*e1c(ij)
    e1f(if) = e1f(if) + w2*e1c(ij)
    e1f(ijf+1-Nf) = e1f(ijf+1-Nf) + w2*e1c(ij)
    e1f(if+1) = e1f(if+1) + w2*e1c(ij)
    e1f(ijf+2-Nf) = e1f(ijf+2-Nf) + w2*e1c(ij)
    e1f(ijf+2) = e1f(ijf+2) + w2*e1c(ij)

    ! least contribution to furthest
    e1f(ijf-1-Nf) = e1f(ijf-1-Nf) + w3*e1c(ij)
    e1f(if-1) = e1f(if-1) + w3*e1c(ij)
    e1f(ijf+2-2*Nf) = e1f(ijf+2-2*Nf) + w3*e1c(ij)
    e1f(if+2-Nf) = e1f(if+2-Nf) + w3*e1c(ij)

    ! largest contribution to nearest
    e2f(ijf) = e2f(ijf) + w1*e2c(ij)
    e2f(ijf+1) = e2f(ijf+1) + w1*e2c(ij)
    e2f(ijf+Nf) = e2f(ijf+Nf) + w1*e2c(ij)
    e2f(ijf+Nf+1) = e2f(ijf+Nf+1) + w1*e2c(ij)

    ! lesser contribution to intermediate
    e2f(ijf-1) = e2f(ijf-1) + w2*e2c(ij)
    e2f(ijf+Nf-1) = e2f(ijf+Nf-1) + w2*e2c(ij)
    e2f(ijf-Nf) = e2f(ijf-Nf) + w2*e2c(ij)
    e2f(if) = e2f(if) + w2*e2c(ij)
    e2f(ijf+1-Nf) = e2f(ijf+1-Nf) + w2*e2c(ij)
    e2f(if+1) = e2f(if+1) + w2*e2c(ij)
    e2f(ijf+2-Nf) = e2f(ijf+2-Nf) + w2*e2c(ij)
    e2f(ijf+2) = e2f(ijf+2) + w2*e2c(ij)

    ! least contribution to furthest
    e2f(ijf-1-Nf) = e2f(ijf-1-Nf) + w3*e2c(ij)
    e2f(if-1) = e2f(if-1) + w3*e2c(ij)
    e2f(ijf+2-2*Nf) = e2f(ijf+2-2*Nf) + w3*e2c(ij)
    e2f(if+2-Nf) = e2f(if+2-Nf) + w3*e2c(ij)
  end subroutine prolongate


  ! TODO: remove, for testing only
  subroutine temp_output_data(b, fname)
    implicit none
    real(dp), pointer, contiguous, intent(in) :: b(:)
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
