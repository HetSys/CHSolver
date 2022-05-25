!> @brief Module defining the subroutines required to run the Pseudo-spectral solver.
!! @details It is important that any arrays passed to FFTW3 subroutine calls are
!! initialised with the FFTW alloc functions as described on the FFTW3
!! <a href="https://www.fftw.org/fftw3_doc/Memory-Allocation.html#Memory-Allocation">website</a>.
module pseudo_spectral_solver
  use, intrinsic :: iso_fortran_env
  use, intrinsic :: iso_c_binding
  use omp_lib
  use fftw3
  use solver_utils
  use globals
  use hdf5_io

  implicit none

contains

  !> @brief Solve the Cahn-Hilliard equation using a pseudo-spectral semi-implicit method.
  !! @details This subroutine solves the Cahn-Hilliard equation \cite{Li2017} using a
  !! semi-implicit (explicit second order adams-bashforth \cite{Li2017} for the non-linear
  !! part, second order implicit for the linear part) pseudo-spectral (spectral
  !! in space, finite difference in time) scheme. The equation can be optionally
  !! stabilised using a term proportional to grad-squared from \cite{Zhang2019}.
  !! @param[in] t0 The time at which to start the solver.
  !! @param[in] A Stabilisation parameter, as seen in \cite{Zhang2019}. Larger values will
  !! significantly disrupt initial dynamics, but will allow the use of large
  !! (up to order 1!) timesteps to reach the steady state solution.
  !! @param[in] Tout An array of times at which the solver will write hdf5 checkpoint files.
  !! @param[in] CH_params An array storing the 6 parameters defining the CH equation in
  !! the project description. These parameters are used to convert the non-dimensionalised
  !! problem solved by the solver to the dimensionalised version saved in the checkpoint file.
  !! @param[in] inarr The input concentration array.
  !! @param[in] eps2 The parameter that controls the behaviour of the non-dimensionalised
  !! CH equation. This can be found with a call to non_dimensionalise().
  !! @param[inout] errors A value to catch internal FFTW3 errors. A value of 1 (strangely) means
  !! FFTW3 has been initialised with openmp correctly, while a value of 0 indicates an error.
  !! @param[in] inarr1 An optional second concentration array, defined to be dt_in
  !! further ahead in time than inarr. If both inarr1 and dt_in are specified,
  !! the solver will begin from the coupled initial condition.
  !! @param[in] dt_in An optional real value, defining the temporal separation between inarr
  !! and inarr1, used when evaluating a system beginning from a paired initial condition.
  subroutine solver_pssi(t0, A, Tout, CH_params, inarr, eps2, errors, inarr1, dt_in)
    implicit none

    real(dp), intent(in)                              :: Tout(:)
    real(dp), intent(in)                              :: eps2, A
    real(dp), intent(in), dimension(6)                :: CH_params
    real(dp), intent(in)                              :: inarr(:,:)
    real(dp), intent(in), optional                    :: inarr1(:,:)
    real(dp), intent(in), optional                    :: dt_in
    real(dp), intent(in)                              :: t0
    integer, intent(inout)                            :: errors
    integer                                           :: N, i, n_threads, j, k, caseflag
    type(C_PTR)                                       :: fwplan, bwplan, pc0, pc1
    real(dp), dimension(:, :), allocatable            :: ksq
    complex(cdc), dimension(:,:), pointer             :: c0, c1
    real(dp), allocatable, dimension(:,:)             :: rc0, rc1
    real(dp)                                          :: dt1, dt2, t, tmax, dt, dt0, t_out, dt_out, kappa, origdt, initdt
    logical                                           :: outflag
    character(len=42)                                 :: msg
    character(1), parameter                           :: newline = NEW_LINE('a')

    ! SETUP=====================================================================
    24 format(A, F7.4)
    N = size(inarr, 1)
    tmax = maxval(Tout)
    i = 1
    t = t0
    kappa = eps2
    outflag = .true.
    caseflag = 0

    if (present(inarr1) .and. present(dt_in)) then
      caseflag = 0
    else if (present(inarr1)) then
      call logger%warning("solver_ps", "c1 was provided while dt isnt specified"// &
        newline//"defaulting to first order solver for initial")
      caseflag = 1
    else if (present(dt_in)) then
      call logger%warning("solver_ps", "dt was specified when c1 was not provided"// &
          newline//"defaulting to first order solver for initial")
      caseflag = 1
    else
      caseflag = 1
    end if

    dt1 = 0.45_dp*kappa
    dt2 = 0.45_dp/real(N, dp)

    if ( dt1 < dt2 ) then
      dt = dt1
    else
      dt = dt2
    end if
    origdt = dt
    initdt = dt*1e-3_dp

    allocate(rc0(N,N))
    allocate(rc1(N,N))

    allocate(ksq(N,N))
    call create_ksq(ksq, N)

    errors = fftw_init_threads()
    n_threads = omp_get_max_threads()
    call fftw_plan_with_nthreads(int(n_threads))

    pc0 = fftw_alloc_complex(int(N * N, C_SIZE_T))
    call c_f_pointer(pc0, c0, [N, N])
    pc1 = fftw_alloc_complex(int(N * N, C_SIZE_T))
    call c_f_pointer(pc1, c1, [N, N])

    fwplan = fftw_plan_dft_2d(N, N, c0, c1, FFTW_FORWARD, FFTW_ESTIMATE)
    bwplan = fftw_plan_dft_2d(N, N, c0, c1, FFTW_BACKWARD, FFTW_ESTIMATE)

    ! START FROM BEGINNING======================================================
    if (caseflag == 1) then
      !$OMP PARALLEL PRIVATE(j, k)
      do j = 1, N
        do k = 1, N
          c0(j,k) = cmplx(inarr(j, k), kind=cdc)
        end do
      end do
      !$OMP END PARALLEL

      if (abs(tout(i)-t) < epsilon(tout(i))) then
        write(msg, 24) "Initial condition output at t=  ", t
        call logger%info("solver_pssi", msg)
        dt_out = dt
        t_out = t
        !$OMP PARALLEL PRIVATE(j, k)
        do j = 1, N
          do k = 1, N
            rc0(j,k) = real(c0(j, k), kind=dp)
            rc1(j,k) = real(c1(j, k), kind=dp)
          end do
        end do
        !$OMP END PARALLEL

        call dimensionalise(CH_params, rc1, t_out)
        call dimensionalise(CH_params, rc0, dt_out)

        call write_to_traj(rc0, rc1, t_out, dt_out, errors)
        i = i + 1
      endif

    ! INITAL TIMESTEP===========================================================
      if ( t + dt + epsilon(t) > Tout(i) ) then
        dt = tout(i) - t
        outflag = .true.
      else
        outflag = .false.
      end if

      t = t + initdt
      dt0 = initdt

      call initial_iteration(initdt, kappa, ksq, c0, c1, fwplan, bwplan)

      dt = origdt

      if (outflag) then
        write(msg, 24) "Output at t=", t
        call logger%info("solver_pssi", msg)

        dt_out = dt
        t_out = t

        !$OMP PARALLEL PRIVATE(j, k)
        do j = 1, N
          do k = 1, N
            rc0(j,k) = real(c0(j, k), kind=dp)
            rc1(j,k) = real(c1(j, k), kind=dp)
          end do
        end do
        !$OMP END PARALLEL

        call dimensionalise(CH_params, rc1, t_out)
        call dimensionalise(CH_params, rc0, dt_out)

        call write_to_traj(rc1, rc0, t_out, dt_out, errors)

        i = i + 1
      end if

    ! START FROM CHECKPOINT=====================================================
    else if (caseflag == 0) then

      dt = dt_in

      !$OMP PARALLEL PRIVATE(j, k)
      do j = 1, N
        do k = 1, N
          c0(j,k) = cmplx(inarr(j, k), kind=cdc)
          c1(j,k) = cmplx(inarr1(j, k), kind=cdc)
        end do
      end do
      !$OMP END PARALLEL

      if (abs(tout(i)-t) < epsilon(tout(i))) then
        write(msg, 24) "Initial from checkpoint output at t=  ", t
        call logger%info("solver_pssi", msg)
        dt_out = dt
        t_out = t

        !$OMP PARALLEL PRIVATE(j, k)
        do j = 1, N
          do k = 1, N
            rc0(j,k) = real(c0(j, k), kind=dp)
            rc1(j,k) = real(c1(j, k), kind=dp)
          end do
        end do
        !$OMP END PARALLEL

        call dimensionalise(CH_params, rc1, t_out)
        call dimensionalise(CH_params, rc0, dt_out)

        call write_to_traj(rc1, rc0, t_out, dt_out, errors)
        i = i + 1
      end if

    end if

    ! REMAINING TIMESTEPS=======================================================
    do while (t < tmax)
      if ( t + dt +epsilon(t) > Tout(i) ) then
        dt = Tout(i) - t
        outflag = .true.
      else
        outflag = .false.
      end if
      t = t + dt
      dt0 = dt

      call iteration(dt, ksq, kappa, A, c0, c1, fwplan, bwplan)
      dt = origdt

      if (outflag) then
        write(msg, 24) "Output at t=", t
        call logger%info("solver_pssi", msg)

        dt_out = dt
        t_out = t

        !$OMP PARALLEL PRIVATE(j, k)
        do j = 1, N
          do k = 1, N
            rc0(j,k) = real(c0(j, k), kind=dp)
            rc1(j,k) = real(c1(j, k), kind=dp)
          end do
        end do
        !$OMP END PARALLEL

        call dimensionalise(CH_params, rc1, t_out)
        call dimensionalise(CH_params, rc0, dt_out)

        call write_to_traj(rc1, rc0, t_out, dt_out, errors)

        i = i + 1
      end if

    end do

    ! CLEANUP===================================================================
    call fftw_destroy_plan(fwplan)
    call fftw_destroy_plan(bwplan)

    call fftw_cleanup_threads()

    call fftw_free(pc0)
    call fftw_free(pc1)
    deallocate(ksq)
    deallocate(rc0)
    deallocate(rc1)

  end subroutine solver_pssi

  !> @brief Calculate the bulk free potential energy.
  !!
  !! @param[in] c The concentration array.
  !! @param[out] res The bulk free energy at that concentration.
  subroutine f(c, res)
    complex(cdc), dimension(:,:), intent(in)   :: c
    complex(cdc), dimension(:,:), intent(out)  :: res
    !$OMP PARALLEL WORKSHARE
    res(:,:) = c(:,:)*c(:,:)*c(:,:) - c(:,:)
    !$OMP END PARALLEL WORKSHARE
  end subroutine f

  !> @brief Takes one second order semi-implicit time step in k-space.
  !! @details None
  !! @param[in] tau The temporal finite difference timestep.
  !! @param[in] ksq The matrix of k-values used to replace the spatial derivatives
  !! in the spectral method.
  !! @param[in] kappa The dimensionless parameter that controls the behaviour of the
  !! CH equation.
  !! @param[in] A Stabilisation parameter, as seen in \cite{Zhang2019}. Larger values will
  !! significantly disrupt initial dynamics, but will allow the use of large
  !! (order 1) timesteps to reach the steady state solution.
  !! @param[in] ft_c1 The fourier transform of c1, the concentration array calculated
  !! 1 time step previously.
  !! @param[in] ft_c0 The fourier transform of c0, the concentration array calculated
  !! 2 time steps previously.
  !! @param[in] ft_fc1 The fourier transform of the bulk potential calculated for c1.
  !! @param[in] ft_fc0 The fourier transform of the bulk potential calculated for c0.
  !! @param[out] res The concentration array at time tau later than c1.
  subroutine step(tau, ksq, kappa, A, ft_c1, ft_c0, ft_fc1, ft_fc0, res)
    complex(cdc), dimension(:,:), intent(in)               :: ft_c1, ft_c0, ft_fc1, ft_fc0
    real(dp), intent(in)                                   :: tau, kappa, A
    real(dp), dimension(:,:), intent(in)                   :: ksq
    complex(cdc), dimension(:,:), intent(out)              :: res
    complex(cdc)                                           :: ctau, ckappa, cA
    complex(cdc), allocatable, dimension(:,:)              :: cksq
    integer                                                :: N
    N = size(ft_c1, 1)
    allocate(cksq(N,N))
    ctau = cmplx(tau, kind=cdc)
    ckappa = cmplx(kappa, kind=cdc)
    cA = cmplx(A, kind=cdc)

    !$OMP PARALLEL WORKSHARE
    cksq(:,:) = cmplx(ksq(:,:), kind=cdc)
    res(:,:) = (cmplx(2.0_dp,kind=cdc)*ctau/(cmplx(3.0_dp,kind=cdc)+cmplx(2.0_dp,kind=cdc)&
    *ctau*ckappa*cksq(:,:)*cksq(:,:)-cmplx(2.0_dp,kind=cdc)*ctau*cA*cksq(:,:)))&
    *(cmplx(2.0_dp,kind=cdc)*cksq(:,:)*ft_fc1(:,:)-cksq(:,:)*ft_fc0(:,:)-&
    cksq(:,:)*cmplx(2.0_dp,kind=cdc)*cA*ft_c1(:,:)+cksq(:,:)*cA*&
    ft_c0(:,:)+(cmplx(4.0_dp,kind=cdc)*ft_c1(:,:)-ft_c0(:,:))/(cmplx(2.0_dp,kind=cdc)*ctau))
    !$OMP END PARALLEL WORKSHARE
    deallocate(cksq)
  end subroutine step

  !> @brief The initial iteration for the ps solver, starting from a single concentration
  !! array c0.
  !! @details This iteration is an explicit first order in time pseudo-spectral step.
  !! As such, the time-step tolerance is much worse than for iteration(). As a rule
  !! of thumb, it should be 100 times smaller.
  !! @param[in] tau The temporal finite difference timestep.
  !! @param[in] kappa The dimensionless parameter that controls the behaviour of the
  !! CH equation.
  !! @param[in] ksq The matrix of k-values used to replace the spatial derivatives in
  !! the spectral method.
  !! @param[inout] c0 The initial concentration array. INOUT to be type-compatible with
  !! FFTW3 plans.
  !! @param[out] res The concentration array at a time tau later than c0.
  !! @param[in] fwplan The FFTW3 plan used to forward fourier transform the concentration
  !! arrays and the bulk potentials.
  !! @param[in] bwplan The FFTW3 plan used to backward fourier transform the concentration
  !! arrays and the bulk potentials.
  subroutine initial_iteration(tau, kappa, ksq, c0, res, fwplan, bwplan)
    complex(cdc), dimension(:,:), intent(inout)           :: c0
    real(dp), intent(in)                                  :: tau, kappa
    real(dp), dimension(:,:), intent(in)                  :: ksq
    complex(cdc), dimension(:,:), intent(out)             :: res
    type(C_PTR), intent(in)                               :: fwplan, bwplan
    integer, dimension(2)                                 :: c_shape
    integer                                               :: N
    real(dp)                                              :: tau1
    complex(cdc)                                          :: ctau1, ckappa
    complex(cdc), allocatable ,dimension(:,:)             :: cksq
    complex(cdc), pointer                                 :: fc0(:,:), ft_c0(:,:), ft_fc0(:,:), ft_c1(:,:)
    type(C_PTR)                                           :: p1, p2, p3, p4


    c_shape = shape(c0)
    N = c_shape(1)
    allocate(cksq(N,N))

    p1 = fftw_alloc_complex(int(N * N, C_SIZE_T))
    call c_f_pointer(p1, fc0, [N, N])
    p2 = fftw_alloc_complex(int(N * N, C_SIZE_T))
    call c_f_pointer(p2, ft_c0, [N, N])
    p3 = fftw_alloc_complex(int(N * N, C_SIZE_T))
    call c_f_pointer(p3, ft_fc0, [N, N])
    p4 = fftw_alloc_complex(int(N * N, C_SIZE_T))
    call c_f_pointer(p4, ft_c1, [N, N])

    call f(c0, fc0)
    call fftw_execute_dft(fwplan, c0, ft_c0)
    call fftw_execute_dft(fwplan, fc0, ft_fc0)
    if (tau < 1.0_dp) then
      tau1 = tau
    else
      tau1 = 1.0_dp
    end if
    tau1 = tau*1e-6_dp

    ctau1 = cmplx(tau1, kind=cdc)
    ckappa = cmplx(kappa, kind=cdc)


    !$OMP PARALLEL WORKSHARE
    cksq(:,:) = cmplx(ksq(:,:), kind=cdc)
    ft_c1(:,:) = ctau1*(-ckappa*cksq(:,:)*cksq(:,:)*ft_c0(:,:)+cksq(:,:)*ft_fc0(:,:))+ft_c0(:,:)
    !$OMP END PARALLEL WORKSHARE

    call fftw_execute_dft(bwplan, ft_c1, res)

    !$OMP PARALLEL WORKSHARE
    res(:,:) = res(:,:)/cmplx(N*N, kind=cdc)
    !$OMP END PARALLEL WORKSHARE

    call fftw_free(p1)
    call fftw_free(p2)
    call fftw_free(p3)
    call fftw_free(p4)
    deallocate(cksq)
  end subroutine initial_iteration

  !> @brief The iteration for the ps solver, starting from a paired initial condition of
  !! c0 and c1 separated by timestep tau.
  !! @details This iteration is a second order in time semi-implicit pseudo-spectral
  !! step. As a rule of thumb, tau should be the lowest of 0.3/N and 0.3*kappa.
  !! At the output of this function, c0 will be overwritten with the input value of c1
  !! while c1 will store the output of the iteration.
  !! @param[in] tau The temporal finite difference timestep.
  !! @param[in] kappa The dimensionless parameter that controls the behaviour of the
  !! CH equation.
  !! @param[in] ksq The matrix of k-values used to replace the spatial derivatives in
  !! the spectral method.
  !! @param[in] A Stabilisation parameter, as seen in \cite{Zhang2019}. Larger values will
  !! significantly disrupt initial dynamics, but will allow the use of large
  !! (order 1) timesteps to reach the steady state solution.
  !! @param[inout] c0 The concentration array at a time tau before c1.
  !! @param[inout] c1 The concentration array at a time tau before the output.
  !! @param[in] fwplan The FFTW3 plan used to forward fourier transform the concentration
  !! arrays and the bulk potentials.
  !! @param[in] bwplan The FFTW3 plan used to backward fourier transform the concentration
  !! arrays and the bulk potentials.
  subroutine iteration(tau, ksq, kappa, A, c0, c1, fwplan, bwplan)
    complex(cdc), dimension(:,:), intent(inout)               :: c1, c0
    real(dp), intent(in)                                      :: tau, kappa, A
    real(dp), dimension(:,:), intent(in)                      :: ksq
    type(C_PTR), intent(in)                                   :: fwplan, bwplan
    integer, dimension(2)                                     :: c_shape
    integer                                                   :: N

    complex(cdc), pointer :: fc1(:,:), fc0(:,:), ft_c1(:,:),&
    ft_fc1(:,:), ft_c0(:,:), ft_fc0(:,:), res(:,:), ft_res(:,:)
    type(C_PTR) :: p1, p2, p3, p4, p5, p6, p7, p8

    c_shape = shape(c0)
    N = c_shape(1)

    p1 = fftw_alloc_complex(int(N * N, C_SIZE_T))
    call c_f_pointer(p1, fc1, [N, N])
    p2 = fftw_alloc_complex(int(N * N, C_SIZE_T))
    call c_f_pointer(p2, fc0, [N, N])
    p3 = fftw_alloc_complex(int(N * N, C_SIZE_T))
    call c_f_pointer(p3, ft_c1, [N, N])
    p4 = fftw_alloc_complex(int(N * N, C_SIZE_T))
    call c_f_pointer(p4, ft_fc1, [N, N])
    p5 = fftw_alloc_complex(int(N * N, C_SIZE_T))
    call c_f_pointer(p5, ft_c0, [N, N])
    p6 = fftw_alloc_complex(int(N * N, C_SIZE_T))
    call c_f_pointer(p6, ft_fc0, [N, N])
    p7 = fftw_alloc_complex(int(N * N, C_SIZE_T))
    call c_f_pointer(p7, res, [N, N])
    p8 = fftw_alloc_complex(int(N * N, C_SIZE_T))
    call c_f_pointer(p8, ft_res, [N, N])

    call f(c1, fc1)
    call f(c0, fc0)

    call fftw_execute_dft(fwplan, c1, ft_c1)
    call fftw_execute_dft(fwplan, fc1, ft_fc1)

    call fftw_execute_dft(fwplan, c0, ft_c0)
    call fftw_execute_dft(fwplan, fc0, ft_fc0)

    call step(tau, ksq, kappa, A, ft_c1, ft_c0, ft_fc1, ft_fc0, ft_res)

    call fftw_execute_dft(bwplan, ft_res, res)
    !$OMP PARALLEL WORKSHARE
    res(:,:) = res(:,:)/cmplx(N*N, kind=cdc)
    !$OMP END PARALLEL WORKSHARE

    !$OMP PARALLEL WORKSHARE
    c0(:,:) = c1(:,:)
    c1(:,:) = res(:,:)
    !$OMP END PARALLEL WORKSHARE

    call fftw_free(p1)
    call fftw_free(p2)
    call fftw_free(p3)
    call fftw_free(p4)
    call fftw_free(p5)
    call fftw_free(p6)
    call fftw_free(p7)
    call fftw_free(p8)
  end subroutine iteration

  !> @brief Create a 2D meshgrid as in numpy's "meshgrid" function.
  !!
  !! @param[in] x, y  The coordinate series for the two dimensions.
  !! @param[inout] x2, y2  The meshgrids.
  subroutine meshgrid(x, y, x2, y2)
    real(dp), intent(in) :: x(:), y(:)
    real(dp), intent(inout) :: x2(:, :), y2(:, :)

    x2 = spread(x, 1, size(y))
    y2 = spread(y, 2, size(x))
  end subroutine meshgrid

  !> @brief Create the fourier (k-) space grid used to calculate real space derivatives.
  !!
  !! @param[inout] ksq  The squre array used to store the k-space grid.
  !! @param[out] N  The size of ksq.
  subroutine create_ksq(ksq, N)
    real(dp), parameter                            :: PI = 3.14159265_dp
    integer, intent(in)                            :: N
    real(dp), dimension(:, :), intent(inout)       :: ksq
    integer                                        :: i
    real(dp), dimension(:, :), allocatable         :: KX, KY
    real(dp), dimension(:), allocatable            :: k

    allocate(KX(N,N))
    allocate(KY(N,N))
    allocate(k(N))

    !$OMP PARALLEL PRIVATE(i)
    do i = 1, int(N/2)
      k(i) = real((i - 1), dp)
    end do
    do i = int(N/2) + 1, N
      k(i) =  real(i - N - 1, dp)
    end do
    !$OMP END PARALLEL

    !$OMP PARALLEL WORKSHARE
    k(:) = 2.0_dp*PI*k(:)!/real(N)
    !$OMP END PARALLEL WORKSHARE

    call meshgrid(k, k, KX, KY)
    !$OMP PARALLEL WORKSHARE
    ksq(:,:) = -(KX(:,:)*KX(:,:) + KY(:,:)*KY(:,:))
    !$OMP END PARALLEL WORKSHARE

    deallocate(KX)
    deallocate(KY)
    deallocate(k)
  end subroutine create_ksq

end module pseudo_spectral_solver
