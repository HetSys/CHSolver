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

  subroutine solver_pssi(A, Tout, CH_params, inarr, eps2, errors, inarr1, dt_in)
    implicit none

    real(dp), intent(in)                              :: Tout(:)
    real(dp), intent(in)                              :: eps2, A
    real(dp), intent(in), dimension(6)                :: CH_params
    real(dp),              intent(in)                 :: inarr(:,:)
    real(dp), intent(in), optional                    :: inarr1(:,:)
    real(dp), intent(in), optional                    :: dt_in
    integer, intent(inout)                            :: errors
    integer                                           :: N, i, n_threads, j, k, caseflag
    type(C_PTR)                                       :: fwplan, bwplan, pc0, pc1
    real(dp), dimension(:, :), allocatable            :: ksq
    complex(cdc), dimension(:,:), pointer             :: c0, c1
    real(dp), allocatable, dimension(:,:)             :: rc0, rc1
    real(dp)                                          :: dt1, dt2, t, tmax, dt, dt0, t_out, dt_out, kappa, origdt
    logical                                           :: outflag
    character(len=48)                                 :: msg
    character(1), parameter                           :: newline = NEW_LINE('a')

    ! SETUP=====================================================================
    24 format(A, F7.4)
    N = size(inarr, 1)
    tmax = maxval(Tout)
    i = 1
    t = 0.0_dp
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

    dt1 = 0.1_dp*kappa
    dt2 = 0.1_dp/real(N, dp)

    if ( dt1 < dt2 ) then
      dt = dt1
    else
      dt = dt2
    end if
    origdt = dt

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

      if (tout(i) < epsilon(tout(i))) then
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

        call write_to_traj(real(rc1,kind=dp), real(rc0,kind=dp), t_out, dt_out, errors)
        i = i + 1
      endif

    ! INITAL TIMESTEP===========================================================
      if ( t + dt + epsilon(t) > Tout(i) ) then
        dt = tout(i) - t
        outflag = .true.
      else
        outflag = .false.
      end if

      t = t + dt
      dt0 = dt

      call initial_iteration(dt, kappa, ksq, c0, c1, fwplan, bwplan)
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

        call write_to_traj(real(rc1,kind=dp), real(rc0,kind=dp), t_out, dt_out, errors)

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

      if (tout(i) < epsilon(tout(i))) then
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

        call write_to_traj(real(rc1,kind=dp), real(rc0,kind=dp), t_out, dt_out, errors)
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

        call write_to_traj(real(rc1,kind=dp), real(rc0,kind=dp), t_out, dt_out, errors)

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
  !! @param[in]
  !! @param[inout]
  subroutine f(c, res)
    complex(cdc), dimension(:,:), intent(in)   :: c
    complex(cdc), dimension(:,:), intent(out)  :: res
    !$OMP PARALLEL WORKSHARE
    res(:,:) = c(:,:)*c(:,:)*c(:,:) - c(:,:)
    !$OMP END PARALLEL WORKSHARE
  end subroutine f

  !> @brief Perform a finite difference time step in k-space.
  !!
  !! @param[in]
  !! @param[inout]
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

  !> @brief Perform the initial iteration of the pseudo spectral method.
  !!
  !! @param[in] tau, kappa, A, ksq, fwplan, bwplan  Stuff here.
  !! @param[inout] c0, res  The concentration grids.
  subroutine initial_iteration(tau, kappa, ksq, c0, res, fwplan, bwplan)
    complex(cdc), dimension(:,:), intent(inout)  :: c0
    real(dp), intent(in)                                  :: tau, kappa
    real(dp), dimension(:,:), intent(in)                  :: ksq
    complex(cdc), dimension(:,:), intent(out)    :: res
    type(C_PTR), intent(in)                                   :: fwplan, bwplan
    integer, dimension(2)                                     :: c_shape
    integer                                                   :: N
    real(dp)                                              :: tau1
    complex(cdc)                                           :: ctau1, ckappa
    complex(cdc), allocatable ,dimension(:,:)              :: cksq
    complex(cdc), pointer :: fc0(:,:), ft_c0(:,:), ft_fc0(:,:), ft_c1(:,:)
    type(C_PTR) :: p1, p2, p3, p4


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

  !> @brief Perform one iteration of the pseudo spectral method.
  !!
  !! @param[in] tau, kappa, A, ksq, fwplan, bwplan  Stuff here.
  !! @param[inout] c0, c1  The concentration grids.
  subroutine iteration(tau, ksq, kappa, A, c0, c1, fwplan, bwplan)
    complex(cdc), dimension(:,:), intent(inout)  :: c1, c0
    real(dp), intent(in)                                  :: tau, kappa, A
    real(dp), dimension(:,:), intent(in)                  :: ksq
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
  !! @param[out] x2, y2  The meshgrids.
  subroutine meshgrid(x, y, x2, y2)
    real(dp), intent(in) :: x(:), y(:)
    real(dp), intent(out) :: x2(:, :), y2(:, :)

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
