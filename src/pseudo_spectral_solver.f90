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
    real(dp), allocatable, intent(in)                 :: inarr(:,:)
    real(dp), allocatable, intent(in), optional       :: inarr1(:,:)
    real(dp), intent(in), optional                    :: dt_in
    integer, intent(inout)                            :: errors
    integer                                           :: N, i, n_threads, j, k, caseflag
    type(C_PTR)                                       :: fwplan, bwplan, pc0, pc1
    real(dp), dimension(:, :), allocatable            :: ksq
    complex(cdc), dimension(:,:), pointer:: c0, c1
    real(dp), allocatable, dimension(:,:)             :: rc0, rc1
    real(dp)                                          :: dt1, dt2, t, tmax, dt, dt0, t_out, dt_out, kappa, origdt
    logical :: outflag
    character(len=48) :: msg
    character(1), parameter :: newline = NEW_LINE('a')

    ! SETUP=====================================================================
    24 format(A, F7.4)
    N = size(inarr, 1)
    tmax = maxval(Tout)
    i = 1
    t = 0.0_dp
    kappa = sqrt(eps2)
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
    if (caseflag == 1) then

      do j = 1, N
        do k = 1, N
          c0(j,k) = cmplx(inarr(j, k), kind=cdc)
        end do
      end do

      if (tout(i) < epsilon(tout(i))) then
        write(msg, 24) "Initial condition output at t=  ", t
        call logger%info("solver_pssi", msg)
        dt_out = dt
        t_out = t

        do j = 1, N
          do k = 1, N
            rc0(j,k) = real(c0(j, k), kind=dp)
            rc1(j,k) = real(c1(j, k), kind=dp)
          end do
        end do

        call dimensionalise(CH_params, rc1, t_out)
        call dimensionalise(CH_params, rc0, dt_out)

        call write_to_traj_2D(real(rc1,kind=dp), real(rc0,kind=dp), t_out, dt_out, errors)
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

        do j = 1, N
          do k = 1, N
            rc0(j,k) = real(c0(j, k), kind=dp)
            rc1(j,k) = real(c1(j, k), kind=dp)
          end do
        end do

        call dimensionalise(CH_params, rc1, t_out)
        call dimensionalise(CH_params, rc0, dt_out)

        call write_to_traj_2D(real(rc1,kind=dp), real(rc0,kind=dp), t_out, dt_out, errors)

        i = i + 1
      end if

    else if (caseflag == 0) then

      dt = dt_in

      do j = 1, N
        do k = 1, N
          c0(j,k) = cmplx(inarr(j, k), kind=cdc)
          c1(j,k) = cmplx(inarr1(j, k), kind=cdc)
        end do
      end do

      if (tout(i) < epsilon(tout(i))) then
        write(msg, 24) "Initial from checkpoint output at t=  ", t
        call logger%info("solver_pssi", msg)
        dt_out = dt
        t_out = t

        do j = 1, N
          do k = 1, N
            rc0(j,k) = real(c0(j, k), kind=dp)
            rc1(j,k) = real(c1(j, k), kind=dp)
          end do
        end do

        call dimensionalise(CH_params, rc1, t_out)
        call dimensionalise(CH_params, rc0, dt_out)

        call write_to_traj_2D(real(rc1,kind=dp), real(rc0,kind=dp), t_out, dt_out, errors)
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
      !print*, t, dt

      call iteration(dt, ksq, kappa, A, c0, c1, fwplan, bwplan)
      dt = origdt

      if (outflag) then
        write(msg, 24) "Output at t=", t
        call logger%info("solver_pssi", msg)

        dt_out = dt
        t_out = t

        do j = 1, N
          do k = 1, N
            rc0(j,k) = real(c0(j, k), kind=dp)
            rc1(j,k) = real(c1(j, k), kind=dp)
          end do
        end do

        call dimensionalise(CH_params, rc1, t_out)
        call dimensionalise(CH_params, rc0, dt_out)

        call write_to_traj_2D(real(rc1,kind=dp), real(rc0,kind=dp), t_out, dt_out, errors)

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


end module pseudo_spectral_solver
