module solver_utils
  use, intrinsic :: iso_fortran_env
  use, intrinsic :: iso_c_binding
  use globals
  use omp_lib
  use fftw3


  implicit none
  save

  contains

  !> @brief converts between dimensional and dimensionless systems
  !!
  !! @param[in] CH_params  array of 6 parameters describing the original system
  !! @param[inout] c       dimensional concentration -> dimensionless
  !! @param[inout] T       dimensional output times -> dimensionless
  !! @param[out] eps2      resulting dimensionless number
  subroutine nondimensionalise(CH_params, c, T, eps2)
    implicit none
    real(dp), intent(in) :: CH_params(6)
    real(dp), intent(inout) :: T(:)
    real(dp), dimension(:,:), allocatable, intent(inout) :: c
    real(dp), intent(out) :: eps2

    real(dp) :: L, A, M, K, p0, p1 ! CH_params contents
    real(dp) :: a0, a1 ! conversion factors

    L = CH_params(1)
    A = CH_params(2)
    M = CH_params(3)
    K = CH_params(4)
    p0 = CH_params(5)
    p1 = CH_params(6)

    a0 = abs(p1-p0)/2.0_dp
    a1 = (p0+p1)/2.0_dp

    ! scale concentrations so {p0,p1} -> {-1,1}
    c = (c-a1)/a0

    ! find dimensionless number
    eps2 = K/(4.0_dp*a0*a0*a0*A*L*L)

    ! scale times
    T = 4.0_dp*a0*a0*a0*A*M/(L*L) * T
  end subroutine nondimensionalise


  !> @brief converts between dimensionless and dimensional systems
  !!
  !! @param[in] CH_params  array of 6 parameters describing the original system
  !! @param[inout] c       dimensionless concentration -> dimensional
  !! @param[inout] t       dimensionless output time -> dimensional
  subroutine dimensionalise(CH_params, c, t)
    implicit none
    real(dp), intent(in) :: CH_params(6)
    real(dp), intent(inout) :: t
    real(dp), dimension(:,:), allocatable, intent(inout) :: c

    real(dp) :: L, A, M, K, p0, p1 ! CH_params contents
    real(dp) :: a0, a1 ! conversion factors

    L = CH_params(1)
    A = CH_params(2)
    M = CH_params(3)
    K = CH_params(4)
    p0 = CH_params(5)
    p1 = CH_params(6)

    a0 = abs(p1-p0)/2.0_dp
    a1 = (p0+p1)/2.0_dp

    ! scale concentrations so {-1,1} -> {p0,p1}
    c = a0*c+a1

    ! scale times
    t = L*L/(4.0_dp*a0*a0*a0*A*M) * t
  end subroutine dimensionalise


  !> @brief integer log2
  !!
  !! @param[in] val   2**k for some positive integer k
  !! @param[out] res  log2(val)
  subroutine ilog2(val, res)
    implicit none
    integer, intent(in) :: val
    integer, intent(out) :: res

    res = bit_size(val)-leadz(val)-1
  end subroutine ilog2

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

end module solver_utils
