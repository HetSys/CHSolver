module solver_utils
  use globals
  use omp_lib


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
    real(dp), dimension(:,:), intent(inout) :: c
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
    real(dp), dimension(:,:), intent(inout) :: c

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
end module solver_utils
