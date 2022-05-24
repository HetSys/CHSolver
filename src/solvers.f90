module solvers
  use globals
  use solver_utils
  use fd_solvers
  use pseudo_spectral_solver

  implicit none
  save

  ! solver codes to select the desired solver
  integer, parameter :: SOLVER_FD2 = 1, SOLVER_PS = 2

  contains

  !> @brief solves the CH equation from a single initial condition
  !!
  !! @param[in] CH_params  array of 6 parameters describing the system
  !! @param[in] c0         initial concentration
  !! @param[in] Tout       output times
  !! @param[in] code       solver-specifier
  subroutine solver_1(Tout, c0, CH_params, code, error)
    implicit none
    integer, intent(out) :: error
    real(dp), intent(in) :: CH_params(6)
    real(dp), dimension(:,:), allocatable, intent(inout) :: c0
    real(dp), allocatable, intent(inout) :: Tout(:)
    integer, intent(in) :: code

    real(dp) :: eps2

    ! non-dimensionalise
    call nondimensionalise(CH_params, c0, Tout, eps2)

    ! call relevant solver
    select case (code)
      case (SOLVER_FD2)
        call logger%info("solver_1", "Solving with fd2")
        call solver_ufds2t2(0.0_dp, Tout, CH_params, c0, eps2, error)
      case (SOLVER_PS)
        call logger%info("solver_1", "Solving with ps")
        call solver_pssi(0.1_dp, Tout, CH_params, c0, eps2, error)
      case default
        call logger%fatal("solver_1", "Invalid solver code")
    endselect
  end subroutine solver_1


  !> @brief solves the CH equation from a paired initial condition
  !!
  !! @param[in] CH_params  array of 6 parameters describing the system
  !! @param[in] c0         initial concentration
  !! @param[in] c1         concentration at t = dt
  !! @param[in] dt         initial timestep
  !! @param[in] Tout       output times
  !! @param[in] code       solver-specifier
  subroutine solver_2(t0, Tout, c0, c1, dt, CH_params, code, error)
    implicit none
    integer :: error
    real(dp), intent(in) :: CH_params(6)
    real(dp) :: t0
    real(dp), allocatable, dimension(:,:), intent(inout) :: c0, c1
    real(dp), intent(in) :: dt
    real(dp), allocatable, intent(inout) :: Tout(:)
    integer, intent(in) :: code

    real(dp) :: eps2
    real(dp) :: dt_(1)
    dt_(1) = dt

    ! non-dimensionalise
    call nondimensionalise(CH_params, c0, Tout, eps2)
    call nondimensionalise(CH_params, c1, dt_, eps2)

    ! call relevant solver
    select case (code)
      case (SOLVER_FD2)
        call logger%info("solver_2", "Solving with fd2")
        call solver_ufds2t2(t0, Tout, CH_params, c0, eps2, error, c1, dt_(1))
      case (SOLVER_PS)
        call logger%info("solver_2", "Solving with ps")
        call solver_pssi(0.1_dp, Tout, CH_params, c0, eps2, error, c1, dt_(1))
      case default
        call logger%fatal("solver_2", "Invalid solver code")
    endselect

  end subroutine solver_2

end module solvers
