module JSON_Parser
  !> Parses JSON file for given run name
  !! Returns the params needed for grid generation, and solving

  use ISO_FORTRAN_ENV
  implicit none
  

  character(*), parameter :: JSON_FILE = "input-data.json"
  contains


  !> @brief Validation of program inputs
  !! @param[in]  CH_params [L, A, M, K, p0, p1]
  !! @param[in]  grid_init Grid initialisation type character
  !! @param[in]  grid_level Controls size of grid
  subroutine validate(CH_params, grid_init, grid_level)
    real(kind=REAL64), intent(in) :: CH_params(6)
    character(1), intent(in) :: grid_init
    integer, intent(in) :: grid_level

  end subroutine

  !> @brief Reads JSON file, and searches for given params
  !! @param[in]  run_name  String key for input values
  !! @param[in]  CH_params [L, A, M, K, p0, p1]
  !! @param[in]  grid_init Grid initialisation type character
  !! @param[in]  grid_level Controls size of grid
  subroutine read_json(run_name, CH_params, grid_init, grid_level)
    character(*) :: run_name
    real(kind=REAL64), intent(out) :: CH_params(6)
    character(1), intent(out) :: grid_init
    integer, intent(out) :: grid_level

    grid_init = "r"
    grid_level = 2
    CH_params = (/1.0_REAL64, 1.0_REAL64, 1.0_REAL64, 1.0_REAL64, 1.0_REAL64, 1.0_REAL64/)
  end subroutine


end module