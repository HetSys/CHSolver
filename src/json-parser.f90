module JSON_Parser
  !> Parses JSON file for given run name
  !! Returns the params needed for grid generation, and solving

  use globals
  use json_module
  implicit none
  

  character(*), parameter :: JSON_FILENAME = "input-data.json"



  interface json_retrieve
    module procedure json_retrieve_real
    module procedure json_retrieve_char
    module procedure json_retrieve_int 
  end interface json_retrieve
  contains


  !> @brief Validation of program inputs
  !! @param[in]  CH_params [L, A, M, K, p0, p1]
  !! @param[in]  grid_init Grid initialisation type character
  !! @param[in]  grid_level Controls size of grid
  subroutine validate_json_params(CH_params, grid_init, grid_level)
    real(kind=dp), intent(in) :: CH_params(6)
    character(1), intent(in) :: grid_init
    integer, intent(in) :: grid_level

  end subroutine

  !> @brief Reads JSON file, and searches for given params
  !! @param[in]  run_name  String key for input values
  !! @param[in]  CH_params [L, A, M, K, p0, p1]
  !! @param[in]  grid_init Grid initialisation type character
  !! @param[in]  grid_level Controls size of grid
  subroutine read_json(run_name, CH_params, grid_init, grid_level)
    character(*), intent(in) :: run_name
    real(kind=dp), intent(out) :: CH_params(6)
    character(1), intent(out) :: grid_init
    integer, intent(out) :: grid_level

    grid_init = "r"
    grid_level = 2
    CH_params = (/1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp/)
  end subroutine


  !> @brief Opens JSON file
  function open_json() result(json)
    type(json_file) :: json

    ! Initialise JSON, allowing for path mode
    ! Nested keys accessed via $["outer key"]["inner key"]
    call json%initialize(path_mode=3)

    call json%load(filename=JSON_FILENAME)
    
    call logger%trivia("open_json", ("Opening "//JSON_FILENAME))
  end function


  !> @brief Grab required params from the JSON file opened by json_handler
  subroutine get_json_params(json, run_name, CH_params, grid_init, grid_level)
    type(json_file), intent(in) ::json
    character(*), intent(in) :: run_name
    real(kind=dp), intent(out) :: CH_params(6)
    character(1), intent(out) :: grid_init
    integer, intent(out) :: grid_level

    logical :: found, all_found

    real(dp) :: L, A, M, K, p0, p1

    character(7) val_string

    all_found = .TRUE.

    ! Search for L
    found = json_retrieve(json, run_name, "L", L)

    if (found) then
      write(val_string, "(f3.3)") L
      call logger%trivia("get_json_params", "Found L value of "// val_string)
    else
      call logger%error("get_json_params", "Input Parameter L not found")
      all_found = .FALSE.
    end if


    ! Search for A
    found = json_retrieve(json, run_name, "A", A)

    if (found) then
      write(val_string, "(f3.3)") A
      call logger%trivia("get_json_params", "Found LAvalue of "// val_string)
    else
      call logger%error("get_json_params", "Input Parameter A not found")
      all_found = .FALSE.
    end if

    ! Search for M
    found = json_retrieve(json, run_name, "M", M)

    if (found) then
      write(val_string, "(f3.3)") M
      call logger%trivia("get_json_params", "Found M value of "// val_string)
    else
      call logger%error("get_json_params", "Input Parameter M not found")
      all_found = .FALSE.
    end if

    ! Search for K
    found = json_retrieve(json, run_name, "K", K)

    if (found) then
      write(val_string, "(f3.3)") K
      call logger%trivia("get_json_params", "Found K value of "// val_string)
    else
      call logger%error("get_json_params", "Input Parameter K not found")
      all_found = .FALSE.
    end if

    ! Search for p0
    found = json_retrieve(json, run_name, "p0", p0)

    if (found) then
      write(val_string, "(f3.3)") p0
      call logger%trivia("get_json_params", "Found p0 value of "// val_string)
    else
      call logger%error("get_json_params", "Input Parameter p0 not found")
      all_found = .FALSE.
    end if

    ! Search for p1
    found = json_retrieve(json, run_name, "p1", p1)

    if (found) then
      write(val_string, "(f3.3)") p1
      call logger%trivia("get_json_params", "Found p1 value of "// val_string)
    else
      call logger%error("get_json_params", "Input Parameter p1 not found")
      all_found = .FALSE.
    end if

    ! Search for level
    found = json_retrieve(json, run_name, "grid_level", grid_level)

    if (found) then
      write(val_string, "(i6)") grid_level
      call logger%trivia("get_json_params", "Found level value of "// trim(val_string))
    else
      call logger%error("get_json_params", "Input Parameter level not found")
      all_found = .FALSE.
    end if


    ! Search for level
    found = json_retrieve(json, run_name, "grid_type", grid_init)

    if (found) then
      write(val_string, "(i6)") grid_level
      call logger%trivia("get_json_params", "Found grid type of "// trim(grid_init))
    else
      call logger%error("get_json_params", "Input Parameter 'grid_type' not found")
      all_found = .FALSE.
    end if
  end subroutine

  !> @brief String formatter for JSON paths
  function json_get_path(run_name, key_name) result(path)
    character(*), intent(in) :: run_name, key_name
    character(128) :: path

    path = "$['" // run_name // "']['" // key_name // "']"
  end function

  function json_retrieve_real(json, run_name, key_name, val) result(found)
    type(json_file), intent(in) ::json
    character(*), intent(in) :: run_name, key_name
    real(dp), intent(out) :: val
    logical :: found

    character(128) :: path

    path = json_get_path(run_name, key_name)

    call json%get(trim(path), val, found)
  end function

  function json_retrieve_char(json, run_name, key_name, val) result(found)
    type(json_file), intent(in) ::json
    character(*), intent(in) :: run_name, key_name
    character(1), intent(out) :: val
    logical :: found

    character(128) :: path

    path = json_get_path(run_name, key_name)

    call json%get(trim(path), val, found)
  end function

  function json_retrieve_int(json, run_name, key_name, val) result(found)
    type(json_file), intent(in) ::json
    character(*), intent(in) :: run_name, key_name
    integer, intent(out) :: val
    logical :: found

    character(128) :: path

    path = json_get_path(run_name, key_name)

    call json%get(trim(path), val, found)
  end function

end module