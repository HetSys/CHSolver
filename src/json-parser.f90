module JSON_Parser
  !> Parses JSON file for given run name
  !! Returns the params needed for grid generation, and solving

  use globals
  use json_module

  implicit none
  

  character(*), parameter :: JSON_FILENAME = "./input-data.json"



  interface json_retrieve
    module procedure json_retrieve_real
    module procedure json_retrieve_char
    module procedure json_retrieve_int 
  end interface json_retrieve
  contains

  !> @brief Reads JSON file, and searches for given params
  !! @param[in]  run_name  String key for input values
  !! @param[in]  CH_params [L, A, M, K, p0, p1]
  !! @param[in]  grid_init Grid initialisation type character
  !! @param[in]  grid_level Controls size of grid
  subroutine read_json(run_name, CH_params, grid_init, grid_level)
    character(*), intent(in) :: run_name
    real(kind=dp), intent(out) :: CH_params(6)
    character(:), allocatable, intent(out) :: grid_init
    integer, intent(out) :: grid_level
    logical :: error

    type(json_file) :: json

    call open_json(json)

    call get_json_params(json, run_name, CH_params, grid_init, grid_level)

    call validate_params(CH_params, grid_init, grid_level, error)

    if (error) then
      call logger%fatal("read_json", "Issues found with input parameters")
      stop
    end if

    call logger%trivia("read_json", "No issues found in input parameters")

    call json%destroy()
  end subroutine


  !> @brief Opens JSON file
  subroutine open_json(json)
    type(json_file), intent(inout) :: json
    character(:), allocatable :: j_string

    ! Initialise JSON, allowing for path mode
    ! Nested keys accessed via $["outer key"]["inner key"]
    call json%initialize(path_mode=3)

    call logger%info("open_json", ("Reading "//JSON_FILENAME))
    call json%load(filename=JSON_FILENAME)

    if (json%failed()) then
      call logger%fatal("open_json", JSON_FILENAME // " could not be opened.")
      stop
    end if
    
    call json%print_to_string(j_string)
    call logger%debug("open_json", "Found:")
    call logger%debug("open_json", j_string)
  end subroutine


  !> @brief Grab required params from the JSON file opened by json_handler
  subroutine get_json_params(json, run_name, CH_params, grid_init, grid_level)
    type(json_file), intent(inout) ::json
    character(*), intent(in) :: run_name
    real(kind=dp), intent(out) :: CH_params(6)
    character(:), allocatable, intent(out) :: grid_init
    integer, intent(out) :: grid_level

    logical :: found, all_found

    real(dp) :: L, A, M, K, p0, p1

    character(100) :: val_string

    ! Verify that run_name exists
    call logger%trivia("get_json_params","Validating run name")

    found = json%valid_path("$['"//run_name//"']")

    if (.NOT. found) then
      call logger%fatal("get_json_params", "Run name '"//run_name//"' not found")
      stop
    end if


    ! Search run_name for input params

    all_found = .TRUE.

    ! Search for L
    found = json_retrieve(json, run_name, "L", L)

    if (found) then
      write(val_string, *) L
      call logger%trivia("get_json_params", "Found L value of "// trim(val_string))
    else
      call logger%error("get_json_params", "Input Parameter L not found")
      all_found = .FALSE.
    end if


    ! Search for A
    found = json_retrieve(json, run_name, "A", A)

    if (found) then
      write(val_string, *) A
      call logger%trivia("get_json_params", "Found A value of "// trim(val_string))
    else
      call logger%error("get_json_params", "Input Parameter A not found")
      all_found = .FALSE.
    end if

    ! Search for M
    found = json_retrieve(json, run_name, "M", M)

    if (found) then
      write(val_string, *) M
      call logger%trivia("get_json_params", "Found M value of "// val_string)
    else
      call logger%error("get_json_params", "Input Parameter M not found")
      all_found = .FALSE.
    end if

    ! Search for K
    found = json_retrieve(json, run_name, "K", K)

    if (found) then
      write(val_string, *) K
      call logger%trivia("get_json_params", "Found K value of "// val_string)
    else
      call logger%error("get_json_params", "Input Parameter K not found")
      all_found = .FALSE.
    end if

    ! Search for p0
    found = json_retrieve(json, run_name, "p0", p0)

    if (found) then
      write(val_string, *) p0
      call logger%trivia("get_json_params", "Found p0 value of "// val_string)
    else
      call logger%error("get_json_params", "Input Parameter p0 not found")
      all_found = .FALSE.
    end if

    ! Search for p1
    found = json_retrieve(json, run_name, "p1", p1)

    if (found) then
      write(val_string, *) p1
      call logger%trivia("get_json_params", "Found p1 value of "// trim(val_string))
    else
      call logger%error("get_json_params", "Input Parameter p1 not found")
      all_found = .FALSE.
    end if

    ! Search for level
    found = json_retrieve(json, run_name, "grid_level", grid_level)

    if (found) then
      write(val_string, "(i4)") grid_level
      call logger%trivia("get_json_params", "Found level value of "// trim(val_string))
    else
      call logger%error("get_json_params", "Input Parameter level not found")
      all_found = .FALSE.
    end if


    ! Search for level
    found = json_retrieve(json, run_name, "grid_type", grid_init)

    if (found) then
      call logger%trivia("get_json_params", "Found grid type of '"// trim(grid_init) // "'")
    else
      call logger%error("get_json_params", "Input Parameter 'grid_type' not found")
      all_found = .FALSE.
    end if


    if (.NOT. all_found) then
      call logger%fatal("get_json_params", "Missing Input Parameters")
      stop
    end if

    ! Fill in CH_params array
    CH_params(1) = L
    CH_params(2) = A
    CH_params(3) = M
    CH_params(4) = K
    CH_params(5) = p0
    CH_params(6) = p1
  end subroutine

  !> @brief String formatter for JSON paths
  function json_get_path(run_name, key_name) result(path)
    character(*), intent(in) :: run_name, key_name
    character(128) :: path

    path = "$['" // run_name // "']['" // key_name // "']"
    call logger%debug("json_get_path", "Searching path "//trim(path))
  end function

  function json_retrieve_real(json, run_name, key_name, val) result(found)
    type(json_file), intent(inout) :: json
    character(*), intent(in) :: run_name, key_name
    real(dp), intent(out) :: val
    logical :: found

    character(128) :: path

    path = json_get_path(run_name, key_name)

    call json%get(trim(path), val, found)
  end function

  function json_retrieve_char(json, run_name, key_name, val) result(found)
    type(json_file), intent(inout) :: json
    character(*), intent(in) :: run_name, key_name
    character(len=:), allocatable, intent(inout) :: val
    logical :: found

    character(128) :: path

    path = json_get_path(run_name, key_name)

    call json%get(trim(path), val, found)
  end function

  function json_retrieve_int(json, run_name, key_name, val) result(found)
    type(json_file), intent(inout) :: json
    character(*), intent(in) :: run_name, key_name
    integer, intent(out) :: val
    logical :: found

    character(128) :: path

    path = json_get_path(run_name, key_name)

    call json%get(trim(path), val, found)
  end function

end module