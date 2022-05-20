!> @brief Parser for input JSON formatted files
!! @details Methods to search a target JSON file for the required input parameters
!! @author Tom Rocke

module JSON_Parser

  use globals
  use json_module

  implicit none

  !> @brief interface for retrieving parameterss from the JSON file
  interface json_retrieve
    module procedure json_retrieve_real
    module procedure json_retrieve_char
    module procedure json_retrieve_int 
    module procedure json_retrieve_realarr
  end interface json_retrieve
  contains

  !> @brief Reads JSON file, and searches for required input parameters
  !! @param[in]  run_name  String key for input values
  !! @param[in]  CH_params [L, A, M, K, p0, p1]
  !! @param[in]  grid_init Grid initialisation type character
  !! @param[in]  grid_level Controls size of grid
  !! @param[in]  Tout Output timestep array
  subroutine read_json(fname, run_name, CH_params, grid_init, grid_level, Tout)
    character(*), intent(in) :: fname, run_name
    real(kind=dp), intent(out) :: CH_params(6)
    character(:), allocatable, intent(out) :: grid_init
    real(dp), allocatable, intent(out), optional :: Tout(:)
    integer, intent(out) :: grid_level
    logical :: error

    type(json_file) :: json

    error = .FALSE.
    
    call open_json(json, fname)

    if (present(Tout)) then
      call get_json_params(json, run_name, CH_params, grid_init, grid_level, Tout, error)
    else
      call get_json_params(json, run_name, CH_params, grid_init, grid_level, error=error)
    end if
    
    if (error) then
      call logger%fatal("read_json", "Issues found fetching JSON params")
      stop
    end if

    call json%destroy()

    if (json%failed()) then
      call logger%fatal("open_json", "Failed to cleanup json object")
      stop
    end if
  end subroutine


  !> @brief Opens JSON file
  !! @param[inout] json json-fortran type(json_file) object to handle the interface
  !! @param[in] fname Filename of the JSON file to be opened
  subroutine open_json(json, fname)
    type(json_file), intent(inout) :: json
    character(*), intent(in) :: fname
    character(:), allocatable :: j_string

    ! Initialise JSON, allowing for path mode
    ! Nested keys accessed via $["outer key"]["inner key"]
    call json%initialize(path_mode=3)

    call logger%info("open_json", ("Reading "//fname))
    call json%load(filename=fname)

    if (json%failed()) then
      call logger%fatal("open_json", fname // " could not be opened.")
      stop
    end if
    
    call json%print_to_string(j_string)
    call logger%debug("open_json", "Found:")
    call logger%debug("open_json", j_string)
  end subroutine

  !> @brief Grab required input parameters from the JSON file open in the json object
  !! @param[inout] json json-fortran interface to the JSON file
  !! @param[in] run_name Name of the run to search for.
  !! @param[out] CH_params Will be populated by the CH parameters (L, A, M, K, p0, p1) 
  !! found in the JSON file
  !! @param[out] grid_init Will be populated with the initialisation character found in the JSON file
  !! @param[out] grid_level Will be populated with the grid level found in the JSON file
  !! @param[out] Tout Will be populated with the time array found in the JSON file
  !! @param[out] error Flag for whether errors occurred (parameters missing, etc)
  subroutine get_json_params(json, run_name, CH_params, grid_init, grid_level, Tout, error)
    type(json_file), intent(inout) :: json
    character(*), intent(in) :: run_name
    real(kind=dp), intent(out) :: CH_params(6)
    character(:), allocatable, intent(out) :: grid_init
    real(dp), allocatable, intent(out), optional :: Tout(:)
    integer, intent(out) :: grid_level

    logical :: found, all_found, error

    real(dp) :: L, A, M, K, p0, p1

    character(100) :: val_string

    error = .FALSE.
    ! Verify that run_name exists
    call logger%trivia("get_json_params","Validating run name")

    found = json%valid_path("$['"//run_name//"']")

    if (.NOT. found) then
      call logger%fatal("get_json_params", "Run name '"//run_name//"' not found")
      error = .TRUE.
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


    ! Search for init
    found = json_retrieve(json, run_name, "grid_type", grid_init)

    if (found) then
      call logger%trivia("get_json_params", "Found grid type of '"// trim(grid_init) // "'")
    else
      call logger%error("get_json_params", "Input Parameter 'grid_type' not found")
      all_found = .FALSE.
    end if

    ! Search for T
    found = .False.
    if (present(Tout)) found = json_retrieve(json, run_name, "T", Tout)

    if (found) then
      call logger%trivia("get_json_params", "Found T array of '"// trim(to_string(Tout)) // "'")
    else
      call logger%warning("get_json_params", "Input Parameter 'grid_type' not found in JSON file.")
    end if


    if (.NOT. all_found) then
      call logger%fatal("get_json_params", "Missing Input Parameters")
      error = .TRUE.
    end if

    ! Fill in CH_params array
    CH_params(1) = L
    CH_params(2) = A
    CH_params(3) = M
    CH_params(4) = K
    CH_params(5) = p0
    CH_params(6) = p1
  end subroutine

  !> @brief String formatter for JSON paths.
  !! @details Converts nested keys run_name and key_name into a 
  !! JSON searchable path from the root of the JSON data structure
  !! @param[in] run_name outermost key in search path
  !! @param[in] key_name nested key to be searched
  !! @result path Path from the JSON root to the value requested by the 
  !! run_name and key_name
  function json_get_path(run_name, key_name) result(path)
    character(*), intent(in) :: run_name, key_name
    character(128) :: path

    path = "$['" // run_name // "']['" // key_name // "']"
    call logger%debug("json_get_path", "Searching path "//trim(path))
  end function


  !> @brief Function to grab a value given by the path from run_name and key_name
  !! @details Searchs the path root -> run_name -> key_name in the JSON data structure held
  !! in the json object. Overwrites val and sets found to true if found.
  !! @param[inout] json JSON data structure object
  !! @param[in] run_name Outermost key in search path
  !! @param[in] key_name Nested key to be searched
  !! @param[out] val value found along path
  !! @result found flag for whether a value was found
  function json_retrieve_real(json, run_name, key_name, val) result(found)
    type(json_file), intent(inout) :: json
    character(*), intent(in) :: run_name, key_name
    real(dp), intent(out) :: val
    logical :: found

    character(128) :: path

    path = json_get_path(run_name, key_name)

    call json%get(trim(path), val, found)
  end function

  !> @brief Function to grab a value given by the path from run_name and key_name
  !! @details Searchs the path root -> run_name -> key_name in the JSON data structure held
  !! in the json object. Overwrites val and sets found to true if found.
  !! @param[inout] json JSON data structure object
  !! @param[in] run_name Outermost key in search path
  !! @param[in] key_name Nested key to be searched
  !! @param[out] val value found along path
  !! @result found flag for whether a value was found
  function json_retrieve_char(json, run_name, key_name, val) result(found)
    type(json_file), intent(inout) :: json
    character(*), intent(in) :: run_name, key_name
    character(len=:), allocatable, intent(inout) :: val
    logical :: found

    character(128) :: path

    path = json_get_path(run_name, key_name)

    call json%get(trim(path), val, found)
  end function

  !> @brief Function to grab a value given by the path from run_name and key_name
  !! @details Searchs the path root -> run_name -> key_name in the JSON data structure held
  !! in the json object. Overwrites val and sets found to true if found.
  !! @param[inout] json JSON data structure object
  !! @param[in] run_name Outermost key in search path
  !! @param[in] key_name Nested key to be searched
  !! @param[out] val value found along path
  !! @result found flag for whether a value was found
  function json_retrieve_int(json, run_name, key_name, val) result(found)
    type(json_file), intent(inout) :: json
    character(*), intent(in) :: run_name, key_name
    integer, intent(out) :: val
    logical :: found

    character(128) :: path

    path = json_get_path(run_name, key_name)

    call json%get(trim(path), val, found)
  end function

  !> @brief Function to grab a value given by the path from run_name and key_name
  !! @details Searchs the path root -> run_name -> key_name in the JSON data structure held
  !! in the json object. Overwrites val and sets found to true if found.
  !! @param[inout] json JSON data structure object
  !! @param[in] run_name Outermost key in search path
  !! @param[in] key_name Nested key to be searched
  !! @param[out] val value found along path
  !! @result found flag for whether a value was found
  function json_retrieve_realarr(json, run_name, key_name, val) result(found)
    type(json_file), intent(inout) :: json
    character(*), intent(in) :: run_name, key_name
    real(dp), allocatable, intent(out) :: val(:)
    logical :: found

    character(128) :: path

    path = json_get_path(run_name, key_name)

    if (allocated(val)) deallocate(val)

    call json%get(trim(path), val, found)
  end function

end module
