!> @brief Globally accessible methods and variables
!! @details Includes methods for string conversion, logging setup,
!! parameter validation, and random seeding setup

module globals
  use iso_fortran_env
  use iso_c_binding
  use logging, only: logger, log_debug => debug, log_trivia => trivia, log_info => info, &
                     log_warning => warning, log_error => error, log_fatal => fatal

  implicit none
  !> @var integer dp
  !! Kind/precision of reals to use
  integer, parameter :: dp = real64
  integer, parameter :: cdc = C_DOUBLE_COMPLEX

  ! Logging defaults

  !> @var integer stderr_threshold
  !! Default threshold for logging to stderr
  integer, parameter :: stderr_threshold = log_error

  !> @var integer stdout_threshold
  !! Default threshold for logging to stdout
  integer, parameter :: stdout_threshold = log_info

  !> @var integer logfile_threshold
  !! Default threshold for logging to logfile
  integer, parameter :: logfile_threshold = log_trivia ! set to debug for more info in the logfile


  !> @brief Interface for converting various types into string representations
  interface to_string
    module procedure real_to_string
    module procedure int_to_string
    module procedure realarr_to_string
  end interface to_string


  contains

  !> @brief Validation of program inputs
  !! @details Performs a series of validation tests (eg non-negative tests)
  !! to ensure that the input parameters are viable
  !! @param[in]  CH_params [L, A, M, K, p0, p1]
  !! @param[in]  grid_init Grid initialisation type character
  !! @param[in]  grid_level Controls size of grid
  !! @param[in]  Tout Output timesteps
  !! @param[out] errors Returns true if any errors have been detected
  subroutine validate_params(CH_params, grid_init, grid_level, Tout, errors)
    real(kind=dp), intent(in) :: CH_params(6)
    character(*), intent(in) :: grid_init
    integer, intent(in) :: grid_level
    real(dp), allocatable, intent(in), optional :: Tout(:)
    logical, intent(out) :: errors
    real(dp) :: L, A, M, K, p0, p1

    real(dp), parameter :: REAL_TOL = 1e-10_DP

    character(*), parameter :: ACCEPTED_INITS = "rcbs"


    errors = .FALSE.

    L = CH_params(1)
    A = CH_params(2)
    M = CH_params(3)
    K = CH_params(4)
    p0 = CH_params(5)
    p1 = CH_params(6)


    ! Grid init char validation
    if (len(grid_init) /=1) then
      call logger%error("validate_params", "Grid init character should be of length 1, got "&
                  // "length "//adjustl(trim(to_string(len(grid_init)))))
      errors = .TRUE.
    end if

    if (verify(grid_init, ACCEPTED_INITS) /= 0) then
      ! Verify() returns 0 if all chars of grid_init are in ACCEPTED_INITS
      call logger%error("validate_params", "Grid init char '"//grid_init//"' not in accepted "&
                  // "list of init chars: "//ACCEPTED_INITS)
      errors = .TRUE.
    end if

    ! Grid level validation
    if (grid_level<0) then
      call logger%error("validate_params", "Grid level cannot be negative")
      call logger%error("validate_params", "Got Grid level ="//adjustl(trim(to_string(grid_level))))
      errors = .TRUE.
    end if

    ! L param validation
    if (L<0) then
      call logger%error("validate_params", "L cannot be negative")
      call logger%error("validate_params", "Got L="//adjustl(trim(to_string(L))))
      errors = .TRUE.
    end if

    ! A param validation
    if (A<=0) then
      call logger%error("validate_params", "A cannot be negative")
      call logger%error("validate_params", "Got A="//adjustl(trim(to_string(A))))
      errors = .TRUE.
    end if

    ! M param validation
    if (M<0) then
      call logger%error("validate_params", "M cannot be negative")
      call logger%error("validate_params", "Got M="//adjustl(trim(to_string(M))))
      errors = .TRUE.
    end if
    if (abs(M) < REAL_TOL) then
      call logger%warning("validate_params", "M=0 implies system will not evolve")
      call logger%warning("validate_params", "abs("//adjustl(trim(to_string(M)))&
                //") was less "//"than tolerance of "//adjustl(trim(to_string(REAL_TOL))))
    end if

    ! K param validation
    if (K<0) then
      call logger%error("validate_params", "K cannot be negative")
      call logger%error("validate_params", "Got K="//adjustl(trim(to_string(K))))
      errors = .TRUE.
    end if

    ! p0 param validation


    ! p1 param validation

    if (abs(p0 - p1) < REAL_TOL) then
      call logger%error("validate_params", "p0 and p1 cannot be equal")

      call logger%error("validate_params", "Got p0="//adjustl(trim(to_string(p0)))//" and p1="&
                      //adjustl(trim(to_string(p1))))

      call logger%error("validate_params", "abs("//adjustl(trim(to_string(p0)))//" - "//adjustl(trim(to_string(p1)))&
                        //") was less "//"than tolerance of "//adjustl(trim(to_string(REAL_TOL))))
      errors = .TRUE.
    end if

    ! Tout validation
    if (present(Tout)) then
      if (.not. allocated(Tout)) then
        call logger%error("validate_params", "Output times not specified")
        errors = .TRUE.
      else if (any(Tout .LT. 0)) then
        call logger%error("validate_params", "Cannot output at negative timesteps")
        errors = .TRUE.
      else if (any(Tout(2:) .LE. Tout(:size(Tout)))) then
        call logger%error("validate_params", "Output times must be strictly increasing")
        errors = .TRUE.
      end if
    end if
  end subroutine


  !> @brief Converter from reals to strings
  !! @param[in] val Real number to be converted
  !! @result str string representation of val
  function real_to_string(val) result(str)
    real(dp), intent(in) :: val
    character(128):: str

    write(str, *) val

  end function real_to_string

  !> @brief Converter from integers to strings
  !! @param[in] val Integer number to be converted
  !! @result str string representation of val
  function int_to_string(val) result(str)
    integer, intent(in) :: val
    character(128) :: str

    write(str, *) val

  end function int_to_string

  !> @brief Converter from an array of reals to a string
  !! @param[in] val Real array to be converted
  !! @result str string representation of val
  function realarr_to_string(val) result(str)
    real(dp), dimension(:), intent(in) :: val
    character(2048) :: str
    write(str, *) val
    str = "{" // trim(str) // "}"
  end function realarr_to_string


  !> @brief Initializes Seed for PRNG
  !! @param[in]  seed Integer seed for PRNG.
  !! If Seed is set to -1, will use datetime for seed
  subroutine set_ranseed(seed)
    integer, intent(in) :: seed

    integer, allocatable :: state(:)
    integer :: state_size
    integer, dimension(8) :: datetime

    call random_seed(size=state_size)
    allocate(state(state_size))

    if (seed == -1) then

      !Current milisecond of the hour
      call date_and_time(values=datetime)
      state = 60000*datetime(6) + 1000*datetime(7) + datetime(8)

    else
      state = seed
    end if

    call random_seed(put=state)

  end subroutine

end module globals
