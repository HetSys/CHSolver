module globals
  use iso_fortran_env
  use logger_mod, only: logger => master_logger, logger_init, debug, trivia, info, &
                                    warning, error, fatal

  implicit none

  ! Kind/precision of reals to use
  integer, parameter :: dp = real64


  ! Logging defaults

  integer, parameter :: stderr_threshold = error

  integer, parameter :: stdout_threshold = info

  integer, parameter :: logfile_threshold = trivia ! set to debug for more info in the logfile

  character(*), parameter :: logfile_prefix = "CH"
  character(*), parameter :: logfolder = "logs/"

  interface to_string
    module procedure real_to_string
    module procedure int_to_string
  end interface to_string


  contains

  !> @Brief Initialise logging
  !! Wrapper for the flogging logger_init() subroutine
  !! Automatically handles logfile creation and naming
  subroutine initialise()
    character(8) :: date 
    character(10) :: time

    character(128) :: logname

    character(4) :: year
    character(2) :: month, day, hr, min, sec

    call date_and_time(date=date, time=time)

    year = date(:5)
    month = date(5:7)
    day = date(7:8)

    hr = time(1:3)
    min = time(3:5)
    sec = time(5:7)

    ! Log filename = "<logfile_prefix>yyyy-mm-dd-hh-mm-ss.log
    logname = logfolder // logfile_prefix // "-" // year // "-" // month // "-" // &
                day // "-" // hr // ":" // min // ":" // sec // ".log"


    ! Initialise master logger
    call logger_init(trim(logname), stderr_threshold, stdout_threshold, &
                      logfile_threshold)
  end subroutine initialise


  !> @brief Validation of program inputs
  !! @param[in]  CH_params [L, A, M, K, p0, p1]
  !! @param[in]  grid_init Grid initialisation type character
  !! @param[in]  grid_level Controls size of grid
  !! @todo grid_level validation
  subroutine validate_params(CH_params, grid_init, grid_level)
    real(kind=dp), intent(in) :: CH_params(6)
    character(*), intent(in) :: grid_init
    integer, intent(in) :: grid_level

    real(dp) :: L, A, M, K, p0, p1

    logical :: errors

    real(dp), parameter :: REAL_TOL = 1e-10_DP

    character(*), parameter :: ACCEPTED_INITS = "r"


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


    ! Check for any errors

    if (errors) then
      call logger%fatal("validate_params", "Issues found with input parameters")
      stop
    end if

    call logger%trivia("validate_params", "No issues found in input parameters")
  end subroutine



  function real_to_string(val) result(str)
    real(dp), intent(in) :: val
    character(128) :: str

    write(str, *) val

  end function real_to_string

  function int_to_string(val) result(str)
    integer, intent(in) :: val
    character(128) :: str

    write(str, *) val

  end function int_to_string

end module globals
