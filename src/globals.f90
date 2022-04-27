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

  character(*), parameter :: logfile_prefix = "ch-log"
  character(*), parameter :: logfolder = "logs/"

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
  !! @todo input validation 
  subroutine validate_params(CH_params, grid_init, grid_level)
    real(kind=dp), intent(in) :: CH_params(6)
    character(1), intent(in) :: grid_init
    integer, intent(in) :: grid_level

  end subroutine

end module globals
