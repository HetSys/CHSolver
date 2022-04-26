module globals
  use iso_fortran_env
  use logger_mod, only: logger => master_logger, logger_init, debug, trivia, info, &
                                    warning, error, fatal

  implicit none

  ! Kind/precision of reals to use
  integer, parameter :: real_kind = real64


  ! Logging defaults

  integer, parameter :: stderr_threshold = error

  integer, parameter :: stdout_threshold = info

  integer, parameter :: logfile_threshold = trivia

  character(*), parameter :: logfile_prefix = "ch-log"

  contains


  subroutine initialise()
    character(8) :: date 
    character(10) :: time

    character(128) :: logname


    call date_and_time(date=date, time=time)


    ! Log filename = "<logfile_prefix>yyyy-mm-dd-hh-mm-ss.log
    logname = logfile_prefix // "-" // date // "-" // time(:7) // ".log"


    ! Initialise master logger
    call logger_init(trim(logname), stderr_threshold, stdout_threshold, &
                      logfile_threshold)
  end subroutine initialise




end module globals
