module globals
  use iso_fortran_env
  use logger_mod, logger => master_logger

  implicit none

  ! Kind/precision of reals to use
  integer, parameter :: real_kind = real64


  ! Logging defaults

  integer, parameter :: stderr_threshold = error

  integer, parameter :: stdout_threshold = info

  integer, parameter :: logfile_threshold = trivia

  char(*), parameter :: logfile_prefix = "ch-log-"

  integer(8) :: dt

  char(128) :: logname

  call date_and_time(values=dt)


  ! Log filename = "<logfile_prefix>yyyy-mm-dd-hh-mm-ss.log
  logname = logfile_prefix + dt(1) + "-" + dt(2) + "-" + dt(3) + "-" &
              dt(5) + "-" + dt(6) + "-" + dt(7) + ".log"


  ! Initialise master logger
  call logger_init(logname, , stderr_threshold, stdout_threshold, &
                    logfile_threshold)




end module globals