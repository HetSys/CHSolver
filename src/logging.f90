module logging
  use iso_fortran_env
  use logger_mod, only: flogger => master_logger, logger_init, debug, trivia, info, &
                                    warning, error, fatal

  implicit none


  type :: logger_type
    logical, private :: log_enabled

    contains
    procedure :: init => log_init
    procedure :: fatal => log_fatal
    procedure :: error => log_error
    procedure :: warning => log_warn
    procedure :: info => log_info
    procedure :: trivia => log_trivia
    procedure :: debug => log_debug
  end type logger_type

  character(*), parameter :: logfile_prefix = "CH"
  character(*), parameter :: logfolder = "logs/"


  type(logger_type) :: logger

  contains

  subroutine log_init(this, err, out, file)
    class(logger_type) :: this
    integer, intent(in) :: err, out, file
    character(100) :: cmd

    this%log_enabled = .TRUE.

    if (COMMAND_ARGUMENT_COUNT()==1)then
      call GET_COMMAND_ARGUMENT(1, cmd)
      if (cmd=="-nolog" .OR. cmd=="-NOLOG") this%log_enabled = .FALSE.
    end if

    if(this%log_enabled) call init_logfile(err, out, file)

  end subroutine

subroutine init_logfile(err, out, file)
  integer :: err, out, file
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


  ! Initialise master logger_type
  call logger_init(trim(logname), err, out, file)
end subroutine
  
  subroutine log_fatal(this, source, msg)
    class(logger_type), intent(in) :: this
    character(*), intent(in) :: source, msg
    
    if (this%log_enabled) then
      call flogger%fatal(source, msg)
    else
      write(ERROR_UNIT, *) "["//source//"]"//"<fatal> "//msg
    end if

  end subroutine

  subroutine log_error(this, source, msg)
    class(logger_type), intent(in) :: this
    character(*), intent(in) :: source, msg
    
    if (this%log_enabled) then
      call flogger%error(source, msg)
    else
      write(ERROR_UNIT, *) "["//source//"]"//"<error> "//msg
    end if
  end subroutine

  subroutine log_warn(this, source, msg)
    class(logger_type), intent(in) :: this
    character(*), intent(in) :: source, msg
    
    if (this%log_enabled) then
      call flogger%warning(source, msg)
    else
      write(ERROR_UNIT, *) "["//source//"]"//"<warning> "//msg
    end if
  end subroutine

  subroutine log_info(this, source, msg)
    class(logger_type), intent(in) :: this
    character(*), intent(in) :: source, msg
    
    if (this%log_enabled) call flogger%info(source, msg)
  end subroutine

  subroutine log_trivia(this, source, msg)
    class(logger_type), intent(in) :: this
    character(*), intent(in) :: source, msg
    
    if (this%log_enabled) call flogger%trivia(source, msg)
  end subroutine

  subroutine log_debug(this, source, msg)
    class(logger_type), intent(in) :: this
    character(*), intent(in) :: source, msg
    
    if (this%log_enabled) call flogger%debug(source, msg)
  end subroutine

end module logging