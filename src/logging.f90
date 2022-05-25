!> @brief Wrapper module to flogging library
!! @details Wrapping allows logging to be disabled and/or simplified to not include logfile
module logging
  use iso_fortran_env
  use logger_mod, only: flogger => master_logger, logger_init, debug, trivia, info, &
                                    warning, error, fatal

  implicit none

  !> @brief logger class - Wrapper to flogging master_logger
  !! @details Provides an interface to flogging logfile, stdout and stderr messages.
  !! Class Methods equivalent to flogging master_logger
  type :: logger_type
    logical :: log_enabled
    logical :: disable_all_logging
    character(128), private :: logdir

    contains
    procedure :: init => log_init
    procedure :: fatal => log_fatal
    procedure :: error => log_error
    procedure :: warning => log_warn
    procedure :: info => log_info
    procedure :: trivia => log_trivia
    procedure :: debug => log_debug
    procedure, private :: init_logfile
  end type logger_type

  character(*), parameter :: logfile_prefix = "CH"
  character(*), parameter :: logfolder = "logs"


  type(logger_type) :: logger

  contains

  !> @brief Initialise logger object
  !! @param this logger object
  !! @param[in] err Threshold for stderr messages
  !! @param[in] out Threshold for stdout messages
  !! @param[in] file Threshold for logfile messages
  !! @param[in] logdir Directory to generate logfile
  subroutine log_init(this, err, out, file, logdir)
    class(logger_type) :: this
    integer, intent(in) :: err, out, file
    character(*), intent(in), optional :: logdir

    if (present(logdir)) then
      this%logdir = logdir
    else
      this%logdir = logfolder
    end if

    if(this%log_enabled) call this%init_logfile(err, out, file)

  end subroutine

  !> @brief Generate logfile with standardised name
  !! @param[in] err Threshold for stderr messages
  !! @param[in] out Threshold for stdout messages
  !! @param[in] file Threshold for logfile messages
  subroutine init_logfile(this, err, out, file)
    class(logger_type) :: this
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
    logname = trim(this%logdir) // "/" // logfile_prefix // "-" // year // "-" // month // "-" // &
                day // "-" // hr // ":" // min // ":" // sec // ".log"


    ! Initialise master logger_type
    call logger_init(trim(logname), err, out, file)
  end subroutine
  
  !> @brief Wrapper for fatal logging messages
  !! @param[in] source Source function of the message
  !! @param[in] msg Message to send to logger
  subroutine log_fatal(this, source, msg)
    class(logger_type), intent(in) :: this
    character(*), intent(in) :: source, msg
    if(.NOT. logger%disable_all_logging) then
      if (this%log_enabled) then
        call flogger%fatal(source, msg)
      else
        write(ERROR_UNIT, *) "["//source//"]"//"<fatal> "//msg
      end if
    end if

  end subroutine


  !> @brief Wrapper for error logging messages
  !! @param[in] source Source function of the message
  !! @param[in] msg Message to send to logger
  subroutine log_error(this, source, msg)
    class(logger_type), intent(in) :: this
    character(*), intent(in) :: source, msg
    if(.NOT. logger%disable_all_logging) then
      if (this%log_enabled) then
        call flogger%error(source, msg)
      else
        write(ERROR_UNIT, *) "["//source//"]"//"<error> "//msg
      end if
    end if
  end subroutine


  !> @brief Wrapper for warning logging messages
  !! @param[in] source Source function of the message
  !! @param[in] msg Message to send to logger
  subroutine log_warn(this, source, msg)
    class(logger_type), intent(in) :: this
    character(*), intent(in) :: source, msg
    if(.NOT. logger%disable_all_logging) then
      if (this%log_enabled) then
        call flogger%warning(source, msg)
      else
        write(ERROR_UNIT, *) "["//source//"]"//"<warning> "//msg
      end if
    end if
  end subroutine


  !> @brief Wrapper for info logging messages
  !! @param[in] source Source function of the message
  !! @param[in] msg Message to send to logger
  subroutine log_info(this, source, msg)
    class(logger_type), intent(in) :: this
    character(*), intent(in) :: source, msg
    
    if(.NOT. logger%disable_all_logging .AND. this%log_enabled) call flogger%info(source, msg)
  end subroutine


  !> @brief Wrapper for trivia logging messages
  !! @param[in] source Source function of the message
  !! @param[in] msg Message to send to logger
  subroutine log_trivia(this, source, msg)
    class(logger_type), intent(in) :: this
    character(*), intent(in) :: source, msg
    
    if(.NOT. logger%disable_all_logging .AND. this%log_enabled) call flogger%trivia(source, msg)
  end subroutine


  !> @brief Wrapper for debug logging messages
  !! @param[in] source Source function of the message
  !! @param[in] msg Message to send to logger
  subroutine log_debug(this, source, msg)
    class(logger_type), intent(in) :: this
    character(*), intent(in) :: source, msg
    
    if(.NOT. logger%disable_all_logging .AND. this%log_enabled) call flogger%debug(source, msg)
  end subroutine

end module logging