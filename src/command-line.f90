module command_line
  use globals
  implicit none

  real(dp), private :: cmd_CH_params(6)
  character(1), private :: cmd_grid_init
  integer, private :: cmd_level
  character, allocatable, private :: cmd_fname, cmd_runname, cmd_outpath, cmd_logdir
  real(dp), allocatable, private :: cmd_timearray(:)

  integer, private :: cmd_stdout_val

  logical, private :: CH_fnd(6), init_fnd, level_fnd
  contains

  !> @Brief Parse command line args, Initialise logging
  !! Wrapper for the flogging logger_init() subroutine
  !! Automatically handles logfile creation and naming
  subroutine initialise(err_threshold, out_threshold, file_threshold, ranseed)
    integer, intent(in), optional :: err_threshold, out_threshold, file_threshold, ranseed

    integer :: err, out, file, seed

    cmd_stdout_val = stdout_threshold

    logger%disable_all_logging = .FALSE.
    logger%log_enabled = .TRUE.

    CH_fnd = .FALSE.
    init_fnd = .FALSE.
    level_fnd = .FALSE.

    call parse_args()

    err = stderr_threshold
    out = cmd_stdout_val
    file = logfile_threshold
    seed = -1

    if (present(err_threshold)) err = err_threshold
    if (present(out_threshold)) out = out_threshold
    if (present(file_threshold)) file = file_threshold
    if (present(ranseed)) seed = ranseed

    if (allocated(cmd_logdir)) then
      call logger%init(err, out, file, cmd_logdir)
    else
      call logger%init(err, out, file)
    end if
  end subroutine initialise




  !! Grab JSON filename, run name, and directory to store outputs from parsed command line args
  !> Will only modify filename, run_name, and output_dir if
  !> relevant command line overrides were found
  subroutine get_io_commands(filename, run_name, output_dir)
    character, allocatable, optional, intent(inout) :: filename, run_name, output_dir

    if (present(filename) .AND. allocated(cmd_fname)) filename = cmd_fname

    if (present(run_name) .AND. allocated(cmd_runname)) run_name = cmd_runname

    if (present(output_dir) .AND. allocated(cmd_fname)) output_dir = cmd_outpath
  end subroutine


  !! Grab CH Params, grid level, grid init, and output times from parsed command line args
  !> Will only modify parts of CH_params, or level, init, time_arr
  !> if relevant overrides were found
  subroutine get_input_commands(CH_params, level, init, time_arr)
    real(dp), intent(inout), optional :: CH_params(6)
    integer, intent(inout), optional :: level
    character(*), intent(inout), optional :: init
    real(dp), allocatable, intent(inout), optional :: time_arr(:)

    integer :: idx

    if (present(CH_params)) then
      do idx=1,6
        if (CH_fnd(idx)) CH_params(idx) = cmd_CH_params(idx)
      end do
    end if

    if (present(level) .AND. level_fnd) level = cmd_level

    if (present(init) .AND. init_fnd) init = cmd_grid_init

    if (present(time_arr) .AND. allocated(cmd_timearray)) then
      if (allocated(time_arr)) deallocate(time_arr)

      allocate(time_arr(size(cmd_timearray)))
      time_arr = cmd_timearray
    end if
    
  end subroutine

  !! Parse command line args, modifying private variables when overriding values found from the command line
  subroutine parse_args()
    integer :: num_args
    character(:), allocatable :: key_arg, val_arg
    character(len=100) :: arg
    integer :: current_arg, len_arg, equals_pos, idx
    num_args = command_argument_count()

    cmd_CH_params = -1.0_dp

    !if (num_args == 0) return ! Nothing to do if no args provided
    current_arg = 1
    do current_arg=1, num_args
      call get_command_argument(current_arg, arg, len_arg)
      ! Ignore arg if it's not a -{key} or --{key}
      if (len_arg < 2 .OR. arg(1:1) /= "-") cycle
      ! Check if arg is short (-{key}), or long (--{key})
      if (arg(1:2) == "--") then
        ! Long arg (--{key}={val})
        equals_pos = index(arg, "=")
        if (equals_pos==0) equals_pos = len_arg + 1 ! No val_arg found
        key_arg = trim(arg(3:equals_pos - 1))
        arg = arg(equals_pos + 1 :)
        if (verify("{", arg) /= 0) then
          ! arg is the start of an array
          val_arg =  trim(parse_array_arg(current_arg, trim(arg)))
        else
          val_arg = trim(arg)
        end if
        call parse_keyval_arg(key_arg, val_arg)

      else if (arg(1:1) == "-") then
        ! Short arg
        key_arg = trim(arg(2:))
        call get_command_argument(current_arg+1, arg)
        if (verify("{", arg) == 0) then
          ! arg is the start of an array
          val_arg =  trim(parse_array_arg(current_arg+2, arg))
        else
          val_arg = trim(arg)
        end if
        do idx=1, len_arg - 1
          ! Loop through all chars in key_arg
          ! EG if -lamk 1.0 specified
          ! Should be equivalent to -l 1.0 -a 1.0...
          call  parse_keyval_arg(key_arg(idx:idx), val_arg)
        end do
      end if
    end do
  end subroutine
  
  !! Find the extent of an array argument, returning the full found array
  function parse_array_arg(current_arg, start_val) result(arr)
    integer, intent(in) :: current_arg
    character(*), intent(in) :: start_val
    character(256) :: arr, arg
    integer :: arg_len, idx
    arr = trim(start_val)
    do idx = current_arg, command_argument_count()
      call get_command_argument(idx, arg, arg_len)
      arr = trim(arr) // " " // trim(arg)
      if (verify("}", arg) == 0) then
        exit ! End of array found
      end if
    end do
  end function





  !! Parses given key arg and val arg
  !> Specified either by -{key} {val} or by 
  !> --{key}={val}
  subroutine parse_keyval_arg(key_arg, val_arg)
    character(*), intent(in) :: key_arg
    character(*), intent(in), optional :: val_arg

    integer :: idx

    select case (key_arg)
      ! HELP
      case ("h", "help")
        call print_help_text()
        stop
      ! LOGGING & STDOUTPUT
      case ("v", "verbose")
        cmd_stdout_val = cmd_stdout_val - 10
      case ("V", "version")
        print *, "CHSolver 1.0"
        stop
      case ("q", "quiet")
        logger%log_enabled = .FALSE.
        cmd_stdout_val = cmd_stdout_val + 10
      case ("s", "silent")
        logger%disable_all_logging = .TRUE.
        logger%log_enabled = .FALSE.
      case ("log_dir")
        cmd_logdir = val_arg
      ! JSON
      case ("j", "json_file")
        if (present(val_arg)) then
          cmd_fname = val_arg
        end if
      case ("r", "run_name")
        if (present(val_arg)) then
          cmd_runname = val_arg
        end if
      ! INPUT PARAM OVERRIDES
      case ("t", "time_array")
        if (present(val_arg)) then
          call allocate_t_array(val_arg)
        end if
      case ("l", "a", "m", "k", "0", "1", "p0", "p1")
        select case (key_arg)
        case ("l")
          idx = 1
        case ("a")
          idx = 2
        case ("m")
          idx = 3
        case ("k")
          idx = 4
        case ("0", "p0")
          idx = 5
        case ("1", "p1")
          idx = 6
        end select
        if (present(val_arg)) then
          cmd_CH_params(idx) = str_to_real(val_arg)
        end if
      case ("i", "init")
        if (present(val_arg)) then
          cmd_grid_init = val_arg
        end if
      case ("L", "level")
        if (present(val_arg)) then
          cmd_level = str_to_int(val_arg)
        end if
      case ("o", "out_dir")
        if (present(val_arg)) then
          cmd_outpath = val_arg
        end if
      case ("p")
      end select
  end subroutine parse_keyval_arg



  subroutine print_help_text()
    character(1), parameter :: newline = NEW_LINE('a')
    print *, "CHSolver Command Line Options", newline
    print *, "Usage:", newline, "  chsolver [options]", newline
    print *, "General Options:"
    print *, "-h, --help                         Show helper message for command line interface and exit."
    print *, "-V, --version                      Print the current version and exit.", newline
    print *, "Verbosity Options:"
    print *, "-v, --verbose                      Increase the amount of information printed to standard output."
    print *, "                                     Can be called twice to add debug messages."
    print *, "-q, --quiet                        Limits standard output messages to warning messages only."
    print *, "                                     Can be called twice to ignore warnings."
    print *, "-s, --silent                       Disable all messages to standard output and error."
    print *, "                                     WARNING: Will silently abort if errors are found.", newline
    print *, "File IO Options:"
    print *, "-j <file>, --json_file=<file>      Sets the path of the JSON input file, defaulting to './input-data.json."
    print *, "-r <name>, --run_name=<name>       Sets the run name to search for input parameters"
    print *, "-o <dir>, --out_dir=<dir>          Sets the path of the directory to save output trajectories to,"
    print *, "                                     defaulting to './output'"
    print *, "--log_dir=<dir>                    Sets the path o the directory to save log files to, defaulting"
    print *, "                                     to './logs'", newline
    print *, "Input Parameter Overrides:"
    print *, "-l <val>, --l=<val>                Sets the 'L' Cahn-Hilliard parameter to <val>"
    print *, "-a <val>, --a=<val>                Sets the 'A' Cahn-Hilliard parameter to <val>"
    print *, "-m <val>, --m=<val>                Sets the 'M' Cahn-Hilliard parameter to <val>"
    print *, "-k <val>, --k=<val>                Sets the 'K' Cahn-Hilliard parameter to <val>"
    print *, "-p0 <val>, -0 <val> , --p0=<val>   Sets the 'p0' Cahn-Hilliard parameter to <val>"
    print *, "-p1 <val>, -1 <val> , --p1=<val>   Sets the 'p1' Cahn-Hilliard parameter to <val>", newline
    print *, "-L <val> --level=<val>             Set the grid level to <val>."
    print *, "                                     The resulting grid of concentrations will be  of shape (2^<val>, 2^<val>)"
    print *, "-i <val>, --init=<val>             Sets the grid initialisation type to <val>"
    print *, "                                     EG: -i r gives an initial grid containing random concentrations"
    print *, "-t <arr>, --time_array=<arr>       Sets the array of output times to <arr>"
    print *, "                                     <arr> is a space separated list inside curly braces"
    print *, "                                     EG: -t {0.0, 1.0, 2.0}", newline
    print *, "Multiple characters in the 'short form' -<key> <val> key will be expanded"
    print *, "  EG: -lap0 1.0 is equivalent to -l 1.0 -a 1.0 -p0 1.0"
  end subroutine


  subroutine allocate_t_array(str_array)
    character(*), intent(in) :: str_array
    integer :: n_elements, idx, last_idx, t_idx
    n_elements = 1
    do idx=2, len(str_array) - 1
      ! Find length of t array to allocate
      if (str_array(idx:idx) == " ") n_elements = n_elements + 1
    end do
    
    allocate(cmd_timearray(n_elements))

    t_idx = 1
    last_idx = 2
    do idx=2, len(str_array) - 1
      ! Find length of t array to allocate
      if (str_array(idx:idx) == " ") then
        cmd_timearray(t_idx) = str_to_real(str_array(last_idx : idx - 1))
        t_idx = t_idx + 1
        last_idx = idx + 1
      end if
    end do
    cmd_timearray(t_idx) = str_to_real(str_array(last_idx : idx - 1))
  end subroutine

  function str_to_real(str_val) result(real_val)
    character(*), intent(in) :: str_val
    real(dp) :: real_val
    read(str_val, *) real_val
  end function

  function str_to_int(str_val) result(int_val)
    character(*), intent(in) :: str_val
    integer :: int_val

    read(str_val, *) int_val
  end function
end module command_line