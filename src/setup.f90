module setup

  use globals
  implicit none

  contains
  
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

  !> @brief Initializes Grid
  !! @param[in] n integer, sets up grid to be divided into 2^n x 2^n squares
  !! @param[in] init char with options [r, c, b, s] 
  !! @param[in] init option r: Uniform randomly initialising a number in (-0.2, 0.2)
  !! @param[in] init option c: A centered circle of radius 0.25 with concentration 1 and -1 everywhere else
  !! @param[in] init option b: A centered bar with dimensions (0.2, 0.5) with concetration 1 and -1 everywhere else
  !! @param[in] init option s: Splits the grid into left half with -1 and right half with 1
  !! @param[out] c allocatable rank 2 real64 array that stores concentration
  subroutine setup_grid(c, n, init)
    integer, intent(in) :: n
    character, intent(in) :: init
    real(dp), intent(out), allocatable, dimension(:,:) :: c

    integer :: i, j, grid

    grid = 2**n

    allocate(c(grid,grid))

    select case(init)

      case ("r") ! Random
        call random_number(c)
        c = 0.4_DP*c - 0.2_DP

      case ("c") ! Circle
        c = -1.0_dp

        do i = 1, grid
          do j = 1, grid

            if( ((i - grid/2 - 0.5_DP))**2 + ((j - grid/2 - 0.5_DP))**2 .lt. grid**2*0.0625_DP) then
              c(i, j) = 1.0_dp
            end if

          end do
        end do

      case ("b") ! Bar

        do i = 1, int(grid*0.4)
          c(i, :) = -1.0_dp
        end do

        do i = int(grid*0.6), grid
          c(i, :) = -1.0_dp
        end do
        
        do i = int(grid*0.4), int(grid*0.6)
          
          do j = 1, int(grid*0.25) 
            c(i, j) = -1.0_dp
          end do
        
          do j = int(grid*0.25), int(grid*0.75)
            c(i,j) = 1.0_dp
          end do
          
          do j = int(grid*0.75), grid 
            c(i, j) = -1.0_dp
          end do

        end do

      case ("s") ! Split

        do i = 1, int(grid*0.5)
          c(i, :) = -1.0_dp
        end do

        do i = int(grid*0.5 + 1), grid
          c(i, :) = 1.0_dp
        end do


      end select


  end subroutine

  !> @brief Linearly interpolate Array of sample times
  !! @param[in] start, the starting time
  !! @param[in] end, the finishing time
  !! @param[in] nsteps, the number of steps between start to finish
  !! @param[out] T, Result Array of Sample times with length nsteps+1
  subroutine lin_tspace(start, end, nsteps, T)
    real(dp), intent(in) :: start, end
    integer, intent(in) :: nsteps
    real(dp), intent(out), allocatable :: T(:)

    integer :: i
    call t_validation(start, end, nsteps)

    allocate(T(0:nsteps))

    T(0) = start

    do i = 1, nsteps 
      T(i) = start + (i)*(end - start)/(nsteps) 
    end do

  end subroutine

  
  !> @brief Loglinearly (base10) interpolate Array of sample times
  !! @param[in] start, the starting time
  !! @param[in] end, the finishing time
  !! @param[in] nsteps, the number of steps between start to finish
  !! @param[out] T, Result Array of Sample times with length nsteps+1
  subroutine log_tspace(start, end, nsteps, T)
    real(dp), intent(in) :: start, end
    integer, intent(in) :: nsteps
    real(dp), intent(out), allocatable :: T(:)
    call t_validation(start, end, nsteps)

    call lin_tspace(0.0_DP, log10(end - start + 1), nsteps, T)
    print *, T
    T = 10**T
    T = T - 1 + start

  end subroutine 


  subroutine t_validation(start, end, nsteps)
    real(dp), intent(in) :: start, end
    integer, intent(in) :: nsteps
    real(dp), parameter :: REAL_TOL = 1e-10_DP
    logical :: error

    error = .FALSE.

    ! Start validation
    if (start < 0) then
      call logger%error("tspace_validation", "Start must be greater than zero")
      error = .TRUE.
    end if

    ! End Validation

    if (end < start) then
      call logger%error("tspace_validation", "End must be greater than Start")
      error = .TRUE.
    end if

    if (abs(start - end) .LT. REAL_TOL) then
      call logger%error("tspace_validation", "Start and End cannot be equal")
      error = .TRUE.
    end if
    ! Nsteps validation
    if (nsteps <= 0) then
      call logger%error("tspace_validation", "Number of steps must be at least 1")
      error = .TRUE.
    end if

    if (error) then
      call logger%fatal("tspace_validation", "Invalid timespace setup")
      stop
    end if
  end subroutine
end module

