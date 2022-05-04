module setup

  use globals
  implicit none

  contains

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
        call logger%info("setup_grid", "Generating Random grid")
        call random_number(c)
        c = 0.4_DP*c - 0.2_DP

      case ("c") ! Circle
        call logger%info("setup_grid", "Generating Circle grid")
        c = -1.0_dp

        do i = 1, grid
          do j = 1, grid

            if( ((real(i, dp) - real(grid, dp)/2 - 0.5_DP))**2 + ((real(j, dp) - real(grid, dp)/2 - 0.5_DP))**2 &
                    .lt. real(grid, dp)**2*0.0625_DP) then
              c(i, j) = 1.0_dp
            end if

          end do
        end do

      case ("b") ! Bar
        call logger%info("setup_grid", "Generating Bar grid")

        do i = 1, int(real(grid, dp)*0.4_dp)
          c(i, :) = -1.0_dp
        end do

        do i = int(real(grid, dp)*0.6_dp), grid
          c(i, :) = -1.0_dp
        end do
        
        do i = int(real(grid, dp)*0.4_dp), int(real(grid, dp)*0.6_dp)
          
          do j = 1, int(real(grid, dp)*0.25_dp) 
            c(i, j) = -1.0_dp
          end do
        
          do j = int(real(grid, dp)*0.25_dp), int(real(grid, dp)*0.75_dp)
            c(i,j) = 1.0_dp
          end do
          
          do j = int(real(grid, dp)*0.75_dp), grid 
            c(i, j) = -1.0_dp
          end do

        end do

      case ("s") ! Split
        call logger%info("setup_grid", "Generating Split grid")
        do i = 1, int(real(grid, dp)*0.5_dp)
          c(i, :) = -1.0_dp
        end do

        do i = int(real(grid, dp)*0.5_dp + 1.0_dp), grid
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
      T(i) = start + real(i, dp)*(end - start)/real(nsteps, dp) 
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

    call lin_tspace(0.0_DP, log10(end - start + 1.0_dp), nsteps, T)
    print *, T
    T = 10.0_dp**T
    T = T - 1.0_dp + start

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

