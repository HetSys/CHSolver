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

      case ("r")
        call random_number(c)
        c = 0.4*c - 0.2

      case ("c")
        c = -1.0_dp

        do i = 1, grid
          do j = 1, grid

            if( ((i - grid/2 - 0.5))**2 + ((j - grid/2 - 0.5))**2 .lt. grid**2*0.0625) then
              c(i, j) = 1.0_dp
            end if

          end do
        end do

      case ("b")

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
        
          do j = int(grid*0.25, grid*0.75)
            c(i,j) = 1.0_dp
          end do
          
          do j = int(grid*0.75), grid 
            c(i, j) = -1.0_dp
          end do

        end do

      case ("s")

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
  subroutine linspace(start, end, nsteps, T)
    real(dp), intent(in) :: start, end
    integer, intent(in) :: nsteps
    real(dp), intent(out), allocatable :: T(:)

    integer :: i

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
  subroutine logspace(start, end, nsteps, T)
    real(dp), intent(in) :: start, end
    integer, intent(in) :: nsteps
    real(dp), intent(out), allocatable :: T(:)

    call linspace(0, log10(end - start), nsteps, T)

    T = 10**T
    T = T - 1 + start

  end subroutine 

end module

