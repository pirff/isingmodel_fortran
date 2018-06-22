program mc1
  implicit none

  integer, parameter :: nmc = 20       ! number of trials
  integer, parameter :: n = 100000      ! sample size
  real, dimension(nmc) :: mu            ! result of each trial
  real :: mean, stdev                   ! mean and standard deviation
  integer :: j
  real temp

  !$OMP PARALLEL DO
  do j = 1, nmc
     print *, 'Experiment', j
     temp = 2.0 + j*0.01
     print *, metropolis(temp)
  end do
  !$OMP END PARALLEL DO

contains

  function metropolis(temperature) result(y)
    real, intent(in) :: temperature
    real :: y, u

    y = temperature

  end function metropolis


  ! Draws n Uniform(0,1) random numbers and returns the sample mean
  function monte_carlo(n) result(y)
    integer, intent(in) :: n
    integer :: i
    real :: y, u

    y = 0.d0
    do i = 1, n
       call random_number(u)            ! draw from Uniform(0,1)
       y = y + u                        ! sum the draws
    end do
    y = y / dble(n)                     ! return the sample mean
  end function monte_carlo

end program mc1
