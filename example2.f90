module rng
  implicit none

  private
  public :: rng_t, rng_seed, rng_uniform

  ! Dimension of the state
  integer, parameter :: ns = 4

  ! Default seed vector
  integer, parameter, dimension(ns) :: default_seed &
       = (/ 521288629, 362436069, 16163801, 1131199299 /)

  ! A data type for storing the state of the RNG
  type :: rng_t
     integer, dimension(ns) :: state = default_seed
  end type rng_t

contains

  ! Seeds the RNG using a single integer and a default seed vector.
  subroutine rng_seed(self, seed)
    type(rng_t), intent(inout) :: self
    integer, intent(in) :: seed
    self%state(1) = seed
    self%state(2:ns) = default_seed(2:ns)
  end subroutine rng_seed

  ! Draws a uniform real number on [0,1].
  function rng_uniform(self) result(u)
    type(rng_t), intent(inout) :: self
    real :: u
    integer :: imz

    imz = self%state(1) - self%state(3)

    if (imz < 0) imz = imz + 2147483579

    self%state(1) = self%state(2)
    self%state(2) = self%state(3)
    self%state(3) = imz
    self%state(4) = 69069 * self%state(4) + 1013904243
    imz = imz + self%state(4)
    u = 0.5d0 + 0.23283064d-9 * imz
  end function rng_uniform

end module rng


program mc1
  use rng
  implicit none

  integer, parameter :: nmc = 100       ! number of trials
  integer, parameter :: n = 100000      ! sample size
  type(rng_t), dimension(nmc) :: rng            ! result of each trial
  real :: mean, stdev                   ! mean and standard deviation
  integer :: j

  !$OMP PARALLEL DO
  do j = 1, nmc
     print *, 'Experiment', j
     call rng_seed(rng(j), 932117 + j)
     mu(j) = monte_carlo(rng,n)
  end do
  !$OMP END PARALLEL DO

  mean = sum(mu) / dble(nmc)
  stdev = sqrt( sum( (mu - mean)**2 ) ) / dble(nmc)

  print *, 'mean', mean
  print *, 'stdev', stdev

contains

  ! Draws n Uniform(0,1) random numbers and returns the sample mean
  function monte_carlo(rng, n) result(y)
    type(rng_t), intent(inout) :: rng
    integer, intent(in) :: n
    integer :: i
    real :: y, u

    y = 0.d0
    do i = 1, n
       u = rng_uniform(rng)           ! draw from Uniform(0,1)
       y = y + u                      ! sum the draws
    end do
    y = y / dble(n)                   ! return the sample mean
  end function monte_carlo

end program mc1
