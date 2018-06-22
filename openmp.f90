program mc1
  ! implicit none

  integer, parameter :: nmc = 10       ! number of trials
  integer, parameter :: n = 100000      ! sample size
  real, dimension(nmc) :: mu            ! result of each trial
  real :: mean, stdev                   ! mean and standard deviation
  integer :: j
  real :: temperature

  !$OMP PARALLEL DO
  do j = 1, nmc
     print *, 'Experiment', j
     temperature = 2.0 + J*0.01
     print *, 'Temperature',temperature
     print *, metropolis(temperature)
  end do
  !$OMP END PARALLEL DO


contains

  integer function nextIndex(i,L) result (r)
    integer i
    integer L
    if (i<L) THEN
      r = i+1
    ELSE
      r = 1
    END if
  END function

  real function beforeIndex(i,L) result (r)
    integer i
    integer L
    if (i>1) THEN
      r = i-1
    ELSE
      r = L
    END if
  END function

  subroutine init_random_seed()

        INTEGER :: i, n, clock
        INTEGER, DIMENSION(:), ALLOCATABLE :: seed

        CALL RANDOM_SEED(size = n)
        ALLOCATE(seed(n))

        CALL SYSTEM_CLOCK(COUNT=clock)

        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        CALL RANDOM_SEED(PUT = seed)

        DEALLOCATE(seed)
  end

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

  function metropolis(temperature) result(y)
    ! integer, intent(in) :: n
    real :: y, u, temperature
    real, parameter :: Kb = 1.0
    real, parameter :: J = 1.0
    integer, parameter :: L = 20 ! dimension
    integer, parameter :: TimeEquil = 1000 ! time to equilibrium
    integer, parameter :: TimeStat = 2000 ! statistics time
    integer, parameter :: WriteToCsv = 1 ! 0 - do not write to csv
    real  :: x,beta, r(L,L)
    double precision :: deltaE, prob, en, magn, en2, magn2, c1,c2, totalE, totalM
    integer :: lattice(L,L), image(L,L), row(L)
    integer :: q,qq,i,k,neighbors,spin,w,z,aaa
    integer size
    character(len=10) :: name

    call init_random_seed()
    lattice = 1
    
    y = temperature
  end function metropolis

end program mc1
