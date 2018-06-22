program mc1
  use omp_lib
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
    double precision :: y(5)

    real, parameter :: Kb = 1.0
    real, parameter :: J = 1.0
    integer, parameter :: L = 20 ! dimension
    integer, parameter :: TimeEquil = 1000 ! time to equilibrium
    integer, parameter :: TimeStat = 5000 ! statistics time
    real  :: x,beta, r(L,L)
    double precision :: deltaE, prob, en, magn, en2, magn2, c1,c2, totalE, totalM
    integer :: lattice(L,L), image(L,L), row(L)
    integer :: q,qq,i,k,neighbors,spin,w,z,aaa

    y = temperature

    lattice = 1
    beta = 1 / temperature

    do w = 1, TimeEquil
      do q = 1, L
        do qq = 1, L
          call random_number(x)
          i = floor(x*L)+1
          call random_number(x)
          k = floor(x*L)+1
          spin = lattice(i,k)
          neighbors = lattice(nextIndex(i,L),k) + lattice(i,nextIndex(k,L))
          aaa = beforeIndex(k,L)
          neighbors = neighbors + lattice(i,aaa)
          aaa = beforeIndex(i,L)
          neighbors = neighbors + lattice(aaa,k)
          deltaE = -2.0 * J * spin * neighbors
          IF (deltaE >= 0) THEN
            lattice(i,k) = -spin
          ELSE
            prob = exp(deltaE * beta)
            call random_number(x)
            IF (x<prob) THEN
              lattice(i,k) = -spin
            END IF
          END IF
        end do
      end do
    end do

    en = 0
    en2 = 0
    magn = 0
    magn2 = 0
    do w = 1, TimeStat
      ! sweep
      do q = 1, L
        do qq = 1, L
          call random_number(x)
          i = floor(x*L) +1
          call random_number(x)
          k = floor(x*L) +1
          spin = lattice(i,k)
          neighbors = lattice(nextIndex(i,L),k) + lattice(i,nextIndex(k,L))
          aaa = beforeIndex(k,L)
          neighbors = neighbors + lattice(i,aaa)
          aaa = beforeIndex(i,L)
          neighbors = neighbors + lattice(aaa,k)
          deltaE = -2.0 * J * spin * neighbors
          IF (deltaE >= 0) THEN
            lattice(i,k) = -spin
          ELSE
            prob = exp(deltaE * beta)
            call random_number(x)
            IF (x<prob) THEN
              lattice(i,k) = -spin
            END IF
          END IF
        end do
      end do
      ! statistics
      totalE = 0
      totalM = 0
      do i=1,l
        do k=1,l
          spin = lattice(i,k)
          neighbors = lattice(nextIndex(i,L),k) + lattice(i,nextIndex(k,L))
          aaa = beforeIndex(k,L)
          neighbors = neighbors + lattice(i,aaa)
          aaa = beforeIndex(i,L)
          neighbors = neighbors + lattice(aaa,k)
          totalE = totalE - J * spin * neighbors
          totalM = totalM + spin
        end do
      end do
      en = en + totalE
      en2 = en2 + totalE*totalE
      magn = magn + totalM
      magn2 = magn2 + totalM*totalM
    end do
    en = en / (L*L)
    en2 = en2 / (L*L*L*L)
    magn = magn / (L*L)
    magn2 = magn2 / (L*L*L*L)
    en = en / TimeStat
    en2 = en2 / TimeStat
    magn = magn / TimeStat
    magn2 = magn2 / TimeStat


    y(2) = en
    y(3) = en2
    y(4) = magn
    y(5) = magn2

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


end program mc1
