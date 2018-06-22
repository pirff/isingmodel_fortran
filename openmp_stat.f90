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

! subroutine to initialize a random seed
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

program data
  real, parameter :: Kb = 1.0
  real, parameter :: J = 1.0
  integer, parameter :: L = 20 ! dimension
  integer, parameter :: TimeEquil = 5000 ! time to equilibrium
  integer, parameter :: TimeStat = 20000 ! statistics time
  integer, parameter :: WriteToCsv = 1 ! 0 - do not write to csv
  real  :: x,beta, temperature, r(L,L)
  double precision :: deltaE, prob, en, magn, en2, magn2, c1,c2, totalE, totalM
  integer :: lattice(L,L), image(L,L), row(L)
  integer :: q,qq,i,k,neighbors,spin,w,z,aaa
  integer size
  character(len=10) :: name
  call init_random_seed()

  call random_number(r)

  ! RANDOM initial
  !lattice = nint(r)*2 - 1

  ! POSITIVE INITIAL
  lattice = 1

  ! NEGATIVE INITIAL
  !lattice = -1


  open(unit=1,file='data9.csv',status='unknown',ACCESS='APPEND')
  !open(unit=1,file='data.csv',status='replace')

  do z=0,2
    !call random_number(r)
    !lattice = nint(r)*2 - 1
    lattice = 1
    temperature = 2.0 + z*0.01
    ! print *, "**************"
    ! print *, temperature
    beta = 1 / temperature
    !****************** making sweep to reach equilibrium
    do w = 1, TimeEquil
      do q = 1, L
        do qq = 1, L
          call random_number(x)
          i = floor(x*L)+1
          call random_number(x)
          k = floor(x*L)+1
          !print *, k
          !print *, i
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
    print *,"   reached equilibrium"

    !************** making sweep for statistics
    !*** defining some coefficients
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
    print *, "  statistics computed"
    if (WriteToCsv>0) THEN
      write (1, "(5(f0.9,',',:))") temperature,en, en2, magn, magn2
    END IF
    print *, temperature,en, en2, magn, magn2
  end do


  close(1)

end program data
