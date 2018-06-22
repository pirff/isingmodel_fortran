program main
  use omp_lib

  implicit none

  integer ( kind = 4 ), parameter :: L = 10
  double precision, parameter :: temperature = 2.00
  integer, parameter :: TimeEquil = 1 ! time to equilibrium
  integer, parameter :: TimeStat = 2000 ! statistics time

  real :: wtime,x
  integer :: lattice(L,L),z,i,k,spin,neighbors
  double precision :: deltaE


  write ( *, '(a,i8,a,i8,a)' ) '  Lattice of size ', l, ' by ', l, ' points.'
  write ( *, '(a,f0.9)' ) &
    '  Temperature = ', temperature
  write ( *, '(a,i8)' ) &
    '  The number of processors available = ', omp_get_num_procs ( )
  write ( *, '(a,i8)' ) &
    '  The number of threads available    = ', omp_get_max_threads ( )

  lattice = 1
  wtime = omp_get_wtime ( )

  !$omp parallel do
  do z = 1, L*L
    call random_number(x)
    i = floor(x*L)+1
    call random_number(x)
    k = floor(x*L)+1
    spin = lattice(i,k)

    !!!!!!!!!! computing sum of neighbors
    neighbors = 0
    if (i<L) THEN
      neighbors = neighbors + lattice(i+1,k)
    ELSE
      neighbors = neighbors + lattice(1,k)
    end if
    if (i>1) THEN
      neighbors = neighbors + lattice(i-1,k)
    ELSE
      neighbors = neighbors + lattice(L,k)
    end if
    if (k<L) THEN
      neighbors = neighbors + lattice(i,k+1)
    ELSE
      neighbors = neighbors + lattice(i,1)
    end if
    if (k>1) THEN
      neighbors = neighbors + lattice(i,k-1)
    ELSE
      neighbors = neighbors + lattice(i,L)
    end if

    !!!!!!!!!!! estimating delta energy
    deltaE = -2.0 * spin * neighbors
  end do
  !$omp end parallel do

  wtime = omp_get_wtime ( ) - wtime
  write ( *, '(a,g14.6)' ) '  Wall clock time = ', wtime
end
