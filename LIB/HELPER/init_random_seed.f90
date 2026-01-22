!> \brief Initialize random seed based on system clock
!> Inspired from GCC example on how to initialize a random seed
subroutine init_random_seed()
  implicit none

  INTEGER :: i, n, clock
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed

  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))

  CALL SYSTEM_CLOCK(COUNT=clock)

  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  CALL RANDOM_SEED(PUT = seed)

  DEALLOCATE(seed)

end subroutine init_random_seed


!> \brief Obtain a gaussian-distributed random number using the Box-Muller transform
!> Instead of giving back two numbers directly, this only gives one for simplicity. However two random numbers are sampled per call
!  _8 just means its double precision
function random_gaussian()
  implicit none

  real(kind=8) r1, r2, random_gaussian
  real(kind=8) :: pi  = 4.0_8 * atan(1.0_8)

  call random_number(r1)
  call random_number(r2)
  random_gaussian = sqrt(-2.0_8*log(r1)) * cos(2.0_8*pi*r2)
  ! in theory we could get a second number with
  ! random_gaussian_2 = sqrt(-2.0*log(u2)) * cos(2.0*pi*u1)

end function