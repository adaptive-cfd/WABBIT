! ********************************
! WABBIT
! --------------------------------
!
! compare 2 arrays
!
! name: array_compare.f90
! date: 29.09.2016
! author: msr
! version: 0.2
!
! ********************************

logical function array_compare(array1, array2, N)

    implicit none

    integer, intent(in)                   :: N
    integer, dimension(N), intent(in)     :: array1, array2

    integer                               :: i

    array_compare = .true.
    do i = 1, N
        if ( array1(i) /= array2(i) ) array_compare = .false.
    enddo

end function array_compare
