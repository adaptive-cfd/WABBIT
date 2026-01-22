!! input:    - 2 size N arrays
!! output:   - .true. if arrays are equal
! ********************************************************************************************

logical function array_compare(array1, array2, N)

    use module_params                                               ! global parameters

    implicit none

    integer(kind=ik), intent(in)        :: N
    integer(kind=ik), intent(in)        :: array1(N), array2(N)
    integer(kind=ik)                    :: i                        ! loop variable

    array_compare = .true.

    do i = 1, N
        if ( array1(i) /= array2(i) ) array_compare = .false.
    enddo

end function array_compare
