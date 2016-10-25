! ********************************
! WABBIT
! --------------------------------
!
! convert a integer i to binary b
! binary return as vector with length N
!
! name: int_to_binary.f90
! date: 25.10.2016
! author: msr
! version: 0.3
!
! ********************************

subroutine int_to_binary(i, N, b)

    implicit none

    integer, intent(in)                     :: i, N
    integer, dimension(N), intent(out)      :: b

    integer                                 :: j, k

    j = 1
    b = 0
    k = i

    do while (k > 0)
        b(j) = mod(k, 2)
        k = int(k/2)
        j = j + 1
    end do

end subroutine int_to_binary
