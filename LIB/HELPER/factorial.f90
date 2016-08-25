! ********************************
! 2D AMR prototype
! --------------------------------
! 
! factorial function
!
! name: factorial.f90
! date: 04.08.2016
! author: msr
! version: 0.1
! 
! ********************************

subroutine factorial(n, erg)

    implicit none

    integer, intent(in)	    :: n
    integer, intent(out)	:: erg

    integer		            :: i

    if (n==0) then
        erg = 1
    else
        erg = 1
        do i = 1, n
            erg = erg*i
        end do
    end if

end subroutine factorial
