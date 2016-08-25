! ********************************
! 2D AMR prototype
! --------------------------------
! 
! flip double precision array a, 
! m dimension of a
!
! name: fliplr.f90
! date: 05.08.2016
! author: msr
! version: 0.1
! 
! ********************************

subroutine fliplr(a, m)

    implicit none

    integer, intent(in)				            :: m
    real(kind=8), dimension(m), intent(inout)	:: a(m)

    integer					                    :: i
    real(kind=8), dimension(m)			        :: a_temp

    do i = 1, m
        a_temp(m + 1 - i) = a(i)
    end do
    a = a_temp
   
end subroutine fliplr
