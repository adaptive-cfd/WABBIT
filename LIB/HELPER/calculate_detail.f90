! ********************************
! 2D AMR prototype
! --------------------------------
!
! calculate the detail (error) between
! two given 2D fields
! [double precision]
!
! name: calculate_detail.f90
! date: 14.09.2016
! author: msr
! version: 0.1
!
! ********************************

subroutine calculate_detail(detail, phi1, phi2, N)

    implicit none

    integer, intent(in)                    :: N

    real(8), dimension(N,N), intent(in)    :: phi1, phi2
    real(8), intent(out)                   :: detail

    integer                                :: i, j

    detail = 0.0_8

    do i = 1, N
        do j = 1, N
            detail = max( detail, sqrt( (phi1(i,j)-phi2(i,j)) * ( phi1(i,j)-phi2(i,j)) ) )
        end do
    end do

end subroutine calculate_detail
