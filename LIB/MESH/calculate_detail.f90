! ********************************
! WABBIT
! --------------------------------
!
! calculate the detail (error) between
! two given 2D fields
! [double precision]
!
! name: calculate_detail.f90
! date: 28.10.2016
! author: msr
! version: 0.3
!
! ********************************

subroutine calculate_detail(detail, phi1, phi2, N)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), intent(in)                 :: N

    real(kind=rk), dimension(N,N), intent(in)    :: phi1, phi2
    real(kind=rk), intent(out)                   :: detail

    integer(kind=ik)                             :: i, j

    detail = 0.0_rk

    do i = 1, N
        do j = 1, N
            detail = max( detail, sqrt( (phi1(i,j)-phi2(i,j)) * ( phi1(i,j)-phi2(i,j)) ) )
        end do
    end do

end subroutine calculate_detail
