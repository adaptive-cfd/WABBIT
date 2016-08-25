! ********************************
! 2D AMR prototype
! --------------------------------
! 
! matrix multiplication
!
! name: matrix_mult.f90
! date: 09.08.2016
! author: msr
! version: 0.1
! 
! ********************************

subroutine matrix_mult(mult_AB, A, B, N)

    implicit none

    integer, intent(in)				            :: N
    real(8), dimension(N, N), intent(in)		:: A, B
    real(8), dimension(N, N), intent(out)	    :: mult_AB

    real(8)					                    :: z
    integer					                    :: i, j, k

    mult_AB = 0.0_8

    do i = 1, N
        do j = 1, N
            z = 0.0_8
            do k = 1, N
                z = z + A(i,k) * B(k,j)
            end do
            mult_AB(i,j) = z
        end do
    end do

end subroutine matrix_mult
