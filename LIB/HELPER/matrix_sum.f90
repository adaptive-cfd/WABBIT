! ********************************
! WABBIT
! --------------------------------
! 
! matrix summation
!
! name: matrix_sum.f90
! date: 29.09.2016
! author: msr
! version: 0.2
! 
! ********************************

subroutine matrix_sum(sum_A, A, N)

   implicit none

   integer, intent(in)				        :: N
   real(8), dimension(N, N), intent(in)		:: A
   real(8), intent(out)				        :: sum_A

   integer					                :: i,j

   sum_A = 0.0_8

   do i = 1, N
      do j = 1, N
         sum_A = sum_A + abs(A(i,j))
      end do
   end do

end
