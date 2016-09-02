! ********************************
! 2D AMR prototype
! --------------------------------
!
! absolute value summation over
! all nodes for all blocks
!
! name: blocks_sum.f90
! date: 16.08.2016
! author: msr
! version: 0.1
!
! ********************************

subroutine blocks_sum(s)

    use module_params
    use module_blocks

    implicit none

    real(kind=rk), intent(inout)    :: s
    s = 0.0_rk
!
!    integer                         :: i, j, k, N, block_num

!    N = size(blocks_params%active_list, dim=1)
!    s = 0.0_rk
!
!    do k = 1, N
!
!        block_num = blocks_params%active_list(k)
!
!        do i = 1, blocks_params%size_block
!            do j = 1, blocks_params%size_block
!                s = s + abs(blocks(block_num)%data1(i,j))
!            end do
!        end do
!
!    end do

end subroutine blocks_sum
