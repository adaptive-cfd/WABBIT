! ********************************
! 2D AMR prototype
! --------------------------------
!
! absolute value summation over all nodes
! for all blocks for data field dF
!
! sum is scaled with block area and number of nodes
!
! name: blocks_sum.f90
! date: 16.08.2016
! author: msr
! version: 0.1
!
! ********************************

subroutine blocks_sum(s, dF)

    use module_params
    use module_blocks

    implicit none

    real(kind=rk), intent(inout)    :: s
    integer(kind=ik), intent(in)    :: dF

    real(kind=rk)                   :: block_s, block_area
    integer                         :: i, j, k, N, block_num, Bs, g

    N           = size(blocks_params%active_list, dim=1)
    Bs          = blocks_params%size_block
    g           = blocks_params%number_ghost_nodes

    s           = 0.0_rk
    block_s     = 0.0_rk
    block_area  = 0.0_rk

    do k = 1, N

        block_num = blocks_params%active_list(k)

        block_area = abs( blocks(block_num)%coord_x(1) - blocks(block_num)%coord_x(Bs) ) * &
                     abs( blocks(block_num)%coord_y(1) - blocks(block_num)%coord_y(Bs) )

        block_s = 0.0_rk
        do i = g+1, Bs+g
            do j = g+1, Bs+g
                block_s = block_s + abs(blocks(block_num)%data_fields(dF)%data_(i,j))
            end do
        end do

        s = s + block_s / real(Bs*Bs, kind=rk) * real(block_area, kind=rk)

    end do

end subroutine blocks_sum
