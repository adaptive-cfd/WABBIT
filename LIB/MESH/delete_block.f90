! ********************************
! WABBIT
! --------------------------------
!
! delete block
!
! name: delete_block.f90
! date: 30.09.2016
! author: msr
! version: 0.2
!
! ********************************

subroutine delete_block(block_num)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), intent(in)    :: block_num

    integer(kind=ik)                :: dF

    ! delete data
    do dF = 1, blocks_params%number_data_fields
        blocks(block_num)%data_fields(dF)%data_  = 0.0_rk
    end do
    blocks(block_num)%coord_x                       = 0.0_rk
    blocks(block_num)%coord_y                       = 0.0_rk
    blocks(block_num)%treecode                      = -1
    blocks(block_num)%level                         = -1
    blocks(block_num)%refinement                    = 0
    blocks(block_num)%neighbor_id                   = -1
    blocks(block_num)%neighbor2_id                  = -1
    blocks(block_num)%neighbor_treecode             = -1
    blocks(block_num)%neighbor2_treecode            = -1
    blocks(block_num)%boundary                      = .false.
    blocks(block_num)%boundary2                     = .false.

    ! update active blocks list
    blocks(block_num)%active                        = .false.

end subroutine delete_block
