! ********************************
! 2D AMR prototype
! --------------------------------
!
! delete block
!
! name: delete_block.f90
! date: 22.08.2016
! author: msr
! version: 0.1
!
! ********************************

subroutine delete_block(block_num)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), intent(in)    :: block_num

    ! delete data
    blocks(block_num)%data1             = 0.0_rk
    blocks(block_num)%data2             = 0.0_rk
    blocks(block_num)%coord_x           = 0.0_rk
    blocks(block_num)%coord_y           = 0.0_rk
    blocks(block_num)%treecode          = -1
    blocks(block_num)%level             = -1
    blocks(block_num)%refinement        = 0
    blocks(block_num)%neighbor_id       = -1
    blocks(block_num)%neighbor2_id      = -1
    blocks(block_num)%neighbor_treecode = -1
    blocks(block_num)%neighbor2_treecode= -1

    ! update active blocks list
    blocks(block_num)%active        = .false.
    call active_blocks_list()

end subroutine delete_block
