! ********************************
! 2D AMR prototype
! --------------------------------
!
! find block id to given treecode
!
! name: find_block_id.f90
! date: 16.08.2016
! author: msr
! version: 0.1
!
! ********************************

subroutine find_block_id(treecode, block_id)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), dimension(10), intent(in)     :: treecode
    integer(kind=ik), intent(out)                   :: block_id

    integer(kind=ik)                                :: k, N, block_num
    logical                                         :: array_compare

    N           = size(blocks_params%active_list, dim=1)
    block_id    = -1

    do k = 1, N

        block_num   = blocks_params%active_list(k)
        if (array_compare(blocks(block_num)%treecode, treecode, 10)) block_id = block_num

    end do

end subroutine find_block_id
