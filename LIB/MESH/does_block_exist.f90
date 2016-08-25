! ********************************
! 2D AMR prototype
! --------------------------------
!
! return true, if block exist
!
! name: update_neighbors.f90
! date: 16.08.2016
! author: msr
! version: 0.1
!
! ********************************

subroutine does_block_exist(treecode, exists)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), dimension(10), intent(in)     :: treecode
    logical, intent(out)                            :: exists

    integer(kind=ik)                                :: k, N, block_num
    logical                                         :: array_compare

    exists      = .false.
    N           = size(blocks_params%active_list, dim=1)

    do k = 1, N

        block_num   = blocks_params%active_list(k)

        if (array_compare(blocks(block_num)%treecode, treecode, 10)) exists = .true.

    end do

end subroutine does_block_exist
