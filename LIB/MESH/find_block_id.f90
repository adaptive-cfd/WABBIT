! ********************************
! WABBIT
! --------------------------------
!
! find block id to given treecode
!
! name: find_block_id.f90
! date: 29.09.2016
! author: msr
! version: 0.2
!
! ********************************

subroutine find_block_id(treecode, block_id)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), dimension(10), intent(in)     :: treecode
    integer(kind=ik), intent(out)                   :: block_id

    integer(kind=ik)                                :: k, N
    logical                                         :: array_compare

    N           = blocks_params%number_max_blocks
    block_id    = -1

    ! loop over all blocks
    do k = 1, N

        if (blocks(k)%active) then

            if (array_compare(blocks(k)%treecode, treecode, 10)) block_id = k

        end if

    end do

end subroutine find_block_id
