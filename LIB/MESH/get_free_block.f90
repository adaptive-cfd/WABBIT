! ********************************
! 2D AMR prototype
! --------------------------------
!
! return free block id
!
! name: get_free_block.f90
! date: 22.08.2016
! author: msr
! version: 0.1
!
! ********************************

subroutine get_free_block(id)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), intent(out)   :: id
    integer(kind=ik)                :: k, N

    N  = blocks_params%number_max_blocks
    id = 0

    ! loop over all blocks
    do k = 1, N

        if (blocks(k)%active .eqv. .false.) then
            id = k
            exit
        end if

    end do

end subroutine get_free_block
