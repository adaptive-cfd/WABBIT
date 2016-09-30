! ********************************
! WABBIT
! --------------------------------
!
! count active blocks
!
! name: block_count.f90
! date: 30.09.2016
! author: msr
! version: 0.2
!
! ********************************

subroutine block_count(active_blocks)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), intent(inout) :: active_blocks
    integer(kind=ik)                :: k, N

    N               = blocks_params%number_max_blocks
    active_blocks   = 0

    ! loop over all blocks
    do k = 1, N
        if (blocks(k)%active) then
            active_blocks = active_blocks + 1
        end if
    end do

end subroutine block_count
