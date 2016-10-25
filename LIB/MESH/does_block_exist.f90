! ********************************
! WABBIT
! --------------------------------
!
! return true, if block exist
!
! name: does_block_exist.f90
! date: 25.10.2016
! author: msr, engels
! version: 0.3
!
! ********************************

subroutine does_block_exist(treecode, exists)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), dimension(10), intent(in)     :: treecode
    logical, intent(out)                            :: exists

    integer(kind=ik)                                :: k, N
    logical                                         :: array_compare

    exists      = .false.
    N           = blocks_params%number_max_blocks

    ! loop over all blocks
    do k = 1, N

        if (blocks(k)%active) then

            if (array_compare(blocks(k)%treecode, treecode, 10)) exists = .true.

        end if

    end do

end subroutine does_block_exist
