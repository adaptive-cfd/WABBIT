! ********************************
! 2D AMR prototype
! --------------------------------
!
! create list of active blocks
!
! name: new_block.f90
! date: 15.08.2016
! author: msr, engels
! version: 0.1
!
! ********************************

subroutine active_blocks_list()

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik)                            :: i, j, allocate_error
    integer(kind=ik), dimension(:), allocatable :: blocks_list

    ! allocate a list to temporarily hold the active blocks. as we do not
    ! know the number of active blocks (that's why we're here), we allocate alonger list
    ! with the maximum number of blocks allowed
    allocate( blocks_list(1:blocks_params%number_max_blocks), stat=allocate_error)

    ! deallocate old list
    deallocate( blocks_params%active_list, stat=allocate_error )

    ! create new list
    ! 1st step: we count how many active blocks we find
    j = 1
    do i = 1, blocks_params%number_max_blocks
        ! is this block actrive?
        if (blocks(i)%active .eqv. .true.) then
            ! and add it to the list of active blocks
            blocks_list(j) = i
            ! yes, count it
            j = j + 1
        end if
    end do
    ! we counted one too far
    j = j-1


    ! allocate memory, saving new list
    allocate( blocks_params%active_list(1:j), stat=allocate_error )
    ! copy the list of active blocks, but only the first part (1:j)
    blocks_params%active_list(1:j) = blocks_list(1:j)

    ! local variable
    deallocate( blocks_list, stat=allocate_error)

end subroutine active_blocks_list
