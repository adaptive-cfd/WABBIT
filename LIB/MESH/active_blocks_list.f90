! ********************************
! 2D AMR prototype
! --------------------------------
!
! create list of active blocks
!
! name: new_block.f90
! date: 15.08.2016
! author: msr
! version: 0.1
!
! ********************************

subroutine active_blocks_list()

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik)                            :: i, j, allocate_error
    integer(kind=ik), dimension(:), allocatable :: blocks_list

    ! local variable
    allocate( blocks_list(blocks_params%number_max_blocks), stat=allocate_error)

    ! deallocate old list
    deallocate( blocks_params%active_list, stat=allocate_error )

    ! create new list
    j = 1
    do i = 1, blocks_params%number_max_blocks
        if (blocks(i)%active .eqv. .true.) then
            blocks_list(j) = i
            j = j + 1
        end if
    end do

    ! allocate memory, saving new list
    allocate( blocks_params%active_list(j-1), stat=allocate_error )
    blocks_params%active_list = blocks_list

    ! local variable
    deallocate( blocks_list, stat=allocate_error)

end subroutine active_blocks_list
