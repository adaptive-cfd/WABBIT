! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: get_heavy_free_block.f90
! version: 0.4
! author: msr, engels
!
! return heavy free block id
!
! input:    - light data list (first column of light data)
!           - size of light data list
! output:   - free heavy data id
!
! = log ======================================================================================
!
! 07/11/16 - switch to v0.4
! ********************************************************************************************

subroutine get_heavy_free_block(id, block_list, N)

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! heavy data id
    integer(kind=ik), intent(out)   :: id
    ! block list length
    integer(kind=ik), intent(in)    :: N
    ! heavy data array - block data
    integer(kind=ik), intent(in)    :: block_list(N)

    ! loop variable, number of blocks (lines) in light data
    integer(kind=ik)                :: k

!---------------------------------------------------------------------------------------------
! variables initialization

    id = -1

!---------------------------------------------------------------------------------------------
! main body

    ! loop over all blocks
    do k = 1, N
        ! check if the block is active. if it is not, then we found a free block
        ! to return
        if (block_list(k) == -1) then
        id = k
        exit
        end if
    end do

    ! error catching: is there no more free blocks on the list?
    if (id == -1) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: We try to fetch a heavy free block ID from the list but all blocks are used."
        stop

    end if

end subroutine get_heavy_free_block
