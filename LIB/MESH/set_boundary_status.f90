! ********************************
! 2D AMR prototype
! --------------------------------
!
! set boundary staus for all blocks
!
! name: set_boundary_status.f90
! date: 16.08.2016
! author: msr
! version: 0.1
!
! ********************************

subroutine set_boundary_status()

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik)                    :: k, N, block_num, Bs

    N               = size(blocks_params%active_list, dim=1)
    Bs              = blocks_params%size_block

    ! loop over all active blocks
    do k = 1, N

        block_num = blocks_params%active_list(k)

        ! reset boundary status
        blocks(block_num)%boundary(:)            = .false.

        if (N==1) then
            ! if only one block, all neighbors are across the  boundary
            blocks(block_num)%boundary(:)               = .true.
        else
            ! more than one block
            ! east
            if ( blocks(block_num)%coord_x(Bs) == params%Lx ) blocks(block_num)%boundary(2) = .true.
            ! west
            if ( blocks(block_num)%coord_x(1) == 0.0_rk ) blocks(block_num)%boundary(4) = .true.
            ! south
            if ( blocks(block_num)%coord_y(Bs) == 0.0_rk ) blocks(block_num)%boundary(3) = .true.
            ! north
            if ( blocks(block_num)%coord_y(1) == params%Lx ) blocks(block_num)%boundary(1) = .true.

        end if

    end do

end subroutine set_boundary_status
