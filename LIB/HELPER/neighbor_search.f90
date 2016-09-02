! ********************************
! 2D AMR prototype
! --------------------------------
!
! test function
!
! ********************************

subroutine neighbor_search()

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik)                    :: k, N, block_num, l

!    N               = size(blocks_params%active_list, dim=1)
!
!    do k = 1, N
!
!        block_num = blocks_params%active_list(k)
!
!        do l = 1, 8
!            if (blocks(block_num)%neighbor_id(l) == -1) then
!                blocks(block_num)%data1 = 7.0_rk
!                print*, 'error: no neighbor found'
!                call save_data(1, 1.0_rk)
!                stop
!            end if
!        end do
!
!    end do

end subroutine neighbor_search
