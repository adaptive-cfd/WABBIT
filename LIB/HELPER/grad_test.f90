! ********************************
! 2D AMR prototype
! --------------------------------
!
! test function
!
! ********************************

subroutine grad_test()

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik) :: k, N, block_num, l, neighbor_id, level, neighbor_level

!    N  = size(blocks_params%active_list, dim=1)
!
!    ! simple test
!    do k = 1, N
!
!        block_num = blocks_params%active_list(k)
!
!        level = blocks(block_num)%level
!
!        if (level < 1) then
!            print*, "error: block has level < 1"
!            stop
!        end if
!
!        if ( blocks(block_num)%refinement /= -2) then
!            level = level + blocks(block_num)%refinement
!        else
!            level = level - 1
!        end if
!
!        do l = 1, 8
!            neighbor_id = blocks(block_num)%neighbor_id(l)
!            if (neighbor_id == -1) then
!                print*, "TEST error 1"
!                blocks(block_num)%data1 = 7.0_rk
!                print*, blocks(block_num)%treecode
!                print*, blocks(block_num)%neighbor_id
!                print*, blocks(block_num)%neighbor2_id
!                call save_data(1, 1.0_rk)
!                stop
!            end if
!            neighbor_level = blocks(neighbor_id)%level
!
!            if (neighbor_level < 1) then
!                print*, "error: neighbor has level < 1"
!                stop
!            end if
!
!            if ( blocks(neighbor_id)%refinement /= -2 ) then
!                neighbor_level = neighbor_level + blocks(neighbor_id)%refinement
!            else
!                neighbor_level = neighbor_level - 1
!            end if
!
!            if ( abs(level - neighbor_level) > 1 ) then
!                print*, "TEST error 2"
!                    blocks(block_num)%data1 = 7.0_rk
!                    print*, blocks(block_num)%treecode
!                    print*, blocks(block_num)%neighbor_id
!                    print*, blocks(block_num)%neighbor2_id
!                    print*, blocks(block_num)%level
!                    print*, blocks(block_num)%refinement
!                    print*, neighbor_id
!                    print*, blocks(neighbor_id)%level
!                    print*, blocks(neighbor_id)%refinement
!                    call save_data(1, 1.0_rk)
!                    stop
!            end if
!        end do
!
!        do l = 1, 4
!            neighbor_id = blocks(block_num)%neighbor2_id(l)
!            if (neighbor_id /= -1) then
!
!                neighbor_level = blocks(neighbor_id)%level
!
!                if ( blocks(neighbor_id)%refinement /= -2 ) then
!                    neighbor_level = neighbor_level + blocks(neighbor_id)%refinement
!                else
!                    neighbor_level = neighbor_level - 1
!                end if
!
!                if ( abs(level - neighbor_level) > 1 ) then
!                    print*, "TEST error 3"
!                    blocks(block_num)%data1 = 7.0_rk
!                    print*, blocks(block_num)%treecode
!                    print*, blocks(block_num)%neighbor_id
!                    print*, blocks(block_num)%neighbor2_id
!                    call save_data(1, 1.0_rk)
!                    stop
!                end if
!            end if
!        end do
!
!    end do

end subroutine grad_test
