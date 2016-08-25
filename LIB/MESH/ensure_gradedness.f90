! ********************************
! 2D AMR prototype
! --------------------------------
!
! check the gradedness after new
! refinement status
!
! name: ensure_gradedness.f90
! date: 18.08.2016
! author: msr
! version: 0.1
!
! ********************************

subroutine ensure_gradedness()

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik) :: k, N, block_num, l, neighbor_id, level, neighbor_level, neighbor2_id

    N  = size(blocks_params%active_list, dim=1)

    ! loop over all active blocks
    do k = 1, N

        block_num = blocks_params%active_list(k)

        ! block want to refine
        if (blocks(block_num)%refinement == 1) then

            ! loop over all directions
            do l = 1, 8

                ! check number of neighbors, first case: exact 1 neighbor
                if (blocks(block_num)%neighbor_number(l) == 1) then

                    ! check neighbor treelevel
                    neighbor_id     = blocks(block_num)%neighbor_id(l)
                    level           = blocks(block_num)%level

                    if (neighbor_id == -1) then
                        print*, 'error: gradedness error, can not refine - neighbor does not exists'
                        stop
                    end if

                    neighbor_level  = blocks(neighbor_id)%level

                    if (level == neighbor_level) then
                        ! neighbor on same level
                        ! neighbor can not coarsen, because block want to refine
                        if (blocks(neighbor_id)%refinement == -1) blocks(neighbor_id)%refinement = 0

                    elseif (level - neighbor_level == 1) then
                        ! neighbor one level lower
                        ! neighbor has to refine, because block want to refine

                        blocks(neighbor_id)%refinement = 1

                    elseif (neighbor_level - level == 1) then
                        ! neighbor one level higher
                        ! nothing to do

                    else
                        ! error
                        print*, 'error: neighbor more than one level up/down in tree'
                        stop

                    end if

                ! second case: 2 neighbors
                elseif (blocks(block_num)%neighbor_number(l) == 2) then
                    ! nothing to do, block can refine

                ! error case
                else
                    print*, "error: block has wrong number of neighbors!"
                    stop

                end if
            end do

        ! block want to coarsen
        elseif (blocks(block_num)%refinement == -1) then

            ! loop over all directions
            do l = 1, 8

                ! check number of neighbors, first case: exact 1 neighbor
                if (blocks(block_num)%neighbor_number(l) == 1) then

                    ! check neighbor treelevel
                    neighbor_id     = blocks(block_num)%neighbor_id(l)
                    level           = blocks(block_num)%level

                    if (neighbor_id == -1) then
                        print*, 'error: gradedness error, can not coarsen - neighbor does not exists'
                        stop
                    end if

                    neighbor_level  = blocks(neighbor_id)%level

                    if (level == neighbor_level) then
                        ! neighbor on same level
                        ! block can not coarsen, because neighbor want to refine
                        if (blocks(neighbor_id)%refinement == 1) blocks(block_num)%refinement = 0

                    elseif (level - neighbor_level == 1) then
                        ! neighbor on lower level
                        ! nothing to do

                    elseif (neighbor_level - level == 1) then
                        ! neighbor on higher level
                        ! neighbor want to refine, ...
                        if (blocks(neighbor_id)%refinement == 1) then
                            ! ... so block have also refine
                            blocks(block_num)%refinement = 1
                        else
                            ! neighbor stays on his level or want to coarsen,
                            ! so block can not coarsen
                            ! TODO: block can coarsen, if neighbor also want to coarsen
                            blocks(block_num)%refinement = 0
                        end if

                    else
                        ! error case
                        print*, "error: block has wrong number of neighbors!"
                        stop

                    end if

                ! second case: 2 neighbors
                elseif (blocks(block_num)%neighbor_number(l) == 2) then
                    ! two neighbor, both are one level higher in tree
                    neighbor_id     = blocks(block_num)%neighbor_id(l)

                    if (neighbor_id == -1) then
                        print*, 'error: gradedness error, can not coarsen - neighbor 1 does not exists'
                        stop
                    end if

                    ! look for second neighbor
                    neighbor2_id    = blocks(block_num)%neighbor2_id(l)

                    if (neighbor_id == -1) then
                        print*, 'error: gradedness error, can not coarsen - neighbor 2 does not exists'
                        stop
                    end if

                    if ( (blocks(neighbor_id)%refinement == 1) .or. (blocks(neighbor2_id)%refinement == 1) ) then
                        ! one of the neighbors want to refine,
                        ! so block have also refine
                        blocks(block_num)%refinement = 1
                    else
                        ! neighbors stay on tree level or want to coarsen,
                        ! so block can not coarsen
                        ! TODO: block can coarsen, if both neighbors also want to coarsen
                        blocks(block_num)%refinement = 0

                    end if

                ! error case
                else
                    print*, "error: block has wrong number of neighbors!"
                    stop
                end if

            end do

        end if

    end do

end subroutine ensure_gradedness
