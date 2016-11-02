! ********************************
! WABBIT
! --------------------------------
!
! update neighbor relations:
! all procs loop over all blocks
! create new com_list after neighbor updating
!
! note: proc_rank and proc_block_id are not stored here
!
! name: update_neighbors.f90
! date: 25.10.2016
! author: msr
! version: 0.3
!
! ********************************

subroutine update_neighbors()

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik)                                    :: N, light_id, l
    character(len=3), dimension(16)                     :: dirs
    logical                                             :: exists

    N               = blocks_params%number_max_blocks
    exists          = .false.

    dirs = (/'__N', '__E', '__S', '__W', '_NE', '_NW', '_SE', '_SW', 'NNE', 'NNW', 'SSE', 'SSW', 'ENE', 'ESE', 'WNW', 'WSW'/)

    ! loop over all blocks
    do light_id = 1, N

        if ( blocks(light_id)%active == .true. ) then

            ! reset neighbor lists
            blocks(light_id)%neighbor_treecode(:,:) = -1
            blocks(light_id)%neighbor_dir(:)        = ""
            blocks(light_id)%neighbor_id(:)         = -1

            if ( blocks_params%size_domain == blocks_params%size_block ) then
                ! only one block, so all neighbors are this block
                ! one block domian has 8 neighbors
                do l = 1, 8

                    blocks(light_id)%neighbor_treecode(l,:)    = blocks(light_id)%treecode
                    blocks(light_id)%neighbor_dir(l)           = dirs(l)
                    blocks(light_id)%neighbor_id(l)            = light_id
                end do

            else
                ! more than one block
                ! find neighbors on edges
                ! north
                call find_neighbor_edge(light_id, '__N')
                ! east
                call find_neighbor_edge(light_id, '__E')
                ! south
                call find_neighbor_edge(light_id, '__S')
                ! west
                call find_neighbor_edge(light_id, '__W')

                ! find neighbors on corners
                ! northeast
                call find_neighbor_corner(light_id, '_NE')
                ! northwest
                call find_neighbor_corner(light_id, '_NW')
                ! southeast
                call find_neighbor_corner(light_id, '_SE')
                ! southwest
                call find_neighbor_corner(light_id, '_SW')

            end if

        end if

    end do

end subroutine update_neighbors
