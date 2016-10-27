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

    integer(kind=ik)                                    :: N, light_id, l, neighbor_light_id, neighbor_number, Bs
    integer(kind=ik), dimension(10)                     :: neighbor, virt1, virt2, neighbor2
    integer(kind=ik), dimension(4,2)                    :: dirsi
    integer(kind=ik), dimension(4)                      :: dirsi2
    character(len=2), dimension(8)                      :: dirs
    logical                                             :: exists

    N               = blocks_params%number_max_blocks
    Bs              = blocks_params%size_block
    dirs            = (/'NO', 'EA', 'SO', 'WE', 'NE', 'NW', 'SE', 'SW'/)
    dirsi           = reshape( (/0, 1, 2, 0, 1, 3, 3, 2/), (/4, 2/) )
    dirsi2          = (/1, 0, 3, 2/)
    exists          = .false.
    neighbor_number = 0

    ! loop over all blocks
    do light_id = 1, N

        if ( blocks(light_id)%active == .true. ) then

            ! reset neighbor lists
            blocks(light_id)%neighbor_treecode(:,:) = -1
            blocks(light_id)%neighbor_dir(:)        = ""
            blocks(light_id)%neighbor_id(:)         = -1
            blocks(light_id)%neighbor2_treecode(:,:)= -1
            blocks(light_id)%neighbor2_dir(:)       = ""
            blocks(light_id)%neighbor2_id(:)        = -1

            if ( blocks_params%size_domain == blocks_params%size_block ) then
                ! only one block, so all neighbors are this block
                do l = 1, 8

                    blocks(light_id)%neighbor_treecode(l,:)    = blocks(light_id)%treecode
                    blocks(light_id)%neighbor_dir(l)           = dirs(l)
                    blocks(light_id)%neighbor_id(l)            = light_id
                end do

            else
                ! more than one block
                ! loop over all directions
                do l = 1, 8

                    neighbor_number = 1

                    ! calculate treecode for neighbor block on same level
                    call adjacent_block(blocks(light_id)%treecode, neighbor, dirs(l))

                    call does_block_exist(neighbor, exists)

                    if (exists) then
                        ! one neighbor block on same level exists
                        call find_block_id(neighbor, neighbor_light_id)

                        blocks(light_id)%neighbor_treecode(l,:)    = neighbor
                        blocks(light_id)%neighbor_dir(l)           = dirs(l)
                        blocks(light_id)%neighbor_id(l)            = neighbor_light_id

                    else
                        ! neighbor block has different level, 1 or 2 neighbors
                        ! check if neighbor is on lower level
                        neighbor(blocks(light_id)%level) = -1
                        call does_block_exist(neighbor, exists)

                        if (exists) then
                            ! on neighbor on lower level exists
                            call find_block_id(neighbor, neighbor_light_id)

                            blocks(light_id)%neighbor_treecode(l,:)    = neighbor
                            blocks(light_id)%neighbor_dir(l)           = dirs(l)
                            blocks(light_id)%neighbor_id(l)            = neighbor_light_id

                        else
                            ! two neighbors on higher level
                            ! diagonal: only one neighbor on higher level
                            if (l <= 4) then
                                ! neighbors on edge
                                neighbor_number                             = 2

                                ! create virtual treecodes and look for their neighbors
                                virt1                                       = blocks(light_id)%treecode
                                virt1(blocks(light_id)%level+1)             = dirsi(l,1)
                                virt2                                       = blocks(light_id)%treecode
                                virt2(blocks(light_id)%level+1)             = dirsi(l,2)

                                ! find neighbors
                                call adjacent_block(virt1, neighbor, dirs(l))
                                call adjacent_block(virt2, neighbor2, dirs(l))

                                call does_block_exist(neighbor, exists)
                                if (exists) then
                                    call find_block_id(neighbor, neighbor_light_id)
                                else
                                   ! error case
                                    print*, 'error: neighbor is missing (1 of 2 neighbors)'
                                    stop
                                end if

                                blocks(light_id)%neighbor_treecode(l,:)    = neighbor
                                blocks(light_id)%neighbor_dir(l)           = dirs(l)
                                blocks(light_id)%neighbor_id(l)            = neighbor_light_id

                                call does_block_exist(neighbor2, exists)
                                if (exists) then
                                    call find_block_id(neighbor2, neighbor_light_id)
                                else
                                  ! error case
                                    print*, 'error: neighbor is missing (2 of 2 neighbors)'
                                    stop
                                end if

                                blocks(light_id)%neighbor2_treecode(l,:)   = neighbor2
                                blocks(light_id)%neighbor2_dir(l)          = dirs(l)
                                blocks(light_id)%neighbor2_id(l)           = neighbor_light_id

                            else
                                ! diagonal neighbor
                                neighbor_number                             = 1

                                ! create virtual treecode and look for the neighbor
                                virt1                                      = blocks(light_id)%treecode
                                virt1(blocks(light_id)%level+1)            = dirsi2(l-4)

                                call adjacent_block(virt1, neighbor, dirs(l))
                                call find_block_id(neighbor, neighbor_light_id)

                                blocks(light_id)%neighbor_treecode(l,:)    = neighbor
                                blocks(light_id)%neighbor_dir(l)           = dirs(l)
                                blocks(light_id)%neighbor_id(l)            = neighbor_light_id

                            end if

                        end if

                    end if

                    ! save number of neighbors in direction l
                    blocks(light_id)%neighbor_number(l)            = neighbor_number

                end do

            end if

        end if

    end do

end subroutine update_neighbors
