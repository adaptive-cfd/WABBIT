! ********************************
! 2D AMR prototype
! --------------------------------
!
! update neighbor relations
!
! name: update_neighbors.f90
! date: 16.08.2016
! author: msr
! version: 0.1
!
! ********************************

subroutine update_neighbors()

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik)                    :: k, N, block_num, l, block_id, neighbor_number
    integer(kind=ik), dimension(10)     :: neighbor, virt1, virt2, neighbor2
    integer(kind=ik), dimension(4,2)    :: dirsi
    integer(kind=ik), dimension(4)      :: dirsi2
    character(len=2), dimension(8)      :: dirs
    logical                             :: exists

    N               = size(blocks_params%active_list, dim=1)
    dirs            = (/'NO', 'EA', 'SO', 'WE', 'NE', 'NW', 'SE', 'SW'/)
    dirsi           = reshape( (/0, 1, 2, 0, 1, 3, 3, 2/), (/4, 2/) )
    dirsi2          = (/1, 0, 3, 2/)
    exists          = .false.
    neighbor_number = 0

    ! loop over all active blocks
    do k = 1, N

        block_num = blocks_params%active_list(k)
        ! reset neighbor lists
        blocks(block_num)%neighbor_treecode(:,:) = -1
        blocks(block_num)%neighbor_dir(:)        = ""
        blocks(block_num)%neighbor_id(:)         = -1
        blocks(block_num)%neighbor2_treecode(:,:)= -1
        blocks(block_num)%neighbor2_dir(:)       = ""
        blocks(block_num)%neighbor2_id(:)        = -1

        ! loop over all directions
        do l = 1, 8

            neighbor_number = 1

            ! calculate treecode for neighbor block on same level
            call adjacent_block(blocks(block_num)%treecode, neighbor, dirs(l))

            call does_block_exist(neighbor, exists)

            if (exists) then
                ! one neighbor block on same level exists
                call find_block_id(neighbor, block_id)

                blocks(block_num)%neighbor_treecode(l,:)    = neighbor
                blocks(block_num)%neighbor_dir(l)           = dirs(l)
                blocks(block_num)%neighbor_id(l)            = block_id

            else
                ! neighbor block has different level, 1 or 2 neighbors
                ! check if neighbor is on lower level
                neighbor(blocks(block_num)%level) = -1
                call does_block_exist(neighbor, exists)

                if (exists) then
                    ! on neighbor on lower level exists
                    call find_block_id(neighbor, block_id)

                    blocks(block_num)%neighbor_treecode(l,:)    = neighbor
                    blocks(block_num)%neighbor_dir(l)           = dirs(l)
                    blocks(block_num)%neighbor_id(l)            = block_id

                else
                    ! two neighbors on higher level
                    ! diagonal: only one neighbor on higher level
                    if (l <= 4) then
                        ! neighbors on edge
                        neighbor_number                             = 2

                        ! create virtual treecodes and look for their neighbors
                        virt1                                       = blocks(block_num)%treecode
                        virt1(blocks(block_num)%level+1)            = dirsi(l,1)
                        virt2                                       = blocks(block_num)%treecode
                        virt2(blocks(block_num)%level+1)            = dirsi(l,2)

                        ! find neighbors
                        call adjacent_block(virt1, neighbor, dirs(l))
                        call adjacent_block(virt2, neighbor2, dirs(l))

                        call does_block_exist(neighbor, exists)
                        if (exists) then
                            call find_block_id(neighbor, block_id)
                        else
                            ! error case
                            print*, 'error: neighbor is missing (1 of 2 neighbors)'
                            stop
                        end if

                        blocks(block_num)%neighbor_treecode(l,:)    = neighbor
                        blocks(block_num)%neighbor_dir(l)           = dirs(l)
                        blocks(block_num)%neighbor_id(l)            = block_id

                        call does_block_exist(neighbor, exists)
                        if (exists) then
                            call find_block_id(neighbor2, block_id)
                        else
                            ! error case
                            print*, 'error: neighbor is missing (2 of 2 neighbors)'
                            stop
                        end if

                        blocks(block_num)%neighbor2_treecode(l,:)   = neighbor2
                        blocks(block_num)%neighbor2_dir(l)          = dirs(l)
                        blocks(block_num)%neighbor2_id(l)           = block_id

                    else
                        ! diagonal neighbor
                        neighbor_number                             = 1

                        ! create virtual treecode and look for the neighbor
                        virt1                                       = blocks(block_num)%treecode
                        virt1(blocks(block_num)%level+1)            = dirsi2(l-4)

                        call adjacent_block(virt1, neighbor, dirs(l))
                        call find_block_id(neighbor, block_id)

                        blocks(block_num)%neighbor_treecode(l,:)    = neighbor
                        blocks(block_num)%neighbor_dir(l)           = dirs(l)
                        blocks(block_num)%neighbor_id(l)            = block_id

                    end if

                end if

            end if

            ! save number of neighbors in direction l
            blocks(block_num)%neighbor_number(l)            = neighbor_number

        end do

    end do

end subroutine update_neighbors
