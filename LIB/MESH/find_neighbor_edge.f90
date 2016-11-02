! ********************************
! WABBIT
! --------------------------------
!
! find neighbor on block edge
! valid cases for edge neighbors:
! i) same level: allways exact one neighbor
! ii) one level down: one neighbor, two possible
! neighbor relations
! iii) one level up: allways two neighbors
!
! name: find_neighbor_edge.f90
! date: 31.10.2016
! author: msr
! version: 0.3
!
! ********************************

subroutine find_neighbor_edge(light_id, dir)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), intent(in)      :: light_id
    character(len=3), intent(in)      :: dir

    integer(kind=ik)                  :: list_id, list_id2, virt_list_id1, virt_list_id2, virt_code1, virt_code2, neighbor_light_id
    integer(kind=ik), dimension(10)   :: neighbor, virt_treecode
    logical                           :: exists
    character(len=3)                  :: dir2, virt_dir1, virt_dir2

    ! -------------------------------------------------------------------------------------------------------------------------
    ! dirs = (/'__N', '__E', '__S', '__W', '_NE', '_NW', '_SE', '_SW', 'NNE', 'NNW', 'SSE', 'SSW', 'ENE', 'ESE', 'WNW', 'WSW'/)
    ! -------------------------------------------------------------------------------------------------------------------------

    ! find id in light data list, set virt_code and lvl_down_neighbor
    select case(dir)
        case('__N')
            list_id    = 1
            ! virtual treecodes, dirs, list_ids for neighbors on higher level
            virt_code1    = 0
            virt_dir1     = 'NNW'
            virt_list_id1 = 10
            virt_code2    = 1
            virt_dir2     = 'NNE'
            virt_list_id2 = 9
            ! dir2 for cases with neighbor one level down
            if (blocks(light_id)%treecode( blocks(light_id)%level ) == 0) then
                dir2     = 'NNW'
                list_id2 = 10
            elseif (blocks(light_id)%treecode( blocks(light_id)%level ) == 1) then
                dir2     = 'NNE'
                list_id2 = 9
            end if

        case('__E')
            list_id = 2
            ! virtual treecodes for neighbors on higher level
            virt_code1    = 1
            virt_dir1     = 'ENE'
            virt_list_id1 = 13
            virt_code2    = 3
            virt_dir2     = 'ESE'
            virt_list_id2 = 14
            ! dir2 for cases with neighbor one level down
            if (blocks(light_id)%treecode( blocks(light_id)%level ) == 1) then
                dir2     = 'ENE'
                list_id2 = 13
            elseif (blocks(light_id)%treecode( blocks(light_id)%level ) == 3) then
                dir2     = 'ESE'
                list_id2 = 14
            end if

        case('__S')
            list_id   = 3
            ! virtual treecodes for neighbors on higher level
            virt_code1    = 2
            virt_dir1     = 'SSW'
            virt_list_id1 = 12
            virt_code2    = 3
            virt_dir2     = 'SSE'
            virt_list_id2 = 11
            ! dir2 for cases with neighbor one level down
            if (blocks(light_id)%treecode( blocks(light_id)%level ) == 3) then
                dir2     = 'SSE'
                list_id2 = 11
            elseif (blocks(light_id)%treecode( blocks(light_id)%level ) == 2) then
                dir2     = 'SSW'
                list_id2 = 12
            end if

        case('__W')
            list_id   = 4
            ! virtual treecodes for neighbors on higher level
            virt_code1    = 0
            virt_dir1     = 'WNW'
            virt_list_id1 = 15
            virt_code2    = 2
            virt_dir2     = 'WSW'
            virt_list_id2 = 16
            ! dir2 for cases with neighbor one level down
            if (blocks(light_id)%treecode( blocks(light_id)%level ) == 0) then
                dir2     = 'WNW'
                list_id2 = 15
            elseif (blocks(light_id)%treecode( blocks(light_id)%level ) == 2) then
                dir2     = 'WSW'
                list_id2 = 16
            end if

    end select

    ! calculate treecode for neighbor on same level
    call adjacent_block(blocks(light_id)%treecode, neighbor, dir)

    ! proof existence of neighbor block
    call does_block_exist(neighbor, exists)

    if (exists) then

        ! neighbor on same level
        ! find light block id
        call find_block_id(neighbor, neighbor_light_id)
        ! write light data
        blocks(light_id)%neighbor_treecode(list_id,:)    = neighbor
        blocks(light_id)%neighbor_dir(list_id)           = dir
        blocks(light_id)%neighbor_id(list_id)            = neighbor_light_id

    else

        ! neighbor could be one level down
        neighbor( blocks(light_id)%level ) = -1
        ! proof existence of neighbor block
        call does_block_exist(neighbor, exists)
        if ( exists ) then
            ! neigbor is one level down
            ! find light block id
            call find_block_id(neighbor, neighbor_light_id)
            ! save dir2, list_id2
            blocks(light_id)%neighbor_treecode(list_id2,:)    = neighbor
            blocks(light_id)%neighbor_dir(list_id2)           = dir2
            blocks(light_id)%neighbor_id(list_id2)            = neighbor_light_id

        elseif ( .not.(exists) ) then
            ! 2 neighbors one level up
            ! first neighbor virtual treecode, one level up
            virt_treecode = blocks(light_id)%treecode
            virt_treecode( blocks(light_id)%level+1 ) = virt_code1
            ! calculate treecode for neighbor on same level (virtual level)
            call adjacent_block(virt_treecode, neighbor, dir)
            ! proof existence of neighbor block
            call does_block_exist(neighbor, exists)

            if (exists) then
                ! neigbor is one level up
                ! find light block id
                call find_block_id(neighbor, neighbor_light_id)
                ! write data
                blocks(light_id)%neighbor_treecode(virt_list_id1,:)    = neighbor
                blocks(light_id)%neighbor_dir(virt_list_id1)           = virt_dir1
                blocks(light_id)%neighbor_id(virt_list_id1)            = neighbor_light_id
            else
                ! error case
                print*, 'error: can not find edge neighbor)'
                stop
            end if

            ! second neighbor virtual treecode, one level up
            virt_treecode = blocks(light_id)%treecode
            virt_treecode( blocks(light_id)%level+1 ) = virt_code2
            ! calculate treecode for neighbor on same level (virtual level)
            call adjacent_block(virt_treecode, neighbor, dir)
            ! proof existence of neighbor block
            call does_block_exist(neighbor, exists)

            if (exists) then
                ! neigbor is one level up
                ! find light block id
                call find_block_id(neighbor, neighbor_light_id)
                ! write data
                blocks(light_id)%neighbor_treecode(virt_list_id2,:)    = neighbor
                blocks(light_id)%neighbor_dir(virt_list_id2)           = virt_dir2
                blocks(light_id)%neighbor_id(virt_list_id2)            = neighbor_light_id
            else
                ! error case
                print*, 'error: can not find edge neighbor)'
                stop
            end if

        end if

    end if


end subroutine find_neighbor_edge
