! ********************************
! WABBIT
! --------------------------------
!
! find neighbor on block corner
! valid cases for corner neighbors:
! i) same level: allways exact one neighbor
! ii) one level down: one neighbor, if this neighbor
! only on block corner (not additional on block side)
! iii) one level up: allways exact one neighbor
!
! name: find_neighbor_corner.f90
! date: 31.10.2016
! author: msr
! version: 0.3
!
! ********************************

subroutine find_neighbor_corner(light_id, dir)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), intent(in)      :: light_id
    character(len=3), intent(in)      :: dir

    integer(kind=ik)                  :: list_id, virt_code, neighbor_light_id
    integer(kind=ik), dimension(10)   :: neighbor, virt_treecode
    logical                           :: exists, lvl_down_neighbor

    ! -------------------------------------------------------------------------------------------------------------------------
    ! dirs = (/'__N', '__E', '__S', '__W', '_NE', '_NW', '_SE', '_SW', 'NNE', 'NNW', 'SSE', 'SSW', 'ENE', 'ESE', 'WNW', 'WSW'/)
    ! -------------------------------------------------------------------------------------------------------------------------

    lvl_down_neighbor = .false.

    ! find id in light data list, set virt_code and lvl_down_neighbor
    select case(dir)
        case('_NE')
            list_id   = 5
            virt_code = 1
            ! only sister block 1, 2 can have valid NE neighbor at one level down
            if ( (blocks(light_id)%treecode( blocks(light_id)%level ) == 1) .or. (blocks(light_id)%treecode( blocks(light_id)%level ) == 2) ) then
                lvl_down_neighbor = .true.
            end if

        case('_NW')
            list_id   = 6
            virt_code = 0
            ! only sister block 0, 3 can have valid NW neighbor at one level down
            if ( (blocks(light_id)%treecode( blocks(light_id)%level ) == 0) .or. (blocks(light_id)%treecode( blocks(light_id)%level ) == 3) ) then
                lvl_down_neighbor = .true.
            end if

        case('_SE')
            list_id   = 7
            virt_code = 3
            ! only sister block 0, 3 can have valid SE neighbor at one level down
            if ( (blocks(light_id)%treecode( blocks(light_id)%level ) == 0) .or. (blocks(light_id)%treecode( blocks(light_id)%level ) == 3) ) then
                lvl_down_neighbor = .true.
            end if

        case('_SW')
            list_id   = 8
            virt_code = 2
            ! only sister block 1, 2 can have valid NE neighbor at one level down
            if ( (blocks(light_id)%treecode( blocks(light_id)%level ) == 1) .or. (blocks(light_id)%treecode( blocks(light_id)%level ) == 2) ) then
                lvl_down_neighbor = .true.
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
        if ( exists .and. lvl_down_neighbor ) then
            ! find light block id
            call find_block_id(neighbor, neighbor_light_id)
            ! neigbor is one level down
            blocks(light_id)%neighbor_treecode(list_id,:)    = neighbor
            blocks(light_id)%neighbor_dir(list_id)           = dir
            blocks(light_id)%neighbor_id(list_id)            = neighbor_light_id

        elseif ( .not.(exists) ) then
            ! neighbor could be on level up
            ! virtual treecode, one level up
            virt_treecode = blocks(light_id)%treecode
            virt_treecode( blocks(light_id)%level+1 ) = virt_code
            ! calculate treecode for neighbor on same level (virtual level)
            call adjacent_block(virt_treecode, neighbor, dir)
            ! proof existence of neighbor block
            call does_block_exist(neighbor, exists)

            if (exists) then
                ! find light block id
                call find_block_id(neighbor, neighbor_light_id)
                ! neigbor is one level up
                blocks(light_id)%neighbor_treecode(list_id,:)    = neighbor
                blocks(light_id)%neighbor_dir(list_id)           = dir
                blocks(light_id)%neighbor_id(list_id)            = neighbor_light_id
            else
                ! error case
                print*, 'error: can not find corner neighbor)'
                stop
            end if

        end if

    end if


end subroutine find_neighbor_corner
