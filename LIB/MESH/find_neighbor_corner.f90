! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: find_neighbor_corner.f90
! version: 0.4
! author: msr
!
! find neighbor on block corner
! valid cases for corner neighbors:
! i) same level: allways exact one neighbor
! ii) one level down: one neighbor, if this neighbor
! only on block corner (not additional on block side)
! iii) one level up: allways exact one neighbor
!
! input:    - heavy and light data id
!           - light data array and max treelevel
!           - direction for neighbor search
! output:   - neighbor list array
!
! -------------------------------------------------------------------------------------------------------------------------
! dirs = (/'__N', '__E', '__S', '__W', '_NE', '_NW', '_SE', '_SW', 'NNE', 'NNW', 'SSE', 'SSW', 'ENE', 'ESE', 'WNW', 'WSW'/)
! -------------------------------------------------------------------------------------------------------------------------
!
! = log ======================================================================================
!
! 08/11/16 - switch to v0.4
! ********************************************************************************************

subroutine find_neighbor_corner(heavy_id, light_id, block_list, max_treelevel, dir, neighbor_list)

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! heavy data id
    integer(kind=ik), intent(in)        :: heavy_id
    ! light data id
    integer(kind=ik), intent(in)        :: light_id
    ! max treelevel
    integer(kind=ik), intent(in)        :: max_treelevel
    ! light data array
    integer(kind=ik), intent(in)        :: block_list(:, :)
    ! direction for neighbor search
    character(len=3), intent(in)        :: dir

    ! heavy data array - neifghbor data
    integer(kind=ik), intent(out)       :: neighbor_list(:)

    ! mesh level
    integer(kind=ik)                    :: level
    ! treecode varaibles
    integer(kind=ik)                    :: my_treecode(max_treelevel), neighbor(max_treelevel), virt_treecode(max_treelevel)
    ! return value from function "does_block_exist"
    logical                             :: exists
    ! variable to show if there is a valid corner neighbor
    logical                             :: lvl_down_neighbor

    ! auxiliary variables
    integer(kind=ik)                    :: list_id, virt_code

    ! neighbor light data id
    integer(kind=ik)                    :: neighbor_light_id

!---------------------------------------------------------------------------------------------
! interfaces

    interface
        subroutine does_block_exist(treecode, block_list, max_treelevel, exists, light_id)
            use module_params
            integer(kind=ik), intent(in)        :: max_treelevel
            integer(kind=ik), intent(in)        :: treecode(max_treelevel)
            integer(kind=ik), intent(in)        :: block_list(:, :)
            logical, intent(out)                :: exists
            integer(kind=ik), intent(out)       :: light_id
        end subroutine does_block_exist

    end interface

!---------------------------------------------------------------------------------------------
! variables initialization

    my_treecode     = block_list( light_id, 1:max_treelevel )
    level           = block_list( light_id, max_treelevel + 1 )

    lvl_down_neighbor = .false.

!---------------------------------------------------------------------------------------------
! main body

    ! find id in light data list, set virt_code and lvl_down_neighbor
    select case(dir)
        case('_NE')
            list_id   = 5
            virt_code = 1
            ! only sister block 1, 2 can have valid NE neighbor at one level down
            if ( (my_treecode( level ) == 1) .or. (my_treecode( level ) == 2) ) then
                lvl_down_neighbor = .true.
            end if

        case('_NW')
            list_id   = 6
            virt_code = 0
            ! only sister block 0, 3 can have valid NW neighbor at one level down
            if ( (my_treecode( level ) == 0) .or. (my_treecode( level ) == 3) ) then
                lvl_down_neighbor = .true.
            end if

        case('_SE')
            list_id   = 7
            virt_code = 3
            ! only sister block 0, 3 can have valid SE neighbor at one level down
            if ( (my_treecode( level ) == 0) .or. (my_treecode( level ) == 3) ) then
                lvl_down_neighbor = .true.
            end if

        case('_SW')
            list_id   = 8
            virt_code = 2
            ! only sister block 1, 2 can have valid NE neighbor at one level down
            if ( (my_treecode( level ) == 1) .or. (my_treecode( level ) == 2) ) then
                lvl_down_neighbor = .true.
            end if

    end select

    ! calculate treecode for neighbor on same level
    call adjacent_block( my_treecode, neighbor, dir, level, max_treelevel)

    ! proof existence of neighbor block
    call does_block_exist(neighbor, block_list, max_treelevel, exists, neighbor_light_id)

    if (exists) then

        ! neighbor on same level
        ! write light data
        neighbor_list( (heavy_id-1)*16 + list_id ) = neighbor_light_id

    else

        ! neighbor could be one level down
        neighbor( level ) = -1
        ! proof existence of neighbor block
        call does_block_exist(neighbor, block_list, max_treelevel, exists, neighbor_light_id)

        if ( exists .and. lvl_down_neighbor ) then
            ! neigbor is one level down
            neighbor_list( (heavy_id-1)*16 + list_id ) = neighbor_light_id

        elseif ( .not.(exists) ) then
            ! neighbor could be on level up
            ! virtual treecode, one level up
            virt_treecode = my_treecode
            virt_treecode( level+1 ) = virt_code

            ! calculate treecode for neighbor on same level (virtual level)
            call adjacent_block( virt_treecode, neighbor, dir, level+1, max_treelevel)
            ! proof existence of neighbor block
            call does_block_exist(neighbor, block_list, max_treelevel, exists, neighbor_light_id)

            if (exists) then
                ! neigbor is one level up
                neighbor_list( (heavy_id-1)*16 + list_id ) = neighbor_light_id

            else
                ! error case
                print*, 'ERROR: can not find corner neighbor)'
                stop
            end if

        end if

    end if

end subroutine find_neighbor_corner
