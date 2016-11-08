! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: find_neighbor_edge.f90
! version: 0.4
! author: msr
!
! find neighbor on block edge
! valid cases for edge neighbors:
! i) same level: allways exact one neighbor
! ii) one level down: one neighbor, two possible
! neighbor relations
! iii) one level up: allways two neighbors
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
! 07/11/16 - switch to v0.4
! ********************************************************************************************

subroutine find_neighbor_edge(heavy_id, light_id, block_list, max_treelevel, dir, neighbor_list)

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

    ! auxiliary variables
    integer(kind=ik)                    :: list_id, virt_code1, virt_list_id1, virt_code2, virt_list_id2, list_id2

    ! mesh level
    integer(kind=ik)                    :: level
    ! treecode varaibles
    integer(kind=ik)                    :: my_treecode(max_treelevel), neighbor(max_treelevel), virt_treecode(max_treelevel)
    ! return value from function "does_block_exist"
    logical                             :: exists
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

!---------------------------------------------------------------------------------------------
! main body

    ! set virt_code and lvl_down_neighbor
    select case(dir)
        case('__N')
            list_id    = 1
            ! virtual treecodes, list_ids for neighbors on higher level
            virt_code1    = 0
            virt_list_id1 = 10
            virt_code2    = 1
            virt_list_id2 = 9
            ! id2 for cases with neighbor one level down
            if ( block_list(light_id, level) == 0) then
                list_id2 = 10
            elseif ( block_list(light_id, level) == 1) then
                list_id2 = 9
            end if

        case('__E')
            list_id = 2
            ! virtual treecodes for neighbors on higher level
            virt_code1    = 1
            virt_list_id1 = 13
            virt_code2    = 3
            virt_list_id2 = 14
            ! id2 for cases with neighbor one level down
            if ( block_list(light_id, level) == 1) then
                list_id2 = 13
            elseif ( block_list(light_id, level) == 3) then
                list_id2 = 14
            end if

        case('__S')
            list_id   = 3
            ! virtual treecodes for neighbors on higher level
            virt_code1    = 2
            virt_list_id1 = 12
            virt_code2    = 3
            virt_list_id2 = 11
            ! id2 for cases with neighbor one level down
            if ( block_list(light_id, level) == 3) then
                list_id2 = 11
            elseif ( block_list(light_id, level) == 2) then
                list_id2 = 12
            end if

        case('__W')
            list_id   = 4
            ! virtual treecodes for neighbors on higher level
            virt_code1    = 0
            virt_list_id1 = 15
            virt_code2    = 2
            virt_list_id2 = 16
            ! id2 for cases with neighbor one level down
            if ( block_list(light_id, level) == 0) then
                list_id2 = 15
            elseif ( block_list(light_id, level) == 2) then
                list_id2 = 16
            end if

    end select

    ! calculate treecode for neighbor on same level
    call adjacent_block( my_treecode, neighbor, dir, level, max_treelevel)

    ! proof existence of neighbor block and find light data id
    call does_block_exist(neighbor, block_list, max_treelevel, exists, neighbor_light_id)

    if (exists) then

        ! neighbor on same level
        ! write neighbor data, 2D: 16 possible neighbor relations
        neighbor_list( (heavy_id-1)*16 + list_id ) = neighbor_light_id

    else

        ! neighbor could be one level down
        neighbor( level ) = -1
        ! proof existence of neighbor block
        call does_block_exist(neighbor, block_list, max_treelevel, exists, neighbor_light_id)
        if ( exists ) then
            ! neigbor is one level down
            ! save list_id2
            neighbor_list( (heavy_id-1)*16 + list_id2 ) = neighbor_light_id

        elseif ( .not.(exists) ) then
            ! 2 neighbors one level up

            ! first neighbor virtual treecode, one level up
            virt_treecode = my_treecode
            virt_treecode( level+1 ) = virt_code1

            ! calculate treecode for neighbor on same level (virtual level)
            call adjacent_block( virt_treecode, neighbor, dir, level+1, max_treelevel)
            ! proof existence of neighbor block
            call does_block_exist(neighbor, block_list, max_treelevel, exists, neighbor_light_id)

            if (exists) then
                ! neigbor is one level up
                ! write data
                neighbor_list( (heavy_id-1)*16 + virt_list_id1 ) = neighbor_light_id

            else
                ! error case
                print*, 'ERROR: can not find edge neighbor)'
                stop
            end if

            ! second neighbor virtual treecode, one level up
            virt_treecode = my_treecode
            virt_treecode( level+1 ) = virt_code2

            ! calculate treecode for neighbor on same level (virtual level)
            call adjacent_block( virt_treecode, neighbor, dir, level+1, max_treelevel)
            ! proof existence of neighbor block
            call does_block_exist(neighbor, block_list, max_treelevel, exists, neighbor_light_id)

            if (exists) then
                ! neigbor is one level up
                ! write data
                neighbor_list( (heavy_id-1)*16 + virt_list_id2 ) = neighbor_light_id

            else
                ! error case
                print*, 'ERROR: can not find edge neighbor)'
                stop
            end if

        end if

    end if

end subroutine find_neighbor_edge
