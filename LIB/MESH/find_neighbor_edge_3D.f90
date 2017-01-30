! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: find_neighbor_edge_3D.f90
! version: 0.5
! author: msr
!
! find neighbor on block edge
!
! input:    - heavy and light data id
!           - light data array and max treelevel
!           - direction for neighbor search
!           - list of active blocks
! output:   - neighbor list array
!
! --------------------------------------------------------------------------------------------
! neighbor codes:
! ---------------
! for imagination:  - 6-sided dice with '1'-side on top, '6'-side on bottom, '2'-side in front
!                   - edge: boundary between two sides - use sides numbers for coding
!                   - corner: between three sides - so use all three sides numbers
!                   - block on higher/lower level: block shares face/edge and one unique corner,
!                     so use this corner code in second part of neighbor code
!
! faces:  '__1/___', '__2/___', '__3/___', '__4/___', '__5/___', '__6/___'
! edges:  '_12/___', '_13/___', '_14/___', '_15/___'
!         '_62/___', '_63/___', '_64/___', '_65/___'
!         '_23/___', '_25/___', '_43/___', '_45/___'
! corner: '123/___', '134/___', '145/___', '152/___'
!         '623/___', '634/___', '645/___', '652/___'
!
! complete neighbor code array, 74 possible neighbor relations
! neighbors = (/'__1/___', '__2/___', '__3/___', '__4/___', '__5/___', '__6/___', '_12/___', '_13/___', '_14/___', '_15/___',
!               '_62/___', '_63/___', '_64/___', '_65/___', '_23/___', '_25/___', '_43/___', '_45/___', '123/___', '134/___',
!               '145/___', '152/___', '623/___', '634/___', '645/___', '652/___', '__1/123', '__1/134', '__1/145', '__1/152',
!               '__2/123', '__2/623', '__2/152', '__2/652', '__3/123', '__3/623', '__3/134', '__3/634', '__4/134', '__4/634',
!               '__4/145', '__4/645', '__5/145', '__5/645', '__5/152', '__5/652', '__6/623', '__6/634', '__6/645', '__6/652',
!               '_12/123', '_12/152', '_13/123', '_13/134', '_14/134', '_14/145', '_15/145', '_15/152', '_62/623', '_62/652',
!               '_63/623', '_63/634', '_64/634', '_64/645', '_65/645', '_65/652', '_23/123', '_23/623', '_25/152', '_25/652',
!               '_43/134', '_43/634', '_45/145', '_45/645' /)
! --------------------------------------------------------------------------------------------
!
! = log ======================================================================================
!
! 30/01/17 - create
!
! ********************************************************************************************

subroutine find_neighbor_edge_3D(heavy_id, light_id, lgt_block, max_treelevel, dir, hvy_neighbor, lgt_active, lgt_n)

!---------------------------------------------------------------------------------------------
! modules

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
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    ! direction for neighbor search
    character(len=7), intent(in)        :: dir
    ! list of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_active(:)
    ! number of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_n

    ! heavy data array - neifghbor data
    integer(kind=ik), intent(out)       :: hvy_neighbor(:,:)

    ! auxiliary variables
    integer(kind=ik)                    :: list_id!, virt_code1, virt_list_id1, virt_code2, virt_list_id2, list_id2

    ! mesh level
    integer(kind=ik)                    :: level
    ! treecode varaibles
    integer(kind=ik)                    :: my_treecode(max_treelevel), neighbor(max_treelevel)!, virt_treecode(max_treelevel)
    ! return value from function "does_block_exist"
    logical                             :: exists
    ! neighbor light data id
    integer(kind=ik)                    :: neighbor_light_id

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    my_treecode     = lgt_block( light_id, 1:max_treelevel )
    level           = lgt_block( light_id, max_treelevel + 1 )

    list_id         = -1
!    virt_code1      = -1
!    virt_list_id1   = -1
!    virt_code2      = -1
!    virt_list_id2   = -1
!    list_id2        = -1

!---------------------------------------------------------------------------------------------
! main body

    ! set auxiliary variables
    select case(dir)

        case('_12/___')
            list_id    = 7

        case('_13/___')
            list_id    = 8

        case('_14/___')
            list_id    = 9

        case('_15/___')
            list_id    = 10

        case('_62/___')
            list_id    = 11

        case('_63/___')
            list_id    = 12

        case('_64/___')
            list_id    = 13

        case('_65/___')
            list_id    = 14

        case('_23/___')
            list_id    = 15

        case('_25/___')
            list_id    = 16

        case('_43/___')
            list_id    = 17

        case('_45/___')
            list_id    = 18

    end select

    ! calculate treecode for neighbor on same level
    call adjacent_block_3D( my_treecode, neighbor, dir, level, max_treelevel)

    ! proof existence of neighbor block and find light data id
    call does_block_exist(neighbor, lgt_block, max_treelevel, exists, neighbor_light_id, lgt_active, lgt_n)

    if (exists) then
        ! neighbor on same level
        hvy_neighbor( heavy_id, list_id ) = neighbor_light_id

    else
        ! to do

    end if

!    ! set virt_code and lvl_down_neighbor
!    select case(dir)
!        case('__N')
!            list_id    = 1
!            ! virtual treecodes, list_ids for neighbors on higher level
!            virt_code1    = 0
!            virt_list_id1 = 10
!            virt_code2    = 1
!            virt_list_id2 = 9
!            ! id2 for cases with neighbor one level down
!            if ( lgt_block(light_id, level) == 0) then
!                list_id2 = 10
!            elseif ( lgt_block(light_id, level) == 1) then
!                list_id2 = 9
!            end if
!
!        case('__E')
!            list_id = 2
!            ! virtual treecodes for neighbors on higher level
!            virt_code1    = 1
!            virt_list_id1 = 13
!            virt_code2    = 3
!            virt_list_id2 = 14
!            ! id2 for cases with neighbor one level down
!            if ( lgt_block(light_id, level) == 1) then
!                list_id2 = 13
!            elseif ( lgt_block(light_id, level) == 3) then
!                list_id2 = 14
!            end if
!
!        case('__S')
!            list_id   = 3
!            ! virtual treecodes for neighbors on higher level
!            virt_code1    = 2
!            virt_list_id1 = 12
!            virt_code2    = 3
!            virt_list_id2 = 11
!            ! id2 for cases with neighbor one level down
!            if ( lgt_block(light_id, level) == 3) then
!                list_id2 = 11
!            elseif ( lgt_block(light_id, level) == 2) then
!                list_id2 = 12
!            end if
!
!        case('__W')
!            list_id   = 4
!            ! virtual treecodes for neighbors on higher level
!            virt_code1    = 0
!            virt_list_id1 = 15
!            virt_code2    = 2
!            virt_list_id2 = 16
!            ! id2 for cases with neighbor one level down
!            if ( lgt_block(light_id, level) == 0) then
!                list_id2 = 15
!            elseif ( lgt_block(light_id, level) == 2) then
!                list_id2 = 16
!            end if
!
!    end select
!
!    ! calculate treecode for neighbor on same level
!    call adjacent_block_2D( my_treecode, neighbor, dir, level, max_treelevel)
!
!    ! proof existence of neighbor block and find light data id
!    call does_block_exist(neighbor, lgt_block, max_treelevel, exists, neighbor_light_id, lgt_active, lgt_n)
!
!    if (exists) then
!
!        ! neighbor on same level
!        ! write neighbor data, 2D: 16 possible neighbor relations
!        hvy_neighbor( heavy_id, list_id ) = neighbor_light_id
!
!    else
!
!        ! neighbor could be one level down
!        neighbor( level ) = -1
!        ! proof existence of neighbor block
!        call does_block_exist(neighbor, lgt_block, max_treelevel, exists, neighbor_light_id, lgt_active, lgt_n)
!        if ( exists ) then
!            ! neigbor is one level down
!            ! save list_id2
!            hvy_neighbor( heavy_id, list_id2 ) = neighbor_light_id
!
!        elseif ( .not.(exists) ) then
!            ! 2 neighbors one level up
!
!            ! first neighbor virtual treecode, one level up
!            virt_treecode = my_treecode
!            virt_treecode( level+1 ) = virt_code1
!
!            ! calculate treecode for neighbor on same level (virtual level)
!            call adjacent_block_2D( virt_treecode, neighbor, dir, level+1, max_treelevel)
!            ! proof existence of neighbor block
!            call does_block_exist(neighbor, lgt_block, max_treelevel, exists, neighbor_light_id, lgt_active, lgt_n)
!
!            if (exists) then
!                ! neigbor is one level up
!                ! write data
!                hvy_neighbor( heavy_id, virt_list_id1 ) = neighbor_light_id
!
!            else
!                ! error case
!                print*, my_treecode
!                print*, neighbor
!                print*, 'ERROR: can not find edge neighbor'
!                stop
!            end if
!
!            ! second neighbor virtual treecode, one level up
!            virt_treecode = my_treecode
!            virt_treecode( level+1 ) = virt_code2
!
!            ! calculate treecode for neighbor on same level (virtual level)
!            call adjacent_block_2D( virt_treecode, neighbor, dir, level+1, max_treelevel)
!            ! proof existence of neighbor block
!            call does_block_exist(neighbor, lgt_block, max_treelevel, exists, neighbor_light_id, lgt_active, lgt_n)
!
!            if (exists) then
!                ! neigbor is one level up
!                ! write data
!                hvy_neighbor( heavy_id, virt_list_id2 ) = neighbor_light_id
!
!            else
!                ! error case
!                print*, my_treecode
!                print*, neighbor
!                print*, 'ERROR: can not find edge neighbor'
!                stop
!            end if
!
!        end if
!
!    end if

end subroutine find_neighbor_edge_3D
