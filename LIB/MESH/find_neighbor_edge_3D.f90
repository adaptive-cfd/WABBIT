!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name find_neighbor_edge_3D.f90
!> \version 0.5
!> \author msr
!
!> \brief find neighbor on block edge
!
!>
!! input:
!!           - heavy and light data id
!!           - light data array and max treelevel
!!           - direction for neighbor search
!!           - list of active blocks
!!
!! output:
!!           - neighbor list array
!!
!! \n
!  --------------------------------------------------------------------------------------------
!> neighbor codes: \n
!  ---------------
!> for imagination:
!!                   - 6-sided dice with '1'-side on top, '6'-side on bottom, '2'-side in front
!!                   - edge: boundary between two sides - use sides numbers for coding
!!                   - corner: between three sides - so use all three sides numbers
!!                   - block on higher/lower level: block shares face/edge and one unique corner,
!!                     so use this corner code in second part of neighbor code
!!
!! faces:  '__1/___', '__2/___', '__3/___', '__4/___', '__5/___', '__6/___' \n
!! edges:  '_12/___', '_13/___', '_14/___', '_15/___'
!!         '_62/___', '_63/___', '_64/___', '_65/___'
!!         '_23/___', '_25/___', '_43/___', '_45/___' \n
!! corner: '123/___', '134/___', '145/___', '152/___'
!!         '623/___', '634/___', '645/___', '652/___' \n
!!
!! complete neighbor code array, 74 possible neighbor relations \n
!! \n
!! neighbors = (/'__1/___', '__2/___', '__3/___', '__4/___', '__5/___', '__6/___', '_12/___', '_13/___', '_14/___', '_15/___',
!!                '_62/___', '_63/___', '_64/___', '_65/___', '_23/___', '_25/___', '_43/___', '_45/___', '123/___', '134/___',
!!                '145/___', '152/___', '623/___', '634/___', '645/___', '652/___', '__1/123', '__1/134', '__1/145', '__1/152',
!!                '__2/123', '__2/623', '__2/152', '__2/652', '__3/123', '__3/623', '__3/134', '__3/634', '__4/134', '__4/634',
!!                '__4/145', '__4/645', '__5/145', '__5/645', '__5/152', '__5/652', '__6/623', '__6/634', '__6/645', '__6/652',
!!                '_12/123', '_12/152', '_13/123', '_13/134', '_14/134', '_14/145', '_15/145', '_15/152', '_62/623', '_62/652',
!!                '_63/623', '_63/634', '_64/634', '_64/645', '_65/645', '_65/652', '_23/123', '_23/623', '_25/152', '_25/652',
!!               '_43/134', '_43/634', '_45/145', '_45/645' /)
! --------------------------------------------------------------------------------------------
!
!> \details
!! = log ======================================================================================
!! \n
!! 30/01/17 - create
!
! ********************************************************************************************

subroutine find_neighbor_edge_3D(heavy_id, lgt_id, lgt_block, max_treelevel, dir, hvy_neighbor, lgt_n, lgt_sortednumlist)

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> heavy data id
    integer(kind=ik), intent(in)        :: heavy_id
    !> light data id
    integer(kind=ik), intent(in)        :: lgt_id
    !> max treelevel
    integer(kind=ik), intent(in)        :: max_treelevel
    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> direction for neighbor search
    character(len=7), intent(in)        :: dir
    !> number of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_n
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), intent(in)   :: lgt_sortednumlist(:,:)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(out)       :: hvy_neighbor(:,:)

    ! auxiliary variables
    integer(kind=ik)                    :: list_id, list_id2, virt_code(2), virt_list_id(2)

    ! mesh level
    integer(kind=ik)                    :: level
    ! treecode varaibles
    integer(kind=ik)                    :: my_treecode(max_treelevel), neighbor(max_treelevel), virt_treecode(max_treelevel)
    ! return value from function "does_block_exist"
    logical                             :: exists
    ! neighbor light data id
    integer(kind=ik)                    :: neighbor_light_id

    ! loop variable
    integer(kind=ik)                    :: k

    ! variable to show if there is a valid edge neighbor
    logical                             :: lvl_down_neighbor

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    my_treecode     = lgt_block( lgt_id, 1:max_treelevel )
    level           = lgt_block( lgt_id, max_treelevel + 1 )

    list_id         = -1
    list_id2        = -1
    virt_code       = -1
    virt_list_id    = -1

    lvl_down_neighbor = .false.

!---------------------------------------------------------------------------------------------
! main body

    ! set auxiliary variables
    select case(dir)

        case('_12/___')

            list_id    = 7

            ! id2 for cases with neighbor one level down
            if ( lgt_block(lgt_id, level) == 4) then
                list_id2 = 52
            elseif ( lgt_block(lgt_id, level) == 6) then
                list_id2 = 51
            end if

            ! virtual treecodes, list_ids for neighbors on higher level
            virt_code(1)    = 4
            virt_list_id(1) = 52
            virt_code(2)    = 6
            virt_list_id(2) = 51

            ! set logical for valid edge neighbors
            if ( (my_treecode( level ) == 4) .or. (my_treecode( level ) == 6) ) then
                lvl_down_neighbor = .true.
            end if

        case('_13/___')

            list_id    = 8

            ! id2 for cases with neighbor one level down
            if ( lgt_block(lgt_id, level) == 6) then
                list_id2 = 53
            elseif ( lgt_block(lgt_id, level) == 7) then
                list_id2 = 54
            end if

            ! virtual treecodes, list_ids for neighbors on higher level
            virt_code(1)    = 6
            virt_list_id(1) = 53
            virt_code(2)    = 7
            virt_list_id(2) = 54

            ! set logical for valid edge neighbors
            if ( (my_treecode( level ) == 6) .or. (my_treecode( level ) == 7) ) then
                lvl_down_neighbor = .true.
            end if

        case('_14/___')

            list_id    = 9

            ! id2 for cases with neighbor one level down
            if ( lgt_block(lgt_id, level) == 5) then
                list_id2 = 56
            elseif ( lgt_block(lgt_id, level) == 7) then
                list_id2 = 55
            end if

            ! virtual treecodes, list_ids for neighbors on higher level
            virt_code(1)    = 5
            virt_list_id(1) = 56
            virt_code(2)    = 7
            virt_list_id(2) = 55

            ! set logical for valid edge neighbors
            if ( (my_treecode( level ) == 5) .or. (my_treecode( level ) == 7) ) then
                lvl_down_neighbor = .true.
            end if

        case('_15/___')

            list_id    = 10

            ! id2 for cases with neighbor one level down
            if ( lgt_block(lgt_id, level) == 4) then
                list_id2 = 58
            elseif ( lgt_block(lgt_id, level) == 5) then
                list_id2 = 57
            end if

            ! virtual treecodes, list_ids for neighbors on higher level
            virt_code(1)    = 4
            virt_list_id(1) = 58
            virt_code(2)    = 5
            virt_list_id(2) = 57

            ! set logical for valid edge neighbors
            if ( (my_treecode( level ) == 4) .or. (my_treecode( level ) == 5) ) then
                lvl_down_neighbor = .true.
            end if

        case('_62/___')

            list_id    = 11

            ! id2 for cases with neighbor one level down
            if ( lgt_block(lgt_id, level) == 0) then
                list_id2 = 60
            elseif ( lgt_block(lgt_id, level) == 2) then
                list_id2 = 59
            end if

            ! virtual treecodes, list_ids for neighbors on higher level
            virt_code(1)    = 0
            virt_list_id(1) = 60
            virt_code(2)    = 2
            virt_list_id(2) = 59

            ! set logical for valid edge neighbors
            if ( (my_treecode( level ) == 0) .or. (my_treecode( level ) == 2) ) then
                lvl_down_neighbor = .true.
            end if

        case('_63/___')

            list_id    = 12

            ! id2 for cases with neighbor one level down
            if ( lgt_block(lgt_id, level) == 2) then
                list_id2 = 61
            elseif ( lgt_block(lgt_id, level) == 3) then
                list_id2 = 62
            end if

            ! virtual treecodes, list_ids for neighbors on higher level
            virt_code(1)    = 2
            virt_list_id(1) = 61
            virt_code(2)    = 3
            virt_list_id(2) = 62

            ! set logical for valid edge neighbors
            if ( (my_treecode( level ) == 2) .or. (my_treecode( level ) == 3) ) then
                lvl_down_neighbor = .true.
            end if

        case('_64/___')

            list_id    = 13

            ! id2 for cases with neighbor one level down
            if ( lgt_block(lgt_id, level) == 1) then
                list_id2 = 64
            elseif ( lgt_block(lgt_id, level) == 3) then
                list_id2 = 63
            end if

            ! virtual treecodes, list_ids for neighbors on higher level
            virt_code(1)    = 1
            virt_list_id(1) = 64
            virt_code(2)    = 3
            virt_list_id(2) = 63

            ! set logical for valid edge neighbors
            if ( (my_treecode( level ) == 1) .or. (my_treecode( level ) == 3) ) then
                lvl_down_neighbor = .true.
            end if

        case('_65/___')

            list_id    = 14

            ! id2 for cases with neighbor one level down
            if ( lgt_block(lgt_id, level) == 0) then
                list_id2 = 66
            elseif ( lgt_block(lgt_id, level) == 1) then
                list_id2 = 65
            end if

            ! virtual treecodes, list_ids for neighbors on higher level
            virt_code(1)    = 0
            virt_list_id(1) = 66
            virt_code(2)    = 1
            virt_list_id(2) = 65

            ! set logical for valid edge neighbors
            if ( (my_treecode( level ) == 0) .or. (my_treecode( level ) == 1) ) then
                lvl_down_neighbor = .true.
            end if

        case('_23/___')

            list_id    = 15

            ! id2 for cases with neighbor one level down
            if ( lgt_block(lgt_id, level) == 2) then
                list_id2 = 68
            elseif ( lgt_block(lgt_id, level) == 6) then
                list_id2 = 67
            end if

            ! virtual treecodes, list_ids for neighbors on higher level
            virt_code(1)    = 2
            virt_list_id(1) = 68
            virt_code(2)    = 6
            virt_list_id(2) = 67

            ! set logical for valid edge neighbors
            if ( (my_treecode( level ) == 2) .or. (my_treecode( level ) == 6) ) then
                lvl_down_neighbor = .true.
            end if

        case('_25/___')

            list_id    = 16

            ! id2 for cases with neighbor one level down
            if ( lgt_block(lgt_id, level) == 0) then
                list_id2 = 70
            elseif ( lgt_block(lgt_id, level) == 4) then
                list_id2 = 69
            end if

            ! virtual treecodes, list_ids for neighbors on higher level
            virt_code(1)    = 0
            virt_list_id(1) = 70
            virt_code(2)    = 4
            virt_list_id(2) = 69

            ! set logical for valid edge neighbors
            if ( (my_treecode( level ) == 0) .or. (my_treecode( level ) == 4) ) then
                lvl_down_neighbor = .true.
            end if

        case('_43/___')

            list_id    = 17

            ! id2 for cases with neighbor one level down
            if ( lgt_block(lgt_id, level) == 3) then
                list_id2 = 72
            elseif ( lgt_block(lgt_id, level) == 7) then
                list_id2 = 71
            end if

            ! virtual treecodes, list_ids for neighbors on higher level
            virt_code(1)    = 3
            virt_list_id(1) = 72
            virt_code(2)    = 7
            virt_list_id(2) = 71

            ! set logical for valid edge neighbors
            if ( (my_treecode( level ) == 3) .or. (my_treecode( level ) == 7) ) then
                lvl_down_neighbor = .true.
            end if

        case('_45/___')

            list_id    = 18

            ! id2 for cases with neighbor one level down
            if ( lgt_block(lgt_id, level) == 1) then
                list_id2 = 74
            elseif ( lgt_block(lgt_id, level) == 5) then
                list_id2 = 73
            end if

            ! virtual treecodes, list_ids for neighbors on higher level
            virt_code(1)    = 1
            virt_list_id(1) = 74
            virt_code(2)    = 5
            virt_list_id(2) = 73

            ! set logical for valid edge neighbors
            if ( (my_treecode( level ) == 1) .or. (my_treecode( level ) == 5) ) then
                lvl_down_neighbor = .true.
            end if

    end select

    ! calculate treecode for neighbor on same level
    call adjacent_block_3D( my_treecode, neighbor, dir, level, max_treelevel)

    ! proof existence of neighbor block and find light data id
    call does_block_exist(neighbor, exists, neighbor_light_id, lgt_sortednumlist, lgt_n)

    if (exists) then
        ! neighbor on same level
        hvy_neighbor( heavy_id, list_id ) = neighbor_light_id

    else

        ! neighbor could be one level down
        neighbor( level ) = -1
        ! proof existence of neighbor block
        call does_block_exist(neighbor, exists, neighbor_light_id, lgt_sortednumlist, lgt_n)
        if ( exists .and. lvl_down_neighbor ) then

            ! neigbor is one level down
            ! save list_id2
            hvy_neighbor( heavy_id, list_id2 ) = neighbor_light_id

        elseif ( .not.(exists) ) then
            ! 2 neighbors one level up
            ! loop over all 2 possible neighbors

            do k = 1, 2

                ! first neighbor virtual treecode, one level up
                virt_treecode = my_treecode
                virt_treecode( level+1 ) = virt_code(k)

                ! calculate treecode for neighbor on same level (virtual level)
                call adjacent_block_3D( virt_treecode, neighbor, dir, level+1, max_treelevel)
                ! proof existence of neighbor block
                call does_block_exist(neighbor, exists, neighbor_light_id, lgt_sortednumlist, lgt_n)
                if (exists) then
                    ! neigbor is one level up
                    ! write data
                    hvy_neighbor( heavy_id, virt_list_id(k) ) = neighbor_light_id

                else
                    ! error case
                    print*, dir
                    print*, my_treecode
                    print*, neighbor
                    print*, 'ERROR: can not find edge neighbor'
                    stop
                end if

            end do

        end if

    end if

end subroutine find_neighbor_edge_3D
