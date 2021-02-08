!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name find_neighbor_corner_2D.f90
!> \version 0.4
!> \author msr
!
!> \brief find neighbor on block corner
!> \details  valid cases for corner neighbors:
!!    1. same level: always exact one neighbor
!!    2. one level down: one neighbor, if this neighbor
!! only on block corner (not additional on block side)
!!    3. one level up: always exact one neighbor
!!
!! input:
!!           - heavy and light data id
!!           - light data array and max treelevel
!!           - direction for neighbor search
!!           - list of active blocks
!!
!!  output:
!!           - neighbor list array
!!
! -------------------------------------------------------------------------------------------------------------------------
!>  dirs = (/'__N', '__E', '__S', '__W', '_NE', '_NW', '_SE', '_SW', 'NNE', 'NNW', 'SSE', 'SSW', 'ENE', 'ESE', 'WNW', 'WSW'/)
! 1 __N
! 2 __E
! 3 __S
! 4 __W
! 5 _NE
! 6 _NW
! 7 _SE
! 8 _SW
! 9 NNE
! 10 NNW
! 11 SSE
! 12 SSW
! 13 ENE
! 14 ESE
! 15 WNW
! 16 WSW
! -------------------------------------------------------------------------------------------------------------------------
!> \details
!! = log ======================================================================================
!! \n
!! 08/11/16 - switch to v0.4
! ********************************************************************************************
!> \image html neighborhood.svg "Neighborhood Relations in 2D" width=400

subroutine find_neighbor_corner_2D(params, heavy_id, light_id, lgt_block, max_treelevel, dir, &
    hvy_neighbor, lgt_n, lgt_sortednumlist, error, n_domain)

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> heavy data id
    integer(kind=ik), intent(in)        :: heavy_id
    !> light data id
    integer(kind=ik), intent(in)        :: light_id
    !> max treelevel
    integer(kind=ik), intent(in)        :: max_treelevel
    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> direction for neighbor search
    character(len=3), intent(in)        :: dir
    !> number of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_n
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), intent(in)     :: lgt_sortednumlist(:,:)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(inout)     :: hvy_neighbor(:,:)
    logical, intent(inout)              :: error
    integer(kind=2), intent(in) :: n_domain(1:3)

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
    integer(kind=ik)                    :: neighbor_light_id, tree_id

    my_treecode = lgt_block( light_id, 1:max_treelevel )
    level       = lgt_block( light_id, max_treelevel + IDX_MESH_LVL )
    tree_id     = lgt_block( light_id, max_treelevel + IDX_TREE_ID )
    lvl_down_neighbor = .false.
    virt_code = -1
    list_id   = -1


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

    !---------------------------------------------------------------------------
    ! Does the neighbor on the same level exist ?
    !---------------------------------------------------------------------------
    ! calculate treecode for neighbor on same level
    call adjacent_block_2D( my_treecode, neighbor, dir, level, max_treelevel)
    ! check existence of neighbor block and find light data id
    call does_block_exist(neighbor, exists, neighbor_light_id, lgt_sortednumlist, lgt_n, tree_id)

    if (exists) then
        ! we found the neighbor on same level in the direction "dir", so we're done.
        hvy_neighbor( heavy_id, list_id ) = neighbor_light_id
    else
        ! neighbor could be one level down
        neighbor( level ) = -1
        ! check existence of neighbor block
        call does_block_exist(neighbor, exists, neighbor_light_id, lgt_sortednumlist, lgt_n, tree_id)

        if ( exists .and. lvl_down_neighbor ) then
            ! neigbor is one level down
            hvy_neighbor( heavy_id, list_id ) = neighbor_light_id

        elseif ( .not.(exists) .and. (level/=params%max_treelevel)) then
            ! neighbor could be on level up
            ! virtual treecode, one level up
            virt_treecode = my_treecode
            virt_treecode( level+1 ) = virt_code

            ! calculate treecode for neighbor on same level (virtual level)
            call adjacent_block_2D( virt_treecode, neighbor, dir, level+1, max_treelevel)
            ! check existence of neighbor block
            call does_block_exist(neighbor, exists, neighbor_light_id, lgt_sortednumlist, lgt_n, tree_id)

            if (exists) then
                ! neigbor is one level up
                hvy_neighbor( heavy_id, list_id ) = neighbor_light_id
            end if

            ! we did not find a neighbor. that may be a bad grid error, or simply, there is none
            ! because symmetry conditions are used.
            if ((.not. exists .and. ALL(params%periodic_BC)).or.(maxval(abs(n_domain))==0.and..not.exists)) then
                write(*,*) "Corner neighbor not found!"
                error = .true.
            endif

        end if

    end if



end subroutine find_neighbor_corner_2D
