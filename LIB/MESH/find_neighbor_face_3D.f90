
!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name find_neighbor_face_3D.f90
!> \version 0.5
!> \author msr
!
!> \brief find neighbor on block face
!
!> \details
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
!!
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
!!faces:  '__1/___', '__2/___', '__3/___', '__4/___', '__5/___', '__6/___' \n
!! edges:  '_12/___', '_13/___', '_14/___', '_15/___' 
!!         '_62/___', '_63/___', '_64/___', '_65/___' 
!!         '_23/___', '_25/___', '_43/___', '_45/___' \n
!! corner: '123/___', '134/___', '145/___', '152/___' 
!!         '623/___', '634/___', '645/___', '652/___' \n
!! \n
!! complete neighbor code array, 74 possible neighbor relations \n
!! neighbors = (/'__1/___', '__2/___', '__3/___', '__4/___', '__5/___', '__6/___', '_12/___', '_13/___', '_14/___', '_15/___',
!!               '_62/___', '_63/___', '_64/___', '_65/___', '_23/___', '_25/___', '_43/___', '_45/___', '123/___', '134/___',
!!               '145/___', '152/___', '623/___', '634/___', '645/___', '652/___', '__1/123', '__1/134', '__1/145', '__1/152',
!!               '__2/123', '__2/623', '__2/152', '__2/652', '__3/123', '__3/623', '__3/134', '__3/634', '__4/134', '__4/634',
!!               '__4/145', '__4/645', '__5/145', '__5/645', '__5/152', '__5/652', '__6/623', '__6/634', '__6/645', '__6/652',
!!               '_12/123', '_12/152', '_13/123', '_13/134', '_14/134', '_14/145', '_15/145', '_15/152', '_62/623', '_62/652',
!!               '_63/623', '_63/634', '_64/634', '_64/645', '_65/645', '_65/652', '_23/123', '_23/623', '_25/152', '_25/652',
!!               '_43/134', '_43/634', '_45/145', '_45/645' /) \n
!
! --------------------------------------------------------------------------------------------
!
!> \details
!! = log ======================================================================================
!! \n
!! 27/01/17 - start
!
! ********************************************************************************************

subroutine find_neighbor_face_3D(heavy_id, lgt_id, lgt_block, max_treelevel, dir, hvy_neighbor, lgt_active, lgt_n)

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
    !> list of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_n

    !> heavy data array - neifghbor data
    integer(kind=ik), intent(out)       :: hvy_neighbor(:,:)

    ! auxiliary variables
    integer(kind=ik)                    :: list_id, list_id2, virt_code(4), virt_list_id(4)

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

!---------------------------------------------------------------------------------------------
! main body

    ! set auxiliary variables
    select case(dir)

        case('__1/___')

            list_id    = 1

            ! id2 for cases with neighbor one level down
            if ( lgt_block(lgt_id, level) == 4) then
                list_id2 = 30
            elseif ( lgt_block(lgt_id, level) == 5) then
                list_id2 = 29
            elseif ( lgt_block(lgt_id, level) == 6) then
                list_id2 = 27
            elseif ( lgt_block(lgt_id, level) == 7) then
                list_id2 = 28
            end if

            ! virtual treecodes, list_ids for neighbors on higher level
            virt_code(1)    = 4
            virt_list_id(1) = 30
            virt_code(2)    = 5
            virt_list_id(2) = 29
            virt_code(3)    = 6
            virt_list_id(3) = 27
            virt_code(4)    = 7
            virt_list_id(4) = 28

        case('__2/___')

            list_id    = 2

            ! id2 for cases with neighbor one level down
            if ( lgt_block(lgt_id, level) == 0) then
                list_id2 = 34
            elseif ( lgt_block(lgt_id, level) == 2) then
                list_id2 = 32
            elseif ( lgt_block(lgt_id, level) == 4) then
                list_id2 = 33
            elseif ( lgt_block(lgt_id, level) == 6) then
                list_id2 = 31
            end if

            ! virtual treecodes, list_ids for neighbors on higher level
            virt_code(1)    = 0
            virt_list_id(1) = 34
            virt_code(2)    = 2
            virt_list_id(2) = 32
            virt_code(3)    = 4
            virt_list_id(3) = 33
            virt_code(4)    = 6
            virt_list_id(4) = 31

        case('__3/___')

            list_id    = 3

            ! id2 for cases with neighbor one level down
            if ( lgt_block(lgt_id, level) == 2) then
                list_id2 = 36
            elseif ( lgt_block(lgt_id, level) == 3) then
                list_id2 = 38
            elseif ( lgt_block(lgt_id, level) == 6) then
                list_id2 = 35
            elseif ( lgt_block(lgt_id, level) == 7) then
                list_id2 = 37
            end if

            ! virtual treecodes, list_ids for neighbors on higher level
            virt_code(1)    = 2
            virt_list_id(1) = 36
            virt_code(2)    = 3
            virt_list_id(2) = 38
            virt_code(3)    = 6
            virt_list_id(3) = 35
            virt_code(4)    = 7
            virt_list_id(4) = 37

        case('__4/___')

            list_id    = 4

            ! id2 for cases with neighbor one level down
            if ( lgt_block(lgt_id, level) == 1) then
                list_id2 = 42
            elseif ( lgt_block(lgt_id, level) == 3) then
                list_id2 = 40
            elseif ( lgt_block(lgt_id, level) == 5) then
                list_id2 = 41
            elseif ( lgt_block(lgt_id, level) == 7) then
                list_id2 = 39
            end if

            ! virtual treecodes, list_ids for neighbors on higher level
            virt_code(1)    = 1
            virt_list_id(1) = 42
            virt_code(2)    = 3
            virt_list_id(2) = 40
            virt_code(3)    = 5
            virt_list_id(3) = 41
            virt_code(4)    = 7
            virt_list_id(4) = 39

        case('__5/___')

            list_id    = 5

            ! id2 for cases with neighbor one level down
            if ( lgt_block(lgt_id, level) == 0) then
                list_id2 = 46
            elseif ( lgt_block(lgt_id, level) == 1) then
                list_id2 = 44
            elseif ( lgt_block(lgt_id, level) == 4) then
                list_id2 = 45
            elseif ( lgt_block(lgt_id, level) == 5) then
                list_id2 = 43
            end if

            ! virtual treecodes, list_ids for neighbors on higher level
            virt_code(1)    = 0
            virt_list_id(1) = 46
            virt_code(2)    = 1
            virt_list_id(2) = 44
            virt_code(3)    = 4
            virt_list_id(3) = 45
            virt_code(4)    = 5
            virt_list_id(4) = 43

        case('__6/___')

            list_id    = 6

            ! id2 for cases with neighbor one level down
            if ( lgt_block(lgt_id, level) == 0) then
                list_id2 = 50
            elseif ( lgt_block(lgt_id, level) == 1) then
                list_id2 = 49
            elseif ( lgt_block(lgt_id, level) == 2) then
                list_id2 = 47
            elseif ( lgt_block(lgt_id, level) == 3) then
                list_id2 = 48
            end if

            ! virtual treecodes, list_ids for neighbors on higher level
            virt_code(1)    = 0
            virt_list_id(1) = 50
            virt_code(2)    = 1
            virt_list_id(2) = 49
            virt_code(3)    = 2
            virt_list_id(3) = 47
            virt_code(4)    = 3
            virt_list_id(4) = 48

    end select

    ! calculate treecode for neighbor on same level
    call adjacent_block_3D( my_treecode, neighbor, dir, level, max_treelevel)

    ! proof existence of neighbor block and find light data id
    call does_block_exist(neighbor, lgt_block, max_treelevel, exists, neighbor_light_id, lgt_active, lgt_n)

    if (exists) then
        ! neighbor on same level
        hvy_neighbor( heavy_id, list_id ) = neighbor_light_id

    else

        ! neighbor could be one level down
        neighbor( level ) = -1
        ! proof existence of neighbor block
        call does_block_exist(neighbor, lgt_block, max_treelevel, exists, neighbor_light_id, lgt_active, lgt_n)
        if ( exists ) then
            ! neigbor is one level down
            ! save list_id2
            hvy_neighbor( heavy_id, list_id2 ) = neighbor_light_id

        elseif ( .not.(exists) ) then
            ! 4 neighbors one level up
            ! loop over all 4 possible neighbors

            do k = 1, 4

                ! first neighbor virtual treecode, one level up
                virt_treecode = my_treecode
                virt_treecode( level+1 ) = virt_code(k)

                ! calculate treecode for neighbor on same level (virtual level)
                call adjacent_block_3D( virt_treecode, neighbor, dir, level+1, max_treelevel)
                ! proof existence of neighbor block
                call does_block_exist(neighbor, lgt_block, max_treelevel, exists, neighbor_light_id, lgt_active, lgt_n)

                if (exists) then
                    ! neigbor is one level up
                    ! write data
                    hvy_neighbor( heavy_id, virt_list_id(k) ) = neighbor_light_id

                else
                    ! error case
                    print*, dir
                    print*, my_treecode
                    print*, neighbor
                    print*, 'ERROR: can not find face neighbor'
                    stop
                end if

            end do

        end if

    end if

end subroutine find_neighbor_face_3D
