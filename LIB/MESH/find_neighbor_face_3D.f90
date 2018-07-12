
!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name find_neighbor_face_3D.f90
!> \version 0.5
!> \author msr
!
!> \brief find neighbor on block face \n
!
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
!! \image html neighborcode.svg "Neighborcode 3D" width=250
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

subroutine find_neighbor_face_3D(params, heavy_id, lgt_id, lgt_block, max_treelevel, dir, hvy_neighbor, &
    lgt_n, lgt_sortednumlist, error)

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none
    !> user defined parameter structure
    type (type_params), intent(in)      :: params
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
    integer(kind=tsize), intent(in)     :: lgt_sortednumlist(:,:)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(inout)     :: hvy_neighbor(:,:)
    logical, intent(inout)              :: error

    ! auxiliary variables
    integer(kind=ik)                    :: neighborID_sameLevel
    integer(kind=ik)                    :: neighborID_coarserLevel, virt_code(4)
    integer(kind=ik)                    :: neighborID_finerLevel(4)

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

    neighborID_sameLevel     = -1
    neighborID_coarserLevel  = -1
    virt_code                = -1
    neighborID_finerLevel    = -1

!---------------------------------------------------------------------------------------------
! main body

    ! set auxiliary variables
    select case(dir)

        case('__1/___')

            neighborID_sameLevel    = 1

            ! If the neighbor is coarser, then we have only one possible block, but
            ! the finer block (me) may be at four positions, which define the neighborhood
            if ( lgt_block(lgt_id, level) == 4) then
                neighborID_coarserLevel = 30
            elseif ( lgt_block(lgt_id, level) == 5) then
                neighborID_coarserLevel = 29
            elseif ( lgt_block(lgt_id, level) == 6) then
                neighborID_coarserLevel = 27
            elseif ( lgt_block(lgt_id, level) == 7) then
                neighborID_coarserLevel = 28
            end if

            ! virtual treecodes, list_ids for neighbors on higher level
            virt_code(1)    = 4
            neighborID_finerLevel(1) = 30
            virt_code(2)    = 5
            neighborID_finerLevel(2) = 29
            virt_code(3)    = 6
            neighborID_finerLevel(3) = 27
            virt_code(4)    = 7
            neighborID_finerLevel(4) = 28

        case('__2/___')

            neighborID_sameLevel    = 2

            ! If the neighbor is coarser, then we have only one possible block, but
            ! the finer block (me) may be at four positions, which define the neighborhood
            if ( lgt_block(lgt_id, level) == 0) then
                neighborID_coarserLevel = 34
            elseif ( lgt_block(lgt_id, level) == 2) then
                neighborID_coarserLevel = 32
            elseif ( lgt_block(lgt_id, level) == 4) then
                neighborID_coarserLevel = 33
            elseif ( lgt_block(lgt_id, level) == 6) then
                neighborID_coarserLevel = 31
            end if

            ! virtual treecodes, list_ids for neighbors on higher level
            virt_code(1)    = 0
            neighborID_finerLevel(1) = 34
            virt_code(2)    = 2
            neighborID_finerLevel(2) = 32
            virt_code(3)    = 4
            neighborID_finerLevel(3) = 33
            virt_code(4)    = 6
            neighborID_finerLevel(4) = 31

        case('__3/___')

            neighborID_sameLevel    = 3

            ! If the neighbor is coarser, then we have only one possible block, but
            ! the finer block (me) may be at four positions, which define the neighborhood
            if ( lgt_block(lgt_id, level) == 2) then
                neighborID_coarserLevel = 36
            elseif ( lgt_block(lgt_id, level) == 3) then
                neighborID_coarserLevel = 38
            elseif ( lgt_block(lgt_id, level) == 6) then
                neighborID_coarserLevel = 35
            elseif ( lgt_block(lgt_id, level) == 7) then
                neighborID_coarserLevel = 37
            end if

            ! virtual treecodes, list_ids for neighbors on higher level
            virt_code(1)    = 2
            neighborID_finerLevel(1) = 36
            virt_code(2)    = 3
            neighborID_finerLevel(2) = 38
            virt_code(3)    = 6
            neighborID_finerLevel(3) = 35
            virt_code(4)    = 7
            neighborID_finerLevel(4) = 37

        case('__4/___')

            neighborID_sameLevel    = 4

            ! If the neighbor is coarser, then we have only one possible block, but
            ! the finer block (me) may be at four positions, which define the neighborhood
            if ( lgt_block(lgt_id, level) == 1) then
                neighborID_coarserLevel = 42
            elseif ( lgt_block(lgt_id, level) == 3) then
                neighborID_coarserLevel = 40
            elseif ( lgt_block(lgt_id, level) == 5) then
                neighborID_coarserLevel = 41
            elseif ( lgt_block(lgt_id, level) == 7) then
                neighborID_coarserLevel = 39
            end if

            ! virtual treecodes, list_ids for neighbors on higher level
            virt_code(1)    = 1
            neighborID_finerLevel(1) = 42
            virt_code(2)    = 3
            neighborID_finerLevel(2) = 40
            virt_code(3)    = 5
            neighborID_finerLevel(3) = 41
            virt_code(4)    = 7
            neighborID_finerLevel(4) = 39

        case('__5/___')

            neighborID_sameLevel    = 5

            ! If the neighbor is coarser, then we have only one possible block, but
            ! the finer block (me) may be at four positions, which define the neighborhood
            if ( lgt_block(lgt_id, level) == 0) then
                neighborID_coarserLevel = 46
            elseif ( lgt_block(lgt_id, level) == 1) then
                neighborID_coarserLevel = 44
            elseif ( lgt_block(lgt_id, level) == 4) then
                neighborID_coarserLevel = 45
            elseif ( lgt_block(lgt_id, level) == 5) then
                neighborID_coarserLevel = 43
            end if

            ! virtual treecodes, list_ids for neighbors on higher level
            virt_code(1)    = 0
            neighborID_finerLevel(1) = 46
            virt_code(2)    = 1
            neighborID_finerLevel(2) = 44
            virt_code(3)    = 4
            neighborID_finerLevel(3) = 45
            virt_code(4)    = 5
            neighborID_finerLevel(4) = 43

        case('__6/___')

            neighborID_sameLevel    = 6

            ! If the neighbor is coarser, then we have only one possible block, but
            ! the finer block (me) may be at four positions, which define the neighborhood
            if ( lgt_block(lgt_id, level) == 0) then
                neighborID_coarserLevel = 50
            elseif ( lgt_block(lgt_id, level) == 1) then
                neighborID_coarserLevel = 49
            elseif ( lgt_block(lgt_id, level) == 2) then
                neighborID_coarserLevel = 47
            elseif ( lgt_block(lgt_id, level) == 3) then
                neighborID_coarserLevel = 48
            end if

            ! virtual treecodes, list_ids for neighbors on higher level
            virt_code(1)    = 0
            neighborID_finerLevel(1) = 50
            virt_code(2)    = 1
            neighborID_finerLevel(2) = 49
            virt_code(3)    = 2
            neighborID_finerLevel(3) = 47
            virt_code(4)    = 3
            neighborID_finerLevel(4) = 48

        case default
            call abort(636300, "well you mustnt")

    end select

    ! calculate treecode for neighbor on same level
    call adjacent_block_3D( my_treecode, neighbor, dir, level, max_treelevel)
    ! check existence of neighbor block and find light data id
    call does_block_exist(neighbor, exists, neighbor_light_id, lgt_sortednumlist, lgt_n)

    if (exists) then
        ! neighbor on same level
        hvy_neighbor( heavy_id, neighborID_sameLevel ) = neighbor_light_id
    else

        ! We did not find the neighbor on the same level, and now check on coarser levels.
        ! Depending on my own treecode, I know what neigbor I am looking for, and I just set
        ! the last index to -1 = I go one level down (coarser)
        neighbor( level ) = -1

        ! check existence of neighbor block
        call does_block_exist(neighbor, exists, neighbor_light_id, lgt_sortednumlist, lgt_n)
        if ( exists ) then
            ! neighbor is one level down (coarser)
            hvy_neighbor( heavy_id, neighborID_coarserLevel ) = neighbor_light_id

        elseif ( .not.(exists) ) then
            ! 4 neighbors one level up
            ! loop over all 4 possible neighbors

            do k = 1, 4

                ! first neighbor virtual treecode, one level up
                virt_treecode = my_treecode
                virt_treecode( level+1 ) = virt_code(k)

                ! calculate treecode for neighbor on same level (virtual level)
                call adjacent_block_3D( virt_treecode, neighbor, dir, level+1, max_treelevel)
                ! check existence of neighbor block
                call does_block_exist(neighbor, exists, neighbor_light_id, lgt_sortednumlist, lgt_n)

                if (exists) then
                    ! neigbor is one level up
                    ! write data
                    hvy_neighbor( heavy_id, neighborID_finerLevel(k) ) = neighbor_light_id
                else
                    ! error case
                    write(*,*) "find_neighbor_face_3D: my treecode", my_treecode, "dir", dir, "neighbor treecode", neighbor
                    error = .true.
                end if

            end do

        end if

    end if

end subroutine find_neighbor_face_3D
