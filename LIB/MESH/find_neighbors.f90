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
! ********************************************************************************************
subroutine find_neighbor(params, hvyID_block, lgtID_block, Jmax, dir, error, n_domain)

    implicit none
    type (type_params), intent(in)      :: params                   !< user defined parameter structure
    integer(kind=ik), intent(in)        :: hvyID_block
    integer(kind=ik), intent(in)        :: lgtID_block
    integer(kind=ik), intent(in)        :: Jmax
    !> direction for neighbor search - number where each digit represents a cardinal direction
    !> 652 -> first 6 (bottom, z-1), then 5 (north, x-1) then 2 (front, y-1) 
    integer(kind=ik), intent(in)        :: dir                      
    logical, intent(inout)              :: error
    integer(kind=2), intent(in)         :: n_domain(1:3)
    integer(kind=ik)                    :: neighborDirCode_sameLevel
    integer(kind=ik)                    :: neighborDirCode_coarserLevel, tcFinerAppendDigit(4)
    integer(kind=ik)                    :: neighborDirCode_finerLevel(4)
    integer(kind=ik)                    :: k, lgtID_neighbor, level, tc_last, tree_ID
    integer(kind=tsize)                 :: tcb_Block, tcb_Neighbor, tcb_Virtual
    logical                             :: exists
    ! variable to show if there is a valid edge neighbor
    logical                             :: lvl_down_neighbor
    logical :: thereMustBeANeighbor

    level      = lgt_block( lgtID_block, IDX_MESH_LVL )
    tree_ID    = lgt_block( lgtID_block, IDX_TREE_ID )
    neighborDirCode_sameLevel     = -1
    neighborDirCode_coarserLevel  = -1
    neighborDirCode_finerLevel    = -1
    tcFinerAppendDigit            = -1

    ! we have to init tcBlock before we set it elsewise fortran doesnt like it
    tcb_Block    = get_tc(lgt_block(lgtID_block, IDX_TC_1 : IDX_TC_2))
    ! last digit is used very often so we only extract it once
    tc_last = tc_get_digit_at_level_b(tcb_Block, dim=params%dim, level=level, max_level=Jmax)

    ! not all blocks can have coarse neighbors.
    ! Consider:
    ! a c E E
    ! b d E E
    ! Then to the right, block b cannot have a coarser neighbor
    lvl_down_neighbor = .false.
    ! For the faces, we insist on finding neighbors (exception: symmetry conditions).
    ! for corners and edges, we do not- in the above example, block d does not
    ! have a neighbor in the top right corner.
    thereMustBeANeighbor = .false.


    ! 2D:
    !   faces: 1. same level: always one neighbor
    !          2. coarser: one neighbor, two possible neighbor codes
    !          3. finer: always two neighbors
    !   edges: 1. same level: always one neighbor
    !          2. coarser: one neighbor, possibly. Exists only if the coarser blocks corner coincides with the finer blocks corner.
    !          3. finer: always two neighbors




    ! set auxiliary variables
    !> direction for neighbor search - number where each digit represents a cardinal direction
    !> 652 -> first 6 (bottom, z-1), then 5 (north, x-1) then 2 (front, y-1) 
    select case(dir)
        case(1)  ! '__1/___'
            neighborDirCode_sameLevel    = 1
            thereMustBeANeighbor = .true.

            ! If the neighbor is coarser, then we have only one possible block, but
            ! the finer block (me) may be at four positions, which define the neighborhood code
            if ( tc_last == 4) then
                neighborDirCode_coarserLevel = 30
            elseif ( tc_last == 5) then
                neighborDirCode_coarserLevel = 29
            elseif ( tc_last == 6) then
                neighborDirCode_coarserLevel = 27
            elseif ( tc_last == 7) then
                neighborDirCode_coarserLevel = 28
            end if
            lvl_down_neighbor = .true.

            ! virtual treecodes, list_ids for neighbors on higher level
            tcFinerAppendDigit(1:4)  = (/ 4, 5, 6, 7 /)
            neighborDirCode_finerLevel(1:4) = (/ 30, 29, 27, 28 /)

        case(2)  ! '__2/___'
            if (params%dim == 2) then
                neighborDirCode_sameLevel = 4
                lvl_down_neighbor = .true.
                thereMustBeANeighbor = .true.
                ! virtual treecodes for neighbors on higher level
                tcFinerAppendDigit(1:2)         = (/ 0, 2 /)
                neighborDirCode_finerLevel(1:2) = (/ 15, 16 /)
                ! neighbor code for coarser neighbors
                if ( tc_last == 0) then
                    neighborDirCode_coarserLevel = 15
                elseif ( tc_last == 2) then
                    neighborDirCode_coarserLevel = 16
                end if

            else
                neighborDirCode_sameLevel    = 2
                thereMustBeANeighbor = .true.

                ! If the neighbor is coarser, then we have only one possible block, but
                ! the finer block (me) may be at four positions, which define the neighborhood code
                if ( tc_last == 0) then
                    neighborDirCode_coarserLevel = 34
                elseif ( tc_last == 2) then
                    neighborDirCode_coarserLevel = 32
                elseif ( tc_last == 4) then
                    neighborDirCode_coarserLevel = 33
                elseif ( tc_last == 6) then
                    neighborDirCode_coarserLevel = 31
                end if
                lvl_down_neighbor = .true.

                ! virtual treecodes, list_ids for neighbors on higher level
                tcFinerAppendDigit(1:4)  = (/ 0, 2, 4, 6 /)
                neighborDirCode_finerLevel(1:4) = (/ 34, 32, 33, 31 /)

            endif

        case(3)  ! '__3/___'  or  '__S'
            if (params%dim == 2) then
                neighborDirCode_sameLevel = 3
                lvl_down_neighbor = .true.
                thereMustBeANeighbor = .true.
                ! virtual treecodes for neighbors on higher level
                tcFinerAppendDigit(1:2)         = (/ 2, 3 /)
                neighborDirCode_finerLevel(1:2) = (/ 12, 11 /)
                ! neighbor code for coarser neighbors
                if ( tc_last == 3) then
                    neighborDirCode_coarserLevel = 11
                elseif ( tc_last == 2) then
                    neighborDirCode_coarserLevel = 12
                end if

            else
                neighborDirCode_sameLevel    = 3
                thereMustBeANeighbor = .true.

                ! If the neighbor is coarser, then we have only one possible block, but
                ! the finer block (me) may be at four positions, which define the neighborhood code
                if ( tc_last == 2) then
                    neighborDirCode_coarserLevel = 36
                elseif ( tc_last == 3) then
                    neighborDirCode_coarserLevel = 38
                elseif ( tc_last == 6) then
                    neighborDirCode_coarserLevel = 35
                elseif ( tc_last == 7) then
                    neighborDirCode_coarserLevel = 37
                end if
                lvl_down_neighbor = .true.

                ! virtual treecodes, list_ids for neighbors on higher level
                tcFinerAppendDigit(1:4)  = (/ 2, 3, 6, 7 /)
                neighborDirCode_finerLevel(1:4) = (/ 36, 38, 35, 37 /)

            endif

        case(4)  ! '__4/___'  or  '__E'
            if (params%dim == 2) then
                neighborDirCode_sameLevel = 2
                lvl_down_neighbor = .true.
                thereMustBeANeighbor = .true.
                ! virtual treecodes for neighbors on higher level
                tcFinerAppendDigit(1:2)         = (/ 1, 3 /)
                neighborDirCode_finerLevel(1:2) = (/ 13, 14 /)
                ! neighbor code for coarser neighbors
                if ( tc_last == 1) then
                    neighborDirCode_coarserLevel = 13
                elseif ( tc_last == 3) then
                    neighborDirCode_coarserLevel = 14
                end if

            else
                neighborDirCode_sameLevel    = 4
                thereMustBeANeighbor = .true.

                ! If the neighbor is coarser, then we have only one possible block, but
                ! the finer block (me) may be at four positions, which define the neighborhood code
                if ( tc_last == 1) then
                    neighborDirCode_coarserLevel = 42
                elseif ( tc_last == 3) then
                    neighborDirCode_coarserLevel = 40
                elseif ( tc_last == 5) then
                    neighborDirCode_coarserLevel = 41
                elseif ( tc_last == 7) then
                    neighborDirCode_coarserLevel = 39
                end if
                lvl_down_neighbor = .true.

                ! virtual treecodes, list_ids for neighbors on higher level
                tcFinerAppendDigit(1:4)  = (/ 1, 3, 5, 7 /)
                neighborDirCode_finerLevel(1:4) = (/ 42, 40, 41, 39 /)

            endif

        case(5)  ! '__5/___'  or  '__N'
            if (params%dim == 2) then
                neighborDirCode_sameLevel = 1
                lvl_down_neighbor = .true.
                thereMustBeANeighbor = .true.
                ! virtual treecodes, list_ids for neighbors on higher level
                tcFinerAppendDigit(1:2)         = (/ 0, 1 /)
                neighborDirCode_finerLevel(1:2) = (/ 10,  9 /)
                ! neighbor code for coarser neighbors
                if ( tc_last == 0) then
                    neighborDirCode_coarserLevel = 10
                elseif ( tc_last == 1) then
                    neighborDirCode_coarserLevel = 9
                end if
            
            else
                neighborDirCode_sameLevel    = 5
                thereMustBeANeighbor = .true.

                ! If the neighbor is coarser, then we have only one possible block, but
                ! the finer block (me) may be at four positions, which define the neighborhood code
                if ( tc_last == 0) then
                    neighborDirCode_coarserLevel = 46
                elseif ( tc_last == 1) then
                    neighborDirCode_coarserLevel = 44
                elseif ( tc_last == 4) then
                    neighborDirCode_coarserLevel = 45
                elseif ( tc_last == 5) then
                    neighborDirCode_coarserLevel = 43
                end if
                lvl_down_neighbor = .true.

                ! virtual treecodes, list_ids for neighbors on higher level
                tcFinerAppendDigit(1:4)  = (/ 0, 1, 4, 5 /)
                neighborDirCode_finerLevel(1:4) = (/ 46, 44, 45, 43 /)

            endif

        case(6)  ! '__6/___'
            neighborDirCode_sameLevel    = 6
            thereMustBeANeighbor = .true.

            ! If the neighbor is coarser, then we have only one possible block, but
            ! the finer block (me) may be at four positions, which define the neighborhood code
            if ( tc_last == 0) then
                neighborDirCode_coarserLevel = 50
            elseif ( tc_last == 1) then
                neighborDirCode_coarserLevel = 49
            elseif ( tc_last == 2) then
                neighborDirCode_coarserLevel = 47
            elseif ( tc_last == 3) then
                neighborDirCode_coarserLevel = 48
            end if
            lvl_down_neighbor = .true.

            ! virtual treecodes, list_ids for neighbors on higher level
            tcFinerAppendDigit(1:4)  = (/ 0, 1, 2, 3 /)
            neighborDirCode_finerLevel(1:4) = (/ 50, 49, 47, 48 /)

        case(12)  ! '_12/___'
            neighborDirCode_sameLevel    = 7

            ! neighbor code for coarser neighbors
            if ( tc_last == 4) then
                neighborDirCode_coarserLevel = 52
            elseif ( tc_last == 6) then
                neighborDirCode_coarserLevel = 51
            end if
            lvl_down_neighbor = ( (tc_last == 4) .or. (tc_last == 6) )

            tcFinerAppendDigit(1:2)         = (/ 4, 6 /)
            neighborDirCode_finerLevel(1:2) = (/ 52, 51 /)

        case(13)  ! '_13/___'
            neighborDirCode_sameLevel    = 8

            ! neighbor code for coarser neighbors
            if ( tc_last == 6) then
                neighborDirCode_coarserLevel = 53
            elseif ( tc_last == 7) then
                neighborDirCode_coarserLevel = 54
            end if
            lvl_down_neighbor = ( (tc_last == 6) .or. (tc_last == 7) )

            tcFinerAppendDigit(1:2)         = (/ 6, 7 /)
            neighborDirCode_finerLevel(1:2) = (/ 53, 54 /)

        case(14)  ! '_14/___'
            neighborDirCode_sameLevel    = 9

            ! neighbor code for coarser neighbors
            if ( tc_last == 5) then
                neighborDirCode_coarserLevel = 56
            elseif ( tc_last == 7) then
                neighborDirCode_coarserLevel = 55
            end if
            lvl_down_neighbor = ( (tc_last == 5) .or. (tc_last == 7) )

            tcFinerAppendDigit(1:2)         = (/ 5, 7 /)
            neighborDirCode_finerLevel(1:2) = (/ 56, 55 /)

        case(15)  ! '_15/___'
            neighborDirCode_sameLevel    = 10

            ! neighbor code for coarser neighbors
            if ( tc_last == 4) then
                neighborDirCode_coarserLevel = 58
            elseif ( tc_last == 5) then
                neighborDirCode_coarserLevel = 57
            end if
            lvl_down_neighbor = ( (tc_last == 4) .or. (tc_last == 5) )

            tcFinerAppendDigit(1:2)         = (/ 4, 5 /)
            neighborDirCode_finerLevel(1:2) = (/ 58, 57 /)

        case(62)  ! '_62/___'
            neighborDirCode_sameLevel    = 11

            ! neighbor code for coarser neighbors
            if ( tc_last == 0) then
                neighborDirCode_coarserLevel = 60
            elseif ( tc_last == 2) then
                neighborDirCode_coarserLevel = 59
            end if
            lvl_down_neighbor = ( (tc_last == 0) .or. (tc_last == 2) )

            tcFinerAppendDigit(1:2)         = (/ 0, 2 /)
            neighborDirCode_finerLevel(1:2) = (/ 60, 59 /)

        case(63)  ! '_63/___'
            neighborDirCode_sameLevel    = 12

            ! neighbor code for coarser neighbors
            if ( tc_last == 2) then
                neighborDirCode_coarserLevel = 61
            elseif ( tc_last == 3) then
                neighborDirCode_coarserLevel = 62
            end if
            lvl_down_neighbor = ( (tc_last == 2) .or. (tc_last == 3) )

            tcFinerAppendDigit(1:2)         = (/ 2, 3 /)
            neighborDirCode_finerLevel(1:2) = (/ 61, 62 /)

        case(64)  ! '_64/___'
            neighborDirCode_sameLevel    = 13

            ! neighbor code for coarser neighbors
            if ( tc_last == 1) then
                neighborDirCode_coarserLevel = 64
            elseif ( tc_last == 3) then
                neighborDirCode_coarserLevel = 63
            end if
            lvl_down_neighbor = ( (tc_last == 1) .or. (tc_last == 3) )

            tcFinerAppendDigit(1:2)         = (/ 1, 3 /)
            neighborDirCode_finerLevel(1:2) = (/ 64, 63 /)

        case(65)  ! '_65/___'
            neighborDirCode_sameLevel    = 14

            ! neighbor code for coarser neighbors
            if ( tc_last == 0) then
                neighborDirCode_coarserLevel = 66
            elseif ( tc_last == 1) then
                neighborDirCode_coarserLevel = 65
            end if
            lvl_down_neighbor = ( (tc_last == 0) .or. (tc_last == 1) )

            tcFinerAppendDigit(1:2)         = (/ 0, 1 /)
            neighborDirCode_finerLevel(1:2) = (/ 66, 65 /)

        case(23)  ! '_23/___'  or  '_SW'
            if (params%dim == 2) then
                neighborDirCode_sameLevel     = 8
                neighborDirCode_coarserLevel  = 8
                neighborDirCode_finerLevel(1) = 8
                tcFinerAppendDigit(1) = 2
                ! only sister block 1, 2 can have valid NE neighbor at one level down
                if ( (tc_last == 1) .or. (tc_last == 2) ) then
                    lvl_down_neighbor = .true.
                end if
            else
                neighborDirCode_sameLevel    = 15

                ! neighbor code for coarser neighbors
                if ( tc_last == 2) then
                    neighborDirCode_coarserLevel = 68
                elseif ( tc_last == 6) then
                    neighborDirCode_coarserLevel = 67
                end if
                lvl_down_neighbor = ( (tc_last == 2) .or. (tc_last == 6) )

                tcFinerAppendDigit(1:2)         = (/ 2, 6 /)
                neighborDirCode_finerLevel(1:2) = (/ 68, 67 /)
            endif

        case(25)  ! '_25/___'  or  '_NW'
            if (params%dim == 2) then
                neighborDirCode_sameLevel     = 6
                neighborDirCode_coarserLevel  = 6
                neighborDirCode_finerLevel(1) = 6
                tcFinerAppendDigit(1) = 0
                ! only sister block 0, 3 can have valid NW neighbor at one level down
                if ( (tc_last == 0) .or. (tc_last == 3) ) then
                    lvl_down_neighbor = .true.
                end if
            else
                neighborDirCode_sameLevel    = 16

                ! neighbor code for coarser neighbors
                if ( tc_last == 0) then
                    neighborDirCode_coarserLevel = 70
                elseif ( tc_last == 4) then
                    neighborDirCode_coarserLevel = 69
                end if
                lvl_down_neighbor = ( (tc_last == 0) .or. (tc_last == 4) )

                tcFinerAppendDigit(1:2)         = (/ 0, 4 /)
                neighborDirCode_finerLevel(1:2) = (/ 70, 69 /)
            endif

        case(43)  ! '_43/___'  or  '_SE'
            if (params%dim == 2) then
                neighborDirCode_sameLevel     = 7
                neighborDirCode_coarserLevel  = 7
                neighborDirCode_finerLevel(1) = 7
                tcFinerAppendDigit(1) = 3
                ! only sister block 0, 3 can have valid SE neighbor at one level down
                if ( (tc_last == 0) .or. (tc_last == 3) ) then
                    lvl_down_neighbor = .true.
                end if
            else
                neighborDirCode_sameLevel    = 17

                ! neighbor code for coarser neighbors
                if ( tc_last == 3) then
                    neighborDirCode_coarserLevel = 72
                elseif ( tc_last == 7) then
                    neighborDirCode_coarserLevel = 71
                end if
                lvl_down_neighbor = ( (tc_last == 3) .or. (tc_last == 7) )

                tcFinerAppendDigit(1:2)         = (/ 3, 7 /)
                neighborDirCode_finerLevel(1:2) = (/ 72, 71 /)
            endif

        case(45)  ! '_45/___'  or  '_NE'
            if (params%dim == 2) then
                neighborDirCode_sameLevel     = 5
                neighborDirCode_coarserLevel  = 5
                neighborDirCode_finerLevel(1) = 5
                tcFinerAppendDigit(1) = 1
                ! only sister block 1, 2 can have valid NE neighbor at one level down
                if ( (tc_last == 1) .or. (tc_last == 2) ) then
                    lvl_down_neighbor = .true.
                end if
            else
                neighborDirCode_sameLevel    = 18

                ! neighbor code for coarser neighbors
                if ( tc_last == 1) then
                    neighborDirCode_coarserLevel = 74
                elseif ( tc_last == 5) then
                    neighborDirCode_coarserLevel = 73
                end if
                lvl_down_neighbor = ( (tc_last == 1) .or. (tc_last == 5) )

                tcFinerAppendDigit(1:2)         = (/ 1, 5 /)
                neighborDirCode_finerLevel(1:2) = (/ 74, 73 /)
            endif

        case(123)  ! '123/___'
            neighborDirCode_sameLevel      = 19
            neighborDirCode_coarserLevel   = 19
            neighborDirCode_finerLevel(1)  = 19
            tcFinerAppendDigit(1) = 6
            lvl_down_neighbor = ( tc_last == 6 )

        case(134)  ! '134/___'
            neighborDirCode_sameLevel      = 20
            neighborDirCode_coarserLevel   = 20
            neighborDirCode_finerLevel(1)  = 20
            tcFinerAppendDigit(1) = 7
            lvl_down_neighbor = ( tc_last == 7 )

        case(145)  ! '145/___'
            neighborDirCode_sameLevel      = 21
            neighborDirCode_coarserLevel   = 21
            neighborDirCode_finerLevel(1)  = 21
            tcFinerAppendDigit(1) = 5
            lvl_down_neighbor = ( tc_last == 5 )

        case(152)  ! '152/___'
            neighborDirCode_sameLevel      = 22
            neighborDirCode_coarserLevel   = 22
            neighborDirCode_finerLevel(1)  = 22
            tcFinerAppendDigit(1) = 4
            lvl_down_neighbor = ( tc_last == 4 )

        case(623)  ! '623/___'
            neighborDirCode_sameLevel      = 23
            neighborDirCode_coarserLevel   = 23
            neighborDirCode_finerLevel(1)  = 23
            tcFinerAppendDigit(1) = 2
            lvl_down_neighbor = ( tc_last == 2 )

        case(634)  ! '634/___'
            neighborDirCode_sameLevel      = 24
            neighborDirCode_coarserLevel   = 24
            neighborDirCode_finerLevel(1)  = 24
            tcFinerAppendDigit(1) = 3
            lvl_down_neighbor = ( tc_last == 3 )

        case(645)  ! '645/___'
            neighborDirCode_sameLevel      = 25
            neighborDirCode_coarserLevel   = 25
            neighborDirCode_finerLevel(1)  = 25
            tcFinerAppendDigit(1) = 1
            lvl_down_neighbor = ( tc_last == 1 )

        case(652)  ! '652/___'
            neighborDirCode_sameLevel      = 26
            neighborDirCode_coarserLevel   = 26
            neighborDirCode_finerLevel(1)  = 26
            tcFinerAppendDigit(1) = 0
            lvl_down_neighbor = ( tc_last == 0 )
        case default
            call abort(636300, "A weird error occured.")

    end select

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 1) Check if we find a neighbor on the SAME LEVEL
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! calculate treecode for neighbor on same level
    ! call adjacent(tcBlock, tcNeighbor, dir, level, Jmax, params%dim)
    call adjacent_wrapper_b(tcb_Block, tcb_Neighbor, dir, level=level, dim=params%dim, max_level=Jmax)

    ! check if (hypothetical) neighbor exists and if so find its lgtID
    call doesBlockExist_tree(tcb_Neighbor, exists, lgtID_neighbor, dim=params%dim, level=level, tree_id=tree_ID, max_level=params%Jmax)

    if (exists) then
        ! we found the neighbor on the same level.
        hvy_neighbor( hvyID_block, neighborDirCode_sameLevel ) = lgtID_neighbor
        return
    endif

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 2) Check if we find a neighbor on the COARSER LEVEL (if that is possible at all)
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (lvl_down_neighbor) then
        ! We did not find the neighbor on the same level, and now check on the coarser level.
        ! Just clear all levels before the lower one
        tcb_Neighbor = tc_clear_until_level_b(tcb_Neighbor, dim=params%dim, level=level-1, max_level=Jmax)
        ! check if (hypothetical) neighbor exists and if so find its lgtID
        call doesBlockExist_tree(tcb_Neighbor, exists, lgtID_neighbor, dim=params%dim, level=level-1, tree_id=tree_ID, max_level=params%Jmax)

        ! call tcb2array(tcb_Neighbor, tc_test, dim=params%dim, level=level, max_level=Jmax)
        ! ! write(*, '("neighbour ", a, " orig=", b64.64, " neig=", b64.64)') dir, tcb_Block, tcb_Neighbor
        ! write(*, '("Coarse neighbour dir=", a, " ta=", 13(i0), " tcb=", 13(i0))') dir, tcNeighbor, tc_test

        if ( exists ) then
            ! neighbor is one level down (coarser)
            hvy_neighbor( hvyID_block, neighborDirCode_coarserLevel ) = lgtID_neighbor
            return
        endif
    endif


    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! 3) Check if we find a neighbor on the FINER LEVEL
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! Note: there are several neighbors possible on finer levels (up to four in 3D!)
    ! loop over all 4 possible neighbors
    if (level < Jmax) then
        do k = 1, 4
            if (tcFinerAppendDigit(k) /= -1) then
                ! first neighbor virtual treecode, one level up
                tcb_Virtual = tc_set_digit_at_level_b(tcb_Block, tcFinerAppendDigit(k), level=level+1, max_level=Jmax, dim=params%dim)

                ! calculate treecode for neighbor on same level (virtual level)
                call adjacent_wrapper_b(tcb_Virtual, tcb_Neighbor, dir, level=level+1, max_level=Jmax, dim=params%dim)
                ! check if (hypothetical) neighbor exists and if so find its lgtID
                call doesBlockExist_tree(tcb_Neighbor, exists, lgtID_neighbor, dim=params%dim, level=level+1, tree_id=tree_ID, max_level=params%Jmax)

                if (exists) then
                    hvy_neighbor( hvyID_block, neighborDirCode_finerLevel(k) ) = lgtID_neighbor
                end if

                ! we did not find a neighbor. that may be a bad grid error, or simply, there is none
                ! because symmetry conditions are used.
                if (thereMustBeANeighbor) then
                    if ((.not. exists .and. ALL(params%periodic_BC)).or.(maxval(abs(n_domain))==0.and..not.exists)) then
                        ! construct print format dynamically after Jmax
                        write(*, '("Dir ", i0, ", lvl_down=", l1, " lvl=", i0, ":", 4(1x, i0), " TC: ", b64.64)') &
                            dir, lvl_down_neighbor, level+1, tcFinerAppendDigit, tcb_Block
                        error = .true.
                    endif
                endif
            endif
        end do
    endif
end subroutine