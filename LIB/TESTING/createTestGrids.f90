! in here are some settings to create some grids that are useful for testing

!---------------------------------------------------------------------------
! One adapted block in +x,-y,-z
!    2 2 3 3
!    2 2 3 3
!    1 1 4 5
!    1 1 6 7
!---------------------------------------------------------------------------
subroutine createGrid_simple_adaptive( params, hvy_block, tree_ID )
    implicit none

    type (type_params), intent(inout)       :: params                     !> user defined parameter structure
    real(kind=rk),  intent(inout)           :: hvy_block(:, :, :, :, :)   !> heavy data array - block data
    integer(kind=ik), intent(in)            :: tree_ID                    !> treeID of grid

    integer(kind=ik)           :: k, lgt_id, hvy_id
    integer(kind=tsize)        :: treecode

    ! init grid
    call createEquidistantGrid_tree( params, hvy_block, 1, .true., tree_ID )

    ! choose block to refine
    do k = 1, hvy_n(tree_ID)
        hvy_ID = hvy_active(k, tree_ID)
        call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)
        treecode = get_tc(lgt_block(lgt_id, IDX_TC_1 : IDX_TC_2))
        ! if (tc_get_digit_at_level_b( treecode, params%dim, 1, params%Jmax) == 3-tc_get_digit_at_level_b( treecode, params%dim, 2, params%Jmax)) then
        if (tc_get_digit_at_level_b( treecode, params%dim, 1, params%Jmax) == 2) then
                lgt_block(lgt_id, IDX_REFINE_STS) = +1
        endif
    enddo

    ! refine and update so that all blocks know of each other
    if ( params%dim == 3 ) then
        call refinementExecute3D_tree( params, hvy_block, tree_ID )
    else
        call refinementExecute2D_tree( params, hvy_block(:,:,1,:,:), tree_ID )
    endif
    call updateMetadata_tree(params, tree_ID)

end subroutine



!---------------------------------------------------------------------------
! Like 1, but with a layer of blocks on lvl 2 around it, adapted block is on lvl 3
! so that all relations can be identified without crossing periodicity border
!    g g f f e e d d
!    g g f f e e d d
!    h h 2 2 3 3 c c
!    h h 2 2 3 3 c c
!    i i 1 1 4 5 b b
!    i i 1 1 6 7 b b
!    8 8 9 9 0 0 a a
!    8 8 9 9 0 0 a a
!---------------------------------------------------------------------------
subroutine createGrid_simple_adaptive_with_borders( params, hvy_block, tree_ID )
    implicit none

    type (type_params), intent(inout)       :: params                     !> user defined parameter structure
    real(kind=rk),  intent(inout)           :: hvy_block(:, :, :, :, :)   !> heavy data array - block data
    integer(kind=ik), intent(in)            :: tree_ID                    !> treeID of grid

    integer(kind=ik)           :: k, lgt_id, hvy_id
    integer(kind=tsize)        :: treecode

    ! init grid
    call createEquidistantGrid_tree( params, hvy_block, 2, .true., tree_ID )

    ! choose block to refine - one block, has TC 1.2
    do k = 1, hvy_n(tree_ID)
        hvy_ID = hvy_active(k, tree_ID)
        call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)
        treecode = get_tc(lgt_block(lgt_id, IDX_TC_1 : IDX_TC_2))
        if (tc_get_digit_at_level_b( treecode, params%dim, 1, params%Jmax) == 1 .and. &
            ((tc_get_digit_at_level_b( treecode, params%dim, 2, params%Jmax) == 2 .and. params%dim == 2) .or. &
             (tc_get_digit_at_level_b( treecode, params%dim, 2, params%Jmax) == 6 .and. params%dim == 3))) then
                lgt_block(lgt_id, IDX_REFINE_STS) = +1
        endif
    enddo

    ! refine and update so that all blocks know of each other
    if ( params%dim == 3 ) then
        call refinementExecute3D_tree( params, hvy_block, tree_ID )
    else
        call refinementExecute2D_tree( params, hvy_block(:,:,1,:,:), tree_ID )
    endif
    call updateMetadata_tree(params, tree_ID)

end subroutine