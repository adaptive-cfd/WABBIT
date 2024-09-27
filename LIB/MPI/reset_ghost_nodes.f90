!> \brief Reset ghosts nodes for all blocks (not just active ones) for debugging
!> Ghost nodes are set to a very large constant.
!> Currently only works if we are synching/overwriting all available ghost-patches,
!! so not for level-wise or ref-wise or any other fancy stuff
!! input:    - params, light and heavy data \n
!! output:   - heavy data array
! ********************************************************************************************
subroutine reset_ghost_nodes(  params, hvy_block, tree_ID)
    implicit none

    type (type_params), intent(in)      :: params                     !< user defined parameter structure
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)   !< heavy data array - block data
    integer(kind=ik), intent(in)        :: tree_ID                    !< which tree to study

    integer(kind=ik)               :: g, hvy_id, k_b, k_n, lgt_id, lgt_id_n, level_me, lvl_diff
    integer(kind=ik)               :: Bs(3), ijk(2,3)

    Bs = params%Bs
    g  = params%g

    do k_b = 1, hvy_n(tree_ID)
        hvy_ID = hvy_active(k_b, tree_ID)
        call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
        level_me = lgt_block( lgt_id, IDX_MESH_LVL )

        ! loop over all neighborhoods, check if this is synched and then wipe it
        do k_n = 1, size(hvy_neighbor, 2)
          lgt_id_n = hvy_neighbor( hvy_id, k_n )
          if ( lgt_id_n /= -1 ) then
            lvl_diff = level_me - lgt_block( lgt_id_n, IDX_MESH_LVL )

            ! get indices of ghost patch that is affected
            ijk = ijkPatches(:,:, k_n, RECVER)
            hvy_block(ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), :, hvy_ID )           = 9.0e9_rk
          endif
        enddo
    enddo
end subroutine reset_ghost_nodes
