!> \brief reset ghosts nodes for all blocks (not just active ones) for debuging
!> ghost nodes are set to a very large constant.
!! input:    - params, light and heavy data \n
!! output:   - heavy data array
! ********************************************************************************************

subroutine reset_ghost_nodes(  params, hvy_block, tree_ID, s_M2M, s_M2C, s_M2F)
    implicit none

    type (type_params), intent(in)      :: params                     !< user defined parameter structure
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)   !< heavy data array - block data
    integer(kind=ik), intent(in)        :: tree_ID                    !< which tree to study
    logical, intent(in), optional  :: s_M2M                         !< Synch from level J   to J
    logical, intent(in), optional  :: s_M2C                         !< Synch from level J   to J-1
    logical, intent(in), optional  :: s_M2F                         !< Synch from level J   to J+1

    integer(kind=ik)               :: g, hvy_id, k_b, k_n, lgt_id, lgt_id_n, level_me, level_diff
    integer(kind=ik)               :: Bs(3), ijk(2,3)
    logical :: SM2M, SM2C, SM2F

    sM2M = .false.
    sM2C = .false.
    sM2F = .false.
    if (present(s_M2M)) sM2M = s_M2M
    if (present(s_M2C)) sM2C = s_M2C
    if (present(s_M2F)) sM2F = s_M2F

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
            level_diff = level_me - lgt_block( lgt_id_n, IDX_MESH_LVL )

            if ((level_diff==0 .and. sM2M) .or. (level_diff==-1 .and. sM2F) .or. (level_diff==+1 .and. sM2C)) then
              ! get indices of ghost patch that is affected
              ijk = ijkPatches(:,:, k_n, level_diff, RECVER)
              ! call get_indices_of_modify_patch(params, k_n, ijk, (/ size(hvy_block, 1), size(hvy_block, 2), size(hvy_block, 3)/), &
              !   (/params%g, params%g, params%g/), (/params%g, params%g, params%g/))
              hvy_block(ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), :, hvy_ID )           = 9.0e9_rk
            endif
          endif
        enddo
        ! !-- x-direction
        ! hvy_block(1:g, :, :, :, hvy_ID )           = 9.0e9_rk
        ! hvy_block(Bs(1)+g+1:Bs(1)+2*g, :, :, :, hvy_ID) = 9.0e9_rk
        ! !-- y-direction
        ! hvy_block(:, 1:g, :, :, hvy_ID)           = 9.0e9_rk
        ! hvy_block(:, Bs(2)+g+1:Bs(2)+2*g, :, :, hvy_ID) = 9.0e9_rk
        ! !-- z-direction
        ! if ( params%dim == 3 ) then
        !   hvy_block(:, :, 1:g, :, hvy_ID)           = 9.0e9_rk
        !   hvy_block(:, :, Bs(3)+g+1:Bs(3)+2*g, :, hvy_ID) = 9.0e9_rk
        ! end if
    enddo
end subroutine reset_ghost_nodes
