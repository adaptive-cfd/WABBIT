!> \brief functional wrapper for the 2D and 3D version of update family
!! input:     - light data array
!!            - params struct
!! output:    - neighbor list array
!!
!> family relations: \n
!  ---------------
!> Index 1:                         Mother, -1 if not existent
!! Index 2-5 (2D) or 2-9 (3D):     Sisters with last TC 0-7
!! Index 6-9 (2D) or 10-17 (3D): Daugthers with last TC 0-7, -1 if not existent
!
!  One could think that for non-CVS computations the mother and daughters are not used,
!  however, for refinement and coarsening they are comfortable because all operations become
!  a send between the individual blocks
! ********************************************************************************************
subroutine updateFamily_tree(params, tree_ID)

    implicit none
    type (type_params), intent(in)      :: params                   !> user defined parameter structure
    integer(kind=ik), intent(in)        :: tree_ID
    integer(kind=ik)                    :: k, lgtID, hvyID


    ! loop over active heavy data blocks
    do k = 1, hvy_n(tree_ID)
        ! delete existing neighbors (they are assumed to be outdated when calling this routine)
        hvyID = hvy_active(k, tree_ID)
        hvy_family(hvyID, :) = -1

        call hvy2lgt( lgtID, hvyID, params%rank, params%number_blocks )

        ! search for family
        call find_mother(params, lgtID, hvy_family(hvyID, 1))
        call find_sisters(params, lgtID, hvy_family(hvyID, 2:1+2**params%dim))
        call find_daughters(params, lgtID, hvy_family(hvyID, 2+2**params%dim:1+2**(params%dim+1)))

        ! write(*, '("3R", i1, " B", i2, " S", i5, " L", i2, " R", i2, " F ", 9(i5, 1x), " TC", 1(b32.32))') &
        ! params%rank, k, lgtID, lgt_block(lgtID, IDX_MESH_LVL), lgt_block(lgtID, IDX_REFINE_STS), hvy_family(hvy_ID, :), lgt_block(lgt_ID, IDX_TC_2)

        ! error check - we should always find atleast one sister, which is the block itself
        if (all(hvy_family(hvyID, 2:1+2**params%dim) == -1)) then
            call abort(2405061, "Block does not find any sister, probably something went wrong")
        endif
        
    end do

end subroutine updateFamily_tree
