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
! ********************************************************************************************
subroutine updateFamily_tree(params, tree_ID)

    implicit none
    type (type_params), intent(in)      :: params                   !> user defined parameter structure
    integer(kind=ik), intent(in)        :: tree_ID
    integer(kind=ik)                    :: k, N, Jmax, lgtID, hvyID
    integer(kind=ik)                    :: lgtID_fam(7)


    N = params%number_blocks
    Jmax = params%Jmax


    ! loop over active heavy data blocks
    do k = 1, hvy_n(tree_ID)
        ! delete existing neighbors (they are assumed to be outdated when calling this routine)
        hvyID = hvy_active(k, tree_ID)
        hvy_family(hvyID, :) = -1

        call hvy2lgt( lgtID, hvyID, params%rank, N )

        ! search for family
        call find_mother(params, lgtID, hvy_family(hvyID, 1))
        call find_sisters(params, lgtID, hvy_family(hvyID, 2:1+2**params%dim))
        call find_daughters(params, lgtID, hvy_family(hvyID, 2+2**params%dim:1+2**(params%dim+1)))

        ! error check - we should always find atleast one sister, which is the block itself
        if (all(hvy_family(hvyID, 2:1+2**params%dim) == -1)) then
            call abort(2405061, "Block does not find any sister, probably something went wrong")
        endif
        
    end do

end subroutine updateFamily_tree
