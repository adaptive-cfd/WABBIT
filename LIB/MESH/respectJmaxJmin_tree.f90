subroutine respectJmaxJmin_tree( params, tree_ID)
    implicit none

    type (type_params), intent(in)      :: params                       !< User defined parameter structure
    integer(kind=ik), intent(in)        :: tree_ID                      !< Tree we look at
    integer(kind=ik)                    :: Jmax, Jmin                   !< Treelevel restrictions
    integer(kind=ik)                    :: k, lgt_id, ref_stat, level  !< Loop variables

    Jmax = params%Jmax
    Jmin = params%Jmin

    ! loop over all active blocks
    do k = 1, lgt_n(tree_ID)

        lgt_id = lgt_active(k, tree_ID)

        ref_stat = lgt_block( lgt_id, IDX_REFINE_STS )
        level = lgt_block( lgt_id, IDX_MESH_LVL )

        if ((ref_stat == +1).and.(level >= Jmax)) then
            ! can not refine (set flag to 0 = stay)
            lgt_block( lgt_id, IDX_REFINE_STS ) = 0
        end if

        if ((ref_stat == -1).and.(level <= Jmin)) then
            ! can not coarsen (set flag to 0 = stay)
            lgt_block( lgt_id, IDX_REFINE_STS ) = 0
        end if

    end do

end subroutine respectJmaxJmin_tree
