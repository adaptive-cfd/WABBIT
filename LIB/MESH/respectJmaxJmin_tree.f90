subroutine respectJmaxJmin_tree( params, tree_ID, Jmin_set, Jmax_set, stay_value )
    implicit none

    type (type_params), intent(in)      :: params                       !< User defined parameter structure
    integer(kind=ik), intent(in)        :: tree_ID                      !< Tree we look at
    integer(kind=ik), intent(in), optional :: Jmin_set, Jmax_set        !< Optional Jmin/Jmax overrides
    integer(kind=ik), intent(in), optional :: stay_value                !< Value to set when refinement/coarsening is not allowed, defaults to 0
    integer(kind=ik)                    :: Jmax, Jmin                   !< Treelevel restrictions
    integer(kind=ik)                    :: k, lgt_id, ref_stat, level   !< Loop variables
    integer(kind=ik)                    :: stay_val

    Jmax = params%Jmax
    Jmin = params%Jmin
    stay_val = 0
    if (present(Jmin_set)) Jmin = Jmin_set
    if (present(Jmax_set)) Jmax = Jmax_set
    if (present(stay_value)) stay_val = stay_value

    ! loop over all active blocks
    do k = 1, lgt_n(tree_ID)

        lgt_id = lgt_active(k, tree_ID)

        ref_stat = lgt_block( lgt_id, IDX_REFINE_STS )
        level = lgt_block( lgt_id, IDX_MESH_LVL )

        if ((ref_stat == +1).and.(level >= Jmax)) then
            ! can not refine (set flag to 0 = stay)
            lgt_block( lgt_id, IDX_REFINE_STS ) = stay_val
        end if

        if ((ref_stat == -1).and.(level <= Jmin)) then
            ! can not coarsen (set flag to 0 = stay)
            lgt_block( lgt_id, IDX_REFINE_STS ) = stay_val
        end if

    end do

end subroutine respectJmaxJmin_tree
