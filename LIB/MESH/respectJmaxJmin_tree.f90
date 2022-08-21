subroutine respectJmaxJmin_tree( params, lgt_block, lgt_active, lgt_n, tree_ID)
    implicit none

    type (type_params), intent(in)      :: params           !> user defined parameter structure
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)  !> light data array
    integer(kind=ik), intent(in)        :: lgt_active(:,:)    !> list of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_n(:)         !> number of active blocks (light data)
    integer(kind=ik), intent(in)        :: tree_ID
    integer(kind=ik)                    :: Jmax, Jmin       ! treelevel restrictions
    integer(kind=ik)                    :: k                ! loop variables

    Jmax = params%max_treelevel
    Jmin = params%min_treelevel

    ! loop over all active blocks
    do k = 1, lgt_n(tree_ID)

        if ((lgt_block( lgt_active(k, tree_ID), Jmax + idx_refine_sts ) == +1_ik).and.(lgt_block( lgt_active(k, tree_ID), Jmax + IDX_MESH_LVL ) >= Jmax)) then
            ! can not refine (set flag to 0 = stay)
            lgt_block( lgt_active(k, tree_ID), Jmax + idx_refine_sts ) = 0_ik
        end if

        if ((lgt_block( lgt_active(k, tree_ID), Jmax + idx_refine_sts ) == -1_ik).and.(lgt_block( lgt_active(k, tree_ID), Jmax + IDX_MESH_LVL ) <= Jmin)) then
            ! can not coarsen
            lgt_block( lgt_active(k, tree_ID), Jmax + idx_refine_sts ) = 0_ik
        end if

    end do

end subroutine respectJmaxJmin_tree
