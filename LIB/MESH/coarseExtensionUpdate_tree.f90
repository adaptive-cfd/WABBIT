subroutine coarseExtensionUpdate_tree( params, hvy_block, hvy_tmp, tree_ID )
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params

    implicit none

    type (type_params), intent(in)      :: params
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik), intent(in)        :: tree_ID

    integer(kind=ik)                    :: Jmax_active, Jmin_active, level

    Jmin_active = minActiveLevel_tree(tree_ID)
    Jmax_active = maxActiveLevel_tree(tree_ID)

    do level = Jmax_active, Jmin_active, -1
        call coarseExtensionUpdate_level( params, lgt_block, hvy_block, hvy_tmp, hvy_neighbor, hvy_active(:,tree_ID), &
        hvy_n(tree_ID),lgt_n(tree_ID), inputDataSynced=.false., level=level, hvy_details=hvy_details )
    enddo

end subroutine