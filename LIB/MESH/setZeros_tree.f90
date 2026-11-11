subroutine setZeros_tree(params, hvy_block, tree_ID)    
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params    
    
    type (type_params), intent(inout)    :: params                      !> user defined parameter structure
    real(kind=rk), intent(inout)         :: hvy_block(:, :, :, :, :)    !> heavy data array - block data
    integer(kind=ik), intent(in)         :: tree_ID

    integer :: k, hvy_id

    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)

        hvy_block(:,:,:,:,hvy_id) = 0.0_rk
    enddo


end subroutine



subroutine copyMask_tree(params, hvy_block, hvy_mask, tree_ID)    
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params    
    
    type (type_params), intent(inout)    :: params                      !> user defined parameter structure
    real(kind=rk), intent(inout)         :: hvy_block(:, :, :, :, :)    !> heavy data array - block data
    real(kind=rk), intent(inout)         :: hvy_mask(:, :, :, :, :)    !> heavy data array - block data
    integer(kind=ik), intent(in)         :: tree_ID

    integer :: k, hvy_id

    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)

        hvy_block(:,:,:,1,hvy_id) = hvy_mask(:,:,:,1,hvy_id)
    enddo


end subroutine