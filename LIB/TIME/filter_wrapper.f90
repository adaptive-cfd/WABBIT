subroutine filter_wrapper(time, params, hvy_block, hvy_tmp, hvy_mask, tree_ID)
   implicit none

    real(kind=rk), intent(in)           :: time
    type (type_params), intent(inout)   :: params                       !> user defined parameter structure, hvy_active
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)     !> heavy data array - block data
    !> heavy temp data: used for saving, filtering, and helper qtys (reaction rate, mask function)
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    !> hvy_mask are qty that depend on the grid and not explicitly on time
    real(kind=rk), intent(inout)        :: hvy_mask(:, :, :, :, :)
    integer(kind=ik), intent(in)        :: tree_ID
    real(kind=rk), dimension(3)         :: dx, x0                       !> spacing and origin of a block
    integer(kind=ik)                    :: k, dF, neqn, lgt_id, i       ! loop variables
    integer(kind=ik)                    :: g                            ! grid parameter, error variable
    integer(kind=ik), dimension(3)      :: Bs
    integer(kind=2)                     :: n_domain(1:3)                ! surface normal
    integer(kind=ik)                    :: level, hvy_id

    Bs = params%Bs
    g  = params%g
    n_domain = 0

    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)

        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

        level = lgt_block(lgt_id, params%Jmax+IDX_MESH_LVL)

        if ((params%filter_only_maxlevel .and. level==params%Jmax) .or. &
            (.not.params%filter_only_maxlevel .and. .not.params%filter_all_except_maxlevel) .or. &
            (params%filter_all_except_maxlevel .and. level/=params%Jmax) ) then

            ! get block spacing for RHS
            call get_block_spacing_origin( params, lgt_id, x0, dx )

            if ( .not. All(params%periodic_BC) ) then
                ! check if block is adjacent to a boundary of the domain, if this is the case we use one sided stencils
                call get_adjacent_boundary_surface_normal( lgt_block(lgt_id, 1:lgt_block(lgt_id,params%Jmax+IDX_MESH_LVL)), &
                params%domain_size, params%Bs, params%dim, n_domain )
            endif

            call filter_meta(params%physics_type, time, hvy_block(:,:,:,:, hvy_id), g, x0, dx, &
            hvy_tmp(:,:,:,:,hvy_id), hvy_mask(:,:,:,:,hvy_id), n_domain)
        endif
    enddo
end subroutine filter_wrapper
