!> \brief Wrapper, that applies boundary conditions on whole tree or only on parts of it
subroutine bound_cond_generic(params, hvy_block, tree_ID, sync_case, spaghetti_form, s_val, edges_only)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params
    
    implicit none

    type (type_params), intent(in) :: params
    real(kind=rk), intent(inout)   :: hvy_block(:, :, :, :, :)      !< heavy data array - block data
    integer(kind=ik), intent(in)   :: tree_ID                       !< which tree to study
    !> String representing which kind of syncing we want to do, varies between restrictions (tree, level or ref)
    character(len=*), intent(in)   :: sync_case
    logical, intent(in)            :: spaghetti_form                !< if true, hvy_block is in spaghetti form
    !> Additional value to be considered for syncing logic, can be level or refinement status to which should be synced, used if sync case includes ref or level
    integer(kind=ik), intent(in), optional  :: s_val
    logical, intent(in), optional  :: edges_only                    !< if true, only set the boundary condition on the edges of the domain, default is false

    integer(kind=ik) :: i_dim, k, hvy_id, lgt_id, sync_id
    real(kind=rk) :: x0(1:3), dx(1:3), tolerance
    integer(kind=2) :: n_domain(1:3)
    logical :: edgesOnly

    edgesOnly = .false.
    if (present(edges_only)) edgesOnly = edges_only

    ! now lets treat the special restrictions, set to the second digit
    if (index(sync_case, "tree") > 0) sync_id = 0
    if (index(sync_case, "level") > 0) sync_id = 1
    if (index(sync_case, "ref") > 0) sync_id = 2
    if (index(sync_case, "leaf") > 0) sync_id = 3
    if (index(sync_case, "root") > 0) sync_id = 4

    if ( .not. All(params%periodic_BC) ) then
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k, tree_ID)
            call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

            if (sync_id == 1 .and. lgt_block(lgt_id, IDX_MESH_LVL) /= s_val) cycle
            if (sync_id == 2 .and. lgt_block(lgt_id, IDX_REFINE_STS) /= s_val) cycle
            if (sync_id == 3 .and. .not. block_is_leaf(params, hvy_id)) cycle
            if (sync_id == 4 .and. .not. block_is_root(params, hvy_id)) cycle

            ! compute block spacing and origin from treecode
            ! call get_block_spacing_origin( params, lgt_id, x0, dx )
            call get_block_spacing_origin_b( get_tc(lgt_block(lgt_id, IDX_TC_1 : IDX_TC_2)), params%domain_size, &
            params%Bs, x0, dx, dim=params%dim, level=lgt_block(lgt_id, IDX_MESH_LVL), max_level=params%Jmax)

            n_domain(:) = 0
            ! ! check if block is adjacent to a boundary of the domain, if this is the case we adapt the ghost patches
            ! call get_adjacent_boundary_surface_normal( params, lgt_id, n_domain )  ! this function relies in MESH module, so we cannot use it here

            tolerance = 1.0e-3_rk * minval(dx(1:params%dim))

            do i_dim = 1, params%dim
                ! check if origin_b = 0
                if (abs(x0(i_dim)-0.0_rk) < tolerance ) then !x_i == 0
                    n_domain(i_dim) = -1
                ! check if origin_b + BS*dx = L
                elseif (abs(x0(i_dim)+dx(i_dim)*real(params%Bs(i_dim),kind=rk) - params%domain_size(i_dim)) < tolerance) then ! x_i == L
                    n_domain(i_dim) = +1
                endif
            end do

            ! set the boundary condition on this block
            ! JB ToDo: Add time support in case we have time-dependent BCs
            call BOUNDCOND_meta(params%physics_type, 0.0_rk, hvy_block(:,:,:,:,hvy_id), params%g, x0, dx, n_domain, spaghetti_form=spaghetti_form, edges_only=edgesOnly)
        enddo
    endif

end subroutine