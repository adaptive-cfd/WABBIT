subroutine set_inicond_blocks(params, lgt_block, hvy_block, hvy_active, hvy_n)

  implicit none

    type (type_params), intent(inout)    :: params                              !> user defined parameter structure
    integer(kind=ik), intent(inout)      :: lgt_block(:, :)                     !> light data array
    real(kind=rk), intent(inout)         :: hvy_block(:, :, :, :, :)            !> heavy data array - block data
    integer(kind=ik), intent(inout)      :: hvy_active(:,:)                     !> list of active blocks (light data)
    integer(kind=ik), intent(inout)      :: hvy_n(:)                            !> number of heavy and light active blocks
    integer(kind=ik)                     :: k, g, Bs(1:3), hvy_id, lgt_id       ! loop variable
    real(kind=rk)                        :: x0(1:3), dx(1:3)                    ! origin and spacing of blocks
    integer(kind=2) :: n_domain(1:3)

    Bs = params%Bs
    g  = params%n_ghosts
    ! The normal on the domain indicates (if non-periodic BC are used), if a block
    ! is at the outer, usually periodic border of the domain ( x,y,z == 0 and x,y,z == L)
    ! Nonzero values indicate this is the case, e.g., n_domain=(/1, 0, -1/) means in x-axis, our block
    ! is way at the back and its boundary normal points in +x, and in z, its at the bottom (z=0), thus
    ! its normal points downwards.
    n_domain = 0


    !---------------------------------------------------------------------------
    ! on the grid, evaluate the initial condition
    !---------------------------------------------------------------------------
    ! loop over my active heavy data
    do k = 1, hvy_n(tree_ID_flow)
        ! hvy_id of the block we're looking at
        hvy_id = hvy_active(k, tree_ID_flow)

        ! light id of this block
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

        ! compute block spacing and origin from treecode
        call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

        ! get block spacing for RHS
        call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

        if ( .not. All(params%periodic_BC) ) then
            ! check if block is adjacent to a boundary of the domain, if this is the case we use one sided stencils
            call get_adjacent_boundary_surface_normal( lgt_block(lgt_id, 1:lgt_block(lgt_id,params%max_treelevel+IDX_MESH_LVL)), &
            params%domain_size, params%Bs, params%dim, n_domain )
        endif

        ! set the initial condition on this block
        call INICOND_meta(params%physics_type, 0.0_rk, hvy_block(:,:,:,:,hvy_id), g, &
        x0, dx, n_domain)
    enddo

end subroutine set_inicond_blocks
