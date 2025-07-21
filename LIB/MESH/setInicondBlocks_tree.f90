! On an existing tree, set the initial condition (via the physics modules) on all blocks
! This routine does not modify the tree; it is simply a wrapper to set the data on the blocks
subroutine setInicondBlocks_tree(params, hvy_block, tree_ID)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params

    implicit none

    type (type_params), intent(inout)    :: params                              !> user defined parameter structure
    real(kind=rk), intent(inout)         :: hvy_block(:, :, :, :, :)            !> heavy data array - block data
    integer(kind=ik), intent(in)         :: tree_ID

    integer(kind=ik)                     :: k, g, Bs(1:3), hvy_id, lgt_id       ! loop variable
    real(kind=rk)                        :: x0(1:3), dx(1:3)                    ! origin and spacing of blocks
    integer(kind=2) :: n_domain(1:3)

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas


    Bs = params%Bs
    g  = params%g
    ! The normal on the domain indicates (if non-periodic BC are used), if a block
    ! is at the outer, usually periodic border of the domain ( x,y,z == 0 and x,y,z == L)
    ! Nonzero values indicate this is the case, e.g., n_domain=(/1, 0, -1/) means in x-axis, our block
    ! is way at the back and its boundary normal points in +x, and in z, its at the bottom (z=0), thus
    ! its normal points downwards.
    n_domain = 0

    !---------------------------------------------------------------------------
    ! on the existig tree, evaluate the initial condition
    !---------------------------------------------------------------------------
    ! loop over my active heavy data
    do k = 1, hvy_n(tree_ID)
        ! hvy_id of the block we're looking at
        hvy_id = hvy_active(k, tree_ID)

        ! light id of this block
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

        ! compute block spacing and origin from treecode
        call get_block_spacing_origin( params, lgt_id, x0, dx )

        if ( .not. All(params%periodic_BC) ) then
            ! check if block is adjacent to a boundary of the domain, if this is the case we use one sided stencils
            call get_adjacent_boundary_surface_normal( params, lgt_id, n_domain )
        endif

        ! set the initial condition on this block
        call INICOND_meta(params%physics_type, 0.0_rk, hvy_block(:,:,:,:,hvy_id), g, &
        x0, dx, n_domain)
    enddo

end subroutine setInicondBlocks_tree
