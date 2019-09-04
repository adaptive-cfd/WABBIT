!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name init_data.f90
!> \version 0.5
!> \author msr
!
!
!
! ********************************************************************************************
subroutine set_inicond_blocks(params, lgt_block, hvy_block, hvy_active, hvy_n)

  implicit none

    !> user defined parameter structure
    type (type_params), intent(inout)    :: params
    !> light data array
    integer(kind=ik), intent(inout)      :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)         :: hvy_block(:, :, :, :, :)
    !> list of active blocks (light data)
    integer(kind=ik), intent(inout)      :: hvy_active(:,:)
    !> number of heavy and light active blocks
    integer(kind=ik), intent(inout)      :: hvy_n(:)

    ! loop variable
    integer(kind=ik)                     :: k, g, Bs(1:3), hvy_id, lgt_id
    ! origin and spacing of blocks
    real(kind=rk)                        :: x0(1:3), dx(1:3)
    integer(kind=2)       :: surface(3)=0

    Bs = params%Bs
    g  = params%n_ghosts


    !---------------------------------------------------------------------------
    ! on the grid, evaluate the initial condition
    !---------------------------------------------------------------------------
    ! loop over my active heavy data
    do k = 1, hvy_n(tree_ID_flow)
        ! hvy_id of the block we're looking at
        hvy_id = hvy_active(k, tree_ID_flow)

        ! light id of this block
        call hvy_id_to_lgt_id( lgt_id, hvy_id, params%rank, params%number_blocks )
        
        ! compute block spacing and origin from treecode
        call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

        ! get block spacing for RHS
        call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )
        if ( .not. All(params%periodic_BC) ) then
            ! check if block is adjacent to a boundary of the domain, if this is the case we use one sided stencils
            call get_adjacent_boundary_surface_normal(params, lgt_id, lgt_block, params%max_treelevel, surface)
        endif
        
        ! set the initial condition on this block
        call INICOND_meta(params%physics_type, 0.0_rk, hvy_block(:,:,:,:,hvy_id), g, &
        x0, dx, boundary_flag=surface)
    enddo

end subroutine set_inicond_blocks
