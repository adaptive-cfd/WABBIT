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
subroutine set_inicond_blocks(params, lgt_block, hvy_block, hvy_active, hvy_n, lgt_active, &
    lgt_n, lgt_sortednumlist, hvy_mask, adapting, hvy_neighbor)

  !---------------------------------------------------------------------------------------------
  ! variables

  implicit none

    !> user defined parameter structure
    type (type_params), intent(inout)    :: params
    !> light data array
    integer(kind=ik), intent(inout)      :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)         :: hvy_block(:, :, :, :, :)
    !> heavy data array - work data
    real(kind=rk), intent(inout)         :: hvy_mask(:, :, :, :, :)
    !> list of active blocks (light data)
    integer(kind=ik), intent(inout)      :: hvy_active(:)
    !> number of heavy and light active blocks
    integer(kind=ik), intent(inout)      :: hvy_n, lgt_n
    !> list of active blocks light data)
    integer(kind=ik), intent(inout)      :: lgt_active(:)
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), intent(inout)   :: lgt_sortednumlist(:,:)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)


    ! loop variable
    integer(kind=ik)                     :: k, g, Bs(1:3), hvy_id, lgt_id
    ! origin and spacing of blocks
    real(kind=rk)                        :: x0(1:3), dx(1:3)
    ! are we still adapting the grid? (for ACM & VPM)
    logical, intent(in)                  :: adapting

    Bs    = params%Bs
    g     = params%n_ghosts

    ! it is often desirable to obtain an initial grid suitable for the penalization
    ! term, i.e. the initial grid is on the finest level where the gradient of the
    ! mask function is nonzero. Even if no mask is used or this behavior is not implemented
    ! in the physics module, the mask is always generated here (even if it is simply cnstant 0.0)
    ! because (a) that is cheap and done only once and (b) it was the default for a long time.
    call create_mask_tree(params, 0.0_rk, lgt_block, hvy_mask, &
    hvy_neighbor, hvy_active, hvy_n, lgt_active, lgt_n, lgt_sortednumlist)

    !---------------------------------------------------------------------------
    ! on the grid, evaluate the initial condition
    !---------------------------------------------------------------------------
    ! loop over my active heavy data
    do k = 1, hvy_n
        ! hvy_id of the block we're looking at
        hvy_id = hvy_active(k)

        ! light id of this block
        call hvy_id_to_lgt_id( lgt_id, hvy_id, params%rank, params%number_blocks )
        ! compute block spacing and origin from treecode
        call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

        ! set the initial condition on this block
        call INICOND_meta(params%physics_type, 0.0_rk, hvy_block(:,:,:,:,hvy_id), g, &
            x0, dx, hvy_mask(:,:,:,:,hvy_id), adapting)
    enddo

end subroutine set_inicond_blocks
