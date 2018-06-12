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
subroutine set_inicond_blocks(params, lgt_block, hvy_block, hvy_active, hvy_n, inicond, hvy_work, adapting)

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
    real(kind=rk), intent(inout)         :: hvy_work(:, :, :, :, :)
    !> list of active blocks (light data)
    integer(kind=ik), intent(inout)      :: hvy_active(:)
    !> number of heavy and light active blocks
    integer(kind=ik), intent(inout)      :: hvy_n
    !> what function to use
    character(len=*), intent(in)         :: inicond
    ! loop variable
    integer(kind=ik)                     :: k
    integer(kind=ik)                     :: Bs, g
    integer(kind=ik)                     :: hvy_id, lgt_id
    ! origin and spacing of blocks
    real(kind=rk)                        :: x0(1:3), dx(1:3)
    ! are we still adapting the grid? (for ACM & VPM)
    logical, intent(in)                  :: adapting

  !---------------------------------------------------------------------------------------------
  ! variables initialization
    Bs    = params%number_block_nodes
    g     = params%number_ghost_nodes

  !---------------------------------------------------------------------------------------------
  ! main body

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
            x0, dx, hvy_work(:,:,:,:,hvy_id), adapting)
    enddo

end subroutine set_inicond_blocks
