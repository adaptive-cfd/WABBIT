!> \brief reset ghosts nodes for all blocks (not just active ones) for debuging
!> ghost nodes are set to a very large constant.
!! input:    - params, light and heavy data \n
!! output:   - heavy data array
! ********************************************************************************************

subroutine reset_ghost_nodes(  params, hvy_block, hvy_active, hvy_n )
    implicit none

    type (type_params), intent(in)      :: params                     !> user defined parameter structure
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)   !> heavy data array - block data
    integer(kind=ik), intent(in)        :: hvy_active(:)              !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n                      !> number of active blocks (heavy data)
    integer(kind=ik)                    :: g                          ! grid parameter
    integer(kind=ik), dimension(3)      :: Bs

    Bs = params%Bs
    g  = params%n_ghosts

    ! delete layer of ghost nodes for all blocks (not just active ones)

    !-- x-direction
    hvy_block(1:g, :, :, :, : )           = 9.0e9_rk
    hvy_block(Bs(1)+g+1:Bs(1)+2*g, :, :, :, : ) = 9.0e9_rk
    !-- y-direction
    hvy_block(:, 1:g, :, :, : )           = 9.0e9_rk
    hvy_block(:, Bs(2)+g+1:Bs(2)+2*g, :, :, : ) = 9.0e9_rk
    !-- z-direction
    if ( params%dim == 3 ) then
      hvy_block(:, :, 1:g, :, : )           = 9.0e9_rk
      hvy_block(:, :, Bs(3)+g+1:Bs(3)+2*g, :, : ) = 9.0e9_rk
    end if
end subroutine reset_ghost_nodes
