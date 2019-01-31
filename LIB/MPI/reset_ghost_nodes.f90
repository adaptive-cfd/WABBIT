!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name reset_ghost_nodes.f90
!> \version 0.5
!> \author msr
!
!> \brief reset ghosts nodes for all blocks (not just active ones) for debuging
!> ghost nodes are set to a very large constant.
!
!>
!! input:    - params, light and heavy data \n
!! output:   - heavy data array
!
!> \details
!! = log ======================================================================================
!! \n
!! 19/05/17 - create
!
! ********************************************************************************************

subroutine reset_ghost_nodes(  params, hvy_block, hvy_active, hvy_n )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    ! grid parameter
    integer(kind=ik)                    :: g
    integer(kind=ik), dimension(3)     :: Bs

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! grid parameter
    Bs = params%Bs
    g  = params%n_ghosts

!---------------------------------------------------------------------------------------------
! main body
    ! delete layer of ghost nodes for all blocks (not just active ones)

    !-- x-direction
    hvy_block(1:g, :, :, :, : )           = 9.0e9_rk
    hvy_block(Bs(1)+g+1:Bs(1)+2*g, :, :, :, : ) = 9.0e9_rk
    !-- y-direction
    hvy_block(:, 1:g, :, :, : )           = 9.0e9_rk
    hvy_block(:, Bs(2)+g+1:Bs(2)+2*g, :, :, : ) = 9.0e9_rk
    !-- z-direction
    if ( params%threeD_case ) then
      hvy_block(:, :, 1:g, :, : )           = 9.0e9_rk
      hvy_block(:, :, Bs(3)+g+1:Bs(3)+2*g, :, : ) = 9.0e9_rk
    end if
end subroutine reset_ghost_nodes
