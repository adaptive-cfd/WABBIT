!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name reset_ghost_nodes.f90
!> \version 0.5
!> \author msr
!
!> \brief reset ghosts nodes for debuging
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

    ! loop variables
    integer(kind=ik)                    :: k, dF

    ! grid parameter
    integer(kind=ik)                    :: g, Bs

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! grid parameter
    Bs = params%number_block_nodes
    g  = params%number_ghost_nodes

!---------------------------------------------------------------------------------------------
! main body

    ! loop over all active blocks
    do k = 1, hvy_n
        ! loop over all datafields
        do dF = 2, params%number_data_fields+1
            ! reset ghost nodes
            if ( params%threeD_case ) then
                ! 3D:
                hvy_block(1:g, :, :, dF, hvy_active(k) )           = 99.0_rk!9.0e9_rk
                hvy_block(Bs+g+1:Bs+2*g, :, :, dF, hvy_active(k) ) = 99.0_rk!9.0e9_rk
                hvy_block(:, 1:g, :, dF, hvy_active(k) )           = 99.0_rk!9.0e9_rk
                hvy_block(:, Bs+g+1:Bs+2*g, :, dF, hvy_active(k) ) = 99.0_rk!9.0e9_rk
                hvy_block(:, :, 1:g, dF, hvy_active(k) )           = 99.0_rk!9.0e9_rk
                hvy_block(:, :, Bs+g+1:Bs+2*g, dF, hvy_active(k) ) = 99.0_rk!9.0e9_rk
            else
                ! 2D:
                hvy_block(1:g, :, 1, 2, hvy_active(k) )           = 9.0e9_rk
                hvy_block(Bs+g+1:Bs+2*g, :, 1, 2, hvy_active(k) ) = 9.0e9_rk
                hvy_block(:, 1:g, 1, 2, hvy_active(k) )           = 9.0e9_rk
                hvy_block(:, Bs+g+1:Bs+2*g, 1, 2, hvy_active(k) ) = 9.0e9_rk
            end if
        end do
    end do


end subroutine reset_ghost_nodes
