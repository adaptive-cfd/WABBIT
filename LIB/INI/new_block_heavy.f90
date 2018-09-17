!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name new_block_heavy.f90
!> \version 0.5
!> \author msr
!
!> \brief write heavy block data in first datafield, reset all other datafields
!
!>
!! input:    
!!           - params
!!           - heavy data id
!!           - data for second datafield
!!           - coordinate vectors (store in first datafield)
!!
!! output:   
!!           - heavy block data array
!!
!! = log ======================================================================================
!! \n
!! 07/11/16 - switch to v0.4 \n
!! 26/01/17 - use 2D/3D hvy data array
!
! ********************************************************************************************

subroutine new_block_heavy( params, hvy_block, heavy_id, phi)

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params

    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)

    !> heavy data id
    integer(kind=ik), intent(in)        :: heavy_id

    !> input data
    real(kind=rk), intent(in)           :: phi(:, :, :)

    ! grid parameter
    integer(kind=ik)                    :: Bs, g
    ! number of datafields
    integer(kind=ik)                    :: dF

    ! loop variable
    integer(kind=ik)                    :: k

!---------------------------------------------------------------------------------------------
! variables initialization

    ! set parameters for readability
    Bs              = params%Bs
    g               = params%nr_ghosts
    dF              = params%n_eqn

!---------------------------------------------------------------------------------------------
! main body

    if ( params%threeD_case ) then
        ! 3D:
        ! save data in first datafield
        hvy_block( g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, 1, heavy_id )      = phi(:, :, :)

    else
        ! 2D:
        ! save data in first datafield
        hvy_block( g+1:Bs+g, g+1:Bs+g, 1, 1, heavy_id )             = phi(:, :, 1)

    end if

    ! reset all other datafields
    do k = 2, dF
        hvy_block( :, :, :, k, heavy_id )                          = 0.0_rk
    end do

end subroutine new_block_heavy
