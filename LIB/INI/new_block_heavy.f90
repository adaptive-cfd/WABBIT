! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: new_block_heavy.f90
! version: 0.5
! author: msr
!
! write heavy block data in first datafield, reset all other datafields
!
! input:    - params
!           - heavy data id
!           - data for second datafield
!           - coordinate vectors (store in first datafield)
! output:   - heavy block data array
!
! = log ======================================================================================
!
! 07/11/16 - switch to v0.4
! 26/01/17 - use 2D/3D hvy data array
!
! ********************************************************************************************

subroutine new_block_heavy( params, hvy_block, heavy_id, phi, coord_x, coord_y, coord_z)

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(in)      :: params

    ! heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)

    ! heavy data id
    integer(kind=ik), intent(in)        :: heavy_id

    ! input data
    real(kind=rk), intent(in)           :: phi(:, :, :)

    ! coordinate vectors
    real(kind=rk), intent(in)           :: coord_x(:), coord_y(:), coord_z(:)

    ! grid parameter
    integer(kind=ik)                    :: Bs, g
    ! number of datafields
    integer(kind=ik)                    :: dF

    ! loop variable
    integer(kind=ik)                    :: k

!---------------------------------------------------------------------------------------------
! variables initialization

    ! set parameters for readability
    Bs              = params%number_block_nodes
    g               = params%number_ghost_nodes
    dF              = params%number_data_fields

!---------------------------------------------------------------------------------------------
! main body

    ! save coordinates in field 1 of heavy data array
    hvy_block( 1, 1:Bs, 1, 1, heavy_id )              = coord_x
    hvy_block( 2, 1:Bs, 1, 1, heavy_id )              = coord_y
    hvy_block( 3, 1:Bs, 1, 1, heavy_id )              = coord_z

    if ( params%threeD_case ) then
        ! 3D:
        ! save data in first datafield (second field)
        hvy_block( g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, 2, heavy_id )      = phi(:, :, :)

    else
        ! 2D:
        ! save data in first datafield (second field)
        hvy_block( g+1:Bs+g, g+1:Bs+g, 1, 2, heavy_id )             = phi(:, :, 1)

    end if

    ! reset all other datafields
    do k = 3, dF
        hvy_block( :, :, :, k, heavy_id )                          = 0.0_rk
    end do

end subroutine new_block_heavy
