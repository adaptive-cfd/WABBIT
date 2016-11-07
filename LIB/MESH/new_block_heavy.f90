! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: new_block_heavy.f90
! version: 0.4
! author: msr
!
! write heavy block data in first datafield, reset all other datafields
!
! input:    - heavy data id
!           - data for second datafield
!           - coordinate vectors (store in first datafield)
!           - grid parameters
!           - number of data fields
! output:   - heavy block data array
!
! = log ======================================================================================
!
! 07/11/16 - switch to v0.4
! ********************************************************************************************

subroutine new_block_heavy(block_data, heavy_id, data_, coord_x, coord_y, Bs, g, dF)

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! heavy data array - block data
    real(kind=rk), intent(inout)        :: block_data(:, :, :, :)
    ! heavy data id
    integer(kind=ik), intent(in)        :: heavy_id
    ! grid parameter
    integer(kind=ik), intent(in)        :: Bs, g
    ! input data
    real(kind=rk), intent(in)           :: data_( Bs, Bs )
    ! coordinate vectors
    real(kind=rk), intent(in)           :: coord_x(Bs), coord_y(Bs)
    ! number of datafields
    integer(kind=ik), intent(in)        :: dF

    ! loop variable
    integer(kind=ik)                    :: k

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    ! save coordinates in field 1 of heavy data array
    block_data( heavy_id, 1, 1:Bs, 1 )              = coord_x
    block_data( heavy_id, 2, 1:Bs, 1 )              = coord_y

    ! save data in first datafield
    block_data( heavy_id, g+1:Bs+g, g+1:Bs+g, 2 )   = data_
    ! reset all other datafields
    do k = 3, dF
        block_data( heavy_id, :, :, k )             = 9.0e9_rk
    end do

end subroutine new_block_heavy
