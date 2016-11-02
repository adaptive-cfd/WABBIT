! ********************************
! WABBIT
! --------------------------------
!
! delete heavy block data
!
! name: delete_block_heavy.f90
! date: 27.10.2016
! author: msr
! version: 0.3
!
! ********************************

subroutine delete_block_heavy(heavy_id)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), intent(in)    :: heavy_id
    integer(kind=ik)                :: dF

    ! delete heavy data
    do dF = 1, blocks_params%number_data_fields
        blocks_data(heavy_id)%data_fields(dF)%data_  = 0.0_rk
    end do
    blocks_data(heavy_id)%coord_x   = 0.0_rk
    blocks_data(heavy_id)%coord_y   = 0.0_rk
    blocks_data(heavy_id)%dx        = 0.0_rk
    blocks_data(heavy_id)%dy        = 0.0_rk
    blocks_data(heavy_id)%block_id  = -1

end subroutine delete_block_heavy
