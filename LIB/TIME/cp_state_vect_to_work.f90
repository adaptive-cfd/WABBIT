!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name final_stage_RK.f90
!> \version 0.5
!> \author sm
!
!> \brief Copy state vector in hvy_block to hvy_work.
!
!>
!! input:
!!          - params
!!          - heavy data
!!
!! output:
!!          - heavy data
!!
!! butcher table, e.g.
!!
!! |   |    |    |   |
!! |---|----|----|---|
!! | 0 | 0  | 0  |  0|
!! |c2 | a21| 0  |  0|
!! |c3 | a31| a32|  0|
!! | 0 | b1 | b2 | b3|
!!
!!
!! = log ======================================================================================
!! \n
!! 23/05/17 - create
!
!**********************************************************************************************

subroutine cp_state_vect_to_work(params, hvy_work, hvy_block, hvy_active, hvy_n)



!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy work data array - block data
    real(kind=rk), intent(inout)        :: hvy_work(:, :, :, :, :)

    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    !> loop variables
    integer(kind=ik)                    :: dF, neq, k

!---------------------------------------------------------------------------------------------
! variables initialization

    neq  = params%number_data_fields

!---------------------------------------------------------------------------------------------
! main body

    ! loop over all active heavy data blocks
    do k = 1, hvy_n
        hvy_work( :, :, :, 1:neq, hvy_active(k) ) = hvy_block( :, :, :, 1:neq, hvy_active(k) )
    end do

end subroutine cp_state_vect_to_work
