!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name get_inicond_from_file.f90
!> \version 0.5
!> \author engels, sm
!
!> \brief read initial condition from a file
!
!>
!! input:
!!           - 
!!           - 
!!           - 
!!           - 
!!           - 
!!
!! output:
!!           -
!!
!!
!! = log ======================================================================================
!! \n
!! 28/09/17 - create
!
! ********************************************************************************************

subroutine get_inicond_from_file( params, )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none
    !> user defined parameter structure
    type (type_params), intent(in)      :: params
        !> light data array
    integer(kind=ik), intent(inout)      :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)         :: hvy_block(:, :, :, :, :)
    !> list of active blocks (light data)
    integer(kind=ik), intent(inout)      :: hvy_active(:)
    !> list of active blocks light data)
    integer(kind=ik), intent(inout)      :: lgt_active(:)
    !> number of heavy and light active blocks
    integer(kind=ik), intent(inout)      :: hvy_n, lgt_n

!---------------------------------------------------------------------------------------------
! variables initialization

    N_dF = params%number_data_fields

!---------------------------------------------------------------------------------------------
! main body

    call fetch_attributes(params%ini_files(1), params, lgt_n, time, iteration, )

    call read_mesh(params%ini_files(1), params, lgt_n, lgt_block, hvy_n )

    !create hvy_active and lgt_active list
    call create_active_and_sorted_lists(params, lgt_block, lgt_active, lgt_n, hvy_active, hvy_n, lgt_sortednumlist, create_sorted_list)

    do dF = 1, N_dF
        call read_field(params%ini_files(dF), dF, params, hvy_block, lgt_n, hvy_n )
    end do


end subroutine get_inicond_from_file
