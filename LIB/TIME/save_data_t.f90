!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name final_stage_RK.f90
!> \version 0.5
!> \author sm
!
!> \brief Save data at time t for Runge-Kutta time step routine
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

subroutine save_data_t(params, hvy_work, hvy_block, hvy_active, hvy_n)



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
    integer(kind=ik)                    :: dF, N_dF, k

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    N_dF  = params%number_data_fields

!---------------------------------------------------------------------------------------------
! main body


    select case (params%physics_type)
        case('2D_convection_diffusion')
            ! loop over all data fields
            do dF = 1, N_dF
                ! loop over all active heavy data blocks
                do k = 1, hvy_n
                    hvy_work( :, :, :, (dF-1)*5+1, hvy_active(k) ) = hvy_block( :, :, :, dF, hvy_active(k) )
                end do
            end do

        case('2D_navier_stokes')
            ! loop over all active heavy data blocks
            do k = 1, hvy_n
                hvy_work( :, :, :, 1:N_dF, hvy_active(k) ) = hvy_block( :, :, :, 1:N_dF, hvy_active(k) )
            end do

        case('3D_convection_diffusion')
            ! loop over all datafields
            do dF = 1, N_dF
                ! loop over all active heavy data blocks
                do k = 1, hvy_n
                     hvy_work( :, :, :, (dF-1)*5+1, hvy_active(k) ) = hvy_block( :, :, :, dF, hvy_active(k) )
                end do
            end do

        case('3D_navier_stokes')
            ! loop over all active heavy data blocks
            do k = 1, hvy_n
                hvy_work( :, :, :, 1:N_dF, hvy_active(k) ) = hvy_block( :, :, :, 1:N_dF, hvy_active(k) )
            end do

        case('2D_advection')
            ! loop over all datafields
            do dF = 1, N_dF
                ! loop over all active heavy data blocks
                do k = 1, hvy_n
                    hvy_work( :, :, :, (dF-1)*5+1, hvy_active(k) ) = hvy_block( :, :, :, dF, hvy_active(k) )
                end do
            end do

        case('2D_acm')
            ! loop over all active heavy data blocks
            do k = 1, hvy_n
                hvy_work( :, :, :, 1:N_dF, hvy_active(k) ) = hvy_block( :, :, :, 1:N_dF, hvy_active(k) )
            end do

        case default
            write(*,'(80("_"))')
            write(*,*) "ERROR: physics type is unknown"
            write(*,*) params%physics_type
            stop
    end select

end subroutine save_data_t
