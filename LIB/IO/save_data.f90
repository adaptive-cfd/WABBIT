!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name save_data.f90
!> \version 0.5
!> \author msr
!
!> \brief save data main function, call write data routine
!
!>
!! input:
!!           - time loop parameter
!!           - parameter array
!!           - light data array
!!           - heavy data array
!!
!! output:
!!           -
!!
!!
!! = log ======================================================================================
!! \n
!! 07/11/16 - switch to v0.4 \n
!! 26/01/17 - switch to 3D, v0.5
!
! ********************************************************************************************

subroutine save_data(iteration, time, params, lgt_block, hvy_block, lgt_active, lgt_n, hvy_n )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> time loop parameters
    real(kind=rk), intent(in)                       :: time
    integer(kind=ik), intent(in)                    :: iteration

    !> user defined parameter structure
    type (type_params), intent(in)                  :: params
    !> light data array
    integer(kind=ik), intent(inout)                 :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)                    :: hvy_block(:, :, :, :, :)
    !> list of active blocks (light data)
    integer(kind=ik), intent(inout)                 :: lgt_active(:)
    !> number of active blocks (light/heavy data)
    integer(kind=ik), intent(inout)                 :: lgt_n, hvy_n

    ! loop variable
    integer(kind=ik)                                :: k
    ! file name
    character(len=80)                               :: fname
    ! cpu time variables for running time calculation
    real(kind=rk)                                   :: t0

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    ! start time
    t0 = MPI_Wtime()

    do k = 1, params%number_data_fields

        ! file name depends on variable names
        select case(params%physics_type)
            case('2D_convection_diffusion')
                ! select corresponding datafield name
                write( fname,'(a, "_", i12.12, ".h5")') trim(adjustl(params%physics%names(k))), nint(time * 1.0e6_rk)
            case('2D_navier_stokes')
                ! select corresponding datafield name
                write( fname,'(a, "_", i12.12, ".h5")') trim(adjustl(params%physics_ns%names(k))), nint(time * 1.0e6_rk)
            case('3D_convection_diffusion')
                ! select corresponding datafield name
                write( fname,'(a, "_", i12.12, ".h5")') trim(adjustl(params%physics%names(k))), nint(time * 1.0e6_rk)
            case('3D_navier_stokes')
                ! select corresponding datafield name
                write( fname,'(a, "_", i12.12, ".h5")') trim(adjustl(params%physics_ns%names(k))), nint(time * 1.0e6_rk)
            case('2D_advection')
                ! select corresponding datafield name
                write( fname,'(a, "_", i12.12, ".h5")') trim(adjustl(params%physics%names(k))), nint(time * 1.0e6_rk)
            case('3D_advection')
                ! select corresponding datafield name
                write( fname,'(a, "_", i12.12, ".h5")') trim(adjustl(params%physics%names(k))), nint(time * 1.0e6_rk)
            case('2D_acm')
                ! select corresponding datafield name
                write( fname,'(a, "_", i12.12, ".h5")') trim(adjustl(params%physics_acm%names(k))), nint(time * 1.0e6_rk)
            case('3D_acm')
                ! select corresponding datafield name
                write( fname,'(a, "_", i12.12, ".h5")') trim(adjustl(params%physics_acm%names(k))), nint(time * 1.0e6_rk)
            case default
                write(*,'(80("_"))')
                write(*,*) "ERROR: physics type is unknown - can not save data"
                stop
        end select

        call write_field( fname, time, iteration, k, params, lgt_block, hvy_block, lgt_active, lgt_n, hvy_n)
    end do

    ! timing
    call toc( params, "save_data", MPI_wtime()-t0 )
end subroutine save_data
