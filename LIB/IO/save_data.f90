! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: save_data.f90
! version: 0.5
! author: msr
!
! save data main function, call write data routine
!
! input:    - time loop parameter
!           - parameter array
!           - light data array
!           - heavy data array
! output:   -
!
! = log ======================================================================================
!
! 07/11/16 - switch to v0.4
! 26/01/17 - switch to 3D, v0.5
!
! ********************************************************************************************

subroutine save_data(iteration, time, params, lgt_block, hvy_block, lgt_active, lgt_n, hvy_n)

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! time loop parameters
    real(kind=rk), intent(in)                       :: time
    integer(kind=ik), intent(in)                    :: iteration

    ! user defined parameter structure
    type (type_params), intent(in)                  :: params
    ! light data array
    integer(kind=ik), intent(in)                    :: lgt_block(:, :)
    ! heavy data array - block data
    real(kind=rk), intent(in)                       :: hvy_block(:, :, :, :, :)
    ! list of active blocks (light data)
    integer(kind=ik), intent(inout)                 :: lgt_active(:)
    ! number of active blocks (light data)
    integer(kind=ik), intent(inout)                 :: lgt_n
    ! number of active blocks (heavy data)
    integer(kind=ik), intent(inout)                 :: hvy_n

    ! loop variable
    integer(kind=ik)                                :: k
    ! file name
    character(len=80)                               :: fname
    ! cpu time variables for running time calculation
    real(kind=rk)                                   :: sub_t0, sub_t1

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    ! start time
    sub_t0 = MPI_Wtime()

    ! FIXME DF
    ! real datafields start at datafield 2
    do k = 2, params%number_data_fields+1

        ! file name depends on variable names
        select case(params%physics_type)
            case('2D_convection_diffusion')
                ! select corresponding datafield name
                write( fname,'(a, "_", i12.12, ".h5")') trim(adjustl(params%physics%names(k-1))), nint(time * 1.0e6_rk)
            case('2D_navier_stokes')
                ! select corresponding datafield name
                write( fname,'(a, "_", i12.12, ".h5")') trim(adjustl(params%physics_ns%names(k-1))), nint(time * 1.0e6_rk)
            case('3D_convection_diffusion')
                ! select corresponding datafield name
                write( fname,'(a, "_", i12.12, ".h5")') trim(adjustl(params%physics%names(k-1))), nint(time * 1.0e6_rk)
            case('3D_navier_stokes')
                ! select corresponding datafield name
                write( fname,'(a, "_", i12.12, ".h5")') trim(adjustl(params%physics_ns%names(k-1))), nint(time * 1.0e6_rk)
            case('2D_advection')
                ! select corresponding datafield name
                write( fname,'(a, "_", i12.12, ".h5")') trim(adjustl(params%physics%names(k-1))), nint(time * 1.0e6_rk)
            case default
                write(*,'(80("_"))')
                write(*,*) "ERROR: physics type is unknown - can not save data"
                stop
        end select

        call write_field( fname, time, iteration, k, params, lgt_block, hvy_block, lgt_active, lgt_n, hvy_n)
    end do

    ! end time
    sub_t1 = MPI_Wtime()

    ! write time
    if ( params%debug ) then
        ! find free or corresponding line
        k = 1
        do while ( debug%name_comp_time(k) /= "---" )
            ! entry for current subroutine exists
            if ( debug%name_comp_time(k) == "save_data" ) exit
            k = k + 1
        end do
        ! write time
        debug%name_comp_time(k) = "save_data"
        debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1
        debug%comp_time(k, 2)   = debug%comp_time(k, 2) + sub_t1 - sub_t0
    end if

end subroutine save_data
