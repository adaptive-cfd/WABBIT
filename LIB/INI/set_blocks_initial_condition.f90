! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: init_data.f90
! version: 0.5
! author: msr
!
! This routine initializes the block data, i.e. it evaluates the initial condition on the grid
!
! input:    - parameter array
!           - light data array
!           - heavy data array
!           - neighbor data array
!           - light and heavy active block list
! output:   - filled user defined data structure for global params
!           - initialized light and heavy data arrays
!
! = log ======================================================================================
!
! 04/11/16 - switch to v0.4, now run complete initialization within these subroutine and return
!            initialized block data to main program
! 07/12/16 - now uses heavy work data array
! 25/01/17 - switch to 3D, v0.5
!
! ********************************************************************************************
subroutine set_blocks_initial_condition(params, lgt_block, hvy_block, hvy_work, hvy_neighbor, lgt_active, hvy_active)

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(inout)                 :: params
    ! light data array
    integer(kind=ik), intent(inout)      :: lgt_block(:, :)
    ! heavy data array - block data
    real(kind=rk), intent(inout)         :: hvy_block(:, :, :, :, :)
    ! heavy work array  )
    real(kind=rk), intent(inout)         :: hvy_work(:, :, :, :, :)
    ! neighbor array (heavy data)
    integer(kind=ik), intent(inout)      :: hvy_neighbor(:,:)
    ! list of active blocks light data)
    integer(kind=ik), intent(inout)      :: lgt_active(:)
    ! list of active blocks (light data)
    integer(kind=ik), intent(inout)      :: hvy_active(:)

    ! loop variable
    integer(kind=ik)                                :: k, allocate_error

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

    ! initial data field
    select case( params%initial_cond )
        case ("gauss-blob","gauss_blob")

            ! 2D: gaussblob, 3D: zero fields
            if ( params%threeD_case ) then
                ! 3D:
                call inicond_sphere( params, lgt_block, hvy_block )
            else
                ! 2D:
                call inicond_gauss_blob( params, lgt_block, hvy_block )
            end if

        case ("vorticity_filaments")
!            call inicond_vorticity_filaments(params, lgt_block, hvy_block)
!            ! set density and pressure
!            hvy_block( :, :, 2, : ) = 1.0_rk
!            hvy_block( :, :, 5, : ) = 1e5_rk

        case ("richtmyer_meshkov")
!            call inicond_richtmyer_meshkov(params, lgt_block, hvy_block)

        case ("shear_layer")
            call inicond_shear_layer(params, lgt_block, hvy_block)

        case default
            write(*,'(80("_"))')
            write(*,*) "ERROR: initial condition is unknown"
            write(*,*) params%initial_cond
            stop

    end select


    ! end time
    sub_t1 = MPI_Wtime()
    ! write time
    if ( params%debug ) then
        ! find free or corresponding line
        k = 1
        do while ( debug%name_comp_time(k) /= "---" )
            ! entry for current subroutine exists
            if ( debug%name_comp_time(k) == "init_data" ) exit
            k = k + 1
        end do
        ! write time
        debug%name_comp_time(k) = "init_data"
        debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1
        debug%comp_time(k, 2)   = debug%comp_time(k, 2) + sub_t1 - sub_t0

    end if

end subroutine set_blocks_initial_condition
