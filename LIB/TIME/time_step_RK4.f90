! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: time_step_RK4.f90
! version: 0.4
! author: msr
!
! time step main function, RK4
!
! input:    - time variable, params, light and heavy data, neighbor list
! output:   - time variable and heavy data array
!
! = log ======================================================================================
!
! 08/11/16 - switch to v0.4
! ********************************************************************************************

subroutine time_step_RK4( time, params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! time varible
    real(kind=rk), intent(inout)        :: time

    ! user defined parameter structure
    type (type_params), intent(in)      :: params
    ! light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    ! heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :)
    ! heavy data array - neifghbor data
    integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)

    ! list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    ! number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    ! grid parameter
    integer(kind=ik)                    :: Bs, g
    ! loop variables
    integer(kind=ik)                    :: k, dF, N_dF

    ! time step, dx
    real(kind=rk)                       :: dt, dx, my_dx

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank

    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0, sub_t1

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    N_dF  = params%number_data_fields

    ! grid parameter
    Bs    = params%number_block_nodes
    g     = params%number_ghost_nodes

    ! reset dx
    my_dx = 9.0e9_rk
    dx    = 9.0e9_rk
    ! reset dt
    dt    = 9.0e9_rk

!---------------------------------------------------------------------------------------------
! main body

    ! start time
    sub_t0 = MPI_Wtime()

    ! determinate process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    ! ----------------------------------------------------------------------------------------
    ! calculate time step
    ! loop over all active blocks (heavy data)
    do k = 1, hvy_n
        my_dx = min(my_dx, hvy_block(1, 2, 1, hvy_active(k) ) - hvy_block(1, 1, 1, hvy_active(k) ) )
    end do

    ! synchronize dx
    call MPI_Allreduce(my_dx, dx, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)

    ! calculate time step, loop over all data fields
    do dF = 2, N_dF+1
        dt = min(dt, params%CFL * dx / norm2( params%u0( (dF-2)*2 + 1 : (dF-2)*2 + 2 ) ) )
    end do

    time = time + dt
    ! last timestep should fit in maximal time
    if (time >= params%time_max) then
        time = time - dt
        dt = params%time_max - time
        time = params%time_max
    end if

    ! loop over all datafields
    do dF = 2, N_dF+1

        !------------------------------
        ! first stage
        ! synchronize ghostnodes
        call synchronize_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

        ! loop over all active heavy data blocks
        do k = 1, hvy_n

            ! save old data
            hvy_block( :, :, N_dF+2, hvy_active(k) ) = hvy_block( :, :, dF, hvy_active(k) )
            ! set k1 step
            hvy_block( :, :, N_dF+3, hvy_active(k) ) = hvy_block( :, :, dF, hvy_active(k) )
            ! RHS
            call RHS_2D_convection_diffusion( hvy_block( :, :, N_dF+3, hvy_active(k) ), &
                                              abs(hvy_block( 1, 2, 1, hvy_active(k) ) - hvy_block( 1, 1, 1, hvy_active(k) )), &
                                              abs(hvy_block( 2, 2, 1, hvy_active(k)) - hvy_block( 2, 1, 1, hvy_active(k) )), &
                                              g, Bs, &
                                              params%u0( (dF-2)*2 + 1 ), params%u0( (dF-2)*2 + 2 ), params%nu( (dF-1) ), &
                                              params%order_discretization  )

        end do

        !------------------------------
        ! second stage
        do k = 1, hvy_n

            ! save old data
            hvy_block( :, :, dF, hvy_active(k) ) = hvy_block( :, :, N_dF+2, hvy_active(k) ) + (0.5_rk * dt) * hvy_block( :, :, N_dF+3, hvy_active(k) )

        end do

        ! synchronize ghostnodes
        call synchronize_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )



        do k = 1, hvy_n

            ! set k2 step
            hvy_block( :, :, N_dF+4, hvy_active(k) ) = hvy_block( :, :, dF, hvy_active(k) )
            ! RHS
            call RHS_2D_convection_diffusion( hvy_block( :, :, N_dF+4, hvy_active(k) ), &
                                              abs(hvy_block( 1, 2, 1, hvy_active(k) ) - hvy_block( 1, 1, 1, hvy_active(k) )), &
                                              abs(hvy_block( 2, 2, 1, hvy_active(k) ) - hvy_block( 2, 1, 1, hvy_active(k) )), &
                                              g, Bs, &
                                              params%u0( (dF-2)*2 + 1 ), params%u0( (dF-2)*2 + 2 ), params%nu( (dF-1) ), &
                                              params%order_discretization  )

        end do

        !------------------------------
        ! third stage
        do k = 1, hvy_n

            ! save old data
            hvy_block( :, :, dF, hvy_active(k) ) = hvy_block( :, :, N_dF+2, hvy_active(k) ) + (0.5_rk * dt) * hvy_block( :, :, N_dF+4, hvy_active(k) )

        end do

        ! synchronize ghostnodes
        call synchronize_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

        do k = 1, hvy_n

            ! set k3 step
            hvy_block( :, :, N_dF+5, hvy_active(k) ) = hvy_block( :, :, dF, hvy_active(k) )
            ! RHS
            call RHS_2D_convection_diffusion( hvy_block( :, :, N_dF+5, hvy_active(k) ), &
                                              abs(hvy_block( 1, 2, 1, hvy_active(k) ) - hvy_block( 1, 1, 1, hvy_active(k) )), &
                                              abs(hvy_block( 2, 2, 1, hvy_active(k) ) - hvy_block( 2, 1, 1, hvy_active(k) )), &
                                              g, Bs, &
                                              params%u0( (dF-2)*2 + 1 ), params%u0( (dF-2)*2 + 2 ), params%nu( (dF-1) ), &
                                              params%order_discretization  )

        end do

        !------------------------------
        ! fourth stage
        do k = 1, hvy_n

            ! save old data
            hvy_block( :, :, dF, hvy_active(k) ) = hvy_block( :, :, N_dF+2, hvy_active(k) ) + (0.5_rk * dt) * hvy_block( :, :, N_dF+5, hvy_active(k) )

        end do

        ! synchronize ghostnodes
        call synchronize_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

        do k = 1, hvy_n

            ! set k4 step
            hvy_block( :, :, N_dF+6, hvy_active(k) ) = hvy_block( :, :, dF, hvy_active(k) )
            ! RHS
            call RHS_2D_convection_diffusion( hvy_block( :, :, N_dF+6, hvy_active(k) ), &
                                              abs(hvy_block( 1, 2, 1, hvy_active(k) ) - hvy_block( 1, 1, 1, hvy_active(k) )), &
                                              abs(hvy_block( 2, 2, 1, hvy_active(k) ) - hvy_block( 2, 1, 1, hvy_active(k) )), &
                                              g, Bs, &
                                              params%u0( (dF-2)*2 + 1 ), params%u0( (dF-2)*2 + 2 ), params%nu( (dF-1) ), &
                                              params%order_discretization  )

        end do

        !------------------------------
        ! final stage
        do k = 1, hvy_n

            ! final step
            hvy_block( :, :, dF, hvy_active(k) ) = hvy_block( :, :, N_dF+2, hvy_active(k) ) &
                                                + (dt/6.0_rk) * ( hvy_block( :, :, N_dF+3, hvy_active(k) ) &
                                                + 2.0_rk * hvy_block( :, :, N_dF+4, hvy_active(k) ) &
                                                + 2.0_rk * hvy_block( :, :, N_dF+5, hvy_active(k) ) &
                                                + hvy_block( :, :, N_dF+6, hvy_active(k) ) )

        end do

    end do

    ! end time
    sub_t1 = MPI_Wtime()
    ! write time
    if ( params%debug ) then
        ! find free or corresponding line
        k = 1
        do while ( debug%name_comp_time(k) /= "---" )
            ! entry for current subroutine exists
            if ( debug%name_comp_time(k) == "time_step" ) exit
            k = k + 1
        end do
        ! write time
        debug%name_comp_time(k) = "time_step"
        debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1
        debug%comp_time(k, 2)   = debug%comp_time(k, 2) + sub_t1 - sub_t0
    end if

end subroutine time_step_RK4
