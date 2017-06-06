!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name time_step_RK4.f90
!> \version 0.5
!> \author msr
!
!> \brief time step main function, RK4
!
!>
!! input:    - time variable, params, light and heavy data, neighbor list \n
!! output:   - time variable and heavy data array \n
!!
!!
!! physics: \n
! --------
!> - convection/diffusion: works only for one datafield \n
!!
!!
!! = log ======================================================================================
!! \n
!! 08/11/16 - switch to v0.4 \n
!! 07/12/16 - now uses heavy work data array and work for different physics \n
!! 31/01/17 - switch to 3D, v0.5 \n
!
! ********************************************************************************************

subroutine time_step_RK4( time, params, lgt_block, hvy_block, hvy_work, hvy_neighbor, hvy_active, hvy_n )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> time varible
    real(kind=rk), intent(inout)        :: time

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy work data array - block data
    real(kind=rk), intent(inout)        :: hvy_work(:, :, :, :, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)

    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    ! grid parameter
    integer(kind=ik)                    :: Bs, g
    ! loop variables
    integer(kind=ik)                    :: k, dF, N_dF, lgt_id, d

    ! time step, dx
    real(kind=rk)                       :: dt, dx, my_dx
    ! spacing and origin of a block
    real(kind=rk)                       :: xx0(1:3), ddx(1:3)
    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank

    ! new time after timestep dt, so current time can be used in RHS
    real(kind=rk)                       :: time_dt

    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0, sub_t1, time_sum

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

    time_sum = 0.0_rk

    ! set MPI parameter
    rank         = params%rank
    d = 2
    if ( params%threeD_case) d = 3

!---------------------------------------------------------------------------------------------
! main body

    ! start time
    sub_t0 = MPI_Wtime()

    ! ----------------------------------------------------------------------------------------
    ! calculate time step
    ! loop over all active blocks (heavy data)
    ! FIXME: you could also look over light data, as ddx is available only from that. no mpi
    do k = 1, hvy_n
        ! light id of this block
        call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )
        ! compute blocks' spacing from treecode
        call get_block_spacing_origin( params, lgt_id, lgt_block, xx0, ddx )
        ! find smallest dx of active blocks
        my_dx = min(my_dx, minval(ddx(1:d)) )

        ! HACK repair first datafield, as we're about to remove it
        hvy_block(:,:,:,1,hvy_active(k)) = 0.0_rk
        hvy_block(1,2,:,1,hvy_active(k)) = ddx(1)
        hvy_block(2,2,:,1,hvy_active(k)) = ddx(2)
        hvy_block(3,2,:,1,hvy_active(k)) = ddx(3)

    end do

    ! find globally smallest dx
    call MPI_Allreduce(my_dx, dx, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)

    ! calculate dt
    call calculate_time_step(params, dx, dt)

    ! calculate value for time after one dt
    time_dt = time + dt
    ! last timestep should fit in maximal time
    if ( time_dt >= params%time_max) then
        dt = params%time_max - time
        time_dt = params%time_max
    end if

    ! time measurement without ghost nodes synchronization
    sub_t1   = MPI_Wtime()
    time_sum = time_sum + (sub_t1 - sub_t0)

    !***************************************************************************
    ! first stage
    !***************************************************************************
    ! synchronize ghostnodes
    call synchronize_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

    ! restart time
    sub_t0 = MPI_Wtime()

    ! RHS, depends on physics
    select case(params%physics_type)
        case('2D_convection_diffusion')
            ! loop over all datafields
            do dF = 2, N_dF+1
                ! loop over all active heavy data blocks
                do k = 1, hvy_n

                    ! save old data
                    hvy_work( :, :, :, (dF-2)*5+1, hvy_active(k) ) = hvy_block( :, :, :, dF, hvy_active(k) )
                    ! set k1 step
                    hvy_work( :, :, :, (dF-2)*5+2, hvy_active(k) ) = hvy_block( :, :, :, dF, hvy_active(k) )
                    ! RHS
                    call RHS_2D_convection_diffusion( hvy_work( :, :, 1, (dF-2)*5+2, hvy_active(k) ), &
                                                      abs(hvy_block( 1, 2, 1, 1, hvy_active(k) ) - hvy_block( 1, 1, 1, 1, hvy_active(k) )), &
                                                      abs(hvy_block( 2, 2, 1, 1, hvy_active(k)) - hvy_block( 2, 1, 1, 1, hvy_active(k) )), &
                                                      g, Bs, &
                                                      params%physics%u0( (dF-2)*2 + 1 ), params%physics%u0( (dF-2)*2 + 2 ), params%physics%nu( (dF-1) ), &
                                                      params%order_discretization  )

                end do
            end do

        case('2D_navier_stokes')
            ! loop over all active heavy data blocks
            do k = 1, hvy_n

                ! save old data
                hvy_work( :, :, :, 1:N_dF, hvy_active(k) ) = hvy_block( :, :, :, 2:N_dF+1, hvy_active(k) )
                ! set k1 step
                hvy_work( :, :, :, N_dF+1:2*N_dF, hvy_active(k) ) = hvy_block( :, :, :, 2:N_dF+1, hvy_active(k) )
                ! RHS
                call RHS_2D_navier_stokes( params%physics_ns, g, Bs, &
                                           abs(hvy_block( 1, 2, 1, 1, hvy_active(k) ) - hvy_block( 1, 1, 1, 1, hvy_active(k) )), &
                                           abs(hvy_block( 2, 2, 1, 1, hvy_active(k)) - hvy_block( 2, 1, 1, 1, hvy_active(k) )), &
                                           N_dF, &
                                           hvy_work( :, :, 1, N_dF+1:2*N_dF, hvy_active(k) ) )

            end do

        case('3D_convection_diffusion')
            ! loop over all datafields
            do dF = 2, N_dF+1
                ! loop over all active heavy data blocks
                do k = 1, hvy_n

                    ! save old data
                    hvy_work( :, :, :, (dF-2)*5+1, hvy_active(k) ) = hvy_block( :, :, :, dF, hvy_active(k) )
                    ! set k1 step
                    hvy_work( :, :, :, (dF-2)*5+2, hvy_active(k) ) = hvy_block( :, :, :, dF, hvy_active(k) )
                    ! RHS
                    call RHS_3D_convection_diffusion( hvy_work( :, :, :, (dF-2)*5+2, hvy_active(k) ), &
                                                      abs(hvy_block( 1, 2, 1, 1, hvy_active(k) ) - hvy_block( 1, 1, 1, 1, hvy_active(k) )), &
                                                      abs(hvy_block( 2, 2, 1, 1, hvy_active(k)) - hvy_block( 2, 1, 1, 1, hvy_active(k) )), &
                                                      abs(hvy_block( 3, 2, 1, 1, hvy_active(k)) - hvy_block( 3, 1, 1, 1, hvy_active(k) )), &
                                                      g, Bs, &
                                                      params%physics%u0( (dF-2)*2 + 1 ), params%physics%u0( (dF-2)*2 + 2 ), params%physics%u0( (dF-2)*2 + 3 ), params%physics%nu( (dF-1) ), &
                                                      params%order_discretization  )

                end do
            end do

        case('3D_navier_stokes')
            ! loop over all active heavy data blocks
            do k = 1, hvy_n

                ! save old data
                hvy_work( :, :, :, 1:N_dF, hvy_active(k) ) = hvy_block( :, :, :, 2:N_dF+1, hvy_active(k) )
                ! set k1 step
                hvy_work( :, :, :, N_dF+1:2*N_dF, hvy_active(k) ) = hvy_block( :, :, :, 2:N_dF+1, hvy_active(k) )
                ! RHS
                call RHS_3D_navier_stokes( params%physics_ns, g, Bs, &
                                           abs(hvy_block( 1, 2, 1, 1, hvy_active(k) ) - hvy_block( 1, 1, 1, 1, hvy_active(k) )), &
                                           abs(hvy_block( 2, 2, 1, 1, hvy_active(k)) - hvy_block( 2, 1, 1, 1, hvy_active(k) )), &
                                           abs(hvy_block( 3, 2, 1, 1, hvy_active(k)) - hvy_block( 3, 1, 1, 1, hvy_active(k) )), &
                                           N_dF, &
                                           hvy_work( :, :, :, N_dF+1:2*N_dF, hvy_active(k) ) )

            end do

        case('2D_advection')
            ! loop over all datafields
            do dF = 2, N_dF+1
                ! loop over all active heavy data blocks
                do k = 1, hvy_n

                    ! light id of this block
                    call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )
                    ! compute blocks' spacing from treecode
                    call get_block_spacing_origin( params, lgt_id, lgt_block, xx0, ddx )

                    ! save old data
                    hvy_work( :, :, :, (dF-2)*5+1, hvy_active(k) ) = hvy_block( :, :, :, dF, hvy_active(k) )
                    ! set k1 step
                    hvy_work( :, :, :, (dF-2)*5+2, hvy_active(k) ) = hvy_block( :, :, :, dF, hvy_active(k) )
                    ! RHS
                    call RHS_2D_advection( hvy_work( :, :, 1, (dF-2)*5+2, hvy_active(k) ), &
                                                      xx0(1:2), &
                                                      ddx(1:2), &
                                                      g, Bs, &
                                                      time, &
                                                      params%order_discretization  )

                end do
            end do

        case default
            write(*,'(80("_"))')
            write(*,*) "ERROR: physics type is unknown"
            write(*,*) params%physics_type
            stop

    end select

    !***************************************************************************
    ! second stage
    !***************************************************************************
    ! loop over all datafields
    select case(params%physics_type)
        case('2D_convection_diffusion')
            do dF = 2, N_dF+1
                do k = 1, hvy_n
                    ! save old data
                    hvy_block( :, :, :, dF, hvy_active(k) ) = hvy_work( :, :, :, (dF-2)*5+1, hvy_active(k) ) + (0.5_rk * dt) * hvy_work( :, :, :, (dF-2)*5+2, hvy_active(k) )
                end do
            end do

        case('2D_navier_stokes')
            do k = 1, hvy_n
                ! save old data
                hvy_block( :, :, :, 2:N_dF+1, hvy_active(k) ) = hvy_work( :, :, :, 1:N_dF, hvy_active(k) ) + (0.5_rk * dt) * hvy_work( :, :, :, N_dF+1:2*N_dF, hvy_active(k) )
            end do

        case('3D_convection_diffusion')
            do dF = 2, N_dF+1
                do k = 1, hvy_n
                    ! save old data
                    hvy_block( :, :, :, dF, hvy_active(k) ) = hvy_work( :, :, :, (dF-2)*5+1, hvy_active(k) ) + (0.5_rk * dt) * hvy_work( :, :, :, (dF-2)*5+2, hvy_active(k) )
                end do
            end do

        case('3D_navier_stokes')
            do k = 1, hvy_n
                ! save old data
                hvy_block( :, :, :, 2:N_dF+1, hvy_active(k) ) = hvy_work( :, :, :, 1:N_dF, hvy_active(k) ) + (0.5_rk * dt) * hvy_work( :, :, :, N_dF+1:2*N_dF, hvy_active(k) )
            end do

        case('2D_advection')
            do dF = 2, N_dF+1
                do k = 1, hvy_n
                    ! save old data
                    hvy_block( :, :, :, dF, hvy_active(k) ) = hvy_work( :, :, :, (dF-2)*5+1, hvy_active(k) ) + (0.5_rk * dt) * hvy_work( :, :, :, (dF-2)*5+2, hvy_active(k) )
                end do
            end do

        case default
            write(*,'(80("_"))')
            write(*,*) "ERROR: physics type is unknown"
            write(*,*) params%physics_type
            stop

    end select

    ! time measurement without ghost nodes synchronization
    sub_t1   = MPI_Wtime()
    time_sum = time_sum + (sub_t1 - sub_t0)

    ! synchronize ghostnodes
    call synchronize_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

    ! restart time
    sub_t0 = MPI_Wtime()

    ! RHS, depends on physics
    select case(params%physics_type)
        case('2D_convection_diffusion')
            ! loop over all datafields
            do dF = 2, N_dF+1
                ! loop over all active heavy data blocks
                do k = 1, hvy_n

                    ! set k2 step
                    hvy_work( :, :, :, (dF-2)*5+3, hvy_active(k) ) = hvy_block( :, :, :, dF, hvy_active(k) )
                    ! RHS
                    call RHS_2D_convection_diffusion( hvy_work( :, :, 1, (dF-2)*5+3, hvy_active(k) ), &
                                                      abs(hvy_block( 1, 2, 1, 1, hvy_active(k) ) - hvy_block( 1, 1, 1, 1, hvy_active(k) )), &
                                                      abs(hvy_block( 2, 2, 1, 1, hvy_active(k) ) - hvy_block( 2, 1, 1, 1, hvy_active(k) )), &
                                                      g, Bs, &
                                                      params%physics%u0( (dF-2)*2 + 1 ), params%physics%u0( (dF-2)*2 + 2 ), params%physics%nu( (dF-1) ), &
                                                      params%order_discretization  )

                end do
            end do

        case('2D_navier_stokes')
            ! loop over all active heavy data blocks
            do k = 1, hvy_n

                ! set k2 step
                hvy_work( :, :, :, 2*N_dF+1:3*N_dF, hvy_active(k) ) = hvy_block( :, :, :, 2:N_dF+1, hvy_active(k) )
                ! RHS
                call RHS_2D_navier_stokes( params%physics_ns, g, Bs, &
                                           abs(hvy_block( 1, 2, 1, 1, hvy_active(k) ) - hvy_block( 1, 1, 1, 1, hvy_active(k) )), &
                                           abs(hvy_block( 2, 2, 1, 1, hvy_active(k)) - hvy_block( 2, 1, 1, 1, hvy_active(k) )), &
                                           N_dF, &
                                           hvy_work( :, :, 1, 2*N_dF+1:3*N_dF, hvy_active(k) ) )

            end do

        case('3D_convection_diffusion')
            ! loop over all datafields
            do dF = 2, N_dF+1
                ! loop over all active heavy data blocks
                do k = 1, hvy_n

                    ! set k2 step
                    hvy_work( :, :, :, (dF-2)*5+3, hvy_active(k) ) = hvy_block( :, :, :, dF, hvy_active(k) )
                    ! RHS
                    call RHS_3D_convection_diffusion( hvy_work( :, :, :, (dF-2)*5+3, hvy_active(k) ), &
                                                      abs(hvy_block( 1, 2, 1, 1, hvy_active(k) ) - hvy_block( 1, 1, 1, 1, hvy_active(k) )), &
                                                      abs(hvy_block( 2, 2, 1, 1, hvy_active(k) ) - hvy_block( 2, 1, 1, 1, hvy_active(k) )), &
                                                      abs(hvy_block( 3, 2, 1, 1, hvy_active(k) ) - hvy_block( 3, 1, 1, 1, hvy_active(k) )), &
                                                      g, Bs, &
                                                      params%physics%u0( (dF-2)*2 + 1 ), params%physics%u0( (dF-2)*2 + 2 ), params%physics%u0( (dF-2)*2 + 3 ), params%physics%nu( (dF-1) ), &
                                                      params%order_discretization  )

                end do
            end do

        case('3D_navier_stokes')
            ! loop over all active heavy data blocks
            do k = 1, hvy_n

                ! set k2 step
                hvy_work( :, :, :, 2*N_dF+1:3*N_dF, hvy_active(k) ) = hvy_block( :, :, :, 2:N_dF+1, hvy_active(k) )
                ! RHS
                call RHS_3D_navier_stokes( params%physics_ns, g, Bs, &
                                           abs(hvy_block( 1, 2, 1, 1, hvy_active(k) ) - hvy_block( 1, 1, 1, 1, hvy_active(k) )), &
                                           abs(hvy_block( 2, 2, 1, 1, hvy_active(k)) - hvy_block( 2, 1, 1, 1, hvy_active(k) )), &
                                           abs(hvy_block( 3, 2, 1, 1, hvy_active(k)) - hvy_block( 3, 1, 1, 1, hvy_active(k) )), &
                                           N_dF, &
                                           hvy_work( :, :, :, 2*N_dF+1:3*N_dF, hvy_active(k) ) )

            end do

        case('2D_advection')
            ! loop over all datafields
            do dF = 2, N_dF+1
                ! loop over all active heavy data blocks
                do k = 1, hvy_n

                    ! light id of this block
                    call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )
                    ! compute blocks' spacing from treecode
                    call get_block_spacing_origin( params, lgt_id, lgt_block, xx0, ddx )

                    ! set k2 step
                    hvy_work( :, :, :, (dF-2)*5+3, hvy_active(k) ) = hvy_block( :, :, :, dF, hvy_active(k) )
                    ! RHS
!                    call RHS_2D_advection( hvy_work( :, :, 1, (dF-2)*5+3, hvy_active(k) ), &
!                                                      xx0(1:2), &
!                                                      ddx(1:2), &
!                                                      g, Bs, &
!                                                      time, &
!                                                      params%order_discretization  )
                    call RHS_2D_advection( hvy_work( :, :, 1, (dF-2)*5+3, hvy_active(k) ), &
                                                      xx0(1:2), &
                                                      ddx(1:2), &
                                                      g, Bs, &
                                                      time+dt/2.0_rk, &
                                                      params%order_discretization  )

                end do
            end do

        case default
            write(*,'(80("_"))')
            write(*,*) "ERROR: physics type is unknown"
            write(*,*) params%physics_type
            stop

    end select

    !***************************************************************************
    ! third stage
    !***************************************************************************
    ! loop over all datafields
    select case(params%physics_type)
        case('2D_convection_diffusion')
            do dF = 2, N_dF+1
                do k = 1, hvy_n
                    ! save old data
                    hvy_block( :, :, :, dF, hvy_active(k) ) = hvy_work( :, :, :, (dF-2)*5+1, hvy_active(k) ) + (0.5_rk * dt) * hvy_work( :, :, :, (dF-2)*5+3, hvy_active(k) )
                end do
            end do

        case('2D_navier_stokes')
            do k = 1, hvy_n
                ! save old data
                hvy_block( :, :, :, 2:N_dF+1, hvy_active(k) ) = hvy_work( :, :, :, 1:N_dF, hvy_active(k) ) + (0.5_rk * dt) * hvy_work( :, :, :, 2*N_dF+1:3*N_dF, hvy_active(k) )
            end do

        case('3D_convection_diffusion')
            do dF = 2, N_dF+1
                do k = 1, hvy_n
                    ! save old data
                    hvy_block( :, :, :, dF, hvy_active(k) ) = hvy_work( :, :, :, (dF-2)*5+1, hvy_active(k) ) + (0.5_rk * dt) * hvy_work( :, :, :, (dF-2)*5+3, hvy_active(k) )
                end do
            end do

        case('3D_navier_stokes')
            do k = 1, hvy_n
                ! save old data
                hvy_block( :, :, :, 2:N_dF+1, hvy_active(k) ) = hvy_work( :, :, :, 1:N_dF, hvy_active(k) ) + (0.5_rk * dt) * hvy_work( :, :, :, 2*N_dF+1:3*N_dF, hvy_active(k) )
            end do

        case('2D_advection')
            do dF = 2, N_dF+1
                do k = 1, hvy_n
                    ! save old data
                    hvy_block( :, :, :, dF, hvy_active(k) ) = hvy_work( :, :, :, (dF-2)*5+1, hvy_active(k) ) + (0.5_rk * dt) * hvy_work( :, :, :, (dF-2)*5+3, hvy_active(k) )
                end do
            end do

        case default
            write(*,'(80("_"))')
            write(*,*) "ERROR: physics type is unknown"
            write(*,*) params%physics_type
            stop

    end select

    ! time measurement without ghost nodes synchronization
    sub_t1   = MPI_Wtime()
    time_sum = time_sum + (sub_t1 - sub_t0)

    ! synchronize ghostnodes
    call synchronize_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

    ! restart time
    sub_t0 = MPI_Wtime()

    ! RHS, depends on physics
    select case(params%physics_type)
        case('2D_convection_diffusion')
            ! loop over all datafields
            do dF = 2, N_dF+1
                ! loop over all active heavy data blocks
                do k = 1, hvy_n

                    ! set k3 step
                    hvy_work( :, :, :, (dF-2)*5+4, hvy_active(k) ) = hvy_block( :, :, :, dF, hvy_active(k) )
                    ! RHS
                    call RHS_2D_convection_diffusion( hvy_work( :, :, 1, (dF-2)*5+4, hvy_active(k) ), &
                                                      abs(hvy_block( 1, 2, 1, 1, hvy_active(k) ) - hvy_block( 1, 1, 1, 1, hvy_active(k) )), &
                                                      abs(hvy_block( 2, 2, 1, 1, hvy_active(k) ) - hvy_block( 2, 1, 1, 1, hvy_active(k) )), &
                                                      g, Bs, &
                                                      params%physics%u0( (dF-2)*2 + 1 ), params%physics%u0( (dF-2)*2 + 2 ), params%physics%nu( (dF-1) ), &
                                                      params%order_discretization  )

                end do
            end do

        case('2D_navier_stokes')
            ! loop over all active heavy data blocks
            do k = 1, hvy_n

                ! set k3 step
                hvy_work( :, :, :, 3*N_dF+1:4*N_dF, hvy_active(k) ) = hvy_block( :, :, :, 2:N_dF+1, hvy_active(k) )
                ! RHS
                call RHS_2D_navier_stokes( params%physics_ns, g, Bs, &
                                           abs(hvy_block( 1, 2, 1, 1, hvy_active(k) ) - hvy_block( 1, 1, 1, 1, hvy_active(k) )), &
                                           abs(hvy_block( 2, 2, 1, 1, hvy_active(k)) - hvy_block( 2, 1, 1, 1, hvy_active(k) )), &
                                           N_dF, &
                                           hvy_work( :, :, 1, 3*N_dF+1:4*N_dF, hvy_active(k) ) )

            end do

        case('3D_convection_diffusion')
            ! loop over all datafields
            do dF = 2, N_dF+1
                ! loop over all active heavy data blocks
                do k = 1, hvy_n

                    ! set k3 step
                    hvy_work( :, :, :, (dF-2)*5+4, hvy_active(k) ) = hvy_block( :, :, :, dF, hvy_active(k) )
                    ! RHS
                    call RHS_3D_convection_diffusion( hvy_work( :, :, :, (dF-2)*5+4, hvy_active(k) ), &
                                                      abs(hvy_block( 1, 2, 1, 1, hvy_active(k) ) - hvy_block( 1, 1, 1, 1, hvy_active(k) )), &
                                                      abs(hvy_block( 2, 2, 1, 1, hvy_active(k) ) - hvy_block( 2, 1, 1, 1, hvy_active(k) )), &
                                                      abs(hvy_block( 3, 2, 1, 1, hvy_active(k) ) - hvy_block( 3, 1, 1, 1, hvy_active(k) )), &
                                                      g, Bs, &
                                                      params%physics%u0( (dF-2)*2 + 1 ), params%physics%u0( (dF-2)*2 + 2 ), params%physics%u0( (dF-2)*2 + 3 ), params%physics%nu( (dF-1) ), &
                                                      params%order_discretization  )

                end do
            end do

        case('3D_navier_stokes')
            ! loop over all active heavy data blocks
            do k = 1, hvy_n

                ! set k3 step
                hvy_work( :, :, :, 3*N_dF+1:4*N_dF, hvy_active(k) ) = hvy_block( :, :, :, 2:N_dF+1, hvy_active(k) )
                ! RHS
                call RHS_3D_navier_stokes( params%physics_ns, g, Bs, &
                                           abs(hvy_block( 1, 2, 1, 1, hvy_active(k) ) - hvy_block( 1, 1, 1, 1, hvy_active(k) )), &
                                           abs(hvy_block( 2, 2, 1, 1, hvy_active(k)) - hvy_block( 2, 1, 1, 1, hvy_active(k) )), &
                                           abs(hvy_block( 3, 2, 1, 1, hvy_active(k)) - hvy_block( 3, 1, 1, 1, hvy_active(k) )), &
                                           N_dF, &
                                           hvy_work( :, :, :, 3*N_dF+1:4*N_dF, hvy_active(k) ) )

            end do

        case('2D_advection')
            ! loop over all datafields
            do dF = 2, N_dF+1
                ! loop over all active heavy data blocks
                do k = 1, hvy_n

                    ! light id of this block
                    call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )
                    ! compute blocks' spacing from treecode
                    call get_block_spacing_origin( params, lgt_id, lgt_block, xx0, ddx )

                    ! set k3 step
                    hvy_work( :, :, :, (dF-2)*5+4, hvy_active(k) ) = hvy_block( :, :, :, dF, hvy_active(k) )
                    ! RHS
!                    call RHS_2D_advection( hvy_work( :, :, 1, (dF-2)*5+4, hvy_active(k) ), &
!                                                      xx0(1:2), &
!                                                      ddx(1:2), &
!                                                      g, Bs, &
!                                                      time, &
!                                                      params%order_discretization  )
                    call RHS_2D_advection( hvy_work( :, :, 1, (dF-2)*5+4, hvy_active(k) ), &
                                                      xx0(1:2), &
                                                      ddx(1:2), &
                                                      g, Bs, &
                                                      time+dt/2.0_rk, &
                                                      params%order_discretization  )

                end do
            end do

        case default
            write(*,'(80("_"))')
            write(*,*) "ERROR: physics type is unknown"
            write(*,*) params%physics_type
            stop

    end select

    !***************************************************************************
    ! fourth stage
    !***************************************************************************
    ! loop over all datafields
    select case(params%physics_type)
        case('2D_convection_diffusion')
            do dF = 2, N_dF+1
                do k = 1, hvy_n
                    ! save old data
                    hvy_block( :, :, :, dF, hvy_active(k) ) = hvy_work( :, :, :, (dF-2)*5+1, hvy_active(k) ) + dt * hvy_work( :, :, :, (dF-2)*5+4, hvy_active(k) )
                end do
            end do

        case('2D_navier_stokes')
            do k = 1, hvy_n
                ! save old data
                hvy_block( :, :, :, 2:N_dF+1, hvy_active(k) ) = hvy_work( :, :, :, 1:N_dF, hvy_active(k) ) + dt * hvy_work( :, :, :, 3*N_dF+1:4*N_dF, hvy_active(k) )
            end do

        case('3D_convection_diffusion')
            do dF = 2, N_dF+1
                do k = 1, hvy_n
                    ! save old data
                    hvy_block( :, :, :, dF, hvy_active(k) ) = hvy_work( :, :, :, (dF-2)*5+1, hvy_active(k) ) + dt * hvy_work( :, :, :, (dF-2)*5+4, hvy_active(k) )
                end do
            end do

        case('3D_navier_stokes')
            do k = 1, hvy_n
                ! save old data
                hvy_block( :, :, :, 2:N_dF+1, hvy_active(k) ) = hvy_work( :, :, :, 1:N_dF, hvy_active(k) ) + dt * hvy_work( :, :, :, 3*N_dF+1:4*N_dF, hvy_active(k) )
            end do

        case('2D_advection')
            do dF = 2, N_dF+1
                do k = 1, hvy_n
                    ! save old data
                    hvy_block( :, :, :, dF, hvy_active(k) ) = hvy_work( :, :, :, (dF-2)*5+1, hvy_active(k) ) + dt * hvy_work( :, :, :, (dF-2)*5+4, hvy_active(k) )
                end do
            end do

        case default
            write(*,'(80("_"))')
            write(*,*) "ERROR: physics type is unknown"
            write(*,*) params%physics_type
            stop

    end select

    ! time measurement without ghost nodes synchronization
    sub_t1   = MPI_Wtime()
    time_sum = time_sum + (sub_t1 - sub_t0)

    ! synchronize ghostnodes
    call synchronize_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

    ! restart time
    sub_t0 = MPI_Wtime()

    ! RHS, depends on physics
    select case(params%physics_type)
        case('2D_convection_diffusion')
            ! loop over all datafields
            do dF = 2, N_dF+1
                ! loop over all active heavy data blocks
                do k = 1, hvy_n

                    ! set k4 step
                    hvy_work( :, :, :, (dF-2)*5+5, hvy_active(k) ) = hvy_block( :, :, :, dF, hvy_active(k) )
                    ! RHS
                    call RHS_2D_convection_diffusion( hvy_work( :, :, 1, (dF-2)*5+5, hvy_active(k) ), &
                                                      abs(hvy_block( 1, 2, 1, 1, hvy_active(k) ) - hvy_block( 1, 1, 1, 1, hvy_active(k) )), &
                                                      abs(hvy_block( 2, 2, 1, 1, hvy_active(k) ) - hvy_block( 2, 1, 1, 1, hvy_active(k) )), &
                                                      g, Bs, &
                                                      params%physics%u0( (dF-2)*2 + 1 ), params%physics%u0( (dF-2)*2 + 2 ), params%physics%nu( (dF-1) ), &
                                                      params%order_discretization  )

                end do
            end do

        case('2D_navier_stokes')
            ! loop over all active heavy data blocks
            do k = 1, hvy_n

                ! set k4 step
                hvy_work( :, :, :, 4*N_dF+1:5*N_dF, hvy_active(k) ) = hvy_block( :, :, :, 2:N_dF+1, hvy_active(k) )
                ! RHS
                call RHS_2D_navier_stokes( params%physics_ns, g, Bs, &
                                           abs(hvy_block( 1, 2, 1, 1, hvy_active(k) ) - hvy_block( 1, 1, 1, 1, hvy_active(k) )), &
                                           abs(hvy_block( 2, 2, 1, 1, hvy_active(k)) - hvy_block( 2, 1, 1, 1, hvy_active(k) )), &
                                           N_dF, &
                                           hvy_work( :, :, 1, 4*N_dF+1:5*N_dF, hvy_active(k) ) )

            end do

        case('3D_convection_diffusion')
            ! loop over all datafields
            do dF = 2, N_dF+1
                ! loop over all active heavy data blocks
                do k = 1, hvy_n

                    ! set k4 step
                    hvy_work( :, :, :, (dF-2)*5+5, hvy_active(k) ) = hvy_block( :, :, :, dF, hvy_active(k) )
                    ! RHS
                    call RHS_3D_convection_diffusion( hvy_work( :, :, :, (dF-2)*5+5, hvy_active(k) ), &
                                                      abs(hvy_block( 1, 2, 1, 1, hvy_active(k) ) - hvy_block( 1, 1, 1, 1, hvy_active(k) )), &
                                                      abs(hvy_block( 2, 2, 1, 1, hvy_active(k) ) - hvy_block( 2, 1, 1, 1, hvy_active(k) )), &
                                                      abs(hvy_block( 3, 2, 1, 1, hvy_active(k) ) - hvy_block( 3, 1, 1, 1, hvy_active(k) )), &
                                                      g, Bs, &
                                                      params%physics%u0( (dF-2)*2 + 1 ), params%physics%u0( (dF-2)*2 + 2 ), params%physics%u0( (dF-2)*2 + 3 ), params%physics%nu( (dF-1) ), &
                                                      params%order_discretization  )

                end do
            end do

        case('3D_navier_stokes')
            ! loop over all active heavy data blocks
            do k = 1, hvy_n

                ! set k4 step
                hvy_work( :, :, :, 4*N_dF+1:5*N_dF, hvy_active(k) ) = hvy_block( :, :, :, 2:N_dF+1, hvy_active(k) )
                ! RHS
                call RHS_3D_navier_stokes( params%physics_ns, g, Bs, &
                                           abs(hvy_block( 1, 2, 1, 1, hvy_active(k) ) - hvy_block( 1, 1, 1, 1, hvy_active(k) )), &
                                           abs(hvy_block( 2, 2, 1, 1, hvy_active(k)) - hvy_block( 2, 1, 1, 1, hvy_active(k) )), &
                                           abs(hvy_block( 3, 2, 1, 1, hvy_active(k)) - hvy_block( 3, 1, 1, 1, hvy_active(k) )), &
                                           N_dF, &
                                           hvy_work( :, :, :, 4*N_dF+1:5*N_dF, hvy_active(k) ) )

            end do

        case('2D_advection')
            ! loop over all datafields
            do dF = 2, N_dF+1
                ! loop over all active heavy data blocks
                do k = 1, hvy_n

                    ! light id of this block
                    call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )
                    ! compute blocks' spacing from treecode
                    call get_block_spacing_origin( params, lgt_id, lgt_block, xx0, ddx )

                    ! set k4 step
                    hvy_work( :, :, :, (dF-2)*5+5, hvy_active(k) ) = hvy_block( :, :, :, dF, hvy_active(k) )
                    ! RHS
!                    call RHS_2D_advection( hvy_work( :, :, 1, (dF-2)*5+5, hvy_active(k) ), &
!                                                      xx0(1:2), &
!                                                      ddx(1:2), &
!                                                      g, Bs, &
!                                                      time, &
!                                                      params%order_discretization  )
                    call RHS_2D_advection( hvy_work( :, :, 1, (dF-2)*5+5, hvy_active(k) ), &
                                                      xx0(1:2), &
                                                      ddx(1:2), &
                                                      g, Bs, &
                                                      time+dt, &
                                                      params%order_discretization  )

                end do
            end do

        case default
            write(*,'(80("_"))')
            write(*,*) "ERROR: physics type is unknown"
            write(*,*) params%physics_type
            stop

    end select

    !***************************************************************************
    ! final stage
    !***************************************************************************
    select case(params%physics_type)
        case('2D_convection_diffusion')
            ! loop over all datafields
            do dF = 2, N_dF+1
                ! loop over all active heavy data blocks
                do k = 1, hvy_n

                    ! final step
                    hvy_block( :, :, :, dF, hvy_active(k) ) = hvy_work( :, :, :, (dF-2)*5+1, hvy_active(k) ) &
                                                            + (dt/6.0_rk) * ( hvy_work( :, :, :, (dF-2)*5+2, hvy_active(k) ) &
                                                            + 2.0_rk * hvy_work( :, :, :, (dF-2)*5+3, hvy_active(k) ) &
                                                            + 2.0_rk * hvy_work( :, :, :, (dF-2)*5+4, hvy_active(k) ) &
                                                            + hvy_work( :, :, :, (dF-2)*5+5, hvy_active(k) ) )

                end do
            end do

        case('2D_navier_stokes')
            ! loop over all active heavy data blocks
            do k = 1, hvy_n

                ! final step
                hvy_block( :, :, :, 2:N_dF+1, hvy_active(k) )   = hvy_work( :, :, :, 1:N_dF, hvy_active(k) ) &
                                                                + (dt/6.0_rk) * ( hvy_work( :, :, :, N_dF+1:2*N_dF, hvy_active(k) ) &
                                                                + 2.0_rk * hvy_work( :, :, :, 2*N_dF+1:3*N_dF, hvy_active(k) ) &
                                                                + 2.0_rk * hvy_work( :, :, :, 3*N_dF+1:4*N_dF, hvy_active(k) ) &
                                                                + hvy_work( :, :, :, 4*N_dF+1:5*N_dF, hvy_active(k) ) )

            end do

        case('3D_convection_diffusion')
            ! loop over all datafields
            do dF = 2, N_dF+1
                ! loop over all active heavy data blocks
                do k = 1, hvy_n

                    ! final step
                    hvy_block( :, :, :, dF, hvy_active(k) ) = hvy_work( :, :, :, (dF-2)*5+1, hvy_active(k) ) &
                                                            + (dt/6.0_rk) * ( hvy_work( :, :, :, (dF-2)*5+2, hvy_active(k) ) &
                                                            + 2.0_rk * hvy_work( :, :, :, (dF-2)*5+3, hvy_active(k) ) &
                                                            + 2.0_rk * hvy_work( :, :, :, (dF-2)*5+4, hvy_active(k) ) &
                                                            + hvy_work( :, :, :, (dF-2)*5+5, hvy_active(k) ) )

                end do
            end do

        case('3D_navier_stokes')
            ! loop over all active heavy data blocks
            do k = 1, hvy_n

                ! final step
                hvy_block( :, :, :, 2:N_dF+1, hvy_active(k) )   = hvy_work( :, :, :, 1:N_dF, hvy_active(k) ) &
                                                                + (dt/6.0_rk) * ( hvy_work( :, :, :, N_dF+1:2*N_dF, hvy_active(k) ) &
                                                                + 2.0_rk * hvy_work( :, :, :, 2*N_dF+1:3*N_dF, hvy_active(k) ) &
                                                                + 2.0_rk * hvy_work( :, :, :, 3*N_dF+1:4*N_dF, hvy_active(k) ) &
                                                                + hvy_work( :, :, :, 4*N_dF+1:5*N_dF, hvy_active(k) ) )

            end do

        case('2D_advection')
            ! loop over all datafields
            do dF = 2, N_dF+1
                ! loop over all active heavy data blocks
                do k = 1, hvy_n

                    ! final step
                    hvy_block( :, :, :, dF, hvy_active(k) ) = hvy_work( :, :, :, (dF-2)*5+1, hvy_active(k) ) &
                                                            + (dt/6.0_rk) * ( hvy_work( :, :, :, (dF-2)*5+2, hvy_active(k) ) &
                                                            + 2.0_rk * hvy_work( :, :, :, (dF-2)*5+3, hvy_active(k) ) &
                                                            + 2.0_rk * hvy_work( :, :, :, (dF-2)*5+4, hvy_active(k) ) &
                                                            + hvy_work( :, :, :, (dF-2)*5+5, hvy_active(k) ) )

                end do
            end do

        case default
            write(*,'(80("_"))')
            write(*,*) "ERROR: physics type is unknown"
            write(*,*) params%physics_type
            stop

    end select

    ! increase time variable after all RHS substeps
    time = time_dt

    ! end time
    sub_t1   = MPI_Wtime()
    time_sum = time_sum + (sub_t1 - sub_t0)
    ! write time
    if ( params%debug ) then
        ! find free or corresponding line
        k = 1
        do while ( debug%name_comp_time(k) /= "---" )
            ! entry for current subroutine exists
            if ( debug%name_comp_time(k) == "time_step (w/o ghost synch.)" ) exit
            k = k + 1
        end do
        ! write time
        debug%name_comp_time(k) = "time_step (w/o ghost synch.)"
        debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1
        debug%comp_time(k, 2)   = debug%comp_time(k, 2) + time_sum
    end if

end subroutine time_step_RK4
