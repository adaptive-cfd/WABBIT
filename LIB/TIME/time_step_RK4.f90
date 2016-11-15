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

subroutine time_step_RK4( time, params, block_list, block_data, neighbor_list )

!---------------------------------------------------------------------------------------------
! modules

    use mpi
    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! time varible
    real(kind=rk), intent(inout)        :: time

    ! user defined parameter structure
    type (type_params), intent(in)      :: params
    ! light data array
    integer(kind=ik), intent(in)        :: block_list(:, :)
    ! heavy data array - block data
    real(kind=rk), intent(inout)        :: block_data(:, :, :, :)
    ! neighbor list
    integer(kind=ik), intent(in)        :: neighbor_list(:)

    ! grid parameter
    integer(kind=ik)                    :: Bs, g
    ! loop variables
    integer(kind=ik)                    :: k, N, dF, N_dF

    ! time step, dx
    real(kind=rk)                       :: dt, dx, my_dx

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank

!---------------------------------------------------------------------------------------------
! interfaces

    interface
        subroutine synchronize_ghosts( params, block_list, block_data, neighbor_list )
            use module_params
            type (type_params), intent(in)              :: params
            integer(kind=ik), intent(in)                :: block_list(:, :)
            real(kind=rk), intent(inout)                :: block_data(:, :, :, :)
            integer(kind=ik), intent(in)                :: neighbor_list(:)
        end subroutine synchronize_ghosts

    end interface

!---------------------------------------------------------------------------------------------
! variables initialization

    N     = params%number_blocks
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

    ! determinate process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    ! ----------------------------------------------------------------------------------------
    ! calculate time step
    ! loop over all blocks (heavy data)
    do k = 1, N
        ! block is active
        if ( block_list(rank*N + k , 1) /= -1 ) then
            my_dx = min(my_dx, block_data(1, 2, 1, k ) - block_data(1, 1, 1, k ) )
        end if
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
        call synchronize_ghosts( params, block_list, block_data, neighbor_list )

!        ! loop over all heavy data blocks
!        do k = 1, N
!            ! block is active
!            if ( block_list( rank*N + k , 1) /= -1 ) then
!                ! save old data
!                block_data( :, :, N_dF+2, k ) = block_data( :, :, dF, k )
!                ! set k1 step
!                block_data( :, :, N_dF+3, k ) = block_data( :, :, dF, k )
!                ! RHS
!                call RHS_2D_convection_diffusion( block_data( :, :, N_dF+3, k ), &
!                                                  abs(block_data( 1, 2, 1, k ) - block_data( 1, 1, 1, k )), &
!                                                  abs(block_data( 2, 2, 1, k ) - block_data( 2, 1, 1, k )), &
!                                                  g, Bs, &
!                                                  params%u0( (dF-2)*2 + 1 ), params%u0( (dF-2)*2 + 2 ), params%nu( (dF-1) ), &
!                                                  params%order_discretization  )
!            end if
!        end do
!
!        !------------------------------
!        ! second stage
!        do k = 1, N
!            ! block is active
!            if ( block_list( rank*N + k , 1) /= -1 ) then
!                ! save old data
!                block_data( :, :, dF, k ) = block_data( :, :, N_dF+2, k ) + (0.5_rk * dt) * block_data( :, :, N_dF+3, k )
!            end if
!        end do
!
!        ! synchronize ghostnodes
!        call synchronize_ghosts( params, block_list, block_data, neighbor_list )
!
!
!
!        do k = 1, N
!            ! block is active
!            if ( block_list( rank*N + k , 1) /= -1 ) then
!                ! set k2 step
!                block_data( :, :, N_dF+4, k ) = block_data( :, :, dF, k )
!                ! RHS
!                call RHS_2D_convection_diffusion( block_data( :, :, N_dF+4, k ), &
!                                                  abs(block_data( 1, 2, 1, k ) - block_data( 1, 1, 1, k )), &
!                                                  abs(block_data( 2, 2, 1, k ) - block_data( 2, 1, 1, k )), &
!                                                  g, Bs, &
!                                                  params%u0( (dF-2)*2 + 1 ), params%u0( (dF-2)*2 + 2 ), params%nu( (dF-1) ), &
!                                                  params%order_discretization  )
!            end if
!        end do
!
!        !------------------------------
!        ! third stage
!        do k = 1, N
!            ! block is active
!            if ( block_list( rank*N + k , 1) /= -1 ) then
!                ! save old data
!                block_data( :, :, dF, k ) = block_data( :, :, N_dF+2, k ) + (0.5_rk * dt) * block_data( :, :, N_dF+4, k )
!            end if
!        end do
!
!        ! synchronize ghostnodes
!        call synchronize_ghosts( params, block_list, block_data, neighbor_list )
!
!        do k = 1, N
!            ! block is active
!            if ( block_list( rank*N + k , 1) /= -1 ) then
!                ! set k3 step
!                block_data( :, :, N_dF+5, k ) = block_data( :, :, dF, k )
!                ! RHS
!                call RHS_2D_convection_diffusion( block_data( :, :, N_dF+5, k ), &
!                                                  abs(block_data( 1, 2, 1, k ) - block_data( 1, 1, 1, k )), &
!                                                  abs(block_data( 2, 2, 1, k ) - block_data( 2, 1, 1, k )), &
!                                                  g, Bs, &
!                                                  params%u0( (dF-2)*2 + 1 ), params%u0( (dF-2)*2 + 2 ), params%nu( (dF-1) ), &
!                                                  params%order_discretization  )
!            end if
!        end do
!
!        !------------------------------
!        ! fourth stage
!        do k = 1, N
!            ! block is active
!            if ( block_list( rank*N + k , 1) /= -1 ) then
!                ! save old data
!                block_data( :, :, dF, k ) = block_data( :, :, N_dF+2, k ) + (0.5_rk * dt) * block_data( :, :, N_dF+5, k )
!            end if
!        end do
!
!        ! synchronize ghostnodes
!        call synchronize_ghosts( params, block_list, block_data, neighbor_list )
!
!        do k = 1, N
!            ! block is active
!            if ( block_list( rank*N + k , 1) /= -1 ) then
!                ! set k4 step
!                block_data( :, :, N_dF+6, k ) = block_data( :, :, dF, k )
!                ! RHS
!                call RHS_2D_convection_diffusion( block_data( :, :, N_dF+6, k ), &
!                                                  abs(block_data( 1, 2, 1, k ) - block_data( 1, 1, 1, k )), &
!                                                  abs(block_data( 2, 2, 1, k ) - block_data( 2, 1, 1, k )), &
!                                                  g, Bs, &
!                                                  params%u0( (dF-2)*2 + 1 ), params%u0( (dF-2)*2 + 2 ), params%nu( (dF-1) ), &
!                                                  params%order_discretization  )
!            end if
!        end do
!
!        !------------------------------
!        ! final stage
!        do k = 1, N
!            ! block is active
!            if ( block_list( rank*N + k , 1) /= -1 ) then
!                ! final step
!                block_data( :, :, dF, k ) = block_data( :, :, N_dF+2, k ) &
!                                          + (dt/6.0_rk) * ( block_data( :, :, N_dF+3, k ) &
!                                          + 2.0_rk * block_data( :, :, N_dF+4, k ) &
!                                          + 2.0_rk * block_data( :, :, N_dF+5, k ) &
!                                          + block_data( :, :, N_dF+6, k ) )
!            end if
!        end do

    end do

end subroutine time_step_RK4
