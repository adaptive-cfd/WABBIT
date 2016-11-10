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
    integer(kind=ik)                    :: k, N, dF

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
            my_dx = min(my_dx, block_data(k, 1, 2, 1 ) - block_data(k, 1, 1, 1 ) )
        end if
    end do

    ! synchronize dx
    call MPI_Allreduce(my_dx, dx, 1, MPI_REAL8, MPI_MIN, MPI_COMM_WORLD, ierr)

    ! calculate time step, loop over all data fields
    do dF = 2, params%number_data_fields+1
        dt = min(dt, params%CFL * dx / norm2( params%u0( (dF-2)*2 + 1 : (dF-2)*2 + 2 ) ) )
    end do

    time = time + dt
    ! last timestep should fit in maximal time
    if (time >= params%time_max) then
        time = time - dt
        dt = params%time_max - time
        time = params%time_max
    end if

    !------------------------------
    ! first stage
    ! synchronize ghostnodes
    call synchronize_ghosts( params, block_list, block_data, neighbor_list )

!    ! loop over all internal blocks
!    do k = 1, N
!        ! if block is active
!        if ( blocks_data(k)%block_id /= -1 ) then
!            ! loop over all data fields
!            do dF = 1, blocks_params%number_data_fields
!                blocks_data(k)%data_fields(dF)%data_old  = blocks_data(k)%data_fields(dF)%data_
!                blocks_data(k)%data_fields(dF)%k1        = blocks_data(k)%data_fields(dF)%data_old
!                ! RHS
!                call RHS_2D_convection_diffusion(blocks_data(k)%data_fields(dF)%k1(:,:), blocks_data(k)%dx, blocks_data(k)%dy, g, Bs, params%u0(1, dF), params%u0(2, dF), params%nu(dF))
!            end do
!        end if
!    end do

!    !------------------------------
!    ! second stage
!    do k = 1, N
!        if ( blocks_data(k)%block_id /= -1 ) then
!            ! loop over all data fields
!            do dF = 1, blocks_params%number_data_fields
!                blocks_data(k)%data_fields(dF)%data_     = blocks_data(k)%data_fields(dF)%data_old + (0.5_rk * dt) * ( blocks_data(k)%data_fields(dF)%k1 )
!            end do
!        end if
!    end do
!    ! synchronize ghostnodes
!    call synchronize_ghosts()
!
!    do k = 1, N
!        if ( blocks_data(k)%block_id /= -1 ) then
!            ! loop over all data fields
!            do dF = 1, blocks_params%number_data_fields
!                blocks_data(k)%data_fields(dF)%k2        = blocks_data(k)%data_fields(dF)%data_
!                ! RHS
!                call RHS_2D_convection_diffusion(blocks_data(k)%data_fields(dF)%k2(:,:), blocks_data(k)%dx, blocks_data(k)%dy, g, Bs, params%u0(1, dF), params%u0(2, dF), params%nu(dF))
!            end do
!        end if
!    end do
!
!    !------------------------------
!    ! third stage
!    do k = 1, N
!        if ( blocks_data(k)%block_id /= -1 ) then
!            ! loop over all data fields
!            do dF = 1, blocks_params%number_data_fields
!                blocks_data(k)%data_fields(dF)%data_     = blocks_data(k)%data_fields(dF)%data_old + (0.5_rk * dt) * ( blocks_data(k)%data_fields(dF)%k2 )
!            end do
!        end if
!    end do
!    ! synchronize ghostnodes
!    call synchronize_ghosts()
!
!    do k = 1, N
!        if ( blocks_data(k)%block_id /= -1 ) then
!            ! loop over all data fields
!            do dF = 1, blocks_params%number_data_fields
!                blocks_data(k)%data_fields(dF)%k3        = blocks_data(k)%data_fields(dF)%data_
!                ! RHS
!                call RHS_2D_convection_diffusion(blocks_data(k)%data_fields(dF)%k3(:,:), blocks_data(k)%dx, blocks_data(k)%dy, g, Bs, params%u0(1, dF), params%u0(2, dF), params%nu(dF))
!            end do
!        end if
!    end do
!
!    !------------------------------
!    ! fourth stage
!    do k = 1, N
!        if ( blocks_data(k)%block_id /= -1 ) then
!            ! loop over all data fields
!            do dF = 1, blocks_params%number_data_fields
!                blocks_data(k)%data_fields(dF)%data_     = blocks_data(k)%data_fields(dF)%data_old + dt * ( blocks_data(k)%data_fields(dF)%k3 )
!            end do
!        end if
!    end do
!    ! synchronize ghostnodes
!    call synchronize_ghosts()
!
!    do k = 1, N
!        if ( blocks_data(k)%block_id /= -1 ) then
!            ! loop over all data fields
!            do dF = 1, blocks_params%number_data_fields
!                blocks_data(k)%data_fields(dF)%k4        = blocks_data(k)%data_fields(dF)%data_
!                ! RHS
!                call RHS_2D_convection_diffusion(blocks_data(k)%data_fields(dF)%k4(:,:), blocks_data(k)%dx, blocks_data(k)%dy, g, Bs, params%u0(1, dF), params%u0(2, dF), params%nu(dF))
!            end do
!        end if
!    end do
!
!    !------------------------------
!    ! final stage
!    do k = 1, N
!        if ( blocks_data(k)%block_id /= -1 ) then
!            ! loop over all data fields
!            do dF = 1, blocks_params%number_data_fields
!                blocks_data(k)%data_fields(dF)%data_     = blocks_data(k)%data_fields(dF)%data_old + (dt/6.0_rk) * ( blocks_data(k)%data_fields(dF)%k1 &
!                                                        + 2.0_rk*blocks_data(k)%data_fields(dF)%k2 &
!                                                        + 2.0_rk*blocks_data(k)%data_fields(dF)%k3 &
!                                                        + blocks_data(k)%data_fields(dF)%k4 )
!            end do
!        end if
!    end do

end subroutine time_step_RK4
