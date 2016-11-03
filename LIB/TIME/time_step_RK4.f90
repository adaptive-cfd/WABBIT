! ********************************
! WABBIT
! --------------------------------
!
! time step main function, RK4
!
! name: time_step_RK4.f90
! date: 26.10.2016
! author: msr
! version: 0.3
!
! ********************************

subroutine time_step_RK4(time)

    use module_params
    use module_blocks

    implicit none

    real(kind=rk), intent(inout)  			    :: time

    integer(kind=ik)				            :: dF, g, Bs, N, k
    real(kind=rk)				                :: dt

    g 					    = blocks_params%number_ghost_nodes
    Bs 					    = blocks_params%size_block

    N                       = blocks_params%number_max_blocks_data

    call calc_dt(dt)

    time 				    = time + dt
    ! last timestep should fit in maximal time
    if (time >= params%time_max) then
        time = time - dt
        dt = params%time_max - time
        time = params%time_max
    end if

    !------------------------------
    ! first stage
    ! synchronize ghostnodes
    call synchronize_ghosts()

    ! loop over all internal blocks
    do k = 1, N
        ! if block is active
        if ( blocks_data(k)%block_id /= -1 ) then
            ! loop over all data fields
            do dF = 1, blocks_params%number_data_fields
                blocks_data(k)%data_fields(dF)%data_old  = blocks_data(k)%data_fields(dF)%data_
                blocks_data(k)%data_fields(dF)%k1        = blocks_data(k)%data_fields(dF)%data_old
                ! RHS
                call RHS_2D_convection_diffusion(blocks_data(k)%data_fields(dF)%k1(:,:), blocks_data(k)%dx, blocks_data(k)%dy, g, Bs, params%u0(1, dF), params%u0(2, dF), params%nu(dF))
            end do
        end if
    end do

    !------------------------------
    ! second stage
    do k = 1, N
        if ( blocks_data(k)%block_id /= -1 ) then
            ! loop over all data fields
            do dF = 1, blocks_params%number_data_fields
                blocks_data(k)%data_fields(dF)%data_     = blocks_data(k)%data_fields(dF)%data_old + (0.5_rk * dt) * ( blocks_data(k)%data_fields(dF)%k1 )
            end do
        end if
    end do
    ! synchronize ghostnodes
    call synchronize_ghosts()

    do k = 1, N
        if ( blocks_data(k)%block_id /= -1 ) then
            ! loop over all data fields
            do dF = 1, blocks_params%number_data_fields
                blocks_data(k)%data_fields(dF)%k2        = blocks_data(k)%data_fields(dF)%data_
                ! RHS
                call RHS_2D_convection_diffusion(blocks_data(k)%data_fields(dF)%k2(:,:), blocks_data(k)%dx, blocks_data(k)%dy, g, Bs, params%u0(1, dF), params%u0(2, dF), params%nu(dF))
            end do
        end if
    end do

    !------------------------------
    ! third stage
    do k = 1, N
        if ( blocks_data(k)%block_id /= -1 ) then
            ! loop over all data fields
            do dF = 1, blocks_params%number_data_fields
                blocks_data(k)%data_fields(dF)%data_     = blocks_data(k)%data_fields(dF)%data_old + (0.5_rk * dt) * ( blocks_data(k)%data_fields(dF)%k2 )
            end do
        end if
    end do
    ! synchronize ghostnodes
    call synchronize_ghosts()

    do k = 1, N
        if ( blocks_data(k)%block_id /= -1 ) then
            ! loop over all data fields
            do dF = 1, blocks_params%number_data_fields
                blocks_data(k)%data_fields(dF)%k3        = blocks_data(k)%data_fields(dF)%data_
                ! RHS
                call RHS_2D_convection_diffusion(blocks_data(k)%data_fields(dF)%k3(:,:), blocks_data(k)%dx, blocks_data(k)%dy, g, Bs, params%u0(1, dF), params%u0(2, dF), params%nu(dF))
            end do
        end if
    end do

    !------------------------------
    ! fourth stage
    do k = 1, N
        if ( blocks_data(k)%block_id /= -1 ) then
            ! loop over all data fields
            do dF = 1, blocks_params%number_data_fields
                blocks_data(k)%data_fields(dF)%data_     = blocks_data(k)%data_fields(dF)%data_old + dt * ( blocks_data(k)%data_fields(dF)%k3 )
            end do
        end if
    end do
    ! synchronize ghostnodes
    call synchronize_ghosts()

    do k = 1, N
        if ( blocks_data(k)%block_id /= -1 ) then
            ! loop over all data fields
            do dF = 1, blocks_params%number_data_fields
                blocks_data(k)%data_fields(dF)%k4        = blocks_data(k)%data_fields(dF)%data_
                ! RHS
                call RHS_2D_convection_diffusion(blocks_data(k)%data_fields(dF)%k4(:,:), blocks_data(k)%dx, blocks_data(k)%dy, g, Bs, params%u0(1, dF), params%u0(2, dF), params%nu(dF))
            end do
        end if
    end do

    !------------------------------
    ! final stage
    do k = 1, N
        if ( blocks_data(k)%block_id /= -1 ) then
            ! loop over all data fields
            do dF = 1, blocks_params%number_data_fields
                blocks_data(k)%data_fields(dF)%data_     = blocks_data(k)%data_fields(dF)%data_old + (dt/6.0_rk) * ( blocks_data(k)%data_fields(dF)%k1 &
                                                        + 2.0_rk*blocks_data(k)%data_fields(dF)%k2 &
                                                        + 2.0_rk*blocks_data(k)%data_fields(dF)%k3 &
                                                        + blocks_data(k)%data_fields(dF)%k4 )
            end do
        end if
    end do

end subroutine time_step_RK4
