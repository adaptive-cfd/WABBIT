! ********************************
! WABBIT
! --------------------------------
!
! time step main function, RK4
!
! name: time_step_RK4.f90
! date: 30.09.2016
! author: msr
! version: 0.2
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

    N                       = blocks_params%number_max_blocks

    call calc_dt(dt)

    time 				    = time + dt
    ! last timestep should fit in maximal time
    if (time >= params%time_max) then
        time = time - dt
        dt = params%time_max - time
        time = params%time_max
    end if

    ! check number of data fields
    if ( blocks_params%number_data_fields == 1 ) then
        ! single data field
        dF = blocks_params%number_data_fields

        !------------------------------
        ! first stage
        ! synchronize ghostnodes
        call synchronize_ghosts()

        do k = 1, N
            if (blocks(k)%active) then
                blocks(k)%data_fields(dF)%data_old  = blocks(k)%data_fields(dF)%data_
                blocks(k)%data_fields(dF)%k1        = blocks(k)%data_fields(dF)%data_old
                ! RHS
                call RHS_2D_convection_diffusion(blocks(k)%data_fields(dF)%k1(:,:), blocks(k)%dx, blocks(k)%dy, g, Bs)
            end if
        end do

        !------------------------------
        ! second stage
        do k = 1, N
            if (blocks(k)%active) then
                blocks(k)%data_fields(dF)%data_     = blocks(k)%data_fields(dF)%data_old + (0.5_rk * dt) * ( blocks(k)%data_fields(dF)%k1 )
            end if
        end do
        ! synchronize ghostnodes
        call synchronize_ghosts()

        do k = 1, N
            if (blocks(k)%active) then
                blocks(k)%data_fields(dF)%k2        = blocks(k)%data_fields(dF)%data_
                ! RHS
                call RHS_2D_convection_diffusion(blocks(k)%data_fields(dF)%k2(:,:), blocks(k)%dx, blocks(k)%dy, g, Bs)
            end if
        end do

        !------------------------------
        ! third stage
        do k = 1, N
            if (blocks(k)%active) then
                blocks(k)%data_fields(dF)%data_     = blocks(k)%data_fields(dF)%data_old + (0.5_rk * dt) * ( blocks(k)%data_fields(dF)%k2 )
            end if
        end do
        ! synchronize ghostnodes
        call synchronize_ghosts()

        do k = 1, N
            if (blocks(k)%active) then
                blocks(k)%data_fields(dF)%k3        = blocks(k)%data_fields(dF)%data_
                ! RHS
                call RHS_2D_convection_diffusion(blocks(k)%data_fields(dF)%k3(:,:), blocks(k)%dx, blocks(k)%dy, g, Bs)
            end if
        end do

        !------------------------------
        ! fourth stage
        do k = 1, N
            if (blocks(k)%active) then
                blocks(k)%data_fields(dF)%data_     = blocks(k)%data_fields(dF)%data_old + dt * ( blocks(k)%data_fields(dF)%k3 )
            end if
        end do
        ! synchronize ghostnodes
        call synchronize_ghosts()

        do k = 1, N
            if (blocks(k)%active) then
                blocks(k)%data_fields(dF)%k4        = blocks(k)%data_fields(dF)%data_
                ! RHS
                call RHS_2D_convection_diffusion(blocks(k)%data_fields(dF)%k4(:,:), blocks(k)%dx, blocks(k)%dy, g, Bs)
            end if
        end do

        !------------------------------
        ! final stage
        do k = 1, N
            if (blocks(k)%active) then
                blocks(k)%data_fields(dF)%data_     = blocks(k)%data_fields(dF)%data_old + (dt/6.0_rk) * ( blocks(k)%data_fields(dF)%k1 &
                                                    + 2.0_rk*blocks(k)%data_fields(dF)%k2 &
                                                    + 2.0_rk*blocks(k)%data_fields(dF)%k3 &
                                                    + blocks(k)%data_fields(dF)%k4 )
            end if
        end do

    else
        ! more than one data field
        ! to do

    end if

end subroutine time_step_RK4
