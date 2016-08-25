! ********************************
! 2D AMR prototype
! --------------------------------
! 
! time step main function, RK4
!
! name: time_step.f90
! date: 02.08.2016
! author: msr
! version: 0.1
! 
! ********************************

subroutine time_step(time)

    use module_params
    use module_blocks

    implicit none

    real(kind=rk), intent(inout)  			    :: time

    integer(kind=ik)				            :: g, N, block_N, k, block_num
    real(kind=rk)				                :: dt

    g 					    = blocks_params%number_ghost_nodes
    N 					    = blocks_params%size_block

    block_N                 = size(blocks_params%active_list, dim=1)

    call calc_dt(dt)
    time 				    = time + dt
    ! last timestep fits in maximal time
    if (time > params%time_max) time = params%time_max

    ! output on screen
    write(*,'(a,i5)') "number of active blocks = ", block_N

    !------------------------------
    ! first stage
    ! synchronize ghostnodes
    call synchronize_ghosts()

    do k = 1, block_N

        block_num                   = blocks_params%active_list(k)
        blocks(block_num)%data_old  = blocks(block_num)%data2
        blocks(block_num)%k1        = blocks(block_num)%data_old
        ! RHS
        call RHS_2D_block(blocks(block_num)%k1, blocks(block_num)%dx, blocks(block_num)%dy, g, N)

    end do

    !------------------------------
    ! second stage
    do k = 1, block_N

        block_num                   = blocks_params%active_list(k)
        blocks(block_num)%data2     = blocks(block_num)%data_old + 0.5_rk * dt * blocks(block_num)%k1
        blocks(block_num)%data1     = blocks(block_num)%data2(g+1:N+g, g+1:N+g)

    end do
    ! synchronize ghostnodes
    call synchronize_ghosts()
    do k = 1, block_N

        block_num                   = blocks_params%active_list(k)
        blocks(block_num)%k2        = blocks(block_num)%data2
        ! RHS
        call RHS_2D_block(blocks(block_num)%k2, blocks(block_num)%dx, blocks(block_num)%dy, g, N)

    end do

    !------------------------------
    ! third stage
    do k = 1, block_N

        block_num                   = blocks_params%active_list(k)
        blocks(block_num)%data2     = blocks(block_num)%data_old + 0.5_rk * dt * blocks(block_num)%k2
        blocks(block_num)%data1     = blocks(block_num)%data2(g+1:N+g, g+1:N+g)

    end do
    ! synchronize ghostnodes
    call synchronize_ghosts()
    do k = 1, block_N

        block_num                   = blocks_params%active_list(k)
        blocks(block_num)%k3        = blocks(block_num)%data2
        ! RHS
        call RHS_2D_block(blocks(block_num)%k3, blocks(block_num)%dx, blocks(block_num)%dy, g, N)

    end do

    !------------------------------
    ! fourth stage
    do k = 1, block_N

        block_num                   = blocks_params%active_list(k)
        blocks(block_num)%data2     = blocks(block_num)%data_old + dt * blocks(block_num)%k3
        blocks(block_num)%data1     = blocks(block_num)%data2(g+1:N+g, g+1:N+g)

    end do
    ! synchronize ghostnodes
    call synchronize_ghosts()
    do k = 1, block_N

        block_num                   = blocks_params%active_list(k)
        blocks(block_num)%k4        = blocks(block_num)%data2
        ! RHS
        call RHS_2D_block(blocks(block_num)%k4, blocks(block_num)%dx, blocks(block_num)%dy, g, N)

    end do

    !------------------------------
    ! final stage
    do k = 1, block_N

        block_num                   = blocks_params%active_list(k)

        blocks(block_num)%data2     = blocks(block_num)%data_old + dt/6.0_rk * ( blocks(block_num)%k1 &
                                                                        + 2.0_rk*blocks(block_num)%k2 &
                                                                        + 2.0_rk*blocks(block_num)%k3 &
                                                                        + blocks(block_num)%k4 )
        blocks(block_num)%data1     = blocks(block_num)%data2(g+1:N+g, g+1:N+g)

    end do

end subroutine time_step
