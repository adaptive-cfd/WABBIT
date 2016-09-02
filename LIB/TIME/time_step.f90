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

    integer(kind=ik)				            :: dF, g, Bs, N, k, block_num
    real(kind=rk)				                :: dt

    g 					    = blocks_params%number_ghost_nodes
    Bs 					    = blocks_params%size_block

    N                       = size(blocks_params%active_list, dim=1)

    call calc_dt(dt)
    time 				    = time + dt
    ! last timestep fits in maximal time
    if (time > params%time_max) time = params%time_max

    ! check number of data fields
    if ( blocks_params%number_data_fields == 1 ) then
        ! single data field
        dF = blocks_params%number_data_fields

        !------------------------------
        ! first stage
        ! synchronize ghostnodes
        call synchronize_ghosts()

        do k = 1, N

            block_num                                   = blocks_params%active_list(k)
            blocks(block_num)%data_fields(dF)%data_old  = blocks(block_num)%data_fields(dF)%data_
            blocks(block_num)%data_fields(dF)%k1        = blocks(block_num)%data_fields(dF)%data_old
            ! RHS
            call RHS_2D_block(blocks(block_num)%data_fields(dF)%k1(:,:), blocks(block_num)%dx, blocks(block_num)%dy, g, Bs)

        end do

        !------------------------------
        ! second stage
        do k = 1, N

            block_num                                   = blocks_params%active_list(k)
            blocks(block_num)%data_fields(dF)%data_     = blocks(block_num)%data_fields(dF)%data_old + 0.5_rk * dt * blocks(block_num)%data_fields(dF)%k1

        end do
        ! synchronize ghostnodes
        call synchronize_ghosts()
        do k = 1, N

            block_num                                   = blocks_params%active_list(k)
            blocks(block_num)%data_fields(dF)%k2        = blocks(block_num)%data_fields(dF)%data_
            ! RHS
            call RHS_2D_block(blocks(block_num)%data_fields(dF)%k2(:,:), blocks(block_num)%dx, blocks(block_num)%dy, g, Bs)

        end do

        !------------------------------
        ! third stage
        do k = 1, N

            block_num                                   = blocks_params%active_list(k)
            blocks(block_num)%data_fields(dF)%data_     = blocks(block_num)%data_fields(dF)%data_old + 0.5_rk * dt * blocks(block_num)%data_fields(dF)%k2

        end do
        ! synchronize ghostnodes
        call synchronize_ghosts()
        do k = 1, N

            block_num                                   = blocks_params%active_list(k)
            blocks(block_num)%data_fields(dF)%k3        = blocks(block_num)%data_fields(dF)%data_
            ! RHS
            call RHS_2D_block(blocks(block_num)%data_fields(dF)%k3(:,:), blocks(block_num)%dx, blocks(block_num)%dy, g, Bs)

        end do

        !------------------------------
        ! fourth stage
        do k = 1, N

            block_num                                   = blocks_params%active_list(k)
            blocks(block_num)%data_fields(dF)%data_     = blocks(block_num)%data_fields(dF)%data_old + dt * blocks(block_num)%data_fields(dF)%k3

        end do
        ! synchronize ghostnodes
        call synchronize_ghosts()
        do k = 1, N

            block_num                                   = blocks_params%active_list(k)
            blocks(block_num)%data_fields(dF)%k4        = blocks(block_num)%data_fields(dF)%data_
            ! RHS
            call RHS_2D_block(blocks(block_num)%data_fields(dF)%k4(:,:), blocks(block_num)%dx, blocks(block_num)%dy, g, Bs)

        end do

        !------------------------------
        ! final stage
        do k = 1, N

            block_num                                   = blocks_params%active_list(k)

            blocks(block_num)%data_fields(dF)%data_     = blocks(block_num)%data_fields(dF)%data_old + dt/6.0_rk * ( blocks(block_num)%data_fields(dF)%k1 &
                                                                            + 2.0_rk*blocks(block_num)%data_fields(dF)%k2 &
                                                                            + 2.0_rk*blocks(block_num)%data_fields(dF)%k3 &
                                                                            + blocks(block_num)%data_fields(dF)%k4 )

        end do

    else
        ! more than one data field
        ! to do

    end if

end subroutine time_step
