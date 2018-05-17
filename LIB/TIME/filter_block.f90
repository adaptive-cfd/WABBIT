!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name filter_block.f90
!> \version 0.5
!> \author msr
!
!> \brief filter heavy data
!
!>
!! input:    - params, heavy data array \n
!! output:   - filtered heavy data array \n
!!
!!
!! = log ======================================================================================
!! \n
!! 27/03/17 - create
!! 02/05/17 - filter (x,y,z)-values separatly
!
! ********************************************************************************************

subroutine filter_block( params, lgt_block, hvy_block, hvy_synch, hvy_neighbor, hvy_active, hvy_n, com_lists, com_matrix, int_send_buffer, int_receive_buffer, real_send_buffer, real_receive_buffer )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params

    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy synch array
    integer(kind=1), intent(inout)      :: hvy_synch(:, :, :, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)

    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    ! communication lists:
    integer(kind=ik), intent(inout)     :: com_lists(:, :, :, :)

    ! communications matrix:
    integer(kind=ik), intent(inout)     :: com_matrix(:,:,:)

    ! send/receive buffer, integer and real
    integer(kind=ik), intent(inout)      :: int_send_buffer(:,:), int_receive_buffer(:,:)
    real(kind=rk), intent(inout)         :: real_send_buffer(:,:), real_receive_buffer(:,:)

    ! loop variables
    integer(kind=ik)                    :: k, i, j, l, dF, N_dF, lgt_id

    ! grid parameter
    integer(kind=ik)                    :: Bs, g

    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0, sub_t1, time_sum

    ! stencil array, note: size is fixed
    real(kind=rk)                       :: stencil(19)

    ! filter position (array postion of value to filter)
    integer(kind=ik)                    :: stencil_size

    ! filtered values and array for old block data
    real(kind=rk)                       :: phi_tilde(3)
    real(kind=rk), allocatable          :: block_old(:, :, :)

    ! spacing and origin of a block
    real(kind=rk)                       :: xx0(1:3), ddx(1:3)


!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    time_sum = 0.0_rk

    ! grid parameter
    Bs    = params%number_block_nodes
    g     = params%number_ghost_nodes

    N_dF  = params%number_data_fields

    ! allocate old data array
    allocate( block_old(Bs+2*g, Bs+2*g, Bs+2*g) )

    stencil_size = 0

    ! set filter pos
    select case(params%filter_type)
        case('explicit_5pt')
            stencil_size = 5
            stencil(1:stencil_size) = (/ -1.0_rk/ 16.0_rk, &
                                          1.0_rk/  4.0_rk, &
                                         -3.0_rk/  8.0_rk, &
                                          1.0_rk/  4.0_rk, &
                                         -1.0_rk/ 16.0_rk/)

        case('explicit_7pt')
            stencil_size = 7
            stencil(1:stencil_size) = (/  1.0_rk/ 64.0_rk, &
                                         -3.0_rk/ 32.0_rk, &
                                         15.0_rk/ 64.0_rk, &
                                         -5.0_rk/ 16.0_rk, &
                                         15.0_rk/ 64.0_rk, &
                                         -3.0_rk/ 32.0_rk, &
                                          1.0_rk/ 64.0_rk/)

        case('explicit_9pt')
            stencil_size = 9
            stencil(1:stencil_size) = (/ -1.0_rk/256.0_rk, &
                                          1.0_rk/ 32.0_rk, &
                                         -7.0_rk/ 64.0_rk, &
                                          7.0_rk/ 32.0_rk, &
                                        -35.0_rk/128.0_rk, &
                                          7.0_rk/ 32.0_rk, &
                                         -7.0_rk/ 64.0_rk, &
                                          1.0_rk/ 32.0_rk, &
                                         -1.0_rk/256.0_rk/)

        case('explicit_11pt')
            stencil_size = 11
            stencil(1:stencil_size) = (/  1.0_rk/1024.0_rk, &
                                         -5.0_rk/ 512.0_rk, &
                                         45.0_rk/1024.0_rk, &
                                        -15.0_rk/ 128.0_rk, &
                                        105.0_rk/ 512.0_rk, &
                                        -63.0_rk/ 256.0_rk, &
                                        105.0_rk/ 512.0_rk, &
                                        -15.0_rk/ 128.0_rk, &
                                         45.0_rk/1024.0_rk, &
                                         -5.0_rk/ 512.0_rk, &
                                          1.0_rk/1024.0_rk/)
            
        case('no_filter')
            ! do nothing..

        case('wavelet')
            ! do nothing..

        case('bogey_shock')
            ! do nothing

        case default
            write(*,'(80("_"))')
            write(*,*) "ERROR: filter type is unknown"
            write(*,*) params%filter_type
            stop

    end select

!---------------------------------------------------------------------------------------------
! main body

    ! check ghost nodes number
    if ( g < ((stencil_size+1)/2-1) ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: number of ghost nodes is too low for filtering"
        stop
    end if

    ! synchronize ghostnodes
    call synchronize_ghosts( params, lgt_block, hvy_block, hvy_synch, hvy_neighbor, hvy_active, hvy_n, com_lists, com_matrix, .false., int_send_buffer, int_receive_buffer, real_send_buffer, real_receive_buffer )
    ! start time
    sub_t0 = MPI_Wtime()

    ! loop over all active heavy data blocks
    do k = 1, hvy_n

        ! loop over all datafields
        do dF = 1, N_dF

            ! switch case
            ! explicit filter -> stencil_size /= 0
            ! wavelet filter -> stencil_size == 0
            if (stencil_size /= 0) then
                ! explicit filter

                ! save old block data
                block_old = hvy_block(:, :, :, dF, hvy_active(k) )

                ! 3D or 2D case
                if ( params%threeD_case ) then
                    ! 3D
                    ! loop over block data
                    do i = g+1, Bs+g
                        do j = g+1, Bs+g
                            do l = g+1, Bs+g

                                ! x direction
                                call filter_1D( block_old(i-( (stencil_size+1)/2-1):i+( (stencil_size+1)/2-1), j, l ), phi_tilde(1), stencil(1:stencil_size) )
                                ! y direction
                                call filter_1D( block_old(i, j-( (stencil_size+1)/2-1):j+( (stencil_size+1)/2-1), l ), phi_tilde(2), stencil(1:stencil_size) )
                                ! z direction
                                call filter_1D( block_old(i, j, l-( (stencil_size+1)/2-1):l+( (stencil_size+1)/2-1) ), phi_tilde(3), stencil(1:stencil_size) )

                                ! filter
                                hvy_block(i, j, l, dF, hvy_active(k) ) = hvy_block(i, j, l, dF, hvy_active(k) ) + phi_tilde(1) + phi_tilde(2) + phi_tilde(3)

                            end do
                        end do
                    end do

                else
                    ! 2D
                    ! loop over block data
                    do i = g+1, Bs+g
                        do j = g+1, Bs+g

                            ! x direction
                            call filter_1D( block_old(i-( (stencil_size+1)/2-1):i+( (stencil_size+1)/2-1), j, 1 ), phi_tilde(1), stencil(1:stencil_size) )
                            ! y direction
                            call filter_1D( block_old(i, j-( (stencil_size+1)/2-1):j+( (stencil_size+1)/2-1), 1 ), phi_tilde(2), stencil(1:stencil_size) )

                            ! filter
                            hvy_block(i, j, 1, dF, hvy_active(k) ) = hvy_block(i, j, 1, dF, hvy_active(k) ) + phi_tilde(1) + phi_tilde(2)

                        end do
                    end do

                end if

            elseif (stencil_size == 0) then

                select case(params%filter_type)

                    case('wavelet')
                        ! wavelet filter
                        call wavelet_filter( params, hvy_block(:, :, :, dF, hvy_active(k) ))

                    case('bogey_shock')
                        ! convert given hvy_id to lgt_id for block spacing routine
                        call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )
                        ! get block spacing for RHS
                        call get_block_spacing_origin( params, lgt_id, lgt_block, xx0, ddx )

                        ! shock filter
                        if ( dF == 1 ) then
                            call bogey_filter( params, params%number_ghost_nodes, params%number_block_nodes, params%number_data_fields, hvy_block(:, :, :, 1:N_dF, hvy_active(k) ), xx0, ddx, lgt_block(lgt_id, params%max_treelevel+1) )
                        end if

                end select

            end if

        end do

    end do

    ! clean up
    deallocate(block_old)

    ! end time
    sub_t1   = MPI_Wtime()
    time_sum = time_sum + (sub_t1 - sub_t0)
    ! write time
    if ( params%debug ) then
        ! find free or corresponding line
        k = 1
        do while ( debug%name_comp_time(k) /= "---" )
            ! entry for current subroutine exists
            if ( debug%name_comp_time(k) == "filter block" ) exit
            k = k + 1
        end do
        ! write time
        debug%name_comp_time(k) = "filter block"
        debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1
        debug%comp_time(k, 2)   = debug%comp_time(k, 2) + time_sum
    end if

end subroutine filter_block
