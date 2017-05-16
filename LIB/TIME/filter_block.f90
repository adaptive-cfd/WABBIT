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

subroutine filter_block( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

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
    !> heavy data array - neighbor data
    integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)

    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    ! loop variables
    integer(kind=ik)                    :: k, i, j, l, dF, N_dF

    ! grid parameter
    integer(kind=ik)                    :: Bs, g

    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0, sub_t1, time_sum

    ! stencil arrays
    real(kind=rk)                       :: stencil_5(5), stencil_7(7), stencil_9(9), stencil_11(11)

    ! filter position (array postion of value to filter)
    integer(kind=ik)                    :: filter_pos

    ! filtered values
    real(kind=rk)                       :: phi_tilde(3)


!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    time_sum = 0.0_rk

    ! grid parameter
    Bs    = params%number_block_nodes
    g     = params%number_ghost_nodes

    N_dF  = params%number_data_fields

    ! initialize stencils
    stencil_5  = (/ -1.0_rk/ 16.0_rk, &
                     1.0_rk/  4.0_rk, &
                    -3.0_rk/  8.0_rk, &
                     1.0_rk/  4.0_rk, &
                    -1.0_rk/ 16.0_rk/)
    stencil_7  = (/  1.0_rk/ 64.0_rk, &
                    -3.0_rk/ 32.0_rk, &
                    15.0_rk/ 64.0_rk, &
                    -5.0_rk/ 16.0_rk, &
                    15.0_rk/ 64.0_rk, &
                    -3.0_rk/ 32.0_rk, &
                     1.0_rk/ 64.0_rk/)
    stencil_9  = (/ -1.0_rk/256.0_rk, &
                     1.0_rk/ 32.0_rk, &
                    -7.0_rk/ 64.0_rk, &
                     7.0_rk/ 32.0_rk, &
                   -35.0_rk/128.0_rk, &
                     7.0_rk/ 32.0_rk, &
                    -7.0_rk/ 64.0_rk, &
                     1.0_rk/ 32.0_rk, &
                    -1.0_rk/256.0_rk/)
    stencil_11 = (/  1.0_rk/1024.0_rk, &
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

    filter_pos = 0

    ! set filter pos
    select case(params%filter_type)
        case('explicit_5pt')
            filter_pos = 3

        case('explicit_7pt')
            filter_pos = 4

        case('explicit_9pt')
            filter_pos = 5

        case('explicit_11pt')
            filter_pos = 6
            
        case('no_filter')
            ! do nothing..

        case default
            write(*,'(80("_"))')
            write(*,*) "ERROR: filter type is unknown"
            write(*,*) params%filter_type
            stop

    end select

!---------------------------------------------------------------------------------------------
! main body

    ! check ghost nodes number
    if ( g < (filter_pos-1) ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: number of ghost nodes is too low for filtering"
        stop
    end if

    ! synchronize ghostnodes
    call synchronize_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

    ! start time
    sub_t0 = MPI_Wtime()

    ! select filter method
    select case(params%filter_type)
        case('explicit_5pt')
            ! loop over all datafields
            do dF = 2, N_dF+1
                ! loop over all active heavy data blocks
                do k = 1, hvy_n

                    if ( params%threeD_case ) then
                        ! 3D:
                        do i = g+1, Bs+g
                            do j = g+1, Bs+g
                                do l = g+1, Bs+g

                                    ! x direction
                                    call filter_1D( hvy_block(i-(filter_pos-1):i+(filter_pos-1), j, l, dF, hvy_active(k) ), phi_tilde(1), stencil_5 )
                                    ! y direction
                                    call filter_1D( hvy_block(i, j-(filter_pos-1):j+(filter_pos-1), l, dF, hvy_active(k) ), phi_tilde(2), stencil_5 )
                                    ! z direction
                                    call filter_1D( hvy_block(i, j, l-(filter_pos-1):l+(filter_pos-1), dF, hvy_active(k) ), phi_tilde(3), stencil_5 )

                                    ! filter
                                    hvy_block(i, j, l, dF, hvy_active(k) ) = hvy_block(i, j, l, dF, hvy_active(k) ) + phi_tilde(1) + phi_tilde(2) + phi_tilde(3)

                                end do
                            end do
                        end do

                    else
                        ! 2D:
                        do i = g+1, Bs+g
                            do j = g+1, Bs+g

                                ! x direction
                                call filter_1D( hvy_block(i-(filter_pos-1):i+(filter_pos-1), j, 1, dF, hvy_active(k) ), phi_tilde(1), stencil_5 )
                                ! y direction
                                call filter_1D( hvy_block(i, j-(filter_pos-1):j+(filter_pos-1), 1, dF, hvy_active(k) ), phi_tilde(2), stencil_5 )

                                ! filter
                                hvy_block(i, j, 1, dF, hvy_active(k) ) = hvy_block(i, j, 1, dF, hvy_active(k) ) + phi_tilde(1) + phi_tilde(2)

                            end do
                        end do

                    end if

                end do
            end do

        case('explicit_7pt')
            ! loop over all datafields
            do dF = 2, N_dF+1
                ! loop over all active heavy data blocks
                do k = 1, hvy_n

                    if ( params%threeD_case ) then
                        ! 3D:
                        do i = g+1, Bs+g
                            do j = g+1, Bs+g
                                do l = g+1, Bs+g

                                    ! x direction
                                    call filter_1D( hvy_block(i-(filter_pos-1):i+(filter_pos-1), j, l, dF, hvy_active(k) ), phi_tilde(1), stencil_7 )
                                    ! y direction
                                    call filter_1D( hvy_block(i, j-(filter_pos-1):j+(filter_pos-1), l, dF, hvy_active(k) ), phi_tilde(2), stencil_7 )
                                    ! z direction
                                    call filter_1D( hvy_block(i, j, l-(filter_pos-1):l+(filter_pos-1), dF, hvy_active(k) ), phi_tilde(3), stencil_7 )

                                    ! filter
                                    hvy_block(i, j, l, dF, hvy_active(k) ) = hvy_block(i, j, l, dF, hvy_active(k) ) + phi_tilde(1) + phi_tilde(2) + phi_tilde(3)

                                end do
                            end do
                        end do

                    else
                        ! 2D:
                        do i = g+1, Bs+g
                            do j = g+1, Bs+g

                                ! x direction
                                call filter_1D( hvy_block(i-(filter_pos-1):i+(filter_pos-1), j, 1, dF, hvy_active(k) ), phi_tilde(1), stencil_7 )
                                ! y direction
                                call filter_1D( hvy_block(i, j-(filter_pos-1):j+(filter_pos-1), 1, dF, hvy_active(k) ), phi_tilde(2), stencil_7 )

                                ! filter
                                hvy_block(i, j, 1, dF, hvy_active(k) ) = hvy_block(i, j, 1, dF, hvy_active(k) ) + phi_tilde(1) + phi_tilde(2)

                            end do
                        end do

                    end if

                end do
            end do

        case('explicit_9pt')
            ! loop over all datafields
            do dF = 2, N_dF+1
                ! loop over all active heavy data blocks
                do k = 1, hvy_n

                    if ( params%threeD_case ) then
                        ! 3D:
                        do i = g+1, Bs+g
                            do j = g+1, Bs+g
                                do l = g+1, Bs+g

                                    ! x direction
                                    call filter_1D( hvy_block(i-(filter_pos-1):i+(filter_pos-1), j, l, dF, hvy_active(k) ), phi_tilde(1), stencil_9 )
                                    ! y direction
                                    call filter_1D( hvy_block(i, j-(filter_pos-1):j+(filter_pos-1), l, dF, hvy_active(k) ), phi_tilde(2), stencil_9 )
                                    ! z direction
                                    call filter_1D( hvy_block(i, j, l-(filter_pos-1):l+(filter_pos-1), dF, hvy_active(k) ), phi_tilde(3), stencil_9 )

                                    ! filter
                                    hvy_block(i, j, l, dF, hvy_active(k) ) = hvy_block(i, j, l, dF, hvy_active(k) ) + phi_tilde(1) + phi_tilde(2) + phi_tilde(3)

                                end do
                            end do
                        end do

                    else
                        ! 2D:
                        do i = g+1, Bs+g
                            do j = g+1, Bs+g

                                ! x direction
                                call filter_1D( hvy_block(i-(filter_pos-1):i+(filter_pos-1), j, 1, dF, hvy_active(k) ), phi_tilde(1), stencil_9 )
                                ! y direction
                                call filter_1D( hvy_block(i, j-(filter_pos-1):j+(filter_pos-1), 1, dF, hvy_active(k) ), phi_tilde(2), stencil_9 )

                                ! filter
                                hvy_block(i, j, 1, dF, hvy_active(k) ) = hvy_block(i, j, 1, dF, hvy_active(k) ) + phi_tilde(1) + phi_tilde(2)

                            end do
                        end do

                    end if

                end do
            end do

        case('explicit_11pt')
            ! loop over all datafields
            do dF = 2, N_dF+1
                ! loop over all active heavy data blocks
                do k = 1, hvy_n

                    if ( params%threeD_case ) then
                        ! 3D:
                        do i = g+1, Bs+g
                            do j = g+1, Bs+g
                                do l = g+1, Bs+g

                                    ! x direction
                                    call filter_1D( hvy_block(i-(filter_pos-1):i+(filter_pos-1), j, l, dF, hvy_active(k) ), phi_tilde(1), stencil_11)
                                    ! y direction
                                    call filter_1D( hvy_block(i, j-(filter_pos-1):j+(filter_pos-1), l, dF, hvy_active(k) ), phi_tilde(2), stencil_11 )
                                    ! z direction
                                    call filter_1D( hvy_block(i, j, l-(filter_pos-1):l+(filter_pos-1), dF, hvy_active(k) ), phi_tilde(3), stencil_11 )

                                    ! filter
                                    hvy_block(i, j, l, dF, hvy_active(k) ) = hvy_block(i, j, l, dF, hvy_active(k) ) + phi_tilde(1) + phi_tilde(2) + phi_tilde(3)

                                end do
                            end do
                        end do

                    else
                        ! 2D:
                        do i = g+1, Bs+g
                            do j = g+1, Bs+g

                                ! x direction
                                call filter_1D( hvy_block(i-(filter_pos-1):i+(filter_pos-1), j, 1, dF, hvy_active(k) ), phi_tilde(1), stencil_11 )
                                ! y direction
                                call filter_1D( hvy_block(i, j-(filter_pos-1):j+(filter_pos-1), 1, dF, hvy_active(k) ), phi_tilde(2), stencil_11 )

                                ! filter
                                hvy_block(i, j, 1, dF, hvy_active(k) ) = hvy_block(i, j, 1, dF, hvy_active(k) ) + phi_tilde(1) + phi_tilde(2)

                            end do
                        end do

                    end if

                end do
            end do

        case('no_filter')
            ! do nothing

        case default
            write(*,'(80("_"))')
            write(*,*) "ERROR: filter type is unknown"
            write(*,*) params%filter_type
            stop

    end select

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
