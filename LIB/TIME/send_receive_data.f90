! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: send_receive_data.f90
! version: 0.4
! author: msr
!
! send and receive data to synchronize ghost nodes
!
! input:    - heavy data array
!           - params struct
!           - communications list
!           - com_list id
!           - number of communications to send/receive
! output:   - heavy data array
!
! -------------------------------------------------------------------------------------------------------------------------
! dirs = (/'__N', '__E', '__S', '__W', '_NE', '_NW', '_SE', '_SW', 'NNE', 'NNW', 'SSE', 'SSW', 'ENE', 'ESE', 'WNW', 'WSW'/)
! -------------------------------------------------------------------------------------------------------------------------
!
! = log ======================================================================================
!
! 09/11/16 - create for v0.4
! ********************************************************************************************

subroutine send_receive_data(params, hvy_block, com_id, com_list, com_number)

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(in)                  :: params
    ! heavy data array - block data
    real(kind=rk), intent(inout)                    :: hvy_block(:, :, :, :)

    ! com list
    integer(kind=ik), intent(in)                    :: com_list(:, :)
    ! com_list id, number of communications
    integer(kind=ik), intent(in)                    :: com_id, com_number

    ! grid parameter
    integer(kind=ik)                                :: Bs, g

    ! interpolation variables
    real(kind=rk), dimension(:,:), allocatable      :: data_corner, data_corner_fine, data_edge, data_edge_fine

    ! allocation error variable
    integer(kind=ik)                                :: allocate_error

    ! MPI error variable
    integer(kind=ik)                                :: ierr
    ! process rank
    integer(kind=ik)                                :: rank
    ! MPI message tag
    integer(kind=ik)                                :: tag
    ! MPI status
    integer(kind=ik)                                :: status(MPI_status_size)

    ! com list elements
    integer(kind=ik)                                :: my_block, neighbor_block, my_dir, level_diff, my_dest

    ! send/receive buffer list
    real(kind=rk)                                   :: send_buff( params%number_data_fields*params%number_blocks*12*(2*params%number_block_nodes)*params%number_ghost_nodes ), &
                                                       recv_buff( params%number_data_fields*params%number_blocks*12*(2*params%number_block_nodes)*params%number_ghost_nodes )
    ! buffer index
    integer(kind=ik)                                :: buffer_i

    ! loop variable
    integer(kind=ik)                                :: k, l, k_shift, k_start, k_end, dF

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! grid parameter
    Bs    = params%number_block_nodes
    g     = params%number_ghost_nodes

    ! determinate process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    allocate( data_corner( g, g), stat=allocate_error )
    allocate( data_corner_fine( 2*g-1, 2*g-1), stat=allocate_error )
    allocate( data_edge( (Bs+1)/2 + g/2, (Bs+1)/2 + g/2), stat=allocate_error )
    allocate( data_edge_fine( Bs+g, Bs+g), stat=allocate_error )

    send_buff = 9.0e9_rk
    recv_buff = 9.0e9_rk
    buffer_i  = 1

    tag = 0

    data_corner      = 9.0e9_rk
    data_corner_fine = 9.0e9_rk
    data_edge        = 9.0e9_rk
    data_edge_fine   = 9.0e9_rk

!---------------------------------------------------------------------------------------------
! main body

    ! check if proc has to send data in first or second part of com list
    ! also find destination proc
    if ( rank == com_list( com_id, 2 ) ) then
        ! proc has to send next data => no shift
        k_shift = 0
        ! destination
        my_dest = com_list( com_id, 3 )
    else
        ! proc has to shift data, so he first receive data
        k_shift = com_number
        ! destination
        my_dest = com_list( com_id, 2 )
    end if

    ! fill send buffer
    do k = 1 + k_shift, com_number + k_shift

        my_block        = com_list( com_id+k-1, 4 )
        neighbor_block  = com_list( com_id+k-1, 5 )
        my_dir          = com_list( com_id+k-1, 6 )
        level_diff      = com_list( com_id+k-1, 8 )

        select case(my_dir)
            ! '__N'
            case(1)
                do dF = 2, params%number_data_fields+1
                    do l = 1, g
                        send_buff(buffer_i:buffer_i+Bs-1)   = hvy_block( g+l+1, g+1:Bs+g, dF, my_block )
                        buffer_i                            = buffer_i + Bs
                    end do
                end do

            ! '__E'
            case(2)
                do dF = 2, params%number_data_fields+1
                    do l = 1, g
                        send_buff(buffer_i:buffer_i+Bs-1)   = hvy_block( g+1:Bs+g, Bs+g-l, dF, my_block )
                        buffer_i                            = buffer_i + Bs
                    end do
                end do

            ! '__S'
            case(3)
                do dF = 2, params%number_data_fields+1
                    do l = 1, g
                        send_buff(buffer_i:buffer_i+Bs-1)   = hvy_block( Bs+g-l, g+1:Bs+g, dF, my_block )
                        buffer_i                            = buffer_i + Bs
                    end do
                end do

            ! '__W'
            case(4)
                do dF = 2, params%number_data_fields+1
                    do l = 1, g
                        send_buff(buffer_i:buffer_i+Bs-1)   = hvy_block( g+1:Bs+g, g+l+1, dF, my_block )
                        buffer_i                            = buffer_i + Bs
                    end do
                end do

            ! '_NE'
            case(5)
                do dF = 2, params%number_data_fields+1
                    if ( level_diff == 0 ) then
                        ! blocks on same level
                        data_corner = hvy_block( g+2:g+1+g, Bs:Bs-1+g, dF, my_block )

                    elseif ( level_diff == -1 ) then
                        ! sender one level down
                        ! interpolate data
                        ! data to refine
                        data_corner = hvy_block( g+1:g+g, Bs+1:Bs+g, dF, my_block )
                        ! interpolate data
                        call prediction_2D( data_corner , data_corner_fine, params%order_predictor)
                        ! data to synchronize
                        data_corner = data_corner_fine(2:g+1, g-1:2*g-2)

                    elseif ( level_diff == 1) then
                        ! sender one level up
                        data_corner = hvy_block( g+3:g+1+g+g:2, Bs-g:Bs-2+g:2, dF, my_block )

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                    ! send data
                    do l = 1, g
                        send_buff(buffer_i:buffer_i+g-1)    = data_corner(l, 1:g)
                        buffer_i                            = buffer_i + g
                    end do
                end do

            ! '_NW'
            case(6)
                do dF = 2, params%number_data_fields+1
                    if ( level_diff == 0 ) then
                        ! blocks on same level
                        ! loop over all datafields
                        data_corner = hvy_block( g+2:g+1+g, g+2:g+1+g, dF, my_block )

                    elseif ( level_diff == -1 ) then
                        ! sender one level down
                        ! interpolate data
                        ! data to refine
                        data_corner = hvy_block( g+1:g+g, g+1:g+g, dF, my_block )
                        ! interpolate data
                        call prediction_2D( data_corner , data_corner_fine, params%order_predictor)
                        ! data to synchronize
                        data_corner = data_corner_fine(2:g+1, 2:g+1)

                    elseif ( level_diff == 1) then
                        ! sender one level up
                        data_corner = hvy_block( g+3:g+1+g+g:2, g+3:g+1+g+g:2, dF, my_block )

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                    ! send data
                    do l = 1, g
                        send_buff(buffer_i:buffer_i+g-1)    = data_corner(l, 1:g)
                        buffer_i                            = buffer_i + g
                    end do
                end do

            ! '_SE'
            case(7)
                do dF = 2, params%number_data_fields+1
                    if ( level_diff == 0 ) then
                        ! blocks on same level
                        data_corner = hvy_block( Bs:Bs-1+g, Bs:Bs-1+g, dF, my_block )

                    elseif ( level_diff == -1 ) then
                        ! sender one level down
                        ! interpolate data
                        ! data to refine
                        data_corner = hvy_block( Bs+1:Bs+g, Bs+1:Bs+g, dF, my_block )
                        ! interpolate data
                        call prediction_2D( data_corner , data_corner_fine, params%order_predictor)
                        ! data to synchronize
                        data_corner = data_corner_fine(g-1:2*g-2, g-1:2*g-2)

                    elseif ( level_diff == 1) then
                        ! sender one level up
                        data_corner = hvy_block( Bs-g:Bs-2+g:2, Bs-g:Bs-2+g:2, dF, my_block )

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                    ! send data
                    do l = 1, g
                        send_buff(buffer_i:buffer_i+g-1)    = data_corner(l, 1:g)
                        buffer_i                            = buffer_i + g
                    end do
                end do

            ! '_SW'
            case(8)
                do dF = 2, params%number_data_fields+1
                    if ( level_diff == 0 ) then
                        ! blocks on same level
                        data_corner = hvy_block( Bs:Bs-1+g, g+2:g+1+g, dF, my_block )

                    elseif ( level_diff == -1 ) then
                        ! sender one level down
                        ! interpolate data
                        ! data to refine
                        data_corner = hvy_block( Bs+1:Bs+g, g+1:g+g, dF, my_block )
                        ! interpolate data
                        call prediction_2D( data_corner , data_corner_fine, params%order_predictor)
                        ! data to synchronize
                        data_corner = data_corner_fine(g-1:2*g-2, 2:g+1)

                    elseif ( level_diff == 1) then
                        ! sender one level up
                        data_corner = hvy_block( Bs-g:Bs-2+g:2, g+3:g+1+g+g:2, dF, my_block )

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                    ! send data
                    do l = 1, g
                        send_buff(buffer_i:buffer_i+g-1)    = data_corner(l, 1:g)
                        buffer_i                            = buffer_i + g
                    end do
                end do

            ! 'NNE'
            case(9)
                do dF = 2, params%number_data_fields+1
                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_edge = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, my_block )
                        ! interpolate data
                        call prediction_2D( data_edge , data_edge_fine, params%order_predictor)

                        ! send data
                        do l = 1, g
                            send_buff(buffer_i:buffer_i+Bs+g-1)    = data_edge_fine(l+1, 1:Bs+g)
                            buffer_i                               = buffer_i + Bs+g
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! send data
                        do l = 1, g
                            send_buff(buffer_i:buffer_i+(Bs+1)/2-1)   = hvy_block( g+(2*l)+1, g+1:Bs+g:2, dF, my_block )
                            buffer_i                                  = buffer_i + (Bs+1)/2
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if
                end do

            ! 'NNW'
            case(10)
                do dF = 2, params%number_data_fields+1
                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_edge = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, dF, my_block )
                        ! interpolate data
                        call prediction_2D( data_edge , data_edge_fine, params%order_predictor)

                        ! send data
                        do l = 1, g
                            send_buff(buffer_i:buffer_i+Bs+g-1)    = data_edge_fine(l+1, 1:Bs+g)
                            buffer_i                               = buffer_i + Bs+g
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! send data
                        do l = 1, g
                            send_buff(buffer_i:buffer_i+(Bs+1)/2-1)   = hvy_block( g+(2*l)+1, g+1:Bs+g:2, dF, my_block )
                            buffer_i                                  = buffer_i + (Bs+1)/2
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if
                end do

            ! 'SSE'
            case(11)
                do dF = 2, params%number_data_fields+1
                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_edge = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, my_block )
                        ! interpolate data
                        call prediction_2D( data_edge , data_edge_fine, params%order_predictor)
                        ! send data
                        do l = 1, g
                            send_buff(buffer_i:buffer_i+Bs+g-1)    = data_edge_fine(Bs+g-l, 1:Bs+g)
                            buffer_i                               = buffer_i + Bs+g
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! send data
                        do l = 1, g
                            send_buff(buffer_i:buffer_i+(Bs+1)/2-1)   = hvy_block( Bs+g-(2*l), g+1:Bs+g:2, dF, my_block )
                            buffer_i                                  = buffer_i + (Bs+1)/2
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if
                end do

            ! 'SSW'
            case(12)
                do dF = 2, params%number_data_fields+1
                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_edge = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, my_block )
                        ! interpolate data
                        call prediction_2D( data_edge , data_edge_fine, params%order_predictor)
                        ! send data
                        do l = 1, g
                            send_buff(buffer_i:buffer_i+Bs+g-1)    = data_edge_fine(Bs+g-l, 1:Bs+g)
                            buffer_i                               = buffer_i + Bs+g
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! send data
                        do l = 1, g
                            send_buff(buffer_i:buffer_i+(Bs+1)/2-1)   = hvy_block( Bs+g-(2*l), g+1:Bs+g:2, dF, my_block )
                            buffer_i                                  = buffer_i + (Bs+1)/2
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if
                end do

            ! 'ENE'
            case(13)
                do dF = 2, params%number_data_fields+1
                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_edge = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, my_block )
                        ! interpolate data
                        call prediction_2D( data_edge , data_edge_fine, params%order_predictor)
                        ! send data
                        do l = 1, g
                            send_buff(buffer_i:buffer_i+Bs+g-1)    = data_edge_fine(1:Bs+g, Bs+l-1)
                            buffer_i                               = buffer_i + Bs+g
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! send data
                        do l = 1, g
                            send_buff(buffer_i:buffer_i+(Bs+1)/2-1)   = hvy_block( g+1:Bs+g:2, Bs-g+2*l-2, dF, my_block )
                            buffer_i                                  = buffer_i + (Bs+1)/2
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if
                end do

            ! 'ESE'
            case(14)
                do dF = 2, params%number_data_fields+1
                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_edge = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, my_block )
                        ! interpolate data
                        call prediction_2D( data_edge , data_edge_fine, params%order_predictor)
                        ! send data
                        do l = 1, g
                            send_buff(buffer_i:buffer_i+Bs+g-1)    = data_edge_fine(1:Bs+g, Bs+l-1)
                            buffer_i                               = buffer_i + Bs+g
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        do l = 1, g
                            send_buff(buffer_i:buffer_i+(Bs+1)/2-1)   = hvy_block( g+1:Bs+g:2, Bs-g+2*l-2, dF, my_block )
                            buffer_i                                  = buffer_i + (Bs+1)/2
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if
                end do

            ! 'WNW'
            case(15)
                do dF = 2, params%number_data_fields+1
                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_edge = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, dF, my_block )
                        ! interpolate data
                        call prediction_2D( data_edge , data_edge_fine, params%order_predictor)
                        ! send data
                        do l = 1, g
                            send_buff(buffer_i:buffer_i+Bs+g-1)    = data_edge_fine(1:Bs+g, l+1)
                            buffer_i                               = buffer_i + Bs+g
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        do l = 1, g
                            send_buff(buffer_i:buffer_i+(Bs+1)/2-1)   = hvy_block( g+1:Bs+g:2, g+(2*l)+1, dF, my_block )
                            buffer_i                                  = buffer_i + (Bs+1)/2
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if
                end do

            ! 'WSW'
            case(16)
                do dF = 2, params%number_data_fields+1
                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_edge = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, my_block )
                        ! interpolate data
                        call prediction_2D( data_edge , data_edge_fine, params%order_predictor)
                        ! send data
                        do l = 1, g
                            send_buff(buffer_i:buffer_i+Bs+g-1)    = data_edge_fine(1:Bs+g, l+1)
                            buffer_i                               = buffer_i + Bs+g
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! send data
                        do l = 1, g
                            send_buff(buffer_i:buffer_i+(Bs+1)/2-1)   = hvy_block( g+1:Bs+g:2, g+(2*l)+1, dF, my_block )
                            buffer_i                                  = buffer_i + (Bs+1)/2
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if
                end do

        end select
    end do

    ! send/receive data
    call MPI_Sendrecv( send_buff, params%number_data_fields*(Bs+g)*g*com_number, MPI_REAL8, my_dest, tag, recv_buff, params%number_data_fields*(Bs+g)*g*com_number, MPI_REAL8, my_dest, tag, MPI_COMM_WORLD, status, ierr)

    k_start = 1+com_number-k_shift
    k_end   = com_number+com_number-k_shift

    ! reset buffer index
    buffer_i  = 1

    ! reset
    data_corner      = 9.0e9_rk
    data_corner_fine = 9.0e9_rk
    data_edge        = 9.0e9_rk
    data_edge_fine   = 9.0e9_rk

    ! write received data in block data
    do k = k_start, k_end

        my_block        = com_list( com_id+k-1, 5 )
        neighbor_block  = com_list( com_id+k-1, 4 )
        my_dir          = com_list( com_id+k-1, 7 )
        level_diff      = com_list( com_id+k-1, 8 )

        select case(my_dir)
            ! '__N'
            case(1)
                do dF = 2, params%number_data_fields+1
                    do l = 1, g
                        hvy_block( Bs+g+l, g+1:Bs+g, dF, my_block )     = recv_buff(buffer_i:buffer_i+Bs-1)
                        buffer_i                                         = buffer_i + Bs
                    end do
                end do

            ! '__E'
            case(2)
                do dF = 2, params%number_data_fields+1
                    do l = 1, g
                        hvy_block( g+1:Bs+g, g+1-l, dF , my_block)      = recv_buff(buffer_i:buffer_i+Bs-1)
                        buffer_i                                         = buffer_i + Bs
                    end do
                end do

            ! '__S'
            case(3)
                do dF = 2, params%number_data_fields+1
                    do l = 1, g
                        hvy_block( g+1-l, g+1:Bs+g, dF, my_block )      = recv_buff(buffer_i:buffer_i+Bs-1)
                        buffer_i                                         = buffer_i + Bs
                    end do
                end do

            ! '__W'
            case(4)
                do dF = 2, params%number_data_fields+1
                    do l = 1, g
                        hvy_block( g+1:Bs+g, Bs+g+l, dF, my_block )     = recv_buff(buffer_i:buffer_i+Bs-1)
                        buffer_i                                         = buffer_i + Bs
                    end do
                end do

            ! '_NE'
            case(5)
                do dF = 2, params%number_data_fields+1
                    ! receive data
                    do l = 1, g
                        data_corner(l, 1:g) = recv_buff(buffer_i:buffer_i+g-1)
                        buffer_i            = buffer_i + g
                    end do
                    ! write data
                    hvy_block( Bs+g+1:Bs+g+g, 1:g, dF, my_block ) = data_corner
                end do

            ! '_NW'
            case(6)
                do dF = 2, params%number_data_fields+1
                    ! receive data
                    do l = 1, g
                        data_corner(l, 1:g) = recv_buff(buffer_i:buffer_i+g-1)
                        buffer_i            = buffer_i + g
                    end do
                    ! write data
                    hvy_block( Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, dF, my_block ) = data_corner
                end do

            ! '_SE'
            case(7)
                do dF = 2, params%number_data_fields+1
                    ! receive data
                    do l = 1, g
                        data_corner(l, 1:g) = recv_buff(buffer_i:buffer_i+g-1)
                        buffer_i            = buffer_i + g
                    end do
                    ! write data
                    hvy_block( 1:g, 1:g, dF, my_block ) = data_corner
                end do

            ! '_SW'
            case(8)
                do dF = 2, params%number_data_fields+1
                    ! receive data
                    do l = 1, g
                        data_corner(l, 1:g) = recv_buff(buffer_i:buffer_i+g-1)
                        buffer_i            = buffer_i + g
                    end do
                    ! write data
                    hvy_block( 1:g, Bs+g+1:Bs+g+g, dF, my_block ) = data_corner
                end do

            ! 'NNE'
            case(9)
                do dF = 2, params%number_data_fields+1
                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! receive data
                        do l = 1, g
                            data_edge_fine(l, 1:Bs+g)         = recv_buff(buffer_i:buffer_i+Bs+g-1)
                            buffer_i                          = buffer_i + Bs+g
                        end do

                        ! write data
                        hvy_block( Bs+g+1:Bs+g+g, 1:Bs+g, dF, my_block ) = data_edge_fine(1:g, 1:Bs+g)

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! receive data
                        do l = 1, g
                            hvy_block( Bs+g+l, g+(Bs+1)/2:Bs+g, dF, my_block )  = recv_buff(buffer_i:buffer_i+(Bs+1)/2-1)
                            buffer_i                                             = buffer_i + (Bs+1)/2
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if
                end do

            ! 'NNW'
            case(10)
                do dF = 2, params%number_data_fields+1
                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! receive data
                        do l = 1, g
                            data_edge_fine(l, 1:Bs+g)         = recv_buff(buffer_i:buffer_i+Bs+g-1)
                            buffer_i                          = buffer_i + Bs+g
                        end do
                        ! write data
                        hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+2*g, dF, my_block ) = data_edge_fine(1:g, 1:Bs+g)

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! receive data
                        do l = 1, g
                            hvy_block( Bs+g+l, g+1:g+(Bs+1)/2, dF, my_block )  = recv_buff(buffer_i:buffer_i+(Bs+1)/2-1)
                            buffer_i                                            = buffer_i + (Bs+1)/2
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if
                end do

            ! 'SSE'
            case(11)
                do dF = 2, params%number_data_fields+1
                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! receive data
                        do l = 1, g
                            data_edge_fine(g-l+1, 1:Bs+g)     = recv_buff(buffer_i:buffer_i+Bs+g-1)
                            buffer_i                          = buffer_i + Bs+g
                        end do
                        ! write data
                        hvy_block( 1:g, 1:Bs+g, dF, my_block ) = data_edge_fine(1:g, 1:Bs+g)

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! receive data
                        do l = 1, g
                            hvy_block( g-l+1, g+(Bs+1)/2:Bs+g, dF, my_block )  = recv_buff(buffer_i:buffer_i+(Bs+1)/2-1)
                            buffer_i                                            = buffer_i + (Bs+1)/2
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if
                end do

            ! 'SSW'
            case(12)
                do dF = 2, params%number_data_fields+1
                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! receive data
                        do l = 1, g
                            data_edge_fine(g-l+1, 1:Bs+g)     = recv_buff(buffer_i:buffer_i+Bs+g-1)
                            buffer_i                          = buffer_i + Bs+g
                        end do
                        ! write data
                        hvy_block( 1:g, g+1:Bs+2*g, dF, my_block ) = data_edge_fine(1:g, 1:Bs+g)

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! receive data
                        do l = 1, g
                            hvy_block( g-l+1, g+1:g+(Bs+1)/2, dF, my_block )  = recv_buff(buffer_i:buffer_i+(Bs+1)/2-1)
                            buffer_i                                           = buffer_i + (Bs+1)/2
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if
                end do

            ! 'ENE'
            case(13)
                do dF = 2, params%number_data_fields+1
                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! receive data
                        do l = 1, g
                            data_edge_fine(1:Bs+g, Bs+l)     = recv_buff(buffer_i:buffer_i+Bs+g-1)
                            buffer_i                         = buffer_i + Bs+g
                        end do
                        ! write data
                        hvy_block( g+1:Bs+2*g, 1:g, dF, my_block ) = data_edge_fine(1:Bs+g, Bs+1:Bs+g)

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! receive data
                        do l = 1, g
                            hvy_block( g+1:g+(Bs+1)/2, l, dF, my_block )  = recv_buff(buffer_i:buffer_i+(Bs+1)/2-1)
                            buffer_i                                       = buffer_i + (Bs+1)/2
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if
                end do

            ! 'ESE'
            case(14)
                do dF = 2, params%number_data_fields+1
                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! receive data
                        do l = 1, g
                            data_edge_fine(1:Bs+g, Bs+l)     = recv_buff(buffer_i:buffer_i+Bs+g-1)
                            buffer_i                          = buffer_i + Bs+g
                        end do
                        ! write data
                        hvy_block( 1:Bs+g, 1:g, dF, my_block ) = data_edge_fine(1:Bs+g, Bs+1:Bs+g)

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! receive data
                        do l = 1, g
                            hvy_block( g+(Bs+1)/2:Bs+g, l, dF, my_block )  = recv_buff(buffer_i:buffer_i+(Bs+1)/2-1)
                            buffer_i                                        = buffer_i + (Bs+1)/2
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if
                end do

            ! 'WNW'
            case(15)
                do dF = 2, params%number_data_fields+1
                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! receive data
                        do l = 1, g
                            data_edge_fine(1:Bs+g, l)     = recv_buff(buffer_i:buffer_i+Bs+g-1)
                            buffer_i                      = buffer_i + Bs+g
                        end do
                        ! write data
                        hvy_block( g+1:Bs+2*g, Bs+g+1:Bs+g+g, dF, my_block ) = data_edge_fine(1:Bs+g, 1:g)

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! receive data
                        do l = 1, g
                            hvy_block( g+1:g+(Bs+1)/2, Bs+g+l, dF, my_block )  = recv_buff(buffer_i:buffer_i+(Bs+1)/2-1)
                            buffer_i                                            = buffer_i + (Bs+1)/2
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if
                end do

            ! 'WSW'
            case(16)
                do dF = 2, params%number_data_fields+1
                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! receive data
                        do l = 1, g
                            data_edge_fine(1:Bs+g, l)     = recv_buff(buffer_i:buffer_i+Bs+g-1)
                            buffer_i                      = buffer_i + Bs+g
                        end do
                        ! write data
                        hvy_block( 1:Bs+g, Bs+g+1:Bs+g+g, dF, my_block ) = data_edge_fine(1:Bs+g, 1:g)

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! receive data
                        do l = 1, g
                            hvy_block( g+(Bs+1)/2:Bs+g, Bs+g+l, dF, my_block )  = recv_buff(buffer_i:buffer_i+(Bs+1)/2-1)
                            buffer_i                                             = buffer_i + (Bs+1)/2
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if
                end do

        end select

    end do

    ! clean up
    deallocate( data_corner, stat=allocate_error )
    deallocate( data_corner_fine, stat=allocate_error )
    deallocate( data_edge, stat=allocate_error )
    deallocate( data_edge_fine, stat=allocate_error )

end subroutine send_receive_data
