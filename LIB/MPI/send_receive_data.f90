! ********************************
! WABBIT
! --------------------------------
!
! send and receive data to synchronize ghost nodes
!
! name: send_receive_data.f90
! date: 26.10.2016
! author: msr
! version: 0.3
!
! ********************************

subroutine send_receive_data(com_id, plan_type, com_number, com_list, com_list_N)

    use mpi
    use module_blocks
    use module_interpolation

    implicit none

    integer(kind=ik), intent(in)                    :: com_list_N

    real(kind=rk) , dimension(1000000)              :: send_buff, recv_buff

    integer, intent(in)                             :: com_id, com_number, plan_type
    integer, dimension(com_list_N, 7), intent(in)   :: com_list
    integer                                         :: Bs, g, k, k_start, k_end, l, buffer_i, dF
    integer                                         :: my_block, neighbor_block, my_dir, tag, ierr, rank, k_shift, my_dest, allocate_error
    integer                                         :: status(MPI_status_size)

    real(kind=rk), dimension(:,:), allocatable      :: data_corner, data_corner_fine, data_edge1, data_edge1_coarse, data_edge1_fine, data_edge2_coarse

    tag      = 0

    send_buff = 9.0e9_rk
    recv_buff = 9.0e9_rk
    buffer_i  = 1

    g  = blocks_params%number_ghost_nodes
    Bs = blocks_params%size_block
    dF = 1

    allocate( data_corner( g, g), stat=allocate_error )
    allocate( data_corner_fine( 2*g-1, 2*g-1), stat=allocate_error )
    allocate( data_edge1( (Bs+1)/2 + g/2, (Bs+1)/2 + g/2), stat=allocate_error )
    allocate( data_edge1_fine( Bs+g, Bs+g), stat=allocate_error )
    allocate( data_edge1_coarse( g, (Bs+1)/2), stat=allocate_error )
    allocate( data_edge2_coarse( (Bs+1)/2, g), stat=allocate_error )

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    ! check if proc has to send data in first or second part of com list
    ! also find destination proc
    if ( rank == com_list( com_id, 2 ) ) then
        ! proc has to send next data => no shift (case also valid for internal communication)
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
    do k = 1+k_shift, com_number+k_shift

        my_block        = com_list( com_id+k-1, 4 )
        neighbor_block  = com_list( com_id+k-1, 5 )
        my_dir          = com_list( com_id+k-1, 6 )

        select case(blocks(my_block)%neighbor_dir(my_dir))
            case('__N') ! north
                do l = 1, g
                    send_buff(buffer_i:buffer_i+Bs-1)   = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g+l+1,g+1:Bs+g)
                    buffer_i                            = buffer_i + Bs
                end do

            case('__W') ! west
                do l = 1, g
                    send_buff(buffer_i:buffer_i+Bs-1)   = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g+1:Bs+g, g+l+1)
                    buffer_i                            = buffer_i + Bs
                end do

            case('__S') ! south
                do l = 1, g
                    send_buff(buffer_i:buffer_i+Bs-1)   = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(Bs+g-l,g+1:Bs+g)
                    buffer_i                            = buffer_i + Bs
                end do

            case('__E') ! east
                do l = 1, g
                    send_buff(buffer_i:buffer_i+Bs-1)   = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g+1:Bs+g, Bs+g-l)
                    buffer_i                            = buffer_i + Bs
                end do

            case('_NE') ! northeast
                if ( blocks(my_block)%level == blocks(neighbor_block)%level ) then
                    ! blocks on same level
                    data_corner = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g+2:g+g+1, Bs:Bs+g-1)

                elseif ( blocks(my_block)%level - blocks(neighbor_block)%level == 1 ) then
                    ! sender on higher level
                    data_corner = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g+3:3*g+1:2, Bs-g:Bs+g-2:2)

                elseif ( blocks(my_block)%level - blocks(neighbor_block)%level == -1 ) then
                    ! sender on lower level
                    ! interpolate data
                    data_corner = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g+1:g+g, Bs+1:Bs+g)
                    call prediction_2D( data_corner , data_corner_fine)
                    data_corner = data_corner_fine(2:g+1, g-1:2*g-2)

                else
                    ! error case
                    print*, "error: mesh not graded, can not send/receive data"
                    stop
                end if
                ! send data
                do l = 1, g
                    send_buff(buffer_i:buffer_i+g-1)    = data_corner(l, 1:g)
                    buffer_i                            = buffer_i + g
                end do

            case('_NW') ! northwest
                if ( blocks(my_block)%level == blocks(neighbor_block)%level ) then
                    ! blocks on same level
                    data_corner = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g+2:g+g+1, g+2:g+g+1)

                elseif ( blocks(my_block)%level - blocks(neighbor_block)%level == 1 ) then
                    ! sender on higher level
                    data_corner = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g+3:3*g+1:2, g+3:3*g+1:2)

                elseif ( blocks(my_block)%level - blocks(neighbor_block)%level == -1 ) then
                    ! sender on lower level
                    ! interpolate data
                    data_corner = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g+1:g+g, g+1:g+g)
                    call prediction_2D( data_corner , data_corner_fine)
                    data_corner = data_corner_fine(2:g+1, 2:g+1)

                else
                    ! error case
                    print*, "error: mesh not graded, can not send/receive data"
                    stop

                end if
                ! send data
                do l = 1, g
                    send_buff(buffer_i:buffer_i+g-1)    = data_corner(l, 1:g)
                    buffer_i                            = buffer_i + g
                end do

            case('_SE') ! southeast
                if ( blocks(my_block)%level == blocks(neighbor_block)%level ) then
                    ! blocks on same level
                    data_corner = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(Bs:Bs+g-1, Bs:Bs+g-1)

                elseif ( blocks(my_block)%level - blocks(neighbor_block)%level == 1 ) then
                    ! sender on higher level
                    data_corner = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(Bs-g:Bs+g-2:2, Bs-g:Bs+g-2:2)

                elseif ( blocks(my_block)%level - blocks(neighbor_block)%level == -1 ) then
                    ! sender on lower level
                    ! interpolate data
                    data_corner = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(Bs+1:Bs+g, Bs+1:Bs+g)
                    call prediction_2D( data_corner , data_corner_fine)
                    data_corner = data_corner_fine(g-1:2*g-2, g-1:2*g-2)

                else
                    ! error case
                    print*, "error: mesh not graded, can not send/receive data"
                    stop

                end if
                ! send data
                do l = 1, g
                    send_buff(buffer_i:buffer_i+g-1)    = data_corner(l, 1:g)
                    buffer_i                            = buffer_i + g
                end do

            case('_SW') ! southwest
                if ( blocks(my_block)%level == blocks(neighbor_block)%level ) then
                    ! blocks on same level
                    data_corner = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(Bs:Bs+g-1, g+2:g+g+1)

                elseif ( blocks(my_block)%level - blocks(neighbor_block)%level == 1 ) then
                    ! sender on higher level
                    data_corner = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(Bs-g:Bs+g-2:2, g+3:3*g+1:2)

                elseif ( blocks(my_block)%level - blocks(neighbor_block)%level == -1 ) then
                    ! sender on lower level
                    ! interpolate data
                    data_corner = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(Bs+1:Bs+g, g+1:g+g)
                    call prediction_2D( data_corner , data_corner_fine)
                    data_corner = data_corner_fine(g-1:2*g-2, 2:g+1)

                else
                    ! error case
                    print*, "error: mesh not graded, can not send/receive data"
                    stop

                end if
                ! send data
                do l = 1, g
                    send_buff(buffer_i:buffer_i+g-1)    = data_corner(l, 1:g)
                    buffer_i                            = buffer_i + g
                end do

            case('NNE') ! northnortheast
                if ( blocks(my_block)%level - blocks(neighbor_block)%level == 1 ) then
                    ! sender on higher level
                    data_edge1_coarse = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g+3:3*g+1:2, g+1:Bs+g:2)
                    ! send data
                    do l = 1, g
                        send_buff(buffer_i:buffer_i+(Bs+1)/2-1)   = data_edge1_coarse(l, 1:(Bs+1)/2)
                        buffer_i                                  = buffer_i + (Bs+1)/2
                    end do

                elseif ( blocks(my_block)%level - blocks(neighbor_block)%level == -1 ) then
                    ! sender on lower level
                    data_edge1 = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g)
                    ! interpolate data
                    call prediction_2D( data_edge1 , data_edge1_fine)

                    ! send data
                    do l = 1, g
                        send_buff(buffer_i:buffer_i+Bs+g-1)    = data_edge1_fine(l+1, 1:Bs+g)
                        buffer_i                               = buffer_i + Bs+g
                    end do

                else
                    ! error case
                    print*, "error: mesh not graded, can not send/receive data"
                    stop

                end if

            case('SSE') ! southsouthheast
                if ( blocks(my_block)%level - blocks(neighbor_block)%level == 1 ) then
                    ! sender on higher level
                    data_edge1_coarse = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(Bs-g:Bs+g-2:2, g+1:Bs+g:2)
                    ! send data
                    do l = 1, g
                        send_buff(buffer_i:buffer_i+(Bs+1)/2-1)   = data_edge1_coarse(l, 1:(Bs+1)/2)
                        buffer_i                                  = buffer_i + (Bs+1)/2
                    end do

                elseif ( blocks(my_block)%level - blocks(neighbor_block)%level == -1 ) then
                    ! sender on lower level
                    data_edge1 = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_((Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g)
                    ! interpolate data
                    call prediction_2D( data_edge1 , data_edge1_fine)

                    ! send data
                    do l = 1, g
                        send_buff(buffer_i:buffer_i+Bs+g-1)    = data_edge1_fine(l+1, 1:Bs+g)
                        buffer_i                               = buffer_i + Bs+g
                    end do

                else
                    ! error case
                    print*, "error: mesh not graded, can not send/receive data"
                    stop

                end if

            case('NNW') ! northnorthwest
                if ( blocks(my_block)%level - blocks(neighbor_block)%level == 1 ) then
                    ! sender on higher level
                    data_edge1_coarse = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g+3:3*g+1:2, g+1:Bs+g:2)
                    ! send data
                    do l = 1, g
                        send_buff(buffer_i:buffer_i+(Bs+1)/2-1)   = data_edge1_coarse(l, 1:(Bs+1)/2)
                        buffer_i                                  = buffer_i + (Bs+1)/2
                    end do

                elseif ( blocks(my_block)%level - blocks(neighbor_block)%level == -1 ) then
                    ! sender on lower level
                    data_edge1 = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g)
                    ! interpolate data
                    call prediction_2D( data_edge1 , data_edge1_fine)

                    ! send data
                    do l = 1, g
                        send_buff(buffer_i:buffer_i+Bs+g-1)    = data_edge1_fine(l+1, 1:Bs+g)
                        buffer_i                               = buffer_i + Bs+g
                    end do

                else
                    ! error case
                    print*, "error: mesh not graded, can not send/receive data"
                    stop

                end if

            case('SSW') ! southsouthwest
                if ( blocks(my_block)%level - blocks(neighbor_block)%level == 1 ) then
                    ! sender on higher level
                    data_edge1_coarse = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(Bs-g:Bs+g-2:2, g+1:Bs+g:2)
                    ! send data
                    do l = 1, g
                        send_buff(buffer_i:buffer_i+(Bs+1)/2-1)   = data_edge1_coarse(l, 1:(Bs+1)/2)
                        buffer_i                                  = buffer_i + (Bs+1)/2
                    end do

                elseif ( blocks(my_block)%level - blocks(neighbor_block)%level == -1 ) then
                    ! sender on lower level
                    data_edge1 = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_((Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g)
                    ! interpolate data
                    call prediction_2D( data_edge1 , data_edge1_fine)

                    ! send data
                    do l = 1, g
                        send_buff(buffer_i:buffer_i+Bs+g-1)    = data_edge1_fine(l+1, 1:Bs+g)
                        buffer_i                               = buffer_i + Bs+g
                    end do

                else
                    ! error case
                    print*, "error: mesh not graded, can not send/receive data"
                    stop

                end if

            case('ENE') ! eastnortheast
                if ( blocks(my_block)%level - blocks(neighbor_block)%level == 1 ) then
                    ! sender on higher level
                    data_edge2_coarse = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g+1:Bs+g:2, Bs-g:Bs+g-2:2)
                    ! send data
                    do l = 1, g
                        send_buff(buffer_i:buffer_i+(Bs+1)/2-1)   = data_edge2_coarse(1:(Bs+1)/2, l)
                        buffer_i                                  = buffer_i + (Bs+1)/2
                    end do

                elseif ( blocks(my_block)%level - blocks(neighbor_block)%level == -1 ) then
                    ! sender on lower level
                    data_edge1 = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g)
                    ! interpolate data
                    call prediction_2D( data_edge1 , data_edge1_fine)

                    ! send data
                    do l = 1, g
                        send_buff(buffer_i:buffer_i+Bs+g-1)    = data_edge1_fine(1:Bs+g, Bs+l)
                        buffer_i                               = buffer_i + Bs+g
                    end do

                else
                    ! error case
                    print*, "error: mesh not graded, can not send/receive data"
                    stop

                end if

            case('ESE') ! eastsoutheast
                if ( blocks(my_block)%level - blocks(neighbor_block)%level == 1 ) then
                    ! sender on higher level
                    data_edge2_coarse = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g+1:Bs+g:2, Bs-g:Bs+g-2:2)
                    ! send data
                    do l = 1, g
                        send_buff(buffer_i:buffer_i+(Bs+1)/2-1)   = data_edge2_coarse(1:(Bs+1)/2, l)
                        buffer_i                                  = buffer_i + (Bs+1)/2
                    end do

                elseif ( blocks(my_block)%level - blocks(neighbor_block)%level == -1 ) then
                    ! sender on lower level
                    data_edge1 = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_((Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g)
                    ! interpolate data
                    call prediction_2D( data_edge1 , data_edge1_fine)

                    ! send data
                    do l = 1, g
                        send_buff(buffer_i:buffer_i+Bs+g-1)    = data_edge1_fine(1:Bs+g, Bs+l)
                        buffer_i                               = buffer_i + Bs+g
                    end do

                else
                    ! error case
                    print*, "error: mesh not graded, can not send/receive data"
                    stop

                end if

            case('WNW') ! westnorthwest
                if ( blocks(my_block)%level - blocks(neighbor_block)%level == 1 ) then
                    ! sender on higher level
                    data_edge2_coarse = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g+1:Bs+g:2, g+3:3*g+1:2)
                    ! send data
                    do l = 1, g
                        send_buff(buffer_i:buffer_i+(Bs+1)/2-1)   = data_edge2_coarse(1:(Bs+1)/2, l)
                        buffer_i                                  = buffer_i + (Bs+1)/2
                    end do

                elseif ( blocks(my_block)%level - blocks(neighbor_block)%level == -1 ) then
                    ! sender on lower level
                    data_edge1 = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g)
                    ! interpolate data
                    call prediction_2D( data_edge1 , data_edge1_fine)

                    ! send data
                    do l = 1, g
                        send_buff(buffer_i:buffer_i+Bs+g-1)    = data_edge1_fine(1:Bs+g, Bs+l)
                        buffer_i                               = buffer_i + Bs+g
                    end do

                else
                    ! error case
                    print*, "error: mesh not graded, can not send/receive data"
                    stop

                end if

            case('WSW') ! westsouthwest
                if ( blocks(my_block)%level - blocks(neighbor_block)%level == 1 ) then
                    ! sender on higher level
                    data_edge2_coarse = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g+1:Bs+g:2, g+3:3*g+1:2)
                    ! send data
                    do l = 1, g
                        send_buff(buffer_i:buffer_i+(Bs+1)/2-1)   = data_edge2_coarse(1:(Bs+1)/2, l)
                        buffer_i                                  = buffer_i + (Bs+1)/2
                    end do

                elseif ( blocks(my_block)%level - blocks(neighbor_block)%level == -1 ) then
                    ! sender on lower level
                    data_edge1 = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_((Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g)
                    ! interpolate data
                    call prediction_2D( data_edge1 , data_edge1_fine)

                    ! send data
                    do l = 1, g
                        send_buff(buffer_i:buffer_i+Bs+g-1)    = data_edge1_fine(1:Bs+g, Bs+l)
                        buffer_i                               = buffer_i + Bs+g
                    end do

                else
                    ! error case
                    print*, "error: mesh not graded, can not send/receive data"
                    stop

                end if

       end select

    end do

    if (plan_type==1) then
        ! internal communication
        recv_buff   = send_buff
        k_start     = 1
        k_end       = com_number
    else
        ! external communication
        ! send/receive data
        call MPI_Sendrecv( send_buff, 1000000, MPI_REAL4, my_dest, tag, recv_buff, 1000000, MPI_REAL4, my_dest, tag, MPI_COMM_WORLD, status, ierr)
        k_start = 1+com_number-k_shift
        k_end   = com_number+com_number-k_shift
    end if

    buffer_i  = 1

    ! write received data in block data
    do k = k_start, k_end

        my_block        = com_list( com_id+k-1, 5 )
        neighbor_block  = com_list( com_id+k-1, 4 )
        my_dir          = com_list( com_id+k-1, 7 )

        select case(blocks(my_block)%neighbor_dir(my_dir))
            case('__S') ! receive data for southern block boundary
                do l = 1, g
                    blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(Bs+g+l,g+1:Bs+g)         = recv_buff(buffer_i:buffer_i+Bs-1)
                    buffer_i                                                                                    = buffer_i + Bs
                end do

            case('__E') ! receive data for eastern block boundary
                do l = 1, g
                    blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g+1:Bs+g, Bs+g+l)        = recv_buff(buffer_i:buffer_i+Bs-1)
                    buffer_i                                                                                    = buffer_i + Bs
                end do

            case('__N') ! receive data for northern block boundary
                do l = 1, g
                    blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g-l+1,g+1:Bs+g)          = recv_buff(buffer_i:buffer_i+Bs-1)
                    buffer_i                                                                                    = buffer_i + Bs
                end do

            case('__W') ! receive data for western block boundary
                do l = 1, g
                    blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g+1:Bs+g, g-l+1)         = recv_buff(buffer_i:buffer_i+Bs-1)
                    buffer_i                                                                                    = buffer_i + Bs
                end do

            case('_SW') ! receive data for southwestern block boundary
                ! receive data
                do l = 1, g
                    data_corner(l, 1:g) = recv_buff(buffer_i:buffer_i+g-1)
                    buffer_i            = buffer_i + g
                end do
                ! write data
                blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(Bs+g+1:Bs+2*g, 1:g) = data_corner

            case('_SE') ! receive data for southeastern block boundary
                ! receive data
                do l = 1, g
                    data_corner(l, 1:g) = recv_buff(buffer_i:buffer_i+g-1)
                    buffer_i            = buffer_i + g
                end do
                ! write data
                blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(Bs+g+1:Bs+2*g, Bs+g+1:Bs+2*g) = data_corner

            case('_NW') ! receive data for northwestern block boundary
                ! receive data
                do l = 1, g
                    data_corner(l, 1:g) = recv_buff(buffer_i:buffer_i+g-1)
                    buffer_i            = buffer_i + g
                end do
                ! write data
                blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(1:g, 1:g) = data_corner

            case('_NE') ! receive data for northeastern block boundary
                ! receive data
                do l = 1, g
                    data_corner(l, 1:g) = recv_buff(buffer_i:buffer_i+g-1)
                    buffer_i            = buffer_i + g
                end do
                ! write data
                blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(1:g, Bs+g+1:Bs+2*g) = data_corner

            case('SSE') ! receive data for southsoutheastern boundary
                if ( blocks(my_block)%level < blocks(neighbor_block)%level ) then
                    ! sender on higher level
                    ! receive data
                    do l = 1, g
                        data_edge1_coarse(l, 1:(Bs+1)/2)   = recv_buff(buffer_i:buffer_i+(Bs+1)/2-1)
                        buffer_i                           = buffer_i + (Bs+1)/2
                    end do
                    ! write data
                    blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(Bs+g+1:Bs+2*g, g+(Bs+1)/2:Bs+g) = data_edge1_coarse
                else
                    ! sender on lower level
                    ! receive data
                    do l = 1, g
                        data_edge1_fine(l, 1:Bs+g)         = recv_buff(buffer_i:buffer_i+Bs+g-1)
                        buffer_i                           = buffer_i + Bs+g
                    end do
                    ! write data
                    blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(Bs+g+1:Bs+2*g, 1:Bs+g) = data_edge1_fine(1:g, :)

                end if

            case('NNE') ! receive data for northnortheastern boundary
                if ( blocks(my_block)%level < blocks(neighbor_block)%level ) then
                    ! sender on higher level
                    ! receive data
                    do l = 1, g
                        data_edge1_coarse(l, 1:(Bs+1)/2)   = recv_buff(buffer_i:buffer_i+(Bs+1)/2-1)
                        buffer_i                           = buffer_i + (Bs+1)/2
                    end do
                    ! write data
                    blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(1:g, g+(Bs+1)/2:Bs+g) = data_edge1_coarse
                else
                    ! sender on lower level
                    ! receive data
                    do l = 1, g
                        data_edge1_fine(l, 1:Bs+g)         = recv_buff(buffer_i:buffer_i+Bs+g-1)
                        buffer_i                           = buffer_i + Bs+g
                    end do
                    ! write data
                    blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(1:g, 1:Bs+g) = data_edge1_fine(1:g, :)

                end if

            case('SSW') ! receive data for southsouthwestern boundary
                if ( blocks(my_block)%level < blocks(neighbor_block)%level ) then
                    ! sender on higher level
                    ! receive data
                    do l = 1, g
                        data_edge1_coarse(l, 1:(Bs+1)/2)   = recv_buff(buffer_i:buffer_i+(Bs+1)/2-1)
                        buffer_i                           = buffer_i + (Bs+1)/2
                    end do
                    ! write data
                    blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(Bs+g+1:Bs+2*g, g+1:g+(Bs+1)/2) = data_edge1_coarse
                else
                    ! sender on lower level
                    ! receive data
                    do l = 1, g
                        data_edge1_fine(l, 1:Bs+g)         = recv_buff(buffer_i:buffer_i+Bs+g-1)
                        buffer_i                           = buffer_i + Bs+g
                    end do
                    ! write data
                    blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(Bs+g+1:Bs+2*g, g+1:Bs+2*g) = data_edge1_fine(1:g, :)

                end if

            case('NNW') ! receive data for northnorthwestern boundary
                if ( blocks(my_block)%level < blocks(neighbor_block)%level ) then
                    ! sender on higher level
                    ! receive data
                    do l = 1, g
                        data_edge1_coarse(l, 1:(Bs+1)/2)   = recv_buff(buffer_i:buffer_i+(Bs+1)/2-1)
                        buffer_i                           = buffer_i + (Bs+1)/2
                    end do
                    ! write data
                    blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(1:g, g+1:g+(Bs+1)/2) = data_edge1_coarse
                else
                    ! sender on lower level
                    ! receive data
                    do l = 1, g
                        data_edge1_fine(l, 1:Bs+g)         = recv_buff(buffer_i:buffer_i+Bs+g-1)
                        buffer_i                           = buffer_i + Bs+g
                    end do
                    ! write data
                    blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(1:g, g+1:Bs+2*g) = data_edge1_fine(1:g, :)

                end if

            case('WNW') ! receive data for westnorthwest boundary
                if ( blocks(my_block)%level < blocks(neighbor_block)%level ) then
                    ! sender on higher level
                    ! receive data
                    do l = 1, g
                        data_edge2_coarse(1:(Bs+1)/2, l)   = recv_buff(buffer_i:buffer_i+(Bs+1)/2-1)
                        buffer_i                           = buffer_i + (Bs+1)/2
                    end do
                    ! write data
                    blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g+1:g+(Bs+1)/2, 1:g) = data_edge2_coarse
                else
                    ! sender on lower level
                    ! receive data
                    do l = 1, g
                        data_edge1_fine(1:Bs+g, l)         = recv_buff(buffer_i:buffer_i+Bs+g-1)
                        buffer_i                           = buffer_i + Bs+g
                    end do
                    ! write data
                    blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g+1:Bs+2*g, 1:g) = data_edge1_fine(:, 1:g)

                end if

            case('WSW') ! receive data for westsouthwest boundary
                if ( blocks(my_block)%level < blocks(neighbor_block)%level ) then
                    ! sender on higher level
                    ! receive data
                    do l = 1, g
                        data_edge2_coarse(1:(Bs+1)/2, l)   = recv_buff(buffer_i:buffer_i+(Bs+1)/2-1)
                        buffer_i                           = buffer_i + (Bs+1)/2
                    end do
                    ! write data
                    blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g+(Bs+1)/2:Bs+g, 1:g) = data_edge2_coarse
                else
                    ! sender on lower level
                    ! receive data
                    do l = 1, g
                        data_edge1_fine(1:Bs+g, l)         = recv_buff(buffer_i:buffer_i+Bs+g-1)
                        buffer_i                           = buffer_i + Bs+g
                    end do
                    ! write data
                    blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(1:Bs+g, 1:g) = data_edge1_fine(:, 1:g)

                end if

            case('ENE') ! receive data for eastnortheast boundary
                if ( blocks(my_block)%level < blocks(neighbor_block)%level ) then
                    ! sender on higher level
                    ! receive data
                    do l = 1, g
                        data_edge2_coarse(1:(Bs+1)/2, l)   = recv_buff(buffer_i:buffer_i+(Bs+1)/2-1)
                        buffer_i                           = buffer_i + (Bs+1)/2
                    end do
                    ! write data
                    blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g+1:g+(Bs+1)/2, Bs+g+1:Bs+2*g) = data_edge2_coarse
                else
                    ! sender on lower level
                    ! receive data
                    do l = 1, g
                        data_edge1_fine(1:Bs+g, l)         = recv_buff(buffer_i:buffer_i+Bs+g-1)
                        buffer_i                           = buffer_i + Bs+g
                    end do
                    ! write data
                    blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g+1:Bs+2*g, Bs+g+1:Bs+2*g) = data_edge1_fine(:, 1:g)

                end if

            case('ESE') ! receive data for eastsoutheast boundary
                if ( blocks(my_block)%level < blocks(neighbor_block)%level ) then
                    ! sender on higher level
                    ! receive data
                    do l = 1, g
                        data_edge2_coarse(1:(Bs+1)/2, l)   = recv_buff(buffer_i:buffer_i+(Bs+1)/2-1)
                        buffer_i                           = buffer_i + (Bs+1)/2
                    end do
                    ! write data
                    blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g+(Bs+1)/2:Bs+g, Bs+g+1:Bs+2*g) = data_edge2_coarse
                else
                    ! sender on lower level
                    ! receive data
                    do l = 1, g
                        data_edge1_fine(1:Bs+g, l)         = recv_buff(buffer_i:buffer_i+Bs+g-1)
                        buffer_i                           = buffer_i + Bs+g
                    end do
                    ! write data
                    blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(1:Bs+g, Bs+g+1:Bs+2*g) = data_edge1_fine(:, 1:g)

                end if

        end select

    end do

    deallocate( data_corner, stat=allocate_error )
    deallocate( data_corner_fine, stat=allocate_error )
    deallocate( data_edge1, stat=allocate_error )
    deallocate( data_edge1_coarse, stat=allocate_error )
    deallocate( data_edge2_coarse, stat=allocate_error )
    deallocate( data_edge1_fine, stat=allocate_error )

end subroutine send_receive_data
