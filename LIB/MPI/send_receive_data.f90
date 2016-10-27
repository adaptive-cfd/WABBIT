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

subroutine send_receive_data(com_id, plan_type, com_number, com_list)

    use mpi
    use module_blocks

    implicit none

    real(kind=rk) , dimension(10000)        :: send_buff, recv_buff

    integer, intent(in)                     :: com_id, com_number, plan_type
    integer, dimension(200, 7), intent(in)  :: com_list
    integer                                 :: Bs, g, k, l, buffer_i, dF
    integer                                 :: my_block, my_dir, tag, ierr, rank, k_shift, my_dest
    integer                                 :: status(MPI_status_size)

    tag      = 0

    send_buff = 9.0e9_rk
    recv_buff = 9.0e9_rk
    buffer_i  = 1

    g  = blocks_params%number_ghost_nodes
    Bs = blocks_params%size_block
    dF = 1

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

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
    do k = 1+k_shift, com_number+k_shift

        my_block        = com_list( com_id+k-1, 4 )
        my_dir          = com_list( com_id+k-1, 6 )

        select case(blocks(my_block)%neighbor_dir(my_dir))
            case('NO') ! north
                do l = 1, g
                    send_buff(buffer_i:buffer_i+Bs-1)   = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g+l+1,g+1:Bs+g)
                    buffer_i                            = buffer_i + Bs
                end do
            case('WE') ! west
                do l = 1, g
                    send_buff(buffer_i:buffer_i+Bs-1)   = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g+1:Bs+g, g+l+1)
                    buffer_i                            = buffer_i + Bs
                end do
            case('SO') ! south
                do l = 1, g
                    send_buff(buffer_i:buffer_i+Bs-1)   = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(Bs+g-l,g+1:Bs+g)
                    buffer_i                            = buffer_i + Bs
                end do
            case('EA') ! east
                do l = 1, g
                    send_buff(buffer_i:buffer_i+Bs-1)   = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g+1:Bs+g, Bs+g-l)
                    buffer_i                            = buffer_i + Bs
                end do
            case('NE') ! northeast
                do l = 1, g
                    send_buff(buffer_i:buffer_i+g-1)    = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g+l+1, Bs:Bs+g-1)
                    buffer_i                            = buffer_i + g
                end do
            case('NW') ! northwest
                do l = 1, g
                    send_buff(buffer_i:buffer_i+g-1)    = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g+l+1, g+2:g+1+g)
                    buffer_i                            = buffer_i + g
                end do
            case('SE') ! southeast
                do l = 1, g
                    send_buff(buffer_i:buffer_i+g-1)    = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(Bs+g-l, Bs:Bs+g-1)
                    buffer_i                            = buffer_i + g
                end do
            case('SW') ! southwest
                do l = 1, g
                    send_buff(buffer_i:buffer_i+g-1)    = blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(Bs+g-l, g+2:g+1+g)
                    buffer_i                            = buffer_i + g
                end do
       end select

    end do

    if (plan_type==1) then
        ! internal communication
        recv_buff = send_buff
    else
        ! external communication
        ! send/receive data
        call MPI_Sendrecv( send_buff, 10000, MPI_REAL4, my_dest, tag, recv_buff, 10000, MPI_REAL4, my_dest, tag, MPI_COMM_WORLD, status, ierr)
    end if

    buffer_i  = 1
    ! write received data in block data
    do k = 1+com_number-k_shift, com_number+com_number-k_shift

        my_block        = com_list( com_id+k-1, 5 )
        my_dir          = com_list( com_id+k-1, 7 )

        select case(blocks(my_block)%neighbor_dir(my_dir))
            case('SO') ! receive data for southern block boundary
                do l = 1, g
                    blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(Bs+g+l,g+1:Bs+g)         = recv_buff(buffer_i:buffer_i+Bs-1)
                    buffer_i                                                                                    = buffer_i + Bs
                end do
            case('EA') ! receive data for eastern block boundary
                do l = 1, g
                    blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g+1:Bs+g, Bs+g+l)        = recv_buff(buffer_i:buffer_i+Bs-1)
                    buffer_i                                                                                    = buffer_i + Bs
                end do
            case('NO') ! receive data for northern block boundary
                do l = 1, g
                    blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g-l+1,g+1:Bs+g)          = recv_buff(buffer_i:buffer_i+Bs-1)
                    buffer_i                                                                                    = buffer_i + Bs
                end do
            case('WE') ! receive data for western block boundary
                do l = 1, g
                    blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g+1:Bs+g, g-l+1)         = recv_buff(buffer_i:buffer_i+Bs-1)
                    buffer_i                                                                                    = buffer_i + Bs
                end do
            case('SW') ! receive data for southwestern block boundary
                do l = 1, g
                    blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(Bs+g+l, 1:g)             = recv_buff(buffer_i:buffer_i+g-1)
                    buffer_i                                                                                    = buffer_i + g
                end do
            case('SE') ! receive data for southeastern block boundary
                do l = 1, g
                    blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(Bs+g+l, Bs+g+1:Bs+g+g)   = recv_buff(buffer_i:buffer_i+g-1)
                    buffer_i                                                                                    = buffer_i + g
                end do
            case('NW') ! receive data for northwestern block boundary
                do l = 1, g
                    blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g-l+1, 1:g)              = recv_buff(buffer_i:buffer_i+g-1)
                    buffer_i                                                                                    = buffer_i + g
                end do
            case('NE') ! receive data for northeastern block boundary
                do l = 1, g
                    blocks_data( blocks(my_block)%proc_data_id )%data_fields(dF)%data_(g-l+1, Bs+g+1:Bs+g+g)    = recv_buff(buffer_i:buffer_i+g-1)
                    buffer_i                                                                                    = buffer_i + g
                end do
        end select

    end do

end subroutine send_receive_data
