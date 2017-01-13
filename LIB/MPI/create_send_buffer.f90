! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: create_send_buffer.f90
! version: 0.4
! author: msr
!
! write send buffer with given com list
!
! input:    - heavy data array
!           - params struct
!           - communications list
!           - number of communications to send/receive
! output:   - send buffer
!
! -------------------------------------------------------------------------------------------------------------------------
! dirs = (/'__N', '__E', '__S', '__W', '_NE', '_NW', '_SE', '_SW', 'NNE', 'NNW', 'SSE', 'SSW', 'ENE', 'ESE', 'WNW', 'WSW'/)
! -------------------------------------------------------------------------------------------------------------------------
!
! = log ======================================================================================
!
! 06/01/17 - create for v0.4
! ********************************************************************************************

subroutine create_send_buffer(params, hvy_block, com_list, com_number, send_buff, buffer_i)

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(in)                  :: params
    ! heavy data array - block data
    real(kind=rk), intent(in)                       :: hvy_block(:, :, :, :)

    ! com list
    integer(kind=ik), intent(in)                    :: com_list(:, :)
    ! com_list id, number of communications
    integer(kind=ik), intent(in)                    :: com_number

    ! send buffer
    real(kind=rk), intent(out)                      :: send_buff(:)

    ! buffer index
    integer(kind=ik), intent(out)                   :: buffer_i

    ! grid parameter
    integer(kind=ik)                                :: Bs, g

    ! interpolation variables
    real(kind=rk), dimension(:,:), allocatable      :: data_corner, data_corner_fine, data_edge, data_edge_fine

    ! allocation error variable
    integer(kind=ik)                                :: allocate_error

    ! com list elements
    integer(kind=ik)                                :: my_block, neighbor_block, my_dir, level_diff

    ! loop variable
    integer(kind=ik)                                :: k, l, dF

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! grid parameter
    Bs    = params%number_block_nodes
    g     = params%number_ghost_nodes

    allocate( data_corner( g, g), stat=allocate_error )
    allocate( data_corner_fine( 2*g-1, 2*g-1), stat=allocate_error )
    allocate( data_edge( (Bs+1)/2 + g/2, (Bs+1)/2 + g/2), stat=allocate_error )
    allocate( data_edge_fine( Bs+g, Bs+g), stat=allocate_error )

    buffer_i         = 1

    data_corner      = 9.0e9_rk
    data_corner_fine = 9.0e9_rk
    data_edge        = 9.0e9_rk
    data_edge_fine   = 9.0e9_rk

    send_buff        = 9.0e9_rk

!---------------------------------------------------------------------------------------------
! main body

    ! fill send buffer
    do k = 1 , com_number

        my_block        = com_list( k, 3 )
        neighbor_block  = com_list( k, 4 )
        my_dir          = com_list( k, 5 )
        level_diff      = com_list( k, 6 )

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

    ! decrease buffer index (for using in other subroutines)
    buffer_i = buffer_i - 1

    ! clean up
    deallocate( data_corner, stat=allocate_error )
    deallocate( data_corner_fine, stat=allocate_error )
    deallocate( data_edge, stat=allocate_error )
    deallocate( data_edge_fine, stat=allocate_error )

end subroutine create_send_buffer
