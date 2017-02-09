! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: create_send_buffer_3D.f90
! version: 0.5
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
! --------------------------------------------------------------------------------------------
! neighbor codes:
! ---------------
! for imagination:  - 6-sided dice with '1'-side on top, '6'-side on bottom, '2'-side in front
!                   - edge: boundary between two sides - use sides numbers for coding
!                   - corner: between three sides - so use all three sides numbers
!                   - block on higher/lower level: block shares face/edge and one unique corner,
!                     so use this corner code in second part of neighbor code
!
! faces:  '__1/___', '__2/___', '__3/___', '__4/___', '__5/___', '__6/___'
! edges:  '_12/___', '_13/___', '_14/___', '_15/___'
!         '_62/___', '_63/___', '_64/___', '_65/___'
!         '_23/___', '_25/___', '_43/___', '_45/___'
! corner: '123/___', '134/___', '145/___', '152/___'
!         '623/___', '634/___', '645/___', '652/___'
!
! complete neighbor code array, 74 possible neighbor relations
! neighbors = (/'__1/___', '__2/___', '__3/___', '__4/___', '__5/___', '__6/___', '_12/___', '_13/___', '_14/___', '_15/___',
!               '_62/___', '_63/___', '_64/___', '_65/___', '_23/___', '_25/___', '_43/___', '_45/___', '123/___', '134/___',
!               '145/___', '152/___', '623/___', '634/___', '645/___', '652/___', '__1/123', '__1/134', '__1/145', '__1/152',
!               '__2/123', '__2/623', '__2/152', '__2/652', '__3/123', '__3/623', '__3/134', '__3/634', '__4/134', '__4/634',
!               '__4/145', '__4/645', '__5/145', '__5/645', '__5/152', '__5/652', '__6/623', '__6/634', '__6/645', '__6/652',
!               '_12/123', '_12/152', '_13/123', '_13/134', '_14/134', '_14/145', '_15/145', '_15/152', '_62/623', '_62/652',
!               '_63/623', '_63/634', '_64/634', '_64/645', '_65/645', '_65/652', '_23/123', '_23/623', '_25/152', '_25/652',
!               '_43/134', '_43/634', '_45/145', '_45/645' /)
! --------------------------------------------------------------------------------------------
!
! = log ======================================================================================
!
! 01/02/17 - create
!
! ********************************************************************************************

subroutine create_send_buffer_3D(params, hvy_block, com_list, com_number, send_buff, buffer_i)

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(in)                  :: params
    ! heavy data array - block data
    real(kind=rk), intent(in)                       :: hvy_block(:, :, :, :, :)

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
    real(kind=rk), dimension(:,:,:), allocatable    :: data_corner, data_corner_fine, data_face, data_face_fine

    ! allocation error variable
    integer(kind=ik)                                :: allocate_error

    ! com list elements
    integer(kind=ik)                                :: my_block, neighbor_block, my_dir, level_diff

    ! loop variable
    integer(kind=ik)                                :: k, i, j, dF

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! grid parameter
    Bs    = params%number_block_nodes
    g     = params%number_ghost_nodes

    allocate( data_corner( g, g, g), stat=allocate_error )
    call check_allocation(allocate_error)
    allocate( data_corner_fine( 2*g-1, 2*g-1, 2*g-1), stat=allocate_error )
    call check_allocation(allocate_error)

    allocate( data_face( (Bs+1)/2 + g/2, (Bs+1)/2 + g/2, (Bs+1)/2 + g/2), stat=allocate_error )
    call check_allocation(allocate_error)
    allocate( data_face_fine( Bs+g, Bs+g, Bs+g), stat=allocate_error )
    call check_allocation(allocate_error)

    buffer_i         = 1

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

            ! '__1/___'
            case(1)
                do dF = 2, params%number_data_fields+1
                    do i = 1, g
                        do j = 1, Bs
                            send_buff(buffer_i:buffer_i+Bs-1)   = hvy_block( g+j, g+1:Bs+g, Bs+g-i, dF, my_block )
                            buffer_i                            = buffer_i + Bs
                        end do
                    end do
                end do

            ! '__2/___'
            case(2)
                do dF = 2, params%number_data_fields+1
                    do i = 1, g
                        do j = 1, Bs
                            send_buff(buffer_i:buffer_i+Bs-1)   = hvy_block( g+j, g+1+i, g+1:Bs+g, dF, my_block )
                            buffer_i                            = buffer_i + Bs
                        end do
                    end do
                end do

            ! '__3/___'
            case(3)
                do dF = 2, params%number_data_fields+1
                    do i = 1, g
                        do j = 1, Bs
                            send_buff(buffer_i:buffer_i+Bs-1)   = hvy_block( Bs+g-i, g+j, g+1:Bs+g, dF, my_block )
                            buffer_i                            = buffer_i + Bs
                        end do
                    end do
                end do

            ! '__4/___'
            case(4)
                do dF = 2, params%number_data_fields+1
                    do i = 1, g
                        do j = 1, Bs
                            send_buff(buffer_i:buffer_i+Bs-1)   = hvy_block( g+j, Bs+g-i, g+1:Bs+g, dF, my_block )
                            buffer_i                            = buffer_i + Bs
                        end do
                    end do
                end do

            ! '__5/___'
            case(5)
                do dF = 2, params%number_data_fields+1
                    do i = 1, g
                        do j = 1, Bs
                            send_buff(buffer_i:buffer_i+Bs-1)   = hvy_block( g+1+i, g+j, g+1:Bs+g, dF, my_block )
                            buffer_i                            = buffer_i + Bs
                        end do
                    end do
                end do

            ! '__6/___'
            case(6)
                do dF = 2, params%number_data_fields+1
                    do i = 1, g
                        do j = 1, Bs
                            send_buff(buffer_i:buffer_i+Bs-1)   = hvy_block( g+j, g+1:Bs+g, g+1+i, dF, my_block )
                            buffer_i                            = buffer_i + Bs
                        end do
                    end do
                end do

            ! '_12/___'
            case(7)
                do dF = 2, params%number_data_fields+1
                    do i = 1, g
                        do j = 1, g
                            send_buff(buffer_i:buffer_i+Bs-1)   = hvy_block( g+1:Bs+g, g+1+i, Bs+g-j, dF, my_block )
                            buffer_i                            = buffer_i + Bs
                        end do
                    end do
                end do

            ! '_13/___'
            case(8)
                do dF = 2, params%number_data_fields+1
                    do i = 1, g
                        do j = 1, g
                            send_buff(buffer_i:buffer_i+Bs-1)   = hvy_block( Bs+g-i, g+1:Bs+g, Bs+g-j, dF, my_block )
                            buffer_i                            = buffer_i + Bs
                        end do
                    end do
                end do

            ! '_14/___'
            case(9)
                do dF = 2, params%number_data_fields+1
                    do i = 1, g
                        do j = 1, g
                            send_buff(buffer_i:buffer_i+Bs-1)   = hvy_block( g+1:Bs+g, Bs+g-i, Bs+g-j, dF, my_block )
                            buffer_i                            = buffer_i + Bs
                        end do
                    end do
                end do

            ! '_15/___'
            case(10)
                do dF = 2, params%number_data_fields+1
                    do i = 1, g
                        do j = 1, g
                            send_buff(buffer_i:buffer_i+Bs-1)   = hvy_block( g+1+i, g+1:Bs+g, Bs+g-j, dF, my_block )
                            buffer_i                            = buffer_i + Bs
                        end do
                    end do
                end do

            ! '_62/___'
            case(11)
                do dF = 2, params%number_data_fields+1
                    do i = 1, g
                        do j = 1, g
                            send_buff(buffer_i:buffer_i+Bs-1)   = hvy_block( g+1:Bs+g, g+1+i, g+1+j, dF, my_block )
                            buffer_i                            = buffer_i + Bs
                        end do
                    end do
                end do

            ! '_63/___'
            case(12)
                do dF = 2, params%number_data_fields+1
                    do i = 1, g
                        do j = 1, g
                            send_buff(buffer_i:buffer_i+Bs-1)   = hvy_block( Bs+g-i, g+1:Bs+g, g+1+j, dF, my_block )
                            buffer_i                            = buffer_i + Bs
                        end do
                    end do
                end do

            ! '_64/___'
            case(13)
                do dF = 2, params%number_data_fields+1
                    do i = 1, g
                        do j = 1, g
                            send_buff(buffer_i:buffer_i+Bs-1)   = hvy_block( g+1:Bs+g, Bs+g-i, g+1+j, dF, my_block )
                            buffer_i                            = buffer_i + Bs
                        end do
                    end do
                end do

            ! '_65/___'
            case(14)
                do dF = 2, params%number_data_fields+1
                    do i = 1, g
                        do j = 1, g
                            send_buff(buffer_i:buffer_i+Bs-1)   = hvy_block( g+1+i, g+1:Bs+g, g+1+j, dF, my_block )
                            buffer_i                            = buffer_i + Bs
                        end do
                    end do
                end do

            ! '_23/___'
            case(15)
                do dF = 2, params%number_data_fields+1
                    do i = 1, g
                        do j = 1, g
                            send_buff(buffer_i:buffer_i+Bs-1)   = hvy_block( Bs+g-i, g+1+j, g+1:Bs+g, dF, my_block )
                            buffer_i                            = buffer_i + Bs
                        end do
                    end do
                end do

            ! '_25/___'
            case(16)
                do dF = 2, params%number_data_fields+1
                    do i = 1, g
                        do j = 1, g
                            send_buff(buffer_i:buffer_i+Bs-1)   = hvy_block( g+1+i, g+1+j, g+1:Bs+g, dF, my_block )
                            buffer_i                            = buffer_i + Bs
                        end do
                    end do
                end do

            ! '_43/___'
            case(17)
                do dF = 2, params%number_data_fields+1
                    do i = 1, g
                        do j = 1, g
                            send_buff(buffer_i:buffer_i+Bs-1)   = hvy_block( Bs+g-i, Bs+g-j, g+1:Bs+g, dF, my_block )
                            buffer_i                            = buffer_i + Bs
                        end do
                    end do
                end do

            ! '_45/___'
            case(18)
                do dF = 2, params%number_data_fields+1
                    do i = 1, g
                        do j = 1, g
                            send_buff(buffer_i:buffer_i+Bs-1)   = hvy_block( g+1+i, Bs+g-j, g+1:Bs+g, dF, my_block )
                            buffer_i                            = buffer_i + Bs
                        end do
                    end do
                end do

            ! '123/___'
            case(19)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == 0 ) then
                        ! blocks on same level
                        data_corner = hvy_block( Bs:Bs-1+g, g+2:g+1+g, Bs:Bs-1+g, dF, my_block )

                    elseif ( level_diff == -1 ) then
                        ! sender one level down
                        ! interpolate data
                        ! data to refine
                        data_corner = hvy_block( Bs+1:Bs+g, g+1:g+g, Bs+1:Bs+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_corner , data_corner_fine, params%order_predictor)
                        ! data to synchronize
                        data_corner = data_corner_fine(g-1:2*g-2, 2:g+1, g-1:2*g-2)

                    elseif ( level_diff == 1) then
                        ! sender one level up
                        data_corner = hvy_block( Bs-g:Bs-2+g:2, g+3:g+1+g+g:2, Bs-g:Bs-2+g:2, dF, my_block )


                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                    ! write buffer
                    do i = 1, g
                        do j = 1, g
                            send_buff(buffer_i:buffer_i+g-1)    = data_corner( i, j, 1:g )
                            buffer_i                            = buffer_i + g
                        end do
                    end do

                end do

            ! '134/___'
            case(20)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == 0 ) then
                        ! blocks on same level
                        data_corner = hvy_block( Bs:Bs-1+g, Bs:Bs-1+g, Bs:Bs-1+g, dF, my_block )

                    elseif ( level_diff == -1 ) then
                        ! sender one level down
                        ! interpolate data
                        ! data to refine
                        data_corner = hvy_block( Bs:Bs-1+g, Bs:Bs-1+g, Bs:Bs-1+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_corner , data_corner_fine, params%order_predictor)
                        ! data to synchronize
                        data_corner = data_corner_fine(g-1:2*g-2, g-1:2*g-2, g-1:2*g-2)

                    elseif ( level_diff == 1) then
                        ! sender one level up
                        data_corner = hvy_block( Bs-g:Bs-2+g:2, Bs-g:Bs-2+g:2, Bs-g:Bs-2+g:2, dF, my_block )

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                    ! write buffer
                    do i = 1, g
                        do j = 1, g
                            send_buff(buffer_i:buffer_i+g-1)    = data_corner( i, j, 1:g )
                            buffer_i                            = buffer_i + g
                        end do
                    end do

                end do

            ! '145/___'
            case(21)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == 0 ) then
                        ! blocks on same level
                        data_corner = hvy_block( g+2:g+1+g, Bs:Bs-1+g, Bs:Bs-1+g, dF, my_block )

                    elseif ( level_diff == -1 ) then
                        ! sender one level down
                        ! interpolate data
                        ! data to refine
                        data_corner = hvy_block( g+1:g+g, Bs+1:Bs+g, Bs+1:Bs+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_corner , data_corner_fine, params%order_predictor)
                        ! data to synchronize
                        data_corner = data_corner_fine(2:g+1, g-1:2*g-2, g-1:2*g-2)

                    elseif ( level_diff == 1) then
                        ! sender one level up
                        data_corner = hvy_block( g+3:g+1+g+g:2, Bs-g:Bs-2+g:2, Bs-g:Bs-2+g:2, dF, my_block )

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                    ! write buffer
                    do i = 1, g
                        do j = 1, g
                            send_buff(buffer_i:buffer_i+g-1)    = data_corner( i, j, 1:g )
                            buffer_i                            = buffer_i + g
                        end do
                    end do

                end do

            ! '152/___'
            case(22)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == 0 ) then
                        ! blocks on same level
                        data_corner = hvy_block( g+2:g+1+g, g+2:g+1+g, Bs:Bs-1+g, dF, my_block )

                    elseif ( level_diff == -1 ) then
                        ! sender one level down
                        ! interpolate data
                        ! data to refine
                        data_corner = hvy_block( g+1:g+g, g+1:g+g, Bs+1:Bs+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_corner , data_corner_fine, params%order_predictor)
                        ! data to synchronize
                        data_corner = data_corner_fine(2:g+1, 2:g+1, g-1:2*g-2)

                    elseif ( level_diff == 1) then
                        ! sender one level up
                        data_corner = hvy_block( g+3:g+1+g+g:2, g+3:g+1+g+g:2, Bs-g:Bs-2+g:2, dF, my_block )

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                    ! write buffer
                    do i = 1, g
                        do j = 1, g
                            send_buff(buffer_i:buffer_i+g-1)    = data_corner( i, j, 1:g )
                            buffer_i                            = buffer_i + g
                        end do
                    end do

                end do

            ! '623/___'
            case(23)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == 0 ) then
                        ! blocks on same level
                        data_corner = hvy_block( Bs:Bs-1+g, g+2:g+1+g, g+2:g+1+g, dF, my_block )

                    elseif ( level_diff == -1 ) then
                        ! sender one level down
                        ! interpolate data
                        ! data to refine
                        data_corner = hvy_block( Bs+1:Bs+g, g+1:g+g, g+1:g+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_corner , data_corner_fine, params%order_predictor)
                        ! data to synchronize
                        data_corner = data_corner_fine(g-1:2*g-2, 2:g+1, 2:g+1)

                    elseif ( level_diff == 1) then
                        ! sender one level up
                        data_corner = hvy_block( Bs-g:Bs-2+g:2, g+3:g+1+g+g:2, g+3:g+1+g+g:2, dF, my_block )

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                    ! write buffer
                    do i = 1, g
                        do j = 1, g
                            send_buff(buffer_i:buffer_i+g-1)    = data_corner( i, j, 1:g )
                            buffer_i                            = buffer_i + g
                        end do
                    end do

                end do

            ! '634/___'
            case(24)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == 0 ) then
                        ! blocks on same level
                        data_corner = hvy_block( Bs:Bs-1+g, Bs:Bs-1+g, g+2:g+1+g, dF, my_block )

                    elseif ( level_diff == -1 ) then
                        ! sender one level down
                        ! interpolate data
                        ! data to refine
                        data_corner = hvy_block( Bs:Bs-1+g, Bs:Bs-1+g, g+1:g+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_corner , data_corner_fine, params%order_predictor)
                        ! data to synchronize
                        data_corner = data_corner_fine(g-1:2*g-2, g-1:2*g-2, 2:g+1)

                    elseif ( level_diff == 1) then
                        ! sender one level up
                        data_corner = hvy_block( Bs-g:Bs-2+g:2, Bs-g:Bs-2+g:2, g+3:g+1+g+g:2, dF, my_block )

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                    ! write buffer
                    do i = 1, g
                        do j = 1, g
                            send_buff(buffer_i:buffer_i+g-1)    = data_corner( i, j, 1:g )
                            buffer_i                            = buffer_i + g
                        end do
                    end do

                end do

            ! '645/___'
            case(25)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == 0 ) then
                        ! blocks on same level
                        data_corner = hvy_block( g+2:g+1+g, Bs:Bs-1+g, g+2:g+1+g, dF, my_block )

                    elseif ( level_diff == -1 ) then
                        ! sender one level down
                        ! interpolate data
                        ! data to refine
                        data_corner = hvy_block( g+1:g+g, Bs+1:Bs+g, g+1:g+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_corner , data_corner_fine, params%order_predictor)
                        ! data to synchronize
                        data_corner = data_corner_fine(2:g+1, g-1:2*g-2, 2:g+1)

                    elseif ( level_diff == 1) then
                        ! sender one level up
                        data_corner = hvy_block( g+3:g+1+g+g:2, Bs-g:Bs-2+g:2, g+3:g+1+g+g:2, dF, my_block )

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                    ! write buffer
                    do i = 1, g
                        do j = 1, g
                            send_buff(buffer_i:buffer_i+g-1)    = data_corner( i, j, 1:g )
                            buffer_i                            = buffer_i + g
                        end do
                    end do

                end do

            ! '652/___'
            case(26)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == 0 ) then
                        ! blocks on same level
                        data_corner = hvy_block( g+2:g+1+g, g+2:g+1+g, g+2:g+1+g, dF, my_block )

                    elseif ( level_diff == -1 ) then
                        ! sender one level down
                        ! interpolate data
                        ! data to refine
                        data_corner = hvy_block( g+1:g+g, g+1:g+g, g+1:g+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_corner , data_corner_fine, params%order_predictor)
                        ! data to synchronize
                        data_corner = data_corner_fine(2:g+1, 2:g+1, 2:g+1)

                    elseif ( level_diff == 1) then
                        ! sender one level up
                        data_corner = hvy_block( g+3:g+1+g+g:2, g+3:g+1+g+g:2, g+3:g+1+g+g:2, dF, my_block )

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                    ! write buffer
                    do i = 1, g
                        do j = 1, g
                            send_buff(buffer_i:buffer_i+g-1)    = data_corner( i, j, 1:g )
                            buffer_i                            = buffer_i + g
                        end do
                    end do

                end do

            ! '__1/123'
            case(27)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, Bs+g
                                send_buff(buffer_i:buffer_i+Bs+g-1)     = data_face_fine(1:Bs+g, j, Bs+i-1)
                                buffer_i                                = buffer_i + Bs+g
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, (Bs+1)/2
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+1:Bs+g:2, g+2*j-1, Bs+g-(2*i), dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '__1/134'
            case(28)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, Bs+g
                                send_buff(buffer_i:buffer_i+Bs+g-1)     = data_face_fine(1:Bs+g, j, Bs+i-1)
                                buffer_i                                = buffer_i + Bs+g
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, (Bs+1)/2
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+1:Bs+g:2, g+2*j-1, Bs+g-(2*i), dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '__1/145'
            case(29)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, Bs+g
                                send_buff(buffer_i:buffer_i+Bs+g-1)     = data_face_fine(1:Bs+g, j, Bs+i-1)
                                buffer_i                                = buffer_i + Bs+g
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, (Bs+1)/2
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+1:Bs+g:2, g+2*j-1, Bs+g-(2*i), dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '__1/152'
            case(30)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, Bs+g
                                send_buff(buffer_i:buffer_i+Bs+g-1)     = data_face_fine(1:Bs+g, j, Bs+i-1)
                                buffer_i                                = buffer_i + Bs+g
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, (Bs+1)/2
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+1:Bs+g:2, g+2*j-1, Bs+g-(2*i), dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '__2/123'
            case(31)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, Bs+g
                                send_buff(buffer_i:buffer_i+Bs+g-1)     = data_face_fine(1:Bs+g, 1+i, j)
                                buffer_i                                = buffer_i + Bs+g
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, (Bs+1)/2
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+1:Bs+g:2, g+(2*i)+1, g+2*j-1, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '__2/623'
            case(32)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, Bs+g
                                send_buff(buffer_i:buffer_i+Bs+g-1)     = data_face_fine(1:Bs+g, 1+i, j)
                                buffer_i                                = buffer_i + Bs+g
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, (Bs+1)/2
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+1:Bs+g:2, g+(2*i)+1, g+2*j-1, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '__2/152'
            case(33)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, Bs+g
                                send_buff(buffer_i:buffer_i+Bs+g-1)     = data_face_fine(1:Bs+g, 1+i, j)
                                buffer_i                                = buffer_i + Bs+g
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, (Bs+1)/2
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+1:Bs+g:2, g+(2*i)+1, g+2*j-1, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '__2/652'
            case(34)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, Bs+g
                                send_buff(buffer_i:buffer_i+Bs+g-1)     = data_face_fine(1:Bs+g, 1+i, j)
                                buffer_i                                = buffer_i + Bs+g
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, (Bs+1)/2
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+1:Bs+g:2, g+(2*i)+1, g+2*j-1, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '__3/123'
            case(35)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, Bs+g
                                send_buff(buffer_i:buffer_i+Bs+g-1)     = data_face_fine(Bs+i-1, 1:Bs+g, j)
                                buffer_i                                = buffer_i + Bs+g
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, (Bs+1)/2
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( Bs+g-(2*i), g+1:Bs+g:2, g+2*j-1, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '__3/623'
            case(36)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, Bs+g
                                send_buff(buffer_i:buffer_i+Bs+g-1)     = data_face_fine(Bs+i-1, 1:Bs+g, j)
                                buffer_i                                = buffer_i + Bs+g
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, (Bs+1)/2
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( Bs+g-(2*i), g+1:Bs+g:2, g+2*j-1, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '__3/134'
            case(37)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, Bs+g
                                send_buff(buffer_i:buffer_i+Bs+g-1)     = data_face_fine(Bs+i-1, 1:Bs+g, j)
                                buffer_i                                = buffer_i + Bs+g
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, (Bs+1)/2
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( Bs+g-(2*i), g+1:Bs+g:2, g+2*j-1, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '__3/634'
            case(38)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, Bs+g
                                send_buff(buffer_i:buffer_i+Bs+g-1)     = data_face_fine(Bs+i-1, 1:Bs+g, j)
                                buffer_i                                = buffer_i + Bs+g
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, (Bs+1)/2
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( Bs+g-(2*i), g+1:Bs+g:2, g+2*j-1, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '__4/134'
            case(39)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, Bs+g
                                send_buff(buffer_i:buffer_i+Bs+g-1)     = data_face_fine(1:Bs+g, Bs+i-1, j)
                                buffer_i                                = buffer_i + Bs+g
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, (Bs+1)/2
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+1:Bs+g:2, Bs+g-(2*i), g+2*j-1, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '__4/634'
            case(40)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, Bs+g
                                send_buff(buffer_i:buffer_i+Bs+g-1)     = data_face_fine(1:Bs+g, Bs+i-1, j)
                                buffer_i                                = buffer_i + Bs+g
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, (Bs+1)/2
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+1:Bs+g:2, Bs+g-(2*i), g+2*j-1, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '__4/145'
            case(41)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, Bs+g
                                send_buff(buffer_i:buffer_i+Bs+g-1)     = data_face_fine(1:Bs+g, Bs+i-1, j)
                                buffer_i                                = buffer_i + Bs+g
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, (Bs+1)/2
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+1:Bs+g:2, Bs+g-(2*i), g+2*j-1, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '__4/645'
            case(42)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, Bs+g
                                send_buff(buffer_i:buffer_i+Bs+g-1)     = data_face_fine(1:Bs+g, Bs+i-1, j)
                                buffer_i                                = buffer_i + Bs+g
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, (Bs+1)/2
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+1:Bs+g:2, Bs+g-(2*i), g+2*j-1, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '__5/145'
            case(43)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, Bs+g
                                send_buff(buffer_i:buffer_i+Bs+g-1)     = data_face_fine(1+i, 1:Bs+g, j)
                                buffer_i                                = buffer_i + Bs+g
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, (Bs+1)/2
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+(2*i)+1, g+1:Bs+g:2, g+2*j-1, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '__5/645'
            case(44)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, Bs+g
                                send_buff(buffer_i:buffer_i+Bs+g-1)     = data_face_fine(1+i, 1:Bs+g, j)
                                buffer_i                                = buffer_i + Bs+g
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, (Bs+1)/2
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+(2*i)+1, g+1:Bs+g:2, g+2*j-1, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '__5/152'
            case(45)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, Bs+g
                                send_buff(buffer_i:buffer_i+Bs+g-1)     = data_face_fine(1+i, 1:Bs+g, j)
                                buffer_i                                = buffer_i + Bs+g
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, (Bs+1)/2
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+(2*i)+1, g+1:Bs+g:2, g+2*j-1, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '__5/652'
            case(46)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, Bs+g
                                send_buff(buffer_i:buffer_i+Bs+g-1)     = data_face_fine(1+i, 1:Bs+g, j)
                                buffer_i                                = buffer_i + Bs+g
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, (Bs+1)/2
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+(2*i)+1, g+1:Bs+g:2, g+2*j-1, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '__6/623'
            case(47)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g , g+1:(Bs+1)/2+g/2+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, Bs+g
                                send_buff(buffer_i:buffer_i+Bs+g-1)     = data_face_fine(1:Bs+g, j, 1+i)
                                buffer_i                                = buffer_i + Bs+g
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, (Bs+1)/2
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+1:Bs+g:2, g+2*j-1, g+(2*i)+1, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '__6/634'
            case(48)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, Bs+g
                                send_buff(buffer_i:buffer_i+Bs+g-1)     = data_face_fine(1:Bs+g, j, 1+i)
                                buffer_i                                = buffer_i + Bs+g
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, (Bs+1)/2
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+1:Bs+g:2, g+2*j-1, g+(2*i)+1, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '__6/645'
            case(49)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, Bs+g
                                send_buff(buffer_i:buffer_i+Bs+g-1)     = data_face_fine(1:Bs+g, j, 1+i)
                                buffer_i                                = buffer_i + Bs+g
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, (Bs+1)/2
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+1:Bs+g:2, g+2*j-1, g+(2*i)+1, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '__6/652'
            case(50)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, Bs+g
                                send_buff(buffer_i:buffer_i+Bs+g-1)     = data_face_fine(1:Bs+g, j, 1+i)
                                buffer_i                                = buffer_i + Bs+g
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, (Bs+1)/2
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+1:Bs+g:2, g+2*j-1, g+(2*i)+1, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '_12/123'
            case(51)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+Bs-1)     = data_face_fine(g+1:Bs+g, 1+i, Bs+j-1)
                                buffer_i                                = buffer_i + Bs
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+1:Bs+g:2, g+(2*i)+1, Bs+g-(2*j), dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '_12/152'
            case(52)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+Bs-1)     = data_face_fine(g+1:Bs+g, 1+i, Bs+j-1)
                                buffer_i                                = buffer_i + Bs
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+1:Bs+g:2, g+(2*i)+1, Bs+g-(2*j), dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '_13/123'
            case(53)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+Bs-1)     = data_face_fine(Bs+i-1, g+1:Bs+g, Bs+j-1)
                                buffer_i                                = buffer_i + Bs
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( Bs+g-(2*i), g+1:Bs+g:2, Bs+g-(2*j), dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '_13/134'
            case(54)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+Bs-1)     = data_face_fine(Bs+i-1, g+1:Bs+g, Bs+j-1)
                                buffer_i                                = buffer_i + Bs
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( Bs+g-(2*i), g+1:Bs+g:2, Bs+g-(2*j), dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '_14/134'
            case(55)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+Bs-1)     = data_face_fine(g+1:Bs+g, Bs+i-1, Bs+j-1)
                                buffer_i                                = buffer_i + Bs
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+1:Bs+g:2, Bs+g-(2*i), Bs+g-(2*j), dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '_14/145'
            case(56)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+Bs-1)     = data_face_fine(g+1:Bs+g, Bs+i-1, Bs+j-1)
                                buffer_i                                = buffer_i + Bs
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+1:Bs+g:2, Bs+g-(2*i), Bs+g-(2*j), dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '_15/145'
            case(57)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+Bs-1)     = data_face_fine(1+i, g+1:Bs+g, Bs+j-1)
                                buffer_i                                = buffer_i + Bs
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+(2*i)+1, g+1:Bs+g:2, Bs+g-(2*j), dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '_15/152'
            case(58)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+Bs-1)     = data_face_fine(1+i, g+1:Bs+g, Bs+j-1)
                                buffer_i                                = buffer_i + Bs
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+(2*i)+1, g+1:Bs+g:2, Bs+g-(2*j), dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

             ! '_62/623'
            case(59)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+Bs-1)     = data_face_fine(g+1:Bs+g, 1+i, 1+j)
                                buffer_i                                = buffer_i + Bs
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+1:Bs+g:2, g+(2*i)+1, g+(2*j)+1, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '_62/652'
            case(60)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+Bs-1)     = data_face_fine(g+1:Bs+g, 1+i, 1+j)
                                buffer_i                                = buffer_i + Bs
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+1:Bs+g:2, g+(2*i)+1, g+(2*j)+1, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '_63/623'
            case(61)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+Bs-1)     = data_face_fine(Bs+i-1, g+1:Bs+g, 1+j)
                                buffer_i                                = buffer_i + Bs
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( Bs+g-(2*i), g+1:Bs+g:2, g+(2*j)+1, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '_63/634'
            case(62)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+Bs-1)     = data_face_fine(Bs+i-1, g+1:Bs+g, 1+j)
                                buffer_i                                = buffer_i + Bs
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( Bs+g-(2*i), g+1:Bs+g:2, g+(2*j)+1, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '_64/634'
            case(63)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+Bs-1)     = data_face_fine(g+1:Bs+g, Bs+i-1, 1+j)
                                buffer_i                                = buffer_i + Bs
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+1:Bs+g:2, Bs+g-(2*i), g+(2*j)+1, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '_64/645'
            case(64)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+Bs-1)     = data_face_fine(g+1:Bs+g, Bs+i-1, 1+j)
                                buffer_i                                = buffer_i + Bs
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+1:Bs+g:2, Bs+g-(2*i), g+(2*j)+1, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '_65/645'
            case(65)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+Bs-1)     = data_face_fine(1+i, g+1:Bs+g, 1+j)
                                buffer_i                                = buffer_i + Bs
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+(2*i)+1, g+1:Bs+g:2, g+(2*j)+1, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '_65/652'
            case(66)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+Bs-1)     = data_face_fine(1+i, g+1:Bs+g, 1+j)
                                buffer_i                                = buffer_i + Bs
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block( g+(2*i)+1, g+1:Bs+g:2, g+(2*j)+1, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '_23/123'
            case(67)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+Bs-1)     = data_face_fine(Bs+i-1, 1+j, g+1:Bs+g)
                                buffer_i                                = buffer_i + Bs
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block(Bs+g-(2*i), g+(2*j)+1, g+1:Bs+g:2, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '_23/623'
            case(68)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+Bs-1)     = data_face_fine(Bs+i-1, 1+j, g+1:Bs+g)
                                buffer_i                                = buffer_i + Bs
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block(Bs+g-(2*i), g+(2*j)+1, g+1:Bs+g:2, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '_25/152'
            case(69)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+Bs-1)     = data_face_fine(1+i, 1+j, g+1:Bs+g)
                                buffer_i                                = buffer_i + Bs
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block(g+(2*i)+1, g+(2*j)+1, g+1:Bs+g:2, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '_25/652'
            case(70)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+Bs-1)     = data_face_fine(1+i, 1+j, g+1:Bs+g)
                                buffer_i                                = buffer_i + Bs
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block(g+(2*i)+1, g+(2*j)+1, g+1:Bs+g:2, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '_43/134'
            case(71)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+Bs-1)     = data_face_fine(Bs+i-1, Bs+j-1, g+1:Bs+g)
                                buffer_i                                = buffer_i + Bs
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block(Bs+g-(2*i), Bs+g-(2*j), g+1:Bs+g:2, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '_43/634'
            case(72)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+Bs-1)     = data_face_fine(Bs+i-1, Bs+j-1, g+1:Bs+g)
                                buffer_i                                = buffer_i + Bs
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block(Bs+g-(2*i), Bs+g-(2*j), g+1:Bs+g:2, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '_45/145'
            case(73)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+Bs-1)     = data_face_fine(1+i, Bs+j-1, g+1:Bs+g)
                                buffer_i                                = buffer_i + Bs
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block(g+(2*i)+1, Bs+g-(2*j), g+1:Bs+g:2, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
                        end do

                    else
                        ! error case
                        write(*,'(80("_"))')
                        write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                        stop
                    end if

                end do

            ! '_45/645'
            case(74)
                do dF = 2, params%number_data_fields+1

                    if ( level_diff == -1 ) then
                        ! sender on lower level
                        ! data to interpolate
                        data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, my_block )
                        ! interpolate data
                        call prediction_3D( data_face , data_face_fine, params%order_predictor)
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+Bs-1)     = data_face_fine(1+i, Bs+j-1, g+1:Bs+g)
                                buffer_i                                = buffer_i + Bs
                            end do
                        end do

                    elseif ( level_diff == 1 ) then
                        ! sender on higher level
                        ! write data
                        do i = 1, g
                            do j = 1, g
                                send_buff(buffer_i:buffer_i+(Bs+1)/2-1)     = hvy_block(g+(2*i)+1, Bs+g-(2*j), g+1:Bs+g:2, dF, my_block )
                                buffer_i                                    = buffer_i + (Bs+1)/2
                            end do
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
    deallocate( data_face, stat=allocate_error )
    deallocate( data_face_fine, stat=allocate_error )

end subroutine create_send_buffer_3D
