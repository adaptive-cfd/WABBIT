!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name coarse_mesh_3D.f90
!> \version 0.5
!> \author msr
!
!> \brief coarse the mesh: \n
!! every proc work on light data array
!
!> \details
!! input:    - params, light and heavy data \n
!! output:   - light and heavy data arrays
!! \n
!> = log ======================================================================================
!! \n
!! 03/02/17 - create
!
! ********************************************************************************************

subroutine coarse_mesh_3D( params, lgt_block, hvy_block, lgt_active, lgt_n, lgt_sortednumlist )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> list of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_n
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), intent(inout)  :: lgt_sortednumlist(:,:)

    ! loop variables
    integer(kind=ik)                    :: k, dF, i, l, N, lgt_id, hvy_id, sister_rank

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank
    ! MPI message tag
    integer(kind=ik)                    :: tag
    ! MPI status
    integer                             :: status(MPI_status_size)

    ! grid parameter
    integer(kind=ik)                    :: Bs, g
    ! data fields for interpolation
    real(kind=rk), allocatable          :: new_data(:,:,:,:), send_receive_data(:,:,:)
    ! new coordinates vectors
    real(kind=rk), allocatable          :: new_coord_x(:), new_coord_y(:), new_coord_z(:), send_receive_coord(:)
    ! treecode varaible
    integer(kind=ik)                    :: me(params%max_treelevel), s1(params%max_treelevel), s2(params%max_treelevel), s3(params%max_treelevel), &
                                        s4(params%max_treelevel), s5(params%max_treelevel), s6(params%max_treelevel), s7(params%max_treelevel)
    ! max treecode level
    integer(kind=ik)                    :: maxtL
    ! mesh level
    integer(kind=ik)                    :: level

    ! list of block ids, proc ranks
    integer(kind=ik)                    :: light_ids(8), proc_rank(8), heavy_ids(8)

    ! sister ids
    integer(kind=ik)                    :: id(7)

    ! variable for block existence
    logical                             :: exists

    ! rank of proc to keep the coarsen data
    integer(kind=ik)                    :: data_rank

    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0, sub_t1

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! set MPI parameter
    rank         = params%rank

    ! start time
    sub_t0 = MPI_Wtime()

    N = params%number_blocks

    maxtL = params%max_treelevel

    id = -1

    tag = 0

    ! grid parameter
    Bs = params%number_block_nodes
    g  = params%number_ghost_nodes

    ! allocate data array
    allocate( new_data(Bs, Bs, Bs, params%number_data_fields) )

    ! new coordinates vectors
    allocate( new_coord_x(Bs) )
    allocate( new_coord_y(Bs) )
    allocate( new_coord_z(Bs) )

    ! send/receive data and coordinates
    allocate( send_receive_data(Bs, Bs, Bs) )
    allocate( send_receive_coord(Bs) )

!---------------------------------------------------------------------------------------------
! main body

    ! loop over all active light data
    do k = 1, lgt_n

        ! only work on light data, if block is active
        ! note: due to previous loops, some light data id are allready
        ! coarsen, but are still in active block list
        ! so: check additional active criteria
        if ( lgt_block( lgt_active(k), 1 ) /= -1 ) then

            ! rank of proc to keep the data
            call lgt_id_to_proc_rank( data_rank, lgt_active(k), N )

            ! block wants to coarsen
            if ( lgt_block( lgt_active(k), params%max_treelevel+2) == -2 ) then

                ! get treecodes for current block and sister-blocks
                me                                                          = lgt_block( lgt_active(k), 1:maxtL)
                light_ids( me( lgt_block( lgt_active(k), maxtL+1) ) + 1 )   = lgt_active(k)

                ! heavy id
                call lgt_id_to_hvy_id( hvy_id, lgt_active(k), data_rank, N )

                heavy_ids( me( lgt_block( lgt_active(k), maxtL+1) ) + 1 )   = hvy_id
                proc_rank( me( lgt_block( lgt_active(k), maxtL+1) ) + 1 )   = data_rank

                ! current block level
                level    = lgt_block( lgt_active(k), maxtL+1 )

                ! reset id array
                id = -1

                ! get sister ids
                i = 0
                do l = 1, 8
                    ! sister treecode differs only on last element
                    if ( lgt_block( lgt_active(k), level ) /= l-1) then

                        i         = i + 1
                        me(level) = l-1
                        ! find block id
                        call does_block_exist(me, exists, lgt_id, lgt_sortednumlist, lgt_n)
                        ! block exists
                        if (exists) then
                            id(i) = lgt_id
                        else
                            ! error case
                            print*, "ERROR: can not find sister id"
                            stop
                        end if

                    end if
                end do

                ! -------------------------------------------------------------------------------------
                ! id1 data
                s1                                                 = lgt_block(id(1), 1:maxtL)
                light_ids( s1( lgt_block(id(1), maxtL+1) ) + 1 )   = id(1)

                ! rank of id1 proc
                call lgt_id_to_proc_rank( sister_rank, id(1), N )
                ! heavy id
                call lgt_id_to_hvy_id( hvy_id, id(1), sister_rank, N )

                heavy_ids( s1( lgt_block(id(1), maxtL+1) ) + 1 )   = hvy_id
                proc_rank( s1( lgt_block(id(1), maxtL+1) ) + 1 )   = sister_rank

                ! -------------------------------------------------------------------------------------
                ! id2 data
                s2                                                 = lgt_block(id(2), 1:maxtL)
                light_ids( s2( lgt_block(id(2), maxtL+1) ) + 1 )   = id(2)

                ! rank of id2 proc
                call lgt_id_to_proc_rank( sister_rank, id(2), N )
                ! heavy id
                call lgt_id_to_hvy_id( hvy_id, id(2), sister_rank, N )

                heavy_ids( s2( lgt_block(id(2), maxtL+1) ) + 1 )   = hvy_id
                proc_rank( s2( lgt_block(id(2), maxtL+1) ) + 1 )   = sister_rank

                ! -------------------------------------------------------------------------------------
                ! id3 data
                s3                                                 = lgt_block(id(3), 1:maxtL)
                light_ids( s3( lgt_block(id(3), maxtL+1) ) + 1 )   = id(3)

                ! rank of id3 proc
                call lgt_id_to_proc_rank( sister_rank, id(3), N )
                ! heavy id
                call lgt_id_to_hvy_id( hvy_id, id(3), sister_rank, N )

                heavy_ids( s3( lgt_block(id(3), maxtL+1) ) + 1 )   = hvy_id
                proc_rank( s3( lgt_block(id(3), maxtL+1) ) + 1 )   = sister_rank

                ! -------------------------------------------------------------------------------------
                ! id4 data
                s4                                                 = lgt_block(id(4), 1:maxtL)
                light_ids( s4( lgt_block(id(4), maxtL+1) ) + 1 )   = id(4)

                ! rank of id4 proc
                call lgt_id_to_proc_rank( sister_rank, id(4), N )
                ! heavy id
                call lgt_id_to_hvy_id( hvy_id, id(4), sister_rank, N )

                heavy_ids( s4( lgt_block(id(4), maxtL+1) ) + 1 )   = hvy_id
                proc_rank( s4( lgt_block(id(4), maxtL+1) ) + 1 )   = sister_rank

                ! -------------------------------------------------------------------------------------
                ! id5 data
                s5                                                 = lgt_block(id(5), 1:maxtL)
                light_ids( s5( lgt_block(id(5), maxtL+1) ) + 1 )   = id(5)

                ! rank of id5 proc
                call lgt_id_to_proc_rank( sister_rank, id(5), N )
                ! heavy id
                call lgt_id_to_hvy_id( hvy_id, id(5), sister_rank, N )

                heavy_ids( s5( lgt_block(id(5), maxtL+1) ) + 1 )   = hvy_id
                proc_rank( s5( lgt_block(id(5), maxtL+1) ) + 1 )   = sister_rank

                ! -------------------------------------------------------------------------------------
                ! id6 data
                s6                                                 = lgt_block(id(6), 1:maxtL)
                light_ids( s6( lgt_block(id(6), maxtL+1) ) + 1 )   = id(6)

                ! rank of id6 proc
                call lgt_id_to_proc_rank( sister_rank, id(6), N )
                ! heavy id
                call lgt_id_to_hvy_id( hvy_id, id(6), sister_rank, N )

                heavy_ids( s6( lgt_block(id(6), maxtL+1) ) + 1 )   = hvy_id
                proc_rank( s6( lgt_block(id(6), maxtL+1) ) + 1 )   = sister_rank

                ! -------------------------------------------------------------------------------------
                ! id7 data
                s7                                                 = lgt_block(id(7), 1:maxtL)
                light_ids( s7( lgt_block(id(7), maxtL+1) ) + 1 )   = id(7)

                ! rank of id7 proc
                call lgt_id_to_proc_rank( sister_rank, id(7), N )
                ! heavy id
                call lgt_id_to_hvy_id( hvy_id, id(7), sister_rank, N )

                heavy_ids( s7( lgt_block(id(7), maxtL+1) ) + 1 )   = hvy_id
                proc_rank( s7( lgt_block(id(7), maxtL+1) ) + 1 )   = sister_rank

                ! -------------------------------------------------------------------------------------
                ! proc with block k keep the data
!                if ( data_rank == rank ) then
!                    ! creating new block data
!                    new_data           = 9.0e9_rk
!                    ! coordinates for new block
!                    new_coord_x        = 9.0e9_rk
!                    new_coord_y        = 9.0e9_rk
!                    new_coord_z        = 9.0e9_rk
!                    ! send_receive buffer
!                    send_receive_coord = 7.0e9_rk
!                    send_receive_data  = 7.0e9_rk
!                end if

                ! loop over proc rank list, proc with light id block data collect new data and coarse block
                ! then save new data on light id, current block_num
                do i = 1, 8

                    ! first: collect coordinates data
                    if ( data_rank == rank ) then
                        ! collect data
                        if ( proc_rank(i) == rank ) then
                            ! no communication needed
                            select case(i)
                                case(1)
                                    ! sister 0
                                    new_coord_x(1:(Bs-1)/2+1)              = hvy_block( 1, 1:Bs:2, 1, 1, heavy_ids(1) )
                                    new_coord_y(1:(Bs-1)/2+1)              = hvy_block( 2, 1:Bs:2, 1, 1, heavy_ids(1) )
                                    new_coord_z(1:(Bs-1)/2+1)              = hvy_block( 3, 1:Bs:2, 1, 1, heavy_ids(1) )
                                case(2)
                                    ! sister 1
                                    new_coord_y((Bs-1)/2+1:Bs)             = hvy_block( 2, 1:Bs:2, 1, 1, heavy_ids(2) )
                                case(3)
                                    ! sister 2
                                    new_coord_x((Bs-1)/2+1:Bs)             = hvy_block( 1, 1:Bs:2, 1, 1, heavy_ids(3) )
                                case(4)
                                    ! sister 3
                                    ! nothing to do
                                case(5)
                                    ! sister 4
                                    new_coord_z((Bs-1)/2+1:Bs)             = hvy_block( 3, 1:Bs:2, 1, 1, heavy_ids(5) )
                                case(6)
                                    ! sister 5
                                    ! nothing to do
                                case(7)
                                    ! sister 6
                                    ! nothing to do
                                case(8)
                                    ! sister 7
                                    ! nothing to do
                            end select
                        else
                            ! receive data from other proc, note: not all blocks have to send coord vectors
                            select case(i)
                                case(1)
                                    ! sister 0
                                    ! receive coords
                                    call MPI_Recv(send_receive_coord, Bs, MPI_REAL8, proc_rank(i), tag, MPI_COMM_WORLD, status, ierr)
                                    new_coord_x(1:(Bs-1)/2+1)              = send_receive_coord(1:Bs:2)
                                    call MPI_Recv(send_receive_coord, Bs, MPI_REAL8, proc_rank(i), tag, MPI_COMM_WORLD, status, ierr)
                                    new_coord_y(1:(Bs-1)/2+1)              = send_receive_coord(1:Bs:2)
                                    call MPI_Recv(send_receive_coord, Bs, MPI_REAL8, proc_rank(i), tag, MPI_COMM_WORLD, status, ierr)
                                    new_coord_z(1:(Bs-1)/2+1)              = send_receive_coord(1:Bs:2)
                                case(2)
                                    ! sister 1
                                    ! receive coords
                                    call MPI_Recv(send_receive_coord, Bs, MPI_REAL8, proc_rank(i), tag, MPI_COMM_WORLD, status, ierr)
                                    new_coord_y((Bs-1)/2+1:Bs)             = send_receive_coord(1:Bs:2)
                                case(3)
                                    ! sister 2
                                    ! receive coords
                                    call MPI_Recv(send_receive_coord, Bs, MPI_REAL8, proc_rank(i), tag, MPI_COMM_WORLD, status, ierr)
                                    new_coord_x((Bs-1)/2+1:Bs)              = send_receive_coord(1:Bs:2)
                                case(4)
                                    ! sister 3
                                    ! nothing to do
                                case(5)
                                    ! sister 4
                                    ! receive coords
                                    call MPI_Recv(send_receive_coord, Bs, MPI_REAL8, proc_rank(i), tag, MPI_COMM_WORLD, status, ierr)
                                    new_coord_z((Bs-1)/2+1:Bs)              = send_receive_coord(1:Bs:2)
                                case(6)
                                    ! sister 5
                                    ! nothing to do
                                case(7)
                                    ! sister 6
                                    ! nothing to do
                                case(8)
                                    ! sister 7
                                    ! nothing to do
                            end select
                        end if

                    elseif ( proc_rank(i) == rank ) then
                        ! send data
                        select case(i)
                            case(1)
                                ! sister 0
                                ! send coord
                                send_receive_coord = hvy_block( 1, 1:Bs, 1, 1, heavy_ids(1) )
                                call MPI_Send( send_receive_coord, Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                                send_receive_coord = hvy_block( 2, 1:Bs, 1, 1, heavy_ids(1) )
                                call MPI_Send( send_receive_coord, Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                                send_receive_coord = hvy_block( 3, 1:Bs, 1, 1, heavy_ids(1) )
                                call MPI_Send( send_receive_coord, Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                            case(2)
                                ! sister 1
                                ! send coord
                                send_receive_coord = hvy_block( 2, 1:Bs, 1, 1, heavy_ids(2) )
                                call MPI_Send( send_receive_coord, Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                            case(3)
                                ! sister 2
                                ! send coord
                                send_receive_coord = hvy_block( 1, 1:Bs, 1, 1, heavy_ids(3) )
                                call MPI_Send( send_receive_coord, Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                            case(4)
                                ! sister 3
                                ! nothing to do
                            case(5)
                                ! sister 4
                                ! send coord
                                send_receive_coord = hvy_block( 3, 1:Bs, 1, 1, heavy_ids(5) )
                                call MPI_Send( send_receive_coord, Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                            case(6)
                                ! sister 5
                                ! nothing to do
                            case(7)
                                ! sister 6
                                ! nothing to do
                            case(8)
                                ! sister 7
                                ! nothing to do
                        end select
                    else
                        ! nothing to do
                    end if

                    ! second: send data, loop over all datafields
                    do dF = 2, params%number_data_fields+1
                        if ( data_rank == rank ) then
                            ! collect data
                            if ( proc_rank(i) == rank ) then
                                ! no communication needed
                                select case(i)

                                    case(1)
                                        ! sister 0
                                        new_data(1:(Bs-1)/2+1, 1:(Bs-1)/2+1, 1:(Bs-1)/2+1, dF-1)   = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, heavy_ids(1) )
                                    case(2)
                                        ! sister 1
                                        new_data(1:(Bs-1)/2+1, (Bs-1)/2+1:Bs, 1:(Bs-1)/2+1, dF-1)  = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, heavy_ids(2) )
                                    case(3)
                                        ! sister 2
                                        new_data((Bs-1)/2+1:Bs, 1:(Bs-1)/2+1, 1:(Bs-1)/2+1, dF-1)  = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, heavy_ids(3) )
                                    case(4)
                                        ! sister 3
                                        new_data((Bs-1)/2+1:Bs, (Bs-1)/2+1:Bs, 1:(Bs-1)/2+1, dF-1) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, heavy_ids(4) )
                                    case(5)
                                        ! sister 4
                                        new_data(1:(Bs-1)/2+1, 1:(Bs-1)/2+1, (Bs-1)/2+1:Bs, dF-1)   = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, heavy_ids(5) )
                                    case(6)
                                        ! sister 5
                                        new_data(1:(Bs-1)/2+1, (Bs-1)/2+1:Bs, (Bs-1)/2+1:Bs, dF-1)  = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, heavy_ids(6) )
                                    case(7)
                                        ! sister 6
                                        new_data((Bs-1)/2+1:Bs, 1:(Bs-1)/2+1, (Bs-1)/2+1:Bs, dF-1)  = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, heavy_ids(7) )
                                    case(8)
                                        ! sister 7
                                        new_data((Bs-1)/2+1:Bs, (Bs-1)/2+1:Bs, (Bs-1)/2+1:Bs, dF-1) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, heavy_ids(8) )

                                end select
                            else
                                ! receive data from other proc, note: not all blocks have to send coord vectors
                                select case(i)
                                    case(1)
                                        ! sister 0
                                        ! receive data
                                        call MPI_Recv(send_receive_data, Bs*Bs*Bs, MPI_REAL8, proc_rank(i), tag, MPI_COMM_WORLD, status, ierr)
                                        new_data(1:(Bs-1)/2+1, 1:(Bs-1)/2+1, 1:(Bs-1)/2+1, dF-1)   = send_receive_data(1:Bs:2, 1:Bs:2, 1:Bs:2)
                                    case(2)
                                        ! sister 1
                                        ! receive data
                                        call MPI_Recv(send_receive_data, Bs*Bs*Bs, MPI_REAL8, proc_rank(i), tag, MPI_COMM_WORLD, status, ierr)
                                        new_data(1:(Bs-1)/2+1, (Bs-1)/2+1:Bs, 1:(Bs-1)/2+1, dF-1)  = send_receive_data(1:Bs:2, 1:Bs:2, 1:Bs:2)
                                    case(3)
                                        ! sister 2
                                        ! receive data
                                        call MPI_Recv(send_receive_data, Bs*Bs*Bs, MPI_REAL8, proc_rank(i), tag, MPI_COMM_WORLD, status, ierr)
                                        new_data((Bs-1)/2+1:Bs, 1:(Bs-1)/2+1, 1:(Bs-1)/2+1, dF-1)  = send_receive_data(1:Bs:2, 1:Bs:2, 1:Bs:2)
                                    case(4)
                                        ! sister 3
                                        ! receive data
                                        call MPI_Recv(send_receive_data, Bs*Bs*Bs, MPI_REAL8, proc_rank(i), tag, MPI_COMM_WORLD, status, ierr)
                                        new_data((Bs-1)/2+1:Bs, (Bs-1)/2+1:Bs, 1:(Bs-1)/2+1, dF-1) = send_receive_data(1:Bs:2, 1:Bs:2, 1:Bs:2)
                                    case(5)
                                        ! sister 4
                                        ! receive data
                                        call MPI_Recv(send_receive_data, Bs*Bs*Bs, MPI_REAL8, proc_rank(i), tag, MPI_COMM_WORLD, status, ierr)
                                        new_data(1:(Bs-1)/2+1, 1:(Bs-1)/2+1, (Bs-1)/2+1:Bs, dF-1)   = send_receive_data(1:Bs:2, 1:Bs:2, 1:Bs:2)
                                    case(6)
                                        ! sister 5
                                        ! receive data
                                        call MPI_Recv(send_receive_data, Bs*Bs*Bs, MPI_REAL8, proc_rank(i), tag, MPI_COMM_WORLD, status, ierr)
                                        new_data(1:(Bs-1)/2+1, (Bs-1)/2+1:Bs, (Bs-1)/2+1:Bs, dF-1)  = send_receive_data(1:Bs:2, 1:Bs:2, 1:Bs:2)
                                    case(7)
                                        ! sister 6
                                        ! receive data
                                        call MPI_Recv(send_receive_data, Bs*Bs*Bs, MPI_REAL8, proc_rank(i), tag, MPI_COMM_WORLD, status, ierr)
                                        new_data((Bs-1)/2+1:Bs, 1:(Bs-1)/2+1, (Bs-1)/2+1:Bs, dF-1)  = send_receive_data(1:Bs:2, 1:Bs:2, 1:Bs:2)
                                    case(8)
                                        ! sister 7
                                        ! receive data
                                        call MPI_Recv(send_receive_data, Bs*Bs*Bs, MPI_REAL8, proc_rank(i), tag, MPI_COMM_WORLD, status, ierr)
                                        new_data((Bs-1)/2+1:Bs, (Bs-1)/2+1:Bs, (Bs-1)/2+1:Bs, dF-1) = send_receive_data(1:Bs:2, 1:Bs:2, 1:Bs:2)
                                end select
                            end if

                        elseif ( proc_rank(i) == rank ) then
                            ! send data
                            select case(i)
                                case(1)
                                    ! sister 0
                                    ! send data
                                    send_receive_data = hvy_block( g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, dF, heavy_ids(1) )
                                    call MPI_Send( send_receive_data, Bs*Bs*Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                                case(2)
                                    ! sister 1
                                    ! send data
                                    send_receive_data = hvy_block( g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, dF, heavy_ids(2) )
                                    call MPI_Send( send_receive_data, Bs*Bs*Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                                case(3)
                                    ! sister 2
                                    ! send data
                                    send_receive_data = hvy_block( g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, dF, heavy_ids(3) )
                                    call MPI_Send( send_receive_data, Bs*Bs*Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                                case(4)
                                    ! sister 3
                                    ! send data
                                    send_receive_data = hvy_block( g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, dF, heavy_ids(4) )
                                    call MPI_Send( send_receive_data, Bs*Bs*Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                                case(5)
                                    ! sister 4
                                    ! send data
                                    send_receive_data = hvy_block( g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, dF, heavy_ids(5) )
                                    call MPI_Send( send_receive_data, Bs*Bs*Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                                case(6)
                                    ! sister 5
                                    ! send data
                                    send_receive_data = hvy_block( g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, dF, heavy_ids(6) )
                                    call MPI_Send( send_receive_data, Bs*Bs*Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                                case(7)
                                    ! sister 6
                                    ! send data
                                    send_receive_data = hvy_block( g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, dF, heavy_ids(7) )
                                    call MPI_Send( send_receive_data, Bs*Bs*Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                                case(8)
                                    ! sister 7
                                    ! send data
                                    send_receive_data = hvy_block( g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, dF, heavy_ids(8) )
                                    call MPI_Send( send_receive_data, Bs*Bs*Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                            end select
                        else
                            ! nothing to do
                        end if
                    end do

                end do

                ! delete all heavy data
                if ( proc_rank(1) == rank ) hvy_block( :, :, :, :, heavy_ids(1) ) = 0.0_rk !9.0e9_rk
                if ( proc_rank(2) == rank ) hvy_block( :, :, :, :, heavy_ids(2) ) = 0.0_rk !9.0e9_rk
                if ( proc_rank(3) == rank ) hvy_block( :, :, :, :, heavy_ids(3) ) = 0.0_rk !9.0e9_rk
                if ( proc_rank(4) == rank ) hvy_block( :, :, :, :, heavy_ids(4) ) = 0.0_rk !9.0e9_rk
                if ( proc_rank(5) == rank ) hvy_block( :, :, :, :, heavy_ids(5) ) = 0.0_rk !9.0e9_rk
                if ( proc_rank(6) == rank ) hvy_block( :, :, :, :, heavy_ids(6) ) = 0.0_rk !9.0e9_rk
                if ( proc_rank(7) == rank ) hvy_block( :, :, :, :, heavy_ids(7) ) = 0.0_rk !9.0e9_rk
                if ( proc_rank(8) == rank ) hvy_block( :, :, :, :, heavy_ids(8) ) = 0.0_rk !9.0e9_rk

                ! delete light data
                lgt_block(light_ids(1), : ) = -1
                lgt_block(light_ids(2), : ) = -1
                lgt_block(light_ids(3), : ) = -1
                lgt_block(light_ids(4), : ) = -1
                lgt_block(light_ids(5), : ) = -1
                lgt_block(light_ids(6), : ) = -1
                lgt_block(light_ids(7), : ) = -1
                lgt_block(light_ids(8), : ) = -1

                ! new treecode, one level down (coarsening)
                me(level) = -1

                ! write new block on current block
                ! write light data
                lgt_block( lgt_active(k), 1:maxtL) = me
                lgt_block( lgt_active(k), maxtL+1) = level-1
                lgt_block( lgt_active(k), maxtL+2) = 0

                ! final steps
                if ( data_rank == rank ) then

                    ! heavy id
                    call lgt_id_to_hvy_id( hvy_id, lgt_active(k), rank, N )

                    ! write coordinates
                    hvy_block( 1, 1:Bs, 1, 1, hvy_id ) = new_coord_x
                    hvy_block( 2, 1:Bs, 1, 1, hvy_id ) = new_coord_y
                    hvy_block( 3, 1:Bs, 1, 1, hvy_id ) = new_coord_z

                    ! loop over all dataFields
                    do dF = 2, params%number_data_fields+1

                        ! write data
                        hvy_block( g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, dF, hvy_id ) = new_data(:,:,:,dF-1)

                    end do

                end if

            end if

        end if

    end do

    ! clean up
    deallocate( new_data )
    deallocate( new_coord_x )
    deallocate( new_coord_y )
    deallocate( send_receive_data )
    deallocate( send_receive_coord )

    ! end time
    sub_t1 = MPI_Wtime()
    ! write time
    if ( params%debug ) then
        ! find free or corresponding line
        k = 1
        do while ( debug%name_comp_time(k) /= "---" )
            ! entry for current subroutine exists
            if ( debug%name_comp_time(k) == "coarse_mesh" ) exit
            k = k + 1
        end do
        ! write time
        debug%name_comp_time(k) = "coarse_mesh"
        debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1
        debug%comp_time(k, 2)   = debug%comp_time(k, 2) + sub_t1 - sub_t0
    end if

end subroutine coarse_mesh_3D
