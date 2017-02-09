! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: coarse_mesh_2D.f90
! version: 0.4
! author: msr
!
! coarse the mesh:
! every proc work on light data array,
!
! input:    - params, light and heavy data
! output:   - light and heavy data arrays
!
! = log ======================================================================================
!
! 08/11/16 - switch to v0.4, split old interpolate_mesh subroutine into two refine/coarsen
!            subroutines
! ********************************************************************************************

subroutine coarse_mesh_2D( params, lgt_block, hvy_block, lgt_active, lgt_n )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(in)      :: params
    ! light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)
    ! heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :)

    ! list of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_active(:)
    ! number of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_n

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
    real(kind=rk), allocatable          :: new_data(:,:,:), send_receive_data(:,:)
    ! new coordinates vectors
    real(kind=rk), allocatable          :: new_coord_x(:), new_coord_y(:), send_receive_coord(:)

    ! allocation error variable
    integer(kind=ik)                    :: allocate_error

    ! treecode varaible
    integer(kind=ik)                    :: me(params%max_treelevel), s1(params%max_treelevel), s2(params%max_treelevel), s3(params%max_treelevel)
    ! max treecode level
    integer(kind=ik)                    :: maxtL
    ! mesh level
    integer(kind=ik)                    :: level

    ! list of block ids, proc ranks
    integer(kind=ik)                    :: light_ids(4), proc_rank(4), heavy_ids(4)

    ! sister ids
    integer(kind=ik)                    :: id(3)

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
    allocate( new_data(Bs, Bs, params%number_data_fields), stat=allocate_error )
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if

    ! new coordinates vectors
    allocate( new_coord_x(Bs), stat=allocate_error )
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if

    allocate( new_coord_y(Bs), stat=allocate_error )
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if

    ! send/receive data and coordinates
    allocate( send_receive_data(Bs, Bs), stat=allocate_error )
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if

    allocate( send_receive_coord(Bs), stat=allocate_error )
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if

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
                do l = 1, 4
                    ! sister treecode differs only on last element
                    if ( lgt_block( lgt_active(k), level ) /= l-1) then

                        i         = i + 1
                        me(level) = l-1
                        ! find block id
                        call does_block_exist(me, lgt_block, maxtL, exists, lgt_id, lgt_active, lgt_n)
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

                ! proc with block k keep the data
                if ( data_rank == rank ) then
                    ! creating new block data
                    new_data           = 9.0e9_rk
                    ! coordinates for new block
                    new_coord_x        = 9.0e9_rk
                    new_coord_y        = 9.0e9_rk
                    ! send_receive buffer
                    send_receive_coord = 9.0e9_rk
                    send_receive_data  = 9.0e9_rk
                end if

                ! loop over proc rank list, proc with light id block data collect new data and coarse block
                ! then save new data on light id, current block_num
                do i = 1, 4

                    ! first: collect coordinates data
                    if ( data_rank == rank ) then
                        ! collect data
                        if ( proc_rank(i) == rank ) then
                            ! no communication needed
                            select case(i)
                                case(1)
                                    ! sister 0
                                    new_coord_x(1:(Bs-1)/2+1)              = hvy_block( 1, 1:Bs:2, 1, heavy_ids(1) )
                                    new_coord_y(1:(Bs-1)/2+1)              = hvy_block( 2, 1:Bs:2, 1, heavy_ids(1) )
                                case(2)
                                    ! sister 1
                                    new_coord_x((Bs-1)/2+1:Bs)             = hvy_block( 1, 1:Bs:2, 1, heavy_ids(2) )
                                case(3)
                                    ! sister 2
                                    new_coord_y((Bs-1)/2+1:Bs)             = hvy_block( 2, 1:Bs:2, 1, heavy_ids(3) )
                                case(4)
                                    ! sister 3
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
                                case(2)
                                    ! sister 1
                                    ! receive coords
                                    call MPI_Recv(send_receive_coord, Bs, MPI_REAL8, proc_rank(i), tag, MPI_COMM_WORLD, status, ierr)
                                    new_coord_x((Bs-1)/2+1:Bs)             = send_receive_coord(1:Bs:2)
                                case(3)
                                    ! sister 2
                                    ! receive coords
                                    call MPI_Recv(send_receive_coord, Bs, MPI_REAL8, proc_rank(i), tag, MPI_COMM_WORLD, status, ierr)
                                    new_coord_y((Bs-1)/2+1:Bs)              = send_receive_coord(1:Bs:2)
                                case(4)
                                    ! sister 3
                                    ! nothing to do
                            end select
                        end if

                    elseif ( proc_rank(i) == rank ) then
                        ! send data
                        select case(i)
                            case(1)
                                ! sister 0
                                ! send coord
                                send_receive_coord = hvy_block( 1, 1:Bs, 1, heavy_ids(1) )
                                call MPI_Send( send_receive_coord, Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                                send_receive_coord = hvy_block( 2, 1:Bs, 1, heavy_ids(1) )
                                call MPI_Send( send_receive_coord, Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                            case(2)
                                ! sister 1
                                ! send coord
                                send_receive_coord = hvy_block( 1, 1:Bs, 1, heavy_ids(2) )
                                call MPI_Send( send_receive_coord, Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                            case(3)
                                ! sister 2
                                ! send coord
                                send_receive_coord = hvy_block( 2, 1:Bs, 1, heavy_ids(3) )
                                call MPI_Send( send_receive_coord, Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                            case(4)
                                ! sister 3
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
                                        new_data(1:(Bs-1)/2+1, 1:(Bs-1)/2+1, dF-1)   = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, dF, heavy_ids(1) )
                                    case(2)
                                        ! sister 1
                                        new_data(1:(Bs-1)/2+1, (Bs-1)/2+1:Bs, dF-1)  = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, dF, heavy_ids(2) )
                                    case(3)
                                        ! sister 2
                                        new_data((Bs-1)/2+1:Bs, 1:(Bs-1)/2+1, dF-1)  = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, dF, heavy_ids(3) )
                                    case(4)
                                        ! sister 3
                                        new_data((Bs-1)/2+1:Bs, (Bs-1)/2+1:Bs, dF-1) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, dF, heavy_ids(4) )
                                end select
                            else
                                ! receive data from other proc, note: not all blocks have to send coord vectors
                                select case(i)
                                    case(1)
                                        ! sister 0
                                        ! receive data
                                        call MPI_Recv(send_receive_data, Bs*Bs, MPI_REAL8, proc_rank(i), tag, MPI_COMM_WORLD, status, ierr)
                                        new_data(1:(Bs-1)/2+1, 1:(Bs-1)/2+1, dF-1)   = send_receive_data(1:Bs:2, 1:Bs:2)
                                    case(2)
                                        ! sister 1
                                        ! receive data
                                        call MPI_Recv(send_receive_data, Bs*Bs, MPI_REAL8, proc_rank(i), tag, MPI_COMM_WORLD, status, ierr)
                                        new_data(1:(Bs-1)/2+1, (Bs-1)/2+1:Bs, dF-1)  = send_receive_data(1:Bs:2, 1:Bs:2)
                                    case(3)
                                        ! sister 2
                                        ! receive data
                                        call MPI_Recv(send_receive_data, Bs*Bs, MPI_REAL8, proc_rank(i), tag, MPI_COMM_WORLD, status, ierr)
                                        new_data((Bs-1)/2+1:Bs, 1:(Bs-1)/2+1, dF-1)  = send_receive_data(1:Bs:2, 1:Bs:2)
                                    case(4)
                                        ! sister 3
                                        ! receive data
                                        call MPI_Recv(send_receive_data, Bs*Bs, MPI_REAL8, proc_rank(i), tag, MPI_COMM_WORLD, status, ierr)
                                        new_data((Bs-1)/2+1:Bs, (Bs-1)/2+1:Bs, dF-1) = send_receive_data(1:Bs:2, 1:Bs:2)
                                end select
                            end if

                        elseif ( proc_rank(i) == rank ) then
                            ! send data
                            select case(i)
                                case(1)
                                    ! sister 0
                                    ! send data
                                    send_receive_data = hvy_block( g+1:Bs+g, g+1:Bs+g, dF, heavy_ids(1) )
                                    call MPI_Send( send_receive_data, Bs*Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                                case(2)
                                    ! sister 1
                                    ! send data
                                    send_receive_data = hvy_block( g+1:Bs+g, g+1:Bs+g, dF, heavy_ids(2) )
                                    call MPI_Send( send_receive_data, Bs*Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                                case(3)
                                    ! sister 2
                                    ! send data
                                    send_receive_data = hvy_block( g+1:Bs+g, g+1:Bs+g, dF, heavy_ids(3) )
                                    call MPI_Send( send_receive_data, Bs*Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                                case(4)
                                    ! sister 3
                                    ! send data
                                    send_receive_data = hvy_block( g+1:Bs+g, g+1:Bs+g, dF, heavy_ids(4) )
                                    call MPI_Send( send_receive_data, Bs*Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                            end select
                        else
                            ! nothing to do
                        end if
                    end do

                end do

                ! delete all heavy data
                if ( proc_rank(1) == rank ) hvy_block( :, :, :, heavy_ids(1) ) = 0.0_rk !9.0e9_rk
                if ( proc_rank(2) == rank ) hvy_block( :, :, :, heavy_ids(2) ) = 0.0_rk !9.0e9_rk
                if ( proc_rank(3) == rank ) hvy_block( :, :, :, heavy_ids(3) ) = 0.0_rk !9.0e9_rk
                if ( proc_rank(4) == rank ) hvy_block( :, :, :, heavy_ids(4) ) = 0.0_rk !9.0e9_rk

                ! delete light data
                lgt_block(light_ids(1), : ) = -1
                lgt_block(light_ids(2), : ) = -1
                lgt_block(light_ids(3), : ) = -1
                lgt_block(light_ids(4), : ) = -1

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
                    hvy_block( 1, 1:Bs, 1, hvy_id ) = new_coord_x
                    hvy_block( 2, 1:Bs, 1, hvy_id ) = new_coord_y

                    ! loop over all dataFields
                    do dF = 2, params%number_data_fields+1

                        ! write data
                        hvy_block( g+1:Bs+g, g+1:Bs+g, dF, hvy_id ) = new_data(:,:,dF-1)

                    end do

                end if

            end if

        end if

    end do

    ! clean up
    deallocate( new_data, stat=allocate_error )
    deallocate( new_coord_x, stat=allocate_error )
    deallocate( new_coord_y, stat=allocate_error )
    deallocate( send_receive_data, stat=allocate_error )
    deallocate( send_receive_coord, stat=allocate_error )

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

end subroutine coarse_mesh_2D
