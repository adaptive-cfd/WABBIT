! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: coarse_mesh.f90
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

subroutine coarse_mesh( params, block_list, block_data )

!---------------------------------------------------------------------------------------------
! modules

    use mpi
    ! global parameters
    use module_params
    ! interpolation routines
    use module_interpolation

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(in)      :: params
    ! light data array
    integer(kind=ik), intent(inout)     :: block_list(:, :)
    ! heavy data array - block data
    real(kind=rk), intent(inout)        :: block_data(:, :, :, :)

    ! loop variables
    integer(kind=ik)                    :: k, N, dF, i, l, j

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

    ! function to compare treecodes
    logical                             :: array_compare

    ! rank of proc to keep the coarsen data
    integer(kind=ik)                    :: data_rank

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    N = size(block_list, 1)

    maxtL = params%max_treelevel

    ! determinate process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    id = 0

    tag = 0

    ! grid parameter
    Bs = params%number_block_nodes
    g  = params%number_ghost_nodes

    ! allocate data array
    allocate( new_data(Bs, Bs, params%number_data_fields), stat=allocate_error )
    ! new coordinates vectors
    allocate( new_coord_x(Bs), stat=allocate_error )
    allocate( new_coord_y(Bs), stat=allocate_error )

    ! send/receive data and coordinates
    allocate( send_receive_data(Bs, Bs), stat=allocate_error )
    allocate( send_receive_coord(Bs), stat=allocate_error )

!---------------------------------------------------------------------------------------------
! main body

    ! loop over all light data
    do k = 1, N

        ! rank of proc to keep the data
        data_rank = (k-1) / params%number_blocks

        ! block is active and wants to coarsen
        if ( (block_list(k, 1) /= -1) .and. (block_list(k, params%max_treelevel+2) == -2) ) then

            ! get treecodes for current block and sister-blocks
            me                                              = block_list(k, 1:maxtL)
            light_ids( me( block_list(k, maxtL+1) ) + 1 )   = k
            heavy_ids( me( block_list(k, maxtL+1) ) + 1 )   = k - ((k-1) / params%number_blocks) * params%number_blocks
            proc_rank( me( block_list(k, maxtL+1) ) + 1 )   = (k-1) / params%number_blocks

            ! get sister ids
            level    = block_list( k, maxtL+1 )

            i = 0
            do l = 1, 4
                ! sister treecode differs only on last element
                if ( block_list( k, level ) /= l-1) then

                    i         = i + 1
                    me(level) = l-1
                    ! find block id
                    do j = 1, N
                        if ( block_list( j, 1 ) /= -1 ) then
                            if (array_compare( block_list( j, 1:maxtL ), me, maxtL)) then
                                id(i) = j
                            end if
                        end if
                    end do

                end if
            end do

            s1                                                  = block_list(id(1), 1:maxtL)
            light_ids( s1( block_list(id(1), maxtL+1) ) + 1 )   = id(1)
            heavy_ids( s1( block_list(id(1), maxtL+1) ) + 1 )   = id(1) - ((id(1)-1) / params%number_blocks) * params%number_blocks
            proc_rank( s1( block_list(id(1), maxtL+1) ) + 1 )   = (id(1)-1) / params%number_blocks

            s2                                                  = block_list(id(2), 1:maxtL)
            light_ids( s2( block_list(id(2), maxtL+1) ) + 1 )   = id(2)
            heavy_ids( s2( block_list(id(2), maxtL+1) ) + 1 )   = id(2) - ((id(2)-1) / params%number_blocks) * params%number_blocks
            proc_rank( s2( block_list(id(2), maxtL+1) ) + 1 )   = (id(2)-1) / params%number_blocks

            s3                                                  = block_list(id(3), 1:maxtL)
            light_ids( s3( block_list(id(3), maxtL+1) ) + 1 )   = id(3)
            heavy_ids( s3( block_list(id(3), maxtL+1) ) + 1 )   = id(3) - ((id(3)-1) / params%number_blocks) * params%number_blocks
            proc_rank( s3( block_list(id(3), maxtL+1) ) + 1 )   = (id(3)-1) / params%number_blocks

            ! proc with block k keep the data
            if ( data_rank == rank ) then
                ! creating new block data
                new_data     = 9.0e9_rk
                ! coordinates for new block
                new_coord_x  = 9.0e9_rk
                new_coord_y  = 9.0e9_rk
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
                                new_coord_x(1:(Bs-1)/2+1)              = block_data( 1, 1:Bs:2, 1, heavy_ids(1) )
                                new_coord_y(1:(Bs-1)/2+1)              = block_data( 2, 1:Bs:2, 1, heavy_ids(1) )
                            case(2)
                                ! sister 1
                                new_coord_x((Bs-1)/2+1:Bs)             = block_data( 1, 1:Bs:2, 1, heavy_ids(2) )
                            case(3)
                                ! sister 2
                                new_coord_y((Bs-1)/2+1:Bs)             = block_data( 2, 1:Bs:2, 1, heavy_ids(3) )
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
                            send_receive_coord = block_data( 1, 1:Bs, 1, heavy_ids(1) )
                            call MPI_Send( send_receive_coord, Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                            send_receive_coord = block_data( 2, 1:Bs, 1, heavy_ids(1) )
                            call MPI_Send( send_receive_coord, Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                        case(2)
                            ! sister 1
                            ! send coord
                            send_receive_coord = block_data( 1, 1:Bs, 1, heavy_ids(2) )
                            call MPI_Send( send_receive_coord, Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                        case(3)
                            ! sister 2
                            ! send coord
                            send_receive_coord = block_data( 2, 1:Bs, 1, heavy_ids(3) )
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
                                    new_data(1:(Bs-1)/2+1, 1:(Bs-1)/2+1, dF-1)   = block_data( g+1:Bs+g:2, g+1:Bs+g:2, dF, heavy_ids(1) )
                                case(2)
                                    ! sister 1
                                    new_data(1:(Bs-1)/2+1, (Bs-1)/2+1:Bs, dF-1)  = block_data( g+1:Bs+g:2, g+1:Bs+g:2, dF, heavy_ids(2) )
                                case(3)
                                    ! sister 2
                                    new_data((Bs-1)/2+1:Bs, 1:(Bs-1)/2+1, dF-1)  = block_data( g+1:Bs+g:2, g+1:Bs+g:2, dF, heavy_ids(3) )
                                case(4)
                                    ! sister 3
                                    new_data((Bs-1)/2+1:Bs, (Bs-1)/2+1:Bs, dF-1) = block_data( g+1:Bs+g:2, g+1:Bs+g:2, dF, heavy_ids(4) )
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
                                send_receive_data = block_data( g+1:Bs+g, g+1:Bs+g, dF, heavy_ids(1) )
                                call MPI_Send( send_receive_data, Bs*Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                            case(2)
                                ! sister 1
                                ! send data
                                send_receive_data = block_data( g+1:Bs+g, g+1:Bs+g, dF, heavy_ids(2) )
                                call MPI_Send( send_receive_data, Bs*Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                            case(3)
                                ! sister 2
                                ! send data
                                send_receive_data = block_data( g+1:Bs+g, g+1:Bs+g, dF, heavy_ids(3) )
                                call MPI_Send( send_receive_data, Bs*Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                            case(4)
                                ! sister 3
                                ! send data
                                send_receive_data = block_data( g+1:Bs+g, g+1:Bs+g, dF, heavy_ids(4) )
                                call MPI_Send( send_receive_data, Bs*Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                        end select
                    else
                        ! nothing to do
                    end if
                end do

            end do

            ! delete all heavy data
            if ( proc_rank(1) == rank ) block_data( :, :, :, heavy_ids(1) ) = 0.0_rk
            if ( proc_rank(2) == rank ) block_data( :, :, :, heavy_ids(2) ) = 0.0_rk
            if ( proc_rank(3) == rank ) block_data( :, :, :, heavy_ids(3) ) = 0.0_rk
            if ( proc_rank(4) == rank ) block_data( :, :, :, heavy_ids(4) ) = 0.0_rk

            ! delete light data
            block_list(light_ids(1), : ) = -1
            block_list(light_ids(2), : ) = -1
            block_list(light_ids(3), : ) = -1
            block_list(light_ids(4), : ) = -1

            ! new treecode, one level down (coarsening)
            me(level) = -1

            ! write new block on current block
            ! write light data
            block_list(k, 1:maxtL) = me
            block_list(k, maxtL+1) = level-1
            block_list(k, maxtL+2) = 0

            ! final steps
            if ( data_rank == rank ) then

                ! write coordinates
                block_data( 1, 1:Bs, 1, k - ((k-1) / params%number_blocks) * params%number_blocks ) = new_coord_x
                block_data( 2, 1:Bs, 1, k - ((k-1) / params%number_blocks) * params%number_blocks ) = new_coord_y

                ! loop over all dataFields
                do dF = 2, params%number_data_fields+1

                    ! write data
                    block_data( g+1:Bs+g, g+1:Bs+g, dF, k - ((k-1) / params%number_blocks) * params%number_blocks ) = new_data(:,:,dF-1)

                end do

            end if

        end if

    end do

    ! clean up
    deallocate( new_data, stat=allocate_error )
    deallocate( new_coord_x, stat=allocate_error )
    deallocate( new_coord_y, stat=allocate_error )
    deallocate( send_receive_data, stat=allocate_error )
    deallocate( send_receive_coord, stat=allocate_error )

end subroutine coarse_mesh
