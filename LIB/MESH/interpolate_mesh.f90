! ********************************
! WABBIT
! --------------------------------
!
! adapt the mesh with interpolation
!
! name: interpolate_mesh.f90
! date: 27.10.2016
! author: msr
! version: 0.3
!
! ********************************

subroutine interpolate_mesh()

    use mpi
    use module_params
    use module_blocks
    use module_interpolation

    implicit none

    integer(kind=ik)                                :: dF, k, N, block_num, id1, id2, id3, allocate_error, Bs, g, i, light_free_id, light_id, heavy_free_id, data_rank
    integer(kind=ik)                                :: rank, ierr, tag
    integer(kind=ik), dimension(10)                 :: me, s1, s2, s3, new_treecode
    integer(kind=ik), dimension(4)                  :: block_ids, proc_rank
    integer                                         :: status(MPI_status_size)

    real(kind=rk), dimension(:,:), allocatable      :: new_data, new_data_w_ghost, data_predict_fine, data_predict_coarse, send_receive_data
    real(kind=rk), dimension(:), allocatable        :: new_coord_x, new_coord_y, send_receive_coord

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    tag = 0

    ! loop over all light data
    N           = blocks_params%number_max_blocks

    Bs = blocks_params%size_block
    g  = blocks_params%number_ghost_nodes

    allocate( new_data(Bs, Bs), stat=allocate_error )
    allocate( new_data_w_ghost(Bs+2*g, Bs+2*g), stat=allocate_error )

    allocate( data_predict_fine(2*Bs-1, 2*Bs-1), stat=allocate_error )
    allocate( data_predict_coarse(Bs, Bs), stat=allocate_error )

    allocate( new_coord_x(Bs), stat=allocate_error )
    allocate( new_coord_y(Bs), stat=allocate_error )

    allocate( send_receive_data(Bs, Bs), stat=allocate_error )
    allocate( send_receive_coord(Bs), stat=allocate_error )

    ! loop over all data fields
    do dF = 1, blocks_params%number_data_fields

        ! loop over all blocks
        do k = 1, N

        if ( blocks(k)%active ) then

            ! current block, here the new data will be stored
            ! heavy block id
            block_num = blocks(k)%proc_data_id
            ! light block id
            light_id = k
            ! current proc rank
            data_rank = blocks(k)%proc_rank

            if (blocks(light_id)%refinement == -2) then
                ! coarsening

                ! get treecodes for current block and sister-blocks
                me                                              = blocks(light_id)%treecode
                block_ids( me( blocks(light_id)%level ) + 1 )   = light_id
                proc_rank( me( blocks(light_id)%level ) + 1 )   = blocks(light_id)%proc_rank

                call get_sister_id(id1, id2, id3, me, blocks(light_id)%level)

                s1                                              = blocks(id1)%treecode
                block_ids( s1( blocks(id1)%level ) + 1 )        = id1
                proc_rank( s1( blocks(id1)%level ) + 1 )        = blocks(id1)%proc_rank

                s2                                              = blocks(id2)%treecode
                block_ids( s2( blocks(id2)%level ) + 1 )        = id2
                proc_rank( s2( blocks(id2)%level ) + 1 )        = blocks(id2)%proc_rank

                s3                                              = blocks(id3)%treecode
                block_ids( s3( blocks(id3)%level ) + 1 )        = id3
                proc_rank( s3( blocks(id3)%level ) + 1 )        = blocks(id3)%proc_rank

                if ( data_rank == rank ) then
                    ! creating new block data
                    new_data                                        = 9.0e9_rk
                    new_data_w_ghost                                = 9.0e9_rk
                    ! coordinates for new block
                    new_coord_x                                     = 9.0e9_rk
                    new_coord_y                                     = 9.0e9_rk
                end if

                ! loop over proc rank list, proc with light id block data collect new data and coarse block
                ! then save new data on light id, current block_num
                                do i = 1, 4

                    if ( data_rank == rank ) then
                        ! collect data
                        if ( proc_rank(i) == rank ) then
                            ! no communication needed
                            select case(i)
                                case(1)
                                    ! sister 0
                                    new_data(1:(Bs-1)/2+1, 1:(Bs-1)/2+1)   = blocks_data( blocks( block_ids(1) )%proc_data_id )%data_fields(dF)%data_(g+1:Bs+g:2, g+1:Bs+g:2)
                                    new_coord_x(1:(Bs-1)/2+1)              = blocks_data( blocks( block_ids(1) )%proc_data_id )%coord_x(1:Bs:2)
                                    new_coord_y(1:(Bs-1)/2+1)              = blocks_data( blocks( block_ids(1) )%proc_data_id )%coord_y(1:Bs:2)
                                case(2)
                                    ! sister 1
                                    new_data(1:(Bs-1)/2+1, (Bs-1)/2+1:Bs)  = blocks_data( blocks( block_ids(2) )%proc_data_id )%data_fields(dF)%data_(g+1:Bs+g:2, g+1:Bs+g:2)
                                    new_coord_x((Bs-1)/2+1:Bs)             = blocks_data( blocks( block_ids(2) )%proc_data_id )%coord_x(1:Bs:2)
                                case(3)
                                    ! sister 2
                                    new_data((Bs-1)/2+1:Bs, 1:(Bs-1)/2+1)  = blocks_data( blocks( block_ids(3) )%proc_data_id )%data_fields(dF)%data_(g+1:Bs+g:2, g+1:Bs+g:2)
                                    new_coord_y((Bs-1)/2+1:Bs)             = blocks_data( blocks( block_ids(3) )%proc_data_id )%coord_y(1:Bs:2)
                                case(4)
                                    ! sister 3
                                    new_data((Bs-1)/2+1:Bs, (Bs-1)/2+1:Bs) = blocks_data( blocks( block_ids(4) )%proc_data_id )%data_fields(dF)%data_(g+1:Bs+g:2, g+1:Bs+g:2)
                            end select
                        else
                            ! receive data from other proc, note: not all blocks have to send coord vectors
                            select case(i)
                                case(1)
                                    ! sister 0
                                    ! receive data
                                    call MPI_Recv(send_receive_data, Bs*Bs, MPI_REAL8, proc_rank(i), tag, MPI_COMM_WORLD, status, ierr)
                                    new_data(1:(Bs-1)/2+1, 1:(Bs-1)/2+1)   = send_receive_data(1:Bs:2, 1:Bs:2)
                                    ! receive coords
                                    call MPI_Recv(send_receive_coord, Bs, MPI_REAL8, proc_rank(i), tag, MPI_COMM_WORLD, status, ierr)
                                    new_coord_x(1:(Bs-1)/2+1)              = send_receive_coord(1:Bs:2)
                                    call MPI_Recv(send_receive_coord, Bs, MPI_REAL8, proc_rank(i), tag, MPI_COMM_WORLD, status, ierr)
                                    new_coord_y(1:(Bs-1)/2+1)              = send_receive_coord(1:Bs:2)
                                case(2)
                                    ! sister 1
                                    ! receive data
                                    call MPI_Recv(send_receive_data, Bs*Bs, MPI_REAL8, proc_rank(i), tag, MPI_COMM_WORLD, status, ierr)
                                    new_data(1:(Bs-1)/2+1, (Bs-1)/2+1:Bs)  = send_receive_data(1:Bs:2, 1:Bs:2)
                                    ! receive coords
                                    call MPI_Recv(send_receive_coord, Bs, MPI_REAL8, proc_rank(i), tag, MPI_COMM_WORLD, status, ierr)
                                    new_coord_x((Bs-1)/2+1:Bs)             = send_receive_coord(1:Bs:2)
                                case(3)
                                    ! sister 2
                                    ! receive data
                                    call MPI_Recv(send_receive_data, Bs*Bs, MPI_REAL8, proc_rank(i), tag, MPI_COMM_WORLD, status, ierr)
                                    new_data((Bs-1)/2+1:Bs, 1:(Bs-1)/2+1)  = send_receive_data(1:Bs:2, 1:Bs:2)
                                    ! receive coords
                                    call MPI_Recv(send_receive_coord, Bs, MPI_REAL8, proc_rank(i), tag, MPI_COMM_WORLD, status, ierr)
                                    new_coord_y((Bs-1)/2+1:Bs)              = send_receive_coord(1:Bs:2)
                                case(4)
                                    ! sister 3
                                    ! receive data
                                    call MPI_Recv(send_receive_data, Bs*Bs, MPI_REAL8, proc_rank(i), tag, MPI_COMM_WORLD, status, ierr)
                                    new_data((Bs-1)/2+1:Bs, (Bs-1)/2+1:Bs) = send_receive_data(1:Bs:2, 1:Bs:2)
                            end select
                        end if

                    elseif ( proc_rank(i) == rank ) then
                        ! send data
                        select case(i)
                            case(1)
                                ! sister 0
                                ! send data
                                send_receive_data = blocks_data( blocks( block_ids(1) )%proc_data_id )%data_fields(dF)%data_(g+1:Bs+g, g+1:Bs+g)
                                call MPI_Send( send_receive_data, Bs*Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                                ! send coord
                                send_receive_coord = blocks_data( blocks( block_ids(1) )%proc_data_id )%coord_x(1:Bs)
                                call MPI_Send( send_receive_coord, Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                                send_receive_coord = blocks_data( blocks( block_ids(1) )%proc_data_id )%coord_y(1:Bs)
                                call MPI_Send( send_receive_coord, Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                            case(2)
                                ! sister 1
                                ! send data
                                send_receive_data = blocks_data( blocks( block_ids(2) )%proc_data_id )%data_fields(dF)%data_(g+1:Bs+g, g+1:Bs+g)
                                call MPI_Send( send_receive_data, Bs*Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                                ! send coord
                                send_receive_coord = blocks_data( blocks( block_ids(2) )%proc_data_id )%coord_x(1:Bs)
                                call MPI_Send( send_receive_coord, Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                            case(3)
                                ! sister 2
                                ! send data
                                send_receive_data = blocks_data( blocks( block_ids(3) )%proc_data_id )%data_fields(dF)%data_(g+1:Bs+g, g+1:Bs+g)
                                call MPI_Send( send_receive_data, Bs*Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                                ! send coord
                                send_receive_coord = blocks_data( blocks( block_ids(3) )%proc_data_id )%coord_y(1:Bs)
                                call MPI_Send( send_receive_coord, Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                            case(4)
                                ! sister 3
                                ! send data
                                send_receive_data = blocks_data( blocks( block_ids(4) )%proc_data_id )%data_fields(dF)%data_(g+1:Bs+g, g+1:Bs+g)
                                call MPI_Send( send_receive_data, Bs*Bs, MPI_REAL8, data_rank, tag, MPI_COMM_WORLD, ierr)
                        end select
                    else
                        ! nothing to do
                    end if

                end do

                ! delete all heavy data
                if ( proc_rank(1) == rank ) call delete_block_heavy( blocks( block_ids(1) )%proc_data_id )
                if ( proc_rank(2) == rank ) call delete_block_heavy( blocks( block_ids(2) )%proc_data_id )
                if ( proc_rank(3) == rank ) call delete_block_heavy( blocks( block_ids(3) )%proc_data_id )
                if ( proc_rank(4) == rank ) call delete_block_heavy( blocks( block_ids(4) )%proc_data_id )

                ! new treecode, one level down (coarsening)
                new_treecode                                    = me
                new_treecode(blocks(light_id)%level)           = -1

                ! write new block on current block
                ! write light data
                call new_block_light(light_id, new_treecode, 10)
                ! set proc_rank and id
                blocks(light_id)%proc_rank = data_rank
                blocks(light_id)%proc_data_id = block_num

                ! final steps, only for sister0-proc
                if ( data_rank == rank ) then

                    ! write block data to data field with ghost nodes
                    new_data_w_ghost(g+1:Bs+g, g+1:Bs+g)            = new_data

                    ! write heavy data
                    call new_block_heavy(block_num, light_id, new_data_w_ghost, new_coord_x, new_coord_y, Bs, g, dF)

                end if

                ! delete light data
                call delete_block_light(id1)
                call delete_block_light(id2)
                call delete_block_light(id3)

            elseif (blocks(light_id)%refinement == 1) then
                ! refinement
                ! new block data, new data, only for proc corresponding block
                if ( blocks(light_id)%proc_rank == rank ) then
                    new_data                                        = 9.0e9_rk
                    new_data_w_ghost                                = 9.0e9_rk
                    data_predict_fine                               = 9.0e9_rk
                    data_predict_coarse                             = 9.0e9_rk

                    ! data prediction
                    data_predict_coarse                             = blocks_data(block_num)%data_fields(dF)%data_(g+1:Bs+g, g+1:Bs+g)
                    call prediction_2D(data_predict_coarse, data_predict_fine)

                end if

                !--------------------------
                ! first new block
                ! find free block id
                light_free_id                                   = 0
                call get_light_free_block(light_free_id)

                ! create new treecode
                new_treecode                                    = blocks(light_id)%treecode
                new_treecode( blocks(light_id)%level + 1 )      = 0

                ! new light block data
                call new_block_light(light_free_id, new_treecode, 10)

                ! new heavy block data, proc corresponding to block do the work
                if ( blocks(light_id)%proc_rank == rank ) then

                    new_data            = data_predict_fine(1:Bs, 1:Bs)

                    ! new coordinates
                    new_coord_x(1:Bs:2) = blocks_data(block_num)%coord_x(1:(Bs-1)/2+1)
                    new_coord_y(1:Bs:2) = blocks_data(block_num)%coord_y(1:(Bs-1)/2+1)

                    do i = 2, Bs, 2
                        new_coord_x(i)  = ( new_coord_x(i-1) + new_coord_x(i+1) ) / 2.0_rk
                        new_coord_y(i)  = ( new_coord_y(i-1) + new_coord_y(i+1) ) / 2.0_rk
                    end do

                    ! write block data to data field with ghost nodes
                    new_data_w_ghost(g+1:Bs+g, g+1:Bs+g)            = new_data

                    ! find free local block id
                    call get_heavy_free_block(heavy_free_id)

                    ! write new block on free id
                    call new_block_heavy(heavy_free_id, light_free_id, new_data_w_ghost, new_coord_x, new_coord_y, Bs, g, dF)

                end if

                ! synchronize new heavy block id
                call MPI_Bcast(heavy_free_id, 1, MPI_INTEGER4, blocks(light_id)%proc_rank, MPI_COMM_WORLD, ierr)

                ! set proc_rank and id
                blocks(light_free_id)%proc_rank = blocks(light_id)%proc_rank
                blocks(light_free_id)%proc_data_id = heavy_free_id

                ! reset refinement status
                blocks(light_free_id)%refinement = 0

                !--------------------------
                ! second new block
                ! find free block id
                light_free_id                                   = 0
                call get_light_free_block(light_free_id)

                ! create new treecode
                new_treecode                                    = blocks(light_id)%treecode
                new_treecode( blocks(light_id)%level + 1 )      = 1

                ! new light block data
                call new_block_light(light_free_id, new_treecode, 10)

                ! new heavy block data, proc corresponding to block do the work
                if ( blocks(light_id)%proc_rank == rank ) then

                    new_data            = data_predict_fine(1:Bs, Bs:2*Bs-1)

                    ! new coordinates
                    new_coord_x(1:Bs:2) = blocks_data(block_num)%coord_x((Bs-1)/2+1:Bs)
                    new_coord_y(1:Bs:2) = blocks_data(block_num)%coord_y(1:(Bs-1)/2+1)

                    do i = 2, Bs, 2
                        new_coord_x(i)  = ( new_coord_x(i-1) + new_coord_x(i+1) ) / 2.0_rk
                        new_coord_y(i)  = ( new_coord_y(i-1) + new_coord_y(i+1) ) / 2.0_rk
                    end do

                    ! write block data to data field with ghost nodes
                    new_data_w_ghost(g+1:Bs+g, g+1:Bs+g)            = new_data

                    ! find free local block id
                    call get_heavy_free_block(heavy_free_id)

                    ! write new block on free id
                    call new_block_heavy(heavy_free_id, light_free_id, new_data_w_ghost, new_coord_x, new_coord_y, Bs, g, dF)

                end if

                ! synchronize new heavy block id
                call MPI_Bcast(heavy_free_id, 1, MPI_INTEGER4, blocks(light_id)%proc_rank, MPI_COMM_WORLD, ierr)

                ! set proc_rank and id
                blocks(light_free_id)%proc_rank = blocks(light_id)%proc_rank
                blocks(light_free_id)%proc_data_id = heavy_free_id

                ! reset refinement status
                blocks(light_free_id)%refinement = 0

                !--------------------------
                ! third new block
                ! find free block id
                light_free_id                                   = 0
                call get_light_free_block(light_free_id)

                ! create new treecode
                new_treecode                                    = blocks(light_id)%treecode
                new_treecode( blocks(light_id)%level + 1 )      = 2

                ! new light block data
                call new_block_light(light_free_id, new_treecode, 10)

                ! new heavy block data, proc corresponding to block do the work
                if ( blocks(light_id)%proc_rank == rank ) then

                    new_data            = data_predict_fine(Bs:2*Bs-1, 1:Bs)

                    ! new coordinates
                    new_coord_x(1:Bs:2) = blocks_data(block_num)%coord_x(1:(Bs-1)/2+1)
                    new_coord_y(1:Bs:2) = blocks_data(block_num)%coord_y((Bs-1)/2+1:Bs)

                    do i = 2, Bs, 2
                        new_coord_x(i)  = ( new_coord_x(i-1) + new_coord_x(i+1) ) / 2.0_rk
                        new_coord_y(i)  = ( new_coord_y(i-1) + new_coord_y(i+1) ) / 2.0_rk
                    end do

                    ! write block data to data field with ghost nodes
                    new_data_w_ghost(g+1:Bs+g, g+1:Bs+g)            = new_data

                    ! find free local block id
                    call get_heavy_free_block(heavy_free_id)

                    ! write new block on free id
                    call new_block_heavy(heavy_free_id, light_free_id, new_data_w_ghost, new_coord_x, new_coord_y, Bs, g, dF)

                end if

                ! synchronize new heavy block id
                call MPI_Bcast(heavy_free_id, 1, MPI_INTEGER4, blocks(light_id)%proc_rank, MPI_COMM_WORLD, ierr)

                ! set proc_rank and id
                blocks(light_free_id)%proc_rank = blocks(light_id)%proc_rank
                blocks(light_free_id)%proc_data_id = heavy_free_id

                ! reset refinement status
                blocks(light_free_id)%refinement = 0

                !--------------------------
                ! fourth new block
                ! write light data on old block

                ! create new treecode
                new_treecode                                    = blocks(light_id)%treecode
                new_treecode( blocks(light_id)%level + 1 )      = 3

                ! new light block data
                call new_block_light(light_id, new_treecode, 10)

                ! new heavy block data, proc corresponding to block do the work
                if ( blocks(light_id)%proc_rank == rank ) then

                    new_data            = data_predict_fine(Bs:2*Bs-1, Bs:2*Bs-1)

                    ! new coordinates
                    new_coord_x(1:Bs:2) = blocks_data(block_num)%coord_x((Bs-1)/2+1:Bs)
                    new_coord_y(1:Bs:2) = blocks_data(block_num)%coord_y((Bs-1)/2+1:Bs)

                    do i = 2, Bs, 2
                        new_coord_x(i)  = ( new_coord_x(i-1) + new_coord_x(i+1) ) / 2.0_rk
                        new_coord_y(i)  = ( new_coord_y(i-1) + new_coord_y(i+1) ) / 2.0_rk
                    end do

                    ! write new block data on old block
                    ! first clear old data
                    call delete_block_heavy(block_num)

                    ! write block data to data field with ghost nodes
                    new_data_w_ghost(g+1:Bs+g, g+1:Bs+g)            = new_data

                    ! find free local block id
                    call get_heavy_free_block(heavy_free_id)

                    ! write new block on free id
                    call new_block_heavy(heavy_free_id, light_id, new_data_w_ghost, new_coord_x, new_coord_y, Bs, g, dF)

                end if

                ! synchronize new heavy block id
                call MPI_Bcast(heavy_free_id, 1, MPI_INTEGER4, blocks(light_id)%proc_rank, MPI_COMM_WORLD, ierr)

                ! set proc heavy id
                blocks(light_id)%proc_data_id   = heavy_free_id

                ! reset refinement status
                blocks(light_id)%refinement = 0

            end if

         end if

        end do

    end do

    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    deallocate( new_data, stat=allocate_error )
    deallocate( new_data_w_ghost, stat=allocate_error )
    deallocate( data_predict_fine, stat=allocate_error )
    deallocate( data_predict_coarse, stat=allocate_error )
    deallocate( new_coord_x, stat=allocate_error )
    deallocate( new_coord_y, stat=allocate_error )
    deallocate( send_receive_data, stat=allocate_error )
    deallocate( send_receive_coord, stat=allocate_error )

end subroutine interpolate_mesh
