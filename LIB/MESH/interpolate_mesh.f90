! ********************************
! 2D AMR prototype
! --------------------------------
!
! adapt the mesh with interpolation
!
! name: interpolate_mesh.f90
! date: 18.08.2016
! author: msr
! version: 0.1
!
! ********************************

subroutine interpolate_mesh()

    use module_params
    use module_blocks
    use module_interpolation

    implicit none

    integer(kind=ik)                                :: k, N, block_num, id1, id2, id3, allocate_error, Bs, g, i, free_id
    integer(kind=ik), dimension(10)                 :: me, s1, s2, s3, new_treecode
    integer(kind=ik), dimension(4)                  :: block_ids

    real(kind=rk), dimension(:,:), allocatable      :: new_data, data_predict_fine, data_predict_coarse
    real(kind=rk), dimension(:), allocatable        :: new_coord_x, new_coord_y

    N  = size(blocks_params%active_list, dim=1)
    Bs = blocks_params%size_block
    g  = blocks_params%number_ghost_nodes

    allocate( new_data(Bs, Bs), stat=allocate_error )
    allocate( data_predict_fine(2*Bs-1, 2*Bs-1), stat=allocate_error )
    allocate( data_predict_coarse(Bs, Bs), stat=allocate_error )

    allocate( new_coord_x(Bs), stat=allocate_error )
    allocate( new_coord_y(Bs), stat=allocate_error )

    ! loop over all blocks
    do k = 1, blocks_params%number_max_blocks

    if (blocks(k)%active) then

        block_num = k

        if (blocks(block_num)%refinement == -2) then
            ! coarsening

            ! get treecodes for current block and sister-blocks
            me                                              = blocks(block_num)%treecode
            block_ids( me( blocks(block_num)%level ) + 1 )  = block_num

            call get_sister_id(id1, id2, id3, me, blocks(block_num)%level)

            s1                                              = blocks(id1)%treecode
            block_ids( s1( blocks(id1)%level ) + 1 )        = id1

            s2                                              = blocks(id2)%treecode
            block_ids( s2( blocks(id2)%level ) + 1 )        = id2

            s3                                              = blocks(id3)%treecode
            block_ids( s3( blocks(id3)%level ) + 1 )        = id3

            ! creating new block data
            new_data                                        = 0.0_rk
            new_data(1:(Bs-1)/2+1, 1:(Bs-1)/2+1)            = blocks( block_ids(1) )%data1(1:Bs:2, 1:Bs:2)
            new_data(1:(Bs-1)/2+1, (Bs-1)/2+1:Bs)           = blocks( block_ids(2) )%data1(1:Bs:2, 1:Bs:2)
            new_data((Bs-1)/2+1:Bs, 1:(Bs-1)/2+1)           = blocks( block_ids(3) )%data1(1:Bs:2, 1:Bs:2)
            new_data((Bs-1)/2+1:Bs, (Bs-1)/2+1:Bs)          = blocks( block_ids(4) )%data1(1:Bs:2, 1:Bs:2)

            ! coordinates for new block
            new_coord_x                                     = 0.0_rk
            new_coord_y                                     = 0.0_rk
            new_coord_x(1:(Bs-1)/2+1)                       = blocks( block_ids(1) )%coord_x(1:Bs:2)
            new_coord_x((Bs-1)/2+1:Bs)                      = blocks( block_ids(2) )%coord_x(1:Bs:2)
            new_coord_y(1:(Bs-1)/2+1)                       = blocks( block_ids(1) )%coord_y(1:Bs:2)
            new_coord_y((Bs-1)/2+1:Bs)                      = blocks( block_ids(3) )%coord_y(1:Bs:2)

            ! new treecode, one level down (coarsening)
            new_treecode                                    = me
            new_treecode(blocks(block_num)%level)           = -1

            ! write new block on actual block
            call new_block(block_num, new_treecode, 10, new_data, new_coord_x, new_coord_y, Bs)

            ! delete sister blocks
            call delete_block(id1)
            call delete_block(id2)
            call delete_block(id3)

        elseif (blocks(block_num)%refinement == 1) then
            ! refinement

            ! new block data, assume synchronized ghost nodes
            new_data                                        = 0.0_rk
            data_predict_fine                               = 0.0_rk
            data_predict_coarse                             = 0.0_rk

            data_predict_coarse                             = blocks(block_num)%data1
            call prediction_2D(data_predict_coarse, data_predict_fine)

            !--------------------------
            ! first new block
            ! find free block id
            free_id                                         = 0
            call get_free_block(free_id)

            ! create new treecode
            new_treecode                                    = blocks(block_num)%treecode
            new_treecode( blocks(block_num)%level + 1 )     = 0

            ! new block data
            new_data                                        = data_predict_fine(1:Bs, 1:Bs)

            ! new coordinates
            new_coord_x(1:Bs:2)                             = blocks(block_num)%coord_x(1:(Bs-1)/2+1)
            new_coord_y(1:Bs:2)                             = blocks(block_num)%coord_y(1:(Bs-1)/2+1)

            do i = 2, Bs, 2
                new_coord_x(i)                              = ( new_coord_x(i-1) + new_coord_x(i+1) ) / 2.0_rk
                new_coord_y(i)                              = ( new_coord_y(i-1) + new_coord_y(i+1) ) / 2.0_rk
            end do

            ! write new block
            call new_block(free_id, new_treecode, 10, new_data, new_coord_x, new_coord_y, Bs)

            !--------------------------
            ! second new block
            ! find free block id
            free_id                                         = 0
            call get_free_block(free_id)

            ! create new treecode
            new_treecode                                    = blocks(block_num)%treecode
            new_treecode( blocks(block_num)%level + 1 )     = 1

            ! new block data
            new_data                                        = data_predict_fine(1:Bs, Bs:2*Bs-1)

            ! new coordinates
            new_coord_x(1:Bs:2)                             = blocks(block_num)%coord_x((Bs-1)/2+1:Bs)
            new_coord_y(1:Bs:2)                             = blocks(block_num)%coord_y(1:(Bs-1)/2+1)

            do i = 2, Bs, 2
                new_coord_x(i)                              = ( new_coord_x(i-1) + new_coord_x(i+1) ) / 2.0_rk
                new_coord_y(i)                              = ( new_coord_y(i-1) + new_coord_y(i+1) ) / 2.0_rk
            end do

            ! write new block
            call new_block(free_id, new_treecode, 10, new_data, new_coord_x, new_coord_y, Bs)

            !--------------------------
            ! third new block
            ! find free block id
            free_id                                         = 0
            call get_free_block(free_id)

            ! create new treecode
            new_treecode                                    = blocks(block_num)%treecode
            new_treecode( blocks(block_num)%level + 1 )     = 2

            ! new block data
            new_data                                        = data_predict_fine(Bs:2*Bs-1, 1:Bs)

            ! new coordinates
            new_coord_x(1:Bs:2)                             = blocks(block_num)%coord_x(1:(Bs-1)/2+1)
            new_coord_y(1:Bs:2)                             = blocks(block_num)%coord_y((Bs-1)/2+1:Bs)

            do i = 2, Bs, 2
                new_coord_x(i)                              = ( new_coord_x(i-1) + new_coord_x(i+1) ) / 2.0_rk
                new_coord_y(i)                              = ( new_coord_y(i-1) + new_coord_y(i+1) ) / 2.0_rk
            end do

            ! write new block
            call new_block(free_id, new_treecode, 10, new_data, new_coord_x, new_coord_y, Bs)

            !--------------------------
            ! fourth new block
            ! find free block id
            free_id                                         = 0

            ! create new treecode
            new_treecode                                    = blocks(block_num)%treecode
            new_treecode( blocks(block_num)%level + 1 )     = 3

            ! new block data
            new_data                                        = data_predict_fine(Bs:2*Bs-1, Bs:2*Bs-1)

            ! new coordinates
            new_coord_x(1:Bs:2)                             = blocks(block_num)%coord_x((Bs-1)/2+1:Bs)
            new_coord_y(1:Bs:2)                             = blocks(block_num)%coord_y((Bs-1)/2+1:Bs)

            do i = 2, Bs, 2
                new_coord_x(i)                              = ( new_coord_x(i-1) + new_coord_x(i+1) ) / 2.0_rk
                new_coord_y(i)                              = ( new_coord_y(i-1) + new_coord_y(i+1) ) / 2.0_rk
            end do

            ! write new block data on old block
            ! first clear old data
            call delete_block(block_num)
            ! should be lock_num
            call get_free_block(free_id)
            ! writing
            call new_block(free_id, new_treecode, 10, new_data, new_coord_x, new_coord_y, Bs)

        end if

       ! reset refinement status
       blocks(block_num)%refinement = 0

    end if

    end do

    deallocate( new_data, stat=allocate_error )
    deallocate( data_predict_fine, stat=allocate_error )
    deallocate( data_predict_coarse, stat=allocate_error )
    deallocate( new_coord_x, stat=allocate_error )
    deallocate( new_coord_y, stat=allocate_error )

end subroutine interpolate_mesh
