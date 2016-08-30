subroutine  allocate_block_memory()
    use module_params
    use module_blocks

    implicit none

    integer :: i, allocate_error
    integer(kind=ik) :: blocksize, ghosts

    blocksize = blocks_params%size_block
    ghosts = blocks_params%number_ghost_nodes
    write(*,'("System is allocating",i7," blocks" )') blocks_params%number_max_blocks

    allocate( blocks(1:blocks_params%number_max_blocks), stat=allocate_error )
    ! dummy allocation
    allocate( blocks_params%active_list(1), stat=allocate_error )

    do i = 1, blocks_params%number_max_blocks
        allocate( blocks(i)%data1(blocksize, blocksize), stat=allocate_error )
        allocate( blocks(i)%data2(blocksize+2*ghosts, blocksize+2*ghosts), stat=allocate_error )
        allocate( blocks(i)%data_old(blocksize+2*ghosts, blocksize+2*ghosts), stat=allocate_error )
        allocate( blocks(i)%k1(blocksize+2*ghosts, blocksize+2*ghosts), stat=allocate_error )
        allocate( blocks(i)%k2(blocksize+2*ghosts, blocksize+2*ghosts), stat=allocate_error )
        allocate( blocks(i)%k3(blocksize+2*ghosts, blocksize+2*ghosts), stat=allocate_error )
        allocate( blocks(i)%k4(blocksize+2*ghosts, blocksize+2*ghosts), stat=allocate_error )
        allocate( blocks(i)%treecode(10), stat=allocate_error )
        allocate( blocks(i)%neighbor_treecode(8,10), stat=allocate_error )
        allocate( blocks(i)%neighbor2_treecode(4,10), stat=allocate_error )

        blocks(i)%data1(:,:)                = 0.0_rk
        blocks(i)%data2(:,:)                = 0.0_rk
        blocks(i)%data_old(:,:)             = 0.0_rk
        blocks(i)%k1(:,:)                   = 0.0_rk
        blocks(i)%k2(:,:)                   = 0.0_rk
        blocks(i)%k3(:,:)                   = 0.0_rk
        blocks(i)%k4(:,:)                   = 0.0_rk
        blocks(i)%active                    = .false.
        blocks(i)%treecode(:)               = -1
        blocks(i)%neighbor_treecode(:,:)    = -1
        blocks(i)%neighbor2_treecode(:,:)   = -1
        blocks(i)%neighbor_id(:)            = -1
        blocks(i)%neighbor2_id(:)           = -1
        blocks(i)%refinement                = 0
        blocks(i)%level                     = -1
        blocks(i)%neighbor_number(:)        = 1

        blocks(i)%dx = params%Lx / real(blocks_params%size_domain,8)
        blocks(i)%dy = params%Lx / real(blocks_params%size_domain,8)

        ! coordinates
        allocate( blocks(i)%coord_x(blocksize) )
        allocate( blocks(i)%coord_y(blocksize) )

        blocks(i)%coord_x(:)    = 0.0_rk
        blocks(i)%coord_y(:)    = 0.0_rk
    end do
end subroutine
