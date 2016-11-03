! ********************************
! WABBIT
! --------------------------------
!
! create block tree from domain at start
!
! name: matrix_to_block_tree.f90
! date: 25.10.2016
! author: msr
! version: 0.3
!
! ********************************

subroutine matrix_to_block_tree()

    use mpi
    use module_params
    use module_blocks

    implicit none

    integer(kind=ik)                                :: num_blocks_x, num_blocks_y, num_blocks
    integer(kind=ik)                                :: Bs, Ds, g, dF, i, j, ib, allocate_error, rank, ierr
    integer(kind=ik), dimension(10)                 :: treecode

    real(kind=rk), dimension(:,:), allocatable      :: data_
    real(kind=rk), dimension(:), allocatable        :: coord_x, coord_y, domain_coord_x, domain_coord_y

    character(len=80)                               :: distribution

    Bs                  = blocks_params%size_block
    g                   = blocks_params%number_ghost_nodes
    Ds                  = blocks_params%size_domain
    dF                  = blocks_params%number_data_fields

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    ! allocate memory for data fields and coordinates vectors
    allocate( data_(Bs+2*g, Bs+2*g), stat=allocate_error )

    allocate( coord_x(Bs), stat=allocate_error )
    allocate( coord_y(Bs), stat=allocate_error )
    allocate( domain_coord_x(Ds), stat=allocate_error )
    allocate( domain_coord_y(Ds), stat=allocate_error )

    ! assign -1 to treecode
    treecode            = -1

    ! print block decomposition information
    ! every block has two more points in a single direction (from his neighbors)
    ! therefore the complete domain has also two additional points
    num_blocks_x        = (Ds-1) / (Bs-1)
    num_blocks_y        = (Ds-1) / (Bs-1)
    num_blocks          = num_blocks_x * num_blocks_y

    ! check given domain and block size
    if ( Ds /= (num_blocks_x-1)*(Bs-1) + Bs ) then
        print*, "error: blocksize do not fit into domain size"
        stop
    end if

    if (rank==0) then
        write(*,'(a,i5,a,i5,a,i5,a,i5,a,i5,a,i5)') "Field with res: ", (Ds), " x", (Ds), " gives: ", &
        num_blocks_x, " x", num_blocks_y, " (", num_blocks, ") blocks of size: ", (Bs)
        write(*,'(80("-"))')
    end if

    ! domain coordinate vectors
    do i = 1, Ds
        domain_coord_x(i) = (i-1) * params%Lx / (Ds-1)
        domain_coord_y(i) = params%Lx - (i-1) * params%Ly / (Ds-1)
    end do

    ! create new blocks, first: only light data

    ib = 1
    do i = 1, num_blocks_x
        do j = 1, num_blocks_y
            ! treecode
            call encoding(treecode, i, j, num_blocks_x, num_blocks_y)
            ! write light data
            call new_block_light(ib, treecode, 10)
            ! block counter
            ib                          = ib + 1
        end do
    end do

    ! distribute blocks to procs
    distribution = "equal"
    call initial_block_distribution(distribution)

    ! create new blocks, second: heavy data
    if (dF == 1) then

        ! exactly one field
        ib = 1
        do i = 1, num_blocks_x
            do j = 1, num_blocks_y
                ! coordinates
                coord_x                     = domain_coord_x( (j-1)*(Bs-1) + 1 : j*(Bs-1) + 1 )
                coord_y                     = domain_coord_y( (i-1)*(Bs-1) + 1 : i*(Bs-1) + 1 )
                ! data
                data_(:,:)                  = 0.0_rk
                data_(g+1:Bs+g,g+1:Bs+g)    = blocks_params%phi( (i-1)*(Bs-1) + 1 : i*(Bs-1) + 1 , (j-1)*(Bs-1) + 1 : j*(Bs-1) + 1 )

                ! write heavy data
                if ( rank == blocks(ib)%proc_rank ) then
                    call new_block_heavy(blocks(ib)%proc_data_id, ib, data_, coord_x, coord_y, Bs, g, dF)
                end if

                ! block counter
                ib                          = ib + 1
            end do
        end do

    else
        ! more than one field
        ! write block data, write start field (phi) data to datafield 1, set all other datafields to zero
        dF = 1
        ib = 1
        do i = 1, num_blocks_x
            do j = 1, num_blocks_y
                ! coordinates
                coord_x                     = domain_coord_x( (j-1)*(Bs-1) + 1 : j*(Bs-1) + 1 )
                coord_y                     = domain_coord_y( (i-1)*(Bs-1) + 1 : i*(Bs-1) + 1 )
                ! data
                data_(:,:)                  = 0.0_rk
                data_(g+1:Bs+g,g+1:Bs+g)    = blocks_params%phi( (i-1)*(Bs-1) + 1 : i*(Bs-1) + 1 , (j-1)*(Bs-1) + 1 : j*(Bs-1) + 1 )

                ! write heavy data
                if ( rank == blocks(ib)%proc_rank ) then
                    call new_block_heavy(blocks(ib)%proc_data_id, ib, data_, coord_x, coord_y, Bs, g, dF)
                end if

                ! block counter
                ib                          = ib + 1
            end do
        end do

        ! set all other datafields to zero
        data_(:,:)                  = 0.0_rk
        do dF = 2, blocks_params%number_data_fields
            ! loop over all heavy data
            do i = 1, blocks_params%number_max_blocks_data
                if ( blocks_data(i)%block_id /= -1 ) then

                    ! block is active
                    call set_heavy_data(i, data_, Bs, g, dF)

                end if
            end do
        end do

    end if

    ! deallocate memory for local variables
    deallocate( data_, stat=allocate_error )
    deallocate( coord_x, stat=allocate_error )
    deallocate( coord_y, stat=allocate_error )
    deallocate( domain_coord_x, stat=allocate_error )
    deallocate( domain_coord_x, stat=allocate_error )

end subroutine matrix_to_block_tree
