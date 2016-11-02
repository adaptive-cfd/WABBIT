! ********************************
! WABBIT
! --------------------------------
!
! allocate memory for block data
! and initialize all data
!
! name: allocate_block_memory.f90
! date: 25.10.2016
! author: engels, msr
! version: 0.3
!
! ********************************

subroutine  allocate_block_memory()

    use mpi
    use module_params
    use module_blocks

    implicit none

    integer                                         :: i, j, allocate_error
    integer(kind=ik)                                :: Bs, g, dF, rank, ierr

    ! grid size parameters
    Bs  = blocks_params%size_block
    g   = blocks_params%number_ghost_nodes
    dF  = blocks_params%number_data_fields

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    ! output
    if (rank==0) then
        write(*,'("System is allocating",i7," blocks for light data" )') blocks_params%number_max_blocks
        write(*,'("System is allocating",i7," blocks for heavy data" )') blocks_params%number_max_blocks_data
        write(*,'(80("-"))')
    end if

    ! allocate light data
    allocate( blocks(1:blocks_params%number_max_blocks), stat=allocate_error )

    do i = 1, blocks_params%number_max_blocks

        allocate( blocks(i)%treecode(10), stat=allocate_error )
        allocate( blocks(i)%neighbor_treecode(16,10), stat=allocate_error )

        blocks(i)%active                    = .false.
        blocks(i)%treecode(:)               = -1
        blocks(i)%neighbor_treecode(:,:)    = -1
        blocks(i)%neighbor_id(:)            = -1
        blocks(i)%refinement                = 0
        blocks(i)%level                     = -1

        blocks(i)%proc_rank                 = -1
        blocks(i)%proc_data_id              = -1

    end do

    ! allocate heavy data
    allocate( blocks_data(1:blocks_params%number_max_blocks_data), stat=allocate_error )

    do i = 1, blocks_params%number_max_blocks_data

        allocate( blocks_data(i)%data_fields(1:dF), stat=allocate_error )

        ! loop over all data fields
        do j = 1, dF

            allocate( blocks_data(i)%data_fields(j)%data_(Bs+2*g, Bs+2*g), stat=allocate_error )
            allocate( blocks_data(i)%data_fields(j)%data_old(Bs+2*g, Bs+2*g), stat=allocate_error )

            allocate( blocks_data(i)%data_fields(j)%k1(Bs+2*g, Bs+2*g), stat=allocate_error )
            allocate( blocks_data(i)%data_fields(j)%k2(Bs+2*g, Bs+2*g), stat=allocate_error )
            allocate( blocks_data(i)%data_fields(j)%k3(Bs+2*g, Bs+2*g), stat=allocate_error )
            allocate( blocks_data(i)%data_fields(j)%k4(Bs+2*g, Bs+2*g), stat=allocate_error )

        end do

        ! loop over all data fields
        do j = 1, dF

            blocks_data(i)%data_fields(j)%data_(:,:)              = 9.0e9_rk
            blocks_data(i)%data_fields(j)%data_old(:,:)           = 9.0e9_rk

            blocks_data(i)%data_fields(j)%k1(:,:)                 = 9.0e9_rk
            blocks_data(i)%data_fields(j)%k2(:,:)                 = 9.0e9_rk
            blocks_data(i)%data_fields(j)%k3(:,:)                 = 9.0e9_rk
            blocks_data(i)%data_fields(j)%k4(:,:)                 = 9.0e9_rk

        end do


        blocks_data(i)%dx = params%Lx / real(blocks_params%size_domain,8)
        blocks_data(i)%dy = params%Lx / real(blocks_params%size_domain,8)

        ! coordinates
        allocate( blocks_data(i)%coord_x(Bs) )
        allocate( blocks_data(i)%coord_y(Bs) )

        blocks_data(i)%coord_x(:)    = 0.0_rk
        blocks_data(i)%coord_y(:)    = 0.0_rk

        ! light data id
        blocks_data(i)%block_id      = -1

    end do


end subroutine allocate_block_memory
