! ********************************
! 2D AMR prototype
! --------------------------------
!
! allocate memory for block data
! and initialize all data
!
! name: allocate_block_memory.f90
! date: 31.08.2016
! author: engels, msr
! version: 0.1
!
! ********************************

subroutine  allocate_block_memory()

    use module_params
    use module_blocks

    implicit none

    integer                                         :: i, j, allocate_error
    integer(kind=ik)                                :: Bs, g, dF

    ! grid size parameters
    Bs  = blocks_params%size_block
    g   = blocks_params%number_ghost_nodes
    dF  = blocks_params%number_data_fields

    write(*,'("System is allocating",i7," blocks" )') blocks_params%number_max_blocks
    write(*,'(80("-"))')

    allocate( blocks(1:blocks_params%number_max_blocks), stat=allocate_error )

    ! dummy allocation
    allocate( blocks_params%active_list(1), stat=allocate_error )

    do i = 1, blocks_params%number_max_blocks

        allocate( blocks(i)%data_fields(1:dF), stat=allocate_error )

        ! loop over all data fields
        do j = 1, dF

            allocate( blocks(i)%data_fields(j)%data_(Bs+2*g, Bs+2*g), stat=allocate_error )
            allocate( blocks(i)%data_fields(j)%data_old(Bs+2*g, Bs+2*g), stat=allocate_error )

            allocate( blocks(i)%data_fields(j)%k1(Bs+2*g, Bs+2*g), stat=allocate_error )
            allocate( blocks(i)%data_fields(j)%k2(Bs+2*g, Bs+2*g), stat=allocate_error )
            allocate( blocks(i)%data_fields(j)%k3(Bs+2*g, Bs+2*g), stat=allocate_error )
            allocate( blocks(i)%data_fields(j)%k4(Bs+2*g, Bs+2*g), stat=allocate_error )

        end do

        allocate( blocks(i)%treecode(10), stat=allocate_error )
        allocate( blocks(i)%neighbor_treecode(8,10), stat=allocate_error )
        allocate( blocks(i)%neighbor2_treecode(4,10), stat=allocate_error )

        ! loop over all data fields
        do j = 1, dF

            blocks(i)%data_fields(j)%data_(:,:)              = 0.0_rk
            blocks(i)%data_fields(j)%data_old(:,:)           = 0.0_rk

            blocks(i)%data_fields(j)%k1(:,:)                 = 0.0_rk
            blocks(i)%data_fields(j)%k2(:,:)                 = 0.0_rk
            blocks(i)%data_fields(j)%k3(:,:)                 = 0.0_rk
            blocks(i)%data_fields(j)%k4(:,:)                 = 0.0_rk

        end do

        blocks(i)%active                    = .false.
        blocks(i)%treecode(:)               = -1
        blocks(i)%neighbor_treecode(:,:)    = -1
        blocks(i)%neighbor2_treecode(:,:)   = -1
        blocks(i)%neighbor_id(:)            = -1
        blocks(i)%neighbor2_id(:)           = -1
        blocks(i)%refinement                = 0
        blocks(i)%level                     = -1
        blocks(i)%neighbor_number(:)        = 1
        blocks(i)%boundary(:)               = .false.
        blocks(i)%boundary2(:)              = .false.

        blocks(i)%dx = params%Lx / real(blocks_params%size_domain,8)
        blocks(i)%dy = params%Lx / real(blocks_params%size_domain,8)

        ! coordinates
        allocate( blocks(i)%coord_x(Bs) )
        allocate( blocks(i)%coord_y(Bs) )

        blocks(i)%coord_x(:)    = 0.0_rk
        blocks(i)%coord_y(:)    = 0.0_rk

    end do

end subroutine allocate_block_memory
