! ********************************
! 2D AMR prototype
! --------------------------------
!
! create block tree from domain at start
!
! name: matrix_to_block_tree.f90
! date: 12.08.2016
! author: msr
! version: 0.1
!
! ********************************

subroutine matrix_to_block_tree()

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik)                                :: num_blocks_x, num_blocks_y, num_blocks
    integer(kind=ik)                                :: Bs, Ds, g, dF, i, j, ib, allocate_error
    integer(kind=ik), dimension(10)                 :: treecode

    real(kind=rk), dimension(:,:), allocatable      :: data
    real(kind=rk), dimension(:), allocatable        :: coord_x, coord_y, domain_coord_x, domain_coord_y

    Bs                  = blocks_params%size_block
    g                   = blocks_params%number_ghost_nodes
    Ds                  = blocks_params%size_domain
    dF                  = blocks_params%number_data_fields

    ! allocate memory for data fields and coordinates vectors
    allocate( data(Bs+2*g, Bs+2*g), stat=allocate_error )

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

    write(*,'(a,i5,a,i5,a,i5,a,i5,a,i5,a,i5)') "Field with res: ", (Ds), " x", (Ds), " gives: ", &
    num_blocks_x, " x", num_blocks_y, " (", num_blocks, ") blocks of size: ", (Bs)
    write(*,'(80("-"))')

    ! domain coordinate vectors
    do i = 1, Ds
        domain_coord_x(i) = (i-1) * params%Lx / (Ds-1)
        domain_coord_y(i) = params%Lx - (i-1) * params%Ly / (Ds-1)
    end do

    ! create new blocks
    if (dF == 1) then

        ! exactly one field
        ib = 1
        do i = 1, num_blocks_x
            do j = 1, num_blocks_y
                ! coordinates
                coord_x                     = domain_coord_x( (j-1)*(Bs-1) + 1 : j*(Bs-1) + 1 )
                coord_y                     = domain_coord_y( (i-1)*(Bs-1) + 1 : i*(Bs-1) + 1 )
                ! data
                data(:,:)                   = 0.0_rk
                data(g+1:Bs+g,g+1:Bs+g)     = blocks_params%phi( (i-1)*(Bs-1) + 1 : i*(Bs-1) + 1 , (j-1)*(Bs-1) + 1 : j*(Bs-1) + 1 )
                ! treecode
                call encoding(treecode, i, j, num_blocks_x, num_blocks_y)
                ! write new block
                call new_block(ib, treecode, 10, data, coord_x, coord_y, Bs, g, dF)
                ! block counter
                ib                          = ib + 1
            end do
        end do

    else
        ! more than one field
        ! to do
    end if

    ! deallocate memory for local variables
    deallocate( data, stat=allocate_error )
    deallocate( coord_x, stat=allocate_error )
    deallocate( coord_y, stat=allocate_error )
    deallocate( domain_coord_x, stat=allocate_error )
    deallocate( domain_coord_x, stat=allocate_error )

end subroutine matrix_to_block_tree
