
subroutine post_dump_neighbors(params)
    use module_precision
    use module_mesh
    use module_params
    use module_IO
    use module_forest
    use module_mpi
    use module_operators

    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=80)      :: file, operator
    real(kind=rk)          :: time
    integer(kind=ik)       :: iteration, k, lgt_id, lgt_n, hvy_n, tc_length
    integer(kind=ik), dimension(3) :: Bs
    character(len=2)       :: order

    integer(kind=ik), allocatable      :: lgt_block(:, :)
    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_work(:, :, :, :, :, :), hvy_tmp(:, :, :, :, :)
    integer(kind=ik), allocatable      :: hvy_neighbor(:,:)
    integer(kind=ik), allocatable      :: lgt_active(:), hvy_active(:)
    integer(kind=tsize), allocatable   :: lgt_sortednumlist(:,:)
    real(kind=rk), dimension(3)        :: dx, x0
    integer(hid_t)                     :: file_id
    real(kind=rk), dimension(3)        :: domain

    call get_command_argument(1, operator)
    call get_command_argument(2, file)

    ! does the user need help?
    if (file=='--help' .or. file=='--h') then
        if (params%rank==0) then
            write(*,*) "-----------------------------------------------------------"
            write(*,*) " Wabbit postprocessing: dump neighbors"
            write(*,*) "-----------------------------------------------------------"
            write(*,*) " Read in a data field (2D or 3D) and computes the neighbor relations"
            write(*,*) " (update_neighbors) for all blocks. then, dumps lgt_block and hvy_neighbor to "
            write(*,*) " ascii file. Useful for development"
            write(*,*) " Best used in serial monoproc: then, hvy_id == lgt_id"
        end if
        return
    endif

    call check_file_exists(trim(file))

    if (params%number_procs /= 1) then
        call abort( 202101281, "use on monoproc")
    endif


    ! get some parameters from one of the files (they should be the same in all of them)
    call read_attributes(file, lgt_n, time, iteration, domain, Bs, tc_length, params%dim, periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)


    ! unused so just fill any value
    params%order_discretization = "FD_2nd_central"
    params%order_predictor = "multiresolution_2nd"
    params%n_ghosts = 2_ik
    params%max_treelevel = tc_length+1
    params%n_eqn = params%dim
    params%domain_size(1) = domain(1)
    params%domain_size(2) = domain(2)
    params%domain_size(3) = domain(3)
    params%Bs = Bs

    allocate(params%butcher_tableau(1,1))
    ! no refinement is made in this postprocessing tool; we therefore allocate about
    ! the number of blocks in the file (and not much more than that)
    params%number_blocks = ceiling(  real(lgt_n)/real(params%number_procs) )

    ! allocate data
    call allocate_grid(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, hvy_active, lgt_sortednumlist, hvy_tmp=hvy_tmp)

    ! read mesh and field
    call read_tree((/file/), 1, params, lgt_n, lgt_block, hvy_block, hvy_tmp, tree_id_optional=1)

    ! create lists of active blocks (light and heavy data)
    ! update list of sorted nunmerical treecodes, used for finding blocks
    call create_active_and_sorted_lists( params, lgt_block, lgt_active, lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_ID=1)

    ! update neighbor relations
    call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n )

    open(14, file="lgt_block.txt", status='replace')
    do k = 1, size(lgt_block,1)
        write(14,'(20(i2,1x))') lgt_block(k,:)
    enddo
    close(14)

    open(14, file="hvy_neighbors.txt", status='replace')
    do k = 1, size(hvy_neighbor,1)
        write(14,'(90(i5,1x))') hvy_neighbor(k,:)
    enddo
    close(14)


    open(14, file="x0_dx.txt", status='replace')
    do k = 1, size(lgt_block,1)
        call get_block_spacing_origin( params, k, lgt_block, x0, dx )
        write(14,'(6(es15.6,1x))') x0, dx
    enddo
    close(14)


end subroutine post_dump_neighbors
