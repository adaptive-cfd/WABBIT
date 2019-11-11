subroutine post_prune_tree(params)
    use module_precision
    use module_mesh
    use module_params
    use module_IO
    use module_mpi
    use module_operators
    use module_physics_metamodule
    use module_time_step
    use module_stl_file_reader
    use module_helpers
    use module_forest

    implicit none

    type (type_params), intent(inout)  :: params
    character(len=80) :: fname_ini, fname1, fname_out

    integer(kind=ik), allocatable      :: lgt_block(:, :)
    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik), allocatable      :: hvy_neighbor(:,:)
    integer(kind=ik), allocatable      :: lgt_active(:,:), hvy_active(:,:)
    integer(kind=tsize), allocatable   :: lgt_sortednumlist(:,:,:)
    integer(kind=ik), allocatable      :: lgt_n(:), hvy_n(:)
    integer :: hvy_id, lgt_id, fsize, j

    integer(kind=ik) :: iteration, Bs(1:3), tc_length1, dim, tc_length2, N1, N2, tree_N
    real(kind=rk) :: time, domain(1:3)

    !-----------------------------------------------------------------------------------------------------
    ! get values from command line (filename and level for interpolation)
    call get_command_argument(2, fname_ini)

    ! does the user need help?
    if (fname_ini=='--help' .or. fname_ini=='--h' .or. fname_ini=='-h') then
        if (params%rank==0) then
            write(*,*) "------------------------------------------------------------------"
            write(*,*) "./wabbit-post --prune-tree input_002.h5 output.h5"
            write(*,*) "------------------------------------------------------------------"
            write(*,*) " Removes all zero blocks from a tree (useful only for the mask function)"
            write(*,*) "------------------------------------------------------------------"
        end if
        return
    endif

    call get_command_argument(2, fname1)
    call check_file_exists(fname1)

    call get_command_argument(3, fname_out)

    N_MAX_COMPONENTS  = 1

    call read_attributes(fname1, N1, time, iteration, domain, params%Bs, tc_length1, params%dim)

    ! just to get some memory:
    params%number_blocks = 20 + N1 / params%number_procs
    params%domain_size = domain
    params%max_treelevel = tc_length1
    params%min_treelevel = 1
    params%n_eqn = 1
    params%n_ghosts = 2
    params%forest_size = 4
    fsize = params%forest_size
    params%order_predictor = "multiresolution_2nd"
    params%block_distribution = "sfc_hilbert"
    ! The ghost nodes will call their own setup on the first call, but for cleaner output
    ! we can also just do it now.
    call init_ghost_nodes( params )

    ! we have to allocate grid if this routine is called for the first time
    call allocate_forest(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, &
    hvy_active, lgt_sortednumlist, hvy_tmp=hvy_tmp, hvy_n=hvy_n, lgt_n=lgt_n)

    hvy_neighbor = -1
    lgt_n = 0 ! reset number of active light blocks
    tree_n = 0 ! reset number of trees in forest

    call read_field2tree(params, (/fname1/), 1, 1, tree_n, &
    lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
    hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor)

    call create_active_and_sorted_lists( params, lgt_block, lgt_active, &
    lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_n)

    call prune_tree( params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
    hvy_block, hvy_active, hvy_n, hvy_neighbor, tree_id=1)

    call create_active_and_sorted_lists( params, lgt_block, lgt_active, &
    lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_n)

    call write_tree_field(fname_out, params, lgt_block, lgt_active, hvy_block, &
    lgt_n, hvy_n, hvy_active, dF=1, tree_id=1, time=time, iteration=iteration )

    ! make a summary of the program parts, which have been profiled using toc(...)
    ! and print it to stdout
    call summarize_profiling( WABBIT_COMM )
end subroutine
