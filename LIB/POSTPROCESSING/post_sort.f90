subroutine post_sort(params)
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_helpers
    use module_forestMetaData

    implicit none

    type (type_params), intent(inout)  :: params
    character(len=clong) :: mode, fname, fname_out, args

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik) :: hvy_id, lgt_id, fsize, j, tree_ID, k_block

    integer(kind=ik) :: iteration, Bs(1:3), N, proc
    real(kind=rk) :: time, domain(1:3)
    integer(kind=ik), allocatable, dimension(:,:) :: tc_map_old
    integer(kind=tsize) :: treecode



    !-----------------------------------------------------------------------------------------------------
    ! get values from command line (filename and level for interpolation)
    call get_command_argument(1, mode)

    ! does the user need help?
    if (mode=='--help' .or. mode=='--h' .or. mode=='-h') then
        if (params%rank==0) then
            write(*,*) "------------------------------------------------------------------"
            write(*,*) "Sort array after treecode. Needed as it is assumed that data for input has same hvy block structure. Doesn't work for full-tree grids."
            write(*,*) "--------------------------------------------------------------"
            write(*,*) ""
            write(*,*) " Call:"
            write(*,*) " ./wabbit-post --sort DATA.h5 output.h5"
            write(*,*) ""
            write(*,*) "--------------------------------------------------------------"
        end if
        return
    endif

    call get_command_argument(2, fname)
    call check_file_exists(fname)

    call get_command_argument(3, fname_out)


    call read_attributes(fname, N, time, iteration, params%domain_size, params%Bs, params%Jmax, params%dim, periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)

    params%number_blocks = ceiling(1.1_rk * dble(N) / dble(params%number_procs)) ! just to get some memory in case not provided, extra space for load balancing

    params%Jmin = 0
    params%n_eqn = 1
    params%g = 1
    tree_ID = tree_ID_flow

    allocate(params%threshold_state_vector_component(params%n_eqn))
    params%threshold_state_vector_component(1:params%n_eqn)=1


    params%wavelet ='CDF20'
    ! We actually don't need the wavelet, but we need to call the setup to get the filter coefficients in case anything is hardcoded
    call setup_wavelet(params)

    params%forest_size = 20
    fsize = params%forest_size
    params%block_distribution = "sfc_z"

    ! we have to allocate grid if this routine is called for the first time
    call allocate_forest(params, hvy_block)

    ! The ghost nodes will call their own setup on the first call, but for cleaner output
    ! we can also just do it now.
    call init_ghost_nodes( params )

    ! read in data
    call readHDF5vct_tree((/fname/), params, hvy_block, tree_ID=1, verbosity=.true.)

    ! do load balancing after z curve to get all necessary data on the right rank
    call balanceLoad_tree( params, hvy_block, tree_ID )
    if (params%rank == 0) write(*, '(A)') "Load balancing finished."

    ! now we need to let each rank sort itself, we build array and then sort it
    allocate(tc_map_old(1:3, 1:hvy_n(tree_ID)))
    do k_block = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k_block, tree_ID)
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

        tc_map_old(1, k_block) = k_block  ! store hvy id
        tc_map_old(2:3, k_block) = lgt_block(lgt_id, IDX_TC_1:IDX_TC_2)  ! store treecode)
    enddo
    call quicksort(tc_map_old, 1, hvy_n(tree_ID), 3)
    ! rewrite indices of tc_map_old so that they are now in sorted order
    do k_block = 1, hvy_n(tree_ID)
        tc_map_old(1, k_block) = k_block  ! store hvy id
    enddo

    ! now let's reorder the hvy arrays according to the new mapping
    call reorder_hvy_arrays(params, hvy_block, tree_ID, tc_map_old, n_leaves=hvy_n(tree_ID))
    if (params%rank == 0) write(*, '(A)') "Reordering finished."

    ! let's save
    call saveHDF5_tree(fname_out, time, iteration, 1, params, hvy_block, tree_ID=tree_ID)

    call deallocate_forest(params, hvy_block, hvy_tmp=hvy_tmp )

    ! make a summary of the program parts, which have been profiled using toc(...)
    ! and print it to stdout
    call summarize_profiling( WABBIT_COMM )

end subroutine