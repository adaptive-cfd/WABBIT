subroutine post_prune_tree(params)
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_operators
    use module_physics_metamodule
    use module_time_step
    use module_stl_file_reader
    use module_helpers
    use module_forestMetaData

    implicit none

    type (type_params), intent(inout)  :: params
    character(len=cshort) :: fname_ini, fname1, fname_out

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_tmp(:, :, :, :, :)
    integer :: hvy_id, lgt_id, fsize, j

    integer(kind=ik) :: iteration, Bs(1:3), tc_length1, dim, tc_length2, N1, N2
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

    call read_attributes(fname1, N1, time, iteration, domain, params%Bs, tc_length1, params%dim, periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)

    ! just to get some memory:
    params%number_blocks = 20 + N1 / params%number_procs
    params%domain_size = domain
    params%Jmax = tc_length1
    params%Jmin = 1
    params%n_eqn = 1
    params%g = 2
    params%forest_size = 4
    fsize = params%forest_size
    params%order_predictor = "multiresolution_2nd"
    params%block_distribution = "sfc_hilbert"
    params%wavelet = "CDF20" ! does not matter here...

    call setup_wavelet(params)

    ! The ghost nodes will call their own setup on the first call, but for cleaner output
    ! we can also just do it now.
    call init_ghost_nodes( params )

    ! we have to allocate grid if this routine is called for the first time
    call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp)

    hvy_neighbor = -1
    lgt_n = 0 ! reset number of active light blocks
    tree_n = 0 ! reset number of trees in forest

    call readHDF5vct_tree((/fname1/), params, hvy_block, tree_ID=1, verbosity=.true.)

    call createActiveSortedLists_forest(params)

    call prune_tree( params, hvy_block, tree_ID=1)

    call createActiveSortedLists_forest(params)

    call saveHDF5_tree(fname_out, time, iteration, 1, params, hvy_block, 1)

    ! make a summary of the program parts, which have been profiled using toc(...)
    ! and print it to stdout
    call summarize_profiling( WABBIT_COMM )
end subroutine
