subroutine post_add_two_masks(params)
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
    character(len=cshort) :: mode, fname1, fname2, fname_out

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_work(:, :, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik) :: hvy_id, lgt_id, fsize, j, tree_ID, k

    integer(kind=ik) :: iteration, Bs(1:3), tc_length1, dim, tc_length2, N1, N2
    real(kind=rk) :: time, domain(1:3), norm

    integer(kind=ik) , save, allocatable :: lgt_active_ref(:,:), lgt_block_ref(:,:)
    integer(kind=ik) , save :: lgt_n_ref(2)=0_ik



    !-----------------------------------------------------------------------------------------------------
    ! get values from command line (filename and level for interpolation)
    call get_command_argument(1, mode)

    ! does the user need help?
    if (mode=='--help' .or. mode=='--h' .or. mode=='-h') then
        if (params%rank==0) then
            write(*,*) "------------------------------------------------------------------"
            write(*,*) "./wabbit-post --add mask_001.h5 mask_002.h5 output.h5"
            write(*,*) "./wabbit-post --subtract mask_001.h5 mask_002.h5 output.h5"
            write(*,*) "./wabbit-post --multiply mask_001.h5 mask_002.h5 output.h5"
            write(*,*) "./wabbit-post --grid1-to-grid2 mask_001.h5 mask_002.h5 output.h5"
            write(*,*) "./wabbit-post --noise-like-grid1 mask_001.h5 mask_002.h5 output.h5"
            write(*,*) "------------------------------------------------------------------"
            write(*,*) ""
            write(*,*) ""
            write(*,*) ""
            write(*,*) ""
            write(*,*) "------------------------------------------------------------------"
        end if
        return
    endif

    call get_command_argument(2, fname1)
    call check_file_exists(fname1)

    call get_command_argument(3, fname2)
    call check_file_exists(fname2)

    call get_command_argument(4, fname_out)

    call read_attributes(fname1, N1, time, iteration, domain, params%Bs, tc_length1, params%dim, periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)
    call read_attributes(fname2, N2, time, iteration, domain, params%Bs, tc_length2, params%dim, periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)


    if (mode=="--test_operations") then
        params%number_blocks = 5*max(N1,N2) ! just to get some memory:
    else
        params%number_blocks = 3*max(N1,N2) ! just to get some memory:
    end if
    params%domain_size = domain
    params%Jmax = max( tc_length2, tc_length1 )
    params%Jmin = 1
    params%n_eqn = 1

    allocate(params%threshold_state_vector_component(params%n_eqn))
    params%threshold_state_vector_component(1:params%n_eqn)=1


    call get_cmd_arg_str( "--wavelet", params%wavelet, default='CDF40' )
    ! initialize wavelet transform
    ! also, set number of ghost nodes params%G to minimal value for this wavelet
    call setup_wavelet(params, params%g)

    params%forest_size = 20
    fsize = params%forest_size
    params%block_distribution = "sfc_hilbert"
    params%time_step_method = 'none'


    ! we have to allocate grid if this routine is called for the first time
    call allocate_forest(params, hvy_block, hvy_work, hvy_tmp=hvy_tmp)

    ! The ghost nodes will call their own setup on the first call, but for cleaner output
    ! we can also just do it now.
    call init_ghost_nodes( params )


    call reset_forest(params)
    lgt_n = 0 ! reset number of active light blocks
    hvy_n = 0
    tree_n = 0 ! reset number of trees in forest

    call readHDF5vct_tree((/fname1/), params, hvy_block, tree_ID=1, verbosity=.false.)
    call readHDF5vct_tree((/fname2/), params, hvy_block, tree_ID=2, verbosity=.false.)

    call createActiveSortedLists_forest(params)


    select case(mode)
    case ("--noise-like-grid1")
        do k = 1, hvy_n(1)
            hvy_id = hvy_active(k, 1)
            call random_data(hvy_block(:,:,:,:,hvy_id))
        enddo

    case ("--add")
        call add_two_trees(params, hvy_block, hvy_tmp, tree_ID1=1, tree_ID2=2, verbosity=.true.)

    case ("--subtract")
        call substract_two_trees(params, hvy_block, hvy_tmp, tree_ID1=1, tree_ID2=2)

    case ("--multiply")
        call multiply_two_trees(params, hvy_block, hvy_tmp, tree_ID1=1, tree_ID2=2)
    
    case ("--grid1-to-grid2")
        ! Store (both) grids (=only treedata) in the *_ref variables
        call store_ref_meshes( lgt_block_ref, lgt_active_ref, lgt_n_ref, tree_ID1=1, tree_ID2=2)
        ! before applying wavelet filters we need to sync
        call sync_ghosts_tree(params, hvy_block, tree_ID=1)
        ! make grid1 equal to grid 1 (by only coarsening).
	! Should refinement be required error will be thrown.
        call coarse_tree_2_reference_mesh(params, lgt_block_ref, lgt_active_ref(:,2), lgt_n_ref(2), &
        hvy_block, hvy_tmp, tree_ID=1, verbosity=.true.)

    case ("--test_operations")
        ! this tests inplace vs out of place operations
        ! z = x*y + x
        ! norm(w)
         call multiply_two_trees(params, hvy_block, hvy_tmp, tree_ID1=1, tree_ID2=2, dest_tree_ID=3, verbosity=.True.)
         call add_two_trees(params, hvy_block, hvy_tmp, tree_ID1=1, tree_ID2=3, dest_tree_ID=4, verbosity=.True.)
        ! y <- x*y + x
         call multiply_two_trees(params, hvy_block, hvy_tmp, tree_ID1=2, tree_ID2=1, verbosity=.True.)
         call add_two_trees(params, hvy_block, hvy_tmp, tree_ID1=2, tree_ID2=1, verbosity=.True.)
         ! w = z-y
         call substract_two_trees(params, hvy_block, hvy_tmp, tree_ID1=2, tree_ID2=4, dest_tree_ID=3)
        !
        norm =  scalar_product_two_trees( params, hvy_block, hvy_tmp, tree_ID1=2, tree_ID2=2)

        if (params%rank==0) then
            write(*,*) "Norm (should be not 0): ", norm
        end if
        norm =  scalar_product_two_trees( params, hvy_block, hvy_tmp ,&
                              tree_ID1=4, tree_ID2=4)
        if (params%rank==0) then
          write(*,*) "Norm (should be not 0): ", norm
        end if


        norm =  scalar_product_two_trees( params, hvy_block, hvy_tmp ,&
                                   tree_ID1=3, tree_ID2=3)
        if (params%rank==0) then
            write(*,*) "Norm (should be 0): ", norm
        end if

    case default
        call abort(301219,"unkown mode...")
    end select

    call createActiveSortedLists_forest(params)

    tree_ID = 1
    call updateNeighbors_tree( params, tree_ID, search_overlapping=.false.)

    tree_ID = 2
    call updateNeighbors_tree( params, tree_ID, search_overlapping=.false.)

    call saveHDF5_tree(fname_out, time, iteration, 1, params, hvy_block, tree_ID=1)

    call deallocate_forest(params, hvy_block, hvy_work, hvy_tmp )

    ! make a summary of the program parts, which have been profiled using toc(...)
    ! and print it to stdout
    call summarize_profiling( WABBIT_COMM )
end subroutine
