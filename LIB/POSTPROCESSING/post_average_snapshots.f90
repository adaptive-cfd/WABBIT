
subroutine post_average_snapshots(params)
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
    character(len=80) :: fname_ini, fname2, fname_out

    character(len=80), allocatable     :: fname_in(:)
    integer(kind=ik), allocatable      :: lgt_block(:, :)
    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_work(:, :, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_tmp(:, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_gridQ(:, :, :, :, :)
    integer(kind=ik), allocatable      :: hvy_neighbor(:,:)
    integer(kind=ik), allocatable      :: lgt_active(:,:), hvy_active(:,:)
    integer(kind=tsize), allocatable   :: lgt_sortednumlist(:,:,:)
    integer(kind=ik), allocatable      :: lgt_n(:), hvy_n(:)
    integer(kind=ik), allocatable      :: Nblocks(:)
    integer :: hvy_id, lgt_id, fsize, i, tree_id, average_tree_id
    integer(kind=ik) :: iteration, Bs(1:3), tc_length1, dim, tc_length2, Nargs, &
                        tree_N, level, N_snapshots
    real(kind=rk) :: time, domain(1:3)
    logical,parameter :: verbosity = .false.
    !-----------------------------------------------------------------------------------------------------
    ! get values from command line (filename and level for interpolation)
    call get_command_argument(2, fname_ini)

    ! does the user need help?
    if (fname_ini=='--help' .or. fname_ini=='--h' .or. fname_ini=='-h') then
        if (params%rank==0) then
            write(*,*) "------------------------------------------------------------------"
            write(*,*) "./wabbit-post --average field_0001.h5 field_0002.h5 ... field_XXXX.h5 output.h5"
            write(*,*) "------------------------------------------------------------------"
            write(*,*) ""
            write(*,*) ""
            write(*,*) ""
            write(*,*) ""
            write(*,*) "------------------------------------------------------------------"
        end if
        return
    endif

    Nargs = command_argument_count()
    N_snapshots = Nargs - 2
    allocate(fname_in(N_snapshots))
    allocate(Nblocks(N_snapshots))

    ! get names and attributes of input files
    do i = 1, N_snapshots

      call get_command_argument(i+1, fname_in(i)) ! i + 1 because we skip the first argument
      call check_file_exists(fname_in(i))

      !! read attributes
      if ( i == 1 ) then
          call read_attributes(fname_in(i), Nblocks(1), time, iteration, &
               params%domain_size, params%Bs, params%max_treelevel, params%dim, periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)
      endif
      call read_attributes(fname_in(i), Nblocks(i), time, iteration, &
           domain, bs, level, dim, periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)

      ! check attributes for consistency
      params%max_treelevel = max(params%max_treelevel, level) ! find the maximal level of all snapshot
      if (any(params%Bs .ne. Bs)) call abort( 203191, " Block size is not consistent ")
      if ( abs(sum(params%domain_size(1:dim) - domain(1:dim))) > 1e-14 ) call abort( 203192, "Domain size is not consistent ")
      if (params%dim .ne. dim) call abort( 203193, "Dimension is not consistent ")
    end do

    ! read name of output file:
    call get_command_argument(Nargs, fname_out)


    params%number_blocks = (N_snapshots+1)*maxval(Nblocks) ! just to get some memory:
    params%min_treelevel = 1
    params%n_eqn = 1
    params%n_ghosts = 4
    params%forest_size = N_snapshots + 1
    fsize = params%forest_size
    params%order_predictor = "multiresolution_4th"
    params%block_distribution = "sfc_hilbert"
    params%time_step_method = 'none'


    ! we have to allocate grid if this routine is called for the first time
    call allocate_forest(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, &
    hvy_active, lgt_sortednumlist, hvy_work, hvy_tmp=hvy_tmp, hvy_n=hvy_n, lgt_n=lgt_n)

    ! The ghost nodes will call their own setup on the first call, but for cleaner output
    ! we can also just do it now.
    call init_ghost_nodes( params )

    hvy_neighbor = -1
    lgt_n = 0 ! reset number of active light blocks
    hvy_n = 0
    tree_n = 0 ! reset number of trees in forest

    do tree_id = 1, N_snapshots
      call read_field2tree(params,fname_in(tree_id) , params%n_eqn, tree_id, &
                  tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, hvy_block, &
                  hvy_active, hvy_n, hvy_tmp, hvy_neighbor)
      if (params%adapt_mesh) then
          ! now, evaluate the refinement criterion on each block, and coarsen the grid where possible.
          ! adapt-mesh also performs neighbor and active lists updates
      !    call adapt_tree_mesh( time(tree_id), params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, &
      !    lgt_sortednumlist, hvy_active, hvy_n, params%coarsening_indicator, hvy_tmp ,tree_id, tree_n)
      endif
    end do


    average_tree_id = N_snapshots + 1
    call average_trees( params, tree_N, lgt_block, lgt_active,&
                        lgt_n, lgt_sortednumlist, hvy_block, hvy_active,&
                        hvy_n, hvy_tmp, hvy_neighbor, (/(i, i= 1, N_snapshots)/), &
                        average_tree_id, verbosity)

    call write_tree_field(fname_out, params, lgt_block, lgt_active, hvy_block, &
    lgt_n, hvy_n, hvy_active, dF=1, tree_id=average_tree_id, time=0.0_rk, iteration=0 )


    call deallocate_forest(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, hvy_active, &
    lgt_sortednumlist, hvy_work, hvy_tmp, hvy_n, lgt_n )

    ! make a summary of the program parts, which have been profiled using toc(...)
    ! and print it to stdout
    call summarize_profiling( WABBIT_COMM )
end subroutine
