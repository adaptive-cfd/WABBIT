subroutine post_dry_run

    use mpi
    use module_helpers
    use module_MPI
    ! global parameters
    use module_params
    ! timing module
    use module_timing
    ! init data module
    use module_initialization
    ! mesh manipulation subroutines
    use module_mesh
    ! IO module
    use module_IO
    ! time step module
    use module_time_step
    ! unit test module
    use module_unit_test
    ! bridge implementation of wabbit
    use module_bridge_interface
    use module_forest
    use module_mask
    ! HACK.We should load only the metamodule, but we require WRITE_INSECT_DATA(time)
    ! to dump kinematics data.
    use module_ACM

    implicit none

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank
    ! number of processes
    integer(kind=ik)                    :: number_procs
    ! cpu time variables for running time calculation
    real(kind=rk)                       :: t0, t1, t2
    ! user defined parameter structure
    type (type_params)                  :: params

    integer(kind=ik), allocatable       :: lgt_block(:, :)
    real(kind=rk), allocatable          :: hvy_mask(:, :, :, :, :)
    integer(kind=ik), allocatable       :: hvy_neighbor(:,:)
    integer(kind=tsize), allocatable    :: lgt_sortednumlist(:, :)
    integer(kind=ik), allocatable       :: lgt_active(:)
    integer(kind=ik), allocatable       :: hvy_active(:)
    integer(kind=ik)       :: lgt_n
    integer(kind=ik)       :: hvy_n
    ! time loop variables
    real(kind=rk)                       :: time
    ! filename of *.ini file used to read parameters
    character(len=80)                   :: filename,fname
    integer(kind=ik) :: k, lgt_id, Bs(1:3), g, tree_n, hvy_id, iter, Jmax, Jmin
    real(kind=rk) :: x0(1:3), dx(1:3)




    !---------------------------------------------------------------------------
    ! Initialize parameters,bridge and grid
    !---------------------------------------------------------------------------
    ! read in the parameter file to setup the case
    ! get the second command line argument: this should be the ini-file name
    call get_command_argument( 2, filename )
    ! read ini-file and save parameters in struct
    call ini_file_to_params( params, filename )

    ! modifications to parameters
    deallocate( params%butcher_tableau )
    allocate( params%butcher_tableau(1,1) )
    params%n_eqn = 6

    Bs = params%Bs
    g  = params%n_ghosts
    Jmax = params%max_treelevel
    Jmin = params%min_treelevel
    tree_n = params%forest_size ! used only for resetting at this point
    time = 0.0_rk

    ! initializes the communicator for Wabbit and creates a bridge if needed
    call initialize_communicator(params)
    ! have the pysics module read their own parameters. They also decide how many grid-qtys
    ! they want
    call init_physics_modules( params, filename, params%N_mask_components )

    ! allocate memory for heavy, light, work and neighbor data
    call allocate_grid(params, lgt_block, hvy_mask, hvy_neighbor, lgt_active, &
    hvy_active, lgt_sortednumlist)

    ! The ghost nodes will call their own setup on the first call, but for cleaner output
    ! we can also just do it now.
    call init_ghost_nodes( params )

    do while ( time < params%time_max )

        ! start with an equidistant grid on coarsest level.
        ! routine also deletes any existing mesh in the tree.
        call create_equidistant_grid( params, lgt_block, hvy_neighbor, &
        lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n, Jmin, verbosity=.true., tree_ID=tree_ID_mask )

        if (params%rank==0) then
            write(*,'("Starting mask generation. Now: Jmax=",i2, " Nb=",i7)') &
            max_active_level(lgt_block, lgt_active, lgt_n), lgt_n
        endif

        ! generate complete mask on the equidistan grid
        do k = 1, hvy_n

            hvy_id = hvy_active(k)

            call hvy_id_to_lgt_id( lgt_id, hvy_id, params%rank, params%number_blocks )
            call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

            call CREATE_MASK_meta( params%physics_type, time, x0, dx, Bs, g, &
            hvy_mask(:,:,:,:,hvy_id), "all-parts" )
        enddo


        ! refine the grid near the interface and re-generate the mask function.
        do iter = 1, Jmax - Jmin
            ! synchronization before refinement (because the interpolation takes place on the extended blocks
            ! including the ghost nodes)
            ! Note: at this point the grid is rather coarse (fewer blocks), and the sync step is rather cheap.
            ! Snyc'ing becomes much more expensive once the grid is refined.
            ! sync possible only before pruning
            call sync_ghosts( params, lgt_block, hvy_mask, hvy_neighbor, hvy_active, hvy_n )


            ! refine the mesh. Note: afterwards, it can happen that two blocks on the same level differ
            ! in their redundant nodes, but the ghost node sync'ing later on will correct these mistakes.
            call refine_mesh( params, lgt_block, hvy_mask, hvy_neighbor, &
            lgt_active, lgt_n, &
            lgt_sortednumlist, hvy_active, &
            hvy_n, "mask-threshold", tree_ID_mask )


            if (params%rank==0) then
                write(*,'("Did one iteration for mask generation. Now: Jmax=",i2, " Nb=",i7)') &
                max_active_level(lgt_block, lgt_active, lgt_n), lgt_n
            endif


            ! call create_mask_tree(params, time, lgt_block, hvy_mask, hvy_tmp, &
            ! hvy_neighbor, hvy_active, hvy_n, lgt_active, lgt_n, lgt_sortednumlist )
            do k = 1, hvy_n
                hvy_id = hvy_active(k)

                call hvy_id_to_lgt_id( lgt_id, hvy_id, params%rank, params%number_blocks )
                call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

                call CREATE_MASK_meta( params%physics_type, time, x0, dx, Bs, g, &
                hvy_mask(:,:,:,:,hvy_id), "all-parts" )
            enddo

        enddo


        call WRITE_INSECT_DATA(time)


        !***********************************************************************
        ! Write fields to HDF5 file
        !***********************************************************************
        ! call save_data( iteration, time, params, lgt_block, hvy_block, lgt_active, &
        ! lgt_n, lgt_sortednumlist, hvy_n, hvy_tmp, hvy_active, hvy_mask, hvy_neighbor )

        ! create filename
        write( fname,'(a, "_", i12.12, ".h5")') "mask", nint(time * 1.0e6_rk)
        call write_field( fname, time, -99, 1, params, lgt_block, hvy_mask, lgt_active, lgt_n, hvy_n, hvy_active)

        !write( fname,'(a, "_", i12.12, ".h5")') "usx", nint(time * 1.0e6_rk)
        !call write_field( fname, time, -99, 2, params, lgt_block, hvy_mask, lgt_active, lgt_n, hvy_n, hvy_active)
	!
        !write( fname,'(a, "_", i12.12, ".h5")') "usy", nint(time * 1.0e6_rk)
        !call write_field( fname, time, -99, 3, params, lgt_block, hvy_mask, lgt_active, lgt_n, hvy_n, hvy_active)
        !
        !write( fname,'(a, "_", i12.12, ".h5")') "usz", nint(time * 1.0e6_rk)
        ! call write_field( fname, time, -99, 4, params, lgt_block, hvy_mask, lgt_active, lgt_n, hvy_n, hvy_active)

        time = time + params%write_time
    end do



end subroutine post_dry_run
