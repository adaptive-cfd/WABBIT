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
    ! number of processes
    integer(kind=ik)                    :: number_procs
    ! cpu time variables for running time calculation
    real(kind=rk)                       :: t0, t1, t2
    ! user defined parameter structure
    type (type_params)                  :: params

    integer(kind=ik), allocatable       :: lgt_block(:, :)
    real(kind=rk), allocatable          :: hvy_mask(:, :, :, :, :)
    integer(kind=ik), allocatable       :: hvy_neighbor(:,:)
    integer(kind=tsize), allocatable    :: lgt_sortednumlist(:, :, :)
    integer(kind=ik), allocatable       :: lgt_active(:,:)
    integer(kind=ik), allocatable       :: hvy_active(:,:)
    integer(kind=ik), allocatable       :: lgt_n(:)
    integer(kind=ik), allocatable       :: hvy_n(:)
    ! time loop variables
    real(kind=rk)                       :: time
    ! filename of *.ini file used to read parameters
    character(len=80)                   :: filename,fname
    integer(kind=ik) :: k, lgt_id, Bs(1:3), g, tree_n, hvy_id, iter, Jmax, Jmin, Jmin_equi
    real(kind=rk) :: x0(1:3), dx(1:3)
    logical :: pruned, help1, help2

    !---------------------------------------------------------------------------
    ! If called with '--help' or '-h', print a help message and exit.
    !---------------------------------------------------------------------------
    call get_cmd_arg( "--help", help1, default=.false. )
    call get_cmd_arg( "-h", help2, default=.false. )

    if ((help1 .or. help2) .and. (params%rank==0)) then
        write(*,*) "--------------------------------------------------------------"
        write(*,*) " Dry run"
        write(*,*) "--------------------------------------------------------------"
        write(*,*) " Just like a normal simulation, this postprocessing routine reads"
        write(*,*) " the parameters from an INI file, but it just generates the mask"
        write(*,*) " function at time intervals specified in the INI file. No fluid "
        write(*,*) " or other eqns are solved."
        write(*,*) "--------------------------------------------------------------"
        write(*,*) " This routine is mostly used for ACM and insects"
        write(*,*) "--------------------------------------------------------------"
        write(*,*) " Call:"
        write(*,*) " ./wabbit-post --dry-run PARAMS.ini --memory=20.0GB"
        write(*,*) ""
        write(*,*) " Other parameters:"
        write(*,*) ""
        write(*,*) " --pruned=1 (or simply --pruned) Uses tree-pruning, i.e. removes"
        write(*,*) "    blocks which do not contain the mask function. This option can"
        write(*,*) "    speed up visualization, but the data is incomplete: you cannot"
        write(*,*) "    read those fields into wabbit."
        write(*,*) ""
        write(*,*) " --Jmin"
        write(*,*) ""
        write(*,*) " [Insects]::smoothing_thickness=local;"
        write(*,*) ""
    endif



    !---------------------------------------------------------------------------
    ! Initialize parameters,bridge and grid
    !---------------------------------------------------------------------------
    ! read in the parameter file to setup the case
    ! get the second command line argument: this should be the ini-file name
    call get_command_argument( 2, filename )
    ! read ini-file and save parameters in struct
    call ini_file_to_params( params, filename )





    call get_cmd_arg( "--pruned", pruned, default=.false. )
    call get_cmd_arg( "--Jmin", Jmin_equi, default=params%min_treelevel )




    ! modifications to parameters
    deallocate( params%butcher_tableau )
    allocate( params%butcher_tableau(1,1) )
    ! mask, usx,usy,usz, color, sponge = 6 components
    params%n_eqn = 6

    Bs = params%Bs
    g  = params%n_ghosts
    Jmax = params%max_treelevel
    Jmin = params%min_treelevel
    tree_n = params%forest_size ! used only for resetting at this point
    time = 0.0_rk

    ! initializes the communicator for Wabbit and creates a bridge if needed
    call initialize_communicator(params)
    ! have the pysics module read their own parameters
    call init_physics_modules( params, filename, params%N_mask_components )

    ! allocate memory for heavy, light, work and neighbor data
    call allocate_forest(params, lgt_block, hvy_mask, hvy_neighbor, lgt_active, &
    hvy_active, lgt_sortednumlist, hvy_n=hvy_n, lgt_n=lgt_n)

    ! The ghost nodes will call their own setup on the first call, but for cleaner output
    ! we can also just do it now.
    call init_ghost_nodes( params )


    !-----------------------------------
    call init_t_file('kinematics.t', .true., (/ &
    "           time", &
    "    xc_body_g_x", &
    "    xc_body_g_y", &
    "    xc_body_g_z", &
    "            psi", &
    "           beta", &
    "          gamma", &
    "     eta_stroke", &
    "        alpha_l", &
    "          phi_l", &
    "        theta_l", &
    "        alpha_r", &
    "          phi_r", &
    "        theta_r", &
    "  rot_rel_l_w_x", &
    "  rot_rel_l_w_y", &
    "  rot_rel_l_w_z", &
    "  rot_rel_r_w_x", &
    "  rot_rel_r_w_y", &
    "  rot_rel_r_w_z", &
    "   rot_dt_l_w_x", &
    "   rot_dt_l_w_y", &
    "   rot_dt_l_w_z", &
    "   rot_dt_r_w_x", &
    "   rot_dt_r_w_y", &
    "   rot_dt_r_w_z"/) )
    !-----------------------------------

    do while ( time < params%time_max )
        ! start with an equidistant grid on coarsest level.
        ! routine also deletes any existing mesh in the tree.
        call create_equidistant_grid( params, lgt_block, hvy_neighbor, &
        lgt_active(:,tree_ID_mask), lgt_n(tree_ID_mask), lgt_sortednumlist(:,:,tree_ID_mask),&
        hvy_active(:,tree_ID_mask), hvy_n(tree_ID_mask), Jmin_equi, verbosity=.true., tree_ID=tree_ID_mask )


        if (params%rank==0) then
            write(*,'("Starting mask generation. Now: Jmax=",i2, " Nb=",i7)') &
            max_active_level(lgt_block, lgt_active(:,tree_ID_mask), lgt_n(tree_ID_mask)), lgt_n(tree_ID_mask)
        endif


        ! generate complete mask on the equidistan grid
        do k = 1, hvy_n(tree_ID_mask)
            hvy_id = hvy_active(k,tree_ID_mask)

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
            call sync_ghosts( params, lgt_block, hvy_mask, hvy_neighbor, hvy_active(:,1), hvy_n(1) )


            ! refine the mesh. Note: afterwards, it can happen that two blocks on the same level differ
            ! in their redundant nodes, but the ghost node sync'ing later on will correct these mistakes.
            call refine_mesh( params, lgt_block, hvy_mask, hvy_neighbor, &
            lgt_active(:,tree_ID_mask), lgt_n(tree_ID_mask), &
            lgt_sortednumlist(:,:,tree_ID_mask), hvy_active(:,tree_ID_mask), &
            hvy_n(tree_ID_mask), "mask-threshold", tree_ID_mask )


            if (params%rank==0) then
                write(*,'("Did one iteration for mask generation. Now: Jmax=",i2, " Nb=",i7)') &
                max_active_level(lgt_block, lgt_active(:,tree_ID_mask), lgt_n(tree_ID_mask)), lgt_n(tree_ID_mask)
            endif


            ! call create_mask_tree(params, time, lgt_block, hvy_mask, hvy_tmp, &
            ! hvy_neighbor, hvy_active, hvy_n, lgt_active, lgt_n, lgt_sortednumlist )
            do k = 1, hvy_n(tree_ID_mask)
                hvy_id = hvy_active(k,tree_ID_mask)

                call hvy_id_to_lgt_id( lgt_id, hvy_id, params%rank, params%number_blocks )
                call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

                call CREATE_MASK_meta( params%physics_type, time, x0, dx, Bs, g, &
                hvy_mask(:,:,:,:,hvy_id), "all-parts" )
            enddo
        enddo


        call WRITE_INSECT_DATA(time)

        if (pruned) then
            if (params%rank==0) write(*,*) "now pruning!"

            call prune_tree( params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
            hvy_mask, hvy_active, hvy_n, hvy_neighbor, tree_id=tree_ID_mask)
        endif

        !***********************************************************************
        ! Write fields to HDF5 file
        !***********************************************************************
        ! call save_data( iteration, time, params, lgt_block, hvy_block, lgt_active, &
        ! lgt_n, lgt_sortednumlist, hvy_n, hvy_tmp, hvy_active, hvy_mask, hvy_neighbor )

        ! create filename
        write( fname,'(a, "_", i12.12, ".h5")') "mask", nint(time * 1.0e6_rk)
        ! call write_field( fname, time, -99, 1, params, lgt_block, hvy_mask, lgt_active, lgt_n, hvy_n, hvy_active)

        call write_tree_field(fname, params, lgt_block, lgt_active, hvy_mask, &
        lgt_n, hvy_n, hvy_active, dF=1, tree_id=tree_ID_mask, time=time, iteration=-1 )

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
