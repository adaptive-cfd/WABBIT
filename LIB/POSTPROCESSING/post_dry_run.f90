subroutine post_dry_run

    use mpi
    use module_helpers
    use module_MPI
    use module_params               ! global parameters
    use module_timing
    use module_initialization       ! init data module
    use module_mesh                 ! mesh manipulation subroutines
    use module_time_step
    use module_unit_test
    use module_bridge_interface     ! bridge implementation of wabbit
    ! HACK.We should load only the metamodule, but we require WRITE_INSECT_DATA(time)
    ! to dump kinematics data.
    use module_ACM
    use module_forestMetaData

    implicit none

    integer(kind=ik)                    :: ierr             ! MPI error variable
    integer(kind=ik)                    :: number_procs     ! number of processes
    real(kind=rk)                       :: t0, t1, t2       ! cpu time variables for running time calculation
    type (type_params)                  :: params           ! user defined parameter structure

    real(kind=rk), allocatable          :: hvy_mask(:, :, :, :, :), hvy_tmp(:, :, :, :, :)
    real(kind=rk)                       :: time             ! time loop variables
    character(len=cshort)               :: filename, fname, grid_list
    integer(kind=ik) :: k, lgt_id, Bs(1:3), g, hvy_id, iter, Jmax, Jmin, Jmin_equi, Jnow, Nmask, io_error, lgt_n_old, lgt_n_new, iteration
    real(kind=rk) :: x0(1:3), dx(1:3)
    logical :: pruned, help1, help2, save_us, iterate, error_OOM, save_color
    type(inifile) :: FILE

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas

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
        write(*,*) " Two modes are available: "
        write(*,*) "    1/ create mask on minimum grid, which may be different from fluid grid"
        write(*,*) "    2/ create mask on a list of given grids, then identical to fluid grid"
        write(*,*) "--------------------------------------------------------------"
        write(*,*) " Call:"
        write(*,*) " ./wabbit-post --dry-run PARAMS.ini --memory=20.0GB --Jmin=1 --grid-list=none --pruned"
        write(*,*) ""
        write(*,*) " Other parameters:"
        write(*,*) ""
        write(*,*) " --pruned=1 (or simply --pruned) Uses tree-pruning, i.e. removes"
        write(*,*) "    blocks which do not contain the mask function. This option can"
        write(*,*) "    speed up visualization, but the data is incomplete: you cannot"
        write(*,*) "    read those fields into wabbit."
        write(*,*) ""
        write(*,*) " --Jmin The minimum tree level, default is taken from INI file."
        write(*,*) ""
        write(*,*) " --save-us Save, in addition to mask_*.h5, also the solid velocity field us (a vector field)"
        write(*,*) ""
        write(*,*) " --grid-list=list.txt"
        write(*,*) "   If given, we read in the specified h5 files and create the mask"
        write(*,*) "   on the same grid. Output is stored in mask_XXXXXXXXXXXX.h5 as usual, "
        write(*,*) "   but the time stamp is taken from the file(s)."
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
    call get_cmd_arg( "--save-color", save_color, default=.false. )
    call get_cmd_arg( "--save-us", save_us, default=.false. )
    call get_cmd_arg( "--Jmin", Jmin_equi, default=params%Jmin )
    call get_cmd_arg( "--grid-list", grid_list, default="none" )
    

    ! for the dry run we dont need to use the fancy wavelets
    params%wavelet = "CDF20"
    call setup_wavelet(params)

    ! modifications to parameters (because we use hvy_block instead of hvy_mask, NEQN set
    ! in ini file is not correct)
    deallocate( params%butcher_tableau )
    allocate( params%butcher_tableau(1,1) )
    ! mask, usx,usy,usz, color, sponge = 6 components
    params%n_eqn = 6
    deallocate(params%threshold_state_vector_component)
    allocate(params%threshold_state_vector_component(1:params%n_eqn))
    params%threshold_state_vector_component = 0
    params%threshold_state_vector_component(1) = 1

    deallocate(params%symmetry_vector_component)
    allocate(params%symmetry_vector_component(1:params%n_eqn))
    params%symmetry_vector_component = "0"

    ! it is generally desired to create the mask on Jmax, which is the finest
    ! level used in the simulation. This is where the RHS is computed. If the dealiasing
    ! switch is .true., blocks on Jmx are then forced to be coarsened to Jmax-1
    params%force_maxlevel_dealiasing = .false.
    ! output this to user, because elsewise it might be confusing. This overwrites the PARAMS-file but is not output to the log-file elsewise
    if (params%force_maxlevel_dealiasing) then
        if (params%rank==0) write(*,'(A)') "Force maxlevel dealiasing set to TRUE"
    else
        if (params%rank==0) write(*,'(A)') "Force maxlevel dealiasing set to FALSE"
    endif

    params%eps = 1.0e-6
    params%coarsening_indicator = "threshold-state-vector"

    if (params%rank==0) then
        write(*,'(A)') "DRY-RUN: creating mask function. Please note that the params"
        write(*,'(A)') "         eps=1.0e-5, force_maxlevel_dealiasing=false and"
        write(*,'(A)') "         threshold_state_vector_component are hard-coded"
        write(*,'(A)') "         and thus NOT read from the INI file."
    endif

    Bs = params%Bs
    g  = params%g
    Jmax = params%Jmax
    Jmin = params%Jmin
    tree_n = params%forest_size ! used only for resetting at this point
    time = 0.0_rk

    ! initializes the communicator for Wabbit and creates a bridge if needed
    call initialize_communicator(params)
    ! have the pysics module read their own parameters
    call init_physics_modules( params, filename, params%N_mask_components )



    ! allocate memory for heavy, light, work and neighbor data
    call allocate_forest(params, hvy_mask, hvy_tmp=hvy_tmp, neqn_hvy_tmp=6)

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

    if (grid_list == "none" ) then
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+
        ! traditional mode: create the grid AND the mask, possibly with pruning.
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~+

        do while ( time < params%time_max )
            ! start with an equidistant grid on coarsest level.
            ! routine also deletes any existing mesh in the tree.
            call createEquidistantGrid_tree( params, hvy_mask, Jmin_equi, verbosity=.true., tree_ID=tree_ID_flow )


            if (params%rank==0) then
                write(*,'("Starting mask generation. Now: Jmax=",i2, " Nb=",i7," time=",g12.4)') &
                maxActiveLevel_tree(tree_ID_flow), lgt_n(tree_ID_flow), time
            endif

            ! generate complete mask on the initial equidistant grid
            call createMask_tree(params, time, hvy_mask, hvy_mask, .false.)


            ! refine the grid near the interface and re-generate the mask function.
            iterate = .true.
            iteration = 0
            do while (iterate)
                lgt_n_old = lgt_n(tree_ID_flow)

                ! synchronization before refinement (because the interpolation takes place on the extended blocks
                ! including the ghost nodes)
                ! Note: at this point the grid is rather coarse (fewer blocks), and the sync step is rather cheap.
                ! Snyc'ing becomes much more expensive once the grid is refined.
                ! Sync possible only before pruning.
                call sync_ghosts_tree( params, hvy_mask, tree_ID_flow )

                ! refine the mesh, but only where the mask is interesting (not everywhere!)
                call refine_tree( params, hvy_mask, hvy_tmp, "mask-threshold", tree_ID_flow, error_OOM )

                if (error_OOM) call abort(2512177,"Refinement failed, out of memory. Try with more memory.")

                ! on new grid, create the mask again
                call createMask_tree(params, time, hvy_mask, hvy_mask, .false.)
                Nmask = lgt_n(tree_ID_flow)

                ! note we do not pass hvy_mask in the last argument, so the switch params%threshold_mask
                ! is effectively ignored. It seems redundant; if we set a small eps (done independent
                ! of the parameter file), this yields the same result
                call adapt_tree( time, params, hvy_mask, tree_ID_flow, params%coarsening_indicator, hvy_tmp, hvy_mask=hvy_mask )
                lgt_n_new = lgt_n(tree_ID_flow)

                ! on new grid, create the mask again
                call createMask_tree(params, time, hvy_mask, hvy_mask, .false.)

                ! current finest level is:
                Jnow = maxActiveLevel_tree(tree_ID_flow)

                if (params%rank==0) then
                    write(*,'("Did iteration ", i0," for mask generation. Mask computed on ",i6," blocks.&
                    & After coarsening: Jmax=",i2, " Nb=",i7)') iteration, Nmask, Jnow, lgt_n(tree_ID_flow)
                endif

                ! We're done once the mask is created on the final level. Relevant only if the start grid is not
                ! created on Jmin, but on Jequi. A second condition is that the grid does not change anymore:
                ! sometimes, the start level is so coarse that the blocks retain not a single point of the mask
                ! function, which can lead to not all parts of the grid being properly generated. If *none* of the
                ! blocks contains a mask function point on Jmin, the only solution is to increase --Jmin in the call.
                if ((Jnow==Jmax) .and. (lgt_n_new==lgt_n_old)) iterate = .false.
                iteration = iteration + 1
            enddo

            ! on new grid, create the mask again
            call createMask_tree(params, time, hvy_mask, hvy_mask, .false.)
            Nmask = lgt_n(tree_ID_flow)

            call WRITE_INSECT_DATA(time)

            ! before (possible) pruning, we sync the ghosts
            call sync_ghosts_tree( params, hvy_mask, tree_ID_flow )

            if (pruned) then
                if (params%rank==0) write(*,*) "now pruning!"

                call prune_tree( params, hvy_mask, tree_ID=tree_ID_flow)
            endif

            !***********************************************************************
            ! Write fields to HDF5 file
            !***********************************************************************
            write( fname,'(a, "_", i12.12, ".h5")') "mask", nint(time * 1.0e6_rk)
            call saveHDF5_tree(fname, time, -1_ik, 1, params, hvy_mask, tree_ID_flow, no_sync=pruned)

            if (save_us) then
                write( fname,'(a, "_", i12.12, ".h5")') "usx", nint(time * 1.0e6_rk)
                call saveHDF5_tree(fname, time, -1_ik, 2, params, hvy_mask, tree_ID_flow, no_sync=pruned)

                write( fname,'(a, "_", i12.12, ".h5")') "usy", nint(time * 1.0e6_rk)
                call saveHDF5_tree(fname, time, -1_ik, 3, params, hvy_mask, tree_ID_flow, no_sync=pruned)

                write( fname,'(a, "_", i12.12, ".h5")') "usz", nint(time * 1.0e6_rk)
                call saveHDF5_tree(fname, time, -1_ik, 4, params, hvy_mask, tree_ID_flow, no_sync=pruned)
            endif

            if (save_color) then
                write( fname,'(a, "_", i12.12, ".h5")') "color", nint(time * 1.0e6_rk)
                call saveHDF5_tree(fname, time, -1_ik, 5, params, hvy_mask, tree_ID_flow, no_sync=pruned)
            endif

            time = time + params%write_time
        end do

    else
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Grid mode: create the mask on a given grid (do not create the grid here)
        ! Given is a list of files to take the grid from.
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        io_error = 0
        open(unit=14,file=grid_list, action='read', status='old')
        do while (io_error==0)
            read (14,'(A)',iostat=io_error) fname

            ! read in the file
            call readHDF5vct_tree( (/fname/), params, hvy_mask, tree_ID_flow, time=time)

            ! create the mask
            call createMask_tree(params, time, hvy_mask, hvy_mask, .false.)

            ! as we intend to create the mask on a given grid, pruning makes no sense

            ! store mask file
            write( fname,'(a, "_", i12.12, ".h5")') "mask", nint(time * 1.0e6_rk)
            call saveHDF5_tree(fname, time, -1_ik, 1, params, hvy_mask, tree_ID_flow, no_sync=.false.)

            if (save_us) then
                write( fname,'(a, "_", i12.12, ".h5")') "usx", nint(time * 1.0e6_rk)
                call saveHDF5_tree(fname, time, -1_ik, 2, params, hvy_mask, tree_ID_flow, no_sync=.false.)

                write( fname,'(a, "_", i12.12, ".h5")') "usy", nint(time * 1.0e6_rk)
                call saveHDF5_tree(fname, time, -1_ik, 3, params, hvy_mask, tree_ID_flow, no_sync=.false.)

                write( fname,'(a, "_", i12.12, ".h5")') "usz", nint(time * 1.0e6_rk)
                call saveHDF5_tree(fname, time, -1_ik, 4, params, hvy_mask, tree_ID_flow, no_sync=.false.)
            endif

            if (save_color) then
                write( fname,'(a, "_", i12.12, ".h5")') "color", nint(time * 1.0e6_rk)
                call saveHDF5_tree(fname, time, -1_ik, 5, params, hvy_mask, tree_ID_flow, no_sync=pruned)
            endif
        enddo
        close (14)
    endif



end subroutine post_dry_run
