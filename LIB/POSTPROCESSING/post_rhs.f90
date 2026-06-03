subroutine post_rhs(params)
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_operators
    use module_physics_metamodule
    use module_time_step
    use module_forestMetaData

    implicit none

    type (type_params), intent(inout)  :: params
    character(len=cshort)  :: fname_ini, dummy, fname
    character(len=cshort), allocatable  :: files(:)
    real(kind=rk)          :: time
    integer(kind=ik)       :: iteration, k, lgt_id, tc_length, g
    integer(kind=ik), dimension(3) :: Bs

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_work(:, :, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_tmp(:, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_mask(:, :, :, :, :)
    real(kind=rk), dimension(3)        :: domain
    integer(kind=ik)                   :: tree_ID=1

    !-----------------------------------------------------------------------------------------------------
    ! get values from command line (filename and level for interpolation)
    call get_command_argument(2, fname_ini)
    ! does the user need help?
    if (fname_ini=='--help' .or. fname_ini=='--h') then
        if (params%rank==0) then
            write(*,*) "------------------------------------------------------------------"
            write(*,*) "wabbit postprocessing routine for computation of RHS"
            write(*,*) "------------------------------------------------------------------"
            write(*,*) " ./wabbit-post --compute-rhs PARAMS.ini state1_0000.h5 [state2_0000.h5]"
            write(*,*) "    Give all statevector components as separate files (inicond section"
            write(*,*) "    in paramsfile is ignored)"
            write(*,*) " Output is written to "
            write(*,*) " rhs1_0000.h5, rhs2...... timestamp matches input timestamp"
            write(*,*) "------------------------------------------------------------------"
        end if
        return
    endif



    ! read ini-file and save parameters in struct
    call ini_file_to_params( params, fname_ini )
    ! have the pysics module read their own parameters
    call init_physics_modules( params, fname_ini, params%N_mask_components  )

    call setup_wavelet(params)

    allocate(files(1:params%n_eqn))

    ! this routine works only on one tree
    allocate( hvy_n(1), lgt_n(1) )

    ! get some parameters from one of the files (they should be the same in all of them)
    call get_command_argument(3, files(1))
    call read_attributes(files(1), lgt_n(tree_ID), time, iteration, domain, Bs, tc_length, params%dim)

    ! in usual parameter files, RK4 (or some other RK) is used an requires a lot of memory
    ! here we do not need that, and hence pretent to use a basic scheme (EE1 maybe)
    deallocate(params%butcher_tableau)
    allocate(params%butcher_tableau(1,1))

    params%number_blocks = ceiling( dble(lgt_n(tree_ID))/dble(params%number_procs) )
    Bs = params%Bs
    g = params%g

    ! allocate data
    call allocate_forest(params, hvy_block, hvy_work, hvy_tmp, hvy_mask)

    ! read in all components of statevector
    do k = 1, params%n_eqn
        call get_command_argument(2+k, files(k))
    enddo

    ! read input data
    call readHDF5vct_tree( files, params, hvy_block, tree_ID)

    ! create lists of active blocks (light and heavy data)
    ! update list of sorted nunmerical treecodes, used for finding blocks
    call updateMetadata_tree(params, tree_ID)

    ! balance the load
    call balanceLoad_tree(params, hvy_block, tree_ID)

    call sync_ghosts_tree(params, hvy_block, tree_ID )

    ! compute right hand side
    call RHS_wrapper(time, params, hvy_block, hvy_work(:,:,:,:,:,1), hvy_mask, hvy_tmp, tree_ID)

    ! save result to disk
    do k = 1, params%n_eqn
        write( fname,'("rhs",i1,"_", i6.6, i6.6, ".h5")') k, int(time+1.0e-12_rk, kind=ik), nint(max((time-int(time+1.0e-12_rk, kind=ik))*1.0e6_rk, 0.0_rk), kind=ik)

        call saveHDF5_tree(fname, time, iteration, k, params, hvy_work(:,:,:,:,:,1), tree_ID)
    enddo

end subroutine post_rhs
