subroutine post_denoising(params)
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_globals
    use module_forestMetaData
    use module_unit_test

    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=cshort)  :: file_in, file_out
    real(kind=rk)          :: time, time_given
    integer(kind=ik)       :: iteration

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_tmp(:, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_work(:, :, :, :, :, :)
    integer(kind=ik)                   :: tree_ID=1, hvy_id, lgtID, hvyID, level, k
    logical                            :: verbose

    ! filename should follow directly after option --denoise
    call get_command_argument(2, file_in)

    ! does the user need help?
    if (file_in=='--help' .or. file_in=='--h' .or. file_in=='-h') then
        if (params%rank==0) then
            write(*,*) "------------------------------------------------------------------"
            write(*,*) "./wabbit-post --denoise FILE_IN --memory=[memory] [options]"
            write(*,*) "------------------------------------------------------------------"
            write(*,*) " This function denoises the input field FILE_IN"
            write(*,*) " Further options are:"
            write(*,*) "    --wavelet, --Jmin"
            write(*,*) "    --verbose : write more fields"
            write(*,*) "------------------------------------------------------------------"
        end if
        return
    endif

    ! this routine works only on one tree
    allocate( hvy_n(1), lgt_n(1) )

    params%cvs = .True.

    call get_cmd_arg( "--wavelet", params%wavelet, default="CDF44" )
    call get_cmd_arg( "--Jmin", params%Jmin, default=1 )
    call get_cmd_arg( "--verbose", verbose, default=.false.)

        ! get some parameters from one of the files (they should be the same in all of them)
    call read_attributes(file_in, lgt_n(tree_ID_flow), time, iteration, params%domain_size, params%Bs, params%Jmax, params%dim, &
    periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)

    ! initialize wavelet transform
    ! also, set number of ghost nodes params%G to minimal value for this wavelet
    params%g = -1
    if (params%g == -1) then
        call setup_wavelet(params, g_wavelet=params%g, g_RHS=params%g_RHS)
    else
        call setup_wavelet(params, g_RHS=params%g_RHS)
    endif

    ! in postprocessing, it is important to be sure that the parameter struct is correctly filled:
    ! most variables are unfortunately not automatically set to reasonable values. In simulations,
    ! the ini files parser takes care of that (by the passed default arguments). But in postprocessing
    ! we do not read an ini file, so defaults may not be set.
    allocate(params%butcher_tableau(1,1))
    ! we read only one datafield in this routine
    params%n_eqn = 1
    params%block_distribution = "sfc_hilbert"
    ! params%Jmin = 1
    params%domain_size = 1.0_rk

    params%eps_normalized = .false.
    params%eps_norm = "L2"
    allocate(params%threshold_state_vector_component(1:params%n_eqn))
    params%threshold_state_vector_component(1:params%n_eqn) = .true.
    params%coarsening_indicator = "threshold-image-denoise"

    call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp, hvy_work=hvy_work, neqn_hvy_tmp=1, nrhs_slots1=2 )

    ! read in data
    call readHDF5vct_tree((/file_in/), params, hvy_block, tree_ID=tree_ID_flow, verbosity=.true.)

    ! save data as full tree
    if (verbose) then
        call saveHDF5_wavelet_decomposed_tree( "field-before-full-tree-WD_000000.h5", 0.0_rk, 0, 1, params, hvy_block, hvy_tmp, tree_ID)
    endif

    ! denoise data
    call adapt_tree( time, params, hvy_block, tree_ID, params%coarsening_indicator, hvy_tmp, hvy_work, ignore_coarsening=.true.)

    call saveHDF5_tree( "denoised_000000.h5", 0.0_rk, 0, 1, params, hvy_block, tree_ID)
    if (verbose) then
        ! Noise will be computed during adapt_tree in hvy_tmp, however this is only correct without load balancing (nProcs=1) and without coarsening
        call saveHDF5_tree( "noise_000000.h5", 0.0_rk, 0, 1, params, hvy_tmp, tree_ID)
    endif

    if (verbose) then
        call saveHDF5_wavelet_decomposed_tree( "field-full-tree-WD_000000.h5", 0.0_rk, 0, 1, params, hvy_block, hvy_tmp, tree_ID)
    endif


    call deallocate_forest(params, hvy_block, hvy_tmp=hvy_tmp)
end subroutine
