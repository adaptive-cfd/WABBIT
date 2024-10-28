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
    real(kind=rk)          :: time, time_given, domain(1:3)

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_tmp(:, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_work(:, :, :, :, :, :)
    integer(kind=ik)                   :: tree_ID=1, hvy_id, lgtID, hvyID, level, k, dim, Bs(3), iteration, lgt_n_tmp
    logical                            :: verbose, denoise_energy, save_noise

    ! help should follow directly after option --denoise
    call get_command_argument(2, file_in)

    ! does the user need help?
    if (file_in=='--help' .or. file_in=='--h' .or. file_in=='-h') then
        if (params%rank==0) then
            write(*,*) "------------------------------------------------------------------"
            write(*,*) "./wabbit-post --denoise --files=""FILES_IN"" --memory=[memory] [options]"
            write(*,*) "------------------------------------------------------------------"
            write(*,*) " This function denoises the input field(s) FILES_IN"
            write(*,*) " Further options with their defaults are:"
            write(*,*) "    --wavelet=CDF44"
            write(*,*) "    --Jmin=1"
            write(*,*) "    --verbose=0     : write more fields"
            write(*,*) "    --noise=0       : save noise field as well"
            write(*,*) "    --energy=0      : denoise fields by energy, needs input fields (ux uy (uz) p, ...)"
            write(*,*) "------------------------------------------------------------------"
        end if
        return
    endif

    ! this routine works only on one tree
    allocate( hvy_n(1), lgt_n(1) )

    call get_cmd_arg( "--wavelet", params%wavelet, default="CDF44" )
    call get_cmd_arg( "--Jmin", params%Jmin, default=1 )
    call get_cmd_arg( "--noise", save_noise, default=.false. )  ! is only correct if loadbalancing does nothing
    call get_cmd_arg( "--verbose", verbose, default=.false.)
    call get_cmd_arg( "--energy", denoise_energy, default=.false. )
    call get_cmd_arg( "--files", params%input_files )

    ! get some parameters from one of the files (they should be the same in all of them)
    call read_attributes(params%input_files(1), lgt_n(tree_ID_flow), time, iteration, params%domain_size, params%Bs, params%Jmax, params%dim, &
    periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)

    ! check if the other files are consistent with this one
    do k = 1, params%n_eqn
        file_in = params%input_files(k)
        call check_file_exists(trim(file_in))
        call read_attributes(file_in, lgt_n_tmp, time, iteration, domain, Bs, level, dim)

        params%Jmax = max(params%Jmax, level) ! find the maximal level of all snapshot

        if (lgt_n(1) .ne. lgt_n_tmp) call abort(241003,"Number of blocks do not agree!")
        if (any(params%Bs .ne. Bs)) call abort( 241003, "Block size is not consistent ")
        if (params%dim .ne. dim) call abort(241003,"Dimensions do not agree!")
        if ( abs(sum(params%domain_size(1:dim) - domain(1:dim))) > 1e-14 ) call abort( 241003, "Domain size is not consistent ")

        ! Concatenate "sparse" with filename
        params%input_files(k) = trim(file_in)
    end do

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
    params%n_eqn = size(params%input_files)
    allocate(params%field_names(params%n_eqn))
    allocate(params%threshold_state_vector_component(params%n_eqn))
    params%block_distribution = "sfc_hilbert"

    params%eps_normalized = .false.
    params%eps_norm = "L2"
    params%eps = 1.0e200_rk  ! practicallz disable it
    params%azzalini_iterations = 100
    params%threshold_wc = .true.
    params%threshold_state_vector_component(1:params%n_eqn) = .true.
    if (.not. denoise_energy) then
        params%coarsening_indicator = "threshold-image-denoise"
    else
        params%coarsening_indicator = "threshold-cvs"
    endif

    call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp, hvy_work=hvy_work, neqn_hvy_tmp=params%n_eqn, nrhs_slots1=2 )

    ! read in data
    call readHDF5vct_tree(params%input_files, params, hvy_block, tree_ID=tree_ID_flow, verbosity=.true.)

    ! save data as full tree
    if (verbose) then
        do k = 1, params%n_eqn
            file_in = trim(params%input_files(k))
            write(file_out, '(A, A, A)') file_in(1:index(file_in, "/", back=.true.)), "full-tree-WD-before-", file_in(index(file_in, "/", back=.true.)+1:len(trim(file_in)))
            call saveHDF5_wavelet_decomposed_tree( file_out, time, iteration, k, params, hvy_block, hvy_tmp, tree_ID)
        enddo
    endif

    ! denoise data
    call adapt_tree( time, params, hvy_block, tree_ID, params%coarsening_indicator, hvy_tmp, hvy_work, ignore_coarsening=.true.)

    do k = 1, params%n_eqn
        file_in = trim(params%input_files(k))
        write(file_out, '(A, A, A)') file_in(1:index(file_in, "/", back=.true.)), "denoised-", file_in(index(file_in, "/", back=.true.)+1:len(trim(file_in)))
        call saveHDF5_tree( file_out, time, iteration, k, params, hvy_block, tree_ID)
    enddo
    if (save_noise) then
        ! Noise will be computed during adapt_tree in hvy_tmp, however this is only correct without load balancing (nProcs=1) and without coarsening
        do k = 1, params%n_eqn
            file_in = trim(params%input_files(k))
            write(file_out, '(A, A, A)') file_in(1:index(file_in, "/", back=.true.)), "noise-", file_in(index(file_in, "/", back=.true.)+1:len(trim(file_in)))
            call saveHDF5_tree( file_out, time, iteration, k, params, hvy_tmp, tree_ID)
        enddo
    endif

    if (verbose) then
        do k = 1, params%n_eqn
            file_in = trim(params%input_files(k))
            write(file_out, '(A, A, A)') file_in(1:index(file_in, "/", back=.true.)), "full-tree-WD-", file_in(index(file_in, "/", back=.true.)+1:len(trim(file_in)))
            call saveHDF5_wavelet_decomposed_tree( file_out, time, iteration, k, params, hvy_block, hvy_tmp, tree_ID)
        enddo
    endif


    call deallocate_forest(params, hvy_block, hvy_tmp=hvy_tmp)
end subroutine
