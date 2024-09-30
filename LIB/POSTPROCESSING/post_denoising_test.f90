subroutine post_denoising_test(params)
    use module_globals
    use module_helpers
    use module_mesh
    use module_params
    use module_mpi
    use module_globals
    use module_forestMetaData
    use module_unit_test

    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=cshort)  :: file_in, file_out, n_type
    real(kind=rk)          :: n_s, n_f, n_i, signal_strength(1:1), r1, r2, t_den, std_est, time

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_tmp(:, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_work(:, :, :, :, :, :)
    integer(kind=ik)                   :: hvy_id, lgtID, hvyID, ix, iy, iz, iteration

    integer(kind=ik)                        :: level, k, tc_length, Bs, Jmin_diff, n_n, i_n, tree_ID_backup=4, random_seed_size
    integer(hid_t)                          :: file_id
    integer(hsize_t), dimension(2)          :: dims_treecode
    real(kind=rk)                           :: norm_L2(1:1), norm_Linfty(1:1)
    logical                                 :: verbose, predictable, save_results

    ! this routine works on two trees - one for copying and one for data to work with
    allocate( hvy_n(1:tree_ID_backup), lgt_n(1:tree_ID_backup) )

    ! filename should follow directly after option --denoise
    call get_command_argument(2, file_in)

    ! does the user need help?
    if (file_in=='--help' .or. file_in=='--h' .or. file_in=='-h') then
        if (params%rank==0) then
            write(*,'(A)') "------------------------------------------------------------------"
            write(*,'(A)') "./wabbit-post --denoising-test FILE_IN --memory=[memory] [options]"
            write(*,'(A)') "------------------------------------------------------------------"
            write(*,'(A)') " This function does denoising study on input field FILE_IN"
            write(*,'(A)') " Further options are:"
            write(*,'(A)') "    --n-s=1e01        - lowest noise to be added and denoised"
            write(*,'(A)') "    --n-f=1e01        - highest noise to be added and denoised"
            write(*,'(A)') "    --n-n=1           - number of noises to be investigated"
            write(*,'(A)') "    --n-type=uniform  - noise type - ""uniform"" or ""gaussian"""
            write(*,'(A)') "    --save=0          - Save results of each iteration"
            write(*,'(A)') "    --wavelet=CDF44   - wavelet used for denoising study"
            write(*,'(A)') "    --Jmin=1          - minimum level for full wavelet transform"
            write(*,'(A)') "    --verbose=0       - write more fields"
            write(*,'(A)') "------------------------------------------------------------------"
        end if
        return
    endif

    ! read in noise test parameters
    call get_cmd_arg( "--n-s", n_s, default=10.0_rk)
    call get_cmd_arg( "--n-f", n_f, default=10.0_rk)
    call get_cmd_arg( "--n-n", n_n, default=1)
    call get_cmd_arg( "--n-type", n_type, default="uniform")
    call get_cmd_arg( "--save", save_results, default=.false.)

    ! read in wabbit parameters
    call get_cmd_arg( "--wavelet", params%wavelet, default="CDF44" )
    call get_cmd_arg( "--Jmin", params%Jmin, default=1 )
    call get_cmd_arg( "--verbose", verbose, default=.false.)

    ! make random values predictable for each noise value or set after system clock
    ! helps for comparing between different wavelets
    predictable = .false.

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

    params%eps_normalized = .false.
    params%eps_norm = "L2"
    allocate(params%threshold_state_vector_component(1:params%n_eqn))
    params%threshold_state_vector_component(1:params%n_eqn) = .true.
    params%coarsening_indicator = "threshold-image-denoise"
    params%forest_size = tree_ID_backup

    call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp, hvy_work=hvy_work, neqn_hvy_tmp=1, nrhs_slots1=3 )

    ! read in data
    call readHDF5vct_tree((/file_in/), params, hvy_block, tree_ID=tree_ID_backup, verbosity=.true.)

    ! save data as full tree
    if (verbose) then
        call saveHDF5_wavelet_decomposed_tree( "image-full-tree-WD_0.h5", 0.0_rk, 0, 1, params, hvy_block, hvy_tmp, tree_ID_backup)
    endif

    ! loop over each noise to be investigated
    do i_n = 0, n_n
        ! noise variance at this iteration, computed in logspace between other noises
        n_i = exp(log(n_s) + (log(n_f) - log(n_s)) * dble(i_n)/dble(n_n))

        ! copy tree
        call copy_tree(params, hvy_block, tree_ID_flow, tree_ID_backup, skip_sync_ghosts=.true.)

        ! add noise
        if (predictable) then
            call random_seed(size=random_seed_size)
            call random_seed(put=int(log(n_i)*1000.0_rk)* (/ (k-1, k=1, random_seed_size) /))
        else
            call init_random_seed()
        endif
        do k = 1, hvy_n(tree_ID_flow)
            hvy_id = hvy_active(k, tree_ID_flow)
            do ix = params%g+1, params%Bs(1)+2*params%g
                do iy = params%g+1, params%Bs(2)+2*params%g
                    do iz = merge(1, params%g+1, params%dim==2), merge(1, params%Bs(3)+2*params%g, params%dim==2)
                        ! add random data as noise
                        call random_number(r1)
                        if (n_type == "uniform") then
                            hvy_block(ix,iy,iz,1,hvy_id) = hvy_block(ix,iy,iz,1,hvy_id) + r1 * n_i * sqrt(12.0_rk)
                        elseif (n_type == "gaussian") then
                            call random_number(r2)
                            hvy_block(ix,iy,iz,1,hvy_id) = hvy_block(ix,iy,iz,1,hvy_id) + sqrt(-2.0_rk*log(r1)) * cos(8.0_rk * atan(1.0_rk)*r2) * n_i
                        else
                            call abort(240917, "What's that noise?")
                        endif
                    enddo
                enddo
            enddo
        enddo
        if (save_results) then
            write(file_out, '(A, i4.4, A)') "image-noise-", i_n, "_0.h5"
            call saveHDF5_tree( file_out, 0.0_rk, 0, 1, params, hvy_block, tree_ID_flow)
        endif

        ! compute signal strength for signal2noise ration
        call componentWiseNorm_tree(params, hvy_block, tree_ID_flow, "CVS", signal_strength(1:1))
        signal_strength(1:1) = signal_strength(1:1) / dble(2_ik**(params%dim*maxActiveLevel_tree(tree_ID_flow))*product(params%bs(1:params%dim)))

        ! denoise data
        t_den = MPI_Wtime()
        call adapt_tree( time, params, hvy_block, tree_ID_flow, params%coarsening_indicator, hvy_tmp, hvy_work, ignore_coarsening=.true., std_est=std_est)
        t_den = MPI_Wtime()-t_den

        if (save_results) then
            call saveHDF5_tree( "denoised_000000.h5", 0.0_rk, 0, 1, params, hvy_block, tree_ID_flow)
            if (verbose) then
                ! Noise will be computed during adapt_tree in hvy_tmp, however this is only correct without load balancing (nProcs=1) and without coarsening
                call saveHDF5_tree( "noise_000000.h5", 0.0_rk, 0, 1, params, hvy_tmp, tree_ID_flow)
            endif
        endif

        ! compare to original data:
        call substract_two_trees(params, hvy_block, hvy_tmp, tree_ID_flow, tree_ID_backup)
        call componentWiseNorm_tree(params, hvy_block, tree_ID_flow, "L2", norm_L2)
        call componentWiseNorm_tree(params, hvy_block, tree_ID_flow, "Linfty", norm_Linfty)

        ! print infos on this iteration
        if (params%rank == 0) then
            write(*, '(A, i4, 6(A, es10.3))') "It", i_n+1, ": Std in / est - ", n_i, " / ", std_est, &
                " , Err L2 / Linfty - ", norm_L2(1:1), " / ", norm_Linfty(1:1), &
                " , S2N - ", signal_strength(1) / n_i**2, " , Time - ", t_den 
        endif

        if (verbose) then
            call saveHDF5_wavelet_decomposed_tree( "image-denoised-full-tree-WD_0.h5", 0.0_rk, 0, 1, params, hvy_block, hvy_tmp, tree_ID_flow)
        endif

    enddo


    call deallocate_forest(params, hvy_block, hvy_tmp=hvy_tmp)
end subroutine
