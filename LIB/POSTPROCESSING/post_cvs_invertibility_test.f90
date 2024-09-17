! JB ToDo

subroutine post_cvs_invertibility_test(params)
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
    character(len=cshort)  :: file_in
    character(len=cshort)  :: file_out
    real(kind=rk)          :: time, time_given, norm_1(1:1), norm_2(1:1), norm_ref1(1:1), norm_ref2(1:1)
    integer(kind=ik)       :: iteration

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_tmp(:, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_work(:, :, :, :, :, :)
    integer(kind=ik)                   :: tree_ID=1, hvy_id, lgtID, hvyID

    integer(kind=ik)                        :: level, k, tc_length
    integer(kind=ik)                        :: Bs, Jmin_diff, Nb_old
    integer(hid_t)                          :: file_id
    character(len=cshort)                   :: order
    character(len=cshort)                   :: fname1
    real(kind=rk), dimension(3)             :: domain
    integer(hsize_t), dimension(2)          :: dims_treecode
    integer(kind=ik)                        :: number_dense_blocks, Nb_file
    logical                                 :: verbose

    ! filename should follow directly after option --denoise
    call get_command_argument(2, fname1)

    ! this routine works only on one tree
    allocate( hvy_n(1), lgt_n(1) )

    call get_cmd_arg( "--wavelet", params%wavelet, default="CDF44" )
    call get_cmd_arg( "--Jmax", params%Jmax, default=9 )
    call get_cmd_arg( "--Jmin", params%Jmin, default=1 )
    call get_cmd_arg( "--dim", params%dim, default=2 )
    call get_cmd_arg( "--Bs", Bs, default=-1 )
    call get_cmd_arg( "--verbose", verbose, default=.false.)

    ! initialize block size dynamically, make it small so that tests dont take too long
    if (Bs == -1) Bs = 20
    params%Bs(1:3) = 1
    params%Bs(1:params%dim) = Bs

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

    call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp, hvy_work=hvy_work, neqn_hvy_tmp=1, nrhs_slots1=4 )

    ! **********************************
    ! first test: not adaptive
    ! read in data
    if (params%rank == 0) then
        write(*, '(A)') "CVS-TEST: First test - not adaptive"
    endif

    call readHDF5vct_tree((/fname1/), params, hvy_block, tree_ID=tree_ID_flow, verbosity=.true.)

    ! save first data
    do k = 1, hvy_n(tree_ID)
        hvy_ID = hvy_active(k, tree_ID)
        hvy_work(:, :, :, :, hvy_id, 3) = hvy_block(:, :, : ,:, hvy_id)
    enddo

    if (params%rank == 0) then
        write(*, '(A)') "CVS-TEST: Starting first adapt_tree loop for first thresholding"
    endif

    ! denoise data
    call adapt_tree( 0.0_rk, params, hvy_block, tree_ID, params%coarsening_indicator, hvy_tmp, hvy_work, ignore_coarsening=.true.)

    ! save second data
    do k = 1, hvy_n(tree_ID)
        hvy_ID = hvy_active(k, tree_ID)
        hvy_work(:, :, :, :, hvy_id, 4) = hvy_block(:, :, : ,:, hvy_id)
    enddo

    if (params%rank == 0) then
        write(*, '(A)') "CVS-TEST: Starting second adapt_tree loop for second thresholding"
    endif

    ! denoise data
    call adapt_tree( 0.0_rk, params, hvy_block, tree_ID, params%coarsening_indicator, hvy_tmp, hvy_work, ignore_coarsening=.true.)

    ! reference norms
    call componentWiseNorm_tree(params, hvy_work(:, :, :, :, :, 3), tree_ID, "L2", norm_ref1)
    call componentWiseNorm_tree(params, hvy_work(:, :, :, :, :, 4), tree_ID, "L2", norm_ref2)

    ! compute errors
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k,tree_ID)
        hvy_work(:, :, :, :, hvy_id, 3) = hvy_work(:, :, :, :, hvy_id, 4) - hvy_work(:, :, :, :, hvy_id, 3)
        hvy_work(:, :, :, :, hvy_id, 4) = hvy_block(:, :, : ,:, hvy_id) - hvy_work(:, :, :, :, hvy_id, 4)
    end do

    ! compare data - compute norm of error
    call componentWiseNorm_tree(params, hvy_work(:, :, :, :, :, 3), tree_ID, "L2", norm_1)
    call componentWiseNorm_tree(params, hvy_work(:, :, :, :, :, 4), tree_ID, "L2", norm_2)

    ! output to screen
    if (params%rank == 0) then
        write(*, '(A, es12.4)') "CVS-TEST: Rel error between initial field and field after first thresholding: ", norm_1/norm_ref1
        write(*, '(A, es12.4)') "CVS-TEST: Rel error between field after first and second thresholding:        ", norm_2/norm_ref2
    endif

    ! **********************************
    ! second test: adaptive
    ! read in data
    if (params%rank == 0) then
        write(*, '(A)') "CVS-TEST: Second test - adaptive"
    endif
    call readHDF5vct_tree((/fname1/), params, hvy_block, tree_ID=tree_ID_flow, verbosity=.true.)

    if (params%rank == 0) then
        write(*, '(A)') "CVS-TEST: Starting first adapt_tree loop for first thresholding and adaption"
    endif

    ! denoise data
    call adapt_tree( 0.0_rk, params, hvy_block, tree_ID, params%coarsening_indicator, hvy_tmp, hvy_work, ignore_coarsening=.false.)

    Nb_old = lgt_n(tree_ID)
    if (params%rank == 0) then
        write(*, '(A, i0)') "CVS-TEST: After first adapt_tree loop, Nb= ", lgt_n(tree_ID)
    endif

    ! save second data
    do k = 1, hvy_n(tree_ID)
        hvy_ID = hvy_active(k, tree_ID)
        hvy_work(:, :, :, :, hvy_id, 4) = hvy_block(:, :, : ,:, hvy_id)
    enddo

    if (params%rank == 0) then
        write(*, '(A)') "CVS-TEST: Starting second adapt_tree loop for second thresholding and adaption"
    endif

    ! denoise data
    call adapt_tree( 0.0_rk, params, hvy_block, tree_ID, params%coarsening_indicator, hvy_tmp, hvy_work, ignore_coarsening=.false.)

    if (params%rank == 0) then
        write(*, '(A, i0)') "CVS-TEST: After second adapt_tree loop, Nb= ", lgt_n(tree_ID)
    endif

    ! reference norms
    call componentWiseNorm_tree(params, hvy_work(:, :, :, :, :, 4), tree_ID, "L2", norm_ref2)

    ! compute errors
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k,tree_ID)
        hvy_work(:, :, :, :, hvy_id, 4) = hvy_block(:, :, : ,:, hvy_id) - hvy_work(:, :, :, :, hvy_id, 4)
    end do

    ! compare data - compute norm of error
    call componentWiseNorm_tree(params, hvy_work(:, :, :, :, :, 4), tree_ID, "L2", norm_2)

    ! output to screen
    if (params%rank == 0) then
        write(*, '(A, i0, " vs. ", i0)') "CVS-TEST: Blocks after first and second iteration: ", Nb_old, lgt_n(tree_ID)
        write(*, '(A, es12.4)') "CVS-TEST: Rel error between field after first and second thresholding with adaption: ", norm_2/norm_ref2
    endif

    call deallocate_forest(params, hvy_block, hvy_tmp=hvy_tmp)
end subroutine
