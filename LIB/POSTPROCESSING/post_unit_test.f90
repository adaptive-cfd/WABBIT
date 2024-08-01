subroutine post_unit_test(params)
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
    real(kind=rk)          :: time, time_given
    integer(kind=ik)       :: iteration

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_tmp(:, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_work(:, :, :, :, :, :)
    integer(kind=ik)                   :: tree_ID=1, hvy_id, lgtID, hvyID

    integer(kind=ik)                        :: level, k, tc_length
    integer(kind=ik)                        :: Bs, Jmin_diff
    integer(hid_t)                          :: file_id
    character(len=cshort)                   :: order
    character(len=cshort)                   :: operator
    real(kind=rk), dimension(3)             :: domain
    integer(hsize_t), dimension(2)          :: dims_treecode
    integer(kind=ik)                        :: number_dense_blocks, Nb_file
    logical                                 :: verbose

    call get_command_argument(1, operator)

    ! this routine works only on one tree
    allocate( hvy_n(1), lgt_n(1) )


    call get_cmd_arg( "--wavelet", params%wavelet, default="CFD44" )
    call get_cmd_arg( "--Jmax", params%Jmax, default=5 )
    call get_cmd_arg( "--Jmin", params%Jmin, default=1 )
    call get_cmd_arg( "--dim", params%dim, default=2 )
    call get_cmd_arg( "--Bs", Bs, default=-1 )
    call get_cmd_arg( "--verbose", verbose, default=.false.)

    ! initialize block size dynamically, make it small so that tests dont take too long
    if (Bs == -1) then
        ! check for X in CDFXY
        if (params%wavelet(4:4) == "2") Bs = 6
        if (params%wavelet(4:4) == "4") Bs = 12
        if (params%wavelet(4:4) == "6") Bs = 20

        ! check for Y in CDFXY
        if (params%wavelet(5:5) == "0") Bs = Bs + 0
        if (params%wavelet(5:5) == "2") Bs = Bs + 4
        if (params%wavelet(5:5) == "4") Bs = Bs + 10
        if (params%wavelet(5:5) == "6") Bs = Bs + 16
        if (params%wavelet(5:5) == "8") Bs = Bs + 22
    endif
    params%Bs(1:3) = 1
    params%Bs(1:params%dim) = Bs

    ! ghost sync test needs at least 32 points over the domain length. We need to ensure that Jmin is fitted accordingly
    if (operator == "--ghost-nodes-test") then
        ! -0.1 to ensure integer cast is done correctly
        Jmin_diff = int(log(32.0/(real(Bs)-0.1)) / log(2.0))+1 - params%Jmin
        ! now increase both Jmin and Jmax accordingly
        params%Jmax = params%Jmax + Jmin_diff
        params%Jmin = params%Jmin + Jmin_diff

        write(*, '(A, i0, A, i0)') "UNIT TEST: Need atleast 32 points over domain length. Adapted Jmin = ", params%Jmin, " and Jmax = ", params%Jmax
    endif

    
    ! initialize wavelet transform
    ! also, set number of ghost nodes params%G to minimal value for this wavelet
    call get_cmd_arg( "--g", params%g, default=-1 )
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

    call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp, hvy_work=hvy_work, neqn_hvy_tmp=1, nrhs_slots1=1 )

    select case(operator)
    case("--ghost-nodes-test")
        call unitTest_ghostSync( params, hvy_block, hvy_work, hvy_tmp, tree_ID, abort_on_fail=.true., verbose=verbose)

    case("--refine-coarsen-test")
        call unitTest_refineCoarsen( params, hvy_block, hvy_work, hvy_tmp, tree_ID, verbose=verbose)

    case("--wavelet-decomposition-unit-test")
        call unitTest_waveletDecomposition( params, hvy_block, hvy_work, hvy_tmp, tree_ID )
    case("--sync-test")
        call unitTest_Sync( params, hvy_block, hvy_work, hvy_tmp, tree_ID, abort_on_fail=.true., verbose=verbose)
    case("--treecode-test")
        call unit_test_treecode( params, hvy_block, hvy_work, hvy_tmp, tree_ID, abort_on_fail=.true.)

    case default
        call abort(202355462,"unknown operator")
    end select

    call deallocate_forest(params, hvy_block, hvy_tmp=hvy_tmp)
end subroutine
