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
    character(len=cshort)                   :: operator, help
    real(kind=rk), dimension(3)             :: domain
    integer(hsize_t), dimension(2)          :: dims_treecode
    integer(kind=ik)                        :: number_dense_blocks, Nb_file
    logical                                 :: verbose

    call get_command_argument(1, operator)
    call get_command_argument(2, help)


    if (params%rank==0) then
        write(*,'("     /¯¯¯\                       /¯¯¯\       /¯¯¯\       /¯¯¯\       /¯¯¯\  ")')
        write(*,'("    /     \     Unit Tests      /     \     /     \     /     \     /     \ ")')
        write(*,'("___/       \___________________/       \___/       \___/       \___/       \")')
        write(*, '("")')
        write(*, '(A)') " Wabbit postprocessing: Unit tests"
        write(*, '("")')
    endif

    ! does the user need help?
    if (help=='--help' .or. help=='--h' .or. help=="-h") then
        if (params%rank==0) then
            write(*,'(20("_/¯\"))')
            write(*, '(A)') " Executes one of the available unit test."
            write(*, '(A)') ""
            write(*, '(A)') "Example:"
            write(*, '(A)') "   ./wabbit-post --ghost-nodes-test --wavelet=CDF42 --memory=8.0GB"
            write(*, '(A)') ""
            write(*,'(20("_/¯\"))')
            write(*, '(A)') " Available unit tests:"
            write(*, '(A)') "   --ghost-nodes-test"
            write(*, '(A)') "   --refine-coarsen-test"
            write(*, '(A)') "   --wavelet-decomposition-unit-test"
            write(*, '(A)') "   --wavelet-decomposition-invertibility-test"
            write(*, '(A)') "   --sync-test"
            write(*, '(A)') "   --treecode-test"
            write(*,'(20("_/¯\"))')
            write(*, '(A)') "Parameters with default value:"
            write(*, '(A)') "   --wavelet=CDF44           - wavelet to be utilized"
            write(*, '(A)') "   --memory=8.0GB            - memory to initialize arrays"
            write(*, '(A)') "   --JMax=6                  - maximum block level"
            write(*, '(A)') "   --JMin=6                  - minimum block level"
            write(*, '(A)') "   --dim=2                   - dimension of test"
            write(*, '(A)') "   --Bs                      - Block size, default adapts to wavelet"
            write(*, '(A)') "   --g                       - Amount of ghost points, default adapts to wavelet"
            write(*, '(A)') "   --verbose=0               - prints and saves more debugging data for tests"
            write(*, '(A)') "   --max_grid_density=0.1    - Percentage of memory utilization to be targeted for random grids"
            write(*,'(20("_/¯\"))')
        end if
        return
    endif

    ! this routine works only on one tree
    allocate( hvy_n(1), lgt_n(1) )


    call get_cmd_arg( "--wavelet", params%wavelet, default="CDF44" )
    call get_cmd_arg( "--Jmax", params%Jmax, default=6 )
    call get_cmd_arg( "--Jmin", params%Jmin, default=1 )
    call get_cmd_arg( "--dim", params%dim, default=2 )
    call get_cmd_arg( "--Bs", Bs, default=-1 )
    call get_cmd_arg( "--verbose", verbose, default=.false.)
    call get_cmd_arg( "--max-grid-density", params%max_grid_density, 0.1_rk)

    ! initialize block size dynamically, make it BSmin for every wavelet
    if (Bs == -1) then
        ! unlifted wavelets, +4 per increase in X of CDFX0
        if (params%wavelet(4:5) == "20") Bs = 4
        if (params%wavelet(4:5) == "40") Bs = 8
        if (params%wavelet(4:5) == "60") Bs = 12

        ! lifted wavelets, +4 per increase in X or Y of CDFXY, commented are the sizes were leaf-first decomposition optimization is used, where we have +6
        if (params%wavelet(4:5) == "22") Bs = 6
        if (params%wavelet(4:5) == "24" .or. params%wavelet(4:5) == "42") Bs = 10  ! 12
        if (params%wavelet(4:5) == "26" .or. params%wavelet(4:5) == "44" .or. params%wavelet(4:5) == "62") Bs = 14  ! 18
        if (params%wavelet(4:5) == "28" .or. params%wavelet(4:5) == "46" .or. params%wavelet(4:5) == "64") Bs = 18  ! 224
        if (params%wavelet(4:5) == "66") Bs = 22  ! 30
    endif
    params%Bs(1:3) = 1
    params%Bs(1:params%dim) = Bs
    
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

    params%eps_normalized = .false.
    allocate(params%threshold_state_vector_component(1:params%n_eqn))
    params%threshold_state_vector_component(1:params%n_eqn) = 1
    params%coarsening_indicator = "threshold-state-vector"

    params%useCoarseExtension = params%isLiftedWavelet
    params%useSecurityZone = params%isLiftedWavelet
    ! params%useCoarseExtension = .true.
    ! params%useSecurityZone = .true.

    call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp, hvy_work=hvy_work, neqn_hvy_tmp=1, nrhs_slots1=1 )

    select case(operator)
    case("--ghost-nodes-test")
        call unit_test_ghostSync( params, hvy_block, hvy_work, hvy_tmp, tree_ID, abort_on_fail=.true., verbose=verbose)

    case("--refine-coarsen-test")
        call unit_test_refineCoarsen( params, hvy_block, hvy_work, hvy_tmp, tree_ID, verbose=verbose)

    case("--wavelet-decomposition-unit-test")
        call unit_test_waveletDecomposition( params, hvy_block, hvy_work, hvy_tmp, tree_ID, verbose=verbose )
    case("--wavelet-decomposition-invertibility-test")
        call unit_test_waveletDecomposition_invertibility( params, hvy_block, hvy_work, hvy_tmp, tree_ID, verbose=verbose )
    case("--sync-test")
        call unit_test_Sync( params, hvy_block, hvy_work, hvy_tmp, tree_ID, abort_on_fail=.true., verbose=verbose)
    case("--treecode-test")
        call unit_test_treecode( params, hvy_block, hvy_work, hvy_tmp, tree_ID, abort_on_fail=.true.)

    case default
        call abort(202355462,"unknown operator")
    end select

    call deallocate_forest(params, hvy_block, hvy_tmp=hvy_tmp)
end subroutine
