subroutine post_unit_test(params)
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_globals
    use module_forestMetaData

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

    integer(kind=ik)                        :: max_neighbors, level, k, tc_length
    integer(kind=ik), dimension(3)          :: Bs
    integer(hid_t)                          :: file_id
    character(len=cshort)                   :: order
    character(len=cshort)                   :: operator
    real(kind=rk), dimension(3)             :: domain
    integer(hsize_t), dimension(2)          :: dims_treecode
    integer(kind=ik)                        :: number_dense_blocks, Nb_file

    ! this routine works only on one tree
    allocate( hvy_n(1), lgt_n(1) )


    call get_cmd_arg( "--wavelet", params%wavelet, default="CFD44" )
    call get_cmd_arg( "--dim", params%dim, default=2 )
    ! initialize wavelet transform
    ! also, set number of ghost nodes params%G to minimal value for this wavelet
    call setup_wavelet(params, params%g)

    ! in postprocessing, it is important to be sure that the parameter struct is correctly filled:
    ! most variables are unfortunately not automatically set to reasonable values. In simulations,
    ! the ini files parser takes care of that (by the passed default arguments). But in postprocessing
    ! we do not read an ini file, so defaults may not be set.
    allocate(params%butcher_tableau(1,1))
    ! we read only one datafield in this routine
    params%n_eqn = 1
    params%block_distribution = "sfc_hilbert"
    params%Jmax = 5
    params%Jmin = 1
    params%Bs = 22
    params%domain_size = 1.0_rk

    call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp, hvy_work=hvy_work, neqn_hvy_tmp=1, nrhs_slots1=1 )

    call get_command_argument(1, operator)

    select case(operator)
    case("--ghost-nodes-test")
        call unitTest_ghostSync( params, hvy_block, hvy_work, hvy_tmp, tree_ID, abort_on_fail=.true.)

    case("--refine-coarsen-test")
        call unitTest_refineCoarsen( params, hvy_block, hvy_work, hvy_tmp, tree_ID )

    case("--wavelet-decomposition-unit-test")
        call unitTest_waveletDecomposition( params, hvy_block, hvy_work, hvy_tmp, tree_ID )

    case default
        call abort(202355462,"unknown operator")
    end select

    call deallocate_forest(params, hvy_block, hvy_tmp=hvy_tmp)
end subroutine
