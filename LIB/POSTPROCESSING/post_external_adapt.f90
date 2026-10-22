subroutine external_adapt(params)
    use module_params
    use module_mesh
    use module_mpi
    use module_globals
    use module_helpers
    use module_forestMetaData

    implicit none

    type(type_params), intent(inout) :: params

    character(len=cshort) :: operator, file_in, file_out
    real(kind=rk) :: time
    integer(kind=ik) :: iteration
    logical :: verbose = .true.
    logical :: error_OOM

    real(kind=rk), allocatable :: hvy_block(:, :, :, :, :), hvy_tmp(:, :, :, :, :)
    integer(kind=ik) :: tree_ID

    ! this routine works only on one tree
    allocate( hvy_n(1), lgt_n(1) )

    ! argument parsing
    call get_command_argument(1, operator)
    call get_command_argument(2, file_in)
    call get_command_argument(3, file_out)

    if (file_in == '--help' .or. file_in == '--h' .or. len_trim(operator) == 0) then
        if (params%rank == 0) then
            write(*,'(A)') "--------------------------------------------------------------"
            write(*,'(A)') "           EXTERNAL ADAPT/REFINE POSTPROCESSING"
            write(*,'(A)') "--------------------------------------------------------------"
            write(*,'(A)') "Use an external refinement status saved in the HDF5 file."
            write(*,'(A)') "Commands:"
            write(*,'(A)') "  ./wabbit-post --external-refine source.h5 target.h5 [--wavelet=CDF40]"
            write(*,'(A)') "  ./wabbit-post --external-coarsen source.h5 target.h5 [--wavelet=CDF40]"
            write(*,'(A)') "Notes: refinement status must be stored in dataset 'refinement_status'"
            write(*,'(A)') "--------------------------------------------------------------"
        endif
        return
    end if

    call check_file_exists(trim(file_in))
    call read_attributes(file_in, lgt_n(1), time, iteration, params%domain_size, params%Bs, params%Jmax, params%dim, &
        periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)
    params%Jmin = 0

    call get_cmd_arg( "--wavelet", params%wavelet, default="CDF40" )
    ! initialize wavelet transform
    ! also, set number of ghost nodes params%G to minimal value for this wavelet
    call setup_wavelet(params, params%g)

    ! minimal forest setup
    params%n_eqn = 1
    params%block_distribution = "sfc_hilbert"
    allocate(params%butcher_tableau(1,1))

    ! allocate forest arrays
    call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp)

    tree_ID = 1

    ! read file and request refinement status to be loaded into light data
    call readHDF5vct_tree( (/trim(file_in)/), params, hvy_block, tree_ID, verbosity=verbose, read_refinement_status=.true. )

    ! ensure metadata and balance
    call balanceLoad_tree(params, hvy_block, tree_ID)

    ! choose action
    if (operator == '--external-refine') then
        ! refine using the external flags set in lgt_block(:, IDX_REFINE_STS)
        call refine_tree(params, hvy_block, 'nothing (external)', tree_ID, error_OOM)
        if (error_OOM) call abort(2512181, "Refinement failed, out of memory. Try with more memory.")

    elseif (operator == '--external-coarsen') then
        ! coarsen according to external flags
        call adapt_tree(0.0_rk, params, hvy_block, tree_ID, 'nothing (external)', hvy_tmp)

    else
        if (params%rank == 0) write(*,'(A)') "ERROR: Unknown operator. Use --external-refine or --external-coarsen"
        call deallocate_forest(params, hvy_block, hvy_tmp=hvy_tmp)
        return
    endif

    ! save result
    if (len_trim(file_out) == 0) then
        call abort(0909191, "You must specify a name for the target! See --external-refine --help")
    endif

    call saveHDF5_tree(file_out, time, iteration, 1, params, hvy_block, tree_ID)

    if (params%rank == 0) then
        write(*,'(A,1x,A)') "Wrote adapted file to", trim(adjustl(file_out))
    endif

    call deallocate_forest(params, hvy_block, hvy_tmp=hvy_tmp)

end subroutine external_adapt
