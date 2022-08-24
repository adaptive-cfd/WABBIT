!##############################################################
!> This function reads a set of data into a specified tree.
!> The set of data can involve multiple files. Each file should
!> contain a quantity at the same snapshot (i.e. same iteration/time :: the same grid!).
!> If no data has been read in before, this function will allocate
!> the heavy data for you.
subroutine read_field2tree(params, fnames, N_files, tree_ID, hvy_block, verbosity)
    implicit none
    !-----------------------------------------------------------------
    type (type_params), intent(inout) :: params           !< params structure
    integer(kind=ik), intent(in)      :: N_files     !< number of fields/quantities to read
    character(len=*), intent(in)      :: fnames(N_files)  !< file names
    integer(kind=ik), intent(in)      :: tree_ID     !< number of the tree
    real(kind=rk), ALLOCATABLE, intent(inout)      :: hvy_block(:, :, :, :, :) !< heavy data array - block data
    logical, intent(in), optional :: verbosity !< if verbosity==True generates log output
    !-------------------------------- ---------------------------------
    integer(kind=ik)  :: iteration, dF, tc_length, dim, i, lgt_n_tmp, &
    rank, level, Bs(3), fsize, number_dense_blocks
    real(kind=rk)     :: time, domain(3)
    logical :: verbose = .true.

    if (present(verbosity)) verbose=verbosity

    ! set MPI parameter
    rank  = params%rank
    fsize = params%forest_size

    if (tree_ID <= tree_n) then
        ! the tree already exists: to overwrite it, we first delete the existing one
        call delete_tree(params, tree_ID)
        call createActiveSortedLists_forest(params)
    endif

    ! From any file (here: 1st in list), we read the essential parameters of the tree:
    ! most importantly Bs, dim, length of treecodes. NOTE: all trees in the forest must
    ! have the same Bs.
    call read_attributes(fnames(1), lgt_n_tmp, time, iteration, domain, Bs, &
    tc_length, dim, verbosity=verbose, periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)

    ! Check if one tree already exists in the forest:
    ! If it doesnt initialize some important parameters like the Block size Bs
    ! and the spatial dimension of the data.
    ! If it does chFeck if the parameters are suitable for the forest
    if (allocated(hvy_block)) then
        if ( any(params%Bs /= Bs) ) then
            ! NOTE: all trees in the forest must have the same Bs.
            call abort(100119, "ERROR: Trying to read a tree (grid) with Bs different from existing trees in forest!")
        endif

        if (params%max_treelevel < tc_length) call abort(10119, "maximal treelevel incompatible")
    else
        params%max_treelevel = tc_length
        params%dim = dim
        params%n_eqn = N_files
        params%N_fields_saved = N_files
        params%Bs = Bs
        params%domain_size = domain

        ! we have to allocate grid if this routine is called for the first time
        call allocate_forest(params, hvy_block)

        hvy_neighbor = -1
        hvy_n = 0
        lgt_n = 0 ! reset number of active light blocks
        tree_n = 0 ! reset number of trees in forest
    endif

    call readHDF5vct_tree(fnames, params, hvy_block, tree_id, verbosity=.true.)


    call createActiveSortedLists_forest(params)
    call updateNeighbors_tree(params, tree_ID)

    call sync_ghosts(params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID) )

end subroutine read_field2tree
