!> \brief Refine the mash (tag / interpolate / update lists)
!
!> \details This routine first sets the refinement flag for all blocks to +1
!! and then executes the refinement directly. Blocks that cannot be refined because they
!! are already on the finest allowed level are unaltered.
!!
!! As the grid changes, active lists and neighbor relations are updated, and load balancing
!! is applied.
!!
!! Note we assume, on call, that active lists / neighbor are up-to-date
!! input:    - params, light and heavy data \n
!! output:   - light and heavy data arrays
! ********************************************************************************************

subroutine refine_tree( params, hvy_block, indicator, tree_ID, error_OOM, check_full_tree)
    use module_indicators

    implicit none

    type (type_params), intent(inout) :: params                   !> user defined parameter structure
    real(kind=rk), intent(inout)   :: hvy_block(:, :, :, :, :)    !> heavy data array - block data
    character(len=*), intent(in)   :: indicator                   !> how to choose blocks for refinement
    integer(kind=ik), intent(in)   :: tree_ID
    logical, intent(out)   :: error_OOM !> Out-of-memory error, causes the main time loop to exit.
    logical, intent(in), optional :: check_full_tree  !> if true, checks if a full tree can be created in deterministic fashion


    ! cpu time variables for running time calculation
    real(kind=rk)                  :: t0, t1, t2, t_misc, test(1:size(hvy_block, 4))
    integer(kind=ik)               :: k, hvy_n_afterRefinement, lgt_id, hvy_id, mpierr
    real(kind=rk) :: norm(1:params%n_eqn)
    logical :: checkFullTree

    checkFullTree = .false.
    if (present(check_full_tree)) checkFullTree = check_full_tree

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas

    error_OOM = .false.

    ! start time
    t0 = MPI_Wtime()
    t_misc = 0.0_rk

    !> (a) loop over the blocks and set their refinement status.
    t1 = MPI_Wtime()
    call refinementIndicator_tree( params, hvy_block, tree_ID, indicator )
    call toc( "refine_tree (refinementIndicator_tree)", 141, MPI_Wtime()-t1 )


    !> (b) remove refinement flag for blocks that are on the finest level and thus
    !! cannot be refined anymore.
    t1 = MPI_Wtime()
    call respectJmaxJmin_tree( params, tree_ID )
    call toc( "refine_tree (respectJmaxJmin_tree)", 142, MPI_Wtime()-t1 )


    !> (c) ensure gradedness of mesh. If the refinement is done everywhere, there is
    !! no way gradedness can be damaged, so we skip the call in this case. However,
    !! in all other indicators, this step is very important.
    t1 = MPI_Wtime()
    if ( indicator /= "everywhere" ) then
      call ensureGradedness_tree( params, tree_ID )
    endif
    call toc( "refine_tree (ensureGradedness_tree)", 143, MPI_Wtime()-t1 )


    ! this (preventive) load balancing step ensures that AFTER the execution of refinement, 
    ! the load is balanced. Can be done only after flagging for refinement, of course. May not
    ! be perfectly balanced (+- 2**D blocks, because we can only distribute a block before refinement
    ! which means 2**D blocks after refinement), but is better than before.
    ! To ensure perfect balancing, a second step is done after the actual refinement.
    if (indicator /= "everywhere") then
        call balanceLoad_tree( params, hvy_block, tree_ID, balanceForRefinement=.true.)
    endif

    !---------------------------------------------------------------------------
    ! check if the refinement step will succeed or if we will run out of memory
    !---------------------------------------------------------------------------
    ! NOTE: in the simulation part, a second check is performed BEFORE this routine
    ! is called. It then aborts the time loop and dumps a backup for late resuming.
    ! If the following check fails, then this will be in postprocessing, and it will kill the code.
    call check_oom(params, tree_id, error_OOM, check_full_tree=.false., check_ref=.true.)
    ! byebye
    if (error_OOM) return
    !---------------------------------------------------------------------------


    !> (d) execute refinement, interpolate the new mesh. Blocks with status=+1 are refined, often
    !! that means all blocks are refined.
    t1 = MPI_Wtime()
    call refinement_execute_tree( params, hvy_block, tree_ID )
    call toc( "refine_tree (refinement_execute)", 144, MPI_Wtime()-t1 )


    !> (e) as the grid changed now with the refinement, we have to update the list of
    !! active blocks so other routines can loop just over these active blocks
    !! and do not have to ensure that the active list is up-to-date
    t2 = MPI_wtime()
    call updateMetadata_tree(params, tree_ID)
    t_misc = MPI_wtime() - t2


    !> (f) At this point the refinement is done. Since not all blocks are refined, namely only those
    !! that were not on Jmax, now, the load distribution of blocks may no longer
    !! be balanced, so we have to balance load now. EXCEPTION: if the flag force_maxlevel_dealiasing is set
    !! then we force blocks on Jmax to coarsen, even if their details are large. Hence, each refinement
    !! step is a true "everywhere". Then, there is no need for balancing, as all mpiranks automatically
    !! hold the same number of blocks, if they started with a balanced distribution (which is the
    !! output of adapt_tree.)
    !! This of course implies any other indicator than "everywhere" requires balancing here.
    if ((params%force_maxlevel_dealiasing .eqv. .false.) .or. (indicator/="everywhere")) then
        t1 = MPI_Wtime()
        call balanceLoad_tree( params, hvy_block, tree_ID )
        call toc( "refine_tree (balanceLoad_tree)", 145, MPI_Wtime()-t1 )
    endif

    ! check once again to see if we can create the full tree
    if (checkFullTree) then
        call check_oom(params, tree_id, error_OOM, check_full_tree=.true., check_ref=.false.)
    endif

    call toc( "refine_tree (lists+neighbors)", 146, t_misc )
    call toc( "refine_tree (TOTAL)", 140, MPI_wtime()-t0 )

end subroutine refine_tree



!> \brief Check if a refinement step will succeed or if we will run out of memory
subroutine check_oom(params, tree_id, error_OOM, check_full_tree, check_ref)
    implicit none

    type (type_params), intent(in) :: params                   !> user defined parameter structure
    integer(kind=ik), intent(in) :: tree_id                    !> tree ID to check
    !> Out-of-memory error, causes the main time loop to exit.
    logical, intent(out) :: error_OOM
    logical, intent(in), optional :: check_full_tree  !> if true, checks if a full tree can be created in deterministic fashion
    logical, intent(in), optional :: check_ref        !> if true, checks for blocks that will be refined

    integer(kind=tsize) :: treecode
    integer(kind=ik) :: hvy_n_afterRefinement, k, i_level, mpierr, lgt_id, hvy_id, level, digit, ref
    logical :: checkFullTree, checkRef
    integer(kind=ik), allocatable :: hvy_n_procs(:)

    checkFullTree = .false.
    if (present(check_full_tree)) checkFullTree = check_full_tree
    checkRef = .false.
    if (present(check_ref)) checkRef = check_ref

    ! Attention: This logic has to be in line with the logic in init_full_tree

    ! ToDo: Change loop over hvy_n to loop over lgt_n and omit allreduce call

    ! check if the refinement step will succeed or if we will run out of memory
    hvy_n_afterRefinement = 0
    do k = 1, hvy_n(tree_ID)
        ! loop over hvy is safer, because we can tell for each mpirank if it fails or not
        ! so we detect also problems arising from load-imbalancing.

        hvy_id = hvy_active(k, tree_ID)
        call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
        ref = lgt_block(lgt_id, IDX_REFINE_STS)
        level = lgt_block(lgt_id, IDX_MESH_LVL)

        ! add block itself
        if ( ref == +1 .and. level /= params%Jmax .and. check_ref) then
            ! this block will be refined, so it creates 2**D new blocks
            hvy_n_afterRefinement = hvy_n_afterRefinement + 2**params%dim
        else
            ! this block will not be refined, so it remains.
            hvy_n_afterRefinement = hvy_n_afterRefinement + 1
        endif

        ! if we check the full tree, we also need to add all blocks that might be spawned by this one
        if (checkFullTree) then
            treecode = get_tc(lgt_block( lgt_ID, IDX_TC_1:IDX_TC_2 ))
            do i_level = level, params%Jmin+1, -1
                digit = tc_get_digit_at_level_b(treecode, params%dim, i_level, params%Jmax)
                if (digit == 2**params%dim-1) then
                    hvy_n_afterRefinement = hvy_n_afterRefinement + 1
                else
                    exit
                endif
            enddo
        endif
    enddo

    ! oh-oh case.
    if ( hvy_n_afterRefinement > params%number_blocks ) then
        if (check_ref) then
            write(*,'("On rank:",i5," hvy_n=",i5," will refine to ",i5," but limit is ",i5)') &
            params%rank, hvy_n(tree_ID), hvy_n_afterRefinement, params%number_blocks
        elseif (checkFullTree) then
            write(*,'("On rank:",i5," hvy_n=",i5," will create full tree with ",i5," blocks but limit is ",i5)') &
            params%rank, hvy_n(tree_ID), hvy_n_afterRefinement, params%number_blocks
        endif
        ! return the error. in simulations, this causes us to jump to backup saving. in postprocessing,
        ! the code shall abort.
        error_OOM = .true.
    else
        error_OOM = .false.
    endif

    call MPI_ALLREDUCE(MPI_IN_PLACE, error_OOM, 1, MPI_LOGICAL, MPI_LOR, WABBIT_COMM, mpierr)

    if (error_OOM) then
        allocate(hvy_n_procs(1:params%number_procs), stat=mpierr)
        ! Gather hvy_n_afterRefinement from all ranks to rank 0
        call MPI_GATHER(hvy_n_afterRefinement, 1, MPI_INTEGER, hvy_n_procs, 1, MPI_INTEGER, 0, WABBIT_COMM, mpierr)
        if (params%rank == 0) then
            open(unit=99, file="oom_hvy_n_per_rank.txt", status="replace")
            write(99,'(A, i0)') "Rank, current hvy_n, expected hvy_n after refinement, limit is: ", params%number_blocks
            do k = 1, params%number_procs
                write(99,'(i0, 1x, i0, 1x, i0)') k, hvy_n(tree_ID), hvy_n_procs(k)
            enddo
            close(99)
        endif
        deallocate(hvy_n_procs)
    endif

end subroutine check_oom