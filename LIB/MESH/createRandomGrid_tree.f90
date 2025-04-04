subroutine createRandomGrid_tree( params, hvy_block, hvy_tmp, level_init, verbosity, iterations, tree_ID )

    implicit none

    type (type_params), intent(inout)   :: params                     !> user defined parameter structure
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)   !> heavy data array - block data
    !> heavy temp data: used for saving, filtering, and helper qtys (reaction rate, mask function)
    real(kind=rk), intent(out)          :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik), intent(in)        :: level_init, iterations           !> what level to initialize?
    !> write output
    logical, intent(in)                 :: verbosity
    integer(kind=ik), intent(in)        :: tree_ID

    logical :: error_OOM 
    integer :: l

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas

    ! set all blocks to free (since if we call inicond twice, all blocks are used in the second call)
    call reset_tree( params, verbosity, tree_ID=tree_ID )

    ! setup the coarsest grid level with some data (we don't care what data, we'll erase it)
    ! Note that active lists + neighbor relations are updated inside this routine as well, as
    ! the grid is modified
    call createEquidistantGrid_tree( params, hvy_block, max(level_init, params%Jmin), verbosity, tree_ID=tree_ID )

    ! second: refine some blocks (random), coarsen some blocks (random)
    do l = 1, iterations
        if (params%rank==0 .and. verbosity) then
            write(*,'("RANDOM GRID GENERATION: iteration ",i2," active=",i9," J=(",i2,"/",i2,")")') &
            l, lgt_n(tree_ID), &
            minActiveLevel_tree(tree_ID), &
            maxActiveLevel_tree(tree_ID)
        endif

        ! randomly refine some blocks
        call refine_tree( params, hvy_block, hvy_tmp, "random", tree_ID=tree_ID, error_OOM=error_OOM   )

        if (error_OOM) call abort(2512111, "Out of memory, refinement failed. Try with more memory.")

        ! randomly coarsen some blocks
        call adapt_tree( 0.0_rk, params, hvy_block, tree_ID, "random", hvy_tmp  )
    end do

    if (params%rank==0 .and. verbosity) then
        write(*,'("DID RANDOM GRID GENERATION: active=",i9," J=(",i2,"/",i2,") density=",es12.4," / ",es12.4)') &
        lgt_n(tree_ID), &
        minActiveLevel_tree(tree_ID), &
        maxActiveLevel_tree(tree_ID), &
        dble(lgt_n(tree_ID)) / dble(size(lgt_block, 1)), params%max_grid_density
    endif

    call balanceLoad_tree( params, hvy_block, tree_ID)
end subroutine createRandomGrid_tree
