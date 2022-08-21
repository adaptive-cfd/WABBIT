subroutine createRandomGrid_tree( params, lgt_block, hvy_block, hvy_tmp, hvy_neighbor, lgt_active, &
    lgt_n, lgt_sortednumlist, hvy_active, hvy_n, Jmin, verbosity, iterations, tree_ID )

    implicit none

    type (type_params), intent(inout)   :: params                     !> user defined parameter structure
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)            !> light data array
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)   !> heavy data array - block data
    !> heavy temp data: used for saving, filtering, and helper qtys (reaction rate, mask function)
    real(kind=rk), intent(out)          :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik), intent(inout)     :: hvy_neighbor(:,:)          !> heavy data array - neighbor data
    integer(kind=ik), intent(inout)     :: lgt_active(:,:)            !> list of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_n(:)                   !> number of active blocks (light data)
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), intent(inout)  :: lgt_sortednumlist(:,:,:)
    integer(kind=ik), intent(inout)     :: hvy_active(:,:)            !> list of active blocks (heavy data)
    integer(kind=ik), intent(inout)     :: hvy_n(:)                   !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: Jmin, iterations           !> what level to initialize?
    !> write output
    logical, intent(in)                 :: verbosity
    integer(kind=ik), intent(in)        :: tree_ID

    integer :: l

    ! set all blocks to free (since if we call inicond twice, all blocks are used in the second call)
    call reset_tree( params, lgt_block, lgt_active, &
    lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true., tree_ID=tree_ID )

    ! setup the coarsest grid level with some data (we don't care what data, we'll erase it)
    ! Note that active lists + neighbor relations are updated inside this routine as well, as
    ! the grid is modified
    call createEquidistantGrid_tree( params, lgt_block, hvy_neighbor, lgt_active, &
    lgt_n, lgt_sortednumlist, hvy_active, hvy_n, 2, .true., tree_ID=tree_ID )

    !---------------------------------------------------------------------------------------------
    ! second: refine some blocks (random), coarsen some blocks (random)
    do l = 1, iterations
        if (params%rank==0 .and. verbosity) then
            write(*,'("RANDOM GRID GENERATION: iteration ",i2," active=",i9," Jmin=",i2," Jmax=",i2)') &
            l, lgt_n, &
            minActiveLevel_tree( lgt_block, tree_ID, lgt_active, lgt_n ), &
            maxActiveLevel_tree( lgt_block, tree_ID, lgt_active, lgt_n )
        endif

        ! randomly refine some blocks
        call refine_tree( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, &
             lgt_sortednumlist, hvy_active, hvy_n, "random", tree_ID=tree_ID   )

        ! randomly coarsen some blocks
        call adapt_tree( 0.0_rk, params, lgt_block, hvy_block, hvy_neighbor, lgt_active, &
             lgt_n, lgt_sortednumlist, hvy_active, hvy_n, tree_ID, "random", hvy_tmp  )
    end do

    if (params%rank==0 .and. verbosity) then
        write(*,'("DID RANDOM GRID GENERATION: active=",i9," Jmin=",i2," Jmax=",i2," density=",es12.4," / ",es12.4)') &
        lgt_n(tree_ID), &
        minActiveLevel_tree( lgt_block, tree_ID, lgt_active, lgt_n ), &
        maxActiveLevel_tree( lgt_block, tree_ID, lgt_active, lgt_n ), &
        dble(lgt_n(tree_ID)) / dble(size(lgt_block, 1)), params%max_grid_density
    endif
end subroutine createRandomGrid_tree
