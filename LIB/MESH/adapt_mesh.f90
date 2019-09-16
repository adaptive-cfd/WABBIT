!! \author   msr, engels
!
!> \brief This routine performs the coarsing of the mesh, where possible. For the given mesh
!! we compute the details-coefficients on all blocks. If four sister blocks have maximum
!! details below the specified tolerance, (so they are insignificant), they are merged to
!! one coarser block one level below. This process is repeated until the grid does not change
!! anymore.
!!
!! As the grid changes, active lists and neighbor relations are updated, and load balancing
!! is applied.
!
!> \note The block thresholding is done with the restriction/prediction operators acting on the
!! entire block, INCLUDING GHOST NODES. Ghost node syncing is performed in threshold_block.
!
!> \note It is well possible to start with a very fine mesh and end up with only one active
!! block after this routine. You do *NOT* have to call it several times.
!
!> \details
!! \date 10/11/16 - switch to v0.4
!> \image html adapt_mesh.svg width=400

subroutine adapt_mesh( time, params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, &
    lgt_sortednumlist, hvy_active, hvy_n, tree_ID, indicator, hvy_tmp, hvy_mask, external_loop)

    implicit none

    real(kind=rk), intent(in)           :: time
    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)
    !> heavy data array
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy work data array - block data.
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    ! mask data. we can use different trees (4est module) to generate time-dependent/indenpedent
    ! mask functions separately. This makes the mask routines tree-level routines (and no longer
    ! block level) so the physics modules have to provide an interface to create the mask at a tree
    ! level. All parts of the mask shall be included: chi, boundary values, sponges.
    ! Optional: if the grid is not adapted to the mask, passing hvy_mask is not required.
    real(kind=rk), intent(inout), optional :: hvy_mask(:, :, :, :, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(inout)     :: hvy_neighbor(:,:)
    !> list of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_n
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), intent(inout)  :: lgt_sortednumlist(:,:)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(inout)     :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(inout)     :: hvy_n
    !> coarsening indicator
    character(len=*), intent(in)        :: indicator
    !> Well, what now. The grid coarsening is an iterative process that runs until no more blocks can be
    !! coarsened. One iteration is not enough. If called without "external_loop", this routine
    !! performs this loop until it is converged. In some situations, this might be undesired, and
    !! the loop needs to be outsourced to the calling routine. This happens currently (07/2019)
    !! only in the initial condition, where the first grid can be so coarse that the inicond is different
    !! only on one point, which is then completely removed (happens for a mask function, for example.)
    !! if external_loop=.true., only one iteration step is performed.
    logical, intent(in), optional       :: external_loop

    integer(kind=ik), intent(in)        :: tree_ID
    ! loop variables
    integer(kind=ik)                    :: lgt_n_old, iteration, k, max_neighbors, lgt_id
    ! cpu time variables for running time calculation
    real(kind=rk)                       :: t0, t1
    ! MPI error variable
    integer(kind=ik)                    :: ierr, k1, hvy_id


    ! start time
    t0 = MPI_Wtime()
    t1 = t0
    lgt_n_old = 0
    iteration = 0

    if ( params%dim == 3 ) then
        max_neighbors = 56
    else
        max_neighbors = 12
    end if


    ! 2D case:
    !       |          |         |
    !   1   |    2     |    3    |  4
    !       |          |         |
    ! -----------------------------------
    !       |                    |
    !   5   |                    |  6
    !       |                    |
    ! ------|       my_rank      |-------
    !       |                    |
    !   7   |                    |  8
    !       |                    |
    ! -----------------------------------
    !       |          |         |
    !   9   |    10    |    11   |  12
    !       |          |         |

    ! To avoid that the incomming hvy_neighbor array and active lists are outdated
    ! we synchronice them.
    t0 = MPI_Wtime()
    call update_grid_metadata(params, lgt_block, hvy_neighbor, lgt_active, lgt_n, &
    lgt_sortednumlist, hvy_active, hvy_n, tree_ID)
    call toc( "adapt_mesh (update neighbors)", MPI_Wtime()-t0 )

    !> we iterate until the number of blocks is constant (note: as only coarsening
    !! is done here, no new blocks arise that could compromise the number of blocks -
    !! if it's constant, its because no more blocks are coarsened)
    do while ( lgt_n_old /= lgt_n )
        lgt_n_old = lgt_n

        !> (a) check where coarsening is possible
        ! ------------------------------------------------------------------------------------
        ! first: synchronize ghost nodes - thresholding on block with ghost nodes
        ! synchronize ghostnodes, grid has changed, not in the first one, but in later loops
        t0 = MPI_Wtime()
        call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )
        call toc( "adapt_mesh (sync_ghosts)", MPI_Wtime()-t0 )

        !! calculate detail on the entire grid. Note this is a wrapper for block_coarsening_indicator, which
        !! acts on a single block only
        t0 = MPI_Wtime()
        if (params%threshold_mask .and. present(hvy_mask)) then
            ! if present, the mask can also be used for thresholding (and not only the state vector). However,
            ! as the grid changes within this routine, the mask will have to be constructed in grid_coarsening_indicator
            call grid_coarsening_indicator( time, params, lgt_block, hvy_block, hvy_tmp, lgt_active, lgt_n, &
            lgt_sortednumlist, hvy_active, hvy_n, indicator, iteration, hvy_neighbor, hvy_mask)
        else
            call grid_coarsening_indicator( time, params, lgt_block, hvy_block, hvy_tmp, lgt_active, lgt_n, &
            lgt_sortednumlist, hvy_active, hvy_n, indicator, iteration, hvy_neighbor)
        endif
        call toc( "adapt_mesh (grid_coarsening_indicator)", MPI_Wtime()-t0 )


        !> (b) check if block has reached maximal level, if so, remove refinement flags
        t0 = MPI_Wtime()
        call respect_min_max_treelevel( params, lgt_block, lgt_active, lgt_n )
        call toc( "adapt_mesh (respect_min_max_treelevel)", MPI_Wtime()-t0 )

        !> (c) unmark blocks that cannot be coarsened due to gradedness and completeness
        t0 = MPI_Wtime()
        call ensure_gradedness( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, &
        lgt_sortednumlist, hvy_active, hvy_n )
        call toc( "adapt_mesh (ensure_gradedness)", MPI_Wtime()-t0 )

        !> (d) adapt the mesh, i.e. actually merge blocks
        t0 = MPI_Wtime()
        if (params%threshold_mask .and. present(hvy_mask)) then
            ! if the mask function is used as secondary coarsening criterion, we also pass the mask data array
            ! the idea is that now coarse-mesh will keep both hvy_block and hvy_mask on the same grid, i.e.
            ! the same coarsening is applied to both. the mask does not have to be re-created here, because
            ! regions with sharp gradients (that's where the mask is interesting) will remain unchanged
            ! by the coarsening function.
            ! This is not entirely true: often, adapt_mesh is called after refine_mesh, which pushes the grid
            ! to Jmax. If the dealiasing function is on (params%force_maxlevel_dealiasing=.true.), the coarsening
            ! has to go down one level. If the mask is on Jmax on input and yet still poorly resolved (say, only one point nonzero)
            ! then it can happen that the mask is gone after coarse_mesh.
            ! There is some overhead involved with keeping both structures the same, in the sense
            ! that MPI-communication is increased, if blocks on different CPU have to merged.
            call coarse_mesh( params, lgt_block, hvy_block, lgt_active, lgt_n, &
            lgt_sortednumlist, hvy_active, hvy_n, tree_ID, hvy_mask )
        else

            call coarse_mesh( params, lgt_block, hvy_block, lgt_active, lgt_n, &
            lgt_sortednumlist, hvy_active, hvy_n, tree_ID)
        endif
        call toc( "adapt_mesh (coarse_mesh)", MPI_Wtime()-t0 )


        ! update grid lists: active list, neighbor relations, etc
        t0 = MPI_Wtime()
        call update_grid_metadata(params, lgt_block, hvy_neighbor, lgt_active, lgt_n, &
        lgt_sortednumlist, hvy_active, hvy_n, tree_ID)
        call toc( "adapt_mesh (update neighbors)", MPI_Wtime()-t0 )

        iteration = iteration + 1

        ! see description above in argument list.
        if (present(external_loop)) then
            if (external_loop) exit ! exit loop
        endif
    end do

    ! The grid adaptation is done now, the blocks that can be coarsened are coarser.
    ! If a block is on Jmax now, we assign it the status +11.
    ! NOTE: Consider two blocks, a coarse on Jmax-1 and a fine on Jmax. If you refine only
    ! the coarse one (Jmax-1 -> Jmax), because you cannot refine the other one anymore
    ! (by defintion of Jmax), then the redundant layer in both blocks is different.
    ! To corrent that, you need to know which of the blocks results from interpolation and
    ! which one has previously been at Jmax. This latter one gets the 11 status.
    do k = 1, lgt_n
        lgt_id = lgt_active(k)
        if ( lgt_block( lgt_id, params%max_treelevel+ IDX_MESH_LVL) == params%max_treelevel ) then
            lgt_block( lgt_id, params%max_treelevel + IDX_REFINE_STS ) = 11
        end if
    end do

    !> At this point the coarsening is done. All blocks that can be coarsened are coarsened
    !! they may have passed several level also. Now, the distribution of blocks may no longer
    !! be balanced, so we have to balance load now
    t0 = MPI_Wtime()
    call balance_load( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, &
    lgt_n, lgt_sortednumlist, hvy_active, hvy_n, tree_ID )
    call toc( "adapt_mesh (balance_load)", MPI_Wtime()-t0 )


    call toc( "adapt_mesh (TOTAL)", MPI_wtime()-t1)
end subroutine adapt_mesh
