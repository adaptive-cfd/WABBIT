! ********************************************************************************************
!> \brief Apply mesh coarsening: Merge tagged blocks into new, coarser blocks
! ********************************************************************************************
subroutine executeCoarsening_tree( params, hvy_block, tree_ID )
    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    integer(kind=ik), intent(in)        :: tree_ID

    ! loop variables
    integer(kind=ik)                    :: k, Jmax, N, j
    ! list of block ids, proc ranks
    integer(kind=ik)                    :: light_ids(1:8), mpirank_owners(1:8)
    integer(kind=ik), allocatable, save :: xfer_list(:,:)
    ! rank of proc to keep the coarsened data
    integer(kind=ik)                    :: data_rank, n_xfer, ierr, lgtID

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas


    Jmax = params%Jmax
    ! number of blocks to merge, 4 or 8
    N = 2**params%dim
    ! at worst every block is on a different rank
    if (.not. allocated(xfer_list)) allocate(xfer_list(size(lgt_block,1),3))

    ! transfer counter
    n_xfer = 0


    !---------------------------------------------------------------------------
    ! first, prepare for xfer (gather information: which blocks are sent where)
    ! COLLECTIVE OPERATION
    !---------------------------------------------------------------------------
    do k = 1, lgt_n(tree_ID)
        ! Check if the block will be coarsened.
        !
        ! FIRST condition: only work on light data, if block is active. Usually, you would do that just with
        ! the active list. NOTE HERE: due to previous iterations, some light data id are already
        ! marked for xfer, (and thus given the -7 status) but they are still in active block list
        ! so: check again if this block is REALLY active (the active list is updated later)
        !
        ! SECOND condition: block wants to coarsen, i.e. it has the status -1. Note the routine
        ! ensureGradedness_tree removes the -1 flag if not all sister blocks share it
        lgtID = lgt_active(k, tree_ID)

        if ( get_tc(lgt_block(lgtID, Jmax+IDX_TC_1 : Jmax+IDX_TC_2)) >= 0 .and. lgt_block(lgtID, Jmax+IDX_REFINE_STS) == -1) then
            ! find all sisters (including the block in question, so four or eight blocks)
            ! their light IDs are in "light_ids" and ordered by their last treecode-digit
            call findSisters_tree( params, lgtID, light_ids(1:N), tree_ID )

            ! figure out on which rank the sisters lie, xfer them if necessary
            do j = 1, N
                call lgt2proc( mpirank_owners(j), light_ids(j), params%number_blocks )
            enddo

            ! The merging will be done on the mpirank which holds the most of the sister blocks
            data_rank = most_common_element( mpirank_owners(1:N) )

            do j = 1, N
                if (mpirank_owners(j) /= data_rank) then
                    ! MPI xfer required. Add the xfer to the list
                    n_xfer = n_xfer + 1
                    xfer_list(n_xfer, 1) = mpirank_owners(j)  ! send from this rank ..
                    xfer_list(n_xfer, 2) = data_rank          ! ... to this rank
                    xfer_list(n_xfer, 3) = light_ids(j)       ! transfer this block
                endif

                ! don't forget: mark all 4/8 sisters as treated here, in order not to trigger this
                ! loop again: we use the temporary status -7
                lgt_block(light_ids(j), Jmax+IDX_REFINE_STS) = -7
            enddo
        endif
    enddo

    ! actual xfer
    call block_xfer( params, xfer_list, n_xfer, hvy_block )

    ! the active lists are outdates after the transfer: we need to create
    ! them or findSisters_tree will not be able to do its job
    call createActiveSortedLists_tree( params, tree_ID )

    ! actual merging
    do k = 1, lgt_n(tree_ID)
        ! FIRST condition: only work on light data, if block is active. Usually, you would do that just with
        ! the active list. NOTE HERE: due to previous loops, some light data id are already
        ! coarsened, (and thus given the -1) but they are still in active block list
        ! so: check again if this block is REALLY active (the active list is updated later)
        !
        ! SECOND condition: block wants to coarsen, i.e. it has the status -1. Note the routine
        ! ensureGradedness_tree removes the -1 flag and sets -7 (temporarily assigned above) flag if not all sister blocks share it
        lgtID = lgt_active(k, tree_ID)
        if ( get_tc(lgt_block(lgtID, Jmax+IDX_TC_1 : Jmax+IDX_TC_2)) >= 0 .and. lgt_block(lgtID, Jmax+IDX_REFINE_STS) == -7) then
            ! merge the four blocks into one new block. Merging is done in two steps,
            ! first for light data (which all CPUS do redundantly, so light data is kept synched)
            ! Then only the responsible rank will perform the heavy data merging.
            call findSisters_tree( params, lgtID, light_ids(1:N), tree_ID )
            ! note the newly merged block has status 0
            call merge_blocks( params, hvy_block, light_ids(1:N) )
        endif
    enddo
end subroutine executeCoarsening_tree
