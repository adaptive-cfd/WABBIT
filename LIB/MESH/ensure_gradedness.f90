!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name ensure_gradedness.f90
!> \version 0.5
!> \author msr, engels
!
!> \brief check the gradedness after new refinement status
!
!> \details This routine is called after all blocks have been tagged whether to refine or coarsen or stay.
!! It now goes through the list of blocks and looks for refinement or coarsening states that would
!! result in an non-graded mesh. These mistakes are corrected, their status -1 or 0 is overwritten.
!! The status +1 is always conserved (recall to call respect_min_Jmax before).
!!
!!
!! Since 04/2017, the new code checks all blocks that want to coarsen or remain, NOT the ones that
!! want to refine, as was done in prototypes. The reason is MPI: I cannot easily set the flags of
!! my neighbors, as they might reside on another proc.
!!
!!
!! = log ======================================================================================
!! \n
!! 10/11/16 - switch to v0.4 \n
!! 23/11/16 - rework complete subroutine: use list of active blocks, procs works now on light data \n
!! 03/02/17 - insert neighbor_num variable to use subroutine for 2D and 3D data \n
!! 05/04/17 - Improvement: Ensure a graded mesh in any case, not only in the coarsen states (which was done before)
!! 09/06/17 - speed up with switching to 8bit integer for working arrays, shorten send/receive buffer
!
! ********************************************************************************************


subroutine ensure_gradedness( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, &
    lgt_sortednumlist, hvy_active, hvy_n )

    !---------------------------------------------------------------------------------------------
    ! modules

    !---------------------------------------------------------------------------------------------
    ! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)
    !> neighbor list
    integer(kind=ik), intent(inout)     :: hvy_neighbor(:, :)
    !> active_block_list (light data)
    integer(kind=ik), intent(inout)     :: lgt_active(:)
    integer(kind=tsize), intent(inout)  :: lgt_sortednumlist(:,:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_n
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(inout)     :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(inout)     :: hvy_n

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank
    ! loop variables
    integer(kind=ik)                    :: k, i, N, mylevel, neighbor_level, &
    counter, hvy_id, neighbor_status, Jmax, proc_id, lgt_id
    ! status of grid changing
    logical                             :: grid_changed, test2

    ! it turned out that in some situations ensure_completeness is very expensive, mostly because
    ! of the find_sisters routine. We therefore at least do this procedure only once and not
    ! in each iteration of the algorithm.
    integer(kind=ik), allocatable, save  :: sisters(:)

    ! number of neighbor relations
    ! 2D: 16, 3D: 74
    integer(kind=ik) :: neighbor_num
    real(kind=rk) :: t0


    N = params%number_blocks
    Jmax = params%max_treelevel

    ! set MPI parameter
    rank = params%rank

    if ( params%threeD_case ) then
        ! 3D:
        neighbor_num = 74
        if (.not.allocated(sisters)) allocate( sisters(8) )
    else
        ! 2D:
        neighbor_num = 16
        if (.not.allocated(sisters)) allocate( sisters(4) )
    end if
    sisters = -1

    ! NOTE: The status +11 is required for ghost nodes synching, if the maxlevel
    ! is reached. It means just the same as 0 (stay) in the context of this routine


    ! we repeat the ensure_gradedness procedure until this flag is .false. since as long
    ! as the grid changes due to gradedness requirements, we have to check it again
    grid_changed = .true. ! set true to trigger the loop
    counter = 0


    do while ( grid_changed )
        ! we hope not to set the flag to .true. again in this iteration
        grid_changed = .false.

        t0 = MPI_wtime()
        ! we loop over heavy data here: parallel execution.
        do k = 1, hvy_n
            hvy_id = hvy_active(k)
            call hvy_id_to_lgt_id(lgt_id, hvy_id, rank, N)

            !-------------------------------------------------------------------
            ! completeness
            !-------------------------------------------------------------------
            ! We first remove the -1 flag from blocks which cannot be coarsened because their sisters
            ! disagree. If 1 of 4 or 1/8 blocks has 0 or +1 status, this cannot be changed. Therefore we first
            ! remove the status -1 from the blocks which have non-1 sisters. This is not only a question of
            ! simplicity. Consider 8 blocks on the same level:
            !      This block has to remove its -1 status as well, as the 4 neighbors to the right cannot coarsen
            !      v
            ! -1  -1    -1   0
            ! -1  -1    -1  -1
            ! It is thus clearly NOT enough to just look at the nearest neighbors in this ensure_gradedness routine.
            if ( lgt_block( lgt_id , Jmax+2 ) == -1) then
                ! find the sisters of this block
                call find_sisters(params, lgt_id, sisters, lgt_block, lgt_n, lgt_sortednumlist)
                ! check if all sisters share the -1 status, remove it if they don't
                call ensure_completeness( params, lgt_block, lgt_id, sisters )
                ! if the flag is removed, then it is removed only on mpiranks that hold at least
                ! one of the blocks, but the removal may have consequences everywhere. hence,
                ! we force the iteration to be executed one more time
                if (lgt_block(lgt_id , Jmax+2) /= -1)  grid_changed = .true.
            endif

            !-----------------------------------------------------------------------
            ! This block (still) wants to coarsen
            !-----------------------------------------------------------------------
            if ( lgt_block( lgt_id , Jmax+2 ) == -1) then
                ! loop over all neighbors
                do i = 1, neighbor_num
                    if ( hvy_neighbor( hvy_id, i ) /= -1 ) then

                        ! check neighbor treelevel
                        mylevel         = lgt_block( lgt_id, Jmax+1 )
                        neighbor_level  = lgt_block( hvy_neighbor( hvy_id, i ) , Jmax+1 )
                        neighbor_status = lgt_block( hvy_neighbor( hvy_id, i ) , Jmax+2 )

                        if (mylevel == neighbor_level) then
                            ! neighbor on same level
                            ! block can not coarsen, if neighbor wants to refine
                            if ( neighbor_status == -1 ) then
                                ! neighbor wants to coarsen, as do I, we're on the same level -> ok
                            elseif ( neighbor_status == 0 .or. neighbor_status == 11 ) then
                                ! neighbor wants to stay, I want to coarsen, we're on the same level -> ok
                            elseif ( neighbor_status == 1 ) then
                                ! neighbor wants to refine, I want to coarsen, we're on the same level -> NOT OK
                                ! I have at least to stay on my level.
                                ! Note we cannot simply set 0 as we could accidentally overwrite a refinement flag
                                if (lgt_block( lgt_id, Jmax+2 )<0) then
                                    lgt_block( lgt_id, Jmax+2 ) = max( 0, lgt_block( lgt_id, Jmax+2 ) )
                                    grid_changed = .true.
                                endif

                            end if
                        elseif (mylevel - neighbor_level == 1) then
                            ! neighbor on lower level
                            if ( neighbor_status == -1 ) then
                                ! neighbor wants to coarsen, as do I, he is one level coarser, -> ok
                            elseif ( neighbor_status == 0 .or. neighbor_status == 11 ) then
                                ! neighbor wants to stay, I want to coarsen, he is one level coarser, -> ok
                            elseif ( neighbor_status == 1 ) then
                                ! neighbor wants to refine, I want to coarsen,  he is one level coarser, -> ok
                            end if
                        elseif (neighbor_level - mylevel == 1) then
                            ! neighbor on higher level
                            ! neighbor wants to refine, ...
                            if ( neighbor_status == +1) then
                                ! ... so I also have to refine (not only can I NOT coarsen, I actually
                                ! have to refine!)
                                if (lgt_block( lgt_id, Jmax+2 )<+1) then
                                    lgt_block( lgt_id, Jmax+2 ) = max( +1, lgt_block( lgt_id, Jmax+2 ) )
                                    grid_changed = .true.
                                endif

                            elseif ( neighbor_status == 0 .or. neighbor_status == 11) then
                                ! neighbor wants to stay and I want to coarsen, but
                                ! I cannot do that (there would be two levels between us)
                                ! Note we cannot simply set 0 as we could accidentally overwrite a refinement flag
                                if (lgt_block( lgt_id, Jmax+2 )<0) then
                                    lgt_block( lgt_id, Jmax+2 ) = max( 0, lgt_block( lgt_id, Jmax+2 ) )
                                    grid_changed = .true.
                                endif

                            elseif ( neighbor_status == -1) then
                                ! neighbor wants to coarsen, which is what I want too,
                                ! so we both would just go up one level together - that's fine
                            end if
                        else
                            call abort(785879, "ERROR: ensure_gradedness: my neighbor does not seem to have -1,0,+1 level diff!")
                        end if
                    end if ! if neighbor exists
                end do ! loop over neighbors

        !-----------------------------------------------------------------------
        ! this block wants to stay on its level
        !-----------------------------------------------------------------------
        elseif (lgt_block( lgt_id , Jmax+2 ) == 0 .or. lgt_block( lgt_id , Jmax+2 ) == 11 ) then
                ! loop over all neighbors
                do i = 1, neighbor_num
                    ! neighbor exists ? If not, this is a bad error
                    if ( hvy_neighbor( hvy_id, i ) /= -1 ) then
                        mylevel     = lgt_block( lgt_id, Jmax+1 )
                        neighbor_level = lgt_block( hvy_neighbor( hvy_id, i ) , Jmax+1 )
                        neighbor_status = lgt_block( hvy_neighbor( hvy_id, i ) , Jmax+2 )

                        if (mylevel == neighbor_level) then
                            ! me and my neighbor are on the same level
                            ! As I'd wish to stay where I am, my neighbor is free to go -1,0,+1
                        elseif (mylevel - neighbor_level == 1) then
                            ! my neighbor is one level coarser
                            ! My neighbor can stay or refine, but not coarsen. This case is however handled above (coarsening inhibited)
                        elseif (neighbor_level - mylevel == 1) then
                            ! my neighbor is one level finer
                            if (neighbor_status == +1) then
                                ! neighbor refines (and we cannot inhibt that) so I HAVE TO do so as well
                                if (lgt_block( lgt_id, Jmax+2 )<+1) then
                                    lgt_block( lgt_id, Jmax+2 ) = max( +1, lgt_block( lgt_id, Jmax+2 ) )
                                    grid_changed = .true.
                                endif
                            end if
                        else
                            call abort(785879, "ERROR: ensure_gradedness: my neighbor does not seem to have -1,0,+1 level diff!")
                        end if
                    end if ! if neighbor exists
                end do
            end if ! refinement status
        end do ! loop over blocks
        call toc( params, "ensure_gradedness (processing part)", MPI_Wtime()-t0 )

        ! since not all mpiranks change something in their light data, but all have to perform
        ! the same iterations, we sync the grid_changed indicator here. Note each mpirank changed
        ! only the blocks it holds, not blocks held by other ranks.
        test2 = grid_changed
        call MPI_Allreduce(test2, grid_changed, 1, MPI_LOGICAL, MPI_LOR, WABBIT_COMM, ierr )


        !> after locally modifying refinement statusses, we need to synchronize light data
        t0 = MPI_wtime()
        call synchronize_lgt_data( params, lgt_block, refinement_status_only=.true. )
        call toc( params, "ensure_gradedness (sync_lgt)", MPI_Wtime()-t0 )

        ! avoid infinite loops
        counter = counter + 1
        if (counter == 10) call abort(785877, "ERROR: unable to build a graded mesh")

    end do ! end do of repeat procedure until grid_changed==.false.

end subroutine ensure_gradedness
