!> \brief check the gradedness after new refinement status
!
!> \details This routine is called after all blocks have been tagged whether to refine or coarsen or stay.
!! It now goes through the list of blocks and looks for refinement or coarsening states that would
!! result in an non-graded mesh. These mistakes are corrected, their status -1 or 0 is overwritten.
!! The status +1 is always conserved (recall to call respect_min_Jmax before).
!!
!! Since 04/2017, the new code checks all blocks that want to coarsen or remain, NOT the ones that
!! want to refine, as was done in prototypes. The reason is MPI: I cannot easily set the flags of
!! my neighbors, as they might reside on another proc.
! ********************************************************************************************

subroutine ensureGradedness_tree( params, tree_ID, mark_TMP_flag, check_daughters )

    implicit none
    type (type_params), intent(in)      :: params
    integer(kind=ik), intent(in)        :: tree_ID
    logical, intent(in), optional :: mark_TMP_flag          !< Set refinement of completeness to 0 or temporary flag
    logical, intent(in), optional :: check_daughters          !< For overfull CVS grids completeness has to check the mother as well

    integer(kind=ik)                    :: ierr, rank
    ! loop variables
    integer(kind=ik)                    :: k, i, N, mylevel, neighbor_level, &
    counter, hvy_id, ref_n, ref_me, Jmax, proc_id, lgt_id, REF_TMP_GRADED_STAY, lgt_id_sisters(2**params%dim)
    logical                             :: grid_changed, grid_changed_tmp, markTMPflag, checkDaughters    ! status of grid changing

    real(kind=rk) :: t0

    ! If we loop leaf-wise we do not want blocks to be set to 0 if not all sisters exist
    ! as we want to keep the wavelet-decomposed values and maybe just have to wait a step longer
    markTMPflag = .false.
    checkDaughters = .false.
    if (present(mark_TMP_flag)) markTMPflag = mark_TMP_flag
    if (present(check_daughters)) checkDaughters = check_daughters
    ! for security zone I have to wait for my finer neighbors to coarsen as their security zone might affect me
    ! as gradedness wants to keep me but completeness wants to coarsen we need another temporary flag 
    ! has to be > 0 to not have infinite loop with ensure_completeness
    REF_TMP_GRADED_STAY = 2

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas


    N = params%number_blocks
    Jmax = params%Jmax
    rank = params%rank

    ! we repeat the ensureGradedness_tree procedure until this flag is .false. since as long
    ! as the grid changes due to gradedness requirements, we have to check it again
    grid_changed = .true. ! set true to trigger the loop
    counter = 0


    do while ( grid_changed )
        ! we hope not to set the flag to .true. again in this iteration
        grid_changed = .false.

        t0 = MPI_wtime()
        ! we loop over heavy data here: parallel execution.
        ! We could make this collectively, only need to change finding the sisters for non-rank blocks without access to hvy_family
        ! However, then we would need to call find_sister over all lgt blocks with -1 (potentially many), so in the end it is sync_lgt against many times doesBlockExist
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k, tree_ID)
            call hvy2lgt(lgt_id, hvy_id, rank, N)
            ref_me = lgt_block( lgt_id , IDX_REFINE_STS )


            !-------------------------------------------------------------------
            ! completeness
            !-------------------------------------------------------------------
            ! We first remove the -1 flag from blocks which cannot be coarsened because their sisters
            ! disagree. If 1/4 or 1/8 blocks has 0 or +1 status, this cannot be changed. Therefore we first
            ! remove the status -1 from the blocks which have non-1 sisters. This is not only a question of
            ! simplicity. Consider 8 blocks on the same level:
            !      This block has to remove its -1 status as well, as the 4 neighbors to the right cannot coarsen
            !      v
            ! -1  -1    -1   0
            ! -1  -1    -1  -1
            ! It is thus clearly NOT enough to just look at the nearest neighbors in this ensureGradedness_tree routine.

            ! With security zone some blocks are waiting for their neighbors to coarsen, we check every first counter if blocks can now coarsen
            ! and they will be elsewise revoked, this is important if all 4/8 sisters share REF_TMP_TREATED_COARSEN
            if ( ref_me == -1 .or. (counter==0 .and. ref_me == REF_TMP_TREATED_COARSEN .and. markTMPflag)) then
                ! check if all sisters want to coarsen, remove the status it if they don't

                call ensure_completeness( params, lgt_id, hvy_family(hvy_ID, 2:1+2**params%dim), mark_TMP_flag=markTMPflag, check_daughters=checkDaughters )
                ! if the flag is removed, then it is removed only on mpiranks that hold at least
                ! one of the blocks, but the removal may have consequences everywhere. hence,
                ! we force the iteration to be executed one more time
                if (lgt_block(lgt_id , IDX_REFINE_STS) /= ref_me) then
                    ref_me = lgt_block(lgt_id , IDX_REFINE_STS)  ! updating variable in case it changed
                    grid_changed = .true.
                endif
            endif


            !-----------------------------------------------------------------------
            ! This block (still) wants to coarsen
            ! Neighbor blocks with REF_TMP_TREATED_COARSEN or REF_TMP_GRADED_STAY want to coarsen in theory but have to stay
            ! atleast for this iteration so they are treated as ref_n=0 and pass on REF_TMP_GRADED_STAY
            !-----------------------------------------------------------------------
            if ( ref_me == -1) then
                ! loop over all neighbors
                do i = 1, size(hvy_neighbor,2)
                    if ( hvy_neighbor( hvy_id, i ) < 0 ) cycle

                    ! check neighbor treelevel
                    mylevel         = lgt_block( lgt_id, IDX_MESH_LVL )
                    neighbor_level  = lgt_block( hvy_neighbor( hvy_id, i ), IDX_MESH_LVL )
                    ref_n = lgt_block( hvy_neighbor( hvy_id, i ), IDX_REFINE_STS )

                    if (mylevel == neighbor_level) then
                        ! neighbor on same level
                        ! block can not coarsen, if neighbor wants to refine
                        if ( ref_n == -1 ) then
                            ! neighbor wants to coarsen, as do I, we're on the same level -> ok
                        elseif ( ref_n == 0 .or. ref_n == REF_TMP_TREATED_COARSEN .or. ref_n == REF_TMP_GRADED_STAY ) then
                            ! neighbor wants to stay, I want to coarsen, we're on the same level -> ok
                        elseif ( ref_n == 1 ) then
                            ! neighbor wants to refine, I want to coarsen, we're on the same level -> NOT OK
                            ! I have at least to stay on my level.
                            lgt_block( lgt_id, IDX_REFINE_STS ) = 0
                            grid_changed = .true.

                        end if
                    elseif (mylevel - neighbor_level == 1) then
                        ! neighbor on lower level
                        if ( ref_n == -1 ) then
                            ! neighbor wants to coarsen, as do I, it is one level coarser, -> ok for me
                        elseif ( ref_n == 0 .or. ref_n == REF_TMP_TREATED_COARSEN .or. ref_n == REF_TMP_GRADED_STAY ) then
                            ! neighbor wants to stay, I want to coarsen, it is one level coarser, -> ok
                        elseif ( ref_n == 1 ) then
                            ! neighbor wants to refine, I want to coarsen, it is one level coarser, -> ok
                        end if
                    elseif (neighbor_level - mylevel == 1) then
                        ! neighbor on higher level
                        ! neighbor wants to refine, ...
                        if ( ref_n == +1) then
                            ! ... so I also have to refine (not only can I NOT coarsen, I actually have to refine!)
                            lgt_block( lgt_id, IDX_REFINE_STS ) = +1
                            grid_changed = .true.

                        elseif ( ref_n == 0 .or. ref_n == REF_TMP_TREATED_COARSEN .or. ref_n == REF_TMP_GRADED_STAY ) then
                            ! neighbor wants to stay and I want to coarsen, but
                            ! I cannot do that (there would be two levels between us)
                            if (ref_n == 0) then
                                lgt_block( lgt_id, IDX_REFINE_STS ) = 0
                            else  ! TMP neighbors pass on TMP flag
                                lgt_block( lgt_id, IDX_REFINE_STS ) = REF_TMP_GRADED_STAY
                            endif
                            grid_changed = .true.

                        elseif ( ref_n == -1 ) then
                            ! neighbor wants to coarsen, which is what I want too,
                            ! so we both would just go down one level together - that's fine
                            ! however, for lifted wavelets this is not fine as our neighbor might want to keep me to stay due to securityZone
                            if (params%useSecurityZone .and. params%isLiftedWavelet .and. markTMPflag) then
                                lgt_block( lgt_id, IDX_REFINE_STS ) = REF_TMP_GRADED_STAY
                                grid_changed = .true.
                            endif
                        end if
                    else
                        call abort(785879, "ERROR: ensureGradedness_tree: my neighbor does not seem to have -1,0,+1 level diff!")
                    end if
                end do ! loop over neighbors

            !-----------------------------------------------------------------------
            ! this block wants to stay on its level
            !-----------------------------------------------------------------------
            elseif (ref_me == 0 .or. ref_me == REF_TMP_GRADED_STAY .or. ref_me == REF_TMP_TREATED_COARSEN) then
                ! loop over all neighbors
                do i = 1, size(hvy_neighbor,2)
                    ! neighbor exists ? If not, we skip it
                    if ( hvy_neighbor( hvy_id, i ) < 0 ) cycle

                    mylevel     = lgt_block( lgt_id, IDX_MESH_LVL )
                    neighbor_level = lgt_block( hvy_neighbor( hvy_id, i ) , IDX_MESH_LVL )
                    ref_n = lgt_block( hvy_neighbor( hvy_id, i ) , IDX_REFINE_STS )

                    if (mylevel == neighbor_level) then
                        ! me and my neighbor are on the same level
                        ! As I'd wish to stay where I am, my neighbor is free to go -1,0,+1
                    elseif (mylevel - neighbor_level == 1) then
                        ! my neighbor is one level coarser
                        ! My neighbor can stay or refine, but not coarsen. This case is however handled above (coarsening inhibited)
                    elseif (neighbor_level - mylevel == 1) then
                        ! my neighbor is one level finer
                        if (ref_n == +1) then
                            ! neighbor refines (and we cannot inhibt that) so I HAVE TO do so as well
                            if (lgt_block( lgt_id, IDX_REFINE_STS )<+1) then
                                lgt_block( lgt_id, IDX_REFINE_STS ) = max( +1, lgt_block( lgt_id, IDX_REFINE_STS ) )
                                grid_changed = .true.
                            endif
                        end if
                    else
                        call abort(785879, "ERROR: ensureGradedness_tree: my neighbor does not seem to have -1,0,+1 level diff!")
                    end if
                end do
            end if ! refinement status

        end do ! loop over blocks
        call toc( "ensureGradedness_tree (processing part)", 130, MPI_Wtime()-t0 )

        ! since not all mpiranks change something in their light data, but all have to perform
        ! the same iterations, we sync the grid_changed indicator here. Note each mpirank changed
        ! only the blocks it holds, not blocks held by other ranks.
        grid_changed_tmp = grid_changed
        call MPI_Allreduce(grid_changed_tmp, grid_changed, 1, MPI_LOGICAL, MPI_LOR, WABBIT_COMM, ierr )

        !> after locally modifying refinement statusses, we need to synchronize light data
        t0 = MPI_wtime()
        call synchronize_lgt_data( params, refinement_status_only=.true. )
        call toc( "ensureGradedness_tree (sync_lgt)", 131, MPI_Wtime()-t0 )

        ! avoid infinite loops
        counter = counter + 1
        if (counter == 10*params%Jmax) call abort(785877, "ERROR: unable to build a graded mesh")

    end do ! end do of repeat procedure until grid_changed==.false.


    ! set REF_TMP_GRADED_STAY back to REF_TMP_TREATED_COARSEN so that they have a fitting flag < 0, collective operation
    do k = 1, lgt_n(tree_ID)
        lgt_id = lgt_active(k, tree_ID)
        if (lgt_block( lgt_id , IDX_REFINE_STS ) == REF_TMP_GRADED_STAY) then
            lgt_block( lgt_id , IDX_REFINE_STS ) = REF_TMP_TREATED_COARSEN
        endif
    enddo

end subroutine ensureGradedness_tree
