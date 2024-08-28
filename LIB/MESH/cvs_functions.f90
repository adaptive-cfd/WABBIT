subroutine cvs_decompose_tree(params, hvy_block, tree_ID, hvy_tmp, log_blocks, log_iterations, verbose_check)
    implicit none

    type (type_params), intent(in)      :: params
    !> heavy data array
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    integer(kind=ik), intent(in)        :: tree_ID
    !> heavy tmp data array - block data.
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    !> some information that we can log so that we know what happened in the loops
    integer, intent(out), optional       :: log_blocks(:), log_iterations
    logical, intent(in), optional       :: verbose_check  !< No matter the value, if this is present we debug

    integer(kind=ik)     :: k, lgt_ID, hvy_ID, lgt_n_old, g_this, level_me, ref_stat, iteration, level
    real(kind=rk)        :: t_block, t_loop
    character(len=clong) :: toc_statement
    logical              :: iterate

    ! ! at first, initialize all mothers without any values yet
    ! t_block = MPI_Wtime()
    ! call cvs_init_mothers_tree(params, tree_ID)
    ! call toc( "adapt_tree (init mothers)", 103, MPI_Wtime()-t_block )
    ! ! Leaf layer starts and gets temporary flag
    ! do k = 1, hvy_n(tree_ID)
    !     hvy_ID = hvy_active(k, tree_ID)
    !     ! check if block is leaf-block / any daughter exists
    !     if (any(hvy_family(hvy_ID, 2+2**params%dim:1+2**(params%dim+1)) /= -1)) then
    !         call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)
    !         lgt_block(lgt_ID, IDX_REFINE_STS) = REF_TMP_UNTREATED
    !     endif
    ! end do
    ! ! sync ref status


    ! In order to not investigate blocks again, we will give everyone a temporary flag that they are not wavelet decomposed yet
    do k = 1, lgt_n(tree_ID)
        lgt_ID = lgt_active(k, tree_ID)
        lgt_block(lgt_ID, IDX_REFINE_STS) = REF_TMP_UNTREATED            
    end do

    ! Idea: One of the very expensive things is updateMetadata_tree for overfull grids, where a block checks all 56*3 possible neighbors
    ! We could create beforehand all mothers for this loop and give them a REF flag of TMP_MOTHER or something so that we do not sync with them
    ! Then, we always only move the leaf layer downwards, WD the blocks and sync values with the mothers
    ! Syncing lgt_block is then also only done for refinement status
    ! JB ToDo, give new refinement status of blocks that have none-sense values and do not allow syncing

    ! iteration for wavelet decomposition - no blocks will be deleted during the loop!
    level       = maxActiveLevel_tree(tree_ID) ! leaf-wise blocks can have this as maximum level
    iterate     = .true.
    iteration   = 0
    do while (iterate)
        t_loop = MPI_Wtime()
        lgt_n_old = lgt_n(tree_ID)

        ! ! I often need to check block values when debugging so I leave it here
        ! do k = 1, hvy_n(tree_ID)
        !     hvy_ID = hvy_active(k, tree_ID)
        !     call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)
        !     if (params%rank == 0 .and. lgt_block(lgt_id, IDX_REFINE_STS) == REF_TMP_UNTREATED .and. iteration >= 2) then
        !     ! if (lgt_block(lgt_id, IDX_MESH_LVL) == 3 .and. lgt_block(lgt_id, IDX_TC_2) == 54525952) then
        !         write(*, '("1 I-", i0, " R-", i0, " B-", i0, " L-", i0, " Ref-", i0, " TC-", i0)') iteration, params%rank, lgt_ID, &
        !             lgt_block(lgt_id, IDX_MESH_LVL), lgt_block(lgt_id, IDX_REFINE_STS), lgt_block(lgt_id, IDX_TC_2)
        !         call write_neighborhood_info(hvy_neighbor(hvy_id, :), params%dim)
        !         ! call dump_block_fancy(hvy_block(:, :, :, 1:1, hvy_id), "block_cvs_9.txt", params%Bs, params%g)
        !     endif
        ! enddo

        ! synchronize ghost nodes - required to apply wavelet filters
        ! block only needs information from medium and fine neighbors as CE will cut dependency to coarse neighbors
        t_block = MPI_Wtime()
        g_this = max(ubound(params%HD,1),ubound(params%GD,1))
        ! for coarse extension we are not dependend on coarser neighbors so lets skip the syncing
        if (params%isLiftedWavelet) then
            call sync_TMP_from_MF( params, hvy_block, tree_ID, REF_TMP_UNTREATED, g_minus=g_this, g_plus=g_this, hvy_tmp=hvy_tmp)
            call toc( "adapt_tree (sync lvl <- MF)", 103, MPI_Wtime()-t_block )

            write(toc_statement, '(A, i0, A)') "adapt_tree (WD it ", iteration, " sync lvl <- MF)"
            call toc( toc_statement, 1100+iteration, MPI_Wtime()-t_block )
        ! unlifted wavelets need coarser neighbor values for their WC so we need to sync them too
        else
            call sync_TMP_from_all( params, hvy_block, tree_ID, REF_TMP_UNTREATED, g_minus=g_this, g_plus=g_this, hvy_tmp=hvy_tmp)
            call toc( "adapt_tree (sync lvl <- all)", 103, MPI_Wtime()-t_block )

            write(toc_statement, '(A, i0, A)') "adapt_tree (WD it ", iteration, " sync lvl <- all)"
            call toc( toc_statement, 1100+iteration, MPI_Wtime()-t_block )
        endif

        ! Wavelet-transform all remaining non-decomposed blocks
        ! From now on until wavelet retransform hvy_block will hold the wavelet decomposed values in spaghetti form
        t_block = MPI_Wtime()
        do k = 1, hvy_n(tree_ID)
            hvy_ID = hvy_active(k, tree_ID)
    
            ! We compute detail coefficients on the fly here, for all blocks
            ! on the level.
            call hvy2lgt( lgt_ID, hvy_ID, params%rank, params%number_blocks )
            ref_stat = lgt_block( lgt_ID, IDX_REFINE_STS )
    
            ! FWT required for a block that is on the level
            if (ref_stat == REF_TMP_UNTREATED) then
                ! hvy_tmp now is a copy with sync'ed ghost points.
                hvy_tmp(:,:,:,1:size(hvy_block, 4),hvy_ID) = hvy_block(:,:,:,1:size(hvy_block, 4),hvy_ID)
                level_me = lgt_block( lgt_ID, IDX_MESH_LVL )

                ! Compute wavelet decomposition
                ! For Jmax if dealiasing, just compute H filter
                ! Data SC/WC now in Spaghetti order
                if (level_me == params%Jmax .and. params%force_maxlevel_dealiasing) then
                    call blockFilterXYZ_vct( params, hvy_tmp(:,:,:,1:size(hvy_block, 4),hvy_ID), hvy_block(:,:,:,1:size(hvy_block, 4),hvy_ID), params%HD, &
                        lbound(params%HD, dim=1), ubound(params%HD, dim=1), do_restriction=.true.)
                else
                    call waveletDecomposition_block(params, hvy_block(:,:,:,:,hvy_ID))
                endif
            endif
        end do
        call toc( "adapt_tree (FWT)", 104, MPI_Wtime()-t_block )

        write(toc_statement, '(A, i0, A)') "adapt_tree (WD it ", iteration, " FWT)"
        call toc( toc_statement, 1200+iteration, MPI_Wtime()-t_block )

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! coarseExtension: remove wavelet coefficients near a fine/coarse interface
        ! on the fine block. Does nothing in the case of CDF60, CDF40 or CDF20.
        ! This is the first coarse extension before removing blocks
        ! As no blocks are deleted, working on newly created blocks is suficient, they are non-leafs anyways and will be skipped
        t_block = MPI_Wtime()
        call coarse_extension_modify(params, hvy_block, hvy_tmp, tree_ID, CE_case="ref", s_val=REF_TMP_UNTREATED)
        call toc( "adapt_tree (coarse_extension_modify)", 105, MPI_Wtime()-t_block )

        ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Mother blocks: we now create mother blocks and fill their values with the SC from their daughters
        ! We do this for every but the last iteration, where no mother blocks can be created anyways
        if (level > params%Jmin) then
            ! Prepare refinement status for updating the mother blocks
            do k = 1, lgt_n(tree_ID)
                lgt_ID = lgt_active(k, tree_ID)
                if (any(lgt_block(lgt_ID, IDX_REFINE_STS) == (/ REF_TMP_UNTREATED , REF_TMP_TREATED_COARSEN, REF_TMP_GRADED_STAY /))) then
                    lgt_block(lgt_ID, IDX_REFINE_STS) = -1
                endif
            end do

            ! respect lower level boundary
            t_block = MPI_Wtime()
            call respectJmaxJmin_tree( params, tree_ID )
            call toc( "adapt_tree (respectJmaxJmin_tree)", 108, MPI_Wtime()-t_block )

            ! when blocks don't have all neighbors available on medium lvl, the newly created mother block will not be able
            ! to find all needed neighbors, so we let the block wait if that is the case
            do k = 1, hvy_n(tree_ID)
                hvy_id = hvy_active(k, tree_ID)
                call ensure_mother_neighbors(params, hvy_id, mark_TMP_flag=.true.)
            enddo
            ! after locally modifying refinement statusses, we need to synchronize light data
            call synchronize_lgt_data( params, refinement_status_only=.true. )
            ! watch for completeness
            do k = 1, hvy_n(tree_ID)
                hvy_id = hvy_active(k, tree_ID)
                call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
                call ensure_completeness( params, lgt_id, hvy_family(hvy_ID, 2:1+2**params%dim), mark_TMP_flag=.true. )
            enddo
            ! Details are not synced between procs, however all procs correctly know which blocks to delete if they have one block in there
            ! The rest will be synced in executeCoarsening
            ! call synchronize_lgt_data( params, refinement_status_only=.true. )

            ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! Create or update mother blocks with SC of daughter blocks
            ! This uses the already wavelet decomposed blocks and copies their SC into the new mother block
            t_block = MPI_Wtime()
            call executeCoarsening_WD_tree( params, hvy_block, tree_ID, mark_TMP_flag=.true., no_deletion=.true.)
            call toc( "adapt_tree (update_mothers)", 111, MPI_Wtime()-t_block )

            write(toc_statement, '(A, i0, A)') "adapt_tree (WD it ", iteration, " update_mothers)"
            call toc( toc_statement, 1650+iteration, MPI_Wtime()-t_block )

            ! In executeCoarsening we already sync list and family but not neighbors, as no blocks are deleted there we only do this then here
            ! This call is surprisingly expensive as for overlapping search we always check all 56*3 possible neighbors for every loop
            t_block = MPI_Wtime()
            call updateMetadata_tree(params, tree_ID, search_overlapping=.true., update_lists=.false., update_neighbors=.true., update_family=.false.)
            call toc( "adapt_tree (updateMetadata_tree)", 101, MPI_Wtime()-t_block )

        endif

        ! iteration counter
        iteration = iteration + 1
        level = level - 1
        ! loop continues until we are on the lowest level.
        iterate = (level >= params%Jmin)

        ! log some statistics - amount of blocks, after iteration it shoul be increased, the difference showing how many blocks we treated
        if (present(log_blocks)) log_blocks(iteration) = lgt_n(tree_ID)
        write(toc_statement, '(A, i0, A)') "adapt_tree (WD it ", iteration, " TOTAL)"
        call toc( toc_statement, 1700+iteration, MPI_Wtime()-t_loop )

        ! write(toc_statement, '(A, i6.6, A)') "test_", iteration, "000000.h5"
        ! call saveHDF5_tree( toc_statement, 0.0 + dble(iteration), iteration, 1, params, hvy_tmp, tree_ID, no_sync=.true.)

        if (present(verbose_check)) then
            write(*, '(A, i0, A, i0, A)') "Loop ", iteration, " with ", lgt_n(tree_ID), " blocks"
        endif
    end do

    ! log some statistics - amount of iterations
    if (present(log_iterations)) log_iterations = iteration
end subroutine


subroutine cvs_reconstruct_tree(params, hvy_block, hvy_tmp, tree_ID)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params
    
    implicit none

    type (type_params), intent(in)      :: params
    !> heavy data array
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy tmp data array - block data.
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik), intent(in)        :: tree_ID

    ! loop variables
    integer(kind=ik)                    :: iteration, k, Jmax_active, Jmin_active, level, hvy_ID, lgt_ID, level_me, Jmin, g_spaghetti
    real(kind=rk)                       :: t_block, t_all, t_loop
    character(len=clong)                :: toc_statement


    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas

    Jmin        = params%Jmin
    if (Jmin<1) call abort(2202243, "Currently, setting Jmin<1 is not possible")

    ! it turns out, when the coefficients are spaghetti-ordered,
    ! we can sync only even numbers of points and save one for odd numbered
    g_spaghetti = params%g/2*2

    Jmax_active = maxActiveLevel_tree(tree_ID)
    iteration = 0
    do level = Jmin, Jmax_active
        t_loop = MPI_Wtime()

        ! synch SC and WC from coarser neighbours and same-level neighbours in order to apply the correct wavelet reconstruction
        ! Attention: This uses hvy_temp for coarser neighbors to predict the data, as we want the correct SC from coarser neighbors
        t_block = MPI_Wtime()
        call sync_SCWC_from_MC( params, hvy_block, tree_ID, hvy_tmp, g_minus=g_spaghetti, g_plus=g_spaghetti, level=level)
        call toc( "adapt_tree (sync all <- MC)", 112, MPI_Wtime()-t_block )

        write(toc_statement, '(A, i0, A)') "adapt_tree (WR it ", iteration, " sync level <- MC)"
        call toc( toc_statement, 1750+iteration, MPI_Wtime()-t_block )

        ! After synching with coarser neighbors, the spaghetti form of the whole ghost patch is filled with values
        ! We now need to delete all values in spots of WC, leaving only the correct SC values from coarse neighbours
        ! This uses the fact that values refined on the coarse grid points with wc=0 are the copied SC values
        ! We do this in order to bundle up the synchronizations in the step before as the modification is really cheap
        t_block = MPI_Wtime()
        call coarse_extension_modify(params, hvy_block, hvy_tmp, tree_ID, CE_case="level", s_val=level, clear_wc_only=.true.)
        call toc( "adapt_tree (coarse_extension_modify)", 105, MPI_Wtime()-t_block )

        ! Wavelet-reconstruct blocks on level
        t_block = MPI_Wtime()
        do k = 1, hvy_n(tree_ID)
            hvy_ID = hvy_active(k, tree_ID)
            call hvy2lgt( lgt_ID, hvy_ID, params%rank, params%number_blocks )
            level_me = lgt_block( lgt_ID, IDX_MESH_LVL )

            ! RWT required for a block that is on the level
            if (level_me == level) then
                ! Compute wavelet reconstruction
                call waveletReconstruction_block(params, hvy_block(:,:,:,:,hvy_ID))

                ! make copy of data in hvy_tmp as we sync from coarser neighbors using that and I do not want to write another logic currently
                hvy_tmp(:,:,:,1:size(hvy_block, 4),hvy_ID) = hvy_block(:,:,:,1:size(hvy_block, 4),hvy_ID)
            endif
        end do
        call toc( "adapt_tree (reset or RWT)", 113, MPI_Wtime()-t_block )

        write(toc_statement, '(A, i0, A)') "adapt_tree (WR it ", iteration, " RWT)"
        call toc( toc_statement, 1850+iteration, MPI_Wtime()-t_block )

        ! now we need to update the SC of all daughters
        t_block = MPI_Wtime()
        call sync_M2D_level(params, hvy_block, tree_ID, level)
        call toc( "adapt_tree (sync level 2 daughters)", 114, MPI_Wtime()-t_block )

        write(toc_statement, '(A, i0, A)') "adapt_tree (WR it ", iteration, " sync level 2 daughters)"
        call toc( toc_statement, 1900+iteration, MPI_Wtime()-t_block )

        write(toc_statement, '(A, i0, A)') "adapt_tree (WR it ", iteration, " TOTAL)"
        call toc( toc_statement, 1950+iteration, MPI_Wtime()-t_loop )

        iteration = iteration + 1
    enddo
end subroutine



! subroutine cvs_init_mothers_tree(params, tree_ID, verbose_check)
!     implicit none

!     type (type_params), intent(in)      :: params
!     integer(kind=ik), intent(in)        :: tree_ID
!     logical, intent(in), optional       :: verbose_check  !< No matter the value, if this is present we debug

!     integer(kind=ik)     :: i_loop, j, k, N, lgt_ID, hvy_ID, lgt_n_old, g_this, level_b, ref_stat, iteration, level, data_rank, lgt_merge_id
!     real(kind=rk)        :: t_block, t_loop
!     character(len=clong) :: toc_statement
!     logical              :: iterate
!     ! list of block ids, proc ranks
!     integer(kind=ik)                    :: lgt_daughters(1:8), rank_daughters(1:8)
!     integer(kind=tsize)                 :: treecode

!     ! Start with all blocks
!     do k = 1, lgt_n(tree_ID)
!         lgt_ID = lgt_active(k, tree_ID)
!         lgt_block(lgt_ID, IDX_REFINE_STS) = REF_TMP_UNTREATED            
!     end do

!     ! loop downwards until all blocks are on Jmin
!     level       = maxActiveLevel_tree(tree_ID) ! leaf-wise blocks can have this as maximum level
!     iteration   = 0
!     do i_loop = level, params%Jmin, -1
!         ! watch for completeness
!         do k = 1, hvy_n(tree_ID)
!             hvy_id = hvy_active(k, tree_ID)
!             call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
!             call ensure_completeness( params, lgt_id, hvy_family(hvy_ID, 2:1+2**params%dim), mark_TMP_flag=.true. )
!         enddo

!         !---------------------------------------------------------------------------
!         ! create new empty blocks on rank with most daughters if it does not exist
!         !---------------------------------------------------------------------------
!         do k = 1, hvy_n(tree_ID)
!             hvy_ID = hvy_active(k, tree_ID)
!             call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)
!             level_b = lgt_block( lgt_ID, IDX_MESH_LVL )

!             ! check if block wants to be refined and does not have a mother yet
!             if ( lgt_block(lgt_ID, IDX_REFINE_STS) == -1 .and. hvy_family(hvy_ID, 1) == -1) then
!                 ! Get all sisters and figure out on which rank the sisters lie
!                 lgt_daughters(1:N) = hvy_family(hvy_ID, 2:1+2**params%dim)
!                 do j = 1, N
!                     call lgt2proc( rank_daughters(j), lgt_daughters(j), params%number_blocks )
!                 enddo
!                 ! The merging will be done on the mpirank which holds the most of the sister blocks
!                 data_rank = most_common_element( rank_daughters(1:N) )
    
!                 ! construct new mother block if on my rank and create light data entry for the new block
!                 if (data_rank == params%rank) then
!                     call get_free_local_light_id(params, data_rank, lgt_merge_id, message="init_mothers")
!                     treecode = get_tc(lgt_block( lgt_daughters(1), IDX_TC_1:IDX_TC_2 ))

!                     ! set meta_data of mother block
!                     lgt_block( lgt_merge_id, : ) = -1
!                     call set_tc(lgt_block( lgt_merge_id, IDX_TC_1:IDX_TC_2), tc_clear_until_level_b(treecode, &
!                         dim=params%dim, level=level_b-1, max_level=params%Jmax))
!                     lgt_block( lgt_merge_id, IDX_MESH_LVL ) = level_b-1
!                     if (level_b-1 > params%Jmin) then ! blocks >Jmin want to coarsen in next iteration
!                         lgt_block( lgt_merge_id, IDX_REFINE_STS ) = -1
!                     else ! blocks on JMin cannot coarsen further
!                         lgt_block( lgt_merge_id, IDX_REFINE_STS ) = 0
!                     endif
!                     lgt_block( lgt_merge_id, IDX_TREE_ID ) = tree_ID

!                     ! update sisters on my rank that they have found their mother and can be skipped
!                     do j = 1, N
!                         if (rank_daughters(j) == params%rank .and. lgt_daughters(j) /= lgt_merge_id) then
!                             call lgt2hvy(hvy_ID, lgt_daughters(j), params%rank, params%number_blocks)
!                             hvy_family(hvy_ID, 1) = lgt_merge_id
!                         endif
!                     enddo

!                     ! write(*, '("Rank ", i0, " created a new block: ", i0, " with TC ", b64.64)') rank, lgt_merge_id, tc_clear_until_level_b(treecode, &
!                     ! dim=params%dim, level=level-1, max_level=params%Jmax)
!                 endif

!                 ! change refinement block that it was treated, we will not need it again
!                 lgt_block( lgt_id, IDX_REFINE_STS ) = 0
!             endif
!         enddo
!         ! the active lists are outdated, so lets resynch
!         call synchronize_lgt_data( params, refinement_status_only=.false.)
!         ! update metadata for overlapping grid, blocks need to find their sisters and mothers!
!         call updateMetadata_tree(params, tree_ID, search_overlapping=.true., update_neighbors=.false.)

!         ! increase iteration counter
!         iteration = iteration+1
!     enddo

!     ! At last, we update all neighbor relations
!     call updateMetadata_tree(params, tree_ID, search_overlapping=.true., update_lists=.false., update_neighbors=.true., update_family=.false.)
! end subroutine



! JB ToDo: This is still kind of work in progress, for 2D L2 thresholding is implemented and for 3D Linfty
subroutine cvs_threshold_tree(params, hvy_block, tree_ID, thresh, JMax_active)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params
    
    implicit none

    type (type_params), intent(in)  :: params
    !> heavy data array, WDed values in spaghetti form
    real(kind=rk), intent(inout)    :: hvy_block(:, :, :, :, :)
    integer(kind=ik), intent(in)    :: tree_ID
    real(kind=rk), intent(in)       :: thresh      !< Threshold for CVS to be applied to all WC
    integer(kind=ik), intent(in)    :: JMax_active !< if not provided it will be computed

    integer(kind=ik) :: hvy_ID, lgt_ID, ix, iy, iz, even_odd, i_b, level_me, JMaxActive
    real(kind=rk)    :: level_fac, wc_fac

    ! is the first point a SC or WC?
    even_odd = mod(params%g + 1, 2)

    if (present(JMax_active)) then
        JMaxActive = JMax_active
    else
        JMaxActive = maxActiveLevel_tree(tree_ID)
    endif

    do i_b = 1, hvy_n(tree_ID)
        hvy_ID = hvy_active(i_b, tree_ID)
        call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)
        level_me = lgt_block(lgt_ID, IDX_MESH_LVL)

        wc_fac = 1.0_rk
        if (params%eps_norm == "L2") then
            level_fac = ( 2.0_rk**(+dble((level_me-JMaxActive)*params%dim)/2.0_rk) )
        elseif (params%eps_norm == "H1") then
            level_fac = ( 2**(-level_me*(params%dim+2.0_rk)*0.5_rk) )
        else
            level_fac = 1.0
        endif
        
        if (params%dim == 2) then
            do iy = 1, params%Bs(2)+2*params%g
                do ix = 1, params%Bs(1)+2*params%g
                    ! skip SC
                    if (mod(ix, 2) == even_odd .and. mod(iy, 2) == even_odd) cycle

                    ! WC need to be renormalized
                    if (params%eps_norm == "L2") then
                        wc_fac = 0.5_rk  ! 2.0_rk**(-(params%dim)/2.0_rk)
                        if (mod(ix, 2) == 1-even_odd) wc_fac = wc_fac * 2.0_rk
                        if (mod(iy, 2) == 1-even_odd) wc_fac = wc_fac * 2.0_rk
                    endif

                    ! apply thresholding
                    if (abs(hvy_block(ix, iy, 1, 1, hvy_ID)) < thresh * level_fac * wc_fac) then
                        hvy_block(ix, iy, 1, 1, hvy_ID) = 0.0
                    endif
                enddo
            enddo
        else
            do iz = 1, params%Bs(3)+2*params%g
                do iy = 1, params%Bs(2)+2*params%g
                    do ix = 1, params%Bs(1)+2*params%g
                        ! skip SC
                        if (mod(ix, 2) == even_odd .and. mod(iy, 2) == even_odd .and. mod(iz, 2) == even_odd) cycle

                        ! WC need to be renormalized
                        if (params%eps_norm == "L2") then
                            wc_fac = 1/(sqrt(2.0_rk)*2.0_rk)  ! 2.0_rk**(-(params%dim)/2.0_rk)
                            if (mod(ix, 2) == 1-even_odd) wc_fac = wc_fac * 2.0_rk
                            if (mod(iy, 2) == 1-even_odd) wc_fac = wc_fac * 2.0_rk
                            if (mod(iz, 2) == 1-even_odd) wc_fac = wc_fac * 2.0_rk
                        endif

                        ! apply thresholding
                        if (abs(hvy_block(ix, iy, iz, 1, hvy_ID)) < thresh * level_fac * wc_fac) then
                            hvy_block(ix, iy, iz, 1, hvy_ID) = 0.0
                        endif
                    enddo
                enddo
            enddo
        endif
    enddo

end subroutine