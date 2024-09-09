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

    ! In order to start with leaf-wise investigation, these will receive the correct ref flag
    do k = 1, lgt_n(tree_ID)
        lgt_ID = lgt_active(k, tree_ID)
        lgt_block(lgt_ID, IDX_REFINE_STS) = REF_TMP_UNTREATED
    end do
    ! do backup here so we need less logic in synching neighbors
    do k = 1, hvy_n(tree_ID)
        hvy_ID = hvy_active(k, tree_ID)
        hvy_tmp(:,:,:,1:size(hvy_block, 4),hvy_ID) = hvy_block(:,:,:,1:size(hvy_block, 4),hvy_ID)
    enddo

    ! at first, initialize all mothers without any values yet
    t_block = MPI_Wtime()
    call cvs_init_mothers_tree(params, tree_ID)
    call toc( "adapt_tree (init mothers)", 102, MPI_Wtime()-t_block )
    ! now, all mothers have been created and share the ref status REF_TMP_EMPTY, leafs share status REF_TMP_UNTREATED

    ! iteration for wavelet decomposition - no blocks will be created or deleted during the loop!
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
        !     ! if (params%rank == 0 .and. lgt_block(lgt_id, IDX_REFINE_STS) == REF_TMP_UNTREATED .and. iteration >= 2) then
        !     ! if (lgt_block(lgt_id, IDX_MESH_LVL) == 3 .and. lgt_block(lgt_id, IDX_TC_2) == 54525952) then
        !         write(*, '("1 I-", i0, " R-", i2, " B-", i6, " L-", i2, " Ref-", i3, " TC-", i0)') iteration, params%rank, lgt_ID, &
        !             lgt_block(lgt_id, IDX_MESH_LVL), lgt_block(lgt_id, IDX_REFINE_STS), lgt_block(lgt_id, IDX_TC_2)
        !         ! call write_neighborhood_info(hvy_neighbor(hvy_id, :), params%dim)
        !         ! call dump_block_fancy(hvy_block(:, :, :, 1:1, hvy_id), "block_cvs_9.txt", params%Bs, params%g)
        !     ! endif
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

        ! Wavelet-transform all blocks which are untreated
        ! From now on until wavelet retransform hvy_block will hold the wavelet decomposed values in spaghetti form for these blocks
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
                ! We actually copy here a second time but this is to ensure we have the correct ghost points saved, and copying is not too expensive anyways
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

            ! ! when blocks don't have all neighbors available on medium lvl, the newly created mother block will not be able
            ! ! to find all needed neighbors, so we let the block wait if that is the case
            ! do k = 1, hvy_n(tree_ID)
            !     hvy_id = hvy_active(k, tree_ID)
            !     call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)
            !     if (lgt_block( lgt_ID, IDX_REFINE_STS) == -1) then
            !         call ensure_mother_neighbors(params, hvy_id, mark_TMP_flag=.true.)
            !     endif
            ! enddo
            ! ! after locally modifying refinement statusses, we need to synchronize light data
            ! call synchronize_lgt_data( params, refinement_status_only=.true. )

            ! watch for completeness
            do k = 1, hvy_n(tree_ID)
                hvy_id = hvy_active(k, tree_ID)
                call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
                call ensure_completeness( params, lgt_id, hvy_family(hvy_ID, 2:1+2**params%dim), mark_TMP_flag=.true. )
            enddo
            ! set mother block status, this is done here so that it will directly be synchronized with next sync_lgt call
            do k = 1, hvy_n(tree_ID)
                hvy_id = hvy_active(k, tree_ID)
                call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
                if (lgt_block( lgt_ID, IDX_REFINE_STS) == -1) then
                    lgt_block(hvy_family(hvy_ID, 1), IDX_REFINE_STS) = REF_TMP_UNTREATED
                endif
            enddo
            ! Details are synced so that blocks send to the correct neighbors
            call synchronize_lgt_data( params, refinement_status_only=.true. )

            ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            ! Create or update mother blocks with SC of daughter blocks
            ! This uses the already wavelet decomposed blocks and copies their SC into the new mother block
            t_block = MPI_Wtime()
            call sync_D2M(params, hvy_block, tree_ID, sync_case="ref", s_val=-1)
            call toc( "adapt_tree (update_mothers)", 111, MPI_Wtime()-t_block )

            ! now the daughter blocks are reset, lgt-loop to avoid synching
            do k = 1, lgt_n(tree_ID)
                lgt_ID = lgt_active(k, tree_ID)
                if (lgt_block( lgt_ID, IDX_REFINE_STS) == -1) lgt_block(lgt_ID, IDX_REFINE_STS) = 0
            enddo

            ! check if mother is happy to exist and can continue or if it wants to wait for all medium neighbors
            ! this also frees it of this burden if it ever happened and it now has all required neighbors
            ! it is the equivalent to ensure_gradedness
            do k = 1, hvy_n(tree_ID)
                hvy_id = hvy_active(k, tree_ID)
                call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
                ! Here we also directly backup all values for the first time
                if (lgt_block( lgt_ID, IDX_REFINE_STS) == REF_TMP_UNTREATED) then
                    hvy_tmp(:,:,:,1:size(hvy_block, 4),hvy_ID) = hvy_block(:,:,:,1:size(hvy_block, 4),hvy_ID)
                endif
                if (any(lgt_block( lgt_ID, IDX_REFINE_STS) == (/REF_TMP_UNTREATED, REF_TMP_UNTREATED_WAIT /) )) then
                    call block_can_WD_safely(params, hvy_id)
                endif
            enddo
            ! Details are synced so that blocks send to the correct neighbors
            call synchronize_lgt_data( params, refinement_status_only=.true. )

            write(toc_statement, '(A, i0, A)') "adapt_tree (WD it ", iteration, " update_mothers)"
            call toc( toc_statement, 1650+iteration, MPI_Wtime()-t_block )
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

        if (present(verbose_check)) then
            write(*, '(A, i0, A, i0, A)') "Loop ", iteration, " with ", lgt_n(tree_ID), " blocks"
        endif

        ! ! for debugging purposes
        ! write( toc_statement,'(A, i6.6, A)') 'TestWD_', iteration, "000000.h5"
        ! call saveHDF5_tree(toc_statement, dble(iteration), iteration, 1, params, hvy_block, tree_ID )
        ! write( toc_statement,'(A, i6.6, A)') 'TestN_', iteration, "000000.h5"
        ! call saveHDF5_tree(toc_statement, dble(iteration), iteration, 1, params, hvy_tmp, tree_ID )
    end do

    ! ! for debugging purposes
    ! write( toc_statement,'(A, i6.6, A)') 'TestWD_', 100, "000000.h5"
    ! call saveHDF5_tree(toc_statement, dble(100), 100, 1, params, hvy_block, tree_ID )
    ! write( toc_statement,'(A, i6.6, A)') 'TestN_', 100, "000000.h5"
    ! call saveHDF5_tree(toc_statement, dble(100), 100, 1, params, hvy_tmp, tree_ID )

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

        ! ! Normal adapt_tree loop does another CE_modify here to copy the SC, but for full WT this is not necessary
        ! ! as interior blocks only get values from their medium neighbors and no SC need to be copied
        ! t_block = MPI_Wtime()
        ! call coarse_extension_modify(params, hvy_block, hvy_tmp, tree_ID, CE_case="level", s_val=level, clear_wc_only=.false.)
        ! call toc( "adapt_tree (coarse_extension_modify)", 105, MPI_Wtime()-t_block )

        ! synch SC and WC from coarser neighbours and same-level neighbours in order to apply the correct wavelet reconstruction
        ! Attention: This uses hvy_temp for coarser neighbors to predict the data, as we want the correct SC from coarser neighbors
        t_block = MPI_Wtime()
        call sync_SCWC_from_MC( params, hvy_block, tree_ID, hvy_tmp, g_minus=g_spaghetti, g_plus=g_spaghetti, level=level)
        call toc( "adapt_tree (sync all <- MC)", 112, MPI_Wtime()-t_block )

        write(toc_statement, '(A, i0, A)') "adapt_tree (WR it ", iteration, " sync level <- MC)"
        call toc( toc_statement, 1750+iteration, MPI_Wtime()-t_block )

        ! In theory, for full WT this CE_modify is only here for two reasons:
        !    1. We need to clear the WC from coarser neighbors, as the values in there are interpolated SC and we need to set those to the correct WC being 0
        !    2. interior WC are deleted so that the SC copy for syncing later is exact with the original values, we assume those WC are not significant
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
        call sync_M2D(params, hvy_block, tree_ID, sync_case="level", s_val=level)
        call toc( "adapt_tree (sync level 2 daughters)", 114, MPI_Wtime()-t_block )

        write(toc_statement, '(A, i0, A)') "adapt_tree (WR it ", iteration, " sync level 2 daughters)"
        call toc( toc_statement, 1900+iteration, MPI_Wtime()-t_block )

        write(toc_statement, '(A, i0, A)') "adapt_tree (WR it ", iteration, " TOTAL)"
        call toc( toc_statement, 1950+iteration, MPI_Wtime()-t_loop )

        iteration = iteration + 1
    enddo
end subroutine



!> \brief This functions inits all mother blocks from a leaf-only grid down to JMin
subroutine cvs_init_mothers_tree(params, tree_ID, verbose_check)
    implicit none

    type (type_params), intent(in)      :: params
    integer(kind=ik), intent(in)        :: tree_ID
    logical, intent(in), optional       :: verbose_check  !< No matter the value, if this is present we debug

    integer(kind=ik)     :: i_level, j, k, N, lgt_ID, hvy_ID, level_b, iteration, level, data_rank, lgt_merge_id, hvy_merge_id, hvy_n_old, digit
    real(kind=rk)        :: t_block, t_loop
    character(len=clong) :: toc_statement
    logical              :: iterate
    ! list of block ids, proc ranks
    integer(kind=ik)                    :: lgt_daughters(1:8), rank_daughters(1:8)
    integer(kind=tsize)                 :: treecode

    ! number of blocks to merge, 4 or 8
    N = 2**params%dim
    ! loop downwards until all blocks are on Jmin
    level       = maxActiveLevel_tree(tree_ID) ! leaf-wise blocks can have this as maximum level
    iteration   = 0

    ! J: I have decided for the following way of creating mother blocks
    !   - main loop goes refine-evolve-coarsen, so first iteration will work on all available block
    !     it is the one with the most data to send and it's nice that we can choose any of the ranks here (as they are all the same)
    !   - after first iteration, not all sister-information are available, I just create the following mothers after a deterministic way
    !     and choose the rank of the block with the highest digit, this ensures some form of spread of ranks, after all I do not really care anymore about the amount of data

    ! just loop until hvy_n does not change anymore, that means all possible mothers have been created
    do while( hvy_n_old /= hvy_n(tree_ID) )
        !---------------------------------------------------------------------------
        ! create new empty blocks on rank with block with digit = 2**dim-1
        !---------------------------------------------------------------------------
        do k = 1, hvy_n(tree_ID)
            hvy_ID = hvy_active(k, tree_ID)
            call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)
            level_b = lgt_block( lgt_ID, IDX_MESH_LVL )

            ! check if block does not have a mother yet and blocks are not on minimum level
            if ( hvy_family(hvy_ID, 1) == -1 .and. level_b /= params%Jmin) then
                ! Get digit
                treecode = get_tc(lgt_block( lgt_ID, IDX_TC_1:IDX_TC_2 ))
                digit = tc_get_digit_at_level_b(treecode, params%dim, level_b, params%Jmax)
                ! The merging will be done on the rank of the block with digit 2**dim-1
                if (digit == 2**params%dim -1) then
                    ! construct new mother block if on my rank and create light data entry for the new block
                    call get_free_local_light_id(params, params%rank, lgt_merge_id, message="init_mothers")
                    call lgt2hvy(hvy_merge_id, lgt_merge_id, params%rank, params%number_blocks)
                    ! set meta_data of mother block
                    lgt_block( lgt_merge_id, : ) = -1
                    call set_tc(lgt_block( lgt_merge_id, IDX_TC_1:IDX_TC_2), tc_clear_until_level_b(treecode, &
                        dim=params%dim, level=level_b-1, max_level=params%Jmax))
                    lgt_block( lgt_merge_id, IDX_MESH_LVL ) = level_b-1
                    lgt_block( lgt_merge_id, IDX_REFINE_STS ) = REF_TMP_EMPTY
                    lgt_block( lgt_merge_id, IDX_TREE_ID ) = tree_ID
                    ! make sure family array for mother block is cleared
                    hvy_family(hvy_merge_id, :) = -1
                    
                    ! update myself that I have created my own mother
                    ! sisters will not be updated but as only one daughter will surely create the mother there will be no doubles
                    hvy_family(hvy_ID, 1) = lgt_merge_id
                endif
            endif
        enddo

        ! we need to update hvy_active locally to be able to loop correctly
        hvy_n_old = hvy_n(tree_ID)
        call createHvyActive_tree(params, tree_ID)

        iteration = iteration + 1

        if (present(verbose_check)) then
            write(*, '(A, i0, A, i0, A, i0)') "R", params%rank, " Deterministic mother addition it ", iteration, ", hvy_n= ", hvy_n(tree_ID)
        endif
    enddo

    ! the active lists are outdated, so lets resynch
    call synchronize_lgt_data( params, refinement_status_only=.false.)
    ! At last, we update all full neighbor relations only once
    call updateMetadata_tree(params, tree_ID, search_overlapping=.true.)
end subroutine



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

    integer(kind=ik) :: hvy_ID, lgt_ID, ix, iy, iz, even_odd, i_b, level_me
    real(kind=rk)    :: level_fac, wc_fac

    ! is the first point a SC or WC?
    even_odd = mod(params%g + 1, 2)

    do i_b = 1, hvy_n(tree_ID)
        hvy_ID = hvy_active(i_b, tree_ID)
        call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)
        level_me = lgt_block(lgt_ID, IDX_MESH_LVL)

        wc_fac = 1.0_rk
        if (params%eps_norm == "L2") then
            level_fac = ( 2.0_rk**(+dble((level_me-JMax_Active)*params%dim)/2.0_rk) )
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