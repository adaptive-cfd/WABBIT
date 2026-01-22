subroutine multigrid_solve(params, hvy_sol, hvy_RHS, hvy_work, tree_ID, init_0, verbose)
    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    real(kind=rk), intent(inout)       :: hvy_sol(:, :, :, :, :), hvy_RHS(:, :, :, :, :), hvy_work(:, :, :, :, :)
    integer(kind=ik), intent(in)       :: tree_ID
    logical, intent(in), optional      :: init_0
    logical, intent(in), optional      :: verbose

    integer(kind=ik)                   :: k_block, lgt_ID, hvy_id, ic, i_cycle
    real(kind=rk)                      :: dx(1:3), x0(1:3), residual(1:4*size(hvy_sol,4)), norm_sol(1:size(hvy_sol,4))
    real(kind=rk)                      :: t_block
    logical                            :: verbose_apply, init0

    verbose_apply = .false.
    if (present(verbose)) verbose_apply = verbose
    init0 = .false.
    if (present(init_0)) init0 = init_0

    ! For the Mehrstellenverfahren, we actually do solve the system Au = Bb
    ! So before doing anything, we apply the matrix B on b, this is a large tensorial matrix, but we only have to apply this operation once
    if (params%poisson_order == "FD_6th_mehrstellen") then
        t_block = MPI_Wtime()
        do k_block = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k_block, tree_ID)
            call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
            call get_block_spacing_origin( params, lgt_ID, x0, dx )
            do ic = 1, size(hvy_RHS,4)
                ! apply the stencil to the RHS
                call GS_compute_Ax(params, hvy_RHS(:,:,:,ic,hvy_id), hvy_work(:,:,:,ic,hvy_id), dx, apply_B_RHS=.true.)
            enddo
            hvy_RHS(:,:,:,1:size(hvy_RHS,4),hvy_id) = hvy_work(:,:,:,1:size(hvy_RHS,4),hvy_id)
        enddo
        call toc( "RHS Preparation Bb", 10000, MPI_Wtime()-t_block )

        t_block = MPI_Wtime()
        call sync_ghosts_tree(params, hvy_RHS(:,:,:,1:size(hvy_RHS,4),:), tree_ID)
        call toc( "Sync Layer", 10010, MPI_Wtime()-t_block )
    endif

    ! compute norm of RHS for relative tolerances
    call componentWiseNorm_tree(params, hvy_RHS(:,:,:,1:size(hvy_RHS,4),:), tree_ID, "Linfty", norm_sol(1:size(hvy_RHS,4)), threshold_state_vector=.false.)

    ! prepare full tree grid, this will be populated along the way
    ! blocks are practically empty, but we fill them along the way so ref flag will be set to 0 to allow synching
    call init_full_tree(params, tree_ID, set_ref=0)

    ! init solution on lower levels as zero, elsewise it is that from last time-step
    do k_block = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k_block, tree_ID)
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
        if (.not. block_is_leaf(params, hvy_ID) .or. init0) then
            hvy_sol(:,:,:,1:size(hvy_RHS,4),hvy_id) = 0.0_rk
        endif
    enddo

    ! compute several v- or f-cycles
    do i_cycle = 1,params%poisson_cycle_max_it

        ! Call the multigrid vcycle function
        call multigrid_vcycle(params, hvy_sol, hvy_RHS, hvy_work, tree_ID, residual_out=residual(1:3*size(hvy_sol,4)), verbose=verbose_apply)
        ! now the residuals are stored as: Linfty, L2, L1, but L2 and L1 are only computed with verbose_apply

        ! if (params%rank == 0 .and. .not. verbose_apply) write(*, '(A, i0, A, i0, A, 10(es10.4, 1x))') "   it ", i_cycle, "/", params%poisson_cycle_it, " Residual Linfty: ", residual(1:size(hvy_sol,4))

        ! choose between different end criteria:
        ! fixed_iterations - do a set number of V-cycles
        ! tolerance - stop when absolute or relative residuals are below a certain threshold or maximum number of iterations is reached
        if (params%poisson_cycle_end_criteria == "fixed_iterations") then
            if (i_cycle >= params%poisson_cycle_it) exit
        elseif (params%poisson_cycle_end_criteria == "tolerance") then
            if (all(residual(1:size(hvy_sol,4)) < params%poisson_cycle_tol_abs)) exit
            if (all(residual(1:size(hvy_sol,4))/norm_sol(1:size(hvy_sol,4)) < params%poisson_cycle_tol_rel) .and. all(norm_sol(1:size(hvy_sol,4)) > 1e-8)) exit
        else
            call abort(250903, "Don't know how to stop! Please choose between 'fixed_iterations' and 'tolerance', thank you :)")
        endif

    enddo
    ! fortran does one last increase and check for do-loops, so I change it back to display correct amount to user
    if (i_cycle > params%poisson_cycle_max_it) i_cycle = params%poisson_cycle_max_it

    ! laplacian is invariant to shifts of constant values
    ! our values are defined with zero mean for comparison
    ! as multigrid might accidently introduce a constant offset, we remove it after all cycles completed
    call componentWiseNorm_tree(params, hvy_sol(:,:,:,1:1,:), tree_ID, "Mean", residual(3*size(hvy_sol,4)+1:4*size(hvy_sol,4)), threshold_state_vector=.false.)
    do k_block = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k_block, tree_ID)
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
        do ic = 1, size(hvy_RHS,4)
            hvy_sol(:,:,:,ic,hvy_id) = hvy_sol(:,:,:,ic,hvy_id) - residual(3*size(hvy_sol,4)+ic)
        enddo
    enddo
    if (params%rank == 0 .and. verbose_apply) write(*, '(A, 10(es10.3, 1x))') "--- Mean value: ", residual(3*size(hvy_sol,4)+1:4*size(hvy_sol,4))

    if (params%rank == 0 .and. .not. verbose_apply) write(*, '(A, i0, A, 10(es10.4, 1x))') "   Final Residual after ", i_cycle, " it, Linfty: ", residual(1:size(hvy_sol,4))

    ! delete all non-leaf blocks with daughters as we for now do not have any use for them
    call prune_fulltree2leafs(params, tree_ID)

end subroutine multigrid_solve



!> \brief Do a multigrid vcycle
!> \details This function executes a multigrid vcycle. The structure is as following:
!>    - It first goes down to the lowest level with downwards sweeps
!>    - It solves the solution on the coarsest level
!>    - It then goes back up to the final level with upwards sweeps
!  A v-cycle for level 4 would look like the following:
!
! d           u
!   d       u
!     d   u
!       c
!
! c - coarsest level solve (FFT, CG or GS)
! u - upwards prolongation and Gauss-Seidel sweeps with it_GS
! d - Decomposition with restriction filter
subroutine multigrid_vcycle(params, hvy_sol, hvy_RHS, hvy_work, tree_ID, verbose, residual_out)
    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    real(kind=rk), intent(inout)       :: hvy_sol(:, :, :, :, :), hvy_RHS(:, :, :, :, :), hvy_work(:, :, :, :, :)
    integer(kind=ik), intent(in)       :: tree_ID
    real(kind=rk), intent(inout)       :: residual_out(:)
    logical, intent(in), optional      :: verbose

    integer(kind=ik)                   :: k_block, lgt_ID, hvy_id
    integer(kind=ik)                   :: i_level, Jmax_a, it_coarsest, i_sweep, nc, ic, mpierr
    logical                            :: sweep_forward

    character(len=cshort)              :: fname
    real(kind=rk)                      :: dx(1:3), x0(1:3), domain(1:3), tol_cg
    integer(kind=2)                    :: n_domain(1:3)
    integer(kind=tsize)                :: treecode
    logical                            :: verbose_apply

    real(kind=rk)        :: t_block, t_loop, t_cycle, t_print(1:1)

    verbose_apply = .false.
    if (present(verbose)) verbose_apply = verbose

    nc = size(hvy_sol,4)

    ! Logic to define the solution layer at each iteration:
    !    Downwards - we simply do a wavelet decomposition of the residual, but only using the scaling filter
    !    Upwards   - we reconstruct similarly to reconstruction, but always keep all leaf-blocks of lower levels to keep a full periodic grid
    !                after each reconstruction we do the Gauss-Seidel sweeps/iterations

    ! ------------------------------------------------------------
    ! Gauss-Seidel Multi-Grid - solve the system Ax = b
    !    When going to lower levels, we compute the residual r = b - Ax and solve on lower levels Ae = r
    !    We decompose the residual until lowest level, solve using FFT, CG or GS there and then prolongate upwards and solve using GS

    t_cycle = MPI_Wtime()

    Jmax_a = maxActiveLevel_tree(tree_id)

    ! sync solution before computing residual
    t_block = MPI_Wtime()
    ! call sync_ghosts_RHS_tree(params, hvy_sol(:,:,:,1:nc,:), tree_ID, g_minus=a, g_plus=a)
    ! call sync_ghosts_tree(params, hvy_sol(:,:,:,1:nc,:), tree_ID, g_minus=a, g_plus=a)
    call sync_ghosts_tree(params, hvy_sol(:,:,:,1:nc,:), tree_ID)
    call toc( "Sync Layer", 10010, MPI_Wtime()-t_block )

    ! JB Comment: I think this is not necessary, we do coarse extension but only keep the block values and not the decomposed ones, so we should reconstruct the original values later on without any problems
    ! ! we need to preserve b
    ! ! Usually we could restore it later with b = r + Ax, however, with coarse extension we have to alter the residual and would not get back the original b
    !     ! downwards - no sweeps are done, we only restrict the residual down to the lowest level
    ! ! compute the residual r = b - Ax to pass it downwards
    ! t_block = MPI_Wtime()
    ! do k_block = 1, hvy_n(tree_ID)
    !     hvy_id = hvy_active(k_block, tree_ID)
    !     call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

    !     ! b only defined on leaf layer
    !     if (.not. block_is_leaf(params, hvy_id)) cycle

    !     ! store b
    !     hvy_work(:,:,:,nc+1:2*nc,hvy_id) = hvy_RHS(:,:,:,1:nc,hvy_id)
    ! enddo
    ! call toc( "GS Downwards - Store b", 10020, MPI_Wtime()-t_block )


    ! downwards - no sweeps are done, we only restrict the residual down to the lowest level
    ! compute the residual r = b - Ax to pass it downwards
    t_block = MPI_Wtime()
    do k_block = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k_block, tree_ID)
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

        ! we start on the leaf layer, here we need to compute residual to pass it downwards
        if (.not. block_is_leaf(params, hvy_id)) cycle

        ! get spacing
        call get_block_spacing_origin( params, lgt_ID, x0, dx )

        ! compute actual residual r = b - Ax, 
        call GS_compute_residual(params, hvy_sol(:,:,:,1:nc,hvy_id), hvy_RHS(:,:,:,1:nc,hvy_id), hvy_work(:,:,:,1:nc,hvy_id), dx)

    enddo
    call toc( "GS Downwards - Compute residual", 10021, MPI_Wtime()-t_block )

    ! we use the decompose function to pass down the residual to all levels
    ! this overwrites b of the leaf layer, but we will recover it later
    ! solution in non-decomposed form is present in hvy_RHS
    t_block = MPI_Wtime()
    call wavelet_decompose_full_tree(params, hvy_work(:,:,:,1:nc,:), tree_ID, hvy_RHS, init_full_tree_grid=.false., compute_SC_only=.true., scaling_filter=params%MGR)
    call toc( "GS Downwards - full sweep", 10002, MPI_Wtime()-t_block )
    if (verbose_apply) then
        t_block = MPI_Wtime()-t_block
        call MPI_ALLREDUCE(MPI_IN_PLACE, t_block, 1, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)
        call append_t_file('multigrid-iteration.t', (/-1.0_rk, -1.0_rk, t_block/))
    endif

    ! usually after decomposition the residual is not synced, when doing multi-step GS sweeps, we need to sync it
    t_block = MPI_Wtime()
    call sync_ghosts_tree(params, hvy_RHS(:,:,:,1:nc,:), tree_ID)
    call toc( "Sync Layer", 10010, MPI_Wtime()-t_block )

    ! for all levels we are not solving Au = b, but rather Ae=r
    ! only on the last upwards step will we iterate on the original problem, so we need to backup u
    ! starting from here, hvy_work contains the backed-up solution
    t_block = MPI_Wtime()
    do k_block = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k_block, tree_ID)
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

        ! only leaf layer has solution
        if (.not. block_is_leaf(params, hvy_id)) cycle
        hvy_work(:,:,:,1:nc,hvy_id) = hvy_sol(:,:,:,1:nc,hvy_id)
    enddo
    call toc( "GS Downwards - Backup solution", 10022, MPI_Wtime()-t_block )

    ! lowest block is treated differently, here we try to solve more thoroughly. Following options exist:
    !    - Conjugent gradient method : only for level 0 for now
    !    - GS sweeps : very inefficient as blocks are usually too large for convergence
    !    - FFT :only for level 0, needs to be checked to solve true solution or discretized laplacian
    call multigrid_coarsest(params, hvy_sol, hvy_RHS, tree_ID, params%Jmin, Jmax_a, verbose_apply)

    ! upwards sweeps - prolong lower level solution, then solve Gauss-Seidel iterations
    ! do i_level = Jmin+1, Jmax_a
    !     call multigrid_upwards(params, hvy_sol, hvy_RHS, hvy_work, tree_ID, i_level, Jmax_a, hvy_depth)
    ! enddo
    call multigrid_upwards(params, hvy_sol, hvy_RHS, hvy_work, tree_ID, params%Jmin, Jmax_a, verbose_apply)

    ! We do not solve Ae=r but Au=b and we need to restore the old solution as well as the RHS b
    ! this is happening for all leaf-blocks
    t_block = MPI_Wtime()
    do k_block = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k_block, tree_ID)
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

        if (block_is_leaf(params, hvy_id)) then
            ! get spacing
            call get_block_spacing_origin( params, lgt_ID, x0, dx )

            ! recompute RHS b = r + Ax, 
            call GS_compute_residual(params, hvy_work(:,:,:,1:nc,hvy_id), hvy_RHS(:,:,:,1:nc,hvy_id), hvy_RHS(:,:,:,1:nc,hvy_id), dx, recompute_b=.true.)
            ! ! restore full RHS b
            ! hvy_RHS(:,:,:,1:nc,hvy_id) = hvy_work(:,:,:,nc+1:2*nc,hvy_id)
    
            ! add previous solution to reconstructed solution
            hvy_sol(:,:,:,1:nc,hvy_id) = hvy_sol(:,:,:,1:nc,hvy_id) + hvy_work(:,:,:,1:nc,hvy_id)
        endif
    enddo

    ! sync RHS for sync_freq > 1
    t_block = MPI_Wtime()
    call sync_ghosts_tree(params, hvy_RHS(:,:,:,1:nc,:), tree_ID)
    call toc( "Sync Layer", 10010, MPI_Wtime()-t_block )
    call toc( "GS - RHS and solution restoration", 10032, MPI_Wtime()-t_block )

    ! sync before computing final residual
    call sync_ghosts_tree(params, hvy_sol(:,:,:,1:nc,:), tree_ID, params%poisson_stencil_size, params%poisson_stencil_size)
    ! call sync_ghosts_tree(params, hvy_sol(:,:,:,1:nc,:), tree_ID)
    call toc( "Sync Layer", 10010, MPI_Wtime()-t_block )

    ! compute residual only on leaf layer
    t_block = MPI_Wtime()
    do k_block = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k_block, tree_ID)
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

        ! get spacing
        call get_block_spacing_origin( params, lgt_ID, x0, dx )

        ! treecode = get_tc(lgt_block(lgt_id, IDX_TC_1 : IDX_TC_2))
        ! write(fname, '(A, i0, A,i0, A,i0, A, i0, A)') "block_dumped_lvl", lgt_block(lgt_id, IDX_MESH_LVL),"_Vcycl", i_cycle ,"_tc=", treecode, ".t"
        ! call dump_block_fancy(hvy_sol(:,:,:,1:nc,hvy_id), fname, params%Bs, params%g, digits=2)

        ! write(*, '("R-", i2, " B-", i6, " L-", i2, " Ref-", i3, " TC-", i0)') params%rank, lgt_ID, &
        !     lgt_block(lgt_id, IDX_MESH_LVL), lgt_block(lgt_id, IDX_REFINE_STS), lgt_block(lgt_id, IDX_TC_2)

        ! skip blocks not on this level
        if (.not. block_is_leaf(params, hvy_id)) cycle

        ! compute actual residual r = b - Ax
        call GS_compute_residual(params, hvy_sol(:,:,:,1:nc,hvy_id), hvy_RHS(:,:,:,1:nc,hvy_id), hvy_work(:,:,:,1:nc,hvy_id), dx)
    enddo
    ! we are subtracting the mean value from the residual in case they have been introduced
    call componentWiseNorm_tree(params, hvy_work(:,:,:,1:nc,:), tree_ID, "Mean", residual_out, threshold_state_vector=.false.)
    if (params%rank == 0 .and. verbose_apply) write(*, '(A, es10.3, A)') "--- Residual Mean: ", residual_out(1), " ---"
    do k_block = 1, hvy_n(tree_ID)
        do ic = 1,nc
            hvy_id = hvy_active(k_block, tree_ID)
            hvy_work(:,:,:,ic,hvy_id) = hvy_work(:,:,:,ic,hvy_id) - residual_out(ic)
        enddo
    enddo

    ! ***************************
    ! compute global residuals
    ! ***************************
    ! we always compute the Linfty norm, and if needed the L1 and L2 norm as well, printing them to file

    call componentWiseNorm_tree(params, hvy_work(:,:,:,1:nc,:), tree_ID, "Linfty", residual_out, threshold_state_vector=.false.)

    if (verbose_apply) then
        call componentWiseNorm_tree(params, hvy_work(:,:,:,1:nc,:), tree_ID, "L2", residual_out(nc+1:2*nc), threshold_state_vector=.false.)
        if (params%rank == 0) write(*, '(A, es10.4, A)') "--- Residual L2: ", residual_out(nc+1), " ---"
        call componentWiseNorm_tree(params, hvy_work(:,:,:,1:nc,:), tree_ID, "L1", residual_out(2*nc+1:3*nc), threshold_state_vector=.false.)
        ! if (params%rank == 0) write(*, '(A, es10.4, A)') "--- Residual L1: ", residual_out(2*nc+1), " ---"

        if (params%rank == 0) write(*, '(A, es10.4, A)') "--- Residual Linfty: ", residual_out(1), " ---"


        t_print = MPI_Wtime()-t_cycle
        call MPI_ALLREDUCE(MPI_IN_PLACE, t_print, 1, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)
        call append_t_file('multigrid-cycle.t', (/residual_out(1), residual_out(nc+1), residual_out(2*nc+1), t_print(1)/))
    endif

    call toc( "Final residual", 10014, MPI_Wtime()-t_block )

    call toc( "V cycle", 10001, MPI_Wtime()-t_cycle )

end subroutine multigrid_vcycle



!> \brief Upwards level iteration for the multigrid solver
!> \details Compute for one level an upwards step
!> Does three things:
!>   - sync the coarser solution from the lower lewel as WD decomposed values, reconstruct the coarser solution
!>   - do GS sweeps on this level
subroutine multigrid_upwards(params, hvy_sol, hvy_RHS, hvy_work, tree_ID, Jmin, Jmax_a, verbose)
    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    real(kind=rk), intent(inout)       :: hvy_sol(:, :, :, :, :), hvy_RHS(:, :, :, :, :), hvy_work(:, :, :, :, :)
    integer(kind=ik), intent(in)       :: tree_ID
    integer(kind=ik), intent(in)       :: Jmin
    integer(kind=ik), intent(in)       :: Jmax_a
    logical, intent(in), optional      :: verbose

    character(len=cshort)              :: fname
    integer(kind=tsize)                :: treecode
    real(kind=rk)                      :: norm(1:2*size(hvy_sol,4))

    logical            :: sweep_forward
    integer(kind=ik)   :: k_block, lgt_ID, hvy_id, nc, ic, i_sweep, mpierr, i_level, level_me, ref_me, sync_freq, i_g(1:3), i_BS(1:3)
    real(kind=rk)      :: t_loop, t_block, t_print(1:1)
    real(kind=rk)      :: dx(1:3), x0(1:3)
    logical            :: verbose_apply

    verbose_apply = .false.
    if (present(verbose)) verbose_apply = verbose

    nc = size(hvy_sol,4)

    i_g = 0
    i_g(1:params%dim) = mod(params%g, 2)
    i_BS = 0
    i_BS(1:params%dim) = mod(params%Bs(1:params%dim), 2)

    ! prepare grid, set all ref stats to EMPTY and only lowest level to 1
    do k_block = 1, lgt_n(tree_ID)
        lgt_ID = lgt_active(k_block, tree_ID)
        level_me = lgt_block( lgt_ID, IDX_MESH_LVL )
        if (level_me == Jmin) then
            lgt_block(lgt_ID, IDX_REFINE_STS) = 1
        else
            lgt_block(lgt_ID, IDX_REFINE_STS) = REF_TMP_EMPTY
        endif
    end do

    ! loop until highest level, we always advance level-wise to get the quickest to max level
    ! however, we need full grids, so leaf blocks will be kept on this level
    do i_level = Jmin+1, Jmax_a

        t_loop = MPI_Wtime()

        if (params%rank == 0 .and. verbose_apply) write(*, '(A, 2(A, i0))') repeat('  ', i_level+1), 'Upwards GS sweep lvl ', i_level, ' it ', params%poisson_GS_it

        ! Preparation for this level:
        !   - ref = -1: last level that is ready to go upwards
        !   - ref =  0: blocks on last or lower levels that stay as they are leaf blocks
        !   - ref =  1: new blocks on this level that will receive values from -1
        !   - ref = EMPTY: all other blocks'/;.l
        t_block = MPI_Wtime()
        do k_block = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k_block, tree_ID)
            call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
            level_me = lgt_block(lgt_id, IDX_MESH_LVL)
            ref_me = lgt_block(lgt_id, IDX_REFINE_STS)

            ! change refinement stats from previous iteration
            !   ref = -1: blocks that are ready to go upwards if not leaf blocks
            !   ref =  0: leaf blocks that stay
            if (ref_me ==  1) then
                if (block_is_leaf(params, hvy_id)) then
                    lgt_block(lgt_id, IDX_REFINE_STS) = 0
                else
                    lgt_block(lgt_id, IDX_REFINE_STS) = -1
                endif
            endif

            ! set new blocks on this level with ref 1
            if (level_me == i_level) then
                lgt_block(lgt_id, IDX_REFINE_STS) = 1
            endif
        enddo
        ! synch the refinement stati
        call synchronize_lgt_data( params, refinement_status_only=.true.)
        call toc( "GS Upwards - Level definition", 10030, MPI_Wtime()-t_block )

        ! do upwards sync of solution field, this then needs to be interpolated
        t_block = MPI_Wtime()
        call sync_M2D(params, hvy_sol(:,:,:,1:nc,:), tree_ID, sync_case="ref", s_val=-1)
        call toc( "Sync M2D", 10012, MPI_Wtime()-t_block )

        ! reordering of refinement flags:
        !   -1 -> REF_TMP_EMPTY : these blocks have synched and are finished
        do k_block = 1, lgt_n(tree_ID)
            lgt_id = lgt_active(k_block, tree_ID)
            if (lgt_block(lgt_id,IDX_REFINE_STS) == -1) then
                lgt_block(lgt_id,IDX_REFINE_STS) = REF_TMP_EMPTY
            endif
        enddo

        ! we need to synch the solution to reconstruct it, this now only synchs on the current solution layer, as the rest is EMPTY
        ! maybe this syncing size can be reduced to params%HR, needs to be checked
        t_block = MPI_Wtime()
        call sync_ghosts_tree(params, hvy_sol(:,:,:,1:nc,:), tree_ID)
        call toc( "Sync Layer", 10010, MPI_Wtime()-t_block )

        t_block = MPI_Wtime()
        do k_block = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k_block, tree_ID)
            call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

            ! only blocks with 1 received new values, so only these are reconstructed
            if (lgt_block(lgt_id,IDX_REFINE_STS) /= 1) cycle

            do ic = 1,nc
                ! We now predict, here we insert the following, example BS=6, g=2:
                !    GS GW IS IW IS IW IS IW GS GW
                !    GS    IS    IS    IS    GS
                ! We use the grid between the outmost GS to interpolate and retrieve the grid with points in between afterwards
                ! Leftmost GS is dependent on g even(1)/odd(2)
                ! Rightmost GS is dependent on BS even g even(N-1)/odd(N) and BS odd g even(N)/odd(N-1)
                if (params%dim==2) then
                    call prediction(hvy_sol(1+i_g(1):params%Bs(1)+2*params%g:2,1+i_g(2):params%Bs(2)+2*params%g:2,:,ic,hvy_id), &
                                    hvy_sol(1+i_g(1):params%Bs(1)+2*params%g-1+mod(i_g(1)+i_BS(1),2),1+i_g(2):params%Bs(2)+2*params%g-1+mod(i_g(2)+i_BS(2),2),:,ic,hvy_id), params%order_predictor)
                else
                    call prediction(hvy_sol(1+i_g(1):params%Bs(1)+2*params%g:2,1+i_g(2):params%Bs(2)+2*params%g:2,1+i_g(3):params%Bs(3)+2*params%g:2,ic,hvy_id), &
                                    hvy_sol(1+i_g(1):params%Bs(1)+2*params%g-1+mod(i_g(1)+i_BS(1),2),1+i_g(2):params%Bs(2)+2*params%g-1+mod(i_g(2)+i_BS(2),2),1+i_g(3):params%Bs(3)+2*params%g-1+mod(i_g(3)+i_BS(3),2),ic,hvy_id), params%order_predictor)
                endif
            enddo

            ! treecode = get_tc(lgt_block(lgt_id, IDX_TC_1 : IDX_TC_2))
            ! write(fname, '(A, i0, A,i0, A,i0, A, i0, A)') "block_dumped_lvl", lgt_block(lgt_id, IDX_MESH_LVL),"_it", i_level ,"_tc=", treecode, ".t"
            ! call dump_block_fancy(hvy_work(:,:,:,hvy_depth(k_block)*nc+1:hvy_depth(k_block)*nc+nc,hvy_id,1), fname, params%Bs, params%g, digits=2, print_ghosts=.false.)
            ! call dump_block_fancy(hvy_sol(:,:,:,1:nc,hvy_id), fname, params%Bs, params%g, digits=2, print_ghosts=.false.)
        enddo
        call toc( "GS Upwards - Solution prolongation", 10031, MPI_Wtime()-t_block )

        ! laplacian is invariant to shifts of constant values
        ! our values are defined with zero mean for comparison
        ! as multigrid might accidently introduce a constant offset, we remove it
        call componentWiseNorm_tree(params, hvy_RHS(:,:,:,1:nc,:), tree_ID, "Mean", norm(1:nc), norm_case="not_empty", threshold_state_vector=.false.)
        if (params%rank == 0 .and. verbose_apply) write(*, '(A, A, es10.3)') repeat('  ', i_level+1), 'RHS mean ', norm(1)
        do k_block = 1, hvy_n(tree_ID)
            do ic = 1,nc
                hvy_id = hvy_active(k_block, tree_ID)
                hvy_RHS(:,:,:,ic,hvy_id) = hvy_RHS(:,:,:,ic,hvy_id) - norm(ic)
            enddo
        enddo

        ! do actual sweeps
        ! if (i_level /= Jmax_a) then
            sync_freq = params%poisson_Sync_it
            ! sync_freq = params%g / params%poisson_stencil_size
            sweep_forward = .true.
            do i_sweep = 1, params%poisson_GS_it
                if (modulo(i_sweep-1, sync_freq) == 0) then
                    ! we need to synch before a sweep
                    t_block = MPI_Wtime()
                    call sync_ghosts_tree(params, hvy_sol(:,:,:,1:nc,:), tree_ID, params%poisson_stencil_size, params%poisson_stencil_size, ignore_Filter=(params%poisson_stencil_size<=ubound(params%HD, dim=1)/2 .or. .not. params%isLiftedWavelet))
                    ! call sync_ghosts_tree(params, hvy_sol(:,:,:,1:nc,:), tree_ID, a*sync_freq, a*sync_freq)
                    call toc( "Sync Layer", 10010, MPI_Wtime()-t_block )
                endif

                ! blocks on this iteration do a GS-sweep, they have refinement status 0 or 1
                call GS_iteration_ref(params, tree_id, (/ 1, 0 /), hvy_sol(:,:,:,1:nc,:), hvy_RHS(:,:,:,1:nc,:), sweep_forward, filter_offset=params%g)
                ! call GS_iteration_ref(params, tree_id, (/ 1, 0 /), hvy_sol(:,:,:,1:nc,:), hvy_RHS(:,:,:,1:nc,:), sweep_forward, filter_offset=max(0,params%g-(sync_freq-1)*params%poisson_stencil_size))

                sweep_forward = .not. sweep_forward
            enddo
        ! endif

        if (verbose_apply) then
            t_print = MPI_Wtime()-t_loop
            call MPI_ALLREDUCE(MPI_IN_PLACE, t_print, 1, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)
            call append_t_file('multigrid-iteration.t', (/1.0_rk, dble(i_level), t_print/))
        endif

        call toc( "GS Upwards - full sweep", 10005, MPI_Wtime()-t_loop )
    enddo

end subroutine multigrid_upwards



!> \brief Solve for the solution on the coarsest level
!> \details Solves for the solution on the coarsest level
!> It solves this with one of the following three modes:
!>   - FFT, but only for level 0
!>   - Conjugent gradient method, (but only for level 0 for now)
!>   - GS sweeps
subroutine multigrid_coarsest(params, hvy_sol, hvy_RHS, tree_ID, i_level, Jmax_a, verbose)
    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    real(kind=rk), intent(inout)       :: hvy_sol(:, :, :, :, :), hvy_RHS(:, :, :, :, :)
    integer(kind=ik), intent(in)       :: tree_ID
    integer(kind=ik), intent(in)       :: i_level
    integer(kind=ik), intent(in)       :: Jmax_a
    logical, intent(in), optional      :: verbose

    character(len=cshort)              :: fname
    integer(kind=tsize)                :: treecode

    logical            :: sweep_forward
    integer(kind=ik)   :: k_block, lgt_ID, hvy_id, nc, i_sweep, mpierr
    real(kind=rk)      :: dx(1:3), x0(1:3)
    real(kind=rk)      :: t_block, t_print(1:1)

    ! should be set in params
    real(kind=rk)      :: tol_cg
    integer(kind=ik)   :: it_sweep

    logical            :: verbose_apply

    verbose_apply = .false.
    if (present(verbose)) verbose_apply = verbose

    tol_cg = 1.0e-4
    it_sweep = 100

    nc = size(hvy_sol,4)

    if (i_level == 0 .and. (params%poisson_coarsest == "FFT" .or. params%poisson_coarsest == "CG")) then
        if (params%poisson_coarsest == "FFT") then
            ! FFT for lowest block if level = 0
            if (params%rank == 0 .and. verbose_apply) write(*, '(A, A, i0)') repeat('  ', i_level+1), 'FFT lvl ', i_level
            t_block = MPI_Wtime()

            ! we want to solve laplacian for finest level, so if FFT has fft precision, we need to adjust it to the finest level
            ! so let's get a dx on the finest level, we could loop over blocks and use get_block_origin_spacing, but maybe the processor with block on Jmin does not have one!
            ! so we aquire it directly by calling the underlying function with TC 0 and the correct level
            call get_block_spacing_origin_b( 0_tsize, params%domain_size, params%Bs, x0, dx, dim=params%dim, level=Jmax_a, max_level=params%Jmax)

            do k_block = 1, hvy_n(tree_ID)
                hvy_id = hvy_active(k_block, tree_ID)
                call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

                ! skip blocks not on this level
                if (lgt_block(lgt_id,IDX_MESH_LVL) /= i_level) cycle

                ! ! get spacing of lowest level block - usually we need that of level max though
                ! call get_block_spacing_origin( params, lgt_ID, x0, dx )

                call fft_solve_poisson(params, hvy_sol(:,:,:,1:nc,hvy_id), hvy_RHS(:,:,:,1:nc,hvy_id), dx)

                ! call dump_block_fancy(hvy_sol(:,:,:,1:nc,hvy_id), 'sol_L0.txt', params%Bs, params%g, digits=2, print_ghosts=.true.)
                ! call dump_block_fancy(hvy_RHS(:,:,:,1:nc,hvy_id), 'RHS_L0.txt', params%Bs, params%g, digits=2, print_ghosts=.true.)

            enddo
            if (verbose_apply) then
                t_print = MPI_Wtime()-t_block
                call MPI_ALLREDUCE(MPI_IN_PLACE, t_print, 1, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)
                call append_t_file('multigrid-iteration.t', (/0.0_rk, dble(i_level), t_print/))
            endif
            call toc( "fft solve poisson", 10004, MPI_Wtime()-t_block )
        elseif (params%poisson_coarsest == "CG") then
            ! conjugent gradient method for lowest block if level = 0
            if (params%rank == 0 .and. verbose_apply) write(*, '(A, 2(A, i0), A, es7.1)') repeat('  ', i_level+1), 'Coarsest CG lvl ', i_level, ' it max ', it_sweep, ' tol ', tol_cg
            t_block = MPI_Wtime()
            do k_block = 1, hvy_n(tree_ID)
                hvy_id = hvy_active(k_block, tree_ID)
                call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

                ! skip blocks not on this level
                if (lgt_block(lgt_id,IDX_MESH_LVL) /= i_level) cycle

                ! get spacing
                call get_block_spacing_origin( params, lgt_ID, x0, dx )
                
                call CG_solve_poisson_level0(params, hvy_sol(:,:,:,1:nc,hvy_id), hvy_RHS(:,:,:,1:nc,hvy_id), dx, tol_cg, it_sweep)
            enddo
            if (verbose_apply) then
                t_print = MPI_Wtime()-t_block
                call MPI_ALLREDUCE(MPI_IN_PLACE, t_print, 1, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)
                call append_t_file('multigrid-iteration.t', (/0.0_rk, dble(i_level), t_print/))
            endif
            call toc( "conjugent gradient", 10004, MPI_Wtime()-t_block )
        endif
    elseif (params%poisson_coarsest == "GS") then
        ! coarsest level sweeps
        if (params%rank == 0 .and. verbose_apply) write(*, '(A, 2(A, i0))') repeat('  ', i_level+1), 'Coarsest GS sweeps lvl ', i_level, ' it ', it_sweep

        sweep_forward = .true.
        do i_sweep = 1, it_sweep
            if (i_sweep == 1) then
                ! init solution as 0
                hvy_sol(:,:,:,1:nc,hvy_id) = 0.0_rk
            else
                ! we need to synch before a sweep
                t_block = MPI_Wtime()
                call sync_level_from_M( params, hvy_sol(:,:,:,1:nc,:), tree_ID, i_level, params%poisson_stencil_size, params%poisson_stencil_size)
                call toc( "Sync Layer", 10010, MPI_Wtime()-t_block )
            endif

            ! blocks on this level do a GS-sweep
            call GS_iteration_level(params, tree_id, i_level, hvy_sol(:,:,:,1:nc,:), hvy_RHS(:,:,:,1:nc,:), sweep_forward)

            sweep_forward = .not. sweep_forward
        enddo
    else
        write(fname, '(i0)') i_level
        call abort(250617, 'Coarsest level solver not implemented for : '//trim(params%poisson_coarsest)//' at level '//trim(fname))
    endif

end subroutine multigrid_coarsest