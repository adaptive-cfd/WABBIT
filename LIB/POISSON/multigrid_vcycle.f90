subroutine multigrid_solve(params, hvy_sol, hvy_RHS, hvy_work, tree_ID, init_0, verbose)
    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    real(kind=rk), intent(inout)       :: hvy_sol(:, :, :, :, :), hvy_RHS(:, :, :, :, :), hvy_work(:, :, :, :, :)
    integer(kind=ik), intent(in)       :: tree_ID
    logical, intent(in), optional      :: init_0  !< Initialize solution to zero
    logical, intent(in), optional      :: verbose

    integer(kind=ik)                   :: k_block, lgt_ID, hvy_id, ic, i_cycle
    real(kind=rk)                      :: dx(1:3), x0(1:3), residual(1:4*size(hvy_sol,4)), norm_sol(1:size(hvy_sol,4))
    real(kind=rk)                      :: t_block
    logical                            :: verbose_apply, init0
    character(len=cshort)              :: fname

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
    call init_full_tree(params, tree_ID, Jmin_set=params%poisson_Jmin, set_ref=0)

    ! init solution on lower levels as zero, leaf is that from input array (reuses last time-step for example) or also set to zero
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

        ! if (params%rank == 0 .and. .not. verbose_apply) write(*, '(A, i0, A, i0, A, 10(es10.3, 1x))') "   it ", i_cycle, "/", params%poisson_cycle_it, " Residual Linfty: ", residual(1:size(hvy_sol,4))

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

    if (params%rank == 0 .and. .not. verbose_apply) then
        write(fname, '(A, i0, A, i0, A)') "(A, i0, A, ", size(hvy_sol,4), "(es10.3, 1x), A, ", size(hvy_sol,4), "(es10.3, 1x))"
        write(*, fname) "   Final Residual after ", i_cycle, " it, Linfty: ", residual(1:size(hvy_sol,4)), " , rel to RHS: ", residual(1:size(hvy_sol,4))/norm_sol(1:size(hvy_sol,4))
    endif

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

    integer(kind=ik)                   :: k_block, lgt_ID, hvy_id, old_hvy_id, n_leaves
    integer(kind=ik)                   :: i_level, Jmax_a, it_coarsest, i_sweep, nc, ic, mpierr
    logical                            :: sweep_forward, tc_found

    character(len=cshort)              :: fname
    real(kind=rk)                      :: dx(1:3), x0(1:3), domain(1:3), tol_cg
    integer(kind=2)                    :: n_domain(1:3)
    integer(kind=tsize)                :: treecode
    logical                            :: verbose_apply
    
    ! Array for treecode mapping: (hvy_id, tc_part1, tc_part2)
    integer(kind=ik), save, allocatable :: tc_map_old(:,:)

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
    ! call toc( "MG Downwards - Store b", 10020, MPI_Wtime()-t_block )


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
    call toc( "MG Downwards - Compute residual", 10021, MPI_Wtime()-t_block )

    ! we use the decompose function to pass down the residual to all levels
    ! this overwrites b of the leaf layer, but we will recover it later
    ! solution in non-decomposed form is present in hvy_RHS
    t_block = MPI_Wtime()
    call wavelet_decompose_full_tree(params, hvy_work(:,:,:,1:nc,:), tree_ID, hvy_RHS, Jmin_set=params%poisson_Jmin, init_full_tree_grid=.false., compute_SC_only=.true., scaling_filter=params%MGR)
    call toc( "MG Downwards - full sweep", 10002, MPI_Wtime()-t_block )
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
    
    ! ========================================================================
    ! Build treecode mapping to handle block reordering after loadbalancing
    ! ========================================================================
    ! Problem: We backup the solution to hvy_work(:,:,:,:,hvy_id) before calling multigrid_upwards.
    !          During upwards, balanceLoad_tree is giving the same deterministic partitioning for the ranks, but the hvy_id of each block may change due to reordering.
    !          Only hvy_sol and hvy_RHS are sorted, so we need to find a way to map hvyt_work to avoid moving 3 arrays for loadBalancing
    !
    ! Solution: Before backup, store (hvy_id, treecode) pairs for all leaf blocks in tc_map_old, sorted by treecode.
    !          After loadbalancing, we search this map by treecode to find the original hvy_id where the backed-up data is stored,
    !          then access hvy_work(:,:,:,:,old_hvy_id) with the current block's new hvy_id.
    !
    ! tc_map_old structure: (hvy_id_old, treecode_part1, treecode_part2) sorted by treecode, similar to how it's done in balanceLoad_tree
    ! ========================================================================
    
    if (params%poisson_balanceLoad) then
        ! Allocate tc_map_old if not yet done (only needed when load balancing is enabled)
        if (.not. allocated(tc_map_old)) then
            allocate(tc_map_old(1:3, 1:params%number_blocks))
        endif
        tc_map_old = -1  ! we have to reset the map to -1 not 0 as 0 is a valid TC
        n_leaves = 0
        
        ! Build treecode mapping before loadbalancing (only for leaf blocks)
        do k_block = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k_block, tree_ID)
            call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

            ! only leaf layer has solution
            if (.not. block_is_leaf(params, hvy_id)) cycle
            
            n_leaves = n_leaves + 1
            tc_map_old(1, n_leaves) = hvy_id  ! store hvy_id
            treecode = get_tc(lgt_block(lgt_id, IDX_TC_1 : IDX_TC_2))
            call set_tc(tc_map_old(2:3, n_leaves), treecode)  ! store treecode
            
            ! backup solution
            hvy_work(:,:,:,1:nc,hvy_id) = hvy_sol(:,:,:,1:nc,hvy_id)

        enddo
        
        ! Sort tc_map_old by treecode for binary search later
        if (n_leaves > 1) then
            call quicksort(tc_map_old, 1, n_leaves, 3)
        endif
    else
        ! No load balancing - simple backup without treecode mapping
        n_leaves = 0
        do k_block = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k_block, tree_ID)
            if (.not. block_is_leaf(params, hvy_id)) cycle
            hvy_work(:,:,:,1:nc,hvy_id) = hvy_sol(:,:,:,1:nc,hvy_id)
        enddo
    endif
    
    call toc( "MG Downwards - Backup solution", 10022, MPI_Wtime()-t_block )

    ! lowest block is treated differently, here we try to solve more thoroughly. Following options exist:
    !    - Conjugent gradient method : only for level 0 for now
    !    - GS sweeps : very inefficient as blocks are usually too large for convergence
    !    - FFT :only for level 0, needs to be checked to solve true solution or discretized laplacian
    call multigrid_coarsest(params, hvy_sol, hvy_RHS, tree_ID, params%poisson_Jmin, Jmax_a, verbose_apply)

    ! upwards sweeps - prolong lower level solution, then solve Gauss-Seidel iterations
    ! do i_level = params%poisson_Jmin+1, Jmax_a
    !     call multigrid_upwards(params, hvy_sol, hvy_RHS, hvy_work, tree_ID, i_level, Jmax_a, hvy_depth)
    ! enddo
    call multigrid_upwards(params, hvy_sol, hvy_RHS, hvy_work, tree_ID, params%poisson_Jmin, Jmax_a, tc_map_old, n_leaves, verbose_apply)

    ! We do not solve Ae=r but Au=b and we need to restore the old solution as well as the RHS b
    ! this is happening for all leaf-blocks
    t_block = MPI_Wtime()
    do k_block = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k_block, tree_ID)
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

        if (block_is_leaf(params, hvy_id)) then
            ! get spacing
            call get_block_spacing_origin( params, lgt_id, x0, dx )
            
            ! After reordering, hvy_work is aligned with current hvy_id, so we can directly index it
            
            ! recompute RHS b = r + Ax, using current hvy_id (data was reordered)
            call GS_compute_residual(params, hvy_work(:,:,:,1:nc,hvy_id), hvy_RHS(:,:,:,1:nc,hvy_id), hvy_RHS(:,:,:,1:nc,hvy_id), dx, recompute_b=.true.)
            ! ! restore full RHS b
            ! hvy_RHS(:,:,:,1:nc,hvy_id) = hvy_work(:,:,:,nc+1:2*nc,hvy_id)
    
            ! add previous solution to reconstructed solution, using current hvy_id (data was reordered)
            hvy_sol(:,:,:,1:nc,hvy_id) = hvy_sol(:,:,:,1:nc,hvy_id) + hvy_work(:,:,:,1:nc,hvy_id)

        endif
    enddo

    ! sync RHS for sync_freq > 1
    t_block = MPI_Wtime()
    call sync_ghosts_tree(params, hvy_RHS(:,:,:,1:nc,:), tree_ID)
    call toc( "Sync Layer", 10010, MPI_Wtime()-t_block )
    call toc( "MG - RHS and solution restoration", 10032, MPI_Wtime()-t_block )

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
        if (params%rank == 0) write(*, '(A, es10.3, A)') "--- Residual L2: ", residual_out(nc+1), " ---"
        call componentWiseNorm_tree(params, hvy_work(:,:,:,1:nc,:), tree_ID, "L1", residual_out(2*nc+1:3*nc), threshold_state_vector=.false.)
        ! if (params%rank == 0) write(*, '(A, es10.3, A)') "--- Residual L1: ", residual_out(2*nc+1), " ---"

        if (params%rank == 0) write(*, '(A, es10.3, A)') "--- Residual Linfty: ", residual_out(1), " ---"


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
subroutine multigrid_upwards(params, hvy_sol, hvy_RHS, hvy_work, tree_ID, Jmin, Jmax_a, tc_map_old, n_leaves, verbose)
    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    real(kind=rk), intent(inout)       :: hvy_sol(:, :, :, :, :), hvy_RHS(:, :, :, :, :), hvy_work(:, :, :, :, :)
    integer(kind=ik), intent(in)       :: tree_ID
    integer(kind=ik), intent(in)       :: Jmin
    integer(kind=ik), intent(in)       :: Jmax_a
    integer(kind=ik), intent(inout)    :: tc_map_old(:, :)
    integer(kind=ik), intent(in)       :: n_leaves
    logical, intent(in), optional      :: verbose

    character(len=cshort)              :: fname
    integer(kind=tsize)                :: treecode
    real(kind=rk)                      :: norm(1:2*size(hvy_sol,4))

    logical            :: sweep_forward
    integer(kind=ik)   :: k_block, lgt_ID, hvy_id, nc, ic, i_sweep, mpierr, i_level, level_me, ref_me, sync_freq, i_g(1:3), i_BS(1:3)
    real(kind=rk)      :: t_loop, t_block, t_print(1:1)
    real(kind=rk)      :: dx(1:3), x0(1:3)
    logical            :: verbose_apply, balance_load_needed
    integer (kind=ik)  :: num_leaf, num_operate
    integer(kind=ik), dimension(:), allocatable, save :: lgt_refinementStatus_backup

    verbose_apply = .false.
    if (present(verbose)) verbose_apply = verbose
    balance_load_needed = .false.

    nc = size(hvy_sol,4)

    ! the refinement_status is used in this routine for things other than the refinement_status
    ! but that may be problematic if we have set the refinement_status externally and like to keep it.
    ! This happens for example in postprocessing (--grid1-to-grid2). In this case, the indicator is
    ! "nothing (external)" and we should not overwrite it. Hence, we make a backup, and copy the original
    ! status (when entering this routine) back to lgt_block
    if (.not. allocated(lgt_refinementStatus_backup)) then
        allocate(lgt_refinementStatus_backup(1:size(lgt_block,1)))
    endif
    lgt_refinementStatus_backup(:) = 0  ! reset back-up to be sure

    i_g = 0
    i_g(1:params%dim) = mod(params%g, 2)
    i_BS = 0
    i_BS(1:params%dim) = mod(params%Bs(1:params%dim), 2)

    ! prepare grid, set all ref stats to EMPTY and only lowest level to 1
    do k_block = 1, lgt_n(tree_ID)
        lgt_ID = lgt_active(k_block, tree_ID)
        level_me = lgt_block( lgt_ID, IDX_MESH_LVL )
        lgt_refinementStatus_backup(lgt_active(k_block, tree_ID)) = lgt_block(lgt_active(k_block, tree_ID), IDX_REFINE_STS)
        if (level_me == Jmin) then
            lgt_block(lgt_ID, IDX_REFINE_STS) = 1
        else
            lgt_block(lgt_ID, IDX_REFINE_STS) = REF_TMP_EMPTY
        endif
    end do

    ! count amount of leaf blocks for load balancing criterion
    num_leaf = 0
    do k_block = 1, hvy_n(tree_ID)
        if (block_is_leaf(params, hvy_active(k_block, tree_ID))) num_leaf = num_leaf + 1
    enddo

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
        call toc( "MG Upwards - Level definition", 10030, MPI_Wtime()-t_block )

        ! do upwards sync of solution field, this then needs to be interpolated
        t_block = MPI_Wtime()
        call sync_M2D(params, hvy_sol(:,:,:,1:nc,:), tree_ID, sync_case="ref", s_ref=-1)
        call toc( "Sync M2D", 10012, MPI_Wtime()-t_block )

        ! reordering of refinement flags:
        !   -1 -> REF_TMP_EMPTY : these blocks have synched and are finished
        do k_block = 1, lgt_n(tree_ID)
            lgt_id = lgt_active(k_block, tree_ID)
            if (lgt_block(lgt_id,IDX_REFINE_STS) == -1) then
                lgt_block(lgt_id,IDX_REFINE_STS) = REF_TMP_EMPTY
            endif
        enddo

        if (params%poisson_balanceLoad) then
            ! we do not always want to do loadbalancing. For example, on the first level, we operate on a maximum of 8 blocks. Even if this is super imbalanced, then the overall load is still very small and we pay more in communication than we gain in balancing.
            ! For a criterion, if any processor exceeds 1/3 of the amount of blocks on leaf layer, we do balancing
            ! Once the threshold is exceeded, we always do balancing on higher levels as well
            ! If no level below the final leaf-layer needs balancing, then we skip balancing on the final leaf-layer as well as this stayed balanced
            if (.not. balance_load_needed .and. i_level /= Jmax_a) then
                num_operate = 0
                do k_block = 1, hvy_n(tree_ID)
                    hvy_id = hvy_active(k_block, tree_ID)
                    call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
                    ! count blocks that will be operated on this level
                    if (lgt_block(lgt_id, IDX_REFINE_STS) == 0 .or. lgt_block(lgt_id, IDX_REFINE_STS) == 1) then
                        num_operate = num_operate + 1
                    endif
                    ! exit criteria
                    if (num_operate > num_leaf/3 .and. num_leaf > 8) then
                        balance_load_needed = .true.
                        exit
                    endif
                enddo
            endif
            ! num_leaf is not the same for all procs, this is why the upper loop has to be local and we need to communicate the result
            call MPI_ALLREDUCE(MPI_IN_PLACE, balance_load_needed, 1, MPI_LOGICAL, MPI_LOR, WABBIT_COMM, mpierr)

            if (balance_load_needed) then

                ! balance load for blocks with ref status 0 and 1 (leaf blocks and new blocks on current level)
                ! transfer both hvy_sol and hvy_RHS together
                write(fname, '(A, i0)') "MultiGrid_upwards_L", i_level
                t_block = MPI_Wtime()
                call balanceLoad_tree(params, hvy_sol(:,:,:,1:nc,:), tree_ID, balanceMode="selective", &
                                    balance_ref=(/ 0_ik, 1_ik /), balance_name=fname, hvy_tmp=hvy_RHS(:,:,:,1:nc,:), Jmin_set=params%poisson_Jmin, full_tree_grid=.true.)
                call toc( "MG Upwards - balanceLoad_tree", 10050, MPI_Wtime()-t_block )

                ! Reorder data arrays after load balancing at the finest level
                ! This ensures consistent ordering across processors and eliminates the need for treecode lookup
                if (i_level == Jmax_a) then
                    t_block = MPI_Wtime()
                    call reorder_hvy_arrays(params, hvy_sol, hvy_RHS, tree_ID, tc_map_old, n_leaves)
                    call toc( "MG upwards - Reorder hvy arrays after balancing", 10051, MPI_Wtime()-t_block )
                endif

            endif
        endif

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
        call toc( "MG Upwards - Solution prolongation", 10031, MPI_Wtime()-t_block )

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
                call GS_iteration_ref(params, tree_id, (/ 1, 0 /), hvy_sol(:,:,:,1:nc,:), hvy_RHS(:,:,:,1:nc,:), sweep_forward, filter_offset=params%g, sweep_number=i_sweep, multigrid_level=i_level)
                ! call GS_iteration_ref(params, tree_id, (/ 1, 0 /), hvy_sol(:,:,:,1:nc,:), hvy_RHS(:,:,:,1:nc,:), sweep_forward, filter_offset=max(0,params%g-(sync_freq-1)*params%poisson_stencil_size))

                sweep_forward = .not. sweep_forward
            enddo
        ! endif

        if (verbose_apply) then
            t_print = MPI_Wtime()-t_loop
            call MPI_ALLREDUCE(MPI_IN_PLACE, t_print, 1, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)
            call append_t_file('multigrid-iteration.t', (/1.0_rk, dble(i_level), t_print/))
        endif

        call toc( "MG Upwards - full sweep", 10005, MPI_Wtime()-t_loop )
    enddo

    ! copy the original refinement_status (when entering this routine) back to lgt_block
    do k_block = 1, lgt_n(tree_ID)
        lgt_block(lgt_active(k_block, tree_ID), IDX_REFINE_STS) = lgt_refinementStatus_backup(lgt_active(k_block, tree_ID))
    enddo

end subroutine multigrid_upwards


!> \brief Reorder heavy data arrays to restore original ordering after load balancing
!> \details Uses cycle-following algorithm to permute data in-place with minimal memory overhead.
!! After load balancing, hvy_sol and hvy_RHS are reordered but hvy_work (backup) is not.
!! This routine reorders hvy_sol and hvy_RHS back to the original ordering (matching hvy_work).
!! 
!! Ordering terminology:
!! - Original ordering: where hvy_work has backed up data (encoded in tc_map_old before load balance)
!! - Current ordering: where hvy_sol/hvy_RHS are after load balancing
!! - Desired ordering: same as original ordering (goal of this routine)
subroutine reorder_hvy_arrays(params, hvy_sol, hvy_RHS, tree_ID, tc_map_old, n_leaves)
    implicit none
    
    type (type_params), intent(inout)  :: params
    real(kind=rk), intent(inout)       :: hvy_sol(:, :, :, :, :), hvy_RHS(:, :, :, :, :)
    integer(kind=ik), intent(in)       :: tree_ID
    integer(kind=ik), intent(inout)    :: tc_map_old(:, :)
    integer(kind=ik), intent(in)       :: n_leaves
    
    integer(kind=ik) :: k_block, hvy_id, lgt_id, original_hvy_id, next_lgt_id, i, j, start_idx, current_idx, next_idx, last_unassigned, nc
    integer(kind=tsize) :: treecode
    logical :: tc_found
    logical, allocatable, save :: visited(:)
    real(kind=rk), allocatable, save :: tmp_block_sol(:,:,:,:), tmp_block_rhs(:,:,:,:)
    integer(kind=ik), allocatable, save :: current_to_original(:), original_to_current(:), tmp_lgt_block(:)
    integer(kind=ik) :: Bs(1:3), g
    
    ! Get grid parameters
    Bs = params%Bs
    g = params%g
    nc = size(hvy_sol, 4)
    
    ! ============================================================================
    ! Build mapping between original and current positions
    ! ============================================================================
    ! We build two arrays:
    ! 1. current_to_original(current_hvy_id) = original_hvy_id
    !    "Block at current position came from which original position"
    ! 2. original_to_current(original_hvy_id) = current_hvy_id
    !    "Data from original position is now at which current position"
    !    This is the inverse of (1) and is used for cycle-following
    !
    ! For unassigned positions (non-leaf or empty blocks):
    !   - current_to_original = 0 AND original_to_current = 0: truly unused
    !   - current_to_original = 0 AND original_to_current ≠ 0: waypoint position
    !   - current_to_original ≠ 0 AND original_to_current = 0: invalid state
    ! ============================================================================
    if (.not. allocated(current_to_original)) allocate(current_to_original(params%number_blocks))
    if (.not. allocated(original_to_current)) allocate(original_to_current(params%number_blocks))
    if (.not. allocated(visited)) allocate(visited(params%number_blocks))
    
    current_to_original = 0
    original_to_current = 0
    visited = .false.
    last_unassigned = 1  ! Optimization: track last used unassigned position
    
    ! Build mapping for leaf blocks only
    do k_block = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k_block, tree_ID)  ! Current position after load balance
        call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
        
        ! Only process leaf blocks
        if (.not. block_is_leaf(params, hvy_id)) cycle
        
        ! Get treecode of current leaf block
        treecode = get_tc(lgt_block(lgt_id, IDX_TC_1 : IDX_TC_2))
        
        ! Find this treecode in tc_map_old to get original hvy_id (before load balance)
        call find_id_by_treecode(tc_map_old, n_leaves, treecode, original_hvy_id, tc_found)
        
        if (.not. tc_found) then
            call abort(260116, "[reorder_hvy_arrays] ERROR: Treecode not found in tc_map_old during reordering")
        endif
        
        ! Build both mappings
        current_to_original(hvy_id) = original_hvy_id
        original_to_current(original_hvy_id) = hvy_id
    enddo

    ! Allocate temporary storage for one block only
    if (allocated(tmp_block_sol)) then
        if (size(tmp_block_sol, 4) /= nc) deallocate(tmp_block_sol)
    endif
    if (.not. allocated(tmp_block_sol)) allocate(tmp_block_sol(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g, nc))
    if (allocated(tmp_block_rhs)) then
        if (size(tmp_block_rhs, 4) /= nc) deallocate(tmp_block_rhs)
    endif
    if (.not. allocated(tmp_block_rhs)) allocate(tmp_block_rhs(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g, nc))
    if (.not. allocated(tmp_lgt_block)) allocate(tmp_lgt_block(size(lgt_block, 2)))

    ! ============================================================================
    ! Cycle-following algorithm for in-place data reordering
    ! ============================================================================
    ! Goal: Move data from current positions back to original (desired) positions
    !       to align hvy_sol/hvy_RHS with hvy_work
    !
    ! Strategy: Use original_to_current mapping to follow permutation cycles.
    !   For each original position:
    !     - Find where its data currently is (original_to_current)
    !     - Pull that data to the desired position
    !     - Continue following the cycle
    !
    ! Example cycle:
    !   Original pos 2 wants data currently at pos 5 (original_to_current(2)=5)
    !   Original pos 5 wants data currently at pos 3 (original_to_current(5)=3)
    !   Original pos 3 wants data currently at pos 2 (original_to_current(3)=2)
    !   This forms cycle: 2→5→3→2
    !
    ! Algorithm:
    !   1. Start at unvisited original position, save its current data
    !   2. Follow chain: next_pos = original_to_current(current_original_pos)
    !   3. Pull data from next_pos to current position
    !   4. If we return to start: close cycle, write saved data
    !   5. If original_to_current(pos) = 0: position is unassigned
    !      - Search for an unassigned waypoint position to continue cycle
    !      - Waypoint criteria: current_to_original=0 AND original_to_current≠0
    !      - Handle inactive blocks (lgt_block < 0) by setting current to -1
    !      - Mark truly unused positions (both mappings = 0) as visited
    !
    ! Special cases:
    !   - Unassigned start positions: not marked visited until cycle completes
    !   - Infinite loop protection: bounded by params%number_blocks iterations
    !
    ! Minimal memory: Only 1 temporary block needed for the start of each cycle
    ! ============================================================================
    
    ! Apply cycle-following with original_to_current mapping
    do k_block = 1, hvy_n(tree_ID)
        start_idx = hvy_active(k_block, tree_ID)
        
        ! Skip conditions:
        if (visited(start_idx)) cycle  ! Already processed
        if (original_to_current(start_idx) == 0 .and. current_to_original(start_idx) == 0) cycle  ! Truly unused, stays in place
        if (original_to_current(start_idx) == start_idx) then  ! Identity mapping, no swap needed
            visited(start_idx) = .true.
            cycle
        endif
        
        ! Start a new cycle/chain at this original position
        current_idx = start_idx  ! Current original position being filled
        next_idx = start_idx
        call hvy2lgt(lgt_id, current_idx, params%rank, params%number_blocks)
        next_lgt_id = lgt_id
        
        ! Save the data currently at the start position (will be overwritten)
        tmp_block_sol = hvy_sol(:,:,:,1:nc,current_idx)
        tmp_block_rhs = hvy_RHS(:,:,:,1:nc,current_idx)
        tmp_lgt_block = lgt_block(lgt_id, :)
        
        ! Follow the cycle: pull data from current locations to original positions
        ! Bounded loop to prevent infinite cycles (safety measure)
        do j = 1, params%number_blocks
            ! Update to next position in cycle
            current_idx = next_idx  ! Original position to fill now
            lgt_id = next_lgt_id
            next_idx = original_to_current(current_idx)  ! Find current location of data for this original pos
            call hvy2lgt(next_lgt_id, next_idx, params%rank, params%number_blocks)

            ! ========================================================================
            ! Handle unassigned positions (original_to_current = 0)
            ! ========================================================================
            ! When the original position has no assigned data (non-leaf or empty),
            ! we need to find a waypoint position to continue the cycle
            if (next_idx == 0) then
                ! Search for an unassigned position to use as temporary storage (waypoint)
                ! Waypoint criteria: current_to_original=0 (not source) AND original_to_current≠0 (is target)
                do i = last_unassigned, params%number_blocks
                    next_idx = i
                    call hvy2lgt(next_lgt_id, next_idx, params%rank, params%number_blocks)
                    
                    ! Check if this position qualifies as a waypoint
                    if (.not. visited(next_idx) .and. current_to_original(next_idx) == 0 .and. original_to_current(next_idx) /= 0) then
                        ! Valid waypoint found
                        if (lgt_block(next_lgt_id, IDX_TC_1) >= 0) then
                            ! Active block: pull data from waypoint
                            if (next_idx == start_idx) then
                                ! Special case: looped back to start, write saved data and finish
                                hvy_sol(:,:,:,1:nc,current_idx) = tmp_block_sol
                                hvy_RHS(:,:,:,1:nc,current_idx) = tmp_block_rhs
                                lgt_block(lgt_id, :) = tmp_lgt_block
                            else
                                ! Normal case: pull from waypoint to current
                                hvy_sol(:,:,:,1:nc,current_idx) = hvy_sol(:,:,:,1:nc,next_idx)
                                hvy_RHS(:,:,:,1:nc,current_idx) = hvy_RHS(:,:,:,1:nc,next_idx)
                                lgt_block(lgt_id, :) = lgt_block(next_lgt_id, :)
                            endif
                        else
                            ! Inactive block: current position becomes empty (propagate empty slot)
                            lgt_block(lgt_id, :) = -1
                        endif
                        visited(current_idx) = .true.
                        last_unassigned = next_idx + 1  ! Optimize next search
                        exit  ! Continue cycle from this waypoint
                    elseif ((current_to_original(next_idx) == 0 .and. original_to_current(next_idx) == 0) .or. lgt_block(next_lgt_id, IDX_TC_1) < 0) then
                        ! Mark truly unused positions as visited to prevent reuse
                        visited(next_idx) = .true.
                    endif
                enddo

                ! Check if we've returned to start (cycle complete via waypoints)
                if (next_idx == start_idx) exit

                ! Safety check: if no waypoint found, abort
                if (next_idx == params%number_blocks .and. .not. (current_to_original(next_idx) == 0 .and. original_to_current(next_idx) /= 0)) then
                    call abort(260115, "[reorder_hvy_arrays] ERROR: Unable to find unvisited unassigned waypoint position")
                endif

            elseif (next_idx == start_idx) then
                ! ====================================================================
                ! Cycle complete - returned to start
                ! ====================================================================
                ! Write the saved data from the start of the cycle to current position
                hvy_sol(:,:,:,1:nc,current_idx) = tmp_block_sol
                hvy_RHS(:,:,:,1:nc,current_idx) = tmp_block_rhs
                lgt_block(lgt_id, :) = tmp_lgt_block
                visited(current_idx) = .true.
                exit  ! Cycle complete
            else
                ! ====================================================================
                ! Continue following the cycle (normal case)
                ! ====================================================================
                ! Pull data from next_idx (current location) to current_idx (original/desired location)
                hvy_sol(:,:,:,1:nc,current_idx) = hvy_sol(:,:,:,1:nc,next_idx)
                hvy_RHS(:,:,:,1:nc,current_idx) = hvy_RHS(:,:,:,1:nc,next_idx)
                lgt_block(lgt_id, :) = lgt_block(next_lgt_id, :)
                
                ! Mark as visited, except for unassigned start positions (handled at cycle completion)
                if (.not. (current_to_original(current_idx) == 0 .and. current_idx == start_idx)) then
                    visited(current_idx) = .true.
                endif
            endif
        enddo
    enddo

    ! update grid metadata is necessary
    call synchronize_lgt_data(params, refinement_status_only=.false.)
    call updateMetadata_tree(params, tree_ID, search_overlapping=.true., Jmin_set=params%poisson_Jmin)
    
    ! Note: Arrays are declared with SAVE attribute, so they persist between calls
    ! and are not deallocated here (only reallocated if nc changes)
    
end subroutine reorder_hvy_arrays


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
            if (params%rank == 0 .and. verbose_apply) write(*, '(A, 2(A, i0), A, es8.1)') repeat('  ', i_level+1), 'Coarsest CG lvl ', i_level, ' it max ', it_sweep, ' tol ', tol_cg
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