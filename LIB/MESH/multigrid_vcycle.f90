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
! u - upwards prolongation and Gauss-Seidel sweeps with it_up
! d - Decomposition with restriction filter
subroutine multigrid_vcycle(params, hvy_sol, hvy_RHS, hvy_work, tree_ID, a, stencil, it_down, it_up, fft_order_discretization)
    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    real(kind=rk), intent(inout)       :: hvy_sol(:, :, :, :, :), hvy_RHS(:, :, :, :, :), hvy_work(:, :, :, :, :)
    integer(kind=ik), intent(in)       :: tree_ID
    integer(kind=ik), intent(in)       :: a
    real(kind=rk), intent(in)          :: stencil(-a:a)  ! FD operator
    integer(kind=ik), intent(in)       :: it_down, it_up
    character(len=cshort), intent(in)  :: fft_order_discretization

    integer(kind=ik)                   :: k_block, lgt_ID, hvy_id
    integer(kind=ik)                   :: i_level, Jmin, Jmax_a, it_coarsest, i_sweep, nc, ic, mpierr
    logical                            :: sweep_forward

    character(len=cshort)              :: fname
    real(kind=rk)                      :: dx(1:3), x0(1:3), domain(1:3), tol_cg
    real(kind=rk), allocatable         :: norm(:)
    integer(kind=2)                    :: n_domain(1:3)
    integer(kind=tsize)                :: treecode

    real(kind=rk)        :: t_block, t_loop, t_cycle, t_print(1:1)

    nc = size(hvy_sol,4)
    allocate(norm(3*nc))

    ! Logic to define the solution layer at each iteration:
    !    Downwards - we simply do a wavelet decomposition of the residual, but only using the scaling filter
    !    Upwards   - we reconstruct similarly to reconstruction, but always keep all leaf-blocks of lower levels to keep a full periodic grid
    !                after each reconstruction we do the Gauss-Seidel sweeps/iterations

    ! ------------------------------------------------------------
    ! Gauss-Seidel Multi-Grid - solve the system Ax = b
    !    When going to lower levels, we compute the residual r = b - Ax and solve on lower levels Ae = r
    !    We decompose the residual until lowest level, solve using FFT, CG or GS there and then prolongate upwards and solve using GS
    Jmin = 0
    it_coarsest = -1  ! JMin=0 : if set to -1, FFT is used, elsewise CG, for JMin>0 GS sweeps are always used
    tol_cg = 1e-4

    t_cycle = MPI_Wtime()

    Jmax_a = maxActiveLevel_tree(tree_id)

    ! ! downwards sweeps - do Gauss-Seidel iterations, compute residual r = b - Ax and restrict residual downwards
    ! do i_level = Jmax_a, Jmin+1, -1
    !     call multigrid_downwards(params, hvy_sol, hvy_RHS, hvy_work, tree_ID, i_level, Jmax_a, a, stencil, it_down, hvy_depth)
    ! enddo

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
        call GS_compute_residual(params, hvy_sol(:,:,:,1:nc,hvy_id), hvy_RHS(:,:,:,1:nc,hvy_id), hvy_work(:,:,:,1:nc,hvy_id), a, stencil, dx)

    enddo
    call toc( "GS Downwards - Compute residual", 10021, MPI_Wtime()-t_block )

    ! we use the decompose function to pass down the residual to all levels
    ! this overwrites b of the leaf layer, but we will recover it later
    ! solution in non-decomposed form is present in hvy_RHS
    t_block = MPI_Wtime()
    call wavelet_decompose_full_tree(params, hvy_work(:,:,:,1:nc,:), tree_ID, hvy_RHS, init_full_tree_grid=.false., compute_SC_only=.true., scaling_filter=params%MGR)
    call toc( "GS Downwards - full sweep", 10002, MPI_Wtime()-t_block )
    t_block = MPI_Wtime()-t_block
    call MPI_ALLREDUCE(MPI_IN_PLACE, t_block, 1, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)
    call append_t_file('multigrid-iteration.t', (/-1.0_rk, -1.0_rk, t_block/))

    ! usually after decomposition the residual is not synced, when doing multi-step GS sweeps, we need to sync it
    t_block = MPI_Wtime()
    call sync_ghosts_tree(params, hvy_RHS(:,:,:,1:nc,:), tree_ID)
    call toc( "Sync Layer", 10010, MPI_Wtime()-t_block )

    ! for all levels we are not solving Au = b, but rather Ae=r
    ! only on the last upwards step will we iterate on the original problem, so we need to backup u
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
    call multigrid_coarsest(params, hvy_sol, hvy_RHS, tree_ID, Jmin, Jmax_a, a, stencil, it_coarsest, tol_cg, fft_order_discretization)

    ! upwards sweeps - prolong lower level solution, then solve Gauss-Seidel iterations
    ! do i_level = Jmin+1, Jmax_a
    !     call multigrid_upwards(params, hvy_sol, hvy_RHS, hvy_work, tree_ID, i_level, Jmax_a, a, stencil, it_up, hvy_depth)
    ! enddo
    call multigrid_upwards_2(params, hvy_sol, hvy_RHS, hvy_work, tree_ID, Jmin, Jmax_a, a, stencil, it_up)

    ! sync before computing final residual
    call sync_ghosts_tree(params, hvy_sol(:,:,:,1:nc,:), tree_ID, a, a)
    call toc( "Sync Layer", 10010, MPI_Wtime()-t_block )

    ! compute and output residual only on leaf layer
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
        call GS_compute_residual(params, hvy_sol(:,:,:,1:nc,hvy_id), hvy_RHS(:,:,:,1:nc,hvy_id), hvy_work(:,:,:,1:nc,hvy_id), a, stencil, dx)
    enddo
    call componentWiseNorm_tree(params, hvy_work(:,:,:,1:nc,:), tree_ID, "L2", norm)
    if (params%rank == 0) write(*, '(A, es10.4, A)') "--- Residual L2: ", norm(1), " ---"
    call componentWiseNorm_tree(params, hvy_work(:,:,:,1:nc,:), tree_ID, "L1", norm(nc+1:2*nc))
    ! if (params%rank == 0) write(*, '(A, es10.4, A)') "--- Residual L1: ", norm(nc+1), " ---"
    call componentWiseNorm_tree(params, hvy_work(:,:,:,1:nc,:), tree_ID, "Linfty", norm(2*nc+1:3*nc))
    if (params%rank == 0) write(*, '(A, es10.4, A)') "--- Residual Linfty: ", norm(nc+2), " ---"

    t_print = MPI_Wtime()-t_cycle
    call MPI_ALLREDUCE(MPI_IN_PLACE, t_print, 1, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)
    call append_t_file('multigrid-cycle.t', (/norm(1), norm(nc+1), norm(2*nc+1), t_print(1)/))

    call toc( "Final residual", 10014, MPI_Wtime()-t_block )

    call toc( "V cycle", 10001, MPI_Wtime()-t_cycle )

    deallocate(norm)

end subroutine multigrid_vcycle



! !> \brief Downwards level iteration for the multigrid solver
! !> \details Compute for one level a downwards step
! !> Does four things:
! !>   - initializes solution
! !>   - do GS sweeps on this level
! !>   - compute the residual r = b - Ax
! !>   - decompose and sync the solution to lower level
! subroutine multigrid_downwards(params, hvy_sol, hvy_RHS, hvy_work, tree_ID, i_level, Jmax_a, a, stencil, it_sweep, hvy_depth)
!     implicit none

!     !> parameter struct
!     type (type_params), intent(inout)  :: params
!     real(kind=rk), intent(inout)       :: hvy_sol(:, :, :, :, :), hvy_RHS(:, :, :, :, :), hvy_work(:, :, :, :, :)
!     integer(kind=ik), intent(in)       :: tree_ID
!     integer(kind=ik), intent(in)       :: i_level
!     integer(kind=ik), intent(in)       :: Jmax_a
!     integer(kind=ik), intent(in)       :: a
!     real(kind=rk), intent(in)          :: stencil(-a:a)
!     integer(kind=ik), intent(in)       :: it_sweep
!     integer(kind=ik), intent(inout)    :: hvy_depth(:)

!     character(len=cshort)              :: fname
!     integer(kind=tsize)                :: treecode
    
!     logical            :: sweep_forward
!     integer(kind=ik)   :: k_block, lgt_ID, hvy_id, nc, i_sweep, mpierr
!     real(kind=rk)      :: dx(1:3), x0(1:3)
!     real(kind=rk)      :: t_loop, t_block, t_print(1:1)
!     nc = size(hvy_sol,4)

!     t_loop = MPI_Wtime()
!     if (params%rank == 0) write(*, '(A, 2(A, i0))') repeat('  ', i_level+1), 'Downwards GS sweep lvl ', i_level, ' it ', it_sweep

!     ! init solution u or e as 0. u is only initialized the first time, while e will always be initialized as 0
!     ! why? Because every V-cycle we compute a new r_i and solve a new system Ae_i = r_i which is different to before
!     t_block = MPI_Wtime()
!     if (i_level /= Jmax_a) then
!         do k_block = 1, hvy_n(tree_ID)
!             hvy_id = hvy_active(k_block, tree_ID)
!             call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

!             ! skip blocks not on this subgrid, defined with Ref flag of -1
!             if (lgt_block(lgt_id,IDX_REFINE_STS) /= -1) cycle
!             hvy_sol(:,:,:,1:nc,hvy_id) = 0.0_rk
!         enddo
!     endif
!     call toc( "GS Downwards - Solution initialization", 10020, MPI_Wtime()-t_block )

!     ! do GS sweeps
!     t_block = MPI_Wtime()
!     do i_sweep = 1, it_sweep
!         ! we need to synch before a sweep - use normal sync as blocks with REF_TMP_EMPTY will not sync, leading in only this iterations layer syncing with itself
!         t_block = MPI_Wtime()
!         call sync_ghosts_tree(params, hvy_sol(:,:,:,1:nc,:), tree_ID, a, a, ignore_Filter=.true.)
!         call toc( "Sync Layer", 10010, MPI_Wtime()-t_block )

!         ! blocks on this subgrid do a GS-sweep
!         call GS_iteration_ref(params, tree_id, (/-1/), hvy_sol(:,:,:,1:nc,:), hvy_RHS(:,:,:,1:nc,:), a, stencil, sweep_forward)
!         sweep_forward = .not. sweep_forward
!     enddo

!     ! we finished all downwards GS sweeps on this level, now we compute the residual r = b - Ax and pass it downwards
!     ! if we do 0 GS-sweeps and are on lower iterations, the solution is 0, and therefore r = b
!     if (.not. (i_level /= Jmax_a .and. it_sweep == 0)) then
!         t_block = MPI_Wtime()
!         do k_block = 1, hvy_n(tree_ID)
!             hvy_id = hvy_active(k_block, tree_ID)
!             call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

!             ! skip blocks not on this subgrid, defined with Ref flag of -1
!             if (lgt_block(lgt_id,IDX_REFINE_STS) /= -1) cycle

!             ! get spacing
!             call get_block_spacing_origin( params, lgt_ID, x0, dx )

!             ! compute actual residual r = b - Ax, overwrite b (we will recover it later)
!             call GS_compute_residual(params, hvy_sol(:,:,:,1:nc,hvy_id), hvy_RHS(:,:,:,1:nc,hvy_id), hvy_RHS(:,:,:,1:nc,hvy_id), a, stencil, dx)

!         enddo
!         call toc( "GS Downwards - Compute residual", 10021, MPI_Wtime()-t_block )
!     endif

!     ! sync residual, so it can be decomposed
!     t_block = MPI_Wtime()
!     call sync_ghosts_tree(params, hvy_RHS(:,:,:,1:nc,:), tree_ID)
!     call toc( "Sync Layer", 10010, MPI_Wtime()-t_block )

!     ! preparation for next level:
!     !   - ref = -1: decompose residual and sync it downwards
!     !   - ref =  0: backup solution of this field by overwriting residual

!     ! now decide which blocks can go down one level
!     do k_block = 1, lgt_n(tree_ID)
!         lgt_id = lgt_active(k_block, tree_ID)
!         ! only blocks on specific level can go down on level, others need to stay
!         if (lgt_block(lgt_id,IDX_REFINE_STS) == -1 .and. lgt_block(lgt_id,IDX_MESH_LVL) /= i_level) then
!             lgt_block(lgt_id,IDX_REFINE_STS) = 0
!         endif
!     enddo

!     !     all: store this solution in hvy_work
!     ! ref = 1: decompose residual
!     ! ref = 0: increase depth counter
!     t_block = MPI_Wtime()
!     do k_block = 1, hvy_n(tree_ID)
!         hvy_id = hvy_active(k_block, tree_ID)
!         call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

!         ! blocks with -1 will decompose residual, others will increase depth counter
!         if (lgt_block(lgt_id,IDX_REFINE_STS) == -1) then
!             ! store r before decomposing so that we can restore it
!             hvy_work(:,:,:,hvy_depth(k_block)*nc+1:hvy_depth(k_block)*nc+nc,hvy_id) = hvy_RHS(:,:,:,1:nc,hvy_id)

!             ! decompose this residual, we only need SC
!             call blockFilterXYZ_vct( params, hvy_RHS(:,:,:,1:nc,hvy_ID), hvy_RHS(:,:,:,1:nc,hvy_ID), params%MGR, &
!             lbound(params%MGR, dim=1), ubound(params%MGR, dim=1), do_restriction=.true.)
!         elseif (lgt_block(lgt_id,IDX_REFINE_STS) == 0) then
!             ! backup solution field
!             hvy_work(:,:,:,hvy_depth(k_block)*nc+1:hvy_depth(k_block)*nc+nc,hvy_id) = hvy_sol(:,:,:,1:nc,hvy_id)
!             hvy_depth(k_block) = hvy_depth(k_block) + 1
!         endif

!         ! activate mother, so that they can receive the residual
!         ! if daughter block exists and any/all have ref=-1: this block will be on subgrid in next iteration, so we set it a new refinement value (0 for now)
!         if (.not. block_is_leaf(params, hvy_id)) then
!             if (any(lgt_block(hvy_family(hvy_ID, 2+2**params%dim:1+2**(params%dim+1)),IDX_REFINE_STS) == -1)) then
!                 lgt_block(lgt_id,IDX_REFINE_STS) = 0
!             endif
!         endif
!     enddo
!     call toc( "GS Downwards - Residual decomposition", 10022, MPI_Wtime()-t_block )

!     ! do downwards sync of decomposed residual r, ready to be used for the next sweeps, only blocks with ref = -1
!     t_block = MPI_Wtime()
!     call sync_D2M(params, hvy_RHS, tree_ID, sync_case="ref", s_val=-1)
!     call toc( "Sync D2M", 10011, MPI_Wtime()-t_block )

!     ! ref = 1: hvy_RHS was decomposed, we restore it and store the solution
!     ! ref = 0: solution was already stored
!     t_block = MPI_Wtime()
!     do k_block = 1, hvy_n(tree_ID)
!         hvy_id = hvy_active(k_block, tree_ID)
!         call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

!         ! restore residual r in RHS, store solution and wipe the buffer (so that we can sync SC into it)
!         if (lgt_block(lgt_id,IDX_REFINE_STS) == -1) then
!             hvy_RHS(:,:,:,1:nc,hvy_id) = hvy_work(:,:,:,hvy_depth(k_block)*nc+1:hvy_depth(k_block)*nc+nc,hvy_id)
!             hvy_work(:,:,:,hvy_depth(k_block)*nc+1:hvy_depth(k_block)*nc+nc,hvy_id) = hvy_sol(:,:,:,1:nc,hvy_id)
!             hvy_sol(:,:,:,1:nc,hvy_id) = 0.0_rk
!         endif
!     enddo
!     ! define next grid: blocks with ref=0 get ref=-1 and blocks with ref=-1 get ref=REF_TMP_EMPTY
!     do k_block = 1, hvy_n(tree_ID)
!         hvy_id = hvy_active(k_block, tree_ID)
!         call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
!         if (lgt_block(lgt_id,IDX_REFINE_STS) == -1) then
!             lgt_block(lgt_id,IDX_REFINE_STS) = REF_TMP_EMPTY
!         elseif (lgt_block(lgt_id,IDX_REFINE_STS) == 0) then
!             lgt_block(lgt_id,IDX_REFINE_STS) = -1
!         endif
!     enddo
!     call toc( "GS Downwards - Copy back", 10023, MPI_Wtime()-t_block )

!     ! update of refinement status needed afterwards
!     call synchronize_lgt_data( params, refinement_status_only=.true.)

!     t_print = MPI_Wtime()-t_loop
!     call MPI_ALLREDUCE(MPI_IN_PLACE, t_print, 1, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)
!     call append_t_file('multigrid-iteration.t', (/-1.0_rk, dble(i_level), t_print(1)/))

!     call toc( "GS Downwards - full sweep", 10002, MPI_Wtime()-t_loop )

! end subroutine multigrid_downwards


! !> \brief Upwards level iteration for the multigrid solver
! !> \details Compute for one level an upwards step
! !> Does three things:
! !>   - sync the residual from the lower lewel as WD decomposed values, reconstruct the residual
! !>   - add previous solution of this level
! !>   - do GS sweeps on this level
! subroutine multigrid_upwards(params, hvy_sol, hvy_RHS, hvy_work, tree_ID, i_level, Jmax_a, a, stencil, it_sweep, hvy_depth)
!     implicit none

!     !> parameter struct
!     type (type_params), intent(inout)  :: params
!     real(kind=rk), intent(inout)       :: hvy_sol(:, :, :, :, :), hvy_RHS(:, :, :, :, :), hvy_work(:, :, :, :, :)
!     integer(kind=ik), intent(in)       :: tree_ID
!     integer(kind=ik), intent(in)       :: i_level
!     integer(kind=ik), intent(in)       :: Jmax_a
!     integer(kind=ik), intent(in)       :: a
!     real(kind=rk), intent(in)          :: stencil(-a:a)
!     integer(kind=ik), intent(in)       :: it_sweep
!     integer(kind=ik), intent(inout)    :: hvy_depth(:)

!     character(len=cshort)              :: fname
!     integer(kind=tsize)                :: treecode

!     logical            :: sweep_forward
!     integer(kind=ik)   :: k_block, lgt_ID, hvy_id, nc, ic, i_sweep, mpierr
!     real(kind=rk)      :: t_loop, t_block, t_print(1:1)
!     real(kind=rk)      :: dx(1:3), x0(1:3)
!     nc = size(hvy_sol,4)

!     t_loop = MPI_Wtime()

!     if (params%rank == 0) write(*, '(A, 2(A, i0))') repeat('  ', i_level+1), 'Upwards GS sweep lvl ', i_level, ' it ', it_sweep

!     ! Preparation for this level, it is the inverse of the downwards steps, we have the following situation:
!     !   - ref = -1: last level that is ready to go upwards
!     !   - ref =  0: last level, that has hvy_depth > 0 and therefore will stay on this level
!     ! We will now go through all blocks and set refinement flag after the mother
!     t_block = MPI_Wtime()
!     do k_block = 1, hvy_n(tree_ID)
!         hvy_id = hvy_active(k_block, tree_ID)
!         call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

!         ! activate daughters, so that they can receive the residual
!         ! if daughter block exists and any/all have ref=-1: this block will be on subgrid in next iteration, so we set it a new refinement value (0 for now)
!         if (.not. block_is_root(params, hvy_id)) then
!             if (lgt_block(hvy_family(hvy_ID, 1),IDX_REFINE_STS) == -1) then
!                 lgt_block(lgt_id,IDX_REFINE_STS) = 0
!                 ! this following is set that we know, that this block is new on this level and we need to reconstruct it after upwards sync
!                 hvy_depth(k_block) = -hvy_depth(k_block)-1
!             endif
!         endif

!         ! we already check for hvy_depth for the next level, in order to synchronize only once per loop
!         ! they will get a temporary flag which will be changed soon
!         if (lgt_block(lgt_id,IDX_REFINE_STS) == 0 .and. (hvy_depth(k_block) > 0 .or. hvy_depth(k_block) < -1)) then
!             lgt_block(lgt_id,IDX_REFINE_STS) = REF_UNSIGNIFICANT_STAY
!         endif
!     enddo
!     ! synch the refinement stati
!     call synchronize_lgt_data( params, refinement_status_only=.true.)
!     call toc( "GS Upwards - Level definition", 10030, MPI_Wtime()-t_block )

!     ! do upwards sync of solution field, this then needs to be interpolated
!     t_block = MPI_Wtime()
!     call sync_M2D(params, hvy_sol(:,:,:,1:nc,:), tree_ID, sync_case="ref", s_val=-1)
!     call toc( "Sync M2D", 10012, MPI_Wtime()-t_block )

!     ! reordering of refinement flags:
!     !                        -1 -> REF_TMP_EMPTY : these blocks have synched and are finished
!     !                         0 -> -1 : these blocks will sync in the next iteration
!     !    REF_UNSIGNIFICANT_STAY ->  0 : these blocks will stay on this level in the next iteration
!     do k_block = 1, lgt_n(tree_ID)
!         lgt_id = lgt_active(k_block, tree_ID)
!         if (lgt_block(lgt_id,IDX_REFINE_STS) == -1) then
!             lgt_block(lgt_id,IDX_REFINE_STS) = REF_TMP_EMPTY
!         elseif (lgt_block(lgt_id,IDX_REFINE_STS) == 0) then
!             lgt_block(lgt_id,IDX_REFINE_STS) = -1
!         elseif (lgt_block(lgt_id,IDX_REFINE_STS) == REF_UNSIGNIFICANT_STAY) then
!             lgt_block(lgt_id,IDX_REFINE_STS) = 0
!         endif
!     enddo

!     ! we need to synch the solution to reconstruct it
!     t_block = MPI_Wtime()
!     call sync_ghosts_tree(params, hvy_sol(:,:,:,1:nc,:), tree_ID)
!     call toc( "Sync Layer", 10010, MPI_Wtime()-t_block )

!     t_block = MPI_Wtime()
!     do k_block = 1, hvy_n(tree_ID)
!         hvy_id = hvy_active(k_block, tree_ID)
!         call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

!         ! skip blocks not active
!         if (all(lgt_block(lgt_id,IDX_REFINE_STS) /= (/-1, 0/))) cycle

!         ! only reconstruct blocks that received values, these were hacked with negative hvy_depth values (to avoid any ref logic or syncing stuff)
!         if (hvy_depth(k_block) < 0) then
!             ! ! now do wavelet reconstruction - we only need interpolation (this is inefficient as it does not use prediction optimization)
!             ! call waveletReconstruction_block(params, hvy_sol(:,:,:,1:nc,hvy_id), SC_reconstruct=.true., WC_reconstruct=.false.)

!             do ic = 1,nc
!                 if (params%dim==2) then
!                     call prediction(hvy_sol(1:params%Bs(1)+2*params%g:2,1:params%Bs(2)+2*params%g:2,:,ic,hvy_id), &
!                                     hvy_sol(1:params%Bs(1)+2*params%g-1,1:params%Bs(2)+2*params%g-1,:,ic,hvy_id), params%order_predictor)
!                 else
!                     call prediction(hvy_sol(1:params%Bs(1)+2*params%g:2,1:params%Bs(2)+2*params%g:2,1:params%Bs(3)+2*params%g:2,ic,hvy_id), &
!                                     hvy_sol(1:params%Bs(1)+2*params%g-1,1:params%Bs(2)+2*params%g-1,1:params%Bs(3)+2*params%g-1,ic,hvy_id), params%order_predictor)
!                 endif
!             enddo
            
!             ! now we need to reset the depth counter to the correct value
!             hvy_depth(k_block) = -hvy_depth(k_block)-1
!         endif

!         ! get spacing
!         call get_block_spacing_origin( params, lgt_ID, x0, dx )

!         ! restore RHS b from residual
!         call GS_compute_residual(params, hvy_work(:,:,:,hvy_depth(k_block)*nc+1:hvy_depth(k_block)*nc+nc,hvy_id), hvy_RHS(:,:,:,1:nc,hvy_id), hvy_RHS(:,:,:,1:nc,hvy_id), a, stencil, dx, recompute_b=.true.)

!         ! treecode = get_tc(lgt_block(lgt_id, IDX_TC_1 : IDX_TC_2))
!         ! write(fname, '(A, i0, A,i0, A,i0, A, i0, A)') "block_dumped_lvl", lgt_block(lgt_id, IDX_MESH_LVL),"_it", i_level ,"_tc=", treecode, ".t"
!         ! call dump_block_fancy(hvy_work(:,:,:,hvy_depth(k_block)*nc+1:hvy_depth(k_block)*nc+nc,hvy_id,1), fname, params%Bs, params%g, digits=2, print_ghosts=.false.)
!         ! call dump_block_fancy(hvy_sol(:,:,:,1:nc,hvy_id), fname, params%Bs, params%g, digits=2, print_ghosts=.false.)

!         ! at last, we can add the interpolated residual or that of this level to the solution
!         hvy_sol(:,:,:,1:nc,hvy_id) = hvy_sol(:,:,:,1:nc,hvy_id) + hvy_work(:,:,:,hvy_depth(k_block)*nc+1:hvy_depth(k_block)*nc+nc,hvy_id)
!     enddo
!     call toc( "GS Upwards - Solution reconstruction", 10031, MPI_Wtime()-t_block )

!     ! do actual sweeps
!     sweep_forward = .true.
!     do i_sweep = 1, it_sweep
!         ! we need to synch before a sweep
!         t_block = MPI_Wtime()
!         call sync_ghosts_tree(params, hvy_sol(:,:,:,1:nc,:), tree_ID, a, a, ignore_Filter=.true.)
!         call toc( "Sync Layer", 10010, MPI_Wtime()-t_block )

!         ! blocks on this iteration do a GS-sweep, they have refinement status 0 or -1
!         call GS_iteration_ref(params, tree_id, (/-1, 0 /), hvy_sol(:,:,:,1:nc,:), hvy_RHS(:,:,:,1:nc,:), a, stencil, sweep_forward)

!         sweep_forward = .not. sweep_forward
!     enddo

!     ! after sweeps, we only need to decrease hvy_depth counters! we already marked blocks with -1 or 0 depending on their hvy_depth status
!     ! this comes from the downwards iterations, where graded grids were always ensured, so we just reverse that order and everything is fine
!     do k_block = 1, hvy_n(tree_ID)
!         hvy_id = hvy_active(k_block, tree_ID)
!         call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
!         if (lgt_block(lgt_id,IDX_REFINE_STS) == 0 .and. hvy_depth(k_block) > 0) hvy_depth(k_block) = hvy_depth(k_block) - 1
!     enddo

!     t_print = MPI_Wtime()-t_loop
!     call MPI_ALLREDUCE(MPI_IN_PLACE, t_print, 1, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)
!     call append_t_file('multigrid-iteration.t', (/1.0_rk, dble(i_level), t_print/))

!     call toc( "GS Upwards - full sweep", 10005, MPI_Wtime()-t_loop )

! end subroutine multigrid_upwards


!> \brief Upwards level iteration for the multigrid solver
!> \details Compute for one level an upwards step
!> Does three things:
!>   - sync the coarser solution from the lower lewel as WD decomposed values, reconstruct the coarser solution
!>   - do GS sweeps on this level
subroutine multigrid_upwards_2(params, hvy_sol, hvy_RHS, hvy_work, tree_ID, Jmin, Jmax_a, a, stencil, it_sweep)
    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    real(kind=rk), intent(inout)       :: hvy_sol(:, :, :, :, :), hvy_RHS(:, :, :, :, :), hvy_work(:, :, :, :, :)
    integer(kind=ik), intent(in)       :: tree_ID
    integer(kind=ik), intent(in)       :: Jmin
    integer(kind=ik), intent(in)       :: Jmax_a
    integer(kind=ik), intent(in)       :: a
    real(kind=rk), intent(in)          :: stencil(-a:a)
    integer(kind=ik), intent(in)       :: it_sweep

    character(len=cshort)              :: fname
    integer(kind=tsize)                :: treecode

    logical            :: sweep_forward
    integer(kind=ik)   :: k_block, lgt_ID, hvy_id, nc, ic, i_sweep, mpierr, i_level, level_me, ref_me, sync_freq, i_g(1:3), i_BS(1:3)
    real(kind=rk)      :: t_loop, t_block, t_print(1:1)
    real(kind=rk)      :: dx(1:3), x0(1:3)
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

        if (params%rank == 0) write(*, '(A, 2(A, i0))') repeat('  ', i_level+1), 'Upwards GS sweep lvl ', i_level, ' it ', it_sweep

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


        ! in case that we hit Jmax_a, then we do not solve Ae=r but Au=b and we need to restore the old solution as well as the RHS b
        ! this is happening for all leaf-blocks
        t_block = MPI_Wtime()
        if (i_level == Jmax_a) then
            do k_block = 1, hvy_n(tree_ID)
                hvy_id = hvy_active(k_block, tree_ID)
                call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

                if (block_is_leaf(params, hvy_id)) then
                    ! get spacing
                    call get_block_spacing_origin( params, lgt_ID, x0, dx )

                    ! recompute RHS b = r + Ax, 
                    call GS_compute_residual(params, hvy_work(:,:,:,1:nc,hvy_id), hvy_RHS(:,:,:,1:nc,hvy_id), hvy_RHS(:,:,:,1:nc,hvy_id), a, stencil, dx, recompute_b=.true.)
            
                    ! add reconstructed solution to residual to previous solution
                    hvy_sol(:,:,:,1:nc,hvy_id) = hvy_sol(:,:,:,1:nc,hvy_id) + hvy_work(:,:,:,1:nc,hvy_id)
                endif
            enddo

            ! sync RHS for sync_freq > 1
            t_block = MPI_Wtime()
            call sync_ghosts_tree(params, hvy_RHS(:,:,:,1:nc,:), tree_ID)
            call toc( "Sync Layer", 10010, MPI_Wtime()-t_block )
        endif
        call toc( "GS Upwards - RHS and solution restoration", 10032, MPI_Wtime()-t_block )

        ! do actual sweeps
        ! if (i_level /= Jmax_a) then
            sync_freq = 1
            ! sync_freq = params%g / a
            sweep_forward = .true.
            do i_sweep = 1, it_sweep
                if (modulo(i_sweep-1, sync_freq) == 0) then
                    ! we need to synch before a sweep
                    t_block = MPI_Wtime()
                    call sync_ghosts_tree(params, hvy_sol(:,:,:,1:nc,:), tree_ID, a, a, ignore_Filter=.true.)
                    ! call sync_ghosts_tree(params, hvy_sol(:,:,:,1:nc,:), tree_ID, a*sync_freq, a*sync_freq)
                    call toc( "Sync Layer", 10010, MPI_Wtime()-t_block )
                endif

                ! blocks on this iteration do a GS-sweep, they have refinement status 0 or 1
                ! call GS_iteration_ref(params, tree_id, (/ 1, 0 /), hvy_sol(:,:,:,1:nc,:), hvy_RHS(:,:,:,1:nc,:), a, stencil, sweep_forward)
                call GS_iteration_ref(params, tree_id, (/ 1, 0 /), hvy_sol(:,:,:,1:nc,:), hvy_RHS(:,:,:,1:nc,:), a, stencil, sweep_forward, filter_offset=params%g-(sync_freq-1)*a)

                sweep_forward = .not. sweep_forward
            enddo
        ! endif

        t_print = MPI_Wtime()-t_loop
        call MPI_ALLREDUCE(MPI_IN_PLACE, t_print, 1, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)
        call append_t_file('multigrid-iteration.t', (/1.0_rk, dble(i_level), t_print/))

        call toc( "GS Upwards - full sweep", 10005, MPI_Wtime()-t_loop )
    enddo

end subroutine multigrid_upwards_2



!> \brief Solve for the solution on the coarsest level
!> \details Solves for the solution on the coarsest level
!> It solves this with one of the following three modes:
!>   - FFT, but only for level 0
!>   - Conjugent gradient method, (but only for level 0 for now)
!>   - GS sweeps
subroutine multigrid_coarsest(params, hvy_sol, hvy_RHS, tree_ID, i_level, Jmax_a, a, stencil, it_sweep, tol_cg, fft_order_discretization)
    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    real(kind=rk), intent(inout)       :: hvy_sol(:, :, :, :, :), hvy_RHS(:, :, :, :, :)
    integer(kind=ik), intent(in)       :: tree_ID
    integer(kind=ik), intent(in)       :: i_level
    integer(kind=ik), intent(in)       :: Jmax_a
    integer(kind=ik), intent(in)       :: a
    real(kind=rk), intent(in)          :: stencil(-a:a)
    integer(kind=ik), intent(in)       :: it_sweep
    real(kind=rk), intent(in)          :: tol_cg
    character(len=cshort), intent(in)  :: fft_order_discretization

    character(len=cshort)              :: fname
    integer(kind=tsize)                :: treecode

    logical            :: sweep_forward
    integer(kind=ik)   :: k_block, lgt_ID, hvy_id, nc, i_sweep, mpierr
    real(kind=rk)      :: dx(1:3), x0(1:3)
    real(kind=rk)      :: t_block, t_print(1:1)
    nc = size(hvy_sol,4)

    if (i_level == 0) then
        if (it_sweep == -1) then
            ! FFT for lowest block if level = 0
            if (params%rank == 0) write(*, '(A, A, i0)') repeat('  ', i_level+1), 'FFT lvl ', i_level
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

                call fft_solve_poisson(params, hvy_sol(:,:,:,1:nc,hvy_id), hvy_RHS(:,:,:,1:nc,hvy_id), fft_order_discretization, dx)

                ! call dump_block_fancy(hvy_sol(:,:,:,1:nc,hvy_id), 'sol_L0.txt', params%Bs, params%g, digits=2, print_ghosts=.true.)
                ! call dump_block_fancy(hvy_RHS(:,:,:,1:nc,hvy_id), 'RHS_L0.txt', params%Bs, params%g, digits=2, print_ghosts=.true.)

            enddo
            t_print = MPI_Wtime()-t_block
            call MPI_ALLREDUCE(MPI_IN_PLACE, t_print, 1, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)
            call append_t_file('multigrid-iteration.t', (/0.0_rk, dble(i_level), t_print/))
            call toc( "fft solve poisson", 10004, MPI_Wtime()-t_block )
        else
            ! conjugent gradient method for lowest block if level = 0
            if (params%rank == 0) write(*, '(A, 2(A, i0), A, es7.1)') repeat('  ', i_level+1), 'Coarsest CG lvl ', i_level, ' it max ', it_sweep, ' tol ', tol_cg
            t_block = MPI_Wtime()
            do k_block = 1, hvy_n(tree_ID)
                hvy_id = hvy_active(k_block, tree_ID)
                call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

                ! skip blocks not on this level
                if (lgt_block(lgt_id,IDX_MESH_LVL) /= i_level) cycle

                ! get spacing
                call get_block_spacing_origin( params, lgt_ID, x0, dx )
                
                call CG_solve_poisson_level0(params, hvy_sol(:,:,:,1:nc,hvy_id), hvy_RHS(:,:,:,1:nc,hvy_id), a, stencil, dx, tol_cg, it_sweep)
            enddo
            t_print = MPI_Wtime()-t_block
            call MPI_ALLREDUCE(MPI_IN_PLACE, t_print, 1, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)
            call append_t_file('multigrid-iteration.t', (/0.0_rk, dble(i_level), t_print/))
            call toc( "conjugent gradient", 10004, MPI_Wtime()-t_block )
        endif
    else
        ! coarsest level sweeps
        if (params%rank == 0) write(*, '(A, 2(A, i0))') repeat('  ', i_level+1), 'Coarsest GS sweeps lvl ', i_level, ' it ', it_sweep

        sweep_forward = .true.
        do i_sweep = 1, it_sweep
            if (i_sweep == 1) then
                ! init solution as 0
                hvy_sol(:,:,:,1:nc,hvy_id) = 0.0_rk
            else
                ! we need to synch before a sweep
                t_block = MPI_Wtime()
                call sync_level_from_M( params, hvy_sol(:,:,:,1:nc,:), tree_ID, i_level, a, a)
                call toc( "Sync Layer", 10010, MPI_Wtime()-t_block )
            endif

            ! blocks on this level do a GS-sweep
            call GS_iteration_level(params, tree_id, i_level, hvy_sol(:,:,:,1:nc,:), hvy_RHS(:,:,:,1:nc,:), a, stencil, sweep_forward)

            sweep_forward = .not. sweep_forward
        enddo
    endif

end subroutine multigrid_coarsest