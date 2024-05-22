subroutine coarseningIndicator_tree( time, params, level_this, hvy_block, hvy_tmp, &
    tree_ID, indicator, iteration, ignore_maxlevel, input_is_WD, hvy_mask)

    use module_indicators

    implicit none
    real(kind=rk), intent(in)           :: time
    type (type_params), intent(in)      :: params                         !> user defined parameter structure
    integer(kind=ik), intent(in)        :: level_this                     !> current level to look at (in the case of biorthogonal wavelets, not in "harten-multiresolution")
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)       !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)         !> heavy work data array - block data.
    !> mask data. we can use different trees (4est module) to generate time-dependent/indenpedent
    !! mask functions separately. This makes the mask routines tree-level routines (and no longer
    !! block level) so the physics modules have to provide an interface to create the mask at a tree
    !! level. All parts of the mask shall be included: chi, boundary values, sponges.
    real(kind=rk), intent(inout), optional :: hvy_mask(:, :, :, :, :)
    integer(kind=ik), intent(in)        :: tree_ID
    character(len=*), intent(in)        :: indicator                      !> how to choose blocks for refinement
    !> coarsening iteration index. coarsening is done until the grid has reached
    !! the steady state; therefore, this routine is called several times during the
    !! mesh adaptation. Random coarsening (used for testing) is done only in the first call.
    integer(kind=ik), intent(in)        :: iteration
    
    logical, intent(in)                 :: input_is_WD                       !< flag if hvy_block is already wavelet decomposed

    ! for the mask generation (time-independent mask) we require the mask on the highest
    ! level so the "force_maxlevel_dealiasing" option needs to be overwritten. Life is difficult, at times.
    logical, intent(in)                 :: ignore_maxlevel

    ! local variables
    integer(kind=ik) :: k, Jmax, neq, lgtID, g, mpierr, hvyID, p, N_thresholding_components, tags, ierr, level
    integer(kind=ik), dimension(3) :: Bs
    ! local block spacing and origin
    real(kind=rk) :: dx(1:3), x0(1:3), crsn_chance, R
    real(kind=rk), allocatable, save :: norm(:)
    logical :: consider_hvy_tmp, inputIsWD
    real(kind=rk) :: t0  !< timing for debugging

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas

    ! in the default case we threshold all statevector components
    N_thresholding_components = params%n_eqn
    consider_hvy_tmp = .false.
    inputIsWD = input_is_WD
    Jmax = params%Jmax
    neq = params%n_eqn
    Bs = params%Bs
    g = params%g

    ! reset refinement status to "stay" on all active blocks
    do k = 1, lgt_n(tree_ID)
        lgtID = lgt_active(k, tree_ID)
        lgt_block( lgtID, IDX_REFINE_STS ) = 0
    enddo

    ! construct mask function, if it is used as secondary criterion. This criterion
    ! ensures that regions with gradients in the mask function (the fluid/solid interface)
    ! are not coarsened (except in the "dealiasing step", because all blocks on Jmax are coarsened)
    if (params%threshold_mask .and. present(hvy_mask) .and. indicator/="everywhere" .and. indicator/="random") then
        t0 = MPI_Wtime()
        ! Note the "all_parts=.false." means that we do not bypass the pruned trees. This functionality should
        ! work as designed, but use it carefully, as it is still developped. If the PARAMS file sets
        ! params%dont_use_pruned_tree_mask=1, it is deactivated anyways.
        call createMask_tree(params, time, hvy_mask, hvy_tmp, all_parts=.false.)
        call toc( "coarseningIndicator (createMask_tree)", MPI_Wtime()-t0 )
    endif


    ! the indicator primary-variables is for compressible Navier-Stokes and first
    ! converts the actual state vector to "traditional" quantities (the primary variables)
    ! and then applies thresholding to that. The result is stored in hvy_tmp
    ! and hence we use this for thresholding
    if ( indicator == "primary-variables" ) then
        consider_hvy_tmp = .true.
    endif

    !---------------------------------------------------------------------------
    !> Preparation for thresholding (if required)
    !---------------------------------------------------------------------------
    !! The idea is that you want to threshold something else than your statevector
    !! whatever reason. For example, vorticity, or in the skew-symmetric case you want
    !! to apply some statevector conversion first.
    !! it is a little unfortunate that we have to do this preparation here and not in
    !! coarseningIndicator_block, where it should be. but after computing it, we have to
    !! sync its ghost nodes in order to apply the detail operator to the entire
    !! derived field (incl gost nodes).
    if (consider_hvy_tmp) then
        call abort(197,"This case currently does not work as Julius needs to sort out the temporary arrays")
        inputIsWD = .false.

        ! t0 = MPI_Wtime()
        ! ! case with derived quantities.
        ! ! loop over my active hvy data:
        ! do k = 1, hvy_n(tree_ID)
        !     hvyID = hvy_active(k, tree_ID)

        !     ! get lgt id of block
        !     call hvy2lgt( lgtID, hvyID, params%rank, params%number_blocks )

        !     ! some indicators may depend on the grid (e.g. the vorticity), hence
        !     ! we pass the spacing and origin of the block
        !     call get_block_spacing_origin( params, lgtID, x0, dx )

        !     ! actual computation of thresholding quantity (vorticity etc)
        !     call PREPARE_THRESHOLDFIELD_meta( params%physics_type, time, hvy_block(:,:,:,:,hvyID), &
        !     g, x0, dx, hvy_tmp(:,:,:,:,hvyID), hvy_mask(:,:,:,:,hvyID), N_thresholding_components )
        ! enddo

        ! ! note here we sync hvy_tmp (=derived qty) and not hvy_block
        ! call sync_ghosts_all( params, lgt_block, hvy_tmp(:,:,:,1:N_thresholding_components,:), hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID) )

        ! if (params%threshold_mask .and. N_thresholding_components /= params%n_eqn) &
        ! call abort(2801191,"your thresholding does not work with threshold-mask.")

        ! call toc( "coarseningIndicator (prepare thresholdfield)", MPI_Wtime()-t0 )
    endif



    !---------------------------------------------------------------------------
    !> Compute normalization for eps, if desired.
    !---------------------------------------------------------------------------
    !! versions <30.05.2018 used fixed eps for all qtys of the state vector, but that is not very smart
    !! as each qty can have different mangitudes. If the switch eps_normalized is on, we compute here
    !! the vector of normalization factors for each qty that adaptivity will be based on (state vector
    !! or vorticity). The norm is specified in params%eps_norm, default is Linfty.

    !! default norm is 1 so in this case eps is an absolute value.
    if (.not. allocated(norm)) allocate(norm(1:N_thresholding_components))
    norm = 1.0_rk

    ! if we coarsen randomly or everywhere, well, why compute the norm ?
    if ( params%eps_normalized .and. indicator/="everywhere" .and. indicator/="random" ) then
        t0 = MPI_Wtime()
        if ( .not. consider_hvy_tmp ) then
            ! Apply thresholding directly to the statevector (hvy_block), not to derived quantities
            call componentWiseNorm_tree(params, hvy_block, tree_ID, params%eps_norm, norm)
        else
            ! use derived qtys instead (hvy_tmp)
            call componentWiseNorm_tree(params, hvy_tmp, tree_ID, params%eps_norm, norm)
        endif

        ! HACK
!        if (params%physics_type == "ACM-new") then
!            if (params%eps_norm == "Linfty") then
!                norm(1:params%dim) = maxval( norm(1:params%dim) )
!            elseif (params%eps_norm == "L2") then
!                norm(1:params%dim) = sqrt( sum(norm(1:params%dim)**2) )
!            endif
        ! endif

        ! avoid division by zero (corresponds to using an absolute eps if the norm is very small)
        do p = 1, N_thresholding_components
            if (norm(p) <= 1.0e-9_rk) norm(p) = 1.0_rk
        enddo

        ! during dev is it useful to know what the normalization is, if that is active
        call append_t_file('eps_norm.t', (/time, norm, params%eps/))

        call toc( "coarseningIndicator (norm)", MPI_Wtime()-t0 )
    endif

    !---------------------------------------------------------------------------
    !> evaluate coarsening criterion on all blocks
    !---------------------------------------------------------------------------
    ! the indicator "random" requires special treatment below (it does not pass via
    ! coarseningIndicator_block)
    select case(indicator)
    case ("everywhere")
        ! coarsen all blocks. Note: it is not always possible (with any grid) to do that!
        ! The flag may be removed again (completeness, Jmin, gradedness)...
        ! Done only in the first iteration of grid adaptation,
        ! because otherwise, the grid would continue changing, and the coarsening
        ! stops only if all blocks are at Jmin.
        if (iteration==0) then
            do k = 1, lgt_n(tree_ID)
                ! flag for coarsening
                lgt_block(lgt_active(k, tree_ID), IDX_REFINE_STS) = -1
            enddo
        endif

    case ("random")
        ! random coarsening. Done only in the first iteration of grid adaptation,
        ! because otherwise, the grid would continue changing, and the coarsening
        ! stops only if all blocks are at Jmin.
        if (iteration==0) then
            ! set random seed
            call init_random_seed()
            ! only root tags blocks (can be messy otherwise)
            if (params%rank == 0) then
                tags = 0
                ! the chance for coarsening: (note to coarsen, all 4 or 8 sisters need to have this status,
                ! hence the probability for a block to actually coarsen is only crsn_chance**(2^D))
                crsn_chance = (0.50_rk)**(1.0_rk / 2.0_rk**params%dim)

                do k = 1, lgt_n(tree_ID)
                    lgtID = lgt_active(k, tree_ID)
                    ! random number
                    call random_number(r)
                    ! set refinement status to coarsen based on random numbers.
                    if ( r <= crsn_chance ) then
                        lgt_block(lgtID, IDX_REFINE_STS) = -1
                        tags = tags + 1
                    endif
                enddo
            endif
            ! sync light data, as only root sets random refinement
            call MPI_BCAST( lgt_block(:, IDX_REFINE_STS), size(lgt_block,1), MPI_INTEGER4, 0, WABBIT_COMM, ierr )
        endif

    case default
        t0 = MPI_Wtime()
        ! Default is wavelet thresholding...
        
        ! NOTE: even if additional mask thresholding is used, passing the mask is optional,
        ! notably because of the ghost nodes unit test, where random refinement / coarsening
        ! is used. hence, checking the flag params%threshold_mask alone is not enough.
        do k = 1, hvy_n(tree_ID)
            hvyID = hvy_active(k, tree_ID)
            call hvy2lgt( lgtID, hvyID, params%rank, params%number_blocks )
            level = lgt_block( lgtID, IDX_MESH_LVL)

            ! level wise coarsening: in the "biorthogonal" case, we start at J_max_active and
            ! iterate down to J_min. Only blocks on the level "level_this" are allowed to coarsen.
            ! this should prevent filtering artifacts at block-block interfaces.
            if (level /= level_this) cycle

            ! force blocks on maximum refinement level to coarsen, if parameter is set.
            ! Note this behavior can be bypassed using the ignore_maxlevel switch.
            if (params%force_maxlevel_dealiasing .and. .not. ignore_maxlevel .and. (level==Jmax)) then
                ! coarsen (no need to evaluate the possibly expensive coarseningIndicator_block)
                lgt_block(lgtID, IDX_REFINE_STS) = -1
            else
                ! evaluate the criterion on this block.
                if (params%threshold_mask .and. present(hvy_mask)) then
                    call coarseningIndicator_block( params, hvy_block(:,:,:,:,hvyID), &
                    hvy_tmp(:,:,:,:,hvyID), indicator, &
                    lgt_block(lgtID, IDX_REFINE_STS), norm, level, inputIsWD, hvy_mask(:,:,:,:,hvyID))
                else
                    call coarseningIndicator_block( params, hvy_block(:,:,:,:,hvyID), &
                    hvy_tmp(:,:,:,:,hvyID), indicator, &
                    lgt_block(lgtID, IDX_REFINE_STS), norm, level, inputIsWD)
                endif
            endif
        enddo
        call toc( "coarseningIndicator (coarseIndicator_block)", MPI_Wtime()-t0 )
    end select

    !> after modifying all refinement flags, we need to synchronize light data
    call synchronize_lgt_data( params,  refinement_status_only=.true. )

end subroutine
