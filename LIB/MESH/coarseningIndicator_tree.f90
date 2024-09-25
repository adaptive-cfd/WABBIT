!> \brief Decides for all blocks if they can stay or want to coarsen.
!> This method goes a bit against the naming convention, as for default wavelet cases it acts
!! level-wise or leaf-wise but for specific indicators (everywhere or random) it acts on the whole tree.
subroutine coarseningIndicator_tree( time, params, hvy_block, hvy_tmp, &
    tree_ID, indicator, ignore_maxlevel, input_is_WD, hvy_mask, norm_inout)

    use module_indicators

    implicit none
    real(kind=rk), intent(in)           :: time                           !< Time used for mask generation
    type (type_params), intent(inout)   :: params                         !< user defined parameter structure
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)       !< heavy data array - block data
    !> heavy work data array - block data, if input_is_wd is true this contains the block values
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    !> mask data. we can use different trees (4est module) to generate time-dependent/indenpedent
    !! mask functions separately. This makes the mask routines tree-level routines (and no longer
    !! block level) so the physics modules have to provide an interface to create the mask at a tree
    !! level. All parts of the mask shall be included: chi, boundary values, sponges.
    real(kind=rk), intent(inout), optional :: hvy_mask(:, :, :, :, :)
    integer(kind=ik), intent(in)        :: tree_ID
    character(len=*), intent(in)        :: indicator                      !< how to choose blocks for refinement
    logical, intent(in)                 :: input_is_WD                    !< flag if hvy_block is already wavelet decomposed
    !> for the mask generation (time-independent mask) we require the mask on the highest
    !! level so the "force_maxlevel_dealiasing" option needs to be overwritten. Life is difficult, at times.
    logical, intent(in)                 :: ignore_maxlevel
    real(kind=rk), intent(inout), optional :: norm_inout(:)  !< We can output the norm as well

    ! local variables
    integer(kind=ik) :: k_b, k_nc, Jmax, n_eqn, lgt_ID, g, mpierr, hvy_ID, p, N_thresholding_components, tags, tags_all, ierr, ref_stat, level
    integer(kind=ik), dimension(3) :: Bs
    ! local block spacing and origin
    real(kind=rk) :: dx(1:3), x0(1:3), crsn_chance, R
    real(kind=rk), allocatable, save :: norm(:)
    logical :: consider_hvy_tmp, inputIsWD, check_mask
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
    n_eqn = params%n_eqn
    Bs = params%Bs
    g = params%g

    ! reset refinement status to "stay" on all active blocks
    ! caution: if the code will ever be adapted again to not work on whole tree then this needs to be adapted
    do k_b = 1, lgt_n(tree_ID)
        lgt_ID = lgt_active(k_b, tree_ID)
        lgt_block( lgt_ID, IDX_REFINE_STS ) = 0
    enddo

    ! construct mask function, if it is used as secondary criterion. This criterion
    ! ensures that regions with gradients in the mask function (the fluid/solid interface)
    ! are not coarsened (except for dealiasing, because all blocks on Jmax are coarsened)
    if (params%threshold_mask .and. present(hvy_mask) .and. indicator/="everywhere" .and. indicator/="random") then

        t0 = MPI_Wtime()
        ! Note the "all_parts=.false." means that we do not bypass the pruned trees. This functionality should
        ! work as designed, but use it carefully, as it is still developped. If the PARAMS file sets
        ! params%dont_use_pruned_tree_mask=1, it is deactivated anyways.
        ! While this uses hvy_tmp as well, we use different tree_ID so there is no clash
        call createMask_tree(params, time, hvy_mask, hvy_tmp, all_parts=.false.)
        call toc( "coarseningIndicator (createMask_tree)", 120, MPI_Wtime()-t0 )
        ! simplify if-condition downwards and check mask for indicator if it was created here
        check_mask = .true.
    else
        check_mask = .false.
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
        !     call hvy2lgt( lgt_ID, hvyID, params%rank, params%number_blocks )

        !     ! some indicators may depend on the grid (e.g. the vorticity), hence
        !     ! we pass the spacing and origin of the block
        !     call get_block_spacing_origin( params, lgt_ID, x0, dx )

        !     ! actual computation of thresholding quantity (vorticity etc)
        !     call PREPARE_THRESHOLDFIELD_meta( params%physics_type, time, hvy_block(:,:,:,:,hvyID), &
        !     g, x0, dx, hvy_tmp(:,:,:,:,hvyID), hvy_mask(:,:,:,:,hvyID), N_thresholding_components )
        ! enddo

        ! ! note here we sync hvy_tmp (=derived qty) and not hvy_block
        ! call sync_ghosts_all( params, hvy_tmp(:,:,:,1:N_thresholding_components,:), tree_ID )

        ! if (params%threshold_mask .and. N_thresholding_components /= params%n_eqn) &
        ! call abort(2801191,"your thresholding does not work with threshold-mask.")

        ! call toc( "coarseningIndicator (prepare thresholdfield)", 121, MPI_Wtime()-t0 )
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
    norm = -99999  ! reset as random value

    ! if we coarsen randomly or everywhere, well, why compute the norm ?
    if ( params%eps_normalized .and. indicator/="everywhere" .and. indicator/="random" .and. indicator/="threshold-cvs" .and. indicator/="threshold-image-denoise") then
        t0 = MPI_Wtime()
        if ( .not. consider_hvy_tmp .and. .not. input_is_WD) then
            ! Apply thresholding directly to the statevector (hvy_block), not to derived quantities
            call componentWiseNorm_tree(params, hvy_block, tree_ID, params%eps_norm, norm)
        else
            ! use derived qtys instead (hvy_tmp) or pre-saved values for leaf-wise loop
            ! For adapt_tree leaf-wise we ensure that at this spot all blocks are wavelet decomposed and original values copied to hvy_tmp
            ! norm is computed from original values, not WD values
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

        call toc( "coarseningIndicator (norm)", 122, MPI_Wtime()-t0 )
    else
        norm = 1.0_rk  ! set to 1
    endif
    if (present(norm_inout)) then
        norm_inout(1:N_thresholding_components) = norm(1:N_thresholding_components)
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
        ! Needs to be a heavy loop as elsewise we don't know which is a leaf-block
        do k_b = 1, hvy_n(tree_ID)
            hvy_ID = hvy_active(k_b, tree_ID)
            call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)

            ! leaf-blocks only are being set
            if (.not. block_is_leaf(params, hvy_ID)) cycle

            ! flag for coarsening
            lgt_block(lgt_ID, IDX_REFINE_STS) = -1
        enddo

    case ("random")
        ! random coarsening, set random seed after cpu clock
        call init_random_seed()

        tags = 0
        ! the chance for coarsening: (note to coarsen, all 4 or 8 sisters need to have this status,
        ! hence the probability for a block to actually coarsen is only crsn_chance**(2^D))
        crsn_chance = (0.50_rk)**(1.0_rk / 2.0_rk**params%dim)

        ! needs to be heavy as elsewise we don't know which is a leaf-block
        do k_b = 1, hvy_n(tree_ID)
            hvy_ID = hvy_active(k_b, tree_ID)
            call hvy2lgt(lgt_ID, hvy_ID, params%rank, params%number_blocks)

            ! leaf-blocks only are being set
            if (.not. block_is_leaf(params, hvy_ID)) cycle

            ! random number
            call random_number(r)
            ! set refinement status to coarsen based on random numbers or let them stay
            if ( r <= crsn_chance ) then
                lgt_block(lgt_ID, IDX_REFINE_STS) = -1
                tags = tags + 1
            else
                lgt_block(lgt_ID, IDX_REFINE_STS) = 0
            endif
        enddo

        ! ! in order to test: sum of tags over processors can be compared
        ! call MPI_ALLREDUCE(MPI_IN_PLACE, tags, 1, MPI_INTEGER4, MPI_SUM, WABBIT_COMM, mpierr)
        ! if (params%rank==0) write(*, '(A, es10.3, A, es10.3)') "Theoretical coarsening chance: ", crsn_chance, " , actual coarsening fraction: ", dble(tags) / dble(lgt_n(tree_ID)) * 2.0_rk**params%dim / (2.0_rk**params%dim - 1.0_rk)

    case default
        t0 = MPI_Wtime()
        ! Default is wavelet thresholding...

        ! if (params%rank == 0) then
        !     write(*, '(A, 10(es10.3, 2x))') "Norm: ", norm(:)
        ! endif
        
        ! NOTE: even if additional mask thresholding is used, passing the mask is optional,
        ! notably because of the ghost nodes unit test, where random refinement / coarsening
        ! is used. hence, checking the flag params%threshold_mask alone is not enough.
        do k_b = 1, hvy_n(tree_ID)
            hvy_ID = hvy_active(k_b, tree_ID)
            call hvy2lgt( lgt_ID, hvy_ID, params%rank, params%number_blocks )
            level = lgt_block( lgt_ID, IDX_MESH_LVL)

            ! force blocks on maximum refinement level to coarsen, if parameter is set.
            ! Note this behavior can be bypassed using the ignore_maxlevel switch.
            if (params%force_maxlevel_dealiasing .and. .not. ignore_maxlevel .and. (level==Jmax)) then
                ! coarsen (no need to evaluate the possibly expensive coarseningIndicator_block)
                lgt_block(lgt_ID, IDX_REFINE_STS) = -1
            else
                ! evaluate the criterion on this block.
                if (check_mask) then
                    call coarseningIndicator_block( params, hvy_block(:,:,:,:,hvy_ID), &
                    hvy_tmp(:,:,:,:,hvy_ID), indicator, &
                    lgt_block(lgt_ID, IDX_REFINE_STS), norm, level, inputIsWD, hvy_mask(:,:,:,:,hvy_ID))
                else
                    call coarseningIndicator_block( params, hvy_block(:,:,:,:,hvy_ID), &
                    hvy_tmp(:,:,:,:,hvy_ID), indicator, &
                    lgt_block(lgt_ID, IDX_REFINE_STS), norm, level, inputIsWD)
                endif
            endif
        enddo
        call toc( "coarseningIndicator (coarseIndicator_block)", 123, MPI_Wtime()-t0 )

    end select

    ! after modifying all refinement flags, we need to synchronize light data
    call synchronize_lgt_data( params,  refinement_status_only=.true. )

end subroutine