
! ============================================================================================
!> \name coarsening_indicator.f90
!> \version 0.5
!> \author engels
!> \brief Set coarsening status for all active blocks, different methods possible
!
!> \details This routine sets the coarsening flag for all blocks. We allow for different
!! mathematical methods (everywhere / random) currently not very complex, but expected to grow
!! in the future.
!! \n
!! ------------------ \n
!! Refinement status: \n
!! ------------------ \n
!! +1 refine \n
!! 0 do nothing \n
!! -1 block wants to coarsen (ignoring other constraints, such as gradedness) \n
!! -2 block will coarsen and be merged with her sisters \n
!! ------------------ \n
!! \n
!! = log ======================================================================================
!! \n
!! 29/05/2018 create
! ********************************************************************************************
subroutine grid_coarsening_indicator( time, params, lgt_block, hvy_block, hvy_tmp, lgt_active, &
    lgt_n, lgt_sortednumlist, hvy_active, hvy_n, indicator, iteration, hvy_neighbor, hvy_mask)

    use module_indicators

    implicit none
    real(kind=rk), intent(in)           :: time
    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy work data array - block data.
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    ! mask data. we can use different trees (4est module) to generate time-dependent/indenpedent
    ! mask functions separately. This makes the mask routines tree-level routines (and no longer
    ! block level) so the physics modules have to provide an interface to create the mask at a tree
    ! level. All parts of the mask shall be included: chi, boundary values, sponges.
    real(kind=rk), intent(inout), optional :: hvy_mask(:, :, :, :, :)
    !> list of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_n
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(inout)     :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(inout)     :: hvy_n
    !> how to choose blocks for refinement
    character(len=*), intent(in)        :: indicator
    !> coarsening iteration index. coarsening is done until the grid has reached
    !! the steady state; therefore, this routine is called several times during the
    !! mesh adaptation. Random coarsening (used for testing) is done only in the first call.
    integer(kind=ik), intent(in)        :: iteration
    !> heavy data array - neighbor data
    integer(kind=ik), intent(inout)     :: hvy_neighbor(:,:)
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), intent(inout)  :: lgt_sortednumlist(:,:)

    ! local variables
    integer(kind=ik) :: k, Jmax, neq, lgt_id, g, mpierr, hvy_id, p, N_thresholding_components, tags, ierr
    integer(kind=ik), dimension(3) :: Bs
    ! local block spacing and origin
    real(kind=rk) :: dx(1:3), x0(1:3), crsn_chance, R
    real(kind=rk), allocatable, save :: norm(:), block_norm(:), tmp(:)
    !> mask term for every grid point in this block
    integer(kind=2), allocatable, save :: mask_color(:,:,:)
    !> velocity of the solid
    real(kind=rk), allocatable, save :: us(:,:,:,:)


    Jmax = params%max_treelevel
    neq = params%n_eqn
    Bs = params%Bs
    g = params%n_ghosts

    !> reset refinement status to "stay" on all blocks
    do k = 1, lgt_n
        lgt_id = lgt_active(k)
        lgt_block( lgt_id, Jmax + IDX_REFINE_STS ) = 0
    enddo

    !---------------------------------------------------------------------------
    !> Preparation for thresholding (if required)
    !---------------------------------------------------------------------------
    !! The idea is that you want to threshold something else than your statevector
    !! whatever reason. For example, vorticity, or in the skew-symmetric case you want
    !! to apply some statevector conversion first.
    !! it is a little unfortunate that we have to do this preparation here and not in
    !! block_coarsening_indicator, where it should be. but after computing it, we have to
    !! sync its ghost nodes in order to apply the detail operator to the entire
    !! derived field (incl gost nodes).
    if (params%coarsening_indicator /= "threshold-state-vector") then
        ! case with derived quantities.

        ! loop over my active hvy data
        do k = 1, hvy_n
            hvy_id = hvy_active(k)

            ! get lgt id of block
            call hvy_id_to_lgt_id( lgt_id, hvy_id, params%rank, params%number_blocks )

            ! some indicators may depend on the grid (e.g. the vorticity), hence
            ! we pass the spacing and origin of the block
            call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

            ! actual computation of thresholding quantity (vorticity etc)
            call PREPARE_THRESHOLDFIELD_meta( params%physics_type, time, hvy_block(:,:,:,:,hvy_id), &
            g, x0, dx, hvy_tmp(:,:,:,:,hvy_id), hvy_mask(:,:,:,:,hvy_id), N_thresholding_components )
        enddo

        ! note here we sync hvy_tmp (=derived qty) and not hvy_block
        call sync_ghosts( params, lgt_block, hvy_tmp(:,:,:,1:N_thresholding_components,:), hvy_neighbor, hvy_active, hvy_n )

        if (params%threshold_mask .and. N_thresholding_components /= params%n_eqn) &
        call abort(2801191,"your thresholding does not work with threshold-mask.")

    else
        ! case without derived quantities. NOTE: while it would be nicer to have
        ! PREPARE_THRESHOLDFIELD_meta simply copy the relevant components to HVY_TMP,
        ! this overhead slows down the code and must be avoided.

        ! in the default case we threshold all statevector components
        N_thresholding_components = params%n_eqn
    endif

    if (.not. allocated(norm)) then
        allocate(norm(1:N_thresholding_components), block_norm(1:N_thresholding_components),&
        tmp(1:N_thresholding_components))
    endif


    !---------------------------------------------------------------------------
    !> Compute normalization for eps, if desired.
    !---------------------------------------------------------------------------
    !! versions <30.05.2018 used fixed eps for all qtys of the state vector, but that is not very smart
    !! as each qty can have different mangitudes. If the switch eps_normalized is on, we compute here
    !! the vector of normalization factors for each qty that adaptivity will bebased on (state vector
    !! or vorticity). We currently use the L_infty norm. I have made bad experience with L_2 norm
    !! (to be checked...)

    !! default norm (e.g. for compressible navier-stokes) is 1 so in this case eps
    !! is an absolute value.
    norm = 1.0_rk

    if (params%eps_normalized) then
        ! normalizing is dependent on the variables you want to threshold!
        ! for customized variables you need custom normalization -> physics module
        norm = 0.0_rk
        block_norm = 0.0_rk

        if ( params%coarsening_indicator == "threshold-state-vector" ) then
            ! Apply thresholding directly to the statevector, not to derived quantities
            do k = 1, hvy_n
                call NORM_THRESHOLDFIELD_meta( params%physics_type, hvy_block(:,:,:,:,hvy_active(k)), block_norm)

                do p = 1, N_thresholding_components
                    norm(p) = max( norm(p), block_norm(p) )
                enddo
            end do
        else
            ! Apply thresholding to derived qtys, such as the vorticity
            do k = 1, hvy_n
                call NORM_THRESHOLDFIELD_meta( params%physics_type, hvy_tmp(:,:,:,:,hvy_active(k)), block_norm)

                do p = 1, N_thresholding_components
                    norm(p) = max( norm(p), block_norm(p) )
                enddo
            end do
        endif

        ! note that a check norm>1.0e-10 is in threshold-block
        tmp = norm
        call MPI_ALLREDUCE(tmp, norm, neq, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)

        ! during dev is it useful to know what the normalization is, if that is active
        call append_t_file('eps_norm.t', (/time, norm, params%eps/))
    endif


    !---------------------------------------------------------------------------
    !> evaluate coarsing criterion on all blocks
    !---------------------------------------------------------------------------
    ! the indicator "random" requires special treatment below (it does not pass via
    ! block_coarsening_indicator)
    if (indicator /= "random") then
        ! NOTE: even if additional mask thresholding is used, passing the mask is optional,
        ! notably because of the ghost nodes unit test, where random refinement / coarsening
        ! is used. hence, checking the flag params%threshold_mask alone is not enough.
        if (params%threshold_mask .and. present(hvy_mask)) then
            ! loop over all my blocks
            do k = 1, hvy_n
                hvy_id = hvy_active(k)
                ! get lgt id of block
                call hvy_id_to_lgt_id( lgt_id, hvy_id, params%rank, params%number_blocks )

                ! some indicators may depend on the grid, hence
                ! we pass the spacing and origin of the block (as we have to compute vorticity
                ! here, this can actually be omitted.)
                call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

                ! evaluate the criterion on this block.
                call block_coarsening_indicator( params, hvy_block(:,:,:,:,hvy_id), &
                hvy_tmp(:,:,:,:,hvy_id), dx, x0, indicator, iteration, &
                lgt_block(lgt_id, Jmax + IDX_REFINE_STS), norm,  hvy_mask(:,:,:,:,hvy_id))
            enddo
        else
            ! loop over all my blocks
            do k = 1, hvy_n
                hvy_id = hvy_active(k)
                ! get lgt id of block
                call hvy_id_to_lgt_id( lgt_id, hvy_id, params%rank, params%number_blocks )

                ! some indicators may depend on the grid, hence
                ! we pass the spacing and origin of the block (as we have to compute vorticity
                ! here, this can actually be omitted.)
                call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

                ! evaluate the criterion on this block.
                call block_coarsening_indicator( params, hvy_block(:,:,:,:,hvy_id), &
                hvy_tmp(:,:,:,:,hvy_id), dx, x0, indicator, iteration, &
                lgt_block(lgt_id, Jmax + IDX_REFINE_STS), norm)
            enddo
        endif
    else
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

                do k = 1, lgt_n
                    lgt_id = lgt_active(k)
                    ! random number
                    call random_number(r)
                    ! set refinement status to coarsen based on random numbers.
                    if ( r <= crsn_chance ) then
                        lgt_block(lgt_id, Jmax + IDX_REFINE_STS) = -1
                        tags = tags + 1
                    endif
                enddo
            endif
            ! sync light data, as only root sets random refinement
            call MPI_BCAST( lgt_block(:, Jmax + IDX_REFINE_STS), size(lgt_block,1), MPI_INTEGER4, 0, WABBIT_COMM, ierr )
        endif
    endif



    !---------------------------------------------------------------------------
    !> force blocks on maximum refinement level to coarsen, if parameter is set
    !---------------------------------------------------------------------------
    if (params%force_maxlevel_dealiasing) then
        do k = 1, lgt_n
            lgt_id = lgt_active(k)
            if (lgt_block(lgt_id, Jmax + IDX_MESH_LVL) == params%max_treelevel) then
                ! force blocks on maxlevel to coarsen
                lgt_block(lgt_id, Jmax + IDX_REFINE_STS) = -1
            endif
        enddo
    endif


    !> after modifying all refinement flags, we need to synchronize light data
    call synchronize_lgt_data( params, lgt_block, refinement_status_only=.true. )

end subroutine grid_coarsening_indicator
