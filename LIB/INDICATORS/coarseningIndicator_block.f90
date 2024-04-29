!> \details This routine sets the coarsening flag for a single blocks. We allow for different
!! mathematical methods (everywhere / random) currently not very complex, but expected to grow
!! in the future.
!! ------------------ \n
!! Refinement status: \n
!! ------------------ \n
!! +1 refine \n
!! 0 do nothing \n
!! -1 block wants to coarsen \n
! ********************************************************************************************

subroutine coarseningIndicator_block( params, block_data, block_work, dx, x0, indicator, &
    iteration, refinement_status, norm, level, detail_precomputed, block_mask)

    implicit none
    type (type_params), intent(in)      :: params
    !> heavy data - this routine is called on one block only, not on the entire grid. hence th 4D array.
    real(kind=rk), intent(inout)        :: block_data(:, :, :, :)
    ! mask data. we can use different trees (4est module) to generate time-dependent/indenpedent
    ! mask functions separately. This makes the mask routines tree-level routines (and no longer
    ! block level) so the physics modules have to provide an interface to create the mask at a tree
    ! level. All parts of the mask shall be included: chi, boundary values, sponges.
    ! On input, the mask array is correctly filled. You cannot create the full mask here.
    ! NOTE: Here, the mask is required only if grid adaptation is also done on the mask.
    real(kind=rk), intent(inout), optional :: block_mask(:, :, :, :)
    !> heavy work data array (expected to hold the VORTICITY if thresholding is applied to vorticity)
    real(kind=rk), intent(inout)        :: block_work(:, :, :, :)
    !> block spacing and origin
    real(kind=rk), intent(in)           :: dx(1:3), x0(1:3)
    !> how to choose blocks for refinement
    character(len=*), intent(in)        :: indicator
    !> coarsening iteration index. coarsening is done until the grid has reached
    !! the steady state; therefore, this routine is called several times for
    integer(kind=ik), intent(in)        :: iteration
    ! If we use L2 or H1 normalization, the threshold eps is level-dependent, hence
    ! we pass the level to this routine
    integer(kind=ik), intent(in)        :: level
    !> output is the refinement_status
    integer(kind=ik), intent(out)       :: refinement_status
    !
    real(kind=rk), intent(inout)        :: norm(1:size(block_data,4))
    real(kind=rk), intent(inout)        :: detail_precomputed(:)

    ! local variables
    integer(kind=ik) :: k, Jmax, d, j, hvy_id, g, refinement_status_mask, tags, ix, iy, iz
    integer(kind=ik), dimension(3) :: Bs
    ! chance for block refinement, random number
    real(kind=rk) :: crsn_chance, r, mask_max, mask_min
    logical :: thresholding_component(1:size(block_data,4))
    real(kind=rk) :: t0  !< timing for debugging

    Jmax = params%Jmax
    Bs = params%Bs
    g = params%g

    !> This routine sets the -1 coarsening flag on a block. it uses different methods to
    !! decide where to coarsen, each acts on one block. Note due to gradedness and completeness
    !! this status may be revoked later in the computation.
    select case (indicator)
    case ("maxval-eps")
        ! debug indicator (useful for grid generation tests)
        ! coarsen a block if the maxval of its first component is smaller 0.9 (for passive scalars)
        if (params%dim==2) then
            if ( maxval(block_data(g+1:Bs(1)+g, g+1:Bs(2)+g, 1, 1) ) <= 0.9_rk ) then
                refinement_status = -1
            endif
        else
            if ( maxval(block_data(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 1) ) <= 0.9_rk ) then
                refinement_status = -1
            endif
        endif

    case ("mask-allzero-noghosts")
        if ( maxval(block_data(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 1) )<1.0e-9_rk ) then
            refinement_status = -1
        endif
        if ( minval(block_data(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 1) )>=1.0_rk- 1.0e-9_rk ) then
            refinement_status = -1
        endif

    case ("threshold-state-vector", "primary-variables")
        t0 = MPI_Wtime()
        !! use wavelet indicator to check where to coarsen. Note here, active components are considered
        !! and the max over all active components results in the coarsening state -1. The components
        !! to be used can be specified in the PARAMS file. default is all componants.
#ifdef DEV
        if (.not. allocated(params%threshold_state_vector_component)) then
            call abort(7363823, "params%threshold_state_vector_component not allocated....")
        endif
#endif

        thresholding_component = params%threshold_state_vector_component
        call threshold_block( params, block_data, thresholding_component, refinement_status, norm, level, detail_precomputed )

        call toc( "coarseningIndicator_block (treshold_block)", MPI_Wtime()-t0 )
    case default
        call abort(151413,"ERROR: unknown coarsening operator: "//trim(adjustl(indicator)))

    end select


    ! mask thresholding on top of regular thresholding?
    ! it can be useful to also use the mask function (if penalization is used) for grid adaptation.
    ! i.e. the grid is always at the finest level on mask interfaces. Careful though: the Penalization
    ! is implemented on physics-module level, i.e. it is not available for all modules.  If it is
    ! not available, the option is useless but can cause errors.
    ! NOTE: since the CDF44 wavelet is expensive, we use an alternative method to detect the gradient.
    if (params%threshold_mask .and. present(block_mask)) then
        t0 = MPI_Wtime()
        ! even if the global eps is very large, we want the fluid/solid (mask interface) to be on the finest level
        refinement_status_mask = -1_ik ! default we coarsen
        mask_max = 0.0_rk
        mask_min = 2.0_rk

        if (params%dim == 3) then
            do iz = g+1, Bs(3)+g
                do iy = g+1, Bs(2)+g
                    do ix = g+1, Bs(1)+g
                        if ((block_mask(ix,iy,iz,1) > 1.0e-9_rk).and.(block_mask(ix,iy,iz,1) < 1.0_rk-1.0e-9_rk)) then
                            ! but if we find intermediate values on the block (i.e. its neither uniformily
                            ! 0 nor 1), the block contains the fluid/solid interface and is not coarsened.
                            refinement_status_mask = 0_ik
                        endif
                        mask_max = max(mask_max, block_mask(ix,iy,iz,1))
                        mask_min = min(mask_min, block_mask(ix,iy,iz,1))
                    enddo
                enddo
            enddo
        else
            do iy = g+1, Bs(2)+g
                do ix = g+1, Bs(1)+g
                    if ((block_mask(ix,iy,1,1) > 1.0e-6_rk).and.(block_mask(ix,iy,1,1) < 1.0_rk-1.0e-6_rk)) then
                        ! but if we find intermediate values on the block (i.e. its neither uniformily
                        ! 0 nor 1), the block contains the fluid/solid interface and is not coarsened.
                        refinement_status_mask = 0_ik
                    endif
                    mask_max = max(mask_max, block_mask(ix,iy,1,1))
                    mask_min = min(mask_min, block_mask(ix,iy,1,1))
                enddo
            enddo
        endif

        ! maybe the resolution is so coarse no point on the smoothing layer exists...
        if ((mask_max-mask_min)>1.0e-6) refinement_status_mask = 0_ik

        ! refinement_status_state: -1 refinement_status_mask: -1 ==>  -1
        ! refinement_status_state: 0  refinement_status_mask: -1 ==>   0
        ! refinement_status_state: 0  refinement_status_mask: 0  ==>   0
        refinement_status = max(refinement_status, refinement_status_mask)

        call toc( "coarseningIndicator_block (mask_comp)", MPI_Wtime()-t0 )
    endif

end subroutine coarseningIndicator_block
