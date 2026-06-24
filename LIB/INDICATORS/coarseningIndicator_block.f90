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

subroutine coarseningIndicator_block( params, block_data, block_work, indicator, &
    refinement_status, norm, level, input_is_WD, indices, verbose_check)
    ! it is not technically required to include the module here, but for VS code it reduces the number of wrong "errors"
    use module_params

    implicit none
    type (type_params), intent(in)      :: params
    !> heavy data - this routine is called on one block only, not on the entire grid. hence th 4D array.
    real(kind=rk), intent(inout)        :: block_data(:, :, :, :)
    !> heavy work data array (expected to hold the VORTICITY if thresholding is applied to vorticity)
    real(kind=rk), intent(inout)        :: block_work(:, :, :, :)
    !> how to choose blocks for refinement
    character(len=*), intent(in)        :: indicator
    ! If we use L2 or H1 normalization, the threshold eps is level-dependent, hence
    ! we pass the level to this routine
    integer(kind=ik), intent(in)        :: level
    !> output is the refinement_status
    integer(kind=ik), intent(out)       :: refinement_status
    !
    real(kind=rk), intent(inout)        :: norm(1:size(block_data,4))
    logical, intent(in)                 :: input_is_WD                       !< flag if hvy_block is already wavelet decomposed
    !> Indices of patch if not the whole interior block should be thresholded, used for securityZone
    integer(kind=ik), intent(in), optional :: indices(1:2, 1:3)
    logical, intent(in), optional       :: verbose_check  !< No matter the value, if this is present we debug

    ! local variables
    integer(kind=ik) :: k, Jmax, d, j, hvy_id, g, refinement_status_mask, tags, ix, iy, iz, idx(2,3)
    integer(kind=ik), dimension(3) :: Bs
    ! chance for block refinement, random number
    real(kind=rk) :: crsn_chance, r, mask_max, mask_min
    real(kind=rk) :: t0  !< timing for debugging

    Jmax = params%Jmax
    Bs = params%Bs
    g = params%g

    ! set the indices we want to threshold
    idx(:, :) = 1
    if (present(indices)) then
        idx(:, :) = indices(:, :)
    else  ! full interior block
        idx(1, 1) = g+1
        idx(2, 1) = Bs(1)+g
        idx(1, 2) = g+1
        idx(2, 2) = Bs(2)+g
        if (params%dim == 3) then
            idx(1, 3) = g+1
            idx(2, 3) = Bs(3)+g
        endif
    endif

    !> This routine sets the -1 coarsening flag on a block. it uses different methods to
    !! decide where to coarsen, each acts on one block. Note due to gradedness and completeness
    !! this status may be revoked later in the computation.
    select case (indicator)
    case ("maxval-eps")
        ! debug indicator (useful for grid generation tests)
        ! coarsen a block if the maxval of its first component is smaller 0.9 (for passive scalars)
        ! merge selects 2D or 3D bounds depending on params%dim
        if ( maxval(block_data(g+1:Bs(1)+g, g+1:Bs(2)+g, merge(1, 1+g, params%dim == 2):merge(1, Bs(3)+g, params%dim == 2), 1) ) <= 0.9_rk ) then
            refinement_status = -1
        endif

    case ("mask-allzero-noghosts")
        ! merge selects 2D or 3D bounds depending on params%dim
        if ( maxval(block_data(g+1:Bs(1)+g, g+1:Bs(2)+g, merge(1, 1+g, params%dim == 2):merge(1, Bs(3)+g, params%dim == 2), 1) )<1.0e-9_rk ) then
            refinement_status = -1
        endif
        if ( minval(block_data(g+1:Bs(1)+g, g+1:Bs(2)+g, merge(1, 1+g, params%dim == 2):merge(1, Bs(3)+g, params%dim == 2), 1) )>=1.0_rk- 1.0e-9_rk ) then
            refinement_status = -1
        endif

    case ("threshold-state-vector", "threshold-cvs", "threshold-image-denoise", "primary-variables")
        ! t0 = MPI_Wtime()
        !! use wavelet indicator to check where to coarsen. Note here, active components are considered
        !! and the max over all active components results in the coarsening state -1. The components
        !! to be used can be specified in the PARAMS file. default is all componants.
#ifdef DEV
        if (.not. allocated(params%threshold_state_vector_component)) then
            call abort(7363823, "params%threshold_state_vector_component not allocated....")
        endif
#endif

        if (indicator == "threshold-cvs" .or. indicator == "threshold-image-denoise") then
            call threshold_block( params, block_data, refinement_status, level, input_is_WD, eps=norm, indices=idx, verbose_check=verbose_check)
        else
            call threshold_block( params, block_data, refinement_status, level, input_is_WD, norm=norm, indices=idx, verbose_check=verbose_check)
        endif

        ! timing for debugging - block based so should not be deployed for productive versions
        ! call toc( "coarseningIndicator_block (threshold_block)", 1000, MPI_Wtime()-t0 )
    case default
        call abort(87455214,"ERROR: unknown coarsening operator: "//trim(adjustl(indicator)))

    end select

end subroutine coarseningIndicator_block
