subroutine threshold_block( params, u, refinement_status, level, input_is_WD, norm, indices, eps, verbose_check)
    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> heavy data - this routine is called on one block only, not on the entire grid. hence th 4D array
    !! When input_is_WD is set they are expected to be already wavelet decomposed in Spaghetti-ordering
    real(kind=rk), intent(inout)        :: u(:, :, :, :)
    !> main output of this routine is the new satus
    integer(kind=ik), intent(out)       :: refinement_status
    !> If we use L2 or H1 normalization, the threshold eps is level-dependent, hence
    !! we pass the level to this routine
    integer(kind=ik), intent(in)        :: level
    logical, intent(in)                 :: input_is_WD                       !< flag if hvy_block is already wavelet decomposed
    real(kind=rk), intent(inout), optional  :: norm( size(u,4) )
    !> if different from the default eps (params%eps), you can pass a different value here. This is optional
    !! and used for example when thresholding the mask function.
    real(kind=rk), intent(in), optional :: eps( size(u,4) )
    !> Indices of patch if not the whole interior block should be thresholded, used for securityZone
    integer(kind=ik), intent(in), optional :: indices(1:2, 1:3)
    logical, intent(in), optional       :: verbose_check  !< No matter the value, if this is present we debug

    integer(kind=ik)                    :: dF, l, p_norm, p, idx(2,3), Bs(1:3), g, i_dim, dim, Jmax, nc
    real(kind=rk)                       :: detail( size(u,4) ), eps_use( size(u,4) )
    real(kind=rk), allocatable, dimension(:,:,:,:), save :: u_wc
    integer, dimension(:), allocatable  :: mask_i, mask

    nc     = size(u, 4)
    Bs     = params%Bs
    g      = params%g
    dim    = params%dim
    Jmax   = params%Jmax
    detail = -1.0_rk

    if (allocated(u_wc)) then
        if (size(u_wc, 4) < nc) deallocate(u_wc)
    endif
    ! ToDo: We should consider reducing the temporary work array size, we just need to take care how to do this with the mask
    if (.not. allocated(u_wc)) allocate(u_wc(1:size(u, 1), 1:size(u, 2), 1:size(u, 3), 1:nc ) )

    ! for threshold state vector components we do some slicing magic
    ! fortran does not allow logical mask slicing, so we build a mask of indices to slice with
    allocate(mask_i(1:nc))
    do l = 1, nc
        mask_i(l) = l
    end do

    ! set the indices we want to threshold
    idx(:, :) = 1
    if (present(indices)) then
        idx(:, :) = indices(:, :)
    else  ! full interior block
        idx(1, 1) = g+1
        idx(2, 1) = Bs(1)+g
        idx(1, 2) = g+1
        idx(2, 2) = Bs(2)+g
        if (dim == 3) then
            idx(1, 3) = g+1
            idx(2, 3) = Bs(3)+g
        endif
    endif

#ifdef DEV
    if (.not. allocated(params%GD)) call abort(1213149, "The cat is angry: Wavelet-setup not yet called?")
    ! if (modulo(Bs(1),2) /= 0) call abort(1213150, "The dog is angry: Block size must be even.")
    ! if (modulo(Bs(2),2) /= 0) call abort(1213150, "The dog is angry: Block size must be even.")
#endif

    if (.not. input_is_WD) then
        u_wc(:, :, :, 1:nc) = u(:, :, :, 1:nc)  ! Wavelet decompose full block
        call waveletDecomposition_block(params, u_wc(:, :, :, 1:nc)) ! data on u (WC/SC) now in Spaghetti order
        call wavelet_renorm_block(params, u_wc, u_wc, level, indices=indices, verbose_check=verbose_check)
    else
        call wavelet_renorm_block(params, u, u_wc, level, indices=indices, verbose_check=verbose_check)
    endif

    ! thresholding of some of the state vector components should be done together, as for example the velocity 
    ! components ux, uy, uz. This means: we do not threshold ux,uy,uz separately, but treat them as one vector.
    ! In this code, this is calle a norm-equivalent. The thresholding is guided by the array params%threshold_state_vector_component
    ! which is set in the params.ini file by the user. If the entry is >1, then this component belongs to a vector and is treated as such
    ! (as opposed to individual thresholding per component). For example, in ACM, we would use params%threshold_state_vector_component=(/2,2,2,1/)
    ! if we were to treat the velocity as one vector, and the pressure as a scalar. If we ever have more than one vector (which 
    ! is not the case currently, 12/2024), we can use integers to distinct between the vector fields. 
    ! So we could set 2 2 3 3 and actually have those treated as independent vector fields, where we compute the norm.
    ! maybe a bit overkill because usually we only have the velocity, but why not have the capacity that we can build up on
    do l = 2, maxval(params%threshold_state_vector_component(:))
        ! convert logical mask to index mask
        mask = pack(mask_i, params%threshold_state_vector_component(:)==l)
        detail(mask) = &
            maxval(sqrt(u_wc(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3):idx(2,3), mask)**2))
    enddo

    ! thresholding of scalar components of state vector (==1)
    do p = 1, nc
        if (params%threshold_state_vector_component(p) == 1) then
            ! if all details are smaller than C_eps, we can coarsen, check interior WC only
            detail(p) = maxval( abs(u_wc(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3):idx(2,3), p)) )
        endif
    enddo

    ! We could disable detail checking for qtys we do not want to consider, (in the sense of not even computing the WC),
    ! but this is more work and selective thresholding is rarely used.
    do p = 1, nc
        if (params%threshold_state_vector_component(p) == 0) detail(p) = 0.0_rk
    enddo

    ! default threshold is the one in the parameter struct
    eps_use(:) = params%eps
    ! but if we pass another one, use that.
    if (present(eps)) eps_use(1:nc) = eps(1:nc)

    ! Normalize threshold with the norm of the component. This corresponds to using a relative vs an absolute EPS.
    if (present(norm)) eps_use(1:nc) = eps_use(1:nc) * norm(1:nc)

    ! evaluate criterion: if this blocks detail is smaller than the prescribed precision,
    ! the block is tagged as "wants to coarsen" by setting the tag -1.
    ! Note gradedness and completeness may prevent it from actually going through with that.
    ! Note the most severe detail (= state vector component) sets the flag for the entire block.
    if ( all(detail(:) <= eps_use(:))) then
        ! coarsen block, -1
        refinement_status = -1
    else
        refinement_status = 0
    end if

    ! debugging info
    if (present(verbose_check) .and. any(detail(:) > eps_use(:) / 2.0_rk)) then
        write(*, '(A, es10.3, A, i2, A, 10(es10.3, 2x))') "Eps: ", eps_use, ", Ref stat: ", refinement_status, ", Details: ", detail(1:nc)
    endif
end subroutine threshold_block
