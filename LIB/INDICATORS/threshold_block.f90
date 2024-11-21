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

    ! we need to know if the first point is a SC or WC for the patch we check and skip it if it is a WC
    !     1 2 3 4 5 6 7 8 9 A B C
    !     G G G S W S W S W G G G
    !                 I I I
    ! 1-C - index numbering in hex format, G - ghost point, S - SC, W - WC, I - point of patch to be checked
    ! Patch I is checked, but we need to know that index 7 has a WC and should be skipped
    ! this is for parity with inflatedMallat version where SC and WC are situated on the SC indices of the spaghetti format
    ! for g=odd, the SC are on even numbers; for g=even, the SC are on odd numbers
    idx(1, 1:params%dim) = idx(1, 1:params%dim) + modulo(g + idx(1, 1:params%dim) + 1, 2)
    ! also, when the last point is a SC, the last point is only partially included but its WC have to be considered
    ! this gives problem if the last point is a SC so we need to handle this special case
    do i_dim = 1, params%dim
        if (idx(2, i_dim) /= size(u, i_dim)) then
            idx(2, i_dim) = idx(2, i_dim) + modulo(idx(2, i_dim) - idx(1, i_dim) + 1, 2)
        endif
    enddo


#ifdef DEV
    if (.not. allocated(params%GD)) call abort(1213149, "The cat is angry: Wavelet-setup not yet called?")
    if (modulo(Bs(1),2) /= 0) call abort(1213150, "The dog is angry: Block size must be even.")
    if (modulo(Bs(2),2) /= 0) call abort(1213150, "The dog is angry: Block size must be even.")
#endif

    if (.not. input_is_WD) then
        u_wc(:, :, :, 1:nc) = u(:, :, :, 1:nc)  ! Wavelet decompose full block
        call waveletDecomposition_block(params, u_wc(:, :, :, 1:nc)) ! data on u (WC/SC) now in Spaghetti order
        call wavelet_renorm_block(params, u_wc, u_wc, level, indices=indices, verbose_check=verbose_check)
    else
        call wavelet_renorm_block(params, u, u_wc, level, indices=indices, verbose_check=verbose_check)
    endif

    ! some thresholding components want to be treated together, as for example the velocity, so we compute a norm-equivalent if wanted
    ! every value >1 in thresholding_component are treated as belonging together
    ! so we could set 2 2 3 3 and actually have those treated as independent vector fields, where we compute the norm
    ! maybe a bit overkill because usually we only have the velocity, but why not have the capacity that we can build up on
    do l = 2, maxval(params%threshold_state_vector_component(:))
        ! convert logical mask to index mask
        mask = pack(mask_i, params%threshold_state_vector_component(:)==l)
        detail(mask) = &
            maxval(sqrt(u_wc(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3):idx(2,3), mask)**2))

        ! p_norm = 0
        ! do p = 1, nc
        !     ! set value onto which we construct the norm
        !     if (params%threshold_state_vector_component(p) == l .and. p_norm == 0) then
        !         p_norm = p
        !         ! first value is set
        !         u_wc(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3):idx(2,3), p_norm) = &
        !         u_wc(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3):idx(2,3), p_norm) **2
        !     else
        !         ! all other values are added on top to where we construct the norm
        !         u_wc(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3):idx(2,3), p_norm) = &
        !         u_wc(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3):idx(2,3), p_norm) + &
        !         u_wc(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3):idx(2,3), p)**2
        !     endif
        ! enddo
        ! ! now compute the square root
        ! u_wc(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3):idx(2,3), p_norm) = sqrt(u_wc(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3):idx(2,3), p_norm))
    enddo

    do p = 1, nc
        if (params%threshold_state_vector_component(p) == 1) then
            ! if all details are smaller than C_eps, we can coarsen, check interior WC only
            detail(p) = maxval( abs(u_wc(idx(1,1):idx(2,1), idx(1,2):idx(2,2), idx(1,3):idx(2,3), p)) )
        endif
    enddo

    ! We could disable detail checking for qtys we do not want to consider,
    ! but this is more work and selective thresholding is rarely used
    do p = 1, nc
        if (params%threshold_state_vector_component(p) == 0) detail(p) = 0.0_rk
    enddo

    ! ! for qtys where we compute norm values we have to disable the rest
    ! if ((params%eps_norm == "L2" .or. params%eps_norm == "H1") .and. params%physics_type == "ACM_new") then
    !     do l = 2, maxval(params%threshold_state_vector_component(:))
    !         p_norm = 0
    !         do p = 1, nc
    !             ! the first value is where we constructed the norm
    !             if (params%threshold_state_vector_component(p) == l .and. p_norm == 0) then
    !                 p_norm = p
    !             ! all others get their details deleted
    !             else
    !                 detail(p) = 0.0_rk
    !             endif
    !         enddo
    !     enddo
    ! endif

    ! default thresholding level is the one in the parameter struct
    eps_use(:) = params%eps
    ! but if we pass another one, use that.
    if (present(eps)) eps_use(1:nc) = eps(1:nc)

    ! Renorm criteria with norm of the field
    if (present(norm)) eps_use(1:nc) = eps_use(1:nc) * norm(1:nc)

    ! evaluate criterion: if this blocks detail is smaller than the prescribed precision,
    ! the block is tagged as "wants to coarsen" by setting the tag -1
    ! note gradedness and completeness may prevent it from actually going through with that
    if ( all(detail(:) < eps_use(:))) then
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
