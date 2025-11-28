subroutine filter_wrapper(time, params, hvy_block, hvy_tmp, hvy_mask, tree_ID)
    implicit none

    real(kind=rk), intent(in)           :: time
    type (type_params), intent(inout)   :: params                       !> user defined parameter structure, hvy_active
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)     !> heavy data array - block data
    !> heavy temp data: used for saving, filtering, and helper qtys (reaction rate, mask function)
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    !> hvy_mask are qty that depend on the grid and not explicitly on time
    real(kind=rk), intent(inout)        :: hvy_mask(:, :, :, :, :)
    integer(kind=ik), intent(in)        :: tree_ID
    real(kind=rk), dimension(3)         :: dx, x0                       !> spacing and origin of a block
    integer(kind=ik)                    :: k, lgt_id, ic       ! loop variables
    integer(kind=ik)                    :: g, a, nc
    integer(kind=ik), dimension(3)      :: Bs
    integer(kind=2)                     :: n_domain(1:3)                ! surface normal
    integer(kind=ik)                    :: lvl, hvy_id
    real(kind=rk), allocatable          :: stencil(:)
    integer(kind=ik)                    :: stencil_size         ! filter position (array postion of value to filter)

    Bs = params%Bs
    g  = params%g
    n_domain = 0
    nc = size(hvy_block, 4)

    if (params%filter_only_maxlevel .and. params%filter_all_except_maxlevel) then
        call abort(251106,"ERROR: Do you want to filter only on max level or all except max level??? Choose one, not both.")
    end if

    select case(params%filter_type)
    case('explicit_3pt', 'superviscosity_2nd')
        call generate_superviscosity_stencil(stencil, a, 2)
    case('explicit_5pt', 'superviscosity_4th')
        call generate_superviscosity_stencil(stencil, a, 4)
        stencil = -stencil  ! filters are defined to be negative at center?
    case('explicit_7pt', 'superviscosity_6th')
        call generate_superviscosity_stencil(stencil, a, 6)
    case('explicit_9pt', 'superviscosity_8th')
        call generate_superviscosity_stencil(stencil, a, 8)
        stencil = -stencil  ! filters are defined to be negative at center?
    case('explicit_11pt', 'superviscosity_10th')
        call generate_superviscosity_stencil(stencil, a, 10)
    case('explicit_13pt', 'superviscosity_12th')
        call generate_superviscosity_stencil(stencil, a, 12)
        stencil = -stencil  ! filters are defined to be negative at center?
    case('explicit_15pt', 'superviscosity_14th')
        call generate_superviscosity_stencil(stencil, a, 14)
    case('explicit_17pt', 'superviscosity_16th')
        call generate_superviscosity_stencil(stencil, a, 16)
        stencil = -stencil  ! filters are defined to be negative at center?
    case('explicit_19pt', 'superviscosity_18th')
        call generate_superviscosity_stencil(stencil, a, 18)
    case('explicit_21pt', 'superviscosity_20th')
        call generate_superviscosity_stencil(stencil, a, 20)
        stencil = -stencil  ! filters are defined to be negative at center?
    case default
        call abort(251107,"ERROR: Filter not known: "//params%filter_type)
    end select

    if (a > params%g) then
        call abort(251108,"ERROR: Nice filter you've selected there, but its stencil size exceeds the ghost layer thickness. Increase g.")
    end if

    ! we sum the filter operation to the flow, so we need to add 1 at the center
    stencil(0) = stencil(0) + 1.0_rk

    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
        lvl = lgt_block(lgt_id, IDX_MESH_LVL)

        ! we have some special filter settings to filter only parts of the domain
        if (params%filter_only_maxlevel .and. lvl < params%Jmax) cycle
        if (params%filter_all_except_maxlevel .and. lvl == params%Jmax) cycle

        ! we filter each component separately, depending on user choice (defaults to all components)
        do ic = 1, params%n_eqn
            if (.not. params%filter_component(ic)) cycle
            call blockFilterXYZ_vct( params, hvy_block(:,:,:,ic:ic, hvy_id), hvy_block(:,:,:,ic:ic, hvy_id), stencil(-a:+a), -a, +a)
        enddo

    enddo
end subroutine filter_wrapper


! creates the explicit filter stencils, which can be used for superviscosity
! these have 2nd order accuracy and are of order 'order'
subroutine generate_superviscosity_stencil(stencil, length, order)
    real(kind=rk), intent(inout), allocatable :: stencil(:)
    integer(kind=ik), intent(in) :: order
    integer(kind=ik), intent(inout) :: length
    integer(kind=ik) :: k
    integer(kind=ik) :: binom

    length = order/2
    if (allocated(stencil)) then
        if (lbound(stencil, dim=1) > -length .or. ubound(stencil, dim=1) < length) then
            deallocate(stencil)
        end if
    end if
    if (.not. allocated(stencil)) allocate( stencil(-length:length) )

    do k = -length, length
        binom = binomial_coeff(length*2, length+k)
        stencil(k) = (-1.0d0)**(k+length) * dble(binom)
    end do
    ! normalize L1 norm to 1
    stencil = stencil / sum(abs(stencil))
end subroutine generate_superviscosity_stencil

! Function to compute binomial coefficient n choose k
integer function binomial_coeff(n, k)
    integer, intent(in) :: n, k
    integer :: i
    if (k < 0 .or. k > n) then
        binomial_coeff = 0
    else
        binomial_coeff = 1
        do i = 1, k
        binomial_coeff = binomial_coeff * (n - i + 1) / i
        end do
    end if
end function binomial_coeff