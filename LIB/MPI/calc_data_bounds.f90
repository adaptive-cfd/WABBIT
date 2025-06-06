logical function patch_crosses_periodic_BC(x0, dx, ijk, dim)
    implicit none
    real(kind=rk), intent(in) :: x0(1:3), dx(1:3)
    integer(kind=ik) , intent(in) :: ijk(2,3), dim
    real(kind=rk) :: x(1:8,1:3)

    ! corners of patch
    x(1,1:3) = x0 + (/ dx(1)*ijk(1,1), dx(2)*ijk(1,2), dx(3)*ijk(1,3) /)
    x(2,1:3) = x0 + (/ dx(1)*ijk(2,1), dx(2)*ijk(2,2), dx(3)*ijk(2,3) /)

    x(3,1:3) = x0 + (/ dx(1)*ijk(2,1), dx(2)*ijk(1,2), dx(3)*ijk(1,3) /)
    x(4,1:3) = x0 + (/ dx(1)*ijk(1,1), dx(2)*ijk(2,2), dx(3)*ijk(1,3) /)
    x(5,1:3) = x0 + (/ dx(1)*ijk(1,1), dx(2)*ijk(1,2), dx(3)*ijk(2,3) /)

    x(6,1:3) = x0 + (/ dx(1)*ijk(2,1), dx(2)*ijk(2,2), dx(3)*ijk(1,3) /)
    x(7,1:3) = x0 + (/ dx(1)*ijk(1,1), dx(2)*ijk(2,2), dx(3)*ijk(2,3) /)
    x(8,1:3) = x0 + (/ dx(1)*ijk(2,1), dx(2)*ijk(1,2), dx(3)*ijk(2,3) /)

    patch_crosses_periodic_BC = .false.
    ! note we normalize xl, yl, zl to 1.0 just in these routines here
    if ( minval(x(:,1:dim)) < 0.0_rk .or. maxval(x(:,1:dim)) > 1.0_rk) then
        patch_crosses_periodic_BC = .true.
    endif

end function


!> \brief Set send bounds and sender buffer for different neighborhood or patch relations
!> This uses get_indices_of_modify_patch for the different relations. For lvl_diff = -1 / +1 restriction or prediction changes
!! patch sizes. For this, the necessary buffers for the operators are computed (in data_bounds) and what will be used from those buffers
!! for sending (in data_buffer)
!
!> neighbor codes: \n
!  ---------------
!>   1- 56 : lvl_diff =  0  (same level)
!>  57-112 : lvl_diff = +1  (coarser neighbor)
!> 113-168 : lvl_diff = -1  (finer   neighbor)
!> For each range, the different 56 entries are:
!> 01-08 : X side (4-,4+)
!> 09-16 : Y-side (4-,4+)
!> 17-24 : Z-side (4-,4+)
!> 25-32 : X-Y edge (2--, 2+-, 2-+, 2++)
!> 33-40 : X-Z edge (2--, 2+-, 2-+, 2++)
!> 41-48 : Y-Z edge (2--, 2+-, 2-+, 2++)
!> 49-56 : corners (---, +--, -+-, ++-, --+, +-+, -++, +++)
subroutine set_send_bounds( params, data_bounds, data_buffer, relation, lvl_diff, gminus, gplus)
    implicit none

    type (type_params), intent(in)                  :: params
    !> data_bounds
    integer(kind=ik), intent(inout)                 :: data_bounds(2,3), data_buffer(2,3)
    !> neighborhood or family relation, id from dirs
    !! -8:-1 is mother/daughter relation
    !! 0 is full block relation, level
    !! 1:74 is neighborhood relation
    integer(kind=ik), intent(in)                    :: relation
    !> difference between block levels
    integer(kind=ik), intent(in)                    :: lvl_diff, gminus, gplus

    integer(kind=ik) :: n(1:3), g(3), Bs(1:3), i_dim, a, min_size, Nsender

    g(:) = params%g
    Bs = params%bs
    ! compute total size of array
    n(:) = 1
    do i_dim = 1,params%dim
        n(i_dim) = Bs(i_dim) + 2*g(i_dim)
    enddo

    ! the number a is how many extra coarse points on the sender side you use to
    ! avoid one-sided interpolation. The actual formula is S = (order-2 )/2 so for
    ! 6th order you would require 2 extra ones.
    if (params%order_predictor == "multiresolution_2nd" ) then
        a = 0
    elseif (params%order_predictor == "multiresolution_4th" ) then
        a = 1
    elseif (params%order_predictor == "multiresolution_6th" ) then
        a = 2
    else
        call abort(2875490, "The predictor method is unknown")
    endif

    ! set 1 and not -1 (or anything else), because 2D bounds ignore 3rd dimension
    ! and thus cycle from 1:1
    data_bounds(:,:) = 1
    data_buffer(:,:) = 1

    if (lvl_diff == 0) then
        call get_indices_of_modify_patch(params%g, params%dim, relation, data_bounds, n, &
            (/gminus, gminus, gminus/), (/gplus, gplus, gplus/), &
            g_m=(/g, g, g/), g_p=(/g, g, g/), lvl_diff=lvl_diff)
    elseif (lvl_diff == +1) then  ! restriction for coarser neighbor
        ! Restriction, example with Bs=8, g=3 for even grids:
        ! Normal grid:              g   g   g   i   i   i   i   i   i   i   i   g   g   g
        ! Restricted:                           l       lr      lr      r
        ! Values have to be odd (coincide with SC from spaghetti form), thats why left and right are treated differently then
        call get_indices_of_modify_patch(params%g, params%dim, relation, data_bounds, n, &
            (/gminus*2-1, gminus*2-1, gminus*2-1/), (/gplus*2-1, gplus*2-1, gplus*2-1/), &
            g_m=(/g, g, g/), g_p=(/g+1, g+1, g+1/), lvl_diff=lvl_diff)
        ! buffer is the same but 1-based
        do i_dim = 1, params%dim
            Nsender = data_bounds(2, i_dim) - data_bounds(1, i_dim) + 1
            data_buffer(2, i_dim) = (Nsender+1)/2
        enddo
    elseif (lvl_diff == -1) then ! prediction for finer neighbor
        ! Prediction, example with Bs=6, g=3, 2nd order for even grids:
        ! Normal grid:              g   g   g   i   i   i   i   i   i   g   g   g
        ! Predicted:                            l l l             r r r
        ! Needed for prediction:                r   r           r   r   r
        call get_indices_of_modify_patch(params%g, params%dim, relation, data_bounds, n, &
            (/gminus/2+1, gminus/2+1, gminus/2+1/), (/(gplus+1)/2+1, (gplus+1)/2+1, (gplus+1)/2+1/), &
            g_m=(/g, g, g/), g_p=(/g-1, g-1, g-1/), lvl_diff=lvl_diff)

        ! enlargen area with stencil size
        do i_dim = 1, params%dim
            data_bounds(1, i_dim) = data_bounds(1, i_dim) - a
            data_bounds(2, i_dim) = data_bounds(2, i_dim) + a
        enddo
        
        ! buffer, where we take the values out of it, those are left-shifted (first point is not interpolated, second one is)
        do i_dim = 1, params%dim
            ! first where do we start? If it is left-bounded then 1+a, else 2+a
            ! first point is always directly after ghost patch, also for even sync all other patches are left-bounded as well
            if (data_bounds(1, i_dim) == params%g+1-a .or. mod(gminus,2) == 0) then
                data_buffer(1, i_dim) = 1+a*2
            else
                data_buffer(1, i_dim) = 2+a*2
            endif
            data_buffer(2, i_dim) = data_buffer(1, i_dim) + gplus-1
            ! for edges the length is increased by Bs
            if (data_bounds(2, i_dim) - data_bounds(1, i_dim) > g(i_dim)+a) then
                data_buffer(2, i_dim) = data_buffer(2, i_dim) + Bs(i_dim)
            endif
        enddo
    endif
end subroutine set_send_bounds


!> \brief Set recv bounds for different neighborhood or patch relations.
!> This is simply a wrapper to get_indices_of_ghost_patch as there is no other special treatment that needs to be done.
!! ATTENTION: Relation is inverted for ghost_patch definition in comparison to interior patches due to sender-receiver logic
!
!> neighbor codes: \n
!  ---------------
!>   1- 56 : lvl_diff =  0  (same level)
!>  57-112 : lvl_diff = +1  (coarser neighbor)
!> 113-168 : lvl_diff = -1  (finer   neighbor)
!> For each range, the different 56 entries are:
!> 01-08 : X side (4-,4+)
!> 09-16 : Y-side (4-,4+)
!> 17-24 : Z-side (4-,4+)
!> 25-32 : X-Y edge (2--, 2+-, 2-+, 2++)
!> 33-40 : X-Z edge (2--, 2+-, 2-+, 2++)
!> 41-48 : Y-Z edge (2--, 2+-, 2-+, 2++)
!> 49-56 : corners (---, +--, -+-, ++-, --+, +-+, -++, +++)
subroutine set_recv_bounds( params, data_bounds, relation, lvl_diff, gminus, gplus)
    implicit none

    type (type_params), intent(in)                  :: params
    !> data_bounds
    integer(kind=ik), intent(inout)                 :: data_bounds(2,3)
    !> neighborhood or family relation, id from dirs
    !! -8:-1 is mother/daughter relation
    !! 0 is full block relation, level
    !! 1:56*3 is neighborhood relation
    integer(kind=ik), intent(in)                    :: relation
    !> difference between block levels
    integer(kind=ik), intent(in)                    :: lvl_diff, gminus, gplus

    integer(kind=ik) :: Bs(1:3), g, i_dim

    Bs = params%Bs
    g = params%g

    ! call function from wavelets
    call get_indices_of_ghost_patch(params%Bs, params%g, params%dim, relation, data_bounds, gminus, gplus, lvl_diff)

end subroutine set_recv_bounds