!> \brief For any block lgt_id this routine computes, from the treecode stored in
!! lgt_block( lgt_id, : ), the block's origin and grid spacing. Note spacing
!! and origin are 3D vectors, the third component being zero in a 2D case. \n

! Function acts as a wrapper, as it is used often we want to shorten calls
subroutine get_block_spacing_origin( params, lgt_id, x0, dx )

  implicit none

  type (type_params), intent(in)             :: params              !> user defined parameter structure
  integer(kind=ik), intent(in)               :: lgt_id              !> the block in question
  real(kind=rk), dimension(1:3), intent(out) :: x0, dx              !> output
  integer(kind=tsize)                        :: treecode
  integer(kind=ik)                           :: J

  call get_block_spacing_origin_b( get_tc(lgt_block(lgt_id, params%Jmax+IDX_TC_1 : params%Jmax+IDX_TC_2)), params%domain_size, &
    params%Bs, x0, dx, dim=params%dim, level=lgt_block(lgt_id, params%Jmax+IDX_MESH_LVL), max_level=params%Jmax)

end subroutine get_block_spacing_origin


!-----------------------------------------------------------------------------
!> \brief Computes the surface normal of the global domain boundary, if the current block is adjacent to the boundary.\n
!> \details
!>   - The surface normal is 0 if the block is not adjacent to the global domain boundary (even if the BC is periodic!!)
!>   - The surface normal is computed from the treecode
! This function is here as it basically is very similar to the upper one
subroutine get_adjacent_boundary_surface_normal(params, lgt_id, n_surface)
  implicit none
  type (type_params), intent(in)             :: params              !> user defined parameter structure
  integer(kind=ik), intent(in)               :: lgt_id              !> the block in question
  ! The normal on the domain indicates (if non-periodic BC are used), if a block
  ! is at the outer, usually periodic border of the domain ( x,y,z == 0 and x,y,z == L)
  ! Nonzero values indicate this is the case, e.g., n_domain=(/1, 0, -1/) means in x-axis, our block
  ! is way at the back and its boundary normal points in +x, and in z, its at the bottom (z=0), thus
  ! its normal points downwards.
  integer(kind=2), intent(out) :: n_surface(3)

  real(kind=rk), dimension(1:3) :: x0, dx
  real(kind=rk) :: tolerance
  integer(kind=ik) :: i_dim

  call get_block_spacing_origin_b( get_tc(lgt_block(lgt_id, params%Jmax+IDX_TC_1 : params%Jmax+IDX_TC_2)), params%domain_size, &
    params%Bs, x0, dx, dim=params%dim, level=lgt_block(lgt_id, params%Jmax+IDX_MESH_LVL), max_level=params%Jmax)

  tolerance = 1.0e-3_rk * minval(dx(1:params%dim))

  n_surface(1:params%dim) = 0
  do i_dim = 1, params%dim
    if (abs(x0(i_dim)-0.0_rk) < tolerance ) then !x_i == 0
      n_surface(i_dim) = -1
    elseif (abs(x0(i_dim)+dx(i_dim)*real(params%Bs(i_dim)-1,kind=rk) - params%domain_size(i_dim)) < tolerance) then ! x_i == L
      n_surface(i_dim) = +1
    endif
  end do

end subroutine get_adjacent_boundary_surface_normal