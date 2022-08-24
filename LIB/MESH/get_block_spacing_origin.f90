!> \brief For any block lgt_id this routine computes, from the treecode stored in
!! lgt_block( lgt_id, : ), the block's origin and grid spacing. Note spacing
!! and origin are 3D vectors, the third component being zero in a 2D case. \n

subroutine get_block_spacing_origin( params, lgt_id, x0, dx )

  implicit none

  type (type_params), intent(in)             :: params              !> user defined parameter structure
  integer(kind=ik), intent(in)               :: lgt_id              !> the block in question
  real(kind=rk), dimension(1:3), intent(out) :: x0, dx              !> output
  integer(kind=ik)                           :: J

  J = lgt_block(lgt_id, params%max_treelevel+IDX_MESH_LVL)
  call get_block_spacing_origin2( lgt_block(lgt_id, 1:J), params%domain_size, params%Bs, params%dim, x0, dx )

end subroutine get_block_spacing_origin
