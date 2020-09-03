!> \file
!> \author engels
!> \brief For any block lgt_id this routine computes, from the treecode stored in
!! lgt_block( lgt_id, : ), the block's origin and grid spacing. Note spacing
!! and origin are 3D vectors, the third component being zero in a 2D case. \n
!
!! \details
!! \date 30/03/17 - create
!
subroutine get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

  implicit none

  !> user defined parameter structure
  type (type_params), intent(in)             :: params
  !> the block in question
  integer(kind=ik), intent(in)               :: lgt_id
  !> output
  real(kind=rk), dimension(1:3), intent(out) :: x0, dx
  !> light data array
  integer(kind=ik), intent(in)               :: lgt_block(:, :)

  integer(kind=ik) :: J

  J = lgt_block(lgt_id, params%max_treelevel+IDX_MESH_LVL)
  call get_block_spacing_origin2( lgt_block(lgt_id, 1:J), params%domain_size, params%Bs, params%dim, x0, dx )

end subroutine get_block_spacing_origin
