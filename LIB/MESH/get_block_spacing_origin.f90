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
  ! loop variables and shortcuts
  integer(kind=ik)                           :: ix,iy,iz,level,bs

  bs = params%number_block_nodes

  ! fetch this blocks level:
  level = lgt_block( lgt_id, params%max_treelevel + idx_mesh_lvl )

  ! compute its coordinates in ijk space
  call decoding( lgt_block( lgt_id, 1:level ), ix,iy,iz, level)

  ! the spacing on a block is the basic spacing Lx/Bs of the coarsest block (if there
  ! is only one block, j=0) divided by 2 for each level, thus the 2^-j factor
  dx = 2.0_rk**(-level) * (/params%domain_size(1) , params%domain_size(2) , params%domain_size(3)/) / real(bs-1, kind=rk)
  ! note zero based indexing:
  x0 = real( ((/ix,iy,iz/) - 1)*(Bs-1) ,kind=rk) * dx

end subroutine get_block_spacing_origin
