! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: get_block_spacing_origin.f90
! version: 0.5
! author: engels
!
! For any block lgt_id this routine computes, from the treecode stored in
! lgt_block( lgt_id, : ), the block's origin and grid spacing. Note spacing
! and origin are 3D vectors, the third component being zero in a 2D case.
!
! input:    params
!           lgt_id (this is the block we look at)
!           lgt_block (this is the list of light block data, which contains the treecode, level (and refinement status, but we don't use that here))
! output:   x0(1:3) vector with the origin of the block
!           dx(1:3) spacing on the block
!
! = log ======================================================================================
!
! 30/03/17 - create
!
subroutine get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )
  implicit none

  ! user defined parameter structure
  type (type_params), intent(in)      :: params
  ! the block in question
  integer(kind=ik)                    :: lgt_id
  ! output
  real(kind=rk), dimension(1:3), intent(out) :: x0, dx
  ! light data array
  integer(kind=ik), intent(inout)     :: lgt_block(:, :)

  integer(kind=ik) :: ix,iy,iz,level,bs

  x0 = 0.0_rk
  dx = 0.0_rk

  bs = params%number_block_nodes

  ! fetch this blocks level:
  level = lgt_block( lgt_id, params%max_treelevel+1 )

  ! compute its coordinates in ijk space
  call decoding( lgt_block( lgt_id, 1:level ), ix,iy,iz, level)

  ! the spacing on a block is the basic spacing Lx/Bs of the coarsest block (if there
  ! is only one block, j=0) divided by 2 for each level, thus the 2^-j factor
  dx = 2.0_rk**(-level) * (/params%Lx , params%Ly , params%Lz/) / real(bs-1, kind=rk)
  ! note zero based indexing:
  x0 = real( ((/ix,iy,iz/) - 1)*(Bs-1) ,kind=rk) * dx
end subroutine get_block_spacing_origin
