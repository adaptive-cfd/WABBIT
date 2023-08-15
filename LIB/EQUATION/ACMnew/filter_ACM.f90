!-----------------------------------------------------------------------------
! main level wrapper to filter a block. Note this is completely
! independent of the grid and any MPI formalism, neighboring relations and the like.
! You just get a block data (e.g. ux, uy, uz, p) and apply your filter to it.
! Ghost nodes are assumed to be sync'ed.
!-----------------------------------------------------------------------------
subroutine filter_ACM( time, u, g, x0, dx, work_array, mask )
  implicit none
  ! it may happen that some source terms have an explicit time-dependency
  ! therefore the general call has to pass time
  real(kind=rk), intent (in) :: time

  ! block data, containg the state vector. In general a 4D field (3 dims+components)
  ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
  real(kind=rk), intent(inout) :: u(1:,1:,1:,1:)

  ! as you are allowed to compute the work_array only in the interior of the field
  ! you also need to know where 'interior' starts: so we pass the number of ghost points
  integer, intent(in) :: g

  ! for each block, you'll need to know where it lies in physical space. The first
  ! non-ghost point has the coordinate x0, from then on its just cartesian with dx spacing
  real(kind=rk), intent(in) :: x0(1:3), dx(1:3)

  ! work array. Note assumed-shape arrays
  real(kind=rk), intent(inout) :: work_array(1:,1:,1:,1:)

  ! penalization mask function
  real(kind=rk), intent(inout) :: mask(1:,1:,1:,1:)

  ! local variables
  integer(kind=ik) :: i
  integer(kind=ik), dimension(3) :: Bs

  ! stencil array, note: size is fixed
  real(kind=rk) :: stencil(19) = 0.0_rk
  ! filter position (array postion of value to filter)
  integer(kind=ik) :: stencil_size
  integer :: j,l
  ! filtered values and array for old block data
  real(kind=rk) :: phi_tilde(3)
  integer :: dF

  ! compute the size of blocks
  Bs(1) = size(u,1) - 2*g
  Bs(2) = size(u,2) - 2*g
  Bs(3) = size(u,3) - 2*g

  if (.not. params_acm%initialized) write(*,*) "WARNING: filter_ACM called but ACM not initialized"

end subroutine
