!-----------------------------------------------------------------------------
! main level wrapper to filter a block. Note this is completely
! independent of the grid and any MPI formalism, neighboring relations and the like.
! You just get a block data (e.g. ux, uy, uz, p) and apply your filter to it.
! Ghost nodes are assumed to be sync'ed.
!-----------------------------------------------------------------------------
subroutine filter_ACM( time, u, g, x0, dx, work_array )
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

  ! output. Note assumed-shape arrays
  real(kind=rk), intent(inout) :: work_array(1:,1:,1:,1:)

  ! local variables
  integer(kind=ik) :: i
  integer(kind=ik), dimension(3) :: Bs

  ! compute the size of blocks
  Bs(1) = size(u,1) - 2*g
  Bs(2) = size(u,2) - 2*g
  Bs(3) = size(u,3) - 2*g

  select case (params_acm%filter_type)
  case ("no_filter")
      ! do nothing (this routine should not be called if "no_filter" is set)

  case ("wavelet_filter")
      ! apply wavelet filter to all components of state vector
      ! NOTE: this is a filter that removes all details.
      do i = 1, size(u,4)
          call wavelet_filter(params_acm%order_predictor, Bs, g, u(:,:,:,i))
      enddo

  end select


end subroutine




!=====================================================================
!  WAVELET FILTER
!=====================================================================
subroutine wavelet_filter(order_predictor, Bs, g, block_data)
    use module_interpolation, only :  restriction_3D, restriction_2D, prediction_2D, prediction_3D

    implicit none
    !> params structure of navier stokes
    character(len=*), intent(in) :: order_predictor
    !> mesh params
    integer(kind=ik), dimension(3), intent(in) :: Bs
    integer(kind=ik), intent(in) :: g
    !> heavy data array - block data
    real(kind=rk), intent(inout) :: block_data(:, :, :)
    real(kind=rk), allocatable, save :: u3(:,:,:)


    if ( size(block_data,3)>1 ) then
        ! ********** 3D **********
        if (.not.allocated(u3)) allocate(u3((Bs(1)+1)/2+g,(Bs(2)+1)/2+g,(Bs(3)+1)/2+g))
        ! now, coarsen array u1 (restriction)
        call restriction_3D( block_data, u3 )  ! fine, coarse
        ! then, re-interpolate to the initial level (prediciton)
        call prediction_3D ( u3, block_data, order_predictor )  ! coarse, fine

    else
        ! ********** 2D **********
        if (.not.allocated(u3)) allocate(u3((Bs(1)+1)/2+g,(Bs(2)+1)/2+g,1))
        ! now, coarsen array u1 (restriction)
        call restriction_2D( block_data(:,:,1), u3(:,:,1) )  ! fine, coarse
        ! then, re-interpolate to the initial level (prediciton)
        call prediction_2D ( u3(:,:,1), block_data(:,:,1), order_predictor )  ! coarse, fine

    end if

end subroutine wavelet_filter
