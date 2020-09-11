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


  select case (params_acm%filter_type)
  case ("explicit_7pt")
      ! stencil_size = 7
     ! stencil(1:stencil_size) = (/ 1.0_rk/64.0_rk, -3.0_rk/32.0_rk, 15.0_rk/64.0_rk, -5.0_rk/16.0_rk, 15.0_rk/64.0_rk, -3.0_rk/32.0_rk, 1.0_rk/64.0_rk/)
     stencil_size = 9
     stencil(1:stencil_size) = (/ -1.0_rk/256.0_rk, &
                                   1.0_rk/ 32.0_rk, &
                                  -7.0_rk/ 64.0_rk, &
                                   7.0_rk/ 32.0_rk, &
                                 -35.0_rk/128.0_rk, &
                                   7.0_rk/ 32.0_rk, &
                                  -7.0_rk/ 64.0_rk, &
                                   1.0_rk/ 32.0_rk, &
                                  -1.0_rk/256.0_rk/)

      ! stencil_size = 11
      ! stencil(1:stencil_size) = (/  1.0_rk/1024.0_rk, &
      !                              -5.0_rk/ 512.0_rk, &
      !                              45.0_rk/1024.0_rk, &
      !                             -15.0_rk/ 128.0_rk, &
      !                             105.0_rk/ 512.0_rk, &
      !                             -63.0_rk/ 256.0_rk, &
      !                             105.0_rk/ 512.0_rk, &
      !                             -15.0_rk/ 128.0_rk, &
      !                              45.0_rk/1024.0_rk, &
      !                              -5.0_rk/ 512.0_rk, &
      !                               1.0_rk/1024.0_rk/)


      work_array(:,:,:,1:params_acm%dim) = u(:,:,:,1:params_acm%dim)

      do dF = 1, params_acm%dim
          ! 3D or 2D case
          if (params_acm%dim==3 ) then
              ! 3D
              do i = g+1, Bs(1)+g
                  do j = g+1, Bs(2)+g
                      do l = g+1, Bs(3)+g
                          ! x direction
                          call filter_1D( work_array(i-( (stencil_size+1)/2-1):i+( (stencil_size+1)/2-1), j, l, dF ), phi_tilde(1), stencil(1:stencil_size) )
                          ! y direction
                          call filter_1D( work_array(i, j-( (stencil_size+1)/2-1):j+( (stencil_size+1)/2-1), l, dF ), phi_tilde(2), stencil(1:stencil_size) )
                          ! z direction
                          call filter_1D( work_array(i, j, l-( (stencil_size+1)/2-1):l+( (stencil_size+1)/2-1), dF ), phi_tilde(3), stencil(1:stencil_size) )
                          ! filter
                          work_array(i, j, l, dF ) = u(i, j, l, dF ) + phi_tilde(1) + phi_tilde(2) + phi_tilde(3)
                      end do
                  end do
              end do
          else
              ! 2D
              do i = g+1, Bs(1)+g
                  do j = g+1, Bs(2)+g
                      ! x direction
                      call filter_1D( work_array(i-( (stencil_size+1)/2-1):i+( (stencil_size+1)/2-1), j, 1, dF ), phi_tilde(1), stencil(1:stencil_size) )
                      ! y direction
                      call filter_1D( work_array(i, j-( (stencil_size+1)/2-1):j+( (stencil_size+1)/2-1), 1, dF ), phi_tilde(2), stencil(1:stencil_size) )
                      ! filter
                      work_array(i, j, :, dF ) = u(i, j, :, dF) + phi_tilde(1) + phi_tilde(2)
                  end do
              end do

          endif
      end do

      u(:,:,:,1:params_acm%dim) = work_array(:,:,:,1:params_acm%dim)

  case ("no_filter")
      ! do nothing (this routine should not be called if "no_filter" is set)


  case ("wavelet_filter")
      ! apply wavelet filter to all components of state vector
      ! NOTE: this is a filter that removes all details.
      do i = 1, size(u,4)
          call wavelet_filter(u(:,:,:,i), Bs, g )
      enddo

  end select


end subroutine




!=====================================================================
!  WAVELET FILTER
!=====================================================================
subroutine wavelet_filter(block_data, Bs, g )
    use module_interpolation

    implicit none
    !> heavy data array - block data
    real(kind=rk), intent(inout) :: block_data(:, :, :)
    !> mesh params
    integer(kind=ik), dimension(3), intent(in) :: Bs
    integer(kind=ik), intent(in) :: g
    real(kind=rk), allocatable, save :: u3(:,:,:)
    real(kind=rk), allocatable, save :: u2(:,:,:)


    if (params_acm%dim == 3) then
        ! ********** 3D **********
        if (.not.allocated(u2)) allocate(u2( 1:Bs(1)+2*g, 1:Bs(2)+2*g, 1:Bs(3)+2*g ))
        if (.not.allocated(u3)) allocate(u3((Bs(1)+1)/2+g,(Bs(2)+1)/2+g,(Bs(3)+1)/2+g))

        if (params_acm%harten_multiresolution) then
            ! now, coarsen array u1 (restriction)
            call restriction_3D( block_data, u3 )  ! fine, coarse
            ! then, re-interpolate to the initial level (prediciton)
            call prediction_3D ( u3, block_data, params_acm%order_predictor )  ! coarse, fine
        else
            ! apply a smoothing filter (low-pass, in wavelet terminology h_tilde)
            call restriction_prefilter_3D( block_data(:, :, :), u2(:,:,:), params_acm%wavelet )
            ! now, coarsen array u1 (restriction)
            call restriction_3D( u2, u3 )  ! fine, coarse
            ! then, re-interpolate to the initial level (prediciton)
            call prediction_3D ( u3, block_data, params_acm%order_predictor )  ! coarse, fine
        endif

    else
        ! ********** 2D **********
        if (.not.allocated(u2)) allocate(u2( 1:Bs(1)+2*g, 1:Bs(2)+2*g, 1 ))
        if (.not.allocated(u3)) allocate(u3((Bs(1)+1)/2+g,(Bs(2)+1)/2+g,1))

        if (params_acm%harten_multiresolution) then
            ! now, coarsen array u1 (restriction)
            call restriction_2D( block_data(:,:,1), u3(:,:,1) )  ! fine, coarse
            ! then, re-interpolate to the initial level (prediciton)
            call prediction_2D ( u3(:,:,1), block_data(:,:,1), params_acm%order_predictor )  ! coarse, fine
        else
            ! apply a smoothing filter (low-pass, in wavelet terminology h_tilde)
            call restriction_prefilter_2D(block_data(:, :, 1), u2(:,:,1), params_acm%wavelet)
            ! now, coarsen array u1 (restriction)
            call restriction_2D( u2(:,:,1), u3(:,:,1) )  ! fine, coarse
            ! then, re-interpolate to the initial level (prediciton)
            call prediction_2D ( u3(:,:,1), block_data(:,:,1), params_acm%order_predictor )  ! coarse, fine
        endif

    end if

end subroutine wavelet_filter




!=====================================================================
!  1D FILTER
!=====================================================================
!> \brief 1D Filter subroutine
!> \details
!> \version 0.5
!> \author msr
!! \date 27/03/17 - create
!! \date 02/05/17 - return filtered value separatly
subroutine filter_1D(phi, phi_tilde, a)

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> datafield
    real(kind=rk), intent(in)           :: phi(:)
    !> filtered value
    real(kind=rk), intent(out)          :: phi_tilde
    !> filter coefficients
    real(kind=rk), intent(in)           :: a(:)

    ! loop variable
    integer(kind=ik)                    :: k

    ! old values
    real(kind=rk)                       :: phi_old(size(phi,1))

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    phi_old   = phi
    phi_tilde = 0.0_rk

!---------------------------------------------------------------------------------------------
! main body

    ! check filter stencil
    if ( size(phi) /= size(a) ) then
        write(*,'(80("_"))')
        print*, phi
        print*, a
        call abort(123980,"ERROR: filter stencil has wrong size")
    end if

    ! filter data
    do k = 1, size(a)
        phi_tilde = phi_tilde + a(k)*phi_old(k)
    end do

end subroutine filter_1D
