!-------------------------------------------------------------------------------
! get plate geometry.
! the plate has an upper and lower surface, which are both 2D arrays of 3D points
! they are described in the plate coordinate system, in which their x,y coordinates
! are given by the beam, and the z coordinate is the height.
!-------------------------------------------------------------------------------
! as of 07 Aug 2014, this subroutine takes fewer arguments; the arrays surfaces
! heights and active_points are global in the solid solver scope. this is a re-
! quirement for running on turing
!-------------------------------------------------------------------------------
! subroutine plate_geometry( beam, surfaces, active_points, nh, M_plate, x0_plate )
subroutine plate_geometry( time )
    use mpi
    implicit none
    real(kind=rk), intent(in) :: time

    integer :: ih, is, isurf, k, nh
    real(kind=rk) :: dh,xf,yf,s,ztop,zbot

    !-----------------------------------------------------------------------------
    !-- decide how many points in the rigid direction
    !-----------------------------------------------------------------------------
    !-- "true" 3D case
    !-- number of interpolation points in rigid direction (span)
    nh = nint( L_span/dx_min  )
    !-- spacing in span direction (which is actually dx_min, you checked)
    dh = L_span/dble(nh)
    !-- allocate "heights" array only once, as nh will not change over time
    if(.not.allocated(heights)) allocate(heights(0:nh))
    !-- fill heights array
    do ih = 0, nh
        heights(ih) = real(ih, kind=rk)*dh - 0.5_rk*L_span
    enddo

    ! this array holds the interpolation points (2 2D arrays of 3D vectors = 4 indices)
    ! allocate arrays only once, as ns,nh will not change over time
    if(.not.allocated(surfaces)) allocate(surfaces(0:ns-1,0:nh,1:2,1:3))
    if(.not.allocated(p_surface)) allocate(p_surface(0:ns-1,0:nh,1:2))

    surfaces = 999.9_rk

    ! ! by default, all points are active
    ! active_points = 1.d0
    ! ! decide which points are active (1.0) and which are not (0.0). The latter
    ! ! will not be used for computing the surface pressure points
    ! do is = 0, ns-1
    !   s = dble(is)*ds
    !   ztop = +z_top(s)
    !   zbot = -z_bottom(s)

    !   do ih=0,nh
    !     if ((heights(ih)>=zbot).and.(heights(ih)<=ztop)) then
    !       active_points(is,ih) = 1.d0
    !     else
    !       active_points(is,ih) = 0.d0
    !     endif
    !   enddo
    ! enddo


    !-----------------------------------------------------------------------------
    ! Compute the two 2D surfaces in 3D space (top and bottom)
    !-----------------------------------------------------------------------------
    do is = 0, ns-1
        do ih = 0, nh
            !-- top surface points
            xf = beams(1)%x(is) - t_beam*dsin(beams(1)%theta(is))
            yf = beams(1)%y(is) + t_beam*dcos(beams(1)%theta(is))
            
            !-- point in global coordinate system
            surfaces(is,ih,1,1:3) = (/ xf, yf, heights(ih) /) + beams(1)%x0_origin

            !-- bottom surface points
            xf = beams(1)%x(is) + t_beam*dsin(beams(1)%theta(is))
            yf = beams(1)%y(is) - t_beam*dcos(beams(1)%theta(is))
            
            !-- point in global coordinate system
            surfaces(is,ih,2,1:3) = (/ xf, yf, heights(ih) /) + beams(1)%x0_origin
        enddo
    enddo

    ! construction of xq_solid for the interpolation routine
    if (.not. allocated(xq_solid)) then
        allocate( xq_solid(1:3, 1:ns*(nh+1)*2) )
    endif

    if (.not. allocated(pq_solid_old)) then
        allocate( pq_solid_old(1, 1:ns*(nh+1)*2) )
    endif

    if (.not. allocated(pq_solid_new)) then
        allocate( pq_solid_new(1, 1:ns*(nh+1)*2) )
    endif

    ! reshaping in the format used by Interpolation routine
    k = 1
    do is = 0, ns-1
        do ih = 0, nh ! note nh+1 values
            do isurf = 1, 2
                xq_solid(1:3, k) = surfaces(is, ih, isurf, 1:3)
                k = k +1
            enddo
        enddo
    enddo

end subroutine plate_geometry


!-------------------------------------------------------------------------------
! shape functions of the plate
! Note the maximum spanwise length is L_span for both parts,
! so each function has to provide z values smaller 0.5*L_span
! (otherwise, there is a conflict with the bounding boxes)
!-------------------------------------------------------------------------------
real(kind=rk) function z_top(s)
    implicit none
    real(kind=rk),intent(in)::s
    real(kind=rk)::R

    select case(plate_shape)
      case("rectangular")
        z_top = 0.5d0*L_span

      ! case("fish")
      !   z_top = 0.25*L_span + 0.25*L_span*s**2

      ! case("triangular_appendix")
      !   if (s <=0.75d0) then
      !       z_top = (0.5d0*L_span-0.05d0)
      !   else
      !       z_top = (0.5d0*L_span-0.05d0) + 0.1d0*((s-0.75d0)/0.25d0)
      !   endif

      ! case("negative_triangular_appendix")
      !   if (s <=0.75d0) then
      !       z_top = (0.5d0*L_span-0.05d0)
      !   else
      !       z_top = (0.5d0*L_span-0.05d0) - 0.1d0*((s-0.75d0)/0.25d0)
      !   endif

      ! case("fish_expanding")
      !   z_top = 0.5d0*0.35d0 + 0.525d0*s**2

      ! case("fish_contracting")
      !   z_top = (0.5d0*1.05d0 - 0.525d0*s**2)+1.0d-3

      ! case("half_circle")
      !   R=0.5d0*L_span
      !   if (s<0.d0) then
      !       z_top = R
      !   elseif (s>=0.d0 .and. s<=1.d0) then
      !       z_top = dsqrt( max(R**2-s**2,0.d0) )
      !   elseif (s > 1.d0) then
      !       z_top = (1.d0-s)
      !   endif
      case default
        call abort(35, "ERROR. plate_shape not defined.."//trim(adjustl(plate_shape)))
    end select
end function



real(kind=rk) function z_bottom(s)
    !use vars
    implicit none
    real(kind=rk),intent(inout)::s
    real(kind=rk)::R

    select case(plate_shape)
      case("rectangular")
        z_bottom = 0.5d0*L_span

        ! case("fish")
        !   z_bottom = 0.25*L_span + 0.25*L_span*s**2

        ! case("fish_expanding")
        !   z_bottom = 0.5d0*0.35d0 + 0.525d0*s**2

        ! case("fish_contracting")
        !   z_bottom = (0.5d0*1.05d0 - 0.525d0*s**2)+1.0d-3

        ! case("triangular_appendix")
        !   if (s <=0.75d0) then
        !     z_bottom = (0.5d0*L_span-0.05d0)
        !   else
        !     z_bottom = (0.5d0*L_span-0.05d0) + 0.1d0*((s-0.75d0)/0.25d0)
        !   endif

        ! case("negative_triangular_appendix")
        !   if (s <=0.75d0) then
        !     z_bottom = (0.5d0*L_span-0.05d0)
        !   else
        !     z_bottom = (0.5d0*L_span-0.05d0) - 0.1d0*((s-0.75d0)/0.25d0)
        !   endif

        ! case("half_circle")
        !   R=0.5d0*L_span
        !   if (s<0.d0) then
        !     z_bottom = R
        !   elseif (s>=0.d0 .and. s<=1.d0) then
        !     z_bottom = dsqrt( max(R**2-s**2,0.d0) )
        !   elseif (s > 1.d0) then
        !     z_bottom = (1.d0-s)
        !   endif

      case default
        call abort(36, "ERROR. plate_shape not defined.."//trim(adjustl(plate_shape)))
    end select
end function
