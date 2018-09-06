 !> \brief reads parameters for mask function from file
subroutine init_triangle(params,FILE)

  implicit none
  !> params structure of navier stokes
   type(type_params_ns),intent(inout)  :: params
 ! character(len=*), intent(in) :: filename
  type(inifile) , intent(inout) :: FILE

  ! read in the angle between the symmetry axis and the triangle side
  call read_param_mpi(FILE, 'VPM', 'angle', triangle%angle, 30.0_rk )
  call read_param_mpi(FILE, 'VPM', 'x_cntr', triangle%x_cntr, (/0.2_rk, 0.5_rk, 0.0_rk/) )
  call read_param_mpi(FILE, 'VPM', 'length', triangle%length, 0.1*params%domain_size(1) )

  ! convert degrees to radians
  if ( triangle%angle>0.0_rk .and. triangle%angle<90.0_rk ) then
    triangle%angle=triangle%angle*PI/180.0_rk
  else
    call abort(45756,"somebody has to go back to preeshool! 0< angle <90")
  end if

  !convert from height to length
  triangle%length=triangle%length/(2*tan(triangle%angle))

  ! a rhombus is only a double triangle
  if ( mask_geometry=="rhombus" ) then
    triangle%rhombus    =.true.
  else
    triangle%rhombus    =.false.
  end if

  if (triangle%x_cntr(1)>1.0_rk .or. triangle%x_cntr(1)<0.0_rk ) then
    triangle%x_cntr(1) =0.2_rk
  end if

  if (triangle%x_cntr(2)>1.0_rk .or. triangle%x_cntr(2)<0.0_rk ) then
      triangle%x_cntr(2) =0.5_rk
  end if
  triangle%x_cntr(1)=triangle%x_cntr(1)*params%domain_size(1)
  triangle%x_cntr(2)=triangle%x_cntr(2)*params%domain_size(2)
  triangle%x_cntr(3)=triangle%x_cntr(3)*params%domain_size(3)

end subroutine init_triangle



!==========================================================================
subroutine draw_triangle(mask, x0, dx, Bs, g )

    implicit none
    ! grid
    integer(kind=ik), intent(in)                              :: Bs, g
    !> mask term for every grid point of this block
    real(kind=rk), dimension(:,:), intent(out)              :: mask
    !> spacing and origin of block
    real(kind=rk), dimension(2), intent(in)                   :: x0, dx

    ! auxiliary variables
    real(kind=rk)                                             :: x, y,h,length,tan_theta,height
    ! loop variables
    integer(kind=ik)                                          :: ix, iy
    !
    logical                                                   :: triangle_is_rhombus

!---------------------------------------------------------------------------------------------
! variables initialization
    if (size(mask,1) /= Bs+2*g) call abort(777109,"wrong array size, there's pirates, captain!")

    ! reset mask array
    mask  = 0.0_rk

    triangle_is_rhombus =triangle%rhombus
    length              =triangle%length
    tan_theta           =tan(triangle%angle)

!---------------------------------------------------------------------------------------------
! main body


    ! parameter for smoothing function (width)
    h = 1.5_rk*max(dx(1), dx(2))

    do ix=1, Bs+2*g
        x = dble(ix-(g+1)) * dx(1) + x0(1)
        x = x-triangle%x_cntr(1)
        !YES a rhombus is made from two isosceles triangles
        if ( triangle_is_rhombus .and. x>length ) then
          ! reflect triangle at the line at x=length parallel to y axis: <|>
          x = 2 * length - x
        end if
        height = tan_theta*x
        do iy=1, Bs+2*g
           y = dble(iy-(g+1)) * dx(2) + x0(2)
           y = abs(y-triangle%x_cntr(2))

           mask(ix,iy) = soft_bump(x,0.0_rk,length+h,h)*smoothstep(y-height,h)

       end do
    end do
end subroutine draw_triangle
!==========================================================================





subroutine add_triangle(penalization, x0, dx, Bs, g ,phi)


    implicit none

    ! grid
    integer(kind=ik), intent(in)                     :: Bs, g
    !> penalization term including mask
    real(kind=rk), dimension(:,:,:), intent(inout)   :: penalization
    !> spacing and origin of block
    real(kind=rk), dimension(2), intent(in)          :: x0, dx
    !> statevector
    real(kind=rk), dimension(:,:,:), intent(in)      :: phi
      !> statevector
    real(kind=rk)                                    :: mask(Bs+2*g,Bs+2*g)

    ! auxiliary variables
    real(kind=rk)                                    :: x, y, r, h,velocity
    ! preasure,density velocities
    real(kind=rk)                                    :: p(Bs+2*g,Bs+2*g),rho(Bs+2*g,Bs+2*g), &
                                                        u(Bs+2*g,Bs+2*g),v(Bs+2*g,Bs+2*g)
    ! loop variables
    integer(kind=ik)                                 :: ix, iy,n
    ! outlets and inlets
    real(kind=rk)                                    :: T0,rho0,u0,p0,Rs

!---------------------------------------------------------------------------------------------
! variables initialization
    if (size(penalization,1) /= Bs+2*g) call abort(777109,"wrong array size, there's pirates, captain!")

    ! reset mask array
    mask    = 0.0_rk
    penalization = 0.0_rk


    T0      = params_ns%initial_temp
    rho0    = params_ns%initial_density
    p0      =  params_ns%initial_pressure
    u0      =  params_ns%initial_velocity(1)
    Rs      =  params_ns%Rs


    ! parameter for smoothing function (width)
    h       = 1.5_rk*max(dx(1), dx(2))

    rho         = phi(:,:,1)**2
    u           = phi(:,:,2)/phi(:,:,1)
    v           = phi(:,:,3)/phi(:,:,1)
    p           = phi(:,:,4)

    ! triangle
    !---------
    call draw_triangle(mask, x0, dx, Bs, g )

    ! x-velocity
    penalization(:,:,2)=C_eta_inv*mask * ( rho*u )
    ! y-velocity
    penalization(:,:,3)=C_eta_inv*mask * ( rho*v )
    ! preasure
    penalization(:,:,4)=C_eta_inv*mask *( p- rho*Rs*T0 )


    ! sponge
    !--------
    call simple_sponge(mask, x0, dx, Bs, g)


    penalization(:,:,1)= penalization(:,:,1) + C_sp_inv*mask * ( rho - rho0 )
    ! x-velocity
    penalization(:,:,2)= penalization(:,:,2) + C_sp_inv*mask * ( rho*u - rho0*u0  )
    ! y-velocity
    penalization(:,:,3)= penalization(:,:,3) + C_sp_inv*mask * ( rho*v )
    ! preasure
    penalization(:,:,4)= penalization(:,:,4) + C_sp_inv*mask *( p - p0 )



          ! sponge-free outlet
          !--------
          call draw_free_outlet_wall(mask, x0, dx, Bs, g )

          penalization(:,:,1)= penalization(:,:,1) + C_sp_inv*mask * ( rho - rho0 )
          ! x-velocity
          !penalization(:,:,2)= penalization(:,:,2) + C_sp_inv*mask * ( rho*u - rho0*u0  )
          ! y-velocity
          penalization(:,:,3)= penalization(:,:,3) + C_sp_inv*mask * ( rho*v )
          ! preasure
          penalization(:,:,4)= penalization(:,:,4) + C_sp_inv*mask *( p - p0 )


end subroutine add_triangle
