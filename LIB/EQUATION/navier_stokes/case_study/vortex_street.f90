!> \brief reads parameters for mask function from file
subroutine init_vortex_street(FILE)

  implicit none

 ! character(len=*), intent(in) :: filename
  type(inifile) , intent(inout) :: FILE


  call read_param_mpi(FILE, 'VPM', 'x_cntr', cyl%x_cntr,(/ 0.25_rk , 0.5_rk/) )
  call read_param_mpi(FILE, 'VPM', 'length', cyl%radius,0.05_rk*min(domain_size(1),domain_size(2)) )
  cyl%radius=cyl%radius*0.5_rk
  cyl%x_cntr(1)=cyl%x_cntr(1)*domain_size(1)
  cyl%x_cntr(2)=cyl%x_cntr(2)*domain_size(2)
end subroutine init_vortex_street



!==========================================================================
subroutine draw_cylinder(mask, x0, dx, Bs, g )


    implicit none

    ! grid
    integer(kind=ik), intent(in)                              :: Bs, g
    !> mask term for every grid point of this block
    real(kind=rk), dimension(:,:), intent(out)     :: mask
    !> spacing and origin of block
    real(kind=rk), dimension(2), intent(in)                   :: x0, dx

    ! auxiliary variables
    real(kind=rk)                                             :: x, y, r, h
    ! loop variables
    integer(kind=ik)                                          :: ix, iy

!---------------------------------------------------------------------------------------------
! variables initialization
    if (size(mask,1) /= Bs+2*g) call abort(777109,"wrong array size, there's pirates, captain!")

    ! reset mask array
    mask = 0.0_rk

!---------------------------------------------------------------------------------------------
! main body


    ! parameter for smoothing function (width)
    h = 1.5_rk*max(dx(1), dx(2))

    do iy=1, Bs+2*g
       y = dble(iy-(g+1)) * dx(2) + x0(2) - cyl%x_cntr(2)
       do ix=1, Bs+2*g
           x = dble(ix-(g+1)) * dx(1) + x0(1) - cyl%x_cntr(1)
           ! distance from center of cylinder
           r = dsqrt(x*x + y*y)
           mask(ix,iy) = smoothstep( r - cyl%radius, h)

       end do
    end do
end subroutine draw_cylinder
!==========================================================================



subroutine add_cylinder(penalization, x0, dx, Bs, g ,phi)


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
    real(kind=rk)                                    :: T0,rho0,u0,p0

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

    ! parameter for smoothing function (width)
    h       = 1.5_rk*max(dx(1), dx(2))

    rho         = phi(:,:,1)**2
    u           = phi(:,:,2)/phi(:,:,1)
    v           = phi(:,:,3)/phi(:,:,1)
    p           = phi(:,:,4)

    ! cylinder
    !---------
    call draw_cylinder(mask, x0, dx, Bs, g )

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


end subroutine add_cylinder
