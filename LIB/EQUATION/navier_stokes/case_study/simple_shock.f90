!> \brief reads parameters for mask function from file
subroutine init_simple_shock(params_ns,FILE)
  use module_navier_stokes_params
  implicit none
  !> params structure of navier stokes
   type(type_params_ns),intent(in)  :: params_ns
 ! character(len=*), intent(in) :: filename
  type(inifile) , intent(inout) :: FILE


  call read_param_mpi(FILE, 'Navier_Stokes', 'initial_pressure', shock_params%p_left,1.0_rk )
  call read_param_mpi(FILE, 'Navier_Stokes', 'initial_velocity', shock_params%u_left, 3.0_rk )
  call read_param_mpi(FILE, 'Navier_Stokes', 'initial_density', shock_params%rho_left, 2.0_rk )

  call shockVals( shock_params%rho_left,shock_params%u_left,shock_params%p_left, &
                  shock_params%rho_right,shock_params%u_right,shock_params%p_right,params_ns%gamma_)


end subroutine init_simple_shock



!==========================================================================
!> \brief Compute mask function of sod shock tube
!> \detail
!>      +For boundary_type=='left'
!>             -> mask will be generatet for 0<x<0.1*Lx
!>      +For boundary_type=='rigth'
!>             -> mask will be generatet for 0.9Lx<x<Lx
!>      +For boundary_type else
!>             -> mask will be generatet for both sides
subroutine draw_simple_shock(mask, x0, dx, Bs, g, boundary_type )


    implicit none

    ! grid
    integer(kind=ik), intent(in)                              :: Bs, g
    !> mask term for every grid point of this block
    real(kind=rk), dimension(:,:), intent(out)     :: mask
    !> spacing and origin of block
    real(kind=rk), dimension(2), intent(in)                   :: x0, dx
    !> boundary_type={'left','rigth'}
     character(len=*),          intent(in)                    :: boundary_type
    ! auxiliary variables
    real(kind=rk)                                             :: x, h, x_boundary(2),boundary_width
    ! loop variables
    integer(kind=ik)                                          :: ix

!---------------------------------------------------------------------------------------------
! variables initialization
    if (size(mask,1) /= Bs+2*g) call abort(777109,"wrong array size, there's pirates, captain!")

    ! reset mask array
    mask = 0.0_rk
    ! left and right boundary of shock tube
    x_boundary(1)   =0.1_rk*domain_size(1)
    x_boundary(2)   =domain_size(1)*0.9_rk
    ! parameter for smoothing function (width)
    h = 1.5_rk*max(dx(1), dx(2))

       do ix=1, Bs+2*g
           x = dble(ix-(g+1)) * dx(1) + x0(1)

              mask(ix,:) = smoothstep(x-x_boundary(1),h) &
                         + smoothstep(x_boundary(2)-x,h)

       end do

end subroutine draw_simple_shock
!==========================================================================



subroutine add_simple_shock(penalization, x0, dx, Bs, g ,phi)


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
    real(kind=rk)                                    :: mask

    ! auxiliary variables
    real(kind=rk)                                    :: x, y, h
    ! preasure,density velocities
    real(kind=rk)                                    :: p(Bs+2*g,Bs+2*g),rho(Bs+2*g,Bs+2*g), &
                                                        u(Bs+2*g,Bs+2*g),v(Bs+2*g,Bs+2*g)
    ! loop variables
    integer(kind=ik)                                 :: ix, iy,n
    ! left and right boundary
    real(kind=rk)                                    :: u_ref,u_R,u_L,rho_R,rho_L,p_L,p_R,x_L,x_R,rho_ref,p_ref

!---------------------------------------------------------------------------------------------
! variables initialization
    if (size(penalization,1) /= Bs+2*g) call abort(777109,"wrong array size, there's pirates, captain!")

    ! reset mask array
    mask         = 0.0_rk
    penalization = 0.0_rk

    ! parameter for smoothing function (width)
    h       = 1.5_rk*max(dx(1), dx(2))

    rho         = phi(:,:,1)**2
    u           = phi(:,:,2)/phi(:,:,1)
    v           = phi(:,:,3)/phi(:,:,1)
    p           = phi(:,:,4)

    ! penalization term
    !-------------------

    ! left boundary

    x_L      =0.1_rk*domain_size(1)
    p_L      = shock_params%p_left
    rho_L    =  shock_params%rho_left
    u_L      =  shock_params%u_left
    ! right boundary

    x_R      =domain_size(1)*0.9_rk
    p_R      = shock_params%p_right
    rho_R    =  shock_params%rho_right
    u_R      =  shock_params%u_right
    ! parameter for smoothing function (width)

    do ix=1, Bs+2*g
      x = dble(ix-(g+1)) * dx(1) + x0(1)
      if (x<domain_size(1)*0.5_rk) then
        mask      = smoothstep(x-x_L,h)
        rho_ref   = rho_L
        p_ref     = p_L
        u_ref   = u_L
      else
        mask      = smoothstep(x_R-x,h)
        rho_ref   = transition(x,0.925_rk*domain_size(1),0.05_rk*domain_size(1),rho_R,rho_L)
        p_ref     = transition(x,0.925_rk*domain_size(1),0.05_rk*domain_size(1),p_R  ,p_L    )
        u_ref     = transition(x,0.925_rk*domain_size(1),0.05_rk*domain_size(1),u_R  ,u_L    )
      endif

      ! density
      penalization(ix,:,1)= penalization(ix,:,1) + C_sp_inv*mask * ( rho(ix,:) - rho_ref )
      ! x-velocity
      penalization(ix,:,2)= penalization(ix,:,2) + C_sp_inv*mask * ( rho(ix,:)*u(ix,:)- rho_ref*u_ref )
      ! y-velocity
      penalization(ix,:,3)= penalization(ix,:,3) + C_sp_inv*mask * ( rho(ix,:)*v(ix,:) )
      ! preasure
      penalization(ix,:,4)= penalization(ix,:,4) + C_sp_inv*mask *( p(ix,:) - p_ref )
    end do



end subroutine add_simple_shock
