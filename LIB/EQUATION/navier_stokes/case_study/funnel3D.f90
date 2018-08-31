!#########################################################################################
!                     3D FUNNEL IMPLEMENTATION
!#########################################################################################
subroutine  funnel_penalization3D(penalization, x0, dx, Bs, g ,phi)
  use module_helpers
  implicit none
    ! -----------------------------------------------------------------
    integer(kind=ik), intent(in)  :: Bs, g          !< grid parameter
    real(kind=rk), intent(in)     :: x0(3), dx(3)   !< coordinates of block and block spacinf
    real(kind=rk), intent(in)     :: phi(:,:,:,:)     !< state vector
    real(kind=rk), intent(inout)  :: penalization(:,:,:,:) !< reference values of penalized volume
    real(kind=rk), allocatable,save :: mask(:,:,:,:)    !< mask function
    integer(kind=2), allocatable,save:: mask_color(:,:,:)!< identifyers of mask parts (plates etc)
    logical                       :: mesh_was_adapted=.true.
    ! -----------------------------------------------------------------
    if (.not. allocated(mask_color))  allocate(mask_color(1:Bs+2*g, 1:Bs+2*g,  1:Bs+2*g))
    if (.not. allocated(mask))        allocate(mask(1:Bs+2*g, 1:Bs+2*g,  1:Bs+2*g, 5))

    !!> todo implement function check_if_mesh_adapted (true/false) in adapt mesh
    if ( mesh_was_adapted .eqv. .true. ) then
      ! reset parameters
      mask        = 0.0_rk
      mask_color  = 0
      call draw_funnel3D(x0, dx, Bs, g, mask, mask_color)
      call draw_sponge3D(x0, dx, Bs, g, mask, mask_color)
      call compute_penal3D(mask_color,mask,phi, x0, dx, Bs, g ,penalization)
    else
      call compute_penal3D(mask_color,mask,phi, x0, dx, Bs, g ,penalization)
    end if

end subroutine  funnel_penalization3D






subroutine compute_penal3D(mask_color,mask,phi, x0, dx, Bs, g ,penalization)
    implicit none
    ! -----------------------------------------------------------------
    integer(kind=ik), intent(in)  :: Bs, g            !< grid parameter
    real(kind=rk), intent(in)     :: x0(3), dx(3)     !< coordinates of block and block spacinf
    integer(kind=2), intent(inout):: mask_color(:,:,:)!< identifyers of mask parts (plates etc)
    real(kind=rk), intent(in)     ::mask(:,:,:,:)     !< mask function
    real(kind=rk), intent(in)     ::phi(:,:,:,:)     !< state vector
    real(kind=rk), intent(inout)  ::penalization(:,:,:,:) !<
    ! -----------------------------------------------------------------
    ! auxiliary variables
    real(kind=rk)     :: x, y, r, z, h,velocity,C_inv
    real(kind=rk)     :: p,rho,u,v,w,chi,v_ref,dq,phi_ref(5)
    ! loop variables
    integer(kind=ik)  :: ix, iy, iz, n
    ! outlets and inlets
    real(kind=rk)     :: velocity_pump,rho_pump,pressure_pump, &
                        rho_capillary,u_capillary,v_capillary,w_capillary,p_capillary, &
                        p_2nd_pump_stage,rho_2nd_pump_stage
    ! smooth width of jet
    real(kind=rk)     ::jet_smooth_width,pump_smooth_width


    u_capillary       =funnel%inlet_velocity(1)
    v_capillary       =funnel%inlet_velocity(2)
    w_capillary       =0.0_rk
    rho_capillary     =funnel%inlet_density
    rho_pump          =funnel%pump_density
    p_capillary       =funnel%inlet_pressure
    velocity_pump     =funnel%pump_speed
    pressure_pump     =funnel%pump_pressure
    p_2nd_pump_stage  =funnel%outlet_pressure
    rho_2nd_pump_stage=funnel%outlet_density
    ! parameter for smoothing function (width)
    h  = 1.5_rk*maxval(dx(1:params_ns%dim))

    if ( 3*minval(dx(2:3)) <= 0.05_rk*funnel%jet_radius ) then
      jet_smooth_width = 0.05_rk*funnel%jet_radius
    else
      jet_smooth_width = 3*minval(dx(2:3))
      !call abort('ERROR [funnel.f90]: discretication constant dy to large')
    endif

    if ( 3*dx(1)<=0.1_rk*funnel%pump_diameter ) then
      pump_smooth_width = 0.025_rk*funnel%pump_diameter
    else
      pump_smooth_width = 3*h
      !call abort('ERROR [funnel.f90]: discretication constant dy to large')
    endif



  do iz=1, Bs+2*g
    z = dble(iz - (g+1)) * dx(3) + x0(3)
    do iy=1, Bs+2*g
       y = dble(iy-(g+1)) * dx(2) + x0(2)
       r = sqrt((y-R_domain)**2+(z-R_domain)**2)
       do ix=1, Bs+2*g
            x = dble(ix-(g+1)) * dx(1) + x0(1)
            rho             = phi(ix,iy,iz,rhoF)**2
            u   = phi(ix,iy,iz,UxF)/phi(ix,iy,iz,rhoF)
            v   = phi(ix,iy,iz,UyF)/phi(ix,iy,iz,rhoF)
            w   = phi(ix,iy,iz,UzF)/phi(ix,iy,iz,rhoF)
            p   = phi(ix,iy,iz,pF)

            C_inv=C_sp_inv
            ! solid obstacles: plates and walls
            ! ---------------------------------
            if ( mask_color(ix,iy,iz) == color_plates &
            .or. mask_color(ix,iy,iz) == color_walls ) then
              Phi_ref(UxF) = 0.0_rk                     ! no velocity in x
              Phi_ref(UyF) = 0.0_rk                     ! no velocity in y
              Phi_ref(UzF) = 0.0_rk                     ! no velocity in z
              Phi_ref(pF) = rho*Rs*funnel%temperatur   ! pressure set according to
              C_inv=C_eta_inv
            endif                                             ! the temperature of the funnel

            ! Outlet flow: PUMPS
            ! ------------------
            if (mask_color(ix,iy,iz) == color_pumps) then
              !compute velocity profile
              Phi_ref(UxF) = 0
              v_ref=velocity_pump*jet_stream(abs(x-funnel%pump_x_center),funnel%pump_diameter*0.5_rk,pump_smooth_width)
              if (y>R_domain) then
                phi_ref(UyF) = rho*v_ref
              else
                phi_ref(UyF) = -rho*v_ref
              endif
            endif
            !energy sink
            if ( mask_color(ix,iy,iz) == color_pumps_sink) then
              Phi_ref(rhoF) = rho_pump
              Phi_ref(pF) = pressure_pump
              C_inv=C_eta_inv
            endif

            if (mask_color(ix,iy,iz) == color_capillary) then
              dq               =jet_stream(r,funnel%jet_radius,jet_smooth_width)
              C_inv=C_sp_inv
              phi_ref(rhoF) =  rho_capillary
              phi_ref(UxF) =  rho_capillary*u_capillary*dq
              phi_ref(UyF) =  rho_capillary*v_capillary
              phi_ref(UzF) =  rho_capillary*w_capillary
              phi_ref(pF) =  p_capillary  !rho*Rs*funnel%temperatur * (1 - dq) + p_capillary * dq
            endif

             ! Outlet flow: Transition to 2pump
             ! ---------------------
              if ( mask_color(ix,iy,iz) == color_outlet) then
                Phi_ref(rhoF) = rho_2nd_pump_stage
                !Phi_ref(2) = 0
                Phi_ref(UyF) = 0
                Phi_ref(UzF) = 0
                Phi_ref(pF) = p_2nd_pump_stage
                C_inv=C_sp_inv
              endif


              ! density
              penalization(ix,iy,iz,rhoF)=C_inv*mask(ix,iy,iz,rhoF)*(rho-  Phi_ref(rhoF) )
              ! x-velocity
              penalization(ix,iy,iz,UxF)=C_inv*mask(ix,iy,iz,UxF)*(rho*u-  Phi_ref(UxF) )
              ! y-velocity
              penalization(ix,iy,iz,UyF)=C_inv*mask(ix,iy,iz,UyF)*(rho*v-  Phi_ref(UyF) )
              ! z-velocity
              penalization(ix,iy,iz,UzF)=C_inv*mask(ix,iy,iz,UzF)*(rho*v-  Phi_ref(UzF) )
              ! preasure
              penalization(ix,iy,iz,pF)=C_inv*mask(ix,iy,iz,pF)*(p- Phi_ref(pF) )

       end do
    end do
  end do
end subroutine compute_penal3D


subroutine draw_funnel3D(x0, dx, Bs, g, mask, mask_color)
    implicit none
    ! -----------------------------------------------------------------
    integer(kind=ik), intent(in)  :: Bs, g                             !< grid parameter
    real(kind=rk), intent(in)     :: x0(3), dx(3)                      !< coordinates of block and block spacinf
    real(kind=rk), intent(inout)  :: mask(:,:,:,:)                     !< mask function
    integer(kind=2), intent(inout):: mask_color(:,:,:)                 !< identifyers of mask parts (plates etc)
    ! -----------------------------------------------------------------
    ! auxiliary variables
    real(kind=rk)                     :: x, y, z, r, h
    ! preasure,density velocities
    real(kind=rk)                     :: chi
    ! loop variables
    integer(kind=ik)                                 :: ix, iy, iz, n


    ! reset mask array
    mask        = 0.0_rk
    mask_color  = 0

    h  = 1.5_rk*maxval(dx(1:params_ns%dim))

  do iz=g+1, Bs+g
    z = dble(iz - (g+1)) * dx(3) + x0(3)
    do iy=g+1, Bs+g
       y = dble(iy-(g+1)) * dx(2) + x0(2)
       r = sqrt((y-R_domain)**2+(z-R_domain)**2)
       do ix=g+1, Bs+g
            x = dble(ix-(g+1)) * dx(1) + x0(1)

            ! plates
            ! ------
            chi = draw_funnel_plates(x,r,funnel,h)
            if (chi>0.0_rk) then
              mask_color(ix,iy,iz)  = color_plates
              mask(ix,iy,iz,2:5)    = mask(ix,iy,iz,2:5) + chi
            endif                                           ! the temperature of the funnel

            ! ! Walls
            ! ! -----
            ! chi = draw_walls(x,r,funnel,h)
            ! if (chi>0.0_rk) then                       ! default values on the walls
            !   mask_color(ix,iy,iz)  = color_walls
            !   mask(ix,iy,iz,2:5)    = mask(ix,iy,iz,2:5) + chi
            ! endif
       end do !ix
    end do !iy
  end do !iz
end subroutine draw_funnel3D


subroutine draw_sponge3D(x0, dx, Bs, g, mask, mask_color)
    implicit none
    ! -----------------------------------------------------------------
    integer(kind=ik), intent(in)  :: Bs, g                             !< grid parameter
    real(kind=rk), intent(in)     :: x0(3), dx(3)                      !< coordinates of block and block spacinf
    real(kind=rk), intent(inout)  :: mask(:,:,:,:)                     !< mask function
    integer(kind=2), intent(inout):: mask_color(:,:,:)                 !< identifyers of mask parts (plates etc)
    ! -----------------------------------------------------------------
    ! auxiliary variables
    real(kind=rk)                     :: x, y, z, r, h
    ! preasure,density velocities
    real(kind=rk)                     :: chi
    ! loop variables
    integer(kind=ik)                                 :: ix, iy, iz, n

    h  = 1.5_rk*maxval(dx(1:params_ns%dim))

  do iz=g+1, Bs+g
    z = dble(iz - (g+1)) * dx(3) + x0(3)
    do iy=g+1, Bs+g
       y = dble(iy-(g+1)) * dx(2) + x0(2)
       r = sqrt((y-R_domain)**2+(z-R_domain)**2)
       do ix=g+1, Bs+g
            x = dble(ix-(g+1)) * dx(1) + x0(1)
            ! Outlet flow: PUMPS
            ! ------------------
            ! pump volume flow
            ! compute mask
            chi=  draw_pumps_volume_flow(x,r,funnel,h)
            if (chi>0.0_rk) then
              mask_color(ix,iy,iz) = color_pumps
              mask(ix,iy,iz,UxF)    = mask(ix,iy,iz,UxF)+chi
              mask(ix,iy,iz,UyF)    = mask(ix,iy,iz,UyF)+chi
            endif
            ! mass and energy sink
            chi=  draw_pumps_sink(x,r,funnel,h)
            if (chi>0.0_rk) then
              mask_color(ix,iy,iz) = color_pumps_sink
              mask(ix,iy,iz,rhoF) = mask(ix,iy,iz,rhoF)+chi
              mask(ix,iy,iz,pF) = mask(ix,iy,iz,pF)+chi
            endif



            ! Inlet flow: Capillary
            ! ---------------------
            chi=  draw_jet(x,r,funnel,h)
            if (chi>0.0_rk) then
              mask_color(ix,iy,iz) = color_capillary
              mask(ix,iy,iz,:)  =  mask(ix,iy,iz,:)+chi
            endif


            ! ! Outlet flow: Transition to 2pump
            ! ! ---------------------
            chi=  draw_outlet(x,r,funnel,h)
              if (chi>0.0_rk) then
                mask_color(ix,iy,iz) = color_outlet
                mask(ix,iy,iz,rhoF) = mask(ix,iy,iz,rhoF)+chi
                mask(ix,iy,iz,pF) = mask(ix,iy,iz,pF)+chi
              endif
       end do !ix
    end do !iy
  end do !iz
end subroutine draw_sponge3D
