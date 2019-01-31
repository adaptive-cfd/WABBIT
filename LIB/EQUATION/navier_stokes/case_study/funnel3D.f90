!#########################################################################################
!                     3D FUNNEL IMPLEMENTATION
!#########################################################################################
subroutine  funnel_penalization3D(Bs, g, x0, dx, phi, mask, phi_ref)
  use module_helpers
  implicit none
    ! -----------------------------------------------------------------
    integer(kind=ik), intent(in)  :: g          !< grid parameter
    integer(kind=ik), dimension(3), intent(in) :: Bs
    real(kind=rk), intent(in)     :: x0(3), dx(3)   !< coordinates of block and block spacinf
    real(kind=rk), intent(in)     :: phi(:,:,:,:)     !< state vector
    real(kind=rk), intent(inout)  :: phi_ref(:,:,:,:) !< reference values of penalized volume
    real(kind=rk), intent(inout)  :: mask(:,:,:,:)    !< mask function
    integer(kind=2), allocatable,save:: mask_color(:,:,:)!< identifyers of mask parts (plates etc)
    logical                       :: mesh_was_adapted=.true.
    ! -----------------------------------------------------------------
    if (.not. allocated(mask_color))  allocate(mask_color(1:Bs(1)+2*g, 1:Bs(2)+2*g,  1:Bs(3)+2*g))
    !!> todo implement function check_if_mesh_adapted (true/false) in adapt mesh
    if ( mesh_was_adapted .eqv. .true. ) then
      ! dont switch the order of draw_funnel3D and draw_sponge3D,
      ! because mask and color are reset in the draw_funnel
      call draw_funnel3D(x0, dx, Bs, g, mask, mask_color)
      call draw_sponge3D(x0, dx, Bs, g, mask, mask_color)
    end if

    call compute_penal3D(mask_color,mask,phi, x0, dx, Bs, g ,phi_ref)

end subroutine  funnel_penalization3D




subroutine draw_funnel3D(x0, dx, Bs, g, mask, mask_color)
    implicit none
    ! -----------------------------------------------------------------
    integer(kind=ik), intent(in)  :: g          !< grid parameter
    integer(kind=ik), dimension(3), intent(in) :: Bs
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

    if (size(mask,1) /= Bs(1)+2*g) call abort(77559,"wrong array size of mask function")

    h  = 1.5_rk*maxval(dx(1:params_ns%dim))

    do iz=g+1, Bs(3)+g
      z = dble(iz - (g+1)) * dx(3) + x0(3)
      do iy=g+1, Bs(2)+g
         y = dble(iy-(g+1)) * dx(2) + x0(2)
         r = sqrt((y-R_domain)**2+(z-R_domain)**2)
         do ix=g+1, Bs(1)+g
              x = dble(ix-(g+1)) * dx(1) + x0(1)

              ! reset the mask function here
              mask_color(ix,iy,iz)=0
              mask(ix,iy,iz,:)=0.0_rk

              ! plates
              ! ------

              chi = draw_funnel_plates(x,r,funnel,h)
              if (chi>0.0_rk) then
                mask_color(ix,iy,iz)  = color_plates
                mask(ix,iy,iz,2:5)    = mask(ix,iy,iz,2:5) + chi
              endif                                           ! the temperature of the funnel

               ! Walls
               ! -----
               chi = draw_walls3D(x,y,z,r,funnel,h)
               if (chi>0.0_rk) then                       ! default values on the walls
                 mask_color(ix,iy,iz)  = color_walls
                 mask(ix,iy,iz,2:5)    = mask(ix,iy,iz,2:5) + chi
               endif
         end do !ix
      end do !iy
    end do !iz
end subroutine draw_funnel3D


subroutine draw_sponge3D(x0, dx, Bs, g, mask, mask_color)
    implicit none
    ! -----------------------------------------------------------------
    integer(kind=ik), intent(in)  :: g          !< grid parameter
    integer(kind=ik), dimension(3), intent(in) :: Bs
    real(kind=rk), intent(in)     :: x0(3), dx(3)      !< coordinates of block and block spacinf
    real(kind=rk), intent(inout)  :: mask(:,:,:,:)     !< mask function
    integer(kind=2), intent(inout):: mask_color(:,:,:) !< identifyers of mask parts (plates etc)
    ! -----------------------------------------------------------------
    ! auxiliary variables
    real(kind=rk)                     :: x, y, z, r, h
    ! preasure,density velocities
    real(kind=rk)                     :: chi
    ! loop variables
    integer(kind=ik)                                 :: ix, iy, iz, n

    h  = 1.5_rk*maxval(dx(1:params_ns%dim))

  do iz=g+1, Bs(3)+g
    z = dble(iz - (g+1)) * dx(3) + x0(3)
    do iy=g+1, Bs(2)+g
       y = dble(iy-(g+1)) * dx(2) + x0(2)
       r = sqrt((y-R_domain)**2+(z-R_domain)**2)
       do ix=g+1, Bs(1)+g
            x = dble(ix-(g+1)) * dx(1) + x0(1)
            ! Outlet flow: PUMPS
            ! ------------------
            ! pump volume flow
            ! compute mask
            chi=  draw_pumps_volume_flow3D(x,y,z,r,funnel,h)
            if (chi>0.0_rk) then
              mask_color(ix,iy,iz) = color_pumps
              mask(ix,iy,iz,UxF)    = mask(ix,iy,iz,UxF)+chi
              mask(ix,iy,iz,UyF)    = mask(ix,iy,iz,UyF)+chi
              mask(ix,iy,iz,UzF)    = mask(ix,iy,iz,UzF)+chi
            endif
            ! mass and energy sink
            chi=  draw_pumps_sink3D(x,y,z,r,funnel,h)
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



subroutine compute_penal3D(mask_color, mask, phi, x0, dx, Bs, g ,phi_ref)
    implicit none
    ! -----------------------------------------------------------------
    integer(kind=ik), intent(in)  :: g          !< grid parameter
    integer(kind=ik), dimension(3), intent(in) :: Bs
    real(kind=rk), intent(in)     :: x0(3), dx(3)   !< coordinates of block and block spacinf
    integer(kind=2), intent(inout):: mask_color(:,:,:)!< identifyers of mask parts (plates etc)
    real(kind=rk), intent(in)     :: phi(:,:,:,:)     !< state vector
    real(kind=rk), intent(inout)  :: mask(:,:,:,:)     !< mask
    real(kind=rk), intent(inout)  :: phi_ref(:,:,:,:)  !< funnel penalization term
    ! -----------------------------------------------------------------
    ! auxiliary variables
    real(kind=rk)     :: x, y, r, z, h,velocity,C_inv
    real(kind=rk)     :: p,rho,u,v,w,chi,v_ref,dq
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



  do iz=g+1, Bs(3)+g
    z = dble(iz - (g+1)) * dx(3) + x0(3)
    do iy=g+1, Bs(2)+g
       y = dble(iy-(g+1)) * dx(2) + x0(2)
       r = sqrt((y-R_domain)**2+(z-R_domain)**2)
       do ix=g+1, Bs(1)+g
            x = dble(ix-(g+1)) * dx(1) + x0(1)

            !reset ref values
            phi_ref(ix,iy,iz,:)=0.0_rk

            ! get primary variables
            rho = phi(ix,iy,iz,rhoF)
            u   = phi(ix,iy,iz,UxF)
            v   = phi(ix,iy,iz,UyF)
            w   = phi(ix,iy,iz,UzF)
            p   = phi(ix,iy,iz,pF)

            C_inv=C_sp_inv
            ! solid obstacles: plates and walls
            ! ---------------------------------
            if ( mask_color(ix,iy,iz) == color_plates &
            .or. mask_color(ix,iy,iz) == color_walls ) then
              phi_ref(ix,iy,iz,UxF) = 0.0_rk                     ! no velocity in x
              phi_ref(ix,iy,iz,UyF) = 0.0_rk                     ! no velocity in y
              phi_ref(ix,iy,iz,UzF) = 0.0_rk                     ! no velocity in z
              phi_ref(ix,iy,iz,pF) = rho*Rs*funnel%temperatur   ! pressure set according to
              C_inv=C_eta_inv
            endif                                             ! the temperature of the funnel

            ! Outlet flow: PUMPS
            ! ------------------
            if (mask_color(ix,iy,iz) == color_pumps) then
              !compute velocity profile
              phi_ref(ix,iy,iz,UxF) = 0
              phi_ref(ix,iy,iz,UzF) = 0
              v_ref=velocity_pump*jet_stream(abs(x-funnel%pump_x_center),funnel%pump_diameter*0.5_rk,pump_smooth_width)
              if (y>R_domain) then
                phi_ref(ix,iy,iz,UyF) = rho*v_ref
              else
                phi_ref(ix,iy,iz,UyF) = -rho*v_ref
              endif
            endif
            !energy sink
            if ( mask_color(ix,iy,iz) == color_pumps_sink) then
              phi_ref(ix,iy,iz,rhoF) = rho_pump
              phi_ref(ix,iy,iz,pF)   = pressure_pump
              C_inv         =C_eta_inv
            endif

            ! Inlet flow: capillary
            ! ------------------
            if (mask_color(ix,iy,iz) == color_capillary) then
              dq               =jet_stream(r,funnel%jet_radius,jet_smooth_width)
              C_inv=C_sp_inv
              phi_ref(ix,iy,iz,rhoF) =  rho_capillary
              phi_ref(ix,iy,iz,UxF) =  rho_capillary*u_capillary*dq
              phi_ref(ix,iy,iz,UyF) =  rho_capillary*v_capillary
              phi_ref(ix,iy,iz,UzF) =  rho_capillary*w_capillary
              phi_ref(ix,iy,iz,pF) =  p_capillary  !rho*Rs*funnel%temperatur * (1 - dq) + p_capillary * dq
            endif

             ! Outlet flow: Transition to 2pump
             ! ---------------------
              if ( mask_color(ix,iy,iz) == color_outlet) then
                phi_ref(ix,iy,iz,rhoF) = rho_2nd_pump_stage
                !phi_ref(ix,iy,iz,2) = 0
                phi_ref(ix,iy,iz,UyF) = 0
                phi_ref(ix,iy,iz,UzF) = 0
                phi_ref(ix,iy,iz,pF) = p_2nd_pump_stage
                C_inv=C_sp_inv
              endif

              ! add the strength parameter to the mask function
              mask(ix,iy,iz,:)=C_inv*mask(ix,iy,iz,:)
       end do
    end do
  end do

end subroutine compute_penal3D





!==========================================================================
  !> Integrates the flow field close to the pump outlet,
  !> and computes the volume of the intagration region.
  subroutine integrate_over_pump_area3D(u,g,Bs,x0,dx,integral,volume)
      implicit none
      ! -----------------------------------------------------------------
      integer(kind=ik), intent(in)  :: g          !< grid parameter
      integer(kind=ik), dimension(3), intent(in) :: Bs
      real(kind=rk), intent(in)            :: u(:,:,:,:)    !< vector of state in pure format
      real(kind=rk), intent(in)            :: x0(3), dx(3)  !< spacing and origin of block
      real(kind=rk),intent(out)            :: integral(params_ns%n_eqn)  !< mean statevector
      real(kind=rk),intent(out)            :: volume        !< volume of the integration domain
      ! -----------------------------------------------------------------
      integer(kind=ik)                   :: ix,iy,iz,neq
      real(kind=rk)                      :: h,r,y,x,z,r0,width
      real(kind=rk),allocatable,save     :: tmp(:)

      neq=params_ns%n_eqn

      if (.not. allocated(tmp) ) allocate(tmp(1:neq+1))
      ! calculate mean density close to the pump
      width =funnel%wall_thickness
      tmp   =  0.0_rk
      r0    =(R_domain-2*funnel%wall_thickness)
      do iz=g+1, Bs(3)+g
       ! relative coordinates z=z-Lz/2
       z = dble(iz-(g+1)) * dx(3) + x0(3) - R_domain
       do iy=g+1, Bs(2)+g
         ! relative coordinates y=y-Ly/2
         y = dble(iy-(g+1)) * dx(2) + x0(2)-R_domain
         ! radius
         r = dsqrt(y**2+z**2)
         do ix=g+1, Bs(1)+g
              !this is the absolut coordinate
              x = dble(ix-(g+1)) * dx(1) + x0(1)
              if (abs(x-funnel%pump_x_center)<= funnel%pump_diameter*0.5_rk .and. &
                  r>r0 .and. r<r0+width) then
                tmp(1:neq)  = tmp(1:neq)+ u(ix,iy,iz,:)
                tmp(neq+1)  = tmp(neq+1)  + 1.0_rk
              endif
          enddo
        enddo
      enddo
        ! integral of all quantities and the volume
        integral  = integral + tmp(1:neq) * product(dx)
        volume    = volume   + tmp(neq+1) * product(dx)

  end subroutine integrate_over_pump_area3D
  !==========================================================================



  !==========================================================================
  function draw_walls3D(x,y,z,r,funnel,h)
    implicit none
    !---------------------------------------------
    real(kind=rk),    intent(in)   :: x, y, z, r, h
    type(type_funnel),intent(in)   ::funnel
    !---------------------------------------------
    real(kind=rk)                  ::  mask, draw_walls3D,r2_rel

    !relative distance (squared) from pump center in the z/x plane
    r2_rel=(x-funnel%pump_x_center)**2+(z-R_domain)**2 ! squaring is cheaper then taking the root

    ! wall in radial direction
    mask = hardstep(R_domain-funnel%wall_thickness-r) &
         * hardstep((funnel%pump_diameter*0.5_rk)**2 - r2_rel)
    mask = mask + hardstep(r2_rel - (funnel%pump_diameter*0.5_rk)**2 ) &
         * hardstep(R_domain-0.333_rk*funnel%wall_thickness-abs(y-R_domain))


    ! wall in east
    !mask=mask+smoothstep(x-funnel%wall_thickness,h)
    if ( r > funnel%jet_radius) then
        mask=mask+hardstep(x-funnel%wall_thickness)
    else
        mask=mask+hardstep(x-funnel%wall_thickness*0.5_rk)
    endif
    ! attach cappilary to wall in EAST
    if (  r > funnel%jet_radius  ) then
        mask=mask+hardstep(x-funnel%plate(1)%x0(1))*hardstep(r-funnel%r_out_cappilary)

           !mask=mask+smoothstep(x-funnel%plate(1)%x0(1),h)*smoothstep(r-funnel%r_out_cappilary,h)
    endif


    ! wall in WEST
    if ( r > funnel%min_inner_diameter*0.5_rk) then
          !  mask=mask+smoothstep(domain_size(1)-x-funnel%wall_thickness,h)

           mask=mask+hardstep(domain_size(1)-x-funnel%wall_thickness)
    else
          mask=mask+hardstep(domain_size(1)-x-funnel%wall_thickness*0.5_rk)
          ! mask=mask+smoothstep(domain_size(1)-x-funnel%wall_thickness*0.5_rk,h)
    endif

     ! is needed because mask off walls overlap
    if (mask>1) then
           mask=1
    endif

    draw_walls3D=mask
  end function draw_walls3D
  !==========================================================================


  !==========================================================================

  function draw_pumps_volume_flow3D(x,y,z,r,funnel,h)
    implicit none
    !-------------------------------------------------
    real(kind=rk),    intent(in)          :: x, y, z, r, h
    type(type_funnel),intent(in)          ::funnel
    !-------------------------------------------------
    real(kind=rk) ::  mask, draw_pumps_volume_flow3D,r0,width,r2_rel,R2_pump

    mask  =0
    width =funnel%wall_thickness*0.333_rk
    r0    =(R_domain-funnel%wall_thickness)
    r2_rel=(x-funnel%pump_x_center)**2+(z-R_domain)**2 ! squaring is cheaper then taking the root
    R2_pump=(funnel%pump_diameter*0.5_rk)**2 !squared pump diameter

    mask = smoothstep(r0-r,h)*smoothstep(r2_rel-(R2_pump-h),h) &
         * smoothstep(abs(y-R_domain)-(R_domain-0.666_rk*funnel%wall_thickness),h)


    draw_pumps_volume_flow3D=mask
  end function draw_pumps_volume_flow3D
  !==========================================================================



  !==========================================================================
  function draw_pumps_sink3D(x,y,z,r,funnel,h)
    implicit none
    !-------------------------------------------------
    real(kind=rk),    intent(in)          :: x, y, z, r, h
    type(type_funnel),intent(in)          ::funnel
    !-------------------------------------------------
    real(kind=rk) ::r0, draw_pumps_sink3D,width,depth,R2_pump,r2_rel

    r0    =(R_domain-funnel%wall_thickness*0.666_rk)
    depth =funnel%wall_thickness*0.333_rk
    r2_rel=(x-funnel%pump_x_center)**2+(z-R_domain)**2 ! squaring is cheaper then taking the root
    R2_pump=(funnel%pump_diameter*0.5_rk)**2 !squared pump diameter

    draw_pumps_sink3D = smoothstep(r2_rel-(R2_pump-h),h) * soft_bump2(abs(y-R_domain),r0,depth,h)
  end function draw_pumps_sink3D
  !==========================================================================
