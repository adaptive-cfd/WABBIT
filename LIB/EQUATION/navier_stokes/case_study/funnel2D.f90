!#########################################################################################
!                     2D FUNNEL IMPLEMENTATION
!#########################################################################################
subroutine  funnel_penalization2D(penalization, x0, dx, Bs, g ,phi)
  use module_helpers
  implicit none
    ! -----------------------------------------------------------------
    integer(kind=ik), intent(in)  :: Bs, g          !< grid parameter
    real(kind=rk), intent(in)     :: x0(3), dx(3)   !< coordinates of block and block spacinf
    real(kind=rk), intent(in)     :: phi(:,:,:)     !< state vector
    real(kind=rk), intent(inout)  :: penalization(:,:,:) !< reference values of penalized volume
    real(kind=rk), allocatable,save :: mask(:,:,:)    !< mask function
    integer(kind=2), allocatable,save:: mask_color(:,:)!< identifyers of mask parts (plates etc)
    logical                       :: mesh_was_adapted=.true.
    ! -----------------------------------------------------------------
    if (.not. allocated(mask_color))  allocate(mask_color(1:Bs+2*g, 1:Bs+2*g))
    if (.not. allocated(mask))        allocate(mask(1:Bs+2*g, 1:Bs+2*g, 4))

    !!> todo implement function check_if_mesh_adapted (true/false) in adapt mesh
    if ( mesh_was_adapted .eqv. .true. ) then
      ! reset parameters
      mask        = 0.0_rk
      mask_color  = 0
      call draw_funnel2D(x0, dx, Bs, g, mask, mask_color)
      call draw_sponge2D(x0, dx, Bs, g, mask, mask_color)
      call compute_penal2D(mask_color,mask,phi, x0, dx, Bs, g ,penalization)
    else
      call compute_penal2D(mask_color,mask,phi, x0, dx, Bs, g ,penalization)
    end if

end subroutine  funnel_penalization2D




subroutine draw_funnel2D(x0, dx, Bs, g, mask, mask_color)
    implicit none
    ! -----------------------------------------------------------------
    integer(kind=ik), intent(in)             :: Bs, g          !< grid parameter
    real(kind=rk), intent(in)                :: x0(3), dx(3)   !< coordinates of block and block spacinf
    real(kind=rk), intent(inout)             :: mask(:,:,:)    !< mask function
    integer(kind=2), intent(inout), optional :: mask_color(:,:)!< identifyers of mask parts (plates etc)
    ! -----------------------------------------------------------------
    real(kind=rk)     :: x, y, r, h
    real(kind=rk)     :: chi
    integer(kind=ik)  :: ix, iy,n ! loop variables
  ! -----------------------------------------------------------------

    ! parameter for smoothing function (width)
    h  = 1.5_rk*max(dx(1), dx(2))

    ! smooth width in x and y direction
    do iy=g+1, Bs+g
       y = dble(iy-(g+1)) * dx(2) + x0(2)
       r = abs(y-domain_size(2)*0.5_rk)
       do ix=g+1, Bs+g
            x = dble(ix-(g+1)) * dx(1) + x0(1)

            ! plates
            ! ------
            chi = draw_funnel_plates(x,r,funnel,h)
            if (chi>0.0_rk) then
              mask_color(ix,iy)  = color_plates
              mask(ix,iy,2:4)    = mask(ix,iy,2:4) + chi
            endif                                           ! the temperature of the funnel

            ! Walls
            ! -----
            chi = draw_walls(x,r,funnel,h)
            if (chi>0.0_rk) then                       ! default values on the walls
              mask_color(ix,iy)  = color_walls
              mask(ix,iy,2:4)    = mask(ix,iy,2:4) + chi
            endif                                           ! the temperature of the funnel

       end do
    end do

end subroutine  draw_funnel2D




subroutine draw_sponge2D(x0, dx, Bs, g, mask, mask_color)
    implicit none
    ! -----------------------------------------------------------------
    integer(kind=ik), intent(in)             :: Bs, g          !< grid parameter
    real(kind=rk), intent(in)                :: x0(3), dx(3)   !< coordinates of block and block spacinf
    real(kind=rk), intent(inout)             :: mask(:,:,:)    !< mask function
    integer(kind=2), intent(inout), optional :: mask_color(:,:)!< identifyers of mask parts (plates etc)
    ! -----------------------------------------------------------------
    real(kind=rk)     :: x, y, r, h
    real(kind=rk)     :: chi
    integer(kind=ik)  :: ix, iy,n ! loop variables
  ! -----------------------------------------------------------------

    ! parameter for smoothing function (width)
    h  = 1.5_rk*max(dx(1), dx(2))

    ! smooth width in x and y direction
    do iy=g+1, Bs+g
       y = dble(iy-(g+1)) * dx(2) + x0(2)
       r = abs(y-domain_size(2)*0.5_rk)
       do ix=g+1, Bs+g
            x = dble(ix-(g+1)) * dx(1) + x0(1)
                                       ! the temperature of the funnel

            ! Outlet flow: PUMPS
            ! ------------------
            ! pump volume flow
            chi=  draw_pumps_volume_flow(x,r,funnel,h)
            if (chi>0.0_rk) then
              mask(ix,iy,2:3)   = mask(ix,iy,2:3)+chi
              mask_color(ix,iy) = color_pumps
            endif
            ! mass and energy sink
            chi=  draw_pumps_sink(x,r,funnel,h)
            if (chi>0.0_rk) then
              mask_color(ix,iy) = color_pumps_sink
              mask(ix,iy,1) = mask(ix,iy,1)+chi
              mask(ix,iy,4) = mask(ix,iy,4)+chi
            endif

            ! Inlet flow: Capillary
            ! ---------------------
            chi=  draw_jet(x,r,funnel,h)
            if (chi>0.0_rk) then
              mask_color(ix,iy) = color_capillary
              mask(ix,iy,1:4)   =  mask(ix,iy,1:4)+chi
            endif

            ! Outlet flow: Transition to 2pump
            ! ---------------------
            chi=  draw_outlet(x,r,funnel,h)
            !   chi=  draw_sink(x,y,funnel,h)
              if (chi>0.0_rk) then
                mask_color(ix,iy) = color_outlet
                mask(ix,iy,1)     = mask(ix,iy,1)+chi
                mask(ix,iy,3:4)   = mask(ix,iy,3:4)+chi
              endif

       end do
    end do

end subroutine  draw_sponge2D





!> Computes the 2D funnel mask with reference values of the penalized system
subroutine  compute_penal2D(mask_color,mask,phi, x0, dx, Bs, g ,penalization)
    implicit none
    ! -----------------------------------------------------------------
    integer(kind=ik), intent(in)  :: Bs, g          !< grid parameter
    real(kind=rk), intent(in)     :: x0(3), dx(3)   !< coordinates of block and block spacinf
    integer(kind=2), intent(inout):: mask_color(:,:)!< identifyers of mask parts (plates etc)
    real(kind=rk), intent(in)     :: phi(:,:,:)     !< state vector
    real(kind=rk), intent(in)     :: mask(:,:,:)     !< state vector
    real(kind=rk), intent(inout)  :: penalization(:,:,:)  !< funnel penalization term
    ! -----------------------------------------------------------------
    real(kind=rk)     :: x, y, r, h,velocity
    real(kind=rk)     :: rho,chi,v_ref,dq,u,v,p,phi_ref(4),C_inv
    integer(kind=ik)  :: ix, iy,n                                    ! loop variables
    real(kind=rk)     :: velocity_pump,rho_pump,pressure_pump, &    ! outlets and inlets
                      rho_capillary,u_capillary,v_capillary,p_capillary, &
                      p_2nd_pump_stage,rho_2nd_pump_stage
    real(kind=rk)     ::jet_smooth_width,pump_smooth_width  ! smooth width of jet and pumpsponges (boundarylayerthickness)
    ! -----------------------------------------------------------------

    ! initialice parameters
    phi_ref     = 0.0_rk
    u_capillary       =funnel%inlet_velocity(1)
    v_capillary       =funnel%inlet_velocity(2)
    rho_capillary     =funnel%inlet_density
    rho_pump          =funnel%pump_density
    p_capillary       =funnel%inlet_pressure
    velocity_pump     =funnel%pump_speed
    pressure_pump     =funnel%pump_pressure
    p_2nd_pump_stage  =funnel%outlet_pressure
    rho_2nd_pump_stage=funnel%outlet_density

    ! parameter for smoothing function (width)
    h  = 1.5_rk*max(dx(1), dx(2))

    if (3*dx(2)<=0.05_rk*funnel%jet_radius) then
      jet_smooth_width = 0.05_rk*funnel%jet_radius
    else
      jet_smooth_width = 3*dx(2)
      !call abort('ERROR [funnel.f90]: discretication constant dy to large')
    endif

    if (3*dx(1)<=0.1_rk*funnel%pump_diameter) then
      pump_smooth_width = 0.025_rk*funnel%pump_diameter
    else
      pump_smooth_width = 3*h
      !call abort('ERROR [funnel.f90]: discretication constant dy to large')
    endif
    ! smooth width in x and y direction
    do iy=1, Bs+2*g
       y = dble(iy-(g+1)) * dx(2) + x0(2)
       r = abs(y-domain_size(2)*0.5_rk)
       do ix=1, Bs+2*g
            x = dble(ix-(g+1)) * dx(1) + x0(1)
            rho = phi(ix,iy,1)**2
            u   = phi(ix,iy,2)/phi(ix,iy,1)
            v   = phi(ix,iy,3)/phi(ix,iy,1)
            p   = phi(ix,iy,4)

            C_inv=C_sp_inv

            !solid obstacles: walls and plates
            ! ------
            if  (mask_color(ix,iy) == color_plates &
            .or. mask_color(ix,iy) == color_walls ) then
              Phi_ref(2) = 0.0_rk                     ! no velocity in x
              Phi_ref(3) = 0.0_rk                     ! no velocity in y
              Phi_ref(4) = rho*Rs*funnel%temperatur   ! pressure set according to
              C_inv=C_eta_inv
            endif                                           ! the temperature of the funnel

            ! Outlet flow: PUMPS
            ! ------------------
            if (mask_color(ix,iy) == color_pumps) then
              v_ref=velocity_pump*jet_stream(abs(x-funnel%pump_x_center), &
                                             funnel%pump_diameter*0.5_rk,pump_smooth_width)
               Phi_ref(2) = 0
               C_inv=C_eta_inv
              if (y>R_domain) then
                Phi_ref(3) = rho*v_ref
              else
                Phi_ref(3) = -rho*v_ref
              endif
            endif
            ! mass and energy sink
            if (mask_color(ix,iy)==color_pumps_sink) then
              Phi_ref(1) = rho_pump
              Phi_ref(4) = pressure_pump
              C_inv=C_eta_inv
            endif

            ! Inlet flow: Capillary
            ! ---------------------
            if (mask_color(ix,iy)==color_capillary) then
              dq               =jet_stream(r,funnel%jet_radius,jet_smooth_width)
              C_inv=C_sp_inv
              Phi_ref(1) =  rho_capillary
              Phi_ref(2) =  rho_capillary*u_capillary*dq
              Phi_ref(3) =  rho_capillary*v_capillary
              Phi_ref(4) =  p_capillary  !rho*Rs*funnel%temperatur * (1 - dq) + p_capillary * dq
            endif

            ! ! Outlet flow: Transition to 2pump
            ! ! ---------------------
              if (mask_color(ix,iy)==color_outlet) then
                Phi_ref(1) = rho_2nd_pump_stage
                !Phi_ref(2) = 0
                Phi_ref(3) = 0
                Phi_ref(4) = p_2nd_pump_stage
                C_inv=C_sp_inv
              endif

            ! density
            penalization(ix,iy,1)=C_inv*mask(ix,iy,1)*(rho-  Phi_ref(1) )
            ! x-velocity
            penalization(ix,iy,2)=C_inv*mask(ix,iy,2)*(rho*u-  Phi_ref(2) )
            ! y-velocity
            penalization(ix,iy,3)=C_inv*mask(ix,iy,3)*(rho*v-  Phi_ref(3) )
            ! preasure
            penalization(ix,iy,4)=C_inv*mask(ix,iy,4)*(p- Phi_ref(4) )
       end do
    end do
end subroutine  compute_penal2D







!> \brief Compute mask term for stacked rings of the ion funnel
function draw_funnel_plates(x,r,funnel,h)

  real(kind=rk)     , intent(in)  :: x, r, h
  type(type_funnel) , intent(in)  :: funnel
  real(kind=rk)                   :: draw_funnel_plates
  integer(kind=ik)                :: n
  ! distance frommin_inner_diameterter of cylinder
  draw_funnel_plates=0.0_rk

  ! loop over all plates
  do n=1,funnel%nr_plates

      draw_funnel_plates=draw_funnel_plates+draw_plate(x,r,funnel%plate(n),h)

  enddo

end function draw_funnel_plates


!> \brief Compute mask term for single ring-plate
function draw_plate(x,r,plate,h)

  real(kind=rk), intent(in)          :: x, r, h
  type(type_funnel_plate),intent(in) ::plate

  real(kind=rk)                      ::draw_plate,delta_r

  delta_r     = plate%r_out-plate%r_in
  !draw_plate  = soft_bump(x,plate%x0(1),plate%width,h)*soft_bump(r,plate%r_in,delta_r,h)
  draw_plate  = hard_bump(x,plate%x0(1),plate%width)*hard_bump(r,plate%r_in,delta_r)

end function draw_plate


function draw_walls(x,r,funnel,h)

  real(kind=rk),    intent(in)          :: x, r, h
  type(type_funnel),intent(in)          ::funnel

  real(kind=rk)                         ::  mask, draw_walls

  mask=0.0_rk

  ! wall in south and north (i.e. in radial direction)
  if (abs(x-funnel%pump_x_center)> funnel%pump_diameter*0.5_rk) then
    ! mask for r>R_domain-wall_thickness   (outer wall)
         !mask=mask+smoothstep(R_domain-funnel%wall_thickness-r,h)
         mask=mask+hardstep(R_domain-funnel%wall_thickness-r)
  else
        ! +h because the sponge domain should not overlap with the walls
        mask=mask+hardstep(R_domain-0.333_rk*funnel%wall_thickness+h-r)

         !mask=mask+smoothstep(R_domain-0.333_rk*funnel%wall_thickness+h-r,h)
  endif

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

  draw_walls=mask
end function draw_walls

function draw_pumps_volume_flow(x,r,funnel,h)

  real(kind=rk),    intent(in)          :: x, r, h
  type(type_funnel),intent(in)          ::funnel

  real(kind=rk)                         ::  mask, draw_pumps_volume_flow,r0,width

  mask  =0
  width =funnel%wall_thickness*0.333_rk
  r0    =(R_domain-funnel%wall_thickness)
  if (abs(x-funnel%pump_x_center)<= funnel%pump_diameter*0.5_rk) then
         !for r0<r<r0+width apply penalization
         mask=soft_bump2(r,r0,width,h)
  endif


  draw_pumps_volume_flow=mask
end function draw_pumps_volume_flow

function draw_pumps_sink(x,r,funnel,h)

  real(kind=rk),    intent(in)          :: x, r, h
  type(type_funnel),intent(in)          ::funnel

  real(kind=rk)                         ::r0, draw_pumps_sink,width,depth,x_lb,x_rb

  draw_pumps_sink  = 0.0_rk
  r0    =(R_domain-funnel%wall_thickness*0.666_rk)
  depth =funnel%wall_thickness*0.3_rk
  x_lb  =funnel%pump_x_center-funnel%pump_diameter*0.5_rk
  x_rb  =funnel%pump_x_center+funnel%pump_diameter*0.5_rk

  draw_pumps_sink=soft_bump2(r,r0,depth,h)*soft_bump2(x,x_lb,x_rb-x_lb,h)

end function draw_pumps_sink

function draw_jet(x,r,funnel,h)

  real(kind=rk),    intent(in)          :: x, r, h
  type(type_funnel),intent(in)          ::funnel

  real(kind=rk)                         ::draw_jet,length_of_jet

  length_of_jet=funnel%plate(1)%x0(1)-funnel%wall_thickness*0.5_rk-h

  !if (r< funnel%jet_radius) then
    ! wall in EAST
    draw_jet=soft_bump2(x,funnel%wall_thickness*0.5_rk+h,length_of_jet,h)*smoothstep(r-funnel%jet_radius+h,h)

    !draw_jet=smoothstep(x-funnel%wall_thickness,h)-smoothstep(x-funnel%wall_thickness/2.0_rk,h)
  !else
  !endif

end function draw_jet


function draw_outlet(x,r,funnel,h)

  real(kind=rk),    intent(in)          :: x, r, h
  type(type_funnel),intent(in)          ::funnel

  real(kind=rk)                         ::draw_outlet

         ! wall in WEST
    !draw_outlet=smoothstep(domain_size(1)-x-funnel%wall_thickness,h)
draw_outlet=soft_bump2(x,domain_size(1)-funnel%wall_thickness,funnel%wall_thickness*0.5_rk,h)&
            *smoothstep(r-(funnel%min_inner_diameter*0.5_rk-h),h)

end function draw_outlet


function draw_sink(x,r,funnel,h)

  real(kind=rk),    intent(in)          :: x, r, h
  type(type_funnel),intent(in)          ::funnel

  real(kind=rk)                         ::draw_sink,radius

  radius     = sqrt((x-domain_size(1)+funnel%wall_thickness*0.6_rk)**2+r**2)
  draw_sink  = smoothstep(r-funnel%min_inner_diameter*0.4_rk,h)


end function draw_sink




!==========================================================================
    !> \brief initialization of all plates in the funnel
    !> \todo insert picture
    subroutine init_plates(funnel)

      implicit none
      !> geometric parameters of ion funnel
      type(type_funnel), intent(inout) :: funnel

      real(kind=rk)                   :: distance,length,length_focus,width
      integer(kind=ik)                :: n
      type(type_funnel_plate)         :: plate

      allocate(funnel%plate(funnel%nr_plates))

      distance    =funnel%plates_distance
      length      =funnel%length
      width       =funnel%plates_thickness
      ! length of focus area in funnel (+distance +wdith because the last plate is in the orifice)
      length_focus           = length +width+distance   - (funnel%max_inner_diameter-funnel%min_inner_diameter)/funnel%slope
      ! origin of funnel
      funnel%offset=(/ domain_size(1)-length-funnel%wall_thickness-distance, &
                       R_domain/)
      if(funnel%offset(1)<funnel%wall_thickness) then
       call abort(13457,'Error [module_mask.f90]: your funnel is to long')
      endif

      ! initialicd all plates
      ! ---------------------
      ! first plate is often fat and ugly (:
      funnel%plate(1)%x0    =funnel%offset                  ! coordinates of the midpoint of the funnel
      funnel%plate(1)%width =funnel%first_plate_thickness   ! as the name says its a fat plate
      funnel%plate(1)%r_in  =funnel%max_inner_diameter/2    ! inner diameter of the ring
      funnel%plate(1)%r_out  =funnel%outer_diameter/2        ! outer diameter of the ring

      ! all the other plates are similar
      do n=2,funnel%nr_plates
        plate%x0(1)     = funnel%plate(n-1)%x0(1)+funnel%plate(n-1)%width+distance
        plate%x0(2)     = funnel%offset(2)
        plate%width     = width
        if (plate%x0(1)-funnel%offset(1)<length_focus) then
           plate%r_in    = funnel%max_inner_diameter/2
        else
           plate%r_in = plate%r_in - funnel%slope*(distance+width)/2
        endif
        plate%r_out   =funnel%outer_diameter/2
        funnel%plate(n)=plate
      enddo

    end subroutine init_plates
!==========================================================================



!==========================================================================
  !> Integrates the flow field close to the pump outlet,
  !> and computes the area of the intagration region.
  subroutine integrate_over_pump_area2D(u,g,Bs,x0,dx,integral,area)
      implicit none

      !> grid parameter (g ghostnotes,Bs Bulk)
      integer(kind=ik), intent(in)                     :: Bs, g
      !> density,pressure
      real(kind=rk), dimension(1:,1:,1:), intent(in)   :: u
      !> spacing and origin of block
      real(kind=rk), dimension(2), intent(in)          :: x0, dx
      !> mean density
      real(kind=rk),intent(out)                      :: integral(4), area

      real(kind=rk)                                    :: h,r,y,x,r0,width
      !temporal data field
      real(kind=rk),dimension(5)                       :: tmp

      integer(kind=ik)                                 :: ix,iy

       h  = 1.5_rk*max(dx(1), dx(2))
      ! calculate mean density close to the pump
      width =funnel%wall_thickness
      tmp   =  0.0_rk
      r0    =(R_domain-2*funnel%wall_thickness)
       do iy=g+1, Bs+g
         y = dble(iy-(g+1)) * dx(2) + x0(2)
         r = abs(y-R_domain)
         do ix=g+1, Bs+g
              x = dble(ix-(g+1)) * dx(1) + x0(1)
              if (abs(x-funnel%pump_x_center)<= funnel%pump_diameter*0.5_rk .and. &
                  r>r0 .and. r<r0+width) then
                tmp(1:4)  = tmp(1:4)+ u(ix,iy,:)
                tmp(5)    = tmp(5)  + 1.0_rk
              endif
          enddo
        enddo
        integral  = integral + tmp(1:4) *dx(1)*dx(2)
        area      = area     + tmp(5)   *dx(1)*dx(2)

  end subroutine integrate_over_pump_area2D
  !==========================================================================
