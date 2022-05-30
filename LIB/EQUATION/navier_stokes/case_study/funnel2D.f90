!#########################################################################################
!                     2D FUNNEL IMPLEMENTATION
!#########################################################################################
subroutine  funnel_penalization2D(Bs, g, x_0, delta_x, phi, mask, phi_ref)
  use module_helpers
  implicit none
    ! -----------------------------------------------------------------
    integer(kind=ik), intent(in)  :: g          !< grid parameter
    integer(kind=ik), dimension(3), intent(in) :: Bs
    real(kind=rk), intent(in)     :: x_0(2), delta_x(2)   !< coordinates of block and block spacinf
    real(kind=rk), intent(in)     :: phi(:,:,:)     !< state vector
    real(kind=rk), intent(inout)  :: phi_ref(:,:,:) !< reference values of penalized volume
    real(kind=rk), intent(inout)  :: mask(:,:,:)    !< mask function
    integer(kind=2), allocatable,save:: mask_color(:,:)!< identifyers of mask parts (plates etc)
    logical                       :: mesh_was_adapted=.true.
    real(kind=rk)                 :: x0, dx, r0, dr
    ! -----------------------------------------------------------------
    call cartesian2cylinder( x_0, delta_x, dx, dr, x0, r0 )

    if (.not. allocated(mask_color))  allocate(mask_color(1:Bs(1)+2*g, 1:Bs(2)+2*g))
    !!> todo implement function check_if_mesh_adapted (true/false) in adapt mesh
    if ( mesh_was_adapted .eqv. .true. ) then
      ! dont switch the order of draw_funnel3D and draw_sponge3D,
      ! because mask and color are reset in the draw_funnel
      call draw_funnel2D(r0, dr, x0, dx, Bs, g, mask, mask_color)
      call draw_sponge2D(r0, dr, x0, dx, Bs, g, mask, mask_color)
    end if

    call compute_penal2D(mask_color, mask, phi, r0, dr, x0, dx, Bs, g ,phi_ref)

end subroutine  funnel_penalization2D



subroutine draw_funnel2D(r0, dr, x0, dx, Bs, g, mask, mask_color)
    implicit none
    ! -----------------------------------------------------------------
    integer(kind=ik), intent(in)             :: g          !< grid parameter
    integer(kind=ik), dimension(3), intent(in) :: Bs
    real(kind=rk), intent(in)                :: r0, x0, dr, dx !< coordinates of block and block spacinf
    real(kind=rk), intent(inout)             :: mask(:,:,:)    !< mask function
    integer(kind=2), intent(inout), optional :: mask_color(:,:)!< identifyers of mask parts (plates etc)
    ! -----------------------------------------------------------------
    real(kind=rk)     :: x, r, h
    real(kind=rk)     :: chi
    integer(kind=ik)  :: ix, ir,n ! loop variables
  ! -----------------------------------------------------------------

    ! parameter for smoothing function (width)
    h  = 1.5_rk*max(dx, dr)
    ! smooth width in x and y direction
    do ir = g+1, Bs(2)+g+ONE_SKIPREDUNDANT
       ! calculate the radial distance from the coordinate origin
       r = abs(dble(ir-(g+1)) * dr + r0)

       do ix = g+1, Bs(1)+g+ONE_SKIPREDUNDANT
            x = dble(ix-(g+1)) * dx + x0
            !=============================
!     /\     reset the mask function!
!    /  \    caution! Reseting should only be done once
!   /    \   for optimal performance
!  / stop \
! +--------+
            mask_color(ix,ir)=0
            mask(ix,ir,:)=0.0_rk
            !=============================

            ! plates
            ! ------
            chi = draw_funnel_plates(x,r,funnel,h)
            if (chi>0.0_rk) then
              mask_color(ix,ir)  = color_plates
              mask(ix,ir,2:4)    = mask(ix,ir,2:4) + chi
            endif                                           ! the temperature of the funnel

            ! Walls
            ! -----
             chi = draw_walls(x,r,funnel,h)
             if (chi>0.0_rk) then                       ! default values on the walls
               mask_color(ix,ir)  = color_walls
               mask(ix,ir,2:4)    = mask(ix,ir,2:4) + chi
             endif                                           ! the temperature of the funnel

       end do
    end do

end subroutine  draw_funnel2D




subroutine draw_sponge2D(r0, dr, x0, dx, Bs, g, mask, mask_color)
    implicit none
    ! -----------------------------------------------------------------
    integer(kind=ik), intent(in)             :: g          !< grid parameter
    integer(kind=ik), dimension(3), intent(in) :: Bs     !< grid parameter
    real(kind=rk), intent(in)                :: r0, x0, dr, dx !< coordinates of block and block spacing
    real(kind=rk), intent(inout)             :: mask(:,:,:)    !< mask function
    integer(kind=2), intent(inout), optional :: mask_color(:,:)!< identifyers of mask parts (plates etc)
    ! -----------------------------------------------------------------
    real(kind=rk)     :: x, r, h
    real(kind=rk)     :: chi
    integer(kind=ik)  :: ix, ir, n ! loop variables
  ! -----------------------------------------------------------------

    ! parameter for smoothing function (width)
    h  = 1.5_rk*max(dx, dr)

    ! smooth width in x and y direction
    do ir=g+1, Bs(2)+g+ONE_SKIPREDUNDANT
      ! calculate the radial distance from the coordinate origin
      r = abs(dble(ir-(g+1)) * dr + r0)

       do ix=g+1, Bs(1)+g+ONE_SKIPREDUNDANT
            x = dble(ix-(g+1)) * dx + x0
                                       ! the temperature of the funnel
            ! Outlet flow: PUMPS
            ! ------------------
            ! pump volume flow
            chi=  draw_pumps_volume_flow(x,r,funnel,h)
            if (chi>0.0_rk) then
              mask(ix,ir,2:3)   = mask(ix,ir,2:3)+chi
              mask_color(ix,ir) = color_pumps
            endif
            ! mass and energy sink
            chi=  draw_pumps_sink(x,r,funnel,h)
            if (chi>0.0_rk) then
              mask_color(ix,ir) = color_pumps_sink
              mask(ix,ir,1) = mask(ix,ir,1)+chi
              mask(ix,ir,4) = mask(ix,ir,4)+chi
            endif

            ! Inlet flow: Capillary
            ! ---------------------
            chi=  draw_jet(x,r,funnel,h)
            if (chi>0.0_rk) then
              mask_color(ix,ir) = color_capillary
              mask(ix,ir,1:4)   =  mask(ix,ir,1:4)+chi
            endif

            ! Outlet flow: Transition to 2pump
            ! ---------------------
            chi=  draw_outlet(x,r,funnel,h)
            !   chi=  draw_sink(x,y,funnel,h)
              if (chi>0.0_rk) then
                mask_color(ix,ir) = color_outlet
                mask(ix,ir,1)     = mask(ix,ir,1)+chi
                mask(ix,ir,3:4)   = mask(ix,ir,3:4)+chi
              endif

       end do
    end do

end subroutine  draw_sponge2D





!> Computes the 2D funnel mask with reference values of the penalized system
subroutine  compute_penal2D(mask_color, mask, phi, r0, dr, x0, dx, Bs, g ,phi_ref)
    implicit none
    ! -----------------------------------------------------------------
    integer(kind=ik), intent(in)  :: g          !< grid parameter
    integer(kind=ik), dimension(3), intent(in) :: Bs
    real(kind=rk), intent(in)     :: r0, x0, dr, dx !< coordinates of block and block spacing
    integer(kind=2), intent(inout):: mask_color(:,:)!< identifyers of mask parts (plates etc)
    real(kind=rk), intent(in)     :: phi(:,:,:)     !< state vector
    real(kind=rk), intent(inout)  :: mask(:,:,:)     !< mask
    real(kind=rk), intent(inout)  :: phi_ref(:,:,:)  !< funnel penalization term
    ! -----------------------------------------------------------------
    real(kind=rk)     :: x, y, r, h,velocity
    real(kind=rk)     :: rho,chi,v_ref,dq,u,v,p,C_inv
    integer(kind=ik)  :: ix, ir,n                                    ! loop variables
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
    h  = 1.5_rk*max(dx, dr)
    if (3*dr<=0.05_rk*funnel%jet_radius) then
      jet_smooth_width = 0.05_rk*funnel%jet_radius
    else
      jet_smooth_width = 3*dr
      !call abort('ERROR [funnel.f90]: discretication constant dy to large')
    endif

    if (3*dx <=0.1_rk*funnel%pump_diameter) then
      pump_smooth_width = 0.025_rk*funnel%pump_diameter
    else
      pump_smooth_width = 3.0_rk*h
      !call abort('ERROR [funnel.f90]: discretication constant dy to large')
    endif
    ! smooth width in x and y direction
    do ir=g+1, Bs(2)+g+ONE_SKIPREDUNDANT
      ! calculate the radial distance from the coordinate origin
       r = abs(dble(ir-(g+1)) * dr + r0)

       do ix=g+1, Bs(1)+g+ONE_SKIPREDUNDANT
            x = dble(ix-(g+1)) * dx + x0
            rho = phi(ix,ir,rhoF)
            u   = phi(ix,ir,UxF)
            v   = phi(ix,ir,UyF)
            p   = phi(ix,ir,pF)

            C_inv=C_sp_inv

            !solid obstacles: walls and plates
            ! ------
            if  (mask_color(ix,ir) == color_plates &
            .or. mask_color(ix,ir) == color_walls ) then
              Phi_ref(ix,ir,2) = 0.0_rk                     ! no velocity in x
              Phi_ref(ix,ir,3) = 0.0_rk                     ! no velocity in y
              Phi_ref(ix,ir,4) = rho*Rs*funnel%temperatur   ! pressure set according to
              C_inv=C_eta_inv
            endif                                           ! the temperature of the funnel

            ! Outlet flow: PUMPS
            ! ------------------
            if (mask_color(ix,ir) == color_pumps) then
              v_ref=velocity_pump*jet_stream(abs(x-funnel%pump_x_center), &
                                             funnel%pump_diameter*0.5_rk,pump_smooth_width)
               Phi_ref(ix,ir,UxF) = 0
               C_inv=C_eta_inv
              if (r0 > 0.0_rk) then
                Phi_ref(ix,ir,UyF) = rho*v_ref
              else
                Phi_ref(ix,ir,UyF) = -rho*v_ref
              endif
            endif
            ! mass and energy sink
            if (mask_color(ix,ir)==color_pumps_sink) then
              Phi_ref(ix,ir,1) = rho_pump
              Phi_ref(ix,ir,4) = pressure_pump
              C_inv=C_eta_inv!C_sp_inv
            endif

            ! Inlet flow: Capillary
            ! ---------------------
            if (mask_color(ix,ir)==color_capillary) then
              dq               =jet_stream(r,funnel%jet_radius,jet_smooth_width)
              C_inv=C_sp_inv
              Phi_ref(ix,ir,1) =  rho_capillary
              Phi_ref(ix,ir,2) =  rho_capillary*u_capillary*dq
              Phi_ref(ix,ir,3) =  rho_capillary*v_capillary
              Phi_ref(ix,ir,4) =  p_capillary  !rho*Rs*funnel%temperatur * (1 - dq) + p_capillary * dq
            endif

            ! Outlet flow: Transition to 2pump
            ! ---------------------
              if (mask_color(ix,ir)==color_outlet) then
                Phi_ref(ix,ir,1) = rho_2nd_pump_stage
                !Phi_ref(ix,ir,2) = 0
                Phi_ref(ix,ir,3) = 0
                Phi_ref(ix,ir,4) = p_2nd_pump_stage
                C_inv=C_sp_inv
              endif
              ! add penalization strength to mask
              mask(ix,ir,:)=C_inv*mask(ix,ir,:)
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

   ! is needed because mask of walls overlap
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
      integer(kind=ik), intent(in)                     :: g
      integer(kind=ik), dimension(3), intent(in) :: Bs
      !> density,pressure
      real(kind=rk), dimension(1:,1:,1:), intent(in)   :: u
      !> spacing and origin of block
      real(kind=rk), dimension(2), intent(in)          :: x0, dx
      !> mean density
      real(kind=rk),intent(out)                      :: integral(4), area

      real(kind=rk)                                    :: h,r,y,x,r0,width
      !temporal data field
      real(kind=rk),dimension(5)                       :: tmp

      integer(kind=ik)                                 :: ix,ir

       h  = 1.5_rk*max(dx(1), dx(2))
      ! calculate mean density close to the pump
      width =funnel%wall_thickness
      tmp   =  0.0_rk
      r0    =(R_domain-2*funnel%wall_thickness)
       do ir=g+1, Bs(2)+g+ONE_SKIPREDUNDANT
         y = dble(ir-(g+1)) * dx(2) + x0(2)
         r = abs(y-R_domain)
         do ix=g+1, Bs(1)+g+ONE_SKIPREDUNDANT
              x = dble(ix-(g+1)) * dx(1) + x0(1)
              if (abs(x-funnel%pump_x_center)<= funnel%pump_diameter*0.5_rk .and. &
                  r>r0 .and. r<r0+width) then
                tmp(1:4)  = tmp(1:4)+ u(ix,ir,:)
                tmp(5)    = tmp(5)  + 1.0_rk
              endif
          enddo
        enddo
        integral  = integral + tmp(1:4) *dx(1)*dx(2)
        area      = area     + tmp(5)   *dx(1)*dx(2)

  end subroutine integrate_over_pump_area2D
  !==========================================================================


!########################################################################
!> Transform computational coordinates to cylindrical coords
subroutine cartesian2cylinder(x_0, delta_x, dx, dr, x0, r0)

  implicit none
  !--------------------------------------------------------
  real(kind=rk), intent(in) :: x_0(2), delta_x(2)     !<cartesian coords
  real(kind=rk), intent(out) :: dx, dr, r0, x0  !<cylindrical coords
  !--------------------------------------------------------
  ! Set the lattice spacing and coordinate origin:
  x0 = x_0(1)
  dx = delta_x(1)
  dr = delta_x(2)
  ! The coordinate origin is dependent on the coordinate system
  if ( params_ns%coordinates=="cylindrical" ) then
    ! The total grid is shifted by R_min, which accounts for the
    ! infinitesimal cylinder centered arround the symmetrie axis.
    ! The origin is therefore shifted to (x,r) = (0, R_min)
    r0 = x_0(2) + params_ns%R_min
  else
    ! The coordinate system is centered at (x,r) = (0, L_r/2) and -L_r/2<r<L_r/2
    r0 = x_0(2) - params_ns%domain_size(2)*0.5_rk
  end if

end subroutine cartesian2cylinder
!########################################################################
