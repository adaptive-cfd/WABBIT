
!> \brief reads parameters for mask function from file
subroutine init_funnel(FILE)

  implicit none

 ! character(len=*), intent(in) :: filename
  type(inifile) , intent(inout) :: FILE

  real(kind=rk)                 :: dmax,dmin
  integer(kind=ik)              :: nr_focus_plates
   ! inifile structure
  !type(inifile) :: FILE
  !call read_ini_file_mpi(FILE, filename, .true.)

  !===========================================================================================================
  ! READ IN geometry
  ! ----------------
  call read_param_mpi(FILE, 'funnel', 'outer_diameter'        , funnel%outer_diameter, R_domain*0.5_rk )
  call read_param_mpi(FILE, 'funnel', 'maximal_inner_diameter', dmax, domain_size(2)/3.0_rk )
  call read_param_mpi(FILE, 'funnel', 'minimal_inner_diameter', dmin, domain_size(2)/4.0_rk )
  call read_param_mpi(FILE, 'funnel', 'Number_of_plates'      , funnel%nr_plates, 30 )
  call read_param_mpi(FILE, 'funnel', 'Number_of_focus_plates', nr_focus_plates, 15)
  call read_param_mpi(FILE, 'funnel', 'Temperatur_of_plates'  , funnel%temperatur, 300.0_rk)
  call read_param_mpi(FILE, 'funnel', 'pump_diameter'         , funnel%pump_diameter, domain_size(1)/5.0_rk)
  call read_param_mpi(FILE, 'funnel', 'jet_diameter'          , funnel%jet_radius, R_domain*0.5_rk)
  call read_param_mpi(FILE, 'funnel', 'pump_x_center'         , funnel%pump_x_center, domain_size(1)*0.5_rk)
  funnel%max_inner_diameter   = dmax
  funnel%min_inner_diameter   = dmin
  funnel%wall_thickness       = 0.05*domain_size(1)
  funnel%length               = domain_size(1)*0.95_rk-funnel%wall_thickness*2.0_rk
  funnel%plates_thickness     = funnel%length/(2.0_rk*funnel%nr_plates)
  funnel%plates_distance      = (funnel%length-funnel%plates_thickness)/(funnel%nr_plates-1)
  funnel%slope                = (dmax - dmin)/((nr_focus_plates-1)*funnel%plates_distance)
  funnel%jet_radius           = funnel%jet_radius/2.0_rk
  funnel%pump_density         = 0
  !===========================================================================================================
  ! READ IN Capillary inlet flow
  ! ----------------------------
  funnel%inlet_velocity=(/ 0.0, 0.0 /)
  call read_param_mpi(FILE, 'funnel', 'inlet_velocity'  , funnel%inlet_velocity ,funnel%inlet_velocity)
  call read_param_mpi(FILE, 'funnel', 'inlet_density'   , funnel%inlet_density  , 1.0_rk )
  call read_param_mpi(FILE, 'funnel', 'inlet_pressure'  , funnel%inlet_pressure, 1.0_rk )
  call read_param_mpi(FILE, 'funnel', 'pump_speed'      , funnel%pump_speed, 30.0_rk )
  call read_param_mpi(FILE, 'funnel', 'outlet_pressure' , funnel%outlet_pressure, 1.0_rk)

  if (funnel%length         >domain_size(1)-2.0_rk*funnel%wall_thickness .or. &
      funnel%outer_diameter >domain_size(2)-2.0_rk*funnel%wall_thickness) then
    call abort(5032,"ERROR [funnel.f90]:funnel is larger then simulation domain!")
  endif

  !initialice geometry of ion funnel plates
  call init_plates(funnel)
end subroutine init_funnel


subroutine add_funnel(penalization, x0, dx, Bs, g ,phi)


    implicit none

    ! grid
    integer(kind=ik), intent(in)                     :: Bs, g
    !> penalization term including mask
    real(kind=rk), dimension(:,:,:), intent(inout)   :: penalization
    !> spacing and origin of block
    real(kind=rk), dimension(2), intent(in)          :: x0, dx
    !> statevector
    real(kind=rk), dimension(:,:,:), intent(in)      :: phi

    !> mask of geometrie
    real(kind=rk)                                    :: mask(Bs+2*g,Bs+2*g,4),phi_ref(Bs+2*g,Bs+2*g,4)
    ! auxiliary variables
    real(kind=rk)                                    :: x, y, r, h,velocity
    ! preasure,density velocities
    real(kind=rk)                                    :: p,rho,u,v,chi,v_ref,dq
    ! loop variables
    integer(kind=ik)                                 :: ix, iy,n
    ! outlets and inlets
    real(kind=rk)                                     :: velocity_pump,rho_pump,pressure_pump, &
                                                         rho_capillary,u_capillary,v_capillary,p_capillary, &
                                                         p_2nd_pump_stage,rho_2nd_pump_stage
    ! smooth width of jet
    real(kind=rk)                                     ::jet_smooth_width,pump_smooth_width
!---------------------------------------------------------------------------------------------
! variables initialization
    if (size(penalization,1) /= Bs+2*g) call abort(777109,"wrong array size, there's pirates, captain!")

    ! reset mask array
    mask = 0.0_rk
    phi_ref=0.0_rk
!---------------------------------------------------------------------------------------------
! main body

    u_capillary       =funnel%inlet_velocity(1)
    v_capillary       =funnel%inlet_velocity(2)
    rho_capillary     =funnel%inlet_density
    rho_pump          =funnel%pump_density
    rho_2nd_pump_stage=funnel%inlet_density/100.0_rk
    p_capillary       =funnel%inlet_pressure
    velocity_pump     =funnel%pump_speed
    pressure_pump     =funnel%pump_pressure
    p_2nd_pump_stage  =funnel%outlet_pressure



    ! parameter for smoothing function (width)
    h  = 1.5_rk*max(dx(1), dx(2))

    if (3*dx(2)<=0.05_rk*funnel%jet_radius) then
      jet_smooth_width = 0.05_rk*funnel%jet_radius
    else
      jet_smooth_width = 2*h
      !call abort('ERROR [funnel.f90]: discretication constant dy to large')
    endif

    if (3*dx(1)<=0.1_rk*funnel%pump_diameter) then
      pump_smooth_width = 0.05_rk*funnel%pump_diameter
    else
      pump_smooth_width = 3*h
      !call abort('ERROR [funnel.f90]: discretication constant dy to large')
    endif



    ! smooth width in x and y direction
    do iy=1, Bs+2*g
       y = dble(iy-(g+1)) * dx(2) + x0(2)
       call continue_periodic(y,domain_size(2))
       r = abs(y-domain_size(2)*0.5_rk)
       do ix=1, Bs+2*g
            x = dble(ix-(g+1)) * dx(1) + x0(1)
            call continue_periodic(x,domain_size(1))
            rho         = phi(ix,iy,1)**2
            u           = phi(ix,iy,2)/phi(ix,iy,1)
            v           = phi(ix,iy,3)/phi(ix,iy,1)
            p           = phi(ix,iy,4)


            ! Funnel
            ! ------
            ! 1. compute mask term:
            chi = draw_funnel_plates(x,r,funnel,h)
            ! 2. set quantities to desired values:
            if (chi>0) then
              mask(ix,iy,2:4) = mask(ix,iy,2:4) + chi
              Phi_ref(ix,iy,2) = 0.0_rk                     ! no velocity in x
              Phi_ref(ix,iy,3) = 0.0_rk                     ! no velocity in y
              Phi_ref(ix,iy,4) = rho*Rs*funnel%temperatur   ! pressure set according to
            endif                                           ! the temperature of the funnel

            ! Walls
            ! -----
            ! draw the walls arround the funnel
            chi = draw_walls(x,r,funnel,h)
            ! set reference values at these region
            if (chi>0) then                       ! default values on the walls
              mask(ix,iy,2:4) = mask(ix,iy,2:4) + chi
              Phi_ref(ix,iy,2) = 0.0_rk                     ! no velocity in x
              Phi_ref(ix,iy,3) = 0.0_rk                     ! no velocity in y
              Phi_ref(ix,iy,4) = rho*Rs*funnel%temperatur   ! pressure set according to
            endif                                           ! the temperature of the funnel

            ! Outlet flow: PUMPS
            ! ------------------
            ! pump volume flow
            ! compute mask
            chi=  draw_pumps_volume_flow(x,r,funnel,h)
            if (chi>0) then
              mask(ix,iy,2:3) = mask(ix,iy,2:3)+chi
              !compute velocity profile
              v_ref=velocity_pump*jet_stream(abs(x-funnel%pump_x_center),funnel%pump_diameter*0.5_rk,pump_smooth_width)
               Phi_ref(ix,iy,2) = 0
              if (y>R_domain) then
                Phi_ref(ix,iy,3) = rho*v_ref
              else
                Phi_ref(ix,iy,3) = -rho*v_ref
              endif
            endif

            ! mass and energy sink
            chi=  draw_pumps_sink(x,r,funnel,h)
            if (chi>0) then
              mask(ix,iy,1) = mask(ix,iy,1)+chi
              mask(ix,iy,4) = mask(ix,iy,4)+chi
              Phi_ref(ix,iy,1) = rho_pump
              Phi_ref(ix,iy,4) = pressure_pump
            endif



            ! Inlet flow: Capillary
            ! ---------------------
            chi=  draw_jet(x,y,funnel,h)
            if (chi>0) then
              dq               =jet_stream(r,funnel%jet_radius,jet_smooth_width)

              mask(ix,iy,1:4)  =  mask(ix,iy,1:4)+chi
              Phi_ref(ix,iy,1) =  rho_capillary
              Phi_ref(ix,iy,2) =  rho_capillary*u_capillary*dq
              Phi_ref(ix,iy,3) =  rho_capillary*v_capillary
              Phi_ref(ix,iy,4) =  p_capillary  !rho*Rs*funnel%temperatur * (1 - dq) + p_capillary * dq
            endif


            ! ! Outlet flow: Transition to 2pump
            ! ! ---------------------
            chi=  draw_outlet(x,y,funnel,h)
            !   chi=  draw_sink(x,y,funnel,h)

              if (chi>0) then
                mask(ix,iy,1) = mask(ix,iy,1)+chi
                mask(ix,iy,4) = mask(ix,iy,4)+chi
                Phi_ref(ix,iy,1) = rho_2nd_pump_stage
                !Phi_ref(ix,iy,2) = 0
                Phi_ref(ix,iy,3) = 0
                Phi_ref(ix,iy,4) = p_2nd_pump_stage
              endif


            ! density
            penalization(ix,iy,1)=C_eta_inv*mask(ix,iy,1)*(rho-  Phi_ref(ix,iy,1) )
            ! x-velocity
            penalization(ix,iy,2)=C_eta_inv*mask(ix,iy,2)*(rho*u-  Phi_ref(ix,iy,2) )
            ! y-velocity
            penalization(ix,iy,3)=C_eta_inv*mask(ix,iy,3)*(rho*v-  Phi_ref(ix,iy,3) )
            ! preasure
            penalization(ix,iy,4)=C_eta_inv*mask(ix,iy,4)*(p-  Phi_ref(ix,iy,4) )

       end do
    end do
end subroutine add_funnel


subroutine integrate_over_pump_area(u,g,Bs,x0,dx,integral,area)


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
    width =funnel%wall_thickness*0.5_rk
    tmp   =  0
    r0    =(R_domain-funnel%wall_thickness)
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

end subroutine integrate_over_pump_area



subroutine mean_quantity(integral,area)
    !> area of taking the mean
    real(kind=rk),intent(in)    :: area
    !> integral over the area
    real(kind=rk),intent(inout) :: integral(1:)

    ! temporary values
    real(kind=rk),allocatable   :: tmp(:)
    real(kind=rk)               :: A
    integer(kind=ik)            :: mpierr,Nq


    Nq = size(integral,1)
    allocate(tmp(Nq))

    tmp=integral

    ! integrate over all procs
    call MPI_ALLREDUCE(tmp  ,integral, Nq , MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
    call MPI_ALLREDUCE(area ,A       , 1  , MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)

    if ( abs(A) <= 1.0e-13_rk) then
      call abort(24636,"Error [funnel.f90]: only chuck norris can devide by zero!!")
    endif

    !devide by the area of the region
    integral = integral / A

    funnel%pump_density = integral(1)/2
    funnel%pump_pressure = integral(4)/2

end subroutine mean_quantity





!==========================================================================
subroutine draw_funnel(mask, x0, dx, Bs, g )
    implicit none
    ! grid
    integer(kind=ik), intent(in)                              :: Bs, g
    !> mask term for every grid point of this block
    real(kind=rk), dimension(:,:), intent(out)                :: mask
    !> spacing and origin of block
    real(kind=rk), dimension(2), intent(in)                   :: x0, dx

    ! auxiliary variables
    real(kind=rk)                                             :: x, y, r, h
    ! loop variables
    integer(kind=ik)                                          :: ix, iy
!-------------------------------------------------------------------------
! variables initialization
    if (size(mask,1) /= Bs+2*g) call abort(777109,"wrong array size, there's pirates, captain!")

    ! reset mask array
    mask = 0.0_rk

    ! parameter for smoothing function (width)
    h = 1.5_rk*max(dx(1), dx(2))
    do iy=1, Bs+2*g
       y = dble(iy-(g+1)) * dx(2) + x0(2)
       r = abs(y-domain_size(2)*0.5_rk)
       do ix=1, Bs+2*g
           x = dble(ix-(g+1)) * dx(1) + x0(1)
            mask(ix,iy) =   draw_walls(x,r,funnel,h) + &
                            draw_funnel_plates(x,r,funnel,h)
       end do
    end do
end subroutine draw_funnel
!==========================================================================





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
  draw_plate  = soft_bump(x,plate%x0(1),plate%width,h)*soft_bump(r,plate%r_in,delta_r,h)
end function draw_plate


function draw_walls(x,r,funnel,h)

  real(kind=rk),    intent(in)          :: x, r, h
  type(type_funnel),intent(in)          ::funnel

  real(kind=rk)                         ::  mask, draw_walls

  mask=0.0_rk

  if (abs(x-funnel%pump_x_center)> funnel%pump_diameter*0.5_rk) then
    ! mask for r>R_domain-wall_thickness   (outer wall)
         mask=mask+smoothstep(R_domain-funnel%wall_thickness-r,h)
  else

         mask=mask+smoothstep(R_domain-0.333_rk*funnel%wall_thickness-r,h)
  endif

  ! wall in EAST
  if (  r > funnel%jet_radius  ) then
         mask=mask+smoothstep(x-funnel%wall_thickness,h)
  else
         mask=mask+smoothstep(x-funnel%wall_thickness*0.5_rk,h)
  endif

  ! wall in WEST
  if ( r > funnel%min_inner_diameter*0.5_rk) then
         mask=mask+smoothstep(domain_size(1)-x-funnel%wall_thickness,h)
  else
         mask=mask+smoothstep(domain_size(1)-x-funnel%wall_thickness*0.5_rk,h)
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

  real(kind=rk)                         ::r0,  mask, draw_pumps_sink,width

  mask  = 0.0_rk
  r0    =(R_domain-funnel%wall_thickness*0.6666_rk)
  width =funnel%wall_thickness*0.333_rk
  if (abs(x-funnel%pump_x_center)<= funnel%pump_diameter/2) then
         !mask=soft_bump2(r,r0,width,h)
         mask=smoothstep(r0-r,h)
  endif


  draw_pumps_sink=mask
end function draw_pumps_sink

function draw_jet(x,y,funnel,h)

  real(kind=rk),    intent(in)          :: x, y, h
  type(type_funnel),intent(in)          ::funnel

  real(kind=rk)                         ::draw_jet

  if (abs(y-R_domain)< funnel%jet_radius) then
    ! wall in EAST
    draw_jet=smoothstep(x-funnel%wall_thickness,h)-smoothstep(x-funnel%wall_thickness/2,h)
  else
    draw_jet=0.0_rk
  endif

end function draw_jet


function draw_outlet(x,y,funnel,h)

  real(kind=rk),    intent(in)          :: x, y, h
  type(type_funnel),intent(in)          ::funnel

  real(kind=rk)                         ::draw_outlet

 if (abs(y-R_domain)< funnel%min_inner_diameter/2) then
         ! wall in WEST
    !draw_outlet=smoothstep(domain_size(1)-x-funnel%wall_thickness,h)
    draw_outlet=soft_bump(x,domain_size(1)-funnel%wall_thickness+h,funnel%wall_thickness*0.5_rk-2*h,h)

  else
    draw_outlet=0.0_rk
  endif

end function draw_outlet


function draw_sink(x,y,funnel,h)

  real(kind=rk),    intent(in)          :: x, y, h
  type(type_funnel),intent(in)          ::funnel

  real(kind=rk)                         ::draw_sink,r

  r=sqrt((x-domain_size(1)+funnel%wall_thickness*0.6_rk)**2+(y-R_domain)**2)
    draw_sink=smoothstep(r-funnel%min_inner_diameter*0.4_rk,h)


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
      ! length of focus area in funnel
      length_focus           = length - (funnel%max_inner_diameter-funnel%min_inner_diameter)/funnel%slope
      ! origin of funnel
      funnel%offset=(/ domain_size(1)-length-funnel%wall_thickness-distance+width, &
                       R_domain/)
      if(funnel%offset(1)<funnel%wall_thickness) then
       call abort(13457,'Error [module_mask.f90]: your funnel is to long')
      endif

      ! initialicd all plates
      do n=1,funnel%nr_plates
        plate%x0(1)     = funnel%offset(1)+(n-1)*distance
        plate%x0(2)     = funnel%offset(2)
        plate%width     = width
        if (plate%x0(1)-funnel%offset(1)<length_focus) then
           plate%r_in    = funnel%max_inner_diameter/2
        else
           plate%r_in = plate%r_in - funnel%slope*distance/2
        endif
        plate%r_out   =funnel%outer_diameter/2
        funnel%plate(n)=plate
      enddo

    end subroutine init_plates
!==========================================================================
