
!> \brief reads parameters for mask function from file
subroutine init_funnel(FILE)

  implicit none
    
 ! character(len=*), intent(in) :: filename
  type(inifile) , intent(inout) :: FILE
   ! inifile structure
  !type(inifile) :: FILE
  !call read_ini_file_mpi(FILE, filename, .true.)

  call read_param_mpi(FILE, 'funnel', 'outer_diameter'        , funnel%outer_diameter, domain_size(2)/2.0_rk )
  call read_param_mpi(FILE, 'funnel', 'maximal_inner_diameter', funnel%max_inner_diameter, domain_size(2)/3.0_rk )
  call read_param_mpi(FILE, 'funnel', 'minimal_inner_diameter', funnel%min_inner_diameter, domain_size(2)/4.0_rk )
  call read_param_mpi(FILE, 'funnel', 'Number_of_plates'      , funnel%nr_plates, 30 )
  call read_param_mpi(FILE, 'funnel', 'diameter_per_plate'    , funnel%slope, domain_size(2)/3.0_rk)     
  call read_param_mpi(FILE, 'funnel', 'plates_thickness'      , funnel%plates_thickness, domain_size(1)/100.0_rk)     
  call read_param_mpi(FILE, 'funnel', 'pump_diameter'         , funnel%pump_diameter, domain_size(1)/5.0_rk)     
  call read_param_mpi(FILE, 'funnel', 'jet_diameter'          , funnel%jet_radius, domain_size(2)/20.0_rk)     
  call read_param_mpi(FILE, 'funnel', 'pump_x_center'         , funnel%pump_x_center, domain_size(1)*0.5_rk)     
  funnel%length               = domain_size(1)*0.85_rk
  funnel%plates_distance      = (funnel%length-funnel%plates_thickness)/(funnel%nr_plates-1)
  funnel%slope                = funnel%slope/(funnel%plates_distance)
  funnel%wall_thickness       = 0.02*domain_size(1)
  funnel%jet_radius           = funnel%jet_radius/2.0_rk
  funnel%temperatur           = 200.0_rk
  write(*,*) "s=", funnel%plates_distance
  write(*,*) "ds=", funnel%slope

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
    real(kind=rk)                                    :: p,rho,u,v,chi
    ! loop variables
    integer(kind=ik)                                 :: ix, iy,n
    ! outlets and inlets 
    real(kind=rk)                   :: velocity_pump,rho_pump, rho_capillary,u_capillary,v_capillary,p_capillary,p_2nd_pump_stage

!---------------------------------------------------------------------------------------------
! variables initialization
    if (size(penalization,1) /= Bs+2*g) call abort(777109,"wrong array size, there's pirates, captain!")

    ! reset mask array
    mask = 0.0_rk
    phi_ref=0.0_rk
!---------------------------------------------------------------------------------------------
! main body
    
    ! test!
    !> \todo put into infile
    velocity_pump   = 50.0_rk
    rho_capillary   = 1.645_rk
    rho_pump        = 1.645_rk
    u_capillary     = 100.0_rk
    v_capillary     = 0.0_rk
    p_capillary     = 101330.0_rk
    p_2nd_pump_stage= 70000.0_rk

    ! parameter for smoothing function (width)
    h = 1.5_rk*max(dx(1), dx(2))

    do iy=1, Bs+2*g
       y = dble(iy-(g+1)) * dx(2) + x0(2)
       do ix=1, Bs+2*g
            x = dble(ix-(g+1)) * dx(1) + x0(1)

            rho         = phi(ix,iy,1)**2
            u           = phi(ix,iy,2)/phi(ix,iy,1)
            v           = phi(ix,iy,3)/phi(ix,iy,1)
            p           = phi(ix,iy,4)
            

            ! Funnel
            ! ------
            ! 1. compute mask term:
            chi = draw_funnel_plates(x,y,funnel,h)
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
            chi = draw_walls(x,y,funnel,h)
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
            chi=  draw_pumps_volume_flow(x,y,funnel,h)
            if (chi>0) then                                 
              mask(ix,iy,2:3) = mask(ix,iy,2:3)+chi
              if (y>domain_size(2)/2) then
                Phi_ref(ix,iy,3) = rho*velocity_pump                                       
              else
                Phi_ref(ix,iy,3) = -rho*velocity_pump                                          
              endif                   
            endif 

            ! mass and energy sink
            chi=  draw_pumps_sink(x,y,funnel,h)
            if (chi>0) then                                 
              mask(ix,iy,1) = mask(ix,iy,1)+chi
              mask(ix,iy,4) = mask(ix,iy,4)+chi
              Phi_ref(ix,iy,1) = rho_pump                                       
              Phi_ref(ix,iy,4) = rho_pump*Rs*funnel%temperatur                                              
            endif                                             



            ! Inlet flow: Capillary
            ! ---------------------
            chi=  draw_jet(x,y,funnel,h)
            if (chi>0) then                                 
              mask(ix,iy,1:4) = mask(ix,iy,1:4)+chi
              Phi_ref(ix,iy,1) = rho_capillary                       
              Phi_ref(ix,iy,2) = rho_capillary*u_capillary  
              Phi_ref(ix,iy,3) = rho_capillary*v_capillary                              
              Phi_ref(ix,iy,4) = p_capillary   
            endif                                              


            ! ! Outlet flow: Transition to 2pump
            ! ! ---------------------
            ! chi=  draw_outlet(x,y,funnel,h)
               chi=  draw_sink(x,y,funnel,h)
              
              if (chi>0) then                                 
                mask(ix,iy,1) = mask(ix,iy,1)+chi
                mask(ix,iy,4) = mask(ix,iy,4)+chi
                Phi_ref(ix,iy,1) = rho_capillary                       
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




! !==========================================================================
! subroutine compute_mask_ref_value(x, y,phi, mask , phi_ref, h  )
!     implicit none
!     ! grid
!     integer(kind=ik), intent(in)                              :: Bs, g
!     !> mask term for every grid point of this block
!     real(kind=rk), dimension(4), intent(out)                  :: mask,phi_ref
!     !> spacing and origin of block
!     real(kind=rk), dimension(2), intent(in)                   :: x,y,phi

!     ! auxiliary variables
!     real(kind=rk)                                             :: x, y, r, h
!     ! loop variables
!     integer(kind=ik)                                          :: ix, iy
! !-------------------------------------------------------------------------
! ! variables initialization
!     ! Funnel
!             ! ------
!             ! 1. compute mask term:
!             chi = draw_funnel_plates(x,y,funnel,h)
!             ! 2. set quantities to desired values:
!             if (chi>0) then
!               mask(2:4) = mask(2:4) + chi
!               Phi_ref(2) = 0.0_rk                     ! no velocity in x
!               Phi_ref(3) = 0.0_rk                     ! no velocity in y
!               Phi_ref(4) = rho*Rs*funnel%temperatur   ! pressure set according to
!             endif                                           ! the temperature of the funnel   

!             ! Walls
!             ! -----
!             ! draw the walls arround the funnel
!             chi = draw_walls(x,y,funnel,h)
!             ! set reference values at these region
!             if (chi>0) then                       ! default values on the walls
!               mask(2:4) = mask(2:4) + chi
!               Phi_ref(2) = 0.0_rk                     ! no velocity in x
!               Phi_ref(3) = 0.0_rk                     ! no velocity in y
!               Phi_ref(4) = rho*Rs*funnel%temperatur   ! pressure set according to
!             endif                                           ! the temperature of the funnel   

!             ! Outlet flow: PUMPS
!             ! ------------------
!             ! pump volume flow
!             chi=  draw_pumps_volume_flow(x,y,funnel,h)
!             if (chi>0) then                                 
!               mask(2:3) = mask(2:3)+chi
!               if (y>domain_size(2)/2) then
!                 Phi_ref(3) = rho*velocity_pump                                       
!               else
!                 Phi_ref(3) = -rho*velocity_pump                                          
!               endif                   
!             endif 

!             ! mass and energy sink
!             chi=  draw_pumps_sink(x,y,funnel,h)
!             if (chi>0) then                                 
!               mask(1) = mask(1)+chi
!               mask(4) = mask(4)+chi
!               Phi_ref(1) = rho_pump                                       
!               Phi_ref(4) = rho_pump*Rs*funnel%temperatur                                              
!             endif                                             



!             ! Inlet flow: Capillary
!             ! ---------------------
!             chi=  draw_jet(x,y,funnel,h)
!             if (chi>0) then                                 
!               mask(1:4) = mask(1:4)+chi
!               Phi_ref(1) = rho_capillary                       
!               Phi_ref(2) = rho_capillary*u_capillary  
!               Phi_ref(3) = rho_capillary*v_capillary                              
!               Phi_ref(4) = p_capillary   
!             endif                                              


!             ! ! Outlet flow: Transition to 2pump
!             ! ! ---------------------
!             ! chi=  draw_outlet(x,y,funnel,h)
!                chi=  draw_sink(x,y,funnel,h)
              
!               if (chi>0) then                                 
!                 mask(1) = mask(1)+chi
!                 mask(4) = mask(4)+chi
!                 Phi_ref(1) = rho_capillary                       
!                 !Phi_ref(2) = 0
!                 Phi_ref(3) = 0                                      
!                 Phi_ref(4) = p_2nd_pump_stage
!               endif    
          
! end subroutine draw_funnel
! !==========================================================================






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
       do ix=1, Bs+2*g
           x = dble(ix-(g+1)) * dx(1) + x0(1) 
            mask(ix,iy) =   draw_walls(x,y,funnel,h) + &
                            draw_funnel_plates(x,y,funnel,h)
       end do
    end do
end subroutine draw_funnel
!==========================================================================





!> \brief Compute mask term for stacked rings of the ion funnel
function draw_funnel_plates(x,y,funnel,h)
 
  real(kind=rk)     , intent(in)  :: x, y, h
  type(type_funnel) , intent(in)  :: funnel
  real(kind=rk)                   :: draw_funnel_plates
  integer(kind=ik)                :: n
  ! distance frommin_inner_diameterter of cylinder
  draw_funnel_plates=0.0_rk
  
  ! loop over all plates
  do n=1,funnel%nr_plates

    if (x>(funnel%plate(n)%x0(1)-3*h) .and. &
        x<(funnel%plate(n)%x0(1)+funnel%plates_thickness+3*h)) then
                    draw_funnel_plates=draw_funnel_plates+draw_plate(x,y,funnel%plate(n),h)
    endif
  
  enddo

end function draw_funnel_plates


!> \brief Compute mask term for single ring-plate
function draw_plate(x,y,plate,h)
 
  real(kind=rk), intent(in)          :: x, y, h
  type(type_funnel_plate),intent(in) ::plate

  real(kind=rk)                      ::draw_plate,delta_r

  delta_r     = plate%r_out-plate%r_in
  draw_plate  = soft_bump(x,plate%x0(1),plate%width,h)*(soft_bump(y,plate%x0(2)+plate%r_in,delta_r,h) &
              + soft_bump(y,plate%x0(2)-(plate%r_out),delta_r,h))

end function draw_plate


function draw_walls(x,y,funnel,h)
 
  real(kind=rk),    intent(in)          :: x, y, h
  type(type_funnel),intent(in)          ::funnel

  real(kind=rk)                         ::  mask, draw_walls

  mask=0.0_rk

  if (abs(x-funnel%pump_x_center)> funnel%pump_diameter/2) then
         ! wall in south
         mask=mask+smoothstep(y-funnel%wall_thickness,h)
         ! wall in north
         mask=mask+smoothstep(domain_size(2)-y-funnel%wall_thickness,h)
  endif

  if (abs(y-domain_size(2)/2)> funnel%jet_radius) then
         ! wall in EAST
         mask=mask+smoothstep(x-funnel%wall_thickness,h)
   else
         mask=mask+smoothstep(x-funnel%wall_thickness/2,h)
  endif

  if (abs(y-domain_size(2)/2)> funnel%min_inner_diameter/2) then
         ! wall in WEST
         mask=mask+smoothstep(domain_size(1)-x-funnel%wall_thickness,h)
  else
         mask=mask+smoothstep(domain_size(1)-x-funnel%wall_thickness/2,h)
  endif

   ! is needed because mask off walls overlap
  if (mask>1) then
         mask=1
  endif

  draw_walls=mask
end function draw_walls

function draw_pumps_volume_flow(x,y,funnel,h)
 
  real(kind=rk),    intent(in)          :: x, y, h
  type(type_funnel),intent(in)          ::funnel

  real(kind=rk)                         ::  mask, draw_pumps_volume_flow,x0,x1,width

  mask  =0
  x0    =funnel%wall_thickness*0.5_rk + h
  width =funnel%wall_thickness*0.5_rk -2.0_rk*h
  x1    =domain_size(2)-funnel%wall_thickness+h 
  if (abs(x-funnel%pump_x_center)<= funnel%pump_diameter/2) then
         ! wall in south
         !mask=smoothstep(y-funnel%wall_thickness,h)
         mask=soft_bump(y,x0,width,h)
         ! wall in north
         mask=mask+soft_bump(y,x1,width,h)
  endif

   ! is needed because mask off walls overlap
  if (mask>1) then
         mask=1
  endif

  draw_pumps_volume_flow=mask
end function draw_pumps_volume_flow

function draw_pumps_sink(x,y,funnel,h)
 
  real(kind=rk),    intent(in)          :: x, y, h
  type(type_funnel),intent(in)          ::funnel

  real(kind=rk)                         ::  mask, draw_pumps_sink

  mask=0

  if (abs(x-funnel%pump_x_center)<= funnel%pump_diameter/2) then
         ! wall in south
         mask=smoothstep(y-funnel%wall_thickness*0.5_rk,h)
         ! wall in north
         mask=mask+smoothstep(domain_size(2)-funnel%wall_thickness*0.5_rk-y,h)
  endif

   ! is needed because mask off walls overlap
  if (mask>1) then
         mask=1
  endif

  draw_pumps_sink=mask
end function draw_pumps_sink

function draw_jet(x,y,funnel,h)
 
  real(kind=rk),    intent(in)          :: x, y, h
  type(type_funnel),intent(in)          ::funnel

  real(kind=rk)                         ::draw_jet

  if (abs(y-domain_size(2)/2)< funnel%jet_radius) then
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

 if (abs(y-domain_size(2)/2)< funnel%min_inner_diameter/2) then
         ! wall in WEST
    draw_outlet=smoothstep(domain_size(1)-x-funnel%wall_thickness,h)
  else
    draw_outlet=0.0_rk
  endif

end function draw_outlet


function draw_sink(x,y,funnel,h)
 
  real(kind=rk),    intent(in)          :: x, y, h
  type(type_funnel),intent(in)          ::funnel

  real(kind=rk)                         ::draw_sink,r

  r=sqrt((x-domain_size(1)+funnel%wall_thickness*0.6_rk)**2+(y-domain_size(2)*0.5_rk)**2)
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
                       domain_size(2)/2/)
      if(funnel%offset(1)<funnel%wall_thickness) then
       call abort(13457,'Error [module_mask.f90]: your funnel is to long')
      endif

      ! initialicd all plates
      do n=1,funnel%nr_plates
        plate%x0(1)     = funnel%offset(1)+(n-1)*distance
        plate%x0(2)     = funnel%offset(2)
        plate%width     = width
        if (plate%x0(1)-funnel%offset(1)<=length_focus) then
           plate%r_in    = funnel%max_inner_diameter/2
        else
           plate%r_in = plate%r_in - funnel%slope*distance/2
        endif
        plate%r_out   =funnel%outer_diameter/2
        funnel%plate(n)=plate   
      enddo
      
    end subroutine init_plates
!==========================================================================

