!-------------------------------------------------------------------------------
! WRAPPER Motion protocol wrapper
!-------------------------------------------------------------------------------
subroutine FlappingMotionWrap ( time, Insect, wingID )
  use module_insects_integration_flusi_wabbit
  implicit none

  real(kind=rk),intent(in) :: time
  real(kind=rk),dimension(0:3)::ep
  type(diptera) :: Insect
  integer(kind=2), intent(in) :: wingID

  select case ( wingID )
  case (1) !("left")
      call FlappingMotion ( time, Insect, Insect%FlappingMotion_left, &
      Insect%phi_l, Insect%alpha_l, Insect%theta_l, Insect%phi_dt_l,&
      Insect%alpha_dt_l, Insect%theta_dt_l, Insect%kine_wing_l, wingID )

  case (2) !("right")
      call FlappingMotion ( time, Insect, Insect%FlappingMotion_right, &
      Insect%phi_r, Insect%alpha_r, Insect%theta_r, Insect%phi_dt_r, &
      Insect%alpha_dt_r, Insect%theta_dt_r, Insect%kine_wing_r, wingID )

  case (3) !("left2")
      call FlappingMotion ( time, Insect, Insect%FlappingMotion_left2, &
      Insect%phi_l2, Insect%alpha_l2, Insect%theta_l2, Insect%phi_dt_l2,&
      Insect%alpha_dt_l2, Insect%theta_dt_l2, Insect%kine_wing_l2, wingID )

  case (4) !("right2")
      call FlappingMotion ( time, Insect, Insect%FlappingMotion_right2, &
      Insect%phi_r2, Insect%alpha_r2, Insect%theta_r2, Insect%phi_dt_r2, &
      Insect%alpha_dt_r2, Insect%theta_dt_r2, Insect%kine_wing_r2, wingID )

  case default
    call abort(77744, "not a valid wing identifier")
  end select

end subroutine FlappingMotionWrap


!-------------------------------------------------------------------------------
! Flapping wing motion protocoll, different choices.
! Input:
!       time (self explanatory)
!       protocoll: string containing what motion you want to have (may be
!                  different for both wings)
! Output:
!       phi: flapping/positonal angle
!       alpha: feathering angle / angle of attack
!       theta: deviation angle / Out-of-stroke-plane
!       phi_dt: flapping time derivative
!       alpha_dt: feathering time derivative
!       theta_dt: deviatory time derivative
!       kine: a Fourier series object, as described in insects.f90
! The actual motion depends on the choices in the parameter file, namely
! Insect%WingMotion, and sub-parameters that may further precise a given motion
! protocoll. Note we allow all four wings to follow a different motion, but they all
! call this routine here.
!-------------------------------------------------------------------------------
subroutine FlappingMotion(time, Insect, protocoll, phi, alpha, theta, phi_dt, &
           alpha_dt, theta_dt, kine, wingID)
  implicit none

  real(kind=rk), intent(in) :: time
  type(diptera), intent(inout) :: Insect
  real(kind=rk), intent(out) :: phi, alpha, theta, phi_dt, alpha_dt, theta_dt
  character (len=*), intent(in) :: protocoll
  type(wingkinematics), intent(inout) :: kine
  ! wingID is used for kineloader (to figure out which columsn to use from the data array)
  integer(kind=2), intent(in) :: wingID

  real(kind=rk) :: phi_max,alpha_max, phase,f
  real(kind=rk) :: bi_alpha_flapper(1:29) ! For comparison with Sane&Dickinson
  real(kind=rk) :: ai_phi_flapper(1:31) ! For comparison with Sane&Dickinson
  real(kind=rk) :: tadv,t0adv ! For comparison with Dickinson, Ramamurti
  real(kind=rk) :: posi,elev,feth,posi_dt,elev_dt,feth_dt,angles ! Comp. w. Maeda
  real(kind=rk) :: dangle_posi,dangle_elev,dangle_feth ! Comp. w. Maeda
  real(kind=rk) :: dangle_posi_dt,dangle_elev_dt,dangle_feth_dt ! Comp. w. Maeda
  real(kind=rk) :: a_posi(1:4),b_posi(1:4),a_elev(1:4),b_elev(1:4),a_feth(1:4),b_feth(1:4)
  real(kind=rk) :: s,c,a0_phi
  real(kind=rk) :: phicdeg
  real(kind=rk) :: alphacdeg, ttau
  integer :: i,mpicode
  character(len=clong) :: dummy
  type(inifile) :: kinefile

  select case ( protocoll )
  case ("from_file")
    !---------------------------------------------------------------------------
    ! Load kinematics for one stroke from file. This can be applied for example
    ! in hovering. if the wing motion varies appreciably between strokes,
    ! the kinematic loader is the method of choice. The file is specified
    ! in the params file and stored in Insect%infile. An example file
    ! (kinematics_fourier_example.ini) is in the git-repository
    !---------------------------------------------------------------------------
    if (index( kine%infile,".ini")==0) then
      call abort(2030,"you're trying to load an old kinematics file,please convert it to &
      &a new *.ini file (see src/insects/kinematics_example.ini&
      & the insects-tools repository can help you do that!")
    endif

    !---------------------------------------------------------------------------
    ! this block is excecuted only once per simulation, on first call of this routine
    !---------------------------------------------------------------------------
    if (.not.kine%initialized) then
        if (root) then
          write(*,'(80("<"))')
          write(*,*) "Initializing wing kinematics!"
          write(*,*) "*.ini file is: "//trim(adjustl(kine%infile))
          write(*,'(80("<"))')
        endif
      ! parse ini file
      call read_ini_file_mpi(kinefile, kine%infile, .true.)

      ! how to interpret numbers: Fourier or Hermite?
      call read_param_mpi(kinefile,"kinematics","type",kine%infile_type,"none")
      ! what units are given, degree or radiant?
      call read_param_mpi(kinefile,"kinematics","units",kine%infile_units,"degree")
      ! what convention/definition does the data follow?
      call read_param_mpi(kinefile,"kinematics","convention",kine%infile_convention,"flusi")

      ! inform about your interpretation
      select case (kine%infile_type)
      case ("Fourier", "fourier", "FOURIER")
        if (root) write(*,*) "The input file is interpreted as FOURIER coefficients"

      case ("Hermite", "hermite", "HERMITE")
        if (root) write(*,*) "The input file is interpreted as HERMITE coefficients"

      case default
        call abort(77771, "kinematics file does not appear to be valid, set type=fourier or type=hermite")

      end select

      !------------------------------------------------------------------------------------
      ! If type is _FOURIER_, we read a0, ai, bi:
      ! The Fourier series evaluation in WABBIT/FLUSI is :
      ! Q = a0_Q / 2 + ( a1_Q*cos(1*2*pi*t) + b1_Q*sin(1*2*pi*t) )
      !              + ( a2_Q*cos(2*2*pi*t) + b2_Q*sin(2*2*pi*t) )
      !              + ....
      ! Note the unfortunate division of a0 by 2, which is an historic artifact.
      !------------------------------------------------------------------------------------
      ! If the type is _HERMITE_, we read function values (ai) and derivatives (bi), 
      ! assuming implicitly an equidistant time vector [0, 1), thus excluding t=1.0
      !------------------------------------------------------------------------------------
      
      call read_param_mpi(kinefile,"kinematics","nfft_phi",kine%nfft_phi,0)
      call read_param_mpi(kinefile,"kinematics","a0_phi",kine%a0_phi,0.0_rk)
      call read_param_mpi(kinefile,"kinematics","ai_phi",kine%ai_phi(1:kine%nfft_phi))
      call read_param_mpi(kinefile,"kinematics","bi_phi",kine%bi_phi(1:kine%nfft_phi))
      
      call read_param_mpi(kinefile,"kinematics","nfft_alpha",kine%nfft_alpha,0)
      call read_param_mpi(kinefile,"kinematics","a0_alpha",kine%a0_alpha,0.0_rk)
      call read_param_mpi(kinefile,"kinematics","ai_alpha",kine%ai_alpha(1:kine%nfft_alpha))
      call read_param_mpi(kinefile,"kinematics","bi_alpha",kine%bi_alpha(1:kine%nfft_alpha))
      
      call read_param_mpi(kinefile,"kinematics","nfft_theta",kine%nfft_theta,0)
      call read_param_mpi(kinefile,"kinematics","a0_theta",kine%a0_theta,0.0_rk)
      call read_param_mpi(kinefile,"kinematics","ai_theta",kine%ai_theta(1:kine%nfft_theta))
      call read_param_mpi(kinefile,"kinematics","bi_theta",kine%bi_theta(1:kine%nfft_theta))

      kine%initialized = .true.
      call clean_ini_file_mpi( kinefile )

      if (root) write(*,'(80(">"))')
    endif

    !---------------------------------------------------------------------------
    ! get actual kinematics from the coefficients
    !---------------------------------------------------------------------------
    ! this block is executed every time. it is called a few times only per time step
    ! so don't worry about performance, we can use the string comparison
    select case (kine%infile_type)
    case ("Fourier","fourier","FOURIER")
      ! evaluate fourier series
      call fseries_eval(time,phi,phi_dt, kine%a0_phi, kine%ai_phi(1:kine%nfft_phi), kine%bi_phi(1:kine%nfft_phi))
      call fseries_eval(time,alpha,alpha_dt, kine%a0_alpha, kine%ai_alpha(1:kine%nfft_alpha), kine%bi_alpha(1:kine%nfft_alpha))
      call fseries_eval(time,theta,theta_dt, kine%a0_theta, kine%ai_theta(1:kine%nfft_theta), kine%bi_theta(1:kine%nfft_theta))

    case ("Hermite","hermite","HERMITE")
      ! evaluate hermite interpolation
      call hermite_eval(time, phi, phi_dt    , kine%ai_phi(1:kine%nfft_phi), kine%bi_phi(1:kine%nfft_phi))
      call hermite_eval(time, alpha, alpha_dt, kine%ai_alpha(1:kine%nfft_alpha), kine%bi_alpha(1:kine%nfft_alpha))
      call hermite_eval(time, theta, theta_dt, kine%ai_theta(1:kine%nfft_theta), kine%bi_theta(1:kine%nfft_theta))

    case default
      call abort(1717,"kinematics file does not appear to be valid, set type=fourier or type=hermite")
    end select

    !---------------------------------------------------------------------------
    ! make sure the output is in the right units (it HAS to be radiants!)
    !---------------------------------------------------------------------------
    select case (kine%infile_units)
    case ("degree","DEGREE","Degree","DEG","deg")
      ! the rest of the code gets radiants, so convert here
      phi    = deg2rad(phi)
      alpha  = deg2rad(alpha)
      theta  = deg2rad(theta)
      phi_dt = deg2rad(phi_dt)
      alpha_dt = deg2rad(alpha_dt)
      theta_dt = deg2rad(theta_dt)
    case ("radian","RADIAN","Radian","radiant","RADIANT","Radiant","rad","RAD")
      ! if the file is already in radiants, do nothing and be happy!
    case default
      call abort(1718,"kinematics file does not appear to be valid, set units=degree or units=radiant")
    end select

    !---------------------------------------------------------------------------
    ! make sure convention / definition of angles is respected
    !---------------------------------------------------------------------------
    select case (kine%infile_convention)
    case ("FLUSI","flusi","Flusi","FluSI")
      ! defintion as used in flusi. all angles are positive in the right hand rule
      ! that means especially that positive deviation puts wing downwards in a
      ! horizontal stroke plane
    case ("SISC","sisc","siam-sisc")
      ! the same as FLUSI, except for the deviation angle. The sign has been changed
      ! to be more in agreement with other people's definitions. however, flusi
      ! always internally works with the right hand rule
      theta    = -theta
      theta_dt = -theta_dt
    case default
      call abort(1719,"kinematics file does not appear to be valid, set convention=flusi or convention=sisc")
    end select


  case ("kinematics_loader")
    !--------------------------------------------------
    ! kinematics loader for non-periodic kinematics
    !--------------------------------------------------
    if (.not. Insect%kineloader_initialized) then
      call load_kine_init(Insect)
    endif

    ! fetch current wingkinematics from the data file, using hermite interpolation.
    ! NOTE: unlike earlier implementations, we assume here the file contains RADIANT and is in FLUSI 
    ! convention, and that time and velocities are in appropriate units (i.e. normalized externally
    ! during the generation of the data file)
    call wing_kine_interp(time, Insect, wingID, phi, alpha, theta, phi_dt, alpha_dt, theta_dt)

  case ("revolving-set1")
    ! revolving wing kinematics, pre-defined set. We fix alpha to 45deg and increase
    ! phi linearily with a short startup conditioner as suggested in [1]. The startup
    ! time is fixed to 0.4, which gives phi=31.35deg at the end of that interval
    ! [1] D. Kolomenskiy, Y. Elimelech and K. Schneider. Leading-edge vortex shedding from rotating wings. Fluid Dyn. Res., 46, 031421, 2014.
    ttau = 0.4
    ! position angle (is directly given in radian)
    ! we use PHI_DOT = 1 as normalization as well (since we have no frequency in this case)
    phi = 1._rk*( ttau*dexp(-time/ttau) + time)
    phi_dt = 1._rk*(1._rk-dexp(-time/ttau))
    ! feathering angle is constant
    alpha = deg2rad(-45._rk)
    alpha_dt = 0.0_rk
    ! elevation angle is always zero
    theta = 0.0_rk
    theta_dt = 0.0_rk

  case ("revolving-anticlock")
    ! revolving wing kinematics. Similar to "revolving-set1", but phi(0)=0
    ttau = 0.4
    ! position angle (is directly given in radian)
    ! we use PHI_DOT = 2*pi as normalization
    phi = 2.0_rk*pi*( ttau*dexp(-time/ttau) - ttau + time)
    phi_dt = 2.0_rk*pi*(1._rk-dexp(-time/ttau))
    ! feathering angle is constant
    alpha = deg2rad(-45._rk)
    alpha_dt = 0.0_rk
    ! elevation angle is always zero
    theta = 0.0_rk
    theta_dt = 0.0_rk

  case ("revolving-clock")
    ! revolving wing kinematics. Opposite direction to "revolving-anticlock"
    ttau = 0.4
    ! position angle (is directly given in radian)
    ! we use PHI_DOT = 2*pi as normalization
    phi = -2.0_rk*pi*( ttau*dexp(-time/ttau) - ttau + time)
    phi_dt = -2.0_rk*pi*(1._rk-dexp(-time/ttau))
    ! feathering angle is constant
    alpha = deg2rad(45._rk)
    alpha_dt = 0.0_rk
    ! elevation angle is always zero
    theta = 0.0_rk
    theta_dt = 0.0_rk

  case ("revolving-set2")
    ! revolving wing kinematics, pre-defined set. We fix alpha and increase
    ! phi linearily with a quadratic startup transient until phi=pi/8,
    ! which is reached at time equal to pi/4
    ! we use PHI_DOT = 1 as normalization
    if (time < pi/4.0_rk) then
      phi = 2.0_rk*time**2/pi
      phi_dt = 4.0_rk*time/pi
    else
      phi = time - pi/8.0_rk
      phi_dt = 1.0_rk
    endif
    ! feathering angle is constant
    alpha = - Insect%init_alpha_phi_theta(1) ! Mind the "-" sign
    alpha_dt = 0.0_rk
    ! elevation angle is always zero
    theta = 0.0_rk
    theta_dt = 0.0_rk

  case ("revolving-set3")
    ! revolving wing kinematics, pre-defined set. We fix alpha and increase
    ! phi linearily with a quadratic startup transient until phi=pi/16,
    ! which is reached at time equal to pi/8
    ! we use PHI_DOT = 1 as normalization
    ! Therefore, the acceleration is twice as fast as in "revolving-set2"
    if (time < pi/8.0_rk) then
      phi = 4.0_rk*time**2/pi
      phi_dt = 8.0_rk*time/pi
    else
      phi = time - pi/16.0_rk
      phi_dt = 1.0_rk
    endif
    ! feathering angle is constant
    alpha = - Insect%init_alpha_phi_theta(1) ! Mind the "-" sign
    alpha_dt = 0.0_rk
    ! elevation angle is always zero
    theta = 0.0_rk
    theta_dt = 0.0_rk

case ("revolving-set3-degree")
    ! revolving wing kinematics, pre-defined set. We fix alpha and increase
    ! phi linearily with a quadratic startup transient until phi=pi/16,
    ! which is reached at time equal to pi/8
    ! we use PHI_DOT = 1 as normalization
    ! Therefore, the acceleration is twice as fast as in "revolving-set2"
    if (time < pi/8.0_rk) then
      phi = 4.0_rk*time**2/pi
      phi_dt = 8.0_rk*time/pi
    else
      phi = time - pi/16.0_rk
      phi_dt = 1.0_rk
    endif
    ! feathering angle is constant
    alpha = - deg2rad(Insect%init_alpha_phi_theta(1)) ! Mind the "-" sign
    alpha_dt = 0.0_rk
    ! elevation angle is always zero
    theta = 0.0_rk
    theta_dt = 0.0_rk

  case ("Drosophila_hovering_fry")
    !---------------------------------------------------------------------------
    ! motion protocoll digitalized from Fry et al JEB 208, 2303-2318 (2005)
    !
    ! fourier coefficients analyzed with matlab
    !---------------------------------------------------------------------------
    if (.not.kine%initialized) then
      kine%nfft_alpha = 10
      kine%nfft_theta = 10
      kine%nfft_phi   = 10

      kine%a0_phi   =25.4649398
      kine%a0_alpha =-0.3056968
      kine%a0_theta =-17.8244658  ! - sign (Dmitry, 10 Nov 2013)

      kine%ai_phi(1:kine%nfft_phi) = (/71.1061858,2.1685448,-0.1986978,0.6095268,-0.0311298,&
                                       -0.1255648,-0.0867778,0.0543518,0.0,0.0/)
      kine%bi_phi(1:kine%nfft_phi) = (/5.4547058,-3.5461688,0.6260698,0.1573728,-0.0360498,-0.0205348,&
                                      -0.0083818,-0.0076848,0.0,0.0/)
      kine%ai_alpha(1:kine%nfft_alpha) = (/3.3288788,0.6303878,-10.9780518,2.1123398,-3.2301198,&
                                          -1.4473158,0.6141758,-0.3071608,0.1458498,0.0848308/)
      kine%bi_alpha(1:kine%nfft_alpha) = (/67.5430838,0.6566888,9.9226018,3.9183988,-2.6882828,0.6433518,&
                                          -0.8792398,-0.4817838,0.0300078,-0.1015118/)
      kine%ai_theta(1:kine%nfft_theta) = (/-3.9750378,-8.2808998,0.0611208,0.3906598,-0.4488778,0.120087,&
                                          0.0717048,-0.0699578,0.0,0.0/)   ! - sign (Dmitry, 10 Nov 2013)
      kine%bi_theta(1:kine%nfft_theta) = (/-2.2839398,-3.5213068,1.9296668,-1.0832488,-0.3011748,0.1786648,&
                                          -0.1228608,0.0004808,0.0,0.0/)   ! - sign (Dmitry, 10 Nov 2013)
    endif
    kine%initialized = .true.

    call fseries_eval(time,phi,phi_dt    ,kine%a0_phi, kine%ai_phi(1:kine%nfft_phi), kine%bi_phi(1:kine%nfft_phi))
    call fseries_eval(time,alpha,alpha_dt,kine%a0_alpha, kine%ai_alpha(1:kine%nfft_alpha), kine%bi_alpha(1:kine%nfft_alpha))
    call fseries_eval(time,theta,theta_dt,kine%a0_theta, kine%ai_theta(1:kine%nfft_theta), kine%bi_theta(1:kine%nfft_theta))

    phi =  deg2rad(phi)
    alpha = deg2rad(alpha)
    theta = deg2rad(theta)

    phi_dt = deg2rad(phi_dt)
    alpha_dt = deg2rad(alpha_dt)
    theta_dt = deg2rad(theta_dt)

  case ("Drosophila_hovering_sun")
    !---------------------------------------------------------------------------
    ! 6th wingbeat from Chen & Sun (2014)
    ! Fourier coefficients analyzed with matlab
    ! Note that it begins with UPSTROKE, unlike Fry's kinematics
    !---------------------------------------------------------------------------
    if (.not.kine%initialized) then
      kine%nfft_alpha = 10
      kine%nfft_theta = 10
      kine%nfft_phi   = 10

      kine%a0_phi   =38.2280144124915
      kine%a0_alpha =1.09156750841542
      kine%a0_theta =-17.0396438317138
      kine%ai_phi(1:kine%nfft_phi)   =(/-60.4766838884452,-5.34194400577534,-2.79100982711466,&
      0.334561891664093,-0.0202655263017256,-0.323801616724358,&
      -0.474435962657283,-0.111655938200879,0.00958151551130752,&
      0.119666224142596/)
      kine%bi_phi(1:kine%nfft_phi)   =(/3.00359559936894,-7.55184535084148,-1.32520563461209,&
      -0.297351445375239,-0.213108013812305,-0.0328282543472566,&
      0.0146299151981855,-0.0385423658663155,-0.512411386850196,&
      0.0785978606299901/)
      kine%ai_alpha(1:kine%nfft_alpha) =(/-7.73228268249956,8.09409174482393,3.98349294858406,&
      6.54460609657175,4.20944598804824,0.138380341939039,&
      1.38813149742271,0.625930107014395,0.607953761451392,&
      -0.688049862096416/)
      kine%bi_alpha(1:kine%nfft_alpha) =(/-52.2064112743351,-3.83568699799253,-14.8200023306913,&
      4.57428035662431,1.01656074807431,-0.113387395332322,&
      0.614350733080735,0.637197524906845,0.878128257565861,&
      0.271333646229075/)
      kine%ai_theta(1:kine%nfft_theta)  =(/-0.0876540586926952,-6.32825811106895,0.461710021840119,&
      -0.196290365124456,0.266829219535534,0.191278837298358,&
      0.0651359677824226,0.132751714936873,0.104342707251547,&
      0.0251049936194829/)
      kine%bi_theta(1:kine%nfft_theta)  =(/5.05944308805575,-0.677547202362310,-1.30945644385840,&
      -0.741962147111828,-0.351209215835472,-0.0374224119537382,&
      -0.0940160961803885,-0.0563224030429001,0.0533369476976694,&
      0.0507212428142968/)
      kine%initialized = .true.
    endif

    call fseries_eval(time,phi,phi_dt    ,kine%a0_phi, kine%ai_phi(1:kine%nfft_phi), kine%bi_phi(1:kine%nfft_phi))
    call fseries_eval(time,alpha,alpha_dt,kine%a0_alpha, kine%ai_alpha(1:kine%nfft_alpha), kine%bi_alpha(1:kine%nfft_alpha))
    call fseries_eval(time,theta,theta_dt,kine%a0_theta, kine%ai_theta(1:kine%nfft_theta), kine%bi_theta(1:kine%nfft_theta))

    phi =  deg2rad(phi)
    alpha = deg2rad(alpha)
    theta = deg2rad(theta)

    phi_dt = deg2rad(phi_dt)
    alpha_dt = deg2rad(alpha_dt)
    theta_dt = deg2rad(theta_dt)
  case ("Drosophila_hovering_maeda")
    !---------------------------------------------------------------------------
    ! Drosophila hovering kinematics protocol
    !
    ! Fourier coefficients provided by Maeda
    ! Diditized from Fry et al.
    !---------------------------------------------------------------------------
    a_posi = (/  0.22700_rk,  1.24020_rk,  0.03610_rk, -0.00360_rk/)
    b_posi = (/  0.00000_rk,  0.08880_rk, -0.07000_rk,  0.01250_rk/)
    a_elev = (/  0.16125_rk,  0.06750_rk,  0.14500_rk,  0.00540_rk/)
    b_elev = (/  0.00000_rk,  0.03670_rk,  0.06840_rk, -0.03390_rk/)
    a_feth = (/ -0.00864_rk, -0.04890_rk, -0.02056_rk,  0.19649_rk/)
    b_feth = (/  0.00000_rk, -1.17586_rk, -0.01216_rk, -0.17590_rk/)

    ! Initialize angles and velocities
    posi = 0.0_rk
    elev = 0.0_rk
    feth = 0.0_rk
    posi_dt = 0.0_rk
    elev_dt = 0.0_rk
    feth_dt = 0.0_rk

    do i=0,3 !! Fourier series
      !! time dependent angle
      angles  = 2.0_rk*dble(i)*pi*time

      selectcase( i )
      case( 0 ) ! Fourier 0th order
        ! mean
        dangle_posi = a_posi(1) ! +shift_mean_posi_
        dangle_elev = a_elev(1) ! +shift_mean_elev_
        dangle_feth = a_feth(1) ! +shift_mean_feth_

        dangle_posi_dt = 0.0_rk
        dangle_elev_dt = 0.0_rk
        dangle_feth_dt = 0.0_rk

      case default !! Fourier n-th orders

        call get_dangle( &
        & angles, &                !! intent(in)
        & i, &                     !! intent(in)
        & a_posi(i+1), & !! intent(in)
        & b_posi(i+1), & !! intent(in)
        & 0.0_rk, &                 !! intent(in)
        & 0.0_rk, &                 !! intent(in)
        & dangle_posi, &      !! intent(out)
        & dangle_posi_dt &   !! intent(out)
        & )

        call get_dangle( &
        & angles, &                !! intent(in
        & i, &                     !! intent(in)
        & a_elev(i+1), & !! intent(in)
        & b_elev(i+1), & !! intent(in)
        & 0.0_rk, &                 !! intent(in)
        & 0.0_rk, &                 !! intent(in)
        & dangle_elev, &      !! intent(out)
        & dangle_elev_dt &   !! intent(out)
        & )

        call get_dangle( &
        & angles, &                !! intent(in
        & i, &                     !! intent(in)
        & a_feth(i+1), & !! intent(in)
        & b_feth(i+1), & !! intent(in)
        & 0.0_rk, &                 !! intent(in)
        & 0.0_rk, &                 !! intent(in)
        & dangle_feth, &      !! intent(out)
        & dangle_feth_dt &   !! intent(out)
        & )

      endselect

      posi = posi +dangle_posi
      elev = elev +dangle_elev
      feth = feth +dangle_feth

      posi_dt = posi_dt +dangle_posi_dt
      elev_dt = elev_dt +dangle_elev_dt
      feth_dt = feth_dt +dangle_feth_dt
    enddo

    ! Convert to FLUSI's variables
    phi = posi
    alpha = -feth
    theta = -elev

    phi_dt = posi_dt
    alpha_dt = -feth_dt
    theta_dt = -elev_dt

  case ("flapper_sane")
    !---------------------------------------------------------------------------
    ! motion protocol from Sane and Dickinson, JEB 204, 2607-2626 (2001)
    !
    ! feathering: fourier coefficients analyzed with matlab, 2nd order
    !             Butterworth filter with cutoff at k=10
    ! positional: similar to above
    ! elevation:  zero
    !
    ! Dmitry, 2 Nov 2013
    !---------------------------------------------------------------------------

    ! *** I. feathering motion ***
    ! Corresponds to Fig. 3D in JEB 204, p. 2613
    ! Note that this is feathering angle measured from the vertical.
    ! This is NOT angle of attack
    bi_alpha_flapper =(/48.807554373967804_rk,&
    0.0_rk,11.14661083909663_rk,0.0_rk,2.242734216805251_rk,&
    0.0_rk,-0.6141899985692184_rk,0.0_rk,-0.7426551158681146_rk,&
    0.0_rk,-0.2329560587573768_rk,0.0_rk,0.038749678276091284_rk,&
    0.0_rk,0.07083462320831221_rk,0.0_rk,0.028982501947490313_rk,&
    0.0_rk,-0.0025202918494477244_rk,0.0_rk,-0.010221019942802941_rk,&
    0.0_rk,-0.005614021318470698_rk,0.0_rk,1.1958884364596903d-6,&
    0.0_rk,0.002186832241254999_rk,0.0_rk,0.0015347995090793172_rk/)

    alpha = 0.0
    alpha_dt = 0.0

    ! frequency factor
    f = 2.0_rk*pi

    ! Fourier series
    do i=1,29
      ! allows the spaces I like with the 80 columns malcolm likes :)
      s = dsin(f*dble(i)*time)
      c = dcos(f*dble(i)*time)
      alpha = alpha + bi_alpha_flapper(i) * s
      alpha_dt = alpha_dt + f*dble(i)* bi_alpha_flapper(i) * c
    enddo

    ! Scale to a given value of max angle in gedrees
    ! alphacdeg is 90deg MINUS alpha of JEB 204 (eg alphedeg=90-50 for Fig 3D)
    alphacdeg = 90.0_rk - 00.0_rk
    alpha = alphacdeg/40.0_rk * alpha
    alpha_dt = alphacdeg/40.0_rk * alpha_dt

    ! convert in radians
    alpha = deg2rad(alpha)
    alpha_dt = deg2rad(alpha_dt)

    ! *** II. position ***
    ai_phi_flapper =(/72.96795908179631_rk,&
    0.0_rk,8.064401876272864_rk,0.0_rk,2.769062401215844_rk,&
    0.0_rk,1.2200252377066352_rk,0.0_rk,0.5584689705779989_rk,&
    0.0_rk,0.2545617536476344_rk,0.0_rk,0.11829515180579572_rk,&
    0.0_rk,0.05754453975774996_rk,0.0_rk,0.02964141751269772_rk,&
    0.0_rk,0.016177705089515895_rk,0.0_rk,0.009315101869467001_rk,&
    0.0_rk,0.005625663922446026_rk,0.0_rk,0.0035424425357352385_rk,&
    0.0_rk,0.0023130422432356247_rk,0.0_rk,0.001558278163264511_rk,&
    0.0_rk,0.001078213692334021_rk/)

    phi = 0.0
    phi_dt = 0.0

    ! frequency factor
    f = 2.0_rk*pi

    ! Fourier series
    do i=1,31
      ! allows the spaces I like with the 80 columns malcolm likes :)
      s = dsin(f*dble(i)*time)
      c = dcos(f*dble(i)*time)
      phi   = phi   + ai_phi_flapper(i) * c
      phi_dt   = phi_dt   + f*dble(i)*(-ai_phi_flapper(i) * s)
    enddo

    ! Scale to a given value of max angle in gedrees
    ! phicdeg is Phi of JEB 204 (eg twice the max value of triangular wave)
    phicdeg = 180.0_rk
    phi = phicdeg/180.0_rk * phi
    phi_dt = phicdeg/180.0_rk * phi_dt

    ! convert in radians
    phi = deg2rad(phi)
    phi_dt = deg2rad(phi_dt)

    ! *** III. elevation ***
    theta = 0.0_rk
    theta_dt = 0.0_rk

  case ("flapper_dickinson")
    !---------------------------------------------------------------------------
    ! motion protocol from Dickinson, Lehmann and Sane, Science (1999)
    !
    ! feathering: fourier coefficients analyzed with matlab,
    !             Gaussian filter that fits fig 3D
    ! positional: similar to above
    ! elevation:  zero
    !
    ! Dmitry, 5 Nov 2013
    !---------------------------------------------------------------------------

    ! *** I. feathering motion ***
    ! Corresponds to Fig. 3D in Science
    ! Note that this is feathering angle measured from the vertical.
    ! This is NOT angle of attack
    bi_alpha_flapper =(/48.23094285611071_rk,&
    0.0_rk,10.224154661301371_rk,0.0_rk,2.1623763046726396_rk,&
    0.0_rk,0.05049394424178093_rk,0.0_rk,-0.17550942623071494_rk,&
    0.0_rk,-0.06634193748204852_rk,0.0_rk,-0.008925020495896451_rk,&
    0.0_rk,0.0011292567942149407_rk,0.0_rk,6.471071566666472d-4,&
    0.0_rk,1.0018757795834964d-4,0.0_rk,3.0105550216312524d-6,&
    0.0_rk,-1.237567150768195d-6,0.0_rk,-1.988004402010933d-7,&
    0.0_rk,-1.10165545174181d-8,0.0_rk,2.4135650975460306d-10/)

    ! Advanced rotation (+ sign) or delayed rotation (- sign)
    tadv = 0.08
    !tadv = - 0.08

    alpha = 0.0
    alpha_dt = 0.0

    ! frequency factor
    f = 2.0_rk*pi

    ! Fourier series
    do i=1,29
      ! allows the spaces I like with the 80 columns malcolm likes :)
      s = dsin(f*dble(i)*(time+tadv))
      c = dcos(f*dble(i)*(time+tadv))
      alpha = alpha + bi_alpha_flapper(i) * s
      alpha_dt = alpha_dt + f*dble(i)* bi_alpha_flapper(i) * c
    enddo

    ! Scale to a given value of max angle in gedrees
    ! alphacdeg is 90deg MINUS alpha of Science (eg alphedeg=90-40 for Fig 3)
    alphacdeg = 90.0_rk - 40.0_rk
    alpha = alphacdeg/40.0_rk * alpha
    alpha_dt = alphacdeg/40.0_rk * alpha_dt

    ! convert in radians
    alpha = deg2rad(alpha)
    alpha_dt = deg2rad(alpha_dt)

    ! *** II. position ***
    ai_phi_flapper =(/63.24528806534019_rk,&
    0.0_rk,5.753991800610726_rk,0.0_rk,1.3887974015525626_rk,&
    0.0_rk,0.3889856512386744_rk,0.0_rk,0.10577402496901325_rk,&
    0.0_rk,0.026061339604144987_rk,0.0_rk,0.005623376646981709_rk,&
    0.0_rk,0.001042285996467963_rk,0.0_rk,1.639611509380189d-4,&
    0.0_rk,2.1716252827442023d-5,0.0_rk,2.408190194815521d-6,&
    0.0_rk,2.2268710288534648d-7,0.0_rk,1.7118916093759426d-8,&
    0.0_rk,1.0914870312823793d-9,0.0_rk,5.76135101855556d-11,&
    0.0_rk,2.513944479978149d-12/)

    phi = 0.0
    phi_dt = 0.0

    ! frequency factor
    f = 2.0_rk*pi

    ! Fourier series
    do i=1,31
      ! allows the spaces I like with the 80 columns malcolm likes :)
      s = dsin(f*dble(i)*time)
      c = dcos(f*dble(i)*time)
      phi   = phi   + ai_phi_flapper(i) * c
      phi_dt   = phi_dt   + f*dble(i)*(-ai_phi_flapper(i) * s)
    enddo

    ! Scale to a given value of max angle in gedrees
    ! phicdeg is Phi of JEB 204 (eg twice the max value of triangular wave)
    phicdeg = 180.0_rk
    phi = phicdeg/180.0_rk * phi
    phi_dt = phicdeg/180.0_rk * phi_dt

    ! convert in radians
    phi = deg2rad(phi)
    phi_dt = deg2rad(phi_dt)

    ! *** III. elevation ***
    theta = 0.0_rk
    theta_dt = 0.0_rk

  case ("flapper_ramamurti")
    !---------------------------------------------------------------------------
    ! motion protocol from Ramamurti and Sandberg, JEB (2002)
    !
    ! feathering: fourier coefficients analyzed with matlab,
    !             Gaussian filter that fits fig 2
    ! positional: similar to above
    ! elevation:  zero
    !
    ! Dmitry, 10 Jun 2013
    !---------------------------------------------------------------------------

    ! Time shift as in JEB
    t0adv = 0.28175_rk

    ! *** I. feathering motion ***
    ! Note that this is feathering angle measured from the vertical.
    ! This is NOT angle of attack
    !    bi_alpha_flapper =(/48.23094285611071_rk,&
    !      0.0_rk,10.224154661301371_rk,0.0_rk,2.1623763046726396_rk,&
    !      0.0_rk,0.05049394424178093_rk,0.0_rk,-0.17550942623071494_rk,&
    !      0.0_rk,-0.06634193748204852_rk,0.0_rk,-0.008925020495896451_rk,&
    !      0.0_rk,0.0011292567942149407_rk,0.0_rk,6.471071566666472d-4,&
    !      0.0_rk,1.0018757795834964d-4,0.0_rk,3.0105550216312524d-6,&
    !      0.0_rk,-1.237567150768195d-6,0.0_rk,-1.988004402010933d-7,&
    !      0.0_rk,-1.10165545174181d-8,0.0_rk,2.4135650975460306d-10/)

    bi_alpha_flapper =(/58.5622945117485_rk,&
    0.0_rk,9.70856389020196_rk,0.0_rk,1.06772463979698_rk,&
    0.0_rk,-0.127455709998572_rk,0.0_rk,-0.0553559839123380_rk,&
    0.0_rk,-0.00500430568361157_rk,0.0_rk,0.000181637429802491_rk,&
    0.0_rk,5.24470001944736d-05,0.0_rk,2.39488192125763d-06,&
    0.0_rk,-1.62462115890220d-08,0.0_rk,-3.60732923972996d-09,&
    0.0_rk,-7.58943924434529d-11,0.0_rk,-3.17275659750908d-16,&
    0.0_rk,1.50381167957580d-14,0.0_rk,1.41146073296091d-16/)

    ! Advanced rotation (+ sign) or delayed rotation (- sign)
    tadv = 0.08_rk
    !tadv = - 0.08_rk

    alpha = 0.0_rk
    alpha_dt = 0.0_rk

    ! frequency factor
    f = 2.0_rk*pi

    ! Fourier series
    do i=1,29
      ! allows the spaces I like with the 80 columns malcolm likes :)
      s = dsin(f*dble(i)*(time+tadv+t0adv))
      c = dcos(f*dble(i)*(time+tadv+t0adv))
      alpha = alpha + bi_alpha_flapper(i) * s
      alpha_dt = alpha_dt + f*dble(i)* bi_alpha_flapper(i) * c
    enddo

    ! convert in radians
    alpha = deg2rad(alpha)
    alpha_dt = deg2rad(alpha_dt)

    ! *** II. position ***
    !    ai_phi_flapper =(/63.24528806534019_rk,&
    !      0.0_rk,5.753991800610726_rk,0.0_rk,1.3887974015525626_rk,&
    !      0.0_rk,0.3889856512386744_rk,0.0_rk,0.10577402496901325_rk,&
    !      0.0_rk,0.026061339604144987_rk,0.0_rk,0.005623376646981709_rk,&
    !      0.0_rk,0.001042285996467963_rk,0.0_rk,1.639611509380189d-4,&
    !      0.0_rk,2.1716252827442023d-5,0.0_rk,2.408190194815521d-6,&
    !      0.0_rk,2.2268710288534648d-7,0.0_rk,1.7118916093759426d-8,&
    !      0.0_rk,1.0914870312823793d-9,0.0_rk,5.76135101855556d-11,&
    !      0.0_rk,2.513944479978149d-12/)

    a0_phi = 20.0_rk
    ai_phi_flapper =(/71.1023524748246_rk,&
    0.0_rk,6.43355277659369_rk,0.0_rk,1.53592591540135_rk,&
    0.0_rk,0.423188982190934_rk,0.0_rk,0.112580047175950_rk,&
    0.0_rk,0.0269873997935609_rk,0.0_rk,0.00563410946130459_rk,&
    0.0_rk,0.00100469248408058_rk,0.0_rk,0.000151188714225520_rk,&
    0.0_rk,1.90437108218227d-05,0.0_rk,1.99627349027768d-06,&
    0.0_rk,1.73399283714958d-07,0.0_rk,1.24381425530736d-08,&
    0.0_rk,7.34711932879597d-10,0.0_rk,3.56493134534579d-11,&
    0.0_rk,1.41757814786810d-12/)

    phi = 0.5_rk * a0_phi
    phi_dt = 0.0_rk

    ! frequency factor
    f = 2.0_rk*pi

    ! Fourier series
    do i=1,31
      ! allows the spaces I like with the 80 columns malcolm likes :)
      s = dsin(f*dble(i)*(time+t0adv))
      c = dcos(f*dble(i)*(time+t0adv))
      phi   = phi   + ai_phi_flapper(i) * c
      phi_dt   = phi_dt   + f*dble(i)*(-ai_phi_flapper(i) * s)
    enddo

    ! convert in radians
    phi = deg2rad(phi)
    phi_dt = deg2rad(phi_dt)

    ! *** III. elevation ***
    theta = 0.0_rk
    theta_dt = 0.0_rk



  case ("suzuki")

    ! frequency
    f = 2.0_rk*pi

    ! amplitudes and other parameters
    phi_max = 80
    alpha_max = 45
    c = 3.3

    ! *** I. position ***
    phi = phi_max * cos(f*time)
    phi_dt = -phi_max*f * sin(f*time)
    phi = deg2rad(phi)
    phi_dt = deg2rad(phi_dt)

    ! *** II. feathering ***
    alpha = alpha_max/tanh(c)*tanh(c*sin(f*time))
    alpha_dt = alpha_max/tanh(c)*f*c*cos(f*time)*(1-tanh(c*sin(f*time))**2)
    alpha = deg2rad(alpha)
    alpha_dt = deg2rad(alpha_dt)

    ! *** III. elevation ***
    theta = 0.0_rk
    theta_dt = 0.0_rk

  case ("simplified")
    !---------------------------------------------------------------------------
    ! simplified motion protocoll
    !
    ! J. comput. Phys. 231 (2012) 1822-1847 "A fluid-structure interaction
    ! model of insect flight with flexible wings"
    !
    ! the pase shift "phase" was my idea
    !---------------------------------------------------------------------------
    phi_max     = deg2rad(80.0_rk)  ! phi is up/down angle (flapping)
    alpha_max   = deg2rad(45._rk)  ! alpha is tethering
    phase       = 0.0_rk! 10.0_rk*pi/180.0_rk  ! phase shift between flapping and tethering
    f = 1._rk*2.0*pi

    phi      = phi_max  *dcos(f*time)
    alpha    = alpha_max*dsin(f*(time+phase))
    theta    = 0.0

    phi_dt   =-phi_max *f *dsin(f*time)
    alpha_dt = alpha_max*f*dcos(f*(time+phase))
    theta_dt = 0.0

  case ("simplified2")
    !---------------------------------------------------------------------------
    ! simplified motion protocoll
    !---------------------------------------------------------------------------
    phi_max     = deg2rad(60.0_rk)  ! phi is up/down angle (flapping)
    alpha_max   = deg2rad(0.0_rk)  ! alpha is tethering
    phase       = 0.0_rk! 10.0_rk*pi/180.0_rk  ! phase shift between flapping and tethering
    f = 1._rk*2.0*pi

    phi      = phi_max  *dcos(f*time)
    alpha    = alpha_max*dsin(f*(time+phase))
    theta    = 0.0

    phi_dt   =-phi_max *f *dsin(f*time)
    alpha_dt = alpha_max*f*dcos(f*(time+phase))
    theta_dt = 0.0

  case ("debug")
    phi      = 0.0
    alpha    = deg2rad(-45._rk)
    theta    = 0.0
    phi_dt   = 0.0
    alpha_dt = 0.0
    theta_dt = 0.0
  case ("debug2")
    phi      = 0.0
    alpha    = deg2rad(+45._rk)
    theta    = 0.0
    phi_dt   = 0.0
    alpha_dt = 0.0
    theta_dt = 0.0
  case ("none")
    phi      = 0.0
    alpha    = deg2rad(45._rk)
    theta    = 0.0
    phi_dt   = 0.0
    alpha_dt = 0.0
    theta_dt = 0.0
  case ("command-line-left")
    ! use values specified in global variables. this is used to draw a single
    ! mask with a dry run, where all 12 parameters are specified in the command line
    if(root) write(*,*) "Reading left wing kinematics (params 10,11,12: phi,alpha,theta)"
    if(root) write(*,*) "note wings do NOT move (no velocity field)"
    call get_command_argument(10,dummy)
    read (dummy,*) phi
    call get_command_argument(11,dummy)
    read (dummy,*) alpha
    call get_command_argument(12,dummy)
    read (dummy,*) theta

    if(root) write(*,'("phi=",g12.4,"° alpha=",g12.4,"° theta=",g12.4,"° ")') phi,alpha,theta
    phi = deg2rad(phi)
    alpha = deg2rad(alpha)
    theta = deg2rad(theta)
    phi_dt   = 0.0
    alpha_dt = 0.0
    theta_dt = 0.0

  case ("command-line-right")
    ! use values specified in global variables. this is used to draw a single
    ! mask with a dry run, where all 12 parameters are specified in the command line
    if(root) write(*,*) "Reading left wing kinematics (params 13,14,15: phi,alpha,theta)"
    if(root) write(*,*) "note wings do NOT move (no velocity field)"
    call get_command_argument(13,dummy)
    read (dummy,*) phi
    call get_command_argument(14,dummy)
    read (dummy,*) alpha
    call get_command_argument(15,dummy)
    read (dummy,*) theta

    if(root) write(*,'("phi=",g12.4,"° alpha=",g12.4,"° theta=",g12.4,"° ")') phi,alpha,theta
    phi = deg2rad(phi)
    alpha = deg2rad(alpha)
    theta = deg2rad(theta)
    phi_dt   = 0.0
    alpha_dt = 0.0
    theta_dt = 0.0
  case default
    write(*,*) "value is: "//trim(adjustl(protocoll))
    call abort(121212,"insects.f90::FlappingMotion: motion case (protocoll) undefined")
  end select

end subroutine FlappingMotion
