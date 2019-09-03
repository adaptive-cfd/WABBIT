!-------------------------------------------------------------------------------
! Body motion protocoll, different choices.
! Input:
!      time (self explanatory)
! Output:
!      Insect% psi:     roll angle
!      Insect%beta:     pitch angle
!      Insect%gamma:    yaw angle
!      Insect%psi_dt:   roll angular velocity
!      Insect%beta_dt:  pitch angular velocity
!      Insect%gamma_dt: yaw angular velocity
!      Insect%xc:       center of gravity coordinate
!      Insect%vc:       translational velocity of the body
! The actual motion depends on the choices in the parameter file, namely
! Insect%BodyMotion, and sub-parameters that may further precise a given motion
! protocoll
! Note that in new versions, all the angles and positions are stored in one
! datastructure, which is then the only output variable of this routine.
!-------------------------------------------------------------------------------
subroutine BodyMotion(time, Insect)
    implicit none

    real(kind=rk), intent(in) :: time
    type(diptera), intent(inout) :: Insect
    real(kind=rk) :: psi, beta, gamma, psi_dt, beta_dt, gamma_dt
    real(kind=rk) :: xc(1:3), vc(1:3), ep(0:3)
    real(kind=rk) :: T,R
    character(len=strlen) :: dummy

    ! the tag body_moves is used to draw the insect's body only once, if the body
    ! does not move (body_moves=="no"). For safety, we initialize the body as moving
    ! so if you forget to specify (body_moves=="no"), the body is drawn every time
    Insect%body_moves = "yes"
    Insect%quaternion_solver_used = .false.

    select case (Insect%BodyMotion)
    case ("command-line")
        ! use values specified in global variables. this is used to draw a single
        ! mask with a dry run, where all 12 parameters are specified in the command line
        if(root) write(*,*) "Reading left wing kinematics (params 4,5,6:x,y,z and 7,8,9:psi,beta,gamma)"
        if(root) write(*,*) "note body does NOT move (no velocity field)"
        call get_command_argument(4,dummy)
        read (dummy,*) xc(1)
        call get_command_argument(5,dummy)
        read (dummy,*) xc(2)
        call get_command_argument(6,dummy)
        read (dummy,*) xc(3)
        call get_command_argument(7,dummy)
        read (dummy,*) psi
        call get_command_argument(8,dummy)
        read (dummy,*) beta
        call get_command_argument(9,dummy)
        read (dummy,*) gamma

        Insect%body_moves = "no"
        psi_dt   = 0.d0
        beta_dt  = 0.d0
        gamma_dt = 0.d0
        vc = (/0.d0, 0.d0, 0.d0/) ! tethered: no velocity

        if(root) write(*,'("x=",g12.4,"y=",g12.4,"z=",g12.4)') xc
        if(root) write(*,'("psi=",g12.4,"beta=",g12.4,"gamma=",g12.4)') psi,beta,gamma

    case ("yawpitchroll")
        psi      = 30.d0*pi/180.d0*sin(2.d0*pi*time)
        beta     = 30.d0*pi/180.d0*sin(2.d0*pi*time) ! pitch
        gamma    = 30.d0*pi/180.d0*sin(2.d0*pi*time)
        psi_dt   = 30.d0*pi/180.d0*cos(2.d0*pi*time)*2.d0*pi
        beta_dt  = 30.d0*pi/180.d0*cos(2.d0*pi*time)*2.d0*pi
        gamma_dt = 30.d0*pi/180.d0*cos(2.d0*pi*time)*2.d0*pi
        xc = Insect%x0 + time*Insect%v0
        vc = Insect%v0
        Insect%body_moves = "yes"

    case ("buffeting")
        psi      = -7.d0*pi/180.d0*sin(2.d0*pi*(23d0/152.d0)*time) ! roll
        beta     = -24.5d0*pi/180.d0 ! pitch
        gamma    = 180.d0*pi/180.d0 ! yaw
        psi_dt   = -7.d0*pi/180.d0*cos(2.d0*pi*(23.d0/152.d0)*time)*2.d0*pi*(23.d0/152.d0)
        beta_dt  = 0.d0
        gamma_dt = 0.d0
        xc = Insect%x0
        xc(2) = xc(2) - 0.45d0/13.2d0*sin(2.d0*pi*(23.d0/152.d0)*time)
        vc = (/0.0d0, -0.45d0/13.2d0*cos(2.d0*pi*(23.d0/152.d0)*time)*2.d0*pi*(23.d0/152.d0), 0.0d0/)
        Insect%body_moves = "yes"

    case ("casting")
        psi      = -20.d0*pi/180.d0*sin(2.d0*pi*(2d0/152.d0)*time) ! roll
        beta     = -24.5d0*pi/180.d0 ! pitch
        gamma    = 180.d0*pi/180.d0 ! yaw
        psi_dt   = -20.d0*pi/180.d0*cos(2.d0*pi*(2.d0/152.d0)*time)*2.d0*pi*(2.d0/152.d0)
        beta_dt  = 0.d0
        gamma_dt = 0.d0
        xc = Insect%x0
        xc(2) = xc(2) + 22.d0/13.2d0*sin(2.d0*pi*(2.d0/152.d0)*time)
        vc = (/0.0d0, 22.d0/13.2d0*cos(2.d0*pi*(2.d0/152.d0)*time)*2.d0*pi*(2.d0/152.d0), 0.0d0/)
        Insect%body_moves = "yes"

    case ("roll_param")
        ! Roll motion for aerodynamic power analysis
        ! Bumblebee model, wingbeat frequency 152 Hz, roll frequency 25 Hz, pitch angle 24.5 deg nose up
        !    psi      = -60.d0  *pi/180.d0*sin(2.d0*pi*(25d0/152.d0)*time) ! roll
        psi      = -60.d0  *pi/180.d0*sin(2.d0*pi*(6.25d0/152.d0)*time) ! roll
        beta     = -24.5d0*pi/180.d0 ! pitch
        gamma    = 180.d0*pi/180.d0 ! yaw
        !    psi_dt   = -60.d0  *pi/180.d0*cos(2.d0*pi*(25.d0/152.d0)*time)*2.d0*pi*(25.d0/152.d0)
        psi_dt   = -60.d0  *pi/180.d0*cos(2.d0*pi*(6.25d0/152.d0)*time)*2.d0*pi*(6.25d0/152.d0)
        beta_dt  = 0.d0
        gamma_dt = 0.d0
        xc = Insect%x0
        vc = (/0.0, 0.0, 0.0/) ! tethered: no velocity
        Insect%body_moves = "yes"

    case ("roll")
        psi      = 30.d0*pi/180.d0*sin(2.d0*pi*time)
        beta     = 0.d0 ! pitch
        gamma    = 0.d0 ! yaw
        psi_dt   = 30.d0*pi/180.d0*cos(2.d0*pi*time)*2.d0*pi
        beta_dt  = 0.d0
        gamma_dt = 0.d0
        xc = Insect%x0
        vc = (/0.0, 0.0, 0.0/) ! tethered: no velocity
        Insect%body_moves = "yes"

    case ("pitch_param")
        ! Pitch motion for aerodynamic power analysis
        psi      = Insect%yawpitchroll_0(3) ! roll
        beta     = Insect%yawpitchroll_0(2) +  15.d0  *pi/180.d0*cos(2.d0*pi*time) ! pitch
        gamma    = Insect%yawpitchroll_0(1) ! yaw
        psi_dt   = 0.d0
        beta_dt  = - 15.d0  *pi/180.d0*sin(2.d0*pi*time)*2.d0*pi
        gamma_dt = 0.d0
        xc = Insect%x0
        vc = (/0.0, 0.0, 0.0/) ! tethered: no velocity
        Insect%body_moves = "yes"

    case ("pitch")
        psi      = 0.d0
        beta     = 30.d0*pi/180.d0*sin(2.d0*pi*time)
        gamma    = 0.d0 ! yaw
        psi_dt   = 0.d0
        beta_dt  = 30.d0*pi/180.d0*cos(2.d0*pi*time)*2.d0*pi
        gamma_dt = 0.d0
        xc = Insect%x0
        vc = (/0.0, 0.0, 0.0/) ! tethered: no velocity
        Insect%body_moves = "yes"

    case ("yaw")
        psi      = 0.d0
        beta     = 0.d0
        gamma    = 30.d0*pi/180.d0*sin(2.d0*pi*time)
        psi_dt   = 0.d0
        beta_dt  = 0.d0
        gamma_dt = 30.d0*pi/180.d0*cos(2.d0*pi*time)*2.d0*pi
        xc = Insect%x0
        vc = (/0.0, 0.0, 0.0/) ! tethered: no velocity
        Insect%body_moves = "yes"

    case ("tethered")
        psi      = Insect%yawpitchroll_0(3) ! roll
        beta     = Insect%yawpitchroll_0(2) ! pitch
        gamma    = Insect%yawpitchroll_0(1) ! yaw
        psi_dt   = 0.d0  ! tethered: angles const
        beta_dt  = 0.d0
        gamma_dt = 0.d0
        xc = Insect%x0
        vc = (/0.0, 0.0, 0.0/) ! tethered: no velocity
        Insect%body_moves = "no" ! tethered: body does not move

        if (Insect%BodyType=="suzuki_thin_rod") then
            Insect%body_moves = "yes"
        endif



    case ("free_flight")
        ! The case "free_flight" is different from the others. The free flight solver
        ! computes the current state of the insect in INSECT%STATE, which is a 13
        ! component vector (6 translation, 4 quaternions, 3 angular velocity)
        ! in this case, the position is dynamically computed, and quaternions are used
        Insect%body_moves = "yes"
        Insect%quaternion_solver_used = .true.

        ! copy data from insect state vector
        xc = Insect%STATE(1:3)
        vc = Insect%STATE(4:6)
        ep = Insect%STATE(7:10)
        Insect%rot_body_b = Insect%STATE(11:13)

        ! compute yaw pitch roll from the quaternion. attention: this is just for
        ! information to dump in the log file (may be useful), but NOT for the rotation
        ! matrix. otherwise, we could just omit the quaternion
        psi      = atan2(2*(ep(2)*ep(3) + ep(0)*ep(1)), ep(0)*ep(0) - ep(1)*ep(1) - ep(2)*ep(2) + ep(3)*ep(3))
        beta     = asin(-2*(ep(1)*ep(3) - ep(0)*ep(2)))
        gamma    = atan2(2*(ep(1)*ep(2) + ep(0)*ep(3)), ep(0)*ep(0) + ep(1)*ep(1) - ep(2)*ep(2) - ep(3)*ep(3))

        ! these values cannot easily be computed, but they are not really necessary
        psi_dt   = 0.0d0
        beta_dt  = 0.0d0
        gamma_dt = 0.0d0

        if (Insect%wing_fsi /= "yes") then
            ! output on screen
            if(root) write(*,'("free--flight: ",f12.4,2x,9(es12.4,1x))') time, xc, vc, Insect%rot_body_b
        endif

    case default
        if (Insect%BodyMotion(1:11) == "from_file::") then
            ! use body state from a pre-existing simulation, which was probably an FSI
            ! run with free-flight. using this case, you can replay the same motion under
            ! different flow conditions
            call read_insect_STATE_from_file(time, Insect, Insect%BodyMotion(12:strlen), verbose=.false.)

            Insect%body_moves = "yes"
            Insect%quaternion_solver_used = .true.

            ! copy data from insect state vector
            xc = Insect%STATE(1:3)
            vc = Insect%STATE(4:6)
            ep = Insect%STATE(7:10)
            Insect%rot_body_b = Insect%STATE(11:13)

            ! compute yaw pitch roll from the quaternion. attention: this is just for
            ! information to dump in the log file (may be useful), but NOT for the rotation
            ! matrix. otherwise, we could just omit the quaternion
            psi      = atan2(2*(ep(2)*ep(3) + ep(0)*ep(1)), ep(0)*ep(0) - ep(1)*ep(1) - ep(2)*ep(2) + ep(3)*ep(3))
            beta     = asin(-2*(ep(1)*ep(3) - ep(0)*ep(2)))
            gamma    = atan2(2*(ep(1)*ep(2) + ep(0)*ep(3)), ep(0)*ep(0) + ep(1)*ep(1) - ep(2)*ep(2) - ep(3)*ep(3))

            ! these values cannot easily be computed, but they are not really necessary
            psi_dt   = 0.0d0
            beta_dt  = 0.0d0
            gamma_dt = 0.0d0

        else
            call abort(220918, "body_motion.f90::BodyMotion: motion case &
            &(Insect%BodyMotion) undefined: "//trim(adjustl(Insect%BodyMotion)))
        endif
    end select




    if ((root).and.(maxval(vc)>0.0d0).and.(Insect%body_moves=="no")) then
        write(*,*) "error in body_motion.f90: I found maxval(vc)>0 but the body_moves"
        write(*,*) "flag is set to no, which means we will draw the body only once"
        write(*,*) "This is probably not intented - you should look into it."
        call abort(10624, "error in body_motion.f90: I found maxval(vc)>0 but the body_moves")
    endif


    ! save above values in the insect
    Insect%psi      = psi
    Insect%beta     = beta
    Insect%gamma    = gamma
    Insect%psi_dt   = psi_dt
    Insect%beta_dt  = beta_dt
    Insect%gamma_dt = gamma_dt
    Insect%xc_body_g  = xc
    Insect%vc_body_g  = vc


    ! for compatibility, we update the x0,y0,z0 also
    ! this is used e.g. for torque computation
    x0 = xc(1)
    y0 = xc(2)
    z0 = xc(3)

end subroutine BodyMotion
