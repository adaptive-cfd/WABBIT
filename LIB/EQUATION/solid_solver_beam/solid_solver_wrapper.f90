!-------------------------------------------------------------------------------
!   SOLID SOLVER WRAPPER
! Input:
!   time: current time level (t^n)
!   dt: time step for new level (t^n+1 - t^n)
!   Beams: array of beams, all at the old time level t^n
! Output:
!   Beams: array of beams, all at the new time level t^n+1
!-------------------------------------------------------------------------------
subroutine SolidSolverWrapper ( time, dt )
    implicit none
    real(kind=rk), intent (in) ::  dt, time
    real(kind=rk) :: t0, soft_startup, dh
    integer :: i, is, ih, isurf, k
    t0 = MPI_wtime()


    dh = dx_min

    ! 1/ unpacking of pq_solid_old to p_surface
    ! The interpolation routine gets a (1:3, 1:N) array, which is the optimized
    ! way to store such data (1st index - fast in fortran). The solid solver uses
    ! a different array ordering, and we first unpack the data into this format.
    k = 1
    do is = 0, ns-1
        do ih = 0, size(p_surface, 2)-1
            do isurf = 1, 2
                p_surface(is, ih, isurf) = pq_solid_old(1, k)
                k = k +1
            enddo
        enddo
    enddo
    ! 2/ Integrate the pressure over the spanwise direction.
    ! The beam model is 1D, so it is flexible only in one direction. In the other, it
    ! is rigid. In this rigid direction, we integrate the pressure using a simple
    ! quadrature formula.
    ! NOTE: original FLUSI code applied masking to points here for non-rect shapes
    do is=0,ns-1
        !-- store result of interpolation at OLD timelevel
        beams(1)%pressure_old(is) = sum(p_surface(is,:,1) - p_surface(is,:,2)) * dh
    enddo

    ! 1/ for new time level
    k = 1
    do is = 0, ns-1
        do ih = 0, size(p_surface, 2)-1
            do isurf = 1, 2
                p_surface(is, ih, isurf) = pq_solid_new(1, k)
                k = k +1
            enddo
        enddo
    enddo
    ! 2/ for new time level
    do is=0,ns-1
        !-- store result of interpolation at NEW timelevel
        beams(1)%pressure_new(is) = sum(p_surface(is,:,1) - p_surface(is,:,2)) * dh
    enddo

    ! viscous tensions, from earler 2D versions, are not used in this implementation right
    ! now.
    beams(1)%tau_new = 0.0_rk
    beams(1)%tau_old = 0.0_rk

    !-----------------------------------------------------------------------------
    ! To avoid startup problems, we can smoothly "turn on" the coupling by
    ! multiplying the pressure with a startup conditioner.
    ! Note this does not apply for gravity or imposed motion.
    !-----------------------------------------------------------------------------
    soft_startup = startup_conditioner(time, T_release, tau)
    beams(1)%pressure_new = beams(1)%pressure_new * soft_startup
    beams(1)%pressure_old = beams(1)%pressure_old * soft_startup

    do i = 1, nBeams
        !------------------------------------------
        ! check if input values are okay
        !------------------------------------------
        if (Vector_isNAN(beams(i)%pressure_new).or.&
            Vector_isNAN(beams(i)%pressure_old).or.&
            Vector_isNAN(beams(i)%tau_new).or.&
            Vector_isNAN(beams(i)%tau_old) ) then
            if (root) write(*,*) "SolidSolver: input values contain NaNs"
            ! time to go..
            call abort(6652,"The input forces for the solid solver contain NaN..abort")
        endif

        if ( time >= T_release) then
            !-------------------------------------------
            ! the beams are released (now: active FSI), call IBES solvers
            !-------------------------------------------
            ! NOTE: IBES contains integrate_position
            select case (TimeMethodSolid)
              case ("CN2","BDF2","EI1")
                ! all implicit solvers are in one subroutine
                call IBES_solver (time, dt, beams(i))
              case ("RK4")
                ! note this solver may be very unstable:
                ! call RK4_wrapper (time, dt, beams(i))
                call abort(99987,"RK4 for beam equation is not implemented")
              case ("EE1")
                ! note this solver may be very unstable:
                ! call EE1_wrapper (time, dt, beams(i))
                call abort(7778,"EE1 for beam equation is not implemented")
              case default
                call abort(723763,"SolidSolverWrapper::invalid value of TimeMethodSolid"//&
                    trim(adjustl(TimeMethodSolid)))
            end select
        else
            !-------------------------------------------
            ! the beams are not yet released, but their leading edges may move
            ! (passive FSI)
            !-------------------------------------------
            call integrate_position (time+dt, beams(i))
        endif

        !-------------------------------------------
        ! compute energies and stuff
        !-------------------------------------------
        call SolidEnergies ( beams(i) )

        !-- check if everything seems okay, if not show beam and abort
        call show_beam_on_error( beams(i) )
    enddo

    call toc("SolidSolver (SolidSolverWrapper)", 223, MPI_wtime() - t0)
end subroutine SolidSolverWrapper


!-------------------------------------------------------------------------------
!   solid solver main entry point
!-------------------------------------------------------------------------------
subroutine OnlySolidSimulation()
    !use vars
    implicit none
    ! type(solid), dimension(1:nBeams) :: beams
    ! real (kind=rk) :: time
    ! integer :: it,nsave

    ! ! in case you forgot to set it
    ! if (itdrag==0) itdrag=10
    ! if (dt_fixed<=1.0d-10) dt_fixed=1.d-3

    ! write (*,*) "*** information: starting OnlySolidSimulation"
    ! time = 0.0
    ! it = 0

    ! call show_solid_model_information

    ! !-- initialization
    ! call init_beams( beams )

    ! !--loop over time steps
    ! do while ((time<=tmax))
    !   !-- external loads for testing purposes
    !   if (debug_pressure==1) then
    !     beams(1)%pressure_new=0.1*dsin(time)
    !     beams(1)%tau_new=0.1*dsin(time)
    !   endif

    !   !-- time stepping
    !   call SolidSolverWrapper( time, dt_fixed , beams )

    !   it   = it+1
    !   time = dble(it)*dt_fixed

    !   call SaveBeamData( time, beams )
    ! enddo

end subroutine


!-------------------------------------------------------------------------------
!   solid solver convergence test
!-------------------------------------------------------------------------------
subroutine SolidModelConvergenceTest()
    ! !use vars
    ! implicit none
    ! type(solid), dimension(1:nBeams) :: beams
    ! type(solid), dimension(1:6) :: solutions
    ! real(kind=rk) :: time, t1
    ! real(kind=rk), dimension(1:6) :: dts, errors_y, errors_vy, errors_th
    ! integer :: it,nsave,i

    ! ! in case you forgot to set it
    ! if (itdrag==0) itdrag=10
    ! if (dt_fixed<=1.0d-10) dt_fixed=1.d-3

    ! write(*,*) "*** information: starting SolidModelConvergenceTest"
    ! dts = (/2.0d-3, 1.0d-3, 5.0d-4, 2.5d-4, 1.0d-4, 5.0d-5/)
    ! write(*,*) "We use the folling time steps: ", dts

    ! !-----------------------------------------------------------------------------
    ! ! loop over the different time steps and perform one simulation for each value of dt
    ! !-----------------------------------------------------------------------------
    ! do i = 1, 6
    !   t1 = MPI_wtime()
    !   dt_fixed = dts(i)
    !   time = 0.0d0
    !   it = 0
    !   call show_solid_model_information
    !   !-- initialization
    !   call init_beams( beams )
    !   !--loop over time steps
    !   do while (time<tmax)
    !     !-- time stepping
    !     call SolidSolverWrapper( time, dt_fixed , beams )
    !     it   = it+1
    !     time = dble(it)*dt_fixed
    !   enddo

    !   !-- run is done, copy the result
    !   solutions(i) = beams(1)
    !   write(*,'(80("-"))')
    !   write(*,*) "run dt=",dts(i)," is done time=", time, "it took t=", MPI_wtime()-t1
    !   write(*,'(80("-"))')
    ! enddo

    ! !-----------------------------------------------------------------------------
    ! ! compute errors, save them in a txt file and show them on screen
    ! !-----------------------------------------------------------------------------
    ! call init_empty_file('time_convergence_result')
    ! open (14, file = 'time_convergence_result', status = 'unknown',position='append')
    ! write(14,'("% dt   error_y    error_vy  error_theta")')
    ! do i=5,1,-1
    !   errors_y(i) = dsqrt(sum(solutions(i)%y - solutions(6)%y)**2) / dsqrt(sum((solutions(6)%y)**2))
    !   errors_th(i) = dsqrt(sum(solutions(i)%theta - solutions(6)%theta)**2) / dsqrt(sum((solutions(6)%theta)**2))
    !   errors_vy(i) = dsqrt(sum(solutions(i)%vy - solutions(6)%vy)**2) / dsqrt(sum((solutions(6)%vy)**2))

    !   write(*,'(4(es15.8))')  dts(i), errors_y(i), errors_vy(i), errors_th(i)
    !   write(14,'(4(es15.8))') dts(i), errors_y(i), errors_vy(i), errors_th(i)
    ! enddo
    ! close(14)

end subroutine
