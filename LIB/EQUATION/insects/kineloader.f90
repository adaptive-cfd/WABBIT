subroutine load_kine_init(kine)
  implicit none

  type(wingkinematics), intent(inout) :: kine
  integer :: j, mpicode, nk
  real (kind=rk) :: t_period, r_wing


  if (root) then
    write(*,*) "Initializing kinematics loader for non-periodic kinematics"
    write(*,*) "file="//trim(adjustl(kine%infile))
    open (10, file = kine%infile, form='formatted', status='old')
    read (10, *) t_period ! stroke period in s, for normalization
    read (10, *) r_wing   ! wing length in mm, for normalization
    read (10, *) kine%nk
    if (kine%nk > nhrmt_max) then
      write(*,*) "WARNING(load_kine_init): adjusting the value kine%nk to nhrmt_max"
      kine%nk = nhrmt_max
    endif
    write(*,'("nk=",i5," t_period=",g12.4," r_wing=",g12.4)') kine%nk, t_period, r_wing
  endif

  call MPI_BCAST(kine%nk,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpicode)
  nk = kine%nk

  ! intitalize vectors as zero (note: not all of the space is actually used)
  kine%vec_t = 0.d0
  kine%vec_vert = 0.d0
  kine%vec_horz = 0.d0
  kine%vec_phi_dt = 0.d0
  kine%vec_alpha_dt = 0.d0
  kine%vec_theta_dt = 0.d0
  kine%vec_pitch_dt = 0.d0
  kine%vec_vert_dt = 0.d0
  kine%vec_horz_dt = 0.d0

  if (root) then
    ! read data from file
    do j = 1, nk
      read (10, *) kine%vec_t(j), &
        kine%vec_phi(j),&
        kine%vec_alpha(j),&
        kine%vec_theta(j),&
        kine%vec_pitch(j),&
        kine%vec_vert(j),&
        kine%vec_horz(j),  &
        kine%vec_phi_dt(j),&
        kine%vec_alpha_dt(j),&
        kine%vec_theta_dt(j),&
        kine%vec_pitch_dt(j),&
        kine%vec_vert_dt(j),&
        kine%vec_horz_dt(j)
    enddo
    close (10)
    print *, "load_kine_init: data read from file, nk=", nk
    print *, "non-dimensionalizing input data:"
    ! non-dimensionalize
    kine%vec_t = kine%vec_t / t_period
    kine%vec_vert = kine%vec_vert / r_wing
    kine%vec_horz = kine%vec_horz / r_wing
    kine%vec_phi_dt = kine%vec_phi_dt * t_period
    kine%vec_alpha_dt = kine%vec_alpha_dt * t_period
    kine%vec_theta_dt = kine%vec_theta_dt * t_period
    kine%vec_pitch_dt = kine%vec_pitch_dt * t_period
    kine%vec_vert_dt = kine%vec_vert_dt * t_period / r_wing
    kine%vec_horz_dt = kine%vec_horz_dt * t_period / r_wing
  endif

  call MPI_BCAST( kine%vec_t, nhrmt_max, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( kine%vec_phi, nhrmt_max, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( kine%vec_alpha, nhrmt_max, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( kine%vec_theta, nhrmt_max, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( kine%vec_pitch, nhrmt_max, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( kine%vec_vert, nhrmt_max, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( kine%vec_horz, nhrmt_max, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( kine%vec_phi_dt, nhrmt_max, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( kine%vec_alpha_dt, nhrmt_max, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( kine%vec_theta_dt, nhrmt_max, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( kine%vec_pitch_dt, nhrmt_max, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( kine%vec_vert_dt, nhrmt_max, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
  call MPI_BCAST( kine%vec_horz_dt, nhrmt_max, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )

  kine%initialized = .true.
  if (root) write(*,*) "done initializing kineoader!"
end subroutine


subroutine load_kine_clean(kine)
  implicit none
  type(wingkinematics), intent(inout) :: kine

  !if(allocated(kine%vec_t)) deallocate(kine%vec_t)
  !if(allocated(kine%vec_phi)) deallocate(kine%vec_phi)
  !if(allocated(kine%vec_phi_dt)) deallocate(kine%vec_phi_dt)
  !if(allocated(kine%vec_alpha)) deallocate(kine%vec_alpha)
  !if(allocated(kine%vec_alpha_dt)) deallocate(kine%vec_alpha_dt)
  !if(allocated(kine%vec_theta)) deallocate(kine%vec_theta)
  !if(allocated(kine%vec_theta_dt)) deallocate(kine%vec_theta_dt)
  !if(allocated(kine%vec_pitch)) deallocate(kine%vec_pitch)
  !if(allocated(kine%vec_pitch_dt)) deallocate(kine%vec_pitch_dt)
  !if(allocated(kine%vec_vert)) deallocate(kine%vec_vert)
  !if(allocated(kine%vec_vert_dt)) deallocate(kine%vec_vert_dt)
  !if(allocated(kine%vec_horz)) deallocate(kine%vec_horz)

end subroutine


subroutine wing_kine_interp(time, kine, phi_i, alpha_i, theta_i, phi_dt_i, alpha_dt_i, theta_dt_i)
  implicit none

  real(kind=rk), intent(in) :: time
  type(wingkinematics), intent(in) :: kine
  real(kind=rk), intent(out) :: phi_i,alpha_i,theta_i,phi_dt_i,alpha_dt_i,theta_dt_i

  if (kine%initialized .eqv. .false.) then
    call abort(1515,"kinematics_loader is not initialized but wing_kine_interp is called")
  endif

  if ((time<0.d0).or.(time>kine%vec_t(kine%nk))) then
    if(root) write(*,'("time=",es15.8)') time
    call abort(1516,"requested time in kineloader out of valid bounds")
  endif

  call hermite1d(kine%nk,kine%vec_t,kine%vec_phi  ,kine%vec_phi_dt,  time,phi_i,phi_dt_i)
  call hermite1d(kine%nk,kine%vec_t,kine%vec_alpha,kine%vec_alpha_dt,time,alpha_i,alpha_dt_i)
  call hermite1d(kine%nk,kine%vec_t,kine%vec_theta,kine%vec_theta_dt,time,theta_i,theta_dt_i)
end subroutine


! subroutine body_kine_interp(t_i,pitch_i,vert_i,horz_i,pitch_dt_i,vert_dt_i,horz_dt_i)
!   use kine
!   implicit none
!
!   real(kind=rk), intent(in) :: t_i
!   real(kind=rk), intent(out) :: pitch_i,vert_i,horz_i,pitch_dt_i,vert_dt_i,horz_dt_i
!
!   call hermite1d(nk,vec_t,vec_pitch,vec_pitch_dt,t_i,pitch_i,pitch_dt_i)
!   call hermite1d(nk,vec_t,vec_vert,vec_vert_dt,t_i,vert_i,vert_dt_i)
!   call hermite1d(nk,vec_t,vec_horz,vec_horz_dt,t_i,horz_i,horz_dt_i)
! end subroutine


subroutine hermite1d(n, xphi, phi, dpdx, xi, phi_interp, dpdx_interp)
  implicit none

  integer, parameter :: rk = 8
  integer :: i0,i1
  integer, intent(in) :: n
  real(kind=rk) :: x,z,ap0,ap1,ax0,ax1,dx,d2pdx2_i0,d2pdx2_i1
  real(kind=rk), intent(in) :: xi
  real(kind=rk), intent(in) :: xphi(1:n),phi(1:n),dpdx(1:n)
  real(kind=rk), intent(out) :: phi_interp,dpdx_interp

  dx = xphi(2) - xphi(1)

  i0 = floor(xi/dx)+1
  i1 = i0+1

  if ((xi<xphi(i0)).or.(xi>xphi(i1))) then
     print *, "hermite1d: not a uniform grid"
     write(*,'("xi=",es12.4," dx=",es12.4)') xi,dx
     write(*,'("xphi(i0)=",es12.4," xphi(i1)=",es12.4)') xphi(i0),xphi(i1)
     call abort(884, "hermite1d: not a uniform grid")
  endif

  x = (xi-xphi(i0))/dx
  z = 1.0d0-x

  ! Phi
  ap0 = 2.0d0*x*x*x-3.0d0*x*x+1.0d0 ! f(x)
  ap1 = 2.0d0*z*z*z-3.0d0*z*z+1.0d0 ! f(1.0d0-x)
  ax0 = x*x*x  -2.0d0*x*x + x !g(x)
  ax1 = -(z*z*z -2.0d0*z*z + z) !-g(1.0d0-x)
  phi_interp = (phi(i0)*ap0+phi(i1)*ap1) + dx*(dpdx(i0)*ax0+dpdx(i1)*ax1)

  ! Phi_x
  if ((i0>1).and.(i1<n)) then
     d2pdx2_i0 = 0.5d0*(dpdx(i1)-dpdx(i0-1))/dx
     d2pdx2_i1 = 0.5d0*(dpdx(i1+1)-dpdx(i0))/dx
     dpdx_interp=(dpdx(i0)*ap0+dpdx(i1)*ap1)+dx*(d2pdx2_i0*ax0+d2pdx2_i1*ax1)
  else
     ap0 = 6.0d0*x*x-6.0d0*x !df(x)
     ap1 = -(6.0d0*z*z-6.0d0*z) !-df(1.0d0-x)
     ax0 = 3.0d0*x*x-4.0d0*x+1.0d0 !dg(x)
     ax1 = 3.0d0*z*z-4.0d0*z+1.0d0 ! dg(1.0d0-x)
     dpdx_interp = (phi(i0)*ap0+phi(i1)*ap1)/dx+(dpdx(i0)*ax0+dpdx(i1)*ax1)
  endif

end subroutine
