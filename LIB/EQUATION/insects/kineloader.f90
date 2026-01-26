subroutine load_kine_init(Insect)
  implicit none

  type(diptera), intent(inout) :: Insect

  integer :: nk, nc

  if (root) then
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,*) "Initializing kinematics loader for non-periodic kinematics"
    write(*,*) "file="//trim(adjustl(Insect%infile_kineloader))
    write(*,*) "Insect%BodyMotion=", trim(adjustl(Insect%BodyMotion))
  endif

  call count_lines_in_ascii_file_mpi(Insect%infile_kineloader, nk, n_header=1)
  call count_cols_in_ascii_file_mpi(Insect%infile_kineloader, nc, n_header=1)

  if (nc /= 43) then
    write(*,*) nc
    call abort(1111251,"The kineloader file should have 43 columns (see f90 ode for naming of indiv. columns)")
  endif

  ! number of time steps read from file
  Insect%nk = nk

  if (root) then
    write(*,*) "nk=", nk
    write(*,*) "nc=", nc
  endif

  ! 1   0   time
  ! 2   1   body_center_g_x
  ! 3   2   body_center_g_x_dt
  ! 4   3   body_center_g_y
  ! 5   4   body_center_g_y_dt
  ! 6   5   body_center_g_z
  ! 7   6   body_center_g_z_dt
  ! 8   7   unused
  ! 9   8   unused
  ! 10  9   unused
  ! 11  10  unused
  ! 12  11  unused
  ! 13  12  unused
  ! 14  13  psi
  ! 15  14  psi_dt
  ! 16  15  beta
  ! 17  16  beta_dt
  ! 18  17  gamma
  ! 19  18  gamma_dt
  ! 20  19  alpha_L
  ! 21  20  alpha_L_dt
  ! 22  21  phi_L
  ! 23  22  phi_L_dt
  ! 24  23  theta_L
  ! 25  24  theta_L_dt
  ! 26  25  alpha_R
  ! 27  26  alpha_R_dt
  ! 28  27  phi_R
  ! 29  28  phi_R_dt
  ! 30  29  theta_R
  ! 31  30  theta_R_dt
  ! 32  31  alpha_L2
  ! 33  32  alpha_L2_dt
  ! 34  33  phi_L2
  ! 35  34  phi_L2_dt
  ! 36  35  theta_L2
  ! 37  36  theta_L2_dt
  ! 38  37  alpha_R2
  ! 39  38  alpha_R2_dt
  ! 40  39  phi_R2
  ! 41  40  phi_R2_dt
  ! 42  41  theta_R2
  ! 43  42  theta_R2_dt
  allocate( Insect%data_kineloader(1:nk, 1:43) )

  call read_array_from_ascii_file_mpi(Insect%infile_kineloader, Insect%data_kineloader, n_header=1)

  if (root) then
    write(*,'("time_max_kineloader=",es15.8)') Insect%data_kineloader(Insect%nk-3:Insect%nk,1)
    write(*,*) "done initializing kineloader!"
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  endif

  Insect%kineloader_initialized = .true.  
end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine wing_kine_interp(time, Insect, wingID, phi, alpha, theta, phi_dt, alpha_dt, theta_dt)
  implicit none

  real(kind=rk), intent(in) :: time
  type(diptera), intent(inout) :: Insect
  ! wingID: 1=left 2=right 3=left2 4=right2
  integer(kind=2), intent(in) :: wingID
  real(kind=rk), intent(out) :: phi, alpha, theta, phi_dt, alpha_dt, theta_dt

  ! column indices for interpolation (see documentaion above), one-based indexing
  integer(kind=ik) :: ialpha, itheta, iphi

  if (.not. Insect%kineloader_initialized) then
    call abort(3010251,"kinematics_loader is not initialized but wing_kine_interp is called")
  endif

  if ( (time<0.0_rk) .or. (time>Insect%data_kineloader(Insect%nk,1)) ) then
    if (root) then
      write(*,'("time=",es15.8)') time
      write(*,'("time_max_kineloader=",es15.8)') Insect%data_kineloader(Insect%nk,1)
    endif
    call abort(3010253,"requested time in kineloader out of valid bounds")
  endif

  select case ( wingID )
  case (1) !("left")
    ialpha = 20
    itheta = 24
    iphi   = 22
  case (2) !("right")
    ialpha = 26
    itheta = 30
    iphi   = 28
  case (3) !("left2")
    ialpha = 32
    itheta = 36
    iphi   = 34
  case (4) !("right2")
    ialpha = 38
    itheta = 42
    iphi   = 40
  case default
    call abort(77747, "not a valid wing identifier")
  end select

  ! apply hermite interpolation on data columns to get wing angles and their time derivatives
  call hermite1d(Insect%nk, Insect%data_kineloader(:,1), Insect%data_kineloader(:,iphi), Insect%data_kineloader(:,iphi+1), time, phi, phi_dt)
  call hermite1d(Insect%nk, Insect%data_kineloader(:,1), Insect%data_kineloader(:,ialpha), Insect%data_kineloader(:,ialpha+1), time, alpha, alpha_dt)
  call hermite1d(Insect%nk, Insect%data_kineloader(:,1), Insect%data_kineloader(:,itheta), Insect%data_kineloader(:,itheta+1), time, theta, theta_dt)
end subroutine


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine body_kine_interp(time, Insect, gamma, beta, psi, gamma_dt, beta_dt, psi_dt, xc, yc, zc, ux, uy, uz )
  implicit none

  real(kind=rk), intent(in) :: time
  type(diptera), intent(inout) :: Insect
  real(kind=rk), intent(out) :: gamma, beta, psi, gamma_dt, beta_dt, psi_dt, xc, yc, zc, ux, uy, uz

  integer(kind=ik) :: ipsi=14, igamma=18, ibeta=16
  integer(kind=ik) :: ix=2, iy=4, iz=6!, iux=8, iuy=10, iuz=12


  if (.not. Insect%kineloader_initialized) then
    call abort(3010251,"kinematics_loader is not initialized but wing_kine_interp is called")
  endif

  if ( (time<0.0_rk) .or. (time>Insect%data_kineloader(Insect%nk,1)) ) then
    if (root) then
      write(*,'("time=",es15.8)') time
      write(*,'("time_max_kineloader=",es15.8)') Insect%data_kineloader(Insect%nk,1)
    endif
    call abort(3010253,"requested time in kineloader out of valid bounds")
  endif


  ! body angles and their time derivatives
  call hermite1d(Insect%nk, Insect%data_kineloader(:,1), Insect%data_kineloader(:,ipsi), Insect%data_kineloader(:,ipsi+1), time, psi, psi_dt)
  call hermite1d(Insect%nk, Insect%data_kineloader(:,1), Insect%data_kineloader(:,ibeta), Insect%data_kineloader(:,ibeta+1), time, beta, beta_dt)
  call hermite1d(Insect%nk, Insect%data_kineloader(:,1), Insect%data_kineloader(:,igamma), Insect%data_kineloader(:,igamma+1), time, gamma, gamma_dt)

  ! body position & velocity
  ! (in global coordinate system)
  call hermite1d(Insect%nk, Insect%data_kineloader(:,1), Insect%data_kineloader(:,ix), Insect%data_kineloader(:,ix+1), time, xc, ux)
  call hermite1d(Insect%nk, Insect%data_kineloader(:,1), Insect%data_kineloader(:,iy), Insect%data_kineloader(:,iy+1), time, yc, uy)
  call hermite1d(Insect%nk, Insect%data_kineloader(:,1), Insect%data_kineloader(:,iz), Insect%data_kineloader(:,iz+1), time, zc, uz)
end subroutine

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

subroutine hermite1d(n, xphi, phi, dpdx, xi, phi_interp, dpdx_interp)
  implicit none

  integer, parameter :: rk = 8
  integer(kind=ik) :: i0,i1
  integer(kind=ik), intent(in) :: n
  real(kind=rk) :: x,z,ap0,ap1,ax0,ax1,dx,d2pdx2_i0,d2pdx2_i1
  real(kind=rk), intent(in) :: xi
  real(kind=rk), intent(in) :: xphi(1:n),phi(1:n),dpdx(1:n)
  real(kind=rk), intent(out) :: phi_interp,dpdx_interp
  real(kind=rk), parameter :: eps = 1.0e-11_rk

  dx = xphi(2) - xphi(1)

  i0 = floor(xi/dx)+1
  i1 = i0+1

  if ((xi<xphi(i0)-eps).or.(xi>xphi(i1)+eps)) then
     print *, "hermite1d: not a uniform grid"
     write(*,'("xi=",es17.9," dx=",es17.9)') xi,dx
     write(*,'("xphi(i0)=",es17.9," xphi(i1)=",es17.9)') xphi(i0),xphi(i1)
     call abort(202505111, "hermite1d: not a uniform grid")
  endif

  x = (xi-xphi(i0))/dx
  z = 1.0_rk - x

  ! Phi
  ap0 = 2.0_rk*x*x*x-3.0_rk*x*x+1.0_rk ! f(x)
  ap1 = 2.0_rk*z*z*z-3.0_rk*z*z+1.0_rk ! f(1.0_rk-x)
  ax0 = x*x*x  -2.0_rk*x*x + x !g(x)
  ax1 = -(z*z*z -2.0_rk*z*z + z) !-g(1.0_rk-x)
  phi_interp = (phi(i0)*ap0+phi(i1)*ap1) + dx*(dpdx(i0)*ax0+dpdx(i1)*ax1)

  ! Phi_x
  if ((i0>1).and.(i1<n)) then
     d2pdx2_i0 = 0.5_rk*(dpdx(i1)-dpdx(i0-1))/dx
     d2pdx2_i1 = 0.5_rk*(dpdx(i1+1)-dpdx(i0))/dx
     dpdx_interp=(dpdx(i0)*ap0+dpdx(i1)*ap1)+dx*(d2pdx2_i0*ax0+d2pdx2_i1*ax1)
  else
     ap0 = 6.0_rk*x*x-6.0_rk*x !df(x)
     ap1 = -(6.0_rk*z*z-6.0_rk*z) !-df(1.0_rk-x)
     ax0 = 3.0_rk*x*x-4.0_rk*x+1.0_rk !dg(x)
     ax1 = 3.0_rk*z*z-4.0_rk*z+1.0_rk ! dg(1.0_rk-x)
     dpdx_interp = (phi(i0)*ap0+phi(i1)*ap1)/dx+(dpdx(i0)*ax0+dpdx(i1)*ax1)
  endif

end subroutine
