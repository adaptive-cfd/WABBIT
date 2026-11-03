!--------------------------------------------
! solves the linear system A*x = b
! assumed shape-array and vector size
!--------------------------------------------
subroutine Solve_linear_system ( A, b, x )
    implicit none
    real(kind=rk),dimension(1:,1:), intent(in) :: A
    real(kind=rk),dimension(1:), intent(inout) :: x
    real(kind=rk),dimension(1:), intent(in) :: b
    real(kind=rk) :: t0
    integer :: error, nn
    integer, allocatable, dimension(:) :: ipiv

    nn = size(b)

    if (size(b) /= size(x) ) call abort(91,'Solve_linear_system: size(b)/=size(x)')

    allocate(ipiv(1:nn))

    t0 = MPI_wtime()
    call dgetrf( nn, nn, A , nn, ipiv, error )
    if (error .ne. 0) then
        write(*,*) "!!! mmCrutial: dgetrf error.", error ,nn
        call abort(92, "Error in solve liner system (dgetrf)")
    endif

    x = b
    call dgetrs( 'N', nn, 1, A, nn, ipiv, x, nn, error )
    if (error .ne. 0) then
        write(*,*) "!!! mmCrutial: dgetrs error.", error ,nn
        call abort(93, "Error in solve liner system (dgetrs)")
    endif


    call toc("SolidSolver (Solve_linear_system)", 222, MPI_wtime()-t0)
end subroutine

! ------------------------------------------------------------------------------

subroutine time_marching_coefs(dt,dt_old,C1,C2,C3,C4)
    implicit none
    real(kind=rk) :: R,dt,dt_old,C1,C2,C3,C4

    if (TimeMethodSolid == "EI1") then
        C1=1.0d0 ! dt factor
        C2=0.d0 ! factor for RHS
        C3=1.0d0 ! factor before the THETA_DOT_N term
        C4=0.d0 ! factor before the THETA_DOT_N-1 term
    elseif (TimeMethodSolid == "CN2") then
        C1=2.0d0 ! dt factor
        C2=1.0d0 ! rhs old factor
        C3=1.0d0 ! factor before the THETA_DOT_N term
        C4=0.d0 ! factor before the THETA_DOT_N-1 term
    elseif (TimeMethodSolid == "BDF2") then
        R  = dt / dt_old
        C1 = (1.d0+2.d0*R)/(1.d0+R)   ! dt factor
        C2 = 0.d0      ! rhs old factor
        C3 = ((1.d0+R)**2)/(1.d0+2.d0*R)   ! factor before the THETA_DOT_N term
        C4 = (-R**2 )/(1.d0+2.d0*R)   ! factor before the THETA_DOT_N-1 term
    endif

    if (C2.ne.0.0) call abort(94,"to use c2/=0 you have to implement the rhs first")
end subroutine



! ------------------------------------------------------------------------------

subroutine Differentiate1D (f, f_derivative, N, dx, order)
    implicit none
    integer, intent (in)         :: order, n
    real(kind=rk), intent (in)            :: dx
    real(kind=rk), dimension(0:N-1)      :: f, f_derivative
    real(kind=rk),dimension(0:N-1, 0:N-1)   :: D1, D2, D3 ,D4

    call create_diff_matrices (D1, D2, D3, D4, N, dx)

    select case (order)
      case(1)
        f_derivative = matmul(D1,f)
      case(2)
        f_derivative = matmul(D2,f)
      case(3)
        f_derivative = matmul(D3,f)
      case(4)
        f_derivative = matmul(D4,f)
      case default
        call abort(95, "Differentiate1D wrong choice")
    end select

end subroutine Differentiate1D

! ------------------------------------------------------------------------------

subroutine create_diff_matrices (D1, D2, D3, D4, nn, dx)
    implicit none
    integer :: i,N
    integer, intent (in) :: nn
    real(kind=rk), intent(in) :: dx
    real(kind=rk),dimension(0:nn-1, 0:nn-1), intent (out) :: D1, D2, D3 ,D4

    !------------------------------------
    ! create diff-matrices.
    ! first derivative: x'=D1 * x
    !------------------------------------

    D1=0.d0
    D2=0.d0
    D3=0.d0
    D4=0.d0

    do i=0,nn-1 !main diag
        D1(i,i) = 0.d0 !main
        D2(i,i) =-2.d0 !main
        D3(i,i) = 0.d0 !main
        D4(i,i) = 6.d0 !main
    enddo

    do i=0,nn-2 !first upper/lower diag
        D1(i,i+1) = 1.d0  !up1
        D1(i+1,i) =-1.d0  !lo1
        D2(i,i+1) = 1.d0  !up1
        D2(i+1,i) = 1.d0  !lo1
        D3(i,i+1) =-2.d0  !up1
        D3(i+1,i) = 2.d0  !lo1
        D4(i,i+1) =-4.d0  !up1
        D4(i+1,i) =-4.d0  !lo1
    enddo

    do i=0,nn-3 !second upper/lower diag
        D3(i,i+2) = 1.d0   !up2
        D3(i+2,i) =-1.d0  !lo2
        D4(i,i+2) = 1.d0   !up2
        D4(i+2,i) = 1.d0   !lo2
    enddo
!   N=ns+2 !dimension of the matrices (we added 3 points for 3 boundary equations
    N = nn-1
!-----------------------------
    D1(0,0)=  -3.d0
    D1(0,1)=   4.d0
    D1(0,2)=  -1.d0
    D1(N,N-2)= 1.d0
    D1(N,N-1)=-4.d0
    D1(N,N)=   3.d0
    D1=D1/(2.0d0* dx)
!-----------------------------
    D2(0,0)=2.d0
    D2(0,1)=-5.d0
    D2(0,2)=4.d0
    D2(0,3)=-1.d0
    D2(N,N-3)=-1.d0
    D2(N,N-2)=4.d0
    D2(N,N-1)=-5.d0
    D2(N,N)=2.d0
    D2=D2/( dx**2 )
!-----------------------------
    D3(0,0)=-17.0d0/2.d0
    D3(0,1)=71.0d0/2.d0
    D3(0,2)=-59.d0
    D3(0,3)=49.d0
    D3(0,4)=-41.0d0/2.d0
    D3(0,5)=7.0d0/2.d0
    D3(1,0)=0.d0
    D3(1,1)=-17.0d0/2.d0
    D3(1,2)=71.0d0/2.d0
    D3(1,3)=-59.d0
    D3(1,4)=49.d0
    D3(1,5)=-41.0/2.d0
    D3(1,6)=7.0d0/2.d0
    D3(N,N-5)=-7.0d0/2.d0
    D3(N,N-4)=41.0d0/2.d0
    D3(N,N-3)=-49.d0
    D3(N,N-2)=59.d0
    D3(N,N-1)=-71.0d0/2.d0
    D3(N,N)=17.0d0/2.d0
    D3(N-1,N-6)=-7.0d0/2.d0
    D3(N-1,N-5)=41.0d0/2.d0
    D3(N-1,N-4)=-49.0d0
    D3(N-1,N-3)=59.0d0
    D3(N-1,N-2)=-71.0/2.d0
    D3(N-1,N-1)=17.0/2.d0
    D3(N-1,N)=0.d0
    D3=D3/(2.0d0* dx**3)
!-----------------------------
    D4(0,0)=35.0d0/6.d0
    D4(0,1)=-31.d0
    D4(0,2)=137.0d0/2.d0
    D4(0,3)=-242.0d0/3.d0
    D4(0,4)=107.0d0/2.d0
    D4(0,5)=-19.0d0
    D4(0,6)=17.0d0/6.d0
    D4(1,0)=0.d0
    D4(1,1)=35.0d0/6.d0
    D4(1,2)=-31.d0
    D4(1,3)=137.0d0/2.d0
    D4(1,4)=-242.0d0/3.0d0
    D4(1,5)=107.0d0/2.d0
    D4(1,6)=-19.0d0
    D4(1,7)=17.0d0/6.0d0
    D4(N,N-6)=17.0d0/6.0d0
    D4(N,N-5)=-19.0d0
    D4(N,N-4)=107.0d0/2.0d0
    D4(N,N-3)=-242.0d0/3.0d0
    D4(N,N-2)=137.0d0/2.0d0
    D4(N,N-1)=-31.0d0
    D4(N,N)=35.0d0/6.0d0
    D4(N-1,N-7)=17.0d0/6.d0
    D4(N-1,N-6)=-19.d0
    D4(N-1,N-5)=+107.0d0/2.d0
    D4(N-1,N-4)=-242.0d0/3.d0
    D4(N-1,N-3)=+137.0d0/2.d0
    D4(N-1,N-2)=-31.0d0
    D4(N-1,N-1)=+35.0d0/6.d0
    D4(N-1,N)=0.d0
    D4=D4/( dx**4)
    !-----------------------------
end subroutine create_diff_matrices

!-------------------------------------------------------------------------------


logical function Vector_isNAN(f)
    real(kind=rk), intent(in) :: f(:)
    integer :: a, i
    a = size(f)
    Vector_isNAN = .false.
    do i=1,a
        if (.not.(f(i).eq.f(i))) then
            Vector_isNAN = .true.
            return
        endif
    enddo
end function


!-------------------------------------------------------------------------------


subroutine show_solid_model_information
    implicit none
    write(*,*) "*****************************************"
    write(*,*) "solid solver debuging: showing all globals"
    write(*,*) "*****************************************"
    write(*,'("mpirank=",i5)') mpirank
    write(*,'("ns=",i3)') ns
    write(*,'("mue0=",es12.4," eta0=",es12.4," grav=",es12.4)') mue0,eta0,grav
    write(*,'("sigma=",es12.4," t_beam=",es12.4," ds=",es12.4)') sigma, t_beam, ds
    write(*,'("T_release=",es12.4," TimeMethodSolid=",A," tau=",es12.4)') &
        T_release, trim(TimeMethodSolid), tau
    write(*,*) "*******************end*******************"
end subroutine


!-------------------------------------------------------------------------------

subroutine show_beam( beam )
    !use vars
    implicit none
    character(len=20):: str
    type(solid),intent(in) :: beam
    if (root) then
        write(str,'("(A,1x,",i3.3,"(f7.3,1x))")') ns
        write(*,str) "beam%x ", beam%x(0:ns-1)
        write(*,str) "beam%y ", beam%y(0:ns-1)
        write(*,str) "beam%vx", beam%vx(0:ns-1)
        write(*,str) "beam%vy", beam%vy(0:ns-1)
        write(*,*) "dt_old", beam%dt_old
        write(*,str) "beam%theta    ", beam%theta(0:ns-1)
        write(*,str) "beam%theta_dot", beam%theta_dot(0:ns-1)
        write(*,str) "beam%pressure_new", beam%pressure_new(0:ns-1)
        write(*,str) "beam%pressure_old", beam%pressure_old(0:ns-1)
        write(*,str) "beam%tau_new", beam%tau_new(0:ns-1)
        write(*,str) "beam%tau_old", beam%tau_old(0:ns-1)
        write(*,*) "-->end beam<--"
    endif
end subroutine show_beam


!-------------------------------------------------------------------------------


subroutine show_beam_on_error( beam )
    !use vars
    implicit none
    type(solid),intent(in) :: beam
    if ( Vector_isNAN(beam%x).or.&
        Vector_isNAN(beam%y).or.&
        Vector_isNAN(beam%vx).or.&
        Vector_isNAN(beam%vy).or.&
        Vector_isNAN(beam%theta).or.&
        Vector_isNAN(beam%theta_dot).or.&
        Vector_isNAN(beam%pressure_new).or.&
        Vector_isNAN(beam%pressure_old).or.&
        Vector_isNAN(beam%tau_new).or.&
        Vector_isNAN(beam%tau_old) &
        ) then

        call show_beam( beam )
        call abort(171717, "Show beam on error lets flusi quit here.")
    endif
end subroutine show_beam_on_error


!-------------------------------------------------------------------------------


subroutine lapack_unit_test()
    !-----------------------------------------------------------------------------
    ! Minimal example for testing the Lapack library
    ! Solves the linear system
    ! I*x = b
    ! where I is the identity matrix, so the solution is trivially x=b
    !-----------------------------------------------------------------------------
    !use vars
    implicit none
    integer, parameter:: n = 20
    real(kind=rk) :: a1(1:n,1:n),a2(1:n,1:n), b(1:n), x(1:n), err
    integer :: i,error
    integer :: ipiv(1:n)
    a2 = 0.d0
    a1 = 0.d0
    b  = 7.d0
    x  = 0.d0

    if (root) write(*,*) "--------------------------------"
    if (root) write(*,*) " Starting LAPACK unit test"

    ! create identity matrix
    do i=1,n
        a1(i,i) = 1.d0
    enddo

    !-----------------------------------------------------------------------------
    ! first step: factorization of the matrix a1 (which is overwritten on output
    ! so we make a copy of it)
    !-----------------------------------------------------------------------------
    a2=a1 ! lapack overwrites a2 with the factorization
    call dgetrf ( n, n, a2, n, ipiv, error )

    !-----------------------------------------------------------------------------
    ! second step: backwards substitution
    !-----------------------------------------------------------------------------
    x = b ! lapack overwrites x with the solution
    call dgetrs( 'N', n, 1, a2, n, ipiv, x, n, error )

    err = maxval( x - b )
    if ( err>1.d-12 ) then
        call abort(96, "LAPACK unit test failed" )
    endif

    if (root) write(*,'(" Done. err=",es15.8)') err
    if (root) write(*,*) "--------------------------------"
end subroutine lapack_unit_test


!-------------------------------------------------------------------------
! dump runtime backup for the solid solver
!-------------------------------------------------------------------------
subroutine dump_solid_backup( time )
    implicit none

    real(kind=rk), intent(in) :: time
    character(len=clong) :: filename
    integer :: i

    !-- only root rank dumps backup
    if (root) then

        write( filename,'(a, "_", a, ".dat")') "solid", trim(adjustl(timestr(time)))

        write(*,'(A)',advance='no') "Backuping solid solver..."
        write(*,'(A)',advance='no') "file="//filename
        write(*,'(" time=",e11.4)',advance='no') time


        open(14,file=filename,status='replace',form='formatted')


        write(14,*) time
        write(14,*) ns, nBeams

        do i = 1, nBeams
            write(14,*) beams(i)%x, beams(i)%y, beams(i)%vx, beams(i)%vy
            write(14,*) beams(i)%theta, beams(i)%theta_dot, beams(i)%ax, beams(i)%ay

            write(14,*) beams(i)%T
            write(14,*) beams(i)%theta_old, beams(i)%theta_oldold
            write(14,*) beams(i)%theta_dot_old, beams(i)%theta_dot_oldold

            write(14,*) beams(i)%pressure_old, beams(i)%pressure_new
            write(14,*) beams(i)%tau_old, beams(i)%tau_new

            write(14,*) beams(i)%dt_old
            write(14,*) beams(i)%StartupStep
        enddo

        close(14)
        write(*,'(A)',advance='yes') "...DONE!"
    endif
end subroutine dump_solid_backup



!-------------------------------------------------------------------------
! read runtime backup for the solid solver
!-------------------------------------------------------------------------
subroutine read_solid_backup( filename )
    implicit none

    integer :: ns_file, nBeams_file, i, mpicode
    character(len=*), intent(in) :: filename
    real(kind=rk) :: time

    ! only root shall read in the file, the data are then sent to the other ranks
    if (root) then
        write(*,'(A)') "Reading in backup of solid solver: "//filename

        open(14,file=filename,status='old',form='formatted',action='read')
        read(14,*) time
        read(14,*) ns_file, nBeams_file

        write(*,'("(time=",e11.4,")")') time

        if ((ns_file/=ns).or.(nBeams_file/=nBeams)) then
            call abort(97,"Cant resume solid backup: resolution ns or beam number doesnt match")
        endif

        do i =1, nBeams
            read(14,*) beams(i)%x, beams(i)%y, beams(i)%vx, beams(i)%vy
            read(14,*) beams(i)%theta, beams(i)%theta_dot, beams(i)%ax, beams(i)%ay

            read(14,*) beams(i)%T
            read(14,*) beams(i)%theta_old, beams(i)%theta_oldold
            read(14,*) beams(i)%theta_dot_old, beams(i)%theta_dot_oldold

            read(14,*) beams(i)%pressure_old, beams(i)%pressure_new
            read(14,*) beams(i)%tau_old, beams(i)%tau_new
            read(14,*) beams(i)%dt_old
            read(14,*) beams(i)%StartupStep
        enddo
        close(14)

        write(*,'(A)') "...done reading solid backup."
    endif

    do i =1, nBeams
        call MPI_Bcast( beams(i)%x, size(beams(i)%x,1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
        call MPI_Bcast( beams(i)%y, size(beams(i)%y,1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
        call MPI_Bcast( beams(i)%vx, size(beams(i)%vx,1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
        call MPI_Bcast( beams(i)%vy, size(beams(i)%vy,1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
        call MPI_Bcast( beams(i)%theta, size(beams(i)%theta,1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
        call MPI_Bcast( beams(i)%theta_dot, size(beams(i)%theta_dot,1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )

        call MPI_Bcast( beams(i)%T, size(beams(i)%T,1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
        call MPI_Bcast( beams(i)%theta_old, size(beams(i)%theta_old,1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
        call MPI_Bcast( beams(i)%theta_oldold, size(beams(i)%theta_oldold,1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
        call MPI_Bcast( beams(i)%theta_dot_old, size(beams(i)%theta_dot_old,1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
        call MPI_Bcast( beams(i)%theta_dot_oldold, size(beams(i)%theta_dot_oldold,1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )

        call MPI_Bcast( beams(i)%ax, size(beams(i)%ax,1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
        call MPI_Bcast( beams(i)%ay, size(beams(i)%ay,1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
        call MPI_Bcast( beams(i)%pressure_old, size(beams(i)%pressure_old,1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
        call MPI_Bcast( beams(i)%pressure_new, size(beams(i)%pressure_new,1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
        call MPI_Bcast( beams(i)%tau_old, size(beams(i)%tau_old,1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
        call MPI_Bcast( beams(i)%tau_new, size(beams(i)%tau_new,1), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )

        call MPI_Bcast( beams(i)%dt_old, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, mpicode )
        call MPI_Bcast( beams(i)%StartupStep, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, mpicode )
    enddo

end subroutine read_solid_backup
