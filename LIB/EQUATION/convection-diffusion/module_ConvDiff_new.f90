!> \brief Implementation of 3d/2d Convection Diffusion Equations
! ********************************************************************************************
module module_convdiff_new

  use module_globals
  ! ini file parser module, used to read parameters. note: in principle, you can also
  ! just use any reader you feel comfortable with, as long as you can read the parameters
  ! from a file.
  use module_ini_files_parser_mpi
  use mpi
  use module_t_files

  implicit none

  ! I usually find it helpful to use the private keyword by itself initially, which specifies
  ! that everything within the module is private unless explicitly marked public.
  PRIVATE

  !**********************************************************************************************
  ! These are the important routines that are visible to WABBIT:
  !**********************************************************************************************
  PUBLIC :: READ_PARAMETERS_convdiff, PREPARE_SAVE_DATA_convdiff, RHS_convdiff, GET_DT_BLOCK_convdiff, &
  INICOND_convdiff, FIELD_NAMES_convdiff, STATISTICS_convdiff, TIME_STATISTICS_convdiff
  !**********************************************************************************************

  ! user defined data structure for time independent parameters, settings, constants
  ! and the like. only visible here.
  type :: type_paramsb
    real(kind=rk) :: CFL, T_end, T_swirl, CFL_nu=0.094, u_const=0.0_rk, gamma, tau
    real(kind=rk) :: domain_size(3)=0.0_rk, scalar_integral=0.0_rk,w0(3)=0.0_rk, scalar_max=0.0_rk
    real(kind=rk), allocatable, dimension(:) :: nu, u0x,u0y,u0z,phi_boundary
    real(kind=rk), allocatable, dimension(:,:) :: blob_width,x0,y0,z0
    integer(kind=ik) :: dim, N_scalars, N_fields_saved, Nblobs
    logical :: time_statistics = .false.
    integer(kind=ik) :: N_time_statistics = 0
    character(len=cshort), allocatable :: names(:), inicond(:), velocity(:), time_statistics_names(:)
    character(len=cshort) :: discretization,boundary_type
    logical,dimension(3):: periodic_BC=(/.true.,.true.,.true./)
  end type type_paramsb

  ! parameters for this module. they should not be seen outside this physics module
  ! in the rest of the code. WABBIT does not need to know them.
  type(type_paramsb), save :: params_convdiff



contains

#include "statistics_convdiff.f90"
#include "time_statistics_convdiff.f90"
#include "rhs_convdiff.f90"

  !-----------------------------------------------------------------------------
  ! main level wrapper routine to read parameters in the physics module. It reads
  ! from the same ini file as wabbit, and it reads all it has to know. note in physics modules
  ! the parameter struct for wabbit is not available.
  subroutine READ_PARAMETERS_convdiff( filename, g )
    implicit none

    character(len=*), intent(in) :: filename
    integer(kind=ik), intent(in) :: g

    real(kind=rk), dimension(3)      :: domain_size=0.0_rk
    real(kind=rk), allocatable, dimension(:,:) :: matr_dummy
    ! inifile structure
    type(inifile) :: FILE
    ! read the file, only process 0 should create output on screen
    call set_lattice_spacing_mpi(1.0d0)
    call read_ini_file_mpi(FILE, filename, .true.)

    call read_param_mpi(FILE, 'ConvectionDiffusion', 'N_scalars', params_convdiff%N_scalars, 1  )
    call read_param_mpi(FILE, 'ConvectionDiffusion', 'Nblobs', params_convdiff%Nblobs, 1 )

    ! nu, and u0 are used for the velocity of the convection, they are per scalar
    allocate( params_convdiff%nu(1:params_convdiff%N_scalars))
    allocate( params_convdiff%u0x(1:params_convdiff%N_scalars))
    allocate( params_convdiff%u0y(1:params_convdiff%N_scalars))
    allocate( params_convdiff%u0z(1:params_convdiff%N_scalars))

    ! x0, y0 and z0 and blob width are used for the initial conditions, they are per scalar and possibly blob in a matrix
    ! as a matrix, they have to be in form like the butcher tableau
    allocate( params_convdiff%x0(1:params_convdiff%N_scalars, 1:params_convdiff%Nblobs))
    allocate( params_convdiff%y0(1:params_convdiff%N_scalars, 1:params_convdiff%Nblobs))
    allocate( params_convdiff%z0(1:params_convdiff%N_scalars, 1:params_convdiff%Nblobs))
    allocate( params_convdiff%blob_width(1:params_convdiff%N_scalars, 1:params_convdiff%Nblobs) )
    allocate( matr_dummy(1:params_convdiff%N_scalars, 1:params_convdiff%Nblobs) )
    matr_dummy = 0.0_rk

    ! following are string names for which condition to choose
    allocate( params_convdiff%inicond(1:params_convdiff%N_scalars))
    allocate( params_convdiff%velocity(1:params_convdiff%N_scalars))



    call read_param_mpi(FILE, 'ConvectionDiffusion', 'nu', params_convdiff%nu )
    call read_param_mpi(FILE, 'ConvectionDiffusion', 'u0x', params_convdiff%u0x )
    call read_param_mpi(FILE, 'ConvectionDiffusion', 'u0y', params_convdiff%u0y )
    call read_param_mpi(FILE, 'ConvectionDiffusion', 'u0z', params_convdiff%u0z )

    call read_param_mpi(FILE, 'ConvectionDiffusion', 'gamma', params_convdiff%gamma, 0.0_rk ) ! reaction constant
    call read_param_mpi(FILE, 'ConvectionDiffusion', 'tau', params_convdiff%tau , 0.0_rk)

    call read_param_mpi(FILE, 'ConvectionDiffusion', 'x0', params_convdiff%x0, matr_dummy )
    call read_param_mpi(FILE, 'ConvectionDiffusion', 'y0', params_convdiff%y0, matr_dummy )
    call read_param_mpi(FILE, 'ConvectionDiffusion', 'z0', params_convdiff%z0, matr_dummy )
    matr_dummy = 1.0_rk  ! avoid 0 division by default for blob_width
    call read_param_mpi(FILE, 'ConvectionDiffusion', 'blob_width', params_convdiff%blob_width, matr_dummy )

    call read_param_mpi(FILE, 'ConvectionDiffusion', 'inicond', params_convdiff%inicond, (/'blob'/) )
    call read_param_mpi(FILE, 'ConvectionDiffusion', 'velocity', params_convdiff%velocity, (/'constant'/) )

    ! allocate( params_convdiff%blob_width(1:params_convdiff%N_scalars))
    ! call read_param_mpi(FILE, 'ConvectionDiffusion', 'Nblobs', params_convdiff%Nblobs, 1 )
    ! allocate( params_convdiff%blobs_width(1:params_convdiff%Nblobs) )
    ! allocate( params_convdiff%blobs_x0(1:params_convdiff%Nblobs) )
    ! allocate( params_convdiff%blobs_y0(1:params_convdiff%Nblobs) )
    ! allocate( params_convdiff%blobs_z0(1:params_convdiff%Nblobs) )
    ! call read_param_mpi(FILE, 'ConvectionDiffusion', 'blobs_width', params_convdiff%blobs_width )
    ! call read_param_mpi(FILE, 'ConvectionDiffusion', 'x0', params_convdiff%blobs_x0 )
    ! call read_param_mpi(FILE, 'ConvectionDiffusion', 'y0', params_convdiff%blobs_y0 )
    ! call read_param_mpi(FILE, 'ConvectionDiffusion', 'z0', params_convdiff%blobs_z0 )

    ! time statistics (averaging or similar)
    call read_param_mpi(FILE, 'Time-Statistics', 'time_statistics', params_convdiff%time_statistics, .false.)
    if (params_convdiff%time_statistics) then
        call read_param_mpi(FILE, 'Time-Statistics', 'N_time_statistics', params_convdiff%N_time_statistics, 1)
        allocate( params_convdiff%time_statistics_names(1:params_convdiff%N_time_statistics) )
        call read_param_mpi(FILE, 'Time-Statistics', 'time_statistics_names', params_convdiff%time_statistics_names, (/"none"/))
    endif

    call read_param_mpi(FILE, 'ConvectionDiffusion', 'u_const', params_convdiff%u_const, 0.0_rk )


    call read_param_mpi(FILE, 'Domain', 'dim', params_convdiff%dim, 2 )
    call read_param_mpi(FILE, 'ConvectionDiffusion', 'w0', params_convdiff%w0(1:params_convdiff%dim) )
    call read_param_mpi(FILE, 'Domain', 'domain_size', params_convdiff%domain_size(1:params_convdiff%dim) )
    call read_param_mpi(FILE, 'Domain', 'periodic_BC', params_convdiff%periodic_BC(1:params_convdiff%dim), &
                                                       params_convdiff%periodic_BC(1:params_convdiff%dim) )
    if ( .not. All(params_convdiff%periodic_BC) ) then
      allocate( params_convdiff%phi_boundary(1:params_convdiff%N_scalars))
      call read_param_mpi(FILE, 'ConvectionDiffusion', 'boundary_type', params_convdiff%boundary_type,'--' )
      call read_param_mpi(FILE, 'ConvectionDiffusion', 'phi_boundary', params_convdiff%phi_boundary,  params_convdiff%phi_boundary )
    endif


    call read_param_mpi(FILE, 'Discretization', 'order_discretization', params_convdiff%discretization, "FD_2nd_central")

    if ( params_convdiff%discretization=='FD_4th_central_optimized' .and. g<=2  ) then
      call abort(91020181, "Number of ghost nodes for this scheme is 3!")
    endif

    call read_param_mpi(FILE, 'Saving', 'N_fields_saved', params_convdiff%N_fields_saved, 1 )
    allocate( params_convdiff%names(1:params_convdiff%N_fields_saved))
    call read_param_mpi(FILE, 'Saving', 'field_names', params_convdiff%names, (/"phi1"/) )


    call read_param_mpi(FILE, 'Time', 'CFL', params_convdiff%CFL, 1.0_rk)
    call read_param_mpi(FILE, 'Time', 'CFL_nu', params_convdiff%CFL_nu, 0.99_rk*2.79_rk/(dble(params_convdiff%dim)*pi**2) )
    call read_param_mpi(FILE, 'Time', 'time_max', params_convdiff%T_end, 1.0_rk)
    call read_param_mpi(FILE, 'ConvectionDiffusion', 'T_swirl', params_convdiff%T_swirl, params_convdiff%T_end)

    call clean_ini_file_mpi( FILE )

    if ( params_convdiff%dim == 2) then
        params_convdiff%u0z=0.0_rk
        params_convdiff%z0=0.0_rk
    endif
  end subroutine READ_PARAMETERS_convdiff


  !-----------------------------------------------------------------------------
  ! save data. Since you might want to save derived data, such as the vorticity,
  ! the divergence etc., which are not in your state vector, this routine has to
  ! copy and compute what you want to save to the work array.
  !
  ! In the main code, save_fields than saves the first N_fields_saved components of the
  ! work array to file.
  !
  ! NOTE that as we have way more work arrays than actual state variables (typically
  ! for a RK4 that would be >= 4*dim), you can compute a lot of stuff, if you want to.
  !-----------------------------------------------------------------------------
  subroutine PREPARE_SAVE_DATA_convdiff( time, u, g, x0, dx, work )
    implicit none
    ! it may happen that some source terms have an explicit time-dependency
    ! therefore the general call has to pass time
    real(kind=rk), intent (in) :: time

    ! block data, containg the state vector. In general a 4D field (3 dims+components)
    ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
    real(kind=rk), intent(in) :: u(1:,1:,1:,1:)

    ! as you are allowed to compute the RHS only in the interior of the field
    ! you also need to know where 'interior' starts: so we pass the number of ghost points
    integer, intent(in) :: g

    ! for each block, you'll need to know where it lies in physical space. The first
    ! non-ghost point has the coordinate x0, from then on its just cartesian with dx spacing
    real(kind=rk), intent(in) :: x0(1:3), dx(1:3)

    ! output in work array.
    real(kind=rk), intent(inout) :: work(1:,1:,1:,1:)

    ! local variables
    integer(kind=ik) :: neqn, nwork, i_var, k, i_time_statistics
    integer(kind=ik), dimension(3) :: Bs
    character(len=cshort) :: name

    Bs(1) = size(u,1) - 2*g
    Bs(2) = size(u,2) - 2*g
    Bs(3) = size(u,3) - 2*g

    do k = 1, size(params_convdiff%names,1)
        name = params_convdiff%names(k)
        
        ! saving of phi, read if it starts with "phi" and has no "-" in there for time statistics
        if (name(1:3) == "phi" .and. index(name, "-") == 0) then
            ! check if the 4th character is actually a digit
            if (.not. (name(4:4) >= '0' .and. name(4:4) <= '9')) then
                call abort(250929, "CONVDIFF: field_name '"//trim(name)//"' has invalid format. Expected 'phiX' where X is a digit (1-9), but found '"//name(4:4)//"' at position 4.")
            end if
            read( name(4:4), * ) i_var
            if (i_var >=1 .and. i_var <= params_convdiff%N_scalars) then
                work(:,:,:,k) = u(:,:,:,i_var)
            else
                call abort(250922, "CONVDIFF: You try to save a variable that does not exist")
            endif
        endif

        ! saving of velocity field
        if (name(1:2) == 'ux' .or. name(1:2) == 'uy' .or. name(1:2) == 'uz') then
            read( name(4:4), * ) i_var
            if (params_convdiff%dim == 2) then
                if (size(work,4) - k < 1) then
                    call abort(250922,"CONVDIFF: Not enough space to compute velocity")
                endif
                if (k > params_convdiff%N_scalars+2) then
                    call abort(250922, "CONVDIFF: You try to save a variable that does not exist")
                endif
                call create_velocity_field_2d( time, g, Bs, dx, x0, work(:,:,:,k), i_var, u(:,:,1,i_var) )
            elseif (params_convdiff%dim == 3) then
                if (size(work,4) - k < 2) then
                    call abort(250922,"CONVDIFF: Not enough space to compute velocity")
                endif
                if (k > params_convdiff%N_scalars+3) then
                    call abort(250922, "CONVDIFF: You try to save a variable that does not exist " // trim(name) )
                endif
                call create_velocity_field_3d( time, g, Bs, dx, x0, work(:,:,:,k:k+2), i_var )
            endif
            if (name(1:2) == 'uy') then
                work(:,:,:,k) = work(:,:,:,k+1)
            elseif (name(1:2) == 'uz') then
                if (params_convdiff%dim == 2) then
                    call abort(250922, "CONVDIFF: You try to save uz in 2d!")
                endif
                work(:,:,:,k) = work(:,:,:,k+2)
            endif
        endif

        ! if any of those endings is in the name, then it is a timestatistics variable
        if (index(name, "-avg") > 0 .or. index(name, "-mean") > 0 .or. index(name, "-var") > 0 &
            .or. index(name, "-minmax") > 0 .or. index(name, "-min") > 0 .or. index(name, "-max") > 0) then
            ! now we have to find the index of it
            do i_time_statistics = 1, params_convdiff%N_time_statistics
                if (name == trim(params_convdiff%time_statistics_names(i_time_statistics))) exit
            end do

            ! safety check if the name was found
            if (i_time_statistics > params_convdiff%N_time_statistics) then
                call abort(250929, "CONV_DIFF: field_name " // trim(name) // " not found in time_statistics_names")
            endif

            work(:,:,:,k) = u(:,:,:,params_convdiff%N_scalars + i_time_statistics)
        endif
    enddo

  end subroutine


  !-----------------------------------------------------------------------------
  ! when savig to disk, WABBIT would like to know how you named you variables.
  ! e.g. u(:,:,:,1) is called "ux"
  !
  ! the main routine save_fields has to know how you label the stuff you want to
  ! store from the work array, and this routine returns those strings
  !-----------------------------------------------------------------------------
  subroutine FIELD_NAMES_convdiff( N, name )
    implicit none
    ! component index
    integer(kind=ik), intent(in) :: N
    ! returns the name
    character(len=cshort), intent(out) :: name

    if (allocated(params_convdiff%names)) then
      name = params_convdiff%names(N)
    else
      call abort(5554,'Something ricked')
    endif

  end subroutine FIELD_NAMES_convdiff



  !-----------------------------------------------------------------------------
  ! subroutine statistics_convdiff()
  !   implicit none
  ! end subroutine


  !-----------------------------------------------------------------------------
  ! setting the time step is very physics-dependent. Sometimes you have a CFL like
  ! condition, sometimes not. So each physic module must be able to decide on its
  ! time step. This routine is called for all blocks, the smallest returned dt is used.
  !-----------------------------------------------------------------------------
  subroutine GET_DT_BLOCK_convdiff( time, u, Bs, g, x0, dx, dt )
    implicit none

    ! it may happen that some source terms have an explicit time-dependency
    ! therefore the general call has to pass time
    real(kind=rk), intent (in) :: time

    ! block data, containg the state vector. In general a 4D field (3 dims+components)
    ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
    real(kind=rk), intent(in) :: u(1:,1:,1:,1:)

    ! as you are allowed to compute the RHS only in the interior of the field
    ! you also need to know where 'interior' starts: so we pass the number of ghost points
    integer, intent(in) :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs

    ! for each block, you'll need to know where it lies in physical space. The first
    ! non-ghost point has the coordinate x0, from then on its just cartesian with dx spacing
    real(kind=rk), intent(in) :: x0(1:3), dx(1:3)

    ! the dt for this block is returned to the caller:
    real(kind=rk), intent(out) :: dt

    ! TODO: make this global and allocatable
    real(kind=rk) :: u0(1:Bs(1)+2*g, 1:Bs(2)+2*g, 1:Bs(3)+2*g, 1:3)
    integer(kind=ik) :: i, ix, iy
    real(kind=rk) :: x,y,unorm,c_react

    dt = 9.9e9_rk


    do i = 1, params_convdiff%N_scalars
        if (params_convdiff%dim == 2) then
            call create_velocity_field_2d( time, g, Bs, dx, x0, u0(:,:,1,1:2), i, u(:,:,1,1) )
            unorm = maxval( u0(:,:,1,1)*u0(:,:,1,1) + u0(:,:,1,2)*u0(:,:,1,2) )
        else
            call create_velocity_field_3d( time, g, Bs, dx, x0, u0, i )
            unorm = maxval( u0(:,:,:,1)*u0(:,:,:,1) + u0(:,:,:,2)*u0(:,:,:,2) + u0(:,:,:,3)*u0(:,:,:,3) )
        endif

        if ( unorm < 1.0e-8_rk) then
            ! if the value of u is very small, which may happen if it is time dependent
            ! we choose some fixed value in order not to miss the instant when u becomes
            ! large again.

            dt = min( 1.0e-3_rk, dt )

            ! constant zero velocity: nothing from CFL
            if (params_convdiff%velocity(i) == "constant") dt = 9e9_rk
        else
            dt = min(params_convdiff%CFL * dx(1) / sqrt(unorm), dt)
        endif

        if (params_convdiff%nu(i) > 1.0e-13_rk) then
          dt = min(dt,  params_convdiff%CFL_nu * minval(dx(1:params_convdiff%dim))**2 / params_convdiff%nu(i))
          if (params_convdiff%gamma>1.0e-13_rk) then
            c_react = 0.15*params_convdiff%gamma
            dt = min(dt,params_convdiff%CFL_nu/c_react * minval(dx(1:params_convdiff%dim)) )
          endif
        endif

    enddo


  end subroutine GET_DT_BLOCK_convdiff

  !-----------------------------------------------------------------------------
  ! main level wrapper for setting the initial condition on a block
  !-----------------------------------------------------------------------------
  subroutine INICOND_convdiff( time, u, g, x0, dx )
    implicit none

    ! it may happen that some source terms have an explicit time-dependency
    ! therefore the general call has to pass time
    real(kind=rk), intent (in) :: time

    ! block data, containg the state vector. In general a 4D field (3 dims+components)
    ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
    real(kind=rk), intent(inout) :: u(1:,1:,1:,1:)

    ! as you are allowed to compute the RHS only in the interior of the field
    ! you also need to know where 'interior' starts: so we pass the number of ghost points
    integer, intent(in) :: g

    ! for each block, you'll need to know where it lies in physical space. The first
    ! non-ghost point has the coordinate x0, from then on its just cartesian with dx spacing
    real(kind=rk), intent(in) :: x0(1:3), dx(1:3)

    integer(kind=ik) :: ix, iy, iz, i, iblob, bx, by, bz
    integer(kind=ik), dimension(3) :: Bs
    real(kind=rk) :: x, y, c0x, c0y, z, c0z, lambd, delta

    ! compute the size of blocks
    Bs(1) = size(u,1) - 2*g
    Bs(2) = size(u,2) - 2*g
    Bs(3) = size(u,3) - 2*g

    u = 0.0_rk

    do i = 1, params_convdiff%N_scalars
    c0x = params_convdiff%x0(i,1)
    c0y = params_convdiff%y0(i,1)
    c0z = params_convdiff%z0(i,1)

    select case (params_convdiff%inicond(i))
        case ("noise")
            do ix = 1, Bs(1)+2*g
                do iy = 1, Bs(2)+2*g
                    u(ix,iy,:,i) = rand_nbr()
                end do
            end do

        case ("zero")
            u(:,:,:,i) = 0.0_rk

        case ("const")
            u(:,:,:,i) = params_convdiff%blob_width(i,1)

        case ("sin")
            if (params_convdiff%dim==2) then
                do ix = 1, Bs(1)+2*g
                    do iy = 1, Bs(2)+2*g
                        ! compute x,y coordinates from spacing and origin
                        x = dble(ix-(g+1)) * dx(1) + x0(1)
                        y = dble(iy-(g+1)) * dx(2) + x0(2)

                        u(ix,iy,:,i) = sin(2.0_rk*pi*x)
                    end do
                end do
            else
                call abort(66273,"this inicond is 2d only..")
            endif

        case ("sin+1")
            if (params_convdiff%dim==2) then
                do ix = 1, Bs(1)+2*g
                    do iy = 1, Bs(2)+2*g
                        ! compute x,y coordinates from spacing and origin
                        x = dble(ix-(g+1)) * dx(1) + x0(1)
                        y = dble(iy-(g+1)) * dx(2) + x0(2)

                        u(ix,iy,:,i) = sin(2.0_rk*pi*x) + 1.0_rk
                    end do
                end do
            else
                call abort(66273,"this inicond is 2d only..")
            endif

        case ("cyclogenesis")
            if (params_convdiff%dim==2) then
                do ix = 1, Bs(1)+2*g
                    do iy = 1, Bs(2)+2*g
                        ! compute x,y coordinates from spacing and origin
                        x = dble(ix-(g+1)) * dx(1) + x0(1) - c0x
                        y = dble(iy-(g+1)) * dx(2) + x0(2) - c0y

                        u(ix,iy,:,i) = -tanh( y / params_convdiff%blob_width(i,1) )
                    end do
                end do
            else
                call abort(66273,"this inicond is 2d only..")
            endif

        case ("circle")
            if (params_convdiff%dim==2) then
                lambd = 0.005 *maxval(params_convdiff%domain_size)
                do ix = 1, Bs(1)+2*g
                    do iy = 1, Bs(2)+2*g
                        ! compute x,y coordinates from spacing and origin
                        x = dble(ix-(g+1)) * dx(1) + x0(1) - c0x
                        y = dble(iy-(g+1)) * dx(2) + x0(2) - c0y
                        !u(ix,iy,:,i) =  1/(1+dexp( (dsqrt(x**2 + y**2) - params_convdiff%blob_width(i,1) ) /lambd))
                        u(ix,iy,:,i) = 0.5*(1-dtanh( (dsqrt(x**2 + y**2) - params_convdiff%blob_width(i,1) ) /lambd))
                    end do
                end do
            else
                call abort(66273,"this inicond is 2d only..")
            endif
        case ("1D-bump")
            if (params_convdiff%dim==2) then
                delta = dsqrt(8/params_convdiff%gamma)
                do ix = 1, Bs(1)+2*g
                    do iy = 1, Bs(2)+2*g
                        ! compute x,y coordinates from spacing and origin
                        x = dble(ix-(g+1)) * dx(1) + x0(1) - 0.5_rk * params_convdiff%domain_size(1)
                        y = dble(iy-(g+1)) * dx(2) + x0(2)
                        !u(ix,iy,:,i) =  1/(1+dexp( (dsqrt(x**2 + y**2) - params_convdiff%blob_width(i) ) /lambd))
                        u(ix,iy,:,i) = 0.5*(1-dtanh( (abs(x) - params_convdiff%blob_width(i,1) ) / delta))
                    end do
                end do
            else
                call abort(66273,"this inicond is 2d only..")
            endif


        case ("blob")
            u = 0.0_rk
            ! loop over all points and compute x,y,z coordinates from spacing and origin
            do iz = merge(1, g+1, params_convdiff%dim==2), merge(1, Bs(3)+g, params_convdiff%dim==2)
                z = 0.0_rk
                if (params_convdiff%dim == 3) z = dble(iz-(g+1)) * dx(3) + x0(3)
                do iy = 1, Bs(2)+2*g
                    y = dble(iy-(g+1)) * dx(2) + x0(2)
                    do ix = 1, Bs(1)+2*g
                        x = dble(ix-(g+1)) * dx(1) + x0(1)

                        ! set actual inicond gauss blob for each blob
                        do iblob = 1, params_convdiff%Nblobs
                            ! if a blob is close to the domain border but has a large width, periodicity effects may impact the flow on the other periodic side
                            ! in order to incorporate this, we add the periodic effect by adding blobs shifted by domain size
                            ! as this initial condition is only executed once, this "more expensive" but more exact approach is chosen
                            do bx = merge(-1, 0, params_convdiff%periodic_BC(1)), merge(1, 0, params_convdiff%periodic_BC(1))
                                do by = merge(-1, 0, params_convdiff%periodic_BC(2)), merge(1, 0, params_convdiff%periodic_BC(2))
                                    do bz = merge(-1, 0, params_convdiff%periodic_BC(3) .and. params_convdiff%dim == 3), &
                                            merge( 1, 0, params_convdiff%periodic_BC(3) .and. params_convdiff%dim == 3)
                                        u(ix,iy,iz,i) = u(ix,iy,iz,i) + dexp( -((x-params_convdiff%x0(i, iblob)+bx*params_convdiff%domain_size(1))**2 &
                                            + (y-params_convdiff%y0(i, iblob)+by*params_convdiff%domain_size(2))**2 &
                                            + (z-params_convdiff%z0(i, iblob)+bz*params_convdiff%domain_size(3))**2) / params_convdiff%blob_width(i,iblob) )
                                    enddo
                                enddo
                            enddo
                        enddo
                    end do
                end do
            enddo
        case default
            call abort(72637,"Error. Inital conditon for conv-diff is unkown: "//trim(adjustl(params_convdiff%inicond(i))))
        end select

    enddo


  end subroutine INICOND_convdiff



end module module_convdiff_new
