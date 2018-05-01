!> \dir
!> \brief
!>Implementation of 3d/2d Navier Stokes Physics

!-----------------------------------------------------------------
!> \file
!> \brief
!! Module of public 2D/3D Navier Stokes equation
!> \details
!!    * reads in params
!!    * sets initial conditions
!!    * calls RHS
!!    * calculates time step
!!
!> \version 23.1.2018
!> \author P.Krah
!-----------------------------------------------------------------

module module_navier_stokes_new

  !---------------------------------------------------------------------------------------------
  ! modules
  use module_navier_stokes_params
  use module_operators, only : compute_vorticity
  use module_ns_penalization

  implicit none

  ! I usually find it helpful to use the private keyword by itself initially, which specifies
  ! that everything within the module is private unless explicitly marked public.
  PRIVATE

  !**********************************************************************************************
  ! These are the important routines that are visible to WABBIT:
  !**********************************************************************************************
  PUBLIC :: READ_PARAMETERS_NSTOKES, PREPARE_SAVE_DATA_NSTOKES, RHS_NSTOKES, GET_DT_BLOCK_NSTOKES, &
            CONVERT_STATEVECTOR2D, PACK_STATEVECTOR2D, INICOND_NSTOKES, FIELD_NAMES_NStokes,&
            STATISTICS_NStokes,FILTER_NSTOKES
  !**********************************************************************************************
  ! parameters for this module. they should not be seen outside this physics module
  ! in the rest of the code. WABBIT does not need to know them.
  type(type_params_ns),save :: params_ns


  ! real(kind=rk)   ,parameter :: rho_init=1.645_rk,p_init=101330.0_rk,u_init=36.4_rk,v0_=0.0_rk,T_init=200!273.15_rk
  real(kind=rk)       ,save :: rho_init=1_rk,p_init=1.0_rk,T_init=1!273.15_rk
  real(kind=rk)       ,save :: u_init(3)=(/1.0_rk,1.0_rk,0.0_rk/)

 
contains


  include "RHS_2D_navier_stokes.f90"
  include "RHS_3D_navier_stokes.f90"
!-----------------------------------------------------------------------------
  !> \brief Reads in parameters of physics module
  !> \details
  !> Main level wrapper routine to read parameters in the physics module. It reads
  !> from the same ini file as wabbit, and it reads all it has to know. note in physics modules
  !> the parameter struct for wabbit is not available.
  subroutine READ_PARAMETERS_NStokes( filename )
    implicit none
    !> name of inifile
    character(len=*), intent(in) :: filename

    ! inifile structure
    type(inifile)               :: FILE
    integer(kind=ik)            :: dF
    integer(kind=ik)           :: mpicode
   


    ! ==================================================================
    ! initialize MPI parameter
    ! ------------------------------------------------------------------
    ! we still need to know about mpirank and mpisize, occasionally
    call MPI_COMM_SIZE (WABBIT_COMM, params_ns%mpisize, mpicode)
    call MPI_COMM_RANK (WABBIT_COMM, params_ns%mpirank, mpicode)

    if (params_ns%mpirank==0) then
      write(*,*)
      write(*,*)
      write(*,'(80("<"))')
      write(*,*) "Initializing Navier Stokes module!"
      write(*,'(80("<"))')
      write(*,*)
      write(*,*)
    endif


    ! read the file, only process 0 should create output on screen
    call set_lattice_spacing_mpi(1.0d0)
    ! open file
    call read_ini_file_mpi(FILE, filename, .true.)
    ! init all parameters used in ns_equations
    call init_navier_stokes_eq(params_ns, FILE)
    ! init all parameters used for penalization
    call init_penalization(    params_ns, FILE)
    ! init all parameters used for the filter
    call init_filter(   params_ns%filter, FILE)
    ! init all params for organisation 
    call init_other_params(params_ns,     FILE )

    ! read in initial conditions
    if (params_ns%mpirank==0) then
      write(*,*)
      write(*,*)
      write(*,*) "PARAMS: initial conditions"
      write(*,'(" ---------------------------")')
    endif
    call read_param_mpi(FILE, 'Navier_Stokes', 'inicond'      , params_ns%inicond, "pressure_blob" )
    call read_param_mpi(FILE, 'Navier_Stokes', 'inicond_width', params_ns%inicond_width, params_ns%Lx*0.1_rk )
    call read_param_mpi(FILE, 'Navier_Stokes', 'initial_pressure' , p_init, p_init )
    call read_param_mpi(FILE, 'Navier_Stokes', 'initial_velocity' , u_init, u_init )
    call read_param_mpi(FILE, 'Navier_Stokes', 'initial_temperature', T_init, T_init )
    call read_param_mpi(FILE, 'Navier_Stokes', 'initial_density', rho_init, rho_init )

    ! set global parameters pF,rohF, UxF etc
    UzF=-1
    do dF = 1, params_ns%number_data_fields
                if ( params_ns%names(dF) == "p" ) pF = dF
                if ( params_ns%names(dF) == "rho" ) rhoF = dF
                if ( params_ns%names(dF) == "Ux" ) UxF = dF
                if ( params_ns%names(dF) == "Uy" ) UyF = dF
                if ( params_ns%names(dF) == "Uz" ) UzF = dF
    end do

    call clean_ini_file_mpi( FILE )

  end subroutine READ_PARAMETERS_NStokes

 




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
  subroutine PREPARE_SAVE_DATA_NStokes( time, u, g, x0, dx, work )
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

    ! output in work array.
    real(kind=rk), allocatable :: vort(:,:,:,:)

    ! local variables
    integer(kind=ik) ::  Bs

    Bs = size(u,1)-2*g



    !
    ! compute vorticity
    allocate(vort(size(u,1),size(u,2),size(u,3),3))

    if (size(u,3)==1) then
          ! ---------------------------------
          ! save all datafields in u
          work(:,:,:,rhoF) = u(:,:,:,1)**2
          work(:,:,:,UxF)  = u(:,:,:,2)/u(:,:,:,1)
          work(:,:,:,UyF)  = u(:,:,:,3)/u(:,:,:,1)
          work(:,:,:,pF)   = u(:,:,:,4)
          ! ---------------------------------
          
        ! only wx,wy (2D - case)
        call compute_vorticity(  u(:,:,:,UxF)/u(:,:,:,rhoF), &
                                 u(:,:,:,UyF)/u(:,:,:,rhoF), &
                                 0*u(:,:,:,rhoF), &
                                 dx, Bs, g, params_ns%discretization, &
                                 vort)

        work(:,:,:,5)=vort(:,:,:,1)

        !write out mask
        if (params_ns%penalization) then
          call get_mask(work(:,:,1,6), x0, dx, Bs, g )
        else
          work(:,:,1,6)=0.0_rk
        endif

    else
      call abort(564567,"Error: [module_navier_stokes.f90] 3D case not implemented")
        ! wx,wy,wz (3D - case)
   !      call compute_vorticity(  u(:,:,:,UxF)/u(:,:,:,rhoF), &
   !                               u(:,:,:,UyF)/u(:,:,:,rhoF), &
   !                               u(:,:,:,UzF)/u(:,:,:,rhoF), &
   !                               dx, Bs, g, params_ns%discretization, &
   !                               vort)
   ! !     work(:,:,:,=vort(:,:,:,1:3)
    endif
    deallocate(vort)

    ! mask

  end subroutine


  !-----------------------------------------------------------------------------
  ! when savig to disk, WABBIT would like to know how you named you variables.
  ! e.g. u(:,:,:,1) is called "ux"
  !
  ! the main routine save_fields has to know how you label the stuff you want to
  ! store from the work array, and this routine returns those strings
  !-----------------------------------------------------------------------------
  subroutine FIELD_NAMES_NStokes( N, name )
    implicit none
    ! component index
    integer(kind=ik), intent(in) :: N
    ! returns the name
    character(len=80), intent(out) :: name

    if (allocated(params_ns%names)) then
      name = params_ns%names(N)
    else
      call abort(5554,'Something ricked')
    endif

  end subroutine FIELD_NAMES_NStokes


  !-----------------------------------------------------------------------------
  ! main level wrapper to set the right hand side on a block. Note this is completely
  ! independent of the grid any an MPI formalism, neighboring relations and the like.
  ! You just get a block data (e.g. ux, uy, uz, p) and compute the right hand side
  ! from that. Ghost nodes are assumed to be sync'ed.
  !-----------------------------------------------------------------------------
  subroutine RHS_NStokes( time, u, g, x0, dx, rhs, stage )
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

    ! output. Note assumed-shape arrays
    real(kind=rk), intent(inout) :: rhs(1:,1:,1:,1:)

    ! stage. there is 3 stages, init_stage, integral_stage and local_stage. If the PDE has
    ! terms that depend on global qtys, such as forces etc, which cannot be computed
    ! from a single block alone, the first stage does that. the second stage can then
    ! use these integral qtys for the actual RHS evaluation.
    character(len=*), intent(in)       :: stage
    ! Area of mean_density
    real(kind=rk)    ,save             :: integral(4),area


    ! local variables
    integer(kind=ik) :: Bs

    ! compute the size of blocks
    Bs = size(u,1) - 2*g

    select case(stage)
    case ("init_stage")
      !-------------------------------------------------------------------------
      ! 1st stage: init_stage.
      !-------------------------------------------------------------------------
      ! this stage is called only once, not for each block.
      ! performs initializations in the RHS module, such as resetting integrals
      integral= 0
      area    = 0

    case ("integral_stage")
      !-------------------------------------------------------------------------
      ! 2nd stage: init_stage.
      !-------------------------------------------------------------------------
      ! For some RHS, the eqn depend not only on local, block based qtys, such as
      ! the state vector, but also on the entire grid, for example to compute a
      ! global forcing term (e.g. in FSI the forces on bodies). As the physics
      ! modules cannot see the grid, (they only see blocks), in order to encapsulate
      ! them nicer, two RHS stages have to be defined: integral / local stage.
      !
      ! called for each block.
      if (params_ns%penalization .and. params_ns%geometry=="funnel") then
        rhs=u
        call convert_statevector2D(rhs(:,:,1,:),'pure_variables')
        call integrate_over_pump_area(rhs(:,:,1,:),g,Bs,x0,dx,integral,area) 
        rhs=0.0_rk
      endif


    case ("post_stage")
      !-------------------------------------------------------------------------
      ! 3rd stage: post_stage.
      !-------------------------------------------------------------------------
      ! this stage is called only once, not for each block.
      if (params_ns%penalization .and. params_ns%geometry=="funnel") then
        ! reduce sum on each block to global sum
        call mean_quant(integral,area)
        
      endif

    case ("local_stage")
      !-------------------------------------------------------------------------
      ! 4th stage: local evaluation of RHS on all blocks
      !-------------------------------------------------------------------------
      ! the second stage then is what you would usually do: evaluate local differential
      ! operators etc.


      ! called for each block.
      if (size(u,3)==1) then
        call  RHS_2D_navier_stokes(g, Bs,x0, (/dx(1),dx(2)/),u(:,:,1,:), rhs(:,:,1,:))
        
        if (params_ns%penalization) then
        ! add volume penalization 
          call add_constraints(rhs(:,:,1,:),Bs, g, x0,(/dx(1),dx(2)/),u(:,:,1,:))
        endif

      else
         call RHS_3D_navier_stokes(g, Bs,x0, (/dx(1),dx(2),dx(3)/), u, rhs)
      endif

    case default
      call abort(7771,"the RHS wrapper requests a stage this physics module cannot handle.")
    end select


  end subroutine RHS_NStokes

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  ! main level wrapper to compute statistics (such as mean flow, global energy,
  ! forces, but potentially also derived stuff such as Integral/Kolmogorov scales)
  ! NOTE: as for the RHS, some terms here depend on the grid as whole, and not just
  ! on individual blocks. This requires one to use the same staging concept as for the RHS.
  !-----------------------------------------------------------------------------
  subroutine STATISTICS_NStokes( time, u, g, x0, dx, stage )
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


    ! stage. there is 3 stages, init_stage, integral_stage and local_stage. If the PDE has
    ! terms that depend on global qtys, such as forces etc, which cannot be computed
    ! from a single block alone, the first stage does that. the second stage can then
    ! use these integral qtys for the actual RHS evaluation.
    character(len=*), intent(in) :: stage

    ! local variables
    integer(kind=ik)            :: Bs, mpierr,ix,iy
    real(kind=rk),save          :: area
    real(kind=rk), allocatable  :: mask(:,:)
    real(kind=rk)               :: eps_inv,tmp(3),y,x,r

    ! compute the size of blocks
    Bs = size(u,1) - 2*g

    select case(stage)
    case ("init_stage")
      !-------------------------------------------------------------------------
      ! 1st stage: init_stage.
      !-------------------------------------------------------------------------
      ! this stage is called only once, NOT for each block.
      ! performs initializations in the RHS module, such as resetting integrals
      params_ns%mean_density  = 0.0_rk
      params_ns%mean_pressure = 0.0_rk
      area                    = 0.0_rk
    case ("integral_stage")
      !-------------------------------------------------------------------------
      ! 2nd stage: integral_stage.
      !-------------------------------------------------------------------------
      ! This stage contains all operations which are running on the blocks
      !
      ! called for each block.

      if (maxval(abs(u))>1.0e5) then
        call abort(6661,"ns fail: very very large values in state vector.")
      endif
      ! compute force
      allocate(mask(Bs+2*g, Bs+2*g))
      call get_mask(mask, x0, dx, Bs, g)


      if (size(u,3)==1) then
        ! compute density and pressure only in physical domain 
        tmp(1:3) =0.0_rk
        
        do iy=g+1, Bs+g
          y = dble(iy-(g+1)) * dx(2) + x0(2)
          do ix=g+1, Bs+g
            x = dble(ix-(g+1)) * dx(1) + x0(1)
            !if (mask(ix,iy)<1e-10) then
                  tmp(1) = tmp(1)   + u(ix,iy, 1, rhoF)**2
                  tmp(2) = tmp(2)   + u(ix,iy, 1, pF)
                  tmp(3) = tmp(3)   + 1.0_rk
            !endif
          enddo
        enddo

        params_ns%mean_density = params_ns%mean_density   + tmp(1)*dx(1)*dx(2)
        params_ns%mean_pressure= params_ns%mean_pressure  + tmp(2)*dx(1)*dx(2)
        area                   = area                     + tmp(3)*dx(1)*dx(2)
      endif ! NOTE: MPI_SUM is perfomed in the post_stage.

    case ("post_stage")
      !-------------------------------------------------------------------------
      ! 3rd stage: post_stage.
      !-------------------------------------------------------------------------
      ! this stage is called only once, NOT for each block.


      tmp(1) = params_ns%mean_density
      call MPI_ALLREDUCE(tmp(1), params_ns%mean_density, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
      tmp(2) = params_ns%mean_pressure
      call MPI_ALLREDUCE(tmp(2), params_ns%mean_pressure, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)
      tmp(3) = area
      call MPI_ALLREDUCE(tmp(3), area                   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpierr)


     
    
       if (params_ns%mpirank == 0) then
         ! write mean flow to disk...
         write(*,*) "area=",area/params_ns%Lx/params_ns%Ly,params_ns%mean_density/area ,params_ns%mean_pressure/area
         open(14,file='meandensity.t',status='unknown',position='append')
         write (14,'(4(es15.8,1x))') time, params_ns%mean_density
         close(14)

         ! write forces to disk...
         open(14,file='meanpressure.t',status='unknown',position='append')
         write (14,'(4(es15.8,1x))') time, params_ns%mean_pressure
         close(14)
       end if

     case default
       call abort(7772,"the STATISTICS wrapper requests a stage this physics module cannot handle.")
     end select


  end subroutine STATISTICS_NStokes



  !-----------------------------------------------------------------------------
  ! setting the time step is very physics-dependent. Sometimes you have a CFL like
  ! condition, sometimes not. So each physic module must be able to decide on its
  ! time step. This routine is called for all blocks, the smallest returned dt is used.
  !-----------------------------------------------------------------------------
  subroutine GET_DT_BLOCK_NStokes( time, u, Bs, g, x0, dx, dt )
    implicit none

    ! it may happen that some source terms have an explicit time-dependency
    ! therefore the general call has to pass time
    real(kind=rk), intent (in) :: time

    ! block data, containg the state vector. In general a 4D field (3 dims+components)
    ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
    real(kind=rk), intent(in) :: u(1:,1:,1:,1:)

    ! as you are allowed to compute the RHS only in the interior of the field
    ! you also need to know where 'interior' starts: so we pass the number of ghost points
    integer, intent(in) :: g, bs

    ! for each block, you'll need to know where it lies in physical space. The first
    ! non-ghost point has the coordinate x0, from then on its just cartesian with dx spacing
    real(kind=rk), intent(in) :: x0(1:3), dx(1:3)

    ! the dt for this block is returned to the caller:
    real(kind=rk), intent(out) :: dt

    ! local variables
    real(kind=rk),allocatable  :: v_physical(:,:,:)
    real(kind=rk) :: deltax

    dt = 9.9e9_rk

    ! get smallest spatial seperation
    if(size(u,3)==1) then
        deltax=minval(dx(1:2))
        allocate(v_physical(2*g+Bs,2*g+Bs,1))
    else
        deltax=minval(dx)
        allocate(v_physical(2*g+Bs,2*g+Bs,2*g+Bs))
    endif

    ! calculate norm of velocity at every spatial point
    if (UzF==-1) then
        v_physical = u(:,:,:,UxF)*u(:,:,:,UxF) + u(:,:,:,UyF)*u(:,:,:,UyF)
    else
        v_physical = u(:,:,:,UxF)*u(:,:,:,UxF) + u(:,:,:,UyF)*u(:,:,:,UyF)+u(:,:,:,UzF)*u(:,:,:,UzF)
    endif

    ! maximal characteristical velocity is u+c where c = sqrt(gamma*p/rho) (speed of sound)
    v_physical = sqrt(v_physical)+sqrt(params_ns%gamma_*u(:,:,:,pF))
    v_physical = v_physical/u(:,:,:,rhoF)

    ! CFL criteria CFL=v_physical/v_numerical where v_numerical=dx/dt
     dt = min(dt, params_ns%CFL * deltax / maxval(v_physical))

    ! penalization requiers dt <= C_eta
    if (params_ns%penalization ) then
        dt=min(dt,params_ns%C_eta)
    endif
    ! penalization requiers dt <= C_eta
    if (params_ns%sponge_layer ) then
        dt=min(dt,params_ns%C_sp)
    endif
    deallocate(v_physical)
  end subroutine GET_DT_BLOCK_NStokes


  !-----------------------------------------------------------------------------
  ! main level wrapper for setting the initial condition on a block
  !-----------------------------------------------------------------------------
  subroutine INICOND_NStokes( time, u, g, x0, dx )
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

    integer(kind=ik)          :: Bs,ix
    real(kind=rk)             :: x


    ! compute the size of blocks
    Bs = size(u,1) - 2*g

    u = 0.0_rk

    if (p_init<=0.0_rk .or. rho_init <=0.0) then
      call abort(6032, "Error [module_navier_stokes_new.f90]: initial pressure and density must be larger then 0")
    endif
    


    select case( params_ns%inicond )
    case ("sinus_2d","sinus2d","sin2d")
      !> \todo implement sinus_2d inicondition
      call abort(7771,"inicond is not implemented yet: "//trim(adjustl(params_ns%inicond)))
    case ("zeros")
      ! add ambient pressure
      u( :, :, :, pF) = p_init
      ! set rho
      u( :, :, :, rhoF) = sqrt(rho_init)
      ! set Ux
      u( :, :, :, UxF) = 0.0_rk
      ! set Uy
      u( :, :, :, UyF) = 0.0_rk

      if (size(u,3).ne.1) then
          ! set Uz to zero
          u( :, :, :, UzF) = 0.0_rk
      endif
    case ("sod_shock_tube")
      ! Sods test case: shock tube
      ! ---------------------------
      !
      ! Test case for shock capturing filter
      ! The initial condition is devided into
      ! Left part x<= Lx/2 
      !
      ! rho=1
      ! p  =1
      ! u  =0
      ! 
      ! Rigth part x> Lx/2
      !
      ! rho=0.125
      ! p  =0.1
      ! u  =0
      do ix=1, Bs+2*g
         x = dble(ix-(g+1)) * dx(1) + x0(1)
         ! left region
         if (x <= params_ns%Lx*0.5_rk) then
           u( ix, :, :, rhoF) = 1.0_rk
           u( ix, :, :, pF)   = 1.0_rk
         else
           u( ix, :, :, rhoF) = sqrt(0.125_rk)
           u( ix, :, :, pF)   = 0.1_rk
         endif
      end do

      ! velocity set to 0
       u( :, :, :, UxF) = 0.0_rk
       u( :, :, :, UyF) = 0.0_rk

    case ("mask")
      ! add ambient pressure
      u( :, :, :, pF) = p_init
      ! set rho
      u( :, :, :, rhoF) = sqrt(rho_init)

      ! set velocity field u(x)=1 for x in mask
      if (size(u,3)==1 .and. params_ns%penalization) then
        call get_mask(u( :, :, 1, UxF), x0, dx, Bs, g )
        call get_mask(u( :, :, 1, UyF), x0, dx, Bs, g )
      endif

      ! u(x)=(1-mask(x))*u0 to make sure that flow is zero at mask values
      u( :, :, :, UxF) = (1-u(:,:,:,UxF))*u_init(1)*sqrt(rho_init) !flow in x
      u( :, :, :, UyF) = (1-u(:,:,:,UyF))*u_init(2)*sqrt(rho_init) !flow in y
    case ("pressure_blob")

        call inicond_gauss_blob( params_ns%inicond_width,Bs,g,(/ params_ns%Lx, params_ns%Ly, params_ns%Lz/), u(:,:,:,pF), x0, dx )
        ! add ambient pressure
        u( :, :, :, pF) = p_init + 1000.0_rk * u( :, :, :, pF)
        ! set rho
        u( :, :, :, rhoF) = sqrt(rho_init)
        ! set Ux
        u( :, :, :, UxF) = 0.0_rk
        ! set Uy
        u( :, :, :, UyF) = 0.0_rk

        if (size(u,3).ne.1) then
            ! set Uz to zero
            u( :, :, :, UzF) = 0.0_rk
        endif
    case default
        call abort(7771,"the initial condition is unkown: "//trim(adjustl(params_ns%inicond)))
    end select

  end subroutine INICOND_NStokes





 !-----------------------------------------------------------------------------
  ! main level wrapper to set the right hand side on a block. Note this is completely
  ! independent of the grid any an MPI formalism, neighboring relations and the like.
  ! You just get a block data (e.g. ux, uy, uz, p) and compute the right hand side
  ! from that. Ghost nodes are assumed to be sync'ed.
  !-----------------------------------------------------------------------------
  subroutine filter_NStokes( time, u, g, x0, dx, work_array )
    implicit none
    ! it may happen that some source terms have an explicit time-dependency
    ! therefore the general call has to pass time
    real(kind=rk), intent (in) :: time

    ! block data, containg the state vector. In general a 4D field (3 dims+components)
    ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
    real(kind=rk), intent(inout) :: u(1:,1:,1:,1:)

    ! as you are allowed to compute the work_array only in the interior of the field
    ! you also need to know where 'interior' starts: so we pass the number of ghost points
    integer, intent(in) :: g

    ! for each block, you'll need to know where it lies in physical space. The first
    ! non-ghost point has the coordinate x0, from then on its just cartesian with dx spacing
    real(kind=rk), intent(in) :: x0(1:3), dx(1:3)

    ! output. Note assumed-shape arrays
    real(kind=rk), intent(inout) :: work_array(1:,1:,1:,1:)


    ! local variables
    integer(kind=ik) :: Bs

    ! compute the size of blocks
    Bs = size(u,1) - 2*g

   
    call filter_block(params_ns%filter, time, u, Bs, g, x0, dx, work_array)
   

  end subroutine filter_NStokes




!> \brief convert the statevector \f$(\sqrt(rho),\sqrt(rho)u,\sqrt(rho)v,p )\f$
!! to the desired format.
subroutine convert_statevector2D(phi,convert2format)

    implicit none
    ! convert to type "conservative","pure_variables"
    character(len=*), intent(in)   :: convert2format
    !phi U=(sqrt(rho),sqrt(rho)u,sqrt(rho)v,sqrt(rho)w,p )
    real(kind=rk), intent(inout)      :: phi(1:,1:,1:)
    ! vector containing the variables in the desired format
    real(kind=rk)                  :: converted_vector(size(phi,1),size(phi,2),size(phi,3))

    select case( convert2format )
    case ("conservative") ! U=(rho, rho u, rho v, rho w, p)
      ! density 
      converted_vector(:,:,1)=phi(:,:,rhoF)**2
      ! rho u 
      converted_vector(:,:,2)=phi(:,:,UxF)*phi(:,:,rhoF)
      ! rho v 
      converted_vector(:,:,3)=phi(:,:,UyF)*phi(:,:,rhoF)
      ! kinetic energie
      converted_vector(:,:,4)=phi(:,:,UyF)**2+phi(:,:,UyF)**2
      converted_vector(:,:,4)=converted_vector(:,:,4)*0.5_rk
      ! e_tot=e_kin+p/(gamma-1)      
      converted_vector(:,:,4)=converted_vector(:,:,4)+phi(:,:,pF)/(params_ns%gamma_-1)
    case ("pure_variables")
      ! add ambient pressure
      !rho
      converted_vector(:,:,1)= phi(:,:,rhoF)**2
      !u
      converted_vector(:,:,2)= phi(:,:, UxF)/phi(:,:,rhoF)
      !v
      converted_vector(:,:,3)= phi(:,:, UyF)/phi(:,:,rhoF)
      !p
      converted_vector(:,:,4)= phi(:,:, pF)    
    case default
        call abort(7771,"the format is unkown: "//trim(adjustl(convert2format)))
    end select

    phi=converted_vector

end subroutine convert_statevector2D



subroutine pack_statevector2D(phi,format)
    implicit none
    ! convert to type "conservative","pure_variables"
    character(len=*), intent(in)   :: format
    !phi U=(sqrt(rho),sqrt(rho)u,sqrt(rho)v,sqrt(rho)w,p )
    real(kind=rk), intent(inout)      :: phi(1:,1:,1:)
    ! vector containing the variables in the desired format
    real(kind=rk)                  :: converted_vector(size(phi,1),size(phi,2),size(phi,3))

    select case( format )
    case ("conservative") ! phi=(rho, rho u, rho v, e_tot)
      ! sqrt(rho) 
      converted_vector(:,:,1)=sqrt(phi(:,:,1))
      ! sqrt(rho) u 
      converted_vector(:,:,2)=phi(:,:,2)/converted_vector(:,:,1)
      ! sqrt(rho) v 
      converted_vector(:,:,3)=phi(:,:,3)/converted_vector(:,:,1)
      ! kinetic energie
      converted_vector(:,:,4)=converted_vector(:,:,2)**2+converted_vector(:,:,3)**2
      converted_vector(:,:,4)=converted_vector(:,:,4)*0.5_rk
      ! p=(e_tot-e_kin)(gamma-1)/rho      
      converted_vector(:,:,4)=(phi(:,:,4)-converted_vector(:,:,4))*(params_ns%gamma_-1)
    case ("pure_variables") !phi=(rho,u,v,p)
      ! add ambient pressure
      ! sqrt(rho) 
      converted_vector(:,:,1)= sqrt(phi(:,:,1))
      ! sqrt(rho) u
      converted_vector(:,:,2)= phi(:,:, 2)*converted_vector(:,:,1)
      ! sqrt(rho)v
      converted_vector(:,:,3)= phi(:,:, 3)*converted_vector(:,:,1)
      !p
      converted_vector(:,:,4)= phi(:,:, 4)    
    case default
        call abort(7771,"the format is unkown: "//trim(adjustl(format)))
    end select

    phi(:,:,rhoF) =converted_vector(:,:,1)
    phi(:,:,UxF)  =converted_vector(:,:,2)
    phi(:,:,UyF)  =converted_vector(:,:,3)
    phi(:,:,pF)   =converted_vector(:,:,4)
end subroutine pack_statevector2D






end module module_navier_stokes_new
