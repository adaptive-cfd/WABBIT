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

!> \brief Implementation of Navier Stokes Physiscs Interface for
!! WABBIT
module module_navier_stokes

  !---------------------------------------------------------------------------------------------
  ! modules
  use module_navier_stokes_params
  use module_helpers, only: block_contains_NaN
  use module_operators, only : compute_vorticity
  use module_ns_penalization
  use module_navier_stokes_cases

  implicit none

  ! I usually find it helpful to use the private keyword by itself initially, which specifies
  ! that everything within the module is private unless explicitly marked public.
  PRIVATE

  !**********************************************************************************************
  ! These are the important routines that are visible to WABBIT:
  !**********************************************************************************************
  PUBLIC :: READ_PARAMETERS_NSTOKES, PREPARE_SAVE_DATA_NSTOKES, RHS_NSTOKES, GET_DT_BLOCK_NSTOKES, &
            CONVERT_STATEVECTOR, PACK_STATEVECTOR, INICOND_NSTOKES, FIELD_NAMES_NStokes,&
            STATISTICS_NStokes,FILTER_NSTOKES,convert2format
  !**********************************************************************************************
  ! parameters for this module. they should not be seen outside this physics module
  ! in the rest of the code. WABBIT does not need to know them.



contains


  include "RHS_2D_navier_stokes.f90"
  include "RHS_3D_navier_stokes.f90"
  include "RHS_2D_cylinder.f90"
  include "filter_block.f90"
  include "inicond_NStokes.f90"
  include "save_data_ns.f90"
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
    integer(kind=ik)            :: mpicode,nx_max


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
    call init_initial_conditions(params_ns,file)
    ! initialice parameters and fields of the specific case study
    call read_case_parameters( params_ns , FILE )
    ! computes initial mach+reynolds number, speed of sound and smallest lattice spacing
    call add_info(params_ns)

    ! set global parameters pF,rohF, UxF etc
    do dF = 1, params_ns%n_eqn
                if ( params_ns%names(dF) == "p" ) pF = dF
                if ( params_ns%names(dF) == "rho" ) rhoF = dF
                if ( params_ns%names(dF) == "Ux" ) UxF = dF
                if ( params_ns%names(dF) == "Uy" ) UyF = dF
                if ( params_ns%names(dF) == "Uz" ) UzF = dF
    end do

    call check_parameters(params_ns)

    call clean_ini_file_mpi( FILE )

  end subroutine READ_PARAMETERS_NStokes








  !-----------------------------------------------------------------------------
  ! main level wrapper to set the right hand side on a block. Note this is completely
  ! independent of the grid and any MPI formalism, neighboring relations and the like.
  ! You just get a block data (e.g. ux, uy, uz, p) and compute the right hand side
  ! from that. Ghost nodes are assumed to be sync'ed.
  !-----------------------------------------------------------------------------
  subroutine RHS_NStokes( time, u, g, x0, dx, rhs, stage, boundary_flag )
    use module_funnel, only:mean_quantity, integrate_over_pump_area
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

    ! when implementing boundary conditions, it is necessary to now if the local field (block)
    ! is adjacent to a boundary, because the stencil has to be modified on the domain boundary.
    ! The boundary_flag tells you if the local field is adjacent to a domain boundary:
    ! boundary_flag(i) can be either 0, 1, -1,
    !  0: no boundary in the direction +/-e_i
    !  1: boundary in the direction +e_i
    ! -1: boundary in the direction - e_i
    ! currently only acessible in the local stage
    integer(kind=2)          , intent(in):: boundary_flag(3)

    ! Area of mean_density
    real(kind=rk)    ,save             :: integral(5),area


    ! local variables
    integer(kind=ik) :: Bs,dF

    ! compute the size of blocks
    Bs = size(u,1) - 2*g

    select case(stage)
    case ("init_stage")
      !-------------------------------------------------------------------------
      ! 1st stage: init_stage.
      !-------------------------------------------------------------------------
      ! this stage is called only once, not for each block.
      ! performs initializations in the RHS module, such as resetting integrals
      integral= 0.0_rk
      area    = 0.0_rk

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
      if (params_ns%case=="funnel") then
        ! since rhs was not computed yet we can use it as a temporary storage
        rhs=u
        call convert_statevector(rhs(:,:,:,:),'pure_variables')
        call integrate_over_pump_area(rhs(:,:,:,:),g,Bs,x0,dx,integral,area)
      endif

    case ("post_stage")
      !-------------------------------------------------------------------------
      ! 3rd stage: post_stage.
      !-------------------------------------------------------------------------
      ! this stage is called only once, not for each block.
      if (params_ns%case=="funnel") then
        ! reduce sum on each block to global sum
        call mean_quantity(integral,area)
      endif

    case ("local_stage")
      !-------------------------------------------------------------------------
      ! 4th stage: local evaluation of RHS on all blocks
      !-------------------------------------------------------------------------
      ! the second stage then is what you would usually do: evaluate local differential
      ! operators etc.

      ! called for each block.
      if (params_ns%dim==2) then
        select case(params_ns%coordinates)
        case ("cartesian")
          call  RHS_2D_navier_stokes(g, Bs,x0, (/dx(1),dx(2)/),u(:,:,1,:), rhs(:,:,1,:),boundary_flag)
        case("cylindrical")
          call RHS_2D_cylinder(g, Bs,x0, (/dx(1),dx(2)/),u(:,:,1,:), rhs(:,:,1,:))
        case default
          call abort(7772,"ERROR [module_navier_stokes]: This coordinate system is not known!")
        end select
        !call  RHS_1D_navier_stokes(g, Bs,x0, (/dx(1),dx(2)/),u(:,:,1,:), rhs(:,:,1,:))
      else
         call RHS_3D_navier_stokes(g, Bs,x0, (/dx(1),dx(2),dx(3)/), u, rhs)
      endif

      if (params_ns%penalization) then
        ! add volume penalization
        call add_constraints(params_ns, rhs, Bs, g, x0, dx, u)
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
    real(kind=rk), allocatable  :: mask(:,:,:)
    real(kind=rk)               :: eta_inv,tmp(5),y,x,r

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
      params_ns%force         = 0.0_rk
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
      ! compute mean density and pressure
      if(.not. allocated(mask)) then
       if (params_ns%dim==2) allocate(mask(Bs+2*g, Bs+2*g, 1))
       if (params_ns%dim==3) allocate(mask(Bs+2*g, Bs+2*g, Bs+2*g))
     endif
      if ( params_ns%penalization ) then
        call get_mask(params_ns, x0, dx, Bs, g , mask)
      else
        mask=0.0_rk
      end if

      eta_inv                 = 1.0_rk/params_ns%C_eta

      if (params_ns%dim==2) then
        ! compute density and pressure only in physical domain
        tmp(1:5) =0.0_rk
        ! we do not want to sum over redudant points so exclude Bs+g!!!
        do iy=g+1, Bs+g-1
          y = dble(iy-(g+1)) * dx(2) + x0(2)
          do ix=g+1, Bs+g-1
            x = dble(ix-(g+1)) * dx(1) + x0(1)
            if (mask(ix,iy,1)<1e-10) then
                  tmp(1) = tmp(1)   + u(ix,iy, 1, rhoF)**2
                  tmp(2) = tmp(2)   + u(ix,iy, 1, pF)
                  tmp(5) = tmp(5)   + 1.0_rk
            endif
            ! force on obstacle (see Boiron)
            !Fx=1/Ceta mask*rho*u
            tmp(3) = tmp(3)   + u(ix,iy, 1, rhoF)*u(ix,iy, 1, UxF)*mask(ix,iy,1)
            !Fy=1/Ceta mask*rho*v
            tmp(4) = tmp(4)   + u(ix,iy, 1, rhoF)*u(ix,iy, 1, UyF)*mask(ix,iy,1)

          enddo
        enddo

        params_ns%mean_density = params_ns%mean_density   + tmp(1)*dx(1)*dx(2)
        params_ns%mean_pressure= params_ns%mean_pressure  + tmp(2)*dx(1)*dx(2)
        params_ns%force(1)     = params_ns%force(1)       + tmp(3)*dx(1)*dx(2)*eta_inv
        params_ns%force(2)     = params_ns%force(2)       + tmp(4)*dx(1)*dx(2)*eta_inv
        params_ns%force(3)     = 0
        area                   = area                     + tmp(5)*dx(1)*dx(2)
      endif ! NOTE: MPI_SUM is perfomed in the post_stage.

    case ("post_stage")
      !-------------------------------------------------------------------------
      ! 3rd stage: post_stage.
      !-------------------------------------------------------------------------
      ! this stage is called only once, NOT for each block.


      tmp(1) = params_ns%mean_density
      call MPI_ALLREDUCE(tmp(1), params_ns%mean_density, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
      tmp(2) = params_ns%mean_pressure
      call MPI_ALLREDUCE(tmp(2), params_ns%mean_pressure, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
      tmp(3) = params_ns%force(1)
      call MPI_ALLREDUCE(tmp(3), params_ns%Force(1)     , 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
      tmp(4) = params_ns%force(2)
      call MPI_ALLREDUCE(tmp(4), params_ns%Force(2)     , 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
      tmp(5) = area
      call MPI_ALLREDUCE(tmp(5), area                   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)




       if (params_ns%mpirank == 0) then
         ! write mean flow to disk...
         write(*,*) 'density=', params_ns%mean_density/area ,&
                    'pressure=',params_ns%mean_pressure/area, &
                    'drag=',params_ns%force(1),&!*2/params_ns%initial_density/params_ns%initial_velocity(1)**2/0.01, &
                    'Fy=',params_ns%force(2)
         open(14,file='meandensity.t',status='unknown',position='append')
         write (14,'(2(es15.8,1x))') time, params_ns%mean_density/area
         close(14)

         ! write mean Force
         open(14,file='Force.t',status='unknown',position='append')
         write (14,'(4(es15.8,1x))') time, params_ns%force
         close(14)

         ! write forces to disk...
         open(14,file='meanpressure.t',status='unknown',position='append')
         write (14,'(2(es15.8,1x))') time, params_ns%mean_pressure/area
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
    real(kind=rk),allocatable,save  :: v_physical(:,:,:)
    real(kind=rk)                   :: dx_min


    dt = 9.9e9_rk
    dx_min=minval(dx(1:params_ns%dim))

    if (maxval(abs(u))>1.0e7 .OR. minval(u(:,:,:,pF))<0 ) then
         call abort(65761,"ERROR [module_navier_stokes.f90]: statevector values out of physical range")
    endif
    if(params_ns%dim==2) then
      if( .not. allocated(v_physical))  allocate(v_physical(2*g+Bs,2*g+Bs,1))
      v_physical = u(:,:,:,UxF)*u(:,:,:,UxF) + u(:,:,:,UyF)*u(:,:,:,UyF)
    else
      if( .not. allocated(v_physical))  allocate(v_physical(2*g+Bs,2*g+Bs,2*g+Bs))
      v_physical = u(:,:,:,UxF)*u(:,:,:,UxF) + u(:,:,:,UyF)*u(:,:,:,UyF)+u(:,:,:,UzF)*u(:,:,:,UzF)
    endif

    v_physical = sqrt(v_physical)+sqrt(params_ns%gamma_*u(:,:,:,pF)) ! v= sqrt(rho u^2) + sqrt(gamma p)
    v_physical = v_physical/u(:,:,:,rhoF)                            ! v= (sqrt(rho u^2) + sqrt (gamma p))/sqrt(rho)

    ! CFL criteria CFL=v_physical/v_numerical where v_numerical=dx/dt
     dt = min(dt, params_ns%CFL * dx_min / maxval(v_physical))
    ! penalization requiers dt <= C_eta
    if (params_ns%penalization ) then
        dt=min(dt,params_ns%C_eta)
    endif
    ! penalization requiers dt <= C_eta
    if (params_ns%sponge_layer ) then
        dt=min(dt,params_ns%C_sp)
    endif
    !deallocate(v_physical)
  end subroutine GET_DT_BLOCK_NStokes






  !-----------------------------------------------------------------------------
  ! main level wrapper to filter a block. Note this is completely
  ! independent of the grid and any MPI formalism, neighboring relations and the like.
  ! You just get a block data (e.g. ux, uy, uz, p) and apply your filter to it.
  ! Ghost nodes are assumed to be sync'ed.
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

    call filter_block(params_ns%filter, time, u, g, Bs, x0, dx, work_array )
    u=work_array(:,:,:,1:params_ns%n_eqn)

  end subroutine filter_NStokes




!> \brief convert the statevector \f$(\sqrt(\rho),\sqrt(\rho)u,\sqrt(\rho)v,p )\f$
!! to the desired format.
subroutine convert_statevector(phi,convert2format)

    implicit none
    ! convert to type "conservative","pure_variables"
    character(len=*), intent(in)   :: convert2format
    !phi U=(sqrt(rho),sqrt(rho)u,sqrt(rho)v,sqrt(rho)w,p )
    real(kind=rk), intent(inout)      :: phi(1:,1:,1:,1:)
    ! vector containing the variables in the desired format
    real(kind=rk)                  :: converted_vector(size(phi,1),size(phi,2),size(phi,3),size(phi,4))


    select case( convert2format )
    case ("conservative") ! U=(rho, rho u, rho v, rho w, p)
      ! density
      converted_vector(:,:,:,rhoF)  =phi(:,:,:,rhoF)**2
      ! rho u
      converted_vector(:,:,:,UxF)   =phi(:,:,:,UxF)*phi(:,:,:,rhoF)
      ! rho v
      converted_vector(:,:,:,UyF)   =phi(:,:,:,UyF)*phi(:,:,:,rhoF)

      if ( params_ns%dim==3 ) then
        ! rho w
        converted_vector(:,:,:,UzF)=phi(:,:,:,UzF)*phi(:,:,:,rhoF)
        !kinetic energie
        converted_vector(:,:,:,pF)=phi(:,:,:,UxF)**2+phi(:,:,:,UyF)**2+phi(:,:,:,UzF)**2
      else
        ! kinetic energie
        converted_vector(:,:,:,pF)=phi(:,:,:,UxF)**2+phi(:,:,:,UyF)**2
      end if
      converted_vector(:,:,:,pF)=converted_vector(:,:,:,pF)*0.5_rk
      ! total energie
      ! e_tot=e_kin+p/(gamma-1)
      converted_vector(:,:,:,pF)=converted_vector(:,:,:,pF)+phi(:,:,:,pF)/(params_ns%gamma_-1)
    case ("pure_variables")
      ! add ambient pressure
      !rho
      converted_vector(:,:,:,rhoF)= phi(:,:,:,rhoF)**2
      !u
      converted_vector(:,:,:,UxF)= phi(:,:,:, UxF)/phi(:,:,:,rhoF)
      !v
      converted_vector(:,:,:,UyF)= phi(:,:,:, UyF)/phi(:,:,:,rhoF)
      !w
      if ( params_ns%dim==3 ) converted_vector(:,:,:,UzF)= phi(:,:,:, UzF)/phi(:,:,:,rhoF)
      !p
      converted_vector(:,:,:,pF)= phi(:,:,:, pF)
    case default
        call abort(7771,"the format is unkown: "//trim(adjustl(convert2format)))
    end select

    phi=converted_vector

end subroutine convert_statevector


!> \brief pack statevector of skewsymetric scheme \f$(\sqrt(\rho),\sqrt(\rho)u,\sqrt(\rho)v,p )\f$ from
!>            + conservative variables \f$(\rho,\rho u,\rho v,e\rho )\f$ or pure variables (rho,u,v,p)
subroutine pack_statevector(phi,format)
    implicit none
    ! convert to type "conservative","pure_variables"
    character(len=*), intent(in)   :: format
    !phi U=(sqrt(rho),sqrt(rho)u,sqrt(rho)v,sqrt(rho)w,p )
    real(kind=rk), intent(inout)   :: phi(1:,1:,1:,1:)
    ! vector containing the variables in the desired format
    real(kind=rk)                  :: converted_vector(size(phi,1),size(phi,2),size(phi,3),size(phi,4))



    select case( format )
    case ("conservative") ! phi=(rho, rho u, rho v, e_tot)
      ! sqrt(rho)
      if ( minval(phi(:,:,:,rhoF))<0 ) then
        write(*,*) "minval=", minval(phi(:,:,:,rhoF))
        call abort(457881,"ERROR [module_navier_stokes.f90]: density smaller then 0!!")
      end if
      converted_vector(:,:,:,rhoF)=sqrt(phi(:,:,:,rhoF))
      ! sqrt(rho) u
      converted_vector(:,:,:,UxF)=phi(:,:,:,UxF)/converted_vector(:,:,:,rhoF)
      ! sqrt(rho) v
      converted_vector(:,:,:,UyF)=phi(:,:,:,UyF)/converted_vector(:,:,:,rhoF)

      if ( params_ns%dim==3 ) then
        converted_vector(:,:,:,UzF)=phi(:,:,:,UzF)/converted_vector(:,:,:,rhoF)
        ! kinetic energie
        converted_vector(:,:,:,pF)=converted_vector(:,:,:,UxF)**2 &
                                  +converted_vector(:,:,:,UyF)**2 &
                                  +converted_vector(:,:,:,UzF)**2
      else
        ! kinetic energie
        converted_vector(:,:,:,pF)=converted_vector(:,:,:,UxF)**2 &
                                  +converted_vector(:,:,:,UyF)**2
      end if
      converted_vector(:,:,:,pF)=converted_vector(:,:,:,pF)*0.5_rk
      ! p=(e_tot-e_kin)(gamma-1)/rho
      converted_vector(:,:,:,pF)=(phi(:,:,:,pF)-converted_vector(:,:,:,pF))*(params_ns%gamma_-1)
    case ("pure_variables") !phi=(rho,u,v,p)
      ! add ambient pressure
      ! sqrt(rho)
      converted_vector(:,:,:,rhoF)= sqrt(phi(:,:,:,rhoF))
      ! sqrt(rho) u
      converted_vector(:,:,:,UxF)= phi(:,:,:,UxF)*converted_vector(:,:,:,rhoF)
      ! sqrt(rho)v
      converted_vector(:,:,:,UyF)= phi(:,:,:,UyF)*converted_vector(:,:,:,rhoF)
      ! sqrt(rho)w
      if ( params_ns%dim==3 ) converted_vector(:,:,:,UzF)= phi(:,:,:,UzF)*converted_vector(:,:,:,rhoF)
      !p
      converted_vector(:,:,:,pF)= phi(:,:,:,pF)
    case default
        call abort(7771,"the format is unkown: "//trim(adjustl(format)))
    end select

    phi(:,:,:,rhoF) =converted_vector(:,:,:,rhoF)
    phi(:,:,:,UxF)  =converted_vector(:,:,:,UxF )
    phi(:,:,:,UyF)  =converted_vector(:,:,:,UyF )
    if ( params_ns%dim==3 ) phi(:,:,:,UzF) = converted_vector(:,:,:,UzF)
    phi(:,:,:,pF)   =converted_vector(:,:,:,pF)
end subroutine pack_statevector



subroutine convert2format(phi_in,format_in,phi_out,format_out)
    implicit none
    ! convert to type "conservative","pure_variables"
    character(len=*), intent(in)   ::  format_in, format_out
    !phi U=(sqrt(rho),sqrt(rho)u,sqrt(rho)v,sqrt(rho)w,p )
    real(kind=rk), intent(in)      :: phi_in(1:,1:,1:,1:)
    ! vector containing the variables in the desired format
    real(kind=rk), intent(inout)   :: phi_out(:,:,:,:)

    ! convert phi_in to skewsymetric variables  \f$(\sqrt(\rho),\sqrt(\rho)u,\sqrt(\rho)v,p )\f$
    if (format_in=="skew") then
      phi_out  =  phi_in
    else
      phi_out  =  phi_in
      call pack_statevector(phi_out,format_in)
    endif

    ! form skewsymetric variables convert to any other scheme
    if (format_out=="skew") then
      !do nothing because format is skew already
    else
      call convert_statevector(phi_out,format_out)
    endif
end subroutine convert2format






end module module_navier_stokes
