!-----------------------------------------------------------------
!> \brief parameters and variables of Navier_Stokes physics
!-----------------------------------------------------------------
!> \details
!! This module contains all parameter types, and functions to init
!! the params structure visible to the navier stokes modules.
!> \details
!!    * params_ns
!!    * params_filter
!> \date 23.1.2018 -creation
!> \author P.Krah
!-----------------------------------------------------------------

module module_navier_stokes_params

  use module_precision
  use module_ini_files_parser_mpi
  use mpi

  implicit none

  !=========================================================
  !               FILTER PARAMETER
  !---------------------------------------------------------
  ! filter structure
  type :: type_params_filter
      ! name of filter
      character(len=80)                  :: name
      ! number or data fields to filter
      integer(kind=ik)                   :: number_data_fields

      ! explicit filters
      ! ----------------
      ! stencil array, note: size is fixed
      real(kind=rk)                      :: stencil(19)
      ! filter position (array postion of value to filter)
      integer(kind=ik)                   :: stencil_size

      !wavelet filter
      !---------------
      ! fine/coarse predictor of wavelets
      character(len=80)                  ::order_predictor
       ! threshold value for wavelets
      real(kind=rk)                     ::eps

      ! bogey shcok filter
      !---------------
      ! bogey shock detector threshold,
      real(kind=rk)                               :: r_th
      !  detector_method (divU,p),
      character(len=80)                           :: detector_method
      ! switch functio (tanh,abs)
      character(len=80)                           :: sigma_switch
      ! save filter strength sigma
      logical                                     :: save_filter_strength=.false.

  end type type_params_filter
  !=========================================================



  !=========================================================
  !               NAVIER STOKES PARAMETER
  !---------------------------------------------------------
  type :: type_params_ns
        ! coordinate system
        character(len=80)                           :: coordinates="cartesian"
        ! Courant-Friedrichs-Lewy
        real(kind=rk)                               :: CFL, T_end
        ! spatial domain%number_data_fields
        real(kind=rk)                               :: R_max, domain_size(3)=0.0_rk
        ! number data fields
        integer(kind=ik)                            :: number_data_fields
        ! number of block nodes
        integer(kind=ik)                            :: Bs
        ! maximal tree level
        integer(kind=ik)                            :: Jmax
        ! dimension
        integer(kind=ik)                            :: dim, N_fields_saved
        ! adiabatic coefficient
        real(kind=rk)                               :: gamma_=0
        ! specific gas constant
        real(kind=rk)                               :: Rs
        ! isochoric heat capacity
        real(kind=rk)                               :: Cv
        ! isobaric heat capacity
        real(kind=rk)                               :: Cp
        ! prandtl number
        real(kind=rk)                               :: Pr
        ! dynamic viscosity
        real(kind=rk)                               :: mu0
        ! dissipation switch
        logical                                     :: dissipation
        ! variable names
        character(len=80), allocatable              :: names(:)
        ! type used for the spatial discretization
        character(len=80)                           :: discretization
        !---------------------------------------------------------------------------------
        ! initial conditions
        !---------------------------------------------------------------------------------
        ! înitial condition
        character(len=80)                           :: inicond
        ! width of initialcond
        real(kind=rk)                               :: inicond_width
        ! width of initialcond
        real(kind=rk)                               :: initial_density
        ! pressure of initialcond
        real(kind=rk)                               :: initial_pressure

        real(kind=rk)                               :: initial_temp
        ! width of initialcond
        real(kind=rk)                               :: initial_velocity(3)
        ! --------------------------------------------------------------------------------
        ! penalization
        ! --------------------------------------------------------------------------------
        logical                                     :: penalization,smooth_mask=.True., sponge_layer
        ! penalization parameter and sponge parameter
        real(kind=rk)                               :: C_eta,C_sp
        ! geometry to display
        character(len=80)                           :: geometry="cylinder",case='--'
        ! geometric parameters for cylinder (x_0,r)
        real(kind=rk)                               :: x_cntr(1:3), R_cyl
        ! mean variables on domain
        real(kind=rk)                               :: mean_density
        real(kind=rk)                               :: mean_pressure
        real(kind=rk)                               :: force(3)
        ! we need to know which mpirank prints output..
        integer(kind=ik)                            :: mpirank, mpisize
        !---------------------------------------------------------
        ! filter structure containing all parameters of the filter
        type(type_params_filter)                    :: filter
        ! --------------------------------------------------------------------------------
        ! dimensionless scales (see function add_info)
        ! --------------------------------------------------------------------------------
        ! reynolds number is computed from intial conditions and Ly as characteristic length scales
        real(kind=rk)                               :: Reynolds
        ! speed of sound computed from inicond density, pressure
        real(kind=rk)                               :: c0
        ! machnumber
        real(kind=rk)                               :: machnumber
        ! smallest lattice spacing
        real(kind=rk)                               :: dx_min

  end type type_params_ns

  ! statevector index
  integer(kind=ik) ,save,public :: rhoF=-1,UxF=-1,UyF=-1,UzF=-1,pF=-1

  type(type_params_ns)          :: params_ns

contains

  include "initial_conditions.f90"


  subroutine init_navier_stokes_eq(params_ns, FILE )
    implicit none
    !> pointer to inifile
    type(inifile) ,intent(inout)     :: FILE
    !> params structure of navier stokes
    type(type_params_ns),intent(inout)  :: params_ns
    real(kind=rk), dimension(3)      :: domain_size=0.0_rk


    if (params_ns%mpirank==0) then
      write(*,*)
      write(*,*)
      write(*,*) "PARAMS: navier stokes equations!"
      write(*,'(" -----------------------------------")')
    endif
    ! read number_data_fields
    call read_param_mpi(FILE, 'Blocks', 'number_data_fields', params_ns%number_data_fields, 1 )
    ! dimension
    call read_param_mpi(FILE, 'Domain', 'dim', params_ns%dim, 2 )
    !
    call read_param_mpi(FILE, 'Navier_Stokes', 'Coordinate_system', params_ns%coordinates, &
                                                                    params_ns%coordinates )
    ! spatial domain size
    call read_param_mpi(FILE, 'Domain', 'dim', params_ns%dim, 2 )
    call read_param_mpi(FILE, 'Domain', 'domain_size', params_ns%domain_size(1:params_ns%dim), (/ 1.0_rk, 1.0_rk, 1.0_rk /) )

    if ( params_ns%coordinates=="cylindrical" ) then
      params_ns%R_max=params_ns%domain_size(2)*0.5_rk
      if ( params_ns%mpirank==0 ) write(*,'("maximal Radius" ,T30,"R_max",T60,"=",TR1, e10.2 )') params_ns%R_max
    end if

    ! read adiabatic coefficient
    call read_param_mpi(FILE, 'Navier_Stokes', 'gamma_', params_ns%gamma_, 0.0_rk )
    ! read specific gas constant
    call read_param_mpi(FILE, 'Navier_Stokes', 'Rs', params_ns%Rs, 0.0_rk )
    ! calculate isochoric heat capacity
    params_ns%Cv = params_ns%Rs/(params_ns%gamma_-1.0_rk)
    ! calculate isobaric heat capacity
    params_ns%Cp = params_ns%Cv*params_ns%gamma_
    ! read prandtl number
    call read_param_mpi(FILE, 'Navier_Stokes', 'Pr', params_ns%Pr, 0.0_rk )
    ! read dynamic viscosity
    call read_param_mpi(FILE, 'Navier_Stokes', 'mu0', params_ns%mu0, 0.0_rk )
    ! read switch to turn on|off dissipation
    call read_param_mpi(FILE, 'Navier_Stokes', 'dissipation', params_ns%dissipation, .true. )

  end subroutine init_navier_stokes_eq



  subroutine init_initial_conditions(params_ns, FILE )
      implicit none
      !> pointer to inifile
      type(inifile) ,intent(inout)        :: FILE
      !> params structure of navier stokes
      type(type_params_ns),intent(inout)  :: params_ns
      ! initial parameters
      real(kind=rk)                       :: rho_init=1,p_init=1,u_init(3)=0,T_init=1,width

      if (params_ns%mpirank==0) then
        write(*,*)
        write(*,*)
        write(*,*) "PARAMS: initial conditions"
        write(*,'(" ---------------------------")')
      endif
      call read_param_mpi(FILE, 'Physics', 'initial_cond', params_ns%inicond, "read_from_files" )
      if ( params_ns%inicond == "read_from_files") then
        if (params_ns%mpirank==0) write(*,'("initial configuration is read from file!")')
        if (params_ns%mpirank==0) write(*,'("we read in (rho,u , v, p) and convert it to skew: (sqrt(rho),sqrt(rho)u, sqrt(rho)v, p)!")')
        return
      end if
      call read_param_mpi(FILE, 'Navier_Stokes', 'inicond'      , params_ns%inicond, "pressure_blob" )
      call read_param_mpi(FILE, 'Navier_Stokes', 'inicond_width',width, params_ns%domain_size(1)*0.1_rk )
      call read_param_mpi(FILE, 'Navier_Stokes', 'initial_pressure' , p_init, p_init )
      call read_param_mpi(FILE, 'Navier_Stokes', 'initial_velocity' , u_init(1:params_ns%dim), u_init(1:params_ns%dim) )
      call read_param_mpi(FILE, 'Navier_Stokes', 'initial_temperature', T_init, T_init )
      call read_param_mpi(FILE, 'Navier_Stokes', 'initial_density', rho_init, rho_init )
      params_ns%initial_density=rho_init
      params_ns%initial_velocity=u_init
      params_ns%initial_pressure=p_init
      params_ns%initial_temp=T_init
       params_ns%inicond_width  =width

    end subroutine init_initial_conditions




subroutine init_other_params(params_ns, FILE )

    use module_helpers , only: list_contains_name
    implicit none
    !> pointer to inifile
    type(inifile) ,intent(inout)        :: FILE
    !> params structure of navier stokes
    type(type_params_ns),intent(inout)  :: params_ns

    if (params_ns%mpirank==0) then
      write(*,*)
      write(*,*)
      write(*,*) "PARAMS: additional data"
      write(*,'(" -----------------------")')
    endif
 ! allocate names list
    call read_param_mpi(FILE, 'Saving', 'N_fields_saved', params_ns%N_fields_saved, 1 )
    allocate( params_ns%names( params_ns%N_fields_saved ) )
    params_ns%names = "---"
    ! read file
    call read_param_mpi(FILE, 'Saving', 'field_names', params_ns%names, params_ns%names )
    if (  list_contains_name(params_ns%names,'sigmax')>0 ) then
      params_ns%filter%save_filter_strength=.true.
    end if

    call read_param_mpi(FILE, 'Blocks', 'number_data_fields', params_ns%number_data_fields, 1 )

    call read_param_mpi(FILE, 'Discretization', 'order_discretization', params_ns%discretization, "FD_2nd_central")

    call read_param_mpi(FILE, 'Time', 'CFL', params_ns%CFL, 1.0_rk   )
    call read_param_mpi(FILE, 'Time', 'time_max', params_ns%T_end, 1.0_rk)

    call read_param_mpi(FILE, 'Blocks', 'max_treelevel', params_ns%Jmax, 1   )
    call read_param_mpi(FILE, 'Blocks', 'number_block_nodes', params_ns%Bs, 1   )


  end subroutine init_other_params



  subroutine continue_periodic(x,L)
        !> position x
        real(kind=rk), intent(inout)     :: x
        !> domain length
        real(kind=rk), intent(in)     :: L

        real(kind=rk)                  :: min_dx

        if ( x>L ) then
          x=x-L
        elseif( x<0 ) then
          ! note it is actually x=L-abs(x) but since x is negative its
          x=L+x
        endif

        min_dx = 2.0_rk**(-params_ns%Jmax) * min(params_ns%domain_size(1),params_ns%domain_size(2))&
                          / real(params_ns%Bs-1, kind=rk)
        ! u(x=0) should be set equal to u(x=L)
        if ( abs(x-L)<min_dx*0.5_rk ) then
          x = 0.0_rk
        end if

  end subroutine continue_periodic


!> \brief Add additional info to navier stokes (initial Mach, Reynolds and speed of sound)
  subroutine add_info(params_ns )
      implicit none
      !> params structure of navier stokes
      type(type_params_ns),intent(inout)  :: params_ns
          integer(kind=ik)                :: nx_max
      ! compute min(dx,dy,dz)
      if ( params_ns%dim==2 ) then
        params_ns%dx_min = 2.0_rk**(-params_ns%Jmax) * min(params_ns%domain_size(1),params_ns%domain_size(2)) &
                                                          / real(params_ns%Bs-1, kind=rk)
      else
        params_ns%dx_min = 2.0_rk**(-params_ns%Jmax) * minval(params_ns%domain_size) &
                                                          / real(params_ns%Bs-1, kind=rk)
      end if

      ! initial speed of sound, Mach number, reynolds number
      params_ns%c0        = sqrt(params_ns%gamma_*params_ns%initial_pressure/params_ns%initial_density)
      params_ns%Machnumber= sqrt(params_ns%initial_velocity(1)**2 &
                                +params_ns%initial_velocity(2)**2 &
                                +params_ns%initial_velocity(3)**2)/params_ns%c0
      params_ns%Reynolds  = params_ns%initial_density*params_ns%domain_size(2)* &
                            params_ns%machnumber*params_ns%c0/params_ns%mu0

      if (params_ns%mpirank==0) then
        write(*,*)
        write(*,*)
        write(*,*) "Additional Information"
        write(*,'(" -----------------------")')
        nx_max = (params_ns%Bs-1) * 2**(params_ns%Jmax)
        write(*,'("minimal lattice spacing:",T40,g12.4)') params_ns%dx_min
        write(*,'("maximal resolution: ",T40,i5," x",i5)') nx_max, nx_max

        if (.not. params_ns%inicond=="read_from_files") then
            write(*,'("initial speed of sound:", T40, f6.2)') params_ns%c0
            write(*,'("initial Machnumber:", T40, f6.2)')     params_ns%machnumber
            write(*,'("Reynolds for Ly:", T40, f12.1)')       params_ns%Reynolds

        endif
      endif

    end subroutine add_info


    !> \brief Add additional info to navier stokes (initial Mach, Reynolds and speed of sound)
      subroutine check_parameters(params_ns )
          implicit none
          !> params structure of navier stokes
          type(type_params_ns),intent(inout)  :: params_ns


          if ( params_ns%dim==3 ) then
            if( params_ns%number_data_fields<5 ) call abort(9898,'Please increase number of data fields (min 5)')
            if( min(pF,UxF,UyF,UzF,rhoF)<0 )  call abort(9898,'Check names of data fields [p,Ux,Uy,Uz,rho]!')
          else
            if( params_ns%number_data_fields<4 ) call abort(9898,'Please increase number of data fields (min 4)')
            if( min(pF,UxF,UyF,rhoF)<0 )      call abort(9898,'Check names of data fields [p,Ux,Uy,rho]!')
          end if


    end subroutine check_parameters


    !> This function allocates a statevector like quantite.
    !> Hence a 2D(3D) array with 4(5) data fields
    subroutine allocate_statevector_ns(data,Bs,g)
        implicit none
        !data=(sqrt(rho),sqrt(rho)u,sqrt(rho)v,sqrt(rho)w,p )
        real(kind=rk), allocatable, intent(inout) :: data(:,:,:,:)
        integer(kind=ik),              intent(in) :: Bs,g

        if (.not.allocated(data)) then
          if (params_ns%dim==3) then
            allocate(data(1:Bs+2*g, 1:Bs+2*g, 1:Bs+2*g, params_ns%number_data_fields))
          else
            allocate(data(1:Bs+2*g, 1:Bs+2*g, 1, params_ns%number_data_fields))
        endif
      endif
    end subroutine allocate_statevector_ns



end module module_navier_stokes_params
