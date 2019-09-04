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
  use module_params, only: read_bs

  implicit none

  !=========================================================
  !               BOUNDARYCONDITIONS Parameter
  !---------------------------------------------------------
  ! filter structure
  type :: type_boundary

      ! Boundary type
      ! -----------------
      ! Type of boundary is an array of dimension 3.
      ! Each boundary direction e_i (i=1,2 in 2D)
      ! has a type:
      !                 + periodic (if ghostnodes in direction are periodically synched)
      !                 + symmetric-open (symmetric in -e_i and open in +e_i)
      !                 + open (open at both sides)
      character(len=80), dimension(3)  :: name

      ! Boundary values
      ! ----------------
      ! possible values of the statevector at the boundaries:
      !  + e_x direction
      real(kind=rk),allocatable  :: phi_xplus(:)
      !  + -e_x direction
      real(kind=rk),allocatable  :: phi_xminus(:)
      !  + e_y direction
      real(kind=rk),allocatable  :: phi_yplus(:)
      !  + -e_y direction
      real(kind=rk),allocatable  :: phi_yminus(:)
      !  + e_z direction
      real(kind=rk),allocatable  :: phi_zplus(:)
      !  + -e_z direction
      real(kind=rk),allocatable  :: phi_zminus(:)

      ! How to initialice Boundaryvalues?
      ! ----------------------------------
      ! Boundary values can be initialiced with parameters red from
      ! parameter file (choose "paramsfile") or taken from the initial condition (choose "inicond")
      ! Currently not supported, but possible: choose new BC in every timestep (choose "time-dependent")
      character(len=80)       :: ref_type

end type type_boundary
  !=========================================================




  !=========================================================
  !               FILTER PARAMETER
  !---------------------------------------------------------
  ! filter structure
  type :: type_params_filter
      ! name of filter
      character(len=80)                  :: name
      ! number or data fields to filter
      integer(kind=ik)                   :: n_eqn

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
        real(kind=rk)                               :: domain_size(3)=0.0_rk
        ! minimal and maximal radius in the case of cylindrical coordinates
        real(kind=rk)                               :: R_max, R_min
        ! number data fields
        integer(kind=ik)                            :: n_eqn
        ! number of block nodes
        integer(kind=ik)                            :: g
        integer(kind=ik), dimension(3)              :: Bs
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
        logical                                     :: read_from_files
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
        real(kind=rk)                               :: C_eta,C_sp,L_sponge=-1.0_rk
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
        integer(kind=ik)                            :: N_mask_components
        !---------------------------------------------------------
        ! filter structure containing all parameters of the filter
        type(type_params_filter)                    :: filter
        !---------------------------------------------------------
        ! filter structure containing all parameters of the filter
        type(type_boundary)                         :: bound
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
        ! -------------------------------------------------------------------------------------
        ! Boundary conditions
        ! -------------------------------------------------------------------------------------
        logical,dimension(3)                        :: periodic_BC=(/.true.,.true.,.true./)
  end type type_params_ns

  ! statevector index
  integer(kind=ik), save, public :: rhoF=-1,UxF=-1,UyF=-1,UzF=-1,pF=-1

  type(type_params_ns), save :: params_ns


  interface convert_statevector
      module procedure convert_statevector2D, convert_statevector3D
  end interface convert_statevector

contains


  subroutine init_navier_stokes_eq( FILE )
    implicit none
    !> pointer to inifile
    type(inifile) ,intent(inout)     :: FILE
    real(kind=rk), dimension(3)      :: domain_size=0.0_rk


    if (params_ns%mpirank==0) then
      write(*,*)
      write(*,*)
      write(*,*) "PARAMS: navier stokes equations!"
      write(*,'(" -----------------------------------")')
    endif
    ! read number_data_fields
    call read_param_mpi(FILE, 'Blocks', 'number_equations', params_ns%n_eqn, 1 )
    ! dimension
    call read_param_mpi(FILE, 'Domain', 'dim', params_ns%dim, 2 )
    call read_param_mpi(FILE, 'Navier_Stokes', 'Coordinate_system', params_ns%coordinates, &
                                                                    params_ns%coordinates )
    ! spatial domain size
    params_ns%domain_size = (/ 1.0_rk, 1.0_rk, 1.0_rk /)
    call read_param_mpi(FILE, 'Domain', 'domain_size', params_ns%domain_size(1:params_ns%dim), &
                        params_ns%domain_size(1:params_ns%dim))

    if ( params_ns%coordinates=="cylindrical" ) then
      params_ns%R_max=params_ns%domain_size(2)*0.5_rk
      if ( params_ns%mpirank==0 ) write(*,'("maximal Radius" ,T30,"R_max",T60,"=",TR1, e10.2 )') &
                                  params_ns%R_max
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
    ! which case is studied in the NStokes module
    call read_param_mpi(FILE, 'Navier_Stokes', 'case', params_ns%case, 'no' )

  end subroutine init_navier_stokes_eq



  subroutine init_initial_conditions( FILE )
      implicit none
      !> pointer to inifile
      type(inifile) ,intent(inout)        :: FILE
      ! initial parameters
      real(kind=rk)                       :: rho_init=-1,p_init=-1,u_init(3)=0,T_init=-1,width
      character(len=*),parameter         :: section='Initial_Values'

      if (params_ns%mpirank==0) then
        write(*,*)
        write(*,*)
        write(*,*) "PARAMS: "//section
        write(*,'(" ---------------------------")')
      endif
      call read_param_mpi(FILE, 'Physics', 'read_from_files', params_ns%read_from_files, .false. )
      if ( params_ns%read_from_files) then
        if (params_ns%mpirank==0) write(*,'("initial configuration is read from file!")')
        if (params_ns%mpirank==0) write(*,'("we read in (rho,u , v, p) and convert it to skew: (sqrt(rho),sqrt(rho)u, sqrt(rho)v, p)!")')
        return
      end if
      call read_param_mpi(FILE, section, 'inicond'      , params_ns%inicond, "no" )
      call read_param_mpi(FILE, section, 'inicond_width',width, params_ns%domain_size(1)*0.1_rk )
      call read_param_mpi(FILE, section, 'initial_pressure' , p_init, p_init )
      call read_param_mpi(FILE, section, 'initial_velocity' , u_init(1:params_ns%dim), u_init(1:params_ns%dim) )
      call read_param_mpi(FILE, section, 'initial_temperature', T_init, T_init )
      call read_param_mpi(FILE, section, 'initial_density', rho_init, rho_init )
      params_ns%initial_density=rho_init
      params_ns%initial_velocity=u_init
      params_ns%initial_pressure=p_init
      params_ns%initial_temp=T_init
       params_ns%inicond_width  =width

    end subroutine init_initial_conditions


subroutine read_boundary_conditions( FILE )
    implicit none
    !------------------------------------------------
    type(inifile) ,intent(inout)        :: FILE !> pointer to inifile
    !------------------------------------------------
    integer :: dim,i

    dim=params_ns%dim

    call read_param_mpi(FILE, 'Domain', 'periodic_BC', params_ns%periodic_BC(1:dim), &
                                                       params_ns%periodic_BC(1:dim) )

    ! only read in parameters if we need them
    if ( ALL(params_ns%periodic_BC) ) then
      params_ns%bound%name="periodic"
      return
    end if

    if (params_ns%mpirank==0) then
      write(*,*)
      write(*,*)
      write(*,*) "PARAMS: Boundary Conditions"
      write(*,'(" ---------------------------")')
    endif

    params_ns%bound%name="---"
    call read_param_mpi(FILE, 'Boundary_Conditions', 'boundary_type', params_ns%bound%name(1:dim), &
                                                                params_ns%bound%name(1:dim) )

    do i = 1, dim
      if ( .not. params_ns%periodic_BC(i) .and. params_ns%bound%name(i)=="periodic" ) then
        call abort(251020181,"ERROR: Try to asign periodic BC to a boundary which is not periodic!"// &
                              "Tip: Change either boundary_type or parameter periodic_BC")
      end if
    end do


    ! ==========================================
    ! set boundary values
    ! ==========================================
    ! call read_param_mpi(FILE, 'Boundary_Conditions', 'reference_type',phi_tmp,phi_tmp)
    ! if ( .not.params_ns%bound%ref_type=='paramsfile' ) then
    !       return
    ! end if

    ! + at the e_x direction boundaries
    if ( .not. params_ns%periodic_BC(1) ) then
      select case(params_ns%bound%name(1))
      case("open")
            allocate(params_ns%bound%phi_xminus(params_ns%n_eqn))
            params_ns%bound%phi_xminus=-222.222_rk
            call read_param_mpi(FILE, 'Boundary_Conditions', 'state_xminus', &
                                    params_ns%bound%phi_xminus,params_ns%bound%phi_xminus)

            allocate(params_ns%bound%phi_xplus(params_ns%n_eqn))
            params_ns%bound%phi_xplus=-222.222_rk
            call read_param_mpi(FILE, 'Boundary_Conditions', 'state_xplus', &
                                    params_ns%bound%phi_xplus,params_ns%bound%phi_xplus)

      case("symmetric-open")
            ! no boundary conditions at yminus since yminus is symetric
            allocate(params_ns%bound%phi_xplus(params_ns%n_eqn))
            params_ns%bound%phi_xplus=-222.222_rk
            call read_param_mpi(FILE, 'Boundary_Conditions', 'state_xplus', &
                                    params_ns%bound%phi_xplus,params_ns%bound%phi_xplus)

      case("wall")
            ! to implement
      case default
            call abort(81020166,"OHHHH no, Unknown Boundary Condition: "// params_ns%bound%name(1))
      end select
    endif

    ! + at the e_y direction boundaries
    if ( .not. params_ns%periodic_BC(2) ) then
      select case(params_ns%bound%name(2))
      case("open")
            allocate(params_ns%bound%phi_yminus(params_ns%n_eqn))
            params_ns%bound%phi_yminus=-222.222_rk
            call read_param_mpi(FILE, 'Boundary_Conditions', 'state_yminus', &
                                    params_ns%bound%phi_yminus,params_ns%bound%phi_yminus)

            allocate(params_ns%bound%phi_yplus(params_ns%n_eqn))
            params_ns%bound%phi_yplus=-222.222_rk
            call read_param_mpi(FILE, 'Boundary_Conditions', 'state_yplus', &
                                    params_ns%bound%phi_yplus,params_ns%bound%phi_yplus)
      case("symmetric-open")
            ! no boundary conditions at yminus since yminus is symetric
            allocate(params_ns%bound%phi_yplus(params_ns%n_eqn))
            params_ns%bound%phi_yplus=-222.222_rk
            call read_param_mpi(FILE, 'Boundary_Conditions', 'state_yplus', &
                                    params_ns%bound%phi_yplus,params_ns%bound%phi_yplus)

      case("symmetryAxis-wall")
           ! to implement
      case default
           call abort(81020162,"OHHHH no, Unknown Boundary Condition: "// params_ns%bound%name(2))
      end select
    endif


end subroutine read_boundary_conditions


subroutine init_other_params( FILE )

    use module_helpers , only: list_contains_name
    implicit none
    !> pointer to inifile
    type(inifile) ,intent(inout)        :: FILE

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

    call read_param_mpi(FILE, 'Discretization', 'order_discretization', params_ns%discretization, "FD_2nd_central")

    call read_param_mpi(FILE, 'Time', 'CFL', params_ns%CFL, 1.0_rk   )
    call read_param_mpi(FILE, 'Time', 'time_max', params_ns%T_end, 1.0_rk)

    call read_param_mpi(FILE, 'Blocks', 'max_treelevel', params_ns%Jmax, 1   )

    params_ns%Bs = read_bs(FILE,'Blocks', 'number_block_nodes', params_ns%Bs, params_ns%dim)

    call read_param_mpi(FILE, 'Blocks', 'number_ghost_nodes', params_ns%g, 1   )


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
                          / real(params_ns%Bs(1)-1, kind=rk)
        ! u(x=0) should be set equal to u(x=L)
        if ( abs(x-L)<min_dx*0.5_rk ) then
          x = 0.0_rk
        end if

  end subroutine continue_periodic


!> \brief Add additional info to navier stokes (initial Mach, Reynolds and speed of sound)
  subroutine add_info(  )
      implicit none
      integer(kind=ik)                :: nx_max

      ! compute min(dx,dy,dz)
      if ( params_ns%dim==2 ) then
        params_ns%dx_min = 2.0_rk**(-params_ns%Jmax) * min(params_ns%domain_size(1),params_ns%domain_size(2)) &
                                                          / real(params_ns%Bs(1)-1, kind=rk)
      else
        params_ns%dx_min = 2.0_rk**(-params_ns%Jmax) * minval(params_ns%domain_size) &
                                                          / real(params_ns%Bs(1)-1, kind=rk)
      end if

      ! initial speed of sound, Mach number, reynolds number
      if ( params_ns%initial_density>0 ) then
        params_ns%c0        = sqrt(params_ns%gamma_*params_ns%initial_pressure/params_ns%initial_density)
        params_ns%Machnumber= sqrt(params_ns%initial_velocity(1)**2 &
                                  +params_ns%initial_velocity(2)**2 &
                                  +params_ns%initial_velocity(3)**2)/params_ns%c0
        params_ns%Reynolds  = params_ns%initial_density*params_ns%domain_size(2)* &
                              params_ns%machnumber*params_ns%c0/params_ns%mu0
      endif

      if (params_ns%mpirank==0) then
        write(*,*)
        write(*,*)
        write(*,*) "Additional Information"
        write(*,'(" -----------------------")')
        nx_max = (params_ns%Bs(1)-1) * 2**(params_ns%Jmax)
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
      subroutine check_parameters()
          implicit none

          if ( params_ns%dim==3 ) then
            if( params_ns%n_eqn<5 ) call abort(9898,'Please increase number of data fields (min 5)')
            if( min(pF,UxF,UyF,UzF,rhoF)<0 )  call abort(9898,'Check names of data fields [p,Ux,Uy,Uz,rho]!')
          else
            if( params_ns%n_eqn<4 ) call abort(9898,'Please increase number of data fields (min 4)')
            if( min(pF,UxF,UyF,rhoF)<0 )      call abort(9898,'Check names of data fields [p,Ux,Uy,rho]!')
          end if


    end subroutine check_parameters


    !> This function allocates a statevector like quantite.
    !> Hence a 2D(3D) array with 4(5) data fields
    subroutine allocate_statevector_ns(data,Bs,g)
        implicit none
        !data=(sqrt(rho),sqrt(rho)u,sqrt(rho)v,sqrt(rho)w,p )
        real(kind=rk), allocatable, intent(inout) :: data(:,:,:,:)
        integer(kind=ik),              intent(in) :: g
        integer(kind=ik), dimension(3), intent(in) :: Bs

        if (.not.allocated(data)) then
          if (params_ns%dim==3) then
            allocate(data(1:Bs(1)+2*g, 1:Bs(2)+2*g, 1:Bs(3)+2*g, params_ns%n_eqn))
          else
            allocate(data(1:Bs(1)+2*g, 1:Bs(2)+2*g, 1, params_ns%n_eqn))
        endif
      endif
    end subroutine allocate_statevector_ns


    !> \brief convert the statevector \f$(\sqrt(\rho),\sqrt(\rho)u,\sqrt(\rho)v,p )\f$
    !! to the desired format.
    subroutine convert_statevector3D(phi,convert2format)

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



    end subroutine convert_statevector3D


        !> \brief convert the statevector \f$(\sqrt(\rho),\sqrt(\rho)u,\sqrt(\rho)v,p )\f$
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
              converted_vector(:,:,rhoF)  =phi(:,:,rhoF)**2
              ! rho u
              converted_vector(:,:,UxF)   =phi(:,:,UxF)*phi(:,:,rhoF)
              ! rho v
              converted_vector(:,:,UyF)   =phi(:,:,UyF)*phi(:,:,rhoF)

                ! kinetic energie
              converted_vector(:,:,pF)=phi(:,:,UxF)**2+phi(:,:,UyF)**2
              converted_vector(:,:,pF)=converted_vector(:,:,pF)*0.5_rk
              ! total energie
              ! e_tot=e_kin+p/(gamma-1)
              converted_vector(:,:,pF)=converted_vector(:,:,pF)+phi(:,:,pF)/(params_ns%gamma_-1)
            case ("pure_variables")
              ! add ambient pressure
              !rho
              converted_vector(:,:,rhoF)= phi(:,:,rhoF)**2
              !u
              converted_vector(:,:,UxF)= phi(:,:, UxF)/phi(:,:,rhoF)
              !v
              converted_vector(:,:,UyF)= phi(:,:, UyF)/phi(:,:,rhoF)
              !p
              converted_vector(:,:,pF)= phi(:,:, pF)
            case default
                call abort(7771,"the format is unkown: "//trim(adjustl(convert2format)))
            end select

            phi=converted_vector

        end subroutine convert_statevector2D


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
          call convert_statevector3D(phi_out,format_out)
        endif
    end subroutine convert2format






end module module_navier_stokes_params
