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
      logical                                     :: save_filter_strength

  end type type_params_filter   
  !=========================================================



  !=========================================================
  !               NAVIER STOKES PARAMETER
  !---------------------------------------------------------
  type :: type_params_ns
        ! Courant-Friedrichs-Lewy
        real(kind=rk)                               :: CFL, T_end
        ! spatial domain
        real(kind=rk)                               :: Lx, Ly, Lz
        ! number data fields
        integer(kind=ik)                            :: number_data_fields
        ! dimension
        integer(kind=ik)                            :: dim, N_fields_saved
        ! adiabatic coefficient
        real(kind=rk)                               :: gamma_
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
        ! width of initialcond
        real(kind=rk)                               :: inicond_width
        ! dissipation switch
        logical                                     :: dissipation
        ! variable names
        character(len=80), allocatable              :: names(:)
        ! înitial condition
        character(len=80)                           :: inicond
        ! type used for the spatial discretization
        character(len=80)                           :: discretization
        ! ------------------------------------------------------------------------------------------
        ! penalization
        logical                                     :: penalization,smooth_mask=.True., sponge_layer
        ! penalization parameter and sponge parameter
        real(kind=rk)                               :: C_eta,C_sp
        ! geometry to display
        character(len=80)                           :: geometry="cylinder"
        ! geometric parameters for cylinder (x_0,r)
        real(kind=rk)                               :: x_cntr(1:3), R_cyl
        ! mean variables on domain
        real(kind=rk)                               :: mean_density
        real(kind=rk)                               :: mean_pressure
        ! we need to know which mpirank prints output..
        integer(kind=ik)                            :: mpirank, mpisize
        !---------------------------------------------------------
        ! filter structure containing all parameters of the filter
        type(type_params_filter)                    :: filter
  end type type_params_ns

  ! statevector index
  integer(kind=ik) ,save,public :: rhoF,UxF,UyF,UzF,pF

contains

  include "initial_conditions.f90"
  

  subroutine init_navier_stokes_eq(params_ns, FILE )
    implicit none
    !> pointer to inifile
    type(inifile) ,intent(inout)     :: FILE
    !> params structure of navier stokes
    type(type_params_ns),intent(inout)  :: params_ns
    
    if (params_ns%mpirank==0) then
      write(*,*)
      write(*,*)
      write(*,*) "PARAMS: navier stokes equations!"
      write(*,'(" -----------------------------------")')
    endif
    ! read number_data_fields
    call read_param_mpi(FILE, 'Blocks', 'number_data_fields', params_ns%number_data_fields, 1 )
    ! dimension
    call read_param_mpi(FILE, 'Dimensionality', 'dim', params_ns%dim, 2 )
    ! spatial domain size
    call read_param_mpi(FILE, 'DomainSize', 'Lx', params_ns%Lx, 1.0_rk )
    call read_param_mpi(FILE, 'DomainSize', 'Ly', params_ns%Ly, 1.0_rk )
    call read_param_mpi(FILE, 'DomainSize', 'Lz', params_ns%Lz, 0.0_rk )
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

  


subroutine init_other_params(params_ns, FILE )
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

    call read_param_mpi(FILE, 'Blocks', 'number_data_fields', params_ns%number_data_fields, 1 )
    call read_param_mpi(FILE, 'Discretization', 'order_discretization', params_ns%discretization, "FD_2nd_central")

    call read_param_mpi(FILE, 'Time', 'CFL', params_ns%CFL, 1.0_rk   )
    call read_param_mpi(FILE, 'Time', 'time_max', params_ns%T_end, 1.0_rk   )
   
    
  end subroutine init_other_params





end module module_navier_stokes_params
