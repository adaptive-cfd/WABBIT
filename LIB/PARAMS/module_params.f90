!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name module_params.f90
!> \version 0.4
!> \author msr
!
!> \brief params data structure, define all constant parameters for global use
!
!> \todo module actually only works for specific RHS, split between RHS parameters and program
!!       parameters in future versions
!
!> \details
!> = log ======================================================================================
!! \n
!! 04/11/16 - switch to v0.4, merge old block_params structure with new structure
! ********************************************************************************************

module module_params

!---------------------------------------------------------------------------------------------
! modules
    use mpi
    ! ini file parser module
    use module_ini_files_parser_mpi
    ! MPI general bridge module
    use module_bridge
!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! global user defined data structure for time independent variables
    type type_params

        ! maximal time for main time loop
        real(kind=rk)                                :: time_max=0.0_rk
        ! CFL criteria for time step calculation
        real(kind=rk)                                :: CFL=0.0_rk
        ! dt
        real(kind=rk)                                :: dt_fixed=0.0_rk, dt_max=0.0_rk
        ! number of allowed time steps
        integer(kind=ik)                             :: nt=99999999, inicond_refinements=0
        ! data writing frequency
        integer(kind=ik)                             :: write_freq=0
        ! data writing frequency
        character(len=80)                            :: write_method="fixed_time"
        ! data writing frequency
        real(kind=rk)                                :: write_time=0.1_rk
        ! data next write time, store here the next time for output data
        real(kind=rk)                                :: next_write_time=0.0_rk

        ! butcher tableau containing coefficients for Runge-Kutta method
        real(kind=rk), dimension(:,:), allocatable   :: butcher_tableau

        ! threshold for wavelet indicator
        real(kind=rk)                                :: eps=0.0_rk
        ! minimal level for blocks in data tree
        integer(kind=ik)                             :: min_treelevel=0
        ! maximal level for blocks in data tree
        integer(kind=ik)                             :: max_treelevel=0

        ! order of refinement predictor
        character(len=80)                            :: order_predictor=""
        ! order of spatial discretization
        character(len=80)                            :: order_discretization=""
        ! boundary condition
        character(len=80)                            :: boundary_cond=""
        ! initial condition
        character(len=80)                            :: initial_cond=""
        ! files we want to read for inital cond.
        character(len=80), dimension(:), allocatable :: input_files

        ! grid parameter
        integer(kind=ik)                             :: number_domain_nodes=0
        integer(kind=ik)                             :: number_block_nodes=0
        integer(kind=ik)                             :: number_ghost_nodes=0

        ! switch for mesh adaption
        logical                                      :: adapt_mesh=.false., adapt_inicond=.false.

        ! number of allocated heavy data fields per process
        integer(kind=ik)                             :: number_blocks=0
        ! number of allocated data fields in heavy data array, number of fields in heavy work data (depend from time step scheme, ...)
        integer(kind=ik)                             :: number_data_fields=0
        integer(kind=ik)                             :: number_fields=0

        ! block distribution for load balancing (also used for start distribution)
        character(len=80)                            :: block_distribution=""

        ! debug flag
        logical                                      :: debug=.false.

        ! use non-uniform mesh correction
        logical                                      :: non_uniform_mesh_correction=.true.

        ! -------------------------------------------------------------------------------------
        ! physics
        ! -------------------------------------------------------------------------------------
        ! physics type
        character(len=80)                            :: physics_type=""

        ! domain length
        real(kind=rk)                                :: Lx=0.0_rk, Ly=0.0_rk, Lz=0.0_rk

        ! use third dimension
        logical                                     :: threeD_case=.false.
        integer(kind=ik)                            :: dim=2 ! can be 2 or 3

        ! -------------------------------------------------------------------------------------
        ! statistics
        ! -------------------------------------------------------------------------------------
        real(kind=rk)    :: tsave_stats=0.0_rk, next_stats_time=0.0_rk
        integer(kind=ik) :: nsave_stats=0

        ! -------------------------------------------------------------------------------------
        ! MPI
        ! -------------------------------------------------------------------------------------
        ! data exchange method
        character(len=80)                           :: mpi_data_exchange=""
        ! process rank
        integer(kind=ik)                            :: rank=-1
        ! number of processes
        integer(kind=ik)                            :: number_procs=-1
        ! WABBIT communicator
        integer(kind=ik)                            :: WABBIT_COMM

        ! -------------------------------------------------------------------------------------
        ! bridge
        ! -------------------------------------------------------------------------------------
        ! bridge for connecting WABBIT to outdoor MPI_WORLD
        type(bridgeMPI)                             :: bridge
        !
        logical                                     :: bridge_exists = .false.
        !--------------------------------------------------------------------------------------
               !! particle connection
        !--------------------------------------------------------------------------------------
        !! - command to use for the particle program (over bridge)
        character(len=100)               :: particleCommand=""
        !! - Usage of a common myWorld_comm
        logical                          :: bridgeCommonMPI
        !! - Consideration of the fluid side as master in case of several myWorld_comms
        logical                          :: bridgeFluidMaster



        ! -------------------------------------------------------------------------------------
        ! saving
        ! -------------------------------------------------------------------------------------
        integer(kind=ik) :: N_fields_saved=0
        character(len=80), allocatable, dimension(:) :: field_names


        ! -------------------------------------------------------------------------------------
        ! unit test
        ! -------------------------------------------------------------------------------------
        ! unit test params struct
        logical                                     :: unit_test=.false.
        ! unit test time_stepper convergence flag
        logical                                     :: test_time_stepper=.false.
        ! unit test spatial convergence flag
        logical                                     :: test_spatial=.false.
        ! unit test wavelet compression
        logical                                     :: test_wavelet_comp=.false.
        ! unit test treecode flag
        logical                                     :: test_treecode=.false., test_ghost_nodes_synch=.false.

        ! -------------------------------------------------------------------------------------
        ! filter
        ! -------------------------------------------------------------------------------------
        ! type
        character(len=80)                           :: filter_type="no_filter"
        ! frequency
        integer(kind=ik)                            :: filter_freq=-1
        ! save filter strength sigma
        logical                                     :: save_filter_strength

    end type type_params

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body
contains

    ! this file reads the ini file and distributes all the parameters to the
    ! various structs holding them
    include "ini_file_to_params.f90"

    ! --------------------------------------------------------------------
    ! currently, wabbit is called ./wabbit 2D params.ini so the decision if we're
    ! running 2d or 3d is done in the command line call. here we figure that out
    ! and save the result in the parameter structure.
    ! --------------------------------------------------------------------
    subroutine decide_if_running_2D_or_3D(params)
      implicit none
      !> user defined parameter structure
      type (type_params), intent(inout) :: params

      character(len=80) :: dim_number

      ! read number of dimensions from command line
      call get_command_argument(1, dim_number)

      ! output dimension number
      if (params%rank==0) then
          write(*,'(80("_"))')
          write(*, '("INIT: running ", a3, " case")') dim_number
      end if

      ! save case dimension in params struct
      select case(dim_number)
          case('2D')
              params%threeD_case = .false.
              params%dim = 2
          case('3D')
              params%threeD_case = .true.
              params%dim = 3
          case('--help')
          case('--h')
          case default
              call abort(1,"ERROR: case dimension is wrong")
      end select

    end subroutine

end module module_params
