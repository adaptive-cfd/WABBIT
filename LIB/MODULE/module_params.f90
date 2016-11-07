! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: module_params.f90
! version: 0.4
! author: msr
!
! params data structure, define all constant parameters for global use
!
! todo: module actually only works for specific RHS, split between RHS parameters and program
!       parameters in future versions
!
! = log ======================================================================================
!
! 04/11/16 - switch to v0.4, merge old block_params structure with new structure
! ********************************************************************************************

module module_params

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! define data precision parameters
    integer, parameter, public   :: sngl_prec=selected_real_kind(4)
    integer, parameter, public   :: dble_prec=selected_real_kind(8)

    integer, parameter, public   :: int_prec=selected_int_kind(8)

    integer, parameter, public   :: rk=dble_prec
    integer, parameter, public   :: ik=int_prec

    ! global user defined data structure for time independent variables
    type type_params

        ! maximal time for main time loop
        real(kind=rk)                               :: time_max
        ! CFL criteria for time step calculation
        real(kind=rk)                               :: CFL
        ! data writing frequency
        integer(kind=ik)                            :: write_freq

        ! threshold for wavelet indicator
        real(kind=rk)                               :: eps
        ! minimal level for blocks in data tree
        integer(kind=ik)                            :: min_treelevel
        ! maximal level for blocks in data tree
        integer(kind=ik)                            :: max_treelevel

        ! order of refinement predictor
        character(len=80)                           :: order_predictor
        ! order of spatial discretization
        character(len=80)                           :: order_discretization
        ! boundary condition
        character(len=80)                           :: boundary_cond
        ! initial condition
        character(len=80)                           :: initial_cond

        ! domain length
        real(kind=rk)                               :: Lx, Ly
        ! grid parameter
        integer(kind=ik)                            :: number_domain_nodes
        integer(kind=ik)                            :: number_block_nodes
        integer(kind=ik)                            :: number_ghost_nodes

        ! switch for mesh adaption
        logical                                     :: adapt_mesh

        ! number of allocated heavy data fields per process
        integer(kind=ik)                            :: number_blocks
        ! number of allocated data fields in heavy data array, number of fields in heavy data (depend from time step scheme, ...)
        integer(kind=ik)                            :: number_data_fields
        integer(kind=ik)                            :: number_fields

        ! each data field can use separat velocity and diffusion coefficient (in 2D convection-diffusion RHS)
        ! stored all velocity components in one 1D vector, same with the diffusion coefficients
        real(kind=rk), dimension(:), allocatable    :: u0
        real(kind=rk), dimension(:), allocatable    :: nu

        ! block distribution for load balancing (also used for start distribution)
        character(len=80)                           :: block_distribution

    end type type_params

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

end module module_params
