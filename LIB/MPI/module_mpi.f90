!>\dir
!> MPI synchronization routines

!> \file
!> \brief File of MPI Module


!> \brief This module implements all MPI synchronization routines
!> \details
!>          * synchronization of ghost nodes
!>          * synchronization of light data
!>          * creation of rank block rank lists
!>          * copy, write send and recieve buffers
!> \version 0.6
!> \author engels, reiss, msr
!! 12/01/17 - create
!! 18/07/2018 - remove all old ghost node routines, build in JR's improved version of MSR's new ghost nodes
module module_MPI

!---------------------------------------------------------------------------------------------
! modules

    use mpi
    ! global parameters
    use module_params
    ! debug module
    use module_debug
    ! interpolation routines
    use module_interpolation

    implicit none
!---------------------------------------------------------------------------------------------
! variables

    SAVE

    ! send/receive buffer, integer and real
    ! allocate in init substep not in synchronize subroutine, to avoid slow down when using
    ! large numbers of processes and blocks per process
    integer(kind=ik), allocatable :: int_send_buffer(:,:), int_receive_buffer(:,:)
    real(kind=rk), allocatable    :: real_send_buffer(:,:), real_receive_buffer(:,:)

    logical :: ghost_nodes_module_ready = .false.

    ! I usually find it helpful to use the private keyword by itself initially, which specifies
    ! that everything within the module is private unless explicitly marked public.
    PRIVATE

    ! ! global variables
    ! data_buffer
    ! res_pre_data
    ! com_matrix
    ! int_pos
    ! data_bounds_names

    PUBLIC :: sync_ghosts, blocks_per_mpirank, synchronize_lgt_data, reset_ghost_nodes
    PUBLIC :: check_redundant_nodes, synchronize_ghosts_generic_sequence


!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

contains

    include "synchronize_ghosts.f90"
    include "blocks_per_mpirank.f90"
    include "synchronize_lgt_data.f90"
    include "reset_ghost_nodes.f90"
    include "check_redundant_nodes.f90"



subroutine init_ghost_nodes( params )
    implicit none
    !> user defined parameter structure
    type (type_params), intent(in) :: params
    ! local variables
    integer(kind=ik) :: buffer_N_int, buffer_N, Bs, g, N_dF, number_blocks, number_procs, rank


    if (.not. ghost_nodes_module_ready) then
        number_blocks   = params%number_blocks
        Bs              = params%number_block_nodes
        g               = params%number_ghost_nodes
        N_dF            = params%number_data_fields
        rank            = params%rank

        ! synchronize buffer length
        ! assume: all blocks are used, all blocks have external neighbors,
        ! max neighbor number: 2D = 12, 3D = 56
        ! max neighborhood size, 2D: (Bs+g+1)*(g+1)
        ! max neighborhood size, 3D: (Bs+g+1)*(g+1)*(g+1)
        if ( params%threeD_case ) then
            buffer_N = number_blocks * 56 * (Bs+g+1)*(g+1)*(g+1) * N_dF
            buffer_N_int = number_blocks * 56 * 3
        else
            buffer_N = number_blocks * 12 * (Bs+g+1)*(g+1) * N_dF
            buffer_N_int = number_blocks * 12 * 3
        end if

        ! allocate synch buffer
        allocate( int_send_buffer( buffer_N_int, params%number_procs) )
        if (rank==0) write(*,'("GHOSTS-INIT: Allocated ",A," shape=",7(i9,1x))') "int_send_buffer", shape(int_send_buffer)

        allocate( int_receive_buffer( buffer_N_int, params%number_procs) )
        if (rank==0) write(*,'("GHOSTS-INIT: Allocated ",A," shape=",7(i9,1x))') "int_receive_buffer", shape(int_receive_buffer)

        allocate( real_send_buffer( buffer_N, params%number_procs) )
        if (rank==0) write(*,'("GHOSTS-INIT: Allocated ",A," shape=",7(i9,1x))') "real_send_buffer", shape(real_send_buffer)

        allocate( real_receive_buffer( buffer_N, params%number_procs) )
        if (rank==0) write(*,'("GHOSTS-INIT: Allocated ",A," shape=",7(i9,1x))') "real_receive_buffer", shape(real_receive_buffer)

        if (rank==0) then
            write(*,'("GHOSTS-INIT: Real buffer size is",g15.3," GB ")') 2.0_rk*size(real_send_buffer)*8.0_rk/1000.0_rk/1000.0_rk/1000.0_rk
            write(*,'("GHOSTS-INIT: Int  buffer size is",g15.3," GB ")') 2.0_rk*size(int_send_buffer)*8.0_rk/1000.0_rk/1000.0_rk/1000.0_rk
        endif

        ghost_nodes_module_ready = .true.
    endif

end subroutine

end module module_MPI
