!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name init_data.f90
!> \version 0.5
!> \author msr
!
!> \brief Allocate grid data (light, heavy, neighbors, active lists etc), initialize
!
!>
!! input:
!!           - parameter array
!!           - light data array
!!           - heavy data array
!!           - neighbor data array
!!           - light and heavy active block list
!!
!! output:
!!           - filled user defined data structure for global params
!!           - initialized light and heavy data arrays
!!
!! = log ======================================================================================
!! \n
!! 04/11/16 - switch to v0.4, now run complete initialization within these subroutine and return
!!            initialized block data to main program \n
!! 07/12/16 - now uses heavy work data array \n
!! 25/01/17 - switch to 3D, v0.5
!
! ********************************************************************************************
subroutine allocate_grid(params, lgt_block, hvy_block, hvy_work, hvy_neighbor, lgt_active, hvy_active, lgt_sortednumlist, int_send_buffer, int_receive_buffer, real_send_buffer, real_receive_buffer)

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(inout)               :: params
    !> light data array
    integer(kind=ik), allocatable, intent(out)      :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), allocatable, intent(out)         :: hvy_block(:, :, :, :, :)
    !> heavy work array  )
    real(kind=rk), allocatable, intent(out)         :: hvy_work(:, :, :, :, :)
    !> neighbor array (heavy data)
    integer(kind=ik), allocatable, intent(out)      :: hvy_neighbor(:,:)
    !> list of active blocks (light data)
    integer(kind=ik), allocatable, intent(out)      :: lgt_active(:)
    !> list of active blocks (light data)
    integer(kind=ik), allocatable, intent(out)      :: hvy_active(:)
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), allocatable, intent(out)   :: lgt_sortednumlist(:,:)
    ! local shortcuts:
    integer(kind=ik)                                :: Bs, g, dF, number_blocks, rank, number_procs, buffer_N

    ! send/receive buffer, integer and real
    integer(kind=ik), allocatable, intent(out)      :: int_send_buffer(:,:), int_receive_buffer(:,:)
    real(kind=rk), allocatable, intent(out)         :: real_send_buffer(:,:), real_receive_buffer(:,:)

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization
    ! set parameters for readability
    rank         = params%rank
    number_blocks   = params%number_blocks
    Bs              = params%number_block_nodes
    g               = params%number_ghost_nodes
    dF              = params%number_data_fields
    number_procs    = params%number_procs

    ! synchronize buffer length
    ! assume: all blocks are used, all blocks have external neighbors,
    ! max neighbor number: 2D = 12, 3D = 56
    ! max neighborhood size, 2D: (Bs+g+1)*(g+1)
    ! max neighborhood size, 3D: (Bs+g+1)*(g+1)*(g+1)
    if ( params%threeD_case ) then
        buffer_N = number_blocks * 56 * (Bs+g+1)*(g+1)*(g+1) * dF
    else
        buffer_N = number_blocks * 12 * (Bs+g+1)*(g+1) * dF
    end if

!---------------------------------------------------------------------------------------------
! main body

    if (rank == 0) then
      write(*,'(80("_"))')
      write(*,'(A)') "INIT: Beginning memory allocation and initialization."
      write(*,'("INIT: mpisize=",i6)') params%number_procs
      write(*,'("INIT: Bs=",i7," blocks-per-rank=",i7," total blocks=", i7)') Bs, number_blocks, number_blocks*number_procs
    endif

    ! allocate memory
    if ( params%threeD_case ) then
        ! 3D:
        ! datafields + 1 -> first field for coordinates, ...
        allocate( hvy_block( Bs+2*g, Bs+2*g, Bs+2*g, dF+1, number_blocks ) )
        ! work data (Runge-Kutta substeps and old time level)
        allocate( hvy_work( Bs+2*g, Bs+2*g, Bs+2*g, dF*5, number_blocks ) )
        ! 3D: maximal 74 neighbors per block
        allocate( hvy_neighbor( params%number_blocks, 74 ) )
    else
        ! 2D:
        ! datafields + 1 -> first field for coordinates, ...
        allocate( hvy_block( Bs+2*g, Bs+2*g, 1, dF+1, number_blocks ) )
        ! work data (Runge-Kutta substeps and old time level)
        allocate( hvy_work( Bs+2*g, Bs+2*g, 1, dF*5, number_blocks ) )
        ! 2D: maximal 16 neighbors per block
        allocate( hvy_neighbor( params%number_blocks, 16 ) )
    end if

    ! allocate memory
    allocate( lgt_block( number_procs*number_blocks, params%max_treelevel+2) )
    allocate( lgt_sortednumlist( size(lgt_block,1), 2) )

    ! allocate synch buffer
    !allocate( int_send_buffer( number_blocks*3+1, number_procs) )
    !allocate( int_receive_buffer( number_blocks*3+1, number_procs) )
    allocate( int_send_buffer( 1000, number_procs) )
    allocate( int_receive_buffer( 1000, number_procs) )
    allocate( real_send_buffer( buffer_N, number_procs) )
    allocate( real_receive_buffer( buffer_N, number_procs) )

    ! reset data:
    ! all blocks are inactive, reset treecode
    lgt_block(:, 1:params%max_treelevel) = -1
    ! all blocks are inactive, reset treelevel
    lgt_block(:, params%max_treelevel+1) = -1
    ! set refinement to 0
    lgt_block(:, params%max_treelevel+2) = 0


    ! reset data
    hvy_block = 9.99e99_rk
    hvy_work = 9.99e99_rk
    hvy_neighbor = -1

    ! allocate active list
    allocate( lgt_active( size(lgt_block, 1) ) )

    ! note: 5th dimension in heavy data is block id
    allocate( hvy_active( size(hvy_block, 5) ) )

    if (rank == 0) then
      ! note we currently use 8byte per real and integer by default, so all the same bytes per point
      write(*,'("INIT: Local memory footprint is ",g15.3,"GB per mpirank")') &
      (size(hvy_block)+size(hvy_work)+size(lgt_block)+size(hvy_neighbor)+size(lgt_active)+size(hvy_active))*8.0_rk/1000.0_rk/1000.0_rk/1000.0_rk
      write(*,'("INIT: TOTAL memory footprint is ",g15.3,"GB")') &
      (size(hvy_block)+size(hvy_work)+size(lgt_block)+size(hvy_neighbor)+size(lgt_active)+size(hvy_active))*8.0_rk*real(number_procs,kind=rk)/1000.0_rk/1000.0_rk/1000.0_rk

      write(*,'("INIT: System is allocating heavy data for ",i7," blocks and ", i3, " fields" )') number_blocks, dF
      write(*,'("INIT: System is allocating light data for ",i7," blocks" )') number_procs*number_blocks
      write(*,'("INIT: System is allocating heavy work data for ",i7," blocks " )') number_blocks
    endif


end subroutine allocate_grid
