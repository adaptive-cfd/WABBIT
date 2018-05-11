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
subroutine allocate_grid(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, hvy_active, &
    lgt_sortednumlist, simulation, hvy_work, hvy_synch, int_send_buffer, int_receive_buffer, &
    real_send_buffer, real_receive_buffer)

    !---------------------------------------------------------------------------------------------
    ! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(inout)                   :: params
    !> light data array
    integer(kind=ik), allocatable, intent(out)          :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), allocatable, intent(out)             :: hvy_block(:, :, :, :, :)
    !> heavy work array
    real(kind=rk), allocatable, optional, intent(out)   :: hvy_work(:, :, :, :, :)
    integer(kind=1), allocatable, optional, intent(out) :: hvy_synch(:, :, :, :)
    !> neighbor array (heavy data)
    integer(kind=ik), allocatable, intent(out)          :: hvy_neighbor(:,:)
    !> list of active blocks (light data)
    integer(kind=ik), allocatable, intent(out)          :: lgt_active(:)
    !> list of active blocks (light data)
    integer(kind=ik), allocatable, intent(out)          :: hvy_active(:)
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), allocatable, intent(out)       :: lgt_sortednumlist(:,:)
    ! local shortcuts:
    integer(kind=ik)                                    :: Bs, g, N_dF, number_blocks,&
                                                      rank, number_procs, buffer_N, buffer_N_int

    !> send/receive buffer, integer and real
    integer(kind=ik), allocatable, optional, intent(out) :: int_send_buffer(:,:), int_receive_buffer(:,:)
    real(kind=rk), allocatable, optional, intent(out)    :: real_send_buffer(:,:), real_receive_buffer(:,:)
    !> do we have to allocate everything?
    logical, intent(in)                                  :: simulation
    integer(kind=ik)                                     :: rk_steps

    !---------------------------------------------------------------------------------------------
    ! interfaces

    !---------------------------------------------------------------------------------------------
    ! variables initialization
    ! set parameters for readability
    rank         = params%rank
    number_blocks   = params%number_blocks
    Bs              = params%number_block_nodes
    g               = params%number_ghost_nodes
    N_dF            = params%number_data_fields
    number_procs    = params%number_procs

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

    !---------------------------------------------------------------------------------------------
    ! main body


    if (rank == 0) then
        write(*,'(80("_"))')
        write(*,'(A)') "INIT: Beginning memory allocation and initialization."
        write(*,'("INIT: mpisize=",i6)') params%number_procs
        write(*,'("INIT: Bs=",i7," blocks-per-rank=",i7," total blocks=", i7)') Bs, number_blocks, number_blocks*number_procs
    endif

    rk_steps = max(size(params%butcher_tableau,1)-1,params%N_fields_saved)
    ! allocate memory
    if ( params%threeD_case ) then
        ! 3D:
        if (rank == 0) write(*,'("INIT: Allocating a 3D case.")')
        ! datafields
        allocate( hvy_block( Bs+2*g, Bs+2*g, Bs+2*g, N_dF, number_blocks ) )
        if (rank==0) write(*,'("INIT: Allocated ",A," shape=",7(i9,1x))') "hvy_block", shape(hvy_block)


        ! work data (Runge-Kutta substeps and old time level)
        if (simulation) then
            allocate( hvy_work( Bs+2*g, Bs+2*g, Bs+2*g, N_dF*(rk_steps+1), number_blocks ) )
            if (rank==0) write(*,'("INIT: Allocated ",A," shape=",7(i9,1x))') "hvy_work", shape(hvy_work)

            ! synch array, use for ghost nodes synchronization
            allocate( hvy_synch( Bs+2*g, Bs+2*g, Bs+2*g, number_blocks ) )
            if (rank==0) write(*,'("INIT: Allocated ",A," shape=",7(i9,1x))') "hvy_synch", shape(hvy_synch)
        end if

        ! 3D: maximal 74 neighbors per block
        allocate( hvy_neighbor( params%number_blocks, 74 ) )
        if (rank==0) write(*,'("INIT: Allocated ",A," shape=",7(i9,1x))') "hvy_neighbor", shape(hvy_neighbor)
    else
        ! 2D:
        if (rank==0) write(*,'("INIT: Allocating a 2D case.")')
        ! datafields
        allocate( hvy_block( Bs+2*g, Bs+2*g, 1, N_dF, number_blocks ) )
        if (rank==0) write(*,'("INIT: Allocated ",A," shape=",7(i9,1x))') "hvy_block", shape(hvy_block)

        ! work data (Runge-Kutta substeps and old time level)
        if (simulation) then
            allocate( hvy_work( Bs+2*g, Bs+2*g, 1, N_dF*(rk_steps+1), number_blocks ) )
            if (rank==0) write(*,'("INIT: Allocated ",A," shape=",5(i9,1x))') "hvy_work", shape(hvy_work)

            ! synch array, use for ghost nodes synchronization
            allocate( hvy_synch( Bs+2*g, Bs+2*g, 1, number_blocks ) )
            if (rank==0) write(*,'("INIT: Allocated ",A," shape=",7(i9,1x))') "hvy_synch", shape(hvy_synch)
        end if

        ! 2D: maximal 16 neighbors per block
        allocate( hvy_neighbor( params%number_blocks, 16 ) )
        if (rank==0) write(*,'("INIT: Allocated ",A," shape=",7(i9,1x))') "hvy_neighbor", shape(hvy_neighbor)
    end if

    ! allocate memory
    allocate( lgt_block( number_procs*number_blocks, params%max_treelevel+2) )
    if (rank==0) write(*,'("INIT: Allocated ",A," shape=",7(i9,1x))') "lgt_block", shape(lgt_block)

    allocate( lgt_sortednumlist( size(lgt_block,1), 2) )
    if (rank==0) write(*,'("INIT: Allocated ",A," shape=",7(i9,1x))') "lgt_sortednumlist", shape(lgt_sortednumlist)

    if (simulation) then
        ! allocate synch buffer
        allocate( int_send_buffer( buffer_N_int, number_procs) )
        if (rank==0) write(*,'("INIT: Allocated ",A," shape=",7(i9,1x))') "int_send_buffer", shape(int_send_buffer)

        allocate( int_receive_buffer( buffer_N_int, number_procs) )
        if (rank==0) write(*,'("INIT: Allocated ",A," shape=",7(i9,1x))') "int_receive_buffer", shape(int_receive_buffer)

        allocate( real_send_buffer( buffer_N, number_procs) )
        if (rank==0) write(*,'("INIT: Allocated ",A," shape=",7(i9,1x))') "real_send_buffer", shape(real_send_buffer)

        allocate( real_receive_buffer( buffer_N, number_procs) )
        if (rank==0) write(*,'("INIT: Allocated ",A," shape=",7(i9,1x))') "real_receive_buffer", shape(real_receive_buffer)
    end if
    ! reset data:
    ! all blocks are inactive, reset treecode
    lgt_block(:, 1:params%max_treelevel) = -1
    ! all blocks are inactive, reset treelevel
    lgt_block(:, params%max_treelevel+1) = -1
    ! set refinement to 0
    lgt_block(:, params%max_treelevel+2) = 0


    ! ! reset data
    ! hvy_block = 9.99e99_rk
    ! if (simulation) then
    !     hvy_work = 9.99e99_rk
    !     hvy_synch = -99
    ! end if
    ! hvy_neighbor = -1

    ! allocate active list
    allocate( lgt_active( size(lgt_block, 1) ) )
    if (rank==0) write(*,'("INIT: Allocated ",A," shape=",7(i9,1x))') "lgt_active", shape(lgt_active)

    ! note: 5th dimension in heavy data is block id
    allocate( hvy_active( size(hvy_block, 5) ) )
    if (rank==0) write(*,'("INIT: Allocated ",A," shape=",7(i9,1x))') "hvy_active", shape(hvy_active)

    if (rank == 0 .and. simulation) then
        write(*,'("INIT: System is allocating heavy data for ",i7," blocks and ", i3, " fields" )') number_blocks, N_dF
        write(*,'("INIT: System is allocating light data for ",i7," blocks" )') number_procs*number_blocks
        write(*,'("INIT: System is allocating heavy work data for ",i7," blocks " )') number_blocks
        write(*,'("INIT: Real buffer size is",g15.3," GB ")') 2.0_rk*size(real_send_buffer)*8.0_rk/1000.0_rk/1000.0_rk/1000.0_rk
        write(*,'("INIT: Int  buffer size is",g15.3," GB ")') 2.0_rk*size(int_send_buffer)*8.0_rk/1000.0_rk/1000.0_rk/1000.0_rk

        ! note we currently use 8byte per real and integer by default, so all the same bytes per point
        write(*,'("INIT: Measured (true) local (on 1 cpu) memory footprint is ",g15.3,"GB per mpirank")') &
        (dble(size(hvy_block)) + dble(size(hvy_work)) + dble(size(lgt_block)) + dble(size(lgt_sortednumlist)) &
        + dble(size(hvy_neighbor)) + dble(size(lgt_active)) + dble(size(hvy_active)) + dble(size(hvy_synch))/8.0 &
        + dble(size(real_send_buffer)) + dble(size(real_receive_buffer)) + dble(size(int_send_buffer)) &
        + dble(size(int_receive_buffer)))*8.0_rk/1000.0_rk/1000.0_rk/1000.0_rk

        write(*,'("INIT: Measured (true) TOTAL (on all CPU) memory footprint is ",g15.3,"GB")') &
        ((dble(size(hvy_block)) + dble(size(hvy_work)) + dble(size(lgt_block)) + dble(size(lgt_sortednumlist)) &
        + dble(size(hvy_neighbor)) + dble(size(lgt_active)) + dble(size(hvy_active)) + dble(size(hvy_synch))/8.0 &
        + dble(size(real_send_buffer)) + dble(size(real_receive_buffer)) + dble(size(int_send_buffer)) &
        + dble(size(int_receive_buffer)))*8.0_rk/1000.0_rk/1000.0_rk/1000.0_rk)*dble(params%number_procs)
    end if


end subroutine allocate_grid





subroutine deallocate_grid(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, hvy_active, &
    lgt_sortednumlist, hvy_work, hvy_synch, int_send_buffer, &
    int_receive_buffer, real_send_buffer, real_receive_buffer)

    !---------------------------------------------------------------------------------------------
    ! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(inout)                   :: params
    !> light data array
    integer(kind=ik), allocatable, intent(out)          :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), allocatable, intent(out)             :: hvy_block(:, :, :, :, :)
    !> heavy work array  )
    real(kind=rk), allocatable, optional, intent(out)   :: hvy_work(:, :, :, :, :)
    integer(kind=1), allocatable, optional, intent(out) :: hvy_synch(:, :, :, :)
    !> neighbor array (heavy data)
    integer(kind=ik), allocatable, intent(out)          :: hvy_neighbor(:,:)
    !> list of active blocks (light data)
    integer(kind=ik), allocatable, intent(out)          :: lgt_active(:)
    !> list of active blocks (light data)
    integer(kind=ik), allocatable, intent(out)          :: hvy_active(:)
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), allocatable, intent(out)       :: lgt_sortednumlist(:,:)
    !> send/receive buffer, integer and real
    integer(kind=ik), allocatable, optional, intent(out) :: int_send_buffer(:,:), int_receive_buffer(:,:)
    real(kind=rk), allocatable, optional, intent(out)    :: real_send_buffer(:,:), real_receive_buffer(:,:)

    if (params%rank == 0) then
        write(*,'(80("-"))')
        write(*,'(A)') "FREE: Beginning freeying of memory."
    endif

    if (allocated(hvy_block)) deallocate( hvy_block )
    if (allocated(hvy_work)) deallocate( hvy_work )
    if (allocated(hvy_synch)) deallocate( hvy_synch )
    if (allocated(hvy_neighbor)) deallocate( hvy_neighbor )
    if (allocated(lgt_block)) deallocate( lgt_block )
    if (allocated(lgt_sortednumlist)) deallocate( lgt_sortednumlist )
    if (allocated(int_send_buffer)) deallocate( int_send_buffer )
    if (allocated(int_receive_buffer)) deallocate( int_receive_buffer )
    if (allocated(real_send_buffer)) deallocate( real_send_buffer )
    if (allocated(real_receive_buffer)) deallocate( real_receive_buffer )
    if (allocated(lgt_active)) deallocate( lgt_active )
    if (allocated(hvy_active)) deallocate( hvy_active )

    if (params%rank == 0) then
        write(*,'(A)') "All memory is cleared!"
        write(*,'(80("-"))')
    endif

end subroutine deallocate_grid
