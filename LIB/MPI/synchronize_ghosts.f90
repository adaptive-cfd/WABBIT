!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name synchronize_ghosts.f90
!> \version 0.5
!> \author msr
!
!> \brief synchronize ghosts nodes
!! synchronization is done in 4 stages:
!!          stage 1: finer blocks give data to coarser blocks (restriction)
!!          stage 2: blocks on same level synch
!!          stage 3: coarser blocks synch with finer blocks (prediction)
!!          stage 4: additional step to correct interpolation errors (needed for stability)
!!
!! subroutine structure:
!! ---------------------
!!      1. create external communication list:
!!              loop over active heavy blocks, if neighbor needs ghost nodes and neighbor is located
!!              on other proc: new com list entry
!!      2. synchronize com matrix:
!!              needed for buffer allocation and communication order
!!      3. fill send buffer:
!!              sender proc loop over his com list and gather all data for all neighbor procs/blocks
!!      4. send/receive data
!!      5. write own and received ghost nodes into own blocks
!!              - loop over own active heavy data (own block has to receive data)
!!              - if neighbor exists: check level difference and synch stage
!!                note: here level diff is neighbor (sender) - me (receiver)
!!              - if something to synchronize: switch neighborhood, because, neighborhood is defined from
!!                senders point of view, but we loop actually over receiver blocks
!!              - if neighbor is internal: copy data, note: correct sender/receiver block ids!
!!              - external neighbor: search for neighbor data in received buffer
!!                note: reiceived buffer includes all data, and write routine can handle this, but
!!                at this point we work onm one specific block, maybe later we can switch back to old
!!                internal/external handling
!!  \todo rework ghost nodes writing routine to avoid receive buffer searching
!!  stage 4 handling:
!!      - external neighbor: check condition in com list creation - then send/receive data as before
!!      - internal neighbor: use new copy_redundant_nodes subroutine
!!  \todo if possible use same copy ghost nodes routine for all 4 stages
!!
!>
!! input:    - params, light and heavy data \n
!! output:   - heavy data array
!
!> \details
!! = log ======================================================================================
!! \n
!! 08/11/16 - switch to v0.4 \n
!! 06/01/17 - use RMA to synchronize data \n
!! 31/01/17 - switch to 3D, v0.5 \n
!! 12/04/17 - redundant ghost nodes workaround
!! 19/05/17 - switch to new synchronization routine (correct redundant nodes handling)
!! 16/06/17 - allocate all send/receive buffer in ini step -> huge performance boost
!
! ********************************************************************************************

subroutine sync_ghosts(  params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n, hvy_synch )
    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n
    integer(kind=1), intent(inout)      :: hvy_synch(:, :, :, :)

    logical :: sync
    character(len=80) :: method
    integer(kind=ik) :: count
    real(kind=rk) :: t0

    if (.not. ghost_nodes_module_ready) then
        call init_ghost_nodes( params )
    endif

    t0 = MPI_wtime()

    count = command_argument_count()
    call get_command_argument(count, method)

    select case (method)
    case ("--staging")
        ! MSR CODE
        sync = .true.
        call check_redundant_nodes( params, lgt_block, hvy_block, hvy_synch, hvy_neighbor, hvy_active, &
        hvy_n, sync, .false., .false. )
    case default
        ! JR version of MSR code above
        call synchronize_ghosts_generic_sequence( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )
    end select

    call toc( params, "WRAPPER: sync ghosts", MPI_wtime()-t0 )
end subroutine sync_ghosts
