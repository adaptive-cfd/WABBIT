! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: init_data.f90
! version: 0.5
! author: msr
!
! This routine initializes the block data, i.e. it evaluates the initial condition on the grid
!
! input:    - parameter array
!           - light data array
!           - heavy data array
!           - neighbor data array
!           - light and heavy active block list
! output:   - filled user defined data structure for global params
!           - initialized light and heavy data arrays
!
! = log ======================================================================================
!
! 04/11/16 - switch to v0.4, now run complete initialization within these subroutine and return
!            initialized block data to main program
! 07/12/16 - now uses heavy work data array
! 25/01/17 - switch to 3D, v0.5
!
! ********************************************************************************************
subroutine allocate_grid(params, lgt_block, hvy_block, hvy_work, hvy_neighbor, lgt_active, hvy_active)

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(inout)               :: params
    ! light data array
    integer(kind=ik), allocatable, intent(out)      :: lgt_block(:, :)
    ! heavy data array - block data
    real(kind=rk), allocatable, intent(out)         :: hvy_block(:, :, :, :, :)
    ! heavy work array  )
    real(kind=rk), allocatable, intent(out)         :: hvy_work(:, :, :, :, :)
    ! neighbor array (heavy data)
    integer(kind=ik), allocatable, intent(out)      :: hvy_neighbor(:,:)
    ! list of active blocks (light data)
    integer(kind=ik), allocatable, intent(out)      :: lgt_active(:)
    ! list of active blocks (light data)
    integer(kind=ik), allocatable, intent(out)      :: hvy_active(:)
    ! allocation error variabel
    integer(kind=ik)                                :: allocate_error
    ! local shortcuts:
    integer(kind=ik)                                :: Bs,g,dF,number_blocks, rank, number_procs

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
        allocate( hvy_block( Bs+2*g, Bs+2*g, Bs+2*g, dF+1, number_blocks ), stat=allocate_error )
        call check_allocation(allocate_error)

        ! work data (Runge-Kutta substeps and old time level)
        allocate( hvy_work( Bs+2*g, Bs+2*g, Bs+2*g, dF*5, number_blocks ), stat=allocate_error )
        call check_allocation(allocate_error)

        ! 3D: maximal 74 neighbors per block
        allocate( hvy_neighbor( params%number_blocks, 74 ), stat=allocate_error )
        call check_allocation(allocate_error)
    else
        ! 2D:
        ! datafields + 1 -> first field for coordinates, ...
        allocate( hvy_block( Bs+2*g, Bs+2*g, 1, dF+1, number_blocks ), stat=allocate_error )
        call check_allocation(allocate_error)

        ! work data (Runge-Kutta substeps and old time level)
        allocate( hvy_work( Bs+2*g, Bs+2*g, 1, dF*5, number_blocks ), stat=allocate_error )
        call check_allocation(allocate_error)

        ! 2D: maximal 16 neighbors per block
        allocate( hvy_neighbor( params%number_blocks, 16 ), stat=allocate_error )
        call check_allocation(allocate_error)
    end if

    ! allocate memory
    allocate( lgt_block( number_procs*number_blocks, params%max_treelevel+2), stat=allocate_error )
    call check_allocation(allocate_error)

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
    allocate( lgt_active( size(lgt_block, 1) ), stat=allocate_error )
    call check_allocation(allocate_error)

    ! note: 5th dimension in heavy data is block id
    allocate( hvy_active( size(hvy_block, 5) ), stat=allocate_error )
    call check_allocation(allocate_error)

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
