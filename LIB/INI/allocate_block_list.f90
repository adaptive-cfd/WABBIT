! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: allocate_block_list.f90
! version: 0.4
! author: msr
!
! allocate and reset block list (global light data)
!
! light data array  -> line number = ( 1 + proc_rank ) * heavy_data_line_number
!                   -> column(1:max_treelevel): block treecode, treecode -1 => block is inactive
!                   -> column(max_treelevel+1): treecode length = mesh level
!                   -> column(max_treelevel+2):   refinement status (-1..coarsen / 0...no change / +1...refine)
!
! input:    - params
! output:   - empty light data array
!
! = log ======================================================================================
!
! 04/11/16 - switch to v0.4
! 26/01/17 - use process rank and number of procs from params struct
!
! ********************************************************************************************

subroutine  allocate_block_list( params, lgt_block )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(in)              :: params

    ! light data array
    integer(kind=ik), allocatable, intent(out)  :: lgt_block(:, :)

    ! number of light and heavy data
    integer(kind=ik)                            :: number_blocks
    ! maximal treelevel = maximal length of treecode
    integer(kind=ik)                            :: max_treelevel

    ! allocation error variable
    integer(kind=ik)                            :: allocate_error

    ! process rank
    integer(kind=ik)                            :: rank
    ! number of processes
    integer(kind=ik)                            :: number_procs

!---------------------------------------------------------------------------------------------
! variables initialization

    ! set MPI parameters
    rank         = params%rank
    number_procs = params%number_procs

    ! set parameters for readability
    number_blocks = params%number_blocks
    max_treelevel = params%max_treelevel

!---------------------------------------------------------------------------------------------
! main body

    ! allocate memory
    allocate( lgt_block( number_procs*number_blocks, max_treelevel+2), stat=allocate_error )
    call check_allocation(allocate_error)

    ! reset data:
    ! all blocks are inactive, reset treecode
    lgt_block(:, 1:max_treelevel) = -1
    ! all blocks are inactive, reset treelevel
    lgt_block(:, max_treelevel+1) = -1
    ! set refinement to 0
    lgt_block(:, max_treelevel+2) = 0

    ! output
    if (rank==0) then
        if ( params%unit_test .eqv. .false. ) then
            write(*,'(80("_"))')
            write(*,'("INIT: System is allocating light data for ",i7," blocks" )') number_procs*number_blocks
        end if
    end if

end subroutine allocate_block_list
