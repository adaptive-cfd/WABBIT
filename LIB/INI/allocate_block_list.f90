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
!                   -> column(max_treelevel+1):   refinement status (-1..coarsen / 0...no change / +1...refine)
!
! input:    - maximal number of blocks per process
!           - maximal treelevel
! output:   - empty light data array
!
! = log ======================================================================================
!
! 04/11/16 - switch to v0.4
! ********************************************************************************************

subroutine  allocate_block_list(block_list, number_blocks, max_treelevel)

!---------------------------------------------------------------------------------------------
! modules

    use mpi
    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! light data array
    integer(kind=ik), allocatable, intent(out)  :: block_list(:, :)
    ! number of light and heavy data
    integer(kind=ik), intent(in)                :: number_blocks
    ! maximal treelevel = maximal length of treecode
    integer(kind=ik), intent(in)                :: max_treelevel

    ! allocation error variable
    integer(kind=ik)                            :: allocate_error

    ! MPI error variable
    integer(kind=ik)                            :: ierr
    ! process rank
    integer(kind=ik)                            :: rank
    ! number of processes
    integer(kind=ik)                            :: number_procs

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    ! determinate process number
    call MPI_Comm_size(MPI_COMM_WORLD, number_procs, ierr)
    ! determinate process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    ! allocate memory
    allocate( block_list( number_procs*number_blocks, max_treelevel+1), stat=allocate_error )

    ! reset data
    !
    ! all blocks are inactive, reset treecode
    block_list(:, 1:max_treelevel) = -1
    ! set refinement to 0
    block_list(:, max_treelevel+1) = 0

    ! output
    if (rank==0) then
        write(*,'(80("_"))')
        write(*,'("INIT: System is allocating light data for ",i7," blocks" )') number_procs*number_blocks
    end if

end subroutine allocate_block_list
