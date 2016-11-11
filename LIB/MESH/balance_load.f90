! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: balance_load.f90
! version: 0.4
! author: msr
!
! balance the load
!
! input:    - params, light and heavy data
! output:   - light and heavy data arrays
!
! = log ======================================================================================
!
! 08/11/16 - switch to v0.4
! ********************************************************************************************

subroutine balance_load( params, block_list, block_data )

!---------------------------------------------------------------------------------------------
! modules

    use mpi
    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(in)      :: params
    ! light data array
    integer(kind=ik), intent(inout)     :: block_list(:, :)
    ! heavy data array - block data
    real(kind=rk), intent(inout)        :: block_data(:, :, :, :)

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank
    ! number of processes
    integer(kind=ik)                    :: number_procs

    ! distribution type
    character(len=80)                   :: distribution

    ! block distribution lists
    integer(kind=ik), allocatable       :: opt_dist_list(:), dist_list(:)

    ! allocation error variable
    integer(kind=ik)                    :: allocate_error

    ! loop variables
    integer(kind=ik)                    :: k, N, num_blocks

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    distribution    = params%block_distribution

    ! determinate process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    ! determinate process number
    call MPI_Comm_size(MPI_COMM_WORLD, number_procs, ierr)

    ! allocate block to proc lists
    allocate( opt_dist_list( number_procs ), stat=allocate_error )
    allocate( dist_list( number_procs ), stat=allocate_error )

    ! number of light data
    N = size(block_list, 1)

    ! reset block count
    num_blocks = 0

    ! reset block to proc list
    dist_list = 0

!---------------------------------------------------------------------------------------------
! main body

    select case(distribution)
        case("equal")
            ! simple uniformly distribution

            !---------------------------------------------------------------------------------
            ! first: calculate optimal block distribution
            ! count number of active blocks and current block distribution
            do k = 1, N
                ! block is active
                if ( block_list(k, 1) /= -1 ) then
                    ! count block for proc corresponing to light data
                    dist_list( (k-1) / params%number_blocks + 1 ) = dist_list( (k-1) / params%number_blocks + 1 ) + 1
                end if
            end do
            ! count blocks
            do k = 1, number_procs
                num_blocks = num_blocks + dist_list(k)
            end do

            ! optimal distribution
            opt_dist_list = (num_blocks - mod(num_blocks, number_procs))/number_procs
            ! distribute remaining blocks
            if (mod(num_blocks, number_procs) > 0) then
                opt_dist_list(1:mod(num_blocks, number_procs)) = (num_blocks - mod(num_blocks, number_procs))/number_procs + 1
            end if

            print*, opt_dist_list - dist_list
            stop


        case("sfc1")

        case default
            write(*,'(80("_"))')
            write(*,*) "ERROR: block distribution scheme is unknown"
            write(*,*) distribution
            stop

    end select

    ! clean up
    deallocate( opt_dist_list, stat=allocate_error )
    deallocate( dist_list, stat=allocate_error )

end subroutine balance_load
