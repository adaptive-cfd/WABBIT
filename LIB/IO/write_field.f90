! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: write_field.f90
! version: 0.4
! author: engels, msr
!
! write data of a single datafield dF at timestep iteration and time t
!
! input:    - time loop parameter
!           - datafield number
!           - parameter array
!           - light data array
!           - heavy data array
! output:   -
!
! = log ======================================================================================
!
! 07/11/16 - switch to v0.4
! ********************************************************************************************

subroutine write_field(time, iteration, dF, params, block_list, block_data, neighbor_list)

!---------------------------------------------------------------------------------------------
! modules

    use mpi
    ! global parameters
    use module_params
    ! hdf5 file wrapper
    use module_hdf5_wrapper

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! time loop parameters
    real(kind=rk), intent(in)           :: time
    integer(kind=ik), intent(in)        :: iteration

    ! datafield number
    integer(kind=ik), intent(in)        :: dF

    ! user defined parameter structure
    type (type_params), intent(in)      :: params
    ! light data array
    integer(kind=ik), intent(in)        :: block_list(:, :)
    ! heavy data array - block data
    real(kind=rk), intent(in)           :: block_data(:, :, :, :)
    ! heavy data array - neifghbor data
    integer(kind=ik), intent(in)        :: neighbor_list(:)

    ! file name
    character(len=80)                   :: fname
    ! set name
    character(len=80)                   :: dsetname
    ! neighbor name list
    character(len=3)                    :: dirs(16)

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank

    ! loop variable
    integer(kind=ik)                    :: k, l, i
    ! grid parameter
    integer(kind=ik)                    :: Bs, g

!---------------------------------------------------------------------------------------------
! variables initialization

    Bs = params%number_block_nodes
    g  = params%number_ghost_nodes

    dirs = (/'__N', '__E', '__S', '__W', '_NE', '_NW', '_SE', '_SW', 'NNE', 'NNW', 'SSE', 'SSW', 'ENE', 'ESE', 'WNW', 'WSW'/)

!---------------------------------------------------------------------------------------------
! main body

    ! determinate process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    ! fiel name
    write( fname,'("data_",i8.8,"_field_",i2.2,".h5")') nint(time * 1.0e4_rk), dF

    ! create the filename, dF == 2 is the first real datafield
    if ( (rank == 0) .and. (dF == 2) ) then

        write(*,'(80("_"))')
        write(*,'("IO: writing data for time = ", f15.8," file = ",A)') time, trim(adjustl(fname))

        ! overwrite the file, if it already exists
        call init_empty_file( fname )
    end if

    l = 1
    ! save block data, loop over all light data
    do k = 1, size(block_list, 1)

        ! calculate proc rank from light data line number
        if (k > l*params%number_blocks) then
            l = l + 1
        end if

        ! block is active, only corresponding proc can work with heavy data
        if ( (block_list(k, 1) /= -1) .and. ( (l-1) == rank) ) then

            ! the name of the block within th hdf5 file
            write(dsetname,'("block_",i8.8)') k

            ! write data field
            !call write_field_hdf5( fname, dsetname, block_data( g+1:Bs+g, g+1:Bs+g, dF, k - (l-1)*params%number_blocks), .false.)
            call write_field_hdf5( fname, dsetname, block_data( 1:Bs+2*g, 1:Bs+2*g, dF, k - (l-1)*params%number_blocks), .false.)

            ! add useful attributes to the block:
            ! write treecode
            call write_attribute( fname, dsetname, "treecode", block_list(k, 1:block_list(k, params%max_treelevel+1) ) )

            call write_attribute( fname, dsetname, "time", (/time/))
            call write_attribute( fname, dsetname, "iteration", (/iteration/))
            call write_attribute( fname, dsetname, "rank", (/rank/))

            ! save coordinates
            call write_attribute( fname, dsetname, "coord_x", block_data( 1, 1:Bs, 1, k - (l-1)*params%number_blocks))
            call write_attribute( fname, dsetname, "coord_y", block_data( 2, 1:Bs, 1, k - (l-1)*params%number_blocks))

            ! save neighbors
            ! loop over all neighbors
            do i = 1, 16
                ! neighbor exists
                if ( neighbor_list( (k - (l-1)*params%number_blocks - 1)*16 + i ) /= -1 ) then
                    call write_attribute( fname, dsetname, dirs(i),  block_list(neighbor_list( (k - (l-1)*params%number_blocks - 1)*16 + i ), 1:block_list(neighbor_list( (k - (l-1)*params%number_blocks - 1)*16 + i ), params%max_treelevel+1) ) )
                end if
            end do

        end if

        ! synchronize procs, to avoid simultaneous file operations
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    end do

end subroutine write_field
