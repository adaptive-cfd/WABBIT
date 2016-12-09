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

subroutine write_field(time, iteration, dF, params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n)

!---------------------------------------------------------------------------------------------
! modules

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
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    ! heavy data array - block data
    real(kind=rk), intent(in)           :: hvy_block(:, :, :, :)
    ! heavy data array - neifghbor data
    integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)
    ! list of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_active(:)
    ! number of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_n

    ! file name
    character(len=80)                   :: fname
    ! set name
    character(len=80)                   :: dsetname
    ! neighbor name list
    character(len=3)                    :: dirs(16)

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank, lgt_rank

    ! loop variable
    integer(kind=ik)                    :: k, l, i, hvy_id
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

    ! file name depends on variable names
    select case(params%physics_type)
        case('2D_convection_diffusion')
            ! select corresponding datafield name
            write( fname,'(a, "_", i12.12, ".h5")') trim(adjustl(params%physics%names(dF-1))), nint(time * 1.0e6_rk)

        case('2D_navier_stokes')
            ! select corresponding datafield name
            write( fname,'(a, "_", i12.12, ".h5")') trim(adjustl(params%physics_ns%names(dF-1))), nint(time * 1.0e6_rk)

    end select

    ! create the filename, dF == 2 is the first real datafield
    if ( (rank == 0) .and. (dF == 2) ) then

        write(*,'(80("_"))')
        write(*,'("IO: writing data for time = ", f15.8," file = ",A)') time, trim(adjustl(fname))

        ! overwrite the file, if it already exists
        call init_empty_file( fname )
    elseif (rank == 0) then
        ! write output
        write(*,'("IO: writing data for time = ", f15.8," file = ",A)') time, trim(adjustl(fname))

    end if

    l = 1
    ! save block data, loop over all active light data
    do k = 1, lgt_n

        ! calculate proc rank from light data line number
        call lgt_id_to_proc_rank( lgt_rank, lgt_active(k), params%number_blocks )

        ! calculate heavy block id corresponding to light id
        call lgt_id_to_hvy_id( hvy_id, lgt_active(k), rank, params%number_blocks )

        ! only proc corresponding to block writes data
        if ( lgt_rank == rank ) then

            ! the name of the block within th hdf5 file
            write(dsetname,'("block_",i8.8)') k

            ! write data field
            call write_field_hdf5( fname, dsetname, hvy_block( g+1:Bs+g, g+1:Bs+g, dF, hvy_id), .false.)

            ! add useful attributes to the block:
            ! write treecode
            call write_attribute( fname, dsetname, "treecode", lgt_block(lgt_active(k), 1:lgt_block(lgt_active(k), params%max_treelevel+1) ) )

            call write_attribute( fname, dsetname, "time", (/time/))
            call write_attribute( fname, dsetname, "iteration", (/iteration/))
            call write_attribute( fname, dsetname, "rank", (/rank/))

            ! save coordinates
            call write_attribute( fname, dsetname, "coord_x", hvy_block( 1, 1:Bs, 1, hvy_id) )
            call write_attribute( fname, dsetname, "coord_y", hvy_block( 2, 1:Bs, 1, hvy_id) )

            ! save neighbors
            ! loop over all neighbors
            do i = 1, 16
                ! neighbor exists
                if ( hvy_neighbor( hvy_id, i ) /= -1 ) then
                    call write_attribute( fname, dsetname, dirs(i),  lgt_block( hvy_neighbor( hvy_id, i ), 1:lgt_block( hvy_neighbor( hvy_id, i ), params%max_treelevel+1) ) )
                end if
            end do

        end if

        ! synchronize procs, to avoid simultaneous file operations
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    end do

end subroutine write_field
