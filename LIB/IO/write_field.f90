! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: write_field.f90
! version: 0.5
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
! 26/01/17 - switch to 3D, v0.5
!          - add dirs_3D array for 3D neighbor codes
!
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
    real(kind=rk), intent(in)           :: hvy_block(:, :, :, :, :)
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
    character(len=3)                    :: dirs_2D(16)
    character(len=7)                    :: dirs_3D(74)

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank, lgt_rank

    ! loop variable
    integer(kind=ik)                    :: k, l, i, hvy_id
    ! grid parameter
    integer(kind=ik)                    :: Bs, g

    ! array for coordinate values, need 6 values for 3D (4 for 2D)
    real(kind=rk)                       :: coords_origin(3), coords_spacing(3)

!---------------------------------------------------------------------------------------------
! variables initialization

    ! set MPI parameters
    rank            = params%rank

    Bs = params%number_block_nodes
    g  = params%number_ghost_nodes

    dirs_2D = (/'__N', '__E', '__S', '__W', '_NE', '_NW', '_SE', '_SW', 'NNE', 'NNW', 'SSE', 'SSW', 'ENE', 'ESE', 'WNW', 'WSW'/)

    dirs_3D = (/'__1/___', '__2/___', '__3/___', '__4/___', '__5/___', '__6/___', '_12/___', '_13/___', '_14/___', '_15/___', &
               '_62/___', '_63/___', '_64/___', '_65/___', '_23/___', '_25/___', '_43/___', '_45/___', '123/___', '134/___', &
               '145/___', '152/___', '623/___', '634/___', '645/___', '652/___', '__1/123', '__1/134', '__1/145', '__1/152', &
               '__2/123', '__2/623', '__2/152', '__2/652', '__3/123', '__3/623', '__3/134', '__3/634', '__4/134', '__4/634', &
               '__4/145', '__4/645', '__5/145', '__5/645', '__5/152', '__5/652', '__6/623', '__6/634', '__6/645', '__6/652', &
               '_12/123', '_12/152', '_13/123', '_13/134', '_14/134', '_14/145', '_15/145', '_15/152', '_62/623', '_62/652', &
               '_63/623', '_63/634', '_64/634', '_64/645', '_65/645', '_65/652', '_23/123', '_23/623', '_25/152', '_25/652', &
               '_43/134', '_43/634', '_45/145', '_45/645' /)

!---------------------------------------------------------------------------------------------
! main body

    ! file name depends on variable names
    select case(params%physics_type)
        case('2D_convection_diffusion')
            ! select corresponding datafield name
            write( fname,'(a, "_", i12.12, ".h5")') trim(adjustl(params%physics%names(dF-1))), nint(time * 1.0e6_rk)

        case('2D_navier_stokes')
            ! select corresponding datafield name
            write( fname,'(a, "_", i12.12, ".h5")') trim(adjustl(params%physics_ns%names(dF-1))), nint(time * 1.0e6_rk)

        case('3D_convection_diffusion')
            ! select corresponding datafield name
            write( fname,'(a, "_", i12.12, ".h5")') trim(adjustl(params%physics%names(dF-1))), nint(time * 1.0e6_rk)

        case('3D_navier_stokes')
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
            if ( params%threeD_case ) then
                ! 3D:
                call write_field_hdf5_3D( fname, dsetname, hvy_block( g+1:Bs+g, g+1:Bs+g, g+1:Bs+g, dF, hvy_id), .false.)
            else
                ! 2D:
                call write_field_hdf5( fname, dsetname, hvy_block( g+1:Bs+g, g+1:Bs+g, 1, dF, hvy_id), .false.)
            end if


            ! add useful attributes to the block:
            ! write treecode
            call write_attribute( fname, dsetname, "treecode", lgt_block(lgt_active(k), 1:lgt_block(lgt_active(k), params%max_treelevel+1) ) )

            call write_attribute( fname, dsetname, "time", (/time/))
            call write_attribute( fname, dsetname, "iteration", (/iteration/))
            call write_attribute( fname, dsetname, "rank", (/rank/))

            ! save coordinates
            call write_attribute( fname, dsetname, "coord_x", hvy_block( 1, 1:Bs, 1, 1, hvy_id) )
            call write_attribute( fname, dsetname, "coord_y", hvy_block( 2, 1:Bs, 1, 1, hvy_id) )
            call write_attribute( fname, dsetname, "coord_z", hvy_block( 3, 1:Bs, 1, 1, hvy_id) )

            ! save neighbors
            if ( params%threeD_case ) then
                ! 3D:
                ! loop over all neighbors
                do i = 1, 74
                    ! neighbor exists
                    if ( hvy_neighbor( hvy_id, i ) /= -1 ) then
                        call write_attribute( fname, dsetname, dirs_3D(i),  lgt_block( hvy_neighbor( hvy_id, i ), 1:lgt_block( hvy_neighbor( hvy_id, i ), params%max_treelevel+1) ) )
                    end if
                end do

            else
                ! 2D:
                ! loop over all neighbors
                do i = 1, 16
                    ! neighbor exists
                    if ( hvy_neighbor( hvy_id, i ) /= -1 ) then
                        call write_attribute( fname, dsetname, dirs_2D(i),  lgt_block( hvy_neighbor( hvy_id, i ), 1:lgt_block( hvy_neighbor( hvy_id, i ), params%max_treelevel+1) ) )
                    end if
                end do

            end if

            ! write coordinates (origin and spacing) in new dataset
            ! dataset name
            write(dsetname,'("block_",i8.8,"_origin")') k
            ! coordinate values
            coords_origin(1) = hvy_block( 1, 1, 1, 1, hvy_id)
            coords_origin(2) = hvy_block( 2, 1, 1, 1, hvy_id)
            coords_origin(3) = hvy_block( 3, 1, 1, 1, hvy_id)
            ! write data
            if ( params%threeD_case ) then
                ! 3D:
                call write_field_hdf5_1D( fname, dsetname, coords_origin(:), .false.)
            else
                ! 2D:
                call write_field_hdf5_1D( fname, dsetname, coords_origin(1:2), .false.)
            end if

            write(dsetname,'("block_",i8.8,"_spacing")') k
            ! coordinate values
            coords_spacing(1) = abs(hvy_block( 1, 1, 1, 1, hvy_id) - hvy_block( 1, 2, 1, 1, hvy_id) )
            coords_spacing(2) = abs(hvy_block( 2, 1, 1, 1, hvy_id) - hvy_block( 2, 2, 1, 1, hvy_id) )
            coords_spacing(3) = abs(hvy_block( 3, 1, 1, 1, hvy_id) - hvy_block( 3, 2, 1, 1, hvy_id) )
            ! write data
            if ( params%threeD_case ) then
                ! 3D:
                call write_field_hdf5_1D( fname, dsetname, coords_spacing(:), .false.)
            else
                ! 2D:
                call write_field_hdf5_1D( fname, dsetname, coords_spacing(1:2), .false.)
            end if

        end if

        ! synchronize procs, to avoid simultaneous file operations
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

    end do

end subroutine write_field
