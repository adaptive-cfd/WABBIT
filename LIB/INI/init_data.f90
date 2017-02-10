! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: init_data.f90
! version: 0.5
! author: msr
!
! initialize all data: read params from ini file, allocate memory, initialize starting condition
! and decompose start matrix into block data
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
subroutine init_data(params, lgt_block, hvy_block, hvy_work, hvy_neighbor, lgt_active, hvy_active)

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(inout)                 :: params

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

    ! inifile name
    character(len=80)                               :: filename

    ! allocation error variabel
    integer(kind=ik)                                :: allocate_error

    ! loop variable
    integer(kind=ik)                                :: k

    ! cpu time variables for running time calculation
    real(kind=rk)                                   :: sub_t0, sub_t1

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    ! start time
    sub_t0 = MPI_Wtime()

    ! get the second command line argument: ini-file name
    call get_command_argument(2, filename)

    ! read ini-file and save parameter
    call ini_file_to_params( params, filename)

    !***************************************************************************
    ! allocate light/heavy data, initialize start field and write block data
    !
    ! allocate block_list
    call allocate_block_list( params, lgt_block )
    ! allocate heavy data
    call allocate_block_data( params, hvy_block )
    ! allocate heavy work data
    call allocate_work_data( params, hvy_work )

    ! initial data field
    select case( params%initial_cond )
        case ("gauss-blob","gauss_blob")

            ! 2D: gaussblob, 3D: zero fields
            if ( params%threeD_case ) then
                ! 3D:
                call inicond_sphere( params, lgt_block, hvy_block )
            else
                ! 2D:
                call inicond_gauss_blob( params, lgt_block, hvy_block )
            end if

        case ("vorticity_filaments")
!            call inicond_vorticity_filaments(params, lgt_block, hvy_block)
!            ! set density and pressure
!            hvy_block( :, :, 2, : ) = 1.0_rk
!            hvy_block( :, :, 5, : ) = 1e5_rk

        case default
            write(*,'(80("_"))')
            write(*,*) "ERROR: initial condition is unknown"
            write(*,*) params%initial_cond
            stop

    end select

    ! allocate active list
    allocate( lgt_active( size(lgt_block, 1) ), stat=allocate_error )
    call check_allocation(allocate_error)

    ! note: 5th dimension in heavy data is block id
    allocate( hvy_active( size(hvy_block, 5) ), stat=allocate_error )
    call check_allocation(allocate_error)

    ! ------------------------------------------------------------------------------------------------------
    ! init neighbor data array
    ! 2D: maximal 16 neighbors per block
    ! 3D: maximal 74 neighbors per block
    if ( params%threeD_case ) then
        ! 3D:
        allocate( hvy_neighbor( params%number_blocks, 74 ), stat=allocate_error )
        call check_allocation(allocate_error)

    else
        ! 2D:
        allocate( hvy_neighbor( params%number_blocks, 16 ), stat=allocate_error )
        call check_allocation(allocate_error)

    end if

    ! reset neighbor data array
    hvy_neighbor = -1

    ! ------------------------------------------------------------------------------------------------------
    ! init debug data
    ! note: fix size of time measurements array
    if ( params%debug ) then

        ! allocate array for time measurements - data
        allocate( debug%comp_time( 20, 4 ), stat=allocate_error )
        call check_allocation(allocate_error)

        ! reset times
        debug%comp_time = 0.0_rk

        ! allocate array for time measurements - names
        allocate( debug%name_comp_time( 20 ), stat=allocate_error )
        call check_allocation(allocate_error)

        ! reset names
        debug%name_comp_time = "---"

    end if

    ! end time
    sub_t1 = MPI_Wtime()
    ! write time
    if ( params%debug ) then
        ! find free or corresponding line
        k = 1
        do while ( debug%name_comp_time(k) /= "---" )
            ! entry for current subroutine exists
            if ( debug%name_comp_time(k) == "init_data" ) exit
            k = k + 1
        end do
        ! write time
        debug%name_comp_time(k) = "init_data"
        debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1
        debug%comp_time(k, 2)   = debug%comp_time(k, 2) + sub_t1 - sub_t0

    end if

end subroutine init_data
