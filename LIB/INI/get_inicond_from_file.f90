!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name get_inicond_from_file.f90
!> \version 0.5
!> \author engels, sm
!
!> \brief read initial condition from a file
!
!>
!! input:
!!           - parameter array
!!
!! output:
!!           - light data array
!!           - heavy data array
!!           - number of active blocks (light and heavy)
!!           - light and heavy active block list
!!           -
!!
!!
!! = log ======================================================================================
!! \n
!! 28/09/17 - create
!
! ********************************************************************************************

subroutine get_inicond_from_file(params, lgt_block, hvy_block, hvy_n, lgt_n, time, iteration)
!---------------------------------------------------------------------------------------------
! modules

    use module_IO

!---------------------------------------------------------------------------------------------
! variables

    implicit none
    !> user defined parameter structure
    type (type_params), intent(in)        :: params
    !> light data array
    integer(kind=ik), intent(inout)       :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)          :: hvy_block(:, :, :, :, :)
    !> list of active blocks (light data)
!    integer(kind=ik), intent(inout)      :: hvy_active(:)
    !> list of active blocks light data)
!    integer(kind=ik), intent(inout)      :: lgt_active(:)
    !> number of heavy and light active blocks
    integer(kind=ik), intent(inout)       :: hvy_n, lgt_n
    !> sorted list of numerical treecodes, used for block finding
!    integer(kind=tsize), intent(inout)   :: lgt_sortednumlist(:,:)
    !> time loop variables
    real(kind=rk), intent(inout)          :: time
    integer(kind=ik), intent(inout)       :: iteration


    ! cpu time variables for running time calculation
    real(kind=rk)                         :: sub_t0, sub_t1
    ! number of files to read from
    integer(kind=ik)                      :: N_files
    ! loop variable
    integer(kind=ik)                      :: k, dF
!---------------------------------------------------------------------------------------------
! variables initialization

    ! number of files to read from
    N_files = params%number_data_fields

!---------------------------------------------------------------------------------------------
! main body

    ! start time
    sub_t0 = MPI_wtime()

    ! call fetch_attributes(params%ini_files(1), params, lgt_n, hvy_n, time, iteration)

    call read_mesh_and_attributes(params%inicond_files(1), params, lgt_n, hvy_n, lgt_block, time, iteration)

!    !create hvy_active and lgt_active list
!    call create_active_and_sorted_lists(params, lgt_block, lgt_active, lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true.)
 
!    ! update neighbor relations
!    call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n )

    do dF = 1, N_files
        call read_field(params%inicond_files(dF), dF, params, hvy_block, lgt_n, hvy_n )
    end do


    ! end time
    sub_t1 = MPI_wtime()

    ! write time
    if ( params%debug ) then
        ! find free or corresponding line
        k = 1
        do while ( debug%name_comp_time(k) /= "---" )
            ! entry for current subroutine exists
            if ( debug%name_comp_time(k) == "read_data" ) exit
            k = k + 1
        end do
        ! write time
        debug%name_comp_time(k) = "read_data"
        debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1
        debug%comp_time(k, 2)   = debug%comp_time(k, 2) + sub_t1 - sub_t0
    end if

end subroutine get_inicond_from_file
