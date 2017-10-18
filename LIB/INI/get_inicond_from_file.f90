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
    !> number of heavy and light active blocks
    integer(kind=ik), intent(inout)       :: hvy_n, lgt_n
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

    ! read treecode, time, iteration and total number of blocks from first input file
    call read_mesh_and_attributes(params%input_files(1), params, lgt_n, hvy_n, lgt_block, time, iteration)

    ! read datafields from files into hvy_block array
    do dF = 1, N_files
        call read_field(params%input_files(dF), dF, params, hvy_block, hvy_n )
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
