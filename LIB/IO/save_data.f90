! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: save_data.f90
! version: 0.4
! author: msr
!
! save data main function, call write data routine
!
! input:    - time loop parameter
!           - parameter array
!           - light data array
!           - heavy data array
! output:   -
!
! = log ======================================================================================
!
! 07/11/16 - switch to v0.4
! ********************************************************************************************

subroutine save_data(iteration, time, params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n)

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! time loop parameters
    real(kind=rk), intent(in)                       :: time
    integer(kind=ik), intent(in)                    :: iteration

    ! user defined parameter structure
    type (type_params), intent(in)                  :: params
    ! light data array
    integer(kind=ik), intent(in)                    :: lgt_block(:, :)
    ! heavy data array - block data
    real(kind=rk), intent(in)                       :: hvy_block(:, :, :, :)
    ! heavy data array - neifghbor data
    integer(kind=ik), intent(in)                    :: hvy_neighbor(:,:)
    ! list of active blocks (light data)
    integer(kind=ik), intent(in)                    :: lgt_active(:)
    ! number of active blocks (light data)
    integer(kind=ik), intent(in)                    :: lgt_n

    ! loop variable
    integer(kind=ik)                                :: k

    ! cpu time variables for running time calculation
    real(kind=rk)                                   :: sub_t0, sub_t1

    ! MPI error variable
    integer(kind=ik)                                :: ierr
    ! process rank
    integer(kind=ik)                                :: rank

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    ! start time
    sub_t0 = MPI_Wtime()

    ! determinate process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    ! real datafields start at datafield 2
    do k = 2, params%number_data_fields+1
        call write_field(time, iteration, k, params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n)
    end do

    ! end time
    sub_t1 = MPI_Wtime()
    ! write time
    if ( params%debug ) then
        ! find first free line
        k = 1
        do while ( debug%name_comp_time(k) /= "---" )
            k = k + 1
        end do
        ! write time
        debug%name_comp_time(k) = "save_data"
        debug%comp_time(rank+1, k) = sub_t1 - sub_t0
    end if

end subroutine save_data
