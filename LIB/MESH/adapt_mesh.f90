! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: adapt_mesh.f90
! version: 0.4
! author: msr
!
! mesh adapting main function
!
! input:    - params, light and heavy data
! output:   - light and heavy data arrays
!
! = log ======================================================================================
!
! 10/11/16 - switch to v0.4
! ********************************************************************************************

subroutine adapt_mesh( params, block_list, block_data, neighbor_list )

!---------------------------------------------------------------------------------------------
! modules

    use mpi
    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(inout)   :: params
    ! light data array
    integer(kind=ik), intent(inout)     :: block_list(:, :)
    ! heavy data array - block data
    real(kind=rk), intent(inout)        :: block_data(:, :, :, :)
    ! neighbor list
    integer(kind=ik), intent(inout)     :: neighbor_list(:)

    ! loop variables
    integer(kind=ik)                    :: k

    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0, sub_t1

!---------------------------------------------------------------------------------------------
! interfaces

    interface
        subroutine threshold_block( params, block_list, block_data, neighbor_list )
            use module_params
            type (type_params), intent(in)              :: params
            integer(kind=ik), intent(inout)             :: block_list(:, :)
            real(kind=rk), intent(inout)                :: block_data(:, :, :, :)
            integer(kind=ik), intent(in)                :: neighbor_list(:)
        end subroutine threshold_block

        subroutine ensure_gradedness( block_list, neighbor_list, N, max_treelevel )
            use module_params
            integer(kind=ik), intent(inout)             :: block_list(:, :)
            integer(kind=ik), intent(in)                :: neighbor_list(:)
            integer(kind=ik), intent(in)                :: N, max_treelevel
        end subroutine ensure_gradedness

        subroutine ensure_completeness( block_list, max_treelevel )
            use module_params
            integer(kind=ik), intent(inout)             :: block_list(:, :)
            integer(kind=ik), intent(in)                :: max_treelevel
        end subroutine ensure_completeness

        subroutine coarse_mesh( params, block_list, block_data )
            use module_params
            type (type_params), intent(in)              :: params
            integer(kind=ik), intent(inout)             :: block_list(:, :)
            real(kind=rk), intent(inout)                :: block_data(:, :, :, :)
        end subroutine coarse_mesh

        subroutine update_neighbors(block_list, neighbor_list, N, max_treelevel)
            use module_params
            integer(kind=ik), intent(in)                :: block_list(:, :)
            integer(kind=ik), intent(out)               :: neighbor_list(:)
            integer(kind=ik), intent(in)                :: N
            integer(kind=ik), intent(in)                :: max_treelevel
        end subroutine update_neighbors

        subroutine balance_load( params, block_list, block_data, neighbor_list )
            use module_params
            type (type_params), intent(in)              :: params
            integer(kind=ik), intent(inout)             :: block_list(:, :)
            real(kind=rk), intent(inout)                :: block_data(:, :, :, :)
            integer(kind=ik), intent(inout)               :: neighbor_list(:)
        end subroutine balance_load

    end interface

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    ! maximal number of loops to coarsen the mesh == one block go down from max_treelevel to min_treelevel
    do k = 1, (params%max_treelevel - params%min_treelevel)

        ! start time
        sub_t0 = MPI_Wtime()
        ! check where to coarsen (refinement done with safety zone)
        call threshold_block( params, block_list, block_data, neighbor_list )
        ! end time
        sub_t1 = MPI_Wtime()
        ! save time diff
        params%comp_time(4 + (k-1)*5) = sub_t1 - sub_t0

        ! start time
        sub_t0 = MPI_Wtime()
        ! unmark blocks that cannot be coarsened due to gradedness
        call ensure_gradedness( block_list, neighbor_list, params%number_blocks, params%max_treelevel )
        ! end time
        sub_t1 = MPI_Wtime()
        ! save time diff
        params%comp_time(5 + (k-1)*5) = sub_t1 - sub_t0

        ! start time
        sub_t0 = MPI_Wtime()
        ! ensure completeness
        call ensure_completeness( block_list, params%max_treelevel )
        ! end time
        sub_t1 = MPI_Wtime()
        ! save time diff
        params%comp_time(6 + (k-1)*5) = sub_t1 - sub_t0

        ! start time
        sub_t0 = MPI_Wtime()
        ! adapt the mesh
        call coarse_mesh( params, block_list, block_data )
        ! end time
        sub_t1 = MPI_Wtime()
        ! save time diff
        params%comp_time(7 + (k-1)*5) = sub_t1 - sub_t0

        ! start time
        sub_t0 = MPI_Wtime()
        ! update neighbor relations
        call update_neighbors( block_list, neighbor_list, params%number_blocks, params%max_treelevel )
        ! end time
        sub_t1 = MPI_Wtime()
        ! save time diff
        params%comp_time(8 + (k-1)*5) = sub_t1 - sub_t0

    end do

    ! start time
    sub_t0 = MPI_Wtime()
    ! balance load
    call balance_load( params, block_list, block_data, neighbor_list )
    ! end time
    sub_t1 = MPI_Wtime()
    ! save time diff
    params%comp_time(24) = sub_t1 - sub_t0

    ! start time
    sub_t0 = MPI_Wtime()
    ! update neighbor relations
    call update_neighbors( block_list, neighbor_list, params%number_blocks, params%max_treelevel )
    ! end time
    sub_t1 = MPI_Wtime()
    ! save time diff
    params%comp_time(25) = sub_t1 - sub_t0

end subroutine adapt_mesh
