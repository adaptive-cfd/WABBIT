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

subroutine adapt_mesh( params, block_list, block_data, neighbor_list, debug )

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
    ! user defined parameter structure
    type (type_debug), intent(inout)    :: debug

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

    ! active_block_list for testing
    integer(kind=ik)                    :: lgt_active( size(block_list,1) )
    ! number of active blocks (light data)
    integer(kind=ik)                    :: n_lgt

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank

!---------------------------------------------------------------------------------------------
! interfaces

    interface
        subroutine threshold_block( params, block_list, block_data, neighbor_list, debug, adapt_count )
            use module_params
            type (type_params), intent(inout)           :: params
            integer(kind=ik), intent(inout)             :: block_list(:, :)
            real(kind=rk), intent(inout)                :: block_data(:, :, :, :)
            integer(kind=ik), intent(in)                :: neighbor_list(:)
            type (type_debug), intent(inout)            :: debug
            integer(kind=ik), intent(in)                :: adapt_count

        end subroutine threshold_block

        subroutine ensure_gradedness( block_list, neighbor_list, N, max_treelevel, lgt_active, n_lgt )
            use module_params
            integer(kind=ik), intent(inout)             :: block_list(:, :)
            integer(kind=ik), intent(in)                :: neighbor_list(:)
            integer(kind=ik), intent(in)                :: N, max_treelevel
            integer(kind=ik), intent(in)                :: lgt_active(:)
            integer(kind=ik), intent(in)                :: n_lgt
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
            integer(kind=ik), intent(inout)             :: neighbor_list(:)
        end subroutine balance_load

        subroutine create_lgt_active_list( block_list, lgt_active, n_lgt )
            use module_params
            integer(kind=ik), intent(in)                :: block_list(:, :)
            integer(kind=ik), intent(out)               :: lgt_active(:)
            integer(kind=ik), intent(out)               :: n_lgt
        end subroutine create_lgt_active_list

    end interface

!---------------------------------------------------------------------------------------------
! variables initialization

    ! determinate process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

!---------------------------------------------------------------------------------------------
! main body

    ! maximal number of loops to coarsen the mesh == one block go down from max_treelevel to min_treelevel
    do k = 1, (params%max_treelevel - params%min_treelevel)

        ! start time
        sub_t0 = MPI_Wtime()

        ! new light active list
        call create_lgt_active_list( block_list, lgt_active, n_lgt )

        ! end time
        sub_t1 = MPI_Wtime()
        if ( params%debug ) then
            debug%name_comp_time(5 + (k-1)* 6) = "lgt_active_list"
            debug%comp_time(rank+1, 5 + (k-1)* 6) = sub_t1 - sub_t0
        end if

        ! check where to coarsen (refinement done with safety zone)
        call threshold_block( params, block_list, block_data, neighbor_list, debug, k )

        ! start time
        sub_t0 = MPI_Wtime()

        ! unmark blocks that cannot be coarsened due to gradedness
        call ensure_gradedness( block_list, neighbor_list, params%number_blocks, params%max_treelevel, lgt_active, n_lgt )

        ! end time
        sub_t1 = MPI_Wtime()
        if ( params%debug ) then
            debug%name_comp_time(8 + (k-1)* 6) = "ensure_gradedness"
            debug%comp_time(rank+1, 8 + (k-1)* 6) = sub_t1 - sub_t0
        end if

        ! start time
        sub_t0 = MPI_Wtime()

        ! ensure completeness
        call ensure_completeness( block_list, params%max_treelevel )

        ! end time
        sub_t1 = MPI_Wtime()
        if ( params%debug ) then
            debug%name_comp_time(9 + (k-1)* 6) = "ensure_completeness"
            debug%comp_time(rank+1, 9 + (k-1)* 6) = sub_t1 - sub_t0
        end if

        ! start time
        sub_t0 = MPI_Wtime()

        ! adapt the mesh
        call coarse_mesh( params, block_list, block_data )

        ! end time
        sub_t1 = MPI_Wtime()
        if ( params%debug ) then
            debug%name_comp_time(10 + (k-1)* 6) = "coarse_mesh"
            debug%comp_time(rank+1, 10 + (k-1)* 6) = sub_t1 - sub_t0
        end if

        ! update neighbor relations
        call update_neighbors( block_list, neighbor_list, params%number_blocks, params%max_treelevel )

    end do

    ! start time
    sub_t0 = MPI_Wtime()

    ! balance load
    call balance_load( params, block_list, block_data, neighbor_list )

    ! end time
    sub_t1 = MPI_Wtime()
    if ( params%debug ) then
        debug%name_comp_time(29) = "balance_load"
        debug%comp_time(rank+1, 29) = sub_t1 - sub_t0
    end if

    ! start time
    sub_t0 = MPI_Wtime()

    ! update neighbor relations
    call update_neighbors( block_list, neighbor_list, params%number_blocks, params%max_treelevel )

    ! end time
    sub_t1 = MPI_Wtime()
    if ( params%debug ) then
        debug%name_comp_time(30) = "update_neighbors"
        debug%comp_time(rank+1, 30) = sub_t1 - sub_t0
    end if

end subroutine adapt_mesh
