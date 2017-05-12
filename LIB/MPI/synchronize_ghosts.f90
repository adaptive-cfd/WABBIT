!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name synchronize_ghosts.f90
!> \version 0.5
!> \author msr
!
!> \brief synchronize ghosts nodes
!
!> \details
!! input:    - params, light and heavy data \n
!! output:   - heavy data array
!
!> \todo change soubroutine, to work only on one datafield, not on all to the same time
!
! -------------------------------------------------------------------------------------------------------------------------
!> dirs = (/'__N', '__E', '__S', '__W', '_NE', '_NW', '_SE', '_SW', 'NNE', 'NNW', 'SSE', 'SSW', 'ENE', 'ESE', 'WNW', 'WSW'/)
! -------------------------------------------------------------------------------------------------------------------------
!> \details
!! = log ======================================================================================
!! \n
!! 08/11/16 - switch to v0.4 \n
!! 06/01/17 - use RMA to synchronize data \n
!! 31/01/17 - switch to 3D, v0.5 \n
!! 12/04/17 - redundant ghost nodes workaround
!
! ********************************************************************************************

subroutine synchronize_ghosts(  params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)

    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    ! loop variables
    integer(kind=ik)                    :: k, N, i, j!, dF

    ! grid parameter
    integer(kind=ik)                    :: g, Bs

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank
    ! number of processes
    integer(kind=ik)                    :: number_procs

    ! communication lists:
    ! dim 1: list elements
    ! dim 2: columns
    !                       1   rank of sender process
    !                       2   rank of receiver process
    !                       3   sender block heavy data id
    !                       4   receiver block heavy data id
    !                       5   sender block neighborhood to receiver (dirs id)
    !                       6   difference between sender-receiver level
    ! dim 3: receiver proc rank
    integer(kind=ik), allocatable       :: com_lists(:, :, :), com_lists2(:, :, :)

    ! allocation error variable
    integer(kind=ik)                    :: allocate_error

    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0, sub_t1

    ! communications matrix:
    ! count the number of communications between procs
    ! row/column number encodes process rank + 1
    ! com matrix pos: position in send buffer
    integer(kind=ik), allocatable       :: com_matrix(:,:), com_matrix_pos(:,:), my_com_matrix(:,:), &
                                           com_matrix2(:,:), com_matrix_pos2(:,:), my_com_matrix2(:,:)

    ! variable for non-uniform mesh correction: remove redundant node between fine->coarse blocks
    integer(kind=ik)                                :: rmv_redundant

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    N = params%number_blocks

    ! grid parameter
    Bs = params%number_block_nodes
    g  = params%number_ghost_nodes

    rmv_redundant = 0

    ! set non-uniform mesh correction
    if ( params%non_uniform_mesh_correction ) then
        rmv_redundant = 1
    else
        rmv_redundant = 0
    end if

    ! set MPI parameter
    rank         = params%rank
    number_procs = params%number_procs

    ! allocate local com_lists
    allocate( com_lists( N*16, 6, number_procs), stat=allocate_error )
    !call check_allocation(allocate_error)
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if

    allocate( com_lists2( N*16, 6, number_procs), stat=allocate_error )
    !call check_allocation(allocate_error)
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if

    ! allocate com matrix
    allocate( com_matrix(number_procs, number_procs), stat=allocate_error )
    !call check_allocation(allocate_error)
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if
    allocate( com_matrix_pos(number_procs, number_procs), stat=allocate_error )
    !call check_allocation(allocate_error)
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if
    allocate( my_com_matrix(number_procs, number_procs), stat=allocate_error )
    !call check_allocation(allocate_error)
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if

    allocate( com_matrix2(number_procs, number_procs), stat=allocate_error )
    !call check_allocation(allocate_error)
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if
    allocate( com_matrix_pos2(number_procs, number_procs), stat=allocate_error )
    !call check_allocation(allocate_error)
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if
    allocate( my_com_matrix2(number_procs, number_procs), stat=allocate_error )
    !call check_allocation(allocate_error)
    if ( allocate_error /= 0 ) then
        write(*,'(80("_"))')
        write(*,*) "ERROR: memory allocation fails"
        stop
    end if

    ! reset com-list, com_plan, com matrix, receiver lists
    com_lists       = -1
    com_lists2      = -1

    com_matrix      =  0
    com_matrix_pos  =  0
    my_com_matrix   =  0

    com_matrix2     =  0
    com_matrix_pos2 =  0
    my_com_matrix2  =  0

!    ! reset ghost nodes for all active blocks
!    ! loop over all active blocks
!    do k = 1, hvy_n
!        ! reset ghost nodes
!        if ( params%threeD_case ) then
!            ! 3D:
!!            hvy_block(1:g, :, :, dF, hvy_active(k) )           = 99.0_rk!9.0e9_rk
!!            hvy_block(Bs+g+1:Bs+2*g, :, :, dF, hvy_active(k) ) = 99.0_rk!9.0e9_rk
!!            hvy_block(:, 1:g, :, dF, hvy_active(k) )           = 99.0_rk!9.0e9_rk
!!            hvy_block(:, Bs+g+1:Bs+2*g, :, dF, hvy_active(k) ) = 99.0_rk!9.0e9_rk
!!            hvy_block(:, :, 1:g, dF, hvy_active(k) )           = 99.0_rk!9.0e9_rk
!!            hvy_block(:, :, Bs+g+1:Bs+2*g, dF, hvy_active(k) ) = 99.0_rk!9.0e9_rk
!        else
!            ! 2D:
!            hvy_block(1:g, :, 1, 2, hvy_active(k) )           = 9.0e9_rk
!            hvy_block(Bs+g+1:Bs+2*g, :, 1, 2, hvy_active(k) ) = 9.0e9_rk
!            hvy_block(:, 1:g, 1, 2, hvy_active(k) )           = 9.0e9_rk
!            hvy_block(:, Bs+g+1:Bs+2*g, 1, 2, hvy_active(k) ) = 9.0e9_rk
!        end if
!    end do

!---------------------------------------------------------------------------------------------
! main body

    ! start time
    sub_t0 = MPI_Wtime()

    ! ----------------------------------------------------------------------------------------
    ! first: synchronize internal ghost nodes, create com_list for external communications

    ! copy internal nodes and create com_matrix/com_lists for external communications
    call synchronize_internal_nodes( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n, my_com_matrix, com_lists, 1 )
    call synchronize_internal_nodes( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n, my_com_matrix2, com_lists2, 2 )

    ! end time
    sub_t1 = MPI_Wtime()
    ! write time
    if ( params%debug ) then
        ! find free or corresponding line
        k = 1
        do while ( debug%name_comp_time(k) /= "---" )
            ! entry for current subroutine exists
            if ( debug%name_comp_time(k) == "synch. ghosts - internal" ) exit
            k = k + 1
        end do
        ! write time
        debug%name_comp_time(k) = "synch. ghosts - internal"
        debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1
        debug%comp_time(k, 2)   = debug%comp_time(k, 2) + sub_t1 - sub_t0
    end if

    ! start time
    sub_t0 = MPI_Wtime()

    ! synchronize com matrix
    call MPI_Allreduce(my_com_matrix, com_matrix, number_procs*number_procs, MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Allreduce(my_com_matrix2, com_matrix2, number_procs*number_procs, MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! end time
    sub_t1 = MPI_Wtime()
    ! write time
    if ( params%debug ) then
        ! find free or corresponding line
        k = 1
        do while ( debug%name_comp_time(k) /= "---" )
            ! entry for current subroutine exists
            if ( debug%name_comp_time(k) == "synch. ghosts - com_matrix" ) exit
            k = k + 1
        end do
        ! write time
        debug%name_comp_time(k) = "synch. ghosts - com_matrix"
        debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1
        debug%comp_time(k, 2)   = debug%comp_time(k, 2) + sub_t1 - sub_t0
    end if

    ! symmetrize com matrix - choose max com number
    do i = 1, number_procs
        do j = i+1, number_procs
            com_matrix(i,j) = max( com_matrix(i,j), com_matrix(j,i) )
            com_matrix(j,i) = com_matrix(i,j)
            com_matrix2(i,j) = max( com_matrix2(i,j), com_matrix2(j,i) )
            com_matrix2(j,i) = com_matrix2(i,j)
        end do
    end do

    ! save com matrix
    if ( params%debug ) then
        call write_com_matrix( com_matrix )
    end if
    if ( params%debug ) then
        call write_com_matrix( com_matrix2 )
    end if

    ! start time
    sub_t0 = MPI_Wtime()

    ! call external nodes synchronization
    ! case 1
    call synchronize_external_nodes(  params, hvy_block, com_lists, com_matrix )
    ! case 2
    call synchronize_external_nodes(  params, hvy_block, com_lists2, com_matrix2 )

    ! workaround: second internal synchronization to overwrite external redundant nodes
    !call synchronize_internal_nodes( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n, my_com_matrix, com_lists )
    call synchronize_internal_nodes( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n, my_com_matrix, com_lists, 1 )
    call synchronize_internal_nodes( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n, my_com_matrix2, com_lists2, 2 )

    ! clean up
    deallocate( com_lists, stat=allocate_error )
    deallocate( com_lists2, stat=allocate_error )

    deallocate( com_matrix, stat=allocate_error )
    deallocate( com_matrix2, stat=allocate_error )
    deallocate( com_matrix_pos, stat=allocate_error )
    deallocate( com_matrix_pos2, stat=allocate_error )
    deallocate( my_com_matrix, stat=allocate_error )
    deallocate( my_com_matrix2, stat=allocate_error )

    ! end time
    sub_t1 = MPI_Wtime()
    ! write time
    if ( params%debug ) then
        ! find free or corresponding line
        k = 1
        do while ( debug%name_comp_time(k) /= "---" )
            ! entry for current subroutine exists
            if ( debug%name_comp_time(k) == "synch. ghosts - write external" ) exit
            k = k + 1
        end do
        ! write time
        debug%name_comp_time(k) = "synch. ghosts - write external"
        debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1
        debug%comp_time(k, 2)   = debug%comp_time(k, 2) + sub_t1 - sub_t0
    end if

end subroutine synchronize_ghosts
