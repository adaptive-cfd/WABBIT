!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name write_future_mesh_lvl.f90
!> \version 0.4
!> \author msr
!
!> \brief debug planned mesh changes \n
!! write file with future mesh level for all blocks \n
!! write also future mesh level for all known neighbor blocks \n
!! 
!!
!! input:    - params, light data \n
!! output:   - status of lgt_block synchronzation \n
!!
!!
!! = log ======================================================================================
!! \n
!! 29/11/16 - create
! ********************************************************************************************

subroutine write_future_mesh_lvl( params, lgt_block, hvy_neighbor, lgt_active, lgt_n )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)
    !> neighbor list
    integer(kind=ik), intent(in)        :: hvy_neighbor(:, :)

    !> active_block_list (light data)
    integer(kind=ik), intent(in)        :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_n

    ! local data array, rows: block id, number of columns: future lvl, 16 possible neighbors
    integer(kind=ik)                    :: my_fut_mesh_lvl( size(lgt_block,1) , 17 ), my_fut_mesh_lvl_0( size(lgt_block,1) , 17 )

    ! MPI error variable
    integer(kind=ik)                    :: ierr
    ! process rank
    integer(kind=ik)                    :: rank

    ! file IO error variable
    integer(kind=ik)                    :: io_error

    ! loop variables
    integer(kind=ik)                    :: k, l, lgt_start, lgt_id, proc_id, hvy_id, ref_status

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! determinate process rank
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

    ! light data start id
    call proc_to_lgt_data_start_id( lgt_start, rank, params%number_blocks )

    ! reset local light data
    my_fut_mesh_lvl   = 0
    my_fut_mesh_lvl_0 = 0

!---------------------------------------------------------------------------------------------
! main body

    ! loop over all active blocks
    do k = 1, lgt_n

        ! light data id
        lgt_id = lgt_active(k)
        ! corresponding proc
        call lgt_id_to_proc_rank( proc_id, lgt_id, params%number_blocks )

        if ( proc_id == rank ) then

            ! set refinement status
            ref_status = lgt_block( lgt_id, params%max_treelevel + 2 )
            if ( ref_status == -2) then
                ref_status = -1
            elseif ( ref_status == -1) then
                ref_status = 0
            end if
            ! write future block lvl
            my_fut_mesh_lvl( lgt_id, 1 ) = lgt_block( lgt_id, params%max_treelevel + 1 ) + ref_status

            ! heavy id
            call lgt_id_to_hvy_id( hvy_id, lgt_id, rank, params%number_blocks )

            ! write future neighbor level
            do l = 1, 16
                ! neighbor is active
                if ( hvy_neighbor(hvy_id, l) /= -1 ) then
                    ! set refinement status
                    ref_status = lgt_block( hvy_neighbor(hvy_id, l), params%max_treelevel + 2 )
                    if ( ref_status == -2) then
                        ref_status = -1
                    elseif ( ref_status == -1) then
                        ref_status = 0
                    end if
                    ! write future level
                    my_fut_mesh_lvl( lgt_id, l+1 ) = lgt_block( hvy_neighbor(hvy_id, l), params%max_treelevel + 1 ) + ref_status
                end if
            end do

        end if

    end do

    ! gather data
    call MPI_Allreduce(my_fut_mesh_lvl, my_fut_mesh_lvl_0, size(my_fut_mesh_lvl,1)*size(my_fut_mesh_lvl,2), MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! write file
    if (rank == 0) then

        ! open file
        open(unit=99,file="future_block_lvl.dat",status='replace',action='write', iostat=io_error)

        ! write file header
        write(99, '("block_id",1x,"treecode",1x,"lvl",1x,"neighbors lvl")', advance='no')
        write(99,*)

        ! write data
        do k = 1, lgt_n

            ! light data id
            lgt_id = lgt_active(k)

            ! write light data id
            write(99, '(i5,6x)', advance='no') lgt_id

            ! write treecode
            do l = 1, params%max_treelevel
                if ( lgt_block(lgt_id, l) /= -1 ) then
                    write(99, '(i1)', advance='no') lgt_block(lgt_id, l)
                else
                    write(99, '(1x)', advance='no')
                end if
            end do

            ! write future mesh level
            write(99, '(1x,i2,3x)', advance='no') my_fut_mesh_lvl_0(lgt_id, 1)

            ! write future neighbor level
            do l = 1, 16
                ! neighbor is active
                if ( my_fut_mesh_lvl_0(lgt_id, l+1) /= 0 ) then
                    write(99, '(i2,1x)', advance='no') my_fut_mesh_lvl_0(lgt_id, l+1)
                end if
            end do

            ! next line
            write(99,*)

        end do

        ! close file
        close(unit=99)

    end if

end subroutine write_future_mesh_lvl
