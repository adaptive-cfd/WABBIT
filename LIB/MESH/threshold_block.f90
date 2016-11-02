! ********************************
! WABBIT
! --------------------------------
!
! thresholding for all blocks
!
! name: threshold_block.f90
! date: 28.10.2016
! author: msr
! version: 0.3
!
! ********************************

subroutine threshold_block()

    use mpi
    use module_params
    use module_blocks
    use module_interpolation

    implicit none

    real(kind=rk), dimension(:,:), allocatable      :: u1, u2, u3
    real(kind=rk)                                   :: detail

    integer(kind=ik)                                :: dF, k, l, buff_i, N_light, N_heavy, Bs, g, allocate_error, max_status
    integer(kind=ik), dimension(:), allocatable     :: dF_status
    integer(kind=ik)                                :: rank, ierr, n_proc
    integer(kind=ik) , dimension(10000)              :: send_buff, recv_buff

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, n_proc, ierr)

    N_light = blocks_params%number_max_blocks
    N_heavy = blocks_params%number_max_blocks_data
    Bs      = blocks_params%size_block
    g       = blocks_params%number_ghost_nodes

    allocate( u1(1:Bs+2*g,1:Bs+2*g), stat=allocate_error )
    allocate( u2(1:Bs+2*g,1:Bs+2*g), stat=allocate_error )
    ! coarsen field are half block size + 1/2
    allocate( u3( 1:(Bs+1)/2 + g , 1:(Bs+1)/2 + g), stat=allocate_error )

    ! allocate memory for array with refinement status for all fields
    allocate( dF_status(blocks_params%number_data_fields), stat=allocate_error )

    dF_status = 0

    ! synchronize ghostnodes
    call synchronize_ghosts()

    ! clear old refinement status for all blocks, work on light data
    do k = 1, N_light
        blocks(k)%refinement = 0
    end do

    ! loop over all blocks to calculate the detail, work on heavy data
    do k = 1, N_heavy

        ! loop over all fields
        do dF = 1, blocks_params%number_data_fields

            if ( blocks_data(k)%block_id /= -1 ) then
                ! block is active

                ! reset old detail
                blocks_data(k)%data_fields(:)%detail = 0.0_rk

                ! reset dF status
                dF_status = 0

                ! reset interpolation fields
                u1        = blocks_data(k)%data_fields(dF)%data_(:,:)
                u2        = 0.0_rk
                u3        = 0.0_rk

                ! wavelet indicator
                call restriction_2D(u1, u3)  ! fine, coarse
                call prediction_2D (u3, u2)  ! coarse, fine
                call calculate_detail(detail, u1, u2, Bs+2*g)

                ! threshold
                if (detail < params%eps) then
                    ! coarsen block, -1
                    dF_status(dF) = -1
                end if

                ! save detail
                blocks_data(k)%data_fields(dF)%detail = detail

            end if

        end do

        ! block is active
        if ( blocks_data(k)%block_id /= -1 ) then

            max_status = -99
            ! loop over all data fields
            do dF = 1, blocks_params%number_data_fields
                ! block refinement status is the maximal status
                max_status = max( max_status, dF_status(dF) )
            end do

            ! save new refinement status in light data, note: light data need to synchronize later!
            blocks( blocks_data(k)%block_id )%refinement = max_status

        end if

    end do

    ! synchronize light data (refinement status)
    ! note: send only light-data corresponding to proc
    ! loop over all processes
    do k = 1, n_proc

        if ( (k-1) == rank ) then
            ! send data
            ! create send buffer
            send_buff = -99
            buff_i = 0
            do l = 1, N_light
                if ( blocks(l)%proc_rank == (k-1) ) then
                    buff_i              = buff_i + 1
                    send_buff(buff_i)   = blocks(l)%refinement
                end if
            end do
            ! send data
            call MPI_Bcast(send_buff, 10000, MPI_INTEGER4, k-1, MPI_COMM_WORLD, ierr)

        else

            ! receive data
            recv_buff = -99
            call MPI_Bcast(recv_buff, 10000, MPI_INTEGER4, k-1, MPI_COMM_WORLD, ierr)
            ! synchronize light data
            buff_i = 0
            do l = 1, N_light
                if ( blocks(l)%proc_rank == (k-1) ) then
                    buff_i               = buff_i + 1
                    blocks(l)%refinement = recv_buff(buff_i)
                end if
            end do

        end if

    end do

    ! check if block has reached minimal level
    call respect_min_max_treelevel()

    deallocate( u1, stat=allocate_error )
    deallocate( u2, stat=allocate_error )
    deallocate( u3, stat=allocate_error )
    deallocate( dF_status, stat=allocate_error )

end subroutine threshold_block
