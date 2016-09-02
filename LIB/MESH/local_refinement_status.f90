! ********************************
! 2D AMR prototype
! --------------------------------
!
! set refinement status for all blocks
!
! name: local_refinement_status.f90
! date: 17.08.2016
! author: msr
! version: 0.1
!
! ********************************

subroutine local_refinement_status()

    use module_params
    use module_blocks
    use module_interpolation

    implicit none

    real(kind=rk), dimension(:,:), allocatable      :: u1, u2, u3
    real(kind=rk)                                   :: error

    integer(kind=ik)                                :: dF, k, N, block_num, Bs, g, allocate_error, i, j, max_status
    integer(kind=ik), dimension(:,:), allocatable   :: dF_status

    N  = size(blocks_params%active_list, dim=1)
    Bs = blocks_params%size_block
    g  = blocks_params%number_ghost_nodes

    allocate( u1(1:Bs+2*g,1:Bs+2*g), stat=allocate_error )
    allocate( u2(1:Bs+2*g,1:Bs+2*g), stat=allocate_error )
    ! coarsen field are half block size + 1/2
    allocate( u3( 1:(Bs+1)/2 + g , 1:(Bs+1)/2 + g), stat=allocate_error )

    ! allocate memory for array with refinement status for all fields
    allocate( dF_status(N, blocks_params%number_data_fields), stat=allocate_error )

    dF_status = 0

    ! synchronize ghostnodes
    call synchronize_ghosts()

    ! clear old refinement status for all blocks
    do k = 1, N
        block_num = blocks_params%active_list(k)
        blocks(block_num)%refinement = 0
    end do

    ! loop over all fields
    do dF = 1, blocks_params%number_data_fields

        ! loop over all active blocks
        do k = 1, N

            block_num = blocks_params%active_list(k)
            u1        = blocks(block_num)%data_fields(dF)%data_(:,:)
            u2        = 0.0_rk
            u3        = 0.0_rk

            ! wavelet indicator
            call restriction_2D(u1, u3)  ! fine, coarse
            call prediction_2D (u3, u2)  ! coarse, fine

            ! error
            error = 0.0_rk
            do i = 1, size(u1,1)
                do j = 1, size(u1,1)
                    error = max( error, sqrt( (u1(i,j)-u2(i,j)) * ( u1(i,j)-u2(i,j)) ) )
                end do
            end do

            if (error < params%eps_coarsen) then
                ! coarsen block, -1
                dF_status(k, dF) = -1
            elseif (error > params%eps_refine) then
                ! refine block, +1
                dF_status(k, dF) = 1
            end if

        end do

    end do

    ! set block refinement status
    do k = 1, N

        block_num = blocks_params%active_list(k)

        max_status = -99
        ! loop over all data fields
        do dF = 1, blocks_params%number_data_fields
            ! block refinement status is the maximal status
            max_status = max( max_status, dF_status(k, dF) )
        end do

        blocks(block_num)%refinement = max_status

    end do

    ! check if block has reached maximal level
    call respect_min_max_treelevel()

    deallocate( u1, stat=allocate_error )
    deallocate( u2, stat=allocate_error )
    deallocate( u3, stat=allocate_error )
    deallocate( dF_status, stat=allocate_error )

end subroutine local_refinement_status
