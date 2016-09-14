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

    real(kind=rk), dimension(:,:), allocatable      :: u1, u2, u3, u5, u6
    real(kind=rk)                                   :: detail1, detail2

    integer(kind=ik)                                :: dF, k, N, block_num, Bs, g, allocate_error, max_status
    integer(kind=ik), dimension(:,:), allocatable   :: dF_status

    N  = size(blocks_params%active_list, dim=1)
    Bs = blocks_params%size_block
    g  = blocks_params%number_ghost_nodes

    allocate( u1(1:Bs+2*g,1:Bs+2*g), stat=allocate_error )
    allocate( u2(1:Bs+2*g,1:Bs+2*g), stat=allocate_error )
    ! coarsen field are half block size + 1/2
    allocate( u3( 1:(Bs+1)/2 + g , 1:(Bs+1)/2 + g), stat=allocate_error )
    ! helper fields to create coarsen field, two level down
    allocate( u5( 1:(Bs+1)/2 + g, 1:(Bs+1)/2 + g ), stat=allocate_error )
    allocate( u6( 1:((Bs+1)/2 + g + 1)/2, 1:((Bs+1)/2 + g + 1)/2 ), stat=allocate_error )

    ! allocate memory for array with refinement status for all fields
    allocate( dF_status(N, blocks_params%number_data_fields), stat=allocate_error )

    dF_status = 0

    ! synchronize ghostnodes
    call synchronize_ghosts()

    ! clear old refinement status for all blocks
    do k = 1, N
        block_num = blocks_params%active_list(k)
        blocks(block_num)%refinement = 0
        blocks(block_num)%data_fields(:)%detail = 0.0_rk
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
            ! calculate wavelet coefficient (detail) for the next two lower grid levels
            ! one level coarser
            call calculate_detail(detail1, u1, u2, Bs+2*g)

            ! two level coarser
            u5 = 0.0_rk
            u6 = 0.0_rk

            call restriction_2D(u3, u6)  ! fine, coarse
            call prediction_2D (u6, u5)  ! coarse, fine
            call calculate_detail(detail2, u3, u5, (Bs+1)/2 + g)

            ! error
            ! case 1: detail1 = coarsen, detail2 = coarsen => coarse block, status -1
            if ( (detail1 < params%eps_coarsen) .and. (detail2 < params%eps_coarsen) ) then
                dF_status(k, dF) = -1
                blocks(block_num)%data_fields(dF)%detail = max(detail1, detail2)
            ! case 2: detail1 = refine, detail2 does not matter => refine block, status +1
            elseif (detail1 > params%eps_coarsen) then
                dF_status(k, dF) = 1
                blocks(block_num)%data_fields(dF)%detail = detail1
            ! case 3: detail1 = coarsen, detail2 = refine => stay, nothing to do
            else
                blocks(block_num)%data_fields(dF)%detail = max(detail1, detail2)
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
    deallocate( u5, stat=allocate_error )
    deallocate( u6, stat=allocate_error )
    deallocate( dF_status, stat=allocate_error )

end subroutine local_refinement_status
