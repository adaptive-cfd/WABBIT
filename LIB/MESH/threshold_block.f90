! ********************************
! WABBIT
! --------------------------------
!
! thresholding for all blocks
!
! name: threshold_block.f90
! date: 30.09.2016
! author: msr
! version: 0.2
!
! ********************************

subroutine threshold_block()

    use module_params
    use module_blocks
    use module_interpolation

    implicit none

    real(kind=rk), dimension(:,:), allocatable      :: u1, u2, u3
    real(kind=rk)                                   :: detail

    integer(kind=ik)                                :: dF, k, N, Bs, g, allocate_error, max_status
    integer(kind=ik), dimension(:), allocatable     :: dF_status

    N  = blocks_params%number_max_blocks
    Bs = blocks_params%size_block
    g  = blocks_params%number_ghost_nodes

    allocate( u1(1:Bs+2*g,1:Bs+2*g), stat=allocate_error )
    allocate( u2(1:Bs+2*g,1:Bs+2*g), stat=allocate_error )
    ! coarsen field are half block size + 1/2
    allocate( u3( 1:(Bs+1)/2 + g , 1:(Bs+1)/2 + g), stat=allocate_error )

    ! allocate memory for array with refinement status for all fields
    allocate( dF_status(blocks_params%number_data_fields), stat=allocate_error )

    dF_status = 0

    ! synchronize ghostnodes
    call synchronize_ghosts()

    ! clear old refinement status for all blocks
    do k = 1, N
        blocks(k)%refinement            = 0
        blocks(k)%data_fields(:)%detail = 0.0_rk
    end do

    ! loop over all blocks
    do k = 1, N

        ! loop over all fields
        do dF = 1, blocks_params%number_data_fields

            if (blocks(k)%active) then

                dF_status = 0

                u1        = blocks(k)%data_fields(dF)%data_(:,:)
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
                blocks(k)%data_fields(dF)%detail = detail

            end if

        end do

        max_status = -99
        ! loop over all data fields
        do dF = 1, blocks_params%number_data_fields
            ! block refinement status is the maximal status
            max_status = max( max_status, dF_status(dF) )
        end do

        blocks(k)%refinement = max_status

    end do

    ! check if block has reached minimal level
    call respect_min_max_treelevel()

    deallocate( u1, stat=allocate_error )
    deallocate( u2, stat=allocate_error )
    deallocate( u3, stat=allocate_error )
    deallocate( dF_status, stat=allocate_error )

end subroutine threshold_block
