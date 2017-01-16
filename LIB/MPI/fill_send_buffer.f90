! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: fill_send_buffer.f90
! version: 0.4
! author: msr
!
! fill send buffer
!
! input:    - params, heavy data
!           - com matrix line
!           - proc rank
! output:   - filled send buffer
!           - second com matrix line with receiver data position (column number) in send buffer
!
! = log ======================================================================================
!
! 13/01/17 - create for v0.4
! ********************************************************************************************

subroutine fill_send_buffer( params, hvy_block, com_lists, com_matrix_line, rank, int_send_buffer, real_send_buffer )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(in)                  :: params

    ! heavy data array - block data
    real(kind=rk), intent(in)                       :: hvy_block(:, :, :, :)

    ! communication lists:
    integer(kind=ik), intent(in)                    :: com_lists(:, :, :)

    ! com matrix line
    integer(kind=ik), intent(in)                    :: com_matrix_line(:)

    ! proc rank
    integer(kind=ik), intent(in)                    :: rank

    ! integer send buffer
    integer(kind=ik), intent(inout)                 :: int_send_buffer(:,:)
    ! real send buffer
    real(kind=rk), intent(inout)                    :: real_send_buffer(:,:)

    ! loop variable
    integer(kind=ik)                                :: k, i

    ! column number of send buffer, position in integer buffer
    integer(kind=ik)                                :: column_pos, int_pos

    ! send buffer for one proc
    real(kind=rk), allocatable                      :: proc_send_buffer(:)

    ! allocation error variable
    integer(kind=ik)                                :: allocate_error

    ! index of send buffer, return from create_send_buffer subroutine
    integer(kind=ik)                                :: buffer_i

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! reset column number
    column_pos = 1

    ! allocate proc send buffer, size = line size of real send buffer
    allocate( proc_send_buffer( size(real_send_buffer,1) ), stat=allocate_error )

!---------------------------------------------------------------------------------------------
! main body

    ! loop over all line elements
    do k = 1, size(com_matrix_line,1)

        ! communication to other proc, do not work with internal communications
        if ( (com_matrix_line(k) /= 0) .and. (k /= rank+1) ) then

            ! first: real data
            ! ----------------

            ! write real send buffer for proc k
            call create_send_buffer(params, hvy_block, com_lists( 1:com_matrix_line(k), :, k), com_matrix_line(k), proc_send_buffer, buffer_i)

            ! real buffer entry
            real_send_buffer( 1 : buffer_i, column_pos ) = proc_send_buffer( 1 : buffer_i )

            ! second: integer data
            ! --------------------

            ! save real buffer length
            int_send_buffer(1, column_pos) = buffer_i

            ! reset position
            int_pos = 2

            ! loop over all communications to this proc
            do i = 1, com_matrix_line(k)

                ! int buffer entry: neighbor block id, neighborhood, level difference
                int_send_buffer( int_pos  , column_pos ) = com_lists( i, 4, k)
                int_send_buffer( int_pos+1, column_pos ) = com_lists( i, 5, k)
                int_send_buffer( int_pos+2, column_pos ) = com_lists( i, 6, k)
                ! increase int buffer position
                int_pos = int_pos + 3

            end do

            ! third: increase column number
            ! ------------------------------------
            column_pos = column_pos + 1

        end if

    end do

    ! clean up
    deallocate( proc_send_buffer, stat=allocate_error )

end subroutine fill_send_buffer
