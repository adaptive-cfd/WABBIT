subroutine check_redundant_nodes_clean( params, lgt_block, hvy_block, hvy_neighbor, &
    hvy_active, hvy_n, stop_status )

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
    ! status of nodes check: if true: stops program
    logical, intent(inout)              :: stop_status

    ! MPI parameter
    integer(kind=ik)                    :: myrank
    ! loop variables
    integer(kind=ik)                    :: N, k, l, neighborhood, level_diff
    ! id integers
    integer(kind=ik)                    :: lgt_id, neighbor_lgt_id, neighbor_rank, hvy_id
    integer(kind=ik), dimension(2,3)    :: ijk, ijk2
    ! data buffer size
    integer(kind=ik)                    :: buffer_size, buffer_position
    ! grid parameter
    integer(kind=ik)                    :: g, i0, ii0, ii1
    integer(kind=ik), dimension(3)      :: Bs
    ! number of datafields
    integer(kind=ik)                    :: NdF
    logical                             :: test2
    integer(kind=ik), parameter         :: stage=1

 !---------------------------------------------------------------------------------------------
! variables initialization

    if (.not. ghost_nodes_module_ready) then
        call init_ghost_nodes( params )
    endif

    ! reset status
    stop_status = .false.

    Bs       = params%Bs
    g        = params%n_ghosts
    NdF      = params%n_eqn
    N        = params%number_blocks
    myrank   = params%rank


    ! the (module-global) communication_counter is the number of neighboring relations
    ! this rank has with all other ranks (it is thus an array of number_procs)
    communication_counter(:, stage) = 0_ik


    int_pos(:, stage) = 1
    real_pos(:, stage) = 0

    ! in this routine we treat internal and external neighbors identically, therefore we need to
    ! count also my own neighorhood relations. this is not highly efficient (involves copying)
    ! but this routine is meant for testing, not production.
    call get_my_sendrecv_amount_with_ranks(params, lgt_block, hvy_neighbor, hvy_active, hvy_n, &
    recv_counter(:, stage), send_counter(:, stage), &
    int_recv_counter(:, stage), int_send_counter(:, stage), INCLUDE_REDUNDANT, .true.)


    ! reset int_send_buffer, but only the parts that will actually be treated.
    do k = 1, params%number_procs
        ii0 = sum(int_recv_counter(0:(k-1)-1, stage)) + 1
        ii1 = ii0 + int_recv_counter(k-1, stage)
        int_send_buffer(ii0:ii1, stage) = -99
    enddo


    do k = 1, hvy_n
        do neighborhood = 1, size(hvy_neighbor, 2)
            ! neighbor exists
            if ( hvy_neighbor( hvy_active(k), neighborhood ) /= -1 ) then

                neighbor_lgt_id = hvy_neighbor( hvy_active(k), neighborhood )
                call lgt_id_to_proc_rank( neighbor_rank, neighbor_lgt_id, N )
                call hvy_id_to_lgt_id( lgt_id, hvy_active(k), myrank, N )
                call lgt_id_to_hvy_id( hvy_id, neighbor_lgt_id, neighbor_rank, N )

                ! define leveldiff: sender - receiver, so +1 means sender on higher level. sender is active block (me)
                level_diff = lgt_block( lgt_id, params%max_treelevel + idx_mesh_lvl ) - lgt_block( neighbor_lgt_id, params%max_treelevel + idx_mesh_lvl )

                ijk = ijkGhosts(:,:, neighborhood, level_diff, INCLUDE_REDUNDANT, SENDER)

                if ( level_diff == 0 ) then
                    !-----------------------------------------------------------
                    ! same level
                    !-----------------------------------------------------------
                    call GhostLayer2Line( params, line_buffer, buffer_size, &
                    hvy_block( ijk(1,1):ijk(2,1), ijk(1,2):ijk(2,2), ijk(1,3):ijk(2,3), :, hvy_active(k)) )

                else
                    !-----------------------------------------------------------
                    ! different level
                    !-----------------------------------------------------------
                    call restrict_predict_data( params, res_pre_data, ijk, neighborhood, level_diff, hvy_block, hvy_active(k))

                    ijk2 = ijkGhosts(1:2, 1:3, neighborhood, level_diff, INCLUDE_REDUNDANT, RESPRE)

                    call GhostLayer2Line( params, line_buffer, buffer_size, res_pre_data( ijk2(1,1):ijk2(2,1), &
                    ijk2(1,2):ijk2(2,2), ijk2(1,3):ijk2(2,3),:) )

                end if


                ! first: fill com matrix, count number of communication to neighboring process, needed for int buffer length
                communication_counter(neighbor_rank+1, stage) = communication_counter(neighbor_rank+1, stage) + 1
                ! active block send data to its neighbor block
                ! fill int/real buffer
                call AppendLineToBuffer( int_send_buffer, new_send_buffer, buffer_size, neighbor_rank+1, line_buffer, &
                hvy_id, neighborhood, level_diff, stage )

            end if
        end do
    end do

    !***********************************************************************
    ! transfer part (send/recv)
    !***********************************************************************
    ! send/receive data
    call isend_irecv_data_2( params, int_send_buffer, new_send_buffer, int_recv_buffer, &
    new_recv_buffer, communication_counter, stage )

    !***********************************************************************
    ! Unpack received data and compare with ghost nodes data
    !***********************************************************************
    ! loop over all mpiranks
    do k = 1, params%number_procs
        if ( communication_counter(k, stage) /= 0 ) then

            ! start index of this mpirank in the int_buffer
            l = sum(int_recv_counter(0:k-1-1, stage)) +1

            ! -99 marks end of data
            do while ( int_recv_buffer(l, stage) > -99 )

                hvy_id          = int_recv_buffer(l,   stage)
                neighborhood    = int_recv_buffer(l+1, stage)
                level_diff      = int_recv_buffer(l+2, stage)
                buffer_position = int_recv_buffer(l+3, stage)
                buffer_size     = int_recv_buffer(l+4, stage)

                !line_buffer(1:buffer_size) = real_receive_buffer( buffer_position : buffer_position-1 + buffer_size, k, 1 )
                i0 = sum(recv_counter(0:k-1-1, stage)) + buffer_position
                line_buffer(1:buffer_size) = new_recv_buffer( i0 : i0+buffer_size-1, stage )

                ! data bounds (2-recv)
                ijk = ijkGhosts(:,:, neighborhood, level_diff, INCLUDE_REDUNDANT, RECVER)

                ! compare data
                call hvy_id_to_lgt_id( lgt_id, hvy_id, myrank, N )
                call compare_hvy_data( params, line_buffer, ijk, hvy_block, hvy_id, stop_status, level_diff, &
                lgt_block(lgt_id, params%max_treelevel + idx_refine_sts), treecode2int( lgt_block(lgt_id, 1:params%max_treelevel) ) )

                l = l + 5
            end do
        end if
    end do

    ! MPI sync the stop status
    test2 = stop_status
    call MPI_Allreduce(test2, stop_status, 1, MPI_LOGICAL, MPI_LOR, WABBIT_COMM, k )
end subroutine check_redundant_nodes_clean




subroutine compare_hvy_data( params, line_buffer, ijk, hvy_block, hvy_id, stop_status, level_diff, my_ref, tc )

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)                  :: params
    !> data buffer
    real(kind=rk), intent(inout)                    :: line_buffer(:)
    !> ijk
    integer(kind=ik), intent(inout)                 :: ijk(2,3)
    !> heavy data array - block data
    real(kind=rk), intent(inout)                    :: hvy_block(:, :, :, :, :)
    !> hvy id
    integer(kind=ik), intent(in)                    :: hvy_id, level_diff, my_ref
    ! status of nodes check: if true: stops program
    logical, intent(inout)              :: stop_status
    integer(kind=tsize)::tc

    ! loop variable
    integer(kind=ik)                                :: i, j, k, dF, buffer_i, oddeven, g
    integer(kind=ik), dimension(3)                  :: Bs

    ! error threshold
    real(kind=rk)                                   :: eps


    ! error norm
    real(kind=rk)       :: error_norm

    Bs = params%Bs
    g = params%n_ghosts

!---------------------------------------------------------------------------------------------
! variables initialization
    buffer_i = 1

    ! NOTE: newer versions do not compare actual data, but fill the hvy_blocks with
    ! their lgt_ids. This makes the comparison much easier, as those values differ by at least 1.0_rk
    eps = 1e-6_rk

    ! reset error norm
    error_norm = 0.0_rk

    ! the first index of the redundant points is (g+1, g+1, g+1)
    ! so if g is even, then we must compare the odd indices i,j,k on the lines
    ! of the redundant points.
    ! if g is odd, then we must compare the even ones
    ! Further note that BS is odd (always), so as odd+even=odd and odd+odd=even
    ! we can simply study the parity of g
    oddeven = mod(params%n_ghosts,2)

!---------------------------------------------------------------------------------------------
! main body
    ! loop over all data fields
    do dF = 1, params%n_eqn
        ! third dimension, note: for 2D cases k is always 1
        do k = ijk(1,3), ijk(2,3)
            ! second dimension
            do j = ijk(1,2), ijk(2,2)
                ! first dimension
                do i = ijk(1,1), ijk(2,1)

                    if (level_diff /= -1) then
                        ! on the same or coarser level, the comparison just takes all points, no odd/even downsampling required.
                        error_norm = max(error_norm, abs(hvy_block( i, j, k, dF, hvy_id ) - line_buffer( buffer_i )))

                        hvy_block_test_err( i, j, k, dF, hvy_id ) = abs(hvy_block( i, j, k, dF, hvy_id ) - line_buffer( buffer_i ))
                        hvy_block_test_val( i, j, k, dF, hvy_id ) = hvy_block( i, j, k, dF, hvy_id )
                        hvy_block_test_interpref( i, j, k, dF, hvy_id ) = line_buffer( buffer_i )

                    else
                        ! if the level diff is -1, I compare with interpolated (upsampled) data. that means every EVEN
                        ! point is the result of interpolation, and not truely redundant.
                        ! Note this routine ALWAYS just compares the redundant nodes, so it will mostly be called
                        ! with a line of points (i.e. one dimension is length one)
                        !
                        ! This routine has been tested:
                        !   - old method (working version): no error found (okay)
                        !   - old method, non_uniform_mesh_correction=0; in params file -> plenty of errors (okay)
                        !   - old method, sync stage 4 deactivated: finds all occurances of "3finer blocks on corner problem" (okay)
                        !   - new method, averaging, no error found (makes sense: okay)
                        if (oddeven==0) then
                            ! even number of ghost nodes -> comparison on ODD points
                             if ( (mod(i,2)/=0) .and. (mod(j,2)/=0) .and. (mod(k,2)/=0) ) then
                                error_norm = max(error_norm, abs(hvy_block( i, j, k, dF, hvy_id ) - line_buffer( buffer_i )))

                                hvy_block_test_err( i, j, k, dF, hvy_id ) = abs(hvy_block( i, j, k, dF, hvy_id ) - line_buffer( buffer_i ))
                                hvy_block_test_val( i, j, k, dF, hvy_id ) = hvy_block( i, j, k, dF, hvy_id )
                                hvy_block_test_interpref( i, j, k, dF, hvy_id ) = line_buffer( buffer_i )
                            endif
                        else
                            ! odd number of ghost nodes -> comparison on EVEN points
                            if ( (mod(i,2)==0) .and. (mod(j,2)==0) .and. (mod(k,2)==0) ) then
                                error_norm = max(error_norm, abs(hvy_block( i, j, k, dF, hvy_id ) - line_buffer( buffer_i )))

                                hvy_block_test_err( i, j, k, dF, hvy_id ) = abs(hvy_block( i, j, k, dF, hvy_id ) - line_buffer( buffer_i ))
                                hvy_block_test_val( i, j, k, dF, hvy_id ) = hvy_block( i, j, k, dF, hvy_id )
                                hvy_block_test_interpref( i, j, k, dF, hvy_id ) = line_buffer( buffer_i )
                            endif
                        endif
                    endif
                    buffer_i = buffer_i + 1
                end do
            end do
        end do
    end do

    if (error_norm > eps)  then
        write(*,'("ERROR: difference in redundant nodes ",es12.4," level_diff=",i2, " hvy_id=",i6,1x,i6," rank=",i5 )') &
        error_norm, level_diff, nint(hvy_block( size(hvy_block,1)/2, size(hvy_block,2)/2, 1, 1, hvy_id )), hvy_id, params%rank
        write(*,*) "refinement status", my_ref, "tc=", tc
        ! stop program
        stop_status = .true.
    end if

end subroutine compare_hvy_data

!############################################################################################################

subroutine write_real5(data_block,hvy_active, hvy_n, fileName, rank )
    ! dump all data including ghost nodes for debugging, eg with matlab:

!    function   [data ] =  read(fileName)

!    fid=fopen(fileName, 'rb');        % Open the file.
!    [dataSize, count ] =fread(fid, 5, 'int32') ;
!    data = zeros(dataSize') ;

!    allVals = zeros(prod(dataSize),1) ;
!    allInd  = zeros(prod(dataSize),5)  ;

!    lineNum = 0 ;
!    while(1)
!     [coord, countC ] =fread(fid, 5, 'int32');
!     [val, countV ] =fread(fid, 1, 'float64') ;

!     if (countC*countV == 0 )
!         disp('all')
!         disp ( countC)
!         disp ( countV)
!         if (prod(dataSize) ~= lineNum)
!             fclose(fid) ;
!             error('wrong linenumberm better check it')
!         end
!         %!disp(lineNum)
!         break
!     end
!     lineNum = lineNum + 1 ;
!     allVals(lineNum)  = val ;
!     allInd(lineNum,:)   = coord';

!      data(coord(1), coord(2) , coord(3), coord(4),coord(5))   = val ;
!    end
!    fclose(fid) ;
!    end

    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
   !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n, rank

    real(kind=rk), intent(in)       :: data_block(:, :, :, :, :)
    character(len=*), intent(in)  :: fileName
    character(len=3) :: rankname
    integer                         :: i,j,k,l,m

write(rankname,'(i3.3)') rank
    open(unit=11, file= fileName//'.'//rankname//'.dat', form='unformatted', status='replace',access='stream')

    write(*,*) "dumping hvy_n=", hvy_n, "rank=", rank
    !write(*,*) size(data_block,1), size(data_block,2), size(data_block,3), size(data_block,4), hvy_n

    write(11) size(data_block,1), size(data_block,2), size(data_block,3), size(data_block,4), hvy_n

    ! loop sequence not very quick, but i prefere this sequence
    do m = 1, hvy_n
        do i = 1, size(data_block,1)
            do j = 1, size(data_block,2)
                do k = 1, size(data_block,3)
                    do l = 1, size(data_block,4)
                        write(11) i, j, k, l, m, data_block(i, j, k, l, hvy_active(m) )
                    end do
                end do
            end do
        end do
    end do

    close(11)
end subroutine
