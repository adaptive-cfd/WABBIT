! ********************************
! 2D AMR prototype
! --------------------------------
!
! synchronize ghosts nodes
!
! name: synchronize_ghosts.f90
! date: 05.08.2016
! author: msr
! version: 0.1
!
! ********************************

subroutine synchronize_ghosts()

    use module_params
    use module_blocks
    use module_interpolation

    implicit none

    integer(kind=ik)	                            :: dF, k, N, block_num, i, neighbor_num, neighbor2_num, g, Bs, allocate_error
    real(kind=rk), dimension(:,:), allocatable      :: data_predict_fine, data_predict_coarse, neighbor_data, neighbor2_data

     N  = size(blocks_params%active_list, dim=1)
     g  = blocks_params%number_ghost_nodes
     Bs = blocks_params%size_block

     allocate( data_predict_fine(2*Bs-1, 2*Bs-1), stat=allocate_error )
     allocate( data_predict_coarse(Bs, Bs), stat=allocate_error )
     allocate( neighbor_data(Bs, Bs), stat=allocate_error )
     allocate( neighbor2_data(Bs, Bs), stat=allocate_error )

     !loop over all data fields
     do dF = 1, blocks_params%number_data_fields

         ! loop over all active blocks
         do k = 1, N

             block_num = blocks_params%active_list(k)

             ! reset ghost nodes
             blocks(block_num)%data_fields(dF)%data_(1:g, : )           = 9.0e9_rk
             blocks(block_num)%data_fields(dF)%data_(Bs+g+1:Bs+2*g, : ) = 9.0e9_rk
             blocks(block_num)%data_fields(dF)%data_(:, 1:g )           = 9.0e9_rk
             blocks(block_num)%data_fields(dF)%data_(:, Bs+g+1:Bs+2*g ) = 9.0e9_rk

             ! loop over all possible neighbors (max=8)
             do i = 1, 8

                 ! neighbor exists, first neighbor list (include neighbors on same level)
                 if (( blocks(block_num)%neighbor_id(i) /= -1 )) then

                     neighbor_num  = blocks(block_num)%neighbor_id(i)
                     ! copy neighbor data for readability
                     neighbor_data = blocks(neighbor_num)%data_fields(dF)%data_(g+1:Bs+g, g+1:Bs+g)

                     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     ! neighbor is on the same level
                     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                     if ( blocks(block_num)%level == blocks(neighbor_num)%level ) then
                         select case(blocks(block_num)%neighbor_dir(i))
                             ! first line/row overlap, so do not copy!
                             !north
                             case('NO')
                             blocks(block_num)%data_fields(dF)%data_(1:g, g+1:Bs+g)                 = neighbor_data(Bs-g:Bs-1, :)

                             !south
                             case('SO')
                             blocks(block_num)%data_fields(dF)%data_(Bs+g+1:Bs+g+g, g+1:Bs+g)       = neighbor_data(2:g+1, :)

                             !east
                             case('EA')
                             blocks(block_num)%data_fields(dF)%data_(g+1:Bs+g, Bs+g+1:Bs+g+g)       = neighbor_data(:, 2:g+1)

                             !west
                             case('WE')
                             blocks(block_num)%data_fields(dF)%data_(g+1:Bs+g, 1:g)                 = neighbor_data(:, Bs-g:Bs-1)

                             !northeast
                             case('NE')
                             blocks(block_num)%data_fields(dF)%data_(1:g, Bs+g+1:Bs+g+g)            = neighbor_data(Bs-g:Bs-1, 2:g+1)

                             !northwest
                             case('NW')
                             blocks(block_num)%data_fields(dF)%data_(1:g, 1:g)                      = neighbor_data(Bs-g:Bs-1, Bs-g:Bs-1)

                             !southeast
                             case('SE')
                             blocks(block_num)%data_fields(dF)%data_(Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g)  = neighbor_data(2:g+1, 2:g+1)

                             !southwest
                             case('SW')
                             blocks(block_num)%data_fields(dF)%data_(Bs+g+1:Bs+g+g, 1:g)            = neighbor_data(2:g+1, Bs-g:Bs-1)
                         end select

                     elseif ( blocks(block_num)%level > blocks(neighbor_num)%level) then
                         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                         ! neighbor is one level coarser
                         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                         ! predict finer data
                         data_predict_fine   = 0.0_rk
                         data_predict_coarse = neighbor_data
                         call prediction_2D(data_predict_coarse, data_predict_fine)

                         select case(blocks(block_num)%neighbor_dir(i))

                             ! first line/row overlap, so do not copy!
                             ! north
                             case('NO')
                             if ( blocks(block_num)%treecode(blocks(block_num)%level) == 0 ) then
                                 ! block is on "0"-side of neighbor
                                 blocks(block_num)%data_fields(dF)%data_(1:g, g+1:Bs+g)                = data_predict_fine(2*Bs-1-g:2*Bs-2, 1:Bs)
                             else
                                 ! block is on "1"-side of neighbor
                                 blocks(block_num)%data_fields(dF)%data_(1:g, g+1:Bs+g)                = data_predict_fine(2*Bs-1-g:2*Bs-2, Bs:2*Bs-1)
                             end if

                             !south
                             case('SO')
                             if ( blocks(block_num)%treecode(blocks(block_num)%level) == 2 ) then
                                 ! block is on "2"-side of neighbor
                                 blocks(block_num)%data_fields(dF)%data_(Bs+g+1:Bs+g+g, g+1:Bs+g)     = data_predict_fine(2:g+1, 1:Bs)
                             else
                                 ! block is on "3"-side of neighbor
                                 blocks(block_num)%data_fields(dF)%data_(Bs+g+1:Bs+g+g, g+1:Bs+g)     = data_predict_fine(2:g+1, Bs:2*Bs-1)
                             end if

                             !east
                             case('EA')
                             if ( blocks(block_num)%treecode(blocks(block_num)%level) == 1 ) then
                                 ! block is on "1"-side of neighbor
                                 blocks(block_num)%data_fields(dF)%data_(g+1:Bs+g, Bs+g+1:Bs+g+g)      = data_predict_fine(1:Bs, 2:g+1)
                             else
                                 ! block is on "3"-side of neighbor
                                 blocks(block_num)%data_fields(dF)%data_(g+1:Bs+g, Bs+g+1:Bs+g+g)      = data_predict_fine(Bs:2*Bs-1, 2:g+1)
                             end if

                             !west
                             case('WE')
                             if ( blocks(block_num)%treecode(blocks(block_num)%level) == 0 ) then
                                 ! block is on "0"-side of neighbor
                                 blocks(block_num)%data_fields(dF)%data_(g+1:Bs+g, 1:g)                = data_predict_fine(1:Bs, 2*Bs-1-g:2*Bs-2)
                             else
                                 ! block is on "2"-side of neighbor
                                 blocks(block_num)%data_fields(dF)%data_(g+1:Bs+g, 1:g)                = data_predict_fine(Bs:2*Bs-1, 2*Bs-1-g:2*Bs-2)
                             end if

                             !northeast
                             case('NE')
                             select case(blocks(block_num)%treecode(blocks(block_num)%level))
                             ! check on which side of coarser neighbor current block is located
                                 case(0)
                                 blocks(block_num)%data_fields(dF)%data_(1:g, Bs+g+1:Bs+g+g)           = data_predict_fine(2*Bs-1-g:2*Bs-2, Bs+1:Bs+g)
                                 case(1)
                                 blocks(block_num)%data_fields(dF)%data_(1:g, Bs+g+1:Bs+g+g)           = data_predict_fine(2*Bs-1-g:2*Bs-2, 2:g+1)
                                 case(3)
                                 blocks(block_num)%data_fields(dF)%data_(1:g, Bs+g+1:Bs+g+g)           = data_predict_fine(Bs-g:Bs-1, 2:g+1)
                             end select

                             !northwest
                             case('NW')
                             select case(blocks(block_num)%treecode(blocks(block_num)%level))
                             ! check on which side of coarser neighbor current block is located
                                 case(0)
                                 blocks(block_num)%data_fields(dF)%data_(1:g, 1:g)                     = data_predict_fine(2*Bs-1-g:2*Bs-2, 2*Bs-1-g:2*Bs-2)
                                 case(1)
                                 blocks(block_num)%data_fields(dF)%data_(1:g, 1:g)                     = data_predict_fine(2*Bs-1-g:2*Bs-2, Bs-g:Bs-1)
                                 case(2)
                                 blocks(block_num)%data_fields(dF)%data_(1:g, 1:g)                     = data_predict_fine(Bs-g:Bs-1, 2*Bs-1-g:2*Bs-2)
                             end select

                             !southeast
                             case('SE')
                             select case(blocks(block_num)%treecode(blocks(block_num)%level))
                             ! check on which side of coarser neighbor current block is located
                                 case(1)
                                 blocks(block_num)%data_fields(dF)%data_(Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g) = data_predict_fine(Bs+1:Bs+g, 2:g+1)
                                 case(2)
                                 blocks(block_num)%data_fields(dF)%data_(Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g) = data_predict_fine(2:g+1, Bs+1:Bs+g)
                                 case(3)
                                 blocks(block_num)%data_fields(dF)%data_(Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g) = data_predict_fine(2:g+1, 2:g+1)
                             end select

                             !southwest
                             case('SW')
                             select case(blocks(block_num)%treecode(blocks(block_num)%level))
                             ! check on which side of coarser neighbor current block is located
                                 case(0)
                                 blocks(block_num)%data_fields(dF)%data_(Bs+g+1:Bs+g+g, 1:g)           = data_predict_fine(Bs+1:Bs+g, 2*Bs-1-g:2*Bs-2)
                                 case(2)
                                 blocks(block_num)%data_fields(dF)%data_(Bs+g+1:Bs+g+g, 1:g)           = data_predict_fine(2:g+1, 2*Bs-1-g:2*Bs-2)
                                 case(3)
                                 blocks(block_num)%data_fields(dF)%data_(Bs+g+1:Bs+g+g, 1:g)           = data_predict_fine(2:g+1, Bs-g:Bs-1)
                             end select

                         end select

                     elseif ( blocks(block_num)%neighbor_number(i) == 1 ) then
                         ! there is exact one neighbor, one level finer, only valid for diagonal neighbors
                         select case(blocks(block_num)%neighbor_dir(i))
                             !northeast
                             case('NE')
                             blocks(block_num)%data_fields(dF)%data_(1:g, Bs+g+1:Bs+g+g)           = neighbor_data(Bs-2*g:Bs-2:2, 3:2*g+1:2)

                             !northwest
                             case('NW')
                             blocks(block_num)%data_fields(dF)%data_(1:g, 1:g)                     = neighbor_data(Bs-2*g:Bs-2:2, Bs-2*g:Bs-2:2)

                             !southeast
                             case('SE')
                             blocks(block_num)%data_fields(dF)%data_(Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g) = neighbor_data(3:2*g+1:2, 3:2*g+1:2)

                             !southwest
                             case('SW')
                             blocks(block_num)%data_fields(dF)%data_(Bs+g+1:Bs+g+g, 1:g)           = neighbor_data(3:2*g+1:2, Bs-2*g:Bs-2:2)

                         end select

                     else
                         ! there are two neighbors, neighbor is located on an edge
                         ! first neighbor id is known
                         ! assume second neighbor id is saved on same element in second neighbor list as first id
                         neighbor2_num  = blocks(block_num)%neighbor2_id(i)
                         neighbor2_data = blocks(neighbor2_num)%data_fields(dF)%data_(g+1:Bs+g, g+1:Bs+g)

                         select case(blocks(block_num)%neighbor_dir(i))
                             !north
                             case('NO')
                             if ( blocks(neighbor_num)%treecode(blocks(neighbor_num)%level) == 2 ) then
                                 ! first neighbor in "2" direction
                                 blocks(block_num)%data_fields(dF)%data_(1:g, g+1:(Bs+1)/2+g)          = neighbor_data(Bs-2*g:Bs-2:2, 1:Bs:2)
                                 ! second neighbor in "3" direction
                                 blocks(block_num)%data_fields(dF)%data_(1:g, (Bs+1)/2+g:Bs+g)         = neighbor2_data(Bs-2*g:Bs-2:2, 1:Bs:2)
                             else
                                 ! first neighbor in "3" direction
                                 blocks(block_num)%data_fields(dF)%data_(1:g, (Bs+1)/2+g:Bs+g)         = neighbor_data(Bs-2*g:Bs-2:2, 1:Bs:2)
                                 ! second neighbor in "2" direction
                                 blocks(block_num)%data_fields(dF)%data_(1:g, g+1:(Bs+1)/2+g)          = neighbor2_data(Bs-2*g:Bs-2:2, 1:Bs:2)
                             end if

                             !south
                             case('SO')
                             if ( blocks(neighbor_num)%treecode(blocks(neighbor_num)%level) == 0 ) then
                                 ! first neighbor in "0" direction
                                 blocks(block_num)%data_fields(dF)%data_(Bs+g+1:Bs+g+g, g+1:(Bs+1)/2+g) = neighbor_data(3:2*g+1:2, 1:Bs:2)
                                 ! second neighbor in "1" direction
                                 blocks(block_num)%data_fields(dF)%data_(Bs+g+1:Bs+g+g, (Bs+1)/2+g:Bs+g)= neighbor2_data(3:2*g+1:2, 1:Bs:2)
                             else
                                 ! first neighbor in "1" direction
                                 blocks(block_num)%data_fields(dF)%data_(Bs+g+1:Bs+g+g, (Bs+1)/2+g:Bs+g)= neighbor_data(3:2*g+1:2, 1:Bs:2)
                                 ! second neighbor in "0" direction
                                 blocks(block_num)%data_fields(dF)%data_(Bs+g+1:Bs+g+g, g+1:(Bs+1)/2+g) = neighbor2_data(3:2*g+1:2, 1:Bs:2)
                             end if

                             !west
                             case('WE')
                             if ( blocks(neighbor_num)%treecode(blocks(neighbor_num)%level) == 1 ) then
                                 ! first neighbor in "1" direction
                                 blocks(block_num)%data_fields(dF)%data_(g+1:(Bs+1)/2+g, 1:g)           = neighbor_data(1:Bs:2, Bs-2*g:Bs-2:2)
                                 ! second neighbor in "3" direction
                                 blocks(block_num)%data_fields(dF)%data_((Bs+1)/2+g:Bs+g, 1:g)          = neighbor2_data(1:Bs:2, Bs-2*g:Bs-2:2)
                             else
                                 ! first neighbor in "3" direction
                                 blocks(block_num)%data_fields(dF)%data_((Bs+1)/2+g:Bs+g, 1:g)          = neighbor_data(1:Bs:2, Bs-2*g:Bs-2:2)
                                 ! second neighbor in "1" direction
                                 blocks(block_num)%data_fields(dF)%data_(g+1:(Bs+1)/2+g, 1:g)           = neighbor2_data(1:Bs:2, Bs-2*g:Bs-2:2)
                             end if

                             !east
                             case('EA')
                             if ( blocks(neighbor_num)%treecode(blocks(neighbor_num)%level) == 0 ) then
                                 ! first neighbor in "0" direction
                                 blocks(block_num)%data_fields(dF)%data_(g+1:(Bs+1)/2+g, Bs+g+1:Bs+g+g) = neighbor_data(1:Bs:2, 3:2*g+1:2)
                                 ! second neighbor in "2" direction
                                 blocks(block_num)%data_fields(dF)%data_((Bs+1)/2+g:Bs+g, Bs+g+1:Bs+g+g)= neighbor2_data(1:Bs:2, 3:2*g+1:2)
                             else
                                 ! first neighbor in "2" direction
                                 blocks(block_num)%data_fields(dF)%data_((Bs+1)/2+g:Bs+g, Bs+g+1:Bs+g+g)= neighbor_data(1:Bs:2, 3:2*g+1:2)
                                 ! second neighbor in "0" direction
                                 blocks(block_num)%data_fields(dF)%data_(g+1:(Bs+1)/2+g, Bs+g+1:Bs+g+g) = neighbor2_data(1:Bs:2, 3:2*g+1:2)
                             end if

                         end select

                     end if

                 end if

             end do

         end do

     end do

     deallocate( data_predict_fine, stat=allocate_error )
     deallocate( data_predict_coarse, stat=allocate_error )
     deallocate( neighbor_data, stat=allocate_error )
     deallocate( neighbor2_data, stat=allocate_error )

end subroutine synchronize_ghosts
