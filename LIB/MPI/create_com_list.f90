! ********************************
! WABBIT
! --------------------------------
!
! creates communication list
!
! name: create_com_list.f90
! date: 26.10.2016
! author: msr
! version: 0.3
!
! ********************************

subroutine create_com_list(com_list, com_list_i, com_list_N)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), intent(in)                            :: com_list_N

    integer(kind=ik), dimension(com_list_N,7), intent(out)  :: com_list
    integer(kind=ik), intent(out)                           :: com_list_i
    integer                                                 :: find_neighborhood

    integer(kind=ik)                                        :: light_id, neighbor_light_id, N, l

    com_list    = -99
    com_list_i  = 1
    N           = blocks_params%number_max_blocks

    ! loop over all blocks to create com-list
    do light_id = 1, N

        if ( blocks(light_id)%active == .true. ) then

            if ( blocks_params%size_domain == blocks_params%size_block ) then
                ! write com_list if only one block
                do l = 1, 8
                    com_list(com_list_i, 1)  = com_list_i
                    com_list(com_list_i, 2)  = blocks(light_id)%proc_rank
                    com_list(com_list_i, 3)  = blocks(light_id)%proc_rank
                    com_list(com_list_i, 4)  = light_id
                    com_list(com_list_i, 5)  = light_id
                    com_list(com_list_i, 6)  = l
                    com_list(com_list_i, 7)  = find_neighborhood(blocks(light_id)%neighbor_dir(l))
                    com_list_i               = com_list_i + 1
                end do

            else
                ! more than one block
                ! loop over all directions
                do l = 1, 16
                    if ( blocks(light_id)%neighbor_id(l) /= -1 ) then
                        ! neighbor light id
                        neighbor_light_id  = blocks(light_id)%neighbor_id(l)
                        ! write list
                        com_list(com_list_i, 1)  = com_list_i
                        com_list(com_list_i, 2)  = blocks(light_id)%proc_rank
                        com_list(com_list_i, 3)  = blocks(neighbor_light_id)%proc_rank
                        com_list(com_list_i, 4)  = light_id
                        com_list(com_list_i, 5)  = neighbor_light_id
                        com_list(com_list_i, 6)  = l
                        com_list(com_list_i, 7)  = find_neighborhood(blocks(light_id)%neighbor_dir(l))
                        com_list_i               = com_list_i + 1
                    end if
                end do

            end if

        end if

    end do

    com_list_i = com_list_i - 1

end subroutine create_com_list
