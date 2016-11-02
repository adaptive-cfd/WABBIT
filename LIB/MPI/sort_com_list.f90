! ********************************
! WABBIT
! --------------------------------
!
! sort com_list
!
! name: sort_com_list.f90
! date: 26.10.2016
! author: msr
! version: 0.3
!
! ********************************

subroutine sort_com_list(com_list, com_plan, n_proc, n_com, com_list_N)

    use module_params
    use module_blocks

    implicit none

    integer(kind=ik), intent(in)                                :: com_list_N

    integer(kind=ik), dimension(com_list_N, 7), intent(inout)   :: com_list
    integer, dimension(2000, 2), intent(out)                    :: com_plan
    integer(kind=ik), intent(in)                                :: n_proc, n_com

    integer                                                     :: i, j, k, allocate_error, com_count, plan_type
    integer, dimension(7)                                       :: com_list_elem

    logical, dimension(:), allocatable                          :: proc_status
    logical                                                     :: com_allowed

    ! reset com_plan
    com_plan    = -99

    allocate( proc_status(n_proc), stat=allocate_error )
    ! at start all procs can communicate
    proc_status = .true.

    ! create communication plan and sort com_list
    i = 1
    j = 1

    do while ( com_list(i,1) /= -99 )

        com_count = 0

        ! first: check if next procs in list can communicate
        if ( proc_status(com_list(i,2)+1) .and. proc_status(com_list(i,3)+1) ) then
            ! both procs can communicate
            ! one additional communication
            com_count = com_count + 1
            ! set status for current proc and target proc to .false.
            proc_status( com_list(i,2)+1 ) = .false.
            proc_status( com_list(i,3)+1 ) = .false.
        else
            ! sort list, switch to next allowed communciation
            do k = (i+1), n_com
                if ( proc_status( com_list(k,2)+1 ) ) then
                    if ( proc_status( com_list(k,3)+1 ) ) then

                        ! sort list
                        com_list_elem       = com_list(k,:)
                        com_list(i+1:k, :)  = com_list(i:k-1, :)
                        com_list(i, :)      = com_list_elem
                        ! increase number of communications
                        com_count = com_count + 1
                        ! set status for current proc and target proc to .false.
                        proc_status( com_list(i,2)+1 ) = .false.
                        proc_status( com_list(i,3)+1 ) = .false.
                        exit

                    end if
                end if
            end do
        end if

        ! second step only if at least one allowed communication
        if ( com_count == 1 ) then
            ! second: sort all similar communications between these procs on top of the list
            do k = (i+1), n_com
                if ( (com_list(k,2) == com_list(i,2)) .and. (com_list(k,3) == com_list(i,3)) ) then
                    if (k == (i+1)) then
                        ! increase number of communications
                        com_count = com_count + 1
                    else
                        ! sort list, switch next communciation with allowed communication
                        com_list_elem       = com_list(k,:)
                        com_list(i+2:k, :)  = com_list(i+1:k-1, :)
                        com_list(i+1, :)    = com_list_elem
                        ! increase number of communications
                        com_count = com_count + 1
                    end if
                end if
            end do
            ! third: sort all reverse communications
            do k = (i+com_count), n_com
                if ( (com_list(k,2) == com_list(i,3)) .and. (com_list(k,3) == com_list(i,2)) ) then
                    if (k == (i+com_count)) then
                        ! nothing to do
                    else
                        ! sort list, switch next communciation with allowed communication
                        com_list_elem                 = com_list(k,:)
                        com_list(i+1+com_count:k, :)  = com_list(i+com_count:k-1, :)
                        com_list(i+com_count, :)      = com_list_elem
                    end if
                end if
            end do
        end if

        ! set type off communication (1: internal, 2: external)
        if ( com_list(i,2) == com_list(i,3) ) then
            plan_type = 1
        else
            plan_type = 2
        end if

        ! move list index
        if (plan_type == 1) then
            ! internal com
            i = i + com_count
        else
            ! external com
            i = i + 2*com_count
        end if

        if ( com_count > 0) then

            if (plan_type == 1) then
                ! set com plan, internal com
                com_plan(j, 1)    = plan_type
                com_plan(j, 2)    = com_count
                ! prepare for next loop
                j = j + 1
            else
                ! set com plan, external com
                com_plan(j, 1)    = plan_type
                com_plan(j, 2)    = com_count
                com_plan(j+1, 1)  = plan_type
                com_plan(j+1, 2)  = com_count
                ! prepare for next loop
                j = j + 2
            end if
        end if

        ! reset proc_status, if necessary
        if ( com_allowed(proc_status, n_proc) == .false. ) then
            proc_status = .true.
        end if

        ! no communications found and not reached end of com_list:
        ! reset proc_status
        if ( (com_count == 0) .and. (i <= n_com) ) then
            proc_status = .true.
        end if

    end do

    deallocate( proc_status, stat=allocate_error )

end subroutine sort_com_list
