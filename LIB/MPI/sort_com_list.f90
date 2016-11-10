! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: sort_com_list.f90
! version: 0.4
! author: msr
!
! sort com_list
!
! input:    - com_list
!           - number of procs
!           - number of communications
! output:   - sorted com_list
!           - com_plan
!
! = log ======================================================================================
!
! 08/11/16 - switch to v0.4
! ********************************************************************************************

subroutine sort_com_list(com_list, com_list_N, com_plan, com_plan_N, n_proc, n_com)

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! com_list length
    integer(kind=ik), intent(in)          :: com_list_N
    ! com list
    integer(kind=ik), intent(inout)       :: com_list(com_list_N, 8)

    ! com_plan length
    integer(kind=ik), intent(in)          :: com_plan_N
    ! com plan
    integer, intent(out)                  :: com_plan(com_plan_N)

    ! number of procs, number of communications
    integer(kind=ik), intent(in)          :: n_proc, n_com

    ! allocation error variable
    integer(kind=ik)                      :: allocate_error

    ! loop variables
    integer(kind=ik)                      :: i, j, k

    ! number of communications for com_plan
    integer(kind=ik)                      :: com_count

    ! line of com_list to copy/sort list
    integer, dimension(8)                 :: com_list_elem

    ! status of proc: .true. - can send/receive, .false. - con not send/receive
    logical, dimension(:), allocatable    :: proc_status
    ! .true. - there are communications between procs left - use to parallelize communication
    logical                               :: com_allowed

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! reset com_plan
    com_plan    = -1

    allocate( proc_status(n_proc), stat=allocate_error )
    ! at start all procs can communicate
    proc_status = .true.

!---------------------------------------------------------------------------------------------
! main body

    ! create communication plan and sort com_list
    i = 1
    j = 1

    do while ( com_list(i,1) /= -1 )

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

        ! move list index
        i = i + 2*com_count

        if ( com_count > 0) then

            ! set com plan, external com
            com_plan(j)    = com_count
            com_plan(j+1)  = com_count
            ! prepare for next loop
            j = j + 2

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
