! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: compute_affinity.f90
! version: 0.4
! author: engels, msr
!
!-------------------------------------------------------------------------------
! compute the affinity for all my blocks to a given rank.
! That means, if a block has many neighbor relations with the target rank, it has a high
! value in the array, otherwise a low value or -10 if the block is not active at all
!
! This list is at the core of heuristic load balancing: we know whom to send to, and using this list
! we decide which blocks will be sent.
!-------------------------------------------------------------------------------
!
! input:    - params, non-synchronized light data list, neighbor list
! output:   - affinity list
!
! = log ======================================================================================
!
! 28/11/16 - create
!
! ********************************************************************************************
subroutine compute_affinity(params, my_block_list, hvy_neighbor, rank, rank_partner, affinity, hvy_active, hvy_n)

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(in)      :: params

    ! light data array
    integer(kind=ik), intent(inout)     :: my_block_list(:, :)

    ! heavy data array - neifghbor data
    integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)

    ! affinity list
    integer(kind=ik), intent(out)       :: affinity(:)

    ! proc status
    integer(kind=ik), intent(in)        :: rank, rank_partner

    ! list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    ! number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    ! loop variables
    integer                             :: heavy_id, light_id, proc_id, q, k

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! affinity list, HEAVY DATA ARRAY
    affinity = 0

!---------------------------------------------------------------------------------------------
! main body

    ! error case
    if( size(affinity) /= params%number_blocks ) then
        write(*,*)"unforeseen"
        stop
    endif

    ! loop over all active heavy
    do k = 1, hvy_n

        ! heavy id
        heavy_id = hvy_active(k)
        ! get corresponding light_id
        call hvy_id_to_lgt_id( light_id, heavy_id, rank, params%number_blocks )

        if (my_block_list( light_id, 1) == -1) then
            ! inactive block have extremely low affinity to be transferred
            affinity(heavy_id) = -1
        else
            ! loop over all directions and count how many neighbor relations I do have with the rank_partner
            do q = 1, 16 ! q is relative direction
                if (hvy_neighbor( heavy_id, q ) > 0) then ! is there a neighbor?
                    ! proc rank of neighbors id
                    call lgt_id_to_proc_rank( proc_id, hvy_neighbor( heavy_id, q ), params%number_blocks )

                    if (proc_id == rank_partner) then !receiver? ZERO BASED
                        ! a shared border with the target rank is a high priority
                        affinity(heavy_id) = affinity(heavy_id)+20
                    elseif (proc_id /= rank) then
                        ! so I dont share this border with the target rank, but it is an mpi border
                        ! when I'm out of good candidates, I should at least send blocks that are not surrounded only by my blocks
                        affinity(heavy_id) = affinity(heavy_id)+1
                    endif
                endif
            enddo
        end if

    end do

end subroutine compute_affinity
