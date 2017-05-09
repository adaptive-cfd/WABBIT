!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name check_ghost_nodes.f90
!> \version 0.4
!> \author msr
!
!> \brief check if ghost nodes synchronization fails
!
!> 
!!first check: if synchroniozation fails, ghost nodes can have values of 9e9, so check if value larger than given value \n
!! 
!!
!! input:    - params, heavy data, list of active heavy blocks, data field number, value for checking \n
!! output:   -       \n
!! 
!!
!! = log ======================================================================================
!! \n
!! 29/11/16 - create
! ********************************************************************************************

subroutine check_ghost_nodes( params, hvy_block, hvy_active, hvy_n, dF, value)

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params

    !> heavy data array - block data
    real(kind=rk), intent(in)           :: hvy_block(:, :, :, :)

    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    !> datafield to check
    integer(kind=ik), intent(in)        :: dF
    !> value to check
    real(kind=rk), intent(in)           :: value

    ! loop variables
    integer(kind=ik)                    :: k, i, j

    ! grid parameter
    integer(kind=ik)                    :: Bs, g

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    Bs = params%number_block_nodes
    g  = params%number_ghost_nodes

!---------------------------------------------------------------------------------------------
! main body

    ! loop over all active blocks
    do k = 1, hvy_n

        ! loop over ghost nodes
        ! north
        do i = 1, g
            do j = 1, Bs+2*g
                if ( hvy_block( i, j, dF, hvy_active(k) ) > value ) then
                    ! error case
                    print*, "ghost nodes synchronization fails"
                    print*, hvy_active(k)
                    stop
                end if
            end do
        end do

        ! east
        do i = 1, Bs+2*g
            do j = Bs+g+1, Bs+2*g
                if ( hvy_block( i, j, dF, hvy_active(k) ) > value ) then
                    ! error case
                    print*, "ghost nodes synchronization fails"
                    print*, hvy_active(k)
                    stop
                end if
            end do
        end do

        ! south
        do i = Bs+g+1, Bs+2*g
            do j = 1, Bs+2*g
                if ( hvy_block( i, j, dF, hvy_active(k) ) > value ) then
                    ! error case
                    print*, "ghost nodes synchronization fails"
                    print*, hvy_active(k)
                    stop
                end if
            end do
        end do

        ! west
        do i = 1, Bs+2*g
            do j = 1, g
                if ( hvy_block( i, j, dF, hvy_active(k) ) > value ) then
                    ! error case
                    print*, "ghost nodes synchronization fails"
                    print*, hvy_active(k)
                    stop
                end if
            end do
        end do

    end do

end subroutine check_ghost_nodes
