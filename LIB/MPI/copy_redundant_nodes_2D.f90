!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name copy_redundant_nodes_2D.f90
!> \version 0.5
!> \author msr
!
!> \brief copy redundant nodes from sender block to receiver block
!> \note works only for block on same process
!
!> \details
!! input:    
!!           - heavy data array
!!           - sender block id
!!           - receiver block id
!!           - neighbor relations between sender/receiver
!!           - level difference between these two blocks
!!
!! output:   
!!           - heavy data array
!!
!! \n
! -------------------------------------------------------------------------------------------------------------------------
!> dirs = (/'__N', '__E', '__S', '__W', '_NE', '_NW', '_SE', '_SW', 'NNE', 'NNW', 'SSE', 'SSW', 'ENE', 'ESE', 'WNW', 'WSW'/) \n
! -------------------------------------------------------------------------------------------------------------------------
!>
!!
!! = log ======================================================================================
!! \n
!! 22/05/17 - create
!!
! ********************************************************************************************

subroutine copy_redundant_nodes_2D( params, hvy_block, sender_id, receiver_id, neighborhood, level_diff, hvy_neighbor)

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)                  :: params
    !> heavy data array - block data
    real(kind=rk), intent(inout)                    :: hvy_block(:, :, :, :)
    !> heavy data id's
    integer(kind=ik), intent(in)                    :: sender_id, receiver_id
    !> neighborhood relation, id from dirs
    integer(kind=ik), intent(in)                    :: neighborhood
    !> difference between block levels
    integer(kind=ik), intent(in)                    :: level_diff

    !> heavy data array - neifghbor data
    integer(kind=ik), intent(in)                    :: hvy_neighbor(:,:)

    ! grid parameter
    integer(kind=ik)                                :: Bs, g
    ! loop variables
    integer(kind=ik)                                :: dF

    ! light data id
    integer(kind=ik)                                :: sender_lgt_id, receiver_lgt_id

!---------------------------------------------------------------------------------------------
! interfaces

    ! grid parameter
    Bs    = params%number_block_nodes
    g     = params%number_ghost_nodes

    ! light id
    call hvy_id_to_lgt_id( sender_lgt_id, sender_id, params%rank, params%number_blocks )
    call hvy_id_to_lgt_id( receiver_lgt_id, receiver_id, params%rank, params%number_blocks )

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    select case(neighborhood)
        ! '__N'
        case(1)

        ! '__E'
        case(2)

        ! '__S'
        case(3)

        ! '__W'
        case(4)

        ! '_NE'
        case(5)
            ! case d)
            ! corner neighbor on same level
            ! side neighbor (at least one of them) one level down
            ! => synchronize with redundant nodes
            ! NE: relevant sides N, E

            ! check if blocks on same level
            if ( level_diff == 0 ) then
                ! check if one side neighbor is one level down
                ! check existence of half-side neighbor
                if ( (hvy_neighbor( sender_id, 9) /= -1) .or. (hvy_neighbor( sender_id, 13) /= -1) ) then
                    ! loop over all datafields
                    do dF = 1, params%number_data_fields
                        hvy_block( Bs+g:Bs+g+g, 1:g+1, dF, receiver_id ) = hvy_block( g+1:g+1+g, Bs:Bs+g, dF, sender_id )
                    end do
                end if
            end if

        ! '_NW'
        case(6)
            ! case d)
            ! corner neighbor on same level
            ! side neighbor (at least one of them) one level down
            ! => synchronize with redundant nodes
            ! NW: relevant sides N, W

            ! check if blocks on same level
            if ( level_diff == 0 ) then
                ! check if one side neighbor is one level down
                ! check existence of half-side neighbor
                if ( (hvy_neighbor( sender_id, 10) /= -1) .or. (hvy_neighbor( sender_id, 15) /= -1) ) then
                    ! loop over all datafields
                    do dF = 1, params%number_data_fields
                        hvy_block( Bs+g:Bs+g+g, Bs+g:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:g+1+g, g+1:g+1+g, dF, sender_id )
                    end do
                end if
            end if

        ! '_SE'
        case(7)
            ! case d)
            ! corner neighbor on same level
            ! side neighbor (at least one of them) one level down
            ! => synchronize with redundant nodes
            ! SE: relevant sides S, E

            ! check if blocks on same level
            if ( level_diff == 0 ) then
                ! check if one side neighbor is one level down
                ! check existence of half-side neighbor
                if ( (hvy_neighbor( sender_id, 11) /= -1) .or. (hvy_neighbor( sender_id, 14) /= -1) ) then
                    ! loop over all datafields
                    do dF = 1, params%number_data_fields
                        hvy_block( 1:g+1, 1:g+1, dF, receiver_id ) = hvy_block( Bs:Bs+g, Bs:Bs+g, dF, sender_id )
                    end do
                end if
            end if

        ! '_SW'
        case(8)
            ! case d)
            ! corner neighbor on same level
            ! side neighbor (at least one of them) one level down
            ! => synchronize with redundant nodes
            ! SW: relevant sides S, W

            ! check if blocks on same level
            if ( level_diff == 0 ) then
                ! check if one side neighbor is one level down
                ! check existence of half-side neighbor
                if ( (hvy_neighbor( sender_id, 12) /= -1) .or. (hvy_neighbor( sender_id, 16) /= -1) ) then
                    ! loop over all datafields
                    do dF = 1, params%number_data_fields
                        hvy_block( 1:g+1, Bs+g:Bs+g+g, dF, receiver_id ) = hvy_block( Bs:Bs+g, g+1:g+1+g, dF, sender_id )
                    end do
                end if
            end if

        ! 'NNE'
        case(9)

        ! 'NNW'
        case(10)

        ! 'SSE'
        case(11)

        ! 'SSW'
        case(12)

        ! 'ENE'
        case(13)

        ! 'ESE'
        case(14)

        ! 'WNW'
        case(15)

        ! 'WSW'
        case(16)

    end select

end subroutine copy_redundant_nodes_2D
