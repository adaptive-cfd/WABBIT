! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: copy_ghost_nodes_3D.f90
! version: 0.5
! author: msr
!
! copy ghost points from sender block to receiver block
! note: works only for block on same process
!
! input:    - heavy data array
!           - sender block id
!           - receiver block id
!           - neighbor relations between sender/receiver
!           - level difference between these two blocks
! output:   - heavy data array
!
! --------------------------------------------------------------------------------------------
! neighbor codes:
! ---------------
! for imagination:  - 6-sided dice with '1'-side on top, '6'-side on bottom, '2'-side in front
!                   - edge: boundary between two sides - use sides numbers for coding
!                   - corner: between three sides - so use all three sides numbers
!                   - block on higher/lower level: block shares face/edge and one unique corner,
!                     so use this corner code in second part of neighbor code
!
! faces:  '__1/___', '__2/___', '__3/___', '__4/___', '__5/___', '__6/___'
! edges:  '_12/___', '_13/___', '_14/___', '_15/___'
!         '_62/___', '_63/___', '_64/___', '_65/___'
!         '_23/___', '_25/___', '_43/___', '_45/___'
! corner: '123/___', '134/___', '145/___', '152/___'
!         '623/___', '634/___', '645/___', '652/___'
!
! complete neighbor code array, 74 possible neighbor relations
! neighbors = (/'__1/___', '__2/___', '__3/___', '__4/___', '__5/___', '__6/___', '_12/___', '_13/___', '_14/___', '_15/___',
!               '_62/___', '_63/___', '_64/___', '_65/___', '_23/___', '_25/___', '_43/___', '_45/___', '123/___', '134/___',
!               '145/___', '152/___', '623/___', '634/___', '645/___', '652/___', '__1/123', '__1/134', '__1/145', '__1/152',
!               '__2/123', '__2/623', '__2/152', '__2/652', '__3/123', '__3/623', '__3/134', '__3/634', '__4/134', '__4/634',
!               '__4/145', '__4/645', '__5/145', '__5/645', '__5/152', '__5/652', '__6/623', '__6/634', '__6/645', '__6/652',
!               '_12/123', '_12/152', '_13/123', '_13/134', '_14/134', '_14/145', '_15/145', '_15/152', '_62/623', '_62/652',
!               '_63/623', '_63/634', '_64/634', '_64/645', '_65/645', '_65/652', '_23/123', '_23/623', '_25/152', '_25/652',
!               '_43/134', '_43/634', '_45/145', '_45/645' /)
! --------------------------------------------------------------------------------------------
!
! = log ======================================================================================
!
! 31/01/17 - create
!
! ********************************************************************************************

subroutine copy_ghost_nodes_3D( params, hvy_block, sender_id, receiver_id, neighborhood, level_diff)

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(in)                  :: params
    ! heavy data array - block data
    real(kind=rk), intent(inout)                    :: hvy_block(:, :, :, :, :)
    ! heavy data id's
    integer(kind=ik), intent(in)                    :: sender_id, receiver_id
    ! neighborhood relation, id from dirs
    integer(kind=ik), intent(in)                    :: neighborhood
    ! difference between block levels
    integer(kind=ik), intent(in)                    :: level_diff

    ! grid parameter
    integer(kind=ik)                                :: Bs, g
    ! loop variables
    integer(kind=ik)                                :: dF

    ! interpolation variables
    !real(kind=rk), dimension(:,:), allocatable      :: data_corner, data_corner_fine, data_edge, data_edge_fine

    ! allocation error variable
    !integer(kind=ik)                                :: allocate_error

!---------------------------------------------------------------------------------------------
! interfaces

    ! grid parameter
    Bs    = params%number_block_nodes
    g     = params%number_ghost_nodes

!---------------------------------------------------------------------------------------------
! variables initialization

!    allocate( data_corner( g, g), stat=allocate_error )
!    call check_allocation(allocate_error)
!    allocate( data_corner_fine( 2*g-1, 2*g-1), stat=allocate_error )
!    call check_allocation(allocate_error)
!    allocate( data_edge( (Bs+1)/2 + g/2, (Bs+1)/2 + g/2), stat=allocate_error )
!    call check_allocation(allocate_error)
!    allocate( data_edge_fine( Bs+g, Bs+g), stat=allocate_error )
!    call check_allocation(allocate_error)

!---------------------------------------------------------------------------------------------
! main body

    select case(neighborhood)
        ! '__1/___'
        case(1)
            if ( level_diff == 0 ) then
                ! sender/receiver on same level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    hvy_block( g+1:Bs+g, g+1:Bs+g, 1:g, dF, receiver_id ) = hvy_block( g+1:Bs+g, g+1:Bs+g, Bs:Bs-1+g, dF, sender_id )
                end do

            else
                ! error case (there should be no level difference between sender/receiver)
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '__2/___'
        case(2)
            if ( level_diff == 0 ) then
                ! sender/receiver on same level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    hvy_block( g+1:Bs+g, Bs+g+1:Bs+g+g, g+1:Bs+g, dF, receiver_id ) = hvy_block( g+1:Bs+g, g+2:g+1+g, g+1:Bs+g, dF, sender_id )
                end do

            else
                ! error case (there should be no level difference between sender/receiver)
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '__3/___'
        case(3)
            if ( level_diff == 0 ) then
                ! sender/receiver on same level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    hvy_block( 1:g, g+1:Bs+g, g+1:Bs+g, dF, receiver_id ) = hvy_block( Bs:Bs-1+g, g+1:Bs+g, g+1:Bs+g, dF, sender_id )
                end do

            else
                ! error case (there should be no level difference between sender/receiver)
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '__4/___'
        case(4)
            if ( level_diff == 0 ) then
                ! sender/receiver on same level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    hvy_block( g+1:Bs+g, 1:g, g+1:Bs+g, dF, receiver_id ) = hvy_block( g+1:Bs+g, Bs:Bs-1+g, g+1:Bs+g, dF, sender_id )
                end do

            else
                ! error case (there should be no level difference between sender/receiver)
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '__5/___'
        case(5)
            if ( level_diff == 0 ) then
                ! sender/receiver on same level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+g, g+1:Bs+g, dF, receiver_id ) = hvy_block( g+2:g+1+g, g+1:Bs+g, g+1:Bs+g, dF, sender_id )
                end do

            else
                ! error case (there should be no level difference between sender/receiver)
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '__6/___'
        case(6)
            if ( level_diff == 0 ) then
                ! sender/receiver on same level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    hvy_block( g+1:Bs+g, g+1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g, g+1:Bs+g, g+2:g+1+g, dF, sender_id )
                end do

            else
                ! error case (there should be no level difference between sender/receiver)
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '_12/___'
        case(7)
            if ( level_diff == 0 ) then
                ! blocks on same level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    hvy_block( g+1:Bs+g, Bs+g+1:Bs+g+g, 1:g, dF, receiver_id ) = hvy_block( g+1:Bs+g, g+2:g+1+g, Bs:Bs-1+g, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '_13/___'
        case(8)
            if ( level_diff == 0 ) then
                ! blocks on same level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    hvy_block( 1:g, g+1:Bs+g, 1:g, dF, receiver_id ) = hvy_block( Bs:Bs-1+g, g+1:Bs+g, Bs:Bs-1+g, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '_14/___'
        case(9)
            if ( level_diff == 0 ) then
                ! blocks on same level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    hvy_block( g+1:Bs+g, 1:g, 1:g, dF, receiver_id ) = hvy_block( g+1:Bs+g, Bs:Bs-1+g, Bs:Bs-1+g, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '_15/___'
        case(10)
            if ( level_diff == 0 ) then
                ! blocks on same level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+g, 1:g, dF, receiver_id ) = hvy_block( g+2:g+1+g, g+1:Bs+g, Bs:Bs-1+g, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '_62/___'
        case(11)
            if ( level_diff == 0 ) then
                ! blocks on same level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    hvy_block( g+1:Bs+g, Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g, g+2:g+1+g, g+2:g+1+g, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '_63/___'
        case(12)
            if ( level_diff == 0 ) then
                ! blocks on same level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    hvy_block( 1:g, g+1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( Bs:Bs-1+g, g+1:Bs+g, g+2:g+1+g, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '_64/___'
        case(13)
            if ( level_diff == 0 ) then
                ! blocks on same level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    hvy_block( g+1:Bs+g, 1:g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g, Bs:Bs-1+g, g+2:g+1+g, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '_65/___'
        case(14)
            if ( level_diff == 0 ) then
                ! blocks on same level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+2:g+1+g, g+1:Bs+g, g+2:g+1+g, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '_23/___'
        case(15)
            if ( level_diff == 0 ) then
                ! blocks on same level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    hvy_block( 1:g, Bs+g+1:Bs+g+g, g+1:Bs+g, dF, receiver_id ) = hvy_block( Bs:Bs-1+g, g+2:g+1+g, g+1:Bs+g, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '_25/___'
        case(16)
            if ( level_diff == 0 ) then
                ! blocks on same level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    hvy_block( Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, g+1:Bs+g, dF, receiver_id ) = hvy_block( g+2:g+1+g, g+2:g+1+g, g+1:Bs+g, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '_43/___'
        case(17)
            if ( level_diff == 0 ) then
                ! blocks on same level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    hvy_block( 1:g, 1:g, g+1:Bs+g, dF, receiver_id ) = hvy_block( Bs:Bs-1+g, Bs:Bs-1+g, g+1:Bs+g, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '_45/___'
        case(18)
            if ( level_diff == 0 ) then
                ! blocks on same level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    hvy_block( Bs+g+1:Bs+g+g, 1:g, g+1:Bs+g, dF, receiver_id ) = hvy_block( g+2:g+1+g, Bs:Bs-1+g, g+1:Bs+g, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '123/___'
        case(19)
            if ( level_diff == 0 ) then
                ! blocks on same level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    hvy_block( 1:g, Bs+g+1:Bs+g+g, 1:g, dF, receiver_id ) = hvy_block( Bs:Bs-1+g, g+2:g+1+g, Bs:Bs-1+g, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '134/___'
        case(20)
            if ( level_diff == 0 ) then
                ! blocks on same level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    hvy_block( 1:g, 1:g, 1:g, dF, receiver_id ) = hvy_block( Bs:Bs-1+g, Bs:Bs-1+g, Bs:Bs-1+g, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '145/___'
        case(21)
            if ( level_diff == 0 ) then
                ! blocks on same level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    hvy_block( Bs+g+1:Bs+g+g, 1:g, 1:g, dF, receiver_id ) = hvy_block( g+2:g+1+g, Bs:Bs-1+g, Bs:Bs-1+g, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '152/___'
        case(22)
            if ( level_diff == 0 ) then
                ! blocks on same level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    hvy_block( Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, 1:g, dF, receiver_id ) = hvy_block( g+2:g+1+g, g+2:g+1+g, Bs:Bs-1+g, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '623/___'
        case(23)
            if ( level_diff == 0 ) then
                ! blocks on same level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    hvy_block( 1:g, Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( Bs:Bs-1+g, g+2:g+1+g, g+2:g+1+g, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '634/___'
        case(24)
            if ( level_diff == 0 ) then
                ! blocks on same level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    hvy_block( 1:g, 1:g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( Bs:Bs-1+g, Bs:Bs-1+g, g+2:g+1+g, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '645/___'
        case(25)
            if ( level_diff == 0 ) then
                ! blocks on same level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    hvy_block( Bs+g+1:Bs+g+g, 1:g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+2:g+1+g, Bs:Bs-1+g, g+2:g+1+g, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '652/___'
        case(26)
            if ( level_diff == 0 ) then
                ! blocks on same level
                ! loop over all datafields
                do dF = 2, params%number_data_fields+1
                    hvy_block( Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+2:g+1+g, g+2:g+1+g, g+2:g+1+g, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

    end select

    ! clean up
!    deallocate( data_corner, stat=allocate_error )
!    deallocate( data_corner_fine, stat=allocate_error )
!    deallocate( data_edge, stat=allocate_error )
!    deallocate( data_edge_fine, stat=allocate_error )

end subroutine copy_ghost_nodes_3D
