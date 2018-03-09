!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name copy_ghost_nodes_3D.f90
!> \version 0.5
!> \author msr
!
!> \brief copy ghost points from sender block to receiver block
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
! --------------------------------------------------------------------------------------------
!> neighbor codes: \n
! ---------------
!> for imagination:  
!!                   - 6-sided dice with '1'-side on top, '6'-side on bottom, '2'-side in front
!!                   - edge: boundary between two sides - use sides numbers for coding
!!                   - corner: between three sides - so use all three sides numbers
!!                   - block on higher/lower level: block shares face/edge and one unique corner,
!!                     so use this corner code in second part of neighbor code
!!
!! faces:  '__1/___', '__2/___', '__3/___', '__4/___', '__5/___', '__6/___' \n
!! edges:  '_12/___', '_13/___', '_14/___', '_15/___'
!!         '_62/___', '_63/___', '_64/___', '_65/___'
!!         '_23/___', '_25/___', '_43/___', '_45/___' \n
!! corner: '123/___', '134/___', '145/___', '152/___'
!!         '623/___', '634/___', '645/___', '652/___' \n
!!
!!
!! complete neighbor code array, 74 possible neighbor relations \n
!!
!! neighbors = (/'__1/___', '__2/___', '__3/___', '__4/___', '__5/___', '__6/___', '_12/___', '_13/___', '_14/___', '_15/___',
!!               '_62/___', '_63/___', '_64/___', '_65/___', '_23/___', '_25/___', '_43/___', '_45/___', '123/___', '134/___',
!!               '145/___', '152/___', '623/___', '634/___', '645/___', '652/___', '__1/123', '__1/134', '__1/145', '__1/152',
!!               '__2/123', '__2/623', '__2/152', '__2/652', '__3/123', '__3/623', '__3/134', '__3/634', '__4/134', '__4/634',
!!               '__4/145', '__4/645', '__5/145', '__5/645', '__5/152', '__5/652', '__6/623', '__6/634', '__6/645', '__6/652',
!!               '_12/123', '_12/152', '_13/123', '_13/134', '_14/134', '_14/145', '_15/145', '_15/152', '_62/623', '_62/652',
!!               '_63/623', '_63/634', '_64/634', '_64/645', '_65/645', '_65/652', '_23/123', '_23/623', '_25/152', '_25/652',
!!               '_43/134', '_43/634', '_45/145', '_45/645' /) \n
! --------------------------------------------------------------------------------------------
!>
!!
!! = log ======================================================================================
!! \n
!! 31/01/17 - create \n
!! 12/04/17 - remove redundant nodes between blocks with meshlevel +1
!
! ********************************************************************************************

subroutine copy_ghost_nodes_3D( params, hvy_block, sender_id, receiver_id, neighborhood, level_diff)

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)                  :: params
    !> heavy data array - block data
    real(kind=rk), intent(inout)                    :: hvy_block(:, :, :, :, :)
    !> heavy data id's
    integer(kind=ik), intent(in)                    :: sender_id, receiver_id
    !> neighborhood relation, id from dirs
    integer(kind=ik), intent(in)                    :: neighborhood
    !> difference between block levels
    integer(kind=ik), intent(in)                    :: level_diff

    ! grid parameter
    integer(kind=ik)                                :: Bs, g
    ! loop variables
    integer(kind=ik)                                :: dF

    ! interpolation variables
    real(kind=rk), dimension(:,:,:), allocatable    :: data_corner, data_corner_fine, data_face, data_face_fine


    ! variable for non-uniform mesh correction: remove redundant node between fine->coarse blocks
    integer(kind=ik)                                :: rmv_redundant

    ! lower and upper block nodes bounds
    integer(kind=ik), dimension(2)                  :: sender_x, sender_y, sender_z, receiver_x, receiver_y, receiver_z, finer_x, finer_y, finer_z

!---------------------------------------------------------------------------------------------
! interfaces

    ! grid parameter
    Bs    = params%number_block_nodes
    g     = params%number_ghost_nodes

    rmv_redundant = 0

    ! set non-uniform mesh correction
    if ( params%non_uniform_mesh_correction ) then
        rmv_redundant = 1
    else
        rmv_redundant = 0
    end if

!---------------------------------------------------------------------------------------------
! variables initialization

    allocate( data_corner( g, g, g)  )
    allocate( data_corner_fine( 2*g-1, 2*g-1, 2*g-1)  )

    !allocate( data_face( (Bs+1)/2 + g/2, (Bs+1)/2 + g/2, (Bs+1)/2 + g/2)  )
    allocate( data_face( (Bs+1)/2 + g, (Bs+1)/2 + g, (Bs+1)/2 + g)  )
    !allocate( data_face_fine( Bs+g, Bs+g, Bs+g)  )
    allocate( data_face_fine( Bs+2*g, Bs+2*g, Bs+2*g)  )

!---------------------------------------------------------------------------------------------
! main body

    select case(neighborhood)
        ! '__1/___'
        case(1)
            if ( level_diff == 0 ) then
                ! sender/receiver on same level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
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
                do dF = 1, params%number_data_fields
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
                do dF = 1, params%number_data_fields
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
                do dF = 1, params%number_data_fields
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
                do dF = 1, params%number_data_fields
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
                do dF = 1, params%number_data_fields
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
                do dF = 1, params%number_data_fields
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
                do dF = 1, params%number_data_fields
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
                do dF = 1, params%number_data_fields
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
                do dF = 1, params%number_data_fields
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
                do dF = 1, params%number_data_fields
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
                do dF = 1, params%number_data_fields
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
                do dF = 1, params%number_data_fields
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
                do dF = 1, params%number_data_fields
                    hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+2:g+1+g, g+1:Bs+g, g+2:g+1+g, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

!!        ! '_62/___'
!!        case(11)
!!            if ( level_diff == 0 ) then
!!                ! blocks on same level
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    hvy_block( g+1:Bs+g, Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g, g+2:g+1+g, g+2:g+1+g, dF, sender_id )
!!                end do
!!
!!            else
!!                ! error case
!!                write(*,'(80("_"))')
!!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!!                stop
!!            end if
!
!        ! '_63/___'
!        case(12)
!            if ( level_diff == 0 ) then
!                ! blocks on same level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    hvy_block( 1:g, g+1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( Bs:Bs-1+g, g+1:Bs+g, g+2:g+1+g, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if
!
!!        ! '_64/___'
!!        case(13)
!!            if ( level_diff == 0 ) then
!!                ! blocks on same level
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    hvy_block( g+1:Bs+g, 1:g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g, Bs:Bs-1+g, g+2:g+1+g, dF, sender_id )
!!                end do
!!
!!            else
!!                ! error case
!!                write(*,'(80("_"))')
!!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!!                stop
!!            end if
!
!        ! '_65/___'
!        case(14)
!            if ( level_diff == 0 ) then
!                ! blocks on same level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+2:g+1+g, g+1:Bs+g, g+2:g+1+g, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if
!
        ! '_23/___'
        case(15)
            if ( level_diff == 0 ) then
                ! blocks on same level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
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
                do dF = 1, params%number_data_fields
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
                do dF = 1, params%number_data_fields
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
                do dF = 1, params%number_data_fields
                    hvy_block( Bs+g+1:Bs+g+g, 1:g, g+1:Bs+g, dF, receiver_id ) = hvy_block( g+2:g+1+g, Bs:Bs-1+g, g+1:Bs+g, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        case(19,20,21,22)
            if ( level_diff == 0 ) then
                ! blocks on same level
                ! z
                sender_z(1)   = Bs
                sender_z(2)   = Bs-1+g
                receiver_z(1) = 1
                receiver_z(2) = g
                ! x,y
                select case(neighborhood)
                    case(19) ! '123/___'
                        sender_x(1)   = Bs
                        sender_x(2)   = Bs-1+g
                        receiver_x(1) = 1
                        receiver_x(2) = g

                        sender_y(1)   = g+2
                        sender_y(2)   = g+1+g
                        receiver_y(1) = Bs+g+1
                        receiver_y(2) = Bs+g+g

                    case(20) ! '134/___'
                        sender_x(1)   = Bs
                        sender_x(2)   = Bs-1+g
                        receiver_x(1) = 1
                        receiver_x(2) = g

                        sender_y(1)   = Bs
                        sender_y(2)   = Bs-1+g
                        receiver_y(1) = 1
                        receiver_y(2) = g

                    case(21) ! '145/___'
                        sender_x(1)   = g+2
                        sender_x(2)   = g+1+g
                        receiver_x(1) = Bs+g+1
                        receiver_x(2) = Bs+g+g

                        sender_y(1)   = Bs
                        sender_y(2)   = Bs-1+g
                        receiver_y(1) = 1
                        receiver_y(2) = g

                    case(22) ! '152/___'
                        sender_x(1)   = g+2
                        sender_x(2)   = g+1+g
                        receiver_x(1) = Bs+g+1
                        receiver_x(2) = Bs+g+g

                        sender_y(1)   = g+2
                        sender_y(2)   = g+1+g
                        receiver_y(1) = Bs+g+1
                        receiver_y(2) = Bs+g+g

                end select

                ! copy data
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! copy data
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = hvy_block( sender_x(1):sender_x(2), sender_y(1):sender_y(2), sender_z(1):sender_z(2), dF, sender_id )
                end do

            elseif ( level_diff == -1 ) then
                ! z
                sender_z(1)   = Bs+1
                sender_z(2)   = Bs+g
                receiver_z(1) = 1
                receiver_z(2) = g
                finer_z(1)    = g-1
                finer_z(2)    = 2*g-2
                ! x,y
                select case(neighborhood)
                    case(19) ! '123/___'
                        sender_x(1)   = Bs+1
                        sender_x(2)   = Bs+g
                        receiver_x(1) = 1
                        receiver_x(2) = g
                        finer_x(1)    = g-1
                        finer_x(2)    = 2*g-2

                        sender_y(1)   = g+1
                        sender_y(2)   = g+g
                        receiver_y(1) = Bs+g+1
                        receiver_y(2) = Bs+g+g
                        finer_y(1)    = 2
                        finer_y(2)    = g+1

                    case(20) ! '134/___'
                        sender_x(1)   = Bs+1
                        sender_x(2)   = Bs+g
                        receiver_x(1) = 1
                        receiver_x(2) = g
                        finer_x(1)    = g-1
                        finer_x(2)    = 2*g-2

                        sender_y(1)   = Bs+1
                        sender_y(2)   = Bs+g
                        receiver_y(1) = 1
                        receiver_y(2) = g
                        finer_y(1)    = g-1
                        finer_y(2)    = 2*g-2

                    case(21) ! '145/___'
                        sender_x(1)   = g+1
                        sender_x(2)   = g+g
                        receiver_x(1) = Bs+g+1
                        receiver_x(2) = Bs+g+g
                        finer_x(1)    = 2
                        finer_x(2)    = g+1

                        sender_y(1)   = Bs+1
                        sender_y(2)   = Bs+g
                        receiver_y(1) = 1
                        receiver_y(2) = g
                        finer_y(1)    = g-1
                        finer_y(2)    = 2*g-2

                    case(22) ! '152/___'
                        sender_x(1)   = g+1
                        sender_x(2)   = g+g
                        receiver_x(1) = Bs+g+1
                        receiver_x(2) = Bs+g+g
                        finer_x(1)    = 2
                        finer_x(2)    = g+1

                        sender_y(1)   = g+1
                        sender_y(2)   = g+g
                        receiver_y(1) = Bs+g+1
                        receiver_y(2) = Bs+g+g
                        finer_y(1)    = 2
                        finer_y(2)    = g+1

                end select

                ! copy data
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    data_corner = hvy_block( sender_x(1):sender_x(2), sender_y(1):sender_y(2), sender_z(1):sender_z(2), dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_corner , data_corner_fine, params%order_predictor)
                    ! copy data
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = data_corner_fine(finer_x(1):finer_x(2), finer_y(1):finer_y(2), finer_z(1):finer_z(2))
                end do

            elseif ( level_diff == 1 ) then
                ! x,y,z
                sender_z(1)   = Bs-g
                sender_z(2)   = Bs+g-2+2*rmv_redundant
                receiver_z(1) = 1
                receiver_z(2) = g+rmv_redundant
                ! x,y
                select case(neighborhood)
                    case(19) ! '123/___'
                        sender_x(1)   = Bs-g
                        sender_x(2)   = Bs+g-2+2*rmv_redundant
                        receiver_x(1) = 1
                        receiver_x(2) = g+rmv_redundant

                        sender_y(1)   = g+3-2*rmv_redundant
                        sender_y(2)   = g+1+g+g
                        receiver_y(1) = Bs+g+1-rmv_redundant
                        receiver_y(2) = Bs+g+g

                    case(20) ! '134/___'
                        sender_x(1)   = Bs-g
                        sender_x(2)   = Bs+g-2+2*rmv_redundant
                        receiver_x(1) = 1
                        receiver_x(2) = g+rmv_redundant

                        sender_y(1)   = Bs-g
                        sender_y(2)   = Bs+g-2+2*rmv_redundant
                        receiver_y(1) = 1
                        receiver_y(2) = g+rmv_redundant

                    case(21) ! '145/___'
                        sender_x(1)   = g+3-2*rmv_redundant
                        sender_x(2)   = g+1+g+g
                        receiver_x(1) = Bs+g+1-rmv_redundant
                        receiver_x(2) = Bs+g+g

                        sender_y(1)   = Bs-g
                        sender_y(2)   = Bs+g-2+2*rmv_redundant
                        receiver_y(1) = 1
                        receiver_y(2) = g+rmv_redundant

                    case(22) ! '152/___'
                        sender_x(1)   = g+3-2*rmv_redundant
                        sender_x(2)   = g+1+g+g
                        receiver_x(1) = Bs+g+1-rmv_redundant
                        receiver_x(2) = Bs+g+g

                        sender_y(1)   = g+3-2*rmv_redundant
                        sender_y(2)   = g+1+g+g
                        receiver_y(1) = Bs+g+1-rmv_redundant
                        receiver_y(2) = Bs+g+g

                end select

                ! copy data
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = hvy_block( sender_x(1):sender_x(2):2, sender_y(1):sender_y(2):2, sender_z(1):sender_z(2):2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

!!        ! '123/___'
!!        case(19)
!!            if ( level_diff == 0 ) then
!!                ! blocks on same level
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    hvy_block( 1:g, Bs+g+1:Bs+g+g, 1:g, dF, receiver_id ) = hvy_block( Bs:Bs-1+g, g+2:g+1+g, Bs:Bs-1+g, dF, sender_id )
!!                end do
!!
!!            elseif ( level_diff == -1 ) then
!!                ! sender one level down
!!                ! interpolate data
!!                do dF = 1, params%number_data_fields
!!                    ! data to refine
!!                    data_corner = hvy_block( Bs+1:Bs+g, g+1:g+g, Bs+1:Bs+g, dF, sender_id )
!!                    ! interpolate data
!!                    call prediction_3D( data_corner , data_corner_fine, params%order_predictor)
!!                    ! data to synchronize
!!                    data_corner = data_corner_fine(g-1:2*g-2, 2:g+1, g-1:2*g-2)
!!                    ! write data
!!                    hvy_block( 1:g, Bs+g+1:Bs+g+g, 1:g, dF, receiver_id ) = data_corner
!!                end do
!!
!!            elseif ( level_diff == 1) then
!!                ! sender one level up
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    !hvy_block( 1:g, Bs+g+1:Bs+g+g, 1:g, dF, receiver_id ) = hvy_block( Bs-g:Bs-2+g:2, g+3:g+1+g+g:2, Bs-g:Bs-2+g:2, dF, sender_id )
!!                    hvy_block( 1:g+rmv_redundant, Bs+g+1-rmv_redundant:Bs+g+g, 1:g+rmv_redundant, dF, receiver_id ) = hvy_block( Bs-g:Bs-2+g+2*rmv_redundant:2, g+3-2*rmv_redundant:g+1+g+g:2, Bs-g:Bs-2+g+2*rmv_redundant:2, dF, sender_id )
!!                end do
!!
!!            else
!!                ! error case
!!                write(*,'(80("_"))')
!!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!!                stop
!!            end if
!        ! '134/___'
!        case(20)
!            if ( level_diff == 0 ) then
!                ! blocks on same level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    hvy_block( 1:g, 1:g, 1:g, dF, receiver_id ) = 17.0_rk!hvy_block( Bs:Bs-1+g, Bs:Bs-1+g, Bs:Bs-1+g, dF, sender_id )
!                end do
!
!            elseif ( level_diff == -1 ) then
!                ! sender one level down
!                ! interpolate data
!                do dF = 1, params%number_data_fields
!                    ! data to refine
!                    data_corner = hvy_block( Bs+1:Bs+g, Bs+1:Bs+g, Bs+1:Bs+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_corner , data_corner_fine, params%order_predictor)
!                    ! data to synchronize
!                    data_corner = data_corner_fine(g-1:2*g-2, g-1:2*g-2, g-1:2*g-2)
!                    ! write data
!                    hvy_block( 1:g, 1:g, 1:g, dF, receiver_id ) = 17.0_rk!data_corner
!                end do
!
!            elseif ( level_diff == 1) then
!                ! sender one level up
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( 1:g, 1:g, 1:g, dF, receiver_id ) = hvy_block( Bs-g:Bs-2+g:2, Bs-g:Bs-2+g:2, Bs-g:Bs-2+g:2, dF, sender_id )
!                    hvy_block( 1:g+rmv_redundant, 1:g+rmv_redundant, 1:g+rmv_redundant, dF, receiver_id ) = 17.0_rk!hvy_block( Bs-g:Bs-2+g+2*rmv_redundant:2, Bs-g:Bs-2+g+2*rmv_redundant:2, Bs-g:Bs-2+g+2*rmv_redundant:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if
!!
!!        ! '145/___'
!!        case(21)
!!            if ( level_diff == 0 ) then
!!                ! blocks on same level
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    hvy_block( Bs+g+1:Bs+g+g, 1:g, 1:g, dF, receiver_id ) = hvy_block( g+2:g+1+g, Bs:Bs-1+g, Bs:Bs-1+g, dF, sender_id )
!!                end do
!!
!!            elseif ( level_diff == -1 ) then
!!                ! sender one level down
!!                ! interpolate data
!!                do dF = 1, params%number_data_fields
!!                    ! data to refine
!!                    data_corner = hvy_block( g+1:g+g, Bs+1:Bs+g, Bs+1:Bs+g, dF, sender_id )
!!                    ! interpolate data
!!                    call prediction_3D( data_corner , data_corner_fine, params%order_predictor)
!!                    ! data to synchronize
!!                    data_corner = data_corner_fine(2:g+1, g-1:2*g-2, g-1:2*g-2)
!!                    ! write data
!!                    hvy_block( Bs+g+1:Bs+g+g, 1:g, 1:g, dF, receiver_id ) = data_corner
!!                end do
!!
!!            elseif ( level_diff == 1) then
!!                ! sender one level up
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    !hvy_block( Bs+g+1:Bs+g+g, 1:g, 1:g, dF, receiver_id ) = hvy_block( g+3:g+1+g+g:2, Bs-g:Bs-2+g:2, Bs-g:Bs-2+g:2, dF, sender_id )
!!                    hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, 1:g+rmv_redundant, 1:g+rmv_redundant, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:g+1+g+g:2, Bs-g:Bs-2+g+2*rmv_redundant:2, Bs-g:Bs-2+g+2*rmv_redundant:2, dF, sender_id )
!!                end do
!!
!!            else
!!                ! error case
!!                write(*,'(80("_"))')
!!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!!                stop
!!            end if
!!        ! '152/___'
!!        case(22)
!!            if ( level_diff == 0 ) then
!!                ! blocks on same level
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    hvy_block( Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, 1:g, dF, receiver_id ) = hvy_block( g+2:g+1+g, g+2:g+1+g, Bs:Bs-1+g, dF, sender_id )
!!                end do
!!
!!            elseif ( level_diff == -1 ) then
!!                ! sender one level down
!!                ! interpolate data
!!                do dF = 1, params%number_data_fields
!!                    ! data to refine
!!                    data_corner = hvy_block( g+1:g+g, g+1:g+g, Bs+1:Bs+g, dF, sender_id )
!!                    ! interpolate data
!!                    call prediction_3D( data_corner , data_corner_fine, params%order_predictor)
!!                    ! data to synchronize
!!                    data_corner = data_corner_fine(2:g+1, 2:g+1, g-1:2*g-2)
!!                    ! write data
!!                    hvy_block( Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, 1:g, dF, receiver_id ) = data_corner
!!                end do
!!
!!            elseif ( level_diff == 1) then
!!                ! sender one level up
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    !hvy_block( Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, 1:g, dF, receiver_id ) = hvy_block( g+3:g+1+g+g:2, g+3:g+1+g+g:2, Bs-g:Bs-2+g:2, dF, sender_id )
!!                    hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, Bs+g+1-rmv_redundant:Bs+g+g, 1:g+rmv_redundant, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:g+1+g+g:2, g+3-2*rmv_redundant:g+1+g+g:2, Bs-g:Bs-2+g+2*rmv_redundant:2, dF, sender_id )
!!                end do
!!
!!            else
!!                ! error case
!!                write(*,'(80("_"))')
!!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!!                stop
!!            end if

          case(23,24,25,26)
            if ( level_diff == 0 ) then
                ! blocks on same level
                ! z
                sender_z(1)   = g+2
                sender_z(2)   = g+1+g
                receiver_z(1) = Bs+g+1
                receiver_z(2) = Bs+g+g
                ! x,y
                select case(neighborhood)
                    case(23) ! '623/___'
                        sender_x(1)   = Bs
                        sender_x(2)   = Bs-1+g
                        receiver_x(1) = 1
                        receiver_x(2) = g

                        sender_y(1)   = g+2
                        sender_y(2)   = g+1+g
                        receiver_y(1) = Bs+g+1
                        receiver_y(2) = Bs+g+g

                    case(24) ! '634/___'
                        sender_x(1)   = Bs
                        sender_x(2)   = Bs-1+g
                        receiver_x(1) = 1
                        receiver_x(2) = g

                        sender_y(1)   = Bs
                        sender_y(2)   = Bs-1+g
                        receiver_y(1) = 1
                        receiver_y(2) = g

                    case(25) ! '645/___'
                        sender_x(1)   = g+2
                        sender_x(2)   = g+1+g
                        receiver_x(1) = Bs+g+1
                        receiver_x(2) = Bs+g+g

                        sender_y(1)   = Bs
                        sender_y(2)   = Bs-1+g
                        receiver_y(1) = 1
                        receiver_y(2) = g

                    case(26) ! '652/___'
                        sender_x(1)   = g+2
                        sender_x(2)   = g+1+g
                        receiver_x(1) = Bs+g+1
                        receiver_x(2) = Bs+g+g

                        sender_y(1)   = g+2
                        sender_y(2)   = g+1+g
                        receiver_y(1) = Bs+g+1
                        receiver_y(2) = Bs+g+g

                end select

                ! copy data
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! copy data
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = hvy_block( sender_x(1):sender_x(2), sender_y(1):sender_y(2), sender_z(1):sender_z(2), dF, sender_id )
                end do

            elseif ( level_diff == -1 ) then
                ! z
                sender_z(1)   = g+1
                sender_z(2)   = g+g
                receiver_z(1) = Bs+g+1
                receiver_z(2) = Bs+g+g
                finer_z(1)    = 2
                finer_z(2)    = g+1
                ! x,y
                select case(neighborhood)
                    case(23) ! '623/___'
                        sender_x(1)   = Bs+1
                        sender_x(2)   = Bs+g
                        receiver_x(1) = 1
                        receiver_x(2) = g
                        finer_x(1)    = g-1
                        finer_x(2)    = 2*g-2

                        sender_y(1)   = g+1
                        sender_y(2)   = g+g
                        receiver_y(1) = Bs+g+1
                        receiver_y(2) = Bs+g+g
                        finer_y(1)    = 2
                        finer_y(2)    = g+1

                    case(24) ! '634/___'
                        sender_x(1)   = Bs+1
                        sender_x(2)   = Bs+g
                        receiver_x(1) = 1
                        receiver_x(2) = g
                        finer_x(1)    = g-1
                        finer_x(2)    = 2*g-2

                        sender_y(1)   = Bs+1
                        sender_y(2)   = Bs+g
                        receiver_y(1) = 1
                        receiver_y(2) = g
                        finer_y(1)    = g-1
                        finer_y(2)    = 2*g-2

                    case(25) ! '645/___'
                        sender_x(1)   = g+1
                        sender_x(2)   = g+g
                        receiver_x(1) = Bs+g+1
                        receiver_x(2) = Bs+g+g
                        finer_x(1)    = 2
                        finer_x(2)    = g+1

                        sender_y(1)   = Bs+1
                        sender_y(2)   = Bs+g
                        receiver_y(1) = 1
                        receiver_y(2) = g
                        finer_y(1)    = g-1
                        finer_y(2)    = 2*g-2

                    case(26) ! '652/___'
                        sender_x(1)   = g+1
                        sender_x(2)   = g+g
                        receiver_x(1) = Bs+g+1
                        receiver_x(2) = Bs+g+g
                        finer_x(1)    = 2
                        finer_x(2)    = g+1

                        sender_y(1)   = g+1
                        sender_y(2)   = g+g
                        receiver_y(1) = Bs+g+1
                        receiver_y(2) = Bs+g+g
                        finer_y(1)    = 2
                        finer_y(2)    = g+1

                end select

                ! copy data
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    data_corner = hvy_block( sender_x(1):sender_x(2), sender_y(1):sender_y(2), sender_z(1):sender_z(2), dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_corner , data_corner_fine, params%order_predictor)
                    ! copy data
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = data_corner_fine(finer_x(1):finer_x(2), finer_y(1):finer_y(2), finer_z(1):finer_z(2))
                end do

            elseif ( level_diff == 1 ) then
                ! x,y,z
                sender_z(1)   = g+3-2*rmv_redundant
                sender_z(2)   = g+1+g+g
                receiver_z(1) = Bs+g+1-rmv_redundant
                receiver_z(2) = Bs+g+g
                ! x,y
                select case(neighborhood)
                    case(23) ! '623/___'
                        sender_x(1)   = Bs-g
                        sender_x(2)   = Bs+g-2+2*rmv_redundant
                        receiver_x(1) = 1
                        receiver_x(2) = g+rmv_redundant

                        sender_y(1)   = g+3-2*rmv_redundant
                        sender_y(2)   = g+1+g+g
                        receiver_y(1) = Bs+g+1-rmv_redundant
                        receiver_y(2) = Bs+g+g

                   case(24) ! '634/___'
                        sender_x(1)   = Bs-g
                        sender_x(2)   = Bs+g-2+2*rmv_redundant
                        receiver_x(1) = 1
                        receiver_x(2) = g+rmv_redundant

                        sender_y(1)   = Bs-g
                        sender_y(2)   = Bs+g-2+2*rmv_redundant
                        receiver_y(1) = 1
                        receiver_y(2) = g+rmv_redundant

                    case(25) ! '645/___'
                        sender_x(1)   = g+3-2*rmv_redundant
                        sender_x(2)   = g+1+g+g
                        receiver_x(1) = Bs+g+1-rmv_redundant
                        receiver_x(2) = Bs+g+g

                        sender_y(1)   = Bs-g
                        sender_y(2)   = Bs+g-2+2*rmv_redundant
                        receiver_y(1) = 1
                        receiver_y(2) = g+rmv_redundant

                    case(26) ! '652/___'
                        sender_x(1)   = g+3-2*rmv_redundant
                        sender_x(2)   = g+1+g+g
                        receiver_x(1) = Bs+g+1-rmv_redundant
                        receiver_x(2) = Bs+g+g

                        sender_y(1)   = g+3-2*rmv_redundant
                        sender_y(2)   = g+1+g+g
                        receiver_y(1) = Bs+g+1-rmv_redundant
                        receiver_y(2) = Bs+g+g

                end select

                ! copy data
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = hvy_block( sender_x(1):sender_x(2):2, sender_y(1):sender_y(2):2, sender_z(1):sender_z(2):2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

!!        ! '623/___'
!!        case(23)
!!            if ( level_diff == 0 ) then
!!                ! blocks on same level
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    hvy_block( 1:g, Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( Bs:Bs-1+g, g+2:g+1+g, g+2:g+1+g, dF, sender_id )
!!                end do
!!
!!            elseif ( level_diff == -1 ) then
!!                ! sender one level down
!!                ! interpolate data
!!                do dF = 1, params%number_data_fields
!!                    ! data to refine
!!                    data_corner = hvy_block( Bs+1:Bs+g, g+1:g+g, g+1:g+g, dF, sender_id )
!!                    ! interpolate data
!!                    call prediction_3D( data_corner , data_corner_fine, params%order_predictor)
!!                    ! data to synchronize
!!                    data_corner = data_corner_fine(g-1:2*g-2, 2:g+1, 2:g+1)
!!                    ! write data
!!                    hvy_block( 1:g, Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_corner
!!                end do
!!
!!            elseif ( level_diff == 1) then
!!                ! sender one level up
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    !hvy_block( 1:g, Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( Bs-g:Bs-2+g:2, g+3:g+1+g+g:2, g+3:g+1+g+g:2, dF, sender_id )
!!                    hvy_block( 1:g+rmv_redundant, Bs+g+1-rmv_redundant:Bs+g+g, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( Bs-g:Bs-2+g+2*rmv_redundant:2, g+3-2*rmv_redundant:g+1+g+g:2, g+3-2*rmv_redundant:g+1+g+g:2, dF, sender_id )
!!                end do
!!
!!            else
!!                ! error case
!!                write(*,'(80("_"))')
!!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!!                stop
!!            end if
!!
!!        ! '634/___'
!!        case(24)
!!            if ( level_diff == 0 ) then
!!                ! blocks on same level
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    hvy_block( 1:g, 1:g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( Bs:Bs-1+g, Bs:Bs-1+g, g+2:g+1+g, dF, sender_id )
!!                end do
!!
!!            elseif ( level_diff == -1 ) then
!!                ! sender one level down
!!                ! interpolate data
!!                do dF = 1, params%number_data_fields
!!                    ! data to refine
!!                    data_corner = hvy_block( Bs+1:Bs+g, Bs+1:Bs+g, g+1:g+g, dF, sender_id )
!!                    ! interpolate data
!!                    call prediction_3D( data_corner , data_corner_fine, params%order_predictor)
!!                    ! data to synchronize
!!                    data_corner = data_corner_fine(g-1:2*g-2, g-1:2*g-2, 2:g+1)
!!                    ! write data
!!                    hvy_block( 1:g, 1:g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_corner
!!                end do
!!
!!            elseif ( level_diff == 1) then
!!                ! sender one level up
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    !hvy_block( 1:g, 1:g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( Bs-g:Bs-2+g:2, Bs-g:Bs-2+g:2, g+3:g+1+g+g:2, dF, sender_id )
!!                    hvy_block( 1:g+rmv_redundant, 1:g+rmv_redundant, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( Bs-g:Bs-2+g+2*rmv_redundant:2, Bs-g:Bs-2+g+2*rmv_redundant:2, g+3-2*rmv_redundant:g+1+g+g:2, dF, sender_id )
!!                end do
!!
!!            else
!!                ! error case
!!                write(*,'(80("_"))')
!!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!!                stop
!!            end if
!!
!!        ! '645/___'
!!        case(25)
!!            if ( level_diff == 0 ) then
!!                ! blocks on same level
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    hvy_block( Bs+g+1:Bs+g+g, 1:g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+2:g+1+g, Bs:Bs-1+g, g+2:g+1+g, dF, sender_id )
!!                end do
!!
!!            elseif ( level_diff == -1 ) then
!!                ! sender one level down
!!                ! interpolate data
!!                do dF = 1, params%number_data_fields
!!                    ! data to refine
!!                    data_corner = hvy_block( g+1:g+g, Bs+1:Bs+g, g+1:g+g, dF, sender_id )
!!                    ! interpolate data
!!                    call prediction_3D( data_corner , data_corner_fine, params%order_predictor)
!!                    ! data to synchronize
!!                    data_corner = data_corner_fine(2:g+1, g-1:2*g-2, 2:g+1)
!!                    ! write data
!!                    hvy_block( Bs+g+1:Bs+g+g, 1:g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_corner
!!                end do
!!
!!            elseif ( level_diff == 1) then
!!                ! sender one level up
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    !hvy_block( Bs+g+1:Bs+g+g, 1:g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+3:g+1+g+g:2, Bs-g:Bs-2+g:2, g+3:g+1+g+g:2, dF, sender_id )
!!                    hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, 1:g+rmv_redundant, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:g+1+g+g:2, Bs-g:Bs-2+g+2*rmv_redundant:2, g+3-2*rmv_redundant:g+1+g+g:2, dF, sender_id )
!!                end do
!!
!!            else
!!                ! error case
!!                write(*,'(80("_"))')
!!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!!                stop
!!            end if
!!
!!        ! '652/___'
!!        case(26)
!!            if ( level_diff == 0 ) then
!!                ! blocks on same level
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    hvy_block( Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+2:g+1+g, g+2:g+1+g, g+2:g+1+g, dF, sender_id )
!!                end do
!!
!!            elseif ( level_diff == -1 ) then
!!                ! sender one level down
!!                ! interpolate data
!!                do dF = 1, params%number_data_fields
!!                    ! data to refine
!!                    data_corner = hvy_block( g+1:g+g, g+1:g+g, g+1:g+g, dF, sender_id )
!!                    ! interpolate data
!!                    call prediction_3D( data_corner , data_corner_fine, params%order_predictor)
!!                    ! data to synchronize
!!                    data_corner = data_corner_fine(2:g+1, 2:g+1, 2:g+1)
!!                    ! write data
!!                    hvy_block( Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g,Bs+g+1:Bs+g+g, dF, receiver_id ) = data_corner
!!                end do
!!
!!            elseif ( level_diff == 1) then
!!                ! sender one level up
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    !hvy_block( Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+3:g+1+g+g:2, g+3:g+1+g+g:2, g+3:g+1+g+g:2, dF, sender_id )
!!                    hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, Bs+g+1-rmv_redundant:Bs+g+g, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:g+1+g+g:2, g+3-2*rmv_redundant:g+1+g+g:2, g+3-2*rmv_redundant:g+1+g+g:2, dF, sender_id )
!!                end do
!!
!!            else
!!                ! error case
!!                write(*,'(80("_"))')
!!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!!                stop
!!            end if

        case(27,28,29,30)
            if ( level_diff == -1 ) then
                ! z
                sender_z(1)   = (Bs+1)/2
                sender_z(2)   = Bs+g
                receiver_z(1) = 1
                receiver_z(2) = g
                finer_z(1)    = Bs+g
                finer_z(2)    = Bs-1+2*g
                ! x,y
                select case(neighborhood)
                    case(27) ! '__1/123'
                        sender_x(1)   = (Bs+1)/2
                        sender_x(2)   = Bs+g
                        receiver_x(1) = 1
                        receiver_x(2) = Bs+g
                        finer_x(1)    = g+1
                        finer_x(2)    = Bs+2*g

                        sender_y(1)   = g+1
                        sender_y(2)   = (Bs+1)/2+g+g
                        receiver_y(1) = g+1
                        receiver_y(2) = Bs+2*g
                        finer_y(1)    = 1
                        finer_y(2)    = Bs+g

                    case(28) ! '__1/134'
                        sender_x(1)   = (Bs+1)/2
                        sender_x(2)   = Bs+g
                        receiver_x(1) = 1
                        receiver_x(2) = Bs+g
                        finer_x(1)    = g+1
                        finer_x(2)    = Bs+2*g

                        sender_y(1)   = (Bs+1)/2
                        sender_y(2)   = Bs+g
                        receiver_y(1) = 1
                        receiver_y(2) = Bs+g
                        finer_y(1)    = g+1
                        finer_y(2)    = Bs+2*g

                    case(29) ! '__1/145'
                        sender_x(1)   = g+1
                        sender_x(2)   = (Bs+1)/2+g+g
                        receiver_x(1) = g+1
                        receiver_x(2) = Bs+2*g
                        finer_x(1)    = 1
                        finer_x(2)    = Bs+g

                        sender_y(1)   = (Bs+1)/2
                        sender_y(2)   = Bs+g
                        receiver_y(1) = 1
                        receiver_y(2) = Bs+g
                        finer_y(1)    = g+1
                        finer_y(2)    = Bs+2*g

                    case(30) ! '__1/152'
                        sender_x(1)   = g+1
                        sender_x(2)   = (Bs+1)/2+g+g
                        receiver_x(1) = g+1
                        receiver_x(2) = Bs+2*g
                        finer_x(1)    = 1
                        finer_x(2)    = Bs+g

                        sender_y(1)   = g+1
                        sender_y(2)   = (Bs+1)/2+g+g
                        receiver_y(1) = g+1
                        receiver_y(2) = Bs+2*g
                        finer_y(1)    = 1
                        finer_y(2)    = Bs+g

                end select

                ! copy data
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    data_face = hvy_block( sender_x(1):sender_x(2), sender_y(1):sender_y(2), sender_z(1):sender_z(2), dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = data_face_fine(finer_x(1):finer_x(2), finer_y(1):finer_y(2), finer_z(1):finer_z(2))
                end do

            elseif ( level_diff == 1 ) then
                ! x,y,z
                sender_z(1)   = Bs-g
                sender_z(2)   = Bs+g-2+2*rmv_redundant
                sender_y(1)   = g+1
                sender_y(2)   = Bs+g
                sender_x(1)   = g+1
                sender_x(2)   = Bs+g

                receiver_z(1) = 1
                receiver_z(2) = g+rmv_redundant
                ! x,y
                select case(neighborhood)
                    case(27) ! '__1/123'
                        receiver_x(1) = g+(Bs+1)/2
                        receiver_x(2) = Bs+g

                        receiver_y(1) = g+1
                        receiver_y(2) = g+(Bs+1)/2

                    case(28) ! '__1/134'
                        receiver_x(1) = g+(Bs+1)/2
                        receiver_x(2) = Bs+g

                        receiver_y(1) = g+(Bs+1)/2
                        receiver_y(2) = Bs+g

                    case(29) ! '__1/145'
                        receiver_x(1) = g+1
                        receiver_x(2) = g+(Bs+1)/2

                        receiver_y(1) = g+(Bs+1)/2
                        receiver_y(2) = Bs+g

                    case(30) ! '__1/152'
                        receiver_x(1) = g+1
                        receiver_x(2) = g+(Bs+1)/2

                        receiver_y(1) = g+1
                        receiver_y(2) = g+(Bs+1)/2

                end select

                ! copy data
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = hvy_block( sender_x(1):sender_x(2):2, sender_y(1):sender_y(2):2, sender_z(1):sender_z(2):2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

!
!        ! '__1/123'
!        case(27)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate
!                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
!                    !data_face = hvy_block( (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, dF, sender_id )
!                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( g+1:Bs+g, g+1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, g+1:Bs+g, Bs:Bs-1+g)
!                    !hvy_block( 1:Bs+g, g+1:Bs+2*g, 1:g, dF, receiver_id ) = data_face_fine(1:Bs+g, 1:Bs+g, Bs:Bs-1+g)
!                    !hvy_block( 1:Bs+g, g+1:Bs+2*g, 1:g, dF, receiver_id ) = data_face_fine(g+1:Bs+2*g, 1:Bs+g, Bs+g:Bs-1+2*g)
!                    hvy_block( g+1:Bs+2*g, 1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(1:Bs+g, g+1:Bs+2*g, Bs+g:Bs-1+2*g)
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( g+(Bs+1)/2:Bs+g, g+1:g+(Bs+1)/2, 1:g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, Bs-g:Bs+g-2:2, dF, sender_id )
!                    !hvy_block( g+(Bs+1)/2:Bs+g, g+1:g+(Bs+1)/2, 1:g+rmv_redundant, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, dF, sender_id )
!                    hvy_block( g+1:g+(Bs+1)/2, g+(Bs+1)/2:Bs+g, 1:g+rmv_redundant, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if
!        ! '__1/134'
!        case(28)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate
!                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
!                    !data_face = hvy_block( (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, dF, sender_id )
!                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( g+1:Bs+g, g+1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, g+1:Bs+g, Bs:Bs-1+g)
!                    !hvy_block( 1:Bs+g, 1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(1:Bs+g, 1:Bs+g, Bs:Bs-1+g)
!                    !hvy_block( 1:Bs+g, 1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(g+1:Bs+2*g, g+1:Bs+2*g, Bs+g:Bs-1+2*g)
!                    hvy_block( g+1:Bs+2*g, g+1:Bs+2*g, 1:g, dF, receiver_id ) = data_face_fine(1:Bs+g, 1:Bs+g, Bs+g:Bs-1+2*g)
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( g+(Bs+1)/2:Bs+g, g+(Bs+1)/2:Bs+g, 1:g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, Bs-g:Bs+g-2:2, dF, sender_id )
!                    !hvy_block( g+(Bs+1)/2:Bs+g, g+(Bs+1)/2:Bs+g, 1:g+rmv_redundant, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, dF, sender_id )
!                    hvy_block( g+1:g+(Bs+1)/2, g+1:g+(Bs+1)/2, 1:g+rmv_redundant, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if
!        ! '__1/145'
!        case(29)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate
!                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
!                    !data_face = hvy_block( g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, dF, sender_id )
!                    data_face = hvy_block( (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( g+1:Bs+g, g+1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, g+1:Bs+g, Bs:Bs-1+g)
!                    !hvy_block( g+1:Bs+2*g, 1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(1:Bs+g, 1:Bs+g, Bs:Bs-1+g)
!                    !hvy_block( g+1:Bs+2*g, 1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(1:Bs+g, g+1:Bs+2*g, Bs+g:Bs-1+2*g)
!                    hvy_block( 1:Bs+g, g+1:Bs+2*g, 1:g, dF, receiver_id ) = data_face_fine(g+1:Bs+2*g, 1:Bs+g, Bs+g:Bs-1+2*g)
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( g+1:g+(Bs+1)/2, g+(Bs+1)/2:Bs+g, 1:g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, Bs-g:Bs+g-2:2, dF, sender_id )
!                    !hvy_block( g+1:g+(Bs+1)/2, g+(Bs+1)/2:Bs+g, 1:g+rmv_redundant, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, dF, sender_id )
!                    hvy_block( g+(Bs+1)/2:Bs+g, g+1:g+(Bs+1)/2, 1:g+rmv_redundant, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if
!        ! '__1/152'
!        case(30)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate
!                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
!                    !data_face = hvy_block( g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, dF, sender_id )
!                    data_face = hvy_block( (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( g+1:Bs+g, g+1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, g+1:Bs+g, Bs:Bs-1+g)
!                    !hvy_block( g+1:Bs+2*g, g+1:Bs+2*g, 1:g, dF, receiver_id ) = data_face_fine(1:Bs+g, 1:Bs+g, Bs:Bs-1+g)
!                    !hvy_block( g+1:Bs+2*g, g+1:Bs+2*g, 1:g, dF, receiver_id ) = data_face_fine(1:Bs+g, 1:Bs+g, Bs+g:Bs-1+2*g)
!                    hvy_block( 1:Bs+g, 1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(g+1:Bs+2*g, g+1:Bs+2*g, Bs+g:Bs-1+2*g)
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( g+1:g+(Bs+1)/2, g+1:g+(Bs+1)/2, 1:g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, Bs-g:Bs+g-2:2, dF, sender_id )
!                    !hvy_block( g+1:g+(Bs+1)/2, g+1:g+(Bs+1)/2, 1:g+rmv_redundant, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, dF, sender_id )
!                    hvy_block( g+(Bs+1)/2:Bs+g, g+(Bs+1)/2:Bs+g, 1:g+rmv_redundant, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if

        case(31,32,33,34)
            if ( level_diff == -1 ) then
                ! y
                sender_y(1)   = g+1
                sender_y(2)   = (Bs+1)/2+g+g
                receiver_y(1) = Bs+g+1
                receiver_y(2) = Bs+g+g
                finer_y(1)    = 2
                finer_y(2)    = g+1
                ! x,z
                select case(neighborhood)
                    case(32) ! '__2/623'
                        sender_x(1)   = (Bs+1)/2
                        sender_x(2)   = Bs+g
                        receiver_x(1) = 1
                        receiver_x(2) = Bs+g
                        finer_x(1)    = g+1
                        finer_x(2)    = Bs+2*g

                        sender_z(1)   = g+1
                        sender_z(2)   = (Bs+1)/2+g+g
                        receiver_z(1) = g+1
                        receiver_z(2) = Bs+2*g
                        finer_z(1)    = 1
                        finer_z(2)    = Bs+g

                    case(31) ! '__2/123'
                        sender_x(1)   = (Bs+1)/2
                        sender_x(2)   = Bs+g
                        receiver_x(1) = 1
                        receiver_x(2) = Bs+g
                        finer_x(1)    = g+1
                        finer_x(2)    = Bs+2*g

                        sender_z(1)   = (Bs+1)/2
                        sender_z(2)   = Bs+g
                        receiver_z(1) = 1
                        receiver_z(2) = Bs+g
                        finer_z(1)    = g+1
                        finer_z(2)    = Bs+2*g

                    case(33) ! '__2/152'
                        sender_x(1)   = g+1
                        sender_x(2)   = (Bs+1)/2+g+g
                        receiver_x(1) = g+1
                        receiver_x(2) = Bs+2*g
                        finer_x(1)    = 1
                        finer_x(2)    = Bs+g

                        sender_z(1)   = (Bs+1)/2
                        sender_z(2)   = Bs+g
                        receiver_z(1) = 1
                        receiver_z(2) = Bs+g
                        finer_z(1)    = g+1
                        finer_z(2)    = Bs+2*g

                    case(34) ! '__2/652'
                        sender_x(1)   = g+1
                        sender_x(2)   = (Bs+1)/2+g+g
                        receiver_x(1) = g+1
                        receiver_x(2) = Bs+2*g
                        finer_x(1)    = 1
                        finer_x(2)    = Bs+g

                        sender_z(1)   = g+1
                        sender_z(2)   = (Bs+1)/2+g+g
                        receiver_z(1) = g+1
                        receiver_z(2) = Bs+2*g
                        finer_z(1)    = 1
                        finer_z(2)    = Bs+g

                end select

                ! copy data
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    data_face = hvy_block( sender_x(1):sender_x(2), sender_y(1):sender_y(2), sender_z(1):sender_z(2), dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = data_face_fine(finer_x(1):finer_x(2), finer_y(1):finer_y(2), finer_z(1):finer_z(2))
                end do

            elseif ( level_diff == 1 ) then
                ! x,y,z
                sender_y(1)   = g+3-2*rmv_redundant
                sender_y(2)   = 3*g+1
                sender_z(1)   = g+1
                sender_z(2)   = Bs+g
                sender_x(1)   = g+1
                sender_x(2)   = Bs+g

                receiver_y(1) = Bs+g+1-rmv_redundant
                receiver_y(2) = Bs+g+g
                ! x,z
                select case(neighborhood)
                    case(32) ! '___2/623'
                        receiver_x(1) = g+(Bs+1)/2
                        receiver_x(2) = Bs+g

                        receiver_z(1) = g+1
                        receiver_z(2) = g+(Bs+1)/2

                    case(31) ! '__2/123'
                        receiver_x(1) = g+(Bs+1)/2
                        receiver_x(2) = Bs+g

                        receiver_z(1) = g+(Bs+1)/2
                        receiver_z(2) = Bs+g

                    case(33) ! '__2/152'
                        receiver_x(1) = g+1
                        receiver_x(2) = g+(Bs+1)/2

                        receiver_z(1) = g+(Bs+1)/2
                        receiver_z(2) = Bs+g

                    case(34) ! '__2/652'
                        receiver_x(1) = g+1
                        receiver_x(2) = g+(Bs+1)/2

                        receiver_z(1) = g+1
                        receiver_z(2) = g+(Bs+1)/2

                end select

                ! copy data
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = hvy_block( sender_x(1):sender_x(2):2, sender_y(1):sender_y(2):2, sender_z(1):sender_z(2):2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

!
!        ! '__2/123'
!        case(31)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate
!                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
!                    !data_face = hvy_block( (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, dF, sender_id )
!                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( g+1:Bs+g, Bs+g+1:Bs+g+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, 2:g+1, g+1:Bs+g)
!                    !hvy_block( 1:Bs+g, Bs+g+1:Bs+g+g, 1:Bs+g, dF, receiver_id ) = data_face_fine(1:Bs+g, 2:g+1, 1:Bs+g)
!                    !hvy_block( 1:Bs+g, Bs+g+1:Bs+g+g, 1:Bs+g, dF, receiver_id ) = data_face_fine(g+1:Bs+2*g, 2:g+1, g+1:Bs+2*g)
!                    hvy_block( g+1:Bs+2*g, Bs+g+1:Bs+g+g, g+1:Bs+2*g, dF, receiver_id ) = data_face_fine(1:Bs+g, 2:g+1, 1:Bs+g)
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( g+(Bs+1)/2:Bs+g, Bs+g+1:Bs+g+g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
!                    !hvy_block( g+(Bs+1)/2:Bs+g, Bs+g+1-rmv_redundant:Bs+g+g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
!                    hvy_block( g+1:g+(Bs+1)/2, Bs+g+1-rmv_redundant:Bs+g+g, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if
!
!        ! '__2/623'
!        case(32)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate
!                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
!                    !data_face = hvy_block( (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( g+1:Bs+g, Bs+g+1:Bs+g+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, 2:g+1, g+1:Bs+g)
!                    !hvy_block( 1:Bs+g, Bs+g+1:Bs+g+g, g+1:Bs+2*g, dF, receiver_id ) = data_face_fine(1:Bs+g, 2:g+1, 1:Bs+g)
!                    hvy_block( g+1:Bs+2*g, Bs+g+1:Bs+g+g, 1:Bs+g, dF, receiver_id ) = data_face_fine(1:Bs+g, 2:g+1, g+1:Bs+2*g)
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( g+(Bs+1)/2:Bs+g, Bs+g+1:Bs+g+g, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
!                    !hvy_block( g+(Bs+1)/2:Bs+g, Bs+g+1-rmv_redundant:Bs+g+g, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
!                    hvy_block( g+1:g+(Bs+1)/2, Bs+g+1-rmv_redundant:Bs+g+g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if
!
!        ! '__2/152'
!        case(33)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate
!                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
!                    !data_face = hvy_block( g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, dF, sender_id )
!                    data_face = hvy_block( (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( g+1:Bs+g, Bs+g+1:Bs+g+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, 2:g+1, g+1:Bs+g)
!                    !hvy_block( g+1:Bs+2*g, Bs+g+1:Bs+g+g, 1:Bs+g, dF, receiver_id ) = data_face_fine(1:Bs+g, 2:g+1, 1:Bs+g)
!                    hvy_block( 1:Bs+g, Bs+g+1:Bs+g+g, g+1:Bs+2*g, dF, receiver_id ) = data_face_fine(g+1:Bs+2*g, 2:g+1, 1:Bs+g)
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( g+1:g+(Bs+1)/2, Bs+g+1:Bs+g+g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
!                    !hvy_block( g+1:g+(Bs+1)/2, Bs+g+1-rmv_redundant:Bs+g+g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
!                    hvy_block( g+(Bs+1)/2:Bs+g, Bs+g+1-rmv_redundant:Bs+g+g, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if
!
!        ! '__2/652'
!        case(34)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate
!                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
!                    !data_face = hvy_block( g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!                    data_face = hvy_block( (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( g+1:Bs+g, Bs+g+1:Bs+g+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, 2:g+1, g+1:Bs+g)
!                    !hvy_block( g+1:Bs+2*g, Bs+g+1:Bs+g+g, g+1:Bs+2*g, dF, receiver_id ) = data_face_fine(1:Bs+g, 2:g+1, 1:Bs+g)
!                    hvy_block( 1:Bs+g, Bs+g+1:Bs+g+g, 1:Bs+g, dF, receiver_id ) = data_face_fine(g+1:Bs+2*g, 2:g+1, g+1:Bs+2*g)
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( g+1:g+(Bs+1)/2, Bs+g+1:Bs+g+g, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
!                    !hvy_block( g+1:g+(Bs+1)/2, Bs+g+1-rmv_redundant:Bs+g+g, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
!                    hvy_block( g+(Bs+1)/2:Bs+g, Bs+g+1-rmv_redundant:Bs+g+g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if

        case(35,36,37,38)
            if ( level_diff == -1 ) then
                ! x
                sender_x(1)   = (Bs+1)/2
                sender_x(2)   = Bs+g
                receiver_x(1) = 1
                receiver_x(2) = g
                finer_x(1)    = Bs+g
                finer_x(2)    = Bs-1+2*g
                ! z,y
                select case(neighborhood)
                    case(35) ! '__3/123'
                        sender_z(1)   = (Bs+1)/2
                        sender_z(2)   = Bs+g
                        receiver_z(1) = 1
                        receiver_z(2) = Bs+g
                        finer_z(1)    = g+1
                        finer_z(2)    = Bs+2*g

                        sender_y(1)   = g+1
                        sender_y(2)   = (Bs+1)/2+g+g
                        receiver_y(1) = g+1
                        receiver_y(2) = Bs+2*g
                        finer_y(1)    = 1
                        finer_y(2)    = Bs+g

                    case(37) ! '__3/134'
                        sender_z(1)   = (Bs+1)/2
                        sender_z(2)   = Bs+g
                        receiver_z(1) = 1
                        receiver_z(2) = Bs+g
                        finer_z(1)    = g+1
                        finer_z(2)    = Bs+2*g

                        sender_y(1)   = (Bs+1)/2
                        sender_y(2)   = Bs+g
                        receiver_y(1) = 1
                        receiver_y(2) = Bs+g
                        finer_y(1)    = g+1
                        finer_y(2)    = Bs+2*g

                    case(38) ! '__3/634'
                        sender_z(1)   = g+1
                        sender_z(2)   = (Bs+1)/2+g+g
                        receiver_z(1) = g+1
                        receiver_z(2) = Bs+2*g
                        finer_z(1)    = 1
                        finer_z(2)    = Bs+g

                        sender_y(1)   = (Bs+1)/2
                        sender_y(2)   = Bs+g
                        receiver_y(1) = 1
                        receiver_y(2) = Bs+g
                        finer_y(1)    = g+1
                        finer_y(2)    = Bs+2*g

                    case(36) ! '__3/623'
                        sender_z(1)   = g+1
                        sender_z(2)   = (Bs+1)/2+g+g
                        receiver_z(1) = g+1
                        receiver_z(2) = Bs+2*g
                        finer_z(1)    = 1
                        finer_z(2)    = Bs+g

                        sender_y(1)   = g+1
                        sender_y(2)   = (Bs+1)/2+g+g
                        receiver_y(1) = g+1
                        receiver_y(2) = Bs+2*g
                        finer_y(1)    = 1
                        finer_y(2)    = Bs+g

                end select

                ! copy data
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    data_face = hvy_block( sender_x(1):sender_x(2), sender_y(1):sender_y(2), sender_z(1):sender_z(2), dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = data_face_fine(finer_x(1):finer_x(2), finer_y(1):finer_y(2), finer_z(1):finer_z(2))
                end do

            elseif ( level_diff == 1 ) then
                ! x,y,z
                sender_x(1)   = Bs-g
                sender_x(2)   = Bs+g-2+2*rmv_redundant
                sender_y(1)   = g+1
                sender_y(2)   = Bs+g
                sender_z(1)   = g+1
                sender_z(2)   = Bs+g

                receiver_x(1) = 1
                receiver_x(2) = g+rmv_redundant
                ! z,y
                select case(neighborhood)
                    case(35) ! '__3/123'
                        receiver_z(1) = g+(Bs+1)/2
                        receiver_z(2) = Bs+g

                        receiver_y(1) = g+1
                        receiver_y(2) = g+(Bs+1)/2

                    case(37) ! '__3/134'
                        receiver_z(1) = g+(Bs+1)/2
                        receiver_z(2) = Bs+g

                        receiver_y(1) = g+(Bs+1)/2
                        receiver_y(2) = Bs+g

                    case(38) ! '__3/634'
                        receiver_z(1) = g+1
                        receiver_z(2) = g+(Bs+1)/2

                        receiver_y(1) = g+(Bs+1)/2
                        receiver_y(2) = Bs+g

                    case(36) ! '__3/623'
                        receiver_z(1) = g+1
                        receiver_z(2) = g+(Bs+1)/2

                        receiver_y(1) = g+1
                        receiver_y(2) = g+(Bs+1)/2

                end select

                ! copy data
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = hvy_block( sender_x(1):sender_x(2):2, sender_y(1):sender_y(2):2, sender_z(1):sender_z(2):2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

!        ! '__3/123'
!        case(35)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate
!                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
!                    !data_face = hvy_block( (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, dF, sender_id )
!                    !data_face = hvy_block( (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!                    data_face = hvy_block( (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( 1:g, g+1:Bs+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, g+1:Bs+g, g+1:Bs+g)
!                    !hvy_block( 1:g, g+1:Bs+2*g, 1:Bs+g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, 1:Bs+g, 1:Bs+g)
!                    !hvy_block( 1:g, 1:Bs+g, g+1:Bs+2*g, dF, receiver_id ) = data_face_fine(Bs+g:Bs-1+2*g, g+1:Bs+2*g, 1:Bs+g)
!                    !hvy_block( 1:g, 1:Bs+g, 1:Bs+g, dF, receiver_id ) = 17.0_rk!data_face_fine(Bs+g:Bs-1+2*g, g+1:Bs+2*g, g+1:Bs+2*g)
!                    hvy_block( 1:Bs+g, g+1:Bs+2*g, Bs+g+1:Bs+g+g, dF, receiver_id ) = 17.0_rk
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( 1:g, g+1:g+(Bs+1)/2, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
!                    !hvy_block( 1:g+rmv_redundant, g+1:g+(Bs+1)/2, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
!                    !hvy_block( 1:g+rmv_redundant, g+(Bs+1)/2:Bs+g, g+1:g+(Bs+1)/2, dF, receiver_id ) = 17.0_rk!hvy_block( Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
!                    hvy_block( g+(Bs+1)/2:Bs+g, g+1:g+(Bs+1)/2, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = 17.0_rk
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if

!        ! '__3/623'
!        case(36)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate
!                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
!                    !data_face = hvy_block( (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!                    data_face = hvy_block( (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( 1:g, g+1:Bs+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, g+1:Bs+g, g+1:Bs+g)
!                    !hvy_block( 1:g, g+1:Bs+2*g, g+1:Bs+2*g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, 1:Bs+g, 1:Bs+g)
!                    !hvy_block( 1:g, g+1:Bs+2*g, g+1:Bs+2*g, dF, receiver_id ) = data_face_fine(Bs+g:Bs-1+2*g, 1:Bs+g, 1:Bs+g)
!                    !hvy_block( 1:g, 1:Bs+g, 1:Bs+g, dF, receiver_id ) = data_face_fine(Bs+g:Bs-1+2*g, g+1:Bs+2*g, g+1:Bs+2*g)
!                    hvy_block( 1:g, 1:Bs+g, g+1:Bs+2*g, dF, receiver_id ) = data_face_fine(Bs+g:Bs-1+2*g, g+1:Bs+2*g, g+1:Bs+2*g)
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( 1:g, g+1:g+(Bs+1)/2, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
!                    !hvy_block( 1:g+rmv_redundant, g+1:g+(Bs+1)/2, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
!                    hvy_block( 1:g+rmv_redundant, g+(Bs+1)/2:Bs+g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if
!        ! '__3/134'
!        case(37)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate
!                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
!                    !data_face = hvy_block( (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, dF, sender_id )
!                    data_face = hvy_block( (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( 1:g, g+1:Bs+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, g+1:Bs+g, g+1:Bs+g)
!                    !hvy_block( 1:g, 1:Bs+g, 1:Bs+g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, 1:Bs+g, 1:Bs+g)
!                    !hvy_block( 1:g, g+1:Bs+2*g, g+1:Bs+2*g, dF, receiver_id ) = data_face_fine(Bs+g:Bs-1+2*g, 1:Bs+g, 1:Bs+g)
!                    hvy_block( 1:g, g+1:Bs+2*g, 1:Bs+g, dF, receiver_id ) = data_face_fine(Bs+g:Bs-1+2*g, 1:Bs+g, 1:Bs+g)
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( 1:g, g+(Bs+1)/2:Bs+g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
!                    !hvy_block( 1:g+rmv_redundant, g+(Bs+1)/2:Bs+g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
!                    hvy_block( 1:g+rmv_redundant, g+1:g+(Bs+1)/2, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if
!        ! '__3/634'
!        case(38)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate
!                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
!                    !data_face = hvy_block( (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!                    data_face = hvy_block( (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( 1:g, g+1:Bs+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, g+1:Bs+g, g+1:Bs+g)
!                    !hvy_block( 1:g, 1:Bs+g, g+1:Bs+2*g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, 1:Bs+g, 1:Bs+g)
!                    !hvy_block( 1:g, 1:Bs+g, g+1:Bs+2*g, dF, receiver_id ) = data_face_fine(Bs+g:Bs-1+2*g, g+1:Bs+2*g, 1:Bs+g)
!                    !hvy_block( 1:g, g+1:Bs+2*g, 1:Bs+g, dF, receiver_id ) = data_face_fine(Bs+g:Bs-1+2*g, 1:Bs+g, g+1:Bs+2*g)
!                    hvy_block( 1:g, g+1:Bs+2*g, g+1:Bs+2*g, dF, receiver_id ) = data_face_fine(Bs+g:Bs-1+2*g, 1:Bs+g, g+1:Bs+2*g)
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( 1:g, g+(Bs+1)/2:Bs+g, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
!                    !hvy_block( 1:g+rmv_redundant, g+(Bs+1)/2:Bs+g, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
!                    hvy_block( 1:g+rmv_redundant, g+1:g+(Bs+1)/2, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if

        case(39,40,41,42)
            if ( level_diff == -1 ) then
                ! y
                sender_y(1)   = (Bs+1)/2
                sender_y(2)   = Bs+g
                receiver_y(1) = 1
                receiver_y(2) = g
                finer_y(1)    = Bs+g
                finer_y(2)    = Bs-1+2*g
                ! x,z
                select case(neighborhood)
                    case(40) ! '__4/634'
                        sender_x(1)   = (Bs+1)/2
                        sender_x(2)   = Bs+g
                        receiver_x(1) = 1
                        receiver_x(2) = Bs+g
                        finer_x(1)    = g+1
                        finer_x(2)    = Bs+2*g

                        sender_z(1)   = g+1
                        sender_z(2)   = (Bs+1)/2+g+g
                        receiver_z(1) = g+1
                        receiver_z(2) = Bs+2*g
                        finer_z(1)    = 1
                        finer_z(2)    = Bs+g

                    case(39) ! '__4/134'
                        sender_x(1)   = (Bs+1)/2
                        sender_x(2)   = Bs+g
                        receiver_x(1) = 1
                        receiver_x(2) = Bs+g
                        finer_x(1)    = g+1
                        finer_x(2)    = Bs+2*g

                        sender_z(1)   = (Bs+1)/2
                        sender_z(2)   = Bs+g
                        receiver_z(1) = 1
                        receiver_z(2) = Bs+g
                        finer_z(1)    = g+1
                        finer_z(2)    = Bs+2*g

                    case(41) ! '__4/145'
                        sender_x(1)   = g+1
                        sender_x(2)   = (Bs+1)/2+g+g
                        receiver_x(1) = g+1
                        receiver_x(2) = Bs+2*g
                        finer_x(1)    = 1
                        finer_x(2)    = Bs+g

                        sender_z(1)   = (Bs+1)/2
                        sender_z(2)   = Bs+g
                        receiver_z(1) = 1
                        receiver_z(2) = Bs+g
                        finer_z(1)    = g+1
                        finer_z(2)    = Bs+2*g

                    case(42) ! '__4/645'
                        sender_x(1)   = g+1
                        sender_x(2)   = (Bs+1)/2+g+g
                        receiver_x(1) = g+1
                        receiver_x(2) = Bs+2*g
                        finer_x(1)    = 1
                        finer_x(2)    = Bs+g

                        sender_z(1)   = g+1
                        sender_z(2)   = (Bs+1)/2+g+g
                        receiver_z(1) = g+1
                        receiver_z(2) = Bs+2*g
                        finer_z(1)    = 1
                        finer_z(2)    = Bs+g

                end select

                ! copy data
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    data_face = hvy_block( sender_x(1):sender_x(2), sender_y(1):sender_y(2), sender_z(1):sender_z(2), dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = data_face_fine(finer_x(1):finer_x(2), finer_y(1):finer_y(2), finer_z(1):finer_z(2))
                end do

            elseif ( level_diff == 1 ) then
                ! x,y,z
                sender_y(1)   = Bs-g
                sender_y(2)   = Bs+g-2+2*rmv_redundant
                sender_z(1)   = g+1
                sender_z(2)   = Bs+g
                sender_x(1)   = g+1
                sender_x(2)   = Bs+g

                receiver_y(1) = 1
                receiver_y(2) = g+rmv_redundant
                ! x,y
                select case(neighborhood)
                    case(40) ! '__4/634'
                        receiver_x(1) = g+(Bs+1)/2
                        receiver_x(2) = Bs+g

                        receiver_z(1) = g+1
                        receiver_z(2) = g+(Bs+1)/2

                    case(39) ! '__4/134'
                        receiver_x(1) = g+(Bs+1)/2
                        receiver_x(2) = Bs+g

                        receiver_z(1) = g+(Bs+1)/2
                        receiver_z(2) = Bs+g

                    case(41) ! '__4/145'
                        receiver_x(1) = g+1
                        receiver_x(2) = g+(Bs+1)/2

                        receiver_z(1) = g+(Bs+1)/2
                        receiver_z(2) = Bs+g

                    case(42) ! '__4/645'
                        receiver_x(1) = g+1
                        receiver_x(2) = g+(Bs+1)/2

                        receiver_z(1) = g+1
                        receiver_z(2) = g+(Bs+1)/2

                end select

                ! copy data
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = hvy_block( sender_x(1):sender_x(2):2, sender_y(1):sender_y(2):2, sender_z(1):sender_z(2):2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

!        ! '__4/134'
!        case(39)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate
!                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
!                    !data_face = hvy_block( (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, dF, sender_id )
!                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( g+1:Bs+g, 1:g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, Bs:Bs-1+g, g+1:Bs+g)
!                    !hvy_block( 1:Bs+g, 1:g, 1:Bs+g, dF, receiver_id ) = data_face_fine(1:Bs+g, Bs:Bs-1+g, 1:Bs+g)
!                    hvy_block( g+1:Bs+2*g, 1:g, g+1:Bs+2*g, dF, receiver_id ) = data_face_fine(1:Bs+g, Bs+g:Bs-1+2*g, 1:Bs+g)
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( g+(Bs+1)/2:Bs+g, 1:g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2:2, g+1:Bs+g:2, dF, sender_id )
!                    !hvy_block( g+(Bs+1)/2:Bs+g, 1:g+rmv_redundant, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, dF, sender_id )
!                    hvy_block( g+1:g+(Bs+1)/2, 1:g+rmv_redundant, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if
!
!        ! '__4/634'
!        case(40)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate
!                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
!                    !data_face = hvy_block( (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( g+1:Bs+g, 1:g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, Bs:Bs-1+g, g+1:Bs+g)
!                    !hvy_block( 1:Bs+g, 1:g, g+1:Bs+2*g, dF, receiver_id ) = data_face_fine(1:Bs+g, Bs:Bs-1+g, 1:Bs+g)
!                    hvy_block( g+1:Bs+2*g, 1:g, 1:Bs+g, dF, receiver_id ) = data_face_fine(1:Bs+g, Bs+g:Bs-1+2*g, g+1:Bs+2*g)
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( g+(Bs+1)/2:Bs+g, 1:g, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2:2, g+1:Bs+g:2, dF, sender_id )
!                    !hvy_block( g+(Bs+1)/2:Bs+g, 1:g+rmv_redundant, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, dF, sender_id )
!                    hvy_block( g+1:g+(Bs+1)/2, 1:g+rmv_redundant, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if
!
!        ! '__4/145'
!        case(41)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate
!                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
!                    !data_face = hvy_block( g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, dF, sender_id )
!                    data_face = hvy_block( (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( g+1:Bs+g, 1:g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, Bs:Bs-1+g, g+1:Bs+g)
!                    !hvy_block( g+1:Bs+2*g, 1:g, 1:Bs+g, dF, receiver_id ) = data_face_fine(1:Bs+g, Bs:Bs-1+g, 1:Bs+g)
!                    hvy_block( 1:Bs+g, 1:g, g+1:Bs+2*g, dF, receiver_id ) = data_face_fine(g+1:Bs+2*g, Bs+g:Bs-1+2*g, 1:Bs+g)
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( g+1:g+(Bs+1)/2, 1:g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2:2, g+1:Bs+g:2, dF, sender_id )
!                    !hvy_block( g+1:g+(Bs+1)/2, 1:g+rmv_redundant, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, dF, sender_id )
!                    hvy_block( g+(Bs+1)/2:Bs+g, 1:g+rmv_redundant, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if
!
!        ! '__4/645'
!        case(42)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate
!                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
!                    !data_face = hvy_block( g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!                    data_face = hvy_block( (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( g+1:Bs+g, 1:g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, Bs:Bs-1+g, g+1:Bs+g)
!                    !hvy_block( g+1:Bs+2*g, 1:g, g+1:Bs+2*g, dF, receiver_id ) = data_face_fine(1:Bs+g, Bs:Bs-1+g, 1:Bs+g)
!                    hvy_block( 1:Bs+g, 1:g, 1:Bs+g, dF, receiver_id ) = data_face_fine(g+1:Bs+2*g, Bs+g:Bs-1+2*g, g+1:Bs+2*g)
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( g+1:g+(Bs+1)/2, 1:g, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2:2, g+1:Bs+g:2, dF, sender_id )
!                    !hvy_block( g+1:g+(Bs+1)/2, 1:g+rmv_redundant, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, dF, sender_id )
!                    hvy_block( g+(Bs+1)/2:Bs+g, 1:g+rmv_redundant, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if

        case(43,44,45,46)
            if ( level_diff == -1 ) then
                ! x
                sender_x(1)   = g+1
                sender_x(2)   = (Bs+1)/2+g+g
                receiver_x(1) = Bs+g+1
                receiver_x(2) = Bs+g+g
                finer_x(1)    = 2
                finer_x(2)    = g+1
                ! z,y
                select case(neighborhood)
                    case(45) ! '__5/152'
                        sender_z(1)   = (Bs+1)/2
                        sender_z(2)   = Bs+g
                        receiver_z(1) = 1
                        receiver_z(2) = Bs+g
                        finer_z(1)    = g+1
                        finer_z(2)    = Bs+2*g

                        sender_y(1)   = g+1
                        sender_y(2)   = (Bs+1)/2+g+g
                        receiver_y(1) = g+1
                        receiver_y(2) = Bs+2*g
                        finer_y(1)    = 1
                        finer_y(2)    = Bs+g

                    case(43) ! '__5/145'
                        sender_z(1)   = (Bs+1)/2
                        sender_z(2)   = Bs+g
                        receiver_z(1) = 1
                        receiver_z(2) = Bs+g
                        finer_z(1)    = g+1
                        finer_z(2)    = Bs+2*g

                        sender_y(1)   = (Bs+1)/2
                        sender_y(2)   = Bs+g
                        receiver_y(1) = 1
                        receiver_y(2) = Bs+g
                        finer_y(1)    = g+1
                        finer_y(2)    = Bs+2*g

                    case(44) ! '__5/645'
                        sender_z(1)   = g+1
                        sender_z(2)   = (Bs+1)/2+g+g
                        receiver_z(1) = g+1
                        receiver_z(2) = Bs+2*g
                        finer_z(1)    = 1
                        finer_z(2)    = Bs+g

                        sender_y(1)   = (Bs+1)/2
                        sender_y(2)   = Bs+g
                        receiver_y(1) = 1
                        receiver_y(2) = Bs+g
                        finer_y(1)    = g+1
                        finer_y(2)    = Bs+2*g

                    case(46) ! '__5/652'
                        sender_z(1)   = g+1
                        sender_z(2)   = (Bs+1)/2+g+g
                        receiver_z(1) = g+1
                        receiver_z(2) = Bs+2*g
                        finer_z(1)    = 1
                        finer_z(2)    = Bs+g

                        sender_y(1)   = g+1
                        sender_y(2)   = (Bs+1)/2+g+g
                        receiver_y(1) = g+1
                        receiver_y(2) = Bs+2*g
                        finer_y(1)    = 1
                        finer_y(2)    = Bs+g

                end select

                ! copy data
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    data_face = hvy_block( sender_x(1):sender_x(2), sender_y(1):sender_y(2), sender_z(1):sender_z(2), dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = data_face_fine(finer_x(1):finer_x(2), finer_y(1):finer_y(2), finer_z(1):finer_z(2))
                end do

            elseif ( level_diff == 1 ) then
                ! x,y,z
                sender_x(1)   = g+3-2*rmv_redundant
                sender_x(2)   = 3*g+1
                sender_y(1)   = g+1
                sender_y(2)   = Bs+g
                sender_z(1)   = g+1
                sender_z(2)   = Bs+g

                receiver_x(1) = Bs+g+1-rmv_redundant
                receiver_x(2) = Bs+g+g
                ! z,y
                select case(neighborhood)
                    case(45) ! '__5/152'
                        receiver_z(1) = g+(Bs+1)/2
                        receiver_z(2) = Bs+g

                        receiver_y(1) = g+1
                        receiver_y(2) = g+(Bs+1)/2

                    case(43) ! '__5/145'
                        receiver_z(1) = g+(Bs+1)/2
                        receiver_z(2) = Bs+g

                        receiver_y(1) = g+(Bs+1)/2
                        receiver_y(2) = Bs+g

                    case(44) ! '__5/645'
                        receiver_z(1) = g+1
                        receiver_z(2) = g+(Bs+1)/2

                        receiver_y(1) = g+(Bs+1)/2
                        receiver_y(2) = Bs+g

                    case(46) ! '__5/652'
                        receiver_z(1) = g+1
                        receiver_z(2) = g+(Bs+1)/2

                        receiver_y(1) = g+1
                        receiver_y(2) = g+(Bs+1)/2

                end select

                ! copy data
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = hvy_block( sender_x(1):sender_x(2):2, sender_y(1):sender_y(2):2, sender_z(1):sender_z(2):2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

!        ! '__5/145'
!        case(43)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate
!                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
!                    !data_face = hvy_block( g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, dF, sender_id )
!                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(2:g+1, g+1:Bs+g, g+1:Bs+g)
!                    !hvy_block( Bs+g+1:Bs+g+g, 1:Bs+g, 1:Bs+g, dF, receiver_id ) = data_face_fine(2:g+1, 1:Bs+g, 1:Bs+g)
!                    hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+2*g, g+1:Bs+2*g, dF, receiver_id ) = data_face_fine(2:g+1, 1:Bs+g, 1:Bs+g)
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( Bs+g+1:Bs+g+g, g+(Bs+1)/2:Bs+g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+3:3*g+1:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
!                    !hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, g+(Bs+1)/2:Bs+g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
!                    hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, g+1:g+(Bs+1)/2, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if
!
!        ! '__5/645'
!        case(44)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate
!                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
!                    !data_face = hvy_block( g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(2:g+1, g+1:Bs+g, g+1:Bs+g)
!                    !hvy_block( Bs+g+1:Bs+g+g, 1:Bs+g, g+1:Bs+2*g, dF, receiver_id ) = data_face_fine(2:g+1, 1:Bs+g, 1:Bs+g)
!                    !hvy_block( Bs+g+1:Bs+g+g, 1:Bs+g, g+1:Bs+2*g, dF, receiver_id ) = data_face_fine(2:g+1, g+1:Bs+2*g, 1:Bs+g)
!                    hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+2*g, 1:Bs+g, dF, receiver_id ) = data_face_fine(2:g+1, 1:Bs+g, g+1:Bs+2*g)
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( Bs+g+1:Bs+g+g, g+(Bs+1)/2:Bs+g, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+3:3*g+1:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
!                    !hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, g+(Bs+1)/2:Bs+g, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
!                    hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, g+1:g+(Bs+1)/2, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if
!
!        ! '__5/152'
!        case(45)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate
!                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
!                    !data_face = hvy_block( g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, dF, sender_id )
!                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(2:g+1, g+1:Bs+g, g+1:Bs+g)
!                    !hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+2*g, 1:Bs+g, dF, receiver_id ) = data_face_fine(2:g+1, 1:Bs+g, 1:Bs+g)
!                    hvy_block( Bs+g+1:Bs+g+g, 1:Bs+g, g+1:Bs+2*g, dF, receiver_id ) = data_face_fine(2:g+1, g+1:Bs+2*g, 1:Bs+g)
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( Bs+g+1:Bs+g+g, g+1:g+(Bs+1)/2, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+3:3*g+1:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
!                    !hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, g+1:g+(Bs+1)/2, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
!                    hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, g+(Bs+1)/2:Bs+g, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if
!        ! '__5/652'
!        case(46)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate
!                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
!                    !data_face = hvy_block( g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(2:g+1, g+1:Bs+g, g+1:Bs+g)
!                    !hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+2*g, g+1:Bs+2*g, dF, receiver_id ) = data_face_fine(2:g+1, 1:Bs+g, 1:Bs+g)
!                    hvy_block( Bs+g+1:Bs+g+g, 1:Bs+g, 1:Bs+g, dF, receiver_id ) = data_face_fine(2:g+1, g+1:Bs+2*g, g+1:Bs+2*g)
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( Bs+g+1:Bs+g+g, g+1:g+(Bs+1)/2, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+3:3*g+1:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
!                    !hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, g+1:g+(Bs+1)/2, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
!                    hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, g+(Bs+1)/2:Bs+g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if

        case(47,48,49,50)
            if ( level_diff == -1 ) then
                ! z
                sender_z(1)   = g+1
                sender_z(2)   = (Bs+1)/2+g+g
                receiver_z(1) = Bs+g+1
                receiver_z(2) = Bs+g+g
                finer_z(1)    = 2
                finer_z(2)    = g+1
                ! x,y
                select case(neighborhood)
                    case(47) ! '__6/623'
                        sender_x(1)   = (Bs+1)/2
                        sender_x(2)   = Bs+g
                        receiver_x(1) = 1
                        receiver_x(2) = Bs+g
                        finer_x(1)    = g+1
                        finer_x(2)    = Bs+2*g

                        sender_y(1)   = g+1
                        sender_y(2)   = (Bs+1)/2+g+g
                        receiver_y(1) = g+1
                        receiver_y(2) = Bs+2*g
                        finer_y(1)    = 1
                        finer_y(2)    = Bs+g

                    case(48) ! '__6/634'
                        sender_x(1)   = (Bs+1)/2
                        sender_x(2)   = Bs+g
                        receiver_x(1) = 1
                        receiver_x(2) = Bs+g
                        finer_x(1)    = g+1
                        finer_x(2)    = Bs+2*g

                        sender_y(1)   = (Bs+1)/2
                        sender_y(2)   = Bs+g
                        receiver_y(1) = 1
                        receiver_y(2) = Bs+g
                        finer_y(1)    = g+1
                        finer_y(2)    = Bs+2*g

                    case(49) ! '__6/645'
                        sender_x(1)   = g+1
                        sender_x(2)   = (Bs+1)/2+g+g
                        receiver_x(1) = g+1
                        receiver_x(2) = Bs+2*g
                        finer_x(1)    = 1
                        finer_x(2)    = Bs+g

                        sender_y(1)   = (Bs+1)/2
                        sender_y(2)   = Bs+g
                        receiver_y(1) = 1
                        receiver_y(2) = Bs+g
                        finer_y(1)    = g+1
                        finer_y(2)    = Bs+2*g

                    case(50) ! '__6/652'
                        sender_x(1)   = g+1
                        sender_x(2)   = (Bs+1)/2+g+g
                        receiver_x(1) = g+1
                        receiver_x(2) = Bs+2*g
                        finer_x(1)    = 1
                        finer_x(2)    = Bs+g

                        sender_y(1)   = g+1
                        sender_y(2)   = (Bs+1)/2+g+g
                        receiver_y(1) = g+1
                        receiver_y(2) = Bs+2*g
                        finer_y(1)    = 1
                        finer_y(2)    = Bs+g

                end select

                ! copy data
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    data_face = hvy_block( sender_x(1):sender_x(2), sender_y(1):sender_y(2), sender_z(1):sender_z(2), dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = data_face_fine(finer_x(1):finer_x(2), finer_y(1):finer_y(2), finer_z(1):finer_z(2))
                end do

            elseif ( level_diff == 1 ) then
                ! x,y,z
                sender_z(1)   = g+3-2*rmv_redundant
                sender_z(2)   = 3*g+1
                sender_y(1)   = g+1
                sender_y(2)   = Bs+g
                sender_x(1)   = g+1
                sender_x(2)   = Bs+g

                receiver_z(1) = Bs+g+1-rmv_redundant
                receiver_z(2) = Bs+g+g
                ! x,y
                select case(neighborhood)
                    case(47) ! '__6/623'
                        receiver_x(1) = g+(Bs+1)/2
                        receiver_x(2) = Bs+g

                        receiver_y(1) = g+1
                        receiver_y(2) = g+(Bs+1)/2

                    case(48) ! '__6/634'
                        receiver_x(1) = g+(Bs+1)/2
                        receiver_x(2) = Bs+g

                        receiver_y(1) = g+(Bs+1)/2
                        receiver_y(2) = Bs+g

                    case(49) ! '__6/645'
                        receiver_x(1) = g+1
                        receiver_x(2) = g+(Bs+1)/2

                        receiver_y(1) = g+(Bs+1)/2
                        receiver_y(2) = Bs+g

                    case(50) ! '__6/652'
                        receiver_x(1) = g+1
                        receiver_x(2) = g+(Bs+1)/2

                        receiver_y(1) = g+1
                        receiver_y(2) = g+(Bs+1)/2

                end select

                ! copy data
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = hvy_block( sender_x(1):sender_x(2):2, sender_y(1):sender_y(2):2, sender_z(1):sender_z(2):2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

!        ! '__6/623'
!        case(47)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate
!                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g , g+1:(Bs+1)/2+g/2+g, dF, sender_id )
!                    data_face = hvy_block( (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( g+1:Bs+g, g+1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, g+1:Bs+g, 2:g+1)
!                    !hvy_block( 1:Bs+g, g+1:Bs+2*g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(1:Bs+g, 1:Bs+g, 2:g+1)
!                    hvy_block( g+1:Bs+2*g, 1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(g+1:Bs+2*g, 1:Bs+g, 2:g+1)
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( g+(Bs+1)/2:Bs+g, g+1:g+(Bs+1)/2, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, g+3:3*g+1:2, dF, sender_id )
!                    !hvy_block( g+(Bs+1)/2:Bs+g, g+1:g+(Bs+1)/2, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, dF, sender_id )
!                    hvy_block( g+1:g+(Bs+1)/2, g+(Bs+1)/2:Bs+g, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if
!
!        ! '__6/634'
!        case(48)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate
!                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
!                    !data_face = hvy_block( (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( g+1:Bs+g, g+1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, g+1:Bs+g, 2:g+1)
!                    !hvy_block( 1:Bs+g, 1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(1:Bs+g, 1:Bs+g, 2:g+1)
!                    !hvy_block( 1:Bs+g, 1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(g+1:Bs+2*g, g+1:Bs+2*g, 2:g+1)
!                    hvy_block( g+1:Bs+2*g, g+1:Bs+2*g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(1:Bs+g, 1:Bs+g, 2:g+1)
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( g+(Bs+1)/2:Bs+g, g+(Bs+1)/2:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, g+3:3*g+1:2, dF, sender_id )
!                    !hvy_block( g+(Bs+1)/2:Bs+g, g+(Bs+1)/2:Bs+g, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, dF, sender_id )
!                    hvy_block( g+1:g+(Bs+1)/2, g+1:g+(Bs+1)/2, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if
!        ! '__6/645'
!        case(49)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate
!                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
!                    !data_face = hvy_block( g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!                    data_face = hvy_block( (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( g+1:Bs+g, g+1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, g+1:Bs+g, 2:g+1)
!                    !hvy_block( g+1:Bs+2*g, 1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(1:Bs+g, 1:Bs+g, 2:g+1)
!                    !hvy_block( g+1:Bs+2*g, 1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(1:Bs+g, g+1:Bs+2*g, 2:g+1)
!                    hvy_block( 1:Bs+g, g+1:Bs+2*g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(g+1:Bs+2*g, 1:Bs+g, 2:g+1)
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( g+1:g+(Bs+1)/2, g+(Bs+1)/2:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, g+3:3*g+1:2, dF, sender_id )
!                    !hvy_block( g+1:g+(Bs+1)/2, g+(Bs+1)/2:Bs+g, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, dF, sender_id )
!                    hvy_block( g+(Bs+1)/2:Bs+g, g+1:g+(Bs+1)/2, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if
!
!        ! '__6/652'
!        case(50)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate
!                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
!                    !data_face = hvy_block( g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!                    data_face = hvy_block( (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( g+1:Bs+g, g+1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, g+1:Bs+g, 2:g+1)
!                    !hvy_block( g+1:Bs+2*g, g+1:Bs+2*g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(1:Bs+g, 1:Bs+g, 2:g+1)
!                    hvy_block( 1:Bs+g, 1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(g+1:Bs+2*g, g+1:Bs+2*g, 2:g+1)
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( g+1:g+(Bs+1)/2, g+1:g+(Bs+1)/2, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, g+3:3*g+1:2, dF, sender_id )
!                    !hvy_block( g+1:g+(Bs+1)/2, g+1:g+(Bs+1)/2, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, dF, sender_id )
!                    hvy_block( g+(Bs+1)/2:Bs+g, g+(Bs+1)/2:Bs+g, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if

        case(51,52)
            if ( level_diff == -1 ) then
                ! z,y
                sender_z(1)   = (Bs+1)/2
                sender_z(2)   = Bs+g
                receiver_z(1) = 1
                receiver_z(2) = g
                finer_z(1)    = Bs+g
                finer_z(2)    = Bs-1+2*g

                sender_y(1)   = g+1
                sender_y(2)   = (Bs+1)/2+g+g
                receiver_y(1) = Bs+g+1
                receiver_y(2) = Bs+2*g
                finer_y(1)    = 2
                finer_y(2)    = g+1
                ! x,y
                select case(neighborhood)
                    case(51) ! '_12/123'
                        sender_x(1)   = (Bs+1)/2
                        sender_x(2)   = Bs+g
                        receiver_x(1) = 1
                        receiver_x(2) = Bs+g
                        finer_x(1)    = g+1
                        finer_x(2)    = Bs+2*g

                    case(52) ! '_12/152'
                        sender_x(1)   = g+1
                        sender_x(2)   = (Bs+1)/2+g+g
                        receiver_x(1) = g+1
                        receiver_x(2) = Bs+2*g
                        finer_x(1)    = 1
                        finer_x(2)    = Bs+g

                end select

                ! copy data
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    data_face = hvy_block( sender_x(1):sender_x(2), sender_y(1):sender_y(2), sender_z(1):sender_z(2), dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = data_face_fine(finer_x(1):finer_x(2), finer_y(1):finer_y(2), finer_z(1):finer_z(2))
                end do

            elseif ( level_diff == 1 ) then
                ! x,y,z
                sender_z(1)   = Bs-g
                sender_z(2)   = Bs+g-2+2*rmv_redundant
                sender_y(1)   = g+3-2*rmv_redundant
                sender_y(2)   = 3*g+1
                sender_x(1)   = g+1
                sender_x(2)   = Bs+g

                receiver_z(1) = 1
                receiver_z(2) = g+rmv_redundant
                receiver_y(1) = Bs+g+1-rmv_redundant
                receiver_y(2) = Bs+g+g
                ! y
                select case(neighborhood)
                    case(51) ! '_12/123'
                        receiver_x(1) = g+(Bs+1)/2
                        receiver_x(2) = Bs+g

                    case(52) ! '_12/152'
                        receiver_x(1) = g+1
                        receiver_x(2) = g+(Bs+1)/2

                end select

                ! copy data
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = hvy_block( sender_x(1):sender_x(2):2, sender_y(1):sender_y(2):2, sender_z(1):sender_z(2):2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

!!        ! '_12/123'
!!        case(51)
!!            if ( level_diff == -1 ) then
!!                ! sender on lower level
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    ! data to interpolate, note: use data_face interpolation variable
!!                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
!!                    data_face = hvy_block( (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, dF, sender_id )
!!                    ! interpolate data
!!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!!                    ! copy data
!!                    !hvy_block( g+1:Bs+g, Bs+g+1:Bs+g+g, 1:g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, 2:g+1, Bs:Bs-1+g)
!!                    hvy_block( g+1:Bs+g, Bs+g+1:Bs+g+g, 1:g, dF, receiver_id ) = data_face_fine(g+g+1:Bs+2*g, 2:g+1, Bs+g:Bs-1+2*g)
!!                end do
!!
!!            elseif ( level_diff == 1 ) then
!!                ! sender on higher level
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    !hvy_block( g+(Bs+1)/2:Bs+g, Bs+g+1:Bs+g+g, 1:g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3:3*g+1:2, Bs-g:Bs+g-2:2, dF, sender_id )
!!                    hvy_block( g+(Bs+1)/2:Bs+g, Bs+g+1-rmv_redundant:Bs+g+g, 1:g+rmv_redundant, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, Bs-g:Bs+g-2+2*rmv_redundant:2, dF, sender_id )
!!                end do
!!
!!            else
!!                ! error case
!!                write(*,'(80("_"))')
!!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!!                stop
!!            end if
!!
!!        ! '_12/152'
!!        case(52)
!!            if ( level_diff == -1 ) then
!!                ! sender on lower level
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    ! data to interpolate, note: use data_face interpolation variable
!!                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
!!                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, dF, sender_id )
!!                    ! interpolate data
!!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!!                    ! copy data
!!                    !hvy_block( g+1:Bs+g, Bs+g+1:Bs+g+g, 1:g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, 2:g+1, Bs:Bs-1+g)
!!                    !hvy_block( g+1:Bs+g, Bs+g+1:Bs+g+g, 1:g, dF, receiver_id ) = data_face_fine(1:Bs, 2:g+1, Bs:Bs-1+g)
!!                    hvy_block( g+1:Bs+g, Bs+g+1:Bs+g+g, 1:g, dF, receiver_id ) = data_face_fine(1:Bs, 2:g+1, Bs+g:Bs-1+2*g)
!!                end do
!!
!!            elseif ( level_diff == 1 ) then
!!                ! sender on higher level
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    !hvy_block( g+1:g+(Bs+1)/2, Bs+g+1:Bs+g+g, 1:g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3:3*g+1:2, Bs-g:Bs+g-2:2, dF, sender_id )
!!                    hvy_block( g+1:g+(Bs+1)/2, Bs+g+1-rmv_redundant:Bs+g+g, 1:g+rmv_redundant, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, Bs-g:Bs+g-2+2*rmv_redundant:2, dF, sender_id )
!!                end do
!!
!!            else
!!                ! error case
!!                write(*,'(80("_"))')
!!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!!                stop
!!            end if
        case(53,54)
            if ( level_diff == -1 ) then
                ! z,x
                sender_z(1)   = (Bs+1)/2
                sender_z(2)   = Bs+g
                receiver_z(1) = 1
                receiver_z(2) = g
                finer_z(1)    = Bs+g
                finer_z(2)    = Bs-1+2*g

                sender_x(1)   = (Bs+1)/2
                sender_x(2)   = Bs+g
                receiver_x(1) = 1
                receiver_x(2) = g
                finer_x(1)    = Bs+g
                finer_x(2)    = Bs-1+2*g
                ! x,y
                select case(neighborhood)
                    case(54) ! '_13/134'
                        sender_y(1)   = (Bs+1)/2
                        sender_y(2)   = Bs+g
                        receiver_y(1) = 1
                        receiver_y(2) = Bs+g
                        finer_y(1)    = g+1
                        finer_y(2)    = Bs+2*g

                    case(53) ! '_13/123'
                        sender_y(1)   = g+1
                        sender_y(2)   = (Bs+1)/2+g+g
                        receiver_y(1) = g+1
                        receiver_y(2) = Bs+2*g
                        finer_y(1)    = 1
                        finer_y(2)    = Bs+g

                end select

                ! copy data
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    data_face = hvy_block( sender_x(1):sender_x(2), sender_y(1):sender_y(2), sender_z(1):sender_z(2), dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = data_face_fine(finer_x(1):finer_x(2), finer_y(1):finer_y(2), finer_z(1):finer_z(2))
                end do

            elseif ( level_diff == 1 ) then
                ! x,y,z
                sender_z(1)   = Bs-g
                sender_z(2)   = Bs+g-2+2*rmv_redundant
                sender_x(1)   = Bs-g
                sender_x(2)   = Bs+g-2+2*rmv_redundant
                sender_y(1)   = g+1
                sender_y(2)   = Bs+g

                receiver_z(1) = 1
                receiver_z(2) = g+rmv_redundant
                receiver_x(1) = 1
                receiver_x(2) = g+rmv_redundant
                ! y
                select case(neighborhood)
                    case(54) ! '_13/134'
                        receiver_y(1) = g+(Bs+1)/2
                        receiver_y(2) = Bs+g

                    case(53) ! '_13/123'
                        receiver_y(1) = g+1
                        receiver_y(2) = g+(Bs+1)/2

                end select

                ! copy data
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = hvy_block( sender_x(1):sender_x(2):2, sender_y(1):sender_y(2):2, sender_z(1):sender_z(2):2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if
!!        ! '_13/123'
!!        case(53)
!!            if ( level_diff == -1 ) then
!!                ! sender on lower level
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    ! data to interpolate, note: use data_face interpolation variable
!!                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
!!                    data_face = hvy_block( (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, dF, sender_id )
!!                    ! interpolate data
!!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!!                    ! copy data
!!                    !hvy_block( 1:g, g+1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, g+1:Bs+g, Bs:Bs-1+g)
!!                    !hvy_block( 1:g, g+1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, 1:Bs, Bs:Bs-1+g)
!!                    hvy_block( 1:g, g+1:Bs+2*g, 1:g, dF, receiver_id ) = data_face_fine(Bs+g:Bs-1+2*g, 1:Bs+g, Bs+g:Bs-1+2*g)
!!                end do
!!
!!            elseif ( level_diff == 1 ) then
!!                ! sender on higher level
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    !hvy_block( 1:g, g+1:g+(Bs+1)/2, 1:g, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2:2, g+1:Bs+g:2, Bs-g:Bs+g-2:2, dF, sender_id )
!!                    hvy_block( 1:g+rmv_redundant, g+1:g+(Bs+1)/2, 1:g+rmv_redundant, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, dF, sender_id )
!!                end do
!!
!!            else
!!                ! error case
!!                write(*,'(80("_"))')
!!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!!                stop
!!            end if
!!        ! '_13/134'
!!        case(54)
!!            if ( level_diff == -1 ) then
!!                ! sender on lower level
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    ! data to interpolate, note: use data_face interpolation variable
!!                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
!!                    data_face = hvy_block( (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, dF, sender_id )
!!                    ! interpolate data
!!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!!                    ! copy data
!!                    !hvy_block( 1:g, g+1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, g+1:Bs+g, Bs:Bs-1+g)
!!                    hvy_block( 1:g, 1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(Bs+g:Bs-1+2*g, g+1:Bs+2*g, Bs+g:Bs-1+2*g)
!!                end do
!!
!!            elseif ( level_diff == 1 ) then
!!                ! sender on higher level
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    !hvy_block( 1:g, g+(Bs+1)/2:Bs+g, 1:g, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2:2, g+1:Bs+g:2, Bs-g:Bs+g-2:2, dF, sender_id )
!!                    hvy_block( 1:g+rmv_redundant, g+(Bs+1)/2:Bs+g, 1:g+rmv_redundant, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, dF, sender_id )
!!                end do
!!
!!            else
!!                ! error case
!!                write(*,'(80("_"))')
!!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!!                stop
!!            end if
        case(55,56)
            if ( level_diff == -1 ) then
                ! z,y
                sender_z(1)   = (Bs+1)/2
                sender_z(2)   = Bs+g
                receiver_z(1) = 1
                receiver_z(2) = g
                finer_z(1)    = Bs+g
                finer_z(2)    = Bs-1+2*g

                sender_y(1)   = (Bs+1)/2
                sender_y(2)   = Bs+g
                receiver_y(1) = 1
                receiver_y(2) = g
                finer_y(1)    = Bs+g
                finer_y(2)    = Bs-1+2*g
                ! x,y
                select case(neighborhood)
                    case(55) ! '_14/134'
                        sender_x(1)   = (Bs+1)/2
                        sender_x(2)   = Bs+g
                        receiver_x(1) = 1
                        receiver_x(2) = Bs+g
                        finer_x(1)    = g+1
                        finer_x(2)    = Bs+2*g

                    case(56) ! '_14/145'
                        sender_x(1)   = g+1
                        sender_x(2)   = (Bs+1)/2+g+g
                        receiver_x(1) = g+1
                        receiver_x(2) = Bs+2*g
                        finer_x(1)    = 1
                        finer_x(2)    = Bs+g

                end select

                ! copy data
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    data_face = hvy_block( sender_x(1):sender_x(2), sender_y(1):sender_y(2), sender_z(1):sender_z(2), dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = data_face_fine(finer_x(1):finer_x(2), finer_y(1):finer_y(2), finer_z(1):finer_z(2))
                end do

            elseif ( level_diff == 1 ) then
                ! x,y,z
                sender_z(1)   = Bs-g
                sender_z(2)   = Bs+g-2+2*rmv_redundant
                sender_y(1)   = Bs-g
                sender_y(2)   = Bs+g-2+2*rmv_redundant
                sender_x(1)   = g+1
                sender_x(2)   = Bs+g

                receiver_z(1) = 1
                receiver_z(2) = g+rmv_redundant
                receiver_y(1) = 1
                receiver_y(2) = g+rmv_redundant
                ! y
                select case(neighborhood)
                    case(55) ! '_14/134'
                        receiver_x(1) = g+(Bs+1)/2
                        receiver_x(2) = Bs+g

                    case(56) ! '_14/145'
                        receiver_x(1) = g+1
                        receiver_x(2) = g+(Bs+1)/2

                end select

                ! copy data
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = hvy_block( sender_x(1):sender_x(2):2, sender_y(1):sender_y(2):2, sender_z(1):sender_z(2):2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if
!!        ! '_14/134'
!!        case(55)
!!            if ( level_diff == -1 ) then
!!                ! sender on lower level
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    ! data to interpolate, note: use data_face interpolation variable
!!                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
!!                    data_face = hvy_block( (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, dF, sender_id )
!!                    ! interpolate data
!!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!!                    ! copy data
!!                    !hvy_block( g+1:Bs+g, 1:g, 1:g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, Bs:Bs-1+g, Bs:Bs-1+g)
!!                    hvy_block( g+1:Bs+g, 1:g, 1:g, dF, receiver_id ) = data_face_fine(g+g+1:Bs+2*g, Bs+g:Bs-1+2*g, Bs+g:Bs-1+2*g)
!!                end do
!!
!!            elseif ( level_diff == 1 ) then
!!                ! sender on higher level
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    !hvy_block( g+(Bs+1)/2:Bs+g, 1:g, 1:g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2:2, Bs-g:Bs+g-2:2, dF, sender_id )
!!                    hvy_block( g+(Bs+1)/2:Bs+g, 1:g+rmv_redundant, 1:g+rmv_redundant, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, Bs-g:Bs+g-2+2*rmv_redundant:2, dF, sender_id )
!!                end do
!!
!!            else
!!                ! error case
!!                write(*,'(80("_"))')
!!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!!                stop
!!            end if
!!
!!        ! '_14/145'
!!        case(56)
!!            if ( level_diff == -1 ) then
!!                ! sender on lower level
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    ! data to interpolate, note: use data_face interpolation variable
!!                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
!!                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, dF, sender_id )
!!                    ! interpolate data
!!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!!                    ! copy data
!!                    !hvy_block( g+1:Bs+g, 1:g, 1:g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, Bs:Bs-1+g, Bs:Bs-1+g)
!!                    !hvy_block( g+1:Bs+g, 1:g, 1:g, dF, receiver_id ) = data_face_fine(1:Bs, Bs:Bs-1+g, Bs:Bs-1+g)
!!                    hvy_block( g+1:Bs+g, 1:g, 1:g, dF, receiver_id ) = data_face_fine(1:Bs, Bs+g:Bs-1+2*g, Bs+g:Bs-1+2*g)
!!                end do
!!
!!            elseif ( level_diff == 1 ) then
!!                ! sender on higher level
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    !hvy_block( g+1:g+(Bs+1)/2, 1:g, 1:g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2:2, Bs-g:Bs+g-2:2, dF, sender_id )
!!                    hvy_block( g+1:g+(Bs+1)/2, 1:g+rmv_redundant, 1:g+rmv_redundant, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, Bs-g:Bs+g-2+2*rmv_redundant:2, dF, sender_id )
!!                end do
!!
!!            else
!!                ! error case
!!                write(*,'(80("_"))')
!!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!!                stop
!!            end if
        case(57,58)
            if ( level_diff == -1 ) then
                ! z,x
                sender_z(1)   = (Bs+1)/2
                sender_z(2)   = Bs+g
                receiver_z(1) = 1
                receiver_z(2) = g
                finer_z(1)    = Bs+g
                finer_z(2)    = Bs-1+2*g

                sender_x(1)   = g+1
                sender_x(2)   = (Bs+1)/2+g+g
                receiver_x(1) = Bs+g+1
                receiver_x(2) = Bs+2*g
                finer_x(1)    = 2
                finer_x(2)    = g+1
                ! x,y
                select case(neighborhood)
                    case(57) ! '_15/145'
                        sender_y(1)   = (Bs+1)/2
                        sender_y(2)   = Bs+g
                        receiver_y(1) = 1
                        receiver_y(2) = Bs+g
                        finer_y(1)    = g+1
                        finer_y(2)    = Bs+2*g

                    case(58) ! '_15/152''
                        sender_y(1)   = g+1
                        sender_y(2)   = (Bs+1)/2+g+g
                        receiver_y(1) = g+1
                        receiver_y(2) = Bs+2*g
                        finer_y(1)    = 1
                        finer_y(2)    = Bs+g

                end select

                ! copy data
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    data_face = hvy_block( sender_x(1):sender_x(2), sender_y(1):sender_y(2), sender_z(1):sender_z(2), dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = data_face_fine(finer_x(1):finer_x(2), finer_y(1):finer_y(2), finer_z(1):finer_z(2))
                end do

            elseif ( level_diff == 1 ) then
                ! x,y,z
                sender_z(1)   = Bs-g
                sender_z(2)   = Bs+g-2+2*rmv_redundant
                sender_x(1)   = g+3-2*rmv_redundant
                sender_x(2)   = 3*g+1
                sender_y(1)   = g+1
                sender_y(2)   = Bs+g

                receiver_z(1) = 1
                receiver_z(2) = g+rmv_redundant
                receiver_x(1) = Bs+g+1-rmv_redundant
                receiver_x(2) = Bs+g+g
                ! y
                select case(neighborhood)
                    case(57) ! '_15/145'
                        receiver_y(1) = g+(Bs+1)/2
                        receiver_y(2) = Bs+g

                    case(58) ! '_15/152'
                        receiver_y(1) = g+1
                        receiver_y(2) = g+(Bs+1)/2

                end select

                ! copy data
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = hvy_block( sender_x(1):sender_x(2):2, sender_y(1):sender_y(2):2, sender_z(1):sender_z(2):2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if
!!        ! '_15/145'
!!        case(57)
!!            if ( level_diff == -1 ) then
!!                ! sender on lower level
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    ! data to interpolate, note: use data_face interpolation variable
!!                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
!!                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, dF, sender_id )
!!                    ! interpolate data
!!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!!                    ! copy data
!!                    !hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(2:g+1, g+1:Bs+g, Bs:Bs-1+g)
!!                    hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(2:g+1, g+g+1:Bs+2*g, Bs+g:Bs-1+2*g)
!!                end do
!!
!!            elseif ( level_diff == 1 ) then
!!                ! sender on higher level
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    !hvy_block( Bs+g+1:Bs+g+g, g+(Bs+1)/2:Bs+g, 1:g, dF, receiver_id ) = hvy_block( g+3:3*g+1:2, g+1:Bs+g:2, Bs-g:Bs+g-2:2, dF, sender_id )
!!                    hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, g+(Bs+1)/2:Bs+g, 1:g+rmv_redundant, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, dF, sender_id )
!!                end do
!!
!!            else
!!                ! error case
!!                write(*,'(80("_"))')
!!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!!                stop
!!            end if
!!
!!        ! '_15/152'
!!        case(58)
!!            if ( level_diff == -1 ) then
!!                ! sender on lower level
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    ! data to interpolate, note: use data_face interpolation variable
!!                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
!!                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, dF, sender_id )
!!                    ! interpolate data
!!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!!                    ! copy data
!!                    !hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(2:g+1, g+1:Bs+g, Bs:Bs-1+g)
!!                    !hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(2:g+1, 1:Bs, Bs:Bs-1+g)
!!                    hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(2:g+1, 1:Bs, Bs+g:Bs-1+2*g)
!!                end do
!!
!!            elseif ( level_diff == 1 ) then
!!                ! sender on higher level
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    !hvy_block( Bs+g+1:Bs+g+g, g+1:g+(Bs+1)/2, 1:g, dF, receiver_id ) = hvy_block( g+3:3*g+1:2, g+1:Bs+g:2, Bs-g:Bs+g-2:2, dF, sender_id )
!!                    hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, g+1:g+(Bs+1)/2, 1:g+rmv_redundant, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, dF, sender_id )
!!                end do
!!
!!            else
!!                ! error case
!!                write(*,'(80("_"))')
!!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!!                stop
!!            end if

          case(59,60)
            if ( level_diff == -1 ) then
                ! z,y
                sender_z(1)   = g+1
                sender_z(2)   = (Bs+1)/2+g+g
                receiver_z(1) = Bs+g+1
                receiver_z(2) = Bs+2*g
                finer_z(1)    = 2
                finer_z(2)    = g+1

                sender_y(1)   = g+1
                sender_y(2)   = (Bs+1)/2+g+g
                receiver_y(1) = Bs+g+1
                receiver_y(2) = Bs+2*g
                finer_y(1)    = 2
                finer_y(2)    = g+1
                ! x,y
                select case(neighborhood)
                    case(59) ! '_62/623'
                        sender_x(1)   = (Bs+1)/2
                        sender_x(2)   = Bs+g
                        receiver_x(1) = 1
                        receiver_x(2) = Bs+g
                        finer_x(1)    = g+1
                        finer_x(2)    = Bs+2*g

                    case(60) ! '_62/652'
                        sender_x(1)   = g+1
                        sender_x(2)   = (Bs+1)/2+g+g
                        receiver_x(1) = g+1
                        receiver_x(2) = Bs+2*g
                        finer_x(1)    = 1
                        finer_x(2)    = Bs+g

                end select

                ! copy data
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    data_face = hvy_block( sender_x(1):sender_x(2), sender_y(1):sender_y(2), sender_z(1):sender_z(2), dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = data_face_fine(finer_x(1):finer_x(2), finer_y(1):finer_y(2), finer_z(1):finer_z(2))
                end do

            elseif ( level_diff == 1 ) then
                ! x,y,z
                sender_z(1)   = g+3-2*rmv_redundant
                sender_z(2)   = 3*g+1
                sender_y(1)   = g+3-2*rmv_redundant
                sender_y(2)   = 3*g+1
                sender_x(1)   = g+1
                sender_x(2)   = Bs+g

                receiver_z(1) = Bs+g+1-rmv_redundant
                receiver_z(2) = Bs+g+g
                receiver_y(1) = Bs+g+1-rmv_redundant
                receiver_y(2) = Bs+g+g
                ! y
                select case(neighborhood)
                    case(59) ! '_62/623'
                        receiver_x(1) = g+(Bs+1)/2
                        receiver_x(2) = Bs+g

                    case(60) ! '_62/652'
                        receiver_x(1) = g+1
                        receiver_x(2) = g+(Bs+1)/2

                end select

                ! copy data
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = hvy_block( sender_x(1):sender_x(2):2, sender_y(1):sender_y(2):2, sender_z(1):sender_z(2):2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

!!         ! '_62/623'
!!        case(59)
!!            if ( level_diff == -1 ) then
!!                ! sender on lower level
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    ! data to interpolate, note: use data_face interpolation variable
!!                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
!!                    data_face = hvy_block( (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!!                    ! interpolate data
!!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!!                    ! copy data
!!                    !hvy_block( g+1:Bs+g, Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, 2:g+1, 2:g+1)
!!                    hvy_block( g+1:Bs+g, Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(g+g+1:Bs+2*g, 2:g+1, 2:g+1)
!!                end do
!!
!!            elseif ( level_diff == 1 ) then
!!                ! sender on higher level
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    !hvy_block( g+(Bs+1)/2:Bs+g, Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3:3*g+1:2, g+3:3*g+1:2, dF, sender_id )
!!                    hvy_block( g+(Bs+1)/2:Bs+g, Bs+g+1-rmv_redundant:Bs+g+g, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, g+3-2*rmv_redundant:3*g+1:2, dF, sender_id )
!!                end do
!!
!!            else
!!                ! error case
!!                write(*,'(80("_"))')
!!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!!                stop
!!            end if
!!
!!        ! '_62/652'
!!        case(60)
!!            if ( level_diff == -1 ) then
!!                ! sender on lower level
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    ! data to interpolate, note: use data_face interpolation variable
!!                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
!!                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!!                    ! interpolate data
!!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!!                    ! copy data
!!                    !hvy_block( g+1:Bs+g, Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, 2:g+1, 2:g+1)
!!                    hvy_block( g+1:Bs+g, Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(1:Bs, 2:g+1, 2:g+1)
!!                end do
!!
!!            elseif ( level_diff == 1 ) then
!!                ! sender on higher level
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    !hvy_block( g+1:g+(Bs+1)/2, Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3:3*g+1:2, g+3:3*g+1:2, dF, sender_id )
!!                    hvy_block( g+1:g+(Bs+1)/2, Bs+g+1-rmv_redundant:Bs+g+g, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, g+3-2*rmv_redundant:3*g+1:2, dF, sender_id )
!!                end do
!!
!!            else
!!                ! error case
!!                write(*,'(80("_"))')
!!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!!                stop
!!            end if

         case(61,62)
            if ( level_diff == -1 ) then
                ! z,x
                sender_z(1)   = g+1
                sender_z(2)   = (Bs+1)/2+g+g
                receiver_z(1) = Bs+g+1
                receiver_z(2) = Bs+2*g
                finer_z(1)    = 2
                finer_z(2)    = g+1

                sender_x(1)   = (Bs+1)/2
                sender_x(2)   = Bs+g
                receiver_x(1) = 1
                receiver_x(2) = g
                finer_x(1)    = Bs+g
                finer_x(2)    = Bs-1+2*g
                ! x,y
                select case(neighborhood)
                    case(62) ! '_63/634'
                        sender_y(1)   = (Bs+1)/2
                        sender_y(2)   = Bs+g
                        receiver_y(1) = 1
                        receiver_y(2) = Bs+g
                        finer_y(1)    = g+1
                        finer_y(2)    = Bs+2*g

                    case(61) ! '_63/623'
                        sender_y(1)   = g+1
                        sender_y(2)   = (Bs+1)/2+g+g
                        receiver_y(1) = g+1
                        receiver_y(2) = Bs+2*g
                        finer_y(1)    = 1
                        finer_y(2)    = Bs+g

                end select

                ! copy data
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    data_face = hvy_block( sender_x(1):sender_x(2), sender_y(1):sender_y(2), sender_z(1):sender_z(2), dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = data_face_fine(finer_x(1):finer_x(2), finer_y(1):finer_y(2), finer_z(1):finer_z(2))
                end do

            elseif ( level_diff == 1 ) then
                ! x,y,z
                sender_z(1)   = g+3-2*rmv_redundant
                sender_z(2)   = 3*g+1
                sender_x(1)   = Bs-g
                sender_x(2)   = Bs+g-2+2*rmv_redundant
                sender_y(1)   = g+1
                sender_y(2)   = Bs+g

                receiver_z(1) = Bs+g+1-rmv_redundant
                receiver_z(2) = Bs+g+g
                receiver_x(1) = 1
                receiver_x(2) = g+rmv_redundant
                ! y
                select case(neighborhood)
                    case(62) ! '_63/634'
                        receiver_y(1) = g+(Bs+1)/2
                        receiver_y(2) = Bs+g

                    case(61) ! '_63/623'
                        receiver_y(1) = g+1
                        receiver_y(2) = g+(Bs+1)/2

                end select

                ! copy data
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = hvy_block( sender_x(1):sender_x(2):2, sender_y(1):sender_y(2):2, sender_z(1):sender_z(2):2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

!        ! '_63/623'
!        case(61)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate, note: use data_face interpolation variable
!                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
!                    !data_face = hvy_block( (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!                    data_face = hvy_block( (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( 1:g, g+1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, g+1:Bs+g, 2:g+1)
!                    !hvy_block( 1:g, g+1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, 1:Bs, 2:g+1)
!                    hvy_block( 1:g, 1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(Bs+g:Bs-1+2*g, g+1:Bs+2*g, 2:g+1)
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( 1:g, g+1:g+(Bs+1)/2, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2:2, g+1:Bs+g:2, g+3:3*g+1:2, dF, sender_id )
!                    !hvy_block( 1:g+rmv_redundant, g+1:g+(Bs+1)/2, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, dF, sender_id )
!                    hvy_block( 1:g+rmv_redundant, g+(Bs+1)/2:Bs+g, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if
!
!        ! '_63/634'
!        case(62)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate, note: use data_face interpolation variable
!                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
!                    !data_face = hvy_block( (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!                    data_face = hvy_block( (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( 1:g, g+1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, g+1:Bs+g, 2:g+1)
!                    !hvy_block( 1:g, g+1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(Bs+g:Bs-1+2*g, g+g+1:Bs+2*g, 2:g+1)
!                    hvy_block( 1:g, g+1:Bs+2*g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(Bs+g:Bs-1+2*g, 1:Bs+g, 2:g+1)
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( 1:g, g+(Bs+1)/2:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2:2, g+1:Bs+g:2, g+3:3*g+1:2, dF, sender_id )
!                    !hvy_block( 1:g+rmv_redundant, g+(Bs+1)/2:Bs+g, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, dF, sender_id )
!                    hvy_block( 1:g+rmv_redundant, g+1:g+(Bs+1)/2, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if

          case(63,64)
            if ( level_diff == -1 ) then
                ! z,y
                sender_z(1)   = g+1
                sender_z(2)   = (Bs+1)/2+g+g
                receiver_z(1) = Bs+g+1
                receiver_z(2) = Bs+2*g
                finer_z(1)    = 2
                finer_z(2)    = g+1

                sender_y(1)   = (Bs+1)/2
                sender_y(2)   = Bs+g
                receiver_y(1) = 1
                receiver_y(2) = g
                finer_y(1)    = Bs+g
                finer_y(2)    = Bs-1+2*g
                ! x,y
                select case(neighborhood)
                    case(63) ! '_64/634'
                        sender_x(1)   = (Bs+1)/2
                        sender_x(2)   = Bs+g
                        receiver_x(1) = 1
                        receiver_x(2) = Bs+g
                        finer_x(1)    = g+1
                        finer_x(2)    = Bs+2*g

                    case(64) ! '_64/645'
                        sender_x(1)   = g+1
                        sender_x(2)   = (Bs+1)/2+g+g
                        receiver_x(1) = g+1
                        receiver_x(2) = Bs+2*g
                        finer_x(1)    = 1
                        finer_x(2)    = Bs+g

                end select

                ! copy data
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    data_face = hvy_block( sender_x(1):sender_x(2), sender_y(1):sender_y(2), sender_z(1):sender_z(2), dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = data_face_fine(finer_x(1):finer_x(2), finer_y(1):finer_y(2), finer_z(1):finer_z(2))
                end do

            elseif ( level_diff == 1 ) then
                ! x,y,z
                sender_z(1)   = g+3-2*rmv_redundant
                sender_z(2)   = 3*g+1
                sender_y(1)   = Bs-g
                sender_y(2)   = Bs+g-2+2*rmv_redundant
                sender_x(1)   = g+1
                sender_x(2)   = Bs+g

                receiver_z(1) = Bs+g+1-rmv_redundant
                receiver_z(2) = Bs+g+g
                receiver_y(1) = 1
                receiver_y(2) = g+rmv_redundant
                ! y
                select case(neighborhood)
                    case(63) ! '_64/634'
                        receiver_x(1) = g+(Bs+1)/2
                        receiver_x(2) = Bs+g

                    case(64) ! '_64/645'
                        receiver_x(1) = g+1
                        receiver_x(2) = g+(Bs+1)/2

                end select

                ! copy data
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = hvy_block( sender_x(1):sender_x(2):2, sender_y(1):sender_y(2):2, sender_z(1):sender_z(2):2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

!!        ! '_64/634'
!!        case(63)
!!            if ( level_diff == -1 ) then
!!                ! sender on lower level
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    ! data to interpolate, note: use data_face interpolation variable
!!                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
!!                    data_face = hvy_block( (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!!                    ! interpolate data
!!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!!                    ! copy data
!!                    !hvy_block( g+1:Bs+g, 1:g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, Bs:Bs-1+g, 2:g+1)
!!                    hvy_block( g+1:Bs+g, 1:g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(g+g+1:Bs+2*g, Bs+g:Bs-1+2*g, 2:g+1)
!!                end do
!!
!!            elseif ( level_diff == 1 ) then
!!                ! sender on higher level
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    !hvy_block( g+(Bs+1)/2:Bs+g, 1:g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2:2, g+3:3*g+1:2, dF, sender_id )
!!                    hvy_block( g+(Bs+1)/2:Bs+g, 1:g+rmv_redundant, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, g+3-2*rmv_redundant:3*g+1:2, dF, sender_id )
!!                end do
!!
!!            else
!!                ! error case
!!                write(*,'(80("_"))')
!!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!!                stop
!!            end if
!!
!!        ! '_64/645'
!!        case(64)
!!            if ( level_diff == -1 ) then
!!                ! sender on lower level
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    ! data to interpolate, note: use data_face interpolation variable
!!                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
!!                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!!                    ! interpolate data
!!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!!                    ! copy data
!!                    !hvy_block( g+1:Bs+g, 1:g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, Bs:Bs-1+g, 2:g+1)
!!                    !hvy_block( g+1:Bs+g, 1:g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(1:Bs, Bs:Bs-1+g, 2:g+1)
!!                    hvy_block( g+1:Bs+g, 1:g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(1:Bs, Bs+g:Bs-1+2*g, 2:g+1)
!!                end do
!!
!!            elseif ( level_diff == 1 ) then
!!                ! sender on higher level
!!                ! loop over all datafields
!!                do dF = 1, params%number_data_fields
!!                    !hvy_block( g+1:g+(Bs+1)/2, 1:g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2:2, g+3:3*g+1:2, dF, sender_id )
!!                    hvy_block( g+1:g+(Bs+1)/2, 1:g+rmv_redundant, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, g+3-2*rmv_redundant:3*g+1:2, dF, sender_id )
!!                end do
!!
!!            else
!!                ! error case
!!                write(*,'(80("_"))')
!!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!!                stop
!!            end if

        case(65,66)
            if ( level_diff == -1 ) then
                ! z,x
                sender_z(1)   = g+1
                sender_z(2)   = (Bs+1)/2+g+g
                receiver_z(1) = Bs+g+1
                receiver_z(2) = Bs+2*g
                finer_z(1)    = 2
                finer_z(2)    = g+1

                sender_x(1)   = g+1
                sender_x(2)   = (Bs+1)/2+g+g
                receiver_x(1) = Bs+g+1
                receiver_x(2) = Bs+2*g
                finer_x(1)    = 2
                finer_x(2)    = g+1
                ! x,y
                select case(neighborhood)
                    case(65) ! '_65/645'
                        sender_y(1)   = (Bs+1)/2
                        sender_y(2)   = Bs+g
                        receiver_y(1) = 1
                        receiver_y(2) = Bs+g
                        finer_y(1)    = g+1
                        finer_y(2)    = Bs+2*g

                    case(66) ! '_65/652''
                        sender_y(1)   = g+1
                        sender_y(2)   = (Bs+1)/2+g+g
                        receiver_y(1) = g+1
                        receiver_y(2) = Bs+2*g
                        finer_y(1)    = 1
                        finer_y(2)    = Bs+g

                end select

                ! copy data
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    data_face = hvy_block( sender_x(1):sender_x(2), sender_y(1):sender_y(2), sender_z(1):sender_z(2), dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = data_face_fine(finer_x(1):finer_x(2), finer_y(1):finer_y(2), finer_z(1):finer_z(2))
                end do

            elseif ( level_diff == 1 ) then
                ! x,y,z
                sender_z(1)   = g+3-2*rmv_redundant
                sender_z(2)   = 3*g+1
                sender_x(1)   = g+3-2*rmv_redundant
                sender_x(2)   = 3*g+1
                sender_y(1)   = g+1
                sender_y(2)   = Bs+g

                receiver_z(1) = Bs+g+1-rmv_redundant
                receiver_z(2) = Bs+g+g
                receiver_x(1) = Bs+g+1-rmv_redundant
                receiver_x(2) = Bs+g+g
                ! y
                select case(neighborhood)
                    case(65) ! '_65/645'
                        receiver_y(1) = g+(Bs+1)/2
                        receiver_y(2) = Bs+g

                    case(66) ! '_65/652'
                        receiver_y(1) = g+1
                        receiver_y(2) = g+(Bs+1)/2

                end select

                ! copy data
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = hvy_block( sender_x(1):sender_x(2):2, sender_y(1):sender_y(2):2, sender_z(1):sender_z(2):2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

!        ! '_65/645'
!        case(65)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate, note: use data_face interpolation variable
!                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
!                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(2:g+1, g+1:Bs+g, 2:g+1)
!                    !hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(2:g+1, g+g+1:Bs+2*g, 2:g+1)
!                    hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+2*g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(2:g+1, 1:Bs+g, 2:g+1)
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( Bs+g+1:Bs+g+g, g+(Bs+1)/2:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+3:3*g+1:2, g+1:Bs+g:2, g+3:3*g+1:2, dF, sender_id )
!                    !hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, g+(Bs+1)/2:Bs+g, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, dF, sender_id )
!                    hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, g+1:g+(Bs+1)/2, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if
!
!        ! '_65/652'
!        case(66)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate, note: use data_face interpolation variable
!                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
!                    !data_face = hvy_block( g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(2:g+1, g+1:Bs+g, 2:g+1)
!                    !hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(2:g+1, 1:Bs, 2:g+1)
!                    hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+2*g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(2:g+1, 1:Bs+g, 2:g+1)
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( Bs+g+1:Bs+g+g, g+1:g+(Bs+1)/2, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+3:3*g+1:2, g+1:Bs+g:2, g+3:3*g+1:2, dF, sender_id )
!                    !hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, g+1:g+(Bs+1)/2, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, dF, sender_id )
!                    hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, g+(Bs+1)/2:Bs+g, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if

         case(67,68)
            if ( level_diff == -1 ) then
                ! y,x
                sender_x(1)   = (Bs+1)/2
                sender_x(2)   = Bs+g
                receiver_x(1) = 1
                receiver_x(2) = g
                finer_x(1)    = Bs+g
                finer_x(2)    = Bs-1+2*g

                sender_y(1)   = g+1
                sender_y(2)   = (Bs+1)/2+g+g
                receiver_y(1) = Bs+g+1
                receiver_y(2) = Bs+2*g
                finer_y(1)    = 2
                finer_y(2)    = g+1

                ! x,y
                select case(neighborhood)
                    case(67) ! '_23/123'
                        sender_z(1)   = (Bs+1)/2
                        sender_z(2)   = Bs+g
                        receiver_z(1) = 1
                        receiver_z(2) = Bs+g
                        finer_z(1)    = g+1
                        finer_z(2)    = Bs+2*g

                    case(68) ! '_23/236''
                        sender_z(1)   = g+1
                        sender_z(2)   = (Bs+1)/2+g+g
                        receiver_z(1) = g+1
                        receiver_z(2) = Bs+2*g
                        finer_z(1)    = 1
                        finer_z(2)    = Bs+g

                end select
               
                ! copy data
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    data_face = hvy_block( sender_x(1):sender_x(2), sender_y(1):sender_y(2), sender_z(1):sender_z(2), dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = data_face_fine(finer_x(1):finer_x(2), finer_y(1):finer_y(2), finer_z(1):finer_z(2))
                end do

            elseif ( level_diff == 1 ) then
                ! x,y,z
                sender_x(1)   = Bs-g
                sender_x(2)   = Bs+g-2+2*rmv_redundant
                sender_y(1)   = g+3-2*rmv_redundant
                sender_y(2)   = 3*g+1
                sender_z(1)   = g+1
                sender_z(2)   = Bs+g

                receiver_x(1) = 1
                receiver_x(2) = g+rmv_redundant
                receiver_y(1) = Bs+g+1-rmv_redundant
                receiver_y(2) = Bs+g+g
                ! y
                select case(neighborhood)
                    case(67) ! '_23/123'
                        receiver_z(1) = g+(Bs+1)/2
                        receiver_z(2) = Bs+g

                    case(68) ! '_23/236'
                        receiver_z(1) = g+1
                        receiver_z(2) = g+(Bs+1)/2

                end select

                ! copy data
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = hvy_block( sender_x(1):sender_x(2):2, sender_y(1):sender_y(2):2, sender_z(1):sender_z(2):2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

!        ! '_23/123'
!        case(67)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate, note: use data_face interpolation variable
!                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
!                    data_face = hvy_block( (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( 1:g, Bs+g+1:Bs+g+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, 2:g+1, g+1:Bs+g)
!                    hvy_block( 1:g, Bs+g+1:Bs+g+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(Bs+g:Bs-1+2*g, 2:g+1, g+g+1:Bs+2*g)
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( 1:g, Bs+g+1:Bs+g+g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2:2, g+3:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
!                    hvy_block( 1:g+rmv_redundant, Bs+g+1-rmv_redundant:Bs+g+g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2+2*rmv_redundant:2, g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if
!
!        ! '_23/623'
!        case(68)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate, note: use data_face interpolation variable
!                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
!                    data_face = hvy_block( (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( 1:g, Bs+g+1:Bs+g+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, 2:g+1, g+1:Bs+g)
!                    !hvy_block( 1:g, Bs+g+1:Bs+g+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, 2:g+1, 1:Bs)
!                    hvy_block( 1:g, Bs+g+1:Bs+g+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(Bs+g:Bs-1+2*g, 2:g+1, 1:Bs)
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( 1:g, Bs+g+1:Bs+g+g, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2:2, g+3:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
!                    hvy_block( 1:g+rmv_redundant, Bs+g+1-rmv_redundant:Bs+g+g, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2+2*rmv_redundant:2, g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if

         case(69,70)
            if ( level_diff == -1 ) then
                ! y,x
                sender_x(1)   = g+1
                sender_x(2)   = (Bs+1)/2+g+g
                receiver_x(1) = Bs+g+1
                receiver_x(2) = Bs+2*g
                finer_x(1)    = 2
                finer_x(2)    = g+1

                sender_y(1)   = g+1
                sender_y(2)   = (Bs+1)/2+g+g
                receiver_y(1) = Bs+g+1
                receiver_y(2) = Bs+2*g
                finer_y(1)    = 2
                finer_y(2)    = g+1

                ! x,y
                select case(neighborhood)
                    case(69) ! '_25/152'
                        sender_z(1)   = (Bs+1)/2
                        sender_z(2)   = Bs+g
                        receiver_z(1) = 1
                        receiver_z(2) = Bs+g
                        finer_z(1)    = g+1
                        finer_z(2)    = Bs+2*g

                    case(70) ! '_25/652''
                        sender_z(1)   = g+1
                        sender_z(2)   = (Bs+1)/2+g+g
                        receiver_z(1) = g+1
                        receiver_z(2) = Bs+2*g
                        finer_z(1)    = 1
                        finer_z(2)    = Bs+g

                end select
               
                ! copy data
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    data_face = hvy_block( sender_x(1):sender_x(2), sender_y(1):sender_y(2), sender_z(1):sender_z(2), dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = data_face_fine(finer_x(1):finer_x(2), finer_y(1):finer_y(2), finer_z(1):finer_z(2))
                end do

            elseif ( level_diff == 1 ) then
                ! x,y,z
                sender_x(1)   = g+3-2*rmv_redundant
                sender_x(2)   = 3*g+1
                sender_y(1)   = g+3-2*rmv_redundant
                sender_y(2)   = 3*g+1
                sender_z(1)   = g+1
                sender_z(2)   = Bs+g

                receiver_x(1) = Bs+g+1-rmv_redundant
                receiver_x(2) = Bs+g+g
                receiver_y(1) = Bs+g+1-rmv_redundant
                receiver_y(2) = Bs+g+g
                ! y
                select case(neighborhood)
                    case(69) ! '_25/152'
                        receiver_z(1) = g+(Bs+1)/2
                        receiver_z(2) = Bs+g

                    case(70) ! '_25/652'
                        receiver_z(1) = g+1
                        receiver_z(2) = g+(Bs+1)/2

                end select

                ! copy data
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = hvy_block( sender_x(1):sender_x(2):2, sender_y(1):sender_y(2):2, sender_z(1):sender_z(2):2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

!        ! '_25/152'
!        case(69)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate, note: use data_face interpolation variable
!                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
!                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(2:g+1, 2:g+1, g+1:Bs+g)
!                    hvy_block( Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(2:g+1, 2:g+1, g+g+1:Bs+2*g)
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+3:3*g+1:2, g+3:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
!                    hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, Bs+g+1-rmv_redundant:Bs+g+g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:3*g+1:2, g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if
!
!        ! '_25/652'
!        case(70)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate, note: use data_face interpolation variable
!                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
!                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(2:g+1, 2:g+1, g+1:Bs+g)
!                    hvy_block( Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(2:g+1, 2:g+1, 1:Bs)
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+3:3*g+1:2, g+3:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
!                    hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, Bs+g+1-rmv_redundant:Bs+g+g, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:3*g+1:2, g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if

         case(71,72)
            if ( level_diff == -1 ) then
                ! y,x
                sender_x(1)   = (Bs+1)/2
                sender_x(2)   = Bs+g
                receiver_x(1) = 1
                receiver_x(2) = g
                finer_x(1)    = Bs+g
                finer_x(2)    = Bs-1+2*g

                sender_y(1)   = (Bs+1)/2
                sender_y(2)   = Bs+g
                receiver_y(1) = 1
                receiver_y(2) = g
                finer_y(1)    = Bs+g
                finer_y(2)    = Bs-1+2*g

                ! x,y
                select case(neighborhood)
                    case(71) ! '_43/134'
                        sender_z(1)   = (Bs+1)/2
                        sender_z(2)   = Bs+g
                        receiver_z(1) = 1
                        receiver_z(2) = Bs+g
                        finer_z(1)    = g+1
                        finer_z(2)    = Bs+2*g

                    case(72) ! '_43/634''
                        sender_z(1)   = g+1
                        sender_z(2)   = (Bs+1)/2+g+g
                        receiver_z(1) = g+1
                        receiver_z(2) = Bs+2*g
                        finer_z(1)    = 1
                        finer_z(2)    = Bs+g

                end select
               
                ! copy data
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    data_face = hvy_block( sender_x(1):sender_x(2), sender_y(1):sender_y(2), sender_z(1):sender_z(2), dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = data_face_fine(finer_x(1):finer_x(2), finer_y(1):finer_y(2), finer_z(1):finer_z(2))
                end do

            elseif ( level_diff == 1 ) then
                ! x,y,z
                sender_x(1)   = Bs-g
                sender_x(2)   = Bs+g-2+2*rmv_redundant
                sender_y(1)   = Bs-g
                sender_y(2)   = Bs+g-2+2*rmv_redundant
                sender_z(1)   = g+1
                sender_z(2)   = Bs+g

                receiver_x(1) = 1
                receiver_x(2) = g+rmv_redundant
                receiver_y(1) = 1
                receiver_y(2) = g+rmv_redundant
                ! y
                select case(neighborhood)
                    case(71) ! '_43/134'
                        receiver_z(1) = g+(Bs+1)/2
                        receiver_z(2) = Bs+g

                    case(72) ! '_43/634'
                        receiver_z(1) = g+1
                        receiver_z(2) = g+(Bs+1)/2

                end select

                ! copy data
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = hvy_block( sender_x(1):sender_x(2):2, sender_y(1):sender_y(2):2, sender_z(1):sender_z(2):2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

!        ! '_43/134'
!        case(71)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate, note: use data_face interpolation variable
!                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
!                    data_face = hvy_block( (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( 1:g, 1:g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, Bs:Bs-1+g, g+1:Bs+g)
!                    hvy_block( 1:g, 1:g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(Bs+g:Bs-1+2*g, Bs+g:Bs-1+2*g, g+g+1:Bs+2*g)
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( 1:g, 1:g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2:2, Bs-g:Bs+g-2:2, g+1:Bs+g:2, dF, sender_id )
!                    hvy_block( 1:g+rmv_redundant, 1:g+rmv_redundant, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2+2*rmv_redundant:2, Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if
!
!        ! '_43/634'
!        case(72)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate, note: use data_face interpolation variable
!                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
!                    data_face = hvy_block( (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( 1:g, 1:g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, Bs:Bs-1+g, g+1:Bs+g)
!                    !hvy_block( 1:g, 1:g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, Bs:Bs-1+g, 1:Bs)
!                    hvy_block( 1:g, 1:g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(Bs+g:Bs-1+2*g, Bs+g:Bs-1+2*g, 1:Bs)
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( 1:g, 1:g, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2:2, Bs-g:Bs+g-2:2, g+1:Bs+g:2, dF, sender_id )
!                    hvy_block( 1:g+rmv_redundant, 1:g+rmv_redundant, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2+2*rmv_redundant:2, Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if

         case(73,74)
            if ( level_diff == -1 ) then
                ! y,x
                sender_y(1)   = (Bs+1)/2
                sender_y(2)   = Bs+g
                receiver_y(1) = 1
                receiver_y(2) = g
                finer_y(1)    = Bs+g
                finer_y(2)    = Bs-1+2*g

                sender_x(1)   = g+1
                sender_x(2)   = (Bs+1)/2+g+g
                receiver_x(1) = Bs+g+1
                receiver_x(2) = Bs+2*g
                finer_x(1)    = 2
                finer_x(2)    = g+1

                ! x,y
                select case(neighborhood)
                    case(73) ! '_45/145'
                        sender_z(1)   = (Bs+1)/2
                        sender_z(2)   = Bs+g
                        receiver_z(1) = 1
                        receiver_z(2) = Bs+g
                        finer_z(1)    = g+1
                        finer_z(2)    = Bs+2*g

                    case(74) ! '_45/645'
                        sender_z(1)   = g+1
                        sender_z(2)   = (Bs+1)/2+g+g
                        receiver_z(1) = g+1
                        receiver_z(2) = Bs+2*g
                        finer_z(1)    = 1
                        finer_z(2)    = Bs+g

                end select
               
                ! copy data
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    data_face = hvy_block( sender_x(1):sender_x(2), sender_y(1):sender_y(2), sender_z(1):sender_z(2), dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = data_face_fine(finer_x(1):finer_x(2), finer_y(1):finer_y(2), finer_z(1):finer_z(2))
                end do

            elseif ( level_diff == 1 ) then
                ! x,y,z
                sender_y(1)   = Bs-g
                sender_y(2)   = Bs+g-2+2*rmv_redundant
                sender_x(1)   = g+3-2*rmv_redundant
                sender_x(2)   = 3*g+1
                sender_z(1)   = g+1
                sender_z(2)   = Bs+g

                receiver_y(1) = 1
                receiver_y(2) = g+rmv_redundant
                receiver_x(1) = Bs+g+1-rmv_redundant
                receiver_x(2) = Bs+g+g
                ! y
                select case(neighborhood)
                    case(73) ! '_45/145'
                        receiver_z(1) = g+(Bs+1)/2
                        receiver_z(2) = Bs+g

                    case(74) ! '_45/645'
                        receiver_z(1) = g+1
                        receiver_z(2) = g+(Bs+1)/2

                end select

                ! copy data
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    hvy_block( receiver_x(1):receiver_x(2), receiver_y(1):receiver_y(2), receiver_z(1):receiver_z(2), dF, receiver_id ) &
                    = hvy_block( sender_x(1):sender_x(2):2, sender_y(1):sender_y(2):2, sender_z(1):sender_z(2):2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

!        ! '_45/145'
!        case(73)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate, note: use data_face interpolation variable
!                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
!                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( Bs+g+1:Bs+g+g, 1:g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(2:g+1, Bs:Bs-1+g, g+1:Bs+g)
!                    hvy_block( Bs+g+1:Bs+g+g, 1:g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(2:g+1, Bs+g:Bs-1+2*g, g+g+1:Bs+2*g)
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( Bs+g+1:Bs+g+g, 1:g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+3:3*g+1:2, Bs-g:Bs+g-2:2, g+1:Bs+g:2, dF, sender_id )
!                    hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, 1:g+rmv_redundant, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:3*g+1:2, Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if
!
!        ! '_45/645'
!        case(74)
!            if ( level_diff == -1 ) then
!                ! sender on lower level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    ! data to interpolate, note: use data_face interpolation variable
!                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
!                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
!                    ! interpolate data
!                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
!                    ! copy data
!                    !hvy_block( Bs+g+1:Bs+g+g, 1:g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(2:g+1, Bs:Bs-1+g, g+1:Bs+g)
!                    !hvy_block( Bs+g+1:Bs+g+g, 1:g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(2:g+1, Bs:Bs-1+g, 1:Bs)
!                    hvy_block( Bs+g+1:Bs+g+g, 1:g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(2:g+1, Bs+g:Bs-1+2*g, 1:Bs)
!                end do
!
!            elseif ( level_diff == 1 ) then
!                ! sender on higher level
!                ! loop over all datafields
!                do dF = 1, params%number_data_fields
!                    !hvy_block( Bs+g+1:Bs+g+g, 1:g, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+3:3*g+1:2, Bs-g:Bs+g-2:2, g+1:Bs+g:2, dF, sender_id )
!                    hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, 1:g+rmv_redundant, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:3*g+1:2, Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, dF, sender_id )
!                end do
!
!            else
!                ! error case
!                write(*,'(80("_"))')
!                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
!                stop
!            end if

    end select

    ! clean up
    deallocate( data_corner  )
    deallocate( data_corner_fine  )
    deallocate( data_face  )
    deallocate( data_face_fine  )

end subroutine copy_ghost_nodes_3D
