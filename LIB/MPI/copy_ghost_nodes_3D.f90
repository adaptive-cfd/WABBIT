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

        ! '123/___'
        case(19)
            if ( level_diff == 0 ) then
                ! blocks on same level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    hvy_block( 1:g, Bs+g+1:Bs+g+g, 1:g, dF, receiver_id ) = hvy_block( Bs:Bs-1+g, g+2:g+1+g, Bs:Bs-1+g, dF, sender_id )
                end do

            elseif ( level_diff == -1 ) then
                ! sender one level down
                ! interpolate data
                do dF = 1, params%number_data_fields
                    ! data to refine
                    data_corner = hvy_block( Bs+1:Bs+g, g+1:g+g, Bs+1:Bs+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_corner , data_corner_fine, params%order_predictor)
                    ! data to synchronize
                    data_corner = data_corner_fine(g-1:2*g-2, 2:g+1, g-1:2*g-2)
                    ! write data
                    hvy_block( 1:g, Bs+g+1:Bs+g+g, 1:g, dF, receiver_id ) = data_corner
                end do

            elseif ( level_diff == 1) then
                ! sender one level up
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( 1:g, Bs+g+1:Bs+g+g, 1:g, dF, receiver_id ) = hvy_block( Bs-g:Bs-2+g:2, g+3:g+1+g+g:2, Bs-g:Bs-2+g:2, dF, sender_id )
                    hvy_block( 1:g+rmv_redundant, Bs+g+1-rmv_redundant:Bs+g+g, 1:g+rmv_redundant, dF, receiver_id ) = hvy_block( Bs-g:Bs-2+g+2*rmv_redundant:2, g+3-2*rmv_redundant:g+1+g+g:2, Bs-g:Bs-2+g+2*rmv_redundant:2, dF, sender_id )
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
                do dF = 1, params%number_data_fields
                    hvy_block( 1:g, 1:g, 1:g, dF, receiver_id ) = hvy_block( Bs:Bs-1+g, Bs:Bs-1+g, Bs:Bs-1+g, dF, sender_id )
                end do

            elseif ( level_diff == -1 ) then
                ! sender one level down
                ! interpolate data
                do dF = 1, params%number_data_fields
                    ! data to refine
                    data_corner = hvy_block( Bs+1:Bs+g, Bs+1:Bs+g, Bs+1:Bs+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_corner , data_corner_fine, params%order_predictor)
                    ! data to synchronize
                    data_corner = data_corner_fine(g-1:2*g-2, g-1:2*g-2, g-1:2*g-2)
                    ! write data
                    hvy_block( 1:g, 1:g, 1:g, dF, receiver_id ) = data_corner
                end do

            elseif ( level_diff == 1) then
                ! sender one level up
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( 1:g, 1:g, 1:g, dF, receiver_id ) = hvy_block( Bs-g:Bs-2+g:2, Bs-g:Bs-2+g:2, Bs-g:Bs-2+g:2, dF, sender_id )
                    hvy_block( 1:g+rmv_redundant, 1:g+rmv_redundant, 1:g+rmv_redundant, dF, receiver_id ) = hvy_block( Bs-g:Bs-2+g+2*rmv_redundant:2, Bs-g:Bs-2+g+2*rmv_redundant:2, Bs-g:Bs-2+g+2*rmv_redundant:2, dF, sender_id )
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
                do dF = 1, params%number_data_fields
                    hvy_block( Bs+g+1:Bs+g+g, 1:g, 1:g, dF, receiver_id ) = hvy_block( g+2:g+1+g, Bs:Bs-1+g, Bs:Bs-1+g, dF, sender_id )
                end do

            elseif ( level_diff == -1 ) then
                ! sender one level down
                ! interpolate data
                do dF = 1, params%number_data_fields
                    ! data to refine
                    data_corner = hvy_block( g+1:g+g, Bs+1:Bs+g, Bs+1:Bs+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_corner , data_corner_fine, params%order_predictor)
                    ! data to synchronize
                    data_corner = data_corner_fine(2:g+1, g-1:2*g-2, g-1:2*g-2)
                    ! write data
                    hvy_block( Bs+g+1:Bs+g+g, 1:g, 1:g, dF, receiver_id ) = data_corner
                end do

            elseif ( level_diff == 1) then
                ! sender one level up
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( Bs+g+1:Bs+g+g, 1:g, 1:g, dF, receiver_id ) = hvy_block( g+3:g+1+g+g:2, Bs-g:Bs-2+g:2, Bs-g:Bs-2+g:2, dF, sender_id )
                    hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, 1:g+rmv_redundant, 1:g+rmv_redundant, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:g+1+g+g:2, Bs-g:Bs-2+g+2*rmv_redundant:2, Bs-g:Bs-2+g+2*rmv_redundant:2, dF, sender_id )
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
                do dF = 1, params%number_data_fields
                    hvy_block( Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, 1:g, dF, receiver_id ) = hvy_block( g+2:g+1+g, g+2:g+1+g, Bs:Bs-1+g, dF, sender_id )
                end do

            elseif ( level_diff == -1 ) then
                ! sender one level down
                ! interpolate data
                do dF = 1, params%number_data_fields
                    ! data to refine
                    data_corner = hvy_block( g+1:g+g, g+1:g+g, Bs+1:Bs+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_corner , data_corner_fine, params%order_predictor)
                    ! data to synchronize
                    data_corner = data_corner_fine(2:g+1, 2:g+1, g-1:2*g-2)
                    ! write data
                    hvy_block( Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, 1:g, dF, receiver_id ) = data_corner
                end do

            elseif ( level_diff == 1) then
                ! sender one level up
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, 1:g, dF, receiver_id ) = hvy_block( g+3:g+1+g+g:2, g+3:g+1+g+g:2, Bs-g:Bs-2+g:2, dF, sender_id )
                    hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, Bs+g+1-rmv_redundant:Bs+g+g, 1:g+rmv_redundant, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:g+1+g+g:2, g+3-2*rmv_redundant:g+1+g+g:2, Bs-g:Bs-2+g+2*rmv_redundant:2, dF, sender_id )
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
                do dF = 1, params%number_data_fields
                    hvy_block( 1:g, Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( Bs:Bs-1+g, g+2:g+1+g, g+2:g+1+g, dF, sender_id )
                end do

            elseif ( level_diff == -1 ) then
                ! sender one level down
                ! interpolate data
                do dF = 1, params%number_data_fields
                    ! data to refine
                    data_corner = hvy_block( Bs+1:Bs+g, g+1:g+g, g+1:g+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_corner , data_corner_fine, params%order_predictor)
                    ! data to synchronize
                    data_corner = data_corner_fine(g-1:2*g-2, 2:g+1, 2:g+1)
                    ! write data
                    hvy_block( 1:g, Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_corner
                end do

            elseif ( level_diff == 1) then
                ! sender one level up
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( 1:g, Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( Bs-g:Bs-2+g:2, g+3:g+1+g+g:2, g+3:g+1+g+g:2, dF, sender_id )
                    hvy_block( 1:g+rmv_redundant, Bs+g+1-rmv_redundant:Bs+g+g, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( Bs-g:Bs-2+g+2*rmv_redundant:2, g+3-2*rmv_redundant:g+1+g+g:2, g+3-2*rmv_redundant:g+1+g+g:2, dF, sender_id )
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
                do dF = 1, params%number_data_fields
                    hvy_block( 1:g, 1:g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( Bs:Bs-1+g, Bs:Bs-1+g, g+2:g+1+g, dF, sender_id )
                end do

            elseif ( level_diff == -1 ) then
                ! sender one level down
                ! interpolate data
                do dF = 1, params%number_data_fields
                    ! data to refine
                    data_corner = hvy_block( Bs+1:Bs+g, Bs+1:Bs+g, g+1:g+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_corner , data_corner_fine, params%order_predictor)
                    ! data to synchronize
                    data_corner = data_corner_fine(g-1:2*g-2, g-1:2*g-2, 2:g+1)
                    ! write data
                    hvy_block( 1:g, 1:g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_corner
                end do

            elseif ( level_diff == 1) then
                ! sender one level up
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( 1:g, 1:g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( Bs-g:Bs-2+g:2, Bs-g:Bs-2+g:2, g+3:g+1+g+g:2, dF, sender_id )
                    hvy_block( 1:g+rmv_redundant, 1:g+rmv_redundant, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( Bs-g:Bs-2+g+2*rmv_redundant:2, Bs-g:Bs-2+g+2*rmv_redundant:2, g+3-2*rmv_redundant:g+1+g+g:2, dF, sender_id )
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
                do dF = 1, params%number_data_fields
                    hvy_block( Bs+g+1:Bs+g+g, 1:g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+2:g+1+g, Bs:Bs-1+g, g+2:g+1+g, dF, sender_id )
                end do

            elseif ( level_diff == -1 ) then
                ! sender one level down
                ! interpolate data
                do dF = 1, params%number_data_fields
                    ! data to refine
                    data_corner = hvy_block( g+1:g+g, Bs+1:Bs+g, g+1:g+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_corner , data_corner_fine, params%order_predictor)
                    ! data to synchronize
                    data_corner = data_corner_fine(2:g+1, g-1:2*g-2, 2:g+1)
                    ! write data
                    hvy_block( Bs+g+1:Bs+g+g, 1:g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_corner
                end do

            elseif ( level_diff == 1) then
                ! sender one level up
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( Bs+g+1:Bs+g+g, 1:g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+3:g+1+g+g:2, Bs-g:Bs-2+g:2, g+3:g+1+g+g:2, dF, sender_id )
                    hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, 1:g+rmv_redundant, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:g+1+g+g:2, Bs-g:Bs-2+g+2*rmv_redundant:2, g+3-2*rmv_redundant:g+1+g+g:2, dF, sender_id )
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
                do dF = 1, params%number_data_fields
                    hvy_block( Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+2:g+1+g, g+2:g+1+g, g+2:g+1+g, dF, sender_id )
                end do

            elseif ( level_diff == -1 ) then
                ! sender one level down
                ! interpolate data
                do dF = 1, params%number_data_fields
                    ! data to refine
                    data_corner = hvy_block( g+1:g+g, g+1:g+g, g+1:g+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_corner , data_corner_fine, params%order_predictor)
                    ! data to synchronize
                    data_corner = data_corner_fine(2:g+1, 2:g+1, 2:g+1)
                    ! write data
                    hvy_block( Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g,Bs+g+1:Bs+g+g, dF, receiver_id ) = data_corner
                end do

            elseif ( level_diff == 1) then
                ! sender one level up
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+3:g+1+g+g:2, g+3:g+1+g+g:2, g+3:g+1+g+g:2, dF, sender_id )
                    hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, Bs+g+1-rmv_redundant:Bs+g+g, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:g+1+g+g:2, g+3-2*rmv_redundant:g+1+g+g:2, g+3-2*rmv_redundant:g+1+g+g:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '__1/123'
        case(27)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
                    data_face = hvy_block( (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( g+1:Bs+g, g+1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, g+1:Bs+g, Bs:Bs-1+g)
                    !hvy_block( 1:Bs+g, g+1:Bs+2*g, 1:g, dF, receiver_id ) = data_face_fine(1:Bs+g, 1:Bs+g, Bs:Bs-1+g)
                    hvy_block( 1:Bs+g, g+1:Bs+2*g, 1:g, dF, receiver_id ) = data_face_fine(g+1:Bs+2*g, 1:Bs+g, Bs+g:Bs-1+2*g)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( g+(Bs+1)/2:Bs+g, g+1:g+(Bs+1)/2, 1:g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, Bs-g:Bs+g-2:2, dF, sender_id )
                    hvy_block( g+(Bs+1)/2:Bs+g, g+1:g+(Bs+1)/2, 1:g+rmv_redundant, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '__1/134'
        case(28)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
                    data_face = hvy_block( (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( g+1:Bs+g, g+1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, g+1:Bs+g, Bs:Bs-1+g)
                    !hvy_block( 1:Bs+g, 1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(1:Bs+g, 1:Bs+g, Bs:Bs-1+g)
                    hvy_block( 1:Bs+g, 1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(g+1:Bs+2*g, g+1:Bs+2*g, Bs+g:Bs-1+2*g)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( g+(Bs+1)/2:Bs+g, g+(Bs+1)/2:Bs+g, 1:g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, Bs-g:Bs+g-2:2, dF, sender_id )
                    hvy_block( g+(Bs+1)/2:Bs+g, g+(Bs+1)/2:Bs+g, 1:g+rmv_redundant, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '__1/145'
        case(29)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( g+1:Bs+g, g+1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, g+1:Bs+g, Bs:Bs-1+g)
                    !hvy_block( g+1:Bs+2*g, 1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(1:Bs+g, 1:Bs+g, Bs:Bs-1+g)
                    hvy_block( g+1:Bs+2*g, 1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(1:Bs+g, g+1:Bs+2*g, Bs+g:Bs-1+2*g)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( g+1:g+(Bs+1)/2, g+(Bs+1)/2:Bs+g, 1:g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, Bs-g:Bs+g-2:2, dF, sender_id )
                    hvy_block( g+1:g+(Bs+1)/2, g+(Bs+1)/2:Bs+g, 1:g+rmv_redundant, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '__1/152'
        case(30)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( g+1:Bs+g, g+1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, g+1:Bs+g, Bs:Bs-1+g)
                    !hvy_block( g+1:Bs+2*g, g+1:Bs+2*g, 1:g, dF, receiver_id ) = data_face_fine(1:Bs+g, 1:Bs+g, Bs:Bs-1+g)
                    hvy_block( g+1:Bs+2*g, g+1:Bs+2*g, 1:g, dF, receiver_id ) = data_face_fine(1:Bs+g, 1:Bs+g, Bs+g:Bs-1+2*g)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( g+1:g+(Bs+1)/2, g+1:g+(Bs+1)/2, 1:g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, Bs-g:Bs+g-2:2, dF, sender_id )
                    hvy_block( g+1:g+(Bs+1)/2, g+1:g+(Bs+1)/2, 1:g+rmv_redundant, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '__2/123'
        case(31)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
                    data_face = hvy_block( (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( g+1:Bs+g, Bs+g+1:Bs+g+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, 2:g+1, g+1:Bs+g)
                    !hvy_block( 1:Bs+g, Bs+g+1:Bs+g+g, 1:Bs+g, dF, receiver_id ) = data_face_fine(1:Bs+g, 2:g+1, 1:Bs+g)
                    hvy_block( 1:Bs+g, Bs+g+1:Bs+g+g, 1:Bs+g, dF, receiver_id ) = data_face_fine(g+1:Bs+2*g, 2:g+1, g+1:Bs+2*g)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( g+(Bs+1)/2:Bs+g, Bs+g+1:Bs+g+g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
                    hvy_block( g+(Bs+1)/2:Bs+g, Bs+g+1-rmv_redundant:Bs+g+g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '__2/623'
        case(32)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
                    data_face = hvy_block( (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( g+1:Bs+g, Bs+g+1:Bs+g+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, 2:g+1, g+1:Bs+g)
                    !hvy_block( 1:Bs+g, Bs+g+1:Bs+g+g, g+1:Bs+2*g, dF, receiver_id ) = data_face_fine(1:Bs+g, 2:g+1, 1:Bs+g)
                    hvy_block( 1:Bs+g, Bs+g+1:Bs+g+g, g+1:Bs+2*g, dF, receiver_id ) = data_face_fine(g+1:Bs+2*g, 2:g+1, 1:Bs+g)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( g+(Bs+1)/2:Bs+g, Bs+g+1:Bs+g+g, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
                    hvy_block( g+(Bs+1)/2:Bs+g, Bs+g+1-rmv_redundant:Bs+g+g, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '__2/152'
        case(33)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( g+1:Bs+g, Bs+g+1:Bs+g+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, 2:g+1, g+1:Bs+g)
                    !hvy_block( g+1:Bs+2*g, Bs+g+1:Bs+g+g, 1:Bs+g, dF, receiver_id ) = data_face_fine(1:Bs+g, 2:g+1, 1:Bs+g)
                    hvy_block( g+1:Bs+2*g, Bs+g+1:Bs+g+g, 1:Bs+g, dF, receiver_id ) = data_face_fine(1:Bs+g, 2:g+1, g+1:Bs+2*g)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( g+1:g+(Bs+1)/2, Bs+g+1:Bs+g+g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
                    hvy_block( g+1:g+(Bs+1)/2, Bs+g+1-rmv_redundant:Bs+g+g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '__2/652'
        case(34)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( g+1:Bs+g, Bs+g+1:Bs+g+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, 2:g+1, g+1:Bs+g)
                    hvy_block( g+1:Bs+2*g, Bs+g+1:Bs+g+g, g+1:Bs+2*g, dF, receiver_id ) = data_face_fine(1:Bs+g, 2:g+1, 1:Bs+g)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( g+1:g+(Bs+1)/2, Bs+g+1:Bs+g+g, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
                    hvy_block( g+1:g+(Bs+1)/2, Bs+g+1-rmv_redundant:Bs+g+g, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '__3/123'
        case(35)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
                    data_face = hvy_block( (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( 1:g, g+1:Bs+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, g+1:Bs+g, g+1:Bs+g)
                    !hvy_block( 1:g, g+1:Bs+2*g, 1:Bs+g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, 1:Bs+g, 1:Bs+g)
                    hvy_block( 1:g, g+1:Bs+2*g, 1:Bs+g, dF, receiver_id ) = data_face_fine(Bs+g:Bs-1+2*g, 1:Bs+g, g+1:Bs+2*g)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( 1:g, g+1:g+(Bs+1)/2, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
                    hvy_block( 1:g+rmv_redundant, g+1:g+(Bs+1)/2, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '__3/623'
        case(36)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
                    data_face = hvy_block( (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( 1:g, g+1:Bs+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, g+1:Bs+g, g+1:Bs+g)
                    !hvy_block( 1:g, g+1:Bs+2*g, g+1:Bs+2*g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, 1:Bs+g, 1:Bs+g)
                    hvy_block( 1:g, g+1:Bs+2*g, g+1:Bs+2*g, dF, receiver_id ) = data_face_fine(Bs+g:Bs-1+2*g, 1:Bs+g, 1:Bs+g)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( 1:g, g+1:g+(Bs+1)/2, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
                    hvy_block( 1:g+rmv_redundant, g+1:g+(Bs+1)/2, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '__3/134'
        case(37)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
                    data_face = hvy_block( (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( 1:g, g+1:Bs+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, g+1:Bs+g, g+1:Bs+g)
                    !hvy_block( 1:g, 1:Bs+g, 1:Bs+g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, 1:Bs+g, 1:Bs+g)
                    hvy_block( 1:g, 1:Bs+g, 1:Bs+g, dF, receiver_id ) = data_face_fine(Bs+g:Bs-1+2*g, g+1:Bs+2*g, g+1:Bs+2*g)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( 1:g, g+(Bs+1)/2:Bs+g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
                    hvy_block( 1:g+rmv_redundant, g+(Bs+1)/2:Bs+g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '__3/634'
        case(38)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
                    data_face = hvy_block( (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( 1:g, g+1:Bs+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, g+1:Bs+g, g+1:Bs+g)
                    !hvy_block( 1:g, 1:Bs+g, g+1:Bs+2*g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, 1:Bs+g, 1:Bs+g)
                    hvy_block( 1:g, 1:Bs+g, g+1:Bs+2*g, dF, receiver_id ) = data_face_fine(Bs+g:Bs-1+2*g, g+1:Bs+2*g, 1:Bs+g)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( 1:g, g+(Bs+1)/2:Bs+g, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
                    hvy_block( 1:g+rmv_redundant, g+(Bs+1)/2:Bs+g, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '__4/134'
        case(39)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
                    data_face = hvy_block( (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( g+1:Bs+g, 1:g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, Bs:Bs-1+g, g+1:Bs+g)
                    !hvy_block( 1:Bs+g, 1:g, 1:Bs+g, dF, receiver_id ) = data_face_fine(1:Bs+g, Bs:Bs-1+g, 1:Bs+g)
                    hvy_block( 1:Bs+g, 1:g, 1:Bs+g, dF, receiver_id ) = data_face_fine(g+1:Bs+2*g, Bs+g:Bs-1+2*g, g+1:Bs+2*g)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( g+(Bs+1)/2:Bs+g, 1:g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2:2, g+1:Bs+g:2, dF, sender_id )
                    hvy_block( g+(Bs+1)/2:Bs+g, 1:g+rmv_redundant, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '__4/634'
        case(40)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
                    data_face = hvy_block( (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( g+1:Bs+g, 1:g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, Bs:Bs-1+g, g+1:Bs+g)
                    !hvy_block( 1:Bs+g, 1:g, g+1:Bs+2*g, dF, receiver_id ) = data_face_fine(1:Bs+g, Bs:Bs-1+g, 1:Bs+g)
                    hvy_block( 1:Bs+g, 1:g, g+1:Bs+2*g, dF, receiver_id ) = data_face_fine(g+1:Bs+2*g, Bs+g:Bs-1+2*g, 1:Bs+g)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( g+(Bs+1)/2:Bs+g, 1:g, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2:2, g+1:Bs+g:2, dF, sender_id )
                    hvy_block( g+(Bs+1)/2:Bs+g, 1:g+rmv_redundant, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '__4/145'
        case(41)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( g+1:Bs+g, 1:g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, Bs:Bs-1+g, g+1:Bs+g)
                    !hvy_block( g+1:Bs+2*g, 1:g, 1:Bs+g, dF, receiver_id ) = data_face_fine(1:Bs+g, Bs:Bs-1+g, 1:Bs+g)
                    hvy_block( g+1:Bs+2*g, 1:g, 1:Bs+g, dF, receiver_id ) = data_face_fine(1:Bs+g, Bs+g:Bs-1+2*g, g+1:Bs+2*g)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( g+1:g+(Bs+1)/2, 1:g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2:2, g+1:Bs+g:2, dF, sender_id )
                    hvy_block( g+1:g+(Bs+1)/2, 1:g+rmv_redundant, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '__4/645'
        case(42)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( g+1:Bs+g, 1:g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, Bs:Bs-1+g, g+1:Bs+g)
                    !hvy_block( g+1:Bs+2*g, 1:g, g+1:Bs+2*g, dF, receiver_id ) = data_face_fine(1:Bs+g, Bs:Bs-1+g, 1:Bs+g)
                    hvy_block( g+1:Bs+2*g, 1:g, g+1:Bs+2*g, dF, receiver_id ) = data_face_fine(1:Bs+g, Bs+g:Bs-1+2*g, 1:Bs+g)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( g+1:g+(Bs+1)/2, 1:g, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2:2, g+1:Bs+g:2, dF, sender_id )
                    hvy_block( g+1:g+(Bs+1)/2, 1:g+rmv_redundant, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '__5/145'
        case(43)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(2:g+1, g+1:Bs+g, g+1:Bs+g)
                    !hvy_block( Bs+g+1:Bs+g+g, 1:Bs+g, 1:Bs+g, dF, receiver_id ) = data_face_fine(2:g+1, 1:Bs+g, 1:Bs+g)
                    hvy_block( Bs+g+1:Bs+g+g, 1:Bs+g, 1:Bs+g, dF, receiver_id ) = data_face_fine(2:g+1, g+1:Bs+2*g, g+1:Bs+2*g)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( Bs+g+1:Bs+g+g, g+(Bs+1)/2:Bs+g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+3:3*g+1:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
                    hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, g+(Bs+1)/2:Bs+g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '__5/645'
        case(44)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(2:g+1, g+1:Bs+g, g+1:Bs+g)
                    !hvy_block( Bs+g+1:Bs+g+g, 1:Bs+g, g+1:Bs+2*g, dF, receiver_id ) = data_face_fine(2:g+1, 1:Bs+g, 1:Bs+g)
                    hvy_block( Bs+g+1:Bs+g+g, 1:Bs+g, g+1:Bs+2*g, dF, receiver_id ) = data_face_fine(2:g+1, g+1:Bs+2*g, 1:Bs+g)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( Bs+g+1:Bs+g+g, g+(Bs+1)/2:Bs+g, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+3:3*g+1:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
                    hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, g+(Bs+1)/2:Bs+g, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '__5/152'
        case(45)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(2:g+1, g+1:Bs+g, g+1:Bs+g)
                    !hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+2*g, 1:Bs+g, dF, receiver_id ) = data_face_fine(2:g+1, 1:Bs+g, 1:Bs+g)
                    hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+2*g, 1:Bs+g, dF, receiver_id ) = data_face_fine(2:g+1, 1:Bs+g, g+1:Bs+2*g)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( Bs+g+1:Bs+g+g, g+1:g+(Bs+1)/2, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+3:3*g+1:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
                    hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, g+1:g+(Bs+1)/2, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '__5/652'
        case(46)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(2:g+1, g+1:Bs+g, g+1:Bs+g)
                    hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+2*g, g+1:Bs+2*g, dF, receiver_id ) = data_face_fine(2:g+1, 1:Bs+g, 1:Bs+g)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( Bs+g+1:Bs+g+g, g+1:g+(Bs+1)/2, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+3:3*g+1:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
                    hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, g+1:g+(Bs+1)/2, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, g+1:Bs+g:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '__6/623'
        case(47)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g , g+1:(Bs+1)/2+g/2+g, dF, sender_id )
                    data_face = hvy_block( (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g , g+1:(Bs+1)/2+g+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( g+1:Bs+g, g+1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, g+1:Bs+g, 2:g+1)
                    !hvy_block( 1:Bs+g, g+1:Bs+2*g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(1:Bs+g, 1:Bs+g, 2:g+1)
                    hvy_block( 1:Bs+g, g+1:Bs+2*g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(g+1:Bs+2*g, 1:Bs+g, 2:g+1)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( g+(Bs+1)/2:Bs+g, g+1:g+(Bs+1)/2, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, g+3:3*g+1:2, dF, sender_id )
                    hvy_block( g+(Bs+1)/2:Bs+g, g+1:g+(Bs+1)/2, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '__6/634'
        case(48)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
                    data_face = hvy_block( (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( g+1:Bs+g, g+1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, g+1:Bs+g, 2:g+1)
                    !hvy_block( 1:Bs+g, 1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(1:Bs+g, 1:Bs+g, 2:g+1)
                    hvy_block( 1:Bs+g, 1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(g+1:Bs+2*g, g+1:Bs+2*g, 2:g+1)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( g+(Bs+1)/2:Bs+g, g+(Bs+1)/2:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, g+3:3*g+1:2, dF, sender_id )
                    hvy_block( g+(Bs+1)/2:Bs+g, g+(Bs+1)/2:Bs+g, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '__6/645'
        case(49)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( g+1:Bs+g, g+1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, g+1:Bs+g, 2:g+1)
                    !hvy_block( g+1:Bs+2*g, 1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(1:Bs+g, 1:Bs+g, 2:g+1)
                    hvy_block( g+1:Bs+2*g, 1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(1:Bs+g, g+1:Bs+2*g, 2:g+1)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( g+1:g+(Bs+1)/2, g+(Bs+1)/2:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, g+3:3*g+1:2, dF, sender_id )
                    hvy_block( g+1:g+(Bs+1)/2, g+(Bs+1)/2:Bs+g, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '__6/652'
        case(50)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate
                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( g+1:Bs+g, g+1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, g+1:Bs+g, 2:g+1)
                    hvy_block( g+1:Bs+2*g, g+1:Bs+2*g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(1:Bs+g, 1:Bs+g, 2:g+1)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( g+1:g+(Bs+1)/2, g+1:g+(Bs+1)/2, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, g+3:3*g+1:2, dF, sender_id )
                    hvy_block( g+1:g+(Bs+1)/2, g+1:g+(Bs+1)/2, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '_12/123'
        case(51)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate, note: use data_face interpolation variable
                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
                    data_face = hvy_block( (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( g+1:Bs+g, Bs+g+1:Bs+g+g, 1:g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, 2:g+1, Bs:Bs-1+g)
                    hvy_block( g+1:Bs+g, Bs+g+1:Bs+g+g, 1:g, dF, receiver_id ) = data_face_fine(g+g+1:Bs+2*g, 2:g+1, Bs+g:Bs-1+2*g)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( g+(Bs+1)/2:Bs+g, Bs+g+1:Bs+g+g, 1:g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3:3*g+1:2, Bs-g:Bs+g-2:2, dF, sender_id )
                    hvy_block( g+(Bs+1)/2:Bs+g, Bs+g+1-rmv_redundant:Bs+g+g, 1:g+rmv_redundant, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, Bs-g:Bs+g-2+2*rmv_redundant:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '_12/152'
        case(52)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate, note: use data_face interpolation variable
                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( g+1:Bs+g, Bs+g+1:Bs+g+g, 1:g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, 2:g+1, Bs:Bs-1+g)
                    !hvy_block( g+1:Bs+g, Bs+g+1:Bs+g+g, 1:g, dF, receiver_id ) = data_face_fine(1:Bs, 2:g+1, Bs:Bs-1+g)
                    hvy_block( g+1:Bs+g, Bs+g+1:Bs+g+g, 1:g, dF, receiver_id ) = data_face_fine(1:Bs, 2:g+1, Bs+g:Bs-1+2*g)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( g+1:g+(Bs+1)/2, Bs+g+1:Bs+g+g, 1:g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3:3*g+1:2, Bs-g:Bs+g-2:2, dF, sender_id )
                    hvy_block( g+1:g+(Bs+1)/2, Bs+g+1-rmv_redundant:Bs+g+g, 1:g+rmv_redundant, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, Bs-g:Bs+g-2+2*rmv_redundant:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '_13/123'
        case(53)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate, note: use data_face interpolation variable
                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
                    data_face = hvy_block( (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( 1:g, g+1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, g+1:Bs+g, Bs:Bs-1+g)
                    !hvy_block( 1:g, g+1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, 1:Bs, Bs:Bs-1+g)
                    hvy_block( 1:g, g+1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(Bs+g:Bs-1+2*g, 1:Bs, Bs+g:Bs-1+2*g)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( 1:g, g+1:g+(Bs+1)/2, 1:g, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2:2, g+1:Bs+g:2, Bs-g:Bs+g-2:2, dF, sender_id )
                    hvy_block( 1:g+rmv_redundant, g+1:g+(Bs+1)/2, 1:g+rmv_redundant, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '_13/134'
        case(54)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate, note: use data_face interpolation variable
                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
                    data_face = hvy_block( (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( 1:g, g+1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, g+1:Bs+g, Bs:Bs-1+g)
                    hvy_block( 1:g, g+1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(Bs+g:Bs-1+2*g, g+g+1:Bs+2*g, Bs+g:Bs-1+2*g)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( 1:g, g+(Bs+1)/2:Bs+g, 1:g, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2:2, g+1:Bs+g:2, Bs-g:Bs+g-2:2, dF, sender_id )
                    hvy_block( 1:g+rmv_redundant, g+(Bs+1)/2:Bs+g, 1:g+rmv_redundant, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '_14/134'
        case(55)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate, note: use data_face interpolation variable
                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
                    data_face = hvy_block( (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( g+1:Bs+g, 1:g, 1:g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, Bs:Bs-1+g, Bs:Bs-1+g)
                    hvy_block( g+1:Bs+g, 1:g, 1:g, dF, receiver_id ) = data_face_fine(g+g+1:Bs+2*g, Bs+g:Bs-1+2*g, Bs+g:Bs-1+2*g)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( g+(Bs+1)/2:Bs+g, 1:g, 1:g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2:2, Bs-g:Bs+g-2:2, dF, sender_id )
                    hvy_block( g+(Bs+1)/2:Bs+g, 1:g+rmv_redundant, 1:g+rmv_redundant, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, Bs-g:Bs+g-2+2*rmv_redundant:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '_14/145'
        case(56)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate, note: use data_face interpolation variable
                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( g+1:Bs+g, 1:g, 1:g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, Bs:Bs-1+g, Bs:Bs-1+g)
                    !hvy_block( g+1:Bs+g, 1:g, 1:g, dF, receiver_id ) = data_face_fine(1:Bs, Bs:Bs-1+g, Bs:Bs-1+g)
                    hvy_block( g+1:Bs+g, 1:g, 1:g, dF, receiver_id ) = data_face_fine(1:Bs, Bs+g:Bs-1+2*g, Bs+g:Bs-1+2*g)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( g+1:g+(Bs+1)/2, 1:g, 1:g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2:2, Bs-g:Bs+g-2:2, dF, sender_id )
                    hvy_block( g+1:g+(Bs+1)/2, 1:g+rmv_redundant, 1:g+rmv_redundant, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, Bs-g:Bs+g-2+2*rmv_redundant:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '_15/145'
        case(57)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate, note: use data_face interpolation variable
                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(2:g+1, g+1:Bs+g, Bs:Bs-1+g)
                    hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(2:g+1, g+g+1:Bs+2*g, Bs+g:Bs-1+2*g)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( Bs+g+1:Bs+g+g, g+(Bs+1)/2:Bs+g, 1:g, dF, receiver_id ) = hvy_block( g+3:3*g+1:2, g+1:Bs+g:2, Bs-g:Bs+g-2:2, dF, sender_id )
                    hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, g+(Bs+1)/2:Bs+g, 1:g+rmv_redundant, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '_15/152'
        case(58)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate, note: use data_face interpolation variable
                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(2:g+1, g+1:Bs+g, Bs:Bs-1+g)
                    !hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(2:g+1, 1:Bs, Bs:Bs-1+g)
                    hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+g, 1:g, dF, receiver_id ) = data_face_fine(2:g+1, 1:Bs, Bs+g:Bs-1+2*g)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( Bs+g+1:Bs+g+g, g+1:g+(Bs+1)/2, 1:g, dF, receiver_id ) = hvy_block( g+3:3*g+1:2, g+1:Bs+g:2, Bs-g:Bs+g-2:2, dF, sender_id )
                    hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, g+1:g+(Bs+1)/2, 1:g+rmv_redundant, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

         ! '_62/623'
        case(59)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate, note: use data_face interpolation variable
                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
                    data_face = hvy_block( (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( g+1:Bs+g, Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, 2:g+1, 2:g+1)
                    hvy_block( g+1:Bs+g, Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(g+g+1:Bs+2*g, 2:g+1, 2:g+1)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( g+(Bs+1)/2:Bs+g, Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3:3*g+1:2, g+3:3*g+1:2, dF, sender_id )
                    hvy_block( g+(Bs+1)/2:Bs+g, Bs+g+1-rmv_redundant:Bs+g+g, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, g+3-2*rmv_redundant:3*g+1:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '_62/652'
        case(60)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate, note: use data_face interpolation variable
                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( g+1:Bs+g, Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, 2:g+1, 2:g+1)
                    hvy_block( g+1:Bs+g, Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(1:Bs, 2:g+1, 2:g+1)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( g+1:g+(Bs+1)/2, Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3:3*g+1:2, g+3:3*g+1:2, dF, sender_id )
                    hvy_block( g+1:g+(Bs+1)/2, Bs+g+1-rmv_redundant:Bs+g+g, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, g+3-2*rmv_redundant:3*g+1:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '_63/623'
        case(61)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate, note: use data_face interpolation variable
                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
                    data_face = hvy_block( (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( 1:g, g+1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, g+1:Bs+g, 2:g+1)
                    !hvy_block( 1:g, g+1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, 1:Bs, 2:g+1)
                    hvy_block( 1:g, g+1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(Bs+g:Bs-1+2*g, 1:Bs, 2:g+1)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( 1:g, g+1:g+(Bs+1)/2, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2:2, g+1:Bs+g:2, g+3:3*g+1:2, dF, sender_id )
                    hvy_block( 1:g+rmv_redundant, g+1:g+(Bs+1)/2, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '_63/634'
        case(62)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate, note: use data_face interpolation variable
                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
                    data_face = hvy_block( (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( 1:g, g+1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, g+1:Bs+g, 2:g+1)
                    hvy_block( 1:g, g+1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(Bs+g:Bs-1+2*g, g+g+1:Bs+2*g, 2:g+1)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( 1:g, g+(Bs+1)/2:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2:2, g+1:Bs+g:2, g+3:3*g+1:2, dF, sender_id )
                    hvy_block( 1:g+rmv_redundant, g+(Bs+1)/2:Bs+g, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '_64/634'
        case(63)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate, note: use data_face interpolation variable
                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
                    data_face = hvy_block( (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( g+1:Bs+g, 1:g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, Bs:Bs-1+g, 2:g+1)
                    hvy_block( g+1:Bs+g, 1:g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(g+g+1:Bs+2*g, Bs+g:Bs-1+2*g, 2:g+1)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( g+(Bs+1)/2:Bs+g, 1:g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2:2, g+3:3*g+1:2, dF, sender_id )
                    hvy_block( g+(Bs+1)/2:Bs+g, 1:g+rmv_redundant, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, g+3-2*rmv_redundant:3*g+1:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '_64/645'
        case(64)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate, note: use data_face interpolation variable
                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( g+1:Bs+g, 1:g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(g+1:Bs+g, Bs:Bs-1+g, 2:g+1)
                    !hvy_block( g+1:Bs+g, 1:g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(1:Bs, Bs:Bs-1+g, 2:g+1)
                    hvy_block( g+1:Bs+g, 1:g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(1:Bs, Bs+g:Bs-1+2*g, 2:g+1)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( g+1:g+(Bs+1)/2, 1:g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2:2, g+3:3*g+1:2, dF, sender_id )
                    hvy_block( g+1:g+(Bs+1)/2, 1:g+rmv_redundant, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( g+1:Bs+g:2, Bs-g:Bs+g-2+2*rmv_redundant:2, g+3-2*rmv_redundant:3*g+1:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '_65/645'
        case(65)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate, note: use data_face interpolation variable
                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(2:g+1, g+1:Bs+g, 2:g+1)
                    hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(2:g+1, g+g+1:Bs+2*g, 2:g+1)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( Bs+g+1:Bs+g+g, g+(Bs+1)/2:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+3:3*g+1:2, g+1:Bs+g:2, g+3:3*g+1:2, dF, sender_id )
                    hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, g+(Bs+1)/2:Bs+g, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '_65/652'
        case(66)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate, note: use data_face interpolation variable
                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(2:g+1, g+1:Bs+g, 2:g+1)
                    hvy_block( Bs+g+1:Bs+g+g, g+1:Bs+g, Bs+g+1:Bs+g+g, dF, receiver_id ) = data_face_fine(2:g+1, 1:Bs, 2:g+1)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( Bs+g+1:Bs+g+g, g+1:g+(Bs+1)/2, Bs+g+1:Bs+g+g, dF, receiver_id ) = hvy_block( g+3:3*g+1:2, g+1:Bs+g:2, g+3:3*g+1:2, dF, sender_id )
                    hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, g+1:g+(Bs+1)/2, Bs+g+1-rmv_redundant:Bs+g+g, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, g+3-2*rmv_redundant:3*g+1:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '_23/123'
        case(67)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate, note: use data_face interpolation variable
                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
                    data_face = hvy_block( (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( 1:g, Bs+g+1:Bs+g+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, 2:g+1, g+1:Bs+g)
                    hvy_block( 1:g, Bs+g+1:Bs+g+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(Bs+g:Bs-1+2*g, 2:g+1, g+g+1:Bs+2*g)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( 1:g, Bs+g+1:Bs+g+g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2:2, g+3:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
                    hvy_block( 1:g+rmv_redundant, Bs+g+1-rmv_redundant:Bs+g+g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2+2*rmv_redundant:2, g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '_23/623'
        case(68)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate, note: use data_face interpolation variable
                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
                    data_face = hvy_block( (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( 1:g, Bs+g+1:Bs+g+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, 2:g+1, g+1:Bs+g)
                    !hvy_block( 1:g, Bs+g+1:Bs+g+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, 2:g+1, 1:Bs)
                    hvy_block( 1:g, Bs+g+1:Bs+g+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(Bs+g:Bs-1+2*g, 2:g+1, 1:Bs)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( 1:g, Bs+g+1:Bs+g+g, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2:2, g+3:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
                    hvy_block( 1:g+rmv_redundant, Bs+g+1-rmv_redundant:Bs+g+g, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2+2*rmv_redundant:2, g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '_25/152'
        case(69)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate, note: use data_face interpolation variable
                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(2:g+1, 2:g+1, g+1:Bs+g)
                    hvy_block( Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(2:g+1, 2:g+1, g+g+1:Bs+2*g)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+3:3*g+1:2, g+3:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
                    hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, Bs+g+1-rmv_redundant:Bs+g+g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:3*g+1:2, g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '_25/652'
        case(70)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate, note: use data_face interpolation variable
                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(2:g+1, 2:g+1, g+1:Bs+g)
                    hvy_block( Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(2:g+1, 2:g+1, 1:Bs)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( Bs+g+1:Bs+g+g, Bs+g+1:Bs+g+g, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+3:3*g+1:2, g+3:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
                    hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, Bs+g+1-rmv_redundant:Bs+g+g, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:3*g+1:2, g+3-2*rmv_redundant:3*g+1:2, g+1:Bs+g:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '_43/134'
        case(71)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate, note: use data_face interpolation variable
                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
                    data_face = hvy_block( (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( 1:g, 1:g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, Bs:Bs-1+g, g+1:Bs+g)
                    hvy_block( 1:g, 1:g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(Bs+g:Bs-1+2*g, Bs+g:Bs-1+2*g, g+g+1:Bs+2*g)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( 1:g, 1:g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2:2, Bs-g:Bs+g-2:2, g+1:Bs+g:2, dF, sender_id )
                    hvy_block( 1:g+rmv_redundant, 1:g+rmv_redundant, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2+2*rmv_redundant:2, Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '_43/634'
        case(72)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate, note: use data_face interpolation variable
                    !data_face = hvy_block( (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
                    data_face = hvy_block( (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( 1:g, 1:g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, Bs:Bs-1+g, g+1:Bs+g)
                    !hvy_block( 1:g, 1:g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(Bs:Bs-1+g, Bs:Bs-1+g, 1:Bs)
                    hvy_block( 1:g, 1:g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(Bs+g:Bs-1+2*g, Bs+g:Bs-1+2*g, 1:Bs)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( 1:g, 1:g, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2:2, Bs-g:Bs+g-2:2, g+1:Bs+g:2, dF, sender_id )
                    hvy_block( 1:g+rmv_redundant, 1:g+rmv_redundant, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( Bs-g:Bs+g-2+2*rmv_redundant:2, Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '_45/145'
        case(73)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate, note: use data_face interpolation variable
                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, (Bs+1)/2+g/2:Bs+g, dF, sender_id )
                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, (Bs+1)/2:Bs+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( Bs+g+1:Bs+g+g, 1:g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(2:g+1, Bs:Bs-1+g, g+1:Bs+g)
                    hvy_block( Bs+g+1:Bs+g+g, 1:g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(2:g+1, Bs+g:Bs-1+2*g, g+g+1:Bs+2*g)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( Bs+g+1:Bs+g+g, 1:g, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+3:3*g+1:2, Bs-g:Bs+g-2:2, g+1:Bs+g:2, dF, sender_id )
                    hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, 1:g+rmv_redundant, g+(Bs+1)/2:Bs+g, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:3*g+1:2, Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

        ! '_45/645'
        case(74)
            if ( level_diff == -1 ) then
                ! sender on lower level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    ! data to interpolate, note: use data_face interpolation variable
                    !data_face = hvy_block( g+1:(Bs+1)/2+g/2+g, (Bs+1)/2+g/2:Bs+g, g+1:(Bs+1)/2+g/2+g, dF, sender_id )
                    data_face = hvy_block( g+1:(Bs+1)/2+g+g, (Bs+1)/2:Bs+g, g+1:(Bs+1)/2+g+g, dF, sender_id )
                    ! interpolate data
                    call prediction_3D( data_face , data_face_fine, params%order_predictor)
                    ! copy data
                    !hvy_block( Bs+g+1:Bs+g+g, 1:g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(2:g+1, Bs:Bs-1+g, g+1:Bs+g)
                    !hvy_block( Bs+g+1:Bs+g+g, 1:g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(2:g+1, Bs:Bs-1+g, 1:Bs)
                    hvy_block( Bs+g+1:Bs+g+g, 1:g, g+1:Bs+g, dF, receiver_id ) = data_face_fine(2:g+1, Bs+g:Bs-1+2*g, 1:Bs)
                end do

            elseif ( level_diff == 1 ) then
                ! sender on higher level
                ! loop over all datafields
                do dF = 1, params%number_data_fields
                    !hvy_block( Bs+g+1:Bs+g+g, 1:g, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+3:3*g+1:2, Bs-g:Bs+g-2:2, g+1:Bs+g:2, dF, sender_id )
                    hvy_block( Bs+g+1-rmv_redundant:Bs+g+g, 1:g+rmv_redundant, g+1:g+(Bs+1)/2, dF, receiver_id ) = hvy_block( g+3-2*rmv_redundant:3*g+1:2, Bs-g:Bs+g-2+2*rmv_redundant:2, g+1:Bs+g:2, dF, sender_id )
                end do

            else
                ! error case
                write(*,'(80("_"))')
                write(*,*) "ERROR: can not synchronize ghost nodes, mesh is not graded"
                stop
            end if

    end select

    ! clean up
    deallocate( data_corner  )
    deallocate( data_corner_fine  )
    deallocate( data_face  )
    deallocate( data_face_fine  )

end subroutine copy_ghost_nodes_3D
