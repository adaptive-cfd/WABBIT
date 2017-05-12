!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name update_neighbors_3D.f90
!> \version 0.5
!> \author msr
!
!>\brief update neighbor relations with light data, store neighbors in neighbor list (heavy data)
!!       3D case
!
!> 
!! input:    
!!           - light data array
!!           - params struct
!!
!! output:   
!!           - neighbor list array
!!
!! \n
!!
!  --------------------------------------------------------------------------------------------
!> neighbor codes: \n
!  --------------- 
!> for imagination:  
!!                   - 6-sided dice with '1'-side on top, '6'-side on bottom, '2'-side in front 
!!                   - edge: boundary between two sides - use sides numbers for coding
!!                   - corner: between three sides - so use all three sides numbers
!!                   - block on higher/lower level: block shares face/edge and one unique corner,
!!                     so use this corner code in second part of neighbor code
!!
!! faces:  '__1/___', '__2/___', '__3/___', '__4/___', '__5/___', '__6/___' \n
!! edges:  '_12/___', '_13/___', '_14/___', '_15/___' \n
!!         '_62/___', '_63/___', '_64/___', '_65/___' \n
!!         '_23/___', '_25/___', '_43/___', '_45/___' \n
!! corner: '123/___', '134/___', '145/___', '152/___' \n
!!         '623/___', '634/___', '645/___', '652/___' \n
!! \n
!! complete neighbor code array, 74 possible neighbor relations \n
!! neighbors = (/'__1/___', '__2/___', '__3/___', '__4/___', '__5/___', '__6/___', '_12/___', '_13/___', '_14/___', '_15/___',
!!               '_62/___', '_63/___', '_64/___', '_65/___', '_23/___', '_25/___', '_43/___', '_45/___', '123/___', '134/___',
!!               '145/___', '152/___', '623/___', '634/___', '645/___', '652/___', '__1/123', '__1/134', '__1/145', '__1/152',
!!               '__2/123', '__2/623', '__2/152', '__2/652', '__3/123', '__3/623', '__3/134', '__3/634', '__4/134', '__4/634',
!!               '__4/145', '__4/645', '__5/145', '__5/645', '__5/152', '__5/652', '__6/623', '__6/634', '__6/645', '__6/652',
!!               '_12/123', '_12/152', '_13/123', '_13/134', '_14/134', '_14/145', '_15/145', '_15/152', '_62/623', '_62/652',
!!               '_63/623', '_63/634', '_64/634', '_64/645', '_65/645', '_65/652', '_23/123', '_23/623', '_25/152', '_25/652',
!!               '_43/134', '_43/634', '_45/145', '_45/645' /) \n
!
!> = log ======================================================================================
!! \n
!! 27/01/17 - create
!
! ********************************************************************************************
!> \image html update_neighbors3d.png "Update Neighbors" width=400

subroutine update_neighbors_3D(params, lgt_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n)

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> heavy data array - neifghbor data
    integer(kind=ik), intent(out)       :: hvy_neighbor(:,:)

    !> list of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_n
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    ! number of blocks per proc
    integer(kind=ik)                    :: N
    ! max treelevel
    integer(kind=ik)                    :: max_treelevel

    ! process rank
    integer(kind=ik)                    :: rank

    ! loop variable
    integer(kind=ik)                    :: k, lgt_id

    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0, sub_t1

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! set MPI parameters
    rank         = params%rank

    ! reset neighbor list
    hvy_neighbor = -1

    ! number of blocks per proc
    N = params%number_blocks

    ! max treelevel
    max_treelevel = params%max_treelevel

!---------------------------------------------------------------------------------------------
! main body

    ! start time
    sub_t0 = MPI_Wtime()

    ! special case:
    ! if there is only one block => all neighbors are this block
    ! one block criteria: size of block_list should be one!
    if ( lgt_n == 1 ) then
        hvy_neighbor(1,1:26) = lgt_active(1)
    end if

    ! loop only if more than one block is active
    if ( lgt_n /= 1 ) then
        ! loop over active heavy data blocks
        do k = 1, hvy_n

            ! light id
            call hvy_id_to_lgt_id( lgt_id, hvy_active(k), rank, N )

            ! faces
            ! '1'-side
            call find_neighbor_face_3D(hvy_active(k), lgt_id, lgt_block, max_treelevel, '__1/___', hvy_neighbor, lgt_active, lgt_n)
            ! '2'-side
            call find_neighbor_face_3D(hvy_active(k), lgt_id, lgt_block, max_treelevel, '__2/___', hvy_neighbor, lgt_active, lgt_n)
            ! '3'-side
            call find_neighbor_face_3D(hvy_active(k), lgt_id, lgt_block, max_treelevel, '__3/___', hvy_neighbor, lgt_active, lgt_n)
            ! '4'-side
            call find_neighbor_face_3D(hvy_active(k), lgt_id, lgt_block, max_treelevel, '__4/___', hvy_neighbor, lgt_active, lgt_n)
            ! '5'-side
            call find_neighbor_face_3D(hvy_active(k), lgt_id, lgt_block, max_treelevel, '__5/___', hvy_neighbor, lgt_active, lgt_n)
            ! '6'-side
            call find_neighbor_face_3D(hvy_active(k), lgt_id, lgt_block, max_treelevel, '__6/___', hvy_neighbor, lgt_active, lgt_n)

            ! edges
            ! '12'-edge
            call find_neighbor_edge_3D(hvy_active(k), lgt_id, lgt_block, max_treelevel, '_12/___', hvy_neighbor, lgt_active, lgt_n)
            ! '13'-edge
            call find_neighbor_edge_3D(hvy_active(k), lgt_id, lgt_block, max_treelevel, '_13/___', hvy_neighbor, lgt_active, lgt_n)
            ! '14'-edge
            call find_neighbor_edge_3D(hvy_active(k), lgt_id, lgt_block, max_treelevel, '_14/___', hvy_neighbor, lgt_active, lgt_n)
            ! '15'-edge
            call find_neighbor_edge_3D(hvy_active(k), lgt_id, lgt_block, max_treelevel, '_15/___', hvy_neighbor, lgt_active, lgt_n)
            ! '62'-edge
            call find_neighbor_edge_3D(hvy_active(k), lgt_id, lgt_block, max_treelevel, '_62/___', hvy_neighbor, lgt_active, lgt_n)
            ! '63'-edge
            call find_neighbor_edge_3D(hvy_active(k), lgt_id, lgt_block, max_treelevel, '_63/___', hvy_neighbor, lgt_active, lgt_n)
            ! '64'-edge
            call find_neighbor_edge_3D(hvy_active(k), lgt_id, lgt_block, max_treelevel, '_64/___', hvy_neighbor, lgt_active, lgt_n)
            ! '65'-edge
            call find_neighbor_edge_3D(hvy_active(k), lgt_id, lgt_block, max_treelevel, '_65/___', hvy_neighbor, lgt_active, lgt_n)
            ! '23'-edge
            call find_neighbor_edge_3D(hvy_active(k), lgt_id, lgt_block, max_treelevel, '_23/___', hvy_neighbor, lgt_active, lgt_n)
            ! '25'-edge
            call find_neighbor_edge_3D(hvy_active(k), lgt_id, lgt_block, max_treelevel, '_25/___', hvy_neighbor, lgt_active, lgt_n)
            ! '43'-edge
            call find_neighbor_edge_3D(hvy_active(k), lgt_id, lgt_block, max_treelevel, '_43/___', hvy_neighbor, lgt_active, lgt_n)
            ! '45'-edge
            call find_neighbor_edge_3D(hvy_active(k), lgt_id, lgt_block, max_treelevel, '_45/___', hvy_neighbor, lgt_active, lgt_n)

            ! corners
            ! '123'-corner
            call find_neighbor_corner_3D(hvy_active(k), lgt_id, lgt_block, max_treelevel, '123/___', hvy_neighbor, lgt_active, lgt_n)
            ! '134'-corner
            call find_neighbor_corner_3D(hvy_active(k), lgt_id, lgt_block, max_treelevel, '134/___', hvy_neighbor, lgt_active, lgt_n)
            ! '145'-corner
            call find_neighbor_corner_3D(hvy_active(k), lgt_id, lgt_block, max_treelevel, '145/___', hvy_neighbor, lgt_active, lgt_n)
            ! '152'-corner
            call find_neighbor_corner_3D(hvy_active(k), lgt_id, lgt_block, max_treelevel, '152/___', hvy_neighbor, lgt_active, lgt_n)
            ! '623'-corner
            call find_neighbor_corner_3D(hvy_active(k), lgt_id, lgt_block, max_treelevel, '623/___', hvy_neighbor, lgt_active, lgt_n)
            ! '634'-corner
            call find_neighbor_corner_3D(hvy_active(k), lgt_id, lgt_block, max_treelevel, '634/___', hvy_neighbor, lgt_active, lgt_n)
            ! '645'-corner
            call find_neighbor_corner_3D(hvy_active(k), lgt_id, lgt_block, max_treelevel, '645/___', hvy_neighbor, lgt_active, lgt_n)
            ! '652'-corner
            call find_neighbor_corner_3D(hvy_active(k), lgt_id, lgt_block, max_treelevel, '652/___', hvy_neighbor, lgt_active, lgt_n)

        end do
    end if

    ! end time
    sub_t1 = MPI_Wtime()
    ! write time
    if ( params%debug ) then
        ! find free or corresponding line
        k = 1
        do while ( debug%name_comp_time(k) /= "---" )
            ! entry for current subroutine exists
            if ( debug%name_comp_time(k) == "update_neighbors" ) exit
            k = k + 1
        end do
        ! write time
        debug%name_comp_time(k) = "update_neighbors"
        debug%comp_time(k, 1)   = debug%comp_time(k, 1) + 1
        debug%comp_time(k, 2)   = debug%comp_time(k, 2) + sub_t1 - sub_t0
    end if

end subroutine update_neighbors_3D
