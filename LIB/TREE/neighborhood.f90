! JB: I would love to slowly collect all things related to neighborhood relations without params in this file, maybe make it a module if it grows
! JB: Yes, this will include calc_data_bounds at some points I hope (or fear?)


!> \brief Convert neighborhood to corresponding patch
!> The selection of ID is taken from the neighborhood relation, the structure goes from -2+, x2z, side2edge2corner
! How I have missed this beast, relations taken from dice-3d-neighbors and neighbor_codes
subroutine neighborhood2patchID(neighborhood, patchID, dim)
    implicit none

    integer(kind=ik), intent(in)   :: neighborhood  !< neighborhood of significant patch
    integer(kind=ik), intent(out)  :: patchID       !< location of neighbor, where is the neighbor physically?
    integer(kind=ik), intent(in)   :: dim           !< params%dim

    if (dim == 2) then
        ! 2D
        select case(neighborhood)
        case (1, 9:10)  ! -x
            patchID = 1
        case (3, 11:12) ! +x
            patchID = 2
        case (4, 15:16) ! -y
            patchID = 3
        case (2, 13:14) ! +y
            patchID = 4
        case(6) ! -x-y
            patchID = 5
        case(5) ! -x+y
            patchID = 6
        case(8) ! +x-y
            patchID = 7
        case(7) ! +x+y
            patchID = 8
        case default
            call abort(197001, "You ended up in the wrong neighborhood")
        end select
    else
        ! 3D
        select case(neighborhood)
        ! ---faces---
        case (5, 43:46) ! -x
            patchID = 1
        case (3, 35:38) ! +x
            patchID = 2
        case (2, 31:34) ! -y
            patchID = 3
        case (4, 39:42) ! +y
            patchID = 4
        case (6, 47:50) ! -z
            patchID = 5
        case (1, 27:30) ! +z
            patchID = 6
        ! ---(partial) edges---
        case (16, 69:70) ! -x-y
            patchID = 7
        case (18, 73:74) ! -x+y
            patchID = 8
        case (15, 67:68) ! +x-y
            patchID = 9
        case (17, 71:72) ! +x+y
            patchID = 10
        case (14, 65:66) ! -x  -z
            patchID = 11
        case (10, 57:58) ! -x  +z
            patchID = 12
        case (12, 61:62) ! +x  -z
            patchID = 13
        case (8,  53:54) ! +x  +z
            patchID = 14
        case (11, 59:60) !   -y-z
            patchID = 15
        case (7,  51:52) !   -y+z
            patchID = 16
        case (13, 63:64) !   +y-z
            patchID = 18
        case (9,  55:56) !   +y+z
            patchID = 19
        ! --- corners ---
        case (26) ! -x-y-z
            patchID = 20
        case (22) ! -x-y+z
            patchID = 21
        case (25) ! -x+y-z
            patchID = 22
        case (21) ! -x+y+z
            patchID = 23
        case (23) ! +x-y-z
            patchID = 24
        case (19) ! +x-y+z
            patchID = 25
        case (24) ! +x+y-z
            patchID = 26
        case (20) ! +x+y+z
            patchID = 14
        case default
            call abort(197001, "You ended up in the wrong neighborhood")
        end select
    endif
end 