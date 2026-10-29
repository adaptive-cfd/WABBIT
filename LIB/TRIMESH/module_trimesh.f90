
module module_trimesh

    use module_globals
    use module_stl_file_reader

    implicit none

    type deformable_superstl
        ! Store information about the surface trimesh. by importing a superstl. Note that we store two versions, the original, providing the "material coordinates" and the deformable version

        ! vertex1        = Insect%body_superSTL_g(i, 1:3)
        ! vertex2        = Insect%body_superSTL_g(i, 4:6)
        ! vertex3        = Insect%body_superSTL_g(i, 7:9)
        ! face_normal    = Insect%body_superSTL_g(i, 10:12)
        ! vertex1_normal = Insect%body_superSTL_g(i, 13:15)
        ! vertex2_normal = Insect%body_superSTL_g(i, 16:18)
        ! vertex3_normal = Insect%body_superSTL_g(i, 19:21)
        ! edge1_normal   = Insect%body_superSTL_g(i, 22:24)
        ! edge2_normal   = Insect%body_superSTL_g(i, 25:27)
        ! edge3_normal   = Insect%body_superSTL_g(i, 28:30)

        real(kind=rk), allocatable, dimension(:,:) :: superstl_mat, superstl_w, verts
        integer, allocatable, dimension(:, :) :: edges   ! n_tri x 3
        integer, allocatable, dimension(:) :: edge_adjacency   ! n_edges x 1
        integer, ALLOCATABLE, DIMENSION(:, :) :: global_verts ! n_verts x RANK(superstl) = n_verts x 2
        integer, ALLOCATABLE, DIMENSION(:, :) :: shared_points ! n_unique x (max. no. of twins)
        integer :: ntri, n_unique
        real(kind=rk) :: scale, origin(3), tolerance
        character(len=cshort) :: fname

    end type deformable_superstl

contains



    subroutine generate_edges_and_points(wing_defSTL)

        type(deformable_superstl), INTENT(inout) :: wing_defSTL

        real(kind=rk) :: scale, origin(3), E1_temp(2, 3), E2_temp(2, 3), p2_temp, norm(2), tolerance
        integer :: i, j, k, m, o, ntri, n_edge, n_uniquepoints, n_adjacentedges, n_adjacentpoints_temp, n_adjacentedges_temp, temp_edge(4), current_edge(4)

        scale = 1.0_rk
        origin = (/0.0_rk, 0.0_rk, 0.0_rk/)
        tolerance=1e-9        
        
        wing_defSTL%superstl_w = wing_defSTL%superstl_mat

        wing_defSTL%origin = origin
        wing_defSTL%scale = scale
        wing_defSTL%tolerance = 1.0e-8_rk

        ntri = size(wing_defSTL%superstl_mat, 1)

        wing_defSTL%ntri = ntri
        ! create an edgelist: n_edges x (RANK(SUPERSTL) * 2) = n_edges x 4
        allocate(wing_defSTL%edges(ntri * 3, 2 * RANK(wing_defSTL%superstl_mat)))
        ! create a vertice list
        ALLOCATE(wing_defSTL%global_verts(ntri * 3, 2))
        ALLOCATE(wing_defSTL%verts(ntri * 3, 2))

        do i = 1, ntri
            wing_defSTL%edges(i * 3 - 2, 1:2) = (/ i, 1/)
            wing_defSTL%edges(i * 3 - 2, 3:4) = (/ i, 4/)

            wing_defSTL%edges(i * 3 - 1, 1:2) = (/ i, 4/)
            wing_defSTL%edges(i * 3 - 1, 3:4) = (/i, 7/)

            wing_defSTL%edges(i * 3, 1:2) = (/i, 7/)
            wing_defSTL%edges(i * 3, 3:4) = (/i, 1/)

            wing_defSTL%global_verts(i * 3 - 2, :) = (/i, 1/)
            wing_defSTL%global_verts(i * 3 - 1, :) = (/i, 4/)
            wing_defSTL%global_verts(i * 3, :) = (/i, 7/)


            wing_defSTL%verts(i * 3 - 2, :) = wing_defSTL%superstl_mat(i, 1:3)
            wing_defSTL%verts(i * 3 - 1, :) = wing_defSTL%superstl_mat(i, 4:6)
            wing_defSTL%verts(i * 3, :) = wing_defSTL%superstl_mat(i, 7:9)
        
        end do

        ! create a list of edge adjacencies: n_edges x 1
        ALLOCATE(wing_defSTL%edge_adjacency(ntri * 3)) 


        wing_defSTL%edge_adjacency=0

        wing_defSTL%n_unique = ntri * 3



        end subroutine generate_edges_and_points

    subroutine find_adjacent_edges(wing_defSTL)

        type(deformable_superstl), INTENT(out) :: wing_defSTL
        real(kind=4), allocatable, dimension(:,:) :: triangles, normals, tmp_verts
        integer :: i, j, k, m, o, ntri, n_edge, n_unique_temp, n_adjacentedges, n_adjacentpoints_temp, n_adjacentedges_temp, temp_edge(4), current_edge(4), P1_idx(2), P2_idx(2)
        real(kind=rk) :: scale, origin(3), E1_temp(2, 3), E2_temp(2, 3), norm(2), P1(3), P2(3), tolerance, norm_1

        scale = 1.0_rk
        origin = (/0.0_rk, 0.0_rk, 0.0_rk/)
        tolerance=1e-9

        ntri = size(wing_defSTL%superstl_mat, 1)
        ! loop through the edge list. check for adjacent edges, in a 2d manifold surface there can only be one other bordering edge per individual edge.

        n_unique_temp = wing_defSTL%n_unique
        wing_defSTL%edge_adjacency(i) = 0

        do i = 1, ntri * 3
            current_edge = wing_defSTL%edges(i, :)

            E1_temp(1, :) = wing_defSTL%superstl_mat(current_edge(1), current_edge(2):(current_edge(2)+2))
            E1_temp(2, :) = wing_defSTL%superstl_mat(current_edge(3), current_edge(4):(current_edge(4) + 2))

            do j = i + 1, ntri * 3
                temp_edge =wing_defSTL%edges(j, :)
                E2_temp(1, :) = wing_defSTL%superstl_mat(temp_edge(1), temp_edge(2):(temp_edge(2)+2))
                E2_temp(2, :) = wing_defSTL%superstl_mat(temp_edge(3), temp_edge(4):(temp_edge(4) + 2))

                norm(1)  = sum((E1_temp(1, :)-E2_temp(1, :)) ** 2)
                norm(2)  = sum((E1_temp(2, :)-E2_temp(2, :)) ** 2)

                if (count(norm < tolerance ** 2) == 2) then
                    wing_defSTL%edge_adjacency(i) = j
                    wing_defSTL%edge_adjacency(j) = i
                    cycle
                end if

                norm(1)  = sum((E1_temp(1, :)-E2_temp(2, :)) ** 2)
                norm(2)  = sum((E1_temp(2, :)-E2_temp(1, :)) ** 2)


                if (count(norm < tolerance ** 2) == 2) then
                    wing_defSTL%edge_adjacency(i) = j
                    wing_defSTL%edge_adjacency(j) = i
                    exit
                end if
            end do
        end do

        ALLOCATE(wing_defSTL%shared_points(n_unique_temp, 8))


    end subroutine find_adjacent_edges


        
        ! ! loop through all vertices, finding identical points (within a tolerance) to count the maximum amount of overlapping points.
        ! ! for a 2D surface without any overlapping triangles, we can only have a maximum of ONE other adjacent edge for each edge, or three adjacent edges per face! This makes it easier to map adjacency in a (n_edges x 2) MAtrix,
        ! ! 

        ! n_adjacentedges = 0
        ! n_adjacentpoints = 0

        !do i = 1, ntri

        !     ! cycle through each face in the list
        !     F1_temp(1, 1:3) = wing_defSTL%superstl_mat(i, 1:3)
        !     F1_temp(2, 1:3) = wing_defSTL%superstl_mat(i, 4:6)
        !     F1_temp(3, 1:3) = wing_defSTL%superstl_mat(i, 7:9)

        !     n_adjacentedges_temp = 0
        !     n_adjacentpoints_temp = 0

        !     do j = 1, i-1

        !         ! cycle through each already visited face

        !         F2_temp(1, 1:3) = wing_defSTL%superstl_mat(j, 1:3)
        !         F2_temp(2, 1:3) = wing_defSTL%superstl_mat(j, 4:6)
        !         F2_temp(3, 1:3) = wing_defSTL%superstl_mat(j, 7:9)

        !         ! calculate distance to current point 
        !         norm(1) = sum((F1_temp(1, 1:3) - F2_temp(1, 1:3)) ** 2)
        !         norm(2) = sum((F1_temp(2, 1:3) - F2_temp(2, 1:3)) ** 2)
        !         norm(3) = sum((F1_temp(3, 1:3) - F2_temp(3, 1:3)) ** 2)

        !         ! if all points are the same, the triangle is the same (should not be the case)
        !         ! if two points are the same, we have an adjacent edge and two indentical points. 
        !         ! if one point is the same, we have no adjcacent edges and one identical point.

        !         n_adjacentedges_temp = count(norm < tolerance)

        !         if (n_adjacentedges_temp == 3 ) then
        !             call abort(2345, "bad STL file gave me overlapping triangles!")
        !         end if



        !     end do


        !end do


end module module_trimesh