
module module_trimesh

    use module_globals

    implicit none

    type surfaceTrimesh
        ! Store information about the surface trimesh.
        ! MUST be manifold, MUST be a surface (closed or open)
        real(kind=rk), allocatable :: vertices(:, :)
        integer, allocatable :: faces(:, :)
        integer, allocatable :: edges(:, :)
        integer, allocatable :: neighbours(:, :)

    end type surfaceTrimesh

contains

    subroutine importGnuplotMesh(filename, trimesh)

        Use, intrinsic :: iso_fortran_env, Only : iostat_end
        
        character(len=*), INTENT(IN) :: filename
        type(surfaceTrimesh), INTENT(OUT) :: trimesh

        real(kind=rk) :: datap1(4, 2)
        integer :: error, iunit, nedges, nelements, nverts
        real(kind=rk) :: x_temp, y_temp
        

        iunit=2
        nedges=0
        nverts = 0
        nelements = 0

        ! read very simple Gnuplot mesh format. 
        ! This consists of a set of triangles, defined by four points (p1 -> p2 -> p3 -> p1).
        ! First we count the lines (how do I do this cleverly?)


        OPEN(unit=iunit, file=filename)

        do 
            read(iunit, *, iostat=error) datap1(1, 1), datap1(1, 2), datap1(2, 1), datap1(2, 2), datap1(3, 1), datap1(3, 2), datap1(4, 1), datap1(4, 2)

            select case(error)
                case(0)
                case(iostat_end)
                    EXIT
                case default
                    write(*, *)
            end select     
            nelements = nelements + 1
            nverts = nverts + 3
            nedges = nedges + 3       
        end do
        write(*, *) 'total ', nelements, ' elements;' 
        write(*, *) 'total ', nverts, ' vertices;' 
        write(*, *) 'total ', nedges, ' edges loaded' 
        CLOSE(iunit)


        ! now we allocate the size of the array and read in the points

        allocate(trimesh%faces(nelements, 3))
        allocate(trimesh%edges(nedges, 2))
        allocate(trimesh%vertices(nverts, 2))

        OPEN(unit=iunit, file=filename)

        nedges=0
        nverts = 0
        nelements = 0

        do 
            read(iunit, *, iostat=error) datap1(1, 1), datap1(1, 2), datap1(2, 1), datap1(2, 2), datap1(3, 1), datap1(3, 2), datap1(4, 1), datap1(4, 2)
            
            select case(error)
                case(0)
                case(iostat_end)
                    EXIT
                case default
                    write(*, *)
            end select
            
            nelements = nelements + 1
            nverts = nverts + 3
            nedges = nedges + 3

            trimesh%vertices(nverts - 2, :) = (/datap1(1, 1), datap1(1, 2)/)
            trimesh%vertices(nverts - 1, :) = (/datap1(2, 1), datap1(2, 2)/)
            trimesh%vertices(nverts, :) = (/datap1(3, 1), datap1(3, 2)/)

            trimesh%edges(nedges - 2, :) = (/nverts - 2, nverts - 1/)
            trimesh%edges(nedges - 1, :) = (/nverts - 1, nverts/)
            trimesh%edges(nedges, :) = (/nverts, nverts - 2/)

            trimesh%faces(nelements, :) = (/nedges - 2, nedges - 1, nedges/)         
            
        end do
        CLOSE(iunit)
        write(*, *) "mesh loaded!"

    end subroutine importGnuplotMesh

    subroutine mergeNeighbours(trimesh, tolerance)
        type(surfaceTrimesh), INTENT(INOUT) :: trimesh
        real(kind=rk), INTENT(IN) :: tolerance

        integer :: i, j, k, n, nelements, nverts, nedges, maxNoNeighbours

        nelements = size(trimesh%faces, 1)
        nverts = size(trimesh%vertices, 1)
        nedges = size(trimesh%edges, 1)

        write(*, *) trimesh%vertices(trimesh%faces(nelements, :), :)

    end subroutine mergeNeighbours


end module module_trimesh
