! Fortran module for reading binary STL files. not suitable
! for ascii versions. adopted from Sukhbinder Singh
! https://sukhbinder.wordpress.com/2011/08/07/stl-files-and-fortran/
module module_stl_file_reader

    use module_globals
    use mpi

    implicit none

    ! I usually find it helpful to use the private keyword by itself initially, which specifies
    ! that everything within the module is private unless explicitly marked public.
    PRIVATE

    PUBLIC :: read_stl_file, normalize_stl_file, pointTriangleDistance, compute_superstl_file

    integer(kind=ik) :: mpirank, mpisize

contains

    ! compute all vertex- and edgenormals, store them in an ascii file
    ! See:
    ! Baerentzen, Aanaes. Generating Signed Distance Fields From Triangle Meshes (IMM-TECHNICAL REPORT-2002-21)
    subroutine compute_superstl_file( triangles, facenormals, ntri, superstl )
        implicit none
        integer, intent(in) :: ntri
        real(kind=4), intent(in) :: triangles(:,:), facenormals(:,:)
        real(kind=rk), allocatable, dimension(:,:), intent(out) :: superstl
        integer :: i, ivertex, last_info, info
        real(kind=rk), allocatable, dimension(:,:) :: points
        real(kind=rk), dimension(1:3) :: p1,p2,p3,n1,n2,n3,e1,e2,e3
        integer :: istart, iend, npercpu, mpirank, mpisize, mpicode

        allocate(superstl(1:ntri, 1:30))
        allocate(points(1:ntri, 1:9))

        ! fill points array
        do i = 1, ntri
            ivertex = 3*i - 2
            points(i, 1:3) = real( triangles(1:3,ivertex  ), kind=rk)
            points(i, 4:6) = real( triangles(1:3,ivertex+1), kind=rk)
            points(i, 7:9) = real( triangles(1:3,ivertex+2), kind=rk)
        enddo

        superstl = -90.e9_rk
        last_info = 0

        call MPI_Comm_rank(MPI_COMM_WORLD, mpirank, mpicode)
        call MPI_Comm_size(MPI_COMM_WORLD, mpisize, mpicode)

        npercpu = ntri / (mpisize)
        istart = mpirank*npercpu + 1
        iend = (mpirank+1)*npercpu

        if (mpirank == mpisize-1) then
            npercpu = ntri - (mpisize-1)*npercpu
            iend = ntri
        endif


        write(*,*) "Rank", mpirank, "triangles:", npercpu, istart, iend
        call MPI_BARRIER(MPI_COMM_WORLD, mpicode)

        ! compute all normals for all triangles.
        do i = istart, iend
            superstl(i,1:9) = points(i,:)
            superstl(i,10:12) = facenormals(1:3,i)

            p1 = points(i, 1:3)
            p2 = points(i, 4:6)
            p3 = points(i, 7:9)

            n1 = vertexnormal( p1, points, real(facenormals, kind=rk) )
            n2 = vertexnormal( p2, points, real(facenormals, kind=rk) )
            n3 = vertexnormal( p3, points, real(facenormals, kind=rk) )

            e1 = edgenormal( p1, p2, i, points, real(facenormals, kind=rk) )
            e2 = edgenormal( p2, p3, i, points, real(facenormals, kind=rk) )
            e3 = edgenormal( p3, p1, i, points, real(facenormals, kind=rk) )

            superstl(i, 13:15) = n1
            superstl(i, 16:18) = n2
            superstl(i, 19:21) = n3

            superstl(i, 22:24) = e1
            superstl(i, 25:27) = e2
            superstl(i, 28:30) = e3

            if (mpirank == mpisize-1) then
                info = int(100.0*real(i-istart) / real(npercpu))
                if (info /= last_info) then
                    last_info = info
                    write(*,'("Parallel processing.. ~",i3,"%")') info
                endif
            endif
        enddo

        if (mpirank==0) write(*,*) "processing part done, gathering up result"
        call MPI_ALLREDUCE(MPI_IN_PLACE, superstl, size(superstl), MPI_DOUBLE_PRECISION,&
         MPI_MAX, MPI_COMM_WORLD, mpicode)
    end subroutine


    function vertexnormal( x, points, facenormals )
        implicit none
        real(kind=rk), intent(in) :: points(:,:), facenormals(:,:), x(1:3)
        real(kind=rk), dimension(1:3) :: vertexnormal
        real(kind=rk), dimension(1:3) :: p1, p2, p3
        integer :: ntri, i,  k
        integer :: faces1(1:20), faces2(1:20), faces3(1:20)
        integer :: N1, N2, N3
        real(kind=rk) :: alpha

        ntri = size(points, 1)
        faces1 = 0
        faces2 = 0
        faces3 = 0
        N1 = 0
        N2 = 0
        N3 = 0
        vertexnormal = 0.0_rk

        do i = 1, ntri
            if (is_same_vector(x, points(i,1:3))) then
                ! faces containing the vertex as 1st point
                N1 = N1 + 1
                faces1(N1) = i
            endif

            if (is_same_vector(x, points(i,4:6))) then
                ! faces containing the vertex as 2nd point
                N2 = N2 + 1
                faces2(N2) = i
            endif

            if (is_same_vector(x, points(i,7:9))) then
                ! faces containing the vertex as 3rd point
                N3 = N3 + 1
                faces3(N3) = i
            endif
        enddo

        do k = 1, N1
            p1 = points( faces1(k), 1:3)
            p2 = points( faces1(k), 4:6)
            p3 = points( faces1(k), 7:9)

            alpha = angle_between_vectors( p2-p1, p3-p1)

            vertexnormal = vertexnormal + alpha * unit_vector(facenormals(1:3,faces1(k)))
        enddo

        do k = 1, N2
            p1 = points( faces2(k), 1:3)
            p2 = points( faces2(k), 4:6)
            p3 = points( faces2(k), 7:9)

            alpha = angle_between_vectors( p3-p2, p1-p2)

            vertexnormal = vertexnormal + alpha * unit_vector(facenormals(1:3,faces2(k)))
        enddo

        do k = 1, N3
            p1 = points( faces3(k), 1:3)
            p2 = points( faces3(k), 4:6)
            p3 = points( faces3(k), 7:9)

            alpha = angle_between_vectors( p1-p3, p2-p3)

            vertexnormal = vertexnormal + alpha * unit_vector(facenormals(1:3,faces3(k)))
        enddo

        if (N1+N2+N3 == 1) then !point only part of a single triangle?
            if (N1==1) then
                vertexnormal = unit_vector( -((p3-p1)+(p2-p1))  )
            elseif (N2==1) then
                vertexnormal = unit_vector(  -((p3-p2)+(p1-p2))   )
            elseif (N3==1) then
                vertexnormal = unit_vector(  -((p1-p3)+(p2-p3))   )
            endif
        else
            vertexnormal = unit_vector(vertexnormal)
        endif
    end function



    function edgenormal( p1, p2, itri, points, facenormals )
        implicit none
        real(kind=rk), intent(in) :: points(:,:), facenormals(:,:)
        real(kind=rk), intent(in) :: p1(1:3), p2(1:3)
        real(kind=rk), dimension(1:3) :: edgenormal
        real(kind=rk) :: n_tmp(1:3), p3(1:3)
        integer, intent(in) :: itri
        integer :: ntri, i
        logical :: found

        ntri = size(points, 1)
        ! given triangle A, we need to find triangle B
        found = .false.

        do i = 1, ntri
            if (i /= itri) then
                if (is_same_vector(p1, points(i,1:3)) .or. is_same_vector(p1, points(i,4:6)) .or. is_same_vector(p1, points(i,7:9))) then
                    ! p1 is on triangle
                    if (is_same_vector(p2, points(i,1:3)) .or. is_same_vector(p2, points(i,4:6)) .or. is_same_vector(p2, points(i,7:9))) then
                        ! p2 as well! found the second triangle
                        found = .true.
                        exit
                    endif
                endif
            endif
        enddo

        if (found) then
            edgenormal = unit_vector( unit_vector(facenormals(1:3,i)) + unit_vector(facenormals(1:3,itri)) )
        else
            ! edgenormal = unit_vector( facenormals(1:3,itri) )
            n_tmp = unit_vector( cross(p1-p2, facenormals(1:3,itri) ) )

            if ((.not. is_same_vector(p1,points(itri,1:3))) .and. (.not. is_same_vector(p2,points(itri,1:3)))) then
                p3 = points(itri, 1:3)
            endif
            if ((.not. is_same_vector(p1,points(itri,4:6))) .and. (.not. is_same_vector(p2,points(itri,4:6)))) then
                p3 = points(itri, 4:6)
            endif
            if ((.not. is_same_vector(p1,points(itri,7:9))) .and. (.not. is_same_vector(p2,points(itri,7:9)))) then
                p3 = points(itri, 7:9)
            endif

            if (dot_product(p3-p1, n_tmp) > 0.0 ) then
                n_tmp = unit_vector( cross(p2-p1, facenormals(1:3,itri) ) )
            endif

            edgenormal = n_tmp
        endif

    end function




    ! read a binary *.stl file from disk. the file contains
    ! ntri triangles, which are returned in "triangles". each
    ! triangle has 3 points and 1 normal, which is returned
    ! in "normals".
    subroutine read_stl_file(filename, ntri, triangles, normals)
        implicit none
        character(len=*), intent(in) :: filename
        integer(kind=4), intent(out) :: ntri
        ! dimensionality is (3,ntri)
        real(kind=4), allocatable, dimension(:,:) :: triangles, normals

        character(len=cshort) :: title
        integer :: iunit, irc, k, i, ierr
        integer*2 :: padding
        real*4 :: n(3),x1(3),x2(3),x3(3)

        call MPI_Comm_rank(MPI_COMM_WORLD, mpirank, ierr)
        call MPI_Comm_size(MPI_COMM_WORLD, mpisize, ierr)

        ! only root reads from file; other procs receive data later
        if (mpirank==0) then
            write(*,*) "--------------------------------------------------------"
            write(*,*) " ____ _____ _                             _           "
            write(*,*) "/ ___|_   _| |         _ __ ___  __ _  __| | ___ _ __ "
            write(*,*) "\___ \ | | | |   _____| '__/ _ \/ _` |/ _` |/ _ \ '__|"
            write(*,*) " ___) || | | |__|_____| | |  __/ (_| | (_| |  __/ |   "
            write(*,*) "|____/ |_| |_____|    |_|  \___|\__,_|\__,_|\___|_|   "
            write(*,*) "--------------------------------------------------------"
            write(*,*) "stlreader is reading file="//trim(adjustl(filename))

            iunit = 13
            irc = 0
            ntri = 0

#ifdef TURING
            ! TURING's IBM XL compiler by default reads binary data in BIG ENDIAN exclusively
            ! while it is supposedly possible to modify this behavior with environment variables
            ! I did not succeed in doing so. the following command ensures reading from unit 13
            ! in little endian.
            ! see also:
            ! http://www.ibm.com/support/knowledgecenter/SS2MB5_14.1.0/com.ibm.xlf141.bg.doc/language_ref/sup-setrteopts.html%23sup-setrteopts
            ! http://www-01.ibm.com/support/docview.wss?uid=swg21243120
            call setrteopts("ufmt_littleendian=13")
#endif

            iunit=13
            open(unit=iunit,file=filename,status='old',access='stream', form='unformatted')

            read(iunit) title
            read(iunit) ntri

            write(*,*) "file header: ", trim(adjustl(title))
            write(*,*) "number of triangles: ", ntri

            allocate(normals(3,ntri))
            allocate(triangles(3,ntri*3))

            k=1
            do i = 1,ntri
                read(iunit) normals(1,i),normals(2,i),normals(3,i)
                read(iunit) triangles(1,k),triangles(2,k),triangles(3,k)
                read(iunit) triangles(1,k+1),triangles(2,k+1),triangles(3,k+1)
                read(iunit) triangles(1,k+2),triangles(2,k+2),triangles(3,k+2)
                read(iunit) padding
                k=k+3
            end do
            close(iunit)


            ! print some diagnostics for the file we just read.
            write(*,*) "----------------"
            write(*,*) "file just read: "//trim(adjustl(filename))
            write(*,*) 'STL-header ',trim(title),' and has ',ntri, ' triangles'
            write(*,*) "dimensions (not normalized, as read from file):"
            write(*,'("Lx=",g15.6)') maxval(triangles(1,:))-minval(triangles(1,:))
            write(*,'("Ly=",g15.6)') maxval(triangles(2,:))-minval(triangles(2,:))
            write(*,'("Lz=",g15.6)') maxval(triangles(3,:))-minval(triangles(3,:))

            write(*,'("xmin=",g15.6," xmax=",g15.6)') minval(triangles(1,:)), maxval(triangles(1,:))
            write(*,'("ymin=",g15.6," ymax=",g15.6)') minval(triangles(2,:)), maxval(triangles(2,:))
            write(*,'("zmin=",g15.6," zmax=",g15.6)') minval(triangles(3,:)), maxval(triangles(3,:))

        endif

        ! only root read from file: now, distribute data to all other procs
        call broadcast_stl_file(ntri, triangles, normals)
    end subroutine


    ! distribute *.stl to all procs, from root. avoids that all cpu have to perform
    ! i/o
    subroutine broadcast_stl_file(ntri, triangles, normals)
        use mpi
        implicit none
        integer, intent(inout) :: ntri
        ! dimensionality is (3,ntri)
        real(kind=4), allocatable, dimension(:,:) :: triangles, normals
        integer :: mpicode, i

        ! check if on root, array is allocated, if not, yell
        if (.not. allocated(triangles) .and. mpirank==0) then
            call abort(2020,"broadcast stl file but file not yet read?")
        endif
        ! broadcast number of triangles and allocate memory on ranks other
        ! than root
        call MPI_BCAST( ntri, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, mpicode )
        if (mpirank/=0) then
            allocate(normals(3,ntri))
            allocate(triangles(3,ntri*3))
        endif

        do i=1,3
            call MPI_BCAST( normals(1,:), ntri, MPI_REAL, 0, MPI_COMM_WORLD, mpicode )
            call MPI_BCAST( normals(2,:), ntri, MPI_REAL, 0, MPI_COMM_WORLD, mpicode )
            call MPI_BCAST( normals(3,:), ntri, MPI_REAL, 0, MPI_COMM_WORLD, mpicode )

            call MPI_BCAST( triangles(1,:), 3*ntri, MPI_REAL, 0, MPI_COMM_WORLD, mpicode )
            call MPI_BCAST( triangles(2,:), 3*ntri, MPI_REAL, 0, MPI_COMM_WORLD, mpicode )
            call MPI_BCAST( triangles(3,:), 3*ntri, MPI_REAL, 0, MPI_COMM_WORLD, mpicode )
        enddo

        ! security barrier: if for some reason not all ranks call MPI_BCAST approx. at the same time,
        ! that causes crashes on some machines.
        call MPI_BARRIER(MPI_COMM_WORLD, mpicode)
    end subroutine


    ! stl files contain non-normalized triangles (i.e. without units). this is not
    ! perturbing since we normalize our code anyways, but it is somewhat tricky to
    ! find an appropriate lengthscale.
    subroutine normalize_stl_file( ntri, triangles, mode, scale, origin )
        implicit none
        integer, intent(in) :: ntri
        real(kind=4), dimension(1:3,3*ntri) :: triangles
        real(kind=rk), dimension(3), intent(in) :: origin
        character(len=*), intent(in) :: mode
        real(kind=rk), intent(inout) :: scale

        select case(mode)
        case ("--lx")
            scale = maxval(triangles(1,:))-minval(triangles(1,:))

        case ("--ly")
            scale = maxval(triangles(2,:))-minval(triangles(2,:))

        case ("--lz")
            scale = maxval(triangles(3,:))-minval(triangles(3,:))

        case ("--scale")
            ! do nothing, scale is given directly

        case default
            call abort(270618,"normalize_stl_file: mode unkown.")

        end select

        if (mpirank==0) then
            write(*,'("normalizing stl-file: scale=",g15.6, "x0=",3(g12.4,1x))') scale, origin
        endif

        ! normalize distances
        triangles = triangles / scale
        ! shift origin to x0,y0,z0
        triangles(1,:) = triangles(1,:) + origin(1)
        triangles(2,:) = triangles(2,:) + origin(2)
        triangles(3,:) = triangles(3,:) + origin(3)

        if (mpirank ==0 ) then
            write(*,*) "dimensions (normalized):"
            write(*,'("Lx=",g15.6)') maxval(triangles(1,:))-minval(triangles(1,:))
            write(*,'("Ly=",g15.6)') maxval(triangles(2,:))-minval(triangles(2,:))
            write(*,'("Lz=",g15.6)') maxval(triangles(3,:))-minval(triangles(3,:))

            write(*,'("xmin=",g15.6," xmax=",g15.6)') minval(triangles(1,:)), maxval(triangles(1,:))
            write(*,'("ymin=",g15.6," ymax=",g15.6)') minval(triangles(2,:)), maxval(triangles(2,:))
            write(*,'("zmin=",g15.6," zmax=",g15.6)') minval(triangles(3,:)), maxval(triangles(3,:))
        endif
    end subroutine normalize_stl_file

    ! same function, but returns sign using face-, vertex- and edgenormals
    real(kind=rk) function pointTriangleDistance(P1,P2,P3,point,normal, n1,n2,n3, e1,e2,e3)
        ! calculate distance between a point and a triangle in 3D
        ! SYNTAX
        !   dist = pointTriangleDistance(TRI,P)
        !   [dist,PP0] = pointTriangleDistance(TRI,P)
        !
        ! DESCRIPTION
        !   Calculate the distance of a given point P from a triangle TRI.
        !   Point P is a row vector of the form 1x3. The triangle is a matrix
        !   formed by three rows of points TRI = [P1P2P3] each of size 1x3.
        !   dist = pointTriangleDistance(TRI,P) returns the distance of the point P
        !   to the triangle TRI.
        !   [dist,PP0] = pointTriangleDistance(TRI,P) additionally returns the
        !   closest point PP0 to P on the triangle TRI.
        !
        ! Author: Gwendolyn Fischer
        ! Release: 1.0
        ! Release date: 09/02/02
        ! Release: 1.1 Fixed Bug because of normalization
        ! Release: 1.2 Fixed Bug because of typo in region 5 20101013
        ! Release: 1.3 Fixed Bug because of typo in region 2 20101014

        ! Possible extention could be a version tailored not to return the distance
        ! and additionally the closest point, but instead return only the closest
        ! point. Could lead to a small speed gain.

        ! Example:
        ! !! The Problem
        ! P0 = [0.5 -0.3 0.5]
        !
        ! P1 = [0 -1 0]
        ! P2 = [1  0 0]
        ! P3 = [0  0 0]
        !
        ! vertices = [P1 P2 P3]
        ! faces = [1 2 3]
        !
        ! !! The Engine
        ! [dist,PP0] = pointTriangleDistance([P1P2P3],P0)
        !
        ! !! Visualization
        ! [x,y,z] = sphere(20)
        ! x = dist*x+P0(1)
        ! y = dist*y+P0(2)
        ! z = dist*z+P0(3)
        !
        ! figure
        ! hold all
        ! patch('Vertices',vertices,'Faces',faces,'FaceColor','r','FaceAlpha',0.8)
        ! plot3(P0(1),P0(2),P0(3),'b*')
        ! plot3(PP0(1),PP0(2),PP0(3),'*g')
        ! surf(x,y,z,'FaceColor','b','FaceAlpha',0.3)
        ! view(3)

        ! The algorithm is based on
        ! "David Eberly, 'Distance Between Point and Triangle in 3D',
        ! Geometric Tools, LLC, (1999)"
        ! http:\\www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf
        !
        !        ^t
        !  \     |
        !   \reg2|
        !    \   |
        !     \  |
        !      \ |
        !       \|
        !        *P2
        !        |\
        !        | \
        !  reg3  |  \ reg1
        !        |   \
        !        |reg0\
        !        |     \
        !        |      \ P1
        ! -------*-------*------->s
        !        |P0      \
        !  reg4  | reg5    \ reg6
        implicit none
        real(kind=rk), dimension(1:3), intent(in) :: P1,P2,P3,normal, n1,n2,n3,e1,e2,e3
        real(kind=rk), dimension(1:3), intent(in) :: point
        real(kind=rk), dimension(1:3) :: BB,EE0,EE1,DD
        real(kind=rk) :: a,b,c,d,e,f,det,s,t,sqrDistance,tmp0,tmp1,numer,denom,invDet
        integer(kind=ik) :: region

        ! rewrite triangle in normal form
        BB = P1
        EE0 = P2-BB
        EE1 = P3-BB


        DD = BB - point
        a = dot_product(EE0,EE0)
        b = dot_product(EE0,EE1)
        c = dot_product(EE1,EE1)
        d = dot_product(EE0,DD)
        e = dot_product(EE1,DD)
        f = dot_product(DD,DD)



        det = a*c - b*b ! do we have to use abs here?
        s   = b*e - c*d
        t   = b*d - a*e

        ! if (abs(det) < 1.0d-16) then
        !   pointTriangleDistance = 9.0d9
        !   return
        ! endif

        region = -1

        ! Terible tree of conditionals to determine in which region of the diagram
        ! shown above the projection of the point into the triangle-plane lies.
        if ((s+t) <= det) then
            if (s < 0.0_rk) then
                    !-------------------------------------------------------
                    if (t < 0.0_rk) then
                        !region4 (vertex)
                        region = 4
                        if (d < 0.0_rk) then
                            t = 0.0_rk
                            if (-d >= a) then
                                s = 1.d0
                                sqrDistance = a + 2.0_rk*d + f
                            else
                                s = -d/a
                                sqrDistance = d*s + f
                            endif
                        else
                            s = 0.0_rk
                            if (e >= 0.0_rk) then
                                t = 0.0_rk
                                sqrDistance = f
                            else
                                if (-e >= c) then
                                    t = 1.d0
                                    sqrDistance = c + 2.0_rk*e + f
                                else
                                    t = -e/c
                                    sqrDistance = e*t + f
                                endif
                            endif
                        endif !of region 4
                    !-------------------------------------------------------
                    else
                    ! region 3 (edge)
                    region = 3
                    s = 0.0_rk
                    if (e >= 0.0_rk) then
                        t = 0.0_rk
                        sqrDistance = f
                    else
                        if (-e >= c) then
                            t = 1.d0
                            sqrDistance = c + 2.0_rk*e +f
                        else
                            t = -e/c
                            sqrDistance = e*t + f
                        endif
                    endif
                endif !of region 3
                !-------------------------------------------------------
            else
                if (t < 0.0_rk) then
                    ! region 5
                    region = 5
                    t = 0.0_rk
                    if (d >= 0.0_rk) then
                        s = 0.0_rk
                        sqrDistance = f
                    else
                        if (-d >= a) then
                            s = 1.d0
                            sqrDistance = a + 2.0_rk*d + f! GF 20101013 fixed typo d*s ->2*d
                        else
                            s = -d/a
                            sqrDistance = d*s + f
                        endif
                    endif
                else
                    ! region 0
                    region = 0
                    invDet = 1.d0/det
                    s = s*invDet
                    t = t*invDet
                    sqrDistance = s*(a*s + b*t + 2.0_rk*d) &
                    + t*(b*s + c*t + 2.0_rk*e) + f
                endif
            endif
        else
            if (s < 0.0_rk) then
                ! region 2
                region = 2
                tmp0 = b + d
                tmp1 = c + e
                if (tmp1 > tmp0) then ! minimum on edge s+t=1
                    numer = tmp1 - tmp0
                    denom = a - 2.0_rk*b + c
                    if (numer >= denom) then
                        s = 1.d0
                        t = 0.0_rk
                        sqrDistance = a + 2.0_rk*d + f ! GF 20101014 fixed typo 2*b -> 2*d
                    else
                        s = numer/denom
                        t = 1.d0-s
                        sqrDistance = s*(a*s + b*t + 2.0_rk*d) &
                        + t*(b*s + c*t + 2.0_rk*e) + f
                    endif
                else          ! minimum on edge s=0
                    s = 0.0_rk
                    if (tmp1 <= 0.0_rk) then
                        t = 1.d0
                        sqrDistance = c + 2.0_rk*e + f
                    else
                        if (e >= 0.0_rk) then
                            t = 0.0_rk
                            sqrDistance = f
                        else
                            t = -e/c
                            sqrDistance = e*t + f
                        endif
                    endif
                endif !of region 2
            else
                if (t < 0.0_rk) then
                    region = 6
                    !region6
                    tmp0 = b + e
                    tmp1 = a + d
                    if (tmp1 > tmp0) then
                        numer = tmp1 - tmp0
                        denom = a-2.0_rk*b+c
                        if (numer >= denom) then
                            t = 1.d0
                            s = 0.0_rk
                            sqrDistance = c + 2.0_rk*e + f
                        else
                            t = numer/denom
                            s = 1.d0 - t
                            sqrDistance = s*(a*s + b*t + 2.0_rk*d) &
                            + t*(b*s + c*t + 2.0_rk*e) + f
                        endif
                    else
                        t = 0.0_rk
                        if (tmp1 <= 0) then
                            s = 1.d0
                            sqrDistance = a + 2.0_rk*d + f
                        else
                            if (d >= 0.0_rk) then
                                s = 0.0_rk
                                sqrDistance = f
                            else
                                s = -d/a
                                sqrDistance = d*s + f
                            endif
                        endif
                    endif
                    !end region 6
                else
                    ! region 1
                    region = 1
                    numer = c + e - b - d
                    if (numer <= 0.0_rk) then
                        s = 0.0_rk
                        t = 1.d0
                        sqrDistance = c + 2.0_rk*e + f
                    else
                        denom = a - 2.0_rk*b + c
                        if (numer >= denom) then
                            s = 1.d0
                            t = 0.0_rk
                            sqrDistance = a + 2.0_rk*d + f
                        else
                            s = numer/denom
                            t = 1-s
                            sqrDistance = s*(a*s + b*t + 2.0_rk*d) &
                            + t*(b*s + c*t + 2.0_rk*e) + f
                        endif
                    endif !of region 1
                endif
            endif
        endif

        ! account for numerical round-off error
        if (sqrDistance < 0.0_rk) then
            sqrDistance = 0.0_rk
        endif



        ! closest point on triangle
        DD = BB + s*EE0 + t*EE1;

        if (is_same_vector(DD, P1)) then
            ! vector from target point to closest point on surface
            DD = point-DD
            t = dot_product(DD, n1)
            region = 11

        elseif (is_same_vector(DD, P2)) then
            ! vector from target point to closest point on surface
            DD = point-DD
            t = dot_product(DD,n2)
            region = 22

        elseif (is_same_vector(DD, P3)) then
            ! vector from target point to closest point on surface
            DD = point-DD
            t = dot_product(DD, n3)
            region = 33
        else
            ! vector from target point to closest point on surface
            DD = point-DD

            select case (region)
            case (0)
                t = dot_product(DD,normal)
            case (1)
                t = dot_product(DD,e2)
            case (3)
                t = dot_product(DD,e3)
            case (5)
                t = dot_product(DD,e1)
            case (2)
                t = dot_product(DD,n3)
            case (4)
                t = dot_product(DD,n1)
            case (6)
                t = dot_product(DD,n2)
            end select
        endif


        if (t >= 0.0_rk) then
            pointTriangleDistance = dsqrt(sqrDistance)
        else
            pointTriangleDistance = -dsqrt(sqrDistance)
        endif

    end function


    subroutine readbin(iunit,a,irc)
        implicit none
        integer*4 ntri,iunit,irc
        real*4 a,b
        integer*2 ib(2)
        equivalence(b,ib)

        read(iunit,rec=irc)ib(1)
        irc=irc+1
        read(iunit,rec=irc)ib(2)
        irc=irc+1

        a=b

    end subroutine

end module module_stl_file_reader
