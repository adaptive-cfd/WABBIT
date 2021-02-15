subroutine draw_fractal_tree(Insect, xx0, ddx, mask, mask_color, us)
    implicit none

    type(diptera),intent(inout) :: Insect
    real(kind=rk),intent(in)    :: xx0(1:3), ddx(1:3)
    real(kind=rk),intent(inout) :: mask(0:,0:,0:)
    real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
    integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)

    real(kind=rk) :: x1(1:3), x2(1:3), R
    integer ::  i

    ! check if the global treedata is already filled = initialization is done
    if (.not.allocated(treedata)) then
        call abort(20200504, "draw_fractal_tree called but fractal_tree_init was not done!")
    endif

    ! 28/01/2019: Thomas. Discovered that this was done block based, i.e. the smoothing layer
    ! had different thickness, if some blocks happened to be at different levels (and still carry
    ! a part of the smoothing layer.) I don't know if that made sense, because the layer shrinks/expands then
    ! and because it might be discontinous. Both options are included now, default is "as before"
    ! Insect%smoothing_thickness=="local"  : smoothing_layer = c_sm * 2**-J * L/(BS-1)
    ! Insect%smoothing_thickness=="global" : smoothing_layer = c_sm * 2**-Jmax * L/(BS-1)
    if (Insect%smoothing_thickness == "local") then
        Insect%smooth = 1.50d0*maxval(ddx)
        Insect%safety = 3.5d0*Insect%smooth
    endif

    !*****************************************************************************
    ! phase 3: all ranks draw the individual cylinders...
    ! please note we assume the root point to be 0 0 0 in the file
    !*****************************************************************************
     do i = 1, size(treedata,1)
        x1 = (treedata(i,1:3) - treedata(1,1:3))*Insect%fractal_tree_scaling + Insect%fractal_tree_x0
        x2 = (treedata(i,4:6) - treedata(1,1:3))*Insect%fractal_tree_scaling + Insect%fractal_tree_x0

        ! the file containes the radius of the cylinder
        R = treedata(i,7)*Insect%fractal_tree_scaling

        call draw_cylinder_new( x1, x2, R, xx0, ddx, mask, mask_color, us, Insect, int(1,kind=2), bounding_box=treedata_boundingbox(i,1:6))
    end do

end subroutine draw_fractal_tree


subroutine fractal_tree_init(Insect)
    implicit none
    type(diptera), intent(inout) :: Insect
    integer :: nlines, i, Nphi, icyl
    real(kind=rk),dimension(1:3) :: e_x, e_r, e_tmp, x1, x2, p2, p1, R
    real(kind=rk),dimension(1:3,1:3) :: M_phi
    real(kind=rk) :: phi, xmin,ymin,zmin,xmax,ymax,zmax

    ! note: smoothing is set in insect_init (which must be called also for fractal trees)

    ! initialization
    if (.not. allocated(treedata)) then
        call check_file_exists(Insect%fractal_tree_file)

        !*****************************************************************************
        ! phase one: read the number of lines (which is the number of rigid cylinders in the tree)
        !*****************************************************************************
        call count_lines_in_ascii_file_mpi(Insect%fractal_tree_file, nlines, 0)

        if (root) write(*,'(80("-"))')
        if (root) write(*,'("Building a fractal tree with ",i5," rigid cylinders")') nlines
        if (root) write(*,'("Fractal tree scaling factor is ",g16.4)') Insect%fractal_tree_scaling
        if (root) write(*,'("Fractal tree scaling root point is ",3(g16.4,1x))') Insect%fractal_tree_x0

        !*****************************************************************************
        ! phase two: read all cylinders from the file into the array and bcast them
        !*****************************************************************************
        allocate( treedata(1:nlines, 1:7) )
        call read_array_from_ascii_file_mpi(Insect%fractal_tree_file, treedata, n_header=0)

        if (root) write(*,'("Done reading ",i5," rigid cylinders")') nlines

        !*****************************************************************************
        ! phase 3: compute bounding boxes of each cylinder (otherwise has to be done for every block)
        !*****************************************************************************
        if (root) write(*,'("Computing static bounding boxes for cylinders in xyz")')
        allocate( treedata_boundingbox(1:nlines, 1:6) )

        do icyl = 1, nlines
            ! cylinder start and end point
            x1 = (treedata(icyl,1:3) - treedata(1,1:3))*Insect%fractal_tree_scaling + Insect%fractal_tree_x0
            x2 = (treedata(icyl,4:6) - treedata(1,1:3))*Insect%fractal_tree_scaling + Insect%fractal_tree_x0

            ! the file containes the radius of the cylinder
            R = treedata(icyl,7)*Insect%fractal_tree_scaling

            ! unit vector in cylinder axis direction
            e_x = x2 - x1
            e_x = e_x / norm2(e_x)

            ! radial unit vector
            ! use a vector perpendicular to e_x, since it is a azimuthal symmetry
            ! it does not really matter which one. however, we must be sure that the vector
            ! we use and the e_x vector are not colinear -- their cross product is the zero vector, if that is the case
            e_r = (/0.d0, 0.d0, 0.d0/)
            do while ( norm2(e_r) <= 1.0d-12 )
                e_r = cross( (/rand_nbr(),rand_nbr(),rand_nbr()/), e_x)
            enddo
            e_r = e_r / norm2(e_r)

            ! initialize bounding box
            xmin = +9999999.0_rk
            ymin = +9999999.0_rk
            zmin = +9999999.0_rk
            xmax = 0.0_rk
            ymax = 0.0_rk
            zmax = 0.0_rk

            ! actual bounding box computation
            ! note for cylinders that are arbitrarily aligned, checking 4 points is not enough
            ! to accurately determine the bounding box
            Nphi = 50
            do i = 0, Nphi
                phi = real(i, kind=rk)/real(Nphi+1, kind=rk) * 2.0_rk * pi
                call Rarb(M_phi, phi, e_x)

                e_tmp = matmul(M_phi, e_r)

                p1 = x1 + R*e_tmp
                p2 = x2 + R*e_tmp

                xmin = min(xmin, minval((/p1(1), p2(1)/)) )
                ymin = min(ymin, minval((/p1(2), p2(2)/)) )
                zmin = min(zmin, minval((/p1(3), p2(3)/)) )

                xmax = max(xmax, maxval((/p1(1), p2(1)/)) )
                ymax = max(ymax, maxval((/p1(2), p2(2)/)) )
                zmax = max(zmax, maxval((/p1(3), p2(3)/)) )
            enddo

            ! save result in array, to be used in draw_cylinder routine
            treedata_boundingbox(icyl, 1:6) = (/xmin, ymin, zmin, xmax, ymax, zmax/)
        enddo

        if (root) write(*,'("Computing static bounding done.")')

    else
        call abort(20200504, "fractal_tree_init seems to be called twice.")
    endif

end subroutine
