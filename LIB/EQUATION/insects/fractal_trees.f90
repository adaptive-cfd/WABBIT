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
    if (Insect%smoothing_thickness=="local") then
        Insect%smooth = 1.5d0*maxval(ddx)
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

        call draw_cylinder_new( x1, x2, R, xx0, ddx, mask, mask_color, us, Insect, int(1,kind=2))
    end do

end subroutine draw_fractal_tree


subroutine fractal_tree_init(Insect)
    implicit none
    type(diptera), intent(inout) :: Insect
    integer :: nlines
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
    else
        call abort(20200504, "fractal_tree_init seems to be called twice.")
    endif

end subroutine
