subroutine draw_fractal_tree(Insect, xx0, ddx, mask, mask_color, us)
    implicit none

    type(diptera),intent(inout) :: Insect
    real(kind=rk),intent(in)    :: xx0(1:3), ddx(1:3)
    real(kind=rk),intent(inout) :: mask(0:,0:,0:)
    real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
    integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)

    real(kind=rk) :: x1(1:3), x2(1:3), R
    integer ::  i


    ! reset everything
    mask = 0.d0
    mask_color = 0
    us = 0.d0

    !*****************************************************************************
    ! phase 3: all ranks draw the individual cylinders...
    ! please note we assume the root point to be 0 0 0 in the file
    !*****************************************************************************
     do i = 1, size(treedata,1)
        x1 = treedata(i,1:3) + (/x0,y0,z0/)
        x2 = treedata(i,4:6) + (/x0,y0,z0/)
        ! the file containes the radius of the cylinder
        R = treedata(i,7)
        call draw_cylinder_new( x1, x2, R, xx0, ddx, mask, mask_color, us, Insect, int(1,kind=2))
    end do

end subroutine draw_fractal_tree


subroutine fractal_tree_init()
    implicit none
    character(len=80) :: file
    integer :: nlines
    ! note: smoothing is set in insect_init (which must be called also for fractal trees)

    ! initialization
    if (.not. allocated(treedata)) then
        file = 'tree_data.in'
        call check_file_exists( file )

        !*****************************************************************************
        ! phase one: read the number of lines (which is the number of rigid cylinders in the tree)
        !*****************************************************************************
        call count_lines_in_ascii_file_mpi(file, nlines, 0)
        if (root) write(*,'("Building a fractal tree with ",i5," rigid cylinders")') nlines

        !*****************************************************************************
        ! phase two: read all cylinders from the file into the array and bcast them
        !*****************************************************************************
        allocate( treedata(1:nlines, 1:7) )
        call read_array_from_ascii_file_mpi(file, treedata, n_header=0)

        if (root) write(*,'("Done reading ",i5," rigid cylinders")') nlines
    endif

end subroutine
