subroutine draw_fractal_tree(Insect, xx0, ddx, mask, mask_color, us)
    implicit none

    type(diptera),intent(inout) :: Insect
    real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
    real(kind=rk),intent(inout) :: mask(0:,0:,0:)
    real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
    integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)

    integer :: io_error, nlines, mpicode,i
    real (kind=rk) :: R, safety, x1(1:3),x2(1:3), t1
    character(len=2048) :: dummy
    character(len=strlen) :: file
    real(kind=rk),allocatable,dimension(:,:) :: treedata

    ! reset everything
    mask = 0.d0
    mask_color = 0
    us = 0.d0

    ! thickness of smoothing layer and safety distance
    Insect%smooth = 1.0_rk*maxval(ddx)
    Insect%safety = 3.5_rk*Insect%smooth

    file = 'tree_data.in'
    call check_file_exists( file )

    !*****************************************************************************
    ! phase one: read the number of lines (which is the number of rigid cylinders in the tree)
    !*****************************************************************************
    call count_lines_in_ascii_file_mpi(file, nlines, 0)
    if (root) then
        write(*,'("Building a fractal tree with ",i5," rigid cylinders")') nlines
    endif

    !*****************************************************************************
    ! phase two: read all cylinders from the file into the array and bcast them
    !*****************************************************************************
    allocate( treedata(1:nlines, 1:7) )

    if (root) then
        open(unit=14,file=file,action='read',status='old')
        do i=1,nlines
            read (14,'(A)',iostat=io_error) dummy
            if (io_error==0) then
                read (dummy,*) treedata(i,:)
                write(*,'("read cylinder ",7(es12.4,1x))') treedata(i,:)
            endif
        enddo
        close (14)
    endif

    call MPI_BCAST(treedata,nlines*7,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpicode)

    !*****************************************************************************
    ! phase 3: all ranks draw the individual cylinders...
    ! please note we assume the root point to be 0 0 0 in the file
    !*****************************************************************************
    t1 = MPI_Wtime()
    do i=1, nlines
        x1 = treedata(i,1:3) + (/x0,y0,z0/)
        x2 = treedata(i,4:6) + (/x0,y0,z0/)
        ! the file containes the radius of the cylinder
        R = treedata(i,7)
        call draw_cylinder_new( x1, x2, R, xx0, ddx, mask, mask_color, us, Insect, int(1,kind=2))
    end do

    if (root) then
        write(*,'(80("-"))')
        write(*,'("done creating fractal tree mask of ",i4," branches")') nlines
        write(*,'("wtime on master process is ",es12.4," secs")') MPI_wtime()-t1
        write(*,'(80("-"))')
    end if

    deallocate(treedata)
end subroutine draw_fractal_tree
