subroutine mask_from_pointcloud(points, normals, xx0, ddx, mask, safety, h_smooth, d0, mask_color, color)
    implicit none
    ! point cloud data (points and normals)
    real(kind=rk),dimension(1:,1:),intent(in):: points, normals
    ! xx0 and spacing of grid
    real(kind=rk),dimension(1:3),intent(in):: xx0, ddx
    real(kind=rk),intent(inout):: mask(0:,0:,0:)
    integer, intent(in) :: safety
    real(kind=rk),intent(in):: d0, h_smooth
    integer(kind=2),intent(in), optional :: color
    integer(kind=2),intent(inout), optional :: mask_color(0:,0:,0:)

    integer, dimension(1:3) :: lbounds, ubounds
    integer :: ix,iy,iz,i,xmin,xmax,ymin,ymax,zmin,zmax,N, mpisize, mpicode
    real(kind=rk) :: x,y,z,tmp,t0,sign

    lbounds = 0
    ubounds = (/size(mask,1), size(mask,2), size(mask,3)/) - 1

    t0 = MPI_wtime()
    ! number of points in cloud
    N = size(points,1)

    ! attention we currently reset EVERYTHING!!!!!
    ! NOTE HACK: in the future, we should not do that but do something
    ! clever instead. As long as were using the interpolation approach, it does not
    ! hurt, though
    mask = 9.0d7


    if (maxval(points(:,1))>xx0(1)+dble(ubounds(1))*ddx(1) .or. minval(points(:,1))<xx0(1)) then
        write(*,*) "WARNING: some points may be outside of the domain! (x) oh-oh"
    endif
    if (maxval(points(:,2))>xx0(2)+dble(ubounds(2))*ddx(2) .or. minval(points(:,2))<xx0(2)) then
        write(*,*) "WARNING: some points may be outside of the domain! (y) oh-oh"
    endif
    if (maxval(points(:,3))>xx0(3)+dble(ubounds(3))*ddx(3) .or. minval(points(:,3))<xx0(3)) then
        write(*,*) "WARNING: some points may be outside of the domain! (z) oh-oh"
    endif


    do i = 1, N
        ! bounding box of the vicinity of the Lagrangian maker point of the cloud.
        xmin = nint((points(i,1)-xx0(1))/ddx(1))-safety
        xmax = nint((points(i,1)-xx0(1))/ddx(1))+safety

        ymin = nint((points(i,2)-xx0(2))/ddx(2))-safety
        ymax = nint((points(i,2)-xx0(2))/ddx(2))+safety

        zmin = nint((points(i,3)-xx0(3))/ddx(3))-safety
        zmax = nint((points(i,3)-xx0(3))/ddx(3))+safety

        do iz = max(zmin,lbounds(3)), min(zmax,ubounds(3))
            z = xx0(3) + ddx(3)*dble(iz)
            do iy = max(ymin,lbounds(2)), min(ymax,ubounds(2))
                y = xx0(2) + ddx(2)*dble(iy)
                do ix = max(xmin,lbounds(1)), min(xmax,ubounds(1))
                    x = xx0(1) + ddx(1)*dble(ix)

                    !-------------------------------------------------------------
                    ! the distance to the current point:
                    ! note this is the square of the distance (cheaper to take sqrt later!)
                    tmp = (x-points(i,1))*(x-points(i,1)) + (y-points(i,2))*(y-points(i,2)) &
                    + (z-points(i,3))*(z-points(i,3))

                    ! if closer (in abs value!) then use this now
                    if ( dabs(tmp) < dabs(mask(ix,iy,iz)) ) then
                        ! take care of sign: exterior / interior
                        sign = (x-points(i,1))*normals(i,1) + (y-points(i,2))*normals(i,2) &
                        + (z-points(i,3))*normals(i,3)
                        ! is interior?
                        if (sign < 0.0d0) then
                            tmp = -tmp
                        endif
                        ! use this value now
                        mask(ix,iy,iz)  = tmp
                    endif
                    !-------------------------------------------------------------
                enddo
            enddo
        enddo

    enddo

    !-----------------------------------------------------------------------------
    ! convert signed distance function to mask function chi
    !-----------------------------------------------------------------------------
    do iz = lbounds(3), ubounds(3)
        do iy = lbounds(2), ubounds(2)
            do ix = lbounds(1), ubounds(1)
                tmp = mask(ix,iy,iz)
                ! exclude area that has not been altered:
                if (tmp<1.0d7) then
                    if (tmp >= 0.0d0) then
                        tmp = dsqrt(tmp)
                    else
                        tmp = -dsqrt(-tmp)
                    endif
                    mask(ix,iy,iz) = steps( tmp-d0, 0.d0, h_smooth )
                else
                    mask(ix,iy,iz) = 0.0d0
                endif
            enddo
        enddo
    enddo

    ! assign color, if that is used
    if (present(mask_color) .and. present(color)) then
        ! reset color. NOTE HACK: in the future, we should not do that but do something
        ! clever instead. As long as were using the interpolation approach, it does not
        ! hurt, though
        mask_color = 0

        do iz = lbounds(3), ubounds(3)
            do iy = lbounds(2), ubounds(2)
                do ix = lbounds(1), ubounds(1)
                    if (mask(ix,iy,iz)>0.0) mask_color(ix,iy,iz)=color
                enddo
            enddo
        enddo
    endif

    t0 = MPI_wtime() -t0
    t0 = mpisum(t0)

    call MPI_COMM_SIZE (MPI_COMM_WORLD,mpisize,mpicode)
    if (root) write(*,*) "elapsed time in pc2chi ", t0/dble(mpisize)

end subroutine mask_from_pointcloud
