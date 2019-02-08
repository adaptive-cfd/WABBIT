subroutine draw_active_grid_winglets(time, Insect, xx0, ddx, mask, mask_color, us)
    implicit none

    real(kind=rk), intent(in) :: time
    type(diptera), intent(inout) :: Insect
    real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
    real(kind=rk),intent(inout)  :: mask(0:,0:,0:)
    real(kind=rk),intent(inout)  :: us(0:,0:,0:,1:)
    integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)

    real(kind=rk), parameter :: smoothing = 1.5d0
    real(kind=rk), parameter :: offset = 4.0d0 / 20.0d0
    real(kind=rk), parameter :: omega(1:4) = (/1.0d0, -1.0d0, +1.0d0, -1.0d0/)
    real(kind=rk), parameter :: alpha0(1:4)=(/0.0d0, 0.0d0, pi/2.0d0, pi/2.0d0/)
    integer(kind=2) :: color_val

    Insect%smooth = 1.0d0*maxval(ddx)
    Insect%safety = 3.5d0*Insect%smooth

    mask = 0.00d0
    mask_color = 0
    us = 0.0d0

    color_val = 1
    call draw_single_winglet( time, (/x0, 1.5d0, 0.0d0/), omega(1), alpha0(1), 'z', color_val, xx0, &
     ddx, mask, mask_color, us, smoothing*maxval(ddx))

    color_val = 2
    call draw_single_winglet( time, (/x0, 0.5d0, 0.0d0/), omega(2), alpha0(2), 'z', color_val, xx0, &
     ddx, mask, mask_color, us, smoothing*maxval(ddx))

    color_val = 3
    call draw_single_winglet( time, (/x0-offset, 0.0d0, 1.5d0/), omega(3), alpha0(3), 'y', color_val, xx0, &
     ddx, mask, mask_color, us, smoothing*maxval(ddx))

    color_val = 4
    call draw_single_winglet( time, (/x0-offset, 0.0d0, 0.5d0/), omega(4), alpha0(4), 'y', color_val, xx0, &
     ddx, mask, mask_color, us, smoothing*maxval(ddx))


end subroutine draw_active_grid_winglets


subroutine draw_single_winglet(time, x0, omega, alpha0, orientation, color_val, xx0, ddx, mask, mask_color, us, sm)
    implicit none
    real(kind=rk), intent(in) :: time, sm
    real(kind=rk), intent(in) :: x0(1:3), omega, alpha0
    character(len=*), intent(in) :: orientation
    integer(kind=2), intent(in) :: color_val
    real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
    real(kind=rk),intent(inout)  :: mask(0:,0:,0:)
    real(kind=rk),intent(inout)  :: us(0:,0:,0:,1:)
    integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)

    real(kind=rk) :: x_glob(1:3), x_wing(1:3), H, h1, h2, h3, M=1.0d0, tmp, tmp2, MM(1:3,1:3), M2(1:3,1:3)
    real(kind=rk) :: alpha, sigma

    real(kind=rk), parameter :: scaling_winglet = 18.5/20.0
    real(kind=rk), parameter :: t_winglet = (3.0/20.0) / 2.0
    real(kind=rk), parameter :: r_axis = (3.0/20.0 ) / 2.0

    integer:: ix,iy,iz

    where (mask_color==color_val)
      mask = 0.d0
      mask_color = 0
    end where

    alpha = alpha0 + 2.0*pi*time*omega

    if (orientation == 'z') then

        call Rz( MM, alpha)
        do iz = g, size(mask,3)-1-g
            do iy = g, size(mask,2)-1-g
                do ix = g, size(mask,1)-1-g
                    ! x_glob is in the global coordinate system
                    ! note origin is shifted to x0
                    x_glob = (/ xx0(1)+dble(ix)*ddx(1), xx0(2)+dble(iy)*ddx(2), xx0(3)+dble(iz)*ddx(3) /) - x0
                    x_glob = periodize_coordinate(x_glob, (/xl,yl,zl/))
                    x_glob = x_glob / M
                    x_wing = matmul(MM, x_glob)

                    ! rotation axis: z-direction
                    if (x_wing(3)<=1.0_pr) then
                        H = -M/2.0_pr + x_wing(3)/M
                    else
                        H =  M/2.0_pr - (x_wing(3)/M -1.0_pr)
                    endif

                    H = abs(H) - (M-M*scaling_winglet)/2.0

                    tmp = max( steps( abs(x_wing(2)), H, sm) * steps( abs(x_wing(1)), t_winglet, sm ), &
                    steps( sqrt(x_wing(2)**2+x_wing(1)**2), r_axis, sm ) )

                    if (tmp >= mask(ix,iy,iz) .and. (tmp>1.0e-10) .and. &
                    (mask_color(ix,iy,iz)==0 .or. mask_color(ix,iy,iz)==color_val)) then
                        mask(ix,iy,iz) = tmp
                        mask_color(ix,iy,iz) = color_val

                        us(ix,iy,iz,1:3) = matmul( transpose(MM), (/-x_wing(2)*omega, 0.0_pr, 0.0_pr/) )
                    endif
                enddo
            enddo
        enddo

    elseif (orientation == 'y') then

        call Ry( MM, alpha)

        do iz = g, size(mask,3)-1-g
            do iy = g, size(mask,2)-1-g
                do ix = g, size(mask,1)-1-g
                    ! x_glob is in the global coordinate system
                    ! note origin is shifted to x0
                    x_glob = (/ xx0(1)+dble(ix)*ddx(1), xx0(2)+dble(iy)*ddx(2), xx0(3)+dble(iz)*ddx(3) /) - x0
                    x_glob = periodize_coordinate(x_glob, (/xl,yl,zl/))
                    x_glob = x_glob / M
                    x_wing = matmul(MM, x_glob)

                    ! rotation axis: z-direction
                    if (x_wing(2)<=1.0_pr) then
                        H =  -M/2.0_pr + x_wing(2)/M
                    else
                        H =  M/2.0_pr - (x_wing(2)/M -1.0_pr)
                    endif

                    H = abs(H) - (M-M*scaling_winglet)/2.0


                    tmp = max( steps( abs(x_wing(3)), H, sm) * steps( abs(x_wing(1)), t_winglet, sm), &
                    steps( sqrt(x_wing(3)**2+x_wing(1)**2), r_axis, sm) )

                    if (tmp >= mask(ix,iy,iz) .and. (tmp>1.0e-10) .and. &
                    (mask_color(ix,iy,iz)==0 .or. mask_color(ix,iy,iz)==color_val) ) then
                        mask(ix,iy,iz) = tmp
                        mask_color(ix,iy,iz) = color_val

                        us(ix,iy,iz,1:3) = matmul( transpose(MM), (/ +x_wing(3)*omega, 0.0_pr, 0.0_pr/) )
                    endif
                enddo
            enddo
        enddo
    endif
end subroutine draw_single_winglet
