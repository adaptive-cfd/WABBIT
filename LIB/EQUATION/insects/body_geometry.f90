
! the old routine (<02/2019) created just the body mask but not its solid velocity field
! now the routine does create the full body mask, including the velocity field. It sets us
! only inside the body (hence the wings set u_wing + u_body)
subroutine draw_insect_body( time, xx0, ddx, mask, mask_color, us, Insect, delete)
    implicit none

    real(kind=rk), intent(in)    :: time
    type(diptera), intent(inout) :: Insect
    real(kind=rk), intent(in)    :: xx0(1:3), ddx(1:3)
    real(kind=rk), intent(inout) :: mask(0:,0:,0:)
    real(kind=rk), intent(inout) :: us(0:,0:,0:,1:)
    integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)
    logical, intent(in)           :: delete

    integer                       :: ix, iy, iz
    real(kind=rk), dimension(1:3) :: x_glob, x_body, v_tmp

    ! 28/01/2019: Thomas. Discovered that this was done block based, i.e. the smoothing layer
    ! had different thickness, if some blocks happened to be at different levels (and still carry
    ! a part of the smoothing layer.) I don't know if that made sense, because the layer shrinks/expands then
    ! and because it might be discontinous. Both options are included now, default is "as before"
    ! Insect%smoothing_thickness=="local"  : smoothing_layer = c_sm * 2**-J * L/(BS-1)
    ! Insect%smoothing_thickness=="global" : smoothing_layer = c_sm * 2**-Jmax * L/(BS-1)
    ! NOTE: for FLUSI, this has no impact! Here, the grid is constant and equidistant.
    if (Insect%smoothing_thickness=="local" .or. .not. grid_time_dependent) then
        Insect%smooth = 1.0d0*maxval(ddx)
        Insect%safety = 3.5d0*Insect%smooth
    endif

    if (size(mask) /= size(mask_color) .or. size(us,4) /= 3) then
        write(*,*) "mask:", shape(mask), "mask_color:", shape(mask_color), "us:", shape(us)
        call abort (08021901,"Insects: arrays have wrong size..")
    endif

    if ((dabs(Insect%time-time)>1.0d-10).and.root) then
        write(*,'("error! time=",es15.8," but Insect%time=",es15.8)') time, Insect%time
        write(*,'("Did you call Update_Insect before draw_insect_body?")')
    endif


    if (delete) then
        if (grid_time_dependent) then
            ! The grid is time-dependent. In this case, the separation between
            ! time-dependent (wings, moving body) and time-independent (fixed body)
            ! is done elsewhere, so deleting means delete entire block
            mask = 0.00_rk
            us(:,:,:,1) = 0.00_rk
            us(:,:,:,2) = 0.00_rk
            us(:,:,:,3) = 0.00_rk
            mask_color = 0
        else
            ! for the fixed-grid codes, delete only the body.
            where (mask_color==Insect%color_body)
                mask = 0.00_rk
                us(:,:,:,1) = 0.00_rk
                us(:,:,:,2) = 0.00_rk
                us(:,:,:,3) = 0.00_rk
                mask_color = 0
            end where
        endif
    endif

    !---------------------------------------------------------------------------
    ! stage I:
    !---------------------------------------------------------------------------
    ! create the body mask, not the solid velocity field.
    select case (Insect%BodyType)
    case ("nobody")
        return
    case ("suzuki_thin_rod")
        call draw_suzuki_thin_rod( xx0, ddx, mask, mask_color, us, Insect)

    case ("superSTL")
        call draw_body_superSTL( xx0, ddx, mask, mask_color, us, Insect)

    case ("jerry","Jerry")
        call draw_body_jerry( xx0, ddx, mask, mask_color, us, Insect)

    case ("hawkmoth","Hawkmoth")
        call draw_body_hawkmoth( xx0, ddx, mask, mask_color, us, Insect)

    case ("platicle")
        call draw_body_platicle( xx0, ddx, mask, mask_color, us, Insect)

    case ("coin")
        call draw_body_coin( xx0, ddx, mask, mask_color, us, Insect)

    case ("sphere","SPHERE","Sphere")
        call draw_body_sphere( xx0, ddx, mask, mask_color, us, Insect)

    case ("drosophila_maeda","drosophila_slim")
        call draw_body_drosophila_maeda( xx0, ddx, mask, mask_color, us, Insect)

    case ("bumblebee")
        call draw_body_bumblebee( xx0, ddx, mask, mask_color, us, Insect)

    case ("paratuposa_simple")
        call draw_body_paratuposa_simple( xx0, ddx, mask, mask_color, us, Insect)

    case ("mosquito_iams")
        call draw_body_mosquito_iams( xx0, ddx, mask, mask_color, us, Insect)

    case ("cone")
        call draw_body_cone( xx0, ddx, mask, mask_color, us, Insect)

    case default
        call abort(10623, "Insect::draw_insect_body::Insect%BodyType unknown..."//trim(adjustl(Insect%BodyType)))

    end select

    !---------------------------------------------------------------------------
    ! stage II:
    !---------------------------------------------------------------------------
    ! if the body does not move, we can skip the creation of us
    if (Insect%body_moves == "no") return

    ! add the solid velocity field to the body mask (i.e. create us)
    do iz = g, size(mask,3)-1-g
        x_glob(3) = xx0(3) + dble(iz)*ddx(3) - Insect%xc_body_g(3)
        do iy = g, size(mask,2)-1-g
            x_glob(2) = xx0(2) + dble(iy)*ddx(2) - Insect%xc_body_g(2)
            do ix = g, size(mask,1)-1-g
                x_glob(1) = xx0(1) + dble(ix)*ddx(1) - Insect%xc_body_g(1)

                ! skip all parts that do not belong to the body (ie they have a different color)
                if ( mask_color(ix,iy,iz) == Insect%color_body .and. mask(ix,iy,iz) > 0.d0 ) then

                    if (periodic_insect) x_glob = periodize_coordinate(x_glob, (/xl,yl,zl/))
                    x_body = matmul(Insect%M_body, x_glob)

                    ! add solid body rotation to the translational velocity field. Note
                    ! that rot_body_b and x_body are in the body reference frame
                    v_tmp(1) = Insect%rot_body_b(2)*x_body(3)-Insect%rot_body_b(3)*x_body(2)
                    v_tmp(2) = Insect%rot_body_b(3)*x_body(1)-Insect%rot_body_b(1)*x_body(3)
                    v_tmp(3) = Insect%rot_body_b(1)*x_body(2)-Insect%rot_body_b(2)*x_body(1)

                    ! the body motion is transformed to the global system, translation is added
                    us(ix,iy,iz,1:3) = matmul( Insect%M_body_inv, v_tmp ) + Insect%vc_body_g
                endif
            enddo
        enddo
    enddo

end subroutine


!-------------------------------------------------------------------------------

! Bumblebee body, BB1 in Dudley & Ellington JEB 1990
subroutine draw_body_bumblebee( xx0, ddx, mask, mask_color, us, Insect)
    implicit none

    type(diptera),intent(inout) :: Insect
    real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
    real(kind=rk),intent(inout) :: mask(0:,0:,0:)
    real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
    integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)

    integer :: ix,iy,iz,j
    real(kind=rk) :: x,y,z,s,s1,a_body,R,R0,R_tmp,x1
    real(kind=rk) :: x_glob(1:3),x_body(1:3),x_head(1:3),xa(1:3),xb(1:3)
    real(kind=rk) :: rbc,thbc1,thbc2,x0bc,z0bc,xcs,zcs
    real(kind=rk) :: xx_head,zz_head,dx_head,dz_head,a_head
    real(kind=rk) :: xl1(5),yl1(5),zl1(5),rl1(4),xl2(5),yl2(5),zl2(5),rl2(4),&
    xl3(5),yl3(5),zl3(5),rl3(4),xf(2),yf(2),zf(2),rf,xan(2),yan(2),zan(2),ran,&
    xmin_bbox,xmax_bbox,ymin_bbox,ymax_bbox,zmin_bbox,zmax_bbox
    real(kind=rk)   :: M_body(1:3,1:3)
    integer(kind=2) :: color_body

    color_body = Insect%color_body
    M_body     = Insect%M_body

    !-----------------------------------------------------------------------------
    ! Body
    !-----------------------------------------------------------------------------
    do iz = g, size(mask,3)-1-g
        x_glob(3) = xx0(3) + dble(iz)*ddx(3) - Insect%xc_body_g(3)
        do iy = g, size(mask,2)-1-g
            x_glob(2) = xx0(2) + dble(iy)*ddx(2) - Insect%xc_body_g(2)
            do ix = g, size(mask,1)-1-g
                x_glob(1) = xx0(1) + dble(ix)*ddx(1) - Insect%xc_body_g(1)
                if (periodic_insect) x_glob = periodize_coordinate(x_glob, (/xl,yl,zl/))
                ! x_body is in the body coordinate system
                x_body = matmul(M_body,x_glob)

                ! ------------------------------------
                ! approximation to mesh
                ! ------------------------------------
                x = x_body(1)
                y = x_body(2)
                z = x_body(3)

                ! symmetry plane is xz
                ! +x direction is forward
                ! body centerline is an arc with center at x0bc,y0bc
                ! radius rbc and angles th counted from negative z
                rbc = 1.3d0
                thbc1 = 112.0d0 *pi/180.0d0
                thbc2 = 53.0d0 *pi/180.0d0
                x0bc = 0.0782301255230126d0
                z0bc = -1.26512552301255d0

                ! chordwise dimensionless coordinate, from head to abdomen
                s = (datan2(z-z0bc,-(x-x0bc))-thbc1)/(thbc2-thbc1)
                ! body center coordinates at s
                xcs = x0bc + (x-x0bc)*rbc/dsqrt((x-x0bc)**2+(z-z0bc)**2)
                zcs = z0bc + (z-z0bc)*rbc/dsqrt((x-x0bc)**2+(z-z0bc)**2)

                ! check if inside body bounds (in s-direction)
                if ( (s>=-Insect%safety) .and. (s<=1.075d0+Insect%safety) ) then
                    R0 = 0.0d0
                    ! round section by default
                    a_body = 1.0d0
                    ! distortion of s
                    s1 = 1.0d0 - ( s + 0.08d0*dtanh(30.0d0*s) ) / (1.0d0+0.08d0*dtanh(30.0d0))
                    s1 = ( s1 + 0.04d0*dtanh(60.0d0*s1) ) / (1.0d0+0.04d0*dtanh(60.0d0))
                    s1 = dsin(1.2d0*s1)/dsin(1.2d0)
                    s1 = sign(abs(s1)**1.25,s1)

                    x1 = 1.075d0 * s1
                    ! compute radius as a function of x1 (counting from the tail on)
                    ! same shape as 'drosophila'
                    if (x1 < 0.6333d0) then
                        ! we're in the ABDOMEN
                        R0 = max( -1.2990d0*x1**2 + 0.9490d0*x1 + 0.0267d0, 0.d0)
                        ! flatten abdomen
                        a_body = 1.0d0-0.07d0*(x1-0.6333d0)*x1/0.0488d0
                    elseif ((x1 >= 0.6333d0) .and. (x1 <=1.075d0 )) then
                        ! we're in the THORAX
                        R0 = max( -2.1667d0*x1**2 + 3.4661d0*x1 - 1.2194d0, 0.d0)
                    endif
                    ! distortion of R0
                    R0 = 1.2d0 * (1.0d0+0.6d0*(1.0d0-s)**2) * R0
                    ! distance to the body center at s
                    R = dsqrt( (x-xcs)**2 + y**2 + (a_body*(z-zcs))**2 )

                    ! smoothing
                    if (( R < R0 + Insect%safety ).and.(R0>0.d0)) then
                        R_tmp = steps(R, R0, Insect%smooth)
                        mask(ix,iy,iz)= max( R_tmp , mask(ix,iy,iz) )
                        mask_color(ix,iy,iz) = color_body
                    endif

                endif
            enddo
        enddo
    enddo

    !-----------------------------------------------------------------------------
    ! Head
    !-----------------------------------------------------------------------------
    a_head = 1.04d0

    ! ellipsoid head, assumes xc_head=0 in .ini file
    xx_head = 0.58125d0
    zz_head = -0.1d0
    dx_head = 0.5d0 * 0.2035d0
    dz_head = 0.5d0 * 0.297d0

    do iz = g, size(mask,3)-1-g
        x_glob(3) = xx0(3) + dble(iz)*ddx(3) - Insect%xc_body_g(3)
        do iy = g, size(mask,2)-1-g
            x_glob(2) = xx0(2) + dble(iy)*ddx(2) - Insect%xc_body_g(2)
            do ix = g, size(mask,1)-1-g
                x_glob(1) = xx0(1) + dble(ix)*ddx(1) - Insect%xc_body_g(1)
                if (periodic_insect) x_glob = periodize_coordinate(x_glob, (/xl,yl,zl/))
                x_body   = matmul(M_body,x_glob)
                x_head   = x_body

                ! check if inside the surrounding box (save comput. time)
                if ( dabs(x_head(2)) <= dz_head + Insect%safety ) then
                    if ( dabs(x_head(3)-zz_head) <= dz_head + Insect%safety ) then
                        ! check for length inside ellipsoid:
                        if ( dabs(x_head(1)-xx_head) < dx_head + Insect%safety ) then

                            R  = dsqrt ( (a_head*x_head(2))**2 + (x_head(3)-zz_head)**2 )
                            ! this gives the R(x) shape
                            if ( ((x_head(1)-xx_head)/dx_head)**2 <= 1.d0) then
                                R0 = dz_head*dsqrt(1.d0- ((x_head(1)-xx_head)/dx_head)**2 )
                                if ( R < R0 + Insect%safety ) then
                                    mask(ix,iy,iz)= max(steps(R,R0, Insect%smooth),mask(ix,iy,iz))
                                    mask_color(ix,iy,iz) = color_body
                                endif
                            endif
                        endif
                    endif
                endif


            enddo
        enddo
    enddo

    !-----------------------------------------------------------------------------
    ! Legs, antennae and proboscis
    !-----------------------------------------------------------------------------
    if ((Insect%HasDetails=="all").or.(Insect%HasDetails=="legs").or.(Insect%HasDetails=="antennae_proboscis")) then
        ! Parameters of legs, antennae and proboscis
        xl1 = (/-0.74,-0.63,-0.4,-0.1,0.1/)
        yl1 = (/0.32,0.32,0.31,0.3,0.12/)
        zl1 = (/-0.35,-0.37,-0.2,-0.1,-0.16/)
        rl1 = (/0.015,0.03,0.04,0.03/)*1.3
        xl2 = (/-0.24,-0.15,0.02,0.17,0.19/)
        yl2 = (/0.33,0.33,0.32,0.3,0.15/)
        zl2 = (/-0.29,-0.28,-0.2,-0.15,-0.19/)
        rl2 = (/0.015,0.03,0.04,0.03/)*1.3
        xl3 = (/0.28,0.35,0.45,0.4,0.35/)
        yl3 = (/0.31,0.30,0.28,0.2,0.15/)
        zl3 = (/-0.3,-0.28,-0.25,-0.18,-0.18/)
        rl3 = (/0.015,0.02,0.03,0.02/)*1.3
        xf = (/0.43,0.6/)
        yf = (/0.0,0.0/)
        zf = (/-0.28,-0.23/)
        rf = 0.017*1.3
        xan = (/0.63,0.8/)
        yan = (/0.05,0.27/)
        zan = (/-0.03,0.1/)
        ran = 0.015*1.3

        ! ----------------------------------------------------------------------------
        ! legs (composed of 3 cylinders)
        ! ----------------------------------------------------------------------------
        do j = 1,  4
          ! transform coordinates to global system. they are defined in the body system
          xa = matmul( transpose(M_body), (/xl1(j)  ,yl1(j)  ,zl1(j)/)  ) + Insect%xc_body_g
          xb = matmul( transpose(M_body), (/xl1(j+1),yl1(j+1),zl1(j+1)/)) + Insect%xc_body_g
          ! note input to draw_cylinder_new is in global coordinates
          call draw_cylinder_new( xa, xb, rl1(j), xx0, ddx, mask, mask_color, us, Insect, color_body)

          xa = matmul( transpose(M_body), (/xl2(j)  ,yl2(j)  ,zl2(j)/)  ) + Insect%xc_body_g
          xb = matmul( transpose(M_body), (/xl2(j+1),yl2(j+1),zl2(j+1)/)) + Insect%xc_body_g
          call draw_cylinder_new( xa, xb, rl2(j), xx0, ddx, mask, mask_color, us, Insect, color_body)

          xa = matmul( transpose(M_body), (/xl3(j)  ,yl3(j)  ,zl3(j)/)  ) + Insect%xc_body_g
          xb = matmul( transpose(M_body), (/xl3(j+1),yl3(j+1),zl3(j+1)/)) + Insect%xc_body_g
          call draw_cylinder_new( xa, xb, rl3(j), xx0, ddx, mask, mask_color, us, Insect, color_body)

          ! right side of body (flip the sign of y)
          xa = matmul( transpose(M_body), (/xl1(j)  ,-yl1(j)  ,zl1(j)/)  ) + Insect%xc_body_g
          xb = matmul( transpose(M_body), (/xl1(j+1),-yl1(j+1),zl1(j+1)/)) + Insect%xc_body_g
          call draw_cylinder_new( xa, xb, rl1(j), xx0, ddx, mask, mask_color, us, Insect, color_body)

          xa = matmul( transpose(M_body), (/xl2(j)  ,-yl2(j)  ,zl2(j)/)  ) + Insect%xc_body_g
          xb = matmul( transpose(M_body), (/xl2(j+1),-yl2(j+1),zl2(j+1)/)) + Insect%xc_body_g
          call draw_cylinder_new( xa, xb, rl2(j), xx0, ddx, mask, mask_color, us, Insect, color_body)

          xa = matmul( transpose(M_body), (/xl3(j)  ,-yl3(j)  ,zl3(j)/)  ) + Insect%xc_body_g
          xb = matmul( transpose(M_body), (/xl3(j+1),-yl3(j+1),zl3(j+1)/)) + Insect%xc_body_g
          call draw_cylinder_new( xa, xb, rl3(j), xx0, ddx, mask, mask_color, us, Insect, color_body)
        enddo

        ! antenna (left)
        xa = matmul( transpose(M_body), (/xan(1),yan(1),zan(1)/) ) + Insect%xc_body_g
        xb = matmul( transpose(M_body), (/xan(2),yan(2),zan(2)/) ) + Insect%xc_body_g
        call draw_cylinder_new( xa, xb, ran, xx0, ddx, mask, mask_color, us, Insect, color_body)

        ! antenna (right)
        xa = matmul( transpose(M_body), (/xan(1),-yan(1),zan(1)/) ) + Insect%xc_body_g
        xb = matmul( transpose(M_body), (/xan(2),-yan(2),zan(2)/) ) + Insect%xc_body_g
        call draw_cylinder_new( xa, xb, ran, xx0, ddx, mask, mask_color, us, Insect, color_body)

        ! proboscis (to drink)
        xa = matmul( transpose(M_body), (/xf(1),yf(1),zf(1)/) ) + Insect%xc_body_g
        xb = matmul( transpose(M_body), (/xf(2),yf(2),zf(2)/) ) + Insect%xc_body_g
        call draw_cylinder_new( xa, xb, rf, xx0, ddx, mask, mask_color, us, Insect, color_body)
    endif

end subroutine draw_body_bumblebee

!------------------------------------------------------------------------------
! A very small bug Paratuposa, highly simplified body shape
subroutine draw_body_paratuposa_simple( xx0, ddx, mask, mask_color, us, Insect)
    implicit none

    type(diptera),intent(inout) :: Insect
    real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
    real(kind=rk),intent(inout) :: mask(0:,0:,0:)
    real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
    integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)

    integer :: ix,iy,iz,j
    real(kind=rk) :: x,y,z,s,a_body,bodylen,R,R0,R_tmp
    real(kind=rk) :: x_glob(1:3),x_body(1:3),x_head(1:3),xa(1:3),xb(1:3)
    real(kind=rk) :: rbc,thbc1,thbc2,x0bc,z0bc,xcs,zcs
    real(kind=rk) :: xx_head,zz_head,dx_head,dz_head,a_head
    real(kind=rk) :: xl1(5),yl1(5),zl1(5),rl1(4),xl2(5),yl2(5),zl2(5),rl2(4),&
    xl3(5),yl3(5),zl3(5),rl3(4),xf(2),yf(2),zf(2),rf,xan(2),yan(2),zan(2),ran,&
    xmin_bbox,xmax_bbox,ymin_bbox,ymax_bbox,zmin_bbox,zmax_bbox
    real(kind=rk)   :: M_body(1:3,1:3)
    integer(kind=2) :: color_body

    ! Body length relative to the wing length
    bodylen = 0.84d0

    color_body = Insect%color_body
    M_body     = Insect%M_body

    do iz = g, size(mask,3)-1-g
        x_glob(3) = xx0(3) + dble(iz)*ddx(3) - Insect%xc_body_g(3)
        do iy = g, size(mask,2)-1-g
            x_glob(2) = xx0(2) + dble(iy)*ddx(2) - Insect%xc_body_g(2)
            do ix = g, size(mask,1)-1-g
                x_glob(1) = xx0(1) + dble(ix)*ddx(1) - Insect%xc_body_g(1)
                if (periodic_insect) x_glob = periodize_coordinate(x_glob, (/xl,yl,zl/))
                ! x_body is in the body coordinate system
                x_body = matmul(M_body,x_glob)

                ! ------------------------------------
                ! approximation to mesh
                ! ------------------------------------
                x = x_body(1)/bodylen
                y = x_body(2)/bodylen
                z = x_body(3)/bodylen

                ! symmetry plane is xz
                ! +x direction is forward
                ! body centerline is an arc with center at x0bc,y0bc
                ! radius rbc and angles th measured from negative z
                rbc = 1.27d0
                thbc1 = 113.0d0 *pi/180.0d0
                thbc2 = 67.0d0 *pi/180.0d0
                x0bc = 0.0d0
                z0bc = -1.17d0

                ! chordwise dimensionless coordinate, from head to abdomen
                s = (datan2(z-z0bc,-(x-x0bc))-thbc1)/(thbc2-thbc1)
                ! body center coordinates at s
                xcs = x0bc + (x-x0bc)*rbc/dsqrt((x-x0bc)**2+(z-z0bc)**2)
                zcs = z0bc + (z-z0bc)*rbc/dsqrt((x-x0bc)**2+(z-z0bc)**2)

                ! check if inside body bounds (in s-direction)
!                if ( (s>=-Insect%safety) .and. (s<=1.0d0+Insect%safety) ) then
                if ( (s>=0.0d0) .and. (s<=1.0d0) ) then

                    ! compute radius as a function of s (counting from the tail on)
                    if (z >= zcs) then
                        R0 = 0.22d0/0.5d0*dsqrt(0.5d0**2-(s-0.5d0)**2)
                        a_body = 0.18d0/0.22d0
                    else
                        R0 = 0.07d0/0.5d0*dsqrt(0.5d0**2-(s-0.5d0)**2)
                        a_body = 0.18d0/0.07d0
                    endif

                    ! distance to the body center at s
                    R = dsqrt( (x-xcs)**2 + (y/a_body)**2 + (z-zcs)**2 )

                    ! smoothing
                    if (( R < R0 + Insect%safety ).and.(R0>0.d0)) then
                        R_tmp = steps(R, R0, Insect%smooth)
                        mask(ix,iy,iz)= max( R_tmp , mask(ix,iy,iz) )
                        mask_color(ix,iy,iz) = color_body
                    endif

                endif
            enddo
        enddo
    enddo

end subroutine draw_body_paratuposa_simple

!------------------------------------------------------------------------------
! Body adapted from Maeda & Liu. It assumes Insect%x_head=0.0
subroutine draw_body_drosophila_maeda( xx0, ddx, mask, mask_color, us, Insect)
    implicit none

    type(diptera),intent(inout) :: Insect
    real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
    real(kind=rk),intent(inout) :: mask(0:,0:,0:)
    real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
    integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)

    integer :: ix,iy,iz
    real(kind=rk) :: x,y,z,s,s1, a_body, R,R0,R_tmp,x1, a_body0
    real(kind=rk) :: x_glob(1:3),x_body(1:3),x_head(1:3)
    real(kind=rk) :: rbc,thbc1,thbc2,x0bc,z0bc,xcs,zcs
    real(kind=rk) :: xx_head,zz_head,dx_head,dz_head,a_head
    real(kind=rk)   :: M_body(1:3,1:3)
    integer(kind=2) :: color_body

    color_body = Insect%color_body
    M_body     = Insect%M_body


    !-----------------------------------------------------------------------------
    ! Body
    !-----------------------------------------------------------------------------
    if (Insect%BodyType == 'drosophila_slim') then
        a_body0 = 1.09d0
    else
        a_body0 = 1.0d0
    endif

    do iz = g, size(mask,3)-1-g
        x_glob(3) = xx0(3) + dble(iz)*ddx(3) - Insect%xc_body_g(3)
        do iy = g, size(mask,2)-1-g
            x_glob(2) = xx0(2) + dble(iy)*ddx(2) - Insect%xc_body_g(2)
            do ix = g, size(mask,1)-1-g
                x_glob(1) = xx0(1) + dble(ix)*ddx(1) - Insect%xc_body_g(1)
                if (periodic_insect) x_glob = periodize_coordinate(x_glob, (/xl,yl,zl/))
                ! x_body is in the body coordinate system
                x_body = matmul(M_body,x_glob)

                ! ------------------------------------
                ! approximation to mesh from Maeda
                ! similar to Aono et al.
                ! ------------------------------------
                x = x_body(1)
                y = x_body(2)
                z = x_body(3)

                ! symmetry plane is xz
                ! +x direction is forward
                ! body centerline is an arc with center at x0bc,y0bc
                ! radius rbc and angles th counted from negative z
                rbc = 0.9464435146443515d0
                thbc1 = 112.0d0 *pi/180.0d0
                thbc2 = 53.0d0 *pi/180.0d0
                x0bc = -0.24476987447698745d0
                z0bc = -0.9301255230125524d0

                ! chordwise dimensionless coordinate, from head to abdomen
                s = (datan2(z-z0bc,-(x-x0bc))-thbc1)/(thbc2-thbc1)
                ! body center coordinates at s
                xcs = x0bc + (x-x0bc)*rbc/dsqrt((x-x0bc)**2+(z-z0bc)**2)
                zcs = z0bc + (z-z0bc)*rbc/dsqrt((x-x0bc)**2+(z-z0bc)**2)

                ! check if inside body bounds (in s-direction)
                if ( (s>=-Insect%safety) .and. (s<=1.075d0+Insect%safety) ) then
                    R0 = 0.0d0
                    a_body = a_body0
                    ! distortion of s
                    s1 = 1.0d0 - ( s + 0.08d0*dtanh(30.0d0*s) ) / (1.0d0+0.08d0*dtanh(30.0d0))
                    s1 = ( s1 + 0.04d0*dtanh(60.0d0*s1) ) / (1.0d0+0.04d0*dtanh(60.0d0))

                    ! s1 = ( dsin(1.2d0*s1)/dsin(1.2d0) )**1.25
                    s1 = dsin(1.2d0*s1)/dsin(1.2d0)
                    s1 = sign(abs(s1)**1.25,s1)

                    x1 = 1.075d0 * s1
                    ! compute radius as a function of x1 (counting from the tail on)
                    ! same shape as 'drosophila'
                    if (x1 < 0.6333d0) then
                        ! we're in the ABDOMEN
                        R0 = max( -1.2990d0*x1**2 + 0.9490d0*x1 + 0.0267d0, 0.d0)
                    elseif ((x1 >= 0.6333d0) .and. (x1 <=1.075d0 )) then
                        ! we're in the THORAX
                        R0 = max( -2.1667d0*x1**2 + 3.4661d0*x1 - 1.2194d0, 0.d0)
                        ! slim body
                        if (Insect%BodyType == 'drosophila_slim') &
                        a_body = 1.09d0-0.19d0*(x1-0.6333d0)*(x1-1.075d0)/0.0488d0
                    endif
                    ! distortion of R0
                    R0 = 0.8158996d0 * (1.0d0+0.6d0*(1.0d0-s)**2) * R0
                    ! distance to the body center at s
                    R = dsqrt( (x-xcs)**2 + (a_body*y)**2 + (z-zcs)**2 )

                    ! smoothing
                    if (( R < R0 + Insect%safety ).and.(R0>0.d0)) then
                        R_tmp = steps(R,R0, Insect%smooth)
                        mask(ix,iy,iz)= max( R_tmp , mask(ix,iy,iz) )
                        mask_color(ix,iy,iz) = color_body
                    endif

                endif
            enddo
        enddo
    enddo

    !-----------------------------------------------------------------------------
    ! Head
    !-----------------------------------------------------------------------------
    if (Insect%BodyType == 'drosophila_slim') then
        a_head = 1.09d0
    else
        a_head = 1.0d0
    endif

    ! ellipsoid head, assumes xc_head=0 in .ini file
    xx_head = 0.17d0
    zz_head = -0.1d0
    dx_head = 0.5d0 * 0.185d0
    dz_head = 0.5d0 * 0.27d0

    do iz = g, size(mask,3)-1-g
        x_glob(3) = xx0(3) + dble(iz)*ddx(3) - Insect%xc_body_g(3)
        do iy = g, size(mask,2)-1-g
            x_glob(2) = xx0(2) + dble(iy)*ddx(2) - Insect%xc_body_g(2)
            do ix = g, size(mask,1)-1-g
                x_glob(1) = xx0(1) + dble(ix)*ddx(1) - Insect%xc_body_g(1)
                if (periodic_insect) x_glob = periodize_coordinate(x_glob, (/xl,yl,zl/))
                ! x_body is in the body coordinate system
                x_body = matmul(M_body,x_glob)
                x_head   = x_body

                ! check if inside the surrounding box (save comput. time)
                if ( dabs(x_head(2)) <= dz_head + Insect%safety ) then
                    if ( dabs(x_head(3)-zz_head) <= dz_head + Insect%safety ) then
                        ! check for length inside ellipsoid:
                        if ( dabs(x_head(1)-xx_head) < dx_head + Insect%safety ) then

                            R  = dsqrt ( (a_head*x_head(2))**2 + (x_head(3)-zz_head)**2 )
                            ! this gives the R(x) shape
                            if ( ((x_head(1)-xx_head)/dx_head)**2 <= 1.d0) then
                                R0 = dz_head*dsqrt(1.d0- ((x_head(1)-xx_head)/dx_head)**2 )
                                if ( R < R0 + Insect%safety ) then
                                    mask(ix,iy,iz)= max(steps(R,R0, Insect%smooth),mask(ix,iy,iz))
                                    mask_color(ix,iy,iz) = color_body
                                endif
                            endif
                        endif
                    endif
                endif


            enddo
        enddo
    enddo
end subroutine draw_body_drosophila_maeda


!-------------------------------------------------------------------------------
subroutine draw_body_jerry( xx0, ddx, mask, mask_color, us, Insect)
    implicit none

    type(diptera),intent(inout) :: Insect
    real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
    real(kind=rk),intent(inout) :: mask(0:,0:,0:)
    real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
    integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)

    real(kind=rk) :: R0,R,a_body
    real(kind=rk) :: x_body(1:3), x_glob(1:3), x_head(1:3), x_eye(1:3)
    real(kind=rk)   :: M_body(1:3,1:3)
    integer(kind=2) :: color_body
    integer :: ix,iy,iz

    color_body = Insect%color_body
    M_body     = Insect%M_body

    ! the following are coordinates of specific points on the insect's body, for
    ! example the position of the head, its size etc. In older versions, these
    ! parameters were set in the *.ini file, which proved to be too much flexibility.
    ! in practice, the insect is created once, while implementing it, and then
    ! no longer changed. For Jerry, we overwrite the coordinates with hard-coded
    ! values here. the advantage is that we can set   BodyType=jerry;  and voilà!
    Insect%R_head = 0.125d0
    Insect%R_eye = 0.0625d0
    Insect%x_pivot_r_b =(/ 0.05d0, -0.2165d0, 0.d0 /)
    Insect%x_pivot_l_b =(/ 0.05d0, +0.2165d0, 0.d0 /)
    Insect%x_pivot_r2_b =(/ 0.d0, -0.d0, 0.d0 /)
    Insect%x_pivot_l2_b =(/ 0.d0, +0.d0, 0.d0 /)
    Insect%b_body = 0.1d0
    Insect%L_body = 1.0d0
    Insect%x_head = (/0.5d0*Insect%L_body,0.d0,0.d0 /)
    Insect%x_eye_r = Insect%x_head+dsin(45.d0*pi/180.d0)*Insect%R_head&
    *0.8d0*(/1.d0,+1.d0,1.d0/)
    Insect%x_eye_l = Insect%x_head+dsin(45.d0*pi/180.d0)*Insect%R_head&
    *0.8d0*(/1.d0,-1.d0,1.d0/)

    a_body = Insect%L_body / 2.d0
    !-----------------------------------------------------------------------------
    ! Jerry's body is an ellipsoid
    !-----------------------------------------------------------------------------
    do iz = g, size(mask,3)-1-g
        x_glob(3) = xx0(3) + dble(iz)*ddx(3) - Insect%xc_body_g(3)
        do iy = g, size(mask,2)-1-g
            x_glob(2) = xx0(2) + dble(iy)*ddx(2) - Insect%xc_body_g(2)
            do ix = g, size(mask,1)-1-g
                x_glob(1) = xx0(1) + dble(ix)*ddx(1) - Insect%xc_body_g(1)
                if (periodic_insect) x_glob = periodize_coordinate(x_glob, (/xl,yl,zl/))

                ! x_body is in the body coordinate system, which is centered at Insect%xc_body_g
                x_body = matmul( M_body, x_glob)

                ! check if inside the surrounding box (save comput. time)
                if ( dabs(x_body(2)) <= Insect%b_body + Insect%safety ) then
                    if ( dabs(x_body(3)) <= Insect%b_body + Insect%safety ) then
                        ! check for length inside ellipsoid:
                        if ( dabs(x_body(1) ) < Insect%L_body/2 + Insect%safety ) then
                            R  = dsqrt ( x_body(2)**2 + x_body(3)**2 )
                            ! this gives the R(x) shape
                            if ( (x_body(1)/a_body)**2 <= 1.d0) then
                                R0 = dsqrt( Insect%b_body**2 *(1.d0- (x_body(1)/a_body)**2 ) )

                                if ( R < R0 + Insect%safety ) then
                                    mask(ix,iy,iz)= max(steps(R,R0, Insect%smooth),mask(ix,iy,iz))
                                    mask_color(ix,iy,iz) = color_body
                                endif
                            endif
                        endif
                    endif
                endif
            enddo
        enddo
    enddo

    !-----------------------------------------------------------------------------
    ! Jerry's head and eyes are spheres
    !-----------------------------------------------------------------------------
    x_head = Insect%xc_body_g + matmul(transpose(M_body),Insect%x_head)
    call drawsphere( x_head,Insect%R_head,xx0, ddx, mask, mask_color, us,Insect,color_body )

    x_eye = Insect%xc_body_g + matmul(transpose(M_body),Insect%x_eye_l)
    call drawsphere( x_eye,Insect%R_eye,xx0, ddx, mask, mask_color, us,Insect,color_body )

    x_eye = Insect%xc_body_g + matmul(transpose(M_body),Insect%x_eye_r)
    call drawsphere( x_eye,Insect%R_eye,xx0, ddx, mask, mask_color, us,Insect,color_body )
end subroutine draw_body_jerry


! a body that is just a sphere of unit diameter. used for particles.
subroutine draw_body_sphere( xx0, ddx, mask, mask_color, us, Insect)
    implicit none

    type(diptera),intent(inout) :: Insect
    real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
    real(kind=rk),intent(inout) :: mask(0:,0:,0:)
    real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
    integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)

    real(kind=rk) :: x,R0,R,R_tmp,x_tmp,a_body
    real(kind=rk) :: corner
    real(kind=rk) :: x_body(1:3), x_glob(1:3), x_head(1:3), x_eye(1:3)
    real(kind=rk)   :: M_body(1:3,1:3)
    integer(kind=2) :: color_body

    color_body = Insect%color_body
    M_body     = Insect%M_body

    x_head = Insect%xc_body_g
    call drawsphere( x_head, Insect%L_body/2.d0, xx0, ddx, mask, mask_color, us, Insect, color_body )

end subroutine draw_body_sphere


! draw a cylinder defined by points (x1,y1,z1), (x2,y2,z2) and radius R0
subroutine draw_cylinder( xp,x1,y1,z1,x2,y2,z2,R0,mask_val,color_val,icolor,safety, h_smooth )
    implicit none

    real(kind=rk), intent(in) :: h_smooth
    real(kind=rk),intent(in)::xp(1:3),x1,x2,y1,y2,z1,z2,R0,safety
    real(kind=rk),intent(inout)::mask_val
    integer(kind=2),intent(in)::icolor
    integer(kind=2),intent(inout)::color_val

    real(kind=rk)::x(1:3),R,xab,yab,zab,xu,yu,zu,xvp,yvp,zvp,&
    cbx,cby,cbz,rbx,rby,rbz

    ! coordinates centered around cylinder mid-point
    cbx = 0.5*(x1+x2) - xp(1)
    cby = 0.5*(y1+y2) - xp(2)
    cbz = 0.5*(z1+z2) - xp(3)

    ! length of cylinder
    rbx = x1-x2
    rby = y1-y2
    rbz = z1-z2

    ! draw cylinder without endpoint treatment
    if ( cbx*cbx + cby*cby + cbz*cbz < 0.25*(rbx*rbx+rby*rby+rbz*rbz) ) then ! the 0.25 is from the 0.5 squared
        xab = xp(1) - x1
        yab = xp(2) - y1
        zab = xp(3) - z1

        ! e_x vector
        xu = x2 - x1
        yu = y2 - y1
        zu = z2 - z1

        xvp = yab*zu - zab*yu
        yvp = zab*xu - xab*zu
        zvp = xab*yu - yab*xu

        R = sqrt( (xvp*xvp + yvp*yvp + zvp*zvp) / (xu*xu + yu*yu + zu*zu) )
        if ( R <= R0+safety ) then
            mask_val = max(steps(R,R0,h_smooth),mask_val)
            color_val = icolor
        endif
    endif

    ! spheres at endpoints
    x = xp - (/x1,y1,z1/)
    if (dabs(x(1)) <= R0+safety) then
        if (dabs(x(2)) <= R0+safety) then
            if (dabs(x(3)) <= R0+safety) then
                R = dsqrt( x(1)*x(1)+x(2)*x(2)+x(3)*x(3) )
                if ( R <= R0+safety ) then
                    mask_val = max(steps(R,R0,h_smooth),mask_val)
                    color_val = icolor
                endif
            endif
        endif
    endif
    x = xp - (/x2,y2,z2/)
    if (dabs(x(1)) <= R0+safety) then
        if (dabs(x(2)) <= R0+safety) then
            if (dabs(x(3)) <= R0+safety) then
                R = dsqrt( x(1)*x(1)+x(2)*x(2)+x(3)*x(3) )
                if ( R <= R0+safety ) then
                    mask_val = max(steps(R,R0,h_smooth),mask_val)
                    color_val = icolor
                endif
            endif
        endif
    endif

end subroutine

!-------------------------------------------------------------------------------
! In our terminology, a macroscopic particle is an insect without wings and no
! flapping motion in free flight. Therefore, the insect module contains nowadays
! also body shapes that are not related to insects. This one is a flat plate of
! size
! Insect%L_span x Insect%L_body x Insect%WingThickness
!-------------------------------------------------------------------------------
subroutine draw_body_platicle( xx0, ddx, mask, mask_color, us, Insect)
    implicit none

    type(diptera),intent(inout) :: Insect
    real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
    real(kind=rk),intent(inout) :: mask(0:,0:,0:)
    real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
    integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)

    real(kind=rk) :: R0,R,a_body, projected_length
    real(kind=rk) :: x_body(1:3), x(1:3), xc(1:3), n_part(1:3)
    integer :: ix,iy,iz,ip, npoints, mpicode, ijk(1:3), box, start,i,j,k
    real(kind=rk)   :: M_body(1:3,1:3)
    integer(kind=2) :: color_body

    color_body = Insect%color_body
    M_body     = Insect%M_body

    do iz = g, size(mask,3)-1-g
        x(3) = xx0(3) + dble(iz)*ddx(3) - Insect%xc_body_g(3)
        do iy = g, size(mask,2)-1-g
            x(2) = xx0(2) + dble(iy)*ddx(2) - Insect%xc_body_g(2)
            do ix = g, size(mask,1)-1-g
                x(1) = xx0(1) + dble(ix)*ddx(1) - Insect%xc_body_g(1)
                if (periodic_insect) x = periodize_coordinate(x, (/xl,yl,zl/))

                ! x_body is in the body coordinate system
                x_body = matmul(M_body,x)

                ! bounding box checks
                if (dabs(x_body(1)) <= Insect%L_span+Insect%safety) then
                    if (dabs(x_body(2)) <= Insect%L_body+Insect%safety) then
                        if (dabs(x_body(3)) <= Insect%WingThickness+Insect%safety) then
                            ! signed distance:
                            R = maxval( (/ dabs(x_body(3))-Insect%WingThickness/2.d0,&
                            dabs(x_body(2))-Insect%L_body/2.d0,&
                            dabs(x_body(1))-Insect%L_span/2.d0 &
                            /) )
                            mask(ix,iy,iz) = max(steps(R,0.d0, Insect%smooth),mask(ix,iy,iz))
                            mask_color(ix,iy,iz) = color_body
                        endif
                    endif
                endif


            enddo
        enddo
    enddo

end subroutine draw_body_platicle


!-------------------------------------------------------------------------------
! In our terminology, a macroscopic particle is an insect without wings and no
! flapping motion in free flight. Therefore, the insect module contains nowadays
! also body shapes that are not related to insects. This one is a flat COIN (D=1)
!-------------------------------------------------------------------------------
subroutine draw_body_coin( xx0, ddx, mask, mask_color, us, Insect)
    implicit none

    type(diptera),intent(inout) :: Insect
    real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
    real(kind=rk),intent(inout) :: mask(0:,0:,0:)
    real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
    integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)

    real(kind=rk) :: R0,R,a_body, projected_length
    real(kind=rk) :: x_body(1:3), x(1:3), xc(1:3), n_part(1:3)
    integer :: ix,iy,iz,ip, npoints, mpicode, ijk(1:3), box, start,i,j,k
    real(kind=rk)   :: M_body(1:3,1:3)
    integer(kind=2) :: color_body

    color_body = Insect%color_body
    M_body     = Insect%M_body

    do iz = g, size(mask,3)-1-g
        x(3) = xx0(3) + dble(iz)*ddx(3) - Insect%xc_body_g(3)
        do iy = g, size(mask,2)-1-g
            x(2) = xx0(2) + dble(iy)*ddx(2) - Insect%xc_body_g(2)
            do ix = g, size(mask,1)-1-g
                x(1) = xx0(1) + dble(ix)*ddx(1) - Insect%xc_body_g(1)
                if (periodic_insect) x = periodize_coordinate(x, (/xl,yl,zl/))
                ! x_body is in the body coordinate system
                x_body = matmul(M_body,x)

                if (dabs(x_body(1)) <= 0.5d0+Insect%safety) then
                    if (dabs(x_body(2)) <= 0.5d0+Insect%safety) then
                        if (dabs(x_body(3)) <= Insect%WingThickness+Insect%safety) then
                            ! signed distance:
                            R = maxval( (/ dabs(x_body(3)) - Insect%WingThickness/2.d0,&
                            dsqrt(x_body(2)**2 + x_body(1)**2)-0.5d0 &
                            /) )
                            mask(ix,iy,iz) = max(steps(R,0.d0, Insect%smooth),mask(ix,iy,iz))
                            mask_color(ix,iy,iz) = color_body
                        endif
                    endif
                endif


            enddo
        enddo
    enddo

end subroutine draw_body_coin



!-------------------------------------------------------------------------------
! Thin rod-like body used in Suzuki et al. JFM 2015 to model a butterfly
!-------------------------------------------------------------------------------
subroutine draw_suzuki_thin_rod( xx0, ddx, mask, mask_color, us, Insect)
    implicit none

    type(diptera),intent(inout) :: Insect
    real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
    real(kind=rk),intent(inout) :: mask(0:,0:,0:)
    real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
    integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)

    real(kind=rk) :: R0,R,a,RR0
    real(kind=rk) :: x_body(1:3), x_glob(1:3), x_head(1:3), x_eye(1:3)
    integer :: ix,iy,iz
    real(kind=rk)   :: M_body(1:3,1:3)
    integer(kind=2) :: color_body

    color_body = Insect%color_body
    M_body     = Insect%M_body

    R0 = ( 0.5d0*Insect%WingThickness + Insect%Safety )**2
    RR0 = 0.5d0*Insect%WingThickness

    do iz = g, size(mask,3)-1-g
        x_glob(3) = xx0(3) + dble(iz)*ddx(3) - Insect%xc_body_g(3)
        do iy = g, size(mask,2)-1-g
            x_glob(2) = xx0(2) + dble(iy)*ddx(2) - Insect%xc_body_g(2)
            do ix = g, size(mask,1)-1-g
                x_glob(1) = xx0(1) + dble(ix)*ddx(1) - Insect%xc_body_g(1)
                if (periodic_insect) x_glob = periodize_coordinate(x_glob, (/xl,yl,zl/))

                ! x_body is in the body coordinate system, which is centered at Insect%xc_body_g
                x_body = matmul( M_body, x_glob)

                if ( dabs(x_body(1))<=0.5d0+Insect%safety) then
                    R = x_body(2)**2 + x_body(3)**2
                    if ( R < R0) then
                        a = steps(dsqrt(R),RR0, Insect%smooth)
                        if (mask(ix,iy,iz)<=a) then
                            mask(ix,iy,iz) = a
                            mask_color(ix,iy,iz) = color_body
                        endif
                    endif
                endif

            enddo
        enddo
    enddo

end subroutine draw_suzuki_thin_rod

!-------------------------------------------------------------------------------
subroutine draw_body_hawkmoth( xx0, ddx, mask, mask_color, us, Insect)
    implicit none

    type(diptera),intent(inout) :: Insect
    real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
    real(kind=rk),intent(inout) :: mask(0:,0:,0:)
    real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
    integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)


    real(kind=rk) :: R0,R,a_body
    real(kind=rk) :: x_body(1:3), x_glob(1:3), x_head(1:3), x_eye(1:3), x_eye_r(1:3), x_eye_l(1:3)
    real(kind=rk), dimension(1:3) :: x1,x2
    real(kind=rk)   :: M_body(1:3,1:3)
    integer(kind=2) :: color_body
    integer :: ix,iy,iz

    color_body = Insect%color_body
    M_body     = Insect%M_body

    Insect%R_head = 0.125d0
    Insect%R_eye = 0.0625d0
    ! Insect%x_pivot_r_b =(/ 0.05d0, -0.2165d0, 0.d0 /)
    ! Insect%x_pivot_l_b =(/ 0.05d0, +0.2165d0, 0.d0 /)
    Insect%b_body = 0.15d0
    Insect%L_body = 1.0d0
    Insect%x_head = (/0.5d0*Insect%L_body,0.d0,0.d0 /)
    Insect%x_eye_r = Insect%x_head+dsin(45.d0*pi/180.d0)*Insect%R_head&
    *0.8d0*(/1.d0,+1.d0,1.d0/)
    Insect%x_eye_l = Insect%x_head+dsin(45.d0*pi/180.d0)*Insect%R_head&
    *0.8d0*(/1.d0,-1.d0,1.d0/)

    a_body = Insect%L_body / 2.d0
    !-----------------------------------------------------------------------------
    ! The body is an ellipsoid
    !-----------------------------------------------------------------------------
    do iz = g, size(mask,3)-1-g
        x_glob(3) = xx0(3) + dble(iz)*ddx(3) - Insect%xc_body_g(3)
        do iy = g, size(mask,2)-1-g
            x_glob(2) = xx0(2) + dble(iy)*ddx(2) - Insect%xc_body_g(2)
            do ix = g, size(mask,1)-1-g
                x_glob(1) = xx0(1) + dble(ix)*ddx(1) - Insect%xc_body_g(1)
                if (periodic_insect) x_glob = periodize_coordinate(x_glob, (/xl,yl,zl/))

                ! x_body is in the body coordinate system, which is centered at Insect%xc_body_g
                x_body = matmul( M_body, x_glob)
                ! check if inside the surrounding box (save comput. time)
                if ( dabs(x_body(2)) <= Insect%b_body + Insect%safety ) then
                    if ( dabs(x_body(3)) <= Insect%b_body + Insect%safety ) then
                        ! check for length inside ellipsoid:
                        if ( dabs(x_body(1) ) < Insect%L_body/2 + Insect%safety ) then
                            R  = dsqrt ( x_body(2)**2 + x_body(3)**2 )
                            ! this gives the R(x) shape
                            if ( (x_body(1)/a_body)**2 <= 1.d0) then
                                R0 = dsqrt( Insect%b_body**2 *(1.d0- (x_body(1)/a_body)**2 ) )

                                if ( R < R0 + Insect%safety ) then
                                    mask(ix,iy,iz)= max(steps(R,R0, Insect%smooth),mask(ix,iy,iz))
                                    mask_color(ix,iy,iz) = color_body
                                endif
                            endif
                        endif
                    endif
                endif
            enddo
        enddo
    enddo

    !-----------------------------------------------------------------------------
    ! Head is a sphere, we add antennae, which are cylinders
    !-----------------------------------------------------------------------------
    x_head = Insect%xc_body_g + matmul(transpose(M_body),Insect%x_head)
    call drawsphere( x_head,Insect%R_head,xx0, ddx, mask, mask_color, us,Insect,color_body )

    ! these guys are in the body system:
    x_eye_r = Insect%x_head+dsin(45.d0*pi/180.d0)*Insect%R_head*0.8d0*(/1.d0,+1.d0,1.d0/)
    x_eye_l = Insect%x_head+dsin(45.d0*pi/180.d0)*Insect%R_head*1.8d0*(/1.d0,+1.d0,1.d0/)
    ! back to global system
    x1 = Insect%xc_body_g + matmul(transpose(M_body),x_eye_l)
    x2 = Insect%xc_body_g + matmul(transpose(M_body),x_eye_r)
    ! draw the cylinder (with spheres at the ends)
    call draw_cylinder_new( x1, x2, 0.015d0*1.3d0, xx0, ddx, mask, mask_color, us, Insect, color_body )


    ! these guys are in the body system:
    x_eye_r = Insect%x_head+dsin(45.d0*pi/180.d0)*Insect%R_head*0.8d0*(/1.d0,-1.d0,1.d0/)
    x_eye_l = Insect%x_head+dsin(45.d0*pi/180.d0)*Insect%R_head*1.8d0*(/1.d0,-1.d0,1.d0/)
    ! back to global system
    x1 = Insect%xc_body_g + matmul(transpose(M_body),x_eye_l)
    x2 = Insect%xc_body_g + matmul(transpose(M_body),x_eye_r)
    ! draw the cylinder (with spheres at the ends)
    call draw_cylinder_new( x1, x2, 0.015d0*1.3d0, xx0, ddx, mask, mask_color, us, Insect, color_body )
end subroutine draw_body_hawkmoth



!-------------------------------------------------------------------------------
! The mosquito is based on the simplified model presented in
! [1] Iams "Flight stability of mosquitos: A reduced model" SIAM J. Appl. Math. 74(5) 1535--1550 (2014)
!-------------------------------------------------------------------------------
subroutine draw_body_mosquito_iams( xx0, ddx, mask, mask_color, us, Insect)
    implicit none

    type(diptera),intent(inout) :: Insect
    real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
    real(kind=rk),intent(inout) :: mask(0:,0:,0:)
    real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
    integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)

    real(kind=rk) :: R0,R,a_body, a,b,c, alpha, Ralpha(1:3,1:3)
    real(kind=rk) :: x_body(1:3), x_glob(1:3), x_head(1:3), x_eye(1:3), x_eye_r(1:3), x_eye_l(1:3)
    real(kind=rk) :: x0_head(1:3), x0_abdomen(1:3), x0_thorax(1:3)
    real(kind=rk), dimension(1:3) :: x1,x2
    real(kind=rk)   :: M_body(1:3,1:3)
    integer(kind=2) :: color_body
    integer :: ix,iy,iz

    color_body = Insect%color_body
    M_body     = Insect%M_body

    ! The mosquito consists of three parts: head, thorax and abdomen (sphere, ellipsoid, ellipsoid)
    ! positions are measured from fig. 1 in [1], we computed also the center of gravity
    ! for this mosquito, Insect%xc_body_g is thus the center of gravity
    x0_head = (/ 0.5652d0, 0.d0, -0.0434d0 /)
    x0_thorax = (/ 0.2579d0, 0.d0, 0.1267d0 /)
    x0_abdomen = (/-0.437d0, 0.d0, -0.2024d0 /)

    !-----------------------------------------------------------------------------
    ! head
    !-----------------------------------------------------------------------------
    ! the head is a simple sphere with radius 0.1154
    R0 = 0.1154d0
    x1 = x0_head + Insect%xc_body_g
    call drawsphere( x1, R0, xx0, ddx, mask, mask_color, us,Insect,color_body )

    !-----------------------------------------------------------------------------
    ! thorax
    !-----------------------------------------------------------------------------
    ! the thorax is a triaxial ellipsiod without rotation
    a = 0.2628d0
    b = 0.1603d0
    c = b ! HACK: for simplicity, assume b=c, otherwise it can be very tough to draw

    do iz = g, size(mask,3)-1-g
        x_glob(3) = xx0(3) + dble(iz)*ddx(3) - Insect%xc_body_g(3)
        do iy = g, size(mask,2)-1-g
            x_glob(2) = xx0(2) + dble(iy)*ddx(2) - Insect%xc_body_g(2)
            do ix = g, size(mask,1)-1-g
                x_glob(1) = xx0(1) + dble(ix)*ddx(1) - Insect%xc_body_g(1)
                if (periodic_insect) x_glob = periodize_coordinate(x_glob, (/xl,yl,zl/))

                ! x_body is in the body coordinate system, which is centered at Insect%xc_body_g
                x_body = matmul( M_body, x_glob)
                ! translate to origin of thorax
                x_body = x_body - x0_thorax

                ! check if inside the surrounding box (save comput. time)
                if ( dabs(x_body(1)) <= b + Insect%safety ) then
                    if ( dabs(x_body(2)) <= b + Insect%safety ) then
                        if ( dabs(x_body(3)) <= a + Insect%safety ) then
                            ! the x-y plane are circles
                            R  = dsqrt ( x_body(1)**2 + x_body(2)**2 )
                            ! this gives the R(x) shape
                            if ( x_body(3)/a <= 1.d0) then
                                R0 = b * dsqrt( 1.d0 - (x_body(3)/a)**2 )
                                if ( R < R0 + Insect%safety ) then
                                    mask(ix,iy,iz)= max(steps(R,R0, Insect%smooth),mask(ix,iy,iz))
                                    mask_color(ix,iy,iz) = color_body
                                endif
                            endif
                        endif
                    endif
                endif
            enddo
        enddo
    enddo

    !-----------------------------------------------------------------------------
    ! abdomen
    !-----------------------------------------------------------------------------
    ! the abdomen is a axi-symmetric ellipsiod inclined by 30.44°
    a = 0.6026d0
    b = 0.1282d0
    ! angle by which the abdomen is tilted (measured from figure 1 in [1])
    alpha = deg2rad(-30.44d0)

    do iz = g, size(mask,3)-1-g
        x_glob(3) = xx0(3) + dble(iz)*ddx(3) - Insect%xc_body_g(3)
        do iy = g, size(mask,2)-1-g
            x_glob(2) = xx0(2) + dble(iy)*ddx(2) - Insect%xc_body_g(2)
            do ix = g, size(mask,1)-1-g
                x_glob(1) = xx0(1) + dble(ix)*ddx(1) - Insect%xc_body_g(1)
                if (periodic_insect) x_glob = periodize_coordinate(x_glob, (/xl,yl,zl/))

                ! x_body is in the body coordinate system, which is centered at Insect%xc_body_g
                x_body = matmul(M_body, x_glob)
                ! translate to origin of abdomen
                x_body = x_body - x0_abdomen
                ! rotate into abdomens principal axis
                call Ry(Ralpha,alpha)
                x_body = matmul(Ralpha, x_body)

                ! check if inside the surrounding box (save comput. time)
                if ( dabs(x_body(1)) <= a + Insect%safety ) then
                    if ( dabs(x_body(2)) <= b + Insect%safety ) then
                        if ( dabs(x_body(3)) <= b + Insect%safety ) then
                            ! the y-z plane are circles
                            R  = dsqrt ( x_body(2)**2 + x_body(3)**2 )
                            ! this gives the R(x) shape
                            if ( x_body(1)/a <= 1.d0) then
                                R0 = b * dsqrt( 1.d0 - (x_body(1)/a)**2 )
                                if ( R < R0 + Insect%safety ) then
                                    mask(ix,iy,iz)= max(steps(R,R0, Insect%smooth),mask(ix,iy,iz))
                                    mask_color(ix,iy,iz) = color_body
                                endif
                            endif
                        endif
                    endif
                endif
            enddo
        enddo
    enddo



end subroutine draw_body_mosquito_iams


!-------------------------------------------------------------------------------
! The flying pyramid, optimistically termed "bug" from the paper
! Liu, Ristroph, Weathers, Childress, Zhang, Intrinsic Stability of a Body hovering in an
! oscillating airflow, Phys. Rev. Lett. 2012
! The HEIGHT is Insect%L_body
! The SIDELENGTH is INsect%b_body
! This is the conical version (so a pyramid with circular base area)
!-------------------------------------------------------------------------------
subroutine draw_body_cone( xx0, ddx, mask, mask_color, us, Insect)
    implicit none

    type(diptera),intent(inout) :: Insect
    real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
    real(kind=rk),intent(inout) :: mask(0:,0:,0:)
    real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
    integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)

    real(kind=rk) :: R0,R,a,H, alpha, thick
    real(kind=rk) :: x_body(1:3), x_glob(1:3)
    integer :: ix,iy,iz
    logical, save :: informed = .false.
    real(kind=rk)   :: M_body(1:3,1:3)
    integer(kind=2) :: color_body

    color_body = Insect%color_body
    M_body     = Insect%M_body

    ! a is the sidelength of the pyramid
    a = Insect%b_body
    ! alpha is HALF the opening angle
    alpha = Insect%eta0
    ! heigh is defined by alpha and a
    H = a * dcos(alpha)
    ! we draw a shell here, and this is its thickness
    thick = Insect%WingThickness

    if (root) then
        if (informed .eqv. .false.) then
            write(*,'("Conical flyer H=",g12.4," alpha=",g12.4," a=",g12.4)') H, alpha*180.d0/pi, a
            informed = .true.
        endif
    endif

    do iz = g, size(mask,3)-1-g
        x_glob(3) = xx0(3) + dble(iz)*ddx(3) - Insect%xc_body_g(3)
        do iy = g, size(mask,2)-1-g
            x_glob(2) = xx0(2) + dble(iy)*ddx(2) - Insect%xc_body_g(2)
            do ix = g, size(mask,1)-1-g
                x_glob(1) = xx0(1) + dble(ix)*ddx(1) - Insect%xc_body_g(1)
                if (periodic_insect) x_glob = periodize_coordinate(x_glob, (/xl,yl,zl/))

                ! x_body is in the body coordinate system, which is centered at Insect%xc_body_g
                x_body = matmul( M_body, x_glob)
                ! shift x body to the center of gravity
                x_body(3) = x_body(3) + H/3.d0

                ! check if inside the surrounding box (save comput. time)
                if ( dabs(x_body(1)) <= a + Insect%safety ) then
                    if ( dabs(x_body(2)) <= a + Insect%safety ) then
                        if ( dabs(x_body(3)) <= H*2.d0/3.d0 + Insect%safety .and. x_body(3)>-Insect%safety-H/3.d0) then
                            ! the x-y plane are circles
                            R  = dsqrt ( x_body(1)**2 + x_body(2)**2 )
                            ! this gives the R(z) shape
                            R0 = max( -a*sin(alpha)*x_body(3) + (2.d0/3.d0)*H*a*sin(alpha) , 0.d0 )
                            ! define the mask. note w shifted the system to the center of gravity
                            ! therefore -H/3 <= z <= 2H/3
                            mask(ix,iy,iz)= max(steps(dabs(R-R0),0.5d0*thick, Insect%smooth)&
                            *steps(x_body(3),H*2.d0/3.d0, Insect%smooth)&
                            *steps(-x_body(3),H/3.d0, Insect%smooth),mask(ix,iy,iz))

                            mask_color(ix,iy,iz) = color_body
                        endif
                    endif
                endif
            enddo
        enddo
    enddo
end subroutine draw_body_cone

!-------------------------------------------------------------------------------
! draw a cylinder defined by GLOBALS points (x1,y1,z1), (x2,y2,z2) and radius R0
! At the start/end point, we add a sphere.
! The solid velocity field us is not touched -- we consider this routine for bodies
! therefore the solid velocity field (which is a solid body rotation around
! insect%xc) is added in the main insect drawing routine.
! The color of the new cylinder will be what you pass in color_val
!-------------------------------------------------------------------------------
subroutine draw_cylinder_new( x1, x2, R0, xx0, ddx, mask, mask_color, us, Insect, color_val)
    implicit none

    real(kind=rk),dimension(1:3),intent(inout )::x1,x2
    real(kind=rk),intent(in)::R0
    type(diptera),intent(inout)::Insect
    real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
    real(kind=rk),intent(inout) :: mask(0:,0:,0:)
    real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
    integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)
    integer(kind=2),intent(in) :: color_val

    real(kind=rk),dimension(1:3) ::  cb, rb, ab, u, vp
    real(kind=rk),dimension(1:3) :: x_glob, e_x, tmp, e_r, e_3
    real(kind=rk)::ceta1, ceta2, ceta3, R, RR0, clength, safety, t
    integer :: ix,iy,iz

    integer, dimension(1:3) :: lbounds, ubounds
    integer :: xmin,xmax,ymin,ymax,zmin,zmax
    integer :: Nsafety

    safety = 1.5*Insect%safety
    Nsafety = ceiling(safety / minval(ddx))

    ! bounds of the current patch of data
    lbounds = g
    ubounds = (/size(mask,1), size(mask,2), size(mask,3)/) - 1 - g

    RR0 = R0 + safety

    ! unit vector in cylinder axis direction and cylinder length
    e_x = x2 - x1
    clength = norm2(e_x)
    e_x = e_x / clength

    ! radial unit vector
    ! use a vector perpendicular to e_x, since it is a azimuthal symmetry
    ! it does not really matter which one. however, we must be sure that the vector
    ! we use and the e_x vector are not colinear -- their cross product is the zero vector, if that is the case
    e_r = (/0.d0,0.d0,0.d0/)
    do while ( norm2(e_r) <= 1.0d-12 )
        e_r = cross( (/rand_nbr(),rand_nbr(),rand_nbr()/), e_x)
    enddo
    e_r = e_r / norm2(e_r)

    ! third (also radial) unit vector, simply the cross product of the others
    e_3 = cross(e_x,e_r)
    e_3 = e_3 / norm2(e_3)

    ! bounding box of the vicinity of the cylinder.
    t = minval( (/x1(1)+RR0*e_r(1), x1(1)-RR0*e_r(1), x1(1)+RR0*e_3(1), x1(1)-RR0*e_3(1), &
                  x2(1)+RR0*e_r(1), x2(1)-RR0*e_r(1), x2(1)+RR0*e_3(1), x2(1)-RR0*e_3(1) /) )
    xmin = nint( (t-xx0(1)) / ddx(1) ) - Nsafety

    t = maxval( (/x1(1)+RR0*e_r(1), x1(1)-RR0*e_r(1), x1(1)+RR0*e_3(1), x1(1)-RR0*e_3(1), &
                  x2(1)+RR0*e_r(1), x2(1)-RR0*e_r(1), x2(1)+RR0*e_3(1), x2(1)-RR0*e_3(1) /) )
    xmax = nint( (t-xx0(1)) / ddx(1) ) + Nsafety

    t = minval( (/x1(2)+RR0*e_r(2), x1(2)-RR0*e_r(2), x1(2)+RR0*e_3(2), x1(2)-RR0*e_3(2), &
                  x2(2)+RR0*e_r(2), x2(2)-RR0*e_r(2), x2(2)+RR0*e_3(2), x2(2)-RR0*e_3(2) /) )
    ymin = nint( (t-xx0(2)) / ddx(2) ) - Nsafety

    t = maxval( (/x1(2)+RR0*e_r(2), x1(2)-RR0*e_r(2), x1(2)+RR0*e_3(2), x1(2)-RR0*e_3(2), &
                  x2(2)+RR0*e_r(2), x2(2)-RR0*e_r(2), x2(2)+RR0*e_3(2), x2(2)-RR0*e_3(2) /) )
    ymax = nint( (t-xx0(2)) / ddx(2) ) + Nsafety

    t = minval( (/x1(3)+RR0*e_r(3), x1(3)-RR0*e_r(3), x1(3)+RR0*e_3(3), x1(3)-RR0*e_3(3), &
                  x2(3)+RR0*e_r(3), x2(3)-RR0*e_r(3), x2(3)+RR0*e_3(3), x2(3)-RR0*e_3(3) /) )
    zmin = nint( (t-xx0(3)) / ddx(3) ) - Nsafety

    t = maxval( (/x1(3)+RR0*e_r(3), x1(3)-RR0*e_r(3), x1(3)+RR0*e_3(3), x1(3)-RR0*e_3(3), &
                  x2(3)+RR0*e_r(3), x2(3)-RR0*e_r(3), x2(3)+RR0*e_3(3), x2(3)-RR0*e_3(3) /) )
    zmax = nint( (t-xx0(3)) / ddx(3) ) + Nsafety


    ! first we draw the cylinder, then the endpoint spheres
    do iz = max(zmin,lbounds(3)), min(zmax,ubounds(3))
        x_glob(3) = xx0(3) + dble(iz)*ddx(3)

        do iy = max(ymin,lbounds(2)), min(ymax,ubounds(2))
            x_glob(2) = xx0(2) + dble(iy)*ddx(2)

            do ix = max(xmin,lbounds(1)), min(xmax,ubounds(1))
                x_glob(1) = xx0(1) + dble(ix)*ddx(1)
                ! if (periodic_insect) x_glob = periodize_coordinate(x_glob, (/xl,yl,zl/))

                ! cb is the distance to the cylinder mid-point
                cb = 0.5d0*(x1+x2) - x_glob
                ! rb is the length of the clinder
                rb = x1 - x2

                ! this is a spherical bounding box, centered around the mid-point
                if ( sum(cb**2) < 0.25*sum(rb**2) ) then ! the 0.25 is from the 0.5 squared
                    ab = x_glob - x1
                    u = x2 - x1

                    vp = cross(ab, u)
                    R = sqrt( sum(vp**2) / sum(u**2) )

                    if (R <= R0+safety) then
                        t = steps(R, R0, Insect%smooth)
                        if (t >= mask(ix,iy,iz)) then
                            mask(ix,iy,iz) = t
                            mask_color(ix,iy,iz) = color_val
                        endif
                    endif
                endif
            enddo
        enddo
    enddo

    !---------------------------------------------------------------------------
    ! endpoint spheres
    !---------------------------------------------------------------------------
    call drawsphere( x1, R0, xx0, ddx, mask, mask_color, us, Insect, color_val )
    call drawsphere( x2, R0, xx0, ddx, mask, mask_color, us, Insect, color_val )
end subroutine draw_cylinder_new



!-------------------------------------------------------------------------------
! draw a sphere with radius R0 centered at the point xc (GLOBAL SYSTEM)
! This routine's intended use is for drawing the insect's body, for example
! the head and eyes of Jerry. The velocity field inside the body is added
! later, thus, the field us is untouched in this routines.
!-------------------------------------------------------------------------------
subroutine drawsphere( xc, R0, xx0, ddx, mask, mask_color, us, Insect, icolor )
    implicit none

    real(kind=rk),intent(inout)::xc(1:3)
    real(kind=rk),intent(in)::R0
    real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
    real(kind=rk),intent(inout) :: mask(0:,0:,0:)
    real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
    integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)
    integer(kind=2),intent(in) :: icolor
    type(diptera),intent(inout) :: Insect

    integer :: ix,iy,iz
    real(kind=rk)::x(1:3),R,tmp
    integer, dimension(1:3) :: lbounds, ubounds
    integer :: xmin,xmax,ymin,ymax,zmin,zmax
    integer :: Nsafety

    ! periodization: if the center point is out of the domain, then correct that
    if (periodic_insect) then
        if (xc(1)<0.0) xc(1)=xc(1)+xl
        if (xc(2)<0.0) xc(2)=xc(2)+yl
        if (xc(3)<0.0) xc(3)=xc(3)+zl

        if (xc(1)>=xl) xc(1)=xc(1)-xl
        if (xc(2)>=yl) xc(2)=xc(2)-yl
        if (xc(3)>=zl) xc(3)=xc(3)-zl
    endif

    Nsafety = nint( (R0+Insect%safety) / minval(ddx))

    ! bounds of the current patch of data
    lbounds = g
    ubounds = (/size(mask,1), size(mask,2), size(mask,3)/) - 1 - g

    ! bounding box of the vicinity of the sphere.
    xmin = nint( (xc(1)-xx0(1)) / ddx(1) ) - (Nsafety)
    xmax = nint( (xc(1)-xx0(1)) / ddx(1) ) + (Nsafety)
    ymin = nint( (xc(2)-xx0(2)) / ddx(2) ) - (Nsafety)
    ymax = nint( (xc(2)-xx0(2)) / ddx(2) ) + (Nsafety)
    zmin = nint( (xc(3)-xx0(3)) / ddx(3) ) - (Nsafety)
    zmax = nint( (xc(3)-xx0(3)) / ddx(3) ) + (Nsafety)


    do iz = max(zmin,lbounds(3)), min(zmax,ubounds(3))
        x(3) = xx0(3) + dble(iz)*ddx(3) - xc(3)

        do iy = max(ymin,lbounds(2)), min(ymax,ubounds(2))
            x(2) = xx0(2) + dble(iy)*ddx(2) - xc(2)

            do ix = max(xmin,lbounds(1)), min(xmax,ubounds(1))
                x(1) = xx0(1) + dble(ix)*ddx(1) - xc(1)
                if (periodic_insect) x = periodize_coordinate(x, (/xl,yl,zl/))

                ! the bounding box check is incorporated in the loop bounds - no if clause!
                ! compute radius
                R = dsqrt( x(1)*x(1)+x(2)*x(2)+x(3)*x(3) )
                if ( R <= R0+Insect%safety ) then
                    tmp = steps(R, R0, Insect%smooth)
                    if (tmp>=mask(ix,iy,iz)) then
                        ! set new value
                        mask(ix,iy,iz) = tmp
                        mask_color(ix,iy,iz) = icolor
                    endif
                endif
            enddo
        enddo
    enddo

end subroutine



!-------------------------------------------------------------------------------
! Draw a body from SUPERSTL file.
! A superstl is just a list of triangles with all face-, edge- and vertex normals precomputed
!
! No scaling or origin shift is applied: we assume you did that when generating
! the superSTL file. The data is thus understood in the body coordinate system.
!
!-------------------------------------------------------------------------------
subroutine draw_body_superSTL(x0, dx, mask, mask_color, us, Insect)
    implicit none

    type(diptera),intent(inout) :: Insect
    real(kind=rk),intent(in)    :: x0(1:3), dx(1:3)
    real(kind=rk),intent(inout) :: mask(0:,0:,0:)
    real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
    integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)

    real(kind=rk), allocatable, save :: tmp_block(:,:,:)

    integer :: safety, i, ntri, Bs(1:3)
    integer :: ix, iy, iz, ivertex, xmin, xmax, ymin, ymax, zmin, zmax
    integer(kind=2) :: color_body

    real(kind=rk), dimension(1:3) :: vertex1, vertex2, vertex3, vertex1_normal
    real(kind=rk), dimension(1:3) :: vertex2_normal, vertex3_normal, face_normal, edge1_normal, edge2_normal, edge3_normal
    real(kind=rk) :: M_body(1:3,1:3), M_body_inv(1:3,1:3)
    real(kind=rk) :: scale, origin(1:3), tmp, x, y, z
    real(kind=rk) :: x_glob(1:3), x_body(1:3)
    character(len=strlen) :: fname_stl

    color_body = Insect%color_body
    M_body     = Insect%M_body
    M_body_inv = transpose(Insect%M_body)
    fname_stl  = Insect%BodySuperSTLfile

    ! we work on a work block because we add the body to existing masks
    if (.not. allocated(tmp_block)) then
        allocate( tmp_block(0:size(mask,1)-1, 0:size(mask,2)-1, 0:size(mask,3)-1) )
    endif

    Bs(1) = size(mask,1) - 2*g
    Bs(2) = size(mask,2) - 2*g
    Bs(3) = size(mask,3) - 2*g

    ! initialize signed distance as very far away
    tmp_block = 9.0e6_rk ! distance as far away
    safety = 6
    ntri = size(xyz_nxnynz, 1)

    ! loop over all triangles
    do i = 1, ntri
        vertex1        = matmul(M_body_inv, xyz_nxnynz(i, 1:3)) + Insect%xc_body_g
        vertex2        = matmul(M_body_inv, xyz_nxnynz(i, 4:6)) + Insect%xc_body_g
        vertex3        = matmul(M_body_inv, xyz_nxnynz(i, 7:9)) + Insect%xc_body_g
        face_normal    = matmul(M_body_inv, xyz_nxnynz(i, 10:12))
        vertex1_normal = matmul(M_body_inv, xyz_nxnynz(i, 13:15))
        vertex2_normal = matmul(M_body_inv, xyz_nxnynz(i, 16:18))
        vertex3_normal = matmul(M_body_inv, xyz_nxnynz(i, 19:21))
        edge1_normal   = matmul(M_body_inv, xyz_nxnynz(i, 22:24))
        edge2_normal   = matmul(M_body_inv, xyz_nxnynz(i, 25:27))
        edge3_normal   = matmul(M_body_inv, xyz_nxnynz(i, 28:30))

        xmin = floor( ( minval((/vertex1(1),vertex2(1),vertex3(1)/)) - x0(1) ) / dx(1)) - safety
        ymin = floor( ( minval((/vertex1(2),vertex2(2),vertex3(2)/)) - x0(2) ) / dx(2)) - safety
        zmin = floor( ( minval((/vertex1(3),vertex2(3),vertex3(3)/)) - x0(3) ) / dx(3)) - safety

        xmax = ceiling( ( maxval((/vertex1(1),vertex2(1),vertex3(1)/)) - x0(1) ) / dx(1)) + safety
        ymax = ceiling( ( maxval((/vertex1(2),vertex2(2),vertex3(2)/)) - x0(2) ) / dx(2)) + safety
        zmax = ceiling( ( maxval((/vertex1(3),vertex2(3),vertex3(3)/)) - x0(3) ) / dx(3)) + safety

        ! note these guys are zero based indexing....
        xmin = max(xmin, 0)
        ymin = max(ymin, 0)
        zmin = max(zmin, 0)

        xmax = min(xmax, size(mask,1)-1)
        ymax = min(ymax, size(mask,2)-1)
        zmax = min(zmax, size(mask,3)-1)


        vertex1        = xyz_nxnynz(i, 1:3)
        vertex2        = xyz_nxnynz(i, 4:6)
        vertex3        = xyz_nxnynz(i, 7:9)
        face_normal    = xyz_nxnynz(i, 10:12)
        vertex1_normal = xyz_nxnynz(i, 13:15)
        vertex2_normal = xyz_nxnynz(i, 16:18)
        vertex3_normal = xyz_nxnynz(i, 19:21)
        edge1_normal   = xyz_nxnynz(i, 22:24)
        edge2_normal   = xyz_nxnynz(i, 25:27)
        edge3_normal   = xyz_nxnynz(i, 28:30)

        do iz = zmin, zmax
            x_glob(3) = x0(3) + dble(iz)*dx(3) - Insect%xc_body_g(3)
            do iy = ymin, ymax
                x_glob(2) = x0(2) + dble(iy)*dx(2) - Insect%xc_body_g(2)
                do ix = xmin, xmax
                    x_glob(1) = x0(1) + dble(ix)*dx(1) - Insect%xc_body_g(1)

                    if (periodic_insect) x_glob = periodize_coordinate(x_glob, (/xl,yl,zl/))
                    ! x_body is in the body coordinate system
                    x_body = matmul(M_body, x_glob)

                    ! the distance to the current triangle:
                    tmp = pointTriangleDistance( vertex1, vertex2, vertex3, x_body, face_normal, &
                    vertex1_normal, vertex2_normal, vertex3_normal, edge1_normal, edge2_normal, edge3_normal)

                    ! if closer (in abs value!) then use this now
                    if ( abs(tmp) < abs(tmp_block(ix,iy,iz)) ) then !.and. abs(tmp)<=dble(safety)*dx(1) ) then
                        tmp_block(ix,iy,iz) = tmp
                    endif
                enddo
            enddo
        enddo
    enddo ! loop over triangles
    !
    ! signed distance to mask function
    do iz = 0, size(mask,3)-1
        do iy = 0, size(mask,2)-1
            do ix = 0, size(mask,1)-1
                tmp = smoothstep( tmp_block(ix,iy,iz), 0.0_rk, 1.5_rk*dx(1) )

                if (mask(ix,iy,iz) <= tmp) then
                    mask(ix,iy,iz) = tmp
                    mask_color(ix,iy,iz) = color_body
                endif
            enddo
        enddo
    enddo

    ! mask = tmp_block

    ! ! store final resut, do not erase existing data in array
    ! where (mask <= tmp_block)
    !     mask = tmp_block
    !     mask_color = color_body
    ! end where
end subroutine
