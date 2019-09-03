
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
        where (mask_color==Insect%color_body)
            mask = 0.00_rk
            us(:,:,:,1) = 0.00_rk
            us(:,:,:,2) = 0.00_rk
            us(:,:,:,3) = 0.00_rk
            mask_color = 0
        end where
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
    case ("jerry","Jerry")
        call draw_body_jerry( xx0, ddx, mask, mask_color, us, Insect)
    case ("hawkmoth","Hawkmoth")
        call draw_body_hawkmoth( xx0, ddx, mask, mask_color, us, Insect)
    case ("particle")
        call draw_body_particle( xx0, ddx, mask, mask_color, us, Insect)
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
    case ("pyramid")
        call draw_body_pyramid( xx0, ddx, mask, mask_color, us, Insect)
    case ("cone")
        call draw_body_cone( xx0, ddx, mask, mask_color, us, Insect)
    case ("birch_seed")
        call draw_birch_seed( xx0, ddx, mask, mask_color, us, Insect)
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
subroutine draw_body_particle( xx0, ddx, mask, mask_color, us, Insect)
    implicit none

    type(diptera),intent(inout) :: Insect
    real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
    real(kind=rk),intent(inout) :: mask(0:,0:,0:)
    real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
    integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)

    real(kind=rk) :: R0,R,a_body, projected_length
    real(kind=rk) :: x_body(1:3), x_glob(1:3), x_part(1:3), n_part(1:3)
    integer :: ix,iy,iz,ip, npoints, mpicode, ijk(1:3), box, start,i,j,k
    real(kind=rk)   :: M_body(1:3,1:3)
    integer(kind=2) :: color_body

    color_body = Insect%color_body
    M_body     = Insect%M_body

    ! NOTE HACK
    ! this routine is deprecated, its functionality is in the STL and pointcloud ideas.

    !   !-----------------------------------------------------------------------------
    !   ! initialization phase, executed only once
    !   !-----------------------------------------------------------------------------
    !   if (.not.allocated(particle_points)) then
    !     ! initialization, this is the first call to the routine.
    !     if (root) then
    !       open(37, file='particle.in', form='formatted', status='old')
    !       read(37,*) npoints
    !       write(*,'("reading particle with ",i7," points")') npoints
    !     endif
    !     call MPI_BCAST( npoints,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpicode )
    !     allocate ( particle_points(1:npoints,1:6) )
    !
    !     if (root) then
    !       do ix = 1,npoints
    !         read(37,*) particle_points(ix,:)
    !       enddo
    !       write(*,*) "done reading particle.in (that is excellent news!)"
    !       close(37)
    !     endif
    !     ! make sure all cpu know the particle well
    !     call MPI_BCAST( particle_points(:,1),npoints,MPI_DOUBLE_PRECISION,0,&
    !     MPI_COMM_WORLD,mpicode )
    !     call MPI_BCAST( particle_points(:,2),npoints,MPI_DOUBLE_PRECISION,0,&
    !     MPI_COMM_WORLD,mpicode )
    !     call MPI_BCAST( particle_points(:,3),npoints,MPI_DOUBLE_PRECISION,0,&
    !     MPI_COMM_WORLD,mpicode )
    !     call MPI_BCAST( particle_points(:,4),npoints,MPI_DOUBLE_PRECISION,0,&
    !     MPI_COMM_WORLD,mpicode )
    !     call MPI_BCAST( particle_points(:,5),npoints,MPI_DOUBLE_PRECISION,0,&
    !     MPI_COMM_WORLD,mpicode )
    !     call MPI_BCAST( particle_points(:,6),npoints,MPI_DOUBLE_PRECISION,0,&
    !     MPI_COMM_WORLD,mpicode )
    !   endif
    !
    !   ! initialize signed distance as very far away
    !   mask = 1.d8
    !
    !   !-----------------------------------------------------------------------------
    !   ! now, we are sure that all CPU know all points on the particle, and we can
    !   ! proceed to draw it
    !   !-----------------------------------------------------------------------------
    !   npoints = size(particle_points,1)
    !   ! loop over the marker points
    !   do ip = 1, npoints
    !     ! coordinate of surface point (in body system)
    !     x_part = particle_points(ip,1:3)
    !     ! normal vector of surface point (in body system)
    !     n_part = particle_points(ip,4:6)
    !
    !     ! go to laboratory frame:
    !     x_glob = matmul( transpose(M_body), x_part) + Insect%xc_body_g
    !     ! periodize:
    !     if (x_glob(1)<0.0) x_glob(1)=x_glob(1)+xl
    !     if (x_glob(2)<0.0) x_glob(2)=x_glob(2)+yl
    !     if (x_glob(3)<0.0) x_glob(3)=x_glob(3)+zl
    !     if (x_glob(1)>=xl) x_glob(1)=x_glob(1)-xl
    !     if (x_glob(2)>=yl) x_glob(2)=x_glob(2)-yl
    !     if (x_glob(3)>=zl) x_glob(3)=x_glob(3)-zl
    !
    !     ! coordinate (global) in integer space:
    !     ijk = nint( x_glob/dx )
    !     ! size of box around markers
    !     box = 1
    !
    !     ! loop over neigborhood of marker
    !     do iz = ijk(3)-box, ijk(3)+box
    !       do iy = ijk(2)-box, ijk(2)+box
    !         do ix = ijk(1)-box, ijk(1)+box
    !           ! check if this point is on my rank. note this also checks implicitly
    !           ! if we're in the domain at all
    ! !          if ( on_proc( (/ix,iy,iz/)  ) ) then
    !             ! x_glob is in the global coordinate system
    !             x_glob = (/ xx0(1)+dble(ix)*ddx(1), xx0(2)+dble(iy)*ddx(2), xx0(3)+dble(iz)*ddx(3) /)
    !             x_glob = periodize_coordinate(x_glob - Insect%xc_body_g, (/xl,yl,zl/))
    !             ! x_body is in the body coordinate system, which is centered at Insect%xc_body_g
    !             x_body = matmul( M_body, x_glob )
    !
    !             ! unsigned distance to point
    !             R = dsqrt(  (x_body(1)-x_part(1))*(x_body(1)-x_part(1)) + &
    !                         (x_body(2)-x_part(2))*(x_body(2)-x_part(2)) + &
    !                         (x_body(3)-x_part(3))*(x_body(3)-x_part(3)) )
    !
    !             if ( R<=abs(mask(ix,iy,iz)) ) then
    !               ! this is closer, so we overwrite
    !               mask(ix,iy,iz) = R
    !               mask_color(ix,iy,iz) = color_body
    !               ! compute scalar product of difference vector and outward pointing normal
    !               projected_length = (x_body(1)-x_part(1))*n_part(1) + &
    !                                  (x_body(2)-x_part(2))*n_part(2) + &
    !                                  (x_body(3)-x_part(3))*n_part(3)
    !
    !               if ( projected_length <= 0.d0  ) then
    !                 ! we're inside the particle
    !                 mask(ix,iy,iz) = -R!mask(ix,iy,iz)
    !               endif
    !             endif
    !
    ! !          endif ! on_proc
    !         enddo ! nx
    !       enddo ! ny
    !     enddo ! nz
    !   enddo ! np
    !
    !
    !   ! do iz = g, size(mask,3)-1-g
    !   !   do iy = g, size(mask,2)-1-g
    !   !     do ix = g, size(mask,1)-1-g
    !   !       if (mask(ix,iy,iz) < 0.d0) then
    !   !         R=0.d0
    !   !         do k=iz-1,iz+1
    !   !           do j=iy-1,iy+1
    !   !             do i=ix-1,ix+1
    !   !               R=R+mask( per(i,nx),per(j,ny),per(k,nz) )
    !   !             enddo
    !   !           enddo
    !   !         enddo
    !   !         R=R-mask(ix,iy,iz)
    !   !         if (R>26.5e8) mask(ix,iy,iz) = 1.0d8
    !   !
    !   !       endif
    !   !     enddo
    !   !   enddo
    !   ! enddo
    !
    !   !-----------------------------------------------------------------------------
    !   ! fill the interior of the particle with "-" signs (the above concentrates
    !   ! on the interface!)
    !   ! we exploit the fact that the x-direction is continuous in memory and not
    !   ! split among processes
    !   !-----------------------------------------------------------------------------
    !   ! we start at a point which is surely not inside
    !   ! the particle. the algorithm cannot start on interior
    !   ! points.
    !   ! start = per(nint( (Insect%xc_body_g(1)-0.5*xl)/dx), nx)
    !   ! if(root) write(*,*) "point is", Insect%xc_body_g(1)-0.5*xl
    ! start =0
    !   do iz = g, size(mask,3)-1-g
    !     do iy = g, size(mask,2)-1-g
    !       do ix = start, start+nx ! we run all points, still
    !         ! is the point not yet touched (i.e. large value)?
    !         if (mask(per(ix,nx),iy,iz) > 99.99e6) then
    !           ! is either of the neighbors inside (e.g. negative)
    !           if (mask(per(ix+1,nx),iy,iz)<0.d0.or.mask(per(ix-1,nx),iy,iz)<0.d0) then
    !             mask(per(ix,nx),iy,iz) = -mask(per(ix,nx),iy,iz)
    !             mask_color(per(ix,nx),iy,iz) = color_body
    !           endif
    !         endif
    !       enddo
    !     enddo
    !   enddo
    !
    !   ! call save_field_hdf5(0.d0,'./mask_00',mask)
    !   ! stop
    !
    !   !-----------------------------------------------------------------------------
    !   ! convert signed distance function to mask function chi
    !   !-----------------------------------------------------------------------------
    !   do iz = g, size(mask,3)-1-g
    !     do iy = g, size(mask,2)-1-g
    !       do ix = g, size(mask,1)-1-g
    !         mask(ix,iy,iz) = steps( mask(ix,iy,iz),0.d0 )
    !       enddo
    !     enddo
    !   enddo


end subroutine draw_body_particle


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


subroutine draw_birch_seed( xx0, ddx, mask, mask_color, us, Insect)
    implicit none

    type(diptera),intent(inout) :: Insect
    real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
    real(kind=rk),intent(inout) :: mask(0:,0:,0:)
    real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
    integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)

    real(kind=rk) :: R0,R,a,H, alpha, thick, a_body
    real(kind=rk) :: x_body(1:3), x_glob(1:3), Id(1:3,1:3)
    integer :: ix,iy,iz
    logical, save :: informed = .false.
    real(kind=rk)   :: M_body(1:3,1:3)
    integer(kind=2) :: color_body

    color_body = Insect%color_body
    M_body     = Insect%M_body

    ! if (root ) then
    !    if (informed .eqv. .false. ) then
    !     write(*,'("Conical flyer H=",g12.4," alpha=",g12.4," a=",g12.4)') H, alpha*180.d0/pi, a
    !     informed = .true.
    !   endif
    ! endif

    Insect%b_body = 0.09d0
    a_body = 0.17d0

    !-----------------------------------------------------------------------------
    ! The seed's core is an ellipsoid
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

    Id = 0.d0
    Id(1,1) = 1.d0; Id(2,2) = 1.d0; Id(3,3) = 1.0d0

    call draw_wing_fourier(xx0, ddx, mask, mask_color, us,Insect,color_body,M_body,Id,(/0.d0,0.d0,0.d0/),(/0.d0,0.d0,0.d0/))

end subroutine draw_birch_seed

!-------------------------------------------------------------------------------
! The flying pyramid, optimistically termed "bug" from the paper
! Liu, Ristroph, Weathers, Childress, Zhang, Intrinsic Stability of a Body hovering in an
! oscillating airflow, Phys. Rev. Lett. 2012
! The HEIGHT is Insect%L_body
! The SIDELENGTH is INsect%b_body
!-------------------------------------------------------------------------------
subroutine draw_body_pyramid( xx0, ddx, mask, mask_color, us, Insect)
    implicit none

    type(diptera),intent(inout) :: Insect
    real(kind=rk),intent(in) :: xx0(1:3), ddx(1:3)
    real(kind=rk),intent(inout) :: mask(0:,0:,0:)
    real(kind=rk),intent(inout) :: us(0:,0:,0:,1:)
    integer(kind=2),intent(inout) :: mask_color(0:,0:,0:)

    real(kind=rk) :: R0,R,a,H,b, alpha, thick, di(1:4)
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
            write(*,'("Pyramid flyer H=",g12.4," alpha=",g12.4," a=",g12.4)') H, alpha*180.d0/pi, a
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
                if ( dabs(x_body(1)) <= a/2.d0 + Insect%safety ) then
                    if ( dabs(x_body(2)) <= a/2.d0 + Insect%safety ) then
                        !  note we shifted the system to the center of gravity therefore -H/3 <= z <= 2H/3
                        if ( dabs(x_body(3)) <= H*2.d0/3.d0 + Insect%safety .and. x_body(3)>-Insect%safety-H/3.d0) then
                            ! realizing that the pyramid consists of 4 triangles, we can use the pointTriangleDistance function from
                            ! the stl_file_reader module to compute, for all points within the bounding box, the distance to the
                            ! pyramid. This is more elegant than the old solution, since it nicely smoothes the surface and edges.
                            ! note we do use the NORMAL only to distinguish between in/out, it's not required to be correct
                            di(1) = pointTriangleDistance2( (/+a/2.d0,-a/2.d0,-H/3.d0/), (/+a/2.d0,+a/2.d0,-H/3.d0/) ,&
                            (/0.d0,0.d0,2.d0*H/3.d0/), x_body, (/1.d0,0.d0,0.d0/) )
                            di(2) = pointTriangleDistance2( (/+a/2.d0,+a/2.d0,-H/3.d0/), (/-a/2.d0,+a/2.d0,-H/3.d0/) ,&
                            (/0.d0,0.d0,2.d0*H/3.d0/), x_body, (/0.d0,1.d0,0.d0/) )
                            di(3) = pointTriangleDistance2( (/-a/2.d0,+a/2.d0,-H/3.d0/), (/-a/2.d0,-a/2.d0,-H/3.d0/) ,&
                            (/0.d0,0.d0,2.d0*H/3.d0/), x_body, (/-1.d0,0.d0,0.d0/) )
                            di(4) = pointTriangleDistance2( (/-a/2.d0,-a/2.d0,-H/3.d0/), (/+a/2.d0,-a/2.d0,-H/3.d0/) ,&
                            (/0.d0,0.d0,2.d0*H/3.d0/), x_body, (/0.d0,-1.d0,0.d0/) )

                            ! note maxval(di) is the union of all the triangles' interiors (thus a solid filled pyramid results),
                            ! and abs(maxval(di))-t/2 is the shell operator which gives us just the shell (negative values are inside the solid, by convention)
                            ! we also directly take the CHI(delta) here.
                            mask(ix,iy,iz)= steps( abs(maxval(di)) -thick/2.d0, 0.d0, Insect%smooth)
                            mask_color(ix,iy,iz) = color_body
                            ! please note that the solid body velocity will be added elsewhere in the code
                        endif
                    endif
                endif
            enddo
        enddo
    enddo
end subroutine draw_body_pyramid


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

    safety = Insect%safety
    Nsafety = nint(safety / minval(ddx))

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
    if (xc(1)<0.0) xc(1)=xc(1)+xl
    if (xc(2)<0.0) xc(2)=xc(2)+yl
    if (xc(3)<0.0) xc(3)=xc(3)+zl

    if (xc(1)>=xl) xc(1)=xc(1)-xl
    if (xc(2)>=yl) xc(2)=xc(2)-yl
    if (xc(3)>=zl) xc(3)=xc(3)-zl

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


! This function is used in body type pyramid
real(kind=rk) function pointTriangleDistance2(tri1,tri2,tri3,point,normal)
    ! calculate distance between a point and a triangle in 3D
    ! SYNTAX
    !   dist = pointTriangleDistance2(TRI,P)
    !   [dist,PP0] = pointTriangleDistance2(TRI,P)
    !
    ! DESCRIPTION
    !   Calculate the distance of a given point P from a triangle TRI.
    !   Point P is a row vector of the form 1x3. The triangle is a matrix
    !   formed by three rows of points TRI = [P1P2P3] each of size 1x3.
    !   dist = pointTriangleDistance2(TRI,P) returns the distance of the point P
    !   to the triangle TRI.
    !   [dist,PP0] = pointTriangleDistance2(TRI,P) additionally returns the
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
    real(kind=rk), dimension(1:3), intent(in) :: tri1,tri2,tri3,point,normal
    real(kind=rk), dimension(1:3) :: BB,EE0,EE1,DD
    real(kind=rk) :: a,b,c,d,e,f,det,s,t,sqrDistance,tmp0,tmp1,numer,denom,invDet

    ! rewrite triangle in normal form
    BB = tri1
    EE0 = tri2-BB
    EE1 = tri3-BB


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

    if (det < 1.0d-11) then
        pointTriangleDistance2 = 9.0d9
        return
    endif


    ! write(*,'(12(es12.4,1x))') tri1,tri2,tri3,point
    ! write(*,'(12(es12.4,1x))') a,b,c,d,e,f,det,s,t

    ! Terible tree of conditionals to determine in which region of the diagram
    ! shown above the projection of the point into the triangle-plane lies.
    if ((s+t) <= det) then
        if (s < 0.d0) then
            if (t < 0.d0) then
                !region4
                if (d < 0.d0) then
                    t = 0.d0
                    if (-d >= a) then
                        s = 1.d0
                        sqrDistance = a + 2.d0*d + f
                    else
                        s = -d/a
                        sqrDistance = d*s + f
                    endif
                else
                    s = 0.d0
                    if (e >= 0.d0) then
                        t = 0.d0
                        sqrDistance = f
                    else
                        if (-e >= c) then
                            t = 1.d0
                            sqrDistance = c + 2.d0*e + f
                        else
                            t = -e/c
                            sqrDistance = e*t + f
                        endif
                    endif
                endif !of region 4
            else
                ! region 3
                s = 0.d0
                if (e >= 0.d0) then
                    t = 0.d0
                    sqrDistance = f
                else
                    if (-e >= c) then
                        t = 1.d0
                        sqrDistance = c + 2.d0*e +f
                    else
                        t = -e/c
                        sqrDistance = e*t + f
                    endif
                endif
            endif !of region 3
        else
            if (t < 0.d0) then
                ! region 5
                t = 0.d0
                if (d >= 0.d0) then
                    s = 0.d0
                    sqrDistance = f
                else
                    if (-d >= a) then
                        s = 1.d0
                        sqrDistance = a + 2.d0*d + f! GF 20101013 fixed typo d*s ->2*d
                    else
                        s = -d/a
                        sqrDistance = d*s + f
                    endif
                endif
            else
                ! region 0
                invDet = 1.d0/det
                s = s*invDet
                t = t*invDet
                sqrDistance = s*(a*s + b*t + 2.d0*d) &
                + t*(b*s + c*t + 2.d0*e) + f
            endif
        endif
    else
        if (s < 0.d0) then
            ! region 2
            tmp0 = b + d
            tmp1 = c + e
            if (tmp1 > tmp0) then ! minimum on edge s+t=1
                numer = tmp1 - tmp0
                denom = a - 2.d0*b + c
                if (numer >= denom) then
                    s = 1.d0
                    t = 0.d0
                    sqrDistance = a + 2.d0*d + f ! GF 20101014 fixed typo 2*b -> 2*d
                else
                    s = numer/denom
                    t = 1.d0-s
                    sqrDistance = s*(a*s + b*t + 2.d0*d) &
                    + t*(b*s + c*t + 2.d0*e) + f
                endif
            else          ! minimum on edge s=0
                s = 0.d0
                if (tmp1 <= 0.d0) then
                    t = 1.d0
                    sqrDistance = c + 2.d0*e + f
                else
                    if (e >= 0.d0) then
                        t = 0.d0
                        sqrDistance = f
                    else
                        t = -e/c
                        sqrDistance = e*t + f
                    endif
                endif
            endif !of region 2
        else
            if (t < 0.d0) then
                !region6
                tmp0 = b + e
                tmp1 = a + d
                if (tmp1 > tmp0) then
                    numer = tmp1 - tmp0
                    denom = a-2.d0*b+c
                    if (numer >= denom) then
                        t = 1.d0
                        s = 0.d0
                        sqrDistance = c + 2.d0*e + f
                    else
                        t = numer/denom
                        s = 1.d0 - t
                        sqrDistance = s*(a*s + b*t + 2.d0*d) &
                        + t*(b*s + c*t + 2.d0*e) + f
                    endif
                else
                    t = 0.d0
                    if (tmp1 <= 0) then
                        s = 1.d0
                        sqrDistance = a + 2.d0*d + f
                    else
                        if (d >= 0.d0) then
                            s = 0.d0
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
                numer = c + e - b - d
                if (numer <= 0.d0) then
                    s = 0.d0
                    t = 1.d0
                    sqrDistance = c + 2.d0*e + f
                else
                    denom = a - 2.d0*b + c
                    if (numer >= denom) then
                        s = 1.d0
                        t = 0.d0
                        sqrDistance = a + 2.d0*d + f
                    else
                        s = numer/denom
                        t = 1-s
                        sqrDistance = s*(a*s + b*t + 2.d0*d) &
                        + t*(b*s + c*t + 2.d0*e) + f
                    endif
                endif !of region 1
            endif
        endif
    endif

    ! account for numerical round-off error
    if (sqrDistance < 0.d0) then
        sqrDistance = 0.d0
    endif



    ! closest point on triangle
    DD = BB + s*EE0 + t*EE1;
    ! vector from target point to closest point on surface
    DD = point-DD
    t = dot_product(DD,normal)
    if (t >= 0.d0) then
        pointTriangleDistance2 = dsqrt(sqrDistance)
    else
        pointTriangleDistance2 = -dsqrt(sqrDistance)
    endif

end function
