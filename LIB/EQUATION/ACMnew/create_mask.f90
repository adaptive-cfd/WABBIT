! mask data. we can use different trees (4est module) to generate time-dependent/indenpedent
! mask functions separately. This makes the mask routines tree-level routines (and no longer
! block level) so the physics modules have to provide an interface to create the mask at a tree
! level. All parts of the mask shall be included: chi, boundary values, sponges.
! This is a block-level wrapper to fill the mask.
subroutine create_mask_3D_ACM( time, x0, dx, Bs, g, mask, stage )
    implicit none

    ! grid
    integer(kind=ik), intent(in) :: Bs(3), g
    !> mask term for every grid point of this block
    real(kind=rk), dimension(:,:,:,:), intent(inout) :: mask
    !     stage == "time-independent-part"
    !     stage == "time-dependent-part"
    !     stage == "all-parts"
    character(len=*), intent(in) :: stage
    !> spacing and origin of block
    real(kind=rk), intent(in) :: x0(1:3), dx(1:3), time

    integer(kind=2), allocatable, save :: mask_color(:,:,:)


    ! usually, the routine should not be called with no penalization, but if it still
    ! happens, do nothing.
    if ( params_acm%penalization .eqv. .false.) then
        mask = 0.0_rk
        return
    endif

    ! check if the array has the right dimension, if not, put money in swear jar.
    if (size(mask,1) /= Bs(1)+2*g .or. size(mask,2) /= Bs(2)+2*g .or. size(mask,3) /= Bs(3)+2*g ) then
        write(*,*) shape(mask)
        call abort(777107, "mask: wrong array size, there's pirates, captain!")
    endif

    if (size(mask,4) < 5 ) then
        write(*,*) shape(mask)
        call abort(777108, "mask: wrong number of components (5), there's pirates, captain!")
    endif


    if (.not. allocated(mask_color)) allocate(mask_color(1:Bs(1)+2*g, 1:Bs(2)+2*g, 1:Bs(3)+2*g))

    if (.not. params_acm%initialized) write(*,*) "WARNING: create_mask_3D_ACM called but ACM not initialized"

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! mask function and boundary values
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    select case (params_acm%geometry)

    case ('sphere-fixed')
        if (stage == "time-independent-part" .or. stage == "all-parts") then
            call draw_fixed_sphere(x0, dx, Bs, g, mask )
        endif

    case ('sphere-free')
        if (stage == "time-dependent-part" .or. stage == "all-parts") then
            call draw_free_sphere(x0, dx, Bs, g, mask )
        endif

    case ('fractal_tree')
        !-----------------------------------------------------------------------
        ! FRACTAL TREE
        !-----------------------------------------------------------------------
        if (stage == "time-independent-part" .or. stage == "all-parts") then
            ! fractal trees are time-independent (they have no time-dependent part)
            call Draw_fractal_tree(Insect, x0-dble(g)*dx, dx, mask(:,:,:,1), mask_color, mask(:,:,:,2:4))

            ! store the mask color array as double
            mask(:,:,:,5) = real(mask_color, kind=rk)
        endif


    case ('active_grid')
        !-----------------------------------------------------------------------
        ! ACTIVE GRID
        !-----------------------------------------------------------------------
        if (stage == "time-dependent-part" .or. stage == "all-parts") then
            ! active grids are always time-dependent
            call draw_active_grid_winglets(time, Insect, x0-dble(g)*dx, dx, &
            mask(:,:,:,1), mask_color, mask(:,:,:,2:4))

            ! store the mask color array as double
            mask(:,:,:,5) = real(mask_color, kind=rk)
        endif


    case ('Insect')
        !-----------------------------------------------------------------------
        ! INSECT MODULE
        !-----------------------------------------------------------------------
        ! the insects require us to determine their state vector before they can be drawn
        ! as this is to do only once, not for all blocks
        ! 18 Feb 2021: deactivated the call here because it is (more efficiently) done in module_mask.f90
        ! This is important as for FSI problems the mask function at TIME may have to be recomputed even if we
        ! already computed it at this time!! Think of RK substeps, where several RHS evaluations are to be done
        ! at the same time level but with different input data. Hence, the check below is NOT sufficient in those
        ! cases.
        ! if (abs(time-Insect%time) >= 1.0e-13_rk) then
        !     call Update_Insect(time, Insect)
        ! endif

        select case(stage)
        case ("time-independent-part")
            ! insect body: note non-tethered-flight is a problem
            if (Insect%body_moves == "no") then
                call draw_insect_body( time, x0-dble(g)*dx, dx, mask(:,:,:,1), &
                mask_color, mask(:,:,:,2:4), Insect, delete=.true.)

                ! store the mask color array as double
                mask(:,:,:,5) = real(mask_color, kind=rk)
            endif

            ! we can also simulate an insect together with a fractal tree as turbulence
            ! generators. This part is time-independent as the tree does not move.
            if (Insect%fractal_tree) then
                call draw_fractal_tree(Insect, x0-dble(g)*dx, dx, mask(:,:,:,1), mask_color, mask(:,:,:,2:4))
                ! store the mask color array as double
                mask(:,:,:,5) = real(mask_color, kind=rk)
            endif

        case ("time-dependent-part")
            if (Insect%body_moves == "no") then
                ! wings
                call draw_insect_wings( time, x0-dble(g)*dx, dx, mask(:,:,:,1), &
                mask_color, mask(:,:,:,2:4), Insect, delete=.true.)
            else
                ! draw entire insect. Note: insect module is ghost-nodes aware, but requires origin shift.
                call Draw_Insect( time, Insect, x0-dble(g)*dx, dx, mask(:,:,:,1), mask_color, mask(:,:,:,2:4) )
            endif

            ! store the mask color array as double
            mask(:,:,:,5) = real(mask_color, kind=rk)

        case ("all-parts")
            ! wings and body
            ! draw entire insect. Note: insect module is ghost-nodes aware, but requires origin shift.
            call Draw_Insect( time, Insect, x0-dble(g)*dx, dx, mask(:,:,:,1), mask_color, mask(:,:,:,2:4) )

            ! we can also simulate an insect together with a fractal tree as turbulence
            ! generators. This part is time-independent as the tree does not move.
            if (Insect%fractal_tree) then
                call draw_fractal_tree(Insect, x0-dble(g)*dx, dx, mask(:,:,:,1), mask_color, mask(:,:,:,2:4))
            endif

            ! store the mask color array as double
            mask(:,:,:,5) = real(mask_color, kind=rk)

        case default
            call abort(16072019, "unknown request to create_mask")

        end select


    case ('none')
        mask = 0.0_rk

    case default
        call abort(120001,"ERROR: geometry for 3d VPM is unknown "//params_acm%geometry)

    end select

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! sponge
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Though time-independent, we treat the sponge mask as if it were time-dependent
    ! because it does not have to be refined to the finest level (unlike the mask
    ! function.)
    if (stage == "time-dependent-part" .or. stage == "all-parts") then
        if (params_acm%use_sponge) then
            call sponge_3D( mask(:,:,:,6), x0, dx, Bs, g)
        endif
    endif

end subroutine create_mask_3D_ACM

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! mask data. we can use different trees (4est module) to generate time-dependent/indenpedent
! mask functions separately. This makes the mask routines tree-level routines (and no longer
! block level) so the physics modules have to provide an interface to create the mask at a tree
! level. All parts of the mask shall be included: chi, boundary values, sponges.
! This is a block-level wrapper to fill the mask.
subroutine create_mask_2D_ACM( time, x0, dx, Bs, g, mask, stage )
    implicit none

    ! grid
    integer(kind=ik), intent(in) :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    !> mask term for every grid point of this block
    real(kind=rk), dimension(:,:,:), intent(inout) :: mask
    !> spacing and origin of block
    real(kind=rk), intent(in) :: x0(1:2), dx(1:2), time
    ! sometimes one has to do preparatory work for the mask function => staging idea.
    character(len=*), intent(in) :: stage

    ! some cheap checks
    if (size(mask,1) /= Bs(1)+2*g .or. size(mask,2) /= Bs(2)+2*g ) then
        call abort(634840, "mask: wrong array size, there's pirates, captain!")
    endif

    mask = 0.0_rk

    ! usually, the routine should not be called with no penalization, but if it still
    ! happens, do nothing.
    if (.not. params_acm%penalization) return

    if (.not. params_acm%initialized) write(*,*) "WARNING: create_mask_2D_ACM called but ACM not initialized"

    !---------------------------------------------------------------------------
    ! Mask function and forcing values
    !---------------------------------------------------------------------------
    select case (params_acm%geometry)
    case ('rotating-rod')
        if (stage == "time-dependent-part" .or. stage == "all-parts") then
            call draw_rotating_rod( time, mask, x0, dx, Bs, g )
        endif

    case ('cylinder')
        if (stage == "time-independent-part" .or. stage == "all-parts") then
            call draw_cylinder( mask, x0, dx, Bs, g )
        endif

    case ('cylinder-free')
        if (stage == "time-dependent-part" .or. stage == "all-parts") then
            call draw_free_cylinder( mask, x0, dx, Bs, g )
        endif

    case ('rotating_cylinder')
        if (stage == "time-dependent-part" .or. stage == "all-parts") then
            call draw_rotating_cylinder( time, mask, x0, dx, Bs, g )
        endif

    case ('two-cylinders')
        if (stage == "time-independent-part" .or. stage == "all-parts") then
            call draw_two_cylinders( mask(:,:,1), x0, dx, Bs, g )
        endif

    case ('two-moving-cylinders')
        if (stage == "time-dependent-part" .or. stage == "all-parts") then
            call draw_two_moving_cylinders( time, mask(:,:,:), x0, dx, Bs, g )
        endif

    case ('flapping-wings')
        if (stage == "time-dependent-part" .or. stage == "all-parts") then
            call draw_2d_flapping_wings( time, mask(:,:,:), x0, dx, Bs, g )
        endif

    case ('cavity')
        if (stage == "time-independent-part" .or. stage == "all-parts") then
            call draw_cavity( mask, x0, dx, Bs, g )
        endif

    case ('2D-wingsection')
        if (stage == "time-dependent-part" .or. stage == "all-parts") then
            call draw_2d_wingsections( time, mask, x0, dx, Bs, g )
        endif

    case ('none')
        mask = 0.0_rk

    case default
        call abort(120002,"ERROR: geometry for 2d VPM is unknown"//params_acm%geometry)

    end select

    !---------------------------------------------------------------------------
    ! sponge
    !---------------------------------------------------------------------------
    ! Though time-independent, we treat the sponge mask as if it were time-dependent
    ! because it does not have to be refined to the finest level (unlike the mask
    ! function.)
    if (stage == "time-dependent-part" .or. stage == "all-parts") then
        if (params_acm%use_sponge) then
            call sponge_2D( mask(:,:,6), x0, dx, Bs, g)
        endif
    endif

end subroutine create_mask_2D_ACM

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

subroutine draw_cylinder(mask, x0, dx, Bs, g )

    use module_params
    use module_precision

    implicit none

    ! grid
    integer(kind=ik), intent(in) :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    !> mask term for every grid point of this block
    real(kind=rk), dimension(:,:,:), intent(out)     :: mask
    !> spacing and origin of block
    real(kind=rk), dimension(2), intent(in) :: x0, dx

    ! auxiliary variables
    real(kind=rk)  :: x, y, r, h, dx_min, tmp
    ! loop variables
    integer(kind=ik) :: ix, iy

    if (size(mask,1) /= Bs(1)+2*g .or. size(mask,2) /= Bs(2)+2*g ) then
        call abort(777107, "mask: wrong array size, there's pirates, captain!")
    endif

    ! reset mask array
    mask = 0.0_rk

    ! parameter for smoothing function (width)
    dx_min = 2.0_rk**(-params_acm%Jmax) * params_acm%domain_size(1) / real(params_acm%Bs(1)-1, kind=rk)
    h = 1.5_rk * dx_min

    ! Note: this basic mask function is set on the ghost nodes as well.
    do iy = 1, Bs(2)+2*g
        y = dble(iy-(g+1)) * dx(2) + x0(2) - params_acm%x_cntr(2)
        do ix = 1, Bs(1)+2*g
            x = dble(ix-(g+1)) * dx(1) + x0(1) - params_acm%x_cntr(1)
            ! distance from center of cylinder
            r = dsqrt(x*x + y*y)

            tmp = smoothstep(r, params_acm%R_cyl, h)
            if (tmp >= mask(ix,iy,1)) then
                ! mask function
                mask(ix,iy,1) = tmp
                ! color
                mask(ix,iy,5) = 1.0_rk
            endif
        end do
    end do

end subroutine draw_cylinder

!-------------------------------------------------------------------------------

subroutine draw_free_cylinder(mask, x0, dx, Bs, g )

    use module_params
    use module_precision

    implicit none

    ! grid
    integer(kind=ik), intent(in) :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    !> mask term for every grid point of this block
    real(kind=rk), dimension(:,:,:), intent(out)     :: mask
    !> spacing and origin of block
    real(kind=rk), dimension(2), intent(in) :: x0, dx

    ! auxiliary variables
    real(kind=rk)  :: x, y, r, h, dx_min, tmp
    ! loop variables
    integer(kind=ik) :: ix, iy

    if (size(mask,1) /= Bs(1)+2*g .or. size(mask,2) /= Bs(2)+2*g ) then
        call abort(777107, "mask: wrong array size, there's pirates, captain!")
    endif

    params_acm%u_vert = Insect%STATE(5)
    params_acm%z_vert = Insect%STATE(2)


    ! reset mask array
    mask = 0.0_rk

    ! parameter for smoothing function (width)
    h = Insect%C_smooth*minval(dx(1:2))

    ! Note: this basic mask function is set on the ghost nodes as well.
    do iy = 1, Bs(2)+2*g
        y = dble(iy-(g+1)) * dx(2) + x0(2) - params_acm%z_vert
        do ix = 1, Bs(1)+2*g
            x = dble(ix-(g+1)) * dx(1) + x0(1) - 0.5*params_acm%domain_size(1)
            ! distance from center of cylinder
            r = dsqrt(x*x + y*y)

            tmp = smoothstep(r, params_acm%R_cyl, h)
            if (tmp >= mask(ix,iy,1)) then
                ! mask function
                mask(ix,iy,1) = tmp
                ! vertical () velocity
                mask(ix,iy,3) = params_acm%u_vert
                ! color
                mask(ix,iy,5) = 1.0_rk
            endif
        end do
    end do

end subroutine draw_free_cylinder

!-------------------------------------------------------------------------------

subroutine draw_free_sphere(x0, dx, Bs, g, mask )

    use module_params
    use module_precision

    implicit none

    ! grid
    integer(kind=ik), intent(in) :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    !> mask term for every grid point of this block
    real(kind=rk), dimension(:,:,:,:), intent(out) :: mask
    !> spacing and origin of block
    real(kind=rk), dimension(1:3), intent(in) :: x0, dx

    ! auxiliary variables
    real(kind=rk)  :: x, y, z, r, h, dx_min, tmp
    ! loop variables
    integer(kind=ik) :: ix, iy, iz

    if (size(mask,1) /= Bs(1)+2*g .or. size(mask,2) /= Bs(2)+2*g ) then
        call abort(777107, "mask: wrong array size, there's pirates, captain!")
    endif

    ! reset mask array
    mask = 0.0_rk

    ! parameter for smoothing function (width)
    h = Insect%C_smooth*minval(dx)

    ! Note: this basic mask function is set on the ghost nodes as well.
    do iz = g+1, Bs(3)+g+ONE_SKIPREDUNDANT
        z = dble(iz-(g+1)) * dx(3) + x0(3) - Insect%STATE(3)
        do iy = g+1, Bs(2)+g+ONE_SKIPREDUNDANT
            y = dble(iy-(g+1)) * dx(2) + x0(2) - Insect%STATE(2)
            do ix = g+1, Bs(1)+g+ONE_SKIPREDUNDANT
                x = dble(ix-(g+1)) * dx(1) + x0(1) - Insect%STATE(1)

                ! distance from center of cylinder
                r = dsqrt(x*x + y*y + z*z)

                mask(ix,iy,iz,1) = smoothstep(r, params_acm%R_cyl, h)
                mask(ix,iy,iz,2) = Insect%STATE(4)
                mask(ix,iy,iz,3) = Insect%STATE(5)
                mask(ix,iy,iz,4) = Insect%STATE(6)
                ! color
                mask(ix,iy,iz,5) = 1.0_rk
            end do
        end do
    end do

end subroutine draw_free_sphere

!-------------------------------------------------------------------------------

subroutine draw_fixed_sphere(x0, dx, Bs, g, mask )

    use module_params
    use module_precision

    implicit none

    ! grid
    integer(kind=ik), intent(in) :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    !> mask term for every grid point of this block
    real(kind=rk), dimension(:,:,:,:), intent(out) :: mask
    !> spacing and origin of block
    real(kind=rk), dimension(1:3), intent(in) :: x0, dx

    ! auxiliary variables
    real(kind=rk)  :: x, y, z, r, h, tmp, dx_min
    ! loop variables
    integer(kind=ik) :: ix, iy, iz

    if (size(mask,1) /= Bs(1)+2*g .or. size(mask,2) /= Bs(2)+2*g ) then
        call abort(777107, "mask: wrong array size, there's pirates, captain!")
    endif

    ! reset mask array
    mask = 0.0_rk

    ! parameter for smoothing function (width)

    !h = 2*minval(dx)
    dx_min = 2.0_rk**(-params_acm%Jmax) * params_acm%domain_size(1) / real(params_acm%Bs(1)-1, kind=rk)
    h = 1.5_rk * dx_min


    ! Note: this basic mask function is set on the ghost nodes as well.
    do iz = g+1, Bs(3)+g+ONE_SKIPREDUNDANT
        z = dble(iz-(g+1)) * dx(3) + x0(3) - params_acm%x_cntr(3)
        do iy = g+1, Bs(2)+g+ONE_SKIPREDUNDANT
            y = dble(iy-(g+1)) * dx(2) + x0(2) - params_acm%x_cntr(2)
            do ix = g+1, Bs(1)+g+ONE_SKIPREDUNDANT
                x = dble(ix-(g+1)) * dx(1) + x0(1) - params_acm%x_cntr(1)

                ! distance from center of cylinder
                r = dsqrt(x*x + y*y + z*z)

                tmp = smoothstep(r, params_acm%R_cyl, h)

                if (tmp >= mask(ix,iy,iz,1)) then
                    mask(ix,iy,iz,1) = tmp
                    ! color
                    mask(ix,iy,iz,5) = 1.0_rk
                endif
            end do
        end do
    end do

end subroutine draw_fixed_sphere

!-------------------------------------------------------------------------------

subroutine draw_cylinderz(x0, dx, Bs, g, mask )

    use module_params
    use module_precision

    implicit none

    ! grid
    integer(kind=ik), intent(in) :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    !> mask term for every grid point of this block
    real(kind=rk), dimension(:,:,:,:), intent(out) :: mask
    !> spacing and origin of block
    real(kind=rk), dimension(1:3), intent(in) :: x0, dx

    ! auxiliary variables
    real(kind=rk)  :: x, y, z, r, h, dx_min, tmp
    ! loop variables
    integer(kind=ik) :: ix, iy, iz

    ! reset mask array
    mask = 0.0_rk

    ! parameter for smoothing function (width)
    h = Insect%C_smooth*minval(dx)

    do iz = g+1, Bs(3)+g+ONE_SKIPREDUNDANT
        z = dble(iz-(g+1)) * dx(3) + x0(3) - params_acm%domain_size(3)/2.0
        do iy = g+1, Bs(2)+g+ONE_SKIPREDUNDANT
            y = dble(iy-(g+1)) * dx(2) + x0(2) - params_acm%domain_size(2)/2.0
            do ix = g+1, Bs(1)+g+ONE_SKIPREDUNDANT
                x = dble(ix-(g+1)) * dx(1) + x0(1) - params_acm%domain_size(1)/2.0

                ! distance from center of cylinder
                r = dsqrt(x*x + y*y)

                mask(ix,iy,iz,1) = smoothstep(r, params_acm%R_cyl, h)
                ! color
                mask(ix,iy,iz,5) = 1.0_rk
            end do
        end do
    end do

end subroutine draw_cylinderz

!-------------------------------------------------------------------------------

subroutine draw_rotating_cylinder(time, mask, x0, dx, Bs, g )

    use module_params
    use module_precision

    implicit none

    real(kind=rk) :: time
    integer(kind=ik), intent(in) :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    !> mask term for every grid point of this block
    real(kind=rk), dimension(:,:,:), intent(out)     :: mask
    !> spacing and origin of block
    real(kind=rk), dimension(2), intent(in) :: x0, dx

    ! auxiliary variables
    real(kind=rk) :: x, y, r, h, dx_min, tmp, x00, y00, radius, frequ, alpha
    ! loop variables
    integer(kind=ik) :: ix, iy

    if (size(mask,1) /= Bs(1)+2*g .or. size(mask,2) /= Bs(2)+2*g ) then
        call abort(777107, "mask: wrong array size, there's pirates, captain!")
    endif

    ! frequency of rotation is unity:
    frequ = 1.0_rk
    ! radius of rotation:
    radius = 1.0_rk
    alpha = 2.0_rk * pi * time * frequ
    ! cylinder mid-point as a function of time:
    x00 = params_acm%x_cntr(1) + dcos(alpha)*radius
    y00 = params_acm%x_cntr(2) + dsin(alpha)*radius


    ! reset mask array
    mask = 0.0_rk

    ! parameter for smoothing function (width)
    dx_min = 2.0_rk**(-params_acm%Jmax) * params_acm%domain_size(1) / real(params_acm%Bs(1)-1, kind=rk)
    h = 1.5_rk * dx_min

    ! Note: this basic mask function is set on the ghost nodes as well.
    do iy = 1, Bs(2)+2*g
        y = dble(iy-(g+1)) * dx(2) + x0(2)
        do ix = 1, Bs(1)+2*g
            x = dble(ix-(g+1)) * dx(1) + x0(1)
            ! distance from center of cylinder
            r = dsqrt( (x-x00)*(x-x00) + (y-y00)*(y-y00) )

            tmp = smoothstep(r, params_acm%R_cyl, h)

            if (tmp >= mask(ix,iy,1) .and. tmp > 0.0_rk) then
                ! mask function
                mask(ix,iy,1) = tmp
                ! usx
                mask(ix,iy,2) = -2.0_rk * pi * frequ * sin(alpha)
                ! usy
                mask(ix,iy,3) = +2.0_rk * pi * frequ * cos(alpha)
                ! color
                mask(ix,iy,5) = 1.0_rk
            endif
        end do
    end do

end subroutine draw_rotating_cylinder

! ------------------------------------------------------------------------------

subroutine draw_cavity(mask, x0, dx, Bs, g )

    use module_params
    use module_precision

    implicit none

    ! grid
    integer(kind=ik), intent(in) :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    !> mask term for every grid point of this block
    real(kind=rk), dimension(:,:,:), intent(out)     :: mask
    !> spacing and origin of block
    real(kind=rk), dimension(2), intent(in) :: x0, dx

    ! auxiliary variables
    real(kind=rk)  :: x, y, length, tmp
    real(kind=rk), parameter :: c_smooth = 2.0_rk
    ! loop variables
    integer(kind=ik) :: ix, iy

    if (size(mask,1) /= Bs(1)+2*g .or. size(mask,2) /= Bs(2)+2*g ) then
        call abort(777107, "mask: wrong array size, there's pirates, captain!")
    endif

    ! reset mask array
    mask = 0.0_rk
    length= params_acm%length

    if (params_acm%dx_min <= 0.0_rk) call abort(0509198,"params_acm%dx_min invalid (not set?)")


    ! Note: this basic mask function is set on the ghost nodes as well.
    ! Discontinuous version
    ! do iy = 1, Bs(2)+2*g
    !     y = dble(iy-(g+1)) * dx(2) + x0(2)
    !     do ix = 1, Bs(1)+2*g
    !         x = dble(ix-(g+1)) * dx(1) + x0(1)
    !
    !         if ( x<length .or. x>params_acm%domain_size(1)-length &
    !             .or. y<length .or. y>params_acm%domain_size(2)-length) then
    !             ! mask function
    !             mask(ix,iy,1) = 1.0_rk
    !             ! color
    !             mask(ix,iy,5) = 1.0_rk
    !         endif
    !     end do
    ! end do

    ! smoothed version
    do iy = 1, Bs(2)+2*g
        y = dble(iy-(g+1)) * dx(2) + x0(2)
        do ix = 1, Bs(1)+2*g
            x = dble(ix-(g+1)) * dx(1) + x0(1)

            ! distance to borders of domain
            tmp = minval( (/x,y,-(x-params_acm%domain_size(1)),-(y-params_acm%domain_size(2))/) )

            mask(ix,iy,1) = smoothstep( tmp, length, c_smooth*params_acm%dx_min)

        end do
    end do

end subroutine draw_cavity

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

subroutine draw_two_cylinders( mask, x0, dx, Bs, g)

    use module_params
    use module_precision

    implicit none

    ! grid
    integer(kind=ik), intent(in) :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    !> mask term for every grid point of this block
    real(kind=rk), dimension(:,:), intent(out)     :: mask
    !> spacing and origin of block
    real(kind=rk), dimension(2), intent(in)        :: x0, dx

    ! auxiliary variables
    real(kind=rk)         :: x1, x2, y1, y2, R, cx1, cx2, cy1,&
    cy2, r_1, r_2, h, mask1, mask2
    ! loop variables
    integer(kind=ik)      :: ix, iy

    !---------------------------------------------------------------------------------------------
    ! variables initialization
    if (size(mask,1) /= Bs(1)+2*g .or. size(mask,2) /= Bs(2)+2*g  ) then
        call abort(777107, "mask: wrong array size, there's pirates, captain!")
    endif

    ! reset mask array
    mask = 0.0_rk
    mask1 = 0.0_rk
    mask2 = 0.0_rk

    !---------------------------------------------------------------------------------------------
    ! main body

    ! center of the first cylinder
    cx1 = 0.5884_rk*params_acm%domain_size(1)
    cy1 = 0.4116_rk*params_acm%domain_size(2)

    ! center of the second cylinder
    cx2 = 0.4116_rk*params_acm%domain_size(1)
    cy2 = 0.5884_rk*params_acm%domain_size(2)

    ! radius of the cylinders
    R = params_acm%R_cyl
    ! parameter for smoothing function (width)
    h = 1.5_rk*max(dx(1), dx(2))

    do iy=1, Bs(2)+2*g
        y1 = dble(iy-(g+1)) * dx(2) + x0(2) - cy1
        y2 = dble(iy-(g+1)) * dx(2) + x0(2) - cy2
        do ix=1, Bs(1)+2*g
            x1 = dble(ix-(g+1)) * dx(1) + x0(1) - cx1
            x2 = dble(ix-(g+1)) * dx(1) + x0(1) - cx2
            ! distance from center of cylinder 1
            r_1 = dsqrt(x1*x1 + y1*y1)
            ! distance from center of cylinder 2
            r_2 = dsqrt(x2*x2 + y2*y2)
            if (params_acm%smooth_mask) then
                mask1 = smoothstep( r_1, R, h)
                mask2 = smoothstep( r_2, R, h)
                mask(ix,iy) = mask1 + mask2
            else
                ! if point is inside one of the cylinders, set mask to 1
                if (r_1 <= R) then
                    mask(ix,iy) = 1.0_rk
                elseif ( r_2 <= R) then
                    mask(ix,iy) = 1.0_rk
                else
                    mask(ix,iy) = 0.0_rk
                end if
            end if
        end do
    end do
end subroutine draw_two_cylinders


subroutine draw_two_moving_cylinders(time, mask, x0, dx, Bs, g)

    use module_params
    use module_precision

    implicit none

    ! grid
    integer(kind=ik), intent(in) :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    !> mask term for every grid point of this block
    real(kind=rk), dimension(:,:,:), intent(out)     :: mask
    !> spacing and origin of block
    real(kind=rk), dimension(2), intent(in)        :: x0, dx
    !> simulation time
    real(kind=rk), intent(in) :: time
    ! auxiliary variables
    real(kind=rk)         :: x1, x2, y1, y2, R1, R2, cx1, cx2, cy1,&
    cy2, r_1, r_2, h, mask1, mask2, freq, vy2
    real(kind=rk), allocatable, save :: mu(:)
    !real(kind=rk)         :: f1, f2, St1, St2 ,Re1, Re2
    ! loop variables
    integer(kind=ik)      :: ix, iy, k, nfft_y0

    !---------------------------------------------------------------------------------------------
    ! variables initialization
    if (size(mask,1) /= Bs(1)+2*g .or. size(mask,2) /= Bs(2)+2*g  ) then
        call abort(777107, "mask: wrong array size, there's pirates, captain!")
    endif
    nfft_y0 = wingsections(1)%nfft_y0
    if (.not. allocated(mu)) allocate(mu(nfft_y0))
    ! reset mask array
    mask = 0.0_rk
    mask1 = 0.0_rk
    mask2 = 0.0_rk
    !---------------------------------------------------------------------------------------------
    ! main body
    ! radius of the cylinders
    R1 = 1.0_rk * params_acm%R_cyl
    R2 = 1.0_rk * params_acm%R_cyl
    ! Here we set the frequency of the 2. cylinder oscillating up and down behind
    ! the 1. cyl. It should be choosen such that the characteristic time scale of
    ! the vortex shedding is larger then the movement!
    ! Reynolds Number:
    ! Re1 = 2 * params_acm%u_mean_set(1) * R1/ params_acm%nu
    ! Re2 = 2 * params_acm%u_mean_set(1) * R2/ params_acm%nu
    ! ! strouhal number: (see: https://en.wikipedia.org/wiki/K%C3%A1rm%C3%A1n_vortex_street)
    ! St1 = 0.198 * (1-19.7/Re1)
    ! St2 = 0.198 * (1-19.7/Re2)
    ! ! vortex shedding frequency
    ! f1 = St1*params_acm%u_mean_set(1)/(2*R1)
    ! f2 = St2*params_acm%u_mean_set(1)/(2*R2)
    ! ! make cylinder movement slow in comparison to vortex shedding frequency:
    ! freq = min(f1,f2) / 10
    freq = 0.2e-2
    mu = wingsections(1)%ai_y0
    ! center of the first cylinder
    cx1 = 0.250_rk * params_acm%domain_size(1)
    cy1 = 0.500_rk * params_acm%domain_size(2)
    ! center of the second cylinder (oscillates behind 1. cylinder)
    cx2 = 0.500_rk * params_acm%domain_size(1)
    cy2 = cy1
    vy2 = 0
    do k = 1, nfft_y0
        cy2 = cy2 + mu(k) * sin(2*k*pi*freq*time)
        ! velocity of moving cylinder
        vy2 = vy2 + mu(k) * 2 * pi * freq * k * cos(2*k*pi*freq*time)
    end do
    ! parameter for smoothing function (width)
    h = 1.5_rk*max(dx(1), dx(2))

    do iy=1, Bs(2)+2*g
        y1 = dble(iy-(g+1)) * dx(2) + x0(2) - cy1
        y2 = dble(iy-(g+1)) * dx(2) + x0(2) - cy2
        do ix=1, Bs(1)+2*g
            x1 = dble(ix-(g+1)) * dx(1) + x0(1) - cx1
            x2 = dble(ix-(g+1)) * dx(1) + x0(1) - cx2
            ! distance from center of cylinder 1
            r_1 = dsqrt(x1*x1 + y1*y1)
            ! distance from center of cylinder 2
            r_2 = dsqrt(x2*x2 + y2*y2)
            if (params_acm%smooth_mask) then
                mask1 = smoothstep( r_1, R1, h)
                mask2 = smoothstep( r_2, R2, h)
                if (mask2>0.0) mask(ix,iy,3) = vy2
                mask(ix,iy,1) = mask1 + mask2
            else
                ! if point is inside one of the cylinders, set mask to 1
                if (r_1 <= R1) then
                    mask(ix,iy,1) = 1.0_rk
                elseif ( r_2 <= R2) then
                    mask(ix,iy,1) = 1.0_rk
                    mask(ix,iy,2) = 0.0_rk
                    mask(ix,iy,3) = vy2
                else
                    mask(ix,iy,:) = 0.0_rk
                end if
            end if
        end do
    end do


end subroutine draw_two_moving_cylinders


subroutine draw_2d_flapping_wings(time, mask, x0, dx, Bs, g)
    ! simple 2D insect
    ! taken from publication:
    !       Fluid Dyn. Res. 44 (2012), Keigo Ota, Kosuke Suzuki
    !       and Takaji Inamuro
    ! Reynolds number is given by:
    !  Re = u_tip*L/nu
    !  time-averaged tip speed u_tip, chord length L, cinematic viscousity nu
    !  u_tip = 4 * L * freq * A   (freq ... frequency, A... max(theta))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !          wing_l     wing_r
    !             \       /
    !              \     /
    !               \   /
    !                \ /
    !                 O        theta is angle between x and wing_r/wing_l
    !                body
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    use module_params
    use module_precision

    implicit none

    ! grid
    integer(kind=ik), intent(in) :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    !> mask term for every grid point of this block
    real(kind=rk), dimension(:,:,:), intent(out)     :: mask
    !> spacing and origin of block
    real(kind=rk), dimension(2), intent(in)        :: x0, dx
    !> simulation time
    real(kind=rk), intent(in) :: time
    ! auxiliary variables
    real(kind=rk)         :: x1, x2, y1, y2, R, cxr, cxl, cyr, cyl, c, K, L, b2,&
                             r_wing_r, r_wing_l,r_body, h, mask_wing_r, x, y, xb, yb, &
                             mask_wing_l, mask_body, freq, x_bodycenter(2), theta, theta_dt, &
                             cos_theta, sin_theta, uwing_x, uwing_y, ubody_x, ubody_y, A
    integer(kind=ik)            :: ix, iy
    integer(kind=2), parameter  :: color_l=2, color_r=3
    integer(kind=2)             :: color
    !---------------------------------------------------------------------------------------------
    ! variables initialization
    if (size(mask,1) /= Bs(1)+2*g .or. size(mask,2) /= Bs(2)+2*g  ) then
        call abort(777107, "mask: wrong array size, there's pirates, captain!")
    endif

    ! reset mask array
    mask = 0.0_rk
    mask_wing_r = 0.0_rk ! right wing
    mask_wing_l = 0.0_rk ! left wing
    mask_body = 0.0_rk
    ! ---------------
    ! body center of cylinder
    x_bodycenter = params_acm%x_cntr(1:2)
    ! size of the body
    R = params_acm%R_cyl
    ! chord length
    c = params_acm%length
    ! the wings are elipsoids, center is half of the chord length + radius of the body
    L = c/2 + R
    !frequency
    freq = params_acm%freq
    ! half axis ration of elipsoid
    b2 = 0.1**2
    ! angle between wing and x axis
    A = pi*0.25
    theta = A * cos(2.0_rk*pi*freq*time)
    ! angular velocity
    theta_dt = -2.0_rk*pi*freq*A*sin(2.0_rk*pi*freq*time)
    ! precompute cosinus and sinus
    cos_theta = cos(theta)
    sin_theta = sin(theta)
    ! velocity of body
    ubody_x = 0.0_rk
    ubody_y = 0.0_rk
    ! centers of wing elipsoids
    ! right wing
    cxr = x_bodycenter(1) + L * cos_theta
    cyr = x_bodycenter(2) + L * sin_theta
    ! left wing
    cxl = x_bodycenter(1) - L * cos_theta
    cyl = x_bodycenter(2) + L * sin_theta
    ! parameter for smoothing function (width)
    h = 1.5_rk*max(dx(1), dx(2))
    do iy=1, Bs(2)+2*g
        y = dble(iy-(g+1)) * dx(2) + x0(2)
        do ix=1, Bs(1)+2*g
            x = dble(ix-(g+1)) * dx(1) + x0(1)
            ! rotate elipsoid of wing right
            x1 =  cos_theta * (x - cxr) + sin_theta * (y - cyr)
            y1 = -sin_theta * (x - cxr) + cos_theta * (y - cyr)
            ! rotate elipsoid of wing left
            x2 =  cos_theta * (x - cxl) - sin_theta * (y - cyl)
            y2 =  sin_theta * (x - cxl) + cos_theta * (y - cyl)
            ! calculate center of body
            xb = x - x_bodycenter(1)
            yb = y - x_bodycenter(2)
            ! distance from center of wing 1
            r_wing_r = dsqrt(x1*x1 + y1*y1/b2)
            ! distance from center of wing 2
            r_wing_l = dsqrt(x2*x2 + y2*y2/b2)
            ! body
            r_body = dsqrt(xb*xb + yb*yb)
            ! compute velocity at the wings
            ! Note that due to symmetry of the movement,
            ! only the x component switches sign from wing_r to wing_l
            uwing_x = -sin_theta * (theta_dt * r_body)
            uwing_y = cos_theta * (theta_dt * r_body)
            ! reset color
            color = 0
            ! draw mask
            if (params_acm%smooth_mask) then
                mask_wing_r = smoothstep( r_wing_r, c*0.5_rk, h)
                mask_wing_l = smoothstep( r_wing_l, c*0.5_rk, h)
                mask_body = smoothstep( r_body, R, h)
                mask(ix,iy,1) = mask_wing_r + mask_wing_l + mask_body
                if (mask_wing_r>0.0_rk) color = color_r
                if (mask_wing_l>0.0_rk) color = color_l
            else
                ! if point is inside one of the cylinders, set mask to 1
                if (r_wing_r <= c*0.5_rk ) then
                    mask(ix,iy,1) = 1.0_rk
                    color = color_r
                elseif ( r_wing_l <= c*0.5_rk ) then
                    mask(ix,iy,1) = 1.0_rk
                    color = color_l
                elseif ( r_body <= R ) then
                    mask(ix,iy,1) = 1.0_rk
                else
                    mask(ix,iy,:) = 0.0_rk
                end if
            end if

            ! set the velocity values inside the rigid body domain
            if (color == color_r) then
                mask(ix,iy,2) = ubody_x + uwing_x
                mask(ix,iy,3) = ubody_y + uwing_y
            elseif (color == color_l) then
                mask(ix,iy,2) = ubody_x - uwing_x ! minus because of oposit wing movement direction
                mask(ix,iy,3) = ubody_y + uwing_y
            else
                mask(ix,iy,2) = ubody_x
                mask(ix,iy,3) = ubody_y
            end if
        end do
    end do


end subroutine draw_2d_flapping_wings

subroutine draw_rotating_rod(time, mask, x0, dx, Bs, g)
    use module_params
    use module_precision

    implicit none

    ! grid
    integer(kind=ik), intent(in) :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    !> mask term for every grid point of this block
    real(kind=rk), dimension(:,:,:), intent(out)     :: mask
    !> spacing and origin of block
    real(kind=rk), dimension(2), intent(in)        :: x0, dx
    !> simulation time
    real(kind=rk), intent(in) :: time

    integer :: ix, iy, iz, mpicode
    real (kind=rk) :: x2, y2, vx2, vy2, vx2t, vy2t, anglez2, omz2, omz2t, x00, y00
    real (kind=rk) :: x, y, xref, yref, xlev, ylev, tmp, N, rref, rmax, hsmth, Am, alpham
    real (kind=rk) :: Af, Sxf, Syf, Jf, forcex, forcey, torquez

    ! reset everything
    mask = 0.0_rk

    N = Insect%C_smooth ! smoothing coefficient
    hsmth = N*minval(dx) ! smoothing layer thickness
    rmax = 0.5d0

    x00 = params_acm%domain_size(1)/2.0_rk
    y00 = params_acm%domain_size(2)/2.0_rk

    ! Flapping parameters
    Am = 1.00d0
    alpham = 0.25d0*pi

    ! Update kinematics
    x2 = x00 + Am * dcos(time/Am)
    y2 = y00
    anglez2 = 0.5d0*pi + alpham * dsin(time/Am)
    vx2 = - dsin(time/Am)
    vy2 = 0.0d0
    omz2 = alpham/Am * dcos(time/Am)
    vx2t = - 1.0/Am * cos(time/Am)
    vy2t = 0.0d0
    omz2t = - alpham/Am**2 * sin(time/Am)

    ! For all grid points of this subdomain
    do iy = g+1, Bs(2)+g+ONE_SKIPREDUNDANT
        y = dble(iy-(g+1)) * dx(2) + x0(2) - y2

        do ix = g+1, Bs(1)+g+ONE_SKIPREDUNDANT
            x = dble(ix-(g+1)) * dx(1) + x0(1) - x2

            xref = x*dcos(anglez2) + y*dsin(anglez2)
            yref = y*dcos(anglez2) - x*dsin(anglez2)
            rref = dsqrt( xref**2 + 4.0d0**2 * yref**2 ) ! Radius in cylindrical coordinates

            tmp = smoothstep(rref, rmax-0.0d0*hsmth, hsmth)
            ! call SmoothStep (tmp, rref, rmax-0.0d0*hsmth, hsmth)


            mask(ix,iy,1) = tmp
            mask(ix,iy,2) = -omz2*y + vx2
            mask(ix,iy,3) = +omz2*x + vy2
            mask(ix,iy,4) = 0.0_rk
            mask(ix,iy,5) = 1.0_rk

        enddo
    enddo


end subroutine draw_rotating_rod
