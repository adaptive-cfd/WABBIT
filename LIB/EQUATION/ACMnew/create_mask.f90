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


    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! mask function and boundary values
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    select case (params_acm%geometry)
    case ('fractal_tree')
        !-----------------------------------------------------------------------
        ! FRACTAL TREE
        !-----------------------------------------------------------------------
        if (stage == "time-independent-part" .or. stage == "all-parts") then
            ! fractal trees are time-independent (they have no time-dependent part)
            call Draw_fractal_tree(Insect, x0-dble(g)*dx, dx, mask(:,:,:,1), mask_color, mask(:,:,:,2:4))
        endif

        ! store the mask color array as double
        mask(:,:,:,5) = real(mask_color, kind=rk)

    case ('Insect')
        !-----------------------------------------------------------------------
        ! INSECT MODULE
        !-----------------------------------------------------------------------
        ! the insects require us to determine their state vector before they can be drawn
        ! as this is to do only once, not for all blocks
        if ( abs(time-Insect%time) >= 1.0e-13_rk) then
            call Update_Insect(time, Insect)
        endif

        select case(stage)
        case ("time-independent-part")
            ! insect body: note non-tethered-flight is a problem
            if (Insect%body_moves == "no") then
                call draw_insect_body( time, x0-dble(g)*dx, dx, mask(:,:,:,1), &
                mask_color, mask(:,:,:,2:4), Insect, delete=.true.)

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

    !---------------------------------------------------------------------------
    ! Mask function and forcing values
    !---------------------------------------------------------------------------
    select case (params_acm%geometry)
    case ('cylinder')
        if (stage == "time-independent-part" .or. stage == "all-parts") then
            call draw_cylinder( mask(:,:,1), x0, dx, Bs, g )
        endif

    case ('two-cylinders')
        if (stage == "time-independent-part" .or. stage == "all-parts") then
            call draw_two_cylinders( mask(:,:,1), x0, dx, Bs, g )
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
    real(kind=rk), dimension(:,:), intent(out)     :: mask
    !> spacing and origin of block
    real(kind=rk), dimension(2), intent(in)                   :: x0, dx

    ! auxiliary variables
    real(kind=rk)                                             :: x, y, r, h, dx_min
    ! loop variables
    integer(kind=ik)                                          :: ix, iy

    if (size(mask,1) /= Bs(1)+2*g .or. size(mask,2) /= Bs(2)+2*g ) then
        call abort(777107, "mask: wrong array size, there's pirates, captain!")
    endif

    ! reset mask array
    mask = 0.0_rk

    ! parameter for smoothing function (width)
    dx_min = 2.0_rk**(-params_acm%Jmax) * params_acm%domain_size(1) / real(params_acm%Bs(1)-1, kind=rk)
    h = 1.5_rk * dx_min

    do iy=1, Bs(2)+2*g
        y = dble(iy-(g+1)) * dx(2) + x0(2) - params_acm%x_cntr(2)
        do ix=1, Bs(1)+2*g
            x = dble(ix-(g+1)) * dx(1) + x0(1) - params_acm%x_cntr(1)
            ! distance from center of cylinder
            r = dsqrt(x*x + y*y)
            if (params_acm%smooth_mask) then
                mask(ix,iy) = smoothstep(r, params_acm%R_cyl, h)
            else
                ! if point is inside the cylinder, set mask to 1
                if (r <= params_acm%R_cyl) then
                    mask(ix,iy) = 1.0_rk
                else
                    mask(ix,iy) = 0.0_rk
                end if
            end if
        end do
    end do

end subroutine draw_cylinder

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
