!! input:    - params, origin and spacing of the block, grid parameters \n
!! output:   - sponge term \n
!*********************************************************************************************
subroutine sponge_2D(sponge, x0, dx, Bs, g)
    implicit none

    integer(kind=ik), intent(in)                   :: g                         ! grid
    integer(kind=ik), dimension(3), intent(in)     :: Bs
    real(kind=rk), dimension(:,:), intent(out)     :: sponge                    !> sponge term for every grid point of this block
    real(kind=rk), dimension(2), intent(in)        :: x0, dx                    !> spacing and origin of block
    real(kind=rk)                                  :: x, y, tmp, p, radius      ! auxiliary variables
    real(kind=rk), dimension(2)                    :: global_min, global_max, center ! bounds of the (possibly cropped) domain
    integer(kind=ik)                               :: ix, iy                    ! loop variables

    if (.not. params_acm%initialized) write(*,*) "WARNING: sponge_2D called but ACM not initialized"

    ! cropped domains can differ from 0 and domain_size, so we compute the actual global min/max here
    global_min = params_acm%domain_cropping_min(1:2) * params_acm%domain_size(1:2)
    global_max = params_acm%domain_cropping_max(1:2) * params_acm%domain_size(1:2)

    if (params_acm%sponge_type == "rect") then
        ! rectangular sponge with 45deg edges
        do iy = g+1, Bs(2)+g
            y = dble(iy-(g+1)) * dx(2) + x0(2)

            do ix = g+1, Bs(1)+g
                x = dble(ix-(g+1)) * dx(1) + x0(1)

                ! distance to borders of domain
                tmp = minval( (/x-global_min(1), y-global_min(2), -(x-global_max(1)), -(y-global_max(2))/) )

                sponge(ix,iy) = step_cosine( tmp, 0.5_rk*params_acm%L_sponge, 0.5_rk*params_acm%L_sponge)
            end do
        end do

    elseif (params_acm%sponge_type == "p-norm") then
        ! p-norm sponge. The shape of the sponge is dictated as the p-norm
        ! https://de.wikipedia.org/wiki/P-Norm
        ! which is a nice and simple way to get a rectangle with round corners.

        if ( maxval(abs((global_max-global_min) - (global_max(1)-global_min(1)))) > 1.0e-10_rk) then
            call abort(1610184,"ERROR: for the p-norm sponge, the (cropped) domain has to be same size in all directions.")
        endif

        p = params_acm%p_sponge
        center = 0.5_rk * (global_min + global_max)
        radius = 0.5_rk * (global_max(1) - global_min(1))

        do iy = g+1, Bs(2)+g
            y = dble(iy-(g+1)) * dx(2) + x0(2) - center(2)
            do ix = g+1, Bs(1)+g
                x = dble(ix-(g+1)) * dx(1) + x0(1) - center(1)

                ! distance to borders of domain
                tmp = -( (x**p + y**p)**(1.0_rk/p) - radius)
                sponge(ix,iy) = step_cosine( tmp, 0.5_rk*params_acm%L_sponge, &
                0.5_rk*params_acm%L_sponge)
            end do
        end do
    else
        call abort(1610180,"Sponge-type is unknown")
    endif

end subroutine sponge_2D


subroutine sponge_3D(sponge, x0, dx, Bs, g)
    implicit none

    ! grid
    integer(kind=ik), intent(in)  :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    !> sponge term for every grid point of this block
    real(kind=rk), dimension(:,:,:), intent(out)     :: sponge
    !> spacing and origin of block
    real(kind=rk), intent(in) :: x0(1:3), dx(1:3)

    ! auxiliary variables
    real(kind=rk)     :: x, y, z, tmp, p, radius, pinv
    real(kind=rk), dimension(3) :: global_min, global_max, center ! bounds of the (possibly cropped) domain
    ! loop variables
    integer(kind=ik)  :: ix, iy, iz

    if (.not. params_acm%initialized) write(*,*) "WARNING: sponge_3D called but ACM not initialized"

    ! cropped domains can differ from 0 and domain_size, so we compute the actual global min/max here
    global_min = params_acm%domain_cropping_min * params_acm%domain_size
    global_max = params_acm%domain_cropping_max * params_acm%domain_size

    if (params_acm%sponge_type == "rect") then

        ! check if this block is in sponge layer
        ! find the point in the block closest to the domain boundary
        tmp = min(minval(x0 - global_min), minval(global_max - (x0 + dx*dble(Bs))))
        if (tmp > params_acm%L_sponge) then
            sponge = 0.0_rk
        else
            ! rectangular sponge with 45deg edges
            do iz = g+1, Bs(3)+g
                z = dble(iz-(g+1)) * dx(3) + x0(3)
                do iy = g+1, Bs(2)+g
                    y = dble(iy-(g+1)) * dx(2) + x0(2)
                    do ix = g+1, Bs(1)+g
                        x = dble(ix-(g+1)) * dx(1) + x0(1)

                        ! distance to borders of domain
                        tmp = minval( (/x-global_min(1),y-global_min(2),z-global_min(3),-(x-global_max(1)),&
                            -(y-global_max(2)),-(z-global_max(3))/) )

                        sponge(ix,iy,iz) = step_cosine( tmp, 0.5_rk*params_acm%L_sponge, 0.5_rk*params_acm%L_sponge)
                    end do
                end do
            end do
        endif

        ! sponge for using with symmetry_BC
        ! insect is supposed to be at y=0
    elseif (params_acm%sponge_type == "rect-symmetry-y") then
        ! rectangular sponge with 45deg edges
        do iz = g+1, Bs(3)+g
            z = dble(iz-(g+1)) * dx(3) + x0(3)
            do iy = g+1, Bs(2)+g
                y = dble(iy-(g+1)) * dx(2) + x0(2)
                do ix = g+1, Bs(1)+g
                    x = dble(ix-(g+1)) * dx(1) + x0(1)

                    ! distance to borders of domain
                    tmp = minval( (/x-global_min(1),z-global_min(3),-(x-global_max(1)),-(y-global_max(2)),-(z-global_max(3))/) )

                    sponge(ix,iy,iz) = step_cosine( tmp, 0.5_rk*params_acm%L_sponge, 0.5_rk*params_acm%L_sponge)
                end do
            end do
        end do

        ! sponge for using with symmetry_BC
        ! insect is supposed to be at y=0 z=0
    elseif (params_acm%sponge_type == "rect-symmetry-yz") then
        ! rectangular sponge with 45deg edges
        do iz = g+1, Bs(3)+g
            z = dble(iz-(g+1)) * dx(3) + x0(3)
            do iy = g+1, Bs(2)+g
                y = dble(iy-(g+1)) * dx(2) + x0(2)
                do ix = g+1, Bs(1)+g
                    x = dble(ix-(g+1)) * dx(1) + x0(1)

                    ! distance to borders of domain
                    tmp = minval( (/x-global_min(1),-(x-global_max(1)),-(y-global_max(2)),-(z-global_max(3))/) )

                    sponge(ix,iy,iz) = step_cosine( tmp, 0.5_rk*params_acm%L_sponge, 0.5_rk*params_acm%L_sponge)
                end do
            end do
        end do

    elseif (params_acm%sponge_type == "inlet-outlet-x") then
        ! outlet sponge in x-direction
        do ix = g+1, Bs(1)+g
            x = dble(ix-(g+1)) * dx(1) + x0(1)

            ! distance to borders of domain
            tmp = minval( (/x-global_min(1),-(x-global_max(1))/) )

            sponge(ix,:,:) = step_cosine( tmp, 0.5_rk*params_acm%L_sponge, 0.5_rk*params_acm%L_sponge)
        end do


    elseif (params_acm%sponge_type == "p-norm") then
        ! p-norm sponge. The shape of the sponge is dictated as the p-norm
        ! https://de.wikipedia.org/wiki/P-Norm
        ! which is a nice and simple way to get a rectangle with round corners.
        ! For 2 it is the eucledian norm (circle), for infinity it is a rectangle and for everything in between a rounded rectangle

        if ( maxval(abs((global_max-global_min) - (global_max(1)-global_min(1)))) > 1.0e-10_rk) then
            call abort(1610184,"ERROR: for the p-norm sponge, the (cropped) domain has to be same size in all directions.")
        endif

        p = params_acm%p_sponge
        pinv = 1.0_rk / p
        center = 0.5_rk * (global_min + global_max)
        radius = 0.5_rk * (global_max(1) - global_min(1))

        ! check if this block is in sponge layer
        ! Find the point in the block furthest from the domain center
        x = max(abs(x0(1)-center(1)), abs(x0(1)+dx(1)*dble(Bs(1))-center(1)))**p
        y = max(abs(x0(2)-center(2)), abs(x0(2)+dx(2)*dble(Bs(2))-center(2)))**p
        z = max(abs(x0(3)-center(3)), abs(x0(3)+dx(3)*dble(Bs(3))-center(3)))**p
        tmp = radius - (x + y + z)**pinv  ! we invert distance to center to distance from boundary
        if (tmp > params_acm%L_sponge) then
            sponge = 0.0_rk
        else
            do iz = g+1, Bs(3)+g
                z = (dble(iz-(g+1)) * dx(3) + x0(3) - center(3))**p
                do iy = g+1, Bs(2)+g
                    y = (dble(iy-(g+1)) * dx(2) + x0(2) - center(2))**p
                    do ix = g+1, Bs(1)+g
                        x = (dble(ix-(g+1)) * dx(1) + x0(1) - center(1))**p

                        ! distance to borders of domain
                        tmp = -( (x + y + z)**pinv - radius)

                        sponge(ix,iy,iz) = step_cosine( tmp, 0.5_rk*params_acm%L_sponge, &
                        0.5_rk*params_acm%L_sponge)
                    end do
                end do
            end do
        endif

    elseif (params_acm%sponge_type == "p-norm-insect-centered") then
        ! p-norm sponge. The shape of the sponge is dictated as the p-norm
        ! https://de.wikipedia.org/wiki/P-Norm
        ! which is a nice and simple way to get a rectangle with round corners.
        ! This sponge type is moving with the insect, which we assume to be in the centre
        ! of it

        ! Attention: For multiple insects, this is always centered around the first insect in the list, as for multiple insects this sponge types makes only little sense anyways

        if ( maxval(abs((global_max-global_min) - (global_max(1)-global_min(1)))) > 1.0e-10_rk) then
            call abort(1610184,"ERROR: for the p-norm sponge, the (cropped) domain has to be same size in all directions.")
        endif

        p = params_acm%p_sponge
        pinv = 1.0_rk / p
        radius = 0.5_rk * (global_max(1) - global_min(1))

        do iz = g+1, Bs(3)+g
            z = (dble(iz-(g+1)) * dx(3) + x0(3) - Insects(1)%xc_body_g(3))**p
            do iy = g+1, Bs(2)+g
                y = (dble(iy-(g+1)) * dx(2) + x0(2) - Insects(1)%xc_body_g(2))**p
                do ix = g+1, Bs(1)+g
                    x = (dble(ix-(g+1)) * dx(1) + x0(1) - Insects(1)%xc_body_g(1))**p

                    ! distance to borders of domain
                    tmp = -( (x + y + z)**pinv - radius)

                    sponge(ix,iy,iz) = step_cosine( tmp, 0.5_rk*params_acm%L_sponge, &
                    0.5_rk*params_acm%L_sponge)
                end do
            end do
        end do

    else
        call abort(1610181,"Sponge-type is unknown")
    endif

end subroutine sponge_3D
