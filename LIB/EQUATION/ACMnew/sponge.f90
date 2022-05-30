!! input:    - params, origin and spacing of the block, grid parameters \n
!! output:   - sponge term \n
!*********************************************************************************************
subroutine sponge_2D(sponge, x0, dx, Bs, g)
    implicit none

    integer(kind=ik), intent(in)                   :: g                         ! grid
    integer(kind=ik), dimension(3), intent(in)     :: Bs
    real(kind=rk), dimension(:,:), intent(out)     :: sponge                    !> sponge term for every grid point of this block
    real(kind=rk), dimension(2), intent(in)        :: x0, dx                    !> spacing and origin of block
    real(kind=rk)                                  :: x, y, tmp, p, offset      ! auxiliary variables
    integer(kind=ik)                               :: ix, iy                    ! loop variables

    if (.not. params_acm%initialized) write(*,*) "WARNING: sponge_2D called but ACM not initialized"

    if (params_acm%sponge_type == "rect") then
        ! rectangular sponge with 45deg edges
        do iy = g+1, Bs(2)+g+ONE_SKIPREDUNDANT
            y = dble(iy-(g+1)) * dx(2) + x0(2)

            do ix = g+1, Bs(1)+g+ONE_SKIPREDUNDANT
                x = dble(ix-(g+1)) * dx(1) + x0(1)

                ! distance to borders of domain
                tmp = minval( (/x,y,-(x-params_acm%domain_size(1)),-(y-params_acm%domain_size(2))/) )

                sponge(ix,iy) = smoothstep( tmp, 0.5_rk*params_acm%L_sponge, 0.5_rk*params_acm%L_sponge)
            end do
        end do

    elseif (params_acm%sponge_type == "p-norm") then
        ! p-norm sponge. The shape of the sponge is dictated as the p-norm
        ! https://de.wikipedia.org/wiki/P-Norm
        ! which is a nice and simple way to get a rectangle with round corners.

        if ( maxval(abs(params_acm%domain_size(1:2)-params_acm%domain_size(1))) > 1.0e-10_rk) then
            call abort(1610184,"ERROR: for the p-norm sponge, the domain has to be same size in all directions.")
        endif

        p = params_acm%p_sponge
        offset = 0.5_rk * params_acm%domain_size(1)

        do iy = g+1, Bs(2)+g+ONE_SKIPREDUNDANT
            y = dble(iy-(g+1)) * dx(2) + x0(2) - offset
            do ix = g+1, Bs(1)+g+ONE_SKIPREDUNDANT
                x = dble(ix-(g+1)) * dx(1) + x0(1) - offset

                ! distance to borders of domain
                tmp = -( (x**p + y**p)**(1.0_rk/p) - offset)
                sponge(ix,iy) = smoothstep( tmp, 0.5_rk*params_acm%L_sponge, &
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
    real(kind=rk)     :: x, y, z, tmp, p, offset, pinv
    ! loop variables
    integer(kind=ik)  :: ix, iy, iz

    if (.not. params_acm%initialized) write(*,*) "WARNING: sponge_3D called but ACM not initialized"


    if (params_acm%sponge_type == "rect") then
        ! rectangular sponge with 45deg edges
        do iz = g+1, Bs(3)+g+ONE_SKIPREDUNDANT
            z = dble(iz-(g+1)) * dx(3) + x0(3)
            do iy = g+1, Bs(2)+g+ONE_SKIPREDUNDANT
                y = dble(iy-(g+1)) * dx(2) + x0(2)
                do ix = g+1, Bs(1)+g+ONE_SKIPREDUNDANT
                    x = dble(ix-(g+1)) * dx(1) + x0(1)

                    ! distance to borders of domain
                    tmp = minval( (/x,y,z,-(x-params_acm%domain_size(1)),&
                         -(y-params_acm%domain_size(2)),-(z-params_acm%domain_size(3))/) )

                    sponge(ix,iy,iz) = smoothstep( tmp, 0.5_rk*params_acm%L_sponge, 0.5_rk*params_acm%L_sponge)
                end do
            end do
        end do

        ! sponge for using with symmetry_BC
        ! insect is supposed to be at y=0
    elseif (params_acm%sponge_type == "rect-symmetry-y") then
        ! rectangular sponge with 45deg edges
        do iz = g+1, Bs(3)+g+ONE_SKIPREDUNDANT
            z = dble(iz-(g+1)) * dx(3) + x0(3)
            do iy = g+1, Bs(2)+g+ONE_SKIPREDUNDANT
                y = dble(iy-(g+1)) * dx(2) + x0(2)
                do ix = g+1, Bs(1)+g+ONE_SKIPREDUNDANT
                    x = dble(ix-(g+1)) * dx(1) + x0(1)

                    ! distance to borders of domain
                    tmp = minval( (/x,z,-(x-params_acm%domain_size(1)),-(y-params_acm%domain_size(2)),-(z-params_acm%domain_size(3))/) )

                    sponge(ix,iy,iz) = smoothstep( tmp, 0.5_rk*params_acm%L_sponge, 0.5_rk*params_acm%L_sponge)
                end do
            end do
        end do

        ! sponge for using with symmetry_BC
        ! insect is supposed to be at y=0 z=0
    elseif (params_acm%sponge_type == "rect-symmetry-yz") then
        ! rectangular sponge with 45deg edges
        do iz = g+1, Bs(3)+g+ONE_SKIPREDUNDANT
            z = dble(iz-(g+1)) * dx(3) + x0(3)
            do iy = g+1, Bs(2)+g+ONE_SKIPREDUNDANT
                y = dble(iy-(g+1)) * dx(2) + x0(2)
                do ix = g+1, Bs(1)+g+ONE_SKIPREDUNDANT
                    x = dble(ix-(g+1)) * dx(1) + x0(1)

                    ! distance to borders of domain
                    tmp = minval( (/x,-(x-params_acm%domain_size(1)),-(y-params_acm%domain_size(2)),-(z-params_acm%domain_size(3))/) )

                    sponge(ix,iy,iz) = smoothstep( tmp, 0.5_rk*params_acm%L_sponge, 0.5_rk*params_acm%L_sponge)
                end do
            end do
        end do

    elseif (params_acm%sponge_type == "inlet-outlet-x") then
        ! outlet sponge in x-direction
        do ix = g+1, Bs(1)+g+ONE_SKIPREDUNDANT
            x = dble(ix-(g+1)) * dx(1) + x0(1)

            ! distance to borders of domain
            tmp = minval( (/x,-(x-params_acm%domain_size(1))/) )

            sponge(ix,:,:) = smoothstep( tmp, 0.5_rk*params_acm%L_sponge, 0.5_rk*params_acm%L_sponge)
        end do


    elseif (params_acm%sponge_type == "p-norm") then
        ! p-norm sponge. The shape of the sponge is dictated as the p-norm
        ! https://de.wikipedia.org/wiki/P-Norm
        ! which is a nice and simple way to get a rectangle with round corners.

        if ( maxval(abs(params_acm%domain_size-params_acm%domain_size(1))) > 1.0e-10_rk) then
            call abort(1610184,"ERROR: for the p-norm sponge, the domain has to be same size in all directions.")
        endif

        p = params_acm%p_sponge
        pinv = 1.0_rk / p
        offset = 0.5_rk * params_acm%domain_size(1)

        do iz = g+1, Bs(3)+g+ONE_SKIPREDUNDANT
            z = (dble(iz-(g+1)) * dx(3) + x0(3) - offset)**p
            do iy = g+1, Bs(2)+g+ONE_SKIPREDUNDANT
                y = (dble(iy-(g+1)) * dx(2) + x0(2) - offset)**p
                do ix = g+1, Bs(1)+g+ONE_SKIPREDUNDANT
                    x = (dble(ix-(g+1)) * dx(1) + x0(1) - offset)**p

                    ! distance to borders of domain
                    tmp = -( (x + y + z)**pinv - offset)

                    sponge(ix,iy,iz) = smoothstep( tmp, 0.5_rk*params_acm%L_sponge, &
                    0.5_rk*params_acm%L_sponge)
                end do
            end do
        end do

    else
        call abort(1610181,"Sponge-type is unknown")
    endif

end subroutine sponge_3D
