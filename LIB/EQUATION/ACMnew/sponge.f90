!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name sponge.f90
!> \version 0.5
!> \author sm
!
!> \brief
!
!>
!! input:    - params, origin and spacing of the block, grid parameters \n
!! output:   - sponge term \n
!!
!!
!! = log ======================================================================================
!! \n
!! 12/18 - create
!*********************************************************************************************
subroutine sponge_2D(sponge, x0, dx, Bs, g)
    implicit none

    ! grid
    integer(kind=ik), intent(in)                   :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    !> sponge term for every grid point of this block
    real(kind=rk), dimension(:,:), intent(out)     :: sponge
    !> spacing and origin of block
    real(kind=rk), dimension(2), intent(in)        :: x0, dx

    ! auxiliary variables
    real(kind=rk)    :: x, y, tmp, p, offset
    ! loop variables
    integer(kind=ik) :: ix, iy

    if (params_acm%sponge_type == "rect") then
        ! rectangular sponge with 45deg edges
        do iy = g+1, Bs(2)+g
            y = dble(iy-(g+1)) * dx(2) + x0(2)

            do ix = g+1, Bs(1)+g
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

        do iy = g+1, Bs(2)+g
            y = dble(iy-(g+1)) * dx(2) + x0(2) - offset
            do ix = g+1, Bs(1)+g
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
    real(kind=rk)     :: x, y, z, tmp, p, offset
    ! loop variables
    integer(kind=ik)  :: ix, iy, iz


    if (params_acm%sponge_type == "rect") then
        ! rectangular sponge with 45deg edges
        do iz = g+1, Bs(3)+g
            z = dble(iz-(g+1)) * dx(3) + x0(3)
            do iy = g+1, Bs(2)+g
                y = dble(iy-(g+1)) * dx(2) + x0(2)
                do ix = g+1, Bs(1)+g
                    x = dble(ix-(g+1)) * dx(1) + x0(1)

                    ! distance to borders of domain
                    tmp = minval( (/x,y,z,-(x-params_acm%domain_size(1)),&
                         -(y-params_acm%domain_size(2)),-(z-params_acm%domain_size(3))/) )

                    sponge(ix,iy,iz) = smoothstep( tmp, 0.5_rk*params_acm%L_sponge, 0.5_rk*params_acm%L_sponge)
                end do
            end do
        end do

    elseif (params_acm%sponge_type == "p-norm") then
        ! p-norm sponge. The shape of the sponge is dictated as the p-norm
        ! https://de.wikipedia.org/wiki/P-Norm
        ! which is a nice and simple way to get a rectangle with round corners.

        if ( maxval(abs(params_acm%domain_size-params_acm%domain_size(1))) > 1.0e-10_rk) then
            call abort(1610184,"ERROR: for the p-norm sponge, the domain has to be same size in all directions.")
        endif

        p = params_acm%p_sponge
        offset = 0.5_rk * params_acm%domain_size(1)

        do iz = g+1, Bs(3)+g
            z = dble(iz-(g+1)) * dx(3) + x0(3) - offset
            do iy = g+1, Bs(2)+g
                y = dble(iy-(g+1)) * dx(2) + x0(2) - offset
                do ix = g+1, Bs(1)+g
                    x = dble(ix-(g+1)) * dx(1) + x0(1) - offset

                    ! distance to borders of domain
                    tmp = -( (x**p + y**p + z**p)**(1.0_rk/p) - offset)

                    sponge(ix,iy,iz) = smoothstep( tmp, 0.5_rk*params_acm%L_sponge, &
                    0.5_rk*params_acm%L_sponge)
                end do
            end do
        end do

    else
        call abort(1610181,"Sponge-type is unknown")
    endif

end subroutine sponge_3D
