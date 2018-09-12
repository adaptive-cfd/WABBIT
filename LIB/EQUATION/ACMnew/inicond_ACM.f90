!-----------------------------------------------------------------------------
! main level wrapper for setting the initial condition on a block
!-----------------------------------------------------------------------------
subroutine INICOND_ACM( time, u, g, x0, dx, work, adapting )
    implicit none

    ! it may happen that some source terms have an explicit time-dependency
    ! therefore the general call has to pass time
    real(kind=rk), intent (in) :: time

    ! block data, containg the state vector. In general a 4D field (3 dims+components)
    ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
    real(kind=rk), intent(inout) :: u(1:,1:,1:,1:)

    ! work data, for mask, vorticity etc. In general a 4D field (3 dims+components)
    ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
    real(kind=rk), intent(inout) :: work(1:,1:,1:,1:)

    ! as you are allowed to compute the RHS only in the interior of the field
    ! you also need to know where 'interior' starts: so we pass the number of ghost points
    integer, intent(in) :: g

    ! for each block, you'll need to know where it lies in physical space. The first
    ! non-ghost point has the coordinate x0, from then on its just cartesian with dx spacing
    real(kind=rk), intent(in) :: x0(1:3), dx(1:3)

    ! if we are still adapting the initial condition, we may use penalization for refinement.
    ! if the initial grid is adapted we set our initial condition without penalization (impulsive start).
    logical, intent(in) :: adapting

    real(kind=rk)    :: x,y,z
    integer(kind=ik) :: Bs, ix, iy, iz, idir
    ! we have quite some of these work arrays in the code, but they are very small,
    ! only one block. They're ngeligible in front of the lgt_block array.
    real(kind=rk), allocatable, save :: mask(:,:,:), us(:,:,:,:)

    ! compute the size of blocks
    Bs = size(u,1) - 2*g

    u = 0.0_rk

    select case (params_acm%inicond)
    case("pressure-blob")
        if (params_acm%dim==2) then
            ! create gauss pulse. Note we loop over the entire block, incl. ghost nodes.
            do iy = 1, Bs+2*g
                do ix = 1, Bs+2*g
                    ! compute x,y coordinates from spacing and origin
                    x = dble(ix-(g+1)) * dx(1) + x0(1) - params_acm%domain_size(1)/2.0_rk
                    y = dble(iy-(g+1)) * dx(2) + x0(2) - params_acm%domain_size(2)/2.0_rk

                    if (x<-params_acm%domain_size(1)/2.0) x = x + params_acm%domain_size(1)
                    if (x>params_acm%domain_size(1)/2.0) x = x - params_acm%domain_size(1)

                    if (y<-params_acm%domain_size(2)/2.0) y = y + params_acm%domain_size(2)
                    if (y>params_acm%domain_size(2)/2.0) y = y - params_acm%domain_size(2)

                    ! set actual inicond gauss blob
                    ! here only for the pressure.
                    u(ix,iy,:,size(u,4)) = dexp( -( (x)**2 + (y)**2 ) / params_acm%beta )
                end do
            end do
        else
            ! create gauss pulse
            do iz = 1, Bs+2*g
                do iy = 1, Bs+2*g
                    do ix = 1, Bs+2*g
                        ! compute x,y coordinates from spacing and origin
                        x = dble(ix-(g+1)) * dx(1) + x0(1) - params_acm%domain_size(1)/2.0_rk
                        y = dble(iy-(g+1)) * dx(2) + x0(2) - params_acm%domain_size(2)/2.0_rk
                        z = dble(iz-(g+1)) * dx(3) + x0(3) - params_acm%domain_size(3)/2.0_rk

                        if (x<-params_acm%domain_size(1)/2.0) x = x + params_acm%domain_size(1)
                        if (x>params_acm%domain_size(1)/2.0) x = x - params_acm%domain_size(1)

                        if (y<-params_acm%domain_size(2)/2.0) y = y + params_acm%domain_size(2)
                        if (y>params_acm%domain_size(2)/2.0) y = y - params_acm%domain_size(2)

                        if (z<-params_acm%domain_size(3)/2.0) z = z + params_acm%domain_size(3)
                        if (z>params_acm%domain_size(3)/2.0) z = z - params_acm%domain_size(3)

                        ! set actual inicond gauss blob
                        u(ix,iy,iz,size(u,4)) = dexp( -( (x)**2 + (y)**2 + (z)**2 ) / params_acm%beta )
                    end do
                end do
            end do
        end if

    case("meanflow")
        do idir = 1, params_acm%dim
            u(:,:,:,idir) = params_acm%u_mean_set(idir)
        enddo

    case("sinewaves-nopress")
        ! some random sine waves, but no pressure imposed.
        if (params_acm%dim == 2) then
            do iy= 1,Bs+2*g
                do ix= 1, Bs+2*g
                    x = x0(1) + dble(ix-g-1)*dx(1)
                    y = x0(2) + dble(iy-g-1)*dx(2)
                    u(ix,iy,1,1) = sin( 2.0_rk*pi*x/params_acm%domain_size(1) ) + 0.5_rk*sin( 10.0_rk*pi*x/params_acm%domain_size(1) )
                    u(ix,iy,1,2) = cos( 2.0_rk*pi*y/params_acm%domain_size(2) )*sin( 2.0_rk*pi*x/params_acm%domain_size(1) )
                enddo
            enddo
        else
            do iz = 1, Bs+2*g
                do iy= 1,Bs+2*g
                    do ix= 1, Bs+2*g
                        x = x0(1) + dble(ix-g-1)*dx(1)
                        y = x0(2) + dble(iy-g-1)*dx(2)
                        z = x0(3) + dble(iz-g-1)*dx(3)

                        u(ix,iy,iz,1) = sin( 2.0_rk*pi*x/params_acm%domain_size(1) ) + 0.5_rk*sin( 10.0_rk*pi*x/params_acm%domain_size(1) )
                        u(ix,iy,iz,2) = cos( 2.0_rk*pi*y/params_acm%domain_size(2) )*sin( 2.0_rk*pi*x/params_acm%domain_size(1) )*sin( 2.0_rk*pi*z/params_acm%domain_size(3) )
                        u(ix,iy,iz,3) = cos( 2.0_rk*pi*y/params_acm%domain_size(2) )*sin( 3.0_rk*pi*x/params_acm%domain_size(1) )*cos( 2.0_rk*pi*z/params_acm%domain_size(3) )
                    enddo
                enddo
            enddo

        endif

    case("taylor_green")
        do iy= 1,Bs+2*g
            do ix= 1, Bs+2*g
                x = x0(1) + dble(ix-g-1)*dx(1)
                y = x0(2) + dble(iy-g-1)*dx(2)

                call continue_periodic(x,params_acm%domain_size(1))
                call continue_periodic(y,params_acm%domain_size(2))

                u(ix,iy,1,1) = params_acm%u_mean_set(1) + dsin(x)*dcos(y)
                u(ix,iy,1,2) = params_acm%u_mean_set(2) - dcos(x)*dsin(y)
                u(ix,iy,1,3) = 0.25_rk*(dcos(2.0_rk*x) + dcos(2.0_rk*y))
            end do
        end do

    case default
        call abort(428764, "ACM inicond: "//trim(adjustl(params_acm%inicond))//" is unkown.")

    end select

    ! if we use volume penalization, the mask is first used for refinement of the grid.
    ! In a second stage, the initial condition without penalization is then applied to the refined grid.
    if (adapting .and. params_acm%penalization) then

        if (params_acm%dim==3) then
            if (.not. allocated(mask)) allocate(mask(1:Bs+2*g, 1:Bs+2*g, 1:Bs+2*g))
            if (.not. allocated(us)) allocate(us(1:Bs+2*g, 1:Bs+2*g, 1:Bs+2*g, 1:3))
        else
            if (.not. allocated(mask)) allocate(mask(1:Bs+2*g, 1:Bs+2*g, 1))
            if (.not. allocated(us)) allocate(us(1:Bs+2*g, 1:Bs+2*g, 1, 1:2))
        endif


        if (params_acm%dim == 2) then
            call create_mask_2D(time, x0, dx, Bs, g, mask(:,:,1), us(:,:,1,1:2))
        else
            call create_mask_3D(time, x0, dx, Bs, g, mask, us)
        endif

        do idir = 1, params_acm%dim
            u(:,:,:,idir) = mask
        enddo
    end if

end subroutine INICOND_ACM
