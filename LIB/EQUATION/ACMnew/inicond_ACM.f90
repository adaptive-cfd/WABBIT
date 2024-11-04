!-----------------------------------------------------------------------------
! main level wrapper for setting the initial condition on a block
!-----------------------------------------------------------------------------
subroutine INICOND_ACM( time, u, g, x0, dx, n_domain )
    implicit none

    ! it may happen that some source terms have an explicit time-dependency
    ! therefore the general call has to pass time
    real(kind=rk), intent (in) :: time

    ! block data, containg the state vector. In general a 4D field (3 dims+components)
    ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
    real(kind=rk), intent(inout) :: u(1:,1:,1:,1:)

    ! as you are allowed to compute the RHS only in the interior of the field
    ! you also need to know where 'interior' starts: so we pass the number of ghost points
    integer, intent(in) :: g

    ! for each block, you'll need to know where it lies in physical space. The first
    ! non-ghost point has the coordinate x0, from then on its just cartesian with dx spacing
    real(kind=rk), intent(in) :: x0(1:3), dx(1:3)

    ! when implementing boundary conditions, it is necessary to know if the local field (block)
    ! is adjacent to a boundary, because the stencil has to be modified on the domain boundary.
    ! The n_domain tells you if the local field is adjacent to a domain boundary:
    ! n_domain(i) can be either 0, 1, -1,
    !  0: no boundary in the direction +/-e_i
    !  1: boundary in the direction +e_i
    ! -1: boundary in the direction - e_i
    ! currently only acessible in the local stage
    ! NOTE: ACM only supports symmetry BC for the moment (which is handled by wabbit and not ACM)
    integer(kind=2), intent(in) :: n_domain(3)

    real(kind=rk)    :: x, y, z
    integer(kind=ik) :: ix, iy, iz, idir, Bs(3), iscalar, ix_global, iy_global

    ! compute the size of blocks
    Bs(1) = size(u,1) - 2*g
    Bs(2) = size(u,2) - 2*g
    Bs(3) = size(u,3) - 2*g

    u = 0.0_rk

    if (.not. params_acm%initialized) write(*,'(A)') "WARNING: INICOND_ACM called but ACM not initialized"

    if (params_acm%dim==2 .and. size(u,4) /= params_acm%dim + 1 + params_acm%N_scalars) then
        call abort(23091801,"ACM: state vector has not the right number of components")
    endif

    select case (params_acm%inicond)
    case ("noise")
        call random_data(u)
        u = u * params_acm%beta

    case ("lamballais")
        if (params_acm%dim /= 2) call abort(1409241, "lamballais is a 2D test case")

        do iy = g+1, Bs(2)+g
            y = dble(iy-(g+1)) * dx(2) + x0(2)
            do ix = g+1, Bs(1)+g
                x = dble(ix-(g+1)) * dx(1) + x0(1)
    
                ix_global = int( x/dx(1) )
                iy_global = int( y/dx(2) )
     
                ! fields are initialized in read_params_ACM
                ! Note inefficiently, each mpirank has the full array (does not matter as is 2D)
                u(ix,iy,:,1) = params_acm%u_lamballais(ix_global, iy_global, 1) ! ux
                u(ix,iy,:,2) = params_acm%u_lamballais(ix_global, iy_global, 2) ! uy
                u(ix,iy,:,3) = params_acm%u_lamballais(ix_global, iy_global, 3) ! p
            end do
        end do

    case("pressure-blob")
        if (params_acm%dim==2) then
            ! create gauss pulse. Note we loop over the entire block, incl. ghost nodes.
            do iy = 1, Bs(2)+2*g
                do ix = 1, Bs(1)+2*g
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
            do iz = 1, Bs(3)+2*g
                do iy = 1, Bs(2)+2*g
                    do ix = 1, Bs(1)+2*g
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

    case ("taylor-green-vanRees2011")
        do iz = 1, Bs(3)+2*g
            do iy = 1, Bs(2)+2*g
                do ix = 1, Bs(1)+2*g
                    ! compute x,y coordinates from spacing and origin
                    x = dble(ix-(g+1)) * dx(1) + x0(1)
                    y = dble(iy-(g+1)) * dx(2) + x0(2)
                    z = dble(iz-(g+1)) * dx(3) + x0(3)

                    ! the initial condition is known analytically. Note in vanrees JCP2011,
                    ! there is a parameter \theta which is set to zero.
                    ! See also: "Problem C3.5 Direct Numerical Simulation of the Taylor-Green Vortex at Re = 1600"
                    ! the constants rho_0 and p_0 are set to 1.0 and 0.0, respectively
                    ! NOTE: domain size has to be 2*pi = 6.283185307179586
                    u(ix, iy, iz, 1) = dsin(x)*dcos(y)*dcos(z)
                    u(ix, iy, iz, 2) =-dcos(x)*dsin(y)*dcos(z)
                    u(ix, iy, iz, 3) = 0.0_rk
                    u(ix, iy, iz, 4) = (dcos(2.0_rk*x) + dcos(2.0_rk*y))*(dcos(2.0_rk*z) + 2.0_rk) / 16.0_rk
                end do
            end do
        end do

    case("jet-x")
        ! a jet with constant (unit) velocity, with a bit of noise to trigger the instability
        do iy = 1, Bs(2)+2*g
            ! compute x,y coordinates from spacing and origin
            y = abs(dble(iy-(g+1)) * dx(2) + x0(2) - params_acm%domain_size(2)/2.0_rk)

            if (y <= 0.25_rk * params_acm%domain_size(2)) then
                ! ux is =1
                u(:,iy,:,1) = 1.0
            endif
        end do

        ! noise in uy
        ! do iz = 1, size(u,3)
        !     do iy = 1, size(u,2)
        !         do ix = 1, size(u,1)
        !             u(ix,iy,iz,2) = u(ix,iy,iz,2) + 1.0e-5_rk * rand_nbr()
        !         enddo
        !     enddo
        ! enddo

    case("velocity-blob")
        if (params_acm%dim==2) then
            ! create gauss pulse. Note we loop over the entire block, incl. ghost nodes.
            do iy = 1, Bs(2)+2*g
                do ix = 1, Bs(1)+2*g
                    ! compute x,y coordinates from spacing and origin
                    x = dble(ix-(g+1)) * dx(1) + x0(1) - params_acm%domain_size(1)/2.0_rk
                    y = dble(iy-(g+1)) * dx(2) + x0(2) - params_acm%domain_size(2)/2.0_rk

                    if (x<-params_acm%domain_size(1)/2.0) x = x + params_acm%domain_size(1)
                    if (x>params_acm%domain_size(1)/2.0) x = x - params_acm%domain_size(1)

                    if (y<-params_acm%domain_size(2)/2.0) y = y + params_acm%domain_size(2)
                    if (y>params_acm%domain_size(2)/2.0) y = y - params_acm%domain_size(2)

                    ! set actual inicond gauss blob
                    ! here only for the pressure.
                    u(ix,iy,:,1:2) = dexp( -( (x)**2 + (y)**2 ) / params_acm%beta )
                end do
            end do
        else
            ! create gauss pulse
            do iz = 1, Bs(3)+2*g
                do iy = 1, Bs(2)+2*g
                    do ix = 1, Bs(1)+2*g
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
                        u(ix,iy,iz,1:3) = dexp( -( (x)**2 + (y)**2 + (z)**2 ) / params_acm%beta )
                    end do
                end do
            end do
        end if

    case("ux-blob")
        if (params_acm%dim==2) then
            ! create gauss pulse. Note we loop over the entire block, incl. ghost nodes.
            do iy = 1, Bs(2)+2*g
                do ix = 1, Bs(1)+2*g
                    ! compute x,y coordinates from spacing and origin
                    x = dble(ix-(g+1)) * dx(1) + x0(1) - params_acm%domain_size(1)/2.0_rk
                    y = dble(iy-(g+1)) * dx(2) + x0(2) - params_acm%domain_size(2)/2.0_rk

                    if (x<-params_acm%domain_size(1)/2.0) x = x + params_acm%domain_size(1)
                    if (x>params_acm%domain_size(1)/2.0) x = x - params_acm%domain_size(1)

                    if (y<-params_acm%domain_size(2)/2.0) y = y + params_acm%domain_size(2)
                    if (y>params_acm%domain_size(2)/2.0) y = y - params_acm%domain_size(2)

                    ! set actual inicond gauss blob
                    ! here only for the pressure.
                    u(ix,iy,:,1) = dexp( -( (x)**2 + (y)**2 ) / params_acm%beta )
                end do
            end do
        else
            ! create gauss pulse
            do iz = 1, Bs(3)+2*g
                do iy = 1, Bs(2)+2*g
                    do ix = 1, Bs(1)+2*g
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
                        u(ix,iy,iz,1) = dexp( -( (x)**2 + (y)**2 + (z)**2 ) / params_acm%beta )
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
            do iy= 1,Bs(2)+2*g
                do ix= 1, Bs(1)+2*g
                    x = x0(1) + dble(ix-g-1)*dx(1)
                    y = x0(2) + dble(iy-g-1)*dx(2)
                    u(ix,iy,1,1) = sin( 2.0_rk*pi*x/params_acm%domain_size(1) ) + 0.5_rk*sin( 10.0_rk*pi*x/params_acm%domain_size(1) )
                    u(ix,iy,1,2) = cos( 2.0_rk*pi*y/params_acm%domain_size(2) )*sin( 2.0_rk*pi*x/params_acm%domain_size(1) )
                enddo
            enddo
        else
            do iz = 1, Bs(3)+2*g
                do iy= 1,Bs(2)+2*g
                    do ix= 1, Bs(1)+2*g
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
        do iy= 1,Bs(2)+2*g
            do ix= 1, Bs(1)+2*g
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

    ! --------------------------------------------------------------------------
    ! initial conditions for passive scalars, if used.
    ! --------------------------------------------------------------------------
    if (params_acm%use_passive_scalar) then
        ! loop over scalars
        do iscalar = 1, params_acm%N_scalars
            select case (params_acm%scalar_inicond(iscalar))
            case ("empty", "none", "zero")
                u(:,:,:,params_acm%dim + 1 + iscalar) = 0.0_rk
            case ("Kadoch2012")
                if (params_acm%dim == 2) then
                    do iy = 1, Bs(2)+2*g
                        do ix = 1, Bs(1)+2*g
                            x = x0(1) + dble(ix-g-1)*dx(1) - 1.0_rk ! domain is -1...1 in kadoch
                            y = x0(2) + dble(iy-g-1)*dx(2) - 1.0_rk

                            x = x - params_acm%length ! finite size of cavity
                            y = y - params_acm%length

                            ! in their original work, they set the initial condition
                            ! everywhere in the domain (even in the penalization layer)
                            u(ix,iy,:,params_acm%dim + 1 + iscalar) = cos(pi*y)*(cos(4.0_rk*pi*x)+cos(pi*x))
                        end do
                    end do
                else
                    call abort(0409191, "Scalar inicond Kadoch2012 is only for 2D")
                endif
            case default
                call abort(0409192, "Unkown scalar inicond")
            end select
        enddo
    endif

end subroutine INICOND_ACM
