!-----------------------------------------------------------------------------
! main level wrapper for setting the initial condition on a block
!
! NOTE: This implementation is for the PRESSURE POISSON formulation.
! Only VELOCITY components are initialized here. Pressure is computed
! separately via the Poisson equation from the velocity field.
!-----------------------------------------------------------------------------
subroutine INICOND_NSPP( time, u, g, x0, dx, n_domain )
    use module_helpers
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
    ! NOTE: NSPP only supports symmetry BC for the moment (which is handled by wabbit and not NSPP)
    integer(kind=2), intent(in) :: n_domain(3)

    real(kind=rk)    :: x, y, z, param_1, param_2
    integer(kind=ik) :: ix, iy, iz, idir, Bs(3), iscalar, ix_global, iy_global

    ! compute the size of blocks
    Bs(1) = size(u,1) - 2*g
    Bs(2) = size(u,2) - 2*g
    Bs(3) = size(u,3) - 2*g

    u = 0.0_rk

    if (.not. params_nspp%initialized) write(*,'(A)') "WARNING: INICOND_NSPP called but NSPP not initialized"

    if (params_nspp%dim==2 .and. size(u,4) /= params_nspp%dim + 1 + params_nspp%N_scalars + params_nspp%N_time_statistics) then
        call abort(23091801,"NSPP: state vector has not the right number of components")
    endif

    select case (params_nspp%inicond)
    case ("noise")
        call random_data(u)
        u = u * params_nspp%beta

    case ("noise2")
        call random_data(u)
        u = 2.0_rk*(u - 0.5_rk) * params_nspp%beta

    case ("noise_unique")
	! Random data, but results in the same random data on a given grid every time it is called
        call random_data_unique(u, x0, dx, (/g,g,g/), params_nspp%domain_size, params_nspp%Jmax, Bs )
        u = 2.0_rk*(u - 0.5_rk) * params_nspp%beta

    case ("lamballais")
        if (params_nspp%dim /= 2) call abort(1409241, "lamballais is a 2D test case")

        do iy = g+1, Bs(2)+g
            y = dble(iy-(g+1)) * dx(2) + x0(2)
            do ix = g+1, Bs(1)+g
                x = dble(ix-(g+1)) * dx(1) + x0(1)
    
                ix_global = int( x/dx(1) )
                iy_global = int( y/dx(2) )
     
                ! fields are initialized in read_params_NSPP
                ! Note inefficiently, each mpirank has the full array (does not matter as is 2D)
                u(ix,iy,:,1) = params_nspp%u_lamballais(ix_global, iy_global, 1) ! ux
                u(ix,iy,:,2) = params_nspp%u_lamballais(ix_global, iy_global, 2) ! uy
                ! Note: pressure component (u_lamballais(:,:,3)) is not used in pressure Poisson approach
                ! Pressure will be computed from velocity field via Poisson equation
            end do
        end do

    case("jet-x")
        ! a jet with constant (unit) velocity, with a bit of noise to trigger the instability
        do iy = 1, Bs(2)+2*g
            ! compute x,y coordinates from spacing and origin
            y = abs(dble(iy-(g+1)) * dx(2) + x0(2) - params_nspp%domain_size(2)/2.0_rk)

            if (y <= 0.25_rk * params_nspp%domain_size(2)) then
                ! ux is =1
                u(:,iy,:,1) = 1.0
            endif
        end do

        ! noise in uy
        do iz = 1, size(u,3)
            do iy = 1, size(u,2)
                do ix = 1, size(u,1)
                    u(ix,iy,iz,2) = u(ix,iy,iz,2) + params_nspp%beta * rand_nbr()
                enddo
            enddo
        enddo

    case("jet-x-sin")
        ! a jet with constant (unit) velocity
        ! no smoothing is applied
        do iy = 1, Bs(2)+2*g
            ! compute x,y coordinates from spacing and origin
            y = abs(dble(iy-(g+1)) * dx(2) + x0(2) - params_nspp%domain_size(2)/2.0_rk)

            if (y <= 0.25_rk * params_nspp%domain_size(2)) then
                ! ux is =1
                u(:,iy,:,1) = 1.0
            endif
        end do

        ! uy is a sin wave (that triggers the instability)
        do ix = 1, Bs(1)+2*g
            ! compute x,y coordinates from spacing and origin
            x = dble(ix-(g+1)) * dx(1) + x0(1)

            u(ix,:,:,2) = 0.1_rk * sin(2.0_rk*pi*x/params_nspp%domain_size(1))
        end do

    case("jet-x-sin-smooth")
        ! a jet with constant (unit) velocity
        ! with smoothing applied
        do iy = 1, Bs(2)+2*g
            ! compute x,y coordinates from spacing and origin
            y = abs(dble(iy-(g+1)) * dx(2) + x0(2) - params_nspp%domain_size(2)/2.0_rk)

            u(:,iy,:,1) = step_cosine( abs(y), 0.10_rk * params_nspp%domain_size(2), &
            params_nspp%beta*params_nspp%domain_size(2) )
        end do
        
        ! uy is a sin wave (that triggers the instability)
        do ix = 1, Bs(1)+2*g
            ! compute x,y coordinates from spacing and origin
            x = dble(ix-(g+1)) * dx(1) + x0(1)

            u(ix,:,:,2) = 0.1_rk * sin(2.0_rk*pi*x/params_nspp%domain_size(1))
        end do


    case("velocity-blob")
        if (params_nspp%dim==2) then
            ! create gauss pulse. Note we loop over the entire block, incl. ghost nodes.
            do iy = 1, Bs(2)+2*g
                do ix = 1, Bs(1)+2*g
                    ! compute x,y coordinates from spacing and origin
                    x = dble(ix-(g+1)) * dx(1) + x0(1) - params_nspp%domain_size(1)/2.0_rk
                    y = dble(iy-(g+1)) * dx(2) + x0(2) - params_nspp%domain_size(2)/2.0_rk

                    if (x<-params_nspp%domain_size(1)/2.0) x = x + params_nspp%domain_size(1)
                    if (x>params_nspp%domain_size(1)/2.0) x = x - params_nspp%domain_size(1)

                    if (y<-params_nspp%domain_size(2)/2.0) y = y + params_nspp%domain_size(2)
                    if (y>params_nspp%domain_size(2)/2.0) y = y - params_nspp%domain_size(2)

                    ! set actual inicond gauss blob
                    ! here only for the pressure.
                    u(ix,iy,:,1:2) = dexp( -( (x)**2 + (y)**2 ) / params_nspp%beta )
                end do
            end do
        else
            ! create gauss pulse
            do iz = 1, Bs(3)+2*g
                do iy = 1, Bs(2)+2*g
                    do ix = 1, Bs(1)+2*g
                        ! compute x,y coordinates from spacing and origin
                        x = dble(ix-(g+1)) * dx(1) + x0(1) - params_nspp%domain_size(1)/2.0_rk
                        y = dble(iy-(g+1)) * dx(2) + x0(2) - params_nspp%domain_size(2)/2.0_rk
                        z = dble(iz-(g+1)) * dx(3) + x0(3) - params_nspp%domain_size(3)/2.0_rk

                        if (x<-params_nspp%domain_size(1)/2.0) x = x + params_nspp%domain_size(1)
                        if (x>params_nspp%domain_size(1)/2.0) x = x - params_nspp%domain_size(1)

                        if (y<-params_nspp%domain_size(2)/2.0) y = y + params_nspp%domain_size(2)
                        if (y>params_nspp%domain_size(2)/2.0) y = y - params_nspp%domain_size(2)

                        if (z<-params_nspp%domain_size(3)/2.0) z = z + params_nspp%domain_size(3)
                        if (z>params_nspp%domain_size(3)/2.0) z = z - params_nspp%domain_size(3)

                        ! set actual inicond gauss blob
                        u(ix,iy,iz,1:3) = dexp( -( (x)**2 + (y)**2 + (z)**2 ) / params_nspp%beta )
                    end do
                end do
            end do
        end if

    case("ux-blob")
        if (params_nspp%dim==2) then
            ! create gauss pulse. Note we loop over the entire block, incl. ghost nodes.
            do iy = 1, Bs(2)+2*g
                do ix = 1, Bs(1)+2*g
                    ! compute x,y coordinates from spacing and origin
                    x = dble(ix-(g+1)) * dx(1) + x0(1) - params_nspp%domain_size(1)/2.0_rk
                    y = dble(iy-(g+1)) * dx(2) + x0(2) - params_nspp%domain_size(2)/2.0_rk

                    if (x<-params_nspp%domain_size(1)/2.0) x = x + params_nspp%domain_size(1)
                    if (x>params_nspp%domain_size(1)/2.0) x = x - params_nspp%domain_size(1)

                    if (y<-params_nspp%domain_size(2)/2.0) y = y + params_nspp%domain_size(2)
                    if (y>params_nspp%domain_size(2)/2.0) y = y - params_nspp%domain_size(2)

                    ! set actual inicond gauss blob
                    ! here only for the pressure.
                    u(ix,iy,:,1) = dexp( -( (x)**2 + (y)**2 ) / params_nspp%beta )
                end do
            end do
        else
            ! create gauss pulse
            do iz = 1, Bs(3)+2*g
                do iy = 1, Bs(2)+2*g
                    do ix = 1, Bs(1)+2*g
                        ! compute x,y coordinates from spacing and origin
                        x = dble(ix-(g+1)) * dx(1) + x0(1) - params_nspp%domain_size(1)/2.0_rk
                        y = dble(iy-(g+1)) * dx(2) + x0(2) - params_nspp%domain_size(2)/2.0_rk
                        z = dble(iz-(g+1)) * dx(3) + x0(3) - params_nspp%domain_size(3)/2.0_rk

                        if (x<-params_nspp%domain_size(1)/2.0) x = x + params_nspp%domain_size(1)
                        if (x>params_nspp%domain_size(1)/2.0) x = x - params_nspp%domain_size(1)

                        if (y<-params_nspp%domain_size(2)/2.0) y = y + params_nspp%domain_size(2)
                        if (y>params_nspp%domain_size(2)/2.0) y = y - params_nspp%domain_size(2)

                        if (z<-params_nspp%domain_size(3)/2.0) z = z + params_nspp%domain_size(3)
                        if (z>params_nspp%domain_size(3)/2.0) z = z - params_nspp%domain_size(3)

                        ! set actual inicond gauss blob
                        u(ix,iy,iz,1) = dexp( -( (x)**2 + (y)**2 + (z)**2 ) / params_nspp%beta )
                    end do
                end do
            end do
        end if

    case("meanflow")
        do idir = 1, params_nspp%dim
            u(:,:,:,idir) = params_nspp%u_mean_set(idir)
        enddo
        
    case("meanflow_channel")
        ! set the initial condition for channel simulations: the ux-flow is 1 inside the channel and 0 inside the wall
        ! the y,z directions have perturbations (sin-waves) with ampltiude beta (set in INI-file)
        if (params_nspp%dim == 2) then
            call abort(2762025, "Initial condition meanflow_channel is valid only for 3D simulations")
        endif

        do iz = 1, Bs(3)+2*g
            do iy = 1, Bs(2)+2*g
                do ix = 1, Bs(1)+2*g
                    ! coordinates
                    x = x0(1) + dble(ix-g-1)*dx(1)
                    y = x0(2) + dble(iy-g-1)*dx(2)
                    z = x0(3) + dble(iz-g-1)*dx(3)

                    ! ux: constant 1 in the fluid, 0 in the solid
                    if (( y>params_nspp%h_channel).and.(y<params_nspp%domain_size(2)-params_nspp%h_channel)) then
                        ! fluid part
                        u(ix,iy,iz,1) = 1.0_rk
                    else
                        ! solid part
                        u(ix,iy,iz,1) = 0.0_rk
                    endif

                    ! uy, uz: sine wave perturbations
                    ! I don't know if this is the smartest choice ( I don't have a reference for it ) -TE
                    u(ix,iy,iz,2) = 0.0_rk + params_nspp%beta * sin( 2.0_rk*pi*x/params_nspp%domain_size(1) )
                    u(ix,iy,iz,3) = 0.0_rk + params_nspp%beta * sin( 2.0_rk*pi*z/params_nspp%domain_size(3) )
                    ! Note: pressure initialized to zero, will be computed via Poisson equation
                end do
            end do
        end do

    case("sinewaves-nopress")
        ! some random sine waves, but no pressure imposed.
        if (params_nspp%dim == 2) then
            do iy= 1,Bs(2)+2*g
                do ix= 1, Bs(1)+2*g
                    x = x0(1) + dble(ix-g-1)*dx(1)
                    y = x0(2) + dble(iy-g-1)*dx(2)
                    u(ix,iy,1,1) = sin( 2.0_rk*pi*x/params_nspp%domain_size(1) ) + 0.5_rk*sin( 10.0_rk*pi*x/params_nspp%domain_size(1) )
                    u(ix,iy,1,2) = cos( 2.0_rk*pi*y/params_nspp%domain_size(2) )*sin( 2.0_rk*pi*x/params_nspp%domain_size(1) )
                enddo
            enddo
        else
            do iz = 1, Bs(3)+2*g
                do iy= 1,Bs(2)+2*g
                    do ix= 1, Bs(1)+2*g
                        x = x0(1) + dble(ix-g-1)*dx(1)
                        y = x0(2) + dble(iy-g-1)*dx(2)
                        z = x0(3) + dble(iz-g-1)*dx(3)

                        u(ix,iy,iz,1) = sin( 2.0_rk*pi*x/params_nspp%domain_size(1) ) + 0.5_rk*sin( 10.0_rk*pi*x/params_nspp%domain_size(1) )
                        u(ix,iy,iz,2) = cos( 2.0_rk*pi*y/params_nspp%domain_size(2) )*sin( 2.0_rk*pi*x/params_nspp%domain_size(1) )*sin( 2.0_rk*pi*z/params_nspp%domain_size(3) )
                        u(ix,iy,iz,3) = cos( 2.0_rk*pi*y/params_nspp%domain_size(2) )*sin( 3.0_rk*pi*x/params_nspp%domain_size(1) )*cos( 2.0_rk*pi*z/params_nspp%domain_size(3) )
                    enddo
                enddo
            enddo

        endif

    case("taylor_green")
        ! this condition is 2D only!
        if (params_nspp%dim==2) then
            do iy= 1,Bs(2)+2*g
                y = x0(2) + dble(iy-g-1)*dx(2)
                call continue_periodic(y,params_nspp%domain_size(2))
                do ix= 1, Bs(1)+2*g
                    x = x0(1) + dble(ix-g-1)*dx(1)

                    call continue_periodic(x,params_nspp%domain_size(1))

                    u(ix,iy,1,1) = params_nspp%u_mean_set(1) - dsin(x)*dcos(y)
                    u(ix,iy,1,2) = params_nspp%u_mean_set(2) + dcos(x)*dsin(y)
                    ! Note: analytical pressure p = 0.25*(cos(2x) + cos(2y)) not initialized
                    ! Pressure will be computed from velocity field via Poisson equation
                end do
            end do
        else
            call abort(250708, "taylor_green is a 2D test case. Use taylor-green-vanRees2011 for 3D case")
        endif
    case ("taylor-green-vanRees2011")
        ! this condition is 3D only!
        if (params_nspp%dim==3) then
            do iz = 1, Bs(3)+2*g
                z = dble(iz-(g+1)) * dx(3) + x0(3)
                do iy = 1, Bs(2)+2*g
                    y = dble(iy-(g+1)) * dx(2) + x0(2)
                    do ix = 1, Bs(1)+2*g
                        x = dble(ix-(g+1)) * dx(1) + x0(1)

                        ! the initial condition is known analytically. Note in vanrees JCP2011,
                        ! there is a parameter \theta which is set to zero.
                        ! See also: "Problem C3.5 Direct Numerical Simulation of the Taylor-Green Vortex at Re = 1600"
                        ! NOTE: domain size has to be 2*pi = 6.283185307179586
                        u(ix, iy, iz, 1) = dsin(x)*dcos(y)*dcos(z)
                        u(ix, iy, iz, 2) =-dcos(x)*dsin(y)*dcos(z)
                        u(ix, iy, iz, 3) = 0.0_rk
                        ! Note: analytical pressure p = (cos(2x) + cos(2y))*(cos(2z) + 2) / 16 not initialized
                        ! Pressure will be computed from velocity field via Poisson equation
                    end do
                end do
            end do
        else
            call abort(250708, "taylor-green-vanRees2011 is a 3D test case. Use taylor-green for 2D case")
        endif
    case("mixing_layer")
        ! random excitement
        call random_data(u)
        u(:,:,:,params_nspp%dim+1:) = 0.0_rk  ! we only need random velocity data

        do iz = merge(1, g+1, params_nspp%dim==2), merge(1, Bs(3)+g, params_nspp%dim==2)
            z = 0.0_rk
            if (params_nspp%dim == 3) then
                z = dble(iz-(g+1)) * dx(3) + x0(3)
                call continue_periodic(z,params_nspp%domain_size(3))
                z = z - params_nspp%domain_size(3)/2.0_rk
            endif
            do iy = 1, Bs(2)+2*g
                y = dble(iy-(g+1)) * dx(2) + x0(2)
                call continue_periodic(y,params_nspp%domain_size(2))
                y = y - params_nspp%domain_size(2)/2.0_rk
                do ix = 1, Bs(1)+2*g
                    x = dble(ix-(g+1)) * dx(1) + x0(1)
                    call continue_periodic(x,params_nspp%domain_size(1))
                    x = x - params_nspp%domain_size(1)/2.0_rk

                    ! ! tanh profile with random variation to trigger 3D instabilities, change last number to set intensity
                    ! ! the exp term focusses the instabilities and fades it out towards the edges to ease the symmetry BC
                    ! u(ix,iy,iz,1) = tanh(z) * 1 + (u(ix,iy,iz,1)-0.5_rk)*2.0_rk*3e-1* exp(-z**2/5**2)
                    ! u(ix,iy,iz,2:params_nspp%dim) = u(ix,iy,iz,2:params_nspp%dim)*0e-6
                    
                    ! ! condition from Roussel & Schneider 2010 for 2D instabilities
                    ! param_1 = 3  ! maximum number of nodes
                    ! u(ix,iy,iz,1) = tanh(z) * 1 + 1/(2.0_rk*param_1*cosh(z)**2.0_rk) * &
                    !     (cos(2.0_rk*pi*(x+y)/params_nspp%domain_size(1)) + cos(2.0_rk*pi*(x-y)/params_nspp%domain_size(1)))
                    ! do idir = 1, param_1-1
                    !     u(ix,iy,iz,1) = u(ix,iy,iz,1) + 1/(2.0_rk*param_1*cosh(z)**2.0_rk) * &
                    !         cos(2**(idir+1)*pi*x/params_nspp%domain_size(1))
                    ! enddo
                    
                    ! condition to provoke 3D instabilities
                    ! amount of modes should change the second number in cosh denominator
                    ! every direction can have a random shift, which is the addition in the cos terms
                    ! for reproducability, these are hardcoded
                    param_1 = 2.0_rk*pi/dble(params_nspp%domain_size(1))  ! maximum number of nodes
                    param_2 = 2.7_rk
                    u(ix,iy,iz,1) = tanh(z) * 1 + 1/(2.0_rk*4.0_rk*cosh(z)**2.0_rk) * &
                        (cos(2.0_rk*(cos(param_2)*x+sin(param_2)*y)*param_1+4.0_rk) + cos(2.0_rk*(cos(param_2)*x-sin(-param_2)*y)*param_1+2.0_rk) + cos(2.0_rk*2.0_rk*pi*(z)/5.0_rk+1.0_rk))
                    param_2 = 1.1_rk
                    u(ix,iy,iz,1) = u(ix,iy,iz,1) * 1 + 1/(2.0_rk*4.0_rk*cosh(z)**2.0_rk) * &
                        (cos(3.0_rk*(cos(param_2)*x+sin(param_2)*y)*param_1+1.5_rk) + cos(3.0_rk*(cos(param_2)*x-sin(-param_2)*y)*param_1+3.7_rk) + cos(3.0_rk*2.0_rk*pi*(z)/5.0_rk+1.8_rk))
                    param_2 = 4.3_rk
                    u(ix,iy,iz,1) = u(ix,iy,iz,1) * 1 + 1/(2.0_rk*4.0_rk*cosh(z)**2.0_rk) * &
                        (cos(5.0_rk*(cos(param_2)*x+sin(param_2)*y)*param_1+5.2_rk) + cos(5.0_rk*(cos(param_2)*x-sin(-param_2)*y)*param_1+1.3_rk) + cos(5.0_rk*2.0_rk*pi*(z)/5.0_rk+3.0_rk))
                    param_2 = 0.5_rk
                    u(ix,iy,iz,1) = u(ix,iy,iz,1) * 1 + 1/(2.0_rk*4.0_rk*cosh(z)**2.0_rk) * &
                        (cos(7.0_rk*(cos(param_2)*x+sin(param_2)*y)*param_1+6.0_rk) + cos(7.0_rk*(cos(param_2)*x-sin(-param_2)*y)*param_1+1.4_rk) + cos(7.0_rk*2.0_rk*pi*(z)/5.0_rk+0.5_rk))

                    ! y- and z- velocity is not set
                    u(ix,iy,iz,2:params_nspp%dim) = u(ix,iy,iz,2:params_nspp%dim)*0e-6

                enddo
            enddo
        enddo

    case("mixing_layer_rollup")
        ! this IC amplifies the most unstable mode in 2D for a 2D roll-up (works in 2D and 3D), domain is symmetric in z
        do iz = merge(1, g+1, params_nspp%dim==2), merge(1, Bs(3)+g, params_nspp%dim==2)
            z = 0.0_rk
            if (params_nspp%dim == 3) then
                z = dble(iz-(g+1)) * dx(3) + x0(3)
                call continue_periodic(z,params_nspp%domain_size(3))
                z = z - params_nspp%domain_size(3)/2.0_rk
            endif
            do iy = 1, Bs(2)+2*g
                y = dble(iy-(g+1)) * dx(2) + x0(2)
                call continue_periodic(y,params_nspp%domain_size(2))
                y = y - params_nspp%domain_size(2)/2.0_rk
                do ix = 1, Bs(1)+2*g
                    x = dble(ix-(g+1)) * dx(1) + x0(1)
                    call continue_periodic(x,params_nspp%domain_size(1))
                    x = x - params_nspp%domain_size(1)/2.0_rk

                    ! tanh profile with specific thickness L/(4*7)
                    u(ix,iy,iz,1) = tanh(2.0_rk*y*28.0_rk/params_nspp%domain_size(2)) + 1.0_rk/(2.0_rk*cosh(y)**2.0_rk) * &
                        (sin(2.0_rk*pi*(x+z)/params_nspp%domain_size(1)) + sin(2.0_rk*pi*(x-z)/params_nspp%domain_size(1)))

                enddo
            enddo
        enddo
    case("mixing_layer_rollup2")
        ! this IC amplifies the most unstable mode in 2D for a 2D roll-up (works in 2D and 3D), domain is periodic in all directions with 2 shear layers
        do iz = merge(1, g+1, params_nspp%dim==2), merge(1, Bs(3)+g, params_nspp%dim==2)
            z = 0.0_rk
            if (params_nspp%dim == 3) then
                z = dble(iz-(g+1)) * dx(3) + x0(3)
                call continue_periodic(z,params_nspp%domain_size(3))
                z = z - params_nspp%domain_size(3)/2.0_rk
            endif
            do iy = 1, Bs(2)+2*g
                y = dble(iy-(g+1)) * dx(2) + x0(2)
                call continue_periodic(y,params_nspp%domain_size(2))
                y = y - params_nspp%domain_size(2)/2.0_rk
                do ix = 1, Bs(1)+2*g
                    x = dble(ix-(g+1)) * dx(1) + x0(1)
                    call continue_periodic(x,params_nspp%domain_size(1))
                    x = x - params_nspp%domain_size(1)/2.0_rk

                    ! tanh profile with specific thickness L/(4*7), but we have two of them to have periodicity
                    u(ix,iy,iz,1) = -1.0_rk - tanh(2.0_rk*(y - params_nspp%domain_size(2)/4.0_rk)*28.0_rk/params_nspp%domain_size(2)*2.0_rk) + &
                        1.0_rk/(2.0_rk*cosh(y - params_nspp%domain_size(2)/4.0_rk)**2.0_rk) * &
                        (sin(2.0_rk*pi*(x+z)/params_nspp%domain_size(1)*2.0_rk) + sin(2.0_rk*pi*(x-z)/params_nspp%domain_size(1)*2.0_rk)) &
                        + tanh(2.0_rk*(y + params_nspp%domain_size(2)/4.0_rk)*28.0_rk/params_nspp%domain_size(2)*2.0_rk) + &
                        1.0_rk/(2.0_rk*cosh(y + params_nspp%domain_size(2)/4.0_rk)**2.0_rk) * &
                        (sin(2.0_rk*pi*(x+z)/params_nspp%domain_size(1)*2.0_rk) + sin(2.0_rk*pi*(x-z)/params_nspp%domain_size(1)*2.0_rk))

                enddo
            enddo
        enddo
    case("mixing_layer_rollup3")
        ! this IC amplifies the most unstable mode in 2D for a 2D roll-up (works in 2D and 3D), domain is periodic in all directions with 2 shear layers
        ! compares to Yasuda2023 and AbdulGafoor2024
        do iz = merge(1, g+1, params_nspp%dim==2), merge(1, Bs(3)+g, params_nspp%dim==2)
            z = 0.0_rk
            if (params_nspp%dim == 3) then
                z = dble(iz-(g+1)) * dx(3) + x0(3)
                call continue_periodic(z,params_nspp%domain_size(3))
                z = z - params_nspp%domain_size(3)/2.0_rk
            endif
            do iy = 1, Bs(2)+2*g
                y = dble(iy-(g+1)) * dx(2) + x0(2)
                call continue_periodic(y,params_nspp%domain_size(2))
                y = y - params_nspp%domain_size(2)/2.0_rk
                do ix = 1, Bs(1)+2*g
                    x = dble(ix-(g+1)) * dx(1) + x0(1)
                    call continue_periodic(x,params_nspp%domain_size(1))
                    x = x - params_nspp%domain_size(1)/2.0_rk

                    ! tanh profile with specific thickness L/(4*7), but we have two of them to have periodicity
                    u(ix,iy,iz,1) = -1.0_rk - tanh((y - params_nspp%domain_size(2)/4.0_rk)*80.0_rk/params_nspp%domain_size(2)*2.0_rk) + &
                        tanh((y + params_nspp%domain_size(2)/4.0_rk)*80.0_rk/params_nspp%domain_size(2)*2.0_rk)    
                    u(ix,iy,iz,2) = 0.05_rk * sin((x - params_nspp%domain_size(2)/4.0_rk)*2.0_rk*pi)

                enddo
            enddo
        enddo

    case ("three-vortices")
        ! Three vortices test case - requires vorticity formulation
        if (params_nspp%dim /= 2) call abort(260202, "ERROR: 'three-vortices' is a 2D flow! set dim=2")
        if (.not. params_nspp%inicond_vorticity_formulation) then
            call abort(260202, "ERROR: 'three-vortices' requires inicond_vorticity_formulation=yes")
        endif

        ! Define vorticity field (will be transformed to velocity later)
        ! Vortex positions are relative to domain size
        do iy = 1, Bs(2)+2*g
            y = dble(iy-(g+1)) * dx(2) + x0(2)
            do ix = 1, Bs(1)+2*g
                x = dble(ix-(g+1)) * dx(1) + x0(1)

                ! Store vorticity in first component (will be transformed to velocity)
                u(ix,iy,:,:) = 0.0_rk

                ! First vortex: gamma = +1.0, center at (0.75/2*Lx, 1.0/2*Ly), width = 1/pi
                u(ix,iy,:,1) = (+1.0_rk/(pi*(1.0_rk/pi)**2)) * &
                        exp(-((x-(0.75_rk/2.0_rk)*params_nspp%domain_size(1))**2 + &
                             (y-(1.0_rk/2.0_rk)*params_nspp%domain_size(2))**2) / &
                            ((1.0_rk/pi)**2))

                ! Second vortex: gamma = +1.0, center at (1.25/2*Lx, 1.0/2*Ly), width = 1/pi
                u(ix,iy,:,1) = u(ix,iy,:,1) + (+1.0_rk/(pi*(1.0_rk/pi)**2)) * &
                        exp(-((x-(1.25_rk/2.0_rk)*params_nspp%domain_size(1))**2 + &
                             (y-(1.0_rk/2.0_rk)*params_nspp%domain_size(2))**2) / &
                            ((1.0_rk/pi)**2))

                ! Third vortex: gamma = -0.5, center at (1.25/2*Lx, (1.0/2 + 1/(4*sqrt(2)))*Ly), width = 1/pi
                u(ix,iy,:,1) = u(ix,iy,:,1) + (-0.5_rk/(pi*(1.0_rk/pi)**2)) * &
                        exp(-((x-(1.25_rk/2.0_rk)*params_nspp%domain_size(1))**2 + &
                             (y-((1.0_rk/2.0_rk + 1.0_rk/(4.0_rk*sqrt(2.0_rk)))*params_nspp%domain_size(2)))**2) / &
                            ((1.0_rk/pi)**2))

            enddo
        enddo

    case default
        call abort(428764, "NSPP inicond: "//trim(adjustl(params_nspp%inicond))//" is unkown.")

    end select

    ! --------------------------------------------------------------------------
    ! initial conditions for passive scalars, if used.
    ! --------------------------------------------------------------------------
    if (params_nspp%use_passive_scalar) then
        ! loop over scalars
        do iscalar = 1, params_nspp%N_scalars
            select case (params_nspp%scalar_inicond(iscalar))
            case ("empty", "none", "zero")
                u(:,:,:,params_nspp%dim + 1 + iscalar) = 0.0_rk
            case ("Kadoch2012")
                if (params_nspp%dim == 2) then
                    do iy = 1, Bs(2)+2*g
                        do ix = 1, Bs(1)+2*g
                            x = x0(1) + dble(ix-g-1)*dx(1) - 1.0_rk ! domain is -1...1 in kadoch
                            y = x0(2) + dble(iy-g-1)*dx(2) - 1.0_rk

                            x = x - params_nspp%length ! finite size of cavity
                            y = y - params_nspp%length

                            ! in their original work, they set the initial condition
                            ! everywhere in the domain (even in the penalization layer)
                            u(ix,iy,:,params_nspp%dim + 1 + iscalar) = cos(pi*y)*(cos(4.0_rk*pi*x)+cos(pi*x))
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

    ! initialize time statistics
    if (params_nspp%N_time_statistics > 0) then
        do iscalar = 1, params_nspp%N_time_statistics
            u(:,:,:, params_nspp%dim + 1 + params_nspp%N_scalars + iscalar) = 0.0_rk
        end do
    endif

end subroutine INICOND_NSPP
