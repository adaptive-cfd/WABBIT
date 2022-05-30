! ============================================================================================
!> The right hand side for the convection diffusion equations is implemented as follows:
!>\f{eqnarray*}{
!! \partial_t \phi &=& -u_0 \cdot \nabla \phi + \nu \nabla^2 \phi
!!\f}
! ********************************************************************************************


  !-----------------------------------------------------------------------------
  ! main level wrapper to set the right hand side on a block. Note this is completely
  ! independent of the grid and any MPI formalism, neighboring relations and the like.
  ! You just get a block data (e.g. ux, uy, uz, p) and compute the right hand side
  ! from that. Ghost nodes are assumed to be sync'ed.
  !-----------------------------------------------------------------------------
  subroutine RHS_convdiff( time, u, g, x0, dx, rhs, stage, boundary_flag  )
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

    ! output. Note assumed-shape arrays
    real(kind=rk), intent(inout) :: rhs(1:,1:,1:,1:)

    ! stage. there is 3 stages, init_stage, integral_stage and local_stage. If the PDE has
    ! terms that depend on global qtys, such as forces etc, which cannot be computed
    ! from a single block alone, the first stage does that. the second stage can then
    ! use these integral qtys for the actual RHS evaluation.
    character(len=*), intent(in) :: stage

    ! when implementing boundary conditions, it is necessary to know if the local field (block)
    ! is adjacent to a boundary, because the stencil has to be modified on the domain boundary.
    ! The boundary_flag tells you if the local field is adjacent to a domain boundary:
    ! boundary_flag(i) can be either 0, 1, -1,
    !  0: no boundary in the direction +/-e_i
    !  1: boundary in the direction +e_i
    ! -1: boundary in the direction - e_i
    ! currently only acessible in the local stage
    integer(kind=2), intent(in):: boundary_flag(3)

    ! local variables
    integer(kind=ik) :: mpierr
    integer(kind=ik), dimension(3) :: Bs
    real(kind=rk) :: tmp(1:3)

    ! compute the size of blocks
    Bs(1) = size(u,1) - 2*g
    Bs(2) = size(u,2) - 2*g
    Bs(3) = size(u,3) - 2*g

    select case(stage)
    case ("init_stage")
      !-------------------------------------------------------------------------
      ! 1st stage: init_stage.
      !-------------------------------------------------------------------------
      ! this stage is called only once, not for each block.
      ! performs initializations in the RHS module, such as resetting integrals

      return

    case ("integral_stage")
      !-------------------------------------------------------------------------
      ! 2nd stage: init_stage.
      !-------------------------------------------------------------------------
      ! For some RHS, the eqn depend not only on local, block based qtys, such as
      ! the state vector, but also on the entire grid, for example to compute a
      ! global forcing term (e.g. in FSI the forces on bodies). As the physics
      ! modules cannot see the grid, (they only see blocks), in order to encapsulate
      ! them nicer, two RHS stages have to be defined: integral / local stage.
      !if abs(params_convdiff%gamma)>0
      ! called for each block.

      return

    case ("post_stage")
      !-------------------------------------------------------------------------
      ! 3rd stage: post_stage.
      !-------------------------------------------------------------------------
      ! this stage is called only once, not for each block.

      return

    case ("local_stage")
      !-------------------------------------------------------------------------
      ! 4th stage: local evaluation of RHS on all blocks
      !-------------------------------------------------------------------------
      ! the second stage then is what you would usually do: evaluate local differential
      ! operators etc.
      !
      ! called for each block.
      call RHS_convdiff_new(time, g, Bs, dx, x0, u, rhs, boundary_flag)


    case default
      call abort(7771,"the RHS wrapper requests a stage this physics module cannot handle.")
    end select


  end subroutine RHS_convdiff




!> Actual implementation of the 2d/3d convdiff RHS
subroutine RHS_convdiff_new(time, g, Bs, dx, x0, phi, rhs, boundary_flag)

    implicit none

    !> time
    real(kind=rk), intent(in)                      :: time
    !> grid parameter
    integer(kind=ik), intent(in)                   :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    !> origin and spacing of the block
    real(kind=rk), intent(in)                      :: x0(:), dx(:)
    !> datafields
    real(kind=rk), intent(inout)                   :: phi(:,:,:,:)
    real(kind=rk), intent(inout)                   :: rhs(:,:,:,:)
    ! when implementing boundary conditions, it is necessary to know if the local field (block)
    ! is adjacent to a boundary, because the stencil has to be modified on the domain boundary.
    ! The boundary_flag tells you if the local field is adjacent to a domain boundary:
    ! boundary_flag(i) can be either 0, 1, -1,
    !  0: no boundary in the direction +/-e_i
    !  1: boundary in the direction +e_i
    ! -1: boundary in the direction - e_i
    ! currently only acessible in the local stage
    integer(kind=2), intent(in):: boundary_flag(3)

    ! real(kind=rk) :: u0(1:Bs(1)+2*g, 1:Bs(2)+2*g, 1:2)
    real(kind=rk) :: u0(1:Bs(1)+2*g, 1:Bs(2)+2*g, 1:Bs(3)+2*g, 1:3)
    real(kind=rk) :: dx_inv, dy_inv, dz_inv, dx2_inv, dy2_inv, dz2_inv, nu, gamma
    real(kind=rk) :: u_dx, u_dy, u_dz
    real(kind=rk) :: u_dxdx, u_dydy, u_dzdz
    real(kind=rk) :: xcb, ycb, beta, x,y, sech,bc
    ! loop variables
    integer(kind=ik) :: ix, iy, iz, i, N, ia1, ia2, ib1, ib2, ia, ib
    ! coefficients for Tam&Webb
    real(kind=rk)                                  :: a(-3:3)
    real(kind=rk)                                  :: b(-2:2)
    ! coefficients for a standard centered 4th order 1st derivative
    real(kind=rk), parameter :: a_FD4(-2:2) = (/1.0_rk/12.0_rk, -2.0_rk/3.0_rk, 0.0_rk, +2.0_rk/3.0_rk, -1.0_rk/12.0_rk/)


    ! set parameters for readability
    N = params_convdiff%N_scalars
    u0 = 0.0_rk
    rhs = 0.0_rk

    if (size(phi,1)/=Bs(1)+2*g .or. size(phi,2)/=Bs(2)+2*g .or. size(phi,4)/=N) then
        call abort(66233,"wrong size. The door rings, I'll leave u to it.")
    endif

    dx_inv = 1.0_rk / dx(1)
    dy_inv = 1.0_rk / dx(2)
    dx2_inv = 1.0_rk / (dx(1)**2)
    dy2_inv = 1.0_rk / (dx(2)**2)

    ! Tam & Webb, 4th order optimized (for first derivative)
    a = (/-0.02651995_rk, +0.18941314_rk, -0.79926643_rk, 0.0_rk, 0.79926643_rk, -0.18941314_rk, 0.02651995_rk/)
    ! Standard 4th order central scheme:
    !a = (/0.0_rk, 1.0_rk/12.0_rk, -2.0_rk/3.0_rk, 0.0_rk, 2.0_rk/3.0_rk, -1.0_rk/12.0_rk, 0.0_rk/)
    ! 4th order coefficients for second derivative
    b = (/ -1.0_rk/12.0_rk, 4.0_rk/3.0_rk, -5.0_rk/2.0_rk, 4.0_rk/3.0_rk, -1.0_rk/12.0_rk /)


    gamma = params_convdiff%gamma
    ! looop over components - they are independent scalars
    do i = 1, N

        ! because p%nu might load the entire params in the cache and thus be slower:
        nu = params_convdiff%nu(i)
        !!!!!!!!!!!!
        ! 2D
        !!!!!!!!!!!!
        if (params_convdiff%dim == 2) then
            ! create the advection velocity field, which may be time and space dependent
            u0 = 0.0_rk
            call create_velocity_field_2D( time, g, Bs, dx, x0, u0(:,:,1,1:2), i, phi(:,:,1,i) )

            select case(params_convdiff%discretization)
            case("FD_2nd_central")
                !-----------------------------------------------------------------------
                ! 2D, 2nd order
                !-----------------------------------------------------------------------
                ! please note the performance penalty associated with the if-clause so
                ! do not put it inside the loop, even if it's tempting.
                if (nu>=1.0e-10) then ! with viscosity
                    do iy = g+1, Bs(2)+g+ONE_SKIPREDUNDANT
                        do ix = g+1, Bs(1)+g+ONE_SKIPREDUNDANT
                            u_dx = (phi(ix+1,iy,1,i)-phi(ix-1,iy,1,i))*dx_inv*0.5_rk
                            u_dy = (phi(ix,iy+1,1,i)-phi(ix,iy-1,1,i))*dy_inv*0.5_rk

                            u_dxdx = (phi(ix-1,iy,1,i)-2.0_rk*phi(ix,iy,1,i)+phi(ix+1,iy,1,i))*dx2_inv
                            u_dydy = (phi(ix,iy-1,1,i)-2.0_rk*phi(ix,iy,1,i)+phi(ix,iy+1,1,i))*dy2_inv

                            rhs(ix,iy,1,i) = -u0(ix,iy,1,1)*u_dx -u0(ix,iy,1,2)*u_dy + nu*(u_dxdx+u_dydy)
                        end do
                    end do
                else !  no viscosity
                    do iy = g+1, Bs(2)+g+ONE_SKIPREDUNDANT
                        do ix = g+1, Bs(1)+g+ONE_SKIPREDUNDANT
                            u_dx = (phi(ix+1,iy,1,i)-phi(ix-1,iy,1,i))*dx_inv*0.5_rk
                            u_dy = (phi(ix,iy+1,1,i)-phi(ix,iy-1,1,i))*dy_inv*0.5_rk
                            rhs(ix,iy,1,i) = -u0(ix,iy,1,1)*u_dx -u0(ix,iy,1,2)*u_dy
                        end do
                    end do
                endif
            case("FD_4th_central_optimized")
                !-----------------------------------------------------------------------
                ! 2D, 4th order
                !-----------------------------------------------------------------------
                if (nu>=1.0e-10) then ! with viscosity
                    do iy = g+1, Bs(2)+g+ONE_SKIPREDUNDANT
                        do ix = g+1, Bs(1)+g+ONE_SKIPREDUNDANT
                            ! gradient
                            u_dx = (a(-3)*phi(ix-3,iy,1,i) + a(-2)*phi(ix-2,iy,1,i) + a(-1)*phi(ix-1,iy,1,i) + a(0)*phi(ix,iy,1,i)&
                            +  a(+3)*phi(ix+3,iy,1,i) + a(+2)*phi(ix+2,iy,1,i) + a(+1)*phi(ix+1,iy,1,i))*dx_inv
                            u_dy = (a(-3)*phi(ix,iy-3,1,i) + a(-2)*phi(ix,iy-2,1,i) + a(-1)*phi(ix,iy-1,1,i) + a(0)*phi(ix,iy,1,i)&
                            +  a(+3)*phi(ix,iy+3,1,i) + a(+2)*phi(ix,iy+2,1,i) + a(+1)*phi(ix,iy+1,1,i))*dy_inv

                            u_dxdx = (b(-2)*phi(ix-2,iy,1,1) + b(-1)*phi(ix-1,iy,1,1) + b(0)*phi(ix,iy,1,1)&
                            + b(+1)*phi(ix+1,iy,1,1) + b(+2)*phi(ix+2,iy,1,1))*dx2_inv
                            u_dydy = (b(-2)*phi(ix,iy-2,1,1) + b(-1)*phi(ix,iy-1,1,1) + b(0)*phi(ix,iy,1,1)&
                            + b(+1)*phi(ix,iy+1,1,1) + b(+2)*phi(ix,iy+2,1,1))*dy2_inv

                            rhs(ix,iy,1,i) = -u0(ix,iy,1,1)*u_dx -u0(ix,iy,1,2)*u_dy + nu*(u_dxdx+u_dydy)
                        end do
                    end do
                else ! no viscosity
                    do iy = g+1, Bs(2)+g+ONE_SKIPREDUNDANT
                        do ix = g+1, Bs(1)+g+ONE_SKIPREDUNDANT
                            ! gradient
                            u_dx = (a(-3)*phi(ix-3,iy,1,i) + a(-2)*phi(ix-2,iy,1,i) + a(-1)*phi(ix-1,iy,1,i) + a(0)*phi(ix,iy,1,i)&
                            +  a(+3)*phi(ix+3,iy,1,i) + a(+2)*phi(ix+2,iy,1,i) + a(+1)*phi(ix+1,iy,1,i))*dx_inv
                            u_dy = (a(-3)*phi(ix,iy-3,1,i) + a(-2)*phi(ix,iy-2,1,i) + a(-1)*phi(ix,iy-1,1,i) + a(0)*phi(ix,iy,1,i)&
                            +  a(+3)*phi(ix,iy+3,1,i) + a(+2)*phi(ix,iy+2,1,i) + a(+1)*phi(ix,iy+1,1,i))*dy_inv

                            rhs(ix,iy,1,i) = -u0(ix,iy,1,1)*u_dx -u0(ix,iy,1,2)*u_dy
                        end do
                    end do
                endif

            case("FD_4th_central")
                !-----------------------------------------------------------------------
                ! 2D, 4th order (standard scheme not the TW)
                !-----------------------------------------------------------------------
                do iy = g+1, Bs(2)+g+ONE_SKIPREDUNDANT
                    do ix = g+1, Bs(1)+g+ONE_SKIPREDUNDANT
                        ! gradient
                        u_dx = (a_FD4(-2)*phi(ix-2,iy,1,i) + a_FD4(-1)*phi(ix-1,iy,1,i) &
                             +  a_FD4(+2)*phi(ix+2,iy,1,i) + a_FD4(+1)*phi(ix+1,iy,1,i))*dx_inv
                        u_dy = (a_FD4(-2)*phi(ix,iy-2,1,i) + a_FD4(-1)*phi(ix,iy-1,1,i) &
                             +  a_FD4(+2)*phi(ix,iy+2,1,i) + a_FD4(+1)*phi(ix,iy+1,1,i))*dy_inv

                        u_dxdx = (b(-2)*phi(ix-2,iy,1,1) + b(-1)*phi(ix-1,iy,1,1) + b(0)*phi(ix,iy,1,1)&
                               +  b(+1)*phi(ix+1,iy,1,1) + b(+2)*phi(ix+2,iy,1,1))*dx2_inv
                        u_dydy = (b(-2)*phi(ix,iy-2,1,1) + b(-1)*phi(ix,iy-1,1,1) + b(0)*phi(ix,iy,1,1)&
                               +  b(+1)*phi(ix,iy+1,1,1) + b(+2)*phi(ix,iy+2,1,1))*dy2_inv

                        rhs(ix,iy,1,i) = -u0(ix,iy,1,1)*u_dx -u0(ix,iy,1,2)*u_dy + nu*(u_dxdx+u_dydy)
                    end do
                end do


            case default
                call abort(442161, params_convdiff%discretization//" discretization unkown, goto hell.")
            end select


            if (params_convdiff%velocity(i)=="cyclogenesis") then
                do iy = g+1, Bs(2)+g+ONE_SKIPREDUNDANT
                    do ix = g+1, Bs(1)+g+ONE_SKIPREDUNDANT
                        y = dble(iy-(g+1)) * dx(2) + x0(2) - params_convdiff%y0(i)
                        y = y / params_convdiff%blob_width(i)
                        ! sech(x) = 1.0 / cosh(x)
                        sech = 1.0_rk / cosh(y)
                        ! source term
                        rhs(ix,iy,1,i) = rhs(ix,iy,1,i) - u0(ix,iy,1,2)*(-sech**2 / params_convdiff%blob_width(i) )
                    end do
                end do
            endif
        !!!!!!!!!!!!!
        ! 3D
        !!!!!!!!!!!!!
        else
            dz_inv = 1.0_rk / dx(3)
            dz2_inv = 1.0_rk / (dx(3)**2)
            ! create the advection velocity field, which may be time and space dependent
            call create_velocity_field_3D( time, g, Bs, dx, x0, u0, i )

            select case(params_convdiff%discretization)
            case("FD_2nd_central")
                !-----------------------------------------------------------------------
                ! 3D, 2nd order
                !-----------------------------------------------------------------------
                if (nu>=1.0e-10) then ! with viscosity
                    do iz = g+1, Bs(3)+g+ONE_SKIPREDUNDANT
                        do iy = g+1, Bs(2)+g+ONE_SKIPREDUNDANT
                            do ix = g+1, Bs(1)+g+ONE_SKIPREDUNDANT
                                u_dx = (phi(ix+1,iy,iz,i) - phi(ix-1,iy,iz,i))*dx_inv*0.5_rk
                                u_dy = (phi(ix,iy+1,iz,i) - phi(ix,iy-1,iz,i))*dy_inv*0.5_rk
                                u_dz = (phi(ix,iy,iz+1,i) - phi(ix,iy,iz-1,i))*dz_inv*0.5_rk

                                u_dxdx = (phi(ix-1,iy,iz,i) -2.0_rk*phi(ix,iy,iz,i) + phi(ix+1,iy,iz,i))*dx2_inv
                                u_dydy = (phi(ix,iy-1,iz,i) -2.0_rk*phi(ix,iy,iz,i) + phi(ix,iy+1,iz,i))*dy2_inv
                                u_dzdz = (phi(ix,iy,iz-1,i) -2.0_rk*phi(ix,iy,iz,i) + phi(ix,iy,iz+1,i))*dz2_inv

                                rhs(ix,iy,iz,i) = -u0(ix,iy,iz,1)*u_dx -u0(ix,iy,iz,2)*u_dy -u0(ix,iy,iz,3)*u_dz &
                                + nu*(u_dxdx+u_dydy+u_dzdz)
                            end do
                        end do
                    end do
                else ! no viscosity
                    do iz = g+1, Bs(3)+g+ONE_SKIPREDUNDANT
                        do iy = g+1, Bs(2)+g+ONE_SKIPREDUNDANT
                            do ix = g+1, Bs(1)+g+ONE_SKIPREDUNDANT
                                u_dx = (phi(ix+1,iy,iz,i) - phi(ix-1,iy,iz,i))*dx_inv*0.5_rk
                                u_dy = (phi(ix,iy+1,iz,i) - phi(ix,iy-1,iz,i))*dy_inv*0.5_rk
                                u_dz = (phi(ix,iy,iz+1,i) - phi(ix,iy,iz-1,i))*dz_inv*0.5_rk

                                rhs(ix,iy,iz,i) = -u0(ix,iy,iz,1)*u_dx -u0(ix,iy,iz,2)*u_dy -u0(ix,iy,iz,3)*u_dz
                            end do
                        end do
                    end do
                endif

            case("FD_4th_central_optimized")
                !-----------------------------------------------------------------------
                ! 3D, 4th order
                !-----------------------------------------------------------------------
                if (nu>=1.0e-10) then ! with viscosity
                    do iz = g+1, Bs(3)+g+ONE_SKIPREDUNDANT
                        do iy = g+1, Bs(2)+g+ONE_SKIPREDUNDANT
                            do ix = g+1, Bs(1)+g+ONE_SKIPREDUNDANT
                                ! gradient
                                u_dx = (a(-3)*phi(ix-3,iy,iz,i) + a(-2)*phi(ix-2,iy,iz,i) + a(-1)*phi(ix-1,iy,iz,i) + a(0)*phi(ix,iy,iz,i)&
                                     +  a(+3)*phi(ix+3,iy,iz,i) + a(+2)*phi(ix+2,iy,iz,i) + a(+1)*phi(ix+1,iy,iz,i))*dx_inv
                                u_dy = (a(-3)*phi(ix,iy-3,iz,i) + a(-2)*phi(ix,iy-2,iz,i) + a(-1)*phi(ix,iy-1,iz,i) + a(0)*phi(ix,iy,iz,i)&
                                     +  a(+3)*phi(ix,iy+3,iz,i) + a(+2)*phi(ix,iy+2,iz,i) + a(+1)*phi(ix,iy+1,iz,i))*dy_inv
                                u_dz = (a(-3)*phi(ix,iy,iz-3,i) + a(-2)*phi(ix,iy,iz-2,i) + a(-1)*phi(ix,iy,iz-1,i) + a(0)*phi(ix,iy,iz,i)&
                                     +  a(+3)*phi(ix,iy,iz+3,i) + a(+2)*phi(ix,iy,iz+2,i) + a(+1)*phi(ix,iy,iz+1,i))*dz_inv

                                u_dxdx = (b(-2)*phi(ix-2,iy,iz,i) + b(-1)*phi(ix-1,iy,iz,i) + b(0)*phi(ix,iy,iz,i)&
                                       +  b(+1)*phi(ix+1,iy,iz,i) + b(+2)*phi(ix+2,iy,iz,i))*dx2_inv
                                u_dydy = (b(-2)*phi(ix,iy-2,iz,i) + b(-1)*phi(ix,iy-1,iz,i) + b(0)*phi(ix,iy,iz,i)&
                                       +  b(+1)*phi(ix,iy+1,iz,i) + b(+2)*phi(ix,iy+2,iz,i))*dy2_inv
                                u_dzdz = (b(-2)*phi(ix,iy,iz-2,i) + b(-1)*phi(ix,iy,iz-1,i) + b(0)*phi(ix,iy,iz,i)&
                                       +  b(+1)*phi(ix,iy,iz+1,i) + b(+2)*phi(ix,iy,iz+2,i))*dy2_inv

                                rhs(ix,iy,iz,i) = -u0(ix,iy,iz,1)*u_dx -u0(ix,iy,iz,2)*u_dy -u0(ix,iy,iz,3)*u_dz + nu*(u_dxdx+u_dydy+u_dzdz)

                            end do
                        end do
                    end do
                else ! no viscosity
                    do iz = g+1, Bs(3)+g+ONE_SKIPREDUNDANT
                        do iy = g+1, Bs(2)+g+ONE_SKIPREDUNDANT
                            do ix = g+1, Bs(1)+g+ONE_SKIPREDUNDANT
                                ! gradient
                                u_dx = (a(-3)*phi(ix-3,iy,iz,i) + a(-2)*phi(ix-2,iy,iz,i) + a(-1)*phi(ix-1,iy,iz,i) + a(0)*phi(ix,iy,iz,i)&
                                     +  a(+3)*phi(ix+3,iy,iz,i) + a(+2)*phi(ix+2,iy,iz,i) + a(+1)*phi(ix+1,iy,iz,i))*dx_inv
                                u_dy = (a(-3)*phi(ix,iy-3,iz,i) + a(-2)*phi(ix,iy-2,iz,i) + a(-1)*phi(ix,iy-1,iz,i) + a(0)*phi(ix,iy,iz,i)&
                                     +  a(+3)*phi(ix,iy+3,iz,i) + a(+2)*phi(ix,iy+2,iz,i) + a(+1)*phi(ix,iy+1,iz,i))*dy_inv
                                u_dz = (a(-3)*phi(ix,iy,iz-3,i) + a(-2)*phi(ix,iy,iz-2,i) + a(-1)*phi(ix,iy,iz-1,i) + a(0)*phi(ix,iy,iz,i)&
                                     +  a(+3)*phi(ix,iy,iz+3,i) + a(+2)*phi(ix,iy,iz+2,i) + a(+1)*phi(ix,iy,iz+1,i))*dz_inv

                                rhs(ix,iy,iz,i) = -u0(ix,iy,iz,1)*u_dx -u0(ix,iy,iz,2)*u_dy -u0(ix,iy,iz,3)*u_dz
                            end do
                        end do
                    end do
                endif

            case default
                call abort(442161, params_convdiff%discretization//" discretization unkown, goto hell.")
            end select
        endif

        if (maxval(phi(:,:,:,i))>2000000.0_rk) call abort(666,"large values in phi. that cannot be good.")
    end do

    ! reaction
    if (gamma>=1.0e-10) then
      rhs = rhs - gamma * (phi-1) * phi**2
    endif

end subroutine RHS_convdiff_new



subroutine create_velocity_field_2D( time, g, Bs, dx, x0, u0, i, u )
    implicit none
    real(kind=rk), intent(in) :: time
    integer(kind=ik), intent(in) :: g, i
    integer(kind=ik), dimension(3), intent(in) :: Bs
    real(kind=rk), intent(in) :: dx(1:2), x0(1:2)
    real(kind=rk), intent(inout) :: u0(:,:,:)
    real(kind=rk), intent(in) :: u(:,:)
    ! note you cannot change these values without recomputing the coefficients
    real(kind=rk), parameter :: tau= 0.30_rk, t0=0.0_rk, t1=0.55_rk, t2=1.0_rk, u1=1.0_rk, u2=-1.2221975311385904
    real(kind=rk) :: u_this

    integer(kind=ik) :: ix,iy
    real(kind=rk) :: x,y,c0x,c0y, T, c0,c1,c2,c3, phi, r, ut
    real(kind=rk) :: x1,y1,x2,y2,r1,r2,r0,w0(2)

    u0 = 0.0_rk
    c0x = 0.5_rk*params_convdiff%domain_size(1)
    c0y = 0.5_rk*params_convdiff%domain_size(2)
    T = params_convdiff%T_swirl


    select case(params_convdiff%velocity(i))
    case ("nonlinear")
        u0(:,:,1) = u+params_convdiff%u_const

    case ("cyclogenesis")
        do iy = 1, Bs(2) + 2*g
            do ix = 1, Bs(1) + 2*g
                x = dble(ix-(g+1)) * dx(1) + x0(1)
                y = dble(iy-(g+1)) * dx(2) + x0(2)
                ! radius
                r = dsqrt( (x-params_convdiff%x0(i))**2 + (y-params_convdiff%y0(i))**2 )
                ! tangential velocity
                ut = 1.0d0 * (1.d0 / (cosh(r))) * (tanh(r))
                ! angle phi
                phi = atan2(y-params_convdiff%y0(i),x-params_convdiff%x0(i))

                u0(ix,iy,1) = -sin(phi) * ut
                u0(ix,iy,2) =  cos(phi) * ut
            enddo
        enddo

    case ("swirl")

        do iy = 1, Bs(2) + 2*g
            do ix = 1, Bs(1) + 2*g
                x = dble(ix-(g+1)) * dx(1) + x0(1)
                y = dble(iy-(g+1)) * dx(2) + x0(2)

                u0(ix,iy,1) = cos((pi*time)/T) * (sin(pi*x))**2 * sin(2*pi*y)
                u0(ix,iy,2) = cos((pi*time)/T) * (sin(pi*y))**2 * (-sin(2*pi*x))
            enddo
        enddo

    case("constant")
        u0(:,:,1) = params_convdiff%u0x(i)
        u0(:,:,2) = params_convdiff%u0y(i)

    case("circular")
          u0(:,:,1) = -params_convdiff%domain_size(1) * pi * 0.5 * sin(2*pi*time)
          u0(:,:,2) = params_convdiff%domain_size(2) * pi * 0.5 * cos(2*pi*time)

    case("vortex-pair")
      x1 = 0.6 -  params_convdiff%u0x(i) * time
      y1 = 0.49 * params_convdiff%domain_size(2)
      x2 = 0.6 -  params_convdiff%u0x(i) * time
      y2 = 0.51 * params_convdiff%domain_size(2)
      w0(1:2) = params_convdiff%w0(1:params_convdiff%dim)
      r0 = 0.005   ! initial size of vortex
      do iy = 1, Bs(2) + 2*g
          do ix = 1, Bs(1) + 2*g
              x = dble(ix-(g+1)) * dx(1) + x0(1)
              y = dble(iy-(g+1)) * dx(2) + x0(2)

              r1 = (x-x1)**2 + (y-y1)**2
              r2 = (x-x2)**2 + (y-y2)**2

              u0(ix,iy,1) = -w0(1) * (y - y1) * dexp(-(r1/r0)**2) + w0(2) * (y - y2) * dexp(-(r2/r0)**2)
              u0(ix,iy,2) =  w0(1) * (x - x1) * dexp(-(r1/r0)**2) - w0(2) * (x - x2) * dexp(-(r2/r0)**2)
          enddo
      enddo
      u0 = u0 * dexp(-(time / params_convdiff%tau) ** 2)

    case default
        call abort(77262,params_convdiff%velocity(i)//' is an unkown velocity field. It is time to go home.')
    end select
end subroutine



subroutine create_velocity_field_3D( time, g, Bs, dx, x0, u0, i )
    implicit none
    real(kind=rk), intent(in) :: time
    integer(kind=ik), intent(in) :: g, i
    integer(kind=ik), dimension(3), intent(in) :: Bs
    real(kind=rk), intent(in) :: dx(1:3), x0(1:3)
    real(kind=rk), intent(inout) :: u0(:,:,:,1:)

    integer(kind=ik) :: ix, iy, iz
    real(kind=rk) :: x, y, z, c0x, c0y, c0z, T

    select case(params_convdiff%velocity(i))
    case ("swirl-helix")
        do iz = 1, Bs(3) + 2*g
            z = dble(iz-(g+1)) * dx(3) + x0(3)
            do iy = 1, Bs(2) + 2*g
                y = dble(iy-(g+1)) * dx(2) + x0(2)
                do ix = 1, Bs(1) + 2*g
                    x = dble(ix-(g+1)) * dx(1) + x0(1)

                    u0(ix,iy,iz,1) = cos((pi*time)/T) * (sin(pi*x))**2 * sin(2*pi*y)
                    u0(ix,iy,iz,2) = cos((pi*time)/T) * (sin(pi*y))**2 * (-sin(2*pi*x))
                    u0(ix,iy,iz,3) = cos((pi*time)/T)
                end do
            end do
        end do

    case("constant")
        u0(:,:,:,1) = params_convdiff%u0x(i)
        u0(:,:,:,2) = params_convdiff%u0y(i)
        u0(:,:,:,3) = params_convdiff%u0z(i)

    case default
        call abort(77262,params_convdiff%velocity(i)//' is an unkown velocity field. It is time to go home.')
    end select
end subroutine
