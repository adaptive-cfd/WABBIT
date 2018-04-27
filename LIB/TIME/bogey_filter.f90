!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name bogey_filter.f90
!> \version 0.5
!> \author msr
!
!> \brief bogey shock filter subroutine
!
!>
!! input:    - params, heavy data,  \n
!! output:   - heavy data  \n
!!
!!
!! = log ======================================================================================
!! \n
!! 21/09/17 - create
! ********************************************************************************************

subroutine bogey_filter( params,Bs,g,N_dF, block_data,xx0,ddx,hvy_WORK)

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: block_data(:, :, :, :)
    !> heavy work    
    real(kind=rk), intent(inout)        :: hvy_work(:, :, :, :)
    ! spacing and origin of a block
    real(kind=rk), intent(in)           :: xx0(1:3), ddx(1:3)


    ! loop parameter
    integer(kind=ik)                    :: i, j, l, dF
!    ! grid parameter
!    integer(kind=ik)                    :: Bs, g, N_dF

    ! filtered values and array for old block data
    real(kind=rk)                       :: phi_tilde(3), r_xyz(3), r_th, eps, c1, c2, c_stencil(4)
!    real(kind=rk), allocatable          :: block_old(:, :, :, :), r_x(:,:), r_y(:,:), sigma_x(:,:), sigma_y(:,:), theta(:,:), u_x(:,:), v_y(:,:), Dtheta_x(:,:), &
!                                           Dtheta_y(:,:), Dthetamag_x(:,:), Dthetamag_y(:,:), c_2(:,:), rho(:,:), u(:,:), v(:,:), p(:,:)
    real(kind=rk)          :: block_old(Bs+2*g, Bs+2*g, Bs+2*g, N_dF), r_x(Bs+2*g, Bs+2*g), r_y(Bs+2*g, Bs+2*g), sigma_x(Bs+2*g, Bs+2*g), sigma_y(Bs+2*g, Bs+2*g), theta(Bs+2*g, Bs+2*g), &
              u_x(Bs+2*g, Bs+2*g),v_y(Bs+2*g, Bs+2*g), Dtheta_x(Bs+2*g, Bs+2*g), Dtheta_y(Bs+2*g, Bs+2*g), Dthetamag_x(Bs+2*g, Bs+2*g), &
              Dthetamag_y(Bs+2*g, Bs+2*g), c_2(Bs+2*g, Bs+2*g), rho(Bs+2*g, Bs+2*g), u(Bs+2*g, Bs+2*g), v(Bs+2*g, Bs+2*g), p(Bs+2*g, Bs+2*g)

    ! stencil array, note: size is fixed
    real(kind=rk)                       :: stencil(3),gamma_,r_rs(2)

    ! filter position (array postion of value to filter)
    integer(kind=ik)                    :: stencil_size

    ! /todo: move to ini file
    character(len=80)                   :: detector_method, sigma_method
    !> grid parameter
    integer(kind=ik)                    :: g, Bs,N_dF

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization
!
    ! grid parameter
    ! Bs  = params%number_block_nodes
    ! g   = params%number_ghost_nodes
    ! N_dF= params%number_data_fields

    ! allocate old data array
!    allocate( block_old(Bs+2*g, Bs+2*g, Bs+2*g, N_dF), r_x(Bs+2*g, Bs+2*g), r_y(Bs+2*g, Bs+2*g), sigma_x(Bs+2*g, Bs+2*g), sigma_y(Bs+2*g, Bs+2*g), theta(Bs+2*g, Bs+2*g), &
!              u_x(Bs+2*g, Bs+2*g),v_y(Bs+2*g, Bs+2*g), Dtheta_x(Bs+2*g, Bs+2*g), Dtheta_y(Bs+2*g, Bs+2*g), Dthetamag_x(Bs+2*g, Bs+2*g), &
!              Dthetamag_y(Bs+2*g, Bs+2*g), c_2(Bs+2*g, Bs+2*g), rho(Bs+2*g, Bs+2*g), u(Bs+2*g, Bs+2*g), v(Bs+2*g, Bs+2*g), p(Bs+2*g, Bs+2*g) )
    sigma_x=0.0_rk
    sigma_y=0.0_rk

    stencil_size = 3
    stencil(1:stencil_size) = (/  1.0_rk/  4.0_rk, &
                                 -1.0_rk/  2.0_rk, &
                                  1.0_rk/  4.0_rk/)

    ! threshold value
    r_th = params%r_th!/( real(level,kind=rk) )**2.0_rk

    detector_method = params%detector_method ! "p", "divU"
    sigma_method    = params%sigma_switch ! "abs", "tanh"

    eps = 1e-16

    ! second order
    ! c1 = -1/4;
    ! c2 = 0;
    ! fourth order
    ! c1 = -3/16;
    ! c2 = 1/16;
    ! fourth order - opt
    c1 = -0.210383_rk;
    c2 = 0.030617_rk;

    c_stencil = (/ -c2, -c1, c1, c2 /);

    ! hack: reset boundary nodes
 !   call set_boundary( Bs, g, block_data(:, :, 1, 1), block_data(:, :, 1, 2), block_data(:, :, 1, 3), block_data(:, :, 1, 4), xx0, ddx, params%Lx, params%Ly )
    !> \todo make filter interface for different physics modules
    if (params%physics_type=="navier_stokes") then
        rho = block_data(:, :, 1, 1)**2.0_rk
        u   = block_data(:, :, 1, 2)/block_data(:, :, 1, 1)
        v   = block_data(:, :, 1, 3)/block_data(:, :, 1, 1)
        p   = block_data(:, :, 1, 4)
        gamma_=1.4
    else
        call abort(13463,'Error [bogey_filter]: only ns supported: everything else not implemented yet!')
    endif
    
!---------------------------------------------------------------------------------------------
! main body

    ! first - shock detector
    ! detector input
    select case (detector_method)
        case("p")
            theta = block_data(:, :, 1, 4)

        case("divU")
            u_x = 0.0_rk
            v_y = 0.0_rk

            do i = g-2,Bs+g+3
                do j = g-2,Bs+g+3

                    !u_x(i,j) = ( u(i-2,j) - 8.0_rk*u(i-1,j) + 8.0_rk*u(i+1,j) - u(i+2,j) ) / (12.0_rk*ddx(1))
                    !v_y(i,j) = ( v(i,j-2) - 8.0_rk*v(i,j-1) + 8.0_rk*v(i,j+1) - v(i,j+2) ) / (12.0_rk*ddx(2))

                    u_x(i,j) = ( u(i+1,j) - u(i-1,j) ) / (2.0_rk*ddx(1))
                    v_y(i,j) = ( v(i,j+1) - v(i,j-1) ) / (2.0_rk*ddx(2))

                end do
            end do

            theta = u_x + v_y

    end select


    ! filter detector values, note: also filter first two ghost nodes
    do i = g-1,Bs+g+2
        do j = g-1,Bs+g+2
            call filter_1D( theta(i-1:i+1, j), Dtheta_x(i,j), stencil(1:stencil_size) )
            call filter_1D( theta(i, j-1:j+1), Dtheta_y(i,j), stencil(1:stencil_size) )
        end do
    end do

    ! magnitude of filtered values
    do i = g,Bs+g+1
        do j = g,Bs+g+1
            Dthetamag_x(i,j) = 0.5_rk * ( ( Dtheta_x(i,j) - Dtheta_x(i+1,j) )**2.0_rk + ( Dtheta_x(i,j) - Dtheta_x(i-1,j) )**2.0_rk )
            Dthetamag_y(i,j) = 0.5_rk * ( ( Dtheta_y(i,j) - Dtheta_y(i,j+1) )**2.0_rk + ( Dtheta_y(i,j) - Dtheta_y(i,j-1) )**2.0_rk )
        end do
    end do

    ! second - shock sensor
    select case(detector_method)
        case("p")
            do i = g,Bs+g+1
                do j = g,Bs+g+1
                    r_x(i,j) = Dthetamag_x(i,j) / p(i,j)**2.0_rk + eps
                    r_y(i,j) = Dthetamag_y(i,j) / p(i,j)**2.0_rk + eps
                end do
            end do

        case("divU")
            !> \todo write physics method interface for speed of sound
            c_2 = gamma_*p/rho
            do i = g,Bs+g+1
                do j = g,Bs+g+1
                    r_x(i,j) = Dthetamag_x(i,j) / (c_2(i,j)/ddx(1)**2.0_rk) + eps
                    r_y(i,j) = Dthetamag_y(i,j) / (c_2(i,j)/ddx(2)**2.0_rk) + eps
                end do
            end do

    end select

    ! third - filter value magnitude
    select case(sigma_method)
        case("abs")
            do i = g,Bs+g+1
                do j = g,Bs+g+1
                    sigma_x(i,j) = 0.5_rk * ( 1.0_rk - r_th/r_x(i,j) + abs( 1.0_rk - r_th/r_x(i,j) ) )
                    sigma_y(i,j) = 0.5_rk * ( 1.0_rk - r_th/r_y(i,j) + abs( 1.0_rk - r_th/r_y(i,j) ) )
                end do
            end do

        case("tanh")
            do i = g,Bs+g+1
                do j = g,Bs+g+1
                     sigma_x(i,j) = 1.0_rk - dtanh( r_th/r_x(i,j)/0.7_rk )
                     sigma_y(i,j) = 1.0_rk - dtanh( r_th/r_y(i,j)/0.7_rk )
                    
                     if (sigma_x(i,j)>1.0_rk .or. sigma_x(i,j)<0.0_rk) then
                        call abort(343674,'something wrong')
                     endif
                    ! r_rs(1)     = r_th/(r_x(i,j)-0.5_rk*r_th)
                    ! r_rs(2)     = r_th/(r_y(i,j)-0.5_rk*r_th)
                    ! if (r_rs(1)<2.0_rk) then
                    !     sigma_x(i,j) = 1.0_rk - dtanh( r_rs(1) )
                    ! else
                    !     sigma_x(i,j) = 0.0_rk
                    ! endif

                    ! if (r_rs(1)<2.0_rk) then
                    !     sigma_y(i,j) = 1.0_rk - dtanh( r_rs(2) )
                    ! else
                    !     sigma_y(i,j)=0.0_rk
                    ! endif
                    
                end do
            end do

    end select

    ! fourth - filter fields
    block_old(:,:,1,:) = block_data(:,:,1,:)
    do dF = 1, N_dF
        do i = g+1,Bs+g
            do j = g+1,Bs+g

                block_data(i,j,1,dF) = block_old(i,j,1,dF)

                block_data(i,j,1,dF) = block_data(i,j,1,dF) - ( 0.5_rk * (sigma_x(i+1,j) + sigma_x(i,j)) &
                                                                       * sum( c_stencil * block_old(i-1:i+2,j,1,dF) ) &
                                                            -   0.5_rk * (sigma_x(i-1,j) + sigma_x(i,j)) &
                                                                       * sum( c_stencil * block_old(i-2:i+1,j,1,dF) ) )

                block_data(i,j,1,dF) = block_data(i,j,1,dF) - ( 0.5_rk * (sigma_y(i,j+1) + sigma_y(i,j)) &
                                                                       * sum( c_stencil * block_old(i,j-1:j+2,1,dF) ) &
                                                            -   0.5_rk * (sigma_y(i,j-1) + sigma_y(i,j)) &
                                                                       * sum( c_stencil * block_old(i,j-2:j+1,1,dF) ) )

            end do
        end do
    end do

    if (params%save_filter_strength) then
        ! save filter strength in x direction
        hvy_WORK(:,:,1,1)=sigma_x(:,:)
        ! save filter strength in y direction
        hvy_WORK(:,:,1,2)=sigma_y(:,:)
    endif

end subroutine bogey_filter

! shock detector - 1D
! -------------------
subroutine shock_detector_1D( phi, phi_old, r )

    ! global parameters
    use module_params

    implicit none

    !> datafield
    real(kind=rk), intent(in)           :: phi(:)
    !> old value
    real(kind=rk), intent(in)           :: phi_old
    !> r value
    real(kind=rk), intent(out)          :: r

    r = 0.5_rk * ( ( phi(2) - phi(3) )**2 + ( phi(2) - phi(1) )**2 ) / ( ( phi_old )**2 + 1.0_rk) + 1.0e-16_rk

end subroutine shock_detector_1D

! -------------------

subroutine set_boundary( Bs, g, rho, u, v, p, xx0, ddx, Lx, Ly )
    use module_params
    integer(kind=ik), intent(in)         :: Bs, g
    real(kind=rk), intent(inout)         :: rho(Bs+2*g, Bs+2*g), u(Bs+2*g, Bs+2*g), v(Bs+2*g, Bs+2*g), p(Bs+2*g, Bs+2*g)
    real(kind=rk), intent(in)            :: xx0(1:2), ddx(1:2), Lx, Ly

    integer(kind=ik)                     :: i

    ! north
    if ( abs( xx0(2) + ddx(2)*real(Bs-1,kind=rk)  - Ly ) < 1e-12_rk ) then
        ! nothing to do
    end if

    ! east
    if ( abs( xx0(1) + ddx(1)*real(Bs-1,kind=rk)  - Lx ) < 1e-12_rk ) then
        ! set boundary
        do i = Bs+g+1, Bs+2*g
            rho(i, :) = rho(Bs+g, :)
            u(i, :)   = u(Bs+g, :)
            v(i, :)   = v(Bs+g, :)
            p(i, :)   = p(Bs+g, :)
        end do
    end if

    ! south
    if ( abs( xx0(2) - 0.0_rk ) < 1e-12_rk ) then
        ! nothing to do
    end if

    ! west
    if ( abs( xx0(1) - 0.0_rk ) < 1e-12_rk ) then
        ! set boundary
        do i = 1, g
            rho(i, :) = rho(g+1, :)
            u(i, :)   = u(g+1, :)
            v(i, :)   = v(g+1, :)
            p(i, :)   = p(g+1, :)
        end do
    end if

end subroutine set_boundary
