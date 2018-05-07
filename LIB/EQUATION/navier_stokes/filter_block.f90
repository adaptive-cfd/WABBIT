!---------------------------------------------------!!!!
!> \file This file includes all filter routines
!---------------------------------------------------!!!!



    ! old (not-filterd) statevector
    ! -----------------------------
    ! + filter values are added to the original statevector
    ! + is initialiced in init_filter
    ! real(kind=rk), allocatable          :: block_old(:, :, :)


!> \brief this function must be called before filter_block!

subroutine init_filter(filter, FILE )
    implicit none
    !> pointer to inifile
    type(inifile)           ,intent(inout)  :: FILE
    !> params structure of navier stokes
    type(type_params_filter),intent(inout)  :: filter

    ! stencil array, note: size is fixed
    real(kind=rk)                      :: stencil(19)
    ! filter position (array postion of value to filter)
    integer(kind=ik)                   :: stencil_size
    ! number of number_ghost_nodes
    integer(kind=ik)                   :: g,mpi_rank,mpicode


    call MPI_COMM_RANK (WABBIT_COMM, mpi_rank, mpicode)
    if (mpi_rank==0) then
      write(*,*)
      write(*,*)
      write(*,*) "PARAMS: filter"
      write(*,'(" ---------------")')
    endif
    ! filter type
    call read_param_mpi(FILE, 'Discretization', 'filter_type', filter%name, "no-filter" )
    ! reset stencil_size
    stencil_size = 0

        ! set filter pos
        select case(filter%name)
            case('explicit_5pt')
                stencil_size = 5
                stencil(1:stencil_size) = (/ -1.0_rk/ 16.0_rk, &
                                              1.0_rk/  4.0_rk, &
                                             -3.0_rk/  8.0_rk, &
                                              1.0_rk/  4.0_rk, &
                                             -1.0_rk/ 16.0_rk/)

            case('explicit_7pt')
                stencil_size = 7
                stencil(1:stencil_size) = (/  1.0_rk/ 64.0_rk, &
                                             -3.0_rk/ 32.0_rk, &
                                             15.0_rk/ 64.0_rk, &
                                             -5.0_rk/ 16.0_rk, &
                                             15.0_rk/ 64.0_rk, &
                                             -3.0_rk/ 32.0_rk, &
                                              1.0_rk/ 64.0_rk/)

            case('explicit_9pt')
                stencil_size = 9
                stencil(1:stencil_size) = (/ -1.0_rk/256.0_rk, &
                                              1.0_rk/ 32.0_rk, &
                                             -7.0_rk/ 64.0_rk, &
                                              7.0_rk/ 32.0_rk, &
                                            -35.0_rk/128.0_rk, &
                                              7.0_rk/ 32.0_rk, &
                                             -7.0_rk/ 64.0_rk, &
                                              1.0_rk/ 32.0_rk, &
                                             -1.0_rk/256.0_rk/)

            case('explicit_11pt')
                stencil_size = 11
                stencil(1:stencil_size) = (/  1.0_rk/1024.0_rk, &
                                             -5.0_rk/ 512.0_rk, &
                                             45.0_rk/1024.0_rk, &
                                            -15.0_rk/ 128.0_rk, &
                                            105.0_rk/ 512.0_rk, &
                                            -63.0_rk/ 256.0_rk, &
                                            105.0_rk/ 512.0_rk, &
                                            -15.0_rk/ 128.0_rk, &
                                             45.0_rk/1024.0_rk, &
                                             -5.0_rk/ 512.0_rk, &
                                              1.0_rk/1024.0_rk/)

            case('no_filter')
                ! do nothing..

            case('wavelet')
                ! order of predictor for refinement
                call read_param_mpi(FILE, 'Discretization', 'order_predictor', filter%order_predictor, "---" )
                ! read threshold value
                call read_param_mpi(FILE, 'Blocks', 'eps', filter%eps, 1e-3_rk )

            case('bogey_shock')
                ! bogey shock detector threshold
                call read_param_mpi(FILE, 'Discretization', 'r_th', filter%r_th, 1e-3_rk )
                ! bogey shock switch tanh
                call read_param_mpi(FILE, 'Discretization', 'switch', filter%sigma_switch, 'tanh' )
                ! bogey shock detection method (p,divU)
                call read_param_mpi(FILE, 'Discretization', 'detector_method', filter%detector_method, 'divU' )
                ! boolean save bogey filter strength
                call read_param_mpi(FILE, 'Discretization', 'save_filter_strength', filter%save_filter_strength, .false. )

            case default
                call abort(4564,"ERROR [filter_block.f90]: filter type is not known!")
        end select

        call read_param_mpi(FILE, 'Blocks', 'number_data_fields', filter%number_data_fields, 1 )
        ! check ghost nodes number
        ! read number_ghost_nodes
        call read_param_mpi(FILE, 'Blocks', 'number_ghost_nodes', g, 1 )
        if ( g < ((stencil_size+1)/2-1) ) then
            call abort(7964,"ERROR: number of ghost nodes is too low for filtering")
        end if

        filter%stencil      =stencil
        filter%stencil_size =stencil_size

end subroutine init_filter

























!> \brief Filter statevector \f$U\f$ in conservative variables on each block
!
!> \details
!> \f{eqnarray*}{
!!   U\mapsto\tilde{U}={\mathrm F}(U)
!!\f}
!! Currently 3 different filters F(U) are available:
!!          + explicit \f$n\f$-point filters (\f$n=5,9,11\f$ )
!>              \f{eqnarray*}{
!!                  F(U)_i  &=U_i+(DU_i)\\
!!                  DU_i    &=\sum\limits_{j=1-n}^n d_j U_{i+j}
!!              \f}
!!              with \f$d_j\f$ being the filter coeffficient
!!          + wavelet filter
!!          + bogey shock filter
subroutine filter_block(filter,time, u,Bs, g, x0, dx, work_array)

    implicit none
     !> params structure of navier stokes
    type(type_params_filter),intent(in) :: filter
    !> time loop parameters
    real(kind=rk), intent(in)           :: time
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: u(:,:,:,:)
    !> heavy work data array - block data
    real(kind=rk), intent(inout)        :: work_array(:, :, :, :)
    !> grid parameter
    integer(kind=ik),intent(in)         :: Bs, g
    ! spacing and origin of a block
    real(kind=rk)                       :: x0(1:3), dx(1:3)


    ! loop variables
    integer(kind=ik)                    :: i, j, l, dF, N_dF, stencil_size


    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0, sub_t1, time_sum


    ! filtered values and array for old block data
    real(kind=rk)                       :: phi_tilde(3)
    real(kind=rk), allocatable          :: block_old(:, :, :)

    stencil_size            =filter%stencil_size
    N_dF                    =filter%number_data_fields

    if (size(u,3)==1) then
        !2D
        call convert_statevector2D(u(:,:,1,:),'conservative')
    else
        !3D
        call abort(5326,"ERROR [filter_block]: 3D case not implemented yet")
    endif
    work_array(:,:,:,1:N_dF)= u

    ! loop over all datafields
    do dF = 1, N_dF
        ! switch case
        ! explicit filter -> stencil_size /= 0
        ! wavelet filter -> stencil_size == 0
        if (stencil_size /= 0) then
            ! explicit filter
            ! save old block data
            !block_old = hvy_block(:, :, :, dF, hvy_active(k) )
            ! 3D or 2D case
            if ( size(u,3) > 1 ) then
                ! 3D
                ! loop over block data
                do i = g+1, Bs+g
                    do j = g+1, Bs+g
                        do l = g+1, Bs+g
                            ! x direction
                            call filter_1D( work_array(i-( (stencil_size+1)/2-1):i+( (stencil_size+1)/2-1), j, l,dF ), phi_tilde(1), filter%stencil(1:stencil_size) )
                            ! y direction
                            call filter_1D( work_array(i, j-( (stencil_size+1)/2-1):j+( (stencil_size+1)/2-1), l,dF ), phi_tilde(2), filter%stencil(1:stencil_size) )
                            ! z direction
                            call filter_1D( work_array(i, j, l-( (stencil_size+1)/2-1):l+( (stencil_size+1)/2-1),dF ), phi_tilde(3), filter%stencil(1:stencil_size) )
                            ! filter
                            u(i, j, l, dF ) = u(i, j, l, dF ) + phi_tilde(1) + phi_tilde(2) + phi_tilde(3)
                        end do
                    end do
                end do
            else
                ! 2D
                ! loop over block data
                do i = g+1, Bs+g
                    do j = g+1, Bs+g
                        ! x direction
                        call filter_1D( work_array(i-( (stencil_size+1)/2-1):i+( (stencil_size+1)/2-1), j, 1, dF ), phi_tilde(1), filter%stencil(1:stencil_size) )
                        ! y direction
                        call filter_1D( work_array(i, j-( (stencil_size+1)/2-1):j+( (stencil_size+1)/2-1), 1, dF ), phi_tilde(2), filter%stencil(1:stencil_size) )
                        ! filter
                        u(i, j, 1, dF ) = u(i, j, 1, dF) + phi_tilde(1) + phi_tilde(2)
                    end do
                end do
            end if
        elseif (stencil_size == 0) then
            select case(filter%name)
                case('wavelet')
                    ! wavelet filter
                    call wavelet_filter(filter,Bs,g, u(:, :, :, dF))
                case('bogey_shock')
                    ! shock filter
                    if ( dF == 1 ) then
                        call bogey_filter(filter, Bs, g, N_dF ,u,x0,dx,work_array)
                    end if
            end select
        end if
    end do

    if (size(u,3)==1) then
        call pack_statevector2D(u(:,:,1,:),'conservative')
    else
        call abort(9820,"3D not implemented yet")
    endif

end subroutine filter_block














!=====================================================================
!  1D FILTER
!=====================================================================


!> \details
!> \name filter_1D.f90
!> \version 0.5
!> \author msr
!
!> \brief 1D Filter subroutine
!
!>
!! input:    - filter stencil, data array, position of value to filter \n
!! output:   - filtered data \n
!!
!!
!! = log ======================================================================================
!! \n
!! 27/03/17 - create
!! 02/05/17 - return filtered value separatly
!
! ********************************************************************************************

subroutine filter_1D(phi, phi_tilde, a)

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> datafield
    real(kind=rk), intent(in)           :: phi(:)
    !> filtered value
    real(kind=rk), intent(out)          :: phi_tilde
    !> filter coefficients
    real(kind=rk), intent(in)           :: a(:)

    ! loop variable
    integer(kind=ik)                    :: k

    ! old values
    real(kind=rk)                       :: phi_old(size(phi,1))

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    phi_old   = phi
    phi_tilde = 0.0_rk

!---------------------------------------------------------------------------------------------
! main body

    ! check filter stencil
    if ( size(phi) /= size(a) ) then
        write(*,'(80("_"))')
        print*, phi
        print*, a
        write(*,*) "ERROR: filter stencil has wrong size"
        stop
    end if

    ! filter data
    do k = 1, size(a)
        phi_tilde = phi_tilde + a(k)*phi_old(k)
    end do

end subroutine filter_1D



















!=====================================================================
!  WAVELET FILTER
!=====================================================================

!> \file
!> \callgraph
!> \name wavelet_filter.f90
!> \version 0.5
!> \author msr
!
!> \brief wavelet filter subroutine
!! \date 24/07/17 - create
!! \date 29/04/18 - update for new physics branch (pKrah)
subroutine wavelet_filter(filter,Bs,g, block_data)
    use module_interpolation, only :    restriction_3D,restriction_2D,&
                                        prediction_2D,prediction_3D

    implicit none
    !> params structure of navier stokes
    type(type_params_filter),intent(in) :: filter
    !> mesh params
    integer(kind=ik),   intent(in)      :: Bs
    integer(kind=ik),   intent(in)      :: g
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: block_data(:, :, :)

    ! loop parameter
    integer(kind=ik)                    :: i, j, l
    ! detail
    real(kind=rk)                       :: detail
    ! grid parameter
    ! interpolation fields
    real(kind=rk)                        ::  u1(Bs+2*g,Bs+2*g,Bs+2*g), &
                                            u2(Bs+2*g,Bs+2*g,Bs+2*g), &
                                            u3((Bs+1)/2+g,(Bs+1)/2+g,(Bs+1)/2+g)


    ! reset detail
    detail = 0.0_rk

    if ( size(block_data,3)>1 ) then
        ! ********** 3D **********
        ! copy block data to array u1
        u1(:,:,:) = block_data( :, :, : )
        ! now, coarsen array u1 (restriction)
        call restriction_3D( u1, u3 )  ! fine, coarse
        ! then, re-interpolate to the initial level (prediciton)
        call prediction_3D ( u3, u2, filter%order_predictor )  ! coarse, fine

        ! Calculate detail by comparing u1 (original data) and u2 (result of predict(restrict(u1)))
        ! NOTE: the error (or detail) is evaluated on the entire block, INCLUDING the ghost nodes layer
        do i = 1, Bs+2*g
            do j = 1, Bs+2*g
                do l = 1, Bs+2*g
                    detail = max( detail, sqrt( (u1(i,j,l)-u2(i,j,l)) * ( u1(i,j,l)-u2(i,j,l)) ) )
                end do
            end do
        end do

        ! evaluate criterion: if this blocks detail is smaller than the prescribed precision,
        ! the block should be filtered, overwrite block data with predicted data
        if (detail < filter%eps) then
            ! wavelet filtering
            !block_data(:,:,:) = u2(:,:,:)
            ! note: do not filter redundant nodes, to avoid instabilities
            block_data(g+2:Bs+g-1,g+2:Bs+g-1,g+2:Bs+g-1) = u2(g+2:Bs+g-1,g+2:Bs+g-1,g+2:Bs+g-1)
        end if

    else
        ! ********** 2D **********
        ! copy block data to array u1
         u1(:,:,1) = block_data( :, :, 1 )
        ! now, coarsen array u1 (restriction)
        call restriction_2D( u1(:,:,1), u3(:,:,1) )  ! fine, coarse
        ! then, re-interpolate to the initial level (prediciton)
        call prediction_2D ( u3(:,:,1), u2(:,:,1), filter%order_predictor )  ! coarse, fine

        ! Calculate detail by comparing u1 (original data) and u2 (result of predict(restrict(u1)))
        ! NOTE: the error (or detail) is evaluated on the entire block, INCLUDING the ghost nodes layer
        ! do i = 1, Bs+2*g
        !     do j = 1, Bs+2*g
        !         detail = max( detail, sqrt( (u1(i,j,1)-u2(i,j,1)) * ( u1(i,j,1)-u2(i,j,1)) ) )
        !     end do
        ! end do

        ! evaluate criterion: if this blocks detail is smaller than the prescribed precision,
        ! the block should be filtered, overwrite block data with predicted data
        ! if (detail < params%eps) then
            ! wavelet filtering
            !block_data(:,:,1) = u2(:,:,1)
            ! note: do not filter redundant nodes, to avoid instabilities
            block_data(g+2:Bs+g-1,g+2:Bs+g-1,1) = u2(g+2:Bs+g-1,g+2:Bs+g-1,1)
        ! end if

    end if

end subroutine wavelet_filter


























!=====================================================================
!  BOGEY FILTER
!=====================================================================
!> \details
!> \callgraph
!> \name bogey_filter.f90
!> \version 0.5
!> \author msr
!
!> \brief bogey shock filter subroutine

!! \date 21/09/17 - create
!! \date 29/04/18 - update for new physics branch (pKrah)

subroutine bogey_filter(filter, Bs, g, N_dF, block_data, xx0, ddx, hvy_work)

!---------------------------------------------------------------------------------------------
! variables

    implicit none
    !> params structure of navier stokes
    type(type_params_filter),intent(in) :: filter
    !> grid parameter
    integer(kind=ik), intent(in)        :: g, Bs, N_dF
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: block_data(:, :, :, :)
    !> heavy work
    real(kind=rk), intent(inout)        :: hvy_work(:, :, :, :)
    ! spacing and origin of a block
    real(kind=rk), intent(in)           :: xx0(1:3), ddx(1:3)


    ! loop parameter
    integer(kind=ik)                    :: i, j, l, dF

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
    r_th = filter%r_th!/( real(level,kind=rk) )**2.0_rk

    detector_method = filter%detector_method ! "p", "divU"
    sigma_method    = filter%sigma_switch ! "abs", "tanh"

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
    if (.true.) then
        !! CAUTION !!
        ! conservative variables are asumed here (rho,rho u, rho v, e)
        call convert2format(block_data(:,:,1,:)    ,'conservative',&
                            hvy_WORK(:,:,1,1:N_dF) ,'pure_variables')

        rho = hvy_WORK(:, :, 1, 1)
        u   = hvy_work(:, :, 1, 2)
        v   = hvy_work(:, :, 1, 3)
        p   = hvy_work(:, :, 1, 4)
        gamma_=params_ns%gamma_

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

    if (filter%save_filter_strength) then
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
