!---------------------------------------------------!!!!
!> \file This file includes all filter routines
!---------------------------------------------------!!!!


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

            case('bogey_shock')
                ! bogey shock detector threshold
                call read_param_mpi(FILE, 'Discretization', 'r_th', filter%r_th, 1e-3_rk )
                ! bogey shock switch tanh
                call read_param_mpi(FILE, 'Discretization', 'switch', filter%sigma_switch, 'tanh' )
                ! bogey shock detection method (p,divU)
                call read_param_mpi(FILE, 'Discretization', 'detector_method', filter%detector_method, 'divU' )

            case('wavelet')
              call read_param_mpi(FILE, 'Discretization', 'order_predictor', filter%order_predictor, "multiresolution_4th")

            case('no_filter')
                ! do nothing..
                return

            case default
                call abort(4564,"ERROR [filter_block.f90]: filter type is not known!")
        end select

        call read_param_mpi(FILE, 'Blocks', 'number_equations', filter%n_eqn, 1 )
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
subroutine filter_block(filter, time, u, g, Bs, x0, dx, work_array)

    implicit none
     !> params structure of navier stokes
    type(type_params_filter),intent(in) :: filter
    !> time loop parameters
    real(kind=rk), intent(in)           :: time
    !> heavy work data array - block data
    real(kind=rk), intent(inout)        :: work_array(:, :, :, :)
    real(kind=rk), intent(inout)           :: u(:, :, :, :)
    !> grid parameter
    integer(kind=ik),intent(in)         :: Bs, g
    ! spacing and origin of a block
    real(kind=rk)                       :: x0(1:3), dx(1:3)


    ! loop variables
    integer(kind=ik)                    :: i, j, l, dF, N_dF,dF_old, stencil_size
    ! cpu time variables for running time calculation
    real(kind=rk)                       :: sub_t0, sub_t1, time_sum

    ! filtered values and array for old block data
    real(kind=rk)                       :: phi_tilde(3)

    stencil_size            =filter%stencil_size
    N_dF                    =filter%n_eqn

    work_array(:,:,:,1:N_dF)=u
    ! use conservative form of statevector for filtering only!!!!
    call convert_statevector(work_array(:,:,:,1:N_dF),'conservative')

    ! loop over all datafields
    work_array(:,:,:,N_dF+1:2*N_dF)=work_array(:,:,:,1:N_dF)

    select case (filter%name)
    case ('bogey_shock')
      if (params_ns%dim==3) then
        call bogey_filter3D(filter, u, g, Bs, N_dF, x0, dx, work_array)
      else
        call bogey_filter2D_(filter, u, g, Bs, N_dF, x0, dx, work_array) !new bogey_filter
      !  call bogey_filter2D(filter, Bs, g, N_dF ,work_array,x0,dx) !old bogey filter
      endif
    case('explicit_5pt','explicit_7pt','explicit_9pt','explicit_11pt') ! explicit filtering
      do dF = 1, N_dF
          dF_old=dF+N_dF
              !block_old = hvy_block(:, :, :, dF, hvy_active(k) )
              ! 3D or 2D case
              if (params_ns%dim==3 ) then
                  ! 3D
                  ! loop over block data
                  do i = g+1, Bs+g
                      do j = g+1, Bs+g
                          do l = g+1, Bs+g
                              ! x direction
                              call filter_1D( work_array(i-( (stencil_size+1)/2-1):i+( (stencil_size+1)/2-1), j, l,dF_old ), phi_tilde(1), filter%stencil(1:stencil_size) )
                              ! y direction
                              call filter_1D( work_array(i, j-( (stencil_size+1)/2-1):j+( (stencil_size+1)/2-1), l,dF_old ), phi_tilde(2), filter%stencil(1:stencil_size) )
                              ! z direction
                              call filter_1D( work_array(i, j, l-( (stencil_size+1)/2-1):l+( (stencil_size+1)/2-1),dF_old ), phi_tilde(3), filter%stencil(1:stencil_size) )
                              ! filter
                              work_array(i, j, l, dF ) = work_array(i, j, l, dF_old ) + phi_tilde(1) + phi_tilde(2) + phi_tilde(3)
                          end do
                      end do
                  end do
              else
                  ! 2D
                  ! loop over block data
                  do i = g+1, Bs+g
                      do j = g+1, Bs+g
                          ! x direction
                          call filter_1D( work_array(i-( (stencil_size+1)/2-1):i+( (stencil_size+1)/2-1), j, 1, dF_old ), phi_tilde(1), filter%stencil(1:stencil_size) )
                          ! y direction
                          call filter_1D( work_array(i, j-( (stencil_size+1)/2-1):j+( (stencil_size+1)/2-1), 1, dF_old ), phi_tilde(2), filter%stencil(1:stencil_size) )
                          ! filter
                          work_array(i, j, 1, dF ) = work_array(i, j, 1, dF_old) + phi_tilde(1) + phi_tilde(2)
                      end do
                  end do

              endif
      end do

    case('wavelet')
      do dF = 1, N_dF
          call wavelet_filter(filter%order_predictor, Bs, g, work_array(:,:,:,dF))
      enddo
    case default
      call abort(100918,"No filter called: "//filter%name //" is known.")

    end select
    ! pack statevector from conservative form !!!
    call pack_statevector(work_array,'conservative')


end subroutine filter_block








!=====================================================================
!  WAVELET FILTER
!=====================================================================
subroutine wavelet_filter(order_predictor, Bs, g, block_data)
    use module_interpolation, only :  restriction_3D, restriction_2D, prediction_2D, prediction_3D

    implicit none
    !> params structure of navier stokes
    character(len=*), intent(in) :: order_predictor
    !> mesh params
    integer(kind=ik), intent(in) :: Bs
    integer(kind=ik), intent(in) :: g
    !> heavy data array - block data
    real(kind=rk), intent(inout) :: block_data(:, :, :)
    real(kind=rk), allocatable, save :: u3(:,:,:)


    if ( size(block_data,3)>1 ) then
        ! ********** 3D **********
        if (.not.allocated(u3)) allocate(u3((Bs+1)/2+g,(Bs+1)/2+g,(Bs+1)/2+g))
        ! now, coarsen array u1 (restriction)
        call restriction_3D( block_data, u3 )  ! fine, coarse
        ! then, re-interpolate to the initial level (prediciton)
        call prediction_3D ( u3, block_data, order_predictor )  ! coarse, fine

    else
        ! ********** 2D **********
        if (.not.allocated(u3)) allocate(u3((Bs+1)/2+g,(Bs+1)/2+g,1))
        ! now, coarsen array u1 (restriction)
        call restriction_2D( block_data(:,:,1), u3(:,:,1) )  ! fine, coarse
        ! then, re-interpolate to the initial level (prediciton)
        call prediction_2D ( u3(:,:,1), block_data(:,:,1), order_predictor )  ! coarse, fine

    end if

end subroutine wavelet_filter







!=====================================================================
!  1D FILTER
!=====================================================================
!> \brief 1D Filter subroutine
!> \details
!> \version 0.5
!> \author msr
!! \date 27/03/17 - create
!! \date 02/05/17 - return filtered value separatly
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
        call abort(123980,"ERROR: filter stencil has wrong size")
    end if

    ! filter data
    do k = 1, size(a)
        phi_tilde = phi_tilde + a(k)*phi_old(k)
    end do

end subroutine filter_1D







!=========================================================================================
!> 3D implementation of the bogey-bailey shock capturing filter
!> \details Publication: <a href="https://doi.org/10.1016/j.jcp.2008.10.042"> DOI: 10.1016/j.jcp.2008.10.042</a>
!> This implementation only supports the div(u) as thresholding indicator and the filtering
!> strength is computed from sigma=1 - tanh( r_th/r/0.7 ).
!> I have rewritten the filter to make it fast and simple!
!> \author P.Krah
subroutine bogey_filter3D(filter, u, g, Bs, N_dF, xx0, ddx, work_array)

  implicit none
  !-----------------------------------------------------------------------
  type(type_params_filter),intent(in) :: filter      !< filter arguments
  integer(kind=ik), intent(in)        :: g, Bs, N_dF       !< grid parameter
  real(kind=rk), intent(inout)        :: work_array(:, :, :, :) !< block data and work arrays
  real(kind=rk), intent(inout)        :: u(:,:,:,:)
  real(kind=rk), intent(in)           :: xx0(1:3), ddx(1:3) !< spacing and origin of a block
  !-----------------------------------------------------------------------

  integer(kind=ik) :: ix,iy,iz,i,dF,d !loop variables
  real(kind=rk)    :: r_th, r, c_stencil(4)
  real(kind=rk)    :: block_old(size(work_array,1), size(work_array,2), size(work_array,3), N_dF)
  integer(kind=ik),parameter :: stencil_size=3, SHIFT(7)=(/-3, -2, -1, 0, 1, 2, 3/),DIM=3
  real(kind=rk), parameter   ::  eps=1e-16, c1=-0.210383_rk, c2 = 0.030617_rk
  real(kind=rk)               :: stencil(3),gamma_,r_rs(2),Dthetamag

  ! PARAMETERS
  real(kind=rk),save    :: divu(size(SHIFT),3),soundspeed2(size(SHIFT),3) &
  ,Dtheta(size(SHIFT),3),sigma(size(SHIFT),3)
  ! /todo: move to ini file
  character(len=80)                   :: detector_method, sigma_method

  stencil(1:stencil_size) = (/  1.0_rk/  4.0_rk, &
  -1.0_rk/  2.0_rk, &
  1.0_rk/  4.0_rk/)

  ! threshold value
  r_th = filter%r_th!/( real(level,kind=rk) )**2.0_rk

  detector_method = filter%detector_method ! "p", "divU"
  sigma_method    = filter%sigma_switch ! "abs", "tanh"


  c_stencil = (/ -c2, -c1, c1, c2 /);
  block_old = work_array(:,:,:,1:N_dF)
  call convert2format( block_old    ,'conservative',&
  work_array(:,:,:,1:N_dF)  ,'pure_variables')
  call convert_statevector(u,'pure_variables')

  gamma_=params_ns%gamma_

  ! This routine is ugly but fucking fast and it saves memory!!
  ! Because filtering is done locally in every step.
  do iz = g+1, Bs+g
    do iy = g+1, Bs+g
      do ix = g+1, Bs+g
        do i=1,7 !shift the stencil loop
          !This loop shifts the stencil in every direction   direction shift(1)=-3,shift(2)=-2,...shift(6)=3

          ! compute the divergence
          ! divergence shifted in x
          divu(i,1)= ( u(ix+1+SHIFT(i),iy,iz,UxF) - u(ix-1+SHIFT(i),iy,iz,UxF) ) / (2.0_rk*ddx(1)) &
          + ( u(ix+SHIFT(i),iy+1,iz,UyF) - u(ix+SHIFT(i),iy-1,iz,UyF) ) / (2.0_rk*ddx(2)) &
          + ( u(ix+SHIFT(i),iy,iz+1,UzF) - u(ix+SHIFT(i),iy,iz-1,UzF) ) / (2.0_rk*ddx(3))
          ! shifted in y
          divu(i,2)= ( u(ix+1,iy+SHIFT(i),iz,UxF) - u(ix-1,iy+SHIFT(i),iz,UxF) ) / (2.0_rk*ddx(1)) &
          + ( u(ix,iy+SHIFT(i)+1,iz,UyF) - u(ix,iy+SHIFT(i)-1,iz,UyF) ) / (2.0_rk*ddx(2)) &
          + ( u(ix,iy+SHIFT(i),iz+1,UzF) - u(ix,iy+SHIFT(i),iz-1,UzF) ) / (2.0_rk*ddx(3))
          ! shifted in z
          divu(i,3)= ( u(ix+1,iy,iz+SHIFT(i),UxF) - u(ix-1,iy,iz+SHIFT(i),UxF) ) / (2.0_rk*ddx(1)) &
          + ( u(ix,iy+1,iz+SHIFT(i),UyF) - u(ix,iy-1,iz+SHIFT(i),UyF) ) / (2.0_rk*ddx(2)) &
          + ( u(ix,iy,iz+SHIFT(i)+1,UzF) - u(ix,iy,iz+SHIFT(i)-1,UzF) ) / (2.0_rk*ddx(3))

          ! speed of sound is actually not needed for i=1,2 and 6,7
          soundspeed2(i,1)= gamma_ * u(ix+SHIFT(i),iy,iz,pF) / u(ix+SHIFT(i),iy,iz,rhoF) ! c^2 shifted in x
          soundspeed2(i,2)= gamma_ * u(ix,iy+SHIFT(i),iz,pF) / u(ix,iy+SHIFT(i),iz,rhoF) ! c^2 shifted in y
          soundspeed2(i,3)= gamma_ * u(ix,iy,iz+SHIFT(i),pF) / u(ix,iy,iz+SHIFT(i),rhoF) ! c^2 shifted in z

        enddo

        ! from the dilatation we can compute the filtersing strength sigma:
        do d=1,DIM
          ! compute a second order filter as the filtering indicator
          do i=2,6
            call filter_1D( divu(i-1:i+1, d), Dtheta(i,d), stencil(1:stencil_size) )
          enddo
          do i=3,5
            ! compute its magnitude
            Dthetamag = 0.5_rk * ( ( Dtheta(i,d) - Dtheta(i+1,d) )**2.0_rk &
            + ( Dtheta(i,d) - Dtheta(i-1,d) )**2.0_rk )
            ! normalize the magnitude to the local speed of sound
            r = Dthetamag / (soundspeed2(i,d)/ddx(d)**2.0_rk) + eps
            ! compute the filtering strength
            sigma(i,d) = 1.0_rk - dtanh( r_th/r/0.7_rk )
          enddo
          ! save the filtering strength
          i=4 ! SHIFT(4)=0 (no shift)
          work_array(ix,iy,iz,N_dF+d)=sigma(i,d)
        enddo ! loop over dimensions

        do dF= 1, N_dF
          work_array(ix,iy,iz,dF) = block_old(ix,iy,iz,dF)
          work_array(ix,iy,iz,dF) = work_array(ix,iy,iz ,dF) &
          !filtering in x
          - ( 0.5_rk * (sigma(i+1,1) + sigma(i,1)) * sum( c_stencil * block_old(ix-1:ix+2,iy,iz,dF) ) &
          -   0.5_rk * (sigma(i-1,1) + sigma(i,1)) * sum( c_stencil * block_old(ix-2:ix+1,iy,iz,dF) ) ) &
          !filtering in y
          - ( 0.5_rk * (sigma(i+1,2) + sigma(i,2)) * sum( c_stencil * block_old(ix,iy-1:iy+2,iz,dF) ) &
          -   0.5_rk * (sigma(i-1,2) + sigma(i,2)) * sum( c_stencil * block_old(ix,iy-2:iy+1,iz,dF) ) ) &
          !filtering in z
          - ( 0.5_rk * (sigma(i+1,3) + sigma(i,3)) * sum( c_stencil * block_old(ix,iy,iz-1:iz+2,dF) ) &
          -   0.5_rk * (sigma(i-1,3) + sigma(i,3)) * sum( c_stencil * block_old(ix,iy,iz-2:iz+1,dF) ) )
        enddo ! loop over degrees of freedom
      enddo !loop over ix
    enddo ! loop over iy
  enddo ! loop over iz



end subroutine bogey_filter3D
!=========================================================================================





!=========================================================================================
!> 2D implementation of the bogey-bailey shock capturing filter
!> \details Publication: <a href="https://doi.org/10.1016/j.jcp.2008.10.042"> DOI: 10.1016/j.jcp.2008.10.042</a>
!> This implementation only supports the div(u) as thresholding indicator and the filtering
!> strength is computed from sigma=1 - tanh( r_th/r/0.7 ).
!> I have rewritten the filter to make it fast and simple!
!> \author P.Krah
subroutine bogey_filter2D_(filter, u, g, Bs, N_dF, xx0, ddx, work_array)

  implicit none
  !-----------------------------------------------------------------------
  type(type_params_filter),intent(in) :: filter      !< filter arguments
  integer(kind=ik), intent(in)        :: g, Bs, N_dF       !< grid parameter
  real(kind=rk), intent(inout)        :: work_array(:, :, :, :) !< block data and work arrays
  real(kind=rk), intent(inout)        :: u(:,:,:,:)
  real(kind=rk), intent(in)           :: xx0(1:3), ddx(1:3) !< spacing and origin of a block
  !-----------------------------------------------------------------------

  integer(kind=ik) :: ix,iy,iz,dF,i,d !loop variables
  real(kind=rk)    :: r_th, r, c_stencil(4)
  real(kind=rk)    :: block_old(size(work_array,1), size(work_array,2), size(work_array,3), N_dF)
  integer(kind=ik),parameter :: stencil_size=3, SHIFT(7)=(/-3, -2, -1, 0, 1, 2, 3/),DIM=2
  real(kind=rk), parameter   ::  eps=1e-16, c1=-0.210383_rk, c2 = 0.030617_rk
  real(kind=rk)               :: stencil(3),gamma_,r_rs(2),Dthetamag

  ! PARAMETERS
  real(kind=rk),save    :: divu(size(SHIFT),3),soundspeed2(size(SHIFT),3) &
  ,Dtheta(size(SHIFT),3),sigma(size(SHIFT),3)
  ! /todo: move to ini file
  character(len=80)                   :: detector_method, sigma_method

  stencil(1:stencil_size) = (/  1.0_rk/  4.0_rk, &
  -1.0_rk/  2.0_rk, &
  1.0_rk/  4.0_rk/)

  ! threshold value
  r_th = filter%r_th!/( real(level,kind=rk) )**2.0_rk

  detector_method = filter%detector_method ! "p", "divU"
  sigma_method    = filter%sigma_switch ! "abs", "tanh"

  c_stencil = (/ -c2, -c1, c1, c2 /);

  ! copy the non filtert data work_array to block_old
  do dF=1,N_dF
      do iy = 1, Bs+2*g, 1
         do ix = 1, Bs+2*g, 1
              block_old(ix,iy,1,df) = work_array(ix,iy,1,dF)
          end do
      end do
  end do
  call convert_statevector(u,'pure_variables')

  gamma_=params_ns%gamma_

  ! This routine is ugly but fucking fast and it saves memory!!
  ! Because filtering is done locally in every step.
  iz=1
  do iy = g+1, Bs+g
    do ix = g+1, Bs+g
      do i=1,7 !shift the stencil loop
        !This loop shifts the stencil in every direction   direction shift(1)=-1,shift(2)=0,shift(3)=1

        ! compute the divergence
        ! divergence shifted in x
        divu(i,1)= ( u(ix+1+SHIFT(i),iy,iz,UxF) - u(ix-1+SHIFT(i),iy,iz,UxF) ) / (2.0_rk*ddx(1)) &
        + ( u(ix+SHIFT(i),iy+1,iz,UyF) - u(ix+SHIFT(i),iy-1,iz,UyF) ) / (2.0_rk*ddx(2))
        ! shifted in y
        divu(i,2)= ( u(ix+1,iy+SHIFT(i),iz,UxF) - u(ix-1,iy+SHIFT(i),iz,UxF) ) / (2.0_rk*ddx(1)) &
        + ( u(ix,iy+SHIFT(i)+1,iz,UyF) - u(ix,iy+SHIFT(i)-1,iz,UyF) ) / (2.0_rk*ddx(2))

        ! speed of sound is actually not needed for i=1,2 and 6,7
        soundspeed2(i,1)= gamma_ * u(ix+SHIFT(i),iy,iz,pF) / u(ix+SHIFT(i),iy,iz,rhoF) ! c^2 shifted in x
        soundspeed2(i,2)= gamma_ * u(ix,iy+SHIFT(i),iz,pF) / u(ix,iy+SHIFT(i),iz,rhoF) ! c^2 shifted in y

      enddo

      ! from the dilatation we can compute the filtersing strength sigma:
      do d=1,DIM
        ! compute a second order filter as the filtering indicator
        do i=2,6
          call filter_1D( divu(i-1:i+1, d), Dtheta(i,d), stencil(1:stencil_size) )
        enddo
        do i=3,5
          ! compute its magnitude
          Dthetamag = 0.5_rk * ( ( Dtheta(i,d) - Dtheta(i+1,d) )**2.0_rk &
          + ( Dtheta(i,d) - Dtheta(i-1,d) )**2.0_rk )
          ! normalize the magnitude to the local speed of sound
          r = Dthetamag / (soundspeed2(i,d)/ddx(d)**2.0_rk) + eps
          ! compute the filtering strength
          sigma(i,d) = 1.0_rk - dtanh( r_th/r/0.7_rk )
        enddo

        i=4 ! SHIFT(4)=0 (no shift)
        work_array(ix,iy,iz,N_dF+d)=sigma(i,d) ! save the filtering strength to save it to file
      enddo ! loop over dimensions

      do dF= 1, N_dF
        work_array(ix,iy,iz,dF) = block_old(ix,iy,iz,dF) &
        !filtering in x
        - ( 0.5_rk * (sigma(i+1,1) + sigma(i,1)) * sum( c_stencil * block_old(ix-1:ix+2,iy,iz,dF) ) &
        -   0.5_rk * (sigma(i-1,1) + sigma(i,1)) * sum( c_stencil * block_old(ix-2:ix+1,iy,iz,dF) ) ) &
        !filtering in y
        - ( 0.5_rk * (sigma(i+1,2) + sigma(i,2)) * sum( c_stencil * block_old(ix,iy-1:iy+2,iz,dF) ) &
        -   0.5_rk * (sigma(i-1,2) + sigma(i,2)) * sum( c_stencil * block_old(ix,iy-2:iy+1,iz,dF) ) )
      enddo ! loop over degrees of freedom
    enddo !loop over ix
  enddo ! loop over iy

end subroutine bogey_filter2D_
!=========================================================================================




!=====================================================================
!  BOGEY FILTER, !!!!!!!! OLD !!!!!!!!!!
!=====================================================================
!> \details
!> \callgraph
!> \name bogey_filter2D.f90
!> \version 0.5
!> \author msr
!
!> \brief bogey shock filter subroutine

!! \date 21/09/17 - create
!! \date 29/04/18 - update for new physics branch (pKrah)

!
! subroutine bogey_filter2D(filter, Bs, g, N_dF ,hvy_work, xx0, ddx )
!
! !---------------------------------------------------------------------------------------------
! ! variables
!
!     implicit none
!     !> params structure of navier stokes
!     type(type_params_filter),intent(in) :: filter
!     !> grid parameter
!     integer(kind=ik), intent(in)        :: g, Bs, N_dF
!     ! !> heavy work
!     real(kind=rk), intent(inout)        :: hvy_work(:, :, :, :)
!     ! spacing and origin of a block
!     real(kind=rk), intent(in)           :: xx0(1:3), ddx(1:3)
!
!
!     ! !> heavy data array - block data
!     ! real(kind=rk),                      :: block_data(:, :, :, :)
!     ! loop parameter
!     integer(kind=ik)                    :: i, j, l, dF
!
!     ! filtered values and array for old block data
!     real(kind=rk)                       :: phi_tilde(3), r_xyz(3), r_th, eps, c1, c2, c_stencil(4)
!     real(kind=rk)                      :: block_old(size(hvy_work,1), size(hvy_work,2), size(hvy_work,3), N_dF)!, r_x(:,:), r_y(:,:), sigma_x(:,:), sigma_y(:,:), theta(:,:), u_x(:,:), v_y(:,:), Dtheta_x(:,:), &
! !                                           Dtheta_y(:,:), Dthetamag_x(:,:), Dthetamag_y(:,:), c_2(:,:), rho(:,:), work_array(:,:), v(:,:), p(:,:)
!     real(kind=rk)          :: r_x(Bs+2*g, Bs+2*g), r_y(Bs+2*g, Bs+2*g), sigma_x(Bs+2*g, Bs+2*g), sigma_y(Bs+2*g, Bs+2*g), theta(Bs+2*g, Bs+2*g), &
!               u_x(Bs+2*g, Bs+2*g),v_y(Bs+2*g, Bs+2*g), Dtheta_x(Bs+2*g, Bs+2*g), Dtheta_y(Bs+2*g, Bs+2*g), Dthetamag_x(Bs+2*g, Bs+2*g), &
!               Dthetamag_y(Bs+2*g, Bs+2*g), c_2(Bs+2*g, Bs+2*g), rho(Bs+2*g, Bs+2*g), u(Bs+2*g, Bs+2*g), v(Bs+2*g, Bs+2*g), p(Bs+2*g, Bs+2*g)
!
!     ! stencil array, note: size is fixed
!     real(kind=rk)                       :: stencil(3),gamma_,r_rs(2)
!
!     ! filter position (array postion of value to filter)
!     integer(kind=ik)                    :: stencil_size
!
!     ! /todo: move to ini file
!     character(len=80)                   :: detector_method, sigma_method
!
! !---------------------------------------------------------------------------------------------
! ! interfaces
!
! !---------------------------------------------------------------------------------------------
! ! variables initialization
! !
!     ! grid parameter
!     ! Bs  = params%Bs
!     ! g   = params%n_ghosts
!     ! N_dF= params%n_eqn
!
!     ! allocate old data array
! !    allocate( block_old(Bs+2*g, Bs+2*g, Bs+2*g, N_dF), r_x(Bs+2*g, Bs+2*g), r_y(Bs+2*g, Bs+2*g), sigma_x(Bs+2*g, Bs+2*g), sigma_y(Bs+2*g, Bs+2*g), theta(Bs+2*g, Bs+2*g), &
! !              u_x(Bs+2*g, Bs+2*g),v_y(Bs+2*g, Bs+2*g), Dtheta_x(Bs+2*g, Bs+2*g), Dtheta_y(Bs+2*g, Bs+2*g), Dthetamag_x(Bs+2*g, Bs+2*g), &
! !              Dthetamag_y(Bs+2*g, Bs+2*g), c_2(Bs+2*g, Bs+2*g), rho(Bs+2*g, Bs+2*g), u(Bs+2*g, Bs+2*g), v(Bs+2*g, Bs+2*g), p(Bs+2*g, Bs+2*g) )
!     sigma_x=0.0_rk
!     sigma_y=0.0_rk
!
!     stencil_size = 3
!     stencil(1:stencil_size) = (/  1.0_rk/  4.0_rk, &
!                                  -1.0_rk/  2.0_rk, &
!                                   1.0_rk/  4.0_rk/)
!
!     ! threshold value
!     r_th = filter%r_th!/( real(level,kind=rk) )**2.0_rk
!
!     detector_method = filter%detector_method ! "p", "divU"
!     sigma_method    = filter%sigma_switch ! "abs", "tanh"
!
!     eps = 1e-16
!
!     ! second order
!     ! c1 = -1/4;
!     ! c2 = 0;
!     ! fourth order
!     ! c1 = -3/16;
!     ! c2 = 1/16;
!     ! fourth order - opt
!     c1 = -0.210383_rk;
!     c2 = 0.030617_rk;
!
!     c_stencil = (/ -c2, -c1, c1, c2 /);
!     block_old=hvy_work(:,:,:,1:N_dF)
!     call convert2format(block_old    ,'conservative',&
!                         hvy_work(:,:,:,1:N_dF)  ,'pure_variables')
!
!
!     rho = hvy_work(:, :, 1, rhoF)
!     u   = hvy_work(:, :, 1, UxF)
!     v   = hvy_work(:, :, 1, UyF)
!     p   = hvy_work(:, :, 1, pF)
!     gamma_=params_ns%gamma_
!
! !---------------------------------------------------------------------------------------------
! ! main body
!
!     ! first - shock detector
!     ! detector input
!     select case (detector_method)
!         case("p")
!             theta = p
!
!         case("divU")
!             u_x = 0.0_rk
!             v_y = 0.0_rk
!
!             do i = g-2,Bs+g+3
!                 do j = g-2,Bs+g+3
!
!                     !u_x(i,j) = ( u(i-2,j) - 8.0_rk*u(i-1,j) + 8.0_rk*u(i+1,j) - u(i+2,j) ) / (12.0_rk*ddx(1))
!                     !v_y(i,j) = ( v(i,j-2) - 8.0_rk*v(i,j-1) + 8.0_rk*v(i,j+1) - v(i,j+2) ) / (12.0_rk*ddx(2))
!
!                     u_x(i,j) = ( u(i+1,j) - u(i-1,j) ) / (2.0_rk*ddx(1))
!                     v_y(i,j) = ( v(i,j+1) - v(i,j-1) ) / (2.0_rk*ddx(2))
!
!                 end do
!             end do
!
!             theta = u_x + v_y
!
!     end select
!
!
!     ! filter detector values, note: also filter first two ghost nodes
!     do i = g-1,Bs+g+2
!         do j = g-1,Bs+g+2
!             call filter_1D( theta(i-1:i+1, j), Dtheta_x(i,j), stencil(1:stencil_size) )
!             call filter_1D( theta(i, j-1:j+1), Dtheta_y(i,j), stencil(1:stencil_size) )
!         end do
!     end do
!
!     ! magnitude of filtered values
!     do i = g,Bs+g+1
!         do j = g,Bs+g+1
!             Dthetamag_x(i,j) = 0.5_rk * ( ( Dtheta_x(i,j) - Dtheta_x(i+1,j) )**2.0_rk + ( Dtheta_x(i,j) - Dtheta_x(i-1,j) )**2.0_rk )
!             Dthetamag_y(i,j) = 0.5_rk * ( ( Dtheta_y(i,j) - Dtheta_y(i,j+1) )**2.0_rk + ( Dtheta_y(i,j) - Dtheta_y(i,j-1) )**2.0_rk )
!         end do
!     end do
!
!     ! second - shock sensor
!     select case(detector_method)
!         case("p")
!             do i = g,Bs+g+1
!                 do j = g,Bs+g+1
!                     r_x(i,j) = Dthetamag_x(i,j) / p(i,j)**2.0_rk + eps
!                     r_y(i,j) = Dthetamag_y(i,j) / p(i,j)**2.0_rk + eps
!                 end do
!             end do
!
!         case("divU")
!             !> \todo write physics method interface for speed of sound
!             c_2 = gamma_*p/rho
!             do i = g,Bs+g+1
!                 do j = g,Bs+g+1
!                     r_x(i,j) = Dthetamag_x(i,j) / (c_2(i,j)/ddx(1)**2.0_rk) + eps
!                     r_y(i,j) = Dthetamag_y(i,j) / (c_2(i,j)/ddx(2)**2.0_rk) + eps
!                 end do
!             end do
!
!     end select
!
!     ! third - filter value magnitude
!     select case(sigma_method)
!         case("abs")
!             do i = g,Bs+g+1
!                 do j = g,Bs+g+1
!                     sigma_x(i,j) = 0.5_rk * ( 1.0_rk - r_th/r_x(i,j) + abs( 1.0_rk - r_th/r_x(i,j) ) )
!                     sigma_y(i,j) = 0.5_rk * ( 1.0_rk - r_th/r_y(i,j) + abs( 1.0_rk - r_th/r_y(i,j) ) )
!                 end do
!             end do
!
!         case("tanh")
!             do i = g,Bs+g+1
!                 do j = g,Bs+g+1
!                      sigma_x(i,j) = 1.0_rk - dtanh( r_th/r_x(i,j)/0.7_rk )
!                      sigma_y(i,j) = 1.0_rk - dtanh( r_th/r_y(i,j)/0.7_rk )
!                     ! r_rs(1)     = r_th/(r_x(i,j)-0.5_rk*r_th)
!                     ! r_rs(2)     = r_th/(r_y(i,j)-0.5_rk*r_th)
!                     ! if (r_rs(1)<2.0_rk) then
!                     !     sigma_x(i,j) = 1.0_rk - dtanh( r_rs(1) )
!                     ! else
!                     !     sigma_x(i,j) = 0.0_rk
!                     ! endif
!
!                     ! if (r_rs(1)<2.0_rk) then
!                     !     sigma_y(i,j) = 1.0_rk - dtanh( r_rs(2) )
!                     ! else
!                     !     sigma_y(i,j)=0.0_rk
!                     ! endif
!
!                 end do
!             end do
!
!     end select
!
!     do dF = 1, N_dF
!         do i = g+1,Bs+g
!             do j = g+1,Bs+g
!                 hvy_work(i,j,1,dF) = block_old(i,j,1,dF)
!
!                 hvy_work(i,j,1,dF) = hvy_work(i,j,1,dF) - ( 0.5_rk * (sigma_x(i+1,j) + sigma_x(i,j)) &
!                                                                        * sum( c_stencil * block_old(i-1:i+2,j,1,dF) ) &
!                                                             -   0.5_rk * (sigma_x(i-1,j) + sigma_x(i,j)) &
!                                                                        * sum( c_stencil * block_old(i-2:i+1,j,1,dF) ) )
!
!                 hvy_work(i,j,1,dF) = hvy_work(i,j,1,dF) - ( 0.5_rk * (sigma_y(i,j+1) + sigma_y(i,j)) &
!                                                                        * sum( c_stencil * block_old(i,j-1:j+2,1,dF) ) &
!                                                             -   0.5_rk * (sigma_y(i,j-1) + sigma_y(i,j)) &
!                                                                        * sum( c_stencil * block_old(i,j-2:j+1,1,dF) ) )
!
!             end do
!         end do
!     end do
!
!     if (filter%save_filter_strength) then
!         ! save filter strength in x direction
!         hvy_WORK(:,:,1,N_dF+1)=sigma_x(:,:)
!         ! save filter strength in y direction
!         hvy_WORK(:,:,1,N_dF+2)=sigma_y(:,:)
!     endif
!
! end subroutine bogey_filter2D
