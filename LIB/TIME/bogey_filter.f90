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

subroutine bogey_filter( params, block_data)

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params

    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: block_data(:, :, :, :)

    ! loop parameter
    integer(kind=ik)                    :: i, j, l, dF
    ! grid parameter
    integer(kind=ik)                    :: Bs, g, N_dF

    ! filtered values and array for old block data
    real(kind=rk)                       :: phi_tilde(3), r_xyz(3), r_th
    real(kind=rk), allocatable          :: block_old(:, :, :, :), r(:,:,:), sigma(:,:,:)

    ! stencil array, note: size is fixed
    real(kind=rk)                       :: stencil(3)

    ! filter position (array postion of value to filter)
    integer(kind=ik)                    :: stencil_size

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization
!
    ! grid parameter
    Bs = params%number_block_nodes
    g  = params%number_ghost_nodes

    N_dF = params%number_data_fields

    ! allocate old data array
    allocate( block_old(Bs+2*g, Bs+2*g, Bs+2*g, N_dF), r(Bs+2*g, Bs+2*g, Bs+2*g), sigma(Bs+2*g, Bs+2*g, Bs+2*g) )

    stencil_size = 3
    stencil(1:stencil_size) = (/  1.0_rk/  4.0_rk, &
                                 -1.0_rk/  2.0_rk, &
                                  1.0_rk/  4.0_rk/)

    ! threshold value
    r_th = params%r_th

    r = 1.0e-16_rk

!---------------------------------------------------------------------------------------------
! main body

    ! first: filter values, note: also filter first ghost node value!
    ! -----------------------------------------
    ! save old block data
    block_old = block_data

    if ( params%threeD_case ) then

        do dF = 1, N_dF
            do i = g-1, Bs+g+2
                do j = g-1, Bs+g+2
                    do l = g-1, Bs+g+2

                        ! x direction
                        call filter_1D( block_old(i-( (stencil_size+1)/2-1):i+( (stencil_size+1)/2-1), j, l, dF ), phi_tilde(1), stencil(1:stencil_size) )
                        ! y direction
                        call filter_1D( block_old(i, j-( (stencil_size+1)/2-1):j+( (stencil_size+1)/2-1), l, dF ), phi_tilde(2), stencil(1:stencil_size) )
                        ! z direction
                        call filter_1D( block_old(i, j, l-( (stencil_size+1)/2-1):l+( (stencil_size+1)/2-1), dF ), phi_tilde(3), stencil(1:stencil_size) )

                        ! filter
                        block_data(i, j, l, dF ) = block_data(i, j, l, dF ) + phi_tilde(1) + phi_tilde(2) + phi_tilde(3)

                    end do
                end do
            end do
        end do

    else

        do dF = 1, N_dF
            do i = g-1, Bs+g+2
                do j = g-1, Bs+g+2

                    ! x direction
                    call filter_1D( block_old(i-( (stencil_size+1)/2-1):i+( (stencil_size+1)/2-1), j, 1, dF ), phi_tilde(1), stencil(1:stencil_size) )
                    ! y direction
                    call filter_1D( block_old(i, j-( (stencil_size+1)/2-1):j+( (stencil_size+1)/2-1), 1, dF ), phi_tilde(2), stencil(1:stencil_size) )

                    ! filter
                    block_data(i, j, 1, dF) = block_data(i, j, 1, dF ) + phi_tilde(1) + phi_tilde(2)

                end do
            end do
        end do

    end if

    ! second: shock detector
    ! -----------------------------------------
    if ( params%threeD_case ) then

        !do dF = 1, N_dF
        dF = 4
            do i = g, Bs+g+1
                do j = g, Bs+g+1
                    do l = g, Bs+g+1

                        ! x direction
                        call shock_detector_1D( block_data(i-1:i+1, j, l, dF), block_old(i, j, l, dF), r_xyz(1) )
                        ! y direction
                        call shock_detector_1D( block_data(i, j-1:j+1, l, dF), block_old(i, j, l, dF), r_xyz(2) )
                        ! z direction
                        call shock_detector_1D( block_data(i, j, l-1:l+1, dF), block_old(i, j, l, dF), r_xyz(3) )
                        ! max r
                        r(i, j, l) = max( r(i, j, l), r_xyz(1), r_xyz(2), r_xyz(3) )
                        ! sigma
                        sigma(i, j, l) = 0.5_rk * ( 1.0_rk - r_th/r(i, j, l) + abs( 1.0_rk - r_th/r(i, j, l)) )

                    end do
                end do
            end do
        !end do

    else

        !do dF = 1, N_dF
        dF = 4
            do i = g, Bs+g+1
                do j = g, Bs+g+1

                        ! x direction
                        call shock_detector_1D( block_data(i-1:i+1, j, 1, dF), block_old(i, j, 1, dF), r_xyz(1) )
                        ! y direction
                        call shock_detector_1D( block_data(i, j-1:j+1, 1, dF), block_old(i, j, 1, dF), r_xyz(2))
                        ! max r
                        r(i, j, 1) = max( r(i, j, 1), r_xyz(1), r_xyz(2) )
                        ! sigma
                        sigma(i, j, 1) = 0.5_rk * ( 1.0_rk - r_th/r(i, j, 1) + abs( 1.0_rk - r_th/r(i, j, 1)) )

                end do
            end do
        !end do

    end if

    ! third: filter values
    ! -----------------------------------------
    if ( params%threeD_case ) then

        do dF = 1, N_dF
            do i = g+1, Bs+g
                do j = g+1, Bs+g
                    do l = g+1, Bs+g

                        ! x direction
                        phi_tilde(1) = 0.5_rk*( sigma(i+1, j, l) + sigma(i, j, l) ) &
                                       * ( block_old(i-1, j, l, dF ) * stencil(1) + block_old(i, j, l, dF ) * stencil(2)/2.0_rk ) &
                                       + 0.5_rk*( sigma(i-1, j, l) + sigma(i, j, l) ) &
                                       * ( block_old(i, j, l, dF ) * stencil(2)/2.0_rk + block_old(i+1, j, l, dF ) * stencil(3))
                        ! y direction
                        phi_tilde(2) = 0.5_rk*( sigma(i, j+1, l) + sigma(i, j, l) ) &
                                       * ( block_old(i, j-1, l, dF ) * stencil(1) + block_old(i, j, l, dF ) * stencil(2)/2.0_rk ) &
                                       + 0.5_rk*( sigma(i, j-1, l) + sigma(i, j, l) ) &
                                       * ( block_old(i, j, l, dF ) * stencil(2)/2.0_rk + block_old(i, j+1, l, dF ) * stencil(3))

                        ! z direction
                        phi_tilde(3) = 0.5_rk*( sigma(i, j, l+1) + sigma(i, j, l) ) &
                                       * ( block_old(i, j, l, dF -1) * stencil(1) + block_old(i, j, l, dF ) * stencil(2)/2.0_rk ) &
                                       + 0.5_rk*( sigma(i, j, l-1) + sigma(i, j, l) ) &
                                       * ( block_old(i, j, l, dF ) * stencil(2)/2.0_rk + block_old(i, j, l+1, dF ) * stencil(3))


                        ! filter
                        block_data(i, j, l, dF ) = block_old(i, j, l, dF ) + phi_tilde(1) + phi_tilde(2) + phi_tilde(3)

                    end do
                end do
            end do
        end do

    else

        do dF = 1, N_dF
            do i = g+1, Bs+g
                do j = g+1, Bs+g

                    ! x direction
                        phi_tilde(1) = 0.5_rk*( sigma(i+1, j, 1) + sigma(i, j, 1) ) &
                                       * ( block_old(i-1, j, 1, dF ) * stencil(1) + block_old(i, j, 1, dF ) * stencil(2)/2.0_rk ) &
                                       + 0.5_rk*( sigma(i-1, j, 1) + sigma(i, j, 1) ) &
                                       * ( block_old(i, j, 1, dF ) * stencil(2)/2.0_rk + block_old(i+1, j, 1, dF ) * stencil(3))
                        ! y direction
                        phi_tilde(2) = 0.5_rk*( sigma(i, j+1, 1) + sigma(i, j, 1) ) &
                                       * ( block_old(i, j-1, 1, dF ) * stencil(1) + block_old(i, j, 1, dF ) * stencil(2)/2.0_rk ) &
                                       + 0.5_rk*( sigma(i, j-1, 1) + sigma(i, j, 1) ) &
                                       * ( block_old(i, j, 1, dF ) * stencil(2)/2.0_rk + block_old(i, j+1, 1, dF ) * stencil(3))

                        ! filter
                        block_data(i, j, 1, dF ) = block_old(i, j, 1, dF ) + phi_tilde(1) + phi_tilde(2)

                end do
            end do
        end do

    end if

    ! clean up
    deallocate(block_old, r, sigma)

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
