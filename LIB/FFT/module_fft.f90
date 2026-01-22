!-------------------------------------------------------------------------------
!> \brief FORTRAN FFT module for WABBIT
!-------------------------------------------------------------------------------
!> \details
!! This module implements some FFT functions, which can be used on a single block only (!!)\n
!! For the case Jmin=0, this block is usually periodic. \n
!! You are correct to thu_hat, that all this is of no great use for WABBIT, but we can utilize it for the lowest level solution in the Poisson solver. \n
!! This module needs FFTW to be installed/loaded. It still compiles without, but FFT functions cannot be called then \n
!! This module is designed as blackbox, users should only use it's function with input and output real data, \n
!! and not need to touch the complex-valued Fourier-coefficients.
!-------------------------------------------------------------------------------
!  Thanks to the FLUSI coders for the inspiration for this code
!-------------------------------------------------------------------------------
module module_FFT

    use mpi
    use module_params
    use module_helpers
    use module_globals
    use, intrinsic :: iso_c_binding

    implicit none

    ! I usually find it helpful to use the private keyword by itself initially, which specifies
    ! that everything within the module is private unless explicitly marked public.
    PRIVATE

    ! everything is save by default
    SAVE

    ! buffer for forward FFT and backward FFT
    type(C_PTR) :: plan_fft_block, plan_ift_block
    complex(kind=rk), allocatable :: var_hat(:,:,:)

    ! scaling for wavevectors, which are often used
    real(kind=rk) :: scalex, scaley, scalez

!---------------------------------------------------------------------------------------------
! public parts of this module

#ifdef FFT_ROOT
    include 'fftw3.f03'
#endif
    PUBLIC :: fft_initialize, fft_destroy, fft_is_installed
    PUBLIC :: fft_solve_poisson

contains


!----------------------------------------------------------------
! fft interface functions. These are the ones, that are available
! from the outside and should contain your fft routines
!----------------------------------------------------------------
subroutine fft_solve_poisson(params, u, f, dx)
    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    real(kind=rk), intent(inout) :: u(:, :, :, :)
    real(kind=rk), intent(in) :: f(:, :, :, :)
    real(kind=rk), intent(in) :: dx(1:3)

    ! Local variables
    integer(kind=ik) :: g(1:3), bs(1:3), ic, nc

    bs = params%bs
    g = 0
    g(1:params%dim) = params%g
    nc = size(u,4)

    ! repeat for every variable
    do ic = 1,nc
        u(g(1)+1:g(1)+bs(1),g(2)+1:g(2)+bs(2),g(3)+1:g(3)+bs(3),ic) = f(g(1)+1:g(1)+bs(1),g(2)+1:g(2)+bs(2),g(3)+1:g(3)+bs(3),ic)  ! Copy f to u as input buffer that could be overwritten

        ! compute forward fft
        call fft_forward(params, u(:,:,:,ic), var_hat)
                
        ! compute inverse laplacian u_hat = f_hat/(-k^2), either spectral or FD accuracy
        if (params%FFT_accuracy == "spectral") then
            call invlaplacian_inplace(params, var_hat)
        elseif (params%FFT_accuracy == "FD") then
            call invlaplacian_inplace_filtered_FD(params, var_hat, dx(1:3))
        else
            call abort(250617, "FFT_accuracy unknown "//params%FFT_accuracy//", please use 'spectral' or 'FD'")
        endif

        ! compute backward fft
        call fft_backward(params, var_hat, u(:,:,:,ic))

    enddo

end subroutine

!----------------------------------------------------------------
! fft transformations, assume that block still contains ghost points
! Attention! These could overwrite the input!
!----------------------------------------------------------------
subroutine fft_forward(params, u, u_hat)
    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    real(kind=rk), intent(inout) :: u(:, :, :)
    complex(kind=rk), intent(inout)    :: u_hat(:,:,:)
    integer(kind=ik)  :: bs(1:3), g(1:3)

    bs = params%bs
    g = 0
    g(1:params%dim) = params%g

#ifdef FFT_ROOT
    if (params%dim == 2) then
        call fftw_execute_dft_r2c(plan_fft_block, u(g(1)+1:g(1)+bs(1),g(2)+1:g(2)+bs(2),1), u_hat(:,:,1))
    else
        call fftw_execute_dft_r2c(plan_fft_block, u(g(1)+1:g(1)+bs(1),g(2)+1:g(2)+bs(2),g(3)+1:g(3)+bs(3)), u_hat)
    endif
#else
    call abort(250418, 'Trying to do forward FFT, but code was compiled without fftw3. please load/install fftw3 and recompile')
#endif
end subroutine
subroutine fft_backward(params, u_hat, u)
    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    complex(kind=rk), intent(inout)    :: u_hat(:,:,:)
    real(kind=rk), intent(inout) :: u(:, :, :)
    integer(kind=ik)  :: bs(1:3), g(1:3)

    bs = params%bs
    g = 0
    g(1:params%dim) = params%g

#ifdef FFT_ROOT
    if (params%dim == 2) then
        call fftw_execute_dft_c2r(plan_ift_block, u_hat(:,:,1), u(g(1)+1:g(1)+bs(1),g(2)+1:g(2)+bs(2),1))
        u(g(1)+1:g(1)+bs(1),g(2)+1:g(2)+bs(2),1) = u(g(1)+1:g(1)+bs(1),g(2)+1:g(2)+bs(2),1) / real(product(bs(1:params%dim)), rk)
    else
        call fftw_execute_dft_c2r(plan_ift_block, u_hat, u(g(1)+1:g(1)+bs(1),g(2)+1:g(2)+bs(2),g(3)+1:g(3)+bs(3)))
        u(g(1)+1:g(1)+bs(1),g(2)+1:g(2)+bs(2),g(3)+1:g(3)+bs(3)) = u(g(1)+1:g(1)+bs(1),g(2)+1:g(2)+bs(2),g(3)+1:g(3)+bs(3)) / real(product(bs(1:params%dim)), rk)
    endif
#else
    call abort(250418, 'Trying to do backward FFT, but code was compiled without fftw3. please load/install fftw3 and recompile')
#endif
end subroutine

!----------------------------------------------------------------
! fft derivation, these include:
!    - derivative
!    - laplacian (spectral and FD2 accuracy)
!    - inverse laplacian (spectral and FD2 accuracy)
! What could still be implemented if needed:
!    - divergence
!    - curl
!    - higher order FD accuracy approximation
!----------------------------------------------------------------

! fft derivative, done in one dimension and with varying order (can be negative for inverse!)
subroutine fft_derivative_inplace(params, u_hat, dim, order)
    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    complex(kind=rk), intent(inout)    :: u_hat(:,:,:)
    integer(kind=ik), intent(in)       :: dim     !< dimension for derivative - 1:x, 2:y, 3:z
    integer(kind=ik), intent(in)       :: order   !< order derivative - can be negative
    integer(kind=ik)                   :: i
    real(kind=rk)                      :: k
    complex(kind=rk) :: imag   ! imaginary unit

    imag = dcmplx(0.d0,1.d0)

    ! x-derivative
    if (dim == 1) then
        do i = 1, params%bs(1)/2+1
            k = wave_x(i, params%bs(1))
            if (order < 0 .and. k == 0) then
                u_hat(i,:,:) = 0
            else
                u_hat(i,:,:) = u_hat(i,:,:) * imag*k**order
            endif
        enddo
    ! y-derivative
    elseif (dim == 2) then
        do i = 1, params%bs(2)
            k = wave_y(i, params%bs(2))
            if (order < 0 .and. k == 0) then
                u_hat(:,i,:) = 0
            else
                u_hat(:,i,:) = u_hat(:,i,:) * imag*k**order
            endif
        enddo
    ! z-derivative
    else
        do i = 1, params%bs(3)
            k = wave_z(i, params%bs(3))
            if (order < 0 .and. k == 0) then
                u_hat(:,:,i) = 0
            else
                u_hat(:,:,i) = u_hat(:,:,i) * imag*k**order
            endif
        enddo
    endif
end subroutine

! computes laplace(u_hat) for a scalar valued field and returns it in the same array
subroutine laplacian_inplace( params, u_hat )
    implicit none
    !> parameter struct
    type (type_params), intent(inout)  :: params
    complex(kind=rk),intent(inout)::u_hat(:,:,:)
  
    integer :: ix,iy,iz
    real(kind=rk) :: kx,ky,kz,k2

    do iz = 1, params%bs(3)
        !-- wavenumber in z-direction
        kz = wave_z(iz, params%bs(3))
        do iy = 1, params%bs(2)
            !-- wavenumber in y-direction
            ky = wave_y(iy, params%bs(2))
            do ix = 1, params%bs(1)/2+1
                !-- wavenumber in x-direction
                kx = wave_x(ix, params%bs(1))
                k2 = kx*kx + ky*ky + kz*kz
                u_hat(ix,iy,iz) = -k2*u_hat(ix,iy,iz)
            enddo
        enddo
    enddo
end subroutine laplacian_inplace
  
  
! computes laplace(u_hat) for a scalar valued field and returns it in the same array
! note wavenumbers are reduced to FD accuracy, taken from canuto with general formula:
! 1/dx^2 * (i_0 + 2*sum_{j=1}^{a} (i_j * cos(j*k*dx)))
subroutine laplacian_inplace_filtered_FD( params, u_hat, dx )
    implicit none
    !> parameter struct
    type (type_params), intent(inout)  :: params
    complex(kind=rk),intent(inout)::u_hat(:,:,:)
    real(kind=rk), intent(in) :: dx(1:3)

    integer :: ix,iy,iz
    real(kind=rk) :: kx, ky, kz, k2

    do iz = 1, params%bs(3)
        !-- wavenumber in z-direction
        if (params%dim < 3) then
            kz = 0.0_rk
        else
            kz = wave_z(iz, params%bs(3))
        endif
        do iy = 1, params%bs(2)
            !-- wavenumber in y-direction
            ky = wave_y(iy, params%bs(2))
            do ix = 1, params%bs(1)/2+1
                !-- wavenumber in x-direction
                kx = wave_x(ix, params%bs(1))
                k2 = wave_k2_FD2(kx,ky,kz,dx,params%poisson_order)
                u_hat(ix,iy,iz) = -k2*u_hat(ix,iy,iz)
            enddo
        enddo
    enddo
end subroutine laplacian_inplace_filtered_FD

! computes laplace(u_hat) for a scalar valued field and returns it in the same array
subroutine invlaplacian_inplace( params, u_hat )
    implicit none
    !> parameter struct
    type (type_params), intent(inout)  :: params
    complex(kind=rk),intent(inout)::u_hat(:,:,:)
  
    integer :: ix,iy,iz
    real(kind=rk) :: kx,ky,kz,k2

    do iz = 1, params%bs(3)
        !-- wavenumber in z-direction
        kz = wave_z(iz, params%bs(3))
        do iy = 1, params%bs(2)
            !-- wavenumber in y-direction
            ky = wave_y(iy, params%bs(2))
            do ix = 1, params%bs(1)/2+1
                !-- wavenumber in x-direction
                kx = wave_x(ix, params%bs(1))
                k2 = kx*kx + ky*ky + kz*kz
                if (k2 .eq. 0.0_rk) then
                    u_hat(ix,iy,iz) = 0.0_rk
                else
                    u_hat(ix,iy,iz) = u_hat(ix,iy,iz) / (-k2)
                endif
            enddo
        enddo
    enddo
end subroutine invlaplacian_inplace


! computes laplace(u_hat) for a scalar valued field and returns it in the same array
! note wavenumbers are reduced to FD accuracy, taken from canuto with general formula:
! 1/dx^2 * (i_0 + 2*sum_{j=1}^{a} (i_j * cos(j*k*dx)))
subroutine invlaplacian_inplace_filtered_FD( params, u_hat, dx )
    implicit none
    !> parameter struct
    type (type_params), intent(inout)  :: params
    complex(kind=rk),intent(inout)::u_hat(:,:,:)
    real(kind=rk), intent(in) :: dx(1:3)

    integer :: ix,iy,iz
    real(kind=rk) :: kx,ky,kz,k2
    
    do iz = 1, params%bs(3)
        !-- wavenumber in z-direction
        if (params%dim < 3) then
            kz = 0.0_rk
        else
            kz = wave_z(iz, params%bs(3))
        endif
        do iy = 1, params%bs(2)
            !-- wavenumber in y-direction
            ky = wave_y(iy, params%bs(2))
            do ix = 1, params%bs(1)/2+1
                !-- wavenumber in x-direction
                kx = wave_x(ix, params%bs(1))
                k2 = wave_k2_FD2(kx,ky,kz,dx,params%poisson_order)
                if (k2 .eq. 0.0_rk) then
                    u_hat(ix,iy,iz) = 0.0_rk
                else
                    u_hat(ix,iy,iz) = u_hat(ix,iy,iz) / (-k2)
                endif
            enddo
        enddo
    enddo
end subroutine invlaplacian_inplace_filtered_FD



!----------------------------------------------------------------
! plan initialization and deallocation, needs to be called once before and after any usage
!----------------------------------------------------------------
! initialize fft plan and buffer
subroutine fft_initialize(params)
    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    real(kind=rk), allocatable :: dummy_in(:, :, :)
    integer(kind=ik)  :: bs(1:3), g(1:3)

    bs = params%bs
    g = 0
    g(1:params%dim) = params%g

    ! compute scales
    scalex = 2.0_rk*pi / params%domain_size(1)
    scaley = 2.0_rk*pi / params%domain_size(2)
    if (params%dim == 3) then
        scalez = 2.0_rk*pi / params%domain_size(3)
    else
        scalez = 0.0_rk
    endif

    ! allocate temporary array
    ! allocate dummy_in(params%bs(1), params%bs(2), params%bs(3))
    ! if (.not. allocated(var_hat)) allocate(var_hat(params%bs(1)/2+1, params%bs(2), params%bs(3)))

#ifdef FFT_ROOT
    ! allocate temporary array
    allocate(dummy_in(params%bs(1), params%bs(2), params%bs(3)))
    if (.not. allocated(var_hat)) allocate(var_hat(params%bs(1)/2+1, params%bs(2), params%bs(3)))

    ! allocate plans
    if (params%dim == 2) then
        plan_fft_block = fftw_plan_dft_r2c_2d(bs(1), bs(2), dummy_in(:,:,1), var_hat(:,:,1), FFTW_MEASURE)
        plan_ift_block = fftw_plan_dft_c2r_2d(bs(1), bs(2), var_hat(:,:,1), dummy_in(:,:,1), FFTW_MEASURE)
    else
        plan_fft_block = fftw_plan_dft_r2c_3d(bs(1), bs(2), bs(3), dummy_in, var_hat, FFTW_MEASURE)
        plan_ift_block = fftw_plan_dft_c2r_3d(bs(1), bs(2), bs(3), var_hat, dummy_in, FFTW_MEASURE)
    endif

    ! deallocate dummy array that was used for fft plan creation
    deallocate(dummy_in)
#else
    call abort(250418, 'Trying to initialize FFT plans, but code was compiled without fftw3. please load/install fftw3 and recompile')
#endif
end subroutine


! free fft plan and buffer
subroutine fft_destroy(params)
    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params

#ifdef FFT_ROOT
    call fftw_destroy_plan(plan_fft_block)
    call fftw_destroy_plan(plan_ift_block)
    call fftw_cleanup()
    deallocate(var_hat)
#else
    call abort(250418, 'Trying to free FFT pans, but code was compiled without fftw3. please load/install fftw3 and recompile')
#endif
end subroutine

logical function fft_is_installed()
    implicit none

#ifdef FFT_ROOT
    fft_is_installed = .true.
#else
    fft_is_installed = .false.
#endif
end function

!----------------------------------------------------------------
! wavenumber functions: return the kx,ky,kz wavenumbers
! as a function of the array index
!----------------------------------------------------------------
real(kind=rk) function wave_x( ix, nx )
  implicit none
  integer, intent (in) :: ix, nx
  wave_x = scalex*dble(ix-1)
end function

real(kind=rk) function wave_y( iy, ny )
  implicit none
  integer, intent (in) :: iy, ny
  wave_y = scaley*dble(modulo(iy-1+ny/2,ny)-ny/2)
end function

real(kind=rk) function wave_z( iz, nz )
  implicit none
  integer, intent (in) :: iz, nz
  wave_z = scalez*dble(modulo(iz-1+nz/2,nz)-nz/2)
end function

! reduce wavenumbers to FD accuracy for first order derivatives with stencil coeffs c going from -a, ..., a
! general formula derived from canuto p142:
! 2/dx * (sum_{j=1}^{a} (c_j * sin(j*k*dx)))
real(kind=rk) function wave_k_FD1( k, dx, order_discretization )
    implicit none
    real(kind=rk), intent (in) :: k, dx
    character(len=cshort), intent (in) :: order_discretization
    real(kind=rk) :: kx

    ! if dx is zero then we probably have a 2D case, we just set the return value to 0 for that to proceed
    if (dx == 0.0_rk) then
        wave_k_FD1 = 0.0_rk
        return
    endif
    
    ! ToDo: This is evaluated point-wise, so we should avoid string matching cases
    select case(order_discretization)
    case("FD_2nd_central")
        wave_k_FD1 = dsin(k*dx) / dx
    case("FD_4th_central")
        wave_k_FD1 = (4.0_rk/3.0_rk * dsin(k*dx) - 1.0_rk/6.0_rk * dsin(2.0_rk*k*dx)) / dx
    case("FD_6th_central")
        wave_k_FD1 = (3.0_rk/2.0_rk * dsin(k*dx) - 3.0_rk/10.0_rk * dsin(2.0_rk*k*dx) + 1.0_rk/30.0_rk * dsin(3.0_rk*k*dx)) / dx
    case("FD_4th_central_optimized")
        wave_k_FD1 = 2.0_rk * (0.79926643_rk * dsin(k*dx) - 0.18941314_rk * dsin(2.0_rk*k*dx) + 0.02651995_rk * dsin(3.0_rk*k*dx)) / dx
    case("FD_4th_central_optimized_revised")
        wave_k_FD1 = 2.0_rk * (0.770882380518_rk * dsin(k*dx) - 0.166705904415_rk * dsin(2.0_rk*k*dx) + 0.020843142770_rk * dsin(3.0_rk*k*dx)) / dx
    case default
        call abort(250430, "Discretization unkown "//order_discretization//", I ll walk into the light now." )
    end select
end function

! reduce wavenumbers to FD accuracy for second order derivatives with stencil coeffs c going from -a, ..., a
! general formula derived from canuto p140:
! 1/dx^2 * (c_0 + 2*sum_{j=1}^{a} (c_j * cos(j*k*dx)))
real(kind=rk) function wave_k2_FD2( kx, ky, kz, dx, order_discretization )
    implicit none
    real(kind=rk), intent (in) :: kx, ky, kz, dx(1:3)
    character(len=cshort), intent (in) :: order_discretization
    
    ! ToDo: This is evaluated point-wise, so we should avoid string matching cases
    select case(order_discretization)
    case("FD_2nd_central")
        wave_k2_FD2 = (2.0_rk - 2.0_rk*dcos(kx*dx(1))) / dx(1)**2 + (2.0_rk - 2.0_rk*dcos(ky*dx(2))) / dx(2)**2
        ! if dx is zero then we probably have a 2D case, we just set the return value to 0 for that to proceed
        if (dx(3) /= 0.0_rk) then
            wave_k2_FD2 = wave_k2_FD2 + (2.0_rk - 2.0_rk*dcos(kz*dx(3))) / dx(3)**2
        endif
    case("FD_4th_central")
        wave_k2_FD2 = (30.0_rk - 32.0_rk * dcos(kx*dx(1)) + 2.0_rk * dcos(2.0_rk*kx*dx(1))) / 12.0_rk / dx(1)**2
        wave_k2_FD2 = wave_k2_FD2 + (30.0_rk - 32.0_rk * dcos(ky*dx(2)) + 2.0_rk * dcos(2.0_rk*ky*dx(2))) / 12.0_rk / dx(2)**2
        ! if dx is zero then we probably have a 2D case, we just set the return value to 0 for that to proceed
        if (dx(3) /= 0.0_rk) then
            wave_k2_FD2 = wave_k2_FD2 + (30.0_rk - 32.0_rk * dcos(kz*dx(3)) + 2.0_rk * dcos(2.0_rk*kz*dx(3))) / 12.0_rk / dx(3)**2
        endif
    case("FD_6th_central")
        wave_k2_FD2 = (490.0_rk - 540.0_rk * dcos(kx*dx(1)) + 54.0_rk * dcos(2.0_rk*kx*dx(1)) - 4.0_rk * dcos(3.0_rk*kx*dx(1))) / 180.0_rk / dx(1)**2
        wave_k2_FD2 = wave_k2_FD2 + (490.0_rk - 540.0_rk * dcos(ky*dx(2)) + 54.0_rk * dcos(2.0_rk*ky*dx(2)) - 4.0_rk * dcos(3.0_rk*ky*dx(2))) / 180.0_rk / dx(2)**2
        ! if dx is zero then we probably have a 2D case, we just set the return value to 0 for that to proceed
        if (dx(3) /= 0.0_rk) then
            wave_k2_FD2 = wave_k2_FD2 + (490.0_rk - 540.0_rk * dcos(kz*dx(3)) + 54.0_rk * dcos(2.0_rk*kz*dx(3)) - 4.0_rk * dcos(3.0_rk*kz*dx(3))) / 180.0_rk / dx(3)**2
        endif
    case("FD_8th_central")
        wave_k2_FD2 = (14350.0_rk - 16128.0_rk * dcos(kx*dx(1)) + 2016.0_rk * dcos(2.0_rk*kx*dx(1)) - 256.0_rk * dcos(3.0_rk*kx*dx(1)) + 18.0_rk * dcos(4.0_rk*kx*dx(1))) / 5040.0_rk / dx(1)**2
        wave_k2_FD2 = wave_k2_FD2 + (14350.0_rk - 16128.0_rk * dcos(ky*dx(2)) + 2016.0_rk * dcos(2.0_rk*ky*dx(2)) - 256.0_rk * dcos(3.0_rk*ky*dx(2)) + 18.0_rk * dcos(4.0_rk*ky*dx(2))) / 5040.0_rk / dx(2)**2
        ! if dx is zero then we probably have a 2D case, we just set the return value to 0 for that to proceed
        if (dx(3) /= 0.0_rk) then
            wave_k2_FD2 = wave_k2_FD2 + (14350.0_rk - 16128.0_rk * dcos(kz*dx(3)) + 2016.0_rk * dcos(2.0_rk*kz*dx(3)) - 256.0_rk * dcos(3.0_rk*kz*dx(3)) + 18.0_rk * dcos(4.0_rk*kz*dx(3))) / 5040.0_rk / dx(3)**2
        endif
    case("FD_4th_comp_0_4")
        wave_k2_FD2 = (4490.0_rk - 7104.0_rk * dcos(kx*dx(1)) + 3552.0_rk * dcos(2.0_rk*kx*dx(1)) - 1088.0_rk * dcos(3.0_rk*kx*dx(1)) + 150.0_rk * dcos(4.0_rk*kx*dx(1))) / 144.0_rk / dx(1)**2
        wave_k2_FD2 = wave_k2_FD2 + (4490.0_rk - 7104.0_rk * dcos(ky*dx(2)) + 3552.0_rk * dcos(2.0_rk*ky*dx(2)) - 1088.0_rk * dcos(3.0_rk*ky*dx(2)) + 150.0_rk * dcos(4.0_rk*ky*dx(2))) / 144.0_rk / dx(2)**2
        ! if dx is zero then we probably have a 2D case, we just set the return value to 0 for that to proceed
        if (dx(3) /= 0.0_rk) then
            wave_k2_FD2 = wave_k2_FD2 + (4490.0_rk - 7104.0_rk * dcos(kz*dx(3)) + 3552.0_rk * dcos(2.0_rk*kz*dx(3)) - 1088.0_rk * dcos(3.0_rk*kz*dx(3)) + 150.0_rk * dcos(4.0_rk*kz*dx(3))) / 144.0_rk / dx(3)**2
        endif
    case("FD_4th_comp_1_3")
        wave_k2_FD2 = (470.0_rk - 528.0_rk * dcos(kx*dx(1)) + 48.0_rk * dcos(2.0_rk*kx*dx(1)) + 16.0_rk * dcos(3.0_rk*kx*dx(1)) - 6.0_rk * dcos(4.0_rk*kx*dx(1))) / 144.0_rk / dx(1)**2
        wave_k2_FD2 = wave_k2_FD2 + (470.0_rk - 528.0_rk * dcos(ky*dx(2)) + 48.0_rk * dcos(2.0_rk*ky*dx(2)) + 16.0_rk * dcos(3.0_rk*ky*dx(2)) - 6.0_rk * dcos(4.0_rk*ky*dx(2))) / 144.0_rk / dx(2)**2
        ! if dx is zero then we probably have a 2D case, we just set the return value to 0 for that to proceed
        if (dx(3) /= 0.0_rk) then
            wave_k2_FD2 = wave_k2_FD2 + (470.0_rk - 528.0_rk * dcos(kz*dx(3)) + 48.0_rk * dcos(2.0_rk*kz*dx(3)) + 16.0_rk * dcos(3.0_rk*kz*dx(3)) - 6.0_rk * dcos(4.0_rk*kz*dx(3))) / 144.0_rk / dx(3)**2
        endif
    case("FD_4th_comp_2_2")
        wave_k2_FD2 = (130.0_rk - 32.0_rk * dcos(kx*dx(1)) - 128.0_rk * dcos(2.0_rk*kx*dx(1)) + 32.0_rk * dcos(3.0_rk*kx*dx(1)) - 2.0_rk * dcos(4.0_rk*kx*dx(1))) / 144.0_rk / dx(1)**2
        wave_k2_FD2 = wave_k2_FD2 + (130.0_rk - 32.0_rk * dcos(ky*dx(2)) - 128.0_rk * dcos(2.0_rk*ky*dx(2)) + 32.0_rk * dcos(3.0_rk*ky*dx(2)) - 2.0_rk * dcos(4.0_rk*ky*dx(2))) / 144.0_rk / dx(2)**2
        ! if dx is zero then we probably have a 2D case, we just set the return value to 0 for that to proceed
        if (dx(3) /= 0.0_rk) then
            wave_k2_FD2 = wave_k2_FD2 + (130.0_rk - 32.0_rk * dcos(kz*dx(3)) - 128.0_rk * dcos(2.0_rk*kz*dx(3)) + 32.0_rk * dcos(3.0_rk*kz*dx(3)) - 2.0_rk * dcos(4.0_rk*kz*dx(3))) / 144.0_rk / dx(3)**2
        endif
    case("FD_6th_comp_2_4")
        wave_k2_FD2 = (9170.0_rk - 9312.0_rk * dcos(kx*dx(1)) - 540.0_rk * dcos(2.0_rk*kx*dx(1)) + 1040.0_rk * dcos(3.0_rk*kx*dx(1)) - 434.0_rk * dcos(4.0_rk*kx*dx(1)) + 80.0_rk * dcos(5.0_rk*kx*dx(1)) - 4.0_rk * dcos(6.0_rk*kx*dx(1))) / 3600.0_rk / dx(1)**2
        wave_k2_FD2 = wave_k2_FD2 + (9170.0_rk - 9312.0_rk * dcos(ky*dx(2)) - 540.0_rk * dcos(2.0_rk*ky*dx(2)) + 1040.0_rk * dcos(3.0_rk*ky*dx(2)) - 434.0_rk * dcos(4.0_rk*ky*dx(2)) + 80.0_rk * dcos(5.0_rk*ky*dx(2)) - 4.0_rk * dcos(6.0_rk*ky*dx(2))) / 3600.0_rk / dx(2)**2
        ! if dx is zero then we probably have a 2D case, we just set the return value to 0 for that to proceed
        if (dx(3) /= 0.0_rk) then
            wave_k2_FD2 = wave_k2_FD2 + (9170.0_rk - 9312.0_rk * dcos(kz*dx(3)) - 540.0_rk * dcos(2.0_rk*kz*dx(3)) + 1040.0_rk * dcos(3.0_rk*kz*dx(3)) - 434.0_rk * dcos(4.0_rk*kz*dx(3)) + 80.0_rk * dcos(5.0_rk*kz*dx(3)) - 4.0_rk * dcos(6.0_rk*kz*dx(3))) / 3600.0_rk / dx(3)**2
        endif
    case("FD_6th_mehrstellen")
        ! MST assumes that dx is the same in all directions
        ! if dx is zero then we probably have a 2D case, we just set the return value to 0 for that to proceed
        if (dx(3) == 0.0_rk) then
            wave_k2_FD2 = (20.0_rk - 8.0_rk * (dcos(kx*dx(1)) + dcos(ky*dx(2))) - 4.0_rk * (dcos(kx*dx(1))*dcos(ky*dx(2)))) / 6.0_rk / dx(1)**2
        else
            wave_k2_FD2 = ( -128.0_rk &
                + 28.0_rk * (dcos(kx * dx(1)) + dcos(ky * dx(2)) + dcos(kz * dx(3))) &
                + 12.0_rk * (dcos(kx * dx(1)) * dcos(ky * dx(2)) &
                            + dcos(kx * dx(1)) * dcos(kz * dx(3)) &
                            + dcos(ky * dx(2)) * dcos(kz * dx(3))) &
                +  8.0_rk * (dcos(kx * dx(1)) * dcos(ky * dx(2)) * dcos(kz * dx(3)))) / 30.0_rk / dx(1)**2
        endif
    case default
        call abort(250430, "Discretization unkown for FFT FD spectral symbol: "//trim(adjustl(order_discretization)))
    end select
end function

end module module_fft
