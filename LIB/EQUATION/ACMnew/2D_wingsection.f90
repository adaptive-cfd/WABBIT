subroutine draw_2d_wingsection(time, mask, x0, dx, Bs, g )

    use module_params
    use module_precision

    implicit none

    real(kind=rk), intent(in) :: time
    integer(kind=ik), intent(in) :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    !> mask term for every grid point of this block
    real(kind=rk), dimension(:,:,:), intent(out)     :: mask
    !> spacing and origin of block
    real(kind=rk), dimension(2), intent(in) :: x0, dx

    ! auxiliary variables
    real(kind=rk)  :: x, y, R, h, dx_min, tmp, xp, yp, omega, u, v
    ! loop variables
    integer(kind=ik) :: ix, iy

    type(inifile) :: ifile
    logical, save :: initialized = .false.
    real(kind=rk), save :: x00, y00, alpha, alpha_dt, u00, v00, time00=0.0_rk, section_thickness
    real(kind=rk), allocatable, save :: ai_x0(:), bi_x0(:)
    real(kind=rk), allocatable, save :: ai_y0(:), bi_y0(:)
    real(kind=rk), allocatable, save :: ai_alpha(:), bi_alpha(:)
    real(kind=rk), save :: a0_x0
    real(kind=rk), save :: a0_y0
    real(kind=rk), save :: a0_alpha
    integer(kind=ik), save :: nfft_x0, nfft_y0, nfft_alpha
    character(len=clong), save :: kinematics_type

    if (size(mask,1) /= Bs(1)+2*g .or. size(mask,2) /= Bs(2)+2*g ) then
        call abort(777107, "mask: wrong array size, there's pirates, captain!")
    endif


    if (.not. initialized) then
        call read_ini_file_mpi(ifile, params_acm%wingsection_inifile, .true.  )

        call read_param_mpi(ifile, "Wingsection", "type", kinematics_type, "Fourier")
        if (kinematics_type /= "Fourier") then
            call abort(2105101, "Only fourier parameterization supported")
        endif

        call read_param_mpi(ifile, "Wingsection", "section_thickness", section_thickness, 0.05_rk)

        call read_param_mpi(ifile, "Wingsection", "nfft_x0", nfft_x0, 0)
        call read_param_mpi(ifile, "Wingsection", "nfft_y0", nfft_y0, 0)
        call read_param_mpi(ifile, "Wingsection", "nfft_alpha", nfft_alpha, 0)

        allocate(ai_x0(1:nfft_x0))
        allocate(bi_x0(1:nfft_x0))
        call read_param_mpi(ifile, "Wingsection", "a0_x0", a0_x0, 2.0_rk*0.5_rk*params_acm%domain_size(1))
        call read_param_mpi(ifile, "Wingsection", "ai_x0", ai_x0)
        call read_param_mpi(ifile, "Wingsection", "bi_x0", bi_x0)

        allocate(ai_y0(1:nfft_y0))
        allocate(bi_y0(1:nfft_y0))
        call read_param_mpi(ifile, "Wingsection", "a0_y0", a0_y0, 2.0_rk*0.5_rk*params_acm%domain_size(2))
        call read_param_mpi(ifile, "Wingsection", "ai_y0", ai_y0)
        call read_param_mpi(ifile, "Wingsection", "bi_y0", bi_y0)

        allocate(ai_alpha(1:nfft_alpha))
        allocate(bi_alpha(1:nfft_alpha))
        call read_param_mpi(ifile, "Wingsection", "a0_alpha", a0_alpha, 0.0_rk)
        call read_param_mpi(ifile, "Wingsection", "ai_alpha", ai_alpha)
        call read_param_mpi(ifile, "Wingsection", "bi_alpha", bi_alpha)
    endif

    ! if (abs(time-time00) > 1.0e-12 .or. .not. initialized) then
        ! evaluate fourier series for all parameters
        if (nfft_x0 > 0) then
            call fseries_eval( time, x00, u00, a0_x0, ai_x0, bi_x0)
        else
            x00 = a0_x0 / 2.0_rk
            u00 = 0.0_rk
        endif

        if (nfft_y0 > 0) then
            call fseries_eval( time, y00, v00, a0_y0, ai_y0, bi_y0)
        else
            y00 = a0_y0 / 2.0_rk
            v00 = 0.0_rk
        endif

        if (nfft_alpha > 0) then
            call fseries_eval( time, alpha, alpha_dt, a0_alpha, ai_alpha, bi_alpha)
        else
            alpha = a0_x0 / 2.0_rk
            alpha_dt = 0.0_rk
        endif
        ! input is in deg
        alpha    = alpha*pi/180.0_rk
        alpha_dt = alpha_dt*pi/180.0_rk
        ! skip this evaluation of fseries for following blocks (until time changes again)
        time00 = time
    ! endif

    initialized = .true.

    ! angular velocity (rad/T)
    omega = 2.0_rk*pi * alpha_dt

    ! reset mask block
    mask = 0.0_rk

    ! parameter for smoothing function (width of smoothing layer)
    dx_min = maxval( (2.0_rk**(-params_acm%Jmax))*params_acm%domain_size(1:2)/real(Bs(1:2)-1,kind=rk) )
    h = params_acm%C_smooth * dx_min

    ! Note: this basic mask function is set on the ghost nodes as well.
    do iy = g+1, Bs(2)+g
        y = dble(iy-(g+1)) * dx(2) + x0(2) - y00

        do ix = g+1, Bs(1)+g
            x = dble(ix-(g+1)) * dx(1) + x0(1) - x00

            xp =  cos(alpha)*x + sin(alpha)*y
            yp = -sin(alpha)*x + cos(alpha)*y

            u = u00 - omega*y
            v = v00 + omega*x

            tmp = smoothstep(abs(xp+0.5_rk), 0.5_rk, h)
            tmp = tmp*smoothstep(abs(yp), 0.5_rk*section_thickness, h)

            if (tmp >= mask(ix,iy,1)) then
                ! mask function
                mask(ix,iy,1) = tmp
                ! solid velocity
                mask(ix,iy,2) = u
                mask(ix,iy,3) = v
                ! color
                mask(ix,iy,5) = 1.0_rk
            endif
        end do
    end do


end subroutine
