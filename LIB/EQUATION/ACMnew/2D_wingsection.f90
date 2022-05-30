subroutine draw_2d_wingsections(time, mask, x0, dx, Bs, g )
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
    integer(kind=ik) :: ix, iy, i, nwings
    real(kind=rk) :: x00, y00, alpha, alpha_dt, u00, v00, time00=0.0_rk
    logical, save :: initialized=.false.

    if (size(mask,1) /= Bs(1)+2*g .or. size(mask,2) /= Bs(2)+2*g ) then
        call abort(777107, "mask: wrong array size, there's pirates, captain!")
    endif

    ! reset mask block
    mask = 0.0_rk
    nwings = 2

    do i = 1, nwings
        ! if (abs(time-wingsections(i)%time) > 1.0e-7 .or. .not. initialized) then
            ! evaluate fourier series for all parameters
            if (wingsections(i)%nfft_x0 > 0) then
                call fseries_eval( time, x00, u00, wingsections(i)%a0_x0, wingsections(i)%ai_x0, wingsections(i)%bi_x0)
            else
                x00 = wingsections(i)%a0_x0 / 2.0_rk
                u00 = 0.0_rk
            endif

            if (wingsections(i)%nfft_y0 > 0) then
                call fseries_eval( time, y00, v00, wingsections(i)%a0_y0, wingsections(i)%ai_y0, wingsections(i)%bi_y0)
            else
                y00 = wingsections(i)%a0_y0 / 2.0_rk
                v00 = 0.0_rk
            endif

            if (wingsections(i)%nfft_alpha > 0) then
                call fseries_eval( time, alpha, alpha_dt, wingsections(i)%a0_alpha, wingsections(i)%ai_alpha, wingsections(i)%bi_alpha)
            else
                alpha = wingsections(i)%a0_alpha / 2.0_rk
                alpha_dt = 0.0_rk
            endif
            ! input is in deg, convert to rad
            alpha    = alpha*pi/180.0_rk
            alpha_dt = alpha_dt*pi/180.0_rk
            ! skip this evaluation of fseries for following blocks (until time changes again)
            wingsections(i)%time = time
        ! endif

        ! angular velocity (rad/T)
        omega = 2.0_rk*pi * alpha_dt

        ! parameter for smoothing function (width of smoothing layer)
        dx_min = maxval( (2.0_rk**(-params_acm%Jmax))*params_acm%domain_size(1:2)/real(Bs(1:2)-1,kind=rk) )
        h = params_acm%C_smooth * dx_min

        ! Note: this basic mask function is set on the ghost nodes as well.
        do iy = g+1, Bs(2)+g+ONE_SKIPREDUNDANT
            y = dble(iy-(g+1)) * dx(2) + x0(2) - y00

            do ix = g+1, Bs(1)+g+ONE_SKIPREDUNDANT
                x = dble(ix-(g+1)) * dx(1) + x0(1) - x00

                xp =  cos(alpha)*x + sin(alpha)*y
                yp = -sin(alpha)*x + cos(alpha)*y

                u = u00 - omega*y
                v = v00 + omega*x

                tmp = smoothstep(abs(xp-0.5_rk), 0.5_rk, h)
                tmp = tmp*smoothstep(abs(yp), 0.5_rk*wingsections(i)%section_thickness, h)

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
    end do ! loop over wing sections

    initialized = .true.
end subroutine



subroutine init_wingsection_from_file(file, theSection, time)
    use module_params

    implicit none

    character(len=80), intent(in) :: file
    real(kind=rk), intent(in) :: time
    type(wingsection), intent(inout) :: theSection
    type(inifile) :: ifile

    call read_ini_file_mpi(ifile, file, .true.  )

    call read_param_mpi(ifile, "Wingsection", "type", theSection%kinematics_type, "Fourier")
    if (theSection%kinematics_type /= "Fourier") then
        call abort(2105101, "Only fourier parameterization supported")
    endif

    call read_param_mpi(ifile, "Wingsection", "section_thickness", theSection%section_thickness, 0.05_rk)

    call read_param_mpi(ifile, "Wingsection", "nfft_x0"  , theSection%nfft_x0, 0)
    call read_param_mpi(ifile, "Wingsection", "nfft_y0"  , theSection%nfft_y0, 0)
    call read_param_mpi(ifile, "Wingsection", "nfft_alpha", theSection%nfft_alpha, 0)

    allocate(theSection%ai_x0(1:theSection%nfft_x0))
    allocate(theSection%bi_x0(1:theSection%nfft_x0))
    call read_param_mpi(ifile, "Wingsection", "a0_x0", theSection%a0_x0, 2.0_rk*0.5_rk*params_acm%domain_size(1))
    call read_param_mpi(ifile, "Wingsection", "ai_x0", theSection%ai_x0)
    call read_param_mpi(ifile, "Wingsection", "bi_x0", theSection%bi_x0)

    allocate(theSection%ai_y0(1:theSection%nfft_y0))
    allocate(theSection%bi_y0(1:theSection%nfft_y0))
    call read_param_mpi(ifile, "Wingsection", "a0_y0", theSection%a0_y0, 2.0_rk*0.5_rk*params_acm%domain_size(2))
    call read_param_mpi(ifile, "Wingsection", "ai_y0", theSection%ai_y0)
    call read_param_mpi(ifile, "Wingsection", "bi_y0", theSection%bi_y0)

    allocate(theSection%ai_alpha(1:theSection%nfft_alpha))
    allocate(theSection%bi_alpha(1:theSection%nfft_alpha))
    call read_param_mpi(ifile, "Wingsection", "a0_alpha", theSection%a0_alpha, 0.0_rk)
    call read_param_mpi(ifile, "Wingsection", "ai_alpha", theSection%ai_alpha)
    call read_param_mpi(ifile, "Wingsection", "bi_alpha", theSection%bi_alpha)

    theSection%time = time
    theSection%initialized = .true.
end subroutine
