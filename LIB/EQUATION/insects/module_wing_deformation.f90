module module_wing_deformation
    use module_globals
    use module_ini_files_parser_mpi
    use generic_functions_module

    implicit none
    ! A wise programmer once told me:
    ! ""I usually find it helpful to use the private keyword by itself initially, which specifies
    ! that everything within the module is private unless explicitly marked public.""
    PRIVATE

    PUBLIC :: add_deformation_velocity, check_deformation_velocity, break_down_into_simple_model, material2wing, wing2material
    PUBLIC :: deformation_model


    type deformation_model
        ! Jonas: Initialise variables for the deformed wing. These Variables are given in a deformation file under [Deformation]

        ! I might want to edit the way I calculate the deformation later, so I'm trying to make this as modular as I can

        character(len=clong) :: model_name
        character(len=clong) :: filename
        character(len=clong) :: units

        ! variables for a model I'm going to call "lin_sin": a spatially linear increase in the deformation angle D_alpha, multiplied by some sin - cos sum

        real(kind=rk) :: D_alpha_slope = 0.0_rk
        real(kind=rk) :: D_alpha_0 = 0.0_rk
        real(kind=rk) :: omega = 0.0_rk

        ! define the axis of the flexion line. Wing deformation occurs around flexion lines, where veins stiffen the membrane (Hsu 2024). 
        ! we define this as the angle in rad between y^(mat)_(right wing) and the flexion line in question, as well as defining a pivot point through which the flexion line passes.
        ! for bluebottles, this corresponds to the claval flexion line (Hsu 2024).  

        real(kind=rk) :: flexion_angle
        real(kind=rk) :: pivot_flexion(1:3)


        ! maybe program a fourier series later

        real(kind=rk) :: ai_D_alpha = 0.0_rk
        real(kind=rk) :: bi_D_alpha = 0.0_rk

        LOGICAL :: deinsect_verbose = .true.

        ! end goal: describe using an arbitrary (continuous) function of time and y coordinate

        ! the deformation angle D_alpha is treated as a separable function of time and space. the time function is A_tilde, the angle distribution is D_alpha_tilde

        type(generic_function1D) :: A_tilde
        type(generic_function1D) :: D_alpha_tilde

        ! a twist offset is problematic when we have an fseries, since the function is then no longer separable. 
        ! to mitigate this, we store the twist offset as a separate array of fourier coefficients, which is added to alpha (feathering angle) in the main flapping motion subroutine.
        ! this also allows us to separate the twisting into pitching of the wing and deformation, and is helpful for exmperimental measurements which describe the twist with an offset.

        real(kind=rk) :: twistOffset_ai(10)
        real(kind=rk) :: twistOffset_bi(10)
        real(kind=rk) :: twistOffset_a0 = 0.0_rk

    end type


    contains

    subroutine break_down_into_simple_model(model)
        ! function to turn derived models into basic models. you only need to call this once when reading in parameters. This also converts those pesky degrees into radians

        type(deformation_model), INTENT(INOUT) :: model
        integer :: i
        real(kind=rk) :: alpha_tip
        real(kind=rk) :: alpha_root
        real(kind=rk) :: R = 1.0_rk
        model%flexion_angle = 3.141592654_rk * model%flexion_angle / 180.0_rk

        select case (model%model_name)

            case ("linearAlpha")

                ! work on the spatial distribution first
                
                alpha_tip = model%D_alpha_tilde%coeffs(2)
                alpha_root = model%D_alpha_tilde%coeffs(1)

                if (model%units == "deg") then
                ! rescale so we can work with degrees

                model%D_alpha_tilde%coeffs(2) = 2 * 3.141592654_rk * ((alpha_tip - alpha_root) / R) / 360.0_rk
                model%D_alpha_tilde%coeffs(1) = 2.0_rk * 3.141592654_rk * (alpha_root / R) / 360.0_rk

                end if


                ! make sure no higher order terms sneak into our spatial distribution (the keyword is LINEAR alpha)

                do i = 3, size(model%D_alpha_tilde%coeffs)
                    model%D_alpha_tilde%coeffs(i) = 0.0_rk
                end do

                model%D_alpha_tilde%function_name = "polynomial"

                if ((model%A_tilde%function_name == "sin") .or. (model%A_tilde%function_name == "square")) then

                    ! total range of amplitudes, also accounting for when the wing twists negatively

                    model%A_tilde%coeffs(1) = 1.0_rk + model%A_tilde%coeffs(4)

                    ! calculate wing reversals per wingbeat

                    model%A_tilde%coeffs(2) = 2.0_rk * 3.141592654_rk * model%A_tilde%coeffs(2)
                    ! calculate reversal phasemodel

                    model%A_tilde%coeffs(3) = 2.0_rk * 3.141592654_rk * model%A_tilde%coeffs(3) / 360.0_rk

                    ! accout for negative twist offset

                    model%A_tilde%coeffs(4) = 1.0_rk - model%A_tilde%coeffs(4)
                endif
                model%model_name = "generic_angle_distribution"

            case("WalkerTaylor")
                ! Import deformation measurements according to 
                ! Walker, S. M., & Taylor, G. K. (n.d.). A semi-empirical model of the aerodynamics of manoeuvring insect flight.
                ! for now, we only consider the periodic components of the local twist. The offset could be added to the feathering angle.
                ! 
                model%D_alpha_tilde%function_name = "polynomial"

                if (model%units == "deg") then
                ! the amplitude is set by the fourier coefficients. note the deg2rad conversion in the polynomial!
                model%D_alpha_tilde%coeffs = (/0.0_rk, 3.141592654_rk/180.0_rk, 0.0_rk, 0.0_rk, 0.0_rk/)

                end if
                model%A_tilde%function_name = "fseries"
                model%model_name = "generic_angle_distribution"

            case("generic_angle_distribution", "lin_sin")

                ! don't try to translate these models

            case default
                call abort(23673, "model not selected!")


        endselect



    end subroutine

    subroutine wing2material(time, x_mat, model, side)

        ! Map material coordinates inherent to the wing to their position in wing coordinates.
        ! x_mat^(mat) -> x_mat^(w)

        type(deformation_model), intent(in) :: model
        real(kind=rk), intent(in) :: time
        real(kind=rk), intent(inout) :: x_mat(1:3)
        character(len=1), intent(in) :: side

        real(kind=rk) :: A_tilde, D_alpha_tilde, dA_tilde_dt, dD_alpha_tilde_dy
        real(kind=rk) :: D_alpha
        real(kind=rk) :: normalisation
        real(kind=rk) :: x_w, y_w, z_w, x_n, y_n, z_n

        real(kind=rk) ::  x_c_mat_w(1:3), x_wing(1:3), e_n_mat_w(1:3), e_x_mat_w(1:3), e_y_mat_w(1:3), de_n_mat_w_dt(1:3), de_x_mat_w_dt(1:3), de_y_mat_w_dt(1:3), R_y(1:3, 1:3), dR_y_dalpha(1:3, 1:3), x_tmp(1:3), A_tmp(1:3, 1:3)
        ! is this wing on the left or right side ?
        select case (model%model_name)
            case("generic_angle_distribution")

                x_tmp = x_mat

                ! translate these coordinates to the pivot of the flexion line, so we can align our rotation axis with the primary flexion line.

                x_wing(1) = x_tmp(1) - model%pivot_flexion(1)
                x_wing(2) = x_tmp(2) - model%pivot_flexion(2)
                x_wing(3) = x_tmp(3) - model%pivot_flexion(3)

                x_tmp = x_wing

                ! align the twisting axis with the flexion line
                x_wing(1) = cos(model%flexion_angle) * x_tmp(1) - sin(model%flexion_angle) * x_tmp(2)
                x_wing(2) = sin(model%flexion_angle) * x_tmp(1) + cos(model%flexion_angle) * x_tmp(2)

                ! more generic case. Here, we assume that the angle distribution follows a prescribed function with a derivative along the y(mat) axis and a time derivative.
                ! We define a material coordinate inherent to the wing, which describes a material point independent of deformation. 
                ! This falls together with the wing coordinate in the rigid case.
                ! x_c_mat_w denotes the points that lie in the centreplane of the wing's material coordinates.

                x_c_mat_w = x_wing
                x_c_mat_w(3) = 0

                ! calculate the deformation angle.

                D_alpha_tilde = x_wing(2)
                A_tilde = time

                ! calculate derivatives

                dD_alpha_tilde_dy = x_wing(2)
                dA_tilde_dt = time

                ! we call the time function and the spatial function, assuming that the angle distribution is a separable function of time and space

                call eval(model%A_tilde, A_tilde, 0) 
                call eval(model%D_alpha_tilde, D_alpha_tilde, 0) 
                call eval(model%A_tilde, dA_tilde_dt, 1) 
                call eval(model%D_alpha_tilde, dD_alpha_tilde_dy, 1) 

                ! depending on the side, invert the deformation angle (we define the angle distribution based on the right wing!)

                if (side == "R") then
                    D_alpha_tilde = -D_alpha_tilde
                    dD_alpha_tilde_dy = -dD_alpha_tilde_dy
                endif

                D_alpha = D_alpha_tilde * A_tilde

                ! factor to rescale our normal to 1

                normalisation = 1 / sqrt((dD_alpha_tilde_dy * A_tilde) ** 2 * x_wing(1) ** 2 + 1)

                ! row wise calculation of the surface normal

                e_n_mat_w(1) = normalisation * sin(D_alpha)
                e_n_mat_w(2) = normalisation * (dD_alpha_tilde_dy * A_tilde) * x_wing(1)
                e_n_mat_w(3) = - normalisation * cos(D_alpha)

                ! section-wise rotation around the flexion axis of the centreplane of the material coordinates
                
                x_c_mat_w(1) = cos(D_alpha)*x_wing(1)
                x_c_mat_w(2) = x_wing(2)
                x_c_mat_w(3) = sin(D_alpha)*x_wing(1)

                ! By definition, the z coordinate is normal to the material x, y plane. Extrude normally to the centreplane, but describe this extrusion in the wing system

                x_wing(1) = x_c_mat_w(1) - x_wing(3) * e_n_mat_w(1) 
                x_wing(2) = x_c_mat_w(2) - x_wing(3) * e_n_mat_w(2) 
                x_wing(3) = x_c_mat_w(3) - x_wing(3) * e_n_mat_w(3) 

                x_tmp = x_wing

                ! align the twisting axis with the flexion line

                x_wing(1) = cos(model%flexion_angle) * x_tmp(1) + sin(model%flexion_angle) * x_tmp(2)
                x_wing(2) = - sin(model%flexion_angle) * x_tmp(1) + cos(model%flexion_angle) * x_tmp(2)

                ! translate these coordinates to the pivot of the flexion line, so we can align our rotation axis with the primary flexion line

                x_wing(1) = x_wing(1) + model%pivot_flexion(1)
                x_wing(2) = x_wing(2) + model%pivot_flexion(2)
                x_wing(3) = x_wing(3) + model%pivot_flexion(3)
                x_mat = x_wing
            case default
                call abort(23673, "model not selected!")
        end select
    end subroutine wing2material

    subroutine material2wing(time, x_mat, deformation, side)

        ! Map material coordinates inherent to the wing to their position in wing coordinates. This is expensive and not exact, which is why we only use this to approximate the bonding box corners!
        ! x_mat^(mat) -> x_mat^(w)

        implicit none
        real(kind=rk), intent(in)  :: time
        real(kind=rk), intent(inout) :: x_mat(1:3)
        real(kind=rk) :: x_w(1:3)
        ! use the same type and intent as in your wing2material interface
        type(deformation_model), intent(in) :: deformation
        character(len=1), intent(in) :: side
        real(kind=rk) :: xt(1:3), res(1:3)
        integer :: it

        x_w = x_mat
        do it = 1, 20
            xt = x_w
            call wing2material(time, xt, deformation, side)
            res = xt - x_mat
            if (maxval(abs(res)) < 1.0e-12_rk) exit
            x_w = x_w - res
        enddo
        x_mat = x_w
    end subroutine material2wing


    subroutine add_deformation_velocity(time, x_wing, u_s_w, model, side)
        type(deformation_model), intent(in) :: model
        real(kind=rk), intent(in) :: time
        real(kind=rk), intent(inout) :: u_s_w(1:3)
        real(kind=rk), intent(in) :: x_wing(1:3)
        character(len=1), intent(in) :: side

        real(kind=rk) :: A_tilde
        real(kind=rk) :: D_alpha, D_alphadot,A_tildedot,normalisation,normalisationdot, D_alpha_tilde, dA_tilde_dt, dD_alpha_tilde_dy
        real(kind=rk) :: x_w, y_w, z_w, x_n, y_n, z_n
        real(kind=rk) :: u_midplane(1:3), u_displaced(1:3)
        real(kind=rk) ::  x_c_mat_w(1:3), x_mat(1:3), e_n_mat_w(1:3), e_x_mat_w(1:3), e_y_mat_w(1:3), de_n_mat_w_dt(1:3), de_x_mat_w_dt(1:3), de_y_mat_w_dt(1:3), R_y(1:3, 1:3), dR_y_dalpha(1:3, 1:3), x_tmp(1:3), A_tmp(1:3, 1:3)

        select case (model%model_name)
            
            case("generic_angle_distribution")


                ! more generic case. Here, we assume that the angle distribution follows a prescribed, separable function with a derivative along the y(mat) axis and a time derivative.
                ! We define a material coordinate inherent to the wing, which describes a material point independent of deformation. This falls together with the wing coordinate in the rigid case.

                x_mat = x_wing

                ! translate these coordinates to the pivot of the flexion line, so we can align our rotation axis with the primary flexion line

                x_mat(1) = x_mat(1) - model%pivot_flexion(1)
                x_mat(2) = x_mat(2) - model%pivot_flexion(2)
                x_mat(3) = x_mat(3) - model%pivot_flexion(3)

                x_tmp= x_mat

                ! align the twisting axis with the flexion line
                x_mat(1) = cos(model%flexion_angle) * x_tmp(1) - sin(model%flexion_angle) * x_tmp(2)
                x_mat(2) = sin(model%flexion_angle) * x_tmp(1) + cos(model%flexion_angle) * x_tmp(2)


                x_c_mat_w = x_mat
                x_c_mat_w(3) = 0


                ! calculate the deformation angle
                D_alpha_tilde = x_mat(2)
                A_tilde = time

                ! calculate derivatives

                dD_alpha_tilde_dy = x_mat(2)
                dA_tilde_dt = time

                ! we call the time function and the spatial function, assuming that the angle distribution is a separable function of time and space

                call eval(model%A_tilde, A_tilde, 0) 
                call eval(model%D_alpha_tilde, D_alpha_tilde, 0) 
                call eval(model%A_tilde, dA_tilde_dt, 1) 
                call eval(model%D_alpha_tilde, dD_alpha_tilde_dy, 1) 

                if (side=="R") then 
                    dD_alpha_tilde_dy = -dD_alpha_tilde_dy
                    D_alpha_tilde = -D_alpha_tilde
                endif

                D_alpha = D_alpha_tilde * A_tilde
                ! factor to rescale our normal to 1
                normalisation = 1.0_rk / sqrt( &
                    1.0_rk + (x_c_mat_w(1) * dD_alpha_tilde_dy * A_tilde)**2 )
                normalisationdot = -normalisation**3 * &
                   x_c_mat_w(1)**2 * &
                   (dD_alpha_tilde_dy * A_tilde) * &
                   (dD_alpha_tilde_dy * dA_tilde_dt)

                ! row wise calculation of the time derivate of the surface normal

                de_n_mat_w_dt(1) = normalisationdot * sin(D_alpha) &
                                    + normalisation * cos(D_alpha) * &
                                    D_alpha_tilde * dA_tilde_dt
                de_n_mat_w_dt(2) = normalisationdot * x_c_mat_w(1) * &
                                        dD_alpha_tilde_dy * A_tilde &
                                    + normalisation * x_c_mat_w(1) * &
                                        dD_alpha_tilde_dy * dA_tilde_dt
                de_n_mat_w_dt(3) = - normalisationdot * cos(D_alpha) &
                    + normalisation * sin(D_alpha) * &
                    D_alpha_tilde * dA_tilde_dt

                ! finally we can calculate the velocity. start by calculating the velocity of the material centre plane and then add the velocity component through the distance normal to the surface

                ! Calculate the velocity of the material centre plane
                ! in the flexion-line frame.

                u_midplane(1) = -x_c_mat_w(1) * sin(D_alpha) * &
                        D_alpha_tilde * dA_tilde_dt

                u_midplane(2) = 0.0_rk

                u_midplane(3) = x_c_mat_w(1) * cos(D_alpha) * &
                        D_alpha_tilde * dA_tilde_dt

                ! Add the velocity component due to the changing normal.

                u_displaced(1) = u_midplane(1) - &
                        x_mat(3) * de_n_mat_w_dt(1)

                u_displaced(2) = u_midplane(2) - &
                        x_mat(3) * de_n_mat_w_dt(2)

                u_displaced(3) = u_midplane(3) - &
                        x_mat(3) * de_n_mat_w_dt(3)

                ! Rotate the deformation velocity back into the wing frame.
                ! Do not rotate the existing u_s_w contribution again.

                x_tmp = u_displaced

                u_displaced(1) = cos(model%flexion_angle) * x_tmp(1) &
                        + sin(model%flexion_angle) * x_tmp(2)

                u_displaced(2) = -sin(model%flexion_angle) * x_tmp(1) &
                        + cos(model%flexion_angle) * x_tmp(2)

                u_displaced(3) = x_tmp(3)

                ! Add deformation velocity to the existing wing velocity.

                u_s_w = u_s_w + u_displaced
                                
                end select


    end subroutine add_deformation_velocity


    subroutine check_deformation_velocity(x_tip_w, us_tip_w, time,M_g2b,M_b2w,x_pivot_b, &
        xc_body_g, vc_body_g, rot_rel_wing_w,rot_body_b, model, side, rigid)
        real(kind=rk), intent(in) :: time
        real(kind=rk),intent(in)::M_g2b(1:3,1:3),M_b2w(1:3,1:3),x_pivot_b(1:3),rot_rel_wing_w(1:3),rot_body_b(1:3), xc_body_g(3), vc_body_g(3)
        type(deformation_model), intent(in) :: model
        character(len=1), intent(in) :: side
        real(kind=rk),intent(inout)::x_tip_w(3), us_tip_w(3)
        real(kind=rk) :: x_tip_mat(3), v_tmp(3)
        logical, INTENT(IN):: rigid

        v_tmp = 0.0_rk
        x_tip_mat = x_tip_w

        ! debugging purposes: Deform the wing coordinates and calculate the deformation velocity. 

        if (.not. rigid) then
            call wing2material(time, x_tip_mat, model, side)
        end if
        
        if (.not. rigid) then
            call add_deformation_velocity(time, x_tip_w, v_tmp, model, side)
        end if
        us_tip_w = v_tmp
        x_tip_w =x_tip_mat
    end subroutine check_deformation_velocity

end module module_wing_deformation