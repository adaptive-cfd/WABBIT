subroutine initialise_wing_deformation(PARAMS, SECTION, Insect)


    type(inifile), INTENT(IN) :: PARAMS
    character(len=*), intent(in) :: SECTION
    ! the actual insect to be initialized
    type(diptera), intent(inout) :: Insect

    real(kind=rk) :: defaultvec(1:3), defaultvec5(5), alpha_tip, alpha_root, reversalPhase, NReversals,nt_frac
    integer(kind=2) :: wingID

        ! which wings are rigid?

    call read_param_mpi(PARAMS, SECTION, "LeftWingRigid", Insect%Wings(1)%rigid, .true.)
    call read_param_mpi(PARAMS, SECTION, "RightWingRigid", Insect%Wings(2)%rigid, .true.)
    call read_param_mpi(PARAMS, SECTION, "LeftWing2Rigid", Insect%Wings(3)%rigid, .true.)
    call read_param_mpi(PARAMS, SECTION, "RightWing2Rigid", Insect%Wings(4)%rigid, .true.)

    do wingID=1,4
        call read_param_mpi(PARAMS, SECTION, 'modelName', Insect%Wings(wingID)%deformation%model_name, "lin_sin")
        ! call read_param_mpi(PARAMS, SECTION, 'offsetDalpha', Insect%Wings(wingID)%deformation%D_alpha_0, 0.0_rk)
        ! call read_param_mpi(PARAMS, SECTION, 'slopeDalpha', Insect%Wings(wingID)%deformation%D_alpha_slope, 0.0_rk)
        ! call read_param_mpi(PARAMS, SECTION, 'ai', Insect%Wings(wingID)%deformation%ai_D_alpha, 0.0_rk)
        ! call read_param_mpi(PARAMS, SECTION, 'bi', Insect%Wings(wingID)%deformation%bi_D_alpha, 0.0_rk)
        ! call read_param_mpi(PARAMS, SECTION, 'omega', Insect%Wings(wingID)%deformation%omega, 0.0_rk)


        ! call read_param_mpi(PARAMS, SECTION, 'timeCoefficients', defaultvec5, (/0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk/))
        ! call read_param_mpi(PARAMS, SECTION, 'timeFunction', Insect%Wings(wingID)%deformation%A_tilde%function_name, "sin")
        ! !search for parameters for the linearAlpha model (temporal)
        ! call read_param_mpi(PARAMS, SECTION, 'NReversals', NReversals, -1e5_rk)
        ! call read_param_mpi(PARAMS, SECTION, 'reversalPhase', reversalPhase, -1e5_rk)

        ! Insect%Wings(wingID)%deformation%A_tilde%coeffs = defaultvec5

        ! call read_param_mpi(PARAMS, SECTION, 'spatialCoefficients', defaultvec5, (/0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk/))
        ! !search for parameters for the linearAlpha model (spatial)
        ! call read_param_mpi(PARAMS, SECTION, 'alphaTip', alpha_tip, -1e5_rk)
        ! call read_param_mpi(PARAMS, SECTION, 'alphaRoot', alpha_root, -1e5_rk)
        ! Insect%Wings(wingID)%deformation%D_alpha_tilde%coeffs = defaultvec5
        ! call read_param_mpi(PARAMS, SECTION, 'spatialFunction', Insect%Wings(wingID)%deformation%D_alpha_tilde%function_name, "polynomial")
        ! call read_param_mpi(PARAMS, SECTION, 'negativeTwistFraction', nt_frac, -1e5_rk)


        ! call read_param_mpi(PARAMS, SECTION, 'flexionLineAngle', Insect%Wings(wingID)%deformation%flexion_angle, 0.0_rk)
        ! call read_param_mpi(PARAMS, SECTION, 'flexionLinePivot', defaultvec, (/0.0_rk, 0.0_rk, 0.0_rk/))
        ! Insect%Wings(wingID)%deformation%pivot_flexion = defaultvec

        ! if ((alpha_tip/=-1e5_rk) .or. (alpha_root/=-1e5_rk) .or. (NReversals/=-1e5_rk) .or. (reversalPhase/=-1e5_rk)) then

        !     Insect%Wings(wingID)%deformation%D_alpha_tilde%coeffs(1) = alpha_root
        !     Insect%Wings(wingID)%deformation%D_alpha_tilde%coeffs(2) = alpha_tip
        !     Insect%Wings(wingID)%deformation%A_tilde%coeffs(2)  = NReversals
        !     Insect%Wings(wingID)%deformation%A_tilde%coeffs(3)  = reversalPhase
        !     Insect%Wings(wingID)%deformation%A_tilde%coeffs(4)  = nt_frac

        ! else
        !     call abort(0000000,"parameter not initialised!")
        ! end if

        ! call break_down_into_simple_model(Insect%Wings(wingID)%deformation)



    end do

    call break_down_into_simple_model(Insect%Wings(wingID)%deformation)

end subroutine