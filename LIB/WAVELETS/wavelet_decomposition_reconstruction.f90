    !-------------------------------------------------------------------------------
    ! Computes a one-level wavelet decomposition of a block.
    ! Data are stored in Spaghetti-order (not Mallat-Order)
    ! u  : input data
    ! u_d: output data (decomposed), can be u to perform inplace
    !
    ! Let's put an example with CDF22, BS=6, g=2:
    !
    !   n n n n n n n n n n           n n h g h g h g n n           n n h g h g h g n n
    !   n n n n n n n n n n           n n h g h g h g n n           n n h g h g h g n n
    !   n n i i i i i i n n           n n h g h g h g n n           n nghggghggghgg n n
    !   n n i i i i i i n n           n n h g h g h g n n           n nhhhghhhghhhg n n
    !   n n i i i i i i n n    ->     n n h g h g h g n n    ->     n nghggghggghgg n n
    !   n n i i i i i i n n           n n h g h g h g n n           n nhhhghhhghhhg n n
    !   n n i i i i i i n n           n n h g h g h g n n           n nghggghggghgg n n
    !   n n i i i i i i n n           n n h g h g h g n n           n nhhhghhhghhhg n n
    !   n n n n n n n n n n           n n h g h g h g n n           n n h g h g h g n n
    !   n n n n n n n n n n           n n h g h g h g n n           n n h g h g h g n n
    !
    ! n = ghost point, i = interior point, h = scaling filter applied, g = wavelet filter applied
    !
    !-------------------------------------------------------------------------------
    subroutine waveletDecomposition_optimized_block(params, u, u_d, SC_decompose, WC_decompose)
        implicit none
        type (type_params), intent(in) :: params
        real(kind=rk), dimension(:,:,:,:), intent(inout) :: u
        real(kind=rk), dimension(:,:,:,:), intent(out) :: u_d
        logical, optional, intent(in) :: SC_decompose, WC_decompose  ! Not yet implemented

        real(kind=rk), allocatable, dimension(:,:,:,:), save :: sc, wc, test, ucopy
        integer(kind=ik) :: nx, ny, nz, nc, g(1:3), Bs(1:3), ix, iy, iz, ic, Hs(1:2), Gs(1:2), maxn, f_size_max(1:3)
        logical :: SC_dec, WC_dec
        real(kind=rk), dimension(:), allocatable, save :: buffer

        ! By default, both decompositions are enabled
        SC_dec = .true.
        WC_dec = .true.
        if (present(SC_decompose)) SC_dec = SC_decompose
        if (present(WC_decompose)) WC_dec = WC_decompose

        nx = size(u,1)
        ny = size(u,2)
        nz = size(u,3)
        nc = size(u,4)
        g(:) = params%g  ! number of ghost points
        if (params%dim == 2) g(3) = 0
        Bs = params%Bs
        Hs = (/ lbound(params%HD,dim=1), ubound(params%HD,dim=1) /)  ! filter sizes
        Gs = (/ lbound(params%GD,dim=1), ubound(params%GD,dim=1) /)  ! filter sizes
        f_size_max(:) = max(maxval(Hs), maxval(Gs))
        if (params%dim == 2) f_size_max(3) = 0

        maxn = maxval((/ nx, ny, nz /))
        if (allocated(buffer)) then
            if (size(buffer, dim=1)<maxn) deallocate(buffer)
        endif
        if (.not.allocated(buffer)) allocate(buffer(1:maxn))

        ! we loop over all individual components one after another
        do ic = 1, nc

            ! we have different algorithms depending on the wavelet, for performance
            ! We then filter along each dimension consecutively
            ! First along X, then Y, then Z

            ! Some optimizations can generally be applied:
            ! we know maximum filter size and can skip some entries of future dimensions to treat
            ! ignore ghost points for dimensions that were already treated

            select case(params%wavelet)
            case ("CDF20")
                ! ~~~~~~~~~~~~~~~~~~~~~~ X ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do iy = g(2)+1-f_size_max(2), Bs(2)+g(2)+f_size_max(2)
                    ! low-pass filter (scaling function) - simple downpwass filter, copy values
                    u_d(g(1)+1:Bs(1)+g(1):2,iy,iz,ic) = u(g(1)+1:Bs(1)+g(1):2,iy,iz,ic)
                    ! high-pass filter (these guys are the details) - indices 0,1
                    do ix = g(1)+2, Bs(1)+g(1), 2
                        u_d(ix,iy,iz,ic) = u(ix-1,iy,iz,ic) * params%GD(-1) + u(ix,iy,iz,ic) * params%GD(0) + u(ix+1,iy,iz,ic) * params%GD( 1)
                    enddo
                enddo; enddo
                ! ~~~~~~~~~~~~~~~~~~~~~~ Y ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do ix = g(1)+1, Bs(1)+g(1)
                    ! low-pass filter (scaling function) - simple downpwass filter, inplace so we do nothing
                    ! high-pass filter (these guys are the details) - indices 0,1
                    do iy = g(2)+2, Bs(2)+g(2), 2
                        u_d(ix,iy,iz,ic) = u_d(ix,iy-1,iz,ic) * params%GD(-1) + u_d(ix,iy,iz,ic) * params%GD(0) + u_d(ix,iy+1,iz,ic) * params%GD( 1)
                    enddo
                enddo; enddo
                ! ~~~~~~~~~~~~~~~~~~~~~~ Z ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if (params%dim == 3) then
                    do iy = g(2)+1, Bs(2)+g(2); do ix = g(1)+1, Bs(1)+g(1)
                        ! low-pass filter (scaling function) - simple downpwass filter, inplace so we do nothing
                        ! high-pass filter (these guys are the details) - indices 0,1
                        do iz = g(3)+2, Bs(3)+g(3), 2
                            u_d(ix,iy,iz,ic) = u_d(ix,iy,iz-1,ic) * params%GD(-1) + u_d(ix,iy,iz,ic) * params%GD(0) + u_d(ix,iy,iz+1,ic) * params%GD( 1)
                        enddo
                    enddo; enddo
                endif
            case ("CDF40")
                ! ~~~~~~~~~~~~~~~~~~~~~~ X ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do iy = g(2)+1-f_size_max(2), Bs(2)+g(2)+f_size_max(2)
                    ! low-pass filter (scaling function) - simple downpwass filter, copy values
                    u_d(g(1)+1:Bs(1)+g(1):2,iy,iz,ic) = u(g(1)+1:Bs(1)+g(1):2,iy,iz,ic)
                    ! high-pass filter (these guys are the details) - indices 0,1,3
                    do ix = g(1)+2, Bs(1)+g(1), 2
                        u_d(ix,iy,iz,ic) = u(ix-3,iy,iz,ic) * params%GD(-3) + u(ix-1,iy,iz,ic) * params%GD(-1) + u(ix,iy,iz,ic) * params%GD(0) &
                                         + u(ix+1,iy,iz,ic) * params%GD( 1) + u(ix+3,iy,iz,ic) * params%GD( 3)
                    enddo
                enddo; enddo

                ! ~~~~~~~~~~~~~~~~~~~~~~ Y ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do ix = g(1)+1, Bs(1)+g(1)
                    ! low-pass filter (scaling function) - simple downpwass filter, inplace so we do nothing
                    ! high-pass filter (these guys are the details) - indices 0,1,3
                    do iy = g(2)+2, Bs(2)+g(2), 2
                        u_d(ix,iy,iz,ic) = u_d(ix,iy-3,iz,ic) * params%GD(-3) + u_d(ix,iy-1,iz,ic) * params%GD(-1) + u_d(ix,iy,iz,ic) * params%GD(0) &
                                         + u_d(ix,iy+1,iz,ic) * params%GD( 1) + u_d(ix,iy+3,iz,ic) * params%GD( 3)
                    enddo
                enddo; enddo

                ! ~~~~~~~~~~~~~~~~~~~~~~ Z ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if (params%dim == 3) then
                    do iy = g(2)+1, Bs(2)+g(2); do ix = g(1)+1, Bs(1)+g(1)
                        ! low-pass filter (scaling function) - simple downpwass filter, inplace so we do nothing
                        ! high-pass filter (these guys are the details) - indices 0,1,3
                        do iz = g(3)+2, Bs(3)+g(3), 2
                            u_d(ix,iy,iz,ic) = u_d(ix,iy,iz-3,ic) * params%GD(-3) + u_d(ix,iy,iz-1,ic) * params%GD(-1) + u_d(ix,iy,iz,ic) * params%GD(0) &
                                             + u_d(ix,iy,iz+1,ic) * params%GD( 1) + u_d(ix,iy,iz+3,ic) * params%GD( 3)
                        enddo
                    enddo; enddo
                endif
            case ("CDF60")
                ! ~~~~~~~~~~~~~~~~~~~~~~ X ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do iy = g(2)+1-f_size_max(2), Bs(2)+g(2)+f_size_max(2)
                    ! low-pass filter (scaling function) - simple downpwass filter, copy values
                    u_d(g(1)+1:Bs(1)+g(1):2,iy,iz,ic) = u(g(1)+1:Bs(1)+g(1):2,iy,iz,ic)
                    ! high-pass filter (these guys are the details) - indices 0,1,3,5
                    do ix = g(1)+2, Bs(1)+g(1), 2
                        u_d(ix,iy,iz,ic) = u(ix-5,iy,iz,ic) * params%GD(-5) + u(ix-3,iy,iz,ic) * params%GD(-3) + u(ix-1,iy,iz,ic) * params%GD(-1) + u(ix,iy,iz,ic) * params%GD(0) &
                                         + u(ix+1,iy,iz,ic) * params%GD( 1) + u(ix+3,iy,iz,ic) * params%GD( 3) + u(ix+5,iy,iz,ic) * params%GD( 5)
                    enddo
                enddo; enddo

                ! ~~~~~~~~~~~~~~~~~~~~~~ Y ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do ix = g(1)+1, Bs(1)+g(1)
                    ! low-pass filter (scaling function) - simple downpwass filter, inplace so we do nothing
                    ! high-pass filter (these guys are the details) - indices 0,1,3,5
                    do iy = g(2)+2, Bs(2)+g(2), 2
                        u_d(ix,iy,iz,ic) = u_d(ix,iy-5,iz,ic) * params%GD(-5) + u_d(ix,iy-3,iz,ic) * params%GD(-3) + u_d(ix,iy-1,iz,ic) * params%GD(-1) + u_d(ix,iy,iz,ic) * params%GD(0) &
                                         + u_d(ix,iy+1,iz,ic) * params%GD( 1) + u_d(ix,iy+3,iz,ic) * params%GD( 3) + u_d(ix,iy+5,iz,ic) * params%GD( 5)
                    enddo
                enddo; enddo

                ! ~~~~~~~~~~~~~~~~~~~~~~ Z ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if (params%dim == 3) then
                    do iy = g(2)+1, Bs(2)+g(2); do ix = g(1)+1, Bs(1)+g(1)
                        ! low-pass filter (scaling function) - simple downpwass filter, inplace so we do nothing
                        ! high-pass filter (these guys are the details) - indices 0,1,3,5
                        do iz = g(3)+2, Bs(3)+g(3), 2
                            u_d(ix,iy,iz,ic) = u_d(ix,iy,iz-5,ic) * params%GD(-5) + u_d(ix,iy,iz-3,ic) * params%GD(-3) + u_d(ix,iy,iz-1,ic) * params%GD(-1) + u_d(ix,iy,iz,ic) * params%GD(0) &
                                             + u_d(ix,iy,iz+1,ic) * params%GD( 1) + u_d(ix,iy,iz+3,ic) * params%GD( 3) + u_d(ix,iy,iz+5,ic) * params%GD( 5)
                        enddo
                    enddo; enddo
                endif
            case ("CDF80")
                ! ~~~~~~~~~~~~~~~~~~~~~~ X ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do iy = g(2)+1-f_size_max(2), Bs(2)+g(2)+f_size_max(2)
                    ! low-pass filter (scaling function) - simple downpwass filter, copy values
                    u_d(g(1)+1:Bs(1)+g(1):2,iy,iz,ic) = u(g(1)+1:Bs(1)+g(1):2,iy,iz,ic)
                    ! high-pass filter (these guys are the details) - indices 0,1,3,5,7
                    do ix = g(1)+2, Bs(1)+g(1), 2
                        u_d(ix,iy,iz,ic) = u(ix-7,iy,iz,ic) * params%GD(-7) + u(ix-5,iy,iz,ic) * params%GD(-5) + u(ix-3,iy,iz,ic) * params%GD(-3) + u(ix-1,iy,iz,ic) * params%GD(-1) + u(ix,iy,iz,ic) * params%GD(0) &
                                         + u(ix+1,iy,iz,ic) * params%GD( 1) + u(ix+3,iy,iz,ic) * params%GD( 3) + u(ix+5,iy,iz,ic) * params%GD( 5) + u(ix+7,iy,iz,ic) * params%GD( 7)
                    enddo
                enddo; enddo

                ! ~~~~~~~~~~~~~~~~~~~~~~ Y ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do ix = g(1)+1, Bs(1)+g(1)
                    ! low-pass filter (scaling function) - simple downpwass filter, inplace so we do nothing
                    ! high-pass filter (these guys are the details) - indices 0,1,3,5,7
                    do iy = g(2)+2, Bs(2)+g(2), 2
                        u_d(ix,iy,iz,ic) = u_d(ix,iy-7,iz,ic) * params%GD(-7) + u_d(ix,iy-5,iz,ic) * params%GD(-5) + u_d(ix,iy-3,iz,ic) * params%GD(-3) + u_d(ix,iy-1,iz,ic) * params%GD(-1) + u_d(ix,iy,iz,ic) * params%GD(0) &
                                         + u_d(ix,iy+1,iz,ic) * params%GD( 1) + u_d(ix,iy+3,iz,ic) * params%GD( 3) + u_d(ix,iy+5,iz,ic) * params%GD( 5) + u_d(ix,iy+7,iz,ic) * params%GD( 7)
                    enddo
                enddo; enddo

                ! ~~~~~~~~~~~~~~~~~~~~~~ Z ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if (params%dim == 3) then
                    do iy = g(2)+1, Bs(2)+g(2); do ix = g(1)+1, Bs(1)+g(1)
                        ! low-pass filter (scaling function) - simple downpwass filter, inplace so we do nothing
                        ! high-pass filter (these guys are the details) - indices 0,1,3,5,7
                        do iz = g(3)+2, Bs(3)+g(3), 2
                            u_d(ix,iy,iz,ic) = u_d(ix,iy,iz-7,ic) * params%GD(-7) + u_d(ix,iy,iz-5,ic) * params%GD(-5) + u_d(ix,iy,iz-3,ic) * params%GD(-3) + u_d(ix,iy,iz-1,ic) * params%GD(-1) + u_d(ix,iy,iz,ic) * params%GD(0) &
                                             + u_d(ix,iy,iz+1,ic) * params%GD( 1) + u_d(ix,iy,iz+3,ic) * params%GD( 3) + u_d(ix,iy,iz+5,ic) * params%GD( 5) + u_d(ix,iy,iz+7,ic) * params%GD( 7)
                        enddo
                    enddo; enddo
                endif
            case ("CDF22")
                ! ~~~~~~~~~~~~~~~~~~~~~~ X ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do iy = g(2)+1-f_size_max(2), Bs(2)+g(2)+f_size_max(2)
                    do ix = g(1)+1, Bs(1)+g(1), 2
                        ! low-pass filter (scaling function) - indices 0,1,2
                        buffer(ix) = u(ix-2,iy,iz,ic) * params%HD(-2) + u(ix-1,iy,iz,ic) * params%HD(-1) + u(ix,iy,iz,ic) * params%HD(0) + u(ix+1,iy,iz,ic) * params%HD( 1) + u(ix+2,iy,iz,ic) * params%HD( 2)
                        ! high-pass filter (these guys are the details) - indices 0,1 and acting on WC
                        buffer(ix+1) = u(ix,iy,iz,ic) * params%GD(-1) + u(ix+1,iy,iz,ic) * params%GD(0) + u(ix+2,iy,iz,ic) * params%GD( 1)
                    enddo
                    u_d(g(1)+1:Bs(1)+g(1),iy,iz,ic) = buffer(g(1)+1:Bs(1)+g(1))  ! insert values from buffer to avoid overwriting
                enddo; enddo
                ! ~~~~~~~~~~~~~~~~~~~~~~ Y ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do ix = g(1)+1, Bs(1)+g(1)
                    do iy = g(2)+1, Bs(2)+g(2), 2
                        ! low-pass filter (scaling function) - indices 0,1,2
                        buffer(iy) = u_d(ix,iy-2,iz,ic) * params%HD(-2) + u_d(ix,iy-1,iz,ic) * params%HD(-1) + u_d(ix,iy,iz,ic) * params%HD(0) + u_d(ix,iy+1,iz,ic) * params%HD( 1) + u_d(ix,iy+2,iz,ic) * params%HD( 2)
                        ! high-pass filter (these guys are the details) - indices 0,1 and acting on WC
                        buffer(iy+1) = u_d(ix,iy,iz,ic) * params%GD(-1) + u_d(ix,iy+1,iz,ic) * params%GD(0) + u_d(ix,iy+2,iz,ic) * params%GD( 1)
                    enddo
                    u_d(ix,g(2)+1:Bs(2)+g(2),iz,ic) = buffer(g(2)+1:Bs(2)+g(2))  ! insert values from buffer to avoid overwriting
                enddo; enddo
                ! ~~~~~~~~~~~~~~~~~~~~~~ Z ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if (params%dim == 3) then
                    do iy = g(2)+1, Bs(2)+g(2); do ix = g(1)+1, Bs(1)+g(1)
                        do iz = g(3)+1, Bs(3)+g(3), 2
                            ! low-pass filter (scaling function) - indices 0,1,2
                            buffer(iz) = u_d(ix,iy,iz-2,ic) * params%HD(-2) + u_d(ix,iy,iz-1,ic) * params%HD(-1) + u_d(ix,iy,iz,ic) * params%HD(0) + u_d(ix,iy,iz+1,ic) * params%HD( 1) + u_d(ix,iy,iz+2,ic) * params%HD( 2)
                            ! high-pass filter (these guys are the details) - indices 0,1 and acting on WC
                            buffer(iz+1) = u_d(ix,iy,iz,ic) * params%GD(-1) + u_d(ix,iy,iz+1,ic) * params%GD(0) + u_d(ix,iy,iz+2,ic) * params%GD( 1)
                        enddo
                        u_d(ix,iy,g(3)+1:Bs(3)+g(3),ic) = buffer(g(3)+1:Bs(3)+g(3))  ! insert values from buffer to avoid overwriting
                    enddo; enddo
                endif
            case ("CDF42")
                ! ~~~~~~~~~~~~~~~~~~~~~~ X ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do iy = g(2)+1-f_size_max(2), Bs(2)+g(2)+f_size_max(2)
                    do ix = g(1)+1, Bs(1)+g(1), 2
                        ! low-pass filter (scaling function) - indices 0,1,2,4
                        buffer(ix) = u(ix-4,iy,iz,ic) * params%HD(-4) + u(ix-2,iy,iz,ic) * params%HD(-2) + u(ix-1,iy,iz,ic) * params%HD(-1) + u(ix,iy,iz,ic) * params%HD(0) &
                                    + u(ix+1,iy,iz,ic) * params%HD( 1) + u(ix+2,iy,iz,ic) * params%HD( 2) + u(ix+4,iy,iz,ic) * params%HD( 4)
                        ! high-pass filter (these guys are the details) - indices 0,1,3 and acting on WC
                        buffer(ix+1) = u(ix-2,iy,iz,ic) * params%GD(-3) + u(ix,iy,iz,ic) * params%GD(-1) + u(ix+1,iy,iz,ic) * params%GD(0) + u(ix+2,iy,iz,ic) * params%GD( 1) + u(ix+4,iy,iz,ic) * params%GD( 3)
                    enddo
                    u_d(g(1)+1:Bs(1)+g(1),iy,iz,ic) = buffer(g(1)+1:Bs(1)+g(1))  ! insert values from buffer to avoid overwriting
                enddo; enddo

                ! ~~~~~~~~~~~~~~~~~~~~~~ Y ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do ix = g(1)+1, Bs(1)+g(1)
                    do iy = g(2)+1, Bs(2)+g(2), 2
                        ! low-pass filter (scaling function) - indices 0,1,2,4
                        buffer(iy) = u_d(ix,iy-4,iz,ic) * params%HD(-4) + u_d(ix,iy-2,iz,ic) * params%HD(-2) + u_d(ix,iy-1,iz,ic) * params%HD(-1) + u_d(ix,iy,iz,ic) * params%HD(0) &
                                    + u_d(ix,iy+1,iz,ic) * params%HD( 1) + u_d(ix,iy+2,iz,ic) * params%HD( 2) + u_d(ix,iy+4,iz,ic) * params%HD( 4)
                        ! high-pass filter (these guys are the details) - indices 0,1,3 and acting on WC
                        buffer(iy+1) = u_d(ix,iy-2,iz,ic) * params%GD(-3) + u_d(ix,iy,iz,ic) * params%GD(-1) + u_d(ix,iy+1,iz,ic) * params%GD(0) + u_d(ix,iy+2,iz,ic) * params%GD( 1) + u_d(ix,iy+4,iz,ic) * params%GD( 3)
                    enddo
                    u_d(ix,g(2)+1:Bs(2)+g(2),iz,ic) = buffer(g(2)+1:Bs(2)+g(2))  ! insert values from buffer to avoid overwriting
                enddo; enddo

                ! ~~~~~~~~~~~~~~~~~~~~~~ Z ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if (params%dim == 3) then
                    do iy = g(2)+1, Bs(2)+g(2); do ix = g(1)+1, Bs(1)+g(1)
                        do iz = g(3)+1, Bs(3)+g(3), 2
                            ! low-pass filter (scaling function) - indices 0,1,2,4
                            buffer(iz) = u_d(ix,iy,iz-4,ic) * params%HD(-4) + u_d(ix,iy,iz-2,ic) * params%HD(-2) + u_d(ix,iy,iz-1,ic) * params%HD(-1) + u_d(ix,iy,iz,ic) * params%HD(0) &
                                        + u_d(ix,iy,iz+1,ic) * params%HD( 1) + u_d(ix,iy,iz+2,ic) * params%HD( 2) + u_d(ix,iy,iz+4,ic) * params%HD( 4)
                            ! high-pass filter (these guys are the details) - indices 0,1,3 and acting on WC
                            buffer(iz+1) = u_d(ix,iy,iz-2,ic) * params%GD(-3) + u_d(ix,iy,iz,ic) * params%GD(-1) + u_d(ix,iy,iz+1,ic) * params%GD(0) + u_d(ix,iy,iz+2,ic) * params%GD( 1) + u_d(ix,iy,iz+4,ic) * params%GD( 3)
                        enddo
                        u_d(ix,iy,g(3)+1:Bs(3)+g(3),ic) = buffer(g(3)+1:Bs(3)+g(3))  ! insert values from buffer to avoid overwriting
                    enddo; enddo
                endif

            case ("CDF44")
                ! ~~~~~~~~~~~~~~~~~~~~~~ X ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do iy = g(2)+1-f_size_max(2), Bs(2)+g(2)+f_size_max(2)
                    do ix = g(1)+1, Bs(1)+g(1), 2
                        ! low-pass filter (scaling function) - indices 0,1,2,3,4,6
                        buffer(ix) = u(ix-6,iy,iz,ic) * params%HD(-6) + u(ix-4,iy,iz,ic) * params%HD(-4) + u(ix-3,iy,iz,ic) * params%HD(-3) + u(ix-2,iy,iz,ic) * params%HD(-2) + u(ix-1,iy,iz,ic) * params%HD(-1) + u(ix,iy,iz,ic) * params%HD(0) &
                                    + u(ix+1,iy,iz,ic) * params%HD( 1) + u(ix+2,iy,iz,ic) * params%HD( 2) + u(ix+3,iy,iz,ic) * params%HD( 3) + u(ix+4,iy,iz,ic) * params%HD( 4) + u(ix+6,iy,iz,ic) * params%HD( 6)
                        ! high-pass filter (these guys are the details) - indices 0,1,3 and acting on WC
                        buffer(ix+1) = u(ix-2,iy,iz,ic) * params%GD(-3) + u(ix,iy,iz,ic) * params%GD(-1) + u(ix+1,iy,iz,ic) * params%GD(0) + u(ix+2,iy,iz,ic) * params%GD( 1) + u(ix+4,iy,iz,ic) * params%GD( 3)
                    enddo
                    u_d(g(1)+1:Bs(1)+g(1),iy,iz,ic) = buffer(g(1)+1:Bs(1)+g(1))  ! insert values from buffer to avoid overwriting
                enddo; enddo
                ! ~~~~~~~~~~~~~~~~~~~~~~ Y ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do ix = g(1)+1, Bs(1)+g(1)
                    do iy = g(2)+1, Bs(2)+g(2), 2
                        ! low-pass filter (scaling function) - indices 0,1,2,3,4,6
                        buffer(iy) = u_d(ix,iy-6,iz,ic) * params%HD(-6) + u_d(ix,iy-4,iz,ic) * params%HD(-4) + u_d(ix,iy-3,iz,ic) * params%HD(-3) + u_d(ix,iy-2,iz,ic) * params%HD(-2) + u_d(ix,iy-1,iz,ic) * params%HD(-1) + u_d(ix,iy,iz,ic) * params%HD(0) &
                                    + u_d(ix,iy+1,iz,ic) * params%HD( 1) + u_d(ix,iy+2,iz,ic) * params%HD( 2) + u_d(ix,iy+3,iz,ic) * params%HD( 3) + u_d(ix,iy+4,iz,ic) * params%HD( 4) + u_d(ix,iy+6,iz,ic) * params%HD( 6)
                        ! high-pass filter (these guys are the details) - indices 0,1,3 and acting on WC
                        buffer(iy+1) = u_d(ix,iy-2,iz,ic) * params%GD(-3) + u_d(ix,iy,iz,ic) * params%GD(-1) + u_d(ix,iy+1,iz,ic) * params%GD(0) + u_d(ix,iy+2,iz,ic) * params%GD( 1) + u_d(ix,iy+4,iz,ic) * params%GD( 3)
                    enddo
                    u_d(ix,g(2)+1:Bs(2)+g(2),iz,ic) = buffer(g(2)+1:Bs(2)+g(2))  ! insert values from buffer to avoid overwriting
                enddo; enddo
                ! ~~~~~~~~~~~~~~~~~~~~~~ Z ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if (params%dim == 3) then
                    do iy = g(2)+1, Bs(2)+g(2); do ix = g(1)+1, Bs(1)+g(1)
                        do iz = g(3)+1, Bs(3)+g(3), 2
                            ! low-pass filter (scaling function) - indices 0,1,2,3,4,6
                            buffer(iz) = u_d(ix,iy,iz-6,ic) * params%HD(-6) + u_d(ix,iy,iz-4,ic) * params%HD(-4) + u_d(ix,iy,iz-3,ic) * params%HD(-3) + u_d(ix,iy,iz-2,ic) * params%HD(-2) + u_d(ix,iy,iz-1,ic) * params%HD(-1) + u_d(ix,iy,iz,ic) * params%HD(0) &
                                        + u_d(ix,iy,iz+1,ic) * params%HD( 1) + u_d(ix,iy,iz+2,ic) * params%HD( 2) + u_d(ix,iy,iz+3,ic) * params%HD( 3) + u_d(ix,iy,iz+4 ,ic) * params%HD( 4) + u_d(ix,iy,iz+6,ic) * params%HD( 6)
                            ! high-pass filter (these guys are the details) - indices 0,1,3 and acting on WC
                            buffer(iz+1) = u_d(ix,iy,iz-2,ic) * params%GD(-3) + u_d(ix,iy,iz,ic) * params%GD(-1) + u_d(ix,iy,iz+1,ic) * params%GD(0) + u_d(ix,iy,iz+2,ic) * params%GD( 1) + u_d(ix,iy,iz+4,ic) * params%GD( 3)
                        enddo
                        u_d(ix,iy,g(3)+1:Bs(3)+g(3),ic) = buffer(g(3)+1:Bs(3)+g(3))  ! insert values from buffer to avoid overwriting
                    enddo; enddo
                endif
            
            case ("CDF62")
                ! ~~~~~~~~~~~~~~~~~~~~~~ X ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do iy = g(2)+1-f_size_max(2), Bs(2)+g(2)+f_size_max(2)
                    do ix = g(1)+1, Bs(1)+g(1), 2
                        ! low-pass filter (scaling function) - indices 0,1,2,4,6
                        buffer(ix) = u(ix-6,iy,iz,ic) * params%HD(-6) + u(ix-4,iy,iz,ic) * params%HD(-4) + u(ix-2,iy,iz,ic) * params%HD(-2) + u(ix-1,iy,iz,ic) * params%HD(-1) + u(ix,iy,iz,ic) * params%HD(0) &
                                    + u(ix+1,iy,iz,ic) * params%HD( 1) + u(ix+2,iy,iz,ic) * params%HD( 2) + u(ix+4,iy,iz,ic) * params%HD( 4) + u(ix+6,iy,iz,ic) * params%HD( 6)
                        ! high-pass filter (these guys are the details) - indices 0,1,3,5 and acting on WC
                        buffer(ix+1) = u(ix-4,iy,iz,ic) * params%GD(-5) + u(ix-2,iy,iz,ic) * params%GD(-3) + u(ix,iy,iz,ic) * params%GD(-1) + u(ix+1,iy,iz,ic) * params%GD(0) + u(ix+2,iy,iz,ic) * params%GD( 1) + u(ix+4,iy,iz,ic) * params%GD( 3) + u(ix+6,iy,iz,ic) * params%GD( 5)
                    enddo
                    u_d(g(1)+1:Bs(1)+g(1),iy,iz,ic) = buffer(g(1)+1:Bs(1)+g(1))  ! insert values from buffer to avoid overwriting
                enddo; enddo
                ! ~~~~~~~~~~~~~~~~~~~~~~ Y ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do ix = g(1)+1, Bs(1)+g(1)
                    do iy = g(2)+1, Bs(2)+g(2), 2
                        ! low-pass filter (scaling function) - indices 0,1,2,4,6
                        buffer(iy) = u_d(ix,iy-6,iz,ic) * params%HD(-6) + u_d(ix,iy-4,iz,ic) * params%HD(-4) + u_d(ix,iy-2,iz,ic) * params%HD(-2) + u_d(ix,iy-1,iz,ic) * params%HD(-1) + u_d(ix,iy,iz,ic) * params%HD(0) &
                                    + u_d(ix,iy+1,iz,ic) * params%HD( 1) + u_d(ix,iy+2,iz,ic) * params%HD( 2) + u_d(ix,iy+4,iz,ic) * params%HD( 4) + u_d(ix,iy+6,iz,ic) * params%HD( 6)
                        ! high-pass filter (these guys are the details) - indices 0,1,3,5 and acting on WC
                        buffer(iy+1) = u_d(ix,iy-4,iz,ic) * params%GD(-5) + u_d(ix,iy-2,iz,ic) * params%GD(-3) + u_d(ix,iy,iz,ic) * params%GD(-1) + u_d(ix,iy+1,iz,ic) * params%GD(0) + u_d(ix,iy+2,iz,ic) * params%GD( 1) + u_d(ix,iy+4,iz,ic) * params%GD( 3) + u_d(ix,iy+6,iz,ic) * params%GD( 5)
                    enddo
                    u_d(ix,g(2)+1:Bs(2)+g(2),iz,ic) = buffer(g(2)+1:Bs(2)+g(2))  ! insert values from buffer to avoid overwriting
                enddo; enddo
                ! ~~~~~~~~~~~~~~~~~~~~~~ Z ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if (params%dim == 3) then
                    do iy = g(2)+1, Bs(2)+g(2); do ix = g(1)+1, Bs(1)+g(1)
                        do iz = g(3)+1, Bs(3)+g(3), 2
                            ! low-pass filter (scaling function) - indices 0,1,2,4,6
                            buffer(iz) = u_d(ix,iy,iz-6,ic) * params%HD(-6) + u_d(ix,iy,iz-4,ic) * params%HD(-4) + u_d(ix,iy,iz-2,ic) * params%HD(-2) + u_d(ix,iy,iz-1,ic) * params%HD(-1) + u_d(ix,iy,iz,ic) * params%HD(0) &
                                        + u_d(ix,iy,iz+1,ic) * params%HD( 1) + u_d(ix,iy,iz+2,ic) * params%HD( 2) + u_d(ix,iy,iz+4,ic) * params%HD( 4) + u_d(ix,iy,iz+6,ic) * params%HD( 6)
                            ! high-pass filter (these guys are the details) - indices 0,1,3,5 and acting on WC
                            buffer(iz+1) = u_d(ix,iy,iz-4,ic) * params%GD(-5) + u_d(ix,iy,iz-2,ic) * params%GD(-3) + u_d(ix,iy,iz,ic) * params%GD(-1) + u_d(ix,iy,iz+1,ic) * params%GD(0) + u_d(ix,iy,iz+2,ic) * params%GD( 1) + u_d(ix,iy,iz+4,ic) * params%GD( 3) + u_d(ix,iy,iz+6,ic) * params%GD( 5)
                        enddo
                        u_d(ix,iy,g(3)+1:Bs(3)+g(3),ic) = buffer(g(3)+1:Bs(3)+g(3))  ! insert values from buffer to avoid overwriting
                    enddo; enddo
                endif

            ! Lazy cases - some wavelets are not often used and have not been hardcoded for performance, hence we use a general algorithm
            case default
                ! ~~~~~~~~~~~~~~~~~~~~~~ X ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do iy = g(2)+1-f_size_max(2), Bs(2)+g(2)+f_size_max(2)
                    ! low-pass filter (scaling function)
                    do ix = g(1)+1, Bs(1)+g(1), 2
                        buffer(ix) = sum(u(ix+Hs(1):ix+Hs(2),iy,iz,ic) * params%HD(:))
                    enddo
                    ! high-pass filter (these guys are the details)
                    do ix = g(1)+2, Bs(1)+g(1), 2
                        buffer(ix) = sum(u(ix+Gs(1):ix+Gs(2),iy,iz,ic) * params%GD(:))
                    enddo
                    u_d(g(1)+1:Bs(1)+g(1),iy,iz,ic) = buffer(g(1)+1:Bs(1)+g(1))  ! insert values from buffer to avoid overwriting
                enddo; enddo

                ! ~~~~~~~~~~~~~~~~~~~~~~ Y ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do ix = g(1)+1, Bs(1)+g(1)
                    ! low-pass filter (scaling function)
                    do iy = g(2)+1, Bs(2)+g(2), 2
                        buffer(iy) = sum(u_d(ix,iy+Hs(1):iy+Hs(2),iz,ic) * params%HD(:))
                    enddo
                    ! high-pass filter (these guys are the details)
                    do iy = g(2)+2, Bs(2)+g(2), 2
                        buffer(iy) = sum(u_d(ix,iy+Gs(1):iy+Gs(2),iz,ic) * params%GD(:))
                    enddo
                    u_d(ix,g(2)+1:Bs(2)+g(2),iz,ic) = buffer(g(2)+1:Bs(2)+g(2))  ! insert values from buffer to avoid overwriting
                enddo; enddo

                ! ~~~~~~~~~~~~~~~~~~~~~~ Z ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if (params%dim == 3) then
                    do iy = g(2)+1, Bs(2)+g(2); do ix = g(1)+1, Bs(1)+g(1)
                        ! low-pass filter (scaling function)
                        do iz = g(3)+1, Bs(3)+g(3), 2
                            buffer(iz) = sum(u_d(ix,iy,iz+Hs(1):iz+Hs(2),ic) * params%HD(:))
                        enddo
                        ! high-pass filter (these guys are the details)
                        do iz = g(3)+2, Bs(3)+g(3), 2
                            buffer(iz) = sum(u_d(ix,iy,iz+Gs(1):iz+Gs(2),ic) * params%GD(:))
                        enddo
                        u_d(ix,iy,g(3)+1:Bs(3)+g(3),ic) = buffer(g(3)+1:Bs(3)+g(3))  ! insert values from buffer to avoid overwriting
                    enddo; enddo
                endif
            end select

        enddo
    end subroutine


    !-----------------------------------------------------------------------------
    ! Reconstruction from low- and high pass filtered coefficients.
    ! Data is first upsampled, then filtered with the reconstruction filters.
    ! Note reconstruction filters are reverse of decomposition filters.
    !
    ! Input: hh gh hh gh hh gh hh gh
    !        hg gg hg gg hg gg hg gg
    !        hh gh hh gh hh gh hh gh
    !        hg gg hg gg hg gg hg gg
    !        hh gh hh gh hh gh hh gh
    !        hg gg hg gg hg gg hg gg
    ! Note: in input in spaghetti ordering is synced
    ! Output: original data, can be u in order to perform inplace
    !
    ! Let's put an example with CDF22, BS=6, g=2:
    !
    !   n n n n n n n n n n           n 0 n 0 n 0 n 0 n 0          0 n 0 n 0 n 0 n 0 n
    !   n n n n n n n n n n           n 0 n 0 n 0 n 0 n 0          0 n 0 n 0 n 0 n 0 n
    !   n nghggghggghgg n n           n 0gh 0gh 0gh 0 n 0          0 n 0gg 0gg 0gg 0 n
    !   n nhhhghhhghhhg n n           n 0hh 0hh 0hh 0 n 0          0 n 0hg 0hg 0hg 0 n
    !   n nghggghggghgg n n    ->     n 0gh 0gh 0gh 0 n 0    +     0 n 0gg 0gg 0gg 0 n
    !   n nhhhghhhghhhg n n           n 0hh 0hh 0hh 0 n 0          0 n 0hg 0hg 0hg 0 n
    !   n nghggghggghgg n n           n 0gh 0gh 0gh 0 n 0          0 n 0gg 0gg 0gg 0 n
    !   n nhhhghhhghhhg n n           n 0hh 0hh 0hh 0 n 0          0 n 0hg 0hg 0hg 0 n
    !   n n n n n n n n n n           n 0 n 0 n 0 n 0 n 0          0 n 0 n 0 n 0 n 0 n
    !   n n n n n n n n n n           n 0 n 0 n 0 n 0 n 0          0 n 0 n 0 n 0 n 0 n
    !
    ! For reconstruction, we separate the two contributions and apply on the first part HR and on the second GR
    ! This is repeated for every dimension
    ! n = ghost point, h = scaling filter applied, g = wavelet filter applied, 0 = zero (upsampled point)
    !
    !-------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    subroutine waveletReconstruction_optimized_block(params, u, u_r, SC_reconstruct, WC_reconstruct)
        implicit none
        type (type_params), intent(in) :: params
        real(kind=rk), dimension(:,:,:,:), intent(inout) :: u
        real(kind=rk), dimension(:,:,:,:), intent(out) :: u_r
        logical, optional, intent(in) :: SC_reconstruct, WC_reconstruct  ! Not yet implemented

        integer(kind=ik) :: nx, ny, nz, nc, g(1:3), Bs(1:3), ix, iy, iz, ic, Hs(1:2), Gs(1:2), maxn, io(1:3), f_size_max(1:3)
        logical :: SC_rec, WC_rec
        real(kind=rk), dimension(:), allocatable, save :: buffer, buffer_out

        ! By default, both reconstructions are enabled
        SC_rec = .true.
        WC_rec = .true.
        if (present(SC_reconstruct)) SC_rec = SC_reconstruct
        if (present(WC_reconstruct)) WC_rec = WC_reconstruct

        nx = size(u,1)
        ny = size(u,2)
        nz = size(u,3)
        nc = size(u,4)
        g(:) = params%g  ! number of ghost points
        if (params%dim == 2) g(3) = 0
        Bs = params%Bs
        Hs = (/ lbound(params%HR,dim=1), ubound(params%HR,dim=1) /)  ! filter sizes
        Gs = (/ lbound(params%GR,dim=1), ubound(params%GR,dim=1) /)  ! filter sizes
        f_size_max = max(maxval(Hs), maxval(Gs))
        if (params%dim == 2) f_size_max(3) = 0
        ! we need to shift if g is even as SCs start from second point
        io = 0
        do ic = 1, params%dim
            if (modulo(g(ic),2)==1) io(ic) = 1
        enddo

        maxn = maxval((/ nx, ny, nz /))
        if (allocated(buffer)) then
            if (size(buffer, dim=1)<maxn) deallocate(buffer)
        endif
        if (.not.allocated(buffer)) allocate(buffer(1:maxn))
        if (allocated(buffer_out)) then
            if (size(buffer_out, dim=1)<maxn) deallocate(buffer_out)
        endif
        if (.not.allocated(buffer_out)) allocate(buffer_out(1:maxn))

        ! we loop over all individual components one after another
        do ic = 1, nc

            ! we have different algorithms depending on the wavelet, for performance
            ! We then filter along each dimension consecutively
            ! First along X, then Y, then Z

            ! Some optimizations can generally be applied:
            ! we know maximum filter size and can skip some entries of future dimensions to treat
            ! ignore ghost points for dimensions that were already treated

            select case(params%wavelet)
            case ("CDF20")
                ! Unlifted case - we copy all values and only add interpolated values from odd indices
                u_r(:,:,:,ic) = u(:,:,:,ic)  ! copy values once, this is like applying center value of HR and full GR
                ! ~~~~~~~~~~~~~~~~~~~~~~ X ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do iy = g(2)+1-f_size_max(2), Bs(2)+g(2)+f_size_max(2)
                    do ix = g(1)+1+io(1), Bs(1)+g(1), 2
                        buffer_out(ix) = u_r(ix,iy,iz,ic) + u(ix-1,iy,iz,ic) * params%HR(-1) + u(ix+1,iy,iz,ic) * params%HR( 1)
                    enddo
                    u_r(g(1)+1+io(1):Bs(1)+g(1):2,iy,iz,ic) = buffer_out(g(1)+1+io(1):Bs(1)+g(1):2)  ! insert values from buffer to avoid overwriting
                enddo; enddo
                ! ~~~~~~~~~~~~~~~~~~~~~~ Y ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do ix = g(1)+1, Bs(1)+g(1)
                    do iy = g(2)+1+io(2), Bs(2)+g(2), 2
                        buffer_out(iy) = u_r(ix,iy,iz,ic) + u_r(ix,iy-1,iz,ic) * params%HR(-1) + u_r(ix,iy+1,iz,ic) * params%HR( 1)
                    enddo
                    u_r(ix,g(2)+1+io(2):Bs(2)+g(2):2,iz,ic) = buffer_out(g(2)+1+io(2):Bs(2)+g(2):2)  ! insert values from buffer to avoid overwriting
                enddo; enddo
                ! ~~~~~~~~~~~~~~~~~~~~~~ Z ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if (params%dim == 3) then
                    do iy = g(2)+1, Bs(2)+g(2); do ix = g(1)+1, Bs(1)+g(1)
                        do iz = g(3)+1+io(3), Bs(3)+g(3), 2
                            buffer_out(iz) = u_r(ix,iy,iz,ic) + u_r(ix,iy,iz-1,ic) * params%HR(-1) + u_r(ix,iy,iz+1,ic) * params%HR( 1)
                        enddo
                        u_r(ix,iy,g(3)+1+io(3):Bs(3)+g(3):2,ic) = buffer_out(g(3)+1+io(3):Bs(3)+g(3):2)  ! insert values from buffer to avoid overwriting
                    enddo; enddo
                endif
            case ("CDF40")
                ! Unlifted case - we copy all values and only add interpolated values from odd indices
                u_r(:,:,:,ic) = u(:,:,:,ic)  ! copy values once, this is like applying center value of HR and full GR
                ! ~~~~~~~~~~~~~~~~~~~~~~ X ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do iy = g(2)+1-f_size_max(2), Bs(2)+g(2)+f_size_max(2)
                    do ix = g(1)+1+io(1), Bs(1)+g(1), 2
                        buffer_out(ix) = u_r(ix,iy,iz,ic) + u(ix-3,iy,iz,ic) * params%HR(-3) + u(ix-1,iy,iz,ic) * params%HR(-1) &
                                                          + u(ix+1,iy,iz,ic) * params%HR( 1) + u(ix+3,iy,iz,ic) * params%HR( 3)
                    enddo
                    u_r(g(1)+1+io(1):Bs(1)+g(1):2,iy,iz,ic) = buffer_out(g(1)+1+io(1):Bs(1)+g(1):2)  ! insert values from buffer to avoid overwriting
                enddo; enddo
                ! ~~~~~~~~~~~~~~~~~~~~~~ Y ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do ix = g(1)+1, Bs(1)+g(1)
                    do iy = g(2)+1+io(2), Bs(2)+g(2), 2
                        buffer_out(iy) = u_r(ix,iy,iz,ic) + u_r(ix,iy-3,iz,ic) * params%HR(-3) + u_r(ix,iy-1,iz,ic) * params%HR(-1) &
                                                          + u_r(ix,iy+1,iz,ic) * params%HR( 1) + u_r(ix,iy+3,iz,ic) * params%HR( 3)
                    enddo
                    u_r(ix,g(2)+1+io(2):Bs(2)+g(2):2,iz,ic) = buffer_out(g(2)+1+io(2):Bs(2)+g(2):2)  ! insert values from buffer to avoid overwriting
                enddo; enddo
                ! ~~~~~~~~~~~~~~~~~~~~~~ Z ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if (params%dim == 3) then
                    do iy = g(2)+1, Bs(2)+g(2); do ix = g(1)+1, Bs(1)+g(1)
                        do iz = g(3)+1+io(3), Bs(3)+g(3), 2
                            buffer_out(iz) = u_r(ix,iy,iz,ic) + u_r(ix,iy,iz-3,ic) * params%HR(-3) + u_r(ix,iy,iz-1,ic) * params%HR(-1) &
                                                              + u_r(ix,iy,iz+1,ic) * params%HR( 1) + u_r(ix,iy,iz+3,ic) * params%HR( 3)
                        enddo
                        u_r(ix,iy,g(3)+1+io(3):Bs(3)+g(3):2,ic) = buffer_out(g(3)+1+io(3):Bs(3)+g(3):2)  ! insert values from buffer to avoid overwriting
                    enddo; enddo
                endif
            case ("CDF60")
                ! Unlifted case - we copy all values and only add interpolated values from odd indices
                u_r(:,:,:,ic) = u(:,:,:,ic)  ! copy values once, this is like applying center value of HR and full GR
                ! ~~~~~~~~~~~~~~~~~~~~~~ X ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do iy = g(2)+1-f_size_max(2), Bs(2)+g(2)+f_size_max(2)
                    do ix = g(1)+1+io(1), Bs(1)+g(1), 2
                        buffer_out(ix) = u_r(ix,iy,iz,ic) + u(ix-5,iy,iz,ic) * params%HR(-5) + u(ix-3,iy,iz,ic) * params%HR(-3) + u(ix-1,iy,iz,ic) * params%HR(-1) &
                                                          + u(ix+1,iy,iz,ic) * params%HR( 1) + u(ix+3,iy,iz,ic) * params%HR( 3) + u(ix+5,iy,iz,ic) * params%HR( 5)
                    enddo
                    u_r(g(1)+1+io(1):Bs(1)+g(1):2,iy,iz,ic) = buffer_out(g(1)+1+io(1):Bs(1)+g(1):2)  ! insert values from buffer to avoid overwriting
                enddo; enddo
                ! ~~~~~~~~~~~~~~~~~~~~~~ Y ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do ix = g(1)+1, Bs(1)+g(1)
                    do iy = g(2)+1+io(2), Bs(2)+g(2), 2
                        buffer_out(iy) = u_r(ix,iy,iz,ic) + u_r(ix,iy-5,iz,ic) * params%HR(-5) + u_r(ix,iy-3,iz,ic) * params%HR(-3) + u_r(ix,iy-1,iz,ic) * params%HR(-1) &
                                                          + u_r(ix,iy+1,iz,ic) * params%HR( 1) + u_r(ix,iy+3,iz,ic) * params%HR( 3) + u_r(ix,iy+5,iz,ic) * params%HR( 5)
                    enddo
                    u_r(ix,g(2)+1+io(2):Bs(2)+g(2):2,iz,ic) = buffer_out(g(2)+1+io(2):Bs(2)+g(2):2)  ! insert values from buffer to avoid overwriting
                enddo; enddo
                ! ~~~~~~~~~~~~~~~~~~~~~~ Z ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if (params%dim == 3) then
                    do iy = g(2)+1, Bs(2)+g(2); do ix = g(1)+1, Bs(1)+g(1)
                        do iz = g(3)+1+io(3), Bs(3)+g(3), 2
                            buffer_out(iz) = u_r(ix,iy,iz,ic) + u_r(ix,iy,iz-5,ic) * params%HR(-5) + u_r(ix,iy,iz-3,ic) * params%HR(-3) + u_r(ix,iy,iz-1,ic) * params%HR(-1) &
                                                              + u_r(ix,iy,iz+1,ic) * params%HR( 1) + u_r(ix,iy,iz+3,ic) * params%HR( 3) + u_r(ix,iy,iz+5,ic) * params%HR( 5)
                        enddo
                        u_r(ix,iy,g(3)+1+io(3):Bs(3)+g(3):2,ic) = buffer_out(g(3)+1+io(3):Bs(3)+g(3):2)  ! insert values from buffer to avoid overwriting
                    enddo; enddo
                endif
            case ("CDF80")
                ! Unlifted case - we copy all values and only add interpolated values from odd indices
                u_r(:,:,:,ic) = u(:,:,:,ic)  ! copy values once, this is like applying center value of HR and full GR
                ! ~~~~~~~~~~~~~~~~~~~~~~ X ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do iy = g(2)+1-f_size_max(2), Bs(2)+g(2)+f_size_max(2)
                    do ix = g(1)+1+io(1), Bs(1)+g(1), 2
                        buffer_out(ix) = u_r(ix,iy,iz,ic) + u(ix-7,iy,iz,ic) * params%HR(-7) + u(ix-5,iy,iz,ic) * params%HR(-5) + u(ix-3,iy,iz,ic) * params%HR(-3) + u(ix-1,iy,iz,ic) * params%HR(-1) &
                                                          + u(ix+1,iy,iz,ic) * params%HR( 1) + u(ix+3,iy,iz,ic) * params%HR( 3) + u(ix+5,iy,iz,ic) * params%HR( 5) + u(ix+7,iy,iz,ic) * params%HR( 7)
                    enddo
                    u_r(g(1)+1+io(1):Bs(1)+g(1):2,iy,iz,ic) = buffer_out(g(1)+1+io(1):Bs(1)+g(1):2)  ! insert values from buffer to avoid overwriting
                enddo; enddo
                ! ~~~~~~~~~~~~~~~~~~~~~~ Y ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do ix = g(1)+1, Bs(1)+g(1)
                    do iy = g(2)+1+io(2), Bs(2)+g(2), 2
                        buffer_out(iy) = u_r(ix,iy,iz,ic) + u_r(ix,iy-7,iz,ic) * params%HR(-7) + u_r(ix,iy-5,iz,ic) * params%HR(-5) + u_r(ix,iy-3,iz,ic) * params%HR(-3) + u_r(ix,iy-1,iz,ic) * params%HR(-1) &
                                                          + u_r(ix,iy+1,iz,ic) * params%HR( 1) + u_r(ix,iy+3,iz,ic) * params%HR( 3) + u_r(ix,iy+5,iz,ic) * params%HR( 5) + u_r(ix,iy+7,iz,ic) * params%HR( 7)
                    enddo
                    u_r(ix,g(2)+1+io(2):Bs(2)+g(2):2,iz,ic) = buffer_out(g(2)+1+io(2):Bs(2)+g(2):2)  ! insert values from buffer to avoid overwriting
                enddo; enddo
                ! ~~~~~~~~~~~~~~~~~~~~~~ Z ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if (params%dim == 3) then
                    do iy = g(2)+1, Bs(2)+g(2); do ix = g(1)+1, Bs(1)+g(1)
                        do iz = g(3)+1+io(3), Bs(3)+g(3), 2
                            buffer_out(iz) = u_r(ix,iy,iz,ic) + u_r(ix,iy,iz-7,ic) * params%HR(-7) + u_r(ix,iy,iz-5,ic) * params%HR(-5) + u_r(ix,iy,iz-3,ic) * params%HR(-3) + u_r(ix,iy,iz-1,ic) * params%HR(-1) &
                                                              + u_r(ix,iy,iz+1,ic) * params%HR( 1) + u_r(ix,iy,iz+3,ic) * params%HR( 3) + u_r(ix,iy,iz+5,ic) * params%HR( 5) + u_r(ix,iy,iz+7,ic) * params%HR( 7)
                        enddo
                        u_r(ix,iy,g(3)+1+io(3):Bs(3)+g(3):2,ic) = buffer_out(g(3)+1+io(3):Bs(3)+g(3):2)  ! insert values from buffer to avoid overwriting
                    enddo; enddo
                endif
            case ("CDF22")
                ! ~~~~~~~~~~~~~~~~~~~~~~ X ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do iy = g(2)+1-f_size_max(2), Bs(2)+g(2)+f_size_max(2)
                    do ix = g(1)+1, Bs(1)+g(1), 2
                        ! interpolation filter - copy SC indices, interpolate on WC indices
                        buffer_out(ix) = u( ix, iy, iz, ic)
                        buffer_out(ix+1) = u( ix, iy, iz, ic) * params%HR(-1) + u( ix+2, iy, iz, ic) * params%HR( 1)
                        ! lifted WR filter - indices 0,1,2, only half is active at any time
                        buffer_out(ix) = buffer_out(ix) + u( ix-1, iy, iz, ic) * params%GR(-1) + u( ix+1, iy, iz, ic) * params%GR( 1)
                        buffer_out(ix+1) = buffer_out(ix+1) + u( ix+1, iy, iz, ic) * params%GR(0) &
                            + u( ix-1, iy, iz, ic) * params%GR(-2) + u( ix+3, iy, iz, ic) * params%GR( 2)
                    enddo
                    u_r(g(1)+1:Bs(1)+g(1),iy,iz,ic) = buffer_out(g(1)+1:Bs(1)+g(1))  ! insert values from buffer to avoid overwriting
                enddo; enddo
                ! ~~~~~~~~~~~~~~~~~~~~~~ Y ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do ix = g(1)+1, Bs(1)+g(1)
                    do iy = g(2)+1, Bs(2)+g(2), 2
                        ! interpolation filter - copy SC indices, interpolate on WC indices
                        buffer_out(iy) = u_r( ix, iy, iz, ic)
                        buffer_out(iy+1) = u_r( ix, iy, iz, ic) * params%HR(-1) + u_r( ix, iy+2, iz, ic) * params%HR( 1)
                        ! lifted WR filter - indices 0,1,2, only half is active at any time
                        buffer_out(iy) = buffer_out(iy) + u_r( ix, iy-1, iz, ic) * params%GR(-1) + u_r( ix, iy+1, iz, ic) * params%GR( 1)
                        buffer_out(iy+1) = buffer_out(iy+1) + u_r( ix, iy+1, iz, ic) * params%GR(0) &
                            + u_r( ix, iy-1, iz, ic) * params%GR(-2) + u_r( ix, iy+3, iz, ic) * params%GR( 2)
                    enddo
                    u_r(ix,g(2)+1:Bs(2)+g(2),iz,ic) = buffer_out(g(2)+1:Bs(2)+g(2))  ! insert values from buffer to avoid overwriting
                enddo; enddo
                ! ~~~~~~~~~~~~~~~~~~~~~~ Z ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if (params%dim == 3) then
                    do iy = g(2)+1, Bs(2)+g(2); do ix = g(1)+1, Bs(1)+g(1)
                        do iz = g(3)+1, Bs(3)+g(3), 2
                            ! interpolation filter - copy SC indices, interpolate on WC indices
                            buffer_out(iz) = u_r( ix, iy, iz, ic)
                            buffer_out(iz+1) = u_r( ix, iy, iz, ic) * params%HR(-1) + u_r( ix, iy, iz+2, ic) * params%HR( 1)
                            ! lifted WR filter - indices 0,1,2, only half is active at any time
                            buffer_out(iz) = buffer_out(iz) + u_r( ix, iy, iz-1, ic) * params%GR(-1) + u_r( ix, iy, iz+1, ic) * params%GR( 1)
                            buffer_out(iz+1) = buffer_out(iz+1) + u_r( ix, iy, iz+1, ic) * params%GR(0) &
                                + u_r( ix, iy, iz-1, ic) * params%GR(-2) + u_r( ix, iy, iz+3, ic) * params%GR( 2)
                        enddo
                        u_r(ix,iy,g(3)+1:Bs(3)+g(3),ic) = buffer_out(g(3)+1:Bs(3)+g(3))  ! insert values from buffer to avoid overwriting
                    enddo; enddo
                endif

            case ("CDF42")
                ! ~~~~~~~~~~~~~~~~~~~~~~ X ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do iy = g(2)+1-f_size_max(2), Bs(2)+g(2)+f_size_max(2)
                    do ix = g(1)+1, Bs(1)+g(1), 2
                        ! interpolation filter - copy SC indices, interpolate on WC indices
                        buffer_out(ix) = u( ix, iy, iz, ic)
                        buffer_out(ix+1) = u( ix-2, iy, iz, ic) * params%HR(-3) + u( ix  , iy, iz, ic) * params%HR(-1) &
                                           + u( ix+2, iy, iz, ic) * params%HR( 1) + u( ix+4, iy, iz, ic) * params%HR( 3)
                        ! lifted WR filter - indices 0,1,2,4, only half is active at any time
                        buffer_out(ix) = buffer_out(ix) + u( ix-1, iy, iz, ic) * params%GR(-1) + u( ix+1, iy, iz, ic) * params%GR( 1)
                        buffer_out(ix+1) = buffer_out(ix+1) + u( ix+1, iy, iz, ic) * params%GR(0) &
                            + u( ix-3, iy, iz, ic) * params%GR(-4) + u( ix-1, iy, iz, ic) * params%GR(-2) + u( ix+3, iy, iz, ic) * params%GR( 2) + u( ix+5, iy, iz, ic) * params%GR( 4)
                    enddo
                    u_r(g(1)+1:Bs(1)+g(1),iy,iz,ic) = buffer_out(g(1)+1:Bs(1)+g(1))  ! insert values from buffer to avoid overwriting
                enddo; enddo
                ! ~~~~~~~~~~~~~~~~~~~~~~ Y ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do ix = g(1)+1, Bs(1)+g(1)
                    do iy = g(2)+1, Bs(2)+g(2), 2
                        ! interpolation filter - copy SC indices, interpolate on WC indices
                        buffer_out(iy) = u_r( ix, iy, iz, ic)
                        buffer_out(iy+1) = u_r( ix, iy-2, iz, ic) * params%HR(-3) + u_r( ix, iy  , iz, ic) * params%HR(-1) &
                                         + u_r( ix, iy+2, iz, ic) * params%HR( 1) + u_r( ix, iy+4, iz, ic) * params%HR( 3)
                        ! lifted WR filter - indices 0,1,2,4, only half is active at any time
                        buffer_out(iy) = buffer_out(iy) + u_r( ix, iy-1, iz, ic) * params%GR(-1) + u_r( ix, iy+1, iz, ic) * params%GR( 1)
                        buffer_out(iy+1) = buffer_out(iy+1) + u_r( ix, iy+1, iz, ic) * params%GR(0) &
                            + u_r( ix, iy-3, iz, ic) * params%GR(-4) + u_r( ix, iy-1, iz, ic) * params%GR(-2) + u_r( ix, iy+3, iz, ic) * params%GR( 2) + u_r( ix, iy+5, iz, ic) * params%GR( 4)
                    enddo
                    u_r(ix,g(2)+1:Bs(2)+g(2),iz,ic) = buffer_out(g(2)+1:Bs(2)+g(2))  ! insert values from buffer to avoid overwriting
                enddo; enddo
                ! ~~~~~~~~~~~~~~~~~~~~~~ Z ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if (params%dim == 3) then
                    do iy = g(2)+1, Bs(2)+g(2); do ix = g(1)+1, Bs(1)+g(1)
                        do iz = g(3)+1, Bs(3)+g(3), 2
                            ! interpolation filter - copy SC indices, interpolate on WC indices
                            buffer_out(iz) = u_r( ix, iy, iz, ic)
                            buffer_out(iz+1) = u_r( ix, iy, iz-2, ic) * params%HR(-3) + u_r( ix, iy, iz  , ic) * params%HR(-1) &
                                             + u_r( ix, iy, iz+2, ic) * params%HR( 1) + u_r( ix, iy, iz+4, ic) * params%HR( 3)
                            ! lifted WR filter - indices 0,1,2,4, only half is active at any time
                            buffer_out(iz) = buffer_out(iz) + u_r( ix, iy, iz-1, ic) * params%GR(-1) + u_r( ix, iy, iz+1, ic) * params%GR( 1)
                            buffer_out(iz+1) = buffer_out(iz+1) + u_r( ix, iy, iz+1, ic) * params%GR(0) &
                                + u_r( ix, iy, iz-3, ic) * params%GR(-4) + u_r( ix, iy, iz-1, ic) * params%GR(-2) + u_r( ix, iy, iz+3, ic) * params%GR( 2) + u_r( ix, iy, iz+5, ic) * params%GR( 4)
                        enddo
                        u_r(ix,iy,g(3)+1:Bs(3)+g(3),ic) = buffer_out(g(3)+1:Bs(3)+g(3))  ! insert values from buffer to avoid overwriting
                    enddo; enddo
                endif
            
            case ("CDF44")
                ! ~~~~~~~~~~~~~~~~~~~~~~ X ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do iy = g(2)+1-f_size_max(2), Bs(2)+g(2)+f_size_max(2)
                    do ix = g(1)+1, Bs(1)+g(1), 2
                        ! interpolation filter - copy SC indices, interpolate on WC indices
                        buffer_out(ix) = u( ix, iy, iz, ic)
                        buffer_out(ix+1) = u( ix-2, iy, iz, ic) * params%HR(-3) + u( ix  , iy, iz, ic) * params%HR(-1) + &
                                             u( ix+2, iy, iz, ic) * params%HR( 1) + u( ix+4, iy, iz, ic) * params%HR( 3)
                        ! lifted WR filter - indices 0,1,2,3,4,6 only half is active at any time
                        buffer_out(ix) = buffer_out(ix) + u( ix-3, iy, iz, ic) * params%GR(-3) + u( ix-1, iy, iz, ic) * params%GR(-1) + u( ix+1, iy, iz, ic) * params%GR( 1) + u( ix+3, iy, iz, ic) * params%GR( 3)
                        buffer_out(ix+1) = buffer_out(ix+1) + u( ix+1, iy, iz, ic) * params%GR(0) &
                            + u( ix-5, iy, iz, ic) * params%GR(-6) + u( ix-3, iy, iz, ic) * params%GR(-4) + u( ix-1, iy, iz, ic) * params%GR(-2) &
                            + u( ix+3, iy, iz, ic) * params%GR( 2) + u( ix+5, iy, iz, ic) * params%GR( 4) + u( ix+7, iy, iz, ic) * params%GR( 6)
                    enddo
                    u_r(g(1)+1:Bs(1)+g(1),iy,iz,ic) = buffer_out(g(1)+1:Bs(1)+g(1))  ! insert values from buffer to avoid overwriting
                enddo; enddo
                ! ~~~~~~~~~~~~~~~~~~~~~~ Y ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do ix = g(1)+1, Bs(1)+g(1)
                    do iy = g(2)+1, Bs(2)+g(2), 2
                        ! interpolation filter - copy SC indices, interpolate on WC indices
                        buffer_out(iy) = u_r( ix, iy, iz, ic)
                        buffer_out(iy+1) = u_r( ix, iy-2, iz, ic) * params%HR(-3) + u_r( ix, iy  , iz, ic) * params%HR(-1) + &
                                           u_r( ix, iy+2, iz, ic) * params%HR( 1) + u_r( ix, iy+4, iz, ic) * params%HR( 3)
                        ! lifted WR filter - indices 0,1,2,3,4,6 only half is active at any time
                        buffer_out(iy) = buffer_out(iy) + u_r( ix, iy-3, iz, ic) * params%GR(-3) + u_r( ix, iy-1, iz, ic) * params%GR(-1) + u_r( ix, iy+1, iz, ic) * params%GR( 1) + u_r( ix, iy+3, iz, ic) * params%GR( 3)
                        buffer_out(iy+1) = buffer_out(iy+1) + u_r( ix, iy+1, iz, ic) * params%GR(0) &
                            + u_r( ix, iy-5, iz, ic) * params%GR(-6) + u_r( ix, iy-3, iz, ic) * params%GR(-4) + u_r( ix, iy-1, iz, ic) * params%GR(-2) &
                            + u_r( ix, iy+3, iz, ic) * params%GR( 2) + u_r( ix, iy+5, iz, ic) * params%GR( 4) + u_r( ix, iy+7, iz, ic) * params%GR( 6)
                    enddo
                    u_r(ix,g(2)+1:Bs(2)+g(2),iz,ic) = buffer_out(g(2)+1:Bs(2)+g(2))  ! insert values from buffer to avoid overwriting
                enddo; enddo
                ! ~~~~~~~~~~~~~~~~~~~~~~ Z ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if (params%dim == 3) then
                    do iy = g(2)+1, Bs(2)+g(2); do ix = g(1)+1, Bs(1)+g(1)
                        do iz = g(3)+1, Bs(3)+g(3), 2
                            ! interpolation filter - copy SC indices, interpolate on WC indices
                            buffer_out(iz) = u_r( ix, iy, iz, ic)
                            buffer_out(iz+1) = u_r( ix, iy, iz-2, ic) * params%HR(-3) + u_r( ix, iy, iz  , ic) * params%HR(-1) + &
                                               u_r( ix, iy, iz+2, ic) * params%HR( 1) + u_r( ix, iy, iz+4, ic) * params%HR( 3)
                            ! lifted WR filter - indices 0,1,2,3,4,6 only half is active at any time
                            buffer_out(iz) = buffer_out(iz) + u_r( ix, iy, iz-3, ic) * params%GR(-3) + u_r( ix, iy, iz-1, ic) * params%GR(-1) + u_r( ix, iy, iz+1, ic) * params%GR( 1) + u_r( ix, iy, iz+3, ic) * params%GR( 3)
                            buffer_out(iz+1) = buffer_out(iz+1) + u_r( ix, iy, iz+1, ic) * params%GR(0) &
                                + u_r( ix, iy, iz-5, ic) * params%GR(-6) + u_r( ix, iy, iz-3, ic) * params%GR(-4) + u_r( ix, iy, iz-1, ic) * params%GR(-2) &
                                + u_r( ix, iy, iz+3, ic) * params%GR( 2) + u_r( ix, iy, iz+5, ic) * params%GR( 4) + u_r( ix, iy, iz+7, ic) * params%GR( 6)
                        enddo
                        u_r(ix,iy,g(3)+1:Bs(3)+g(3),ic) = buffer_out(g(3)+1:Bs(3)+g(3))  ! insert values from buffer to avoid overwriting
                    enddo; enddo
                endif

            case ("CDF62")
                ! ~~~~~~~~~~~~~~~~~~~~~~ X ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do iy = g(2)+1-f_size_max(2), Bs(2)+g(2)+f_size_max(2)
                    do ix = g(1)+1, Bs(1)+g(1), 2
                        ! interpolation filter - copy SC indices, interpolate on WC indices
                        buffer_out(ix) = u( ix, iy, iz, ic)
                        buffer_out(ix+1) = u( ix-4, iy, iz, ic) * params%HR(-5) + u( ix-2, iy, iz, ic) * params%HR(-3) + u( ix  , iy, iz, ic) * params%HR(-1) &
                                           + u( ix+2, iy, iz, ic) * params%HR( 1) + u( ix+4, iy, iz, ic) * params%HR( 3) + u( ix+6, iy, iz, ic) * params%HR( 5)
                        ! lifted WR filter - indices 0,1,2,4,6 only half is active at any time
                        buffer_out(ix) = buffer_out(ix) + u( ix-1, iy, iz, ic) * params%GR(-1) + u( ix+1, iy, iz, ic) * params%GR( 1)
                        buffer_out(ix+1) = buffer_out(ix+1) + u( ix+1, iy, iz, ic) * params%GR(0) &
                            + u( ix-5, iy, iz, ic) * params%GR(-6) + u( ix-3, iy, iz, ic) * params%GR(-4) + u( ix-1, iy, iz, ic) * params%GR(-2) &
                            + u( ix+3, iy, iz, ic) * params%GR( 2) + u( ix+5, iy, iz, ic) * params%GR( 4) + u( ix+7, iy, iz, ic) * params%GR( 6)
                    enddo
                    u_r(g(1)+1:Bs(1)+g(1),iy,iz,ic) = buffer_out(g(1)+1:Bs(1)+g(1))  ! insert values from buffer to avoid overwriting
                enddo; enddo
                ! ~~~~~~~~~~~~~~~~~~~~~~ Y ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do ix = g(1)+1, Bs(1)+g(1)
                    do iy = g(2)+1, Bs(2)+g(2), 2
                        ! interpolation filter - copy SC indices, interpolate on WC indices
                        buffer_out(iy) = u_r( ix, iy, iz, ic)
                        buffer_out(iy+1) = u_r( ix, iy-4, iz, ic) * params%HR(-5) + u_r( ix, iy-2, iz, ic) * params%HR(-3) + u_r( ix, iy  , iz, ic) * params%HR(-1) &
                                         + u_r( ix, iy+2, iz, ic) * params%HR( 1) + u_r( ix, iy+4, iz, ic) * params%HR( 3) + u_r( ix, iy+6, iz, ic) * params%HR( 5)
                        ! lifted WR filter - indices 0,1,2,4,6 only half is active at any time
                        buffer_out(iy) = buffer_out(iy) + u_r( ix, iy-1, iz, ic) * params%GR(-1) + u_r( ix, iy+1, iz, ic) * params%GR( 1)
                        buffer_out(iy+1) = buffer_out(iy+1) + u_r( ix, iy+1, iz, ic) * params%GR(0) &
                            + u_r( ix, iy-5, iz, ic) * params%GR(-6) + u_r( ix, iy-3, iz, ic) * params%GR(-4) + u_r( ix, iy-1, iz, ic) * params%GR(-2) &
                            + u_r( ix, iy+3, iz, ic) * params%GR( 2) + u_r( ix, iy+5, iz, ic) * params%GR( 4) + u_r( ix, iy+7, iz, ic) * params%GR( 6)
                    enddo
                    u_r(ix,g(2)+1:Bs(2)+g(2),iz,ic) = buffer_out(g(2)+1:Bs(2)+g(2))  ! insert values from buffer to avoid overwriting
                enddo; enddo
                ! ~~~~~~~~~~~~~~~~~~~~~~ Z ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if (params%dim == 3) then
                    do iy = g(2)+1, Bs(2)+g(2); do ix = g(1)+1, Bs(1)+g(1)
                        do iz = g(3)+1, Bs(3)+g(3), 2
                            ! interpolation filter - copy SC indices, interpolate on WC indices
                            buffer_out(iz) = u_r( ix, iy, iz, ic)
                            buffer_out(iz+1) = u_r( ix, iy, iz-4, ic) * params%HR(-5) + u_r( ix, iy, iz-2, ic) * params%HR(-3) + u_r( ix, iy, iz  , ic) * params%HR(-1) &
                                             + u_r( ix, iy, iz+2, ic) * params%HR( 1) + u_r( ix, iy, iz+4, ic) * params%HR( 3) + u_r( ix, iy, iz+6, ic) * params%HR( 5)
                            ! lifted WR filter - indices 0,1,2,4,6 only half is active at any time
                            buffer_out(iz) = buffer_out(iz) + u_r( ix, iy, iz-1, ic) * params%GR(-1) + u_r( ix, iy, iz+1, ic) * params%GR( 1)
                            buffer_out(iz+1) = buffer_out(iz+1) + u_r( ix, iy, iz+1, ic) * params%GR(0) &
                                + u_r( ix, iy, iz-5, ic) * params%GR(-6) + u_r( ix, iy, iz-3, ic) * params%GR(-4) + u_r( ix, iy, iz-1, ic) * params%GR(-2) &
                                + u_r( ix, iy, iz+3, ic) * params%GR( 2) + u_r( ix, iy, iz+5, ic) * params%GR( 4) + u_r( ix, iy, iz+7, ic) * params%GR( 6)
                        enddo
                        u_r(ix,iy,g(3)+1:Bs(3)+g(3),ic) = buffer_out(g(3)+1:Bs(3)+g(3))  ! insert values from buffer to avoid overwriting
                    enddo; enddo
                endif

            ! Lazy cases - some wavelets are not often used and have not been hardcoded for performance, hence we use a general algorithm
            case default
                ! ~~~~~~~~~~~~~~~~~~~~~~ X ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do iy = g(2)+1-f_size_max(2), Bs(2)+g(2)+f_size_max(2)
                    ! fill upsampling buffer for low-pass filter: every second point
                    buffer = 0.0_rk
                    buffer(1+io(1):nx:2) = u(1+io(1):Bs(1)+2*g(1):2, iy, iz, ic) ! SC
                    do ix = g(1)+1, Bs(1)+g(1)
                        buffer_out(ix) = sum(buffer(ix+Hs(1):ix+Hs(2)) * params%HR(:))
                    enddo
                    ! fill upsampling buffer for high-pass filter: every second point
                    buffer = 0.0_rk
                    buffer(2-io(1):nx:2) = u(2-io(1):Bs(1)+2*g(1):2, iy, iz, ic) ! WC
                    do ix = g(1)+1, Bs(1)+g(1)
                        u_r(ix, iy, iz, ic) = buffer_out(ix) + sum(buffer(ix+Gs(1):ix+Gs(2)) * params%GR(:))
                    enddo
                enddo; enddo

                ! ~~~~~~~~~~~~~~~~~~~~~~ Y ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                ! ignore ghost points for dimensions that were already treated
                do iz = g(3)+1-f_size_max(3), Bs(3)+g(3)+f_size_max(3); do ix = g(1)+1, Bs(1)+g(1)
                    ! fill upsampling buffer for low-pass filter: every second point
                    buffer = 0.0_rk
                    buffer(1+io(2):ny:2) = u_r(ix, 1+io(2):Bs(2)+2*g(2):2, iz, ic) ! SC
                    do iy = g(2)+1, Bs(2)+g(2)
                        buffer_out(iy) = sum(buffer(iy+Hs(1):iy+Hs(2)) * params%HR(:))
                    enddo
                    ! fill upsampling buffer for high-pass filter: every second point
                    buffer = 0.0_rk
                    buffer(2-io(2):ny:2) = u_r(ix, 2-io(2):Bs(2)+2*g(2):2, iz, ic) ! WC
                    do iy = g(2)+1, Bs(2)+g(2)
                        u_r(ix,iy,iz,ic) = buffer_out(iy) + sum(buffer(iy+Gs(1):iy+Gs(2)) * params%GR(:))
                    enddo
                enddo; enddo


                ! ~~~~~~~~~~~~~~~~~~~~~~ Z ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                if (params%dim == 3) then
                    ! ignore ghost points for dimensions that were already treated
                    do iy = g(2)+1, Bs(2)+g(2); do ix = g(1)+1, Bs(1)+g(1)
                        ! fill upsampling buffer for low-pass filter: every second point
                        buffer = 0.0_rk
                        buffer(1+io(3):nz:2) = u_r(ix, iy, 1+io(3):Bs(3)+2*g(3):2, ic) ! SC
                        do iz = g(3)+1, Bs(3)+g(3)
                            buffer_out(iz) = sum(buffer(iz+Hs(1):iz+Hs(2)) * params%HR(:))
                        enddo
                        ! fill upsampling buffer for high-pass filter: every second point
                        buffer = 0.0_rk
                        buffer(2-io(3):nz:2) = u_r(ix, iy, 2-io(3):Bs(3)+2*g(3):2, ic) ! WC
                        do iz = g(3)+1, Bs(3)+g(3)
                            u_r(ix,iy,iz,ic) = buffer_out(iz) + sum(buffer(iz+Gs(1):iz+Gs(2)) * params%GR(:))
                        enddo
                    enddo; enddo
                endif

            end select

        enddo

    end subroutine