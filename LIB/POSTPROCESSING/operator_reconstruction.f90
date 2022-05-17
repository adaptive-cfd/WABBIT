subroutine operator_reconstruction(params)
    use module_precision
    use module_mesh
    use module_params
    use module_IO
    use module_forest
    use module_mpi
    use module_acm
    use module_time_step, only: filter_wrapper

    implicit none

    type (type_params), intent(inout)  :: params
    character(len=cshort) :: file, infile
    real(kind=rk) :: time, x, y, dx_fine, u_dx, u_dxdx, dx_inv, val, x2, y2, nu, x_in, y_in
    integer(kind=ik) :: iteration, k, lgt_id, tc_length, tree_N, iblock, ix, iy, &
    g, lgt_n, hvy_n, iz, a1, b1, a2, b2, level
    integer(kind=ik) :: ixx, iyy, ix2, iy2, nx_fine, ixx2,iyy2, n_nonzero
    integer(kind=ik), dimension(3) :: Bs
    character(len=2)       :: order

    integer(kind=ik), allocatable      :: lgt_block(:, :)
    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_work(:, :, :, :, :, :), hvy_tmp(:, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_mask(:, :, :, :, :)
    real(kind=rk), allocatable         :: stencil1(:),stencil2(:)
    integer(kind=ik), allocatable      :: hvy_neighbor(:,:)
    integer(kind=ik), allocatable      :: lgt_active(:), hvy_active(:)
    integer(kind=tsize), allocatable   :: lgt_sortednumlist(:,:)
    character(len=cshort)              :: fname
    real(kind=rk), dimension(3)        :: dx, x0
    integer(hid_t)                     :: file_id
    real(kind=rk), dimension(3)        :: domain
    character(len=1) :: dir
    logical :: refine, notAllPointsAreZero

    ! Tam & Webb, 4th order optimized (for first derivative)
    ! a = (/-0.02651995_rk, +0.18941314_rk, -0.79926643_rk, 0.0_rk, 0.79926643_rk, -0.18941314_rk, 0.02651995_rk/)
    ! standard 4th central FD stencil
    ! a = (/0.0_rk , 1.0_rk/12.0_rk, -2.0_rk/3.0_rk, 0.0_rk, +2.0_rk/3.0_rk, -1.0_rk/12.0_rk, 0.0_rk/)
    ! 4th order coefficients for second derivative
    ! b = (/ -1.0_rk/12.0_rk, 4.0_rk/3.0_rk, -5.0_rk/2.0_rk, 4.0_rk/3.0_rk, -1.0_rk/12.0_rk /)


    if (params%number_procs>1) call abort(2205121, "OperatorReconstruction is a serial routine...")

    call get_command_argument(2, file)
    call check_file_exists(file)

    ! get some parameters from the grid file
    call read_attributes(file, lgt_n, time, iteration, domain, Bs, tc_length, params%dim, &
    periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)


    params%max_treelevel = tc_length+2 ! to allow refinement
    params%n_eqn = 2
    params%domain_size = domain
    params%Bs = Bs
    allocate(params%butcher_tableau(1,1))
    allocate(params%symmetry_vector_component(1:params%n_eqn))
    params%symmetry_vector_component = "0"
    params%number_blocks = 8*lgt_n ! to allow refinement

    ! Note:
    ! When comparing with the basic hand-made matlab operator script, keep in mind
    ! that the coarseWins solution overwrites the fine with coarse data on the interface.
    ! The basic matlab script does not keep these points, but wabbit does. So at the interface, two
    ! new points are added in wabbit


    ! discretization
    call get_cmd_arg( "--discretization", params%order_discretization, default="FD_4th_central" )
    call get_cmd_arg( "--predictor", params%order_predictor, default="multiresolution_4th" )
    ! viscosity
    call get_cmd_arg( "--viscosity", nu, default=0.0_rk )
    ! coarseWins or fineWins
    call get_cmd_arg( "--coarse-wins", params%ghost_nodes_redundant_point_coarseWins, default=.false. )
    call get_cmd_arg( "--refine-coarsen", refine, default=.false. )
    call get_cmd_arg( "--wavelet", params%wavelet, default="CDF40" )

    !---------------------------------------------------------------------------
    ! Adjustable PARAMETERS
    !---------------------------------------------------------------------------
    dir = "x"
    if (params%wavelet=="CDF40") then
        params%wavelet_transform_type = "harten-multiresolution"
    else
        params%wavelet_transform_type = "biorthogonal"
    endif
    params%iter_ghosts = .false.
    !---------------------------------------------------------------------------


    select case(params%order_discretization)
    case("FD_4th_central_optimized")
        ! Tam & Webb, 4th order optimized (for first derivative)
        allocate(stencil1(-3:+3))
        stencil1 = (/-0.02651995_rk, +0.18941314_rk, -0.79926643_rk, 0.0_rk, 0.79926643_rk, -0.18941314_rk, 0.02651995_rk/)
        ! 2nd derivative
        allocate(stencil2(-2:+2))
        stencil2 = (/-1.0_rk/12.0_rk, 4.0_rk/3.0_rk, -5.0_rk/2.0_rk, 4.0_rk/3.0_rk, -1.0_rk/12.0_rk/)

        params%n_ghosts = 4_ik

    case("FD_4th_central")
        ! standard 4th central FD stencil
        allocate(stencil1(-2:+2))
        stencil1 = (/1.0_rk/12.0_rk, -2.0_rk/3.0_rk, 0.0_rk, +2.0_rk/3.0_rk, -1.0_rk/12.0_rk/)
        ! 2nd derivative
        allocate(stencil2(-2:+2))
        stencil2 = (/-1.0_rk/12.0_rk, 4.0_rk/3.0_rk, -5.0_rk/2.0_rk, 4.0_rk/3.0_rk, -1.0_rk/12.0_rk/)

        params%n_ghosts = 2_ik

    case("FD_2nd_central")
        ! Tam & Webb, 4th order optimized (for first derivative)
        allocate(stencil1(-1:+1))
        stencil1 = (/-0.5_rk, 0.0_rk, 0.5_rk/)
        ! 2nd derivative
        allocate(stencil2(-1:+1))
        stencil2 = (/1.0_rk, -2.0_rk, 1.0_rk/)

        params%n_ghosts = 2_ik

    case default
        call abort(1919191222,"unknown discretization set?!")

    end select

    if (params%wavelet=="CDF44") then
        params%n_ghosts = 6_ik
    endif

    open(17, file=trim(adjustl(file))//'.info.txt', status='replace')
    write(17,'(A,1x,A,1x,A," g=",i1," Bs=",i2, " coarseWins=",L1," nu=",es12.4," refineCoarsen=",L1,1x,A)') trim(params%order_discretization), &
    trim(params%order_predictor), dir, params%n_ghosts, params%Bs(1), params%ghost_nodes_redundant_point_coarseWins, nu, refine, params%wavelet
    close(17)

    !---------------------------------------------------------------------------

    if ((params%order_discretization == "FD_4th_central_optimized").and.(params%n_ghosts<4)) then
        call abort(33,"not enough g")
    endif


    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------


    g = params%n_ghosts
    a1 = lbound(stencil1, dim=1)
    b1 = ubound(stencil1, dim=1)
    a2 = lbound(stencil2, dim=1)
    b2 = ubound(stencil2, dim=1)
    iz = 1

    call allocate_grid(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, hvy_active, lgt_sortednumlist, hvy_tmp=hvy_tmp)

    call init_ghost_nodes( params )

    call read_mesh(file, params, lgt_n, hvy_n, lgt_block)

    call create_active_and_sorted_lists(params, lgt_block, lgt_active, &
    lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_ID=1)

    call update_neighbors(params, lgt_block, hvy_neighbor, lgt_active, &
    lgt_n, lgt_sortednumlist, hvy_active, hvy_n)

    dx_fine = (2.0_rk**-max_active_level(lgt_block, lgt_active, lgt_n))*domain(2)/real((Bs(2)-1), kind=rk)
    nx_fine = nint(domain(2)/dx_fine)

    write(*,*) "nx_fine=", nx_fine
    write(*,*) "nblocks=", lgt_n, "bs=", Bs, "npoints (op. matrix size!)=", lgt_n*bs(1)*bs(2)

    ! this hack ensures, on mono_CPU, that later on, refine+coarsening always ends up in the same order in hvy_actve
    if (refine) then
        call sync_ghosts(params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n)
        call refine_mesh( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, &
        lgt_sortednumlist, hvy_active, hvy_n, "everywhere", tree_ID=1 )

        call sync_ghosts(params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n)
        call adapt_mesh( time, params, lgt_block, hvy_block, hvy_neighbor, lgt_active, &
        lgt_n, lgt_sortednumlist, hvy_active, hvy_n, tree_ID_flow, "everywhere", hvy_tmp, external_loop=.true. )
    endif

    !---------------------------------------------------------------------------
    ! save the grid (for plotting in python)
    !---------------------------------------------------------------------------
    open(19, file=trim(adjustl(file))//'.operator_grid_points.txt', status='replace')
    do iblock = 1, hvy_n
        call hvy2lgt(lgt_id, hvy_active(iblock), params%rank, params%number_blocks)
        call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )
        level = lgt_block(lgt_id, params%max_treelevel+IDX_MESH_LVL)
        do ix = g+1, Bs(1)+g
            do iy = g+1, Bs(2)+g
                x = dble(ix-(g+1)) * dx(1) + x0(1)
                y = dble(iy-(g+1)) * dx(2) + x0(2)

                ixx = nint(x/dx_fine)+1
                iyy = nint(y/dx_fine)+1

                write(19,*) ixx, iyy, x, y, level
            enddo
        enddo
    enddo
    close(19)

    !---------------------------------------------------------------------------
    ! compute operator matrix
    !---------------------------------------------------------------------------
    open(17, file=trim(adjustl(file))//'.operator_matrix.txt', status='replace')
    ! open(18, file=trim(adjustl(file))//'.operator_matrix2.txt', status='replace')

    ! loop over points on the fine level (expensive - loops over a lot of points that do not exist...)
    do iyy = 1, nx_fine!+1 ! skip periodic up/right
        do ixx = 1, nx_fine!+1
    ! do ixx = 2, nx_fine ! skip awful periodic points
    !     do iyy = 2, nx_fine

            ! coordinates of the point to be set to one
            x_in = dble(ixx-1) * dx_fine
            y_in = dble(iyy-1) * dx_fine

            ! if (abs((x_in-domain(1))) <= 1.0e-12) x_in = 0.0_rk
            ! if (abs((y_in-domain(2))) <= 1.0e-12) y_in = 0.0_rk

            notAllPointsAreZero = .false.

            ! reset entire grid to zeros (we just set one point to 1.0)
            do iblock = 1, hvy_n
                hvy_block(:, :, :, :, hvy_active(iblock)) = 0.0_rk
            enddo

            ! set this one point we're looking at to 1 (on all blocks!)
            do iblock = 1, hvy_n
                call hvy2lgt(lgt_id, hvy_active(iblock), params%rank, params%number_blocks)
                call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

                if ((x_in >= x0(1)) .and. (x_in <= x0(1)+dble(Bs(1)-1)*dx(1) )) then
                    if ((y_in >= x0(2)) .and. (y_in <= x0(2)+dble(Bs(2)-1)*dx(2) )) then
                        ! the "one" lies on this block
                        ix = nint((x_in-x0(1))/dx(1)) + g+1
                        iy = nint((y_in-x0(2))/dx(2)) + g+1

                        x = dble(ix-(g+1)) * dx(1) + x0(1)
                        y = dble(iy-(g+1)) * dx(2) + x0(2)

                        ! is this *really* the right point? Think of two levels, fine & coarse. you strive to set 1 on the finest grid point,
                        ! but on the coarse one, it lies exactly between two coarser points. neither one is correct ! so, check again, if (ix,iy)
                        ! really correspond to (x_in,y_in)
                        if ((abs(x-x_in)<1.0e-13_rk) .and. (abs(y-y_in)<1.0e-13_rk)) then
                            hvy_block(ix, iy, iz, :, hvy_active(iblock)) = 1.0_rk
                            notAllPointsAreZero = .true.
                            ! write(*,*) x_in, y_in, x, y, "--", ix,iy,iblock
                        endif
                    endif
                endif

                if ((x_in+1.0_rk >= x0(1)) .and. (x_in+1.0_rk <= x0(1)+dble(Bs(1)-1)*dx(1) )) then
                    if ((y_in >= x0(2)) .and. (y_in <= x0(2)+dble(Bs(2)-1)*dx(2) )) then
                        ! the "one" lies on this block
                        ix = nint((x_in+1.0_rk-x0(1))/dx(1)) + g+1
                        iy = nint((y_in-x0(2))/dx(2)) + g+1

                        x = dble(ix-(g+1)) * dx(1) + x0(1)
                        y = dble(iy-(g+1)) * dx(2) + x0(2)

                        ! is this *really* the right point? Think of two levels, fine & coarse. you strive to set 1 on the finest grid point,
                        ! but on the coarse one, it lies exactly between two coarser points. neither one is correct ! so, check again, if (ix,iy)
                        ! really correspond to (x_in,y_in)
                        if ((abs(x-(x_in+1.0_rk))<1.0e-13_rk) .and. (abs(y-y_in)<1.0e-13_rk)) then
                            hvy_block(ix, iy, iz, :, hvy_active(iblock)) = 1.0_rk
                            notAllPointsAreZero = .true.
                            ! write(*,*) x_in+1.0_rk, y_in, x, y, "--", ix,iy,iblock
                        endif
                    endif
                endif
                if ((x_in >= x0(1)) .and. (x_in <= x0(1)+dble(Bs(1)-1)*dx(1) )) then
                    if ((y_in+1.0_rk >= x0(2)) .and. (y_in+1.0_rk <= x0(2)+dble(Bs(2)-1)*dx(2) )) then
                        ! the "one" lies on this block
                        ix = nint((x_in-x0(1))/dx(1)) + g+1
                        iy = nint((y_in+1.0_rk-x0(2))/dx(2)) + g+1

                        x = dble(ix-(g+1)) * dx(1) + x0(1)
                        y = dble(iy-(g+1)) * dx(2) + x0(2)

                        ! is this *really* the right point? Think of two levels, fine & coarse. you strive to set 1 on the finest grid point,
                        ! but on the coarse one, it lies exactly between two coarser points. neither one is correct ! so, check again, if (ix,iy)
                        ! really correspond to (x_in,y_in)
                        if ((abs(x-x_in)<1.0e-13_rk) .and. (abs(y-(y_in+1.0_rk))<1.0e-13_rk)) then
                            hvy_block(ix, iy, iz, :, hvy_active(iblock)) = 1.0_rk
                            notAllPointsAreZero = .true.
                            ! write(*,*) x_in, y_in+1.0_rk, x, y, "--", ix,iy,iblock
                        endif
                    endif
                endif
                if ((x_in+1.0_rk >= x0(1)) .and. (x_in+1.0_rk <= x0(1)+dble(Bs(1)-1)*dx(1) )) then
                    if ((y_in+1.0_rk >= x0(2)) .and. (y_in+1.0_rk <= x0(2)+dble(Bs(2)-1)*dx(2) )) then
                        ! the "one" lies on this block
                        ix = nint((x_in+1.0_rk-x0(1))/dx(1)) + g+1
                        iy = nint((y_in+1.0_rk-x0(2))/dx(2)) + g+1

                        x = dble(ix-(g+1)) * dx(1) + x0(1)
                        y = dble(iy-(g+1)) * dx(2) + x0(2)

                        ! is this *really* the right point? Think of two levels, fine & coarse. you strive to set 1 on the finest grid point,
                        ! but on the coarse one, it lies exactly between two coarser points. neither one is correct ! so, check again, if (ix,iy)
                        ! really correspond to (x_in,y_in)
                        if ((abs(x-(x_in+1.0_rk))<1.0e-13_rk) .and. (abs(y-(y_in+1.0_rk))<1.0e-13_rk)) then
                            hvy_block(ix, iy, iz, :, hvy_active(iblock)) = 1.0_rk
                            notAllPointsAreZero = .true.
                            ! write(*,*) x_in+1.0_rk, y_in+1.0_rk, x, y, "--", ix,iy,iblock
                        endif
                    endif
                endif
            enddo

            ! many points on the finest level do not exist - if the entire grid is zeros,
            ! then we can cycle here.
            if (.not. notAllPointsAreZero) cycle

            !---------------------------------------------------------------
            ! synchronize ghosts (important! if e.g. coarseWins is active and you happen to set the redundant value of a refined block, its overwritten to be zero again)
            ! Note: this also applies to coarse block bordering on a coarse block, if its ID is lower.
            ! In fact, each point is then computed only once. Note: if you set the point on a high lgt_id, then
            ! it will be "sync'ed down" to lower light IDs, so you can find the point more than once
            ! NOTE: with the new version, where the outer loops are ixx and iyy (so we loop over the grid point on
            ! the finest active level), this sync is no longer required: we set ALL points (x_in,y_in) to one, on ALL blocks
            ! that have those points
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call sync_ghosts(params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            if (refine) then
                call refine_mesh( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, &
                lgt_sortednumlist, hvy_active, hvy_n, "everywhere", tree_ID=1 )
            endif

            call sync_ghosts(params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n)

            !---------------------------------------------------------------
            ! Now compute the derivative
            do k = 1, hvy_n
                call hvy2lgt(lgt_id, hvy_active(k), params%rank, params%number_blocks)
                call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

                dx_inv = 1.0_rk / dx(1)

                if (dir=="x") then
                    do iy2 = g+1, Bs(2)+g
                        do ix2 = g+1, Bs(1)+g
                            u_dx   = sum( stencil1*hvy_block(ix2+a1:ix2+b1, iy2, iz, 1, hvy_active(k)) )*dx_inv
                            u_dxdx = sum( stencil2*hvy_block(ix2+a2:ix2+b2, iy2, iz, 1, hvy_active(k)) )*dx_inv**2
                            hvy_block(ix2,iy2,iz,2,hvy_active(k)) = u_dx + nu*u_dxdx
                        end do
                    end do
                elseif (dir=='y') then
                    do iy2 = g+1, Bs(2)+g
                        do ix2 = g+1, Bs(1)+g
                            u_dx   = sum( stencil1*hvy_block(ix2, iy2+a1:iy2+b1, iz, 1, hvy_active(k)) )*dx_inv
                            u_dxdx = sum( stencil2*hvy_block(ix2, iy2+a2:iy2+b2, iz, 1, hvy_active(k)) )*dx_inv**2
                            hvy_block(ix2,iy2,iz,2,hvy_active(k)) = u_dx + nu*u_dxdx
                        end do
                    end do
                else
                    call abort(123,'X or Y baby, nothing else.')
                endif
            end do

            ! This second sync step also synchronizes the derivative we computed previously
            ! Note that on a coarse/fine interface, wabbit computes two values for the derivative
            ! on the coarse and fine level. This synchronizing step lets us keep only either of those,
            ! depending on fineWins or coarseWins
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call sync_ghosts(params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            if (refine) then
                call adapt_mesh( time, params, lgt_block, hvy_block, hvy_neighbor, lgt_active, &
                lgt_n, lgt_sortednumlist, hvy_active, hvy_n, tree_ID_flow, "everywhere", hvy_tmp, external_loop=.true. )
            endif

            call sync_ghosts(params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n)

            ! save operator line to text file.
            ! note: unfortunately, we use the index on the finest level, i.e., we temporarily
            ! create a matrix N_max**2 by N_max**2, where N_max is Npoints on the finest level.
            ! Many points do not exist; they are on coarse levels. however, this is a problem
            ! for the python script, because it first reads the entie matrix, then removes zero cols/rows.
            do k = 1, hvy_n
                call hvy2lgt(lgt_id, hvy_active(k), params%rank, params%number_blocks)
                call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

                do iy2 = g+1, Bs(2)+g
                    do ix2 = g+1, Bs(1)+g
                        x = dble(ix2-(g+1)) * dx(1) + x0(1)
                        y = dble(iy2-(g+1)) * dx(2) + x0(2)

                        if (abs((x-domain(1))) <=1.0e-9) x = 0.0_rk
                        if (abs((y-domain(2))) <=1.0e-9) y = 0.0_rk

                        ixx2 = nint(x/dx_fine)+1
                        iyy2 = nint(y/dx_fine)+1

                        val = hvy_block(ix2, iy2, iz, 2, hvy_active(k)) ! u_dx

                        if (abs(val) > 1.0e-13) then
                            ! this point is a nonzero value
                            write(17,'(i6,1x,i6,1x,es15.8,1x,4(es15.7,1x))') ixx+(iyy-1)*nx_fine, ixx2+(iyy2-1)*nx_fine, val, x_in, y_in, x, y
                            ! write(17,'(i6,1x,i6,1x,es15.8)') ixx+(iyy-1)*nx_fine, ixx2+(iyy2-1)*nx_fine, val
                            ! write(*,*) "python col=", ixx2+(iyy2-1)*nx_fine -1 , "row=", ixx+(iyy-1)*nx_fine -1, val, "xy=",x,y, "block", k, ix,iy, Bs(1)+g, Bs(2)+g
                            ! write(*,*) x_in, y_in, "---", x, y, "---", val
                        endif
                    end do
                end do
            end do


        enddo
    enddo
    close(17)
    ! close(18)
end subroutine
