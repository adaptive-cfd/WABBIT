subroutine operator_reconstruction(params)
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_acm
    use module_time_step, only: filter_wrapper
    use module_forestMetaData

    implicit none

    type (type_params), intent(inout)  :: params
    character(len=cshort) :: file, infile
    real(kind=rk) :: time, x, y, dx_fine, u_dx, u_dxdx, dx_inv, val, x2, y2, nu, x_in, y_in
    integer(kind=ik) :: iteration, k, lgt_id, tc_length, iblock, ix, iy, &
    g, iz, a1, b1, a2, b2, level, j
    integer(kind=ik) :: ixx, iyy, ix2, iy2, nx_fine, ixx2,iyy2, n_nonzero
    integer(kind=ik), dimension(3) :: Bs
    character(len=2)       :: order

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_work(:, :, :, :, :, :), hvy_tmp(:, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_mask(:, :, :, :, :)
    real(kind=rk), allocatable         :: stencil1(:),stencil2(:)

    integer(kind=ik)                   :: tree_ID=1, hvy_id

    character(len=cshort)              :: fname
    real(kind=rk), dimension(3)        :: dx, x0
    integer(hid_t)                     :: file_id
    real(kind=rk), dimension(3)        :: domain
    character(len=1) :: dir
    logical :: refine, notAllPointsAreZero
    real(kind=rk), dimension(4), PARAMETER :: permut_x=(/0.0_rk, 1.0_rk, 0.0_rk, 1.0_rk/), permut_y=(/0.0_rk, 0.0_rk, 1.0_rk, 1.0_rk/)
    real(kind=rk) :: x_in2, y_in2

    if (params%number_procs>1) call abort(2205121, "OperatorReconstruction is a serial routine...")

    ! this routine works only on one tree
    allocate( hvy_n(1), lgt_n(1) )

    call get_command_argument(2, file)
    call check_file_exists(file)

    ! get some parameters from the grid file
    call read_attributes(file, lgt_n(tree_ID), time, iteration, domain, Bs, tc_length, params%dim, &
    periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)


    params%Jmax = tc_length+2 ! to allow refinement
    params%n_eqn = 2
    params%domain_size = domain
    params%Bs = Bs
    allocate(params%butcher_tableau(1,1))
    allocate(params%symmetry_vector_component(1:params%n_eqn))
    params%symmetry_vector_component = "0"
    params%number_blocks = ceiling(8.0_rk*lgt_n(tree_ID) * 2.0_rk**params%dim / (2.0_rk**params%dim - 1.0_rk)) ! to allow refinement and adapt_tree

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
    call get_cmd_arg( "--refine-coarsen", refine, default=.false. )
    call get_cmd_arg( "--wavelet", params%wavelet, default="CDF40" )

    !---------------------------------------------------------------------------
    ! Adjustable PARAMETERS
    !---------------------------------------------------------------------------
    dir = "x"
    !---------------------------------------------------------------------------


    select case(params%order_discretization)
    case("FD_4th_central_optimized")
        ! Tam & Webb, 4th order optimized (for first derivative)
        allocate(stencil1(-3:+3))
        stencil1 = (/-0.02651995_rk, +0.18941314_rk, -0.79926643_rk, 0.0_rk, 0.79926643_rk, -0.18941314_rk, 0.02651995_rk/)
        ! 2nd derivative
        allocate(stencil2(-2:+2))
        stencil2 = (/-1.0_rk/12.0_rk, 4.0_rk/3.0_rk, -5.0_rk/2.0_rk, 4.0_rk/3.0_rk, -1.0_rk/12.0_rk/)

        params%g = 4_ik

    case("FD_4th_central")
        ! standard 4th central FD stencil
        allocate(stencil1(-2:+2))
        stencil1 = (/1.0_rk/12.0_rk, -2.0_rk/3.0_rk, 0.0_rk, +2.0_rk/3.0_rk, -1.0_rk/12.0_rk/)
        ! 2nd derivative
        allocate(stencil2(-2:+2))
        stencil2 = (/-1.0_rk/12.0_rk, 4.0_rk/3.0_rk, -5.0_rk/2.0_rk, 4.0_rk/3.0_rk, -1.0_rk/12.0_rk/)

        params%g = 2_ik

    case("FD_2nd_central")
        ! Tam & Webb, 4th order optimized (for first derivative)
        allocate(stencil1(-1:+1))
        stencil1 = (/-0.5_rk, 0.0_rk, 0.5_rk/)
        ! 2nd derivative
        allocate(stencil2(-1:+1))
        stencil2 = (/1.0_rk, -2.0_rk, 1.0_rk/)

        params%g = 2_ik

    case default
        call abort(1919191222,"unknown discretization set?!")

    end select

    if (params%wavelet=="CDF44") then
        params%g = 6_ik
    endif

    open(17, file=trim(adjustl(file))//'.info.txt', status='replace')
    write(17,'(A,1x,A,1x,A," g=",i1," Bs=",i2," nu=",es12.4," refineCoarsen=",L1,1x,A)') trim(params%order_discretization), &
    trim(params%order_predictor), dir, params%g, params%Bs(1), nu, refine, params%wavelet
    close(17)

    !---------------------------------------------------------------------------

    if ((params%order_discretization == "FD_4th_central_optimized").and.(params%g<4)) then
        call abort(33,"not enough g")
    endif


    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------
    !---------------------------------------------------------------------------


    g = params%g
    a1 = lbound(stencil1, dim=1)
    b1 = ubound(stencil1, dim=1)
    a2 = lbound(stencil2, dim=1)
    b2 = ubound(stencil2, dim=1)
    iz = 1

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas
    call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp)

    call init_ghost_nodes( params )

    ! read data
    call readHDF5vct_tree( (/file/), params, hvy_block, tree_ID)

    call updateMetadata_tree(params, tree_ID)

    dx_fine = (2.0_rk**-maxActiveLevel_tree(tree_ID))*domain(2)/real((Bs(2)-1), kind=rk)
    nx_fine = nint(domain(2)/dx_fine)

    write(*,*) "nx_fine=", nx_fine
    write(*,*) "nblocks=", lgt_n(tree_ID), "bs=", Bs, "npoints (op. matrix size!)=", lgt_n(tree_ID)*bs(1)*bs(2)

    !---------------------------------------------------------------------------
    ! save the grid (for plotting in python)
    !---------------------------------------------------------------------------
    open(19, file=trim(adjustl(file))//'.operator_grid_points.txt', status='replace')
    do iblock = 1, hvy_n(tree_ID)

        call hvy2lgt(lgt_id, hvy_active(iblock, tree_ID), params%rank, params%number_blocks)
        call get_block_spacing_origin( params, lgt_id, x0, dx )

        level = lgt_block(lgt_id, IDX_MESH_LVL)

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

            ! coordinates of the point to be set to one
            x_in = dble(ixx-1) * dx_fine
            y_in = dble(iyy-1) * dx_fine

            ! if (abs((x_in-domain(1))) <=1.0e-9) x_in = 0.0_rk
            ! if (abs((y_in-domain(2))) <=1.0e-9) y_in = 0.0_rk

            notAllPointsAreZero = .false.

            ! reset entire grid to zeros (we just set one point to 1.0)
            do iblock = 1, hvy_n(tree_ID)
                hvy_block(:, :, :, :, hvy_active(iblock, tree_ID)) = 0.0_rk
            enddo

            ! set this one point we're looking at to 1 (on all blocks!)
            do iblock = 1, hvy_n(tree_ID)
                call hvy2lgt(lgt_id, hvy_active(iblock, tree_ID), params%rank, params%number_blocks)
                call get_block_spacing_origin( params, lgt_id, x0, dx )

                do j = 1, 4
                    x_in2 = x_in + permut_x(j)*domain(1)
                    y_in2 = y_in + permut_y(j)*domain(2)

                    if ((x_in2 >= x0(1)) .and. (x_in2 <= x0(1)+dble(Bs(1)-1)*dx(1) )) then
                        if ((y_in2 >= x0(2)) .and. (y_in2 <= x0(2)+dble(Bs(2)-1)*dx(2) )) then
                            ! the "one" lies on this block
                            ix = nint((x_in2-x0(1))/dx(1)) + g+1
                            iy = nint((y_in2-x0(2))/dx(2)) + g+1

                            x = dble(ix-(g+1)) * dx(1) + x0(1)
                            y = dble(iy-(g+1)) * dx(2) + x0(2)

                            ! is this *really* the right point? Think of two levels, fine & coarse. you strive to set 1 on the finest grid point,
                            ! but on the coarse one, it lies exactly between two coarser points. neither one is correct ! so, check again, if (ix,iy)
                            ! really correspond to (x_in2,y_in2)
                            if ((abs(x-x_in2)<1.0e-13_rk) .and. (abs(y-y_in2)<1.0e-13_rk)) then
                                hvy_block(ix, iy, iz, :, hvy_active(iblock,tree_ID)) = 1.0_rk
                                notAllPointsAreZero = .true.
                                ! write(*,*) x_in, y_in, x, y, "--", ix,iy,iblock
                            endif
                        endif
                    endif
                enddo
            enddo

            ! many points on the finest level do not exist - if the entire grid is zeros,
            ! then we can cycle here.
            if (.not. notAllPointsAreZero) cycle

            !---------------------------------------------------------------
            ! synchronize ghosts. This was done for the redundantGrid -> check if still required for the uniqueGrid (as of 2023)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            call sync_ghosts_tree(params, hvy_block, tree_ID)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            if (refine) then
call abort(99999, "need to adapt refine_tree call to include hvy_tmp")
!                call refine_tree( params, hvy_block, "everywhere", tree_ID )

                call sync_ghosts_tree(params, hvy_block, tree_ID)
            endif


            !---------------------------------------------------------------
            ! Now compute the derivative
            do k = 1, hvy_n(tree_ID)
                hvy_id = hvy_active(k, tree_ID)
                call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
                call get_block_spacing_origin( params, lgt_id, x0, dx )

                dx_inv = 1.0_rk / dx(1)

                if (dir=="x") then
                    do iy2 = g+1, Bs(2)+g
                        do ix2 = g+1, Bs(1)+g
                            u_dx   = sum( stencil1*hvy_block(ix2+a1:ix2+b1, iy2, iz, 1, hvy_id) )*dx_inv
                            u_dxdx = sum( stencil2*hvy_block(ix2+a2:ix2+b2, iy2, iz, 1, hvy_id) )*dx_inv**2
                            hvy_block(ix2,iy2,iz,2,hvy_id) = u_dx + nu*u_dxdx
                        end do
                    end do
                elseif (dir=='y') then
                    do iy2 = g+1, Bs(2)+g
                        do ix2 = g+1, Bs(1)+g
                            u_dx   = sum( stencil1*hvy_block(ix2, iy2+a1:iy2+b1, iz, 1, hvy_id) )*dx_inv
                            u_dxdx = sum( stencil2*hvy_block(ix2, iy2+a2:iy2+b2, iz, 1, hvy_id) )*dx_inv**2
                            hvy_block(ix2,iy2,iz,2,hvy_id) = u_dx + nu*u_dxdx
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
            call sync_ghosts_tree(params, hvy_block, tree_ID)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            if (refine) then
                call adapt_tree( time, params, hvy_block, tree_ID_flow, "everywhere", hvy_tmp )

                call sync_ghosts_tree(params, hvy_block, tree_ID)
            endif


            ! save operator line to text file.
            ! note: unfortunately, we use the index on the finest level, i.e., we temporarily
            ! create a matrix N_max**2 by N_max**2, where N_max is Npoints on the finest level.
            ! Many points do not exist; they are on coarse levels
            do k = 1, hvy_n(tree_ID)
                hvy_id = hvy_active(k, tree_ID)

                call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
                call get_block_spacing_origin( params, lgt_id, x0, dx )

                do iy2 = g+1, Bs(2)+g
                    do ix2 = g+1, Bs(1)+g
                        x = dble(ix2-(g+1)) * dx(1) + x0(1)
                        y = dble(iy2-(g+1)) * dx(2) + x0(2)

                        if (abs((x-domain(1))) <=1.0e-9) x = 0.0_rk
                        if (abs((y-domain(2))) <=1.0e-9) y = 0.0_rk

                        ixx2 = nint(x/dx_fine)+1
                        iyy2 = nint(y/dx_fine)+1

                        val = hvy_block(ix2, iy2, iz, 2, hvy_id) ! u_dx

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
