subroutine operator_reconstruction(params)
    use module_precision
    use module_mesh
    use module_params
    use module_IO
    use module_forest
    use module_mpi
    use module_acm
    use module_bridge_interface, only : initialize_communicator
    use module_time_step, only: filter_wrapper

    implicit none

    type (type_params), intent(inout)  :: params
    character(len=cshort) :: file, infile
    real(kind=rk) :: time, x, y, dx_fine, u_dx, dx_inv, val, x2, y2
    integer(kind=ik) :: iteration, k, lgt_id, tc_length, tree_N, iblock, ix, iy, &
    g, lgt_n, hvy_n, iz, a, b
    integer(kind=ik) :: ixx, iyy, ix2, iy2, nx_fine, ixx2,iyy2, n_nonzero
    integer(kind=ik), dimension(3) :: Bs
    character(len=2)       :: order

    integer(kind=ik), allocatable      :: lgt_block(:, :)
    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_work(:, :, :, :, :, :), hvy_tmp(:, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_mask(:, :, :, :, :)
    real(kind=rk), allocatable         :: stencil(:)
    integer(kind=ik), allocatable      :: hvy_neighbor(:,:)
    integer(kind=ik), allocatable      :: lgt_active(:), hvy_active(:)
    integer(kind=tsize), allocatable   :: lgt_sortednumlist(:,:)
    character(len=cshort)                  :: fname
    real(kind=rk), dimension(3)        :: dx, x0
    integer(hid_t)                     :: file_id
    real(kind=rk), dimension(3)        :: domain
    character(len=1) :: dir

    ! Tam & Webb, 4th order optimized (for first derivative)
    ! a = (/-0.02651995_rk, +0.18941314_rk, -0.79926643_rk, 0.0_rk, 0.79926643_rk, -0.18941314_rk, 0.02651995_rk/)
    ! standard 4th central FD stencil
    ! a = (/0.0_rk , 1.0_rk/12.0_rk, -2.0_rk/3.0_rk, 0.0_rk, +2.0_rk/3.0_rk, -1.0_rk/12.0_rk, 0.0_rk/)
    ! 4th order coefficients for second derivative
    ! b = (/ -1.0_rk/12.0_rk, 4.0_rk/3.0_rk, -5.0_rk/2.0_rk, 4.0_rk/3.0_rk, -1.0_rk/12.0_rk /)


    call get_command_argument(2, file)
    call check_file_exists(file)

    ! get some parameters from one of the files (they should be the same in all of them)
    call read_attributes(file, lgt_n, time, iteration, domain, Bs, tc_length, params%dim, &
    periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)


    params%max_treelevel = tc_length
    params%n_eqn = 2
    params%domain_size = domain
    params%Bs = Bs
    allocate(params%butcher_tableau(1,1))
    allocate(params%symmetry_vector_component(1:params%n_eqn))
    params%symmetry_vector_component = "0"
    params%number_blocks = lgt_n

    ! Note:
    ! When comparing with the basic hand-made matlab operator script, keep in mind
    ! that the coarseWins solution overwrites the fine with coarse data on the interface.
    ! The basic matlab script does not keep these points, but wabbit does. So at the interface, two
    ! new points are added in wabbit



    !---------------------------------------------------------------------------
    ! Adjustable PARAMETERS
    !---------------------------------------------------------------------------
    params%n_ghosts = 6_ik
    ! coarseWins or fineWins
    params%ghost_nodes_redundant_point_coarseWins = .false.
    dir = "x"

    ! params%order_discretization = "FD_2nd_central"
    ! params%order_predictor = "multiresolution_2nd"

    params%order_discretization = "FD_4th_central_optimized"
    params%order_predictor = "multiresolution_4th"
    params%wavelet = "CDF44"
    params%wavelet_transform_type = "biorthogonal"!"harten-multiresolution"

    params%iter_ghosts = .true.
    !---------------------------------------------------------------------------


    if (params%order_discretization == "FD_4th_central_optimized") then
        allocate(stencil(-3:+3))
        ! Tam & Webb, 4th order optimized (for first derivative)
        stencil = (/-0.02651995_rk, +0.18941314_rk, -0.79926643_rk, 0.0_rk, 0.79926643_rk, -0.18941314_rk, 0.02651995_rk/)
    endif

    if (params%order_discretization == "FD_4th_central") then
        allocate(stencil(-2:+2))
        ! standard 4th central FD stencil
        stencil = (/1.0_rk/12.0_rk, -2.0_rk/3.0_rk, 0.0_rk, +2.0_rk/3.0_rk, -1.0_rk/12.0_rk/)
    endif

    if (params%order_discretization == "FD_2nd_central") then
        allocate(stencil(-1:+1))
        ! Tam & Webb, 4th order optimized (for first derivative)
        stencil = (/-0.5_rk, 0.0_rk, 0.5_rk/)
    endif

    if (params%order_discretization == "CDF44") then
        allocate( stencil(-6:6) )
        stencil = (/ -2.0d0**(-9.d0), 0.0d0,  9.0d0*2.0d0**(-8.d0), -2.0d0**(-5.d0),  -63.0d0*2.0d0**(-9.d0),  9.0d0*2.0d0**(-5.d0), &
        87.0d0*2.0d0**(-7.d0), &
        9.0d0*2.0d0**(-5.d0), -63.0d0*2.0d0**(-9.d0), -2.0d0**(-5.d0), 9.0d0*2.0d0**(-8.d0), 0.0d0, -2.0d0**(-9.d0)/) ! H TILDE

        write(*,'(13(es15.8,1x))') stencil
    endif


    open(17, file=trim(adjustl(file))//'.info.txt', status='replace')
    write(17,'(A,1x,A,1x,A," g=",i1," Bs=",i2, " coarseWins=",L1)') trim(params%order_discretization), &
    trim(params%order_predictor), dir, params%n_ghosts, params%Bs(1), params%ghost_nodes_redundant_point_coarseWins
    close(17)

    !---------------------------------------------------------------------------

    if ((params%order_discretization == "FD_4th_central_optimized").and.(params%n_ghosts<4)) then
        call abort(33,"no enough g")
    endif


    g = params%n_ghosts
    a = lbound(stencil, dim=1)
    b = ubound(stencil, dim=1)
    iz = 1

    call allocate_grid(params, lgt_block, hvy_block, hvy_neighbor, &
    lgt_active, hvy_active, lgt_sortednumlist, hvy_tmp=hvy_tmp)

    call init_ghost_nodes( params )

    call read_mesh(file, params, lgt_n, hvy_n, lgt_block)

    ! call create_equidistant_grid( params, lgt_block, hvy_neighbor, lgt_active, &
    ! lgt_n, lgt_sortednumlist, hvy_active, hvy_n, 1, .true., 1 )

    call create_active_and_sorted_lists(params, lgt_block, lgt_active, &
    lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_ID=1)

    call update_neighbors(params, lgt_block, hvy_neighbor, lgt_active, &
    lgt_n, lgt_sortednumlist, hvy_active, hvy_n)

    dx_fine = (2.0_rk**-max_active_level(lgt_block, lgt_active, lgt_n))*domain(2)/real((Bs(2)-1), kind=rk)
    nx_fine = nint(domain(2)/dx_fine)

    write(*,*) "nx_fine=", nx_fine

    !---------------------------------------------------------------------------
    ! save the grid (for plotting in python)
    !---------------------------------------------------------------------------
    open(19, file=trim(adjustl(file))//'.operator_grid_points.txt', status='replace')
    do iblock = 1, hvy_n
        call hvy2lgt(lgt_id, hvy_active(iblock), params%rank, params%number_blocks)
        call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )
        do ix = g+1, Bs(1)+g
            do iy = g+1, Bs(2)+g
                x = dble(ix-(g+1)) * dx(1) + x0(1)
                y = dble(iy-(g+1)) * dx(2) + x0(2)

                ixx = nint(x/dx_fine)+1
                iyy = nint(y/dx_fine)+1

                write(19,*) ixx, iyy, x, y
            enddo
        enddo
    enddo
    close(19)


    !---------------------------------------------------------------------------
    ! compute operator matrix
    !---------------------------------------------------------------------------
    open(17, file=trim(adjustl(file))//'.operator_matrix.txt', status='replace')
    open(18, file=trim(adjustl(file))//'.operator_matrix_nosync.txt', status='replace')

    do iblock = 1, hvy_n
        do ix = g+1, Bs(1)+g
            do iy = g+1, Bs(2)+g
                write(*,*) "--------------------point---------------------------"
                !---------------------------------------------------------------
                ! reset entire grid to zeros (do not care about performance, just reset all)
                hvy_block = 0.0_rk

                !---------------------------------------------------------------
                ! set this one point we're looking at to 1
                call hvy2lgt(lgt_id, hvy_active(iblock), params%rank, params%number_blocks)
                call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

                x = dble(ix-(g+1)) * dx(1) + x0(1)
                y = dble(iy-(g+1)) * dx(2) + x0(2)

                if (abs((x-domain(1))) <=1.0e-9) x = 0.0_rk
                if (abs((y-domain(2))) <=1.0e-9) y = 0.0_rk

                ! save its indices on the fine grid
                ixx = nint(x/dx_fine)+1
                iyy = nint(y/dx_fine)+1

                hvy_block(ix, iy, iz, 1, hvy_active(iblock)) = 1.0_rk

                !---------------------------------------------------------------
                ! synchronize ghosts (important! if e.g. coarseWins is active and you happen to set the redundant value of a refined block, its overwritten to be zero again)
                ! Note: this also applies to coarse block bordering on a coarse block, if its ID is lower.
                ! In fact, each point is then computed only once. Note: if you set the point on a high lgt_id, then
                ! it will be "sync'ed down" to lower light IDs, so you can find the point more than once
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                call sync_ghosts(params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n)
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                params_acm%filter_type = "wavelet_filter"
                params_acm%order_predictor = params%order_predictor

                do k = 1, hvy_n
                    call hvy2lgt(lgt_id, hvy_active(k), params%rank, params%number_blocks)
                    call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

                    call filter_ACM( 0.0_rk, hvy_block(:,:,:,:,hvy_active(k)), g, x0, dx, &
                    hvy_block(:,:,:,:,hvy_active(k)), hvy_block(:,:,:,:,hvy_active(k)) )
                end do


                !---------------------------------------------------------------
                ! Now compute the derivative
                do k = 1, hvy_n
                    call hvy2lgt(lgt_id, hvy_active(k), params%rank, params%number_blocks)
                    call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

                    ! if (lgt_block(lgt_id, params%max_treelevel+IDX_MESH_LVL) == 1) then
                    !     ! coarse
                    !     dx_inv = 0.5_rk
                    ! else
                    !     ! fine
                    !     dx_inv = 1.0_rk !/ dx(1)
                    ! endif
                    dx_inv = 1.0_rk / dx(1)

                    do iy2 = g+1, Bs(2)+g
                        do ix2 = g+1, Bs(1)+g

                            if (dir=="x") then
                                u_dx = sum( stencil*hvy_block(ix2+a:ix2+b, iy2, iz, 1, hvy_active(k)) )*dx_inv
                            elseif (dir=='y') then
                                u_dx = sum( stencil*hvy_block(ix2, iy2+a:iy2+b, iz, 1, hvy_active(k)) )*dx_inv
                            else
                                call abort(123,'X or Y baby, nothing else.')
                            endif

                            hvy_block(ix2,iy2,iz,2,hvy_active(k)) = u_dx
                        end do
                    end do
                end do


                ! This second sync step also synchronizes the derivative we computed previously
                ! Note that on a coarse/fine interface, wabbit computes two values for the derivative
                ! on the coarse and fine level. This synchronizing step lets us keep only either of those,
                ! depending on fineWins or coarseWins
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                call sync_ghosts(params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n)
                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






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
                                write(17,'(i6,1x,i6,1x,es15.8)') ixx+(iyy-1)*nx_fine, ixx2+(iyy2-1)*nx_fine, val
                                write(*,*) "python col=", ixx2+(iyy2-1)*nx_fine -1 , "row=", ixx+(iyy-1)*nx_fine -1, val, "xy=",x,y, "block", k
                            endif
                        end do
                    end do
                end do



            enddo
        enddo
    enddo
    close(17)
    close(18)

    !---------------------------------------------------------------------------
    ! Now, we set a sine wave on the grid and compute its derivative.
    ! We store it in 1D indexing format.
    ! It can be compared directly with the above operator:
    ! (python) error = np.max(np.abs( u_dx - np.matmul(D.transpose(),u) ))
    ! With the right ghost nodes sync'ing (also after the derivative) the process
    ! works like a charm, illustrating that it IS indeed the right operator matrix.
    !---------------------------------------------------------------------------


    ! !---------------------------------------------------------------
    ! ! reset entire grid to zeros (do not care about performance, just reset all)
    ! hvy_block = 0.0_rk
    !
    ! !---------------------------------------------------------------
    ! ! set sine-wave
    ! do k = 1, hvy_n
    !     call hvy2lgt(lgt_id, hvy_active(k), params%rank, params%number_blocks)
    !     call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )
    !
    !     do iy2 = g+1, Bs(2)+g
    !         do ix2 = g+1, Bs(1)+g
    !             x2 = dble(ix2-(g+1)) * dx(1) + x0(1)
    !             y2 = dble(iy2-(g+1)) * dx(2) + x0(2)
    !
    !             if (abs((x2-domain(1))) <=1.0e-9) x2 = 0.0_rk ! not needed strictly speaking
    !             if (abs((y2-domain(2))) <=1.0e-9) y2 = 0.0_rk
    !
    !             hvy_block(ix2, iy2, iz, 1, hvy_active(k)) = sin(2.0_rk*pi*x2)
    !         enddo
    !     enddo
    ! enddo
    !
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! call sync_ghosts(params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n)
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !
    ! !---------------------------------------------------------------
    ! ! Now compute the derivative
    ! do k = 1, hvy_n
    !     call hvy2lgt(lgt_id, hvy_active(k), params%rank, params%number_blocks)
    !     call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )
    !
    !     dx_inv = 1.0_rk/dx(1)
    !
    !     do iy2 = g+1, Bs(2)+g
    !         do ix2 = g+1, Bs(1)+g
    !             if (params%order_discretization == "FD_2nd_central") then
    !                 if (dir=="x") then
    !                     u_dx = (hvy_block(ix2+1,iy2,iz,1,hvy_active(k)) - hvy_block(ix2-1,iy2,iz,1,hvy_active(k)))*dx_inv*0.5_rk
    !                 elseif (dir=='y') then
    !                     u_dx = (hvy_block(ix2,iy2+1,iz,1,hvy_active(k)) - hvy_block(ix2,iy2-1,iz,1,hvy_active(k)))*dx_inv*0.5_rk
    !                 else
    !                     call abort(123,'XY')
    !                 endif
    !
    !
    !             elseif (params%order_discretization == 'FD_4th_central_optimized')then
    !                 if (dir=="x") then
    !                     u_dx = (  a(-3)*hvy_block(ix2-3, iy2, iz, 1, hvy_active(k)) &
    !                     + a(-2)*hvy_block(ix2-2, iy2, iz, 1, hvy_active(k)) &
    !                     + a(-1)*hvy_block(ix2-1, iy2, iz, 1, hvy_active(k)) &
    !                     + a(0 )*hvy_block(ix2  , iy2, iz, 1, hvy_active(k)) &
    !                     + a(+1)*hvy_block(ix2+1, iy2, iz, 1, hvy_active(k)) &
    !                     + a(+2)*hvy_block(ix2+2, iy2, iz, 1, hvy_active(k)) &
    !                     + a(+3)*hvy_block(ix2+3, iy2, iz, 1, hvy_active(k)))*dx_inv
    !                 elseif (dir=='y') then
    !                     u_dx = (  a(-3)*hvy_block(ix, iy2-3, iz, 1, hvy_active(k)) &
    !                     + a(-2)*hvy_block(ix, iy2-2, iz, 1, hvy_active(k)) &
    !                     + a(-1)*hvy_block(ix, iy2-1, iz, 1, hvy_active(k)) &
    !                     + a(0 )*hvy_block(ix, iy2  , iz, 1, hvy_active(k)) &
    !                     + a(+1)*hvy_block(ix, iy2+1, iz, 1, hvy_active(k)) &
    !                     + a(+2)*hvy_block(ix, iy2+2, iz, 1, hvy_active(k)) &
    !                     + a(+3)*hvy_block(ix, iy2+3, iz, 1, hvy_active(k)))*dx_inv
    !                 else
    !                     call abort(123,'XY')
    !                 endif
    !             else
    !                 call abort(771272,"unknown")
    !             endif
    !
    !             hvy_block(ix2,iy2,iz,2,hvy_active(k)) = u_dx
    !         end do
    !     end do
    ! end do
    !
    !
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! call sync_ghosts(params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n)
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    !
    ! open(20, file=trim(adjustl(file))//'.u_derivative.txt', status='replace')
    ! open(21, file=trim(adjustl(file))//'.u_function.txt', status='replace')
    ! do k = 1, hvy_n
    !     call hvy2lgt(lgt_id, hvy_active(k), params%rank, params%number_blocks)
    !     call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )
    !
    !     do iy2 = g+1, Bs(2)+g
    !         do ix2 = g+1, Bs(1)+g
    !             x = dble(ix2-(g+1)) * dx(1) + x0(1)
    !             y = dble(iy2-(g+1)) * dx(2) + x0(2)
    !
    !             if (abs((x-domain(1))) <=1.0e-9) x = 0.0_rk
    !             if (abs((y-domain(2))) <=1.0e-9) y = 0.0_rk
    !
    !             ixx2 = nint(x/dx_fine)+1
    !             iyy2 = nint(y/dx_fine)+1
    !
    !             write(20,'(i6,1x,es15.8)') ixx2+(iyy2-1)*nx_fine, hvy_block(ix2, iy2, iz, 2, hvy_active(k))
    !             write(21,'(i6,1x,es15.8)') ixx2+(iyy2-1)*nx_fine, hvy_block(ix2, iy2, iz, 1, hvy_active(k))
    !         end do
    !     end do
    ! end do
    ! close(20)
    ! close(21)
    !
    ! call write_field("ww_000.h5", time, 1, 1, params, lgt_block, hvy_block, lgt_active, lgt_n, hvy_n, hvy_active )
    ! call write_field("wwdx_000.h5", time, 1, 2, params, lgt_block, hvy_block, lgt_active, lgt_n, hvy_n, hvy_active )

end subroutine
