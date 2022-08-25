subroutine post_generate_forest(params)
    use module_precision
    use module_mesh
    use module_params
    use module_mpi
    use module_operators
    use module_physics_metamodule
    use module_time_step
    use module_helpers
    use module_forestMetaData

    implicit none

    type (type_params), intent(inout)  :: params
    character(len=cshort) ::  fname_out,mode

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_work(:, :, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_tmp(:, :, :, :, :)

    integer :: hvy_id, lgt_id, fsize, j, tree_ID, N_treeGenerate
    integer(kind=ik) :: it,ix,iy,iz,k,p,n,m, Bs(1:3)=17, g=4,i,f
    real(kind=rk) :: time, domain(1:3), norm, x,y,z,xrel,yrel,zrel,x0(1:3),dx(1:3),dt, r
    integer(kind=ik) :: freq(15**2) = (/(i, i=1,15**2)/)
    !-----------------------------------------------------------------------------------------------------
    ! get values from command line (filename and level for interpolation)
    call get_command_argument(2, mode)

    ! does the user need help?
    if (mode=='--help' .or. mode=='--h' .or. mode=='-h') then
        if (params%rank==0) then
            write(*,*) "------------------------------------------------------------------"
            write(*,*) "./wabbit-post --generate_forest --Jmax=[5] --Ntrees=[20] --dim=[2,3]"
            write(*,*) "------------------------------------------------------------------"
            write(*,*) ""
            write(*,*) ""
            write(*,*) ""
            write(*,*) ""
            write(*,*) "------------------------------------------------------------------"
        end if
        return
    endif

    call get_cmd_arg( "--Jmax", params%max_treelevel, default=5_ik )
    call get_cmd_arg( "--Ntrees", tree_N, 20 )
    call get_cmd_arg( "--dim", params%dim, 2 )

    params%number_blocks = tree_N*2**(params%dim*params%max_treelevel)  ! just to get some memory:
    params%domain_size = (/ 30, 30 ,30 /)
    params%Bs = Bs
    params%min_treelevel = 1
    params%n_eqn = 1
    params%n_ghosts = g
    params%forest_size = tree_N+2
    fsize = params%forest_size
    params%order_predictor = "multiresolution_4th"
    params%block_distribution = "sfc_hilbert"
    params%time_step_method = 'none'

    N_treeGenerate = tree_N ! note tree_N is overwritten elsewhere

    N_MAX_COMPONENTS = params%n_eqn ! used for ghost node sync'ing (buffer allocation)


    ! we have to allocate grid if this routine is called for the first time
    call allocate_forest(params, hvy_block, hvy_work, hvy_tmp=hvy_tmp)

    call reset_forest(params)
    ! The ghost nodes will call their own setup on the first call, but for cleaner output
    ! we can also just do it now.
    call init_ghost_nodes( params )

    hvy_neighbor = -1
    lgt_n = 0 ! reset number of active light blocks
    hvy_n = 0

    dt = 2*pi/real(N_treeGenerate+1,kind=rk)

    ! generate random shuffle of frequencies
    do i = 15**2, 1, -1
        r=rand_nbr()
        j = modulo(i**2, 15**2)+1
        k = freq(i)
        freq(i) = freq(j)
        freq(j) = k
    end do

    do tree_ID = 1, N_treeGenerate
        call createEquidistantGrid_tree( params, params%max_treelevel, .false., tree_ID )

        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k, tree_ID)

            ! convert given hvy_id to lgt_id for block spacing routine
            call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
            ! get block spacing for RHS
            call get_block_spacing_origin( params, lgt_id, x0, dx )

            if (params%dim == 2) then
                do ix = g+1,Bs(1)+g
                    do iy = g+1,Bs(2)+g
                        hvy_block(ix,iy,1,1,hvy_id) = 0.0_rk
                        ! compute x,y coordinates from spacing and origin
                        x = dble(ix-(g+1)) * dx(1) + x0(1)
                        y = dble(iy-(g+1)) * dx(2) + x0(2)
                        do p = 1, N_treeGenerate
                            f = freq(p)
                            m = modulo(f,15)
                            n = (f-m)/15
                            xrel = x - ( 2* m - 1 )
                            yrel = y - ( 2* n - 1 )
                            ! set actual inicond gauss blob
                            hvy_block(ix,iy,1,1,hvy_id) =hvy_block(ix,iy,1,1,hvy_id)+ &
                            exp(-f/100.0_rk)*sin(pi*f*dt*tree_ID)* bump(sqrt(xrel**2 + yrel**2))
                        end do
                    end do
                end do
            else
                do ix = 1,Bs(1)+2*g
                    do iy = 1,Bs(2)+2*g
                        do iz = 1,Bs(3)+2*g
                            hvy_block(ix,iy,iz,1,hvy_id) = 0.0_rk
                            ! compute x,y coordinates from spacing and origin
                            x = dble(ix-(g+1)) * dx(1) + x0(1)
                            y = dble(iy-(g+1)) * dx(2) + x0(2)
                            z = dble(iz-(g+1)) * dx(3) + x0(3)
                            do p = 1, N_treeGenerate
                                f = freq(p)
                                m = modulo(f,15)
                                n = (f-m)/15
                                xrel = x - ( 2* m - 1 )
                                yrel = y - ( 2* n - 1 )
                                zrel = z - ( 2* n - 1 )
                                ! set actual inicond gauss blob
                                hvy_block(ix,iy,iz,1,hvy_id) =hvy_block(ix,iy,iz,1,hvy_id)+ &
                                exp(-f/3.0_rk)*sin(pi*f*dt*tree_ID)* bump(sqrt(xrel**2 + yrel**2 + zrel**2))
                            end do
                        end do
                    end do
                end do
            end if

        end do

        write( fname_out, '("u_", i12.12, ".h5")')  tree_ID

        call saveHDF5_tree(fname_out, real(tree_ID*dt,kind=rk), 0_ik, 1, params, hvy_block, tree_ID)
    end do

    call deallocate_forest(params, hvy_block, hvy_work, hvy_tmp )

    ! make a summary of the program parts, which have been profiled using toc(...)
    ! and print it to stdout
    call summarize_profiling( WABBIT_COMM )


contains

    function bump(x) result (res)
        use module_precision

        implicit none
        real(kind=rk), intent (in) :: x
        real(kind=rk) :: res

        if (abs(x)<1) then
            res = exp(-1/(1-x**2))
        else
            res = 0.0_rk
        end if
    end function bump

end subroutine
