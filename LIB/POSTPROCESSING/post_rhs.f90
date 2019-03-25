subroutine post_rhs(params)
    use module_precision
    use module_mesh
    use module_params
    use module_IO
    use module_mpi
    use module_operators
    use module_physics_metamodule
    use module_time_step

    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=80)      :: fname_ini, fname_input, dummy
    real(kind=rk)          :: time
    integer(kind=ik)       :: iteration, k, lgt_id, lgt_n, hvy_n, tc_length, g
    integer(kind=ik), dimension(3) :: Bs
    character(len=2)       :: order

    integer(kind=ik), allocatable      :: lgt_block(:, :)
    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_work(:, :, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_tmp(:, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_gridQ(:, :, :, :, :)
    integer(kind=ik), allocatable      :: hvy_neighbor(:,:)
    integer(kind=ik), allocatable      :: lgt_active(:), hvy_active(:)
    integer(kind=tsize), allocatable   :: lgt_sortednumlist(:,:)
    character(len=80)                  :: fname
    real(kind=rk), dimension(3)        :: dx, x0
    real(kind=rk), allocatable         :: us(:,:,:,:)
    integer(hid_t)                     :: file_id
    real(kind=rk), dimension(3)        :: domain
    logical                            :: adaptive
    integer                            :: i

    !-----------------------------------------------------------------------------------------------------
    ! get values from command line (filename and level for interpolation)
    call get_command_argument(2, fname_ini)
    ! does the user need help?
    if (fname_ini=='--help' .or. fname_ini=='--h') then
        if (params%rank==0) then
            write(*,*) "------------------------------------------------------------------"
            write(*,*) "wabbit postprocessing routine for computation of RHS"
            write(*,*) "------------------------------------------------------------------"
            write(*,*) " ./wabbit-post --compute-rhs PARAMS.ini state1_0000.h5 [state2_0000.h5]"
            write(*,*) "    Give all statevector components as separate files (inicond section"
            write(*,*) "    in paramsfile is ignored)"
            write(*,*) " Output is written to "
            write(*,*) " rhs1_0000.h5, rhs2...... timestamp matches input timestamp"
            write(*,*) "------------------------------------------------------------------"
        end if
        return
    endif

    ! read ini-file and save parameters in struct
    call ini_file_to_params( params, fname_ini )
    ! have the pysics module read their own parameters
    call init_physics_modules( params, fname_ini, params%n_gridQ  )


    ! get some parameters from one of the files (they should be the same in all of them)
    call get_command_argument(3, fname_input)
    call read_attributes(fname_input, lgt_n, time, iteration, domain, Bs, tc_length, params%dim)


    adaptive = .false.
    do i = 1, command_argument_count()
        call get_command_argument(i, dummy)
        if (index(dummy,'--adaptive') /= 0) adaptive = .true.
    enddo

    ! in usual parameter files, RK4 (or some other RK) is used an requires a lot of memory
    ! here we do not need that, and hence pretent to use a basic scheme (EE1 maybe)
    deallocate(params%butcher_tableau)
    allocate(params%butcher_tableau(1,1))

    ! only (4* , for safety) lgt_n/number_procs blocks necessary (since we do not want to refine)
    !> \todo change that for 3d case
    params%number_blocks = 4*lgt_n/params%number_procs
    Bs = params%Bs
    g = params%n_ghosts

    ! allocate data
    call allocate_grid(params, lgt_block, hvy_block, hvy_neighbor, &
    lgt_active, hvy_active, lgt_sortednumlist, hvy_work, hvy_tmp, hvy_gridQ)

    ! read mesh and field
    call read_mesh(fname_input, params, lgt_n, hvy_n, lgt_block)


    ! read in all components of statevector
    do k = 1, params%n_eqn
        call get_command_argument(2+k, fname_input)
        call read_field(fname_input, k, params, hvy_block, hvy_n)
    enddo

    ! create lists of active blocks (light and heavy data)
    ! update list of sorted nunmerical treecodes, used for finding blocks
    call create_active_and_sorted_lists( params, lgt_block, lgt_active, &
    lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true. )

    ! update neighbor relations
    call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, &
    lgt_n, lgt_sortednumlist, hvy_active, hvy_n )



    if ( adaptive ) then
        call adapt_mesh( time, params, lgt_block, hvy_block, hvy_neighbor, lgt_active, &
        lgt_n, lgt_sortednumlist, hvy_active, hvy_n, params%coarsening_indicator, hvy_tmp, hvy_gridQ )
    endif

    ! compute right hand side
    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )
    call RHS_wrapper(time, params, hvy_block, hvy_work(:,:,:,:,:,1), hvy_gridQ, lgt_block, &
    hvy_active, hvy_n )

    ! if (params%filter_type /= "no_filter") then
    !     call sync_ghosts( params, lgt_block, hvy_work(:,:,:,:,:,1), hvy_neighbor, hvy_active, hvy_n )
    !     call filter_wrapper(time, params, hvy_work(:,:,:,:,:,1), hvy_tmp, lgt_block, hvy_active, hvy_n)
    ! end if


    ! save result to disk
    do k = 1, params%n_eqn
        write( fname,'("rhs",i1,"_", i12.12, ".h5")') k, nint(time * 1.0e6_rk)

        call write_field(fname, time, iteration, k, params, lgt_block,&
        hvy_work(:,:,:,:,:,1), lgt_active, lgt_n, hvy_n, hvy_active )
    enddo

end subroutine post_rhs
