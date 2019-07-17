
!> \file
!> \name sparse_to_dense.f90
!> \author PKrah
!> \brief postprocessing for sparsing data from a dense wabbit field
!> \date 29.03.2019 creation
!-----------------------------------------------------------------------------------------------------

subroutine dense_to_sparse(params)
    use module_precision
    use module_mesh
    use module_params
    use module_IO
    use module_mpi
    use module_initialization

    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=80)      :: indicator="threshold-state-vector", file_in, args
    character(len=80)      :: tail_string
    real(kind=rk)          :: time, eps=-1.0_rk,maxmem=-1.0_rk
    integer(kind=ik)       :: iteration
    character(len=80), allocatable :: file_out(:)
    integer(kind=ik), allocatable           :: lgt_block(:, :)
    real(kind=rk), allocatable              :: hvy_block(:, :, :, :, :), hvy_work(:, :, :, :, :, :)
    real(kind=rk), allocatable              :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik), allocatable           :: hvy_neighbor(:,:)
    integer(kind=ik), allocatable           :: lgt_active(:,:), hvy_active(:,:), hvy_n(:), lgt_n(:)
    integer(kind=tsize), allocatable        :: lgt_sortednumlist(:,:,:)
    integer(kind=ik)                        :: max_neighbors, level, k, tc_length, lgt_n_tmp
    integer(kind=ik), dimension(3)          :: Bs
    integer(hid_t)                          :: file_id
    character(len=2)                        :: level_in, order
    real(kind=rk), dimension(3)             :: domain
    integer(hsize_t), dimension(2)          :: dims_treecode
    integer(kind=ik)                        :: treecode_size, number_dense_blocks, &
    n_opt_args, i, l, dim
    !-----------------------------------------------------------------------------------------------------

    call get_command_argument(2, file_in)
    if (file_in == '--help' .or. file_in == '--h') then
        if ( params%rank==0 ) then
            write(*,*) "postprocessing subroutine to refine/coarse mesh to a uniform grid (up and downsampling ensured). command line:"
            write(*,*) "mpi_command -n number_procs ./wabbit-post --dense-to-sparse source.h5 target.h5 target_treelevel order-predictor(2 or 4)"
        end if
        return
    end if
    !----------------------------------
    ! read predefined params
    !----------------------------------
    n_opt_args = 1 ! counting all extra arguments, which are not *h5 files

    do i = 1, command_argument_count()
        call get_command_argument(i,args)
        !-------------------------------
        ! order of predictor
        if ( index(args,"--order=")==1 ) then
            read(args(9:len_trim(args)),* ) order
            n_opt_args = n_opt_args + 1
        end if
        !-------------------------------
        ! Threshold indicator [threshold-vorticity, threshold-state-vector (default)]
        if ( index(args,"--indicator=")==1 ) then
            read(args(13:len_trim(args)),* ) indicator
            n_opt_args = n_opt_args + 1
        end if
        !-------------------------------
        ! Threshold for keeping wavelet coefficients
        if ( index(args,"--eps=")==1 ) then
            read(args(7:len_trim(args)),* ) eps
            n_opt_args = n_opt_args + 1
        end if
        !-------------------------------
        ! MEMORY AVAILABLE
        if ( index(args,"--memory=")==1 ) then
            read(args(10:len_trim(args)-2),* ) maxmem
            n_opt_args = n_opt_args + 1
        endif

    end do

    ! Check parameters for correct inputs:
    if (order == "4") then
        params%order_predictor = "multiresolution_4th"
        params%n_ghosts = 4_ik
    else
        params%order_predictor = "multiresolution_2nd"
        params%n_ghosts = 2_ik
    end if

    if (eps > 0) then
        params%eps=eps
    else
        call abort(2303191,"You must specify the threshold value --eps")
    endif
    params%coarsening_indicator=indicator
    params%forest_size = 1

    params%n_eqn = command_argument_count() - n_opt_args
    if (params%n_eqn <= 0 ) call abort(2603191,"No files are specified for sparsing!!!")

    allocate(params%input_files(params%n_eqn))
    allocate(params%field_names(params%n_eqn))
    allocate(file_out(params%n_eqn))
    allocate(params%threshold_state_vector_component(params%n_eqn))
    params%eps_normalized = .true.
    params%physics_type = "POD"
    params%threshold_state_vector_component = .true.

    !-------------------------------------------
    ! check and find common params in all h5-files
    !-------------------------------------------
    call get_command_argument(n_opt_args+1, params%input_files(1))
    call read_attributes(params%input_files(1), lgt_n_tmp, time, iteration, params%domain_size, &
    params%Bs,params%max_treelevel, params%dim)

    do i = 1, params%n_eqn
        call get_command_argument(n_opt_args+i, file_in)
        call check_file_exists(trim(file_in))
        call read_attributes(file_in, lgt_n_tmp, time, iteration, domain, Bs, level, dim)

        params%min_treelevel = 1
        params%max_treelevel = max(params%max_treelevel, level) ! find the maximal level of all snapshot

        if (any(params%Bs .ne. Bs)) call abort( 203192, " Block size is not consistent ")
        if (params%dim .ne. dim) call abort(243191,"Dimensions do not agree!")
        if ( abs(sum(params%domain_size(1:dim) - domain(1:dim))) > 1e-14 ) call abort( 203195, "Domain size is not consistent ")

        ! Concatenate "sparse" with filename
        params%input_files(i) = trim(file_in)
        file_out(i) = trim(file_in)
    end do



    ! in postprocessing, it is important to be sure that the parameter struct is correctly filled:
    ! most variables are unfortunately not automatically set to reasonable values. In simulations,
    ! the ini files parser takes care of that (by the passed default arguments). But in postprocessing
    ! we do not read an ini file, so defaults may not be set.
    allocate(params%butcher_tableau(1,1))


    params%block_distribution="sfc_hilbert"

    ! read attributes from file. This is especially important for the number of
    ! blocks the file contains: this will be the number of active blocks right
    ! after reading.
    if (params%dim==3) then
        ! how many blocks do we need for the desired level?
        number_dense_blocks = 8_ik**level
        max_neighbors = 74
    else
        number_dense_blocks = 4_ik**level
        max_neighbors = 12
    end if

    if (params%rank==0) then
        write(*,'(80("-"))')
        write(*,*) "Wabbit dense-to-sparse."
        do i = 1, params%n_eqn
            write(*,'(A20,1x,A80)') "Reading file:", params%input_files(i)
            write(*,'(A20,1x,A80)') "Writing to file:", file_out(i)
        end do
        write(*,'(A20,1x,A80)') "Predictor used:", params%order_predictor
        write(*,'(A20,1x,es9.3)') "eps:", params%eps
        write(*,'(A20,1x,A80)')"indicator:", params%coarsening_indicator
        write(*,'(80("-"))')
    endif

    ! is lgt_n > number_dense_blocks (downsampling)? if true, allocate lgt_n blocks
    !> \todo change that for 3d case
    params%number_blocks = ceiling( 4.0*dble(max(lgt_n_tmp, number_dense_blocks)) / dble(params%number_procs) )

    if (params%rank==0) then
        write(*,'("Data dimension: ",i1,"D")') params%dim
        write(*,'("File contains Nb=",i6," blocks of size Bs=",i4," x ",i4," x ",i4)') lgt_n_tmp, Bs(1),Bs(2),Bs(3)
        write(*,'("Domain size is ",3(g12.4,1x))') domain
        write(*,'("Time=",g12.4," it=",i9)') time, iteration
        write(*,'("Length of treecodes in file=",i3," in memory=",i3)') level, params%max_treelevel
        write(*,'("NCPU=",i6)') params%number_procs
        write(*,'("File   Nb=",i6," blocks")') lgt_n_tmp
        write(*,'("Memory Nb=",i6)') params%number_blocks
        write(*,'("Dense  Nb=",i6)') number_dense_blocks
    endif

    !----------------------------------
    ! allocate data and reset grid
    !----------------------------------
    call allocate_grid(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, hvy_active, &
        lgt_sortednumlist, hvy_work=hvy_work, hvy_tmp=hvy_tmp, hvy_n=hvy_n, lgt_n=lgt_n)

    ! reset the grid: all blocks are inactive and empty
    call reset_grid( params, lgt_block, hvy_block, hvy_work, hvy_tmp, hvy_neighbor, lgt_active(:,tree_ID_flow), &
    lgt_n(tree_ID_flow), hvy_active(:,tree_ID_flow), hvy_n(tree_ID_flow), lgt_sortednumlist(:,:,tree_ID_flow), .true. )

    ! The ghost nodes will call their own setup on the first call, but for cleaner output
    ! we can also just do it now.
    call init_ghost_nodes( params )

    !----------------------------------
    ! READ Grid and coarse if possible
    !----------------------------------
    params%adapt_mesh=.true.
    params%adapt_inicond=.true.
    params%read_from_files=.true.
    call set_initial_grid( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, hvy_active, &
    lgt_n, hvy_n, lgt_sortednumlist, params%adapt_inicond, time, iteration, hvy_tmp=hvy_tmp )

    !----------------------------------
    ! Write sparse files
    !----------------------------------
    do i = 1, params%n_eqn
        call write_field(file_out(i), time, iteration, i, params, lgt_block, &
        hvy_block, lgt_active(:,tree_ID_flow), lgt_n(tree_ID_flow), hvy_n(tree_ID_flow), &
        hvy_active(:,tree_ID_flow))
    enddo

    call deallocate_grid(params, lgt_block, hvy_block, hvy_neighbor, lgt_active,&
    hvy_active, lgt_sortednumlist, hvy_work, hvy_tmp=hvy_tmp, hvy_n=hvy_n , lgt_n=lgt_n)
end subroutine dense_to_sparse
