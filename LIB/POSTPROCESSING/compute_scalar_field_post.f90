!> \brief postprocessing routine for differential operators on a scalar field fld saved in .h5 files
!-----------------------------------------------------------------------------------------------------

subroutine compute_scalar_field_post(params)
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_operators
    use module_forestMetaData

    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=cshort)              :: file_fld, operator
    real(kind=rk)                      :: time
    integer(kind=ik)                   :: iteration, k, lgt_id, tc_length
    integer(kind=ik), dimension(3)     :: Bs
    character(len=2)                   :: order

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_work(:, :, :, :, :, :), hvy_tmp(:, :, :, :, :)
    integer(kind=ik)                   :: tree_ID=1, hvy_id

    character(len=cshort)              :: fname
    real(kind=rk), dimension(3)        :: dx, x0
    integer(hid_t)                     :: file_id
    real(kind=rk), dimension(3)        :: domain

    ! this routine works only on one tree
    allocate( hvy_n(1), lgt_n(1) )

    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas


    !-----------------------------------------------------------------------------------------------------
    ! get values from command line (filename and level for interpolation)
    call get_command_argument(1, operator)
    call get_command_argument(2, file_fld)

    ! does the user need help?
    if (file_fld=='--help' .or. file_fld=='--h') then
        if (params%rank==0) then
            write(*,*) "-----------------------------------------------------------"
            write(*,*) " Wabbit postprocessing: gradient"
            write(*,*) "-----------------------------------------------------------"
            write(*,*) " Computes differential operator of a scalar field. Output is stored"
            write(*,*) " in predefined files."
            write(*,*) "-----------------------------------------------------------"
            write(*,*) " --gradient"
            write(*,*) "./wabbit-post --gradient source.h5 [ORDER]"
            write(*,*) " Computes divergence of a scalar field, saves in "
            write(*,*) " gradx_*.h5 grady_*.h5 [gradz_*.h5]"
            write(*,*) "-----------------------------------------------------------"
        end if
        return
    endif

    call check_file_exists(trim(file_fld))

    ! get some parameters from one of the files (they should be the same in all of them)
    call read_attributes(file_fld, lgt_n(tree_ID), time, iteration, domain, Bs, tc_length, params%dim, &
    periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)

    call get_command_argument(3, order)

    ! decide which order for discretization and predictor is used. Note predictor
    ! is used in ghost nodes sync'ing
    if (order == "4") then
        params%order_discretization = "FD_4th_central_optimized"
        params%order_predictor = "multiresolution_4th"
        params%g = 4_ik

    elseif (order == "2") then
        params%order_discretization = "FD_2nd_central"
        params%order_predictor = "multiresolution_2nd"
        params%g = 2_ik

    else
        call abort(8765,"chosen discretization order invalid or not (yet) implemented. choose between 4 (FD_4th_central_optimized) and 2 (FD_2nd_central)")

    end if

    params%Jmax = tc_length
    params%n_eqn = params%dim
    params%domain_size(1) = domain(1)
    params%domain_size(2) = domain(2)
    params%domain_size(3) = domain(3)
    params%Bs = Bs
    allocate(params%butcher_tableau(1,1))
    allocate(params%symmetry_vector_component(1))
    params%symmetry_vector_component(1) = "0" ! scalar
    ! no refinement is made in this postprocessing tool; we therefore allocate about
    ! the number of blocks in the file (and not much more than that)
    params%number_blocks = ceiling(  real(lgt_n(tree_ID))/real(params%number_procs) )

    ! allocate data
    call allocate_forest(params, hvy_block, hvy_tmp=hvy_tmp)

    ! read mesh and field
    call readHDF5vct_tree( (/file_fld/), params, hvy_block, tree_ID)

    ! create lists of active blocks (light and heavy data)
    ! update list of sorted nunmerical treecodes, used for finding blocks
    call createActiveSortedLists_tree( params, tree_ID)

    ! update neighbor relations
    call updateNeighbors_tree(params, tree_ID, search_overlapping=.false.)

    call sync_ghosts_tree( params, hvy_block, tree_ID )

    ! calculate vorticity from velocities
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k,tree_ID)

        call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
        call get_block_spacing_origin( params, lgt_id, x0, dx )

        if (operator == "--gradient") then
            if (params%dim == 3) then
                call compute_gradient( hvy_block(:,:,:,1,hvy_id), dx, params%Bs, params%g, &
                params%order_discretization, hvy_tmp(:,:,:,1:3,hvy_id))
            else
                call compute_gradient( hvy_block(:,:,:,1,hvy_id), dx, params%Bs, params%g, &
                params%order_discretization, hvy_tmp(:,:,:,1:2,hvy_id))
            end if

        else
            call abort(1812017,"operator is not --gradient")
        endif
    end do

    if (operator=="--gradient") then
        write( fname,'(a, "_", i12.12, ".h5")') 'gradx', nint(time * 1.0e6_rk)
        call saveHDF5_tree(fname, time, iteration, 1, params, hvy_tmp, tree_ID )

        write( fname,'(a, "_", i12.12, ".h5")') 'grady', nint(time * 1.0e6_rk)
        call saveHDF5_tree(fname, time, iteration, 2, params, hvy_tmp, tree_ID )

        if (params%dim == 3) then
            write( fname,'(a, "_", i12.12, ".h5")') 'gradz', nint(time * 1.0e6_rk)
            call saveHDF5_tree(fname, time, iteration, 3, params, hvy_tmp, tree_ID )
        end if
    endif
end subroutine compute_scalar_field_post
