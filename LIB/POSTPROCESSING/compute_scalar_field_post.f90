    !> \file
! WABBIT
!> \name compute_scalar_field_post.f90
!> \version 0.1
!> \author dk
!
!> \brief postprocessing routine for differential operators on a scalar field fld saved in .h5 files
! = log ======================================================================================
!
!-----------------------------------------------------------------------------------------------------

subroutine compute_scalar_field_post(params)
    use module_precision
    use module_mesh
    use module_params
    use module_IO
    use module_mpi
    use module_operators

    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=80)      :: file_fld, operator
    real(kind=rk)          :: time
    integer(kind=ik)       :: iteration, k, lgt_id, lgt_n, hvy_n, tc_length
    integer(kind=ik), dimension(3) :: Bs
    character(len=2)       :: order

    integer(kind=ik), allocatable      :: lgt_block(:, :)
    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_work(:, :, :, :, :, :), hvy_tmp(:, :, :, :, :)
    integer(kind=ik), allocatable      :: hvy_neighbor(:,:)
    integer(kind=ik), allocatable      :: lgt_active(:), hvy_active(:)
    integer(kind=tsize), allocatable   :: lgt_sortednumlist(:,:)
    character(len=80)                  :: fname
    real(kind=rk), dimension(3)        :: dx, x0
    integer(hid_t)                     :: file_id
    real(kind=rk), dimension(3)        :: domain

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
    call read_attributes(file_fld, lgt_n, time, iteration, domain, Bs, tc_length, params%dim)

    call get_command_argument(3, order)

    ! decide which order for discretization and predictor is used. Note predictor
    ! is used in ghost nodes sync'ing
    if (order == "4") then
        params%order_discretization = "FD_4th_central_optimized"
        params%order_predictor = "multiresolution_4th"
        params%n_ghosts = 4_ik

    elseif (order == "2") then
        params%order_discretization = "FD_2nd_central"
        params%order_predictor = "multiresolution_2nd"
        params%n_ghosts = 2_ik

    else
        call abort(8765,"chosen discretization order invalid or not (yet) implemented. choose between 4 (FD_4th_central_optimized) and 2 (FD_2nd_central)")

    end if

    params%max_treelevel = tc_length
    params%n_eqn = params%dim
    params%domain_size(1) = domain(1)
    params%domain_size(2) = domain(2)
    params%domain_size(3) = domain(3)
    params%Bs = Bs
    allocate(params%butcher_tableau(1,1))
    ! no refinement is made in this postprocessing tool; we therefore allocate about
    ! the number of blocks in the file (and not much more than that)
    params%number_blocks = ceiling(  real(lgt_n)/real(params%number_procs) )

    ! allocate data
    call allocate_grid(params, lgt_block, hvy_block, hvy_neighbor, &
    lgt_active, hvy_active, lgt_sortednumlist, hvy_tmp=hvy_tmp)

    ! read mesh and field
    call read_mesh(file_fld, params, lgt_n, hvy_n, lgt_block)
    call read_field(file_fld, 1, params, hvy_block, hvy_n)

    ! create lists of active blocks (light and heavy data)
    ! update list of sorted nunmerical treecodes, used for finding blocks
    call create_active_and_sorted_lists( params, lgt_block, lgt_active, &
    lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_ID=1)
    ! update neighbor relations
    call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, &
    lgt_n, lgt_sortednumlist, hvy_active, hvy_n )

    call sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

    ! calculate vorticity from velocities
    do k = 1, hvy_n
        call hvy_id_to_lgt_id(lgt_id, hvy_active(k), params%rank, params%number_blocks)
        call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

        if (operator == "--gradient") then
            if (params%dim == 3) then
                call gradient( hvy_block(:,:,:,1,hvy_active(k)), &
                dx, params%Bs, params%n_ghosts, &
                params%order_discretization, hvy_tmp(:,:,:,1:3,hvy_active(k)))
            else
                call gradient( hvy_block(:,:,:,1,hvy_active(k)), &
                dx, params%Bs, params%n_ghosts, &
                params%order_discretization, hvy_tmp(:,:,:,1:2,hvy_active(k)))
            end if

        else
            call abort(1812017,"operator is not --gradient")
        endif
    end do

    if (operator=="--gradient") then
        write( fname,'(a, "_", i12.12, ".h5")') 'gradx', nint(time * 1.0e6_rk)
        call write_field(fname, time, iteration, 1, params, lgt_block,&
        hvy_tmp, lgt_active, lgt_n, hvy_n, hvy_active )
        write( fname,'(a, "_", i12.12, ".h5")') 'grady', nint(time * 1.0e6_rk)
        call write_field(fname, time, iteration, 2, params, lgt_block,&
        hvy_tmp, lgt_active, lgt_n, hvy_n, hvy_active )
        if (params%dim == 3) then
            write( fname,'(a, "_", i12.12, ".h5")') 'gradz', nint(time * 1.0e6_rk)
            call write_field(fname, time, iteration, 3, params, lgt_block, &
            hvy_tmp, lgt_active, lgt_n, hvy_n, hvy_active )
        end if
    endif
end subroutine compute_scalar_field_post
