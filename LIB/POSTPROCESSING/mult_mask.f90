subroutine mult_mask(params)
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
    character(len=80)      :: fname_input, fname_mask, fname_result
    real(kind=rk)          :: time
    integer(kind=ik)       :: iteration, k, lgt_id, lgt_n, hvy_n, Bs, tc_length, g
    character(len=2)       :: order

    integer(kind=ik), allocatable      :: lgt_block(:, :)
    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_work(:, :, :, :, :, :)
    integer(kind=ik), allocatable      :: hvy_neighbor(:,:)
    integer(kind=ik), allocatable      :: lgt_active(:), hvy_active(:)
    integer(kind=tsize), allocatable   :: lgt_sortednumlist(:,:)
    character(len=80)                  :: fname
    real(kind=rk), dimension(3)        :: dx, x0
    real(kind=rk), allocatable :: us(:,:,:,:)
    integer(hid_t)                     :: file_id
    real(kind=rk), dimension(3)        :: domain

    !-----------------------------------------------------------------------------------------------------
    ! does the user need help?
    call get_command_argument(2, fname_input)

    if (fname_input=='--help' .or. fname_input=='--h') then
        if (params%rank==0) then
            write(*,*) "------------------------------------------------------------------"
            write(*,*) "wabbit postprocessing routine for mutliplication with mask"
            write(*,*) "------------------------------------------------------------------"
            write(*,*) " ./wabbit-post --mult-mask input_0000.h5 mask_0000.h5 result_0000.h5"
            write(*,*) "------------------------------------------------------------------"
        end if
        return
    endif

    call check_file_exists(trim(fname_input))

    call get_command_argument(3, fname_mask)
    call check_file_exists(trim(fname_mask))
    call get_command_argument(4, fname_result)

    ! get some parameters from one of the files (they should be the same in all of them)
    call read_attributes(fname_input, lgt_n, time, iteration, domain, Bs, tc_length, params%dim)

    if (params%rank==0) then
        write(*,*) "------------------------------------------------------------------"
        write(*,*) "Mutliplying field with 1-mask: result = (1.0-mask)*input"
        write(*,*) "------------------------------------------------------------------"
        write(*,*) "input= ", trim(adjustl(fname_input))
        write(*,*) "mask=  ", trim(adjustl(fname_mask))
        write(*,*) "result=", trim(adjustl(fname_result))
        write(*,*) "------------------------------------------------------------------"
    endif

    if (params%dim==3) then
        params%threeD_case = .true.
    else
        params%threeD_case = .false.
    end if

    params%max_treelevel = tc_length
    params%n_eqn = 2
    params%domain_size(1) = domain(1)
    params%domain_size(2) = domain(2)
    params%domain_size(3) = domain(3)
    params%Bs = Bs

    allocate(params%butcher_tableau(1,1))

    ! only (4* , for safety) lgt_n/number_procs blocks necessary
    !> \todo change that for 3d case
    params%number_blocks = 4_ik*lgt_n/params%number_procs
    if (params%rank==0) params%number_blocks = params%number_blocks + mod(lgt_n, params%number_procs)

    ! allocate data
    call allocate_grid(params, lgt_block, hvy_block, hvy_neighbor, &
    lgt_active, hvy_active, lgt_sortednumlist, hvy_work=hvy_work)

    ! read mesh and field
    call read_mesh(fname_input, params, lgt_n, hvy_n, lgt_block)
    call read_field(fname_input, 1, params, hvy_block, hvy_n)
    call read_field(fname_mask , 2, params, hvy_block, hvy_n)


    ! create lists of active blocks (light and heavy data)
    ! update list of sorted nunmerical treecodes, used for finding blocks
    call create_active_and_sorted_lists( params, lgt_block, lgt_active, &
    lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true. )

    ! update neighbor relations
    call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, &
    lgt_n, lgt_sortednumlist, hvy_active, hvy_n )


    ! calculate vorticity from velocities
    do k = 1, hvy_n
        hvy_block(:,:,:,1,hvy_active(k)) = hvy_block(:,:,:,1,hvy_active(k)) * (1.0_rk - hvy_block(:,:,:,2,hvy_active(k)))
    end do

    call write_field(fname_result, time, iteration, 1, params, lgt_block, &
    hvy_block, lgt_active, lgt_n, hvy_n, hvy_active )

end subroutine mult_mask
