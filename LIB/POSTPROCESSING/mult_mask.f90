subroutine mult_mask(params)
    use module_globals
    use module_mesh
    use module_params
    use module_mpi
    use module_operators
    use module_physics_metamodule
    use module_time_step
    use module_forestMetaData

    implicit none

    !> parameter struct
    type (type_params), intent(inout)  :: params
    character(len=cshort)      :: fname_input, fname_mask, fname_result, operation
    real(kind=rk)          :: time
    integer(kind=ik)       :: iteration, k, lgt_id, tc_length
    integer(kind=ik), dimension(3) :: Bs
    character(len=2)       :: order

    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :), hvy_work(:, :, :, :, :, :)
    integer(kind=ik)                   :: tree_ID=1, hvy_id

    character(len=cshort)              :: fname
    real(kind=rk), dimension(3)        :: dx, x0
    real(kind=rk), allocatable :: us(:,:,:,:)
    integer(hid_t)                     :: file_id
    real(kind=rk), dimension(3)        :: domain
    ! NOTE: after 24/08/2022, the arrays lgt_active/lgt_n hvy_active/hvy_n as well as lgt_sortednumlist,
    ! hvy_neighbors, tree_N and lgt_block are global variables included via the module_forestMetaData. This is not
    ! the ideal solution, as it is trickier to see what does in/out of a routine. But it drastically shortenes
    ! the subroutine calls, and it is easier to include new variables (without having to pass them through from main
    ! to the last subroutine.)  -Thomas


    ! this routine works only on one tree
    allocate( hvy_n(1), lgt_n(1) )

    !-----------------------------------------------------------------------------------------------------
    ! operation name
    call get_command_argument(1, operation)

    ! does the user need help?
    call get_command_argument(2, fname_input)

    if (fname_input=='--help' .or. fname_input=='--h') then
        if (params%rank==0) then
            write(*,*) "------------------------------------------------------------------"
            write(*,*) "wabbit postprocessing routine for mutliplication with mask"
            write(*,*) "------------------------------------------------------------------"
            write(*,*) " ./wabbit-post --mult-mask input_0000.h5 mask_0000.h5 result_0000.h5"
            write(*,*) "------------------------------------------------------------------"
            write(*,*) " ./wabbit-post --mult-mask-direct input_0000.h5 mask_0000.h5 result_0000.h5"
            write(*,*) "------------------------------------------------------------------"
            write(*,*) " ./wabbit-post --mult-mask-inverse input_0000.h5 mask_0000.h5 result_0000.h5"
            write(*,*) "------------------------------------------------------------------"
        end if
        return
    endif

    call check_file_exists(trim(fname_input))

    call get_command_argument(3, fname_mask)
    call check_file_exists(trim(fname_mask))
    call get_command_argument(4, fname_result)

    ! get some parameters from one of the files (they should be the same in all of them)
    call read_attributes(fname_input, lgt_n(tree_ID), time, iteration, domain, Bs, tc_length, &
    params%dim, periodic_BC=params%periodic_BC, symmetry_BC=params%symmetry_BC)

    if (params%rank==0) then
        if (operation == "--mult-mask" .or. operation == "--mult-mask-inverse") then
            write(*,*) "------------------------------------------------------------------"
            write(*,*) "Mutliplying field with 1-mask: result = (1.0-mask)*input"
            write(*,*) "------------------------------------------------------------------"
            write(*,*) "input= ", trim(adjustl(fname_input))
            write(*,*) "mask=  ", trim(adjustl(fname_mask))
            write(*,*) "result=", trim(adjustl(fname_result))
            write(*,*) "------------------------------------------------------------------"
        elseif (operation == "--mult-mask-direct") then
            write(*,*) "------------------------------------------------------------------"
            write(*,*) "Mutliplying field with mask: result = mask*input"
            write(*,*) "------------------------------------------------------------------"
            write(*,*) "input= ", trim(adjustl(fname_input))
            write(*,*) "mask=  ", trim(adjustl(fname_mask))
            write(*,*) "result=", trim(adjustl(fname_result))
            write(*,*) "------------------------------------------------------------------"
        end if
    endif

    params%Jmax = tc_length
    params%n_eqn = 2
    params%domain_size(1) = domain(1)
    params%domain_size(2) = domain(2)
    params%domain_size(3) = domain(3)
    params%Bs = Bs
    params%g = 2_ik
    params%order_predictor = "multiresolution_2nd"

    allocate(params%butcher_tableau(1,1))

    ! only (4* , for safety) lgt_n/number_procs blocks necessary
    !> \todo change that for 3d case
    params%number_blocks = 4_ik*lgt_n(tree_ID)/params%number_procs

    ! allocate data
    call allocate_forest(params, hvy_block, hvy_work=hvy_work)

    ! read data
    call readHDF5vct_tree( (/fname_input, fname_mask/), params, hvy_block, tree_ID)

    ! create lists of active blocks (light and heavy data)
    ! update list of sorted nunmerical treecodes, used for finding blocks
    call updateMetadata_tree( params, tree_ID)

    call sync_ghosts_all( params, lgt_block, hvy_block, hvy_neighbor, hvy_active(:,tree_ID), hvy_n(tree_ID) )

    ! calculate product
    if (operation == "--mult-mask" .or. operation == "--mult-mask-inverse") then
        do k = 1, hvy_n(tree_ID)
            hvy_block(:,:,:,1,hvy_active(k,tree_ID)) = hvy_block(:,:,:,1,hvy_active(k,tree_ID)) * (1.0_rk - hvy_block(:,:,:,2,hvy_active(k,tree_ID)))
        end do
    elseif (operation == "--mult-mask-direct") then
        do k = 1, hvy_n(tree_ID)
            hvy_block(:,:,:,1,hvy_active(k,tree_ID)) = hvy_block(:,:,:,1,hvy_active(k,tree_ID)) * hvy_block(:,:,:,2,hvy_active(k,tree_ID))
        end do
    end if

    call saveHDF5_tree(fname_result, time, iteration, 1, params, hvy_block, tree_ID)

end subroutine mult_mask
