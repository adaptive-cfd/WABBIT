!> \brief save data main function, call write data routine
!! input:    - time loop parameter
!!           - parameter array
!!           - light data array
!!           - heavy data array
! ********************************************************************************************

subroutine save_data(iteration, time, params, hvy_block, hvy_tmp, hvy_mask, tree_ID)

    implicit none
    !> time loop parameters
    real(kind=rk), intent(in)                       :: time
    integer(kind=ik), intent(in)                    :: iteration
    !> user defined parameter structure
    type (type_params), intent(in)                  :: params
    integer(kind=ik), intent(in)                    :: tree_ID
    !> heavy data array - block data
    real(kind=rk), intent(inout)                    :: hvy_block(:, :, :, :, :)
    !> heavy temp data: used for saving, filtering, and helper qtys (reaction rate, mask function)
    real(kind=rk), intent(inout)                    :: hvy_tmp(:, :, :, :, :)
    ! mask data. we can use different trees (4est module) to generate time-dependent/indenpedent
    ! mask functions separately. This makes the mask routines tree-level routines (and no longer
    ! block level) so the physics modules have to provide an interface to create the mask at a tree
    ! level. All parts of the mask shall be included: chi, boundary values, sponges.
    ! On input, the mask array is correctly filled. You cannot create the full mask here.
    real(kind=rk), intent(inout)                    :: hvy_mask(:, :, :, :, :)

    ! loop variable
    integer(kind=ik)      :: k, lgt_id, hvy_id
    ! file name
    character(len=cshort) :: fname, tmp
    ! cpu time variables for running time calculation
    real(kind=rk)         :: t0, x0(1:3), dx(1:3)
    integer(kind=2)       :: n_domain(1:3)

    t0 = MPI_Wtime()
    if (params%rank == 0) then
        write(*,'("IO: Saving data triggered, time=",g15.8)')  time
    endif

    n_domain = 0

    ! create mask function at current time. (this routine is rarely called and thus
    ! the overhead of calling createMask_tree if the mask is not stored is supposed
    ! to be small)
    call createMask_tree(params, time, hvy_mask, hvy_tmp)


    ! preparatory step. The physics modules have to copy everything they want to
    ! save to disk to the work array. missing qty's shall be computed.
    do k = 1, hvy_n(tree_ID_flow)

        hvy_id = hvy_active(k,tree_ID_flow)

        ! convert given hvy_id to lgt_id for block spacing routine
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

        ! get block spacing for RHS
        call get_block_spacing_origin( params, lgt_id, x0, dx )

        if ( .not. All(params%periodic_BC) ) then
            ! check if block is adjacent to a boundary of the domain, if this is the case we use one sided stencils
            call get_adjacent_boundary_surface_normal( lgt_block(lgt_id, 1:lgt_block(lgt_id,params%max_treelevel+IDX_MESH_LVL)), &
            params%domain_size, params%Bs, params%dim, n_domain )
        endif

        ! call preparatory routines. this routine saves the variables to be stored
        ! to disk in the work array. This way, we can also store derived variables
        ! such as the vorticity. Note in most cases, this copies just the state vector
        ! to work.
        call PREPARE_SAVE_DATA_meta(params%physics_type, time, hvy_block(:,:,:,:,hvy_id), &
        params%n_ghosts, x0, dx, hvy_tmp(:,:,:,:,hvy_id), hvy_mask(:,:,:,:,hvy_id), n_domain)

    enddo

    ! actual saving step. one file per component.
    ! loop over components/qty's:
    do k = 1, params%N_fields_saved
        ! physics modules shall provide an interface for wabbit to know how to label
        ! the components to be stored to hard disk (in the work array)
        call FIELD_NAMES_meta(params%physics_type, k, tmp)

        ! create filename
        if (params%use_iteration_as_fileid) then
          write( fname,'(a, "_", i12.12, ".h5")') trim(adjustl(tmp)), iteration
        else
          write( fname,'(a, "_", i12.12, ".h5")') trim(adjustl(tmp)), nint(time * 1.0e6_rk)
        endif

        ! actual writing
        call saveHDF5_tree( fname, time, iteration, k, params, hvy_tmp, tree_ID)
    enddo

    call toc( "save_data", MPI_wtime()-t0 )
end subroutine save_data
