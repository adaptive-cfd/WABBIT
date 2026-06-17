subroutine save_data(iteration, time, params, hvy_block, hvy_tmp, hvy_mask, tree_ID)

    implicit none
    !> time loop parameters
    real(kind=rk), intent(in)                       :: time
    integer(kind=ik), intent(in)                    :: iteration
    !> user defined parameter structure
    type (type_params), intent(inout)               :: params
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
    character(len=clong) :: fname, tmp
    ! cpu time variables for running time calculation
    real(kind=rk)         :: t0, x0(1:3), dx(1:3)
    integer(kind=2)       :: n_domain(1:3)

    t0 = MPI_Wtime()
    if (params%rank == 0) then
        write(*,'("IO: Saving data triggered, time=",g15.8)')  time
    endif

    n_domain = 0

    ! we need to sync ghost nodes in order to compute the vorticity, if it is used and stored.
    call sync_ghosts_tree( params, hvy_block, tree_ID_flow )

    ! create mask function at current time. (this routine is rarely called and thus
    ! the overhead of calling createMask_tree if the mask is not stored is supposed
    ! to be small)
    call createMask_tree(params, time, hvy_mask, hvy_tmp)

    ! extra preparatory step. Some things have to be prepared on grid-level and not on block-level. This is especially important when deriving any quantities by solving the poisson-equation
    if (params%physics_type == "NSPP") then
        call pressure_from_velocity(params, time, hvy_block, hvy_tmp, hvy_mask, tree_ID_flow)
    endif


    ! preparatory step. The physics modules have to copy everything they want to
    ! save to disk to the work array. missing qty's shall be computed.
    do k = 1, hvy_n(tree_ID_flow)
        hvy_id = hvy_active(k,tree_ID_flow)
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
        ! get block spacing for RHS
        call get_block_spacing_origin( params, lgt_id, x0, dx )
        ! check if block is adjacent to a boundary of the domain, if this is the case we use one sided stencils
        if ( .not. All(params%periodic_BC) ) call get_adjacent_boundary_surface_normal( params, lgt_id, n_domain )

        ! call preparatory routines. this routine saves the variables to be stored
        ! to disk in the work array. This way, we can also store derived variables
        ! such as the vorticity. Note in most cases, this copies just the state vector
        ! to work.
        ! In case that hvy_mask is not used (which means all its sizes are one), we do some fortran magic to avoid errors
        call PREPARE_SAVE_DATA_meta(params%physics_type, time, hvy_block(:,:,:,:,hvy_id), &
        params%g, x0, dx, hvy_tmp(:,:,:,:,hvy_id), hvy_mask(:,:,:,:,merge(1, hvy_id, size(hvy_mask,5) == 1)), n_domain)

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
            write( fname,'(a, "_", a, ".h5")') trim(adjustl(tmp)), trim(adjustl(timestr(time)))
        endif

        ! actual writing
        call saveHDF5_tree( fname, time, iteration, k, params, hvy_tmp, tree_ID, no_sync=.false.)
    enddo

    call toc( "TOPLEVEL: save_data", 16, MPI_wtime()-t0 )
end subroutine save_data


subroutine save_backup( iteration, time, params, hvy_block, hvy_tmp, hvy_mask, tree_ID )

implicit none
    !> time loop parameters
    real(kind=rk), intent(in)                       :: time
    integer(kind=ik), intent(in)                    :: iteration
    !> user defined parameter structure
    type (type_params), intent(inout)               :: params
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
    character(len=clong) :: fname, tmp
    ! cpu time variables for running time calculation
    real(kind=rk)         :: t0, x0(1:3), dx(1:3), obsolete_time, obsolete_iteration
    integer(kind=2)       :: n_domain(1:3)
    ! variables for file deletion
    logical :: exists
    integer :: iu, ios

    t0 = MPI_Wtime()
    if (params%rank == 0) then
        write(*,'("IO: Saving data triggered, time=",g15.8)')  time
    endif

    obsolete_time = params%write_backup_time_list(1)
    do k = 1, params%write_backup_n-1
        params%write_backup_time_list(k) = params%write_backup_time_list(k+1)
    end do
    params%write_backup_time_list(params%write_backup_n) = time
    obsolete_iteration = params%write_backup_iteration_list(1)
    do k = 1, params%write_backup_n-1
        params%write_backup_iteration_list(k) = params%write_backup_iteration_list(k+1)
    end do
    params%write_backup_iteration_list(params%write_backup_n) = iteration

    n_domain = 0

    ! we need to sync ghost nodes in order to compute the vorticity, if it is used and stored.
    call sync_ghosts_tree( params, hvy_block, tree_ID_flow )

    ! create mask function at current time. (this routine is rarely called and thus
    ! the overhead of calling createMask_tree if the mask is not stored is supposed
    ! to be small)
    call createMask_tree(params, time, hvy_mask, hvy_tmp)

    ! extra preparatory step. Some things have to be prepared on grid-level and not on block-level. This is especially important when deriving any quantities by solving the poisson-equation
    if (params%physics_type == "NSPP") then
        call pressure_from_velocity(params, time, hvy_block, hvy_tmp, hvy_mask, tree_ID_flow)
    endif


    ! preparatory step. The physics modules have to copy everything they want to
    ! save to disk to the work array. missing qty's shall be computed.
    do k = 1, hvy_n(tree_ID_flow)
        hvy_id = hvy_active(k,tree_ID_flow)
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
        ! get block spacing for RHS
        call get_block_spacing_origin( params, lgt_id, x0, dx )
        ! check if block is adjacent to a boundary of the domain, if this is the case we use one sided stencils
        if ( .not. All(params%periodic_BC) ) call get_adjacent_boundary_surface_normal( params, lgt_id, n_domain )

        ! call preparatory routines. this routine saves the variables to be stored
        ! to disk in the work array. This way, we can also store derived variables
        ! such as the vorticity. Note in most cases, this copies just the state vector
        ! to work.
        ! In case that hvy_mask is not used (which means all its sizes are one), we do some fortran magic to avoid errors
        call PREPARE_SAVE_DATA_meta(params%physics_type, time, hvy_block(:,:,:,:,hvy_id), &
        params%g, x0, dx, hvy_tmp(:,:,:,:,hvy_id), hvy_mask(:,:,:,:,merge(1, hvy_id, size(hvy_mask,5) == 1)), n_domain)

    enddo

    ! actual saving step. one file per component.
    ! loop over components/qty's:
    do k = 1, params%N_fields_saved
        ! physics modules shall provide an interface for wabbit to know how to label
        ! the components to be stored to hard disk (in the work array)
        call FIELD_NAMES_meta(params%physics_type, k, tmp)

        ! create filename
        if (params%use_iteration_as_fileid) then
            write( fname,'( "backup-", a, "_", i12.12, ".h5")') trim(adjustl(tmp)), iteration
        else
            write( fname,'( "backup-", a, "_", a, ".h5")') trim(adjustl(tmp)), trim(adjustl(timestr(time)))
        endif

        ! actual writing
        call saveHDF5_tree( fname, time, iteration, k, params, hvy_tmp, tree_ID, no_sync=.false.)
    enddo

    ! now we need to delete the obsolete backups, only done on rank 0
    if (params%rank == 0 .and. (obsolete_time >= 0.0_rk .or. obsolete_iteration >= 0_ik)) then
        do k = 1, params%N_fields_saved
            call FIELD_NAMES_meta(params%physics_type, k, tmp)
            if (params%use_iteration_as_fileid) then
                write( fname,'( "backup-", a, "_", i12.12, ".h5")') trim(adjustl(tmp)), obsolete_iteration
            else
                write( fname,'( "backup-", a, "_", a, ".h5")') trim(adjustl(tmp)), trim(adjustl(timestr(obsolete_time)))
            endif
            ! file is deleted here
            inquire(file=fname, exist=exists)
            if (exists) then
                open(newunit=iu, file=fname, status='old', action='read', iostat=ios)
                if (ios == 0) then
                    close(iu, status='delete', iostat=ios)
                    write(*,'("IO: Deleted old backup file = ", a)') trim(adjustl(fname))
                end if
            end if
        end do
    endif

    call toc( "TOPLEVEL: save_data", 16, MPI_wtime()-t0 )

end subroutine save_backup



!> Delete all current backups
subroutine delete_backups( params )

    implicit none
    type (type_params), intent(inout)               :: params
    integer :: k, k_backups, iu, ios
    character(len=clong) :: fname, tmp
    logical :: exists

    do k_backups = 1, params%write_backup_n
        if (params%rank == 0 .and. (params%write_backup_time_list(k_backups) >= 0.0 .or. params%write_backup_iteration_list(k_backups) >= 0_ik)) then
            do k = 1, params%N_fields_saved
                call FIELD_NAMES_meta(params%physics_type, k, tmp)
                if (params%use_iteration_as_fileid) then
                    write( fname,'( "backup-", a, "_", i12.12, ".h5")') trim(adjustl(tmp)), params%write_backup_iteration_list(k_backups)
                else
                    write( fname,'( "backup-", a, "_", i6.6, i6.6, ".h5")') trim(adjustl(tmp)), int(params%write_backup_time_list(k_backups)+1.0e-12_rk, kind=ik), nint(max((params%write_backup_time_list(k_backups)-int(params%write_backup_time_list(k_backups)+1.0e-12_rk, kind=ik))*1.0e6_rk, 0.0_rk), kind=ik)
                endif
                inquire(file=fname, exist=exists)
                if (exists) then
                    open(newunit=iu, file=fname, status='old', action='read', iostat=ios)
                    if (ios == 0) then
                        close(iu, status='delete', iostat=ios)
                        write(*,'("IO: Deleted old backup file = ", a)') trim(adjustl(fname))
                    end if
                end if
            end do
        end if
    enddo

end subroutine delete_backups



!> Determine, if it is time to save data and/or backup. This is determined by the parameters in the params structure, which can be based on time, iteration, or walltime. The decision is made on rank 0 and then broadcasted to all other ranks to ensure synchronization.
subroutine is_it_time_to_save_data( params, time, iteration, tstart, it_is_time_to_save_data, it_is_time_to_save_backup )

    implicit none
    type (type_params), intent(inout)       :: params
    real(kind=rk), intent(in)               :: time
    integer(kind=ik), intent(in)            :: iteration
    real(kind=rk), intent(in)               :: tstart  !< cpu time at the start of the simulation, used for walltime-based saving
    logical, intent(out)                    :: it_is_time_to_save_data
    logical, intent(out)                    :: it_is_time_to_save_backup

    integer(kind=ik) :: mpicode

    ! ************
    ! determine if it is time to save data
    ! ************
    it_is_time_to_save_data = .false.

    if ((params%write_method=='fixed_freq' .and. modulo(iteration, params%write_freq)==0) .or. &
        (params%write_method=='fixed_time' .and. (abs(mod(time, params%write_time))<1.0e-12_rk .or. abs(mod(time, params%write_time)-params%write_time)<1.0e-12_rk))) then
        it_is_time_to_save_data = .true.
    endif

    ! do not save any output before this time (so maybe revoke the previous decision)
    if (time+1e-12_rk<params%write_time_first) then
        it_is_time_to_save_data = .false.
    endif
    ! save after walltime unit is not affected by write_time_first
    if ((MPI_wtime()-tstart) - params%walltime_last_write > params%walltime_write*3600.0_rk) then
        params%walltime_last_write = MPI_wtime()-tstart
        it_is_time_to_save_data = .true.
    endif
    ! it can rarely happen that not all proc arrive at the same time at the above condition, then some decide to
    ! save data and others do not. this is a rare but severe problem, to solve it, synchronize:
    call MPI_BCAST( it_is_time_to_save_data, 1, MPI_LOGICAL, 0, WABBIT_COMM, mpicode )

    ! ************
    ! determine if it is time to save backup
    ! ************
    it_is_time_to_save_backup = .false.

    if (params%write_backup_n <= 0) return  ! if no backup is wanted, we can skip the rest of the checks

    ! backup if specific time is reached
    if (abs(mod(time, params%write_backup_time))<1.0e-12_rk .or. abs(mod(time, params%write_backup_time)-params%write_backup_time)<1.0e-12_rk) then
        it_is_time_to_save_backup = .true.
    endif
    ! save backup after walltime
    if ((MPI_wtime()-tstart) - params%write_backup_walltime > params%write_backup_walltime*3600.0_rk) then
        params%write_backup_walltime = MPI_wtime()-tstart
        it_is_time_to_save_backup = .true.
    endif
    ! it can rarely happen that not all proc arrive at the same time at the above condition, then some decide to
    ! save data and others do not. this is a rare but severe problem, to solve it, synchronize:
    call MPI_BCAST( it_is_time_to_save_backup, 1, MPI_LOGICAL, 0, WABBIT_COMM, mpicode )

end subroutine is_it_time_to_save_data