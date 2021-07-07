!> \brief functional wrapper for the 2D and 3D version of update neighbors.
!! input:     - light data array
!!            - params struct
!! output:    - neighbor list array
! ********************************************************************************************
subroutine update_neighbors(params, lgt_block, hvy_neighbor, lgt_active, lgt_n, &
    lgt_sortednumlist, hvy_active, hvy_n, skip_diagonal_neighbors)

    implicit none
    type (type_params), intent(in)      :: params                   !> user defined parameter structure
    integer(kind=ik), intent(in)        :: lgt_block(:, :)          !> light data array
    integer(kind=ik), intent(out)       :: hvy_neighbor(:,:)        !> heavy data array - neighbor data
    integer(kind=ik), intent(in)        :: lgt_active(:)            !> list of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_n                    !> number of active blocks (light data)
    integer(kind=tsize), intent(in)     :: lgt_sortednumlist(:,:)   !> sorted list of numerical treecodes, used for block finding
    integer(kind=ik), intent(in)        :: hvy_active(:)            !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n                    !> number of active blocks (heavy data)
    logical, intent(in), optional       :: skip_diagonal_neighbors  ! currently not working (Thomas, 02-2021) [unused]
    logical                             :: error = .false., error2 = .false.
    integer(kind=ik )                   :: mpierror, k
    real(kind=rk)                       :: x0(1:3), dx(1:3)

    if ( params%dim == 3 ) then
        if (present(skip_diagonal_neighbors)) then
            call update_neighbors_3D(params, lgt_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n, error, skip_diagonal_neighbors)
        else
            call update_neighbors_3D(params, lgt_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n, error)
        endif
    else
        if (present(skip_diagonal_neighbors)) then
            call update_neighbors_2D(params, lgt_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n, error, skip_diagonal_neighbors)
        else
            call update_neighbors_2D(params, lgt_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n, error)
        endif
    end if

    if (error) then
        ! open(14, file="lgt_block.txt", status='replace')
        ! do k = 1, lgt_n
        !     write(14,'(20(i5,1x))') lgt_block(lgt_active(k),:)
        ! enddo
        ! close(14)
        !
        ! open(14, file="hvy_neighbors.txt", status='replace')
        ! do k = 1, lgt_n
        !     write(74,'(90(i5,1x))') hvy_neighbor(lgt_active(k),:)
        ! enddo
        ! close(14)
        !
        ! open(14, file="x0_dx.txt", status='replace')
        ! do k = 1, lgt_n
        !     write(*,*) "dumping", k, lgt_active(k)
        !     call get_block_spacing_origin( params,lgt_active(k) , lgt_block, x0, dx )
        !     write(14,'(6(es15.6,1x))') x0, dx
        ! enddo
        ! close(14)

        call abort(71737, "Grid error: we did not find a neighboring block. Either your data is corrupt (reading from file) or you found a bug.")
    endif

    ! Is there any non-periodic boundary ?
    if ( .not. All(params%periodic_BC) ) then
        call remove_nonperiodic_neighbors(params, lgt_block, hvy_neighbor, hvy_active, hvy_n)
    endif

end subroutine update_neighbors
