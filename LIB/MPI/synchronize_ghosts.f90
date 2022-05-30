subroutine sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )
    implicit none
    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n
    integer :: k, hvy_id, lgt_id
    integer(kind=2) :: n_domain(1:3)
    real(kind=rk) :: t0

    t0 = MPI_wtime()

    if (REDUNDANT_GRID) then

        ! call reset_ghost_nodes(  params, hvy_block, hvy_active, hvy_n )
        if (params%wavelet_transform_type == "harten-multiresolution") then
            filter = .false.
            call sync_ghosts_redundantGrid( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )
        else
            filter = .false.
            call sync_ghosts_redundantGrid( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

            if (params%iter_ghosts) then
                do k = 1, NiterationsGhosts
                    filter = .true.
                    call sync_ghosts_redundantGrid( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )
                enddo
            endif
        endif

    else
        call sync_ghosts_uniqueGrid( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

    endif

    call toc( "WRAPPER: sync ghosts", MPI_wtime()-t0 )
end subroutine sync_ghosts
