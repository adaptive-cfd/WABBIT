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

    real(kind=rk) :: t0

    t0 = MPI_wtime()

    call synchronize_ghosts_generic_sequence( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

    call toc( "WRAPPER: sync ghosts", MPI_wtime()-t0 )
end subroutine sync_ghosts
