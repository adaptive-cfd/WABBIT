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

    ! call reset_ghost_nodes(  params, hvy_block, hvy_active, hvy_n )


    call synchronize_ghosts_generic_sequence( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

    ! do k = 1, hvy_n
    !     hvy_id = hvy_active(k)
    !     call hvy_id_to_lgt_id(lgt_id, hvy_id, params%rank, params%number_blocks)
    !
    !     call get_adjacent_boundary_surface_normal( lgt_block(lgt_id, 1:lgt_block(lgt_id,params%max_treelevel+IDX_MESH_LVL)), &
    !     params%domain_size, params%Bs, params%dim, n_domain )
    !
    !
    !     ! if (((maxval(hvy_block(:,:,:,1, hvy_id))>1.0e3).or.(minval(hvy_block(:,:,:,1, hvy_id))<-1.0e3)).and.(maxval(abs(n_domain))>0)) then
    !     if (((maxval(hvy_block(:,:,:,1, hvy_id))>1.0e3).or.(minval(hvy_block(:,:,:,1, hvy_id))<-1.0e3))) then
    !         write(*,'(74(i4,1x))') hvy_neighbor(hvy_id,:)
    !         write(*,*) "n", n_domain
    !         write(*,'(25(es12.4,1x))') hvy_block(:,1,1,1, hvy_id)
    !         write(*,*) "max", maxloc(hvy_block(:,:,:,1, hvy_id)), maxval(hvy_block(:,:,:,1, hvy_id)), hvy_id
    !         write(*,*) "min", minloc(hvy_block(:,:,:,1, hvy_id)), minval(hvy_block(:,:,:,1, hvy_id))
    !         call abort(12345, " maybe some ghost node patches are not filled")
    !     endif
    ! enddo

    call toc( "WRAPPER: sync ghosts", MPI_wtime()-t0 )
end subroutine sync_ghosts
