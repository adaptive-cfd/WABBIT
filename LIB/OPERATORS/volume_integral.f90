!> \brief computation of the velocity volume integral
!! input:    - heavy block data
!!           - params
!!           - hvy_active, lgt_block list
!! output:   - volume integral
! ********************************************************************************************

! subroutine volume_integral(volume_int, hvy_block, params, hvy_active, hvy_n, lgt_block)
!
!     implicit none
!     real(kind=rk), dimension(3), intent(out) :: volume_int                !> global volume integral
!     real(kind=rk), intent(in)                :: hvy_block(:, :, :, :, :)  !> actual block data
!     type (type_params), intent(in)           :: params                    !> physics parameter structure
!     integer(kind=ik), intent(in)             :: hvy_active(:)             !> list of active blocks (heavy data)
!     integer(kind=ik), intent(in)             :: lgt_block(:,:)
!     integer(kind=ik), intent(in)             :: hvy_n
!
!     real(kind=rk), dimension(3)              :: dx, x0                    !> origin and spacing of the block
!     real(kind=rk), dimension(3)              :: int_block, int_local      !> volume integral of one block and of all blocks on my process
!     integer(kind=ik)                         :: Bs, g                     !> grid parameter
!     integer(kind=ik)                         :: k, lgt_id                 !> loop variables
!     integer                                  :: mpi_err
!
!
!    int_block  = 0.0_rk
!    int_local  = 0.0_rk
!    volume_int = 0.0_rk
!
!    Bs = params%Bs
!    g  = params%n_ghosts
!
!
!    do k = 1, hvy_n
!
!        ! convert given hvy_id to lgt_id for block spacing routine
!        call hvy2lgt( lgt_id, hvy_active(k), params%rank, params%number_blocks )
!
!        ! get block spacing for RHS
!        call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )
!
!        if (params%threeD_case) then
!            int_block(1) = sum(hvy_block(g+1:Bs+g-1,g+1:Bs+g-1,g+1:Bs+g-1,1,hvy_active(k)))&
!                *dx(1)*dx(2)*dx(3)
!            int_block(2) = sum(hvy_block(g+1:Bs+g-1,g+1:Bs+g-1,g+1:Bs+g-1,2,hvy_active(k)))&
!                *dx(1)*dx(2)*dx(3)
!            int_block(3) = sum(hvy_block(g+1:Bs+g-1,g+1:Bs+g-1,g+1:Bs+g-1,3,hvy_active(k)))&
!                *dx(1)*dx(2)*dx(3)
!        else
!            int_block(1) = sum(hvy_block(g+1:Bs+g-1,g+1:Bs+g-1,1,1,hvy_active(k)))&
!                *dx(1)*dx(2)
!            int_block(2) = sum(hvy_block(g+1:Bs+g-1,g+1:Bs+g-1,1,2,hvy_active(k)))&
!                *dx(1)*dx(2)
!        end if
!
!        int_local = int_local + int_block
!
!     end do
!
!
!    call MPI_ALLREDUCE(int_local, volume_int, 3, MPI_DOUBLE_PRECISION, MPI_SUM, &
!        params%WABBIT_COMM, mpi_err)
!
! end subroutine volume_integral
