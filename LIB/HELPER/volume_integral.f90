!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name volume_integral.f90
!> \version 0.5
!> \author sm
!
!> \brief computation of the velocity volume integral 
!
!>
!! input:    -  \n
!! output:   -  \n
!!
!!
!! = log ======================================================================================
!! \n
!! 27/07/17 - create
! ********************************************************************************************

subroutine volume_integral(volume_int, hvy_block, params, hvy_active, hvy_n, lgt_block)

!----------------------------------------------------------------------------------------------
! modules

!----------------------------------------------------------------------------------------------
! variables

    implicit none

    !> global volume integral
    real(kind=rk), dimension(3), intent(out) :: volume_int
    !> actual block data
    real(kind=rk), intent(in)                :: hvy_block(:, :, :, :, :)
    !> physics parameter structure
    type (type_params), intent(in)           :: params
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)             :: hvy_active(:)
    integer(kind=ik), intent(in)             :: lgt_block(:,:)
    integer(kind=ik), intent(in)             :: hvy_n


    !> origin and spacing of the block
    real(kind=rk), dimension(3)              :: dx, x0
    !> volume integral of one block and of all blocks on my process
    real(kind=rk), dimension(3)              :: int_block, int_local
    !> grid parameter
    integer(kind=ik)                         :: Bs, g
    !> loop variables
    integer(kind=ik)                         :: k, lgt_id
    !> MPI error
    logical                                  :: mpi_err

   
!---------------------------------------------------------------------------------------------
! variables initialization

   int_block  = 0.0_rk
   int_local  = 0.0_rk
   volume_int = 0.0_rk

   Bs = params%number_block_nodes
   g  = params%number_ghost_nodes

!---------------------------------------------------------------------------------------------
! main body

   do k = 1, hvy_n

       ! convert given hvy_id to lgt_id for block spacing routine
       call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )

       ! get block spacing for RHS
       call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )
       
       if (params%threeD_case) then 
           int_block(1) = sum(hvy_block(g+1:Bs+g,g+1:Bs+g,g+1:Bs+g,1,hvy_active(k)))*dx(1)*dx(2)*dx(3)
           int_block(2) = sum(hvy_block(g+1:Bs+g,g+1:Bs+g,g+1:Bs+g,2,hvy_active(k)))*dx(1)*dx(2)*dx(3)
           int_block(3) = sum(hvy_block(g+1:Bs+g,g+1:Bs+g,g+1:Bs+g,3,hvy_active(k)))*dx(1)*dx(2)*dx(3)
       else
           int_block(1) = sum(hvy_block(g+1:Bs+g,g+1:Bs+g,1,1,hvy_active(k)))*dx(1)*dx(2)
           int_block(2) = sum(hvy_block(g+1:Bs+g,g+1:Bs+g,1,2,hvy_active(k)))*dx(1)*dx(2)
       end if

       int_local = int_local + int_block
   
    end do


   call MPI_ALLREDUCE(int_local, volume_int, 3, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, mpi_err)


   
end subroutine volume_integral
