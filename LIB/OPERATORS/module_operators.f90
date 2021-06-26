!> \brief module for all operator routines
! *********************************************************************************************

module module_operators

use mpi
use module_params     ! global parameters
use module_timing
use module_treelib

implicit none

PRIVATE
!**********************************************************************************************
! These are the important routines that are visible to WABBIT:
!**********************************************************************************************
PUBLIC :: compute_vorticity, compute_vorticity_abs, divergence, gradient, compute_Qcriterion, component_wise_tree_norm

contains

#include "compute_Qcriterion.f90"
#include "compute_vorticity.f90"
#include "divergence.f90"
#include "gradient.f90"


subroutine component_wise_tree_norm(params, lgt_block, hvy_block, hvy_active, hvy_n, which_norm, norm)
    implicit none

    type (type_params), intent(in)      :: params                               !> user defined parameter structure
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)             !> heavy data array - block data
    integer(kind=ik), intent(inout)     :: hvy_active(:)                        !> list of active blocks (heavy data)
    integer(kind=ik), intent(inout)     :: hvy_n                                !> number of active blocks (heavy data)
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)
    character(len=*), intent(in)        :: which_norm                           !> which norm to use ? "L2", "Linfty"
    real(kind=rk), intent(inout)        :: norm(:)                              !> the computed norm for each component of the vector
    real(kind=rk)                       :: x0(1:3), dx(1:3)
    integer(kind=ik) :: k, hvy_id, n_eqn, Bs(1:3), g, p, mpierr, lgt_id, D, J

    ! note: if norm and hvy_block components are of different size, we use the smaller one.
    n_eqn = min( size(norm, 1), size(hvy_block,4) )
    Bs = params%Bs
    g = params%n_ghosts
    D = params%dim

    select case (which_norm)
    case ("L2")
        norm = 0.0_rk
        if (params%dim == 2) then
            do k = 1, hvy_n
                hvy_id = hvy_active(k)
                call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

                J = lgt_block(lgt_id, params%max_treelevel+IDX_MESH_LVL)
                call get_block_spacing_origin2( lgt_block(lgt_id, 1:J), params%domain_size, Bs, D, x0, dx )

                do p = 1, n_eqn
                    ! attention: exclude redundant point here
                    norm(p) = norm(p) + dx(1)*dx(2)*sum( hvy_block(g+1:Bs(1)+g-1, g+1:Bs(2)+g-1, 1, p, hvy_id )**2 )
                enddo
            enddo
        else
            do k = 1, hvy_n
                hvy_id = hvy_active(k)
                call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

                J = lgt_block(lgt_id, params%max_treelevel+IDX_MESH_LVL)
                call get_block_spacing_origin2( lgt_block(lgt_id, 1:J), params%domain_size, Bs, D, x0, dx )

                do p = 1, n_eqn
                    ! attention: exclude redundant point here
                    norm(p) = norm(p) + dx(1)*dx(2)*dx(3)*sum( hvy_block(g+1:Bs(1)+g-1, g+1:Bs(2)+g-1, g+1:Bs(3)+g-1, p, hvy_id )**2 )
                enddo
            enddo
        endif

        call MPI_ALLREDUCE(MPI_IN_PLACE, norm, n_eqn, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
        ! norm = sqrt( norm )
        norm = sqrt( norm )

    case ("Linfty")
        norm = -1.0_rk
        do k = 1, hvy_n
            hvy_id = hvy_active(k)
            do p = 1, n_eqn
                ! yes, we can include the ghost nodes: it does not matter for the infty
                ! norm.
                norm(p) = max( norm(p), maxval( abs(hvy_block(:,:,:,p,hvy_id))) )
            enddo
        enddo

        call MPI_ALLREDUCE(MPI_IN_PLACE, norm, n_eqn, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)

    case default
        write(*,*) which_norm
        call abort(20030201, "The tree norm you desire is not implemented. How dare you.")

    end select

end subroutine


end module module_operators
