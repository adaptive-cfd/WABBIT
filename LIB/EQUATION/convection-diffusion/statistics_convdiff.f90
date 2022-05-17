!-----------------------------------------------------------------------------
! main level wrapper to compute statistics (such as mean flow, global energy,
! forces, but potentially also derived stuff such as Integral/Kolmogorov scales)
! NOTE: as for the RHS, some terms here depend on the grid as whole, and not just
! on individual blocks. This requires one to use the same staging concept as for the RHS.
!-----------------------------------------------------------------------------
subroutine STATISTICS_convdiff( time, dt, u, g, x0, dx, stage)
    implicit none

    ! it may happen that some source terms have an explicit time-dependency
    ! therefore the general call has to pass time
    real(kind=rk), intent (in) :: time, dt

    ! block data, containg the state vector. In general a 4D field (3 dims+components)
    ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
    real(kind=rk), intent(inout) :: u(1:,1:,1:,1:)

    ! as you are allowed to compute the RHS only in the interior of the field
    ! you also need to know where 'interior' starts: so we pass the number of ghost points
    integer, intent(in) :: g

    ! for each block, you'll need to know where it lies in physical space. The first
    ! non-ghost point has the coordinate x0, from then on its just cartesian with dx spacing
    real(kind=rk), intent(in) :: x0(1:3), dx(1:3)

    ! stage. there is 3 stages, init_stage, integral_stage and local_stage. If the PDE has
    ! terms that depend on global qtys, such as forces etc, which cannot be computed
    ! from a single block alone, the first stage does that. the second stage can then
    ! use these integral qtys for the actual RHS evaluation.
    character(len=*), intent(in) :: stage

    ! local variables
    integer(kind=ik) :: mpierr, ix, iy, iz, k
    integer(kind=ik), dimension(3) :: Bs
    real(kind=rk) :: scalar_integral, scalar_max
    real(kind=rk), save :: umag, umax, dx_min


    ! compute the size of blocks
    Bs(1) = size(u,1) - 2*g
    Bs(2) = size(u,2) - 2*g
    Bs(3) = size(u,3) - 2*g

    select case(stage)
    case ("init_stage")
        !-------------------------------------------------------------------------
        ! 1st stage: init_stage.
        !-------------------------------------------------------------------------
        ! this stage is called only once, NOT for each block.
        ! performs initializations in the RHS module, such as resetting integrals
        params_convdiff%scalar_integral = 0.0_rk
        params_convdiff%scalar_max = 0.0_rk

    case ("integral_stage")
        !-------------------------------------------------------------------------
        ! 2nd stage: integral_stage.
        !-------------------------------------------------------------------------
        ! This stage contains all operations which are running on the blocks
        ! called for each block.

        ! tmp values for computing the current block only
        scalar_integral = 0.0_rk
        scalar_max = -9.0e9_rk

        if (params_convdiff%dim == 2) then
            ! --- 2D --- --- 2D --- --- 2D --- --- 2D --- --- 2D --- --- 2D ---
            do iy = g+1, Bs(2)+g-1 ! Note: loops skip redundant points
                do ix = g+1, Bs(1)+g-1
                    scalar_integral = scalar_integral + u(ix,iy,1,1)
                enddo
            enddo
        else
            ! --- 3D --- --- 3D --- --- 3D --- --- 3D --- --- 3D --- --- 3D ---
            do iz = g+1, Bs(3)+g-1 ! Note: loops skip redundant points
                do iy = g+1, Bs(2)+g-1
                    do ix = g+1, Bs(1)+g-1
                        scalar_integral = scalar_integral + u(ix,iy,iz,1)
                    enddo
                enddo
            enddo
        endif

        scalar_max = maxval(u(:,:,:,1))

        ! we just computed the values on the current block, which we now add to the
        ! existing blocks in the variables (recall normalization by dV)
        params_convdiff%scalar_integral = params_convdiff%scalar_integral + scalar_integral * product(dx(1:params_convdiff%dim))
        params_convdiff%scalar_max = max(params_convdiff%scalar_max, scalar_max)
    case ("post_stage")
        !-------------------------------------------------------------------------
        ! 3rd stage: post_stage.
        !-------------------------------------------------------------------------
        ! this stage is called only once, NOT for each block.

        call MPI_ALLREDUCE(MPI_IN_PLACE, params_convdiff%scalar_integral, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
        call MPI_ALLREDUCE(MPI_IN_PLACE, params_convdiff%scalar_max     , 1, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)

        call append_t_file( 'scalar_integral.t', (/time, dt, params_convdiff%scalar_integral, params_convdiff%scalar_max /) )

    case default
        call abort(7772,"the STATISTICS wrapper requests a stage this physics module cannot handle.")

    end select


end subroutine
