!-------------------------------------------------------------------------------
!> \brief FORTRAN FFT module for WABBIT
!-------------------------------------------------------------------------------
!> \details
!! This module implements some FFT functions, which can be used on a single block only (!!)\n
!! For the case Jmin=0, this block is usually periodic. \n
!! You are correct to thu_hat, that all this is of no great use for WABBIT, but we can utilize it for the lowest level solution in the Poisson solver. \n
!! This module needs FFTW to be installed. \n
!! This module is designed as blackbox, users should only use it's function with input and output real data, \n
!! and not need to touch the complex-valued Fourier-coefficients.
!-------------------------------------------------------------------------------
!  Thanks to the FLUSI coders for the inspiration for this code
!-------------------------------------------------------------------------------
module module_poisson

    use mpi
    use module_globals
    use module_params
    use module_helpers
    use module_forestMetaData
    use module_timing               ! debug module
    use module_treelib              ! module with evrything related to treecodes (encoding, decoding, neighbors, etc)

    implicit none

    ! I usually find it helpful to use the private keyword by itself initially, which specifies
    ! that everything within the module is private unless explicitly marked public.
    PRIVATE

    ! everything is save by default
    SAVE

!---------------------------------------------------------------------------------------------
! public parts of this module
    
    PUBLIC :: GS_iteration_level, GS_iteration_ref, GS_compute_residual, CG_solve_poisson_level0


contains

    subroutine GS_iteration_level(params, tree_id, level, u, b, a, stencil, sweep_forward)
        implicit none

        !> parameter struct
        type (type_params), intent(inout)  :: params
        integer(kind=ik), intent(in)       :: tree_id
        real(kind=rk), intent(inout)       :: u(:, :, :, :, :)
        real(kind=rk), intent(in)          :: b(:, :, :, :, :)
        integer(kind=ik), intent(in)       :: a
        real(kind=rk), intent(in)          :: stencil(-a:a)
        logical, intent(in)                :: sweep_forward
        integer(kind=ik), intent(in)       :: level

        integer(kind=ik)     :: k, hvy_id, lgt_id
        real(kind=rk)        :: t_block, x0(1:3), dx(1:3)

        t_block = MPI_Wtime()
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k, tree_ID)
            call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

            ! skip blocks not on this level
            if (lgt_block(lgt_id,IDX_MESH_LVL) /= level) cycle

            ! get spacing
            call get_block_spacing_origin_b( get_tc(lgt_block(lgt_id, IDX_TC_1 : IDX_TC_2)), params%domain_size, &
            params%Bs, x0, dx, dim=params%dim, level=lgt_block(lgt_id, IDX_MESH_LVL), max_level=params%Jmax)
            
            ! blocks on this level do a GS-sweep
            call GS_iteration(params, u(:,:,:,:,hvy_id), b(:,:,:,:,hvy_id), a, stencil, dx, sweep_forward)
        enddo
        call toc( "Gauss-Seidel iteration", 10006, MPI_Wtime()-t_block )
    end subroutine


    !> \brief Gauss-Seidel iteration for a given tree_id and refinement level
    subroutine GS_iteration_ref(params, tree_id, ref, u, b, a, stencil, sweep_forward, filter_offset)
        implicit none

        !> parameter struct
        type (type_params), intent(inout)  :: params
        integer(kind=ik), intent(in)       :: tree_id
        integer(kind=ik), intent(in)       :: ref(:)  ! can be several, input as (/ VAL1, VAL2, ... /)
        real(kind=rk), intent(inout)       :: u(:, :, :, :, :)
        real(kind=rk), intent(in)          :: b(:, :, :, :, :)
        integer(kind=ik), intent(in)       :: a
        real(kind=rk), intent(in)          :: stencil(-a:a)
        logical, intent(in)                :: sweep_forward
        integer(kind=ik), intent(in), optional  :: filter_offset  ! where to apply the filter, default is params%g resulting in only interior points

        integer(kind=ik)     :: k, hvy_id, lgt_id, filterOffset
        real(kind=rk)        :: t_block, x0(1:3), dx(1:3)

        filterOffset = params%g
        if (present(filter_offset)) filterOffset = filter_offset

        t_block = MPI_Wtime()
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k, tree_ID)
            call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )

            ! skip blocks not on this level
            if (all(lgt_block(lgt_id,IDX_REFINE_STS) /= ref)) cycle

            ! get spacing
            call get_block_spacing_origin_b( get_tc(lgt_block(lgt_id, IDX_TC_1 : IDX_TC_2)), params%domain_size, &
            params%Bs, x0, dx, dim=params%dim, level=lgt_block(lgt_id, IDX_MESH_LVL), max_level=params%Jmax)
            
            ! blocks on this level do a GS-sweep
            call GS_iteration(params, u(:,:,:,:,hvy_id), b(:,:,:,:,hvy_id), a, stencil, dx, sweep_forward, filterOffset)
        enddo
        call toc( "Gauss-Seidel iteration", 10006, MPI_Wtime()-t_block )
    end subroutine

    ! do one Gauss Seidel iteration, either in backwards or forwards fashion, to solve Ax=b
    ! A is supposed to be periodic and sparse with given stencil
    subroutine GS_iteration(params, u, b, a, stencil, dx, sweep_forward, filter_offset)
        implicit none

        !> parameter struct
        type (type_params), intent(inout)  :: params
        real(kind=rk), intent(inout)       :: u(:, :, :, :)
        real(kind=rk), intent(in)          :: b(:, :, :, :)
        integer(kind=ik), intent(in)       :: a
        real(kind=rk), intent(in)          :: stencil(-a:a)
        real(kind=rk), intent(in)          :: dx(:)           ! needed to weight b, stencil is assumed to not be adapted by dx
        logical, intent(in)                :: sweep_forward
        integer(kind=ik), intent(in), optional :: filter_offset  ! where to apply the filter, default is params%g resulting in only interior points

        real(kind=rk)                      :: stencil_GS(-a:a), fac_GS
        integer(kind=ik)                   :: ix,iy,iz, ic, filterOffset

        filterOffset = params%g
        if (present(filter_offset)) filterOffset = filter_offset

        ! prepare GS stencil
        fac_GS = stencil(0)*params%dim            ! diagonal factor that we need to divide b with
        stencil_GS = -stencil / fac_GS  ! inverted stencil for GS application
        stencil_GS(0) = 0.0_rk         ! main diagonal needs to be 0

        do ic = 1, size(u,4)
            ! forward sweep
            if (sweep_forward) then
                if (params%dim == 3) then
                    do iz = filterOffset+1, params%bs(3)+2*params%g-filterOffset, 1
                        do iy = filterOffset+1, params%bs(2)+2*params%g-filterOffset, 1
                            do ix = filterOffset+1, params%bs(1)+2*params%g-filterOffset, 1
                                u(ix,iy,iz,ic) = sum(stencil_GS * u(ix-a:ix+a,iy,iz,ic)) + sum(stencil_GS * u(ix,iy-a:iy+a,iz,ic)) + sum(stencil_GS * u(ix,iy,iz-a:iz+a,ic)) + b(ix,iy,iz,ic) / fac_GS * dx(1)*dx(1)
                            enddo
                        enddo
                    enddo
                else
                    do iy = filterOffset+1, params%bs(2)+2*params%g-filterOffset, 1
                        do ix = filterOffset+1, params%bs(1)+2*params%g-filterOffset, 1
                            u(ix,iy,1,ic) = sum(stencil_GS * u(ix-a:ix+a,iy,1,ic)) + sum(stencil_GS * u(ix,iy-a:iy+a,1,ic)) + b(ix,iy,1,ic) / fac_GS * dx(1)*dx(1)
                        enddo
                    enddo
                endif
            else
                if (params%dim == 3) then
                    do iz = params%bs(3)+2*params%g-filterOffset, filterOffset+1, -1
                        do iy = params%bs(2)+2*params%g-filterOffset, filterOffset+1, -1
                            do ix = params%bs(1)+2*params%g-filterOffset, filterOffset+1, -1
                                u(ix,iy,iz,ic) = sum(stencil_GS * u(ix-a:ix+a,iy,iz,ic)) + sum(stencil_GS * u(ix,iy-a:iy+a,iz,ic)) + sum(stencil_GS * u(ix,iy,iz-a:iz+a,ic)) + b(ix,iy,iz,ic) / fac_GS * dx(1)*dx(1)
                            enddo
                        enddo
                    enddo
                else
                    do iy = params%bs(2)+2*params%g-filterOffset, filterOffset+1, -1
                        do ix = params%bs(1)+2*params%g-filterOffset, filterOffset+1, -1
                            u(ix,iy,1,ic) = sum(stencil_GS * u(ix-a:ix+a,iy,1,ic)) + sum(stencil_GS * u(ix,iy-a:iy+a,1,ic)) + b(ix,iy,1,ic) / fac_GS * dx(1)*dx(1)
                        enddo
                    enddo
                endif
            endif
        enddo

    end subroutine

    !> \brief compute residual r = b - Ax
    !> \details Compute the residual r = b - Ax, where A is the operator defined by the stencil.
    !> Can also be used to recompute b as b = r + Ax
    subroutine GS_compute_residual(params, u, b, r, a, stencil, dx, recompute_b)
        implicit none

        !> parameter struct
        type (type_params), intent(inout)  :: params
        real(kind=rk), intent(in)         :: u(:, :, :, :)
        real(kind=rk), intent(in)         :: b(:, :, :, :)  !< RHS b
        real(kind=rk), intent(out)        :: r(:, :, :, :)  !< residual r, can be inplace with b
        integer(kind=ik), intent(in)      :: a
        real(kind=rk), intent(in)         :: stencil(-a:a)
        real(kind=rk), intent(in)         :: dx(:)
        logical, intent(in), optional     :: recompute_b    !< if .true. then b = r + Ax

        integer :: Ax_factor

        real(kind=rk)                        :: stencil_x(-a:a)
        real(kind=rk)                        :: stencil_y(-a:a)
        real(kind=rk)                        :: stencil_z(-a:a)
        integer(kind=ik)                     :: ix,iy,iz,ic

        ! factor to multiply Ax with, decides if we compute r = b - Ax or recompute b = r + Ax
        Ax_factor = -1
        if (present(recompute_b)) then
            if (recompute_b) Ax_factor = 1
        endif

        stencil_x = stencil / (dx(1)*dx(1))
        stencil_y = stencil / (dx(2)*dx(2))
        if (params%dim==3) then
            stencil_z = stencil / (dx(3)*dx(3))
        endif

        do ic = 1,size(u,4)
            if (params%dim == 3) then
                do iz = params%g+1, params%g+params%bs(3)
                    do iy = params%g+1, params%g+params%bs(2)
                        do ix = params%g+1, params%g+params%bs(1)
                            r(ix,iy,iz,ic) = b(ix,iy,iz,ic) + Ax_factor * (sum(stencil_x * u(ix-a:ix+a,iy,iz,ic)) \
                                + sum(stencil_y * u(ix,iy-a:iy+a,iz,ic)) + sum(stencil_z * u(ix,iy,iz-a:iz+a,ic)))
                        enddo
                    enddo
                enddo
            else
                do iy = params%g+1, params%g+params%bs(2)
                    do ix = params%g+1, params%g+params%bs(1)
                        r(ix,iy,1,ic) = b(ix,iy,1,ic) + Ax_factor * (sum(stencil_x * u(ix-a:ix+a,iy,1,ic)) \
                            + sum(stencil_y * u(ix,iy-a:iy+a,1,ic)))
                    enddo
                enddo
            endif
        enddo
    end subroutine


    ! apply laplacian and compute Ax
    subroutine GS_compute_laplacian(params, u, ddu, a, stencil, dx)
        implicit none

        !> parameter struct
        type (type_params), intent(inout) :: params
        real(kind=rk), intent(in)         :: u(:, :, :)
        real(kind=rk), intent(inout)      :: ddu(:, :, :)
        integer(kind=ik), intent(in)      :: a
        real(kind=rk), intent(in)         :: stencil(-a:a)
        real(kind=rk), intent(in)         :: dx(:)

        real(kind=rk)                     :: stencil_x(-a:a)
        real(kind=rk)                     :: stencil_y(-a:a)
        real(kind=rk)                     :: stencil_z(-a:a)
        integer(kind=ik)                  :: ix,iy,iz

        stencil_x = stencil / (dx(1)*dx(1))
        stencil_y = stencil / (dx(2)*dx(2))
        if (params%dim==3) then
            stencil_z = stencil / (dx(3)*dx(3))
        endif

        if (params%dim == 3) then
            do iz = params%g+1, params%g+params%bs(3)
                do iy = params%g+1, params%g+params%bs(2)
                    do ix = params%g+1, params%g+params%bs(1)
                        ddu(ix,iy,iz) = sum(stencil_x * u(ix-a:ix+a,iy,iz)) \
                            + sum(stencil_y * u(ix,iy-a:iy+a,iz)) + sum(stencil_z * u(ix,iy,iz-a:iz+a))
                    enddo
                enddo
            enddo
        else
            do iy = params%g+1, params%g+params%bs(2)
                do ix = params%g+1, params%g+params%bs(1)
                    ddu(ix,iy,1) =  sum(stencil_x * u(ix-a:ix+a,iy,1)) \
                        + sum(stencil_y * u(ix,iy-a:iy+a,1))
                enddo
            enddo
        endif
    end subroutine



    subroutine CG_solve_poisson_level0(params, u, f, a, stencil, dx, tol, max_iter)
        use module_params
        implicit none

        ! Arguments
        type (type_params), intent(inout)  :: params
        real(kind=rk), intent(inout) :: u(:, :, :, :)
        real(kind=rk), intent(in) :: f(:, :, :, :)
        integer(kind=ik), intent(in)      :: a
        real(kind=rk), intent(in)         :: stencil(-a:a)
        real(kind=rk), intent(in)         :: dx(:)
        real(kind=rk), intent(in) :: tol
        integer, intent(in) :: max_iter

        ! Locals
        real(kind=rk), allocatable, save :: r(:,:,:), p(:,:,:), Ap(:,:,:)
        real(kind=rk) :: alpha, beta, rsold, rsnew
        integer(kind=ik) :: ic, it, nx, ny, nz, nc, g(1:3), bs(1:3)

        character(len=cshort)              :: fname

        nx = size(u,1)
        ny = size(u,2)
        nz = size(u,3)
        nc = size(u,4)
        g = 0
        g(1:params%dim) = params%g
        bs = params%bs

        if (.not. allocated(r)) allocate(r(nx, ny, nz))
        if (.not. allocated(p)) allocate(p(nx, ny, nz))
        if (.not. allocated(Ap)) allocate(Ap(nx, ny, nz))

        do ic = 1,nc
            ! sync u and to start, maybe they have not been synched before
            ! f does not need to be synched, as no stencil will ever be applied on it
            ! assume that we have only one periodic block here and laplacian is only a + shaped stencil
            call sync_block_level0(params, u(:,:,:,ic), sync_sides_only=.true.)

            ! write(fname, '(A, i0, A,i0, A)') "block_dumped_it", 0, ".t"
            ! call dump_block_fancy(u(:,:,:,:), fname, params%Bs, params%g, digits=2)

            ! r = f - A u
            call GS_compute_laplacian(params, u(:,:,:,ic), Ap, a, stencil, dx)
            r(g(1)+1:bs(1)+g(1),g(2)+1:bs(2)+g(2),g(3)+1:bs(3)+g(3)) = f(g(1)+1:bs(1)+g(1),g(2)+1:bs(2)+g(2),g(3)+1:bs(3)+g(3),ic) - Ap(g(1)+1:bs(1)+g(1),g(2)+1:bs(2)+g(2),g(3)+1:bs(3)+g(3))
            p(g(1)+1:bs(1)+g(1),g(2)+1:bs(2)+g(2),g(3)+1:bs(3)+g(3)) = r(g(1)+1:bs(1)+g(1),g(2)+1:bs(2)+g(2),g(3)+1:bs(3)+g(3))
            rsold = sum(r(g(1)+1:bs(1)+g(1),g(2)+1:bs(2)+g(2),g(3)+1:bs(3)+g(3))**2)
            ! write(*, "(A, i4, A, es9.2)") "      CG it ", 0, " res= ", rsold
            do it = 1, max_iter
                ! sync only p to compute Ap
                ! assume that we have only one periodic block here and laplacian is only a + shaped stencil
                call sync_block_level0(params, p, sync_sides_only=.true.)

                call GS_compute_laplacian(params, p, Ap, a, stencil, dx)
                alpha = rsold / sum(p(g(1)+1:bs(1)+g(1),g(2)+1:bs(2)+g(2),g(3)+1:bs(3)+g(3)) * Ap(g(1)+1:bs(1)+g(1),g(2)+1:bs(2)+g(2),g(3)+1:bs(3)+g(3)))
                u(g(1)+1:bs(1)+g(1),g(2)+1:bs(2)+g(2),g(3)+1:bs(3)+g(3),ic) = u(g(1)+1:bs(1)+g(1),g(2)+1:bs(2)+g(2),g(3)+1:bs(3)+g(3),ic) + alpha * p(g(1)+1:bs(1)+g(1),g(2)+1:bs(2)+g(2),g(3)+1:bs(3)+g(3))
                r(g(1)+1:bs(1)+g(1),g(2)+1:bs(2)+g(2),g(3)+1:bs(3)+g(3)) = r(g(1)+1:bs(1)+g(1),g(2)+1:bs(2)+g(2),g(3)+1:bs(3)+g(3)) - alpha * Ap(g(1)+1:bs(1)+g(1),g(2)+1:bs(2)+g(2),g(3)+1:bs(3)+g(3))
                rsnew = sum(r(g(1)+1:bs(1)+g(1),g(2)+1:bs(2)+g(2),g(3)+1:bs(3)+g(3))**2)
        
                ! write(*, "(A, i4, A, es9.2)") "      CG it ", it, " res= ", sqrt(rsnew)

                if (sqrt(rsnew) < tol) exit
        
                beta = rsnew / rsold
                p(g(1)+1:bs(1)+g(1),g(2)+1:bs(2)+g(2),g(3)+1:bs(3)+g(3)) = r(g(1)+1:bs(1)+g(1),g(2)+1:bs(2)+g(2),g(3)+1:bs(3)+g(3)) + beta * p(g(1)+1:bs(1)+g(1),g(2)+1:bs(2)+g(2),g(3)+1:bs(3)+g(3))
                rsold = rsnew

                ! write(fname, '(A, i0, A,i0, A)') "block_dumped_it", it, ".t"
                ! call dump_block_fancy(u(:,:,:,:), fname, params%Bs, params%g, digits=2)
            enddo

            ! sync only u to finish
            call sync_block_level0(params, u(:,:,:,ic), sync_sides_only=.true.)
        enddo

    end subroutine CG_solve_poisson_level0

    ! synch routine for block on lowest level
    subroutine sync_block_level0(params, u, sync_sides_only)
        use module_params
        implicit none

        ! Arguments
        type (type_params), intent(inout)  :: params
        real(kind=rk), intent(inout) :: u(:, :, :)
        logical, intent(in), optional :: sync_sides_only
        logical :: sync_sides_only_flag
        integer(kind=ik) :: is(1:3)=1, ie(1:3)

        ! select indices, for laplacian we have a +-shaped stencil, so we do not need to sync corners and edges
        sync_sides_only_flag = .false.
        if (present(sync_sides_only)) sync_sides_only_flag = sync_sides_only
        if (sync_sides_only_flag) then
            is(1:params%dim) = params%g+1
            ie(1:params%dim) = params%g+params%bs(1:params%dim)
        else
            ie =params%bs+2*params%g
            if (params%dim == 2) ie(3) = 1
        endif
        
        ! sync and assume that we have only one periodic block here
        u(1:params%g                                     ,is(2):ie(2),is(3):ie(3)) = u(1+params%bs(1):params%g+params%bs(1),is(2):ie(2),is(3):ie(3))
        u(params%g+params%bs(1)+1:params%g*2+params%bs(1),is(2):ie(2),is(3):ie(3)) = u(params%g+1:2*params%g               ,is(2):ie(2),is(3):ie(3))
        u(is(1):ie(1),1:params%g                                     ,is(3):ie(3)) = u(is(1):ie(1),1+params%bs(2):params%g+params%bs(2),is(3):ie(3))
        u(is(1):ie(1),params%g+params%bs(2)+1:params%g*2+params%bs(2),is(3):ie(3)) = u(is(1):ie(1),params%g+1:2*params%g               ,is(3):ie(3))
        if (params%dim == 3) then
            u(is(1):ie(1),is(2):ie(2),1:params%g)                                      = u(is(1):ie(1),is(2):ie(2),1+params%bs(3):params%g+params%bs(3))
            u(is(1):ie(1),is(2):ie(2),params%g+params%bs(3)+1:params%g*2+params%bs(3)) = u(is(1):ie(1),is(2):ie(2),params%g+1:2*params%g)
        endif

    end subroutine


end module module_poisson