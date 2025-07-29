!-------------------------------------------------------------------------------
!> \brief FORTRAN poisson solver module for WABBIT
!-------------------------------------------------------------------------------
!> \details
!! ToDo: Explanations
!! ATTENTION! This module might assume equidistant grid spacing, which might not always be the case. Needs to be checked!
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

        ! we want to solve Ax = Bb where
        !       A is the discrete Laplacian operator
        !       B is an auxiliary operator for the RHS, used in the Mehrstellenverfahren
        !       u is the solution vector
        !       b is the right-hand side
        ! Here we define the different stencils for each method, which then will later be applied

        ! different stencils, either for one-dimensional cross-stencil form or tensorial form
        real(kind=rk), allocatable  :: stencil(:), stencil_tensor_2D(:,:), stencil_tensor_3D(:,:,:), stencil_RHS_tensor_2D(:,:), stencil_RHS_tensor_3D(:,:,:)
        real(kind=rk)  :: stencil_RHS
        integer(kind=ik) :: stencil_size, stencil_RHS_size
        logical :: use_tensor



!---------------------------------------------------------------------------------------------
! public parts of this module
    
    PUBLIC :: GS_iteration_level, GS_iteration_ref, GS_compute_residual, GS_compute_Ax, CG_solve_poisson_level0, setup_Laplacian_stencils


contains

    subroutine setup_Laplacian_stencils(params, g)
        implicit none
        !> parameter struct
        type (type_params), intent(inout)  :: params
        integer(kind=ik), intent(inout), optional :: g  ! ghost points, if present then we set it if we require more points

        ! prepare stencil
        if (allocated(stencil)) deallocate(stencil)
        if (allocated(stencil_tensor_2D)) deallocate(stencil_tensor_2D)
        if (allocated(stencil_tensor_3D)) deallocate(stencil_tensor_3D)
        if (allocated(stencil_RHS_tensor_2D)) deallocate(stencil_RHS_tensor_2D)
        if (allocated(stencil_RHS_tensor_3D)) deallocate(stencil_RHS_tensor_3D)

        if (params%laplacian_order == "FD_2nd_central") then
            allocate(stencil(-1:1))
            stencil = (/  1.0_rk, -2.0_rk,    1.0_rk /)
            stencil_RHS = 1.0_rk
            params%laplacian_stencil_size = 1
            stencil_size = 1
            stencil_RHS_size = 0
            use_tensor = .false.
        elseif (params%laplacian_order == "FD_4th_central") then
            allocate(stencil(-2:2))
            stencil = (/ -1.0_rk,  16.0_rk,  -30.0_rk,   16.0_rk,   -1.0_rk /) / 12.0_rk
            stencil_RHS = 1.0_rk
            params%laplacian_stencil_size = 2
            stencil_size = 2
            stencil_RHS_size = 0
            use_tensor = .false.
        elseif (params%laplacian_order == "FD_6th_central") then
            allocate(stencil(-3:3))
            stencil = (/  2.0_rk, -27.0_rk,   270.0_rk, -490.0_rk,   270.0_rk,  -27.0_rk,    2.0_rk/) / 180.0_rk
            stencil_RHS = 1.0_rk
            params%laplacian_stencil_size = 3
            stencil_size = 3
            stencil_RHS_size = 0
            use_tensor = .false.
        elseif (params%laplacian_order == "FD_8th_central") then
            allocate(stencil(-4:4))
            stencil = (/ -9.0_rk,  128.0_rk, -1008.0_rk, 8064.0_rk, -14350.0_rk, 8064.0_rk, -1008.0_rk, 128.0_rk, -9.0_rk /) / 5040_rk
            stencil_RHS = 1.0_rk
            params%laplacian_stencil_size = 4
            stencil_size = 4
            stencil_RHS_size = 0
            use_tensor = .false.
        elseif (params%laplacian_order == "FD_4th_comp_0_4") then
            allocate(stencil(-4:4))
            stencil = (/  -75.0_rk, 544.0_rk, -1776.0_rk, 3552.0_rk, -4490.0_rk, 3552.0_rk, -1776.0_rk, 544.0_rk, -75.0_rk /) / 144.0_rk
            stencil_RHS = 1.0_rk
            params%laplacian_stencil_size = 4
            stencil_size = 4
            stencil_RHS_size = 0
            use_tensor = .false.
        elseif (params%laplacian_order == "FD_4th_comp_2_2") then
            allocate(stencil(-4:4))
            stencil = (/  1.0_rk, -16.0_rk, 64.0_rk, 16.0_rk, -130.0_rk, 16.0_rk, 64.0_rk, -16.0_rk, 1.0_rk /) / 144.0_rk
            stencil_RHS = 1.0_rk
            params%laplacian_stencil_size = 4
            stencil_size = 4
            stencil_RHS_size = 0
            use_tensor = .false.
        elseif (params%laplacian_order == "FD_4th_comp_1_3") then
            allocate(stencil(-4:4))
            stencil = (/  3.0_rk, -8.0_rk, -24.0_rk, 264.0_rk, -470.0_rk, 264.0_rk, -24.0_rk, -8.0_rk, 3.0_rk /) / 144.0_rk
            stencil_RHS = 1.0_rk
            params%laplacian_stencil_size = 4
            stencil_size = 4
            stencil_RHS_size = 0
            use_tensor = .false.
        elseif (params%laplacian_order == "FD_6th_comp_3_3") then
            allocate(stencil(-6:6))
            stencil = (/  1.0_rk, -18.0_rk, 171.0_rk, -810.0_rk, 1935.0_rk, 828.0_rk, -4214.0_rk, 828.0_rk, 1935.0_rk, -810.0_rk, 171.0_rk, -18.0_rk, 1.0_rk /) / 3600.0_rk
            stencil_RHS = 1.0_rk
            params%laplacian_stencil_size = 6
            stencil_size = 6
            stencil_RHS_size = 0
            use_tensor = .false.
        elseif (params%laplacian_order == "FD_6th_comp_2_4") then
            allocate(stencil(-6:6))
            stencil = (/  2.0_rk, -40.0_rk, 217.0_rk, -520.0_rk, 270.0_rk, 4656.0_rk, -9170.0_rk, 4656.0_rk, 270.0_rk, -520.0_rk, 217.0_rk, -40.0_rk, 2.0_rk /) / 3600.0_rk
            stencil_RHS = 1.0_rk
            params%laplacian_stencil_size = 6
            stencil_size = 6
            stencil_RHS_size = 0
            use_tensor = .false.
        elseif (params%laplacian_order == "FD_6th_comp_1_5") then
            allocate(stencil(-6:6))
            stencil = (/  20.0_rk, 4.0_rk, -955.0_rk, 5300.0_rk, -15300.0_rk, 31560.0_rk, -41258.0_rk, 31560.0_rk, -15300.0_rk, 5300.0_rk, -955.0_rk, 4.0_rk, 20.0_rk /) / 3600.0_rk
            stencil_RHS = 1.0_rk
            params%laplacian_stencil_size = 6
            stencil_size = 6
            stencil_RHS_size = 0
            use_tensor = .false.
        elseif (params%laplacian_order == "FD_6th_comp_0_6") then
            allocate(stencil(-6:6))
            stencil = (/  -1470.0_rk, 14184.0_rk, -63495.0_rk, 176200.0_rk, -342450.0_rk, 501840.0_rk, -569618.0_rk, 501840.0_rk, -342450.0_rk, 176200.0_rk, -63495.0_rk, 14184.0_rk, -1470.0_rk /) / 3600.0_rk
            stencil_RHS = 1.0_rk
            params%laplacian_stencil_size = 6
            stencil_size = 6
            stencil_RHS_size = 0
            use_tensor = .false.
        elseif (params%laplacian_order == "FD_6th_mehrstellen") then
            if (params%dim == 2) then
                allocate(stencil_tensor_2D(-1:1, -1:1))
                allocate(stencil_RHS_tensor_2D(-2:2, -2:2))
                stencil_tensor_2D = reshape( (/ &
                    1.0_rk,   4.0_rk, 1.0_rk, &
                    4.0_rk, -20.0_rk, 4.0_rk, &
                    1.0_rk,   4.0_rk, 1.0_rk /), (/3,3/) ) / 6.0_rk
                stencil_RHS_tensor_2D = reshape( (/ &
                     0.0_rk,  0.0_rk,  -3.0_rk,  0.0_rk,  0.0_rk, &
                     0.0_rk,  8.0_rk,  56.0_rk,  8.0_rk,  0.0_rk, &
                    -3.0_rk, 56.0_rk, 476.0_rk, 56.0_rk, -3.0_rk, &
                     0.0_rk,  8.0_rk,  56.0_rk,  8.0_rk,  0.0_rk, &
                     0.0_rk,  0.0_rk,  -3.0_rk,  0.0_rk,  0.0_rk /), (/5,5/) ) / 720.0_rk
            else
                allocate(stencil_tensor_3D(-1:1, -1:1, -1:1))
                allocate(stencil_RHS_tensor_3D(-2:2, -2:2, -2:2))
                stencil_tensor_3D = reshape( (/ &
                    1.0_rk,  3.0_rk, 1.0_rk,  3.0_rk,   14.0_rk,  3.0_rk, 1.0_rk,  3.0_rk, 1.0_rk, &
                    3.0_rk, 14.0_rk, 3.0_rk, 14.0_rk, -128.0_rk, 14.0_rk, 3.0_rk, 14.0_rk, 3.0_rk, &
                    1.0_rk,  3.0_rk, 1.0_rk,  3.0_rk,   14.0_rk,  3.0_rk, 1.0_rk,  3.0_rk, 1.0_rk /), (/3,3,3/) ) / 30.0_rk
                stencil_RHS_tensor_3D = reshape( (/ &
                0.0_rk,  0.0_rk,   0.0_rk,  0.0_rk,  0.0_rk, &
                0.0_rk,  0.0_rk,   0.0_rk,  0.0_rk,  0.0_rk, &
                0.0_rk,  0.0_rk,  -3.0_rk,  0.0_rk,  0.0_rk, &
                0.0_rk,  0.0_rk,   0.0_rk,  0.0_rk,  0.0_rk, &
                0.0_rk,  0.0_rk,   0.0_rk,  0.0_rk,  0.0_rk, &

                0.0_rk,  0.0_rk,   0.0_rk,  0.0_rk,  0.0_rk, &
                0.0_rk,  0.0_rk,   8.0_rk,  0.0_rk,  0.0_rk, &
                0.0_rk,  8.0_rk,  40.0_rk,  8.0_rk,  0.0_rk, &
                0.0_rk,  0.0_rk,   8.0_rk,  0.0_rk,  0.0_rk, &
                0.0_rk,  0.0_rk,   0.0_rk,  0.0_rk,  0.0_rk, &

                0.0_rk,  0.0_rk,   -3.0_rk,  0.0_rk,  0.0_rk, &
                0.0_rk,  8.0_rk,   40.0_rk,  8.0_rk,  0.0_rk, &
               -3.0_rk, 40.0_rk,  402.0_rk, 40.0_rk, -3.0_rk, &
                0.0_rk,  8.0_rk,   40.0_rk,  8.0_rk,  0.0_rk, &
                0.0_rk,  0.0_rk,   -3.0_rk,  0.0_rk,  0.0_rk, &

                0.0_rk,  0.0_rk,   0.0_rk,  0.0_rk,  0.0_rk, &
                0.0_rk,  0.0_rk,   8.0_rk,  0.0_rk,  0.0_rk, &
                0.0_rk,  8.0_rk,  40.0_rk,  8.0_rk,  0.0_rk, &
                0.0_rk,  0.0_rk,   8.0_rk,  0.0_rk,  0.0_rk, &
                0.0_rk,  0.0_rk,   0.0_rk,  0.0_rk,  0.0_rk, &

                0.0_rk,  0.0_rk,   0.0_rk,  0.0_rk,  0.0_rk, &
                0.0_rk,  0.0_rk,   0.0_rk,  0.0_rk,  0.0_rk, &
                0.0_rk,  0.0_rk,  -3.0_rk,  0.0_rk,  0.0_rk, &
                0.0_rk,  0.0_rk,   0.0_rk,  0.0_rk,  0.0_rk, &
                0.0_rk,  0.0_rk,   0.0_rk,  0.0_rk,  0.0_rk /), (/5,5,5/) ) / 720.0_rk
            endif
            params%laplacian_stencil_size = 1
            stencil_size = 1
            stencil_RHS_size = 2
            use_tensor = .true.
        else
            call abort(250612, "I don't know about this laplacian discretization order: "//trim(adjustl(params%laplacian_order))//" in GS_compute_residual")
        endif

        if (present(g)) then
            g = max(g, params%laplacian_stencil_size)  ! we need at least the stencil size
            g = max(g, stencil_RHS_size)  ! also take RHS into account
        endif
        
    end subroutine setup_Laplacian_stencils

    subroutine GS_iteration_level(params, tree_id, level, u, b, sweep_forward)
        implicit none

        !> parameter struct
        type (type_params), intent(inout)  :: params
        integer(kind=ik), intent(in)       :: tree_id
        real(kind=rk), intent(inout)       :: u(:, :, :, :, :)
        real(kind=rk), intent(in)          :: b(:, :, :, :, :)
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
            call GS_iteration(params, u(:,:,:,:,hvy_id), b(:,:,:,:,hvy_id), dx, sweep_forward)
        enddo
        call toc( "Gauss-Seidel iteration", 10006, MPI_Wtime()-t_block )
    end subroutine


    !> \brief Gauss-Seidel iteration for a given tree_id and refinement level
    subroutine GS_iteration_ref(params, tree_id, ref, u, b, sweep_forward, filter_offset)
        implicit none

        !> parameter struct
        type (type_params), intent(inout)  :: params
        integer(kind=ik), intent(in)       :: tree_id
        integer(kind=ik), intent(in)       :: ref(:)  ! can be several, input as (/ VAL1, VAL2, ... /)
        real(kind=rk), intent(inout)       :: u(:, :, :, :, :)
        real(kind=rk), intent(in)          :: b(:, :, :, :, :)
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
            call GS_iteration(params, u(:,:,:,:,hvy_id), b(:,:,:,:,hvy_id), dx, sweep_forward, filterOffset)
        enddo
        call toc( "Gauss-Seidel iteration", 10006, MPI_Wtime()-t_block )
    end subroutine

    ! do one Gauss Seidel iteration, either in backwards or forwards fashion, to solve Ax=b
    ! A is supposed to be periodic and sparse with given stencil
    subroutine GS_iteration(params, u, b, dx, sweep_forward, filter_offset)
        implicit none

        !> parameter struct
        type (type_params), intent(inout)  :: params
        real(kind=rk), intent(inout)       :: u(:, :, :, :)
        real(kind=rk), intent(in)          :: b(:, :, :, :)
        real(kind=rk), intent(in)          :: dx(:)           ! needed to weight b, stencil is assumed to not be adapted by dx
        logical, intent(in)                :: sweep_forward
        integer(kind=ik), intent(in), optional :: filter_offset  ! where to apply the filter, default is params%g resulting in only interior points

        integer(kind=ik)                   :: ix,iy,iz, ic, filterOffset, is(1:3), ie(1:3), dir, a

        filterOffset = params%g
        if (present(filter_offset)) filterOffset = filter_offset

        a = params%laplacian_stencil_size

        is(:) = 1
        ie(:) = 1
        if (sweep_forward) then
            is(1:params%dim) = filterOffset + 1
            ie(1:params%dim) = params%bs(1:params%dim) + 2*params%g - filterOffset
            dir = 1
        else
            is(1:params%dim) = params%bs(1:params%dim) + 2*params%g - filterOffset
            ie(1:params%dim) = filterOffset + 1
            dir = -1
        endif

        do ic = 1, size(u,4)
            if (.not. use_tensor) then
                if (params%dim == 2) then
                    do iy = is(2), ie(2), dir
                        do ix = is(1), ie(1), dir
                            u(ix,iy,1,ic) = (-sum(stencil * u(ix-a:ix+a,iy,1,ic)) - sum(stencil * u(ix,iy-a:iy+a,1,ic)) + b(ix,iy,1,ic)*dx(1)*dx(1) + params%dim*stencil(0)*u(ix,iy,1,ic)) / (stencil(0)*params%dim)
                        enddo
                    enddo
                else
                    do iz = is(3), ie(3), dir
                        do iy = is(2), ie(2), dir
                            do ix = is(1), ie(1), dir
                                u(ix,iy,iz,ic) = (-sum(stencil * u(ix-a:ix+a,iy,iz,ic)) - sum(stencil * u(ix,iy-a:iy+a,iz,ic)) - sum(stencil * u(ix,iy,iz-a:iz+a,ic)) + b(ix,iy,iz,ic)*dx(1)*dx(1) + params%dim*stencil(0)*u(ix,iy,iz,ic)) / (stencil(0)*params%dim)
                            enddo
                        enddo
                    enddo
                endif
            else
                if (params%dim == 2) then
                    do iy = is(2), ie(2), dir
                        do ix = is(1), ie(1), dir
                            u(ix,iy,1,ic) = (-sum(stencil_tensor_2D * u(ix-a:ix+a,iy-a:iy+a,1,ic)) + b(ix,iy,1,ic)*dx(1)*dx(1) + stencil_tensor_2D(0,0)*u(ix,iy,1,ic)) / stencil_tensor_2D(0,0)
                        enddo
                    enddo
                else
                    do iz = is(3), ie(3), dir
                        do iy = is(2), ie(2), dir
                            do ix = is(1), ie(1), dir
                                u(ix,iy,iz,ic) = (-sum(stencil_tensor_3D * u(ix-a:ix+a,iy-a:iy+a,iz-a:iz+a,ic)) + b(ix,iy,iz,ic)*dx(1)*dx(1) + stencil_tensor_3D(0,0,0)*u(ix,iy,iz,ic)) / stencil_tensor_3D(0,0,0)
                            enddo
                        enddo
                    enddo
                endif
            endif
        enddo

    end subroutine

    !> \brief compute residual r = b - Ax
    !> \details Compute the residual r = b - Ax, where A is the operator defined by the stencil.
    !> Can also be used to recompute b as b = r + Ax
    subroutine GS_compute_residual(params, u, b, r, dx, recompute_b)
        implicit none

        !> parameter struct
        type (type_params), intent(inout)  :: params
        real(kind=rk), intent(in)         :: u(:, :, :, :)
        real(kind=rk), intent(in)         :: b(:, :, :, :)  !< RHS b
        real(kind=rk), intent(out)        :: r(:, :, :, :)  !< residual r, can be inplace with b
        real(kind=rk), intent(in)         :: dx(:)
        logical, intent(in), optional     :: recompute_b    !< if .true. then b = r + Ax

        integer :: Ax_factor

        integer(kind=ik)                  :: ix,iy,iz,ic, a

        ! factor to multiply Ax with, decides if we compute r = b - Ax or recompute b = r + Ax
        Ax_factor = -1
        if (present(recompute_b)) then
            if (recompute_b) Ax_factor = 1
        endif

        a = params%laplacian_stencil_size

        do ic = 1,size(u,4)
            if (.not. use_tensor) then
                if (params%dim == 2) then
                    do iy = params%g+1, params%g+params%bs(2)
                        do ix = params%g+1, params%g+params%bs(1)
                            r(ix,iy,1,ic) = b(ix,iy,1,ic) + Ax_factor * (sum(stencil * u(ix-a:ix+a,iy,1,ic)) / (dx(1)*dx(1)) \
                                + sum(stencil * u(ix,iy-a:iy+a,1,ic))  / (dx(2)*dx(2)))
                        enddo
                    enddo
                else
                    do iz = params%g+1, params%g+params%bs(3)
                        do iy = params%g+1, params%g+params%bs(2)
                            do ix = params%g+1, params%g+params%bs(1)
                                r(ix,iy,iz,ic) = b(ix,iy,iz,ic) + Ax_factor * (sum(stencil * u(ix-a:ix+a,iy,iz,ic)) / (dx(1)*dx(1)) \
                                    + sum(stencil * u(ix,iy-a:iy+a,iz,ic)) / (dx(2)*dx(2)) + sum(stencil * u(ix,iy,iz-a:iz+a,ic)) / (dx(3)*dx(3)))
                            enddo
                        enddo
                    enddo
                endif
            else
                if (params%dim == 2) then
                    do iy = params%g+1, params%g+params%bs(2)
                        do ix = params%g+1, params%g+params%bs(1)
                            r(ix,iy,1,ic) = b(ix,iy,1,ic) + Ax_factor * sum(stencil_tensor_2D * u(ix-a:ix+a,iy-a:iy+a,1,ic)) / (dx(1)*dx(1))
                        enddo
                    enddo
                else
                    do iz = params%g+1, params%g+params%bs(3)
                        do iy = params%g+1, params%g+params%bs(2)
                            do ix = params%g+1, params%g+params%bs(1)
                                r(ix,iy,iz,ic) = b(ix,iy,iz,ic) + Ax_factor * sum(stencil_tensor_3D * u(ix-a:ix+a,iy-a:iy+a,iz-a:iz+a,ic)) / (dx(1)*dx(1))
                            enddo
                        enddo
                    enddo
                endif
            endif
        enddo
    end subroutine


    ! apply laplacian and compute Ax, used for CG
    subroutine GS_compute_Ax(params, u, ddu, dx, apply_B_RHS)
        implicit none

        !> parameter struct
        type (type_params), intent(inout) :: params
        real(kind=rk), intent(in)         :: u(:, :, :)
        real(kind=rk), intent(inout)      :: ddu(:, :, :)
        real(kind=rk), intent(in)         :: dx(:)
        logical, intent(in), optional     :: apply_B_RHS  !< if .true. then apply the B operator stencil_RHS instead of stencil

        integer(kind=ik)                  :: ix,iy,iz, a, is(1:3), ie(1:3)
        logical :: use_stencil_RHS

        ! get indices over which we loop
        is(:) = 1
        ie(:) = 1
        is(1:params%dim) = params%g + 1
        ie(1:params%dim) = params%g + params%bs(1:params%dim)

        use_stencil_RHS = .false.
        if (present(apply_B_RHS)) use_stencil_RHS = apply_B_RHS

        if (.not. apply_B_RHS) then
            a = params%laplacian_stencil_size
            if (.not. use_tensor) then
                if (params%dim == 2) then
                    do iy = params%g+1, params%g+params%bs(2)
                        do ix = params%g+1, params%g+params%bs(1)
                            ddu(ix,iy,1) =  sum(stencil * u(ix-a:ix+a,iy,1)) / (dx(1)*dx(1)) \
                                + sum(stencil * u(ix,iy-a:iy+a,1)) / (dx(2)*dx(2))
                        enddo
                    enddo
                else
                    do iz = params%g+1, params%g+params%bs(3)
                        do iy = params%g+1, params%g+params%bs(2)
                            do ix = params%g+1, params%g+params%bs(1)
                                ddu(ix,iy,iz) = sum(stencil * u(ix-a:ix+a,iy,iz)) / (dx(1)*dx(1)) \
                                    + sum(stencil * u(ix,iy-a:iy+a,iz)) / (dx(2)*dx(2)) + sum(stencil * u(ix,iy,iz-a:iz+a)) / (dx(3)*dx(3))
                            enddo
                        enddo
                    enddo
                endif
            else
                if (params%dim == 2) then
                    do iy = params%g+1, params%g+params%bs(2)
                        do ix = params%g+1, params%g+params%bs(1)
                            ddu(ix,iy,1) = sum(stencil_tensor_2D * u(ix-a:ix+a,iy-a:iy+a,1)) / (dx(1)*dx(1))
                        enddo
                    enddo
                else
                    do iz = params%g+1, params%g+params%bs(3)
                        do iy = params%g+1, params%g+params%bs(2)
                            do ix = params%g+1, params%g+params%bs(1)
                                ddu(ix,iy,iz) = sum(stencil_tensor_3D * u(ix-a:ix+a,iy-a:iy+a,iz-a:iz+a)) / (dx(1)*dx(1))
                            enddo
                        enddo
                    enddo
                endif
            endif
        else
            if (.not. use_tensor) then
                ddu(:,:,:) = u(:,:,:)
            else
                a = stencil_RHS_size
                ddu(is(1):ie(1),is(2):ie(2),is(3):ie(3)) = 0.0_rk
                if (params%dim == 2) then
                    do iy = -a,a
                        do ix = -a,a
                            if (stencil_RHS_tensor_2D(ix,iy) /= 0.0_rk) then
                                ddu(is(1):ie(1),is(2):ie(2),1) = ddu(is(1):ie(1),is(2):ie(2),1) + stencil_RHS_tensor_2D(ix,iy) * u(is(1)+ix:ie(1)+ix,is(2)+iy:ie(2)+iy,1)
                            endif
                        enddo
                    enddo
                else
                    do iz = -a,a
                        do iy = -a,a
                            do ix = -a,a
                                if (stencil_RHS_tensor_3D(ix,iy,iz) /= 0.0_rk) then
                                    ddu(is(1):ie(1),is(2):ie(2),is(3):ie(3)) = ddu(is(1):ie(1),is(2):ie(2),is(3):ie(3)) + stencil_RHS_tensor_3D(ix,iy,iz) * u(is(1)+ix:ie(1)+ix,is(2)+iy:ie(2)+iy,is(3)+iz:ie(3)+iz)
                                endif
                            enddo
                        enddo
                    enddo
                endif
            endif
        endif
    end subroutine



    subroutine CG_solve_poisson_level0(params, u, f, dx, tol, max_iter)
        use module_params
        implicit none

        ! Arguments
        type (type_params), intent(inout)  :: params
        real(kind=rk), intent(inout) :: u(:, :, :, :)
        real(kind=rk), intent(in) :: f(:, :, :, :)
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
            call GS_compute_Ax(params, u(:,:,:,ic), Ap, dx)
            r(g(1)+1:bs(1)+g(1),g(2)+1:bs(2)+g(2),g(3)+1:bs(3)+g(3)) = f(g(1)+1:bs(1)+g(1),g(2)+1:bs(2)+g(2),g(3)+1:bs(3)+g(3),ic) - Ap(g(1)+1:bs(1)+g(1),g(2)+1:bs(2)+g(2),g(3)+1:bs(3)+g(3))
            p(g(1)+1:bs(1)+g(1),g(2)+1:bs(2)+g(2),g(3)+1:bs(3)+g(3)) = r(g(1)+1:bs(1)+g(1),g(2)+1:bs(2)+g(2),g(3)+1:bs(3)+g(3))
            rsold = sum(r(g(1)+1:bs(1)+g(1),g(2)+1:bs(2)+g(2),g(3)+1:bs(3)+g(3))**2)
            ! write(*, "(A, i4, A, es9.2)") "      CG it ", 0, " res= ", rsold
            do it = 1, max_iter
                ! sync only p to compute Ap
                ! assume that we have only one periodic block here and laplacian is only a + shaped stencil
                call sync_block_level0(params, p, sync_sides_only=.true.)

                call GS_compute_Ax(params, p, Ap, dx)
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