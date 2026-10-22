!> \file
!> \brief Poisson-related operations for WABBIT
!! \details
!! This file contains routines for various operations related to the Poisson equation:
!! - Helmholtz projection (velocity correction)
!! - Computing pressure from velocity
!! - Computing velocity from vorticity
!! etc.

!-------------------------------------------------------------------------------
!> \brief Helmholtz projection: correct velocity field to be divergence-free
!> \details
!! Performs Helmholtz projection: solve Laplacian(phi) = div(u), then u := u - grad(phi)
!-------------------------------------------------------------------------------
subroutine helmholtz_projection(params, time, hvy_block, hvy_tmp, tree_ID, diagnostics, force_convergence)
    implicit none

    type(type_params), intent(inout)   :: params
    real(kind=rk), intent(in)          :: time
    real(kind=rk), intent(inout)       :: hvy_block(:, :, :, :, :)
    real(kind=rk), intent(inout)       :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik), intent(in)       :: tree_ID
    logical, intent(in), optional      :: diagnostics
    logical, intent(in), optional      :: force_convergence  !< enforce stricter convergence by using fixed iterations

    integer(kind=ik) :: k, hvy_id, lgt_id, Bs(1:3), g, g_RHS
    real(kind=rk) :: x0(1:3), dx(1:3), div_max, div_rms, div_sum
    integer(kind=ik) :: mpicode
    logical :: do_diagnostics, do_force_convergence
    character(len=cshort) :: saved_cycle_end_criteria
    integer(kind=ik) :: saved_cycle_it

    Bs = params%Bs
    g  = params%g
    g_RHS = params%g_RHS
    do_diagnostics = .false.
    if (present(diagnostics)) do_diagnostics = diagnostics
    
    do_force_convergence = .false.
    if (present(force_convergence)) do_force_convergence = force_convergence
    
    ! Temporarily enforce stricter convergence if requested
    if (do_force_convergence) then
        saved_cycle_end_criteria = params%poisson_cycle_end_criteria
        saved_cycle_it = params%poisson_cycle_it
        params%poisson_cycle_end_criteria = "fixed_iterations"
        params%poisson_cycle_it = max(15, params%poisson_cycle_it)
    endif

    !---------------------------------------------------------------------------
    ! Step 1: Compute divergence -> hvy_tmp(:,:,:,1:,:)
    !---------------------------------------------------------------------------
    call sync_ghosts_tree(params, hvy_block(:,:,:,1:params%dim,:), tree_ID)
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)
        call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
        call get_block_spacing_origin(params, lgt_id, x0, dx)
        call compute_divergence(hvy_block(:,:,:,1:params%dim,hvy_id), dx, Bs, g, &
                                params%poisson_order, hvy_tmp(:,:,:,1,hvy_id))
    enddo

    if (do_diagnostics) then
        div_max = 0.0_rk
        div_sum = 0.0_rk
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k, tree_ID)
            call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
            call get_block_spacing_origin(params, lgt_id, x0, dx)
            div_max = max(div_max, maxval(abs(hvy_tmp(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 1, hvy_id))))
            div_sum = div_sum + sum(hvy_tmp(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 1, hvy_id)**2) * product(dx(1:params%dim))
        enddo
        call MPI_ALLREDUCE(MPI_IN_PLACE, div_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpicode)
        call MPI_ALLREDUCE(MPI_IN_PLACE, div_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpicode)
        div_rms = sqrt(div_sum / product(params%domain_size(1:params%dim)))
        if (params%rank == 0) write(*,'("Helmholtz: before projection max|div(u)|=", es12.4, " RMS=", es12.4)') div_max, div_rms
    endif

    !---------------------------------------------------------------------------
    ! Step 2: Solve Laplacian(phi) = div(u)  ->  hvy_tmp(:,:,:,1,:)
    ! Multigrid solver does synching before and after solve
    !---------------------------------------------------------------------------
    call multigrid_solve(params, hvy_tmp(:,:,:,2:2,:), hvy_tmp(:,:,:,1:1,:), hvy_tmp(:,:,:,3:3,:), tree_ID, init_0=.true., hvy_full=hvy_tmp, verbose=.false., time=time)

    !---------------------------------------------------------------------------
    ! Step 3: Correct velocity: u := u - grad(phi)
    !---------------------------------------------------------------------------
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)
        call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
        call get_block_spacing_origin(params, lgt_id, x0, dx)
        call subtract_scalar_gradient(hvy_block(:,:,:,1:params%dim,hvy_id), dx, Bs, g, &
                                      params%poisson_order, hvy_tmp(:,:,:,2,hvy_id))
    enddo

    if (do_diagnostics) then
        call sync_ghosts_tree(params, hvy_block(:,:,:,1:params%dim,:), tree_ID)
        div_max = 0.0_rk
        div_sum = 0.0_rk
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k, tree_ID)
            call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
            call get_block_spacing_origin(params, lgt_id, x0, dx)
            call compute_divergence(hvy_block(:,:,:,1:params%dim,hvy_id), dx, Bs, g, &
                                    params%poisson_order, hvy_tmp(:,:,:,1,hvy_id))
            div_max = max(div_max, maxval(abs(hvy_tmp(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 1, hvy_id))))
            div_sum = div_sum + sum(hvy_tmp(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 1, hvy_id)**2) * product(dx(1:params%dim))
        enddo
        call MPI_ALLREDUCE(MPI_IN_PLACE, div_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpicode)
        call MPI_ALLREDUCE(MPI_IN_PLACE, div_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpicode)
        div_rms = sqrt(div_sum / product(params%domain_size(1:params%dim)))
        if (params%rank == 0) write(*,'("Helmholtz: after projection  max|div(u)|=", es12.4, " RMS=", es12.4)') div_max, div_rms
    endif
    
    ! Restore original convergence settings
    if (do_force_convergence) then
        params%poisson_cycle_end_criteria = saved_cycle_end_criteria
        params%poisson_cycle_it = saved_cycle_it
    endif

end subroutine helmholtz_projection


!-------------------------------------------------------------------------------
!> \brief Compute pressure from velocity field
!> \details
!! This routine computes the pressure field from a given velocity field by:
!!   1. Computing the velocity RHS (nonlinear + viscous terms)
!!   2. Computing divergence of the RHS: div(dU/dt)
!!   3. Solving the pressure Poisson equation: Laplacian(p) = div(dU/dt)
!!   4. Storing pressure in hvy_block(:,:,:,dim+1,:)
!!
!! This is useful for post-processing and saving pressure alongside velocity.
!-------------------------------------------------------------------------------
subroutine pressure_from_velocity(params, time, hvy_block, hvy_tmp, hvy_mask, tree_ID, force_convergence)
    implicit none

    type(type_params), intent(inout)   :: params
    real(kind=rk), intent(in)          :: time
    real(kind=rk), intent(inout)       :: hvy_block(:, :, :, :, :)
    real(kind=rk), intent(inout)       :: hvy_tmp(:, :, :, :, :)
    real(kind=rk), intent(in)          :: hvy_mask(:, :, :, :, :)
    integer(kind=ik), intent(in)       :: tree_ID
    logical, intent(in), optional      :: force_convergence  !< enforce stricter convergence by using fixed iterations

    integer(kind=ik) :: k, hvy_id, lgt_id, Bs(1:3), g, g_RHS
    real(kind=rk) :: x0(1:3), dx(1:3)
    integer(kind=2), dimension(3) :: n_domain
    character(len=cshort) :: format_string, saved_cycle_end_criteria
    integer(kind=ik) :: saved_cycle_it
    logical :: do_force_convergence

    Bs = params%Bs
    g  = params%g
    g_RHS = params%g_RHS  ! could optimize synching, but I'm not sure if g_RHS is always correct for post functions
    
    do_force_convergence = .false.
    if (present(force_convergence)) do_force_convergence = force_convergence
    
    ! Temporarily enforce stricter convergence if requested
    if (do_force_convergence) then
        saved_cycle_end_criteria = params%poisson_cycle_end_criteria
        saved_cycle_it = params%poisson_cycle_it
        params%poisson_cycle_end_criteria = "fixed_iterations"
        params%poisson_cycle_it = max(15, params%poisson_cycle_it)
    endif

    !---------------------------------------------------------------------------
    ! Step 1: Compute velocity RHS (nonlinear + viscous terms, no pressure)
    ! Use the local_stage from the physics module
    !---------------------------------------------------------------------------
    call sync_ghosts_tree(params, hvy_block(:,:,:,1:params%dim,:), tree_ID)
    
    ! Zero out pressure field to avoid including pressure gradient in RHS
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)
        hvy_block(:,:,:,params%dim+1,hvy_id) = 0.0_rk
    enddo

    ! compute RHS of momentum equation, needed for RHS of pressure-poisson equation
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)
        call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
        call get_block_spacing_origin(params, lgt_id, x0, dx)

        if ( .not. All(params%periodic_BC) ) then
            ! check if block is adjacent to a boundary of the domain, if this is the case we use one sided stencils
            call get_adjacent_boundary_surface_normal( params, lgt_id, n_domain )
        endif
        
        ! Call the physics RHS with "local_stage" to compute dU/dt without pressure gradient
        ! we'd ideally just love to compute the velocity equations of the RHS for ACM or NSPP but ACM will always compute pressure RHS as well - we'll simply ignore it
        call RHS_meta(params%physics_type, time, hvy_block(:,:,:,1:params%n_eqn_rhs,hvy_id), g, x0, dx, &
                 hvy_tmp(:,:,:,1:params%n_eqn_rhs,hvy_id), hvy_mask(:,:,:,:,hvy_id), &
                 "local_stage", n_domain, params%poisson_order)
    enddo

    ! ! debugging output
    ! call saveHDF5_tree("RHSx_0.h5", 0.0_rk, 0, 1, params, hvy_tmp, tree_ID)
    ! call saveHDF5_tree("RHSy_0.h5", 0.0_rk, 0, 2, params, hvy_tmp, tree_ID)

    !---------------------------------------------------------------------------
    ! Step 2: Compute divergence of velocity RHS -> hvy_tmp(:,:,:,params%dim+1,:)
    !---------------------------------------------------------------------------
    call sync_ghosts_tree(params, hvy_tmp(:,:,:,1:params%dim,:), tree_ID)
    
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)
        call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
        call get_block_spacing_origin(params, lgt_id, x0, dx)
        call compute_divergence(hvy_tmp(:,:,:,1:params%dim,hvy_id), dx, Bs, g, &
                                params%poisson_order, hvy_tmp(:,:,:,params%dim+1,hvy_id))
    enddo

    ! ! debugging output
    ! call saveHDF5_tree("divRHS_0.h5", 0.0_rk, 0, params%dim+1, params, hvy_tmp, tree_ID)

    !---------------------------------------------------------------------------
    ! Step 3: Solve Laplacian(p) = div(dU/dt)  ->  hvy_tmp(:,:,:,1,:)
    ! Multigrid solver does synching before and after solve
    !---------------------------------------------------------------------------
    call multigrid_solve(params, hvy_tmp(:,:,:,1:1,:), hvy_tmp(:,:,:,params%dim+1:params%dim+1,:), &
                         hvy_tmp(:,:,:,2:2,:), tree_ID, &
                         init_0=.true., verbose=.false., hvy_full=hvy_tmp, time=time)

    !---------------------------------------------------------------------------
    ! Step 4: Copy pressure to state vector hvy_block(:,:,:,dim+1,:)
    !---------------------------------------------------------------------------
    do k = 1, hvy_n(tree_ID)
        hvy_id = hvy_active(k, tree_ID)
        if (params%dim == 2) then
            ! 2D case
            hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, 1, params%dim+1, hvy_id) = &
            hvy_tmp(g+1:Bs(1)+g, g+1:Bs(2)+g, 1, 1, hvy_id)
        else
            ! 3D case
            hvy_block(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, params%dim+1, hvy_id) = &
            hvy_tmp(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 1, hvy_id)
        end if
    enddo
    
    ! Restore original convergence settings
    if (do_force_convergence) then
        params%poisson_cycle_end_criteria = saved_cycle_end_criteria
        params%poisson_cycle_it = saved_cycle_it
    endif

end subroutine pressure_from_velocity


!-------------------------------------------------------------------------------
!> \brief Reconstruct velocity from vorticity field
!> \details
!! This routine reconstructs a divergence-free velocity field from a given vorticity:
!!
!! 2D: Given omega (scalar), solve for stream function psi:
!!     Laplacian(psi) = -omega
!!     Then: u = dpsi/dy, v = -dpsi/dx
!!
!! 3D: Given omega_vector, solve for vector potential A using ∇·A = 0 gauge:
!!     Laplacian(A) = -omega  (3 independent Poisson equations)
!!     Then: u = ∇ × A
!!
!! Input:  vorticity in hvy_block(:,:,:,1,:) for 2D, hvy_block(:,:,:,1:3,:) for 3D
!! Output: velocity in hvy_block(:,:,:,1:2,:) for 2D, hvy_block(:,:,:,1:3,:) for 3D
!-------------------------------------------------------------------------------
subroutine velocity_from_vorticity(params, hvy_block, hvy_tmp, tree_ID, force_convergence)
    implicit none

    type(type_params), intent(inout)   :: params
    real(kind=rk), intent(inout)       :: hvy_block(:, :, :, :, :)
    real(kind=rk), intent(inout)       :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik), intent(in)       :: tree_ID
    logical, intent(in), optional      :: force_convergence  !< enforce stricter convergence by using fixed iterations

    integer(kind=ik) :: k, hvy_id, lgt_id, Bs(1:3), g, n_hvy_tmp, icomp
    real(kind=rk) :: x0(1:3), dx(1:3)
    character(len=cshort) :: saved_cycle_end_criteria
    integer(kind=ik) :: saved_cycle_it
    logical :: do_force_convergence

    Bs = params%Bs
    g  = params%g
    n_hvy_tmp = size(hvy_tmp, 4)

    do_force_convergence = .false.
    if (present(force_convergence)) do_force_convergence = force_convergence
    
    ! Temporarily enforce stricter convergence if requested
    if (do_force_convergence) then
        saved_cycle_end_criteria = params%poisson_cycle_end_criteria
        saved_cycle_it = params%poisson_cycle_it
        params%poisson_cycle_end_criteria = "fixed_iterations"
        params%poisson_cycle_it = max(15, params%poisson_cycle_it)
    endif

    if (params%dim == 2) then
        !-----------------------------------------------------------------------
        ! 2D case: solve for stream function psi from vorticity omega
        !-----------------------------------------------------------------------
        ! Input: omega in hvy_block(:,:,:,1,:)
        ! Step 1: Setup RHS = -omega in hvy_block(:,:,:,1,:) (already negative)
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k, tree_ID)
            hvy_block(:,:,:,1,hvy_id) = -hvy_block(:,:,:,1,hvy_id)
        enddo
        
        ! Step 2: Solve Laplacian(psi) = -omega
        ! Solution goes to hvy_block(:,:,:,2,:), work array in hvy_block(:,:,:,3,:)
        call multigrid_solve(params, hvy_block(:,:,:,3:3,:), hvy_block(:,:,:,1:1,:), &
                             hvy_block(:,:,:,2:2,:), tree_ID, init_0=.true., verbose=.false., hvy_full=hvy_block)
        
        ! Step 3: Compute velocity from stream function: u = dpsi/dy, v = -dpsi/dx
        call sync_ghosts_tree(params, hvy_block(:,:,:,3:3,:), tree_ID)
        do k = 1, hvy_n(tree_ID)
            hvy_id = hvy_active(k, tree_ID)
            call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
            call get_block_spacing_origin(params, lgt_id, x0, dx)
            
            ! Compute gradient of psi
            call compute_gradient(hvy_block(:,:,:,3,hvy_id), dx, Bs, g, &
                                 params%poisson_order, hvy_block(:,:,:,1:2,hvy_id))
            
            ! Rotate: (dpsi/dx, dpsi/dy) -> (u, v) = (dpsi/dy, -dpsi/dx)
            ! Store temporarily in component 3
            hvy_block(:,:,:,3,hvy_id) = hvy_block(:,:,:,1,hvy_id)  ! save dpsi/dx
            hvy_block(:,:,:,1,hvy_id) = hvy_block(:,:,:,2,hvy_id)   ! u = dpsi/dy
            hvy_block(:,:,:,2,hvy_id) = -hvy_block(:,:,:,3,hvy_id)  ! v = -dpsi/dx
        enddo

    else
        !-----------------------------------------------------------------------
        ! 3D case: solve for vector potential A, then compute u = curl(A)
        !-----------------------------------------------------------------------
        ! Input: (omega_x, omega_y, omega_z) in hvy_block(:,:,:,1:3,:)
        ! Solve: Laplacian(A) = -omega for all 3 components
        ! Then: u = ∇ × A
        
        if (n_hvy_tmp >= 9) then
            !-------------------------------------------------------------------
            ! Option A: Sufficient space - use hvy_tmp(:,:,:,1:9,:)
            ! Slots 1:3 for solution A, slots 4:6 for RHS, slots 7:9 for work
            !-------------------------------------------------------------------
            ! Step 1: Copy RHS = -omega to hvy_tmp(:,:,:,4:6,:)
            do k = 1, hvy_n(tree_ID)
                hvy_id = hvy_active(k, tree_ID)
                hvy_tmp(:,:,:,4:6,hvy_id) = -hvy_block(:,:,:,1:3,hvy_id)
            enddo
            
            ! Step 2: Solve 3 Poisson equations: Laplacian(A_i) = -omega_i
            call multigrid_solve(params, hvy_tmp(:,:,:,1:3,:), hvy_tmp(:,:,:,4:6,:), hvy_tmp(:,:,:,7:9,:), tree_ID, &
                                init_0=.true., verbose=.false., hvy_full=hvy_tmp)
            ! Vector potential A is now in hvy_tmp(:,:,:,1:3,:)
            
            ! Step 3: Compute velocity u = curl(A)
            call sync_ghosts_tree(params, hvy_tmp(:,:,:,1:3,:), tree_ID)
            do k = 1, hvy_n(tree_ID)
                hvy_id = hvy_active(k, tree_ID)
                call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
                call get_block_spacing_origin(params, lgt_id, x0, dx)
                
                ! Compute curl of vector potential A
                call compute_vorticity(hvy_tmp(:,:,:,1:3,hvy_id), dx, Bs, g, &
                                      params%poisson_order, hvy_block(:,:,:,1:3,hvy_id))
            enddo
            
        elseif (n_hvy_tmp >= 5) then
            !-------------------------------------------------------------------
            ! Option B: Limited space - solve sequentially using hvy_block
            ! hvy_block has dim+1=4 components, use slots 1-3 for solving
            ! Need to preserve input vorticity: backup components 2-3 to hvy_tmp(5:6)
            ! and store solutions A in hvy_tmp(1:3)
            !-------------------------------------------------------------------
            ! Backup vorticity components 2-3 (component 1 used first, doesn't need backup)
            do k = 1, hvy_n(tree_ID)
                hvy_id = hvy_active(k, tree_ID)
                hvy_tmp(:,:,:,4:5,hvy_id) = hvy_block(:,:,:,2:3,hvy_id)
            enddo
            
            do icomp = 1, 3
                ! Step 1: Setup RHS = -omega_icomp in hvy_block(:,:,:,1,:)
                do k = 1, hvy_n(tree_ID)
                    hvy_id = hvy_active(k, tree_ID)
                    if (icomp == 1) then
                        hvy_block(:,:,:,1,hvy_id) = -hvy_block(:,:,:,1,hvy_id)
                    elseif (icomp == 2) then
                        hvy_block(:,:,:,1,hvy_id) = -hvy_tmp(:,:,:,4,hvy_id)  ! from backup
                    else ! icomp == 3
                        hvy_block(:,:,:,1,hvy_id) = -hvy_tmp(:,:,:,5,hvy_id)  ! from backup
                    endif
                enddo
                
                ! Step 2: Solve Laplacian(A_icomp) = -omega_icomp
                ! Solution in hvy_block(:,:,:,2,:), work in hvy_block(:,:,:,3,:)
                call multigrid_solve(params, hvy_block(:,:,:,2:2,:), hvy_block(:,:,:,1:1,:), &
                                   hvy_block(:,:,:,3:3,:), tree_ID, init_0=.true., verbose=.false., hvy_full=hvy_block)
                
                ! Step 3: Store A_icomp in hvy_tmp for later curl computation
                do k = 1, hvy_n(tree_ID)
                    hvy_id = hvy_active(k, tree_ID)
                    hvy_tmp(:,:,:,icomp,hvy_id) = hvy_block(:,:,:,2,hvy_id)
                enddo
            enddo
            
            ! Step 4: Compute velocity u = curl(A) where A is now in hvy_tmp(:,:,:,1:3,:)
            call sync_ghosts_tree(params, hvy_tmp(:,:,:,1:3,:), tree_ID)
            do k = 1, hvy_n(tree_ID)
                hvy_id = hvy_active(k, tree_ID)
                call hvy2lgt(lgt_id, hvy_id, params%rank, params%number_blocks)
                call get_block_spacing_origin(params, lgt_id, x0, dx)
                
                ! Compute curl of vector potential A
                call compute_vorticity(hvy_tmp(:,:,:,1:3,hvy_id), dx, Bs, g, &
                                      params%poisson_order, hvy_block(:,:,:,1:3,hvy_id))
            enddo
        else
            call abort(260130, "velocity_from_vorticity: need at least 5 slots in hvy_tmp for 3D")
        endif
    endif    
    ! Restore original convergence settings
    if (do_force_convergence) then
        params%poisson_cycle_end_criteria = saved_cycle_end_criteria
        params%poisson_cycle_it = saved_cycle_it
    endif
end subroutine velocity_from_vorticity
