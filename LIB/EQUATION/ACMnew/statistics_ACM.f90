!-----------------------------------------------------------------------------
! main level wrapper to compute statistics (such as mean flow, global energy,
! forces, but potentially also derived stuff such as Integral/Kolmogorov scales)
! NOTE: as for the RHS, some terms here depend on the grid as whole, and not just
! on individual blocks. This requires one to use the same staging concept as for the RHS.
!-----------------------------------------------------------------------------
subroutine STATISTICS_ACM( time, u, g, x0, dx, stage, work )
    implicit none

    ! it may happen that some source terms have an explicit time-dependency
    ! therefore the general call has to pass time
    real(kind=rk), intent (in) :: time

    ! block data, containg the state vector. In general a 4D field (3 dims+components)
    ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
    real(kind=rk), intent(inout) :: u(1:,1:,1:,1:)

    ! work data, for mask, vorticity etc. In general a 4D field (3 dims+components)
    ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
    real(kind=rk), intent(inout) :: work(1:,1:,1:,1:)

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
    integer(kind=ik) :: Bs, mpierr, ix, iy, idir
    real(kind=rk) :: tmp(1:6)
    real(kind=rk) :: x, y
    real(kind=rk) :: eps_inv
    ! we have quite some of these work arrays in the code, but they are very small,
    ! only one block. They're ngeligible in front of the lgt_block array.
    real(kind=rk), allocatable, save :: mask(:,:,:), us(:,:,:,:)

    ! compute the size of blocks
    Bs = size(u,1) - 2*g

    if (params_acm%dim==3) then
        if (.not. allocated(mask)) allocate(mask(1:Bs+2*g, 1:Bs+2*g, 1:Bs+2*g))
        if (.not. allocated(us)) allocate(us(1:Bs+2*g, 1:Bs+2*g, 1:Bs+2*g, 1:3))
    else
        if (.not. allocated(mask)) allocate(mask(1:Bs+2*g, 1:Bs+2*g, 1))
        if (.not. allocated(us)) allocate(us(1:Bs+2*g, 1:Bs+2*g, 1, 1:2))
    endif


    select case(stage)
    case ("init_stage")
        !-------------------------------------------------------------------------
        ! 1st stage: init_stage.
        !-------------------------------------------------------------------------
        ! this stage is called only once, NOT for each block.
        ! performs initializations in the RHS module, such as resetting integrals
        params_acm%mean_flow = 0.0_rk
        if (params_acm%forcing_type(1) .eq. "taylor_green") params_acm%error = 0.0_rk
        params_acm%force = 0.0_rk
        params_acm%e_kin = 0.0_rk
        params_acm%enstrophy = 0.0_rk

    case ("integral_stage")
        !-------------------------------------------------------------------------
        ! 2nd stage: integral_stage.
        !-------------------------------------------------------------------------
        ! This stage contains all operations which are running on the blocks
        !
        ! called for each block.

        if (maxval(abs(u))>1.0e5) then
            call abort(6661,"ACM fail: very very large values in state vector.")
        endif

        !-------------------------------------------------------------------------
        ! compute mean flow for output in statistics
        if (params_acm%dim == 2) then
            params_acm%mean_flow(1) = params_acm%mean_flow(1) + sum(u(g+1:Bs+g-1, g+1:Bs+g-1, 1, 1))*dx(1)*dx(2)
            params_acm%mean_flow(2) = params_acm%mean_flow(2) + sum(u(g+1:Bs+g-1, g+1:Bs+g-1, 1, 2))*dx(1)*dx(2)
        else
            params_acm%mean_flow(1) = params_acm%mean_flow(1) + sum(u(g+1:Bs+g-1, g+1:Bs+g-1, g+1:Bs+g-1, 1))*dx(1)*dx(2)*dx(3)
            params_acm%mean_flow(2) = params_acm%mean_flow(2) + sum(u(g+1:Bs+g-1, g+1:Bs+g-1, g+1:Bs+g-1, 2))*dx(1)*dx(2)*dx(3)
            params_acm%mean_flow(3) = params_acm%mean_flow(3) + sum(u(g+1:Bs+g-1, g+1:Bs+g-1, g+1:Bs+g-1, 3))*dx(1)*dx(2)*dx(3)
        endif ! NOTE: MPI_SUM is perfomed in the post_stage.

        !-------------------------------------------------------------------------
        ! if the forcing is taylor-green, then we know the exact solution in time. Therefore
        ! we compute the error w.r.t. this solution heres
        if (params_acm%forcing_type(1) .eq. "taylor_green") then
            do iy = g+1,Bs+g
                do ix = g+1, Bs+g
                    x = x0(1) + dble(ix-g-1)*dx(1)
                    y = x0(2) + dble(iy-g-1)*dx(2)
                    tmp(1) = params_acm%u_mean_set(1) + dsin(x-params_acm%u_mean_set(1)*time)*&
                    dcos(y-params_acm%u_mean_set(2)*time)*dcos(time)
                    tmp(2) = params_acm%u_mean_set(2) - dcos(x-params_acm%u_mean_set(1)*time)*&
                    dsin(y-params_acm%u_mean_set(2)*time)*dcos(time)
                    tmp(3) = 0.25_rk*(dcos(2.0_rk*(x-params_acm%u_mean_set(1)*time)) +&
                    dcos(2.0_rk*(y-params_acm%u_mean_set(2)*time)))*dcos(time)**2
                    params_acm%error(1:3) = params_acm%error(1:3) + abs(u(ix,iy,1,:)-&
                    tmp(1:3))
                    params_acm%error(4:6) = params_acm%error(4:6) + sqrt(tmp(1:3)**2)
                end do
            end do
            params_acm%error = params_acm%error*dx(1)*dx(2)
        end if

        !-------------------------------------------------------------------------
        ! compute fluid force on penalized obstacle. The force can be computed by
        ! volume integration (which is much easier than surface integration), see
        ! Angot et al. 1999
        if (params_acm%dim == 2) then
            call create_mask_2D(time, x0, dx, Bs, g, mask(:,:,1), us(:,:,1,1:2))
        else
            call create_mask_3D(time, x0, dx, Bs, g, mask, us)
        endif

        eps_inv = 1.0_rk / params_acm%C_eta
        mask = mask * eps_inv

        if (params_acm%dim == 2) then
            params_acm%force(1) = params_acm%force(1) + sum( &
            (u(g+1:Bs+g-1, g+1:Bs+g-1, 1, 1)-us(g+1:Bs+g-1, g+1:Bs+g-1, 1, 1))*mask(g+1:Bs+g-1, g+1:Bs+g-1,1))*dx(1)*dx(2)

            params_acm%force(2) = params_acm%force(2) + sum(&
            (u(g+1:Bs+g-1, g+1:Bs+g-1, 1, 2)-us(g+1:Bs+g-1, g+1:Bs+g-1, 1, 2))*mask(g+1:Bs+g-1, g+1:Bs+g-1,1))*dx(1)*dx(2)

            params_acm%force(3) = 0.d0
        else
            do idir = 1, 3
                params_acm%force(idir) = params_acm%force(idir) + sum( &
                ( u(g+1:Bs+g-1, g+1:Bs+g-1, g+1:Bs+g-1, idir) - us(g+1:Bs+g-1, g+1:Bs+g-1, g+1:Bs+g-1, idir))&
                * mask(g+1:Bs+g-1, g+1:Bs+g-1, g+1:Bs+g-1))*dx(1)*dx(2)*dx(3)
            enddo
        endif

        !-------------------------------------------------------------------------
        ! compute kinetic energy in the whole domain (including penalized regions)
        if (params_acm%dim == 2) then
            params_acm%e_kin = params_acm%e_kin + 0.5_rk*sum(u(g+1:Bs+g-1, g+1:Bs+g-1, 1, 1:2)**2)*dx(1)*dx(2)
        else
            params_acm%e_kin = params_acm%e_kin + 0.5_rk*sum(u(g+1:Bs+g-1,g+1:Bs+g-1,g+1:Bs+g-1, 1:3)**2)*dx(1)*dx(2)*dx(3)
        end if

        !-------------------------------------------------------------------------
        ! compute enstrophy in the whole domain (including penalized regions)
        call compute_vorticity(u(:,:,:,1), u(:,:,:,2), work(:,:,:,2), dx, Bs, g, params_acm%discretization, work(:,:,:,:))
        if (params_acm%dim ==2) then
            params_acm%enstrophy = params_acm%enstrophy + sum(work(g+1:Bs+g-1,g+1:Bs+g-1,1,1)**2)*dx(1)*dx(2)
        else
            call abort(6661,"ACM 3D not implemented.")
        end if

    case ("post_stage")
        !-------------------------------------------------------------------------
        ! 3rd stage: post_stage.
        !-------------------------------------------------------------------------
        ! this stage is called only once, NOT for each block.


        !-------------------------------------------------------------------------
        ! mean flow
        tmp(1:3) = params_acm%mean_flow
        call MPI_ALLREDUCE(tmp(1:3), params_acm%mean_flow, 3, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
        if (params_acm%dim == 2) then
            params_acm%mean_flow = params_acm%mean_flow / (params_acm%Lx*params_acm%Ly)
        else
            params_acm%mean_flow = params_acm%mean_flow / (params_acm%Lx*params_acm%Ly*params_acm%Lz)
        endif

        !-------------------------------------------------------------------------
        ! force
        tmp(1:3) = params_acm%force
        call MPI_ALLREDUCE(tmp(1:3), params_acm%force, 3, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)

        !-------------------------------------------------------------------------
        ! kinetic energy
        tmp(1) = params_acm%e_kin
        call MPI_ALLREDUCE(tmp(1), params_acm%e_kin, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)

        !-------------------------------------------------------------------------
        ! kinetic enstrophy
        tmp(1)= params_acm%enstrophy
        call MPI_ALLREDUCE(tmp(1), params_acm%enstrophy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)

        !-------------------------------------------------------------------------
        ! write statistics to ascii files.
        if (params_acm%mpirank == 0) then
            ! write mean flow to disk...
            open(14,file='meanflow.t',status='unknown',position='append')
            write (14,'(4(es15.8,1x))') time, params_acm%mean_flow
            close(14)

            ! write forces to disk...
            open(14,file='forces.t',status='unknown',position='append')
            write (14,'(4(es15.8,1x))') time, params_acm%force
            close(14)

            ! write kinetic energy to disk...
            open(14,file='e_kin.t',status='unknown',position='append')
            write (14,'(2(es15.8,1x))') time, params_acm%e_kin
            close(14)

            ! write enstrophy to disk...
            open(14,file='enstrophy.t',status='unknown',position='append')
            write (14,'(2(es15.8,1x))') time, params_acm%enstrophy
            close(14)
        end if

        if (params_acm%forcing_type(1) .eq. "taylor_green") then
            tmp = params_acm%error
            call MPI_REDUCE(tmp, params_acm%error, 6, MPI_DOUBLE_PRECISION, MPI_SUM, 0, WABBIT_COMM,mpierr)
            !params_acm%error(1:3) = params_acm%error(1:3)/params_acm%error(4:6)
            params_acm%error(1:3) = params_acm%error(1:3)/(params_acm%Lx*params_acm%Ly)

            if (params_acm%mpirank == 0) then
                ! write error to disk...
                open(15,file='error_taylor_green.t',status='unknown',position='append')
                write (15,'(4(es15.8,1x))') time, params_acm%error(1:3)
                close(15)
            end if

        end if
    case default
        call abort(7772,"the STATISTICS wrapper requests a stage this physics module cannot handle.")
    end select


end subroutine STATISTICS_ACM
