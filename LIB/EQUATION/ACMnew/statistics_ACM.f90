!-----------------------------------------------------------------------------
! main level wrapper to compute statistics (such as mean flow, global energy,
! forces, but potentially also derived stuff such as Integral/Kolmogorov scales)
! NOTE: as for the RHS, some terms here depend on the grid as whole, and not just
! on individual blocks. This requires one to use the same staging concept as for the RHS.
!-----------------------------------------------------------------------------
subroutine STATISTICS_ACM( time, dt, u, g, x0, dx, stage, work, grid_qty )
    implicit none

    ! it may happen that some source terms have an explicit time-dependency
    ! therefore the general call has to pass time
    real(kind=rk), intent (in) :: time, dt

    ! block data, containg the state vector. In general a 4D field (3 dims+components)
    ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
    real(kind=rk), intent(inout) :: u(1:,1:,1:,1:)

    ! work data, for mask, vorticity etc. In general a 4D field (3 dims+components)
    ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
    real(kind=rk), intent(inout) :: work(1:,1:,1:,1:)
    real(kind=rk), intent(inout) :: grid_qty(1:,1:,1:,1:)

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
    integer(kind=ik) :: Bs, mpierr, ix, iy, iz
    real(kind=rk) :: tmp(1:6), tmp_meanflow(1:3), tmp_force(1:3), tmp_residual(1:3), tmp_ekin, tmp_volume
    real(kind=rk) :: x, y, CFL, CFL_eta, CFL_nu
    real(kind=rk) :: eps_inv, dV, dx_min
    real(kind=rk), save :: umag
    ! we have quite some of these work arrays in the code, but they are very small,
    ! only one block. They're ngeligible in front of the lgt_block array.
    real(kind=rk), allocatable, save :: mask(:,:,:), us(:,:,:,:), div(:,:,:)
    integer(kind=2), allocatable, save :: mask_color(:,:,:)

    ! compute the size of blocks
    Bs = size(u,1) - 2*g

    if (params_acm%dim==3) then
        if (.not. allocated(mask_color)) allocate(mask_color(1:Bs+2*g, 1:Bs+2*g, 1:Bs+2*g))
        if (.not. allocated(mask)) allocate(mask(1:Bs+2*g, 1:Bs+2*g, 1:Bs+2*g))
        if (.not. allocated(div)) allocate(div(1:Bs+2*g, 1:Bs+2*g, 1:Bs+2*g))
        if (.not. allocated(us)) allocate(us(1:Bs+2*g, 1:Bs+2*g, 1:Bs+2*g, 1:3))
        dV = dx(1)*dx(2)*dx(3)
    else
        if (.not. allocated(mask_color)) allocate(mask_color(1:Bs+2*g, 1:Bs+2*g, 1))
        if (.not. allocated(mask)) allocate(mask(1:Bs+2*g, 1:Bs+2*g, 1))
        if (.not. allocated(div)) allocate(div(1:Bs+2*g, 1:Bs+2*g, 1))
        if (.not. allocated(us)) allocate(us(1:Bs+2*g, 1:Bs+2*g, 1, 1:2))
        dV = dx(1)*dx(2)
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
        params_acm%mask_volume = 0.0_rk
        params_acm%u_residual = 0.0_rk
        params_acm%div_max = 0.0_rk
        params_acm%div_min = 0.0_rk
        umag = 0.0_rk

        if (params_acm%geometry == "Insect") call Update_Insect(time, Insect)

    case ("integral_stage")
        !-------------------------------------------------------------------------
        ! 2nd stage: integral_stage.
        !-------------------------------------------------------------------------
        ! This stage contains all operations which are running on the blocks
        !
        ! called for each block.

        if (maxval(abs(u))>1.0e5) then
            call abort(16102018, "ACM fail: very very large values in state vector.")
        endif

        !-------------------------------------------------------------------------
        ! if the forcing is taylor-green, then we know the exact solution in time. Therefore
        ! we compute the error w.r.t. this solution here
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

        ! tmp values for computing the current block only
        tmp_meanflow = 0.0_rk
        tmp_force = 0.0_rk
        tmp_residual = 0.0_rk
        tmp_ekin = 0.0_rk
        tmp_volume = 0.0_rk

        if (params_acm%dim == 2) then
            ! --- 2D --- --- 2D --- --- 2D --- --- 2D --- --- 2D --- --- 2D ---
            call create_mask_2D( time, x0, dx, Bs, g, mask(:,:,1), us(:,:,1,1:2) )
            eps_inv = 1.0_rk / params_acm%C_eta

            ! note in 2D case, uz is ignored, so we pass p just for fun.
            call divergence( u(:,:,:,1), u(:,:,:,2), u(:,:,:,3), dx, Bs, g, params_acm%discretization, div)

            ! mask divergence inside the solid body
            where (mask>0.0_rk)
                div = 0.00_rk
            end where

            do iy = g+1, Bs+g-1 ! Note: loops skip redundant points
            do ix = g+1, Bs+g-1
                ! compute mean flow for output in statistics
                tmp_meanflow(1) = tmp_meanflow(1) + u(ix,iy,1,1)
                tmp_meanflow(2) = tmp_meanflow(2) + u(ix,iy,1,2)

                ! volume of mask (useful to see if it is properly generated)
                tmp_volume = tmp_volume + mask(ix,iy,1)

                ! for the penalization term, we need to divide by C_eta
                mask(ix,iy,1) = mask(ix,iy,1) * eps_inv

                ! forces acting on body
                tmp_force(1) = tmp_force(1) + (u(ix,iy,1,1)-us(ix,iy,1,1))*mask(ix,iy,1)
                tmp_force(2) = tmp_force(2) + (u(ix,iy,1,2)-us(ix,iy,1,2))*mask(ix,iy,1)

                ! residual velocity in the solid domain
                tmp_residual(1) = max( tmp_residual(1), (u(ix,iy,1,1)-us(ix,iy,1,1)) * mask(ix,iy,1))
                tmp_residual(2) = max( tmp_residual(2), (u(ix,iy,1,2)-us(ix,iy,1,2)) * mask(ix,iy,1))

                ! kinetic energy
                tmp_ekin = tmp_ekin + 0.5_rk*sum( u(ix,iy,1,1:2)**2 )

                ! maximum of velocity in the field
                umag = max( umag, u(ix,iy,1,1)*u(ix,iy,1,1) + u(ix,iy,1,2)*u(ix,iy,1,2) )

                ! maximum/min divergence in velocity field
                params_acm%div_max = max( params_acm%div_max, div(ix,iy,1) )
                params_acm%div_min = min( params_acm%div_min, div(ix,iy,1) )
            enddo
            enddo
        else
            ! --- 3D --- --- 3D --- --- 3D --- --- 3D --- --- 3D --- --- 3D ---
            call create_mask_3D( time, x0, dx, Bs, g, mask, mask_color, us, grid_qty=grid_qty )
            eps_inv = 1.0_rk / params_acm%C_eta

            ! compute divergence on this block
            call divergence( u(:,:,:,1), u(:,:,:,2), u(:,:,:,3), dx, Bs, g, params_acm%discretization, div)

            ! mask divergence inside the solid body
            where (mask>0.0_rk)
                div = 0.00_rk
            end where

            do iz = g+1, Bs+g-1 ! Note: loops skip redundant points
            do iy = g+1, Bs+g-1
            do ix = g+1, Bs+g-1
                ! compute mean flow for output in statistics
                tmp_meanflow(1) = tmp_meanflow(1) + u(ix,iy,iz,1)
                tmp_meanflow(2) = tmp_meanflow(2) + u(ix,iy,iz,2)
                tmp_meanflow(3) = tmp_meanflow(3) + u(ix,iy,iz,3)

                ! volume of mask (useful to see if it is properly generated)
                tmp_volume = tmp_volume + mask(ix, iy, iz)

                ! for the penalization term, we need to divide by C_eta
                mask(ix,iy,iz) = mask(ix,iy,iz) * eps_inv

                ! forces acting on body
                tmp_force(1) = tmp_force(1) + (u(ix,iy,iz,1) - us(ix,iy,iz,1)) * mask(ix,iy,iz)
                tmp_force(2) = tmp_force(2) + (u(ix,iy,iz,2) - us(ix,iy,iz,2)) * mask(ix,iy,iz)
                tmp_force(3) = tmp_force(3) + (u(ix,iy,iz,3) - us(ix,iy,iz,3)) * mask(ix,iy,iz)

                ! residual velocity in the solid domain
                tmp_residual(1) = max( tmp_residual(1), (u(ix,iy,iz,1)-us(ix,iy,iz,1))*mask(ix,iy,iz) )
                tmp_residual(2) = max( tmp_residual(2), (u(ix,iy,iz,2)-us(ix,iy,iz,2))*mask(ix,iy,iz) )
                tmp_residual(3) = max( tmp_residual(3), (u(ix,iy,iz,3)-us(ix,iy,iz,3))*mask(ix,iy,iz) )

                tmp_ekin = tmp_ekin + 0.5_rk*sum( u(ix,iy,iz,1:3)**2 )

                ! maximum of velocity in the field
                umag = max( umag, u(ix,iy,iz,1)*u(ix,iy,iz,1) + u(ix,iy,iz,2)*u(ix,iy,iz,2) + u(ix,iy,iz,3)*u(ix,iy,iz,3) )

                ! maximum/min divergence in velocity field
                params_acm%div_max = max( params_acm%div_max, div(ix,iy,iz) )
                params_acm%div_min = min( params_acm%div_min, div(ix,iy,iz) )
            enddo
            enddo
            enddo
        endif

        ! we just computed the values on the current block, which we now add to the
        ! existing blocks in the variables (recall normalization by dV)
        params_acm%u_residual = params_acm%u_residual + tmp_residual * dV
        params_acm%mean_flow = params_acm%mean_flow + tmp_meanflow * dV
        params_acm%mask_volume = params_acm%mask_volume + tmp_volume * dV
        params_acm%force = params_acm%force + tmp_force * dV
        params_acm%e_kin = params_acm%e_kin + tmp_ekin * dV

        !-------------------------------------------------------------------------
        ! compute enstrophy in the whole domain (including penalized regions)
        call compute_vorticity(u(:,:,:,1), u(:,:,:,2), work(:,:,:,2), dx, Bs, g, params_acm%discretization, work(:,:,:,:))

        if (params_acm%dim ==2) then
            params_acm%enstrophy = params_acm%enstrophy + sum(work(g+1:Bs+g-1,g+1:Bs+g-1,1,1)**2)*dx(1)*dx(2)
        else
            params_acm%enstrophy = 0.0_rk
            ! call abort(6661,"ACM 3D not implemented.")
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
            params_acm%mean_flow = params_acm%mean_flow / (params_acm%domain_size(1)*params_acm%domain_size(2))
        else
            params_acm%mean_flow = params_acm%mean_flow / (params_acm%domain_size(1)*params_acm%domain_size(2)*params_acm%domain_size(3))
        endif

        !-------------------------------------------------------------------------
        ! force
        tmp(1:3) = params_acm%force
        call MPI_ALLREDUCE(tmp(1:3), params_acm%force, 3, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)

        !-------------------------------------------------------------------------
        ! residual velocity in solid domain
        tmp(1:3) = params_acm%u_residual
        call MPI_ALLREDUCE(tmp(1:3), params_acm%u_residual, 3, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)

        !-------------------------------------------------------------------------
        ! volume of mask (useful to see if it is properly generated)
        tmp(1) = params_acm%mask_volume
        call MPI_ALLREDUCE(tmp(1), params_acm%mask_volume, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)

        !-------------------------------------------------------------------------
        ! kinetic energy
        tmp(1) = params_acm%e_kin
        call MPI_ALLREDUCE(tmp(1), params_acm%e_kin, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)

        !-------------------------------------------------------------------------
        ! divergence
        tmp(1) = params_acm%div_min
        call MPI_ALLREDUCE(tmp(1), params_acm%div_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, WABBIT_COMM, mpierr)

        tmp(1) = params_acm%div_max
        call MPI_ALLREDUCE(tmp(1), params_acm%div_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)

        !-------------------------------------------------------------------------
        ! kinetic enstrophy
        tmp(1)= params_acm%enstrophy
        call MPI_ALLREDUCE(tmp(1), params_acm%enstrophy, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)

        tmp(1) = umag
        call MPI_ALLREDUCE(tmp(1), umag, 1, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)

        !-------------------------------------------------------------------------
        ! write statistics to ascii files.
        if (params_acm%mpirank == 0) then
            open(14,file='umag.t',status='unknown',position='append')
            write(14,'(5(es15.8,1x))') time, sqrt(umag), params_acm%c_0, &
            params_acm%c_0/sqrt(umag), sqrt(umag) + sqrt(params_acm%c_0**2 + umag)
            close(14)

            dx_min = 2.0_rk**(-params_acm%Jmax) * params_acm%domain_size(1) / real(params_acm%Bs-1, kind=rk)
            CFL   = dt * (sqrt(umag) + sqrt(params_acm%c_0**2 + umag)) / dx_min
            CFL_nu = dt * params_acm%nu / dx_min**2
            CFL_eta = dt / params_acm%C_eta

            open(14,file='CFL.t',status='unknown',position='append')
            write(14,'(4(es15.8,1x))') time, CFL, CFL_nu, CFL_eta
            close(14)

            ! write mean flow to disk...
            open(14,file='meanflow.t',status='unknown',position='append')
            write(14,'(4(es15.8,1x))') time, params_acm%mean_flow
            close(14)

            ! write divergence to disk...
            open(14,file='div.t',status='unknown',position='append')
            write(14,'(3(es15.8,1x))') time, params_acm%div_max, params_acm%div_min
            close(14)

            ! write forces to disk...
            open(14,file='forces.t',status='unknown',position='append')
            write(14,'(4(es15.8,1x))') time, params_acm%force
            close(14)

            ! write kinetic energy to disk...
            open(14,file='e_kin.t',status='unknown',position='append')
            write(14,'(2(es15.8,1x))') time, params_acm%e_kin
            close(14)

            ! write enstrophy to disk...
            open(14,file='enstrophy.t',status='unknown',position='append')
            write(14,'(2(es15.8,1x))') time, params_acm%enstrophy
            close(14)

            ! write mask_volume to disk...
            open(14,file='mask_volume.t',status='unknown',position='append')
            write(14,'(2(es15.8,1x))') time, params_acm%mask_volume
            close(14)

            ! write residual velocity to disk...
            open(14,file='u_residual.t',status='unknown',position='append')
            write(14,'(4(es15.8,1x))') time, params_acm%u_residual
            close(14)
        end if

        if (params_acm%forcing_type(1) .eq. "taylor_green") then
            tmp = params_acm%error
            call MPI_REDUCE(tmp, params_acm%error, 6, MPI_DOUBLE_PRECISION, MPI_SUM, 0, WABBIT_COMM,mpierr)
            !params_acm%error(1:3) = params_acm%error(1:3)/params_acm%error(4:6)
            params_acm%error(1:3) = params_acm%error(1:3)/(params_acm%domain_size(1)*params_acm%domain_size(2))

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
