!> \brief main level wrapper for setting non-periodic boundary conditions in ghost patches
subroutine BOUNDCOND_ACM( time, u, g, x0, dx, n_domain, spaghetti_form, edges_only)
    implicit none

    ! it may happen that some source terms have an explicit time-dependency
    ! therefore the general call has to pass time
    real(kind=rk), intent (in) :: time

    ! block data, containg the state vector. In general a 4D field (3 dims+components)
    ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
    real(kind=rk), intent(inout) :: u(1:,1:,1:,1:)

    ! as you are allowed to compute the RHS only in the interior of the field
    ! you also need to know where 'interior' starts: so we pass the number of ghost points
    integer, intent(in) :: g

    ! for each block, you'll need to know where it lies in physical space. The first
    ! non-ghost point has the coordinate x0, from then on its just cartesian with dx spacing
    real(kind=rk), intent(in) :: x0(1:3), dx(1:3)

    ! when implementing boundary conditions, it is necessary to know if the local field (block)
    ! is adjacent to a boundary, because the stencil has to be modified on the domain boundary.
    ! The n_domain tells you if the local field is adjacent to a domain boundary:
    ! n_domain(i) can be either 0, 1, -1,
    !  0: no boundary in the direction +/-e_i
    !  1: boundary in the direction +e_i
    ! -1: boundary in the direction - e_i
    ! currently only acessible in the local stage
    ! NOTE: ACM only supports symmetry BC for the moment (which is handled by wabbit and not ACM)
    integer(kind=2), intent(in) :: n_domain(3)

    ! sometimes we sync data that is decomposed in spaghetti form of SC WC, this can alter the way we set the BCs
    logical , intent(in) :: spaghetti_form

    logical, intent(in), optional  :: edges_only                    !< if true, only set the boundary condition on the edges of the domain, default is false

    real(kind=rk)    :: x, y, z
    integer(kind=ik) :: ix, iy, iz, idir, Bs(3), is, dim
    logical :: disable_antisymmetric, edgesOnly

    edgesOnly = .false.
    if (present(edges_only)) edgesOnly = edges_only

    ! compute the size of blocks
    Bs(1) = size(u,1) - 2*g
    Bs(2) = size(u,2) - 2*g
    Bs(3) = size(u,3) - 2*g
    ! dim from physics module
    dim = params_acm%dim

    if (.not. params_acm%initialized) write(*,'(A)') "WARNING: BOUNDCOND_ACM called but ACM not initialized"

    ! we do not sync only for state vector, if this is the case, then everything is treated as symmetric
    disable_antisymmetric = size(u,4) /= dim + 1 + params_acm%N_scalars

    ! check if this block has any border
    if (maxval(abs(n_domain(1:dim)))==0) return

    ! currently, only symmetric BC is implemented, but here would be the position for more of them
    ! values are mirrored (or anti-symmetric for normal velocities) and boundary values set to 0 or kept
    ! for x=L, the points are in the ghost node layer so we have to approximate them as they do not exist
    ! I central average, computed using taylor expansion, to estimate those values. Central averaging coefficients:
    !    | i    | $\pm5$ | $\pm4$ | $\pm3$ | $\pm2$ | $\pm1$ |
    !    | ---- | ------ | ------ | ------ | ------ | ------ |
    !    | 2nd  |        |        |        |        | 1/2    |
    !    | 4th  |        |        |        | -1/6   | 2/3    |
    !    | 6th  |        |        | 1/20   | -3/10  | 3/4    |
    !    | 8th  |        | -1/70  | 4/35   | -4/10  | 4/5    |
    !    | 10th | 1/252  | -5/126 | 5/28   | -10/21 | 5/6    |
    ! Sometimes we sync for decomposed values, so I just double dx in that case to sync SC and WC correctly

    ! ATTENTION : This is experimental. It works smoothly for 1 symmetry plane in negative direction however there are some problems:
    !   - in positive direction something is yet missing, maybe my ansatz with averaging is just not clever or something with the wavelets
    !   - I did not check for several planes of symmetry - maybe there are some corner problems but I don't think so - to be checked!

    ! ToDo: implement edgesOnly so that we copy less values

    if (n_domain(1) == -1 .and. params_acm%symmetry_BC(1)) then
        do is = 1,size(u,4)
            if (is == 1 .and. .not. disable_antisymmetric) then  ! antisymmetric components
                u(g+1:-1,:,:,is) = 0.0_rk
                u(g:1:-1,:,:,is) = - u(g+2:2*g+1,:,:,is)
            else  ! symmetric components
                ! values on symmetry line are untouched
                u(g:1:-1,:,:,is) = + u(g+2:2*g+1,:,:,is)
            endif
        enddo
    elseif (n_domain(1) == +1 .and. params_acm%symmetry_BC(1)) then
        do is = 1,size(u,4)
            if (is == 1 .and. .not. disable_antisymmetric) then  ! antisymmetric components
                u(Bs(1)+g+1,:,:,is) = 0.0_rk
                u(Bs(1)+g+2:Bs(1)+2*g,:,:,is) = - u(Bs(1)+g:Bs(1)+2:-1,:,:,is)
            else  ! symmetric components
                ! for the approximation of the mirror line, we need to use only SC for the spaghetti form
                if (.not. spaghetti_form) then
                    u(Bs(2)+g+1,:,:,is) = 2.0_rk*(4.0_rk/5.0_rk * u(Bs(2)+g-1,:,:,is) - 4.0_rk/10.0_rk * u(Bs(2)+g-2,:,:,is) + 4.0_rk/35.0_rk * u(Bs(2)+g-3,:,:,is) - 1.0_rk/70.0_rk * u(Bs(2)+g-4,:,:,is))
                else
                    u(Bs(2)+g+1,:,:,is) = 2.0_rk*(4.0_rk/5.0_rk * u(Bs(2)+g-1,:,:,is) - 4.0_rk/10.0_rk * u(Bs(2)+g-3,:,:,is) + 4.0_rk/35.0_rk * u(Bs(2)+g-5,:,:,is) - 1.0_rk/70.0_rk * u(Bs(2)+g-7,:,:,is))
                endif
                u(Bs(1)+g+2:Bs(1)+2*g,:,:,is) = + u(Bs(1)+g:Bs(1)+2:-1,:,:,is)
            endif
        enddo
    elseif (n_domain(2) == -1 .and. params_acm%symmetry_BC(2)) then
        do is = 1,size(u,4)
            if (is == 2 .and. .not. disable_antisymmetric) then  ! antisymmetric components
                u(:,g+1:-1,:,is) = 0.0_rk
                u(:,g:1:-1,:,is) = - u(:,g+2:2*g+1,:,is)
            else  ! symmetric components
                ! values on symmetry line are untouched
                u(:,g:1:-1,:,is) = + u(:,g+2:2*g+1,:,is)
            endif
        enddo
    elseif (n_domain(2) == +1 .and. params_acm%symmetry_BC(2)) then
        do is = 1,size(u,4)
            if (is == 2 .and. .not. disable_antisymmetric) then  ! antisymmetric components
                u(:,Bs(2)+g+1,:,is) = 0.0_rk
                u(:,Bs(2)+g+2:Bs(2)+2*g,:,is) = - u(:,Bs(2)+g:Bs(2)+2:-1,:,is)
            else  ! symmetric components
                if (.not. spaghetti_form) then
                    u(:,Bs(2)+g+1,:,is) = 2.0_rk*(4.0_rk/5.0_rk * u(:,Bs(2)+g-1,:,is) - 4.0_rk/10.0_rk * u(:,Bs(2)+g-2,:,is) + 4.0_rk/35.0_rk * u(:,Bs(2)+g-3,:,is) - 1.0_rk/70.0_rk * u(:,Bs(2)+g-4,:,is))
                else
                    u(:,Bs(2)+g+1,:,is) = 2.0_rk*(4.0_rk/5.0_rk * u(:,Bs(2)+g-1,:,is) - 4.0_rk/10.0_rk * u(:,Bs(2)+g-3,:,is) + 4.0_rk/35.0_rk * u(:,Bs(2)+g-5,:,is) - 1.0_rk/70.0_rk * u(:,Bs(2)+g-7,:,is))
                endif
                u(:,Bs(2)+g+2:Bs(2)+2*g,:,is) = + u(:,Bs(2)+g:Bs(2)+2:-1,:,is)
            endif
        enddo
    endif
    if (dim == 3) then
        if (n_domain(3) == -1 .and. params_acm%symmetry_BC(3)) then
            do is = 1,size(u,4)
                if (is == 3 .and. .not. disable_antisymmetric) then  ! antisymmetric components
                    u(:,:,g+1:-1,is) = 0.0_rk
                    u(:,:,g:1:-1,is) = - u(:,:,g+2:2*g+1,is)
                else  ! symmetric components
                    ! values on symmetry line are untouched
                    u(:,:,g:1:-1,is) = + u(:,:,g+2:2*g+1,is)
                endif
            enddo
        elseif (n_domain(3) == +1 .and. params_acm%symmetry_BC(3)) then
            do is = 1,size(u,4)
                if (is == 3 .and. .not. disable_antisymmetric) then  ! antisymmetric components
                    u(:,:,Bs(3)+g+1,is) = 0.0_rk
                    u(:,:,Bs(3)+g+2:Bs(3)+2*g,is) = - u(:,:,Bs(3)+g:Bs(3)+2:-1,is)
                else  ! symmetric components
                    if (.not. spaghetti_form) then
                        u(:,:,Bs(3)+g+1,is) = 2.0_rk*(4.0_rk/5.0_rk * u(:,:,Bs(3)+g-1,is) - 4.0_rk/10.0_rk * u(:,:,Bs(3)+g-2,is) + 4.0_rk/35.0_rk * u(:,:,Bs(3)+g-3,is) - 1.0_rk/70.0_rk * u(:,:,Bs(3)+g-4,is))
                    else
                        u(:,:,Bs(3)+g+1,is) = 2.0_rk*(4.0_rk/5.0_rk * u(:,:,Bs(3)+g-1,is) - 4.0_rk/10.0_rk * u(:,:,Bs(3)+g-3,is) + 4.0_rk/35.0_rk * u(:,:,Bs(3)+g-5,is) - 1.0_rk/70.0_rk * u(:,:,Bs(3)+g-7,is))
                    endif
                    u(:,:,Bs(3)+g+2:Bs(3)+2*g,is) = + u(:,:,Bs(3)+g:Bs(3)+2:-1,is)
                endif
            enddo
        endif
    endif

    end subroutine BOUNDCOND_ACM