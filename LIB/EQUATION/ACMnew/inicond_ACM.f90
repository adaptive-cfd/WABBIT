!-----------------------------------------------------------------------------
! main level wrapper for setting the initial condition on a block
!-----------------------------------------------------------------------------
subroutine INICOND_ACM( time, u, g, x0, dx, work, adapting )
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

  ! if we are still adapting the initial condition, we may use penalization for refinement.
  ! if the initial grid is adapted we set our initial condition without penalization (impulsive start).
  logical, intent(in) :: adapting

  real(kind=rk)    :: x,y
  integer(kind=ik) :: Bs, ix, iy

  ! compute the size of blocks
  Bs = size(u,1) - 2*g

  select case (params_acm%inicond)
  case("meanflow")
    u = 0.0_rk
    u(:,:,:,1) = params_acm%u_mean_set(1)
    u(:,:,:,2) = params_acm%u_mean_set(2)
    if (params_acm%dim == 3) then
      u(:,:,:,3) = params_acm%u_mean_set(3)
    endif
  case("taylor_green")
    do iy= 1,Bs+2*g
      do ix= 1, Bs+2*g
        x = x0(1) + dble(ix-g-1)*dx(1)
        y = x0(2) + dble(iy-g-1)*dx(2)
        call continue_periodic(x,params_acm%Lx)
        call continue_periodic(y,params_acm%Ly)
        u(ix,iy,1,1) = params_acm%u_mean_set(1) + dsin(x)*dcos(y)
        u(ix,iy,1,2) = params_acm%u_mean_set(2) - dcos(x)*dsin(y)
        u(ix,iy,1,3) = 0.25_rk*(dcos(2.0_rk*x) + dcos(2.0_rk*y))
      end do
    end do
  case default
    write(*,*) "errorrroororor"
  end select
  ! if we use volume penalization, the mask is first used for refinement of the grid.
  ! In a second stage, the initial condition without penalization is then applied to the refined grid.
  if (adapting .and. params_acm%penalization) then
      call create_mask_2D(work(:,:,1,1), x0, dx, Bs, g )
      u(:,:,:,1) = (1.0_rk-work(:,:,:,1))*u(:,:,:,1)
      u(:,:,:,2) = (1.0_rk-work(:,:,:,1))*u(:,:,:,2)
      if (params_acm%dim == 3) then
          u(:,:,:,3) = (1.0_rk-work(:,:,:,1))*u(:,:,:,3)
      end if
  end if

end subroutine INICOND_ACM
