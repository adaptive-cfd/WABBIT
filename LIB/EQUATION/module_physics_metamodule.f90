module module_physics_metamodule

    use module_precision
    ! at this point, you bind all physics modules into one metamodule, so in the rest
    ! of the code, we just load that. as all other physics modules, it provides some
    ! public routines, at which the corresponding actual physics modules are called
    use module_ConvDiff_new
    use module_ACM_new
    use module_navier_stokes_new

    implicit none

    ! I usually find it helpful to use the private keyword by itself initially, which specifies
    ! that everything within the module is private unless explicitly marked public.
    PRIVATE

    !**********************************************************************************************
    ! These are the important routines that are visible to WABBIT:
    !**********************************************************************************************
    PUBLIC :: READ_PARAMETERS, PREPARE_SAVE_DATA, RHS_meta, GET_DT_BLOCK, INICOND_meta, FIELD_NAMES
    !**********************************************************************************************

contains

 !-----------------------------------------------------------------------------
 ! main level wrapper routine to read parameters in the physics module. It reads
 ! from the same ini file as wabbit, and it reads all it has to know. note in physics modules
 ! the parameter struct for wabbit is not available.
 subroutine READ_PARAMETERS( physics, filename )
   implicit none
   character(len=*), intent(in) :: physics
   character(len=*), intent(in) :: filename

   select case ( physics )
   case ('ACM-new')
     call READ_PARAMETERS_ACM( filename )

   case ('ConvDiff-new')
     call READ_PARAMETERS_convdiff( filename )

   case ('navier_stokes')
     call READ_PARAMETERS_NStokes( filename )

   case default
     call abort(1212,'unknown physics...say whaaat?')

   end select

 end subroutine READ_PARAMETERS


 !-----------------------------------------------------------------------------
 ! save data. Since you might want to save derived data, such as the vorticity,
 ! the divergence etc., which are not in your state vector, this routine has to
 ! copy and compute what you want to save to the work array.
 !
 ! In the main code, save_fields than saves the first N_fields_saved components of the
 ! work array to file.
 !
 ! NOTE that as we have way more work arrays than actual state variables (typically
 ! for a RK4 that would be >= 4*dim), you can compute a lot of stuff, if you want to.
 !-----------------------------------------------------------------------------
 subroutine PREPARE_SAVE_DATA( physics, time, u, g, x0, dx, work )
   implicit none
   character(len=*), intent(in) :: physics

   ! it may happen that some source terms have an explicit time-dependency
   ! therefore the general call has to pass time
   real(kind=rk), intent (in) :: time

   ! block data, containg the state vector. In general a 4D field (3 dims+components)
   ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
   real(kind=rk), intent(in) :: u(1:,1:,1:,1:)

   ! as you are allowed to compute the RHS only in the interior of the field
   ! you also need to know where 'interior' starts: so we pass the number of ghost points
   integer, intent(in) :: g

   ! for each block, you'll need to know where it lies in physical space. The first
   ! non-ghost point has the coordinate x0, from then on its just cartesian with dx spacing
   real(kind=rk), intent(in) :: x0(1:3), dx(1:3)

   ! output in work array.
   real(kind=rk), intent(inout) :: work(1:,1:,1:,1:)

   select case(physics)
   case ('ACM-new')
     call PREPARE_SAVE_DATA_ACM( time, u, g, x0, dx, work )

   case ('ConvDiff-new')
     call PREPARE_SAVE_DATA_convdiff( time, u, g, x0, dx, work )

   case ('navier_stokes')
     call PREPARE_SAVE_DATA_NStokes( time, u, g, x0, dx, work )

   case default
     call abort(88119, "[PREPARE_SAVE_DATA (metamodule)] unknown physics....")

   end select

 end subroutine


 !-----------------------------------------------------------------------------
 ! when savig to disk, WABBIT would like to know how you named you variables.
 ! e.g. u(:,:,:,1) is called "ux"
 !
 ! the main routine save_fields has to know how you label the stuff you want to
 ! store from the work array, and this routine returns those strings
 !-----------------------------------------------------------------------------
 subroutine FIELD_NAMES( physics, N, name )
   implicit none
   character(len=*), intent(in) :: physics
   ! component index
   integer(kind=ik), intent(in) :: N
   ! returns the name
   character(len=80), intent(out) :: name

   select case(physics)
   case ('ACM-new')
     call FIELD_NAMES_ACM(N, name)

   case ('ConvDiff-new')
     call FIELD_NAMES_convdiff(N, name)

    case ('navier_stokes')
     call FIELD_NAMES_NStokes(N, name)

   case default
     call abort(88119, "[FIELD_NAMES (metamodule):] unknown physics....")

   end select

 end subroutine FIELD_NAMES


 !-----------------------------------------------------------------------------
 ! main level wrapper to set the right hand side on a block. Note this is completely
 ! independent of the grid any an MPI formalism, neighboring relations and the like.
 ! You just get a block data (e.g. ux, uy, uz, p) and compute the right hand side
 ! from that. Ghost nodes are assumed to be sync'ed.
 !-----------------------------------------------------------------------------
 subroutine RHS_meta( physics, time, u, g, x0, dx, rhs, stage )
   implicit none

   character(len=*), intent(in) :: physics
   ! it may happen that some source terms have an explicit time-dependency
   ! therefore the general call has to pass time
   real(kind=rk), intent (in) :: time

   ! block data, containg the state vector. In general a 4D field (3 dims+components)
   ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
   real(kind=rk), intent(in) :: u(1:,1:,1:,1:)

   ! as you are allowed to compute the RHS only in the interior of the field
   ! you also need to know where 'interior' starts: so we pass the number of ghost points
   integer, intent(in) :: g

   ! for each block, you'll need to know where it lies in physical space. The first
   ! non-ghost point has the coordinate x0, from then on its just cartesian with dx spacing
   real(kind=rk), intent(in) :: x0(1:3), dx(1:3)

   ! output. Note assumed-shape arrays
   real(kind=rk), intent(inout) :: rhs(1:,1:,1:,1:)

   ! stage. there is 3 stages, init_stage, integral_stage and local_stage. If the PDE has
   ! terms that depend on global qtys, such as forces etc, which cannot be computed
   ! from a single block alone, the first stage does that. the second stage can then
   ! use these integral qtys for the actuall RHS evaluation.
   character(len=*), intent(in) :: stage

   select case(physics)
   case ("ACM-new")
     ! this call is not done for all blocks, but only once, globally.
     call RHS_ACM( time, u, g, x0, dx,  rhs, stage )

   case ("ConvDiff-new")
     ! this call is not done for all blocks, but only once, globally.
     call RHS_convdiff( time, u, g, x0, dx, rhs, stage )

   case ("navier_stokes")
     ! this call is not done for all blocks, but only once, globally.
     call RHS_NStokes( time, u, g, x0, dx, rhs, stage )

   case default
     call abort(2152000, "[RHS_wrapper.f90]: physics_type is unknown"//physics)

   end select

 end subroutine RHS_meta

 !-----------------------------------------------------------------------------
 ! subroutine statistics()
 !   implicit none
 ! end subroutine


 !-----------------------------------------------------------------------------
 ! setting the time step is very physics-dependent. Sometimes you have a CFL like
 ! condition, sometimes not. So each physic module must be able to decide on its
 ! time step. This routine is called for all blocks, the smallest returned dt is used.
 !-----------------------------------------------------------------------------
 subroutine GET_DT_BLOCK( physics, time, u, Bs, g, x0, dx, dt )
   implicit none
   character(len=*), intent(in) :: physics
   ! it may happen that some source terms have an explicit time-dependency
   ! therefore the general call has to pass time
   real(kind=rk), intent (in) :: time

   ! block data, containg the state vector. In general a 4D field (3 dims+components)
   ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
   real(kind=rk), intent(in) :: u(1:,1:,1:,1:)

   ! as you are allowed to compute the RHS only in the interior of the field
   ! you also need to know where 'interior' starts: so we pass the number of ghost points
   integer, intent(in) :: g, bs

   ! for each block, you'll need to know where it lies in physical space. The first
   ! non-ghost point has the coordinate x0, from then on its just cartesian with dx spacing
   real(kind=rk), intent(in) :: x0(1:3), dx(1:3)

   ! the dt for this block is returned to the caller:
   real(kind=rk), intent(out) :: dt

   select case (physics)
   case ('ACM-new')
     ! artificial compressibility
     call GET_DT_BLOCK_ACM( time, u, Bs, g, x0, dx, dt )

   case ('ConvDiff-new')
     ! convection-diffusion
     call GET_DT_BLOCK_convdiff( time, u, Bs, g, x0, dx, dt )

   case ('navier_stokes')
     ! navier stokes
     call GET_DT_BLOCK_NStokes( time, u, Bs, g, x0, dx, dt )

   case default
     call abort('phycics module unkown.')

   end select

 end subroutine GET_DT_BLOCK


 !-----------------------------------------------------------------------------
 ! main level wrapper for setting the initial condition on a block
 !-----------------------------------------------------------------------------
 subroutine INICOND_meta( physics, time, u, g, x0, dx )
   implicit none

   character(len=*), intent(in) :: physics
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

   select case (physics)
   case ("ACM-new")
     call INICOND_ACM( time, u, g, x0, dx )

   case ("ConvDiff-new")
     call INICOND_ConvDiff( time, u, g, x0, dx )

   case ("navier_stokes")
     call INICOND_NStokes( time, u, g, x0, dx )

   case default
     call abort(999,"[INICOND (metamodule):] unkown physics. Its getting hard to find qualified personel.")

   end select

 end subroutine INICOND_meta



end module module_physics_metamodule
