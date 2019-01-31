
!-----------------------------------------------------------------
!> Implementation of zylinder pipe flow
!> \details
!> \version 23.2.2018
!> \author P.Krah
!-----------------------------------------------------------------

module module_pipe_flow

  use module_navier_stokes_params
  use module_precision
  use module_ini_files_parser_mpi
  use module_ns_penalization
  use module_helpers


  implicit none

  !**********************************************************************************************
  ! only this functions are visible outside this module
  PUBLIC :: read_params_pipe_flow,pipe_flow_penalization2D, draw_pipe_sponges
  !**********************************************************************************************
  ! make everything private if not explicitly marked public
  PRIVATE
  !**********************************************************************************************
  real(kind=rk),save  :: p_in, p_out

contains


!=========================================================================================
! INITIALIZATIONs
!=========================================================================================

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> \brief reads parameters for mask function from file
subroutine read_params_pipe_flow( params,FILE )
    implicit none
    !--------------------------------------------
    !> pointer to inifile
    type(inifile) ,intent(inout)       :: FILE
   !> params structure of navier stokes
    type(type_params_ns),intent(inout)  :: params
    !---------------------------------------------
    ! add sponge parameters for the in and outflow of the domain
    call init_simple_sponge(FILE)

    if (params%mpirank==0) then
     write(*,*)
     write(*,*)
     write(*,*) "PARAMS: pipe flow!"
     write(*,'(" ------------------------")')
   endif

   ! geometry of the solid obstacle
    call read_param_mpi(FILE, 'Pipe_flow','p_in', p_in, params%initial_pressure )
    ! free outlet sponge
    call read_param_mpi(FILE, 'Pipe_flow','p_out', p_out, params%initial_pressure )


end subroutine read_params_pipe_flow
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!==========================================================================
subroutine draw_pipe_sponges(mask, x0, dx, Bs, g )
    implicit none
    ! grid
    integer(kind=ik), intent(in) :: g
    integer(kind=ik), dimension(3), intent(in) :: Bs
    !> mask term for every grid point of this block
    real(kind=rk), dimension(:,:,:), intent(out)     :: mask
    !> spacing and origin of block
    real(kind=rk), dimension(2), intent(in)                   :: x0, dx

    call sponge_2D(mask(:,:,1), x0, dx, Bs, g)


end subroutine draw_pipe_sponges
!==========================================================================




!==========================================================================
!> This function adds a penalization term
!> to navier stokes equations
subroutine pipe_flow_penalization2D(Bs, g, x0, dx, mask, phi_ref)
      implicit none
      ! -----------------------------------------------------------------
      integer(kind=ik), intent(in) :: g
      integer(kind=ik), dimension(3), intent(in) :: Bs!< grid parameter
      real(kind=rk), intent(in)     :: x0(2), dx(2)   !< coordinates of block and block spacinf
      real(kind=rk), intent(inout)  :: phi_ref(:,:,:) !< reference values of penalized volume
      real(kind=rk), intent(inout)  :: mask(:,:,:)    !< mask function
      ! -----------------------------------------------------------------
      real(kind=rk) :: z
      real(kind=rk),save,allocatable :: tmp_mask(:,:)    !< mask function
      integer(kind=ik):: ix, ir

      if (.not. allocated(tmp_mask)) allocate(tmp_mask(Bs(1)+2*g,Bs(2)+2*g))
      if (size(mask,1) /= Bs(1)+2*g) call abort(7109,"wrong array size, there's pirates, captain!")
      if (size(phi_ref,1) /= Bs(1)+2*g) call abort(777109,"wrong array size, there's pirates, captain!")
      ! reset mask array
      mask    = 0.0_rk
      phi_ref = 0.0_rk

      ! sponge for in and outflow
      !---------------------------
      call sponge_2D(tmp_mask, x0, dx, Bs, g)

      do ir = g+1, Bs(2) + g
        do ix = g+1, Bs(1) + g
          ! this if is necessary to not overwrite the privious values
          z  = dble(ix-(g+1)) * dx(1) + x0(1)
          if (tmp_mask(ix,ir) > 0.0_rk) then
            ! mask of the inlet and outlet sponge
          !  mask(ix,ir,rhoF) = C_sp_inv*tmp_mask(ix,ir)
          !  mask(ix,ir,UxF ) = C_sp_inv*tmp_mask(ix,ir)
          !  mask(ix,ir,UyF ) = C_sp_inv*tmp_mask(ix,ir)
            mask(ix,ir,pF  ) = C_sp_inv*tmp_mask(ix,ir)
            ! values of the sponge inlet
            !phi_ref(ix,ir,rhoF)= rho0
          !  phi_ref(ix,ir,UxF) = 0.0_rk
          !  phi_ref(ix,ir,UyF) = 0.0_rk
            if (z > params_ns%domain_size(1)*0.5_rk) then
              phi_ref(ix,ir,pF) = p_out
            else
              phi_ref(ix,ir,pF) = p_in
            end if
          end if
        end do
      end do

end subroutine pipe_flow_penalization2D
!==========================================================================


end module module_pipe_flow
