!> \brief Wrapper to call the parameter-initialization routines for the physics module in use.
!! NOTE: this is *NOT* a part of module_params.f90 in order to avoid circulars in the makefile
!! input:    - params struct to know what module is used.
!!           - filename of the ini file to be read
!! output:   - parameters are stored in the modules themselves.
!**********************************************************************************************

subroutine init_physics_modules( params, filename, N_mask_components)
  ! of course, we have to load all available physics modules here.
  use module_physics_metamodule
  ! NOTE: this is *NOT* a part of module_params.f90 in order to avoid circulars in the makefile
  ! therefore load also params module.
  use module_params

  implicit none

  type (type_params), intent(in) :: params
  ! number of grid-dependent (and not time-dependend qtys) is decided by the physics modules
  integer(kind=ik), intent(out) :: N_mask_components
  character(len=*), intent(in) :: filename

  if (params%rank==0) then
    write(*,'(80("─"))')
    write(*,*) "Initializing physics modules"
    write(*,'(80("─"))')
  endif

  N_mask_components = 0

  ! call the initialization routines for the physics module that is in use.
  ! WABBIT decides how many ghost nodes we have (because the versions >=2024 determine G 
  ! automatically depending on the wavelet). Therefore, we pass the number.
  call READ_PARAMETERS_meta( params%physics_type, filename, N_mask_components, params%g )
  

end subroutine
