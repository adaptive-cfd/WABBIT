!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name init_physics_modules.f90
!> \version 0.5
!> \author engels
!
!> \brief Wrapper to call the parameter-initialization routines for the physics module in use.
!!
!! NOTE: this is *NOT* a part of module_params.f90 in order to avoid circulars in the makefile
!!
!! input:
!!           - params struct to know what module is used.
!!           - filename of the ini file to be read
!! output:
!!           - parameters are stored in the modules themselves.
!!
!! = log ======================================================================================
!! \n
!! 17/12/17 - create
!
!**********************************************************************************************

subroutine init_physics_modules( params, filename )
  ! of course, we have to load all available physics modules here.
  use module_ACM_new
  use module_convdiff_new
  use module_convection_diffusion
  use module_navier_stokes
  ! NOTE: this is *NOT* a part of module_params.f90 in order to avoid circulars in the makefile
  ! therefore load also params module.
  use module_params

  implicit none

  type (type_params), intent(in) :: params
  character(len=*), intent(in) :: filename

  ! call the initialization routines for the physics module that is in use
  select case ( params%physics_type )
  case ('ACM-new')
    call READ_PARAMETERS_ACM( filename )

  case ('ConvDiff-new')
    call READ_PARAMETERS_convdiff( filename )

  case default
    call abort(1212,'unknown physics...say whaaat?')

  end select

end subroutine
