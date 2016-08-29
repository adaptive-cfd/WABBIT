module module_interpolation
  use module_params
  use module_blocks

  implicit none

contains

  include "prediction_2D.f90"
  include "restriction_2D.f90"

end module
