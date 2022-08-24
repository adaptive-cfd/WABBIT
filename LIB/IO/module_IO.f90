module module_IO
    use mpi
    use hdf5
    use module_params
    use module_timing
    use module_hdf5_wrapper
    use module_mpi
    ! use mesh module, since we want to compute origin/spacing of blocks
    use module_mesh
    ! use module operators for computation of the vorticity field
    use module_operators, only: compute_vorticity
    ! use physics modules to save the data
    use module_physics_metamodule
    ! for reading parameters and storing them in each hdf file
    use module_ini_files_parser_mpi
    use module_helpers, only : check_file_exists, block_contains_NaN
    use module_treelib
    use module_forestMetaData

    implicit none

contains

#include "saveHDF5_tree.f90"
#include "readHDF5vct_tree.f90"
#include "read_attributes.f90"
#include "read_file_flusi.f90"
#include "forest_IO.f90"


    ! returns the number of blocks stored in an HDF5 file
    function getNumberBlocksH5File(fname)
        character(len=*), intent(in) :: fname
        integer(kind=ik) :: getNumberBlocksH5File

        integer(hid_t) :: file_id

        ! open the file
        call open_file_hdf5( trim(adjustl(fname)), file_id, .false.)

        ! read number of blocks...
        call read_attribute(file_id, "blocks", "total_number_blocks", getNumberBlocksH5File)

        ! close file and HDF5 library
        call close_file_hdf5(file_id)
    end function

end module module_IO
