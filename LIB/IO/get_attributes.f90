!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name get_attributes.f90
!> \version 0.5
!> \author sm
!
!> \brief read attributes saved in a hdf5-file
!
!>
!! input:
!!           - parameter array
!!           - name of the file we want to read from
!!
!! output:
!!           - time, iteration and domain size
!!           - number of active blocks (light and heavy)
!!
!!
!! = log ======================================================================================
!! \n
!! 02/02/18 - create
!

subroutine get_attributes(fname, lgt_n, time, iteration, domain)

    implicit none
    !> file name
    character(len=*), intent(in)                  :: fname
    !> number of active blocks
    integer(kind=ik), intent(out)                 :: lgt_n
    !> time (to be read from file)
    real(kind=rk), intent(out)                    :: time
    !> iteration (to be read from file)
    integer(kind=ik), intent(out)                 :: iteration
    !> domain size
    real(kind=rk), dimension(3), intent(out)      :: domain
    integer(kind=ik), dimension(1)                :: iiteration, number_blocks
    real(kind=rk), dimension(1)                   :: ttime
    integer(hid_t)                                :: file_id


    call check_file_exists(fname)
    ! open the file
    call open_file_hdf5( trim(adjustl(fname)), file_id, .false.)
    ! read attributes
    call read_attribute(file_id, "blocks", "domain-size", domain)
    call read_attribute(file_id, "blocks", "time", ttime)
    time = ttime(1)
    call read_attribute(file_id, "blocks", "iteration", iiteration)
    iteration = iiteration(1)
    call read_attribute(file_id, "blocks", "total_number_blocks", number_blocks)
    lgt_n = number_blocks(1)

    ! close file and HDF5 library
    call close_file_hdf5(file_id)
end subroutine get_attributes
