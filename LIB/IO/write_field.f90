! ********************************
! 2D AMR prototype
! --------------------------------
!
! write data of a single 2D field phi
! at timestep n and time t
!
! name: write_field.f90
! date: 30.08.2016
! author: engels
! version: 0.2
!
! ********************************

subroutine write_field(iteration, time, error, dF)

    use module_params
    use module_blocks
    use hdf5_wrapper

    implicit none

    real(kind=rk), intent(in) 	    :: time, error
    integer(kind=ik), intent(in)    :: iteration

    character(len=80)               :: fname, dsetname
    integer(kind=ik)                :: k

    integer(kind=ik), intent(in)    :: dF

    ! create the filename
    write( fname,'("data_",i8.8,".h5")') nint(time * 1.0e4_rk)

    write(*,'(80("*"))')
    write(*,'("Writing data... time=",f15.8," fname=",A)') time, trim(adjustl(fname))
    write(*,'(80("*"))')

    ! overwrite the file, if it already exists
    call init_empty_file( fname )

    ! save block data
    do k = 1, blocks_params%number_max_blocks

        if (blocks(k)%active) then

          ! the name of the block within th hdf5 file
          write(dsetname,'("block_",i8.8)') k

          ! actual writing of block data to file:
          call write_field_hdf5( fname, dsetname, blocks(k)%data_fields(dF)%data_(:,:), .false.)

          ! add useful attributes to the block:
          call write_attribute( fname, dsetname, "treecode", blocks(k)%treecode)
          call write_attribute( fname, dsetname, "time", (/time/))
          call write_attribute( fname, dsetname, "iteration", (/iteration/))
          call write_attribute( fname, dsetname, "coord_x", blocks(k)%coord_x)
          call write_attribute( fname, dsetname, "coord_y", blocks(k)%coord_y)

          call write_attribute( fname, dsetname, "neighbor-id1", (/blocks(k)%neighbor_id(1)/) )
          call write_attribute( fname, dsetname, "neighbor-treecode1", blocks(k)%neighbor_treecode(1,:))

          call write_attribute( fname, dsetname, "neighbor-id2", (/blocks(k)%neighbor_id(2)/) )
          call write_attribute( fname, dsetname, "neighbor-treecode2", blocks(k)%neighbor_treecode(2,:))

          call write_attribute( fname, dsetname, "neighbor-id3", (/blocks(k)%neighbor_id(3)/) )
          call write_attribute( fname, dsetname, "neighbor-treecode3", blocks(k)%neighbor_treecode(3,:))

          call write_attribute( fname, dsetname, "neighbor-id4", (/blocks(k)%neighbor_id(4)/) )
          call write_attribute( fname, dsetname, "neighbor-treecode4", blocks(k)%neighbor_treecode(4,:))

          call write_attribute( fname, dsetname, "neighbor-id5", (/blocks(k)%neighbor_id(5)/) )
          call write_attribute( fname, dsetname, "neighbor-treecode5", blocks(k)%neighbor_treecode(5,:))

          call write_attribute( fname, dsetname, "neighbor-id6", (/blocks(k)%neighbor_id(6)/) )
          call write_attribute( fname, dsetname, "neighbor-treecode6", blocks(k)%neighbor_treecode(6,:))

          call write_attribute( fname, dsetname, "neighbor-id7", (/blocks(k)%neighbor_id(7)/) )
          call write_attribute( fname, dsetname, "neighbor-treecode7", blocks(k)%neighbor_treecode(7,:))

          call write_attribute( fname, dsetname, "neighbor-id8", (/blocks(k)%neighbor_id(8)/) )
          call write_attribute( fname, dsetname, "neighbor-treecode8", blocks(k)%neighbor_treecode(8,:))

          call write_attribute( fname, dsetname, "errors", (/error/))

          call write_attribute( fname, dsetname, "detail", (/blocks(k)%data_fields(dF)%detail/))

        endif

    end do

end subroutine write_field
