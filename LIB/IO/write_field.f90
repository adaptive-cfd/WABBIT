! ********************************
! 2D AMR prototype
! --------------------------------
!
! write data of a single 2D field phi
! at timestep n and time t
!
! name: write_field.f90
! date: 03.08.2016
! author: msr
! version: 0.1
!
! ********************************

subroutine write_field(iteration, time)

    use module_params
    use module_blocks

    implicit none

    real(kind=rk), intent(in) 	    :: time
    integer(kind=ik), intent(in)    :: iteration

    character(len=20)		        :: str_helper
    character(len=265) 		        :: name_file, name_file_xy

    integer(kind=ik)		        :: io_error, i, j, nx, k, N, block_num

    N = size(blocks_params%active_list, dim=1)

    write(*,'("Writing data... time=",f15.8," iteration",i8, " N_active=",i8)') time, iteration, N

    ! save block data
    do k = 1, N

        block_num = blocks_params%active_list(k)

        ! create file name
        write(str_helper, '(i5.5)') iteration
        name_file = trim(params%name_workdir) // trim(params%name_case) // '/' // trim(adjustl(trim(str_helper)))
        write(str_helper, '(f6.2)') time
        name_file = trim(name_file) // '_' // adjustl(trim(str_helper))
        write(str_helper, '(i4.4)') block_num
        name_file = trim(name_file) // '/block_' // adjustl(trim(str_helper))
        write(str_helper, '(i5.5)') iteration
        name_file = trim(name_file) // '_iteration_' // trim(str_helper)
        write(str_helper, '(f6.2)') time

        name_file_xy = trim(name_file) // '_time_' // trim(adjustl(trim(str_helper))) // '.xy'
        name_file = trim(name_file) // '_time_' // trim(adjustl(trim(str_helper))) // '.dat'

        nx = blocks_params%size_block

        ! write coords
        open(unit=99,file=name_file_xy,status='new',action='write', iostat=io_error)

        if (io_error == 0) then
            do i = 1, nx
                write(99, '(f10.4,1x)', advance='no') blocks(block_num)%coord_x(i)
                write(99, '(f10.4,1x)', advance='no') blocks(block_num)%coord_y(i)
                write(99, *)
            end do
        end if
        close(unit=99)

        ! write data
        open(unit=99,file=name_file,status='new',action='write', iostat=io_error)

        if (io_error == 0) then
            do i = 1, nx
                do j = 1, nx
                    write(99, '(f10.4,1x)', advance='no') blocks(block_num)%data1(i,j)
                end do
                write(99, *)
            end do
        end if
        close(unit=99)

    end do

    ! save treecode

    ! create file name
    write(str_helper, '(i5.5)') iteration
    name_file = trim(params%name_workdir) // trim(params%name_case) // '/' // trim(adjustl(trim(str_helper)))
    write(str_helper, '(f6.2)') time
    name_file = trim(name_file) // '_' // adjustl(trim(str_helper))

    name_file = trim(name_file) // '/treecode.dat'

    open(unit=1000,file=name_file,status='new',action='write', iostat=io_error)

    do k = 1, N

        block_num = blocks_params%active_list(k)
        ! write treecode
        if (io_error == 0) then
            write(1000, '(i6,5x)', advance='no') block_num
            do i = 1, blocks(block_num)%level
                write(1000, '(i2,1x)', advance='no') blocks(block_num)%treecode(i)
            end do
            write(1000, *)
        end if

    end do

    close(unit=1000)

end subroutine write_field
