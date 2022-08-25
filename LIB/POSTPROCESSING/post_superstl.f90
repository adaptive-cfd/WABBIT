subroutine post_superstl(params)
    use module_precision
    use module_mesh
    use module_params
    use module_mpi
    use module_operators
    use module_physics_metamodule
    use module_time_step
    use module_stl_file_reader
    use module_helpers
    use module_ini_files_parser_mpi
    use module_t_files

    implicit none

    type (type_params), intent(inout)  :: params
    real(kind=4), allocatable, dimension(:,:) :: triangles, normals
    real(kind=rk), allocatable, dimension(:,:) :: superstl
    real(kind=rk) :: scale, origin(3)
    character(len=cshort) :: fname_stl, fname_out, dummy
    integer :: i, ntri

    call get_command_argument(2, dummy)

    ! does the user need help?
    if (dummy=='--help' .or. dummy=='--h' .or. dummy=='-h') then
        if (params%rank==0) then
            write(*,*) "------------------------------------------------------------------"
            write(*,*) "./wabbit-post --superstl --x0 2.0,3.0,4.0 --o vertexnormals.sstl --stl input.stl"
            write(*,*) "------------------------------------------------------------------"
            write(*,*) ""
            write(*,*) ""
            write(*,*) ""
            write(*,*) ""
            write(*,*) "------------------------------------------------------------------"
        end if
        return
    endif

    ! defaults
    scale = 1.0_rk
    origin = 0.0_rk
    N_MAX_COMPONENTS = 1

    ! fetch parameters from command line call
    do i = 1, COMMAND_ARGUMENT_COUNT()
        call get_command_argument(i,dummy)

        select case (dummy)
        case ("--stl")
            call get_command_argument(i+1, fname_stl)
            call check_file_exists( fname_stl )
        case ("-o")
            call get_command_argument(i+1, fname_out)
        case ("--scale")
            call get_command_argument(i+1, dummy)
            read(dummy,*) scale
        case ("--x0")
            call get_command_argument(i+1, dummy)
            read(dummy,*) origin
        end select
    enddo

    ! if (params%number_procs>1) call abort(2610191,"routine not yet parallelized")

    if (params%rank==0) then
        write(*,'("super-STL file is ",A)') trim(adjustl(fname_stl))
        write(*,'("STL scaling factor is scale=",g12.3)') scale
        write(*,'("origin shift=",3(g12.3,1x))') origin
    endif

    call read_stl_file(fname_stl, ntri, triangles, normals)
    call normalize_stl_file( ntri, triangles, "--scale", scale, origin )
    ! compute all vertex- and edgenormals, store them in an ascii file
    ! See:
    ! Baerentzen, Aanaes. Generating Signed Distance Fields From Triangle Meshes. (IMM-TECHNICAL REPORT-2002-21)
    call compute_superstl_file( triangles, normals, ntri, superstl )

    if (params%rank == 0) then
        write(*,*) "Writing output to: "//trim(adjustl(fname_out))

        open (14,file = fname_out, status='replace')
        do i = 1, size(superstl, 1)
            write(14,'(30(es15.8,1x))') superstl(i,:)
        enddo
        close(14)
    endif
end subroutine
