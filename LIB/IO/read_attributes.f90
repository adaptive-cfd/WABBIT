!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name read_attributes.f90
!> \version 0.5
!> \author sm
!
!> \brief read attributes saved in a hdf5-file
! = log ======================================================================================
!
!> \date 02/02/18 - create
!

subroutine read_attributes(fname, lgt_n, time, iteration, domain, bs, tc_length, dim, periodic_BC, symmetry_BC, verbosity)

    implicit none
    !> file name
    character(len=*), intent(in)                  :: fname
    !> number of active blocks (required to allocate light data, prior to reading)
    integer(kind=ik), intent(out)                 :: lgt_n
    !> time (to be read from file)
    real(kind=rk), intent(out)                    :: time
    !> iteration (to be read from file)
    integer(kind=ik), intent(out)                 :: iteration
    !> blocksize in the file (required to allocate light data, prior to reading)
    integer(kind=ik), dimension(3), intent(out)   :: Bs
    !> length of treecodes in the file (required to allocate light data, prior to reading)
    integer(kind=ik), intent(out)                 :: tc_length
    !> data dimensionality (2 or 3)
    integer(kind=ik), intent(out)                 :: dim
    !> domain size
    real(kind=rk), dimension(3), intent(out)      :: domain
    logical, intent(in), optional                 :: verbosity !< if verbosity==True generates log output
    logical, intent(inout), optional              :: periodic_BC(1:3), symmetry_BC(1:3)

    integer(kind=ik), dimension(1)                :: iiteration, number_blocks, version
    real(kind=rk), dimension(1)                   :: ttime
    integer(hid_t)                                :: file_id
    integer(kind=ik)                              :: datarank, Nb, rank, ierr, tmp1(1:3), tmp2(1:3)
    integer(kind=hsize_t)                         :: size_field(1:4)
    integer(hsize_t), dimension(2)                :: dims_treecode
    logical :: verbose = .true.

    if (present(verbosity)) verbose=verbosity

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call check_file_exists(fname)

    ! open the file
    call open_file_hdf5( trim(adjustl(fname)), file_id, .false.)


    !----------------------------------------------------------------------------------
    ! version check
    !----------------------------------------------------------------------------------
    ! Files created using the newGhostNodes branch (after 08 April 2020) contain a version number.
    ! if the number is not found, version=0.
    call read_attribute( file_id, "blocks", "version", version)

    if (version(1) == 20200408 .and. rank==0) then
        write(*,*) "--------------------------------------------------------------------"
        write(*,*) "-----WARNING----------WARNING----------WARNING----------WARNING-----"
        write(*,*) "--------------------------------------------------------------------"
        write(*,*) "The file we are trying to read is generated with an intermediate version"
        write(*,*) "of wabbit (after 08 April 2020). In this the file, the grid"
        write(*,*) "definition does not include a redundant point, i.e., a block is defined"
        write(*,*) "with spacing dx = L*2^-J / Bs. This definition was a dead-end, as it lead"
        write(*,*) "to instabilities and other problems. Current versions of WABBIT include a redundant point again."
        write(*,*) ""
        write(*,*) "The newGhostNodes branch still STORED the redundant point for visualization."
        write(*,*) "This was simply the first ghost node. If Bs was odd, this lead to an even number"
        write(*,*) "of points, and this cannot be read with present code versions."
        write(*,*) ""
        write(*,*) "A workaround must be done in preprocessing: upsampling to equidistant resolution and"
        write(*,*) "re-gridding is required. I am truely sorry for this."
        write(*,*) ""
        write(*,*) "If bs was even, the resulting data size is odd, and the file can be read. note however"
        write(*,*) "that the resolution changes slightly, and results cannot be perfectly identical to what would"
        write(*,*) "have been obtained with the newGhostNodes branch"
        write(*,*) "--------------------------------------------------------------------"
        write(*,*) "-----WARNING----------WARNING----------WARNING----------WARNING-----"
        write(*,*) "--------------------------------------------------------------------"
    endif

    if (present(periodic_BC)) then
        call read_attribute(file_id, "blocks", "periodic_BC", tmp1, (/1_ik, 1_ik, 1_ik/))

        periodic_BC = .false.
        where (tmp1 > 0_ik)
            periodic_BC = .true.
        end where
    endif
    if (present(symmetry_BC)) then
        call read_attribute(file_id, "blocks", "symmetry_BC", tmp2, (/0_ik, 0_ik, 0_ik/))

        symmetry_BC = .false.
        where (tmp2 > 0_ik)
            symmetry_BC = .true.
        end where
    endif

    call read_attribute(file_id, "blocks", "domain-size", domain)
    call read_attribute(file_id, "blocks", "time", ttime)
    time = ttime(1)
    call read_attribute(file_id, "blocks", "iteration", iiteration)
    iteration = iiteration(1)
    call read_attribute(file_id, "blocks", "total_number_blocks", number_blocks)
    lgt_n = number_blocks(1)

    !---------------------------------------------------------------------------
    ! Number of blocks and blocksize
    !---------------------------------------------------------------------------
    ! check if we deal with 2D or 3D data
    call get_rank_datafield(file_id, "blocks", datarank)
    if (datarank == 3) then
        ! 2D data
        call get_size_datafield( datarank, file_id, "blocks", size_field(1:datarank))
        Bs(1) = int( size_field(1), kind=ik)
        Bs(2) = int( size_field(2), kind=ik)
        Bs(3) = 1
        Nb = int( size_field(3), kind=ik)
        domain(3) = 0.0_rk
        dim = 2

    elseif (datarank == 4) then
        ! 3D data
        call get_size_datafield( datarank, file_id, "blocks", size_field(1:datarank))
        Bs(1) = int( size_field(1), kind=ik)
        Bs(2) = int( size_field(2), kind=ik)
        Bs(3) = int( size_field(3), kind=ik)
        Nb = int( size_field(4), kind=ik)
        dim = 3

    else
        ! crazy data
        call abort(33321, "Datarank neither 2d nor 3d..that is unusual.")

    endif

    if (modulo(Bs(1),2) == 0) then
        call abort(202009021, "Blocksize Bs(1) is an even number, which this code version cannot handle.")
    endif
    if (modulo(Bs(2),2) == 0) then
        call abort(202009021, "Blocksize Bs(2) is an even number, which this code version cannot handle.")
    endif
    if (modulo(Bs(3),2) == 0) then
        call abort(202009021, "Blocksize Bs(3) is an even number, which this code version cannot handle.")
    endif

    if ( Nb /= lgt_n ) then
        ! the number of blocks stored in metadata and the dimensionality of the
        ! array do not match.
        write(*,*) "Nb= ", Nb, "lgt_n= ", lgt_n
        call abort(333139, "the number of blocks stored in metadata and the dimensionality of the array do not match.")
    endif

    !---------------------------------------------------------------------------
    ! length of treecodes in file
    !---------------------------------------------------------------------------
    ! NOTE: we do store only the treecode, not the level or refinement status
    ! so the length of this array is indeed the treecode length, and not treecode_length+2
    call get_size_datafield(2, file_id, "block_treecode", dims_treecode)
    tc_length = int(dims_treecode(1), kind=ik)



    ! close file and HDF5 library
    call close_file_hdf5(file_id)

    if (rank == 0 .and. verbose) then
        write(*,'(80("~"))')
        write(*,*) "read_attributes.f90: Read important numbers from a file."
        write(*,'(80("~"))')
        write(*,*) "These numbers are used either for initializiation/allocation (in post-"
        write(*,*) "processing) or to check if a file we try to load matches the"
        write(*,*) "specification of the current simulation."
        write(*,*) "We read:"
        write(*,'(" file         = ",A)') trim(adjustl(fname))
        write(*,'(" Bs           = ",3(i3,1x))') Bs
        write(*,'(" domain       = ",3(g15.8,1x))') domain
        write(*,'(" tc_length    = ",i3)') tc_length
        write(*,'(" lgt_n        = ",i8)') lgt_n
        write(*,'(" dim          = ",i8)') dim
        if (present(periodic_BC)) write(*,'(" periodic_BC  = ",3(L1))') periodic_BC
        if (present(symmetry_BC)) write(*,'(" symmetry_BC  = ",3(L1))') symmetry_BC
        write(*,'(80("~"))')
    endif

end subroutine read_attributes
