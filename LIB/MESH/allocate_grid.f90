subroutine allocate_tree(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, hvy_active, &
    lgt_sortednumlist, hvy_work, hvy_tmp, hvy_gridQ)

    implicit none

    !> user defined parameter structure
    type (type_params), intent(inout)                   :: params
    !> light data array
    integer(kind=ik), allocatable, intent(out)          :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), allocatable, intent(out)             :: hvy_block(:, :, :, :, :)
    !> heavy temp data: used for saving, filtering etc (work array)
    real(kind=rk), allocatable, optional, intent(out)   :: hvy_tmp(:, :, :, :, :)
    !> grid qtys (qtys that do depend only on the grid, i.e. mask funciton and geometry factors)
    real(kind=rk), allocatable, optional, intent(out)   :: hvy_gridQ(:, :, :, :, :)
    !> heavy work array: used for RHS evaluation in multistep methods (like RK4: 00, k1, k2 etc)
    real(kind=rk), allocatable, optional, intent(out)   :: hvy_work(:, :, :, :, :, :)
    !> neighbor array (heavy data)
    integer(kind=ik), allocatable, intent(out)          :: hvy_neighbor(:,:)
    !> list of active blocks (light data)
    integer(kind=ik), allocatable, intent(out)          :: lgt_active(:)
    !> list of active blocks (light data)
    integer(kind=ik), allocatable, intent(out)          :: hvy_active(:)
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), allocatable, intent(out)       :: lgt_sortednumlist(:,:)
    ! local shortcuts:
    integer(kind=ik)                                    :: g, Neqn, number_blocks,&
    rank, number_procs, dim
    integer(kind=ik), dimension(3)                      :: Bs
    integer(kind=ik)    :: rk_steps
    real(kind=rk)       :: effective_memory
    integer             :: status, nrhs_slots, nwork, nx, ny, nz, max_neighbors

    real(kind=rk)      :: maxmem, mem_per_block
    real(kind=rk), parameter ::  nstages = 2.0_rk ! stages for ghost node synching
    character(len=80)  :: memstring
    integer(kind=ik)   :: i



    ! set parameters for readability
    rank            = params%rank
    Bs              = params%Bs
    g               = params%n_ghosts
    Neqn            = params%n_eqn
    number_procs    = params%number_procs
    nx = Bs(1)+2*g
    ny = Bs(2)+2*g
    nz = Bs(3)+2*g

    if (params%dim==2) then
        nz = 1
        max_neighbors = 16
    else
        max_neighbors   = 74
    endif

    ! 19 oct 2018: The work array hvy_work is modified to be used in "register-form"
    ! that means one rhs is stored in a 5D subset of a 6D array.
    ! Hence, nrhs_slots is number of slots for RHS saving:
    if (params%time_step_method == "RungeKuttaGeneric") then
        nrhs_slots = size(params%butcher_tableau,1)
    elseif (params%time_step_method == "Krylov") then
        nrhs_slots = params%M_krylov +3
    elseif ((params%time_step_method == 'none').or.(params%time_step_method == 'no')) then
        nrhs_slots = 0
    else
        call abort(191018161, "time_step_method is unkown: "//trim(adjustl(params%time_step_method)))
    endif

    nwork = max( 2*Neqn, params%N_fields_saved)

    if (rank == 0) then
        write(*,'(80("_"))')
        write(*,'(A)') "INIT: Beginning memory allocation and initialization."
        write(*,'(A)') "ALLOCATING A TREE"
        write(*,'(80("_"))')
    endif

    !Automatic memory management. If specified --memory=0.3GB in the call line,
    if (params%number_blocks < 1) then
        !---------------------------------------------------------------------------
        ! Automatic memory management. If specified --memory=0.3GB in the call line,
        ! wabbit will automatically select the number of blocks per rank to be allocated
        ! to consume this amount of memory. helpful on local machines.
        !---------------------------------------------------------------------------
        ! loop over all arguments until you find the string "--memory" in them (or not)
        ! this ensures compatibility with future versions, as the argument may be anywhere in the call.
        do i = 1, command_argument_count()
            call get_command_argument(i, memstring)
            ! is memory limitation used?
            if ( index(memstring,"--memory=")==1 ) then
                if (params%rank==0) write(*,'("INIT: automatic selection of blocks per rank is active!")')

                ! read memory from command line (in GB)
                read(memstring(10:len_trim(memstring)-2),* ) maxmem

                if (params%rank==0) write(*,'("INIT: total memory: ",f9.4,"GB")') maxmem

                ! memory per MPIRANK (in GB)
                maxmem = maxmem / dble(params%number_procs)

                if (params%rank==0) write(*,'("INIT: memory-per-rank: ",f9.4,"GB")') maxmem

                if (params%dim==3) then
                    mem_per_block = real(Neqn) * real((Bs(1)+2*g)*(Bs(2)+2*g)*(Bs(3)+2*g)) & ! hvy_block
                    + 2.0 * nstages * real(Neqn) * real((Bs(1)+2*g)*(Bs(2)+2*g)*(Bs(3)+2*g) - ((Bs(1))*(Bs(2))*(Bs(3)))) &  ! real buffer ghosts
                    + 2.0 * nstages * real(max_neighbors) * 5 / 2.0 ! int bufer (4byte hence /2)

                    ! hvy_gridQ
                    if ( present(hvy_gridQ) ) then
                        mem_per_block = mem_per_block + real(params%n_gridQ) * real((Bs(1)+2*g)*(Bs(2)+2*g)*(Bs(3)+2*g))
                    endif

                    ! hvy_tmp
                    if ( present(hvy_tmp) ) then
                        mem_per_block = mem_per_block + real(max( 2*Neqn, params%N_fields_saved)) * real((Bs(1)+2*g)*(Bs(2)+2*g)*(Bs(3)+2*g))
                    endif

                    ! hvy_work
                    if ( present(hvy_work) ) then
                        mem_per_block = mem_per_block + real(Neqn) * real(nrhs_slots) * real((Bs(1)+2*g)*(Bs(2)+2*g)*(Bs(3)+2*g))
                    endif


                else
                    mem_per_block = real(Neqn) * real((Bs(1)+2*g)*(Bs(2)+2*g)) & ! hvy_block
                    + 2.0 * nstages * real(Neqn) * real((Bs(1)+2*g)*(Bs(2)+2*g) - (Bs(1)*Bs(2))) &  ! real buffer ghosts
                    + 2.0 * nstages * real(max_neighbors) * 5 / 2.0 ! int bufer (4byte hence /2)

                    ! hvy_gridQ
                    if ( present(hvy_gridQ) ) then
                        mem_per_block = mem_per_block + real(params%n_gridQ) * real((Bs(1)+2*g)*(Bs(2)+2*g))
                    endif

                    ! hvy_tmp
                    if ( present(hvy_tmp) ) then
                        mem_per_block = mem_per_block + real(max(2*Neqn, params%N_fields_saved)) * real((Bs(1)+2*g)*(Bs(2)+2*g))
                    endif

                    ! hvy_work
                    if ( present(hvy_work) ) then
                        mem_per_block = mem_per_block + real(Neqn) * real(nrhs_slots) * real((Bs(1)+2*g)*(Bs(2)+2*g))
                    endif

                endif

                ! in GB:
                mem_per_block = mem_per_block * 8.0e-9

                params%number_blocks = nint( maxmem / mem_per_block)

                if (params%rank==0) then
                    write(*,'("INIT: for the desired memory we can allocate ",i8," blocks per rank")') params%number_blocks
                    write(*,'("INIT: we allocated ",i8," blocks per rank (total: ",i8," blocks) ")') params%number_blocks, &
                    params%number_blocks*params%number_procs
                endif

                if ((params%adapt_mesh .or. params%adapt_inicond) .and. (params%number_blocks<2**params%dim)) then
                    call abort(1909181740,"[ini_file_to_params.f90]: The number of blocks for the --memory specification is lower than 2^d&
                    & and comp is adaptive: we cannot fetch all 2^d blocks on one CPU for merging. Use more memory or less cores.")
                endif
            endif
        end do
    endif


    if (rank == 0) then
        write(*,'("INIT: mpisize=",i6)') params%number_procs
        write(*,'("INIT: nwork=",i6)') nwork
        write(*,'("INIT: Bs(1)=",i7," Bs(2)=",i7," Bs(3)=",i7," blocks-per-rank=",i7," total blocks=", i7)') &
        Bs(1),Bs(2),Bs(3), params%number_blocks, params%number_blocks*params%number_procs
        write(*,'("INIT: Allocating a ",i1,"D case.")') params%dim
    endif


    !---------------------------------------------------------------------------
    allocate( hvy_block( nx, ny, nz, Neqn, params%number_blocks ) )
    if (rank==0) then
        write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
        "hvy_block", product(real(shape(hvy_block)))*8.0e-9, shape(hvy_block)
    endif

    !---------------------------------------------------------------------------
    ! work data (Runge-Kutta substeps and old time level)
    if (present(hvy_work)) then
        allocate( hvy_work( nx, ny, nz, Neqn, params%number_blocks, nrhs_slots ) )
        if (rank==0) then
            write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
            "hvy_work", product(real(shape(hvy_work)))*8.0e-9, shape(hvy_work)
        endif
    end if

    if ( present(hvy_tmp) ) then
        allocate( hvy_tmp( nx, ny, nz, nwork, params%number_blocks )  )
        if (rank==0) then
            write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
            "hvy_tmp", product(real(shape(hvy_tmp)))*8.0e-9, shape(hvy_tmp)
        endif
    endif

    if ( present(hvy_gridQ) .and. params%n_gridQ > 0 ) then
        allocate( hvy_gridQ( nx, ny, nz, params%n_gridQ, params%number_blocks )  )
        if (rank==0) then
            write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
            "hvy_gridQ", product(real(shape(hvy_gridQ)))*8.0e-9, shape(hvy_gridQ)
        endif
    endif

    !---------------------------------------------------------------------------
    allocate( hvy_neighbor( params%number_blocks, max_neighbors ) )
    if (rank==0) then
        write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
        "hvy_neighbor", product(real(shape(hvy_neighbor)))*8.0e-9, shape(hvy_neighbor)
    endif

    !---------------------------------------------------------------------------)
    allocate( lgt_block( number_procs*params%number_blocks, params%max_treelevel+EXTRA_LGT_FIELDS) )
    if (rank==0) then
        write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
        "lgt_block", product(real(shape(lgt_block)))*4.0e-9, shape(lgt_block)
    endif

    !---------------------------------------------------------------------------
    allocate( lgt_sortednumlist( size(lgt_block,1), 2) )
    if (rank==0) then
        write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
        "lgt_sortednumlist", product(real(shape(lgt_sortednumlist)))*4.0e-9, shape(lgt_sortednumlist)
    endif

    !---------------------------------------------------------------------------
    allocate( lgt_active( size(lgt_block, 1) ) )
    if (rank==0) then
        write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
        "lgt_active", product(real(shape(lgt_active)))*4.0e-9, shape(lgt_active)
    endif

    !---------------------------------------------------------------------------
    ! note: 5th dimension in heavy data is block id
    allocate( hvy_active( size(hvy_block, 5) ) )
    if (rank==0) then
        write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
        "hvy_active", product(real(shape(hvy_active)))*4.0e-9, shape(hvy_active)
    endif

    ! setting -1 is required to avoid "ERROR: We try to fetch a light free block ID from the list but all blocks are used on this CPU"
    ! because it looks for a -1 block.
    lgt_block(:,1) = -1


    if (rank == 0) then
        write(*,'("INIT: System is allocating heavy data for ",i7," blocks and ", i3, " fields" )') params%number_blocks, Neqn
        write(*,'("INIT: System is allocating light data for ",i7," blocks" )') number_procs*params%number_blocks
        write(*,'("INIT: System is allocating heavy work data for ",i7," blocks " )') params%number_blocks

        effective_memory = (dble(size(hvy_block)) + & ! real data
        dble(size(lgt_block)+size(lgt_sortednumlist)+size(hvy_neighbor)+size(lgt_active)+size(hvy_active))/2.0 & ! integer (hence /2)
        )*8.0e-9 ! in GB

        if (present(hvy_tmp)) effective_memory = effective_memory + dble(size(hvy_tmp))*8.0e-9 ! in GB
        if (present(hvy_work)) effective_memory = effective_memory + dble(size(hvy_work))*8.0e-9 ! in GB

        ! note we currently use 8byte per real and integer by default, so all the same bytes per point
        write(*,'("INIT: Measured (true) local (on 1 cpu) memory (hvy_block+hvy_work+lgt_block no ghosts!) is ",g15.3,"GB per mpirank")') &
        effective_memory

        write(*,'("INIT-GLOBAL: Measured (true) TOTAL (on all CPU) memory (hvy_block+hvy_work+lgt_block no ghosts!) is ",g15.3,"GB")') &
        effective_memory*dble(params%number_procs)
    end if

end subroutine allocate_tree

!-------------------------------------------------------------------------------

subroutine allocate_forest(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, hvy_active, &
    lgt_sortednumlist, hvy_work, hvy_tmp, hvy_gridQ, hvy_n, lgt_n)
    implicit none

    !> user defined parameter structure
    type (type_params), intent(inout)                   :: params
    !> light data array
    integer(kind=ik), allocatable, intent(out)          :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), allocatable, intent(out)             :: hvy_block(:, :, :, :, :)
    !> heavy temp data: used for saving, filtering etc (work array)
    real(kind=rk), allocatable, optional, intent(out)   :: hvy_tmp(:, :, :, :, :)
    !> grid qtys (qtys that do depend only on the grid, i.e. mask funciton and geometry factors)
    real(kind=rk), allocatable, optional, intent(out)   :: hvy_gridQ(:, :, :, :, :)
    !> heavy work array: used for RHS evaluation in multistep methods (like RK4: 00, k1, k2 etc)
    real(kind=rk), allocatable, optional, intent(out)   :: hvy_work(:, :, :, :, :, :)
    !> neighbor array (heavy data)
    integer(kind=ik), allocatable, intent(out)          :: hvy_neighbor(:,:)
    !> list of active blocks (light data)
    integer(kind=ik), allocatable, intent(out)          :: lgt_active(:,:)
    !> list of active blocks (light data)
    integer(kind=ik), allocatable, intent(out)          :: hvy_active(:,:)
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), allocatable, intent(out)       :: lgt_sortednumlist(:,:,:)
    integer(kind=ik), allocatable, intent(out)          :: hvy_n(:), lgt_n(:)
    ! local shortcuts:
    integer(kind=ik)                                    :: g, Neqn, number_blocks,&
    rank, number_procs,  dim
    integer(kind=ik), dimension(3)                      :: Bs
    integer(kind=ik)    :: rk_steps
    real(kind=rk)       :: effective_memory
    integer             :: status, nrhs_slots, nwork, nx, ny, nz, max_neighbors

    real(kind=rk)      :: maxmem, mem_per_block
    real(kind=rk), parameter ::  nstages = 2.0_rk ! stages for ghost node synching
    character(len=80)  :: memstring
    integer(kind=ik)   :: i

    rank            = params%rank
    Bs              = params%Bs
    g               = params%n_ghosts
    Neqn            = params%n_eqn
    number_procs    = params%number_procs
    nx = Bs(1)+2*g
    ny = Bs(2)+2*g
    nz = Bs(3)+2*g

    if (params%dim==2) then
        nz = 1
        max_neighbors = 16
    else
        max_neighbors   = 74
    endif

    ! 19 oct 2018: The work array hvy_work is modified to be used in "register-form"
    ! that means one rhs is stored in a 5D subset of a 6D array.
    ! Hence, nrhs_slots is number of slots for RHS saving:
    if (params%time_step_method == "RungeKuttaGeneric") then
        nrhs_slots = size(params%butcher_tableau,1)
    elseif (params%time_step_method == "Krylov") then
        nrhs_slots = params%M_krylov +3
    elseif ((params%time_step_method == 'none').or.(params%time_step_method == 'no')) then
        nrhs_slots = 0
    else
        call abort(191018161, "time_step_method is unkown: "//trim(adjustl(params%time_step_method)))
    endif

    nwork = max( 2*Neqn, params%N_fields_saved)

    if (rank == 0) then
        write(*,'(80("_"))')
        write(*,'(A)') "INIT: Beginning memory allocation and initialization."
        write(*,'(A)') "ALLOCATING A FOREST"
        write(*,'(80("_"))')
    endif

    ! Automatic memory management. If specified --memory=0.3GB in the call line,
    if (params%number_blocks < 1) then
        !---------------------------------------------------------------------------
        ! Automatic memory management. If specified --memory=0.3GB in the call line,
        ! wabbit will automatically select the number of blocks per rank to be allocated
        ! to consume this amount of memory. helpful on local machines.
        !---------------------------------------------------------------------------
        ! loop over all arguments until you find the string "--memory" in them (or not)
        ! this ensures compatibility with future versions, as the argument may be anywhere in the call.
        do i = 1, command_argument_count()
            call get_command_argument(i, memstring)
            ! is memory limitation used?
            if ( index(memstring,"--memory=")==1 ) then
                if (params%rank==0) write(*,'("INIT: automatic selection of blocks per rank is active!")')

                ! read memory from command line (in GB)
                read(memstring(10:len_trim(memstring)-2),* ) maxmem

                if (params%rank==0) write(*,'("INIT: total memory: ",f9.4,"GB")') maxmem

                ! memory per MPIRANK (in GB)
                maxmem = maxmem / dble(params%number_procs)

                if (params%rank==0) write(*,'("INIT: memory-per-rank: ",f9.4,"GB")') maxmem

                if (params%dim==3) then
                    mem_per_block = real(Neqn) * real((Bs(1)+2*g)*(Bs(2)+2*g)*(Bs(3)+2*g)) & ! hvy_block
                    + 2.0 * nstages * real(Neqn) * real((Bs(1)+2*g)*(Bs(2)+2*g)*(Bs(3)+2*g) - ((Bs(1))*(Bs(2))*(Bs(3)))) &  ! real buffer ghosts
                    + 2.0 * nstages * real(max_neighbors) * 5 / 2.0 ! int bufer (4byte hence /2)

                    ! hvy_gridQ
                    if ( present(hvy_gridQ) ) then
                        mem_per_block = mem_per_block + real(params%n_gridQ) * real((Bs(1)+2*g)*(Bs(2)+2*g)*(Bs(3)+2*g))
                    endif

                    ! hvy_tmp
                    if ( present(hvy_tmp) ) then
                        mem_per_block = mem_per_block + real(max( 2*Neqn, params%N_fields_saved)) * real((Bs(1)+2*g)*(Bs(2)+2*g)*(Bs(3)+2*g))
                    endif

                    ! hvy_work
                    if ( present(hvy_work) ) then
                        mem_per_block = mem_per_block + real(Neqn) * real(nrhs_slots) * real((Bs(1)+2*g)*(Bs(2)+2*g)*(Bs(3)+2*g))
                    endif


                else
                    mem_per_block = real(Neqn) * real((Bs(1)+2*g)*(Bs(2)+2*g)) & ! hvy_block
                    + 2.0 * nstages * real(Neqn) * real((Bs(1)+2*g)*(Bs(2)+2*g) - (Bs(1)*Bs(2))) &  ! real buffer ghosts
                    + 2.0 * nstages * real(max_neighbors) * 5 / 2.0 ! int bufer (4byte hence /2)

                    ! hvy_gridQ
                    if ( present(hvy_gridQ) ) then
                        mem_per_block = mem_per_block + real(params%n_gridQ) * real((Bs(1)+2*g)*(Bs(2)+2*g))
                    endif

                    ! hvy_tmp
                    if ( present(hvy_tmp) ) then
                        mem_per_block = mem_per_block + real(max(2*Neqn, params%N_fields_saved)) * real((Bs(1)+2*g)*(Bs(2)+2*g))
                    endif

                    ! hvy_work
                    if ( present(hvy_work) ) then
                        mem_per_block = mem_per_block + real(Neqn) * real(nrhs_slots) * real((Bs(1)+2*g)*(Bs(2)+2*g))
                    endif

                endif

                ! in GB:
                mem_per_block = mem_per_block * 8.0e-9

                params%number_blocks = nint( maxmem / mem_per_block)

                if (params%rank==0) then
                    write(*,'("INIT: for the desired memory we can allocate ",i8," blocks per rank")') params%number_blocks
                    write(*,'("INIT: we allocated ",i8," blocks per rank (total: ",i8," blocks) ")') params%number_blocks, &
                    params%number_blocks*params%number_procs
                endif

                if ((params%adapt_mesh .or. params%adapt_inicond) .and. (params%number_blocks<2**params%dim)) then
                    call abort(1909181740,"[ini_file_to_params.f90]: The number of blocks for the --memory specification is lower than 2^d&
                    & and comp is adaptive: we cannot fetch all 2^d blocks on one CPU for merging. Use more memory or less cores.")
                endif
            endif
        end do
    endif


    if (rank==0) then
        write(*,'("INIT: mpisize=",i6)') params%number_procs
        write(*,'("INIT: nwork=",i6)') nwork
        write(*,'("INIT: Bs(1)=",i7," Bs(2)=",i7," Bs(3)=",i7," blocks-per-rank=",i7," total blocks=", i7)') &
        Bs(1),Bs(2),Bs(3), params%number_blocks, params%number_blocks*params%number_procs
        write(*,'("INIT: Allocating a ",i1,"D case.")') params%dim
    endif


    !---------------------------------------------------------------------------
    allocate( hvy_block( nx, ny, nz, Neqn, params%number_blocks ) )
    if (rank==0) then
        write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
        "hvy_block", product(real(shape(hvy_block)))*8.0e-9, shape(hvy_block)
    endif

    !---------------------------------------------------------------------------
    ! work data (Runge-Kutta substeps and old time level)
    if (present(hvy_work)) then
        allocate( hvy_work( nx, ny, nz, Neqn, params%number_blocks, nrhs_slots ) )
        if (rank==0) then
            write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
            "hvy_work", product(real(shape(hvy_work)))*8.0e-9, shape(hvy_work)
        endif
    end if

    if ( present(hvy_tmp) ) then
        allocate( hvy_tmp( nx, ny, nz, nwork, params%number_blocks )  )
        if (rank==0) then
            write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
            "hvy_tmp", product(real(shape(hvy_tmp)))*8.0e-9, shape(hvy_tmp)
        endif
    endif

    if ( present(hvy_gridQ) .and. params%n_gridQ > 0 ) then
        allocate( hvy_gridQ( nx, ny, nz, params%n_gridQ, params%number_blocks )  )
        if (rank==0) then
            write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
            "hvy_gridQ", product(real(shape(hvy_gridQ)))*8.0e-9, shape(hvy_gridQ)
        endif
    endif

    !---------------------------------------------------------------------------
    allocate( hvy_neighbor( params%number_blocks, max_neighbors ) )
    if (rank==0) then
        write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
        "hvy_neighbor", product(real(shape(hvy_neighbor)))*8.0e-9, shape(hvy_neighbor)
    endif

    !---------------------------------------------------------------------------)
    allocate( lgt_block( number_procs*params%number_blocks, params%max_treelevel+EXTRA_LGT_FIELDS) )
    if (rank==0) then
        write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
        "lgt_block", product(real(shape(lgt_block)))*4.0e-9, shape(lgt_block)
    endif

    !---------------------------------------------------------------------------
    allocate( lgt_sortednumlist( size(lgt_block,1), 2, params%forest_size) )
    if (rank==0) then
        write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
        "lgt_sortednumlist", product(real(shape(lgt_sortednumlist)))*4.0e-9, shape(lgt_sortednumlist)
    endif

    !---------------------------------------------------------------------------
    allocate( lgt_active( size(lgt_block,1), params%forest_size ) )
    if (rank==0) then
        write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
        "lgt_active", product(real(shape(lgt_active)))*4.0e-9, shape(lgt_active)
    endif

    allocate( lgt_n(params%forest_size ) )
    if (rank==0) then
        write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
        "lgt_n", product(real(shape(lgt_n)))*4.0e-9, shape(lgt_n)
    endif

    allocate( hvy_n(params%forest_size ) )
    if (rank==0) then
        write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
        "hvy_n", product(real(shape(hvy_n)))*4.0e-9, shape(hvy_n)
    endif
    !---------------------------------------------------------------------------
    ! note: 5th dimension in heavy data is block id
    allocate( hvy_active( size(hvy_block, 5), params%forest_size ) )
    if (rank==0) then
        write(*,'("INIT: ALLOCATED ",A19," MEM=",f8.4,"GB SHAPE=",7(i9,1x))') &
        "hvy_active", product(real(shape(hvy_active)))*4.0e-9, shape(hvy_active)
    endif

    lgt_n = 0
    hvy_n = 0
    ! setting -1 is required to avoid "ERROR: We try to fetch a light free block ID from the list but all blocks are used on this CPU"
    ! because it looks for a -1 block.
    lgt_block(:,1) = -1
    if (rank == 0) then
        write(*,'("INIT: System is allocating heavy data for ",i7," blocks and ", i3, " fields" )') params%number_blocks, Neqn
        write(*,'("INIT: System is allocating light data for ",i7," blocks" )') number_procs*params%number_blocks
        write(*,'("INIT: System is allocating heavy work data for ",i7," blocks " )') params%number_blocks

        effective_memory = (dble(size(hvy_block)) + & ! real data
        dble(size(lgt_block)+size(lgt_sortednumlist)+size(hvy_neighbor)+size(lgt_active)+size(hvy_active))/2.0 & ! integer (hence /2)
        )*8.0e-9 ! in GB

        if (present(hvy_tmp)) effective_memory = effective_memory + dble(size(hvy_tmp))*8.0e-9 ! in GB
        if (present(hvy_work)) effective_memory = effective_memory + dble(size(hvy_work))*8.0e-9 ! in GB

        ! note we currently use 8byte per real and integer by default, so all the same bytes per point
        write(*,'("INIT: Measured (true) local (on 1 cpu) memory (hvy_block+hvy_work+lgt_block no ghosts!) is ",g15.3,"GB per mpirank")') &
        effective_memory

        write(*,'("INIT-GLOBAL: Measured (true) TOTAL (on all CPU) memory (hvy_block+hvy_work+lgt_block no ghosts!) is ",g15.3,"GB")') &
        effective_memory*dble(params%number_procs)
    end if

end subroutine allocate_forest





subroutine deallocate_tree(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, hvy_active, &
    lgt_sortednumlist, hvy_work, hvy_tmp )

    implicit none

    !> user defined parameter structure
    type (type_params), intent(inout)                   :: params
    !> light data array
    integer(kind=ik), allocatable, intent(out)          :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), allocatable, intent(out)             :: hvy_block(:, :, :, :, :)
    real(kind=rk), allocatable, intent(out)             :: hvy_tmp(:, :, :, :, :)
    !> heavy work array
    real(kind=rk), allocatable, optional, intent(out)   :: hvy_work(:, :, :, :, :, :)
    !> neighbor array (heavy data)
    integer(kind=ik), allocatable, intent(out)          :: hvy_neighbor(:,:)
    !> list of active blocks (light data)
    integer(kind=ik), allocatable, intent(out)          :: lgt_active(:)
    !> list of active blocks (light data)
    integer(kind=ik), allocatable, intent(out)          :: hvy_active(:)
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), allocatable, intent(out)       :: lgt_sortednumlist(:,:)

    if (params%rank == 0) then
        write(*,'(80("-"))')
        write(*,'(A)') "FREE: Beginning freeying of memory."
    endif

    if (allocated(hvy_block)) deallocate( hvy_block )
    if (allocated(hvy_work)) deallocate( hvy_work )
    if (allocated(hvy_tmp)) deallocate( hvy_tmp )
    if (allocated(hvy_neighbor)) deallocate( hvy_neighbor )
    if (allocated(lgt_block)) deallocate( lgt_block )
    if (allocated(lgt_sortednumlist)) deallocate( lgt_sortednumlist )
    if (allocated(lgt_active)) deallocate( lgt_active )
    if (allocated(hvy_active)) deallocate( hvy_active )

    if (allocated(params%threshold_state_vector_component)) deallocate( params%threshold_state_vector_component )
    if (allocated(params%input_files)) deallocate( params%input_files )

    if (params%rank == 0) then
        write(*,'(A)') "All memory is cleared!"
        write(*,'(80("-"))')
    endif

end subroutine deallocate_tree



subroutine deallocate_forest(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, hvy_active, &
    lgt_sortednumlist, hvy_work, hvy_tmp, hvy_n, lgt_n )

    implicit none

    !> user defined parameter structure
    type (type_params), intent(inout)                   :: params
    !> light data array
    integer(kind=ik), allocatable, intent(out)          :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), allocatable, intent(out)             :: hvy_block(:, :, :, :, :)
    real(kind=rk), allocatable, intent(out)             :: hvy_tmp(:, :, :, :, :)
    !> heavy work array
    real(kind=rk), allocatable, optional, intent(out)   :: hvy_work(:, :, :, :, :, :)
    !> neighbor array (heavy data)
    integer(kind=ik), allocatable, intent(out)          :: hvy_neighbor(:,:)
    !> list of active blocks (light data)
    integer(kind=ik), allocatable, intent(out)          :: lgt_active(:,:)
    !> list of active blocks (light data)
    integer(kind=ik), allocatable, intent(out)          :: hvy_active(:,:)
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), allocatable, intent(out)       :: lgt_sortednumlist(:,:,:)
    integer(kind=ik), allocatable, intent(out)          :: hvy_n(:), lgt_n(:)

    if (params%rank == 0) then
        write(*,'(80("-"))')
        write(*,'(A)') "FREE: Beginning freeying of memory."
    endif

    if (allocated(hvy_block)) deallocate( hvy_block )
    if (allocated(hvy_work)) deallocate( hvy_work )
    if (allocated(hvy_tmp)) deallocate( hvy_tmp )
    if (allocated(hvy_neighbor)) deallocate( hvy_neighbor )
    if (allocated(lgt_block)) deallocate( lgt_block )
    if (allocated(lgt_sortednumlist)) deallocate( lgt_sortednumlist )
    if (allocated(lgt_active)) deallocate( lgt_active )
    if (allocated(hvy_active)) deallocate( hvy_active )
    if (allocated(hvy_n)) deallocate( hvy_n )
    if (allocated(lgt_n)) deallocate( lgt_n )

    if (allocated(params%threshold_state_vector_component)) deallocate( params%threshold_state_vector_component )
    if (allocated(params%input_files)) deallocate( params%input_files )

    if (params%rank == 0) then
        write(*,'(A)') "All memory is cleared!"
        write(*,'(80("-"))')
    endif

end subroutine deallocate_forest
