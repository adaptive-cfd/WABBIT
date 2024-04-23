!===========================================================================
!> this module implements the MPI_bridge between fluid and particle on the side
!> of the fluid
!> \details
!> 1. it initializes the bridge see init_bridge()
!> 2. it provides a routine which sends the treecode and further information
!>    to the particle side of the bridge. This is needed to map the
!>    particle coordinates onto the right processors subdomain
!> 3. interpolates the velocity field \f$(u_\alpha)_{ijk}\f$ to a requested position \f$x \in R^3\f$
!>    and sends it back to the requesting processor
!-------------------------------------------------------------------------!!!!
module module_bridge_interface

    use mpi
    use module_bridge
    use module_treelib
    use module_mesh
    use module_ini_files_parser_mpi
    ! use constants
    ! use methods
    ! use geometry
    ! use parametersParticle
    ! use initializationParticle
    ! use commonParticle
    ! use chargedParticle
    ! use bridge
    use module_forestMetaData

    use module_params

    implicit none
    private

    ! ! Declarations

    ! !! Variable relatives to the other world (fluid domain)
    ! !! - array corresponding to the parameters of the grid (number of points in each direction)
    ! !!   1. total number of grid points
    ! !!   2. number of grid points in one process
    ! integer         , dimension(3,2), private  :: nbGridPointsFluid
    !! - array corresponding to the number of fluid processes in each direction (1,1,1 if only one process)
    integer         , dimension(3)  , private  :: processesParticle
    ! !! - streching parameter to modify the coordinates - WARNING: only one to begin
    ! double precision,                 private  :: stretchingParameter

    ! !! Useful subroutines
    ! public :: getParticleProcesses
    ! public :: getProcessWithPosition
    ! public :: exchangeParametersParticle
    public :: init_bridge
    public :: initialize_communicator
    public :: send_lgt_data
    public :: serve_data_request
    public :: position_to_lgt_id
    ! public :: signalEndCalculation
    ! public :: particleDistributionOverFluidProcesses
    ! public :: exchangeDataParticle




contains

    !===========================================================================
    !> \brief read out from the command line the amount of fluid processes to start
    subroutine getParticleProcesses(params)

        !> structure containing the parameters
        type(type_params), intent(in)         :: params

        ! (total) number of fluid processes
        character(10)                            :: nbParticleProcesses

        ! Read out the total number of particle processes
        !to start from the given command line
        call get_command_argument(2, nbParticleProcesses)

        ! Assign this number in the first value of processesParticle (1 otherwise)
        processesParticle = 1
        read(nbParticleProcesses,'(i10)') processesParticle(1)
        if ( params%rank==0 ) then
            write(*,*) 'Number of Particle Processes: ', processesParticle(1)
        endif
    end subroutine getParticleProcesses
    !===========================================================================


    !===========================================================================
    subroutine send_fixed_params (myBridge, params)
        !> \brief this routine sends all fixed parameters (e.g. \f$L_x,L_y)\f$ to the particle world
        !! Subroutine-declarations
        type(bridgeMPI)  , intent(in)                :: myBridge           !< type bridge on the particle side
        type(type_params), intent(in)                :: params             !> structure containing the parameters
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer                                      :: ierr               ! MPI communication error
        integer         , dimension(MPI_STATUS_SIZE) :: status             ! MPI communication status
        double precision, dimension(3,2)             :: Domain     ! received domain parameters to transform
        integer         , dimension(5)             :: discretizationParams ! Number of grid points and processes
        character(len=cshort)                            :: geometry

        ! Only the root particle process receive the parameters from the root fluid process
        if (myBridge%myWorldRank == 0) then                              ! check if the process is the root process

            if (sum(params%Bs(1:params%dim) - params%Bs(1)) .ne. 0) then
                call abort(24010149, "You tried to send something with different Bs in different directions, but this is not (yet) implemented here.")
            endif
            ! - receive the discretization parameters
            discretizationParams=(/ params%Bs(1),  &
            params%number_blocks,       &
            params%Jmax,       &
            params%dim,                 &
            params%number_procs/)
            call MPI_send(discretizationParams, 5, MPI_integer, myBridge%minOtherWorldRank, &
            parameters_delivery, myBridge%commonWorld,  ierr)

            !! - receive the domain parameters
            Domain(:,1)=(/params%domain_size(1),params%domain_size(2),params%domain_size(3)/) !> Domain size
            Domain(:,2)=(/0, 0, 0/)                         !> Domain origin
            call MPI_send(Domain, 6, MPI_double_precision, myBridge%minOtherWorldRank, &
            parameters_delivery, myBridge%commonWorld, ierr)

            !! - receive the coordinate system type parameters (overwrite the one given in the particle parameters)
            geometry='cartesian' ! 'cartesian' / 'cylindrical'
            call MPI_send(geometry, 80, MPI_character, myBridge%minOtherWorldRank, &
            parameters_delivery, myBridge%commonWorld, ierr)
        end if
    end subroutine send_fixed_params
    !===========================================================================



    !===========================================================================
    subroutine send_lgt_data (params)
        !> \brief this routine sends all fixed parameters (e.g. \f$L_x,L_y)\f$ to the particle world
        ! Subroutine-declarations
        type(type_params), intent(in)         :: params             !< structure containing the parameters

        !--------------------------------------------------------------------------------
        integer                              :: ierr               ! MPI communication error
        integer , dimension(MPI_STATUS_SIZE) :: status             ! MPI communication status
        integer                              :: n,m
        integer :: tree_ID = 1

        ! Only the root particle process receive the parameters from the root fluid process
        if ( params%bridge%myWorldRank == 0 ) then                              ! check if the process is the root process
            ! number of blocks
            n = size(lgt_block,1)

            ! number of columns in matrix
            m = EXTRA_LGT_FIELDS

            ! send number of active and maximal number of blocks
            call MPI_send((/lgt_n/), 1, MPI_integer, &
            params%bridge%minOtherWorldRank, parameters_delivery, &
            params%bridge%commonWorld,  ierr )

            ! send list of active blocks
            call MPI_send(lgt_active, lgt_n(tree_ID), MPI_integer, &
            params%bridge%minOtherWorldRank, parameters_delivery, &
            params%bridge%commonWorld,  ierr )

            ! send light data
            call MPI_send(lgt_block(:,1:m), n*m, MPI_integer, &
            params%bridge%minOtherWorldRank, parameters_delivery, &
            params%bridge%commonWorld, ierr )
        end if
    end subroutine send_lgt_data
    !===========================================================================




    !===========================================================================
    subroutine serve_data_request(hvy_block, hvy_work, params)
        !> \brief this routine sends all fixed parameters (e.g. \f$L_x,L_y)\f$ to the particle world
        ! Subroutine-declarations
        type(type_params), intent(in)       :: params             !< structure containing the parameters
        !> heavy data array - block data
        real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
        !> heavy work data array - block data
        real(kind=rk), intent(inout)        :: hvy_work(:, :, :, :, :)
        !--------------------------------------------------------------------------------
        integer                              :: ierr ,k,maxpoints            ! MPI communication error
        integer , dimension(MPI_STATUS_SIZE) :: status             ! MPI communication status
        integer                              :: Ndat,request,sta
        integer         , allocatable        :: statuses(:,:)   ! MPI communication status array
        integer         , allocatable        :: requests(:)   !  array
        integer                              :: send_id           !! rank of sender
        double precision, allocatable        :: distributedParticles(:,:) ! stores the distributed value
        double precision, allocatable        :: u_inter(:,:) ! stores the distributed value

        character(1)                         :: buf ! Message sent to the fluid world

        !> \detail
        !! \c serve_data_request is
        !!    - if \c MPI_TAG is \c data_positionDelivery
        !!          -> interpolate fluid data to requestet position
        !!          -> send fluid data back to the sender (i.e. particle rank)
        !!    - if \c MPI_TAG is \c end_communication stop waiting for requests
        maxpoints=100000
        k=0

        allocate(distributedParticles(4,maxpoints))
        distributedParticles=-99
        allocate(requests(params%bridge%otherWorldSize+1))
        do
            call MPI_probe(MPI_ANY_SOURCE,MPI_ANY_TAG,params%bridge%otherWorld,status,ierr)

            !=========================================
            select case(status(MPI_TAG))
                !======= case 1: DATA delivery
            case(data_positionDelivery)
                !=======
                send_id=status(MPI_SOURCE)
                !  write(*,*) 'send_id:',send_id
                k=k+1
                call MPI_irecv(distributedParticles,4*maxpoints , MPI_double_precision, send_id, &
                data_positionDelivery,params%bridge%otherWorld,request,ierr)
                call MPI_Get_count( status,MPI_double_precision , Ndat , ierr)
                !this is necessary to make sure that the data is actually writen to the local proc memory
                call MPI_wait( request, status, ierr)
                ! divide Ndat to get the particle Number
                Ndat=Ndat/4

                !write(*,"('wabbit rank ',i6,' - recieved data from pig rank',i6,' x1=(', f6.3, f6.3')')")  params%Bridge%myWorldRank, &
                !send_id,  distributedParticles(1:2,1)

                ! write(*,*) "my wabbit rank= ",params%bridge%commonWorldRank
                ! write(*,*) "number of positions to interpolate=",Ndat
                ! write(*,*)  "WABBIT x1=" ,distributedParticles(1:2,1)
                ! write(*,*)  "WABBIT xend=" ,distributedParticles(1:2,Ndat)
                call interpolate_data(hvy_block, hvy_work, params,distributedParticles(:,1:Ndat),u_inter)

                call MPI_isend(u_inter, 6*Ndat, MPI_double_precision, send_id, data_fieldDelivery, &
                params%bridge%otherWorld,requests(k),  ierr)
                if (allocated(u_inter)) deallocate(u_inter)
                ! write(*,*) 'number of particles:',Ndat/4,'positiondata=', distributedParticles(:,1)
                !======= case 2: end communication
            case (end_communication)
                !=======
                call MPI_recv(buf, 1, MPI_character, params%bridge%minOtherWorldRank,end_communication,params%bridge%otherWorld, status,ierr)
                allocate(statuses(MPI_STATUS_SIZE,k))
                ! write(*,*)"k=",k
                call MPI_waitall(k, requests, statuses, ierr)
                !write(*,*) "fluidrank",params%bridge%myWorldRank,"processed data request"
                exit
                !  case ()
            end select
            !=========================================
        enddo

        deallocate(distributedParticles)
        deallocate(requests)


    end subroutine serve_data_request
    !===========================================================================







    !> \brief Initialize global communicator
    !> \details There are two possible scenarios:
    !!      + Bridge Modus:
    !!            - 2 communicators for the two MPI-WORLDS are created
    !!      + Default Modus:
    !!            - the MPI_COMM_WORLD is duplicated to WABBIT_COMM:
    !!            - thus it creates a new communicator that has a new communication
    !!              context but contains the same group of processes as the
    !!              input communicator
    subroutine initialize_communicator(params)

        !> given structure containing the parameters
        type(type_params), intent(inout), optional  :: params
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer                                     :: ierr

        ! Initialize the particle MPI world (according to the given kind of hierarchy)
        if (params%bridge_exists) then                             ! if there is a single WABBIT_COMM

            if (params%rank==0) then
                write(*,*) ''
                write(*,*) '-----------------'
                write(*,*) 'Modus: BRIDGE'
                write(*,*) '-----------------'
            end if
            params%bridge%myWorld=WABBIT_COMM
            call init_bridge(params%bridge,params)                              ! initialization using split
        else                                                              ! if there are different WABBIT_COMM
            call MPI_comm_dup(WABBIT_COMM, params%bridge%myWorld, ierr)
            call MPI_comm_size(params%bridge%myWorld, params%bridge%myWorldSize, ierr)
            call MPI_comm_rank(params%bridge%myWorld, params%bridge%myWorldRank, ierr)
        end if
        ! after communicator is created set it to the global communicator of WABBIT
        WABBIT_COMM         = params%bridge%myWorld
        params%rank         = params%bridge%myWorldRank
        params%number_procs = params%bridge%myWorldSize
    end subroutine initialize_communicator

    !===========================================================================


    !===========================================================================

    subroutine init_bridge(myBridge,params)
        !! call every necessary subroutine / function to initialize the fluid world
        !! Subroutine-declarations
        type(bridgeMPI)  , intent(out)              :: myBridge            ! type bridge on the fluid side
        type(type_params), intent(inout)            :: params         ! given structure containing the parameters
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! integer                                     :: unitNum = 1         ! file handle
        integer                                     :: ierr                ! MPI communication error

        !    paramsInit = params                                       ! associate the effective parameters structure to the given one

        ! Get the total number of fluid processes
        call getParticleProcesses(params)
        ! Initialize the particle MPI world (according to the given kind of hierarchy)
        if (params%bridgeCommonMPI) then                             ! if there is a single WABBIT_COMM
            call createMPIWorlds(myBridge, 1)                              ! initialization using split
        else                                                             ! if there are different WABBIT_COMM
            if (params%bridgeFluidMaster) then
                !      write(*,*) "write particlecommand: ", params%particleCommand                   ! if the fluid side is the master
                call createMPIWorldsMaster(myBridge, processesParticle(1), params%particleCommand) ! initialization using spawn (slave)
            else                                                           ! if the fluid side is the slave
                call createMPIWorldsSlave(myBridge)                          ! initialization using spawn (master)
            end if                                                         ! end condition regarding the role of the particle side
        end if                                                           ! end condition regarding the WABBIT_COMM structure
        !   write(*,*) 'Fluid side - MPI worlds are created'
        call send_fixed_params(myBridge,params)
        !   write(*,*) 'fixed params send to particle world'
        ! ! Make sure the given number of process for the particle world match the number of processes in the parameters file
        ! if (myBridge%myWorldSize /= (paramsInit%NumberProcesses)) then   ! condition on the total number of processes
        !   if (myBridge%myWorldRank == 0) then                            ! only display the error message if this is the root particle process
        !     write(*,*) "The number of processes given in mpiexec for the particles does't match the &
        !                 &number of processes given in the parameters file. Stop!"
        !     write(*,*) "Given number of particle processes: ", myBridge%myWorldSize
        !     write(*,*) "Expected number of particle processes: ", paramsInit%NumberProcesses
        !     call MPI_abort(myBridge%commonWorld, 99, ierr)               ! Abort the execution
        !   end if                                                         ! end condition on the rank in the fluid world
        ! end if                                                           ! end condition on the total number of processes

        !   ! Get the total fluid domain values from the fluid world
        ! !   write(*,*) 'Particle side - Try to communicate parameters with fluid'
        !   call exchangeParametersParticle(myBridge, paramsInit)
        ! !   if (myBridge%myWorldRank == 0) then
        ! !   write(*,*) 'Particle side - I got from NSF'
        ! !   write(*,*) 'Number of discretization points: ', nbGridPointsFluid
        ! !   write(*,*) 'Number of NSF processes: ', processesFluid
        ! !   end if

    end subroutine init_bridge
    !===========================================================================




    ! !===========================================================================
    ! !> \brief function returns the lgt_id to the corresponding treecode
    ! integer(kind=tsize) function treecode2lgt_id(maxtreelevel,treecode)

    !     ! input:
    !     implicit none
    !     integer(kind=ik), intent(in)    :: treecode(:)      !< treecode
    !     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !     integer         , allocatable           :: meshlevel(:)        ! mesh level of the blocks
    !     integer                                 :: nblock,i,j,k    ! Number of blocks
    !     integer                                 :: maxtreelevel, tree_ID=1


    !     nblock= size(lgt_block,1)
    !     allocate(meshlevel(lgt_n(tree_ID)))

    !     meshlevel(:) = lgt_block(lgt_active(1:lgt_n(tree_ID), tree_ID), maxtreelevel+1)

    !     do k=1,lgt_n(tree_ID) !loop over blocks k
    !         ! only loop over active blocks
    !         j = lgt_active(k, tree_ID)
    !         i = 1
    !         do while(lgt_block(j,i)==treecode(maxtreelevel-i) .and. i<=meshlevel(k))
    !             i=i+1
    !         enddo

    !         if ((i-1)==meshlevel(k)) then
    !             treecode2lgt_id=j
    !             ! write(*,*) 'light_id=',treecode2lgt_id

    !             !return
    !         endif

    !     enddo

    ! end function
    ! !===========================================================================



    !===========================================================================

    !> \brief return the light_id for the given position
    integer(kind=tsize) function position_to_lgt_id (position, params)
        ! Function-declarations
        double precision, dimension(3)  , intent(in)    :: position         !< position vector \f$ 0<position(i)\le 1\f$
        type(type_params),                intent(in)    :: params           !< params of infile
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        integer, dimension(3)           :: blockLength     ! number of grid points in each block
        integer, dimension(4,8)         :: procIndex       ! index of the processes in each direction for the 8 possibilities and activity flag
        integer                         :: d               ! count integer for dimension
        double precision, dimension(3)  :: dx,x0,x1              ! seperation
        integer                         :: k, tree_ID=1
        !> \detail
        position_to_lgt_id=-999999

        do k = 1, lgt_n(tree_ID) ! loop over all active blocks

            call get_block_spacing_origin( params, lgt_active(k, tree_ID), x0, dx )
            x1=x0(:)+dx(:)*(params%Bs(1)-1)

            ! check if position is inside the block with origin x0 and size dx
            ! if yes:   found the right block with its lgt_id
            ! if not:   go to the next block
            !                    x0(1) -------------------------- x1(1)
            !                         |                           |
            !                         |         lgt_id            |
            !                         |               *           |
            !                         |               (position)  |
            !                         |                           |
            !                    x0(2) -------------------------- x1(2)

            if (min(position(1)-x0(1),position(2)-x0(2))>=0 .and. max(position(1)-x1(1),position(2)-x1(2))<=0) then
                if (params%dim == 3) then ! in 3d we have to check the 3. component as well
                    ! 3D case
                    if (position(3)>x0(3) .and. position(3)<x1(3)) then
                        position_to_lgt_id=lgt_active(k, tree_ID)
                        return
                    endif
                else
                    ! 2D case
                    position_to_lgt_id=lgt_active(k, tree_ID)
                    return
                endif ! end 3dcase
            endif


        enddo
        write(*,'("[bridgefluid.f90:] No block found for position= (", f6.3,",",f6.3,") STOP!")') position(1:2)
        write(*,*) position
        call abort(272372,'STOP!')
    end function position_to_lgt_id

    !===========================================================================







    !===========================================================================

    subroutine interpolate_data(hvy_block, hvy_work, params,positions,u_inter)
        !! Send a message to the other world to signal the end of the communication session
        !! Subroutine-declarations

        !> heavy data array - block data
        real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
        !> heavy work data array - block data
        real(kind=rk), intent(inout)        :: hvy_work(:, :, :, :, :)
        !>
        real(kind=rk), allocatable, intent(inout)   :: u_inter(:,:)

        real(kind=rk), intent(inout)        :: positions(:,:)

        type(type_params),                intent(in)    :: params           !< params of infile

        !================================================================

        integer(kind=ik)                    ::k,Nr_particle,lgt_id
        ! block index
        integer(kind=ik)                    ::ibx, iby, ibz

        integer(kind=ik), allocatable       ::particle_id(:)

        real(kind=rk)                      ::u,v,w,rho,p

        real(kind=rk), dimension(1:3)       :: x0, dx

        !================================================================

        Nr_particle=size(positions,2)

        ! allocate state vector
        ! state vector consists of 6 components:
        !   - space components of velocity
        !   - preasure
        !   - density
        if (.not. allocated(u_inter)) allocate(u_inter(6,Nr_particle))
        !allocate particle ids
        allocate(particle_id(Nr_particle))


        ! convert position to hvy_id

        !write(*,*) "dimensions=",size(hvy_block,1),size(hvy_block,2),size(hvy_block,3),size(hvy_block,4  ),size(hvy_block,5)
        do k=1,Nr_particle
            if ( minval(positions(1:2,k))<1e-10 ) then
                write(*,*)positions(1:2,k),k
            endif

            particle_id(k)   = position_to_lgt_id(positions(:,k),params)
            !! lgt_id
            lgt_id=particle_id(k)

            !! hvy_id
            particle_id(k)   = particle_id(k)-params%bridge%myWorldRank*params%number_blocks
            !  write(*,'("hvy_id=",i6,i6,i6)')particle_id(k),params%bridge%myWorldRank,params%number_blocks
            call get_block_spacing_origin( params, lgt_id, x0, dx )
            !!! calculate grid point ibx, iby, ibz
            ibx   = int( (positions(1,k)-x0(1))/dx(1) ) + params%g +1
            iby   = int( (positions(2,k)-x0(2))/dx(2) ) + params%g +1
            if (params%dim == 3) then
                ibz = int( (positions(3,k)-x0(3))/dx(3) ) + params%g +1
            else
                ibz = 1
            endif
            ! if ( positions(1,k)>1.0_rk .and. params%rank==0 ) then
            !   write(*,'("hvy_id =", i6, " x=", f6.3," [xmin,xmax]=[",f6.3,",",f6.3,"]")')particle_id(k),positions(1,k),x0(1), x0(1)+dx(1)*(params%Bs(1)-1)
            !   write(*,'("hvy_id =", i6, " y=", f6.3," [ymin,ymax]=[",f6.3,",",f6.3,"]")')particle_id(k),positions(2,k),x0(2), x0(2)+dx(2)*(params%Bs(2)-1)
            ! endif
            ! remark: the first 3 dimensions of hvy_block project the grid structure into the array
            ! this means:

            !        * *  * * * * * * *  * *
            !        * *  * * * * * * *  * *
            !            --------------
            !        * *| + + + + + + +| * *
            !        * *| + + + + + + +| * *  -->
            !        * *| + + + + + + +| * *  -->  hvy(ibx,iby,ibz,...)
            !        * *| + + + + + + +| * *  -->
            !        * *| + + + + + + +| * *
            !        * *| + + + + + + +| * *      * ... ghost_node
            !        * *| + + + + + + +| * *      + ... inner grid point
            !            --------------           | ... indication of boundary
            !     ^  * *  * * * * * * *  * *
            ! iby |  * *  * * * * * * *  * *
            !        1 2  3 4 5 6 7 8 9  10 11
            !     -->
            !    ibx

            ! write(*,*)"ibx =",ibx,"iby =",iby,"ibz =",ibz

            !interpolate:
            !> \todo write interpolation
            rho = hvy_block(ibx, iby, ibz,1, particle_id(k) )**2

            u   =   hvy_block(ibx, iby, ibz,2, particle_id(k) ) &
            / hvy_block(ibx, iby, ibz,1, particle_id(k) )
            v   =   hvy_block(ibx, iby, ibz,3, particle_id(k) ) &
            / hvy_block(ibx, iby, ibz,1, particle_id(k) )
            if (params%dim == 3) then
                if (params%n_eqn /= 5) then
                    call abort(333990,"[bridgefluid] number of data fields is less then 5, Stop")
                else
                    w   =   hvy_block(ibx, iby, ibz,4, particle_id(k) ) &
                    / hvy_block(ibx, iby, ibz,1, particle_id(k) )
                    p   =   hvy_block(ibx, iby, ibz,5, particle_id(k) )
                    u_inter(:,k)   = hvy_block(ibx, iby, ibz,1:params%n_eqn, particle_id(k) )
                endif
            else
                w   = 0
                p   =   hvy_block(ibx, iby, ibz,4, particle_id(k))
            endif

            ! u_inter is the interpolated data which will be send back to PIG
            u_inter(1,k)   = rho
            u_inter(2,k)   = u
            u_inter(3,k)   = v
            u_inter(4,k)   = w
            u_inter(5,k)   = p
            u_inter(6,k)   = positions(4,k)
            !write(*,'("Nr:",i6,"   rho=",f9.2,"  ux=",f9.2)') int(u_inter(6,k)),u_inter(1,k),u_inter(2,k)
        enddo
        deallocate(particle_id)

    end subroutine interpolate_data

    !===========================================================================

    subroutine signalEndCalculation (myBridge)
        !! Send a message to the other world to signal the end of the calculation session
        !! Subroutine-declarations
        type(bridgeMPI), intent(in)                     :: myBridge        ! type bridge on the particle side
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        character(1)                                    :: endCalculationMessage ! Message sent to the fluid world
        integer                                         :: ierr            ! MPI communication error
        integer                                         :: nProc           ! loop index (fluid process rank)

        ! Wait for all particle processes
        call MPI_barrier(myBridge%myWorld, ierr)

        ! Only the root particle process sent the message
        ! NOTE: one should try to implement this via an intercommunicator broadcast
        if (myBridge%myWorldRank == 0) then                              ! only the root process sends the message
            !! - assign the message to send (but this value is not important)
            endCalculationMessage = 'b'
            !! - send the message to the fluid processes
            do nProc = 0, myBridge%otherWorldSize - 1                      ! loop over all fluid processes
                call MPI_send(endCalculationMessage, 1, MPI_character, nProc, end_calculation, &
                myBridge%otherWorld, ierr)
            end do                                                         ! end loop over all fluid processes
        end if                                                           ! end condition regarding the rank

    end subroutine signalEndCalculation



end module module_bridge_interface
