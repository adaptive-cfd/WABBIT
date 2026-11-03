subroutine SaveBeamData( time )
    !use vars
    implicit none
    real (kind=rk), intent (in) :: time
    character(len=16) :: format_ns1
    character(len=3)  :: ns1_string
    character(len=1)  :: beamstr
    real (kind=rk) :: alpha, alpha_t, alpha_tt
    real (kind=rk), dimension(1:6) :: LeadingEdge !LeadingEdge: x, y, vx, vy, ax, ay (Array)
    integer :: n,step,i

    ! set up formats
    write(ns1_string, '(I3)') ns+1
    format_ns1 = '('//ns1_string//'(es15.8,1x))'

    !-- loop over beams
    do i = 1, nBeams
        !-- for naming files..
        write (beamstr,'(i1)') i

        !-- save trailing edge data
        open  (14, file = 'beam_data'//beamstr//'.t', status = 'unknown',position='append')
        write (14, '(21(es15.8,1x))') &
            time,&
            beams(i)%pressure_new(ns-1),&
            beams(i)%pressure_new(3*ns/4),&
            beams(i)%pressure_new(ns/2),&
            beams(i)%x(ns-1),&
            beams(i)%y(ns-1),&
            beams(i)%vx(ns-1),&
            beams(i)%vy(ns-1),&
            beams(i)%theta(ns-1),&
            beams(i)%theta_dot(ns-1),&
            beams(i)%E_kinetic,&
            beams(i)%E_elastic,&
            beams(i)%E_pot
        close (14)

        !-- save leading edge motions
        call mouvement(time, alpha, alpha_t, alpha_tt, LeadingEdge, beams(i) )
        open (14, file = 'mouvement'//beamstr//'.t', status = 'unknown',position='append')
        write (14, '(10(es15.8,1x))') time, alpha, alpha_t, alpha_tt, (LeadingEdge(n), n=1,6)
        close (14)

        open (14, file = 'beam_x'//beamstr//'.t', status = 'unknown',position='append')
        write (14, format_ns1) time, beams(i)%x(0:ns-1)
        close (14)

        open (14, file = 'beam_y'//beamstr//'.t', status = 'unknown',position='append')
        write (14, format_ns1) time, beams(i)%y(0:ns-1)
        close (14)

        open (14, file = 'beam_vx'//beamstr//'.t', status = 'unknown',position='append')
        write (14, format_ns1) time, beams(i)%vx(0:ns-1)
        close (14)

        open (14, file = 'beam_vy'//beamstr//'.t', status = 'unknown',position='append')
        write (14, format_ns1) time, beams(i)%vy(0:ns-1)
        close (14)

        open (14, file = 'beam_p'//beamstr//'.t', status = 'unknown',position='append')
        write (14, format_ns1) time, beams(i)%pressure_new(0:ns-1)
        close (14)

        open (14, file = 'beam_theta'//beamstr//'.t', status = 'unknown',position='append')
        write (14, format_ns1) time, beams(i)%theta(0:ns-1)
        close (14)
    enddo


end subroutine SaveBeamData
