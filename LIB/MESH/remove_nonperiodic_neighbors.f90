! If some BC are non-periodic, we have to remove their neighborhood relationship
! by setting it to -1. Otherwise, the sync_ghosts routine will overwrite the redundant points,
! which is an error in this case.
subroutine remove_nonperiodic_neighbors(params, tree_ID)
    implicit none

    type (type_params), intent(in)      :: params               !> user defined parameter structure
    integer(kind=ik), intent(in)        :: tree_ID

    integer(kind=ik) :: k, hvy_id, lgt_id, lgt_id_neighbor, J1, J2, a
    logical                             :: remove
    integer(kind=2)                     :: n_domain(1:3)

    do k = 1, hvy_n(tree_ID)
        ! the block w're looking at ...
        hvy_id = hvy_active(k, tree_ID)
        call hvy2lgt( lgt_id, hvy_id, params%rank, params%number_blocks )
        ! ... and its level
        J1 = lgt_block(lgt_id, params%max_treelevel + IDX_MESH_LVL)

        call get_adjacent_boundary_surface_normal( lgt_block(lgt_id, 1:J1), params%domain_size, params%Bs, params%dim, n_domain )

        ! is this an interior block ? (note: this is completely equivalent to checking if a neighborhood crosses the periodic domain,
        ! because only in this case, n_domain has nonzero value)
        if (maxval(abs(n_domain(1:params%dim)))==0) cycle

        do a = 1, size(hvy_neighbor,2)
            ! if this neighborhood is unused (there is no neighbor in this direction): skip it
            if (hvy_neighbor(hvy_id, a) == -1) then
                cycle
            endif

            ! note the hvy_neighbor stores light IDs
            lgt_id_neighbor = hvy_neighbor(hvy_id, a)
            J2 = lgt_block(lgt_id_neighbor, params%max_treelevel + IDX_MESH_LVL)

            ! Terrible conditionals to determine whether this BC is affected by non-periodic
            ! conditions or not.
            if (params%dim==2) then
                ! sender (neighborhood-codes in hvy_neighbors)
                remove =              (((a==3).or.(a==7).or.(a==8).or.(a==11).or.(a==12)) .and. (.not.params%periodic_BC(1) .and. n_domain(1)==+1) )
                remove = (remove .or. (((a==5).or.(a==1).or.(a==6).or.(a== 9).or.(a==10)) .and. (.not.params%periodic_BC(1) .and. n_domain(1)==-1) ))
                remove = (remove .or. (((a==5).or.(a==2).or.(a==7).or.(a==13).or.(a==14)) .and. (.not.params%periodic_BC(2) .and. n_domain(2)==+1) ))
                remove = (remove .or. (((a==6).or.(a==4).or.(a==8).or.(a==15).or.(a==16)) .and. (.not.params%periodic_BC(2) .and. n_domain(2)==-1) ))
            else
                remove =              (((a==3).or.(a==15).or.(a==17).or.(a== 8).or.(a==12).or.(a==19).or.(a==20).or.(a==23).or.(a==24)) .and. (.not.params%periodic_BC(1) .and. n_domain(1)==+1) )
                remove = (remove .or. (((a==5).or.(a==14).or.(a==10).or.(a==16).or.(a==18).or.(a==21).or.(a==22).or.(a==25).or.(a==26)) .and. (.not.params%periodic_BC(1) .and. n_domain(1)==-1) ))

                remove = (remove .or. (((a==4).or.(a==17).or.(a==18).or.(a==09).or.(a==13).or.(a==21).or.(a==20).or.(a==24).or.(a==25)) .and. (.not.params%periodic_BC(2) .and. n_domain(2)==+1) ))
                remove = (remove .or. (((a==2).or.(a==11).or.(a==15).or.(a==07).or.(a==16).or.(a==22).or.(a==19).or.(a==26).or.(a==23)) .and. (.not.params%periodic_BC(2) .and. n_domain(2)==-1) ))

                remove = (remove .or. (((a==1).or.(a==10).or.(a==07).or.(a==08).or.(a==09).or.(a==21).or.(a==22).or.(a==20).or.(a==19)) .and. (.not.params%periodic_BC(3) .and. n_domain(3)==+1) ))
                remove = (remove .or. (((a==6).or.(a==12).or.(a==14).or.(a==11).or.(a==13).or.(a==25).or.(a==26).or.(a==23).or.(a==24)) .and. (.not.params%periodic_BC(3) .and. n_domain(3)==-1) ))


                remove = (remove .or. (((a==35).or.(a==36).or.(a==37).or.(a==38).or.(a==67).or.(a==68).or.(a==61).or.(a==62).or.(a==71).or.(a==72).or.(a==53).or.(a==54)) .and. (.not.params%periodic_BC(1) .and. n_domain(1)==+1) ))
                remove = (remove .or. (((a==43).or.(a==45).or.(a==44).or.(a==46).or.(a==73).or.(a==74).or.(a==69).or.(a==70).or.(a==65).or.(a==66).or.(a==57).or.(a==58)) .and. (.not.params%periodic_BC(1) .and. n_domain(1)==-1) ))

                remove = (remove .or. (((a==42).or.(a==40).or.(a==41).or.(a==39).or.(a==64).or.(a==63).or.(a==71).or.(a==72).or.(a==56).or.(a==55).or.(a==73).or.(a==74)) .and. (.not.params%periodic_BC(2) .and. n_domain(2)==+1) ))
                remove = (remove .or. (((a==33).or.(a==31).or.(a==34).or.(a==32).or.(a==52).or.(a==51).or.(a==67).or.(a==68).or.(a==59).or.(a==60).or.(a==69).or.(a==70)) .and. (.not.params%periodic_BC(2) .and. n_domain(2)==-1) ))

                remove = (remove .or. (((a==29).or.(a==28).or.(a==27).or.(a==30).or.(a==56).or.(a==55).or.(a==54).or.(a==53).or.(a==51).or.(a==52).or.(a==58).or.(a==57)) .and. (.not.params%periodic_BC(3) .and. n_domain(3)==+1) ))
                remove = (remove .or. (((a==50).or.(a==47).or.(a==49).or.(a==48).or.(a==60).or.(a==59).or.(a==61).or.(a==62).or.(a==64).or.(a==63).or.(a==65).or.(a==66)) .and. (.not.params%periodic_BC(3) .and. n_domain(3)==-1) ))
            endif

            if (remove) then
                ! deactivate this neighborhood relation
                hvy_neighbor(hvy_id, a) = -1
            endif
        enddo
    enddo
end subroutine
