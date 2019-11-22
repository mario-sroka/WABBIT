subroutine sync_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )
    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(in)        :: hvy_neighbor(:,:)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    real(kind=rk) :: t0

    ! loop variable
    integer(kind=ik)                    :: k, i
    ! dummy variables
    integer(kind=ik)                    :: lgt_id, level, ix, iy, iz, Bs(3), g, d
    real(kind=rk)                       :: dx(3), x0(3)

    t0 = MPI_wtime()

    ! readability
    Bs = params%Bs
    g  = params%n_ghosts
    d  = params%dim

    call synchronize_ghosts_generic_sequence( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n )

    ! ----------------------------------------------------------------------------------------------------------------------
    ! HACK: with non-periodic boundaries some ghost nodes can not be synchronized
    ! for these cases we need to set the ghost nodes with meaningful values
    ! otherwise a lot of subroutines need proper one-sided mechanism (filtering, thresholding, statistics, conversions, ...)
    ! \todo: remove or rework?
    ! \todo: use block coordinates instead of block origins?
    
    ! ghost nodes are overwritten with redundant node value (\todo: actually need to overwrite redundant node?)
    ! check if block is located at domain boundary and if domain boundary is not periodic

    ! work only with reactive navier stokes physics
    if ( params%physics_type == 'reactive_navier_stokes' ) then

        ! loop over all active blocks
        do k = 1, hvy_n

            ! (I): light id
            call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )

            ! (II): mesh level
            level   = lgt_block( lgt_id, params%max_treelevel + IDX_MESH_LVL )

            ! (III): spacing
            dx(1:d) = 2.0_rk**(-level) * params%domain_size(1:d) / real(Bs(1:d)-1, kind=rk)

            ! (IV): origin
            call decoding( lgt_block( lgt_id, 1:level), ix, iy, iz, level)
            x0 = real( ((/ix,iy,iz/) - 1)*(Bs-1) ,kind=rk) * dx

            ! (V): check different directions
            ! -------------------------------

            ! x+
            if ( ( ( x0(1) + dx(1)*real(Bs(1)-1,kind=rk)  - params%domain_size(1) ) + 1e-14_rk ) > 0.0_rk ) then
                ! check non periodicity
                if ( .NOT.(params%periodic_BC(1)) ) then
                    do i = Bs(1)+g, Bs(1)+2*g
                        hvy_block(i, :, :, :, hvy_active(k)) = hvy_block(Bs(1)+g-1, :, :, :, hvy_active(k))
                    end do
                end if
            end if

            ! x-
            if ( ( ( x0(1) - dx(1) ) - 1e-16_rk ) < 0.0_rk ) then
                ! check non periodicity
                if ( .NOT.(params%periodic_BC(1)) ) then
                    do i = 1, g+1
                        hvy_block(i, :, :, :, hvy_active(k)) = hvy_block(g+2, :, :, :, hvy_active(k))
                    end do
                end if
            end if

        end do
        ! ----------------------------------------------------------------------------------------------------------------------

    end if

    call toc( "WRAPPER: sync ghosts", MPI_wtime()-t0 )
end subroutine sync_ghosts
