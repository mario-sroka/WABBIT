
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \file
!> \callgraph
!> \brief wrapper for filter
!> \version 0.5
!> \author Pkrah
!! \date 30/04/18 - create
!!
!
!**********************************************************************************************

subroutine filter_wrapper(time, params, hvy_block, hvy_tmp, lgt_block, hvy_active, hvy_n)

!----------------------------------------------------------------------------------------------
! modules

!----------------------------------------------------------------------------------------------
! variables

   implicit none

    !> time variable
    real(kind=rk), intent(in)           :: time
    !> user defined parameter structure, hvy_active
    type (type_params), intent(in)      :: params
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy temp data: used for saving, filtering, and helper qtys (reaction rate, mask function)
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    !> global integral
    real(kind=rk), dimension(3)         :: volume_int

    !> spacing and origin of a block
    real(kind=rk), dimension(3)         :: dx, x0
    ! loop variables
    integer(kind=ik)                    :: k, dF, neqn, lgt_id
    ! grid parameter, error variable
    integer(kind=ik)                    :: g
    integer(kind=ik), dimension(3)      :: Bs
    !  surface normal
    integer(kind=2) :: surface(3)=0
    integer         :: level


    ! grid parameter
    Bs    = params%Bs
    g     = params%n_ghosts

    if ( params%filter_type == 'spectral' .and. params%physics_type == 'reactive_navier_stokes' ) then

        ! notes: reactive navier stokes spectral filtering needs all data blocks
        call filter_meta(params%physics_type, time, hvy_block(:,:,:,:,:), g, x0, dx, hvy_tmp(:,:,:,:,:), surface) 

    else

        ! traditional filtering
        do k = 1, hvy_n
          ! convert given hvy_id to lgt_id for block spacing routine
          call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )

          ! level of the block:
          level = lgt_block(lgt_id, params%max_treelevel+IDX_MESH_LVL)

          if ((params%filter_only_maxlevel .and. level==params%max_treelevel) .or. .not. params%filter_only_maxlevel) then

              ! get block spacing for RHS
              call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )
    
              if ( .not. All(params%periodic_BC) ) then
                ! check if block is adjacent to a boundary of the domain, if this is the case we use one sided stencils
                call get_adjacent_boundary_surface_normal(params, lgt_id, lgt_block, params%max_treelevel, surface)
              endif

              call filter_meta(params%physics_type, time, hvy_block(:,:,:,:, hvy_active(k)), g, x0, dx,&
                  hvy_tmp(:,:,:,:,hvy_active(k)),surface)
          endif
        enddo

    end if

end subroutine filter_wrapper
