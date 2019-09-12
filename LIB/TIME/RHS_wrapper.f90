
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \file
!> \callgraph
!> \brief wrapper for RHS call in time step function, computes RHS in work array
!! (inplace)
!> \version 0.5
!> \author sm
!! \date 23/05/17 - create
!!
!
!>\details
!! calls RHS depending on physics
!!
!! butcher table, e.g.
!!
!! |   |    |    |   |
!! |---|----|----|---|
!! | 0 | 0  | 0  |  0|
!! |c2 | a21| 0  |  0|
!! |c3 | a31| a32|  0|
!! | 0 | b1 | b2 | b3|
!!
!!
!
!**********************************************************************************************

subroutine RHS_wrapper(time, params, hvy_block, hvy_rhs, hvy_mask, hvy_tmp, lgt_block, &
    lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n, hvy_neighbor)
   implicit none

    !> time variable
    real(kind=rk), intent(in)           :: time
    !> user defined parameter structure, hvy_active
    type (type_params), intent(in)      :: params
    !> heavy work data array - block data
    real(kind=rk), intent(inout)        :: hvy_rhs(:, :, :, :, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> hvy_mask are qtys that depend on grid and not explicitly on time
    real(kind=rk), intent(inout)        :: hvy_mask(:, :, :, :, :)
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    !> light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(inout)     :: hvy_active(:,:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(inout)     :: hvy_n(:)
    !> list of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_active(:,:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_n(:)
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), intent(inout)  :: lgt_sortednumlist(:,:,:)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(inout)     :: hvy_neighbor(:,:)

    !> global integral
    real(kind=rk), dimension(3)         :: volume_int
    !> spacing and origin of a block
    real(kind=rk), dimension(3)         :: dx, x0
    ! loop variables
    integer(kind=ik)                    :: k, dF, neqn, lgt_id, hvy_id
    ! grid parameter, error variable
    integer(kind=ik)                    :: g
    integer(kind=ik), dimension(3)      :: Bs
    integer(kind=2)                     :: surface(3)=0

    ! grid parameter
    Bs = params%Bs
    g  = params%n_ghosts

    !-------------------------------------------------------------------------
    ! create mask function at current time
    !-------------------------------------------------------------------------
    call create_mask_tree(params, time, lgt_block, hvy_mask, hvy_tmp, &
        hvy_neighbor, hvy_active, hvy_n, lgt_active, lgt_n, lgt_sortednumlist)


    !-------------------------------------------------------------------------
    ! 1st stage: init_stage. (called once, not for all blocks)
    !-------------------------------------------------------------------------
    ! performs initializations in the RHS module, such as resetting integrals
    hvy_id = hvy_active(1, tree_ID_flow) ! for this stage, just pass any block (result does not depend on block)
    call RHS_meta( params%physics_type, time, hvy_block(:,:,:,:,hvy_id), g, x0, dx, &
         hvy_rhs(:,:,:,:,hvy_id), hvy_mask(:,:,:,:,hvy_id), "init_stage" )

    !-------------------------------------------------------------------------
    ! 2nd stage: integral_stage. (called for all blocks)
    !-------------------------------------------------------------------------
    ! For some RHS, the eqn depend not only on local, block based qtys, such as
    ! the state vector, but also on the entire grid, for example to compute a
    ! global forcing term (e.g. in FSI the forces on bodies). As the physics
    ! modules cannot see the grid, (they only see blocks), in order to encapsulate
    ! them nicer, two RHS stages have to be defined: integral / local stage.
    do k = 1, hvy_n(tree_ID_flow)
        hvy_id = hvy_active(k, tree_ID_flow)
        ! convert given hvy_id to lgt_id for block spacing routine
        call hvy_id_to_lgt_id( lgt_id, hvy_id, params%rank, params%number_blocks )
        ! get block spacing for RHS
        call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )
        if ( .not. All(params%periodic_BC) ) then
            ! check if block is adjacent to a boundary of the domain, if this is the case we use one sided stencils
            call get_adjacent_boundary_surface_normal(params, lgt_id, lgt_block, params%max_treelevel, surface)
        endif
        call RHS_meta( params%physics_type, time, hvy_block(:,:,:,:, hvy_id), g, x0, dx,&
        hvy_rhs(:,:,:,:,hvy_id), hvy_mask(:,:,:,:,hvy_id), "integral_stage", boundary_flag=surface )
    enddo


    !-------------------------------------------------------------------------
    ! 3rd stage: post integral stage. (called once, not for all blocks)
    !-------------------------------------------------------------------------
    ! in rhs module, used ror example for MPI_REDUCES
    hvy_id = hvy_active(1, tree_ID_flow) ! for this stage, just pass any block (result does not depend on block)
    call RHS_meta( params%physics_type, time, hvy_block(:,:,:,:, hvy_id), g, x0, dx, &
    hvy_rhs(:,:,:,:,hvy_id), hvy_mask(:,:,:,:,hvy_id), "post_stage" )


    !-------------------------------------------------------------------------
    ! 4th stage: local evaluation of RHS on all blocks (called for all blocks)
    !-------------------------------------------------------------------------
    ! the second stage then is what you would usually do: evaluate local differential
    ! operators etc.

    ! I think the idea would be to add mask generation here; take it out of
    ! physics modules rhs.f90

    do k = 1, hvy_n(tree_ID_flow)
        hvy_id = hvy_active(k, tree_ID_flow)
        ! convert given hvy_id to lgt_id for block spacing routine
        call hvy_id_to_lgt_id( lgt_id, hvy_id, params%rank, params%number_blocks )
        ! get block spacing for RHS
        call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )
        if ( .not. All(params%periodic_BC) ) then
            ! check if block is adjacent to a boundary of the domain, if this is the case we use one sided stencils
            call get_adjacent_boundary_surface_normal(params, lgt_id, lgt_block, params%max_treelevel, surface)
        endif

        call RHS_meta( params%physics_type, time, hvy_block(:,:,:,:, hvy_id), g, &
        x0, dx, hvy_rhs(:,:,:,:, hvy_id), hvy_mask(:,:,:,:, hvy_id), "local_stage", &
        boundary_flag=surface )
    enddo

end subroutine RHS_wrapper
