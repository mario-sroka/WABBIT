!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name RHS_wrapper_reactive_ns.f90
!> \version 0.5
!> \author msr
!
!> \brief wrapper for RHS of reactive navier stokes module
!>
!! input:    - params_physics, data for one block, spacing parameter \n
!! output:   - computed rhs \n
!!
!!
!! = log ======================================================================================
!! \n
!! 10/10/18 - create
!
! ********************************************************************************************

subroutine RHS_wrapper_reactive_ns( params_physics, time, phi, phi_work, x0, dx, rhs, stage, gas )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> parameter struct
    type(type_params_rns), intent(inout)    :: params_physics

    ! it may happen that some source terms have an explicit time-dependency
    ! therefore the general call has to pass time
    real(kind=rk), intent (in)              :: time
    !> state vector for one block
    real(kind=rk), intent(inout)            :: phi(:, :, :, :), phi_work(:, :, :, :)
    !> spacing and origin of block
    real(kind=rk), intent(in)               :: x0(1:3), dx(1:3)

    !> rhs output
    real(kind=rk), intent(inout)            :: rhs(:, :, :, :)

    ! stage. there is 3 stages, init_stage, integral_stage and local_stage. If the PDE has
    ! terms that depend on global qtys, such as forces etc, which cannot be computed
    ! from a single block alone, the first stage does that. the second stage can then
    ! use these integral qtys for the actual RHS evaluation.
    character(len=*), intent(in)            :: stage

    !> Cantera gas mixture struct
    type(phase_t), intent(inout)            :: gas

    ! periodic or non-periodic boundary 
    logical                                 :: periodic

    ! readability
    integer(kind=ik)                        :: g, Bs(3), NdF, Ds(3)

!---------------------------------------------------------------------------------------------
! variables initialization

    periodic = params_physics%periodic_BC(1) .and. params_physics%periodic_BC(2) .and. params_physics%periodic_BC(3)

    ! readability
    g   = params_physics%g
    Bs  = params_physics%Bs
    Ds  = params_physics%Bs + 2*params_physics%g
    NdF = params_physics%NdF

!---------------------------------------------------------------------------------------------
! main body

    select case(stage)
        case ("init_stage")
        !-------------------------------------------------------------------------
        ! 1st stage: init_stage.
        !-------------------------------------------------------------------------
        ! this stage is called only once, not for each block.
        ! performs initializations in the RHS module, such as resetting integrals

    case ("integral_stage")
        !-------------------------------------------------------------------------
        ! 2nd stage: init_stage.
        !-------------------------------------------------------------------------
        ! For some RHS, the eqn depend not only on local, block based qtys, such as
        ! the state vector, but also on the entire grid, for example to compute a
        ! global forcing term (e.g. in FSI the forces on bodies). As the physics
        ! modules cannot see the grid, (they only see blocks), in order to encapsulate
        ! them nicer, two RHS stages have to be defined: integral / local stage.

    case ("post_stage")
        !-------------------------------------------------------------------------
        ! 3rd stage: post_stage.
        !-------------------------------------------------------------------------
        ! this stage is called only once, not for each block.

    case ("local_stage")
        !-------------------------------------------------------------------------
        ! 4th stage: local evaluation of RHS on all blocks
        !-------------------------------------------------------------------------
        ! the second stage then is what you would usually do: evaluate local differential
        ! operators etc.

        if     ( params_physics%chemistry_model == "cantera" ) then

                 !-------------------------------------------------------------------------
                 if (params_physics%d == 2) then
                     
                     if (periodic) then
                         ! periodic RHS
                         call RHS_2D_CANTERA_navier_stokes_reactive_periodicBC( params_physics, Bs, &
                         g, NdF, x0, dx, phi(:,:,:,:), rhs(:,:,:,:), gas )
                     else
                         ! /todo
                     end if

                 else

                     if (periodic) then
                         ! periodic RHS
                         call RHS_3D_CANTERA_navier_stokes_reactive_periodicBC( params_physics, Bs, &
                         g, NdF, x0, dx, phi(:,:,:,:), rhs(:,:,:,:), gas )
                     else
                         ! non-periodic RHS
                         call RHS_3D_CANTERA_navier_stokes_reactive_non_periodicBC( params_physics, Bs, &
                         g, NdF, x0, dx, phi(:,:,:,:), rhs(:,:,:,:), gas )
                     end if
                
                 end if

        elseif ( params_physics%chemistry_model == "inert" .and. &
                 periodic ) then

                 !-------------------------------------------------------------------------
                 if (params_physics%d == 2) then
                     ! /todo
                 else
                     call RHS_3D_navier_stokes_non_reactive_periodicBC(params_physics, Ds, &
                     Bs, g, NdF, x0, dx, phi(:,:,:,:), rhs(:,:,:,:), time)

                 end if

        else
                !-------------------------------------------------------------------------
                ! error case
                call abort(101018001,"reactive ns RHS wrapper: can not find proper RHS subroutine.")

        end if

    case default
        call abort(101018002,"the RHS wrapper requests a stage this physics module cannot handle.")

    end select

end subroutine RHS_wrapper_reactive_ns
