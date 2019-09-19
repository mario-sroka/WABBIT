!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name filter_reactive_ns.f90
!> \version 0.5
!> \author msr
!
!> \brief filter data
!>
!! input:    - params_physics, data for one block, spacing parameter, work_array \n
!! output:   - filtered data \n
!!
!!
!! = log ======================================================================================
!! \n
!! 10/10/18 - create
!
! ********************************************************************************************

subroutine filter_reactive_ns( params_physics, phi, phi_work, x0, dx, stage, filter_type )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> parameter struct
    type(type_params_rns), intent(inout)    :: params_physics

    !> state vector for one block, work array
    real(kind=rk), intent(inout)            :: phi(:, :, :, :), phi_work(:, :, :, :)

    !> spacing and origin of block
    real(kind=rk), intent(in)               :: x0(1:3), dx(1:3)

    !> filter stage: init, filter or post
    character(len=*), intent(in)            :: stage

    !> filter type, needed 
    character(len=*), intent(in)            :: filter_type

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    select case (filter_type)

        case('explicit_3pt', 'explicit_5pt', 'explicit_7pt', 'explicit_9pt', 'explicit_11pt')
            ! no init or post stage
            if ( stage=='filter_stage' ) call explicit_filter( params_physics, phi, phi_work, x0, dx )

        case('spectral')
            call spectral_filter( params_physics, phi, phi_work, x0, dx, stage )

        case('no_filter')
            ! nothing to do

        case default
            call abort(101018003," ERROR: in filter reactive_ns, filter "// params_physics%filter_type //" is unknown.")

    end select

end subroutine filter_reactive_ns
