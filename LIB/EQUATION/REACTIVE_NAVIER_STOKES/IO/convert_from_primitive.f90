!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name convert_from_primitive.f90
!> \version 0.5
!> \author msr
!
!> \brief convert to skew symmetric state vector from primitive state vector
!>
!! input:    - params_physics, state vector \n
!! output:   - work array \n
!!
!!
!! = log ======================================================================================
!! \n
!! 12/08/19 - create
!
! ********************************************************************************************

subroutine convert_from_primitive( params_physics, phi, phi_work )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> parameter struct
    type(type_params_rns), intent(inout)    :: params_physics

    !> state vector for one block, work array
    real(kind=rk), intent(inout)            :: phi(:, :, :, :), phi_work(:, :, :, :)

    ! field indexes
    integer(kind=ik)                        :: rhoF, UxF, UyF, UzF, EF, YF

    ! loop variables
    integer(kind=ik)                        :: k

!---------------------------------------------------------------------------------------------
! variables initialization

    ! field indexes from params
    rhoF = params_physics%rhoF
    UxF  = params_physics%UxF
    UyF  = params_physics%UyF
    UzF  = params_physics%UzF
    EF   = params_physics%EF
    YF   = params_physics%YF

!---------------------------------------------------------------------------------------------
! main body

    ! species
    !-----------------------------------------------------------------------------------------
    ! check chemistry model: 
    select case(params_physics%chemistry_model)

        case('inert')
            ! nothing to do

        case('one_step')
            phi_work(:, :, :, YF) = phi(:,:,:,YF) * phi_work(:, :, :, rhoF)

        case('cantera')
            ! compute inert specie (assume specie is stored in last datafield), to ensure mass conservation
            phi_work(:, :, :, YF+params_physics%species-1) = 1.0_rk
            do k = 1, params_physics%species-1
                 phi_work(:, :, :, YF+k-1) = phi(:,:,:,YF+k-1) * phi_work(:, :, :, rhoF)
                 phi_work(:, :, :, YF+params_physics%species-1) = phi_work(:, :, :, YF+params_physics%species-1) - phi_work(:, :, :, YF+k-1)
            end do

        case default
            call abort(120819002,"ERROR: can not convert to primitive variables, no chemistry model specified")

    end select

    ! energy field
    !-----------------------------------------------------------------------------------------
    ! check chemistry model: 
    ! inert and one step chemistry use pressure in energy equation
    ! cantera chemistry models use sensible internal energy in energy equation
    select case(params_physics%chemistry_model)

        case('inert', 'onestep')
            phi_work(:, :, :, EF)   = phi(:,:,:,EF)

        case('cantera')
            phi_work(:, :, :, EF)   = phi(:,:,:,EF) * phi_work(:, :, :, rhoF)

        case default
            call abort(120819003,"ERROR: can not convert to primitive variables, no chemistry model specified")

    end select

    ! rho
    !-----------------------------------------------------------------------------------------
    phi_work(:, :, :, rhoF) = dsqrt(phi(:,:,:,rhoF))

    ! Ux
    !-----------------------------------------------------------------------------------------
    phi_work(:, :, :, UxF)  = phi(:,:,:,UxF) * phi_work(:, :, :, rhoF)

    ! Uy
    !-----------------------------------------------------------------------------------------
    phi_work(:, :, :, UyF)  = phi(:,:,:,UyF) * phi_work(:, :, :, rhoF) 

    ! Uz
    !-----------------------------------------------------------------------------------------
    ! 3D?
    if ( params_physics%d == 3 ) then
        phi_work(:, :, :, UzF) = phi(:,:,:,UzF) * phi_work(:, :, :, rhoF)
    end if

end subroutine convert_from_primitive
