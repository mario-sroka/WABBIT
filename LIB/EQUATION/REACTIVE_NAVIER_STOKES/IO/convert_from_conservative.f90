!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name convert_from_conservative.f90
!> \version 0.5
!> \author msr
!
!> \brief convert to skew symmetric state vector from conservative state vector
!> q_conservative = (rho, rho*u, rho*v, rho*w, p)
!>
!! input:    - params_physics, state vector \n
!! output:   - work array \n
!!
!!
!! = log ======================================================================================
!! \n
!! 09/08/19 - create
!
! ********************************************************************************************

subroutine convert_from_conservative( params_physics, phi, phi_work )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> parameter struct
    type(type_params_rns), intent(inout)    :: params_physics

    !> state vector for one block, work array
    real(kind=rk), intent(inout)            :: phi(:, :, :, :), phi_work(:, :, :, :)

    ! field indexes
    integer(kind=ik)                        :: rhoF, UxF, UyF, UzF, EF

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

!---------------------------------------------------------------------------------------------
! main body

    ! rho
    !-----------------------------------------------------------------------------------------
    phi_work(:, :, :, rhoF) = dsqrt(phi(:, :, :, rhoF))

    ! Ux
    !-----------------------------------------------------------------------------------------
    phi_work(:, :, :, UxF)  = phi(:, :, :, UxF) / phi_work(:, :, :, rhoF)

    ! Uy
    !-----------------------------------------------------------------------------------------
    phi_work(:, :, :, UyF)  = phi(:, :, :, UyF) / phi_work(:, :, :, rhoF) 

    ! Uz
    !-----------------------------------------------------------------------------------------
    ! 3D?
    if ( params_physics%d == 3 ) then
        phi_work(:, :, :, UzF) = phi(:, :, :, UzF) / phi_work(:, :, :, rhoF)
    end if

    ! energy field
    !-----------------------------------------------------------------------------------------
    ! check chemistry model: 
    ! inert and one step chemistry use pressure in energy equation
    ! cantera chemistry models use sensible internal energy in energy equation
    select case(params_physics%chemistry_model)

        case('inert', 'onestep')
            ! kinetic energy
            if ( params_physics%d == 3 ) then
                phi_work(:, :, :, EF)   = phi_work(:,:,:,UxF)**2.0_rk + phi_work(:,:,:,UyF)**2.0_rk + phi_work(:,:,:,UzF)**2.0_rk
            else
                phi_work(:, :, :, EF)   = phi_work(:,:,:,UxF)**2.0_rk + phi_work(:,:,:,UyF)**2.0_rk
            end if
            phi_work(:, :, :, EF) = phi_work(:, :, :, EF) * 0.5_rk

            ! p = (gamma-1)*(e_tot-e_kin)
            phi_work(:, :, :, EF) = (phi(:, :, :, EF) - phi_work(:, :, :, EF)) * (params_physics%gamma_ - 1.0_rk)

        case default
            call abort(080620003,"ERROR: can not convert from conservative variables, wrong chemistry model specified")

    end select

end subroutine convert_from_conservative
