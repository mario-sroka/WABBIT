!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name convert_to_conservative.f90
!> \version 0.5
!> \author msr
!
!> \brief convert skew symmetric state vector to conservative state vector
!> q = (rho, rho*u, rho*v, rho*w, p)
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

subroutine convert_to_conservative( params_physics, phi, phi_work )

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
    phi_work(:, :, :, rhoF) = phi(:,:,:,rhoF) * phi(:,:,:,rhoF)

    ! Ux
    !-----------------------------------------------------------------------------------------
    phi_work(:, :, :, UxF)  = phi(:,:,:,UxF) * phi(:,:,:,rhoF)

    ! Uy
    !-----------------------------------------------------------------------------------------
    phi_work(:, :, :, UyF)  = phi(:,:,:,UyF) * phi(:,:,:,rhoF) 

    ! Uz
    !-----------------------------------------------------------------------------------------
    ! 3D?
    if ( params_physics%d == 3 ) then
        phi_work(:, :, :, UzF) = phi(:,:,:,UzF) * phi(:,:,:,rhoF)
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
                phi_work(:, :, :, EF)   = phi(:,:,:,UxF)**2.0_rk + phi(:,:,:,UyF)**2.0_rk + phi(:,:,:,UzF)**2.0_rk
            else
                phi_work(:, :, :, EF)   = phi(:,:,:,UxF)**2.0_rk + phi(:,:,:,UyF)**2.0_rk
            end if
            phi_work(:, :, :, EF) = phi_work(:, :, :, EF) * 0.5_rk

            ! total energy
            ! e_tot=e_kin+p/(gamma-1)
            phi_work(:, :, :, EF)   = phi_work(:, :, :, EF) + phi(:,:,:,EF)/(params_physics%gamma_ - 1.0_rk)

        case default
            call abort(080620002,"ERROR: can not convert to conservative variables, wrong chemistry model specified")

    end select

end subroutine convert_to_conservative
