!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name prepare_saved_data_reactive_ns.f90
!> \version 0.5
!> \author msr
!
!> \brief write data to work array
!>
!! input:    - params_physics, state vector \n
!! output:   - work array \n
!!
!!
!! = log ======================================================================================
!! \n
!! 10/10/18 - create
!
! ********************************************************************************************

subroutine prepare_saved_data_reactive_ns( params_physics, phi, phi_work, time, dx, x0, gas )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> parameter struct
    type(type_params_rns), intent(inout)    :: params_physics

    !> state vector for one block, work array
    real(kind=rk), intent(inout)            :: phi(:, :, :, :), phi_work(:, :, :, :)

    ! time of simulation
    real(kind=rk), intent (in)              :: time

    !> spacing and origin of block
    real(kind=rk), intent(in)               :: dx(1:3), x0(1:3)

    !> Cantera gas mixture struct
    type(phase_t), intent(inout)            :: gas

    ! field indexes
    integer(kind=ik)                        :: rhoF, UxF, UyF, UzF, EF, YF
    ! loop variables
    integer(kind=ik)                        :: k !, l, i, j, n
    ! temp array
    real(kind=rk), allocatable, save        :: tmp(:,:,:,:)
    ! grid parameter
    integer(kind=ik)                        :: Bs(3), g

!---------------------------------------------------------------------------------------------
! variables initialization

    ! field indexes from params
    rhoF = params_physics%rhoF
    UxF  = params_physics%UxF
    UyF  = params_physics%UyF
    UzF  = params_physics%UzF
    EF   = params_physics%EF
    YF   = params_physics%YF

    ! grid parameter
    Bs = params_physics%Bs
    g  = params_physics%g

    ! allocate temp array
    if ( .not. allocated(tmp) ) allocate(tmp(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g, 6))

!---------------------------------------------------------------------------------------------
! main body

    !-----------------------------------------------------------------------------------------
    ! check if primitive variable IO is required
    !-----------------------------------------------------------------------------------------
    if (params_physics%save_primitive) then

        ! convert from skew symmetric form
        call convert_to_primitive(params_physics, phi, phi_work(:, :, :, 1:params_physics%NdF) )

    else
        ! nothing to do
        phi_work(:, :, :, 1:params_physics%NdF) = phi

    end if
    
    !-----------------------------------------------------------------------------------------
    ! other fields
    !-----------------------------------------------------------------------------------------

    do k = params_physics%NdF+1, params_physics%N_fields_saved

        select case(params_physics%names_saved(k))
            !---------------------------------------------------------------------------------
            case("T")

                if ( params_physics%chemistry_model == 'cantera' ) then

                    ! cantera chemistry
                    call compute_temperature( params_physics, phi, phi_work(:,:,:,k), gas )

                else

                    ! non cantera chemistry
                    phi_work(:, :, :, k) = phi(:,:,:,EF) / phi(:,:,:,rhoF)**2.0_rk / params_physics%Rs

                end if

            !---------------------------------------------------------------------------------
            case("p")

                if ( params_physics%chemistry_model == 'cantera' ) then

                    ! cantera chemistry
                    call compute_pressure( params_physics, phi, phi_work(:,:,:,k), gas )

                else

                    ! non cantera chemistry
                    phi_work(:, :, :, k) = phi(:,:,:,EF)

                end if

            !---------------------------------------------------------------------------------
            case("wk")

                if ( params_physics%chemistry_model == 'cantera' ) then

                    ! cantera chemistry
                    call compute_reaction_rate( params_physics, phi, phi_work(:,:,:,k), gas )

                else

                    ! non cantera chemistry
                    phi_work(:, :, :, k) = phi(:,:,:,YF) * &
                    params_physics%A * dexp( -params_physics%Ta / &
                    (phi(:,:,:,EF) / phi(:,:,:,rhoF)**2.0_rk / params_physics%Rs) )

                end if

            !---------------------------------------------------------------------------------
            case("vort")
                if (params_physics%d == 3) then
                    call compute_vorticity(phi(:,:,:,UxF)/phi(:,:,:,rhoF), phi(:,:,:,UyF)/phi(:,:,:,rhoF), &
                                           phi(:,:,:,UzF)/phi(:,:,:,rhoF), &
                                           dx, Bs, g, &
                                           'FD_4th_central_optimized', tmp(:,:,:,1:3))

                    phi_work(:, :, :, k) = sqrt(tmp(:, :, :, 1)**2.0_rk + tmp(:, :, :, 2)**2.0_rk + tmp(:, :, :, 3)**2.0_rk)

                else
                    ! 2D case, third velocity component is a dummy argument
                    call compute_vorticity(phi(:,:,:,UxF)/phi(:,:,:,rhoF), phi(:,:,:,UyF)/phi(:,:,:,rhoF), &
                                           phi(:,:,:,UxF)/phi(:,:,:,rhoF), &
                                           dx, Bs, g, &
                                           'FD_4th_central_optimized', tmp(:,:,:,1:2))

                    phi_work(:, :, 1, k) = dsqrt( tmp(:, :, 1, 1)**2.0_rk + tmp(:, :, 1, 2)**2.0_rk )

                end if

            !---------------------------------------------------------------------------------
            case default
                call abort(121018004,"ERROR: have no data to save field: "//trim(adjustl(params_physics%names_saved(k))))

        end select

    end do

end subroutine prepare_saved_data_reactive_ns
