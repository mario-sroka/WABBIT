!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name ini_reactive_ns.f90
!> \version 0.5
!> \author msr
!
!> \brief initialize all data for reactive navier stokes
!>
!! input:    - params_physics, state vector \n
!! output:   - \n
!!
!!
!! = log ======================================================================================
!! \n
!! 09/10/18 - create
!
! ********************************************************************************************

subroutine ini_reactive_ns( params_physics, phi, x0, dx, gas )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> parameter struct
    type(type_params_rns), intent(inout)    :: params_physics

    !> state vector for one block
    real(kind=rk), intent(inout)            :: phi(:, :, :, :)

    !> spacing and origin of block
    real(kind=rk), intent(in)               :: x0(1:3), dx(1:3)

    !> Cantera gas mixture struct
    type(phase_t), intent(inout)            :: gas

    ! field indexes
    integer(kind=ik)                        :: rhoF, UxF, UyF, UzF, EF, YF
    ! loop variable
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

    ! case specific ini conditions
    select case (params_physics%inicond_name)

        case ("from_file")
            !---------------------------------------------------------------------------------------------
            ! compute skew symmetric variables if data from file were saved in primitive form
            ! note: if save_primitive is set, then we assume the data is stored in primitive from
            ! note: data is actually read in wabbit core, so if you wish to use this ini condition, 
            ! ensure read_from_files option is set accordingly
            if (params_physics%save_primitive) then                
                call  convert_from_primitive( params_physics, phi, phi )
            end if

        case ("cantera_spark")
            !---------------------------------------------------------------------------------------------
            ! check chemistry
            if (params_physics%chemistry_model /= 'cantera') then
                call abort(090819002,"ERROR: you can not use current chemistry model with cantera spark initial conditions!.")
            end if

            ! set velocity to zero
            call inicond_zero_velocity( params_physics, phi )

            ! insert spark into domain           
            call inicond_spark( params_physics, phi, x0, dx, gas )

        case ("taylor_green")
            !---------------------------------------------------------------------------------------------
            ! set velocity to taylor green vortices, density and pressure to inicond values
            call inicond_taylor_green( params_physics, phi, x0, dx )

        case ("pressure_blob")
            !---------------------------------------------------------------------------------------------
            ! check chemistry
            ! cantera physics need special initialization, because p is not element of the state vector
            if (params_physics%chemistry_model == 'cantera') then
                call abort(090819002,"ERROR: pressure blob is not implemented for cantera chemistry model.")
            end if

            ! set velocity to zero
            call inicond_zero_velocity( params_physics, phi )

            ! set density to initial density
            ! note: set up skew symmetric form here
            phi(:, :, :, rhoF) = dsqrt(params_physics%inicond_rho)

            ! set up density blob     
            ! note: a few parameters are given explicitly to the subroutine in order to have a more generic 
            ! call (at other positions in wabbit)      
            call inicond_blob( params_physics, phi(:, :, :, EF) , x0, dx, params_physics%d, &
                               params_physics%inicond_position(:), (/1.0_rk, params_physics%inicond_scales(2)/) )
            ! magnitude of blob is between [0,1]
            ! so: set up here pressure blob magnitude
            phi(:, :, :, EF) = params_physics%inicond_p + phi(:, :, :, EF) * params_physics%inicond_scales(1)

        case default
            call abort(091018002,"ERROR: unknown ini condition for reactive navier stokes.")

    end select

end subroutine ini_reactive_ns
