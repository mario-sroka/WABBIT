!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name get_dt_reactive_ns.f90
!> \version 0.5
!> \author msr
!
!> \brief compute dt
!>
!! input:    - params_physics, data for one block, spacing parameter \n
!! output:   - \n
!!
!!
!! = log ======================================================================================
!! \n
!! 28/03/19 - create
!
! ********************************************************************************************

subroutine get_dt_reactive_ns( params_physics, time, phi, x0, dx, dt )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> parameter struct
    type(type_params_rns), intent(inout)    :: params_physics

    ! it may happen that some source terms have an explicit time-dependency
    ! therefore the general call has to pass time
    real(kind=rk), intent (in)              :: time
    !> state vector for one block
    real(kind=rk), intent(inout)            :: phi(:, :, :, :)
    !> spacing and origin of block
    real(kind=rk), intent(in)               :: x0(1:3), dx(1:3)

    ! the dt for this block is returned to the caller:
    real(kind=rk), intent(inout)            :: dt

    ! grid parameter
    integer(kind=ik)                        :: Bs(3), g
    ! loop parameter
    integer(kind=ik)                        :: i, j, k
    ! field indexes
    integer(kind=ik)                        :: rhoF, UxF, UyF, UzF, EF, YF
    ! dummy variables
    real(kind=rk)                           :: dx_min
    ! physical velocity
    real(kind=rk), allocatable, save        :: u_physical(:,:,:)

!---------------------------------------------------------------------------------------------
! variables initialization

    Bs = params_physics%Bs
    g  = params_physics%g

    ! field indexes from params
    rhoF = params_physics%rhoF
    UxF  = params_physics%UxF
    UyF  = params_physics%UyF
    UzF  = params_physics%UzF
    EF   = params_physics%EF
    YF   = params_physics%YF

    ! temp array
    if ( params_physics%d == 2 ) then
        if ( .not. allocated(u_physical) ) allocate(u_physical(Bs(1)+2*g, Bs(2)+2*g, 1))
        u_physical = phi(:,:,:,UxF) * phi(:,:,:,UxF) + phi(:,:,:,UyF) * phi(:,:,:,UyF)
    else
        if ( .not. allocated(u_physical) ) allocate(u_physical(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g))
        u_physical = phi(:,:,:,UxF) * phi(:,:,:,UxF) + phi(:,:,:,UyF) * phi(:,:,:,UyF) + phi(:,:,:,UzF) * phi(:,:,:,UzF)
    end if

    ! init dt
    dt = 9.9e9_rk
    dx_min=minval(dx(1:params_physics%d))

!---------------------------------------------------------------------------------------------
! main body

    if (params_physics%chemistry_model /= 'inert') then
        ! CFL criteria not working for reactive computations
        call abort(120819001,"ERROR: current chemistry model should not use adaptive time stepping.")

    else
        ! assume skew symmetric RHS, so EF is actually pressure p
        u_physical = dsqrt(u_physical) + dsqrt(params_physics%gamma_ * phi(:,:,:,EF)) 
        u_physical = u_physical/phi(:,:,:,rhoF)     

        ! CFL criteria CFL=u_physical/u_numerical where u_numerical=dx/dt
        dt = min(dt, params_physics%CFL * dx_min / maxval(u_physical))

    end if

end subroutine get_dt_reactive_ns
