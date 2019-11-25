!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name RHS_3D_CANTERA_navier_stokes_reactive_periodicBC.f90
!> \version 0.5
!> \author msr
!
!> \brief RHS for 3D navier stokes equation with cantera chemistry and periodic boundaries
!
!>
!! input:    - datafields, grid parameter, chemistry mixture variable \n
!! output:   - RHS(datafield) \n
!!
!!
!! = log ======================================================================================
!! \n
!! 14/02/19 - create
!
! ********************************************************************************************

subroutine RHS_3D_CANTERA_navier_stokes_reactive_periodicBC(params_physics, Bs, g, NdF, x0, delta_x, phi, rhs, gas)

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> navier stokes params struct, note: contails all parameters needed by RHS
    type(type_params_rns), intent(inout)                    :: params_physics
    !> grid parameter
    integer(kind=ik), intent(inout)                         :: g, Bs(3), NdF
    !> rhs parameter
    real(kind=rk), dimension(3), intent(in)                 :: x0, delta_x
    !> datafields
    real(kind=rk), intent(inout)                            :: phi(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g, NdF)
    ! rhs array
    real(kind=rk),intent(inout)                             :: rhs(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g, NdF)
    !> Cantera gas mixture struct
    type(phase_t), intent(inout)                            :: gas

    ! specific gas constant
    real(kind=rk)                                           :: Rs, Wk, Wk_inv, Wk_species_inv(params_physics%species)
    ! dynamic viscosity
    real(kind=rk)                                           :: mu_d, two_three, mu0

    ! spacing
    real(kind=rk)                                           :: dx, dy, dz
    ! dissipation switch
    logical                                                 :: dissipation

    ! variables
    real(kind=rk)                                           :: rho(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), u(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), v(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), &
                                                               w(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), p(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), T(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), &
                                                               tau11(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), tau22(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), tau33(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), &
                                                               tau12(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), tau13(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), tau23(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), &
                                                               mu(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), lambda(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), &
                                                               es(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), Y(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g, params_physics%species), &
                                                               X(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g, params_physics%species), &
                                                               D(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g, params_physics%species)

    ! dummy field
    real(kind=rk)                                           :: dummy(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), dummy2(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), dummy3(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), &
                                                               dummy4(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), dummy5(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), dummy6(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), &
                                                               q(params_physics%species)


    ! inverse sqrt(rho) field 
    real(kind=rk)                                           :: phi1_inv(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

    ! loop variables
    integer(kind=ik)                                        :: i, j, k, n

    ! array indexes
    integer(kind=ik)                                        :: gp, gm, Bp(3), Bm(3)

    ! chemistry parameter
    ! reaction heat
    real(kind=rk)                                           :: dh(params_physics%species)

    ! field indexes
    integer(kind=ik)                                        :: rhoF, UxF, UyF, UzF, EF, YF

    ! switch to enabled one sided derivatives
    logical                                                 :: onesided(2,3)

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! compute array indexes
    gp   = g+1
    gm   = g-1
    Bp   = Bs(:) + g
    Bm   = Bs(:) + g + 2

    ! periodic boundary here, so:
    onesided = .false.

    ! inverse parameters
    do n = 1, params_physics%species
        Wk_species_inv(n) = 1.0_rk / params_physics%Wk(n)
    end do
    two_three = 2.0_rk / 3.0_rk

    ! field indexes from params
    rhoF = params_physics%rhoF
    UxF  = params_physics%UxF
    UyF  = params_physics%UyF
    UzF  = params_physics%UzF
    EF   = params_physics%EF
    YF   = params_physics%YF

    ! set physics parameters for readability
    dissipation = params_physics%dissipation

    ! compute 1/rho for better performance
    phi1_inv(:,:,:)  = 1.0_rk / phi(:,:,:,rhoF)

    ! primitive variables
    ! use rhs as dummy fields
    call convert_to_primitive( params_physics, phi, rhs, phi1_inv )
 
    rho = rhs(:,:,:,rhoF)
    u   = rhs(:,:,:,UxF)
    v   = rhs(:,:,:,UyF)
    w   = rhs(:,:,:,UzF)
    es  = rhs(:,:,:,EF)
    do n = 1, params_physics%species
        Y(:,:,:,n) = rhs(:,:,:,YF+n-1)
    end do

    ! discretization constant
    dx = delta_x(1)
    dy = delta_x(2)
    dz = delta_x(3)

    ! chemistry
    dh  = params_physics%dh

    !#########################################################################################
    ! CANTERA common code
    include "RHS/RHS_3D_CANTERA_common_code.f90"
    !#########################################################################################

!---------------------------------------------------------------------------------------------
! main body

    !#########################################################################################
    ! inlcude RHS for periodic BC
    include "RHS/RHS_3D_common_code.f90"
    !#########################################################################################

    !#########################################################################################
    ! CANTERA common code
    include "RHS/RHS_3D_CANTERA_common_code_equations.f90"
    !#########################################################################################    

    !#########################################################################################
    ! forcing
    !#########################################################################################
    if ( params_physics%forcing ) then
        ! /todo
    end if

    !#########################################################################################
    ! penalization -> nothing to do here
    !#########################################################################################

end subroutine RHS_3D_CANTERA_navier_stokes_reactive_periodicBC
