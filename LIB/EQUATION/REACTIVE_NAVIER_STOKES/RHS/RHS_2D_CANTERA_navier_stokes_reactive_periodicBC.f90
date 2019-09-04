!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name RHS_2D_CANTERA_navier_stokes_reactive_periodicBC.f90
!> \version 0.5
!> \author msr
!
!> \brief RHS for 2D navier stokes equation with cantera chemistry and periodic boundaries
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

subroutine RHS_2D_CANTERA_navier_stokes_reactive_periodicBC(params_physics, Bs, g, NdF, x0, delta_x, phi, rhs, gas)

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
    real(kind=rk), intent(inout)                            :: phi(Bs(1)+2*g, Bs(2)+2*g, 1, NdF)
    ! rhs array
    real(kind=rk),intent(inout)                             :: rhs(Bs(1)+2*g, Bs(2)+2*g, 1, NdF)
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
    real(kind=rk)                                           :: rho(Bs(1)+2*g, Bs(2)+2*g), u(Bs(1)+2*g, Bs(2)+2*g), v(Bs(1)+2*g, Bs(2)+2*g), &
                                                               w(Bs(1)+2*g, Bs(2)+2*g), p(Bs(1)+2*g, Bs(2)+2*g), T(Bs(1)+2*g, Bs(2)+2*g), &
                                                               tau11(Bs(1)+2*g, Bs(2)+2*g), tau22(Bs(1)+2*g, Bs(2)+2*g), tau33(Bs(1)+2*g, Bs(2)+2*g), &
                                                               tau12(Bs(1)+2*g, Bs(2)+2*g), tau13(Bs(1)+2*g, Bs(2)+2*g), tau23(Bs(1)+2*g, Bs(2)+2*g), &
                                                               mu(Bs(1)+2*g, Bs(2)+2*g), lambda(Bs(1)+2*g, Bs(2)+2*g), &
                                                               es(Bs(1)+2*g, Bs(2)+2*g), Y(Bs(1)+2*g, Bs(2)+2*g, params_physics%species), &
                                                               X(Bs(1)+2*g, Bs(2)+2*g, params_physics%species), &
                                                               D(Bs(1)+2*g, Bs(2)+2*g, params_physics%species)

    ! dummy field
    real(kind=rk)                                           :: dummy(Bs(1)+2*g, Bs(2)+2*g), dummy2(Bs(1)+2*g, Bs(2)+2*g), dummy3(Bs(1)+2*g, Bs(2)+2*g), &
                                                               dummy4(Bs(1)+2*g, Bs(2)+2*g), q(params_physics%species)


    ! inverse sqrt(rho) field 
    real(kind=rk)                                           :: phi1_inv(Bs(1)+2*g, Bs(2)+2*g)

    ! loop variables
    integer(kind=ik)                                        :: i, j, k, n

    ! chemistry parameter
    ! reaction heat
    real(kind=rk)                                           :: dh(params_physics%species)

    ! field indexes
    integer(kind=ik)                                        :: rhoF, UxF, UyF, EF, YF

    ! switch to enabled one sided derivatives
    logical                                                 :: onesided(2,3)

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

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
    EF   = params_physics%EF
    YF   = params_physics%YF

    ! set physics parameters for readability
    dissipation = params_physics%dissipation

    ! primitive variables
    ! use rhs as dummy fields
    call convert_to_primitive( params_physics, phi, rhs )
 
    rho = rhs(:,:,1,rhoF)
    u   = rhs(:,:,1,UxF)
    v   = rhs(:,:,1,UyF)
    es  = rhs(:,:,1,EF)
    do n = 1, params_physics%species
        Y(:,:,n) = rhs(:,:,1,YF+n-1)
    end do

    ! compute 1/rho for better performance
    phi1_inv(:,:)  = 1.0_rk / phi(:,:,1,rhoF)

    ! discretization constant
    dx = delta_x(1)
    dy = delta_x(2)
    dz = delta_x(3)

    ! chemistry
    dh  = params_physics%dh

    !#########################################################################################
    ! CANTERA common code
    include "RHS/RHS_2D_CANTERA_common_code.f90"
    !#########################################################################################

!---------------------------------------------------------------------------------------------
! main body

    !#########################################################################################
    ! inlcude RHS for periodic BC
    include "RHS/RHS_2D_common_code.f90"
    !#########################################################################################

    !#########################################################################################
    ! CANTERA common code
    include "RHS/RHS_2D_CANTERA_common_code_equations.f90"
    !#########################################################################################    

    !#########################################################################################
    ! penalization -> nothing to do here
    !#########################################################################################

end subroutine RHS_2D_CANTERA_navier_stokes_reactive_periodicBC
