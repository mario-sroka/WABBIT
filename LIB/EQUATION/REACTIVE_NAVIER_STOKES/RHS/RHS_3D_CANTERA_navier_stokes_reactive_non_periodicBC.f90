!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name RHS_3D_CANTERA_navier_stokes_reactive_non_periodicBC.f90
!> \version 0.5
!> \author msr
!
!> \brief RHS for 3D navier stokes equation with cantera chemistry and non-periodic boundaries
!
!>
!! input:    - datafields, grid parameter, chemistry mixture variable \n
!! output:   - RHS(datafield) \n
!!
!!
!! = log ======================================================================================
!! \n
!! 30/08/19 - create
!
! ********************************************************************************************

subroutine RHS_3D_CANTERA_navier_stokes_reactive_non_periodicBC(params_physics, Bs, g, NdF, x0, delta_x, phi, rhs, gas)

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

    ! reference side for sponge computation
    integer(kind=ik)                                        :: ref_side

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! compute array indexes
    gp   = g+1
    gm   = g-1
    Bp   = Bs(:) + g
    Bm   = Bs(:) + g + 2

    ! set one sided boundary
    ! note: boundaries are set blockwise, so we need to calculate this here for every RHS call!
    call boundaries_xyz( onesided, x0, delta_x, params_physics%L, Bs, g, params_physics%periodic_BC )

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

        ! lundgren forcing in physical space
        ! pure solenoidal forcing, /todo maybe implement other forcing methods
        ! --------------------------------------------------------------------
        ! first: compute IDFT to get the fourier filtered velocity field
        call compute_IDFT( params_physics, params_physics%phi_hat_s(:,:,:,1), dummy(g+1:Bs(1)+g,  g+1:Bs(2)+g, g+1:Bs(3)+g), x0, delta_x )
        call compute_IDFT( params_physics, params_physics%phi_hat_s(:,:,:,2), dummy2(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g), x0, delta_x )
        call compute_IDFT( params_physics, params_physics%phi_hat_s(:,:,:,3), dummy3(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g), x0, delta_x )

        ! second: forcing strength
        ! use mu0 as dummy
        mu0 = - params_physics%target_force * ( params_physics%epsilon_s_mean - params_physics%eps_s_target ) / params_physics%eps_s_target

        ! third: add forcing to rhs
        rhs(g+1:Bs(1)+g,  g+1:Bs(2)+g, g+1:Bs(3)+g, UxF) = rhs(g+1:Bs(1)+g,  g+1:Bs(2)+g, g+1:Bs(3)+g, UxF) + mu0 * rho(g+1:Bs(1)+g,  g+1:Bs(2)+g, g+1:Bs(3)+g) &
                                                         * dummy(g+1:Bs(1)+g,  g+1:Bs(2)+g, g+1:Bs(3)+g)  * phi1_inv(g+1:Bs(1)+g,  g+1:Bs(2)+g, g+1:Bs(3)+g)
        rhs(g+1:Bs(1)+g,  g+1:Bs(2)+g, g+1:Bs(3)+g, UyF) = rhs(g+1:Bs(1)+g,  g+1:Bs(2)+g, g+1:Bs(3)+g, UyF) + mu0 * rho(g+1:Bs(1)+g,  g+1:Bs(2)+g, g+1:Bs(3)+g) &
                                                         * dummy2(g+1:Bs(1)+g,  g+1:Bs(2)+g, g+1:Bs(3)+g) * phi1_inv(g+1:Bs(1)+g,  g+1:Bs(2)+g, g+1:Bs(3)+g)
        rhs(g+1:Bs(1)+g,  g+1:Bs(2)+g, g+1:Bs(3)+g, UzF) = rhs(g+1:Bs(1)+g,  g+1:Bs(2)+g, g+1:Bs(3)+g, UzF) + mu0 * rho(g+1:Bs(1)+g,  g+1:Bs(2)+g, g+1:Bs(3)+g) &
                                                         * dummy3(g+1:Bs(1)+g,  g+1:Bs(2)+g, g+1:Bs(3)+g) * phi1_inv(g+1:Bs(1)+g,  g+1:Bs(2)+g, g+1:Bs(3)+g)

        rhs(g+1:Bs(1)+g,  g+1:Bs(2)+g, g+1:Bs(3)+g, EF)  = rhs(g+1:Bs(1)+g,  g+1:Bs(2)+g, g+1:Bs(3)+g, EF) - mu0 * rho(g+1:Bs(1)+g,  g+1:Bs(2)+g, g+1:Bs(3)+g) &
                                                         * ( u(g+1:Bs(1)+g,  g+1:Bs(2)+g, g+1:Bs(3)+g) * dummy(g+1:Bs(1)+g,  g+1:Bs(2)+g, g+1:Bs(3)+g)  &
                                                           + v(g+1:Bs(1)+g,  g+1:Bs(2)+g, g+1:Bs(3)+g) * dummy2(g+1:Bs(1)+g,  g+1:Bs(2)+g, g+1:Bs(3)+g) &
                                                           + w(g+1:Bs(1)+g,  g+1:Bs(2)+g, g+1:Bs(3)+g) * dummy3(g+1:Bs(1)+g,  g+1:Bs(2)+g, g+1:Bs(3)+g) )

    end if

    !#########################################################################################
    ! penalization
    !#########################################################################################
    ! set sponge field, use dummy array, note: /todo: store sponge field globally
    call set_penalization( params_physics, Bs, g, x0, delta_x, dummy )

    ! set sponge at non-periodic BCs
    ! /todo: better coding
    if ( .NOT.(params_physics%periodic_BC(1)) ) then
        ! compute boundary side: x+ or x-, note: split domain at domain center
        if ( x0(1) < (params_physics%L(1)/2.0_rk + 1e-12_rk) ) then
            ref_side = 1
        else
            ref_side = 2
        end if

        do k = g+1, Bs(3)+g
            do j = g+1, Bs(2)+g
                do i = g+1, Bs(1)+g

                    rhs(i,j,k,rhoF) = rhs(i,j,k,rhoF) - ( rho(i,j,k) - params_physics%rho_ref(ref_side) ) * dummy(i, j, k)

                    rhs(i,j,k,UxF)  = rhs(i,j,k,UxF)  - ( u(i,j,k)   - params_physics%u_ref(ref_side)   ) * dummy(i, j, k)

                    rhs(i,j,k,UyF)  = rhs(i,j,k,UyF)  - ( v(i,j,k)   - params_physics%v_ref(ref_side)   ) * dummy(i, j, k)

                    rhs(i,j,k,UzF)  = rhs(i,j,k,UzF)  - ( w(i,j,k)   - params_physics%w_ref(ref_side)   ) * dummy(i, j, k)

                    rhs(i,j,k,EF)   = rhs(i,j,k,EF)   - ( es(i,j,k)  - params_physics%es_ref(ref_side)  ) * dummy(i, j, k)

!                    do n=1, params_physics%species 
!                        rhs(i,j,k,YF+n-1)   = rhs(i,j,k,YF+n-1)   - ( Y(i,j,k,n)  - params_physics%Y_ref(n,ref_side)  ) * dummy(i, j, k)
!                    end do

                end do
            end do
        end do

    end if

    if ( .NOT.(params_physics%periodic_BC(2)) ) then
        ! compute boundary side: x+ or x-, note: split domain at domain center
        if ( x0(2) < (params_physics%L(2)/2.0_rk + 1e-12_rk) ) then
            ref_side = 3
        else
            ref_side = 4
        end if

        do k = g+1, Bs(3)+g
            do j = g+1, Bs(2)+g
                do i = g+1, Bs(1)+g

                    rhs(i,j,k,rhoF) = rhs(i,j,k,rhoF) - ( rho(i,j,k) - params_physics%rho_ref(ref_side) ) * dummy(i, j, k)

                    rhs(i,j,k,UxF)  = rhs(i,j,k,UxF)  - ( u(i,j,k)   - params_physics%u_ref(ref_side)   ) * dummy(i, j, k)

                    rhs(i,j,k,UyF)  = rhs(i,j,k,UyF)  - ( v(i,j,k)   - params_physics%v_ref(ref_side)   ) * dummy(i, j, k)

                    rhs(i,j,k,UzF)  = rhs(i,j,k,UzF)  - ( w(i,j,k)   - params_physics%w_ref(ref_side)   ) * dummy(i, j, k)

                    rhs(i,j,k,EF)   = rhs(i,j,k,EF)   - ( es(i,j,k)  - params_physics%es_ref(ref_side)  ) * dummy(i, j, k)

!                    do n=1, params_physics%species 
!                        rhs(i,j,k,YF+n-1)   = rhs(i,j,k,YF+n-1)   - ( Y(i,j,k,n)  - params_physics%Y_ref(n,ref_side)  ) * dummy(i, j, k)
!                    end do

                end do
            end do
        end do

    end if

    if ( .NOT.(params_physics%periodic_BC(3)) ) then
        ! compute boundary side: x+ or x-, note: split domain at domain center
        if ( x0(3) < (params_physics%L(3)/2.0_rk + 1e-12_rk) ) then
            ref_side = 5
        else
            ref_side = 6
        end if

        do k = g+1, Bs(3)+g
            do j = g+1, Bs(2)+g
                do i = g+1, Bs(1)+g

                    rhs(i,j,k,rhoF) = rhs(i,j,k,rhoF) - ( rho(i,j,k) - params_physics%rho_ref(ref_side) ) * dummy(i, j, k)

                    rhs(i,j,k,UxF)  = rhs(i,j,k,UxF)  - ( u(i,j,k)   - params_physics%u_ref(ref_side)   ) * dummy(i, j, k)

                    rhs(i,j,k,UyF)  = rhs(i,j,k,UyF)  - ( v(i,j,k)   - params_physics%v_ref(ref_side)   ) * dummy(i, j, k)

                    rhs(i,j,k,UzF)  = rhs(i,j,k,UzF)  - ( w(i,j,k)   - params_physics%w_ref(ref_side)   ) * dummy(i, j, k)

                    rhs(i,j,k,EF)   = rhs(i,j,k,EF)   - ( es(i,j,k)  - params_physics%es_ref(ref_side)  ) * dummy(i, j, k)

!                    do n=1, params_physics%species 
!                        rhs(i,j,k,YF+n-1)   = rhs(i,j,k,YF+n-1)   - ( Y(i,j,k,n)  - params_physics%Y_ref(n,ref_side)  ) * dummy(i, j, k)
!                    end do

                end do
            end do
        end do

    end if

end subroutine RHS_3D_CANTERA_navier_stokes_reactive_non_periodicBC
