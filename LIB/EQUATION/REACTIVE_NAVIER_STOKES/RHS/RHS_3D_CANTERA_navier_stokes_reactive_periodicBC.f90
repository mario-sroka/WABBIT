!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name RHS_3D_navier_stokes.f90
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

    ! primitive variables
    ! use rhs as dummy fields
    call convert_to_primitive( params_physics, phi, rhs )
 
    rho = rhs(:,:,:,rhoF)
    u   = rhs(:,:,:,UxF)
    v   = rhs(:,:,:,UyF)
    w   = rhs(:,:,:,UzF)
    es  = rhs(:,:,:,EF)
    do n = 1, params_physics%species
        Y(:,:,:,n) = rhs(:,:,:,YF+n-1)
    end do

    ! compute 1/rho for better performance
    phi1_inv(:,:,:)  = 1.0_rk / phi(:,:,:,rhoF)

    ! Compute Rs, W, X, D
    if (dissipation) then

        do k = 1, Bs(3)+2*g
            do j = 1, Bs(2)+2*g
                do i = 1, Bs(1)+2*g

                    ! mean molecular weight
                    Wk = 0.0_rk
                    do n = 1, params_physics%species
                        Wk = Wk + Wk_species_inv(n) * Y(i,j,k,n)
                    end do
                    ! note: W is actually 1/W !!
                    Wk_inv = 1.0_rk / Wk
                    ! gas constant
                    Rs = Wk * params_physics%gas_constant

                    ! mole fraction
                    do n = 1, params_physics%species
                        X(i,j,k,n)         = Wk_inv * Wk_species_inv(n) * Y(i,j,k,n)
                    end do

                    ! diffusion coefficient
                    do n = 1, params_physics%species
                        D(i,j,k,n) = params_physics%Wk(n) * Wk
                    end do

                end do
            end do
        end do

        mu_d = 0.0_rk

    end if

    ! discretization constant
    dx = delta_x(1)
    dy = delta_x(2)
    dz = delta_x(3)

    ! chemistry
    dh  = params_physics%dh

    !#########################################################################################
    ! CANTERA
    do k = 1, Bs(3)+2*g
        do j = 1, Bs(2)+g+g
            do i = 1, Bs(1)+g+g

                if ( (k<g-1) .or. (j<g-1) .or. (i<g-1) .or. (k>Bs(3)+g+2) .or. (j>Bs(2)+g+2) .or. (i>Bs(1)+g+2) ) then

                    ! setState call depends on choosen mechanism
                    ! OHN
                    call setMoleFractions(gas, X(i,j,k,:))
                    call setState_UV(gas, es(i,j,k) + sum( dh(:) * Y(i,j,k,:) ), phi1_inv(i,j,k)*phi1_inv(i,j,k)  )

                    ! temperature and pressure
                    T(i,j,k) = temperature(gas)
                    p(i,j,k) = pressure(gas)

                elseif ( (k<g+1) .or. (j<g+1) .or. (i<g+1) .or. (k>Bs(3)+g) .or. (j>Bs(2)+g) .or. (i>Bs(1)+g) ) then

                    ! setState call depends on choosen mechanism
                    ! OHN
                    call setMoleFractions(gas, X(i,j,k,:))
                    call setState_UV(gas, es(i,j,k) + sum( dh(:) * Y(i,j,k,:) ), phi1_inv(i,j,k)*phi1_inv(i,j,k)  )

                    ! temperature and pressure
                    T(i,j,k) = temperature(gas)
                    p(i,j,k) = pressure(gas)

                    ! viscosity
                    mu(i,j,k) = viscosity(gas)

                    ! thermal conductivity
                    lambda(i,j,k) = thermalConductivity(gas)    
                
                    ! diffusion
                    call getMixDiffCoeffs(gas, q)

                    do n = 1, params_physics%species
                        D(i,j,k,n) = D(i,j,k,n)  * q(n) 
                    end do

                else

                    ! setState call depends on choosen mechanism
                    ! OHN
                    call setMoleFractions(gas, X(i,j,k,:))
                    call setState_UV(gas, es(i,j,k) + sum( dh(:) * Y(i,j,k,:) ), phi1_inv(i,j,k)*phi1_inv(i,j,k)  )

                    ! temperature and pressure
                    T(i,j,k) = temperature(gas)
                    p(i,j,k) = pressure(gas)

                    ! viscosity
                    mu(i,j,k) = viscosity(gas)

                    ! thermal conductivity
                    lambda(i,j,k) = thermalConductivity(gas)

                    ! diffusion
                    call getMixDiffCoeffs(gas, q)

                    do n = 1, params_physics%species
                        D(i,j,k,n) = D(i,j,k,n)  * q(n) 
                    end do

                    ! net production srates
                    call getNetProductionRates(gas, q)

                    ! correct reaction rate
                    do n = 1, params_physics%species
                        rhs(i,j,k,YF+n-1) = q(n) * params_physics%Wk(n)
                    end do

                    rhs(i,j,k,EF) = - sum( dh(:) * rhs(i,j,k,YF:YF+params_physics%species-1) )

                end if
    
            end do
        end do
    end do
    !#########################################################################################

!---------------------------------------------------------------------------------------------
! main body

    !#########################################################################################
    ! inlcude RHS for periodic BC
    include "RHS/RHS_3D_common_code.f90"
    !#########################################################################################

    ! RHS of energy equation:  p_t = -gamma*div(U_tilde p) + gamm1 *U x grad(p)
    do k = g-1, Bs(3)+g+2
        do j = g-1, Bs(2)+g+2
            do i = g-1, Bs(1)+g+2
                dummy(i,j,k)  = u(i,j,k) * phi(i,j,k,EF) + u(i,j,k) * p(i,j,k)
                dummy2(i,j,k) = v(i,j,k) * phi(i,j,k,EF) + v(i,j,k) * p(i,j,k)
                dummy3(i,j,k) = w(i,j,k) * phi(i,j,k,EF) + w(i,j,k) * p(i,j,k)
            end do
        end do
    end do

    call diff_wrapper_3D( Bs, g, dummy,  'periodic_u_x', onesided, dx=dx, dudx=dummy4)
    call diff_wrapper_3D( Bs, g, dummy2, 'periodic_u_y', onesided, dy=dy, dudy=dummy5)
    call diff_wrapper_3D( Bs, g, dummy3, 'periodic_u_z', onesided, dz=dz, dudz=dummy6)

    do k = g+1, Bs(3)+g
        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                rhs(i,j,k,EF) = rhs(i,j,k,EF) - ( dummy4(i,j,k) + dummy5(i,j,k) + dummy6(i,j,k) )
            end do
        end do 
    end do 

    !##############################################################################################

    ! diffusion
    ! (rho * D_k * W_k / W * (X_k)_xi)_xi
    !---------------------------------------------------------------------------------------------
    do n = 1, params_physics%species-1

        call diff_wrapper_3D( Bs, g, X(:,:,:,n), 'periodic_u_xyz', onesided, dx=dx, dy=dy, dz=dz, dudx=dummy, dudy=dummy2, dudz=dummy3)

        do k = g-1, Bs(3)+g+2
            do j = g-1, Bs(2)+g+2
                do i = g-1, Bs(1)+g+2

                    ! diffusion
                    dummy(i,j,k)  = rho(i,j,k) * D(i,j,k,n) * dummy(i,j,k)
                    dummy2(i,j,k) = rho(i,j,k) * D(i,j,k,n) * dummy2(i,j,k)
                    dummy3(i,j,k) = rho(i,j,k) * D(i,j,k,n) * dummy3(i,j,k)

                end do
            end do
        end do

        call diff_wrapper_3D( Bs, g, dummy,  'periodic_u_x', onesided, dx=dx, dudx=dummy4)
        call diff_wrapper_3D( Bs, g, dummy2, 'periodic_u_y', onesided, dy=dy, dudy=dummy5)
        call diff_wrapper_3D( Bs, g, dummy3, 'periodic_u_z', onesided, dz=dz, dudz=dummy6)

        do k = g+1, Bs(3)+g
            do j = g+1, Bs(2)+g
                do i = g+1, Bs(1)+g
                    rhs(i,j,k,YF+n-1) = rhs(i,j,k,YF+n-1) + ( dummy4(i,j,k) + dummy5(i,j,k) + dummy6(i,j,k) )
                end do
            end do
        end do

    end do

    ! transport
    ! - (rho * (u_i + Vc_i) * Y_k)_xi
    ! note: YF contains rho*Y_k
    ! -------------------------
    do n = 1, params_physics%species-1

        do k = g-1, Bs(3)+g+2
            do j = g-1, Bs(2)+g+2
                do i = g-1, Bs(1)+g+2
                    dummy(i,j,k)  = ( u(i,j,k) ) * phi(i,j,k,YF+n-1)
                    dummy2(i,j,k) = ( v(i,j,k) ) * phi(i,j,k,YF+n-1)
                    dummy3(i,j,k) = ( w(i,j,k) ) * phi(i,j,k,YF+n-1)
                end do
            end do
        end do

        call diff_wrapper_3D( Bs, g, dummy,  'periodic_u_x', onesided, dx=dx, dudx=dummy4)
        call diff_wrapper_3D( Bs, g, dummy2, 'periodic_u_y', onesided, dy=dy, dudy=dummy5)
        call diff_wrapper_3D( Bs, g, dummy3, 'periodic_u_z', onesided, dz=dz, dudz=dummy6)

        do k = g+1, Bs(3)+g
            do j = g+1, Bs(2)+g
                do i = g+1, Bs(1)+g
                    rhs(i,j,k,YF+n-1) = rhs(i,j,k,YF+n-1) - ( dummy4(i,j,k) + dummy5(i,j,k) + dummy6(i,j,k) )
                end do
            end do
        end do

    end do

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
