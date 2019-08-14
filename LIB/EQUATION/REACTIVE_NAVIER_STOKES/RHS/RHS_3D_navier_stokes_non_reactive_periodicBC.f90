!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name RHS_3D_navier_stokes.f90
!> \version 0.5
!> \author msr
!
!> \brief RHS for 3D navier stokes equation
!
!>
!! input:    - datafield, grid parameter, velocity, diffusion coefficient, derivative order \n
!! output:   - RHS(datafield) \n
!!
!!
!! = log ======================================================================================
!! \n
!! 14/02/17 - create
!
! ********************************************************************************************

subroutine RHS_3D_navier_stokes_non_reactive_periodicBC(params_physics, Bs, g, NdF, x0, delta_x, phi, rhs, time)

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
    ! it may happen that some source terms have an explicit time-dependency
    ! therefore the general call has to pass time
    real(kind=rk), intent (in)                              :: time

    ! adiabatic coefficien t
    real(kind=rk)                                           :: gamma_
    ! specific gas constant
    real(kind=rk)                                           :: Rs
    ! isochoric heat capacity
    real(kind=rk)                                           :: Cv
    ! isobaric heat capacity
    real(kind=rk)                                           :: Cp
    ! prandtl number
    real(kind=rk)                                           :: Pr
    ! dynamic viscosity
    real(kind=rk)                                           :: mu0, mu_d, T0, two_three

    ! spacing
    real(kind=rk)                                           :: dx, dy, dz

    ! dissipation switch
    logical                                                 :: dissipation

    ! variables
    real(kind=rk)                                           :: rho(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), u(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), &
                                                               v(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), w(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), &
                                                               p(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), T(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), &
                                                               tau11(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), tau22(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), &
                                                               tau33(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), tau12(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), &
                                                               tau13(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), tau23(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), &
                                                               mu(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), lambda(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

    ! dummy field
    real(kind=rk)                                           :: dummy(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), dummy2(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), &
                                                               dummy3(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), dummy4(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), &
                                                               dummy5(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), dummy6(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

    ! inverse sqrt(rho) field 
    real(kind=rk)                                           :: phi1_inv(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

    ! field indexes
    integer(kind=ik)                                        :: rhoF, UxF, UyF, UzF, EF

    ! loop variables
    integer(kind=ik)                                        :: i, j, k

    ! switch to enabled one sided derivatives
    logical                                                 :: onesided(2,3)

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    ! periodic boundary here, so:
    onesided = .false.

    ! dummy for better performance
    two_three = 2.0_rk / 3.0_rk

    ! field indexes from params
    rhoF = params_physics%rhoF
    UxF  = params_physics%UxF
    UyF  = params_physics%UyF
    UzF  = params_physics%UzF
    EF   = params_physics%EF

    ! set physics parameters for readability
    gamma_      = params_physics%gamma_
    ! only use 1/Rs
    Rs          = 1.0_rk/params_physics%Rs
    Cv          = params_physics%Cv
    Cp          = params_physics%Cp
    Pr          = params_physics%Pr
    mu0         = params_physics%mu0
    dissipation = params_physics%dissipation
    T0          = params_physics%T0

    ! primitive variables
    ! use rhs as dummy fields
    call convert_to_primitive( params_physics, phi, rhs )
 
    rho = rhs(:,:,:,rhoF)
    u   = rhs(:,:,:,UxF)
    v   = rhs(:,:,:,UyF)
    w   = rhs(:,:,:,UzF)
    p   = rhs(:,:,:,EF)

    ! compute 1/rho for better performance
    phi1_inv(:,:,:)  = 1.0_rk / phi(:,:,:,rhoF)

    ! Compute mu and T
    if (dissipation) then
        do k = 1, Bs(3)+2*g
            do j = 1, Bs(2)+2*g
                do i = 1, Bs(1)+2*g
                    T(i,j,k) = p(i,j,k) * phi1_inv(i,j,k) * phi1_inv(i,j,k) * Rs
                end do
            end do
        end do

        ! viscosity
        select case(params_physics%viscosity_model)
            case('constant')
                mu   = mu0
            case('sutherland')

                do k = 1, Bs(3)+2*g
                    do j = 1, Bs(2)+2*g
                        do i = 1, Bs(1)+2*g
                            mu(i,j,k) = mu0 * ( T(i,j,k) / T0)**1.5_rk * (T0 + 110.4_rk) / (T(i,j,k) + 110.4_rk)
                        end do
                    end do
                end do

            case default
                call abort(161018002,"ERROR in RHS: unknown viscosity model")

        end select

        mu_d = 0.0_rk

        ! thermal conductivity
        lambda= Cp * mu/Pr
    
    end if

    ! discretization constant
    dx = delta_x(1)
    dy = delta_x(2)
    dz = delta_x(3)

    ! necessary to use common code with cantera RHS
    rhs(:,:,:,EF) = 0.0_rk

    !#########################################################################################
    ! inlcude RHS for periodic BC
    include "RHS/RHS_3D_common_code.f90"

    !#########################################################################################
    ! RHS of energy equation:  p_t = -gamma*div(U_tilde p) + gamm1 *U x grad(p)
    do k = g-1, Bs(3)+g+2
        do j = g-1, Bs(2)+g+2
            do i = g-1, Bs(1)+g+2
                dummy(i,j,k)  = u(i,j,k)*p(i,j,k)
                dummy2(i,j,k) = v(i,j,k)*p(i,j,k)
                dummy3(i,j,k) = w(i,j,k)*p(i,j,k)
            end do
        end do
    end do
    call diff_wrapper_3D( Bs, g, dummy,  'periodic_u_x', onesided, dx=dx, dudx=dummy4)
    call diff_wrapper_3D( Bs, g, dummy2, 'periodic_u_y', onesided, dy=dy, dudy=dummy5)
    call diff_wrapper_3D( Bs, g, dummy3, 'periodic_u_z', onesided, dz=dz, dudz=dummy6)

    do k = g+1, Bs(3)+g
        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                rhs(i,j,k,EF) = ( gamma_ - 1.0_rk ) * rhs(i,j,k,EF) - gamma_ * ( dummy4(i,j,k) + dummy5(i,j,k) + dummy6(i,j,k) )
            end do
        end do 
    end do 

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
                                                           + w(g+1:Bs(1)+g,  g+1:Bs(2)+g, g+1:Bs(3)+g) * dummy3(g+1:Bs(1)+g,  g+1:Bs(2)+g, g+1:Bs(3)+g) ) * ( gamma_ - 1.0_rk )

    end if

end subroutine RHS_3D_navier_stokes_non_reactive_periodicBC
