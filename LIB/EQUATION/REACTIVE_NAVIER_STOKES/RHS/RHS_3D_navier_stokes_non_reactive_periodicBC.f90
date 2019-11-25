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

subroutine RHS_3D_navier_stokes_non_reactive_periodicBC(params_physics, Ds, Bs, g, NdF, x0, delta_x, phi, rhs, time)

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> navier stokes params struct, note: contails all parameters needed by RHS
    type(type_params_rns), intent(inout)                    :: params_physics
    !> grid parameter
    integer(kind=ik), intent(inout)                         :: g, Bs(3), NdF, Ds(3)
    !> rhs parameter
    real(kind=rk), dimension(3), intent(in)                 :: x0, delta_x
    !> datafields
    real(kind=rk), intent(inout)                            :: phi(Ds(1), Ds(2), Ds(3), NdF)
    ! rhs array
    real(kind=rk),intent(inout)                             :: rhs(Ds(1), Ds(2), Ds(3), NdF)
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
    real(kind=rk), allocatable, save                        :: rho(:,:,:), u(:,:,:), v(:,:,:), w(:,:,:), p(:,:,:), T(:,:,:), &
                                                               tau11(:,:,:), tau22(:,:,:), tau33(:,:,:), tau12(:,:,:), &
                                                               tau13(:,:,:), tau23(:,:,:), mu(:,:,:), lambda(:,:,:)

    ! dummy field
    real(kind=rk), allocatable, save                        :: dummy(:,:,:), dummy2(:,:,:), dummy3(:,:,:), dummy4(:,:,:), &
                                                               dummy5(:,:,:), dummy6(:,:,:)

    ! inverse sqrt(rho) field 
    real(kind=rk), allocatable, save                        :: phi1_inv(:,:,:)

    ! field indexes
    integer(kind=ik)                                        :: rhoF, UxF, UyF, UzF, EF

    ! loop variables
    integer(kind=ik)                                        :: i, j, k

    ! array indexes
    integer(kind=ik)                                        :: gp, gm, Bp(3), Bm(3)

    ! switch to enabled one sided derivatives
    logical                                                 :: onesided(2,3)

!---------------------------------------------------------------------------------------------
! interfaces

! REMOVE
mu_d = 0.0_rk

!---------------------------------------------------------------------------------------------
! variables initialization

    ! compute array indexes
    gp   = g+1
    gm   = g-1
    Bp   = Bs(:) + g
    Bm   = Bs(:) + g + 2

    ! allocate dummy fields
    if ( .not. allocated(dummy) ) allocate( dummy(Ds(1), Ds(2), Ds(3)), dummy2(Ds(1), Ds(2), Ds(3)), &
                                            dummy3(Ds(1), Ds(2), Ds(3)), dummy4(Ds(1), Ds(2), Ds(3)), &
                                            dummy5(Ds(1), Ds(2), Ds(3)), dummy6(Ds(1), Ds(2), Ds(3)), &
                                            phi1_inv(Ds(1), Ds(2), Ds(3)), rho(Ds(1), Ds(2), Ds(3)), &
                                            u(Ds(1), Ds(2),Ds(3)), v(Ds(1), Ds(2),Ds(3)), w(Ds(1), Ds(2), Ds(3)), &
                                            p(Ds(1), Ds(2),Ds(3)), T(Ds(1), Ds(2), Ds(3)), &
                                            tau11(Ds(1),Ds(2), Ds(3)), tau22(Ds(1), Ds(2), Ds(3)), &
                                            tau33(Ds(1),Ds(2), Ds(3)), tau12(Ds(1), Ds(2), Ds(3)), &
                                            tau13(Ds(1),Ds(2), Ds(3)), tau23(Ds(1), Ds(2), Ds(3)), &
                                            mu(Ds(1), Ds(2),Ds(3)), lambda(Ds(1), Ds(2), Ds(3)) )

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

    ! compute 1/rho for better performance
    phi1_inv(:,:,:)  = 1.0_rk / phi(:,:,:,rhoF)

    ! primitive variables
    ! use rhs as dummy fields
    call convert_to_primitive( params_physics, phi, rhs, phi1_inv )
 
    rho = rhs(:,:,:,rhoF)
    u   = rhs(:,:,:,UxF)
    v   = rhs(:,:,:,UyF)
    w   = rhs(:,:,:,UzF)
    p   = rhs(:,:,:,EF)  

    ! Compute mu and T
    if (dissipation) then

        T(:, :, :) = p(:, :, :) * phi1_inv(:, :, :) * phi1_inv(:, :, :) * Rs

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

        ! mu_d - 2/3 mu, mu_d = 0
        dummy6(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = - two_three * mu(gm:Bm(1), gm:Bm(2), gm:Bm(3))

        ! thermal conductivity
        lambda(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = Cp * mu(gm:Bm(1), gm:Bm(2), gm:Bm(3)) / Pr
    
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
    dummy( gm:Bm(1), gm:Bm(2), gm:Bm(3)) = u(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * p(gm:Bm(1), gm:Bm(2), gm:Bm(3))
    dummy2(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = v(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * p(gm:Bm(1), gm:Bm(2), gm:Bm(3))
    dummy3(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = w(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * p(gm:Bm(1), gm:Bm(2), gm:Bm(3))

    call diff_wrapper_3D( Bs, g, dummy,  'diff_x_c_gp', onesided, dx=dx, dudx=dummy4)
    call diff_wrapper_3D( Bs, g, dummy2, 'diff_y_c_gp', onesided, dy=dy, dudy=dummy5)
    call diff_wrapper_3D( Bs, g, dummy3, 'diff_z_c_gp', onesided, dz=dz, dudz=dummy6)

    dummy4(gp:Bp(1), gp:Bp(2), gp:Bp(3)) = dummy4( gp:Bp(1), gp:Bp(2), gp:Bp(3)) &
                                         + dummy5( gp:Bp(1), gp:Bp(2), gp:Bp(3)) &
                                         + dummy6( gp:Bp(1), gp:Bp(2), gp:Bp(3))

    rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), EF) = rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), EF) * ( gamma_ - 1.0_rk ) &
                                          - dummy4(gp:Bp(1), gp:Bp(2), gp:Bp(3)) * gamma_

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
