!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name RHS_2D_navier_stokes_non_reactive_periodicBC.f90
!> \version 0.5
!> \author msr
!
!> \brief RHS for 2D navier stokes equation
!
!>
!! input:    - datafield, grid parameter, velocity, diffusion coefficient, derivative order \n
!! output:   - RHS(datafield) \n
!!
!!
!! = log ======================================================================================
!! \n
!! 09/06/20 - create
!
! ********************************************************************************************

subroutine RHS_2D_navier_stokes_non_reactive_periodicBC(params_physics, Ds, Bs, g, NdF, x0, delta_x, phi, rhs, time)

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
    real(kind=rk), intent(inout)                            :: phi(Ds(1), Ds(2), 1, NdF)
    ! rhs array
    real(kind=rk),intent(inout)                             :: rhs(Ds(1), Ds(2), 1, NdF)
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
    real(kind=rk)                                           :: mu0, T0, two_three, mu_d

    ! spacing
    real(kind=rk)                                           :: dx, dy, dz

    ! dissipation switch
    logical                                                 :: dissipation

    ! variables
    real(kind=rk), allocatable, save                        :: rho(:,:), u(:,:), v(:,:), p(:,:), T(:,:), &
                                                               tau11(:,:), tau22(:,:), tau12(:,:), tau33(:,:),&
                                                               mu(:,:), lambda(:,:)

    ! dummy field
    real(kind=rk), allocatable, save                        :: dummy(:,:), dummy2(:,:), dummy3(:,:), dummy4(:,:), &
                                                               dummy5(:,:), dummy6(:,:)

    ! inverse sqrt(rho) field 
    real(kind=rk), allocatable, save                        :: phi1_inv(:,:), phi1_DUMMY(:,:,:)

    ! field indexes
    integer(kind=ik)                                        :: rhoF, UxF, UyF, EF

    ! loop variables
    integer(kind=ik)                                        :: i, j, k

    ! array indexes
    integer(kind=ik)                                        :: gp, gm, Bp(3), Bm(3)

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

    ! allocate dummy fields
    if ( .not. allocated(dummy) ) allocate( dummy(Ds(1), Ds(2)), dummy2(Ds(1), Ds(2)), &
                                            dummy3(Ds(1), Ds(2)), dummy4(Ds(1), Ds(2)), &
                                            dummy5(Ds(1), Ds(2)), dummy6(Ds(1), Ds(2)), &
                                            phi1_inv(Ds(1), Ds(2)), phi1_DUMMY(Ds(1), Ds(2), Ds(3)), &
                                            rho(Ds(1), Ds(2)), &
                                            u(Ds(1), Ds(2)), v(Ds(1), Ds(2)), &
                                            p(Ds(1), Ds(2)), T(Ds(1), Ds(2)), &
                                            tau11(Ds(1),Ds(2)), tau22(Ds(1), Ds(2)), &
                                            tau33(Ds(1),Ds(2)), tau12(Ds(1), Ds(2)), &
                                            mu(Ds(1), Ds(2)), lambda(Ds(1), Ds(2)) )

    ! periodic boundary here, so:
    onesided = .false.

    ! dummy for better performance
    two_three = 2.0_rk / 3.0_rk

    ! field indexes from params
    rhoF = params_physics%rhoF
    UxF  = params_physics%UxF
    UyF  = params_physics%UyF
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
    phi1_inv(:,:)  = 1.0_rk / phi(:,:,1,rhoF)

    ! primitive variables
    ! use rhs as dummy fields
    ! use a special DUMMY
    phi1_DUMMY(:,:,1) = phi1_inv
    call convert_to_primitive( params_physics, phi, rhs, phi1_DUMMY(:,:,1:1) )
 
    rho = rhs(:,:,1,rhoF)
    u   = rhs(:,:,1,UxF)
    v   = rhs(:,:,1,UyF)
    p   = rhs(:,:,1,EF)  

    ! Compute mu and T
    if (dissipation) then

        T(:, :) = p(:, :) * phi1_inv(:, :) * phi1_inv(:, :) * Rs

        ! viscosity
        select case(params_physics%viscosity_model)
            case('constant')
                mu   = mu0
            case('sutherland')

                do j = 1, Bs(2)+2*g
                    do i = 1, Bs(1)+2*g
                        mu(i,j) = mu0 * ( T(i,j) / T0)**1.5_rk * (T0 + 110.4_rk) / (T(i,j) + 110.4_rk)
                    end do
                end do

            case default
                call abort(161018002,"ERROR in RHS: unknown viscosity model")

        end select

        ! mu_d - 2/3 mu, mu_d = 0
        mu_d = 0.0_rk

        ! thermal conductivity
        lambda(gm:Bm(1), gm:Bm(2)) = Cp * mu(gm:Bm(1), gm:Bm(2)) / Pr
    
    end if

    ! discretization constant
    dx = delta_x(1)
    dy = delta_x(2)

    ! necessary to use common code with cantera RHS
    rhs(:,:,1,EF) = 0.0_rk

    !#########################################################################################
    ! inlcude RHS for periodic BC
    include "RHS/RHS_2D_common_code.f90"

    !#########################################################################################
    ! RHS of energy equation:  p_t = -gamma*div(U_tilde p) + gamm1 *U x grad(p)
    dummy( gm:Bm(1), gm:Bm(2)) = u(gm:Bm(1), gm:Bm(2)) * p(gm:Bm(1), gm:Bm(2))
    dummy2(gm:Bm(1), gm:Bm(2)) = v(gm:Bm(1), gm:Bm(2)) * p(gm:Bm(1), gm:Bm(2))

    call diff_wrapper_2D( Bs, g, dummy, 'periodic_u_x', onesided, dx=dx, dudx=dummy4)
    call diff_wrapper_2D( Bs, g, dummy2, 'periodic_u_y', onesided, dy=dy, dudy=dummy5)

    dummy4(gp:Bp(1), gp:Bp(2)) = dummy4( gp:Bp(1), gp:Bp(2)) + dummy5( gp:Bp(1), gp:Bp(2))

    rhs(gp:Bp(1), gp:Bp(2), 1, EF) = rhs(gp:Bp(1), gp:Bp(2), 1, EF) * ( gamma_ - 1.0_rk ) &
                                   - dummy4(gp:Bp(1), gp:Bp(2)) * gamma_



end subroutine RHS_2D_navier_stokes_non_reactive_periodicBC
