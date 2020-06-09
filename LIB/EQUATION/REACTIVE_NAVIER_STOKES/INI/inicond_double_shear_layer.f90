!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name inicond_double_shear_layer.f90
!> \version 0.5
!> \author msr
!
!> \brief initialize 2D double shear layer setup
!>
!! input:    - params, block coordinates \n
!! output:   - state vector \n
!!
!!
!! = log ======================================================================================
!! \n
!! 05/06/20 - create
!
! ********************************************************************************************

subroutine inicond_double_shear_layer( params_physics, phi, x0, dx )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> parameter struct
    type(type_params_rns), intent(inout)    :: params_physics

    !> state vector
    real(kind=rk), intent(inout)            :: phi(:,:,:,:)
    
    !> spacing and origin of block
    real(kind=rk), intent(in)               :: x0(1:3), dx(1:3)

    ! grid parameter
    integer(kind=ik)                        :: Bs(3), g
    ! coordinates, domain size and scale parameter
    real(kind=rk)                           :: x, y, z, L(3), width(2), xPos(2), yPos
    ! loop variables
    integer(kind=ik)                        :: i, j
    ! field indexes
    integer(kind=ik)                        :: rhoF, UxF, UyF, EF
    ! initial values
    real(kind=rk)                           :: p0, rho0, u0

!---------------------------------------------------------------------------------------------
! variables initialization

    ! grid parameter
    Bs  = params_physics%Bs
    g   = params_physics%g

    ! domain size
    L   = params_physics%L

    ! field indexes from params
    rhoF = params_physics%rhoF
    UxF  = params_physics%UxF
    UyF  = params_physics%UyF
    EF   = params_physics%EF

    ! scale parameter, readability
    xPos(1)    = params_physics%inicond_position(1)
    xPos(2)    = params_physics%inicond_position(1)
    yPos       = params_physics%inicond_position(2)
    width(1:2) = params_physics%inicond_scales(1:2)

    ! initial values
    rho0 = params_physics%inicond_rho
    p0   = params_physics%inicond_p
    u0   = 0.0_rk ! placeholder

!---------------------------------------------------------------------------------------------
! main body

    ! place layer
    xPos(1) = xPos(1) - width(1)/2.0_rk
    xPos(2) = xPos(2) + width(1)/2.0_rk
    yPos    = yPos * L(2)

    if (params_physics%d == 2) then
        ! 2D case
        ! create shear layer
        do i = 1, Bs(1)+2*g
            do j = 1, Bs(2)+2*g

                ! compute x,y coordinates from spacing and origin
                x = dble(i-(g+1)) * dx(1) + x0(1)
                y = dble(j-(g+1)) * dx(2) + x0(2)

                ! Uy, rho
                if ( x <= 0.5_rk*L(1) ) then
                    phi(i, j, 1, UyF) = dtanh( width(2)/L(1) * ( x - xPos(1) ) ) + u0
                    phi(i, j, 1, rhoF) = ( rho0 + dtanh( width(2) * ( x - xPos(1)  ) ) ) / 2.0_rk + rho0
                else
                    phi(i, j, 1, UyF) = dtanh( width(2)/L(1) * ( xPos(2) - x ) ) + u0
                    phi(i, j, 1, rhoF) = ( rho0 + dtanh( width(2) * ( xPos(2) - x ) ) ) / 2.0_rk + rho0
                end if


                ! Ux
                phi(i, j, 1, UxF) = 0.01_rk * dsin( 2.0_rk * pi * ( y - yPos  ) )

                ! p
                phi(i, j, 1, EF) = p0

            end do
        end do

    else
        ! 3D case
        call abort(140819001,"ERROR: can not use double shear layer initial condition for a 3D computation.")

    endif

    ! note: set up skew symmetric form here
    ! set density
    phi(:, :, :, rhoF) = dsqrt(phi( :, :, :, rhoF))
    ! skew symmetric velocity
    phi( :, :, :, UxF) = phi( :, :, :, UxF) * phi( :, :, :, rhoF)
    phi( :, :, :, UyF) = phi( :, :, :, UyF) * phi( :, :, :, rhoF)  

end subroutine inicond_double_shear_layer

