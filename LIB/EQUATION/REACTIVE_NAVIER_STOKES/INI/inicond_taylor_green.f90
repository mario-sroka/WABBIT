!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name inicond_taylor_green.f90
!> \version 0.5
!> \author msr
!
!> \brief initialize taylor green vortices for 3D cases
!>
!! input:    - params, state vector \n
!! output:   - velocity in state vector \n
!!
!!
!! = log ======================================================================================
!! \n
!! 14/08/19 - create
!
! ********************************************************************************************

subroutine inicond_taylor_green( params_physics, phi, x0, dx )

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
    real(kind=rk)                           :: x, y, z, L(3), V0, N
    ! loop variables
    integer(kind=ik)                        :: i, j, k
    ! field indexes
    integer(kind=ik)                        :: rhoF, UxF, UyF, UzF, EF

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
    UzF  = params_physics%UzF
    EF   = params_physics%EF

    ! scale parameter, readability
    V0  = params_physics%inicond_scales(1)
    N   = params_physics%inicond_scales(2)

!---------------------------------------------------------------------------------------------
! main body

    ! check chemistry
    if (params_physics%chemistry_model == 'cantera') then
        call abort(090819002,"ERROR: initial condition taylor green fro cantera chemistry not implemented yet.")
    end if    

    if (params_physics%d == 2) then
        ! 2D case
        call abort(140819001,"ERROR: can not use taylor green initial condition for a 2D computation.")

    else
        ! 3D case
        ! create vortices
        do i = g+1, Bs(1)+g
            do j = g+1, Bs(2)+g
                do k = g+1, Bs(3)+g

                    ! coordinates
                    x = dble(i-(g+1)) * dx(1) + x0(1)
                    y = dble(j-(g+1)) * dx(2) + x0(2)
                    z = dble(k-(g+1)) * dx(3) + x0(3)

                    ! shift coordinates to enforce periodicity
                    x = x - 0.5_rk * L(1)
                    y = y - 0.5_rk * L(2)
                    z = z - 0.5_rk * L(3)

                    ! set velocity
                    phi(i, j, k, UxF) =  V0 * dsin(N*2*pi*x/L(1)) * dcos(N*2*pi*y/L(2)) * dcos(N*2*pi*z/L(3))
                    phi(i, j, k, UyF) = -V0 * dcos(N*2*pi*x/L(1)) * dsin(N*2*pi*y/L(2)) * dcos(N*2*pi*z/L(3))
                    phi(i, j, k, UzF) = 0.0_rk

                    ! set pressure
                    phi(i, j, k, EF)  = params_physics%inicond_p + params_physics%inicond_rho * V0**2.0_rk / 16.0_rk &
                                      * ( dcos(N*4*pi*x/L(1)) + dcos(N*4*pi*y/L(2)) ) * ( dcos(N*4*pi*z/L(3)) + 2.0_rk )

                end do
            end do
        end do

    endif

    ! set density
    phi(:, :, :, rhoF) = dsqrt(params_physics%inicond_rho)
    
    ! skew symmetric velocity
    phi( :, :, :, UxF) = phi( :, :, :, UxF) * phi( :, :, :, rhoF)
    phi( :, :, :, UyF) = phi( :, :, :, UyF) * phi( :, :, :, rhoF)
    if (params_physics%d == 3) phi( :, :, :, UzF) = phi( :, :, :, UzF) * phi( :, :, :, rhoF)    

end subroutine inicond_taylor_green

