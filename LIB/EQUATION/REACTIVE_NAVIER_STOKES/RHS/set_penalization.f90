!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name set_penalization.f90
!> \version 0.5
!> \author msr
!
!> \brief set penalization sponge 
!!
!! input:    - params struct, grid and domain parameter, periodic/non-periodic BCs \n
!! output:   - sponge field \n
!!
!!
!! = log ======================================================================================
!! \n
!! 04/07/19 - create
!
! ********************************************************************************************

subroutine set_penalization( params_physics, Bs, g, x0, dx, sponge )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> navier stokes params struct, note: contails all parameters needed by RHS
    type(type_params_rns), intent(inout) :: params_physics

    !> grid parameter
    integer(kind=ik), intent(in)         :: Bs(3), g

    !> block origin, spacing
    real(kind=rk), intent(in)            :: x0(3), dx(3)

    !> sponge field
    real(kind=rk), intent(out)           :: sponge(:, :, :)

    !> loop variables
    integer(kind=ik)                     :: i, j, k
    !> dummy variables
    real(kind=rk)                        :: x, y, z, L(3), sponge_width, &
                                            sponge_strength, wall_distance, ref_size
    !> periodic BC switch from wabbit core
    logical                              :: periodic_BC(3)

!---------------------------------------------------------------------------------------------
! variables initialization

    ! readability
    L               = params_physics%L
    sponge_width    = params_physics%sponge_width
    sponge_strength = params_physics%sponge_strength

    ! reference size, choose minimum domain size
    if ( params_physics%d == 2 ) then
        ref_size    = min(L(1), L(2))
    else
        ref_size    = min(L(1), L(2), L(3))
    end if

    ! periodic switch from wabbit core
    periodic_BC     = params_physics%periodic_BC

!---------------------------------------------------------------------------------------------
! main body

    if ( params_physics%d == 2 ) then

        ! 2D case
        do i = g+1, Bs(1)+g
            do j = g+1, Bs(2)+g

                ! min distance to wall
                x = dble(i-(g+1)) * dx(1) + x0(1)
                y = dble(j-(g+1)) * dx(2) + x0(2)

                ! wall distance is only set for non-periodic boundaries
                ! otherwise wall distance is very large, so tanh function return zero
                wall_distance = 9999.99_rk
                if ( .NOT.(periodic_BC(1)) ) wall_distance = min( x, L(1)-x, wall_distance )
                if ( .NOT.(periodic_BC(2)) ) wall_distance = min( y, L(2)-y, wall_distance )

                sponge(i,j,1) = sponge_strength * (-dtanh( 60.0_rk/ref_size * ( wall_distance - sponge_width*ref_size ) ) + 1.0_rk)/2.0_rk

            end do
        end do


    else

        ! 3D case
        do i = g+1, Bs(1)+g
            do j = g+1, Bs(2)+g
                do k = g+1, Bs(3)+g

                    ! min distance to wall
                    x = dble(i-(g+1)) * dx(1) + x0(1)
                    y = dble(j-(g+1)) * dx(2) + x0(2)
                    z = dble(k-(g+1)) * dx(3) + x0(3)

                    ! wall distance is only set for non-periodic boundaries
                    ! otherwise wall distance is very large, so tanh function return zero
                    wall_distance = 9999.99_rk
                    if ( .NOT.(periodic_BC(1)) ) wall_distance = min( x, L(1)-x, wall_distance )
                    if ( .NOT.(periodic_BC(2)) ) wall_distance = min( y, L(2)-y, wall_distance )
                    if ( .NOT.(periodic_BC(3)) ) wall_distance = min( z, L(3)-z, wall_distance )

                    sponge(i,j,k) = sponge_strength * (-dtanh( 60.0_rk/ref_size * ( wall_distance - sponge_width*ref_size ) ) + 1.0_rk)/2.0_rk

                end do
            end do
        end do

    end if

end subroutine set_penalization
