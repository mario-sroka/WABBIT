!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name inicond_vortex.f90
!> \version 0.5
!> \author msr
!
!> \brief initialize velocity field with vortex structure, note: works only on velocity field,
!> returns primitive variables
!>
!! input:    - params, grid and spacing parameter\n
!! output:   - block data \n
!!
!!
!! = log ======================================================================================
!! \n
!! 11/10/18 - create
!
! ********************************************************************************************

subroutine inicond_vortex( params_physics, phi, x0, dx, d )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> parameter struct
    type(type_params_rns), intent(inout)    :: params_physics

    !> actual block data (note this routine acts only on one block)
    real(kind=rk), intent(inout)            :: phi(:,:,:,:)

    !> spacing and origin of block
    real(kind=rk), intent(in)               :: x0(1:3), dx(1:3)

    !> dimension
    integer(kind=ik), intent(in)            :: d

    ! grid parameter
    integer(kind=ik)                        :: Bs(d), g
    ! coordinates and domain parameter
    real(kind=rk)                           :: x, y, z, L(d), x2, y2
    ! gauss blob center coordinates, width, size, radius
    real(kind=rk)                           :: mux(d), sigma(2), omega, r, r2
    ! loop variables
    integer(kind=ik)                        :: i, j
    ! field indexes
    integer(kind=ik)                        :: UxF, UyF
    ! dummy fields
    real(kind=rk), allocatable              :: dummy(:, :), dummy2(:, :)

!---------------------------------------------------------------------------------------------
! variables initialization

    Bs      = params_physics%Bs(1:d)
    g       = params_physics%g
    L       = params_physics%L(1:d)

    ! use inicond scales from ini file
    sigma(1)= params_physics%inicond_scales(1)
    sigma(2)= params_physics%inicond_scales(2)
    omega   = params_physics%inicond_scales(3)

    ! center point of the pressure blob
    mux = params_physics%inicond_position(1:d)

    ! field indexes from params
    UxF  = params_physics%UxF
    UyF  = params_physics%UyF

    ! dummy fields
    allocate( dummy(1:Bs(1)+2*g, 1:Bs(2)+2*g), dummy2(1:Bs(1)+2*g, 1:Bs(2)+2*g) )

!---------------------------------------------------------------------------------------------
! main body

    ! check dimension, works only in 2D!
    if ( d /= 2 ) call abort(040919002,"ERROR: inicond_vortex works only in 2 dimensions.")

    do i = g+1, Bs(1)+g
        do j = g+1, Bs(2)+g

            ! compute x,y coordinates from spacing and origin
            x = dble(i-(g+1)) * dx(1) + x0(1)
            y = dble(j-(g+1)) * dx(2) + x0(2)

            ! coordinate shift
            x2 = x - mux(1)
            y2 = y - mux(2) - 6.0_rk*sigma(1)*L(2)

            x  = x - mux(1)
            y  = y - mux(2) + 6.0_rk*sigma(2)*L(2)

            ! polar coordinates
            r  = dsqrt( x**2.0_rk + y**2.0_rk )
            r2 = dsqrt( x2**2.0_rk + y2**2.0_rk )

            if ( r  < 1.0e-12_rk ) r  = 1.0e-12_rk
            if ( r2 < 1.0e-12_rk ) r2 = 1.0e-12_rk

            ! set dummies
            dummy(i,j)  = dexp( - r**2.0_rk  / (sigma(1))**2.0_rk )
            dummy2(i,j) = dexp( - r2**2.0_rk / (sigma(1))**2.0_rk )

            ! set velocity
            phi(i,j,1,UxF) = dummy(i,j) * ( omega/(2.0_rk*pi) * y/r   * &
                           ( 1 - dexp( - x**2.0_rk / sigma(1)**2.0_rk - y**2.0_rk / sigma(2)**2.0_rk ) ) ) &
                           - dummy2(i,j) * ( omega/(2.0_rk*pi) * y2/r2 * &
                           ( 1 - dexp( - x2**2.0_rk / sigma(1)**2.0_rk - y2**2.0_rk / sigma(2)**2.0_rk ) ) )

            phi(i,j,1,UyF) = - dummy(i,j) * ( omega/(2.0_rk*pi) * x/r   * &
                           ( 1 - dexp( - x**2.0_rk / sigma(1)**2.0_rk - y**2.0_rk / sigma(2)**2.0_rk ) ) ) &
                           + dummy2(i,j) * ( omega/(2.0_rk*pi) * x2/r2 * &
                           ( 1 - dexp( - x2**2.0_rk / sigma(1)**2.0_rk - y2**2.0_rk / sigma(2)**2.0_rk ) ) )

        end do
    end do

    ! clean up
    deallocate(dummy, dummy2)

end subroutine inicond_vortex


