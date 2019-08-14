!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name compute_dilatational_dissipation.f90
!> \version 0.5
!> \author msr
!
!> \brief compute turbulent dissipation, dilatational part
!>
!! input:    - velocity components, spacing vector, grid parameter \n
!! output:   - tensor components \n
!!
!!
!! = log ======================================================================================
!! \n
!! 13/08/19 - create
!
! ********************************************************************************************

subroutine compute_dilatational_dissipation(u, v, w, dx, Bs, g, d, dissipation)

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> local datafields, velocity components
    real(kind=rk), dimension(:,:,:), intent(in)    :: u, v, w

    !> spacing of the block
    real(kind=rk), dimension(3), intent(in)        :: dx

    !> grid parameters, dimension
    integer(kind=ik), intent(in)                   :: Bs(3), g, d

    !> dissipation tensor, should be given with 3 or 6 components
    ! order of components
    ! 2D: S11, S12, S22
    ! 3D: S11, S12, S13, S22, S23, S33
    real(kind=rk), dimension(:,:,:,:), intent(out) :: dissipation


    !> derivatives
    real(kind=rk)                                  :: u_dx, v_dy, w_dz
    !> inverse of dx, dy, dz
    real(kind=rk)                                  :: dx_inv, dy_inv, dz_inv
    ! loop variables
    integer(kind=ik)                               :: i, j, k
    ! coefficients for Tam&Webb
    real(kind=rk)                                  :: a(-3:3)
!---------------------------------------------------------------------------------------------
! variables initialization

    ! Tam & Webb, 4th order optimized (for first derivative)
    a = (/-0.02651995_rk, +0.18941314_rk, -0.79926643_rk, 0.0_rk, &
        0.79926643_rk, -0.18941314_rk, 0.02651995_rk/)

    dx_inv = 1.0_rk / dx(1)
    dy_inv = 1.0_rk / dx(2)
    dz_inv = 1.0_rk / dx(3)

!---------------------------------------------------------------------------------------------
! main body

    if (d==2) then
        ! 2D
        call abort(130819001,"ERROR: computation for 2D dissipation not implemented.")

    else 
        ! 3D
        do k = g+1, Bs(3)+g
            do j = g+1, Bs(2)+g
                do i = g+1, Bs(1)+g

                    ! compute derivatives
                    u_dx = (a(-3)*u(i-3, j, k) + a(-2)*u(i-2, j, k) &
                         +  a(-1)*u(i-1, j, k) + a(0) *u(i,   j, k)    &
                         +  a(+1)*u(i+1, j, k) + a(+2)*u(i+2, j, k) &
                         +  a(+3)*u(i+3, j, k))*dx_inv

                    v_dy = (a(-3)*v(i, j-3, k) + a(-2)*v(i, j-2, k) &
                         +  a(-1)*v(i, j-1, k) + a(0) *v(i, j,   k) &
                         +  a(+1)*v(i, j+1, k) + a(+2)*v(i, j+2, k) &
                         +  a(+3)*v(i, j+3, k))*dy_inv

                    w_dz = (a(-3)*w(i, j, k-3) + a(-2)*w(i, j, k-2) &
                         +  a(-1)*w(i, j, k-1) + a(0) *w(i, j, k  ) &
                         +  a(+1)*w(i, j, k+1) + a(+2)*w(i, j, k+2) &
                         +  a(+3)*w(i, j, k+3))*dz_inv


                    ! dissipation tensor
                    ! S11 
                    dissipation(i, j, k, 1) = u_dx**2.0_rk 
                    ! S12 
                    dissipation(i, j, k, 2) = 2.0_rk * u_dx * v_dy
                    ! S13 
                    dissipation(i, j, k, 3) = 2.0_rk * u_dx * w_dz
                    ! S22 
                    dissipation(i, j, k, 4) = v_dy**2.0_rk  
                    ! S23 
                    dissipation(i, j, k, 5) = 2.0_rk * v_dy * w_dz
                    ! S33 = w_dz + w_dz
                    dissipation(i, j, k, 6) = w_dz**2.0_rk 

                end do
            end do
        end do

    end if

end subroutine compute_dilatational_dissipation
