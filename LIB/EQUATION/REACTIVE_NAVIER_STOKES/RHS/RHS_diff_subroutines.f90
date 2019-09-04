!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name RHS_diff_subroutines.f90
!> \version 0.5
!> \author msr
!
!> \brief file contains all derivations
!
! ********************************************************************************************
!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------

subroutine  diffxy_c_opt( Bs, g, dx, dy, u, dudx, dudy)

    integer(kind=ik), intent(in)    :: g, Bs(3)
    real(kind=rk), intent(in)       :: dx, dy
    real(kind=rk), intent(in)       :: u(Bs(1)+2*g, Bs(2)+2*g)
    real(kind=rk), intent(out)      :: dudx(Bs(1)+2*g, Bs(2)+2*g), dudy(Bs(1)+2*g, Bs(2)+2*g)

    integer                         :: i, j
    real(kind=rk)                   :: dx_inv, dy_inv

    ! - do not use all ghost nodes (note: use only two ghost nodes to get correct second derivatives)
    ! - no one sided stencils necessary
    ! - write loops explicitly,
    ! - use multiplication for dx
    ! - access array in column-major order

    dx_inv = 1.0_rk/(12.0_rk*dx)
    dy_inv = 1.0_rk/(12.0_rk*dy)

    do j = g-1, Bs(2)+g+2
        do i = g-1, Bs(1)+g+2
            dudx(i,j) = ( u(i-2,j) - 8.0_rk*u(i-1,j) + 8.0_rk*u(i+1,j) - u(i+2,j) ) * dx_inv
            dudy(i,j) = ( u(i,j-2) - 8.0_rk*u(i,j-1) + 8.0_rk*u(i,j+1) - u(i,j+2) ) * dy_inv
        end do
    end do

end subroutine diffxy_c_opt

subroutine  diffx_c_opt( Bs, g, dx, u, dudx)

    integer(kind=ik), intent(in)    :: g, Bs(3)
    real(kind=rk), intent(in)       :: dx
    real(kind=rk), intent(in)       :: u(Bs(1)+2*g, Bs(2)+2*g)
    real(kind=rk), intent(out)      :: dudx(Bs(1)+2*g, Bs(2)+2*g)

    integer                         :: i, j
    real(kind=rk)                   :: dx_inv

    ! - do not use ghost nodes
    ! - no one sided stencils necessary
    ! - write loops explicitly,
    ! - use multiplication for dx
    ! - access array in column-major order

    dx_inv = 1.0_rk/(12.0_rk*dx)

    do j = g+1, Bs(2)+g
        do i = g+1, Bs(1)+g
            dudx(i,j) = ( u(i-2,j) - 8.0_rk*u(i-1,j) + 8.0_rk*u(i+1,j) - u(i+2,j) ) * dx_inv
        end do
    end do

end subroutine diffx_c_opt

subroutine  diffy_c_opt( Bs, g, dy, u, dudy)

    integer(kind=ik), intent(in)    :: g, Bs(3)
    real(kind=rk), intent(in)       :: dy
    real(kind=rk), intent(in)       :: u(Bs(1)+2*g, Bs(2)+2*g)
    real(kind=rk), intent(out)      :: dudy(Bs(1)+2*g, Bs(2)+2*g)

    integer                         :: i, j
    real(kind=rk)                   :: dy_inv

    ! - do not use ghost nodes
    ! - no one sided stencils necessary
    ! - write loops explicitly,
    ! - use multiplication for dx
    ! - access array in column-major order

    dy_inv = 1.0_rk/(12.0_rk*dy)

    do j = g+1, Bs(2)+g
        do i = g+1, Bs(1)+g
            dudy(i,j) = ( u(i,j-2) - 8.0_rk*u(i,j-1) + 8.0_rk*u(i,j+1) - u(i,j+2) ) * dy_inv
        end do
    end do

end subroutine diffy_c_opt


! --------
!    3D 
! --------

subroutine  diffx_oc_3D_opt( Bs, g, dx, u, dudx)

    integer(kind=ik), intent(in)    :: g, Bs(3)
    real(kind=rk), intent(in)       :: dx
    real(kind=rk), intent(in)       :: u(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)
    real(kind=rk), intent(out)      :: dudx(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

    integer                         :: i, j, k
    real(kind=rk)                   :: dx_inv

    ! - do not use ghost nodes
    ! - write loops explicitly,
    ! - use multiplication for dx
    ! - access array in column-major order

    dx_inv = 1.0_rk/(12.0_rk*dx)

    do k = g+1, Bs(3)+g
        do j = g+1, Bs(2)+g
            dudx(g+1,j,k) = ( - 18.0_rk*u(g+1,j,k) + 24.0_rk*u(g+2,j,k) - 6.0_rk*u(g+3,j,k) ) * dx_inv
            dudx(g+2,j,k) = ( - 6.0_rk*u(g+1,j,k)                       + 6.0_rk*u(g+3,j,k) ) * dx_inv
            do i = g+3, Bs(1)+g
                dudx(i,j,k) = ( u(i-2,j,k) - 8.0_rk*u(i-1,j,k) + 8.0_rk*u(i+1,j,k) - u(i+2,j,k) ) * dx_inv
            end do

        end do
    end do

end subroutine diffx_oc_3D_opt

subroutine  diffx_co_3D_opt( Bs, g, dx, u, dudx)

    integer(kind=ik), intent(in)    :: g, Bs(3)
    real(kind=rk), intent(in)       :: dx
    real(kind=rk), intent(in)       :: u(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)
    real(kind=rk), intent(out)      :: dudx(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

    integer                         :: i, j, k
    real(kind=rk)                   :: dx_inv

    ! - do not use ghost nodes
    ! - write loops explicitly,
    ! - use multiplication for dx
    ! - access array in column-major order

    dx_inv = 1.0_rk/(12.0_rk*dx)

    do k = g+1, Bs(3)+g
        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g-2
                dudx(i,j,k) = ( u(i-2,j,k) - 8.0_rk*u(i-1,j,k) + 8.0_rk*u(i+1,j,k) - u(i+2,j,k) ) * dx_inv
            end do
            dudx(Bs(1)+g-1,j,k) = ( - 6.0_rk*u(Bs(1)+g-2,j,k)                            + 6.0_rk*u(Bs(1)+g,j,k) ) * dx_inv
            dudx(Bs(1)+g,j,k)   = ( + 6.0_rk*u(Bs(1)+g-2,j,k) - 24.0_rk*u(Bs(1)+g-1,j,k) + 18.0_rk*u(Bs(1)+g,j,k) ) * dx_inv
        end do
    end do

end subroutine diffx_co_3D_opt


subroutine  diffx_oc_3D_opt_2( Bs, g, dx, u, dudx)

    integer(kind=ik), intent(in)    :: g, Bs(3)
    real(kind=rk), intent(in)       :: dx
    real(kind=rk), intent(in)       :: u(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)
    real(kind=rk), intent(out)      :: dudx(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

    integer                         :: i, j, k
    real(kind=rk)                   :: dx_inv

    ! - do not use ghost nodes
    ! - write loops explicitly,
    ! - use multiplication for dx
    ! - access array in column-major order

    dx_inv = 1.0_rk/(12.0_rk*dx)

    do k = g-1, Bs(3)+g+2
        do j = g-1, Bs(2)+g+2
            dudx(g+1,j,k) = ( - 18.0_rk*u(g+1,j,k) + 24.0_rk*u(g+2,j,k) - 6.0_rk*u(g+3,j,k) ) * dx_inv
            dudx(g+2,j,k) = ( - 6.0_rk*u(g+1,j,k)                       + 6.0_rk*u(g+3,j,k) ) * dx_inv
            do i = g+3, Bs(1)+g+2
                dudx(i,j,k) = ( u(i-2,j,k) - 8.0_rk*u(i-1,j,k) + 8.0_rk*u(i+1,j,k) - u(i+2,j,k) ) * dx_inv
            end do

        end do
    end do

end subroutine diffx_oc_3D_opt_2

subroutine  diffx_co_3D_opt_2( Bs, g, dx, u, dudx)

    integer(kind=ik), intent(in)    :: g, Bs(3)
    real(kind=rk), intent(in)       :: dx
    real(kind=rk), intent(in)       :: u(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)
    real(kind=rk), intent(out)      :: dudx(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

    integer                         :: i, j, k
    real(kind=rk)                   :: dx_inv

    ! - do not use ghost nodes
    ! - write loops explicitly,
    ! - use multiplication for dx
    ! - access array in column-major order

    dx_inv = 1.0_rk/(12.0_rk*dx)

    do k = g-1, Bs(3)+g+2
        do j = g-1, Bs(2)+g+2
            do i = g-1, Bs(1)+g-2
                dudx(i,j,k) = ( u(i-2,j,k) - 8.0_rk*u(i-1,j,k) + 8.0_rk*u(i+1,j,k) - u(i+2,j,k) ) * dx_inv
            end do
            dudx(Bs(1)+g-1,j,k) = ( - 6.0_rk*u(Bs(1)+g-2,j,k)                            + 6.0_rk*u(Bs(1)+g,j,k) ) * dx_inv
            dudx(Bs(1)+g,j,k)   = ( + 6.0_rk*u(Bs(1)+g-2,j,k) - 24.0_rk*u(Bs(1)+g-1,j,k) + 18.0_rk*u(Bs(1)+g,j,k) ) * dx_inv
        end do
    end do

end subroutine diffx_co_3D_opt_2

subroutine  diffy_oc_3D_opt( Bs, g, dy, u, dudy)

    integer(kind=ik), intent(in)    :: g, Bs(3)
    real(kind=rk), intent(in)       :: dy
    real(kind=rk), intent(in)       :: u(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)
    real(kind=rk), intent(out)      :: dudy(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

    integer                         :: i, j, k
    real(kind=rk)                   :: dy_inv

    ! - do not use ghost nodes
    ! - write loops explicitly,
    ! - use multiplication for dx
    ! - access array in column-major order

    dy_inv = 1.0_rk/(12.0_rk*dy)

    do k = g+1, Bs(3)+g
        do i = g+1, Bs(1)+g
            dudy(i,g+1,k) = ( - 18.0_rk*u(i,g+1,k) + 24.0_rk*u(i,g+2,k) - 6.0_rk*u(i,g+3,k) ) * dy_inv
            dudy(i,g+2,k) = ( - 6.0_rk*u(i,g+1,k)                       + 6.0_rk*u(i,g+3,k) ) * dy_inv
        end do
        do j = g+3, Bs(2)+g
            do i = g+1, Bs(1)+g
                dudy(i,j,k) = ( u(i,j-2,k) - 8.0_rk*u(i,j-1,k) + 8.0_rk*u(i,j+1,k) - u(i,j+2,k) ) * dy_inv
            end do
        end do
    end do

end subroutine diffy_oc_3D_opt

subroutine  diffy_co_3D_opt( Bs, g, dy, u, dudy)

    integer(kind=ik), intent(in)    :: g, Bs(3)
    real(kind=rk), intent(in)       :: dy
    real(kind=rk), intent(in)       :: u(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)
    real(kind=rk), intent(out)      :: dudy(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

    integer                         :: i, j, k
    real(kind=rk)                   :: dy_inv

    ! - do not use ghost nodes
    ! - write loops explicitly,
    ! - use multiplication for dx
    ! - access array in column-major order

    dy_inv = 1.0_rk/(12.0_rk*dy)

    do k = g+1, Bs(3)+g
        do j = g+1, Bs(2)+g-2
            do i = g+1, Bs(1)+g
                dudy(i,j,k) = ( u(i,j-2,k) - 8.0_rk*u(i,j-1,k) + 8.0_rk*u(i,j+1,k) - u(i,j+2,k) ) * dy_inv
            end do
        end do
        do i = g+1, Bs(1)+g
            dudy(i,Bs(2)+g-1,k) = ( - 6.0_rk*u(i,Bs(2)+g-2,k)                            + 6.0_rk*u(i,Bs(2)+g,k) ) * dy_inv
            dudy(i,Bs(2)+g,k)   = ( + 6.0_rk*u(i,Bs(2)+g-2,k) - 24.0_rk*u(i,Bs(2)+g-1,k) + 18.0_rk*u(i,Bs(2)+g,k) ) * dy_inv
        end do
    end do

end subroutine diffy_co_3D_opt

subroutine  diffy_oc_3D_opt_2( Bs, g, dy, u, dudy)

    integer(kind=ik), intent(in)    :: g, Bs(3)
    real(kind=rk), intent(in)       :: dy
    real(kind=rk), intent(in)       :: u(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)
    real(kind=rk), intent(out)      :: dudy(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

    integer                         :: i, j, k
    real(kind=rk)                   :: dy_inv

    ! - do not use ghost nodes
    ! - write loops explicitly,
    ! - use multiplication for dx
    ! - access array in column-major order

    dy_inv = 1.0_rk/(12.0_rk*dy)

    do k = g-1, Bs(3)+g+2
        do i = g-1, Bs(1)+g+2
            dudy(i,g+1,k) = ( - 18.0_rk*u(i,g+1,k) + 24.0_rk*u(i,g+2,k) - 6.0_rk*u(i,g+3,k) ) * dy_inv
            dudy(i,g+2,k) = ( - 6.0_rk*u(i,g+1,k)                       + 6.0_rk*u(i,g+3,k) ) * dy_inv
        end do
        do j = g+3, Bs(2)+g+2
            do i = g-1, Bs(1)+g+2
                dudy(i,j,k) = ( u(i,j-2,k) - 8.0_rk*u(i,j-1,k) + 8.0_rk*u(i,j+1,k) - u(i,j+2,k) ) * dy_inv
            end do
        end do
    end do

end subroutine diffy_oc_3D_opt_2

subroutine  diffy_co_3D_opt_2( Bs, g, dy, u, dudy)

    integer(kind=ik), intent(in)    :: g, Bs(3)
    real(kind=rk), intent(in)       :: dy
    real(kind=rk), intent(in)       :: u(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)
    real(kind=rk), intent(out)      :: dudy(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

    integer                         :: i, j, k
    real(kind=rk)                   :: dy_inv

    ! - do not use ghost nodes
    ! - write loops explicitly,
    ! - use multiplication for dx
    ! - access array in column-major order

    dy_inv = 1.0_rk/(12.0_rk*dy)

    do k = g-1, Bs(3)+g+2
        do j = g-1, Bs(2)+g-2
            do i = g-1, Bs(1)+g+2
                dudy(i,j,k) = ( u(i,j-2,k) - 8.0_rk*u(i,j-1,k) + 8.0_rk*u(i,j+1,k) - u(i,j+2,k) ) * dy_inv
            end do
        end do
        do i = g-1, Bs(1)+g+2
            dudy(i,Bs(2)+g-1,k) = ( - 6.0_rk*u(i,Bs(2)+g-2,k)                            + 6.0_rk*u(i,Bs(2)+g,k) ) * dy_inv
            dudy(i,Bs(2)+g,k)   = ( + 6.0_rk*u(i,Bs(2)+g-2,k) - 24.0_rk*u(i,Bs(2)+g-1,k) + 18.0_rk*u(i,Bs(2)+g,k) ) * dy_inv
        end do
    end do

end subroutine diffy_co_3D_opt_2


subroutine  diffz_oc_3D_opt( Bs, g, dz, u, dudz)

    integer(kind=ik), intent(in)    :: g, Bs(3)
    real(kind=rk), intent(in)       :: dz
    real(kind=rk), intent(in)       :: u(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)
    real(kind=rk), intent(out)      :: dudz(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

    integer                         :: i, j, k
    real(kind=rk)                   :: dz_inv

    ! - do not use ghost nodes
    ! - write loops explicitly,
    ! - use multiplication for dx
    ! - access array in column-major order

    dz_inv = 1.0_rk/(12.0_rk*dz)

    do j = g+1, Bs(2)+g
        do i = g+1, Bs(1)+g
            dudz(i,j,g+1) = ( - 18.0_rk*u(i,j,g+1) + 24.0_rk*u(i,j,g+2) - 6.0_rk*u(i,j,g+3) ) * dz_inv
            dudz(i,j,g+2) = ( - 6.0_rk*u(i,j,g+1)                       + 6.0_rk*u(i,j,g+3) ) * dz_inv
        end do
    end do
    do k = g+3, Bs(3)+g
        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                dudz(i,j,k) = ( u(i,j,k-2) - 8.0_rk*u(i,j,k-1) + 8.0_rk*u(i,j,k+1) - u(i,j,k+2) ) * dz_inv
            end do
        end do
    end do

end subroutine diffz_oc_3D_opt

subroutine  diffz_co_3D_opt( Bs, g, dz, u, dudz)

    integer(kind=ik), intent(in)    :: g, Bs(3)
    real(kind=rk), intent(in)       :: dz
    real(kind=rk), intent(in)       :: u(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)
    real(kind=rk), intent(out)      :: dudz(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

    integer                         :: i, j, k
    real(kind=rk)                   :: dz_inv

    ! - do not use ghost nodes
    ! - write loops explicitly,
    ! - use multiplication for dx
    ! - access array in column-major order

    dz_inv = 1.0_rk/(12.0_rk*dz)

    do k = g+1, Bs(3)+g-2
        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                dudz(i,j,k) = ( u(i,j,k-2) - 8.0_rk*u(i,j,k-1) + 8.0_rk*u(i,j,k+1) - u(i,j,k+2) ) * dz_inv
            end do
        end do
    end do
    do j = g+1, Bs(2)+g
        do i = g+1, Bs(1)+g
            dudz(i,j,Bs(3)+g-1) = ( - 6.0_rk*u(i,j,Bs(3)+g-2)                            + 6.0_rk*u(i,j,Bs(3)+g) ) * dz_inv
            dudz(i,j,Bs(3)+g)   = ( + 6.0_rk*u(i,j,Bs(3)+g-2) - 24.0_rk*u(i,j,Bs(3)+g-1) + 18.0_rk*u(i,j,Bs(3)+g) ) * dz_inv
        end do
    end do

end subroutine diffz_co_3D_opt

subroutine  diffz_oc_3D_opt_2( Bs, g, dz, u, dudz)

    integer(kind=ik), intent(in)    :: g, Bs(3)
    real(kind=rk), intent(in)       :: dz
    real(kind=rk), intent(in)       :: u(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)
    real(kind=rk), intent(out)      :: dudz(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

    integer                         :: i, j, k
    real(kind=rk)                   :: dz_inv

    ! - do not use ghost nodes
    ! - write loops explicitly,
    ! - use multiplication for dx
    ! - access array in column-major order

    dz_inv = 1.0_rk/(12.0_rk*dz)

    do j = g-1, Bs(2)+g+2
        do i = g-1, Bs(1)+g+2
            dudz(i,j,g+1) = ( - 18.0_rk*u(i,j,g+1) + 24.0_rk*u(i,j,g+2) - 6.0_rk*u(i,j,g+3) ) * dz_inv
            dudz(i,j,g+2) = ( - 6.0_rk*u(i,j,g+1)                       + 6.0_rk*u(i,j,g+3) ) * dz_inv
        end do
    end do
    do k = g+3, Bs(3)+g+2
        do j = g-1, Bs(2)+g+2
            do i = g-1, Bs(1)+g+2
                dudz(i,j,k) = ( u(i,j,k-2) - 8.0_rk*u(i,j,k-1) + 8.0_rk*u(i,j,k+1) - u(i,j,k+2) ) * dz_inv
            end do
        end do
    end do

end subroutine diffz_oc_3D_opt_2

subroutine  diffz_co_3D_opt_2( Bs, g, dz, u, dudz)

    integer(kind=ik), intent(in)    :: g, Bs(3)
    real(kind=rk), intent(in)       :: dz
    real(kind=rk), intent(in)       :: u(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)
    real(kind=rk), intent(out)      :: dudz(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

    integer                         :: i, j, k
    real(kind=rk)                   :: dz_inv

    ! - do not use ghost nodes
    ! - write loops explicitly,
    ! - use multiplication for dx
    ! - access array in column-major order

    dz_inv = 1.0_rk/(12.0_rk*dz)

    do k = g-1, Bs(3)+g-2
        do j = g-1, Bs(2)+g+2
            do i = g-1, Bs(1)+g+2
                dudz(i,j,k) = ( u(i,j,k-2) - 8.0_rk*u(i,j,k-1) + 8.0_rk*u(i,j,k+1) - u(i,j,k+2) ) * dz_inv
            end do
        end do
    end do
    do j = g-1, Bs(2)+g+2
        do i = g-1, Bs(1)+g+2
            dudz(i,j,Bs(3)+g-1) = ( - 6.0_rk*u(i,j,Bs(3)+g-2)                            + 6.0_rk*u(i,j,Bs(3)+g) ) * dz_inv
            dudz(i,j,Bs(3)+g)   = ( + 6.0_rk*u(i,j,Bs(3)+g-2) - 24.0_rk*u(i,j,Bs(3)+g-1) + 18.0_rk*u(i,j,Bs(3)+g) ) * dz_inv
        end do
    end do

end subroutine diffz_co_3D_opt_2

subroutine  diffxyz_c_3D_opt( Bs, g, dx, dy, dz, u, dudx, dudy, dudz)

    integer(kind=ik), intent(in)    :: g, Bs(3)
    real(kind=rk), intent(in)       :: dx, dy, dz
    real(kind=rk), intent(in)       :: u(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)
    real(kind=rk), intent(out)      :: dudx(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), dudy(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), dudz(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

    integer                         :: i, j, k
    real(kind=rk)                   :: dx_inv, dy_inv, dz_inv

    ! - do not use all ghost nodes (note: use only two ghost nodes to get correct second derivatives)
    ! - no one sided stencils necessary 
    ! - write loops explicitly, 
    ! - use multiplication for dx
    ! - access array in column-major order

    dx_inv = 1.0_rk/(12.0_rk*dx)
    dy_inv = 1.0_rk/(12.0_rk*dy)
    dz_inv = 1.0_rk/(12.0_rk*dz)

    do k = g-1, Bs(3)+g+2
        do j = g-1, Bs(2)+g+2
            do i = g-1, Bs(1)+g+2
                dudx(i,j,k) = ( u(i-2,j,k) - 8.0_rk*u(i-1,j,k) + 8.0_rk*u(i+1,j,k) - u(i+2,j,k) ) * dx_inv
                dudy(i,j,k) = ( u(i,j-2,k) - 8.0_rk*u(i,j-1,k) + 8.0_rk*u(i,j+1,k) - u(i,j+2,k) ) * dy_inv
                dudz(i,j,k) = ( u(i,j,k-2) - 8.0_rk*u(i,j,k-1) + 8.0_rk*u(i,j,k+1) - u(i,j,k+2) ) * dz_inv
            end do
        end do
    end do

end subroutine diffxyz_c_3D_opt

subroutine  diffx_c_3D_opt( Bs, g, dx, u, dudx)

    integer(kind=ik), intent(in)    :: g, Bs(3)
    real(kind=rk), intent(in)       :: dx
    real(kind=rk), intent(in)       :: u(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)
    real(kind=rk), intent(out)      :: dudx(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

    integer                         :: i, j, k
    real(kind=rk)                   :: dx_inv

    ! - do not use ghost nodes
    ! - no one sided stencils necessary 
    ! - write loops explicitly, 
    ! - use multiplication for dx
    ! - access array in column-major order

    dx_inv = 1.0_rk/(12.0_rk*dx)

    do k = g+1, Bs(3)+g
        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                dudx(i,j,k) = ( u(i-2,j,k) - 8.0_rk*u(i-1,j,k) + 8.0_rk*u(i+1,j,k) - u(i+2,j,k) ) * dx_inv
            end do
        end do
    end do

end subroutine diffx_c_3D_opt

subroutine  diffy_c_3D_opt( Bs, g, dy, u, dudy)

    integer(kind=ik), intent(in)    :: g, Bs(3)
    real(kind=rk), intent(in)       :: dy
    real(kind=rk), intent(in)       :: u(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)
    real(kind=rk), intent(out)      :: dudy(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

    integer                         :: i, j, k
    real(kind=rk)                   :: dy_inv

    ! - do not use ghost nodes
    ! - no one sided stencils necessary 
    ! - write loops explicitly, 
    ! - use multiplication for dx
    ! - access array in column-major order

    dy_inv = 1.0_rk/(12.0_rk*dy)

    do k = g+1, Bs(3)+g
        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                dudy(i,j,k) = ( u(i,j-2,k) - 8.0_rk*u(i,j-1,k) + 8.0_rk*u(i,j+1,k) - u(i,j+2,k) ) * dy_inv
            end do
        end do
    end do

end subroutine diffy_c_3D_opt

subroutine  diffz_c_3D_opt( Bs, g, dz, u, dudz)

    integer(kind=ik), intent(in)    :: g, Bs(3)
    real(kind=rk), intent(in)       :: dz
    real(kind=rk), intent(in)       :: u(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)
    real(kind=rk), intent(out)      :: dudz(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

    integer                         :: i, j, k
    real(kind=rk)                   :: dz_inv

    ! - do not use ghost nodes
    ! - no one sided stencils necessary 
    ! - write loops explicitly, 
    ! - use multiplication for dx
    ! - access array in column-major order

    dz_inv = 1.0_rk/(12.0_rk*dz)

    do k = g+1, Bs(3)+g
        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                dudz(i,j,k) = ( u(i,j,k-2) - 8.0_rk*u(i,j,k-1) + 8.0_rk*u(i,j,k+1) - u(i,j,k+2) ) * dz_inv
            end do
        end do
    end do

end subroutine diffz_c_3D_opt
