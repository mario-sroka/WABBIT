!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name RHS_2D_common_code.f90
!> \version 0.5
!> \author msr
!
!> \brief source code for 2D non reactive ns, prevent double code
!
! ********************************************************************************************

!---------------------------------------------------------------------------------------------
! main body

    ! derivatives
    ! u_x, u_y
    !---------------------------------------------------------------------------------------------
    call diff_wrapper_2D( Bs, g, u, 'periodic_u_xy', onesided, dx=dx, dy=dy, dudx=dummy, dudy=dummy2)

    do j = g+1, Bs(2)+g
        do i = g+1, Bs(1)+g
            rhs(i,j,1,UxF) = - 0.5_rk * rho(i,j) * ( u(i,j) * dummy(i,j) + v(i,j) * dummy2(i,j) )
        end do
    end do

    if (dissipation) then
        do j = g-1, Bs(2)+g+2
            do i = g-1, Bs(1)+g+2
                ! u_x
                tau11(i,j) = ( mu(i,j) * 2.0_rk +  mu_d - two_three * mu(i,j) ) * dummy(i,j)
                tau22(i,j) = ( mu_d - two_three * mu(i,j) ) * dummy(i,j)
                ! u_y
                tau12(i,j) = mu(i,j) * dummy2(i,j)
            end do
        end do
    end if

    ! v_x, v_y
    !---------------------------------------------------------------------------------------------
    call diff_wrapper_2D( Bs, g, v, 'periodic_u_xy', onesided, dx=dx, dy=dy, dudx=dummy, dudy=dummy2)

    do j = g+1, Bs(2)+g
        do i = g+1, Bs(1)+g
            rhs(i,j,1,UyF) = - 0.5_rk * rho(i,j) * ( u(i,j) * dummy(i,j) + v(i,j) * dummy2(i,j) )
        end do
    end do

    if (dissipation) then
        do j = g-1, Bs(2)+g+2
            do i = g-1, Bs(1)+g+2
                ! v_x
                tau12(i,j) = tau12(i,j) + mu(i,j) * dummy(i,j)
                ! v_y
                tau11(i,j) = tau11(i,j) + ( mu_d - two_three * mu(i,j) ) * dummy2(i,j)
                tau22(i,j) = tau22(i,j) + ( mu(i,j) * 2.0_rk + mu_d - two_three * mu(i,j) ) * dummy2(i,j)
            end do
        end do
    end if

    ! p_x, p_y
    !---------------------------------------------------------------------------------------------
    call diff_wrapper_2D( Bs, g, p, 'periodic_u_x_p', onesided, dx=dx, dudx=dummy)
    call diff_wrapper_2D( Bs, g, p, 'periodic_u_y_p', onesided, dy=dy, dudy=dummy2)

    do j = g+1, Bs(2)+g
        do i = g+1, Bs(1)+g
            rhs(i,j,1,UxF) = rhs(i,j,1,UxF) - dummy(i,j)
            rhs(i,j,1,UyF) = rhs(i,j,1,UyF) - dummy2(i,j)
            rhs(i,j,1,EF) =  rhs(i,j,1,EF) + u(i,j) * dummy(i,j) + v(i,j) * dummy2(i,j)
        end do
    end do

    ! friction
    if (dissipation) then

        ! Friction terms for Momentum equation = div(tau_i*)/(J*srho)
        ! tau11_x
        !---------------------------------------------------------------------------------------------
        call diff_wrapper_2D( Bs, g, tau11, 'periodic_u_x', onesided, dx=dx, dudx=dummy)
        
        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                rhs(i,j,1,UxF) = rhs(i,j,1,UxF) + dummy(i,j)
                rhs(i,j,1,EF) = rhs(i,j,1,EF) - u(i,j) * dummy(i,j)
            end do
        end do

        ! tau12_x
        !---------------------------------------------------------------------------------------------
        call diff_wrapper_2D( Bs, g, tau12, 'periodic_u_x', onesided, dx=dx, dudx=dummy)

        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                rhs(i,j,1,UyF) = rhs(i,j,1,UyF) + dummy(i,j)
                rhs(i,j,1,EF) = rhs(i,j,1,EF) - v(i,j) * dummy(i,j)
            end do
        end do

        ! tau12_y
        !---------------------------------------------------------------------------------------------
        call diff_wrapper_2D( Bs, g, tau12, 'periodic_u_y', onesided, dy=dy, dudy=dummy)

        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                rhs(i,j,1,UxF) = rhs(i,j,1,UxF) + dummy(i,j)
                rhs(i,j,1,EF) = rhs(i,j,1,EF) - u(i,j) * dummy(i,j)
            end do
        end do

        ! tau22_y
        !---------------------------------------------------------------------------------------------
        call diff_wrapper_2D( Bs, g, tau22, 'periodic_u_y', onesided, dy=dy, dudy=dummy)

        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                rhs(i,j,1,UyF) = rhs(i,j,1,UyF) + dummy(i,j)
                rhs(i,j,1,EF) = rhs(i,j,1,EF) - v(i,j) * dummy(i,j)
            end do
        end do

        ! Friction terms for the energy equation
        ! Heat Flux
        call diff_wrapper_2D( Bs, g, T, 'periodic_u_xy', onesided, dx=dx, dy=dy, dudx=dummy, dudy=dummy2)

        do j = g-1, Bs(2)+g+2
            do i = g-1, Bs(1)+g+2
                dummy3(i,j) = u(i,j)*tau11(i,j) + v(i,j)*tau12(i,j) + lambda(i,j) * dummy(i,j)
                dummy4(i,j) = u(i,j)*tau12(i,j) + v(i,j)*tau22(i,j) + lambda(i,j) * dummy2(i,j)
            end do
        end do

        call diff_wrapper_2D( Bs, g, dummy3, 'periodic_u_x', onesided, dx=dx, dudx=dummy)
        call diff_wrapper_2D( Bs, g, dummy4, 'periodic_u_y', onesided, dy=dy, dudy=dummy2)

        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                rhs(i,j,1,EF) = rhs(i,j,1,EF) + ( dummy(i,j) + dummy2(i,j) )
            end do
        end do

    end if

    ! EQUATIONS
    ! --------------------------------------------------------------------------------------------------------------
    ! RHS of equation of mass: J*srho*2 * srho_t = -div(rho*U_tilde)
    do j = g-1, Bs(2)+g+2
        do i = g-1, Bs(1)+g+2
            dummy(i,j)  = rho(i,j)*u(i,j)
            dummy2(i,j) = rho(i,j)*v(i,j)
        end do
    end do

    call diff_wrapper_2D( Bs, g, dummy,  'periodic_u_x', onesided, dx=dx, dudx=dummy3)
    call diff_wrapper_2D( Bs, g, dummy2, 'periodic_u_y', onesided, dy=dy, dudy=dummy4)

    do j = g+1, Bs(2)+g
        do i = g+1, Bs(1)+g
            rhs(i,j,1,rhoF) = (-dummy3(i,j) - dummy4(i,j) ) * 0.5_rk * phi1_inv(i,j)
        end do
    end do

    ! RHS of  momentum equation for u: sru_t = -1/2 * div(rho U_tilde u ) - 1/2 * (rho*U_tilde)*Du - Dp
    do j = g-1, Bs(2)+g+2
        do i = g-1, Bs(1)+g+2
            dummy(i,j)  = u(i,j)*rho(i,j)*u(i,j)
            dummy2(i,j) = v(i,j)*rho(i,j)*u(i,j)
        end do
    end do

    call diff_wrapper_2D( Bs, g, dummy,  'periodic_u_x', onesided, dx=dx, dudx=dummy3)
    call diff_wrapper_2D( Bs, g, dummy2, 'periodic_u_y', onesided, dy=dy, dudy=dummy4)

    do j = g+1, Bs(2)+g
        do i = g+1, Bs(1)+g
            rhs(i,j,1,UxF) = ( rhs(i,j,1,UxF) - 0.5_rk * ( dummy3(i,j) + dummy4(i,j) ) ) * phi1_inv(i,j)
        end do
    end do

    ! RHS of  momentum equation for v
    do j = g-1, Bs(2)+g+2
        do i = g-1, Bs(1)+g+2
            dummy(i,j)  = u(i,j)*rho(i,j)*v(i,j)
            dummy2(i,j) = v(i,j)*rho(i,j)*v(i,j)
        end do
    end do

    call diff_wrapper_2D( Bs, g, dummy,  'periodic_u_x', onesided, dx=dx, dudx=dummy3)
    call diff_wrapper_2D( Bs, g, dummy2, 'periodic_u_y', onesided, dy=dy, dudy=dummy4)

    do j = g+1, Bs(2)+g
        do i = g+1, Bs(1)+g
            rhs(i,j,1,UyF) = ( rhs(i,j,1,UyF) - 0.5_rk * ( dummy3(i,j) + dummy4(i,j) ) ) * phi1_inv(i,j)
        end do
    end do
