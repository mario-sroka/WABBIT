!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name RHS_3D_common_code.f90
!> \version 0.5
!> \author msr
!
!> \brief source code for 3D non reactive ns, prevent double code
!
! ********************************************************************************************

!---------------------------------------------------------------------------------------------
! main body

    ! derivatives
    ! u_x, u_y, u_z
    !---------------------------------------------------------------------------------------------
    call diff_wrapper_3D( Bs, g, u, 'periodic_u_xyz', onesided, dx=dx, dy=dy, dz=dz, dudx=dummy, dudy=dummy2, dudz=dummy3)

    do k = g+1, Bs(3)+g
        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                rhs(i,j,k,UxF) = - 0.5_rk * rho(i,j,k) * ( u(i,j,k) * dummy(i,j,k) + v(i,j,k) * dummy2(i,j,k) + w(i,j,k) * dummy3(i,j,k))
            end do
        end do
    end do

    if (dissipation) then
        do k = g-1, Bs(3)+g+2
            do j = g-1, Bs(2)+g+2
                do i = g-1, Bs(1)+g+2
                    ! u_x
                    tau11(i,j,k) = ( mu(i,j,k) * 2.0_rk +  mu_d - two_three * mu(i,j,k) ) * dummy(i,j,k)
                    tau22(i,j,k) = ( mu_d - two_three * mu(i,j,k) ) * dummy(i,j,k)
                    tau33(i,j,k) = ( mu_d - two_three * mu(i,j,k) ) * dummy(i,j,k)
                    ! u_y
                    tau12(i,j,k) = mu(i,j,k) * dummy2(i,j,k)
                    ! u_z
                    tau13(i,j,k) = mu(i,j,k) * dummy3(i,j,k)
                end do
            end do
        end do

    end if

    ! v_x, v_y, v_z
    !---------------------------------------------------------------------------------------------
    call diff_wrapper_3D( Bs, g, v, 'periodic_u_xyz', onesided, dx=dx, dy=dy, dz=dz, dudx=dummy, dudy=dummy2, dudz=dummy3)

    do k = g+1, Bs(3)+g
        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                rhs(i,j,k,UyF) = - 0.5_rk * rho(i,j,k) * ( u(i,j,k) * dummy(i,j,k) + v(i,j,k) * dummy2(i,j,k) + w(i,j,k) * dummy3(i,j,k))
            end do
        end do
    end do

    if (dissipation) then
        do k = g-1, Bs(3)+g+2
            do j = g-1, Bs(2)+g+2
                do i = g-1, Bs(1)+g+2
                    ! v_x
                    tau12(i,j,k) = tau12(i,j,k) + mu(i,j,k) * dummy(i,j,k)
                    ! v_y
                    tau11(i,j,k) = tau11(i,j,k) + ( mu_d - two_three * mu(i,j,k) ) * dummy2(i,j,k)
                    tau22(i,j,k) = tau22(i,j,k) + ( mu(i,j,k) * 2.0_rk + mu_d - two_three * mu(i,j,k) ) * dummy2(i,j,k)
                    tau33(i,j,k) = tau33(i,j,k) + ( mu_d - two_three * mu(i,j,k) ) * dummy2(i,j,k)
                    ! v_z
                    tau23(i,j,k) = mu(i,j,k) * dummy3(i,j,k)
                end do
            end do
        end do
    end if

    ! w_x, w_y, w_z
    !---------------------------------------------------------------------------------------------
    call diff_wrapper_3D( Bs, g, w, 'periodic_u_xyz', onesided, dx=dx, dy=dy, dz=dz, dudx=dummy, dudy=dummy2, dudz=dummy3)

    do k = g+1, Bs(3)+g
        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                rhs(i,j,k,UzF) = - 0.5_rk * rho(i,j,k) * ( u(i,j,k) * dummy(i,j,k) + v(i,j,k) * dummy2(i,j,k) + w(i,j,k) * dummy3(i,j,k))
            end do
        end do
    end do

    if (dissipation) then
        do k = g-1, Bs(3)+g+2
            do j = g-1, Bs(2)+g+2
                do i = g-1, Bs(1)+g+2
                    ! w_x
                    tau13(i,j,k) = tau13(i,j,k) + mu(i,j,k) * dummy(i,j,k)
                    ! w_z
                    tau11(i,j,k) = tau11(i,j,k) + ( mu_d - two_three * mu(i,j,k) ) * dummy3(i,j,k)
                    tau22(i,j,k) = tau22(i,j,k) + ( mu_d - two_three * mu(i,j,k) ) * dummy3(i,j,k)
                    tau33(i,j,k) = tau33(i,j,k) + ( mu(i,j,k) * 2.0_rk + mu_d - two_three * mu(i,j,k) ) * dummy3(i,j,k)
                    ! w_y
                    tau23(i,j,k) = tau23(i,j,k) + mu(i,j,k) * dummy2(i,j,k)
                end do
            end do
        end do
    end if

    ! p_x, p_y, p_z
    !---------------------------------------------------------------------------------------------
    call diff_wrapper_3D( Bs, g, p, 'periodic_u_x_p', onesided, dx=dx, dudx=dummy)
    call diff_wrapper_3D( Bs, g, p, 'periodic_u_y_p', onesided, dy=dy, dudy=dummy2)
    call diff_wrapper_3D( Bs, g, p, 'periodic_u_z_p', onesided, dz=dz, dudz=dummy3)

    do k = g+1, Bs(3)+g
        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                rhs(i,j,k,UxF) = rhs(i,j,k,UxF) - dummy(i,j,k)
                rhs(i,j,k,UyF) = rhs(i,j,k,UyF) - dummy2(i,j,k)
                rhs(i,j,k,UzF) = rhs(i,j,k,UzF) - dummy3(i,j,k)
                rhs(i,j,k,EF) =  rhs(i,j,k,EF) + u(i,j,k) * dummy(i,j,k) + v(i,j,k) * dummy2(i,j,k) + w(i,j,k) * dummy3(i,j,k)
            end do
        end do
    end do

    ! friction
    if (dissipation) then

        ! Friction terms for Momentum equation = div(tau_i*)/(J*srho)
        ! tau11_x
        !---------------------------------------------------------------------------------------------
        call diff_wrapper_3D( Bs, g, tau11, 'periodic_u_x', onesided, dx=dx, dudx=dummy)
        
        do k = g+1, Bs(3)+g
            do j = g+1, Bs(2)+g
                do i = g+1, Bs(1)+g
                    rhs(i,j,k,UxF) = rhs(i,j,k,UxF) + dummy(i,j,k)
                    rhs(i,j,k,EF) = rhs(i,j,k,EF) - u(i,j,k) * dummy(i,j,k)
                end do
            end do
        end do

        ! tau12_x
        !---------------------------------------------------------------------------------------------
        call diff_wrapper_3D( Bs, g, tau12, 'periodic_u_x', onesided, dx=dx, dudx=dummy)

        do k = g+1, Bs(3)+g
            do j = g+1, Bs(2)+g
                do i = g+1, Bs(1)+g
                    rhs(i,j,k,UyF) = rhs(i,j,k,UyF) + dummy(i,j,k)
                    rhs(i,j,k,EF) = rhs(i,j,k,EF) - v(i,j,k) * dummy(i,j,k)
                end do
            end do
        end do

        ! tau12_y
        !---------------------------------------------------------------------------------------------
        call diff_wrapper_3D( Bs, g, tau12, 'periodic_u_y', onesided, dy=dy, dudy=dummy)

        do k = g+1, Bs(3)+g
            do j = g+1, Bs(2)+g
                do i = g+1, Bs(1)+g
                    rhs(i,j,k,UxF) = rhs(i,j,k,UxF) + dummy(i,j,k)
                    rhs(i,j,k,EF) = rhs(i,j,k,EF) - u(i,j,k) * dummy(i,j,k)
                end do
            end do
        end do

        ! tau13_x
        !---------------------------------------------------------------------------------------------
        call diff_wrapper_3D( Bs, g, tau13, 'periodic_u_x', onesided, dx=dx, dudx=dummy)

        do k = g+1, Bs(3)+g
            do j = g+1, Bs(2)+g
                do i = g+1, Bs(1)+g
                    rhs(i,j,k,UzF) = rhs(i,j,k,UzF) + dummy(i,j,k)
                    rhs(i,j,k,EF) = rhs(i,j,k,EF) - w(i,j,k) * dummy(i,j,k)
                end do
            end do
        end do

        ! tau13_z
        !---------------------------------------------------------------------------------------------
        call diff_wrapper_3D( Bs, g, tau13, 'periodic_u_z', onesided, dz=dz, dudz=dummy)

        do k = g+1, Bs(3)+g
            do j = g+1, Bs(2)+g
                do i = g+1, Bs(1)+g
                    rhs(i,j,k,UxF) = rhs(i,j,k,UxF) + dummy(i,j,k)
                    rhs(i,j,k,EF) = rhs(i,j,k,EF) - u(i,j,k) * dummy(i,j,k)
                end do
            end do
        end do

        ! tau22_y
        !---------------------------------------------------------------------------------------------
        call diff_wrapper_3D( Bs, g, tau22, 'periodic_u_y', onesided, dy=dy, dudy=dummy)

        do k = g+1, Bs(3)+g
            do j = g+1, Bs(2)+g
                do i = g+1, Bs(1)+g
                    rhs(i,j,k,UyF) = rhs(i,j,k,UyF) + dummy(i,j,k)
                    rhs(i,j,k,EF) = rhs(i,j,k,EF) - v(i,j,k) * dummy(i,j,k)
                end do
            end do
        end do

        ! tau23_y
        !---------------------------------------------------------------------------------------------
        call diff_wrapper_3D( Bs, g, tau23, 'periodic_u_y', onesided, dy=dy, dudy=dummy)

        do k = g+1, Bs(3)+g
            do j = g+1, Bs(2)+g
                do i = g+1, Bs(1)+g
                    rhs(i,j,k,UzF) = rhs(i,j,k,UzF) + dummy(i,j,k)
                    rhs(i,j,k,EF) = rhs(i,j,k,EF) - w(i,j,k) * dummy(i,j,k)
                end do
            end do
        end do

        ! tau23_z
        !---------------------------------------------------------------------------------------------
        call diff_wrapper_3D( Bs, g, tau23, 'periodic_u_z', onesided, dz=dz, dudz=dummy)

        do k = g+1, Bs(3)+g
            do j = g+1, Bs(2)+g
                do i = g+1, Bs(1)+g
                    rhs(i,j,k,UyF) = rhs(i,j,k,UyF) + dummy(i,j,k)
                    rhs(i,j,k,EF) = rhs(i,j,k,EF) - v(i,j,k) * dummy(i,j,k)
                end do
            end do
        end do

        ! tau33_z
        !---------------------------------------------------------------------------------------------
        call diff_wrapper_3D( Bs, g, tau33, 'periodic_u_z', onesided, dz=dz, dudz=dummy)

        do k = g+1, Bs(3)+g
            do j = g+1, Bs(2)+g
                do i = g+1, Bs(1)+g
                    rhs(i,j,k,UzF) = rhs(i,j,k,UzF) + dummy(i,j,k)
                    rhs(i,j,k,EF) = rhs(i,j,k,EF) - w(i,j,k) * dummy(i,j,k)
                end do
            end do
        end do

        ! Friction terms for the energy equation
        ! Heat Flux
        call diff_wrapper_3D( Bs, g, T, 'periodic_u_xyz', onesided, dx=dx, dy=dy, dz=dz, dudx=dummy, dudy=dummy2, dudz=dummy3)

        do k = g-1, Bs(3)+g+2
            do j = g-1, Bs(2)+g+2
                do i = g-1, Bs(1)+g+2
                    dummy4(i,j,k) = u(i,j,k)*tau11(i,j,k) + v(i,j,k)*tau12(i,j,k) + w(i,j,k)*tau13(i,j,k) + lambda(i,j,k) * dummy(i,j,k)
                    dummy5(i,j,k) = u(i,j,k)*tau12(i,j,k) + v(i,j,k)*tau22(i,j,k) + w(i,j,k)*tau23(i,j,k) + lambda(i,j,k) * dummy2(i,j,k)
                    dummy6(i,j,k) = u(i,j,k)*tau13(i,j,k) + v(i,j,k)*tau23(i,j,k) + w(i,j,k)*tau33(i,j,k) + lambda(i,j,k) * dummy3(i,j,k)
                end do
            end do
        end do

        call diff_wrapper_3D( Bs, g, dummy4, 'periodic_u_x', onesided, dx=dx,  dudx=dummy)
        call diff_wrapper_3D( Bs, g, dummy5, 'periodic_u_y', onesided, dy=dy, dudy=dummy2)
        call diff_wrapper_3D( Bs, g, dummy6, 'periodic_u_z', onesided, dz=dz, dudz=dummy3)

        do k = g+1, Bs(3)+g
            do j = g+1, Bs(2)+g
                do i = g+1, Bs(1)+g
                    rhs(i,j,k,EF) = rhs(i,j,k,EF) + ( dummy(i,j,k) + dummy2(i,j,k) + dummy3(i,j,k) )
                end do
            end do
        end do

    end if

    ! EQUATIONS
    ! --------------------------------------------------------------------------------------------------------------
    ! RHS of equation of mass: J*srho*2 * srho_t = -div(rho*U_tilde)
    do k = g-1, Bs(3)+g+2
        do j = g-1, Bs(2)+g+2
            do i = g-1, Bs(1)+g+2
                dummy(i,j,k)  = rho(i,j,k)*u(i,j,k)
                dummy2(i,j,k) = rho(i,j,k)*v(i,j,k)
                dummy3(i,j,k) = rho(i,j,k)*w(i,j,k)
            end do
        end do
    end do
    call diff_wrapper_3D( Bs, g, dummy,  'periodic_u_x', onesided, dx=dx, dudx=dummy4)
    call diff_wrapper_3D( Bs, g, dummy2, 'periodic_u_y', onesided, dy=dy, dudy=dummy5)
    call diff_wrapper_3D( Bs, g, dummy3, 'periodic_u_z', onesided, dz=dz, dudz=dummy6)

    do k = g+1, Bs(3)+g
        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                rhs(i,j,k,rhoF) = (-dummy4(i,j,k) - dummy5(i,j,k) - dummy6(i,j,k)) * 0.5_rk * phi1_inv(i,j,k)
            end do
        end do
    end do

    ! RHS of  momentum equation for u: sru_t = -1/2 * div(rho U_tilde u ) - 1/2 * (rho*U_tilde)*Du - Dp
    do k = g-1, Bs(3)+g+2
        do j = g-1, Bs(2)+g+2
            do i = g-1, Bs(1)+g+2
                dummy(i,j,k)  = u(i,j,k)*rho(i,j,k)*u(i,j,k)
                dummy2(i,j,k) = v(i,j,k)*rho(i,j,k)*u(i,j,k)
                dummy3(i,j,k) = w(i,j,k)*rho(i,j,k)*u(i,j,k)
            end do
        end do
    end do
    call diff_wrapper_3D( Bs, g, dummy,  'periodic_u_x', onesided, dx=dx, dudx=dummy4)
    call diff_wrapper_3D( Bs, g, dummy2, 'periodic_u_y', onesided, dy=dy, dudy=dummy5)
    call diff_wrapper_3D( Bs, g, dummy3, 'periodic_u_z', onesided, dz=dz, dudz=dummy6)

    do k = g+1, Bs(3)+g
        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                rhs(i,j,k,UxF) = ( rhs(i,j,k,UxF) - 0.5_rk * ( dummy4(i,j,k) + dummy5(i,j,k) + dummy6(i,j,k) ) ) * phi1_inv(i,j,k)
            end do
        end do
    end do

    ! RHS of  momentum equation for v
    do k = g-1, Bs(3)+g+2
        do j = g-1, Bs(2)+g+2
            do i = g-1, Bs(1)+g+2
                dummy(i,j,k)  = u(i,j,k)*rho(i,j,k)*v(i,j,k)
                dummy2(i,j,k) = v(i,j,k)*rho(i,j,k)*v(i,j,k)
                dummy3(i,j,k) = w(i,j,k)*rho(i,j,k)*v(i,j,k)
            end do
        end do
    end do
    call diff_wrapper_3D( Bs, g, dummy,  'periodic_u_x', onesided, dx=dx, dudx=dummy4)
    call diff_wrapper_3D( Bs, g, dummy2, 'periodic_u_y', onesided, dy=dy, dudy=dummy5)
    call diff_wrapper_3D( Bs, g, dummy3, 'periodic_u_z', onesided, dz=dz, dudz=dummy6)

    do k = g+1, Bs(3)+g
        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                rhs(i,j,k,UyF) = ( rhs(i,j,k,UyF) - 0.5_rk * ( dummy4(i,j,k) + dummy5(i,j,k) + dummy6(i,j,k)) ) * phi1_inv(i,j,k)
            end do
        end do
    end do

    ! RHS of  momentum equation for w
    do k = g-1, Bs(3)+g+2
        do j = g-1, Bs(2)+g+2
            do i = g-1, Bs(1)+g+2
                dummy(i,j,k)  = u(i,j,k)*rho(i,j,k)*w(i,j,k)
                dummy2(i,j,k) = v(i,j,k)*rho(i,j,k)*w(i,j,k)
                dummy3(i,j,k) = w(i,j,k)*rho(i,j,k)*w(i,j,k)
            end do
        end do
    end do
    call diff_wrapper_3D( Bs, g, dummy,  'periodic_u_x', onesided, dx=dx, dudx=dummy4)
    call diff_wrapper_3D( Bs, g, dummy2, 'periodic_u_y', onesided, dy=dy, dudy=dummy5)
    call diff_wrapper_3D( Bs, g, dummy3, 'periodic_u_z', onesided, dz=dz, dudz=dummy6)

    do k = g+1, Bs(3)+g
        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                rhs(i,j,k,UzF) = ( rhs(i,j,k,UzF) - 0.5_rk * ( dummy4(i,j,k) + dummy5(i,j,k) + dummy6(i,j,k)) ) * phi1_inv(i,j,k)
            end do
        end do
    end do
