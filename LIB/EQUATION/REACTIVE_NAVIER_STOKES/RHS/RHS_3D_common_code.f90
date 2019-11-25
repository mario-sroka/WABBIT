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
    call diff_wrapper_3D( Bs, g, u, 'diff_x_c_gm', onesided, dx=dx, dudx=dummy)
    call diff_wrapper_3D( Bs, g, u, 'diff_y_c_gm', onesided, dy=dy, dudy=dummy2)
    call diff_wrapper_3D( Bs, g, u, 'diff_z_c_gm', onesided, dz=dz, dudz=dummy3)

    ! save dummy
    dummy5(gp:Bp(1), gp:Bp(2), gp:Bp(3)) = dummy(gp:Bp(1), gp:Bp(2), gp:Bp(3)) * u(gp:Bp(1), gp:Bp(2), gp:Bp(3))

    if (dissipation) then

        ! tauij = mu * uk
        tau11(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = 2.0_rk * mu(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * dummy(gm:Bm(1), gm:Bm(2), gm:Bm(3))
        tau12(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = mu(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * dummy2(gm:Bm(1), gm:Bm(2), gm:Bm(3))
        tau13(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = mu(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * dummy3(gm:Bm(1), gm:Bm(2), gm:Bm(3))

        ! (mu_d - 2/3 mu) * u_x
        dummy(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = dummy6(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * dummy(gm:Bm(1), gm:Bm(2), gm:Bm(3))
        ! tauii
        tau11(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = tau11(gm:Bm(1), gm:Bm(2), gm:Bm(3)) + dummy(gm:Bm(1), gm:Bm(2), gm:Bm(3))
        tau22(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = dummy(gm:Bm(1), gm:Bm(2), gm:Bm(3))
        tau33(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = dummy(gm:Bm(1), gm:Bm(2), gm:Bm(3))

    end if

    ! Ux 
    dummy2(gp:Bp(1), gp:Bp(2), gp:Bp(3)) = dummy5(gp:Bp(1), gp:Bp(2), gp:Bp(3)) + dummy2(gp:Bp(1), gp:Bp(2), gp:Bp(3)) * v(gp:Bp(1), gp:Bp(2), gp:Bp(3))
    dummy3(gp:Bp(1), gp:Bp(2), gp:Bp(3)) = dummy2(gp:Bp(1), gp:Bp(2), gp:Bp(3)) + dummy3(gp:Bp(1), gp:Bp(2), gp:Bp(3)) * w(gp:Bp(1), gp:Bp(2), gp:Bp(3))

    rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UxF) = - 0.5_rk * rho(gp:Bp(1), gp:Bp(2), gp:Bp(3)) * dummy3(gp:Bp(1), gp:Bp(2), gp:Bp(3))

    ! v_x, v_y, v_z
    !---------------------------------------------------------------------------------------------
    call diff_wrapper_3D( Bs, g, v, 'diff_x_c_gm', onesided, dx=dx, dudx=dummy)
    call diff_wrapper_3D( Bs, g, v, 'diff_y_c_gm', onesided, dy=dy, dudy=dummy2)
    call diff_wrapper_3D( Bs, g, v, 'diff_z_c_gm', onesided, dz=dz, dudz=dummy3)

    ! save dummy2
    dummy5(gp:Bp(1), gp:Bp(2), gp:Bp(3)) = dummy2(gp:Bp(1), gp:Bp(2), gp:Bp(3)) * v(gp:Bp(1), gp:Bp(2), gp:Bp(3))

    if (dissipation) then

        ! tauij = mu * uk
        tau12(gm:Bm(1), gm:Bm(2), gm:Bm(3))  = tau12(gm:Bm(1), gm:Bm(2), gm:Bm(3)) + mu(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * dummy(gm:Bm(1), gm:Bm(2), gm:Bm(3))
        tau22(gm:Bm(1), gm:Bm(2), gm:Bm(3))  = tau22(gm:Bm(1), gm:Bm(2), gm:Bm(3)) + 2.0_rk * mu(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * dummy2(gm:Bm(1), gm:Bm(2), gm:Bm(3))
        tau23(gm:Bm(1), gm:Bm(2), gm:Bm(3))  = mu(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * dummy3(gm:Bm(1), gm:Bm(2), gm:Bm(3))

        ! (mu_d - 2/3 mu) * v_y
        dummy2(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = dummy6(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * dummy2(gm:Bm(1), gm:Bm(2), gm:Bm(3))
        ! tauij
        tau11(gm:Bm(1), gm:Bm(2), gm:Bm(3))  = tau11(gm:Bm(1), gm:Bm(2), gm:Bm(3)) + dummy2(gm:Bm(1), gm:Bm(2), gm:Bm(3))
        tau22(gm:Bm(1), gm:Bm(2), gm:Bm(3))  = tau22(gm:Bm(1), gm:Bm(2), gm:Bm(3)) + dummy2(gm:Bm(1), gm:Bm(2), gm:Bm(3))
        tau33(gm:Bm(1), gm:Bm(2), gm:Bm(3))  = tau33(gm:Bm(1), gm:Bm(2), gm:Bm(3)) + dummy2(gm:Bm(1), gm:Bm(2), gm:Bm(3))

    end if

    ! Uy 
    dummy( gp:Bp(1), gp:Bp(2), gp:Bp(3)) = dummy5(gp:Bp(1), gp:Bp(2), gp:Bp(3)) + dummy( gp:Bp(1), gp:Bp(2), gp:Bp(3)) * u(gp:Bp(1), gp:Bp(2), gp:Bp(3))
    dummy3(gp:Bp(1), gp:Bp(2), gp:Bp(3)) = dummy( gp:Bp(1), gp:Bp(2), gp:Bp(3)) + dummy3(gp:Bp(1), gp:Bp(2), gp:Bp(3)) * w(gp:Bp(1), gp:Bp(2), gp:Bp(3))

    rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UyF) = - 0.5_rk * rho(gp:Bp(1), gp:Bp(2), gp:Bp(3)) * dummy3(gp:Bp(1), gp:Bp(2), gp:Bp(3))


    ! w_x, w_y, w_z
    !---------------------------------------------------------------------------------------------
    call diff_wrapper_3D( Bs, g, w, 'diff_x_c_gm', onesided, dx=dx, dudx=dummy)
    call diff_wrapper_3D( Bs, g, w, 'diff_y_c_gm', onesided, dy=dy, dudy=dummy2)
    call diff_wrapper_3D( Bs, g, w, 'diff_z_c_gm', onesided, dz=dz, dudz=dummy3)

    ! save dummy3
    dummy5(gp:Bp(1), gp:Bp(2), gp:Bp(3)) = dummy3(gp:Bp(1), gp:Bp(2), gp:Bp(3)) * w(gp:Bp(1), gp:Bp(2), gp:Bp(3))

    if (dissipation) then

        ! tauij = mu * uk
        tau13(gm:Bm(1), gm:Bm(2), gm:Bm(3))  = tau13(gm:Bm(1), gm:Bm(2), gm:Bm(3)) + mu(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * dummy(gm:Bm(1), gm:Bm(2), gm:Bm(3))
        tau23(gm:Bm(1), gm:Bm(2), gm:Bm(3))  = tau23(gm:Bm(1), gm:Bm(2), gm:Bm(3)) + mu(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * dummy2(gm:Bm(1), gm:Bm(2), gm:Bm(3))
        tau33(gm:Bm(1), gm:Bm(2), gm:Bm(3))  = tau33(gm:Bm(1), gm:Bm(2), gm:Bm(3)) + 2.0_rk * mu(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * dummy3(gm:Bm(1), gm:Bm(2), gm:Bm(3))

        ! (mu_d - 2/3 mu) * w_z
        dummy3(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = dummy6(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * dummy3(gm:Bm(1), gm:Bm(2), gm:Bm(3))
        ! tau
        tau11(gm:Bm(1), gm:Bm(2), gm:Bm(3))  = tau11(gm:Bm(1), gm:Bm(2), gm:Bm(3)) + dummy3(gm:Bm(1), gm:Bm(2), gm:Bm(3))
        tau22(gm:Bm(1), gm:Bm(2), gm:Bm(3))  = tau22(gm:Bm(1), gm:Bm(2), gm:Bm(3)) + dummy3(gm:Bm(1), gm:Bm(2), gm:Bm(3))
        tau33(gm:Bm(1), gm:Bm(2), gm:Bm(3))  = tau33(gm:Bm(1), gm:Bm(2), gm:Bm(3)) + dummy3(gm:Bm(1), gm:Bm(2), gm:Bm(3))

    end if

    ! Uz 
    dummy( gp:Bp(1), gp:Bp(2), gp:Bp(3)) = dummy5(gp:Bp(1), gp:Bp(2), gp:Bp(3)) + dummy( gp:Bp(1), gp:Bp(2), gp:Bp(3)) * u(gp:Bp(1), gp:Bp(2), gp:Bp(3))
    dummy2(gp:Bp(1), gp:Bp(2), gp:Bp(3)) = dummy( gp:Bp(1), gp:Bp(2), gp:Bp(3)) + dummy2(gp:Bp(1), gp:Bp(2), gp:Bp(3)) * v(gp:Bp(1), gp:Bp(2), gp:Bp(3))

    rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UzF) = - 0.5_rk * rho(gp:Bp(1), gp:Bp(2), gp:Bp(3)) * dummy2(gp:Bp(1), gp:Bp(2), gp:Bp(3))

    ! p_x, p_y, p_z
    !---------------------------------------------------------------------------------------------
    call diff_wrapper_3D( Bs, g, p, 'diff_x_c_gp', onesided, dx=dx, dudx=dummy)
    call diff_wrapper_3D( Bs, g, p, 'diff_y_c_gp', onesided, dy=dy, dudy=dummy2)
    call diff_wrapper_3D( Bs, g, p, 'diff_z_c_gp', onesided, dz=dz, dudz=dummy3)

    ! momentum
    rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UxF) = rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UxF) - dummy( gp:Bp(1), gp:Bp(2), gp:Bp(3))
    rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UyF) = rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UyF) - dummy2(gp:Bp(1), gp:Bp(2), gp:Bp(3))
    rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UzF) = rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UzF) - dummy3(gp:Bp(1), gp:Bp(2), gp:Bp(3))   

    ! px*u + py*v + pz*w
    dummy( gp:Bp(1), gp:Bp(2), gp:Bp(3)) = dummy( gp:Bp(1), gp:Bp(2), gp:Bp(3)) * u(gp:Bp(1), gp:Bp(2), gp:Bp(3)) 
    dummy2(gp:Bp(1), gp:Bp(2), gp:Bp(3)) = dummy( gp:Bp(1), gp:Bp(2), gp:Bp(3)) + dummy2(gp:Bp(1), gp:Bp(2), gp:Bp(3)) * v(gp:Bp(1), gp:Bp(2), gp:Bp(3))
    dummy3(gp:Bp(1), gp:Bp(2), gp:Bp(3)) = dummy2(gp:Bp(1), gp:Bp(2), gp:Bp(3)) + dummy3(gp:Bp(1), gp:Bp(2), gp:Bp(3)) * w(gp:Bp(1), gp:Bp(2), gp:Bp(3))

    ! energy
    rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), EF) = rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), EF) + dummy3(gp:Bp(1), gp:Bp(2), gp:Bp(3))

    ! friction
    if (dissipation) then

        ! Friction terms for Momentum equation = div(tau_i*)/(J*srho)
        ! tau11_x
        !---------------------------------------------------------------------------------------------
        call diff_wrapper_3D( Bs, g, tau11, 'diff_x_c_gp', onesided, dx=dx, dudx=dummy)

        rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UxF) = rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UxF) + dummy( gp:Bp(1), gp:Bp(2), gp:Bp(3))
        rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), EF)  = rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), EF)  - dummy(gp:Bp(1), gp:Bp(2), gp:Bp(3)) * u(gp:Bp(1), gp:Bp(2), gp:Bp(3))

        ! tau12_x
        !---------------------------------------------------------------------------------------------
        call diff_wrapper_3D( Bs, g, tau12, 'diff_x_c_gp', onesided, dx=dx, dudx=dummy)

        rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UyF) = rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UyF) + dummy( gp:Bp(1), gp:Bp(2), gp:Bp(3))
        rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), EF)  = rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), EF)  - dummy(gp:Bp(1), gp:Bp(2), gp:Bp(3)) * v(gp:Bp(1), gp:Bp(2), gp:Bp(3))

        ! tau12_y
        !---------------------------------------------------------------------------------------------
        call diff_wrapper_3D( Bs, g, tau12, 'diff_y_c_gp', onesided, dy=dy, dudy=dummy)

        rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UxF) = rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UxF) + dummy( gp:Bp(1), gp:Bp(2), gp:Bp(3))
        rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), EF)  = rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), EF)  - dummy(gp:Bp(1), gp:Bp(2), gp:Bp(3)) * u(gp:Bp(1), gp:Bp(2), gp:Bp(3))

        ! tau13_x
        !---------------------------------------------------------------------------------------------
        call diff_wrapper_3D( Bs, g, tau13, 'diff_x_c_gp', onesided, dx=dx, dudx=dummy)

        rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UzF) = rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UzF) + dummy( gp:Bp(1), gp:Bp(2), gp:Bp(3))
        rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), EF)  = rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), EF)  - dummy(gp:Bp(1), gp:Bp(2), gp:Bp(3)) * w(gp:Bp(1), gp:Bp(2), gp:Bp(3))

        ! tau13_z
        !---------------------------------------------------------------------------------------------
        call diff_wrapper_3D( Bs, g, tau13, 'diff_z_c_gp', onesided, dz=dz, dudz=dummy)

        rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UxF) = rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UxF) + dummy( gp:Bp(1), gp:Bp(2), gp:Bp(3))
        rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), EF)  = rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), EF)  - dummy(gp:Bp(1), gp:Bp(2), gp:Bp(3)) * u(gp:Bp(1), gp:Bp(2), gp:Bp(3))

        ! tau22_y
        !---------------------------------------------------------------------------------------------
        call diff_wrapper_3D( Bs, g, tau22, 'diff_y_c_gp', onesided, dy=dy, dudy=dummy)

        rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UyF) = rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UyF) + dummy( gp:Bp(1), gp:Bp(2), gp:Bp(3))
        rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), EF)  = rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), EF)  - dummy(gp:Bp(1), gp:Bp(2), gp:Bp(3)) * v(gp:Bp(1), gp:Bp(2), gp:Bp(3))

        ! tau23_y
        !---------------------------------------------------------------------------------------------
        call diff_wrapper_3D( Bs, g, tau23, 'diff_y_c_gp', onesided, dy=dy, dudy=dummy)

        rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UzF) = rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UzF) + dummy( gp:Bp(1), gp:Bp(2), gp:Bp(3))
        rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), EF)  = rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), EF)  - dummy(gp:Bp(1), gp:Bp(2), gp:Bp(3)) * w(gp:Bp(1), gp:Bp(2), gp:Bp(3))

        ! tau23_z
        !---------------------------------------------------------------------------------------------
        call diff_wrapper_3D( Bs, g, tau23, 'diff_z_c_gp', onesided, dz=dz, dudz=dummy)

        rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UyF) = rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UyF) + dummy( gp:Bp(1), gp:Bp(2), gp:Bp(3))
        rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), EF)  = rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), EF)  - dummy(gp:Bp(1), gp:Bp(2), gp:Bp(3)) * v(gp:Bp(1), gp:Bp(2), gp:Bp(3))

        ! tau33_z
        !---------------------------------------------------------------------------------------------
        call diff_wrapper_3D( Bs, g, tau33, 'diff_z_c_gp', onesided, dz=dz, dudz=dummy)

        rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UzF) = rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UzF) + dummy( gp:Bp(1), gp:Bp(2), gp:Bp(3))
        rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), EF)  = rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), EF)  - dummy(gp:Bp(1), gp:Bp(2), gp:Bp(3)) * w(gp:Bp(1), gp:Bp(2), gp:Bp(3))

        ! Friction terms for the energy equation
        ! Heat Flux
        !---------------------------------------------------------------------------------------------
        call diff_wrapper_3D( Bs, g, T, 'diff_x_c_gm', onesided, dx=dx, dudx=dummy)
        call diff_wrapper_3D( Bs, g, T, 'diff_y_c_gm', onesided, dy=dy, dudy=dummy2)
        call diff_wrapper_3D( Bs, g, T, 'diff_z_c_gm', onesided, dz=dz, dudz=dummy3)

        tau11(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = u(     gm:Bm(1), gm:Bm(2), gm:Bm(3)) * tau11(gm:Bm(1), gm:Bm(2), gm:Bm(3))
        tau11(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = tau11( gm:Bm(1), gm:Bm(2), gm:Bm(3)) + lambda(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * dummy(gm:Bm(1), gm:Bm(2), gm:Bm(3))
        tau11(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = tau11( gm:Bm(1), gm:Bm(2), gm:Bm(3)) + v(     gm:Bm(1), gm:Bm(2), gm:Bm(3)) * tau12(gm:Bm(1), gm:Bm(2), gm:Bm(3))
        tau11(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = tau11( gm:Bm(1), gm:Bm(2), gm:Bm(3)) + w(     gm:Bm(1), gm:Bm(2), gm:Bm(3)) * tau13(gm:Bm(1), gm:Bm(2), gm:Bm(3))

        dummy4(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = tau11( gm:Bm(1), gm:Bm(2), gm:Bm(3))

        !-----------------------------------------------------------------------------------------
        tau22(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = v(    gm:Bm(1), gm:Bm(2), gm:Bm(3)) * tau22( gm:Bm(1), gm:Bm(2), gm:Bm(3))
        tau22(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = tau22(gm:Bm(1), gm:Bm(2), gm:Bm(3)) + lambda(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * dummy2(gm:Bm(1), gm:Bm(2), gm:Bm(3))
        tau22(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = tau22(gm:Bm(1), gm:Bm(2), gm:Bm(3)) + u(     gm:Bm(1), gm:Bm(2), gm:Bm(3)) * tau12( gm:Bm(1), gm:Bm(2), gm:Bm(3))
        tau22(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = tau22(gm:Bm(1), gm:Bm(2), gm:Bm(3)) + w(     gm:Bm(1), gm:Bm(2), gm:Bm(3)) * tau23( gm:Bm(1), gm:Bm(2), gm:Bm(3))

        dummy5(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = tau22( gm:Bm(1), gm:Bm(2), gm:Bm(3)) 

        !-----------------------------------------------------------------------------------------
        tau33(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = w(    gm:Bm(1), gm:Bm(2), gm:Bm(3)) * tau33( gm:Bm(1), gm:Bm(2), gm:Bm(3))
        tau33(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = tau33(gm:Bm(1), gm:Bm(2), gm:Bm(3)) + lambda(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * dummy3(gm:Bm(1), gm:Bm(2), gm:Bm(3))
        tau33(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = tau33(gm:Bm(1), gm:Bm(2), gm:Bm(3)) + u(     gm:Bm(1), gm:Bm(2), gm:Bm(3)) * tau13( gm:Bm(1), gm:Bm(2), gm:Bm(3))
        tau33(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = tau33(gm:Bm(1), gm:Bm(2), gm:Bm(3)) + v(     gm:Bm(1), gm:Bm(2), gm:Bm(3)) * tau23( gm:Bm(1), gm:Bm(2), gm:Bm(3))

        dummy6(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = tau33( gm:Bm(1), gm:Bm(2), gm:Bm(3))

        !-----------------------------------------------------------------------------------------
        call diff_wrapper_3D( Bs, g, dummy4, 'diff_x_c_gp', onesided, dx=dx, dudx=dummy)
        call diff_wrapper_3D( Bs, g, dummy5, 'diff_y_c_gp', onesided, dy=dy, dudy=dummy2)
        call diff_wrapper_3D( Bs, g, dummy6, 'diff_z_c_gp', onesided, dz=dz, dudz=dummy3)

        rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), EF) = rhs(   gp:Bp(1), gp:Bp(2), gp:Bp(3), EF) &
                                              + dummy( gp:Bp(1), gp:Bp(2), gp:Bp(3)) &
                                              + dummy2(gp:Bp(1), gp:Bp(2), gp:Bp(3)) &
                                              + dummy3(gp:Bp(1), gp:Bp(2), gp:Bp(3)) 

    end if

    ! EQUATIONS
    ! --------------------------------------------------------------------------------------------------------------
    ! RHS of equation of mass: J*srho*2 * srho_t = -div(rho*U_tilde)
    ! --------------------------------------------------------------------------------------------------------------
    ! use tau11, tau22, tau33 as dummy fields
    tau11(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = rho(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * u(gm:Bm(1), gm:Bm(2), gm:Bm(3))
    tau22(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = rho(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * v(gm:Bm(1), gm:Bm(2), gm:Bm(3))
    tau33(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = rho(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * w(gm:Bm(1), gm:Bm(2), gm:Bm(3))

    call diff_wrapper_3D( Bs, g, tau11, 'diff_x_c_gp', onesided, dx=dx, dudx=dummy4)
    call diff_wrapper_3D( Bs, g, tau22, 'diff_y_c_gp', onesided, dy=dy, dudy=dummy5)
    call diff_wrapper_3D( Bs, g, tau33, 'diff_z_c_gp', onesided, dz=dz, dudz=dummy6)

    rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), rhoF) = dummy4( gp:Bp(1), gp:Bp(2), gp:Bp(3)) &
                                            + dummy5( gp:Bp(1), gp:Bp(2), gp:Bp(3)) &
                                            + dummy6( gp:Bp(1), gp:Bp(2), gp:Bp(3))

    rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), rhoF) = - rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), rhoF) * 0.5_rk * phi1_inv(gp:Bp(1), gp:Bp(2), gp:Bp(3))

    ! RHS of  momentum equation for u: sru_t = -1/2 * div(rho U_tilde u ) - 1/2 * (rho*U_tilde)*Du - Dp
    ! --------------------------------------------------------------------------------------------------------------
    dummy( gm:Bm(1), gm:Bm(2), gm:Bm(3)) = u(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * tau11(gm:Bm(1), gm:Bm(2), gm:Bm(3))
    dummy2(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = v(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * tau11(gm:Bm(1), gm:Bm(2), gm:Bm(3))
    dummy3(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = w(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * tau11(gm:Bm(1), gm:Bm(2), gm:Bm(3))

    call diff_wrapper_3D( Bs, g, dummy,  'diff_x_c_gp', onesided, dx=dx, dudx=dummy4)
    call diff_wrapper_3D( Bs, g, dummy2, 'diff_y_c_gp', onesided, dy=dy, dudy=dummy5)
    call diff_wrapper_3D( Bs, g, dummy3, 'diff_z_c_gp', onesided, dz=dz, dudz=dummy6)

    dummy4(gp:Bp(1), gp:Bp(2), gp:Bp(3)) = dummy4( gp:Bp(1), gp:Bp(2), gp:Bp(3)) &
                                         + dummy5( gp:Bp(1), gp:Bp(2), gp:Bp(3)) &
                                         + dummy6( gp:Bp(1), gp:Bp(2), gp:Bp(3))

    rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UxF) = rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UxF) &
                                           - dummy4(gp:Bp(1), gp:Bp(2), gp:Bp(3)) * 0.5_rk

    rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UxF) = rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UxF) &
                                           * phi1_inv(gp:Bp(1), gp:Bp(2), gp:Bp(3))

    ! RHS of  momentum equation for v
    ! --------------------------------------------------------------------------------------------------------------
    dummy( gm:Bm(1), gm:Bm(2), gm:Bm(3)) = u(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * tau22(gm:Bm(1), gm:Bm(2), gm:Bm(3))
    dummy2(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = v(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * tau22(gm:Bm(1), gm:Bm(2), gm:Bm(3))
    dummy3(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = w(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * tau22(gm:Bm(1), gm:Bm(2), gm:Bm(3))

    call diff_wrapper_3D( Bs, g, dummy,  'diff_x_c_gp', onesided, dx=dx, dudx=dummy4)
    call diff_wrapper_3D( Bs, g, dummy2, 'diff_y_c_gp', onesided, dy=dy, dudy=dummy5)
    call diff_wrapper_3D( Bs, g, dummy3, 'diff_z_c_gp', onesided, dz=dz, dudz=dummy6)

    dummy4(gp:Bp(1), gp:Bp(2), gp:Bp(3)) = dummy4( gp:Bp(1), gp:Bp(2), gp:Bp(3)) &
                                         + dummy5( gp:Bp(1), gp:Bp(2), gp:Bp(3)) &
                                         + dummy6( gp:Bp(1), gp:Bp(2), gp:Bp(3))

    rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UyF) = rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UyF) &
                                           - dummy4(gp:Bp(1), gp:Bp(2), gp:Bp(3)) * 0.5_rk

    rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UyF) = rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UyF) &
                                           * phi1_inv(gp:Bp(1), gp:Bp(2), gp:Bp(3))

    ! RHS of  momentum equation for w
    ! --------------------------------------------------------------------------------------------------------------
    dummy( gm:Bm(1), gm:Bm(2), gm:Bm(3)) = u(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * tau33(gm:Bm(1), gm:Bm(2), gm:Bm(3))
    dummy2(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = v(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * tau33(gm:Bm(1), gm:Bm(2), gm:Bm(3))
    dummy3(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = w(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * tau33(gm:Bm(1), gm:Bm(2), gm:Bm(3))

    call diff_wrapper_3D( Bs, g, dummy,  'diff_x_c_gp', onesided, dx=dx, dudx=dummy4)
    call diff_wrapper_3D( Bs, g, dummy2, 'diff_y_c_gp', onesided, dy=dy, dudy=dummy5)
    call diff_wrapper_3D( Bs, g, dummy3, 'diff_z_c_gp', onesided, dz=dz, dudz=dummy6)

    dummy4(gp:Bp(1), gp:Bp(2), gp:Bp(3)) = dummy4( gp:Bp(1), gp:Bp(2), gp:Bp(3)) &
                                         + dummy5( gp:Bp(1), gp:Bp(2), gp:Bp(3)) &
                                         + dummy6( gp:Bp(1), gp:Bp(2), gp:Bp(3))

    rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UzF) = rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UzF) &
                                           - dummy4(gp:Bp(1), gp:Bp(2), gp:Bp(3)) * 0.5_rk

    rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UzF) = rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), UzF) &
                                           * phi1_inv(gp:Bp(1), gp:Bp(2), gp:Bp(3))
