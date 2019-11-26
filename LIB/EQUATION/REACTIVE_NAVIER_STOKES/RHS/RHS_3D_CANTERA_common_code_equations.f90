!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name RHS_3D_CANTERA_common_code_equations.f90
!> \version 0.5
!> \author msr
!
!> \brief source code for 3D CANTERA, prevent double code
!> energy and species equation
!
! ********************************************************************************************

!---------------------------------------------------------------------------------------------
! main body

    !#########################################################################################
    ! energy equation
    ! RHS of energy equation:  p_t = -gamma*div(U_tilde p) + gamm1 *U x grad(p)
    dummy( gm:Bm(1), gm:Bm(2), gm:Bm(3)) = u(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * p(gm:Bm(1), gm:Bm(2), gm:Bm(3))
    dummy2(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = v(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * p(gm:Bm(1), gm:Bm(2), gm:Bm(3))
    dummy3(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = w(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * p(gm:Bm(1), gm:Bm(2), gm:Bm(3))

    dummy( gm:Bm(1), gm:Bm(2), gm:Bm(3)) = dummy( gm:Bm(1), gm:Bm(2), gm:Bm(3)) + u(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * phi(gm:Bm(1), gm:Bm(2), gm:Bm(3), EF)
    dummy2(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = dummy2(gm:Bm(1), gm:Bm(2), gm:Bm(3)) + v(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * phi(gm:Bm(1), gm:Bm(2), gm:Bm(3), EF)
    dummy3(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = dummy3(gm:Bm(1), gm:Bm(2), gm:Bm(3)) + w(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * phi(gm:Bm(1), gm:Bm(2), gm:Bm(3), EF)


    call diff_wrapper_3D( Bs, g, dummy,  'diff_x_c_gp', onesided, dx=dx, dudx=dummy4)
    call diff_wrapper_3D( Bs, g, dummy2, 'diff_y_c_gp', onesided, dy=dy, dudy=dummy5)
    call diff_wrapper_3D( Bs, g, dummy3, 'diff_z_c_gp', onesided, dz=dz, dudz=dummy6)

    dummy4(gp:Bp(1), gp:Bp(2), gp:Bp(3)) = dummy4( gp:Bp(1), gp:Bp(2), gp:Bp(3)) &
                                         + dummy5( gp:Bp(1), gp:Bp(2), gp:Bp(3)) &
                                         + dummy6( gp:Bp(1), gp:Bp(2), gp:Bp(3))

    rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), EF) = rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), EF) &
                                          - dummy4(gp:Bp(1), gp:Bp(2), gp:Bp(3))  

    !#########################################################################################
    ! species equation 
    ! diffusion
    ! (rho * D_k * W_k / W * (X_k)_xi)_xi
    !---------------------------------------------------------------------------------------------
    do n = 1, params_physics%species-1

        call diff_wrapper_3D( Bs, g, X(:,:,:,n), 'diff_x_c_gm', onesided, dx=dx, dudx=dummy )
        call diff_wrapper_3D( Bs, g, X(:,:,:,n), 'diff_y_c_gm', onesided, dy=dy, dudy=dummy2)
        call diff_wrapper_3D( Bs, g, X(:,:,:,n), 'diff_z_c_gm', onesided, dz=dz, dudz=dummy3)

        dummy4( gm:Bm(1), gm:Bm(2), gm:Bm(3)) = rho(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * D(gm:Bm(1), gm:Bm(2), gm:Bm(3), n)

        dummy( gm:Bm(1), gm:Bm(2), gm:Bm(3)) = dummy( gm:Bm(1), gm:Bm(2), gm:Bm(3)) * dummy4(gm:Bm(1), gm:Bm(2), gm:Bm(3))
        dummy2(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = dummy2(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * dummy4(gm:Bm(1), gm:Bm(2), gm:Bm(3))
        dummy3(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = dummy3(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * dummy4(gm:Bm(1), gm:Bm(2), gm:Bm(3))

        call diff_wrapper_3D( Bs, g, dummy,  'diff_x_c_gp', onesided, dx=dx, dudx=dummy4)
        call diff_wrapper_3D( Bs, g, dummy2, 'diff_y_c_gp', onesided, dy=dy, dudy=dummy5)
        call diff_wrapper_3D( Bs, g, dummy3, 'diff_z_c_gp', onesided, dz=dz, dudz=dummy6)

        dummy4(gp:Bp(1), gp:Bp(2), gp:Bp(3)) = dummy4( gp:Bp(1), gp:Bp(2), gp:Bp(3)) &
                                             + dummy5( gp:Bp(1), gp:Bp(2), gp:Bp(3)) &
                                             + dummy6( gp:Bp(1), gp:Bp(2), gp:Bp(3))

        rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), YF+n-1) = rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), YF+n-1) &
                                                  + dummy4(gp:Bp(1), gp:Bp(2), gp:Bp(3)) 

    end do

    ! transport
    ! - (rho * (u_i + Vc_i) * Y_k)_xi
    ! note: YF contains rho*Y_k
    ! -------------------------
    do n = 1, params_physics%species-1

        dummy( gm:Bm(1), gm:Bm(2), gm:Bm(3)) = u(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * phi(gm:Bm(1), gm:Bm(2), gm:Bm(3), YF+n-1)
        dummy2(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = v(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * phi(gm:Bm(1), gm:Bm(2), gm:Bm(3), YF+n-1)
        dummy3(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = w(gm:Bm(1), gm:Bm(2), gm:Bm(3)) * phi(gm:Bm(1), gm:Bm(2), gm:Bm(3), YF+n-1)       

        call diff_wrapper_3D( Bs, g, dummy,  'diff_x_c_gp', onesided, dx=dx, dudx=dummy4)
        call diff_wrapper_3D( Bs, g, dummy2, 'diff_y_c_gp', onesided, dy=dy, dudy=dummy5)
        call diff_wrapper_3D( Bs, g, dummy3, 'diff_z_c_gp', onesided, dz=dz, dudz=dummy6)

        dummy4(gp:Bp(1), gp:Bp(2), gp:Bp(3)) = dummy4( gp:Bp(1), gp:Bp(2), gp:Bp(3)) &
                                             + dummy5( gp:Bp(1), gp:Bp(2), gp:Bp(3)) &
                                             + dummy6( gp:Bp(1), gp:Bp(2), gp:Bp(3))

        rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), YF+n-1) = rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), YF+n-1) &
                                                  - dummy4(gp:Bp(1), gp:Bp(2), gp:Bp(3)) 

    end do
