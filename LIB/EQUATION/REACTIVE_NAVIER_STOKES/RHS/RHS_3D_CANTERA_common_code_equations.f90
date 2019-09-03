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
    do k = g-1, Bs(3)+g+2
        do j = g-1, Bs(2)+g+2
            do i = g-1, Bs(1)+g+2
                dummy(i,j,k)  = u(i,j,k) * phi(i,j,k,EF) + u(i,j,k) * p(i,j,k)
                dummy2(i,j,k) = v(i,j,k) * phi(i,j,k,EF) + v(i,j,k) * p(i,j,k)
                dummy3(i,j,k) = w(i,j,k) * phi(i,j,k,EF) + w(i,j,k) * p(i,j,k)
            end do
        end do
    end do

    call diff_wrapper_3D( Bs, g, dummy,  'periodic_u_x', onesided, dx=dx, dudx=dummy4)
    call diff_wrapper_3D( Bs, g, dummy2, 'periodic_u_y', onesided, dy=dy, dudy=dummy5)
    call diff_wrapper_3D( Bs, g, dummy3, 'periodic_u_z', onesided, dz=dz, dudz=dummy6)

    do k = g+1, Bs(3)+g
        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                rhs(i,j,k,EF) = rhs(i,j,k,EF) - ( dummy4(i,j,k) + dummy5(i,j,k) + dummy6(i,j,k) )
            end do
        end do 
    end do    

    !#########################################################################################
    ! species equation 
    ! diffusion
    ! (rho * D_k * W_k / W * (X_k)_xi)_xi
    !---------------------------------------------------------------------------------------------
    do n = 1, params_physics%species-1

        call diff_wrapper_3D( Bs, g, X(:,:,:,n), 'periodic_u_xyz', onesided, dx=dx, dy=dy, dz=dz, dudx=dummy, dudy=dummy2, dudz=dummy3)

        do k = g-1, Bs(3)+g+2
            do j = g-1, Bs(2)+g+2
                do i = g-1, Bs(1)+g+2

                    ! diffusion
                    dummy(i,j,k)  = rho(i,j,k) * D(i,j,k,n) * dummy(i,j,k)
                    dummy2(i,j,k) = rho(i,j,k) * D(i,j,k,n) * dummy2(i,j,k)
                    dummy3(i,j,k) = rho(i,j,k) * D(i,j,k,n) * dummy3(i,j,k)

                end do
            end do
        end do

        call diff_wrapper_3D( Bs, g, dummy,  'periodic_u_x', onesided, dx=dx, dudx=dummy4)
        call diff_wrapper_3D( Bs, g, dummy2, 'periodic_u_y', onesided, dy=dy, dudy=dummy5)
        call diff_wrapper_3D( Bs, g, dummy3, 'periodic_u_z', onesided, dz=dz, dudz=dummy6)

        do k = g+1, Bs(3)+g
            do j = g+1, Bs(2)+g
                do i = g+1, Bs(1)+g
                    rhs(i,j,k,YF+n-1) = rhs(i,j,k,YF+n-1) + ( dummy4(i,j,k) + dummy5(i,j,k) + dummy6(i,j,k) )
                end do
            end do
        end do

    end do

    ! transport
    ! - (rho * (u_i + Vc_i) * Y_k)_xi
    ! note: YF contains rho*Y_k
    ! -------------------------
    do n = 1, params_physics%species-1

        do k = g-1, Bs(3)+g+2
            do j = g-1, Bs(2)+g+2
                do i = g-1, Bs(1)+g+2
                    dummy(i,j,k)  = ( u(i,j,k) ) * phi(i,j,k,YF+n-1)
                    dummy2(i,j,k) = ( v(i,j,k) ) * phi(i,j,k,YF+n-1)
                    dummy3(i,j,k) = ( w(i,j,k) ) * phi(i,j,k,YF+n-1)
                end do
            end do
        end do

        call diff_wrapper_3D( Bs, g, dummy,  'periodic_u_x', onesided, dx=dx, dudx=dummy4)
        call diff_wrapper_3D( Bs, g, dummy2, 'periodic_u_y', onesided, dy=dy, dudy=dummy5)
        call diff_wrapper_3D( Bs, g, dummy3, 'periodic_u_z', onesided, dz=dz, dudz=dummy6)

        do k = g+1, Bs(3)+g
            do j = g+1, Bs(2)+g
                do i = g+1, Bs(1)+g
                    rhs(i,j,k,YF+n-1) = rhs(i,j,k,YF+n-1) - ( dummy4(i,j,k) + dummy5(i,j,k) + dummy6(i,j,k) )
                end do
            end do
        end do

    end do
