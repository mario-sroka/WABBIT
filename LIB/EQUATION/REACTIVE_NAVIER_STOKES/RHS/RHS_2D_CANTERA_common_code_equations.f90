!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name RHS_2D_CANTERA_common_code_equations.f90
!> \version 0.5
!> \author msr
!
!> \brief source code for 2D CANTERA, prevent double code
!> energy and species equation
!
! ********************************************************************************************

!---------------------------------------------------------------------------------------------
! main body

    !#########################################################################################
    ! energy equation
    ! RHS of energy equation:  p_t = -gamma*div(U_tilde p) + gamm1 *U x grad(p)
    do j = g-1, Bs(2)+g+2
        do i = g-1, Bs(1)+g+2
            dummy(i,j)  = u(i,j) * phi(i,j,1,EF) + u(i,j) * p(i,j)
            dummy2(i,j) = v(i,j) * phi(i,j,1,EF) + v(i,j) * p(i,j)
        end do
    end do

    call diff_wrapper_2D( Bs, g, dummy,  'periodic_u_x', onesided, dx=dx, dudx=dummy3)
    call diff_wrapper_2D( Bs, g, dummy2, 'periodic_u_y', onesided, dy=dy, dudy=dummy4)

    do j = g+1, Bs(2)+g
        do i = g+1, Bs(1)+g
            rhs(i,j,1,EF) = rhs(i,j,1,EF) - ( dummy3(i,j) + dummy4(i,j) )
        end do
    end do 
 

    !#########################################################################################
    ! species equation 
    ! diffusion
    ! (rho * D_k * W_k / W * (X_k)_xi)_xi
    !---------------------------------------------------------------------------------------------
    do n = 1, params_physics%species-1

        call diff_wrapper_2D( Bs, g, X(:,:,n), 'periodic_u_xy', onesided, dx=dx, dy=dy, dudx=dummy, dudy=dummy2)

        do j = g-1, Bs(2)+g+2
            do i = g-1, Bs(1)+g+2

                ! diffusion
                dummy(i,j)  = rho(i,j) * D(i,j,n) * dummy(i,j)
                dummy2(i,j) = rho(i,j) * D(i,j,n) * dummy2(i,j)

            end do
        end do

        call diff_wrapper_2D( Bs, g, dummy,  'periodic_u_x', onesided, dx=dx, dudx=dummy3)
        call diff_wrapper_2D( Bs, g, dummy2, 'periodic_u_y', onesided, dy=dy, dudy=dummy4)

        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                rhs(i,j,1,YF+n-1) = rhs(i,j,1,YF+n-1) + ( dummy3(i,j) + dummy4(i,j) )
            end do
        end do

    end do

    ! transport
    ! - (rho * (u_i + Vc_i) * Y_k)_xi
    ! note: YF contains rho*Y_k
    ! -------------------------
    do n = 1, params_physics%species-1

        do j = g-1, Bs(2)+g+2
            do i = g-1, Bs(1)+g+2
                dummy(i,j)  = ( u(i,j) ) * phi(i,j,1,YF+n-1)
                dummy2(i,j) = ( v(i,j) ) * phi(i,j,1,YF+n-1)
            end do
        end do

        call diff_wrapper_2D( Bs, g, dummy,  'periodic_u_x', onesided, dx=dx, dudx=dummy3)
        call diff_wrapper_2D( Bs, g, dummy2, 'periodic_u_y', onesided, dy=dy, dudy=dummy4)

        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g
                rhs(i,j,1,YF+n-1) = rhs(i,j,1,YF+n-1) - ( dummy3(i,j) + dummy4(i,j) )
            end do
        end do

    end do
