!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name RHS_3D_CANTERA_common_code.f90
!> \version 0.5
!> \author msr
!
!> \brief source code for 3D CANTERA ns, prevent double code
!
! ********************************************************************************************

!---------------------------------------------------------------------------------------------
! main body

    !#########################################################################################
    ! Compute Rs, W, X, D
    if (dissipation) then

        ! use X as dummy
        do n = 1, params_physics%species
            X(:,:,:,n) = Wk_species_inv(n) * Y(:,:,:,n)
        end do

        ! compute W
        dummy(:,:,:) = X(:,:,:,1)
        do n = 2, params_physics%species
            dummy(:,:,:) = dummy(:,:,:) + X(:,:,:,n)
        end do

        ! compute X
        do n = 1, params_physics%species
            X(:,:,:,n) = X(:,:,:,n) / dummy(:,:,:)
        end do

        ! compute D
        do n = 1, params_physics%species
            D(:,:,:,n) = params_physics%Wk(n) * dummy(:,:,:)
        end do

    end if    

    !#########################################################################################
    ! dummy energy for set state calls
    dummy(:,:,:) = es(:,:,:)
    do n = 1, params_physics%species
        dummy(:,:,:) = dummy(:,:,:) + dh(n) * Y(:,:,:, n)
    end do

    !#########################################################################################
    ! compute thermodynamics for outer ghost nodes
    do k = 1, Ds(3)
        do j = 1, Ds(2)

            do index3 = 1, 4
                i = outer_bound(index3, 1)

                    ! setState call depends on choosen mechanism
                    ! OHN
                    call setMoleFractions(gas, X(i,j,k,:))
                    call setState_UV(gas, dummy(i,j,k), phi1_inv(i,j,k)*phi1_inv(i,j,k)  )

                    ! temperature and pressure
                    T(i,j,k) = temperature(gas)
                    p(i,j,k) = pressure(gas)
  
            end do
        end do
    end do

    do k = 1, Ds(3)

        do index2 = 1, 4
            j = outer_bound(index2, 2)

            do i = gm, Bm(1)

                    ! setState call depends on choosen mechanism
                    ! OHN
                    call setMoleFractions(gas, X(i,j,k,:))
                    call setState_UV(gas, dummy(i,j,k), phi1_inv(i,j,k)*phi1_inv(i,j,k)  )

                    ! temperature and pressure
                    T(i,j,k) = temperature(gas)
                    p(i,j,k) = pressure(gas)
  
            end do
        end do
    end do

    do index1 = 1, 4
        k = outer_bound(index1, 3)

        do j = gm, Bm(2)
            do i = gm, Bm(1)

                    ! setState call depends on choosen mechanism
                    ! OHN
                    call setMoleFractions(gas, X(i,j,k,:))
                    call setState_UV(gas, dummy(i,j,k), phi1_inv(i,j,k)*phi1_inv(i,j,k)  )

                    ! temperature and pressure
                    T(i,j,k) = temperature(gas)
                    p(i,j,k) = pressure(gas)
  
            end do
        end do
    end do

    !#########################################################################################
    ! compute thermodynamics for inner ghost nodes
    do k = gm, Bm(3)
        do j = gm, Bm(2)

            do index3 = 1, 4
                i = inner_bound(index3, 1)

                    ! setState call depends on choosen mechanism
                    ! OHN
                    call setMoleFractions(gas, X(i,j,k,:))
                    call setState_UV(gas, dummy(i,j,k), phi1_inv(i,j,k)*phi1_inv(i,j,k)  )

                    ! temperature and pressure
                    T(i,j,k) = temperature(gas)
                    p(i,j,k) = pressure(gas)

                    ! viscosity
                    mu(i,j,k) = viscosity(gas)

                    ! thermal conductivity
                    lambda(i,j,k) = thermalConductivity(gas)    
                
                    ! diffusion
                    call getMixDiffCoeffs(gas, q)

                    do n = 1, params_physics%species
                        D(i,j,k,n) = D(i,j,k,n)  * q(n) 
                    end do
  
            end do
        end do
    end do

    do k = gm, Bm(3)

        do index2 = 1, 4
            j = inner_bound(index2, 2)

            do i = gp, Bp(1)

                    ! setState call depends on choosen mechanism
                    ! OHN
                    call setMoleFractions(gas, X(i,j,k,:))
                    call setState_UV(gas, dummy(i,j,k), phi1_inv(i,j,k)*phi1_inv(i,j,k)  )

                    ! temperature and pressure
                    T(i,j,k) = temperature(gas)
                    p(i,j,k) = pressure(gas)

                    ! viscosity
                    mu(i,j,k) = viscosity(gas)

                    ! thermal conductivity
                    lambda(i,j,k) = thermalConductivity(gas)    
                
                    ! diffusion
                    call getMixDiffCoeffs(gas, q)

                    do n = 1, params_physics%species
                        D(i,j,k,n) = D(i,j,k,n)  * q(n) 
                    end do
  
            end do
        end do
    end do

    do index1 = 1, 4
        k = inner_bound(index1, 3)

        do j = gp, Bp(2)
            do i = gp, Bp(1)

                    ! setState call depends on choosen mechanism
                    ! OHN
                    call setMoleFractions(gas, X(i,j,k,:))
                    call setState_UV(gas, dummy(i,j,k), phi1_inv(i,j,k)*phi1_inv(i,j,k)  )

                    ! temperature and pressure
                    T(i,j,k) = temperature(gas)
                    p(i,j,k) = pressure(gas)

                    ! viscosity
                    mu(i,j,k) = viscosity(gas)

                    ! thermal conductivity
                    lambda(i,j,k) = thermalConductivity(gas)    
                
                    ! diffusion
                    call getMixDiffCoeffs(gas, q)

                    do n = 1, params_physics%species
                        D(i,j,k,n) = D(i,j,k,n)  * q(n) 
                    end do
  
            end do
        end do
    end do

    !#########################################################################################
    ! compute thermodynamics for inner block nodes
    do k = gp, Bp(3)
        do j = gp, Bp(2)
            do i = gp, Bp(1)

                    ! setState call depends on choosen mechanism
                    ! OHN
                    call setMoleFractions(gas, X(i,j,k,:))
                    call setState_UV(gas, dummy(i,j,k), phi1_inv(i,j,k)*phi1_inv(i,j,k)  )

                    ! temperature and pressure
                    T(i,j,k) = temperature(gas)
                    p(i,j,k) = pressure(gas)

                    ! viscosity
                    mu(i,j,k) = viscosity(gas)

                    ! thermal conductivity
                    lambda(i,j,k) = thermalConductivity(gas)

                    ! diffusion
                    call getMixDiffCoeffs(gas, q)

                    do n = 1, params_physics%species
                        D(i,j,k,n) = D(i,j,k,n)  * q(n) 
                    end do

                    ! net production srates
                    call getNetProductionRates(gas, q)

                    ! correct reaction rate
                    do n = 1, params_physics%species
                        rhs(i,j,k,YF+n-1) = q(n) * params_physics%Wk(n)
                    end do
    
            end do
        end do
    end do


    !#########################################################################################
    ! insert reaction rate in energy equation
    rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), EF) = - dh(1) * rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), YF)
    do n = 2, params_physics%species
       rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), EF) = rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), EF) &
                                             - dh(n) * rhs(gp:Bp(1), gp:Bp(2), gp:Bp(3), YF+n-1)
    end do

    !#########################################################################################
    ! viscosity dummy field
    ! mu_d - 2/3 mu, mu_d = 0
    dummy6(gm:Bm(1), gm:Bm(2), gm:Bm(3)) = - two_three * mu(gm:Bm(1), gm:Bm(2), gm:Bm(3))
