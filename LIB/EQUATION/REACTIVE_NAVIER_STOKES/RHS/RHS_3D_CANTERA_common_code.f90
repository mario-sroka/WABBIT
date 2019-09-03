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

        do k = 1, Bs(3)+2*g
            do j = 1, Bs(2)+2*g
                do i = 1, Bs(1)+2*g

                    ! mean molecular weight
                    Wk = 0.0_rk
                    do n = 1, params_physics%species
                        Wk = Wk + Wk_species_inv(n) * Y(i,j,k,n)
                    end do
                    ! note: W is actually 1/W !!
                    Wk_inv = 1.0_rk / Wk
                    ! gas constant
                    Rs = Wk * params_physics%gas_constant

                    ! mole fraction
                    do n = 1, params_physics%species
                        X(i,j,k,n)         = Wk_inv * Wk_species_inv(n) * Y(i,j,k,n)
                    end do

                    ! diffusion coefficient
                    do n = 1, params_physics%species
                        D(i,j,k,n) = params_physics%Wk(n) * Wk
                    end do

                end do
            end do
        end do

        mu_d = 0.0_rk

    end if    

    !#########################################################################################
    ! CANTERA calls
    ! - principle: set thermodynamic state for each node, compute all necessary values with cantera
    ! note: if-condition inside loop, should ensure the minimal number of cantera calls 
    do k = 1, Bs(3)+2*g
        do j = 1, Bs(2)+g+g
            do i = 1, Bs(1)+g+g

                if ( (k<g-1) .or. (j<g-1) .or. (i<g-1) .or. (k>Bs(3)+g+2) .or. (j>Bs(2)+g+2) .or. (i>Bs(1)+g+2) ) then

                    ! setState call depends on choosen mechanism
                    ! OHN
                    call setMoleFractions(gas, X(i,j,k,:))
                    call setState_UV(gas, es(i,j,k) + sum( dh(:) * Y(i,j,k,:) ), phi1_inv(i,j,k)*phi1_inv(i,j,k)  )

                    ! temperature and pressure
                    T(i,j,k) = temperature(gas)
                    p(i,j,k) = pressure(gas)

                elseif ( (k<g+1) .or. (j<g+1) .or. (i<g+1) .or. (k>Bs(3)+g) .or. (j>Bs(2)+g) .or. (i>Bs(1)+g) ) then

                    ! setState call depends on choosen mechanism
                    ! OHN
                    call setMoleFractions(gas, X(i,j,k,:))
                    call setState_UV(gas, es(i,j,k) + sum( dh(:) * Y(i,j,k,:) ), phi1_inv(i,j,k)*phi1_inv(i,j,k)  )

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

                else

                    ! setState call depends on choosen mechanism
                    ! OHN
                    call setMoleFractions(gas, X(i,j,k,:))
                    call setState_UV(gas, es(i,j,k) + sum( dh(:) * Y(i,j,k,:) ), phi1_inv(i,j,k)*phi1_inv(i,j,k)  )

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

                    rhs(i,j,k,EF) = - sum( dh(:) * rhs(i,j,k,YF:YF+params_physics%species-1) )

                end if
    
            end do
        end do
    end do
