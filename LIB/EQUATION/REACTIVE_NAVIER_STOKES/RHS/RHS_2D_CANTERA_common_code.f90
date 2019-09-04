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

        do j = 1, Bs(2)+2*g
            do i = 1, Bs(1)+2*g

                ! mean molecular weight
                Wk = 0.0_rk
                do n = 1, params_physics%species
                    Wk = Wk + Wk_species_inv(n) * Y(i,j,n)
                end do
                ! note: W is actually 1/W !!
                Wk_inv = 1.0_rk / Wk
                ! gas constant
                Rs = Wk * params_physics%gas_constant

                ! mole fraction
                do n = 1, params_physics%species
                    X(i,j,n)         = Wk_inv * Wk_species_inv(n) * Y(i,j,n)
                end do

                ! diffusion coefficient
                do n = 1, params_physics%species
                    D(i,j,n) = params_physics%Wk(n) * Wk
                end do

                end do
            end do

        mu_d = 0.0_rk

    end if    

    !#########################################################################################
    ! CANTERA calls
    ! - principle: set thermodynamic state for each node, compute all necessary values with cantera
    ! note: if-condition inside loop, should ensure the minimal number of cantera calls 
    do j = 1, Bs(2)+g+g
        do i = 1, Bs(1)+g+g

            if ( (j<g-1) .or. (i<g-1) .or. (j>Bs(2)+g+2) .or. (i>Bs(1)+g+2) ) then

                ! setState call depends on choosen mechanism
                ! OHN
                call setMoleFractions(gas, X(i,j,:))
                call setState_UV(gas, es(i,j) + sum( dh(:) * Y(i,j,:) ), phi1_inv(i,j)*phi1_inv(i,j)  )

                ! temperature and pressure
                T(i,j) = temperature(gas)
                p(i,j) = pressure(gas)

            elseif ( (j<g+1) .or. (i<g+1) .or. (j>Bs(2)+g) .or. (i>Bs(1)+g) ) then

                ! setState call depends on choosen mechanism
                ! OHN
                call setMoleFractions(gas, X(i,j,:))
                call setState_UV(gas, es(i,j) + sum( dh(:) * Y(i,j,:) ), phi1_inv(i,j)*phi1_inv(i,j)  )

                ! temperature and pressure
                T(i,j) = temperature(gas)
                p(i,j) = pressure(gas)

                ! viscosity
                mu(i,j) = viscosity(gas)

                ! thermal conductivity
                lambda(i,j) = thermalConductivity(gas)    
                
                ! diffusion
                call getMixDiffCoeffs(gas, q)

                do n = 1, params_physics%species
                    D(i,j,n) = D(i,j,n)  * q(n) 
                end do

            else

                ! setState call depends on choosen mechanism
                ! OHN
                call setMoleFractions(gas, X(i,j,:))
                call setState_UV(gas, es(i,j) + sum( dh(:) * Y(i,j,:) ), phi1_inv(i,j)*phi1_inv(i,j)  )

                ! temperature and pressure
                T(i,j) = temperature(gas)
                p(i,j) = pressure(gas)

                ! viscosity
                mu(i,j) = viscosity(gas)

                ! thermal conductivity
                lambda(i,j) = thermalConductivity(gas)

                ! diffusion
                call getMixDiffCoeffs(gas, q)

                do n = 1, params_physics%species
                    D(i,j,n) = D(i,j,n)  * q(n) 
                end do

                ! net production srates
                call getNetProductionRates(gas, q)

                ! correct reaction rate
                do n = 1, params_physics%species
                    rhs(i,j,1,YF+n-1) = q(n) * params_physics%Wk(n)
                end do

                rhs(i,j,1,EF) = - sum( dh(:) * rhs(i,j,1,YF:YF+params_physics%species-1) )

            end if
    
        end do
    end do
