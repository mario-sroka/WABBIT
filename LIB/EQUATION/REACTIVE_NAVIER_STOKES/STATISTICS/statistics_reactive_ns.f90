!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name statistics_reactive_ns.f90
!> \version 0.5
!> \author msr
!
!> \brief compute statistics (e.g. total energy, ...)
!>
!! input:    - params_physics, data for one block, spacing parameter \n
!! output:   - \n
!!
!!
!! = log ======================================================================================
!! \n
!! 19/10/18 - create
!
! ********************************************************************************************

subroutine statistics_reactive_ns( params_physics, time, phi, phi_work, x0, dx, stage, lgt_n, gas )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> parameter struct
    type(type_params_rns), intent(inout)    :: params_physics

    ! it may happen that some source terms have an explicit time-dependency
    ! therefore the general call has to pass time
    real(kind=rk), intent (in)              :: time
    !> state vector for one block
    real(kind=rk), intent(in)               :: phi(:, :, :, :)
    real(kind=rk), intent(inout)            :: phi_work(:, :, :, :)
    !> spacing and origin of block
    real(kind=rk), intent(in)               :: x0(1:3), dx(1:3)

    ! stage. there is 3 stages, init_stage, integral_stage and local_stage. If the PDE has
    ! terms that depend on global qtys, such as forces etc, which cannot be computed
    ! from a single block alone, the first stage does that. the second stage can then
    ! use these integral qtys for the actual RHS evaluation.
    character(len=*), intent(in)            :: stage

    !> number of active blocks (light data)
    integer(kind=ik), intent(inout)         :: lgt_n

    !> Cantera gas mixture struct
    type(phase_t), intent(inout)            :: gas

    ! grid parameter
    integer(kind=ik)                        :: Bs(3), g, gridNodes
    ! loop parameter
    integer(kind=ik)                        :: i, j, k, n
    ! field indexes
    integer(kind=ik)                        :: rhoF, UxF, UyF, UzF, EF, YF
    ! dummy variables
    real(kind=rk)                           :: dummy(4), Y(params_physics%species), X(params_physics%species)
    ! temp array
    real(kind=rk), allocatable, save        :: tmp(:,:,:,:)
    ! MPI error variable
    integer(kind=ik)                        :: ierr

    ! forcing:
    ! --------
    ! wavenumber limit, wavenumber
    real(kind=rk)                           :: k0, wn(3)
    ! DFT loop limit
    integer(kind=ik)                        :: kmax
    ! fourier coeficients
    complex(kind=rk), allocatable, save     :: phi_hat(:,:,:,:)
    ! dummy variable
    complex(kind=rk)                        :: cdummy, wn2(3)

!---------------------------------------------------------------------------------------------
! variables initialization

    ! grid parameter
    Bs = params_physics%Bs
    g  = params_physics%g

    ! field indexes from params
    rhoF = params_physics%rhoF
    UxF  = params_physics%UxF
    UyF  = params_physics%UyF
    UzF  = params_physics%UzF
    EF   = params_physics%EF
    YF   = params_physics%YF

    ! temp array
    if ( .not. allocated(tmp) ) allocate(tmp(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g, 6))
    
    ! forcing parameter
    k0   = params_physics%k0
    kmax = params_physics%kmax
    if ( .not. allocated(phi_hat) ) allocate(phi_hat(kmax, kmax, kmax, 3))

!---------------------------------------------------------------------------------------------
! main body

    select case(stage)
        case ("init_stage")
        !-------------------------------------------------------------------------
        ! 1st stage: init_stage.
        !-------------------------------------------------------------------------
        ! this stage is called only once, not for each block.
        ! performs initializations, such as resetting integrals

        ! energy mean values
        params_physics%e_tot_mean       = 0.0_rk
        params_physics%e_kin_mean       = 0.0_rk
        params_physics%e_int_mean       = 0.0_rk

        ! turbulence properties mean
        params_physics%UU_mean          = 0.0_rk

        params_physics%k_mean           = 0.0_rk
        params_physics%k_d_mean         = 0.0_rk
        params_physics%k_s_mean         = 0.0_rk
               
        params_physics%epsilon_mean     = 0.0_rk
        params_physics%epsilon_d_mean   = 0.0_rk
        params_physics%epsilon_s_mean   = 0.0_rk
        
        ! fluid properties mean
        params_physics%rho_mean         = 0.0_rk
        params_physics%u_mean           = 0.0_rk
        params_physics%v_mean           = 0.0_rk
        params_physics%w_mean           = 0.0_rk
        params_physics%mu_mean          = 0.0_rk

        ! fourier coefficients
        params_physics%phi_hat          = 0.0_rk

    case ("integral_stage")
        !-------------------------------------------------------------------------
        ! 2nd stage: integral_stage.
        !-------------------------------------------------------------------------
        ! This stage contains all operations which are running on the blocks
        ! called for each block

        !-------------------------------------------------------------------------
        ! first: compute mean values for dissipation and turbulent kinetic energy
        !-------------------------------------------------------------------------

        !-------------------------------
        ! mean dilatational dissipation
        !-------------------------------
        call compute_dilatational_dissipation(phi(:,:,:,UxF)/phi(:,:,:,rhoF), &
                                              phi(:,:,:,UyF)/phi(:,:,:,rhoF), &
                                              phi(:,:,:,UzF)/phi(:,:,:,rhoF), &
                                              dx, Bs, g, params_physics%d, tmp(:,:,:,1:6))

        if ( params_physics%d == 2 ) then
            ! /todo: 2D statistics

        else
            ! 3D statistics
            do k = g+1, Bs(3)+g-1
                do j = g+1, Bs(2)+g-1
                    do i = g+1, Bs(1)+g-1

                        ! readability
                        dummy(1) = phi(i, j, k, rhoF)**2.0_rk
                       
                        ! note: favre mean values, compute .../rho_mean and mu_mean in post-stage!
                        params_physics%epsilon_d_mean = params_physics%epsilon_d_mean &
                                                      + dummy(1) * sum(tmp(i, j, k, 1:6))


                    end do
                end do
            end do

        end if

        !-------------------------------
        ! mean solenodial dissipation
        !-------------------------------
        call compute_solenoidal_dissipation(phi(:,:,:,UxF)/phi(:,:,:,rhoF), &
                                            phi(:,:,:,UyF)/phi(:,:,:,rhoF), &
                                            phi(:,:,:,UzF)/phi(:,:,:,rhoF), &
                                            dx, Bs, g, params_physics%d, tmp(:,:,:,1:6))

        if ( params_physics%d == 2 ) then
            ! /todo: 2D statistics

        else
            ! 3D statistics
            do k = g+1, Bs(3)+g-1
                do j = g+1, Bs(2)+g-1
                    do i = g+1, Bs(1)+g-1

                        ! readability
                        dummy(1) = phi(i, j, k, rhoF)**2.0_rk
                       
                        ! note: favre mean values, compute .../rho_mean and mu_mean in post-stage!
                        params_physics%epsilon_s_mean = params_physics%epsilon_s_mean &
                                                      + dummy(1) * sum(tmp(i, j, k, 1:6))

                    end do
                end do
            end do

        end if

        !-------------------------------
        ! mean dilatational turbulent kinetic energy
        !-------------------------------
        ! note: no 2D case
        if (params_physics%d == 3) then

            ! disable solenodidal/dilatational decomposition for cases without forcing
            if (params_physics%forcing) then

                ! compute IDFT, dilatational
                call compute_IDFT( params_physics, params_physics%phi_hat_d(:,:,:,1), tmp(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 1), x0, dx )
                call compute_IDFT( params_physics, params_physics%phi_hat_d(:,:,:,2), tmp(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 2), x0, dx )
                call compute_IDFT( params_physics, params_physics%phi_hat_d(:,:,:,3), tmp(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 3), x0, dx )

                do k = g+1, Bs(3)+g-1
                    do j = g+1, Bs(2)+g-1
                        do i = g+1, Bs(1)+g-1

                            ! readability
                            dummy(1) = phi(i, j, k, rhoF)**2.0_rk
                            dummy(2) = phi(i, j, k, UxF) / phi(i, j, k, rhoF)
                            dummy(3) = phi(i, j, k, UyF) / phi(i, j, k, rhoF)
                            dummy(4) = phi(i, j, k, UzF) / phi(i, j, k, rhoF)
                       
                            ! note: favre mean values, compute .../rho_mean and mu_mean in post-stage!
                            params_physics%k_d_mean = params_physics%k_d_mean + &
                                                     0.5_rk * dummy(1) * ( dummy(2) * tmp(i, j, k, 1) &
                                                            + dummy(3) * tmp(i, j, k, 2) &
                                                            + dummy(4) * tmp(i, j, k, 3) )

                        end do
                    end do
                end do

            end if

        end if

        !-------------------------------------------------------------------------
        ! second: compute all other mean values in one big loop over all block nodes
        ! note: assume mean value variable in params struct have set to zero in init stage 
        ! note: do not double sum over redundant nodes!
        !-------------------------------------------------------------------------

        if ( params_physics%d == 2 ) then
            ! /todo: 2D statistics

        else
            ! 3D statistics
            do k = g+1, Bs(3)+g-1
                do j = g+1, Bs(2)+g-1
                    do i = g+1, Bs(1)+g-1

                        ! avoid double computations
                        dummy(1) = phi(i, j, k, rhoF)**2.0_rk
                        dummy(2) = phi(i, j, k, UxF) / phi(i, j, k, rhoF)
                        dummy(3) = phi(i, j, k, UyF) / phi(i, j, k, rhoF)
                        dummy(4) = phi(i, j, k, UzF) / phi(i, j, k, rhoF)

                        !-----------------------
                        ! mean fluid properties
                        !-----------------------
                        params_physics%rho_mean = params_physics%rho_mean + dummy(1)
                        params_physics%u_mean   = params_physics%u_mean   + dummy(2)
                        params_physics%v_mean   = params_physics%v_mean   + dummy(3)
                        params_physics%w_mean   = params_physics%w_mean   + dummy(4)

                        !-----------------------
                        ! mean energies
                        !-----------------------
                        ! note: total energy assume pressure as energy field
                        ! for other chemistry (e.g. cantera) use e_kin + e_int as total energy in postprocessing!
                        ! /todo: add chemical energy?
                        params_physics%e_kin_mean = params_physics%e_kin_mean + &
                                           0.5_rk * ( dummy(2)**2.0_rk + dummy(3)**2.0_rk + dummy(4)**2.0_rk )

                        params_physics%e_int_mean = params_physics%e_int_mean + phi(i, j, k, EF)


                        params_physics%e_tot_mean = params_physics%e_tot_mean + &
                                           0.5_rk * ( dummy(2)**2.0_rk + dummy(3)**2.0_rk + dummy(4)**2.0_rk ) + &
                                           1.0_rk / ( params_physics%gamma_ - 1.0_rk ) * phi(i, j, k, EF) / dummy(1)
                        
                        !-----------------------
                        ! turbulence properties
                        !-----------------------
                        ! note: favre mean values, compute .../rho_mean in post-stage!
                        params_physics%k_mean = params_physics%k_mean + &
                                       0.5_rk * dummy(1) * ( dummy(2)**2.0_rk + dummy(3)**2.0_rk + dummy(4)**2.0_rk )

                        params_physics%UU_mean = params_physics%UU_mean + dummy(2)**2.0_rk + dummy(3)**2.0_rk + dummy(4)**2.0_rk

                    end do
                end do
            end do

        end if

        !-------------------------------------------------------------------------
        ! third: fourier coefficents
        !------------------------------------------------------------------------- 

        if ( params_physics%d == 3 ) then
            ! disable DFT for cases without forcing
            if (params_physics%forcing) then

                ! only 3D, remove mean velocities
                call compute_DFT( params_physics, phi(g+1:Bs(1)+g,g+1:Bs(2)+g,g+1:Bs(3)+g,UxF) &
                                                / phi(g+1:Bs(1)+g,g+1:Bs(2)+g,g+1:Bs(3)+g,rhoF) - params_physics%u_mean_0 &
                                                , phi_hat(:,:,:,1), x0, dx )
                call compute_DFT( params_physics, phi(g+1:Bs(1)+g,g+1:Bs(2)+g,g+1:Bs(3)+g,UyF) &
                                                / phi(g+1:Bs(1)+g,g+1:Bs(2)+g,g+1:Bs(3)+g,rhoF) - params_physics%v_mean_0 &
                                                , phi_hat(:,:,:,2), x0, dx )
                call compute_DFT( params_physics, phi(g+1:Bs(1)+g,g+1:Bs(2)+g,g+1:Bs(3)+g,UzF) &
                                                / phi(g+1:Bs(1)+g,g+1:Bs(2)+g,g+1:Bs(3)+g,rhoF) - params_physics%w_mean_0 &
                                                , phi_hat(:,:,:,3), x0, dx )

                params_physics%phi_hat = params_physics%phi_hat + phi_hat

            end if

        end if

        !-------------------------------------------------------------------------
        ! fourth: viscosity mean
        !-------------------------------------------------------------------------
        if ( params_physics%chemistry_model == "cantera" ) then 

            ! cantera chemistry
            if ( params_physics%d == 2 ) then
                ! /todo: 2D statistics

            else
                ! 3D statistics
                do k = g+1, Bs(3)+g-1
                    do j = g+1, Bs(2)+g-1
                        do i = g+1, Bs(1)+g-1

                            ! for readability
                            dummy(1) = phi(i,j,k,rhoF) * phi(i,j,k,rhoF)

                            ! mass fraction in primitive form
                            do n = 1, params_physics%species-1
                                Y(n) = phi(i,j,k,YF+n-1) / dummy(1)
                            end do
                            Y(params_physics%species) = 1.0_rk - sum(Y(1:params_physics%species-1))

                            ! mean molecular weigth
                            dummy(2) = 0.0_rk
                            do n = 1, params_physics%species
                                dummy(2) = dummy(2) + 1.0_rk / params_physics%Wk(n) * Y(n)
                            end do

                            ! mole fractions
                            do n = 1, params_physics%species
                                X(n) = 1.0_rk / dummy(2) / params_physics%Wk(n) * Y(n)
                            end do

                            ! set thermodynamic state
                            dummy(3) = phi(i,j,k,EF) / dummy(1) + sum( params_physics%dh(:) * Y(:) )
                            call setMoleFractions(gas, X)
                            call setState_UV(gas, dummy(3), 1.0_rk/dummy(1)  )

                            ! viscosity
                            params_physics%mu_mean = params_physics%mu_mean +  viscosity(gas)

                        end do
                    end do
                end do

            end if

        else
         
            ! non cantera chemistry -> use mu0 /todo: does not work in 2D!
            params_physics%mu_mean = params_physics%mu_mean + params_physics%mu0 * (Bs(1)-1) * (Bs(2)-1) * (Bs(3)-1)

        end if




    case ("post_stage")
        !-------------------------------------------------------------------------
        ! 3rd stage: post_stage.
        !-------------------------------------------------------------------------
        ! this stage is called only once, not for each block.

        !-------------------------------------------------------------------------
        ! forcing, solenodidal and dilatational fourier coefficients
        !-------------------------------------------------------------------------
        ! disable solenoidal/dilatational decomposition for cases without forcing
        if (params_physics%forcing) then

            ! sum local fourier coefficients
            phi_hat = params_physics%phi_hat
            call MPI_Allreduce(phi_hat, params_physics%phi_hat, kmax*kmax*kmax*3, MPI_DOUBLE_COMPLEX, MPI_SUM, WABBIT_COMM, ierr)
   
            ! compute solenoidal and dilatational components
            do i = 1, kmax
                do j = 1, kmax
                    do k = 1, kmax

                        ! set first mode to zero, condition from [Petersen2010]
                        if ( (i==1) .and. (j==1) .and. (k==1) ) then

                            params_physics%phi_hat(i,j,k,:)   = 0.0_rk
                            params_physics%phi_hat_d(i,j,k,:) = 0.0_rk
                            params_physics%phi_hat_s(i,j,k,:) = 0.0_rk

                        else

                            ! compute wave number and magnitude
                            wn(1) = 2.0_rk * pi / params_physics%L(1) * real(i-1, kind=rk)
                            wn(2) = 2.0_rk * pi / params_physics%L(2) * real(j-1, kind=rk)
                            wn(3) = 2.0_rk * pi / params_physics%L(3) * real(k-1, kind=rk)
                            dummy(1) = dsqrt( wn(1)**2.0_rk + wn(2)**2.0_rk + wn(3)**2.0_rk)

                            ! dot product between wavenumber and coeffcient vector
                            cdummy = wn(1) * params_physics%phi_hat(i,j,k,1) &
                                   + wn(2) * params_physics%phi_hat(i,j,k,2) &
                                   + wn(3) * params_physics%phi_hat(i,j,k,3)

                            ! dilatational component
                            params_physics%phi_hat_d(i,j,k,1) = wn(1) * cdummy / dummy(1)**2.0_rk
                            params_physics%phi_hat_d(i,j,k,2) = wn(2) * cdummy / dummy(1)**2.0_rk
                            params_physics%phi_hat_d(i,j,k,3) = wn(3) * cdummy / dummy(1)**2.0_rk

                            ! solenoidal component
                            wn2(1) = wn(3) * params_physics%phi_hat(i,j,k,2) - wn(2) * params_physics%phi_hat(i,j,k,3)
                            wn2(2) = wn(1) * params_physics%phi_hat(i,j,k,3) - wn(3) * params_physics%phi_hat(i,j,k,1)
                            wn2(3) = wn(2) * params_physics%phi_hat(i,j,k,1) - wn(1) * params_physics%phi_hat(i,j,k,2)

                            params_physics%phi_hat_s(i,j,k,1) = ( wn2(3) * wn(2) - wn2(2) * wn(3) ) / dummy(1)**2.0_rk
                            params_physics%phi_hat_s(i,j,k,2) = ( wn2(1) * wn(3) - wn2(3) * wn(1) ) / dummy(1)**2.0_rk
                            params_physics%phi_hat_s(i,j,k,3) = ( wn2(2) * wn(1) - wn2(1) * wn(2) ) / dummy(1)**2.0_rk

                        end if

                    end do
                end do
            end do

        end if

        !-------------------------------------------------------------------------
        ! domain size, note: all summations should exclude redundant nodes
        !-------------------------------------------------------------------------
        gridNodes = (Bs(1)-1) * (Bs(2)-1) * (Bs(3)-1) * lgt_n

        !-------------------------------------------------------------------------
        ! fluid properties
        !-------------------------------------------------------------------------
        ! density mean
        dummy(1) = params_physics%rho_mean / real(gridNodes, kind=rk)
        call MPI_Allreduce(dummy(1), params_physics%rho_mean, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, ierr)

        ! velocity means
        dummy(1) = params_physics%u_mean / real(gridNodes, kind=rk)
        call MPI_Allreduce(dummy(1), params_physics%u_mean, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, ierr)
        dummy(1) = params_physics%v_mean / real(gridNodes, kind=rk)
        call MPI_Allreduce(dummy(1), params_physics%v_mean, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, ierr)
        dummy(1) = params_physics%w_mean / real(gridNodes, kind=rk)
        call MPI_Allreduce(dummy(1), params_physics%w_mean, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, ierr)

        ! viscosity
        dummy(1) = params_physics%mu_mean / real(gridNodes, kind=rk)
        call MPI_Allreduce(dummy(1), params_physics%mu_mean, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, ierr)

        !-------------------------------------------------------------------------
        ! turbulent kinetic energy
        !-------------------------------------------------------------------------
        dummy(1) = params_physics%k_mean / real(gridNodes, kind=rk)
        call MPI_Allreduce(dummy(1), params_physics%k_mean, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, ierr)
        ! favre mean
        params_physics%k_mean = params_physics%k_mean / params_physics%rho_mean

        dummy(1) = params_physics%k_d_mean / real(gridNodes, kind=rk)
        call MPI_Allreduce(dummy(1), params_physics%k_d_mean, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, ierr)
        ! favre mean
        params_physics%k_d_mean = params_physics%k_d_mean / params_physics%rho_mean

        ! use k = k_d + k_s to compute k_s
        params_physics%k_s_mean = params_physics%k_mean - params_physics%k_d_mean

        !-----------------------
        ! turbulent dissipation
        !-----------------------
        dummy(1) = 4.0_rk/3.0_rk * params_physics%mu_mean * params_physics%epsilon_d_mean / real(gridNodes, kind=rk)
        call MPI_Allreduce(dummy(1), params_physics%epsilon_d_mean, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, ierr)
        params_physics%epsilon_d_mean = params_physics%epsilon_d_mean / params_physics%rho_mean

        dummy(1) = params_physics%mu_mean * params_physics%epsilon_s_mean / real(gridNodes, kind=rk)
        call MPI_Allreduce(dummy(1), params_physics%epsilon_s_mean, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, ierr)
        params_physics%epsilon_s_mean = params_physics%epsilon_s_mean / params_physics%rho_mean

        ! use eps = eps_d + eps_s to compute dissipation
        params_physics%epsilon_mean = params_physics%epsilon_d_mean + params_physics%epsilon_s_mean

        !-------------------------------------------------------------------------
        ! energy mean
        !-------------------------------------------------------------------------
        dummy(1) = params_physics%e_tot_mean
        call MPI_Allreduce(dummy(1), params_physics%e_tot_mean, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, ierr)

        dummy(1) = params_physics%e_kin_mean
        call MPI_Allreduce(dummy(1), params_physics%e_kin_mean, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, ierr)

        dummy(1) = params_physics%e_int_mean
        call MPI_Allreduce(dummy(1), params_physics%e_int_mean, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, ierr)

        !-------------------------------------------------------------------------
        ! u_rms square
        !-------------------------------------------------------------------------
        dummy(1) = params_physics%UU_mean / real(gridNodes, kind=rk)
        call MPI_Allreduce(dummy(1), params_physics%UU_mean, 1, MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, ierr)
       


        !-------------------------------------------------------------------------
        ! write statistics to ascii files.
        if ( params_physics%rank == 0) then

            ! write total energy to disk...
            open(14,file='e_tot.t',status='unknown',position='append')
            write (14,'(2(es15.8,1x))') time, params_physics%e_tot_mean
            close(14)

            open(14,file='e_kin_mean.t',status='unknown',position='append')
            write (14,'(2(es15.8,1x))') time, params_physics%e_kin_mean
            close(14)

            open(14,file='e_int_mean.t',status='unknown',position='append')
            write (14,'(2(es15.8,1x))') time, params_physics%e_int_mean
            close(14)

            ! write turbulent kinetic energy to disk...
            open(14,file='k_mean.t',status='unknown',position='append')
            write (14,'(2(es15.8,1x))') time, params_physics%k_mean
            close(14)
            ! write turbulent kinetic energy to disk...
            open(14,file='k_d.t',status='unknown',position='append')
            write (14,'(2(es15.8,1x))') time, params_physics%k_d_mean
            close(14)
            ! write turbulent kinetic energy to disk...
            open(14,file='k_s.t',status='unknown',position='append')
            write (14,'(2(es15.8,1x))') time, params_physics%k_s_mean
            close(14)

            ! write turbulent dissipation to disk...
            open(14,file='epsilon_mean.t',status='unknown',position='append')
            write (14,'(2(es15.8,1x))') time, params_physics%epsilon_mean
            close(14)
            ! write turbulent dissipation to disk...
            open(14,file='epsilon_d.t',status='unknown',position='append')
            write (14,'(2(es15.8,1x))') time, params_physics%epsilon_d_mean
            close(14)
            ! write turbulent dissipation to disk...
            open(14,file='epsilon_s.t',status='unknown',position='append')
            write (14,'(2(es15.8,1x))') time, params_physics%epsilon_s_mean
            close(14)

            ! write U square mean to disk...
            open(14,file='UU_mean.t',status='unknown',position='append')
            write (14,'(2(es15.8,1x))') time, params_physics%UU_mean
            close(14)

            ! write mu mean to disk...
            open(14,file='mu_mean.t',status='unknown',position='append')
            write (14,'(2(es15.8,1x))') time, params_physics%mu_mean
            close(14)

        end if

        ! save velocity means
        params_physics%u_mean_0 =  params_physics%u_mean
        params_physics%v_mean_0 =  params_physics%v_mean
        params_physics%w_mean_0 =  params_physics%w_mean

    case default
        call abort(181018001,"the statistics wrapper requests a stage this physics module cannot handle.")

    end select

end subroutine statistics_reactive_ns
