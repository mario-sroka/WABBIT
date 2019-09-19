!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name spectral_filter.f90
!> \version 0.5
!> \author msr
!
!> \brief spectral filtering
!> 3 stages: init   - reset fourier coefficients
!>           filter - compute fourier coefficients blockwise, delete specific coefficients
!>           post   - sum all coefficients with MPI, inverse fourier transformation blockwise 
!>
!! input:    - params_physics, data for one block, spacing parameter \n
!! output:   - filtered data fields\n
!!
!!
!! = log ======================================================================================
!! \n
!! 12/09/19 - create
!
! ********************************************************************************************

subroutine spectral_filter( params_physics, phi, phi_work, x0, dx, stage )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> parameter struct
    type(type_params_rns), intent(inout)    :: params_physics

    !> state vector for one block
    real(kind=rk), intent(inout)            :: phi(:, :, :, :)
    real(kind=rk), intent(inout)            :: phi_work(:, :, :, :)

    !> spacing and origin of block
    real(kind=rk), intent(in)               :: x0(1:3), dx(1:3)

    ! stage. there is 3 stages, init_stage, filter_stage and local_stage. If the PDE has
    ! terms that depend on global qtys, such as forces etc, which cannot be computed
    ! from a single block alone, the first stage does that. the second stage can then
    ! use these integral qtys for the actual RHS evaluation.
    character(len=*), intent(in)            :: stage

    ! grid parameter
    integer(kind=ik)                        :: Bs(3), g, Ds(3)
    ! field indexes
    integer(kind=ik)                        :: rhoF, UxF, UyF, UzF, EF
    ! fourier coeficients
    complex(kind=rk), allocatable, save     :: phi_hat(:,:,:,:)
    ! wavenumber limit, wavenumber
    real(kind=rk)                           :: w0, w
    ! loop parameter
    integer(kind=ik)                        :: i, j, k
    ! MPI error
    integer(kind=ik)                        :: ierr

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

    ! domain sizes
    Ds(1) = size(params_physics%phi_hat, 1 )
    Ds(2) = size(params_physics%phi_hat, 2 )
    Ds(3) = size(params_physics%phi_hat, 3 )

    ! local coefficients array
    if ( .not. allocated(phi_hat) ) allocate(phi_hat(Ds(1), Ds(2), Ds(3), params_physics%NdF))

!---------------------------------------------------------------------------------------------
! main body

    select case(stage)
        case ("init_stage")
        !-------------------------------------------------------------------------
        ! 1st stage: init_stage.
        !-------------------------------------------------------------------------
        ! this stage is called only once, not for each block.
        ! performs initializations, such as resetting integrals
       
        ! reset fourier coefficients
        ! note: make sure, old values are obsolete!
        params_physics%phi_hat          = 0.0_rk

    case ("filter_stage")
        !-------------------------------------------------------------------------
        ! 2nd stage: integral_stage.
        !-------------------------------------------------------------------------
        ! This stage contains all operations which are running on the blocks
        ! called for each block
  
        ! compute DFT, use DFT_full subroutine to compute all coefficients, /todo: remove?
        if ( params_physics%d == 3 ) then
            ! only 3D, remove mean velocities
            ! *** rho ***
            call compute_DFT_full( params_physics, phi(g+1:Bs(1)+g,g+1:Bs(2)+g,g+1:Bs(3)+g,rhoF)**2.0_rk &
                                                 , phi_hat(:,:,:,rhoF), x0, dx )
            ! *** Ux ***
            call compute_DFT_full( params_physics, phi(g+1:Bs(1)+g,g+1:Bs(2)+g,g+1:Bs(3)+g,UxF) &
                                                 / phi(g+1:Bs(1)+g,g+1:Bs(2)+g,g+1:Bs(3)+g,rhoF) - params_physics%u_mean_0 &
                                                 , phi_hat(:,:,:,UxF), x0, dx )
            ! *** Uy ***
            call compute_DFT_full( params_physics, phi(g+1:Bs(1)+g,g+1:Bs(2)+g,g+1:Bs(3)+g,UyF) &
                                                 / phi(g+1:Bs(1)+g,g+1:Bs(2)+g,g+1:Bs(3)+g,rhoF) - params_physics%v_mean_0 &
                                                 , phi_hat(:,:,:,UyF), x0, dx )
            ! *** Uz ***
            call compute_DFT_full( params_physics, phi(g+1:Bs(1)+g,g+1:Bs(2)+g,g+1:Bs(3)+g,UzF) &
                                                 / phi(g+1:Bs(1)+g,g+1:Bs(2)+g,g+1:Bs(3)+g,rhoF) - params_physics%w_mean_0 &
                                                 , phi_hat(:,:,:,UzF), x0, dx )
            ! *** EF ***
            call compute_DFT_full( params_physics, phi(g+1:Bs(1)+g,g+1:Bs(2)+g,g+1:Bs(3)+g,EF) &
                                                 / phi(g+1:Bs(1)+g,g+1:Bs(2)+g,g+1:Bs(3)+g,rhoF)**2.0_rk &
                                                 , phi_hat(:,:,:,EF), x0, dx )
            ! /todo species filtering?
    
            params_physics%phi_hat = params_physics%phi_hat + phi_hat

        end if

    case ("post_stage")
        !-------------------------------------------------------------------------
        ! 3rd stage: post_stage.
        !-------------------------------------------------------------------------
        ! this stage is called for each block.

        ! sum all coefficients
        phi_hat = params_physics%phi_hat
        call MPI_Allreduce(phi_hat, params_physics%phi_hat, Ds(1)*Ds(2)*Ds(3)*3, MPI_DOUBLE_COMPLEX, MPI_SUM, WABBIT_COMM, ierr)

        ! filter phi hat, wavenumber limit
        w0 = params_physics%spectral_w_limit
        do i = 1, Ds(1)
            do j = 1, Ds(2)
                do k = 1, Ds(3)

                    ! compute wave number
                    w = dsqrt( (2.0_rk * pi / params_physics%L(1) * real(i-1, kind=rk))**2.0_rk &
                             + (2.0_rk * pi / params_physics%L(2) * real(j-1, kind=rk))**2.0_rk &
                             + (2.0_rk * pi / params_physics%L(3) * real(k-1, kind=rk))**2.0_rk)

                    if ( w > w0  ) then
                        phi_hat(i,j,k,:) = 0.0_rk
                    end if

                end do
            end do
        end do
 
        ! IDFT
        call compute_IDFT_full( params_physics, params_physics%phi_hat(:,:,:,rhoF), phi_work(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 1), x0, dx )
        call compute_IDFT_full( params_physics, params_physics%phi_hat(:,:,:,UxF),  phi_work(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 2), x0, dx )
        call compute_IDFT_full( params_physics, params_physics%phi_hat(:,:,:,UyF),  phi_work(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 3), x0, dx )
        call compute_IDFT_full( params_physics, params_physics%phi_hat(:,:,:,UzF),  phi_work(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 2), x0, dx )
        call compute_IDFT_full( params_physics, params_physics%phi_hat(:,:,:,EF),   phi_work(g+1:Bs(1)+g, g+1:Bs(2)+g, g+1:Bs(3)+g, 3), x0, dx )

        ! skew symmetric
        call  convert_from_primitive( params_physics, phi_work, phi )

    case default
        call abort(181018001,"the filter wrapper requests a stage this physics module cannot handle.")

    end select

end subroutine spectral_filter
