!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name inicond_spark.f90
!> \version 0.5
!> \author msr
!>
!> \brief initialize rho, es and species with given premixed gas conditions, places a spark 
!> (temperature blob) and burned gas into domain
!> note: subroutine returns skew symmetric values!
!>
!! input:    - params_physics, state vector, gas mixture \n
!! output:   - state vector \n
!!
!!
!! = log ======================================================================================
!! \n
!! 25/01/19 - create
!
! ********************************************************************************************

subroutine inicond_spark( params_physics, phi, x0, dx, gas )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> parameter struct
    type(type_params_rns), intent(inout)    :: params_physics

    !> actual block data (note this routine acts only on one block)
    real(kind=rk), intent(inout)            :: phi(:,:,:,:)

    !> spacing and origin of block
    real(kind=rk), intent(in)               :: x0(1:3), dx(1:3)

    !> Cantera gas mixture struct
    type(phase_t), intent(inout)            :: gas

    ! grid parameter
    integer(kind=ik)                        :: Bs(3), g
    ! coordinates and domain parameter
    real(kind=rk)                           :: x, y, z, L(params_physics%d)
    ! gauss blob center coordinates, width, size, radius
    real(kind=rk)                           :: mux(params_physics%d), sigma, omega, r
    ! loop variables
    integer(kind=ik)                        :: i, j, k, n
    ! field indexes
    integer(kind=ik)                        :: rhoF, UxF, UyF, UzF, EF, YF

    ! dummy fields and variables
    real(kind=rk), allocatable              :: blob_field(:, :, :), blob_field2(:, :, :)
    real(kind=rk)                           :: burned(params_physics%NdF-params_physics%d), unburned(params_physics%NdF-params_physics%d)

!---------------------------------------------------------------------------------------------
! variables initialization

    Bs      = params_physics%Bs
    g       = params_physics%g
    L       = params_physics%L(1:params_physics%d)

    ! use inicond width from ini file
    sigma   = params_physics%inicond_scales(1)
    omega   = params_physics%inicond_scales(2)

    ! field indexes from params
    rhoF = params_physics%rhoF
    UxF  = params_physics%UxF
    UyF  = params_physics%UyF
    UzF  = params_physics%UzF
    EF   = params_physics%EF
    YF   = params_physics%YF

    ! dummy fields
    allocate( blob_field(1:Bs(1)+2*g, 1:Bs(2)+2*g, 1:Bs(3)+2*g), blob_field2(1:Bs(1)+2*g, 1:Bs(2)+2*g, 1:Bs(3)+2*g) )

!---------------------------------------------------------------------------------------------
! main body

    ! center point of the pressure blob
    mux(:) = params_physics%inicond_position(:)

    ! create blob fields
    call inicond_blob( params_physics, blob_field , x0, dx, params_physics%d, mux, (/ 2.0_rk, omega /) )
    call inicond_blob( params_physics, blob_field2, x0, dx, params_physics%d, mux, (/ sigma , omega /) )

    ! blobs all over the world
    ! ------------------------
    ! unburned conditions
    call setState_TPX(gas, params_physics%T0, params_physics%inicond_p, params_physics%inicond_X )

    unburned(1) = dsqrt(density(gas))
    unburned(2) = ( intEnergy_mass(gas) - sum( params_physics%dh(:) * params_physics%inicond_Y(:) ) ) * density(gas)
    do n = 1, params_physics%species
        unburned(2+n) = params_physics%inicond_Y(n) * density(gas)
    end do

    ! burned conditions
    call setState_TPX(gas, params_physics%inicond_T, params_physics%inicond_p, params_physics%inicond_X )
    call equilibrate(gas, 'TP')
    call getMassFractions(gas, burned(3:size(burned)))

    burned(1) = dsqrt(density(gas))
    burned(2) = ( intEnergy_mass(gas) - sum( params_physics%dh(:) * burned(3:size(burned)) ) ) * density(gas)
    do n = 1, params_physics%species
        burned(2+n) = burned(2+n) * density(gas)
    end do

    ! blobing
    phi( :, :, :, rhoF) = unburned(1) + (burned(1) - unburned(1)) * blob_field(:, :, :)
    phi( :, :, :, EF)   = unburned(2) + (burned(2) - unburned(2)) * blob_field(:, :, :)
    do n = 1, params_physics%species
        phi( :, :, :, YF+n-1) = unburned(2+n) + (burned(2+n) - unburned(2+n)) * blob_field2(:, :, :)
    end do

    ! skew symmetric velocity
    phi( :, :, :, UxF) = phi( :, :, :, UxF) * phi( :, :, :, rhoF)
    phi( :, :, :, UyF) = phi( :, :, :, UyF) * phi( :, :, :, rhoF)
    if (params_physics%d == 3) phi( :, :, :, UzF) = phi( :, :, :, UzF) * phi( :, :, :, rhoF)

    ! clean up
    deallocate(blob_field, blob_field2)

end subroutine inicond_spark

