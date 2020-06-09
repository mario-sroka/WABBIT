!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name bogey_filter.f90
!> \version 0.5
!> \author msr
!
!> \brief filter subroutine for bogey filtering
!
!! input:    - params_physics, data for one block, work_array \n
!! output:   - filtered data \n
!!
!!
!! = log ======================================================================================
!! \n
!! 08/06/20 - fresh copy from navier stokes module 
!
! ********************************************************************************************

subroutine bogey_filter( params_physics, phi, phi_work, x0, dx )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> parameter struct
    type(type_params_rns), intent(inout)    :: params_physics

    !> state vector for one block, work array
    real(kind=rk), intent(inout)            :: phi(:, :, :, :), phi_work(:, :, :, :)
    !> spacing and origin of block
    real(kind=rk), intent(in)               :: x0(1:3), dx(1:3)

    !> grid parameter
    integer(kind=ik)                        :: Bs(3), g, NdF
    ! loop variables
    integer(kind=ik)                        :: i, j, k, l, d, dF
    ! field indexes
    integer(kind=ik)                        :: rhoF, UxF, UyF, UzF, EF, YF

    ! stencils
    integer(kind=ik)                        :: stencil_size
    real(kind=rk)                           :: stencil(3), c_stencil(4), c1=-0.210383_rk, c2 = 0.030617_rk

    ! filter parameter: threshold 
    real(kind=rk)                           :: r_th
    character(len=80)                       :: detector_method, sigma_method

    ! dummy field for old block data
    real(kind=rk), allocatable, save        :: dummy(:,:,:,:)

    ! isotropic coefficient for speed of sound computation
    real(kind=rk)                           :: gamma_

    ! parameter
    integer(kind=ik), parameter             :: SHIFT(7)=(/-3, -2, -1, 0, 1, 2, 3/)
    real(kind=rk), save                     :: divu(size(SHIFT),3), soundspeed2(size(SHIFT),3), &
                                               Dtheta(size(SHIFT),3), sigma(size(SHIFT),3)
    real(kind=rk)                           :: r, Dthetamag

!---------------------------------------------------------------------------------------------
! variables initialization

    Bs  = params_physics%Bs
    g   = params_physics%g
    NdF = params_physics%NdF

    ! field indexes from params
    rhoF = params_physics%rhoF
    UxF  = params_physics%UxF
    UyF  = params_physics%UyF
    UzF  = params_physics%UzF
    EF   = params_physics%EF

    stencil_size = 3
    stencil(1:stencil_size) = (/  1.0_rk/  4.0_rk, -1.0_rk/  2.0_rk, 1.0_rk/  4.0_rk/)
    c_stencil = (/ -c2, -c1, c1, c2 /)

    r_th = params_physics%r_th
    detector_method = params_physics%bogey_detector ! "p", "divU"
    sigma_method    = params_physics%bogey_sigma ! "abs", "tanh"

    gamma_ = params_physics%gamma_

    ! dummy array
    if (params_physics%d == 3 ) then
        if ( .not. allocated(dummy) ) allocate(dummy(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g, params_physics%NdF))
    else
        if ( .not. allocated(dummy) ) allocate(dummy(Bs(1)+2*g, Bs(2)+2*g, 1, params_physics%NdF))
    end if

!---------------------------------------------------------------------------------------------
! main body

    ! convert from skew symmetric form
    call convert_to_primitive(params_physics, phi, phi_work(:, :, :, 1:params_physics%NdF), 1.0_rk/phi(:,:,:,rhoF) )
    ! convert to conservative and store in dummy field
    call convert_to_conservative( params_physics, phi, dummy )

    ! 3D or 2D case
    if (params_physics%d == 3 ) then

        call abort(080620001," ERROR: in filter reactive_ns, 3D bogey filter not implemented yet.")

    else

        k = 1
        do j = g+1, Bs(2)+g
            do i = g+1, Bs(1)+g

                !This loop shifts the stencil in every direction   direction shift(1)=-1,shift(2)=0,shift(3)=1
                do l = 1, 7 !shift the stencil loop

                    ! compute the divergence
                    ! divergence shifted in x
                    divu(l,1) = ( phi_work(i+1+SHIFT(l), j, k, UxF) - phi_work(i-1+SHIFT(l), j, k, UxF) ) / (2.0_rk*dx(1)) &
                              + ( phi_work(i+SHIFT(l), j+1, k, UyF) - phi_work(i+SHIFT(l), j-1, k, UyF) ) / (2.0_rk*dx(2))
                    ! shifted in y
                    divu(l,2) = ( phi_work(i+1, j+SHIFT(l), k, UxF) - phi_work(i-1, j+SHIFT(l), k, UxF) ) / (2.0_rk*dx(1)) &
                              + ( phi_work(i, j+SHIFT(l)+1, k, UyF) - phi_work(i, j+SHIFT(l)-1, k, UyF) ) / (2.0_rk*dx(2))
    
                    ! speed of sound is actually not needed for i=1,2 and 6,7
                    soundspeed2(l,1)= gamma_ * phi_work(i+SHIFT(l), j, k, EF) / phi_work(i+SHIFT(l), j, k, rhoF) ! c^2 shifted in x
                    soundspeed2(l,2)= gamma_ * phi_work(i, j+SHIFT(l), k, EF) / phi_work(i, j+SHIFT(l), k, rhoF) ! c^2 shifted in y

                end do

                ! from the dilatation we can compute the filtersing strength sigma:
                do d = 1, 2
    
                    ! compute a second order filter as the filtering indicator
                    do l = 2, 6
                        call filter_1D( divu(l-1:l+1, d), Dtheta(l,d), stencil(1:stencil_size) )
                    end do    

                    do l = 3, 5
                        ! compute its magnitude
                        Dthetamag = 0.5_rk * ( ( Dtheta(l,d) - Dtheta(l+1,d) )**2.0_rk &
                                  + ( Dtheta(l,d) - Dtheta(l-1,d) )**2.0_rk )
                        ! normalize the magnitude to the local speed of sound
                        r = Dthetamag / (soundspeed2(l,d)/dx(d)**2.0_rk) + 1e-16_rk
                        ! compute the filtering strength
                        sigma(l,d) = 1.0_rk - dtanh( r_th/r/0.7_rk )
                    end do

                    l = 4 ! SHIFT(4)=0 (no shift)

                end do ! loop over dimensions

                ! filter all datafields
                do dF = 1, NdF
                    phi(i, j, k, dF) = dummy(i, j, k, dF) &
                    !filtering in x
                    - ( 0.5_rk * (sigma(l+1,1) + sigma(l,1)) * sum( c_stencil * dummy(i-1:i+2, j, k, dF) ) &
                    -   0.5_rk * (sigma(l-1,1) + sigma(l,1)) * sum( c_stencil * dummy(i-2:i+1, j, k, dF) ) ) &
                    !filtering in y
                    - ( 0.5_rk * (sigma(l+1,2) + sigma(l,2)) * sum( c_stencil * dummy(i, j-1:j+2, k, dF) ) &
                    -   0.5_rk * (sigma(l-1,2) + sigma(l,2)) * sum( c_stencil * dummy(i, j-2:j+1, k, dF) ) )
                end do 

            end do
        end do

    endif

    ! convert back to skew symmetric form and write filtered values back to state vector
    call convert_from_conservative( params_physics, phi, phi_work )
    phi(:,:,:,:) = phi_work(:,:,:,1:params_physics%NdF)

end subroutine bogey_filter
