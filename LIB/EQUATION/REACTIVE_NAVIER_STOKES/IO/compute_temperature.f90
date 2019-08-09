!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name compute_temperature.f90
!> \version 0.5
!> \author msr
!
!> \brief compute temperature with cantera function 
!>
!! input:    - params_physics, state vector \n
!! output:   - work array \n
!!
!!
!! = log ======================================================================================
!! \n
!! 09/08/19 - create
!
! ********************************************************************************************

subroutine compute_temperature( params_physics, phi, phi_work, gas )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> parameter struct
    type(type_params_rns), intent(inout)    :: params_physics

    !> state vector for one block, work array
    real(kind=rk), intent(inout)            :: phi(:, :, :, :), phi_work(:, :, :)

    !> Cantera gas mixture struct
    type(phase_t), intent(inout)            :: gas

    ! field indexes
    integer(kind=ik)                        :: rhoF, EF, YF

    ! grid parameter
    integer(kind=ik)                        :: Bs(3), g
    ! loop variables
    integer(kind=ik)                        :: i, j, k, l
    ! dummy 
    real(kind=rk)                           :: dummy, dummy2(params_physics%species)

!---------------------------------------------------------------------------------------------
! variables initialization

    ! field indexes from params
    rhoF = params_physics%rhoF
    EF   = params_physics%EF
    YF   = params_physics%YF

    ! grid parameter
    Bs = params_physics%Bs
    g  = params_physics%g

!---------------------------------------------------------------------------------------------
! main body

     if (params_physics%d == 2) then

         ! loop over all nodes and use cantera function
         do j = 1, Bs(1)+g+g
             do i = 1, Bs(2)+g+g

                 dummy = 0.0_rk
                 dummy2(params_physics%species) = 1.0_rk

                 ! compute primitive mass fractions
                 do l = 1, params_physics%species-1
                     dummy2(l) = phi(i,j,1,YF+l-1) / phi(i,j,1,rhoF)**2.0_rk
                     dummy2(params_physics%species) = dummy2(params_physics%species) - dummy2(l)
                 end do

                 ! set gas mixture and thermodynamic state
                 call setMassFractions(gas, dummy2)
                 call setState_UV(gas, phi(i,j,1,EF)/phi(i,j,1,rhoF)**2.0_rk + sum( params_physics%dh(:) * dummy2(:) ), &
                                  1.0_rk/phi(i,j,1,rhoF)**2.0_rk )

                 ! temperature from cantera
                 phi_work(i, j, 1) = temperature(gas)
               
             end do
         end do

     else

         ! loop over all nodes and use cantera function
         do j = 1, Bs(1)+g+g
             do i = 1, Bs(2)+g+g
                 do k = 1, Bs(3)+g+g

                     dummy = 0.0_rk
                     dummy2(params_physics%species) = 1.0_rk

                     ! compute primitive mass fractions
                     do l = 1, params_physics%species-1
                         dummy2(l) = phi(i,j,k,YF+l-1) / phi(i,j,k,rhoF)**2.0_rk
                         dummy2(params_physics%species) = dummy2(params_physics%species) - dummy2(l)
                     end do

                     ! set gas mixture and thermodynamic state
                     call setMassFractions(gas, dummy2)
                     call setState_UV(gas, phi(i,j,k,EF)/phi(i,j,k,rhoF)**2.0_rk + sum( params_physics%dh(:) * dummy2(:) ), &
                                      1.0_rk/phi(i,j,k,rhoF)**2.0_rk )

                     ! temperature from cantera
                     phi_work(i, j, k) = temperature(gas)
               
                 end do
             end do
         end do

     end if

end subroutine compute_temperature
