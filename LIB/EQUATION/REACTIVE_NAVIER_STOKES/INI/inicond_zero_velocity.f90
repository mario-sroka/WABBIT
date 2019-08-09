!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name inicond_zero_velocity.f90
!> \version 0.5
!> \author msr
!>
!> \brief set velocity fields to zero
!>
!! input:    - params_physics, state vector \n
!! output:   - state vector \n
!!
!! = log ======================================================================================
!! \n
!! 09/08/19 - create
!
! ********************************************************************************************

subroutine inicond_zero_velocity( params_physics, phi )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> parameter struct
    type(type_params_rns), intent(inout)    :: params_physics

    !> actual block data (note this routine acts only on one block)
    real(kind=rk), intent(inout)            :: phi(:,:,:,:)

    ! field indexes
    integer(kind=ik)                        :: UxF, UyF, UzF

!---------------------------------------------------------------------------------------------
! variables initialization

    ! field indexes from params
    UxF  = params_physics%UxF
    UyF  = params_physics%UyF
    UzF  = params_physics%UzF

!---------------------------------------------------------------------------------------------
! main body
    
    ! 3D or 2D?
    if (params_physics%d == 3) then
        phi( :, :, :, UxF) = 0.0_rk
        phi( :, :, :, UyF) = 0.0_rk
        phi( :, :, :, UzF) = 0.0_rk
    else
        phi( :, :, :, UxF) = 0.0_rk
        phi( :, :, :, UyF) = 0.0_rk
    endif

end subroutine inicond_zero_velocity

