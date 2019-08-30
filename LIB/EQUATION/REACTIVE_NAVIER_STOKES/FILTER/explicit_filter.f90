!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name explicit_filter.f90
!> \version 0.5
!> \author msr
!
!> \brief filter subroutine for explicit filtering
!
!! input:    - params_physics, data for one block, work_array \n
!! output:   - filtered data \n
!!
!!
!! = log ======================================================================================
!! \n
!! 15/10/18 - create
!
! ********************************************************************************************

subroutine explicit_filter( params_physics, phi, phi_work, x0, dx )

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
    integer(kind=ik)                        :: i, j, k, dF
    ! filter size
    integer(kind=ik)                        :: fSize

    ! filtered values and array for old block data
    real(kind=rk)                           :: phi_tilde(3)

!---------------------------------------------------------------------------------------------
! variables initialization

    fSize = (size(params_physics%filter_stencil,1)+1)/2-1

    Bs  = params_physics%Bs
    g   = params_physics%g
    NdF = params_physics%NdF

!---------------------------------------------------------------------------------------------
! main body

    do dF = 1, NdF
        ! 3D or 2D case
        if (params_physics%d == 3 ) then

            ! 3D
            ! loop over block data
            do i = g+1, Bs(1)+g
                do j = g+1, Bs(2)+g
                    do k = g+1, Bs(3)+g
                        ! x direction
                        call filter_1D( phi(i-fSize:i+fSize, j, k, dF ), phi_tilde(1), params_physics%filter_stencil )
                        ! y direction
                        call filter_1D( phi(i, j-fSize:j+fSize, k, dF ), phi_tilde(2), params_physics%filter_stencil )
                        ! z direction
                        call filter_1D( phi(i, j, k-fSize:k+fSize, dF ), phi_tilde(3), params_physics%filter_stencil )
                        ! filter
                        phi_work(i, j, k, dF ) = phi(i, j, k, dF ) + 1.0_rk/3.0_rk*(phi_tilde(1) + phi_tilde(2) + phi_tilde(3))
                    end do
                end do
            end do

            ! write filtered values back to state vector
            do j = g+1, Bs(2)+g
                do i = g+1, Bs(1)+g
                    do k = g+1, Bs(3)+g
                        phi(i, j, k, dF)   = phi_work(i, j, k, dF)
                    end do
                end do
            end do

        else

            ! 2D
            ! loop over block data
            do i = g+1, Bs(1)+g
                do j = g+1, Bs(2)+g
                    ! x direction
                    call filter_1D( phi(i-fSize:i+fSize, j, 1, dF ), phi_tilde(1), params_physics%filter_stencil )
                    ! y direction
                    call filter_1D( phi(i, j-fSize:j+fSize, 1, dF ), phi_tilde(2), params_physics%filter_stencil )
                    ! filter
                    phi_work(i, j, 1, dF ) = phi(i, j, 1, dF) + 1.0_rk/2.0_rk*(phi_tilde(1) + phi_tilde(2))
                end do
            end do

            ! write filtered values back to state vector
            do j = g+1, Bs(2)+g
                do i = g+1, Bs(1)+g
                    phi(i, j, 1, dF)   = phi_work(i, j, 1, dF)
                end do
            end do

        endif
    end do

end subroutine explicit_filter
