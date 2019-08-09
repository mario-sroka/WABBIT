!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name inicond_blob.f90
!> \version 0.5
!> \author msr
!
!> \brief initialize blob for 2D/3D navier stokes physics
!> note: blob is located in domain at given coordinates and with given scales
!>
!> note: in order to use the subroutine more generic, same parameters from the parameter struct 
!> are given explicitly
!>
!! input:    - params, grid and spacing parameter, blob field \n
!! output:   - block data \n
!!
!!
!! = log ======================================================================================
!! \n
!! 11/10/18 - create
!
! ********************************************************************************************

subroutine inicond_blob( params_physics, blob_field, x0, dx, d, location, scales )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> parameter struct
    type(type_params_rns), intent(inout)    :: params_physics

    !> blob field
    real(kind=rk), intent(inout)            :: blob_field(:,:,:)

    !> spacing and origin of block
    real(kind=rk), intent(in)               :: x0(3), dx(3)

    !> dimension
    integer(kind=ik), intent(in)            :: d

    !> location and blob scales
    real(kind=rk), intent(in)               :: location(3), scales(2)

    ! grid parameter
    integer(kind=ik)                        :: Bs(3), g
    ! coordinates and domain parameter
    real(kind=rk)                           :: x, y, z
    ! loop variables
    integer(kind=ik)                        :: i, j, k

!---------------------------------------------------------------------------------------------
! variables initialization

    Bs      = params_physics%Bs
    g       = params_physics%g

!---------------------------------------------------------------------------------------------
! main body

    if (d == 2) then
        ! 2D case
        ! create gauss pulse
        do i = g+1, Bs(1)+g
            do j = g+1, Bs(2)+g

                ! compute x,y coordinates from spacing and origin
                x = dble(i-(g+1)) * dx(1) + x0(1)
                y = dble(j-(g+1)) * dx(2) + x0(2)

                ! set gauss blob
                blob_field(i, j, 1) = dexp( -( (x-location(1))**2.0_rk + (y-location(2))**2.0_rk ) / (scales(1)*scales(2)**2.0_rk) )

            end do
        end do

    else
        ! 3D case
        ! create gauss pulse
        do i = g+1, Bs(1)+g
            do j = g+1, Bs(2)+g
                do k = g+1, Bs(3)+g

                    x = dble(i-(g+1)) * dx(1) + x0(1)
                    y = dble(j-(g+1)) * dx(2) + x0(2)
                    z = dble(k-(g+1)) * dx(3) + x0(3)

                    ! set gauss blob
                    blob_field(i, j, k) = dexp( -( (x-location(1))**2.0_rk + (y-location(2))**2.0_rk + &
                                                   (z-location(3))**2.0_rk ) / (scales(1)*scales(2)**2.0_rk) )

                end do
            end do
        end do

    endif

end subroutine inicond_blob

