!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name filter_1D.f90
!> \version 0.5
!> \author msr
!
!> \brief 1D filtering
!
!! input:    - data for one stencil, stencil \n
!! output:   - filtered data \n
!!
!!
!! = log ======================================================================================
!! \n
!! 15/10/18 - create
!
! ********************************************************************************************

subroutine filter_1D(phi, phi_tilde, a)

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> datafield
    real(kind=rk), intent(in)           :: phi(:)
    !> filtered value
    real(kind=rk), intent(out)          :: phi_tilde
    !> filter coefficients
    real(kind=rk), intent(in)           :: a(:)

    ! loop variable
    integer(kind=ik)                    :: k

    ! old values
    real(kind=rk)                       :: phi_old(size(phi,1))

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    phi_old   = phi
    phi_tilde = 0.0_rk

!---------------------------------------------------------------------------------------------
! main body

    ! check filter stencil
    if ( size(phi) /= size(a) ) then
        write(*,'(80("_"))')
        print*, phi
        print*, a
        call abort(123980,"ERROR: filter stencil has wrong size")
    end if

    ! filter data
    do k = 1, size(a)
        phi_tilde = phi_tilde + a(k)*phi_old(k)
    end do

end subroutine filter_1D
