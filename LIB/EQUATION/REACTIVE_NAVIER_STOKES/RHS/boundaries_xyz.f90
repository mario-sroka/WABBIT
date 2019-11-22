!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name boundaries_xyz.f90
!> \version 0.5
!> \author msr
!
!> \brief set boundary switch for non periodic boundaries, subroutine works for on single block 
!> /todo: maybe store this switch to prevent computation in every RHS call
!!
!! input:    - grid and domain parameter \n
!! output:   - boundary switch \n
!!
!!
!! = log ======================================================================================
!! \n
!! 04/07/19 - create
!
! ********************************************************************************************

subroutine boundaries_xyz( onesided, x0, dx, L, Bs, g, periodic_BC)

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> boundary switch: FALSE: periodic boundary, TRUE: non-periodic boundary
    logical, intent(inout)               :: onesided(2,3)

    !> block origin, spacing, domain size
    real(kind=rk), intent(in)            :: x0(3), dx(3), L(3)

    !> grid parameter
    integer(kind=ik), intent(in)         :: Bs(3), g

    !> boundary condition from wabbit core - note: TRUE: periodic BC, FALSE: non-periodic
    logical, intent(in)                  :: periodic_BC(3)

!---------------------------------------------------------------------------------------------
! variables initialization

    ! reset switch
    onesided = .false.

!---------------------------------------------------------------------------------------------
! main body
  
    ! x+
    if ( abs( x0(1) + dx(1)*real(Bs(1)-1,kind=rk)  - L(1) ) < 1e-12_rk ) then
        onesided(2,1) = .NOT.(periodic_BC(1))
    end if
    ! x-
    if ( ( ( x0(1) - dx(1) ) - 1e-16_rk ) < 0.0_rk ) then
        onesided(1,1) = .NOT.(periodic_BC(1))
    end if

    ! y+
    if ( abs( x0(2) + dx(2)*real(Bs(2)-1,kind=rk)  - L(2) ) < 1e-12_rk ) then
        onesided(2,2) = .NOT.(periodic_BC(2))
    end if
    ! y-
    if ( ( ( x0(2) - dx(2) ) - 1e-16_rk ) < 0.0_rk ) then
        onesided(1,2) = .NOT.(periodic_BC(2))
    end if

    ! z+
    if ( abs( x0(3) + dx(3)*real(Bs(3)-1,kind=rk)  - L(3) ) < 1e-12_rk ) then
        onesided(2,3) = .NOT.(periodic_BC(3))
    end if
    ! z-
    if ( ( ( x0(3) - dx(3) ) - 1e-16_rk ) < 0.0_rk ) then
        onesided(1,3) = .NOT.(periodic_BC(3))
    end if

end subroutine boundaries_xyz
