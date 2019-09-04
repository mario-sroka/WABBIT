!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name diff_wrapper_2D.f90
!> \version 0.5
!> \author msr
!
!> \brief wrapper for derivatives, note: pure periodic derivatives
!>
!! input:    - params_physics, data for one block, spacing parameter \n
!! output:   - computed rhs \n
!!
!!
!! = log ======================================================================================
!! \n
!! 10/10/18 - create
!
! ********************************************************************************************

subroutine  diff_wrapper_2D( Bs, g, u, deriv_type, onesided, dx, dy, dudx, dudy)

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters

!---------------------------------------------------------------------------------------------
! variables

    !> grid parameter
    integer(kind=ik), intent(inout)             :: g, Bs(3)

    !> datafield
    real(kind=rk), intent(inout)                :: u(Bs(1)+2*g, Bs(2)+2*g)

    !> derivative type
    character(len=*), intent(in)                :: deriv_type

    !> swithc for one sided derivatives
    logical, intent(in)                         :: onesided(2,3)

    !> spacing parameter
    real(kind=rk), intent(inout), optional      :: dx, dy
    
    !> results arrays
    real(kind=rk), intent(out), optional        :: dudx(Bs(1)+2*g, Bs(2)+2*g), dudy(Bs(1)+2*g, Bs(2)+2*g)

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    if ( onesided(1,1) .OR. onesided(2,1) .OR. onesided(1,2) .OR. onesided(2,2) ) then
        call abort(040919003,"ERROR: one sided derivation not implemented in 2D.")
    end if

!---------------------------------------------------------------------------------------------
! main body

    select case(deriv_type)

        ! u_x, u_y periodic
        case ("periodic_u_xy")
            ! periodic
            call diffxy_c_opt( Bs, g, dx, dy, u, dudx, dudy)

        ! special p
        case ("periodic_u_x_p")
            ! periodic
            call diffx_c_opt( Bs, g, dx, u, dudx)

        case ("periodic_u_y_p")
            ! periodic
            call diffy_c_opt( Bs, g, dy, u, dudy)

        ! u_x periodic
        case ("periodic_u_x")
            ! periodic
            call diffx_c_opt( Bs, g, dx, u, dudx)

        ! u_y periodic
        case ("periodic_u_y")
            ! periodic
            call diffy_c_opt( Bs, g, dy, u, dudy)

    end select

end subroutine diff_wrapper_2D
