!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name diff_wrapper_3D.f90
!> \version 0.5
!> \author msr
!
!> \brief wrapper for derivatives 
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

subroutine  diff_wrapper_3D( Bs, g, u, deriv_type, onesided, dx, dy, dz, dudx, dudy, dudz)

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters

!---------------------------------------------------------------------------------------------
! variables

    !> grid parameter
    integer(kind=ik), intent(inout)             :: g, Bs(3)

    !> datafield
    real(kind=rk), intent(inout)                :: u(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

    !> derivative type
    character(len=*), intent(in)                :: deriv_type

    !> swithc for one sided derivatives
    logical, intent(in)                         :: onesided(2,3)

    !> spacing parameter
    real(kind=rk), intent(inout), optional      :: dx, dy, dz
    
    !> results arrays
    real(kind=rk), intent(out), optional        :: dudx(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), &
                                                   dudy(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g), dudz(Bs(1)+2*g, Bs(2)+2*g, Bs(3)+2*g)

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    select case(deriv_type)

        ! u_x, u_y, u_z periodic
        case ("periodic_u_xyz")
            ! periodic
            call diffxyz_c_3D_opt( Bs, g, dx, dy, dz, u, dudx, dudy, dudz)

            if ( onesided(1,1) ) then
                ! onesided x-
                call diffx_oc_3D_opt_2( Bs, g, dx, u, dudx)
            end if

            if ( onesided(2,1) ) then
                ! onesided x+
                call diffx_co_3D_opt_2( Bs, g, dx, u, dudx)
            end if

            if ( onesided(1,2) ) then
                ! onesided y-
                call diffy_oc_3D_opt_2( Bs, g, dy, u, dudy)
            end if

            if ( onesided(2,2) ) then
                ! onesided y+
                call diffy_co_3D_opt_2( Bs, g, dy, u, dudy)
            end if

            if ( onesided(1,3) ) then
                ! onesided z-
                call diffz_oc_3D_opt_2( Bs, g, dz, u, dudz)
            end if

            if ( onesided(2,3) ) then
                ! onesided z+
                call diffz_co_3D_opt_2( Bs, g, dz, u, dudz)
            end if

        ! special p
        case ("periodic_u_x_p")
            if ( onesided(1,1) ) then
                ! onesided x-
                call diffx_oc_3D_opt_2( Bs, g, dx, u, dudx)
            elseif ( onesided(2,1) ) then
                ! onesided x+
                call diffx_co_3D_opt_2( Bs, g, dx, u, dudx)
            else
                call diffx_c_3D_opt( Bs, g, dx, u, dudx)
            end if

        case ("periodic_u_y_p")
            if ( onesided(1,2) ) then
                ! onesided y-
                call diffy_oc_3D_opt_2( Bs, g, dy, u, dudy)
            elseif ( onesided(2,2) ) then
                ! onesided y+
                call diffy_co_3D_opt_2( Bs, g, dy, u, dudy)
            else
                call diffy_c_3D_opt( Bs, g, dy, u, dudy)
            end if

        case ("periodic_u_z_p")

            if ( onesided(1,3) ) then
                ! onesided z-
                call diffz_oc_3D_opt_2( Bs, g, dz, u, dudz)
            elseif ( onesided(2,3) ) then
                ! onesided z+
                call diffz_co_3D_opt_2( Bs, g, dz, u, dudz)
            else
                call diffz_c_3D_opt( Bs, g, dz, u, dudz)
            end if

        ! u_x periodic
        case ("periodic_u_x")
            if ( onesided(1,1) ) then
                ! onesided x-
                call diffx_oc_3D_opt( Bs, g, dx, u, dudx)

            elseif ( onesided(2,1) ) then
                ! onesided x+
                call diffx_co_3D_opt( Bs, g, dx, u, dudx)

            else
                ! periodic
                call diffx_c_3D_opt( Bs, g, dx, u, dudx)

            end if

        ! u_y periodic
        case ("periodic_u_y")
            if ( onesided(1,2) ) then
                ! onesided y-
                call diffy_oc_3D_opt( Bs, g, dy, u, dudy)

            elseif ( onesided(2,2) ) then
                ! onesided y+
                call diffy_co_3D_opt( Bs, g, dy, u, dudy)

            else
                ! periodic
                call diffy_c_3D_opt( Bs, g, dy, u, dudy)

            end if

        ! u_z periodic
        case ("periodic_u_z")
            if ( onesided(1,3) ) then
                ! onesided z-
                call diffz_oc_3D_opt( Bs, g, dz, u, dudz)

            elseif ( onesided(2,3) ) then
                ! onesided z+
                call diffz_co_3D_opt( Bs, g, dz, u, dudz)

            else
                ! periodic
                call diffz_c_3D_opt( Bs, g, dz, u, dudz)

            end if

    end select

end subroutine diff_wrapper_3D
