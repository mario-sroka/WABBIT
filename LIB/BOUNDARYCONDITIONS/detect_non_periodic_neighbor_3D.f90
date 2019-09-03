!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name detect_non_periodic_neighbor_3D.f90
!> \version 0.5
!> \author msr
!
!> \brief returns TRUE/FALSE if neighboring treecode is over a periodic/non_periodic boundary
!> note: subrouitne works only for 3D 
!>
!! input:    - block treecode, neighbor treecode, domain periodicty \n
!! output:   - TRUE/FALSE \n
!!
!!
!! = log ======================================================================================
!! \n
!! 02/09/19 - create
!
! ********************************************************************************************

subroutine detect_non_periodic_neighbor_3D( my_treecode, neighbor_treecode, periodic_BC, level, L, Bs, data_bounds, non_periodic )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! treecode varaibles
    integer(kind=ik), intent(in)         :: my_treecode(:), neighbor_treecode(:)

    !> boundary condition wabbit core - note: TRUE: periodic BC, FALSE: non-periodic
    logical, intent(in)                  :: periodic_BC(3)

    !> mesh level 
    integer(kind=ik), intent(inout)      :: level

    !> domain sizes 
    real(kind=rk), intent(in)            :: L(3)

    !> block sizes
    integer(kind=ik), intent(in)         :: Bs(3)

    !> data_bounds array
    integer(kind=ik), intent(in)         :: data_bounds(2,3)

    !> neighbor over domain boundary, TRUE: non-periodic, FALSE: periodic
    logical, intent(out)                 :: non_periodic

    ! dummy variables
    integer(kind=ik)                     :: ix, iy, iz
    real(kind=rk)                        :: x01(3), x02(3), dx(3)

!---------------------------------------------------------------------------------------------
! variables initialization

    ! reset output variable
    non_periodic = .FALSE.

    ! zero: compute block origins and spacing
    ! /todo: use block coordinates instead of block origins?
    !-----------------------------------------------------------------------------------------
    dx = 2.0_rk**(-level) * L / real(Bs-1, kind=rk)

    call decoding( my_treecode(1:level), ix, iy, iz, level)
    x01 = real( ((/ix,iy,iz/) - 1)*(Bs-1) ,kind=rk) * dx

    call decoding( neighbor_treecode(1:level), ix, iy, iz, level)
    x02 = real( ((/ix,iy,iz/) - 1)*(Bs-1) ,kind=rk) * dx
    

!---------------------------------------------------------------------------------------------
! main body

    ! principle: for meshlevel 2 or higher compute origins (x01,x02) for active block and 
    ! neighbor block, then add Bs*dx to the largest x0i, compute abs(x01-x02) and if the 
    ! result is larger than the domain size, we have found a neighbor relation above the
    ! domain boundary
    ! for meshlevel 1 we need to check the data bounds additionally, because two blocks
    ! have allways two neighborhood relations: inside and above the domain boundary
    !-----------------------------------------------------------------------------------------
    
    ! first: test mesh level
    !-----------------------------------------------------------------------------------------
    if ( level > 1 ) then

        ! second: check non periodicity
        ! x direction
        !-------------------------------------------------------------------------------------
        ! non periodic boundary
        if ( .NOT.(periodic_BC(1)) ) then

            ! third: check block origins
            !---------------------------------------------------------------------------------
            call check_block_origins(x01(1), x02(1), dx(1), L(1), Bs(1), non_periodic)

        end if

        !-------------------------------------------------------------------------------------

        ! second: check non periodicity
        ! y direction
        !-------------------------------------------------------------------------------------
        ! non periodic boundary
        if ( .NOT.(periodic_BC(2)) ) then

            ! third: check block origins
            !---------------------------------------------------------------------------------
            call check_block_origins(x01(2), x02(2), dx(2), L(2), Bs(2), non_periodic)

        end if
        !-------------------------------------------------------------------------------------

        ! second: check non periodicity
        ! z direction
        !-------------------------------------------------------------------------------------
        ! non periodic boundary
        if ( .NOT.(periodic_BC(3)) ) then

            ! third: check block origins
            !---------------------------------------------------------------------------------
            call check_block_origins(x01(3), x02(3), dx(3), L(3), Bs(3), non_periodic)

        end if
        !-------------------------------------------------------------------------------------

    else
        
        ! second: check non periodicity
        ! x direction
        !-------------------------------------------------------------------------------------
        ! non periodic boundary
        if ( .NOT.(periodic_BC(1)) ) then

            ! third: check block origins and data bounds
            !---------------------------------------------------------------------------------
            ! x02 is receiver origin
            call check_databounds(x01(1), x02(1), data_bounds(2,1), Bs(1), non_periodic)

        end if

        !-------------------------------------------------------------------------------------

        ! second: check non periodicity
        ! y direction
        !-------------------------------------------------------------------------------------
        ! non periodic boundary
        if ( .NOT.(periodic_BC(2)) ) then

            ! third: check block origins
            !---------------------------------------------------------------------------------
            ! x02 is receiver origin
            call check_databounds(x01(2), x02(2), data_bounds(2,2), Bs(2), non_periodic)

        end if
        !-------------------------------------------------------------------------------------

        ! second: check non periodicity
        ! z direction
        !-------------------------------------------------------------------------------------
        ! non periodic boundary
        if ( .NOT.(periodic_BC(3)) ) then

            ! third: check block origins
            !---------------------------------------------------------------------------------
            ! x02 is receiver origin
            call check_databounds(x01(3), x02(3), data_bounds(2,3), Bs(3), non_periodic)

        end if
        !-------------------------------------------------------------------------------------
    end if

end subroutine detect_non_periodic_neighbor_3D

!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
! for readability
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------
subroutine check_block_origins(x01, x02, dx, L, Bs, non_periodic)
!-------------------------------------------------------------------------------------
    implicit none

    ! blocksizes, domain sizes, origins and spacing
    integer(kind=ik), intent(in)    :: Bs
    real(kind=rk), intent(inout)    :: x01, x02
    real(kind=rk), intent(in)       :: dx, L
    ! non periodicity
    logical, intent(inout)          :: non_periodic
!-------------------------------------------------------------------------------------

    ! check block origins
    !---------------------------------------------------------------------------------
    if ( x01 > x02 ) then
                
        ! x01 is closer to L
        x01 = x01 + real(Bs+1, kind=rk)*dx
    else

        ! x02 is closer to L
        x02 = x02 + real(Bs+1, kind=rk)*dx
    end if

    ! check distance between blocks
    !---------------------------------------------------------------------------------
    if ( abs(x01-x02) > L ) non_periodic = .TRUE.

end subroutine check_block_origins
!-------------------------------------------------------------------------------------

subroutine check_databounds(x01, x02, data_bound, Bs, non_periodic)
!-------------------------------------------------------------------------------------
    implicit none

    ! blocksizes, domain sizes, origins and spacing
    integer(kind=ik), intent(in)    :: Bs
    real(kind=rk), intent(inout)    :: x01, x02
    integer(kind=ik), intent(in)    :: data_bound
    ! non periodicity
    logical, intent(inout)          :: non_periodic
!-------------------------------------------------------------------------------------

    !---------------------------------------------------------------------------------
    ! x02 is receiver origin
    if ( x02 > x01 ) then

        ! check databounds
        ! x02 larger, so if receiver has Bs+g as databound, we have found the 
        ! domain boundary
        if ( data_bound > Bs ) non_periodic = .TRUE. 

    else

        ! check databounds
        ! x02 smaller, so if receiver has g+1 as databound, we have found the 
        ! domain boundary
        if ( data_bound < Bs ) non_periodic = .TRUE.

    end if

end subroutine check_databounds
!-------------------------------------------------------------------------------------
