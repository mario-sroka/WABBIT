!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name field_names_reactive_ns.f90
!> \version 0.5
!> \author msr
!
!> \brief return datafield label from physics parameter struct
!>
!! input:    - params_physics, datafield index \n
!! output:   - name of datafield \n
!!
!!
!! = log ======================================================================================
!! \n
!! 10/10/18 - create
!
! ********************************************************************************************

subroutine field_names_reactive_ns( params_physics, dF, dFname )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> parameter struct
    type(type_params_rns), intent(inout)    :: params_physics
    ! field index
    integer(kind=ik), intent(in)            :: dF
    ! returns the name of the field, corresponds to datafield index from ini file
    character(len=80), intent(out)          :: dFname

    ! loop parameter
    integer(kind=ik)                        :: k

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    ! use label from ini file
    if (dF <= params_physics%N_fields_saved) then
        dFname = params_physics%names_saved(dF)
    else
        ! error case
        call abort(121018001,"ERROR: can not find data field name for saving.")
    end if

end subroutine field_names_reactive_ns
