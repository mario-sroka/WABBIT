!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name module_reactive_navier_stokes.f90
!> \version 0.5
!> \author msr
!
!> \brief module for combustion in 2D/3D reactive navier stokes equations:
!> structure: 
!> - include all subroutines from separat files
!> - public interface subroutine is called from wabbit-core and works as a 
!>   simple wrapper reactive navier-stokes subroutines
!> - physics specific params must be used as global data, because they have removed 
!>   from wabbit-core, same for cantera gas structure
!>
!!
!!
!! = log ======================================================================================
!! \n
!! 09/10/18 - create
!
! ********************************************************************************************

module module_reactive_navier_stokes

!---------------------------------------------------------------------------------------------
! modules

    use module_precision
    use module_ini_files_parser_mpi
    use module_reactive_navier_stokes_params
    use module_operators
#ifdef CANTERA
    use cantera
#endif

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !**********************************************************************************************
    ! make everything private if not explicitly marked public
    PRIVATE

    !**********************************************************************************************
    ! only this functions are visible outside this module
    PUBLIC :: interface_reactive_navier_stokes

    !**********************************************************************************************

!---------------------------------------------------------------------------------------------
! main body

    contains

    !----- IO -----------------------------------------
    include "IO/read_parameter_wabbit_core.f90"
    include "IO/field_names_reactive_ns.f90"
    include "IO/read_parameter_combustion.f90"
    include "IO/read_parameter_chemistry.f90"
    include "IO/prepare_saved_data_reactive_ns.f90"
    include "IO/convert_to_primitive.f90"
    include "IO/convert_from_primitive.f90"
    include "IO/compute_temperature.f90"
    include "IO/compute_reaction_rate.f90"
    !--------------------------------------------------

    !----- INI ----------------------------------------
    include "INI/ini_reactive_ns.f90"
    include "INI/inicond_zero_velocity.f90"
    include "INI/inicond_spark.f90"
    include "INI/inicond_blob.f90"
    include "INI/inicond_taylor_green.f90"
    !--------------------------------------------------

    !----- FILTER -------------------------------------
    include "FILTER/filter_reactive_ns.f90"
    include "FILTER/explicit_filter.f90"
    include "FILTER/filter_1D.f90"
    !--------------------------------------------------

    !----- RHS ----------------------------------------
    include "RHS/RHS_wrapper_reactive_ns.f90"
    include "RHS/RHS_3D_CANTERA_navier_stokes_reactive_periodicBC.f90"
    include "RHS/RHS_diff_subroutines.f90"
    include "RHS/diff_wrapper_3D.f90"
    include "RHS/RHS_3D_navier_stokes_non_reactive_periodicBC.f90"
    !--------------------------------------------------

    !----- TIME ---------------------------------------
    include "TIME/get_dt_reactive_ns.f90"
    !--------------------------------------------------

    !----- STATISTICS ---------------------------------
    include "STATISTICS/statistics_reactive_ns.f90"
    include "STATISTICS/compute_dilatational_dissipation.f90"
    include "STATISTICS/compute_solenoidal_dissipation.f90"
    include "STATISTICS/compute_DFT.f90"
    include "STATISTICS/compute_IDFT.f90"
    !--------------------------------------------------

    ! ********************************************************************************************
    ! interface
    ! note: parameter struct is definded and saved here, and all subroutine calls start here
    ! ********************************************************************************************
    subroutine interface_reactive_navier_stokes( interface_switch, filename, phi, phi_work, x0, dx, dF, dFname, rhs, stage, time, dt, filter_type, lgt_n )

    !---------------------------------------------------------------------------------------------
    ! variables

    implicit none

    !> type of specific interface: ini, read, ...
    character(len=*), intent(in)                    :: interface_switch

    ! optional parameter:
    ! -------------------
    !> filename
    character(len=*), intent(in), optional          :: filename
    !> statevector and work array
    real(kind=rk), intent(inout), optional          :: phi(:, :, :, :), phi_work(:, :, :, :)
    !> spacing and origin of block
    real(kind=rk), intent(in), optional             :: x0(1:3), dx(1:3)
    ! field index
    integer(kind=ik), intent(in), optional          :: dF
    ! returns the name of the field, corresponds to datafield index from ini file
    character(len=80), intent(out), optional        :: dFname
    !> rhs output
    real(kind=rk), intent(inout), optional          :: rhs(:, :, :, :)
    ! stage. there is 3 stages, init_stage, integral_stage and local_stage. If the PDE has
    ! terms that depend on global qtys, such as forces etc, which cannot be computed
    ! from a single block alone, the first stage does that. the second stage can then
    ! use these integral qtys for the actual RHS evaluation.
    character(len=*), intent(in), optional          :: stage
    ! it may happen that some source terms have an explicit time-dependency
    ! therefore the general call has to pass time
    real(kind=rk), intent (in), optional            :: time
    ! calculated time step
    real(kind=rk), intent (inout), optional         :: dt
    ! filter type
    character(len=*), intent(inout), optional       :: filter_type
    !> number of active blocks (light data)
    integer(kind=ik), intent(inout), optional       :: lgt_n

    ! physics datatype:
    ! note: use save to keep values in memory, when leave this subroutine to wabbit
    ! -----------------------------------------------------------------------------
    type(type_params_rns), save                     :: params_physics

    ! Cantera gas datatype:
    ! note: use save to keep values in memory, when leave this subroutine to wabbit
    ! -----------------------------------------------------------------------------
    type(phase_t), save                             :: gas

    !---------------------------------------------------------------------------------------------
    ! main body

        select case (interface_switch)

            ! ---------------------------------------
            ! read ini file and set all other params
            ! ---------------------------------------
            case('read_ini_file')
                call read_parameter_wabbit_core( params_physics, filename )
                call read_parameter_combustion( params_physics, params_physics%reactive_ns_file, gas )
                call read_parameter_chemistry( params_physics, params_physics%chemistry_file, gas )

            ! ---------------------------------------
            ! initialize data fields
            ! ---------------------------------------
            case('initialize_data')
                call ini_reactive_ns( params_physics, phi, x0, dx, gas )

            ! ---------------------------------------
            ! prepare save data
            ! ---------------------------------------
            case('prepare_save_data')
                call prepare_saved_data_reactive_ns( params_physics, phi, phi_work, time, dx, x0, gas )

            ! ---------------------------------------
            ! field names
            ! ---------------------------------------
            case('field_names')
                call field_names_reactive_ns( params_physics, dF, dFname )

            ! ---------------------------------------
            ! RHS
            ! ---------------------------------------
            case('rhs')
                call RHS_wrapper_reactive_ns( params_physics, time, phi, phi_work, x0, dx, rhs, stage, gas )

            ! ---------------------------------------
            ! filter
            ! ---------------------------------------
            case('filter_data')
                call filter_reactive_ns( params_physics, phi, phi_work, x0, dx, params_physics%filter_type )

            ! ---------------------------------------
            ! statistics
            ! ---------------------------------------
            case('statistics')
                call statistics_reactive_ns( params_physics, time, phi, phi_work, x0, dx, stage, lgt_n, gas )

            ! ---------------------------------------
            ! time step
            ! ---------------------------------------
            case('get_dt')
                call get_dt_reactive_ns( params_physics, time, phi, x0, dx, dt )

            ! ---------------------------------------
            ! error case
            ! ---------------------------------------
            case default
                call abort(091018001,"ERROR: unknown switch type in reactive navier stokes physics module interface.")

        end select

    end subroutine interface_reactive_navier_stokes

end module module_reactive_navier_stokes
