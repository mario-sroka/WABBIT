!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name read_parameter_chemistry.f90
!> \version 0.5
!> \author msr
!
!> \brief read chemistry parameter from additional file
!>
!! input:    - params_physics, filename \n
!! output:   - \n
!!
!!
!! = log ======================================================================================
!! \n
!! 09/10/18 - create
!
! ********************************************************************************************

subroutine read_parameter_chemistry( params_physics, filename, gas )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> parameter struct
    type(type_params_rns), intent(inout)    :: params_physics
    !> filename
    character(len=*), intent(in)            :: filename
    !> Cantera gas mixture struct
    type(phase_t), intent(inout)            :: gas

    ! inifile structure
    type(inifile)                           :: FILE

    ! loop variable, dummy
    integer(kind=ik)                        :: k
    real(kind=rk)                           :: dummy

!---------------------------------------------------------------------------------------------
! variables initialization

    ! open file
    call read_ini_file_mpi(FILE, filename, .true.)

!---------------------------------------------------------------------------------------------
! main body

    !-----------------------------------------------------------------------------------------
    ! INI
    !-----------------------------------------------------------------------------------------
    ! read field indexes
    ! note: YF marks first species datafield, we assume all other species follows after  
    call read_param_mpi(FILE, 'INI', 'YF', params_physics%YF, -1 )

    ! number of species
    call read_param_mpi(FILE, 'INI', 'species', params_physics%species, 1 )
    ! additionally ini conditions
    call read_param_mpi(FILE, 'INI', 'T0', params_physics%T0, 0.0_rk )
    call read_param_mpi(FILE, 'INI', 'inicond_T', params_physics%inicond_T, 0.0_rk )

    allocate( params_physics%inicond_Y(params_physics%species) )
    call read_param_mpi(FILE, 'INI', 'inicond_Y', params_physics%inicond_Y )
    allocate( params_physics%inicond_X(params_physics%species) )

    !-----------------------------------------------------------------------------------------
    ! FLUID
    !-----------------------------------------------------------------------------------------
    ! read fluid properties
!    call read_param_mpi(FILE, 'FLUID', 'gamma_', params_physics%gamma_, 0.0_rk )
!    call read_param_mpi(FILE, 'FLUID', 'Rs', params_physics%Rs, 0.0_rk )
!    params_physics%Cv = params_physics%Rs/(params_physics%gamma_-1.0_rk)
!    params_physics%Cp = params_physics%Cv*params_physics%gamma_
!    call read_param_mpi(FILE, 'FLUID', 'Pr', params_physics%Pr, 0.0_rk )
!    call read_param_mpi(FILE, 'FLUID', 'mu0', params_physics%mu0, 0.0_rk )
    call read_param_mpi(FILE, 'FLUID', 'dissipation', params_physics%dissipation, .true. )
!    call read_param_mpi(FILE, 'FLUID', 'viscosity_model', params_physics%viscosity_model, '---' )

    !-----------------------------------------------------------------------------------------
    ! CHEM
    !-----------------------------------------------------------------------------------------
    call read_param_mpi(FILE, 'CHEM', 'chemistry_model', params_physics%chemistry_model, '---' )

    allocate( params_physics%dh(params_physics%species) )
    ! /todo read dh values from file for non cantera chemistries
    allocate( params_physics%Wk(params_physics%species) )

    !-----------------------------------------------------------------------------------------
    ! CANTERA
    !-----------------------------------------------------------------------------------------
    ! name of kinetics file and mixture
    call read_param_mpi(FILE, 'INI', 'kinetics_file', params_physics%kinetics_file, '---' )
    call read_param_mpi(FILE, 'INI', 'mixture_scheme', params_physics%mixture_scheme, '---' )

    ! set kinetics for cantera chemistry
    if (params_physics%chemistry_model == 'cantera') then

        ! set gas mixture
        gas = importPhase( trim(params_physics%kinetics_file), trim(params_physics%mixture_scheme) )

        ! set gas constant
        params_physics%gas_constant = 8.3144621_rk

        ! get molecular weights from cantera
        call getMolecularWeights(gas, params_physics%Wk)

        ! mole fractions inicond
        dummy = 0.0_rk
        do k = 1, params_physics%species
            dummy = dummy + params_physics%inicond_Y(k)/params_physics%Wk(k)
        end do
        dummy = 1.0_rk/dummy
        do k = 1, params_physics%species
            params_physics%inicond_X(k) = dummy/params_physics%Wk(k) * params_physics%inicond_Y(k)
        end do

        ! get enthalpies of formation
        call setState_TPX(gas, params_physics%T0, params_physics%inicond_p, params_physics%inicond_X )
        call getEnthalpies_RT(gas, params_physics%dh)
        params_physics%dh = params_physics%dh * params_physics%gas_constant * params_physics%T0
        ! correction 
        do k = 1, params_physics%species
            params_physics%dh(k) = params_physics%dh(k) / (params_physics%Wk(k) * 1e-3_rk) 
        end do

    end if

    ! clean up
    call clean_ini_file_mpi(FILE)

end subroutine read_parameter_chemistry
