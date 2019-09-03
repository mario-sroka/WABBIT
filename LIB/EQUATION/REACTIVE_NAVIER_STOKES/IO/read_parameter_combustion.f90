!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name read_parameter_combustion.f90
!> \version 0.5
!> \author msr
!
!> \brief read parameter for reactive navier stokes physics
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

subroutine read_parameter_combustion( params_physics, filename, gas )

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

    ! dummy variable
    integer(kind=ik)                        :: dummy(3), k
    real(kind=rk)                           :: dummy2(6)

!---------------------------------------------------------------------------------------------
! variables initialization

    ! open file
    call read_ini_file_mpi(FILE, filename, .true.)

!---------------------------------------------------------------------------------------------
! main body

    ! reactive navier stokes parameter
    !-----------------------------------------------------------------------------------------
    ! get inicond name
    call read_param_mpi(FILE, 'Combustion', 'inicond', params_physics%inicond_name, '---' )

    ! inicond scales and position, note: vector elements have different meaning depending on specific initial condition
    call read_param_mpi(FILE, 'Combustion', 'inicond_scales', params_physics%inicond_scales, params_physics%inicond_scales )
    call read_param_mpi(FILE, 'Combustion', 'inicond_position', params_physics%inicond_position, params_physics%inicond_position )
    params_physics%inicond_position = params_physics%inicond_position * params_physics%L

    ! constant inicond values
    call read_param_mpi(FILE, 'Combustion', 'inicond_p', params_physics%inicond_p, 0.0_rk )
    call read_param_mpi(FILE, 'Combustion', 'inicond_rho', params_physics%inicond_rho, 0.0_rk )

    ! read field indexes
    ! note: the given names are used for initializing and rhs computation,
    ! but the do not correspond to data io
    call read_param_mpi(FILE, 'Combustion', 'rhoF', params_physics%rhoF, -1 )
    call read_param_mpi(FILE, 'Combustion', 'UxF', params_physics%UxF, -1 )
    call read_param_mpi(FILE, 'Combustion', 'UyF', params_physics%UyF, -1 )
    call read_param_mpi(FILE, 'Combustion', 'UzF', params_physics%UzF, -1 )
    call read_param_mpi(FILE, 'Combustion', 'EF', params_physics%EF, -1 )

    ! IO parameter
    call read_param_mpi(FILE, 'Combustion', 'save_primitive', params_physics%save_primitive, .false. )

    ! forcing parameters
    call read_param_mpi(FILE, 'Combustion', 'forcing', params_physics%forcing, .false. )
    call read_param_mpi(FILE, 'Combustion', 'k0', params_physics%k0, 0.0_rk )
    call read_param_mpi(FILE, 'Combustion', 'kmax', params_physics%kmax, 1 )
    call read_param_mpi(FILE, 'Combustion', 'eps_s_target', params_physics%eps_s_target, 0.0_rk )
    call read_param_mpi(FILE, 'Combustion', 'target_force', params_physics%target_force, 0.0_rk )

    ! allocate fourier coeffcients array, only if forcing enabled
    if ( params_physics%forcing ) then
        allocate( params_physics%phi_hat( params_physics%kmax, params_physics%kmax, params_physics%kmax, 3 ) )
        allocate( params_physics%phi_hat_d( params_physics%kmax, params_physics%kmax, params_physics%kmax, 3 ) )
        allocate( params_physics%phi_hat_s( params_physics%kmax, params_physics%kmax, params_physics%kmax, 3 ) )

        params_physics%phi_hat_d = 0.0_rk
        params_physics%phi_hat_s = 0.0_rk

        ! complex roots of unity
        ! Domainsizes
        dummy(1) = (params_physics%Bs(1)-1) * 2**params_physics%maxLvl
        dummy(2) = (params_physics%Bs(2)-1) * 2**params_physics%maxLvl
        dummy(3) = (params_physics%Bs(3)-1) * 2**params_physics%maxLvl

        ! allocate
        allocate( params_physics%rootsX( dummy(1) * params_physics%kmax  ) )
        allocate( params_physics%rootsY( dummy(2) * params_physics%kmax  ) )
        allocate( params_physics%rootsZ( dummy(3) * params_physics%kmax  ) )

        ! unity roots
        do k = 1, dummy(1) * params_physics%kmax
            params_physics%rootsX(k) = exp( cmplx( 0.0_rk, 2.0_rk*pi/real(dummy(1), kind=rk), kind=rk) ) &
                                     **real(k-1, kind=rk)
        end do
        do k = 1, dummy(2) * params_physics%kmax
            params_physics%rootsY(k) = exp( cmplx( 0.0_rk, 2.0_rk*pi/real(dummy(2), kind=rk), kind=rk) ) &
                                     **real(k-1, kind=rk)
        end do
        do k = 1, dummy(3) * params_physics%kmax
            params_physics%rootsZ(k) = exp( cmplx( 0.0_rk, 2.0_rk*pi/real(dummy(3), kind=rk), kind=rk) ) &
                                     **real(k-1, kind=rk)
        end do

    end if   

    ! sponge parameter
    call read_param_mpi(FILE, 'Combustion', 'sponge_width', params_physics%sponge_width, 0.0_rk )
    call read_param_mpi(FILE, 'Combustion', 'sponge_strength', params_physics%sponge_strength, 0.0_rk )

    ! sponge reference values
    dummy2 = 0.0_rk
    call read_param_mpi(FILE, 'Combustion', 'rho_ref', params_physics%rho_ref, dummy2 )
    call read_param_mpi(FILE, 'Combustion', 'u_ref', params_physics%rho_ref, dummy2 )
    call read_param_mpi(FILE, 'Combustion', 'v_ref', params_physics%rho_ref, dummy2 )
    call read_param_mpi(FILE, 'Combustion', 'w_ref', params_physics%rho_ref, dummy2 )
    call read_param_mpi(FILE, 'Combustion', 'es_ref', params_physics%rho_ref, dummy2 )

    ! clean up
    call clean_ini_file_mpi(FILE)

end subroutine read_parameter_combustion
