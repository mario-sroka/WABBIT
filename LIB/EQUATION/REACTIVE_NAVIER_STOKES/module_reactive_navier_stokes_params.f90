!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name module_reactive_navier_stokes_params.f90
!> \version 0.5
!> \author msr
!
!> \brief params struct for reactive navier stokes equation
!>
!!
!!
!! = log ======================================================================================
!! \n
!! 09/10/18 - create
!
! ********************************************************************************************

module module_reactive_navier_stokes_params

!---------------------------------------------------------------------------------------------
! modules

    use module_precision

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    type :: type_params_rns

        ! parameter from WABBIT core
        ! --------------------------
        ! grid parameter
        integer(kind=ik)    :: Bs(3)            = 0      ! number of block nodes
        integer(kind=ik)    :: g                = 0      ! number of ghost nodes
        logical             :: periodic_BC(3)            ! boundaries

        ! domain parameter
        integer(kind=ik)    :: d                = 2      ! dimension
        real(kind=rk)       :: L(3)             = 0.0_rk ! domain size

        ! equation paramter
        integer(kind=ik)    :: NdF              = 2      ! number of datafields

        ! level parameter
        integer(kind=ik)    :: maxLvl           = 1      ! max level

        ! filter parameter
        character(len=80)   :: filter_type      = '---'
        real(kind=rk)       :: r_th             = 0.0_rk ! filter threshold
        character(len=80)   :: bogey_detector   = '---'  ! detector method
        character(len=80)   :: bogey_sigma      = '---'  ! sigma computation method
        real(kind=rk), &
        allocatable         :: filter_stencil(:)         ! filter stencil

        ! saving
        integer(kind=ik)    :: N_fields_saved   = -1     ! number of field to save
        character(len=80), &
        allocatable         :: names_saved(:)            ! names of saved field

        ! MPI parameter
        integer(kind=ik)    :: rank                      ! MPI rank

        ! init
        logical             :: read_from_files           ! start from input files
        character(len=80), &
        allocatable         :: input_files(:)            ! files we want to read for inital cond.
        real(kind=rk)       :: input_time       = -1.0_rk! time we want to read for inital cond. 		

        ! time stepping
        real(kind=rk)       :: CFL                       ! CFL number

        ! parameter for reactive navier stokes physics
        ! --------------------------------------------
        ! parameter file
        character(len=80)   :: reactive_ns_file = '---'  ! parameter file name
        ! inicondition
        character(len=80)   :: inicond_name     = '---'
        real(kind=rk)       :: inicond_scales(2) &
                                                = 0.0_rk ! scales for different ini conditions
        real(kind=rk)       :: inicond_position(3) &
                                                = 0.0_rk ! position for different ini conditions

        real(kind=rk)       :: T0               = 0.0_rk ! reference temperature
        real(kind=rk)       :: inicond_rho      = 0.0_rk ! constant initial values, maybe not used by given inicond
        real(kind=rk)       :: inicond_p        = 0.0_rk
        real(kind=rk)       :: inicond_u(3)     = 0.0_rk
        real(kind=rk)       :: inicond_T        = 0.0_rk
        real(kind=rk),&
        allocatable         :: inicond_Y(:)              ! mass fractions iniconditions
        real(kind=rk),&
        allocatable         :: inicond_X(:)              ! mole fractions iniconditions

        integer(kind=ik)    :: rhoF             = -1, &
                               UxF              = -1, &
                               UyF              = -1, &
                               UzF              = -1, &
                               EF               = -1, &
                               YF               = -1     ! data field indexes

        ! equation parameter, set initial value to -1 to force errors from start
        real(kind=rk)       :: gamma_           = 0.0_rk ! adiabatic coefficient
        real(kind=rk)       :: Rs               = 0.0_rk ! specific gas constant
        real(kind=rk)       :: Cv               = 0.0_rk ! isochoric heat capacity
        real(kind=rk)       :: Cp               = 0.0_rk ! isobaric heat capacity
        real(kind=rk)       :: Pr               = 0.0_rk ! prandtl number
        real(kind=rk)       :: mu0              = 0.0_rk ! dynamic viscosity
        logical             :: dissipation      = .true. ! dissipation switch
        character(len=80)   :: viscosity_model  = '---'  ! viscosity model switch

        ! IO
        logical             :: save_primitive   = .false.! IO switch to save non skew symmetric data
        character(len=80)   :: data_file_name   = '---'  ! additionally data file

        ! chmemistry parameter
        integer(kind=ik)    :: species          = 1      ! number of species
        character(len=80)   :: chemistry_file   = '---'  ! parameter file name
        character(len=80)   :: chemistry_model  = '---'  ! chemistry model
        character(len=80)   :: kinetics_file    = '---'  ! kinetics file name
        character(len=80)   :: mixture_scheme   = '---'  ! gas mixture

        real(kind=rk)       :: A                = 0.0_rk ! preexponential constant
        real(kind=rk)       :: Ta               = 0.0_rk ! activation temperature
        real(kind=rk)       ,&
        allocatable         :: dh(:)                     ! reaction heat
        real(kind=rk)       ,&
        allocatable         :: Wk(:)                     ! molecular weights
        real(kind=rk)       :: Lew              = 1.0_rk ! lewis number
        real(kind=rk)       :: F                = 1.0_rk ! flame thickness
        real(kind=rk)       :: gas_constant     = 0.0_rk ! gas constant

        ! statistics parameter
        real(kind=rk)       :: e_tot_mean       = 0.0_rk ! total energy mean
        real(kind=rk)       :: k_mean           = 0.0_rk ! turbulent kinetic energy
        real(kind=rk)       :: k_d_mean         = 0.0_rk ! solenoidal, dilatational component
        real(kind=rk)       :: k_s_mean         = 0.0_rk
        real(kind=rk)       :: epsilon_mean     = 0.0_rk ! turbulent dissipation
        real(kind=rk)       :: epsilon_d_mean   = 0.0_rk ! solenoidal, dilatational component
        real(kind=rk)       :: epsilon_s_mean   = 0.0_rk
        real(kind=rk)       :: UU_mean          = 0.0_rk ! U square mean
        real(kind=rk)       :: e_kin_mean       = 0.0_rk ! kinetic energy mean
        real(kind=rk)       :: e_int_mean       = 0.0_rk ! internal energy mean

        real(kind=rk)       :: rho_mean         = 0.0_rk ! density mean
        real(kind=rk)       :: u_mean           = 0.0_rk ! velocity mean
        real(kind=rk)       :: v_mean           = 0.0_rk
        real(kind=rk)       :: w_mean           = 0.0_rk
        real(kind=rk)       :: u_mean_0         = 0.0_rk ! save mean values from last time step
        real(kind=rk)       :: v_mean_0         = 0.0_rk
        real(kind=rk)       :: w_mean_0         = 0.0_rk

        real(kind=rk)       :: mu_mean          = 0.0_rk ! viscosity mean value

        ! forcing parameter
        logical             :: forcing          = .false.! enable/disable forcing
        real(kind=rk)       :: k0               = 0.0_rk ! wavenumber limit
        integer(kind=ik)    :: kmax             = 1      ! loop limit

        complex(kind=rk)    ,&
        allocatable         :: phi_hat(:,:,:,:)          ! fourier coeffcients
        complex(kind=rk)    ,&
        allocatable         :: phi_hat_d(:,:,:,:)        ! dilatational fourier coeffcients
        complex(kind=rk)    ,&
        allocatable         :: phi_hat_s(:,:,:,:)        ! solenoidal fourier coeffcients
        complex(kind=rk)    ,&
        allocatable         :: rootsX(:), rootsY(:), &   ! complex roots, note: need 3 arrays, for
                               rootsZ(:)                 ! cases with different domain sizes                                                         

        real(kind=rk)       :: target_force     = 0.0_rk ! target values
        real(kind=rk)       :: eps_s_target     = 0.0_rk !

        ! sponge
        real(kind=rk)       :: rho_ref(6)
        real(kind=rk)       :: u_ref(6)
        real(kind=rk)       :: v_ref(6)
        real(kind=rk)       :: w_ref(6)
        real(kind=rk)       :: es_ref(6)
        real(kind=rk)       ,&
        allocatable         :: Y_ref(:,:)

    end type type_params_rns

!---------------------------------------------------------------------------------------------
! main body

    contains

end module module_reactive_navier_stokes_params
