!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name read_parameter_wabbit_core.f90
!> \version 0.5
!> \author msr
!
!> \brief read parameter from wabbit core ini file
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

subroutine read_parameter_wabbit_core( params_physics, filename )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> parameter struct
    type(type_params_rns), intent(inout)    :: params_physics
    !> filename
    character(len=*), intent(in)            :: filename

    ! inifile structure
    type(inifile)                           :: FILE

    ! MPI error variable
    integer(kind=ik)                        :: ierr

    ! loop variable
    integer(kind=ik)                        :: i

    ! dummy variable
    logical                                 :: dummy(3)

!---------------------------------------------------------------------------------------------
! variables initialization

    ! open file
    call read_ini_file_mpi(FILE, filename, .true.)

!---------------------------------------------------------------------------------------------
! main body

    ! WABBIT parameter
    !-----------------------------------------------------------------------------------------
    ! dimension
    call read_param_mpi(FILE, 'Domain', 'dim', params_physics%d, 2 )

    ! read number_block_nodes
    if ( params_physics%d == 2 ) then
        call read_param_mpi(FILE, 'Blocks', 'number_block_nodes', params_physics%Bs(1:2), (/0, 0/) )
    else
        call read_param_mpi(FILE, 'Blocks', 'number_block_nodes', params_physics%Bs(1:3), (/0, 0, 0/) )
    end if
    ! read number_ghost_nodes
    call read_param_mpi(FILE, 'Blocks', 'number_ghost_nodes', params_physics%g, 1 )
    ! boundaries
    params_physics%periodic_BC = .true.
    call read_param_mpi(FILE, 'Domain', 'periodic_BC', params_physics%periodic_BC(1:params_physics%d), &
                                                       params_physics%periodic_BC(1:params_physics%d) )

    ! domain size
    params_physics%L = (/ 1.0_rk, 1.0_rk, 0.0_rk /) !default
    call read_param_mpi(FILE, 'Domain', 'domain_size', params_physics%L(1:params_physics%d), &
                                                       params_physics%L(1:params_physics%d) )

    ! number of datafields
    call read_param_mpi(FILE, 'Blocks', 'number_equations', params_physics%NdF, 1 )

    ! max block level
    call read_param_mpi(FILE, 'Blocks', 'max_treelevel', params_physics%maxLvl, 1 )

    ! filter parameter
    call read_param_mpi(FILE, 'Discretization', 'filter_type', params_physics%filter_type, '---' )
    call read_param_mpi(FILE, 'Discretization', 'r_th', params_physics%r_th, 0.0_rk )
    call read_param_mpi(FILE, 'Discretization', 'bogey_detector', params_physics%bogey_detector, '---' )
    call read_param_mpi(FILE, 'Discretization', 'bogey_sigma', params_physics%bogey_sigma, '---' )

    ! init filter
    select case(params_physics%filter_type)
        case('explicit_3pt')
            allocate( params_physics%filter_stencil(3) )
            params_physics%filter_stencil = (/  1.0_rk/  4.0_rk, &
                                               -1.0_rk/  2.0_rk, &
                                                1.0_rk/  4.0_rk/)
        case('explicit_5pt')
            allocate( params_physics%filter_stencil(5) )
            params_physics%filter_stencil = (/ -1.0_rk/ 16.0_rk, &
                                                1.0_rk/  4.0_rk, &
                                               -3.0_rk/  8.0_rk, &
                                                1.0_rk/  4.0_rk, &
                                               -1.0_rk/ 16.0_rk/)
        case('explicit_7pt')
            allocate( params_physics%filter_stencil(7) )
            params_physics%filter_stencil = (/  1.0_rk/ 64.0_rk, &
                                               -3.0_rk/ 32.0_rk, &
                                               15.0_rk/ 64.0_rk, &
                                               -5.0_rk/ 16.0_rk, &
                                               15.0_rk/ 64.0_rk, &
                                               -3.0_rk/ 32.0_rk, &
                                                1.0_rk/ 64.0_rk/)
        case('explicit_9pt')
            allocate( params_physics%filter_stencil(9) )
            params_physics%filter_stencil = (/ -1.0_rk/256.0_rk, &
                                                1.0_rk/ 32.0_rk, &
                                               -7.0_rk/ 64.0_rk, &
                                                7.0_rk/ 32.0_rk, &
                                              -35.0_rk/128.0_rk, &
                                                7.0_rk/ 32.0_rk, &
                                               -7.0_rk/ 64.0_rk, &
                                                1.0_rk/ 32.0_rk, &
                                               -1.0_rk/256.0_rk/)
        case('explicit_11pt')
            allocate( params_physics%filter_stencil(11) )
            params_physics%filter_stencil = (/  1.0_rk/1024.0_rk, &
                                               -5.0_rk/ 512.0_rk, &
                                               45.0_rk/1024.0_rk, &
                                              -15.0_rk/ 128.0_rk, &
                                              105.0_rk/ 512.0_rk, &
                                              -63.0_rk/ 256.0_rk, &
                                              105.0_rk/ 512.0_rk, &
                                              -15.0_rk/ 128.0_rk, &
                                              45.0_rk/1024.0_rk, &
                                               -5.0_rk/ 512.0_rk, &
                                                1.0_rk/1024.0_rk/)
        case('bogey_shock', 'wavelet', 'no_filter')
            ! nothing to do
        case default
                call abort(151018005,"ERROR in read parameter: filter type is unknown!")
    end select

    ! saving
    call read_param_mpi(FILE, 'Saving', 'N_fields_saved', params_physics%N_fields_saved, 1 )
    allocate( params_physics%names_saved( params_physics%N_fields_saved ) )
    params_physics%names_saved = "---"
    call read_param_mpi(FILE, 'Saving', 'field_names', params_physics%names_saved, params_physics%names_saved)

    ! MPI
    call MPI_Comm_rank (WABBIT_COMM, params_physics%rank, ierr)

    ! init
    call read_param_mpi(FILE, 'Physics', 'read_from_files', params_physics%read_from_files, .false. )
    ! read variable names
    allocate( params_physics%input_files( params_physics%NdF ) )

    ! alternatively to explicitly file names, you can start from a given input time
    ! if a input time >= 0 is read, then the file names are computed, not read
    ! so: if you want to use file names, set time to any value lower than 0    
    ! --------------------------------------------------------------------------------
    call read_param_mpi(FILE, 'Physics', 'input_time', params_physics%input_time, -1.0_rk )

    if ( params_physics%input_time < 0.0_rk ) then
        ! read file names
        params_physics%input_files = "---"
        call read_param_mpi(FILE, 'Physics', 'input_files', params_physics%input_files, params_physics%input_files)
        
    else

        ! use saving file names, assume first file names correspond to datafields 
        ! (as it would be in a restarted computation)
        
        ! compute file names, /todo: time factor should synchronize with file saving or read from ini file
        do i = 1, params_physics%NdF
            write( params_physics%input_files(i) ,'(a, "_", i12.12, ".h5")') trim(adjustl(params_physics%names_saved(i))), nint(params_physics%input_time * 1.0e9_rk)
        end do

    end if

    ! time stepping
    call read_param_mpi(FILE, 'Time', 'CFL', params_physics%CFL, 1.0_rk )

    ! read names for reactive navier stokes ini files from wabbit core ini file
    call read_param_mpi(FILE, 'Physics', 'reactive_parameter_file', params_physics%reactive_ns_file, '---')
    call read_param_mpi(FILE, 'Physics', 'chemistry_parameter_file', params_physics%chemistry_file, '---')

    ! clean up
    call clean_ini_file_mpi(FILE)

end subroutine read_parameter_wabbit_core
