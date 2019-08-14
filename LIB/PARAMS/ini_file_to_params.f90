!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name ini_file_to_params.f90
!> \version 0.5
!> \author msr
!
!> \brief distribute blocks at start => create light data array
!
!>
!! input:    - filename \n
!! output:   - filled parameter struct \n
!!
!!
!! = log ======================================================================================
!! \n
!! 25/01/17    - create \n
!! 29/01/17    - add filter parameter \n
!! 30/01/17    - add automatic memory management
!
! ********************************************************************************************

subroutine ini_file_to_params( params, filename )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(inout)               :: params

    character(len=*), intent(in)                   :: filename

    ! process rank
    integer(kind=ik)                                :: rank
    ! number of processes
    integer(kind=ik)                                :: number_procs
    ! inifile structure
    type(inifile)                                   :: FILE
    ! maximum memory available on all cpus
    real(kind=rk)                                   :: maxmem, mem_per_block, max_neighbors, nstages
    ! string read from command line call
    character(len=80)                               :: memstring
    !
    integer(kind=ik)                                :: d,i, Nblocks_Jmax, g, Neqn, Nrk
    integer(kind=ik), dimension(3)                  :: Bs

    ! dummy file name array
    character(len=80), allocatable                  :: dummy_name(:)

!---------------------------------------------------------------------------------------------
! variables initialization

    ! set MPI parameter
    rank         = params%rank
    number_procs = params%number_procs

!---------------------------------------------------------------------------------------------
! main body

    ! read the file, only process 0 should create output on screen
    call set_lattice_spacing_mpi(1.0d0)
    call read_ini_file_mpi(FILE, filename, .true.)

    ! which physics module is used? (note that the initialization of different parameters takes
    ! place in those modules, i.e., they are not read here.)
    call read_param_mpi(FILE, 'Physics', 'physics_type', params%physics_type, "---" )

    call ini_domain(params, FILE )
    call ini_blocks(params,FILE)
    call ini_time(params,FILE)


    !**************************************************************************
    ! read INITIAL CONDITION parameters

    ! which physics module is used? (note that the initialization of different parameters takes
    ! place in those modules, i.e., they are not read here.)
    call read_param_mpi(FILE, 'Physics', 'physics_type', params%physics_type, "---" )

    ! if the initial condition is read from file, it is handled by wabbit itself, i.e. not
    ! by the physics modules. the pyhsics modules cannot do this, because they just see 'blocks'
    ! and never the entire grid as such.
    call read_param_mpi(FILE, 'Physics', 'read_from_files', params%read_from_files, .false. )

    ! wabbit does need to know how many fiels are written to disk when saving is triggered.
    ! e.g. saving ux, uy and p would mean 3. The names of these files as well as their contents
    ! are defined by the physics modules.
    call read_param_mpi(FILE, 'Saving', 'N_fields_saved', params%N_fields_saved, 3 )

    if (params%read_from_files ) then
        ! read variable names
        allocate( params%input_files( params%n_eqn ) )

        ! alternatively to explicitly file names, you can start from a given input time
        ! if a input time >= 0 is read, then the file names are computed, not read
        ! so: if you want to use file names, set time to any value lower than 0    
        ! --------------------------------------------------------------------------------
        call read_param_mpi(FILE, 'Physics', 'input_time', params%input_time, -1.0_rk )

        if ( params%input_time < 0.0_rk ) then
            ! read file names
            params%input_files = "---"
            call read_param_mpi(FILE, 'Physics', 'input_files', params%input_files, params%input_files)
        
        else

            ! read file names to dummy array, assume first file names correspond to datafields 
            ! (as it would be in a restarted computation)
            allocate( dummy_name( params%N_fields_saved ) )
            dummy_name = '---'
            call read_param_mpi(FILE, 'Saving', 'field_names', dummy_name, dummy_name)
           
            ! compute file names, /todo: time factor should synchronize with file saving or read from ini file
            do i = 1, params%n_eqn
                write( params%input_files(i) ,'(a, "_", i12.12, ".h5")') trim(adjustl(dummy_name(i))), nint(params%input_time * 1.0e9_rk)
            end do

        end if

    end if

    !***************************************************************************
    ! read DISCRETIZATION parameters
    !
    ! discretization order
    call read_param_mpi(FILE, 'Discretization', 'order_discretization', params%order_discretization, "---" )
    ! order of predictor for refinement
    call read_param_mpi(FILE, 'Discretization', 'order_predictor', params%order_predictor, "---" )
    ! filter frequency
    call read_param_mpi(FILE, 'Discretization', 'filter_type', params%filter_type, "no_filter" )
    call read_param_mpi(FILE, 'Discretization', 'filter_only_maxlevel', params%filter_only_maxlevel, .false. )
    if (params%filter_type /= "no_filter") then
        call read_param_mpi(FILE, 'Discretization', 'filter_freq', params%filter_freq, -1 )
    endif

    !***************************************************************************
    ! read statistics parameters
    call read_param_mpi(FILE, 'Statistics', 'nsave_stats', params%nsave_stats, 99999999_ik )
    call read_param_mpi(FILE, 'Statistics', 'tsave_stats', params%tsave_stats, 9999999.9_rk )
    !> assume start at time 0.0 /todo change if start with reloaded data
    params%next_stats_time = 0.0_rk + params%tsave_stats

    !***************************************************************************
    ! WABBIT needs to know about the mask function (if penalization is used): does it contain
    ! a time-dependent-part (e.g. moving obstacles, time-dependent forcing)? does it contain
    ! a time-independent part (fixed walls, homogeneous forcing)? or both? WABBIT needs to know
    ! that since we try to create the time-independent mask function only once, but the time-dependent
    ! part of course in every time step.
    call read_param_mpi(FILE, 'VPM', 'penalization', params%penalization, .false.)
    call read_param_mpi(FILE, 'VPM', 'mask_time_dependent_part', params%mask_time_dependent_part, .true.)
    call read_param_mpi(FILE, 'VPM', 'mask_time_independent_part', params%mask_time_independent_part, .true.)


    !***************************************************************************
    ! read DEBUG parameters
    !
    ! unit test treecode flag
    call read_param_mpi(FILE, 'Debug', 'test_treecode', params%test_treecode, .false.)
    call read_param_mpi(FILE, 'Debug', 'test_ghost_nodes_synch', params%test_ghost_nodes_synch, .false.)
    call read_param_mpi(FILE, 'Debug', 'check_redundant_nodes', params%check_redundant_nodes, .false.)

    !***************************************************************************
    ! read MPI parameters
    !
    ! data exchange method
    call ini_MPI(params, FILE )

    ! clean up
    if (params%rank==0) write(*,'("INIT: cleaning ini file")')
    call clean_ini_file_mpi(FILE)


    ! check ghost nodes number
    if (params%rank==0) write(*,'("INIT: checking if g and predictor work together")')
    ! if ( (params%n_ghosts < 4) .and. (params%order_predictor == 'multiresolution_4th') ) then
    !     call abort("ERROR: need more ghost nodes for given refinement order")
    ! end if
    if ( (params%n_ghosts < 2) .and. (params%order_predictor == 'multiresolution_2nd') ) then
        call abort("ERROR: need more ghost nodes for given refinement order")
    end if
    if ( (params%n_ghosts < 2) .and. (params%order_discretization == 'FD_4th_central_optimized') ) then
        call abort("ERROR: need more ghost nodes for given derivative order")
    end if

end subroutine ini_file_to_params




!> @brief     reads parameters for initializing a bridge from file
  subroutine ini_MPI(params, FILE )
    implicit none
    !> pointer to inifile
    type(inifile) ,intent(inout)     :: FILE
    !> params structure of WABBIT
    type(type_params),intent(inout)  :: params

    if (params%rank==0) then
      write(*,*)
      write(*,*)
      write(*,*) "PARAMS: MPI Communication!"
      write(*,'(" ----------------------------")')
    endif

    ! READ Bridge Parameter
    ! ----------------------
    ! decide if we need a bridge
    call read_param_mpi(FILE, 'BRIDGE', 'connect_with_bridge', params%bridge_exists, .false.)
    if (params%bridge_exists) then               ! if a bridge structure is required
      call read_param_mpi(FILE, 'BRIDGE', 'bridgeCommonMPI', params%bridgeCommonMPI, .false. )
      call read_param_mpi(FILE, 'BRIDGE', 'bridgeFluidMaster', params%bridgeFluidMaster, .false. )
      call read_param_mpi(FILE, 'BRIDGE', 'particleCommand', params%particleCommand, "---" )
    endif

  end subroutine ini_MPI

!-------------------------------------------------------------------------!!!!


!> @brief     reads parameters for initializing grid parameters
  subroutine ini_domain(params, FILE )
    implicit none
    !> pointer to inifile
    type(inifile) ,intent(inout)     :: FILE
    !> params structure of WABBIT
    type(type_params),intent(inout)  :: params


    if (params%rank==0) then
      write(*,*)
      write(*,*)
      write(*,*) "PARAMS: Domain"
      write(*,'(" -----------------")')
    endif


    call read_param_mpi(FILE, 'Domain', 'dim', params%dim, 2 )
    if ( .not.(params%dim==2 .or. params%dim==3) ) then
         call abort(234534,"Hawking: ERROR! The idea of 10 dimensions might sound exciting, &
         & but they would cause real problems if you forget where you parked your car. Tip: &
         & Try dim=2 or dim=3 ")
    endif

    params%domain_size=(/ 1.0_rk, 1.0_rk, 0.0_rk /) !default
    call read_param_mpi(FILE, 'Domain', 'domain_size', params%domain_size(1:params%dim), &
                                                       params%domain_size(1:params%dim) )

    params%periodic_BC = .true.
    call read_param_mpi(FILE, 'Domain', 'periodic_BC', params%periodic_BC(1:params%dim), &
                                                       params%periodic_BC(1:params%dim) )
  end subroutine ini_domain


!> @brief     reads parameters for initializing grid parameters
  subroutine ini_blocks(params, FILE )
    implicit none
    !> pointer to inifile
    type(inifile) ,intent(inout)     :: FILE
    !> params structure of WABBIT
    type(type_params),intent(inout)  :: params
    !> power used for dimensionality (d=2 or d=3)
    integer(kind=ik) :: i
    real(kind=rk), dimension(:), allocatable  :: tmp
    integer(kind=ik), dimension(:), allocatable  :: tmp_int
    if (params%rank==0) then
      write(*,*)
      write(*,*)
      write(*,*) "PARAMS: Blocks"
      write(*,'(" -----------------")')
    endif

    ! read number_block_nodes
    params%Bs =read_Bs(FILE, 'Blocks', 'number_block_nodes', params%Bs,params%dim)

    call read_param_mpi(FILE, 'Blocks', 'max_forest_size', params%forest_size, 3 )
    call read_param_mpi(FILE, 'Blocks', 'number_ghost_nodes', params%n_ghosts, 1 )
    call read_param_mpi(FILE, 'Blocks', 'number_blocks', params%number_blocks, -1 )
    call read_param_mpi(FILE, 'Blocks', 'number_equations', params%n_eqn, 1 )
    call read_param_mpi(FILE, 'Blocks', 'eps', params%eps, 1e-3_rk )
    call read_param_mpi(FILE, 'Blocks', 'eps_normalized', params%eps_normalized, .false. )
    call read_param_mpi(FILE, 'Blocks', 'max_treelevel', params%max_treelevel, 5 )
    call read_param_mpi(FILE, 'Blocks', 'min_treelevel', params%min_treelevel, 1 )

    if ( params%max_treelevel < params%min_treelevel ) then
        call abort(2609181,"Error: Minimal Treelevel cant be larger then Max Treelevel! ")
    end if

    if ( params%max_treelevel > 18 ) then
        ! as we internally convert the treecode to a single integer number, the number of digits is
        ! limited by that type. The largest 64-bit integer is 9 223 372 036 854 775 807
        ! which is 19 digits, but the 18th digit cannot be arbitrarily set. Therefore, 18 refinement levels
        ! are the maximum this code can currently perform.
        call abort(170619,"Error: Max treelevel cannot be larger 18 (64bit long integer problem) ")
    end if

    ! read switch to turn on|off mesh refinement
    call read_param_mpi(FILE, 'Blocks', 'adapt_mesh', params%adapt_mesh, .true. )
    call read_param_mpi(FILE, 'Blocks', 'adapt_inicond', params%adapt_inicond, params%adapt_mesh )
    call read_param_mpi(FILE, 'Blocks', 'inicond_refinements', params%inicond_refinements, 0 )
    call read_param_mpi(FILE, 'Blocks', 'block_dist', params%block_distribution, "---" )
    call read_param_mpi(FILE, 'Blocks', 'coarsening_indicator', params%coarsening_indicator, "threshold-state-vector" )
    call read_param_mpi(FILE, 'Blocks', 'threshold_mask', params%threshold_mask, .false. )
    call read_param_mpi(FILE, 'Blocks', 'force_maxlevel_dealiasing', params%force_maxlevel_dealiasing, .false. )
    call read_param_mpi(FILE, 'Blocks', 'N_dt_per_grid', params%N_dt_per_grid, 1_ik )

    ! Which components of the state vector (if indicator is "threshold-state-vector") shall we
    select case(params%coarsening_indicator)
    
        case('threshold-state-vector')
            ! use? in ACM, it can be good NOT to apply it to the pressure.
            allocate(tmp(1:params%n_eqn))
            allocate(params%threshold_state_vector_component(1:params%n_eqn))
            ! as default, use ones (all components used for indicator)
            tmp = 1.0_rk
            call read_param_mpi(FILE, 'Blocks', 'threshold_state_vector_component',  tmp, tmp )
            do i = 1, params%n_eqn
                 if (tmp(i)>0.0_rk) then
                     params%threshold_state_vector_component(i) = .true.
                 else
                     params%threshold_state_vector_component(i) = .false.
                 endif
            enddo
            deallocate(tmp)

         case('threshold-state-vector-rns')
             ! threshold state vector for reactive navier stokes: need list of state vector indexes

             ! read number of components for thresholding
             call read_param_mpi(FILE, 'Blocks', 'number_of_threshold_state_vector_components', params%number_of_threshold_state_vector_components, params%n_eqn )

             ! allocate dummy fields
             allocate(tmp_int(1:params%number_of_threshold_state_vector_components))
             allocate(params%threshold_state_vector_component(1:params%number_of_threshold_state_vector_components))

             ! default values are .false.
             params%threshold_state_vector_component = .false.

             ! read values from ini file
             call read_param_mpi(FILE, 'Blocks', 'threshold_state_vector_component', tmp_int, tmp_int )
             do i = 1, params%number_of_threshold_state_vector_components
                 params%threshold_state_vector_component(tmp_int(i)) = .true.
             end do
             
             ! clean up
             deallocate(tmp_int)

     end select

    ! number of maximal components for ghost nodes synchronizing
    ! note: do not change default value!
    call read_param_mpi(FILE, 'Blocks', 'N_max_components', params%N_max_components, 6 )

  end subroutine ini_blocks

  !-------------------------------------------------------------------------!!!!


  !> @brief     reads parameters for time stepping
    subroutine ini_time(params, FILE )
      implicit none
      !> pointer to inifile
      type(inifile) ,intent(inout)     :: FILE
      !> params structure of WABBIT
      type(type_params),intent(inout)  :: params
      real(kind=rk) :: butcher_RK4(1:5,1:5)

      if (params%rank==0) then
        write(*,*)
        write(*,*)
        write(*,*) "PARAMS: Time"
        write(*,'(" --------------")')
      endif

      ! time to reach in simulation
      call read_param_mpi(FILE, 'Time', 'time_max', params%time_max, 1.0_rk )
      ! maximum walltime before ending job
      call read_param_mpi(FILE, 'Time', 'walltime_max', params%walltime_max, 24.0_rk*7-0_rk )
      ! number of time steps to be performed. default value is very large, so if not set
      ! the limit will not be reached
      call read_param_mpi(FILE, 'Time', 'nt', params%nt, 99999999_ik )
      call read_param_mpi(FILE, 'Time', 'time_step_method', params%time_step_method, "RungeKuttaGeneric" )
      call read_param_mpi(FILE, 'Time', 'M_krylov', params%M_krylov, 12 )
      call read_param_mpi(FILE, 'Time', 'krylov_err_threshold', params%krylov_err_threshold, 1.0e-3_rk )
      call read_param_mpi(FILE, 'Time', 'krylov_subspace_dimension', params%krylov_subspace_dimension, "fixed" )
      ! read output write method
      call read_param_mpi(FILE, 'Time', 'write_method', params%write_method, "fixed_freq" )
      ! read output write frequency
      call read_param_mpi(FILE, 'Time', 'write_freq', params%write_freq, 25 )
      ! read output write frequency
      call read_param_mpi(FILE, 'Time', 'write_time', params%write_time, 1.0_rk )
      ! assume start at time 0.0
      params%next_write_time = 0.0_rk + params%write_time
      ! read value of fixed time step
      call read_param_mpi(FILE, 'Time', 'dt_fixed', params%dt_fixed, 0.0_rk )
      ! read value of fixed time step
      call read_param_mpi(FILE, 'Time', 'dt_max', params%dt_max, 0.0_rk )
      ! read CFL number
      call read_param_mpi(FILE, 'Time', 'CFL', params%CFL, 0.5_rk )

      ! read butcher tableau (set default value to RK4)
      butcher_RK4(1,1:5) = (/0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk, 0.0_rk/)
      butcher_RK4(2,1:5) = (/0.5_rk, 0.5_rk, 0.0_rk, 0.0_rk, 0.0_rk/)
      butcher_RK4(3,1:5) = (/0.5_rk, 0.0_rk, 0.5_rk, 0.0_rk, 0.0_rk/)
      butcher_RK4(4,1:5) = (/1.0_rk, 0.0_rk, 0.0_rk, 1.0_rk, 0.0_rk/)
      butcher_RK4(5,1:5) = (/0.0_rk, 1.0_rk/6.0_rk, 1.0_rk/3.0_rk, 1.0_rk/3.0_rk, 1.0_rk/6.0_rk/)

      call read_param_mpi(FILE, 'Time', 'butcher_tableau', params%butcher_tableau, butcher_RK4)
    end subroutine ini_time

    !-------------------------------------------------------------------------!
    !> @brief remove (multiple) blancs as seperators in a string
    subroutine merge_blancs(string_merge)
      ! this routine removes blanks at the beginning and end of an string
      ! and multiple blanks which are right next to each other

      implicit none
      character(len=*), intent(inout) :: string_merge
      integer(kind=ik) :: i, j, len_str, count

      len_str = len(string_merge)
      count = 0

      string_merge = string_merge
      do i=1,len_str-1
        if (string_merge(i:i)==" " .and. string_merge(i+1:i+1)==" ") then
          count = count + 1
          string_merge(i+1:len_str-1) = string_merge(i+2:len_str)
        end if
      end do

      string_merge = adjustl(string_merge)

    end subroutine merge_blancs


    !-------------------------------------------------------------------------!
    !> @brief Read Bs from inifile for unknown number of Bs in inifile
function read_Bs(FILE, section, keyword, default_Bs, dims) result(Bs)
  type(inifile) ,intent(inout)     :: FILE
  character(len=*), intent(in)    :: section ! What section do you look for? for example [Resolution]
  character(len=*), intent(in)    :: keyword ! what keyword do you
  integer(kind=ik), intent(in)    :: default_Bs(:)
  integer(kind=ik), intent(in)    :: dims !number of dimensions
  character(len=:), allocatable :: output_trim
  integer(kind=ik):: Bs(3)
  integer(kind=ik):: i, n_entries
  character(len=80):: output

  Bs = 1
  ! read number_block_nodes
  call read_param_mpi(FILE, section, keyword, output, "empty")
  if (trim(output) .eq. "empty") then
      write(*,'("Warning!! ", A, "[",A,"] is empty! Using default! ")') keyword, section
      Bs=default_Bs
  else
    call merge_blancs(output)
    output_trim=trim(output)
    call count_entries(output_trim, " ", n_entries)
    ! check if the number of entries is valid
    if (n_entries > dims) call abort(10519,"Dimensions and number of Bs entries dissagree!")
    ! Cast the output string into the integer
    read(output_trim,*) (Bs(i), i=1,n_entries)
    ! If only one Bs is given in the ini file, we duplicate it
    ! for the rest of the Bs array:
    if (n_entries==1) then
      Bs(1:dims) = Bs(1)
    endif
  endif
end function

    !-------------------------------------------------------------------------!
    !> @brief count number of vector elements in a string
    subroutine count_entries(string_cnt, seperator, n_entries)
      ! only to be used after merged blaks
      ! this routine counts the seperators and gives back this value +1

      implicit none
      character(len=1), intent(in) :: seperator
      character(len=*), intent(in) :: string_cnt
      integer(kind=ik), intent(out) :: n_entries
      integer(kind=ik) :: count_seperator, i, l_string
      character(len=Len_trim(adjustl(string_cnt))):: string_trim

      count_seperator = 0
      string_trim = trim(adjustl(string_cnt))
      l_string = LEN(string_trim)

      do i=1,l_string
        if (string_trim(i:i)==seperator) then
          count_seperator = count_seperator + 1
        end if
      end do

      n_entries = count_seperator + 1

    end subroutine count_entries
