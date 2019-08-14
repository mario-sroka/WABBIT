!----------------------------------------------------------------
!> Interface between WABBIT and Physics Modules
!> This module contains all functions which WABBIT provides to
!> implement your physics module
!> \details
!> \version 0.5
!> \author engels
!----------------------------------------------------------------
module module_physics_metamodule

    use module_globals
    ! at this point, you bind all physics modules into one metamodule, so in the rest
    ! of the code, we just load that. as all other physics modules, it provides some
    ! public routines, at which the corresponding actual physics modules are called
    use module_ConvDiff_new
    use module_acm
    use module_navier_stokes
    use module_reactive_navier_stokes

    implicit none

    ! I usually find it helpful to use the private keyword by itself initially, which specifies
    ! that everything within the module is private unless explicitly marked public.
    PRIVATE

    !**********************************************************************************************
    ! These are the important routines that are visible to WABBIT:
    !**********************************************************************************************
    PUBLIC :: READ_PARAMETERS_meta, PREPARE_SAVE_DATA_meta, RHS_meta, GET_DT_BLOCK_meta, &
    INICOND_meta, FIELD_NAMES_meta,&
    STATISTICS_meta, FILTER_meta, CREATE_MASK_meta
    !**********************************************************************************************

contains

    ! this is a public routine to create the mask function for the penalization term.
    ! not all physics modules use it
    subroutine CREATE_MASK_meta( physics, time, x0, dx, Bs, g, mask, stage )
        implicit none
        character(len=*), intent(in) :: physics
        real(kind=rk), intent(in)    :: x0(1:3), dx(1:3), time
        integer(kind=ik), intent(in) :: Bs(1:3), g
        !     stage == "time-dependent-part"
        !     stage == "time-independent-part"
        !     stage == "all-parts"
        character(len=*), intent(in) :: stage
        ! mask data. we can use different trees (4est module) to generate time-dependent/indenpedent
        ! mask functions separately. This makes the mask routines tree-level routines (and no longer
        ! block level) so the physics modules have to provide an interface to create the mask at a tree
        ! level. All parts of the mask shall be included: chi, boundary values, sponges.
        ! On input, the mask array is correctly filled. You cannot create the full mask here.
        real(kind=rk), dimension(:,:,:,:), intent(inout) :: mask

        mask = 0.0_rk

        select case ( physics )
        case ('ACM-new')
            if (size(mask,3)==1) then
                call create_mask_2D_ACM( time, x0, dx, Bs, g, mask(:,:,1,:), stage )
            else
                call create_mask_3D_ACM( time, x0, dx, Bs, g, mask, stage)
            endif

        case ('ConvDiff-new')
            ! not implemented yet (this module does not use penalization)

        case ('navier_stokes')
            ! not implemented yet

        case ('reactive_navier_stokes')
            ! /todo implement create mask for reactive navier stokes physics in physics metamodule

        case default
            call abort(1212,'unknown physics...say whaaat?')

        end select

    end subroutine


    !-----------------------------------------------------------------------------
    !> \brief main level wrapper routine to read parameters in the physics module. It reads
    !> from the same ini file as wabbit, and it reads all it has to know. note in physics modules
    !> the parameter struct for wabbit is not available.
    subroutine READ_PARAMETERS_meta( physics, filename, N_mask_components )
        implicit none
        character(len=*), intent(in) :: physics
        ! number of grid-dependent (and not time-dependend qtys) is decided by the physics modules
        integer(kind=ik), intent(out) :: N_mask_components
        character(len=*), intent(in) :: filename

        ! default is none (for navier-stokes and convection-diffusion)
        N_mask_components = 0

        select case ( physics )
        case ('ACM-new')
            call READ_PARAMETERS_ACM( filename, N_mask_components )

        case ('ConvDiff-new')
            call READ_PARAMETERS_convdiff( filename )

        case ('navier_stokes')
            call READ_PARAMETERS_NStokes( filename )

        case ('reactive_navier_stokes')
            call interface_reactive_navier_stokes("read_ini_file", filename=filename)

        case default
            call abort(1212,'unknown physics...say whaaat?')

        end select

    end subroutine


    !-----------------------------------------------------------------------------
    ! save data. Since you might want to save derived data, such as the vorticity,
    ! the divergence etc., which are not in your state vector, this routine has to
    ! copy and compute what you want to save to the work array.
    !
    ! In the main code, save_fields than saves the first N_fields_saved components of the
    ! work array to file.
    !
    ! NOTE that as we have way more work arrays than actual state variables (typically
    ! for a RK4 that would be >= 4*dim), you can compute a lot of stuff, if you want to.
    !-----------------------------------------------------------------------------
    subroutine PREPARE_SAVE_DATA_meta( physics, time, u, g, x0, dx, work, mask )
        implicit none
        character(len=*), intent(in) :: physics

        ! it may happen that some source terms have an explicit time-dependency
        ! therefore the general call has to pass time
        real(kind=rk), intent (in) :: time

        ! block data, containg the state vector. In general a 4D field (3 dims+components)
        ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
        real(kind=rk), intent(inout) :: u(1:,1:,1:,1:)

        ! as you are allowed to compute the RHS only in the interior of the field
        ! you also need to know where 'interior' starts: so we pass the number of ghost points
        integer, intent(in) :: g

        ! for each block, you'll need to know where it lies in physical space. The first
        ! non-ghost point has the coordinate x0, from then on its just cartesian with dx spacing
        real(kind=rk), intent(in) :: x0(1:3), dx(1:3)

        ! output in work array.
        real(kind=rk), intent(inout) :: work(1:,1:,1:,1:)
        ! mask data. we can use different trees (4est module) to generate time-dependent/indenpedent
        ! mask functions separately. This makes the mask routines tree-level routines (and no longer
        ! block level) so the physics modules have to provide an interface to create the mask at a tree
        ! level. All parts of the mask shall be included: chi, boundary values, sponges.
        ! On input, the mask array is correctly filled. You cannot create the full mask here.
        real(kind=rk), intent(inout) :: mask(1:,1:,1:,1:)

        select case(physics)
        case ('ACM-new')
            call PREPARE_SAVE_DATA_ACM( time, u, g, x0, dx, work, mask )

        case ('ConvDiff-new')
            call PREPARE_SAVE_DATA_convdiff( time, u, g, x0, dx, work )

        case ('navier_stokes')
            call PREPARE_SAVE_DATA_NStokes( time, u, g, x0, dx, work )

        case ('reactive_navier_stokes')
            call interface_reactive_navier_stokes("prepare_save_data", phi=u, phi_work=work, time=time, dx=dx, x0=x0)

        case default
            call abort(88119, "[PREPARE_SAVE_DATA (metamodule)] unknown physics....")

        end select

    end subroutine


    !-----------------------------------------------------------------------------
    ! when savig to disk, WABBIT would like to know how you named you variables.
    ! e.g. u(:,:,:,1) is called "ux"
    !
    ! the main routine save_fields has to know how you label the stuff you want to
    ! store from the work array, and this routine returns those strings
    !-----------------------------------------------------------------------------
    subroutine FIELD_NAMES_meta( physics, N, name )
        implicit none
        character(len=*), intent(in) :: physics
        ! component index
        integer(kind=ik), intent(in) :: N
        ! returns the name
        character(len=80), intent(out) :: name

        select case(physics)
        case ('ACM-new')
            call FIELD_NAMES_ACM(N, name)

        case ('ConvDiff-new')
            call FIELD_NAMES_convdiff(N, name)

        case ('navier_stokes')
            call FIELD_NAMES_NStokes(N, name)

        case ('reactive_navier_stokes')
            call interface_reactive_navier_stokes("field_names", dF=N, dFname=name)

        case default
            call abort(88119, "[FIELD_NAMES (metamodule):] unknown physics....")

        end select

    end subroutine FIELD_NAMES_meta


    !-----------------------------------------------------------------------------
    ! main level wrapper to set the right hand side on a block. Note this is completely
    ! independent of the grid and any MPI formalism, neighboring relations and the like.
    ! You just get a block data (e.g. ux, uy, uz, p) and compute the right hand side
    ! from that. Ghost nodes are assumed to be sync'ed.
    !-----------------------------------------------------------------------------
    subroutine RHS_meta( physics, time, u, g, x0, dx, rhs, mask, stage, boundary_flag)
        implicit none

        character(len=*), intent(in) :: physics
        ! it may happen that some source terms have an explicit time-dependency
        ! therefore the general call has to pass time
        real(kind=rk), intent (in) :: time

        ! block data, containg the state vector. In general a 4D field (3 dims+components)
        ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
        real(kind=rk), intent(inout) :: u(1:,1:,1:,1:)

        ! as you are allowed to compute the RHS only in the interior of the field
        ! you also need to know where 'interior' starts: so we pass the number of ghost points
        integer, intent(in) :: g

        ! for each block, you'll need to know where it lies in physical space. The first
        ! non-ghost point has the coordinate x0, from then on its just cartesian with dx spacing
        real(kind=rk), intent(in) :: x0(1:3), dx(1:3)

        ! output. Note assumed-shape arrays
        real(kind=rk), intent(inout) :: rhs(1:,1:,1:,1:)

        ! mask data. we can use different trees (4est module) to generate time-dependent/indenpedent
        ! mask functions separately. This makes the mask routines tree-level routines (and no longer
        ! block level) so the physics modules have to provide an interface to create the mask at a tree
        ! level. All parts of the mask shall be included: chi, boundary values, sponges.
        ! On input, the mask array is correctly filled. You cannot create the full mask here.
        real(kind=rk), intent(in) :: mask(1:,1:,1:,1:)


        ! stage. there is 3 stages, init_stage, integral_stage and local_stage. If the PDE has
        ! terms that depend on global qtys, such as forces etc, which cannot be computed
        ! from a single block alone, the first stage does that. the second stage can then
        ! use these integral qtys for the actual RHS evaluation.
        character(len=*), intent(in) :: stage

        ! when implementing boundary conditions, it is necessary to now if the local field (block)
        ! is adjacent to a boundary, because the stencil has to be modified on the domain boundary.
        ! The boundary_flag tells you if the local field is adjacent to a domain boundary:
        ! boundary_flag(i) can be either 0, 1, -1,
        !  0: no boundary in the direction +/-e_i
        !  1: boundary in the direction +e_i
        ! -1: boundary in the direction - e_i
        ! currently only acessible in the local stage
        integer(kind=2), optional, intent(in):: boundary_flag(3)

        select case(physics)
        case ("ACM-new")
            call RHS_ACM( time, u, g, x0, dx,  rhs, mask, stage )

        case ("ConvDiff-new")
            call RHS_convdiff( time, u, g, x0, dx, rhs, stage, boundary_flag )

        case ("navier_stokes")
            call RHS_NStokes( time, u, g, x0, dx, rhs, stage, boundary_flag )

        case ('reactive_navier_stokes')
            call interface_reactive_navier_stokes("rhs", time=time, phi=u, x0=x0, dx=dx, rhs=rhs, &
                 stage=stage)

        case default
            call abort(2152000, "[RHS_wrapper.f90]: physics_type is unknown"//physics)

        end select

    end subroutine RHS_meta


    !-----------------------------------------------------------------------------
    ! main level wrapper to compute statistics (such as mean flow, global energy,
    ! forces, but potentially also derived stuff such as Integral/Kolmogorov scales)
    ! NOTE: as for the RHS, some terms here depend on the grid as whole, and not just
    ! on individual blocks. This requires one to use the same staging concept as for the RHS.
    !-----------------------------------------------------------------------------
    subroutine STATISTICS_meta( physics, time, dt, u, g, x0, dx, rhs, stage, mask, lgt_n )
        implicit none

        character(len=*), intent(in) :: physics
        ! it may happen that some source terms have an explicit time-dependency
        ! therefore the general call has to pass time
        real(kind=rk), intent (in) :: time, dt

        ! block data, containg the state vector. In general a 4D field (3 dims+components)
        ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
        real(kind=rk), intent(inout) :: u(1:,1:,1:,1:)

        ! as you are allowed to compute the RHS only in the interior of the field
        ! you also need to know where 'interior' starts: so we pass the number of ghost points
        integer, intent(in) :: g

        ! for each block, you'll need to know where it lies in physical space. The first
        ! non-ghost point has the coordinate x0, from then on its just cartesian with dx spacing
        real(kind=rk), intent(in) :: x0(1:3), dx(1:3)

        ! output. Note assumed-shape arrays
        real(kind=rk), intent(inout) :: rhs(1:,1:,1:,1:)

        ! mask data. we can use different trees (4est module) to generate time-dependent/indenpedent
        ! mask functions separately. This makes the mask routines tree-level routines (and no longer
        ! block level) so the physics modules have to provide an interface to create the mask at a tree
        ! level. All parts of the mask shall be included: chi, boundary values, sponges.
        ! On input, the mask array is correctly filled. You cannot create the full mask here.
        real(kind=rk), intent(inout) :: mask(1:,1:,1:,1:)


        ! stage. there is 3 stages, init_stage, integral_stage and post_stage. The 1st and 3rd
        ! stages are called only once, not for every block.
        character(len=*), intent(in) :: stage

        !> number of active blocks (light data)
        integer(kind=ik), intent(inout)     :: lgt_n

        select case(physics)
        case ("ACM-new")
            call STATISTICS_ACM( time, dt, u, g, x0, dx, stage, rhs, mask )

        case ("ConvDiff-new")
            ! not implemented yet.
            !  call STATISTICS_convdiff( time, u, g, x0, dx, rhs, stage )

        case ("navier_stokes")
            call STATISTICS_NStokes( time, u, g, x0, dx, stage )

        case ('reactive_navier_stokes')
            call interface_reactive_navier_stokes("statistics", time=time, phi=u, phi_work=rhs, &
                 x0=x0, dx=dx, stage=stage, lgt_n=lgt_n)

        case default
            call abort(2152000, "[RHS_wrapper.f90]: physics_type is unknown"//physics)

        end select

    end subroutine STATISTICS_meta


    !-----------------------------------------------------------------------------
    ! setting the time step is very physics-dependent. Sometimes you have a CFL like
    ! condition, sometimes not. So each physic module must be able to decide on its
    ! time step. This routine is called for all blocks, the smallest returned dt is used.
    !-----------------------------------------------------------------------------
    subroutine GET_DT_BLOCK_meta( physics, time, u, Bs, g, x0, dx, dt )
        implicit none
        character(len=*), intent(in) :: physics
        ! it may happen that some source terms have an explicit time-dependency
        ! therefore the general call has to pass time
        real(kind=rk), intent (in) :: time

        ! block data, containg the state vector. In general a 4D field (3 dims+components)
        ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
        real(kind=rk), intent(inout) :: u(1:,1:,1:,1:)

        ! as you are allowed to compute the RHS only in the interior of the field
        ! you also need to know where 'interior' starts: so we pass the number of ghost points
        integer, intent(in) :: g
        integer, dimension(3), intent(in) :: Bs

        ! for each block, you'll need to know where it lies in physical space. The first
        ! non-ghost point has the coordinate x0, from then on its just cartesian with dx spacing
        real(kind=rk), intent(in) :: x0(1:3), dx(1:3)

        ! the dt for this block is returned to the caller:
        real(kind=rk), intent(out) :: dt

        select case (physics)
        case ('ACM-new')
            ! artificial compressibility
            call GET_DT_BLOCK_ACM( time, u, Bs, g, x0, dx, dt )

        case ('ConvDiff-new')
            ! convection-diffusion
            call GET_DT_BLOCK_convdiff( time, u, Bs, g, x0, dx, dt )

        case ('navier_stokes')
            ! navier stokes
            call GET_DT_BLOCK_NStokes( time, u, Bs, g, x0, dx, dt )

        case ('reactive_navier_stokes')
            call interface_reactive_navier_stokes("get_dt", time=time, phi=u, x0=x0, dx=dx, dt=dt)

        case default
            call abort('phycics module unkown.')

        end select

    end subroutine


    !-----------------------------------------------------------------------------
    ! main level wrapper for setting the initial condition on a block
    !-----------------------------------------------------------------------------
    subroutine INICOND_meta( physics, time, u, g, x0, dx)
        implicit none

        character(len=*), intent(in) :: physics
        ! it may happen that some source terms have an explicit time-dependency
        ! therefore the general call has to pass time
        real(kind=rk), intent (in) :: time

        ! block data, containg the state vector. In general a 4D field (3 dims+components)
        ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
        real(kind=rk), intent(inout) :: u(1:,1:,1:,1:)

        ! as you are allowed to compute the RHS only in the interior of the field
        ! you also need to know where 'interior' starts: so we pass the number of ghost points
        integer, intent(in) :: g

        ! for each block, you'll need to know where it lies in physical space. The first
        ! non-ghost point has the coordinate x0, from then on its just cartesian with dx spacing
        real(kind=rk), intent(in) :: x0(1:3), dx(1:3)


        select case (physics)
        case ("ACM-new")
            call INICOND_ACM( time, u, g, x0, dx)

        case ("ConvDiff-new")
            call INICOND_ConvDiff( time, u, g, x0, dx )

        case ("navier_stokes")
            call INICOND_NStokes( time, u, g, x0, dx )

        case ("reactive_navier_stokes")
            call interface_reactive_navier_stokes("initialize_data", phi=u, x0=x0, dx=dx)

        case default
            call abort(999,"[INICOND (metamodule):] unkown physics. Its getting hard to find qualified personel.")

        end select

    end subroutine INICOND_meta



    !-----------------------------------------------------------------------------
    ! wrapper for filter u -> u_tilde
    ! Note this function is completely
    ! independent of the grid and any MPI formalism, neighboring relations and the like.
    ! You just get a block data (e.g. ux, uy, uz, p) and apply your filter to it.
    ! Ghost nodes are assumed to be sync'ed.
    !-----------------------------------------------------------------------------
    subroutine FILTER_meta(physics, time, u, g, x0, dx, work_array, boundary_flag)
        implicit none
        !> physics type
        character(len=*), intent(in) :: physics
        !> time in physical units
        real(kind=rk), intent (in) :: time

        ! block data, containg the state vector. In general a 4D field (3 dims+components)
        ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
        real(kind=rk), intent(inout) :: u(1:,1:,1:,1:)

        !> number of ghost nodes
        integer, intent(in) :: g

        ! for each block, you'll need to know where it lies in physical space. The first
        ! non-ghost point has the coordinate x0, from then on its just cartesian with dx spacing
        real(kind=rk), intent(in) :: x0(1:3), dx(1:3)

        ! the work array is an additional array which can be used to store temporal
        ! values of the statevector field
        real(kind=rk), intent(inout) :: work_array(1:,1:,1:,1:)
        ! when implementing boundary conditions, it is necessary to now if the local field (block)
        ! is adjacent to a boundary, because the stencil has to be modified on the domain boundary.
        ! The boundary_flag tells you if the local field is adjacent to a domain boundary:
        ! boundary_flag(i) can be either 0, 1, -1,
        !  0: no boundary in the direction +/-e_i
        !  1: boundary in the direction +e_i
        ! -1: boundary in the direction - e_i
        ! currently only acessible in the local stage
        integer(kind=2), optional, intent(in):: boundary_flag(3)

        select case(physics)
        case ("ACM-new")
            call filter_ACM( time, u, g, x0, dx,  work_array)

        case ("ConvDiff-new")
            call abort(1009181817, "filter not implemented for convection-diffusion.")

        case ("navier_stokes")
            call filter_NStokes( time, u, g, x0, dx, work_array, boundary_flag)

        case ("reactive_navier_stokes")
            call interface_reactive_navier_stokes("filter_data", phi=u, phi_work=work_array, x0=x0, dx=dx)

        case default
            call abort(2152001, "ERROR [filter_wrapper.f90]: physics_type is unknown "//trim(adjustl(physics)))

        end select

    end subroutine FILTER_meta



end module module_physics_metamodule
