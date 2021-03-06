subroutine post_add_two_masks(params)
    use module_precision
    use module_mesh
    use module_params
    use module_IO
    use module_mpi
    use module_operators
    use module_physics_metamodule
    use module_time_step
    use module_stl_file_reader
    use module_helpers
    use module_forest

    implicit none

    type (type_params), intent(inout)  :: params
    character(len=80) :: fname_ini, fname1, fname2, fname_out

    integer(kind=ik), allocatable      :: lgt_block(:, :)
    real(kind=rk), allocatable         :: hvy_block(:, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_work(:, :, :, :, :, :)
    real(kind=rk), allocatable         :: hvy_tmp(:, :, :, :, :)
    integer(kind=ik), allocatable      :: hvy_neighbor(:,:)
    integer(kind=ik), allocatable      :: lgt_active(:,:), hvy_active(:,:)
    integer(kind=tsize), allocatable   :: lgt_sortednumlist(:,:,:)
    integer(kind=ik), allocatable      :: lgt_n(:), hvy_n(:)
    integer :: hvy_id, lgt_id, fsize, j, tree_id

    integer(kind=ik) :: iteration, Bs(1:3), tc_length1, dim, tc_length2, N1, N2, tree_N
    real(kind=rk) :: time, domain(1:3)

    !-----------------------------------------------------------------------------------------------------
    ! get values from command line (filename and level for interpolation)
    call get_command_argument(2, fname_ini)

    ! does the user need help?
    if (fname_ini=='--help' .or. fname_ini=='--h' .or. fname_ini=='-h') then
        if (params%rank==0) then
            write(*,*) "------------------------------------------------------------------"
            write(*,*) "./wabbit-post --add-two-masks mask_001.h5 mask_002.h5 output.h5"
            write(*,*) "------------------------------------------------------------------"
            write(*,*) ""
            write(*,*) ""
            write(*,*) ""
            write(*,*) ""
            write(*,*) "------------------------------------------------------------------"
        end if
        return
    endif

    call get_command_argument(2, fname1)
    call check_file_exists(fname1)

    call get_command_argument(3, fname2)
    call check_file_exists(fname2)

    call get_command_argument(4, fname_out)

    call read_attributes(fname1, N1, time, iteration, domain, params%Bs, tc_length1, params%dim)
    call read_attributes(fname2, N2, time, iteration, domain, params%Bs, tc_length2, params%dim)


    params%number_blocks = 2*max(N1,N2) ! just to get some memory:
    params%domain_size = domain
    params%max_treelevel = max( tc_length2, tc_length1 )
    params%min_treelevel = 1
    params%n_eqn = 1
    params%n_ghosts = 4
    params%forest_size = 20
    fsize = params%forest_size
    params%order_predictor = "multiresolution_4th"
    params%block_distribution = "sfc_hilbert"
    params%time_step_method = 'none'


    ! we have to allocate grid if this routine is called for the first time
    call allocate_forest(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, &
    hvy_active, lgt_sortednumlist, hvy_work, hvy_tmp=hvy_tmp, hvy_n=hvy_n, lgt_n=lgt_n)

    ! The ghost nodes will call their own setup on the first call, but for cleaner output
    ! we can also just do it now.
    call init_ghost_nodes( params )

    hvy_neighbor = -1
    lgt_n = 0 ! reset number of active light blocks
    hvy_n = 0
    tree_n = 0 ! reset number of trees in forest

    call read_field2tree(params, (/fname1/), 1, 1, tree_n, &
    lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
    hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor)

    call read_field2tree(params, (/fname2/), 1, 2, tree_n, &
    lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
    hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor)



    ! call copy_tree(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
    ! hvy_block, hvy_active, hvy_n, hvy_neighbor, tree_id_dest=2, tree_id_source=1)
    !
    ! call create_active_and_sorted_lists( params, lgt_block, lgt_active, &
    ! lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_n)
    !
    ! call prune_tree( params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
    ! hvy_block, hvy_active, hvy_n, hvy_neighbor, tree_id=2)

    call create_active_and_sorted_lists( params, lgt_block, lgt_active, &
    lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_n)




    !     tree_id = 1
    ! N1 =  max_active_level(lgt_block,lgt_active(:,tree_id),lgt_n(tree_id))
    !     call refine_tree(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
    !         hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor, tree_id)
    !
    !         N2= max_active_level(lgt_block,lgt_active(:,tree_id),lgt_n(tree_id))
    ! write(*,*) N1, n2
    ! call create_active_and_sorted_lists( params, lgt_block, lgt_active, &
    ! lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_n)

    ! j = 2
    ! call balance_load( params, lgt_block, hvy_block,  hvy_neighbor, &
    ! lgt_active(:, j), lgt_n(j), lgt_sortednumlist(:,:,j), hvy_active(:, j), hvy_n(j), hvy_tmp )
    !
    !     call create_active_and_sorted_lists( params, lgt_block, lgt_active, &
    !     lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_n)
    ! !
    !     call add_pruned_to_full_tree( params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
    !     hvy_block, hvy_active, hvy_n, hvy_neighbor, tree_id_pruned=2, tree_id_full=1)
    !
    !         call create_active_and_sorted_lists( params, lgt_block, lgt_active, &
    !         lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_n)

    call add_two_trees(params, tree_n, lgt_block, lgt_active, lgt_n, lgt_sortednumlist, &
    hvy_block, hvy_active, hvy_n, hvy_tmp, hvy_neighbor, tree_id1=1, tree_id2=2)

    call create_active_and_sorted_lists( params, lgt_block, lgt_active, &
    lgt_n, hvy_active, hvy_n, lgt_sortednumlist, tree_n)

    tree_id = 1
    call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active(:,tree_id),&
    lgt_n(tree_id), lgt_sortednumlist(:,:,tree_id), hvy_active(:,tree_id), hvy_n(tree_id) )

    tree_id = 2
    call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active(:,tree_id),&
    lgt_n(tree_id), lgt_sortednumlist(:,:,tree_id), hvy_active(:,tree_id), hvy_n(tree_id) )

    call write_tree_field(fname_out, params, lgt_block, lgt_active, hvy_block, &
    lgt_n, hvy_n, hvy_active, dF=1, tree_id=1, time=time, iteration=iteration )


    call deallocate_forest(params, lgt_block, hvy_block, hvy_neighbor, lgt_active, hvy_active, &
    lgt_sortednumlist, hvy_work, hvy_tmp, hvy_n, lgt_n )

    ! make a summary of the program parts, which have been profiled using toc(...)
    ! and print it to stdout
    call summarize_profiling( WABBIT_COMM )
end subroutine
