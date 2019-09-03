!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name find_neighbor_corner_3D.f90
!> \version 0.5
!> \author msr
!
!> \brief find neighbor on block corner
!
!> \details  input:
!!                   - heavy and light data id
!!                   - light data array and max treelevel
!!                   - direction for neighbor search
!!                   - list of active blocks
!!
!!            output:
!!                   - neighbor list array
!!
! --------------------------------------------------------------------------------------------
!> \details  neighbor codes: \n
! ---------------
!> for imagination:
!!                   - 6-sided dice with '1'-side on top, '6'-side on bottom, '2'-side in front
!!                   - edge: boundary between two sides - use sides numbers for coding
!!                   - corner: between three sides - so use all three sides numbers
!!                   - block on higher/lower level: block shares face/edge and one unique corner,
!!                     so use this corner code in second part of neighbor code
!!
!! \image html neighborcode.svg "Neighborcode 3D" width=250
!!
!! faces:  '__1/___', '__2/___', '__3/___', '__4/___', '__5/___', '__6/___' \n
!! edges:  '_12/___', '_13/___', '_14/___', '_15/___' \n
!!         '_62/___', '_63/___', '_64/___', '_65/___' \n
!!         '_23/___', '_25/___', '_43/___', '_45/___' \n
!! corner: '123/___', '134/___', '145/___', '152/___' \n
!!         '623/___', '634/___', '645/___', '652/___' \n
!!\n
!! complete neighbor code array, 74 possible neighbor relations \n
!! neighbors = (/'__1/___', '__2/___', '__3/___', '__4/___', '__5/___', '__6/___', '_12/___', '_13/___', '_14/___', '_15/___',
!!                '_62/___', '_63/___', '_64/___', '_65/___', '_23/___', '_25/___', '_43/___', '_45/___', '123/___', '134/___',
!!                '145/___', '152/___', '623/___', '634/___', '645/___', '652/___', '__1/123', '__1/134', '__1/145', '__1/152',
!!               '__2/123', '__2/623', '__2/152', '__2/652', '__3/123', '__3/623', '__3/134', '__3/634', '__4/134', '__4/634',
!!               '__4/145', '__4/645', '__5/145', '__5/645', '__5/152', '__5/652', '__6/623', '__6/634', '__6/645', '__6/652',
!!               '_12/123', '_12/152', '_13/123', '_13/134', '_14/134', '_14/145', '_15/145', '_15/152', '_62/623', '_62/652',
!!               '_63/623', '_63/634', '_64/634', '_64/645', '_65/645', '_65/652', '_23/123', '_23/623', '_25/152', '_25/652',
!!               '_43/134', '_43/634', '_45/145', '_45/645' \)
! --------------------------------------------------------------------------------------------
!> \details
!! = log ======================================================================================
!! \n
!! 30/01/17 - create
!
! ********************************************************************************************

subroutine find_neighbor_corner_3D(params, heavy_id, light_id, lgt_block, max_treelevel, dir, hvy_neighbor, &
    lgt_n, lgt_sortednumlist, error)

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none
    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> heavy data id
    integer(kind=ik), intent(in)        :: heavy_id
    !> light data id
    integer(kind=ik), intent(in)        :: light_id
    !> max treelevel
    integer(kind=ik), intent(in)        :: max_treelevel
    !> light data array
    integer(kind=ik), intent(in)        :: lgt_block(:, :)
    !> direction for neighbor search
    character(len=7), intent(in)        :: dir
    !> number of active blocks (light data)
    integer(kind=ik), intent(in)        :: lgt_n
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), intent(in)     :: lgt_sortednumlist(:,:)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(inout)     :: hvy_neighbor(:,:)
    logical, intent(inout)              :: error

    ! mesh level
    integer(kind=ik)                    :: level
    ! treecode varaibles
    integer(kind=ik)                    :: my_treecode(max_treelevel), neighbor(max_treelevel), virt_treecode(max_treelevel)

    ! return value from function "does_block_exist"
    logical                             :: exists

    ! variable to show if there is a valid corner neighbor
    logical                             :: lvl_down_neighbor

    ! auxiliary variables
    integer(kind=ik)                    :: list_id, virt_code

    ! neighbor light data id
    integer(kind=ik)                    :: neighbor_light_id, tree_id

    ! variable for non-peridoic neighbor block detection
    logical                             :: non_periodic

    ! data_bounds array
    integer(kind=ik)                    :: data_bounds(2,3)

!---------------------------------------------------------------------------------------------
! interfaces

!---------------------------------------------------------------------------------------------
! variables initialization

    my_treecode= lgt_block( light_id, 1:max_treelevel )
    level      = lgt_block( light_id, max_treelevel + IDX_MESH_LVL )
    tree_id    = lgt_block( light_id, max_treelevel + IDX_TREE_ID )

    ! it is not always possible to have a corner neighbor on a coarser level, because
    ! the 4/8 sister blocks are complete. That means, a block cannot have all neighbors
    ! coarser. but has instead to have some on the same level.
    lvl_down_neighbor = .false.

    virt_code = -1
    list_id   = -1


!---------------------------------------------------------------------------------------------
! main body

    ! set auxiliary variables
    select case(dir)

        case('123/___')
            list_id    = 19

            ! virtual treecodes for neighbor search on higher level
            virt_code = 6

            ! is it possible to have a coarser neighbor?
            if ( my_treecode( level ) == 6 ) then
                lvl_down_neighbor = .true.
            end if

        case('134/___')
            list_id    = 20

            ! virtual treecodes for neighbor search on higher level
            virt_code = 7

            ! is it possible to have a coarser neighbor?
            if ( my_treecode( level ) == 7 ) then
                lvl_down_neighbor = .true.
            end if

        case('145/___')
            list_id    = 21

            ! virtual treecodes for neighbor search on higher level
            virt_code = 5

            ! is it possible to have a coarser neighbor?
            if ( my_treecode( level ) == 5 ) then
                lvl_down_neighbor = .true.
            end if

        case('152/___')
            list_id    = 22

            ! virtual treecodes for neighbor search on higher level
            virt_code = 4

            ! is it possible to have a coarser neighbor?
            if ( my_treecode( level ) == 4 ) then
                lvl_down_neighbor = .true.
            end if

        case('623/___')
            list_id    = 23

            ! virtual treecodes for neighbor search on higher level
            virt_code = 2

            ! is it possible to have a coarser neighbor?
            if ( my_treecode( level ) == 2 ) then
                lvl_down_neighbor = .true.
            end if

        case('634/___')
            list_id    = 24

            ! virtual treecodes for neighbor search on higher level
            virt_code = 3

            ! is it possible to have a coarser neighbor?
            if ( my_treecode( level ) == 3 ) then
                lvl_down_neighbor = .true.
            end if

        case('645/___')
            list_id    = 25

            ! virtual treecodes for neighbor search on higher level
            virt_code = 1

            ! is it possible to have a coarser neighbor?
            if ( my_treecode( level ) == 1 ) then
                lvl_down_neighbor = .true.
            end if

        case('652/___')
            list_id    = 26

            ! virtual treecodes for neighbor search on higher level
            virt_code = 0

            ! is it possible to have a coarser neighbor?
            if ( my_treecode( level ) == 0 ) then
                lvl_down_neighbor = .true.
            end if

    end select

    ! calculate treecode for neighbor on same level
    call adjacent_block_3D( my_treecode, neighbor, dir, level, max_treelevel)
    ! check existence of neighbor block and find light data id
    call does_block_exist(neighbor, exists, neighbor_light_id, &
                                lgt_sortednumlist, lgt_n, tree_id)

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! non periodic boundary
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    ! compute receiver bounds, need for lvl 1 meshes in detect non periodic neighbor subrioutine
    ! note: use level difference of zero and neighborhood ids at same level, include_redundant id = 2
    call set_recv_bounds( params, data_bounds, list_id, 0, 2, 'receiver')

    ! detect, if a neighbor is over a non-periodic domain boundary
    call detect_non_periodic_neighbor_3D( my_treecode, neighbor, params%periodic_BC, level, &
                                          params%domain_size, params%Bs, data_bounds, non_periodic )

    ! reset neighbor_light_id for non-periodic cases
    if (non_periodic) neighbor_light_id = -1

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    if (exists) then
        ! neighbor on same level
        hvy_neighbor( heavy_id, list_id ) = neighbor_light_id

    else

        ! neighbor could be one level down
        neighbor( level ) = -1
        ! check existence of neighbor block
        call does_block_exist(neighbor, exists, neighbor_light_id, &
                                lgt_sortednumlist, lgt_n, tree_id)

        if ( exists .and. lvl_down_neighbor ) then
            ! neigbor is one level down
            hvy_neighbor( heavy_id, list_id ) = neighbor_light_id

        elseif ( .not.(exists) ) then
            ! neighbor could be on level up
            ! virtual treecode, one level up
            virt_treecode = my_treecode
            virt_treecode( level+1 ) = virt_code

            ! calculate treecode for neighbor on same level (virtual level)
            call adjacent_block_3D( virt_treecode, neighbor, dir, level+1, max_treelevel)
            ! check existence of neighbor block
            call does_block_exist(neighbor, exists, neighbor_light_id, &
                                lgt_sortednumlist, lgt_n, tree_id)


            if (exists) then
                ! neigbor is one level up
                hvy_neighbor( heavy_id, list_id ) = neighbor_light_id

            else
                ! error case
                write(*,*) "find_neighbor_corner_3D: my treecode", my_treecode, "dir", dir, "neighbor treecode", virt_treecode
                error = .true.
            end if

        end if

    end if

end subroutine find_neighbor_corner_3D
