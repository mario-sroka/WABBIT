!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name compute_DFT.f90
!> \version 0.5
!> \author msr
!
!> \brief compute inverse fourier transformation with a given fourier coefficients
!> note: works only for a fixed (equidistant) grid!
!>
!! input:    - fourier coefficients, grid parameter \n
!! output:   - velocity component \n
!!
!!
!! = log ======================================================================================
!! \n
!! 13/08/19 - create
!
! ********************************************************************************************

subroutine compute_IDFT( params_physics, phi_hat, phi, x0, dx )

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> parameter struct
    type(type_params_rns), intent(inout)    :: params_physics

    !> velocity component
    real(kind=rk), intent(inout)            :: phi(:, :, :)

    !> array of fourier coefficients
    complex(kind=rk), intent(inout)         :: phi_hat(:, :, :)

    !> spacing and origin of block
    real(kind=rk), intent(in)               :: x0(:), dx(:)

    ! grid parameter
    integer(kind=ik)                        :: Bs(3), Ds(3), Lvl
    ! loop parameter
    integer(kind=ik)                        :: i, j, k, l
    ! wavenumber limit, wavenumber
    real(kind=rk)                           :: k0, w
    ! DFT loop limit, dummy array size
    integer(kind=ik)                        :: kmax, dummyN(3)
    ! dummy variable
    complex(kind=rk), allocatable, save     :: dummy(:,:,:), dummy2(:,:,:)
    ! start indexes
    integer(kind=ik)                        :: start_i(3)

!---------------------------------------------------------------------------------------------
! variables initialization

    ! note: DFT is computed on Bs-1 grid, to remove redundant nodes
    ! assume Bs here is the blocksize with redundant node included
    Bs = params_physics%Bs

    Lvl = params_physics%maxLvl
    ! note: domain size is calculated with max tree level, 
    ! so this value should correspond to the real treelevel
    Ds  = (Bs-1) * 2**Lvl

    ! DFT parameters
    k0   = params_physics%k0
    kmax = params_physics%kmax

    ! dummy array allocation
    ! at least a size of Bs is needed
    if (kmax < Bs(1)) then
        dummyN(1) = Bs(1)
    else
        dummyN(1) = kmax
    end if
    if (kmax < Bs(2)) then
        dummyN(2) = Bs(2)
    else
        dummyN(2) = kmax
    end if
    if (kmax < Bs(3)) then
        dummyN(3) = Bs(3)
    else
        dummyN(3) = kmax
    end if

    if ( .not. allocated(dummy) ) allocate(dummy(dummyN(1), dummyN(2), dummyN(3)), dummy2(dummyN(1), dummyN(2), dummyN(3)))

    ! start indexes for complex roots array
    start_i = int( x0/dx+1.0_rk , kind=ik)

!---------------------------------------------------------------------------------------------
! main body

    ! IDFT in all three dimensions
    ! use kmax as loop limit, to reduce computations

    do l = 1, kmax
        do k = 1, kmax
            do i = 1, Bs(1)
                dummy2(i,l,k) = phi_hat(1,l,k)

                do j = 2, kmax

                    dummy2(i,l,k) = dummy2(i,l,k) &
                                  + phi_hat(j,l,k) &
                                  * dconjg(params_physics%rootsX( (start_i(1)+i-2)*(j-1) + 1))

                end do
            end do
        end do
    end do

    do l = 1, Bs(1)
        do k = 1, kmax
            do i = 1, Bs(3)
                dummy(l,k,i) = dummy2(l,k,1)

                do j = 2, kmax

                    dummy(l,k,i) = dummy(l,k,i) &
                                 + dummy2(l,k,j) &
                                 * dconjg(params_physics%rootsZ( (start_i(3)+i-2)*(j-1) + 1))

                end do
            end do
        end do
    end do

    do l = 1, Bs(1)
        do k = 1, Bs(3)
            do i = 1, Bs(2)
                dummy2(l,i,k) = dummy(l,1,k)

                do j = 2, kmax

                    dummy2(l,i,k) = dummy2(l,i,k) &
                                  + dummy(l,j,k) &
                                  * dconjg(params_physics%rootsY( (start_i(2)+i-2)*(j-1) + 1))

                end do
            end do
        end do
    end do

    ! return real valued velocity field
    phi = real(dummy2(1:Bs(1),1:Bs(2),1:Bs(3)), kind=rk)

end subroutine compute_IDFT
