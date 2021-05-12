! =============================================================================
!                       Test ellipse multi merge
!
!         This unit test checks the symmetry by performing a mirrored
!         multi-merge (mirror axis x).
! =============================================================================
program test_ellipse_multi_merge_symmetry
    use unit_test
    use constants, only : pi, one, two, four
    use parcel_container
    use parcel_interpl, only : par2grid_elliptic_symmetry_check
    use parcel_merge, only : merge_ellipses
    use options, only : parcel_info, grid
    use parameters, only : update_parameters, extent
    use ellipse
    implicit none

    double precision :: error
    double precision :: volg(-1:3, 0:1, 1)

    grid = (/3, 3/)

    call update_parameters()

    call parcel_alloc(6)

    !
    ! muti-geometric merging
    !

    call parcel_setup

    ! geometric merge
    parcel_info%lambda = 5.0
    parcel_info%merge_type = 'multi-geometric'
    parcel_info%vfraction = 3

    call merge_ellipses(parcels)

    ! check result
    error = eval_max_error()

    call print_result_dp('Test ellipse multi-merge symmetry (geometric)', error)

    !
    ! muti-optimal merging
    !

    call parcel_setup

    ! optimal merge
    parcel_info%lambda = 5.0
    parcel_info%merge_type = 'multi-optimal'
    parcel_info%vfraction = 3

    call merge_ellipses(parcels)

    ! check result
    error = eval_max_error()

    call print_result_dp('Test ellipse multi-merge symmetry (optimal)', error)

    call parcel_dealloc

    contains

        subroutine parcel_setup
            double precision :: d, a1b1, a2b2
            integer :: n

            a1b1 = 1.44d0
            a2b2 = 0.25d0

            n_parcels = 6
            parcels%position(1, 1) = -0.5d0
            parcels%position(1, 2) = 0.2d0
            parcels%volume(1, 1) = a1b1 * pi
            parcels%B(1, 1) = 1.2d0 * a1b1
            parcels%B(1, 2) = -0.4d0

            parcels%position(2, 1) = -0.6d0
            parcels%position(2, 2) = 0.3d0
            parcels%volume(2, 1) = a2b2 * pi
            parcels%B(2, 1) = 0.8d0 * a2b2
            parcels%B(2, 2) = 0.5d0

            parcels%position(3, 1) = -0.3d0
            parcels%position(3, 2) = -0.1d0
            parcels%volume(3, 1) = a2b2 * pi
            parcels%B(3, 1) = 0.9d0 * a2b2
            parcels%B(3, 2) = -0.1d0

            !
            ! mirrored parcels
            !
            do n = 1, 3
                parcels%position(3+n, 1) = -parcels%position(n, 1)
                parcels%position(3+n, 2) =  parcels%position(n, 2)
                parcels%volume(3+n, 1) = parcels%volume(n, 1)
                parcels%B(3+n, 1) =  parcels%B(n, 1)
                parcels%B(3+n, 2) = -parcels%B(n, 2)
            enddo

        end subroutine parcel_setup

        function eval_max_error() result(max_err)
            double precision :: max_err

            max_err = 0.0d0
            max_err = max(max_err, abs(dble(n_parcels - 2)))

            call par2grid_elliptic_symmetry_check(parcels, parcels%volume, volg)

            max_err = max(max_err, maxval(abs(volg)))

        end function eval_max_error

end program test_ellipse_multi_merge_symmetry