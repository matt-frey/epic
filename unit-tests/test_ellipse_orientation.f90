! =============================================================================
!                       Test ellipse orientation
!
!       This unit test checks the ellipse orientation computed using the
!       B matrix values, i.e., B11 and B12. The angle computed via the
!       B matrix is within [-pi/2, pi/2].
! =============================================================================
program test_ellipse_orientation
    use unit_test
    use constants, only : pi
    use parcel_container, only : parcels, n_parcels, parcel_alloc, parcel_dealloc
    use ellipse
    implicit none

    double precision :: origin(2) = (/0.0, 0.0/)
    double precision :: extent(2) =  (/0.2, 0.2/)
    integer :: iter
    integer :: grid(2) = (/2, 2/)
    double precision :: angle, B11, B12, V, a2, b2
    double precision, parameter :: lam = 3.0
    logical :: failed = .false.

    n_parcels = 1
    call parcel_alloc(1)

    parcels%position = 0.0
    parcels%volume = 0.25 * product(extent / (grid - 1))

    V = parcels%volume(1, 1)

    a2 = (lam * V) / pi
    b2 = V / (pi * lam)

    do iter = 0, 360

        angle = dble(iter) * pi / 180.0d0

        B11 = a2 * cos(angle) ** 2 + b2 * sin(angle) ** 2

        B12 = (a2 - b2) * sin(angle) * cos(angle)

        parcels%B(1, 1) = B11

        parcels%B(1, 2) = B12

        ! get_angle computes the angle in the first and fourth quadrant, i.e.,
        ! -pi/2 <= get_angle <= pi/2
        if (angle > pi / 2 .and. angle <= 3.0 * pi / 2) then
            angle = angle - pi
        else if (angle > 3.0 * pi / 2) then
            angle = angle - 2 * pi
        endif

        failed = (failed .or. abs(angle - get_angle(B11, B12, V)) > 1.0e-13)

    enddo

    call parcel_dealloc()

    call print_result_logical('Test ellipse orientation', failed)

end program test_ellipse_orientation
