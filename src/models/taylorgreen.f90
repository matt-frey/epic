! =============================================================================
!                           2D TAYLOR-GREEN FLOW
!
!                   u(x, y) = A * cos(ax + d) * sin(by + e)
!                   v(x, y) = B * sin(ax + d) * cos(by + e)
!
!                   The code uses following variable names:
!                           flow%freq  = (/a, b/)
!                           flow%amp   = (/A, B/)
!                           flow%phase = (/d, e/)
! =============================================================================
module taylorgreen
    use options, only : flow
    use constants, only : zero
    use parcel_container, only : parcels, n_parcels
    use fields
    implicit none

    contains
        subroutine taylorgreen_init
            double precision :: pos(2)
            integer          :: i, j, n

            do n = 1, n_parcels
                parcels%velocity(n, :) = get_flow_velocity(parcels%position(n, :))
            enddo

            do i = 0, nx-1
                do j = -1, nz+1
                    call get_position(i, j, pos)

                    velog(j, i, :) = get_flow_velocity(pos)

                    velgradg(j, i, :) = get_flow_gradient(pos)

                    if ((j > -1) .and. (j < nz+1)) then
                        vortg(j, i) = get_flow_vorticity(pos)
                    endif
                enddo
            enddo
        end subroutine taylorgreen_init


        function get_flow_velocity(pos) result(vel)
            double precision, intent(in) :: pos(2)
            double precision             :: xx, zz
            double precision             :: vel(2)

            xx = flow%freq(1) * pos(1) + flow%phase(1)
            zz = flow%freq(2) * pos(2) + flow%phase(2)

            vel(1) = flow%amp(1) * cos(xx) * sin(zz)
            vel(2) = flow%amp(2) * sin(xx) * cos(zz)
        end function get_flow_velocity

        ! grad ordering : dudx, dudy, dvdx, dvdy
        function get_flow_gradient(pos) result(grad)
            double precision, intent(in) :: pos(2)
            double precision             :: xx, zz
            double precision             :: grad(4)

            xx = flow%freq(1) * pos(1) + flow%phase(1)
            zz = flow%freq(2) * pos(2) + flow%phase(2)

            ! du/dx = - a * A * sin(xx) * sin(zz)
            grad(1) = - flow%freq(1) * flow%amp(1) * sin(xx) * sin(zz)

            ! du/dy = b * A * cos(xx) * cos(zz)
            grad(2) = flow%freq(2) * flow%amp(1) * cos(xx) * cos(zz)

            ! dv/dx = a * B * cos(xx) * np.cos(zz)
            grad(3) = flow%freq(1) * flow%amp(2) * cos(xx) * cos(zz)

            ! dv/dy = - b * B * sin(xx) * sin(zz)
            grad(4) = - flow%freq(2) * flow%amp(2) * sin(xx) * sin(zz)

        end function get_flow_gradient

        function get_flow_vorticity(pos) result(omega)
            double precision, intent(in) :: pos(2)
            double precision             :: xx, zz
            double precision             :: omega

            xx = flow%freq(1) * pos(1) + flow%phase(1)
            zz = flow%freq(2) * pos(2) + flow%phase(2)

            omega = (flow%amp(2) * flow%freq(1)     &
                   - flow%amp(1) * flow%freq(2))    &
                   * cos(xx) * cos(zz)
        end function get_flow_vorticity

end module taylorgreen
