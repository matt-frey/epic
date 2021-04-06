module taylorgreen
    use options, only : flow
    implicit none

    contains
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
