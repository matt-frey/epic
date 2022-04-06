! =============================================================================
!                             Field diagnostics
! =============================================================================
module field_diagnostics
    use parameters, only : vcell, vcelli, nx, nz, ngridi, ncelli
    use fields
    use timer, only : start_timer, stop_timer
    implicit none

    integer :: field_stats_timer

    double precision :: rms_v,      &       ! rms volume error
                        abserr_v,   &       ! max absolute normalised volume error
                        max_npar,   &       ! max num parcels per cell
                        min_npar,   &       ! min num parcels per cell
                        avg_npar,   &       ! average num parcels per cell
                        avg_nspar,  &       ! average num small parcels per cell
                        max_dudx,   &       ! maximum strain du/dx value in magnitude
                        max_dudy,   &       ! maximum strain du/dy value in magnitude
                        max_dudz,   &       ! maximum strain du/dz value in magnitude
                        max_dvdx,   &       ! maximum strain dv/dx value in magnitude
                        max_dvdy,   &       ! maximum strain dv/dy value in magnitude
                        max_dvdz,   &       ! maximum strain dv/dz value in magnitude
                        max_dwdx,   &       ! maximum strain dw/dx value in magnitude
                        max_dwdy,   &       ! maximum strain dw/dy value in magnitude
                        max_dwdz,   &       ! maximum strain dw/dz value in magnitude
                        max_xvortg, &       ! maximum x-vorticity value in magnitude
                        max_yvortg, &       ! maximum y-vorticity value in magnitude
                        max_zvortg          ! maximum z-vorticity value in magnitude
    contains

        subroutine calculate_field_diagnostics
            double precision :: sqerrsum

            call start_timer(field_stats_timer)

            ! do not take halo cells into account
            sqerrsum = sum((volg(0:nz, :, :) - vcell) ** 2)
            rms_v = dsqrt(sqerrsum * ngridi) * vcelli

            abserr_v = maxval(abs(volg(0:nz, :, :)  - vcell)) * vcelli

            max_npar = maxval(nparg(0:nz-1, :, :))

            min_npar = minval(nparg(0:nz-1, :, :))

            avg_npar = sum(nparg(0:nz-1, :, :)) * ncelli

            avg_nspar = sum(nsparg(0:nz-1, :, :)) * ncelli

            max_dudx = maxval(abs(velgradg(0:nz, :, :, 1)))
            max_dudy = maxval(abs(velgradg(0:nz, :, :, 2)))
            max_dudz = maxval(abs(velgradg(0:nz, :, :, 4) + vortg(0:nz, :, :, 2)))
            max_dvdx = maxval(abs(velgradg(0:nz, :, :, 2) + vortg(0:nz, :, :, 3)))
            max_dvdy = maxval(abs(velgradg(0:nz, :, :, 3)))
            max_dvdz = maxval(abs(velgradg(0:nz, :, :, 5) - vortg(0:nz, :, :, 1)))
            max_dwdx = maxval(abs(velgradg(0:nz, :, :, 4)))
            max_dwdy = maxval(abs(velgradg(0:nz, :, :, 5)))
            max_dwdz = maxval(abs(velgradg(0:nz, :, :, 1) + velgradg(0:nz, :, :, 3)))

            max_xvortg = maxval(abs(vortg(0:nz, :, :, 1)))
            max_yvortg = maxval(abs(vortg(0:nz, :, :, 2)))
            max_zvortg = maxval(abs(vortg(0:nz, :, :, 3)))


            call stop_timer(field_stats_timer)

        end subroutine calculate_field_diagnostics

end module field_diagnostics
