! =============================================================================
!                       EPIC - Elliptical Parcel-in-Cell
! =============================================================================
program epic
    use constants, only : max_num_parcels, zero
    use timer
    use field_diagnostics
    use parser, only : read_config_file
    use parcel_container
    use parcel_bc
    use parcel_split, only : split_ellipses, split_timer
    use parcel_point, only : point_split, point_merge
    use parcel_merge, only : merge_ellipses, merge_timer
    use parcel_correction, only : init_parcel_correction, &
                                  apply_laplace,          &
                                  apply_gradient,         &
                                  lapl_corr_timer,        &
                                  grad_corr_timer
    use parcel_diagnostics
    use parcel_hdf5
    use fields
    use field_hdf5
    use tri_inversion, only : init_inversion, vor2vel_timer, vtend_timer
    use parcel_interpl
    use parcel_init, only : init_parcels, init_timer
    use ls_rk4
    use h5_utils, only : initialise_hdf5, finalise_hdf5
    implicit none

    integer :: epic_timer

    ! Read command line (verbose, filename, etc.)
    call parse_command_line

    ! Create the model
    call pre_run

    ! Run the model
    call run

    ! Deallocate memory
    call post_run

    contains

        subroutine pre_run
            use options, only : field_file, field_tol, output

            call register_timer('epic', epic_timer)
            call register_timer('vol2grid', vol2grid_timer)
            call register_timer('par2grid', par2grid_timer)
            call register_timer('grid2par', grid2par_timer)
            call register_timer('parcel split', split_timer)
            call register_timer('parcel merge', merge_timer)
            call register_timer('laplace correction', lapl_corr_timer)
            call register_timer('gradient correction', grad_corr_timer)
            call register_timer('parcel init', init_timer)
            call register_timer('parcel hdf5', hdf5_parcel_timer)
            call register_timer('field hdf5', hdf5_field_timer)
            call register_timer('vor2vel', vor2vel_timer)
            call register_timer('vorticity tendency', vtend_timer)
            call register_timer('parcel push', rk4_timer)
#ifndef NDEBUG
            call register_timer('symmetric vol2grid', sym_vol2grid_timer)
#endif

            call start_timer(epic_timer)

            call initialise_hdf5

            ! parse the config file
            call read_config_file

            call parcel_alloc(max_num_parcels)

            call init_parcels(field_file, field_tol)

            call ls_rk4_alloc(max_num_parcels)

            call init_inversion

            call init_parcel_correction

            call init_parcel_diagnostics

            call field_default

            call par2grid

            if (output%h5_write_fields) then
                call create_h5_field_file(trim(output%h5_basename), output%h5_overwrite)
            endif

            if (output%h5_write_parcels) then
                call create_h5_parcel_file(trim(output%h5_basename), output%h5_overwrite)
            endif

        end subroutine


        subroutine run
            use options, only : time, output, verbose, parcel
            double precision :: t    = zero ! current time
            double precision :: dt   = zero ! time step
            integer          :: iter = 1    ! simulation iteration
            integer          :: nfw  = 0    ! number of field writes to h5
            integer          :: npw  = 0    ! number of parcel writes to h5
            integer          :: cor_iter    ! iterator for parcel correction

            do while (t <= time%limit)

                dt = get_time_step()

#ifdef ENABLE_VERBOSE
                if (verbose) then
                    print "(a15, f0.4)", "time:          ", t
                    print "(a15, i0)", "iteration:     ", iter
                endif
#endif

                ! make sure we always write initial setup
                if (output%h5_write_fields .and. &
                    (mod(iter - 1, output%h5_field_freq) == 0)) then
#ifndef NDEBUG
                    call vol2grid_symmetry_error
#endif
                    call write_h5_field_step(nfw, t, dt)
                endif

                if (output%h5_write_parcels .and. &
                    (mod(iter - 1, output%h5_parcel_freq) == 0)) then
                    call write_h5_parcel_step(npw, t, dt)
                endif

                call ls_rk4_step(dt)

                if (mod(iter, parcel%merge_freq) == 0) then
                    if (parcel%is_elliptic) then
                        call merge_ellipses(parcels)
                    else
                        call point_merge(parcel%vfraction)
                    endif
                endif

                if (mod(iter, parcel%split_freq) == 0) then
                    if (parcel%is_elliptic) then
                        call split_ellipses(parcels, parcel%lambda, parcel%vmaxfraction)
                    else
                        call point_split(parcel%lambda, parcel%prefactor)
                    endif
                endif

                if (mod(iter, parcel%correction_freq) == 0) then
                    call vol2grid
                    do cor_iter=1,parcel%correction_iters
                        if (parcel%apply_laplace) then
                            call apply_laplace(volg)
                            call vol2grid
                        endif
                        if (parcel%apply_gradient) then
                            call apply_gradient(volg,parcel%gradient_pref,parcel%max_compression)
                            call vol2grid
                        end if
                    end do
                 endif


                t = t + dt
                iter = iter + 1
            end do

            call par2grid

            ! write final step
            if (output%h5_write_fields) then
#ifndef NDEBUG
                call vol2grid_symmetry_error
#endif
                call write_h5_field_step(nfw, t, dt)
            endif

            if (output%h5_write_parcels) then
                call write_h5_parcel_step(npw, t, dt)
            endif

        end subroutine run

        subroutine post_run
            call parcel_dealloc
            call ls_rk4_dealloc
            call finalise_hdf5

            call stop_timer(epic_timer)

            call write_time_to_csv(output%h5_basename)
            call print_timer
        end subroutine


        function get_time_step() result(dt)
            use options, only : time
            double precision :: dt
            double precision :: max_vorticity
            double precision :: H, S11, S12, S21, S22, gmax
            integer          :: i, j

            H = epsilon(zero)
            if (parcel%is_elliptic .and. time%is_adaptive) then
                do i = 0, nx-1
                    do j = 0, nz
                        S11 = velgradg(j, i, 1)
                        S12 = velgradg(j, i, 2)
                        S21 = velgradg(j, i, 3)
                        S22 = velgradg(j, i, 4)
                        H = max(H, (S11 - S22) ** 2 + (S12 + S21) ** 2)
                    enddo
                enddo

                gmax = f12 * dsqrt(H)
                dt = min(time%dt_max, time%alpha / gmax)

            else if (time%is_adaptive) then
                    ! adaptive time stepping according to
                    ! https://doi.org/10.1002/qj.3319
                    max_vorticity = maxval(abs(vortg))
                    dt = min(time%alpha / max_vorticity, time%dt)
            else
                dt = time%dt
            endif

            if (dt <= zero) then
                print "(a10, f0.2, a6)", "Time step ", dt, " <= 0!"
                stop
            endif
        end function get_time_step


    ! Get the file name provided via the command line
    subroutine parse_command_line
        use options, only : filename, verbose
        integer                          :: i
        character(len=512)               :: arg

        filename = ''
        i = 0
        do
            call get_command_argument(i, arg)
            if (len_trim(arg) == 0) then
                exit
            endif

            if (arg == '--config') then
                i = i + 1
                call get_command_argument(i, arg)
                filename = trim(arg)
            else if (arg == '--help') then
                print *, 'Run code with "./epic --config [config file]"'
                stop
#ifdef ENABLE_VERBOSE
            else if (arg == '--verbose') then
                verbose = .true.
#endif
            endif
            i = i+1
        end do

        if (filename == '') then
            print *, 'No configuration file provided. Run code with "./epic --config [config file]"'
            stop
        endif

#ifdef ENABLE_VERBOSE
        ! This is the main application of EPIC
        if (verbose) then
            print *, 'Running EPIC with "', trim(filename), '"'
        endif
#endif
    end subroutine parse_command_line
end program epic
