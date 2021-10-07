program genspec
    use sta2dfft
    use constants, only : pi, twopi, f14, f12, zero, one, two
    use h5_reader
    use h5_writer, only : get_step_group_name
    implicit none

    ! Grid dimensions:
    integer :: nx, nz

    ! Width and height of the domain:
    double precision :: extent(2)

    double precision :: dx, dz

    ! Array to contain data:
    double precision, allocatable :: pp(:, :)

    ! Its Fourier transform:
    double precision, allocatable :: ss(:, :)

    ! The spectrum:
    double precision, allocatable :: spec(:)

    ! x and z wavenumbers:
    double precision, allocatable :: rkx(:), hrkx(:), rkz(:)

    ! Generic arrays needed for the FFTs:
    double precision, allocatable :: xtrig(:), ztrig(:)
    integer :: xfactors(5), zfactors(5)

    ! Wavenumber magnitude used for the spectrum:
    integer, allocatable :: kmag(:, :)

    ! Other work variables:
    double precision :: scx, rkxmax, scz, rkzmax, delk, delki, snorm
    double precision :: pnorm, pnorms
    integer :: kxc, kmax, kx, kz, k

    character(len=512) :: filename
    character(len=64)  :: dset
    integer            :: step

    call initialise_hdf5

    call parse_command_line

    call get_domain

    dx = extent(1) / dble(nx)
    dz = extent(2) / dble(nz)

    print *, 'Field: ' // trim(dset)
    print '(a23, i5, a6, i5)', 'Grid dimensions: nx = ', nx, ' nz = ', nz

    call alloc_arrays

    ! read data into array pp:
    call read_data

    pnorm = dx * dz * (f12 * (sum(pp(0, :) ** 2 + pp(nz, :) ** 2)) + sum(pp(1:nz-1, :) ** 2))
    write(*,*) ' Integral of b^2 dx dy = ', pnorm

    !---------------------------------------------------------------------
    ! Set up FFTs:
    call init2dfft(nx, nz, extent(1), extent(2), xfactors, zfactors, xtrig, ztrig, hrkx, rkz)

    !Define x wavenumbers:
    rkx(0) = zero
    do kx = 1, nx / 2 - 1
        kxc = nx - kx
        rkx(kx ) = hrkx(2 * kx)
        rkx(kxc) = hrkx(2 * kx)
    enddo
    rkx(nx / 2) = hrkx(nx)

    !Initialise arrays for computing the spectrum:
    scx = twopi / extent(1)
    rkxmax = scx * dble(nx / 2)
    scz = pi / extent(2)
    rkzmax = scz * dble(nz)
    delk = sqrt(scx ** 2 + scz **2)
    delki = one / delk
    kmax = nint(dsqrt(rkxmax ** 2 + rkzmax ** 2) * delki)
    do kz = 1, nz
        do kx = 0, nx - 1
            kmag(kx, kz) = nint(dsqrt(rkx(kx) ** 2 + rkz(kz) ** 2) * delki)
        enddo
    enddo
    do kx = 0, nx - 1
        kmag(kx, 0) = nint(rkx(kx) * delki)
    enddo

    !Compute spectrum multiplication factor (snorm) so that the sum
    !of the spectrum is equal to the L2 norm of the original field:
    snorm = two * dx * dz * delki

    !---------------------------------------------------------------------
    !Compute spectrum:

    !Transform data in pp to spectral space:
    call ptospc_fc(nx, nz, pp, ss, xfactors, zfactors, xtrig, ztrig)

    do k = 0, kmax
        spec(k) = zero
    enddo

    !x and y-independent mode:
    k = kmag(0, 0)
    spec(k) = spec(k) + f14 * ss(0, 0) ** 2

    !y-independent mode:
    do kx = 1, nx - 1
        k = kmag(kx, 0)
        spec(k) = spec(k) + f12 * ss(kx, 0) ** 2
    enddo

    !x-independent mode:
    do kz = 1, nz
        k = kmag(0, kz)
        spec(k) = spec(k) + f12 * ss(0, kz) ** 2
    enddo

    !All other modes:
    do kz = 1, nz
        do kx = 1, nx - 1
            k = kmag(kx, kz)
            spec(k) = spec(k) + ss(kx, kz) ** 2
        enddo
    enddo

     !Normalise:
    do k = 0, kmax
        spec(k) = snorm * spec(k)
    enddo

    !---------------------------------------------------------------------
    !Write spectrum contained in spec(k):
    call write_spectrum

    call dealloc_arrays

    call finalise_hdf5

    contains

        subroutine get_domain
            integer(hid_t)   :: h5handle
            double precision :: lower(2)
            ! read domain dimensions
            call open_h5_file(trim(filename), H5F_ACC_RDONLY_F, h5handle)
            call read_h5_box(h5handle, nx, nz, extent, lower)
            call close_h5_file(h5handle)
        end subroutine get_domain

        subroutine alloc_arrays
            allocate(pp(0:nz, 0:nx - 1))
            allocate(ss(0:nx - 1, 0:nz))
            allocate(spec(0:max(nx, nz)))
            allocate(rkx(0:nx - 1))
            allocate(hrkx(nx))
            allocate(rkz(nz))
            allocate(xtrig(2 * nx))
            allocate(ztrig(2 * nz))
            allocate(kmag(0:nx - 1, 0:nz))
        end subroutine alloc_arrays

        subroutine dealloc_arrays
            deallocate(pp)
            deallocate(ss)
            deallocate(spec)
            deallocate(rkx)
            deallocate(hrkx)
            deallocate(rkz)
            deallocate(xtrig)
            deallocate(ztrig)
            deallocate(kmag)
        end subroutine dealloc_arrays

        subroutine read_data
            double precision, allocatable :: buffer_2d(:, :)
            integer(hid_t)                :: h5handle, group
            character(:), allocatable     :: grn

            call open_h5_file(trim(filename), H5F_ACC_RDONLY_F, h5handle)

            grn = trim(get_step_group_name(step))

            call open_h5_group(h5handle, grn, group)

            if (has_dataset(group, trim(dset))) then
                call read_h5_dataset_2d(group, trim(dset), buffer_2d)
                call fill_field_from_buffer_2d(buffer_2d, pp)
                deallocate(buffer_2d)
            else
                print *, "Error: No dataset '" // trim(dset) // "' in the file."
                stop
            endif

            call close_h5_group(group)
            call close_h5_file(h5handle)
        end subroutine read_data

        ! After reading the H5 dataset into the buffer, copy
        ! the data to a field container
        ! @pre field and buffer must be of rank 2
        subroutine fill_field_from_buffer_2d(buffer, field)
            double precision, allocatable :: buffer(:, :)
            double precision              :: field(0:nz, 0:nx-1)
            integer                       :: dims(2), bdims(2), i, j

            dims = (/nz+1, nx/)

            bdims = shape(buffer)
            if (.not. sum(dims - bdims) == 0) then
                print "(a32, i4, a1, i4, a6, i4, a1, i4, a1)", &
                      "Field dimensions do not agree: (", dims(1), ",", &
                      dims(2), ") != (", bdims(1), ",", bdims(2), ")"
                stop
            endif

            do j = 0, nz
                do i = 0, nx-1
                    field(j, i) = buffer(j, i)
                enddo
            enddo
        end subroutine fill_field_from_buffer_2d

        subroutine write_spectrum
            logical                   :: exists = .false.
            character(:), allocatable :: fname
            integer                   :: pos, kx, kz

            ! 1 October 2021
            ! https://stackoverflow.com/questions/36731707/fortran-how-to-remove-file-extension-from-character
            pos = scan(trim(filename), '.', back=.true.)

            if (pos > 0) then
                fname = filename(1:pos-1) // '_spectrum.asc'
            else
                print *, "Error in reading the filename. File extension not found."
            endif

            inquire(file=fname, exist=exists)
            if (exists) then
                print *, "Error: File '" // trim(fname) // "' already exists."
                stop
            else
                open(unit=1235, file=fname, status='replace')
                write(1235, *) '# The power spectrum of the ' // trim(dset) // ' field.'
                write(1235, *) '# The first column is the wavenumber, the second column the spectrum.'
                write(1235, *) '#         k   P(k)'
            endif

            do kz = 0, nz
                do kx = 0, nx - 1
                    k = kmag(kx, kz)
                    write(1235, *) k * delk, spec(k)
                enddo
            enddo

            close(1235)


            !*** For testing (delete after):
            !Compute P(k)*dk:
            pnorms = delk * sum(spec(0:kmax))
            write(*,*) ' Integral of P(k) dk = ', pnorms
            write(*,*)
            write(*,*) ' Ratio of spectral/physical integrals = ', pnorms / pnorm
            write(*,*)
            write(*,*) ' delk = ', delk

        end subroutine write_spectrum


        ! Get the file name provided via the command line
        subroutine parse_command_line
            integer            :: i, stat
            character(len=512) :: arg
            logical            :: exists

            step = -1
            filename = ''
            dset = 'total buoyancy'
            i = 0
            do
                call get_command_argument(i, arg)
                if (len_trim(arg) == 0) then
                    exit
                endif

                if (arg == '--filename') then
                    i = i + 1
                    call get_command_argument(i, arg)
                    filename = trim(arg)
                else if (arg == '--dset') then
                    i = i + 1
                    call get_command_argument(i, arg)
                    dset = trim(arg)
                else if (arg == '--step') then
                    i = i + 1
                    call get_command_argument(i, arg)

                    ! 1 October 2021
                    ! https://stackoverflow.com/questions/24071722/converting-a-string-to-an-integer-in-fortran-90
                    read(arg, *, iostat=stat)  step
                    if (stat .ne. 0) then
                        print *, 'Error conversion failed.'
                        stop
                    endif
                else if (arg == '--help') then
                    print *, 'This program computes the power spectrum and writes it to file.'
                    print *, 'An EPIC field output must be provided with the step number to analyse.'
                    print *, 'Run code with "genspec --filename [field file]" --step [step number] --dset [dataset]'
                    stop
                endif
                i = i+1
            enddo

            if ((filename == '') .or. (step == -1) ) then
                print *, 'No file or step provided. Run code with "genspec --help"'
                stop
            endif

            ! check if correct file is passed
            stat = index(trim(filename), '_fields.hdf5', back=.true.)
            if (stat == 0) then
                print *, "Error: No EPIC field output file provided."
                stop
            endif

            ! check if file exsits
            inquire(file=trim(filename), exist=exists)
            if (.not. exists) then
                print *, "Error: File '" // trim(filename) // "' does not exist."
                stop
            endif
        end subroutine parse_command_line

end program genspec
