module sta3dfft
    use constants, only : zero, f12
    use stafft
    use sta2dfft
    use deriv1d
    implicit none

    ! Spectral filter:
    double precision, allocatable :: filt(:, :, :)

    ! Wavenumbers::
    double precision, allocatable :: rkx(:), hrkx(:), rky(:), hrky(:), rkz(:)

    ! Quantities needed in FFTs:
    double precision, allocatable :: xtrig(:), ytrig(:), ztrig(:)
    integer                       :: xfactors(5), yfactors(5), zfactors(5)


    double precision, allocatable :: skx(:), sky(:), skz(:)

    private :: xtrig, ytrig, xfactors, yfactors,   &
               rkx, hrkx, rky, hrky, rkz

    contains

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
        subroutine init3dfft(nx, ny, nz, extent)
            integer,          intent(in) :: nx, ny, nz
            double precision, intent(in) :: extent(3)
            double precision             :: kxmaxi, kymaxi, kzmaxi
!             double precision             :: kmax, kc !cfilt,
            integer                      :: nwx, nwy, kx, ky, kz

            if (.not. allocated(filt)) then
                allocate(filt(0:nz, nx, ny))
                allocate(rkx(nx))
                allocate(hrkx(nx))
                allocate(rky(ny))
                allocate(hrky(ny))
                allocate(rkz(nz))
                allocate(xtrig(2*nx))
                allocate(ytrig(2*ny))
                allocate(ztrig(2*nz))
                allocate(skx(nx))
                allocate(sky(ny))
                allocate(skz(0:nz))
            endif

            nwx = nx / 2
            nwy = ny / 2

            !----------------------------------------------------------
            ! Set up FFTs:
            ! Initialise FFTs and wavenumber arrays:
            call init2dfft(nx, ny, extent(1), extent(2), xfactors, yfactors, xtrig, ytrig, hrkx, hrky)
            call initfft(nz, zfactors, ztrig)

            !Define x wavenumbers:
            rkx(1) = zero
            do kx = 1, nwx-1
                rkx(kx+1)    = hrkx(2*kx)
                rkx(nx+1-kx) = hrkx(2*kx)
            enddo
            rkx(nwx+1) = hrkx(nx)

            !Define y wavenumbers:
            rky(1) = zero
            do ky = 1, nwy-1
                rky(ky+1)    = hrky(2*ky)
                rky(ny+1-ky) = hrky(2*ky)
            enddo
            rky(nwy+1) = hrky(ny)

            !Define z wavenumbers:
            call init_deriv(nz, extent(3), rkz)


            !----------------------------------------------------------
            !Define Hou and Li filter:
            kxmaxi = one / maxval(rkx)
            skx = -36.d0 * (kxmaxi * rkx) ** 36
            kymaxi = one/maxval(rky)
            sky = -36.d0 * (kymaxi * rky) ** 36
            kzmaxi = one / maxval(rkz)
            skz(0) = zero
            skz(1:nz)=-36.d0 * (kzmaxi * rkz) ** 36
            do ky = 1, ny
                do kx = 1, nx
                    do kz = 0, nz
                        filt(kz, kx, ky) = dexp(skx(kx) + sky(ky) + skz(kz))
                    enddo
                enddo
            enddo

!             !Define filter:
!
!             !kz > 0:
!             do ky = 1, ny
!                 do kx = 1, nx
!                     do kz = 1, nz
!                         filt(kz, kx, ky) = rkx(kx)**2 + rky(ky)**2 + rkz(kz)**2
!                     enddo
!                 enddo
!             enddo
!
!             !kz = 0:
!             do ky = 1, ny
!                 do kx = 1, nx
!                     filt(0, kx, ky) = rkx(kx)**2 + rky(ky)**2
!                 enddo
!             enddo
!
!             kmax = dsqrt(maxval(filt))
! !             write(*,*)
! !             write(*,'(a,f9.3)') '  Note, kx_max = ',rkx(nwx+1)
! !             write(*,'(a,f9.3)') '        ky_max = ',rky(nwy+1)
! !             write(*,'(a,f9.3)') '        kz_max = ',rkz(nz)
! !             write(*,'(a,f9.3)') '   and  k_max = ',kmax
! !             write(*,*)
! !             write(*,*) ' Enter k_c/k_max:'
! !             read(*,*) kc
!
!             kc =  0.827133988d0 !2.0d0/3.0d0 !0.389d0  !FIXME
!             kc = kc * kmax
! !             cfilt = -one / kc**2
!
!             print *, "kmax = ", kmax
!             print *, "kc * kmax = ", kc
!
!              filt = 0.0d0
!             do ky = 1, ny
!                 do kx = 1, nx
!                     do kz = 0, nz
!                         if (filt(kz, kx, ky) < kc ** 2) then
!                              filt(kz, kx, ky) = 1.0d0
!                          endif
! !                          filt(kz, kx, ky) = dexp(cfilt * filt(kz, kx, ky))
!                     enddo
!                 enddo
!             enddo
! !             print *, "max. cfilt*filt", maxval(filt * filt)
        end subroutine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Computes a 3D FFT of an array fp in physical space and
        !returns the result as fs in spectral space.  It is assumed that
        !fp is generally non-zero at the z boundaries (so a cosine
        !transform is used in z).
        !*** fp is destroyed upon exit ***
        subroutine fftczp2s(fp, fs, nx, ny, nz)
            integer,          intent(in)    :: nx, ny, nz
            double precision, intent(inout) :: fp(0:nz, ny, nx)  !Physical
            double precision, intent(out)   :: fs(0:nz, nx, ny)  !Spectral
            integer                         :: kx, ky, iy

            !Carry out a full x transform first:
            call forfft((nz+1)*ny, nx, fp, xtrig, xfactors)

            !Transpose array:
            do kx = 1, nx
                do iy = 1, ny
                    fs(:, kx, iy) = fp(:, iy, kx)
                enddo
            enddo

            !Carry out a full y transform on transposed array:
            call forfft((nz+1)*nx, ny, fs, ytrig, yfactors)

            !Carry out z FFT for each kx and ky:
            do ky = 1, ny
                do kx = 1, nx
                    call dct(1, nz, fs(:, kx, ky), ztrig, zfactors)
                enddo
            enddo
        end subroutine

        !::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

        !Computes an *inverse* 3D FFT of an array fs in spectral space and
        !returns the result as fp in physical space.  It is assumed that
        !fp is generally non-zero at the z boundaries (so a cosine
        !transform is used in z).
        !*** fs is destroyed upon exit ***
        subroutine fftczs2p(fs, fp, nx, ny, nz)
            integer,          intent(in)    :: nx, ny, nz
            double precision, intent(inout) :: fs(0:nz, nx, ny)  !Spectral
            double precision, intent(out)   :: fp(0:nz, ny, nx)  !Physical
            integer                         :: ix, iy, kx

            !Carry out a full inverse y transform first:
            call revfft((nz+1)*nx,ny,fs,ytrig,yfactors)

            !Transpose array:
            do kx = 1, nx
                do iy = 1, ny
                    fp(:, iy, kx) = fs(:, kx, iy)
                enddo
            enddo

            !Carry out a full inverse x transform:
            call revfft((nz+1)*ny, nx, fp, xtrig, xfactors)

            !Carry out z FFT for each ix and iy:
            do ix = 1, nx
                do iy = 1, ny
                    call dct(1, nz, fp(:, iy, ix), ztrig, zfactors)
                enddo
            enddo
        end subroutine

end module sta3dfft
