! =============================================================================
! Module to compute the eigenvalues and eigenvectors of 3x3 real symmetric
! matrices. It uses the qr algorithm of
! https://doi.org/10.1142/S0129183108012303
! (or https://arxiv.org/abs/physics/0610206)
! =============================================================================
module qr
    use constants
    implicit none
    external dsyevq3

    private
        integer, parameter :: n = 3
        logical :: do_evec

    public :: qr_diagonalise


    contains

    subroutine apply_rotation(A, V, i, j)
            double precision, intent(inout)           :: A(n, n)
            double precision, intent(inout), optional :: V(n, n)
            integer,          intent(in)              :: i, j
            double precision                          :: theta, c, s, t, tau
            integer                                   :: k
            double precision                          :: g, h, aij

            ! compute the rotation angle theta
            ! Reference:    Rutishauser, H. The Jacobi method for real symmetric matrices.
            !               Numer. Math. 9, 1-10 (1966). https://doi.org/10.1007/BF02165223
            g = hundred * dabs(A(i, j))
            h = A(j, j) - A(i, i)
            if (dabs(h) + g == dabs(h)) then
                t = A(i, j) / h
            else
                theta = f12 * h / A(i, j)
                t = one / (dabs(theta) + dsqrt(one + theta ** 2))
                if (theta < zero) then
                    t = -t
                endif
            endif

            ! c = cos(theta)
            c = one / dsqrt(one + t ** 2)

            ! s = sin(theta)
            s = t * c

            tau = s / (one + c)

            ! accumulate eigenvector
            if (do_evec) then
                do k = 1, n
                    g = V(i, k)
                    h = V(j, k)
                    V(i, k) = c * g - s * h
                    V(j, k) = s * g + c * h
                enddo
            endif

            !
            ! Apply Givens rotation to matrix
            !

            aij = A(i, j)

            ! update off-diagonal entries
            do k = 1, n
                g = A(k, i)
                h = A(k, j)
                if (.not. (k == i)) then
                    A(k, i) = g - s * (h + tau * g)
                    A(i, k) = A(k, i)
                endif

                if (.not. (k == j)) then
                    A(k, j) = h + s * (g - tau * h)
                    A(j, k) = A(k, j)
                endif
            enddo

            ! update diagonal entries
            A(i, i) = A(i, i) - t * aij
            A(j, j) = A(j, j) + t * aij

            ! set A(i, j) = A(j, i) explicitly to zero
            A(i, j) = zero
            A(j, i) = zero

        end subroutine apply_rotation

        subroutine givens(xi, xj, c, s)
            double precision, intent(in)    :: xi, xj
            double precision, intent(out)   :: c, s
            double precision                :: r


            r = one / dsqrt(xi ** 2 + xj ** 2)

            c = xi * r
            s = -xj * r

!             double precision                :: theta, t
!             integer                         :: k
!             double precision                :: g, h



!             g = hundred * dabs(A(i, j))
!             h = A(j, j) - A(i, i)
!             if (dabs(h) + g == dabs(h)) then
!                 t = A(i, j) / h
!             else
!                 theta = f12 * h / A(i, j)
!                 t = one / (dabs(theta) + dsqrt(one + theta ** 2))
!                 if (theta < zero) then
!                     t = -t
!                 endif
!             endif
!
!             ! c = cos(theta)
!             c = one / dsqrt(one + t ** 2)
!
!             ! s = sin(theta)
!             s = t * c
        end subroutine givens

        ! Diagonalise a real symmetric 3x3 matrix A = V^T * D * V
        ! where D is a diagonal matrix with the eigenvalues and
        ! V is the eigenvector matrix. The eigenvalues are in
        ! descending order.
        ! The input matrix is overwritten and will be the diagonal
        ! matrix storing the eigenvalues, V contains the eigenvectors.
        ! The eigenvector V(:, i) belongs to eigenvalue A(i, i).
        ! @param[inout] A real symmetric 3x3 matrix
        ! @param[out] V eigenvector matrix (optional)
        subroutine qr_diagonalise(A, V)
            double precision, intent(inout)           :: A(3, 3)
            double precision, intent(out), optional   :: V(3, 3)
            double precision                          :: P(3, 3), G(2, 2), rho, c, s
            double precision                          :: rtol, a11, U(3, 3), sig
            integer                                   :: k, i, m
            double precision                          :: aa(3), bb(3), d, x, y, w, z
            double precision :: uu(2), xx(2), rn, cc(2), ss(2), ww(3)

            do_evec = .true. !present(V)


!             call dsyevq3(A, P, ww)

!             print *, ww
!             print *, P
!             stop

!             ! Householder transform to bring A to Hessenberg form
!             P = zero
!             P(1, 1) = one
!             P(2, 2) = one
!             P(3, 3) = one
!
!             rho = sign(one, A(2, 1))
!
!             xx = A(2:3, 1)
!
!             rn = rho * norm2(xx)
!
!             uu(1) = xx(1) - rn
!             uu(2) = xx(2)
!
!             uu = uu / norm2(xx - rn * (/1.0d0, 0.0d0/))
!
!             P(2, 2:3) = P(2, 2:3) - two * uu(1) * uu
!             P(3, 2:3) = P(3, 2:3) - two * uu(2) * uu
!
! !                 print *, P(1, :)
! !                 print *, P(2, :)
! !                 print *, P(3, :)
!
            P = zero
            P(1, 1) = one
            P(2, 2) = one
            P(3, 3) = one
            call apply_rotation(A, P, 3, 1)
!             A = matmul(matmul(P, A), P)
!             U = matmul(matmul(P, U), P)

            U(1, :) = P(:, 1)
            U(2, :) = P(:, 2)
            U(3, :) = P(:, 3)
            print *, A(1, :)
            print *, A(2, :)
            print *, A(3, :)
            print *, ""
            print *, U(1, :)
            print *, U(2, :)
            print *, U(3, :)
!             stop
!             stop

            print *, "--------------------"
            print *, "Hessenberg"

            ! A tridiagonal

            print *, A(1, :)
            print *, A(2, :)
            print *, A(3, :)

            m = n

            aa(1) = A(1, 1)
            aa(2) = A(2, 2)
            aa(3) = A(3, 3)
            bb(2) = A(2, 1)
            bb(3) = A(3, 2)

            do while (m > 1)
                d = (aa(m-1) - aa(m)) / two
                if (d == 0) then
                    s = aa(m) - dabs(bb(m))
                else
                    s = aa(m) - bb(m) ** 2 / (d + sign(1.0d0, d) * dsqrt(d ** 2 + bb(m) ** 2))
                endif

                x = aa(1) - s
                y = bb(2)

                do k = 1, m-1
                    if (m > 2) then
                        call givens(x, y, c, s)
                    else
                        call givens(aa(1), bb(2), c, s)
                    endif

                    w = c * x - s * y
                    d = aa(k) - aa(k+1)
                    z = (two * c * bb(k+1) + d * s ) * s
                    aa(k) = aa(k) - z
                    aa(k+1) = aa(k+1) + z
                    bb(k+1) = d * c * s + (c ** 2 - s ** 2) * bb(k+1)
                    x = bb(k+1)

                    if (k > 1) then
                        bb(k) = w
                    endif

                    if (k < m - 1) then
                        y = -s * bb(k+2)
                        bb(k+2) = c * bb(k+2)
                    endif
                    G(1, 1) = c
                    G(1, 2) = s
                    G(2, 1) = -s
                    G(2, 2) = c
                    U(1:n,k:k+1) = matmul(U(1:n,k:k+1), G)
                enddo

!                 stop

                if (dabs(bb(m)) < 1.0e-15 * (dabs(aa(m-1)) + dabs(aa(m)))) then
                    m = m - 1
                endif
            enddo

!                 print *, "aa", aa
!                 print *, "bb", bb


!             rtol = 1000
!
!             do while (rtol > 1.0e-15)
!
!                     print *, "rtol", rtol
!
!                     a11 = A(1, 1)
!
!                     do k = 1, 2
!                         call givens(A(k, k), A(k+1, k), cc(k), ss(k))
!
!     !                     print *, "c", cc(k)
!     !                     print *, "s", ss(k)
!
!
!                         G(1, 1) = cc(k)
!                         G(1, 2) = -ss(k)
!                         G(2, 1) = ss(k)
!                         G(2, 2) = cc(k)
!
!                         A(k:k+1,k:3) = matmul(G, A(k:k+1, k:3))
!
!
!
!                     enddo
!
!                     do k = 1, 2
!
!                         G(1, 1) = cc(k)
!                         G(1, 2) = ss(k)
!                         G(2, 1) = -ss(k)
!                         G(2, 2) = cc(k)
!
!                         A(1:k+1, k:k+1) = matmul(A(1:k+1, k:k+1), G)
!
!                         U(1:k+1,k:k+1) = matmul(U(1:k+1,k:k+1), G)
!                     enddo
!
!                     A(1, 1) = A(1, 1) + sig
!                     A(2, 2) = A(2, 2) + sig
!                     A(3, 3) = A(3, 3) + sig
!
!     !                     print *, "--------------------"
!     !                     print *, A(1, :)
!     !                     print *, A(2, :)
!     !                     print *, A(3, :)
!
!                     rtol = (A(1, 1) - a11) / a11
!             enddo

!             print *, A(1, :)
!             print *, A(2, :)
!             print *, A(3, :)
!             print *, "-------------------"
            print *, aa
            print *, U(:, 1)
            print *, U(:, 2)
            print *, U(:, 3)

!             call givens(A, 2, 3, c, s)
!
!             G = zero
!             G(1, 1) = one
!             G(2, 2) = one
!             G(3, 3) = one
!
!
!             G(2, 2) = c
!             G(2, 3) = s
!             G(3, 2) = -s
!             G(3, 3) = c
!
!             A = matmul(G, A)
!
!             print *, "--------------------"
!             print *, A(1, :)
!             print *, A(2, :)
!             print *, A(3, :)

!             print *, c, s


!             call sort_descending(A, V)


        end subroutine qr_diagonalise

        subroutine sort_descending(D, V)
            double precision, intent(inout)           :: D(n, n)
            double precision, intent(inout), optional :: V(n, n)
            double precision                          :: teval, tevec(n)

            if (D(2, 2) > D(1, 1)) then
                teval = D(1, 1)
                D(1, 1) = D(2, 2)
                D(2, 2) = teval
                if (do_evec) then
                    tevec = V(1, :)
                    V(1, :) = V(2, :)
                    V(2, :) = tevec
                endif
            endif

            if (D(3, 3) > D(2, 2)) then
                teval = D(2, 2)
                D(2, 2) = D(3, 3)
                D(3, 3) = teval
                if (do_evec) then
                    tevec = V(2, :)
                    V(2, :) = V(3, :)
                    V(3, :) = tevec
                endif
            endif

            if (D(2, 2) > D(1, 1)) then
                teval = D(1, 1)
                D(1, 1) = D(2, 2)
                D(2, 2) = teval
                if (do_evec) then
                    tevec = V(1, :)
                    V(1, :) = V(2, :)
                    V(2, :) = tevec
                endif
            endif
        end subroutine sort_descending

end module qr
