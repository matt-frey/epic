! =============================================================================
! Module to compute the eigenvalues and eigenvectors of 3x3 real symmetric
! matrices. It uses the qr algorithm of
! https://doi.org/10.1142/S0129183108012303
! (or https://arxiv.org/abs/physics/0610206)
! =============================================================================
module qr
    use constants
    implicit none

    private
        integer, parameter :: n = 3
        logical :: do_evec

    public :: qr_diagonalise


    contains

    pure subroutine apply_rotation(A, V, i, j)
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
                    g = V(k, i)
                    h = V(k, j)
                    V(k, i) = c * g - s * h
                    V(k, j) = s * g + c * h
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

        pure subroutine givens(xi, xj, c, s)
            double precision, intent(in)    :: xi, xj
            double precision, intent(out)   :: c, s
            double precision                :: r


            r = one / dsqrt(xi ** 2 + xj ** 2)

            c = xi * r
            s = -xj * r

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
            double precision                          :: G(2, 2), b2, c, s
            integer                                   :: k, m
            double precision                          :: bb(3), d, x, y, w, z ! aa(3),

            do_evec = present(V)

            V(1, 1) = one
            V(1, 2) = zero
            V(1, 3) = zero
            V(2, 2) = one
            V(2, 1) = zero
            V(2, 3) = zero
            V(3, 1) = zero
            V(3, 2) = zero
            V(3, 3) = one

            if (dabs(A(2, 1)) + dabs(A(3, 1)) + dabs(A(3, 2)) < 1.0e-12) then
                ! diagonal matrix
                return
            endif

            call apply_rotation(A, V, 3, 1)

            m = n

            bb(2) = A(2, 1)
            bb(3) = A(3, 2)

            do while (m > 1)
                d = f12 * (A(m-1, m-1) - A(m, m))
                if (d == 0) then
                    s = A(m, m) - dabs(bb(m))
                else
                    b2 = bb(m) ** 2
                    s = A(m, m) - b2 / (d + sign(one, d) * dsqrt(d ** 2 + b2))
                endif

                x = A(1, 1) - s
                y = bb(2)

                do k = 1, m-1
                    if (m > 2) then
                        call givens(x, y, c, s)
                    else
                        call givens(A(1, 1), bb(2), c, s)
                    endif

                    w = c * x - s * y
                    d = A(k, k) - A(k+1,k+1)
                    z = (two * c * bb(k+1) + d * s ) * s
                    A(k, k) = A(k, k) - z
                    A(k+1,k+1) = A(k+1,k+1) + z
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
                    V(1:n,k:k+1) = matmul(V(1:n,k:k+1), G)
                enddo

                if (dabs(bb(m)) < 1.0e-15 * (dabs(A(m-1, m-1)) + dabs(A(m, m)))) then
                    m = m - 1
                endif
            enddo

            call sort_descending(A, V)

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
                    tevec = V(:, 1)
                    V(:, 1) = V(:, 2)
                    V(:, 2) = tevec
                endif
            endif

            if (D(3, 3) > D(2, 2)) then
                teval = D(2, 2)
                D(2, 2) = D(3, 3)
                D(3, 3) = teval
                if (do_evec) then
                    tevec = V(:, 2)
                    V(:, 2) = V(:, 3)
                    V(:, 3) = tevec
                endif
            endif

            if (D(2, 2) > D(1, 1)) then
                teval = D(1, 1)
                D(1, 1) = D(2, 2)
                D(2, 2) = teval
                if (do_evec) then
                    tevec = V(:, 1)
                    V(:, 1) = V(:, 2)
                    V(:, 2) = tevec
                endif
            endif
        end subroutine sort_descending

end module qr
