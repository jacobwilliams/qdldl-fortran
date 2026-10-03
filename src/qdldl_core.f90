!*****************************************************************************************
!>
!  The one-to-one port of QDLDL's C routines: the elimination tree, the numeric
!  \(LDL^T\) factorization, and the triangular solves, for a quasi-definite matrix
!  given by its upper triangle in compressed sparse column (CSC) form.
!
!  Changes from the C code (`qdldl.c`):
!
!  * Indexing is 1-based: column pointers `Ap(1:n+1)` with `Ap(1) = 1`, row indices
!    `1..n`. The "no parent" marker of the elimination tree is `0` (`qdldl_unknown`)
!    instead of `-1`.
!  * Return values are named constants (`qdldl_error_*`). An empty column now has its
!    own code (`qdldl_error_empty_column`); upstream returns `-1` for it too.
!  * The factorization loop starts from column 1 instead of treating the first
!    column separately (column 1 of an upper triangle holds only its diagonal, so the
!    general loop body does exactly what upstream's special case does).
!  * [[qdldl_factor_ext]] adds a zero-pivot tolerance and dynamic regularization of
!    the pivots (as in the Rust port used by Clarabel). [[qdldl_factor]] keeps
!    upstream's behaviour.
!
!### License
!
!  This file is a derivative of QDLDL (<https://github.com/osqp/qdldl>),
!  Copyright 2018, Paul Goulart, Bartolomeo Stellato, Goran Banjac, Ian McInerney,
!  The OSQP developers. Licensed under the Apache License, Version 2.0
!  (SPDX-License-Identifier: Apache-2.0). Changed from the original: translated to
!  Fortran, as described above.

    module qdldl_core

    use qdldl_kinds, only: wp, ip

    implicit none

    private

    integer(ip),parameter,public :: qdldl_unknown = 0_ip  !! "No parent" in the elimination tree (`QDLDL_UNKNOWN`)

    ! status codes (`0` is success; every error is negative)
    integer,parameter,public :: qdldl_success             =  0  !! Success
    integer,parameter,public :: qdldl_error_not_upper     = -1  !! An entry below the diagonal in input that must be upper triangular
    integer,parameter,public :: qdldl_error_overflow      = -2  !! The nonzeros of \(L\) overflow the integer kind (build with `INT64`)
    integer,parameter,public :: qdldl_error_empty_column  = -3  !! A column of the matrix has no entries
    integer,parameter,public :: qdldl_error_zero_pivot    = -4  !! A pivot is zero (or within the zero-pivot tolerance)
    integer,parameter,public :: qdldl_error_invalid_input = -5  !! Invalid input (an index out of range, `n < 0`, wrong sizes, ...)
    integer,parameter,public :: qdldl_error_not_analyzed  = -6  !! The pattern has not been analyzed
    integer,parameter,public :: qdldl_error_not_factored  = -7  !! The matrix has not been factored (successfully)
    integer,parameter,public :: qdldl_error_out_of_memory = -8  !! An allocation failed
    integer,parameter,public :: qdldl_error_not_finite    = -9  !! A value (matrix, right-hand side, or solution) is not finite

    public :: qdldl_etree
    public :: qdldl_factor
    public :: qdldl_factor_ext
    public :: qdldl_solve
    public :: qdldl_lsolve
    public :: qdldl_ltsolve

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Computes the elimination tree of a quasi-definite matrix, and the number of
!  nonzeros of each column of \(L\) (strictly below the diagonal), from the upper
!  triangle of the matrix in CSC form, without duplicate entries.
!
!  A column with entries but a zero on its diagonal is not an error here (the
!  factorization may still succeed), but a missing diagonal entry makes a zero
!  pivot unless fill-in reaches that position.

    function qdldl_etree(n, Ap, Ai, work, Lnz, etree) result(sumLnz)

    integer(ip),intent(in)  :: n              !! number of columns of `A` (assumed square)
    integer(ip),intent(in)  :: Ap(n+1)        !! column pointers of `A` (`Ap(1) = 1`)
    integer(ip),intent(in)  :: Ai(Ap(n+1)-1)  !! row indices of `A`
    integer(ip),intent(out) :: work(n)        !! work array (no meaning on return)
    integer(ip),intent(out) :: Lnz(n)         !! number of nonzeros of each column of `L` below the diagonal
    integer(ip),intent(out) :: etree(n)       !! the elimination tree: the parent of each column, or `qdldl_unknown`
    integer(ip)             :: sumLnz         !! total nonzeros of `L` below the diagonal, or a negative error code:
                                              !! `qdldl_error_not_upper`, `qdldl_error_empty_column`,
                                              !! `qdldl_error_invalid_input` (a row index below 1),
                                              !! or `qdldl_error_overflow` (`sumLnz + 1` would overflow)

    integer(ip) :: i, j, p

    do i = 1, n
        ! zero out Lnz and work, and set all etree values to unknown
        work(i)  = 0
        Lnz(i)   = 0
        etree(i) = qdldl_unknown
        ! abort if a column of A has no entry
        if (Ap(i) == Ap(i+1)) then
            sumLnz = qdldl_error_empty_column
            return
        end if
    end do

    do j = 1, n
        work(j) = j
        do p = Ap(j), Ap(j+1) - 1
            i = Ai(p)
            ! abort if entries on the lower triangle (or out of range)
            if (i > j) then
                sumLnz = qdldl_error_not_upper
                return
            else if (i < 1) then
                sumLnz = qdldl_error_invalid_input
                return
            end if
            do while (work(i) /= j)
                if (etree(i) == qdldl_unknown) etree(i) = j
                Lnz(i) = Lnz(i) + 1     ! nonzeros in this column
                work(i) = j
                i = etree(i)
            end do
        end do
    end do

    ! the total nonzeros in L: the space required for Li and Lx. It must leave room
    ! for Lp(n+1) = sumLnz + 1 (upstream, 0-based, only needs sumLnz itself to fit).
    sumLnz = 0
    do i = 1, n
        if (sumLnz > huge(sumLnz) - 1_ip - Lnz(i)) then
            sumLnz = qdldl_error_overflow
            return
        end if
        sumLnz = sumLnz + Lnz(i)
    end do

    end function qdldl_etree
!*****************************************************************************************

!*****************************************************************************************
!>
!  Computes the \(LDL^T\) factorization of a quasi-definite matrix from the upper
!  triangle of the matrix in CSC form, without duplicate entries: \(L\) in CSC form
!  (unit diagonal not stored), \(D\), and \(D^{-1}\).
!
!  `Li` and `Lx` must have room for the number of nonzeros returned by
!  [[qdldl_etree]], whose `Lnz` and `etree` are inputs here. Nothing is allocated.
!
!  This is upstream's `QDLDL_factor`: it stops at the first pivot that is exactly
!  zero. See [[qdldl_factor_ext]] for a tolerance and regularization.

    function qdldl_factor(n, Ap, Ai, Ax, Lp, Li, Lx, D, Dinv, Lnz, etree, bwork, iwork, fwork, &
                          zero_col) result(positiveValuesInD)

    integer(ip),intent(in)            :: n                  !! number of columns of `L` and `A` (both square)
    integer(ip),intent(in)            :: Ap(n+1)            !! column pointers of `A`
    integer(ip),intent(in)            :: Ai(Ap(n+1)-1)      !! row indices of `A`
    real(wp),intent(in)               :: Ax(Ap(n+1)-1)      !! values of `A`
    integer(ip),intent(out)           :: Lp(n+1)            !! column pointers of `L`
    real(wp),intent(out)              :: D(n)               !! the diagonal factor `D`
    real(wp),intent(out)              :: Dinv(n)            !! `1/D`
    integer(ip),intent(in)            :: Lnz(n)             !! nonzeros of each column of `L`, from [[qdldl_etree]]
    integer(ip),intent(in)            :: etree(n)           !! the elimination tree, from [[qdldl_etree]]
    integer(ip),intent(out)           :: Li(sum(Lnz))       !! row indices of `L` (size: the return value of [[qdldl_etree]])
    real(wp),intent(out)              :: Lx(sum(Lnz))       !! values of `L` (same size as `Li`)
    logical,intent(out)               :: bwork(n)           !! work array
    integer(ip),intent(out)           :: iwork(3*n)         !! work array
    real(wp),intent(out)              :: fwork(n)           !! work array
    integer(ip),intent(out),optional  :: zero_col           !! the column of the zero pivot, if any (else 0)
    integer(ip)                       :: positiveValuesInD  !! the number of positive entries of `D`, or
                                                            !! `qdldl_error_zero_pivot` if a pivot is exactly zero

    integer(ip) :: npos, nneg, nzero, nreg, zcol
    integer :: status

    status = factor_kernel(n, Ap, Ai, Ax, Lp, Li, Lx, D, Dinv, Lnz, etree, &
                           bwork, iwork(1:n), iwork(n+1:2*n), iwork(2*n+1:3*n), fwork, &
                           0.0_wp, .false., 0.0_wp, 0.0_wp, npos, nneg, nzero, nreg, zcol)
    if (status == qdldl_success) then
        positiveValuesInD = npos
    else
        positiveValuesInD = status
    end if
    if (present(zero_col)) zero_col = zcol

    end function qdldl_factor
!*****************************************************************************************

!*****************************************************************************************
!>
!  [[qdldl_factor]] with a zero-pivot tolerance, dynamic regularization, and the
!  inertia of `D`.
!
!  Each pivot \(d_k\) is classified, as computed (before any replacement), as zero
!  (\(|d_k| \le\) `zero_tol`), positive, or negative. Then:
!
!  * Without regularization, a zero pivot stops the factorization
!    (`qdldl_error_zero_pivot`, with its column in `zero_col`), as upstream does
!    (with `zero_tol = 0`, only an exact zero).
!  * With regularization, a pivot with \(s_k d_k <\) `reg_eps`, or a zero pivot, is
!    replaced by \(s_k\) `reg_delta`, where \(s_k\) is the expected sign of the pivot:
!    `signs(k)` if it is given and nonzero, else the sign of \(d_k\) itself
!    (\(+1\) for zero). The factorization then always completes; `n_regularized`
!    counts the replacements.

    function qdldl_factor_ext(n, Ap, Ai, Ax, Lp, Li, Lx, D, Dinv, Lnz, etree, bwork, iwork, fwork, &
                              zero_tol, regularize, reg_eps, reg_delta, &
                              n_positive, n_negative, n_zero, n_regularized, zero_col, signs) result(status)

    integer(ip),intent(in)          :: n              !! number of columns of `L` and `A` (both square)
    integer(ip),intent(in)          :: Ap(n+1)        !! column pointers of `A`
    integer(ip),intent(in)          :: Ai(Ap(n+1)-1)  !! row indices of `A`
    real(wp),intent(in)             :: Ax(Ap(n+1)-1)  !! values of `A`
    integer(ip),intent(out)         :: Lp(n+1)        !! column pointers of `L`
    real(wp),intent(out)            :: D(n)           !! the diagonal factor `D` (after any regularization)
    real(wp),intent(out)            :: Dinv(n)        !! `1/D`
    integer(ip),intent(in)          :: Lnz(n)         !! nonzeros of each column of `L`, from [[qdldl_etree]]
    integer(ip),intent(in)          :: etree(n)       !! the elimination tree, from [[qdldl_etree]]
    integer(ip),intent(out)         :: Li(sum(Lnz))   !! row indices of `L` (size: the return value of [[qdldl_etree]])
    real(wp),intent(out)            :: Lx(sum(Lnz))   !! values of `L` (same size as `Li`)
    logical,intent(out)             :: bwork(n)       !! work array
    integer(ip),intent(out)         :: iwork(3*n)     !! work array
    real(wp),intent(out)            :: fwork(n)       !! work array
    real(wp),intent(in)             :: zero_tol       !! a pivot with `abs(d) <= zero_tol` is zero (`0`: exactly zero only)
    logical,intent(in)              :: regularize     !! replace small or wrong-signed pivots instead of stopping
    real(wp),intent(in)             :: reg_eps        !! a pivot with `s*d < reg_eps` is regularized
    real(wp),intent(in)             :: reg_delta      !! the magnitude of a regularized pivot (must be positive)
    integer(ip),intent(out)         :: n_positive     !! number of positive pivots (as computed)
    integer(ip),intent(out)         :: n_negative     !! number of negative pivots (as computed)
    integer(ip),intent(out)         :: n_zero         !! number of zero pivots (as computed: `abs(d) <= zero_tol`)
    integer(ip),intent(out)         :: n_regularized  !! number of pivots replaced by regularization
    integer(ip),intent(out)         :: zero_col       !! the column of the zero pivot that stopped the factorization (else 0)
    real(wp),intent(in),optional    :: signs(n)       !! expected sign of each pivot (`+1`, `-1`, or `0` for unknown)
    integer                         :: status         !! `qdldl_success` or `qdldl_error_zero_pivot`

    status = factor_kernel(n, Ap, Ai, Ax, Lp, Li, Lx, D, Dinv, Lnz, etree, &
                           bwork, iwork(1:n), iwork(n+1:2*n), iwork(2*n+1:3*n), fwork, &
                           zero_tol, regularize, reg_eps, reg_delta, &
                           n_positive, n_negative, n_zero, n_regularized, zero_col, signs)

    end function qdldl_factor_ext
!*****************************************************************************************

!*****************************************************************************************
!>
!  The numeric factorization shared by [[qdldl_factor]] and [[qdldl_factor_ext]],
!  with the integer work array already partitioned (as upstream does with pointers).

    function factor_kernel(n, Ap, Ai, Ax, Lp, Li, Lx, D, Dinv, Lnz, etree, &
                           yMarkers, yIdx, elimBuffer, LNextSpaceInCol, yVals, &
                           zero_tol, regularize, reg_eps, reg_delta, &
                           npos, nneg, nzero, nreg, zero_col, signs) result(status)

    integer(ip),intent(in)       :: n                   !! number of columns
    integer(ip),intent(in)       :: Ap(n+1)             !! column pointers of `A`
    integer(ip),intent(in)       :: Ai(Ap(n+1)-1)       !! row indices of `A`
    real(wp),intent(in)          :: Ax(Ap(n+1)-1)       !! values of `A`
    integer(ip),intent(out)      :: Lp(n+1)             !! column pointers of `L`
    real(wp),intent(out)         :: D(n)                !! `D`
    real(wp),intent(out)         :: Dinv(n)             !! `1/D`
    integer(ip),intent(in)       :: Lnz(n)              !! nonzeros of each column of `L`
    integer(ip),intent(in)       :: etree(n)            !! the elimination tree
    integer(ip),intent(out)      :: Li(sum(Lnz))        !! row indices of `L`
    real(wp),intent(out)         :: Lx(sum(Lnz))        !! values of `L`
    logical,intent(out)          :: yMarkers(n)         !! which entries of `y` have been visited
    integer(ip),intent(out)      :: yIdx(n)             !! the nonzero pattern of `y` (the current row of `L`)
    integer(ip),intent(out)      :: elimBuffer(n)       !! an elimination path, in reverse
    integer(ip),intent(out)      :: LNextSpaceInCol(n)  !! the next free position in each column of `L`
    real(wp),intent(out)         :: yVals(n)            !! the values of `y`
    real(wp),intent(in)          :: zero_tol            !! zero-pivot tolerance
    logical,intent(in)           :: regularize          !! regularize pivots instead of stopping
    real(wp),intent(in)          :: reg_eps             !! regularization threshold
    real(wp),intent(in)          :: reg_delta           !! regularized pivot magnitude
    integer(ip),intent(out)      :: npos                !! number of positive pivots
    integer(ip),intent(out)      :: nneg                !! number of negative pivots
    integer(ip),intent(out)      :: nzero               !! number of zero pivots
    integer(ip),intent(out)      :: nreg                !! number of regularized pivots
    integer(ip),intent(out)      :: zero_col            !! the column of a zero pivot that stopped the factorization
    real(wp),intent(in),optional :: signs(n)            !! expected signs of the pivots
    integer                      :: status              !! `qdldl_success` or `qdldl_error_zero_pivot`

    integer(ip) :: i, j, k, nnzY, bidx, cidx, nextIdx, nnzE, tmpIdx
    real(wp) :: yVals_cidx, dk, s

    status   = qdldl_success
    npos     = 0
    nneg     = 0
    nzero    = 0
    nreg     = 0
    zero_col = 0

    Lp(1) = 1   ! first column starts at index one

    do i = 1, n
        ! compute L column indices
        Lp(i+1) = Lp(i) + Lnz(i)    ! cumsum, total at the end
        ! set all yIdx to be 'unused' initially. In each column of L, the next
        ! available space to start is just the first space in the column
        yMarkers(i) = .false.
        yVals(i) = 0.0_wp
        D(i) = 0.0_wp
        LNextSpaceInCol(i) = Lp(i)
    end do

    do k = 1, n

        ! For each k, we compute a solution to y = L(1:k-1,1:k-1)\b, where b is the
        ! kth column of A that sits above the diagonal. The solution y is then the
        ! kth row of L, with an implied '1' at the diagonal entry.

        nnzY = 0   ! number of nonzeros in this row of L

        ! this loop determines where nonzeros will go in the kth row of L,
        ! but doesn't compute the actual values
        do i = Ap(k), Ap(k+1) - 1

            bidx = Ai(i)   ! we are working on this element of b

            ! Initialize D(k) as the element of this column corresponding to the
            ! diagonal place. Don't use this element as part of the elimination step
            ! that computes the kth row of L
            if (bidx == k) then
                D(k) = Ax(i)
                cycle
            end if

            yVals(bidx) = Ax(i)   ! initialise y(bidx) = b(bidx)

            ! use the forward elimination tree to figure out which elements
            ! must be eliminated after this element of b
            nextIdx = bidx

            if (.not. yMarkers(nextIdx)) then   ! this y term not already visited

                yMarkers(nextIdx) = .true.  ! I touched this one
                elimBuffer(1) = nextIdx     ! it goes at the start of the current list
                nnzE = 1                    ! length of unvisited elimination path from here

                nextIdx = etree(bidx)

                do while (nextIdx /= qdldl_unknown .and. nextIdx < k)
                    if (yMarkers(nextIdx)) exit
                    yMarkers(nextIdx) = .true.      ! I touched this one
                    nnzE = nnzE + 1                 ! the list is one longer than before
                    elimBuffer(nnzE) = nextIdx      ! it goes in the current list
                    nextIdx = etree(nextIdx)        ! one step further along tree
                end do

                ! now I put the buffered elimination list into
                ! my current ordering in reverse order
                do while (nnzE > 0)
                    nnzY = nnzY + 1
                    yIdx(nnzY) = elimBuffer(nnzE)
                    nnzE = nnzE - 1
                end do

            end if

        end do

        ! this loop places nonzero values in the kth row
        do i = nnzY, 1, -1

            cidx = yIdx(i)    ! which column are we working on?

            ! loop along the elements in this column of L and subtract to solve to y
            tmpIdx = LNextSpaceInCol(cidx)
            yVals_cidx = yVals(cidx)
            do j = Lp(cidx), tmpIdx - 1
                yVals(Li(j)) = yVals(Li(j)) - Lx(j) * yVals_cidx
            end do

            ! now I have the cidx-th element of y = L\b, so compute the corresponding
            ! element of this row of L and put it into the right place
            Li(tmpIdx) = k
            Lx(tmpIdx) = yVals_cidx * Dinv(cidx)

            ! D(k) -= yVals(cidx)*yVals(cidx)*Dinv(cidx)
            D(k) = D(k) - yVals_cidx * Lx(tmpIdx)
            LNextSpaceInCol(cidx) = LNextSpaceInCol(cidx) + 1

            ! reset the y values and markers once I'm done with them
            yVals(cidx) = 0.0_wp
            yMarkers(cidx) = .false.

        end do

        ! Classify the pivot. If it is zero, we can't factor this matrix:
        ! abort, unless regularizing.
        dk = D(k)
        if (abs(dk) <= zero_tol) then
            if (.not. regularize) then
                zero_col = k
                status = qdldl_error_zero_pivot
                return
            end if
            nzero = nzero + 1
        else if (dk > 0.0_wp) then
            npos = npos + 1
        else
            nneg = nneg + 1
        end if

        if (regularize) then
            ! the expected sign of this pivot
            s = 0.0_wp
            if (present(signs)) s = signs(k)
            if (s == 0.0_wp) then
                if (dk < 0.0_wp) then
                    s = -1.0_wp
                else
                    s = 1.0_wp
                end if
            else
                s = sign(1.0_wp, s)
            end if
            if (s*dk < reg_eps .or. abs(dk) <= zero_tol) then
                D(k) = s * reg_delta
                nreg = nreg + 1
            end if
        end if

        ! compute the inverse of the diagonal
        Dinv(k) = 1.0_wp / D(k)

    end do

    end function factor_kernel
!*****************************************************************************************

!*****************************************************************************************
!>
!  Solves \((L+I)x = b\) in place.

    pure subroutine qdldl_lsolve(n, Lp, Li, Lx, x)

    integer(ip),intent(in) :: n              !! number of columns of `L`
    integer(ip),intent(in) :: Lp(n+1)        !! column pointers of `L`
    integer(ip),intent(in) :: Li(Lp(n+1)-1)  !! row indices of `L`
    real(wp),intent(in)    :: Lx(Lp(n+1)-1)  !! values of `L`
    real(wp),intent(inout) :: x(n)           !! on entry `b`, on exit `x`

    integer(ip) :: i, j
    real(wp) :: val

    do i = 1, n
        val = x(i)
        do j = Lp(i), Lp(i+1) - 1
            x(Li(j)) = x(Li(j)) - Lx(j) * val
        end do
    end do

    end subroutine qdldl_lsolve
!*****************************************************************************************

!*****************************************************************************************
!>
!  Solves \((L+I)^T x = b\) in place.

    pure subroutine qdldl_ltsolve(n, Lp, Li, Lx, x)

    integer(ip),intent(in) :: n              !! number of columns of `L`
    integer(ip),intent(in) :: Lp(n+1)        !! column pointers of `L`
    integer(ip),intent(in) :: Li(Lp(n+1)-1)  !! row indices of `L`
    real(wp),intent(in)    :: Lx(Lp(n+1)-1)  !! values of `L`
    real(wp),intent(inout) :: x(n)           !! on entry `b`, on exit `x`

    integer(ip) :: i, j
    real(wp) :: val

    do i = n, 1, -1
        val = x(i)
        do j = Lp(i), Lp(i+1) - 1
            val = val - Lx(j) * x(Li(j))
        end do
        x(i) = val
    end do

    end subroutine qdldl_ltsolve
!*****************************************************************************************

!*****************************************************************************************
!>
!  Solves \(LDL^Tx = b\) in place, given the factors from [[qdldl_factor]].

    pure subroutine qdldl_solve(n, Lp, Li, Lx, Dinv, x)

    integer(ip),intent(in) :: n              !! number of columns of `L`
    integer(ip),intent(in) :: Lp(n+1)        !! column pointers of `L`
    integer(ip),intent(in) :: Li(Lp(n+1)-1)  !! row indices of `L`
    real(wp),intent(in)    :: Lx(Lp(n+1)-1)  !! values of `L`
    real(wp),intent(in)    :: Dinv(n)        !! `1/D`
    real(wp),intent(inout) :: x(n)           !! on entry `b`, on exit `x`

    integer(ip) :: i

    call qdldl_lsolve(n, Lp, Li, Lx, x)
    do i = 1, n
        x(i) = x(i) * Dinv(i)
    end do
    call qdldl_ltsolve(n, Lp, Li, Lx, x)

    end subroutine qdldl_solve
!*****************************************************************************************

!*****************************************************************************************
    end module qdldl_core
!*****************************************************************************************
