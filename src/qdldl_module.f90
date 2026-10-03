!*****************************************************************************************
!>
!  The public interface of the QDLDL library: [[qdldl_type]], a sparse \(LDL^T\)
!  solver for quasi-definite matrices, with a fill-reducing ordering, and the
!  low-level routines of [[qdldl_core]].
!
!  Typical use:
!
!```fortran
!  type(qdldl_type) :: ldl
!  call ldl%analyze(n, irow, icol, istat)   ! the pattern, once
!  call ldl%factor(val, istat)              ! values in the pattern's order
!  call ldl%solve(b, istat)                 ! b is overwritten by x
!```
!
!  The factorization does not pivot. It is stable for quasi-definite matrices
!  \(\begin{bmatrix} H & B^T \\ B & -C \end{bmatrix}\) with \(H\) and \(C\) positive
!  definite, in any symmetric order. For other matrices it may run, but can be
!  unstable; see the regularization options of [[qdldl_type]].
!
!### License
!
!  Apache License 2.0. This library is a derivative of QDLDL
!  (<https://github.com/osqp/qdldl>, Copyright 2018, Paul Goulart, Bartolomeo
!  Stellato, Goran Banjac, Ian McInerney, The OSQP developers), and includes a
!  modernized version of AMD (BSD-3-clause; see [[qdldl_amd]]).

    module qdldl_module

    use qdldl_kinds,    only: wp, ip
    use qdldl_core,     only: qdldl_unknown, qdldl_success, qdldl_error_not_upper, qdldl_error_overflow, &
                              qdldl_error_empty_column, qdldl_error_zero_pivot, qdldl_error_invalid_input, &
                              qdldl_error_not_analyzed, qdldl_error_not_factored, qdldl_error_out_of_memory, &
                              qdldl_error_not_finite, qdldl_etree, qdldl_factor, qdldl_factor_ext, qdldl_solve, &
                              qdldl_lsolve, qdldl_ltsolve
    use qdldl_ordering, only: qdldl_order_natural, qdldl_order_user, qdldl_order_rcm, qdldl_order_amd, &
                              qdldl_order_default, qdldl_symmetric_pattern, qdldl_order_natural_perm, qdldl_rcm, &
                              qdldl_amd_order, qdldl_is_permutation

    implicit none

    private

    ! kinds
    integer,parameter,public :: qdldl_wp = wp  !! real kind
    integer,parameter,public :: qdldl_ip = ip  !! integer kind of indices and counts

    ! orderings
    public :: qdldl_order_natural, qdldl_order_user, qdldl_order_rcm, qdldl_order_amd, qdldl_order_default

    ! status codes
    public :: qdldl_success, qdldl_error_not_upper, qdldl_error_overflow, qdldl_error_empty_column, &
              qdldl_error_zero_pivot, qdldl_error_invalid_input, qdldl_error_not_analyzed, &
              qdldl_error_not_factored, qdldl_error_out_of_memory, qdldl_error_not_finite
    public :: qdldl_status_message

    ! low-level routines
    public :: qdldl_unknown
    public :: qdldl_etree, qdldl_factor, qdldl_factor_ext, qdldl_solve, qdldl_lsolve, qdldl_ltsolve
    public :: qdldl_rcm, qdldl_amd_order, qdldl_symmetric_pattern, qdldl_is_permutation

    type,public :: qdldl_type

        !! A sparse \(LDL^T\) factorization of a symmetric quasi-definite matrix
        !! \(A\), as \(PAP^T = LDL^T\) with a fill-reducing permutation \(P\).
        !!
        !! Give the pattern once ([[qdldl_type:analyze]] or [[qdldl_type:analyze_csc]]),
        !! then [[qdldl_type:factor]] and [[qdldl_type:solve]] as often as needed.
        !! Neither of those allocates.
        !!
        !! All components are allocatable (no pointers), so intrinsic assignment
        !! (`ldl2 = ldl1`) makes an independent copy of the whole state, including a
        !! factorization.
        !!
        !! The options are public components with defaults; set them before
        !! [[qdldl_type:analyze]] (`ordering`) or [[qdldl_type:factor]] (the others).
        !! The results are public components too, but read-only by convention.

        private

        ! options
        integer,public  :: ordering = qdldl_order_default
            !! the fill-reducing ordering (`qdldl_order_default` is AMD), used by
            !! [[qdldl_type:analyze]] unless a permutation or ordering is passed to it
        real(wp),public :: zero_pivot_tol = 0.0_wp
            !! a pivot with \(|d_k| \le\) `zero_pivot_tol` \(\max_{ij} |A_{ij}|\) is zero.
            !! `0`: only an exact zero (upstream's behaviour)
        logical,public  :: regularize = .false.
            !! dynamic regularization: replace a zero pivot, or one with
            !! \(s_k d_k <\) `reg_eps` (\(s_k\) its expected sign, see
            !! [[qdldl_type:set_signs]]), by \(s_k\) `reg_delta`, and continue.
            !! Off: a zero pivot stops the factorization (upstream's behaviour)
        real(wp),public :: reg_eps = epsilon(1.0_wp)**0.8_wp
            !! threshold of dynamic regularization (absolute)
        real(wp),public :: reg_delta = sqrt(epsilon(1.0_wp))
            !! the magnitude of a dynamically regularized pivot (absolute, positive)
        real(wp),public :: static_reg = 0.0_wp
            !! static regularization: `static_reg` \(s_k\) is added to each diagonal entry
            !! before factoring (\(s_k\): the expected sign of row \(k\), or the sign of
            !! its diagonal entry if unknown, \(+1\) for zero). Iterative refinement
            !! in [[qdldl_type:solve]] is against the unregularized matrix. `0`: off
        integer,public  :: max_refine = 0
            !! the maximum number of steps of iterative refinement in
            !! [[qdldl_type:solve]] (unless its `refine` argument is given). `0`: none

        ! results
        integer(ip),public :: n = 0                  !! order of the matrix
        integer(ip),public :: nnz_a = 0              !! entries of the upper triangle of \(A\) (duplicates merged, diagonal included)
        integer(ip),public :: nnz_l = 0              !! nonzeros of \(L\) below the diagonal (from the analysis)
        integer(ip),public :: n_positive = 0         !! number of positive pivots (as computed, before regularization)
        integer(ip),public :: n_negative = 0         !! number of negative pivots (as computed, before regularization)
        integer(ip),public :: n_zero = 0             !! number of zero pivots (within the tolerance; nonzero only when regularizing)
        integer(ip),public :: n_regularized = 0      !! number of pivots replaced by dynamic regularization
        integer(ip),public :: zero_pivot_column = 0  !! after `qdldl_error_zero_pivot`: the row/column of \(A\) whose pivot was zero
        real(wp),public    :: max_abs_l = 0.0_wp     !! \(\max |L_{ij}|\) (large: an unstable factorization)
        real(wp),public    :: pivot_ratio = 0.0_wp   !! \(\max_k |d_k| / \min_k |d_k|\) (after regularization)
        integer,public     :: refine_steps = 0       !! steps of iterative refinement taken by the last [[qdldl_type:solve]]
        real(wp),public    :: residual = 0.0_wp      !! the scaled residual \(\lVert Ax-b \rVert_\infty /
                                                     !! (\lVert A \rVert_\infty \lVert x \rVert_\infty + \lVert b \rVert_\infty)\)
                                                     !! of the last [[qdldl_type:solve]] with refinement (else 0)

        ! state
        logical :: analyzed = .false.      !! the pattern has been analyzed
        logical :: has_values = .false.    !! values have been given to [[qdldl_type:factor]] (for [[qdldl_type:multiply]])
        logical :: factored = .false.      !! the last [[qdldl_type:factor]] succeeded
        integer(ip) :: nz = 0              !! number of entries given to the analysis (the size of `val`)
        real(wp) :: anorm = 0.0_wp         !! \(\lVert A \rVert_\infty\) of the last values

        ! the matrix: the upper triangle of P A P^T in CSC form
        integer(ip),dimension(:),allocatable :: Ap         !! column pointers (`n+1`)
        integer(ip),dimension(:),allocatable :: Ai         !! row indices (`nnz_a`)
        real(wp),dimension(:),allocatable    :: Ax         !! values (`nnz_a`), unregularized
        integer(ip),dimension(:),allocatable :: map        !! position in `Ax` of each input entry (`nz`)
        integer(ip),dimension(:),allocatable :: diag_pos   !! position in `Ax` of each diagonal entry (`n`)
        integer(ip),dimension(:),allocatable :: perm       !! the permutation: row `perm(k)` of A is row `k` of P A P^T (`n`)
        real(wp),dimension(:),allocatable    :: signs      !! expected sign of each pivot, permuted (`n`; unallocated: unknown)

        ! the factors
        integer(ip),dimension(:),allocatable :: etree      !! elimination tree (`n`)
        integer(ip),dimension(:),allocatable :: Lnz        !! nonzeros in each column of L (`n`)
        integer(ip),dimension(:),allocatable :: Lp         !! column pointers of L (`n+1`)
        integer(ip),dimension(:),allocatable :: Li         !! row indices of L (`nnz_l`)
        real(wp),dimension(:),allocatable    :: Lx         !! values of L (`nnz_l`)
        real(wp),dimension(:),allocatable    :: D          !! D (`n`)
        real(wp),dimension(:),allocatable    :: Dinv       !! 1/D (`n`)

        ! work arrays
        logical,dimension(:),allocatable     :: bwork      !! factorization (`n`)
        integer(ip),dimension(:),allocatable :: iwork      !! factorization (`3n`)
        real(wp),dimension(:),allocatable    :: fwork      !! factorization (`n`)
        real(wp),dimension(:),allocatable    :: wb         !! permuted right-hand side (`n`)
        real(wp),dimension(:),allocatable    :: wx         !! permuted solution (`n`)
        real(wp),dimension(:),allocatable    :: wr         !! residual and correction (`n`)
        real(wp),dimension(:),allocatable    :: wy         !! trial solution (`n`)

    contains

        procedure,public :: analyze
        procedure,public :: analyze_csc
        procedure,public :: set_signs
        procedure,public :: factor
        procedure,public :: solve
        procedure,public :: multiply
        procedure,public :: inertia
        procedure,public :: get_permutation
        procedure,public :: destroy
        procedure,public :: is_analyzed
        procedure,public :: is_factored

    end type qdldl_type

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  A short description of a status code.

    pure function qdldl_status_message(istat) result(msg)

    integer,intent(in)            :: istat  !! a status code
    character(len=:),allocatable  :: msg    !! its description

    select case (istat)
    case (qdldl_success);             msg = 'success'
    case (qdldl_error_not_upper);     msg = 'an entry below the diagonal in upper-triangular input'
    case (qdldl_error_overflow);      msg = 'the nonzeros overflow the integer kind (build with INT64)'
    case (qdldl_error_empty_column);  msg = 'a column of the matrix is empty'
    case (qdldl_error_zero_pivot);    msg = 'zero pivot (the matrix is singular, or not quasi-definite)'
    case (qdldl_error_invalid_input); msg = 'invalid input'
    case (qdldl_error_not_analyzed);  msg = 'the pattern has not been analyzed'
    case (qdldl_error_not_factored);  msg = 'the matrix has not been factored'
    case (qdldl_error_out_of_memory); msg = 'out of memory'
    case (qdldl_error_not_finite);    msg = 'a value is not finite'
    case default;                     msg = 'unknown status'
    end select

    end function qdldl_status_message
!*****************************************************************************************

!*****************************************************************************************
!>
!  Analyzes the sparsity pattern of the matrix, given in coordinate form: computes
!  the ordering, the elimination tree, and the nonzeros of \(L\), and allocates
!  everything that [[qdldl_type:factor]] and [[qdldl_type:solve]] need.
!
!  The entries `(irow(k), icol(k))` may be in either triangle, or in both: an
!  off-diagonal entry \((i,j)\) stands for both \((i,j)\) and \((j,i)\), and
!  duplicates (including an entry given in both triangles) are **added
!  together**. So for a matrix given by both triangles, give the values of only
!  one of them (and zeros for the other), or each off-diagonal value is doubled.
!  Missing diagonal entries are added as explicit zeros.
!
!  The ordering is: `perm` if given (row `perm(k)` of \(A\) becomes row `k`), else
!  the `ordering` argument if given, else the `ordering` component.
!
!  Any previous analysis and factorization are discarded (the options are kept,
!  but the expected signs of [[qdldl_type:set_signs]] must be set again).

    subroutine analyze(me, n, irow, icol, istat, perm, ordering)

    class(qdldl_type),intent(inout)   :: me
    integer(ip),intent(in)            :: n         !! order of the matrix (`n >= 0`)
    integer(ip),intent(in)            :: irow(:)   !! row indices of the entries (`1..n`)
    integer(ip),intent(in)            :: icol(:)   !! column indices of the entries (`1..n`, same size as `irow`)
    integer,intent(out)               :: istat     !! `qdldl_success`, `qdldl_error_invalid_input`,
                                                   !! `qdldl_error_out_of_memory`, or `qdldl_error_overflow`
    integer(ip),intent(in),optional   :: perm(:)   !! a permutation of `1..n` (the ordering `qdldl_order_user`)
    integer,intent(in),optional       :: ordering  !! the ordering (`qdldl_order_*`); overrides the component

    integer :: order
    integer(ip) :: k

    order = me%ordering
    if (present(ordering)) order = ordering
    if (present(perm)) order = qdldl_order_user

    call me%destroy()

    ! check the input
    if (n < 0 .or. n >= huge(n) - 2_ip .or. size(irow, kind=ip) /= size(icol, kind=ip)) then
        istat = qdldl_error_invalid_input
        return
    end if
    if (size(irow, kind=ip) > huge(n) - n) then
        istat = qdldl_error_overflow
        return
    end if
    do k = 1, size(irow, kind=ip)
        if (irow(k) < 1 .or. irow(k) > n .or. icol(k) < 1 .or. icol(k) > n) then
            istat = qdldl_error_invalid_input
            return
        end if
    end do
    select case (order)
    case (qdldl_order_default, qdldl_order_natural, qdldl_order_rcm, qdldl_order_amd)
    case (qdldl_order_user)
        if (.not. present(perm)) then
            istat = qdldl_error_invalid_input
            return
        end if
        if (.not. qdldl_is_permutation(n, perm)) then
            istat = qdldl_error_invalid_input
            return
        end if
    case default
        istat = qdldl_error_invalid_input
        return
    end select

    call analyze_pattern(me, n, irow, icol, order, istat, perm)
    if (istat /= qdldl_success) call me%destroy()

    end subroutine analyze
!*****************************************************************************************

!*****************************************************************************************
!>
!  [[qdldl_type:analyze]] for a matrix given in compressed sparse column (CSC) form:
!  the entries of column `j` are `Ai(Ap(j):Ap(j+1)-1)`, with `Ap(1) = 1`. Normally the
!  upper triangle (as for the low-level routines), but, as for coordinate input,
!  entries of either triangle are accepted and duplicates are added. The values
!  given to [[qdldl_type:factor]] are then in the order of `Ai`.

    subroutine analyze_csc(me, n, Ap, Ai, istat, perm, ordering)

    class(qdldl_type),intent(inout)   :: me
    integer(ip),intent(in)            :: n         !! order of the matrix (`n >= 0`)
    integer(ip),intent(in)            :: Ap(:)     !! column pointers (size `n+1`, `Ap(1) = 1`, nondecreasing)
    integer(ip),intent(in)            :: Ai(:)     !! row indices (size at least `Ap(n+1)-1`)
    integer,intent(out)               :: istat     !! as for [[qdldl_type:analyze]]
    integer(ip),intent(in),optional   :: perm(:)   !! a permutation of `1..n`
    integer,intent(in),optional       :: ordering  !! the ordering (`qdldl_order_*`)

    integer(ip),allocatable :: icol(:), irow(:)
    integer(ip) :: j, p, nnz
    integer :: stat

    call me%destroy()
    if (n < 0 .or. size(Ap, kind=ip) /= n + 1) then
        istat = qdldl_error_invalid_input
        return
    end if
    if (Ap(1) /= 1) then
        istat = qdldl_error_invalid_input
        return
    end if
    do j = 1, n
        if (Ap(j+1) < Ap(j)) then
            istat = qdldl_error_invalid_input
            return
        end if
    end do
    nnz = Ap(n+1) - 1
    if (size(Ai, kind=ip) < nnz) then
        istat = qdldl_error_invalid_input
        return
    end if

    allocate(irow(nnz), icol(nnz), stat=stat)
    if (stat /= 0) then
        istat = qdldl_error_out_of_memory
        return
    end if
    do j = 1, n
        do p = Ap(j), Ap(j+1) - 1
            irow(p) = Ai(p)
            icol(p) = j
        end do
    end do
    call me%analyze(n, irow, icol, istat, perm, ordering)

    end subroutine analyze_csc
!*****************************************************************************************

!*****************************************************************************************
!>
!  The analysis proper, after the input has been checked.

    subroutine analyze_pattern(me, n, irow, icol, order, istat, user_perm)

    class(qdldl_type),intent(inout)   :: me
    integer(ip),intent(in)            :: n             !! order of the matrix
    integer(ip),intent(in)            :: irow(:)       !! row indices
    integer(ip),intent(in)            :: icol(:)       !! column indices
    integer,intent(in)                :: order         !! the ordering
    integer,intent(out)               :: istat         !! status code
    integer(ip),intent(in),optional   :: user_perm(:)  !! the user's permutation

    integer(ip),allocatable :: Bp(:), Bi(:), map1(:), map2(:), bucket(:), next(:), marker(:), iperm(:)
    integer(ip),allocatable :: xadj(:), adj(:)
    integer(ip) :: nz, ntot, k, i, j, c, r, p, q, start, nnzb, sumLnz
    integer :: stat

    istat = qdldl_success
    nz = size(irow, kind=ip)
    ntot = nz + n            ! the entries, then a (zero) diagonal entry for each column

    allocate(Bp(n+1), next(n+1), bucket(ntot), Bi(ntot), map1(ntot), marker(n), stat=stat)
    if (stat /= 0) then
        istat = qdldl_error_out_of_memory
        return
    end if

    !---------------------------------------------------------------------------
    ! the upper triangle of A in CSC form, duplicates merged, with every diagonal
    !---------------------------------------------------------------------------

    ! bucket the entries by column of the upper triangle (counting sort)
    Bp = 0
    do k = 1, ntot
        call entry(k, r, c)
        Bp(c) = Bp(c) + 1
    end do
    q = 1
    do j = 1, n
        p = Bp(j)
        Bp(j) = q
        next(j) = q
        q = q + p
    end do
    Bp(n+1) = q
    do k = 1, ntot
        call entry(k, r, c)
        bucket(next(c)) = k
        next(c) = next(c) + 1
    end do

    ! merge duplicates within each column: map1(k) is the position of entry k
    marker = 0
    q = 1
    do j = 1, n
        start = q
        do p = Bp(j), Bp(j+1) - 1
            k = bucket(p)
            call entry(k, r, c)
            if (marker(r) >= start) then
                map1(k) = marker(r)     ! a duplicate
            else
                marker(r) = q
                Bi(q) = r
                map1(k) = q
                q = q + 1
            end if
        end do
        Bp(j) = start
    end do
    Bp(n+1) = q
    nnzb = q - 1
    deallocate(bucket, marker)

    !---------------------------------------------------------------------------
    ! the ordering
    !---------------------------------------------------------------------------

    allocate(me%perm(n), iperm(n), stat=stat)
    if (stat /= 0) then
        istat = qdldl_error_out_of_memory
        return
    end if
    select case (order)
    case (qdldl_order_natural)
        call qdldl_order_natural_perm(n, me%perm)
    case (qdldl_order_user)
        me%perm = user_perm
    case (qdldl_order_rcm, qdldl_order_amd, qdldl_order_default)
        call qdldl_symmetric_pattern(n, Bp, Bi, xadj, adj, istat)
        if (istat /= qdldl_success) return
        if (order == qdldl_order_rcm) then
            call qdldl_rcm(n, xadj, adj, me%perm, istat)
        else
            call qdldl_amd_order(n, xadj, adj, me%perm, istat)
        end if
        if (istat /= qdldl_success) return
        deallocate(xadj, adj)
    case default
        istat = qdldl_error_invalid_input   ! (checked by the caller)
        return
    end select
    do k = 1, n
        iperm(me%perm(k)) = k
    end do

    !---------------------------------------------------------------------------
    ! the upper triangle of P A P^T in CSC form
    !---------------------------------------------------------------------------

    allocate(me%Ap(n+1), me%Ai(max(1_ip, nnzb)), me%Ax(max(1_ip, nnzb)), map2(max(1_ip, nnzb)), &
             me%map(max(1_ip, nz)), me%diag_pos(n), stat=stat)
    if (stat /= 0) then
        istat = qdldl_error_out_of_memory
        return
    end if

    me%Ap = 0
    do j = 1, n
        do p = Bp(j), Bp(j+1) - 1
            c = max(iperm(Bi(p)), iperm(j))
            me%Ap(c) = me%Ap(c) + 1
        end do
    end do
    q = 1
    do j = 1, n
        p = me%Ap(j)
        me%Ap(j) = q
        next(j) = q
        q = q + p
    end do
    me%Ap(n+1) = q
    do j = 1, n
        do p = Bp(j), Bp(j+1) - 1
            i = iperm(Bi(p))
            c = max(i, iperm(j))
            r = min(i, iperm(j))
            me%Ai(next(c)) = r
            map2(p) = next(c)
            if (r == c) me%diag_pos(c) = next(c)
            next(c) = next(c) + 1
        end do
    end do
    do k = 1, nz
        me%map(k) = map2(map1(k))
    end do
    deallocate(Bp, Bi, map1, map2, next, iperm)

    !---------------------------------------------------------------------------
    ! the elimination tree, and the factors' storage
    !---------------------------------------------------------------------------

    allocate(me%etree(n), me%Lnz(n), me%Lp(n+1), me%D(n), me%Dinv(n), &
             me%bwork(n), me%iwork(3*n), me%fwork(n), &
             me%wb(n), me%wx(n), me%wr(n), me%wy(n), stat=stat)
    if (stat /= 0) then
        istat = qdldl_error_out_of_memory
        return
    end if

    sumLnz = qdldl_etree(n, me%Ap, me%Ai, me%iwork, me%Lnz, me%etree)
    if (sumLnz < 0) then
        istat = int(sumLnz)     ! only overflow is possible here (the input has every diagonal)
        return
    end if

    allocate(me%Li(max(1_ip, sumLnz)), me%Lx(max(1_ip, sumLnz)), stat=stat)
    if (stat /= 0) then
        istat = qdldl_error_out_of_memory
        return
    end if

    me%n        = n
    me%nz       = nz
    me%nnz_a    = nnzb
    me%nnz_l    = sumLnz
    me%analyzed = .true.

    contains

        subroutine entry(kk, row, col)
        !! Entry `kk` in the upper triangle: an input entry, or (after them) a diagonal.
        integer(ip),intent(in)  :: kk    !! entry number (`1..nz+n`)
        integer(ip),intent(out) :: row   !! its row in the upper triangle
        integer(ip),intent(out) :: col   !! its column in the upper triangle
        if (kk <= nz) then
            row = min(irow(kk), icol(kk))
            col = max(irow(kk), icol(kk))
        else
            row = kk - nz
            col = row
        end if
        end subroutine entry

    end subroutine analyze_pattern
!*****************************************************************************************

!*****************************************************************************************
!>
!  Sets the expected sign of each pivot (in the original order of the rows of
!  \(A\)), for regularization: `+1` for a row of the positive definite block (\(H\)),
!  `-1` for one of the negative definite block (\(-C\)), `0` if unknown (the sign is
!  then inferred from the values). Must be called after [[qdldl_type:analyze]] (which
!  discards it).

    subroutine set_signs(me, signs, istat)

    class(qdldl_type),intent(inout) :: me
    integer,intent(in)              :: signs(:)  !! expected signs (`-1`, `0`, or `+1`; size `n`)
    integer,intent(out)             :: istat     !! `qdldl_success`, `qdldl_error_not_analyzed`,
                                                 !! `qdldl_error_invalid_input`, or `qdldl_error_out_of_memory`

    integer(ip) :: k
    integer :: stat

    if (.not. me%analyzed) then
        istat = qdldl_error_not_analyzed
        return
    end if
    if (size(signs, kind=ip) /= me%n) then
        istat = qdldl_error_invalid_input
        return
    end if
    do k = 1, me%n
        if (abs(signs(k)) > 1) then
            istat = qdldl_error_invalid_input
            return
        end if
    end do
    if (.not. allocated(me%signs)) then
        allocate(me%signs(me%n), stat=stat)
        if (stat /= 0) then
            istat = qdldl_error_out_of_memory
            return
        end if
    end if
    do k = 1, me%n
        me%signs(k) = real(signs(me%perm(k)), wp)
    end do
    istat = qdldl_success

    end subroutine set_signs
!*****************************************************************************************

!*****************************************************************************************
!>
!  Factors the matrix with new values, in the order of the entries given to
!  [[qdldl_type:analyze]] (duplicates are added). Allocates nothing.
!
!  On success, the inertia (`n_positive`, `n_negative`, `n_zero`), `n_regularized`,
!  `max_abs_l`, and `pivot_ratio` are set. On a zero pivot, `zero_pivot_column`
!  is the row of \(A\) whose pivot was zero.

    subroutine factor(me, val, istat)

    class(qdldl_type),intent(inout) :: me
    real(wp),intent(in)             :: val(:)  !! the values of the entries (same size as `irow`)
    integer,intent(out)             :: istat   !! `qdldl_success`, `qdldl_error_not_analyzed`,
                                               !! `qdldl_error_invalid_input` (wrong size, or
                                               !! `regularize` with `reg_delta <= 0`),
                                               !! `qdldl_error_not_finite`, or `qdldl_error_zero_pivot`

    integer(ip) :: k, j, p, i, zcol
    real(wp) :: amax, tau, s, a, dmin, dmax

    me%factored          = .false.
    me%n_positive        = 0
    me%n_negative        = 0
    me%n_zero            = 0
    me%n_regularized     = 0
    me%zero_pivot_column = 0
    me%max_abs_l         = 0.0_wp
    me%pivot_ratio       = 0.0_wp

    if (.not. me%analyzed) then
        istat = qdldl_error_not_analyzed
        return
    end if
    if (size(val, kind=ip) /= me%nz .or. (me%regularize .and. .not. me%reg_delta > 0.0_wp)) then
        istat = qdldl_error_invalid_input
        return
    end if
    do k = 1, me%nz
        if (.not. is_finite(val(k))) then
            istat = qdldl_error_not_finite
            return
        end if
    end do

    ! scatter the values into the permuted upper triangle
    do p = 1, me%nnz_a
        me%Ax(p) = 0.0_wp
    end do
    do k = 1, me%nz
        me%Ax(me%map(k)) = me%Ax(me%map(k)) + val(k)
    end do
    me%has_values = .true.

    ! the largest entry, and the infinity norm (row sums, in fwork)
    amax = 0.0_wp
    do k = 1, me%n
        me%fwork(k) = 0.0_wp
    end do
    do j = 1, me%n
        do p = me%Ap(j), me%Ap(j+1) - 1
            i = me%Ai(p)
            a = abs(me%Ax(p))
            amax = max(amax, a)
            me%fwork(i) = me%fwork(i) + a
            if (i /= j) me%fwork(j) = me%fwork(j) + a
        end do
    end do
    me%anorm = 0.0_wp
    do k = 1, me%n
        me%anorm = max(me%anorm, me%fwork(k))
    end do
    tau = me%zero_pivot_tol * amax

    ! static regularization (the diagonal is saved in wy, and restored after)
    if (me%static_reg > 0.0_wp) then
        do k = 1, me%n
            p = me%diag_pos(k)
            me%wy(k) = me%Ax(p)
            s = 0.0_wp
            if (allocated(me%signs)) s = me%signs(k)
            if (s == 0.0_wp) then
                if (me%Ax(p) < 0.0_wp) then
                    s = -1.0_wp
                else
                    s = 1.0_wp
                end if
            end if
            me%Ax(p) = me%Ax(p) + s * me%static_reg
        end do
    end if

    istat = qdldl_factor_ext(me%n, me%Ap, me%Ai, me%Ax, me%Lp, me%Li, me%Lx, me%D, me%Dinv, &
                             me%Lnz, me%etree, me%bwork, me%iwork, me%fwork, &
                             tau, me%regularize, me%reg_eps, me%reg_delta, &
                             me%n_positive, me%n_negative, me%n_zero, me%n_regularized, zcol, me%signs)

    if (me%static_reg > 0.0_wp) then
        do k = 1, me%n
            me%Ax(me%diag_pos(k)) = me%wy(k)
        end do
    end if

    if (istat /= qdldl_success) then
        me%zero_pivot_column = me%perm(zcol)
        return
    end if

    ! check the factors, and the growth statistics
    dmin = huge(1.0_wp)
    dmax = 0.0_wp
    do k = 1, me%n
        if (.not. is_finite(me%D(k))) then
            istat = qdldl_error_not_finite
            return
        end if
        dmin = min(dmin, abs(me%D(k)))
        dmax = max(dmax, abs(me%D(k)))
    end do
    do p = 1, me%nnz_l
        if (.not. is_finite(me%Lx(p))) then
            istat = qdldl_error_not_finite
            return
        end if
        me%max_abs_l = max(me%max_abs_l, abs(me%Lx(p)))
    end do
    if (me%n > 0) me%pivot_ratio = dmax / dmin

    me%factored = .true.

    end subroutine factor
!*****************************************************************************************

!*****************************************************************************************
!>
!  Solves \(Ax = b\) in place with the factorization, with optional iterative
!  refinement against the (unregularized) matrix: steps of \(r = b - Ax\),
!  \(x \leftarrow x + \text{solve}(r)\), while the scaled residual is above
!  round-off and each step at least halves the residual (a step that doesn't
!  reduce it is discarded). Allocates nothing.

    subroutine solve(me, b, istat, refine)

    class(qdldl_type),intent(inout) :: me
    real(wp),intent(inout)          :: b(:)    !! on input the right-hand side, on output the solution (size `n`)
    integer,intent(out)             :: istat   !! `qdldl_success`, `qdldl_error_not_factored`,
                                               !! `qdldl_error_invalid_input` (wrong size), or
                                               !! `qdldl_error_not_finite` (in `b` or the solution;
                                               !! `b` is then unchanged)
    integer,intent(in),optional     :: refine  !! maximum number of refinement steps (default: `max_refine`)

    integer(ip) :: k, n
    integer :: nref, it
    real(wp) :: rn, rnew, bn, tol

    me%refine_steps = 0
    me%residual = 0.0_wp
    if (.not. me%factored) then
        istat = qdldl_error_not_factored
        return
    end if
    n = me%n
    if (size(b, kind=ip) /= n) then
        istat = qdldl_error_invalid_input
        return
    end if
    do k = 1, n
        if (.not. is_finite(b(k))) then
            istat = qdldl_error_not_finite
            return
        end if
    end do
    nref = me%max_refine
    if (present(refine)) nref = refine

    ! permute, and solve
    do k = 1, n
        me%wb(k) = b(me%perm(k))
        me%wx(k) = me%wb(k)
    end do
    call qdldl_solve(n, me%Lp, me%Li, me%Lx, me%Dinv, me%wx)

    if (nref > 0) then
        bn = max_abs(me%wb)
        call residual(me%wx, me%wr, rn)
        do it = 1, nref
            tol = epsilon(1.0_wp) * (me%anorm * max_abs(me%wx) + bn)
            if (rn <= tol .or. .not. is_finite(rn)) exit
            ! the correction, in wr, and the trial solution, in wy
            call qdldl_solve(n, me%Lp, me%Li, me%Lx, me%Dinv, me%wr)
            do k = 1, n
                me%wy(k) = me%wx(k) + me%wr(k)
            end do
            call residual(me%wy, me%wr, rnew)
            if (.not. (rnew < rn)) exit   ! no better: keep x
            do k = 1, n
                me%wx(k) = me%wy(k)
            end do
            me%refine_steps = it
            if (rnew > 0.5_wp * rn) then
                rn = rnew
                exit     ! stagnating
            end if
            rn = rnew
        end do
        tol = me%anorm * max_abs(me%wx) + bn
        if (tol > 0.0_wp) me%residual = rn / tol
    end if

    do k = 1, n
        if (.not. is_finite(me%wx(k))) then
            istat = qdldl_error_not_finite
            return
        end if
    end do
    do k = 1, n
        b(me%perm(k)) = me%wx(k)
    end do
    istat = qdldl_success

    contains

        subroutine residual(x, r, rnorm)
        !! \(r = b - Ax\) (permuted), and its infinity norm.
        real(wp),intent(in)  :: x(me%n)  !! the (permuted) solution
        real(wp),intent(out) :: r(me%n)  !! the residual
        real(wp),intent(out) :: rnorm   !! its infinity norm
        integer(ip) :: kk
        call sym_multiply(me%n, me%Ap, me%Ai, me%Ax, x, r)
        do kk = 1, me%n
            r(kk) = me%wb(kk) - r(kk)
        end do
        rnorm = max_abs(r)
        end subroutine residual

    end subroutine solve
!*****************************************************************************************

!*****************************************************************************************
!>
!  \(y = Ax\) with the last values given to [[qdldl_type:factor]] (unregularized),
!  whether or not that factorization succeeded.

    subroutine multiply(me, x, y, istat)

    class(qdldl_type),intent(in)  :: me
    real(wp),intent(in)           :: x(:)   !! the vector (size `n`)
    real(wp),intent(out)          :: y(:)   !! the product (size `n`)
    integer,intent(out),optional  :: istat  !! `qdldl_success`, `qdldl_error_not_factored` (no values yet),
                                            !! or `qdldl_error_invalid_input` (wrong sizes); `y = 0` on error

    integer(ip) :: j, p, i, pi, pj
    real(wp) :: a

    y = 0.0_wp
    if (.not. me%has_values) then
        if (present(istat)) istat = qdldl_error_not_factored
        return
    end if
    if (size(x, kind=ip) /= me%n .or. size(y, kind=ip) /= me%n) then
        if (present(istat)) istat = qdldl_error_invalid_input
        return
    end if
    do j = 1, me%n
        pj = me%perm(j)
        do p = me%Ap(j), me%Ap(j+1) - 1
            i = me%Ai(p)
            pi = me%perm(i)
            a = me%Ax(p)
            y(pi) = y(pi) + a * x(pj)
            if (i /= j) y(pj) = y(pj) + a * x(pi)
        end do
    end do
    if (present(istat)) istat = qdldl_success

    end subroutine multiply
!*****************************************************************************************

!*****************************************************************************************
!>
!  The inertia of the last successful factorization: the numbers of positive,
!  negative, and zero pivots (as computed, before any regularization).
!
!  By Sylvester's law of inertia, these are the numbers of positive, negative, and
!  zero eigenvalues of \(A\), exactly, if \(A\) is quasi-definite and no pivot was
!  regularized (`n_regularized = 0`). Otherwise they are only as reliable as the
!  pivots: they describe the matrix that was factored, after any regularization of
!  earlier pivots. All zero if the matrix has not been factored.

    subroutine inertia(me, n_positive, n_negative, n_zero, n_regularized)

    class(qdldl_type),intent(in)       :: me
    integer(ip),intent(out)            :: n_positive     !! number of positive pivots
    integer(ip),intent(out)            :: n_negative     !! number of negative pivots
    integer(ip),intent(out)            :: n_zero         !! number of zero pivots (within the tolerance)
    integer(ip),intent(out),optional   :: n_regularized  !! number of pivots replaced by dynamic regularization

    if (me%factored) then
        n_positive = me%n_positive
        n_negative = me%n_negative
        n_zero     = me%n_zero
        if (present(n_regularized)) n_regularized = me%n_regularized
    else
        n_positive = 0
        n_negative = 0
        n_zero     = 0
        if (present(n_regularized)) n_regularized = 0
    end if

    end subroutine inertia
!*****************************************************************************************

!*****************************************************************************************
!>
!  The permutation of the analysis: row `perm(k)` of \(A\) is row `k` of \(PAP^T\).
!  Unallocated if the pattern has not been analyzed.

    subroutine get_permutation(me, perm)

    class(qdldl_type),intent(in)                 :: me
    integer(ip),dimension(:),allocatable,intent(out) :: perm  !! the permutation (size `n`)

    if (me%analyzed) perm = me%perm

    end subroutine get_permutation
!*****************************************************************************************

!*****************************************************************************************
!>
!  Whether the pattern has been analyzed.

    pure logical function is_analyzed(me)

    class(qdldl_type),intent(in) :: me

    is_analyzed = me%analyzed

    end function is_analyzed
!*****************************************************************************************

!*****************************************************************************************
!>
!  Whether the last [[qdldl_type:factor]] succeeded (so [[qdldl_type:solve]] can be used).

    pure logical function is_factored(me)

    class(qdldl_type),intent(in) :: me

    is_factored = me%factored

    end function is_factored
!*****************************************************************************************

!*****************************************************************************************
!>
!  Frees all the storage, and resets the results. The options are kept.

    subroutine destroy(me)

    class(qdldl_type),intent(inout) :: me

    if (allocated(me%Ap))       deallocate(me%Ap)
    if (allocated(me%Ai))       deallocate(me%Ai)
    if (allocated(me%Ax))       deallocate(me%Ax)
    if (allocated(me%map))      deallocate(me%map)
    if (allocated(me%diag_pos)) deallocate(me%diag_pos)
    if (allocated(me%perm))     deallocate(me%perm)
    if (allocated(me%signs))    deallocate(me%signs)
    if (allocated(me%etree))    deallocate(me%etree)
    if (allocated(me%Lnz))      deallocate(me%Lnz)
    if (allocated(me%Lp))       deallocate(me%Lp)
    if (allocated(me%Li))       deallocate(me%Li)
    if (allocated(me%Lx))       deallocate(me%Lx)
    if (allocated(me%D))        deallocate(me%D)
    if (allocated(me%Dinv))     deallocate(me%Dinv)
    if (allocated(me%bwork))    deallocate(me%bwork)
    if (allocated(me%iwork))    deallocate(me%iwork)
    if (allocated(me%fwork))    deallocate(me%fwork)
    if (allocated(me%wb))       deallocate(me%wb)
    if (allocated(me%wx))       deallocate(me%wx)
    if (allocated(me%wr))       deallocate(me%wr)
    if (allocated(me%wy))       deallocate(me%wy)

    me%analyzed          = .false.
    me%has_values        = .false.
    me%factored          = .false.
    me%nz                = 0
    me%anorm             = 0.0_wp
    me%n                 = 0
    me%nnz_a             = 0
    me%nnz_l             = 0
    me%n_positive        = 0
    me%n_negative        = 0
    me%n_zero            = 0
    me%n_regularized     = 0
    me%zero_pivot_column = 0
    me%max_abs_l         = 0.0_wp
    me%pivot_ratio       = 0.0_wp
    me%refine_steps      = 0
    me%residual          = 0.0_wp

    end subroutine destroy
!*****************************************************************************************

!*****************************************************************************************
!>
!  \(y = Ax\) for a symmetric matrix given by its upper triangle in CSC form.

    pure subroutine sym_multiply(n, Ap, Ai, Ax, x, y)

    integer(ip),intent(in) :: n              !! order of the matrix
    integer(ip),intent(in) :: Ap(n+1)        !! column pointers
    integer(ip),intent(in) :: Ai(Ap(n+1)-1)  !! row indices
    real(wp),intent(in)    :: Ax(Ap(n+1)-1)  !! values
    real(wp),intent(in)    :: x(n)           !! the vector
    real(wp),intent(out)   :: y(n)           !! the product

    integer(ip) :: j, p, i

    do j = 1, n
        y(j) = 0.0_wp
    end do
    do j = 1, n
        do p = Ap(j), Ap(j+1) - 1
            i = Ai(p)
            y(i) = y(i) + Ax(p) * x(j)
            if (i /= j) y(j) = y(j) + Ax(p) * x(i)
        end do
    end do

    end subroutine sym_multiply
!*****************************************************************************************

!*****************************************************************************************
!>
!  \(\max_i |x_i|\) (0 for an empty vector).

    pure function max_abs(x) result(m)

    real(wp),intent(in) :: x(:)  !! the vector
    real(wp)            :: m     !! its infinity norm

    integer(ip) :: i

    m = 0.0_wp
    do i = 1, size(x, kind=ip)
        m = max(m, abs(x(i)))
    end do

    end function max_abs
!*****************************************************************************************

!*****************************************************************************************
!>
!  Whether a real is finite (neither infinite nor NaN).

    elemental logical function is_finite(x)

    real(wp),intent(in) :: x  !! the value

    is_finite = abs(x) <= huge(x)

    end function is_finite
!*****************************************************************************************

!*****************************************************************************************
    end module qdldl_module
!*****************************************************************************************
