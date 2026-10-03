!*****************************************************************************************
!>
!  Stage 2: the object-oriented interface, [[qdldl_type]]: upstream's matrices
!  through the type, coordinate input (duplicates, both triangles, missing
!  diagonals), CSC input, refactoring, copies, `multiply`, and the status codes.

    program test_type

    use qdldl_kinds, only: wp, ip
    use qdldl_module
    use qdldl_test_utils

    implicit none

    real(wp),parameter :: res_tol = 100.0_wp * epsilon(1.0_wp)   !! scaled residual tolerance

    write(*,'(A)') 'test_type'

    call test_upstream_matrices()
    call test_coordinate_input()
    call test_refactor_and_copy()
    call test_multiply()
    call test_status_codes()
    call test_empty()

    call finish_tests('test_type')

    contains

    !*************************************************************************************
    !>
    !  Upstream's solvable matrices (CSC, 0-based in upstream) through the type, with
    !  every ordering; the singular one gives a zero pivot in the right column.

    subroutine test_upstream_matrices()

    integer :: o
    integer,parameter :: orders(3) = [qdldl_order_natural, qdldl_order_rcm, qdldl_order_amd]
    character(len=*),parameter :: oname(3) = ['natural', 'rcm    ', 'amd    ']

    write(*,'(A)') ' upstream matrices'
    do o = 1, 3
        call csc_case('basic/'//trim(oname(o)), orders(o), 10_ip, &
            [0,1,2,4,5,6,8,10,12,14,17], [0,1,1,2,3,4,1,5,0,6,3,7,6,8,1,2,9], &
            [1.0_wp, 0.460641_wp, -0.121189_wp, 0.417928_wp, 0.177828_wp, 0.1_wp, &
             -0.0290058_wp, -1.0_wp, 0.350321_wp, -0.441092_wp, -0.0845395_wp, -0.316228_wp, &
             0.178663_wp, -0.299077_wp, 0.182452_wp, -1.56506_wp, -0.1_wp], &
            [1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp, 7.0_wp, 8.0_wp, 9.0_wp, 10.0_wp])
        call csc_case('osqp_kkt/'//trim(oname(o)), orders(o), 7_ip, &
            [0,1,2,5,6,7,8,12], [0,1,2,1,0,3,4,5,5,6,4,3], &
            [-0.25_wp, -0.25_wp, 1.0_wp, 0.513578_wp, 0.529142_wp, -0.25_wp, &
             -0.25_wp, 1.10274_wp, 0.15538_wp, 1.25883_wp, 0.13458_wp, 0.621134_wp], &
            [-0.595598_wp, -0.0193715_wp, -0.576156_wp, -0.168746_wp, 0.61543_wp, 0.419073_wp, 1.31087_wp])
        call csc_case('two_by_two/'//trim(oname(o)), orders(o), 2_ip, [0,1,3], [0,0,1], &
            [1.0_wp, 1.0_wp, -1.0_wp], [2.0_wp, 4.0_wp])
        call csc_case('identity/'//trim(oname(o)), orders(o), 4_ip, [0,1,2,3,4], [0,1,2,3], &
            [1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp], [2.0_wp, 2.0_wp, 2.0_wp, 2.0_wp])
        call csc_case('singleton/'//trim(oname(o)), orders(o), 1_ip, [0,1], [0], [0.2_wp], [2.0_wp])
        ! upstream's sym_structure and tril_structure are errors in the low-level
        ! routines; the type accepts them (duplicates added: the off-diagonal is 2)
        call csc_case('sym_structure/'//trim(oname(o)), orders(o), 2_ip, [0,2,4], [0,1,0,1], &
            [5.0_wp, 1.0_wp, 1.0_wp, 5.0_wp], [1.0_wp, 1.0_wp])
        call csc_case('tril_structure/'//trim(oname(o)), orders(o), 2_ip, [0,2,3], [0,1,1], &
            [5.0_wp, 1.0_wp, 5.0_wp], [1.0_wp, 1.0_wp])
    end do
    ! the zero on the diagonal is removed by fill-in only in the natural order
    call csc_case('zero_on_diag/natural', qdldl_order_natural, 3_ip, [0,1,2,5], [0,0,0,1,2], &
        [4.0_wp, 1.0_wp, 2.0_wp, 1.0_wp, -3.0_wp], [6.0_wp, 9.0_wp, 12.0_wp])
    ! rank deficient: [1 1; 1 1]
    call csc_case('rank_deficient', qdldl_order_natural, 2_ip, [0,1,3], [0,0,1], &
        [1.0_wp, 1.0_wp, 1.0_wp], [1.0_wp, 1.0_wp], expected=qdldl_error_zero_pivot, zero_col=2_ip)

    end subroutine test_upstream_matrices
    !*************************************************************************************

    !*************************************************************************************
    !>
    !  One matrix in CSC form (0-based) through [[qdldl_type:analyze_csc]]: the
    !  solution against a dense solve, and the scaled residual.

    subroutine csc_case(name, order, n, Ap0, Ai0, Ax, b, expected, zero_col)

    character(len=*),intent(in)     :: name       !! name of the case
    integer,intent(in)              :: order      !! the ordering
    integer(ip),intent(in)          :: n          !! order of the matrix
    integer,intent(in)              :: Ap0(:)     !! column pointers (0-based)
    integer,intent(in)              :: Ai0(:)     !! row indices (0-based)
    real(wp),intent(in)             :: Ax(:)      !! values
    real(wp),intent(in)             :: b(:)       !! right-hand side
    integer,intent(in),optional     :: expected   !! expected error of the factorization
    integer(ip),intent(in),optional :: zero_col   !! expected column of the zero pivot

    type(qdldl_type) :: ldl
    integer(ip),allocatable :: Ap(:), Ai(:), irow(:), icol(:)
    real(wp),allocatable :: x(:), xd(:), a(:,:)
    integer(ip) :: j, p
    integer :: istat

    Ap = int(Ap0, ip) + 1_ip
    Ai = int(Ai0, ip) + 1_ip
    call ldl%analyze_csc(n, Ap, Ai, istat, ordering=order)
    call check(istat == qdldl_success, name//': analyze_csc')
    call ldl%factor(Ax, istat)
    if (present(expected)) then
        call check(istat == expected, name//': factor error code')
        if (present(zero_col)) call check(ldl%zero_pivot_column == zero_col, name//': zero pivot column')
        return
    end if
    call check(istat == qdldl_success, name//': factor')
    x = b
    call ldl%solve(x, istat)
    call check(istat == qdldl_success, name//': solve')

    allocate(irow(size(Ai)), icol(size(Ai)))
    do j = 1, n
        do p = Ap(j), Ap(j+1) - 1
            irow(p) = Ai(p)
            icol(p) = j
        end do
    end do
    call coo_to_dense(n, irow, icol, Ax, a)
    call dense_solve(a, b, xd)
    call check(max_abs_diff(x, xd) <= 1.0e3_wp * epsilon(1.0_wp) * maxval(abs(xd)), name//': against a dense solve')
    call check(residual_norm(irow, icol, Ax, x, b) < res_tol, name//': scaled residual')

    end subroutine csc_case
    !*************************************************************************************

    !*************************************************************************************
    !>
    !  Coordinate input with entries in both triangles, duplicates (split values),
    !  and a missing diagonal, against the same matrix given cleanly.

    subroutine test_coordinate_input()

    ! a 5x5 KKT matrix [H B^T; B 0] (rows 1-3: H, rows 4-5: constraints, zero (2,2) block)
    integer(ip),parameter :: n = 5
    integer(ip),parameter :: r1(*) = [1, 1, 2, 2, 3, 1, 2, 3, 3]
    integer(ip),parameter :: c1(*) = [1, 2, 2, 3, 3, 4, 4, 5, 4]
    real(wp),parameter :: v1(*) = [4.0_wp, 1.0_wp, 5.0_wp, -1.0_wp, 6.0_wp, 1.0_wp, 2.0_wp, 1.0_wp, -1.0_wp]
    ! the same, with (1,2) in the lower triangle, (2,2) split in two, (3,4) split
    ! between both triangles, (2,3) given twice as (3,2), and (4,4), (5,5) explicit zeros
    integer(ip),parameter :: r2(*) = [1, 2, 2, 2, 3, 3, 3, 1, 2, 3, 3, 4, 4, 5]
    integer(ip),parameter :: c2(*) = [1, 1, 2, 2, 2, 2, 3, 4, 4, 5, 4, 3, 4, 5]
    real(wp),parameter :: v2(*) = [4.0_wp, 1.0_wp, 2.0_wp, 3.0_wp, -0.5_wp, -0.5_wp, 6.0_wp, &
                                   1.0_wp, 2.0_wp, 1.0_wp, -0.25_wp, -0.75_wp, 0.0_wp, 0.0_wp]
    type(qdldl_type) :: ldl1, ldl2
    real(wp) :: b(n), x1(n), x2(n)
    real(wp),allocatable :: a1(:,:), a2(:,:), xd(:)
    integer :: istat

    write(*,'(A)') ' coordinate input'
    call coo_to_dense(n, r1, c1, v1, a1)
    call coo_to_dense(n, r2, c2, v2, a2)
    call check(all(a1 == a2), 'the two inputs describe the same matrix')

    b = [1.0_wp, -2.0_wp, 3.0_wp, 0.5_wp, 1.5_wp]
    call ldl1%analyze(n, r1, c1, istat, ordering=qdldl_order_natural)
    call check(istat == qdldl_success .and. ldl1%nnz_a == 11, 'clean input: analyze (missing diagonals added)')
    call ldl1%factor(v1, istat)
    call check(istat == qdldl_success, 'clean input: factor')
    x1 = b
    call ldl1%solve(x1, istat)

    call ldl2%analyze(n, r2, c2, istat, ordering=qdldl_order_natural)
    call check(istat == qdldl_success .and. ldl2%nnz_a == 11, 'messy input: analyze (duplicates merged)')
    call ldl2%factor(v2, istat)
    call check(istat == qdldl_success, 'messy input: factor')
    x2 = b
    call ldl2%solve(x2, istat)
    call dense_solve(a1, b, xd)
    call check(max_abs_diff(x1, xd) <= 1.0e3_wp * epsilon(1.0_wp) * maxval(abs(xd)), 'clean input: solution')
    call check(max_abs_diff(x2, xd) <= 1.0e3_wp * epsilon(1.0_wp) * maxval(abs(xd)), 'messy input: solution')
    call check(ldl2%n_positive == 3 .and. ldl2%n_negative == 2, 'inertia of the KKT matrix')

    end subroutine test_coordinate_input
    !*************************************************************************************

    !*************************************************************************************
    !>
    !  Refactoring with new values gives what a fresh analysis gives; a copy of the
    !  object solves independently of the original.

    subroutine test_refactor_and_copy()

    type(qdldl_type) :: ldl, fresh, copy
    integer(ip),allocatable :: irow(:), icol(:)
    real(wp),allocatable :: val(:), val2(:), b(:), x(:), xf(:), xc(:)
    integer,allocatable :: s(:)
    integer(ip) :: n, k
    integer :: istat

    write(*,'(A)') ' refactor and copy'
    call lcg_init(7)
    n = 60
    call random_qd(n, 0.08_wp, s, irow, icol, val)
    allocate(val2(size(val)), b(n))
    do k = 1, size(val, kind=ip)
        val2(k) = val(k) * (1.0_wp + 0.3_wp * lcg_uniform())
    end do
    do k = 1, n
        b(k) = lcg_uniform()
    end do

    call ldl%analyze(n, irow, icol, istat)
    call ldl%factor(val, istat)
    call check(istat == qdldl_success, 'first factorization')
    copy = ldl                            ! a copy of the first factorization
    call ldl%factor(val2, istat)          ! refactor the original with new values
    call check(istat == qdldl_success, 'refactorization')
    x = b
    call ldl%solve(x, istat)

    call fresh%analyze(n, irow, icol, istat)
    call fresh%factor(val2, istat)
    xf = b
    call fresh%solve(xf, istat)
    call check(all(x == xf), 'refactor gives the same solution as a fresh analysis')

    xc = b
    call copy%solve(xc, istat)
    call check(istat == qdldl_success, 'the copy solves')
    call check(residual_norm(irow, icol, val, xc, b) < res_tol, 'the copy solves with the first values')
    call check(residual_norm(irow, icol, val2, x, b) < res_tol, 'the original solves with the new values')
    call copy%destroy()
    x = b
    call ldl%solve(x, istat)
    call check(istat == qdldl_success .and. all(x == xf), 'destroying the copy leaves the original intact')

    end subroutine test_refactor_and_copy
    !*************************************************************************************

    !*************************************************************************************
    !>
    !  `multiply` against the product with the coordinate form.

    subroutine test_multiply()

    type(qdldl_type) :: ldl
    integer(ip),allocatable :: irow(:), icol(:)
    real(wp),allocatable :: val(:), x(:), y(:), yref(:)
    integer,allocatable :: s(:)
    integer(ip) :: n, k
    integer :: istat

    write(*,'(A)') ' multiply'
    call lcg_init(11)
    n = 40
    call random_qd(n, 0.1_wp, s, irow, icol, val)
    ! add some entries in the lower triangle and duplicates
    irow = [irow, 3_ip, 7_ip, 7_ip]
    icol = [icol, 1_ip, 2_ip, 2_ip]
    val  = [val, 0.5_wp, -0.25_wp, 0.125_wp]
    allocate(x(n), y(n), yref(n))
    do k = 1, n
        x(k) = lcg_uniform()
    end do
    call ldl%analyze(n, irow, icol, istat)
    call ldl%multiply(x, y, istat)
    call check(istat == qdldl_error_not_factored, 'multiply before factor')
    call ldl%factor(val, istat)
    call ldl%multiply(x, y, istat)
    call coo_multiply(irow, icol, val, x, yref)
    call check(istat == qdldl_success .and. max_abs_diff(y, yref) <= 100*epsilon(1.0_wp)*maxval(abs(yref)), &
               'multiply against the coordinate form')

    end subroutine test_multiply
    !*************************************************************************************

    !*************************************************************************************
    !>
    !  Every status code of the interface.

    subroutine test_status_codes()

    type(qdldl_type) :: ldl
    real(wp) :: b(2), nan, inf
    integer :: istat

    write(*,'(A)') ' status codes'
    nan = 0.0_wp
    nan = nan / nan
    inf = huge(1.0_wp)
    inf = inf * 2

    call ldl%factor([1.0_wp], istat)
    call check(istat == qdldl_error_not_analyzed, 'factor before analyze')
    b = 1.0_wp
    call ldl%solve(b, istat)
    call check(istat == qdldl_error_not_factored, 'solve before factor')
    call ldl%set_signs([1, -1], istat)
    call check(istat == qdldl_error_not_analyzed, 'set_signs before analyze')

    call ldl%analyze(-1_ip, [integer(ip) ::], [integer(ip) ::], istat)
    call check(istat == qdldl_error_invalid_input, 'n < 0')
    call ldl%analyze(2_ip, [1_ip, 3_ip], [1_ip, 2_ip], istat)
    call check(istat == qdldl_error_invalid_input, 'index out of range')
    call ldl%analyze(2_ip, [1_ip, 0_ip], [1_ip, 2_ip], istat)
    call check(istat == qdldl_error_invalid_input, 'index zero')
    call ldl%analyze(2_ip, [1_ip, 2_ip], [1_ip], istat)
    call check(istat == qdldl_error_invalid_input, 'irow and icol of different sizes')
    call ldl%analyze(2_ip, [1_ip, 2_ip], [1_ip, 2_ip], istat, ordering=qdldl_order_user)
    call check(istat == qdldl_error_invalid_input, 'user ordering without perm')
    call ldl%analyze(2_ip, [1_ip, 2_ip], [1_ip, 2_ip], istat, perm=[1_ip, 1_ip])
    call check(istat == qdldl_error_invalid_input, 'perm not a permutation')
    call ldl%analyze(2_ip, [1_ip, 2_ip], [1_ip, 2_ip], istat, ordering=99)
    call check(istat == qdldl_error_invalid_input, 'unknown ordering')
    call ldl%analyze_csc(2_ip, [1_ip, 2_ip], [1_ip], istat)
    call check(istat == qdldl_error_invalid_input, 'CSC: Ap of the wrong size')
    call ldl%analyze_csc(2_ip, [1_ip, 3_ip, 2_ip], [1_ip, 2_ip], istat)
    call check(istat == qdldl_error_invalid_input, 'CSC: Ap decreasing')

    call ldl%analyze(2_ip, [1_ip, 1_ip, 2_ip], [1_ip, 2_ip, 2_ip], istat, perm=[2_ip, 1_ip])
    call check(istat == qdldl_success, 'analyze with a user permutation')
    call ldl%factor([1.0_wp, 2.0_wp], istat)
    call check(istat == qdldl_error_invalid_input, 'factor: val of the wrong size')
    call ldl%factor([1.0_wp, nan, 3.0_wp], istat)
    call check(istat == qdldl_error_not_finite, 'factor: NaN in val')
    call ldl%factor([1.0_wp, inf, 3.0_wp], istat)
    call check(istat == qdldl_error_not_finite, 'factor: Inf in val')
    call ldl%set_signs([1, 2], istat)
    call check(istat == qdldl_error_invalid_input, 'set_signs: a sign of 2')
    call ldl%set_signs([1], istat)
    call check(istat == qdldl_error_invalid_input, 'set_signs: wrong size')
    ldl%regularize = .true.
    ldl%reg_delta = 0.0_wp
    call ldl%factor([1.0_wp, 2.0_wp, 3.0_wp], istat)
    call check(istat == qdldl_error_invalid_input, 'factor: regularize with reg_delta = 0')
    ldl%regularize = .false.
    call ldl%factor([1.0_wp, 2.0_wp, 3.0_wp], istat)
    call check(istat == qdldl_success .and. ldl%is_factored(), 'factor [1 2; 2 3]')
    call ldl%solve(b(1:1), istat)
    call check(istat == qdldl_error_invalid_input, 'solve: b of the wrong size')
    b = [1.0_wp, nan]
    call ldl%solve(b, istat)
    call check(istat == qdldl_error_not_finite .and. b(1) == 1.0_wp, 'solve: NaN in b (b unchanged)')
    b = [1.0_wp, 1.0_wp]
    call ldl%solve(b, istat)
    call check(istat == qdldl_success .and. max_abs_diff(b, [-1.0_wp, 1.0_wp]) <= 10*epsilon(1.0_wp), &
               'solve [1 2; 2 3] x = [1; 1]')
    call ldl%factor([1.0_wp, 2.0_wp, 4.0_wp], istat)
    call check(istat == qdldl_error_zero_pivot .and. .not. ldl%is_factored(), 'zero pivot: not factored')
    call check(ldl%zero_pivot_column == 1, 'zero pivot column (row 1 of A, pivoted second)')
    call ldl%solve(b, istat)
    call check(istat == qdldl_error_not_factored, 'solve after a failed factorization')
    call check(len(qdldl_status_message(qdldl_error_zero_pivot)) > 0, 'status message')

    end subroutine test_status_codes
    !*************************************************************************************

    !*************************************************************************************
    !>
    !  An empty (0 by 0) matrix.

    subroutine test_empty()

    type(qdldl_type) :: ldl
    real(wp) :: b(0)
    integer :: istat

    write(*,'(A)') ' empty matrix'
    call ldl%analyze(0_ip, [integer(ip) ::], [integer(ip) ::], istat)
    call check(istat == qdldl_success, 'analyze n = 0')
    call ldl%factor([real(wp) ::], istat)
    call check(istat == qdldl_success, 'factor n = 0')
    call ldl%solve(b, istat)
    call check(istat == qdldl_success, 'solve n = 0')

    end subroutine test_empty
    !*************************************************************************************

    end program test_type
!*****************************************************************************************
