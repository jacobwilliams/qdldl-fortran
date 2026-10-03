!*****************************************************************************************
!>
!  Stage 1: upstream QDLDL's unit tests (`tests/test_*.h`), through the low-level
!  routines of [[qdldl_core]]. The matrices are upstream's, with indices converted
!  to 1-based. Each solvable case checks the solution against upstream's (to
!  upstream's tolerance, 1e-4, widened in single precision) and the scaled residual;
!  each error case checks its error code.

    program test_core

    use qdldl_kinds, only: wp, ip
    use qdldl_core
    use qdldl_test_utils

    implicit none

    real(wp),parameter :: res_tol = 100.0_wp * epsilon(1.0_wp)   !! scaled residual tolerance

    write(*,'(A)') 'test_core'

    ! test_basic
    call run_case('basic', 10_ip, &
        [0,1,2,4,5,6,8,10,12,14,17], &
        [0,1,1,2,3,4,1,5,0,6,3,7,6,8,1,2,9], &
        [1.0_wp, 0.460641_wp, -0.121189_wp, 0.417928_wp, 0.177828_wp, 0.1_wp, &
         -0.0290058_wp, -1.0_wp, 0.350321_wp, -0.441092_wp, -0.0845395_wp, -0.316228_wp, &
         0.178663_wp, -0.299077_wp, 0.182452_wp, -1.56506_wp, -0.1_wp], &
        [1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp, 6.0_wp, 7.0_wp, 8.0_wp, 9.0_wp, 10.0_wp], &
        [10.2171_wp, 3.9416_wp, -5.69096_wp, 9.28661_wp, 50.0_wp, &
         -6.11433_wp, -26.3104_wp, -27.7809_wp, -45.8099_wp, -3.74178_wp])

    ! test_identity
    call run_case('identity', 4_ip, [0,1,2,3,4], [0,1,2,3], &
        [1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp], [2.0_wp, 2.0_wp, 2.0_wp, 2.0_wp], [2.0_wp, 2.0_wp, 2.0_wp, 2.0_wp])

    ! test_rank_deficient: [1 1; 1 1]
    call run_case('rank_deficient', 2_ip, [0,1,3], [0,0,1], [1.0_wp, 1.0_wp, 1.0_wp], &
        [1.0_wp, 1.0_wp], expected=qdldl_error_zero_pivot)

    ! test_singleton
    call run_case('singleton', 1_ip, [0,1], [0], [0.2_wp], [2.0_wp], [10.0_wp])

    ! test_sym_structure: both triangles given
    call run_case('sym_structure', 2_ip, [0,2,4], [0,1,0,1], [5.0_wp, 1.0_wp, 1.0_wp, 5.0_wp], &
        [1.0_wp, 1.0_wp], expected=qdldl_error_not_upper)

    ! test_tril_structure: the lower triangle given
    call run_case('tril_structure', 2_ip, [0,2,3], [0,1,1], [5.0_wp, 1.0_wp, 5.0_wp], &
        [1.0_wp, 1.0_wp], expected=qdldl_error_not_upper)

    ! test_two_by_two
    call run_case('two_by_two', 2_ip, [0,1,3], [0,0,1], [1.0_wp, 1.0_wp, -1.0_wp], &
        [2.0_wp, 4.0_wp], [3.0_wp, -1.0_wp])

    ! test_zero_on_diag: solvable, with a zero on the diagonal that fill-in removes
    call run_case('zero_on_diag', 3_ip, [0,1,2,5], [0,0,0,1,2], &
        [4.0_wp, 1.0_wp, 2.0_wp, 1.0_wp, -3.0_wp], &
        [6.0_wp, 9.0_wp, 12.0_wp], [17.0_wp, -46.0_wp, -8.0_wp])

    ! test_osqp_kkt: row indices unordered within columns
    call run_case('osqp_kkt', 7_ip, [0,1,2,5,6,7,8,12], [0,1,2,1,0,3,4,5,5,6,4,3], &
        [-0.25_wp, -0.25_wp, 1.0_wp, 0.513578_wp, 0.529142_wp, -0.25_wp, &
         -0.25_wp, 1.10274_wp, 0.15538_wp, 1.25883_wp, 0.13458_wp, 0.621134_wp], &
        [-0.595598_wp, -0.0193715_wp, -0.576156_wp, -0.168746_wp, 0.61543_wp, 0.419073_wp, 1.31087_wp], &
        [1.13141_wp, -1.1367_wp, -0.591044_wp, 1.68867_wp, -2.24209_wp, 0.32254_wp, 0.407998_wp])

    ! an empty column (upstream returns -1 for it, as for a lower-triangle entry)
    call run_case('empty_column', 2_ip, [0,1,1], [0], [1.0_wp], [1.0_wp, 1.0_wp], &
        expected=qdldl_error_empty_column)

    call test_separate_solves()
    call test_zero_col()

    call finish_tests('test_core')

    contains

    !*************************************************************************************
    !>
    !  Upstream's `ldl_factor_solve` and the checks of one test: the etree, the
    !  factorization, and the solve, with 0-based input converted to 1-based.

    subroutine run_case(name, n, Ap0, Ai0, Ax, b, xsol, expected)

    character(len=*),intent(in)  :: name          !! the name of the upstream test
    integer(ip),intent(in)       :: n             !! matrix order
    integer,intent(in)           :: Ap0(:)        !! column pointers (0-based, as in upstream)
    integer,intent(in)           :: Ai0(:)        !! row indices (0-based, as in upstream)
    real(wp),intent(in)          :: Ax(:)         !! values
    real(wp),intent(in)          :: b(:)          !! right-hand side
    real(wp),intent(in),optional :: xsol(:)       !! upstream's solution (absent for the error cases)
    integer,intent(in),optional  :: expected      !! the expected error code (absent if the case must succeed)

    integer(ip),allocatable :: Ap(:), Ai(:), Lp(:), Li(:), Lnz(:), etree(:), iwork(:)
    integer(ip),allocatable :: irow(:), icol(:)
    real(wp),allocatable :: Lx(:), D(:), Dinv(:), fwork(:), x(:)
    logical,allocatable :: bwork(:)
    integer(ip) :: sumLnz, status, j, p
    real(wp) :: tol

    write(*,'(A)') ' '//name
    Ap = int(Ap0, ip) + 1_ip
    Ai = int(Ai0, ip) + 1_ip
    allocate(Lp(n+1), Lnz(n), etree(n), iwork(3*n), D(n), Dinv(n), fwork(n), bwork(n))

    sumLnz = qdldl_etree(n, Ap, Ai, iwork, Lnz, etree)
    if (sumLnz < 0) then
        status = sumLnz
    else
        allocate(Li(max(1_ip, sumLnz)), Lx(max(1_ip, sumLnz)))
        status = qdldl_factor(n, Ap, Ai, Ax, Lp, Li, Lx, D, Dinv, Lnz, etree, bwork, iwork, fwork)
    end if

    if (present(expected)) then
        call check(status == expected, name//': error code')
        return
    end if

    call check(status >= 0, name//': factorisation')
    if (status < 0) return
    x = b
    call qdldl_solve(n, Lp, Li, Lx, Dinv, x)
    tol = max(1.0e-4_wp, 1.0e4_wp * epsilon(1.0_wp)) * max(1.0_wp, maxval(abs(xsol)))
    call check(max_abs_diff(x, xsol) < tol, name//': solution against upstream''s')

    ! the scaled residual, with the matrix in coordinate form
    allocate(irow(size(Ai)), icol(size(Ai)))
    do j = 1, n
        do p = Ap(j), Ap(j+1) - 1
            irow(p) = Ai(p)
            icol(p) = j
        end do
    end do
    call check(residual_norm(irow, icol, Ax, x, b) < res_tol, name//': scaled residual')

    end subroutine run_case
    !*************************************************************************************

    !*************************************************************************************
    !>
    !  The two triangular solves and the diagonal, called separately, give what
    !  [[qdldl_solve]] gives; and \(L\) and \(D\) of the 2x2 case are the exact ones.

    subroutine test_separate_solves()

    integer(ip),parameter :: n = 2
    integer(ip) :: Ap(3), Ai(3), Lp(3), Li(1), Lnz(2), etree(2), iwork(6), sumLnz, npos
    real(wp) :: Ax(3), Lx(1), D(2), Dinv(2), fwork(2), x(2), y(2)
    logical :: bwork(2)

    write(*,'(A)') ' separate solves'
    Ap = [1, 2, 4]; Ai = [1, 1, 2]; Ax = [1.0_wp, 1.0_wp, -1.0_wp]
    sumLnz = qdldl_etree(n, Ap, Ai, iwork, Lnz, etree)
    call check(sumLnz == 1 .and. etree(1) == 2 .and. etree(2) == qdldl_unknown, 'etree of the 2x2 case')
    npos = qdldl_factor(n, Ap, Ai, Ax, Lp, Li, Lx, D, Dinv, Lnz, etree, bwork, iwork, fwork)
    call check(npos == 1, 'one positive pivot')
    call check(Lp(1) == 1 .and. Lp(3) == 2 .and. Li(1) == 2 .and. Lx(1) == 1.0_wp, 'L of the 2x2 case')
    call check(D(1) == 1.0_wp .and. D(2) == -2.0_wp, 'D of the 2x2 case')
    x = [2.0_wp, 4.0_wp]
    y = x
    call qdldl_solve(n, Lp, Li, Lx, Dinv, x)
    call qdldl_lsolve(n, Lp, Li, Lx, y)
    y = y * Dinv
    call qdldl_ltsolve(n, Lp, Li, Lx, y)
    call check(all(x == y) .and. max_abs_diff(x, [3.0_wp, -1.0_wp]) <= 4*epsilon(1.0_wp), &
               'Lsolve, D, and Ltsolve together')

    end subroutine test_separate_solves
    !*************************************************************************************

    !*************************************************************************************
    !>
    !  The column of a zero pivot is reported, and the tolerance of [[qdldl_factor_ext]]
    !  treats a small pivot as zero.

    subroutine test_zero_col()

    integer(ip),parameter :: n = 3
    integer(ip) :: Ap(4), Ai(5), Lp(4), Li(3), Lnz(3), etree(3), iwork(9), sumLnz, zc, npos, nneg, nzero, nreg
    real(wp) :: Ax(5), Lx(3), D(3), Dinv(3), fwork(3)
    logical :: bwork(3)
    integer(ip) :: r
    integer :: status

    write(*,'(A)') ' zero pivot column'
    ! [1 1 0; 1 1 0; 0 0 1] + 1e-10 in (2,2): the second pivot is (almost) zero
    Ap = [1, 2, 4, 5]; Ai = [1, 1, 2, 3, 0]; Ax = [1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 0.0_wp]
    sumLnz = qdldl_etree(n, Ap, Ai, iwork, Lnz, etree)
    r = qdldl_factor(n, Ap, Ai, Ax, Lp, Li, Lx, D, Dinv, Lnz, etree, bwork, iwork, fwork, zero_col=zc)
    call check(r == qdldl_error_zero_pivot .and. zc == 2, 'zero pivot in column 2')

    Ax(3) = 1.0_wp + 1.0e-3_wp
    status = qdldl_factor_ext(n, Ap, Ai, Ax, Lp, Li, Lx, D, Dinv, Lnz, etree, bwork, iwork, fwork, &
                              0.0_wp, .false., 0.0_wp, 0.0_wp, npos, nneg, nzero, nreg, zc)
    call check(status == qdldl_success .and. npos == 3, 'small pivot accepted with no tolerance')
    status = qdldl_factor_ext(n, Ap, Ai, Ax, Lp, Li, Lx, D, Dinv, Lnz, etree, bwork, iwork, fwork, &
                              1.0e-2_wp, .false., 0.0_wp, 0.0_wp, npos, nneg, nzero, nreg, zc)
    call check(status == qdldl_error_zero_pivot .and. zc == 2, 'small pivot is zero within the tolerance')
    status = qdldl_factor_ext(n, Ap, Ai, Ax, Lp, Li, Lx, D, Dinv, Lnz, etree, bwork, iwork, fwork, &
                              1.0e-2_wp, .true., 0.0_wp, 0.5_wp, npos, nneg, nzero, nreg, zc)
    call check(status == qdldl_success .and. nzero == 1 .and. nreg == 1 .and. npos == 2 .and. D(2) == 0.5_wp, &
               'small pivot regularized')

    end subroutine test_zero_col
    !*************************************************************************************

    end program test_core
!*****************************************************************************************
