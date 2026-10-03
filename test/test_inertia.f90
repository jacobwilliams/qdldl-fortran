!*****************************************************************************************
!>
!  Stage 4: inertia, zero pivots, regularization, iterative refinement, and the
!  pivot growth statistics.

    program test_inertia

    use qdldl_kinds, only: wp, ip
    use qdldl_module
    use qdldl_test_utils

    implicit none

    real(wp),parameter :: res_tol = 100.0_wp * epsilon(1.0_wp)   !! scaled residual tolerance

    write(*,'(A)') 'test_inertia'

    call test_quasi_definite_inertia()
    call test_singular_h()
    call test_kkt_static_regularization()
    call test_not_quasi_definite()
    call test_zero_pivot_tolerance()
    call test_growth()

    call finish_tests('test_inertia')

    contains

    !*************************************************************************************
    !>
    !  A quasi-definite KKT matrix with known block sizes: the inertia is exact with
    !  every ordering.

    subroutine test_quasi_definite_inertia()

    integer,parameter :: orders(3) = [qdldl_order_natural, qdldl_order_rcm, qdldl_order_amd]
    type(qdldl_type) :: ldl
    integer(ip),allocatable :: irow(:), icol(:)
    real(wp),allocatable :: val(:)
    integer(ip) :: npos, nneg, nzero, nreg
    integer :: o, istat

    write(*,'(A)') ' quasi-definite inertia'
    call kkt_from_laplacian(12_ip, 50_ip, 1.0e-3_wp, irow, icol, val)
    do o = 1, 3
        call ldl%analyze(194_ip, irow, icol, istat, ordering=orders(o))
        call ldl%factor(val, istat)
        call ldl%inertia(npos, nneg, nzero, nreg)
        call check(istat == qdldl_success .and. npos == 144 .and. nneg == 50 .and. nzero == 0 .and. nreg == 0, &
                   'KKT matrix: 144 positive, 50 negative')
    end do

    end subroutine test_quasi_definite_inertia
    !*************************************************************************************

    !*************************************************************************************
    !>
    !  A KKT matrix \([H\ J^T; J\ 0]\) with a singular \(H\) (but a nonsingular KKT
    !  matrix): the zero pivot is reported, or regularized; static regularization
    !  plus refinement reaches the unregularized solution.

    subroutine test_singular_h()

    ! H = diag(2, 1, 0), J = [1 1 1; 0 0 1]
    integer(ip),parameter :: n = 5
    integer(ip),parameter :: irow(*) = [1, 2, 3, 1, 2, 3, 3]
    integer(ip),parameter :: icol(*) = [1, 2, 3, 4, 4, 4, 5]
    real(wp),parameter :: val(*) = [2.0_wp, 1.0_wp, 0.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp]
    integer,parameter :: signs(*) = [1, 1, 1, -1, -1]
    type(qdldl_type) :: ldl
    real(wp) :: b(n), x(n)
    real(wp),allocatable :: a(:,:), xd(:)
    integer(ip) :: npos, nneg, nzero, nreg
    integer :: istat

    write(*,'(A)') ' singular H'
    b = [1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp]
    call coo_to_dense(n, irow, icol, val, a)
    call dense_solve(a, b, xd)

    ldl%ordering = qdldl_order_natural
    call ldl%analyze(n, irow, icol, istat)
    call ldl%factor(val, istat)
    call check(istat == qdldl_error_zero_pivot .and. ldl%zero_pivot_column == 3, 'zero pivot reported in column 3')

    ! dynamic regularization: completes, one regularized pivot
    ldl%regularize = .true.
    call ldl%set_signs(signs, istat)
    call ldl%factor(val, istat)
    call ldl%inertia(npos, nneg, nzero, nreg)
    call check(istat == qdldl_success .and. nzero == 1 .and. nreg == 1 .and. npos == 2 .and. nneg == 2, &
               'dynamic regularization: completes, one zero pivot regularized')
    x = b
    call ldl%solve(x, istat, refine=10)
    call check(istat == qdldl_success .and. ldl%refine_steps > 0, 'refinement steps taken')
    call check(max_abs_diff(x, xd) <= 1.0e3_wp*epsilon(1.0_wp)*maxval(abs(xd)), &
               'dynamic regularization + refinement: unregularized solution')

    ! static regularization: quasi-definite, so no pivot is zero; refinement removes it
    ldl%regularize = .false.
    ldl%static_reg = sqrt(epsilon(1.0_wp))
    call ldl%factor(val, istat)
    call ldl%inertia(npos, nneg, nzero, nreg)
    call check(istat == qdldl_success .and. npos == 3 .and. nneg == 2 .and. nreg == 0, &
               'static regularization: factors, with the inertia of the KKT matrix')
    x = b
    call ldl%solve(x, istat)
    call check(max_abs_diff(x, xd) > 1.0e3_wp*epsilon(1.0_wp)*maxval(abs(xd)), &
               'static regularization without refinement: perturbed solution')
    ldl%max_refine = 10
    x = b
    call ldl%solve(x, istat)
    call check(istat == qdldl_success .and. max_abs_diff(x, xd) <= 1.0e3_wp*epsilon(1.0_wp)*maxval(abs(xd)), &
               'static regularization + refinement: unregularized solution')
    call check(ldl%residual <= 10*epsilon(1.0_wp), 'refined scaled residual at round-off')

    end subroutine test_singular_h
    !*************************************************************************************

    !*************************************************************************************
    !>
    !  The direct QP case of sqpopt: \([H\ J^T; J\ 0]\) with \(H\) positive definite and
    !  a zero (2,2) block, ordered by AMD. Static regularization makes it
    !  quasi-definite; refinement against the unregularized matrix gives its solution.

    subroutine test_kkt_static_regularization()

    type(qdldl_type) :: ldl
    integer(ip),allocatable :: irow(:), icol(:)
    real(wp),allocatable :: val(:), b(:), x(:)
    integer,allocatable :: signs(:)
    integer(ip) :: n, k, npos, nneg, nzero
    integer :: istat

    write(*,'(A)') ' KKT with a zero (2,2) block'
    call kkt_from_laplacian(15_ip, 60_ip, 0.0_wp, irow, icol, val)
    n = 225 + 60
    allocate(b(n), signs(n))
    do k = 1, n
        b(k) = real(mod(3*k, 11_ip), wp) - 5.0_wp
    end do
    signs(1:225) = 1
    signs(226:n) = -1

    call ldl%analyze(n, irow, icol, istat)
    call ldl%set_signs(signs, istat)
    ldl%static_reg = sqrt(epsilon(1.0_wp))   ! (1.5e-8 in double precision)
    ldl%max_refine = 10
    call ldl%factor(val, istat)
    call ldl%inertia(npos, nneg, nzero)
    call check(istat == qdldl_success .and. npos == 225 .and. nneg == 60, 'factors, with the right inertia')
    x = b
    call ldl%solve(x, istat)
    call check(istat == qdldl_success .and. residual_norm(irow, icol, val, x, b) < res_tol, &
               'refined solution of the unregularized system')

    end subroutine test_kkt_static_regularization
    !*************************************************************************************

    !*************************************************************************************
    !>
    !  An indefinite matrix that is not quasi-definite (an indefinite \(H\) block):
    !  dynamic regularization with the expected signs completes, and counts the
    !  replacement.

    subroutine test_not_quasi_definite()

    ! H = [1 2; 2 1] (eigenvalues 3, -1), C = [1]: [H B^T; B -C] with B = [1 0]
    integer(ip),parameter :: n = 3
    integer(ip),parameter :: irow(*) = [1, 1, 2, 1, 3]
    integer(ip),parameter :: icol(*) = [1, 2, 2, 3, 3]
    real(wp),parameter :: val(*) = [1.0_wp, 2.0_wp, 1.0_wp, 1.0_wp, -1.0_wp]
    type(qdldl_type) :: ldl
    integer(ip) :: npos, nneg, nzero, nreg
    integer :: istat

    write(*,'(A)') ' not quasi-definite'
    ldl%ordering = qdldl_order_natural
    call ldl%analyze(n, irow, icol, istat)
    call ldl%factor(val, istat)
    call ldl%inertia(npos, nneg, nzero, nreg)
    ! without regularization it factors (no zero pivot); the inertia is A's
    call check(istat == qdldl_success .and. npos == 1 .and. nneg == 2 .and. nreg == 0, &
               'factors without regularization (inertia 1, 2)')
    ldl%regularize = .true.
    call ldl%set_signs([1, 1, -1], istat)
    call ldl%factor(val, istat)
    call ldl%inertia(npos, nneg, nzero, nreg)
    call check(istat == qdldl_success .and. nreg == 1, 'with the expected signs: the wrong-signed pivot is regularized')
    call check(ldl%n_negative == 2, 'the counts are of the pivots as computed')

    end subroutine test_not_quasi_definite
    !*************************************************************************************

    !*************************************************************************************
    !>
    !  The zero-pivot tolerance is relative to the largest entry of the matrix.

    subroutine test_zero_pivot_tolerance()

    ! [1 1; 1 1+1e-3] (scaled by 100): second pivot 0.1, largest entry 100.1
    integer(ip),parameter :: irow(*) = [1, 1, 2]
    integer(ip),parameter :: icol(*) = [1, 2, 2]
    real(wp),parameter :: val(*) = [100.0_wp, 100.0_wp, 100.1_wp]
    type(qdldl_type) :: ldl
    integer :: istat

    write(*,'(A)') ' zero-pivot tolerance'
    call ldl%analyze(2_ip, irow, icol, istat, ordering=qdldl_order_natural)
    ldl%zero_pivot_tol = 1.0e-4_wp
    call ldl%factor(val, istat)
    call check(istat == qdldl_success, 'pivot 1e-3 relative: not zero with tolerance 1e-4')
    ldl%zero_pivot_tol = 1.0e-2_wp
    call ldl%factor(val, istat)
    call check(istat == qdldl_error_zero_pivot .and. ldl%zero_pivot_column == 2, &
               'pivot 1e-3 relative: zero with tolerance 1e-2')
    ldl%regularize = .true.
    call ldl%factor(val, istat)
    call check(istat == qdldl_success .and. ldl%n_zero == 1 .and. ldl%n_regularized == 1, &
               'regularized within the tolerance')

    end subroutine test_zero_pivot_tolerance
    !*************************************************************************************

    !*************************************************************************************
    !>
    !  A tiny pivot that isn't zero shows up as large growth.

    subroutine test_growth()

    ! [d 1; 1 0] with a tiny d: L(2,1) = 1/d
    integer(ip),parameter :: irow(*) = [1, 1, 2]
    integer(ip),parameter :: icol(*) = [1, 2, 2]
    type(qdldl_type) :: ldl
    integer :: istat

    write(*,'(A)') ' growth'
    call ldl%analyze(2_ip, irow, icol, istat, ordering=qdldl_order_natural)
    call ldl%factor([1.0e-3_wp, 1.0_wp, 0.0_wp], istat)
    call check(istat == qdldl_success .and. abs(ldl%max_abs_l - 1.0e3_wp) <= 1.0e3_wp*1.0e3_wp*epsilon(1.0_wp), &
               'max |L| = 1/d')
    call check(abs(ldl%pivot_ratio - 1.0e6_wp) <= 1.0e6_wp*1.0e3_wp*epsilon(1.0_wp), 'pivot ratio = 1/d^2')
    call ldl%factor([1.0_wp, 0.5_wp, 2.0_wp], istat)
    call check(istat == qdldl_success .and. ldl%max_abs_l == 0.5_wp, 'max |L| of a well-conditioned matrix')

    end subroutine test_growth
    !*************************************************************************************

    end program test_inertia
!*****************************************************************************************
