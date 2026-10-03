!*****************************************************************************************
!>
!  Factor and solve a small KKT system with [[qdldl_type]]:
!
!  \[ \begin{bmatrix} H & J^T \\ J & 0 \end{bmatrix}
!     \begin{bmatrix} x \\ y \end{bmatrix} =
!     \begin{bmatrix} g \\ c \end{bmatrix} \]
!
!  with \(H\) positive definite (3x3) and \(J\) of full rank (2x3). The zero (2,2)
!  block makes the matrix not quite quasi-definite, so a small static
!  regularization (with the expected signs of the pivots) is added, and iterative
!  refinement removes its effect from the solution.

    program example_kkt

    use qdldl_module, only: qdldl_type, qdldl_wp, qdldl_ip, qdldl_success, qdldl_status_message

    implicit none

    integer,parameter :: wp = qdldl_wp
    integer,parameter :: ip = qdldl_ip
    integer(ip),parameter :: n = 5

    ! the upper triangle, in coordinate form (any order; duplicates would be added)
    integer(ip),parameter :: irow(*) = [1, 1, 2, 2, 3,   1, 2, 3, 3]
    integer(ip),parameter :: icol(*) = [1, 2, 2, 3, 3,   4, 4, 4, 5]
    real(wp),parameter    :: val(*)  = [4.0_wp, 1.0_wp, 5.0_wp, -1.0_wp, 6.0_wp, &   ! H
                                        1.0_wp, 2.0_wp, -1.0_wp, 1.0_wp]            ! J^T

    type(qdldl_type) :: ldl
    real(wp) :: rhs(n), x(n), ax(n)
    integer(ip) :: npos, nneg, nzero
    integer :: istat

    ! the pattern, once (the ordering is AMD by default)
    call ldl%analyze(n, irow, icol, istat)
    if (istat /= qdldl_success) error stop qdldl_status_message(istat)
    write(*,'(A,I0)') 'nonzeros of L: ', ldl%nnz_l

    ! regularization: +eps on the H rows, -eps on the constraint rows, and refinement
    call ldl%set_signs([1, 1, 1, -1, -1], istat)
    ldl%static_reg = 1.0e-8_wp
    ldl%max_refine = 3

    ! the values (as often as needed, with the same pattern)
    call ldl%factor(val, istat)
    if (istat /= qdldl_success) error stop qdldl_status_message(istat)
    call ldl%inertia(npos, nneg, nzero)
    write(*,'(A,I0,A,I0,A,I0,A)') 'inertia: (', npos, ', ', nneg, ', ', nzero, ')'

    ! solve (in place)
    rhs = [1.0_wp, 2.0_wp, 3.0_wp, 0.5_wp, -1.0_wp]
    x = rhs
    call ldl%solve(x, istat)
    if (istat /= qdldl_success) error stop qdldl_status_message(istat)
    write(*,'(A,*(F12.8))') 'x =', x
    write(*,'(A,I0)') 'refinement steps: ', ldl%refine_steps

    ! check: the residual with the unregularized matrix
    call ldl%multiply(x, ax)
    write(*,'(A,ES10.3)') 'max |Ax - b| = ', maxval(abs(ax - rhs))

    end program example_kkt
!*****************************************************************************************
