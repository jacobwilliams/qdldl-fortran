!*****************************************************************************************
!>
!  Helpers shared by the tests: pass/fail bookkeeping, products with a matrix in
!  coordinate form, norms, and generators of test matrices.

    module qdldl_test_utils

    use qdldl_kinds, only: wp, ip

    implicit none

    private

    integer,public :: n_failed = 0   !! number of failed checks so far
    integer,public :: n_checks = 0   !! number of checks so far

    public :: check, finish_tests
    public :: coo_multiply, residual_norm, max_abs_diff
    public :: lcg_init, lcg_uniform
    public :: laplacian_2d, kkt_from_laplacian
    public :: coo_to_dense, dense_solve, random_qd

    integer(8) :: lcg_state = 12345_8   !! state of the random number generator

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  Records one check, and prints it if it fails (or always, with `verbose`).

    subroutine check(ok, name)

    logical,intent(in)          :: ok     !! the condition that must hold
    character(len=*),intent(in) :: name   !! what is checked

    n_checks = n_checks + 1
    if (ok) then
        write(*,'(A)') '  pass: '//name
    else
        n_failed = n_failed + 1
        write(*,'(A)') '  FAIL: '//name
    end if

    end subroutine check
!*****************************************************************************************

!*****************************************************************************************
!>
!  Prints a summary, and stops with a nonzero code if anything failed.

    subroutine finish_tests(name)

    character(len=*),intent(in) :: name   !! name of the test program

    write(*,'(A,I0,A,I0,A)') trim(name)//': ', n_checks - n_failed, ' of ', n_checks, ' checks passed'
    if (n_failed > 0) error stop 1

    end subroutine finish_tests
!*****************************************************************************************

!*****************************************************************************************
!>
!  \(y = Ax\) for a symmetric matrix given by entries of one or both triangles in
!  coordinate form (an off-diagonal entry \((i,j)\) stands for both \((i,j)\) and \((j,i)\);
!  duplicates are added).

    subroutine coo_multiply(irow, icol, val, x, y)

    integer(ip),intent(in) :: irow(:)   !! row indices
    integer(ip),intent(in) :: icol(:)   !! column indices
    real(wp),intent(in)    :: val(:)    !! values
    real(wp),intent(in)    :: x(:)      !! the vector
    real(wp),intent(out)   :: y(:)      !! the product

    integer(ip) :: k

    y = 0.0_wp
    do k = 1, size(irow, kind=ip)
        y(irow(k)) = y(irow(k)) + val(k) * x(icol(k))
        if (irow(k) /= icol(k)) y(icol(k)) = y(icol(k)) + val(k) * x(irow(k))
    end do

    end subroutine coo_multiply
!*****************************************************************************************

!*****************************************************************************************
!>
!  The scaled residual \(\lVert Ax-b \rVert_\infty / (\lVert A \rVert_\infty \lVert x \rVert_\infty
!  + \lVert b \rVert_\infty)\), for a matrix in coordinate form.

    function residual_norm(irow, icol, val, x, b) result(r)

    integer(ip),intent(in) :: irow(:)   !! row indices
    integer(ip),intent(in) :: icol(:)   !! column indices
    real(wp),intent(in)    :: val(:)    !! values
    real(wp),intent(in)    :: x(:)      !! the solution
    real(wp),intent(in)    :: b(:)      !! the right-hand side
    real(wp)               :: r         !! the scaled residual

    real(wp),allocatable :: y(:), rowsum(:)
    integer(ip) :: k
    real(wp) :: anorm

    allocate(y(size(x)), rowsum(size(x)))
    call coo_multiply(irow, icol, val, x, y)
    rowsum = 0.0_wp
    do k = 1, size(irow, kind=ip)
        rowsum(irow(k)) = rowsum(irow(k)) + abs(val(k))
        if (irow(k) /= icol(k)) rowsum(icol(k)) = rowsum(icol(k)) + abs(val(k))
    end do
    anorm = maxval(rowsum)
    r = maxval(abs(y - b)) / (anorm * maxval(abs(x)) + maxval(abs(b)))

    end function residual_norm
!*****************************************************************************************

!*****************************************************************************************
!>
!  \(\max_i |x_i - y_i|\).

    pure function max_abs_diff(x, y) result(d)

    real(wp),intent(in) :: x(:)   !! a vector
    real(wp),intent(in) :: y(:)   !! another vector of the same size
    real(wp)            :: d      !! the largest difference

    d = maxval(abs(x - y))

    end function max_abs_diff
!*****************************************************************************************

!*****************************************************************************************
!>
!  Seeds the (portable, reproducible) random number generator.

    subroutine lcg_init(seed)

    integer,intent(in) :: seed   !! the seed

    lcg_state = int(max(1, seed), 8)

    end subroutine lcg_init
!*****************************************************************************************

!*****************************************************************************************
!>
!  A uniform random number in \((-1, 1)\), from the Park-Miller generator
!  (the same sequence in every real kind).

    function lcg_uniform() result(r)

    real(wp) :: r   !! the random number

    lcg_state = mod(16807_8 * lcg_state, 2147483647_8)   ! Park-Miller (no overflow in 64 bits)
    r = real(lcg_state, wp) / 2147483647.0_wp
    r = 2.0_wp * r - 1.0_wp

    end function lcg_uniform
!*****************************************************************************************

!*****************************************************************************************
!>
!  The 5-point Laplacian of an `m` by `m` grid (upper triangle, coordinate form):
!  positive definite, of order `m**2`.

    subroutine laplacian_2d(m, irow, icol, val)

    integer(ip),intent(in)                :: m        !! grid size
    integer(ip),allocatable,intent(out)   :: irow(:)  !! row indices
    integer(ip),allocatable,intent(out)   :: icol(:)  !! column indices
    real(wp),allocatable,intent(out)      :: val(:)   !! values

    integer(ip) :: i, j, k, nz

    allocate(irow(5*m*m), icol(5*m*m), val(5*m*m))
    nz = 0
    do j = 1, m
        do i = 1, m
            k = i + (j-1)*m
            nz = nz + 1; irow(nz) = k; icol(nz) = k; val(nz) = 4.0_wp
            if (i < m) then
                nz = nz + 1; irow(nz) = k; icol(nz) = k + 1; val(nz) = -1.0_wp
            end if
            if (j < m) then
                nz = nz + 1; irow(nz) = k; icol(nz) = k + m; val(nz) = -1.0_wp
            end if
        end do
    end do
    irow = irow(1:nz); icol = icol(1:nz); val = val(1:nz)

    end subroutine laplacian_2d
!*****************************************************************************************

!*****************************************************************************************
!>
!  A quasi-definite KKT matrix \(\begin{bmatrix} H & B^T \\ B & -\epsilon I \end{bmatrix}\)
!  with \(H\) the 2-D Laplacian of an `m` by `m` grid and \(B\) the first `nc` rows of the
!  same Laplacian (so the matrix has `m**2` positive and `nc` negative eigenvalues).
!  Upper triangle, coordinate form. Rows `1..m**2` are \(H\); rows `m**2+1..m**2+nc` the constraints.

    subroutine kkt_from_laplacian(m, nc, eps, irow, icol, val)

    integer(ip),intent(in)                :: m        !! grid size
    integer(ip),intent(in)                :: nc       !! number of constraints (`<= m**2`)
    real(wp),intent(in)                   :: eps      !! the (positive) regularization of the (2,2) block
    integer(ip),allocatable,intent(out)   :: irow(:)  !! row indices
    integer(ip),allocatable,intent(out)   :: icol(:)  !! column indices
    real(wp),allocatable,intent(out)      :: val(:)   !! values

    integer(ip),allocatable :: hr(:), hc(:)
    real(wp),allocatable :: hv(:)
    integer(ip) :: k, nh, nz, nhz

    call laplacian_2d(m, hr, hc, hv)
    nh = m*m
    nhz = size(hr, kind=ip)
    allocate(irow(3*nhz + nc), icol(3*nhz + nc), val(3*nhz + nc))
    nz = 0
    ! H
    do k = 1, nhz
        nz = nz + 1; irow(nz) = hr(k); icol(nz) = hc(k); val(nz) = hv(k)
    end do
    ! B (rows of the Laplacian with index <= nc), stored as B^T in the upper triangle
    do k = 1, nhz
        if (hr(k) <= nc) then
            nz = nz + 1; irow(nz) = hc(k); icol(nz) = nh + hr(k); val(nz) = 0.5_wp * hv(k)
        end if
        if (hc(k) <= nc .and. hr(k) /= hc(k)) then
            nz = nz + 1; irow(nz) = hr(k); icol(nz) = nh + hc(k); val(nz) = 0.5_wp * hv(k)
        end if
    end do
    ! -eps I
    do k = 1, nc
        nz = nz + 1; irow(nz) = nh + k; icol(nz) = nh + k; val(nz) = -eps
    end do
    irow = irow(1:nz); icol = icol(1:nz); val = val(1:nz)

    end subroutine kkt_from_laplacian
!*****************************************************************************************

!*****************************************************************************************
!>
!  The dense symmetric matrix of entries in coordinate form (as [[coo_multiply]]).

    subroutine coo_to_dense(n, irow, icol, val, a)

    integer(ip),intent(in)            :: n         !! order
    integer(ip),intent(in)            :: irow(:)   !! row indices
    integer(ip),intent(in)            :: icol(:)   !! column indices
    real(wp),intent(in)               :: val(:)    !! values
    real(wp),allocatable,intent(out)  :: a(:,:)    !! the dense matrix

    integer(ip) :: k

    allocate(a(n,n), source=0.0_wp)
    do k = 1, size(irow, kind=ip)
        a(irow(k), icol(k)) = a(irow(k), icol(k)) + val(k)
        if (irow(k) /= icol(k)) a(icol(k), irow(k)) = a(icol(k), irow(k)) + val(k)
    end do

    end subroutine coo_to_dense
!*****************************************************************************************

!*****************************************************************************************
!>
!  Solves a small dense system by Gaussian elimination with partial pivoting.

    subroutine dense_solve(a, b, x)

    real(wp),intent(in)               :: a(:,:)   !! the matrix
    real(wp),intent(in)               :: b(:)     !! the right-hand side
    real(wp),allocatable,intent(out)  :: x(:)     !! the solution

    real(wp),allocatable :: m(:,:), row(:)
    integer :: n, i, k, piv
    real(wp) :: t

    n = size(b)
    m = a
    x = b
    do k = 1, n
        piv = k - 1 + maxloc(abs(m(k:n, k)), dim=1)
        if (piv /= k) then
            row = m(k,:); m(k,:) = m(piv,:); m(piv,:) = row
            t = x(k); x(k) = x(piv); x(piv) = t
        end if
        do i = k + 1, n
            t = m(i,k) / m(k,k)
            m(i,k:n) = m(i,k:n) - t * m(k,k:n)
            x(i) = x(i) - t * x(k)
        end do
    end do
    do k = n, 1, -1
        x(k) = (x(k) - dot_product(m(k,k+1:n), x(k+1:n))) / m(k,k)
    end do

    end subroutine dense_solve
!*****************************************************************************************

!*****************************************************************************************
!>
!  A random sparse quasi-definite matrix (upper triangle, coordinate form): each
!  index `i` gets a random sign `s(i)`; entries between indices of the same sign
!  form a diagonally dominant block of that sign, and entries between indices of
!  different signs are arbitrary. So the matrix has exactly `count(s > 0)` positive
!  and `count(s < 0)` negative eigenvalues.

    subroutine random_qd(n, density, s, irow, icol, val)

    integer(ip),intent(in)                :: n         !! order
    real(wp),intent(in)                   :: density   !! probability of an off-diagonal entry
    integer,allocatable,intent(out)       :: s(:)      !! the sign of each index
    integer(ip),allocatable,intent(out)   :: irow(:)   !! row indices
    integer(ip),allocatable,intent(out)   :: icol(:)   !! column indices
    real(wp),allocatable,intent(out)      :: val(:)    !! values

    real(wp),allocatable :: rowsum(:)
    integer(ip) :: i, j, nz, cap
    real(wp) :: v

    allocate(s(n), rowsum(n))
    do i = 1, n
        s(i) = merge(1, -1, lcg_uniform() > -0.2_wp)
    end do
    cap = n + int(density * real(n, wp) * real(n, wp)) + 10*n
    allocate(irow(cap), icol(cap), val(cap))
    rowsum = 0.0_wp
    nz = 0
    do j = 1, n
        do i = 1, j - 1
            if (0.5_wp * (lcg_uniform() + 1.0_wp) < density .and. nz < cap - n) then
                v = lcg_uniform()
                nz = nz + 1; irow(nz) = i; icol(nz) = j; val(nz) = v
                if (s(i) == s(j)) then
                    rowsum(i) = rowsum(i) + abs(v)
                    rowsum(j) = rowsum(j) + abs(v)
                end if
            end if
        end do
    end do
    do i = 1, n
        nz = nz + 1; irow(nz) = i; icol(nz) = i
        val(nz) = s(i) * (rowsum(i) + 0.1_wp + 0.5_wp * (lcg_uniform() + 1.0_wp))
    end do
    irow = irow(1:nz); icol = icol(1:nz); val = val(1:nz)

    end subroutine random_qd
!*****************************************************************************************

!*****************************************************************************************
    end module qdldl_test_utils
!*****************************************************************************************
