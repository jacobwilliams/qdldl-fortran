!*****************************************************************************************
!>
!  Timings on scalable matrices: 2-D and 3-D Laplacians, and the quasi-definite
!  KKT matrices built from them, \(\begin{bmatrix} I & A^T \\ A & -\epsilon I
!  \end{bmatrix}\) (with \(A\) the Laplacian). For each ordering: the nonzeros of
!  \(L\), and the times of the analysis, the factorization, and a solve.
!
!  Run with `fpm run --example benchmark --profile release`, optionally with a
!  scale factor: `-- 2` doubles the grid sizes.

    program benchmark

    use qdldl_module, only: qdldl_type, qdldl_wp, qdldl_ip, qdldl_success, qdldl_status_message, &
                            qdldl_order_natural, qdldl_order_rcm, qdldl_order_amd
    use, intrinsic :: iso_fortran_env, only: int64

    implicit none

    integer,parameter :: wp = qdldl_wp
    integer,parameter :: ip = qdldl_ip
    integer,parameter :: orders(3) = [qdldl_order_natural, qdldl_order_rcm, qdldl_order_amd]
    character(len=*),parameter :: oname(3) = ['natural', 'rcm    ', 'amd    ']

    integer(ip),allocatable :: irow(:), icol(:)
    real(wp),allocatable :: val(:)
    integer(ip) :: n, scale
    integer :: c, o
    character(len=32) :: arg, name

    scale = 1
    if (command_argument_count() > 0) then
        call get_command_argument(1, arg)
        read(arg, *) scale
    end if

    write(*,'(A)') 'matrix                    n       nnz(A)  ordering        nnz(L)   analyze    factor     solve  (s)'
    do c = 1, 4
        select case (c)
        case (1)
            call laplacian(2, 100*scale, n, irow, icol, val, .false.)
            write(name,'(A,I0,A,I0)') '2-D Laplacian ', 100*scale, '^2'
        case (2)
            call laplacian(2, 100*scale, n, irow, icol, val, .true.)
            write(name,'(A,I0,A,I0)') '2-D KKT ', 100*scale, '^2'
        case (3)
            call laplacian(3, 20*scale, n, irow, icol, val, .false.)
            write(name,'(A,I0,A,I0)') '3-D Laplacian ', 20*scale, '^3'
        case (4)
            call laplacian(3, 20*scale, n, irow, icol, val, .true.)
            write(name,'(A,I0,A,I0)') '3-D KKT ', 20*scale, '^3'
        case default
            error stop 'unknown case'
        end select
        do o = 1, 3
            call run(name, orders(o), oname(o))
        end do
    end do

    contains

    subroutine run(name, order, oname)
    !! Analyze, factor, and solve one matrix with one ordering, and print the times.
    character(len=*),intent(in) :: name    !! the matrix
    integer,intent(in)          :: order   !! the ordering
    character(len=*),intent(in) :: oname   !! its name
    type(qdldl_type) :: ldl
    real(wp),allocatable :: b(:)
    real(wp) :: t0, t1, t2, t3
    integer :: istat
    t0 = now()
    call ldl%analyze(n, irow, icol, istat, ordering=order)
    if (istat /= qdldl_success) then
        write(*,'(A24,I9,I13,2X,A8,A)') name, n, size(irow), oname, '  '//qdldl_status_message(istat)
        return
    end if
    if (ldl%nnz_l > 200000000_ip) then   ! (natural order on the large 3-D matrices)
        write(*,'(A24,I9,I13,2X,A8,I14,A)') name, n, size(irow), oname, ldl%nnz_l, '  (skipped)'
        return
    end if
    t1 = now()
    call ldl%factor(val, istat)
    t2 = now()
    allocate(b(n))
    b = 1.0_wp
    call ldl%solve(b, istat)
    t3 = now()
    write(*,'(A24,I9,I13,2X,A8,I14,3F10.4)') name, n, size(irow), oname, ldl%nnz_l, t1 - t0, t2 - t1, t3 - t2
    end subroutine run

    function now() result(t)
    !! Wall-clock time in seconds.
    real(wp) :: t   !! the time
    integer(int64) :: ticks, rate
    call system_clock(ticks, rate)
    t = real(ticks, wp) / real(rate, wp)
    end function now

    subroutine laplacian(dim, m, n, irow, icol, val, kkt)
    !! The 2-D (5-point) or 3-D (7-point) Laplacian of a grid with `m` points per side
    !! (upper triangle, coordinate form), or the KKT matrix [I A^T; A -eps I] built from it.
    integer,intent(in)                    :: dim       !! 2 or 3
    integer(ip),intent(in)                :: m         !! grid points per side
    integer(ip),intent(out)               :: n         !! order of the matrix
    integer(ip),allocatable,intent(out)   :: irow(:)   !! row indices
    integer(ip),allocatable,intent(out)   :: icol(:)   !! column indices
    real(wp),allocatable,intent(out)      :: val(:)    !! values
    logical,intent(in)                    :: kkt       !! the KKT matrix instead of the Laplacian
    integer(ip) :: i, j, k, p, nl, nz, nlap, stride(3)
    nl = m**dim
    n = merge(2*nl, nl, kkt)
    allocate(irow(2*(dim+1)*nl + 2*nl), icol(2*(dim+1)*nl + 2*nl), val(2*(dim+1)*nl + 2*nl))
    stride = [1_ip, m, m*m]
    ! the Laplacian (upper triangle)
    nz = 0
    do p = 1, nl
        nz = nz + 1; irow(nz) = p; icol(nz) = p; val(nz) = real(2*dim, wp)
        k = p - 1
        do i = 1, dim
            j = mod(k / stride(i), m)     ! coordinate i of point p (0-based)
            if (j < m - 1) then
                nz = nz + 1; irow(nz) = p; icol(nz) = p + stride(i); val(nz) = -1.0_wp
            end if
        end do
    end do
    if (kkt) then
        ! the (1,2) block is A^T = A, all of it (both triangles), in rows 1..nl and
        ! columns nl+1..2nl; then I in the (1,1) block and -eps I in the (2,2) block
        nlap = nz
        do p = 1, nlap
            if (irow(p) /= icol(p)) then
                nz = nz + 1; irow(nz) = icol(p); icol(nz) = irow(p) + nl; val(nz) = val(p)
            end if
            icol(p) = icol(p) + nl
        end do
        do p = 1, nl
            nz = nz + 1; irow(nz) = p; icol(nz) = p; val(nz) = 1.0_wp
            nz = nz + 1; irow(nz) = p + nl; icol(nz) = p + nl; val(nz) = -1.0e-4_wp
        end do
    end if
    irow = irow(1:nz); icol = icol(1:nz); val = val(1:nz)
    end subroutine laplacian

    end program benchmark
!*****************************************************************************************
