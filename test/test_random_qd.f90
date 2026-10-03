!*****************************************************************************************
!>
!  Random sparse quasi-definite matrices of several sizes and densities, with
!  every ordering: the scaled residual is at round-off, and the inertia is exact.

    program test_random_qd

    use qdldl_kinds, only: wp, ip
    use qdldl_module
    use qdldl_test_utils

    implicit none

    real(wp),parameter :: res_tol = 100.0_wp * epsilon(1.0_wp)   !! scaled residual tolerance
    integer,parameter :: orders(3) = [qdldl_order_natural, qdldl_order_rcm, qdldl_order_amd]
    character(len=*),parameter :: oname(3) = ['natural', 'rcm    ', 'amd    ']
    integer(ip),parameter :: sizes(5) = [2_ip, 10_ip, 50_ip, 200_ip, 500_ip]
    real(wp),parameter :: densities(3) = [0.01_wp, 0.05_wp, 0.2_wp]

    type(qdldl_type) :: ldl
    integer(ip),allocatable :: irow(:), icol(:)
    real(wp),allocatable :: val(:), b(:), x(:)
    integer,allocatable :: s(:)
    integer(ip) :: n, k, npos, nneg, nzero
    integer :: i, j, o, istat, seed
    logical :: ok_res, ok_inertia, ok_status
    character(len=80) :: name

    write(*,'(A)') 'test_random_qd'
    seed = 0
    do i = 1, size(sizes)
        do j = 1, size(densities)
            seed = seed + 1
            call lcg_init(1000 + seed)
            n = sizes(i)
            call random_qd(n, densities(j), s, irow, icol, val)
            allocate(b(n))
            do k = 1, n
                b(k) = lcg_uniform()
            end do
            do o = 1, 3
                write(name,'(A,I0,A,F4.2,A)') 'n = ', n, ', density ', densities(j), ', '//trim(oname(o))
                call ldl%analyze(n, irow, icol, istat, ordering=orders(o))
                ok_status = istat == qdldl_success
                call ldl%factor(val, istat)
                ok_status = ok_status .and. istat == qdldl_success
                x = b
                call ldl%solve(x, istat)
                ok_status = ok_status .and. istat == qdldl_success
                call check(ok_status, trim(name)//': analyze, factor, solve')
                ok_res = residual_norm(irow, icol, val, x, b) < res_tol
                call check(ok_res, trim(name)//': scaled residual')
                call ldl%inertia(npos, nneg, nzero)
                ok_inertia = npos == count(s > 0) .and. nneg == count(s < 0) .and. nzero == 0
                call check(ok_inertia, trim(name)//': inertia')
            end do
            deallocate(b)
        end do
    end do

    call finish_tests('test_random_qd')

    end program test_random_qd
!*****************************************************************************************
