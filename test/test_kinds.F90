!*****************************************************************************************
!>
!  Reports the real and integer kinds of the build, and checks a factorization and
!  solve in them. CI runs the whole test suite in each kind (`REAL32`, `REAL64`,
!  `REAL128`, and `INT64`); every test scales its tolerances with `epsilon(1.0_wp)`.

    program test_kinds

    use qdldl_kinds, only: wp, ip
    use qdldl_module
    use qdldl_test_utils

    implicit none

    type(qdldl_type) :: ldl
    integer(ip),allocatable :: irow(:), icol(:)
    real(wp),allocatable :: val(:), b(:), x(:)
    integer(ip) :: n, k
    integer :: istat

    write(*,'(A)') 'test_kinds'
    write(*,'(A,I0,A,I0,A,ES10.3)') '  real kind: ', qdldl_wp, ' (', storage_size(1.0_qdldl_wp), &
                                    ' bits), epsilon ', epsilon(1.0_qdldl_wp)
    write(*,'(A,I0,A,I0,A)') '  integer kind: ', qdldl_ip, ' (', storage_size(1_qdldl_ip), ' bits)'
    call check(qdldl_wp == wp .and. qdldl_ip == ip, 'the public kinds are the library''s')
#ifdef REAL32
    call check(storage_size(1.0_wp) == 32, 'REAL32 selects 32-bit reals')
#endif
#ifdef REAL128
    call check(storage_size(1.0_wp) == 128, 'REAL128 selects 128-bit reals')
#endif
#ifdef INT64
    call check(storage_size(1_ip) == 64, 'INT64 selects 64-bit integers')
#else
    call check(storage_size(1_ip) == 32, 'default integers are 32-bit')
#endif

    call kkt_from_laplacian(20_ip, 100_ip, 1.0e-2_wp, irow, icol, val)
    n = 500
    allocate(b(n))
    do k = 1, n
        b(k) = 1.0_wp / real(k, wp)
    end do
    call ldl%analyze(n, irow, icol, istat)
    call ldl%factor(val, istat)
    x = b
    call ldl%solve(x, istat)
    call check(istat == qdldl_success .and. residual_norm(irow, icol, val, x, b) < 100*epsilon(1.0_wp), &
               'KKT solve at round-off of this kind')
    call check(ldl%n_positive == 400 .and. ldl%n_negative == 100, 'KKT inertia')

    call finish_tests('test_kinds')

    end program test_kinds
!*****************************************************************************************
