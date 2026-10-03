    program benchmark

    use qdldl_kinds, only: wp, ip
    use qdldl_core, only: qdldl_etree, qdldl_factor, qdldl_solve
    use, intrinsic :: iso_c_binding, only: c_int, c_int8_t, c_double
    use, intrinsic :: iso_fortran_env, only: int64

    implicit none

    integer,parameter :: nsamples = 5
    integer,parameter :: grid_sizes(3) = [20, 50, 100]
    integer,parameter :: etree_repeats = 20
    integer,parameter :: factor_repeats = 5
    integer,parameter :: solve_repeats = 100

    interface
        function c_etree(n, Ap, Ai, work, Lnz, etree) bind(c, name='benchmark_c_etree') result(sumLnz)
            import :: c_int
            integer(c_int),value :: n
            integer(c_int),intent(in) :: Ap(*), Ai(*)
            integer(c_int),intent(out) :: work(*), Lnz(*), etree(*)
            integer(c_int) :: sumLnz
        end function c_etree

        function c_factor(n, Ap, Ai, Ax, Lp, Li, Lx, D, Dinv, Lnz, etree, &
                          bwork, iwork, fwork) bind(c, name='benchmark_c_factor') result(npos)
            import :: c_int, c_int8_t, c_double
            integer(c_int),value :: n
            integer(c_int),intent(in) :: Ap(*), Ai(*)
            real(c_double),intent(in) :: Ax(*)
            integer(c_int),intent(out) :: Lp(*), Li(*)
            real(c_double),intent(out) :: Lx(*), D(*), Dinv(*), fwork(*)
            integer(c_int),intent(in) :: Lnz(*), etree(*)
            integer(c_int8_t),intent(out) :: bwork(*)
            integer(c_int),intent(out) :: iwork(*)
            integer(c_int) :: npos
        end function c_factor

        subroutine c_solve(n, Lp, Li, Lx, Dinv, x) bind(c, name='benchmark_c_solve')
            import :: c_int, c_double
            integer(c_int),value :: n
            integer(c_int),intent(in) :: Lp(*), Li(*)
            real(c_double),intent(in) :: Lx(*), Dinv(*)
            real(c_double),intent(inout) :: x(*)
        end subroutine c_solve
    end interface

    integer(ip),allocatable :: Ap(:), Ai(:), Ap_c(:), Ai_c(:), Lnz_f(:), etree_f(:), work_f(:)
    integer(ip),allocatable :: Lnz_c(:), etree_c(:), work_c(:)
    integer(ip),allocatable :: Lp_f(:), Li_f(:), Lp_c(:), Li_c(:)
    integer(ip),allocatable :: iwork_f(:), iwork_c(:)
    integer(c_int8_t),allocatable :: bwork_c(:)
    logical,allocatable :: bwork_f(:)
    real(wp),allocatable :: Ax(:), Lx_f(:), D_f(:), Dinv_f(:), fwork_f(:)
    real(wp),allocatable :: Lx_c(:), D_c(:), Dinv_c(:), fwork_c(:)
    real(wp),allocatable :: rhs(:,:), xf(:,:), xc(:,:)
    real(wp) :: etree_f_time(nsamples), etree_c_time(nsamples)
    real(wp) :: factor_f_time(nsamples), factor_c_time(nsamples)
    real(wp) :: solve_f_time(nsamples), solve_c_time(nsamples)
    integer(ip) :: n, nnz, nnz_l, result_f, result_c
    integer :: m, sample, i
    character(len=32) :: matrix_name

    if (storage_size(0_ip) /= storage_size(0_c_int) .or. storage_size(0.0_wp) /= storage_size(0.0_c_double)) then
        error stop 'benchmark requires default QDLDL kinds (32-bit integers and 64-bit reals)'
    end if

    write(*,'(A)') '2-D five-point Laplacian, natural ordering; median time per call'
    write(*,'(A)') 'matrix           n      nnz(A)   nnz(L)  stage       Fortran (us)       C (us)   C/Fortran'
    do i = 1, size(grid_sizes)
        m = grid_sizes(i)
        write(matrix_name,'(I0,"x",I0)') m, m
        call make_laplacian(m, n, Ap, Ai, Ax)
        nnz = size(Ai, kind=ip)
        allocate(Ap_c(n+1), Ai_c(nnz))
        Ap_c = Ap - 1
        Ai_c = Ai - 1
        allocate(Lnz_f(n), etree_f(n), work_f(n), Lnz_c(n), etree_c(n), work_c(n))

        result_f = qdldl_etree(n, Ap, Ai, work_f, Lnz_f, etree_f)
        result_c = c_etree(int(n, c_int), Ap_c, Ai_c, work_c, Lnz_c, etree_c)
        if (result_f <= 0 .or. result_f /= result_c .or. any(Lnz_f /= Lnz_c) .or. &
            any(etree_f - 1_ip /= etree_c)) then
            error stop 'C and Fortran elimination trees differ'
        end if
        nnz_l = result_f

        allocate(Lp_f(n+1), Li_f(nnz_l), D_f(n), Dinv_f(n), Lx_f(nnz_l), &
                 Lp_c(n+1), Li_c(nnz_l), D_c(n), Dinv_c(n), Lx_c(nnz_l), &
                 bwork_f(n), bwork_c(n), iwork_f(3*n), iwork_c(3*n), &
                 fwork_f(n), fwork_c(n), rhs(n,solve_repeats), xf(n,solve_repeats), xc(n,solve_repeats))
        rhs = 1.0_wp

        result_f = qdldl_factor(n, Ap, Ai, Ax, Lp_f, Li_f, Lx_f, D_f, Dinv_f, &
                                Lnz_f, etree_f, bwork_f, iwork_f, fwork_f)
        result_c = c_factor(int(n,c_int), Ap_c, Ai_c, Ax, Lp_c, Li_c, Lx_c, D_c, Dinv_c, &
                            Lnz_c, etree_c, bwork_c, iwork_c, fwork_c)
        if (result_f /= result_c .or. result_f < 0) error stop 'C and Fortran factorization failed'
        if (any(Lp_f - 1_ip /= Lp_c) .or. any(Li_f - 1_ip /= Li_c)) &
            error stop 'C and Fortran L structures differ'
        if (.not. close_enough(Lx_f, Lx_c) .or. .not. close_enough(D_f, D_c)) &
            error stop 'C and Fortran factor values differ'

        do sample = 1, nsamples
            if (mod(sample,2) == 1) then
                call time_etree_f(etree_f_time(sample))
                call time_etree_c(etree_c_time(sample))
                call time_factor_f(factor_f_time(sample))
                call time_factor_c(factor_c_time(sample))
                call time_solve_f(solve_f_time(sample))
                call time_solve_c(solve_c_time(sample))
            else
                call time_etree_c(etree_c_time(sample))
                call time_etree_f(etree_f_time(sample))
                call time_factor_c(factor_c_time(sample))
                call time_factor_f(factor_f_time(sample))
                call time_solve_c(solve_c_time(sample))
                call time_solve_f(solve_f_time(sample))
            end if
        end do

        call print_row(matrix_name, n, nnz, nnz_l, 'etree  ', median(etree_f_time), median(etree_c_time))
        call print_row(matrix_name, n, nnz, nnz_l, 'factor ', median(factor_f_time), median(factor_c_time))
        call print_row(matrix_name, n, nnz, nnz_l, 'solve  ', median(solve_f_time), median(solve_c_time))
        deallocate(Ap, Ai, Ax, Lnz_f, etree_f, work_f, Lnz_c, etree_c, work_c, &
                   Lp_f, Li_f, D_f, Dinv_f, Lx_f, Lp_c, Li_c, D_c, Dinv_c, Lx_c, &
                   bwork_f, bwork_c, iwork_f, iwork_c, fwork_f, fwork_c, rhs, xf, xc, Ap_c, Ai_c)
    end do

    contains

    subroutine make_laplacian(m, n, Ap, Ai, Ax)
    integer,intent(in) :: m
    integer(ip),intent(out) :: n
    integer(ip),allocatable,intent(out) :: Ap(:), Ai(:)
    real(wp),allocatable,intent(out) :: Ax(:)
    integer(ip) :: j, p
    n = int(m,ip)**2
    allocate(Ap(n+1), Ai(n + 2*int(m,ip)*(m-1)), Ax(n + 2*int(m,ip)*(m-1)))
    p = 1
    do j = 1, n
        Ap(j) = p
        if (j > m) then
            Ai(p) = j - m
            Ax(p) = -1.0_wp
            p = p + 1
        end if
        if (mod(j-1,m) /= 0) then
            Ai(p) = j - 1
            Ax(p) = -1.0_wp
            p = p + 1
        end if
        Ai(p) = j
        Ax(p) = 4.0_wp
        p = p + 1
    end do
    Ap(n+1) = p
    end subroutine make_laplacian

    subroutine time_etree_f(seconds)
    real(wp),intent(out) :: seconds
    integer(int64) :: t0, t1, rate
    integer :: k
    call system_clock(count_rate=rate)
    call system_clock(t0)
    do k = 1, etree_repeats
        result_f = qdldl_etree(n, Ap, Ai, work_f, Lnz_f, etree_f)
    end do
    call system_clock(t1)
    if (result_f /= nnz_l) error stop 'Fortran etree timed run failed'
    seconds = real(t1-t0,wp) / real(rate*etree_repeats,wp)
    end subroutine time_etree_f

    subroutine time_etree_c(seconds)
    real(wp),intent(out) :: seconds
    integer(int64) :: t0, t1, rate
    integer :: k
    call system_clock(count_rate=rate)
    call system_clock(t0)
    do k = 1, etree_repeats
        result_c = c_etree(int(n,c_int), Ap_c, Ai_c, work_c, Lnz_c, etree_c)
    end do
    call system_clock(t1)
    if (result_c /= nnz_l) error stop 'C etree timed run failed'
    seconds = real(t1-t0,wp) / real(rate*etree_repeats,wp)
    end subroutine time_etree_c

    subroutine time_factor_f(seconds)
    real(wp),intent(out) :: seconds
    integer(int64) :: t0, t1, rate
    integer :: k
    call system_clock(count_rate=rate)
    call system_clock(t0)
    do k = 1, factor_repeats
        result_f = qdldl_factor(n, Ap, Ai, Ax, Lp_f, Li_f, Lx_f, D_f, Dinv_f, &
                                Lnz_f, etree_f, bwork_f, iwork_f, fwork_f)
    end do
    call system_clock(t1)
    if (result_f < 0) error stop 'Fortran factor timed run failed'
    seconds = real(t1-t0,wp) / real(rate*factor_repeats,wp)
    end subroutine time_factor_f

    subroutine time_factor_c(seconds)
    real(wp),intent(out) :: seconds
    integer(int64) :: t0, t1, rate
    integer :: k
    call system_clock(count_rate=rate)
    call system_clock(t0)
    do k = 1, factor_repeats
        result_c = c_factor(int(n,c_int), Ap_c, Ai_c, Ax, Lp_c, Li_c, Lx_c, D_c, Dinv_c, &
                            Lnz_c, etree_c, bwork_c, iwork_c, fwork_c)
    end do
    call system_clock(t1)
    if (result_c < 0) error stop 'C factor timed run failed'
    seconds = real(t1-t0,wp) / real(rate*factor_repeats,wp)
    end subroutine time_factor_c

    subroutine time_solve_f(seconds)
    real(wp),intent(out) :: seconds
    integer(int64) :: t0, t1, rate
    integer :: k
    xf = rhs
    call system_clock(count_rate=rate)
    call system_clock(t0)
    do k = 1, solve_repeats
        call qdldl_solve(n, Lp_f, Li_f, Lx_f, Dinv_f, xf(:,k))
    end do
    call system_clock(t1)
    seconds = real(t1-t0,wp) / real(rate*solve_repeats,wp)
    end subroutine time_solve_f

    subroutine time_solve_c(seconds)
    real(wp),intent(out) :: seconds
    integer(int64) :: t0, t1, rate
    integer :: k
    xc = rhs
    call system_clock(count_rate=rate)
    call system_clock(t0)
    do k = 1, solve_repeats
        call c_solve(int(n,c_int), Lp_c, Li_c, Lx_c, Dinv_c, xc(:,k))
    end do
    call system_clock(t1)
    seconds = real(t1-t0,wp) / real(rate*solve_repeats,wp)
    end subroutine time_solve_c

    function close_enough(a, b) result(ok)
    real(wp),intent(in) :: a(:), b(:)
    logical :: ok
    ok = maxval(abs(a-b)) <= 1.0e-12_wp * max(1.0_wp, maxval(abs(a)))
    end function close_enough

    function median(values) result(value)
    real(wp),intent(in) :: values(:)
    real(wp) :: value, sorted(size(values)), item
    integer :: a, b
    sorted = values
    do a = 2, size(sorted)
        item = sorted(a)
        b = a - 1
        do while (b >= 1)
            if (sorted(b) <= item) exit
            sorted(b+1) = sorted(b)
            b = b - 1
        end do
        sorted(b+1) = item
    end do
    value = sorted((size(sorted)+1)/2)
    end function median

    subroutine print_row(name, n, nnz, nnz_l, stage, tf, tc)
    character(len=*),intent(in) :: name, stage
    integer(ip),intent(in) :: n, nnz, nnz_l
    real(wp),intent(in) :: tf, tc
    write(*,'(A14,3I10,2X,A7,2F17.3,F12.2)') name, n, nnz, nnz_l, stage, &
        tf*1.0e6_wp, tc*1.0e6_wp, tc/tf
    end subroutine print_row

    end program benchmark
