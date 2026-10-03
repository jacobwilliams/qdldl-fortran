!*****************************************************************************************
!>
!  The Fortran side of the cross-check with the C QDLDL (see `run.sh`): writes a set
!  of test matrices (upper triangle, CSC, 0-based, natural order) to `matrices.txt`,
!  factors and solves them with [[qdldl_core]], and writes `D`, `L`, and the
!  solutions to `fortran_out.txt`, in the format of `crosscheck.c`.

    program crosscheck

    use qdldl_kinds, only: wp, ip
    use qdldl_core, only: qdldl_etree, qdldl_factor, qdldl_solve
    use, intrinsic :: iso_fortran_env, only: int64

    implicit none

    integer,parameter :: nmat = 6
    integer(ip),allocatable :: Ap(:), Ai(:), Lp(:), Li(:), Lnz(:), etree(:), iwork(:)
    real(wp),allocatable :: Ax(:), b(:), Lx(:), D(:), Dinv(:), fwork(:), x(:)
    logical,allocatable :: bwork(:)
    integer(ip) :: n, sumLnz, r, k
    integer :: im, um, uo
    integer(int64) :: state

    state = 4242
    open(newunit=um, file='matrices.txt', status='replace')
    open(newunit=uo, file='fortran_out.txt', status='replace')
    write(um,'(I0)') nmat
    do im = 1, nmat
        select case (im)
        case (1); call random_qd(5_ip, 0.5_wp)
        case (2); call random_qd(40_ip, 0.1_wp)
        case (3); call random_qd(200_ip, 0.02_wp)
        case (4); call random_qd(300_ip, 0.05_wp)
        case (5); call random_qd(1000_ip, 0.003_wp)
        case (6); call random_qd(2000_ip, 0.002_wp)
        case default; error stop 'unknown matrix'
        end select
        write(um,'(I0,1X,I0)') n, Ap(n+1) - 1
        write(um,'(*(I0,:,1X))') Ap - 1
        write(um,'(*(I0,:,1X))') Ai - 1
        write(um,'(*(ES26.17E3,:,1X))') Ax
        write(um,'(*(ES26.17E3,:,1X))') b

        allocate(Lp(n+1), Lnz(n), etree(n), iwork(3*n), D(n), Dinv(n), fwork(n), bwork(n))
        sumLnz = qdldl_etree(n, Ap, Ai, iwork, Lnz, etree)
        allocate(Li(max(1_ip, sumLnz)), Lx(max(1_ip, sumLnz)))
        r = qdldl_factor(n, Ap, Ai, Ax, Lp, Li, Lx, D, Dinv, Lnz, etree, bwork, iwork, fwork)
        x = b
        call qdldl_solve(n, Lp, Li, Lx, Dinv, x)
        write(uo,'(I0,1X,I0,1X,I0)') n, sumLnz, r
        write(uo,'(*(I0,:,1X))') (etree(k) - 1, k = 1, n)   ! (no parent: -1, as in C)
        write(uo,'(*(I0,:,1X))') Lp - 1
        write(uo,'(*(I0,:,1X))') Li(1:sumLnz) - 1
        write(uo,'(*(ES26.17E3,:,1X))') Lx(1:sumLnz)
        write(uo,'(*(ES26.17E3,:,1X))') D
        write(uo,'(*(ES26.17E3,:,1X))') x
        deallocate(Lp, Li, Lx, Lnz, etree, iwork, D, Dinv, fwork, bwork, x)
    end do
    close(um)
    close(uo)

    contains

    function uniform() result(u)
    !! A uniform random number in (-1, 1) (Park-Miller).
    real(wp) :: u   !! the number
    state = mod(16807_int64 * state, 2147483647_int64)
    u = 2.0_wp * real(state, wp) / 2147483647.0_wp - 1.0_wp
    end function uniform

    subroutine random_qd(nn, density)
    !! A random quasi-definite matrix (upper CSC, by columns, rows sorted), and a right-hand side.
    integer(ip),intent(in) :: nn        !! order
    real(wp),intent(in)    :: density   !! probability of an off-diagonal entry
    integer,allocatable :: s(:)
    real(wp),allocatable :: rowsum(:)
    integer(ip) :: i, j, q
    real(wp) :: v
    n = nn
    if (allocated(Ap)) deallocate(Ap, Ai, Ax, b)
    allocate(s(n), rowsum(n), Ap(n+1), Ai(n*n), Ax(n*n), b(n))
    do i = 1, n
        s(i) = merge(1, -1, uniform() > -0.2_wp)
        b(i) = uniform()
    end do
    rowsum = 0.0_wp
    q = 1
    do j = 1, n
        Ap(j) = q
        do i = 1, j - 1
            if (0.5_wp*(uniform() + 1.0_wp) < density) then
                v = uniform()
                Ai(q) = i; Ax(q) = v; q = q + 1
                if (s(i) == s(j)) then
                    rowsum(i) = rowsum(i) + abs(v); rowsum(j) = rowsum(j) + abs(v)
                end if
            end if
        end do
        Ai(q) = j; Ax(q) = 0.0_wp; q = q + 1      ! the diagonal (set below)
    end do
    Ap(n+1) = q
    ! diagonals: the row sums are complete only now
    do j = 1, n
        Ax(Ap(j+1)-1) = s(j) * (rowsum(j) + 0.1_wp + 0.5_wp*(uniform() + 1.0_wp))
    end do
    Ai = Ai(1:q-1); Ax = Ax(1:q-1)
    end subroutine random_qd

    end program crosscheck
!*****************************************************************************************
