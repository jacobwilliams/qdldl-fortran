!*****************************************************************************************
!>
!  Stage 3: the orderings. Each gives a valid permutation and the same solution;
!  on 2-D grid matrices, AMD's and RCM's fill is well below the natural order's
!  (the counts are recorded as a regression check). AMD gives the same ordering
!  whether or not it has to compress its work array.

    program test_orderings

    use qdldl_kinds, only: wp, ip
    use qdldl_module
    use qdldl_amd, only: amd
    use qdldl_test_utils

    implicit none

    real(wp),parameter :: res_tol = 100.0_wp * epsilon(1.0_wp)   !! scaled residual tolerance

    integer(ip) :: nnz_l(4)

    write(*,'(A)') 'test_orderings'

    ! the 2-D Laplacian of a 30 x 30 grid (n = 900)
    call grid_case('laplacian 30x30', 30_ip, 0_ip, nnz_l)
    ! regression: the nonzeros of L, by ordering (natural, rcm, amd, default)
    call check(nnz_l(1) == 26129, 'laplacian: natural fill')
    call check(nnz_l(2) == 18415, 'laplacian: RCM fill')
    call check(nnz_l(3) == 9331, 'laplacian: AMD fill')
    call check(nnz_l(2) < nnz_l(1), 'laplacian: RCM fill below natural')
    call check(nnz_l(3) < nnz_l(1) / 2, 'laplacian: AMD fill well below natural')
    call check(nnz_l(4) == nnz_l(3), 'laplacian: the default is AMD')

    ! a KKT matrix from it, with 300 constraints (n = 1200)
    call grid_case('kkt 30x30, 300 constraints', 30_ip, 300_ip, nnz_l)
    call check(nnz_l(1) == 304258, 'kkt: natural fill')
    call check(nnz_l(2) == 32872, 'kkt: RCM fill')
    call check(nnz_l(3) == 15533, 'kkt: AMD fill')
    call check(nnz_l(2) < nnz_l(1), 'kkt: RCM fill below natural')
    call check(nnz_l(3) < nnz_l(1) / 2, 'kkt: AMD fill well below natural')

    call test_amd_compression()
    call test_special_patterns()

    call finish_tests('test_orderings')

    contains

    !*************************************************************************************
    !>
    !  A grid matrix with every ordering, and a user ordering (the reverse of AMD's):
    !  valid permutations, and solutions with small residuals that agree.

    subroutine grid_case(name, m, nc, nnz)

    character(len=*),intent(in) :: name     !! name of the case
    integer(ip),intent(in)      :: m        !! grid size
    integer(ip),intent(in)      :: nc       !! number of constraints (0: the Laplacian alone)
    integer(ip),intent(out)     :: nnz(4)   !! the nonzeros of L (natural, rcm, amd, default)

    integer,parameter :: orders(4) = [qdldl_order_natural, qdldl_order_rcm, qdldl_order_amd, qdldl_order_default]
    character(len=*),parameter :: oname(4) = ['natural', 'rcm    ', 'amd    ', 'default']
    type(qdldl_type) :: ldl
    integer(ip),allocatable :: irow(:), icol(:), perm(:), rperm(:)
    real(wp),allocatable :: val(:), b(:), x(:), x0(:)
    integer(ip) :: n, k
    integer :: o, istat

    write(*,'(A)') ' '//name
    if (nc == 0) then
        call laplacian_2d(m, irow, icol, val)
    else
        call kkt_from_laplacian(m, nc, 1.0e-2_wp, irow, icol, val)
    end if
    n = m*m + nc
    allocate(b(n))
    do k = 1, n
        b(k) = real(mod(k, 7_ip), wp) - 3.0_wp
    end do

    do o = 1, 4
        call ldl%analyze(n, irow, icol, istat, ordering=orders(o))
        call check(istat == qdldl_success, name//': analyze '//trim(oname(o)))
        call ldl%get_permutation(perm)
        call check(qdldl_is_permutation(n, perm), name//': valid permutation '//trim(oname(o)))
        call ldl%factor(val, istat)
        call check(istat == qdldl_success, name//': factor '//trim(oname(o)))
        x = b
        call ldl%solve(x, istat)
        call check(residual_norm(irow, icol, val, x, b) < res_tol, name//': residual '//trim(oname(o)))
        if (o == 1) then
            x0 = x
        else
            call check(max_abs_diff(x, x0) <= 1.0e4_wp*epsilon(1.0_wp)*maxval(abs(x0)), &
                       name//': same solution '//trim(oname(o)))
        end if
        nnz(o) = ldl%nnz_l
        write(*,'(A,A8,A,I0)') '    nnz(L) ', trim(oname(o)), ': ', ldl%nnz_l
    end do

    ! a user ordering: the reverse of AMD's
    call ldl%analyze(n, irow, icol, istat, ordering=qdldl_order_amd)
    call ldl%get_permutation(perm)
    rperm = perm(n:1:-1)
    call ldl%analyze(n, irow, icol, istat, perm=rperm)
    call ldl%get_permutation(perm)
    call check(istat == qdldl_success .and. all(perm == rperm), name//': user permutation kept')
    call ldl%factor(val, istat)
    x = b
    call ldl%solve(x, istat)
    call check(istat == qdldl_success .and. residual_norm(irow, icol, val, x, b) < res_tol, &
               name//': residual with the user permutation')

    end subroutine grid_case
    !*************************************************************************************

    !*************************************************************************************
    !>
    !  [[amd]] with a work array too short for the elimination (so it compresses)
    !  gives the same permutation as with plenty of room.

    subroutine test_amd_compression()

    integer(ip),allocatable :: irow(:), icol(:), Ap(:), Ai(:), xadj(:), adj(:), cnt(:)
    integer(ip),allocatable :: iw(:), pe(:), len(:), nv(:), next(:), last(:), head(:), elen(:), degree(:), w(:)
    integer(ip),allocatable :: perm_room(:)
    real(wp),allocatable :: val(:)
    integer(ip) :: n, nadj, k, j, iwlen, pfree, ncmpa, trial
    integer :: istat

    write(*,'(A)') ' amd compression'
    call laplacian_2d(20_ip, irow, icol, val)
    n = 400
    ! upper CSC (the generator's entries are upper triangular, without duplicates)
    allocate(Ap(n+1), Ai(size(irow)), cnt(n+1))
    cnt = 0
    do k = 1, size(irow, kind=ip)
        cnt(icol(k)) = cnt(icol(k)) + 1
    end do
    Ap(1) = 1
    do j = 1, n
        Ap(j+1) = Ap(j) + cnt(j)
    end do
    cnt(1:n) = Ap(1:n)
    do k = 1, size(irow, kind=ip)
        Ai(cnt(icol(k))) = irow(k)
        cnt(icol(k)) = cnt(icol(k)) + 1
    end do
    call qdldl_symmetric_pattern(n, Ap, Ai, xadj, adj, istat)
    nadj = xadj(n+1) - 1

    allocate(pe(n), len(n), nv(n), next(n), last(n), head(n), elen(n), degree(n), w(n))
    do trial = 1, 2
        if (trial == 1) then
            iwlen = 3*nadj + 10*n    ! plenty of room
        else
            iwlen = nadj + 1         ! (almost) none
        end if
        if (allocated(iw)) deallocate(iw)
        allocate(iw(iwlen))
        iw(1:nadj) = adj(1:nadj)
        pe = xadj(1:n)
        len = xadj(2:n+1) - xadj(1:n)
        pfree = nadj + 1
        call amd(n, pe, iw, len, iwlen, pfree, nv, next, last, head, elen, degree, ncmpa, w)
        call check(qdldl_is_permutation(n, last), 'amd: valid permutation')
        call check(all(last(elen) == [(k, k = 1, n)]), 'amd: elen is the inverse permutation')
        if (trial == 1) then
            perm_room = last
            call check(ncmpa == 0, 'amd: no compression with room')
        else
            call check(ncmpa > 0, 'amd: compressions without room')
            call check(all(last == perm_room), 'amd: the same permutation either way')
        end if
    end do

    end subroutine test_amd_compression
    !*************************************************************************************

    !*************************************************************************************
    !>
    !  Patterns that test the edges of the orderings: a diagonal matrix (no edges, `n`
    !  components), a 1x1 matrix, an arrow matrix (one dense row and column), and two
    !  disconnected blocks.

    subroutine test_special_patterns()

    integer,parameter :: orders(2) = [qdldl_order_rcm, qdldl_order_amd]
    character(len=*),parameter :: oname(2) = ['rcm', 'amd']
    type(qdldl_type) :: ldl
    integer(ip),allocatable :: irow(:), icol(:), perm(:)
    real(wp),allocatable :: val(:), b(:), x(:)
    integer(ip) :: n, k
    integer :: o, c, istat
    character(len=:),allocatable :: name

    write(*,'(A)') ' special patterns'
    do c = 1, 4
        select case (c)
        case (1)  ! diagonal
            name = 'diagonal'
            n = 50
            irow = [(k, k = 1, n)]
            icol = irow
            val = [(real(k, wp), k = 1, n)]
        case (2)  ! 1x1
            name = '1x1'
            n = 1
            irow = [1_ip]; icol = [1_ip]; val = [-3.0_wp]
        case (3)  ! arrow: dense last row/column (natural order gives no fill; AMD must find it)
            name = 'arrow'
            n = 40
            irow = [(k, k = 1, n), (k, k = 1, n-1)]
            icol = [(k, k = 1, n), (n, k = 1, n-1)]
            val = [(4.0_wp, k = 1, n-1), real(2*n, wp), (1.0_wp, k = 1, n-1)]
        case (4)  ! two disconnected tridiagonal blocks
            name = 'two blocks'
            n = 30
            irow = [(k, k = 1, n), (k, k = 1, 14), (k, k = 16, n-1)]
            icol = [(k, k = 1, n), (k+1, k = 1, 14), (k+1, k = 16, n-1)]
            val = [(3.0_wp, k = 1, n), (-1.0_wp, k = 1, 14), (-1.0_wp, k = 16, n-1)]
        end select
        allocate(b(n))
        b = 1.0_wp
        do o = 1, 2
            call ldl%analyze(n, irow, icol, istat, ordering=orders(o))
            call ldl%get_permutation(perm)
            call check(istat == qdldl_success .and. qdldl_is_permutation(n, perm), &
                       name//': valid permutation '//oname(o))
            call ldl%factor(val, istat)
            x = b
            call ldl%solve(x, istat)
            call check(istat == qdldl_success .and. residual_norm(irow, icol, val, x, b) < res_tol, &
                       name//': residual '//oname(o))
            if (c == 3 .and. o == 2) call check(ldl%nnz_l == n - 1, 'arrow: AMD gives no fill')
        end do
        deallocate(b)
    end do

    end subroutine test_special_patterns
    !*************************************************************************************

    end program test_orderings
!*****************************************************************************************
