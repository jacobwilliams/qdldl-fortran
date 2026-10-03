!*****************************************************************************************
!>
!  Fill-reducing orderings for the \(LDL^T\) factorization: natural (none),
!  user-given, reverse Cuthill-McKee (RCM), and approximate minimum degree (AMD,
!  through the modernized SuiteSparse routine in [[qdldl_amd]]).
!
!  A permutation `perm` lists the original indices in pivot order: row `perm(k)`
!  of \(A\) is row `k` of \(PAP^T\). Its inverse `iperm` gives the position of
!  each original index: `iperm(perm(k)) = k`.
!
!  The orderings work on the full symmetric adjacency structure of the matrix
!  (both triangles, no diagonal, no duplicates), which [[qdldl_symmetric_pattern]]
!  builds from the upper triangle in CSC form. They allocate their own work
!  arrays (they run once per pattern, in the analysis).

    module qdldl_ordering

    use qdldl_kinds, only: ip
    use qdldl_core,  only: qdldl_success, qdldl_error_out_of_memory, qdldl_error_invalid_input, &
                           qdldl_error_overflow
    use qdldl_amd,   only: amd

    implicit none

    private

    ! orderings
    integer,parameter,public :: qdldl_order_natural = 0   !! no permutation
    integer,parameter,public :: qdldl_order_user    = 1   !! a permutation given by the caller
    integer,parameter,public :: qdldl_order_rcm     = 2   !! reverse Cuthill-McKee
    integer,parameter,public :: qdldl_order_amd     = 3   !! approximate minimum degree
    integer,parameter,public :: qdldl_order_default = -1  !! the default: AMD

    public :: qdldl_symmetric_pattern
    public :: qdldl_order_natural_perm
    public :: qdldl_rcm
    public :: qdldl_amd_order
    public :: qdldl_is_permutation

    contains
!*****************************************************************************************

!*****************************************************************************************
!>
!  The full symmetric adjacency structure (both triangles, without the diagonal) of
!  a matrix given by its upper triangle in CSC form without duplicates: the
!  neighbours of node `i` are `adj(xadj(i):xadj(i+1)-1)`.

    subroutine qdldl_symmetric_pattern(n, Ap, Ai, xadj, adj, istat)

    integer(ip),intent(in)               :: n              !! matrix order
    integer(ip),intent(in)               :: Ap(n+1)        !! column pointers of the upper triangle
    integer(ip),intent(in)               :: Ai(Ap(n+1)-1)  !! row indices of the upper triangle
    integer(ip),allocatable,intent(out)  :: xadj(:)        !! pointers into `adj` (size `n+1`)
    integer(ip),allocatable,intent(out)  :: adj(:)         !! the neighbours of each node
    integer,intent(out)                  :: istat          !! `qdldl_success`, `qdldl_error_out_of_memory`,
                                                           !! or `qdldl_error_overflow`

    integer(ip) :: i, j, p, q, noff
    integer(ip),allocatable :: next(:)

    istat = qdldl_success
    allocate(xadj(n+1), next(n), stat=i)
    if (i /= 0) then
        istat = qdldl_error_out_of_memory
        return
    end if

    ! count the neighbours of each node
    xadj = 0
    noff = 0
    do j = 1, n
        do p = Ap(j), Ap(j+1) - 1
            i = Ai(p)
            if (i /= j) then
                xadj(i) = xadj(i) + 1
                xadj(j) = xadj(j) + 1
                noff = noff + 1
            end if
        end do
    end do
    if (noff > (huge(noff) - 1_ip) / 2_ip) then
        istat = qdldl_error_overflow
        return
    end if

    ! pointers (cumulative sum)
    q = 1
    do i = 1, n
        p = xadj(i)
        xadj(i) = q
        next(i) = q
        q = q + p
    end do
    xadj(n+1) = q

    allocate(adj(max(1_ip, 2_ip*noff)), stat=i)
    if (i /= 0) then
        istat = qdldl_error_out_of_memory
        return
    end if
    do j = 1, n
        do p = Ap(j), Ap(j+1) - 1
            i = Ai(p)
            if (i /= j) then
                adj(next(i)) = j
                next(i) = next(i) + 1
                adj(next(j)) = i
                next(j) = next(j) + 1
            end if
        end do
    end do

    end subroutine qdldl_symmetric_pattern
!*****************************************************************************************

!*****************************************************************************************
!>
!  The natural ordering: `perm(k) = k`.

    pure subroutine qdldl_order_natural_perm(n, perm)

    integer(ip),intent(in)  :: n        !! matrix order
    integer(ip),intent(out) :: perm(n)  !! the permutation

    integer(ip) :: k

    do k = 1, n
        perm(k) = k
    end do

    end subroutine qdldl_order_natural_perm
!*****************************************************************************************

!*****************************************************************************************
!>
!  Whether `perm` is a permutation of `1..n`.

    function qdldl_is_permutation(n, perm) result(ok)

    integer(ip),intent(in) :: n        !! matrix order
    integer(ip),intent(in) :: perm(:)  !! the candidate permutation
    logical                :: ok       !! true if `perm` has size `n` and holds each of `1..n` once
                                       !! (false also if the work array can't be allocated)

    logical,allocatable :: seen(:)
    integer(ip) :: k
    integer :: stat

    ok = size(perm, kind=ip) == n
    if (.not. ok) return
    allocate(seen(n), stat=stat)
    if (stat /= 0) then
        ok = .false.
        return
    end if
    seen = .false.
    do k = 1, n
        if (perm(k) < 1 .or. perm(k) > n) then
            ok = .false.
            return
        end if
        if (seen(perm(k))) then
            ok = .false.
            return
        end if
        seen(perm(k)) = .true.
    end do

    end function qdldl_is_permutation
!*****************************************************************************************

!*****************************************************************************************
!>
!  The reverse Cuthill-McKee ordering, which reduces the bandwidth (and the profile)
!  of the matrix.
!
!  Each connected component is ordered by a breadth-first search from a
!  pseudo-peripheral node (found as in George and Liu's `FNROOT`, starting from a
!  node of minimum degree), visiting the neighbours of each node in order of
!  increasing degree. The whole ordering is then reversed.

    subroutine qdldl_rcm(n, xadj, adj, perm, istat)

    integer(ip),intent(in)  :: n                 !! matrix order
    integer(ip),intent(in)  :: xadj(n+1)         !! pointers into `adj`
    integer(ip),intent(in)  :: adj(xadj(n+1)-1)  !! the neighbours of each node (no diagonal, no duplicates)
    integer(ip),intent(out) :: perm(n)           !! the permutation
    integer,intent(out)     :: istat             !! `qdldl_success` or `qdldl_error_out_of_memory`

    integer(ip),allocatable :: deg(:), mark(:), queue(:), bydeg(:), cnt(:)
    logical,allocatable :: done(:)
    integer(ip) :: i, j, k, p, root, num, head, tail, first, stamp, scan, maxdeg
    integer(ip) :: nlvl, nlvl_new, last_start, ncomp, v, tmp
    integer :: stat

    istat = qdldl_success
    if (n == 0) return
    allocate(deg(n), mark(n), queue(n), bydeg(n), done(n), stat=stat)
    if (stat /= 0) then
        istat = qdldl_error_out_of_memory
        return
    end if

    maxdeg = 0
    do i = 1, n
        deg(i) = xadj(i+1) - xadj(i)
        maxdeg = max(maxdeg, deg(i))
    end do

    ! the nodes sorted by degree (counting sort), to pick a start node of
    ! minimum degree in each component in O(n) overall
    allocate(cnt(0:maxdeg+1), stat=stat)
    if (stat /= 0) then
        istat = qdldl_error_out_of_memory
        return
    end if
    cnt = 0
    do i = 1, n
        cnt(deg(i)+1) = cnt(deg(i)+1) + 1
    end do
    cnt(0) = 1
    do k = 1, maxdeg + 1
        cnt(k) = cnt(k) + cnt(k-1)
    end do
    do i = 1, n
        bydeg(cnt(deg(i))) = i
        cnt(deg(i)) = cnt(deg(i)) + 1
    end do

    done = .false.
    mark = 0
    stamp = 0
    num = 0
    scan = 1

    do while (num < n)

        ! an unnumbered node of minimum degree
        do while (done(bydeg(scan)))
            scan = scan + 1
        end do
        root = bydeg(scan)

        ! find a pseudo-peripheral node: the root of a level structure of
        ! (locally) maximum depth
        call level_structure(root, nlvl, last_start, ncomp)
        do
            if (ncomp == 1) exit
            ! a node of minimum degree in the last level
            v = queue(last_start)
            do k = last_start + 1, ncomp
                if (deg(queue(k)) < deg(v)) v = queue(k)
            end do
            call level_structure(v, nlvl_new, last_start, ncomp)
            if (nlvl_new <= nlvl) exit
            root = v
            nlvl = nlvl_new
        end do

        ! Cuthill-McKee from the root, numbering into perm(num+1:...)
        first = num + 1
        num = num + 1
        perm(num) = root
        done(root) = .true.
        head = first
        do while (head <= num)
            i = perm(head)
            head = head + 1
            tail = num
            do p = xadj(i), xadj(i+1) - 1
                j = adj(p)
                if (.not. done(j)) then
                    done(j) = .true.
                    num = num + 1
                    perm(num) = j
                end if
            end do
            ! sort the new nodes by increasing degree (insertion sort: the lists are short)
            do k = tail + 2, num
                tmp = perm(k)
                p = k - 1
                do while (p > tail)
                    if (deg(perm(p)) <= deg(tmp)) exit
                    perm(p+1) = perm(p)
                    p = p - 1
                end do
                perm(p+1) = tmp
            end do
        end do

    end do

    ! reverse
    do k = 1, n / 2
        tmp = perm(k)
        perm(k) = perm(n+1-k)
        perm(n+1-k) = tmp
    end do

    contains

        subroutine level_structure(r, nlevels, last_level_start, nreach)
        !! The level structure rooted at `r` of the unnumbered nodes reachable from it
        !! (breadth-first search, in `queue(1:nreach)`).
        integer(ip),intent(in)  :: r                  !! the root
        integer(ip),intent(out) :: nlevels            !! the number of levels
        integer(ip),intent(out) :: last_level_start   !! where the last level starts in `queue`
        integer(ip),intent(out) :: nreach             !! the number of nodes reached
        integer(ip) :: lstart, lend, q, u, w, pp
        stamp = stamp + 1
        queue(1) = r
        mark(r) = stamp
        nreach = 1
        lstart = 1
        nlevels = 0
        do
            lend = nreach
            nlevels = nlevels + 1
            last_level_start = lstart
            do q = lstart, lend
                u = queue(q)
                do pp = xadj(u), xadj(u+1) - 1
                    w = adj(pp)
                    if (.not. done(w) .and. mark(w) /= stamp) then
                        mark(w) = stamp
                        nreach = nreach + 1
                        queue(nreach) = w
                    end if
                end do
            end do
            if (nreach == lend) exit
            lstart = lend + 1
        end do
        end subroutine level_structure

    end subroutine qdldl_rcm
!*****************************************************************************************

!*****************************************************************************************
!>
!  The approximate minimum degree ordering, by [[amd]], with elbow room of 20% plus
!  `n` in its work array (so that compressions are rare).
!
!  As in SuiteSparse's C AMD, dense nodes (more than \(\max(16, 10\sqrt{n})\)
!  neighbours) are removed from the graph before ordering and placed last, in their
!  original order. Otherwise a dense row and column (common in KKT matrices) can make
!  the ordering take \(O(n^2)\) time.

    subroutine qdldl_amd_order(n, xadj, adj, perm, istat)

    integer(ip),intent(in)  :: n                 !! matrix order
    integer(ip),intent(in)  :: xadj(n+1)         !! pointers into `adj`
    integer(ip),intent(in)  :: adj(xadj(n+1)-1)  !! the neighbours of each node (no diagonal, no duplicates)
    integer(ip),intent(out) :: perm(n)           !! the permutation
    integer,intent(out)     :: istat             !! `qdldl_success`, `qdldl_error_out_of_memory`,
                                                 !! or `qdldl_error_overflow` (the work array's size)

    integer(ip),allocatable :: iw(:), pe(:), len(:), nv(:), next(:), head(:), elen(:), degree(:), w(:)
    integer(ip),allocatable :: last(:), newidx(:), oldidx(:)
    integer(ip) :: i, j, p, k, nadj, iwlen, pfree, ncmpa, dense, nr, q
    integer :: stat

    istat = qdldl_success
    if (n == 0) return

    nadj = xadj(n+1) - 1
    if (nadj > huge(nadj) - n - 1) then
        istat = qdldl_error_overflow
        return
    end if
    ! iwlen = 1.2*nadj + n, or as much as fits
    if (nadj / 5_ip > huge(nadj) - nadj - n - 1) then
        iwlen = nadj + n + 1
    else
        iwlen = nadj + nadj / 5_ip + n + 1
    end if

    allocate(iw(iwlen), pe(n), len(n), nv(n), next(n), head(n), elen(n), degree(n), w(n), &
             last(n), newidx(n), oldidx(n), stat=stat)
    if (stat /= 0) then
        istat = qdldl_error_out_of_memory
        return
    end if

    ! number the nodes that are not dense (newidx = 0 for a dense node)
    dense = max(16_ip, int(10.0*sqrt(real(n)), ip))
    nr = 0
    do i = 1, n
        if (xadj(i+1) - xadj(i) > dense) then
            newidx(i) = 0
        else
            nr = nr + 1
            newidx(i) = nr
            oldidx(nr) = i
        end if
    end do

    ! the graph without the dense nodes (it fits: it has at most nadj entries)
    q = 0
    do k = 1, nr
        i = oldidx(k)
        pe(k) = q + 1
        do p = xadj(i), xadj(i+1) - 1
            j = newidx(adj(p))
            if (j /= 0) then
                q = q + 1
                iw(q) = j
            end if
        end do
        len(k) = q + 1 - pe(k)
    end do
    pfree = q + 1

    if (nr > 0) call amd(nr, pe, iw, len, iwlen, pfree, nv, next, last, head, elen, degree, ncmpa, w)

    do k = 1, nr
        perm(k) = oldidx(last(k))
    end do
    k = nr
    do i = 1, n
        if (newidx(i) == 0) then
            k = k + 1
            perm(k) = i
        end if
    end do

    end subroutine qdldl_amd_order
!*****************************************************************************************

!*****************************************************************************************
    end module qdldl_ordering
!*****************************************************************************************
